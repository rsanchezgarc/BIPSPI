from __future__ import print_function
import os, gc
import numpy as np
import json
import pandas as pd
import sys
from joblib import Parallel, delayed
from joblib import load as joblib_load
from joblib import dump as joblib_save
from sklearn.model_selection import train_test_split
from sklearn.utils import shuffle as shuffle_matrices
from sklearn.model_selection import GroupKFold

from Config import Configuration
from utils import getItemsFromList
from .crossValidationSplitter import subSplitFolds, mergeSplitFolds
from .processOneFold import trainAndTestOneFold, getResultsOutname
from .classifiers.gridSearch import tuneClassifier
from .getScopeGroups import getScopeGroups

from codifyComplexes.ComplexCodified import ComplexCodified, ComplexSeqStructCodified #required for joblib.load
pd.set_option('precision', 4)


RANDOM_STATE= 121
SKIP_EVALUATION=False

SKIP_LOWER_PREDICTION=True

SAMPLE_TRAINING_PAIRS= 0.7

class TrainAndTestWorker(Configuration):
  '''
    Class that performs train and test
  '''
  def __init__(self, trainDataPath, testPath, outputPath=None, nFolds=None, isLastStep=False,
               saveModelFname=None, verbose=True, numProc=1):
    '''
      builder
       
       :param trainDataPath: str. Path to a dir where training data files are stored
       :param testPath: str. Path to a dir where testing data files are stored
       :param outputPath: str. Path to a dir where predictions will be stored. If None, results will not be saved
                               and just performance evaluation will be carried out
       :param nFolds: int. Number of folds for k-fold cross-validation. If -1, leave-one-out will be performed.
                           If 0, testing set will be used as if it where independent. Cannot be 1

       :param isLastStep: bool. True if this train is the second step of a two steps workflow or the first one in one step workflow
       :param saveModelFname: str. A path where the final model, trained with all data will be saved. If None,
                                  model won't be saved              
       :param verbose: boolean. Whether or not print to stdout info
       :param numProc: int. Number of processes to use in parallel
    '''
    Configuration.__init__(self)  # Load configuration parameters such as path to programs
    
    parentPath, __ = os.path.split(trainDataPath)
    parentPath, stepName= os.path.split(parentPath)
    parentPath, __= os.path.split(parentPath)

    self.outputPath= outputPath
    self.saveModelFname= saveModelFname
    self.verbose= verbose
    self.numProc= numProc
    self.nFolds= nFolds
    self.isLastStep= isLastStep
    self.trainPath= trainDataPath
    
    trainFilesNames= sorted(os.listdir(self.trainPath))
    self.trainPrefixes=   sorted( set([ fname.split(".")[0] for fname in trainFilesNames]))

    self.testPath= testPath
    if not self.testPath is None:
      testFilesNames= sorted(os.listdir(self.testPath))
      self.testPrefixes= sorted(set([ fname.split(".")[0] for fname in testFilesNames ] ))
    else:
      self.testPrefixes= []

    self.data, self.prefixesUsedInModel= None, None # self.loadTrainingData(sharedMemoryPath=None) will be executed latter

    if self.verbose: 
      print("%d train complexes loaded."%(len(self.trainPrefixes)))   
    self.numTestComplexes= 0 if self.testPrefixes==None else len(self.testPrefixes)

  def loadTrainingData(self, sharedMemoryPath=None):
    '''
      :param sharedMemoryPath: str. Path to shared memory to store data matrices or None if regular memory will be used
                                    e.g. sharedMemoryPath = "shm://bipspi_train_cplxId_label_feats"
      :return data rank=2  tensor, second dimensions columns are : [complexesNumId, trainLabels, trainData0, trainData1]
            e.g. [12, 1, rasaL, rasaR, pssmL, ...]
              12 maps to self.trainPrefixes[12]
              1 is positive label
    '''
    if self.data is not None and self.prefixesUsedInModel is not None:
      return self.data, self.prefixesUsedInModel
    paralVerbose= 0
    if self.verbose: 
      print("Loading train data")
      paralVerbose= 1 
      
    readingNThreads= self.numProc
    if self.numProc > 12:
      readingNThreads= 12
    if sharedMemoryPath:
      import SharedArray as sa
      try:
        data = sa.attach(sharedMemoryPath )
        return data
      except OSError:
        pass

    trainComplexesLoaded= Parallel(n_jobs= readingNThreads, verbose=paralVerbose, backend="multiprocessing")(
                                                  delayed(loadOneTrainPrefixFun)( prefix, self.trainPath) for
                                                                                prefix in self.trainPrefixes)
    prefixesUsedInModel={}
    dataList= []
    labels= []
    prefixes= []
    complexesNumId=[]
    totalRows=0
    for complexNum, ppiComplex in enumerate(trainComplexesLoaded):
#      print(ppiComplex.prefix)
      data_d,data_t= ppiComplex.getData() #Data_d, data_t are either pairs L->R, R->L or seqL-structR, seqR-structL
      current_labels= ppiComplex.getLabels().astype(np.int32)
      if SAMPLE_TRAINING_PAIRS:
        data_idxs= np.random.randint(0, current_labels.shape[0], int(current_labels.shape[0]*SAMPLE_TRAINING_PAIRS) )
        data_d= None if isinstance(data_d, type(None) ) else data_d[data_idxs]
        data_t= None if isinstance(data_t, type(None) ) else data_t[data_idxs]
        current_labels = getItemsFromList(data_idxs, current_labels)

      prefixesUsedInModel[ppiComplex.prefix]= ppiComplex.prefixesInvolvedInCoding
      prefixes.append(ppiComplex.getPrefix())
      if not isinstance(data_d, type(None) ):
#        print( "d", data_d.shape)
        dataList.append(  data_d.astype(np.float32))
        labels.append( current_labels )
        complexesNumId+= [complexNum]* data_d.shape[0]
        totalRows+= data_d.shape[0]
      if not isinstance(data_t, type(None) ):
#        print( "t", data_t.shape)
        dataList.append( data_t.astype(np.float32))
        labels.append( current_labels )
        complexesNumId+= [complexNum]* data_t.shape[0]
        totalRows+= data_t.shape[0]
        
    nrows, ncols= ppiComplex.shape
    allDataMat= np.empty( (totalRows, 2+ncols), dtype=np.float32)
    allDataMat[:, 0]= complexesNumId
    del complexesNumId; gc.collect()
    allDataMat[:, 1]= np.concatenate(labels)
    del labels; gc.collect()    
    currentRow=0
    for mat in dataList:
      allDataMat[currentRow: currentRow+mat.shape[0], 2:]= mat 
      currentRow+= mat.shape[0]
    del dataList; gc.collect()
    if sharedMemoryPath:
      data = sa.create( sharedMemoryPath, allDataMat.shape )
      data[ ... ]= allDataMat
      del allDataMat, trainComplexesLoaded
    else:
      data= allDataMat
    return data, prefixesUsedInModel
    
  def doGridSearch(self):
    '''
      Performs a grid search. It will print out the parameters that achieved best score .

    '''
    self.data, self.prefixesUsedInModel = self.loadTrainingData(sharedMemoryPath=None)
    #TODO refactoring
    raise ValueError("Not re-implemented")
    
    trainData= self.data[:,2:]
    trainLabels= self.data[:,1]
    trainGroups= self.data[:,0]
    
    trainData, trainLabels, trainGroups= shuffle_matrices(trainData, trainLabels, trainGroups)
    tuneClassifier( trainData, trainLabels, trainGroups, nFolds= self.nFolds, ncpu=self.numProc,  
                    verboseLevel=1 if self.verbose else 0)
      
  def computeTrainAndTest(self):
    '''
      Performs train and testing. Returns a pandas.DataFrame which is a summary of performance evaluation
      :return summary: Pandas.DataFrame whose columns are:
          pdb  auc_pair  prec_50  reca_50  prec_100  reca_100  prec_500  reca_500   auc_l  prec_l  reca_l
          mcc_l   auc_r  prec_r  reca_r   mcc_r

    '''
    self.data, self.prefixesUsedInModel = self.loadTrainingData(sharedMemoryPath=None)

    if self.verbose: 
      print("Num pdbs to evaluate %d"%( self.numTestComplexes) )
    nTrain= len(self.trainPrefixes)

    if isinstance(self.nFolds, str) or (self.nFolds!=0 and self.nFolds>=-1):
      if len(set(self.trainPrefixes).intersection(self.testPrefixes)) !=nTrain:
        raise ValueError("Error. There are a different number of train and test data files. They must be equal to perform cross-validation")
        
      if isinstance(self.nFolds, str):
        if hasattr(self, "scopeFamiliesFname"): raise ValueError("Incompatible parameters in config file, if scopeFamiliesFname "+
                                                                 "defined, no folds file can be provided, just cross validation")
        trainPrefixes_testPrefixes_ori_list= self.splitAccordigFoldsFile(self.nFolds)
      else:
        trainPrefixes_testPrefixes_ori_list= self.splitKFoldCV()
        
      if not self.isLastStep:
        testIdxsAndTrainIdxs= subSplitFolds(trainPrefixes_testPrefixes_ori_list, self.trainPrefixes, random_state=RANDOM_STATE)
      else:
        testIdxsAndTrainIdxs= mergeSplitFolds(self.testPrefixes, trainPrefixes_testPrefixes_ori_list, self.prefixesUsedInModel)
    else:
      if hasattr(self, "scopeFamiliesFname"): raise ValueError("Incompatible parameters in config file, if train and test "+
                                                           "defined, no scope file should be provided")
      print("Performing testing with testing set (independency must be ensure by user")
      testIdxsAndTrainIdxs= [ (range(len(self.testPrefixes)), [range(nTrain)]) ]      
    summary = None
    if not SKIP_EVALUATION:
      listOfResultsSummary=[]
      testIdxsAndTrainIdxs= list(testIdxsAndTrainIdxs)
      nTrains=len(testIdxsAndTrainIdxs)
      for i, (testIdxs, trainIdxs_list) in enumerate(testIdxsAndTrainIdxs):
        if self.verbose==True and self.nFolds!=0: print("working on fold number %d / %d"%( (i +1), nTrains))
        listOfResultsSummary+= self.recursive_computeOneFold( testIdxs, trainIdxs_list)
      
      summary= pd.concat(listOfResultsSummary, ignore_index=True)
      means= summary.mean(axis=0)
      summary= summary.append( summary.ix[summary.shape[0]-1,:],ignore_index=True )
      summary.ix[summary.shape[0]-1,0]=  "mean"
      summary.ix[summary.shape[0]-1, 1:]=  means
      if self.verbose:
        print("Results:")
        print(summary.to_string(index=False))
    if not self.saveModelFname is None and not os.path.isfile(self.saveModelFname):
      if self.verbose:
        print("Training model with all data")
      model = self.computeOneFold( range(nTrain),  [], returnModel=True)         
      if self.verbose:    
        print("saving model at %s"%self.saveModelFname)
      joblib_save(model, self.saveModelFname)
    return summary

  def recursive_computeOneFold(self, testIdxs, trainIdxs_nested_list, i=0):   
    if  isinstance(trainIdxs_nested_list, tuple) or isinstance(trainIdxs_nested_list[0],int):
      oneResult_list = self.computeOneFold( trainIdxs_nested_list,  testIdxs, trainSubsetN= i, returnModel=False)
      return [oneResult.getPerformanceSummary() for oneResult in oneResult_list]
    else:
      listOfResultsSummary=[]
      for new_i, nested_list in enumerate(trainIdxs_nested_list):
          listOfResultsSummary+=self.recursive_computeOneFold( testIdxs, nested_list, i= (i, new_i))
      return listOfResultsSummary
  

  def computeOneFold(self, trainIdx,  testIdx, trainSubsetN=None, returnModel=False, evaluateJustFirstInTest=False):
    '''
      Trains and tests one fold of cross validation. Returns either a pandas.DataFrame which is a summary
      of performance evaluation or a modelObject[e.g.sklearn.ensemble.RandomForest ] object
      :param trainIdx: int[]. Indices of trainingComplexesLoaded list (self.trainComplexesLoaded) to be used in this fold
      :param testIdx:  int[]. Indices of testingComplexesPrefixes list (self.testPrefixes) to be test in this fold
      :param trainSubsetN: int Tuple. The numerical ids of the training split.
      :param returnModel: boolean. If True, returns modelObject[e.g.sklearn.ensemble.RandomForest ] object trained. Otherwise it returns
                                   a pandas.DataFrame which is a summary of performance evaluation
      :param evaluateJustFirstInTest: boolean. If true, evaluate just the first element in testIdx
      :return summary: Pandas.DataFrame  or modelObject[e.g.sklearn.ensemble.RandomForest ].
                    summary whose columns are:
                          pdb  auc_pair  prec_50  reca_50  prec_100  reca_100  prec_500  reca_500   auc_l  prec_l  reca_l
                          mcc_l   auc_r  prec_r  reca_r   mcc_r
                    modelObject[e.g.sklearn.ensemble.RandomForest ]. Trained model
    '''
    if self.verbose: print("Training with subset %s"%( str(trainSubsetN)))
    print(trainIdx,  testIdx)
    if self.isLastStep:
      trainSubsetN=None
    labelsAndTrainData= self.data[ np.isin(self.data[:, 0], trainIdx), 1: ]
    if returnModel:
      testPrefixes=[]
    else:
      if len(testIdx)!=0:
        testPrefixes= getItemsFromList( testIdx , self.testPrefixes)
        for testPrefix in testPrefixes:
          trainingMembers= getItemsFromList(trainIdx , self.trainPrefixes)
          crossValInfoFname= getResultsOutname(self.outputPath, testPrefix, trainSubsetN, suffix=".crossValInfo.json")
          with open(crossValInfoFname, 'w') as outfile:
            json.dump(trainingMembers, outfile)          
      else:
        return None
    
    if evaluateJustFirstInTest and not returnModel:
      testPrefixes= [ testPrefixes[0] ]

    if SKIP_LOWER_PREDICTION and self.isLastStep:
      testPrefixes= [ elem for elem in testPrefixes if elem[:4].isupper() ]

    results, model= trainAndTestOneFold( labelsAndTrainData, testPrefixes, trainSubsetN, self.testPath, self.outputPath,
                                         self.verbose, ncpu= self.numProc)
    if returnModel:                                  
      return model
    else:
      return results
    
      
  def splitAccordigFoldsFile(self,  foldsFname):
    '''
      :param foldsFname: str. A fname that contains the folds as a json
      [{"train":[2pr3,1a2k,1atn;1acb,1ahw], "test":[1A2K, 1ACB]}, {"train":[2pr1,1aak,1asn;1atb], "test":[1atn]} ]
    '''
    
    assert len( set(self.trainPrefixes).intersection( self.testPrefixes)) == len(self.testPrefixes)
    print("using custom folds")
    testIncluded=set([])
    trainPrefixes_testPrefixes_ori_list=[]
    with open( os.path.expanduser(foldsFname)) as f:
      data= json.load(f)
      for fold  in data:
        testPrefixes, trainPrefixes= fold["test"], fold["train"]
        trainPrefixes_testPrefixes_ori_list.append( (trainPrefixes, testPrefixes) )
        for testPrefix in testPrefixes: testIncluded.add(testPrefix )
      originalPrefixes= set([prefix.split("@")[0] for prefix in  self.testPrefixes])
      if not self.isLastStep:
        assert len(testIncluded.intersection( originalPrefixes
                  ))==len( originalPrefixes), "Error, in order to perform a second step"+\
                                        " all prefixes used for training should also be used for testing"+\
                                        "\n"+str(originalPrefixes.difference(testIncluded ))
      return trainPrefixes_testPrefixes_ori_list

  def checkIfIsMixedProtocol(self):
    return self.trainPrefixes[0].split("@")[0][-3:] in ["-#sr", "#sl"]
  
  def getGroupsForMixedProtocol(self, originalPrefixes):
    i=0
    groups_mixed=[]
    groups_mixed_dict={} 
    for prefix in originalPrefixes:
      if prefix[:-3] not in groups_mixed_dict:
        groups_mixed_dict[prefix[:-3]]=i
        groups_mixed.append(i)
        i+=1
      else:
        groups_mixed.append(groups_mixed_dict[prefix[:-3]])
    return groups_mixed
    
  def splitKFoldCV(self):

    originalPrefixes= sorted(set([prefix.split("@")[0] for prefix in  self.testPrefixes]))
    if hasattr(self, "scopeFamiliesFname") and os.path.isfile(self.scopeFamiliesFname):
      scopeGroups= list(getScopeGroups(originalPrefixes, self.scopeFamiliesFname,
                                       os.path.join(self.computedFeatsRootDir,
                                                  "seqStep/computedFeatures/seqStep/extractedSeqs/seqsData/")))
    else:
      scopeGroups=None

    if scopeGroups:
      groups= scopeGroups
    else:
      if self.checkIfIsMixedProtocol():
        groups= self.getGroupsForMixedProtocol( originalPrefixes)
      else:
        groups= range(len(originalPrefixes))
              
    if self.nFolds==-1:
      nFolds= len(set(groups))
      if self.verbose: print("Performing leave-one-out cross-validation")
    else:
      nFolds= self.nFolds
      if self.verbose: print("Performing %d-fold cross-validation"%self.nFolds)
    trainTest_idx_tuple_list_original= list(GroupKFold(n_splits= nFolds ).split(originalPrefixes, groups=groups))
    trainPrefixes_testPrefixes_ori_list=[]
    for trainIdxs, testIdxs in trainTest_idx_tuple_list_original:
      trainPrefixes= getItemsFromList(trainIdxs,  originalPrefixes)
      testPrefixes= getItemsFromList(testIdxs,  originalPrefixes)
      trainPrefixes_testPrefixes_ori_list.append( (trainPrefixes, testPrefixes) )
#    print(self.trainPrefixes)
#    print(trainPrefixes_testPrefixes_ori_list); raw_input("enter")
    return trainPrefixes_testPrefixes_ori_list
  
  @staticmethod  
  def loadPrefixFile( prefix, filesPath):
    '''
      Load a data file whose name startswith prefix and it is contained in filesPath

      :param prefix:  str. The prefix of the filename to be loaded. E.x. "1A2K"
      :param filesPath: str. The path where data files are contained
      :return complex_data: codifyComplexes.ComplexCodified.ComplexCodified class
      
    '''
    complexChunks= []
    for fname in sorted(os.listdir(filesPath)):
      if fname.split(".")[0]==prefix:
        complexChunks.append(joblib_load(os.path.join(filesPath, fname) ))
    assert len(complexChunks)==1
    return complexChunks[0]
    
  @staticmethod  
  def loadPrefixFilesIterator( prefix, filesPath):
    '''
      Load all data files whose name startswith prefix and it is contained in filesPath.
      Works as an iterator

      :param prefix:  str. The prefix of the filename to be loaded. E.x. "1A2K"
      :param filesPath: str. The path where data files are contained
      @yields complex_data: codifyComplexes.ComplexCodified.ComplexCodified class
      
    '''
    for fname in sorted(os.listdir(filesPath)):
      if fname.split(".")[0]==prefix:
        yield joblib_load(os.path.join(filesPath, fname) )
        
    
def loadOneTrainPrefixFun(prefix, filesPath):
  '''
    Load a data file whose name startswith prefix and it is contained in filesPath.
    It does the same that TrainAndTestWorker.loadPrefixFile() but can be used
    with Parallel delayed

    :param prefix:  str. The prefix of the filename to be loaded. E.x. "1A2K"
    :param filesPath: str. The path where data files are contained
    :return complex_data: codifyComplexes.ComplexCodified.ComplexCodified class
  '''
  return TrainAndTestWorker.loadPrefixFile( prefix, filesPath )


if __name__=="__main__":
  CONSIDER_AUGMENT=3
  trainDataPath, testPath, nFolds, dirToSave=sys.argv[1:]
  nFolds= int(nFolds)
  tt= TrainAndTestWorker( trainDataPath, testPath, nFolds = nFolds,  verbose = False)
  proposedFolds= tt.splitKFoldCV()

  for i,(fold_train, fold_test) in enumerate(proposedFolds):
    fold_test = [ x.split("_")[0] for x in fold_test if x[:4].isupper() ]
    fold_train= map(lambda x: x.split("_")[0] , fold_train)

    if CONSIDER_AUGMENT>0:
      fold_train_=[]
      for x in fold_train:
        fold_train_.append(x)
        if x[:4].islower():
          for j in range(1,1+CONSIDER_AUGMENT):
            fold_train_.append( x+"-%d"%j)
      fold_train= fold_train_
    foldData= {"train": list(fold_train), "test":list(fold_test)}
    with open(os.path.join(dirToSave, "fold_%d.json"%i), "w") as f:
      json.dump(foldData, f)
