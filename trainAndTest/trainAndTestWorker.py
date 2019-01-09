from __future__ import print_function
import os, sys, copy
import numpy as np
import random
import pandas as pd
from joblib import Parallel, delayed
from joblib import load as joblib_load
from joblib import dump as joblib_save
from sklearn.model_selection import KFold, GroupKFold
from sklearn.utils import shuffle as shuffle_matrices
from operator import itemgetter

from Config import Configuration
from .processOneFold import trainAndTestOneFold
from .classifiers.gridSearch import tuneClassifier
from codifyComplexes.ComplexCodified import ComplexCodified #requiered for joblib.load
pd.set_option('precision', 4)

PARALLEL_AT_COMPLEX_LEVEL= False #If True, each complex will be computed in a different process
                                 #otherwise training will be done using multiple threads (but not results evaluation).
                                 #True recommended when lots of memory available (8 gb per cpu aprox)

RANDOM_STATE= 121
TRAIN_WITH_BENCH= False
class TrainAndTestWorker(Configuration):
  '''
    Class that performs train and test
  '''
  def __init__(self, trainDataPath, testPath, outputPath=None, nFolds=None,
               saveModelFname=None, verbose=True, numProc=1):
    '''
      builder
       
       @param trainDataPath: str. Path to a dir where training data files are stored
       @param testPath: str. Path to a dir where testing data files are stored
       @param outputPath: str. Path to a dir where predictions will be stored. If None, results will not be saved
                               and just performance evaluation will be carried out
       @param nFolds: int. Number of folds for k-fold cross-validation. If -1, leave-one-out will be performed.
                           If 0, testing set will be used as if it where independent. Cannot be 1
       @param saveModelFname: str. A path where the final model, trained with all data will be saved. If None,
                                  model won't be saved                           
       @param verbose: boolean. Whether or not print to stdout info
       @param numProc: int. Number of processes to use in parallel
    '''
    Configuration.__init__(self)  # Load configuration parameters such as path to programs
    
    parentPath, __ = os.path.split(trainDataPath)
    parentPath, stepName= os.path.split(parentPath)
    parentPath, codifiedName= os.path.split(parentPath)
    parentPath, inputDataName= os.path.split(parentPath)
          
    self.outputPath= outputPath
    self.saveModelFname= saveModelFname
    self.verbose= verbose
    self.numProc= numProc
    self.nFolds= nFolds
    
    self.trainPath= trainDataPath
    trainFilesNames= sorted(os.listdir(self.trainPath))
    self.trainPrefixes=   sorted( set([ fname.split(".")[0] for fname in trainFilesNames]))

    self.testPath= testPath
    if not self.testPath is None:
      testFilesNames= sorted(os.listdir(self.testPath))
      self.testPrefixes= sorted(set([ fname.split(".")[0] for fname in testFilesNames ] ))
    else:
      self.testPrefixes= []

    self.data= self.loadTrainingData(sharedMemoryPath=None)

    if self.verbose: 
      print("%d train complexes loaded."%(len(self.trainPrefixes)))   
    self.numTestComplexes= 0 if self.testPrefixes==None else len(self.testPrefixes)

  def loadTrainingData(self, sharedMemoryPath=None):
  
    '''
      @param sharedMemoryPath: str. Path to shared memory to store data matrices or None if regular memory will be used
                                    e.g. sharedMemoryPath = "shm://bipspi_train_cplxId_label_feats"
      @return data rank=2  tensor, second dimensions columns are : [complexesNumId, trainLabels, trainData0, trainData1]
            e.g. [12, 1, rasaL, rasaR, pssmL, ...]
              12 maps to self.trainPrefixes[12]
              1 is positive label
    '''
    paralVerbose= 0
    if self.verbose: 
      print("Loading train data")
      paralVerbose= 1 
      
    readingNThreads= self.numProc
    if self.numProc > 18:
      readingNThreads= 18
    if sharedMemoryPath:
      import SharedArray as sa
      try:
        data = sa.attach(sharedMemoryPath )
        return data
      except OSError:
        pass
        
    trainComplexesLoaded= Parallel(n_jobs= readingNThreads, verbose=paralVerbose)(
                                                  delayed(loadOneTrainPrefixFun)( prefix, self.trainPath) for
                                                                                prefix in self.trainPrefixes)
    dataDir= []
    dataTrans= []
    labels= []
    prefixes= []
    complexesNumId=[]
    for complexNum, ppiComplex in enumerate(trainComplexesLoaded):
      data_d,data_t= ppiComplex.getData() #Data_d, data_t are either pairs L->R, R->L or seqL-structR, seqR-structL
      dataDir.append(  data_d)
      dataTrans.append( data_t)
      labels.append( ppiComplex.getLabels())
      prefixes.append(ppiComplex.getPrefix())
      complexesNumId+= [complexNum]* data_d.shape[0]

    dataDir= np.concatenate(dataDir).astype(np.float32)
    dataTrans= np.concatenate(dataTrans).astype(np.float32)
    labels= np.concatenate(labels)     
    trainData= np.concatenate([dataDir,dataTrans])
    del dataDir, dataTrans
    trainLabels= np.expand_dims(np.concatenate([labels,labels]), axis=-1)
    complexesNumId= np.array(complexesNumId, dtype= np.int64)
    complexesNumId= np.expand_dims( np.concatenate( [complexesNumId, complexesNumId]), axis=-1 )
    allDataMat= np.concatenate([complexesNumId, trainLabels, trainData], axis=-1)
    del trainData, complexesNumId, trainLabels
    if sharedMemoryPath:
#      import SharedArray as sa
      data = sa.create( sharedMemoryPath, allDataMat.shape )      
      data[ ... ]= allDataMat
      del allDataMat, trainComplexesLoaded
    else:
      data= allDataMat
    return data
    
  def doGridSearch(self):
    '''
      Performs a grid search. It will print out the parameters that achieved best score .

    '''
    raise ValueError("Not reimplemented")
    
    trainData= self.data[:,2:]
    trainLabels= np.concatenate([self.data[:,1], self.data[:,1]])
    trainGroups= np.concatenate([self.data[:,0],self.data[:,0]])
    
    trainData, trainLabels, trainGroups= shuffle_matrices(trainData, trainLabels, trainGroups)
    tuneClassifier( trainData, trainLabels, trainGroups, nFolds= self.nFolds, ncpu=self.numProc,  
                    verboseLevel=1 if self.verbose else 0)
      
  def computeTrainAndTest(self):
    '''
      Performs train and testing. Returns a pandas.DataFrame which is a summary of performance evaluation
      @return summary: Pandas.DataFrame whose columns are:
          pdb  auc_pair  prec_50  reca_50  prec_100  reca_100  prec_500  reca_500   auc_l  prec_l  reca_l
          mcc_l   auc_r  prec_r  reca_r   mcc_r

    '''
    if self.verbose: 
      print("Num pdbs to evaluate %d"%( self.numTestComplexes) )

    nTrain= len(self.trainPrefixes)
    nTest= len(self.testPrefixes)
    '''
      trainTest_idx_tuple_list is a list of tuples in which first element is a list of prefixes indices of self.trainPrefixes
      and the second element is a list of prefixes indices of self.testPrefixes
    '''
    if self.nFolds!=0 :#cross Validation:
      if len(set(self.trainPrefixes).intersection(self.testPrefixes)) !=nTrain:
        print("Error. There are a different number of train and test data files. They must be equal "+
              "to perform cross-validation")
        sys.exit(1)
      if self.nFolds==-1:
        trainTest_idx_tuple_list= KFold(n_splits= len(self.testPrefixes), random_state=RANDOM_STATE ).split(self.testPrefixes)
        print("Performing leave-one-out cross-validation")
      else:
        trainTest_idx_tuple_list= KFold(n_splits= self.nFolds, shuffle=True, random_state=RANDOM_STATE).split(self.testPrefixes)
        print("Performing %d-fold cross-validation"%self.nFolds)        
    else: #training and testing will be used independently
      if hasattr(self, "leaveExamplesOutFname") and os.path.isfile(self.leaveExamplesOutFname):
        trainTest_idx_tuple_list= list(self.splitAccordigFoldsFile( self.leaveExamplesOutFname ))      
      elif hasattr(self, "scopeFamiliesFname") and os.path.isfile(self.scopeFamiliesFname):
        trainTest_idx_tuple_list= list(self.splitAccordigScopePairs( self.scopeFamiliesFname ))
  #    trainTest_idx_tuple_list= list(self.splitAccordigScopePairs( self.scopeFamiliesFname ))
      else:
        print("Performing testing with testing set (independency must be ensure by user")
        trainTest_idx_tuple_list= [ (range(nTrain), range(nTest)) ]
        
#    trainTest_idx_tuple_list= reversed(list(trainTest_idx_tuple_list))
    listOfResults=[]    
    if PARALLEL_AT_COMPLEX_LEVEL:  
      argsList= [ ( trainIdx, testIdx, self.trainComplexesLoaded, self.testPrefixes,
                    self.testPath, self.outputPath, self.verbose) for (trainIdx, testIdx) 
                              in trainTest_idx_tuple_list] 

      listOfResults= Parallel(n_jobs=self.numProc)(delayed(computeOneFoldForParallel)(*args) for args in argsList)
      listOfResults= [ elem.getPerformanceSummary() for oneResult in listOfResults 
                                                    for elem in oneResult if not oneResult is None]
    else:
      trainTest_idx_tuple_list= list(trainTest_idx_tuple_list)
      nTrains= len( trainTest_idx_tuple_list)
      for i, (trainIdx, testIdx) in enumerate(trainTest_idx_tuple_list):
        if self.verbose==True and self.nFolds!=0:
          print("working on train number %d / %d"%( (i +1), nTrains))
        oneResult = self.computeOneFold( trainIdx,  testIdx, returnModel=False) 
        if oneResult is not None:
          listOfResults+= [elem.getPerformanceSummary() for elem in oneResult]
    summary= pd.concat(listOfResults, ignore_index=True)
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

  def splitAccordigScopePairs(self, scopeFname, ligandReceptorOrder=True):
    '''

    scope independency will be carried out as pairs of families, e.g (A-B) (A-B') are independent (A-B) (A'-B') not
    
    given a path to a file that contains scope info as bellow, splits train and test prefixes
    as leave-one-scope-family-out

    <<<< scopeFname example >>>
    1qfw    A      C      b.1.1.1:l.1.1.1;b.1.1.1 g.17.1.4;g.17.1.4       2
    2jel    B      A       b.1.1.1:b.1.1.2;b.1.1.1:b.1.1.2 d.94.1.1        2
    1avx    A      B       b.47.1.2        b.42.4.1        2

    @param ligandReceptorOrder: True if first scope column is for ligand and second is for receptor. 
                                if False first scope column is for receptor and second is for ligand.

    '''
    import itertools
    assert self.trainPrefixes== self.testPrefixes
    print("ensuring scope independency")
    prefixToScope={}
    scopeToPrefix={}
    with open(scopeFname) as f:
      for line in f:
        lineArray= line.split()
#        print(lineArray); raw_input("enter")
        if ligandReceptorOrder:
          prefix, chainsL, chainsR, scopesL, scopesR= lineArray[:5]
        else:
          prefix, chainsR, chainsL, scopesR, scopesL= lineArray[:5]
        prefix= prefix.upper()
        scopesL= set( itertools.chain.from_iterable([elem.split(":") for elem in scopesL.split(";")] ))
        scopesR= set( itertools.chain.from_iterable([elem.split(":") for elem in scopesR.split(";")] ))
        prefixToScope[prefix]= ( scopesL, scopesR )
    for i, prefix_i in enumerate(self.trainPrefixes):
      test_ix=[]
      train_ix=[]
      prefix_i=  prefix_i.split(":")[0].split("-")[0].upper()
      if prefix_i in prefixToScope:
        scope_L, scope_R= prefixToScope[prefix_i]
      for j, prefix_other in enumerate(self.trainPrefixes):
        prefix_other=  prefix_other.split(":")[0].split("-")[0].upper()
        if prefix_other in prefixToScope:
          scope_L_other, scope_R_other= prefixToScope[prefix_other]
          if (scope_L.intersection(scope_L_other) and scope_R.intersection(scope_R_other) or
              scope_L.intersection(scope_R_other) and scope_R.intersection(scope_L_other) ):
            test_ix+=[j]
          else:
            train_ix+=[j]
      yield ( np.array(train_ix ), np.array(test_ix ) )

      
  def splitAccordigScopeMonomers(self, scopeFname):
    '''
    scope independency will be carried out as nonomers of families, e.g (A-B) (A-B') are not independent
    
      scopeFname example: (Receptor, ligand)
    1qfw    IM      AB      b.1.1.1:l.1.1.1;b.1.1.1 g.17.1.4;g.17.1.4       2
    2jel    HL      P       b.1.1.1:b.1.1.2;b.1.1.1:b.1.1.2 d.94.1.1        2
    1avx    A       B       b.47.1.2        b.42.4.1        2

    '''
    import itertools
    assert self.trainPrefixes== self.testPrefixes    
    print("ensuring scope independency")
    prefixToScope={}
    scopeToPrefix={}
    with open(scopeFname) as f:
      for line in f:
        lineArray= line.split()
        prefix, chainsR, chainsL, scopesR, scopesL= lineArray[:5]
        prefix= prefix.upper()
        scopesL= set( itertools.chain.from_iterable([elem.split(":") for elem in scopesL.split(";")] ))
        scopesR= set( itertools.chain.from_iterable([elem.split(":") for elem in scopesR.split(";")] ))
        prefixToScope[prefix]= set([ elem for elem in scopesL.union( scopesR ) if elem!="NA" ])

    for i, prefix_i in enumerate(self.trainPrefixes):
      test_ix=[]
      train_ix=[]
      prefix_i= prefix_i.split(":")[0].split("-")[0].upper()
      if prefix_i in prefixToScope:
        scopes= prefixToScope[prefix_i]
      for j, prefix_other in enumerate(self.trainPrefixes):
        prefix_other= prefix_other.split(":")[0].split("-")[0].upper()
        if prefix_other in prefixToScope:
          scopes_other= prefixToScope[prefix_other]
          if scopes_other.intersection(scopes):
            test_ix+=[j]
          else:
            train_ix+=[j]
      yield ( np.array(train_ix ), np.array(test_ix ) )
      
  def splitAccordigFoldsFile(self,  foldsFname):
    '''
      foldsFname example:
      ------------------------
      2pr3 2pr3-A:B  
      2uuy 2uuy-A:B  2uuy-B:A 
      1m9x 1m9x-C:B  4dge-A:C  1ak4-D:A 
      1stf 1stf-I:E  3kfq-A:D  1yvb-I:A  1stf-E:I 
      1he1 1he1-A:C  1g4u-R:S  1he1-C:A 
      2ccl 2ccl-B:A  2b59-B:A  2ccl-B:A 
      1dev 1dev-A:B  1mk2-A:B 
      2e2d 2e2d-C:A  3v96-B:A  1uea-A:B  2e2d-A:C  3ma2-D:C  2j0t-D:A  1bqq-M:T  1gxd-A:C  4ilw-A:D 
      1izn 1izn-A:B  3aaa-A:B 
      ------------------------
    '''
    assert len( set(self.trainPrefixes).intersection( self.testPrefixes)) == len(self.testPrefixes)
    print("using custom folds")
    testAndLeaveOutList= []
    with open( os.path.expanduser(foldsFname)) as f:
      for line in f:
        lineArray= line.split()
        testAndLeaveOutList.append( lineArray[1:] )
    for i, oneFoldList in enumerate( testAndLeaveOutList ):
      test_ix= [ self.trainPrefixes.index(oneFoldList[0]) ]     
      for leaveOutExample in oneFoldList[1:]:
        try:
          idx= self.trainPrefixes.index(leaveOutExample)
          if idx not in test_ix:
            test_ix+=[ idx ]
        except ValueError:
          continue
      train_ix= sorted( set(range(len(self.trainPrefixes))).difference(test_ix))
      yield ( np.array(train_ix ), np.array(test_ix ) )


  def computeOneFold(self, trainIdx,  testIdx, returnModel=False, evaluateJustFirstInTest=False):
    '''
      Trains and tests one fold of cross validation. Returns either a pandas.DataFrame which is a summary
      of performance evaluation or a sklearn.ensemble.RandomForest object
      @param trainIdx: int[]. Indices of trainingComplexesLoaded list (self.trainComplexesLoaded) to be used in this fold
      @param testIdx:  int[]. Indices of testingComplexesPrefixes list (self.testPrefixes) to be test in this fold
      @param returnModel: boolean. If True, returns sklearn.ensemble.RandomForest object trained. Otherwise it returns 
                                   a pandas.DataFrame which is a summary of performance evaluation
      @param evaluateJustFirstInTest: boolean. If true, evaluate just the first element in testIdx
      @return summary: Pandas.DataFrame  or sklearn.ensemble.RandomForest.
                    summary whose columns are:
                          pdb  auc_pair  prec_50  reca_50  prec_100  reca_100  prec_500  reca_500   auc_l  prec_l  reca_l
                          mcc_l   auc_r  prec_r  reca_r   mcc_r
                    sklearn.ensemble.RandomForest. Trained model
    '''
    new_trainIdx= []

    if returnModel:
      testPrefixes=[]
    else:
      if len(testIdx)!=0:
        testPrefixes= itemgetter( *testIdx )( self.testPrefixes)
      else:
        return None

    print("Putative test prefixes %s"%(str(testPrefixes)))
    if isinstance(testPrefixes, str): testPrefixes= [testPrefixes]
    if not TRAIN_WITH_BENCH:
      testPrefixes=[ prefix for prefix in testPrefixes if prefix[1].isupper()]
      if len(testPrefixes)==0:
        return None


    if not TRAIN_WITH_BENCH:
      for i in trainIdx:
        prefix = self.trainPrefixes[i]
        if not prefix[1:3].isupper():
          new_trainIdx.append(i)
      trainIdx= new_trainIdx
      assert len(trainIdx)>2, "Error, TRAIN_WITH_BENCH=False, but all complexes belong to benchmark (capital letters)"
    labelsAndTrainData= self.data[ np.isin(self.data[:, 0], trainIdx), 1: ]
#    print( itemgetter(* trainIdx)(self.trainPrefixes), len(trainIdx) )
#    print( set(self.trainPrefixes).difference(set(itemgetter(* trainIdx)(self.trainPrefixes))) )
#    print(labelsAndTrainData.shape)
#    raw_input("enter")

    if evaluateJustFirstInTest and not returnModel:
      testPrefixes= [ testPrefixes[0] ]

    results, model= trainAndTestOneFold( labelsAndTrainData, testPrefixes, self.testPath, self.outputPath, 
                                         self.verbose, ncpu= self.numProc)
    if returnModel:
      return model
    else:
      return results

    
  @staticmethod  
  def loadPrefixFile( prefix, filesPath):
    '''
      Load a data file whose name startswith prefix and it is contained in filesPath

      @param prefix:  str. The prefix of the filename to be loaded. E.x. "1A2K"
      @param filesPath: str. The path where data files are contained
      @return complex_data: codifyComplexes.ComplexCodified.ComplexCodified class
      
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

      @param prefix:  str. The prefix of the filename to be loaded. E.x. "1A2K"
      @param filesPath: str. The path where data files are contained
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

    @param prefix:  str. The prefix of the filename to be loaded. E.x. "1A2K"
    @param filesPath: str. The path where data files are contained
    @return complex_data: codifyComplexes.ComplexCodified.ComplexCodified class
  '''
  return TrainAndTestWorker.loadPrefixFile( prefix, filesPath )

def computeOneFoldForParallel(trainIdx,  testIdx, trainComplexesLoaded, testPrefixes, testPath,
                              outputPath, verbose):
  '''
    Trains and tests one fold of cross validation. Returns a pandas.DataFrame which is a summary
    of performance evaluation
    @param trainIdx: int[]. Indices of trainComplexesLoaded list  to be used in this fold
    @param testIdx:  int[]. Indices of testPrefixes list to be test in this fold
    @param trainComplexesLoaded: codifyComplexes.ComplexCodified[] 
    @param testPrefixes: str[]. A list of pdb identifiers
    @param testPath: str. A path where testing pdb files are located
    @param outputPath: str. A path to the directory where results will be saved
    @param verbose: boolean. Whether or not print to stdout info
    @return summary: Pandas.DataFrame which is a performance summary whose columns are:
                        pdb  auc_pair  prec_50  reca_50  prec_100  reca_100  prec_500  reca_500   auc_l  prec_l  reca_l
                        mcc_l   auc_r  prec_r  reca_r   mcc_r
  '''
  trainComplexes= itemgetter( *trainIdx )(trainComplexesLoaded)

  if len(testIdx)!=0 or not returnModel:
    testPrefixes= itemgetter( *testIdx )( testPrefixes)
  else:
    testPrefixes= []
    
  if isinstance(testPrefixes, str): testPrefixes= [testPrefixes]

  results, __= trainAndTestOneFold( trainComplexes, testPrefixes, testPath, outputPath, 
                                    verbose, ncpu= 1)

  return results
