
from __future__ import print_function
import os, sys, copy
import numpy as np
import random
import pandas as pd
from joblib import Parallel, delayed
from joblib import load as joblib_load
from joblib import dump as joblib_save
from sklearn.model_selection import KFold
from operator import itemgetter

from .processOneFold import trainAndTestOneFold, getDataForClassifierFromComplexes
from .classifiers.gridSearch import  tuneClassifier
 
from codifyComplexes.ComplexCodified import ComplexCodified
pd.set_option('precision', 4)

PARALLEL_AT_COMPLEX_LEVEL= False #If True, each complex will be computed in a different one thread process
                                 #otherwise training will be done using multiple threads (but not results evaluation).
                                 #True recommended when lots of memory available (8 gb per cpu aprox)
                                 
RANDOM_STATE= 121
                                 
class TrainAndTestWorker(object):
  '''
    Class that performs train and test
  '''
  def __init__(self, trainDataPath, testPath, outputPath=None, nFolds=None,
               saveModelFname=None, verbose=True, numProc=1):
    '''
      builder
       
       @param trainDataPath: str. Path to a dir where training data files are stored
       @param testPath: str. Path to a dir where testing data files are stored
       @param outputPath: str. Path to a dir where predictions will be stored
       @param nFolds: int. Number of folds for k-fold cross-validation. If -1, leave-one-out will be performed.
                           If 0, testing set will be used as if it where independent. Cannot be 1
       @param saveModelFname: str. A path where the final model, trained with all data will be saved. If None,
                                  model won't be saved                           
       @param verbose: boolean. Whether or not print to stdout info
       @param numProc: int. Number of processes to use in parallel
    '''
    
    prefixPath, __ = os.path.split(trainDataPath)
    prefixPath, stepName= os.path.split(prefixPath)
    prefixPath, codifiedName= os.path.split(prefixPath)
    prefixPath, inputDataName= os.path.split(prefixPath)
          
    self.outputPath= outputPath
    self.saveModelFname= saveModelFname
    self.verbose= verbose
    self.numProc= numProc
    self.nFolds= nFolds
    
    self.trainPath= trainDataPath
    trainFilesNames= sorted(os.listdir(self.trainPath))
    self.trainPrefixes=   sorted(
    set([ fname.split(".")[0] for fname in trainFilesNames]))

    self.testPath= testPath
    if not self.testPath is None:
      testFilesNames= sorted(os.listdir(self.testPath))
      self.testPrefixes= sorted(set([ fname.split(".")[0] for fname in testFilesNames ] ))
    else:
      self.testPrefixes= []
    paralVerbose= 0
    if self.verbose: 
      print("Loading train data")
      paralVerbose= 1 
      
    readingNThreads= self.numProc
    if self.numProc > 12:
      readingNThreads= 12
         
    self.trainComplexesLoaded= Parallel(n_jobs= readingNThreads, verbose=paralVerbose)(
                                                delayed(loadOneTrainPrefixFun)( prefix, self.trainPath) for
                                                                              prefix in self.trainPrefixes)
    if self.verbose: 
      print("%d train complexes loaded."%(len(self.trainComplexesLoaded)))   
    self.numTestComplexes= 0 if self.testPrefixes==None else len(self.testPrefixes)

  def doGridSearch(self):
    '''
      Performs a grid search. It will print out the parameters that achieved best score .

    '''
    trainComplexes= copy.copy(self.trainComplexesLoaded)

    random.seed(121)
    random.shuffle(trainComplexes)
    random.seed()
    getDataForClassifierFromComplexes(trainComplexes)    
    dataDir,dataTrans, labels, groups = getDataForClassifierFromComplexes(trainComplexes)
    trainData= np.concatenate([dataDir,dataTrans])
    trainLabels= np.concatenate([labels,labels])
    trainGroups= np.concatenate([groups,groups])
    dataDir,dataTrans, labels, groups = (None, None, None, None)
    randIndices= np.random.choice(trainData.shape[0],trainData.shape[0], False)
    trainData= trainData[randIndices, ...]
    trainLabels= trainLabels[randIndices, ...]
    trainGroups= trainGroups[randIndices, ...]
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
      trainTest_idx_tuple_list= [ (range(nTrain), range(nTest)) ]
      print("Performing testing with testing set (independency must be ensure by user")
    
#    scopeFamFile= os.path.expanduser("~/Tesis/rriPredMethod/data/bench5Data/newCodeData/DB5.scope.tsv")
#    trainTest_idx_tuple_list= self.splitAccordigScopePairs( scopeFamFile )
#    trainTest_idx_tuple_list= self.splitAccordigScopeMonomers( scopeFamFile)

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

  def splitAccordigScopePairs(self, scopeFname):
    '''
      scopeFname example: (Receptor, ligand)
    1qfw    IM      AB      b.1.1.1:l.1.1.1;b.1.1.1 g.17.1.4;g.17.1.4       2
    2jel    HL      P       b.1.1.1:b.1.1.2;b.1.1.1:b.1.1.2 d.94.1.1        2
    1avx    A       B       b.47.1.2        b.42.4.1        2

    '''
    import itertools
    print("ensuring scope independency")
    prefixToScope={}
    scopeToPrefix={}
    with open(scopeFname) as f:
      for line in f:
        prefix, chainsR, chainsL, scopesR, scopesL, benchVersion= line.split()
        prefix= prefix.upper()
        scopesL= set( itertools.chain.from_iterable([elem.split(":") for elem in scopesL.split(";")] ))
        scopesR= set( itertools.chain.from_iterable([elem.split(":") for elem in scopesR.split(";")] ))
        prefixToScope[prefix]= ( scopesL, scopesR )
        
    for i, prefix_i in enumerate(self.trainPrefixes):
      test_ix=[]
      train_ix=[]
      if prefix_i in prefixToScope:
        scope_L, scope_R= prefixToScope[prefix_i]
      for j, prefix_other in enumerate(self.trainPrefixes):
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
      scopeFname example: (Receptor, ligand)
    1qfw    IM      AB      b.1.1.1:l.1.1.1;b.1.1.1 g.17.1.4;g.17.1.4       2
    2jel    HL      P       b.1.1.1:b.1.1.2;b.1.1.1:b.1.1.2 d.94.1.1        2
    1avx    A       B       b.47.1.2        b.42.4.1        2

    '''
    import itertools
    print("ensuring scope independency")
    prefixToScope={}
    scopeToPrefix={}
    with open(scopeFname) as f:
      for line in f:
        prefix, chainsR, chainsL, scopesR, scopesL, benchVersion= line.split()
        prefix= prefix.upper()
        scopesL= set( itertools.chain.from_iterable([elem.split(":") for elem in scopesL.split(";")] ))
        scopesR= set( itertools.chain.from_iterable([elem.split(":") for elem in scopesR.split(";")] ))
        prefixToScope[prefix]= set([ elem for elem in scopesL.union( scopesR ) if elem!="NA" ])
    for i, prefix_i in enumerate(self.trainPrefixes):
      test_ix=[]
      train_ix=[]
      if prefix_i in prefixToScope:
        scopes= prefixToScope[prefix_i]
      for j, prefix_other in enumerate(self.trainPrefixes):
        if prefix_other in prefixToScope:
          scopes_other= prefixToScope[prefix_other]
          if scopes_other.intersection(scopes):
            test_ix+=[j]
          else:
            train_ix+=[j]
      yield ( np.array(train_ix ), np.array(test_ix ) )
      
  def computeOneFold(self, trainIdx,  testIdx, returnModel=False):
    '''
      Trains and tests one fold of cross validation. Returns either a pandas.DataFrame which is a summary
      of performance evaluation or a sklearn.ensemble.RandomForest object
      @param trainIdx: int[]. Indices of trainingComplexesLoaded list (self.trainComplexesLoaded) to be used in this fold
      @param testIdx:  int[]. Indices of testingComplexesPrefixes list (self.testPrefixes) to be test in this fold
      @param returnModel: boolean. If True, returns sklearn.ensemble.RandomForest object trained. Otherwise it returns 
                                   a pandas.DataFrame which is a summary of performance evaluation
      @return summary: Pandas.DataFrame  or sklearn.ensemble.RandomForest.
                    summary whose columns are:
                          pdb  auc_pair  prec_50  reca_50  prec_100  reca_100  prec_500  reca_500   auc_l  prec_l  reca_l
                          mcc_l   auc_r  prec_r  reca_r   mcc_r
                    sklearn.ensemble.RandomForest. Trained model
    '''
#    print(trainIdx, testIdx, self.testPrefixes)
    trainComplexes= itemgetter( *trainIdx )(self.trainComplexesLoaded)
    if returnModel:
      testPrefixes=[]
    else:
      if len(testIdx)!=0:
        testPrefixes= itemgetter( *testIdx )( self.testPrefixes)
      else:
        return None
      
    if isinstance(testPrefixes, str): testPrefixes= [testPrefixes]

    results, model= trainAndTestOneFold( trainComplexes, testPrefixes, self.testPath, self.outputPath, 
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
