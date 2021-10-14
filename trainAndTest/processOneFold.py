from __future__ import print_function
import os
import gc
import hashlib
import random

from Config import Configuration
from joblib import load as joblib_load, dump as joblib_save, Parallel, delayed
import numpy as np

from utils import checkFreeMemory, getFileSize, getTotalMemory, tryToRemove
from .resultsManager import ResultsManager
from codifyComplexes.ComplexCodified import ComplexSeqStructCodified

# from .classifiers.randomForest import  trainMethod, predictMethod
from .classifiers.xgBoost import trainMethod, predictMethod

SAMPLE_TRAIN_EXAMPLES=False
MAX_SAMPLING_PAIRS=500000

conf=Configuration()

def findFullTestPPIName(testPrefix, testPath):
  for fname in sorted(os.listdir(testPath)):
    if fname.startswith(testPrefix):
      return os.path.join(testPath, fname)
  return None

def getDataForTestFromPrefix(testPrefix, testPath):
  '''
    Load a data file whose name startswith testPrefix and it is contained in testPath.
    Returns a tuple with all data needed to perform predictions and testing

    :param prefix:  str. The prefix of the filename to be loaded. E.g. "1A2K"
    :param filesPath: str. The path where data files are contained
    :return (data_d, data_t, ppiComplex.getLabels(), ppiComplex.getIds())
          data_d: np.array (n,m). A np.array that can be feed to the classifier. Each row represents
                                  a pair of amino acids in direct form (first ligand aa second receptor aa)
          data_l: np.array (n,m). A np.array that can be feed to the classifier. Each row represents
                                  a pair of amino acids in transpose form (first receptor aa second ligand aa)
          ppiComplex.getLabels(): np.array which contains the labels (-1, 1 ) of each row (pair of amino acids)
          ppiComplex.getIds(): pandas.DataFrame whose columns are:
                    chainIdL resIdL resNameL chainIdR resIdR resNameR categ
  '''

  ppiComplex = joblib_load( findFullTestPPIName(testPrefix, testPath))
  isSeqStruct = isinstance(ppiComplex, ComplexSeqStructCodified)
  data_d, data_t = ppiComplex.getData()

  labels= ppiComplex.getLabels()
  ids=  ppiComplex.getIds()
  if SAMPLE_TRAIN_EXAMPLES and testPrefix[:4].islower():

    gc.collect()
    condition= labels>0
    posLabelsIdx= np.where(condition)[0]
    negIdexIdx= np.where(~ condition)[0]

    nToSample= min(MAX_SAMPLING_PAIRS, len(negIdexIdx) )

    if nToSample==MAX_SAMPLING_PAIRS:
      print("Random sampling for %s"%(testPrefix))
      random_state = abs(hash(testPrefix.split("@")[0].split("#s")[0]))  # ensure that all complexes with the same prefix are equally sampled. Required for results average
      random_state = random_state // 2 ** 32 - 1
      random.seed(random_state)
      np.random.seed(random_state)

      negIdexIdx= np.random.choice(negIdexIdx, size= nToSample, replace=False)
      selectIdxs= np.concatenate([posLabelsIdx, negIdexIdx])
      selectIdxs= np.array( [  int(l) for l in selectIdxs if not np.isnan(l) and not np.isnan(labels[l]) ],  dtype=np.int32)
      assert len(selectIdxs)>0, "Error, empty selectIdx"+str(testPrefix)+"\n"+str(ids.head())
      data_d= data_d[selectIdxs,:] if data_d is not None else None
      data_t= data_t[selectIdxs,:] if data_t is not None else None
      labels= labels[selectIdxs]
      ids= ids.iloc[selectIdxs, :].reset_index()
      gc.collect()
      random.seed(None)
      np.random.seed(None)

  return isSeqStruct, (data_d, data_t, labels, ids)


def getResultsOutname(outputPath, testPrefix, trainSubsetN=None, suffix=".res.tab.gz"):
  testPrefix = testPrefix.split("@")[0].split(".")[0].split("_")[0]
  if trainSubsetN is None:  # In second step just one file is wanted to be written
    return os.path.join(outputPath, testPrefix + suffix)
  else:  # In first step out of two, one file per subset will be written
    trainSubsetN_to_str = "@" + prettifyNestedTuple(trainSubsetN)[:-1] if trainSubsetN else ""
    return os.path.join(outputPath, testPrefix + trainSubsetN_to_str + suffix)


def prettifyNestedTuple(nestT):
  if not isinstance(nestT, tuple):
    return str(nestT)
  else:
    result = ""
    for elem in nestT:
      result += prettifyNestedTuple(elem)
    result += "@"
  return result


def getOriginalToActualPrefixs(newPrefixes):
  originalPrefixes = [prefix.split("@")[0] if "@" in prefix else prefix for prefix in newPrefixes]
  newPrefixToOriginalPrefix = {}
  originalPrefixToNewPrefix = {}
  for newId, originalId in zip(newPrefixes, originalPrefixes):
    newPrefixToOriginalPrefix[newId] = originalId
    try:
      originalPrefixToNewPrefix[originalId].add(newId)
    except KeyError:
      originalPrefixToNewPrefix[originalId] = set([newId])

  return originalPrefixToNewPrefix, newPrefixToOriginalPrefix


def estimateRequiredMemoryPerComplex(prefixesList, path, diskToMemoryFactor=8000):
  if len(prefixesList)==0: return 0.5*diskToMemoryFactor
  memSize=0
  for i, prefix in enumerate( prefixesList):
    fname= findFullTestPPIName(prefix, path)
    # print("fname for memory estimation", fname, getFileSize(fname))
    memSize+= getFileSize(fname)*diskToMemoryFactor
  return memSize/float(i+1)


def trainAndTestOneFold(trainData, testPrefixes, trainSubsetN, testPath, outputPath, verbose=False, ncpu=1):
  '''
    Trains and tests one fold
     
     :param trainData: a numpy array for training with first column labels and the others are features
     :param testPrefixes: str[]. A list that contains prefixes for all complexes to be tested
     :param trainSubsetN: int Tuple. The numerical ids of the training split.
     :param testPath: str. Path to a dir where testing data files are stored
     :param outputPath: str. Path to a dir where predictions will be stored. None if results will not be saved
     :param verbose: boolean. Whether or not print to stdout info
     :param ncpu: int. Number of cpu's to use in parallel
  '''

  testPrefixesNotEvaluated = []
  originalTestPrefixToNewPrefix, __ = getOriginalToActualPrefixs(testPrefixes)
  alreadyComputedPrefixes_and_outnames= []
  for testPrefix in originalTestPrefixToNewPrefix:
    if outputPath is not None:
      outName = getResultsOutname(outputPath, testPrefix, trainSubsetN)
      if verbose and os.path.isfile(outName):
        print("Complex already computed: %s" % (outName))
        alreadyComputedPrefixes_and_outnames.append(  (testPrefix, outName) )
      else:
        testPrefixesNotEvaluated.append((testPrefix, outName))
    else:
      testPrefixesNotEvaluated.append((testPrefix, None))

  modelo = None
  modelFname= os.path.join(conf.tmp, hashlib.md5("".join(sorted(testPrefixes))).hexdigest()+str(trainSubsetN)+"bipspi2.pckl")

  resultsForEvaluation_list=[]
  if len(testPrefixesNotEvaluated) > 0 or len(testPrefixes) == 0:
    if verbose:
      print("Testing:", [ x[0] for x in testPrefixesNotEvaluated])
      verboseLevel = 1
    else:
      verboseLevel = 0

    if os.path.exists(modelFname):
      print("Loading classifier")
      modelo= joblib_load(modelFname)
    else:
      print("Training classifier")
      modelo = trainMethod(trainData[:, 1:], trainData[:, 0], verboseLevel=verboseLevel, ncpu=ncpu)
      joblib_save(modelo, modelFname)
    del trainData
    gc.collect()
    if verbose: print("Classifier fitted.")
    
    expectedSize= estimateRequiredMemoryPerComplex(testPrefixesNotEvaluated, testPath)
    freeMem= checkFreeMemory()
    nJobs= int(max(1, min(ncpu, freeMem/expectedSize, len(testPrefixesNotEvaluated))))
    print("Free memory for predictOnePrefix: %s GB. Njobs: %s (%s expected size)"%(freeMem, nJobs, expectedSize))

    resultsForEvaluation_list= Parallel(n_jobs=nJobs)(delayed(predictOnePrefix)(originalTestPrefixToNewPrefix[testPrefix],
                                                                      modelo, outName, testPath)
                                      for testPrefix, outName in testPrefixesNotEvaluated )
    gc.collect()

  expectedSize= estimateRequiredMemoryPerComplex(alreadyComputedPrefixes_and_outnames, testPath)
  freeMem= checkFreeMemory()
  nJobs= int(max(1, min(ncpu, freeMem/expectedSize, len(alreadyComputedPrefixes_and_outnames))))     
  resultsForEvaluation_list+= Parallel(n_jobs=nJobs)(delayed(loadExistingResults)( testPrefix, outName,)
                                    for testPrefix, outName in alreadyComputedPrefixes_and_outnames )
    
  if len(resultsForEvaluation_list)>0:
    freeMem = checkFreeMemory()
    totMem= getTotalMemory()
    usedMem= totMem-freeMem
    nJobs = int(max(1, min(ncpu, freeMem / (usedMem/(1+len(resultsForEvaluation_list))))))
    print("Free memory for evaluateOneResultObj: %s GB. Njobs: %s" % (freeMem, nJobs))
    Parallel(n_jobs=nJobs)(delayed(evaluateOneResultObj)(testPrefix, resultObj, False)
                           for testPrefix, resultObj in resultsForEvaluation_list)
    finalResults= zip(*resultsForEvaluation_list)[1]
  else:
    finalResults=[]
  del resultsForEvaluation_list
  tryToRemove(modelFname)
  return finalResults, modelo

def loadExistingResults(prefix, outName):
  return prefix, ResultsManager.loadExistingResults(outName)
    
def predictOnePrefix(testPrefixes, modelo, outName, testingDataPath):
  '''
  testPrefixes: Prefixes of a given complex. E.g. 1ACB@1 1ACB@2 for 1ACB, or directly 1ACB
  '''
  prob_predictionsDir = []
  prob_predictionsTrans = []
  initPrefix = list(testPrefixes)[0].split("@")[0].split("_")[0].split("#s")[0]
  for testPrefix in testPrefixes:
    print("predicting %s"%(testPrefix))
    currentPrefix = testPrefix.split("@")[0].split("_")[0].split("#s")[0]
    assert initPrefix == currentPrefix, "Error, all test prefixes provided to predictAndEvaluateOnePrefix should share the same root prefix"
    isSeqStruct, ppiComplexData = getDataForTestFromPrefix(testPrefix, testingDataPath)
    testDataDirect, testDataTrans, testlabels, testPairsIds = ppiComplexData
    if not isinstance(testDataDirect, type(None)):
      preds = predictMethod(modelo, testDataDirect)
      del testDataDirect
      gc.collect()
      prob_predictionsDir += [preds]
    #      print(ResultsManager(testPrefix, preds, None, testPairsIds, doAverageScore=not isSeqStruct).getPerformanceSummary())
    if not isinstance(testDataTrans, type(None)):
      preds = predictMethod(modelo, testDataTrans)
      del testDataTrans
      gc.collect()
      prob_predictionsTrans += [preds]
  #      print(ResultsManager(testPrefix, preds, None, testPairsIds, doAverageScore=not isSeqStruct).getPerformanceSummary())
  if len(prob_predictionsDir) > 0:
    prob_predictionsDir = np.stack(prob_predictionsDir, axis=1)
    prob_predictionsDir = np.mean(prob_predictionsDir, axis=1)
  else:
    prob_predictionsDir = None
  if len(prob_predictionsTrans) > 0:
    prob_predictionsTrans = np.stack(prob_predictionsTrans, axis=1)
    prob_predictionsTrans = np.mean(prob_predictionsTrans, axis=1)
  else:
    prob_predictionsTrans = None
  testPrefix = os.path.basename(outName).split(".")[0].split("_")[0]

  if isSeqStruct:
    prob_predictionsDir = prob_predictionsDir if not isinstance(prob_predictionsDir,
                                                                type(None)) else prob_predictionsTrans
    prob_predictionsTrans = None
  resultObj = ResultsManager(testPrefix, prob_predictionsDir, prob_predictionsTrans, testPairsIds, doAverageScore=not isSeqStruct)
  del testlabels, testPairsIds

  if not os.path.isfile(outName):
    resultObj.writeResults(outName)

  del prob_predictionsDir, prob_predictionsTrans
  gc.collect()
  return testPrefix, resultObj

def evaluateOneResultObj(testPrefix, resultObj, verbose):
  if verbose is True: print("Evaluating predictions for %s" % (testPrefix))
  resultObj.getFullEvaluation()
  if verbose is True: print(resultObj)
