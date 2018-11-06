from __future__ import print_function
import itertools
import sys, os
import inspect
import numpy as np
from joblib import load as joblib_load

from .resultsManager import ResultsManager
#from .classifiers.randomForest import  trainMethod, predictMethod
from .classifiers.xgBoost import  trainMethod, predictMethod

def getDataForTestFromPrefix( testPrefix, testPath ):
  '''
    Load a data file whose name startswith testPrefix and it is contained in testPath.
    Returns a tuple with all data needed to perform predictions and testing

    @param prefix:  str. The prefix of the filename to be loaded. E.g. "1A2K"
    @param filesPath: str. The path where data files are contained
    @return (data_d, data_t, ppiComplex.getLabels(), ppiComplex.getIds())
          data_d: np.array (n,m). A np.array that can be feed to the classifier. Each row represents
                                  a pair of amino acids in direct form (first ligand aa second receptor aa)
          data_l: np.array (n,m). A np.array that can be feed to the classifier. Each row represents
                                  a pair of amino acids in transpose form (first receptor aa second ligand aa)
          ppiComplex.getLabels(): np.array which contains the labels (-1, 1 ) of each row (pair of amino acids)
          ppiComplex.getIds(): pandas.DataFrame whose columns are:
                    chainIdL resIdL resNameL chainIdR resIdR resNameR categ            
  '''
  for fname in sorted(os.listdir(testPath)):
    if fname.startswith(testPrefix):
      ppiComplex= joblib_load(os.path.join(testPath, fname) )
      data_d,data_t= ppiComplex.getData()      
      return (data_d, data_t, ppiComplex.getLabels(), ppiComplex.getIds())
  
def trainAndTestOneFold(trainData, testPrefixes, testPath, outputPath, verbose=False, ncpu= 1):
  '''
    Trains and tests one fold
     
     @param trainData: a numpy array for training with first column labels and the others are features
     @param testPrefixes: str[]. A list that contains prefixes for all complexes to be tested
     @param testPath: str. Path to a dir where testing data files are stored
     @param outputPath: str. Path to a dir where predictions will be stored
     @param verbose: boolean. Whether or not print to stdout info
     @param ncpu: int. Number of cpu's to use in parallel
  ''' 
  resultsForEvaluation_list= []
  testPrefixesNotEvaluated=[]
  finalResults=[]
  for testPrefix in testPrefixes:
    if outputPath is not None:
      outName= os.path.join( outputPath, testPrefix+".res.tab") 
      if verbose and os.path.isfile(outName):
        print("Complex already computed: %s"%(outName))
        resultsForEvaluation_list.append( (testPrefix, ResultsManager.loadExistingResults(outName) ) )
      else:
        testPrefixesNotEvaluated.append( testPrefix )
    else:
      testPrefixesNotEvaluated.append( testPrefix )
  modelo=None
  if len(testPrefixesNotEvaluated)> 0 or len(testPrefixes)==0:
    if verbose:  
      print("Testing:", testPrefixesNotEvaluated)
      print("Training classifier")
      verboseLevel=1
    else: 
      verboseLevel=0
    
    modelo= trainMethod(trainData[:, 1:], trainData[:, 0], verboseLevel= verboseLevel, ncpu= ncpu)
    if verbose==True: print ("Classifier fitted.")

    for testPrefix in testPrefixesNotEvaluated:
      prob_predictionsDir_list= []
      prob_predictionsTrans_list=[]
      testlabels_list=[]
      testPairsIds_list=[]
      if verbose==True: print("Computing predictions for %s"%(testPrefix))
      testDataDirect, testDataTrans, testlabels, testPairsIds= getDataForTestFromPrefix( testPrefix, testPath )
      prob_predictionsDir= predictMethod(modelo, testDataDirect)
      prob_predictionsTrans= predictMethod(modelo,testDataTrans)
      resultEval= ResultsManager(testPrefix, prob_predictionsDir, prob_predictionsTrans, testPairsIds)
      if verbose==True: print("Evaluating predictions of %s"%(testPrefix))
      resultEval.getFullEvaluation()
      if verbose==True: print(resultEval)
##      raw_input("press enter")
      finalResults.append( resultEval )
      if not outputPath is None:
        outName= os.path.join(outputPath, testPrefix+".res.tab") 
        if not os.path.isfile(outName):
          if verbose==True: print("Saving results at %s"%(outName))
          resultEval.writeResults(outName)
      
  for testPrefix, resultEval in resultsForEvaluation_list:
    if verbose==True: print("Evaluating predictions for %s"%(testPrefix))
    resultEval.getFullEvaluation()
    if verbose==True: print(resultEval)
    finalResults.append( resultEval )

  return finalResults, modelo
  
  
