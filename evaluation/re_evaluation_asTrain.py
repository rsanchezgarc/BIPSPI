import sys, os
from sklearn.metrics import roc_auc_score, recall_score, precision_score
import pandas as pd
import numpy as np
from evaluation.evaluateScoresList import evaluateScoresLists
from trainAndTest.resultsManager import ResultsManager

SKIP_LOWER=True


def reEvaluate(resultsPath, newResultsPath=None):

  allScores=[]
  allLabels=[]
  listDf=[]
  roc_aucs_list= []
  for fname in sorted(os.listdir(resultsPath)):
    if not (fname.endswith(".res.tab") or fname.endswith(".res.tab.gz")): continue
    if SKIP_LOWER and fname[:4].islower(): continue
    res= ResultsManager.loadExistingResults ( os.path.join(resultsPath, fname) )
    oneComplexDf= res.getFullEvaluation()
    scores, labels= res.bindingSiteScoresAndLabels()
    allScores+= scores 
    allLabels+= labels
    roc_aucs_list.append( roc_auc_score( labels, scores))
    print("%s: %f auc binding site"%(fname, roc_aucs_list[-1]))
    listDf.append(oneComplexDf)
    if not newResultsPath is None:
      outName= os.path.join(newResultsPath, fname)
      res.writeResults(outName)
  resDf= pd.concat( listDf, axis=0, ignore_index=True)
  means= resDf.mean(axis=0)
  resDf= resDf.append( resDf.ix[resDf.shape[0]-1,:],ignore_index=True )
  resDf.ix[resDf.shape[0]-1,0]=  "mean"
  resDf.ix[resDf.shape[0]-1, 1:]=  means
  print("Results:")
  print(resDf.to_string(index=False))
  evaluateScoresLists(allLabels, allScores, resDf, np.mean(roc_aucs_list), None)  
  
if __name__=="__main__":  
  resultsPath= os.path.expanduser("~/Tesis/rriPredMethod/data/bench5Data/newCodeData/results/mixed_2")
  outPath=None
  if len(sys.argv)==2:
    resultsPath= os.path.expanduser(sys.argv[1])
  elif len(sys.argv)==3:
    resultsPath= os.path.expanduser(sys.argv[1])
    outPath= os.path.expanduser(sys.argv[2])
  reEvaluate(resultsPath, outPath)

