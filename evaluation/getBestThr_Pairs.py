import sys, os
from sklearn.metrics import roc_auc_score, recall_score, precision_score, roc_curve
import pandas as pd
import numpy as np
from evaluation.evaluateScoresList import evaluateScoresLists

DO_ROC_CURVE= False

def loadResults( resultsPath, fnameResults):
  prefix= fnameResults.split("_")[0].split(".")[0]
    
  scoresDf= pd.read_table(os.path.join(resultsPath, fnameResults), comment="#", sep="\s+", 
              dtype={"structResIdL":str, "structResIdR":str, "chainIdL":str, "chainIdR":str})
  return scoresDf
  
def get_pairs_statistics(prefix, labels, scores):
  EVAL_PAIRS_AT= [ 2**2, 2**3]
  precisionAt=[]
  recallAt=[]
  scores= np.array(scores)
  labels= np.array(labels)
  try:
    roc_complex= roc_auc_score(labels, scores)
  except ValueError:
    roc_complex= np.nan
  probability_sorted_indexes = scores.argsort(axis=0)
  probability_sorted_indexes = probability_sorted_indexes[::-1]
  for evalPoint in EVAL_PAIRS_AT:
    try:
      label_predictions= np.ones(scores.shape[0])* np.min( labels)
      label_predictions[probability_sorted_indexes[0 : evalPoint]]= np.repeat(1, evalPoint)
      precisionAt.append( precision_score(labels[probability_sorted_indexes[0 : evalPoint]],
                                        label_predictions[probability_sorted_indexes[0 : evalPoint]]))
      recallAt.append( recall_score(labels, label_predictions))
#      print(sum(labels==1), sum(label_predictions==1),precisionAt[-1], recallAt[-1])
    except IndexError:
      precisionAt.append( 0.0)
      recallAt.append( 0.0)
  summary= pd.DataFrame({"pdb":[prefix]})
  summary["auc_pairs"]= [roc_complex] 
  for evalPoint, precisionAt, recallAt in zip(EVAL_PAIRS_AT,precisionAt, recallAt):
    summary["prec_%d"%evalPoint]= [precisionAt]
    summary["reca_%d"%evalPoint]= [recallAt]
  
  rocCurve= None
  if DO_ROC_CURVE:
    fpr, tpr, __= roc_curve(labels, scores)
    rocCurve= (fpr, tpr, roc_complex)
  return summary, rocCurve
  
def getOptimThr(resultsPath):
  allScores=[]
  allLabels=[]
  perComplexSummaries=[]
  rocCurves= []
  for fname in sorted(os.listdir(resultsPath)):
    if fname.endswith(".tab"):
      print(fname)
      results= loadResults( resultsPath, fname)
      if results is None: continue
      scores= list(results["prediction"].values)
      labels= list(results["categ"].values)
      summary, rocCurve= get_pairs_statistics(fname, labels, scores)
      if rocCurve: rocCurves.append(rocCurve)
      perComplexSummaries.append(summary)
#      print("%s %f"%(fname,roc_complex))
      allScores+= scores
      allLabels+= labels

  summary= pd.concat(perComplexSummaries, ignore_index=True)
  means= summary.mean(axis=0)
  summary= summary.append( summary.ix[summary.shape[0]-1,:],ignore_index=True )
  summary.ix[summary.shape[0]-1,0]=  "mean"
  summary.ix[summary.shape[0]-1, 1:]=  means
  evaluateScoresLists(allLabels, allScores, summary, summary.iloc[-1,1], None if not DO_ROC_CURVE else rocCurves)
        
if __name__=="__main__":
  '''
python -m evaluation.getBestThr_Pairs  ~/Tesis/rriPredMethod/data/bench5Data/newCodeData/results/mixed_2/
  '''
  resultsPath= "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/newCodeData/results_xgb/mixed_2/"
  if len(sys.argv)==2:
    resultsPath= sys.argv[1]
  getOptimThr(resultsPath)
  
