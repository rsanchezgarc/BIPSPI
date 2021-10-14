import sys, os
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, matthews_corrcoef
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json

def loadResults( resultsPath, fnameResults):
  prefix= fnameResults.split("_")[0].split(".")[0]
  if ".lig" in fnameResults:
    chainType="l"
  else:
    chainType="r"
    
  scoresDf= pd.read_table(os.path.join(resultsPath, fnameResults), comment="#", sep="\s+", dtype={"resId":str, "chainId":str})
  return scoresDf
  
def get_single_chain_statistics(prefix, labels, scores, thr):
  scores= np.array(scores)
  labels= np.array(labels)
  try:
    roc_complex= roc_auc_score(labels, scores)
  except ValueError:
    roc_complex= np.nan
  probability_sorted_indexes = scores.argsort(axis=0)
  probability_sorted_indexes = probability_sorted_indexes[::-1]
  evalPoint= 0
  closestDist=9999
  for i in probability_sorted_indexes:
    if abs(scores[i] - thr)< closestDist:
      closestDist= abs(scores[i] - thr)	
      evalPoint=i
  print(evalPoint, scores[evalPoint])
  try:
    label_predictions= np.ones(scores.shape[0], dtype=np.int32)* np.min( labels)
    label_predictions[probability_sorted_indexes[0 : evalPoint]]= np.repeat(1, evalPoint)
    precision=  precision_score(labels[probability_sorted_indexes[0 : evalPoint]],
                                        label_predictions[probability_sorted_indexes[0 : evalPoint]])
    recall= recall_score(labels, label_predictions)
    mcc= matthews_corrcoef(labels, label_predictions)
#      print(sum(labels==1), sum(label_predictions==1),precisionAt[-1], recallAt[-1])
  except IndexError:
    precision= 0.0
    recall= 0.0
  summary= pd.DataFrame({"pdb":[prefix]})
  summary["auc_chains"]= [roc_complex] 
  summary["mcc"]= mcc
  summary["prec"]= precision
  summary["reca"]= recall
  return summary
  

def analyzeData(thr, resultsPath):
  prefixDict={}
  for fname in sorted(os.listdir(resultsPath)):
    if fname.split(".gz")[0].endswith(".rec") or fname.split(".gz")[0].endswith(".lig"):
      prefix= fname.split(".")[0]
      if not prefix in prefixDict:
        prefixDict[prefix]=[]
      prefixDict[prefix].append(fname)
  summaryList=[]
  labels=[]
  scores=[]
  for prefix in prefixDict:
    print(prefix)
    results0= loadResults(resultsPath, prefixDict[prefix][0])
    results1= loadResults(resultsPath, prefixDict[prefix][1])
    results= pd.concat([results0, results1], axis=0)
    labels.append(results["categ"])
    scores.append(results["prediction"])
    summary=  get_single_chain_statistics(prefix, results["categ"], results["prediction"], thr)
    summaryList.append(summary)
#    print(summary)
  summary= pd.concat(summaryList, axis=0, ignore_index=True)
  summary.reset_index(inplace=True)
  print(summary)
  print(summary.mean(axis=0))
  labels= np.concatenate(labels)
  scores= np.concatenate(scores)
  auc= roc_auc_score(labels, scores)
  scores[scores<thr]=-1
  scores[scores>=thr]=1
  recall= recall_score(labels, scores)
  precision= precision_score(labels, scores)
  mcc= matthews_corrcoef(labels, scores)
  print( "auc, mcc, precision, recall" )
  print(auc, mcc, precision, recall)
if __name__=="__main__":
  resultsPath="~/Tesis/rriPredMethod/data/ispred4/xgbWd/results/capri_mixed_2/"
  thr= float(sys.argv[1])
  if len(sys.argv)==3:
    resultsPath= sys.argv[2]
  resultsPath= os.path.expanduser(resultsPath)
  analyzeData(thr, resultsPath)
  
