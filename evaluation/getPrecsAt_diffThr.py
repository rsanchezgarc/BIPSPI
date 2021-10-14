import sys, os
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, matthews_corrcoef
import pandas as pd
import numpy as np
import json

EVALUATE_CAPITAL_ONLY=True

DO_SEQ_AVERAGING=False

#BINDING_CMAPS_PATH="/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/newCodeData/computedFeatures/common/contactMapsBinding"
BINDING_CMAPS_PATH=None # if None, contact maps contained in the results file will be used, otherwise, new contacts will be added from the files in the path

def averageScores(scoresDf):
  labels_list=[]
  scores_list=[]
  for chain in scoresDf["chainId"].unique():
    df= scoresDf.loc[scoresDf["chainId"]==chain,:]
    result= [ ((int(elem), ""), (score, label)) if elem[-1].isdigit() else ((int(elem[:-1]), elem[-1]), (score, label)) 
                                for elem,score,label in zip(df["resId"], df["prediction"], df["categ"]) ]
    result.sort(key= lambda x: x[0])
    result= zip(* result)[1]
    scores, labels= zip(* result)
    scores= list(np.convolve(scores, np.array([1, 3, 1])/5.0, mode='same')+ np.array(scores))
    labels_list+= labels
    scores_list+= scores
  return scores_list, labels_list
  
def loadResults( resultsPath, fnameResults, cMapsPath=BINDING_CMAPS_PATH):
  prefix= fnameResults.split("_")[0].split(".")[0]
  if ".lig" in fnameResults:
    chainType="l"
  else:
    chainType="r"
    
  scoresDf= pd.read_table(os.path.join(resultsPath, fnameResults), comment="#", sep="\s+", dtype={"resId":str, "chainId":str})
  if not cMapsPath is None:
    newCmapSet=set([])
    for fname in os.listdir(cMapsPath):
      if ((chainType=="l" and "_l_" in fname) or (chainType=="r" and "_r_" in fname)) and fname.startswith(prefix):
        df= pd.read_table(os.path.join(cMapsPath,fname),sep='\s+', header='infer', comment="#", 
                          dtype= {"chainIdL":str, "chainIdR":str, "resIdL":str, "resIdR":str,
                                        "chainId":str, "resId":str,  "resId":str})
        for i in range(df.shape[0]):
          chainId, resId, categ= df.iloc[i,:]
          if categ==1:
            newCmapSet.add((chainId, resId))
    for chainId, resId in newCmapSet:
      scoresDf.loc[(scoresDf["chainId"]==chainId) & (scoresDf["resId"]==resId),"categ"]=1
  return scoresDf
  
def get_single_chain_statistics(prefix, labels, scores):
  EVAL_PAIRS_AT= [ 2**3, 2**4]
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
      label_predictions= np.ones(scores.shape[0], dtype=np.int32)* np.min( labels)
      
      label_predictions[probability_sorted_indexes[0 : evalPoint]]= np.repeat(1, evalPoint)
      precisionAt.append( precision_score(labels[probability_sorted_indexes[0 : evalPoint]],
                                        label_predictions[probability_sorted_indexes[0 : evalPoint]]))
      recallAt.append( recall_score(labels, label_predictions))
#      print(sum(labels==1), sum(label_predictions==1),precisionAt[-1], recallAt[-1])
    except IndexError:
      precisionAt.append( 0.0)
      recallAt.append( 0.0)
  summary= pd.DataFrame({"pdb":[prefix]})
  summary["auc_chains"]= [roc_complex] 
  for evalPoint, precisionAt, recallAt in zip(EVAL_PAIRS_AT,precisionAt, recallAt):
    summary["prec_%d"%evalPoint]= [precisionAt]
    summary["reca_%d"%evalPoint]= [recallAt]
  return summary
  
def getRecallAtPrec(resultsPath, targetPrec=0.666, jsonFname=None, useSeqAvera= DO_SEQ_AVERAGING):
  allScores=[]
  allLabels=[]
  perComplexSummaries=[]
  for fname in sorted(os.listdir(resultsPath)):
    if EVALUATE_CAPITAL_ONLY and fname[1:4].islower(): continue
    print(fname)
    if ".rec" in fname or ".lig" in fname:
      results= loadResults( resultsPath, fname)
      if results is None: continue
      if useSeqAvera:
        scores, labels= averageScores(results)
      else:
        scores= list(results["prediction"].values)
        labels= list(results["categ"].values)

      summary= get_single_chain_statistics(fname, labels, scores)
      perComplexSummaries.append(summary)
#      print("%s %f"%(fname,roc_complex))
      allScores+= scores
      allLabels+= labels

  summary= pd.concat(perComplexSummaries, ignore_index=True)
  means= summary.mean(axis=0)
  summary= summary.append( summary.ix[summary.shape[0]-1,:],ignore_index=True )
  summary.ix[summary.shape[0]-1,0]=  "mean"
  summary.ix[summary.shape[0]-1, 1:]=  means
    
  roc_auc= roc_auc_score(allLabels, allScores)
  
  allScores= np.array(allScores)
  allScoresUnique= np.unique(np.round(allScores, decimals=4))
  indices= np.argsort(allScoresUnique)[1:-1]
  # import matplotlib.pyplot as plt
  #  plt.hist(allScoresUnique[:-2], bins=50)
  #  plt.show()
  bestThr=-1
  bestMcc= 0
  bestPrec= -1
  bestRec= -1
  thr_2_prec={"thr":[], "prec":[]}
  nonClassLabel= int(np.min(allLabels))
  for i in indices:
    tmpThr= allScoresUnique[i]
    binaryScores= np.array([1.0 if score>tmpThr else nonClassLabel for score in allScores])

    currentPrec= precision_score(allLabels, binaryScores)
    thr_2_prec["thr"].append(tmpThr)
    thr_2_prec["prec"].append(currentPrec)
    if abs(currentPrec- targetPrec) < abs(bestPrec- targetPrec):
      bestMcc= matthews_corrcoef(allLabels, binaryScores)
      bestThr= tmpThr
      bestPrec= precision_score(allLabels, binaryScores)
      bestRec= recall_score(allLabels, binaryScores)
  print(thr_2_prec)
  if jsonFname:
    with open(jsonFname,"w") as f:
      json.dump(thr_2_prec, f)
  print(summary.to_string(index=False))
  print("roc auc(mixed/mean): %f/%f Thr: %f mcc: %f prec: %f rec: %f"%(roc_auc,summary.iloc[-1,1], bestThr, bestMcc, bestPrec, bestRec))
  return bestThr

def countIfAboveThr(bestThr, resultsPath):
  totalNum=0
  numAbove=0
  for fname in sorted(os.listdir(resultsPath)):
    #print(fname)
    if ".rec" in fname or ".lig" in fname:
      results= loadResults( resultsPath, fname)
      totalNum+=1
      if results is None: 
        continue
      else:
        scores= np.array(results["prediction"].values)
      if np.sum(scores>=bestThr):
        numAbove+=1
  print("Above thr %d / %d "%(numAbove,totalNum))
  
if __name__=="__main__":
  '''
   python evaluation/getPrecsAt_diffThr.py ~/Tesis/rriPredMethod/data/bench5Data/newCodeData/results/mixed_2/ ~/Tesis/rriPredMethod/data/bench5Data/newCodeData/results/thr_2_prec_b5_mixed_2.json

  '''
  if len(sys.argv)==3:
    resultsPath= os.path.expanduser(sys.argv[1])
    jsonFname= os.path.expanduser(sys.argv[2])
  else:
    raise ValueError("Bad arguments.")
  bestThr= getRecallAtPrec(resultsPath, jsonFname=jsonFname)
  countIfAboveThr(bestThr, resultsPath)
  
