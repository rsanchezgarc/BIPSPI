import sys, os

import re
from sklearn.metrics import roc_auc_score, recall_score, precision_score, roc_curve
import pandas as pd
import numpy as np

from evaluation.evaluateScoresList import evaluateScoresLists
from computeFeatures.featuresComputer import getPrefix

SKIP_LOWER=True

SEQ_ONLY_IN_MIXED=False
STRUCT_ONLY_IN_MIXED=False

DO_ROC_CURVE= False

DO_SEQ_AVERAGING=False # True

#BINDING_CMAPS_PATH="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/newCodeData/computedFeatures/common/contactMapsBinding"
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
  if fnameResults.replace(".gz","").endswith(".lig"):
    chainType="l"
  else:
    chainType="r"

  mixedFlagForSeqOrNone= re.match("[a-zA-Z0-9]+\#s([lr])", prefix)
  if mixedFlagForSeqOrNone:
    mixedFlagForSeqOrNone= mixedFlagForSeqOrNone.group(1)
    if SEQ_ONLY_IN_MIXED:
      if chainType!=mixedFlagForSeqOrNone:
        return None
    if STRUCT_ONLY_IN_MIXED:
      if chainType == mixedFlagForSeqOrNone:
        return None

  scoresDf= pd.read_table(os.path.join(resultsPath, fnameResults), comment="#", sep="\s+", dtype={"resId":str, "chainId":str})
  if not cMapsPath is None:
    newCmapSet=set([])
    for fname in os.listdir(cMapsPath):
      #print(fnameResults, fname, prefix, os.path.join(cMapsPath,fname))
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
#  EVAL_PAIRS_AT= [ 2**3, 2**4]
  EVAL_PAIRS_AT= [ 0.01, 0.05, 0.1]
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
    if evalPoint<1:
      evalPoint= int( scores.shape[0] * evalPoint)
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
  summary["auc_chains"]= [roc_complex] 
  for evalPoint, precisionAt, recallAt in zip(EVAL_PAIRS_AT,precisionAt, recallAt):
    summary["prec_%s"%evalPoint]= [precisionAt]
    summary["reca_%s"%evalPoint]= [recallAt]
    
  rocCurve= None
  if DO_ROC_CURVE:
    fpr, tpr, __= roc_curve(labels, scores)
    rocCurve= (fpr, tpr, roc_complex)
  return summary, rocCurve

def readTxt(fname):
  print("reading ids from %s"%fname)
  ids=[]
  with open(fname) as f:
    for line in f:
      ids.append( getPrefixInMixed(line.split()[0]) )
  return ids

def getPrefixInMixed(fname):
  fname= fname.replace("#sl","")
  fname = fname.replace("#sr","")
  # fname= fname[:4]
  return getPrefix(fname)

def getOptimThr(resultsPath, fileOfRestrictedIds=None, invertSelection=False, useSeqAvera= DO_SEQ_AVERAGING):
  allScores=[]
  allLabels=[]
  perComplexSummaries=[]
  rocCurves= []
  allFnames= sorted( [ fname for fname in os.listdir(resultsPath)
                              if ((".rec" in fname or ".lig" in fname) and not ".json" in fname )  ])
  if invertSelection:
    selectionFunction= lambda fname: not getPrefixInMixed(fname) in selectedIds
  else:
    selectionFunction= lambda fname: getPrefixInMixed(fname) in selectedIds

  if fileOfRestrictedIds:
    selectedIds= set(readTxt(fileOfRestrictedIds))
    allFnames= [fname for fname in allFnames if selectionFunction(fname)  ]


  for fname in allFnames:
    if SKIP_LOWER and fname[:4].islower(): continue
    print(fname)
    results= loadResults( resultsPath, fname)
    if results is None: print("skip"); continue
    if useSeqAvera:
      scores, labels= averageScores(results)
    else:
      scores= list(results["prediction"].values)
      labels= list(results["categ"].values)

    summary, rocCurve= get_single_chain_statistics(fname, labels, scores)
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
python -m evaluation.getBestThr_bindingSite  ~/Tesis/rriPredMethod/data/bench5Data/newCodeData/results/mixed_2/
  '''
  resultsPath= os.path.expanduser("~/Tesis/rriPredMethod/data/bench5Data/newCodeData/results/mixed_2/")
  fileOfRestrictedIds=None
  invertSelection=False
  if len(sys.argv)>=2:
    resultsPath= os.path.expanduser(sys.argv[1])
  if len(sys.argv)>=3:
    fileOfRestrictedIds= os.path.expanduser(sys.argv[2])
  if len(sys.argv)==4:
    invertSelection=True
  getOptimThr(resultsPath, fileOfRestrictedIds= fileOfRestrictedIds, invertSelection= invertSelection)
  
