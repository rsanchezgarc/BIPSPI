import sys, os
from sklearn.metrics import roc_auc_score, recall_score, precision_score, roc_curve
import pandas as pd
import numpy as np
from evaluation.evaluateScoresList import evaluateScoresLists
  
DO_ROC_CURVE= False

# if INVERT_SELECTION==True, evaluate chains that are not in the file
INVERT_SELECTION= False 
# if EVAL_AT_COMPLEX_LEVEL==True, select at complex level, otherwise select at chains levels
EVAL_AT_COMPLEX_LEVEL=True

pathDiffComplex= "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/benchmark5/difficultComplexes.txt"
difficultComplexes= set([])
with open(pathDiffComplex) as f:
  for line in f:
    prefix= line.split("_")[0]
    difficultComplexes.add(prefix)
pathMediumComplex= "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/benchmark5/mediumComplexes.txt"
with open(pathMediumComplex) as f:
  for line in f:
    prefix= line.split("_")[0]
    difficultComplexes.add(prefix)
        
def loadResults( resultsPath, fnameResults, selectChains, invertSelChains):
  prefix= fnameResults.split("_")[0].split(".")[0]
  if fnameResults.endswith(".lig"):
    chainType="l"
  else:
    chainType="r"
    
  scoresDf= pd.read_table(os.path.join(resultsPath, fnameResults), comment="#", sep="\s+", dtype={"resId":str, "chainId":str})      
  if selectChains:
    if invertSelChains:
      scoresDf= scoresDf[~ scoresDf["chainId"].isin(selectChains) ]
    else:
      scoresDf= scoresDf[ scoresDf["chainId"].isin(selectChains) ]
  nPos= np.sum(scoresDf["categ"]==1 )
  if scoresDf.shape[0]==0:
    fraccPos= np.nan
  else:
    fraccPos= nPos/ float(scoresDf.shape[0])
  print(scoresDf.shape)
  return scoresDf, fraccPos
  
def get_single_chain_statistics(prefix, labels, scores):
  EVAL_PAIRS_AT= [ 2**3, 2**4, 2**5]
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
  summary["auc_chains"]= [roc_complex] 
  for evalPoint, precisionAt, recallAt in zip(EVAL_PAIRS_AT,precisionAt, recallAt):
    summary["prec_%d"%evalPoint]= [precisionAt]
    summary["reca_%d"%evalPoint]= [recallAt]
    
  rocCurve= None
  if DO_ROC_CURVE:
    fpr, tpr, __= roc_curve(labels, scores)
    rocCurve= (fpr, tpr, roc_complex)
  return summary, rocCurve
  
def loadChainsToEval(chainsToEvalFname):
  prefixToChains={}
  with open(chainsToEvalFname) as f:
    for line in f:
      lineArray= line.split()
      for elem in lineArray[1:]:
        prefix, chain= elem.split("_")
        if prefix not in prefixToChains:
          prefixToChains[prefix]=[]
        prefixToChains[prefix].append(chain)
  return prefixToChains
  
def getOptimThr(resultsPath, chainsToEvalFname):
  allScores=[]
  allLabels=[]
  perComplexSummaries=[]
  rocCurves= []
  prefixToChains= loadChainsToEval(chainsToEvalFname)
  prefixes= set([])
  nDifficult= 0
  fraccPos_list=[]
  for fname in sorted(os.listdir(resultsPath)):
    print(fname)
    if fname.endswith(".rec") or fname.endswith(".lig"):
      prefix = fname.split(".")[0]
      if not EVAL_AT_COMPLEX_LEVEL:
        selectChains= []
        if prefix in prefixToChains:
          selectChains= prefixToChains[prefix]
          if len(selectChains)==0:
            if not INVERT_SELECTION:
              continue
      else:
        selectChains= None
        if prefix in prefixToChains and INVERT_SELECTION:
          continue
        if not prefix in prefixToChains and not INVERT_SELECTION:
          continue
        
      results, fraccPos= loadResults( resultsPath, fname, selectChains, INVERT_SELECTION)
      if results is None or results.shape[0]==0: print("skip"); continue
      fraccPos_list.append(fraccPos)
      prefixes.add(prefix)
      if prefix in difficultComplexes:
        nDifficult+=1
        difficultComplexes.remove(prefix)
               
      scores= list(results["prediction"].values)
      labels= list(results["categ"].values)

      summary, rocCurve= get_single_chain_statistics(fname, labels, scores)
      if rocCurve: rocCurves.append(rocCurve)
      perComplexSummaries.append(summary)
#      print("%s %f"%(fname,roc_complex))
      allScores+= scores
      allLabels+= labels

  summary= pd.concat(perComplexSummaries, ignore_index=True)
  summary= summary.dropna(axis=0)
  means= summary.mean(axis=0)
  summary= summary.append( summary.ix[summary.shape[0]-1,:],ignore_index=True )
  summary.ix[summary.shape[0]-1,0]=  "mean"
  summary.ix[summary.shape[0]-1, 1:]=  means
  print("%d %d -> %f (difficult/total)"%( nDifficult, len(prefixes), float(nDifficult) / len(prefixes)) )
  print("fraccPos %f"%( np.nanmean(fraccPos_list) ))

  evaluateScoresLists(allLabels, allScores, summary, summary.iloc[-1,1], None if not DO_ROC_CURVE else rocCurves)

if __name__=="__main__":
  '''
python -m evaluation.evaluateChainsInList ~/Tesis/rriPredMethod/data/bench5Data/newCodeData/results_xgb/mixed_2/ ~/Tesis/rriPredMethod/data/partnerSpecificity/equivalentChainsB5.txt 
  '''
  if len(sys.argv)!=3:
    raise ValueError("Error, bad number of arguments.usage:\npython -m evaluation.evaluateChainsInList pathToResults pathToChainsToEvaluate")
  resultsPath= os.path.expanduser(sys.argv[1])
  chainsPath= os.path.expanduser(sys.argv[2])
  getOptimThr(resultsPath, chainsPath)
  
