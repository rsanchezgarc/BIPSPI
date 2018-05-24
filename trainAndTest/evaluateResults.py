from __future__ import print_function
import os
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, matthews_corrcoef
from collections import Counter
from Config import Configuration
pd.set_option('precision', 4)

EVAL_PAIRS_AT= [50,100,500]

PSAIA_PATH= os.path.join(Configuration().computedFeatsRootDir,"structStep/PSAIA/procPSAIA")

  
def computeAUC(testLabels, predictions):
  '''
    returns ROC's auc or 0.5 if it was not possible to compute it
    @param testLabels: int[]. List of labels (-1 for negative class and 1 for positive class)
    @param testLabels: float[]. List of scores predicted
    @return auc_score: float
  '''
  try:
    return roc_auc_score(testLabels, predictions)
  except ValueError:
    return 0.5
        
def evaluatePairs(prefix, resDf):
  '''
    Computes performance evaluation at pairs level.
    gets auc_score and precision and recall at EVAL_PAIRS_AT thresholds (the n==EVAL_PAIRS_AT[i] highest 
    predictions will be regarded as positive and all others as negatives).

    @param prefix: str. An id for the complex
    @param resDf: pandas.DataFrame. A DataFrame that represents the pairs predictions. It has the following columns:
                          chainIdL structResIdL resNameL chainIdR structResIdR resNameR categ prediction
      
    @return summary: pandas.DataFrame. A DataFrame that has the following columns:
                          pdb  auc_pair  prec_50  reca_50  prec_100  reca_100  prec_500  reca_500 
                     This DataFrame will have just one row. All evaluation metrics are computing at pair level (rows)
  '''
  
  meanPairsScore= resDf["prediction"]
  testLabels= resDf["categ"]
  auc_val = computeAUC( testLabels, meanPairsScore)
  #compute precision and recall at EVAL_PAIRS_AT highest score pair
  precisionAt=[]
  recallAt=[]
  probability_sorted_indexes = meanPairsScore.argsort(axis=0)
  probability_sorted_indexes = probability_sorted_indexes[::-1]
  for evalPoint in EVAL_PAIRS_AT:
    try:
      label_predictions= np.ones(meanPairsScore.shape[0])* np.min( testLabels)
      label_predictions[probability_sorted_indexes[0 : evalPoint]]= np.repeat(1, evalPoint)

      precisionAt.append( precision_score(testLabels[probability_sorted_indexes[0 : evalPoint]],
                                        label_predictions[probability_sorted_indexes[0 : evalPoint]]))
      recallAt.append( recall_score(testLabels, label_predictions))
    except IndexError:
      precisionAt.append( 0.0)
      recallAt.append( 0.0)
  summary= pd.DataFrame({"pdb":[prefix]})
  summary["auc_pair"]= [auc_val] 
  for evalPoint, precisionAt, recallAt in zip(EVAL_PAIRS_AT,precisionAt, recallAt):
    summary["prec_%d"%evalPoint]= [precisionAt]
    summary["reca_%d"%evalPoint]= [recallAt]
  return summary
  
def evaluateBindingSite(prefix, res_pos, res_all, scoresDict, chainType="l"):
  '''
    Computes performance evaluation at interface level.
    gets auc_score and matthews_corrcoef, precision and recall at a threshold computed 
    with getThr()

    @param prefix: str. An id for the complex
    @param res_pos: set(str[])). A set of resIds that are part of the interface. E.x. set(["A_32","*_13"])
    @param res_all: set(str[])). A set of all resIds that are part of the protein. E.x. set(["A_32","*_13"])
    @param scoresDict:{str:float} {ChainId_resId: interface_score}. The scores for each of the amino acids
    @param chainType: str. "l" if residues are from ligand or "r" if they are from receptor
    @return summary: pandas.DataFrame. A DataFrame that has the following columns:
                          auc_%s prec_%s  reca_%s mcc_%s where %s is L or R
                     This DataFrame will have just one row. All evaluation metrics are computing at binding site level
  '''

  scoresList_all, labelsList_all, idsList_all= getClassPredsOnChain(scoresDict, res_all, res_pos)
  try:
    accessSet,nonAcessSet= loadAccesibility(prefix,chainType)
    scoreThr= getThr(scoresList_all, idsList_all, accessSet )
  except OSError:
    scoreThr= 0.04

#  scoreThr= 0.050700
    
  (auc, precision, recall, mcc), __= getInterfaceStatis(scoresList_all, labelsList_all, idsList_all, thr= scoreThr)
  summary= pd.DataFrame({"auc_%s"%chainType:[auc] })
  summary["prec_%s"%chainType]= [precision]
  summary["reca_%s"%chainType]= [recall]
  summary["mcc_%s"%chainType]= [mcc]
  return summary
  
  
def evaluateBothBindingSites(prefix, res_posL, res_allL, scoresDictL, res_posR, res_allR, scoresDictR, scoreThr= 2.39698):
  '''
    Computes performance evaluation at interface level.
    gets auc_score and matthews_corrcoef, precision and recall at a threshold computed 
    with getThr()

    @param prefix: str. An id for the complex
    @param res_posL: set(str[])). A set of resIds that are part of the ligand interface. E.x. set(["A_32","*_13"])
    @param res_allL: set(str[])). A set of all resIds that are part of the ligand. E.x. set(["A_32","*_13"])
    @param scoresDictL:{str:float} {ChainId_resId: interface_score}. The scores for each of the amino acids of the ligand
    @param  scoreThr: float.  scoreThr= 0.050700 for mixed_2 best mcc.  scoreThr= 0.036700 for seq best mcc 
    @return summary: pandas.DataFrame. A DataFrame that has the following columns:
                          auc_bs prec_bs reca_bs mcc_bs
                     This DataFrame will have just one row. All evaluation metrics are computing at binding site level
  '''

  scoresList_allL, labelsList_allL, idsList_allL= getClassPredsOnChain(scoresDictL, res_allL, res_posL)
  scoresList_allR, labelsList_allR, idsList_allR= getClassPredsOnChain(scoresDictR, res_allR, res_posR)
  scoresList_all= scoresList_allL+ scoresList_allR
  labelsList_all= labelsList_allL+ labelsList_allR
  idsList_all= idsList_allL+ idsList_allR
  (auc, precision, recall, mcc), __= getInterfaceStatis(scoresList_all, labelsList_all, idsList_all, thr= scoreThr)
  summary= pd.DataFrame({"auc_bs":[auc] })
  summary["prec_bs"]= [precision]
  summary["reca_bs"]= [recall]
  summary["mcc_bs"]= [mcc]
  return summary
  
def getThr(scoresList_all, idsList_all, accessSet):
  '''
    Gets a threshold for binding site score binarization. It employs the  Howook Hwang et al.
    Protein Sci 2016 Jan 25(1) 159-165. criterium
      n=6.1N**0.3 where n is the expected number of interface residues and N is the number of
      accesible residues
      
    @param scoresList_all: float[]. A list of all the binding-site scores of residues that belongs to 
                                    the receptor or the ligand.
    @param idsList_all: str[]. A set of resIds that are part of the receptor or the ligand. 
                                     E.x. set(["A_32","*_13"])
                                     
    @param accessSet: set(str[])). A set of resIds (belonging to receptor or ligand) that are accesible
                                     E.x. set(["A_32","*_13"])                                     

    @return thr: float. The score value such that residues with greater score than thr will be considered as positive
                        and residues with less score that thr will be considered as negative.
  '''
  nAcces= len(accessSet)
  expectedNum= int(np.round(6.1*(nAcces**0.3)))
  scores_ids_list= sorted(zip(scoresList_all,idsList_all), reverse=True)
  if len(scores_ids_list)<= expectedNum:
    thr= scores_ids_list[-1][0]
  else:
    thr= scores_ids_list[expectedNum][0]
  return thr



def getInterfaceStatis(scoresOrig, labelsOrig, resIdsOrig, resIdsSubset=None, thr=0.03):
  '''
    Computes statistics (auc, precision, recall, mcc) for interface prediction evaluation at a given threshold

    @param scoresOrig: float[]. A list of the binding-site scores of all residues 
    @param labelsOrig: int[].   A list of the labels of all residues. (Same order than scoresOrig)
    @param resIdsOrig: str[].   A list of resIds of all residues. (Same order than scoresOrig). E.x. ["A_32", "*_13"]
    
    @param resIdsSubset: set(str[]). A set of resIds that are wanted to evaluate (resIdsSubset is a subset of 
                                     resIdsOrig). If None, all residues in resIdsOrig will be evaluated.
    @param thr: float. A threshold to decide whether or not an amino acid belongs to the binding-site. precision, 
                       recall and mcc will be calculated using this threshold.
    
    @return (auc, precision, recall, mcc), ids_list.
              auc: float.
              precision: float.
              recall: float.
              mcc: float.
              ids_list: str[]. A list of resIds that, accordint to the threshold belongs to the interface (and are
                               included in resIdsSubset if it is not None)
  '''

  
  scores= []
  scoresBin= []
  labels= []
  ids_list= []
  for  score, label, resId  in zip(scoresOrig, labelsOrig, resIdsOrig):
    if score>= thr  and ( resIdsSubset is None or (resId in resIdsSubset)):
      scores.append(score)
      labels.append(label)
      scoresBin.append(1)
      ids_list.append(resId)
    else:
      scores.append(score)
      labels.append(label)
      scoresBin.append(0)
  try:
    auc=roc_auc_score(labels, scores)
  except ValueError:
    auc= np.nan
  precision=precision_score(labels, scoresBin)
  recall=recall_score(labels, scoresBin)
  mcc= matthews_corrcoef(labels, scoresBin)
#  print("auc/prec/recall/mcc: %.4f %.4f %.4f %.4f"%(auc, precision, recall, mcc))
  return (auc, precision, recall, mcc), ids_list
    
def getClassPredsOnChain( predIdsDict, allIds, posIds):
  '''
    Extracts from the inputs a list of scores, a list of labels and a list of resIds, keeping the
    same order

    @param predIdsDict: {resId:score } A dict with resIds as keys and binding-site scores as values
    @param allIds: set(str[]).         A set of resIds of all residues. E.x. set(["A_32", "*_13"])
    @param posIds: set(str[]).         A set of resIds of  residues that belongs to the interface. 
                                       E.x. set(["A_32", "*_13"])
        
    @return (scores, labels, resIds)
            scores. float[]
            labels. int[]
            resIds. str[]
  '''
  
  labels= np.zeros(len(allIds) ).astype(np.int64)
  scores= np.zeros( labels.shape)
  resIds= sorted(allIds)
  i=0
  for res in resIds:
    if res in posIds:
      labels[i]=1
    if res in predIdsDict:
      scores[i]= predIdsDict[res]
    else:
      scores[i]=0
    i+=1  
  return list(scores), list(labels), resIds
  
def loadAccesibility(pdbId, chainType="l", rasaThr=10.0):
  '''
    Loads psaia files for por a given pdbId and returns a set of accesible
    resIds and non-accesible resIds.

    @param pdbId: str. The identifier for pdb file
    @param chainType: str. "l" for ligan and "r" for receptor
    @param rasaThr: float. A threshold of relative asa to decide whether or not a residue is accesible or not
        
    @return (accesibleSet, nonAccesibleSet)
        accesibleSet: set(str[]). Set of resIds of residues that are accesible according to PSAIA and the threshold
        nonAccesibleSet: set(str[]). Set of resIds of residues that are non-accesible according to PSAIA and the threshold
  '''
  accesibleSet=set([])
  nonAccesibleSet=set([])
  for fname in os.listdir(PSAIA_PATH):
    if fname.startswith(pdbId+"_"+chainType) and fname.endswith(".psaia.tab"):
      with open(os.path.join(PSAIA_PATH, fname)) as f_:
        f_.readline()
        for line in f_:
          lineArray= line.split()
          chainId= lineArray[0]
          resId= lineArray[1]
          res_full_id= chainId+"_"+resId
          if float(lineArray[8])> rasaThr:
            accesibleSet.add(res_full_id)
          else:
            nonAccesibleSet.add(res_full_id)
  return accesibleSet, nonAccesibleSet
  
if __name__=="__main__":
  acces,nonAcess= loadAccesibility("1A2K","r")
  print(len(acces), len(nonAcess))
