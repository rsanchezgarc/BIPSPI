import sys, os
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, average_precision_score, auc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib import rcParams
from scipy import interp

fontpath = '/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf'
prop = font_manager.FontProperties(fname = fontpath)
rcParams['font.family'] = prop.get_name()
rcParams['font.serif'] = ['Times New Roman']
rcParams['text.usetex'] = True


EVAL_PAIRS_SCORE=False
AVERAGE_PER_COMPLEX=False

SAVE_FIG_NAME=None
SAVE_FIG_NAME="ROC-and-precision-recall_DImS_struct_differentScoring_BIPSPI.png"
#SAVE_FIG_NAME="ROC-and-precision-recall_DBv5_Pairs.png"

def loadResults( resultsPath, prefix, evalPairsScores):
  if EVAL_PAIRS_SCORE:
    fname= os.path.join(resultsPath, prefix+".res.tab")
    scoresDf= pd.read_table(fname, comment="#", sep="\s+", dtype={"structResIdL":str, 
    "chainIdL":str, "structResIdR":str, "chainIdR":str})
  else:
    fname= os.path.join(resultsPath, prefix+".res.tab.rec")  
    scoresDf1= pd.read_table(fname, comment="#", sep="\s+", dtype={"resId":str, "chainId":str})
    fname= os.path.join(resultsPath, prefix+".res.tab.lig")  
    scoresDf2= pd.read_table(fname, comment="#", sep="\s+", dtype={"resId":str, "chainId":str})
    scoresDf= pd.concat( [scoresDf1, scoresDf2])
  return list(scoresDf["categ"]), list(scoresDf["prediction"])
    
def plotSeveralPathsROC(dataSetName, resultsPath_list):
  rocCurves= []
  precRecallCurves=[]
  for resultsPath in resultsPath_list:
    resultsPath, name= resultsPath.split(":")
    resultsPath= os.path.expanduser(resultsPath)
    print("%s %s"%(name, resultsPath) )
    allLabels=[]
    allScores=[]
    allRocAucs=[]
    allPRAucs=[]
    allAvPR=[]
    mean_tpr= 0.0 #will be numpy array
    mean_fpr = np.linspace(0, 1, 100)  
    reversed_mean_precision = 0.0
    mean_recall = np.linspace(0, 1, 100)
    nComplexes=0  #will be numpy array
    prefixes= sorted(set([ fname.split(".")[0] for fname in os.listdir(resultsPath)]))
    for prefix in sorted(prefixes):
      print(prefix)
      results= loadResults( resultsPath, prefix, EVAL_PAIRS_SCORE)
      if results is None: continue
      labels, scores= results
      nComplexes+=1        
      if AVERAGE_PER_COMPLEX:
        fpr, tpr, __= roc_curve(labels, scores)
        mean_tpr+= interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        roc_auc_val= roc_auc_score(labels, scores)
        prec,rec, __= precision_recall_curve(labels, scores)
        reversed_recall = np.fliplr([rec])[0]
        reversed_precision = np.fliplr([prec])[0]
        reversed_mean_precision += interp(mean_recall, reversed_recall, reversed_precision)
        reversed_mean_precision[0] = 0.0          
        average_precision= average_precision_score(labels, scores)
        auc_precision_recall= auc(rec, prec)
        allRocAucs.append( roc_auc_val)
        allPRAucs.append(auc_precision_recall)
        allAvPR.append(average_precision) 

      allScores+= scores
      allLabels+= labels
        
    if AVERAGE_PER_COMPLEX:
      mean_tpr /= nComplexes
      tpr= mean_tpr
      fpr= mean_fpr
      reversed_mean_precision /= nComplexes
      reversed_mean_precision[0] = 1
      prec = reversed_mean_precision
      rec= mean_recall   
      roc_auc_val= np.mean( allRocAucs)
      auc_precision_recall= np.mean( allPRAucs)
      average_precision= np.mean( allAvPR)
    else:      
      fpr, tpr, __= roc_curve(allLabels, allScores)
      roc_auc_val= roc_auc_score(allLabels, allScores)
      prec,rec, __= precision_recall_curve(allLabels, allScores)
      prec,rec = zip(* [ (p,r) for p,r in zip(prec,rec) if r>0.01] )
      average_precision= average_precision_score(allLabels, allScores)
      auc_precision_recall= auc(rec, prec)
      
    positiveFraction=( np.sum(np.array(allLabels)==1)/float(len(allLabels)))      
    rocCurves.append(( (fpr, tpr), roc_auc_val, name) )
    precRecallCurves.append(( (prec,rec), (average_precision, auc_precision_recall), name) )


#  plt.rc('font', family='font.serif', serif='Times New Roman') 
  plt.rc('xtick', labelsize=8)
  plt.rc('ytick', labelsize=8)
  plt.rc('axes', labelsize=8)
  
  fig, axArray = plt.subplots(1, 2 )
  tprs=[]
  aucs=[]
  predType= "binding site" if not EVAL_PAIRS_SCORE else "residue-residue contact"
  ax=axArray[0]
  for rocCurve, roc_auc, name in rocCurves:
#    print(name)
#    raw_input("enter")
    fpr, tpr= rocCurve
    ax.plot(fpr, tpr, label="%s (auc= %0.4f)"%(name,roc_auc))
    print("%s (%s ROC auc= %0.4f)"%( name,"mean" if AVERAGE_PER_COMPLEX else "",roc_auc))    
  ax.plot([0, 1], [0, 1],'k--', label="baseline")
  ax.set_xlim([-0.01, 1.01])
  ax.set_ylim([-0.01, 1.01])
  ax.set_xlabel('False Positive Rate')
  ax.set_ylabel('True Positive Rate')
  ax.set_title("%s ROC %s"%(dataSetName, predType))
  ax.legend(loc="lower right", prop={'size': 8})

  ax=axArray[1]
  for precRecCur, (average_prec,auc_precision_recall), name in precRecallCurves:
#    print(name)
#    raw_input("enter")
    prec, reca= precRecCur
#    ax.plot(reca, prec, label="%s (auc= %0.4f; avg precision= %0.4f)"%(name,auc_precision_recall, average_prec))
    ax.plot(reca, prec, label="%s (auc= %0.4f)"%(name,auc_precision_recall))
    print("%s (%s PR auc= %0.4f)"%(name, "mean" if AVERAGE_PER_COMPLEX else "", auc_precision_recall))
  ax.plot([0, 1], [positiveFraction, positiveFraction],'k--', label="baseline")      

  ax.set_xlim([-0.01, 1.01])
  ax.set_ylim([-0.01, 1.01])
  ax.set_xlabel('Recall')
  ax.set_ylabel('Precision')
  ax.set_title("%s Precision-Recall %s"%(dataSetName, predType))
  ax.legend(loc="best", prop={'size': 8})
  fig.set_size_inches(17.8/2.54,8.6/2.54)
  if SAVE_FIG_NAME:  
    fig.savefig(SAVE_FIG_NAME, format="png", dpi=350)
  plt.show()


if __name__=="__main__":
  '''
  example:
python evaluation/getRocCurvesSeveralResults.py DBv5 ~/Tesis/rriPredMethod/data/bench5Data/newCodeData/results_xgb/mixed:struct-1step ~/Tesis/rriPredMethod/data/bench5Data/newCodeData/results_xgb/mixed_2:struct-2steps ~/Tesis/rriPredMethod/data/bench5Data/newCodeData/results_xgb/seq:seq 

python evaluation/getRocCurvesSeveralResults.py DImS ~/Tesis/rriPredMethod/data/joanDimers/results/mixed:struct-1step ~/Tesis/rriPredMethod/data/joanDimers/results/mixed_2:struct-2step ~/Tesis/rriPredMethod/data/joanDimers/results/seq/:seq

python evaluation/getRocCurvesSeveralResults.py DBv3 ~/Tesis/rriPredMethod/data/bench3Data/newCodeData/results/xgb_results/mixed/:BIPSPI-struct-1-step ~/Tesis/rriPredMethod/data/bench3Data/newCodeData/results/xgb_results/mixed_2:BIPSPI-struct-2-step  ~/Tesis/rriPredMethod/data/bench3Data/newCodeData/results/xgb_results/seq/:BIPSPI-seq  ~/Tesis/rriPredMethod/data/pairPred/b3/b3PairPredBindingSiteScores_max/:PAIRpred-d ~/Tesis/rriPredMethod/data/pairPred/b3/b3PairPredBindingSiteScores_ourFun/:PAIRpred-p

  '''
  if len(sys.argv)<3: raise ValueError("Error, bad number of arguments")
  dataSetName= sys.argv[1]
  resultsPath_list= sys.argv[2:]
  plotSeveralPathsROC(dataSetName, resultsPath_list)
  
