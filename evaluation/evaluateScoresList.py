from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, matthews_corrcoef, roc_curve
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import interp
'''
This is a helper function used by other components.
'''
def evaluateScoresLists(allLabels, allScores, summaryDf, meanAuc, rocCurves=None):
  roc_auc= roc_auc_score(allLabels, allScores)
  if rocCurves:
    fig, ax = plt.subplots()
    tprs=[]
    aucs=[]
    base_fpr = np.linspace(0, 1, 101)
    for rocCurve in rocCurves:
      fpr, tpr, auc_val= rocCurve
      aucs.append(auc_val)
      ax.plot(fpr, tpr, alpha=0.1)
      tpr = interp(base_fpr, fpr, tpr)
      tpr[0] = 0.0
      tprs.append(tpr)
    tprs = np.array(tprs)
    mean_tprs = tprs.mean(axis=0)
    std = tprs.std(axis=0)
    tprs_upper = np.minimum(mean_tprs + std, 1)
    tprs_lower = mean_tprs - std
    ax.plot(base_fpr, mean_tprs, 'b', label="Mean ROC curve leave-one-out (auc= %0.3f)"%np.mean(aucs))
    ax.fill_between(base_fpr, tprs_lower, tprs_upper, color='grey', alpha=0.2)
    ax.plot([0, 1], [0, 1],'k--')
    ax.set_xlim([-0.01, 1.01])
    ax.set_ylim([-0.01, 1.01])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title("model structure no feedback")
    
    fpr, tpr, __= roc_curve(allLabels, allScores)
    ax.plot(fpr, tpr, 'r', label="ROC curve all scores (auc= %0.3f)"%roc_auc )
    ax.legend(loc="lower right")
    
    fig.savefig("ROC-bindingSite-struct-no-feed.png", format="png", dpi=350)
    plt.show()
    
  allScores= np.array(allScores)
  allScoresUnique= np.unique(np.round(allScores, decimals=4))
  indices= np.argsort(allScoresUnique)[1:-1]
#  plt.hist(allScoresUnique[:-2], bins=50)
#  plt.show()
  bestThr=-1
  bestMcc= 0
  bestPrec= -1
  bestRec= -1
  bestAcc= -1
  bestFpr= -1
  bestSpc= -1
  bestNPV= -1
  nonClassLabel= int(np.min(allLabels))
  print(summaryDf.to_string(index=False))
  for i in indices:
    tmpThr= allScoresUnique[i]
    binaryScores= np.where(allScores<=tmpThr,nonClassLabel,1)
    mcc= matthews_corrcoef(allLabels, binaryScores)
    if bestMcc<mcc:
      bestMcc= mcc
      bestThr= tmpThr
      bestPrec= precision_score(allLabels, binaryScores)
      bestRec= recall_score(allLabels, binaryScores)
      bestAcc= accuracy_score(allLabels, binaryScores)
      fp= sum([1 if  score==1 and label!=1 else 0 for label, score in zip(allLabels,binaryScores)])
      tn= sum([1 if  score!=1 and label!=1 else 0 for label, score in zip(allLabels,binaryScores)])
      fn= sum([1 if  score!=1 and label==1 else 0 for label, score in zip(allLabels,binaryScores)])
      bestFpr= fp/ float( fp+ tn)
      bestSpc= tn/float(tn+fp)
      bestNPV= tn/float( tn+fn)
  print(summaryDf.to_string(index=False))
  print( ("roc auc(mixed/mean): %f/%f Thr: %f mcc: %f precision(ppv): %f recal(tpr): %f acc: %f fpr: %f spc(tnr): %f "
          " NPV: %f")%(
          roc_auc, meanAuc, bestThr, bestMcc, bestPrec, bestRec, bestAcc, bestFpr, bestSpc, bestNPV))
  return bestThr
