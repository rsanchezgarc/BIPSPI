import sys, os
import re
from numpy import array
from joblib import load
import xgboost as xgb
from codifyComplexes.ComplexCodified import ComplexCodified
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

SAVE_FIGURES=False
modelPath="/home/rsanchez/Tesis/rriPredMethod/data/develData/modelsComputed/model.struct_2"
featurePath="/home/rsanchez/Tesis/rriPredMethod/data/develData/codifiedInput/struct_2/sampledInputs/1ACB@00.train.pkl.gz"
if len(sys.argv)>=2:
  modelPath= os.path.abspath(os.path.expanduser(sys.argv[1]))

if len(sys.argv)==3:
  featurePath= os.path.abspath(os.path.expanduser(sys.argv[2]))
if SAVE_FIGURES:
  if modelPath.endswith("seq"):
    FIGNAME= 'Importance of features in seq classifier'
  elif modelPath[-1]=="2":
    FIGNAME= 'Importance of features in 2-steps classifier'
  else:
    FIGNAME= 'Importance of features in 1-step classifier'

else:
  FIGNAME=None
  
def get_xgb_imp(model, feat_names):
  imp_vals = model.feature_importances_
  imp_dict = {feat_names[i]:imp_vals[i] for i in range(len(feat_names))}
  total = array(imp_dict.values()).sum()
  return sorted([(k,v/total) for k,v in imp_dict.items()], key=lambda x: x[1])
  
  
def displayImportance(modelPath, featurePath):
  model = load(modelPath)
  feats = load(featurePath)
  featsNames= list(feats.pairsDirect.columns)[7:]
  featImp= get_xgb_imp(model, featsNames)
  print(featImp)
  #featImp= [ elem for elem in featImp if "Aggr" not in elem[0] and "." not in elem[0] and "_dummy" not in elem[0]]
  print("---------------------------------------\n\n")
  df= pd.DataFrame(featImp)
  df.columns= ["featName", "featImp"]

  selRows=[]
  for i in range(df.shape[0]):
    if "ASA" in df.featName[i]:
      selRows.append(i)
  print(df.iloc[selRows,:])
  print("Most important feature index:%d value: %f"%(df.shape[0]-1, df["featImp"][df.shape[0]-1]))
  print(feats.pairsDirect.shape)

def importancePerGroup(modelPath, featurePath, figName= FIGNAME):
  model = load(modelPath)
  feats = load(featurePath)
  featsNames= list(feats.pairsDirect.columns)[7:]

#  print(list(feats.pairsDirect.columns)[:15])
  
#  print( sorted(set([ feat.split(".")[0] for feat in featsNames])))
  print( sorted(set([ feat for feat in featsNames])))
  print(len(featsNames))
#  raw_input("press enter")

  featImp= get_xgb_imp(model, featsNames)
  df= pd.DataFrame(featImp)
  featsNames= df[0]
  df.columns= ["featName", "featImp"]

  typesOfFeatures= {"Conservation":[".*pssm.*", ".*psfm.*", ".*information.*", ".*al2coScore.*"], 
  "AA Symbol":[".*aaWin.*"],
  "HS exposure":[".*Exposure.*"],
  "2nd structure":[".*2ndStruct.*", ".*score_P.*"], "Accesibility":[".*ASA.*", ".*score_asa.*"], 
  "Hydrophobicity":[".*Hydrophobicity.*"], "Depth index":[".*DPX.*"], "Protrusion index":[".*CX.*"],
  "Prev step score":[".*prediction.*"], "SeqLen": [".*sequence length.*"]}

  
#  typesOfFeatures= {"environment-struct":["^.*Aggr.*"],"environment-seq":["^aaWin.*", ".*psfmWin.*", "^pssmWin\..*", "^informationWin\..*"]}
  
  summary={"feat_type":[], "num_variables":[], "mean_importance":[], "sum_importance":[], "max_importance":[]}
  featsNamesSet= set(featsNames)
  involvedNames=[]
  for featType in sorted(typesOfFeatures):
    involvedRows=[]
    for nameFragment in typesOfFeatures[featType]:
      matchesRows=[ i for i,name in enumerate(featsNames) if re.match(nameFragment, name)]
      matchesNames=[ name for name in featsNames if re.match(nameFragment, name)]
      if len(matchesRows)==0:
        print(nameFragment)
      involvedRows+= matchesRows
      involvedNames+= matchesNames
    importanceVals= df["featImp"][involvedRows]
    if featType=="prevStep":
      print(df.iloc[involvedRows])
    summary["feat_type"].append(featType)
    summary["num_variables"].append( len(involvedRows))
    summary["mean_importance"].append( np.mean(importanceVals))
    summary["sum_importance"].append( np.sum(importanceVals))
    summary["max_importance"].append( np.max(importanceVals))
    
  nonInvolvedNames= featsNamesSet.difference( involvedNames)
  involvedRows=[ i for i,name in enumerate(featsNames) if name in nonInvolvedNames]
#  featType= "single residue/pair"
  typesOfFeatures[featType]=None
  if len(matchesRows)!=0:
    importanceVals= df["featImp"][involvedRows]
  summary["feat_type"].append(featType)
  summary["num_variables"].append( len(involvedRows))
  summary["mean_importance"].append( np.mean(importanceVals))
  summary["sum_importance"].append( np.sum(importanceVals))
  summary["max_importance"].append( np.max(importanceVals))
    
      
  summary= pd.DataFrame(summary)
  summary= summary[["feat_type", "num_variables", "mean_importance", "sum_importance", "max_importance"]]
  summary.columns= ["feat_type", "num_variables", "mean_importance", "sum_importance", "max_importance"]
  summary= summary.sort_values(by="mean_importance")
  summary.dropna(inplace=True)
  print(summary.to_string())
  print(summary.sum(axis=0))
  
  plt.rc('font', family='serif', serif='Times New Roman')
  plt.rc('xtick', labelsize=8)
  plt.rc('ytick', labelsize=8)
  plt.rc('axes', labelsize=8)
  

  fig = plt.figure(figsize=(7.2,(3/5.0)*7.2))
  fig.suptitle(figName, fontsize=10, fontweight='bold')
  #THIS IS WHERE OTHERS IS ADDED
  namesOfTypesOfFeatures= list(summary["feat_type"])+["others"]
  importanceValuesForPie= list(summary["sum_importance"])+[1- np.sum(summary["sum_importance"])]
  
  colorsDict= { name: color for name, color in zip(namesOfTypesOfFeatures, plt.cm.Paired(np.linspace(0., 1., len(namesOfTypesOfFeatures))))}  
  colors=[ colorsDict[name] for name in namesOfTypesOfFeatures]


  gridspec.GridSpec(6,6)
  ax= plt.subplot2grid((6,6), (0,0), colspan=3, rowspan=5)
  ax.set_title('global', fontsize=9, fontweight='bold')

  patches, texts = ax.pie( importanceValuesForPie, startangle=100, colors= colors)
  ax.legend( patches, labels=['%s, %2.2f %%' % (l, s*100.0) 
                    for l, s in zip(namesOfTypesOfFeatures, importanceValuesForPie)] , 
             loc="lower left", bbox_to_anchor=(.0,-.3))
                            
  ax= plt.subplot2grid((6,6), (0,4), colspan=2, rowspan=5)

  ax.set_title('mean', fontsize=9, fontweight='bold')
#  ax.yaxis.tick_right()
  ind= 2*np.arange(1, summary.shape[0]+1)
  ax.bar(ind,summary[ "mean_importance"]*100., color= colors)
  ax.set_xticks(ind)
  ax.set_xticklabels(list(summary["feat_type"]), rotation=75) #50
  ax.spines["top"].set_visible(False)
  ax.spines["right"].set_visible(False)
  ax.set_ylabel('mean importance per variable (%)')

  fig.subplots_adjust(left=.03, bottom=.1, right=.95, top=.90)
  if figName:
    fig.savefig(figName.replace(" ","_")+".png", format="png")
  plt.show()
  return None
  
if __name__=="__main__":

#  displayImportance(modelPath, featurePath)
  importancePerGroup(modelPath, featurePath)

