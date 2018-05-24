import sys, time, os
import pandas as pd
import getopt
from .workers.showPymolPath import showPDB_patches_all

RESULTS_THR= 0.3

#PDBs_PATH="/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/b5structures"
#RESULTS_PATH="~/Tesis/rriPredMethod/data/bench5Data/newCodeData/results_xgb/mixed_2/"

PDBs_PATH="/home/rsanchez/Tesis/rriPredMethod/data/develData/pdbFiles/"
RESULTS_PATH="/home/rsanchez/Tesis/rriPredMethod/data/develData/results/mixed_2/"

#PDBs_PATH="/home/rsanchez/Tesis/rriPredMethod/data/partnerSpecificity/trimers/splitPDBs"
#RESULTS_PATH="/home/rsanchez/Tesis/rriPredMethod/data/partnerSpecificity/trimers/results/preds/mixed_2"

#BINDING_CMAPS_PATH="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/newCodeData/computedFeatures/common/contactMapsBinding"
BINDING_CMAPS_PATH=None

def loadDfsCMap(prefix, chainType, cMapsPath= BINDING_CMAPS_PATH):
  resultDF_list= []
  for fname in os.listdir(cMapsPath):
#    print(fname, prefix, os.path.join(cMapsPath,fname))
    if ((chainType=="l" and "_l_" in fname) or (chainType=="r" and "_r_" in fname)) and fname.startswith(prefix):
      df= pd.read_table(os.path.join(cMapsPath,fname),sep='\s+', header='infer', comment="#", 
                        dtype= {"chainIdL":str, "chainIdR":str, "structResIdL":str, "structResIdR":str,
                                      "chainId":str, "structResId":str, "chain":str,  "resIds":str})
      resultDF_list.append( df )
  df= pd.concat(resultDF_list)
  categColNum= list(df.columns).index("categ")
  return fromDfToResIds(df, categColNum, thr=0.5)

def fromDfToResIds(df, colNum, thr):
  interfaceList=[]
  for i in range(df.shape[0]):
    if df.iloc[i, colNum]>=thr:
#      print(df.iloc[i,:])
      interfaceList.append( (df.iloc[i, 0], df.iloc[i, 1]))
  return interfaceList
  
def loadResults(prefix, resultsPath, filterBy="prediction", thr=RESULTS_THR):
  fnameRoot= os.path.join(resultsPath, prefix+".res.tab")
  df_l= pd.read_table(fnameRoot+".lig",sep='\s+', header='infer', comment="#", 
                        dtype= {"chainIdL":str, "chainIdR":str, "structResIdL":str, "structResIdR":str,
                                 "chainId":str, "structResId":str, "chain":str,  "resIds":str})
  df_r= pd.read_table(fnameRoot+".rec",sep='\s+', header='infer', comment="#", 
                        dtype= {"chainIdL":str, "chainIdR":str, "structResIdL":str, "structResIdR":str,
                                 "chainId":str, "structResId":str, "chain":str,  "resIds":str})
  pickColNum= list(df_l.columns).index(filterBy)
  return fromDfToResIds(df_l, pickColNum, thr), fromDfToResIds(df_r, pickColNum, thr)
  
def showOneComplex(prefix, thr= RESULTS_THR, resultsPath= os.path.expanduser(RESULTS_PATH), 
                    cMapsPath= BINDING_CMAPS_PATH, pdbsPath= os.path.expanduser(PDBs_PATH)):
  res_pred_l, res_pred_r= loadResults(prefix, resultsPath, filterBy="prediction", thr=thr)
  if not cMapsPath is None:
    res_true_l= loadDfsCMap(prefix, "l", cMapsPath= cMapsPath)
    res_true_r= loadDfsCMap(prefix, "r", cMapsPath= cMapsPath)
  else:
    res_true_l, res_true_r= loadResults(prefix, resultsPath, filterBy="categ", thr=thr)
  assert len(res_true_l)>0 and len( res_true_r)>0, "Error, there is no residue at interface"
  showPDB_patches_all(pdbsPath, prefix, res_pred_l, res_pred_r, res_true_l, res_true_r)
  

if __name__=="__main__":
    
  options, remainder= getopt.getopt(sys.argv[1:], 'p:i:r:t',
                                     ['prefix=',
                                      'inputFiles=',
                                      'resultsFiles=',
                                      'threshold='
                                     ])
  prefix=None
  thr=RESULTS_THR
  inputFiles= PDBs_PATH
  resultsFiles=RESULTS_PATH
  for opt, arg in options:
    if opt in   ('-i', '--inputFiles'):
      inputFiles = os.path.abspath(os.path.expanduser(arg))
    elif opt in ('-r', '--resultsFiles'):
      resultsFiles = os.path.abspath(os.path.expanduser(arg))
    elif opt in ('-t', '--threshold'):
      print(opt, arg)
      thr = float(arg)
    elif opt in ('-p', '--prefix'):
      prefix = arg
    else:
      raise ValueError("Bad arguments")
  if prefix is None:
    raise ValueError("No prefix provided")
  showOneComplex(prefix, thr, pdbsPath=inputFiles, resultsPath=resultsFiles)
  
