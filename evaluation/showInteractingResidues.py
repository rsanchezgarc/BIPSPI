import sys, time, os
import pandas as pd
import getopt
from .workers.showPymolPath import showPDB_interactions, showPDB_patches_all

RESULTS_THR= 0.95

#PDBs_PATH="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3"
##RESULTS_PATH="~/Tesis/rriPredMethod/data/ppdockingBenchData/newCodeData/results/struct_2"
##RESULTS_PATH="~/Tesis/rriPredMethod/data/ppdockingBenchData/newCodeData/computedFeatures/common/contactMaps_backup" 
#RESULTS_PATH="~/Tesis/rriPredMethod/data/ppdockingBenchData/newCodeData/computedFeatures/common/contactMaps" 

PDBs_PATH="/home/rsanchez/Tesis/rriPredMethod/data/sppider/allPdbsDownloaded"
#RESULTS_PATH="~/Tesis/rriPredMethod/data/sppider/computedFeatures/common/contactMaps" 
RESULTS_PATH="~/Tesis/rriPredMethod/data/sppider/results/struct"

#FILTER_BY="categ"  # "prediction"
FILTER_BY="prediction"
DISPLAY_PAIRS= False # If True edges will be shown, otherwise binding sites

def fromDfToResIds(df, colNum, thr):
  interfaceList=[]
  values= df.iloc[:, colNum].values
  ids= df.iloc[:, :6].values
  for i in range(df.shape[0]):
    if values[i]>=thr:
      interfaceList.append( ((ids[i, 0], ids[i, 1], ids[i, 2]),(ids[i, 3], ids[i, 4], ids[i, 5])))
#  print(len(interfaceList))
#  raw_input("enter")  
  return interfaceList
  
def loadResults(prefix, resultsPath, filterBy=FILTER_BY, thr=RESULTS_THR):
  fname= os.path.join(resultsPath, prefix+".res.tab")
  if not os.path.isfile(fname):
    fname= os.path.join(resultsPath, prefix+".cMap.tab")  
  df_pair= pd.read_table(fname,sep='\s+', header='infer', comment="#", 
                        dtype= {"chainIdL":str, "chainIdR":str, "structResIdL":str, "structResIdR":str,
                                 "chainId":str, "structResId":str, "chain":str,  "resIds":str})
  pickColNum= list(df_pair.columns).index(filterBy)

  return fromDfToResIds(df_pair, pickColNum, thr)
  
def showOneComplex(prefix, resultsPath= os.path.expanduser(RESULTS_PATH), pdbsPath= PDBs_PATH):
  l_r_pairsres= loadResults(prefix, resultsPath)
  print(l_r_pairsres, len(l_r_pairsres))
  showPDB_interactions(pdbsPath, prefix, l_r_pairsres, showPairs=DISPLAY_PAIRS )
  
if __name__=="__main__":


  options, remainder= getopt.getopt(sys.argv[1:], 'p:i:r:',
                                     ['prefix=',
                                      'inputFiles=',
                                      'resultsFiles='
                                     ])
  prefix=None
  inputFiles= PDBs_PATH
  resultsFiles=RESULTS_PATH
  for opt, arg in options:
    if opt in   ('-i', '--inputFiles'):
      inputFiles = os.path.abspath(os.path.expanduser(arg))
    elif opt in ('-r', '--resultsFiles'):
      resultsFiles = os.path.abspath(os.path.expanduser(arg))
    elif opt in ('-p', '--prefix'):
      prefix = arg
    else:
      raise ValueError("Bad arguments")
      
  if prefix is None:
    print(options)
    raise ValueError("No prefix provided")
  showOneComplex(prefix, pdbsPath=inputFiles, resultsPath=resultsFiles)

  
