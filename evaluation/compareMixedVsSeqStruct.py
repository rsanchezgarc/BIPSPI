import sys, os
import pandas as pd
from joblib import Parallel
from joblib import delayed
from sklearn.metrics import roc_auc_score

N_JOBS=8

def getPrefixes(resultsPath):
  prefixes= set([])
  for fname in os.listdir( os.path.expanduser(resultsPath) ):
    prefix= fname.split("#")[0].split(".")[0]
    prefixes.add(prefix)
  return prefixes


def loadResult(fname):
  df = pd.read_table(fname, comment="#", sep="\s+", dtype={"resId":str, "chainId":str})
  return df

def getAUC(fname):
  df= loadResult(fname)
  roc_auc= roc_auc_score(df["categ"], df["prediction"])
  return roc_auc


def processOnePrefix(prefix, ligOrRec, seqResultsPath, mixedResultsPath, structResultsPath):
  fnameSeq= os.path.join(seqResultsPath, prefix+".res.tab."+ligOrRec+".gz")
  fnameStruct= os.path.join(structResultsPath, prefix+".res.tab."+ligOrRec+".gz")
  fnameMixed= os.path.join(mixedResultsPath, prefix+"#s"+ligOrRec[0]+".res.tab."+ligOrRec+".gz")

  try:
    rocSeq= getAUC(fnameSeq)
    rocStruct = getAUC(fnameStruct)
    rocMixed= getAUC(fnameMixed)

    out= abs(rocMixed-rocStruct), rocSeq, rocMixed, rocStruct, ligOrRec, prefix
    if N_JOBS<2: print( out )

    return out

  except ValueError:
    return None

seqResultsPath, mixedResultsPath, structResultsPath= sys.argv[1:]

# seqResultsPath=   "~/Tesis/rriPredMethod/data/bench_joan_2411/bipspi_v2/bench5_and_bench_2411/results/seq_2/"
# mixedResultsPath= "~/Tesis/rriPredMethod/data/bench_joan_2411/bipspi_v2/bench5_and_bench_2411/results/mixed_2/"
# structResultsPath="~/Tesis/rriPredMethod/data/bench_joan_2411/bipspi_v2/bench5_and_bench_2411/results/struct_2/"



seqPrefixes= getPrefixes(seqResultsPath)
mixedPrefixes= getPrefixes(mixedResultsPath)
structPrefixes= getPrefixes(structResultsPath)

sharedPrefixes= seqPrefixes.union(mixedPrefixes).union(structPrefixes)

# sharedPrefixes=list(sharedPrefixes)[:10]

results=Parallel(n_jobs=N_JOBS)( delayed(processOnePrefix)(prefix, ligOrRec, seqResultsPath, mixedResultsPath, structResultsPath) for prefix in sharedPrefixes
                                                                                                                for ligOrRec in ["lig", "rec"] )

results= [elem for elem in results if elem is not None ]
results= pd.DataFrame(results)
results.columns= [ "deltaRoc", "rocSeq", "rocMixed", "rocStruct", "ligOrRec", "prefix"]
results.sort_values(by= ["deltaRoc"], ascending=True, inplace=True)

results = results[results["rocMixed"]>0.7]

print(results)
# results.to_csv("seqVsMixVsStruct.csv", sep="\t", index=False)


#fname="seqVsMixVsStruct.csv"
#x=  pd.read_csv(fname, header="infer", sep="\t")
#x[ (0.8<x.rocStruct) & (x.rocStruct<0.95) & (x.rocSeq<0.8) & (x.rocSeq< x.rocMixed)].head()
