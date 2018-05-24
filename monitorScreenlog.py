import sys, os
import pandas as pd

pd.set_option('precision', 4)

BENCH3_PATH= os.path.expanduser("~/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/")
EVAL_JUST_BENCH3= False

def main(fname):
  if EVAL_JUST_BENCH3:
    bench3Prefixes= set([ elem.split(".")[0].split("_")[0] for elem in os.listdir(BENCH3_PATH)])
  pdbsScoresDict={}
  nextLineIsGood=False
  colNames= None
  with open(fname) as f:
    for line in f:
      if line.startswith("pdb  auc_pair  "):
        colNames= line.split()
        nextLineIsGood=True
        continue
      if nextLineIsGood:
        dataLine= line.split()
        pdbId= dataLine[0].split(".")[0]
        if EVAL_JUST_BENCH3 and not pdbId in bench3Prefixes: continue
        dataLine= [dataLine[0]]+ [float(elem) for elem in dataLine[1:]]
        pdbsScoresDict[dataLine[0]]= dataLine
        nextLineIsGood=False
  resDf= pd.DataFrame.from_records(pdbsScoresDict.values())
  resDf.columns= colNames
  resDf.sort_values("pdb", inplace=True)

  means= resDf.mean(axis=0)
  resDf= resDf.append( resDf.ix[resDf.shape[0]-1,:],ignore_index=True )
  resDf.ix[resDf.shape[0]-1,0]=  "mean"
  resDf.ix[resDf.shape[0]-1, 1:]=  means
  print(resDf.to_string(index=False))
  
if __name__=="__main__":
  if len(sys.argv)==2:
    fname=sys.argv[1]
  else:
    fname= "screenlog.0"
  main(fname)
