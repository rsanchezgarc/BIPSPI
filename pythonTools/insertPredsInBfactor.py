import sys, os
import numpy as np
import pandas as pd
from Bio import PDB

try:
  from pythonTools.myPDBParser import myPDBParser as PDBParser
  from computeFeatures.toolManagerGeneric import UNKNOWN_CHAIN_SYMB
except ImportError:
  from Bio.PDB.PDBParser import PDBParser
  UNKNOWN_CHAIN_SYMB="+"

SCALING_RANGE=50
def insertPredsInBfactor(pdbFnameIn, scoresFnameIn, pdfFnameOut, scale=False, min_max_norm=False, verbose=True):
  parser = PDBParser(QUIET=True)
  struct = parser.get_structure(os.path.basename(pdbFnameIn), pdbFnameIn)
  scores= pd.read_table(scoresFnameIn,sep='\s+', header='infer', comment="#", 
                          dtype= {"chainIdL":str, "chainIdR":str, "resIdL":str, "resIdR":str,
                                   "chainId":str, "resId":str,  "resId":str})

  if scale: scores["prediction"]= SCALING_RANGE* scores["prediction"]
  if min_max_norm:
    scores["prediction"]=  ( scores["prediction"] - np.min(scores["prediction"]) ) / ( np.max(scores["prediction"]) - np.min(scores["prediction"])  )

  scoresDict= {}
  for i in range(scores.shape[0]):                
    scoresDict[(scores["chainId"][i], scores["resId"][i])]= scores["prediction"][i]
#  print( sorted([ (key, scoresDict[key]) for key in scoresDict]))    
  for chain in struct[0]:
    chainId=chain.get_id()
    if chainId==" ": chainId=UNKNOWN_CHAIN_SYMB
    for res in chain:
      resId= res.get_id()
      strResId= (str(resId[1])+resId[2]).strip()
      if (chainId, strResId) in scoresDict:
        predVal= scoresDict[(chainId, strResId)]
      else:
        predVal= 0.0
      if verbose: print(chainId, strResId, (chainId, strResId) in scoresDict, predVal)
      for atom in res:
        atom.set_bfactor(predVal)
                                
  writer = PDB.PDBIO()
  writer.set_structure(struct)
  writer.save(pdfFnameOut)
  
if __name__=="__main__":
#  pdbFnameIn="/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/b5structures/1ACB_l_u.pdb"
#  scoresIn="/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/newCodeData/results_xgb/mixed_2/1ACB.res.tab.lig"
#  pdfFnameOut="/home/rsanchez/tmp/b_factor_pdb_trial.pdb"
#  insertPredsInBfactor(pdbFnameIn, scoresIn, pdfFnameOut)
  
  if len(sys.argv)==4:
    pdbFnameIn= sys.argv[1]
    scoresIn= sys.argv[2]
    pdfFnameOut= sys.argv[3]
    insertPredsInBfactor(pdbFnameIn, scoresIn, pdfFnameOut, min_max_norm=True)
  else:
    raise ValueError("Bad number of arguments")  
