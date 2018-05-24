import sys, os
import pandas as pd
from Bio import PDB
try:
  from computeFeatures.structStep.myPDBParser import myPDBParser as PDBParser
except ImportError:
  from Bio.PDB.PDBParser import PDBParser
  
def insertPredsInBfactor(pdbFnameIn, scoresFnameIn, pdfFnameOut):
  parser = PDBParser(QUIET=True)
  struct = parser.get_structure(os.path.basename(pdbFnameIn), pdbFnameIn)
  scores= pd.read_table(scoresFnameIn,sep='\s+', header='infer', comment="#", 
                          dtype= {"chainIdL":str, "chainIdR":str, "structResIdL":str, "structResIdR":str,
                                   "chainId":str, "structResId":str,  "resId":str})
  scoresDict= {}
  for i in range(scores.shape[0]):                
    scoresDict[(scores["chainId"][i], scores["resId"][i])]= scores["prediction"][i]
#  print( sorted([ (key, scoresDict[key]) for key in scoresDict]))    
  for chain in struct[0]:
    chainId=chain.get_id()
    if chainId==" ": chainId="*"
    for res in chain:
      resId= res.get_id()
      strResId= (str(resId[1])+resId[2]).strip()
      print(strResId, (chainId, strResId) in scoresDict)
      if (chainId, strResId) in scoresDict:
        predVal= scoresDict[(chainId, strResId)]
      else:
        predVal= 0.0
#      print(predVal)
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
    insertPredsInBfactor(pdbFnameIn, scoresIn, pdfFnameOut)    
  else:
    raise ValueError("Bad number of arguments")  
