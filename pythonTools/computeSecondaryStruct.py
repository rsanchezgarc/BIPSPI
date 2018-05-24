import sys, os
from Bio.PDB.DSSP import DSSP

DSSP_SYMBOLS= "HBEGITS-"

def getDSSP(struct, fname, dsspPath=os.path.expanduser("~/Tesis/rriPredMethod/dependencies/bioinformaticTools/dssp/mkdssp")):
  dssp = DSSP(struct[0], fname, dssp=dsspPath)
  chains= struct[0].child_list
  dsspDict= { chain.get_id():{symbol:[] for symbol in DSSP_SYMBOLS} for chain in chains}
  for chainId, resId in dssp.keys():
    secStruct= dssp[(chainId, resId)][2]
    dsspDict[chainId][secStruct].append( resId)

  return dsspDict
    
if __name__ =="__main__":
  from Bio.PDB.PDBParser import PDBParser
  struct= PDBParser().get_structure(sys.argv[1], os.path.expanduser(sys.argv[1]))
  dsspDict= getDSSP(struct, os.path.expanduser(sys.argv[1]))
  print(dsspDict)
