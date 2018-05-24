import sys, os
from Bio.PDB.NeighborSearch import NeighborSearch

def getInterface(struct, useCA=True, res2res_dist=8):
  if useCA:
    atomList = [atom for atom in struct[0].get_atoms() if atom.name.startswith("CA")]
  else:
    atomList = [atom for atom in struct[0].get_atoms() if not atom.name.startswith("H")]
  chains= struct[0].child_list
  searcher= NeighborSearch(atomList)
  allNeigs= searcher.search_all(res2res_dist, level= "R")
  residuesBindingSitePerChain={chain.get_id():{"bindingSite":[]} for chain in chains}
  for res1,res2 in allNeigs:
    pdbId1, modelId1, chainId1, resId1 = res1.get_full_id()
    pdbId2, modelId2, chainId2, resId2 = res2.get_full_id()
    if chainId1 !=chainId2:
      residuesBindingSitePerChain[chainId1]["bindingSite"].append( res1.get_id())
      residuesBindingSitePerChain[chainId2]["bindingSite"].append( res2.get_id())

  return residuesBindingSitePerChain
      
      
if __name__ =="__main__":
  from Bio.PDB.PDBParser import PDBParser
  struct= PDBParser().get_structure(sys.argv[1], os.path.expanduser(sys.argv[1]))
  residues= getInterface(struct)
  print(residues)
