import os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide  import is_aa
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.NeighborSearch import NeighborSearch 
from itertools import combinations

MIN_SEQ_LEN= 30
MAX_SEQ_LEN= 800

parser= PDBParser(QUIET= True)
def splitOnePDB(fname, outPath):

  try:
    s= parser.get_structure(fname, fname)
  except Exception:
    print ("Error loading pdb")
    return 0
  banLenChains=[]    
  try:
    for chain in s[0]:
      badResInChain=0
      for res in  chain.get_list():
        if not is_aa(res,standard=True):
          badResInChain+=1
      chainLen= sum(1 for res in chain if "CA" in res) - badResInChain
      if chainLen < MIN_SEQ_LEN or chainLen > MAX_SEQ_LEN:
        print(chainLen)
        banLenChains.append(chain.get_id())
  except KeyError:
    print ("Not good model")
    return 0  
  for badChainId in banLenChains:
    s[0].detach_child(badChainId)

  receptorChainList= []
  ligandChainList= []
  if len( s[0].get_list())<2:
    print(s)
    print( s[0].get_list())
    print("Not enough good chains")
    return 0
  for chain1 in s[0]:

    tmpReceptorList=[]
    for chain2 in s[0]:
      if chain1!= chain2:
        tmpReceptorList.append(chain2)
    if len(tmpReceptorList)>1 or not tmpReceptorList[0] in ligandChainList:   
      ligandChainList.append(chain1)
      receptorChainList.append(tmpReceptorList)
    
  prefix= os.path.basename(fname).split(".")[0]
  for i, (ligandChain, receptorChains) in enumerate(zip(ligandChainList, receptorChainList)):
    io=PDBIO()
    ligandStruct= Structure(prefix+"ligand")
    ligandStruct.add(Model(0))
    ligandChain.set_parent(ligandStruct[0])
    ligandStruct[0].add(ligandChain)
    io.set_structure(ligandStruct)
    io.save(os.path.join(outPath,prefix+"-"+str(i)+"_l_u.pdb"))

    io=PDBIO()
    receptorStruct= Structure(prefix+"receptor")
    receptorStruct.add(Model(0))
    for receptorChain in receptorChains:
      receptorChain.set_parent(receptorStruct[0])    
      receptorStruct[0].add(receptorChain)
    io.set_structure(receptorStruct)
    io.save(os.path.join(outPath,prefix+"-"+str(i)+"_r_u.pdb"))
    print( "ligand:", ligandChain, "receptor:",receptorChains )
    

def splitAllInDir(dirIn, dirOut):
  for fname in os.listdir(dirIn):
    splitOnePDB(os.path.join(dirIn,fname), dirOut)
if __name__=="__main__":
  inPath= "/home/rsanchez/Tesis/rriPredMethod/data/joanDimers/pdbFiles/rawPDBs"
  outPath= "/home/rsanchez/Tesis/rriPredMethod/data/joanDimers/pdbFiles/pdbsForTrain"
  splitAllInDir(inPath, outPath)
