'''

'''
import os, sys

from Bio import PDB
from string import ascii_uppercase
from Bio.PDB.Polypeptide import is_aa
from computeFeatures.toolManagerGeneric import fromRes2ChainResIdAndName

try:
  from myPDBParser import myPDBParser as PDBParser
except ImportError:
  from pythonTools.myPDBParser import myPDBParser as PDBParser

def get_unpacked_list(self):
  """
  Dirty fix for Bio.PDB bug https://github.com/biopython/biopython/issues/455
  Returns all atoms from the residue,
  in case of disordered, keep only first alt loc and remove the alt-loc tag
  """
  atom_list = self.get_list()
  undisordered_atom_list = []
  for atom in atom_list:
     if atom.is_disordered():
         atom.altloc=" "
         undisordered_atom_list.append(atom)
     else:
         undisordered_atom_list.append(atom)
  return undisordered_atom_list
PDB.Residue.Residue.get_unpacked_list = get_unpacked_list


def combinePDBs(fnamesOrStructs_list, outName, fixResIds= True, override=False, ignoreNoStandard=True):
  '''
  :param. fnamesOrStructs_list: list of fnames to pdb files or Bio.PDB.Structure.Structure
  
  :return Bio.PBB.Structure.Structure
  :return chainsMap: { chainIdNew: (chainIdOld,fnameOld, pdbNum) } pdbNum The number of the pdb in the fnamesOrStructs_list
  :return resMap: { "b2u":{chainId1Bound: {(resId1Bound, resName1Bound):  ((chainId1Unbound, resId1Unbound), i_struct) ...} } ,
                    "u2b":{i_struct: {chainId1Unbound: {(resId1Unbound, resName1Unbound):  (chainId1Bound, resId1Bound) ...} } ,
                          } 
                  }
                       if fixResIds==True else resMap= None
  
  
      resIds are tuples of 3 elements as used in docking project (chainId, resId, resName)
  '''
  if not override and os.path.isfile(outName):
    raise ValueError("Error %s file already exits"%outName)
  print("combining %s to %s"%(fnamesOrStructs_list, outName))
  parser = PDBParser(QUIET=True)
  writer = PDB.PDBIO()
  struct_list=[]
  chains_list=[]
  autoRename=False
  nChains=0
  chainMap={}
  for fnameOrStruct in fnamesOrStructs_list:
    if isinstance(fnameOrStruct, PDB.Structure.Structure):
      struct= fnameOrStruct.copy()
    else:
      struct = parser.get_structure( fnameOrStruct)

    struct_list.append(struct)
    for chain in struct[0]:
      chainId= chain.get_id()
      nChains+=1
      if chainId in chains_list:
        autoRename=True
      else:
        chains_list.append(chainId)
  if nChains> len(ascii_uppercase):
    raise ValueError("Maximun limit of chains %d / %d"%(nChains, len(ascii_uppercase)))
  else:
    letterIter= ( c for c in ascii_uppercase)
  boundStruct= PDB.Structure.Structure("joined")
  boundStruct.add( PDB.Model.Model(0))
  
  resMap_b2u={}
  resMap_u2b={}
  for i_struct, struct in enumerate(struct_list):
    resMap_u2b[i_struct]={}
    for chain in struct[0]:
      unboundChainId= chain.get_id()
      chain._reset_full_id()
      chain.detach_parent()
      if autoRename:
        chain.id= letterIter.next()
      boundChainId= chain.get_id()
      
      if boundChainId==" ":
        boundChainId="*"
      if unboundChainId==" ":
        unboundChainId="*"
        
      chainMap[ boundChainId ]= (unboundChainId, fnamesOrStructs_list[i_struct], i_struct)
      boundStruct[0].add(chain)
      if fixResIds:
        boundResList=[]
        unboundResIdsList=[]
        for resIdNum_new, res in enumerate(chain):
          het, resNum, flag= res.get_id()
          __, resIdUnbound, resName = fromRes2ChainResIdAndName(res)
          unboundResIdsList.append( (unboundChainId, resIdUnbound, resName) )
          res.detach_parent()
          res._reset_full_id()
#          res.id= (het, resIdNum_new, flag) #Remove flag for cases such as 1EAW
          res.id= (het, resIdNum_new, " ")
          boundResList.append(res)
        boundStruct[0].detach_child(chain.get_id())
        chain= PDB.Chain.Chain(boundChainId)
        boundStruct[0].add(chain)
        for res, oldResId_with_format  in zip(boundResList, unboundResIdsList):
          chainIdOld, resIdOld, resNameOld= oldResId_with_format
          if ignoreNoStandard and not is_aa(res, standard=True): continue
          chain.add(res)
          res.set_parent(chain)
          res._reset_full_id()
          newResId_with_format= fromRes2ChainResIdAndName(res)
          chainIdNew, resIdNew, resNameNew=  newResId_with_format
          if chainIdNew not in resMap_b2u:
            resMap_b2u[chainIdNew]= {}
          resMap_b2u[chainIdNew][(resIdNew, resNameNew)]= (oldResId_with_format, i_struct)
          if chainIdOld not in resMap_u2b[i_struct]:
            resMap_u2b[i_struct][chainIdOld]= {}
          resMap_u2b[i_struct][chainIdOld][(resIdOld, resNameOld)]= newResId_with_format

#  print(resMap_u2b)
#  raw_input("enter")
  resMap={"b2u": resMap_b2u,  "u2b": resMap_u2b }
  writer.set_structure(boundStruct)
  writer.save(outName)
  if not fixResIds: resMap={}
  return boundStruct, chainMap, resMap


def test():
  structureL="/home/rsanchez/Tesis/dockScoring/data/develData/decoys/1A2K/1A2K_F23_l_.pdb.gz"
  structureR="/home/rsanchez/Tesis/dockScoring/data/develData/decoys/1A2K/1A2K_F23_r_.pdb.gz"
  boundComplexFname="/home/rsanchez/tmp/myTmp.pdb"
  boundComplex, chainsMap, resMap = combinePDBs( [structureL, structureR], boundComplexFname, override=True,
                                                ignoreNoStandard=False )
  return 

if __name__ == "__main__":
  '''Parses PDB id's desired chains, and creates new PDB structures.'''

  test()
  if len(sys.argv) < 3:
    print("Usage: $ python %s 'pdbOutName' 'pdbIn_0' ['pdbIn_i']" % __file__)
    print("Example:\npython pythonTools/combinePDBs.py pdbPrueba.pdb  "+
          "~/Tesis/rriPredMethod/data/develData/pdbFiles/1A2K_l_b.pdb "+
          "~/Tesis/rriPredMethod/data/develData/pdbFiles/1A2K_r_b.pdb")
    sys.exit()
  else:
    pdbOutName= os.path.expanduser(sys.argv[1])
    pdbNames_list= [ os.path.expanduser(elem) for elem in sys.argv[2:]]
    combinePDBs(pdbNames_list, pdbOutName)


