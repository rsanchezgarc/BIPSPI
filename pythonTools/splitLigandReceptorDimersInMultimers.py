import os
import random

import sys

from .myPDBParser import myPDBParser as PDBParser
from Bio.PDB.Polypeptide  import is_aa
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.NeighborSearch import NeighborSearch 
from joblib import Parallel, delayed

'''
Extracts one chain from a pdb as ligand and the others as ligands
'''
MIN_SEQ_LEN= 30
MAX_SEQ_LEN= 900
N_JOBS=16

parser= PDBParser(QUIET= True)

def findNeigChains(struct, chainIdL, chainIdR, res2res_dist=6, minContacts=20):

  searcher = NeighborSearch( [ atom for atom in struct[0].get_atoms() if atom is not None and not atom.name.startswith("H") ] )
  allNeigs = searcher.search_all(res2res_dist, level="C")
  # print(allNeigs)
  chainL= struct[0][chainIdL]
  chainR= struct[0][chainIdR]
  ligandChains=  set([ chainL ])
  receptorChains= set([ chainR ])
  addedChains= ligandChains.union(receptorChains)
  for neigsGroup in allNeigs:
    searcher = NeighborSearch([atom for chain in neigsGroup  for atom in chain.get_atoms() if atom is not None and
                                                    not atom.name.startswith("H") and atom.get_parent().resname!="HOH"])
    resNeigs = searcher.search_all(res2res_dist, level="R")
    numContacts=0
    for r1,r2 in resNeigs:
      if r1.get_parent().get_id()!= r2.get_parent().get_id():
        numContacts+=1
    # print(neigsGroup, numContacts)
    if numContacts<minContacts:
      continue
    if chainL in neigsGroup or chainR in neigsGroup:
      for chain in neigsGroup:
        if chain not in addedChains:
          if chainR in neigsGroup and chainL in neigsGroup:
            if bool(random.getrandbits(1)):
              receptorChains.add(chain)
            else:
              ligandChains.add(chain)
            addedChains.add(chain)
          elif chainR in neigsGroup:
            receptorChains.add(chain)
            addedChains.add(chain)
          elif chainL in neigsGroup:
            ligandChains.add(chain)
            addedChains.add(chain)

  return ligandChains, receptorChains


def splitOnePDB(fname, chainIdL, chainIdR, outPath):
  print(os.path.basename(fname))
  try:
    s= parser.get_structure(os.path.basename(fname), fname)
  except Exception:
    print ("Error loading pdb")
    return 0

  banLenChains=[]
  try:
    for chain in s[0]:
      badResInChain=0
      for res in  chain.get_list():
        if not is_aa(res,standard=True) and res.resname!="HOH":
          badResInChain+=1
      # for res in chain: print(res)
      chainLen= sum(1 for res in chain if "CA" in res) - badResInChain
      if chainLen < MIN_SEQ_LEN or chainLen > MAX_SEQ_LEN:
        print(chain, chainLen)
        banLenChains.append(chain.get_id())
  except KeyError:
    print ("Not good model")
    return 0

  # print(banLenChains)
  if len( s[0].get_list())-len(banLenChains)<2:
    print(s)
    print( s[0].get_list())
    print("Not enough good chains")
    return 0

  ligandChains, receptorChains= findNeigChains(s, chainIdL, chainIdR)
  print( "ligand:", ligandChains, "receptor:",receptorChains )

  prefix= os.path.basename(fname).split(".")[0]

  io = PDBIO()
  ligandStruct = Structure(prefix + "ligand")
  ligandStruct.add(Model(0))

  for ligandChain in ligandChains:
    ligandChain.set_parent(ligandStruct[0])
    ligandStruct[0].add(ligandChain)
  io.set_structure(ligandStruct)
  io.save(os.path.join(outPath, prefix +"-" + chainIdL + chainIdR + "_l_u.pdb"))

  io=PDBIO()
  receptorStruct= Structure(prefix+"receptor")
  receptorStruct.add(Model(0))
  for receptorChain in receptorChains:
    receptorChain.set_parent(receptorStruct[0])
    receptorStruct[0].add(receptorChain)
  io.set_structure(receptorStruct)
  io.save(os.path.join(outPath, prefix +"-" + chainIdL + chainIdR + "_r_u.pdb"))


def splitFromFile(fname, dirInName, dirOutName):
  args=[]
  with open(fname) as f:
    for line in f:
      lineArray= line.split()
      if len(lineArray)<3:
        continue
      pdbId, chainIdL, chainIdR = lineArray[:3]
      # if pdbId.isupper(): continue
      pdbFname= os.path.join( dirInName, pdbId+".pdb" )
      expectedOutput= os.path.join(dirOutName, pdbId+"-"+chainIdL+chainIdR+"_r_u.pdb")
      if os.path.isfile(pdbFname) and not os.path.isfile(expectedOutput):
        args.append( (pdbFname, chainIdL, chainIdR) )

  Parallel(n_jobs=N_JOBS)( delayed(splitOnePDB)(*(arg+(dirOutName,))) for arg in args)

def test():
  splitOnePDB("/home/ruben/tmp/1kxt.pdb", "B", "A", "/home/ruben/tmp")

if __name__=="__main__":
  # test(); sys.exit();

  __, fname, dirInName, dirOutName= sys.argv
  splitFromFile(fname, dirInName, dirOutName)

