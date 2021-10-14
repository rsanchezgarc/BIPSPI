from __future__ import absolute_import
import os
from computeFeatures.toolManagerGeneric import loadPdbIfIsPath
from Bio.PDB.PDBIO import PDBIO
from ..featuresComputer import FeaturesComputer, FeatureComputerException
from utils import myMakeDir, tryToRemove
from bigExceptions import BadNumberOfResidues

from .structTools.psaiaManager import PSAIAComputer
from .structTools.halfSphereManager import HalfSphereComputer
from .structTools.dsspManager import DsspComputer

class StructFeatComputer(FeaturesComputer):
  '''
  Intended to be used for computing all structural features for one complex
  '''
  def __init__(self, prefix, computedFeatsRootDir= None, statusManager=None):
    '''
      @prefix. An id for a given complex. E.g. 1A2K
      :param computedFeatsRootDir: str. path where features will be stored
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''
    FeaturesComputer.__init__(self, prefix=prefix, computedFeatsRootDir=computedFeatsRootDir, statusManager=statusManager)
    self.computedFeatsRootDir= myMakeDir(self.computedFeatsRootDir, "structStep")
    
    self.psaiaComputer= PSAIAComputer(self.computedFeatsRootDir, statusManager= self.statusManager)
    self.hsComputer= HalfSphereComputer(self.computedFeatsRootDir, statusManager= self.statusManager)
    self.dsspComputer= DsspComputer(self.computedFeatsRootDir, statusManager= self.statusManager)
    self.allComputers= [self.psaiaComputer, self.hsComputer, self.dsspComputer]

    
  def computeComplex(self, fnameL, fnameR, structureL, structureR):
    '''
    Applies self.computeOneFile method to both ligand and receptor pdbFiles.
    :param fnameOrStructL: str or Bio.PDB.Structure. fname or structure to pdbfile or fname to fasta of ligand
    :param fnameOrStructR: str or Bio.PDB.Structure. fname or structure to pdbfile or fname to fasta of receptor
    :param structureL: Bio.PDB.Structure.Structure. Structure of ligand protein (contained in fnameL). None if fasta
    :param structureR: Bio.PDB.Structure.Structure. Structure of receptor protein (contained in fnameR). None if fasta
    '''
    try:
      for computer in self.allComputers:
        computer.computeComplex( fnameL, fnameR, structureL, structureR  )
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing structural features for %s"%(self.prefix))
      raise
    return 0
      
  def getChainsIdsFromStruct(self, structure):
    chainsId= [ chain.get_id()  for chain in structure[0] ]
    chainsId= set([ chainId if chainId!=" " else type(self).UNKNOWN_CHAIN_SYMB for chainId in chainsId ])
    return chainsId
      
def moveAndWriteAsPDBIfMmcif(fnameIn, fnameOut, removeInput=False):
  from Config import Configuration
  conf= Configuration()
  minNumResidues, maxNumResidues= conf.minNumResiduesPartner, conf.maxNumResiduesPartner
  try:
    struct, __= loadPdbIfIsPath(fnameIn)
    totalNumRes=0
    for chain in struct[0]:
      nResInChain= len(chain.get_list())
      totalNumRes+= nResInChain
    if not ( minNumResidues < totalNumRes < maxNumResidues):
      raise BadNumberOfResidues(totalNumRes)
    else:
      writter=PDBIO()
      writter.set_structure(struct)
      writter.save(fnameOut)
      if removeInput: os.remove(fnameIn)
      return True
  except Exception as e:
    print("Error in moveAndWriteAsPDBIfMmcif !!!", e)
    return False
    
