from __future__ import absolute_import
import os,sys
from Bio.PDB.Polypeptide import three_to_one
from computeFeatures.structStep.myPDBParser import myPDBParser as PDBParser
from Bio.PDB.PDBIO import PDBIO
from ..FeaturesComputer import FeaturesComputer, FeatureComputerException
from utils import myMakeDir, tryToRemove #utils is at the root of the package
from Config import Configuration
from bigExceptions import BadNumberOfResidues

class StructFeatComputer(FeaturesComputer):
  '''
  Abstract class. It will be extended with different features computer classes, p.e, PSAIA computer
  Intended to be used for computing one type of features
  '''
  def __init__(self, rFname, lFname, computedFeatsRootDir= None, statusManager=None):
    '''
      @param rFname: str. path to receptor pdb file
      @param lFname: str. path to ligand pdb file      
      @param computedFeatsRootDir: str. path where features will be stored
      @param statusManager: class that implements .setStatus(msg) to communicate
    '''
    FeaturesComputer.__init__(self, rFname, lFname, computedFeatsRootDir, statusManager=statusManager)
    self.computeFun= self.computeOneComplex
    self.computedFeatsRootDir= myMakeDir(self.computedFeatsRootDir, "structStep")

  def threeLetterAA_to_one(self, aa3Letters):
    '''
      Return the one letter code for a given 3 letters code
      @param aa3Letters: str
      return oneLetterCode: str
      
    '''
    try:
      return three_to_one(aa3Letters)
    except KeyError:
      return "X"
  def computeOneComplex(self):
    '''
      Applies self.computeOneFile method to both ligand and receptor pdbFiles.
    '''
    self.computeOneFile( self.rFname)
    self.computeOneFile( self.lFname)
    return 0
    
  def computeOneFile(self, fileName):
    '''
      abstract method
    '''
    return None
    

def moveAndWriteAsPDBIfMmcif(fnameIn, fnameOut, removeInput=False):
  conf= Configuration()
  minNumResidues, maxNumResidues= conf.minNumResiduesPartner, conf.maxNumResiduesPartner
  try:
    parser= PDBParser(QUIET=True)
    struct= parser.get_structure("pdbStruct", fnameIn)
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
    
