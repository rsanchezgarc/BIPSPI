from __future__ import absolute_import

from ...toolManagerGeneric import ToolManager
from ...featuresComputer import   FeatureComputerException

class StructToolManager(ToolManager):
  '''
  Abastract class that will be used to implement psaia, half sphere ... calculators
  
  '''
           
  def __init__(self, computedFeatsRootDir, statusManager=None):
    '''
      :param computedFeatsRootDir: str. root path where results will be saved
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''   
    ToolManager.__init__(self, computedFeatsRootDir, statusManager)
      
      
  def computeComplex(self, fnameL, fnameR, structureL, structureR):
    '''
    Applies self.computeOneFile method to both ligand and receptor pdbFiles.
    :param fnameL: str. fname to pdbfile or fname to fasta of ligand
    :param fnameR: str. fname to pdbfile or fname to fasta of receptor
    :param structureL: Bio.PDB.Structure.Structure. Structure of ligand protein (contained in fnameL). None if fasta
    :param structureR: Bio.PDB.Structure.Structure. Structure of receptor protein (contained in fnameR). None if fasta
    '''
    if structureL:
      self.computeOneFile(fnameL, structureL)
    if structureR:
      self.computeOneFile(fnameR, structureR)
    
  def computeOneFile(self, fnameOrStruct):
    '''
      Computes PSAIA, halfsphere,... for a given pdb file
      :param fnameOrStruct: str or Bio.PDB.Structure. fname or structure to pdbfile or fname to fasta of receptor
      :param struct: ignored
    '''
    raise ValueError("Not implemented")
    

  
      
if __name__=="__main__":
  print("Done")
