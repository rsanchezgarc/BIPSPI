from __future__ import absolute_import
import sys, os

import numpy as np
from Bio.PDB.Polypeptide import d1_to_index

from ..seqStep.SeqFeatComputer import SeqFeatComputer, FeatureComputerException

from utils import myMakeDir, myMakeDirUnique, tryToRemove, tryToCopy #utils is at the root of the package
from bigExceptions import BadNumberOfResidues

class SeqInputPreproceser(SeqFeatComputer):
  '''
    Extends SeqFeatComputer class. Creates a putative contact map given 2 fasta sequences.
  '''
  def __init__(self, rFname, lFname, computedFeatsRootDir= None, statusManager=None):
    '''
      @param rFname: str. path to receptor fasta file
      @param lFname: str. path to ligand fasta file     
      @param computedFeatsRootDir: str. path where features will be stored
      @param statusManager: class that implements .setStatus(msg) to communicate
    '''
    SeqFeatComputer.__init__(self, rFname, lFname, computedFeatsRootDir, statusManager=statusManager)
#    at this point: computedFeatsRootDir= /path_to_computedFeatsRootDir/
#                   self.computedFeatsRootDir= /path_to_computedFeatsRootDir/seqStep This will be overrode

    self.outPath= myMakeDir(computedFeatsRootDir, "common/contactMaps")
    self.prefixR= os.path.split(rFname)[1].split(".")[0].split("_")[0]
    self.prefixL= os.path.split(lFname)[1].split(".")[0].split("_")[0] 
    if self.prefixR == self.prefixL:
      self.prefix= self.prefixR
    else:
      if "<" in self.prefixL:
        raise FeatureComputerException("Error. Partner 1 fasta file %s must not contain '<' or '>' character"%lFname)
      if ">" in self.prefixR:
        raise FeatureComputerException("Error. Partner 1 fasta file %s must not contain '<' or'>' character"%rFname)
      self.prefixR= self.getExtendedPrefix(rFname)
      self.prefixL= self.getExtendedPrefix(lFname)      
        
      self.prefix= self.prefixL+"<->"+self.prefixR

    self.outName= os.path.join(self.outPath,self.prefix+".cMap.tab")    
    self.computeFun= self.contactMapOneComplex
    
  def contactMapOneComplex(self):
    '''
      Computes the contact map of a complex. Initial input for complex codification. Contact map is a file written at
      self.computedFeatsRootDir/common/contactMaps/ with name prefix.cMap.tab where prefix is either the common name of
      ligand and receptor pdb files or the concatenation of ligand and receptor names.
      1A2K_l_u.pdb and 1A2K_r_u.pdb  --> 1A2K.cMap.tab
      1A2K_l_u.pdb and 1A22.pdb  --> 1A2K-1A22.cMap.tab
      
    '''    
    outName= self.outName
    print (outName)
    if os.path.isfile(outName):
      print ('Already computed contact map')
      return 0

    seqL =  self.parseFasta(self.lFname, inputNumber="1")
    seqR =  self.parseFasta( self.rFname, inputNumber="2")
#    print(repr(seqL))
#    print(repr(seqR))
    nResiduesL= len(seqL)
    nResiduesR= len(seqR)
    if not (self.minNumResiduesPartner< nResiduesL < self.maxNumResiduesPartner):
      raise BadNumberOfResidues(nResiduesL, "1")
    if not (self.minNumResiduesPartner< nResiduesR < self.maxNumResiduesPartner):
      raise BadNumberOfResidues(nResiduesL, "2")
      
    with open(outName,"w") as outFile:
      outFile.write("chainIdL structResIdL resNameL chainIdR structResIdR resNameR categ\n")
      try:
        for ixL, resnameL in enumerate(seqL):
          chainIdL="L"
          resIdL= str(ixL)
          if not resnameL in d1_to_index: continue
          for ixR, resnameR in enumerate(seqR):
            if not resnameR in d1_to_index: continue          
            chainIdR="R"
            resIdR= str(ixR)
            categ= np.nan
#            print("%s %s %s %s %s %s %s\n" %(chainIdL, resIdL, resnameL, chainIdR, resIdR, resnameR, categ))
            outFile.write("%s %s %s %s %s %s %s\n" %(chainIdL, resIdL, resnameL, chainIdR, resIdR, resnameR, categ))
      except (KeyboardInterrupt, Exception):
        print("Exception happend computing %s"%outName)
        tryToRemove(outName)    
        raise

def testModulePredict():
  
  lFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/inputFasta/seq1_l_u.fasta"
  rFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/inputFasta/seq1_r_u.fasta"
  
  computedFeatsRootDir="/home/rsanchez/tmp/computedFeats"

  comput= SeqInputPreproceser(rFname, lFname, computedFeatsRootDir)
  comput.computeFun()
  raw_input("test_done")
  raise FeatureComputerException("Debug Mode. Comment out testModule() line 105")
      
            
if __name__=="__main__":

  testModulePredict()

