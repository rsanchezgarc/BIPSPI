from __future__ import absolute_import
import os,sys
from ..FeaturesComputer import FeaturesComputer, FeatureComputerException
from .seqToolManagers.conservationTools.windowPSSM import AA_CODE_ELEMENTS
from utils import myMakeDir, myMakeDirUnique, tryToRemove #utils is at the root of the package
from bigExceptions import BadSequence
import re

class SeqFeatComputer(FeaturesComputer):
  '''
  Abstract class. It will be extended with different features computer classes, p.e, PSSM's.
  Intended to be used for computing one type of features
  '''
  cleanSeqRegEx= re.compile(r"\n|\s|\r")
  def __init__(self, rFname, lFname, computedFeatsRootDir= None, statusManager=None):
    '''
      @param rFname: str. path to receptor pdb file
      @param lFname: str. path to ligand pdb file      
      @param computedFeatsRootDir: str. path where features will be stored. If None they will be stored
                                        at default path (assigned in ../Config.py)
      @param statusManager: class that implements .setStatus(msg) to communicate
    '''
    FeaturesComputer.__init__(self, rFname, lFname, computedFeatsRootDir, statusManager=statusManager)
    self.computeFun= self.computeOneComplex
    self.computedFeatsRootDir= myMakeDir(self.computedFeatsRootDir, "seqStep")

  def computeOneComplex(self):
    '''
      Applies self.computeOneFile method to both ligand and receptor pdbFiles.
    '''
    self.computeOneFile( self.rFname, chainType="r")
    self.computeOneFile( self.lFname, chainType="l")
    return 0

  def computeOneFile(self, fileName):
    '''
      abstract method
    '''
    return None
    
  def checkIfIsFasta(self, fileName):
    with open(fileName) as f:
      head= f.readline()
      if head.startswith(">"):
        for line in f:
          for char in line:
            if char.isdigit():
              return False
        return True
      else:
        return False
        
        
  def parseFasta(self, fname, inputNumber="1"):
    with open(fname) as f:
      head= f.readline()
      if not head.startswith(">"):
        raise BadSequence("Error, wrong fasta format for partner %s or several sequences"%inputNumber)
      seq=""
      for line in f:
        if ">" in line:
          raise BadSequence("Error for partner %s, Just one sequnce allowed in fasta files"%inputNumber)
        else:
          seq+= re.sub(SeqFeatComputer.cleanSeqRegEx,"",line.strip().upper())
    seq= parseSeq(seq)
    if re.match(r".*(;|\||&)", seq):
      raise BadSequence("Error for partner %s, bad character in sequence"%inputNumber)
    return seq


def parseSeq(seq):

  if seq[0]==">":
    seq= "".join(seq.split("\n")[1:])
    
  seq= re.sub(SeqFeatComputer.cleanSeqRegEx,"",seq.strip().upper())
  for letter in seq:
    if letter not in AA_CODE_ELEMENTS:
      raise BadSequence("Error unexpected character '%s' found in %s"%(letter, seq))
  return seq
