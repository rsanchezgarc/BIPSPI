from __future__ import absolute_import

import os
import re
import subprocess
import gzip
from Bio.PDB.Polypeptide import three_to_one, one_to_three
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from itertools import chain as IterChainer

from .featuresComputer import getExtendedPrefix, splitExtendedPrefix
from Config import Configuration
from bigExceptions import BadSequence
from pythonTools.myPDBParser import myPDBParser
from utils import openForReadingFnameOrGz, tryToRemove, tryToMove, tryToSymlink

FILTER_OUT_LABELS=False
IGNORE_NO_STANDARD= False
UNKNOWN_CHAIN_SYMB="+"
AA_CODE_ELEMENTS= ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","Z","X"] #Z stands for no residue, X non standard

class ToolManager(Configuration):
  '''
  Abastract class that will be used to implement psiblast, psaia, etc computers
  '''
  HEADER_SINGLE_CHAIN= "chainId resId resName %(features)s\n"
  HEADER_PAIRWISE= "chainIdL resIdL resNameL chainIdR resIdR resNameR %(features)s\n"
  UNKNOWN_CHAIN_SYMB= UNKNOWN_CHAIN_SYMB
  AA_CODE_ELEMENTS= AA_CODE_ELEMENTS
  def __init__(self, computedFeatsRootDir, statusManager=None):
    '''
      :param computedFeatsRootDir: str. root path where results will be saved
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''    
    Configuration.__init__(self)
    self.computedFeatsRootDir= computedFeatsRootDir
    self.filterOutLabels= FILTER_OUT_LABELS
    self.filterOutNoStandard= IGNORE_NO_STANDARD
    self.statusManager= statusManager

  def computeComplex(self, fnameL, fnameR, structureL, structureR):
    '''
    :param fnameL: str. fname to pdbfile or fname to fasta of ligand
    :param fnameR: str. fname to pdbfile or fname to fasta of receptor
    :param structureL: Bio.PDB.Structure.Structure. Structure of ligand protein (contained in fnameL). None if fasta
    :param structureR: Bio.PDB.Structure.Structure. Structure of receptor protein (contained in fnameR). None if fasta
    '''
    raise ValueError("Not implemented")
     
  def getFNames(self, prefixExtended):
    '''
    abstract method
    Returns a list that contains the fnames that will contain the computed features.
    :param prefixExtended. prefix for output fnames.
    :return list of fnames: [ fname1, fnam2, ...]
    '''
    raise ValueError("Not implemented")
    
  def makePrefixExtendedPairwise(self, prefix):
    return prefix+"_"

  def makePrefixExtendedSingle(self, prefix, chainType):
    return prefix+"_"+chainType+"_"
    
  def loadPdbIfIsPath(self, fnameOrStruct):
    return loadPdbIfIsPath( fnameOrStruct)
      
  def getNLines(self, fname):
    if fname.endswith(".gz"):
      return int(subprocess.check_output('zcat {} | wc -l'.format(fname), shell=True).split()[0])
    else:
      return int(subprocess.check_output('wc -l {}'.format(fname), shell=True).split()[0])
      
  def checkAlreayComputed(self, prefixExtended):
    '''
    checks if output files have been computed
    :param prefixExtended. prefix for output fnames.
    :return boolean. True if prefixExtended sequences have been already computed, False otherwise
    '''  

    for fname in self.getFNames( prefixExtended):
      if os.path.isfile(fname):
        nlines= self.getNLines(fname)
      elif os.path.isfile(fname+".gz"):
        nlines= self.getNLines(fname+".gz")
      else:
        return False
      if nlines< self.minNumResiduesPartner and "joblib" not in fname: #Too small file to be correct
        print("fname %s does not exists"%(fname))
        return False
    return True


  def tryToRemoveAllFnames(self, prefixExtended):
    '''
    try to remove all fnames returned by getFNames (useful to clean if some exception happens
    :param prefixExtended. prefix for output fnames.
    '''  
    for fname in self.getFNames( prefixExtended):
      if os.path.isfile(fname):
        tryToRemove(fname)
    

  def fromRes2ChainResIdAndName(self, resObj):
    return fromRes2ChainResIdAndName(resObj)
    
  def fromResObj2ResIdStr(self, resObj):
    return fromResObj2ResIdStr(resObj)
    
  def makeStrResId(self, resId):
    return makeStrResId(resId)
    
  def threeLetterAA_to_one( self, candidateResName):
    return threeLetterAA_to_one(candidateResName)
    
  def getExtendedPrefix(self, fname):
    '''
      Given a filename, obtains its unambiguous id
      :param fname: str. A filename. pe. "/path/to/file/1A2K_l_.pdb"
      :return unambiguous id: str. pe "1A2K_l_"  This id will be the prefix of output names
    '''
    return getExtendedPrefix(fname)
    
  def splitExtendedPrefix(self, extendedPrefix):
    '''
      Given an extendedPrefix, splits it into: prefix, chainType ,[chainId]
      :param fname: str. An extendedPrefix. pe. "1A2K_l_"
      :return    prefix,  chainType,
             or  prefix,  chainType, chainId
    '''
    return splitExtendedPrefix(extendedPrefix)
    
  def openForReadingFnameOrGz(self, fname):
    return openForReadingFnameOrGz(fname)

  def checkIfIsFasta(self, fileName):    
    return checkIfIsFasta(fileName)

  def parseFasta( self, fname, inputNumber="1"):
    return parseFasta(fname, inputNumber)
    
  def filterOutSmallDataDict(self, dataDict):
    '''
    Remove pdb chains that are too small to be considered
    :param dataDict:{chainId:[listOfDataRows]}
    '''
    for chainId in list(dataDict.keys()):
      if len(dataDict[chainId])<= self.minNumResiduesPartner: #Too small chains will not be considered
        del dataDict[chainId]
    return dataDict
    
  def getStructFromFasta(self, fname, chainType):
    '''
    Creates a Bio.PDB.Structure object from a fasta file contained in fname. Atoms are not filled
    and thus no coordiantes availables. Implements from Structure to Residue hierarchy.
    :param fname: str. path to fasta file
    @chainType: str. "l" or "r"
    '''
    
    seq= self.parseFasta( fname, inputNumber="1" if chainType=="l" else "2") #inpuNumber is used to report which partner fails if error
    prefix= self.splitExtendedPrefix(self.getExtendedPrefix(fname))[0]
    chainId= chainType.upper()
    residues= []
    struct= Structure(prefix)
    model= Model(0)
    struct.add(model)
    chain= Chain(chainId)
    model.add(chain)
    for i, aa in enumerate(seq):
      try:
        resname= one_to_three(aa)
      except KeyError:
        resname= "UNK"
      res= Residue( (' ', i, ' '), resname, prefix)
      chain.add(res)
    return struct
    
  def writeResultsFromDataDictSingleChain(self, dataDict, outName, featuresNames=None, categoricalLevels=None):
    '''
    Used to uniformly write features that describe a single residue in a chain
    :param dataDict: {chainId: [list_of_str_that_represents_a_row ] }
                 e.g {"A":["A 123 L 0.1 0.8 -1", ...]  #3 first fields are chainId resId and resName
    :param outName: str: path where results will be saved
    :param featuresNames: [ str]. A list of names for features. If None, type(self).VAR_LIST will be used
    :param categoricalLevels: { ("categ1A, categ2B..,): [featName1, featName2...]}
        featNames must be included in featuresNames
                e.g.    { ("A","V","L","H"...): ["aaSymbol-1", "aaSymbol0", "aaSymbol+1"],
                           ("H","L","T","-"): ["secondStruc"]}
    '''
    dataDict= self.filterOutSmallDataDict(dataDict)
    listOfRecords= list(IterChainer.from_iterable( [sorted(dataDict[chainId]) for chainId in dataDict ] ))
    featuresNames=  featuresNames if featuresNames else type(self).VAR_LIST
    headerStr=""
    if categoricalLevels:
      headerStr=  "#Levels: "
      for tupleOfFeatNames in categoricalLevels:
        headerStr+= ";".join(tupleOfFeatNames)+":"+";".join( categoricalLevels[tupleOfFeatNames] )+" "
      headerStr= headerStr[:-1]+"\n"
    headerStr+= ToolManager.HEADER_SINGLE_CHAIN%{"features": " ".join( featuresNames) }
    self.writeGzResults(outName, headerStr, listOfRecords)
      
  def writeResultsFromDataDictPairL2R(self, dataList, outName, featuresNames=None):
    '''
    Used to uniformly write features that describe a pair of residues, each from a different partner
    :param dataList: a list of str that represents the rows of the dataframe
                 e.g [ "A 123 L B 2 I 0.1 0.8 -1", ... ] #6 first fields are chainIdL resIdL resNameL 
                                                                          chainIdR resIdR and resNameR
    :param outName: str: path where results will be saved
    :param featuresNames: [ str]. A list of names for features. If None, type(self).VAR_LIST will be used
    '''
    featuresNames=  featuresNames if featuresNames else type(self).VAR_LIST
    headerStr= ToolManager.HEADER_PAIRWISE%{"features": " ".join( featuresNames)  }
    self.writeGzResults(outName, headerStr, dataList)
    
  def writeGzResults(self, outName, headerStr, listOfRecords):
    '''
    Used to write a list of records as a .gz file

    :param outName: str. path where results will be saved
    :param headerStr: str. string that will be written at the beginning of the file
    :param dataList: a list of str that represents the rows of the dataframe
                 e.g [ "A 123 L B 2 I 0.1 0.8 -1", ...]
    '''
    dirName, baseName = os.path.split(outName)
    tmpOutName= os.path.join(dirName, "tmp-"+baseName)
    try:
      with gzip.open(tmpOutName,"w") as outFile:
        outFile.write(headerStr)
        outFile.write("\n".join(listOfRecords))
      tryToMove(tmpOutName, outName)
    except (KeyboardInterrupt, Exception):
      tryToRemove(outName)
      tryToRemove(tmpOutName)      
      raise    
      
  def uncompressFile(self, fileName, destPath=None):
    return uncompressFile(fileName, destPath)
    
def fromRes2ChainResIdAndName( resObj ):
  '''
  Given a Residue object,  obtain its  chainId resId and resName
  as strings. resId is the concatenation of resId number and flag.
  :param resObj: Bio.PDB.Residue
  returns chainId:str, resId:str, resName:str
  '''
  pdbId, model, chainId, resId= resObj.get_full_id()
  if chainId==" ": chainId=UNKNOWN_CHAIN_SYMB
  return chainId, makeStrResId(resId), threeLetterAA_to_one(resObj.resname)
  
def fromResObj2ResIdStr( resObj):
  return makeStrResId(resObj.get_id())
  
def makeStrResId( resId):
  '''
  :param resId: (str1, int, str2)
  resId is the concatenation of resId number and flag.-> str(int)+str2
  '''
  valList= [ str(elem) for elem in resId[1:]]
  finalId= "".join(valList).strip()
  return finalId

def threeLetterAA_to_one( candidateResName):
  '''
    Return the one letter code for a given 3 letters code
    :param candidateResName: str
    return oneLetterCode: str
    
  '''
  try:
    return three_to_one( candidateResName)
  except KeyError:
    if candidateResName in ["H1S", "H2S", "H3S"]:
      return "H"
    elif candidateResName in ["MSE"]:
      return "M"
    elif candidateResName in ["PTR", "TPO"]:
      return "T"
    elif candidateResName in ["SEP"]:
      return "S"
    elif candidateResName in ["HYP"]:
      return "P"
    elif candidateResName in ["TYS"]:
      return "T"
    elif candidateResName in ["MLZ", "MLY", "M3L"]:
      return "K"
    elif candidateResName in ["HOH"]:
      return None
    else:
      return "X"
              
def loadPdbIfIsPath( fnameOrStruct):
  '''
  given a fname or a Bio.PDB.Structure
  
  return None, fname if fname is a fasta file
        Bio.PDB.Structure, fname otherwise
  '''
  if not isinstance(fnameOrStruct, Structure):
    if checkIfIsFasta(fnameOrStruct):
      return None, fnameOrStruct
    else:
      try:
        parser= myPDBParser(QUIET=True)
        struct= parser.get_structure(fnameOrStruct)
        return parser.get_structure(fnameOrStruct),  fnameOrStruct
      except Exception:
        raise
  else:
    return fnameOrStruct, fnameOrStruct.get_id()
    

def checkIfIsPdb(fname):
  try:
    for line in openForReadingFnameOrGz(fname):
      if line.startswith("ATOM"):
        if line[42]=="." and line[50]==".":
          return True
    return False
  except Exception:
    return False
    
CLEAN_SEQ_REGEX= re.compile(r"\n|\s|\r")

def checkIfIsFasta(fname):
  try:
    parseFasta(fname)
  except BadSequence:
    return False
  return True
  
def parseFasta( fname, inputNumber="1"):
  '''
  Returns the sequence contained in fasta or fasta.gz as a str
  @fname: the path to a fasta file
  :param inputNumber: str. Id of the file. Used to report which partners has failed if
                           parsing was unsucessful
  '''
  if fname.endswith(".gz"):
    openFun= gzip.open
  else:
    openFun= open
  with openFun(fname) as f:
    head= f.readline()
    if not head.startswith(">"):
      raise BadSequence("Error, wrong fasta format for partner %s or several sequences"%inputNumber)
    seq=""
    for line in f:
      if ">" in line:
        raise BadSequence("Error for partner %s, Just one sequnce allowed in fasta files"%inputNumber)
      else:
        seq+= re.sub(CLEAN_SEQ_REGEX,"",line.strip().upper())
  seq= parseSeq(seq)
  if re.match(r".*(;|\||&)", seq):
    raise BadSequence("Error for partner %s, bad character in sequence"%inputNumber)
  return seq

def parseSeq(seq):

  if seq[0]==">":
    seq= "".join(seq.split("\n")[1:])
    
  seq= re.sub( CLEAN_SEQ_REGEX,"",seq.strip().upper())
  for letter in seq:
    if letter not in AA_CODE_ELEMENTS:
      raise BadSequence("Error unexpected character '%s' found in %s"%(letter, seq))
  return seq  
  
def uncompressFile(fileName, destPath=None):
  if destPath:
    __, fname_base= os.path.split(fileName)
  else:
    destPath, fname_base= os.path.split(fileName)
  if fname_base.endswith(".gz"):
    uncompressFileName= os.path.join(destPath, fname_base.split(".gz")[0])
    command= " gzip -dcf %s  > %s"%( fileName, uncompressFileName)
    subprocess.check_call(command, shell=True)
  else:
    uncompressFileName= os.path.join(destPath, fname_base)
    tryToSymlink(fileName, uncompressFileName)
  return uncompressFileName
  
def compressFile( fileName):
  if fileName.endswith(".gz"):
    return fileName
  command= " gzip %s"%( fileName)
  subprocess.check_call(command, shell=True)
  return fileName+".gz"
  

if __name__=="__main__":
  print("Done")
