from __future__ import absolute_import, print_function

import gzip
import os, shutil
import pandas as pd

from .featuresComputer import FeaturesComputer
from .toolManagerGeneric import checkIfIsFasta, loadPdbIfIsPath

from .common.computeContactMap import ContactMapper
from .seqStep.seqFeatsComputer import SeqFeatComputer
from .structStep.structFeatsComputer import StructFeatComputer

from pythonTools.myPDBParser import myPDBParser
from utils import myMakeDir
from bigExceptions import NoValidPDBFile

class OneComplexFeatComputer(FeaturesComputer):

  def __init__(self, prefix, computedFeatsRootDir= None, methodProtocol="struct", areForTrainAndTest=False, 
                    boundAvailable=False, statusManager= None):
    '''
    :param prefix: str. Id for a given complex
    :param computedFeatsRootDir: str. Path where features files will be saved. By default it uses
                                Config.py DEFAULT_PARAMETERS["computedFeatsRootDir"] will be used as out_path
    :param methodProtocol: str. "seq" if just sequential features will be used; "struct" if sequential and
                                structural features will be used. "mixed" behaves equally than struct
    :param areForTrainAndTest: boolean. True if ligand and receptor are in interacting coordinates and thus,
                                contact maps are needed to be computed in order to perform evaluation.
                                False if ligand and receptor coordinates are not related and thus, 
                                evaluation does not makes sense.
    :param: boundAvailable. True if there is a bound and unbound pdb for each complex. False otherwise
    :param statusManager: class that implements .setStatus(msg) to communicate
    '''
    FeaturesComputer.__init__(self, prefix, computedFeatsRootDir, statusManager= statusManager)
    if methodProtocol=="mixed": methodProtocol="struct"
    self. methodProtocol= methodProtocol
    self.areForTrainAndTest= areForTrainAndTest
    self.boundAvailable= boundAvailable
    assert self.methodProtocol in ["seq", "struct"], 'Error methodProtocol in computeFeaturesAllPdbsOneDir ' +\
                                                     ' must be "seq", "struct", or "mixed" found -> %s'%(methodProtocol)
    self.computedFeatsRootDir= os.path.expanduser(self.computedFeatsRootDir) #Creates root path where features will be saved
    self.computedFeatsRootDir=  myMakeDir(self.computedFeatsRootDir)
    

  def computeFeaturesOneComplex(self, fnameL, fnameR, lPdbId=None, rPdbId=None, isHomoComplex=False):
    '''
      Computes all features needed for complex codification.
      :param fnameL: str. Path to the the pdb or fasta file of the ligand of the complex
      :param fnameR: str. Path to the the pdb or fasta file of the receptor of the complex
      :param lPdbId: str. pdbId of receptor in order to query 3dCons. If None, psiblast will be launched directly
      :param rPdbId: str. pdbId of receptor in order to query 3dCons. If None, psiblast will be launched directly
      :param isHomoComplex: boolean. If True, just lPartner is provided and both partners are equal.
    '''

    prefixL, chainTypeL = self.splitExtendedPrefix(self.getExtendedPrefix(fnameL))[:2]
    prefixR, chainTypeR = self.splitExtendedPrefix(self.getExtendedPrefix(fnameR))[:2]
    assert prefixL==prefixR , "Error, prefixes are different for %s - %s"%(fnameL, fnameR)
    
    structureL, __ =  loadPdbIfIsPath( fnameL ) #structureL will be none if fnameL is a fasta file
    structureR, __ =  loadPdbIfIsPath( fnameR ) #structureR will be none if fnameL is a fasta file
    if (structureL is None or structureR is None) and self.areForTrainAndTest:
      raise NoValidPDBFile("Error, Training requeries both inputs to be pdbs. \n"+
                            "For prediction they can be either pdbs or sequences. You are training")
    
    cMapper= ContactMapper(self.prefix, computedFeatsRootDir= self.computedFeatsRootDir, 
                           areForTrainAndTest= self.areForTrainAndTest, boundAvailable= self.boundAvailable, 
                           statusManager= self.statusManager)

    if isHomoComplex and not self.areForTrainAndTest:
      resDictL, resDictR= cMapper.computeComplex( fnameL, None, structureL, None)
    else:
      resDictL, resDictR= cMapper.computeComplex( fnameL, fnameR, structureL, structureR, isHomoLR= isHomoComplex)

    seqFeatComputer= SeqFeatComputer(self.prefix, computedFeatsRootDir= self.computedFeatsRootDir, statusManager= self.statusManager)

    if not self.methodProtocol.startswith("seq"):
      structComputer= StructFeatComputer(self.prefix, computedFeatsRootDir= self.computedFeatsRootDir, statusManager=self.statusManager)
      if isHomoComplex and not self.areForTrainAndTest:
        structComputer.computeComplex(fnameL, None, structureL, None)
      else:
        structComputer.computeComplex(fnameL, fnameR, structureL, structureR)


    seqFeatComputer.computeComplex(fnameL, fnameR, structureL, structureR, lPdbId=lPdbId, rPdbId=rPdbId,
                                   areLRequivalentProteins= isHomoComplex and not self.areForTrainAndTest)

    if isHomoComplex and not self.areForTrainAndTest:
      self.copyFeaturesToReceptor(fnameL, fnameR, isSeqInput=self.methodProtocol.startswith("seq"))

  def copyFeaturesToReceptor(self, fnameL, fnameR, isSeqInput=False):
    assert fnameL== fnameR, "error, in homoComplex mode both fnames must be equal"
    dirName, baseName= os.path.split(fnameL)
    shutil.copyfile(os.path.join(dirName, baseName), os.path.join(dirName, baseName.replace("_l_", "_r_") )) 
    for root, dirs, files in os.walk(self.computedFeatsRootDir):
      for fName in files:
        if "_l_" in fName:
          if not isSeqInput:
            shutil.copyfile(os.path.join(root, fName), os.path.join(root, fName.replace("_l_", "_r_")))
          else:
            newName= fName.replace("_l_L", "_r_R")
            try:
              data= pd.read_csv(os.path.join(root, fName),  sep='\s+', comment="#", dtype={"resIdL":str, "resIdR":str,
                                                                                           "resId":str, "chainId":str,
                                                                                           "chainIdL":str, "chainIdR":str})
              comments= ""
              with gzip.open(os.path.join(root, fName)) as f:
                for line in f:
                  if line.startswith("#"):
                    comments+=line
                  else:
                    break
              if "chainId" in data:
                data["chainId"]="R"
              else:
                data["chainIdR"] = "R"
              if len(comments)==0:
                data.to_csv(os.path.join(root, newName), index=False, sep=" ", compression="gzip")
              else:
                with gzip.open(os.path.join(root, newName),"w") as f:
                  f.write(comments)
                  data.to_csv(f, index=False, sep=" ") #, compression="gzip")

            except (pd.errors.ParserError, IOError):
              continue
      
def computeFeaturesOneComplex(  lFname, rFname, prefix=None, lPdbId=None, rPdbId= None, computedFeatsRootDir= None,
                                   methodProtocol=None, areForTrainAndTest=None, isHomoComplex=False, statusManager= None):
  if prefix is None:
    if lPdbId is None or rPdbId is None:
      p_l = os.path.basename(lFname).split("_")[0].split(":")[0].split("@")[0]
      p_r = os.path.basename(rFname).split("_")[0].split(":")[0].split("@")[0]
      if p_l == p_r:
        prefix= p_l
      else:
        raise ValueError(" fnames mismatch %s %s"%(lFname,rFname))
    else:
      prefix= lPdbId+":"+rPdbId
  occ= OneComplexFeatComputer(prefix, computedFeatsRootDir= computedFeatsRootDir, methodProtocol=methodProtocol,
                              areForTrainAndTest=areForTrainAndTest,statusManager= statusManager)                            
  occ.computeFeaturesOneComplex(  lFname, rFname, lPdbId=lPdbId, rPdbId=rPdbId, isHomoComplex=isHomoComplex)      
  
  
