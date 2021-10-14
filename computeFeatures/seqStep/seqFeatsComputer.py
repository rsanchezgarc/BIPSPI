from __future__ import absolute_import
import os
import gzip
from ..featuresComputer import FeaturesComputer
from utils import myMakeDir, tryToSymlink #utils is at the root of the package
from bigExceptions import BadSequence
import re

from .seqToolManagers.seqExtraction.seqStructMapper import SeqStructMapper
from .seqToolManagers.aaTypeTools.windowSeqManager import WindowSeqManager
from .seqToolManagers.conservationTools.PsiBlastManager import PsiBlastManager
from .seqToolManagers.conservationTools.Al2coManager import Al2coManager

from .seqToolManagers.conservationTools.HHblitsManager import HHBlitsManager
#from .seqToolManagers.conservationTools.PSICOVManager import PsicovManager as CorrMutManager, HHBLITS_CMD_TEMPLATE
from .seqToolManagers.conservationTools.CCMPredManger import CCMPredManager as CorrMutManager, HHBLITS_CMD_TEMPLATE

from .seqToolManagers.asaPredictions.spider2Manager import Spider2Manager

WIN_SIZE= 11

class SeqFeatComputer(FeaturesComputer):
  '''
  Intended to be used for computing all sequential features for one complex
  '''
  CLEAN_SEQ_REGEX= re.compile(r"\n|\s|\r")
  
  def __init__(self, prefix, computedFeatsRootDir= None, statusManager=None):
    '''
      @prefix. An id for a complex. Pe. 1A2K. Extended prefix in seq step is prefix+chainType+chainId
      :param computedFeatsRootDir: str. path where features will be stored
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''
    FeaturesComputer.__init__(self, prefix=prefix, computedFeatsRootDir=computedFeatsRootDir, statusManager=statusManager)
    self.computedFeatsRootDir= myMakeDir(self.computedFeatsRootDir, "seqStep")

    self.SeqStructMapper= SeqStructMapper( computedFeatsRootDir= self.computedFeatsRootDir, 
                                   statusManager= self.statusManager)
    self.seqsWorkingDir= self.SeqStructMapper.getSeqsOutDir()
    self.seqWindower= WindowSeqManager( self.computedFeatsRootDir, winSize= WIN_SIZE) #Object to launch and get sliding window of seqs
    self.conservationOutPath= myMakeDir(self.computedFeatsRootDir, "conservation")
    self.psiBlastManager= PsiBlastManager( self.conservationOutPath, winSize= WIN_SIZE) #Object to launch and process blast

    self.al2coManager= Al2coManager( self.conservationOutPath, winSize= None) #Object to launch and process al2co
    
    self.spider2Manager= Spider2Manager( self.computedFeatsRootDir, winSize= None) #Object to launch and process spyder2
    self.singleChainOutPaths=[ self.seqsWorkingDir, self.psiBlastManager.getFinalPath(), self.al2coManager.getFinalPath(),
                              self.seqWindower.getFinalPath()] + [self.spider2Manager.getFinalPath()]
    self.pairwiseOutPaths=[]
    if self.useCorrMut:
      self.hHBlitsManager= HHBlitsManager( self.conservationOutPath, winSize= None, 
                                           hhBlitsCMD_template=HHBLITS_CMD_TEMPLATE) #Object to launch and process hhblits
      self.singleChainOutPaths.append( self.hHBlitsManager.getFinalPath())
      self.corrMutManager= CorrMutManager( self.conservationOutPath) #Object to launch and process psicov
      self.pairwiseOutPaths.append( self.corrMutManager.getFinalPath())
      
    self.singleChainOutPaths= filter(None, self.singleChainOutPaths)
    self.pairwiseOutPaths= filter(None, self.pairwiseOutPaths)

  def computeComplex(self, fnameL, fnameR, structureL, structureR, lPdbId=None, rPdbId=None, areLRequivalentProteins=False):
    '''
    Applies self.computeOneFile method to both ligand and receptor pdbFiles.

    :param fnameL: str. fname or structure to pdbfile or fname to fasta of ligand
    :param fnameR: str. fname or structure to pdbfile or fname to fasta of receptor
    :param structureL: Bio.PDB.Structure.Structure. Structure of ligand protein (contained in fnameL). None if fasta
    :param structureR: Bio.PDB.Structure.Structure. Structure of receptor protein (contained in fnameR). None if fasta
    
    :param lPdbId: str. pdbId of receptor in order to query 3dCons. If None, psiblast will be launched directly
    :param rPdbId: str. pdbId of receptor in order to query 3dCons. If None, psiblast will be launched directly

    '''

    self.SeqStructMapper.computeComplex(fnameL, fnameR, structureL, structureR)
    if areLRequivalentProteins:
      fnameR,structureR = None, None
    computedSeqs= self.getAlreadyComputedSeqs()
    if fnameL and fnameR:
      n_seqs= len(list(self.SeqStructMapper.enumSeqs("l"))) + len(list(self.SeqStructMapper.enumSeqs("r")))
      HHblitsAligsDict= {"l":{}, "r":{}}
      chainTypes= ["l","r"]
    else:
      n_seqs= len(list(self.SeqStructMapper.enumSeqs("l")))
      HHblitsAligsDict= {"l":{}}
      chainTypes= ["l"]
    iterNum=0
    for chainType in chainTypes:
      for chainType, chainId in  self.SeqStructMapper.enumSeqs(chainType):
        seqStr, fastaFname= self.SeqStructMapper.getSeq(chainType, chainId)
        extendedPrefix= self.getExtendedPrefix( fastaFname) #Extended prefix in seq step is prefix+chainType+chainId
        iterNum+=1
        self.reportStatus("\n... Computing sequence features for partner %d chain %s (%d/%d)"%(1 if chainType=="l" else 2, 
                                              chainId, iterNum, n_seqs) )
        if seqStr in computedSeqs:
          if extendedPrefix in computedSeqs[seqStr]:
            print("seq %s already computed"%(extendedPrefix))
          else:
            computedExtenPrefix= computedSeqs[seqStr].pop()
            print("seq %s is the same than %s"%(extendedPrefix, computedExtenPrefix))
            computedSeqs[seqStr].add(computedExtenPrefix)
            self.copySameSeq( computedExtenPrefix, extendedPrefix) 
        pdbCode= lPdbId if chainType=="l" else rPdbId
        psiblastOutName, pssmOutNameRaw = self.psiBlastManager.computeFromSeqStructMapper( self.SeqStructMapper, 
                                                                                       extendedPrefix, pdbCode)
        self.seqWindower.computeFromSeqStructMapper( self.SeqStructMapper, extendedPrefix )
        self.al2coManager.computeFromSeqStructMapper( self.SeqStructMapper, extendedPrefix, psiblastOutName, pssmOutNameRaw )
        self.spider2Manager.computeFromSeqStructMapper( self.SeqStructMapper, extendedPrefix, pssmOutNameRaw )
        if self.useCorrMut:
          aligsName, __, __ =self.hHBlitsManager.computeFromSeqStructMapper(self.SeqStructMapper, extendedPrefix) # Must be executed after psiBlastManager execution  
          HHblitsAligsDict[chainType][chainId]= aligsName
        self.reportStatus("..... Done" )

    if self.useCorrMut:    
      self.reportStatus("\n Computing correlated mutations" )
      self.corrMutManager.computeFromSeqStructMapper(self.SeqStructMapper, extendedPrefix, HHblitsAligsDict)

    #compress raw files
    for chainType in chainTypes:
      for chainType, chainId in  self.SeqStructMapper.enumSeqs(chainType):
        seqStr, fastaFname= self.SeqStructMapper.getSeq(chainType, chainId)
        extendedPrefix= self.getExtendedPrefix( fastaFname) #Extended prefix in seq step is prefix+chainType+chainId
        print(extendedPrefix)
        self.psiBlastManager.compressRawData(extendedPrefix)
        if self.useCorrMut:    self.hHBlitsManager.compressRawData(extendedPrefix)
    self.reportStatus("..... Done" )

  def getAlreadyComputedSeqs(self):
    '''
      Gets the sequences that have been already computed.
      return alreadyComputedSeqs: dict
                {"seq": set([unambiguous_id_i])
                e.g.:
                {"CTARRTCAAYALLPYV": set(["1A2K_r_A_T","1A2K_r_B_T1"])}
    '''
    extendedPrefixList_of_Sets=[]
    for path in self.singleChainOutPaths:
      extendedPrefixList_of_Sets.append( set([  self.getExtendedPrefix(elem) for elem in os.listdir(  path ) ]))
      
    prefixesFullyDone= reduce((lambda x,y: x.intersection(y)), extendedPrefixList_of_Sets)
    
    if self.useCorrMut:      
      pairwisePrefixList_of_Sets=[]
      for path in self.pairwiseOutPaths:
        pairwisePrefixList_of_Sets.append( set([  self.getExtendedPrefix(elem).split("_")[0]  for elem in os.listdir(  path ) ]))
      
      corrMutPrefixes= reduce((lambda x,y: x.intersection(y)), extendedPrefixList_of_Sets)
    alreadyComputedSeqs={}    
    for seqFname in os.listdir(  self.seqsWorkingDir ):
      extendedPrefix= self.getExtendedPrefix(seqFname)
      regularPrefix= extendedPrefix.split("_")[0]
      if not extendedPrefix in prefixesFullyDone or (self.useCorrMut and not regularPrefix in corrMutPrefixes): 
        continue
      if seqFname.endswith(".gz"):
        openFun= gzip.open
      else:
        openFun= open
      with openFun( os.path.join(self.seqsWorkingDir, seqFname)) as f:
        f.readline()
        seq= f.readline().strip()
        if seq in alreadyComputedSeqs:
          alreadyComputedSeqs[seq].add( extendedPrefix )
        else:
          alreadyComputedSeqs[seq]= set([ extendedPrefix ])
          
    return alreadyComputedSeqs
  
  
  def copySameSeq(self, oldExtendedPrefix, newExtendedPrefix):
    '''
      Given an unambiguous id that is known to be already computed (oldExtendedPrefix), and a non-computed
      unambiguous id (newExtendedPrefix), copy results files from oldExtendedPrefix to newExtendedPrefix
      as their sequences are the same.
      :param oldExtendedPrefix: unambiguous id of a sequence that has already been computed
      :param newExtendedPrefix: unambiguous id of a sequence that has not been computed yet and is equal
                                to the sequence of oldExtendedPrefix
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''  
    fnamesToCopy=[]
    oldPrefix, oldChainType, oldChainId = self.splitExtendedPrefix(oldExtendedPrefix)[:3]
    newPrefix, newChainType, newChainId = self.splitExtendedPrefix(newExtendedPrefix)[:3]
    if oldPrefix==newPrefix and oldChainType== newChainType and oldChainId== newChainId:
      fnamesToCopy= self.psiBlastManager.getFNames(oldExtendedPrefix)
      fnamesToCopy+= list(self.al2coManager.getFNames(oldExtendedPrefix))
      if self.useCorrMut:
        fnamesToCopy+= list(self.hHBlitsManager.getFNames(oldExtendedPrefix))

    for oriName in fnamesToCopy:
      destName= oriName.replace(oldExtendedPrefix, newExtendedPrefix)
      tryToSymlink(oriName, destName)
      
    if oldPrefix==newPrefix and oldChainType== newChainType and oldChainId== newChainId and self.useCorrMut:
      fnamesToCopy= list(self.corrMutManager.getFNames(oldExtendedPrefix))
      for oriName in fnamesToCopy:
        destName= oriName.replace(oldPrefix, newPrefix)
        tryToSymlink(oriName, destName)

        
 
