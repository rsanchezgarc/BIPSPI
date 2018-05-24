from __future__ import absolute_import
import os,sys
from subprocess import Popen, PIPE

from .manageSeqs.seqsManager import SeqsManager
from .seqToolManagers.conservationTools.PsiBlastManager import PsiBlastManager
from .seqToolManagers.conservationTools.Al2coManager import Al2coManager

from .seqToolManagers.conservationTools.HHblitsManager import HHBlitsManager
#from .seqToolManagers.conservationTools.PSICOVManager import PsicovManager as CorrMutManager, HHBLITS_CMD_TEMPLATE
from .seqToolManagers.conservationTools.CCMPredManger import CCMPredManager as CorrMutManager, HHBLITS_CMD_TEMPLATE
from .seqToolManagers.ASA_SS_PREDS.spider2Manager import Spider2Manager

from .SeqFeatComputer import SeqFeatComputer, FeatureComputerException
from utils import myMakeDir, myMakeDirUnique, tryToRemove, tryToSymlink #utils is at the root of the package

class SeqFeaturesCalculator(SeqFeatComputer):
  '''
  Extends SeqFeatComputer class. Computes PSSMs, al2co and possibly correlated mutations. 
  '''
  def __init__(self, rFname, lFname, computedFeatsRootDir= None, winSize= 11, rPdbId=None, lPdbId=None, 
                     useCorrMut=False, statusManager=None):
    '''
      @param rFname: str. path to receptor pdb file or fasta file
      @param lFname: str. path to ligand pdb file or fasta file
      @param computedFeatsRootDir: str. path where features will be stored. If None they will be stored
                                        at default path (assigned in ../Config.py)
      @param winSize: int. The size of the window that will be computed with the sequential features.
      @param rPdbId: str. pdbId of receptor in order to query 3dCons. If None, psiblast will be launched directly
      @param lPdbId: str. pdbId of receptor in order to query 3dCons. If None, psiblast will be launched directly
      @param useCorrMut: boolean. True if corrMuts wanted to be used
      @param statusManager: class that implements .setStatus(msg) to communicate
    '''
    
    SeqFeatComputer.__init__(self, rFname, lFname, computedFeatsRootDir, statusManager=statusManager)
    '''
        computedFeatsRootDir=/path/to...computedFeatsRootDir/
       at this point self.computedFeatsRootDir=/path/to...computedFeatsRootDir/seqStep
    '''

    conservationOutPath= myMakeDir(self.computedFeatsRootDir, "conservation")    
    self.seqsManager= SeqsManager( rFname, lFname, computedFeatsRootDir)
    self.seqsManager.computeFun() #Extracts sequences from pdb files
    self.seqsWorkingDir= self.seqsManager.getSeqsOutDir()
    self.useCorrMut= useCorrMut
    self.lPdbId= lPdbId
    self.rPdbId= rPdbId

    self.psiBlastManager= PsiBlastManager( self.seqsManager, conservationOutPath, winSize) #Object to launch and process blast
    self.al2coManager= Al2coManager( self.seqsManager, conservationOutPath, None) #Object to launch and process al2co
    
    if self.useCorrMut:
      self.hHBlitsManager= HHBlitsManager( self.seqsManager, conservationOutPath, winSize, HHBLITS_CMD_TEMPLATE) #Object to launch and process hhblits
      self.corrMutManager= CorrMutManager( self.seqsManager, conservationOutPath) #Object to launch and process psicov
    
    self.spider2Manager= Spider2Manager( self.seqsManager, self.computedFeatsRootDir, None) #Object to launch and process spider2
    
  def getAlreadyConputedSeqs(self):
    '''
      Gets the sequences that have been already computed.
      return alreadyComputedSeqs: dict
                {"seq": set([unambiguous_id_i])
                e.g.:
                {"CTARRTCAAYALLPYV": set(["1A2K_r_A_u","1A2K_r_B_u"])}
    '''  
    extendedPrefixesPssm=  set([  self.getExtendedPrefix(elem) for elem in os.listdir(  self.psiBlastManager.getFinalPath()) ])
    extendedPrefixesAl2Co= set([  self.getExtendedPrefix(elem) for elem in  os.listdir(  self.al2coManager.getFinalPath()) ])
    extendedPrefixesSpider2= set([  self.getExtendedPrefix(elem) for elem in  os.listdir(  self.spider2Manager.getFinalPath()) ])    

    extendedPrefixesSeqs=  set([  self.getExtendedPrefix(elem) for elem in os.listdir(  self.seqsWorkingDir ) ])

    prefixesFullyDone= extendedPrefixesPssm.intersection(extendedPrefixesAl2Co).intersection(extendedPrefixesSeqs)
    prefixesFullyDone= prefixesFullyDone.intersection(extendedPrefixesSpider2)
    if self.useCorrMut:
      extendedPrefixesHHblits= set([  self.getExtendedPrefix(elem) for elem in  os.listdir(  self.hHBlitsManager.getFinalPath()) ]) 
      corrMutPrefix= set([ elem.split(".")[0].split("_")[0] for elem in  os.listdir( self.corrMutManager.getFinalPath()) ])
      prefixesFullyDone= prefixesFullyDone.intersection(extendedPrefixesHHblits)
    
    alreadyComputedSeqs={}    
    for seqFname in os.listdir(  self.seqsWorkingDir ):
      extendedPrefix= self.getExtendedPrefix(seqFname)
      regularPrefix= extendedPrefix.split("_")[0]
      if not extendedPrefix in prefixesFullyDone or (self.useCorrMut and not regularPrefix in corrMutPrefix): 
        continue
      with open( os.path.join(self.seqsWorkingDir, seqFname)) as f:
        f.readline()
        seq= f.readline().strip()
        if seq in alreadyComputedSeqs:
          alreadyComputedSeqs[seq].add( extendedPrefix )
        else:
          alreadyComputedSeqs[seq]= set([ extendedPrefix ])
          
    return alreadyComputedSeqs
        
  def copySameSeq(self, oldExtenPrefix, newExtendedPrefix):
    '''
      Given an unambiguous id that is known to be already computed (oldExtenPrefix), and a non-computed
      unambiguous id (newExtendedPrefix), copy results files from oldExtenPrefix to newExtendedPrefix
      as their sequences are the same.
      @param oldExtenPrefix: unambiguous id of a sequence that has already been computed
      @param newExtendedPrefix: unambiguous id of a sequence that has not been computed yet and is equal
                                to the sequence of oldExtenPrefix
    '''  
    dictOfNames= self.psiBlastManager.getFNames(oldExtenPrefix)
    fnamesToCopy= [dictOfNames["psiblast"][0], dictOfNames["pssm"][0] ] #Do not copy processed results as they are edited for each prefix
    for oriName in fnamesToCopy:
      destName= oriName.replace(oldExtenPrefix,newExtendedPrefix)
      tryToSymlink(oriName, destName)

  def computeOneComplex(self):
    '''
      Computes conservation features for a given complex
      Overrides computeOneComplex of parent class SeqFeatComputer
    '''
    computedSeqs= self.getAlreadyConputedSeqs()
    n_seqs= len(list(self.seqsManager.enumSeqs("l"))) + len(list(self.seqsManager.enumSeqs("r")))
    iterNum=0
    HHblitsAligsDict= {"l":{}, "r":{}}
    for chainType in ["l","r"]:
      allSeqsInvolved={}
      for chainType, chainId in  self.seqsManager.enumSeqs(chainType):
        seqStr, fastaFname= self.seqsManager.getSeq(chainType, chainId)
        extendedPrefix= self.getExtendedPrefix( fastaFname)
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
        pdbCode= self.lPdbId if chainType=="l" else self.rPdbId
        psiblastOutName, pssmOutNameRaw = self.psiBlastManager.compute( extendedPrefix, pdbCode) #Uses 3dcons if possible
        self.al2coManager.compute(extendedPrefix, psiblastOutName, pssmOutNameRaw ) # Must be executed after psiBlastManager execution  
        self.spider2Manager.compute(extendedPrefix, pssmOutNameRaw )  # Must be executed after psiBlastManager execution
        if self.useCorrMut:
          aligsName, __, __ =self.hHBlitsManager.compute(extendedPrefix) # Must be executed after psiBlastManager execution  
          HHblitsAligsDict[chainType][chainId]= aligsName
        self.reportStatus("..... Done" )
        
    if self.useCorrMut:
      self.reportStatus("\n Computing correlated mutations. It may take a while" )
      self.corrMutManager.compute(HHblitsAligsDict, extendedPrefix.split("_")[0])
      self.reportStatus("..... Done" )
    
  
if __name__=="__main__":
  pdbsIndir= "/home/rsanchez/Tesis/rriPredMethod/data/develData/pdbFiles"
  computedFeatsRootDir= "/home/rsanchez/Tesis/rriPredMethod/data/develData/computedFeatures"
  SeqFeatComputer.computeFeaturesAllComplexes(SeqFeaturesCalculator,pdbsIndir= pdbsIndir ,
                                              computedFeatsRootDir= computedFeatsRootDir, ncpu=1 )
