# -*- coding: utf-8 -*-

from __future__ import absolute_import, print_function
import os
from ...seqToolManager import SeqToolManager
from utils import myMakeDir  #utils is at the root of the package

class WindowSeqManager(SeqToolManager):
  '''
    Computes sliding window from amino acids
  '''
  NO_AA_LETTER='Z'
  def __init__(self, computedFeatsRootDir, winSize, statusManager=None):
    '''
      :param computedFeatsRootDir: str. root path where results will be saved
      :param winSize: int>=1 or None. The size of the windows for sliding window if desired
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''
    assert not winSize is None and winSize>=1, "Error, winSize must be greater or equal to one for WindowSeqManager"
    SeqToolManager.__init__(self, computedFeatsRootDir, winSize)
    self.outPath= myMakeDir(self.computedFeatsRootDir,"slidingWinSeq"+str(winSize))

    
  def getFinalPath(self):
    '''
      returns paths where final results are saved
      :return outPath: str
    '''
    return self.outPath
    
  def getFNames(self, prefixExtended):
    '''
      Returns a lsit that contains the fnames that will be used  by WindowSeqManager
      :param prefixExtended. prefix for output fnames.
      :return List of fname
    '''
    windowedOutName= os.path.join(self.outPath, prefixExtended+".wsize"+str(self.winSize)+".seq.gz")
    return [windowedOutName]
    
    
  def computeFromSeqStructMapper(self, seqStructMap, prefixExtended):
    '''
      Computes sliding window for the sequence seqStr contained in seqStructMap under prefixExtended
      unambiguous id
      :param seqStructMap: computeFeatures.seqStep.seqToolManagers.seqExtraction.SeqStructMapper
      :param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
    '''
    __, chainType, chainId= self.splitExtendedPrefix(prefixExtended)[:3]
    seqStr, fastaFname= seqStructMap.getSeq(chainType, chainId) # repeat as psiBlastManager can modify seqs     
    seqStructMap.setCurrentSeq(seqStr, chainType, chainId)
    if self.checkAlreayComputed(prefixExtended):
      print("winSeq already computed for %s"%prefixExtended)
      return 0
    fNames= self.getFNames(prefixExtended)
    print("launching winSeq over %s"%prefixExtended)      
    
    dataList= []
    for i, letter in enumerate(seqStr):
      structIndex= seqStructMap.seqToStructIndex(chainType, chainId, i, asString= True)
      if structIndex:
        if self.filterOutLabels and structIndex[-1].isalpha():
          continue
      else:
        structIndex=str(i)+"?"
      full_resId_tuple= (chainId, structIndex, letter)
      dataList.append( (full_resId_tuple, ([letter])) )
    self.makeWindowed( dataList, ["aa"], [WindowSeqManager.NO_AA_LETTER], [type(self).AA_CODE_ELEMENTS], 
                         outName= fNames[0])
      
  
