from __future__ import absolute_import

import sys, os
import subprocess

from ..SeqFeatComputer import SeqFeatComputer, FeatureComputerException
from Config import Configuration

FILTER_OUT_LABELS=False  
class seqToolManager(Configuration):
  '''
  Abastract class that will be used to implement psiblast and al2co managers
  '''
  def __init__(self, seqsManager, outPath, winSize):
    '''
      @param seqsManager: ..manageSeqs.seqsManager.SeqsManager 
      @param outPath: str. root path where psiblast and al2co scores will be saved
      @param winSize: int. The size of sliding window 
    '''    
    Configuration.__init__(self)

    self.seqsManager= seqsManager
    self.seqsWorkingDir= self.seqsManager.getSeqsOutDir()
    self.outPath= outPath
    self.winSize= winSize
    self.filterOutLabels= FILTER_OUT_LABELS
  def getFNames(self, prefixExtended):
    '''
    abstract method
    Returns a dict that contains the fnames that will be used by psiblast and other tools
    @param prefixExtended. prefix for output fnames. They are form as follows: prefix+chainType+chainId
    @return Dict   { featName: (fname1,fnam2...) ...}
    '''
    return None

  def checkAlreayComputed(self, prefixExtended):
    '''
    checks if output files have been computed
    @param prefixExtended. prefix for output fnames. They are form as follows: prefix+chainType+chainId
    @return boolean. True if prefixExtended sequences have been already computed, False otherwise
    '''  
    dictOfNames= self.getFNames(prefixExtended)
    for fnameTuple in dictOfNames.values():
      for fname in fnameTuple:
#        print( "Checking fname: %s"%fname)
        if not os.path.isfile(fname):
          return False
        else:
          nlines=int(subprocess.check_output('wc -l {}'.format(fname), shell=True).split()[0])
          if nlines<11: #Too small file to be correct
            return False
    return True

if __name__=="__main__":
  print("Done")
