from __future__ import absolute_import, print_function
import sys, os
from subprocess import Popen, PIPE, check_output
from Config import Configuration
from utils import myMakeDir, tryToRemove #utils is at the root of the package
import re
import time

from ..seqToolManager import  FeatureComputerException
from .corrMutGeneric import CorrMutGeneric
from utils import myMakeDir, tryToRemove #utils is at the root of the package

HHBLITS_CMD_TEMPLATE=("%(hhblitsBin)s -i %(fastaInFname)s -n 4 -d %(hhblitsDB)s "+
                "-oa3m %(aligsName)s -cpu %(psiBlastNThrs)d -ohhm %(profileNameRaw)s")
class PsicovManager(CorrMutGeneric):
  '''
    Computes corrMut and processes their outputs. Extends class CorrMutGeneric
  '''
  def __init__(self, seqsManager, outPath):
    '''
      @param seqsManager: ..manageSeqs.seqsManager.SeqsManager 
      @param outPath: str. path where corrMut scores will be saved
    '''  
    CorrMutGeneric.__init__(self,seqsManager, outPath)
    
    self.seqsManager= seqsManager
    self.corrMutOutPath= myMakeDir(outPath,"corrMut")
    self.featName="psicov"
 
  def lauchCorrMutProgram(self, aligFormatedName):
    cmdArray=[self.psicovBin,"-p", "-z", str(self.psiBlastNThrs),"-d", "0.03", "-r", "0.001","-o", 
              "-j", "0", aligFormatedName ]
    print(" ".join(cmdArray))
    process= Popen(cmdArray, stdout=PIPE, stderr=PIPE)
    processOut= process.communicate()
    try:
      iterOfCorrelatedRows= self.processOutput(processOut)
    except ValueError:
      iterOfCorrelatedRows=None
    return iterOfCorrelatedRows
          
  def processOutput(self, processOut):
    '''
      @param processOut (stdin, stderr)
    '''
    if len(processOut[1])>0 or  processOut[0]=="" or "*** Sorry" in processOut[0]: #Error happend
      print(processOut)
      raise ValueError("Error processing data psicov")
    for line in processOut[0].split("\n")[1:]:
      if len(line)==0: break
      i, j, score= line.split()
      i, j=int(i)-1, int(j)-1
      yield [i,j, float(score)]
