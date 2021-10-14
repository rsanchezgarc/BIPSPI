from __future__ import absolute_import, print_function
from subprocess import Popen, PIPE
from utils import myMakeDir #utils is at the root of the package
from .corrMutGeneric import CorrMutGeneric

HHBLITS_CMD_TEMPLATE=("%(hhBlitsBinPath)s/hhblits -i %(fastaInFname)s -n 4 -d %(hhblitsDB)s "+
                "-oa3m %(aligsName)s -cpu %(hhBlitsNThrs)d -ohhm %(profileNameRaw)s -o /dev/null")
class PsicovManager(CorrMutGeneric):
  '''
    Computes corrMut and processes their outputs. Extends class CorrMutGeneric
  '''
  def __init__(self, computedFeatsRootDir):
    '''
      :param outPath: str. path where corrMut scores will be saved
    '''  
    CorrMutGeneric.__init__(self, computedFeatsRootDir)
    self.corrMutOutPath= myMakeDir(computedFeatsRootDir,"corrMut")
    self.featName="psicov"
 
  def lauchCorrMutProgram(self, aligFormatedName):
    cmdArray=[self.psicovBin,"-p", "-z", str(self.corrMutNThrs),"-d", "0.03", "-r", "0.001","-o", 
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
      :param processOut (stdin, stderr)
    '''
    if len(processOut[1])>0 or  processOut[0]=="" or "*** Sorry" in processOut[0]: #Error happend
      print(processOut)
      raise ValueError("Error processing data psicov")
    for line in processOut[0].split("\n")[1:]:
      if len(line)==0: break
      i, j, score= line.split()
      i, j=int(i)-1, int(j)-1
      yield [i,j, float(score)]
