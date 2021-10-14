from __future__ import absolute_import, print_function
import  os
from subprocess import Popen, PIPE
from utils import myMakeDir
import time
import random
import GPUtil
import numpy as np

from .corrMutGeneric import CorrMutGeneric


HHBLITS_CMD_TEMPLATE=("%(hhBlitsBinPath)s/hhblits -i %(fastaInFname)s -d %(hhblitsDB)s -maxmem 8 "+
                "-oa3m %(aligsName)s -cpu %(hhBlitsNThrs)d -ohhm %(profileNameRaw)s -maxfilt 100000 "+
                " -realign_max 100000 -all -B 100000 -Z 100000 -n 3 -e 0.001  -o /dev/null && "+
                "%(hhBlitsBinPath)s/hhfilter -id 90 -i %(aligsName)s -o %(aligsName)s")
                
class CCMPredManager(CorrMutGeneric):
  '''
    Computes corrMut and processes their outputs. Extends class CorrMutGeneric
  '''
  def __init__(self, computedFeatsRootDir, statusManager=None):
    '''
      :param computedFeatsRootDir: str. path where corrMut scores will be saved
    '''
    CorrMutGeneric.__init__(self, computedFeatsRootDir, statusManager)
    
    self.corrMutOutPath= myMakeDir(computedFeatsRootDir,"corrMut")
    self.featName="ccmPred"
 
  def getAligsDims(self, aligFormatedName):
    ncols=None
    nrows=0
    with open(aligFormatedName) as f:
      for line in f:
        nrows+=1
        ncols_tmp= len(line.strip())
        if not ncols is None:
          assert ncols_tmp == ncols, "Bad formated alig %s"%aligFormatedName
        else:
          ncols=ncols_tmp
    return nrows, ncols
    
  def getTotalGPUMemory(self, gpuNumber):
    time.sleep(random.random())
    print("Getting total gpuMemory for %s"%gpuNumber)
    if gpuNumber==None:
      return 0
    try:
      process= Popen(["nvidia-smi", "--query-gpu=memory.free", "--format=csv"], stdout=PIPE, stderr=PIPE)
      processOut= process.communicate()
      if processOut[1]!="" :
        return 0
      else:
        array_lines= processOut[0].split("\n")
        freeMem, memUnit= array_lines[1+gpuNumber].split()
        freeMem= float(freeMem)
        if memUnit=="GiB":
          freeMem*= 2**30
        elif memUnit=="MiB":
          freeMem*= 2**20
        elif memUnit=="KiB":
          freeMem*= 2**10
        return freeMem
    except Exception:
      return 0
     
  def lauchCorrMutProgram(self, aligFormatedName, ignoreGPU=False):
    nrows, ncols= self.getAligsDims( aligFormatedName)
    memoryRequired= 4*(4*(ncols*ncols*32*21 + ncols*20) + 23*nrows*ncols + nrows + ncols*ncols) + 2*nrows*ncols + 1024
    tmpResults= os.path.basename(aligFormatedName)
    tmpResults= os.path.join(self.tmp, tmpResults)
    
    try:
      gpuNumber= GPUtil.getFirstAvailable(order = 'first', maxLoad=0.3, maxMemory=0.3, attempts=2, interval=3)[0]
    except (RuntimeError, OSError, ValueError):
      gpuNumber=None      

    wasRunOnGPU=False
    if not ignoreGPU and memoryRequired*1.1 < self.getTotalGPUMemory( gpuNumber): #*1.1 as a margin of tolerance 
      cmdArray=[self.ccmPredBin,"-R", "-d", str(gpuNumber), aligFormatedName, tmpResults ]
      wasRunOnGPU= True
    else:
      cmdArray=[self.ccmPredBin,"-R", "-t", str(self.corrMutNThrs), aligFormatedName, tmpResults ]
    print(" ".join(cmdArray))
    process= Popen(cmdArray, stdout=PIPE, stderr=PIPE)
    processOut= process.communicate()
    try:
      iterOfCorrelatedRows= self.processOutput(processOut, tmpResults)
    except ValueError as e:
      print(e)
      iterOfCorrelatedRows=None
    
    if iterOfCorrelatedRows is None and wasRunOnGPU==True:
      print("Error running ccmpred on gpu, trying cpu")
      iterOfCorrelatedRows= self.lauchCorrMutProgram(aligFormatedName, ignoreGPU=True)
    return iterOfCorrelatedRows
          
  def processOutput(self, processOut, tmpResults):
    '''
      :param processOut (stdin, stderr)
      :param tmpResults: str. Path to temporal results
    '''
    
    stdout_list= processOut[0].split("\n")
    if (len(stdout_list)<=1 or stdout_list[-7]!="Done with status code 0 - Success!"):
      raise ValueError("Error processing data ccmpred:", processOut)
    return self.processOutputHelper( processOut, tmpResults)

#    results= np.loadtxt(tmpResults)
#    for i in range(results.shape[0]):
#      for j in range(results.shape[1]):
#        yield  [i,j,results[i,j]]
    return self.processOutputHelper( processOut, tmpResults)

  def processOutputHelper(self, processOut, tmpResults):
    results= np.loadtxt(tmpResults)
    for i in range(results.shape[0]):
      for j in range(results.shape[1]):
        yield  [i,j,results[i,j]]
