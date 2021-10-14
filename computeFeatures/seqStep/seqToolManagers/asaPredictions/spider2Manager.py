from __future__ import absolute_import, print_function
import os
import numpy as np
from subprocess import Popen, PIPE
from Bio.PDB.Polypeptide import aa1 as AA_STANDARD
from ...seqToolManager import SeqToolManager
from utils import myMakeDir, tryToRemove #utils is at the root of the package

class Spider2Manager(SeqToolManager):
  '''
    Computes asa and secondary structure predictions and processes their outputs. 
    Extends class seqToolManager
  '''
  VAR_LIST= ["pred_asa", "P_C", "P_E", "P_H"]
  BAD_SCORE_PREDS = ["-1048576",  "-1.0",  "-1.0", "-1.0" ] #Something went wrong tag
  def __init__(self, computedFeatsRootDir, winSize=None, statusManager=None):
    '''
      :param computedFeatsRootDir: str. root path where results will be saved
      :param winSize: int>=1 or None. The size of the windows for sliding window if desired
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''
    SeqToolManager.__init__(self, computedFeatsRootDir, winSize)
    # self.spider2PyScript inherited from Config.py
    self.spider2OutPath= myMakeDir(self.computedFeatsRootDir,"SPIDER2_predAsaSS")
    if winSize:
      self.spider2PathWindowed= myMakeDir(self.computedFeatsRootDir,"spyder2_wSize"+str(winSize))
    else:
      self.spider2PathWindowed= None
    
  def getFinalPath(self):
    '''
      returns path where results are saved
    '''
    return self.spider2OutPath
    
  def getFNames(self, prefixExtended):
    '''
      Returns a dict that contains the fnames that will be used  by spyder2
      :param prefixExtended. prefix for output fnames.
      :return list of fnames: [ fname1, fnam2, ...]
    '''    
    spyderProc= os.path.join(self.spider2OutPath, prefixExtended+".spyder2.gz")
    fNames=[spyderProc]
    if not self.winSize is None:
      spyder2WindowedOutName= os.path.join(self.spider2PathWindowed, prefixExtended+".wsize"+str(self.winSize)+".spyder2.gz")
      fNames+= [spyder2WindowedOutName]
    return fNames
    
    
  def computeFromSeqStructMapper(self, seqStructMap, prefixExtended, pssmOutNameRaw):
    '''
      Computes spyder2 for the sequence seqStr, that is contained at fastaInFname. This sequence is
      associated with prefixExtended as an unambiguous id
      :param seqStructMap: computeFeatures.seqStep.seqToolManagers.seqExtraction.SeqStructMapper
      :param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      :param pssmOutNameRaw: str. Path to psiblast pssms results
    '''
    if pssmOutNameRaw is not None:
      if not os.path.isfile(pssmOutNameRaw):
        pssmOutNameRaw= pssmOutNameRaw+".gz"
      if os.path.isfile(pssmOutNameRaw):
        uncompressFileName= self.uncompressFile(pssmOutNameRaw, self.tmp)
      else:
        uncompressFileName = None
    else:
      uncompressFileName= None
    try:
      prefix, chainType, chainId= self.splitExtendedPrefix(prefixExtended)[:3]
      seqStr, fastaFname= seqStructMap.getSeq(chainType, chainId) # repeat as psiBlastManager can modify seqs
      seqStructMap.setCurrentSeq(seqStr, chainType, chainId)
      if self.checkAlreayComputed(prefixExtended):
        print("spyder2 already computed for %s"%prefixExtended)
        return 0
      fNames= self.getFNames(prefixExtended)
      spider2ProcName= fNames[0]
      spider2RawName= os.path.join(self.spider2OutPath, prefixExtended+".spd3")
      print("launching spyder2 over %s"%prefixExtended)
      curWd= os.getcwd()
      os.chdir(self.spider2OutPath)
      if uncompressFileName is not None:
        cmd= ["python", self.spider2PyScript, uncompressFileName]
        process= Popen(cmd, stdout=PIPE, stderr=PIPE)
        processOut= process.communicate()
        os.chdir(curWd)
        if len(processOut[1])>0:
          print("Error computing spider2. Caught stdin/stderr:\n",processOut[0],processOut[1])
      else:
        spider2RawName= None
      dataList= self.processSpider2(seqStr, seqStructMap, prefixExtended, spider2RawName, spider2ProcName)
          
      if self.winSize:
        self.makeWindowed( dataList, ["asa", "P_C", "P_E", "P_H"], Spider2Manager.BAD_SCORE_PREDS, [None]*4,
                           fNames[1])
    except (Exception, KeyboardInterrupt):
      self.tryToRemoveAllFnames(prefixExtended)
      raise
    finally:
      if uncompressFileName is not None:
        tryToRemove(uncompressFileName)
  
  def processSpider2(self, seq, seqStructMap, prefixExtended, spider2Raw, spider2Proc):
    '''
      Reads spider2 output file and writes another one with tabulated format, headers and
      some error checking.
      :param: seq: str. Sequence of the chain
      :param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      :param spider2Raw: str. Path to spyder2 results
      :param spider2Proc: str. Path where formated results will be saved.
    '''

    if spider2Raw is not None:
      try:
        predictionsData = self.loadSpider2(spider2Raw)
      except IOError:
        predictionsData= [ (letter, type(self).BAD_SCORE_PREDS ) for letter in seq]
    else:
      predictionsData = [(letter, type(self).BAD_SCORE_PREDS) for letter in seq]

    prefix, chainType, chainId= self.splitExtendedPrefix(prefixExtended)[:3]
    try:
      spyderIx=0
      seqIx=0
      seqLen= len(seq)
      spyderLen= len(predictionsData)
      dataList=[]
      listOfRowsToPrint=[]

      while seqIx<seqLen and spyderIx<spyderLen:
        letter= seq[seqIx]
        letterSpider, predsVal= predictionsData[spyderIx]
        if letterSpider== letter:
          structIndex= seqStructMap.seqToStructIndex(chainType, chainId, seqIx, asString= True)
          if structIndex:
            if self.filterOutLabels and structIndex[-1].isalpha():
              continue
          else:
            structIndex=str(seqIx)+"?"
          dataList.append( ( (chainId, structIndex,letter), ( [predsVal[0]], predsVal[1:]) ) )
          listOfRowsToPrint.append( "%s %s %s %s %s %s %s"%(( chainId, structIndex, letter)+tuple(predsVal)) )
          spyderIx+=1
          seqIx+=1
        elif not letter in AA_STANDARD and letterSpider=="-":
          spyderIx+=1
          seqIx+=1
        elif letterSpider=="-":
          spyderIx+=1
        else:
          print(predictionsData)
          print(spyderIx, seqIx)
          raise ValueError("Spyder mismatch %s %s "%(letterSpider, letter))

      self.writeResultsFromDataDictSingleChain( {chainId: listOfRowsToPrint }, outName= spider2Proc)
      return dataList
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing %s"%spider2Proc)
      tryToRemove(spider2Proc)
      raise
    finally:
      if spider2Raw is not None:
        tryToRemove(spider2Raw)
      pass
    
  def loadSpider2(self, fileName):
    '''
      Loads an spider2 file
      :param fileName: str. Path to spider2 file.

              #       AA      SS      ASA     Phi     Psi     Theta(i-1=>i+1) Tau(i-2=>i+1)   P(C)    P(E)    P(H)
        1       E       C       135.0    -89.3   131.8   116.5  -113.0  0.977   0.010   0.011
      
      :return  list of strings. ["row0_spider2","row1spider2"...]
    '''  
    dataList= []
    if not os.path.isfile(fileName):
      fileName+=".gz"
    with self.openForReadingFnameOrGz(fileName) as f:
      f.readline()
      for line in f:
        lineArray=line.split()
        if lineArray[0][0].isdigit():
          dataList.append( (lineArray[1],  tuple([lineArray[3]] + lineArray[8:] )))
        else:
          break
    return dataList  
       
