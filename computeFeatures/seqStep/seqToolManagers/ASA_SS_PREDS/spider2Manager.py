from __future__ import absolute_import, print_function
import sys, os
from subprocess import Popen, PIPE, check_output
from ..seqToolManager import seqToolManager, FeatureComputerException
from utils import myMakeDir, tryToRemove #utils is at the root of the package

class Spider2Manager(seqToolManager):
  '''
    Computes asa and secondary structure predictions and processes their outputs. 
    Extends class seqToolManager
  '''
  BAD_SCORE_CONSERVATION = "-1048576"  #Something went wrong tag
  def __init__(self, seqsManager, outPath, winSize):
    '''
      @param seqsManager: ..manageSeqs.seqsManager.SeqsManager 
      @param outPath: str. path where predictions will be saved
      @param winSize: int. The size of sliding window  NOT USED
    '''
    seqToolManager.__init__(self, seqsManager, outPath, winSize)
    
    # self.spider2PyScript inherited in Config.py
    self.spider2OutPath= myMakeDir(self.outPath,"SPIDER2")
    
  def getFinalPath(self):
    '''
      returns path where results are saved
      @return spider2OutPath: str
    '''
    return self.spider2OutPath
    
  def getFNames(self, prefixExtended):
    '''
      Returns a dict that contains the fnames that will be used spider2
      @param prefixExtended. prefix for output fnames. They are formed as follows: prefix+chainType+chainId
      @return Dict   {"spider2":(spider2Raw, spider2Proc) } Final scores will be saved at spider2Proc
    '''
    spider2Raw= os.path.join(self.spider2OutPath, prefixExtended+".spd3")
    spider2Proc= os.path.join(self.spider2OutPath, prefixExtended+".spider2")
    return {"spider2":(spider2Raw, spider2Proc) }
    
  def checkAlreayComputed(self, prefixExtended):
    '''
      Overrides checkAlreayComputed of seqToolManager
      checks if output files have been computed
      @param prefixExtended. prefix for output fnames. They are formed as follows: prefix+chainType+chainId+b/u
      @return boolean. True if prefixExtended sequences have been already computed, False otherwise
    '''
    spider2RawName, spider2ProcName= self.getFNames(prefixExtended)["spider2"]
    if not os.path.isfile(spider2ProcName):
      return False
    else:
      nlines=int( check_output('wc -l {}'.format(spider2ProcName), shell=True).split()[0])
      if nlines<18: #Too small file to be correct
        return False
    return True
    
  def compute(self,  prefixExtended, pssmOutNameRaw):
    '''
      Computes spider2 for the sequence seqStr, whose pssm is contained at pssmOutNameRaw. This sequence is
      associated with prefixExtended as an unambiguous id
      @param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      @param pssmOutNameRaw: str. Path to psiblast pssms results    
    '''
    chainType, chainId= prefixExtended.split("_")[1:3]
    seqStr, fastaFname= self.seqsManager.getSeq(chainType, chainId) # repeat as psiBlastManager can modify seqs 
    if self.checkAlreayComputed(prefixExtended):
      print("spider2 already computed for %s"%prefixExtended)
      return 0
    print("lauching spider2 over %s"%prefixExtended)      
    spider2RawName, spider2ProcName= self.getFNames(prefixExtended)["spider2"]
    
    curWd= os.getcwd()
    os.chdir(self.spider2OutPath)
    process= Popen(["python", self.spider2PyScript, pssmOutNameRaw], stdout=PIPE, stderr=PIPE)
    processOut= process.communicate()
    os.chdir(curWd)
    if len(processOut[1])>0:
      print("Error computing spider2. Caught stdin/stderr:\n",processOut[0],processOut[1])
    self.processspider2(seqStr, prefixExtended, spider2RawName, spider2ProcName)
    
  def processspider2(self, seq, prefixExtended, spider2Raw, spider2Proc):
    '''
      Reads spider2 output file and writes another one with tabulated format, headers and
      some error checking.
      @param: seq: str. Sequence of the chain
      @param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      @param spider2Raw: str. Path to spider2 results
      @param spider2Proc: str. Path where formated results will be saved.
        head spider2Proc
        chainId seqIndex structResId resName score_asa score_Pc score_Pe score_Ph
        A 0 6 E 146.4 0.982 0.011 0.006
        A 1 7 P 100.2 0.977 0.012 0.012
    '''
    try:
      predictionsData = self.loadspider2(spider2Raw)
    except IOError:
      predictionsData= [ (letter, (tuple([Spider2Manager.BAD_SCORE_CONSERVATION]*4)) ) for letter in seq]
    prefix, chainType, chainId, __= prefixExtended.split("_")
    try:
      outFile= open(spider2Proc,"w")
      outFile.write("chainId seqIndex structResId resName score_asa score_Pc score_Pe score_Ph\n")

      predsIx=0
      seqIx=0
      seqLen= len(seq)
      alcoLen= len(predictionsData)
      while seqIx<seqLen and predsIx<alcoLen:
        letter= seq[seqIx]
        letterspider2, consValTuple= predictionsData[predsIx]
        if letterspider2== letter:
           structIndex= self.seqsManager.seqToStructIndex(chainType, chainId, seqIx, asString= True)
           if self.filterOutLabels and structIndex[-1].isalpha():
             continue
           outFile.write("%s %d %s %s %s %s %s %s\n"%((chainId, seqIx, structIndex, letter)+ consValTuple))
           predsIx+=1
           seqIx+=1
        elif letter=="X" and letterspider2=="-":
           predsIx+=1
           seqIx+=1
        elif letterspider2=="-":
          predsIx+=1
        else:
          print(predictionsData)
          print(seq)
          print(predsIx, seqIx)
          raise ValueError("spider2 mismatch %s %s "%(letterspider2, letter))
#      for i, (letter, (consVal,letterspider2)) in enumerate(zip(seq, predictionsData)):
#        if letter!="X" and  letterspider2!= letter: continue
#        structIndex= self.seqsManager.seqToStructIndex(chainType, chainId, i, asString= True)
#        if self.filterOutLabels and structIndex[-1].isalpha():
#          continue
#        outFile.write("%s %d %s %s %s\n"%(chainId, i, structIndex, letter, consVal))

      outFile.close()
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing %s"%spider2Proc)
      tryToRemove(spider2Proc)
      raise
    finally:
      tryToRemove(spider2Raw)
      pass
      
  def loadspider2(self, filename):
    '''
      Loads an spider2 file
      @param filename: str. Path to spider2 file.
        head    filename   
              #       AA      SS      ASA     Phi     Psi     Theta(i-1=>i+1) Tau(i-2=>i+1)   P(C)    P(E)    P(H)
        1       E       C       135.0    -89.3   131.8   116.5  -113.0  0.977   0.010   0.011
      
      @return  list of strings. ["row0_spider2","row1spider2"...]
    '''  
    dataList= []
    with open(filename) as f:
      f.readline()
      for line in f:
        lineArray=line.split()
        if lineArray[0][0].isdigit():
          dataList.append( (lineArray[1],  tuple([lineArray[3]] + lineArray[8:] )))
        else:
          break
    return dataList     
