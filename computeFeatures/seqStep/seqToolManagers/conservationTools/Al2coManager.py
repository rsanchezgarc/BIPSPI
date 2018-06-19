from __future__ import absolute_import, print_function
import sys, os
from subprocess import Popen, PIPE, check_output
from Bio.PDB.Polypeptide import aa1 as AA_STANDARD
from ..seqToolManager import seqToolManager, FeatureComputerException
from utils import myMakeDir, tryToRemove #utils is at the root of the package

class Al2coManager(seqToolManager):
  '''
    Computes al2co and processes their outputs. Extends class seqToolManager
  '''
  BAD_SCORE_CONSERVATION = "-1048576"  #Something went wrong tag
  def __init__(self, seqsManager, outPath, winSize):
    '''
      @param seqsManager: ..manageSeqs.seqsManager.SeqsManager 
      @param outPath: str. path where al2co scores will be saved
      @param winSize: int. The size of sliding window  NOT USED
    '''
    seqToolManager.__init__(self, seqsManager, outPath, winSize)
    
    self.al2coPerlScript= os.path.join(self.al2coRootDir, "get_al2co_scores.pl")
    self.al2coOutPath= myMakeDir(self.outPath,"al2co")
    
  def getFinalPath(self):
    '''
      returns path where results are saved
      @return al2coOutPath: str
    '''
    return self.al2coOutPath
    
  def getFNames(self, prefixExtended):
    '''
      Returns a dict that contains the fnames that will be used al2co
      @param prefixExtended. prefix for output fnames. They are formed as follows: prefix+chainType+chainId+b/u
      @return Dict   {"al2co":(al2coRaw, al2coProc) } Final scores will be saved at al2coProc
    '''
    al2coRaw= os.path.join(self.al2coOutPath, prefixExtended+".fasta.csv")
    al2coProc= os.path.join(self.al2coOutPath, prefixExtended+".al2co")
    return {"al2co":(al2coRaw, al2coProc) }
    
  def checkAlreayComputed(self, prefixExtended):
    '''
      Overrides checkAlreayComputed of seqToolManager
      checks if output files have been computed
      @param prefixExtended. prefix for output fnames. They are formed as follows: prefix+chainType+chainId+b/u
      @return boolean. True if prefixExtended sequences have been already computed, False otherwise
    '''
    al2coRawName, al2coProcName= self.getFNames(prefixExtended)["al2co"]
    if not os.path.isfile(al2coProcName):
      return False
    else:
      nlines=int( check_output('wc -l {}'.format(al2coProcName), shell=True).split()[0])
      if nlines<18: #Too small file to be correct
        return False
    return True
    
  def compute(self, prefixExtended, psiblastOutName, pssmOutNameRaw):
    '''
      Computes al2co for the sequence seqStr, that is contained at fastaInFname. This sequence is
      associated with prefixExtended as an unambiguous id
      @param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      @param psiblastOutName: str. Path to psiblast aligments results
      @param pssmOutNameRaw: str. Path to psiblast pssms results
    '''
    chainType, chainId= prefixExtended.split("_")[1:3]
    seqStr, fastaFname= self.seqsManager.getSeq(chainType, chainId) # repeat as psiBlastManager can modify seqs     
    if self.checkAlreayComputed(prefixExtended):
      print("Al2co already computed for %s"%prefixExtended)
      return 0
    print("lauching al2co over %s"%prefixExtended)      
    al2coRawName, al2coProcName= self.getFNames(prefixExtended)["al2co"]

    process= Popen(["perl", self.al2coPerlScript, fastaFname, self.al2coRootDir,
                    psiblastOutName, pssmOutNameRaw,
                    self.al2coOutPath], stdout=PIPE, stderr=PIPE)
    processOut= process.communicate()
    if len(processOut[1])>0:
      print("Error computing al2co. Caught stdin/stderr:\n",processOut[0],processOut[1])
#      raise FeatureComputerException("Al2co was not able to compute scores")
    self.processAl2co(seqStr, prefixExtended, al2coRawName, al2coProcName)
    
  def processAl2co(self, seq, prefixExtended, al2coRaw, al2coProc):
    '''
      Reads al2co output file and writes another one with tabulated format, headers and
      some error checking.
      @param: seq: str. Sequence of the chain
      @param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      @param al2coRaw: str. Path to al2co results
      @param al2coProc: str. Path where formated results will be saved.
    '''
    try:
      conserData = self.loadAl2co(al2coRaw)
    except IOError:
      conserData= [ (letter, Al2coManager.BAD_SCORE_CONSERVATION) for letter in seq]
    prefix, chainType, chainId, __= prefixExtended.split("_")
    try:
      outFile= open(al2coProc,"w")
      outFile.write("chainId seqIndex structResId resName score\n")

#      print(conserData)
      alcoIx=0
      seqIx=0
      seqLen= len(seq)
      alcoLen= len(conserData)
      while seqIx<seqLen and alcoIx<alcoLen:
        letter= seq[seqIx]
        letterAl2co, consVal= conserData[alcoIx]
        if letterAl2co== letter:
           structIndex= self.seqsManager.seqToStructIndex(chainType, chainId, seqIx, asString= True)
           if self.filterOutLabels and structIndex[-1].isalpha():
             continue
           outFile.write("%s %d %s %s %s\n"%(chainId, seqIx, structIndex, letter, consVal))
           alcoIx+=1
           seqIx+=1
        elif not letter in AA_STANDARD and letterAl2co=="-":
           alcoIx+=1
           seqIx+=1
        elif letterAl2co=="-":
          alcoIx+=1
        else:
          print(conserData)
          print(seq)
          print(alcoIx, seqIx)
          raise ValueError("Al2co mismatch %s %s "%(letterAl2co, letter))
      outFile.close()
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing %s"%al2coProc)
      tryToRemove(al2coProc)
      raise
    finally:
      tryToRemove(al2coRaw)
      pass
  def loadAl2co(self, filename):
    '''
      Loads an al2co file
      @param fname: str. Path to al2co file.
      @return  list of strings. ["row0_Al2co","row1Al2co"...]
    '''  
    conserv= []
#    print(filename)
    for line in open(filename):
      lineArray=line.split()
#      print(lineArray)
      if lineArray[0][0].isdigit():
        conserv.append(lineArray[1:3])
      else:
        break
    return conserv     
