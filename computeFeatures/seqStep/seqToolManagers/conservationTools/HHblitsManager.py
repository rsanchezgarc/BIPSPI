from __future__ import absolute_import
import sys, os
from subprocess import Popen, PIPE, check_output
import socket

from ..seqToolManager import seqToolManager, FeatureComputerException
from .windowHHblits import WindowHHblits
from utils import myMakeDir, tryToRemove #utils is at the root of the package

class HHBlitsManager(seqToolManager):
  BAD_SCORE_CONSERVATION="-1048576"
  def __init__(self, seqsManager, outPath, winSize, 
                     hhBlitsCMD_template="%(hhBlitsBinPath)s/hhblits -i %(fastaInFname)s -n 4 -d %(hhblitsDB)s "+
                     "-oa3m %(aligsName)s -cpu %(psiBlastNThrs)d -ohhm %(profileNameRaw)s -o /dev/null"):
    seqToolManager.__init__(self, seqsManager, outPath, winSize)
    '''
      @param seqsManager: ..manageSeqs.seqsManager.SeqsManager 
      @param outPath: str. path where hhblits results will be saved
      @param winSize: int. The size of sliding window
    '''
    
#    self.hhBlits should be inherited from
#    self.hhBlitsBinPath Inherited from ../../Config
#    self.hhBlitsDB Inherited from ../../Config
#    self.psiBlastNThrs Inherited from ../../Config
    
    self.hhBlitsOut=  myMakeDir(self.outPath,"hhBlits")
    self.hhBlitsRaw=  myMakeDir(self.hhBlitsOut,"rawHhBlits")
    self.hhBlitsProc= myMakeDir(self.hhBlitsOut,"procHhBlits")
    self.hhBlitsPathWindowed= myMakeDir(self.hhBlitsOut,"windowedHhBlits/wSize"+str(winSize))
    self.hhBlitsCMD_template= hhBlitsCMD_template
  def getFinalPath(self):
    '''
      returns path where final results (win hhBlits) are saved
      @return self.hhBlitsPathWindowed: str
    '''  
    return self.hhBlitsPathWindowed
        
  def getFNames(self, prefixExtended):
    '''
    Returns a dict that contains the fnames that will be used by hhblits
    @param prefixExtended. prefix for output fnames. They are formed as follows: prefix+chainType+chainId+b/u
    @return Dict   {"psiblast":(psiblastOutName, ), "pssm":(pssmOutNameRaw, pssmOutNameProc), "pssmWindow":(pssmWindowedOutName,)}
      Processed pssm and pssmWindow are the ones that will be used for classification.
    '''
    hhBlitsAligName= os.path.join( self.hhBlitsRaw,  prefixExtended+".a3m")
    rawHhblits= os.path.join( self.hhBlitsRaw,  prefixExtended+".ohhm")
    procHhblits= os.path.join(self.hhBlitsProc, prefixExtended+".ohhm") 
    hhblitsWindowedOutName= os.path.join(self.hhBlitsPathWindowed, prefixExtended+".wsize"+str(self.winSize)+".ohhm")

    fNames= { "hhBlitsAligName":(hhBlitsAligName,),"hhBlitsProfiles":(rawHhblits, procHhblits), 
              "hhBlitsProfilesWindow":(hhblitsWindowedOutName,)}
    return fNames

  def compute(self, prefixExtended):
    '''
      Computes hhblits for the sequence associated with prefixExtended as an unambiguous id and included in
                                  self.seqsManager
      @param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names. Must be included
                                  in self.seqsManager
      @return (profileNameProc, winProfileOutName)
                  profileNameProc: str
                  winProfileOutName: str    
    '''
    prefix, chainType, chainId, __= prefixExtended.split("_")
    seqStr, fastaInFname= self.seqsManager.getSeq(chainType, chainId) 
    fNames= self.getFNames( prefixExtended)
    aligsName=  fNames["hhBlitsAligName"][0]
    profileNameRaw, profileNameProc= fNames["hhBlitsProfiles"]
    winProfileOutName= fNames["hhBlitsProfilesWindow"][0]
    if self.checkAlreayComputed(prefixExtended):
      print("hhblits already computed for %s"%prefixExtended)
      return aligsName, profileNameRaw, profileNameProc
    # run hhblits
    print("lauching hhblits over %s"%prefixExtended)
    self.launchHhblits( fastaInFname, aligsName, profileNameRaw)
    #Process psi-blast
    self.processHhblits( seqStr, prefixExtended, profileNameRaw, profileNameProc)
    #Compute windows    
    self.makeWindowedPSSMHhblits( profileNameProc, winProfileOutName)
    return  aligsName, profileNameProc, winProfileOutName
    
  def launchHhblits(self, fastaInFname, aligsName, profileNameRaw):
    '''
      Launches hhblits command with fastaInFname as input file, aligsName as the output file that will
      contain the aligments and profileNameRaw as the output file that will contain the profile.
      @param fastaInFname: str. Path to fasta file where sequence is saved
      @param aligsName: str. Path to results file where aligments will be saved
      @param profileNameRaw: str. Path to results file where profile will be saved
    '''
    if os.path.isfile(profileNameRaw) and int(check_output('wc -l {}'.format(profileNameRaw), shell=True).split()[0])> 11:
      print("hhblits raw files alredy computed")
      return
        
    hhBlitsBinPath=  self.hhBlitsBinPath
    hhblitsDB=   self.hhBlitsDB

    psiBlastNThrs= self.psiBlastNThrs if socket.gethostname()!="servet" else 1
    hhblitsCMD = self.hhBlitsCMD_template%locals()
    hhblitsCMD.replace("_*", "_\\*")
    print(hhblitsCMD)
    process= Popen( hhblitsCMD, shell=True, stdout=PIPE, stderr=PIPE)
    processOut= process.communicate()
    #Check for failure
    if len(processOut[1])>0 and "Error" in processOut[1]: #No hits found will be dealt at processResults
      print("Error computing hhblits. Caught stdin/stderr:\n",processOut[0],processOut[1])
      raise FeatureComputerException("hhblits was not able to compute profile")
      
  def processHhblits(self, seq, prefixExtended, profileNameRaw, profileNameProc):
    '''
      Reads hhblits profile output file and writes another one with tabulated format, headers and
      some error checking.
      @param: seq: str. Sequence of the chain
      @param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      @param profileNameRaw: str.  Path to profiles results
      @param profileNameProc: str. Path where formated results will be saved.
    '''
    try:
      hhBlitsData = self.loadHhblits(profileNameRaw)
    except IOError:
      hhBlitsData= [ " ".join([HHBlitsManager.BAD_SCORE_CONSERVATION for i in range(31)]) for i in range(len(seq))]
    prefix, chainType, chainId, __= prefixExtended.split("_")
    try:    
      outFile= open(profileNameProc,"w")
      outFile.write("chainId seqIndex structResId resName "+ "hhblits "*31+"\n")
      for i, (hhBlitsArrayJoined,letter) in enumerate(zip(hhBlitsData,seq)):
        structIndex= self.seqsManager.seqToStructIndex(chainType, chainId, i, asString= True)
        if self.filterOutLabels and structIndex[-1].isalpha():
          continue
        outFile.write("%s %d %s %s "%(chainId, i, structIndex, letter)+ hhBlitsArrayJoined +"\n")
      outFile.close()
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing %s"%profileNameProc)
      tryToRemove(profileNameProc)    
      raise

  def loadHhblits(self, fname):
    '''
      Loads a hhblits profile file
      @param fname: str. Path to hhblits profile file.
      @return  list of strings. ["row0_hhblits_values","row1_hhblits_values"...]
    '''
    scores=[]
    begin=False
    count=0    
    with open(fname) as f:
      for line in f:
        if line.startswith("#"):
          begin=True
          continue
        if begin==True:
          count+=1
        if count==4:
          break          
      for i,line in enumerate(f):
        if line.startswith("//"):
          break
        lineArray= line.split()
        nElems= len(lineArray)
        if i%3==0:
          lineArray= lineArray[2:]
          if nElems!=23:
            raise ValueError("Bad format in hhblits file %s"%fname)
          scores.append(lineArray)
        elif i%3==1:
          scores[-1]+= lineArray
          scores[-1]= " ".join([ elem if elem!="*" else "-1" for elem in scores[-1] ])
    return scores
    
  def makeWindowedPSSMHhblits(self, profileNameProc, winProfileOutName):
    '''
      Computes sliding windows for a given profileNameProc.
      @param profileNameProc: str. Path to processed hhblits profile file
      @param winProfileOutName: str. Path to windowed results.
    '''  
    try:
      WindowHHblits(self.winSize).compute(profileNameProc, winProfileOutName)
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing %s"%winProfileOutName)
      tryToRemove(winProfileOutName)    
      raise     

def test():
  fname="/home/rsanchez/Tesis/rriPredMethod/dependencies/bioinformaticTools/hh-Tools/seqExample.ohhm"
  from computeFeatures.seqStep.manageSeqs.seqsManager import SeqsManager
  seqManag= SeqsManager("rFname", "lFname", computedFeatsRootDir= ".")
  hhblitsObj= HHBlitsManager( seqsManager= seqManag, outPath=".", winSize=11)
  hhblitsObj.loadHhblits(fname)

if __name__=="__main__":
  test()
