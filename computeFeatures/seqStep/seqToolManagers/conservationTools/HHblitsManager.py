from __future__ import absolute_import
import os
from subprocess import Popen, PIPE, check_output
import socket
from ...seqToolManager import SeqToolManager
from ....featuresComputer import FeatureComputerException
from utils import myMakeDir #utils is at the root of the package

class HHBlitsManager(SeqToolManager):
  VAR_LIST= ["hhblits%d"%i for i in range(31)]
  BAD_SCORE_CONSERVATION="-1024"
  def __init__(self,  computedFeatsRootDir, winSize, 
                     hhBlitsCMD_template="%(hhBlitsBinPath)s/hhblits -i %(fastaInFname)s -n 4 -d %(hhblitsDB)s "+
                     "-oa3m %(aligsName)s -cpu %(hhBlitsNThrs)d -ohhm %(profileNameRaw)s -o /dev/null", statusManager=None):
    
    '''
      :param computedFeatsRootDir: str. root path where results will be saved
      :param winSize: int>=1 or None. The size of the windows for sliding window if desired
      :param hhBlitsCMD_template: str. The shell command that will be used when computing hhblits in order to get the
                                       best correlated mutation scores
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''  
    SeqToolManager.__init__(self, computedFeatsRootDir, winSize)
    
#    self.hhBlitsBinPath Inherited from ../../Config
#    self.hhBlitsDB Inherited from ../../Config
#    self.hhBlitsNThrs Inherited from ../../Config
    
    self.hhBlitsOut=  myMakeDir(computedFeatsRootDir,"hhBlits")
    self.hhBlitsRaw=  myMakeDir(self.hhBlitsOut,"rawHhBlits")
    self.hhBlitsProc= myMakeDir(self.hhBlitsOut,"procHhBlits")
    if winSize:
      self.hhBlitsPathWindowed= myMakeDir(self.hhBlitsOut,"windowedHhBlits/wSize"+str(winSize))
    else:
      self.hhBlitsPathWindowed= None
    self.hhBlitsCMD_template= hhBlitsCMD_template
    
  def getFinalPath(self):
    '''
      returns path where final results (win hhBlits) are saved
      :return self.hhBlitsPathWindowed: str
    '''  
    return self.hhBlitsProc
        
  def getFNames(self, prefixExtended):
    '''
    Returns a list that contains the fnames that will be produced by hhblits
    :param prefixExtended. prefix for output fnames. They are formed as follows: prefix+chainType+chainId+
    :return list of fnames: [ fname1, fnam2, ...]
    '''
    hhBlitsAligName= os.path.join( self.hhBlitsRaw,  prefixExtended+".a3m")
    procHhblits= os.path.join(self.hhBlitsProc, prefixExtended+".ohhm.gz")
    fNames= [ hhBlitsAligName, procHhblits]
    if self.winSize:
      fNames+= [os.path.join(self.hhBlitsPathWindowed, prefixExtended+".wsize"+str(self.winSize)+".ohhm.gz")]
    return fNames

  def computeFromSeqStructMapper(self, seqStructMap, prefixExtended):    
    '''
      Computes hhblits for the sequence associated with prefixExtended as an unambiguous id and included in
                                  self.seqStructMap
      :param seqStructMap: computeFeatures.seqStep.seqToolManagers.seqExtraction.SeqStructMapper
      :param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names. Must be included
                                  in self.seqStructMap
      :return (profileNameProc, winProfileOutName)
                  profileNameProc: str
                  winProfileOutName: str    
    '''

    __, chainType, chainId = self.splitExtendedPrefix( prefixExtended)[:3]
    seqStr, fastaInFname= seqStructMap.getSeq(chainType, chainId)
    seqStructMap.setCurrentSeq(seqStr, chainType, chainId)
    fNames= self.getFNames( prefixExtended)
    aligsName=  fNames[0]
    profileNameProc= fNames[1]
    profileNameRaw= os.path.join( self.hhBlitsRaw,  prefixExtended+".ohhm")
    if self.checkAlreayComputed(prefixExtended):
      print("hhblits already computed for %s"%prefixExtended)
      return aligsName, profileNameRaw, profileNameProc
    # run hhblits
    print("launching hhblits over %s"%prefixExtended)
    self.launchHhblits( fastaInFname, aligsName, profileNameRaw)
    #Process psi-blast
    dataList= self.processHhblits( seqStr, seqStructMap, prefixExtended, profileNameRaw, profileNameProc)
    if self.winSize:
      #Compute windows
      winProfileOutName= fNames[-1]
      self.makeWindowed( dataList, ["hhBlitsProfile"], [ HHBlitsManager.BAD_SCORE_CONSERVATION], [None], winProfileOutName)
    else:
      winProfileOutName= None
    return  aligsName, profileNameProc, winProfileOutName
    
  def compressRawData(self, prefixExtended):
    for fname in os.listdir(self.hhBlitsRaw):
      if fname.startswith(prefixExtended) and not fname.endswith(".gz"):
        if os.path.islink(os.path.join(self.hhBlitsRaw, fname)):
          os.remove(os.path.join(self.hhBlitsRaw, fname))
        else:
          cmd= ["gzip", "-f", os.path.join(self.hhBlitsRaw, fname)]
          check_output(cmd)
        
  def launchHhblits(self, fastaInFname, aligsName, profileNameRaw):
    '''
      Launches hhblits command with fastaInFname as input file, aligsName as the output file that will
      contain the aligments and profileNameRaw as the output file that will contain the profile.
      :param fastaInFname: str. Path to fasta file where sequence is saved
      :param aligsName: str. Path to results file where aligments will be saved
      :param profileNameRaw: str. Path to results file where profile will be saved
    '''
    if os.path.isfile(profileNameRaw) and int(check_output('wc -l {}'.format(profileNameRaw), 
                                              shell=True).split()[0])> self.minNumResiduesPartner :
      print("hhblits raw files alredy computed")
      return
        
    hhBlitsBinPath=  self.hhBlitsBinPath
    hhblitsDB=   self.hhBlitsDB
    hhBlitsNThrs= self.hhBlitsNThrs if socket.gethostname()!="servet" else 1
    hhblitsCMD = self.hhBlitsCMD_template%locals()
    hhblitsCMD.replace("_*", "_\\*")
    print(hhblitsCMD)
    process= Popen( hhblitsCMD, shell=True, stdout=PIPE, stderr=PIPE)
    processOut= process.communicate()
    #Check for failure
    if len(processOut[1])>0 and "Error" in processOut[1]: #No hits found will be dealt at processResults
      print("Error computing hhblits. Caught stdin/stderr:\n",processOut[0],processOut[1])
      raise FeatureComputerException("hhblits was not able to compute profile")
      
  def processHhblits(self, seq, seqStructMap, prefixExtended, profileNameRaw, profileNameProc):
    '''
      Reads hhblits profile output file and writes another one with tabulated format, headers and
      some error checking.
      :param: seq: str. Sequence of the chain
      :param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      :param profileNameRaw: str.  Path to profiles results
      :param profileNameProc: str. Path where formated results will be saved.
    '''
    try:
      hhBlitsData = self.loadHhblits(profileNameRaw)
    except IOError:
      hhBlitsData= [ [HHBlitsManager.BAD_SCORE_CONSERVATION for i in range(31)] for i in range(len(seq))]
    __, chainType, chainId = self.splitExtendedPrefix(prefixExtended)[:3]
    dataList=[]
    listOfRowsToPrint=[]
#    headerStr= "chainId resId resName "+ " ".join(["hhblits%d"%i for i in range(31)])+"\n"

    for i, (hhBlitsArray,letter) in enumerate(zip(hhBlitsData,seq)):
      structIndex= seqStructMap.seqToStructIndex(chainType, chainId, i, asString= True)
      if structIndex:
        if self.filterOutLabels and structIndex[-1].isalpha():
          continue
      else:
        structIndex=str(i)+"?"
      listOfRowsToPrint.append( " ".join([ chainId, structIndex, letter]+hhBlitsArray) )
      dataList.append( ((chainId, structIndex, letter), (hhBlitsArray,) ) )
    self.writeResultsFromDataDictSingleChain( {chainId: listOfRowsToPrint }, outName= profileNameProc)
    return dataList
    
  def loadHhblits(self, fname):
    '''
      Loads a hhblits profile file
      :param fname: str. Path to hhblits profile file.
      :return  list of strings. ["row0_hhblits_values","row1_hhblits_values"...]
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
          scores[-1]= [ elem if elem!="*" else "-1" for elem in scores[-1] ]
    return scores
    

