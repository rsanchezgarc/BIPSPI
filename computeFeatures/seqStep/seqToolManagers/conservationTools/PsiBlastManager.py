from __future__ import absolute_import
import sys, os
from subprocess import Popen, PIPE, check_output
import joblib

from ..seqToolManager import seqToolManager, FeatureComputerException
from .consDbManager import ConsDbManager
from .windowPSSM import WindowPSSM
from utils import myMakeDir, tryToRemove #utils is at the root of the package

INCLUDE_PSSM= True
INCLUDE_PSFM= True
class PsiBlastManager(seqToolManager):
  BAD_SCORE_CONSERVATION="-1048576"
  def __init__(self, seqsManager, outPath, winSize):
    seqToolManager.__init__(self, seqsManager, outPath, winSize)
    '''
      @param seqsManager: ..manageSeqs.seqsManager.SeqsManager 
      @param outPath: str. root path where pssm results will be saved (./procPssms, ./rawPssms ./windowedPSSMs/winSize%d)
      @param winSize: int. The size of sliding window
    '''
#    self.psiBlastBin Inherited from ../../../Config
#    self.psiblastDB Inherited from  ../../..Config
#    self.blastNThrs Inherited from  ../../../Config

    self.consDbManager= ConsDbManager()
    self.psiblastPathOutPath= myMakeDir(self.outPath,"psiblast")
    self.pssmsPathRaw= myMakeDir(self.outPath,"pssms/rawPssms")
    self.pssmsPathProc= myMakeDir(self.outPath,"pssms/procPssms")
    self.pssmsPathWindowed= myMakeDir(self.outPath,"pssms/windowedPSSMs/wSize"+str(winSize))

  def getFinalPath(self):
    '''
      returns path where final results (win pssms) are saved
      @return pssmsPathWindowed: str
    '''
    return self.pssmsPathWindowed
 
  def getFNames(self, prefixExtended):
    '''
    Returns a dict that contains the fnames that will be used by psiblast
    @param prefixExtended. prefix for output fnames. They are formed as follows: prefix+chainType+chainId+b/u
    @return Dict   {"psiblast":(psiblastOutName, ), "pssm":(pssmOutNameRaw, pssmOutNameProc), "pssmWindow":(pssmWindowedOutName,)}
      Processed pssm and pssmWindow are the ones that will be used for classification.
    '''
    psiblastOutName= os.path.join(self.psiblastPathOutPath, prefixExtended+".psiblast")

    pssmOutNameRaw= os.path.join(self.pssmsPathRaw, prefixExtended+".pssm")
    pssmOutNameProc= os.path.join(self.pssmsPathProc, prefixExtended+".pssm")

    pssmWindowedOutName= os.path.join(self.pssmsPathWindowed, prefixExtended+".wsize"+str(self.winSize)+".pssm")

    fNames= {"psiblast":(psiblastOutName, ), "pssm":(pssmOutNameRaw, pssmOutNameProc), "pssmWindow":(pssmWindowedOutName,)}
    return fNames

  def compute(self, prefixExtended, pdbCode=None):
    '''
      Computes psiblast for the sequence seqStr, that is contained at fastaInFname (located in 
                                                      computedFeatsRootDir/seqStep/extractedSeqs/seqsData/).
      This sequence is associated with prefixExtended as an unambiguous id

      @param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      @param pdbCode: str pdbId in order to query 3dCons. If None, psiblast will be launched directly
      @return psioutName, pssmNameRaw
              psioutName: str
              pssmNameRaw: str
    '''
    chainType, chainId= prefixExtended.split("_")[1:3]
    seqStr, fastaInFname= self.seqsManager.getSeq(chainType, chainId) # repeat as psiBlastManager can modify seqs
    fNames= self.getFNames( prefixExtended)
    psioutName= fNames["psiblast"][0]
    pssmNameRaw, pssmNameProc= fNames["pssm"]
    winPssmOutName= fNames["pssmWindow"][0]

    if self.checkAlreayComputed(prefixExtended):
      print("Blast already computed for %s"%prefixExtended)
      return psioutName, pssmNameRaw
    # run psi-blast
    print("getting psiblast results over %s"%prefixExtended)
    areSeqIdsMapped= self.getPsiBlast( fastaInFname, psioutName, pssmNameRaw, seqStr, prefixExtended, pdbCode)
    #Process psi-blast
    self.processPSSM( seqStr, prefixExtended, pssmNameRaw, pssmNameProc, areSeqIdsMapped)
    #Compute windows
    self.makeWindowedPSSM( pssmNameProc, winPssmOutName)
    return psioutName, pssmNameRaw

  def getPsiBlastFrom3DCons(self, psioutName, pssmNameRaw, seqStr, pdbCode, chainId):
    '''
      @param psioutName: str. Path to results file where aligments will be saved
      @param pssmNameRaw: str. Path to results file where pssms will be saved
      @param seqStr: str. Sequence of the chain
      @param pdbCode: str pdbId in order to query 3dCons. If None and no perfect match with seqStr, 
                          psiblast will be launched directly
      @param chainId: str chainId to recover

      @return 0 if seq was available as it is in 3dCons (seq serch, res mapping needed), 
              1 if pdbId_chainId available in 3dCons (aligment of seqId struct id included in pssms but not in aligments file),
              -1 if not available in 3dCons in any form
              2 if was not possible to compute it
    '''
    
    if not self.consDbManager.consDbIsAvailable():
      return -1
    if not pdbCode is None:
      returnVal= self.consDbManager.retrieve3DConsFromPDBChain(pdbCode, chainId, pssmNameRaw, psioutName)
      print("3dcons search status %d"%returnVal)
      if returnVal==0:
        return 1
    returnVal= self.consDbManager.retrieve3DConsFromSeq(seqStr, pssmNameRaw, psioutName)
    print("3dcons search status %d"%returnVal)
    if returnVal==0:
      return 0
    elif returnVal in [-1, -2,-3]:
      return -1
    else:
      return 2

  def getPsiBlast(self, fastaInFname, psioutName, pssmNameRaw, seqStr, prefixExtended, pdbCode):
    '''
      Retrieves psiblast output (pssm and aligments). First, it tries to download from
      3dconsDb. If not available, launches psiblast command with fastaInFname as input file,
      psioutName as the output file that will save the aligments and pssmNameRaw as the output
      file that will save the pssms.
      @param fastaInFname: str. Path to fasta file where sequence is saved
      @param psioutName: str. Path to results file where aligments will be saved
      @param pssmNameRaw: str. Path to results file where pssms will be saved
      @param seqStr: str. Sequence of the chain
      @param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      @param pdbCode: str pdbId in order to query 3dCons. If None and no perfect match with seqStr, 
                      psiblast will be launched directly

      @return seqIdsMapped: boolean. True if psiblast output is obtained from 3dConsDb and thus, 
                structIds are included in first col of pssm files instead of seqIds, false if
                first column are seqIds.
    '''
    
    prefix, chainType, chainId= prefixExtended.split("_")[:-1]
#    print(pdbCode)
#    raw_input("press enter to continue")
    returnVal= self.getPsiBlastFrom3DCons( psioutName, pssmNameRaw, seqStr, pdbCode, chainId)
    if returnVal==-1:
      if os.path.isfile(psioutName) and int(check_output('wc -l {}'.format(psioutName.replace("_*", "_\\*")), 
         shell=True).split()[0])> 11 and os.path.isfile(pssmNameRaw):
        print("psiblast raw files alredy computed")
      else:
        print("Not found in 3dconds db")
        self.launchBlast( fastaInFname, psioutName, pssmNameRaw, seqStr)
      return False
    elif returnVal==1:
      return True
    else:
      return False

  def launchBlast(self, fastaInFname, psioutName, pssmNameRaw, seqStr):
    '''
      Launches psiblast command with fastaInFname as input file, psioutName as the output file that will
      save the aligments and pssmNameRaw as the output file that will save the pssms.
      @param fastaInFname: str. Path to fasta file where sequence is saved
      @param psioutName: str. Path to results file where aligments will be saved
      @param pssmNameRaw: str. Path to results file where pssms will be saved
      @param seqStr: str. Sequence of the chain
    '''
    blastbin= self.psiBlastBin
    blastDB=  self.psiBlastDB
    numThr=   int(self.psiBlastNThrs)
    blastCMD = ("%(blastbin)s -query %(fastaInFname)s -num_iterations 3 -inclusion_ethresh 0.0001 -show_gis -db %(blastDB)s "+
                "-out %(psioutName)s -num_threads %(numThr)d -out_ascii_pssm %(pssmNameRaw)s -num_alignments 250")%locals()
    blastCMD= blastCMD.replace("_*", "_\\*")
    print(blastCMD)
    process= Popen( blastCMD, shell=True, stdout=PIPE, stderr=PIPE)
    processOut= process.communicate()
    #Check for failure
    if len(processOut[1])>0 and processOut[1].startswith(b"BLAST engine error"): #No hits found will be dealt at processResults
      print("Error computing blast. Caught stdin/stderr:\n",processOut[0],processOut[1])
      with open(psioutName,"w") as f:
        f.write(processOut[0]+"\n"+processOut[1])
      self.makeFakePSSM( pssmNameRaw, seqStr)
      
  def makeFakePSSM( self, pssmNameRaw, seqStr):
    header='''
Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts
            A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
'''
    oneFakeLine='''    %d %s     -1024  -1024  -1024  -1024  -1024  -1024  -1024  -1024  -1024   -1024  -1024  -1024  -1024  -1024  -1024   -1024   -1024  -1024  -1024   -1024    -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024  -1024  -1024   -1024   -1024  -1024  -1024.0 -1024.0'''
    with open(pssmNameRaw, "w") as f:
      f.write(header)
      for i, resName in enumerate(seqStr):
        f.write(oneFakeLine%(i, resName))
    return None

  def processPSSM(self, seq, prefixExtended, pssmNameRaw, pssmNameProc, areSeqIdsMapped):
    '''
      Reads psiblast pssms output file and writes another one with tabulated format, headers and
      some error checking.
      @param seq: str. Sequence of the chain
      @param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      @param pssmNameRaw: str. Path to psiblast aligments results
      @param pssmNameProc: str. Path where formated results will be saved.
      @param areSeqIdsMapped: boolean. True if psiblast output is obtained from 3dConsDb and thus, 
                structIds are included in first col of pssm files instead of seqIds, false if
                first column are seqIds.
    '''
    try:
      pssmData, pssmResIds, pssmSeq = self.loadPSSM(pssmNameRaw)
      if areSeqIdsMapped:
        seq= pssmSeq
        prefix, chainType, chainId, __= prefixExtended.split("_")
        self.seqsManager.addResiduesToSeqToStructMap( chainType, chainId, pssmSeq, pssmResIds)
      else:
        assert pssmSeq== seq
    except IOError:
      print("Pssm was not computed. Default value inserted instead")
      pssmSeq= seq
      pssmData= [ " ".join([PsiBlastManager.BAD_SCORE_CONSERVATION for i in range(42)]) for i in range(len(seq))]
    prefix, chainType, chainId, __= prefixExtended.split("_")
    try:
      outFile= open(pssmNameProc,"w")
      outFile.write("chainId seqIndex structResId resName "+ "pssm "*20+ "psfm "*20+ "score "*2+"\n")
      assert len(pssmData) == len(seq)
      for i, (pssmArrayJoined,letter) in enumerate(zip(pssmData,seq)):
        structIndex= self.seqsManager.seqToStructIndex(chainType, chainId, i, asString= True)

        if self.filterOutLabels and structIndex[-1].isalpha():
          continue
        outFile.write("%s %d %s %s "%(chainId, i, structIndex, letter)+ pssmArrayJoined +"\n")
      outFile.close()
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing %s"%pssmNameProc)
      tryToRemove(pssmNameProc)
      raise
    return pssmSeq

  def loadPSSM(self, fname):
    '''
      Loads a pssm file (psiblast format)
      @param fname: str. Path to pssm file.
      @return  list of strings. ["row0_PSSM_values","row1_PSSM_values"...]
    '''
    scores=[]
    indices=[]
    seq=""
    skip=2
    for line in open(fname):
      skip-=1
      if line.startswith("Last position-specific"):
        continue
      if skip==1:
        skip+=1
        continue
      if skip<0:
        if line.startswith("\n"):
          break
        lineArray=line.split()
#        print(lineArray)
        scores.append(" ".join(lineArray[2:]))
        indices.append(lineArray[0])
        seq+= lineArray[1]
    return scores, indices, seq
    
  def makeWindowedPSSM(self, pssmNameProc, winPssmOutName):
    '''
      Computes sliding windows for a given pssmFile. Windows will include aa code and pssm features
      @param pssmNameProc: str. Path to processed pssm file (my format)
      @param winPssmOutName: str. Path to windowed results.
    '''  
    try:
      WindowPSSM(self.winSize, True, INCLUDE_PSSM, INCLUDE_PSFM).compute(pssmNameProc, winPssmOutName)
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing %s"%winPssmOutName)
      tryToRemove(winPssmOutName)    
      raise     


def test():
  from ...manageSeqs.seqsManager  import SeqsManager
  rFname="/home/rsanchez/tmp/tmpRRI/wdir/in/1ACB_r_u.pdb"
  lFname="/home/rsanchez/tmp/tmpRRI/wdir/in/1ACB_l_u.pdb"
  wd="/home/rsanchez/tmp/tmpRRI/wdir/computedFeatures"
  seqsManager= SeqsManager( rFname, lFname, computedFeatsRootDir= wd)
  seqsManager.computeOneFile(rFname, "r")
  seqsManager.computeOneFile(lFname, "l")
  psiManger= PsiBlastManager(seqsManager, outPath="/home/rsanchez/tmp/tmpRRI/wdir/out", winSize=11)
  prefixExtended="1ACB_l_*_u"
  psiManger.compute( prefixExtended, pdbCode=None)
  return None
  
def test():
  from ...manageSeqs.seqsManager  import SeqsManager
  rFname="/home/rsanchez/tmp/tmpRRI/wdir/in/1ACB_r_b.pdb"
  lFname="/home/rsanchez/tmp/tmpRRI/wdir/in/1ACB_l_b.pdb"
  wd="/home/rsanchez/tmp/tmpRRI/wdir/computedFeatures"
  seqsManager= SeqsManager( rFname, lFname, computedFeatsRootDir= wd)
  seqsManager.computeOneFile(rFname, "r")
  seqsManager.computeOneFile(lFname, "l")
  psiManger= PsiBlastManager(seqsManager, outPath="/home/rsanchez/tmp/tmpRRI/wdir/out", winSize=11)
  prefixExtended="1ACB_l_I_u"
#  psiManger.compute( prefixExtended, pdbCode=None)
  psiManger.compute( prefixExtended, pdbCode="1ACB")
  return None
  
if __name__=="__main__":
  test()
