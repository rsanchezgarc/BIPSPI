from __future__ import absolute_import
import  os
from subprocess import Popen, PIPE, check_output
import gzip
import numpy as np
from ...seqToolManager import SeqToolManager
from .consDbManager import ConsDbManager
from utils import myMakeDir, getItemsFromList 

INCLUDE_PSSM= True
INCLUDE_PSFM= True
NORMALIZE_IPP_AND_RWGRMP=True

class PsiBlastManager(SeqToolManager):
  VAR_LIST= ["pssm%d"%i for i in range(20)]+["psfm%d"%i for i in range(20)]+["IPP", "RWGRMP"]
  if NORMALIZE_IPP_AND_RWGRMP:
    VAR_LIST+= ["IPP_norm", "RWGRMP_norm"]  
  BAD_SCORE_CONSERVATION="-1024"
  CONS_DB_MANAGER= ConsDbManager()
  def __init__(self, computedFeatsRootDir, winSize):
    '''
      :param computedFeatsRootDir: str. root path where results will be saved
      :param winSize: int>=1 or None. The size of the windows for sliding window if desired
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''  
    SeqToolManager.__init__(self, computedFeatsRootDir, winSize)

#    self.psiBlastBin Inherited from ../../../Config
#    self.psiblastDB Inherited from  ../../..Config
#    self.blastNThrs Inherited from  ../../../Config

    self.consDbManager= PsiBlastManager.CONS_DB_MANAGER
    self.psiblastPathOutPath= myMakeDir(self.computedFeatsRootDir,"psiblast")
    self.pssmsPathRaw= myMakeDir(self.computedFeatsRootDir,"pssms/rawPssms")
    self.pssmsPathProc= myMakeDir(self.computedFeatsRootDir,"pssms/procPssms")
    if winSize:
      self.pssmsPathWindowed= myMakeDir(self.computedFeatsRootDir,"pssms/windowedPSSMs/wSize"+str(winSize))
    else:
      self.pssmsPathWindowed= None

  def getFinalPath(self):
    '''
      returns path where final results (win pssms) are saved
      :return pssmsPathWindowed: str
    '''
    return self.pssmsPathWindowed
 
  def getFNames(self, prefixExtended):
    '''
    Returns a list that contains the fnames that will be produced by psiblast
    :param prefixExtended. prefix for output fnames. They are formed as follows: prefix+chainType+chainId+
    :return list of fnames: [ fname1, fnam2, ...]
    '''
    psiblastOutName= os.path.join(self.psiblastPathOutPath, prefixExtended+".psiblast")
    pssmOutNameProc= os.path.join(self.pssmsPathProc, prefixExtended+".pssm.gz")
    fNames= [pssmOutNameProc, psiblastOutName   ]
    if self.winSize:
      pssmWindowedOutName= os.path.join(self.pssmsPathWindowed, prefixExtended+".wsize"+str(self.winSize)+".pssm.gz")
      fNames+= [ pssmWindowedOutName ]
    return fNames

  def computeFromSeqStructMapper(self, seqStructMap, prefixExtended, pdbCode=None):
    '''
      Computes psiblast for the sequence seqStr, that is contained at fastaInFname (located in 
      computedFeatsRootDir/seqStep/extractedSeqs/seqsData/).
      This sequence is associated with prefixExtended as an unambiguous id. 
      WARNING: If the sequence is contained in 3dcons-db, the seqStructMap object my be modified if some discrepancies exists.
      :param seqStructMap: computeFeatures.seqStep.seqToolManagers.seqExtraction.SeqStructMapper
      :param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      :param pdbCode: str pdbId in order to query 3dCons. If None, psiblast will be launched directly
      :return psioutName, pssmNameRaw
              psioutName: str
              pssmNameRaw: str
    '''
    __, chainType, chainId = self.splitExtendedPrefix(prefixExtended)[:3]
    seqStr, fastaInFname= seqStructMap.getSeq(chainType, chainId) 
    fNames= self.getFNames( prefixExtended)
    pssmNameProc= fNames[0]
    psioutName= fNames[1]
    pssmNameRaw= os.path.join(self.pssmsPathRaw, prefixExtended+".pssm")
    if self.checkAlreayComputed(prefixExtended):
      print("Blast already computed for %s"%prefixExtended)
      return psioutName, pssmNameRaw
    # run psi-blast
    print("getting psiblast results over %s"%prefixExtended)
    areSeqIdsMapped= self.getPsiBlast( fastaInFname, psioutName, pssmNameRaw, seqStr, prefixExtended, pdbCode)
    #Process psi-blast
    seqStr, fastaInFname= seqStructMap.getSeq(chainType, chainId) 
    seqStructMap.setCurrentSeq(seqStr, chainType, chainId)
    __, dataList= self.processPSSM( seqStr, seqStructMap, prefixExtended, pssmNameRaw, pssmNameProc, areSeqIdsMapped)

    if self.winSize:
      #Compute windows
      winPssmOutName= fNames[-1]
      badConsScore= PsiBlastManager.BAD_SCORE_CONSERVATION
      self.makeWindowed( dataList, ["pssm", "psfm","information"], [badConsScore]*3, [None]*3, winPssmOutName)
    return psioutName, pssmNameRaw

  def getPsiBlastFrom3DCons(self, psioutName, pssmNameRaw, seqStr, pdbCode, chainId):
    '''
      :param psioutName: str. Path to results file where aligments will be saved
      :param pssmNameRaw: str. Path to results file where pssms will be saved
      :param seqStr: str. Sequence of the chain
      :param pdbCode: str pdbId in order to query 3dCons. If None and no perfect match with seqStr,
                          psiblast will be launched directly
      :param chainId: str chainId to recover

      :return 0 if seq was available as it is in 3dCons (seq serch, res mapping needed),
              1 if pdbId_chainId available in 3dCons (aligment of seqId struct id included in pssms but not in aligments file),
              -1 if not available in 3dCons in any form
              2 if was not possible to compute it
    '''
    if not self.consDbManager.consDbIsAvailable():
      return -1    
    if not pdbCode is None:
      returnVal= self.consDbManager.retrieve3DConsFromPDBChain(pdbCode, chainId, pssmNameRaw, psioutName, uncompress= True)
      print("3dcons search status chain %s %d"%(chainId, returnVal))
      if returnVal==0:
        return 1
    returnVal= self.consDbManager.retrieve3DConsFromSeq(seqStr, pssmNameRaw, psioutName, uncompress= True)
    print("3dcons search status chain %s %d"%(chainId, returnVal))
    if returnVal==0:
      return 0
    elif returnVal in [-1, -2,-3]:
      return -1
    else:
      return 2

  def checkIfRawFilesComputed(self, psioutName, pssmNameRaw):
    psioutName= psioutName.replace("_*", "_\\*")
    pssmNameRaw= pssmNameRaw.replace("_*", "_\\*")
    cmd= "zgrep -e 'PSI Gapped' %s"%pssmNameRaw
    process= Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE)
    processOut= process.communicate()
    if len(processOut[1])>0 or len(processOut[0])==0:
      return False
      
    cmd= "zgrep -e 'Window for multiple hits:' %s"%psioutName
    process= Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE)
    processOut= process.communicate()
    if len(processOut[1])>0 or len(processOut[0])==0:
      return False
    return True

      
  def getPsiBlast(self, fastaInFname, psioutName, pssmNameRaw, seqStr, prefixExtended, pdbCode):
    '''
      Retrieves psiblast output (pssm and aligments). First, it tries to download from
      3dconsDb. If not available, launches psiblast command with fastaInFname as input file,
      psioutName as the output file that will save the aligments and pssmNameRaw as the output
      file that will save the pssms.
      :param fastaInFname: str. Path to fasta file where sequence is saved
      :param psioutName: str. Path to results file where aligments will be saved
      :param pssmNameRaw: str. Path to results file where pssms will be saved
      :param seqStr: str. Sequence of the chain
      :param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      :param pdbCode: str pdbId in order to query 3dCons. If None and no perfect match with seqStr,
                      psiblast will be launched directly

      :return seqIdsMapped: boolean. True if psiblast output is obtained from 3dConsDb and thus,
                structIds are included in first col of pssm files instead of seqIds, false if
                first column are seqIds.
    '''

    __, __, chainId = self.splitExtendedPrefix(prefixExtended)[:3]
    returnVal= self.getPsiBlastFrom3DCons( psioutName, pssmNameRaw, seqStr, pdbCode, chainId)
    if returnVal==-1:
      if self.checkIfRawFilesComputed(psioutName, pssmNameRaw):
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
      :param fastaInFname: str. Path to fasta file where sequence is saved
      :param psioutName: str. Path to results file where aligments will be saved
      :param pssmNameRaw: str. Path to results file where pssms will be saved
      :param seqStr: str. Sequence of the chain
    '''
    blastbin= self.psiBlastBin
    blastDB=  self.psiBlastDB
    assert os.path.isfile( self.psiBlastDB+".pal"), "Error, psiblast DB not found %s"%(self.psiBlastDB)
    numThr=   int(self.psiBlastNThrs)
    with open(pssmNameRaw, 'w') as f: #create pssm file in case error or no match happen. It will be override
      f.write("!ERROR\n"*(self.maxNumResiduesPartner+1))
    blastCMD = ("%(blastbin)s -query %(fastaInFname)s -num_iterations 3 -inclusion_ethresh 0.0001 -show_gis -db %(blastDB)s "+
                "-out %(psioutName)s -num_threads %(numThr)d -out_ascii_pssm %(pssmNameRaw)s -num_alignments 250  ")%locals()
    blastCMD= blastCMD.replace("_*", "_\\*")
    print(blastCMD)
    process= Popen( blastCMD, shell=True, stdout=PIPE, stderr=PIPE)
    processOut= process.communicate()
#    processOut=("","")
    #Check for failure
    with open(pssmNameRaw, 'r') as f:
      firstLine= f.readline()
    if len(processOut[1])>0 and processOut[1].startswith(b"BLAST engine error") or firstLine.startswith("!ERROR"): #No hits found will be dealt at processResults
      print("Error computing blast. Caught stdin/stderr:\n",processOut[0],processOut[1])
      with open(psioutName,"a") as f:
        f.write("!ERROR\n"+processOut[0]+"\n"+processOut[1])
      self.makeFakePSSM( pssmNameRaw, seqStr)
      
  def makeFakePSSM( self, pssmNameRaw, seqStr):
    header='''
Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts
            A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
'''
    oneFakeLine='''    %d %s     -1024  -1024  -1024  -1024  -1024  -1024  -1024  -1024  -1024   -1024  -1024  -1024  -1024  -1024  -1024   -1024   -1024  -1024  -1024   -1024    -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024   -1024  -1024  -1024   -1024   -1024  -1024  -1024.0 -1024.0\n'''
    with open(pssmNameRaw, "w") as f:
      f.write(header)
      for i, resName in enumerate(seqStr):
        f.write(oneFakeLine%(i, resName))
    return None

  def processPSSM(self, seq, seqStructMap, prefixExtended, pssmNameRaw, pssmNameProc, areSeqIdsMapped):
    '''
      Reads psiblast pssms output file and writes another one with tabulated format, headers and
      some error checking.
      :param seq: str. Sequence of the chain
      :param seqStructMap: ..manageSeqs.seqStructMap.seqStructMap
      :param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      :param pssmNameRaw: str. Path to psiblast aligments results
      :param pssmNameProc: str. Path where formated results will be saved.
      :param areSeqIdsMapped: boolean. True if psiblast output is obtained from 3dConsDb and thus,
                structIds are included in first col of pssm files instead of seqIds, false if
                first column are seqIds.
      :return pssmSeq: str
      :return dataList: [ (res1_Feats), (res2_Feats), ...]
            resI_Feats: ( (chainId, resId, resName), ([feat1], [feat2]) )
    '''
    prefix, chainType, chainId = self.splitExtendedPrefix(prefixExtended)[:3]
    try:
      pssmData, pssmResIds, pssmSeq = self.loadPSSM(pssmNameRaw)
      if areSeqIdsMapped: #3dcons db seqs are generally bigger as they have not 3d solved residues. Add them to compute windows
        seq= pssmSeq
        seqStructMap.addResiduesToSeqToStructMap( chainType, chainId, pssmSeq, pssmResIds)
      else:
        assert pssmSeq== seq, "Mismatch in 3DconsDB?\n%s\n%s"%(pssmSeq, seq)
    except IOError:
      print("Pssm %s was not computed. Default value inserted instead"%(pssmNameRaw))
      pssmSeq= seq
      pssmData= [ [PsiBlastManager.BAD_SCORE_CONSERVATION for i in range(42)] for i in range(len(seq))]
    dataList=[]
    listOfRowsToPrint=[]
    assert len(pssmData) == len(seq), "Error, mismtach in sequence downloaded from 3dcons"
    if NORMALIZE_IPP_AND_RWGRMP:
      pssmData= self.applyNormalization( pssmData, norm_vars=[40,41])
    for i, (pssmArray,letter) in enumerate(zip(pssmData,seq)):
      structIndex= seqStructMap.seqToStructIndex(chainType, chainId, i, asString= True)
      if structIndex:
        if self.filterOutLabels and structIndex[-1].isalpha():
          continue
      else:
        structIndex=str(i)+"?"
      full_resId_tuple= (chainId, structIndex, letter)
      rowRecord= list(full_resId_tuple)+ pssmArray
      dataList.append( (full_resId_tuple, (pssmArray[:20], pssmArray[20:40], pssmArray[40:])) )
      listOfRowsToPrint.append(" ".join(rowRecord))
      
    self.writeResultsFromDataDictSingleChain( {chainId: listOfRowsToPrint }, outName= self.getFNames(prefixExtended)[0])
    return pssmSeq, dataList

  def loadPSSM(self, fname):
    '''
      Loads a pssm file (psiblast format)
      :param fname: str. Path to pssm file.
      :return  list of lists. [ [row0_PSSM_values], [row1_PSSM_values]...]
    '''
    scores=[]
    indices=[]
    seq=""
    skip=2
    openFun= open
    if not os.path.isfile(fname) and os.path.isfile(fname+".gz"):
      fname+=".gz"
      openFun= gzip.open
    with openFun(fname) as f:
      for line in f:
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
          scores.append( lineArray[2:])
          indices.append(lineArray[0])
          seq+= lineArray[1]
    return scores, indices, seq
    

  def compressRawData(self, prefixExtended):
    for fname in os.listdir(self.pssmsPathRaw):
      if fname.startswith(prefixExtended) and not fname.endswith(".gz"):
        if os.path.islink(os.path.join(self.pssmsPathRaw, fname)):
          os.remove(os.path.join(self.pssmsPathRaw, fname))
        else:
          cmd= ["gzip", "-f", os.path.join(self.pssmsPathRaw, fname)]
          check_output(cmd)

#    for fname in os.listdir(self.psiblastPathOutPath): #al2co not working if using this
#      if fname.startswith(prefixExtended) and not fname.endswith(".gz"):
#        if os.path.islink(os.path.join(self.psiblastPathOutPath, fname)):
#          os.remove(os.path.join(self.psiblastPathOutPath, fname))
#        else:
#          cmd= ["gzip", "-f", os.path.join(self.psiblastPathOutPath, fname)]
#          check_output(cmd)

  def applyNormalization(self, listOfRows, norm_vars=[40,41]):
    '''
    :param listOfRows: [ residue1Vars, residue2Vars, ...]  residue_i_Vars: [0.3, 0.32, -1, ...]
    :param norm_vars: int[]. The variables to normalize
    :return listOfRows+[NormalizedVars]
    '''
    values= np.array([ [float(elem) for elem in getItemsFromList(norm_vars, row)] for row in listOfRows])
    meanVal = np.mean(values, axis=0)
    stdVal = np.std(values, axis=0)
    for row_i in range(len(listOfRows)):
      normFeats=[]
      for i,featNum in enumerate(norm_vars):
        if stdVal[i]==0:
          stdVal[i]= 1e-8
        normFeats+= [( float(listOfRows[row_i][featNum])-meanVal[i]) /stdVal[i]]
      listOfRows[row_i]+= [str(round(elem,5)) for elem in normFeats]
    return listOfRows
if __name__=="__main__":
  print("done")
