from __future__ import absolute_import, print_function
import os
import numpy as np
from subprocess import Popen, PIPE
from Bio.PDB.Polypeptide import aa1 as AA_STANDARD

from ....featuresComputer import FeatureComputerException
from ...seqToolManager import SeqToolManager
from .al2coWorkers.parsePsiBlast import parsePsiBlast
from utils import myMakeDir, tryToRemove

class Al2coManager(SeqToolManager):
  '''
    Computes al2co and processes their outputs. Extends class seqToolManager
  '''
  VAR_LIST= ["al2coScore", "al2coScoreNorm"]
  BAD_SCORE_CONSERVATION = "-1048576"  #Something went wrong tag
  def __init__(self, computedFeatsRootDir, winSize=None, statusManager=None):
    '''
      :param computedFeatsRootDir: str. root path where results will be saved
      :param winSize: int>=1 or None. The size of the windows for sliding window if desired
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''
    SeqToolManager.__init__(self, computedFeatsRootDir, winSize)
    
    self.al2coOutPath= myMakeDir(self.computedFeatsRootDir,"al2co")
    if winSize:
      self.al2coPathWindowed= myMakeDir(self.computedFeatsRootDir,"al2co_wSize"+str(winSize))
    else:
      self.al2coPathWindowed= None
    
  def getFinalPath(self):
    '''
      returns path where results are saved
      :return al2coOutPath: str
    '''
    return self.al2coOutPath
    
  def getFNames(self, prefixExtended):
    '''
      Returns a dict that contains the fnames that will be used  by al2co
      :param prefixExtended. prefix for output fnames.
      :return list of fnames: [ fname1, fnam2, ...]
    '''    
    al2coProc= os.path.join(self.al2coOutPath, prefixExtended+".al2co.gz")
    fNames=[al2coProc]
    if not self.winSize is None:
      al2coWindowedOutName= os.path.join(self.al2coPathWindowed, prefixExtended+".wsize"+str(self.winSize)+".al2co.gz")
      fNames+= [al2coWindowedOutName]
    return fNames
    
    
  def computeFromSeqStructMapper(self, seqStructMap, prefixExtended, psiblastOutName, pssmOutNameRaw):
    '''
      Computes al2co for the sequence seqStr, that is contained at fastaInFname. This sequence is
      associated with prefixExtended as an unambiguous id
      :param seqStructMap: computeFeatures.seqStep.seqToolManagers.seqExtraction.SeqStructMapper
      :param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      :param psiblastOutName: str. Path to psiblast aligments results
      :param pssmOutNameRaw: str. Path to psiblast pssms results
    '''
    msaFname= None

    prefix, chainType, chainId= self.splitExtendedPrefix(prefixExtended)[:3]
    seqStr, fastaFname= seqStructMap.getSeq(chainType, chainId) # repeat as psiBlastManager can modify seqs
    seqStructMap.setCurrentSeq(seqStr, chainType, chainId)
    if self.checkAlreayComputed(prefixExtended):
      print("Al2co already computed for %s"%prefixExtended)
      return 0
    fNames= self.getFNames(prefixExtended)
    print("launching al2co over %s"%prefixExtended)      
    al2coProcName= fNames[0]
    al2coRawName= os.path.join(self.al2coOutPath, prefixExtended+".fasta.csv")
    try:
      if os.path.isfile(psiblastOutName):
        alignedSeqsDict= parsePsiBlast( inputSeq=seqStr, psiBlastOut=psiblastOutName)

        filteredSeqsFname= self.runCdHit(alignedSeqsDict, inputSeq=seqStr, psiBlastOut=psiblastOutName)
        msaFname= self.runClustalW(filteredSeqsFname, psiBlastOut=psiblastOutName)

        cmd= [self.al2coBin, "-i", msaFname,"-m", "0", "-f", "2", "-a", "F", "-b", "50",
              "-g", "0.50", "-w", "1", "-c", "0", "-o", al2coRawName, "-t", al2coProcName]

        print(" ".join(cmd))
        process= Popen(cmd, stdout=PIPE, stderr=PIPE)
        processOut= process.communicate()
        if len(processOut[1])>0:
          print("Error computing al2co. Caught stdin/stderr:\n",processOut[0],processOut[1])
      else:
        print("Error computing al2co. Psiout does not exists for %s"%(prefixExtended))
        al2coRawName=None

      dataList= self.processAl2co(seqStr, seqStructMap, prefixExtended, al2coRawName, al2coProcName)
      if self.winSize:
        self.makeWindowed( dataList, ["al2co", "al2coNorm"], [Al2coManager.BAD_SCORE_CONSERVATION]*2, [None]*2, 
                           fNames[1])
    except (Exception, KeyboardInterrupt):
      self.tryToRemoveAllFnames(prefixExtended)
      raise
    finally:
      if msaFname: tryToRemove(msaFname)

  def processAl2co(self, seq, seqStructMap, prefixExtended, al2coRaw, al2coProc):
    '''
      Reads al2co output file and writes another one with tabulated format, headers and
      some error checking.
      :param: seq: str. Sequence of the chain
      :param prefixExtended: str. unambiguous id of the sequence that will be the prefix of output names
      :param al2coRaw: str. Path to al2co results
      :param al2coProc: str. Path where formatted results will be saved.
    '''
    if al2coRaw is None:
      conserData = [(letter, Al2coManager.BAD_SCORE_CONSERVATION) for letter in seq]
    else:
      try:
        conserData = self.loadRawAl2co(al2coRaw)
      except IOError:
        conserData= [ (letter, Al2coManager.BAD_SCORE_CONSERVATION) for letter in seq]

    prefix, chainType, chainId= self.splitExtendedPrefix(prefixExtended)[:3]
#    print(len(conserData)); raw_input("enter")
    try:
      alcoIx=0
      seqIx=0
      seqLen= len(seq)
      letters, conserVals = zip(* conserData)
      conserVals= [float(elem) for elem in conserVals]
      alcoLen= len(conserData)
      dataList=[]
      listOfRowsToPrint=[]
      mean_val= np.mean(conserVals)
      std_val= np.std(conserVals)

      while seqIx<seqLen and alcoIx<alcoLen:
        letter= seq[seqIx]
        letterAl2co, consVal= conserData[alcoIx]
        if letterAl2co== letter or (letterAl2co=="-" and letter=="X"):
          structIndex= seqStructMap.seqToStructIndex(chainType, chainId, seqIx, asString= True)
#          print(seqIx, letter, alcoIx, structIndex)
          if structIndex:
            if self.filterOutLabels and structIndex[-1].isalpha():
              continue
          else:
            structIndex=str(seqIx)+"?"
          if std_val!=0:
           consValNormalized= (float(consVal)- mean_val)/std_val
          else:
           consValNormalized=float(consVal)
          dataList.append( ( (chainId, structIndex,letter), ( [consVal], [str(consValNormalized)],) ) )
          listOfRowsToPrint.append( "%s %s %s %s %s"%( chainId, structIndex, letter, consVal, consValNormalized) )
          alcoIx+=1
          seqIx+=1
        elif not letter in AA_STANDARD and letterAl2co=="-":
          alcoIx+=1
          seqIx+=1
        elif letterAl2co=="-":
          alcoIx+=1
        else:
          print(conserData)
          print(alcoIx, seqIx)
          raise ValueError("Al2co mismatch %s %s "%(letterAl2co, letter))
#      print(len(listOfRowsToPrint)); raw_input("enter to continue")
      self.writeResultsFromDataDictSingleChain( {chainId: listOfRowsToPrint }, outName= al2coProc)
      return dataList
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing %s"%al2coProc)
      tryToRemove(al2coProc)
      raise
    finally:
      if al2coRaw is not None:
        tryToRemove(al2coRaw)
      pass
    
  def loadRawAl2co(self, filename):
    '''
      Loads an al2co file
      :param fname: str. Path to al2co file.
      :return  list of strings. ["row0_Al2co","row1Al2co"...]
    '''  
    conserv= []
    for line in open(filename):
      lineArray=line.split()
      if lineArray[0][0].isdigit():
        conserv.append(lineArray[1:3])
      else:
        break
    return conserv


  def runCdHit(self, allHits, inputSeq, psiBlastOut, pairSeqIdThr=0.95):
    tmpName= os.path.basename(psiBlastOut).split(".")[0]
    tmpName= os.path.join(self.tmp, tmpName)
    cdhitInName= tmpName+".in-cdhit"
    cdhitOutName= tmpName+".out-cdhit"
    try:
      with open(cdhitInName, "w") as f:
        for hit in allHits:
          f.write("> %s\n"%(hit["target_full_id"]))
          f.write("%s\n"%(hit["targetSeq"].replace("-","")) )

      if(pairSeqIdThr > .70 and pairSeqIdThr <= 1.00): n=5
      elif (pairSeqIdThr <= .70 and pairSeqIdThr >= .55): n=4
      elif (pairSeqIdThr < .55 and pairSeqIdThr >= .50): n=3
      elif (pairSeqIdThr < .50 and pairSeqIdThr >= .40): n=2
      else: raise ValueError("Error, just .4<=pairSeqIdThr<=1.00 allowed")
      
      cdhitCmd= [self.cdHitBin, "-i", cdhitInName, "-o", cdhitOutName, "-n", str(n), 
                 "-c", str(pairSeqIdThr), "-T", str(self.psiBlastNThrs)]
      print(" ".join(cdhitCmd))
      proc = Popen(cdhitCmd, stdin= PIPE, stdout=PIPE, stderr=PIPE)
      output=  proc.communicate()
      if output== None or output[1]!="" or "There was an error cd-hit psiblast" in output[0]:
        print(output)
        print ("Error when parsing %s for al2Co"%psiBlastOut)
        raise FeatureComputerException("Error when cd-hit %s for al2Co"%psiBlastOut)
        
      with open(cdhitOutName, "r+") as f:
        fileData = f.read()
        f.seek(0, 0)
        f.write("> InputSeq\n")
        f.write("%s\n"%(inputSeq.replace("-","")) )
        f.write(fileData+"\n")
      return cdhitOutName
    except (Exception, KeyboardInterrupt):
      tryToRemove(cdhitOutName)
      raise
    finally:
      tryToRemove(cdhitInName)
      
  def runClustalW(self, filteredSeqsFname, psiBlastOut, clustalWOutName=None):
    tmpFnameCommon= ".".join(filteredSeqsFname.split(".")[:-1])
    if clustalWOutName is None:
      clustalWOutName= tmpFnameCommon+".clustalw"
    clustalCommand=[self.clustalW, "-infile=%s"%filteredSeqsFname, "-outfile=%s"%clustalWOutName, "-outorder=INPUT"]
    print(" ".join(clustalCommand))
    try :
      proc = Popen(clustalCommand, stdin= PIPE, stdout=PIPE, stderr=PIPE)
      output=  proc.communicate()
      if output== None or output[1]!="" or "There was an error parsing psiblast, clustalw" in output[0]:
        print(output)
        print ("Error when clustalw %s for al2Co"%psiBlastOut)
        raise FeatureComputerException("Error when clustalw %s for al2Co"%psiBlastOut)
      return clustalWOutName
    except (Exception, KeyboardInterrupt):
      tryToRemove(clustalWOutName)
      raise
    finally:
      tryToRemove(filteredSeqsFname)
      tryToRemove(filteredSeqsFname+".clstr")
      tryToRemove( tmpFnameCommon+".dnd")
