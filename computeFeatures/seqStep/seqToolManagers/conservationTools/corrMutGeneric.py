from __future__ import absolute_import, print_function
import os
from utils import myMakeDir, tryToRemove #utils is at the root of the package
import re
import time

from ...seqToolManager import SeqToolManager
                
class CorrMutGeneric(SeqToolManager):
  '''
    Computes corrMut and processes their outputs. Extends class ConservationManager. This is an
    abstract class
  '''
  BAD_CASE_SCORE = -1024  #Something went wrong tag
  SEQUENCES_SEPARATOR="-"*5
  MIN_N_SEQS_MSA= 10
  
  def __init__(self, computedFeatsRootDir, statusManager=None):
    '''
      :param computedFeatsRootDir: str. root path where results will be saved
      :param winSize: int>=1 or None. The size of the windows for sliding window if desired
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''  
    SeqToolManager.__init__(self, computedFeatsRootDir, None, statusManager)
    
    self.corrMutOutPath= myMakeDir(computedFeatsRootDir,"corrMut")
    self.chainsL={}
    self.chainsR={}
    
  def getFinalPath(self):
    '''
      returns path where results are saved
      :return corrMutOutPath: str
    '''
    return self.corrMutOutPath
    
  def getFNames(self, prefixExtended):
    '''
    Returns a dict that contains the fnames that will be used by corrMutManager
    :param prefixExtended. prefix for output fnames.
    :return Dict   {"corrMut":(corrMutProc, ), }
    '''
    assert len(self.chainsL)>0 and len(self.chainsR)>0
    
    prefix = self.splitExtendedPrefix(prefixExtended)[0]
    listOfCorrMutNames=[]
    for chainIdL in self.chainsL:
      for chainIdR in self.chainsR:
        listOfCorrMutNames.append(self.generateOutName(prefix, chainIdL, chainIdR) )
    return listOfCorrMutNames
    
  def generateOutName(self, prefix, chainIdL, chainIdR):
    return os.path.join(self.corrMutOutPath, prefix+".l-"+chainIdL+"_r-"+chainIdR+"_"+"u.corrMut.gz")
  
  def loadOneAligFile(self, aligFile):
  
    taxaRegex1= re.compile(".* OS=(([]'[\w]+.?\s){1,2})")
    taxaRegex2= re.compile(".* \[((\w+\s){1,2})")
    seqsByLine=set([])
    seqsByTaxa={}
#    print(len(seqsByTaxa), aligFile)
    with self.openForReadingFnameOrGz(aligFile) as f:
      f.readline()
      seq1= f.readline().strip()
      taxa= "targetProtein"
      seqsByLine.add( (seq1, taxa) )
      for line in f:
        if line.startswith(">"):
          matchObj= re.match(taxaRegex1, line)
          if matchObj:
            taxa= re.sub("[\[\]']","",matchObj.group(1).strip())
          else:
            matchObj= re.match(taxaRegex2, line)
            if matchObj:
              taxa= re.sub("[\[\]]","",matchObj.group(1).strip())
            else:
#              print(line)
#              print("taxa not found in aligFile %s"%aligFile)
              taxa=None
        else:
          line= line.strip()
          if not line in seqsByLine:
            seqsByLine.add( (line,taxa) )

    for seq, taxa in seqsByLine:
      aligQuality= len(seq)- seq.count("-")- sum( (aa.islower() for aa in seq))
      seq= re.sub(r"[a-z]", "", seq)
      if taxa not in seqsByTaxa:
        seqsByTaxa[taxa]= (seq, aligQuality)
      else:
        prevSeq, prevQuality= seqsByTaxa[taxa]
        if aligQuality> prevQuality:
          seqsByTaxa[taxa]= (seq, aligQuality)

#    print(len(seqsByTaxa), aligFile)
#    raw_input("enter")
    return seqsByTaxa
    
  def createPairedAlignmet(self, aligsDictL, aligsDictR, aligFormatedName):
    pairedAlig=""
    nPairs= 0

    if "targetProtein" not in aligsDictL or "targetProtein" not in aligsDictR: 
      return None
    seqL_initial, scoreL = aligsDictL["targetProtein"]
    seqR_initial, scoreR = aligsDictR["targetProtein"]
    pairedAlig+=("%s%s%s\n"%(seqL_initial, CorrMutGeneric.SEQUENCES_SEPARATOR, seqR_initial))
    for taxaL in aligsDictL:
#      print(taxaL)
      if taxaL!="targetProtein" and taxaL in aligsDictR:
        seqL, scoreL = aligsDictL[taxaL]
        seqR, scoreR = aligsDictR[taxaL]
        pairedAlig+=("%s%s%s\n"%(seqL, CorrMutGeneric.SEQUENCES_SEPARATOR, seqR))
        nPairs+=1
    with open(aligFormatedName,"w") as f:
      f.write(pairedAlig)
    return aligFormatedName, pairedAlig, nPairs, seqL_initial, seqR_initial
    
  def computeFromSeqStructMapper(self, seqStructMap, extendedPrefix, HHBlitsFnamesDict):
    '''
      Computes corrMut for the Multiple Sequence aligment hhBlitsOut after pairing it by taxa. If more than 2 sequences
      are found for one taxa, just best match is choosen
      :param seqStructMap: computeFeatures.seqStep.seqToolManagers.seqExtraction.SeqStructMapper
      :param HHBlitsFnamesDict: {"l":{"A":"1A2K_F0_l_C_.a3m"}, "r":{"B":"1A2K_F0_r_A_.a3m", "C":"1A2K_F0_r_B_.a3m"}}
    '''
    self.chainsL=set(HHBlitsFnamesDict["l"].keys())
    self.chainsR=set(HHBlitsFnamesDict["r"].keys())
    if self.checkAlreayComputed(extendedPrefix): 
      print("%s already computed correlated mutations for "%extendedPrefix)
      return

    prefix, __, chainType, chainId = self.splitExtendedPrefix(extendedPrefix)
    aligsDict= {chainType_:{ chainId_: self.loadOneAligFile(HHBlitsFnamesDict[chainType_][chainId_]) 
                    for chainId_ in HHBlitsFnamesDict[chainType_]} for chainType_ in HHBlitsFnamesDict}

    for chainIdL in aligsDict["l"]:
      for chainIdR in aligsDict["r"]:
        aligFormatedName= os.path.join(self.corrMutOutPath, "tmp_"+prefix+"_l-"+chainIdL+"-r-"+chainIdR+"_"+"u.ali")
        try:
          corrMutOutName= self.generateOutName(prefix, chainIdL, chainIdR)          
          if os.path.isfile(corrMutOutName) and self.getNLines(corrMutOutName)>self.minNumResiduesPartner: 
            print("%s already computed"%corrMutOutName)
            continue
          print("launching corrMut over chains %s - %s"%(chainIdL, chainIdR))
          aligOut= self.createPairedAlignmet(aligsDict["l"][chainIdL], aligsDict["r"][chainIdR], 
                                                                      aligFormatedName)
          if aligOut:
            __, __, nAlig, seqL, seqR= aligOut
          else:
            nAlig=0
          if nAlig> CorrMutGeneric.MIN_N_SEQS_MSA:
            startTime= time.time()
            iterOfCorrelatedRows= self.lauchCorrMutProgram(aligFormatedName)
            print("Time CorrMut", time.time()- startTime)
          else:
            iterOfCorrelatedRows= None #( "*** Sorry", "Error, not enough sequences in MSA") 
          self.processResults(iterOfCorrelatedRows, seqStructMap, chainIdL, chainIdR, corrMutOutName, nAlig)
        except (KeyboardInterrupt, Exception):
          print("Exception happend computing corrMut for %s over chains %s - %s"%(prefix, chainIdL, chainIdR))
          tryToRemove(corrMutOutName)
          raise
        finally:
          tryToRemove(aligFormatedName)
          pass
      
  def lauchCorrMutProgram(self, aligFormatedName):
    '''
    abstract method
    returns an iterator of [i,j, float(score)] where i and j are the indices of the residues whose score is included
    :return iterOfCorrelatedRows: iterator of elements as [res_i, res_j, corrMuScore] ] res_i and res_j are 0 based

    '''
    raise ValueError("Not implemented")
    
  def processResults(self, iterOfCorrelatedRows, seqStructMap, chainIdL, chainIdR, corrMutOutName, nAlig):

    '''
      Reads corrMut output file and writes another one with tabulated format, headers and
      some error checking.
      :param iterOfCorrelatedRows: iterator of [i,j, float(score)] where i and j are the indices of the residues whose score is included
      :param seqStructMap: ..manageSeqs.seqStructMap.seqStructMap
      :param corrMutOutName: str. Fname where formated results will be saved.
      :param iterOfCorrelatedRows: iterator of elements as [res_i, res_j, corrMuScore] ] res_i and res_j are 0 based
      :param chainIdL:str. The chain Id for the ligand
      :param chainIdR:str. The chain Id for the receptor
      :param nAlig: int. The number of rows of MSA
    '''

    seqL, __= seqStructMap.getSeq("l", chainIdL)
    seqR, __= seqStructMap.getSeq("r", chainIdR)
    seqStructMap.setCurrentSeq(seqL, "l", chainIdL)
    seqStructMap.setCurrentSeq(seqR, "r", chainIdR)
    lenSeqL= len(seqL)
    print("N sequences MSA: %d"%nAlig)
    corrMutQuality=  float(nAlig)/ (lenSeqL+len(seqR))
    addedI_J= set([])
    listOfRowsToPrint=[]

    if iterOfCorrelatedRows==None:
      notIncludedInResults= CorrMutGeneric.BAD_CASE_SCORE
    else:
      notIncludedInResults= 0.0
      lenSeparator= len(CorrMutGeneric.SEQUENCES_SEPARATOR)
      for line in iterOfCorrelatedRows:
        i, j, score= line
        if i>=lenSeqL or j <(lenSeqL+lenSeparator): continue
        j= j-lenSeqL-lenSeparator
        assert j>=0
        addedI_J.add((i,j))
        letterL= seqL[i]
        letterR= seqR[j]
        score= float(score)
        structIndexL= seqStructMap.seqToStructIndex("l", chainIdL, i, asString= True) 
        structIndexR= seqStructMap.seqToStructIndex("r", chainIdR, j, asString= True)
        if structIndexR is None or (self.filterOutLabels and structIndexR[-1].isalpha()): continue
        if structIndexL is None or (self.filterOutLabels and structIndexL[-1].isalpha()): continue
        listOfRowsToPrint.append("%s %s %s %s %s %s %f %f"%(chainIdL, structIndexL, letterL, chainIdR, structIndexR,
                                                              letterR, score, corrMutQuality))
    for i,letterL in enumerate(seqL):
      structIndexL= seqStructMap.seqToStructIndex("l", chainIdL, i, asString= True)
      if structIndexL is None or (self.filterOutLabels and structIndexL[-1].isalpha()): continue
      for j,letterR in enumerate(seqR):
        if not (i,j) in addedI_J:
          structIndexR= seqStructMap.seqToStructIndex("r", chainIdR, j, asString= True)
          if structIndexR is None or (self.filterOutLabels and structIndexR[-1].isalpha()): continue
          listOfRowsToPrint.append("%s %s %s %s %s %s %f %f"%(chainIdL, structIndexL, letterL, chainIdR, structIndexR,
                                                                letterR, notIncludedInResults, corrMutQuality))
    self.writeResultsFromDataDictPairL2R( listOfRowsToPrint, corrMutOutName, 
                                          featuresNames=[self.featName, self.featName+"Quality"])


