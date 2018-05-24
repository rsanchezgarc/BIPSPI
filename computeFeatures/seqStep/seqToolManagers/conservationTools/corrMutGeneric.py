from __future__ import absolute_import, print_function
import sys, os
from subprocess import Popen, PIPE, check_output
from Config import Configuration
from utils import myMakeDir, tryToRemove #utils is at the root of the package
import re
import time

from ..seqToolManager import seqToolManager, FeatureComputerException
from utils import myMakeDir, tryToRemove #utils is at the root of the package
                
class CorrMutGeneric(seqToolManager):
  '''
    Computes corrMut and processes their outputs. Extends class ConservationManager. This is an
    abstract class
  '''
  BAD_SCORE_CONSERVATION = -1048576  #Something went wrong tag
  SEQUENCES_SEPARATOR="-"*5
  MIN_N_SEQS_MSA= 10
  def __init__(self, seqsManager, outPath):
    '''
      @param seqsManager: ..manageSeqs.seqsManager.SeqsManager 
      @param outPath: str. path where corrMut scores will be saved
    '''  
    seqToolManager.__init__(self,seqsManager, outPath, None)
    
    self.seqsManager= seqsManager
    self.corrMutOutPath= myMakeDir(outPath,"corrMut")
    
  def getFinalPath(self):
    '''
      returns path where results are saved
      @return corrMutOutPath: str
    '''
    return self.corrMutOutPath
    
    
  def checkAlreayComputed(self, corrMutOutName):
    '''
    
      checks if output files have been computed
      @param corrMutOutName. Output file
      @return boolean. True if Output file sequences have been computed, False otherwise
    '''    

    if not os.path.isfile(corrMutOutName):
      return False
    else:
      nlines=int( check_output('wc -l {}'.format(corrMutOutName), shell=True).split()[0])
      if nlines<12: #Too small file to be correct
        return False
    return True

  def loadOneAligFile(self, aligFile):
  
    taxaRegex1= re.compile(".* OS=(([]'[\w]+.?\s){1,2})")
    taxaRegex2= re.compile(".* \[((\w+\s){1,2})")
    seqsByLine=set([])
    seqsByTaxa={}
    if not os.path.isfile(aligFile):
      return seqsByTaxa
    with open(aligFile) as f:
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
              print(line)
              print("taxa not found in aligFile %s"%aligFile)
              taxa=None
#              raw_input("continue?")
#              raise ValueError("taxa not found in aligFile %s"%aligFile)
#          print(taxa)
#          raw_input("press enter to continue")
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
    return seqsByTaxa
    
  def createPairedAlignmet(self, aligsDictL, aligsDictR, aligFormatedName):
    pairedAlig=""
    nPairs= 0
#    print(sorted(aligsDictL.keys()))
#    raw_input("press enter")
    if "targetProtein" not in aligsDictL or "targetProtein" not in aligsDictR: 
      return None
    seqL_initial, scoreL = aligsDictL["targetProtein"]
    seqR_initial, scoreR = aligsDictR["targetProtein"]
    pairedAlig+=("%s%s%s\n"%(seqL_initial, CorrMutGeneric.SEQUENCES_SEPARATOR, seqR_initial))
    for taxaL in aligsDictL:
      if taxaL!="targetProtein" and taxaL in aligsDictR:
        seqL, scoreL = aligsDictL[taxaL]
        seqR, scoreR = aligsDictR[taxaL]
        pairedAlig+=("%s%s%s\n"%(seqL, CorrMutGeneric.SEQUENCES_SEPARATOR, seqR))
        nPairs+=1
    with open(aligFormatedName,"w") as f:
      f.write(pairedAlig)
      
    return aligFormatedName, pairedAlig, nPairs, seqL_initial, seqR_initial
    
  def compute(self, HHBlitsFnamesDict, prefix):
    '''
      Computes corrMut for the Multiple Sequence aligment hhBlitsOut after pairing it by taxa. If more than 2 sequences
      are found for one taxa, just best match is choosen
      @param HHBlitsFnamesDict: {"l":{"A":"1A2K_l_A_u.a3m"}, "r":{"B":"1A2K_r_B_u.a3m", "C":"1A2K_r_C_u.a3m"}}
      @param prefix: str. The prefix of the complex, p.e. 1A2K
    '''
    
    aligsDict= {chainType:{ chainId: self.loadOneAligFile(HHBlitsFnamesDict[chainType][chainId]) 
                    for chainId in HHBlitsFnamesDict[chainType]} for chainType in HHBlitsFnamesDict}
    for chainIdL in aligsDict["l"]:
      for chainIdR in aligsDict["r"]:
        print("launching corrMut over chains %s - %s"%(chainIdL, chainIdR))
#        raw_input("press enter to procced")
        aligFormatedName= os.path.join(self.corrMutOutPath, "tmp_"+prefix+"_l-"+chainIdL+"-r-"+chainIdR+"_"+"u.ali")
        try:
          corrMutOutName= os.path.join(self.corrMutOutPath, prefix+"_l-"+chainIdL+"_r-"+chainIdR+"_"+"u.corrMut")
          if self.checkAlreayComputed(corrMutOutName): 
            print("%s already computed"%corrMutOutName)
            continue
          aligOut= self.createPairedAlignmet(aligsDict["l"][chainIdL], aligsDict["r"][chainIdR], 
                                                                      aligFormatedName)
          if aligOut:
            __, __, nAlig, seqL, seqR= aligOut
          else:
            nAlig=0
            seqL, __= self.seqsManager.getSeq("l", chainIdL)
            seqR, __= self.seqsManager.getSeq("r", chainIdR)

          if nAlig> CorrMutGeneric.MIN_N_SEQS_MSA:
            startTime= time.time()
            iterOfCorrelatedRows= self.lauchCorrMutProgram(aligFormatedName)
            print("Time CorrMut", time.time()- startTime)
          else:
            iterOfCorrelatedRows= None #( "*** Sorry", "Error, not enough sequences in MSA") 
#          if len(processOut[1])>0:
#            print("Error computing corrMut. Caught stdin/stderr:\n",processOut[0],processOut[1])
          self.saveProcResults(seqL, seqR, corrMutOutName, iterOfCorrelatedRows, chainIdL, chainIdR, nAlig)
        except (KeyboardInterrupt, Exception):
          print("Exception happend computing corrMut for %s over chains %s - %s"%(prefix, chainIdL, chainIdR))
          tryToRemove(corrMutOutName)
          raise
        finally:
          tryToRemove(aligFormatedName)
          pass
      
  def lauchCorrMutProgram(self, aligFormatedName):
    #abstract Method
    return None
    
  def writeHeader(self, fHandler):
    fHandler.write("chainIdL structResIdL resNameL chainIdR structResIdR resNameR %s %s\n"%(self.featName, 
                  self.featName+"Quality"))
                  
  def saveProcResults(self, seqL, seqR, corrMutOutName, iterOfCorrelatedRows, chainIdL, chainIdR, nAlig):
    '''
      Reads corrMut output file and writes another one with tabulated format, headers and
      some error checking.
      @param: seqL: str. Sequence of the ligand chain
      @param: seqR: str. Sequence of the receptor chain
      @param corrMutOutName: str. Fname where formated results will be saved.
      @param iterOfCorrelatedRows: iterator of elements as [res_i, res_j, corrMuScore] ] res_i and res_j are 0 based
      @param chainIdL:str. The chain Id for the ligand
      @param chainIdR:str. The chain Id for the receptor
      @param nAlig: int. The number of rows of MSA
    '''

    corrMutQuality=  float(nAlig)/ (len(seqL)+len(seqR))
    if iterOfCorrelatedRows==None:
      self.makeFakeFile( seqL, seqR, corrMutOutName, corrMutQuality, chainIdL, chainIdR)
      return 1
    else:
      try:
        with open(corrMutOutName,"w") as outFile:
          self.writeHeader(outFile)
          scoresDict={}
          lenSeqL= len(seqL)
          lenSeparator= len(CorrMutGeneric.SEQUENCES_SEPARATOR)
          addedI_J= set([])
#          for line in corrMutOut.split("\n")[1:]:
          for line in iterOfCorrelatedRows:
            i, j, score= line
#            i, j=int(i)-1, int(j)-1
            if i>=lenSeqL or j <(lenSeqL+lenSeparator): continue
            j= j-lenSeqL-lenSeparator
            assert j>=0
            addedI_J.add((i,j))
            letterL= seqL[i]
            letterR= seqR[j]
            score= float(score)
            structIndexL= self.seqsManager.seqToStructIndex("l", chainIdL, i, asString= True) 
            structIndexR= self.seqsManager.seqToStructIndex("r", chainIdR, j, asString= True)
            if structIndexR is None or (self.filterOutLabels and structIndexR[-1].isalpha()): continue
            if structIndexL is None or (self.filterOutLabels and structIndexL[-1].isalpha()): continue
            outFile.write("%s %s %s %s %s %s %f %f\n"%(chainIdL, structIndexL, letterL, chainIdR, structIndexR, letterR, 
                                                  score, corrMutQuality))
          for i in range(len(seqL)):
              letterL= seqL[i]            
              for j in range(len(seqR)):
                if not (i,j) in addedI_J:
                  letterR= seqR[j]
                  structIndexL= self.seqsManager.seqToStructIndex("l", chainIdL, i, asString= True) 
                  structIndexR= self.seqsManager.seqToStructIndex("r", chainIdR, j, asString= True)
                  if structIndexR is None or (self.filterOutLabels and structIndexR[-1].isalpha()): continue
                  if structIndexL is None or (self.filterOutLabels and structIndexL[-1].isalpha()): continue
                  outFile.write("%s %s %s %s %s %s %f %f\n"%(chainIdL, structIndexL, letterL, chainIdR, structIndexR, letterR, 
                                                    0.0, corrMutQuality))
          return 0
      except (KeyboardInterrupt, Exception) as e:
        print(e)
        print("Exception happend computing %s"%corrMutOutName)
        tryToRemove(corrMutOutName)    
        raise

  def makeFakeFile(self, seqL, seqR, corrMutOutName, corrMutQuality, chainIdL, chainIdR):
    try:
      with open(corrMutOutName,"w") as outFile:
        self.writeHeader(outFile)
        for i,letterL in enumerate(seqL):
          structIndexL= self.seqsManager.seqToStructIndex("l", chainIdL, i, asString= True)
          if structIndexL is None or (self.filterOutLabels and structIndexL[-1].isalpha()): continue
          for j,letterR in enumerate(seqR):
            structIndexR= self.seqsManager.seqToStructIndex("r", chainIdR, j, asString= True)
            if structIndexR is None or (self.filterOutLabels and structIndexR[-1].isalpha()): continue
            outFile.write("%s %s %s %s %s %s %f %f\n"%(chainIdL, structIndexL, letterL, chainIdR, structIndexR, letterR, 
                                                  CorrMutGeneric.BAD_SCORE_CONSERVATION, corrMutQuality))
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing %s"%corrMutOutName)
      tryToRemove(corrMutOutName)    
      raise

