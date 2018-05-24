from __future__ import absolute_import
import os
import numpy
import itertools

AA_CODE_ELEMENTS= ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","Z","X"] #Z stands for no residue, X non standard
class WindowPSSM(object):
  '''
  This class computes sliding windows at the sequence level using as features psiblast derived ones and aminoacid codes.
  '''
  AA_CODE_ELEMENTS= AA_CODE_ELEMENTS
  NonAminoacidLetter='Z'
  NonAminoacidValuePSSM= -1024
  NonAminoacidValueEntropy= 0.0

  def __init__(self, winSize, normalizeSeqEntr= True, includePSSM= True, includePSFM= True):
    '''
      @param winSize: int. The size that sliding window will have
      @param normalizeSeqEntr: booelan. whether seq entropy will be normalized at chain level or not
      @param includePSSM: booelan. whether PSSM will be included or not
      @param includePSFM: booelan. whether PDFM will be included or not  
    '''  
    self.winSize= winSize
    self.centralPos= winSize//2
    self.normalize= normalizeSeqEntr
    
    self.profileStarts=0
    if includePSSM:
      self.numProfileElem= 20
    else:
      self.numProfileElem= 0
    if includePSFM:
      self.numProfileElem+= 20
      if not includePSSM:
        self.profileStarts= 20
      
    self.numSeqEntro= 2
      
  def compute(self,pssmFname,outWindowsFile):
    '''
      Computes sliding windows for one sequence
      @param pssmFname: str. path to pssms (processed) results for one sequence
      @param outWindowsFile: str. name of output file
    '''
    prefix, chainType, f_chain, rest = os.path.split(pssmFname)[-1].split("_")
    winSizeStep= int(self.winSize/2)

    pssmProfile=[]
    entropyProfile=[]
    letras=[]
    ids=[]
    noAA_PSSMprofile=[WindowPSSM.NonAminoacidValuePSSM for i in range(self.numProfileElem)] 
    noAA_SeqEntropyProfile=[WindowPSSM.NonAminoacidValueEntropy for i in range(self.numSeqEntro)]
    f=open(pssmFname)
    header=f.readline().split()
    resIdIndex= header.index("structResId")
    letterIndex= header.index("resName")
    firstValueIndex= header.index("pssm") + self.profileStarts
    for line in f:
      lineArray= line.split()
      if len(lineArray)!=0:
        resInd=lineArray[resIdIndex]
        letras.append(lineArray[letterIndex])
        pssmProfile.append([int(val) for val in lineArray[firstValueIndex:(firstValueIndex+self.numProfileElem)]]) 
        entropyProfile.append([float(val) for val in lineArray[-2:]])
        ids.append((f_chain,resInd))
    f.close()

    numElements= len(ids)
    myLetras={}
    myPSSMs={}
    myEntropy={}
    for resNum,resId in enumerate(ids):

      myLetras[resId]=[  letras[i] 
                  if (i>=0 and i<= numElements-1) 
                  else WindowPSSM.NonAminoacidLetter
                  for i in range(resNum-winSizeStep, resNum+winSizeStep+1)]

      myPSSMs[resId]=[  pssmProfile[i] if (i>=0 and i<= numElements-1) else noAA_PSSMprofile
                for i in range(resNum-winSizeStep, resNum+winSizeStep+1)]


      myEntropy[resId]=[  entropyProfile[i] if (i>=0 and i<= numElements-1) else noAA_SeqEntropyProfile
                for i in range(resNum-winSizeStep, resNum+winSizeStep+1)]

    if self.normalize:
      myEntropy= self.applyNormalization(myEntropy)

    LettersStartPoint= 3
    outFile= open(outWindowsFile,"w")
    for resId in sorted(ids):

      levelsTags= [str(LettersStartPoint+ i)+ ":"+";".join(WindowPSSM.AA_CODE_ELEMENTS) for i in range(len(myLetras[resId]))]

      outFile.write( "#Levels: " + " ".join(levelsTags)+"\n")

      header= ["chainId structResId resName"]+["resNameWin" for i in myLetras[resId]]+ \
              ["seqEntropy" for i in  itertools.chain.from_iterable(myEntropy[resId])]+ \
              ["pssm" for i in itertools.chain.from_iterable(myPSSMs[resId])]

      outFile.write( " ".join(header)+"\n")

      out= list(resId)+list(myLetras[resId][self.centralPos])+myLetras[resId]+list(itertools.chain.from_iterable(myEntropy[resId]))+ list(itertools.chain.from_iterable(myPSSMs[resId]))
      out= [str(val) for val in out]
      outFile.write( " ".join(out)+"\n")
      break

    for resId in sorted(ids)[1:]:
      out= list(resId)+list(myLetras[resId][self.centralPos])+myLetras[resId]+list(itertools.chain.from_iterable(myEntropy[resId]))+ list(itertools.chain.from_iterable(myPSSMs[resId]))
      out= [str(val) for val in out]
      outFile.write( " ".join(out)+"\n")
    outFile.close()



  def applyNormalization(self,res):
    '''
      Normalize seq entropy
      @param res: dict of list  {resId: [ [0entropyVal_0, 1entropyVal_0], [0entropyVal_1, 1entropyVal_1]...]
      return  Same dict but mean std normalized 
    '''
    values= numpy.array([val for valueProfiles in res.values() for val in valueProfiles ])

    meanVal = numpy.mean(values, axis=0)
    stdVal = numpy.std(values, axis=0)

    for featNum in range(len(stdVal)):
      if stdVal[featNum]==0:
        stdVal[featNum]= 1e-8
      for chainResId in res:
        for i in range(self.winSize):
          res[chainResId][i] =  res[chainResId][i] + [(res[chainResId][i][featNum]-meanVal[featNum]) /stdVal[featNum]]

    return res

 
def windowAllPSSMs(inputDir, outPutDir, wsize, normalizeSeqEntr= True, includePSSM= True, includePSFM= False):
  '''
  Computes windowed pssms for all pssms files in path imnputDir. Pssms must be in processed format
  '''
  inputDir= os.path.expanduser(inputDir)
  outPutDir= os.path.expanduser(outPutDir)
  if not os.path.isdir(outPutDir):
    os.makedirs(outPutDir)
  for fname in os.listdir(inputDir):
    print( fname)
    pssmFname= os.path.join(inputDir, fname)
    prefixAndChainInfo= fname.split(".")[0]
    prefixAndChainInfo+= ".wsize"+str(wsize)+".pssm"
    outWindowsFile= os.path.join(outPutDir, prefixAndChainInfo)
    WindowPSSM(wsize, normalizeSeqEntr,includePSSM, includePSFM).compute(pssmFname,outWindowsFile)
    
  return None

def test():
#  pssmFname= "/home/rsanchez/Desktop/computedFeatures/seqStep/conservation/pssms/procPssms/1A2K_l_A_u.pssm"
#  outWindowsFile= "/home/rsanchez/Tesis/rriPredMethod/results/resultsPrueba/1A2K_l_A_u.pssmWinPrueba"
#  WindowPSSM(5, normalize= True).compute(pssmFname,outWindowsFile)

  inputDir= "~/Tesis/rriPredMethod/data/ppdockingBenchData/computedFeatures/seqStep/conservation/pssms/procPssms"
  outPutDir= "~/Tesis/rriPredMethod/data/ppdockingBenchData/computedFeatures/seqStep/conservation/pssms/windowedPSSMs/wSize11"

  windowAllPSSMs(inputDir, outPutDir, wsize= 11, normalizeSeqEntr= True, includePSSM= True, includePSFM= True)

if __name__== "__main__":

  test()

