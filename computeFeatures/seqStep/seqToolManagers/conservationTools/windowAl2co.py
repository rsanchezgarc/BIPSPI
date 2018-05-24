'''
This module is no longer used
'''
from __future__ import absolute_import
import os
import numpy as np
import itertools

FEATURES_TO_INCLUDE=["score"]

class windowAl2co(object):

  def __init__(self, winSize):
    self.winSize= winSize
    self.noProfileWin=[ 0.0 for feat in FEATURES_TO_INCLUDE]
      
  def compute(self,inName,outName):

    prefix, chainType, f_chain, rest = os.path.split(inName)[-1].split("_")
    winSizeStep= int(self.winSize/2)

    profile=[]
    letras=[]
    ids=[]

    f=open(inName)
    header=f.readline().split()
    resIdIndex= header.index("structResId")
    letterIndex= header.index("resName")
    featuresIndices= []
    for featName in FEATURES_TO_INCLUDE:
      featuresIndices.append(header.index(featName))
      
    for line in f:
##          print( [line])
      lineArray= line.split()
      if len(lineArray)!=0:
        resInd=lineArray[resIdIndex]
        letras.append(lineArray[letterIndex])
#        print(lineArray, len(lineArray), firstValueIndex, firstValueIndex+self.numProfileElem)
#        print (lineArray[firstValueIndex:(firstValueIndex+self.numProfileElem)])
        profile.append([float(lineArray[i]) for i in featuresIndices])
        ids.append((f_chain,resInd))
    f.close()
        

    numElements= len(ids)
    myProfile={}
    myLetters={}
    mean_of_profile={}
    for resNum,resId in enumerate(ids):
  
      myProfile[resId]=[  profile[i] if (i>=0 and i<= numElements-1) else self.noProfileWin
                for i in range(resNum-winSizeStep, resNum+winSizeStep+1)]
      myLetters[resId]= [letras[resNum]]
      mean_of_profile[resId]= np.mean(myProfile[resId])


    mean_of_means_of_profile= np.mean(list(mean_of_profile.values()))
    sigma_of_means_of_profile= np.std(list(mean_of_profile.values()))
    
    outFile= open(outName,"w")
    for resId in sorted(ids):

      header= ["chainId structResId resName"] + list(itertools.chain.from_iterable
                                                              ([FEATURES_TO_INCLUDE for i in range(self.winSize)])) + \
                                            ["winAverage", "winAverageNormalized"]

      outFile.write( " ".join(header)+"\n")

      out= list(resId)+myLetters[resId]+list(itertools.chain.from_iterable(myProfile[resId])) + [mean_of_profile[resId]] + \
            [ (mean_of_profile[resId] - mean_of_means_of_profile)/ sigma_of_means_of_profile ]
      out= [str(val) for val in out]
      outFile.write( " ".join(out)+"\n")
      break

    for resId in sorted(ids):
      out= list(resId)+myLetters[resId]+list(itertools.chain.from_iterable(myProfile[resId])) + [mean_of_profile[resId]] + \
            [ (mean_of_profile[resId] - mean_of_means_of_profile)/ sigma_of_means_of_profile ]
      out= [str(val) for val in out]
      outFile.write( " ".join(out)+"\n")
    outFile.close()

 
def windowAllAl2co(inputDir, outPutDir, wsize):
  inputDir= os.path.expanduser(inputDir)
  outPutDir= os.path.expanduser(outPutDir)
  if not os.path.isdir(outPutDir):
    os.makedirs(outPutDir)
  for fname in os.listdir(inputDir):
    print( fname)
    inName= os.path.join(inputDir, fname)
    prefixAndChainInfo= fname.split(".")[0]
    prefixAndChainInfo+= ".wsize"+str(wsize)+".winAl2co.tab"
    outName= os.path.join(outPutDir, prefixAndChainInfo)
    windowAl2co(wsize).compute(inName,outName)
    
  return None

def test():

#  inName= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/computedFeatures/seqStep/conservation/al2co/1A2K_l_A_u.al2co"
#  outName= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/computedFeatures/seqStep/conservation/winAl2co/1A2K_l_A_u.al2coWinPrueba"
#  windowAl2co(5).compute(inName,outName)

  inputDir= "~/Tesis/rriPredMethod/data/ppdockingBenchData/computedFeatures/seqStep/conservation/al2co/"
  outPutDir= "~/Tesis/rriPredMethod/data/ppdockingBenchData/computedFeatures/seqStep/conservation/winAl2co"

  windowAllAl2co(inputDir, outPutDir, wsize= 11)

if __name__== "__main__":

  test()

