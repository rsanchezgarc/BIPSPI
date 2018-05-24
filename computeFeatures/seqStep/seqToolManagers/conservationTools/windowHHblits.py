from __future__ import absolute_import
import os
import numpy
import itertools

class WindowHHblits(object):
  '''
  This class computes sliding windows at the sequence level using as features psiblast derived ones and aminoacid codes.
  '''
  NonAminoacidValueProfile= -1024
  numProfileElem= 31
  def __init__(self, winSize):
    '''
      @param winSize: int. The size that sliding window will have
    '''  
    self.winSize= winSize
    self.centralPos= winSize//2
      
  def compute(self, profileFname, outWindowsFile):
    '''
      Computes sliding windows for one sequence
      @param profileFname: str. path to profiles (processed) results for one sequence
      @param outWindowsFile: str. name of output file
    '''
    prefix, chainType, f_chain, rest = os.path.split(profileFname)[-1].split("_")
    winSizeStep= int(self.winSize/2)

    profileList=[]
    ids=[]
    resNames={}
    noAA_profile=[WindowHHblits.NonAminoacidValueProfile for i in range(WindowHHblits.numProfileElem)] 

    f=open(profileFname)
    header=f.readline().split()
    resIdIndex= header.index("structResId")
    letterIndex= header.index("resName")
    firstValueIndex= header.index("hhblits")
    for line in f:
      lineArray= line.split()
      if len(lineArray)!=0:
        resInd=lineArray[resIdIndex]
        profileList.append([int(val) for val in lineArray[firstValueIndex:(firstValueIndex+WindowHHblits.numProfileElem)]])
        ids.append((f_chain,resInd))
        resNames[ids[-1]]= lineArray[letterIndex]
    f.close()

    numElements= len(ids)
    myProfile={}
    for resNum,resId in enumerate(ids):


      myProfile[resId]=[  profileList[i] if (i>=0 and i<= numElements-1) else noAA_profile
                for i in range(resNum-winSizeStep, resNum+winSizeStep+1)]


    LettersStartPoint= 3
    outFile= open(outWindowsFile,"w")
    for resId in sorted(ids):

      header= ["chainId structResId resName"] + ["hhblits" for i in itertools.chain.from_iterable(myProfile[resId])]

      outFile.write( " ".join(header)+"\n")

      out= list(resId)+ [resNames[resId]]+ list(itertools.chain.from_iterable(myProfile[resId]))
      out= [str(val) for val in out]
      outFile.write( " ".join(out)+"\n")
      break

    for resId in sorted(ids)[1:]:
      out= list(resId)+ [resNames[resId]] + list(itertools.chain.from_iterable(myProfile[resId]))
      out= [str(val) for val in out]
      outFile.write( " ".join(out)+"\n")
    outFile.close()



