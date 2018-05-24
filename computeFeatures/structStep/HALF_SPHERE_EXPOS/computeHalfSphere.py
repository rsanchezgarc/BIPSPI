from __future__ import absolute_import
import os,sys
import numpy as np
from Bio.PDB.Polypeptide import is_aa 
from computeFeatures.structStep.myPDBParser import myPDBParser as PDBParser
#from Bio.PDB.HSExposure import HSExposureCA, HSExposureCB, ExposureCN
from .HSEExposureFixed import HSExposureCA, HSExposureCB, ExposureCN

from ..StructFeatComputer import StructFeatComputer
from utils import myMakeDir, tryToRemove #utils is at the root of the package


class HalfSphereComputer(StructFeatComputer):
  '''
  Extends StructFeatComputer class. Computes half exposure for each residue
  '''
  def __init__(self, rFname, lFname, computedFeatsRootDir=None, statusManager=None):
    '''
        @param statusManager: class that implements .setStatus(msg) to communicate
    '''

    StructFeatComputer.__init__(self, rFname, lFname, computedFeatsRootDir, statusManager=statusManager)

    self.outPath= myMakeDir(self.computedFeatsRootDir, "halfSphereExpos")
    self.parser= PDBParser(QUIET=True)


  def computeOneFile(self, fileName):
    '''
      Computes DSSP for a given pdb file
      @param fileName: str. fname to pdb file
    '''
    
    prefixAndChainTypeId = (fileName.split("/")[-1]).split(".pdb")[0]
    prefixAndChainTypeId= "_".join( prefixAndChainTypeId.split("_")[:2])
    structure=  self.parser.get_structure(prefixAndChainTypeId, fileName)
    model = structure[0]
    outNames= {}
    for chain in structure[0]:
      chainId= chain.get_id()
      nResidues= sum( (1 for res in chain if is_aa(res)) )
      if chainId==" ": chainId="*"
      outName= os.path.join(self.outPath,prefixAndChainTypeId+"_"+chainId+"_u.hse")
      if nResidues>5 and not os.path.isfile(outName):
        outNames[chainId]= outName
    if len(outNames)==0:
      print( "HalfSphere already computed")
      return 0
    
    featuresDict= {}
    hse = HSExposureCA(model)
    for aa, feat in hse:
      featuresDict[aa]= [elem if elem!= None else -1 for elem in feat]
    hseDict={ aa: [elem if elem!= None else -1 for elem in feat] for aa,feat in HSExposureCB(model) }
##    print(len(hseDict))
##    raw_input("press enter")
    for aa in set(featuresDict.keys()).union(set(hseDict.keys())):
      try:
        prevFeatures= featuresDict[aa]
      except KeyError:
        prevFeatures= [-1, -1, -1.0]
      try:
        newFeatures= hseDict[aa]
      except KeyError:
        newFeatures= [-1, -1, -1.0]      
      featuresDict[aa]= prevFeatures + newFeatures
      
    hseDict={ aa: [feat] if feat!= None else -1 for aa,feat in ExposureCN(model) }      
    
    for aa in set(featuresDict.keys()).union(set(hseDict.keys())):
      try:
        prevFeatures= featuresDict[aa]
      except KeyError:
        prevFeatures= [-1, -1, -1.0, -1, -1, -1.0]
      try:
        newFeatures= hseDict[aa]
      except KeyError:
        newFeatures= [-1 ]      
      featuresDict[aa]= prevFeatures + newFeatures
      
    filesHandlers={ chainId:open( outNames[chainId],"w") for chainId in outNames.keys()}
    for fHand in filesHandlers.values():
      fHand.write("chainId structResId resName HSExposureCA1 HSExposureCA2 HSExposureCA3"+
                " HSExposureCB1 HSExposureCB2 HSExposureCB3 ExposureCN\n")

    resuisduesList= [res for chain in structure[0] for res in chain if is_aa(res)]
    badExample= [-1 for elem in list(featuresDict.values())[0]]
    for res in resuisduesList:
#    for res in featuresDict:
#      print(res.get_full_id())
#      print(filesHandlers, res, is_aa(res))
#      raw_input("press enter to continue")
#      
      structId,modelId,chainId, resId= res.get_full_id()
      resId= list(resId)
      resId[1]=str(resId[1])
      resId= "".join(resId[1:])
      resName= self.threeLetterAA_to_one(res.resname)

      if chainId==" ":
        chainId="*"
      try:        
        if chainId not in filesHandlers: continue #small chains
        try:
          valuesForRes= featuresDict[res]
        except KeyError:
          valuesForRes= badExample
        filesHandlers[chainId].write(chainId+" "+resId+" "+resName+" "+" ".join([str(val) for val in valuesForRes])+"\n")
      except (KeyboardInterrupt, Exception):
        for outName in outNames.values():
          print("Exception happend computing %s"%outName)        
          tryToRemove(outName)    
        raise 
    for outFile in filesHandlers.values():
      outFile.close()
 
    return 0

def testModule():
  rFname="/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/b5structures/1GRN_r_u.pdb"
  lFname="/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/b5structures/1GRN_l_u.pdb"
  
#  rFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1A22.pdb"
#  lFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1ABR.pdb"
  
#  computedFeatsRootDir="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/tmpComputedFeatures"
  computedFeatsRootDir="/home/rsanchez/tmp/computedFeats"
  hfsComp= HalfSphereComputer(rFname, lFname, computedFeatsRootDir)
  hfsComp.computeFun()
  raise ValueError("Debug Mode. Comment out testModule() line 128")                
                    

if __name__=="__main__":
  testModule()
  pdbsIndir= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3"
  computedFeatsRootDir= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/computedFeatures"


#  pdbsIndir= "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/b5structures/"
#  computedFeatsRootDir= "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/computedFeatures/structStep/halfSphereExpos"

  StructFeatComputer.computeFeaturesAllComplexes(HalfSphereComputer,pdbsIndir= pdbsIndir ,computedFeatsRootDir= computedFeatsRootDir)


