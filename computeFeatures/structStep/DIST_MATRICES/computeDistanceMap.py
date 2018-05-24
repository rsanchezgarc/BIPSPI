from __future__ import absolute_import
import os,sys
import numpy as np
from Bio.PDB.Polypeptide import is_aa 
from computeFeatures.structStep.myPDBParser import myPDBParser as PDBParser

from ..StructFeatComputer import StructFeatComputer
from utils import myMakeDir, tryToRemove #utils is at the root of the package

class DistMapComputer(StructFeatComputer):
  '''
  Extends StructFeatComputer class. Computes pairwise distance between all amino acids of receptor and ligand independently
  '''
  def __init__(self, rFname, lFname, computedFeatsRootDir=None, statusManager=None):
    '''
      @param rFname: str. path to receptor pdb file
      @param lFname: str. path to ligand pdb file      
      @param computedFeatsRootDir: str. path where features will be stored
      @param statusManager: class that implements .setStatus(msg) to communicate
    '''
    StructFeatComputer.__init__(self, rFname, lFname, computedFeatsRootDir, statusManager=statusManager)

    self.outPath= myMakeDir(self.computedFeatsRootDir, "distanceMatricesData")
    self.parser= PDBParser(QUIET=True)


  def computeOneFile(self, fileName):
    '''
      Computes distance for each pair of aminoacids for a given pdb file
      @param fileName: str. fname to pdb file
    '''
    prefixAndChainTypeId = (fileName.split("/")[-1]).split(".pdb")[0]
    outName= os.path.join(self.outPath,prefixAndChainTypeId+".distMat")
    if os.path.isfile(outName): 
      print("Already computed Distance Maps")
      return 0
    structure=  self.parser.get_structure(prefixAndChainTypeId, fileName)
    structCenterMass= self.getStructCenterMass(structure)

    try:
      outFile= open(outName,"w")
      outFile.write("chainId1 structResId1 chainId2 structResId2 distance angle_to_protCM\n")
      for res1 in structure[0].get_residues():
        if is_aa(res1,standard=True):
    ##        print res, res.get_full_id()
          structId1,modelId1,chainId1, resId1= res1.get_full_id()
          resId1= list(resId1)
          resId1[1]=str(resId1[1])
          resId1= "".join(resId1[1:])
          if chainId1==" ":
            chainId1="*"
          for res2 in structure[0].get_residues():
            if is_aa(res2,standard=True):
      ##        print( res, res.get_full_id())
              structId2,modelId2,chainId2, resId2= res2.get_full_id()
              resId2= list(resId2)
              resId2[1]=str(resId2[1])
              resId2= "".join(resId2[1:])
              if chainId2==" ":
                chainId2="*"
              magnitude= self.getMagnitude(res1, res2, structCenterMass)
#              print( chainId1, resId1, chainId2, resId2, magnitude)
#              a= raw_input()
              outFile.write(chainId1+" "+resId1+" "+chainId2+" "+resId2+" "+" ".join([str(val) for val in magnitude])+"\n")
      outFile.close()
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing %s"%outName)
      tryToRemove(outName)    
      raise
    return 0

  def getStructCenterMass(self, struct):
    '''
      Computes center of mass for a given structure
      @param pdbStruct: Bio.PDB.Structure. Structure of the pdb that is being analyzed
      @return np.array len=3
    '''  
    resList= struct[0].get_residues()
    coords= np.zeros((len(list(resList)),3))
    for i,res in enumerate(resList):
      try:
        coords[i,:]= res["CA"].get_coord()
      except KeyError:
        try:
          coords[i,:]= res["CB"].get_coord()
        except KeyError:
          try:
            coords[i,:]= res["N"].get_coord()
          except KeyError:
            try:
              coords[i,:]= res["O"].get_coord()
            except KeyError:
              coords[i,:]= res.get_list()[0].get_coord()
    return np.mean(coords,axis=0)


  def computeAngle(self, res1_atom,  res2_atom, structCenterMass):
    '''
      Computes angle between 2 atoms with respect one point
      @param res1_atom: Bio.PDB.Atom. Atom
      @param res2_atom: Bio.PDB.Atom. Atom     
      @param structCenterMass: np.array len=3. Point where vector origin will be set
      @return float. The angle
    '''    
    vect1= res1_atom.get_coord() - structCenterMass
    vect2= res2_atom.get_coord() - structCenterMass
    vect1 = vect1 /np.linalg.norm(vect1)
    vect2 = vect2 /np.linalg.norm(vect2)
    return np.arccos(np.clip(np.dot(vect1, vect2), -1.0, 1.0))

  def getMagnitude(self, res1, res2, structCenterMass):
    '''
      Computes both distance and angle with respect protein center of mass of 2 residues
      @param res1: Bio.PDB.Residue. res1
      @param res2: Bio.PDB.Residue. res2     
      @param structCenterMass: np.array len=3. Point where vector origin will be set
      @return Tuple (float, float). (Distance, angle)
    '''      
    try:
       res1_atom= res1["CA"]
    except KeyError:
      try:
         res1_atom= res1["CB"]
      except KeyError:
        try:
           res1_atom=  res1["CA"]
        except KeyError:
           res1_atom=  res1.get_list()[0]
    try:
       res2_atom= res2["CA"]
    except KeyError:
      try:
         res2_atom= res2["CB"]
      except KeyError:
        try:
           res2_atom=  res2["CA"]
        except KeyError:
           res2_atom=  res2.get_list()[0]

    resDist= round( res1_atom - res2_atom,4)
    resAngl= round( self.computeAngle( res1_atom,  res2_atom, structCenterMass), 4)
    return resDist,resAngl

def testModule():
#  rFname="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/1A2K_r_u.pdb"
#  lFname="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/1A2K_l_u.pdb"
  
  rFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1A22.pdb"
  lFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1ABR.pdb"
  
  computedFeatsRootDir="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/tmpComputedFeatures"
  distComp= DistMapComputer(rFname, lFname, computedFeatsRootDir)
  distComp.computeFun()
  raise FeatureComputerException("Debug Mode. Comment out testModule() line 128")                
                    

if __name__=="__main__":

#  testModule()

  pdbsIndir= "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/b5structures/"
  computedFeatsRootDir= "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/newCodeData/computedFeatures/structStep/DIST_MAT"

  StructFeatComputer.computeFeaturesAllComplexes(DistMapComputer,pdbsIndir= pdbsIndir ,computedFeatsRootDir= computedFeatsRootDir)

