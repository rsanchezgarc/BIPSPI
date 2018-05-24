'''
This is used for training feature calculations
'''
from __future__ import absolute_import, print_function
import sys, os
from multiprocessing import cpu_count
from Config import Configuration
from .FeaturesComputer import FeaturesComputer
from .common.computeContactMap import ContactMapper
from .common.seqInputPreproceser import SeqInputPreproceser

from .seqStep.getSeqFeatures import SeqFeaturesCalculator

from .structStep.PSAIA.computePSAIA import PSAIAComputer
from .structStep.VORONOI.computeVoronoi import VORONOIComputer
from .structStep.DSSP.computeDssp import DsspComputer
from .structStep.HALF_SPHERE_EXPOS.computeHalfSphere import HalfSphereComputer


#Default parameters
conf= Configuration()
pdbsIndirDefDefault= conf.pdbsIndir
computedFeatsRootDirDefault= conf.computedFeatsRootDir
useCorrMut= conf.useCorrMut

featuresComputers={ "mixed":[("ContactMapper", (ContactMapper, {"boundAvailable": True}) ),
                              ("PSAIAComputer", (PSAIAComputer,{})),
                              ("VORONOIComputer30", (VORONOIComputer, {"maxDist": 30})),
                              ("DsspComputer", (DsspComputer,{})),
                              ("HalfSphereComputer", (HalfSphereComputer,{})),
                              ("SeqFeaturesCalculator", (SeqFeaturesCalculator, {"useCorrMut": useCorrMut}))
                             ],
                    "seq_train":[("ContactMapper", (ContactMapper, {"boundAvailable": True}) ),
                                 ("SeqFeaturesCalculator", (SeqFeaturesCalculator, {"useCorrMut": useCorrMut}))
                          ],
                    "seq_pred":[("SeqInputPreproceser", (SeqInputPreproceser, {}) ),
                                ("SeqFeaturesCalculator", (SeqFeaturesCalculator, {"useCorrMut": useCorrMut}))
                               ]
                   }

def computeFeaturesAllPdbsOneDir(pdbsIndir= pdbsIndirDefDefault, computedFeatsRootDir= computedFeatsRootDirDefault,
                                 methodProtocol="mixed", ncpu=2):
  '''
    Computes all features needed for complex codification for all complexes in pdbsIndir. Used for training
    @param pdbsIndir: str. Path to the directory where pdb files are located. Must be named as follows:
                                pdbId_l_u.pdb and pdbId_r_u.pdb. E.x. 1A2k_l_u.pdb, 1A2k_r_u.pdb. By default,
                                it uses as pdbsIndir Config.py DEFAULT_PARAMETERS["pdbsIndir"] 
    @param computedFeatsRootDir: str. Path where features files will be saved. By default it uses
                                Config.py DEFAULT_PARAMETERS["computedFeatsRootDir"] will be used as out_path
    @param methodProtocol: str. "mixed" if structural and sequential features will be used. "seq" if just 
                                sequential features used
    @param ncpu: int. Number of cpu's to use. If -1, all cpu's will be used
  '''
  assert methodProtocol=="mixed" or methodProtocol=="seq_train", "Error methodProtocol in computeFeaturesAllPdbsOneDir "+\
                                                              "must be 'mixed' or 'seq_train'"
  if ncpu<1:
    ncpu= cpu_count()
  for name, (featCompClass, classArgs) in featuresComputers[methodProtocol]:
    print ("Applying %s"% str(featCompClass))
    FeaturesComputer.computeFeaturesAllComplexes( featCompClass,pdbsIndir= pdbsIndir,
                      computedFeatsRootDir= computedFeatsRootDir, classArgs= classArgs, ncpu= ncpu)
    print(" %s Done"%str(featCompClass))

def computeFeaturesOneComplex(lFname, rFname, lPdbId=None, rPdbId=None,computedFeatsRootDir= computedFeatsRootDirDefault, 
                              methodProtocol="mixed", areForTrainAndTest=False, statusManager= None):
  '''
    Computes all features needed for complex codification.
    @param lFname: str. Path to the the pdb file or the fasta file of the ligand of the complex
    @param rFname: str. Path to the the pdb file or the fasta file of the receptor of the complex
    @param rPdbId: str. pdbId of receptor in order to query 3dCons. If None, psiblast will be launched directly
    @param lPdbId: str. pdbId of receptor in order to query 3dCons. If None, psiblast will be launched directly
    @param computedFeatsRootDir: str. Path where features files will be saved. By default it uses
                                Config.py DEFAULT_PARAMETERS["computedFeatsRootDir"] will be used as out_path
    @param methodProtocol: str. "mixed" if structural and sequential features will be used. "seq_pred" if just 
                                sequential features used
    @param areForTrainAndTest: boolean. True if ligand and receptor are in interacting coordinates and thus,
                                contact maps are needed to be computed in order to perform evaluation.
                                False if ligand and receptor coordinates are not related and thus, 
                                evaluation does not makes sense.
    @param statusManager: class that implements .setStatus(msg) to communicate
  '''
  assert methodProtocol=="mixed" or methodProtocol=="seq_pred", "Error methodProtocol in computeFeaturesAllPdbsOneDir "+\
                                                    "must be 'mixed' or 'seq_pred'"
  for name, (featCompClass, classArgs) in featuresComputers[methodProtocol]:
    print ("Applying %s"% str(featCompClass))
    if name== "ContactMapper": 
      classArgs= classArgs.copy()
      classArgs.update({"boundAvailable": areForTrainAndTest, "isForPrediction": not areForTrainAndTest})
    elif name== "SeqFeaturesCalculator": 
      if not areForTrainAndTest: #if they are not for benchmark train and test but for de novo evaluation
        classArgs.update({"rPdbId": rPdbId, "lPdbId": lPdbId})
    classArgs.update({"statusManager": statusManager})
    featCompClass(rFname, lFname, computedFeatsRootDir, **classArgs).computeFun()
    print(" %s Done"%str(featCompClass))

def testModulePredictMixed():

#  rFname="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/1A2K_r_u.pdb"
#  lFname="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/1A2K_l_u.pdb"
  lFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/4YGA.pdb"
  rFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1ABR.pdb"
  
  computedFeatsRootDir="/home/rsanchez/tmp/computedFeats"
  computeFeaturesOneComplex( lFname,rFname,  computedFeatsRootDir)
  raise ValueError("Debug Mode. Comment out testModule() line 102")


def testModulePredictSeq():
  
  lFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/inputFasta/seq1_l_u.fasta"
  rFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/inputFasta/seq1_r_u.fasta"
  
  computedFeatsRootDir="/home/rsanchez/tmp/computedFeats"
  computeFeaturesOneComplex(lFname, rFname, computedFeatsRootDir, methodProtocol="seq_pred")  

  raise ValueError("Debug Mode. Comment out testModule() line 105")
  
if __name__ == "__main__":
#  testModulePredictMixed(); print("test done"); sys.exit(1)
#  testModulePredictSeq(); print("test done"); sys.exit(1)
  
  if len(sys.argv)>1:
    pdbsIndir= sys.argv[1]
  if len(sys.argv)>2:
    computedFeatsRootDir= sys.argv[2]
  if len(sys.argv)>3:
    ncpu= int(sys.argv[3])
    
  computeFeaturesAllPdbsOneDir(pdbsIndir, computedFeatsRootDir, ncpu)
  

