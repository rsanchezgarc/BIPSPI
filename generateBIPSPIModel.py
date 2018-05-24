'''
This scripts is used to generate models.

Models are trained from pdb files which are contained in the directory pdbsIndir
that has to be indicated in ./configFiles/configFile.cfg, (p.e. pdbsIndir ~/Tesis/rriPredMethod/data/develData/pdbFiles)
  For each pdb in dataset, there must be 4 files included in pdbsIndir directory, and they have to follow the following
  name schema:
    prefix_l_u.pdb File for the ligand in unbound conformation. Features are computed from this file
    prefix_r_u.pdb File for the receptor in unbound conformation. Features are computed from this file
    prefix_l_b.pdb File for the ligand in bound conformation. Residue contacts are computed from this file
    prefix_r_b.pdb File for the receptor in bound conformation. Residue contacts are computed from this file
  prefix are generally PDBIDs, but can be any string provided it is shared by ligand and receptor.
  
Features computed are hierachically stored in the computedFeatsRootDir directory that has to be indicated in 
./configFiles/configFile.cfg (e.g. computedFeatsRootDir ~/Tesis/rriPredMethod/data/develData/computedFeatures)

Ready to train and predict complexes are stored in the codifiedDataRootDir directory that has to be indicated in 
./configFiles/configFile.cfg (e.g. codifiedDataRootDir  ~/Tesis/rriPredMethod/data/develData/codifiedInput).
Complexes are saved as joblib.dump files of the Class ./codifyComplex/ComplexCodified. See that class for more
info.

Results obtained in cross validation are stored in the resultsRootDir directory that has to be indicated in 
./configFiles/configFile.cfg (e.g. resultsRootDir  ~/Tesis/rriPredMethod/data/develData/results).
Three files per complex are generated:
  prefix.res.tab      Residue-residue contact predictions
  prefix.res.tab.lig  binding site predictions for ligand
  prefix.res.tab.lig  binding site predictions for receptor
  

The models trained are strored in the savedModelsPath that has to be indicated in 
./configFiles/configFile.cfg (e.g. savedModelsPath ~/Tesis/rriPredMethod/data/develData/modelsComputed).
Models are saved as joblib.dump files of the Class xgboost.XGBClassifier

To run this script, edit config file
load conda environment
  source activate xgbpred
and execute
  python generateMixedModel

'''
from __future__ import print_function
import os, sys
from psutil import virtual_memory

import computeFeatures.computeFeatsForPdbs as pComCode
import codifyComplexes.codifyPDBsForTraining as pCodifyAll
import trainAndTest.trainAndTest as pTrainTest
from utils import myMakeDir
from Config import Configuration

memMax_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') 
memMax_gib = memMax_bytes/(1024.**3)

conf= Configuration()

def computeFeatures(methodProtocol= conf.modelType):
  pComCode.computeFeaturesAllPdbsOneDir(ncpu= 1+ conf.ncpu//conf.psiBlastNThrs, methodProtocol= methodProtocol)

def codifyStep(stepType= conf.modelType, feedbackPaths=None):
  benchCod= pCodifyAll.BenchmarkCodificator( feedback_paths= feedbackPaths, environType=stepType, 
                                             ncpu=min( int(1+ memMax_gib//8), conf.ncpu), overridePrevComp= False)
  cMapsPath= os.path.join(os.path.expanduser(benchCod.dataRootPath),"common","contactMaps")
  prefixes= [ elem.split(".")[0] for elem in os.listdir(cMapsPath)]
  prefixes.sort()
#  skipComplexesList= prefixes[4:]
  skipComplexesList= []
  codifiedPath= benchCod.codifyAll( skipComplexesList=skipComplexesList, samplingFold= 3)
  return codifiedPath
  
def trainAndTest(inputRoot, outputRoot, stepType= conf.modelType, saveModelPath= None):
  trainAndTester= pTrainTest.TrainAndTest()
  sampledComplexesPath= os.path.join(inputRoot, "sampledInputs")
  wholeComplexesPath= os.path.join(inputRoot, "allInputs")
  predictOutputPath= os.path.join(outputRoot, stepType)
  myMakeDir(predictOutputPath)
  
  numResults= sum([1 for elem in os.listdir(predictOutputPath)])
  if numResults!=0:
    print("Warning: predictions outpath is not empty %s"%predictOutputPath)
  trainAndTester.setParameters(sampledComplexesPath, wholeComplexesPath, predictOutputPath, numProc=conf.ncpu,
                      nFolds=conf.N_KFOLD, verbose=True, saveModelFname=saveModelPath)
  trainAndTester.main()                    
  return predictOutputPath
  
def main(resultsRoot= None, saveModelName=None):
  if resultsRoot is None:
    resultsRoot= os.path.expanduser(conf.resultsRootDir)
  computeFeatures()
  outpathCodif= codifyStep()
  saveModelNameTmp= None if saveModelName is None else saveModelName+"."+conf.modelType
  predictOutPath_1= trainAndTest( outpathCodif, resultsRoot, stepType= conf.modelType, saveModelPath= saveModelNameTmp)
  if conf.modelType in ["struct", "mixed"]:
    outpathCodif= codifyStep(stepType=conf.modelType+"_2", feedbackPaths= predictOutPath_1)
    saveModelNameTmp= None if saveModelName is None else saveModelName+"."+conf.modelType+"_2"  
    predictOutPath_2= trainAndTest( outpathCodif, resultsRoot, stepType= conf.modelType+"_2", 
                                    saveModelPath= saveModelNameTmp)
if __name__=="__main__":
  main(saveModelName=os.path.join(conf.savedModelsPath,"model"))


