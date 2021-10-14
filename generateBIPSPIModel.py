'''
This scripts is used to generate BIPSPI models.

Models are trained from pdb files which are contained in the directory pdbsIndir
that has to be indicated in ./configFiles/configFile.cfg, (p.e. pdbsIndir ~/Tesis/rriPredMethod/data/develData/pdbFiles)
  For each pdb in dataset, there must be 4 files included in pdbsIndir directory, and they have to follow the following
  name standard:
    prefix_l_u.pdb File for the ligand in unbound conformation. Features are computed from this file
    prefix_r_u.pdb File for the receptor in unbound conformation. Features are computed from this file
    prefix_l_b.pdb File for the ligand in bound conformation. Residue contacts are computed from this file
    prefix_r_b.pdb File for the receptor in bound conformation. Residue contacts are computed from this file
  prefix are generally PDBIDs, but can be any string provided it is shared by ligand and receptor.
  if no unbound pdbs available, use symlinks  prefix_r_b.pdb -->  prefix_r_u.pdb; prefix_l_b.pdb -->  prefix_l_u.pdb
  
Computed features are hierachically stored in the computedFeatsRootDir directory that has to be indicated in 
./configFiles/configFile.cfg (e.g. computedFeatsRootDir ~/Tesis/rriPredMethod/data/develData/computedFeatures)

Ready to train and predict complexes are stored in the codifiedDataRootDir directory that has to be indicated in 
./configFiles/configFile.cfg (e.g. codifiedDataRootDir  ~/Tesis/rriPredMethod/data/develData/codifiedInput).
Complexes are saved as joblib.dump files of the class ./codifyComplex/ComplexCodified. See that class for more
info.

Results obtained in cross validation are stored in the resultsRootDir directory that has to be indicated in 
./configFiles/configFile.cfg (e.g. resultsRootDir  ~/Tesis/rriPredMethod/data/develData/results).
Three files per complex are generated:
  prefix.res.tab.gz      Residue-residue contact predictions
  prefix.res.tab.lig.gz  binding site predictions for ligand
  prefix.res.tab.lig.gz  binding site predictions for receptor
  

The models trained are strored in the savedModelsPath that has to be indicated in 
./configFiles/configFile.cfg (e.g. savedModelsPath ~/Tesis/rriPredMethod/data/develData/modelsComputed).
Models are saved as joblib.dump files of the Class xgboost.XGBClassifier

To run this script, edit config file
load conda environment
  source activate xgbpred
and execute
  python generateBIPSPIModel

'''
from __future__ import print_function
import os

import computeFeatures.computeFeatsForPdbs as pComCode
import codifyComplexes.codifyPDBsForTraining as pCodifyAll
import trainAndTest.trainAndTest as pTrainTest
from utils import myMakeDir, checkFreeMemory
from Config import Configuration


GiB_PER_PROC=16

conf= Configuration()

USE_2_STEPS= True

def computeFeatures(methodProtocol= conf.modelType, isHomeSet= conf.checkHomoInteractionInTraining):
  pComCode.computeFeaturesAllPdbsOneDir(ncpu= 1+ conf.ncpu//conf.psiBlastNThrs, methodProtocol= methodProtocol, isHomeSet=isHomeSet)

def codifyStep(methodProtocol= conf.modelType, feedbackPaths=None):
  benchCod= pCodifyAll.BenchmarkCodificator( feedback_paths= feedbackPaths, environType=methodProtocol, 
                                             ncpu=min( int(1+ checkFreeMemory()//GiB_PER_PROC), conf.ncpu), overridePrevComp= False)
  if not benchCod.checkIfAllCodified():
  #  skipComplexesList= benchCod.prefixes[4:]
    skipComplexesList= []
    codifiedPath= benchCod.codifyAll( skipComplexesList=skipComplexesList, samplingFold= 3)
  else:
    codifiedPath= benchCod.getCodifiedPath()
    print("All complexes already codified")
  return codifiedPath
  
def trainAndTest(inputRoot, outputRoot, methodProtocol= conf.modelType, saveModelPath= None, isLastStep=False):

  sampledComplexesPath= os.path.join(inputRoot, "sampledInputs")
  wholeComplexesPath= os.path.join(inputRoot, "allInputs")
  predictOutputPath= os.path.join(outputRoot, methodProtocol)
  myMakeDir(predictOutputPath)
  
  numResults= len(os.listdir(predictOutputPath))
  if numResults!=0:
    print("Warning: predictions outpath is not empty %s"%predictOutputPath)
  
  trainAndTester= pTrainTest.TrainAndTestWorker(trainDataPath=sampledComplexesPath, testPath=wholeComplexesPath, 
                              outputPath=predictOutputPath, nFolds=conf.N_KFOLD, isLastStep=isLastStep,
                              saveModelFname=saveModelPath, verbose=True, numProc=conf.ncpu)
  trainAndTester.computeTrainAndTest()
  return predictOutputPath
  
def main(resultsRoot= None, saveModelName=None):
  if resultsRoot is None:
    resultsRoot= os.path.expanduser(conf.resultsRootDir)
  computeFeatures()
  outpathCodif= codifyStep()
  myMakeDir(saveModelName)
  saveModelNameTmp= None if saveModelName is None else saveModelName+"."+conf.modelType
  isLastStep= True
  if USE_2_STEPS:
    isLastStep=False
  predictOutPath_1= trainAndTest( outpathCodif, resultsRoot, methodProtocol= conf.modelType, saveModelPath= saveModelNameTmp, isLastStep=isLastStep)
  if USE_2_STEPS:
    outpathCodif= codifyStep(methodProtocol=conf.modelType+"_2", feedbackPaths= predictOutPath_1)
    saveModelNameTmp= None if saveModelName is None else saveModelName+"."+conf.modelType+"_2"

    predictOutPath_2= trainAndTest( outpathCodif, resultsRoot, methodProtocol= conf.modelType+"_2", 
                                    saveModelPath= saveModelNameTmp, isLastStep=True)

if __name__=="__main__":
  parser = Configuration.getArgParser()
  parser.modify_field("pdbsIndir", help="Directory where training pdbs are located", _type= Configuration.file_path)
  parser.modify_field("wdir", help="Directory where partial results and final results will be saved", _type= Configuration.file_path)
  parser.modify_field("tmp", help="Temporary directory", _type= Configuration.file_path)

  parser.modify_field("ncpu", help="Number of cpus for trainng. Each complex in a cross-validation fold is computed in an indepented worker. NCPU workers are computed in parallel")
  parser.modify_field("modelType", help="The type of model to train depending on input options. Struct,sequence or mixed (one seq and one struct)", choices=["struct", "mixed", "seq"])
  parser.modify_field("checkHomoInteractionInTraining", help="Corrects contact maps to consider all homologous residues as contact if one does. For heter-datasets impact is minorSet it to 'True' "
                                                             "if training homo-complexes. For hetero complexes, it has little impact", _type=bool)
  parser.modify_field("N_KFOLD", help="Type of cross validation. -1 for leave-one-complex out, positive values for k= N_KFOLD cross-validation.", _type=Configuration.int_or_filePath)
  parser.modify_field("scopeFamiliesFname", help="Filename containing the familes of the protein chains of ligand and receptor", _type= Configuration.file_path)

  parser.parse_args()

  main(saveModelName=os.path.join(conf.savedModelsPath,"model"))


