'''
This module allows for the prediction of Residue-Residue Contacts and Binding Sites given 2 protein partners
'''
from __future__ import print_function

import argparse
import gzip
import os
import re
import shutil
import sys
import traceback

import codifyComplexes.codifyOneComplex as pCodifyOne
from bigExceptions import NoAvailableForDownloadPDB, NoValidPDBFile
from computeFeatures.computeFeatsOneComplex import computeFeaturesOneComplex
from computeFeatures.structStep.structFeatsComputer import moveAndWriteAsPDBIfMmcif
from computeFeatures.toolManagerGeneric import checkIfIsFasta, parseFasta, parseSeq
from joblib import Parallel, delayed
from pythonTools.downloadPdb import downloadPDB
from trainAndTest.predictOneCodifiedComplex import ComplexPredictor
from utils import myMakeDir

DEFAULT_PREFIX="PRED"
COMPUTE_FEATS_ONLY=False

def prepareInput(allRootDir, idName, pdbId, fname, sequence, chainType="l", removeInputs=False):
  '''
  :param allRootDir: str. Path where files will be generated
  :param idName: an id for the complex, typically a pdbId
  :param pdbId: A pdbId, potentially augmented with chainId and bioUnit info.
  :param fname: str. A path to a fasta or pdb file. Ignored if pdbId provided
  :param sequence: str. A  sequence of amino acids. Ignored if pdbId or fname provided
  :param chainType: "l" or "r"
  :param removeInputs: bool. Remove input files
  '''
  assert chainType in ["l","r"], "Error, chain type must be either 'l' or 'r'"
  finalInputPath= myMakeDir(allRootDir, "finalInput")

  newFname= os.path.join(finalInputPath,"%s_%s_u.pdb"%(idName, chainType))

  if pdbId and fname is None  and sequence is None:
    print("%s_%s"%(pdbId, chainType))
    matchObj= re.match(r"(^\d[a-zA-Z0-9]{3})(:[a-zA-Z0-9])?(_\d*)?$", pdbId)
    if matchObj:
      pdbId= matchObj.group(1)
      chainId= matchObj.group(2)
      biounit= matchObj.group(3)
      isThisParnerSeq=False
      if not chainId is None: chainId= chainId[1:]
      if biounit is None or len(biounit)==1: biounit=None
      else:
        biounit= int(biounit[1:])
        if biounit==0: biounit=None
      try:
        fname = downloadPDB( pdbId, finalInputPath, chainId= chainId, bioUnit=biounit, removeOtherChains=removeInputs)
        print(type(fname), fname)
      except ValueError as e:
        traceback.print_exc()
        raise NoAvailableForDownloadPDB( str(e)+"-> %s"%(pdbId))
      try:
        if not removeInputs:
          shutil.copyfile(fname, newFname)
        else:
          os.rename(fname, newFname)
      except OSError:
        NoAvailableForDownloadPDB("Error. pdb %s was not recover from wwpdb.org"%pdbId)
    else:
      raise NoAvailableForDownloadPDB("Error: bad pdb provided %s"%(pdbId))
  elif fname:
    if checkIfIsFasta(fname):
      sequence= parseFasta(fname)
      if removeInputs:
        os.remove(fname)
      isThisParnerSeq=True
    else: #check if fname is valid pdb
      success= moveAndWriteAsPDBIfMmcif(fname, newFname, removeInput= removeInputs)
      if not success:
        raise NoValidPDBFile("Error. It was not possible to parse your pdb file %d"%( 1 if chainType=="l" else 2))
      isThisParnerSeq=False
  elif sequence: 
    sequence= parseSeq(sequence)
    isThisParnerSeq=True
  else:    
    raise ValueError("One of the following arguments should be provided pdbId, fname, sequence")
  if sequence: #Save sequence
    newFname= os.path.join(finalInputPath,"%s_%s_u.fasta"%(idName, chainType))
    with open(newFname,"w") as f:
      f.write(">"+ os.path.basename(newFname)+"\n"+sequence)
  return isThisParnerSeq, newFname, pdbId
  
def computeFeatures(allRootDir, idName, lFname, rFname, lPdbId=None, rPdbId=None, methodProtocol="struct", 
                    areForTrainAndTest=False, isHomoComplex=False, statusManager=None):
  '''

  :param allRootDir: str. A path where all files will be saved
  :param idName: An id that will be used for files prefix
  :param lFname: A path to a fasta or pdb file
  :param rFname: A path to a fasta or pdb file
  :param lPdbId: A pdbId (4 letters) or None. Used to query 3D-cons
  :param rPdbId: A pdbId (4 letters) or None. Used to query 3D-cons
  :param methodProtocol: str. "struct" "mixed" or "seq_pred"
  :param areForTrainAndTest: boolean. True if ligand and receptor are in interacting coordinates and thus,
                            contact maps are needed to be computed in order to perform evaluation.
                            False if ligand and receptor coordinates are not related and thus,
                            evaluation does not makes sense.
  :param isHomoComplex: : boolean. If True, just lPartner is provided and both partners are equal. It will use model.[inputType][_2].homo
  :param statusManager: class that implements .setStatus(msg) to communicate
  :return:
  '''
  
  computedFeatsRootDir= myMakeDir(allRootDir, "computedFeatures")
  computeFeaturesOneComplex( lFname=lFname, rFname=rFname, prefix=idName, lPdbId=lPdbId, rPdbId= rPdbId,
                            computedFeatsRootDir= computedFeatsRootDir, methodProtocol=methodProtocol,
                            areForTrainAndTest=areForTrainAndTest, isHomoComplex=isHomoComplex, statusManager= statusManager)
  return computedFeatsRootDir

def codifyStep(stepType, computedFeatsRootDir, feedbackPaths=None, idName=None):

    
#  wholeComplexOutPath= os.path.join(allRootDir, "codifiedInput")
  wholeComplexOutPath= None  #We don't want the ready to used data to be saved
  oneComplexCod= pCodifyOne.OneComplexCodifier( dataRootPath= computedFeatsRootDir, sampledOutPath=None,
                                                wholeComplexOutPath=wholeComplexOutPath,
                                                environType=stepType, feedback_paths= feedbackPaths)

  codifiedObject= oneComplexCod.codifyComplex( idName)
  print(codifiedObject, codifiedObject.shape)
  return codifiedObject
  
def predict(allRootDir, stepType, complexCodifiedObject, isHomoComplex, areForTrainAndTest=False, isLastStep=False):
  predictOutputPath= myMakeDir(allRootDir, stepType)
  predictor= ComplexPredictor(stepType, isHomoComplex, averageLRscores= isHomoComplex and not areForTrainAndTest) # If it is a homo_complex but not benchmark, then L and R are the same
  isMixedProto= True if stepType.startswith("mixed") else False
  outPath, resultsDfs= predictor.predictOneComplex(complexCodifiedObject, 
                                            os.path.join(predictOutputPath,"%s.res.tab"%complexCodifiedObject.prefix),
                                            isMixedProto=isMixedProto, isLastStep=isLastStep)
  return predictOutputPath, resultsDfs
  

def cleanDirectory(dirName):
  print("cleaning %s"%dirName )
  for root, dirs, files in os.walk(dirName, topdown=False):
    for name in files:
      fname= os.path.join(root, name)
      if not fname.endswith(".gz"):
        print(fname)
        with open(fname, 'rb') as f_in:
          with open(fname+'.gz', 'wb') as f_out:
            with gzip.GzipFile(fname, 'wb', fileobj=f_out) as f_out2:
                shutil.copyfileobj(f_in, f_out2)
        os.remove(fname)
    for name in dirs:
      dirName= os.path.join(root, name)
      if name.startswith("computedFeatures"):
        shutil.rmtree(dirName)

def predictOneComplex(allRootDir, lPdbId=None, rPdbId=None, lFname=None, rFname=None,
           lSequence=None, rSequence=None, idName=DEFAULT_PREFIX, cleanForDisplay=False,
           areForTrainAndTest=False, removeInputs=False, isHomoComplex=False, statusManager=None):
  '''
    obtains Residue-Residue Contact predictions and binding site predictions for the sequences or structures contained in
            the files lFname, rFname  or indicated in lpartnerStr, rpartnerStr (as a sequence string or as a pdbId)

    :param allRootDir: str. A path where all files will be saved
    :param lPdbId: str. A pdbId (4 letters) for partner l. Can contain info about chain and/or biounit: pdb[:chain][_bioUnit]
    :param lFname: str. A path to a fasta or pdb file for partner r. Ignored if lPdbId provided
    :param lSequence: str. A  sequnce of amino acids for partner l. Ignored if lPdbId or lFname provided

    :param rPdbId: str. A pdbId (4 letters) for partner r. Can contain info about chain and/or biounit: pdb[:chain][_bioUnit]
    :param rFname: str. A path to a fasta or pdb file for partner r. Ignored if rPdbId provided
    :param rSequence: str. A  sequnce of amino acids for partner r. Ignored if rPdbId or rFname provided

    :param idName: str. An id for the complex that will be used for files prefix
    :param cleanForDisplay: Boolean. If true, categ will be removed from output Dataframes and l will be substituted by
                                      partner 1 and r by parnter 2.
    :param areForTrainAndTest: boolean. True if ligand and receptor are in interacting coordinates and thus,
                            contact maps are needed to be computed in order to perform evaluation.
                            False if ligand and receptor coordinates are not related and thus,
                            evaluation does not makes sense.
    :param removeInputs: boolen. If True, inputs will be removed
    :param isHomoComplex: boolean. If True, just lPartner is provided and both partners are equal. It will use model.[inputType][_2].homo
    :param statusManager: class that implements .setStatus(msg) to communicate status
    :return predictOutPath, resultsDfs
               predictOutPath: path were results were saved
               resultsDfs: (resultsPairs, resultsL, resultsR).
  '''
  assert lFname is not None or lPdbId is not None or lSequence is not None, "Error, input not provided for partner l"
  if not isHomoComplex:
    assert rFname is not None or rPdbId is not None or rSequence is not None, "Error, input not provided for partner r"

  if not lFname is None: lFname= os.path.expanduser(lFname)
  if not rFname is None: rFname= os.path.expanduser(rFname)

  if areForTrainAndTest:
    if not lFname is None:
      prefixL= os.path.basename(lFname).split(".")[0].split("_")[0]
      lPdbId= prefixL
    else:
      prefixL= lPdbId.split(":")[0]
    if not rFname is None:
      prefixR= os.path.basename(rFname).split(".")[0].split("_")[0]
      rPdbId= prefixR
    else:
      prefixR= rPdbId.split(":")[0]
    assert prefixL==prefixR, "Error, different prefixes when areForTrainAndTest=True. %s %s"%(prefixL, prefixR)

  isLigandSeq, lFname, lPdbId= prepareInput(allRootDir, idName, lPdbId, lFname, lSequence, chainType="l",
                            removeInputs=removeInputs and not areForTrainAndTest)

  if isHomoComplex:
    isReceptSeq, rFname, rPdbId = isLigandSeq, lFname, lPdbId
  else:
    isReceptSeq, rFname, rPdbId= prepareInput(allRootDir, idName, rPdbId, rFname, rSequence, chainType="r",
                              removeInputs=removeInputs and not areForTrainAndTest)
  print("Input is sequence? l: %s r: %s"%(isLigandSeq, isReceptSeq) )

  numTotalSteps= 5
  if isLigandSeq and isReceptSeq:
    protocolType="seq"
  elif not( isLigandSeq or isReceptSeq):
    if hasattr(list, "force_protocol") and conf.force_protocol is not None:
      protocolType = conf.force_protocol
    else:
      protocolType="struct"
  else:
    protocolType="mixed"

  if not statusManager is None: statusManager.setStatus("STEP 1/%d: Computing features"%(numTotalSteps))
  print("Protocol %s"%(protocolType))
  computedFeatsRootDir= computeFeatures(allRootDir, idName, lFname=lFname, lPdbId=lPdbId, rFname=rFname, rPdbId=rPdbId,
                                        methodProtocol=protocolType, areForTrainAndTest=areForTrainAndTest,
                                        isHomoComplex= isHomoComplex, statusManager= statusManager)
  if COMPUTE_FEATS_ONLY:
    print("MODE: COMPUTE_FEATS_ONLY.\nDONE")
    return
  predsRootDir= myMakeDir(allRootDir, "preds")
  if protocolType.startswith("seq") or protocolType.startswith("mixed") or protocolType.startswith("struct"):
    if not statusManager is None: statusManager.setStatus("STEP 2/%d: Data encoding"%(numTotalSteps))
    codifiedObject= codifyStep(protocolType, computedFeatsRootDir, idName=idName)

    if not statusManager is None: statusManager.setStatus("STEP 3/%d: Predictions calculation"%(numTotalSteps))
    predictOutPathMixed, __= predict(predsRootDir, protocolType, codifiedObject, isHomoComplex, areForTrainAndTest)

    if not statusManager is None: statusManager.setStatus("STEP 4/%d: Feedback data encoding"%(numTotalSteps))
    codifiedObject= codifyStep(protocolType + "_2", computedFeatsRootDir, [predictOutPathMixed], idName=idName)

    if not statusManager is None: statusManager.setStatus("STEP 5/%d: Feedback predictions calculation"%(numTotalSteps))
    predictOutPath, resultsDfs= predict(predsRootDir, protocolType+"_2", codifiedObject, isHomoComplex, areForTrainAndTest,
                                        isLastStep=cleanForDisplay)
  else:
    raise ValueError("protocolType must be mixed or seq or struct")

  if cleanForDisplay==True:
    for df in resultsDfs:
      if "categ" in df.columns:
        df= df.drop("categ",axis=1)
      for colname in df.columns:
        if colname.startswith("chainId"):
          df[colname]= df[colname].map(lambda x: x if  x!="*" else " ")
    cleanDirectory(allRootDir)
  return predictOutPath, resultsDfs


def predictAllComplexesInDir(inputDir, outDir, isHomo=False, areForTrainAndTest=True, ncpu=1):
  prefixes2fname= {}
  for fname in sorted(os.listdir(inputDir)):
    fname_split= fname.split(".")
    prefix, chainType, boundState= fname_split[0].split("_")
    if not prefix in prefixes2fname:
      prefixes2fname[prefix]=[None, None]
    if boundState=="u":      
      prefixes2fname[prefix][0 if chainType=="l" else 1]= os.path.join(inputDir, fname)
  oneFname= list(prefixes2fname.values())[0][0]
  assert oneFname.endswith(".pdb") or oneFname.endswith(".fasta"), "Error, files in directory must be .pdb or .fasta"

  argsList=[]
  for prefix in sorted(prefixes2fname):
    assert not None in prefixes2fname[prefix], "Error, pdb %s has only one partner"%(prefix)
    if len(prefixes2fname[prefix])==2:
      argsList.append({"allRootDir":outDir, "lFname":prefixes2fname[prefix][0],
                       "rFname":prefixes2fname[prefix][1], "idName":prefix, "cleanForDisplay":False, 
                       "areForTrainAndTest": areForTrainAndTest, "removeInputs":False, "isHomoComplex":isHomo})
  Parallel(n_jobs= ncpu, backend="multiprocessing")(delayed(predictOneComplex)(**params) for params in argsList)

  
def predictAllComplexesInList(fnameListFile, outDir, isHomo=False, areForTrainAndTest=True, ncpu=1):
  '''
  fnameListFile
  
  1acb:E,1acb:I
  '''
  
  lPdbIds={}
  rPdbIds={}

  if isinstance(fnameListFile, file):
    f = fnameListFile
  else:
    f = open(fnameListFile)

  for line in f:
    if line.startswith(">") or line.startswith("#"): continue
    prefix1_chains1, prefix2_chains2=  line.split(",")
    prefix1, chains1 = prefix1_chains1.split(":")
    prefix2, chains2 = prefix2_chains2.split(":")

    lPdbIds[prefix1]= prefix1+":"+chains1
    rPdbIds[prefix2]= prefix2+":"+chains2

  argsList=[]
  for prefix in sorted(lPdbIds):
    argsList.append({"allRootDir":outDir, "lFname":None, "rFname":None,
                     "lPdbId":lPdbIds[prefix], "rPdbId":rPdbIds[prefix], "idName":prefix,
                     "cleanForDisplay":False, "areForTrainAndTest": areForTrainAndTest,
                     "removeInputs":False, "isHomoComplex":isHomo})
  f.close()
  Parallel(n_jobs= ncpu, backend="multiprocessing")(delayed(predictOneComplex)(**params) for params in argsList)

if __name__=="__main__":

  from Config import Configuration

  Configuration.update_config_dict(project_config="./configFiles/cmdTool/configFile_pred.cfg", other_dict=dict(force_protocol= None, wdir="./wdir"))
  conf = Configuration()
  parser =  conf.getArgParser()


  command_group = parser.add_mutually_exclusive_group()
  command_group.add_argument("-i", "--inputDir", help="Directory where input files are located. They can be .pdb or .fasta, but they have to"
                                                      "follow Benchmark 5 naming convenctions PREFIX_[lr]_[bu].(pdb|fasta)", type=Configuration.file_path, default=None)
  command_group.add_argument("-f", "--inputFile", help="A .csv file with PDB ids and chains to download", type=argparse.FileType('r'), default=None) #sys.stdin)

  parser.modify_field("wdir", help="Directory where partial results and final results will be saved")
  parser.modify_field("tmp", help="Temporary directory")

  parser.modify_field("ncpu", help="Number of cpus for trainng. Each complex in a cross-validation fold is computed in an indepented worker. NCPU workers are computed in parallel")
  parser.modify_field("N_KFOLD", help="Type of cross validation. -1 for leave-one-complex out, positive values for k= N_KFOLD cross-validation.", _type=int)
  parser.modify_field("scopeFamiliesFname", help="Filename containing the familes of the protein chains of ligand and receptor")

  command_group.add_argument( "--areForTrainAndTest", action="store_true", help="Compute predictions for evaluation purposes. As a requisit, the two partners have to be in bound coordinates")

  args = parser.parse_args()

  assert conf.wdir, "Error, --wdir required"
  assert os.path.isdir( conf.savedModelsPath), "Error, savedModelsPath (%s) does not exists" % conf.savedModelsPath

  if args.get("inputFile"):
    infile = args.get("inputFile")
    assert os.path.isfile(infile.name) or isinstance(infile, file), "Error, %s is not a file"%infile.name
    predictAllComplexesInList(fnameListFile=args.get("inputFile"),outDir=conf.wdir,
                              areForTrainAndTest=args.get("areForTrainAndTest"), ncpu=conf.ncpu)
  elif args.get("inputDir"):
    indir = args.get("inputDir")
    assert os.path.isdir(indir), "Error, %s is not a directory"%indir

    predictAllComplexesInDir(indir, outDir= conf.wdir, areForTrainAndTest=args.get("areForTrainAndTest"), ncpu= conf.ncpu)
  
  
