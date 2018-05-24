'''
This module allows for the prediction of Residue-Residue Contacts and Binding Sites given 2 protein partners
'''
from __future__ import print_function
import os, sys
import shutil
import computeFeatures.computeFeatsForPdbs as pComCode
import computeFeatures.common.computeContactMap as pCMAP
import codifyComplexes.codifyOneComplex as pCodifyOne
import trainAndTest.predictOneCodifiedComplex as pPredict
import gzip
import re
import traceback
from joblib import Parallel, delayed
from utils import myMakeDir
from pythonTools.downloadPdb import downloadPDB
from computeFeatures.seqStep.SeqFeatComputer import parseSeq
from computeFeatures.structStep.StructFeatComputer import moveAndWriteAsPDBIfMmcif
from bigExceptions import NoAvailableForDownloadPDB, NoValidPDBFile, MyException
from Config import Configuration


DEFAULT_PREFIX="PRED"


def prepareInput(fname, partnerStr, finalInputPath, chainType, idName, methodProtocol, removeInputs=False):
  '''
    Takes a pdb or fasta file in fname 
    or 
    a sequence or pdbId in partnerStr 
    and preprocess input (and download it in the case pdbId is provided but no pdbFile )
    
    @param fname: str. A path to a pdb file or a fasta file
    @param partnerStr: str. A  sequnce of amino acids or pdbId that can contain 
                            info about chain and/or biounit: pdb[:chain][_bioUnit]. e.g. 1A2K, 1A2K:A 1A2K_1, 1A2K:B_1
    @param finalInputPath: str. Where preprocessed input will be saved
    @param chainType: str. "l" or "r"
    @param idName: str. An id that will be used for files prefix
    @param methodProtocol str: "mixed" for mixed model or "seq_pred" for sequnce
    @param removeInputs: boolean. True if input files want to be removed after preprocessed versions was generated
  '''

  pdbId= None
  if methodProtocol.startswith("mixed"):
    newFname= os.path.join(finalInputPath,"%s_%s_u.pdb"%(idName, chainType))
    if fname==None and partnerStr!= None:
      pdbId, chainId, biounit = [None]*3
      print(partnerStr)
      matchObj= re.match(r"(^\d[a-zA-Z0-9]{3})(:[a-zA-Z0-9])?(_\d*)?$", partnerStr)
      if matchObj:
        pdbId= matchObj.group(1)
        chainId= matchObj.group(2)
        biounit= matchObj.group(3)
        if not chainId is None: chainId= chainId[1:]
        if biounit is None or len(biounit)==1: biounit=None
        else:
          biounit= int(biounit[1:])
          if biounit==0: biounit==None
      else:
        raise NoAvailableForDownloadPDB("Error: bad pdb provided %s"%(partnerStr.split("_")[0]))
      try:
        fname = downloadPDB( pdbId, finalInputPath, chainId= chainId, bioUnit=biounit)
      except ValueError as e:
        traceback.print_exc()
        raise NoAvailableForDownloadPDB( str(e)+"-> %s"%(partnerStr))
      try:
        if not removeInputs:
          shutil.copyfile(fname, newFname)
        else:
          os.rename(fname, newFname)
      except OSError:
        NoAvailableForDownloadPDB("Error. pdb %s was not recover from wwpdb.org"%partnerStr)
    else:
      if not moveAndWriteAsPDBIfMmcif(fname, newFname, removeInput= removeInputs): 
        raise NoValidPDBFile("Error. It was not possible to parse your pdb file %d"%( 1 if chainType=="l" else 2))
  else:
    newFname= os.path.join(finalInputPath,"%s_%s_u.fasta"%(idName, chainType))
    print(fname, newFname)
    if partnerStr!= None and len(partnerStr)>0:
      seq= parseSeq(partnerStr)
      with open(newFname, "w") as f:
        f.write(">%s\n%s"%(chainType, seq))
    else:
      if not removeInputs:
        shutil.copyfile(fname, newFname)
      else:
        os.rename(fname, newFname)
  return newFname, pdbId
  
def computeFeatures(allRootDir, lFname=None, rFname=None, lpartnerStr=None, rpartnerStr=None, idName=None,
                    methodProtocol="mixed", areForTrainAndTest=False, removeInputs=False, statusManager=None):
  '''
    @param allRootDir: str. A path where all files will be saved
    @param predictionType str: "m" for mixed model or "s" for sequnce
    @param lFname: str. A path to a fasta or pdb file
    @param rFname: str. A path to a fasta or pdb file      
    @param lpartnerStr: str. A pdbId (4 letters) or a sequnce of amino acids
    @param rpartnerStr: str. A pdbId (4 letters) or a sequnce of amino acids 
    @param idName: str. An id that will be used for files prefix
    @param methodProtocol: str. "mixed" or "seq_pred"
    @param areForTrainAndTest: boolean. True if ligand and receptor are in interacting coordinates and thus,
                            contact maps are needed to be computed in order to perform evaluation.
                            False if ligand and receptor coordinates are not related and thus, 
                            evaluation does not makes sense.
    @param removeInputs: boolean. True if input files want to be removed after preprocessed versions was generated
    @param statusManager: class that implements .setStatus(msg) to communicate 
  '''
  
  computedFeatsRootDir= myMakeDir(allRootDir, "computedFeatures")
  pdbsPath= myMakeDir(allRootDir, "finalInput")
  if not lFname is None and lFname==rFname:
    raise MyException("partner1 file and partner2 file must be different")
  newLFname, pdbIdL= prepareInput(lFname, lpartnerStr, pdbsPath, "l", idName, methodProtocol, 
                                    removeInputs= removeInputs and not areForTrainAndTest)
  newRFname, pdbIdR= prepareInput(rFname, rpartnerStr, pdbsPath, "r", idName, methodProtocol,
                                    removeInputs= removeInputs and not areForTrainAndTest)

  pComCode.computeFeaturesOneComplex( lFname= newLFname, rFname=newRFname, lPdbId=pdbIdL, rPdbId= pdbIdR,
                                      computedFeatsRootDir= computedFeatsRootDir, methodProtocol=methodProtocol,
                                      areForTrainAndTest=areForTrainAndTest, statusManager= statusManager)
  return computedFeatsRootDir

def codifyStep(allRootDir, stepType, computedFeatsRootDir, feedbackPaths=None, idName=None):

    
#  wholeComplexOutPath= os.path.join(allRootDir, "codifiedInput")
  wholeComplexOutPath= None  #We don't want the ready to used data to be saved
  oneComplexCod= pCodifyOne.OneComplexCodifier( dataRootPath= computedFeatsRootDir, sampledOutPath=None,
                                                wholeComplexOutPath=wholeComplexOutPath,
                                                environType=stepType, feedback_paths= feedbackPaths)

  codifiedObject= oneComplexCod.codifyComplex( idName)
  print(codifiedObject, codifiedObject.pairsDirect.shape)
  return codifiedObject
  
def predict(allRootDir, stepType, complexCodifiedObject, isLastStep=False):
  predictOutputPath= myMakeDir(allRootDir, stepType)
  predictor= pPredict.ComplexPredictor(stepType)
  outPath, resultsDfs= predictor.predictOneComplex(complexCodifiedObject, 
                                            os.path.join(predictOutputPath,"%s.res.tab"%complexCodifiedObject.prefix),
                                            isLastStep)
  return predictOutputPath, resultsDfs
  

def cleanDirectory(dirName):
  print("cleaning %s"%dirName )
  for root, dirs, files in os.walk(dirName, topdown=False):
    for name in files:
      fname= os.path.join(root, name)
      with open(fname, 'rb') as f_in:
        with open(fname+'.gz', 'wb') as f_out:
          with gzip.GzipFile(fname, 'wb', fileobj=f_out) as f_out:
              shutil.copyfileobj(f_in, f_out)
      os.remove(fname)
    for name in dirs:
      dirName= os.path.join(root, name)
      if name.startswith("computedFeatures"):
        shutil.rmtree(dirName)
       
def predictOneComplex(predictionType, allRootDir, lFname=None, rFname=None, lpartnerStr=None, rpartnerStr=None, 
          idName=DEFAULT_PREFIX, cleanForDisplay=False, areForTrainAndTest=False, removeInputs=False, statusManager=None):
  '''
    obtains Residue-Residue Contact predictions and binding site predictions for the sequences or structures contained in
            the files lFname, rFname  or indicated in lpartnerStr, rpartnerStr (as a sequence string or as a pdbId)
            
    @param predictionType str: "m" for mixed model or "s" for sequnce
    @param allRootDir: str. A path where all files will be saved
    @param lFname: str. A path to a fasta or pdb file.
    @param rFname: str. A path to a fasta or pdb file.
    @param lpartnerStr: str. A  sequnce of amino acids for partner l or pdbId (4 letters). Can contain info about 
                             chain and/or biounit: pdb[:chain][_bioUnit]
    @param rpartnerStr: str. A  sequnce of amino acids for partner r or pdbId (4 letters). Can contain info about 
                             chain and/or biounit: pdb[:chain][_bioUnit]
    @param idName: str. An id that will be used for files prefix
    @param cleanForDisplay: Boolean. If true categ will be removed from output Dataframes and l will be substituted by
                                      partner 1 and r by parnter 2.
    @param areForTrainAndTest: boolean. True if ligand and receptor are in interacting coordinates and thus,
                            contact maps are needed to be computed in order to perform evaluation.
                            False if ligand and receptor coordinates are not related and thus, 
                            evaluation does not makes sense. If False 3dcons will be searched by seq. Else 3dcons 
                            will be serched first by pdb id and then by seq
    @param removeInputs: boolen. If True, inputs will be removed        
    @param statusManager: class that implements .setStatus(msg) to communicate status
    @return predictOutPath, resultsDfs
             predictOutPath: path were results were saved
             resultsDfs: (resultsPairs, resultsL, resultsR).
  '''
          
  assert (lFname!=None or lpartnerStr!=None) and (rFname!=None or rpartnerStr!=None) and not (lFname==None and lpartnerStr==None) \
          and not (rFname==None and rpartnerStr==None)
  if predictionType=="m":
    predictionType="mixed"
  elif predictionType=="s":
    predictionType="seq_pred"
  else:
    raise ValueError("bad predictionType selected. It can be m or s")
    
  if not lFname is None: lFname= os.path.expanduser(lFname)
  if not rFname is None: rFname= os.path.expanduser(rFname)
  
  if areForTrainAndTest:
    if not lFname is None:
      prefixL= os.path.basename(lFname).split(".")[0].split("_")[0]
      lpartnerStr= prefixL
    else:
      prefixL= lpartnerStr.split(":")[0]
    if not rFname is None:
      prefixR= os.path.basename(rFname).split(".")[0].split("_")[0]
      rpartnerStr= prefixR
    else:
      prefixR= rpartnerStr.split(":")[0]
      
    assert prefixL==prefixR, "Error, different prefixes when areForTrainAndTest=False. %s %s"%(prefixL, prefixR)
    
  numTotalSteps= 5 if predictionType[0]=="m" else 3
  if not statusManager is None: statusManager.setStatus("STEP 1/%d: Computing features"%(numTotalSteps))
  computedFeatsRootDir= computeFeatures(allRootDir, lFname, rFname, lpartnerStr, rpartnerStr, idName, 
                                      methodProtocol=predictionType, areForTrainAndTest=areForTrainAndTest,
                                      removeInputs= removeInputs, statusManager= statusManager)
  computedFeatsRootDir= os.path.join(allRootDir, "computedFeatures")
  predsRootDir= myMakeDir(allRootDir, "preds")
  if predictionType=="seq_pred":
    if not statusManager is None: statusManager.setStatus("STEP 2/%d: Data encoding"%(numTotalSteps))
    codifiedObject= codifyStep(allRootDir, "seq", computedFeatsRootDir, idName=idName,)
    
    if not statusManager is None: statusManager.setStatus("STEP 3/%d: Predictions calculation"%(numTotalSteps))
    predictOutPath, resultsDfs= predict(predsRootDir, "seq_train", codifiedObject, isLastStep=cleanForDisplay)
    print("seq step done.")
  elif predictionType=="mixed":
    if not statusManager is None: statusManager.setStatus("STEP 2/%d: Data encoding"%(numTotalSteps))
    codifiedObject= codifyStep(allRootDir, "mixed", computedFeatsRootDir, idName=idName)
    
    if not statusManager is None: statusManager.setStatus("STEP 3/%d: Predictions calculation"%(numTotalSteps))
    predictOutPathMixed, __= predict(predsRootDir, "mixed", codifiedObject)
    
    if not statusManager is None: statusManager.setStatus("STEP 4/%d: Feedback data encoding"%(numTotalSteps))
    codifiedObject= codifyStep(allRootDir, "mixed_2", computedFeatsRootDir, 
                               [ predictOutPathMixed], idName=idName)
                               
    if not statusManager is None: statusManager.setStatus("STEP 5/%d: Feedback predictions calculation"%(numTotalSteps))
    predictOutPath, resultsDfs= predict(predsRootDir, "mixed_2", codifiedObject, isLastStep=cleanForDisplay)
  else:
    raise ValueError("predictionType must be m or s")
  
  if cleanForDisplay==True:
    for df in resultsDfs: 
      if "categ" in df.columns:
        df= df.drop("categ",axis=1)
      for colname in df.columns:
        if colname.startswith("chainId"):
          df[colname]= df[colname].map(lambda x: x if  x!="*" else " ")
    cleanDirectory(allRootDir)
  return predictOutPath, resultsDfs
  
def predictAllComplexesInDir(inDir, outDir, areForTrainAndTest=True, ncpu=1):
  prefixes= {}
  fileExtensions =set([])
  for fname in sorted(os.listdir(inDir)):
    fname_split= fname.split(".")
    prefix, chainType, boundState= fname_split[0].split("_")
    fileExtensions.add( fname_split[-1] )
    if not prefix in prefixes:
      prefixes[prefix]=[None, None]
    if boundState=="u":      
      prefixes[prefix][0 if chainType=="l" else 1]= os.path.join(inDir, fname)
  oneFname= list(prefixes.values())[0][0]
  assert oneFname.endswith(".pdb") or oneFname.endswith(".fasta"), "Error, files in directory must be .pdb or .fasta"
  assert len(fileExtensions)==1, "Several extensions files in input directory, all files must be .pdb or .fasta"
  if oneFname.endswith(".fasta"):
    areForTrainAndTest=False
    methodType= "s"
  else:
    methodType="m"
  argsList=[]
  for prefix in sorted(prefixes):
    assert not None in prefixes[prefix], "Error, pdb %s has only one partner"%(prefix)
    if len(prefixes[prefix])==2:
      argsList.append({"predictionType": methodType, "allRootDir":outDir, "lFname":prefixes[prefix][0],
                       "rFname":prefixes[prefix][1], "lpartnerStr":None, "rpartnerStr":None,
                        "idName":prefix, "cleanForDisplay":False, "areForTrainAndTest": areForTrainAndTest,
                        "removeInputs":False})
  
  Parallel(n_jobs= ncpu)(delayed(predictOneComplex)(**params) for params in argsList)
  
def testPredictOne():
  lFname="~/Tesis/rriPredMethod/data/develData/pdbFiles/1ACB_l_u.pdb"
  rFname="~/Tesis/rriPredMethod/data/develData/pdbFiles/1ACB_r_u.pdb"
#  lFname="~/Tesis/rriPredMethod/data/develData/pdbFiles/1ACB_l_u.pdb"
#  rFname="~/Tesis/rriPredMethod/data/develData/pdbFiles/1A2K_l_u.pdb"

#  rFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/inputFasta/seq1_r_u.fasta"
#  lFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/inputFasta/seq1_l_u.fasta"
#  rFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/inputFasta/1ACB_r_u.fasta"
#  lFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/inputFasta/1ACB_l_u.fasta"
  
#  lpartnerStr="4OV6:E"
#  rpartnerStr="4OV6:D"
  
#  lFname=None  
#  rFname=None
  lpartnerStr=None
  rpartnerStr=None

  methodType="m"
  outName, (p,l,r)= predictOneComplex(predictionType= methodType, allRootDir="/home/rsanchez/tmp/tmpRRI/wdir/trial_ori_1", 
                                        lFname=lFname, rFname=rFname, lpartnerStr= lpartnerStr, rpartnerStr=rpartnerStr, 
                                        areForTrainAndTest=False, removeInputs=False)
  
if __name__=="__main__":
  '''
  usage example: python predictComplexes.py ~/Tesis/rriPredMethod/data/develData/inputFasta/ ~/tmp/tmpRRI/wdir/
  '''
#  testPredictOne(); sys.exit()
  conf= Configuration()
  usage= ("Error. usage: python predictComplexes.py inDir outDir\n"+
          "ne.g. python predictComplexes.py ~/my/path/to/L-R-inFiles ~/my/path/to/results\n"+
          "inDir must contain:\n\tone pdb file for ligand and one pdb file for receptor named \n"+
          "\t\tprefix_l_u.pdb, prefix_r_u.pdb, where prefix is a string (id) for each complex\n"+
          "\tone fasta file for ligand and one fasta file for receptor named \n"+
          "\t\tprefix_l_u.fasta, prefix_r_u.fasta, where prefix is a string (id) for each complex\n")
  print(sys.argv)
  if len(sys.argv)<3:
    print(usage)
    sys.exit(1)
  inDir= sys.argv[1]
  outDir= sys.argv[2]
  if len(sys.argv)>3:
    print(usage)
    sys.exit(1)
  predictAllComplexesInDir(inDir, outDir , areForTrainAndTest=False, ncpu= conf.ncpu)
  
  
