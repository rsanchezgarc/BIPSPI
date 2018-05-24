from __future__ import absolute_import
import os, sys
from subprocess import Popen, PIPE
from Bio.PDB.Polypeptide import is_aa
from computeFeatures.structStep.myPDBParser import myPDBParser as PDBParser

from ..StructFeatComputer import StructFeatComputer
from utils import myMakeDir, tryToRemove #utils is at the root of the package
from bigExceptions import NoValidPDBFile

def decodeFun(oneStr):
  if sys.version_info >= (3, 0):
    return oneStr.decode("utf-8")
  else:
    return oneStr

class DsspComputer(StructFeatComputer):
  '''
  Extends StructFeatComputer class. Computes Secondary Structure with DSSP
  '''
  DSSP_HEADER= "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA \n"
  def __init__(self, rFname, lFname, computedFeatsRootDir=None, statusManager=None):
    '''
      @param rFname: str. path to receptor pdb file
      @param lFname: str. path to ligand pdb file      
      @param computedFeatsRootDir: str. path where features will be stored
      @param statusManager: class that implements .setStatus(msg) to communicate
    '''
    
    StructFeatComputer.__init__(self, rFname, lFname, computedFeatsRootDir, statusManager=statusManager)

    self.outPathRaw=  myMakeDir(self.computedFeatsRootDir, "DSSP/rawDSSP")
    self.outPathProc= myMakeDir(self.computedFeatsRootDir, "DSSP/procDSSP")

    self.dsspBinPath= os.path.join(self.dsspRootDir,"mkdssp")


  def processDSSP(self, prefixAndChainTypeId, pdbStruct, dsspFName):
    '''
      Parses raw output from DSSP and creates results in tab format.
      @param prefixAndChainTypeId: str. fname prefix of raw psaia results (not taking into account 
                                  the full path but the name itself)
      @param pdbStruct: Bio.PDB.Structure. Structure of the pdb that is being analyzed
      @param dsspFName: str. fname of raw Dssp output
    '''
    struct= pdbStruct
    resDict= {}
    for chain in pdbStruct[0]:
      chainId=chain.get_id()
      if chainId==" ": chainId="*"
      resDict[chainId]= set([])
      for res in chain:
        if not is_aa(res, standard=True): continue
        resId= str(res.get_id()[1])
        if res.get_id()[2]!=" ":
          resId+= res.get_id()[2]
        try:
          resDict[chainId].add( (resId, self.threeLetterAA_to_one(res.resname)) )
        except KeyError: continue
  
    chain=None
    stringDict={}
    computed= False
    header=("#Levels: 3:H;B;E;G;I;T;S;Z\n"+ 
            "chainId structResId resName 2ndStruct\n")

    f=open(dsspFName)
    error= True
    for line in f:
      if line.startswith('  #  RESIDUE '):
        error= False
        break
    if error== True:
      raise Exception("Error processing DSSP file "+dsspFName)

    for line in f:
      secStruc=line[16]
      if secStruc==" ":
        secStruc="Z"
      lineArray= line.split()
      chain= line[11]
      resLetter= line[13]
      if resLetter.islower(): resLetter="C"
      resInd=lineArray[1]
      if chain==" ": chain="*"
      if resInd.startswith("!"):  #Chain changed
        continue
      resDict[chain].discard( (resInd, resLetter) )
      out="\t".join([chain, resInd, resLetter, secStruc])
      try:
        stringDict[chain].append(out)
      except KeyError:
        stringDict[chain]=[out]
    #To fill residues for which we do not have enough information

    for chainId in resDict:
      for resId,resLetter in sorted(resDict[chainId]):
        out="\t".join([chainId, resId, resLetter, "Z"])
        try:
          stringDict[chain].append(out)
        except KeyError:
          stringDict[chain]=[out]
    f.close()
    splitName= prefixAndChainTypeId.split("_")    
    for chainId in stringDict:
      if len(splitName)==3:
        prefix, chainType, unbound = splitName
        outName= os.path.join(self.outPathProc, prefix+"_"+chainType+"_"+chainId+"_u.dssp.tab")
      else:
        outName= os.path.join(self.outPathProc, prefixAndChainTypeId+"_"+chainId+"_u.dssp.tab")      
##      print(outName)
      try:      
        outFile=open(outName,"w")
        outFile.write(header)
        outFile.write("\n".join(stringDict[chainId]))
        outFile.close()
      except (KeyboardInterrupt, Exception):
        print("Exception happend computing %s"%outName)
        tryToRemove(outName)    
        raise  

  def computeOneFile(self, pdbFName):
    '''
      Computes DSSP for a given pdb file
      @param pdbFName: str. fname to pdb file
    '''
    parser= PDBParser(QUIET=True)
    struct= parser.get_structure("pdbStruct", pdbFName)
    prefixAndChainTypeId = self.getExtendedPrefix(pdbFName)
    rawDsspOutName= os.path.join(self.outPathRaw, prefixAndChainTypeId+".dssp.tab")

    proc = Popen([self.dsspBinPath,'-i', pdbFName ,'-o',rawDsspOutName],
                   stdin= PIPE, stdout=PIPE, stderr=PIPE)
    output=  proc.communicate()
    if output== None or decodeFun(output[1])!="":
##    no atoms read before TER record
      print("Error when computing DSSP: %s"% pdbFName)
      print( output)
      ## 'no atoms read before TER record \nTER    
      if not decodeFun(output[1]).startswith('no atoms read before TER record'):
        self.createFileForError(struct, rawDsspOutName)
    self.processDSSP(prefixAndChainTypeId, struct, rawDsspOutName)
#    raw_input("enter to continue")
    return None

  def createFileForError(self, pdbStruct, outName):
    '''
      Creates a fake DSSP raw output generated when DSSP fails. All residues will be assigned secStruc= Z
      @param pdbStruct: Bio.PDB.Structure. Structure of the psb that is being analyzed
      @param outName: str. output fname
    '''
    oneResLine= "%5d%5d%2s%2s %2s\n"
    try:
      f= open(outName,"w")
      f.write( DsspComputer.DSSP_HEADER)
      if len(pdbStruct)==0:
        raise NoValidPDBFile("No valid pdb File. There are no models contained")
      for chain in pdbStruct[0]:
        for i, res in enumerate(chain):
          if not is_aa(res): continue
  ##        print i,res,res.get_id()
          seqIndex= i+1
          structIndex = res.get_id()[1]
          letter= self.threeLetterAA_to_one(res.resname)
          fakeSecStruct="Z"
          fakeCharacters1= tuple("f"*7)
          fakeDigits1= tuple([0,0,"f",0])
          fakeStrs= tuple("f"*4)
          fakeFloats= tuple(elem + 0.0 for elem in range(8))
  ##        print ((seqIndex,structIndex, chain.get_id(), letter, fakeSecStruct)+
  ##                            fakeCharacters1+fakeDigits1+fakeStrs+fakeFloats)
  ##        print oneResLine%( (seqIndex,structIndex, chain.get_id(), letter, fakeSecStruct)+
  ##                            fakeCharacters1+fakeDigits1+fakeStrs+fakeFloats)

          f.write(oneResLine%(seqIndex,structIndex, chain.get_id(), letter, fakeSecStruct))
      f.close()
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing %s"%outName)
      tryToRemove(outName)    
      raise        
    return 0


def testModule():
#  rFname="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/1A2K_r_u.pdb"
#  lFname="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/1A2K_l_u.pdb"
  
  rFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1A22.pdb"
  lFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1ABR.pdb"
  
  computedFeatsRootDir="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/tmpComputedFeatures"
  dsspComp= DsspComputer(rFname, lFname, computedFeatsRootDir)
  dsspComp.computeFun()
  raise FeatureComputerException("Debug Mode. Comment out testModule() line 128")                
                    
if __name__=="__main__":
  testModule()
  
##  pdbsIndir= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3"
##  computedFeatsRootDir= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/computedFeatures"

  pdbsIndir= "~/Tesis/rriPredMethod/data/bench5Data/b5structures/"
  computedFeatsRootDir= "~/Tesis/rriPredMethod/data/bench5Data/computedFeatures"


  StructFeatComputer.computeFeaturesAllComplexes(DsspComputer,pdbsIndir= pdbsIndir ,computedFeatsRootDir= computedFeatsRootDir)


