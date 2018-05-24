import os,sys
from subprocess import Popen, PIPE

from ..StructFeatComputer import StructFeatComputer, FeatureComputerException
from utils import myMakeDir, tryToRemove #utils is at the root of the package

class PSAIAComputer(StructFeatComputer):
  '''
  Extends StructFeatComputer class. Computes PSAIA (asa, hidrofobicity...) for a given complex
  '''
  def __init__(self, rFname, lFname, computedFeatsRootDir= None, statusManager=None):
    '''
      @param rFname: str. path to receptor pdb file
      @param lFname: str. path to ligand pdb file      
      @param computedFeatsRootDir: str. path where features will be stored. If None they will be stored
                                        at default path (assigned in ../Config.py)
      @param statusManager: class that implements .setStatus(msg) to communicate
    '''
    StructFeatComputer.__init__(self, rFname, lFname, computedFeatsRootDir, statusManager=statusManager)

    self.outPathRaw=  myMakeDir(self.computedFeatsRootDir, "PSAIA/rawPSAIA")
    self.outPathProc= myMakeDir(self.computedFeatsRootDir, "PSAIA/procPSAIA")
    
#    self.listFileNameForPSAIA = os.path.join(self.psaiaRootDir, "%s_inputFile.fls")
#    self.configFileName= os.path.join(self.psaiaRootDir, "%s_psa_config.cfg")

    self.listFileNameForPSAIA = os.path.join(self.outPathRaw, "%s_inputFile.fls")
    self.configFileName= os.path.join(self.outPathRaw, "%s_psa_config.cfg")

    self.configForPSAIATemplate='''analyze_bound:	1
analyze_unbound:	0
calc_asa:	1
z_slice:	0.25
r_solvent:	1.4
write_asa:	1
calc_rasa:	1
standard_asa:	%s
calc_dpx:	1
calc_cx:	1
cx_threshold:	10
cx_volume:	20.1
calc_hydro:	1
hydro_file: %s
radii_filename:	%s
write_xml:	0
write_table:	1
output_dir:	%s
''' %(
      os.path.join(self.psaiaRootDir, "amac_data/natural_asa.asa"),
      os.path.join(self.psaiaRootDir, "amac_data/hydrophobicity.hpb"),
      os.path.join(self.psaiaRootDir, "amac_data/chothia.radii"),
      self.outPathRaw
    )

  def cleanAllPSAIA(self):
    '''
      Removes all temporal files generated for PSAIA. Warning: If run when other processes
      are computing PSAIA, it will interfer with them
    '''

    for fname in os.listdir(self.outPathRaw):
      if fname.endswith(".fls") or fname.endswith(".cfg"):
        os.remove(os.path.join(self.outPathRaw,fname))

  def processPSAIA(self, prefixAndChainTypeId):
    '''
      Parses raw output from PSAIA and creates results in tab format.
      @param prefixAndChainTypeId: str. fname prefix of raw psaia results (not taking into account 
                                  the full path but the name itself)
    '''
    stringDict={}
    computed= False
    header=("chainId structResId resName total_ASA b-bone_ASA s-chain_ASA polar_ASA n-polar_ASA total_RASA "+
        "b-bone_RASA s-chain_RASA polar_RASA n-polar_RASA average_DPX s_avg_DPX s-ch_avg_DPX s-ch_s_avg_DPX "
        "max_DPX min_DPX average_CX s_avg_CX s-ch_avg_CX s-ch_s_avg_CX max_CX min_CX Hydrophobicity\n")

    for fname in os.listdir(self.outPathRaw):
      if fname.endswith(".tbl"):
        if fname.startswith(prefixAndChainTypeId) and not computed:
          f=open(os.path.join(self.outPathRaw,fname))
          for i in range(8):
            f.readline()
          for line in f:
            arrayLine=line.split()
            out=[arrayLine[0]]+[arrayLine[6]]+[self.threeLetterAA_to_one(arrayLine[7])]+arrayLine[8:]
            out="\t".join(out)
            try:
              stringDict[arrayLine[0]].append(out)
            except KeyError:
              stringDict[arrayLine[0]]=[out]
          f.close()
          computed=True

    splitName= prefixAndChainTypeId.split("_")
    outNames= []    
    try:

      for chainId in stringDict:

        if len(splitName)==3:
          prefix, chainType, unbound = splitName
          outName= os.path.join(self.outPathProc, prefix+"_"+chainType+"_"+chainId+"_u.psaia.tab")
        else:
          outName= os.path.join(self.outPathProc, prefixAndChainTypeId+"_"+chainId+"_u.psaia.tab")
        outNames.append(outName)
        outFile=open(outName,"w")
        outFile.write(header)
        outFile.write("\n".join(stringDict[chainId]))
        outFile.close()
    except (KeyboardInterrupt, Exception):
      for outName in outNames:
        print("Exception happend computing %s"%outName)
        tryToRemove(outName)    
      raise            


  def checkIfAlreadyComputed(self, prefixAndChainTypeId):
    sharedPrefix= "_".join( prefixAndChainTypeId.split("_")[:2])
    for fname in sorted(os.listdir(self.outPathProc)): #remove old psaia runs
      if fname.startswith(sharedPrefix):
        return True
    return False
    
  def computeOneFile(self, fileName):
    '''
      Computes PSAIA for a given pdb file
      @param fileName: str. fname to pdb file
    '''
    try:
      prefixAndChainTypeId = (os.path.split(fileName)[-1]).split(".pdb")[0]
      for fname in os.listdir(self.outPathRaw): #remove old psaia runs
        if fname.startswith(prefixAndChainTypeId):
          os.remove( os.path.join(self.outPathRaw,fname))
      if self.checkIfAlreadyComputed( prefixAndChainTypeId):
        print("PSAIA already computed")
        return 0
      f= open(self.listFileNameForPSAIA%prefixAndChainTypeId,"w")
      f.write(fileName)
      f.close()

      f= open(self.configFileName%prefixAndChainTypeId,"w")
      f.write(self.configForPSAIATemplate)
      f.close()

      proc = Popen([os.path.join(self.psaiaRootDir,"psa"),(self.configFileName%prefixAndChainTypeId),
                    self.listFileNameForPSAIA%prefixAndChainTypeId],
                    stdin= PIPE, stdout=PIPE, stderr=PIPE)
      output=  proc.communicate(input="y\n")
      if output== None or output[1]!="" or "There was an error in PDB" in output[0]:
        print(output)
        print ("Error when computing PSAIA for %s"%fileName)
        raise FeatureComputerException("Error when computing PSAIA for %s"%fileName)
      else:      
        self.processPSAIA(prefixAndChainTypeId)
      return 0
    finally:
      tryToRemove(self.listFileNameForPSAIA%prefixAndChainTypeId)
      tryToRemove(self.configFileName%prefixAndChainTypeId)

def testModule():
  rFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/pdbFiles/1A2K_r_u.pdb"
  lFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/pdbFiles/1A2K_l_u.pdb"
  
  computedFeatsRootDir="/home/rsanchez/tmp/tmpRRI"
  psaiaComp= PSAIAComputer(rFname, lFname, computedFeatsRootDir)
  psaiaComp.computeFun()
  raise FeatureComputerException("Debug Mode. Comment out testModule() line 128")
  
if __name__=="__main__":
  testModule()
  
  pdbsIndir= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3"
  computedFeatsRootDir= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/computedFeatures"

#  pdbsIndir= "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/b5structures/"
#  computedFeatsRootDir= "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/computedFeatures"

#  pdbsIndir= "~/Tesis/rriPredMethod/data/morePDBsDimers/benchPDBs/newBenchPDB_splited_by_chains"
#  computedFeatsRootDir= "~/Tesis/rriPredMethod/data/morePDBsDimers/computedFeatures"

  StructFeatComputer.computeFeaturesAllComplexes(PSAIAComputer,pdbsIndir= pdbsIndir ,computedFeatsRootDir= computedFeatsRootDir)


