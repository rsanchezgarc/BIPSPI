from __future__ import absolute_import, print_function
import os

from subprocess import Popen, PIPE
from .structToolManager import StructToolManager, FeatureComputerException
from utils import myMakeDir, tryToRemove, tryToCleanDir

class PSAIAComputer(StructToolManager):
  '''
  Extends StructToolManager class. Computes PSAIA (asa, hidrofobicity...) for a given complex
  '''
  
  VAR_LIST= ['total_ASA', 'b-bone_ASA', 's-chain_ASA', 'polar_ASA', 'n-polar_ASA', 'total_RASA', 
             'b-bone_RASA', 's-chain_RASA', 'polar_RASA', 'n-polar_RASA', 'average_DPX', 's_avg_DPX', 
             's-ch_avg_DPX', 's-ch_s_avg_DPX', 'max_DPX', 'min_DPX', 'average_CX', 's_avg_CX',
             's-ch_avg_CX', 's-ch_s_avg_CX', 'max_CX', 'min_CX', 'Hydrophobicity']
                 
  def __init__(self, computedFeatsRootDir, statusManager= None):
    '''
      :param computedFeatsRootDir: str. path where features will be stored. If None they will be stored
                                        at default path (assigned in ../Config.py)
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''

    StructToolManager.__init__(self, computedFeatsRootDir, statusManager)

    self.outPathRaw=  myMakeDir(self.computedFeatsRootDir, "PSAIA/rawPSAIA")
    self.outPath= myMakeDir(self.computedFeatsRootDir, "PSAIA/procPSAIA")
    
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

  def getFNames(self, extendedPrefix):
    '''
    Returns a list of feature files obtained from psaia
    :param extendedPrefix. prefix for output fnames.
    :return list of fnames: [ fname1, fname2, ...]
    '''
    prefix, chainType  = self.splitExtendedPrefix(extendedPrefix )[:2]
    return [os.path.join(self.outPath, prefix+"_"+chainType+"_.psaia.tab.gz") ]
  
  def computeOneFile(self, pdbFname, struct):
    '''
      Computes PSAIA for a given pdb file
      :param pdbFname: str. fname to pdb file
      :param struct: ignored
           
    '''
    assert isinstance(pdbFname, str), "Error, PSAIA computeOneFile first argument is a path to pdb file (str). given %s"%pdbFname
    prefixExtended= self.getExtendedPrefix(pdbFname)
    prefix, chainType= self.splitExtendedPrefix(prefixExtended)[:2]

    if self.checkAlreayComputed(prefixExtended):
      print("Psaia already computed for %s"%prefixExtended)
      return 0
    print("launching PSAIA over %s"%prefixExtended)

    uncompressFileName= self.uncompressFile(pdbFname, self.tmp)
    try:
      with open(self.listFileNameForPSAIA%prefixExtended,"w") as f:
        f.write(uncompressFileName)

      with open(self.configFileName%prefixExtended,"w") as f:
        f.write(self.configForPSAIATemplate)

      proc = Popen([os.path.join(self.psaiaRootDir,"psa"),(self.configFileName%prefixExtended),
                    self.listFileNameForPSAIA%prefixExtended],
                    stdin= PIPE, stdout=PIPE, stderr=PIPE)
      output=  proc.communicate(input="y\n")
      if output== None or output[1]!="" or "There was an error in PDB" in output[0]:
        print(output)
        print ("Error when computing PSAIA for %s"%pdbFname)
        raise FeatureComputerException("Error when computing PSAIA for %s"%pdbFname)
      else:      
        self.processPSAIA(prefixExtended)
      return 0
    except (Exception, KeyboardInterrupt):
      self.tryToRemoveAllFnames(prefixExtended)
      raise
    finally:
      tryToRemove(self.listFileNameForPSAIA%prefixExtended)
      tryToRemove(self.configFileName%prefixExtended)
      tryToRemove(uncompressFileName)
      tryToCleanDir(self.outPathRaw, prefixExtended, rootDataDir=self.computedFeatsRootDir)
      
  def processPSAIA(self, prefixExtended):
    '''
      Parses raw output from PSAIA and creates results in tab format.
      :param prefixExtended: str. fname prefix of raw psaia results (not taking into account
                                  the full path but the name itself). E.g. 

    '''
    prefix, chainType= self.splitExtendedPrefix(prefixExtended)[:2]
    dataDict= {}    
    for fname in os.listdir(self.outPathRaw):
      if fname.endswith(".tbl"):
        if fname.startswith(prefixExtended):
          with open(os.path.join(self.outPathRaw,fname)) as f:
            for i in range(8):
              f.readline()
            for line in f:
              arrayLine=line.split()
              chainId= arrayLine[0]
              if chainId=="*": chainId= type(self).UNKNOWN_CHAIN_SYMB
              chainId_resIdStr_resName=  (chainId, arrayLine[6], self.threeLetterAA_to_one(arrayLine[7]))
              if chainId_resIdStr_resName[-1] is None: continue
              record= list(chainId_resIdStr_resName) + arrayLine[8:]
              record= " ".join(record)
              try:
                dataDict[chainId_resIdStr_resName[0]].append(record)
              except KeyError:
                dataDict[chainId_resIdStr_resName[0]]=[record]
              
          break
    if len(dataDict)==0: raise ValueError("Error, psaia  was not correctly processed")    
    self.writeResultsFromDataDictSingleChain(dataDict, outName= self.getFNames(prefixExtended)[0])


