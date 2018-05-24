from __future__ import absolute_import, print_function

import os
from joblib import delayed, Parallel

from Config import Configuration
from utils import myMakeDir, tryToRemove #utils is at the root of the package

class FeaturesComputer(Configuration):
  '''
  Abstract class. It will be extended with different features computer classes, p.e, PSAIA computer, VORONOI computer...
  Intended to be used for computing one type of features each
  '''
  def __init__(self, rFname, lFname, computedFeatsRootDir= None, statusManager=None):
    '''
      @param rFname: str. path to receptor pdb file
      @param lFname: str. path to ligand pdb file
      @param computedFeatsRootDir: str. path where features will be stored
      @param statusManager: class that implements .setStatus(msg) to communicate
    '''

    Configuration.__init__(self)  # Load configuration parameters such as path to programs
    self.statusManager= statusManager
    if computedFeatsRootDir!= None:
      self.computedFeatsRootDir= computedFeatsRootDir
    self.computedFeatsRootDir= os.path.expanduser(self.computedFeatsRootDir) #Creates root path where features will be saved
    myMakeDir(self.computedFeatsRootDir)
    self.rFname= rFname
    self.lFname= lFname

    if not os.path.isdir(self.computedFeatsRootDir):
      os.mkdir(self.computedFeatsRootDir)

  def reportStatus(self, msg):
    if not self.statusManager is None:
      self.statusManager.appendToStatus(msg)
    
  def getExtendedPrefix(self, fname):
    '''
      Given a filename, obtains its unambiguous id
      @param fname: str. A filename. pe. "/path/to/file/1A2K_l_u.pdb"
      @return unambiguous id: str. pe "1A2K_l_u"  This id will be the prefix of output names
    '''
    return os.path.split(fname)[-1].split(".")[0]
    
  def computeFun(self):
    '''
      abstract method.
      This method will be used to compute the features of a given complex (form by files self.rFname and self.lFname)
    '''
    return None


  @staticmethod 
  def computeFeaturesAllComplexes(OneFeaturesComputerClass, pdbsIndir, computedFeatsRootDir, classArgs={}, ncpu=1):
    '''
      Computes one type of feature over all complexes that are inside pdbsIndir.
      
      @param OneFeaturesComputerClass: FeaturesComputer. class to use for compute one kind of features
      @param pdbsIndir: str. path where pdb files to be computed are located. There must be 2 pdb files per complex
                             To distinguish them _l_ or _r_ infixed are used. P.e: "1A2K_l_u.pdb"  and "1A2K_r_u.pdb".
                             pdb files ended with "b.pdb" will be skipped.
                             
      @param computedFeatsRootDir: str. path where features will be stored
      @param classArgs: Dict. The arguments that will be passed to OneFeaturesComputerClass()
      @param ncpu: the number of subprocess to use in parallel (parallelism at complex level)
    '''

    pdbsIndir= os.path.expanduser(pdbsIndir)
    computedFeatsRootDir= os.path.expanduser(computedFeatsRootDir)
    ConfigObject= Configuration()
    if pdbsIndir== None:
      pdbsIndir= ConfigObject.pdbsIndir

    if computedFeatsRootDir== None:
      computedFeatsRootDir= ConfigObject.computedFeatsRootDir

    fnames= {}
    for fname in sorted(os.listdir(pdbsIndir)):
      if not fname.endswith(".pdb"): continue  #skip no pdb files
      if not fname.endswith("b.pdb"):
        prefix= fname.split("_")[0]      
        if "_r_" in fname or "_l_" in fname:
          if prefix not in fnames:
            fnames[ prefix]= [None, None]
          if "_r_" in fname:
            fnames[ prefix][0]= os.path.join(pdbsIndir, fname)
          if "_l_" in fname:
            fnames[ prefix][1]= os.path.join(pdbsIndir, fname)
        else:
          fnames[ prefix].append( os.path.join(pdbsIndir, fname))
          
    if len(fnames)<1: raise ValueError("There are not files to be processed")
    for prefix in fnames: # check for errors
      if len( fnames[ prefix])> 2 or sum([1 for elem in fnames[prefix] if elem is None]):
        raise ValueError("There must be just 2 pdb files for each complex to be predicted")
      else:
        fnames[ prefix]= tuple(fnames[ prefix])

    Parallel( n_jobs= ncpu, backend="multiprocessing", batch_size=2)(delayed(computeFunction)(
                            OneFeaturesComputerClass, fnames[ prefix][0], fnames[ prefix][1], computedFeatsRootDir, classArgs) 
                            for prefix in sorted(fnames))
      
  
def computeFunction(OneFeaturesComputerClass, rFname, lFname, computedFeatsRootDir, classArgs):
  '''
    Computes one type of features (the one specified by OneFeaturesComputerClass) for one complex whose pdbFiles are 
    the 2 first elements of the tuple classArgs.
    
    @param OneFeaturesComputerClass: FeaturesComputer. class to use for compute one kind of features
    @param rFname: str. path to receptor pdb file
    @param lFname: str. path to ligand pdb file      
    @param computedFeatsRootDir: str. path where features will be stored
    @param classArgs: dict. The arguments that will be passed to OneFeaturesComputerClass()
  '''
  pdbsPath, rSuffix= os.path.split(rFname)
  pdbsPath, lSuffix= os.path.split(lFname)  
  print(">>Computing:", str(OneFeaturesComputerClass), rSuffix, lSuffix)
  featCompObj= OneFeaturesComputerClass( rFname, lFname, computedFeatsRootDir, ** classArgs )
  featCompObj.computeFun()
  return 0
  
class FeatureComputerException(Exception):
  pass
  
def testModule():
  fnameL= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/1A2K_l_u.pdb"
  fnameR= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/1A2K_r_u.pdb"
  fetComp= FeaturesComputer(fnameR, fnameL)
  print(fetComp.__dict__)
  pdbsIndir= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/"
  computedFeatsRootDir= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/computedFeatures"
  FeaturesComputer.computeFeaturesAllComplexes(FeaturesComputer, pdbsIndir, computedFeatsRootDir)
  
if __name__=="__main__":
  testModule()

