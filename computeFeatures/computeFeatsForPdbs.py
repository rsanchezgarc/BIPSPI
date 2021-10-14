'''
This module is used for training feature calculations
'''
from __future__ import absolute_import, print_function
import sys, os
from collections import OrderedDict

from multiprocessing import cpu_count
from Config import Configuration
from joblib import Parallel, delayed
from .featuresComputer import  getExtendedPrefix, splitExtendedPrefix, FeatureComputerException
from .computeFeatsOneComplex import OneComplexFeatComputer



def computeFeaturesAllPdbsOneDir(pdbsIndir= None, computedFeatsRootDir= None,
                                 methodProtocol="struct", isHomeSet=False, ncpu=2):
  '''
    Computes all features needed for complex codification for all complexes in pdbsIndir. Used for training
    :param pdbsIndir: str. Path to the directory where pdb files are located. Must be named as follows:
                    path/to/pdbsIndir/
                                      prefix1_[chainType]_u.pdb or pdb.gz
                                      prefix2_[chainType]_u.pdb or pdb
                            By default, it uses as pdbsIndir Config.py DEFAULT_PARAMETERS["pdbsIndir"] 
    :param computedFeatsRootDir: str. Path where features files will be saved. By default it uses
                                Config.py DEFAULT_PARAMETERS["computedFeatsRootDir"] will be used as out_path
    :param methodProtocol: str. "seq" if just sequential features will be used; "struct" if sequential and
                                structural features will be used. "mixed" behaves as "struct"
    :param isHomeSet: True if homo-dataset is to be used, False for hetero dataset or unknown
    :param ncpu: int. Number of cpu's to use. If -1, all cpu's will be used to parallelize at complex level
  '''
  assert methodProtocol in ["seq", "struct", "mixed"], "Error methodProtocol in computeFeaturesAllPdbsOneDir must " + \
                                                   "be 'seq' or 'struct' or 'mixed'->"+str(methodProtocol)

  if pdbsIndir is None or computedFeatsRootDir is None:
    # Default parameters
    conf = Configuration()
    pdbsIndir= pdbsIndir if pdbsIndir else conf.pdbsIndir
    computedFeatsRootDir = computedFeatsRootDir if computedFeatsRootDir else conf.computedFeatsRootDir


  allCodifiedFname= os.path.join(computedFeatsRootDir, "allFeaturesComputed.txt")
  if os.path.exists(allCodifiedFname):
    with open(allCodifiedFname) as f:
      print(f.read())
    return

  if ncpu<1:
    ncpu= cpu_count()

  if pdbsIndir== None:
    pdbsIndir= conf.pdbsIndir
  if computedFeatsRootDir== None:
    computedFeatsRootDir= conf.computedFeatsRootDir
  pdbsIndir= os.path.expanduser(pdbsIndir)
  computedFeatsRootDir= os.path.expanduser(computedFeatsRootDir)
  fnames= OrderedDict({})
  fnamesOther= OrderedDict({})
  for fname in sorted(os.listdir(pdbsIndir)):
    if fname.endswith("_u.pdb.gz") or fname.endswith("_u.pdb"):  #skip no pdb files
      prefix, chainType= splitExtendedPrefix(getExtendedPrefix(fname, splitTag="_u."))
      if prefix not in fnames:
        fnames[prefix]=  [None, None]
      if chainType=="r":
        fnames[ prefix][1]= os.path.join(pdbsIndir, fname)
      elif chainType=="l":
        fnames[ prefix][0]= os.path.join(pdbsIndir, fname)
      else:
        raise FeatureComputerException("Error in filename %s"%fname)
    elif fname.endswith("_b.pdb.gz") or fname.endswith("_b.pdb"):
      prefix, chainType= splitExtendedPrefix(getExtendedPrefix(fname, splitTag="_b."))
      if prefix not in fnamesOther:
        fnamesOther[prefix]=  [None, None]
      if chainType=="r":
        fnamesOther[ prefix][1]= os.path.join(pdbsIndir, fname)
      elif chainType=="l":
        fnamesOther[ prefix][0]= os.path.join(pdbsIndir, fname)
      else:
        raise FeatureComputerException("Error in filename %s"%fname)

  if len(fnames) ==len(fnamesOther):
    boundAvailable=True
  else:
    boundAvailable=False
       
  if len(fnames)<1: raise ValueError("There are not files to be processed")
  for prefix in fnames: # check for errors
    if sum([1 for elem in fnames[prefix] if elem is None])>0:
      print(fnames[ prefix])
      raise ValueError("There must be just 2 pdb files for each complex to be processed")
    if boundAvailable:
      if sum([1 for elem in fnamesOther[prefix] if elem is None])>0:
        raise ValueError("There is no bound structure for some of your pdb files")
    
  print("Bound available:", boundAvailable)
  Parallel( n_jobs= ncpu, backend="multiprocessing", batch_size=2)(delayed(launchComputeFeaturesOneComplex)(
                              fnames[ prefix], prefix, computedFeatsRootDir= computedFeatsRootDir, 
                              boundAvailable= boundAvailable, methodProtocol= methodProtocol,
                              checkIfLRHomo=isHomeSet)
                          for prefix in sorted(fnames))

  with open(allCodifiedFname, "w") as f:
    f.write("All features computed for: %d"%len(fnames) )

def launchComputeFeaturesOneComplex(ligAndRecFnames, prefix, computedFeatsRootDir, boundAvailable, methodProtocol,
                                    checkIfLRHomo ):
  '''
    :param ligAndRecFnames: [pathToLigandPdb, pathToReceptorPDB]
    :param: boundAvailable. True if there is a bound and unbound pdb for each complex. False otherwise
  '''    
  featComp= OneComplexFeatComputer(prefix, computedFeatsRootDir, methodProtocol, areForTrainAndTest=True, 
                                   boundAvailable= boundAvailable, statusManager= None)
  featComp.computeFeaturesOneComplex( * ligAndRecFnames, **{"isHomoComplex":checkIfLRHomo})
  
  
def test():
  '''
  python -m computeFeatures.computeFeatsForPdbs
  '''
  computeFeaturesAllPdbsOneDir(ncpu=1)
  pass
    
if __name__ == "__main__":
  '''
  python -m computeFeatures.computeFeatsForPdbs
  '''
  test(); sys.exit(0)
  proto="mixed"
  if len(sys.argv)>1:
    pdbsIndir= sys.argv[1]
  if len(sys.argv)>2:
    computedFeatsRootDir= sys.argv[2]
  if len(sys.argv)>3:
    ncpu= int(sys.argv[3])
    
  computeFeaturesAllPdbsOneDir(pdbsIndir, computedFeatsRootDir, methodProtocol=proto, ncpu=ncpu)
  

