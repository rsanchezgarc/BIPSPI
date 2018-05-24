'''
Codify all benchmark complexes from previously computed features files
(possibly scores from previous steps)
For features description and file requirements see codifyComplexes.codifyProtocols
'''
from __future__ import absolute_import, print_function

import os, sys, shutil
import multiprocessing
from joblib import Parallel, delayed  

from utils import myMakeDir, myMakeDirUnique, tryToRemove #utils is at the root of the package
from .codifyOneComplex import OneComplexCodifier
from .CodifyComplexException import CodifyComplexException
from Config import Configuration

class BenchmarkCodificator(Configuration):
  '''
    Class that allows for all benchmark complexes codification
  '''
  def __init__(self, features_path=None, out_Codified_path=None, feedback_paths= None, environType=None, ncpu=1, 
                      overridePrevComp= False, verbose=False):
    '''
      @param features_path: str. A path to the computedFeatures directory that contains needed features.
                            If None, Config.py DEFAULT_PARAMETERS["computedFeatsRootDir"] will be used
                            Example:
                            features_path/
                              common/
                                contactMaps/
                              seqStep/
                                conservation/
                                ...
                              structStep/
                                PSAIA/
                                VORONOI/
                                 ...

      @param out_Codified_path: str. Root directory where codified complexes will be saved. Files will be saved at directory
                            out_Codified_path/seq[_n] if sequential enviroment protocol will be used or at
                            out_Codified_path/struct[_n] if structural enviroment protocol will be used. (attribute self.environType)
                            If more than one step of the same type, the path will end with "_#", p.e:
                              1) path/to/outpath/seq/
                              2) path/to/outpath/struct_0/
                              3) path/to/outpath/struct_1/
                            If None, Config.py DEFAULT_PARAMETERS["codifiedDataRootDir"] will be used as out_Codified_path                                                

      @param feedback_paths: str or str[]. A path to a previous results files directory. Contact maps for evaluation will 
                                          be obtained from this file. If None, contact maps will be loaded from contactMaps
                                          files contained at features_path/common/contactMaps/
                                 
      @param environType: str. "seq" if sequential enviroment protocol want to be used (sliding window of pssms...)
                               "struct" if VORONOI neighbours enviroment protocol want to be used (mean, min, max, sum and
                                count for neighbour residues and their properties).

      @param ncpu: int. Number of processes to use in parallel (each process will process one complex)
      
      @param overridePrevComp: boolean. If True and there are complexes at out_Codified_path, those complexes will be overrided.
                             If False, already computed complexes will be kept and codification will continue
                             with non computed complexes 
          
    '''
    Configuration.__init__(self)
    if not (environType.startswith("seq") or environType.startswith("struct") or environType.startswith("mixed")):
      raise CodifyComplexException("environType must be 'seq' or 'struct' or 'mixed'")
      
    self.environType= environType
    if features_path is None:
      features_path= self.computedFeatsRootDir
    self.dataRootPath= os.path.realpath(os.path.expanduser( features_path))
    if out_Codified_path is None:
      out_Codified_path= self.codifiedDataRootDir
    self.out_Codified_path= myMakeDir( os.path.realpath(os.path.expanduser( out_Codified_path)) )
    self.feedback_paths=  feedback_paths  #Either a path or None
    self.overridePrevComp= overridePrevComp
    self.verbose= verbose

    self.out_Codified_path=  myMakeDir(self.out_Codified_path, environType)
    
    self.ncpu= ncpu   
    self.testingDataPath= os.path.join( self.out_Codified_path,"allInputs")
    self.trainingDataPath=os.path.join(self.out_Codified_path,"sampledInputs")
    
    if self.ncpu>multiprocessing.cpu_count() or self.ncpu==-1:
      self.ncpu=multiprocessing.cpu_count()
    elif self.ncpu<1:
      self.ncpu= 1

  def inspectOldDir(self, dirname):
    '''
      check if dirname contain files that look like previously codified complexes and:
        deletes them if self.overridePrevComp==True  or 
        continue execution
      if some filename does not match expected extension, a  CodifyComplexException will be raised 
      @param dirname: str. A path to codified complexes

    '''

    areNotAllPickles= sum( (1 for fname in os.listdir(dirname) if not fname.endswith(".pkl.gz") ) ) #Sanity check for wrong paths
    if areNotAllPickles>0:
      errorMsg="There are files that does not match rriPredMethod format. Check them out and consider to use"\
            +":\nrm -r %s/*"%( dirname ) 
      print(errorMsg+"\naborting")
      raise CodifyComplexException(errorMsg)    
    if self.overridePrevComp: #remove  previously codified files
      print("Deleting all previously codified files at %s"%( dirname ))
      for fname in os.listdir(dirname):
        if fname.endswith(".pkl.gz"):
          os.remove( os.path.join(dirname, fname))
    else:
      print("Resuming benchmark codification at %s"%dirname)
      
  def setResultsDirs(self):
    '''
      check if directories where codified complexes will be saved already exists. If don't create 
      them. If they do exist, clean if self.overridePrevComp==True otherwise continue execution
      continue execution
    '''  
    testingDataPath= self.testingDataPath
    trainingDataPath= self.trainingDataPath
    
    print("Coded input will be save at", self.out_Codified_path)

    if not os.path.exists( testingDataPath):
      os.mkdir(testingDataPath)
    else:
      self.inspectOldDir( testingDataPath)
      
    if not os.path.exists( trainingDataPath):
      os.mkdir(trainingDataPath)
    else:
      self.inspectOldDir( trainingDataPath)          
          
  def codifyAll(self, skipComplexesList=[], samplingFold=2):
    '''
      Codify all complexes for which there are perviousResultsFile or contactMap files.
                
      @param skipComplexesList: str[]. A list of complexes (prefixes) that are included in features_path but are not wanted 
                                         to be codified. Prefixes will be define as follows: 
                                            1A2K_l_r.u.pdb --> 1A2K

      @param samplingFold: int.  Number of times the number of negative sampled pairs is bigger than positive pairs
                                 numNegativePairs= samplingFold *numPositivePairs
                                 (Dealing with imbalanced data sets)
      @return self.out_Codified_path: str. The path where complexes where codified
    '''

    self.setResultsDirs()
    
    if self.feedback_paths== None:  # if no feedback path, then use contactMap info for interacting residues ground truth
      cMapPath= os.path.join( self.dataRootPath, "common", "contactMaps")
    elif isinstance(self.feedback_paths,list): #if several previous step results use last to get interacting residues ground truth
      cMapPath= os.path.realpath(os.path.expanduser(self.feedback_paths[-1]))
    else:
      cMapPath= os.path.realpath(os.path.expanduser(self.feedback_paths)) # use  previous step results to get interacting residues ground truth

    skipComplexesSet= set(skipComplexesList)
    argsToCodifyOneComplex=[]
    totalNum=0
    alreadyComputedNum=0
    for fName in sorted( os.listdir( cMapPath)):
      if fName.endswith(".res.tab") or fName.endswith(".cMap.tab"):
        prefix=fName.split(".")[0]
        if prefix in skipComplexesSet:
#          print("Skipping %s"%prefix )
          continue
        totalNum+= 1
        saveNameTrain=    os.path.join( self.testingDataPath,  prefix+".predict.pkl.gz")
        saveNameTest= os.path.join( self.trainingDataPath, prefix+".train.pkl.gz")
        if (os.path.isfile(saveNameTrain) and os.path.isfile(saveNameTest) and
            os.stat(saveNameTrain).st_size > 0 and os.stat(saveNameTest).st_size > 0):
          alreadyComputedNum+=1
          if self.verbose: print("%s already computed"%prefix)
          continue  
        argsToCodifyOneComplex.append((prefix, self.dataRootPath, self.testingDataPath, self.environType, 
                                        self.feedback_paths, self.trainingDataPath, samplingFold, self.verbose))
    nComplexes= len(argsToCodifyOneComplex)
    if nComplexes==0 and self.overridePrevComp:
      raise CodifyComplexException("No complexes to compute. There are no a single cMap or prevStep result in %s"%cMapPath)
    print ("Number of remaining complexes to codify: %d. Already computed: %d"%(nComplexes, alreadyComputedNum) )
    #Parallel execution
    checksum=Parallel(n_jobs=self.ncpu)(delayed(launchCodifyOneComplex)(inp) for inp in argsToCodifyOneComplex) 
    return self.out_Codified_path
              
def launchCodifyOneComplex( funArgs ):
  '''
    Launches codification of onw complex. This method was declared at global scope to be used with Parallel
              
    @param funArgs: Tuple
                         funArgs[0] : str Prefix. Prefixes will be define as follows: 
                                            1A2K_l_r.u.pdb --> 1A2K
                         funArgs[1...]: Tuple. Arguments to OneComplexCodifier 
                                                      ( dataRootPath, wholeComplexOutPath=None, environType,
                                                         feedback_paths, sampledOutPath, samplingFold=)
             
  '''
  prefix= funArgs[0]
  builderArgs= funArgs[1:]
  print ("Working on complex %s"%funArgs[0])
  OneComplexCodifier( *builderArgs).codifyComplex(prefix)
  return 0              
    
if __name__=="__main__":

  print ("Running default")
  features_path= "~/Tesis/rriPredMethod/data/ppdockingBenchData/tmpComputedFeatures"
#  features_path= "~/Tesis/rriPredMethod/data/ppdockingBenchData/computedFeatures"
  out_Codified_pathRoot= "~/Tesis/rriPredMethod/data/ppdockingBenchData/tmpCodifiedInput"
  feedback_paths= None

  prefixes= [ elem.split(".")[0] for elem in os.listdir(os.path.join(os.path.expanduser(features_path),"common","contactMaps"))]
  prefixes.sort()

  
  skipComplexesList= [ "2PCC"]
  skipComplexesList= prefixes[4:]
#  skipComplexesList= [elem for elem in prefixes if elem== "1A2K"]

  samplingFold=2
  ncpu= 1
  overridePrevComp= False
  environType="struct"
  benchCod= BenchmarkCodificator( features_path, out_Codified_pathRoot, feedback_paths, environType,ncpu, overridePrevComp)
  benchCod.codifyAll( skipComplexesList=skipComplexesList, samplingFold= samplingFold)                                   

