from __future__ import absolute_import, print_function
import os, sys
import random
import joblib
import numpy as np
import pandas as pd

from utils import tryToRemove #utils is at the root of the package
from .CodifyComplexException import CodifyComplexException
from .ComplexCodified import ComplexCodified
from .codifyProtocols.SeqProtocol import SeqProtocol
from .codifyProtocols.StructProtocol import StructProtocol
from .codifyProtocols.MixedProtocol import MixedProtocol

class OneComplexCodifier(object):
  def __init__(self, dataRootPath, wholeComplexOutPath=None, environType="seq",
                     feedback_paths=None, sampledOutPath=None, samplingFold=2, verbose=False):
    '''
      @param dataRootPath: str. A path to computedFeatures directory that contains needed features. Example:
                computedFeatures/
                  common/
                    contactMaps/
                  seqStep/
                    conservation/
                    ...
                  structStep/
                    PSAIA/
                    VORONOI/
                    ...    
                  
      @param wholeComplexOutPath: str. A path where the complex codified (all positive and negative pairs) will be saved.
                                       If None, no file will be saved

      @param environType: str. "seq" if sequential enviroment protocol want to be used (sliding window of pssms...)
                               "struct" if VORONOI neighbours enviroment protocol want to be used (mean, min, max and
                               count for residues and their properties).

      @param feedback_paths: str or str[]. A path to previous results files directory. If it is None, contactMaps files will be used
                                 to define which residue pairs are in contact. Can also be a str[] if multiple feedback_paths
                                 wanted
                                
      @param sampledOutPath: str. Path where positive and sampled negative pairs will be saved.   
                                  If None, no file will be saved
        
      @param samplingFold: int.  Number of times the number of negative sampled pairs is bigger than positive pairs. 
                                 numNegativePairs= samplingFold*numPositivePairs. (dealing with imbalanced data sets)
    '''
    self.CodProtocol=None
    if environType.startswith("seq"):
      self.CodProtocol= SeqProtocol #Now self.CodProtocol is a class, later on self.CodProtocol will be an object
    elif environType.startswith("struct"):
      self.CodProtocol= StructProtocol #Now self.CodProtocol is a class, later on self.CodProtocol will be an object
    elif environType.startswith("mixed"):
      self.CodProtocol= MixedProtocol #Now self.CodProtocol is a class, later on self.CodProtocol will be an object
    else:
      raise CodifyComplexException("environType must startwith 'seq' or 'struct' or 'mixed'")
       
    if samplingFold<=0:
      raise CodifyComplexException("samplingFold must be greater than zero")
      
    self.verbose= verbose
    if self.verbose:
      print( ("Params: wholeComplexOutPath: %s, sampledOutPath: %s, dataRootPath: %s, feedback_paths: %s, "+\
                     "environType: %s, samplingFold: %s")%(  str(wholeComplexOutPath),  str(sampledOutPath),
                                    str(dataRootPath), feedback_paths, str(environType), str(samplingFold)))

    self.dataRootPath= dataRootPath                
    self.wholeComplexOutPath= wholeComplexOutPath
    self.environType= environType
    self.feedback_paths= feedback_paths    
    self.sampledOutPath= sampledOutPath
    self.samplingFold= samplingFold
    self.verbose= verbose
    self.cMapPath= None
    if self.feedback_paths== None:  # if no feedback path, then use contactMap info for interacting residues ground truth
      self.cMapPath= os.path.join( self.dataRootPath, "common", "contactMaps")
    elif isinstance(feedback_paths,list): #if several previous step results use last to get interacting residues ground truth
      self.cMapPath= os.path.realpath(os.path.expanduser(feedback_paths[-1]))
    else:
      self.cMapPath= os.path.realpath(os.path.expanduser(feedback_paths)) # use  previous step results to get interacting residues ground truth
    self.CodProtocol= self.CodProtocol(self.dataRootPath, self.cMapPath, self.feedback_paths, verbose=self.verbose) 

    
  def codifyComplex(self,prefix):
    '''
      Codifies one complex whose identifier is prefix. The features of the complex must have been computed previously
      and they must be located at self.dataRootPath path.
      @param prefix: str. A pdb id for a complex. If complex is formed by 2 pdb files, then prefix 
            looks like "receptorId<->ligandId" e.g. "1da2<->1lla"
      @return wholeComplexObject: ComplexCodified.ComplexCodified. A ComplexCodified object containing all putative
                                  pairs
    '''
    
    if self.verbose: print( "Codifying %s"%prefix )
    if "<->" in prefix:
      prefixR, prefixL= prefix.split("<->")
    else:
      prefixR, prefixL= (prefix,prefix)
    prefixR= prefix+"_r"
    prefixL= prefix+"_l"
    
    pairsCodifiedDir= self.CodProtocol.applyProtocol(prefix, prefixL, prefixR)
    wholeComplexObject= ComplexCodified(prefix, pairsCodifiedDir)
    if not self.sampledOutPath is None:
      if self.verbose: print("Sampling %s"%prefix)
      sampledComplexObject= wholeComplexObject.getSampledVersion(self.samplingFold)
      outName= os.path.join(self.sampledOutPath, prefix+".train.pkl.gz")      
      try:
        joblib.dump(sampledComplexObject, outName , compress= 5, protocol= 2)
      except (KeyboardInterrupt, Exception):
        print("Exception happened computing %s"%outName)
        tryToRemove(outName)
        raise
      
    if not self.wholeComplexOutPath is None:
      outName= os.path.join(self.wholeComplexOutPath, prefix+".predict.pkl.gz")
      try:
        if self.verbose: print("Writing results to disk")
        joblib.dump(wholeComplexObject, outName , compress= 5, protocol= 2)
      except (KeyboardInterrupt, Exception):
        print("Exception happened computing %s"%outName)
        tryToRemove(outName)    
        raise
    if self.verbose: print("%s succesfully codified"%(prefix))
##    a= raw_input("press enter to continue")
    return wholeComplexObject
    

