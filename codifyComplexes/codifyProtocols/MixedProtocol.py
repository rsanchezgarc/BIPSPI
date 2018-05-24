from __future__ import absolute_import, print_function
import os
import pandas as pd
import numpy as np
import copy
import time
from itertools import chain as itertoolchain 
from itertools import product as cartesianProd
from .StructProtocol import StructProtocol
from codifyComplexes.CodifyComplexException import CodifyComplexException
from computeFeatures.seqStep.seqToolManagers.conservationTools.windowPSSM import AA_CODE_ELEMENTS

#(feature_name, path_to_dir, columns, namedColsDict). If columns==None, all columns will be used
#Structural features must come first as sequential single chains features muy contain more aminoacids
# (non solved residues)
FEATURES_TO_INCLUDE_CHAIN= [
  ("psaia", ("structStep/PSAIA/procPSAIA", None, {"total_RASA":8})),
  ("halfSphereExpos", ("structStep/halfSphereExpos", None, {})),
  ("dssp", ("structStep/DSSP/procDSSP", [3],{})),
  ("al2co", ("seqStep/conservation/al2co",[4],{})),
  ("winSeqAndConservation", ("seqStep/conservation/pssms/windowedPSSMs/wSize11", None,{}))
]
FEATURES_TO_INCLUDE_PAIR= [
  ("corrMut", ("seqStep/conservation/corrMut", None, {})),
]

class MixedProtocol(StructProtocol):
  '''
    This class implements structural voronoi environment codification
  '''
  def __init__(self, dataRootPath, cMapPath, prevStepPaths=None, verbose=False):
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
                  
      @param cMapPath: str. A path to a dir that contains the contact map of the protein complex

      @param prevStepPaths: str or str[]. A path to previous results files directory. If it is None, contactMaps files will be used
                                 to define which residue pairs are in contact. Can also be a str[] if multiple feedback_path's
                                 wanted
                                
      @param environ_path: str. relative path (from dataRootPath) where Voronoi neighbours files are contained
    '''    
    StructProtocol.__init__(self, dataRootPath, cMapPath, prevStepPaths,
                                  singleChainfeatsToInclude=FEATURES_TO_INCLUDE_CHAIN, 
                                  pairfeatsToInclude= FEATURES_TO_INCLUDE_PAIR, verbose= verbose)

  def loadSingleChainFeatures(self, prefixOneChainType, chainType):
    '''
      @overrides StructProtocol method
      Loads all features files computed for ligand or receptor chains. Returns a pandas.DataFrame 
      that contains in each row all features from all files for each amino acid. Just amino acids
      that appears in each file will be included. Others will be ruled out (intersection)
      @param prefixOneChainType: str. A prefixOneChainType that identifies the receptor or ligand
      @param chainType: str. "l" for ligand and "r" for receptor
      @return df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      one amino acid
                      Column names are:
                      'chainId%s', 'structResId%s', 'resName%s', [properties] #no defined order for properties
                      %s is L if chainType=="l" and R if  chainType=="r"
    '''
    singleChainFeats= super(StructProtocol,self).loadSingleChainFeatures( prefixOneChainType, chainType)
        #add resName as feature
    chainType= chainType.upper()
    nAAsInWindow= 1+ len(set([ int(elem.split(".")[-1].split("_")[0]) 
                    for elem in singleChainFeats.columns if elem.startswith("resNameWin.")] ))
    centralRes= nAAsInWindow//2

    selectedSeqEntr= set([ 'seqEntropy.%d%s'%(i, chainType) for i in range(centralRes*4 , centralRes*4 +4)])
    selectedPssm= set([ 'pssm.%d%s'%(i, chainType) for i in range(centralRes*40,centralRes*40+40)])
    selectedWinAA= set([ 'resNameWin.%d_dummy_%s%s'%(centralRes,letter, chainType) for letter in AA_CODE_ELEMENTS])
    centralResCols= selectedSeqEntr.union(selectedPssm).union(selectedWinAA)
    winCols= set([col for col in singleChainFeats.columns if not "ggr" in col and(
                                                         "Win" in col or "pssm" in col or "seqEntropy" in col)])
    allWinButCentralCols= winCols.difference(centralResCols)

    allButWinData= singleChainFeats[ [col for col in singleChainFeats.columns if not col in allWinButCentralCols] ]
    winData= singleChainFeats[ list(singleChainFeats.columns[:3])+[col for col in singleChainFeats.columns if col in allWinButCentralCols] ]
    
    singleChainFeats= self.addSingleChainAggregation(allButWinData, chainType)
    mergeOn= [ elem%chainType.upper() for elem in ["chainId%s", "structResId%s", "resName%s"] ]
    singleChainFeats= pd.merge(singleChainFeats, winData, how='inner', on=mergeOn) 

    return singleChainFeats
    
    

