from __future__ import absolute_import, print_function
import os
import pandas as pd
import numpy as np

from .BaseStructProtocol import BaseStructProtocol
from codifyComplexes.CodifyComplexException import CodifyComplexException
from computeFeatures.seqStep.seqToolManager import SeqToolManager
AA_CODE_ELEMENTS= SeqToolManager.AA_CODE_ELEMENTS
'''
(feature_name, path_to_dir, columns ). If columns==None, all columns will be used
Structural features must come first as sequential single chains features muy contain more aminoacids
(e.g. non 3D-solved residues)
'''
FEATURES_TO_INCLUDE_CHAIN= [
  ("psaia", ("structStep/PSAIA/procPSAIA", None)),
  ("halfSphereExpos", ("structStep/halfSphereExpos", None)),
  ("dssp", ("structStep/DSSP", [3])),
  ("al2co", ("seqStep/conservation/al2co",None)),
  ("winPssms", ("seqStep/conservation/pssms/windowedPSSMs/wSize11", None)),
  ("winSeq", ("seqStep/slidingWinSeq11", None))
]
FEATURES_TO_INCLUDE_PAIR= [
  ("corrMut", ("seqStep/conservation/corrMut", None)),
]

class StructProtocol(BaseStructProtocol):
  '''
    This class implements structural voronoi environment codification
  '''
  def __init__(self, dataRootPath, cMapPath, prevStepPaths=None, singleChainfeatsToInclude= FEATURES_TO_INCLUDE_CHAIN, 
                     pairfeatsToInclude=FEATURES_TO_INCLUDE_PAIR, verbose=False):
    '''
      :param dataRootPath: str. A path to computedFeatures directory that contains needed features. Example:
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
                  
      :param cMapPath: str. A path to a dir that contains the contact map of the protein complex

      :param prevStepPaths: str or str[]. A path to previous results files directory. If it is None, contactMaps files will be used
                                 to define which residue pairs are in contact. Can also be a str[] if multiple feedback_path's
                                 wanted
    '''
    BaseStructProtocol.__init__(self, dataRootPath, cMapPath, prevStepPaths,
                                  singleChainfeatsToInclude=FEATURES_TO_INCLUDE_CHAIN, 
                                  pairfeatsToInclude= FEATURES_TO_INCLUDE_PAIR, verbose= verbose)

  def loadSingleChainFeatures(self, prefixOneChainType, chainType):
    '''
      @overrides BaseStructProtocol method to make use of sequence profiles (loaded directly) and struct
                  neighbour but not computing struct neighbours on non central residue features of sliding window
      Loads all features files computed for ligand or receptor chains. Returns a pandas.DataFrame 
      that contains in each row all features from all files for each amino acid. Just amino acids
      that appears in each file will be included. Others will be ruled out (intersection)
      :param prefixOneChainType: str. A prefixOneChainType that identifies the receptor or ligand
      :param chainType: str. "l" for ligand and "r" for receptor
      :return df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      one amino acid
                      Column names are:
                      'chainId%s', 'resId%s', 'resName%s', [properties] #no defined order for properties
                      %s is L if chainType=="l" and R if  chainType=="r"
    '''
    #super (BaseStructProtocol,self) is AbstractProtocol
    singleChainFeats= super(BaseStructProtocol,self).loadSingleChainFeatures( prefixOneChainType, chainType) #Load with no aggregation

    chainType= chainType.upper()
    winSize= max([ int(elem.split(".")[-1][:-1]) for elem in singleChainFeats.columns if elem.startswith("pssmWin") ])+1
    centralRes= winSize//2
    #find variables that will not be considered for structural aggreation: sliding window features of non central amino acids
  
    selectedSeqEntr= set([ 'informationWin.%d.%d%s'%(i, centralRes, chainType) for i in range(2)])
    selectedPssm= set([    'pssmWin.%d.%d%s'%(i, centralRes, chainType) for i in range(20)])
    selectedPsfm= set([    'psfmWin.%d.%d%s'%(i, centralRes, chainType) for i in range(20)])
    selectedWinAA= set([   'aaWin.0.%d_dummy_%s%s'%(centralRes,letter, chainType) for letter in AA_CODE_ELEMENTS ])

    #this variables will be aggregated
    centralResCols= selectedSeqEntr.union(selectedPssm).union(selectedPsfm).union(selectedWinAA)

    winCols= set([col for col in singleChainFeats.columns if not "ggr" in col and "Win" in col ])
    #this variables will not be aggreaged
    allWinButCentralCols= winCols.difference(centralResCols)
    
    allButWinData= singleChainFeats[ [col for col in singleChainFeats.columns if not col in allWinButCentralCols] ]
    winData= singleChainFeats[ list(singleChainFeats.columns[:3])+[col for col in singleChainFeats.columns if col in allWinButCentralCols] ]
#    print( list( allButWinData.columns) );raw_input("enter")
    singleChainFeats= self.addSingleChainAggregation(allButWinData, chainType)
    mergeOn= [ elem%chainType.upper() for elem in ["chainId%s", "resId%s", "resName%s"] ]
    singleChainFeats= pd.merge(singleChainFeats, winData, how='inner', on=mergeOn) 
    return singleChainFeats


# uncomment to use product terms
#  def addProductTerms(self, df):

#    winSize= max([ int(elem.split(".")[-1][:-1]) for elem in df.columns if elem.startswith("pssmWin") ])+1
#    centralRes= winSize//2

#    selectedColsL= sorted(['pssmWin.%d.%d%s'%(i, centralRes, "L") for i in range(20)] +
#                          [ 'total_RASAL', 'average_DPXL', 'HydrophobicityL'])
#    selectedColsR= sorted(['pssmWin.%d.%d%s'%(i, centralRes, "R") for i in range(20)] +
#                          [ 'total_RASAR', 'average_DPXR', 'HydrophobicityR'])
#    print(selectedColsL)
#    for colL in selectedColsL:
#      for colR in selectedColsR:
#        df[ colL+colR+"_P"]= df[colL]*df[colR]
#    return df

#  def applyProtocol( self, prefixComplex, prefixL, prefixR):
#    '''
#      This method is the basic skeleton for applyProtocol of subclasses
#      Given a prefix that identifies the complex and prefixes that identifies
#      the ligand and the receptor, this method integrates the information that
#      is contained in self.dataRootPath and is described in self.singleChainfeatsToInclude
#      
#      :param prefixComplex: str. A prefix that identifies a complex
#      :param prefixL: str. A prefix that identifies the ligand of the complex
#      :param prefixR: str. A prefix that identifies the receptor of the complex
#      :return df: pandas.DataFrame. A pandas.Dataframe in which each row represents
#                      a pair of amino acids in direct form (L to R).
#                      Column names are:
#                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ'
#                       [ propertiesL .... propertiesR .... propertiesP]
#    '''
#    allPairsCodified= super(StructProtocol,self).applyProtocol( prefixComplex, prefixL, prefixR)
#    allPairsCodified= self.addProductTerms(allPairsCodified)
#    allPairsCodified= self.reorderColumns(allPairsCodified)
#    return allPairsCodified

