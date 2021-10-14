from __future__ import absolute_import, print_function
import os
import pandas as pd
import numpy as np
from itertools import chain as itertoolchain 

from .AbstractProtocol import AbstractProtocol, computeNumericAggr, computeFactorAggr
from codifyComplexes.CodifyComplexException import CodifyComplexException

#(feature_name, path_to_dir, columns). If columns==None, all columns will be used. colums=[12,13] Use columns 12 and 13
FEATURES_TO_INCLUDE_CHAIN= [
  ("winPssms", ("seqStep/conservation/pssms/windowedPSSMs/wSize11", None)),
  ("winSeq", ("seqStep/slidingWinSeq11", None)),
  ("al2co", ("seqStep/conservation/al2co",None)),
  ("predAsaAndSS", ("seqStep/SPIDER2_predAsaSS",None))
]
FEATURES_TO_INCLUDE_PAIR= [
  ("corrMut", ("seqStep/conservation/corrMut", None)),
]

class SeqProtocol(AbstractProtocol):
  '''
    This class implements sequential environment codification (sliding window)
  '''
  IGNORE_FOR_AGGREGATION=["ccmPredQuality_P", "psicovQuality_P"]
  def __init__(self, dataRootPath, cMapPath, prevStepPaths=None, singleChainfeatsToInclude=FEATURES_TO_INCLUDE_CHAIN, 
                     pairfeatsToInclude=FEATURES_TO_INCLUDE_PAIR, verbose=False):
    '''
      :param dataRootPath: str. A path to computedFeatures directory that contains needed features. It follows the
                                following standard:
                computedFeatures/
                  common/
                    contactMaps/
                  seqStep/
                    conservation/
                    ...
                  structStep/
                    PSAIA/
                    ...    
                  
      :param cMapPath: str. A path to a dir that contains the contact map of the protein complex

      :param prevStepPaths: str or str[]. A path to previous results files directory. If it is None, contactMaps files will be used
                                 to define which residue pairs are in contact. Can also be a str[] if multiple feedback_paths
                                 wanted
      :param verbose: bool.
    '''
    AbstractProtocol.__init__(self, dataRootPath, cMapPath, prevStepPaths, 
                                  singleChainfeatsToInclude= singleChainfeatsToInclude, 
                                  pairfeatsToInclude=pairfeatsToInclude, verbose=verbose)
    
    
  def applyProtocol( self, prefixComplex):
    '''
      This method is the basic skeleton for applyProtocol of subclasses
      Given a prefix that identifies the complex and prefixes that identifies
      the ligand and the receptor, this method integrates the information that
      is contained in self.dataRootPath and is described in self.singleChainfeatsToInclude
      
      :param prefixComplex: str. A prefix that identifies a complex e.g. 1ACB
      :return df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R).
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ'
                       [propertiesP .... propertiesL     .... propertiesR] #no defined order for properties
    '''

    allPairsCodified= AbstractProtocol.applyProtocol(self,  prefixComplex)
    #adding seq length
    allPairsCodified= self.addSeqLen(allPairsCodified, "l")
    allPairsCodified= self.addSeqLen(allPairsCodified, "r")
    allPairsCodified= self.addPairwiseAggregation(allPairsCodified )
    allPairsCodified= self.reorderColumns(allPairsCodified)

    return allPairsCodified
     
  def addProductTerms(self, df):
    winSize= max([ int(elem.split(".")[-1][:-1]) for elem in df.columns if elem.startswith("pssmWin") ])+1
    selectedPssmL= sorted([ 'pssmWin.%d.%d%s'%(i, winSize//2, "L") for i in range(20)])
    selectedPssmR= sorted([ 'pssmWin.%d.%d%s'%(i, winSize//2, "R") for i in range(20)])
    for colL in selectedPssmL:
      for colR in selectedPssmR:
        df[ colL+colR+"_P"]= df[colL]*df[colR]
    return df
    
  def addSeqLen(self, df, chainType):
    '''
    '''
    chainType= chainType.upper()
    chainsList= df['chainId%s'%chainType].unique()
    seqLenDict={}
    for chainId in chainsList:
      seqLenDict[chainId]= len(df['resId%s'%chainType][
                                                df['chainId%s'%chainType]==chainId].unique())
    df['seqLen%s'%chainType]= -1* np.ones(df.shape[0])
    for chainId in seqLenDict:
      df.loc[df['chainId%s'%chainType]==chainId, 'seqLen%s'%chainType]= seqLenDict[chainId]   
    return df
    
    
  def _getSeqNeigs(self, dataFull, chainType, wsize=3):
  
    chainType= chainType.upper()
    categIndex= list(dataFull.columns).index("categ")
    dataIds=  dataFull.iloc[:, :(categIndex+1)].reset_index()
    neigsRowsFromIds={}
    rowsFromId={}
    for chainId in dataIds["chainId%s"%chainType].unique():
      data= dataIds.loc[dataIds["chainId%s"%chainType]==chainId, :]
      allRowsByResInChain={}
      for i, resId in zip(data.index.values, data["resId%s"%chainType].values):
        resId= (int(resId),"") if resId[-1].isdigit() else (int(resId[:-1]), resId[-1])
        if not resId in allRowsByResInChain:
          allRowsByResInChain[resId]=[]
        allRowsByResInChain[resId].append(i)

      resIdsNames= sorted(allRowsByResInChain)
      nResidues= len(resIdsNames)
      
      for i, resId in enumerate(resIdsNames):
        neigsRowsForRes= set(itertoolchain.from_iterable( [allRowsByResInChain[resIdsNames[resIx]] 
                          for resIx in range(i-wsize//2, i+wsize//2 +1) if resIx>=0 and resIx< nResidues and resIx!= i]))           
        neigsRowsFromIds[(chainId, resId)]= neigsRowsForRes
        rowsFromId[(chainId, resId)]= set(allRowsByResInChain[resId])
    return rowsFromId, neigsRowsFromIds
   
  def prepareDataForPairwiseAggregat(self, df):
    '''
      resturns different data to make it easier to compute aggregation of pairwise features
    
      :param df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R).
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ' 
                      [properties_P .... propertiesL     .... propertiesR] #no defined order for properties

      :return  pairwiseDf, ids2RowL, ids2RowR, neigsids2rowL, neigsids2rowR
          pairwiseDf: the Pandas.DataFrame of features after aggregation
          ids2RowL: a dict that maps from resIdL to all rows of df that contain resIdL
          ids2RowL: a dict that maps from resIdR to all rows of df that contain resIdR
          neigsids2rowL: a dict that maps from resIdL to all rows of df that contain a neighbour of resIdL
          neigsids2rowR: a dict that maps from resIdR to all rows of df that contain a neighbour of resIdR
    '''
#    categIndex= list(df.columns).index("categ")
#    df_ids= df.iloc[:, :(categIndex+1) ]
    pairwiseDf= df[ [elem for elem in df.columns if elem.endswith("P") and elem not in SeqProtocol.IGNORE_FOR_AGGREGATION] ]
    ids2RowL, neigsids2rowL= self._getSeqNeigs( df, chainType="l")
    ids2RowR, neigsids2rowR= self._getSeqNeigs( df, chainType="r")
        
    return pairwiseDf, ids2RowL, ids2RowR, neigsids2rowL, neigsids2rowR  

  def addPairwiseAggregation(self, dataFull):
    '''
      Overrides AbstractProtocol.addPairwiseAggregation
    '''
    varlist= [ elem for elem in dataFull.columns if elem.endswith("_P")]
    if len(varlist)==0:
      return dataFull      
    
    return AbstractProtocol.addPairwiseAggregation(self, dataFull)
 
