from __future__ import absolute_import, print_function
import os
import pandas as pd
import numpy as np
import time
from itertools import chain as itertoolchain 

from .AbstractProtocol import AbstractProtocol, computeNumericAggr, computeFactorAggr
from codifyComplexes.CodifyComplexException import CodifyComplexException
from computeFeatures.seqStep.seqToolManagers.conservationTools.windowPSSM import AA_CODE_ELEMENTS

#(feature_name, path_to_dir, columns, namedColsDict). If columns==None, all columns will be used
FEATURES_TO_INCLUDE_CHAIN= [
  ("winSeqAndConservation", ("seqStep/conservation/pssms/windowedPSSMs/wSize11", None,{})),
  ("al2co", ("seqStep/conservation/al2co",[4],{})),
  ("predAsaAndSS", ("seqStep/SPIDER2",None, {}))
]
FEATURES_TO_INCLUDE_PAIR= [
  ("corrMut", ("seqStep/conservation/corrMut", None, {})),
]

class SeqProtocol(AbstractProtocol):
  '''
    This class implements sequential environment codification (sliding window)
  '''
  IGNORE_FOR_AGGREGATION=["ccmPredQuality_P", "psicovQuality_P"]
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
    '''
    AbstractProtocol.__init__(self, dataRootPath, cMapPath, prevStepPaths, 
                                  singleChainfeatsToInclude=FEATURES_TO_INCLUDE_CHAIN, 
                                  pairfeatsToInclude= FEATURES_TO_INCLUDE_PAIR, verbose=verbose)
    
    
  def applyProtocol( self, prefixComplex, prefixL, prefixR):
    '''
      This method is the basic skeleton for applyProtocol of subclasses
      Given a prefix that identifies the complex and prefixes that identifies
      the ligand and the receptor, this method integrates the information that
      is contained in self.dataRootPath and is described in self.singleChainfeatsToInclude
      
      @param prefixComplex: str. A prefix that identifies a complex
      @param prefixL: str. A prefix that identifies the ligand of the complex
      @param prefixR: str. A prefix that identifies the receptor of the complex
      @return df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R).
                      Column names are:
                      'chainIdL', 'structResIdL', 'resNameL', 'chainIdR', 'structResIdR', 'resNameR', 'categ'
                       [propertiesP .... propertiesL     .... propertiesR] #no defined order for properties
    '''
    s=time.time()
    allPairsCodified= super(SeqProtocol,self).applyProtocol( prefixComplex, prefixL, prefixR)
    #adding seq length
    allPairsCodified= self.addSeqLen(allPairsCodified, "l")
    allPairsCodified= self.addSeqLen(allPairsCodified, "r")
    allPairsCodified= self.addPairwiseAggregation(allPairsCodified )
    allPairsCodified= self.reorderColumns(allPairsCodified)
    print("Time for %s codification:"%prefixComplex, time.time() -s)
    return allPairsCodified
     
  def addProductTerms(self, df):
    selectedPssmL= sorted([ 'pssm.%d%s'%(i, "L") for i in range(200,220)])
    selectedPssmR= sorted([ 'pssm.%d%s'%(i, "R") for i in range(200,220)])
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
      seqLenDict[chainId]= len(df['structResId%s'%chainType][
                                                df['chainId%s'%chainType]==chainId].unique())
    df['seqLen%s'%chainType]= -1* np.ones(df.shape[0])
    for chainId in seqLenDict:
      df.loc[df['chainId%s'%chainType]==chainId, 'seqLen%s'%chainType]= seqLenDict[chainId]   
    return df
    
    
  def _getSeqNeigs(self, dataFull, chainType, wsize=3):
  
    chainType= chainType.upper()
    dataIds=  dataFull.iloc[:,:6].reset_index()
    neigsRowsFromIds={}
    rowsFromId={}
    for chainId in dataIds["chainId%s"%chainType].unique():
      data= dataIds.loc[dataIds["chainId%s"%chainType]==chainId, :]
      allRowsByResInChain={}
      for i, resId in zip(data.index.values, data["structResId%s"%chainType].values):
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
      resturns featuresDicts to make it easier to compute aggregation of pairwise features
    
      @param df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R).
                      Column names are:
                      'chainIdL', 'structResIdL', 'resNameL', 'chainIdR', 'structResIdR', 'resNameR', 'categ' 
                      [properties_P .... propertiesL     .... propertiesR] #no defined order for properties
      @return df
                         
    '''
    df_ids= df.iloc[:, [0,1,3,4 ]]
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
 
