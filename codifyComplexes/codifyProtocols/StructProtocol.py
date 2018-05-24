from __future__ import absolute_import, print_function
import os
import pandas as pd
import numpy as np
import copy
import time
from itertools import chain as itertoolchain 
from itertools import product as cartesianProd
from .PairwiseAAPots import PairwiseAAindex
from .AbstractProtocol import AbstractProtocol, computeNumericAggr, computeFactorAggr
from .NeigbManager import NeigbManager
from codifyComplexes.CodifyComplexException import CodifyComplexException
from computeFeatures.seqStep.seqToolManagers.conservationTools.windowPSSM import AA_CODE_ELEMENTS


#(feature_name, path_to_dir, columns, namedColsDict). If columns==None, all columns will be used
#Structural features must come first as sequential single chains features muy contain more aminoacids
# (non solved residues)
FEATURES_TO_INCLUDE_CHAIN= [
  ("psaia", ("structStep/PSAIA/procPSAIA", None, {"total_RASA":8})),
  ("halfSphereExpos", ("structStep/halfSphereExpos", None, {})),
  ("dssp", ("structStep/DSSP/procDSSP", [3],{})),
  ("psiBlastPSSMAndSeqEntropy", ("seqStep/conservation/pssms/procPssms", list(range(4,4+20))+[44,45],{})),
  ("al2co", ("seqStep/conservation/al2co",[4],{})),
]

FEATURES_TO_INCLUDE_PAIR= [
  ("corrMut", ("seqStep/conservation/corrMut", None, {})),
]


class StructProtocol(AbstractProtocol, NeigbManager):
  '''
    This class implements structural voronoi environment codification
  '''
  IGNORE_FOR_AGGREGATION=["ccmPredQuality_P", "psicovQuality_P"]
  def __init__(self, dataRootPath, cMapPath, prevStepPaths=None, verbose=False, 
                    singleChainfeatsToInclude= FEATURES_TO_INCLUDE_CHAIN,
                    pairfeatsToInclude=FEATURES_TO_INCLUDE_PAIR):
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
    
    NeigbManager.__init__(self, dataRootPath, singleChainfeatsToInclude) #This must be applied first to not override pairfeatsToInclude
    AbstractProtocol.__init__(self, dataRootPath, cMapPath, prevStepPaths, singleChainfeatsToInclude=singleChainfeatsToInclude,
                                    pairfeatsToInclude= pairfeatsToInclude, verbose= verbose)

#    self.aaIndexManager= PairwiseAAindex()

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
                       [propertiesP .... propertiesL     .... propertiesR]
    '''
    s=time.time()
    self.loadEnvirons(prefixL, prefixR)
    allPairsCodified= super(StructProtocol,self).applyProtocol( prefixComplex, prefixL, prefixR)
    allPairsCodified= self.addPairwiseAggregation(allPairsCodified)
    allPairsCodified= self.reorderColumns(allPairsCodified)    
    print("Time for %s codification:"%prefixComplex, time.time() -s)
    return allPairsCodified

  def addProductTerms(self, df):
    selectedColsL= sorted(['pssm%s'%("L")]+ [ 'pssm.%d%s'%(i, "L") for i in range(1,20)] +
                          ['score_xL', 'total_RASAL', 'average_DPXL', 'HydrophobicityL'])
    selectedColsR= sorted(['pssm%s'%("R")]+ [ 'pssm.%d%s'%(i, "R") for i in range(1,20)]+
                          ['score_xR', 'total_RASAR', 'average_DPXR', 'HydrophobicityR'])
    for colL in selectedColsL:
      for colR in selectedColsR:
        df[ colL+colR+"_P"]= df[colL]*df[colR]
    return df
    
  def loadSingleChainFeatures(self, prefixOneChainType, chainType):
    '''
      @overrides AbstractProtocol method
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
    aaLevels= AA_CODE_ELEMENTS
    chainType= chainType.upper()

    singleChainFeats["resName%s"%chainType]=pd.Categorical(singleChainFeats["resName%s"%chainType], 
                                                          categories= aaLevels, ordered=False) 
    singleChainFeats= pd.concat([singleChainFeats, pd.get_dummies(singleChainFeats[["resName%s"%chainType]], 
                                            prefix_sep='_dummy_', columns= ["resName%s"%chainType])], axis=1)
    newColNames= []
    for colname in singleChainFeats.columns:
      if "resName%s_dummy"%chainType in colname:
        newColNames.append( colname+chainType) 
      else:
        newColNames.append(colname)
    singleChainFeats.columns= newColNames

    singleChainFeats= self.addSingleChainAggregation(singleChainFeats, chainType)

    return singleChainFeats
    
    
  def addSingleChainAggregation(self, df, chainType):
    '''
      Adds environment single chain features to a df that represent single chain features.
      @param df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      an amino acid of ligand or receptor.
                      Column names are:
                      'chainId%s', 'structResId%s', 'resName%s',
                      [properties_P .... propertiesL .... propertiesR] #no defined order for properties
      @param chainType: str. "l" for ligand and "r" for receptor
      
      @return newDf: pandas.DataFrame. A pandas.Dataframe in which each row represents an amino acid of 
                      ligand or receptor. New columns will be added, all then
                      named "%sAggr%d%s"%(num/factor, numberOfNewFeature, ligand/receptor)
                      Column names are:
                      'chainId%s', 'structResId%s', 'resName%s',
                      [ properties%s ] #no defined order for properties    
                      
                         
    '''
    chainType= chainType.upper()
    factorCols=[]
    numericCols=[]
    for i,colName in enumerate(df.columns[3:]): #2: To skip chainId, resId, resName
      if "_dummy_" in colName:
        factorCols.append(3+i)
      else:
        numericCols.append(3+i)

    numericData= df.iloc[:, numericCols]
    factorData=  df.iloc[:, factorCols ]

    

    ids= zip( df.iloc[:, 0].values, df.iloc[:, 1].values)
    name2Id={ fullId: i for i, fullId in enumerate(ids)}
    id2Name={ i: fullId for fullId,i in name2Id.iteritems()}

    involvedRows= [[] for row in range(numericData.shape[0]) ]
    for fullId in name2Id:
      ix1= name2Id[fullId]
      neigsSet= self.getNeigs(fullId[0],fullId[1],chainType)
#      involvedRows[ix1].append( ix1)
      if not neigsSet is None:
        for neig in neigsSet:
          if neig in name2Id:
            involvedRows[ix1].append( name2Id[neig])

    numericAggre= computeNumericAggr(numericData, involvedRows)
    factorAggre= computeFactorAggr(factorData, involvedRows)
    aggrResults= pd.concat([df,  numericAggre, factorAggre], axis=1)
    return aggrResults
  

  def prepareDataForPairwiseAggregat(self, df):
    '''
      resturns featuresDicts to make it easier to compute aggregation of pairwise features
    
      @param df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R).
                      Column names are:
                      'chainIdL', 'structResIdL', 'resNameL', 'chainIdR', 'structResIdR', 'resNameR', 'categ' 
                      [properties_P .... propertiesL     .... propertiesR] #no defined order for properties
      @return #TODO write what this returns
                         
    '''
    df_ids= df.iloc[:, [0,1,3,4 ]].values
    pairwiseDf= df[ [elem for elem in df.columns if elem.endswith("P") and elem not in StructProtocol.IGNORE_FOR_AGGREGATION] ]

    ids2RowL={}
    ids2RowR={}
    for rowNum in range(df_ids.shape[0]):
      chainIdL, resIdL, chainIdR, resIdR= df_ids[rowNum,:]
      if (chainIdL,resIdL) not in ids2RowL:
        ids2RowL[(chainIdL,resIdL)]=set([])
      ids2RowL[(chainIdL,resIdL)].add(rowNum)
      
      if (chainIdR,resIdR) not in ids2RowR:
        ids2RowR[(chainIdR,resIdR)]=set([])
      ids2RowR[(chainIdR,resIdR)].add(rowNum)
     
    neigsids2rowL={}
    neigsids2rowR={}
    for resId in ids2RowL:
      neigsList= self.getNeigs(resId[0], resId[1], chainType="l")
      if neigsList is None:
        neigsids2rowL[resId]= set([])
      else:
        neigsids2rowL[resId]= set( itertoolchain.from_iterable([ ids2RowL[neig] for neig in neigsList if neig in ids2RowL]))

    for resId in ids2RowR:
      neigsList= self.getNeigs(resId[0], resId[1], chainType="r")
      if neigsList is None:
        neigsids2rowR[resId]= set([])
      else:
        neigsids2rowR[resId]= set( itertoolchain.from_iterable([ ids2RowR[neig] for neig in neigsList if neig in ids2RowR]))
        
    return pairwiseDf, ids2RowL, ids2RowR, neigsids2rowL, neigsids2rowR  



