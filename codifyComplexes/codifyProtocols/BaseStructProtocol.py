from __future__ import absolute_import, print_function
import os
import pandas as pd
import numpy as np

from itertools import chain as itertoolchain

from .AbstractProtocol import AbstractProtocol, computeNumericAggr, computeFactorAggr
from .NeigbManager import NeigbManager
from codifyComplexes.CodifyComplexException import CodifyComplexException
from computeFeatures.seqStep.seqToolManager import SeqToolManager

AA_CODE_ELEMENTS= SeqToolManager.AA_CODE_ELEMENTS


class BaseStructProtocol(AbstractProtocol, NeigbManager):
  '''
    This class implements structural voronoi environment codification
  '''
  IGNORE_FOR_AGGREGATION=["ccmPredQuality_P", "psicovQuality_P"]
  def __init__(self, dataRootPath, cMapPath, prevStepPaths=None, singleChainfeatsToInclude=None, 
                     pairfeatsToInclude=None, verbose=False ):
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
      :param singleChainfeatsToInclude: List that contains the paths where single chain features needed for complex
                                      codification are located. Must have the following format:
                                      ["featName":(relativePath_from_dataRootPath, listOfColumnNumbers)]
                                         
      :param pairfeatsToInclude: List that contains the paths where pair features needed for complex
                                      codification are located. Must have the following format:
                                      ["featName":(relativePath_from_dataRootPath, listOfColumnNumbers)]
      :param verbose: bool.
    
    '''
    if singleChainfeatsToInclude is None: raise ValueError("Error, singleChainfeatsToInclude must be provided")
    NeigbManager.__init__(self, dataRootPath) 
    AbstractProtocol.__init__(self, dataRootPath, cMapPath, prevStepPaths, singleChainfeatsToInclude=singleChainfeatsToInclude,
                                    pairfeatsToInclude= pairfeatsToInclude, verbose= verbose)

  def applyProtocol( self, prefixComplex):
    '''
      This method is the basic skeleton for applyProtocol of subclasses
      Given a prefix that identifies the complex and prefixes that identifies
      the ligand and the receptor, this method integrates the information that
      is contained in self.dataRootPath and is described in self.singleChainfeatsToInclude
      
      :param prefixComplex: str. A prefix that identifies a complex, e.g. 1ACB
      :return df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R).
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ'
                       [ propertiesL .... propertiesR .... propertiesP]
    '''
    self.loadEnvirons(prefixComplex )
    allPairsCodified= super(BaseStructProtocol,self).applyProtocol( prefixComplex )
    allPairsCodified= self.addPairwiseAggregation(allPairsCodified)
    allPairsCodified= self.reorderColumns(allPairsCodified)    
    return allPairsCodified

    
    
  def loadSingleChainFeatures(self, prefix, chainType):
    '''
      @overrides AbstractProtocol method. Computes structural neighbour aggregation on every feature (if sliding window
                 features included, the final number of variables will be huge, use StructProtocol instead)
      Loads all features files computed for ligand or receptor chains. Returns a pandas.DataFrame 
      that contains in each row all features from all files for each amino acid. Just amino acids
      that appears in each file will be included. Others will be ruled out (intersection)
      :param prefix: str. A prefixOneChainType that identifies the receptor or ligand. e.g. 1A2K
      :param chainType: str. "l" for ligand and "r" for receptor
      :return df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      one amino acid
                      Column names are:
                      'chainId%s', 'resId%s', 'resName%s', [properties]
                      %s is L if chainType=="l" and R if  chainType=="r"
    '''
    singleChainFeats= super(BaseStructProtocol,self).loadSingleChainFeatures( prefix, chainType)
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
      :param df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      an amino acid of ligand or receptor.
                      Column names are:
                      'chainId%s', 'resId%s', 'resName%s',
                      [propertiesL .... propertiesR .... properties_P] #no defined order for properties
      :param chainType: str. "l" for ligand and "r" for receptor
      
      :return newDf: pandas.DataFrame. A pandas.Dataframe in which each row represents an amino acid of
                      ligand or receptor. New columns will be added, all then
                      named "%sAggr%d%s"%(num/factor, numberOfNewFeature, ligand/receptor)
                      Column names are:
                      'chainId%s', 'resId%s', 'resName%s',
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
#    id2Name={ i: fullId for fullId,i in name2Id.iteritems()}

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
      resturns different data to make it easier to compute aggregation of pairwise features
    
      :param df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R).
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ' 
                      [ propertiesL .... propertiesR .... properties_P] #no defined order for properties

      :return  pairwiseDf, ids2RowL, ids2RowR, neigsids2rowL, neigsids2rowR
          pairwiseDf: the Pandas.DataFrame of features after aggregation
          ids2RowL: a dict that maps from resIdL to all rows of df that contain resIdL
          ids2RowL: a dict that maps from resIdR to all rows of df that contain resIdR
          neigsids2rowL: a dict that maps from resIdL to all rows of df that contain a neighbour of resIdL
          neigsids2rowR: a dict that maps from resIdR to all rows of df that contain a neighbour of resIdR
    '''
    categIndex= list(df.columns).index("categ")
    df_ids= df.iloc[:, :(categIndex+1) ].values
    pairwiseDf= df[ [elem for elem in df.columns if elem.endswith("P") and elem not in BaseStructProtocol.IGNORE_FOR_AGGREGATION] ]

    ids2RowL={}
    ids2RowR={}
    for rowNum in range(df_ids.shape[0]):
      chainIdL, resIdL, __, chainIdR, resIdR, __= df_ids[rowNum,:6]
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
#      print(self.getNeigs(resId[0], resId[1], chainType="r")); raw_input("enter")
      neigsList= self.getNeigs(resId[0], resId[1], chainType="r")
      if neigsList is None:
        neigsids2rowR[resId]= set([])
      else:
        neigsids2rowR[resId]= set( itertoolchain.from_iterable([ ids2RowR[neig] for neig in neigsList if neig in ids2RowR]))
        
    return pairwiseDf, ids2RowL, ids2RowR, neigsids2rowL, neigsids2rowR  



