from __future__ import absolute_import, print_function

import pandas as pd
import numpy as np

from .StructProtocol import StructProtocol
from .SeqProtocol import SeqProtocol
from codifyComplexes.CodifyComplexException import CodifyComplexException
from .AbstractProtocol import AbstractProtocol, computeNumericAggr, computeFactorAggr
from codifyComplexes.CodifyComplexException import CodifyComplexException
from .NeigbManager import NeigbManager
from computeFeatures.seqStep.seqToolManager import SeqToolManager

AA_CODE_ELEMENTS= SeqToolManager.AA_CODE_ELEMENTS
'''
(feature_name, path_to_dir, columns ). If columns==None, all columns will be used
Structural features must come first as sequential single chains features muy contain more aminoacids
(e.g. non 3d-solved residues)
'''
FEATURES_TO_INCLUDE_CHAIN_STRUCT= [
  ("psaia", ("structStep/PSAIA/procPSAIA", None)),
  ("halfSphereExpos", ("structStep/halfSphereExpos", None)),
  ("dssp", ("structStep/DSSP", [3])),
  ("al2co", ("seqStep/conservation/al2co",[4])),
  ("winPssms", ("seqStep/conservation/pssms/windowedPSSMs/wSize11", None)),
  ("winSeq", ("seqStep/slidingWinSeq11", None))
]


FEATURES_TO_INCLUDE_CHAIN_SEQ= [
  ("winPssms", ("seqStep/conservation/pssms/windowedPSSMs/wSize11", None)),
  ("winSeq", ("seqStep/slidingWinSeq11", None)),
  ("al2co", ("seqStep/conservation/al2co",[4])),
  ("predAsaAndSS", ("seqStep/SPIDER2_predAsaSS",None))
]

FEATURES_TO_INCLUDE_PAIR= [
  ("corrMut", ("seqStep/conservation/corrMut", None)),
]

class MixedProtocol(SeqProtocol, StructProtocol):
  '''
    This class implement sequentail environment codification for one partner and sequence and voronoi environment 
    codification for the other partner.
    Used when one of the partners is a sequence and the other is an structural model
  '''
  DESIRED_ORDER=["Q", "T", "P"] #sequenceBased, structBased and Pairwise
  IGNORE_FOR_AGGREGATION=["ccmPredQuality_P", "psicovQuality_P"]
  def __init__(self, dataRootPath, cMapPath, prevStepPaths=None, verbose=False):
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
                                 to define which residue pairs are in contact. Can also be a str[] if multiple feedback_paths
                                 wanted
      :param singleChainfeatsToInclude: List that contains the paths where single chain features needed for complex
                                      codification are located. Must have the following format:
                                      ["featName":(relativePath_from_dataRootPath, listOfColumnNumbers)]
                                         
      :param pairfeatsToInclude: List that contains the paths where pair features needed for complex
                                      codification are located. Must have the following format:
                                      ["featName":(relativePath_from_dataRootPath, listOfColumnNumbers)]
      :param verbose: bool.
    '''
    NeigbManager.__init__(self, dataRootPath ) 

    SeqProtocol.__init__(self, dataRootPath, cMapPath, prevStepPaths,
                                singleChainfeatsToInclude= FEATURES_TO_INCLUDE_CHAIN_SEQ, 
                                pairfeatsToInclude=FEATURES_TO_INCLUDE_PAIR, verbose= verbose)
    #self.singleChainfeatsToInclude are the features that describe one single chain and will be loaded 

  def checkWhoIsSequenceOnly(self, prefixComplex):

    prefixComplex= prefixComplex.split("@")[0]
    if prefixComplex[-3:]=="#sl":
      return ["l"]
    elif prefixComplex[-3:]=="#sr":
      return ["r"]
    else:
      raw_prefix= prefixComplex.split("@")[0].split("#")[0]
      isSeqOnly=[]
      __, (fnamesIteratorL, selectedCols)= self.getParamsForLoadingFile( raw_prefix+"_l", type(self).RASA_FEAT_DESCR)
      if len(list(fnamesIteratorL))==0:
        isSeqOnly.append("l")
      __, (fnamesIteratorR, selectedCols)= self.getParamsForLoadingFile( raw_prefix+"_r", type(self).RASA_FEAT_DESCR)
      if len(list(fnamesIteratorR))==0:
        isSeqOnly.append("r")
      return isSeqOnly
    
  def applyProtocol( self, prefixComplex, sequenceOnly="l"):
    '''
      @overrides SeqProtocol applyProtocol
      This method is the basic skeleton for applyProtocol of subclasses
      Given a prefix that identifies the complex and prefixes that identifies
      the ligand and the receptor, this method integrates the information that
      is contained in self.dataRootPath and is described in self.singleChainfeatsToInclude
      
      :param prefixComplex: str. A prefix that identifies a complex
      :param sequenceOnly: str. "l" or "r". Indicates whether the ligand or the receptor has no structural information
      :return df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R).
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ'
                       [ propertiesL .... propertiesR .... propertiesP]
    '''
    assert sequenceOnly in ["l","r"], "Error, sequence only must be 'l' or 'r'"
    raw_prefix= prefixComplex.split("@")[0].split("#")[0]
    self.loadEnvirons(raw_prefix, ligandHas= False if sequenceOnly=="l" else True, 
                                   receptorHas= False if sequenceOnly=="r" else True)
      
    #load sequence features as Sequential protocol does
    allPairsCodified= SeqProtocol.applyProtocol( self, prefixComplex)
    
    chainTypeWithStruct= "l" if sequenceOnly=="r" else "r"
    structFeatures= self.loadSingleChainFeatures( prefixComplex, chainType=chainTypeWithStruct, loadAsSeqProto=False)
    structFeatures_ids= structFeatures.iloc[:, :3]
    structFeatures= structFeatures.iloc[:,3:]
    
    seqFeatNames= list(allPairsCodified.columns)
    structFeatsNamesToAdd= list(set(structFeatures.columns).difference(set(seqFeatNames)))
    structFeatures= structFeatures[structFeatsNamesToAdd]
    categIndex= seqFeatNames.index("categ")
    allPairsCodified.columns= seqFeatNames[:categIndex+1]+[name if name.endswith("_P") else name+"_seQ" for name in seqFeatNames[categIndex+1:] ]
    
    structFeatures.columns= [name+"_strucT" for name in structFeatsNamesToAdd ] 
    structDf= pd.concat([structFeatures_ids, structFeatures], axis=1)

    allPairsCodified= pd.merge(allPairsCodified,   structDf, how='inner', on=None)
    self.addPairwiseMixedAggregation( allPairsCodified, chainTypeWithStruct)
    allPairsCodified= self.reorderColumns(allPairsCodified, orderSeqStruct=True)
    return allPairsCodified
    
  def loadSingleChainFeatures(self, prefix, chainType, loadAsSeqProto=True):
    '''
      @overrides SequenceProtocol method to make use of sequence profiles (loaded directly) and struct
                  neighbour but not computing struct neighbours on non central residue features of sliding window
      Loads all features files computed for ligand or receptor chains. Returns a pandas.DataFrame 
      that contains in each row all features from all files for each amino acid. Just amino acids
      that appears in each file will be included. Others will be ruled out (intersection)
      :param prefix: str. A prefixOneChainType that identifies the receptor or ligand. e.g. 1A2K
      :param chainType: str. "l" for ligand and "r" for receptor
      :return df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      one amino acid
                      Column names are:
                      'chainId%s', 'resId%s', 'resName%s', [properties] #no defined order for properties
                      %s is L if chainType=="l" and R if  chainType=="r"
    '''
    raw_prefix= prefix.split("@")[0].split("#")[0]
    if loadAsSeqProto:
      singleChainFeats= SeqProtocol.loadSingleChainFeatures( self, raw_prefix, chainType)
    else:
      singleChainfeatsToInclude= self.singleChainfeatsToInclude
      self.singleChainfeatsToInclude= FEATURES_TO_INCLUDE_CHAIN_STRUCT #Use struct features 
      singleChainFeats= StructProtocol.loadSingleChainFeatures( self, raw_prefix, chainType)
      self.singleChainfeatsToInclude= singleChainfeatsToInclude #Set seq features as it was previously
      
    return singleChainFeats
    

  def reorderColumns(self, allPairsCodified, orderSeqStruct=False):

    if not orderSeqStruct:
      return SeqProtocol.reorderColumns(self, allPairsCodified)
    else:
      colNames= list(allPairsCodified.columns)
      categIndex= colNames.index("categ")# categ is the last non-feature column. All previous columns are ids
      seqFeatNames= [elem for elem in colNames[(categIndex+1):] if elem[-1]=="Q"]
      structFeatNames= [elem for elem in colNames[(categIndex+1):] if elem[-1]=="T"]
      pairwiseFeatNames= [elem for elem in colNames[(categIndex+1):] if elem.endswith("_P")]
      colOrder= list(colNames[:(categIndex+1)]) #first columns are pair ids and label

      for featType in MixedProtocol.DESIRED_ORDER:
        if featType=="Q":
          colOrder+=seqFeatNames
        elif featType=="T":
          colOrder+=structFeatNames
        elif featType=="P":
          colOrder+=pairwiseFeatNames
        else:
          raise ValueError("just Q,T,P allowed in MixedProtocol.DESIRED_ORDER")
      allPairsCodified= allPairsCodified[ colOrder ]
      allPairsCodified.columns= colOrder
      return allPairsCodified


  def detectPartnerWithStruct(self,df):
    colNames= list(df.columns)
    chainTypes= set([])
    for colName in colNames:
      if colName.endswith("_strucT"):
        chainTypes.add( colName.split("_")[-2][-1].lower() )
    if len(chainTypes)<1:
      return None
    elif len(chainTypes)==1:
      return chainTypes.pop()
    else:
      return "Error, in mixed protocol just one chain can contain structural information"
    
  def addPairwiseMixedAggregation(self, df,  chainTypeWithStruct):
    '''
      Adds environment pairwise features to a df that contains pairwise features (named featName_P)
      :param df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (ligand, receptor)
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ' 
                      [ propertiesL .... propertiesR .... properties_P] #no defined order for properties
      
      :param  chainTypeWithStruct: str. "l" or "r" chainType that has structural and sequence-based info

      :return newDf: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R). New columns will have been added, all then
                      named %sAggr%d%s(numAggr/factorAggr, numberOfNewFeature, ligand/receptor)
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ' 
                      [  propertiesL propertiesR properties_P ] #no defined order for properties
    '''
    

    hasPairwise= bool( sum( (1   for colName in df.columns if colName.endswith("_P"))))
    if not hasPairwise:
      return df
    pairwiseDf, ids2RowL, ids2RowR, neigsids2rowL, neigsids2rowR = StructProtocol.prepareDataForPairwiseAggregat(self, df)
        
    l_2_r_neigs= [[] for i in range(df.shape[0])]
    r_2_l_neigs= [[] for i in range(df.shape[0])]
    l_neigs_2_r_neigs= [[] for i in range(df.shape[0])]
    for idL in ids2RowL:
      rowsInvolvingL= ids2RowL[idL]
      rowsNeigsL= neigsids2rowL[idL]
      for idR in ids2RowR:
        rowsInvolvingR= ids2RowR[idR]
        rowsNeigsR= neigsids2rowR[idR]
        row_Index= tuple(rowsInvolvingL.intersection(rowsInvolvingR))[0]
        l_2_r_neigs[row_Index]= tuple(rowsInvolvingL.intersection(rowsNeigsR))
        r_2_l_neigs[row_Index]= tuple(rowsInvolvingR.intersection(rowsNeigsL))
        l_neigs_2_r_neigs[row_Index]= tuple(rowsNeigsR.intersection(rowsNeigsL))

    if chainTypeWithStruct=="l":
      numericAggre= self.computeNumericAggr(pairwiseDf, r_2_l_neigs)
      numericAggre.columns= [ "r2l-pair"+elem for elem in numericAggre.columns]
    elif chainTypeWithStruct=="r":
      numericAggre= self.computeNumericAggr(pairwiseDf, l_2_r_neigs)
      numericAggre.columns= [ "l2r-pair"+elem for elem in numericAggre.columns]
    else:
      raise ValueError("chainTypeWithStruct must be l or r")

    df= pd.concat([df,  numericAggre ], axis=1)
    return df
          
