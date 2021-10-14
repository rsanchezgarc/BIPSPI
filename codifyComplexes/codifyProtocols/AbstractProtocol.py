import os
import pandas as pd
import numpy as np
from itertools import chain

from codifyComplexes.CodifyComplexException import CodifyComplexException
from .DataLoaderClass import  DataLoader

#TODO: REWRITE PROTOCOLs CLASSES TO REDUCE COUPLING. E.G. Pairwise agregation should be generic for both seq and struct, just neigs different
FEATURES_MISMATH_TOLERANCE=0.1 #Fraction of missing residues between different features that will trigger error
class AbstractProtocol(DataLoader):
  '''
    This class is a template for codification protocols
  '''
  DESIRED_ORDER=["L", "R", "P"]
  def __init__(self, dataRootPath, cMapPath, prevStepPaths, singleChainfeatsToInclude, 
                     pairfeatsToInclude=None, verbose=False):
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
    if not hasattr(self, "dataRootPath"):
      DataLoader.__init__(self, dataRootPath, verbose)

    self.cMapPath= cMapPath
    self.singleChainfeatsToInclude= singleChainfeatsToInclude
    self.pairfeatsToInclude=  pairfeatsToInclude
    if prevStepPaths is None:
      self.prevFnamesList=None
    else:
      self.prevFnamesList= prevStepPaths if isinstance(prevStepPaths,list) else [prevStepPaths]
      self.prevFnamesList= [os.path.join(onePath, fname) for onePath in self.prevFnamesList for fname in os.listdir(onePath)
                                          if fname.endswith(".res.tab.gz")]
    if self.pairfeatsToInclude and not self.useCorrMut:
      self.pairfeatsToInclude = [ elem for elem in self.pairfeatsToInclude if elem[0]!="corrMut"]
      if len(self.pairfeatsToInclude)==0:
        self.pairfeatsToInclude = None
    
  def applyProtocol( self,  prefixComplex):
    '''
      This method is the basic skeleton for applyProtocol of subclasses
      Given a prefix that identifies the complex and prefixes that identifies
      the ligand and the receptor, this method integrates the information that
      is contained in self.dataRootPath and is described in self.singleChainfeatsToInclude
      
      :param prefixComplex: str. A prefix that identifies a complex
      :return df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R).
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ'
                       [propertiesP .... propertiesL     .... propertiesR]
    '''
    #load contact maps or previous results to define pairs of residues ligand to receptor
    raw_prefix= prefixComplex.split("@")[0].split("#")[0]
    if self.prevFnamesList is None:
      cmapNameList= list(self.getFullNamesIterInPath( raw_prefix, self.cMapPath))
      if len(cmapNameList)>1:
        raise ValueError("There are more than 1 Contact map for %s in %s path"%(prefixComplex,self.cMapPath))
      allPairsCodified= self.loadDataFile(cmapNameList)
    else:
      allPairsCodified= self.loadPreviousResults(prefixComplex)
    if not self.pairfeatsToInclude is None: #add pairwise features if there are available
      allPairsCodified= self.addPairFeatures(raw_prefix, allPairsCodified)

    lFeats= self.loadSingleChainFeatures( raw_prefix, chainType="l")
    rFeats= self.loadSingleChainFeatures( raw_prefix, chainType="r")
    #add single chain features to contact map (or pairwise features)
    allPairsCodified= self.combinePairwiseAndSingleChainFeats(allPairsCodified,lFeats, rFeats)
    assert allPairsCodified.shape[0]>1, "Error, %s dataset is empty"%prefixComplex
    #Reorder columns to DESIRED_ORDER order
    allPairsCodified= self.reorderColumns(allPairsCodified)
    return allPairsCodified
    
  def reorderColumns(self, allPairsCodified):
    colNames= list(allPairsCodified.columns)
    categIndex= colNames.index("categ")# categ is the last non-feature column. All previous columns are ids
    lFeatNames= [elem for elem in colNames[(categIndex+1):] if elem[-1]=="L"]
    rFeatNames= [elem for elem in colNames[(categIndex+1):] if elem[-1]=="R"]
    pairwiseFeatNames= [elem for elem in colNames[(categIndex+1):] if elem.endswith("_P")]
    colOrder= list(colNames[:(categIndex+1)]) #first columns are pair ids and label
    for featType in AbstractProtocol.DESIRED_ORDER:
      if featType=="L":
        colOrder+=lFeatNames
      elif featType=="R":
        colOrder+=rFeatNames
      elif featType=="P":
        colOrder+=pairwiseFeatNames
      else:
        raise ValueError("just L, R or P allowed in AbstractProtocol.DESIRED_ORDER")

    allPairsCodified= allPairsCodified[ colOrder ]
    allPairsCodified.columns= colOrder
    return allPairsCodified
    
  def loadPreviousResults(self, prefixComplex):
    '''
      Loads previous results. Returns a pandas.DataFrame that contains in each row
      the scores of previous steps for a given pair of amino acids.
      :param prefixComplex: str. A prefix (a complex id) that identifies the receptor or ligand

      :return df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      one amino acid
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR',
                      'categ', [prev_step_scores]
    '''
    previousResultsList= [fname for fname in self.prevFnamesList if fname.endswith(prefixComplex+".res.tab.gz")]
    assert len(previousResultsList)>0, "No previous results"
    for fname in previousResultsList:
      if self.verbose: print("loading previous resuls for %s"%(prefixComplex))
      prevResults= self.loadDataFile( iter([fname]) )
      prevResults.loc[:,"prediction_norm"]= (prevResults["prediction"] - np.mean(prevResults["prediction"])) / np.std(prevResults["prediction"]) 
      break # break to add just one type of previous predictions
    if len(previousResultsList)>1:
      for fname in previousResultsList[1:]: #Load the remaining previous predictions
        if self.verbose: print("loading previous resuls for %s"%(prefixComplex))

        prevNrow= prevResults.shape[0]
        newData= self.loadDataFile( iter([fname] ))
        newData.loc[:,"prediction_norm"]= (newData["prediction"] - np.mean(newData["prediction"])) / np.std(newData["prediction"]) 
        merge_on_cols= ["chainIdL",  "resIdL","resNameL","chainIdR", "resIdR","resNameR", "categ"]
        prevResults= pd.merge(prevResults, newData, how='inner', on=merge_on_cols) 
        curNrow= prevResults.shape[0]
        if prevNrow<1:
          raise CodifyComplexException(("Error merging previous results in %s. There are 0 rows "+
                                        "in previous features")%(prefixComplex ))
        elif (abs(float(prevNrow- curNrow)) / prevNrow) > FEATURES_MISMATH_TOLERANCE:
          raise CodifyComplexException(("Error merging previous results  in %s. There are a different number of residues "+
                                        "compared to previous file\nNrows previously/now %d/%d %s")%(prefixComplex,
                                                                                                  prevNrow, curNrow,fname))
    return prevResults
  
  def loadSingleChainFeatures(self, prefix, chainType):
    '''
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
    assert chainType=="l" or chainType=="r"
    prefixOneChainType= prefix+"_"+chainType
#    print(prefixOneChainType); raw_input("enter")
    oneChainTypeFeats= None
    featName, params= self.getParamsForLoadingFile(prefixOneChainType, self.singleChainfeatsToInclude[0]) #Load just first single chain feature
    if self.verbose: print("loading %s for %s"%(featName,prefixOneChainType))
    oneChainTypeFeats= self.loadDataFile(*params)
      
    if len(self.singleChainfeatsToInclude)>1:
      for featTuple in self.singleChainfeatsToInclude[1:]: #Load the remaining single chain features
        featName, params= self.getParamsForLoadingFile( prefixOneChainType, featTuple)
        if self.verbose: print("loading %s for %s"%(featName, prefixOneChainType))
        prevNrow= oneChainTypeFeats.shape[0]
        # params= ( list(params[0]), params[1] )
        # print(featTuple, list(params[0]))
        newData= self.loadDataFile(*params)
        oneChainTypeFeats_prev= oneChainTypeFeats
        oneChainTypeFeats= pd.merge(oneChainTypeFeats, newData, how='inner', on=["chainId",  "resId","resName"])      
        curNrow= oneChainTypeFeats.shape[0]
        if prevNrow<1:
          raise CodifyComplexException(("Error merging previous single chain feature %s in %s. There are 0 rows "+
                                        "in previous feature to %s")%(featName, prefixOneChainType, featName))
        elif (abs(float(prevNrow- curNrow)) / prevNrow) > FEATURES_MISMATH_TOLERANCE:
          if prevNrow>curNrow:
            df_diff= getDifRows_pd(oneChainTypeFeats_prev, oneChainTypeFeats, isSingleChain=True)
          else:
            df_diff= getDifRows_pd(oneChainTypeFeats, oneChainTypeFeats_prev, isSingleChain=True)
#          print(df_diff);raw_input()
#          print(oneChainTypeFeats_prev);raw_input()
#          print(oneChainTypeFeats);raw_input()
          errorMsg= str(df_diff)+"\n%s Nrows previously/now %d/%d  %s"%(prefixOneChainType, prevNrow, curNrow,featName)    
          raise CodifyComplexException((errorMsg+"\nError merging single chain feature %s in %s. There are a different number of residues "+
                                "in %s compared to previous features")%(featName, prefixOneChainType, featName))
    chainType= chainType.upper()
    oneChainTypeFeats.rename(columns={elem:elem+chainType for elem in list(oneChainTypeFeats.columns.values)}, inplace=True)
    return oneChainTypeFeats
    
  def addPairFeatures(self, prefix, allPairs):
    '''
      Loads all pairwise features files and adds them to the pairs of residues contained in allPairs df
      Returns a pandas.DataFrame that contains in each row all pairwise features from all files for each pair of amino acids. 
      Just amino acid pairs that appears in each file will be included. Others will be ruled out (intersection)
      :param prefix: str. A prefix that identifies the complex
      :param allPairs: pandas.DataFrame. A pandas.Dataframe in which each row represents one amino acid pair
                        Column names are:
                        "chainIdL", "resIdL","resNameL","chainIdR",  "resIdR","resNameR", "categ", [previousScores]
      :return df: pandas.DataFrame. A pandas.Dataframe in which each row represents one amino acid pair
                      Column names are:
                      "chainIdL",  "resIdL","resNameL","chainIdR",  "resIdR","resNameR", "categ", [previousScores] [pairFeats]
    '''

    pairTypeFeats= None

    featName, params= self.getParamsForLoadingFile( prefix, self.pairfeatsToInclude[0])
    if self.verbose: print("loading %s for %s"%(featName,prefix))
    pairTypeFeats= self.loadDataFile(*params)

    if len(self.pairfeatsToInclude)>1:
      for featTuple in self.pairfeatsToInclude[1:]: #Load the remaining single chain feature
        featName, params= self.getParamsForLoadingFile( prefix, featTuple)
        if self.verbose: print("loading %s for %s"%(featName, prefix))
        prevNrow= pairTypeFeats.shape[0]
        pairTypeFeats_prev= pairTypeFeats
        newData= self.loadDataFile(*params)
        pairTypeFeats= pd.merge(pairTypeFeats, newData, how='inner', on=["chainIdL",  "resIdL","resNameL",
                                                                         "chainIdR",  "resIdR","resNameR"])     
        curNrow= pairTypeFeats.shape[0]
        if prevNrow<1:
          raise CodifyComplexException(("Error merging previous pair feature %s in %s. There are 0 rows "+
                                        "in previous feature to %s")%(featName, prefix, featName))
        elif (abs(float(prevNrow- curNrow)) / prevNrow) > FEATURES_MISMATH_TOLERANCE**2:
          df_diff= getDifRows_pd(pairTypeFeats_prev, pairTypeFeats, isSingleChain=False)
          errorMsg= str(df_diff)+"\n%s Nrows previously/now %d/%d  %s"%(prefixOneChainType, prevNrow, curNrow,featName)    
          raise CodifyComplexException((errorMsg+"\nError merging pair feature %s in %s. There are a different number of residues "+
                                        "in %s compared to previous features")%(featName, prefixOneChainType, featName))

    #check if _P has already been assigned as mark to allPairs pairwise features
    nPs= sum( (1 for elem in pairTypeFeats.columns.values if elem.endswith("_P")))
    if nPs>0: #if so add _P to pairwise features
      pairTypeFeats.rename(columns={elem:elem+"_P" for elem in set(pairTypeFeats.columns.values
                                                       ).difference(AbstractProtocol.ARE_STR_TYPE_COLUMNS)}, inplace=True)
    prevNrow= allPairs.shape[0]
    allPairs_prev= allPairs
    allPairs= pd.merge(allPairs, pairTypeFeats, how='inner', on=["chainIdL",  "resIdL","resNameL",
                                                                       "chainIdR",  "resIdR","resNameR"])
    curNrow= allPairs.shape[0]
    if prevNrow<1:
      raise CodifyComplexException(("Error merging previous pair feature %s in %s. There are 0 rows "+
                                    "in previous feature to %s")%(featName, prefix, "allPairs"))
    elif (abs(float(prevNrow- curNrow)) / prevNrow) > FEATURES_MISMATH_TOLERANCE**2:
      df_diff= getDifRows_pd(allPairs_prev, allPairs, isSingleChain=False)
      errorMsg= str(df_diff)+"\n%s Nrows previously/now %d/%d  %s"%(prefixOneChainType, prevNrow, curNrow, "allPairs")    
      raise CodifyComplexException((errorMsg+"\nError merging pair feature %s in %s. There are a different number of residues "+
                                        "in %s compared to previous features")%(featName, prefixOneChainType, "allPairs"))                                                             
    return allPairs
    
  def combinePairwiseAndSingleChainFeats(self, pairFeats, singleFeatsL, singleFeatsR):
    '''
      Merges pairFeats pandas.DataFrame with singleFeatsL and singleFeatsR dataFrames.
      singleFeatsL has n rows (as many as ligand residues) and singleFeatsR has m rows
      (as many as ligand residues), and pairsFeats has ~ n*m rows
      (some amino acids pairs might not be considered due to several reasons).
    
      :param pairFeats: pandas.DataFrame. A pandas.Dataframe in which each row represents the properties of
                        a residue of the receptor
                        Column names are:
                        'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ' 
                         [propertiesP]    

      :param singleFeatsL: pandas.DataFrame. A pandas.Dataframe in which each row represents the properties of
                           a residue of the ligand
                           Column names are:
                           'chainIdL', 'resIdL', 'resNameL',[propertiesL] 
                           
      :param singleFeatsR: pandas.DataFrame. A pandas.Dataframe in which each row represents the properties of
                           a residue of the receptor
                           Column names are:
                           'chainIdR', 'resIdR', 'resNameR',[propertiesR]

      :return directPairs: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R).
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ' 
                      [ propertiesL .... propertiesR .... properties_P] 
                         
    '''
    if singleFeatsL.shape[1]!= singleFeatsR.shape[1]:
      print( "singleFeatsL and singleFeatsR have different number of variables")
      featsLNo_chain= set( (elem[:-1] for elem in singleFeatsL.columns))
      featsRNo_chain= set( (elem[:-1] for elem in singleFeatsR.columns))
      errorMsg= "\n".join(["L %d R %d"%(len(singleFeatsL.columns),len(singleFeatsR.columns)) ,"L dif R",
            str(sorted(featsLNo_chain.difference(featsRNo_chain))), "R diff L",
            str(sorted(featsRNo_chain.difference(featsLNo_chain))) ])
      raise ValueError( errorMsg+"\nsingleFeatsL and singleFeatsR have different number of variables")
    otherPairColumns= set(pairFeats.columns.values).difference(set(AbstractProtocol.ARE_STR_TYPE_COLUMNS+["categ"]))
    if len(otherPairColumns)>0: #add _P suffix to all pairwise features
      pairFeats.rename(columns={elem:elem+"_P" for elem in otherPairColumns}, inplace=True)      
    directPairs= pd.merge(pairFeats,   singleFeatsL, how='inner', on=None)   
    directPairs= pd.merge(directPairs, singleFeatsR, how='inner', on=None)
    return directPairs

        
  def prepareDataForPairwiseAggregat(self, df):
    '''
      abstract method
      :param   df: the Pandas.DataFrame of features before aggregation
      :return  pairwiseDf, ids2RowL, ids2RowR, neigsids2rowL, neigsids2rowR
          pairwiseDf: the Pandas.DataFrame of features after aggregation
          ids2RowL: a dict that maps from resIdL to all rows of df that contain resIdL
          ids2RowL: a dict that maps from resIdR to all rows of df that contain resIdR
          neigsids2rowL: a dict that maps from resIdL to all rows of df that contain a neighbour of resIdL
          neigsids2rowR: a dict that maps from resIdR to all rows of df that contain a neighbour of resIdR
    '''
    raise ValueError("Not implemented")
    

  def computeNumericAggr(self, df, selectedRows):
    return computeNumericAggr( df, selectedRows)
    
  def addPairwiseAggregation(self, df):
    '''
      Adds environment pairwise features to a df that contains pairwise features (named featName_P)
      :param df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (ligand, receptor)
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ' 
                      [ propertiesL .... propertiesR .... properties_P] #no defined order for properties
      
      :return newDf: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R). New columns will have been added, all then
                      named %sAggr%d%s(numAggr/factorAggr, numberOfNewFeature, ligand/receptor)
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ' 
                      [  propertiesL propertiesR properties_P ] #no defined order for properties
    '''

    pairwiseDf, ids2RowL, ids2RowR, neigsids2rowL, neigsids2rowR = self.prepareDataForPairwiseAggregat(df)
    nElems= df.shape[0]

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

    numericAggreL2r= self.computeNumericAggr(pairwiseDf, l_2_r_neigs)
    numericAggreR2l= self.computeNumericAggr(pairwiseDf, r_2_l_neigs)
    numericAggreN2N= self.computeNumericAggr(pairwiseDf, l_neigs_2_r_neigs)
    numericAggreL2r.columns= [ "l2r-pair"+elem for elem in numericAggreL2r.columns]
    numericAggreR2l.columns= [ "r2l-pair"+elem for elem in numericAggreR2l.columns]
    numericAggreN2N.columns= [ "n2n-pair"+elem for elem in numericAggreN2N.columns]
    df= pd.concat([df,  numericAggreL2r, numericAggreR2l, numericAggreN2N ], axis=1)
    return df

def computeNumericAggr(df, selectedRows):
  '''
    Compute aggregation functions (min, max, mean and sum) for the rows of data. Each row is aggregated
    over the rows included in selectedRows
    :param data: pandas.DataFrame(nrow, nfeatures). The numeric data to aggregate. Each row i is averaged, sum... with
                                          the rows such that selectedRows[i,:]==True 
    :param selectedRows: [[]]: int (nrow, variableLength)): Each row i is averaged, with the rows included in
                                                    selectedRows[i]    
    :return np.array (nrow, 4*nfeatures)
  '''
  data= df.values
  nVars= data.shape[1]
  aggregatedResults= -(2**10)* np.ones( (data.shape[0], nVars* 4))
  splitPoint2= 2*nVars
  splitPoint3= 3*nVars
  splitPoint4= 4*nVars
  for i in range(data.shape[0]):
    dataToAggregate= data[selectedRows[i], :]

    if dataToAggregate.shape[0]>0:
      aggregatedResults[i, 0:nVars]=                np.mean( dataToAggregate, axis=0)
      aggregatedResults[i, nVars:splitPoint2]=      np.max(  dataToAggregate, axis=0)
      aggregatedResults[i, splitPoint2:splitPoint3]= np.min( dataToAggregate, axis=0)
      aggregatedResults[i, splitPoint3:splitPoint4]= np.sum( dataToAggregate, axis=0)
      

  aggregatedResults= pd.DataFrame(aggregatedResults)
  aggregatedResults.columns= list(chain.from_iterable([[ "numAggr-"+oper+"-"+name for name in df.columns] for 
                                                      oper in ["mean", "max", "min", "sum"] ]))  
  return aggregatedResults


def computeFactorAggr(df, selectedRows):
  '''
    Compute aggregation function (sum) for the rows of data. Each row is aggregated
    over the rows included in selectedRows. This is equivalent to count the number of neighbours of each category
    :param data: pandas.DataFrame(nrow, nfeatures). The numeric data to aggregate. Each row i is averaged,  with
                                          the rows such that selectedRows[i,:]==True 
    :param selectedRows: [[]]: int (nrow, variableLength)): Each row i is added with the rows included in
                                                    selectedRows[i]
                           
    :return np.array (nrow, nfeatures)
  '''
  data= df.values
  aggregatedResults= -(2**10)* np.ones( (data.shape[0], data.shape[1]))
  for i in range(data.shape[0]):
    dataToAggregate= data[selectedRows[i],:]
    if dataToAggregate.shape[0]>0:
      aggregatedResults[i, :]= np.sum( dataToAggregate, axis=0)
  aggregatedResults= pd.DataFrame(aggregatedResults)
  aggregatedResults.columns= [ "factorAggr-"+name for name in df.columns] 
  return aggregatedResults


def getDifRows_pd(df1, df2, isSingleChain=True):
  '''
    detects the different rows in two dataframes. Used for debugging
  '''
  print("shapes:", df1.shape, df2.shape)
  nCols= 3 if isSingleChain else 6
  df1= df1.iloc[:,:nCols]
  df2= df2.iloc[:,:nCols]
#  df1.to_csv("/home/rsanchez/tmp/error1.txt")
#  df2.to_csv("/home/rsanchez/tmp/error2.txt")
  common = df1.merge(df2)
  dif_set= pd.concat([df1,df2]).drop_duplicates(keep=False)
  print(dif_set.shape)
  print(dif_set.head())
  print("-----------")
  return dif_set


