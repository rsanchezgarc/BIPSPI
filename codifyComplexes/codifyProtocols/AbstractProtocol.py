import os
import pandas as pd
import numpy as np
from itertools import chain

from codifyComplexes.CodifyComplexException import CodifyComplexException
from .DataLoaderClass import  DataLoader

class AbstractProtocol(DataLoader):
  '''
    This class is a template for codification protocols
  '''

  DESIRED_ORDER=["L", "R", "P"]
  def __init__(self, dataRootPath, cMapPath, prevStepPaths, singleChainfeatsToInclude, pairfeatsToInclude=None, verbose=False):
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
    
      @param singleChainfeatsToInclude: List that contains the paths where single chain features needed for complex 
                                      codification are located. Must have the following format:
                                      ["featName":(relativePath_from_dataRootPath, listOfColumnNumbers, dictForNamedColums)]
                                          dictForNamedColums= {"myFeatName":colNumber} or Empty
      @param pairfeatsToInclude: List that contains the paths where pair features needed for complex 
                                      codification are located. Must have the following format:
                                      ["featName":(relativePath_from_dataRootPath, listOfColumnNumbers, dictForNamedColums)]
                                          dictForNamedColums= {"myFeatName":colNumber} or {}
      @param verbose: bool.
    '''
    DataLoader.__init__(self, dataRootPath, singleChainfeatsToInclude, pairfeatsToInclude, verbose)
    self.dataRootPath= dataRootPath
    self.cMapPath= cMapPath
    self.verbose= verbose
    if prevStepPaths is None:
      self.prevFnamesList=None
    else:
      self.prevFnamesList= prevStepPaths if isinstance(prevStepPaths,list) else [prevStepPaths]
      self.prevFnamesList= [os.path.join(onePath, fname) for onePath in self.prevFnamesList for fname in os.listdir(onePath)
                                          if fname.endswith(".res.tab")]
    if not self.useCorrMut:
      self.pairfeatsToInclude = [ elem for elem in self.pairfeatsToInclude if elem[0]!="corrMut"]
      if len(self.pairfeatsToInclude)==0:
        self.pairfeatsToInclude = None
    
  def applyProtocol( self,  prefixComplex, prefixL, prefixR):
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
 
    if self.prevFnamesList is None:
      cmapNameList= list(self.getFullNamesIterInPath( prefixComplex, self.cMapPath))
      if len(cmapNameList)>1:
        raise ValueError("There are more than 1 Contact map for %s in %s path"%(prefixComplex,self.cMapPath))
      allPairsCodified= self.loadDataFile(cmapNameList)
    else:
      allPairsCodified= self.loadPreviousResults(prefixComplex)
    if not self.pairfeatsToInclude is None:
      allPairsCodified= self.addPairFeatures(prefixComplex, allPairsCodified)
    lFeats= self.loadSingleChainFeatures( prefixL, "l")
    rFeats= self.loadSingleChainFeatures( prefixR, "r")
#    print(lFeats["chainIdL"].unique(), rFeats["chainIdR"].unique(), allPairsCodified["chainIdL"].unique(), allPairsCodified["chainIdR"].unique())
    allPairsCodified= self.combinePairwiseAndSingleChainFeats(allPairsCodified,lFeats, rFeats)
    assert allPairsCodified.shape[0]>1, "Error, %s dataset is empty"%prefixComplex
    #Reorder columns
    allPairsCodified= self.reorderColumns(allPairsCodified)
    return allPairsCodified
    
  def reorderColumns(self, allPairsCodified):
    colNames= list(allPairsCodified.columns)
    categIndex= colNames.index("categ")
    lFeatNames= [elem for elem in colNames[(categIndex+1):] if elem[-1]=="L"]
    rFeatNames= [elem for elem in colNames[(categIndex+1):] if elem[-1]=="R"]
    pairwiseFeatNames= [elem for elem in colNames[(categIndex+1):] if elem.endswith("_P")]
    colOrder= list(colNames[:(categIndex+1)])
    for featType in AbstractProtocol.DESIRED_ORDER:
      if featType=="L":
        colOrder+=lFeatNames
      elif featType=="R":
        colOrder+=rFeatNames
      elif featType=="P":
        colOrder+=pairwiseFeatNames
      else:
        raise ValueError("just L, R or P allowed in idAbstractProtocol.DESIRED_ORDER")

    allPairsCodified= allPairsCodified[ colOrder ]
    allPairsCodified.columns= colOrder

    return allPairsCodified
    
  def loadPreviousResults(self, prefixComplex):
    '''
      Loads previous results. Returns a pandas.DataFrame that contains in each row
      the scores of previous steps for a given pair of amino acids.
      @param prefixComplex: str. A prefixOneChainType that identifies the receptor or ligand

      @return df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      one amino acid
                      Column names are:
                      'chainIdL', 'structResIdL', 'resNameL', 'chainIdR', 'structResIdR', 'resNameR',
                      'categ', [prev_step_scores]
    '''
    previousResultsList= [fname for fname in self.prevFnamesList if fname.endswith(prefixComplex+".res.tab")]
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
        
        prevResults= pd.merge(prevResults, newData, how='inner', on=["chainIdL",  "structResIdL","resNameL",
                                                                     "chainIdR",  "structResIdR","resNameR","categ"]) 
        curNrow= prevResults.shape[0]
        if prevNrow<1:
          raise CodifyComplexException(("Error merging previous results in %s. There are 0 rows "+
                                        "in previous features")%(prefixComplex ))
        elif (abs(float(prevNrow- curNrow)) / prevNrow) > 0.2:
          print("%s Nrows previously/now %d/%d  %s"%(prefixComplex, prevNrow, curNrow,fname))
          raise CodifyComplexException(("Error merging previous results  in %s. There are a different number of residues "+
                                        "compared to previous file")%(prefixComplex))
    return prevResults
  
  def loadSingleChainFeatures(self, prefixOneChainType, chainType):
    '''
      Loads all features files computed for ligand or receptor chains. Returns a pandas.DataFrame 
      that contains in each row all features from all files for each amino acid. Just amino acids
      that appears in each file will be included. Others will be ruled out (intersection)
      @param prefixOneChainType: str. A prefixOneChainType that identifies the receptor or ligand
      @param chainType: str. "l" for ligand and "r" for receptor
      @return df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      one amino acid
                      Column names are:
                      'chainId%s', 'structResId%s', 'resName%s', [properties]
                      %s is L if chainType=="l" and R if  chainType=="r"
    '''

    assert chainType=="l" or chainType=="r"
    chainType= chainType.upper()
    oneChainTypeFeats= None
#    self.verbose=True
    for featNum in range(len(self.singleChainfeatsToInclude)): #Load just one single chain feature
      featName, params= self.getParamsForLoadingFile( prefixOneChainType, featNum)
      if self.verbose: print("loading %s for %s"%(featName,prefixOneChainType))
      oneChainTypeFeats= self.loadDataFile(*params)
      break # break to add just one type of feature. Next will be added in following loop
    if len(self.singleChainfeatsToInclude)>1:
      for featNum in range(len(self.singleChainfeatsToInclude))[1:]: #Load the remaining single chain feature
        featName, params= self.getParamsForLoadingFile( prefixOneChainType, featNum)
        if self.verbose: print("loading %s for %s"%(featName, prefixOneChainType))
        prevNrow= oneChainTypeFeats.shape[0]
        newData= self.loadDataFile(*params)
#        print(newData.head())
#        print(newData.shape, newData["chainId"].unique())
#        print(prefixOneChainType, chainType, featName)
#        raw_input("enter to continue")
#        print(oneChainTypeFeats.head())
#        print(oneChainTypeFeats.shape, oneChainTypeFeats["chainId"].unique())
#        raw_input("enter to continue")
#        print(list(oneChainTypeFeats.columns), list(newData.columns),oneChainTypeFeats.shape)
        oneChainTypeFeats= pd.merge(oneChainTypeFeats, newData, how='inner', on=["chainId",  "structResId","resName"]) 
#        print(oneChainTypeFeats.columns, oneChainTypeFeats.shape)
#        raw_input()        
        curNrow= oneChainTypeFeats.shape[0]
        if prevNrow<1:
          raise CodifyComplexException(("Error merging previous single chain feature %s in %s. There are 0 rows "+
                                        "in previous feature to %s")%(featName, prefixOneChainType, featName))
        elif (abs(float(prevNrow- curNrow)) / prevNrow) > 0.2:
          print("Nrows previously/now", prevNrow, curNrow)
          print("%s Nrows previously/now %d/%d  %s"%(prefixOneChainType, prevNrow, curNrow,featName))          
          raise CodifyComplexException(("Error merging single chain feature %s in %s. There are a different number of residues "+
                                        "in %s compared to previous features")%(featName, prefixOneChainType, featName))

    oneChainTypeFeats.rename(columns={elem:elem+chainType for elem in list(oneChainTypeFeats.columns.values)}, inplace=True)
    return oneChainTypeFeats
    
  def addPairFeatures(self, prefix, allPairs):
    '''
      Loads all files computed for pairwise features and adds them to pairs residue df contained in allPairs. 
      Returns a pandas.DataFrame that contains in each row all pairwise features from all files for each pair of amino acids. 
      Just amino acid pairs that appears in each file will be included. Others will be ruled out (intersection)
      @param prefix: str. A prefix that identifies the complex
      @param allPairs: pandas.DataFrame. A pandas.Dataframe in which each row represents one amino acid pair
                        Column names are:
                        "chainIdL",  "structResIdL","resNameL","chainIdR",  "structResIdR","resNameR", "categ", [previousScores]
      @return df: pandas.DataFrame. A pandas.Dataframe in which each row represents one amino acid pair
                      Column names are:
                      "chainIdL",  "structResIdL","resNameL","chainIdR",  "structResIdR","resNameR", "categ", [properties]
    '''

    pairTypeFeats= None
    self.verbose=True
    for featNum in range(len(self.pairfeatsToInclude)): #Load just one single chain feature
      featName, params= self.getParamsForLoadingFile( prefix, featNum, lookForPairFeats=True)
      if self.verbose: print("loading %s for %s"%(featName,prefix))
      pairTypeFeats= self.loadDataFile(*params)
      break # break to add just one type of feature. Next will be added in following loop
    if len(self.pairfeatsToInclude)>1:
      for featNum in range(len(self.pairfeatsToInclude))[1:]: #Load the remaining single chain feature
        featName, params= self.getParamsForLoadingFile( prefix, featNum, lookForPairFeats=True)
        if self.verbose: print("loading %s for %s"%(featName, prefix))
        prevNrow= pairTypeFeats.shape[0]
        newData= self.loadDataFile(*params)
        pairTypeFeats= pd.merge(pairTypeFeats, newData, how='inner', on=["chainIdL",  "structResIdL","resNameL",
                                                                         "chainIdR",  "structResIdR","resNameR"])     
        curNrow= pairTypeFeats.shape[0]
        if prevNrow<1:
          raise CodifyComplexException(("Error merging previous pair feature %s in %s. There are 0 rows "+
                                        "in previous feature to %s")%(featName, prefix, featName))
        elif (abs(float(prevNrow- curNrow)) / prevNrow) > 0.2:
          print("Nrows previously/now", prevNrow, curNrow)
          print("%s Nrows previously/now %d/%d  %s"%(prefix, prevNrow, curNrow,featName))          
          raise CodifyComplexException(("Error merging pair feature %s in %s. There are a different number of residues "+
                                        "in %s compared to previous features")%(featName, prefix, featName))

    #check if _P has being assigned to allPairs pairwise features
    nPs= sum( (1 for elem in pairTypeFeats.columns.values if elem.endswith("_P")))
    if nPs>0: #if so add _P to pairwise features
      pairTypeFeats.rename(columns={elem:elem+"_P" for elem in set(pairTypeFeats.columns.values
                                                           ).difference(AbstractProtocol.ARE_STR_TYPE_COLUMNS)}, inplace=True)
#    print(allPairs.shape)
    prevNrow= allPairs.shape[0]
    allPairs= pd.merge(allPairs, pairTypeFeats, how='inner', on=["chainIdL",  "structResIdL","resNameL",
                                                                       "chainIdR",  "structResIdR","resNameR"])

#    print(allPairs.head())
#    print(allPairs.shape)
#    raw_input("enter")
    curNrow= allPairs.shape[0]
    if abs(float(prevNrow- curNrow)) / prevNrow > 0.05:
      print("Nrows previously/now", prevNrow, curNrow)
      print("%s Nrows previously/now %d/%d  %s"%(prefix, prevNrow, curNrow,featName))          
      raise CodifyComplexException(("Error merging pair feature %s in %s. There are a different number of residues "+
                                    "in %s compared to previous features")%(featName, prefix, featName))                                                                  
    return allPairs
    
  def combinePairwiseAndSingleChainFeats(self, pairFeats, singleFeatsL, singleFeatsR):
    '''
      Merges pairFeats pandas.DataFrame with singleFeatsL and singleFeatsR dataFrames.
      singleFeatsL has n rows (as many as ligand residues) and singleFeatsR has m rows
      (as many as ligand residues), and pairsFeats has ~ n*m rows
      (some amino acids pairs might not be considered due to several reasons).
    
      @param pairFeats: pandas.DataFrame. A pandas.Dataframe in which each row represents the properties of 
                        a residue of the receptor
                        Column names are:
                        'chainIdL', 'structResIdL', 'resNameL', 'chainIdR', 'structResIdR', 'resNameR', 'categ' 
                         [propertiesP]    

      @param singleFeatsL: pandas.DataFrame. A pandas.Dataframe in which each row represents the properties of 
                           a residue of the ligand
                           Column names are:
                           'chainIdL', 'structResIdL', 'resNameL',[propertiesL] 
                           
      @param singleFeatsR: pandas.DataFrame. A pandas.Dataframe in which each row represents the properties of 
                           a residue of the receptor
                           Column names are:
                           'chainIdR', 'structResIdR', 'resNameR',[propertiesR]

      @return directPairs: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R).
                      Column names are:
                      'chainIdL', 'structResIdL', 'resNameL', 'chainIdR', 'structResIdR', 'resNameR', 'categ' 
                      [properties_P .... propertiesL     .... propertiesR] 
                         
    '''
    if singleFeatsL.shape[1]!= singleFeatsR.shape[1]:
      print( "singleFeatsL and singleFeatsR have different number of variables")
      featsLNo_chain= set( (elem[:-1] for elem in singleFeatsL.columns))
      featsRNo_chain= set( (elem[:-1] for elem in singleFeatsR.columns))
      print("L %d R %d"%(len(singleFeatsL.columns),len(singleFeatsR.columns)) ,"L dif R", 
            sorted(featsLNo_chain.difference(featsRNo_chain)), "R diff L",
            sorted(featsRNo_chain.difference(featsLNo_chain)) )
      raise ValueError( "singleFeatsL and singleFeatsR have different number of variables")
    otherPairColumns= set(pairFeats.columns.values).difference(set(AbstractProtocol.ARE_STR_TYPE_COLUMNS+["categ"]))
    if len(otherPairColumns)>0: #add _P suffix to all pairwise features
      pairFeats.rename(columns={elem:elem+"_P" for elem in otherPairColumns}, inplace=True)      
    directPairs= pd.merge(pairFeats,   singleFeatsL, how='inner', on=None)   
    directPairs= pd.merge(directPairs, singleFeatsR, how='inner', on=None)
    return directPairs

        
  def prepareDataForPairwiseAggregat(self, df):
    '''
      abstract method
    '''
    return None
    
  def addPairwiseAggregation(self, df):
    '''
      Adds environment pairwise features to a df that contains pairwise features (named featName_P)
      @param df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (ligand, receptor)
                      Column names are:
                      'chainIdL', 'structResIdL', 'resNameL', 'chainIdR', 'structResIdR', 'resNameR', 'categ' 
                      [properties_P .... propertiesL .... propertiesR] #no defined order for properties
      
      @return newDf: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R). New columns will have been added, all then
                      named %sAggr%d%s(numAggr/factorAggr, numberOfNewFeature, ligand/receptor)
                      Column names are:
                      'chainIdL', 'structResIdL', 'resNameL', 'chainIdR', 'structResIdR', 'resNameR', 'categ' 
                      [ properties_P propertiesR propertiesL ] #no defined order for properties
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

    numericAggreL2r= computeNumericAggr(pairwiseDf, l_2_r_neigs)
    numericAggreR2l= computeNumericAggr(pairwiseDf, r_2_l_neigs)
    numericAggreN2N= computeNumericAggr(pairwiseDf, l_neigs_2_r_neigs)
    numericAggreL2r.columns= [ "l2r-pair"+elem for elem in numericAggreL2r.columns]
    numericAggreR2l.columns= [ "r2l-pair"+elem for elem in numericAggreR2l.columns]
    numericAggreN2N.columns= [ "n2n-pair"+elem for elem in numericAggreN2N.columns]
    df= pd.concat([df,  numericAggreL2r, numericAggreR2l, numericAggreN2N ], axis=1)
    return df

def computeNumericAggr(df, selectedRows):
  '''
    Compute aggregation functions (min, max, mean and sum) for the rows of data. Each row is aggregated
    over the rows included in selectedRows
    @param data: pandas.DataFrame(nrow, nfeatures). The numeric data to aggregate. Each row i is averaged, sum... with
                                          the rows such that selectedRows[i,:]==True 
    @param selectedRows: [[]]: int (nrow, variableLength)): Each row i is averaged, with the rows included in 
                                                    selectedRows[i]    
    @return np.array (nrow, 4*nfeatures)  
  '''
#  s= time.time()
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
    over the rows included in selectedRows
    @param data: pandas.DataFrame(nrow, nfeatures). The numeric data to aggregate. Each row i is averaged,  with
                                          the rows such that selectedRows[i,:]==True 
    @param selectedRows: [[]]: int (nrow, variableLength)): Each row i is added with the rows included in 
                                                    selectedRows[i]
                           
    @return np.array (nrow, nfeatures)  
  '''
  data= df.values
#  aggregatedResults= -(2**16)* np.ones( (data.shape[0], data.shape[1]))
  aggregatedResults= -(2**10)* np.ones( (data.shape[0], data.shape[1]))
  for i in range(data.shape[0]):
    dataToAggregate= data[selectedRows[i],:]
    if dataToAggregate.shape[0]>0:
      aggregatedResults[i, :]= np.sum( dataToAggregate, axis=0)
  aggregatedResults= pd.DataFrame(aggregatedResults)
  aggregatedResults.columns= [ "factorAggr-"+name for name in df.columns] 
  return aggregatedResults
     
