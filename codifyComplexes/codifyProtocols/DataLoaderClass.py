import os
import pandas as pd
import numpy as np

from codifyComplexes.CodifyComplexException import CodifyComplexException
from Config import Configuration

class DataLoader(Configuration):
  '''
    This class implements loading file functions
  '''
  
  ARE_STR_TYPE_COLUMNS= ["chainId", "structResId", "chainIdL", "structResIdL", "resNameL", "chainIdR", 
                            "structResIdR", "resNameR", "resName" ]
  ALWAYS_SKIP_COLUMNS=["seqIndex"]
  IGNORE_X_AA= False
  
  def __init__(self, dataRootPath, singleChainfeatsToInclude, pairfeatsToInclude=None, verbose=False):
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

      @param singleChainfeatsToInclude: dict. Dictionary that contains the paths where features needed for complex codification 
                                        that describre single chain properties are
                                      located. Must have the following format:
                                      {"featName":(relativePath_from_dataRootPath, listOfColumnNumbers, dictForNamedColums)}
                                      dictForNamedColums= {"myFeatName":colNumber}
      @param pairfeatsToInclude: dict. Dictionary that contains the paths where features needed for complex codification 
                                        that describre single chain properties are
                                      located. Must have the following format:
                                      {"featName":(relativePath_from_dataRootPath, listOfColumnNumbers, dictForNamedColums)}
                                      dictForNamedColums= {"myFeatName":colNumber}
      @param verbose: bool.
    '''  
    Configuration.__init__(self)
    self.dataRootPath= dataRootPath
    self.verbose= verbose
    self.singleChainfeatsToInclude= singleChainfeatsToInclude
    self.pairfeatsToInclude= None if pairfeatsToInclude is None else pairfeatsToInclude[:]

  def getParamsForLoadingFile(self, prefixOneChainType, featNum, useNameColum=None, lookForPairFeats=False):
    '''
      Given a ligand or receptor id (prefixOneChainType) and a feature_name(featName) which is a key 
      in self.singleChainfeatsToInclude dictionary, returns an iterator over all files involved in that 
      feature and the column numbers that where selected. By default, selected columns are stored
      in self.singleChainfeatsToInclude, but if useNameColumn is specified, selectedCols will be a one 
      element list with the index of the column whose name matches useNameColumn
      @param prefixOneChainType: str. A prefix that identifies the receptor or the ligand of a complex
      @param featNum: int. indice of singleChainfeatsToInclude list that contains information of features
      @param useNameColum: str. If None, selectedCols will be obtained from self.singleChainfeatsToInclude.
                                Otherwise, selectedCols= [ colNum ], where colNum is the column number
                                of the column whose name matches useNameColum
      @param lookForPairFeats: bool. Whether to load pairwise features or single chain features otherwise 
      @return (fnamesIterator, selectedCols)
            fnamesIterator: Iterator str. Fnames that contains feature featName for ligand or receptor id
                                           prefixOneChainType
            selectedCols: int[]. Columns that will to be selected                                           
    '''
##    print(self.singleChainfeatsToInclude, featNum)
    if lookForPairFeats:
      featName, (relativePathToFeats, selectedCols, nameColDict) = self.pairfeatsToInclude[featNum]
    else:
      featName, (relativePathToFeats, selectedCols, nameColDict) = self.singleChainfeatsToInclude[featNum]
    if not useNameColum is None: selectedCols= [nameColDict[useNameColum]]
    featuresPath= os.path.join(self.dataRootPath, relativePathToFeats)
    return featName, (self.getFullNamesIterInPath(prefixOneChainType, featuresPath), selectedCols)
        
  def loadDataFile(self, fnamesIterator, selectedCols=None):
    '''
      Returns the pandas.DataFrame resulted from all fnames files concatenation
      contained at fnamesIterator
      @param fname: str. Iterator of file names containing features to be loaded. 
      @param selectedCols: int[]. List of colums to be selected. If None
                                  all colums will be selected
      @return  df: pandas.DataFrame. A pandas.Dataframe in which each row represents the feature of 
                           a residue that is contained in one of the file names of fnamesIterator
                           Column names are specified in the header of files

    '''
    resultDF_list= []
    factorColumns={}
    fnamesList= list(fnamesIterator)
    for fname in fnamesList:
      try:
        df= pd.read_table(fname,sep='\s+', header='infer', comment="#", 
                          dtype= {elem:str for elem in DataLoader.ARE_STR_TYPE_COLUMNS})
        if DataLoader.IGNORE_X_AA:
          if "resName" in df:
            df= df[df["resName"]!="X"]
          if "resNameL" in df:
            df= df[(df["resNameL"]!="X") & (df["resNameR"]!="X")]
      except pd.io.common.CParserError:
        print("Error reading %s"%fname)
        raise
      df_toCheck= df if not "categ" in df else df.drop("categ",axis=1)
##      print(df_toCheck.head())
      if df_toCheck.isnull().values.any()==True:
        raise ValueError("There are missing values or nans in %s"%fname)
      with open(fname) as f:
        firstLine= f.readline()
        if firstLine.startswith("#Levels"):
          factorColumns=  self._parseLevels(firstLine, df.columns)
      for colName in factorColumns:
          df[colName]=  pd.Categorical(df[colName], categories=factorColumns[colName], ordered=False) 
      #Select desired columns    
      if not selectedCols is None:
        ids_cols=[]
        for i,colum in enumerate(df.columns):
          if colum in DataLoader.ARE_STR_TYPE_COLUMNS:
            ids_cols.append(i)
        assert len(set(ids_cols).intersection(set(selectedCols))) ==0, ("problems loading file %s. Ids_cols and "+
                                                                  "selectedCols overlap")%fname                                                
        df= df.ix[:, ids_cols+ selectedCols]
      for colName in DataLoader.ALWAYS_SKIP_COLUMNS:
        if colName in df:
          del df[colName]

      df= pd.get_dummies(df, prefix_sep='_dummy_', columns= list(factorColumns.keys()))
      resultDF_list.append( df )
    df= pd.concat(resultDF_list)
    assert df.shape[0]> 1, "Error loading files %s there are no rows"%str(fnamesList)
    return df
 
  def _parseLevels(self, firstLine, colNames):
    '''
      Helper method for self.loadDataFile. Itparses first line of a file
      looking for factor variables levels
      @param firstLine: str. The first line of a file
      @param colNames: str[]. The column names of the pandas.DataFrame
      @return colNumDict: { colName:[level0, level1 ...]}
    '''  
    colNumDict= [ factorDes.split(":") for factorDes in firstLine.split()[1:] ]
    colNumDict= { colNames[int(factorDes[0])]:factorDes[1].split(";") for factorDes in colNumDict }
    return colNumDict
      
  def getFullNamesIterInPath(self, prefix, dirname):
    '''
      returns a full filename that startswith given prefix an belongs to 
      directory dirname
      @param prefix: str. the prefixes with which file names must start to be considered
      @param dirname: str. The directory where files will be looked for
      @yields fullName. The full name of a file that startswith prefix and is in dirname
    '''
    fullName= None
    for fname in os.listdir(dirname):
      if fname.startswith(prefix):
        fullName= os.path.join(dirname, fname)
        yield fullName

