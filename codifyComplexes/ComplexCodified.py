import random
import pandas as pd

class ComplexCodified(object):
  '''
    This class contains pairs of amino acids (a,b) such that a belongs
    to the ligand of the complex and b to the receptor of the complex.
    Each pair of amino acids is codified as a numerical vector. Factorial
    variables are codified as one hot encoding (dummy) variables
    All the pairs are represented as pandas.DataFrame and stored in self.pairsDirect
    for ligand features first  and self.pairsTranspose for receptor features first.
  '''
  def __init__(self, prefix, pairsDf):
    '''
      @param prefix: str. A prefix that identifies a protein complex
 
      @param pairsDf: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids.
                      Column names are:
                      'chainIdL', 'structResIdL', 'resNameL', 'chainIdR', 'structResIdR', 'resNameR', 'categ' ...
                       [ propertiesL     .... propertiesR ... propertiesP] #no defined order for properties

    '''
    self.prefix= prefix
    self.colNames= pairsDf.columns.values.tolist()
    self.pairsDirect= pairsDf
    self.categIndex= self.colNames.index("categ")    
    self.pairsTranspose= self._obtainTranspose(pairsDf)
    assert self.categIndex == self.pairsTranspose.columns.values.tolist().index("categ"), "Error direct and transpose categIndex are differents"
    assert list(self.pairsTranspose.columns)== list(self.pairsDirect.columns), "Error, direct pairs codification.columns!= transpose.columns"
    
  def getPrefix(self):
    '''
      returns the prefix that identifies the complex
      @return self.prefix: str
    '''
    return self.prefix
    
  def getLabels(self):
    '''
      returns a 1D-np.array with the labels (1->contact; -1 or 0 -> no contact)
      @return self.pairsDirect["categ"].values:  np.array (n,) where n is the number of pairs
    '''
    return self.pairsDirect["categ"].values

  def getIds(self):
    '''
      @return pd.DataFrame which the next columns:
            'chainIdL', 'structResIdL', 'resNameL', 'chainIdR', 'structResIdR', 'resNameR', 'categ' 
      categ values:  (1->contact; -1 or 0 -> no contact)            
    '''
    return self.pairsDirect.iloc[:,0:(self.categIndex+1)]
    
  def getData(self):
    '''
      returns codified pairs in direct form (L to R) and transpose form (R to L)
      Both them are returned as 2D-np.array
      @return (PairsDirect, PairsTranspose): (np.array (n,m), np.array (n,m)) where n is the number of pairs 
              and m is the number of features for each pair
    '''
    return (self.pairsDirect.iloc[:,(self.categIndex+1):].values, self.pairsTranspose.iloc[:,(self.categIndex+1):].values)
  
  def _obtainTranspose(self,df):
    '''
      given a pandas.DataFrame that contains the amino acid pairs in direct form (L to R),
      returns an equivalent pandas.DataFrame in transpose form (R to L)
      @param df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R).
                      Column names are:
                      'chainIdL', 'structResIdL', 'resNameL', 'chainIdR', 'structResIdR', 'resNameR', 'categ' ...
                       [propertiesP .... propertiesL     .... propertiesR]
      @return transposeDf: np.array (n,m) where n is the number of pairs and m is the number of features for each pair
    '''
    colNames= self.colNames
    categIndex= self.categIndex
    idsNames= colNames[:(categIndex+1)]
    lFeatNames= [elem for elem in colNames[(categIndex+1):] if elem[-1]=="L"]
    rFeatNames= [elem for elem in colNames[(categIndex+1):] if elem[-1]=="R"]
    pairwiseFeatNames= [elem for elem in colNames[(categIndex+1):] if elem.endswith("_P")]
    colNamesNew= idsNames+ rFeatNames+ lFeatNames+ pairwiseFeatNames
    transposeDf= df[ colNamesNew ].copy()
    transposeDf.columns= colNames
    return transposeDf
    
  def getSampledVersion(self, samplingFold= 1.0):
    '''
      return a sampled version of self. Thus, sampled version will have just n= num_posPairs*(1+ samplingFold)
      pairs: all the positive pairs and num_posPairs*(samplingFold) randomly selected negative pairs.
      @param samplingFold: Float. Number of times the number of negative sampled pairs is bigger than positive pairs. 
                                 numNegativePairs= samplingFold*numPositivePairs. (dealing with imbalanced data sets)
      @return sampledComplex: ComplexCodified object which contains just sampled pairs
    ''' 

    df= self.pairsDirect
#    df= df[ (df["total_RASAL"] >= 5.0) &  (df["total_RASAR"] >= 5.0) ] #samplig just accesible residues
        
    posIndices = [i for i, x in enumerate(df['categ']) if x == 1]
    negIndices = [i for i, x in enumerate(df['categ']) if x != 1]
    numNegToPick= int(samplingFold*len(posIndices))
    if numNegToPick ==0:
      ValueError("There are no positive pairs in data for %s"%(self.prefix))
    if len(negIndices) > numNegToPick:
      negIndices= random.sample(negIndices, numNegToPick)
    posDf= df.iloc[posIndices,:]
    negDf= df.iloc[negIndices,:]
    sampledDf= pd.concat([posDf, negDf])
    sampledDf = sampledDf.sample(frac=1, replace=False)
    sampledDf.reset_index(inplace=True, drop=True)
    sampledComplex= type(self)(self.prefix, sampledDf)
    return sampledComplex
    
  def __str__(self):
    return "ComplexCodified object: %s"%self.prefix
    

    
    
