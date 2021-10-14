import random
import pandas as pd
from codifyComplexes.codifyProtocols.AbstractProtocol import FEATURES_MISMATH_TOLERANCE
class ComplexCodified(object):
  '''
    This class contains pairs of amino acids (a,b) such that a belongs
    to the ligand of the complex and b to the receptor of the complex.
    Each pair of amino acids is codified as a numerical vector. Factorial
    variables are codified as one hot encoded (dummy) variables
    All the pairs are represented in a pandas.DataFrame and stored in self.pairsDirect
    for ligand features first and self.pairsTranspose for receptor features first
  '''
  def __init__(self, prefix, pairsDf, prefixesInvolvedInCoding):
    '''
      :param prefix: str. A prefix that identifies a protein complex
      :param pairsDf: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids.
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ' ...
                       [ propertiesL     .... propertiesR ... propertiesP] #no defined order for properties
      :param prefixesInvolvedInCoding: str[]. All prefixes that indirectly (feedback scores) have been involved
                                              in the codification of this protein.
    '''
    self.prefix= prefix
    self.colNames= pairsDf.columns.values.tolist()
    self.pairsDirect= pairsDf
    self.categIndex= self.colNames.index("categ")
    self.shape= self.pairsDirect.iloc[:, self.categIndex+1:].shape
    self.pairsTranspose= self._obtainTranspose(pairsDf)
    assert self.categIndex == self.pairsTranspose.columns.values.tolist().index("categ"), "Error direct and transpose categIndex are differents"
    assert list(self.pairsTranspose.columns)== list(self.pairsDirect.columns), "Error, direct pairs codification.columns!= transpose.columns"
    self.prefixesInvolvedInCoding= prefixesInvolvedInCoding



  def getPrefix(self):
    '''
      returns the prefix that identifies the complex
      :return self.prefix: str
    '''
    return self.prefix
    
  def getLabels(self):
    '''
      returns a 1D-np.array with the labels (1->contact; -1 or 0 -> no contact)
      :return self.pairsDirect["categ"].values:  np.array (n,) where n is the number of pairs
    '''
    return self.pairsDirect["categ"].values

  def getIds(self):
    '''
      :return pd.DataFrame which the next columns:
            'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ' 
      categ values:  (1->contact; -1 or 0 -> no contact)            
    '''
    return self.pairsDirect.iloc[:,0:(self.categIndex+1)]
    
  def getData(self):
    '''
      returns codified pairs in direct form (L to R) and transpose form (R to L)
      Both them are returned as 2D-np.array
      :return (PairsDirect, PairsTranspose): (np.array (n,m), np.array (n,m)) where n is the number of pairs
              and m is the number of features for each pair
    '''
    return (self.pairsDirect.iloc[:,(self.categIndex+1):].values, self.pairsTranspose.iloc[:,(self.categIndex+1):].values)
  
  def _obtainTranspose(self,df):
    '''
      given a pandas.DataFrame that contains the amino acid pairs in direct form (L to R),
      returns an equivalent pandas.DataFrame in transpose form (R to L)
      :param df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acids in direct form (L to R).
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ' ...
                       [propertiesP .... propertiesL     .... propertiesR]
      :return transposeDf: np.array (n,m) where n is the number of pairs and m is the number of features for each pair
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
      :param samplingFold: Float. Number of times the number of negative sampled pairs is bigger than positive pairs.
                                 numNegativePairs= samplingFold*numPositivePairs. (dealing with imbalanced data sets)
      :return sampledComplex: ComplexCodified object which contains just sampled pairs
    ''' 
    return type(self)(self.prefix, self._sample_oneDf( self.pairsDirect, samplingFold), self.prefixesInvolvedInCoding )


  def _sample_oneDf(self,df, samplingFold):
    '''
      :param df: pandas.Dataframe which contains all pairs of a complex
      :param samplingFold: Float. Number of times the number of negative sampled pairs is bigger than positive pairs.
                                 numNegativePairs= samplingFold*numPositivePairs. (dealing with imbalanced data sets)
      :return sampledDf: pandas.Dataframe which contains just sampled pairs
    '''
#    df= df[ (df["total_RASAL"] >= 5.0) &  (df["total_RASAR"] >= 5.0) ] #samplig just accesible residues   
    posIndices = [i for i, x in enumerate(df['categ']) if x == 1]
    negIndices = [i for i, x in enumerate(df['categ']) if x != 1]
    numNegToPick= int(samplingFold*len(posIndices))
    
    random_state= abs(hash(self.prefix.split("@")[0].split("#s")[0])) #ensure that all complexes with the same prefix are equally sampled. Required for results average
    if "@" in self.prefix: #sample different pairs in second step to get more variability
      random_state+= 121
    random_state= random_state // 2**32 -1
    
    if numNegToPick ==0:
      ValueError("There are no positive pairs in data for %s"%(self.prefix))
      
    if len(negIndices) > numNegToPick:
      random.seed(random_state)
      negIndices= random.sample(negIndices, numNegToPick)
      random.seed(None)
    posDf= df.iloc[posIndices,:]
    negDf= df.iloc[negIndices,:]
    sampledDf= pd.concat([posDf, negDf])
    sampledDf = sampledDf.sample(frac=1, replace=False, random_state= random_state)
    sampledDf.reset_index(inplace=True, drop=True)
    return sampledDf
    
  def __str__(self):
    return "ComplexCodified object: %s"%self.prefix
    


class ComplexSeqStructCodified( ComplexCodified):
  '''
    This class contains pairs of amino acids (a,b) such that a belongs
    to the ligand of the complex and b to the receptor of the complex.
    Each pair of amino acids is codified as a numerical vector. Factorial
    variables are codified as one hot encoded (dummy) variables
    All the pairs are represented in a pandas.DataFrame and stored in self.pairsSeqStruct
  '''
  
  def __init__(self, prefix, pairsDf_seqL=None, pairsDf_seqR=None, prefixesInvolvedInCoding= []):
    '''
      :param prefix: str. A prefix that identifies a protein complex
      :param pairsDf_seqL: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acid. If predicting, pairsDf_seqL is None if pairsDf_seqR not is None 
                      Column names are:
                      'chainIdL', 'resIdL', 'resNameL', 'chainIdR', 'resIdR', 'resNameR', 'categ' ...
                       [ propertiesSeq     .... propertiesStruct ... propertiesP] #no defined order for properties
      :param pairsDf_seqR: pandas.DataFrame. A pandas.Dataframe in which each row represents
                      a pair of amino acid. If predicting, pairsDf_seqL is None if pairsDf_seqR not is None 
                      [ propertiesSeq     .... propertiesStruct ... propertiesP] #no defined order for properties                      
      :param prefixesInvolvedInCoding: str[]. All prefixes that indirectly (feedback scores) have been involved
                                              in the codification of this protein.
    '''
    assert not (pairsDf_seqL is None and pairsDf_seqR is None), "Error, at least one should be provided"
    self.prefix= prefix
    if not isinstance(pairsDf_seqL, type(None)):
      self.colNames= pairsDf_seqL.columns.values.tolist()
      self.pairsDirect= pairsDf_seqL
    else:
      self.colNames= pairsDf_seqR.columns.values.tolist()
      self.pairsDirect= pairsDf_seqR
    self.pairsDf_seqL= pairsDf_seqL
    self.categIndex= self.colNames.index("categ")
    self.pairsDf_seqR= pairsDf_seqR
    self.shape= self.pairsDirect.iloc[:, self.categIndex+1:].shape
    self.prefixesInvolvedInCoding= prefixesInvolvedInCoding


    #For the case of mixed protocol training, make sure that the two dataframes have equal shape
    if not isinstance(self.pairsDf_seqL, type(None)) and not isinstance(self.pairsDf_seqR, type(None)):
      if self.pairsDf_seqL.shape != self.pairsDf_seqR.shape:
        idsL= self.pairsDf_seqL.iloc[:, :(self.categIndex+1)]
        idsR= self.pairsDf_seqR.iloc[:, :(self.categIndex+1)]
        idsL['indexL'] = idsL.index
        idsR['indexR'] = idsR.index
        shared_ids= pd.merge(  idsL, idsR, how="outer", indicator=True,
                              on=["chainIdL",  "resIdL","resNameL","chainIdR","resIdR","resNameR", "categ"])
        nonSharedIds= shared_ids[ shared_ids["_merge"]!="both"]

        mismatchFrac= float(nonSharedIds.shape[0]) / shared_ids.shape[0]
        del shared_ids
        if mismatchFrac > FEATURES_MISMATH_TOLERANCE:
          raise ValueError("Too much mismatch between sequence only codification in receptor and ligand for %s"%self.prefix)
        remove_fromL= nonSharedIds[nonSharedIds["_merge"]=="left_only"]["indexL"].values
        self.pairsDf_seqL= self.pairsDf_seqL.drop(remove_fromL, axis=0)
        remove_fromR= nonSharedIds[nonSharedIds["_merge"]=="right_only"]["indexR"].values
        self.pairsDf_seqR= self.pairsDf_seqR.drop(remove_fromR, axis=0)

  def getData(self):
    '''
      returns codified pairs in ligand_is_seqOnly) and  receptor_is_seqOnly form
      Both them are returned as 2D-np.array
      :return (ligandIsSeqOnly, receptorIsSeqOnly): (np.array (n,m), np.array (n,m)) where n is the number of pairs
              and m is the number of features for each pair
    '''
    pairs_seqL= self.pairsDf_seqL.iloc[:,(self.categIndex+1):].values if not isinstance(self.pairsDf_seqL, type(None)) else None
    pairs_seqR= self.pairsDf_seqR.iloc[:,(self.categIndex+1):].values if not isinstance(self.pairsDf_seqR, type(None)) else None
    return (pairs_seqL, pairs_seqR)

  def getSampledVersion(self, samplingFold= 1.0):
    '''
      return a sampled version of self. Thus, sampled version will have just n= num_posPairs*(1+ samplingFold)
      pairs: all the positive pairs and num_posPairs*(samplingFold) randomly selected negative pairs.
      :param samplingFold: Float. Number of times the number of negative sampled pairs is bigger than positive pairs.
                                 numNegativePairs= samplingFold*numPositivePairs. (dealing with imbalanced data sets)
      :return sampledComplex: ComplexCodified object which contains just sampled pairs
    ''' 
    sampled_seqL= None if isinstance(self.pairsDf_seqL, type(None)) else self._sample_oneDf(self.pairsDf_seqL, samplingFold)
    sampled_seqR= None if isinstance(self.pairsDf_seqR, type(None)) else self._sample_oneDf(self.pairsDf_seqR, samplingFold)
    return type(self)(self.prefix, sampled_seqL, sampled_seqR, self.prefixesInvolvedInCoding )
    
  def __str__(self):
    return "ComplexSeqStructCodified object: %s"%self.prefix    
    
