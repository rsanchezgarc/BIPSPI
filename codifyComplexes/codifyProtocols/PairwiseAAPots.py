import sys, os
import numpy as np
from Config import Configuration

class PairwiseAAindex(Configuration):
  '''
    This class deals with AAIndex pairwise potentials and makes easy to get 
    the aaIndex potential involving 2 amino acids by giving their one letter
    code
  '''
  badValue= 1024
  def __init__(self, data_path= None):
    '''
      @param data_path: str: Path where AAIndex files are located
    '''
    Configuration.__init__(self)  # Load configuration parameters such as path to programs

    self.protein_proteinIndexes=["KESO980101","KESO980102","MOOG990101"]
    if data_path is None:
      self.data_path= self.AAindexPath
    else:
      self.data_path=data_path
    self.data=self.load()
    
  def load(self):
    '''
      loads all aaIndex potentials contained in self.protein_proteinIndexes and
      returns them as a dictionary
      return res: Dict.   res[aaCodeL][aaCodeR]= [pot0, pot1 ...]
    '''  
    res={}
    for indexName in self.protein_proteinIndexes:
      fname=os.path.join(self.data_path,indexName)
      fname= os.path.realpath(fname)
      f=open(fname)
      aaTypes=f.readline()
      aaTypes=aaTypes.split()[1]
      aaTypes=list(aaTypes)

      for iLetter,line in zip(aaTypes,f):
        lineArray= line.split()
        if iLetter not in res:
          res[iLetter]={}
        for jLetter,score in zip(aaTypes,lineArray):
          if jLetter not in res[iLetter]:
            res[iLetter][jLetter]=[]
          res[iLetter][jLetter].append(float(score))
    
      f.close()
    return res
  def getScore(self,aa1,aa2):
    '''
      given aa1 and aa2, that are one letter amino acid codes, return the 
      associated aaIndex values.
      @param aa1: str. one letter amino acid code
      @param aa2: str. one letter amino acid code      
      
      @return self.data[aa1][aa2]: list containing the associated values
    '''
    try:
      values = self.data[aa1][aa2]
    except KeyError:
#      print ("Non standard aa")
      values= [PairwiseAAindex.BadValue for i in xrange(len(self.protein_proteinIndexes))]
    return values
    
  def addAAIndexToDF(self, df):
    '''
      adds to a pandas.DataFrame that represents amino acid pairs codification, new columns added
      using the aaIndex values.
      @param df: pandas.DataFrame. A pandas.Dataframe in which each row represents
                a pair of amino acids in direct form (L to R).
                Column names are:
                'chainIdL', 'structResIdL', 'resNameL', 'chainIdR', 'structResIdR', 'resNameR', 'categ' ...
                 [propertiesP .... propertiesL     .... propertiesR] #no defined order for properties
      @return df: pandas.DataFrame. The same pandas.DataFrame given as input but with new columns added
                                    whose name is "aaIndex%d"%i 
    '''

    getScore= self.getScore
    nPots= len(self.protein_proteinIndexes)

    resNames= zip(df["resNameL"],df["resNameR"])
    scoresNp= np.zeros((df.shape[0], 3))
    for i in range(df.shape[0]):
      scoresNp[i,:]= getScore( *resNames[i])
    for i in range(nPots):
      df.loc[:,"aaIndex%d_P"%i]= scoresNp[:, i]
    return df
    
if __name__== "__main__":
  mapa=PairwiseAAindex().data
  print mapa["F"]["S"]

