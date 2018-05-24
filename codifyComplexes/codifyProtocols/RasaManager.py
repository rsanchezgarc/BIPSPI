import sys, os
from Config import Configuration
from .DataLoaderClass import  DataLoader

FEATURES_TO_INCLUDE= [
  ("psaia", ("structStep/PSAIA/procPSAIA", [8], {"total_RASA":8}))
]

class RasaManager(DataLoader):
  def __init__(self, dataRootPath=None, featuresToInclude=None):
    if dataRootPath is None:
      dataRootPath= Configuration().computedFeatsRootDir
    if featuresToInclude is None:
      featuresToInclude= FEATURES_TO_INCLUDE
    else:
      if "psaia" not in zip(* featuresToInclude)[0]:
        featuresToInclude.insert(0, ("psaia", ("structStep/PSAIA/procPSAIA", [8], {"total_RASA":8})) )
    self.psaiaIndex= zip(* featuresToInclude)[0].index("psaia" )
    DataLoader.__init__(self, dataRootPath, featuresToInclude)
    self.rasaLDict={}
    self.rasaRDict={}    
  def loadRasa(self, prefixL, prefixR):
    '''
      loads relative asa of both receptor and ligand. To do so, it fills self.rasaLDict and
      self.rasaRDict dictionaries
      @param prefixL:str. Identifier of the ligand (1A2K_l)
      @param prefixR:str. Identifier of the receptor (1A2K_r)
    '''  
    __, (fnamesIterator, selectedCols)= self.getParamsForLoadingFile( prefixL, self.psaiaIndex, useNameColum="total_RASA")
    rasaData= self.loadDataFile( fnamesIterator, selectedCols=selectedCols).to_dict(orient="list")
    self.rasaLDict= { (chainId, resId): rasa for chainId, resId, rasa in 
                       zip(rasaData["chainId"], rasaData["structResId"], rasaData["total_RASA"])}

    __, (fnamesIterator, selectedCols)= self.getParamsForLoadingFile( prefixR, self.psaiaIndex, useNameColum="total_RASA")
    rasaData= self.loadDataFile( fnamesIterator, selectedCols=selectedCols).to_dict(orient="list")
    self.rasaRDict= { (chainId, resId): rasa for chainId, resId, rasa in 
                       zip(rasaData["chainId"], rasaData["structResId"], rasaData["total_RASA"])}
    
  def getRasaL(self, chain, resId):
    '''
      returns relative asa of residue with chain and resId in ligand
      @param chain:str. Chain id
      @param resId:str. resId (as found in pdb file)
      @return None if no rasa available for asked residue or float otherwise
    '''  
    if (chain, resId) in self.rasaLDict:
      return self.rasaLDict[(chain, resId)]
    else:
      return None
      
  def getRasaR(self, chain,resId):
    '''
      returns relative asa of residue with chain and resId in receptor
      @param chain:str. Chain id
      @param resId:str. resId (as found in pdb file)
      @return None if no rasa available for asked residue or float otherwise
    '''
    if (chain, resId) in self.rasaRDict:
      return self.rasaRDict[(chain, resId)]
    else:
      return None
