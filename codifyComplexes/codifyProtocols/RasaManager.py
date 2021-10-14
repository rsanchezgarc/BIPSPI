import sys, os
from Config import Configuration
from .DataLoaderClass import  DataLoader

class RasaManager(DataLoader):
  ''' 
  This class is used to query the relative asa of the residues
  '''
  
  RASA_FEAT_DESCR= ("psaia", ("structStep/PSAIA/procPSAIA", [8])) 
  def __init__(self, dataRootPath=None, verbose=False):


    if not hasattr(self, "dataRootPath"):
      DataLoader.__init__(self, dataRootPath, verbose)
    else:
      DataLoader.__init__(self, self.dataRootPath, verbose)

    self.rasaLDict = {}
    self.rasaRDict = {}

  def loadRasa(self, prefixL, prefixR):
    '''
      loads relative asa of both receptor and ligand. To do so, it fills self.rasaLDict and
      self.rasaRDict dictionaries
      :param prefixL:str. Identifier of the ligand (1A2K_l) or None if no structural information available
      :param prefixR:str. Identifier of the receptor (1A2K_r) or None if no structural information available
    '''
    if prefixL:
      __, (fnamesIterator, selectedCols)= self.getParamsForLoadingFile( prefixL, RasaManager.RASA_FEAT_DESCR, useNameColum="total_RASA")
      rasaData= self.loadDataFile( fnamesIterator, selectedCols=selectedCols).to_dict(orient="list")
      self.rasaLDict= { (chainId, resId): rasa for chainId, resId, rasa in 
                         zip(rasaData["chainId"], rasaData["resId"], rasaData["total_RASA"])}

    if prefixR:
      __, (fnamesIterator, selectedCols)= self.getParamsForLoadingFile( prefixR, RasaManager.RASA_FEAT_DESCR, useNameColum="total_RASA")
      rasaData= self.loadDataFile( fnamesIterator, selectedCols=selectedCols).to_dict(orient="list")
      self.rasaRDict= { (chainId, resId): rasa for chainId, resId, rasa in 
                         zip(rasaData["chainId"], rasaData["resId"], rasaData["total_RASA"])}

  def getRasaL(self, chain, resId):
    '''
      returns relative asa of residue with chain and resId in ligand
      :param chain:str. Chain id
      :param resId:str. resId
      :return None if no rasa available for asked residue or float otherwise
    '''  
    if (chain, resId) in self.rasaLDict:
      return self.rasaLDict[(chain, resId)]
    else:
      return None
      
  def getRasaR(self, chain,resId):
    '''
      returns relative asa of residue with chain and resId in receptor
      :param chain:str. Chain id
      :param resId:str. resId
      :return None if no rasa available for asked residue or float otherwise
    '''
    if (chain, resId) in self.rasaRDict:
      return self.rasaRDict[(chain, resId)]
    else:
      return None
