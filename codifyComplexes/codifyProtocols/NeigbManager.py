import os
from .RasaManager import RasaManager
from joblib import load

NEIGS_PATH="common/voroNeigs"
rASA_THR= 5

class NeigbManager( RasaManager):
  def __init__(self, dataRootPath=None, neigsPath= NEIGS_PATH):

    RasaManager.__init__(self, dataRootPath)

    self.neigsPath=  os.path.join(self.dataRootPath, neigsPath)
    self.neigsL={}
    self.neigsR={}
    
  def loadEnvirons(self, prefix, ligandHas=True, receptorHas=True, asaThr= rASA_THR):
    '''
      Loads neighbour residues from files contained at self.neigsPath. Updates self.neigsL
      and self.neigsR
      :param prefix: str. A prefix that identifies a complex
      :param ligandHas: bool. True if there is structural neighbourhood for ligand
      :param receptorHas: bool.  True if there is structural neighbourhood for receptor
      :param asaThr: float. relative asa threshold. If any of the aminoacids contained in
                     each row of the file are less accesible (relative) than asaThr, they won't
                     be cosidered as neighbours
    '''
    
    __, (fnameIter, __)=  self.getParamsForLoadingFile( prefix, featTuple=("voro", (self.neigsPath, None)),
                                               useNameColum=None )
    fname= list(fnameIter)[0]
    neigsDicts= load(fname)
    chainsWithNeigs= []
    if ligandHas: 
      chainsWithNeigs.append("l")
      prefixL= prefix+"_l"
    else:
      prefixL=None
      
    if receptorHas: 
      chainsWithNeigs.append("r")
      prefixR= prefix+"_r"
    else:
      prefixR=None

    self.neigsL= neigsDicts["neigsDictL"]
    self.neigsR= neigsDicts["neigsDictR"]

    if self.rasaLDict=={} or self.rasaLDict=={}:
      self.loadRasa(prefixL, prefixR)


    for chainType in chainsWithNeigs:
      if chainType=="l":
        getRasa= self.getRasaL
        neigsDict= self.neigsL
      elif chainType=="r":
        getRasa= self.getRasaR
        neigsDict= self.neigsR     
      else:
        raise ValueError("Error, chainType %s must be either 'l' or 'r'"%chainType)
      if neigsDict is None: continue
      idsToRemoveKeys=[]
      idsToRemoveValues= {}
      for chain_resId_1 in neigsDict:
        asaRes1 = getRasa( * chain_resId_1)
        if asaRes1 is None or asaRes1<  asaThr:
          idsToRemoveKeys.append( chain_resId_1 )
          continue
        for chain_resId_2 in neigsDict[chain_resId_1]:
          asaRes2 = getRasa( * chain_resId_2)
          if asaRes2 is None or asaRes2<  asaThr:
            try:
              idsToRemoveValues[chain_resId_1].add( chain_resId_2)
            except KeyError:
              idsToRemoveValues[chain_resId_1]=set([chain_resId_2])
              
      for idKey in  idsToRemoveKeys:
        del neigsDict[idKey]
      for idKey in idsToRemoveValues:  
        for idVal in idsToRemoveValues[idKey]:
          neigsDict[idKey].remove(idVal)
          
  def getNeigs2(self, chainL, resIdL, chainR, resIdR):
    '''
      Returns residues that are neighbours to a given ligand residue and to
      a given receptor residue
      returns relative asa of residue with chain and resId in receptor
      :param chainL:str. Chain id of ligand residue
      :param resIdL:str. resId (as found in pdb file) of ligand residue
      :param chainR:str. Chain id of receptor residue
      :param resIdR:str. resId (as found in pdb file) of receptor residue
      :return setL, setR. Set of residues that are neighbours. Each neighbour
                        residue is represented as tuple(chainId, resId)

    '''
    neigsL={}
    neigsR={}
    try:
      neigsL= self.neigsL[ (chainL, resIdL) ]
    except (KeyError, TypeError):
      pass
    try:
      neigsR= self.neigsR[ (chainR, resIdR) ]
    except (KeyError, TypeError):
      pass
    return neigsL, neigsR


  def getNeigs(self, chainId, resId, chainType):
    '''
      Returns residues that are neighbours to a given ligand residue and to
      a given receptor residue
      returns relative asa of residue with chain and resId in receptor
      :param chainL:str. Chain id of residue
      :param resIdL:str. resId (as found in pdb file) of residue
      :param chainType: str. "l" for ligand "r" for receptor"

    '''
    if chainType=="l" or chainType=="L":
      try:
        return self.neigsL[ (chainId, resId) ]
      except (KeyError, TypeError):
        return None
    else:
      try:
        return self.neigsR[ (chainId, resId) ]
      except (KeyError, TypeError):
        return None


