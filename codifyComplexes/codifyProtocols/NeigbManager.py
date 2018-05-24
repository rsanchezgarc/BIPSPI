import sys, os
from Config import Configuration
from .RasaManager import RasaManager


ENVIRONMENT_MAP="structStep/VORONOI/VORONOI_30"
#ENVIRONMENT_MAP="structStep/VORONOI/VORONOI_10"
rASA_THR= 5

class NeigbManager( RasaManager):
  def __init__(self, dataRootPath=None, featuresToInclude=None, environPath= ENVIRONMENT_MAP):
    if not dataRootPath is None:
      self.computedFeatsRootDir= dataRootPath
    self.environ_path=environPath
    RasaManager.__init__(self, self.computedFeatsRootDir, featuresToInclude)
    
  def loadEnvirons(self, prefixL, prefixR, asaThr= rASA_THR):
    '''
      Loads neighbour residues from files contained at self.environ_path. Updates self.environL
      and self.environR
      @param prefixL:str. Identifier of the ligand
      @param prefixR:str. Identifier of the receptor.
      @param asaThr: float. relative asa threshold. If any of the aminoacids contained in
                     each row of the file are less accesible than asaThr, they won't
                     be cosidered as neighbours
    '''
    if self.rasaLDict=={} or self.rasaLDict=={}:
      self.loadRasa(prefixL, prefixR)
    for chainType in ["l", "r"]:
      if chainType=="l":
        getRasa= self.getRasaL
        prefix= prefixL
      elif chainType=="r":
        getRasa= self.getRasaR
        prefix= prefixR        
      else:
        raise ValueError("Error, chainType %s must be either 'l' or 'r'"%chainType)
      fname=os.path.join(self.dataRootPath, self.environ_path, prefix+"_u.voro")
      res={}  
      with open(fname) as f:
        for line in f:
          lineArray= line.split()
          resInd1, chain1 = lineArray[0].split("_")
          resInd2, chain2 = lineArray[1].split("_")

          asaRes1 = getRasa(chain1,resInd1)
          asaRes2 = getRasa(chain2, resInd2)

          if asaRes1==None or asaRes2==None:
            if self.verbose: print("No rasa available for ",resInd1, chain1,resInd2, chain2)
            continue
          if chain1 not in res:
            res[chain1]={}
          if chain2 not in res:
            res[chain2]={}
          if resInd1 not in res[chain1]:
            res[chain1][resInd1]=set([])
          if resInd2 not in res[chain2]:
            res[chain2][resInd2]=set([])
          if asaRes2 > asaThr:
            res[chain1][resInd1].add((chain2,resInd2))
          if asaRes1 > asaThr:
            res[chain2][resInd2].add((chain1,resInd1))
      if chainType=="l":          
        self.environL= res
      else:
        self.environR= res

    
  def getNeigs2(self, chainL, resIdL, chainR, resIdR):
    '''
      Returns residues that are neighbours to a given ligand residue and to
      a given receptor residue
      returns relative asa of residue with chain and resId in receptor
      @param chainL:str. Chain id of ligand residue
      @param resIdL:str. resId (as found in pdb file) of ligand residue
      @param chainR:str. Chain id of receptor residue
      @param resIdR:str. resId (as found in pdb file) of receptor residue      
      @return setL, setR. Set of residues that are neighbours. Each neighbour
                        residue is represented as tuple(chainId, resId)

    '''
    neigsL={}
    neigsR={}
    try:
      neigsL= self.environL[chainL][resIdL]
    except KeyError:
      pass
    try:
      neigsR= self.environR[chainR][resIdR]
    except KeyError:
      pass
    return neigsL, neigsR


  def getNeigs(self, chainId, resId, chainType):
    '''
      Returns residues that are neighbours to a given ligand residue and to
      a given receptor residue
      returns relative asa of residue with chain and resId in receptor
      @param chainL:str. Chain id of residue
      @param resIdL:str. resId (as found in pdb file) of residue
      @param chainType: str. "l" for ligand "r" for receptor"

    '''
    if chainType=="l" or chainType=="L":
      try:
        return self.environL[chainId][resId]
      except KeyError:
        return None
    else:
      try:
        return self.environR[chainId][resId]
      except KeyError:
        return None
