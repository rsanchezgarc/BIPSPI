from __future__ import absolute_import
import os, sys
from Bio.PDB.DSSP import DSSP
from .structToolManager import StructToolManager
from utils import myMakeDir, tryToRemove

class DsspComputer(StructToolManager):
  '''
  Extends StructToolManager class. Computes DSSP for a given complex
  '''
  VAR_LIST= ['2ndStruct']
        
  def __init__(self, computedFeatsRootDir, statusManager= None):
    '''
      :param computedFeatsRootDir: str. path where features will be stored. If None they will be stored
                                        at default path (assigned in ../Config.py)
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''

    StructToolManager.__init__(self, computedFeatsRootDir, statusManager)
    self.dsspBinPath= os.path.join(self.dsspRootDir,"mkdssp" )
    self.outPath=  myMakeDir(self.computedFeatsRootDir, "DSSP")

  def getFNames(self, extendedPrefix):
    '''
      Returns a dict that contains the fnames that will be used for hse
      :param extendedPrefix. prefix for output fnames.
      :return list of fnames: [ fname1, fnam2, ...]
    '''
    prefix, chainType = self.splitExtendedPrefix(extendedPrefix )[:2]
    return [os.path.join(self.outPath, prefix+"_"+chainType+"_.dssp.tab.gz")]
                                          

  def computeOneFile(self, pdbFname, struct):
    '''
      Computes DSSP for a given pdb file
      :param pdbFname: str. fname to pdb file
      :param struct: Bio.PDB.Structure
    '''
    
    allResidues= set([])
    for chain in struct[0]:
      residues= chain.get_residues()
      allResidues= allResidues.union(set(residues))
      
    prefixExtended= self.getExtendedPrefix(pdbFname)
    prefix, chainType= self.splitExtendedPrefix(prefixExtended)[:2]
    if self.checkAlreayComputed(prefixExtended):
      print("Dssp already computed for %s"%prefixExtended)
      return 0
    print("launching Dssp over %s"%prefixExtended)
    try:
    
      featuresDict= {}
      try:
        dssp_out = DSSP(struct[0], pdbFname, dssp=self.dsspBinPath)
      except Exception as e:
        if "DSSP failed to produce an output" in e.message:
          dssp_out={}
        else:
          print(e)
          raise e
      for chainId, resId in dssp_out.keys():
        secStruct= dssp_out[(chainId, resId)][2]
        if secStruct=="-": secStruct="Z"
        try:
          featuresDict[ struct[0][chainId][resId] ]= secStruct
        except KeyError:
          continue

      dataDict={}
      for aa in allResidues:
        chainId_resIdStr_resName= self.fromRes2ChainResIdAndName(aa)
        chainId, resIdStr, resName= chainId_resIdStr_resName
        if resName is None: continue
        if aa in featuresDict:
          values= [featuresDict[aa]]
        else:
          values= ["Z"]

        record= [ chainId, resIdStr, resName] + values
        record= " ".join(record)

        try:
          dataDict[chainId].append(record)
        except KeyError:
          dataDict[chainId]=[record]
      categoricalLevels= {("H","B","E","G","I","T","S","Z"):("2ndStruct",)}
      self.writeResultsFromDataDictSingleChain( dataDict, outName= self.getFNames(prefixExtended)[0], 
                                                categoricalLevels=categoricalLevels)
        
    except (Exception, KeyboardInterrupt):
      self.tryToRemoveAllFnames(prefixExtended)
      raise


