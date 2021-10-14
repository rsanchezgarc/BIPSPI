from __future__ import absolute_import, print_function
import os
from .structToolManager import StructToolManager
from .hSEExposureFixed import HSExposureCA, HSExposureCB, ExposureCN
from utils import myMakeDir

class HalfSphereComputer(StructToolManager):
  '''
  Extends StructToolManager class. Computes half sphere exposure for a given complex
  '''
  VAR_LIST= ['HSExposureCA1', 'HSExposureCA2', 'HSExposureCA3', 'HSExposureCB1', 'HSExposureCB2',
             'HSExposureCB3', 'ExposureCN']
  
                 
  def __init__(self, computedFeatsRootDir, statusManager= None):
    '''
      :param computedFeatsRootDir: str. path where features will be stored. If None they will be stored
                                        at default path (assigned in ../Config.py)
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''

    StructToolManager.__init__(self, computedFeatsRootDir, statusManager)
    self.outPath=  myMakeDir(self.computedFeatsRootDir, "halfSphereExpos")

  def getFNames(self, extendedPrefix):
    '''
      Returns a dict that contains the fnames that will be used for hse
      :param extendedPrefix. prefix for output fnames.
      :return list of fnames: [ fname1, fnam2, ...]
    '''
    prefix, chainType = self.splitExtendedPrefix(extendedPrefix )[:2]
    return [os.path.join(self.outPath, prefix+"_"+chainType+"_.hse.tab.gz")]
                                          

  def computeOneFile(self, pdbFname, struct):
    '''
      Computes HalfSphere for a given pdb file
      :param pdbFname: str. fname to pdb file
      :param struct: ignored
    '''
    
    allResidues= set([])
    for chain in struct[0]:
      residues= chain.get_residues()
      allResidues= allResidues.union(set(residues))
      
    prefixExtended= self.getExtendedPrefix(pdbFname)
    prefix, chainType= self.splitExtendedPrefix(prefixExtended)[:2]
    if self.checkAlreayComputed(prefixExtended):
      print("HalfSphere already computed for %s"%prefixExtended)
      return 0
    print("launching HalfSphere over %s"%prefixExtended)
    try:
      model = struct[0]
      featuresDict= {}
      for aa, feat in HSExposureCA(model):
        featuresDict[aa]= [elem if elem!= None else -1 for elem in feat]
      hseDict={ aa: [elem if elem!= None else -1 for elem in feat] for aa,feat in HSExposureCB(model) }

      for aa in set(featuresDict.keys()).union(set(hseDict.keys())):
        try:
          prevFeatures= featuresDict[aa]
        except KeyError:
          prevFeatures= [-1, -1, -1.0]
        try:
          newFeatures= hseDict[aa]
        except KeyError:
          newFeatures= [-1, -1, -1.0]      
        featuresDict[aa]= prevFeatures + newFeatures  
      hseDict={ aa: [feat] if feat!= None else -1 for aa,feat in ExposureCN(model) }

      for aa in set(featuresDict.keys()).union(set(hseDict.keys())):
        try:
          prevFeatures= featuresDict[aa]
        except KeyError:
          prevFeatures= [-1, -1, -1.0, -1, -1, -1.0]
        try:
          newFeatures= hseDict[aa]
        except KeyError:
          newFeatures= [-1 ]      
        featuresDict[aa]= prevFeatures + newFeatures
        
      dataDict={}
      for aa in allResidues:
        chainId_resIdStr_resName= self.fromRes2ChainResIdAndName(aa)
        chainId, resIdStr, resName= chainId_resIdStr_resName
        if resName is None: continue
        if aa in featuresDict:
          values= [ str(val) for val in featuresDict[aa] ]
        else:
          values= ["-1", "-1", "-1.0", "-1", "-1", "-1.0","-1"]

        record= [ chainId, resIdStr, resName] + values
        record= " ".join(record)
        try:
          dataDict[chainId].append(record)
        except KeyError:
          dataDict[chainId]=[record]
          
      self.writeResultsFromDataDictSingleChain(dataDict, outName= self.getFNames(prefixExtended)[0])
        
    except (Exception, KeyboardInterrupt):
      self.tryToRemoveAllFnames(prefixExtended)
      raise


