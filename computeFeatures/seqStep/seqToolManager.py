from __future__ import absolute_import
from ..toolManagerGeneric import ToolManager
from itertools import chain as IterChainer

class SeqToolManager(ToolManager):
  '''
  Abastract class that will be used to implement psiblast, al2co ... 
  '''
  def __init__(self, computedFeatsRootDir, winSize=None, statusManager=None):
    '''
      :param computedFeatsRootDir: str. root path where results will be saved
      :param winSize: int>=1 or None. The size of the windows for sliding window if desired or Nones
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''   

    ToolManager.__init__(self, computedFeatsRootDir, statusManager)
    self.winSize= winSize
    assert winSize is None or  winSize%2==1, "Error, winSize must be odd"
    
    
  def getFNames(self, prefixExtended):
    '''
    abstract method
    Returns a list that contains the fnames that will be output by psiblast and other tools
    :param prefixExtended. prefix for output fnames.
    :return a list   [fname1,fnam2...]
    '''
    raise ValueError("Not implemented")
    
  def computeFromSeqStructMapper(self, seqsManager, prefixExtended):
    '''
    :param seqsManager: computeFeatures.seqStep.seqToolManagers.seqExtraction.SeqStructMapper
    :param prefixExtended. prefix for output fnames. They are form as follows: prefix+chainType+chainId
    '''
    raise ValueError("Not implemented")
    
  def makeWindowed(self, dataList, featsNames, featsDefault, featsLevels, outName):
    '''
    :param dataList: [ (res1_Feats), (res2_Feats), ...]
            resI_Feats: ( (chainId, resId, resName), ([feats1], [feast2]) )  e.g. feats1:[pssm0, pssm1], feats2:[entropy1, entropy2])
    :param featsNames: [str] Name for feats1 feats2...  e.g. ["pssms", "entropy"]
    :param featsDefault: integer or character e.g. [-1, 0.3]
    :param featsLevels: [ [str] or None ]. For instance, the levels for pssm and aaCode will be  [ None, ["A","V","L"...]]
    :param outName: The path wehere the windowed version will be saved
    '''
    featsBase_Names= featsNames
    featsNames= [ elem+"Win" for elem in featsNames]
    numElements= len(dataList)
    nFeats= len(featsNames)
    exampleResId, exampleFeatsList= dataList[0]
    chainId= exampleResId[0]
    winSizeStep= int(self.winSize/2)
    noAA_FeatsProfile= [ [featsDefault[i]]*len(exampleFeatsList[i]) for i in range(nFeats)] 
    resIds, feats= zip(* dataList )
    list_of_profilesList= zip(* feats)
    resultsList=[]
    for resNum, full_resId  in enumerate(resIds):
      oneResProfile= [list(full_resId)]
      for featNum in range(nFeats):
        oneResProfile+= [  list_of_profilesList[featNum][i] if (i>=0 and i<= numElements-1) 
                                                                  else noAA_FeatsProfile[featNum]
                                                    for i in range(resNum-winSizeStep, resNum+winSizeStep+1)]

      resultsList.append( " ".join(list(IterChainer.from_iterable(oneResProfile))))

    level_idxs= []
    currI= len(exampleResId)
    assert currI == 3, "Error, resIds are tuples of 3 elements"
    featsNames_=[]
    for featNum,featName in enumerate(featsNames):
      for varNum in range(len(exampleFeatsList[featNum])):
        for i in range(self.winSize):
          featsNames_.append( featName+".%d.%d"%(varNum, i) )
    featsNames= featsNames_
    categoricalLevels={}
    for featNamePrefix, featLevels in zip(featsBase_Names, featsLevels):
      if not featLevels is None:
        categoricalLevels[ tuple(featLevels) ]= [ featName for featName in featsNames if featName.startswith(featNamePrefix)]
        
    if len(categoricalLevels)==0: categoricalLevels=None
    dataDict= {chainId: resultsList }
    self.writeResultsFromDataDictSingleChain(dataDict, featuresNames=featsNames, outName= outName, 
                                             categoricalLevels=categoricalLevels)
      

  
def testMakeWindow():
  computedFeatsRootDir= "/home/rsanchez/tmp"
  winSize=3
  stm= SeqToolManager(computedFeatsRootDir, winSize)
  dataList=[ [("chain%s"%i, "res%d"%i, "resname%d"%i),( [str(i), "b%s"%i],[i%2],[i+1,i-1])] for i in range(7)]
  featsNames= ["featLetras", "featOneNum", "feat2Nums"]
  featsDefault= [ "caca", -1, -1024]
  featsLevels=[["a","b"], ["0","1","2"], None]
  outName= None
  stm.makeWindowed(dataList, featsNames, featsDefault, featsLevels, outName)
  
if __name__=="__main__":
  testMakeWindow()
  print("Done")
  
