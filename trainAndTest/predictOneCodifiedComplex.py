from __future__ import print_function
import sys, os
from joblib import load as joblib_load
from .resultsManager import ResultsManager
from .classifiers.xgBoost import  predictMethod

from Config import Configuration

class ComplexPredictor(Configuration):
  '''
    This class is used to predict new pdbs once models have already been computed.
  '''
  def __init__(self, stepName, isHomoComplex, savedModelsPath=None, averageLRscores=False):
    '''

      :param stepName: str. Must startswith seq_train or struct or mixed (seq_train, mixed_2, structX, seq_train1... are also valid)
      :param isHomoComplex: boolean. Is the target complex homo or hetero
      :param savedModelsPath: str. A path to the directory where models have been saved. If None,
                                   it will used the path indicated in Config
      :param averageLRscores: True if Ligand and receptor are the same protein and thus, binding site prediction should be averaged
    '''
    Configuration.__init__(self)

    self.isHomoComplex = isHomoComplex
    self.stepName= stepName
    self.averageLRscores= averageLRscores
    if not savedModelsPath is None:
      self.savedModelsPath= savedModelsPath

    self.model=None
    print(stepName)
    self.savedModelsPath= os.path.join(self.savedModelsPath, "homo" if self.isHomoComplex else "hetero")
    for fname in os.listdir(self.savedModelsPath):
      if fname.endswith(stepName):
        print("Loading model %s %s"%("homo" if isHomoComplex else "hetero",fname))
        self.model= joblib_load(os.path.join(self.savedModelsPath, fname))
    assert not self.model is None, "Error, there is no valid model in %s for step %s"%(self.savedModelsPath, self.stepName)

  def predictOneComplex(self, complexCodifiedObject, outName, isMixedProto, isLastStep=False):
    '''
      computes pairwise RR interaction and binding-site predictions for the codified complex passed as argument
 
      :param complexCodifiedObject: A codifyComplexes.ComplexCodified object that represents a protein protein
                                    interaction
      :param outName: str.  The name that RR interaction results file will have. Ligand binding-site results
                            name will be outName+'lig' and Receptor binding-site results
                            name will be outName+'rec'
      :param isMixedProto: True if mixed protocol
      :param outName: boolean. True if no more feedback steps will be done
    '''
    data_d, data_t= complexCodifiedObject.getData()
    data_ids= complexCodifiedObject.getIds()
    prefix= complexCodifiedObject.getPrefix()
    data_d_preds= None if isinstance(data_d, type(None)) else predictMethod(self.model, data_d)
    data_t_preds= None if isinstance(data_t, type(None)) else predictMethod(self.model, data_t)

    if isMixedProto:
      prob_predictionsDir= data_d_preds if not isinstance(data_d_preds, type(None)) else data_t_preds
      resultObj= ResultsManager(prefix, prob_predictionsDir, prob_predictionsTrans=None, ids=data_ids,
                                averageLRscores= self.averageLRscores, doAverageScore=False)
    else:
      resultObj= ResultsManager(prefix, data_d_preds, data_t_preds, data_ids, averageLRscores= self.averageLRscores,
                                doAverageScore=True)
    p,l,r= resultObj.getResultsPdDf("p"), resultObj.getResultsPdDf("l"), resultObj.getResultsPdDf("r")
    resultObj.writeResults(outName, isLastStep)
    return  outName, (p,l,r)
