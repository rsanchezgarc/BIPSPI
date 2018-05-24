from __future__ import print_function
import itertools
import sys, os
import inspect
import numpy as np
from joblib import load as joblib_load
from .resultsManager import ResultsManager
from .classifiers.randomForest import  predictMethod

from Config import Configuration

class ComplexPredictor(Configuration):
  '''
    This class is used to predict new pdbs once models have already been computed.
  '''
  def __init__(self, stepName, savedModelsPath=None):
    '''

      @param stepName: str. Must startswith seq_train or struct or mixed (seq_train, mixed_2, structX, seq_train1... are also valid)
      @param savedModelsPath: str. A path to the directory where models have been saved. If None, 
                                   it will used the path indicated in Config
    '''
    Configuration.__init__(self)

    self.stepName= stepName
    if not savedModelsPath is None:
      self.savedModelsPath= savedModelsPath

    self.model=None
    print(stepName)
    for fname in os.listdir(self.savedModelsPath):
      if fname.endswith(stepName):
        print("Loading model %s"%(fname))
        self.model= joblib_load(os.path.join(self.savedModelsPath, fname))
    assert not self.model is None, "Error, there is no valid model in %s for step %s"%(self.savedModelsPath, self.stepName)

  def predictOneComplex(self, complexCodifiedObject, outName, isLastStep=False):
    '''
      computes pairwise RR interaction and binding-site predictions for the codified complex passed as argument
 
      @param complexCodifiedObject: A codifyComplexes.ComplexCodified object that represents a protein protein
                                    interaction
      @param outName: str.  The name that RR interaction results file will have. Ligand binding-site results
                            name will be outName+'lig' and Receptor binding-site results
                            name will be outName+'rec'
      @param outName: boolean. True if no more feedback steps will be done
    '''
    data_d, data_t= complexCodifiedObject.getData()
    data_ids= complexCodifiedObject.getIds()
    prefix= complexCodifiedObject.getPrefix()
    prob_predictionsDir= predictMethod(self.model, data_d)
    prob_predictionsTrans= predictMethod(self.model, data_t)
    resultObj= ResultsManager(prefix, prob_predictionsDir, prob_predictionsTrans, data_ids)
    p,l,r= resultObj.getResultsPdDf("p"), resultObj.getResultsPdDf("l"), resultObj.getResultsPdDf("r")
    resultObj.writeResults(outName, isLastStep)
    return  outName, (p,l,r)

