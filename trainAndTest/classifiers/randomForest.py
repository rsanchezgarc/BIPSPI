from sklearn.ensemble import RandomForestClassifier as EnsembleClassifier

import numpy as np

def trainMethod(trainData, trainLabels, verboseLevel=0, ncpu= 1):
  print("Training shape %s"%(str(trainData.shape)))
  if trainData.shape[1]>1000:
    params= {"n_estimators": 1000, 'max_features': 100, "max_depth":None, "n_jobs": ncpu, "verbose": verboseLevel}
  else:
    params= {"n_estimators": 800, 'max_features': 80, "max_depth":None, "n_jobs": ncpu, "verbose": verboseLevel}
  modelo= EnsembleClassifier( **params)
  modelo.fit(trainData, trainLabels)
  return modelo
  
def predictMethod(modelo, testData):
  return modelo.predict_proba(testData)[:,1]
