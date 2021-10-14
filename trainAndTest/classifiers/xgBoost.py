import xgboost as xgb
import numpy as np

DEBUG=False
USE_GPU= False
def trainMethod(trainData, trainLabels, verboseLevel=0, ncpu= 1, useGpu=USE_GPU):
  print("Training shape %s"%(str(trainData.shape)))

# gpu params 'tree_method':"gpu_hist" max_bin':32, 'n_gpus':1

#  params= {'objective':'binary:logistic', 'colsample_bytree': 0.9, 'learning_rate': 0.1, 'min_child_weight': 1,
#           'n_estimators': 2000, 'subsample': 0.9, 'reg_lambda': 10.0, 'max_depth': 12, 'gamma': 0, 'nthread': ncpu}

  params= {'objective':'binary:logistic', 'colsample_bytree': 0.9, 'learning_rate': 0.1, 'min_child_weight': 1,
           'n_estimators': 2000, 'subsample': 0.9, 'reg_lambda': 10.0, 'max_depth': 12, 'gamma': 0, 'nthread': ncpu,
           'tree_method':"hist", 'max_bin':32}  #scale_pos_weight= 3./1.

  if useGpu:
    params['tree_method']= "gpu_hist"
    params['max_bin']= 32
    params['n_gpus']= 1
    params['predictor']= "cpu_predictor"
  if DEBUG:
    params[ 'n_estimators']= 20
    print("WARNING: DEBUG MODE ")
  print(params)
  modelo= xgb.XGBClassifier( **params)
  modelo.fit(trainData, trainLabels)
  return modelo

def predictMethod(modelo, testData):
  print("Predict shape %s"%(str(testData.shape)))  
  return modelo.predict_proba(testData)[:,1]
