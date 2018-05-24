import xgboost as xgb
import numpy as np

def trainMethod(trainData, trainLabels, verboseLevel=0, ncpu= 1):
  print("Training shape %s"%(str(trainData.shape)))
#  if trainData.shape[1]>1300:
#    params= {'objective':'binary:logistic', 'n_estimators': 2000, 'subsample': 0.9, 'learning_rate': 0.1, 'max_depth': 12,
#             'min_child_weight': 1, 'nthread': ncpu, 'silent': verboseLevel}
#  else:
#    params= {'objective':'binary:logistic','colsample_bytree': 0.9, 'learning_rate': 0.1, 'n_estimators': 1000,
#             'subsample': 0.9, 'max_depth': 16, 'gamma': 0, 'nthread': ncpu, "silent": verboseLevel}

  params= {'objective':'binary:logistic', 'colsample_bytree': 0.9, 'learning_rate': 0.1, 'min_child_weight': 1,
           'n_estimators': 2000, 'subsample': 0.9, 'reg_lambda': 10.0, 'max_depth': 12, 'gamma': 0, 'nthread': ncpu}
  modelo= xgb.XGBClassifier( **params)
  modelo.fit(trainData, trainLabels)
  return modelo
  
def predictMethod(modelo, testData):
  return modelo.predict_proba(testData)[:,1]
