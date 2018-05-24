
from sklearn.model_selection import GridSearchCV, cross_val_score, learning_curve
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
import sys, os
from time import time
import numpy as np

def getModel0(verboseLevel):
  from sklearn.ensemble import RandomForestClassifier
  model= RandomForestClassifier( n_estimators= 800, verbose=verboseLevel, n_jobs=-1)
  # {'max_features': 70, 'n_estimators': 1000, 'max_depth': None} struct train= 1.000000  test= 0.879251
  # {'max_features': 80, 'n_estimators': 800, 'max_depth': None}  struct train= 1.000000  test= 0.880443
  # {'n_estimators'= 1000, 'max_features': 100, 'max_depth':None} seq train= 1.000000  test= 0.848445
  param_grid = {"n_estimators":[ 800], "max_features":[ 60, 70, 80, 90], "max_depth": [None]}
  return model, param_grid

  
def getModel1(verboseLevel):
  import xgboost as xgb
  model= xgb.XGBClassifier( objective="binary:logistic", nthread=-1, silent=verboseLevel)
  #seq 0.8787 {'n_estimators': 2000, 'subsample': 0.9, 'learning_rate': 0.1, 'max_depth': 12,'min_child_weight': 1}
  #seq 0.8802 {'colsample_bytree': 0.9, 'learning_rate': 0.1, 'min_child_weight': 1, 'n_estimators': 2000, 'subsample': 0.9, 'reg_lambda': 10.0, 'max_depth': 12, 'gamma': 0}
  #seq 0.8804 {'colsample_bytree': 0.9, 'learning_rate': 0.1, 'min_child_weight': 1, 'n_estimators': 2000, 'subsample': 0.9, 'reg_lambda': 100.0, 'max_depth': 12, 'gamma': 0}

  #best scores struct: train= 1.000000  test= 0.925382  {'colsample_bytree': 0.9, 'learning_rate': 0.1, 'n_estimators': 1000, 'subsample': 0.9, 'max_depth': 16, 'gamma': 0})
  
  #mixed 0.8979 Parameters: {'colsample_bytree': 0.9, 'learning_rate': 0.1, 'min_child_weight': 1, 'n_estimators': 2000, 'subsample': 0.9, 'max_depth': 12, 'gamma': 0}

  #mixed 0.8992 Parameters: {'colsample_bytree': 0.9, 'learning_rate': 0.1, 'min_child_weight': 1, 'n_estimators': 2000, 'subsample': 0.9, 'reg_lambda': 10.0, 'max_depth': 12, 'gamma': 0}


  param_grid = {"n_estimators":[2000], "learning_rate":[0.1], "max_depth":[12, 13],  'min_child_weight':[1],
                 "subsample":[ 0.9], "gamma":[0], 'colsample_bytree':[0.9], 'reg_lambda':[ 1e1]}
  return model, param_grid

getModel= getModel1


getModel= getModel0
def tuneClassifier(trainData, trainLabels, trainGroups, nFolds= 10, ncpu=1, verboseLevel=0):
  print("Training shape", trainData.shape)
  model, param_grid= getModel(verboseLevel)
  
  # run grid search
  print("Running grid search with n_jobs: %d. Verbose level: %d"%(ncpu,verboseLevel))
  if nFolds==-1:
    crossValIter= LeaveOneGroupOut()
    print("Leave one out will be performed")
  elif nFolds>=2:
    crossValIter= GroupKFold(n_splits= nFolds)
    print("k=%d fold cross validation will be performed"%(nFolds))    
  else:  
    raise ValueError("nFolds must be >=2 or -1 (for Leave-one-out)")
    
  grid_search = GridSearchCV(model, cv=crossValIter,  param_grid=param_grid, scoring= "roc_auc", n_jobs=ncpu, verbose= verboseLevel)
  start = time()
  grid_search.fit(trainData, trainLabels, trainGroups)
  print("GridSearchCV took %.2f seconds for %d candidate parameter settings."
        % (time() - start, len(grid_search.cv_results_['params'])))
  bestParams= reportBestParams(grid_search.cv_results_)
  model, __= getModel(verboseLevel)
  model.set_params(**bestParams)
  plot_learning_curve(model, "learning_curve", trainData, trainLabels)
  print("Best params:", bestParams)
  
  sys.exit(0)

def reportBestParams(results, n_top=5):
  ''' 
    Function to report best scores 
  '''
  bestCandidate=None
  for i in range(1, n_top + 1):
    candidates = np.flatnonzero(results['rank_test_score'] == i)
    for candidate in candidates:
      print("Model with rank: {0}".format(i))
      print("Mean validation score: {0:.4f} (std: {1:.4f})".format(
            results['mean_test_score'][candidate],
            results['std_test_score'][candidate]))
      print("Parameters: {0}".format(results['params'][candidate]))
      if bestCandidate==None:
        bestCandidate= results['params'][candidate]
      print("")
  return bestCandidate
      
def plot_learning_curve(estimator, title, X, y, ylim=None, cv=None,
                        n_jobs=1, train_sizes=np.linspace(.1, 1.0, 10), doPlot=True):
    import matplotlib.pyplot as plt
    """
    Generate a simple plot of the test and training learning curve.

    Parameters
    ----------
    estimator : object type that implements the "fit" and "predict" methods
        An object of that type which is cloned for each validation.

    title : string
        Title for the chart.

    X : array-like, shape (n_samples, n_features)
        Training vector, where n_samples is the number of samples and
        n_features is the number of features.

    y : array-like, shape (n_samples) or (n_samples, n_features), optional
        Target relative to X for classification or regression;
        None for unsupervised learning.

    ylim : tuple, shape (ymin, ymax), optional
        Defines minimum and maximum yvalues plotted.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
          - None, to use the default 3-fold cross-validation,
          - integer, to specify the number of folds.
          - An object to be used as a cross-validation generator.
          - An iterable yielding train/test splits.

        For integer/None inputs, if ``y`` is binary or multiclass,
        :class:`StratifiedKFold` used. If the estimator is not a classifier
        or if ``y`` is neither binary nor multiclass, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validators that can be used here.

    n_jobs : integer, optional
        Number of jobs to run in parallel (default 1).
    """
    plt.figure()
    plt.title(title)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=cv, n_jobs=n_jobs, train_sizes=train_sizes)
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)
    plt.grid()

    plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1,
                     color="r")
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r",
             label="Training score")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g",
             label="Cross-validation score")

    plt.legend(loc="best")
    print("estimator parms", estimator.get_params())
    print("best scores: train= %f  test= %f"%(train_scores_mean[-1], test_scores_mean[-1]))
    if doPlot:
      plt.show()
    return plt
    
