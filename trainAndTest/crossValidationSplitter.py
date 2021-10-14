# -*- coding: utf-8 -*-
from sklearn.model_selection import KFold, GroupKFold, StratifiedKFold
from utils import getItemsFromList

N_SUB_FOLDS= 3

def subSplitFolds(trainPrefixes_testPrefixes_list, prefixes, random_state=None):
  trainPrefixes_testPrefixes_list= list( trainPrefixes_testPrefixes_list )
  groups= {}
  for i,(trainPrefixes, testPrefixes) in enumerate(trainPrefixes_testPrefixes_list):
    for prefix in testPrefixes:
      groups[prefix]= i
  testIdxsToTrainIdxs= []
#  print(trainPrefixes_testPrefixes_list);print(groups); raw_input("enter")
  for trainPrefixes, testPrefixes in trainPrefixes_testPrefixes_list:
    train_list_tmp=[]
    current_groups=  [ groups[trainPrefix] for trainPrefix in trainPrefixes]
    try:
      splits= list( GroupKFold(N_SUB_FOLDS).split(trainPrefixes, groups= current_groups))
    except ValueError as e:
      msg = (" Error, the PARAM N_SUB_FOLDS is larger than the different groups  contained in the train split "+
            "%s. This training fold should be enriched with additional training classes or the number of folds increased. Groups %s"%(str(trainPrefixes),
                                                                                                     str(current_groups)))
      raise ValueError(e.message + msg)
    for (idxList_i, __) in splits:
      trainPrefixes_subset= getItemsFromList(idxList_i, trainPrefixes)
      trainIdxs_subset= [ prefixes.index(prefix) for prefix in  trainPrefixes_subset]
      train_list_tmp.append( trainIdxs_subset )

    testIdxs= [ prefixes.index(prefix) for prefix in  testPrefixes]
    testIdxsToTrainIdxs.append( (testIdxs, train_list_tmp) )
  testIdxsToTrainIdxs= sorted( testIdxsToTrainIdxs, key= lambda x: x[0][0] )
  return testIdxsToTrainIdxs


def mergeSplitFolds(dataIds_multistep, trainPrefixes_testPrefixes_list, prefixesUsedInModel ):
  '''

  :param dataIds_multistep:
  :param trainPrefixes_testPrefixes_list:
  :param prefixesUsedInModel: Dict {currentPrefix:[prefixes used for training the model that predicted currentPrefix] }
  :return: testIdxsToTrainIdxs
  '''
  fromOriPrefix_to_prefixMult= {}
  for i, prefix_mult in enumerate( dataIds_multistep):
    ori_prefix= prefix_mult.split("@")[0]
    if not ori_prefix in fromOriPrefix_to_prefixMult:
      fromOriPrefix_to_prefixMult[ori_prefix]= set([])
    fromOriPrefix_to_prefixMult[ori_prefix].add( (prefix_mult, i) )
  testIdxsToTrainIdxs= []
  for foldNum, (trainPrefixes_ori, testPrefixes_ori) in enumerate(trainPrefixes_testPrefixes_list):
    # testPrefixes_ori= [elem for elem in testPrefixes_ori if elem[:4].isupper() ]
    testPrefixes_idx= [ i for i, prefix in enumerate(dataIds_multistep) if prefix.split("@")[0] in testPrefixes_ori]    
#    print(trainPrefixes_ori, testPrefixes_ori)
    trainPrefixes_idx=[]
    for ori_trainPrefix in trainPrefixes_ori:
      for new_trainPrefix, idx in fromOriPrefix_to_prefixMult[ori_trainPrefix]:
        if new_trainPrefix[-3:] in ["#sl", "#sr"]:
          new_trainPrefix= new_trainPrefix[:-3] 
        involvedPrefixes= set(prefixesUsedInModel[new_trainPrefix])
        # print("%s ------------> %s"%(new_trainPrefix, filter(lambda x: x[:4].isupper(),involvedPrefixes)))
        if not involvedPrefixes.intersection(testPrefixes_ori):
          trainPrefixes_idx.append( idx)
    if len(trainPrefixes_idx)==0: raise ValueError("Error, the train/test split for fold %d is not orthogonal for %s/%s"%(
      foldNum, trainPrefixes_ori, testPrefixes_ori) +". You may try to increase number of subfolds or k-fold number")
    testIdxsToTrainIdxs.append( (testPrefixes_idx, [trainPrefixes_idx]) )
#  print(trainPrefixes_testPrefixes_list)
#  print(testIdxsToTrainIdxs)
  return testIdxsToTrainIdxs


if __name__=="__main__":
  import joblib
  dataIds_multistep= joblib.load("/home/ruben/tmp/BIPSPI/crossValSplitter/dataIds_multistep.joblib")
  prefixesUsedInModel = joblib.load("/home/ruben/tmp/BIPSPI/crossValSplitter/prefixesUsedInModel.joblib")
  trainPrefixes_testPrefixes_list = joblib.load("/home/ruben/tmp/BIPSPI/crossValSplitter/trainPrefixes_testPrefixes_list.joblib")
  mergeSplitFolds(dataIds_multistep, trainPrefixes_testPrefixes_list, prefixesUsedInModel )