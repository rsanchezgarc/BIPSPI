from __future__ import print_function
import os
import numpy as np
import pandas as pd
from collections import Counter
from .evaluateResults import evaluateBindingSite, evaluateBothBindingSites, evaluatePairs, loadAccesibility
pd.set_option('precision', 4)


DO_SEQ_AVER_BIND_SITE=True
EVALUATE_L_R_INDEPENDTLY= False

class ResultsManager(object):
  '''
    This class will store predictions for one complex. It will also allow for results evaluations and
    results writting to disk

  '''
  def __init__(self, prefix, prob_predictionsDir, prob_predictionsTrans, ids, doAverageScore=True):
    '''
      builder
 
       @param prefix: str. An identifier of the pdb file.
       @param prob_predictionsDir: float[]. scores for Residue Residue interaction for pairs in direct mode
                                            (first ligand amino acid, second receptor amino acid)
       @param prob_predictionsTrans: float[]. scores for Residue Residue interaction for pairs in transpose mode
                                            (first receptor amino acid, second ligand amino acid). Same order
                                            than prob_predictionsDir
       @param ids: str[]. ids of the pairs whose scores direct and transpose are prob_predictionsDir and 
                          prob_predictionsTrans. Same order than prob_predictionsDir
       @param doAverageScore: boolean. If True:
                                          Final pair score will be mean of prob_predictionsDir and  prob_predictionsTrans
                                       else:
                                          Final pair score will be prob_predictionsDir                                      
    '''
    
    self.getInterfacesScores= self.getInterfacesScoresMiCount
#    self.getInterfacesScores= self.getInterfacesScoresPairPred
        
    self.prefix= prefix
    self.prob_predictionsDir= prob_predictionsDir
    self.prob_predictionsTrans= prob_predictionsTrans
    self.ids= ids
    self.testLabels= ids["categ"].values
    
    self.performance_summary= None

    if doAverageScore == True:
      meanPairsScore= np.mean(np.concatenate( [[self.prob_predictionsDir],
                                               [self.prob_predictionsTrans]], axis=0),axis=0)
    else:
      meanPairsScore= self.prob_predictionsDir

    self.resultsPairs= self.ids
    self.resultsPairs["prediction"]= meanPairsScore
    self.resultsLigand, self.resultsRecep= self.getInterfacesScores(self.resultsPairs)
    self.l_residues_pos= None
    self.r_residues_pos= None

    
  def getIdsForInteface(self, oneDf):
    '''
      Given a pandas.DataFrame that represents the pairwise predicitions, obtains
      the resIds that belongs to the ligand and to the receptor
       
      @param oneDf: pandas.DataFrame. A data frame that has to have at least the following columns:
                                        chainIdL structResIdL chainIdR structResIdR
      @return l_res, r_res
              l_res: {resId: Number of times resId appears in a row}. For ligand amino acids
              r_res: {resId: Number of times resId appears in a row}. For receptor amino acids
                Residues ids are form by concatenation of chainId and resId using as separator _
                 ("A", "123") --> "A_123"
    '''
            
    l_res= oneDf["chainIdL"].map(str) + "_" +oneDf["structResIdL"].map(str)
    l_res= dict(Counter(l_res))
    r_res= oneDf["chainIdR"].map(str) + "_" +oneDf["structResIdR"].map(str)
    r_res= dict(Counter(r_res))
    return l_res, r_res


  def getResultsPdDf(self, chainType="p", removeCateg=True):
    '''
      returns a pd.DataFrame that contains results of pair predictions when (chainType=="p") or
      single chain predictions (when chainType=="r" or chainType=="l")
      @param chainType: str. "p" for residue-residue predictions, "r" for receptor predictions and "l" for
                             ligand predictions
      @param removeCateg: boolean.  If removeCateg==False files for DataFrames will contain categ column,
    '''
    if chainType=="p":
      df= self.resultsPairs
    elif chainType in ["r","l"]:
      df= self.resultsLigand if chainType=="l" else self.resultsRecep
      df= pd.DataFrame.from_dict(df, orient="index")
      chainIds_ResIds= df.index.map(lambda x: x.split("_"))
      chainIds, resIds= zip(* chainIds_ResIds)
      df= pd.DataFrame( {"chainId": chainIds, "resId":resIds, "prediction": df.values[:,-1]})
      df= df.reindex( ["chainId", "resId", "prediction"], axis=1)
    else:
      raise ValueError("chainType must be 'p', 'r' or 'l'")
    if removeCateg and "categ" in df:
      df= df.drop("categ", axis=1)
    df= df.sort_values("prediction", ascending=False) 
    return df
    
  def getSeqAverScore(self, scoresDict):
    '''
    Averages interface scores over sequence window of size 3
    @param scoresDict: {str:float} {chainAndResId: binding-siteScore}
    @return newScoresDict: {str:float} {chainAndResId: NewbindingSiteScore}
    '''
    if not DO_SEQ_AVER_BIND_SITE:
      return scoresDict
    scores_list=[]
    ids_list=[]
    chainsDict={}
    for elem in scoresDict:
      chain, resId = elem.split("_")
      if chain not in chainsDict:
        chainsDict[chain]=[]
      score= scoresDict[elem]
      chainsDict[chain].append( ((int(resId), ""), score) if resId[-1].isdigit() else ((int(resId[:-1]), resId[-1]), score ) )
    for chain in chainsDict:
      result= chainsDict[chain]
      result.sort(key= lambda x: x[0])
      ids, scores= zip(* result)
      scores= list(np.convolve(scores, np.array([1, 3, 1])/5.0, mode='same')+ np.array(scores))
      scores_list+= scores
      ids_list+= [ chain+"_"+str(elem[0])+elem[1] for elem in ids]
    scores_dict= dict( zip(ids_list, scores_list))
    return scores_dict
  
  def getInterfacesScoresMiCount(self, resDf):
    '''
    Computes binding-site scores from pairwise results contained in resDf
    To do so, it counts the number of times one ligand (receptor) residue
    appears in the 2**i highest score pair predicitions and divides this
    count by 2**i. i ranges from 3 to 11. The sum of all those values will
    be the final score.
    
    @param resDf: pandas.DataFrame. A DataFrame that represents the pairs predictions. It has the following columns:
                          chainIdL structResIdL resNameL chainIdR structResIdR resNameR categ prediction
    @return scores_l, scores_r
            scores_l: {resId: binding-site score}. Binding site score for each residue of the ligand
            scores_l: {resId: binding-site score}. Binding site score for each residue of the ligandeceptor
              Residues ids are form by concatenation of chainId and resId with _
               ("A", "123") --> "A_123"
    '''
    predDf= resDf.sort_values(by="prediction", ascending=False).reset_index(drop=True)
    allIds_l, allIds_r= self.getIdsForInteface(predDf)
    scores_l= { elem:0.0 for elem in allIds_l}
    scores_r= { elem:0.0 for elem in allIds_r}
    for i in reversed(range(3,12)):
      eval_at=2**i
      predDf= predDf[:eval_at]  # If eval_at>=predDf.shape[0] it behaves as eval_at= predDf.shape[0], so no border concerns
      pred_l_res, pred_r_res= self.getIdsForInteface(predDf) # pred_l_res and pred_r_res are {id_: number of times it appear}
      eval_at= float(eval_at)
      for elem in pred_l_res:
        scores_l[elem]+= pred_l_res[elem]/eval_at
      for elem in pred_r_res:
        scores_r[elem]+= pred_r_res[elem]/eval_at
    scores_l= self.getSeqAverScore(scores_l)
    scores_r= self.getSeqAverScore(scores_r)
    
#    scores_l= self.correctScoresWithNumResAcces(scores_l, chainType="l")
#    scores_r= self.correctScoresWithNumResAcces(scores_r, chainType="r")
    
    return scores_l, scores_r

  def correctScoresWithNumResAcces(self, scoresDict, chainType):
  
    accessSet,nonAcessSet= loadAccesibility(self.prefix, chainType)
    accesByChain= {}
    for resId in accessSet:
      chainId=resId[0]
      if not chainId in accesByChain:
        accesByChain[chainId]=1.0
      else:
        accesByChain[chainId]+=1.0
    for resId in scoresDict:
#      scoresDict[resId]= scoresDict[resId]/ accesByChain[resId[0]]
      scoresDict[resId]= scoresDict[resId]/ float(len(accessSet))
    return scoresDict
    
  def getInterfacesScoresPairPred(self, resDf):
    '''
    Computes binding-site scores from pairwise results contained in resDf
    To do so it employed pairPred strategy, f(a) = max(f_pair(a,bj))
    
    @param resDf: pandas.DataFrame. A DataFrame that represents the pairs predictions. It has the following columns:
                          chainIdL structResIdL resNameL chainIdR structResIdR resNameR categ prediction
    @return scores_l, scores_r
            scores_l: {resId: binding-site score}. Binding site score for each residue of the ligand
            scores_l: {resId: binding-site score}. Binding site score for each residue of the ligandeceptor
              Residues ids are form by concatenation of chainId and resId with _
               ("A", "123") --> "A_123"
    '''
    allIds_l, allIds_r= self.getIdsForInteface(resDf)
    scores_l= { elem:0.0 for elem in allIds_l}
    scores_r= { elem:0.0 for elem in allIds_r}

    for chainId_resId in scores_l:
      chainId, resId= chainId_resId.split("_")
      scores_l[chainId_resId]=np.max( resDf.loc[ (resDf["chainIdL"]==chainId) & (resDf["structResIdL"]==resId) , "prediction"])
    for chainId_resId in scores_r:
      chainId, resId= chainId_resId.split("_")
      scores_r[chainId_resId]=np.max( resDf.loc[ (resDf["chainIdR"]==chainId) & (resDf["structResIdR"]==resId) , "prediction"])
    return scores_l, scores_r
    
  def bindingSiteScoresAndLabels(self):

    scores= []
    categ=[]
    for chainType in ["l", "r"]:
      residuesAtInterface= self.l_residues_pos if chainType=="l" else self.r_residues_pos
      residuesAtInterface= residuesAtInterface if not residuesAtInterface is None else set({})
      resultsOneProt= self.resultsLigand if chainType=="l" else self.resultsRecep
      for chainId_resId in sorted(resultsOneProt):
        chain, resId= chainId_resId.split("_")
        scores.append(resultsOneProt[ chainId_resId ])
        if chainId_resId in residuesAtInterface:
          categ.append(1)
        else:
          categ.append(-1 )
    return scores, categ
 
  def writeInterfResults(self, outName, resultsOneProt, chainType="l", removeCateg=False):
    '''
    Writes binding sites results to disk (ligand or receptor)
    
    @param outName: str.  Ligand binding-site results
                          name will be outName+'lig' and Receptor binding-site results
                          name will be outName+'rec'
    @param resultsOneProt: {str: float} {resId: binding_siteScore}. The binding-site score for all
                           residues of ligand or receptor.
                           
    @param chainType: str. "l" for ligand and "r" for receptor
    
    @param removeCateg: boolean.  If removeCateg==False files for DataFrames will contain categ column,

    '''

    fname= outName+".lig" if chainType=="l" else outName+".rec"
    residuesAtInterface= self.l_residues_pos if chainType=="l" else self.r_residues_pos
    
    if residuesAtInterface is None:
      posDf= self.resultsPairs[self.resultsPairs["categ"]==1] # all pairs whose residues are at interface
      if posDf.shape[0]>0:
        l_res_pos, r_res_pos= self.getIdsForInteface(posDf)  #{'A_39': 3, 'A_44': 2}
        self.l_residues_pos= l_res_pos
        self.r_residues_pos= r_res_pos
        residuesAtInterface= self.l_residues_pos if chainType=="l" else self.r_residues_pos
    residuesAtInterface= residuesAtInterface if not residuesAtInterface is None else set({})
    trueCategLabel, falseCategLabel= 1,-1
    if self.resultsPairs["categ"].isnull().values.all():
      trueCategLabel, falseCategLabel= np.nan, np.nan
    with open(fname,"w") as f:
      chainIds= []
      resIds= []
      scores= []
      categ=[]
      for chainId_resId in sorted(resultsOneProt):
        chain, resId= chainId_resId.split("_")
        chainIds.append(chain)
        resIds.append(resId)
        scores.append(resultsOneProt[ chainId_resId ])
        if chainId_resId in residuesAtInterface:
          categ.append(trueCategLabel)
        else:
          categ.append(falseCategLabel )
          
      results= pd.DataFrame()
      results["chainId"]= chainIds
      results["resId"]= resIds
      if not removeCateg:
        results["categ"]= categ
      if not self.performance_summary is None and EVALUATE_L_R_INDEPENDTLY:
        f.write("#auc= %2.5f\n"%(self.performance_summary["auc_%s"%chainType]))
 
      results["prediction"]= scores
      results.to_csv(f, sep=" ", index=False, na_rep="NaN")
      
  def writeResults(self, outName, removeCateg=False):
    '''
    Writes both pairwise and binding-sites results to disk
    
    @param outName: str.  The name that pairwise results file will have. Ligand binding-site results
                          name will be outName+'lig' and Receptor binding-site results
                          name will be outName+'rec'
    @param removeCateg: boolean.  If removeCateg==False files for DataFrames will contain categ column,

    '''
    results= self.resultsPairs.copy()
    with open(outName,"w") as f:
      if not self.performance_summary is None:
        f.write("#auc= %2.5f\n"%(self.performance_summary["auc_pair"]))
      if removeCateg:
        results.drop("categ", axis=1).to_csv(f, sep=" ", index=False, na_rep="NaN")
      else:
        results.to_csv(f, sep=" ", index=False, na_rep="NaN")
        
    self.writeInterfResults(outName, self.resultsLigand, chainType="l", removeCateg=removeCateg)
    self.writeInterfResults(outName, self.resultsRecep,  chainType="r", removeCateg=removeCateg)

  def getFullEvaluation(self):
    '''
    Performs full evaluation of predictions (pairwise and binding-site level)
    Results will be summarized in a pandas.DataFrame, which is stored at the
    attribute self.performance_summary and also returns it
    
    @param return self.performance_summary: pandas.DataFrame. summary of scores evaluation which have the 
                                            following columns:
                                      pdb  auc_pair  prec_50  reca_50  prec_100  reca_100  prec_500  reca_500   auc_l  prec_l
                                      reca_l   mcc_l   auc_r  prec_r  reca_r   mcc_r
    '''
    resDf=self.resultsPairs
    self.performance_summary= evaluatePairs(self.prefix, resDf)

    posDf= resDf[resDf["categ"]==1] # all pairs whose residues are at interface
    l_res_pos, r_res_pos= self.getIdsForInteface(posDf)  #{'A_39': 3, 'A_44': 2}
    l_res_all, r_res_all= self.getIdsForInteface(resDf)  #{'A_12': 3, 'A_27': 2}
    self.l_residues_pos= l_res_pos
    self.r_residues_pos= r_res_pos
    
    if EVALUATE_L_R_INDEPENDTLY:
      bindingSiteSummaryL= evaluateBindingSite(self.prefix, l_res_pos, l_res_all, self.resultsLigand, chainType="l")
      bindingSiteSummaryR= evaluateBindingSite(self.prefix, r_res_pos, r_res_all, self.resultsRecep,  chainType="r")
      self.performance_summary= pd.concat([self.performance_summary, bindingSiteSummaryL, bindingSiteSummaryR],axis=1 )
    else:
      bindingSiteSummary= evaluateBothBindingSites(self.prefix, l_res_pos, l_res_all, self.resultsLigand, 
                                                   r_res_pos, r_res_all, self.resultsRecep)
      self.performance_summary= pd.concat([self.performance_summary, bindingSiteSummary],axis=1 )
    return self.performance_summary
    
  def getPerformanceSummary(self):
    '''
    returns the attribute self.performance_summary. If self.performance_summary is None, computes it first
    
    @param return self.performance_summary: pandas.DataFrame. summary of scores evaluation which have the 
                                            following columns:
                                      pdb  auc_pair  prec_50  reca_50  prec_100  reca_100  prec_500  reca_500   auc_l  prec_l
                                      reca_l   mcc_l   auc_r  prec_r  reca_r   mcc_r
    '''
    if self.performance_summary is None:
      self.getFullEvaluation()
    return self.performance_summary
    
  def __str__(self):
    if self.performance_summary is None:
      return self.resultsPairs.head().to_string(index=False)
    else:
      return self.performance_summary.to_string(index=False)

    
  @staticmethod
  def loadExistingResults(fname):
    '''
    loads an existing pairwise result file (tabular data) with name fname
    @param fname: str. The name of the already computed result
    @return prevResults: ResultsManager(). An instance of ResultsManager with the previous
                         results loaded
    '''
    prefix= (os.path.split(fname)[-1]).split(".")[0]
    result= pd.read_table(fname, sep='\s+', comment="#", dtype={"structResIdL":str, "structResIdR":str,
                                                                "chainIdL":str, "chainIdR":str})
    ids= result[["chainIdL", "structResIdL", "resNameL", "chainIdR", "structResIdR", "resNameR", "categ"]]
    prob_predictions= result["prediction"].values
    return ResultsManager(prefix, prob_predictions, prob_predictions, ids )

