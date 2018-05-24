suppressMessages(library(pROC))

getRankFirstPos <- function(scoreDf){
  x<-scoreDf[order(scoreDf$prediction,decreasing = T),]
  return(which(x$categ==1)[1])
}

getNumHits <- function(scoresDf,numPairs=500){
  scoresDf<- scoresDf[order(scoresDf$prediction,decreasing = T),]
  scoresDf<- scoresDf[1:numPairs,"categ"]
  scoresDf[scoresDf==-1]<-0
  return(sum(scoresDf))
}

getPrecisionTopPairs <- function(scoresDf,numPairs=500){
  scoresDf<- scoresDf[order(scoresDf$prediction,decreasing = T),]
  categOfTopPairs<- scoresDf[1:numPairs,"categ"]
  categOfTopPairs[categOfTopPairs==-1]<-0 #In case -1 is used as tag
  return(sum(categOfTopPairs)/numPairs)
}

getAUC_ROC <- function(scoresDf){
  return(roc(scoresDf$categ,scoresDf$prediction,direction = "<")$auc)
}


getFullEvaluation <- function(scoresDf,numPairs=500){
  return(data.frame(RankFirstPos= getRankFirstPos(scoresDf),
                    PrecisionTopPairs= getPrecisionTopPairs(scoresDf,numPairs=500),
                    AUC_ROC= getAUC_ROC(scoresDf,numPairs=500)
                  )
         )
}

getFullComparation <- function(scoresDf1,scoresDf2,
                               numPairs=500,numPairs2=numPairs){
  return(data.frame(RankFirstPos1= getRankFirstPos(scoresDf1),
                    RankFirstPos2= getRankFirstPos(scoresDf2),
                    Precision1= getPrecisionTopPairs(scoresDf1,numPairs=numPairs),
                    Precision2= getPrecisionTopPairs(scoresDf2,numPairs=numPairs2),
                    AUC_ROC1= getAUC_ROC(scoresDf1,numPairs=numPairs),
                    AUC_ROC2= getAUC_ROC(scoresDf2,numPairs=numPairs2)
  )
  )
}