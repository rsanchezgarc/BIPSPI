# Se encarga de obtener una lista de data frame con las diversas interacciones 
# contenidas en un data frame
preproc<-function(scoresDf){
  
  removeBadAA<- function(scoreDf){
    # print("fixing scoreDf")
    i<-1
    numElem<-c()
    # print (sort(unique(scoreDf$resIdL)))
    for(lIndex in sort(unique(scoreDf$resIdL))){
      numElem[i]<-nrow(scoreDf[scoreDf$resIdL==lIndex,])
      i<-i+1
    }
    # print(numElem)
    numDifAA<-unique(numElem)
    # print(numDifAA)
    if(length(numDifAA)>1){
      # count<-table(numDifAA)
      count<-table(numElem)
      print("Some mismatch")
      # print(count)
#       print(which.min(count))
      # print(names(which.min(count)))
      badIndexes<- which(numElem==as.integer(names(which.min(count))))
      badLIndex<- sort(unique(scoreDf$resIdL))[badIndexes]
      for(i in badLIndex){
        scoreDf<-scoreDf[!scoreDf$resIdL==i,]
      }

    }
#     print(sort(unique(scoreDf$resIdL)))
#     print(sort(unique(scoreDf$resIdR)))
    return(scoreDf)
  }
  
  chainsR<-unique(scoresDf$chainIdR)
  chainsL<-unique(scoresDf$chainIdL)
  dataFrameList<-list()
  k<-1
  for(r in chainsR){
    # print(r)
    for(l in chainsL){
      # print(l)
      temDf<- scoresDf[scoresDf$chainIdL==l & scoresDf$chainIdR==r,]
#       print("antes")
#       print(dim(temDf))
      temDf<- removeBadAA(temDf)
#       print("despues")
#       print(dim(temDf))
      # if (sum(temDf$categ==1)==0){next}
      dataFrameList[[k]]<- temDf[order(temDf$resIdL,
                                       temDf$resIdR),]
      
      k<- k +1
    }
  }
  # print("preproc done")
  return(dataFrameList)
}

getScoreMatrixFromDf<-function(scoreDf){
  scoreDf<- scoreDf[order(scoreDf$resIdL,
                    scoreDf$resIdR),]
  lenSeq1<-length(unique(scoreDf$resIdL))
  scoreMat<-matrix(scoreDf$prediction,
                   nrow=lenSeq1,byrow=TRUE)
  return(scoreMat)
  
}

getActualMatrixFromDf<-function(scoreDf){
  scoreDf<- scoreDf[order(scoreDf$resIdL,
                          scoreDf$resIdR),]
  lenSeq1<-length(unique(scoreDf$resIdL))
  actualMat<-matrix(scoreDf$categ,
                    nrow=lenSeq1,byrow=TRUE)
  return(actualMat)
  
}

getDfFromScoreMatrix<-function(scoreMatrix,origDf){
  
#   print("Computing again df")
#   print(dim(scoreMatrix))
#   print(nrow(origDf))
  origDf<- origDf[order(origDf$resIdL,
                          origDf$resIdR),]
  prediction<-c(t(scoreMatrix))
  #Esto es una prueba
  # origDf$prediction<-prediction
  origDf$prediction<- origDf$prediction + prediction
  return(origDf)
  
}

