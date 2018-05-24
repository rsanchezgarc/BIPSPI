
args <- commandArgs(trailingOnly = TRUE)
if(length(args)==0){
  resPath<-"~/Tesis/rriPredMethod/data/bench5Data/newCodeData/results_xgb/mixed/"
}else{
  resPath<-args[1]
}

if(length(args)==2){
  precisionAt<- as.integer(args[2])
}else{
  precisionAt<- 500
}

numMaxRows<-1024

if(precisionAt> numMaxRows){
  numMaxRows<- precisionAt
}

fnames<- list.files(resPath,full.names = TRUE, pattern="tab$")

listOfDf<-lapply( fnames,
                read.table,header=T,colClasses=c(rep("NULL",6),"integer","numeric"))

numberContacts<-sapply(listOfDf,function(x){sum(x$categ==1)})

minNumRowsInDf<- min(sapply(listOfDf,nrow))
if (numMaxRows>minNumRowsInDf){
  numMaxRows<-minNumRowsInDf
}

listOfCategVects<-lapply(listOfDf,
                  function(x){
                    x<-x[order(x$prediction,decreasing = T),]
                    if(numMaxRows<= nrow(x)){
                      categs<-x$categ[1:numMaxRows]  
                    }else{
                      categs<-x$categ
                    }
                    
                    categs[categs != 1]<-0
                    return(categs)
                  }
                )
numComplexes<-length(listOfDf)
categMatrix<-do.call(cbind,listOfCategVects)
desiredRffpPs<- c(10, 25, 50, 75, 90)
desiredRffpDif<- rep(99, length(desiredRffpPs))
rffpP_val<- rep(-1, length(desiredRffpPs))

currNumPairs<- 1
precision<- -1
cat("RFPP(f)\tf(%)\n")
for(row in currNumPairs:numMaxRows){
  subMatrix<- categMatrix[1:row,, drop=F]
  positiveHit<- colSums(subMatrix)
  meanNumHits<-mean(positiveHit)
  positiveHit<- colSums(subMatrix) >0
  percentAtI<-(sum(positiveHit)/numComplexes)*100
  if ( currNumPairs  <= row){
    cat(row," -->",percentAtI,"%. Mean num hits",meanNumHits,"\n")
    currNumPairs<-2*currNumPairs
  }
  if(row== precisionAt){
    precision<- colSums(subMatrix) / precisionAt
    precisionOverContacts<-precision / numberContacts
  }
  for( i in 1:length(desiredRffpPs)){
    dif_val<- percentAtI- desiredRffpPs[i]
    if ( abs(dif_val) < abs(desiredRffpDif[i]) ){
      rffpP_val[i]<- row
      desiredRffpDif[i]<- dif_val
    }
  }
}


if(row<precisionAt){
  precisionAt<- row
  precision<- colSums(subMatrix) / precisionAt
  precisionOverContacts<-precision / numberContacts
}

cat("mean (Precision(best",precisionAt,"))",mean(precision),"\n")
cat("mean (Precision(best",precisionAt,")/numberOfcontacts)", mean(precisionOverContacts),"\n")

cat("rfpp(p)\n")
cat(desiredRffpPs+desiredRffpDif)
cat("\n---------------------\n")
cat(rffpP_val)
cat("\n")


