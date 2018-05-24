suppressMessages(library(lattice))
source("workers/manageDfToMatrices.R")
source("workers/evaluatePerformance.R")


dataPath<- "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/productVersionData/results/mixed_2"
fotosPath<- "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/productVersionData/results/_fotos_"


plotOneMatrix<-function(x,mess=""){
  plot<-levelplot(t(x),
                   main=mess,aspect="fill",
                   colorkey=NULL)
  return(plot)
}

plotActualVsPredMatrix<-function(xDf, pairsNum= 256,addMess=""){
  
  x<-xDf[order(xDf$prediction,decreasing = T),]
  if(nrow(x)<= pairsNum){
    pairsNum<-nrow(x)
  }
  x[(x$categ==0 | x$categ==-1),"prediction"]<-0
  x[x$categ==1,"prediction"]<-1
  x[1:pairsNum,"prediction"]<-2
  precision<-getPrecisionTopPairs(x,pairsNum)
  x[which(x[1:pairsNum,"categ"]==1),"prediction"]<-3
  predMatrix<-getScoreMatrixFromDf(x)
  if(precision==0){
    myColours<-colorRampPalette(c("white","blue","lightgrey"))
  }else{
    myColours<-colorRampPalette(c("white","blue","lightgrey","green"))  
  }


  plot1<-levelplot(t(predMatrix),
                   main=paste(addMess,"Prediction matrix. Pecision(1st ",pairsNum,")= ",precision,sep=""),
                   xlab="receptor seq",
                   ylab="ligand seq",
                   aspect="fill",
                   col.regions = myColours,
                   colorkey=NULL
  )

  print(plot1, split = c(1, 1, 1, 1), more = FALSE)
  
}

plotearUnDfWithDensities <- function(xDf,pairsNum=500){
  scatterBarNorm <- function(x,realX,xlimits,ylimits, dcol="blue",lhist=20, num.dnorm=5*lhist, ...){
    ## check input

    stopifnot(ncol(x)==2)
    ## set up layout and graphical parameters
    layMat <- matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
    layout(layMat, widths=c(5/7, 2/7), heights=c(2/7, 5/7))
    ospc <- 0.5 # outer space
    pext <- 4 # par extension down and to the left
    bspc <- 1 # space between scatter plot and bar plots
    par. <- par(mar=c(pext, pext, bspc, bspc),
                oma=rep(ospc, 4)) # plot parameters
    ## scatter plot
    plot(x, xlim=xlimits, ylim=ylimits, ...)
    points(realX, xlim=xlimits, ylim=ylimits, col="red",...)
    ## 3) determine barplot and height parameter
    ## histogram (for barplot-ting the density)
    xhist <- hist(x[,1], plot=FALSE, breaks=seq(from=xlimits[1], to=xlimits[2],
                                                length.out=lhist))
    yhist <- hist(x[,2], plot=FALSE, breaks=seq(from=ylimits[1], to=ylimits[2],
                                                length.out=lhist)) # note: this uses probability=TRUE
    ## determine the plot range and all the things needed for the barplots and lines
    xx <- seq(min(x[,1]), max(x[,1]), length.out=num.dnorm) # evaluation points for the overlaid density
    xy <- dnorm(xx, mean=mean(x[,1]), sd=sd(x[,1])) # density points
    yx <- seq(min(x[,2]), max(x[,2]), length.out=num.dnorm)
    yy <- dnorm(yx, mean=mean(x[,2]), sd=sd(x[,2]))
    ## barplot and line for x (top)
    par(mar=c(0, pext, 0, 0))
    barplot(xhist$density, axes=FALSE, ylim=c(0, max(xhist$density, xy)),
            space=0) # barplot
   
    ## barplot and line for y (right)
    par(mar=c(pext , 0, 0+1, 0))
    barplot(yhist$density, axes=FALSE, xlim=c(0, max(yhist$density, yy)),
            space=0, horiz=TRUE) # barplot

    ## restore parameters
    par(par.)
  }
  
  numPosContacts<-sum(xDf$categ==1)
  x<-xDf[order(xDf$prediction,decreasing = T),]
  x[(x$categ==0 | x$categ==-1),"prediction"]<-0
  x[x$categ==1,"prediction"]<-1
  x[1:pairsNum,"prediction"]<-2
  x[which(x[1:pairsNum,"categ"]==1),"prediction"]<-3
  predMatrix<-getScoreMatrixFromDf(x)
  predictedIndexes<-which(predMatrix==2,arr.ind = T)
  # x11()
  realIndexes<-which(predMatrix==1 | predMatrix==3 ,arr.ind = T)
  xlimits<-c(1,(dim(predMatrix)[1]))
  ylimits<-c(1,(dim(predMatrix)[2]))
  scatterBarNorm(predictedIndexes,realIndexes,xlimits,ylimits,lhist = min(dim(predMatrix)))
  
  # points(realIndexes, pch=2,col= "red")
}


plotContactMatricesFromDirectory<-function(pathData=dataPath,outPath=fotosPath){
  fileNames<- list.files(pathData)
  fileFullNames<- list.files(pathData,full.names = T)
  resList<-list()
  k<-1
  for (fname in fileFullNames){
    print(fname)
    xDf<-read.table(fname,header=T)
    dataFramesList<-preproc(xDf)
    j<-1
    for(df in dataFramesList){
      png(filename = file.path(outPath,paste(fileNames[k],j,".png",sep="")),
          width= 1024,height = 768)
      plotActualVsPredMatrix(df)
      dev.off()
      j<-j+1
    }
    k<-k+1
  }
  
}


 
plotContactMatricesFromDirectory()
