#' @export

fullchronologyC= function(folder,DataP,DataC,resolution=200,supp_type=1,Sample_year=2017,cc=1,ccpb=0,
                         memory_shape=4., memory_mean=.7,acc_shape=1.5,acc_mean=20,
                         fi_mean=50,fi_acc=2,As_mean=20,As_acc=2,w_plot=.8,main1 =""){
  par(mfrow=c(1,1))
  Ages=read.table(paste(folder,"Results/dates.csv",sep=""),sep=" ")
  Ages=Ages+1950-Sample_year
  intervals=read.table(paste(folder,"Results/intervals.csv",sep=""),sep=",")
  intervals=intervals+1950-Sample_year
  Depths=as.numeric(read.table(paste(folder,"Results/depths.csv",sep=""),sep=",") )
  Output=read.table(paste(folder,"Results/Results_output.csv",sep=""),sep=",")
  Lead=read.table(paste(folder,DataP,sep=""),sep=",")
  Carbon=read.table(paste(folder,DataC,sep=""),sep=",")
  Plotval=read.table(paste(folder,"Results/Graphs.csv",sep=""),sep=",")
  Slopes=read.table(paste(folder,"Results/Slopes.csv",sep=""),sep=",")
  num_var=length(Output[0,])
  if(length(Lead[1,])==5){
    Ran=1  } else if(length(Lead[1,])==7){
      if(supp_type==1){Ran=1}else{
      Ran=length(Lead[,1])}}

  maxA=max(Ages[,length(Ages[1,])])+.10
  ageSeq=seq(from=0,to=maxA,maxA/resolution)
  deptsSeq=seq(from=0,to=Depths[length(Depths)],Depths[length(Depths)]/resolution)
  deptsSeq=deptsSeq
  diffSep=(deptsSeq[2]-deptsSeq[1])/2
  TotSeq=length(Ages[,1])

  layout(matrix(c(1,1,2,2,3,3,4,4,
                  1,1,2,2,3,3,4,4,
                  5,5,5,5,5,5,5,5,
                  5,5,5,5,5,5,5,5,
                  5,5,5,5,5,5,5,5,
                  5,5,5,5,5,5,5,5,
                  5,5,5,5,5,5,5,5), 7, 8, byrow = TRUE))


  d <- density(as.numeric(Output[-1,(Ran+2)]))
  plot(d, xlab="",main="Memory",ylab = "",xlim=c(0,1),xaxs="i",yaxs="i")
  polygon(d, col=gray(.6))
  memory_shape2=(memory_shape*(1-memory_mean) )/memory_mean
  lines(seq(0,1,.01),dbeta(seq(0,1,.01),memory_shape,memory_shape2),col="green")

  d <- density(as.numeric(unlist(Output[-1,-c(c(1:(Ran+2)),num_var)])))
  plot(d,xlab="",main="Acc",ylab = "",xaxs="i",yaxs="i")
  polygon(d, col=gray(.6))
  lines(seq(0,100,.5),dgamma(seq(0,100,.5),shape=acc_shape,scale=acc_mean/acc_shape),col="green")

  if(Ran==1){
  d <- density(as.numeric(Output[-1,2]))
  plot(d,xlab="",main="Supported Act",ylab = "",xaxs="i",yaxs="i")
  polygon(d, col=gray(.6))
  lines(seq(0,100,.05),dgamma(seq(0,100,.05),shape=As_acc,scale=As_mean/As_acc),col="green")
  }else {
    min1=min(as.numeric(unlist(Output[-1,2:(Ran+1)])))+.1
    max1=max(as.numeric(unlist(Output[-1,2:(Ran+1)])))+.1
    plot(-10,-10,xlim=c(Lead[1,1],Lead[Ran,1]),ylim=c(min1,max1),xlab="Depth (cm)",ylab="",main="Supported 210Pb",xaxs="i",yaxs="i")
    for (k in 1:Ran) {
      points(rep(Lead[k,1],length(Output[-1,1+k])), Output[-1,1+k],pch=19,col=rgb(0,0,0,.03) )
      points(Lead[k,1],Lead[k,6],col="red",pch=18)

    }
  }

  d <- density(as.numeric(Output[-1,1]))
  plot(d,xlab="",main="Supply of 210Pb",ylab = "",xaxs="i",yaxs="i")
  polygon(d, col=gray(.6))
  lines(seq(0,350,.05),dgamma(seq(0,350,.05),fi_acc,scale=fi_mean/fi_acc),col="green")

  chronologylinesC(folder,DataP,DataC,Sample_year,cc,ccpb,w = w_plot,main1 = main1)


  par(mfrow=c(1,1))
}





#' @export
chronologylinesC= function(folder,DataP,DataC,Sample_year=2017,cc=1,ccpb=0,main1=T,w=.8){
  Ages=read.table(paste(folder,"Results/dates.csv",sep=""),sep=" ")
  Ages=Ages+1950-Sample_year
  intervals=read.table(paste(folder,"Results/intervals.csv",sep=""),sep=",")
  intervals=intervals+1950-Sample_year
  Depths=as.numeric(read.table(paste(folder,"Results/depths.csv",sep=""),sep=",") )
  Output=read.table(paste(folder,"Results/Results_output.csv",sep=""),sep=",")
  Lead=read.table(paste(folder,DataP,sep=""),sep=",")
  Carbon=read.table(paste(folder,DataC,sep=""),sep=",")
  Plotval=read.table(paste(folder,"Results/Graphs.csv",sep=""),sep=",")
  Slopes=read.table(paste(folder,"Results/Slopes.csv",sep=""),sep=",")
  num_var=length(Output[0,])
  iterations=length(Ages[,1])
  
  if(main1==T){main1=gsub("\\..*","",DataP)}


plot(Depths,c(1950-Sample_year,Ages[2,]),type="l",col=rgb(0,0,0,.01), lwd=2,
     ylim = c(1950-Sample_year,max(Ages[,length(Ages[1,])])),
     xlab = "Depth (cm)",ylab="cal yr BP",main= main1)
for (i in 1:(iterations-1)){
  lines(Depths,c(1950-Sample_year,Ages[i,]),type="l",col=rgb(0,0,0,.01), lwd=2)
}
lines(Depths,c(1950-Sample_year,intervals[,2]),type="l", lty=2, lwd=1,col="red")
lines(Depths,(c(1950-Sample_year,intervals[,4])),type="l", lty=2, lwd=1,col="red")
lines(Depths,(c(1950-Sample_year,intervals[,3])),type="l", lty=2, lwd=1,col="red")

for (i in 1:length(Lead[,1])){
#  print(Lead[i,1]+Lead[i,5])
  rect(xleft = Lead[i,1]-Lead[i,5],ybottom =1950-Sample_year-10, xright = Lead[i,1],ytop =  1950-Sample_year,col=rgb(0,0,1,.8))
  
  }



for (i in 1:length(Carbon[,1])){
  plot14C(cdate = as.numeric(Carbon[i,]),cc = cc,ccpb = ccpb,S_year = Sample_year,w = w)
}


}



#' @export
slopes= function(folder,Data){
    Depths=as.numeric(read.table(paste(folder,"Results/depths.csv",sep=""),sep=",") )
    Slopes=read.table(paste(folder,"Results/Slopes.csv",sep=""),sep=",")
    iterations=length(Slopes[,1])
  maxS=max(Slopes)+.10
  plot(-1000,-1000, xlim=c(Depths[length(Depths)]),ylim = c(0, maxS),xlab = "Depth (cm)",ylab="Slopes (Accumulations)" )
  for (i in 1:(iterations-1)){
    lines(Depths,as.numeric(c(Slopes[i,])),type="l",col=rgb(0,0,0,.01), lwd=2)
  }
  for (i in 1:length(Lead[,1])){
    rug(Lead[i,1], lwd = Lead[i,5],col = "blue")
  }

}

#' @export
ageof=function(x,interval=.95,folder,Data){
  if(x!=0){
    Ages=read.table(paste(folder,"Results/dates.csv",sep=""),sep=" ")
    Depths=as.numeric(read.table(paste(folder,"Results/depths.csv",sep=""),sep=",") )
    Slopes=read.table(paste(folder,"Results/Slopes.csv",sep=""),sep=",")
    depfix=which(Depths<x)
    depfix=depfix[length(depfix)]
    m2=Slopes[,depfix]
    sumages=c()

    if(depfix!=1){
      for (i in 1:(length(Ages[,1])-1)){
        sumages=c(sumages, Ages[i,(depfix-1)] + Slopes[i,depfix]* (x-Depths[depfix]) )

      }
    }else{sumages= Slopes[,depfix]* (x)}
    n=length(Ages[,1])
    mean1=mean(sumages)
    inter=(1-interval)/2
    lim1=sort(sumages)[as.integer(n*inter)]
    lim2=sort(sumages)[as.integer(n*(inter+interval))]
    d<- density(sumages)

    plot(d, main=paste("Age at ",x, "cm"),ylab = "",xlab = "Age",xaxs="i",yaxs="i")
    polygon(d, col=gray(.6))
    points(x=lim1,y=(max(d$y)*.015),col="red",pch="(")
    points(x=lim2,y=(max(d$y)*.015),col="red",pch=")")
    points(x=mean1,y=(max(d$y)*.015),col="red",pch=16)
    print(paste("the age of depth",x,"is between") )
    print(c(lim1,lim2))
    print(paste("with a ", interval, "%"," condifence interval and a mean of:",sep = ""))
    print(mean1)


  }else{
    print("For depth 0, the age is equal to the collection date.")
  }

  return(c(lim1,mean1,lim2))

}


require(CPlum)
plot14C <- function(cdate,cc,ccpb,S_year,w=.8){
  print(paste("Plotting",cdate[1],"at depth",cdate[3]) )
  library(rPython)
  modirec=path.package("CPlum", quiet = T)
  ccdir=paste(modirec,"/Calibration Curves/",sep="")
  python.load(paste(modirec,"/","CaPb.py",sep="") )
  Xs=c(0,0)
  if(cdate[1]<=0){
    if(ccpb!=5){postlim=-59.62}else{postlim=-61.19}
    Xs[1]=postlim
    Xs[2]=-.01
  }else{Xs=python.call( "invlookup",cdate[1],cdate[2],cc,ccpb,ccdir)
  Xs[1]=0.01
  }
  
  Xs=seq(Xs[1],Xs[2],length.out = 150)
  Ys=c()
  for (ic in Xs){
    Ys=c(Ys,python.call( "Calibrate",ic,cdate,cc,ccpb,paste0(modirec,"/"),TRUE) )
    #Calibrate(ic,cdate[-3],cc,ccpb))
  }
#  Ys=Ys#-(min(Ys))
  Ys=(Ys/max(Ys))*w
  
  #S_year=1950-S_year
  
  Ys0=((cdate[3]- Ys))
  Ys1=((cdate[3]+ Ys))
  
  polygon(x=c(Ys0,rev(Ys1)),c(Xs,rev(Xs)),  col = rgb(0,0,1,.5), lty = 2, lwd = 2,border=FALSE)
  #polygon((cdate[3]-Ys+Ys[1]),Xs ,  col = rgb(0,0,1,.5), lty = 2, lwd = 2,border=FALSE)
 # return(Ys0)
}



