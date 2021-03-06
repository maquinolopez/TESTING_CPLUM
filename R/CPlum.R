#' MCMC to obtain the posterior distribution of the Chronology.
#'
#' It takes the Data of eat each depth and creates samples from the
#' posterior distribution using the Twalk methodology.
#'
#' @param itera is the number of iterations that Twalk will do (1e+6 is the recomended)
#' @param Datos is the activity at each depth and organize starting from the surface.
#' @param sdDatos is the standar deviation of each activity.
#' @param Conti is a vector showing which data was sample discontinuosly (1 would mean the data was sample discontinusly and 0 means the data was sampled continously).
#' @param Bqkg is True when the data is in the form Bq/kg, if it False the data would be consider as ...

#' @export

runCPlum=function(folder,DataPb,DataC,iterations=2e+3,by=5.0,number_supported=FALSE,detection_limit=.05,
                memory_shape=4., memory_mean=.4,fi_mean=50,fi_acc=2,filediv=1,
                As_mean=20,As_acc=2,resolution=200,burnin=1000,thi=10,Ratype=F,
                 acc_shape=1.5,acc_mean=20,cc=1,ccpb=0,Sample_year=2017,seeds=12345678){
  library(rPython)
  folder=paste(normalizePath(folder),"/",sep="")
  Lead=read.table(paste(folder,DataPb,sep=""),sep=",")

  
  print("the length")
  print(by)
  if (as.character(by)==TRUE){
    by=(Lead[length(Lead[,1]),1])/10#25
  }
  print(by)
  
  
  
  
  
  if(Ratype==F){
    if(number_supported==FALSE){
      if(length(Lead[1,])==5){
        n.check=check.equi(Lead)
        print(n.check[1])
        usrresp=1
        while(!(usrresp=="Yes"||usrresp=="No"||usrresp=="no"||usrresp=="yes")){
          cat("Are you sure these data represent the supported 210Pb for this site?")
          usrresp=readline( "Yes or No \n ")
          if(usrresp=="Yes"||usrresp=="No"||usrresp=="no"||usrresp=="yes"){
            if(usrresp=="Yes"||usrresp=="yes"){
              number_supported=n.check[1]}
            if(usrresp=="No"||usrresp=="no"){
              checksupp="notoy"
              while(typeof(checksupp)!="double"){
                cat("Please indicate how many data points whould be use for estimating the supported 210Pb")
                number_supported=scan("",n=1)
                }
            }
          }
        }
        usemod=1
      }else if(length(Lead[1,])==7){
        cat("You have 226Ra data. \n")
        plot(Lead[,1],Lead[,6],pch=16,ylim=c(min(Lead[,6]-Lead[,7]),max(Lead[,6]+Lead[,7])), ylab="Concentration of 226Ra", xlab="Depth (cm)")
        segments(Lead[,1], Lead[,6]-Lead[,7], x1 = Lead[,1], y1 = Lead[,6]+Lead[,7])
        cat("Plum can assum to have a constant supported 210Pb and use the 226Ra data to infer this one value\n")
        cat("Plum can also assum individual supporeted 210Pb per data point.\n
            It is important to consider that this will greatly increses the computing time and it should only be use when clear parters are observed in the 226Ra data.\n")
        cat("\n If you want to use the constant supported 210Pb press 1, if you want to use the individual 210Pb press 2\n")
        usemod=0
        while(!(usemod==1||usemod==2)){
          usemod=scan("",n=1)
        }
        }
    }else{usemod=1}
  }else {if (Ratype==2){
    usemod=2
  }
    if(Ratype==1){
      usemo=1
      if(number_supported==FALSE){
        cat("Please indicate how many data points whould be use for estimating the supported 210Pb")
        number_supported=scan("",n=1)
      }
    }
    }




modirec=path.package("CPlum", quiet = T)
ccdir=paste(modirec,"/","Calibration Curves/",sep="")
if (usemod==1){
  MCMC=paste(modirec,"/","CaPb.py",sep="")
}else if(usemod==2){
  MCMC=paste(modirec,"/","CaPb_Ra.py",sep="")
  number_supported=0
}

print(MCMC)

python.load(MCMC)
dir.create(paste(folder,"Results",sep = ""))

python.call("runmod",folder,DataPb,DataC,ccdir,FALSE,TRUE, Sample_year,   number_supported   ,    detection_limit   ,
            iterations,  by ,memory_shape     ,memory_mean    ,acc_shape       ,acc_mean,
            fi_mean,fi_acc, As_mean,As_acc,   cc,ccpb,resolution,seeds,burnin,thi)

                     #dirt,plomo,carbon,Dircc,T_mod,T_mod_C,S_year,num_sup,det_lim,
#           iterations, by,shape1_m,          mean_m,       shape_acc,          mean_acc,
#           fi_mean,fi_acc,As_mean,As_acc,cc,ccpb,resolution,seeds
##############


Ages=read.table(paste(folder,"Results/dates.csv",sep=""),sep=" ")
intervals=read.table(paste(folder,"Results/intervals.csv",sep=""),sep=",")
Depths=as.numeric(read.table(paste(folder,"Results/depths.csv",sep=""),sep=",") )
Output=read.table(paste(folder,"Results/Results_output.csv",sep=""),sep=",")
Plotval=read.table(paste(folder,"Results/Graphs.csv",sep=""),sep=",")
Slopes=read.table(paste(folder,"Results/Slopes.csv",sep=""),sep=",")
Carbon=read.table(paste(folder,DataC,sep=""),sep=",")

num_var=length(Output[0,])

maxA=max(Ages[,length(Ages[1,])])+.10
ageSeq=seq(from=0,to=maxA,maxA/resolution)
deptsSeq=seq(from=0,to=Depths[length(Depths)],Depths[length(Depths)]/resolution)
deptsSeq=deptsSeq
diffSep=(deptsSeq[2]-deptsSeq[1])/2
TotSeq=length(Ages[,1])


pdf(paste(folder,'Fi.pdf',sep=""))
plot(as.numeric(Output[-1,1]),type="l",main="fi",xlab="",ylab="")
dev.off()
pdf(paste(folder,'Supported.pdf',sep=""))
plot(as.numeric(Output[-1,2]),type="l",main="Supported Act",xlab="",ylab="")
dev.off()
pdf(paste(folder,'Energy.pdf',sep=""))
plot(as.numeric(Output[-1,num_var]),type="l",main="Energy",xlab="",ylab="")
dev.off()
pdf(paste(folder,'Memory.pdf',sep=""))
plot(as.numeric(Output[-1,3]),type="l",main="Memory",xlab="",ylab="")
dev.off()

pdf(paste(folder,'Fi hist.pdf',sep=""))
hist(as.numeric(Output[-1,1]),breaks=50,probability=T,main="fi",xlab="")
dev.off()
pdf(paste(folder,'Supported hist.pdf',sep=""))
hist(as.numeric(Output[-1,2]),breaks=50,probability=T,main="Supported Act",xlab="")
dev.off()
pdf(paste(folder,'Energy hist.pdf',sep=""))
hist(as.numeric(Output[-1,num_var]),breaks=50,main="Energy",xlab="")
dev.off()
pdf(paste(folder,'Memory hist.pdf',sep=""))
hist(as.numeric(Output[-1,3]),breaks=50,main="Memory",probability=T,xlim=c(0,1),xlab="")
lines(seq(0,100,.01),dbeta(seq(0,100,.01),4,1.5),col="red")
dev.off()


pdf(paste(folder,'acc hist.pdf',sep=""))
hist(as.numeric(unlist(Output[-1,-c(1,2,3,num_var)])),breaks=50,main="Acc",probability=T,xlab="")
lines(seq(0,100,.5),dgamma(seq(0,100,.5),1.5,scale=20/1.5),col="red")
dev.off()




pdf(paste(folder,'Chronology.pdf',sep=""))
fullchronologyC(folder,DataPb,DataC,Sample_year = Sample_year,cc = cc,ccpb = ccpb,
               memory_shape=memory_shape,memory_mean=memory_mean,acc_mean=acc_mean,
               acc_shape=acc_shape,supp_type = usemod)
dev.off()


par(mfrow=c(1,1))
fullchronologyC(folder,DataPb,DataC,Sample_year = Sample_year,cc = cc,ccpb = ccpb,
               memory_shape=memory_shape,memory_mean=memory_mean,acc_mean=acc_mean,
               acc_shape=acc_shape,supp_type = usemod)


intervalfile=by_cm(folder,filediv)
write.csv(x = intervalfile,file = paste0(folder,"ages.csv"))


}

#################Check equilibrium #########################
#' @export
check.equi = function (rawdat){
  rawdata=rawdat[,3]
  rawsd=rawdat[,4]
  deps=rawdat[,1]
  cat("Because linear regression needs at least three data points, this function will start checking for equilibrium starting from the last three data points\n")
  lendat=length(rawdata)
  numdat=as.integer(.5*length(rawdata))
  usedat=rawdata[(lendat-3):lendat]
  usesd=rawsd[(lendat-3):lendat]
  usex=1:length((lendat-3):lendat)
  usereg= lm(usedat ~ usex, weights=1/(usesd^2))
  reg=coef(summary(usereg))[2,4]
  est=coef(summary(usereg))[1,1]
  coe=3
  for (i in 1:numdat){
    usedat=rawdata[(lendat-3-i):lendat]
    usesd=rawsd[(lendat-3-i):lendat]
    usex=1:length((lendat-3-i):lendat)
    usereg= lm(usedat ~ usex, weights=1/(usesd^2))
    reg1=coef(summary(usereg))[2,4]
    est1=coef(summary(usereg))[1,1]
    if(reg1>reg){ reg=reg1;coe=(3+i);est=est1 }
  }

  cat(paste( "The regression process proposes the use of the last", as.integer(coe), " data points as estimates of the supported activity, with a p-value of" ,round((reg),3)) )
  cat(".\n")
    plot(deps,rawdata,pch=16)
  points(deps[(lendat-coe+1):lendat],rawdata[(lendat-coe+1):lendat],col="red",pch=16,xlab="Depth (cm)",ylab="210Pb concentration")
  abline(h=est,lty=2) #mean(rawdata[(lendat-coe+1):lendat])
  return (c(coe,reg1))
}

#' @export
check.qui.Ra = function(rawdat){
  `226Ra`=rawdat[,6]
  rawsd=rawdat[,7]
  deps=rawdat[,1]
  tests=shapiro.test(`226Ra`)
  print(tests)
  if(tests$p.value<=.05){
    cat("Because the p-value of the normality test is smaller than .05 it cannot be assumed that the 226Ra data comes from a single distribution.\n")
    return(2)
  }else{
    cat("Because the p-value of the normality test is bigger than .05 it can be assumed that the 226Ra data comes from a single distribution.\n")
    return(1)
  }
}


#' @export
Calibrate =function(x,cdate,cc,ccpb){
  library(rPython)
  modirec=path.package("CPlum", quiet = T)
  ccdir=paste(modirec,"/",sep="")
  MCMC=paste(modirec,"/","CaPb.py",sep="")
  python.load(MCMC)
  inc=python.call( "incallookup2",x,cc,ccpb,ccdir)
  sigm=(inc[2]^2 + cdate[2]^2)
  mu=inc[1]
  utest=((4. + ((cdate[1]-mu)^2.)/((2.*sigm)) )^(-7./2)  )#/sqrt(sigm)
  return (utest)
}



#' @export
by_cm=function(folder,filediv){
  intv=read.table(paste(folder,"Results/intervals.csv",sep=""),sep=",")
  dpts= seq(filediv,max(intv[,1]),by=filediv)
  agedmodl=approx(c(0,intv[,1]),c(0,intv[,2]),dpts)$y
  agedmodm=approx(c(0,intv[,1]),c(0,intv[,3]),dpts)$y
  agedmodu=approx(c(0,intv[,1]),c(0,intv[,4]),dpts)$y
  Age_max=matrix(c(c(1:max(intv[,1])),agedmodl,agedmodm,agedmodu),ncol = 4,byrow=F)
  colnames(Age_max) <- c("depth","min","mean","max")
  return(Age_max)
}
  