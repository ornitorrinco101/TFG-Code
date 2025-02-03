
#nd280data-numu-cc0pi-xs-on-c-2015

library(latex2exp)

#nd280data-numu-cc0pi-xs-on-c-2015

rm=function(a){return(a<0.5)}



data1=read.table("/Users/mcytron/TFG-CODE/DATA-EXP/cross-section_analysisI.txt",header=TRUE)
data2=read.table("/Users/mcytron/TFG-CODE/DATA-EXP/rps_crossSection_analysis2.txt",header=TRUE)

data1=cbind(data1,apply(cbind(data1[,2],data1[,3]),1,mean),apply(cbind(data1[,4],data1[,5]),1,mean))
data2=cbind(data2,apply(cbind(data2[,2],data2[,3]),1,mean),apply(cbind(data2[,4],data2[,5]),1,mean))

names(data1)[7:8]=c("costheta","p")
names(data2)[7:8]=c("costheta","p")

c1=data1$costheta
cmin1=data1$cosThetamin
cmax1=data1$cosThetamax

x1=data1$xsec

p1=data1$p
pmin1=data1$momentummin
pmax1=data1$momentummax

c2=data2$costheta
cmin2=data2$cosThetamin
cmax2=data2$cosThetamax

x2=data2$xsec

p2=data2$p
pmin2=data2$momentummin
pmax2=data2$momentummax

err11=read.table("/Users/mcytron/TFG-CODE/DATA-EXP/covariance_fluxNormalizationSystematics_analysisI.txt",header=TRUE)
err11.f=factor(err11$i)
err11=tapply(err11$cov,err11.f,sum)

err12=read.table("/Users/mcytron/TFG-CODE/DATA-EXP/covariance_shapeSystematics_analysisI.txt",header=TRUE)
err12.f=factor(err12$i)
err12=tapply(err12$cov,err12.f,sum)

err13=read.table("/Users/mcytron/TFG-CODE/DATA-EXP/covariance_statisticUncertainty_analysisI.txt",header=TRUE)
err13.f=factor(err13$i)
err13=tapply(err13$cov,err13.f,sum)

err1=err11+err12+err13

plotz1=function(deg,min,max,color){
  j=0
  dsig=numeric(80)
  kf=numeric(80)
  for(i in deg){
    fn=paste0("/Users/mcytron/TFG-CODE/RFG_QE_c++/con-flujo/th",
              i,".out")
    d=read.table(fn,header=FALSE)
    names(d)=c("kf","kfp","d3sig")
    dsig=dsig+d$d3sig
    kf=d$kf
    j=j+1
  }

  u=(cmax1==max&cmin1==min&(pmax1-pmin1)<2)
  
  kf=kf/(1000)
  dsig=dsig/j
  m=paste(min,TeX("$\\leq\\cos\\theta\\leq$"),max)
  
  plot(kf,dsig/10,type="l",lty="dotted",main=TeX(paste(min,"$\\leq\\cos\\theta\\leq$",max)),
       xlab="",ylab="",cex.main=2) #10-39cm^2 ->10-38 cm^2
  title(xlab=TeX("$k_{f}\\;(GeV)$"),cex.lab=1.5)
  title(ylab=TeX("$d\\sigma/dk_{f} d\\cos \\theta\\; (10^{-38} cm^{2}\\cdot (GeV\\cdot Nucleon)^{-1})$"),
        line=2, cex.lab=1.5)
  points(p1[u],x1[u],col=color,pch=20,cex=0.5)
  
  legend(x="topright",legend=c("RFG","T2K"),
         col=c("black",color),lty=c("dotted",NA),
         pch=c(NA,20))
  
  arrows(pmin1[u], x1[u], pmax1[u], x1[u], code=3,angle=90, length=0.02, col=color)
  arrows(p1[u],x1[u]-err1[u],p1[u],x1[u]+err1[u], code=3,angle=90, length=0.02, col=color)
}
par(mfrow=c(2,3))

#2)
#plotz1(72,0,0.6,"hotpink3")


#3)
plotz1(c(49),0.6,0.7,"hotpink3")
#4)
plotz1(c(41),0.7,0.8,"hotpink3")
#5)
plotz1(34,0.8,0.85,"hotpink3")
#6)
plotz1(28,0.85,0.9,"hotpink3")
#7)
plotz1(23,0.9,0.94,"hotpink3")
#8)
plotz1(16,0.94,0.98,"hotpink3")
#9)
#plotz1(8,0.98,1,"hotpink3")




err21=read.table("/Users/mcytron/TFG-CODE/DATA-EXP/rps_fluxNormCov_analysis2.txt",header=TRUE)
err21.f=factor(err21$i)
err21=tapply(err21$cov,err21.f,sum)

# err22=read.table("/Users/mcytron/TFG-CODE/DATA-EXP/rps_statsCov_analysis2.txt",header=TRUE)
# err22.f=factor(err22$i)
# err22=tapply(err22$cov,err22.f,sum)

# err23=read.table("/Users/mcytron/TFG-CODE/DATA-EXP/rps_systCov_analysis2.txt",header=TRUE)
# err23.f=factor(err23$i)
# err23=tapply(err23$cov,err23.f,sum)

# err2=err21+err22+err23

err2=err21

plotz2=function(deg,min,max,color){
  fn=paste0("/Users/mcytron/TFG-CODE/RFG_QE_c++/con-flujo/th",
            deg,".out")
  dataC=read.table(fn,header=FALSE)
  names(dataC)=c("kf","kfp","d3sig")
  u=(cmax2==max&cmin2==min&(pmax2-pmin2)<2)
  kf=dataC$kfp/1000
  m=paste(min,TeX("$\\leq\\cos\\theta\\leq$"),max)
  
  plot(kf,dataC$d3sig/10,type="l",lty="dotted",main=TeX(paste(min,"$\\leq\\cos\\theta\\leq$",max)),
       xlab=TeX("$k_{f}\\;(GeV)$"),ylab="",cex.lab=1.5,cex.main=2) #10-39cm^2 ->10-38 cm^2
  title(ylab=TeX("$d\\sigma/dk_{f} d\\cos \\theta\\; (10^{-38} cm^{2}\\cdot (GeV\\cdot Nucleon)^{-1})$"),
        line=2, cex.lab=1.5)
  points(p2[u],x2[u],col=color,pch=20,cex=0.5)
  
  legend(x="topright",legend=c("RFG","T2K"),
         col=c("black",color),lty=c("dotted",NA),
         pch=c(NA,20))
  
  arrows(pmin2[u], x2[u], pmax2[u], x2[u], code=3,angle=90, length=0.02, col=color)
  arrows(p2[u],x2[u]-err2[u],p2[u],x2[u]+err2[u], code=3,angle=90, length=0.02, col=color)
}



par(mfrow=c(2,3))


#2)
plotz2(41,0.7,0.8,"royalblue")
#3)
plotz2(34,0.8,0.85,"royalblue")
#4)
plotz2(28,0.85,0.9,"royalblue")
#5)
plotz2(24,0.9,0.925,"royalblue")
#6)
plotz2(20,0.925,0.95,"royalblue")
#7)
plotz2(15,0.95,0.975,"royalblue")
#8)
#plotz2(9,0.975,1,"royalblue")


