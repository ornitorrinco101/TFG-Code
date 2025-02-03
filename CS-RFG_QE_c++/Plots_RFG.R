library(latex2exp)



dataRFG1=read.table("/Users/mcytron/TFG-CODE/RFG_QE_c++/Ei1000_C/Ei1000_th15.out",header=FALSE)
dataRFG2=read.table("/Users/mcytron/TFG-CODE/RFG_QE_c++/Ei1000_C/Ei1000_th90.out",header=FALSE)

names(dataRFG1)=c("omega","kf","d3sig")
names(dataRFG2)=c("omega","kf","d3sig")


#dev.new(width=4, height=4, unit="in")

#par(mfrow=c(1,3))
fixlooks=function(k,sigma){
  l=length(k)
  if(k[1]>k[l]){
    extend=seq(min(k),0)
    ktot=c(k,extend)
    stot=c(sigma,numeric(length(extend)))
    return(cbind(ktot,stot))
  }else{
    extend=seq(0,min(k))
    ktot=c(extend,k)
    stot=c(numeric(length(extend)),sigma)
    return(cbind(ktot,stot))
  }
}

vals=fixlooks(dataRFG1$kf,dataRFG1$d3sig)

kftot=vals[,1]
stot=vals[,2]


plot(kftot/1000,stot,
     xlab=TeX("$k_{f}\\;(GeV)$"),
     ylab=TeX("$d^{3}\\sigma/dk_{f}\\,d\\Omega_{f}\\;\\;(fm^{2}\\cdot MeV^{-1}\\cdot sr^{-1})$"),
     cex.lab=0.75,cex.axis=0.75,
     type="l",
     col=2)

plotRFG=function(first,last,step,color){
  n_vals=seq(first,last,step)
  data=list()
  i=1
  
  for (n in n_vals) {#CARBON
    # Construct the file name
    fn=paste0("/Users/mcytron/TFG-CODE/RFG_QE_c++/Ei1000_C/Ei1000_th", n, ".out")
    
    # Read the data from the file (adjust the path as needed)
    data[[n]]=read.table(fn, header = FALSE)
    u=data[[n]]
    names(u)=c("omega","kf","d3sig")
    vals=fixlooks(u$kf,u$d3sig)
    kftot=vals[,1]
    stot=vals[,2]
    lines(kftot/1000,stot,col=color[i])
    i=i+1
  }
}

plotRFG(30,90,15,color = 3:7)

legend(x="topleft",legend=c(TeX("$\\theta=15^{o}$"),TeX("$\\theta=30^{o}$"),TeX("$\\theta=45^{o}$"),TeX("$\\theta=60^{o}$"),TeX("$\\theta=75^{o}$"),TeX("$\\theta=90^{o}$")),lty=c("solid","solid","solid","solid","solid","solid"),col=c(2:7))
mtext(TeX("$Sección\\;eficaz\\;diferencial\\;del\\;proceso\\;a\\;\\epsilon_{i}=1\\;GeV"),side=3,las=1,line=0.8,cex=1)


plotit=function(element,N,first,last,step,color){
  n_vals=seq(first,last,step)
  data=list()
  d2sig=numeric(length(n_vals))
  i=1
  
  for (n in n_vals) {#CARBON
    # Construct the file name
    fn=paste0("/Users/mcytron/TFG-CODE/RFG_QE_c++/Ei1000_",element,"/Ei1000_th", n, ".out")
    
    # Read the data from the file (adjust the path as needed)
    data[[n]]=read.table(fn, header = FALSE)
    u=data[[n]]
    names(u)=c("omega","kf","d3sig")
    d2sig[i]=(10**12)*sum(u$d3sig)*2*pi/N
    i=i+1
  }
  lines(n_vals,d2sig,col=color)
}


data1000=read.table("/Users/mcytron/TFG-CODE/lepton-nucleon-c++/Ei1000.out",header=FALSE)
names(data1000)=c("omega","K_f","P_n","Q2","dsigQ2","costheta_f","dsig")


deg=acos(data1000$costheta_f)*180/pi
good=deg<=90

plot(deg[good],data1000$dsig[good]*10**12,
     xlab=TeX("$\\theta\\;(deg)$"),
     ylab=" ",
     cex.lab=1,
     type="l",
     col="royalblue")

plotit("C",6,1,30,1,"hotpink3")
plotit("C",6,50,90,5,"hotpink3")
plotit("C",6,30,50,2,"hotpink3")

plotit("D",1,20,90,5,"firebrick3")
plotit("D",1,1,20,1,"firebrick3")

plotit("3He",1,20,90,5,"coral1")
plotit("3He",1,1,20,1,"coral1")

plotit("Ar",22,30,90,5,"darkolivegreen")
plotit("Ar",22,1,30,1,"darkolivegreen")

legend(x="topright",
       legend=c("Dispersión elástica",TeX("$QE\\;RFG\\phantom{o}^{2}H$"),TeX("$QE\\;RFG\\phantom{o}^{3}He$"),TeX("$QE\\;RFG\\phantom{o}^{12}C$"),TeX("$QE\\;RFG\\phantom{o}^{56}Fe$")),
       col=c("royalblue","firebrick3","coral1","hotpink3","darkolivegreen"),lty=c("solid","solid","solid","solid","solid"),cex=0.75)

mtext(TeX("$d\\sigma/d\\cos(\\theta)\\;\\;(10^{-38}\\;cm^{2}/Nucleón)$"),side=2,line=2,cex=1)
#mtext(TeX("$Sección\\;eficaz\\;diferencial\\;del\\;proceso\\;a\\;\\epsilon_{i}=1\\;GeV"),side=3,las=1,line=0.8,cex=1)


