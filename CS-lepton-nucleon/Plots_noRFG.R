library(latex2exp)



data1=read.table("/Users/mcytron/TFG-CODE/lepton-nucleon-c++/Ei1000.out",header=FALSE)
data2=read.table("/Users/mcytron/TFG-CODE/lepton-nucleon-c++/Ei500.out",header=FALSE)
names(data1)=c("omega","K_f","P_n","Q2","dsigQ2","costheta_f","dsig")
names(data2)=c("omega","K_f","P_n","Q2","dsigQ2","costheta_f","dsig")
#data$Q es omega

#dev.new(width=4, height=4, unit="in")

#par(mfrow=c(1,3))
plot(data1$costheta_f,data1$dsig,
     xlab=TeX("$\\cos(\\theta)$"),
     ylab=" ",
     cex.lab=0.75,cex.axis=0.75,
     type="l",
     col="darkorchid")

lines(data2$costheta_f,data2$dsig,
      type="l",
      col="deeppink")

legend(x="topleft",
       legend=c(TeX("$\\epsilon_{i}=1000\\;MeV$"),TeX("$\\epsilon_{i}=500\\;MeV$")),
       col=c("darkorchid","deeppink"),lty=c("solid","solid"),cex=0.6)

mtext(TeX("$d\\sigma/d\\cos(\\theta)\\;\\;(fm^{2})$"),side=2,line=2,cex=0.75)
mtext(TeX("$Sección\\;eficaz\\;diferencial\\;del\\;proceso\\;\\nu_{\\mu}+n\\to\\mu^{-}+p$"),side=3,las=1,line=0.8,cex=1)

#mtext(TeX("$\\frac{d\\sigma(\\theta)}{d\\cos(\\theta)}$"),side=2,las=1,line=2,cex=0.75)
#mtext("\n\n\n\n (fm²)",side=2,las=1,line=2,cex=0.75)


plot(data1$Q2,data1$dsigQ2,
     xlab="",
     ylab="",
     cex.lab=0.75,cex.axis=0.75,
     type="l",
     col="darkorchid")

lines(data2$Q2,data2$dsigQ2,
     type="l",
     col="deeppink")

mtext(TeX("$Q^{2}\\;(GeV^{2})$"),side=1,line=2,cex=0.75)
mtext(TeX("$d\\sigma(Q^{2})/dQ^{2}\\;(fm^{2}/GeV^{2})$"),side=2,line=2,cex=0.75)
mtext(TeX("$Sección\\;eficaz\\;diferencial\\;del\\;proceso\\;\\nu_{\\mu}+n\\to\\mu^{-}+p$"),side=3,las=1,line=1,cex=1)

legend(x="topright",
       legend=c(TeX("$\\epsilon_{i}=1000\\;MeV$"),TeX("$\\epsilon_{i}=500\\;MeV$")),
       col=c("darkorchid","deeppink"),lty=c("solid","solid"),cex=0.6)





