library(latex2exp)

flux=read.table("/Users/mcytron/TFG-CODE/RFG_QE_c++/con-flujo/t2kflux_numu_2016.dat",header=FALSE)
names(flux)=c("enu","flux")


plot(flux$enu/1000,1000*flux$flux,type="l",lty="dotted",
     xlab=TeX("$\\epsilon_{i}\\;(GeV)$"),ylab="",cex.lab=1)
title(ylab=TeX("$\\phi(\\epsilon_{i})\\;(GeV^{-1})$"),line=2,cex.lab=1)


