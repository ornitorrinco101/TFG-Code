c23456789
       program Single_Nucleon
       implicit real*8(a-h,p-z)
       common/param/rmN,rmPr,rmNe,pi,GF,hb
       real*8 xq2x(10000),xtnx(10000)
       real*8 enerI(10000), flux(10000), estep(10000)
       character*40 ipivot

!        i,j,k,l,m,n integers

!        open(356,file='q2.in',status='old')
!        jmax = 0
!        do 43 i=1,10000
!        read(356,*,end=45)xtnx(i),xq2x(i)
!        jmax=i
! 43       continue
! 45       continue
        
       jmax=3000
       do i=1,jmax
       xq2x(i) = 0.001d0*i
       enddo

       open(233,file='flujo_interpolado.in')
       read(233,*)ipivot
       nflux = 0
       totalflux = 0.d0
       totalf = 0.d0
       enerI(1) = 0.d0
       flux(1) = 0.d0
       estep(1) = 0.d0
       do 9879 iflux = 2,10000,1
       read(233,*,end=9889)enerI(iflux),flux(iflux)
       estep(iflux) = enerI(iflux)-enerI(iflux-1)
       totalf = flux(iflux)*20.d0*estep(iflux) + totalf
       nflux = iflux
9879       continue
9889       continue
       totalflux = totalf
       print*,'totalflux=',totalflux
       print*,'nflux=',nflux


       iparam = 1

       hb = 0.19733
       pi = datan(1.d0)*4.d0
       GF = 1.16639E-5
       rmPr = 938.27231d0/1000.d0
       rmNe = 939.5739d0/1000.d0
c       rmN = (rmPr + rmNe)/2.d0
       rmN = rmPr

!        xni = 71.890871/1000.
!        xfi = 870./1000.
!        step = 35.945436/1000.

       do 10 j = 100,jmax
       Q2 = xq2x(j)
       tau = Q2/4.d0/rmN/rmN

c ----------------------- FORM FACTORS ----------------------------------
       if(iparam.eq.1)then
       call weakff(Q2,GEPROT,GMPROT,GENEUT,GMNEUT,
     1                  GEVPROT,GMVPROT,GEVNEUT,GMVNEUT,
     1                      gaprot,ganeut)
       EF1 = 0.5D0*(GEVPROT+TAU*GMVPROT)/(1.D0+TAU)
       EF2 = 0.5D0*(GMVPROT-GEVPROT)/(1.D0+TAU)
       endif

       if(iparam.eq.2)then
       call bodek(Q2,GEPROT,GMPROT,GENEUT,GMNEUT,
     1               f1vp,f1vn,f2vp,f2vn)
       ef1 = f1vp
       ef2 = f2vp
       endif

       if(iparam.eq.3.)then
c --  probe: same definition given in paper PRD 82, 092005 , Apendice B -
       call galster(Q2,GEPROT,GMPROT,GENEUT,GMNEUT)
c       call bodek(Q2,GEPROT,GMPROT,GENEUT,GMNEUT,
c     1               f1vp,f1vn,f2vp,f2vn)
       sin2w = 0.23122d0

       ef1EMp = (GEPROT + tau*GMPROT)/(1.d0+tau)
       ef2EMp = (GMPROT - GEPROT)/(1.d0+tau)
       ef1EMn = (GENEUT + tau*GMNEUT)/(1.d0+tau)
       ef2EMn = (GMNEUT - GENEUT)/(1.d0+tau)

       ef1 = (0.5-sin2w)*(ef1EMp - ef1EMn)-sin2w*(ef1EMp+ef1EMn)
     1      -0.5*ef1s
       ef2 = (0.5-sin2w)*(ef2EMp - ef2EMn)-sin2w*(ef2EMp+ef2EMn)
     1      -0.5*ef2s
       endif
c -----------------------------

cccccccccccc ---- Axial ---- ccccccccccccccccccccccccccccccccccccccc
       gA = 1.2695
       xMA = 1.032
       gas = 0.d0
c       FAz = 0.5*gA/(1.+Q2/xMA**2)**2 - 0.5*FAs
       FAz = (0.5*gA - 0.5*gas)/(1.0+tau*3.32)**2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if(iparam.eq.1) then
       FAz = -0.5d0 * gaprot
       endif

c -------------------
       Aq2 = 0.25*(FAz**2*(1.+tau)-(ef1**2-tau*ef2**2)*(1.-tau) 
     1              + 4.*tau*ef1*ef2)
       Bq2 = 0.25*FAz*(ef1 + ef2)
       Cq2 = rmN**2/16./Q2*(FAz**2 + ef1**2 + tau*ef2**2)
c -----------------------------------------------------------------------
       suma = 0.d0
       do 20 iiflux = 1,nflux
       enu = enerI(iiflux)
c       print*,iiflux

       TnMAX = enu*(1.-rmN/(2.*enu+rmN))
       Q2max = 2.*rmN*TnMAX
cccccccccccccc
       tn = Q2/2.d0/rmN
       if(tn.gt.enu)then
       Suma = 0.d0
       goto 20
       endif
cccccccccccccc
c       print*,'Q2max=',Q2max,'   ','enerI=',enu
       if(Q2.ge.Q2max)then
       Suma = 0.d0
       goto 20
       endif

       W = 4.*enu/rmN - Q2/rmN**2
       factor = GF**2*Q2/(2.d0*pi*enu**2)
       xSig = factor * (Aq2 + Bq2*W + Cq2*W**2)
c unidades: Sig ~ GeV⁻⁴

       uni = 1.E13*hb**2
       Sig = xSig*uni
c unidades: Sig ~ 10⁻¹³fm²/GeV² = 10⁻³⁹cm²/GeV²

       Suma = Suma + (Sig*flux(iiflux)*20.d0)*estep(iiflux)
20       continue

       SigFinal = (Suma/totalflux)
       write(25,*)Q2,SigFinal

10       continue

       stop
       end 




C +++++++++++  WEAK FORM FACTORS OF THE NUCLEON ++++++++++++++++++++++++
      SUBROUTINE WEAKFF(Q2,GEPROT,GMPROT,GENEUT,GMNEUT,
     1                  GEVPROT,GMVPROT,GEVNEUT,GMVNEUT,
     1                      gaprot,ganeut)
C +++++++++++++++++++++++++++++++++++++++++
       IMPLICIT REAL*8(A-H,O-Z)
       common/param/rmN,rmPr,rmNe,pi,GF,hb

        open(1001,file='strangeness.in')
! ccccc PARAMETROS de EXTRAÑEZA  cccccccccc
        read(1001,*)
        read(1001,*)smu
       read(1001,*)
        read(1001,*)gas0
        read(1001,*)
        read(1001,*)rhos
        read(1001,*)
        read(1001,*)xMAx
c       print*,smu,rhos,gas0,xMAx
cccccccccccccccccccccccccccccccccccccccc
        close(1001)

       tau = Q2/4./rmN/rmN
!        smu=-0.020
!        gas0 = 0.0
!        rhos = 0.59
!        xMAx = 1230.0

      sinw=0.23122d0
      xiv0=-1.d0
      xivt0=-4.d0*sinw
      xivt1=2.d0-4.d0*sinw
      xia0=1.d0
      xiat0=0.d0
      xiat1=-2.d0

      call lomon(Q2,GEPROT,GMPROT,GENEUT,GMNEUT)


      get1=0.5*(geprot-geneut)
      get0=0.5*(geprot+geneut)
      gmt1=0.5*(gmprot-gmneut)
      gmt0=0.5*(gmprot+gmneut)
      


c factor de forma de extragneza vectorial
      gvd=1.d0/(1.+tau*4.97)**2
      gest=rhos*tau*gvd
      gmst=smu*gvd
c       gest=0.d0
c       gmst=0.d0

      gevprot=xivt1*get1+xivt0*get0+xiv0*gest
      gmvprot=xivt1*gmt1+xivt0*gmt0+xiv0*gmst
      gevneut=-xivt1*get1+xivt0*get0+xiv0*gest
      gmvneut=-xivt1*gmt1+xivt0*gmt0+xiv0*gmst

c factores de forma axial isovector e isoescalar

      dd=1.2695/1.64
      ff=0.64*dd

      xMAx = xMAx/1000.d0
      xLx = 4.d0*rmN**2/xMAx**2

      gad=1.0/(1.0+tau*xLx)**2
      ga3=0.5*gad*(dd+ff)
      ga8=(3.d0*ff-dd)*gad*0.5/dsqrt(3.d0)
      
c factor de forma de extragneza axial
      gas=gas0*gad

      gaprot=xiat1*ga3+xiat0*ga8+xia0*gas
      ganeut=-xiat1*ga3+xiat0*ga8+xia0*gas

      RETURN
      END




C------- GALSTER ----------------
C
       subroutine galster(Q2,GEPROT,GMPROT,GENEUT,GMNEUT)
             IMPLICIT REAL*8 (A-H,O-Z)
       common/param/rmN,rmPr,rmNe,pi,GF,hb
C
       amed = rmN 
       TAU=Q2/4.D0/AMED/AMED
       RMUP=2.79285D0
       RMUN=-1.91304D0
       RGVD=1.D0/((1.D0+4.97D0*TAU)**2.D0)
       RXIN=1.D0/(1.D0+5.6D0*TAU)
       GEPROT=RGVD
       GMPROT=RMUP*RGVD
       GENEUT=-RMUN*TAU*RGVD*RXIN
       GMNEUT=RMUN*RGVD
       return
       end


      subroutine bodek(Q2,GEPROT,GMPROT,GENEUT,GMNEUT,
     1               f1vp,f1vn,f2vp,f2vn)
      implicit real*8(a-h,o-z)
C
C CALCOLO DEI FATTORI DI FORMA F1, F2, GA e FP.
C
C SI CONSIDERA ANCHE IL CONTRIBUTO DELLA STRANEZZA
C
C PARAMETRIZZAZIONE: H.Budd, A.Bodek, J.Arrington, hep-ex/0308005
C
C Q2G in GeV**2
      parameter(a2e=3.253,a4e=1.422,a6e=.08582,a8e=.3318,a10e=-.09371)
      parameter(a12e=.01076)
      parameter(a2m=3.104,a4m=1.428,a6m=.1112,a8m=-.006981)
      parameter(a10m=.0003705,a12m=-0.7063E-05)
      parameter(a2mn=3.043,a4mn=0.8548,a6mn=0.6806,a8mn=-0.1287)
      parameter(a10mn=0.008912)
      parameter(ane=0.942,bne=4.61,umn=-1.913,ump=2.793,p=0.9382723)
      parameter(vm2=0.71)
      parameter(stw2=0.23143,rgaa=1.267)
*
      real*8 gep,gmp,gen,gmn,f1n,f2n,f1p,f2p,q2,tau,gd
      real*8 f1vp,f1vn,f2vp,f2vn,f1s,f2s,gap,gan,fp
*     
c ----- raul -----------------------
      hb=197.3d0

      axm = 1.032d0
      rgas = 0.d0
      rhos = 0.d0
      rmus = 0.d0
c ----------------------------
      axm2=axm*axm
      gda=(1.+q2/axm2)**(-2)
      gd=1./(1.+q2/vm2)/(1.+q2/vm2)
      tau=0.25*q2/p/p
      q4=q2*q2
      q6=q2*q2*q2
      q8=q2*q2*q2*q2
      q10=q2*q2*q2*q2*q2
      q12=q2*q2*q2*q2*q2*q2
      gep=(1.+a2e*q2+a4e*q4+a6e*q6+a8e*q8+a10e*q10+a12e*q12)**(-1)
      gen=-umn*ane*tau*gd/(1.+tau*bne)
      gmp=ump/(a2m*q2+a4m*q4+a6m*q6+a8m*q8+a10m*q10+a12m*q12+1.)
      gmn=umn/(a2mn*q2+a4mn*q4+a6mn*q6+a8mn*q8+a10mn*q10+1.)
      
             GEPROT = GEP
       GENEUT = GEN
       GMPROT = GMP
       GMNEUT = GMN
*
* fattori di forma elettromagnetici
*
      f1p=(gep+tau*gmp)/(1.+tau)
      f2p=(gmp-gep)/(1.+tau)
      f1n=(gen+tau*gmn)/(1.+tau)
      f2n=(gmn-gen)/(1.+tau)
*
* fattori di forma strani
*
      f1s=(rhos+rmus)*tau*gd/(1.+tau)
      f2s=(rmus-tau*rhos)*gd/(1.+tau)             
*
* fattori di forma isovettoriali deboli
*
       f1vp=-(2.*stw2-0.5)*f1p-0.5*f1n-0.5*f1s
       f2vp=-(2.*stw2-0.5)*f2p-0.5*f2n-0.5*f2s
       f1vn=-(2.*stw2-0.5)*f1n-0.5*f1p-0.5*f1s
       f2vn=-(2.*stw2-0.5)*f2n-0.5*f2p-0.5*f2s
*
* fattore di forma assiale
*
       gap=-0.5*(rgaa-rgas)*gda
       gan=0.5*(rgaa+rgas)*gda  
       fp=0.0  
      return
      end


       subroutine lomon(Q2,GEPROT,GMPROT,GENEUT,GMNEUT)  !GKex

       implicit real*8 (a-h,o-z)
C
C----------------------------------------------------------------------------
C                           All masses, etc. in GeV
C
       q2t=q2
              xmnucleon=0.939D0
C                    nucleon mass
              xmrho=0.776D0
C                    rho mass
C
C                           irho=0 : no rho width
C                           irho=1 : with rho width
C
       irho=1 
              xmrho1=xmrho-irho*0.03465D0
C                    mass rho1
              xmrho2=xmrho-irho*0.04374D0
C                    mass rho2
              alpha1=irho*0.0781808D0
C                    parameter alpha1
              alpha2=irho*0.0632907D0
C                    parameter alpha2
              qq1=0.3176D0
C                    parameter q1**2
              qq2=0.1422D0
C                    parameter q2**2
              xmrhop=1.45D0
C                    rho-prime mass
              xmomega=0.784D0
C                    omega mass
              xmomegap=1.419D0
C                    omega-prime mass
              xmphi=1.019D0
C                    phi mass
C
C----------------------------------------------------------------
C
              xkaps=-0.12D0
C                    isoscalar anomalous mag. mom.
              xkapv=3.706D0        
C                    isovector anomalous mag. mom.
              xkaprho=5.51564D0
C                    rho anomalous moment
              xkaprhop=12.0D0
C                    rho-prime anomalous moment
              xkapomega=0.4027D0
C                    omega anomalous moment
              xkapomegap=-2.973D0
C                    omega-prime anomalous moment
              xkapphi=0.01D0
C                    phi anomalous moment
C
C----------------------------------------------------------------
C
              fgrho=0.5596D0
C                    rho coupling
              fgrhop=0.007208D0
C                    rho-prime coupling
              fgomega=0.7021D0
C                    omega coupling
              fgomegap=0.164D0
C                    omega-prime coupling
              fgphi=-0.1711D0
C                    phi coupling
C
C----------------------------------------------------------------
C
              xlam1=0.93088D0
C                    parameter lambda1
              xlam2=2.6115D0
C                    parameter lambda2
              xlamd=1.181D0
C                    parameter lambdaD
              xlamqcd=0.15D0
C                    parameter lambdaQCD
              arat= DLOG(xlamd**2/xlamqcd**2)
              xmuphi=0.2D0
C                    parameter mu-phi
C
C----------------------------------------------------------------------------

              tau=q2t/(4.0D0*xmnucleon**2)
C                    q2 is 4-momentum q**2; tau is the usual
              femrho1=xmrho1**2/(xmrho1**2+q2t)
              femrho2=xmrho2**2/(xmrho2**2+q2t)
              femrhop=xmrhop**2/(xmrhop**2+q2t)
              femomega=xmomega**2/(xmomega**2+q2t)
              femomegap=xmomegap**2/(xmomegap**2+q2t)
              femphi=xmphi**2/(xmphi**2+q2t)
C
C             In the following the notation should be obvious: 
C             e.g., f1somega means F1, isoscalar (s), omega piece, etc.
C
              
              q2tilda=q2t*DLOG((xlamd**2+q2t)/xlamqcd**2)/arat
              f1=xlam1**2/(xlam1**2+q2tilda)
              f2=xlam2**2/(xlam2**2+q2tilda)
              fhad1=f1*f2
              fhad2=f1*f2**2
              fhad1s=fhad1*(q2t/(xlam1**2+q2t))**1.5
              fmuphi=(xmuphi**2+q2t)/xmuphi**2
              fhad2s=fhad2*(fmuphi*xlam1**2/(xlam1**2+q2t))**1.5
              fqcd=xlamd**2/(xlamd**2+q2tilda)
              fhad1qcd=fqcd*f2
              fhad2qcd=fqcd*f2**2

C
              zzf1s=1-fgomega-fgomegap
       zzf2s=xkaps-xkapomega*fgomega-xkapomegap*fgomegap-xkapphi*fgphi
              zzf1v=1-fgrho-fgrhop
              zzf2v=xkapv-xkaprho*fgrho-xkaprhop*fgrhop

c

              f1somega=fgomega*femomega*fhad1
              f1somegap=fgomegap*femomegap*fhad1
              f1sphi= fgphi*femphi*fhad1s
              f1spqcd= zzf1s*fhad1qcd
              f2somega=xkapomega*fgomega*femomega*fhad2
              f2somegap= xkapomegap*fgomegap*femomegap*fhad2
              f2sphi= xkapphi*fgphi*femphi*fhad2s
              f2spqcd= zzf2s*fhad2qcd
              width1=1-alpha1+alpha1/(1+q2t/qq1)**2
              width2=1-alpha2+alpha2/(1+q2t/qq2)
              f1vrho=fgrho*femrho1*fhad1*width1
              f1vrhop= fgrhop*femrhop*fhad1
              f1vpqcd= zzf1v*fhad1qcd
              f2vrho=xkaprho*fgrho*femrho2*fhad2*width2
              f2vrhop= xkaprhop*fgrhop*femrhop*fhad2
              f2vpqcd= zzf2v*fhad2qcd

C
       
              f1s=f1somega+f1somegap+f1sphi+f1spqcd
              f2s=f2somega+f2somegap+f2sphi+f2spqcd
              f1v=f1vrho+f1vrhop+f1vpqcd
              f2v=f2vrho+f2vrhop+f2vpqcd

C
              
              f1prho=0.5D0*(f1vrho)
              f1prhop= 0.5D0*(f1vrhop)
              f1pomega= 0.5D0*(f1somega)
              f1pomegap= 0.5D0*(f1somegap)
              f1pphi= 0.5D0*(f1sphi)
              f1ppqcd= 0.5D0*(f1spqcd+f1vpqcd)
              f1p=f1prho+f1prhop+f1pomega+f1pomegap+f1pphi+f1ppqcd
C
              f1nrho=0.5D0*(-f1vrho)
              f1nrhop= 0.5D0*(-f1vrhop)
              f1nomega= 0.5D0*(f1somega)
              f1nomegap= 0.5D0*(f1somegap)
              f1nphi= 0.5D0*(f1sphi)
              f1npqcd= 0.5D0*(f1spqcd-f1vpqcd)
              f1n=f1nrho+f1nrhop+f1nomega+f1nomegap+f1nphi+f1npqcd
C
              f2prho=0.5D0*(f2vrho)
              f2prhop= 0.5D0*(f2vrhop)
              f2pomega= 0.5D0*(f2somega)
              f2pomegap= 0.5D0*(f2somegap)
              f2pphi= 0.5D0*(f2sphi)
              f2ppqcd= 0.5D0*(f2spqcd+f2vpqcd)
              f2p=f2prho+f2prhop+f2pomega+f2pomegap+f2pphi+f2ppqcd
C
              f2nrho=0.5D0*(-f2vrho)
              f2nrhop= 0.5D0*(-f2vrhop)
              f2nomega= 0.5D0*(f2somega)
              f2nomegap= 0.5D0*(f2somegap)
              f2nphi= 0.5D0*(f2sphi)
              f2npqcd= 0.5D0*(f2spqcd-f2vpqcd)
              f2n=f2nrho+f2nrhop+f2nomega+f2nomegap+f2nphi+f2npqcd
C
              
              gevrho=f1vrho-tau*f2vrho
              gevrhop=f1vrhop-tau*f2vrhop
              gesomega=f1somega-tau*f2somega
              gesomegap=f1somegap-tau*f2somegap
              gesphi=f1sphi-tau*f2sphi
              gespqcd=f1spqcd-tau*f2spqcd
              gevpqcd=f1vpqcd-tau*f2vpqcd
              ges=gesomega+gesomegap+gesphi+gespqcd
              gev=gevrho+gevrhop+gevpqcd
C
              gmvrho=f1vrho+f2vrho
              gmvrhop=f1vrhop+f2vrhop
              gmsomega=f1somega+f2somega
              gmsomegap=f1somegap+f2somegap
              gmsphi=f1sphi+f2sphi
              gmspqcd=f1spqcd+f2spqcd
              gmvpqcd=f1vpqcd+f2vpqcd
              gms=gmsomega+gmsomegap+gmsphi+gmspqcd
              gmv=gmvrho+gmvrhop+gmvpqcd
C
              geprho= 0.5D0*(gevrho)
              geprhop= 0.5D0*(gevrhop)
              gepomega= 0.5D0*(gesomega)
              gepomegap= 0.5D0*(gesomegap)
              gepphi= 0.5D0*(gesphi)
              geppqcd= 0.5D0*(gespqcd+gevpqcd)
              xgep=geprho+geprhop+gepomega+gepomegap+gepphi+geppqcd
C
              genrho= 0.5D0*(-gevrho)
              genrhop= 0.5D0*(-gevrhop)
              genomega= 0.5D0*(gesomega)
              genomegap= 0.5D0*(gesomegap)
              genphi= 0.5D0*(gesphi)
              genpqcd= 0.5D0*(gespqcd-gevpqcd)
              xgen=genrho+genrhop+genomega+genomegap+genphi+genpqcd
C
              gmprho= 0.5D0*(gmvrho)
              gmprhop= 0.5D0*(gmvrhop)
              gmpomega= 0.5D0*(gmsomega)
              gmpomegap= 0.5D0*(gmsomegap)
              gmpphi= 0.5D0*(gmsphi)
              gmppqcd= 0.5D0*(gmspqcd+gmvpqcd)
              xgmp=gmprho+gmprhop+gmpomega+gmpomegap+gmpphi+gmppqcd
C
              gmnrho= 0.5D0*(-gmvrho)
              gmnrhop= 0.5D0*(-gmvrhop)
              gmnomega= 0.5D0*(gmsomega)
              gmnomegap= 0.5D0*(gmsomegap)
              gmnphi= 0.5D0*(gmsphi)
              gmnpqcd= 0.5D0*(gmspqcd-gmvpqcd)
              xgmn=gmnrho+gmnrhop+gmnomega+gmnomegap+gmnphi+gmnpqcd

C
       GEPROT=xgep
       GENEUT=xgen
       GMPROT=xgmp
       GMNEUT=xgmn
              return
              end
