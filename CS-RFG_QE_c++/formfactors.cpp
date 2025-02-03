#include "formfactors.h"

void formfactors( double Qsq, double &F1p, double &F2p, double &F1n, double &F2n, double &GA, double &FP)
{
  
  double QsqGeV, xmu, tau;
  QsqGeV = Qsq/1.E6;
  xmu = 2*MN;
  tau = Qsq/(xmu*xmu);

  double DipV, DipA, MA2, MV2;
  MA2 = 1.065024; //MA=1.032 GeV
  MV2 = 0.710649; //MV=0.843 GeV
  DipV = 1./pow( 1+QsqGeV/MV2, 2 ); 
  DipA = 1./pow( 1+QsqGeV/MA2, 2 ); 

  double rmup=2.79285;
  double rmun=-1.91304;

  double GEp, GEn, GMp, GMn;
  
  int iff;
//     iff = 0; // Dipole
  iff = 1; // Kelly
  
  if(iff == 0){ //Dipole parametrization

    GEp = DipV;
    GMp = rmup*GEp;
    GEn = -rmun*tau*GEp/(1.+5.6*tau);
    GMn = rmun*GEp;

  }else if(iff == 1){ //parametrization by Kelly Phys. Rev. C 70, 068202 (2004)

    double As=1.70;
    double Bs=3.30;

    double a0=1.0;

    double a1Gep=-0.24;
    double b1Gep=10.98;
    double b2Gep=12.82;
    double b3Gep=21.97;

    double a1Gmp=0.12;
    double b1Gmp=10.97;
    double b2Gmp=18.86;
    double b3Gmp=6.55;
    
    double a1Gmn=2.33;
    double b1Gmn=14.72;
    double b2Gmn=24.20;
    double b3Gmn=84.1;

    GEn = As*tau/(1.+Bs*tau)*DipV;
    GEp = (a0 + a1Gep*tau)/(1. + b1Gep*tau + b2Gep*pow(tau,2) + b3Gep*pow(tau,3));
    GMp = rmup*(a0 + a1Gmp*tau)/(1. + b1Gmp*tau + b2Gmp*pow(tau,2) + b3Gmp*pow(tau,3));
    GMn = rmun*(a0 + a1Gmn*tau)/(1. + b1Gmn*tau + b2Gmn*pow(tau,2) + b3Gmn*pow(tau,3));

  }
    
  F1p = (GEp+tau*GMp)/(1.+tau);
  F2p = (GMp-GEp)/(1.+tau);
  F1n = (GEn+tau*GMn)/(1.+tau);
  F2n = (GMn-GEn)/(1.+tau);

  GA = -gA * DipA ;
  FP = GA*2*MN/(Qsq + Mpi2);

}

