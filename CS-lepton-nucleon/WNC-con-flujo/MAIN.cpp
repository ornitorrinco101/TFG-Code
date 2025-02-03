#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <cstdarg>

using namespace std;

#include "Constants.h"
#include "Leptonic.h"
#include "Hadronic.h"
#include "flux.h"

int main()
{

  // An input file is opened   
  ifstream in("Input.txt"); // Input file 

  int process;
  in >> process; // 0 EM; 1 CC; 2 WNC
  in.ignore(1000, '\n');  // line skipper
    
  int Helicity;  
  in >> Helicity; // 1 antineutrinos; -1 neutrinos
  in.ignore(1000, '\n');

  in.close();
  // // // // // // // // // 

  double enu_flux[201]; double flux[201]; int iflux_max;
  readtheflux( enu_flux, flux, iflux_max );

  double mf = 0.;
  if(process == 1){ mf = muonmass;}
  
  
  // // // building the 4-vectors of incoming particles and defining the reference frame // // // 
  double P[4] = {MN,0,0,0}; //target nucleon  
  // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //   
  
  // // // OUTPUT FILE // // //
  ofstream out; //name of the output unit in the code 
  out.open("dsigma-dQ2.out"); //name of the output file
  out.precision(4); //number of digits in the printed numbers 
  // // // // // // // // // //
  
  for(double Q2GeV = 0.001; Q2GeV <= 2.; Q2GeV += 0.01)
  {
    double Q2 = -Q2GeV*1.E6; //MeV^2
//     cout << "Q2GeV= " << Q2GeV << endl;
    
    double TN=-Q2/(2*MN);
    double EN = TN + MN;
    
    double dsig[2]={};
    for(int nucleon=0; nucleon<=1; nucleon++){

      double Efstep=1.; //MeV
      for(double Ef=mf; Ef<=3000.; Ef+=Efstep)
      {
        double Ei = Ef - MN + EN;   if(Ei<Ef){continue;}
        double kf = sqrt(Ef*Ef - mf*mf);      

        double TN_max = Ei*(1.-MN/(2.*Ei+MN)); if(TN>TN_max){continue;}
        
        double flujo=0;
        int stop=0;
        findtheflux(Ei, enu_flux, flux, iflux_max, flujo, stop);  if(stop != 0){continue;}
  //       cout << Ei << " " << flujo << endl;

        double ki=Ei; //3-mom of initial lepton (ultrarrelativistic approximation)
        double Ki[4] = {Ei,0,0,ki}; //initial lepton 
        double Ki_c[4] = {Ei,0,0,-ki};

        double w = Ei - Ef;
        double q2 = w*w - Q2;   if(q2<w*w){continue;}
        double q = sqrt(q2);
        
        double costheta_f = (ki*ki + kf*kf - q2)/(2*ki*kf);  if( abs(costheta_f)>1 ){continue;}    
        double sintheta_f = sqrt(1-costheta_f*costheta_f);   
        
        double Kf[4] = {Ef, kf*sintheta_f, 0, kf*costheta_f}; //4-mom of final lepton
        double Kf_c[4] = {Kf[0], -Kf[1], -Kf[2], -Kf[3]};
        // // // // // // // // 
        

        // // exchanged boson // //
        double Q[4];
        for(int i=0; i<4; i++)
        {
          Q[i] = Ki[i] - Kf[i];
        }
        // // // // // // // // 
        
        // // final nucleon // //
        double PN[4];
        for(int i=0; i<4; i++)
        {
          PN[i] = P[i] + Q[i];
        }
        if(PN[0] < MN){continue;} //the energy must be larger than the nucleon mass
        
    //     cout << "sqrt(PN^2) = " << sqrt(PN[0]*PN[0] - PN[1]*PN[1] - PN[2]*PN[2] - PN[3]*PN[3]) << ", MN=" << MN << endl; //checking that everything is OK 
        // // // // // // // // 

        // // Lepton tensor // // 
        complex<double> Lmn[4][4]={}; //~MeV^2
        Lmunu(process, Helicity, Ki_c, Kf_c, Lmn); //output is Lmn (with indices down), the rest is input 
        // // // // // // // // // // //

        // // Hadron tensor // // 
        complex<double> Hmn[4][4]={}; //~dimensionless
        Hmunu(process, nucleon, P, PN, Hmn); //output is Lmn, the rest is input
        // // // // // // // // // // //

        // // Contraction // //
        complex<double> LmnHmn=0;
        for(int i=0; i<4; i++)
        {
          for(int j=0; j<4; j++)
          {
            
            LmnHmn += Lmn[i][j]*Hmn[i][j]; //~MeV^2
            
          }
        }    
        // // // // // // // // // 
        
        double FF; // factor from the exchanged boson (photon or W)
        if(process == 2){ FF = pow(G_Fermi/sqrt(2.),2); } //WNC interaction
        else if(process == 1){ FF = pow(G_Fermi*Cabibbo/sqrt(2.),2); } //Weak Charged-Current interaction g^4/(2^6 M_W^4) * (cosCabibbo)^2, ~MeV^-4
        else if(process == 0){ FF = pow(4*Pi*Alpha/Q2,2); } //EM interaction e^4/Q^4, ~MeV^-4
        
        
        // // computing the cross section // //
        double Ki_P = Ki[0]*P[0] - Ki[1]*P[1] - Ki[2]*P[2] - Ki[3]*P[3];
        double KK = MN2/(abs(Ki_P) * Kf[0]*PN[0]) * FF; //~MeV^-6
        
          // // recoil factor
        double frec = 1 + Kf[0]*(kf-ki*costheta_f)/(kf*PN[0]); // expression valid in general  //~dimensionless

          // // // // // // // //     
        
        // d3sigma/(dcosthetaf dphif dEf) ~ MeV^-3
        double d3sig =  Ef*kf/frec * flujo * KK /pow(2*Pi,2) * abs(LmnHmn); //~MeV^2 * MeV^-1 * MeV^-6 * MeV^2 = MeV^-3
        
        // d2sigma/(dcosthetaf dEf) ~ fm^2/MeV
        double d2sig = 2*Pi*pow(hbc,2) * d3sig; // 2pi to integrate over phi_f, hbc^2 to get the cross section in fm^2/MeV

        // // Jacobian dQ^2/dcosthetaf.  
        // // Using Q2 = w^2 - q^2 = 2 Ei Ef (costheta_f-1), one gets dQ^2/dcosthetaf = 2Ei ( dEf/dcosthetaf (costheta_f-1) + Ef).  On the other hand, we have that  Ef = Ei / (1 + Ei/MN*(1 - costheta_f)), so dEf/dcosthetaf = ... = -Ef^2/MN
        double jacQ2 = 2*Ei*Ef*(1. + Ef/MN*(costheta_f-1));

        // // d2sigma/(dQ2 dEf) ~ fm^2/MeV^3
        double d2sigQ2 = d2sig/abs(jacQ2); 
        // // // // // // // // // // // // 
        
  //       // // printing the output in a file (and in the screen) // //
  //       out << thetaf << " " << Ei << " " << Kf[0] << " " << kf << " " << d2sig << endl;
  //           
  //       cout << thetaf << " " << Ei << " " << Kf[0] << " " << kf << " " << d2sig << endl;
      
        // // d2sigma/dQ2 ~ fm^2/MeV^2      
        dsig[nucleon] += d2sigQ2 * Efstep; //integrating over dEf
      } //Ef
      
      dsig[nucleon] = dsig[nucleon]*1.E19 ; // 1.E16 to go from fm^2/MeV^2 to 10^-39 cm^2/GeV^2

    }//nucleon


    // // printing the output in a file (and in the screen) // //
    out << -Q2*1.E-6 << " " << dsig[0] /*neutron*/ << " " << dsig[1] /*proton*/ << endl;
        
    cout << "Q2= " << -Q2*1.E-6 << ", dsig/dQ2(neutron)= " << dsig[0] << ", dsig/dQ2(proton)= " << dsig[1] << endl;

    // // // // // // // // // // // // // // 
    
  } //Q2
  
  return 0; 

} //main closed



