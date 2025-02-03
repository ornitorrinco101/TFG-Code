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
    
  int nucleon; // 0 neutron; 1 proton 
  in >> nucleon;
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
  out.open("d2sigma.out"); //name of the output file
  out.precision(4); //number of digits in the printed numbers 
  // // // // // // // // // //
  
  for(double thetaf_deg=34.4115; thetaf_deg<=34.4115; thetaf_deg +=2)
  {

//     cout << "thetaf_deg= " << thetaf_deg << endl;
    double thetaf = thetaf_deg*Pi/180.,    costheta_f = cos(thetaf),    sintheta_f=sin(thetaf);
    
    
    double Eistep=10; //MeV
    for(double Ei=enu_flux[0]; Ei<=enu_flux[200]; Ei+=Eistep)
    {
      if(Ei<mf){continue;} //if the initial energy is lower than the mass of the final lepton, then we're out if physical p.s. 
      
      double ki=Ei; //3-mom of initial lepton (ultrarrelativistic approximation)
      double Ki[4] = {Ei,0,0,ki}; //initial lepton 
      double Ki_c[4] = {Ei,0,0,-ki};

      // // final lepton // // 
      // // in the ultrarrelativistic limit: 
//       double Ef_1 = Ei / (1 + Ei/MN*(1 - costheta_f));
      // // // // // // 
            
      double flujo=0;
      int stop=0;
      findtheflux(Ei, enu_flux, flux, iflux_max, flujo, stop);  if(stop != 0){continue;}
//       cout << Ei << " " << flujo << endl;
      
      // // considering the lepton masses:
      double AA, BB, CC, determ, mi=0;
      AA = Ei*MN + (mi*mi + mf*mf)/2;   BB = Ei + MN;  CC = ki*costheta_f;    
      determ = pow(AA*BB,2) - (AA*AA + CC*CC*mf*mf)*(BB*BB-CC*CC);   if(determ<0){continue;}
      double Ef;
      if(costheta_f<0)
      {
        Ef = (AA*BB - sqrt(determ))/(BB*BB-CC*CC);
      }else
      {
        Ef = (AA*BB + sqrt(determ))/(BB*BB-CC*CC);
      }
//       cout << Ef << endl;
      // // // // // // // // //     
            
      if(Ef >= Ei || Ef < mf){continue;}
      double kf = sqrt(Ef*Ef - mf*mf);
      double Kf[4] = {Ef,kf*sintheta_f,0,kf*costheta_f};
      double Kf_c[4] = {Kf[0],-Kf[1],-Kf[2],-Kf[3]};
      // // // // // // // // 

      // // exchanged boson // //
      double Q[4];
      for(int i=0; i<4; i++)
      {
        Q[i] = Ki[i] - Kf[i];
      }
      if(Q[0] < 0){continue;} //the energy must be positive
      double Q2 = Q[0]*Q[0] - Q[1]*Q[1] - Q[2]*Q[2] -Q[3]*Q[3];   if(Q2 >= 0){continue;} //Q2 must be negative
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
        double frec = abs(-1 + (ki-kf*costheta_f)*Ki[0]/(ki*PN[0]) ); // expression valid in general  //~dimensionless

        // // // // // // // //     
      
      // d3sigma/(dcosthetaf dphif dkf)
      double d3sig =  kf*kf/frec * flujo * KK /pow(2*Pi,2) * abs(LmnHmn); // ~MeV^2 * MeV^-1 * MeV^-6 * MeV^2 = MeV^-3
      
      double d2sig = 2*Pi*pow(hbc,2) * d3sig; // 2pi to integrate over phi_f, hbc^2 to get the cross section in fm^2, 
      d2sig = d2sig * 1.E16; //10^13 to go from fm^2/MeV to 10^-39 cm^2/GeV,
      
      // // printing the output in a file (and in the screen) // //
      out << thetaf << " " << Ei << " " << Kf[0] << " " << kf << " " << d2sig << endl;
          
      cout << thetaf << " " << Ei << " " << Kf[0] << " " << kf << " " << d2sig << endl;
    
    } //Ei

    // // // // // // // // // // // // // // 
    
  } //thetaf_deg
  
  return 0; 

} //main closed



