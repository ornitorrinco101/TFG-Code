#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <cstdarg>

using namespace std;

#include "Constants.h"
#include "Leptonic.h"
#include "Hadronic.h"

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

  double Ei; // incoming lepton energy (MeV)
  in >> Ei;   
  in.ignore(1000, '\n');

  int Helicity;  
  in >> Helicity; // 1 antineutrinos; -1 neutrinos
  in.ignore(1000, '\n');

  in.close();
  // // // // // // // // // 

  double mf = electronmass;
  if(process == 1){ mf = muonmass;}
  if(process == 2){ mf = 0;}
  
  // // // building the 4-vectors of incoming particles and defining the reference frame // // // 
  double P[4] = {MN,0,0,0}; //target nucleon
  
  double ki=Ei; //3-mom of initial lepton (ultrarrelativistic approximation)
  double Ki[4] = {Ei,0,0,ki}; //initial lepton 
  double Ki_c[4] = {Ei,0,0,-ki}; //initial lepton covariant contravariant?
  // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //   
  
  // // // OUTPUT FILE // // //
  int iEi=Ei;
  ofstream out; //name of the output unit in the code 
  out.open("Ei"+to_string(iEi)+".out"); //name of the output file 
  out.precision(4); //number of digits in the printed numbers 
  // // // // // // // // // //
  
  for(double thetaf_deg=1; thetaf_deg<=180; thetaf_deg +=2)
  {

//     cout << "thetaf_deg= " << thetaf_deg << endl;
    double thetaf = thetaf_deg*Pi/180.,    costheta_f = cos(thetaf),    sintheta_f=sin(thetaf);
    
    double d2sig=0; //declaration and initialization of the cross section at zero
    
    // // final lepton // // 
// //       // // in the ultrarrelativistic limit: 
//       double Ef_1 = Ei / (1 + Ei/MN*(1 - costheta_f));
// //       // // // // // // 
      
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
      // // // // // // // // //     
            
      if(Ef >= Ei || Ef < 0){continue;}
      double kf = sqrt(Ef*Ef-mf*mf);
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
//       double frec_1 = 1 + (Kf[0]-Ki[0]*costheta_f)/PN[0]; // expression valid in the ultrarrelativistic limit for the leptons //~dimensionless
      double frec_2 = 1 + Kf[0]*(kf-ki*costheta_f)/(kf*PN[0]); // expression valid in general  //~dimensionless

      double frec = frec_2;
//       cout << "frec_1= " << frec_1 << ", frec_2= " << frec_2 << endl;
      // // // // // // // //     
    
    // d2sigma/dcosthetaf
    d2sig = KK * pow(kf/(2*Pi),2)/abs(frec) * abs(LmnHmn); // ~MeV^-6 * MeV^2 * MeV^2 = MeV^-2 
    
    double dsig = 2*Pi*pow(hbc,2)*d2sig; // 2pi to integrate over phi_f, hbc^2 to get the cross section in fm^2, 
//     dsig = dsig * 1.E13; //10^13 to go from fm^2 to 10^-39 cm^2,

    // // Jacobian dQ^2/dcosthetaf.  Using Q2 = w^2 - q^2 = 2 Ei Ef (costheta_f-1), one gets dQ^2/dcosthetaf = 2Ei ( dEf/dcosthetaf (costheta_f-1) + Ef).  On the other hand, we have that  Ef = Ei / (1 + Ei/MN*(1 - costheta_f)), so dEf/dcosthetaf = ... = -Ef^2/MN
    double jacQ2 = 2*Ei*Ef*(1. + Ef/MN*(costheta_f-1));
    double dsigQ2 = dsig/abs(jacQ2)  * 1.E6; // 10^6 to go from MeV^-2 to GeV^-2
    // // // // // // // // // // // // 
    

    // // printing the output in a file (and in the screen) // //
    out << Q[0] << " " << Kf[0] << " " << PN[0] << " " << -Q2/1.E6 /*4*/ << " " << dsigQ2 << " " << costheta_f /*6*/ <<" "<< dsig << endl;
        
    cout << "thetaf_deg= " << thetaf_deg << ", w= " << Q[0] << ", Ef= " << Kf[0] << ", EN= " << PN[0] << ", |Q2|= " << -Q2/1.E6 << ", dsig/dQ2= " << dsigQ2 << ", dsig/dcosth_f= " << dsig << endl;

    // // // // // // // // // // // // // // 
    
  } //thetaf_deg
  
  return 0; 

} //main closed



