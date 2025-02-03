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

  int process; // 0 EM; 1 CC
  in >> process;
  in.ignore(1000, '\n');  // line skipper
    
  int A, Z, pF; // A: mass number; Z: atomic number; Fermi momentum (MeV)
  in >> A >> Z >> pF;
  in.ignore(1000, '\n');  // line skipper

  double thetaf_deg; // Scattering angle (degrees)
  in >> thetaf_deg;
  in.ignore(1000, '\n');
  
    double thetaf = thetaf_deg*Pi/180.,    costheta_f = cos(thetaf),    sintheta_f=sin(thetaf);

  int Helicity; // 1 antineutrinos; -1 neutrinos
  in >> Helicity; 
  in.ignore(1000, '\n');

  in.close();
  // // // // // // // // // 

  double enu_flux[201]; double flux[201]; int iflux_max;
  readtheflux( enu_flux, flux, iflux_max );
  
//   // // // building the 4-vectors of incoming lepton and defining the reference frame // // //   
//   double ki=Ei; //3-mom of initial lepton (ultrarrelativistic approximation)
//   double Ki[4] = {Ei,0,0,ki}; //initial lepton 
//   double Ki_c[4] = {Ei,0,0,-ki};
//   // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //   
  
  // // // OUTPUT FILE // // //
  int itheta=thetaf_deg;
  ofstream out; //name of the output unit in the code 
  out.open("th"+to_string(itheta)+".out"); //name of the output file
  out.precision(4); //number of digits in the printed numbers 
  // // // // // // // // // //
  
  for(double kf=1; kf<=2000.; kf+=25)
  {

    cout << "kf= " << kf << endl;
    
    
    // // final lepton // //
    double Ef = sqrt(kf*kf+muonmass*muonmass);       //if(Ef >= Ei || Ef < 0){continue;}
    double kf2 = pow(Ef,2) - pow(muonmass,2);                 if(kf2 < 0){continue;}  // "continue" breaks the flow and sends it back to the next value of the loop (w in our case)  
//     double kf=sqrt(kf2);
    double Kf[4] = {Ef,kf*sintheta_f,0,kf*costheta_f};
    double Kf_c[4] = {Kf[0],-Kf[1],-Kf[2],-Kf[3]};
    // // // // // // // // 

    double d3sig_ave=0; //declaration and initialization of the cross section at zero

    double Eistep=10; //MeV
    for(double Ei=enu_flux[0]; Ei<=enu_flux[200]; Ei+=Eistep)
    {

      // // // building the 4-vectors of incoming lepton and defining the reference frame // // //   
      double ki=Ei; //3-mom of initial lepton (ultrarrelativistic approximation)
      double Ki[4] = {Ei,0,0,ki}; //initial lepton 
      double Ki_c[4] = {Ei,0,0,-ki};
      // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //   
      
      if(Ei <= muonmass){continue;} //cross section is zero if Ei<muonmass

      if(Ef > Ei){continue;}
      
      double w = Ei - Ef; 
      
      double flujo=0;
      int stop=0;
      findtheflux(Ei, enu_flux, flux, iflux_max, flujo, stop);  if(stop != 0){continue;}
      
      // // exchanged boson // //
      double Q[4];
      for(int i=0; i<4; i++)
      {
        Q[i] = Ki[i] - Kf[i];
      }
      if(Q[0] < 0){continue;} //the energy must be positive
      
      double Q2 = Q[0]*Q[0] - Q[1]*Q[1] - Q[2]*Q[2] -Q[3]*Q[3];   if(Q2 >= 0){continue;} //Q2 must be negative
      
      double q2 = w*w - Q2;  //magnitude of the 3-momentum transfer (squared)  
      if( q2 <= 0 ){ continue; }
      double q = sqrt(q2);
      // // // // // // // // 



      
      
      double d3sig=0; //declaration and initialization of the cross section at zero
      
      double thetaN_deg_step=1.;
      for(double thetaN_deg=0; thetaN_deg<=180; thetaN_deg+=thetaN_deg_step)
      {
        double thetaN=thetaN_deg*Pi/180;  
          if(thetaN_deg==0){thetaN=0.00001;}
          if(thetaN_deg==180){thetaN=Pi-0.00001;}
        
        double thetaN_step=thetaN_deg_step*Pi/180;
        double costhetaN = cos(thetaN),  sinthetaN = sin(thetaN);
        
        
        double phiN_deg_step=1.;
        double d4sig=0;
        for(double phiN_deg=0; phiN_deg<=360; phiN_deg+=phiN_deg_step)
        {
          double phiN=phiN_deg*Pi/180;
          double phiN_step=phiN_deg_step*Pi/180;
          double cosphiN = cos(phiN),  sinphiN = sin(phiN);

          
          // // // // // // // // //
          // // final nucleon // //

          // // direction of the final nucleon // //
          double hpNx = sinthetaN*cosphiN, 
                hpNy = sinthetaN*sinphiN, 
                hpNz = costhetaN;
          // // // // // // // // // // // // // // 
          
          double costheta_qN = (Q[1]*hpNx + Q[2]*hpNy + Q[3]*hpNz)/q;
          
          double AA = Q2/(2*w),   BB = q*costheta_qN/w; 
          double deter = AA*AA + (BB*BB - 1)*MN2;               if(deter<0){continue;}
          
          double pN = ( -AA*BB + sqrt(deter) ) / ( BB*BB - 1 );    if(pN<0){continue;}
          
          // // Pauli blocking // //
          if(pN < pF){continue;}
          // // // // // // // // //
          
          double PN[4]={ sqrt(pN*pN + MN2), pN*hpNx, pN*hpNy, pN*hpNz };
          if(PN[0] < MN){continue;} //the energy must be larger than the nucleon mass
          // // // // // // // // // // 
          // // // // // // // // // //
          
          
          // // initial nucleon // //
          double P[4];
          for(int i=0; i<4; i++)
          {
            P[i] = PN[i] - Q[i];
          }
          double p = sqrt(pow(P[1],2) + pow(P[2],2) + pow(P[3],2));

          if( p > pF ){continue;} // above the Fermi level, then not a bound nucleon in the RFG model 
          // // // // // // // // // //

          // // Lepton tensor // // 
          complex<double> Lmn[4][4]; //~MeV^2
          Lmunu(process, Helicity, Ki_c, Kf_c, Lmn); //output is Lmn (with indices down), the rest is input 
          // // // // // // // // // // //

          // // Contraction // //
          complex<double> LmnHmn=0;
          
          if(process==0){
            // // Hadron tensor // // 
            complex<double> Hmn_p[4][4]={}; //~dimensionless
            Hmunu(process, 1, P, PN, Hmn_p); //output is Hmn, the rest is input
            complex<double> Hmn_n[4][4]={}; //~dimensionless
            Hmunu(process, 0, P, PN, Hmn_n); //output is Hmn, the rest is input
            // // // // // // // // // // //
            
            for(int i=0; i<4; i++)
            {
              for(int j=0; j<4; j++)
              {
                
                LmnHmn += Lmn[i][j] * ( double(Z)*Hmn_p[i][j] + double(A-Z)*Hmn_n[i][j] ); //~MeV^2
                
              }
            }  
          }
          else if(process==1)
          { 
            // // Hadron tensor // // 
            complex<double> Hmn[4][4]={}; //~dimensionless
            Hmunu(process, 0, P, PN, Hmn); //output is Lmn, the rest is input
            // // // // // // // // // // //
            
            double NN;
            if(Helicity==-1){NN=double(A-Z);} //neutrinos interact with neutrons
            else if(Helicity==1){NN=double(Z);} //antineutrinos interact with protons
            
            for(int i=0; i<4; i++)
            {
              for(int j=0; j<4; j++)
              {
                
                LmnHmn += Lmn[i][j] * NN *Hmn[i][j] ; //~MeV^2
                
              }
            }            
          }
          // // // // // // // // // 
          
          double FF; // factor from the exchanged boson (photon or W)
          if(process == 1){ FF = pow(G_Fermi*Cabibbo/sqrt(2.),2); } //Weak Charged-Current interaction g^4/(2^6 M_W^4) * (cosCabibbo)^2, ~MeV^-4
          else if(process == 0){ FF = pow(4*Pi*Alpha/Q2,2); } //EM interaction e^4/Q^4, ~MeV^-4
          
          // // computing the cross section // //
          double Ki_P = Ki[0]*P[0] - Ki[1]*P[1] - Ki[2]*P[2] - Ki[3]*P[3];
          double KK = MN2/(abs(Ki_P) * Kf[0]*PN[0]) * FF; //~MeV^-6
        
          double frec = abs(pN/PN[0] - (pN-q*costheta_qN)/P[0]); //~dimensionless

          // d5sig/(dkf dOmegaf dOmegaN)
          double d5sig = 3/(4*Pi*pow(pF,3)) * KK * pow(pN*kf/(2*Pi),2) / frec * abs(LmnHmn); // ~MeV^-3 
          
          // // // // // // // // // // // // 
          
          // // integral over phiN // // 
          if(phiN_deg==0 || phiN_deg==360){phiN_step=phiN_step/2;}
          d4sig += d5sig * phiN_step;
          // // // // // // // // // // 

        } //phiN
        
        // // integral over phiN // //
        if(thetaN_deg==0 || thetaN_deg==180){thetaN_step=thetaN_step/2;}
        d3sig += d4sig * sinthetaN*thetaN_step;
        // // // // // // // // // // 
        
      }//thetaN
    
      d3sig_ave += d3sig*Eistep*flujo;
    
    } //Ei

    double d2sig=0;
    
    // // printing the output in a file (and in the screen) // //
    if(d3sig_ave>0)
    {
      // // so far  d3sig is in 1/(MeV^3 srad)

      d2sig = 2*Pi*d3sig_ave*pow(hbc,2);// 2Pi to integrate over phil; hbc^2 to get the cross section in fm^2/(MeV srad)
      
      d2sig = d2sig *(1.E16)/A; // 1.E-26 (fm^2 --> cm^2) * 1.E3 (1/MeV --> 1/GeV) /A (per nucleon) * 1.E39 (to have it in 10^-39 cm^2)
      
      double kfprime = sqrt( pow(Kf[0]-20.,2) - pow(muonmass,2) );
      out << kf << " " << kfprime << " " << d2sig << endl; //d3s/(dkf dOmegaf) ; units~(fm^2/MeV srad)
    }
    
    cout << kf << " " << d2sig << endl;
    // // // // // // // // // // // // // // 
    
  } //kf
  
  return 0; 

} //main closed



