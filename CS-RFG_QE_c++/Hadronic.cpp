#include "Hadronic.h"

void Hmunu( int process, int nucleon, double Pm[], double PN[], complex<double> Hadron[][4])
{
  
  double Q[4]={ PN[0]-Pm[0], PN[1]-Pm[1], PN[2]-Pm[2], PN[3]-Pm[3] };
  double Qsq = -(pow(Q[0],2) - pow(Q[1],2) - pow(Q[2],2) - pow(Q[3],2)); // defined positive here Qsq > 0 

  double F1p, F1n, F2p, F2n, GA, FP;  
  formfactors( Qsq, F1p, F2p, F1n, F2n, GA, FP );

  if(process==0)
  {
    double PmPN = Pm[0]*PN[0] - Pm[1]*PN[1] - Pm[2]*PN[2] - Pm[3]*PN[3];
    
    double F1, F2;
    if(nucleon == 1){F1 = F1p; F2=F2p;}
    else{F1 = F1n; F2=F2n;}
    
    for(int i=0; i<4; i++){ 
      for(int j=0; j<4; j++){
    
        Hadron[i][j] = 1/(2*MN2) * 
        ( pow((F1+F2),2) * (Pm[i]*PN[j] + Pm[j]*PN[i] + (MN2-PmPN)*gmunu[i][j])
        +( pow(F2/(2*MN),2) *( PmPN + MN2 ) - F2*(F1+F2) ) * (Pm[i]+PN[i])*(Pm[j]+PN[j]) );
        
      }
    }
        
  }
  else if(process==1){
    
    double Pm_c[4] = {Pm[0], -Pm[1], -Pm[2], -Pm[3]};
    double Q_c[4] = {Q[0], -Q[1], -Q[2], -Q[3]};
    
    double F1V = (F1p - F1n);
    double F2V = (F2p - F2n);    

    double GA2 = GA*GA;
    double FP2 = FP*FP; //MeV^-2
    double tau = Qsq/(4.*MN2);
    
    double W1 = tau * ( pow(F1V+F2V,2) + GA2 ) + GA2;
    double W2 = pow(F1V,2) + tau * pow(F2V,2) + GA2;
    double W3 = 2.*GA*(F1V+F2V);
    double W4 = pow(F2V,2)/4. * (tau-1.) - F1V*F2V/2. - FP*GA*MN + MN2*tau*FP2;
    double W5 = W2;

    
// // // // // Symmetric tensor // // // // 
    double Hadron_s[4][4]={};
    for(int i=0; i<4; i++){ 
      for(int j=0; j<4; j++){
        Hadron_s[i][j] = -W1*MN2 * gmunu[i][j] 
                    + W2 * Pm[i]*Pm[j] 
                    + W4 * Q[i]*Q[j] 
                    + W5/(2.) * ( Pm[i]*Q[j] + Q[i]*Pm[j] );        
      }
    }
// // // // // // // // // // // // // // //   
    
// // // // // Antisymmetric tensor // // // // 
    double Hadron_a[4][4]={};
//     for(int i=0; i<4; i++){ 
//       for(int j=0; j<4; j++){
//         for(int al=0; al<4; al++){ 
//           for(int be=0; be<4; be++){
//             Hadron_a[i][j] += epsilon[i][j][al][be] * Pm_c[al]*Q_c[be];
//           }
//         }
//       }
//     }
//        CONVENTION:           epsilon^(0123) = -1
      Hadron_a[0][1] = (-1 * Pm_c[2]*Q_c[3] + 1 * Pm_c[3]*Q_c[2]) * W3/2.;
      Hadron_a[1][0] = -Hadron_a[0][1];
      
      Hadron_a[0][2] = (1 * Pm_c[1]*Q_c[3] + -1 * Pm_c[3]*Q_c[1]) * W3/2.;
      Hadron_a[2][0] = -Hadron_a[0][2];
      
      Hadron_a[0][3] = (-1 * Pm_c[1]*Q_c[2] + 1 * Pm_c[2]*Q_c[1]) * W3/2.;
      Hadron_a[3][0] = -Hadron_a[0][3];
      
      Hadron_a[1][2] = (-1 * Pm_c[0]*Q_c[3] + 1 * Pm_c[3]*Q_c[0]) * W3/2.;
      Hadron_a[2][1] = -Hadron_a[1][2];
      
      Hadron_a[1][3] = (1 * Pm_c[0]*Q_c[2] + -1 * Pm_c[2]*Q_c[0]) * W3/2.;
      Hadron_a[3][1] = -Hadron_a[1][3];
      
      Hadron_a[2][3] = (-1 * Pm_c[0]*Q_c[1] + 1 * Pm_c[1]*Q_c[0]) * W3/2.;
      Hadron_a[3][2] = -Hadron_a[2][3];

    for(int i=0; i<4; i++){ 
      for(int j=0; j<4; j++){
        
        Hadron[i][j] = ( Hadron_s[i][j] + I * Hadron_a[i][j] )/MN2;
        
      }
    }
    

// // // // // // // // // // // // Alternative way of computing the Antisymmetric tensor (computationally more demanding) // // // // // // // // // // // // // // // // // // // //    
//   int epsi[4][4][4][4]={};
//   epsi[0][1][2][3] = -1; epsi[0][1][3][2] = 1;
//   epsi[0][2][1][3] = 1;  epsi[0][2][3][1] = -1;
//   epsi[0][3][2][1] = 1;  epsi[0][3][1][2] = -1;
// 
//   epsi[1][0][2][3] = 1;  epsi[1][0][3][2] = -1;
//   epsi[1][2][0][3] = -1; epsi[1][2][3][0] = 1;
//   epsi[1][3][2][0] = -1; epsi[1][3][0][2] = 1;
// 
//   epsi[2][1][0][3] = 1;  epsi[2][1][3][0] = -1;
//   epsi[2][0][1][3] = -1; epsi[2][0][3][1] = 1;
//   epsi[2][3][0][1] = -1; epsi[2][3][1][0] = 1;
// 
//   epsi[3][0][2][1] = -1; epsi[3][0][1][2] = 1;
//   epsi[3][2][0][1] = 1;  epsi[3][2][1][0] = -1;
//   epsi[3][1][2][0] = 1;  epsi[3][1][0][2] = -1;
// // // // // // // // // // // // // // // // // // // //
//     
//     double Hadron_a[4][4]={};
//     for(int i=0; i<4; i++){
//       for(int j=0; j<4; j++){
// 
//         for(int al=0; al<4; al++){
//           for(int be=0; be<4; be++){
//             Hadron_a[i][j] += epsi[i][j][al][be] * Pm_c[al]*Q_c[be];
//           }
//         }
// 
//         Hadron[i][j] = ( Hadron_s[i][j] + I * W3/2. * Hadron_a[i][j] )/MN2;
//       }
//     }

    
  }//process ==1 
  
  // // // // // // // // // // // // // // // // // // // // // 
}
