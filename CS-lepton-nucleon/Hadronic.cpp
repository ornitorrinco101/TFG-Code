#include "Hadronic.h"

void Hmunu( int process, int nucleon, double Pm[], double PN[], complex<double> Hadron[][4])
{
  
  double Q[4]={ PN[0]-Pm[0], PN[1]-Pm[1], PN[2]-Pm[2], PN[3]-Pm[3] };
  double Qsq = -(pow(Q[0],2) - pow(Q[1],2) - pow(Q[2],2) - pow(Q[3],2)); // defined positive here Qsq > 0 

  double F1p, F1n, F2p, F2n, GA, FP, F1s, F2s, GAs;  
  formfactors( Qsq, F1p, F2p, F1n, F2n, GA, FP, F1s, F2s, GAs );

  if(process==0) //EM interaction
  {
    double PmPN = Pm[0]*PN[0] - Pm[1]*PN[1] - Pm[2]*PN[2] - Pm[3]*PN[3];
    
    double F1, F2;
    if(nucleon == 1){F1 = F1p; F2=F2p;}
    else{F1 = F1n; F2=F2n;}
    
    for(int i=0; i<4; i++){ 
      for(int j=0; j<4; j++){
    
        Hadron[i][j] = 1/(2*MN2) * 
        ( pow(F1+F2,2) * (Pm[i]*PN[j] + Pm[j]*PN[i] + (MN2-PmPN)*gmunu[i][j])
        +( pow(F2/(2*MN),2) *( PmPN + MN2 ) - F2*(F1+F2) ) * (Pm[i]+PN[i])*(Pm[j]+PN[j]) );
        
      }
    }
        
  }//process==0
  else if(process == 1)
  {
    
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
  else if(process == 2)
  {
    
    double Pm_c[4] = {Pm[0], -Pm[1], -Pm[2], -Pm[3]};
    double Q_c[4] = {Q[0], -Q[1], -Q[2], -Q[3]};
    
    // // // // // Symmetric tensor // // // // 
    double PmPN = Pm[0]*PN[0] - Pm[1]*PN[1] - Pm[2]*PN[2] - Pm[3]*PN[3];
    
    // // Form factors // //     
    int tau3=1;  if(nucleon == 0){tau3=-1;} //proton=1, neutron=0
    
    double wF1 = (0.5-Sin2W)*(F1p-F1n)*tau3 - Sin2W*(F1p+F1n) - 0.5*F1s;
    double wF2 = (0.5-Sin2W)*(F2p-F2n)*tau3 - Sin2W*(F2p+F2n) - 0.5*F2s;

    double wGA = 0.5*(tau3*GA - GAs); 
    // // // // // // // // 
    
    
    complex<double> Hadron_VV[4][4], Hadron_AA[4][4], Hadron_VA[4][4]={};
    
    // // Symmetric tensor // // 
    for(int i=0; i<4; i++){ 
      for(int j=0; j<4; j++){
    
        // // vector-vector contribution
        Hadron_VV[i][j] = 1/(2*MN2) * 
        ( pow(wF1+wF2,2) * (Pm[i]*PN[j] + Pm[j]*PN[i] + (MN2-PmPN)*gmunu[i][j])
        +( pow(wF2/(2*MN),2) *( PmPN + MN2 ) - wF2*(wF1+wF2) ) * (Pm[i]+PN[i])*(Pm[j]+PN[j]) );
        
        // // axial-axial contribution        
        double tau=Qsq/(4*MN2);
        Hadron_AA[i][j] =  ((-1)*gmunu[i][j] * (tau  + 1 ) + Pm[i]*Pm[j]/MN2 + (Pm[i]*Q[j]+Pm[j]*Q[i])/(2*MN2) ) * wGA*wGA ;
      }
    }
    // // // // // // // // // // // // // // //   
    
    // Antisymmetric tensor (vector-axial contribution) // // 
    double PN_c[4] = {PN[0], -PN[1], -PN[2], -PN[3]};
    complex<double> temp = 2.*I*wGA*(wF1+wF2)/(2*MN2);
    //        CONVENTION:           epsilon^(0123) = -1
    Hadron_VA[0][1] = temp*(-1 * Pm_c[2]*PN_c[3] + 1 * Pm_c[3]*PN_c[2]);
    Hadron_VA[1][0] = -Hadron_VA[0][1];
    
    Hadron_VA[0][2] = temp*(1 * Pm_c[1]*PN_c[3] + -1 * Pm_c[3]*PN_c[1]);
    Hadron_VA[2][0] = -Hadron_VA[0][2];
    
    Hadron_VA[0][3] = temp*(-1 * Pm_c[1]*PN_c[2] + 1 * Pm_c[2]*PN_c[1]);
    Hadron_VA[3][0] = -Hadron_VA[0][3];
    
    Hadron_VA[1][2] = temp*(-1 * Pm_c[0]*PN_c[3] + 1 * Pm_c[3]*PN_c[0]);
    Hadron_VA[2][1] = -Hadron_VA[1][2];
    
    Hadron_VA[1][3] = temp*(1 * Pm_c[0]*PN_c[2] + -1 * Pm_c[2]*PN_c[0]);
    Hadron_VA[3][1] = -Hadron_VA[1][3];
    
    Hadron_VA[2][3] = temp*(-1 * Pm_c[0]*PN_c[1] + 1 * Pm_c[1]*PN_c[0]);
    Hadron_VA[3][2] = -Hadron_VA[2][3];
    // // // // // // // // // // // // // // //   

// // Equivalent expression in terms of Q instead of PN
//     //        CONVENTION:           epsilon^(0123) = -1
//     Hadron_VA[0][1] = temp*(-1 * Pm_c[2]*Q_c[3] + 1 * Pm_c[3]*Q_c[2]);
//     Hadron_VA[1][0] = -Hadron_VA[0][1];
//     
//     Hadron_VA[0][2] = temp*(1 * Pm_c[1]*Q_c[3] + -1 * Pm_c[3]*Q_c[1]);
//     Hadron_VA[2][0] = -Hadron_VA[0][2];
//     
//     Hadron_VA[0][3] = temp*(-1 * Pm_c[1]*Q_c[2] + 1 * Pm_c[2]*Q_c[1]);
//     Hadron_VA[3][0] = -Hadron_VA[0][3];
//     
//     Hadron_VA[1][2] = temp*(-1 * Pm_c[0]*Q_c[3] + 1 * Pm_c[3]*Q_c[0]);
//     Hadron_VA[2][1] = -Hadron_VA[1][2];
//     
//     Hadron_VA[1][3] = temp*(1 * Pm_c[0]*Q_c[2] + -1 * Pm_c[2]*Q_c[0]);
//     Hadron_VA[3][1] = -Hadron_VA[1][3];
//     
//     Hadron_VA[2][3] = temp*(-1 * Pm_c[0]*Q_c[1] + 1 * Pm_c[1]*Q_c[0]);
//     Hadron_VA[3][2] = -Hadron_VA[2][3];
//     // // // // // // // // // // // // // // //   
    
    for(int i=0; i<4; i++){ 
      for(int j=0; j<4; j++){
        
        Hadron[i][j] = Hadron_VV[i][j] + Hadron_AA[i][j] + Hadron_VA[i][j];
        
      }
    }

  }
  
    
  // // // // // // // // // // // // // // // // // // // // // 
}


















