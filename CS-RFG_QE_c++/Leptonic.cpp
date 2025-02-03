#include "Leptonic.h"

/***************************************

    *     THE LEPTON TENSOR     *
// // down indices
***************************************/

void Lmunu( int process, int Helicity, double Ki[], double Kf[], complex<double> Lepton[][4] )
{

  double Ki_Kf = Ki[0]*Kf[0] - Ki[1]*Kf[1] - Ki[2]*Kf[2] - Ki[3]*Kf[3];

  complex<double> Lepton_S[4][4];
    
  for(int i=0; i<4; i++)
  {
    for(int j=0; j<4; j++)
    {
      Lepton_S[i][j] = 0.5*( Ki[i]*Kf[j] + Ki[j]*Kf[i] - gmunu[i][j]*Ki_Kf );
    }
  }
  
  for(int i=0; i<4; i++)
  {
    for(int j=0; j<4; j++)
    {
      Lepton[i][j] = Lepton_S[i][j];
    }
  }
  
  if( process != 0 ){
    
    complex<double> Lepton_A[4][4];
    
    //        CONVENTION:           epsilon_(0123) = 1
    
    // // L_ij = (-I) * epsilon_ijkm Ki^k Kf^m  
    
    Lepton_A[0][0] = 0.;
    Lepton_A[0][1] = (-I)*(Ki[2]*Kf[3] - Kf[2]*Ki[3]); // zero
    Lepton_A[0][2] = (-I)*(-Ki[1]*Kf[3] + Kf[1]*Ki[3]); // Kf*k_l_inc*sin(thetal) --> 0 (for small angles) 
    Lepton_A[0][3] = (-I)*(Ki[1]*Kf[2] - Kf[1]*Ki[2]); // zero
    
    Lepton_A[1][0] = -Lepton_A[0][1];
    Lepton_A[1][1] = 0.;
    Lepton_A[1][2] = (I)*(Ki[0]*Kf[3] - Kf[0]*Ki[3]); // k_l_inc*Kf*cos(thetal) - El*k_l_inc --> 0 (for small angles) 
    Lepton_A[1][3] = (I)*(-Ki[0]*Kf[2] + Kf[0]*Ki[2]); // zero
    
    Lepton_A[2][0] = -Lepton_A[0][2]; 
    Lepton_A[2][1] = -Lepton_A[1][2];
    Lepton_A[2][2] = 0.;
    Lepton_A[2][3] = (I)*(Ki[0]*Kf[1] - Kf[0]*Ki[1]); // k_l_inc*Kf*sin(thetal) --> 0 (for small angles)
    
    Lepton_A[3][0] = -Lepton_A[0][3]; // zero
    Lepton_A[3][1] = -Lepton_A[1][3]; // zero
    Lepton_A[3][2] = -Lepton_A[2][3]; 
    Lepton_A[3][3] = 0.;

      
    for(int i=0; i<4; i++)
    {
      for(int j=0; j<4; j++)
      {
        
        Lepton[i][j] = 4.*Lepton_S[i][j] + 2.*(-Helicity)*Lepton_A[i][j];  //adding the symmetric and antisymmetric parts
        
      }
    }
    
  }//process != 0


}
