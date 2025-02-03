#include "flux.h"
// #include "parameters.h"

// Outputs: flux[], enu_flux[], iflux_max
void readtheflux( double enu_flux[], double flux[], int &iflux_max )
{
  
  ifstream influx;
  influx.open("t2kflux_numu_2016.dat");
  if (!influx) {
      cout << "flux.cpp: Unable to open flux input file" << endl;
      exit(1);   // call system to stop
  }
  
 
  //prepared to the T2K flux (t2kflux_numu_2016.dat), the following lines must be adapted to each particular input file
  iflux_max=0;
  int i=0;
  
  while(i<201){
    influx >> enu_flux[i] >> flux[i]; // enu_flux[i] in MeV // now flux[i] in MeV^-1
    influx.ignore(1000, '\n');  // line skipper
    iflux_max = i;
//       cout << iflux_max << ", enu (MeV)=" << enu_flux[i] << ", flux= " << flux[i] << endl;
    i++;
  }
  // // // // // // // // // // // 
  
  influx.close();  
  cout << "flux.cpp: iflux_max= " << iflux_max << ", enu_flux[iflux_max]= " << enu_flux[iflux_max] << " " << flux[iflux_max] << endl;
  
}


//Inputs: Ei, enu_flux[], flux[], iflux_max. Outputs: flujo, saltador1
void findtheflux(double Ei, double enu_flux[], double flux[], int iflux_max, double &flujo, int &stop)
{
  if( (Ei > 0) && (Ei < enu_flux[iflux_max]) )
  { 
    
    if( enu_flux[0] >= Ei ){ //Enu is small so we take this as the value
      flujo = flux[0];            
    }else{
      int i=0;
      while( enu_flux[i] < Ei ){i++;}            
      if( enu_flux[i] == Ei ){
        
        flujo = flux[i];
        
      }else{ // linear interpolation
        
        double enu_low = enu_flux[i-1]; 
        double enu_up = enu_flux[i];
        double flujo_low = flux[i-1];
        double flujo_up = flux[i];   

        double fac = (Ei - enu_low)/(enu_up - enu_low);
        flujo = flujo_low + (flujo_up - flujo_low)*fac; //flujo ~ MeV^-1
            
//         if(flujo < 0 ){
//             cout << Ei << " " << enu_low << " " << enu_up << " " << flujo_low << " " << flujo_up << endl;
//         }

        }
    }
      
  }else{ stop = 1;}
  
}
// // // // // // // // // // // // // // // // // // // // // // // 
