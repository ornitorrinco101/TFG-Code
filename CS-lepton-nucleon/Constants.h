#ifndef CONSTANTS_H
#define CONSTANTS_H

/*
The "Constants.h" header file, when included, provides some useful constant objects that will be used throughout the code.
*/

#include <complex>
using namespace std;

// The complex number unit
const complex<double> I = complex<double> (0, 1);

// Some real constants
const double Mp = 938.27203;  // proton mass (MeV)
// const double Mn = 939.56536;  // neutron mass (MeV)
// const double MN = 938.918695;  // average
const double MN = 938.918695;
// const double Mp = 938.918695;
const double MN2 = 881568.315821;//MN^2

// const double Mpi_chrgd = 139.57018;  // charged-pion mass (MeV)
// const double Mpi_ntrl = 134.9766;  // neutral pion mass (MeV)
const double Mpi = 138.0389867; // average (2*Mpi_chrgd+Mpi_ntrl)/3
const double Mpi2 = 19054.761849; //Mpi^2

const double electronmass = 0.511; // MeV
const double muonmass = 105.658369;  // muon mass (MeV)

const double Pi = 3.141592654;  // \pi
const double hbc = 197.3270;  // hbar*c (MeV*fm)

// Electroweak constants
const double M_W = 80385.;  // W-boson mass (MeV)
const double M_Z = 91187.6; // Z-boson mass (MeV)
const double Cabibbo = 0.974;  // cosine Cabibbo angle
const double G_Fermi = 1.16637e-11;  // Fermi constant (MeV^{-2})
const double Alpha = 0.007297353;  // fine structure constant
const double Sin2W = 0.23122; 
const double QWeak = 0.07512; //weak charge of the proton = 1-4*Sin2W
const double gA = 1.2695;
const double fpi = 93.;
// // // // // // // // // // // // // 

// // // Metric // // // 
const double gmunu[4][4] = { {1,0,0,0}, {0,-1,0,0}, {0,0,-1,0}, {0,0,0,-1} };


#endif
