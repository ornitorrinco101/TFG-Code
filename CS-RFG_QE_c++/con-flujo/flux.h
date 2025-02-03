#ifndef FLUX_H
#define FLUX_H

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void readtheflux( double enu_flux[], double flux[], int &iflux_max );

void findtheflux( double Ei, double enu_flux[], double flux[], int iflux_max, double &flujo, int &stop);

#endif
