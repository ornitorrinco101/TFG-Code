#ifndef LEPTONIC_H
#define LEPTONIC_H

#include "Constants.h"

/*****************************************************************************

*                            THE LEPTON TENSOR                               *

*****************************************************************************/

void Lmunu( int process, int Helicity, double Ki[], double Kf[], complex<double> Lepton[][4] );

#endif
