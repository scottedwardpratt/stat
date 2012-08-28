# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;
# include "latin_center.h"

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    LATIN_CENTER_PRB tests the Latin Center Square routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int seed,NDIM=2,NPTS=10;
  //seed = get_seed ( );
  seed=1234;
  cout << "\n";
  cout << "LATIN_CENTER_PRB, C++ version\n";
  double *x=new double[NDIM*NPTS];
  latin_center ( NDIM, NPTS, &seed, x );
  int i,j,k,kk;
  k = 0;
  for ( j = 0; j < NPTS; j++ )
  {
    kk = k;
    for ( i = 0; i < NDIM; i++ )
    {
      cout << setw(10) << x[kk] << "  ";
      kk = kk + NPTS;
    }
    cout << "\n";
    k = k + 1;
  }
  
  return 0;
}

