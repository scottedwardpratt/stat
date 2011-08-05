/* 
Lapack++ Testing Example

Written 2011 Hal Canary <hal@cs.unc.edu> for MADAI <http://madai.us/>

How to compile:
  Under Linux, I installed the packages 
	blas
	blas-devel
	lapack
	lapack-devel
  Downloaded the lapackpp-2.5.4.tar.gz package and unzipped it:
	PREFIX=/playpen2/local
	./configure --prefix=$PREFIX
	make
	make install
  Then in this directory, one can
	make
	make run
*/
/************************Begin*Makefile******************************
PREFIX = /playpen2/local
CFLAGS = -I $(PREFIX)/include/lapackpp -L $(PREFIX)/lib -llapackpp -g -Wall
CXX = g++
lapackpp_example: lapackpp_example.cxx
	$(CXX) $(CFLAGS) $< -o $@
run: lapackpp_example
	LD_LIBRARY_PATH=$(PREFIX)/lib ./lapackpp_example
**************************End*Makefile******************************/
#include <blas3pp.h>
#include <gmd.h>
#include <lavd.h> 
#include <laslv.h>
#include <lavli.h>
#include <laexcp.h>

#include <time.h>
#include <iostream>
#include <vector>

using namespace std;

void print_vec_double(const std::vector<double> & vector_double) {
  int vector_size = vector_double.size();
  if (vector_size == 0) {
    cout << endl;
    return;
  }
  cout << '<';
  for (int i = 0; i < vector_size - 1 ; i++)
    cout << vector_double[i] << " ";
  cout << vector_double[vector_size - 1] << '>' << endl;
}

/** Returns a pointer to a newly allocated inverse matrix on the heap.
    If, for some reason, we can not compute an inverse matrix (id est
    it does not exist or it is not computable using LU factorization)
    then we return NULL.  Caller should check for that. */
LaGenMatDouble * invert(LaGenMatDouble & mat) {
  int N = mat.size(0);
  assert(mat.size(1) == N);
  LaGenMatDouble * mat2 = new LaGenMatDouble(mat);
  LaVectorLongInt pivot = LaVectorLongInt(N);
  try {
    LUFactorizeIP(*mat2,pivot);
    LaLUInverseIP(*mat2,pivot);
  } catch (LaException& e) {
    //cout << e.what() << endl;
    delete mat2;
    return NULL;
  }
  return mat2;
}

/** multiplies a vector by a  */

void transform(LaGenMatDouble & trans,
	       std::vector<double> & ivec,
	       std::vector<double> & ovec) {
  int N = static_cast<int>(ovec.size());
  assert(static_cast<int>(ivec.size()) == N);
  assert(trans.size(0) == N);
  assert(trans.size(1) == N);
  for (int row = 0; row < N; row++) {
    double sum = 0.0;
    for (int column = 0; column < N; column++) {
      sum += (trans(row,column) * ivec[column]);
    }
    ovec[row] = sum;
  }
  return;
}

/** Checks to see that the two matrices multiply to identity. */
bool isInverse(LaGenMatDouble & m1, LaGenMatDouble & m2){
  int N = m1.size(0);
  assert(m1.size(1) == N);
  assert(m2.size(0) == N);
  assert(m2.size(1) == N);
  for (int row = 0; row < N; row++) {
    for (int column = 0; column < N; column++) { 
      double sum = 0.0;
      for (int k = 0; k < N; k++)
	sum += (m1(row,k) * m2(k,column));
      if (row == column)
	sum -= 1.0;
      if (abs(sum) > 1.0e-10) 
	return false;
    }
  }
  return true;
}

double test_equality(const std::vector<double> & dv1,
	       const std::vector<double> & dv2) {
  double biggest_diff = 0.0;
  unsigned int N = dv1.size();
  assert(dv2.size() == N);
  for (unsigned int i = 0; i< N; i++) {
    double x = abs(dv1[i] - dv2[i]);
    if (x > biggest_diff)
      biggest_diff = x;
  }
  return biggest_diff;
}

double test(int N) {
  vector<double> dv1(N);
  vector<double> dv2(N);
  vector<double> dv3(N);
  srand ( time(NULL) );
  LaGenMatDouble mat = LaGenMatDouble::rand(N,N,-100.0,100.0);
  LaGenMatDouble * matinv = invert(mat);

  if (matinv == NULL) {
    cerr << "Can't find inverse.\n";
    return -1.0;
  }
  for (int i = 0; i< N; i++)
    dv1[i] = i + 1;

  transform(mat, dv1, dv2);
  transform(*matinv, dv2, dv3);
  return test_equality(dv1, dv3);
}

int main (int argc, char *argv[]) {
  int arg;
  int N = 1536;
  if (argc > 1)
    if ((arg = atoi(argv[1])) > 0)
      N = arg;
  cout << "N = " << N << endl;
  cout << "biggest diff (small is good) = " << test(N) << endl;
}
