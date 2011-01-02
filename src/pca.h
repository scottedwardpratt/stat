#ifndef __PCA_H__
#define __PCA_H__
#include "coral.h"
#include "qualifier.h"
using namespace std;

class CPCA{
public:
	string yname[1000];
	double sigmay[1000];
	int nruns,ny;
	CQualifiers qualifiers;
	CPCA(int nruns_set);
	double *ybar,*value,**spread;
	void ReadResults();
	double *eigenval;
	double **evec;
	CGSLMatrix_Real *gslmatrix;
	void Calc();
};

#endif

