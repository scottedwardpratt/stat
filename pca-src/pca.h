#ifndef __PCA_H__
#define __PCA_H__
#include "coral.h"
#include "qualifier.h"
using namespace std;

class CPCA{
public:
	string yname[10000];
	double sigmay[10000];
	int nruns,ny,nnames;
	CQualifiers qualifiers;
	CPCA(int nruns_set);
	double *ybar,*value,**spread;
	void ReadResults();
	double *eigenval;
	double **evec;
	CGSLMatrix_Real *gslmatrix;
	bool namecheck(string varname);
	string pcaname[10000];
	void Calc();
};

#endif

