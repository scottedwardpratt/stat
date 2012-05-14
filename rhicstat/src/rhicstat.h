#ifndef __PCA_H__
#define __PCA_H__
#include "coralutils.h"
#include "qualifier.h"
using namespace std;

class CPCA{
public:
	string yname[100];
	double sigmay[100];
	int nruns,ny,nnames;
	CQualifiers qualifiers;
	CPCA(int nruns_set);
	double *ybar,*value,**spread;
	void ReadResults();
	double *eigenval;
	double **evec;
	CGSLMatrix_Real *gslmatrix;
	bool namecheck(string varname);
	void Calc();
};

class CRunInfo{  /** this is info that pertains to a specific run */
public:
	CRunInfo(int NX,int NY);
	double *x,*w;
	double *y,*z;
	double *sigmay;
	double *ylinear,*zlinear,*xlinear;
};

class CRHICStat{
public:
	CRHICStat(int NRUNS);
	string *yname;
	string *xname;
	double *xmin,*xmax;
	int NX,NY,NRUNS;
	double *xbar,*ybar,*sigmaybar;
	double **sigmaxx,**sigmayy,**dxdz,**dxdz_inv;
	double **Uytoz,**Uytoz_inv,**Uxtow,**Uxtow_inv;
	double *eigenvalyy,*eigenvalxx;
	double *uncertainty;
	CRunInfo **runinfo;
	CRunInfo *expinfo,*fitinfo;
	CGSLMatrix_Real *gslmatrix_NY;
	CGSLMatrix_Real *gslmatrix_NX;
	void FitExpData();
	void ScaleXY();
	void InitArrays();
	void InitX();
	void InitY();
	void ReadAllX();
	void ReadX(string filename,CRunInfo *runinfo);
	void ReadAllY();
	void ReadY(string filename,CRunInfo *runinfo);
	void CalcSensitivity();
	void GetZFromY(CRunInfo *runinfo);
	void GetYFromZ(CRunInfo *runinfo);
	void GetXlinearFromZ(CRunInfo *runinfo);
	void GetZlinearFromX(CRunInfo *runinfo);
	void GetZlinearFromXlinear(CRunInfo *runinfo);
	void GetXlinearFromY(CRunInfo *runinfo);
	void GetYlinearFromX(CRunInfo *runinfo);
	void GetXlinearFromZlinear(CRunInfo *runinfo);
	void GetYlinearFromZlinear(CRunInfo *runinfo);
	void GetWFromXlinear(CRunInfo *runinfo);
	void GetXlinearFromW(CRunInfo *runinfo);
	void PrintXlinear(CRunInfo *runinfo);
	void PrintX(CRunInfo *runinfo);
	void PrintY(CRunInfo *runinfo);
	void PrintZ(CRunInfo *runinfo);
	void PrintZlinear(CRunInfo *runinfo);
	void PrintYlinear(CRunInfo *runinfo);
	void PlotZvsX();
};



#endif

