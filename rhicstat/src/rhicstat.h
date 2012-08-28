#ifndef __RHICSTAT_H__
#define __RHICSTAT_H__
#include "coralutils.h"
#include "qualifier.h"
using namespace std;

class CRunInfo{  /** this is info that pertains to a specific run */
public:
	CRunInfo(int NX,int NY);
	double *x,*w;
	double *y,*z;
	double *sigmay;
	double *ylinear,*zlinear,*xlinear,*zquad,*yquad,*xquad;
	double netdiff_exp,netdiff_quad,netdiff_quadexp;
	bool good;
};

class CRHICStat{
public:
	CRHICStat(int NRUNS,int NTESTRUNS);
	string *yname;
	string *xname;
	double *xmin,*xmax;
	int NX,NY,NRUNS,NTESTRUNS,NGOODRUNS;
	double *xbar,*ybar,*sigmaybar;
	double **sigmaxx,**sigmayy,**dxdz,**dxdz_inv;
	double **Uytoz,**Uytoz_inv,**Uxtow,**Uxtow_inv;
	double *eigenvalyy,*eigenvalxx;
	double *uncertainty;
	double ***Aquad,**Bquad,*Cquad;
	CRunInfo **runinfo,**testinfo;
	CRunInfo *expinfo,*fitinfo,*bestinfo;
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
	void QuadFit();
	void GetZquad(double *x,double *z);
	void GetZFromY(CRunInfo *runinfo);
	void GetYFromZ(CRunInfo *runinfo);
	void GetXlinearFromZ(CRunInfo *runinfo);
	void GetZlinearFromX(CRunInfo *runinfo);
	void GetZlinearFromXlinear(CRunInfo *runinfo);
	void GetXlinearFromY(CRunInfo *runinfo);
	void GetYlinearFromX(CRunInfo *runinfo);
	void GetXlinearFromZlinear(CRunInfo *runinfo);
	void GetYlinearFromZlinear(CRunInfo *runinfo);
	void GetYquadFromZquad(CRunInfo *runinfo);
	void GetWFromXlinear(CRunInfo *runinfo);
	void GetXlinearFromW(CRunInfo *runinfo);
	void PrintXlinear(CRunInfo *runinfo);
	void PrintX(CRunInfo *runinfo);
	void PrintY(CRunInfo *runinfo);
	void PrintZ(CRunInfo *runinfo);
	void PrintZlinear(CRunInfo *runinfo);
	void PrintYlinear(CRunInfo *runinfo);
	void PrintZquad(CRunInfo *runinfo);
	void PrintYquad(CRunInfo *runinfo);
	void CheckTestRuns();
	void PlotZvsX();
	void PrintQuadCode();
	double GetLL(double *x);
	void Metropolis(unsigned int nburn,unsigned int nsample);
	void PrintCoefficients();
	void CalcError(CRunInfo*);
	void CalcNetDiffExp(CRunInfo *runinfo);
	void CalcNetDiffQuad(CRunInfo *runinfo);
	void CalcNetDiffQuadExp(CRunInfo *runinfo);
};



#endif

