#ifndef __oscarmaker_h__
#define __oscarmaker_h__
#define __JOSH_FORMAT__

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "coral.h"
#include "eqofst.h"

using namespace std;
class CRing;

class COSCARmaker{
public:
	string oscarfilename,inputfilename;
	// lambdafact is the ratio s.t., lambda_ij = pi_ij / lambdafact,
	// where pi_ij is the shear tensor in matter frame and lambda_ij tells how the momenta are scaled 
	// during collision time  p_i = ptilde_i + lambda_ij ptilde_j, where ptilde is generated isotropically
	double T,ETA_MAX; // will generate particles with -eta_max < eta < eta_max
	double lambdafact,ptbar,ptbar_norm;
	int MakeEvent();
	COSCARmaker(double T,double eta_max,string inputfilename_in,string oscarfilename_in);
protected:
	int MC_NWrite;
	double MC_Ntarget,MC_Nbar,MC_NSample,Ncheck;
	int GenerateParticles();
	CRandom *randy;
	CSpecies *species;
	double *density;
	string tmpfilename;
	FILE *tmpfile,*oscarfile,*input;
	CRing *earlier,*later;
	int nevents;
};

// This fills out Pi[4][4] given Pi_xx, Pi_yy, Pi_xy, u_x and u_y (all specified in lab frame)
void FillOutPi(double **Pi,double pixx,double piyy,double pixy,double ux,double uy);

// Stores info about a ring that breaks up at specific tau
class CRing{
public:
	CRing();
	void Read_VH2(FILE *fptr);
	void Read_Josh(FILE *fptr);
	void ReadHeader_Josh(FILE *fptr);
	void FillOutArrays_Josh(double *x,double *y,double *uxi,double *uyi,double *lambda_xx,double *lambda_xy,double *lambda_yy);
	void FillOutArrays_VH2(double *x,double *y,double *uxi,double *uyi,double *lambda_xx,double *lambda_xy,double *lambda_yy);
	int nphi;
	int nread;
	double tau;
	double *r,*ux,*uy;
	double ***lambda;
};

#endif