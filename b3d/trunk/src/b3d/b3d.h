#ifndef __B3D_H__
#define __B3D_H__

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <complex>
#include <cstdio>
#include <sys/stat.h>
#include "coral.h"
#include "H5Cpp.h"
#include "resonances.h"
#include "bjmaker.h"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;

class CPart; class CAction; class CB3D; class CCell; class CPartH5; class CHYDROtoB3D;

typedef multimap<int,CPart *> CPartMap;
typedef multimap<double,CAction *> CActionMap;
typedef pair<int,CPart*> CPartPair;
typedef pair<double,CAction*> CActionPair;

class CB3D{
public:
	parameterMap parmap;
	CPartMap DeadPartMap,PartMap,FinalPartMap;
	CActionMap ActionMap,DeadActionMap;
	CResList *reslist;
	CPart **part;
	CAction **action;

	int NXY,NETA; // will make meshes of size (2NXmesh+1,2NYmesh+1,2NETAmesh+1)
	int NRINGSMAX;
	double XYMAX,ETAMAX,DXY,DETA;
	CHYDROtoB3D *hydrotob3d;
	CBjMaker bjmaker;
	bool BJORKEN,COLLISIONS,VIZWRITE;
	CCell ****cell;

	CB3D(string run_name_set);
	~CB3D();
	double tau,TAUCOLLMAX;
	int itau;
	int ievent_write,ievent_read;
	//
	// READ IN FROM PARAMETER FILE
	int NACTIONSMAX;
	int NPARTSMAX;
	double SIGMAMAX,SIGMADEFAULT; // cross sections in sq. fm
	string input_dataroot;
	string output_dataroot;
	string run_name,qualifier;
	CompType *ptype;

	H5File *h5outfile, *h5infile;// *h5vizfile;
	int NACTIONS;
	int NSAMPLE;
	//
	void SetQualifier(string qualifier_set);
	int ReadDataH5(int ievent);
	int WriteDataH5();
	void FindAllCollisions();
	void FindAllCellExits();
	void PerformAllActions();
	void Reset();
	void KillAllActions();
	void KillAllParts();

	void AddAction_Activate(CPart *part);
	void AddAction_Decay(CPart *part);
	void AddAction_Collision(CPart *part1,CPart *part2,double tau);
	void AddAction_ResetCollisions(double taureset);
	//void AddAction_SwallowParticles(double tau_breakup);
	void AddAction_ExitCell(CPart *part);
	void AddAction_VizWrite(double tauwrite);

	void ListFutureCollisions();
	void PrintPartList();

	bool FindCollision(CPart *part1,CPart *part2,double &taucoll);
	void Decay(CPart *&mother,int &nbodies, CPart **&daughter);

	CRandom *randy;

	void PrintActionMap(CActionMap *actionmap);

	double GetPiBsquared(CPart *part1,CPart *part2);
	int Collide(CPart *part1,CPart *part2,double scompare);// will collide if sigma>scompare
	void Scatter(CPart *part1,CPart *part2);
	void Merge(CPart *part1,CPart *part2,CResInfo *resinfo);

	void CheckActions();
	bool ERROR_PRINT;

	long long int nscatter,ndecay,nmerge,nswallow,npass,nexit,nactivate,ncheck;
	long long int nactions;

	// These are used for VizWrite
	hid_t viz_file_id;
	//

	void freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt);

};

class CPart{
public:
	CPart();
	CCell *cell,*nextcell;
	double tau0,tau_lastint,tauexit;
	double y,eta;
	double p[4],r[4],mass;
	int listid;
	int actionmother; //refers to action from which particle was created
	CResInfo *resinfo;

	void Propagate(double tau);
	void FindCellExit();
	void Init(int ID,double x,double y,double tau,double eta,double px,double py,double mass,double rapidity);
	void Init_NoColls(int ID,double x,double y,double tau,double eta,double px,double py,double mass,double rapidity);
	void InitH5(CPartH5*);
	void Setp0();
	double GetMass();
	void CyclicReset();
	int key;
	void SetInitialKey();

	void Kill();
	void Reset();
	void SubtractAction(CAction *actionptr);
	void AddAction(CAction *actionptr);
	void KillActions();
	void Print();
	void CheckMap(CPartMap *expectedpartmap);
	void ChangeMap(CPartMap *newmap);
	//~CPart();

	// These are the actions involving these particles
	CActionMap actionmap;

	CPartMap *currentmap; // PartList for a Cell, or b3d->DeadPartList
	CCell *FindCell();


	static CB3D *b3d;
	double GetEta(double tau);
	double GetMT();
	void SetY();
	void SetEta(double neweta);
	void FindCollisions();

	CPartMap::iterator GetPos(CPartMap *pmap);
	CPartMap::iterator DeleteFromMap(CPartMap *partmap);
	void ChangePartMap(CPartMap *newmap);
	CPartMap::iterator DeleteFromCurrentMap();
	void AddToMap(CPartMap *newmap);
	void AddToMap(CPartMap::iterator guess,CPartMap *newmap);
	bool active;
};

class CAction{
public:
	double tau;
	int listid;
	double key;
	void SetKey();
	int type; // =0 for activation, 1 for collision, 2 for decay, ....  6 for ExitCell
	// These are the particles in the action
	CPartMap partmap;

	void Kill();
	void CleanPartMap();
	void AddPart(CPart *partptr);
	void Print();

	void Perform();
	void PerformActivate();
	void PerformExitCell();
	void PerformDecay();
	void PerformCollide();
	void PerformResetCollisions();
	void PerformVizWrite();
	//void PerformSwallowParticles();
	CAction(CB3D *b3dset);
	~CAction();

	static CB3D *b3d;

	CActionMap::iterator GetPos(CActionMap *actionmap);
	void ChangeMap(CActionMap *newmap);
	CActionMap::iterator DeleteFromCurrentMap();
	void AddToMap(CActionMap *newmap);
	void AddToMap(CActionMap::iterator guess,CActionMap *newmap);
	void CheckPartList();
	CActionMap *currentmap;
};

class CCell{
public:
	class CCell *neighbor[3][3][3];
	int ix,iy,ieta;
	class CCell *creflection;
	int ireflection;
	double xmin,xmax,ymin,ymax,etamin,etamax;
	CPartMap partmap;
	void PrintPartMap(CPartMap *partmap);
	void KillAllParts();
	void ReKeyAllParts();
	void Print();

	CCell(double xmin,double xmax,double ymin,double ymax,double etamin,double etamax);
	static CB3D *b3d;
};

class CPartH5{
public:
	int listid,ID;
	double px,py,rapidity,mass,x,y,eta,tau;
};

class CRing; class CResList;

class CHYDROtoB3D{
public:
	//CHYDROtoB3D();
	// lambdafact is the ratio s.t., lambda_ij = pi_ij / lambdafact,
	// where pi_ij is the shear tensor in matter frame and lambda_ij tells how the momenta are scaled 
	// during collision time  p_i = ptilde_i + lambda_ij ptilde_j, where ptilde is generated isotropically
	double T,ETAMAX;
	bool initialization;
	// will generate particles with -eta_max < eta < eta_max
	CResList *reslist;
	double epsilon,P,lambdafact;
	int MakeEvent();
	int GenerateParticles(int iring1,int iring2);
	void freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt);
	void Init();
	void ReadInput();
	CB3D *b3d;
	double GetLambdaFact();
protected:
	int MC_NWrite;
	double MC_Ntarget,MC_Nbar,nsample,Ncheck;
	CRandom *randy;
	int nres,nrings;
	double *density,*ID;
	string tmpfilename;
	FILE *tmpfile,*oscarfile,*input;
	CRing *ring;
	void GetRingInfo();
	void ReadHeader(FILE *);
};

// This fills out Pi[4][4] given Pi_xx, Pi_yy, Pi_xy, u_x and u_y (all specified in lab frame)
void FillOutPi(double **Pi,double pixx,double piyy,double pixy,double ux,double uy);

// Stores info about a ring that breaks up at specific tau
class CRing{
public:
	CRing();
	void Read(FILE *fptr);
	void FillOutArrays(double *x,double *y,double *uxi,double *uyi,double *lambda_xx,double *lambda_xy,double *lambda_yy);
	int nphi;
	int nread;
	double tau;
	double *r,*ux,*uy;
	double ***lambda;
};

#endif
