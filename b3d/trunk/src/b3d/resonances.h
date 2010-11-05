#ifndef __RESONANCES_H__
#define __RESONANCES_H__

#include "b3d.h"

using namespace std;

class CResInfo;

class CB3D;

class CBranchInfo{
public:
	CBranchInfo *nextbptr;
	int nbodies;
	CResInfo *resinfoptr[5];
	double branching;
	CBranchInfo();
};

class CMerge{
public:
	CMerge(CResInfo *resinfo,double branching);
	CResInfo *resinfo;
	double branching;
	CMerge *next;
};

class CResInfo{
public:
	int ires;
	double mass,spin,width,minmass;
	string name;
	int code,charge,strange,baryon,count;
	bool decay;
	CResInfo *nextResInfoptr;
	CBranchInfo *firstbptr,*bptr_minmass;
	void Print();
	void DecayGetResInfoptr(int &nbodies,CResInfo **&daughterresinfoptr);
	void DecayGetResInfoptr_minmass(int &nbodies,CResInfo **&daughterresinfoptr);
	CResInfo();
	static CRandom *ranptr;
};

class CResList{
public:
	int NResonances;
	CResList();
	CResInfo *GfirstResInfoptr;
	void GetResInfoptr(int code,CResInfo *&resinfoptr);
	void ReadResInfo();
	void PrintYields();
	void CalcEoS(double T_0,double T_F,double delT);
	void freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt);
	parameterMap *parmap;
	CMerge ***MergeArray;
	static CB3D *b3d;
};


#endif
