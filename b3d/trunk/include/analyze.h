#ifndef __ANALYZE_H__
#define __ANALYZE_H__

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <sys/stat.h>
#include "coral.h"
#include "H5Cpp.h"
#include "b3d.h"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;

/*
class CPartH5{
public:
	int listid,ID;
	double px,py,rapidity,mass,x,y,eta,tau;
};
*/

class CAnalyze{
public:
	string qualifier;
	int npartsmax,neventsmax,nsample;
	parameterMap parmap;
	CompType *ptype;
	//CompType ptype(sizeof(CPartH5));
	CAnalyze(string qualifier,string parsfilename);
	string input_dataroot,output_dataroot;
	string h5_infilename,outfilename;
	H5File *h5file;
	CPartH5 *partH5;
	bool CALCGARRAYS;
	bool STAR_ACCEPTANCE;
	int ReadDataH5(int ievent); // returns nparts for given event
	void CalcSpectra();
	void CalcV2();
	void CalcHBT();
	void CalcGamma();
	void CalcGamma_BothSigns();
	void CalcBalance();
};



#endif
