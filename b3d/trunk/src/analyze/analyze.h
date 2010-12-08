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

class CAnalyze{
public:
	string qualifier,run_name;
	int npartsmax,neventsmax,nsample;
	parameterMap parmap;
	CompType *ptype;
	//CompType ptype(sizeof(CPartH5));
	CAnalyze(string run_name);
	void SetQualifier(string qualname);
	string input_dataroot,output_dataroot;
	string h5_infilename,outfilename,vizfilename;
	H5File *h5file,*vizfile;
	CPartH5 *partH5;
	bool CALCGARRAYS;
	bool STAR_ACCEPTANCE;
	int ReadDataH5(int ievent); // returns nparts for given event
	int ReadVizData(double tau);
	void CalcSpectra();
	void CalcV2();
	void CalcHBT();
	void GetHBTpars(CPartH5 *pH5,double &tau,double &rout,double &rside,double &rlong);
	void CalcGamma();
	void CalcGamma_BothSigns();
	void CalcBalance();
};



#endif
