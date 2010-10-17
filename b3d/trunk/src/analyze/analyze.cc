#ifndef __ANALYZE_CC__
#define __ANALYZE_CC__

#include "analyze.h"

CAnalyze::CAnalyze(string qualifier_set,string analpars_filename){
	ptype=new CompType(sizeof(CPartH5));
	qualifier=qualifier_set;
	parameter::ReadParsFromFile(parmap,analpars_filename);
	neventsmax=parameter::getI(parmap,"B3D_NEVENTS",40000);
	npartsmax=parameter::getI(parmap,"B3D_NPARTSMAX",3000);
	nsample=parameter::getI(parmap,"B3D_NSAMPLE",1);
	npartsmax*=nsample;
	input_dataroot=parameter::getS(parmap,"B3D_INPUT_DATAROOT","data/b3d");
	output_dataroot=parameter::getS(parmap,"B3D_OUTPUT_DATAROOT","data/b3d");
	h5_infilename=parameter::getS(parmap,"B3D_H5_INFILENAME","b3d.h5");	
	partH5=new CPartH5[npartsmax];
	STAR_ACCEPTANCE=parameter::getB(parmap,"B3D_STAR_ACCEPTANCE","false");
	CALCGARRAYS=false;
	//CompType ptype( sizeof(CPart) );
	ptype->insertMember("listid", HOFFSET(CPartH5,listid), PredType::NATIVE_INT);
	ptype->insertMember("ID", HOFFSET(CPartH5,ID), PredType::NATIVE_INT);
	ptype->insertMember("x", HOFFSET(CPartH5,x), PredType::NATIVE_DOUBLE);
	ptype->insertMember("y", HOFFSET(CPartH5,y), PredType::NATIVE_DOUBLE);
	ptype->insertMember("eta", HOFFSET(CPartH5,eta), PredType::NATIVE_DOUBLE);
	ptype->insertMember("tau", HOFFSET(CPartH5,tau), PredType::NATIVE_DOUBLE);
	ptype->insertMember("px", HOFFSET(CPartH5,px), PredType::NATIVE_DOUBLE);
	ptype->insertMember("py", HOFFSET(CPartH5,py), PredType::NATIVE_DOUBLE);
	ptype->insertMember("rapidity", HOFFSET(CPartH5,rapidity), PredType::NATIVE_DOUBLE);
	ptype->insertMember("mass", HOFFSET(CPartH5,mass), PredType::NATIVE_DOUBLE);
}

int CAnalyze::ReadDataH5(int ievent){
	int nparts=0;

	char eventno[20];
	H5D_space_status_t status;
	sprintf(eventno,"event%d",ievent);
	hsize_t dim[1];
	DataSet *dataset = new DataSet (h5file->openDataSet(eventno));
	//dataset->getSpaceStatus(status);
	//hsize_t datasetsize=dataset->getStorageSize();
	//printf("ievent=%d, status=%d, size=%d\n",ievent,int(status),int(datasetsize));
	DataSpace filespace = dataset->getSpace();
	int rank=filespace.getSimpleExtentDims(dim);
	nparts=dim[0];
	if(nparts>npartsmax){
		printf("Increase NPARTSMAX, nparts=%d\n",nparts);
		exit(1);
	}
	dataset->read(partH5,*ptype);
	delete dataset;
	//printf("READ IN %d PARTS\n",nparts);
	return nparts;
}

#endif