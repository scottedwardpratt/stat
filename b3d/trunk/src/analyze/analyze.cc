#ifndef __ANALYZE_CC__
#define __ANALYZE_CC__

#include "analyze.h"

CAnalyze::CAnalyze(string run_name_set){
	run_name=run_name_set;
	string parsfilename="parameters/"+run_name+"/fixed.param";
	parameter::ReadParsFromFile(parmap,parsfilename);
	parsfilename="parameters/"+run_name+"/stats.param";
	parameter::ReadParsFromFile(parmap,parsfilename);

	ptype=new CompType(sizeof(CPartH5));
	
	neventsmax=parameter::getI(parmap,"B3D_NEVENTS",40000);
	npartsmax=parameter::getI(parmap,"B3D_NPARTSMAX",3000);
	nsample=parameter::getI(parmap,"B3D_NSAMPLE",1);
	npartsmax*=nsample;
	
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
	
	partH5=new CPartH5[npartsmax];
	STAR_ACCEPTANCE=parameter::getB(parmap,"B3D_STAR_ACCEPTANCE","false");
	CALCGARRAYS=false;
	//CompType ptype( sizeof(CPart) );

}

void CAnalyze::SetQualifier(string qualifiername){
	input_dataroot="output/"+run_name+"/"+qualifiername;
	output_dataroot="analysis/"+run_name+"/"+qualifiername;
	string command="mkdir -p "+output_dataroot;
	system(command.c_str());
	h5_infilename=input_dataroot+"/b3d.h5";
	//if(h5file!=NULL) delete h5file;
	//h5_infilename=parameter::getS(parmap,"B3D_H5_INFILENAME","b3d.h5");	
	//vizfilename=parameter::getS(parmap,"B3D_VIZ_INFILENAME","b3dviz.h5");
	//vizfile = new H5File(infilename,H5F_ACC_RDONLY);
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

int CAnalyze::ReadVizData(double tau){
	int nparts,ipart,rank;
	string infilename=input_dataroot+"/"+qualifier+"/"+vizfilename;
	printf("will read %s for viz data\n",infilename.c_str());
	vizfile = new H5File(infilename,H5F_ACC_RDONLY);

	char setname[40];
	H5D_space_status_t status;
	sprintf(setname,"px_tau%g",tau);
	hsize_t pdim[1];
	DataSet *dataset = new DataSet (vizfile->openDataSet(setname));
	//dataset->getSpaceStatus(status);
	//hsize_t datasetsize=dataset->getStorageSize();
	//printf("ievent=%d, status=%d, size=%d\n",ievent,int(status),int(datasetsize));
	DataSpace filespace = dataset->getSpace();
	rank=filespace.getSimpleExtentDims(pdim);
	nparts=pdim[0];
	printf("For px: nparts=%d\n",nparts);
	double *px=new double[nparts];
	if(nparts>npartsmax){
		printf("Increase NPARTSMAX, nparts=%d\n",nparts);
		exit(1);
	}
	dataset->read(px,PredType::NATIVE_DOUBLE);
	delete dataset;
	
	
	printf("------------------\n");
	sprintf(setname,"xyz_tau%g",tau);
	DataSet *xyzdataset = new DataSet (vizfile->openDataSet(setname));
	//dataset->getSpaceStatus(status);
	//hsize_t datasetsize=dataset->getStorageSize();
	//printf("ievent=%d, status=%d, size=%d\n",ievent,int(status),int(datasetsize));
	hsize_t xyzdim[]={nparts,3};
	DataSpace xyzfilespace = xyzdataset->getSpace();
	rank=xyzfilespace.getSimpleExtentDims(xyzdim);
	printf("xyz rank=%d\n",rank);
	nparts=xyzdim[0];
	printf("For Reading set %s, nparts=%d, dimension=%d\n",setname,nparts,int(xyzdim[1]));
	double (*xyz)[3]=new double[nparts][3];
	if(nparts>npartsmax){
		printf("Increase NPARTSMAX, nparts=%d\n",nparts);
		exit(1);
	}
	xyzdataset->read(xyz,PredType::NATIVE_DOUBLE);
	delete xyzdataset;
	
	
	printf("READ IN %d PARTS\n",nparts);
	for(ipart=0;ipart<nparts;ipart++){
		printf("ipart=%d: px=%g, xyz=(%g,%g,%g)\n",ipart,px[ipart],xyz[ipart][0],xyz[ipart][1],xyz[ipart][2]);
	}
	
	delete vizfile;
	delete [] px;
	delete [] xyz;
	return nparts;
}

#endif