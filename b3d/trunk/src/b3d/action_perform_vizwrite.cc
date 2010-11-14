#ifndef __B3DVIZ_CC__
#define __B3DVIZ_CC__

#include"b3d.h"

void CAction::PerformVizWrite(){
	double v,pperp,eperp,t,tauwrite=tau,eta;
	int ipart;
	CPart *part;
	CPartMap::iterator ppos;
	CCell *c;
	int ix,iy,ieta;
	int nparts,npartsmax=b3d->NPARTSMAX;
	double (*xyz)[3]=new double[npartsmax][3];
	int *listid=new int[npartsmax];
	int *ID=new int[npartsmax];
	double *px=new double[npartsmax];
	double *py=new double[npartsmax];
	double *rapidity=new double[npartsmax];
	double *mass=new double[npartsmax];
	
	nparts=0;
	for(ix=0;ix<2*b3d->NXY;ix++){
		for(iy=0;iy<2*b3d->NXY;iy++){
			for(ieta=0;ieta<2*b3d->NETA;ieta++){
				c=b3d->cell[ix][iy][ieta];
				ppos=c->partmap.begin();
				while(ppos!=c->partmap.end()){
					part=ppos->second;
					pperp=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
					eperp=sqrt(pperp*pperp+part->resinfo->mass*part->resinfo->mass);
					v=pperp/eperp;
					eta=part->GetEta(tauwrite);
					t=double(tauwrite*cosh(eta));
					listid[nparts]=part->listid;
					ID[nparts]=part->resinfo->code;
					xyz[nparts][0]=part->r[1]+(part->p[1]/part->p[0])*(t-part->r[0]);
					xyz[nparts][1]=part->r[2]+(part->p[2]/part->p[0])*(t-part->r[0]);
					xyz[nparts][2]=tauwrite*eta;
					px[nparts]=part->p[1];
					py[nparts]=part->p[2];
					rapidity[nparts]=part->y;
					mass[nparts]=part->GetMass();
					++ppos;
					nparts+=1;
				}
			}
		}
	}
	ppos=b3d->FinalPartMap.begin();
	while(ppos!=b3d->FinalPartMap.end()){
		part=ppos->second;
		pperp=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
		eperp=sqrt(pperp*pperp+part->resinfo->mass*part->resinfo->mass);
		v=pperp/eperp;
		eta=part->GetEta(tauwrite);
		t=double(tauwrite*cosh(eta));
		listid[nparts]=part->listid;
		ID[nparts]=part->resinfo->code;
		xyz[nparts][0]=part->r[1]+(part->p[1]/part->p[0])*(t-part->r[0]);
		xyz[nparts][1]=part->r[2]+(part->p[2]/part->p[0])*(t-part->r[0]);
		xyz[nparts][2]=tauwrite*eta;
		px[nparts]=part->p[1];
		py[nparts]=part->p[2];
		rapidity[nparts]=part->y;
		mass[nparts]=part->GetMass();
		++ppos;
		nparts+=1;
	}
	

	printf("Writing Viz Data for tau=%g, nparts=%d\n",tau,nparts);
	DataSet *dataset;
	char setname[40];
	
	sprintf(setname,"xyz_tau%g",tauwrite);
	hsize_t xyzdim[]={nparts,3};
	DataSpace xyzspace(2,xyzdim);
	dataset=new DataSet(b3d->h5vizfile->createDataSet(setname,PredType::NATIVE_DOUBLE,xyzspace));
	dataset->write(xyz,PredType::NATIVE_DOUBLE);
	delete dataset;
	
	sprintf(setname,"px_tau%g",tauwrite);
	hsize_t pdim[]={nparts};
	DataSpace pspace(1,pdim);
	dataset=new DataSet(b3d->h5vizfile->createDataSet(setname,PredType::NATIVE_DOUBLE,pspace));
	dataset->write(px,PredType::NATIVE_DOUBLE);
	delete dataset;
	
	sprintf(setname,"py_tau%g",tauwrite);
	dataset=new DataSet(b3d->h5vizfile->createDataSet(setname,PredType::NATIVE_DOUBLE,pspace));
	dataset->write(py,PredType::NATIVE_DOUBLE);
	delete dataset;
	
	sprintf(setname,"rapidity_tau%g",tauwrite);
	dataset=new DataSet(b3d->h5vizfile->createDataSet(setname,PredType::NATIVE_DOUBLE,pspace));
	dataset->write(rapidity,PredType::NATIVE_DOUBLE);
	delete dataset;
	
	sprintf(setname,"ID_tau%g",tauwrite);
	dataset=new DataSet(b3d->h5vizfile->createDataSet(setname,PredType::NATIVE_INT,pspace));
	dataset->write(ID,PredType::NATIVE_INT);
	delete dataset;
	
	delete [] xyz;
	delete [] px;
	delete [] py;
	delete [] ID;
	delete [] listid;
	delete [] mass;
	delete [] rapidity;
		
}

#endif
