#ifndef __B3D_CC__
#define __B3D_CC__

#include "b3d.h"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif
using namespace std;

CB3D::CB3D(string run_name_set){
	run_name=run_name_set;
	string parsfilename,dirname;
	dirname="parameters/"+run_name;
	parsfilename=dirname+"/fixed.param";
	printf("reading %s\n",parsfilename.c_str());
	parameter::ReadParsFromFile(parmap,parsfilename);
	parsfilename=dirname+"/stats.param";
	printf("reading %s\n",parsfilename.c_str());
	parameter::ReadParsFromFile(parmap,parsfilename);
	NACTIONSMAX=parameter::getI(parmap,"B3D_NACTIONSMAX",5000);
	NPARTSMAX=parameter::getI(parmap,"B3D_NPARTSMAX",2000);
	TAUCOLLMAX=parameter::getD(parmap,"B3D_TAUCOLLMAX",50.0);
	VIZWRITE=parameter::getB(parmap,"B3D_VIZWRITE",false);
	input_dataroot=parameter::getS(parmap,"B3D_INPUT_DATAROOT","data/b3d");
	output_dataroot=parameter::getS(parmap,"B3D_OUTPUT_DATAROOT","data/b3d");
	//NxExVxExNxTxS=parameter::getI(parmap,"B3D_NACTIONS",1);
	NSAMPLE=parameter::getI(parmap,"B3D_NSAMPLE",1);
	ERROR_PRINT=parameter::getB(parmap,"B3D_ERROR_PRINT",true);
	XYMAX=parameter::getD(parmap,"B3D_XYMAX",15);
	ETAMAX=parameter::getD(parmap,"B3D_ETAMAX",1.0);
	NETA=parameter::getI(parmap,"B3D_NETA",10);
	NXY=parameter::getI(parmap,"B3D_NXY",10);
	SIGMAMAX=parameter::getD(parmap,"B3D_SIGMAMAX",30);
	NRINGSMAX=parameter::getI(parmap,"B3D_NRINGSMAX",200);
	COLLISIONS=parameter::getB(parmap,"B3D_COLLISIONS",true);
	BJORKEN=parameter::getB(parmap,"B3D_BJORKEN",false);
	SIGMADEFAULT=parameter::getD(parmap,"B3D_SIGMADEFAULT",1.5);

	SIGMAMAX=SIGMAMAX/double(NSAMPLE);
	NPARTSMAX*=NSAMPLE;
	NACTIONSMAX*=NSAMPLE;
	string command="mkdir -p output/"+run_name;
	system(command.c_str());
	double xmin,xmax,ymin,ymax,etamin,etamax;
	DXY=XYMAX/double(NXY);
	DETA=ETAMAX/double(NETA);
	int jx,jy,jeta;
	CCell *c;
	CCell::b3d=this;
	CPart::b3d=this;
	hydrotob3d=new CHYDROtoB3D();
	hydrotob3d->b3d=this;
	hydrotob3d->initialization=false;
	bjmaker.b3d=this;
	CResList::b3d=this;
	tau=0.0;
	itau=0;
	DeadPartMap.clear();
	PartMap.clear();
	FinalPartMap.clear();
	ActionMap.clear();
	DeadActionMap.clear();
	cell=new CCell***[2*NXY];
	int ix,iy,ieta;
	for(ix=0;ix<2*NXY;ix++){
		xmin=-XYMAX+ix*DXY;
		xmax=xmin+DXY;
		cell[ix]=new CCell**[2*NXY];
		for(iy=0;iy<2*NXY;iy++){
			ymin=-XYMAX+iy*DXY;
			ymax=ymin+DXY;
			cell[ix][iy]=new CCell*[2*NETA];
			for(ieta=0;ieta<2*NETA;ieta++){
				etamin=-ETAMAX+DETA*ieta;
				etamax=etamin+DETA;
				cell[ix][iy][ieta]=new CCell(xmin,xmax,ymin,ymax,etamin,etamax);
			}
		}
	}
	for(ix=0;ix<2*NXY;ix++){
		for(iy=0;iy<2*NXY;iy++){
			for(ieta=0;ieta<2*NETA;ieta++){
				c=cell[ix][iy][ieta];
				c->ix=ix; c->iy=iy; c->ieta=ieta;
				c->creflection=NULL;
				for(jx=ix-1;jx<=ix+1;jx++){
					for(jy=iy-1;jy<=iy+1;jy++){
						for(jeta=ieta-1;jeta<=ieta+1;jeta++){
							if(jx>=0 && jy>=0 && jeta>=0 && jx<2*NXY && jy<2*NXY && jeta<2*NETA){
								c->neighbor[1+jx-ix][1+jy-iy][1+jeta-ieta]=cell[jx][jy][jeta];
							}
							else c->neighbor[1+jx-ix][1+jy-iy][1+jeta-ieta]=NULL;
						}
					}
				}
			}
		}
	}
	if(BJORKEN){
		for(ix=0;ix<2*NXY;ix++){
			for(iy=0;iy<2*NXY;iy++){
				c=cell[ix][iy][0];
				c->ireflection=-1;
				c->creflection=cell[ix][iy][2*NETA-1];
				for(jx=ix-1;jx<=ix+1;jx++){
					for(jy=iy-1;jy<=iy+1;jy++){
						if(jx>=0 && jy>=0 && jx<2*NXY && jy<2*NXY){
							c->neighbor[1+jx-ix][1+jy-iy][0]=cell[jx][jy][2*NETA-1];
						}
					}
				}				
				c=cell[ix][iy][2*NETA-1];
				c->ireflection=1;
				c->creflection=cell[ix][iy][0];
				for(jx=ix-1;jx<=ix+1;jx++){
					for(jy=iy-1;jy<=iy+1;jy++){
						if(jx>=0 && jy>=0 && jx<2*NXY && jy<2*NXY){
							c->neighbor[1+jx-ix][1+jy-iy][2]=cell[jx][jy][0];
						}
					}
				}
			}
		}
	}
	randy=new CRandom(-1234);
	ievent_write=ievent_read=0;
// create particle objects, and put into deadpartmap
	int ipart;
	part=new CPart *[NPARTSMAX];
	for(ipart=0;ipart<NPARTSMAX;ipart++){
		part[ipart]=new CPart();
		part[ipart]->listid=ipart;
		part[ipart]->key=ipart;
		DeadPartMap.insert(CPartPair(part[ipart]->key,part[ipart]));
	}
// create array of action objects
	int iaction;
	action=new CAction *[NACTIONSMAX];
	for(iaction=0;iaction<NACTIONSMAX;iaction++){
		action[iaction]=new CAction(this);
		action[iaction]->listid=iaction;
		action[iaction]->SetKey();
		DeadActionMap.insert(CActionPair(action[iaction]->key,action[iaction]));
	}

//
	reslist=new CResList();

	h5infile=NULL;
	h5outfile=NULL;
	viz_file_id=-1;

	ptype=new CompType(sizeof(CPartH5));
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

void CB3D::SetQualifier(string qualifier_set){
	qualifier=qualifier_set;
	if(h5outfile!=NULL) delete h5outfile;
	if(h5infile!=NULL) delete h5infile;
	if(VIZWRITE && viz_file_id>0){
		H5Fflush(viz_file_id, H5F_SCOPE_LOCAL);//redundant with H5Fclose()
		H5Fclose(viz_file_id);
				//delete h5vizfile;
	}

	string outfilename="output/"+run_name+"/"+qualifier+"/b3d.h5";
	printf("will write to %s\n",outfilename.c_str());
	string command="mkdir -p output/"+run_name+"/"+qualifier;
	system(command.c_str());
	h5outfile = new H5File(outfilename,H5F_ACC_TRUNC);
	if(VIZWRITE){
		string vizfilename="output/"+run_name+"/"+qualifier+"/b3dviz.h5";
		printf("vizfilename=%s\n",vizfilename.c_str());
			//Retrieve current default error printing info.
			//herr_t (*err_stack_traverse_func)(void*);
		H5E_auto2_t err_stack_traverse_func;
		void *err_stack_client_data;
		H5Eget_auto2(H5E_DEFAULT, &err_stack_traverse_func, &err_stack_client_data);
			//Turn off automatic error printing
		if(H5Eset_auto2(H5E_DEFAULT, NULL, NULL)<0){
			printf("(-) Cannot turn off automatic HDF5 error printing.\n");
		}

			//----Create a new HDF5 file----
		hid_t pList_H5Faccess_id = H5P_DEFAULT;
		hid_t pList_H5Dwrite_id = H5P_DEFAULT;
		herr_t   status;
		viz_file_id = H5Fcreate(vizfilename.c_str(), H5F_ACC_DEBUG | H5F_ACC_TRUNC, H5P_DEFAULT, pList_H5Faccess_id);
		if(viz_file_id<0){
			printf("(X) Error in H5Fcreate(%s). Aborting...\n", vizfilename.c_str());
			exit(-1);
		}
			//----Create a new HDF5 file----

		//H5Fflush(viz_file_id, H5F_SCOPE_LOCAL);//redundant with H5Fclose()
		//H5Fclose(viz_file_id);
			//----Write to HDF5 format----

			//Restore the default automatic error trav
		//if(H5Eset_auto2(H5E_DEFAULT, err_stack_traverse_func, err_stack_client_data)<0){
		//	printf("(-) Cannot restore automatic HDF5 error printing. Abort...\n");
		//}
	}

	h5infile=NULL;
}

int CB3D::ReadDataH5(int ievent){
	KillAllActions();
	nactions=0;
	KillAllParts();
	if(h5infile==NULL){
		string infilename="output/"+run_name+"/"+qualifier+"/hydro.h5";
		h5infile = new H5File(infilename,H5F_ACC_RDONLY);
	}

	int nparts,ipart;
	char eventno[20];
	H5D_space_status_t status;
	sprintf(eventno,"event%d",ievent);
	hsize_t dim[1];
	DataSet *dataset = new DataSet (h5infile->openDataSet(eventno));
	//dataset->getSpaceStatus(status);
	//hsize_t datasetsize=dataset->getStorageSize();
	//printf("ievent=%d, status=%d, size=%d\n",ievent,int(status),int(datasetsize));
	DataSpace filespace = dataset->getSpace();
	int rank=filespace.getSimpleExtentDims(dim);
	nparts=dim[0];
	CPartH5 *partH5=new CPartH5[nparts];
	if(nparts>NPARTSMAX){
		printf("Increase NPARTSMAX, nparts=%d\n",nparts);
		exit(1);
	}
	dataset->read(partH5,*ptype);
	delete dataset;
	for(ipart=0;ipart<nparts;ipart++){
		part[ipart]->Init(partH5[ipart].ID,partH5[ipart].x,partH5[ipart].y,partH5[ipart].tau,partH5[ipart].eta,
			partH5[ipart].px,partH5[ipart].py,partH5[ipart].mass,partH5[ipart].rapidity);
	}
	delete [] partH5;
	//printf("READ IN %d PARTS\n",nparts);
	return nparts;
}

int CB3D::WriteDataH5(){
	double v,pperp,eperp,dnchdeta=0.0,t,tauwrite,eta;
	int dnchdy=0,ipart;
	CPart *part;
	CPartMap::iterator ppos;
	CCell *c;
	int ix,iy,ieta;
	for(ix=0;ix<2*NXY;ix++){
		for(iy=0;iy<2*NXY;iy++){
			for(ieta=0;ieta<2*NETA;ieta++){
				c=cell[ix][iy][ieta];
				while(c->partmap.size()>0){
					ppos=c->partmap.begin();
					part=ppos->second;
					part->DeleteFromCurrentMap();
					part->AddToMap(&FinalPartMap);
					//part->Print();
				}
			}
		}
	}

	ievent_write+=1;
	int nparts=int(FinalPartMap.size());
	CPartH5 *partH5=new CPartH5[nparts];
	hsize_t dim[] = {nparts};   /* Dataspace dimensions */
	DataSpace space(1,dim );

	ppos=FinalPartMap.begin();
	ipart=0;
	while(ppos!=FinalPartMap.end()){
		part=ppos->second;
		pperp=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
		eperp=sqrt(pperp*pperp+part->resinfo->mass*part->resinfo->mass);
		v=pperp/eperp;
		if(part->resinfo->charge!=0) dnchdy+=1.0;
		if(part->resinfo->charge!=0 && v>0.0 && v<1.0) dnchdeta+=v;
		tauwrite=part->tau_lastint;
		eta=part->GetEta(tauwrite);
		t=double(tauwrite*cosh(eta));
		partH5[ipart].listid=part->listid;
		partH5[ipart].ID=part->resinfo->code;
		partH5[ipart].x=part->r[1]+(part->p[1]/part->p[0])*(t-part->r[0]);
		partH5[ipart].y=part->r[2]+(part->p[2]/part->p[0])*(t-part->r[0]);
		partH5[ipart].tau=tauwrite;
		partH5[ipart].eta=eta;
		partH5[ipart].px=part->p[1];
		partH5[ipart].py=part->p[2];
		partH5[ipart].rapidity=part->y;
		partH5[ipart].mass=part->GetMass();
		++ppos;
		ipart+=1;
	}

/*
//CompType ptype( sizeof(CPartH5) );
//ptype.insertMember("listid", HOFFSET(CPartH5,listid), PredType::NATIVE_INT);
//ptype.insertMember("ID", HOFFSET(CPartH5,ID), PredType::NATIVE_INT);
//ptype.insertMember("px", HOFFSET(CPartH5,px), PredType::NATIVE_DOUBLE);
//ptype.insertMember("py", HOFFSET(CPartH5,py), PredType::NATIVE_DOUBLE);
//ptype.insertMember("rapidity", HOFFSET(CPartH5,rapidity), PredType::NATIVE_DOUBLE);
//ptype.insertMember("mass", HOFFSET(CPartH5,mass), PredType::NATIVE_DOUBLE);
//ptype.insertMember("x", HOFFSET(CPartH5,x), PredType::NATIVE_DOUBLE);
//ptype.insertMember("y", HOFFSET(CPartH5,y), PredType::NATIVE_DOUBLE);
//ptype.insertMember("eta", HOFFSET(CPartH5,eta), PredType::NATIVE_DOUBLE);
//ptype.insertMember("tau", HOFFSET(CPartH5,tau), PredType::NATIVE_DOUBLE);
*/

	DataSet* dataset;
	char event_number[20];
	sprintf(event_number,"event%d",ievent_write);
	printf("ievent_write=%d, writing data set %s\n",ievent_write,event_number);
	dataset = new DataSet(h5outfile->createDataSet(event_number,*ptype, space));
	dataset->write(partH5,*ptype);

	delete dataset;
	delete [] partH5;

	return dnchdy/(2.0*ETAMAX);
}

void CB3D::PerformAllActions(){
	CActionMap::iterator epos=ActionMap.begin();
	CAction *action;
	nscatter=ndecay=npass=nmerge=nswallow=npass=nexit=nactivate=ncheck=0;
	tau=0.0;
	nactions=0;

	while(epos!=ActionMap.end()){
		action=epos->second;
		action->Perform();
		epos=ActionMap.begin();
	}
}

void CB3D::KillAllActions(){
	CAction *action;
	CActionMap::iterator epos=ActionMap.begin();
	while(epos!=ActionMap.end()){
		action=epos->second;
		action->Kill();
		epos=ActionMap.begin();
	}
	nactions=0;
}

void CB3D::KillAllParts(){
	CPartMap *partmap=&PartMap;
	CPart *part;
	CPartMap::iterator ppos=partmap->begin();
	while(ppos!=partmap->end()){
		part=ppos->second;
		part->Kill();
		ppos=partmap->begin();
	}
	partmap=&FinalPartMap;
	ppos=partmap->begin();
	while(ppos!=partmap->end()){
		part=ppos->second;
		part->Kill();
		ppos=partmap->begin();
	}
	int ix,iy,ieta;
	for(ix=0;ix<2*NXY;ix++){
		for(iy=0;iy<2*NXY;iy++){
			for(ieta=0;ieta<2*NETA;ieta++){
				partmap=&(cell[ix][iy][ieta]->partmap);
				ppos=partmap->begin();
				while(ppos!=partmap->end()){
					part=ppos->second;
					part->Kill();
					ppos=partmap->begin();
				}
			}
		}
	}
}

void CB3D::AddAction_Activate(CPart *part){
	part->active=false;
	CAction *action;
	if(DeadActionMap.size()==0){
		printf("MUST INCREASE NACTIONS_MAX\n");
		exit(1);
	}
	if(BJORKEN && fabs(part->eta)>ETAMAX){
		printf("CB3D::AddAction_Activate, eta out of bounds, =%g\n",fabs(part->eta));
		exit(1);
	}
	CActionMap::iterator epos=DeadActionMap.begin();
	action=epos->second;
	action->DeleteFromCurrentMap();
	action->CleanPartMap(); // shouldn't be necessary, check to see if can be deleted
	action->type=0;
	action->tau=part->tau0;
	action->partmap.insert(CPartPair(part->key,part));
	action->AddToMap(&ActionMap);
	if(action->tau<tau){
		printf("trying to AddAction_Activate at earler time!!! action->tau=%g, tau=%g\n",action->tau,tau);
		exit(1);
	}
	part->actionmap.insert(CActionPair(action->key,action));
}

void CB3D::AddAction_VizWrite(double tauwrite){
	CAction *action;
	if(DeadActionMap.size()==0){
		printf("MUST INCREASE NACTIONS_MAX\n");
		exit(1);
	}
	CActionMap::iterator epos=DeadActionMap.begin();
	action=epos->second;
	action->DeleteFromCurrentMap();
	action->type=3;
	action->tau=tauwrite;
	action->AddToMap(&ActionMap);
	action->CleanPartMap();
	if(action->tau<tau){
		printf("trying to AddAction_VizWrite at earler time!!! action->tau=%g, tau=%g\n",action->tau,tau);
		exit(1);
	}
}

void CB3D::AddAction_Decay(CPart *part){
	CAction *action;
	double t,gamma,vz,newt,newz,taudecay;
	t=HBARC/part->resinfo->width;
	gamma=part->p[0]/part->GetMass();
	t=-t*gamma*log(randy->ran());
	vz=part->p[3]/part->p[0];
	newt=part->r[0]+t;
	newz=part->r[3]+vz*t;
	taudecay=sqrt(newt*newt-newz*newz);
	//printf("taudecay=%g, ID=%d, mR=%g\n",taudecay,part->resinfo->code,part->resinfo->mass);
	if(taudecay<part->tauexit || part->cell==NULL){
		if(part->tau0<tau-0.00001){
			printf("CB3D::AddAction_Decay, trying to add action for particle behind the times, tau=%g\n",tau);
			part->Print();
			exit(1);
		}
		if(DeadActionMap.size()==0){
			printf("MUST INCREASE NACTIONS MAX\n");
			exit(1);
		}
		CActionMap::iterator epos=DeadActionMap.begin();
		action=epos->second;
		action->DeleteFromCurrentMap();
		action->CleanPartMap();
		action->type=1;
		action->tau=taudecay;
		action->SetKey();
	//printf("added action at tau=%g, key=%lld\n",action->tau,action->key);
		action->partmap.insert(CPartPair(part->key,part));
		action->AddToMap(&ActionMap);
		part->actionmap.insert(CActionPair(action->key,action));
		if(action->tau<tau){
			printf("CB3D::AddAction_Decay, trying to AddAction_Decay at earler time!!! action->tau=%g, tau=%g\n",action->tau,tau);
			printf("t=%g, vz=%g\n",t,vz);
			printf("oldtau=%g, decaytau=%g\n",sqrt(pow(part->r[0],2)-pow(part->r[3],2)),sqrt(newt*newt-newz*newz));
			part->Print();
			part->cell->Print();
			exit(1);
		}
	}
}

void CB3D::AddAction_ExitCell(CPart *part){
	CAction *action;
	if(part->tauexit<TAUCOLLMAX){
		if(DeadActionMap.size()==0){
			printf("MUST INCREASE NACTIONS MAX\n");
			exit(1);
		}
		CActionMap::iterator epos=DeadActionMap.begin();
		action=epos->second;
		action->DeleteFromCurrentMap();
		action->CleanPartMap();
		action->type=6;
		action->tau=part->tauexit;
		action->SetKey();

	//printf("added action at tau=%g, key=%lld\n",action->tau,action->key);
		action->partmap.insert(CPartPair(part->key,part));
		action->AddToMap(&ActionMap);
		part->actionmap.insert(CActionPair(action->key,action));
		if(action->tau<tau){
			printf("CB3D::AddAction_ExitCell, trying to AddAction_ExitCell at earler time!!! action->tau=%g, tau=%g\n",action->tau,tau);
			part->Print();
			part->cell->Print();
			exit(1);
		}
	}
}

void CB3D::AddAction_Collision(CPart *part1,CPart *part2,double taucoll){
	CAction *action;
	if(DeadActionMap.size()==0){
		printf("MUST INCREASE NACTIONSMAX\n");
		exit(1);
	}
	CActionMap::iterator epos=DeadActionMap.begin();
	action=epos->second;
	action->DeleteFromCurrentMap();
	action->CleanPartMap(); // shouldn't be necessary, check to see if can be deleted
	action->type=2;
	action->tau=taucoll;
	action->SetKey();
	if(action->tau<tau){
		printf("trying to AddAction_Collision at earler time!!!  tau=%g\n",tau);
		action->Print();
		exit(1);
	}

	action->partmap.insert(CPartPair(part1->key,part1));
	action->partmap.insert(CPartPair(part2->key,part2));
	action->AddToMap(&ActionMap);

	part1->actionmap.insert(CActionPair(action->key,action));
	part2->actionmap.insert(CActionPair(action->key,action));
}

void CB3D::PrintActionMap(CActionMap *actionmap){
	CActionMap::iterator epos;
	CAction *action;
	int iaction=0;
	printf("_________________ ACTIONMAP %d actions _________________________\n",int(actionmap->size()));
	for(epos=actionmap->begin();epos!=actionmap->end();++epos){
		iaction+=1;
		action=epos->second;
		printf("iaction=%d : ",iaction);
		action->Print();
	}
}

void CB3D::FindAllCollisions(){
	double taucoll;
	CPartMap::iterator ppos1,ppos2;
	CPart *part1,*part2;
	CActionMap::iterator epos;
	CAction *action;
	int nbefore=ActionMap.size();
	//printf("CB3D::FindAllCollisions, Resetting Collisions\n");

	for(ppos1=PartMap.begin();ppos1!=PartMap.end();++ppos1){
		part1=ppos1->second;
		part1->KillActions();
	}

	for(epos=ActionMap.begin();epos!=ActionMap.end();++epos){
		action=epos->second;
		if(action->type==2){
			printf("CB3D::FindAllCollisions, expected all type-2 actions to be dead\n");
			exit(1);
		}
	}

	ppos1=PartMap.begin();
	part1=ppos1->second;
	ppos2=ppos1; ++ppos2;
	while(ppos2!=PartMap.end()){
		part1=ppos1->second; part2=ppos2->second;
		FindCollision(part1,part2,taucoll);
		ppos1=ppos2;
		++ppos2;
	}
	int nafter=ActionMap.size();
	if(nbefore!=nafter){
		printf("CB3D::FindAllCollisions, nbefore=%d, nafter=%d\n",nbefore,nafter);
		Misc::Pause();
	}

}

void CB3D::PrintPartList(){
	CPartMap::iterator ppos2,ppos1=FinalPartMap.begin();
	while(ppos1!=FinalPartMap.end()){
		printf("%d ",ppos1->second->listid);
		ppos2=ppos1; ++ppos2;
		if(ppos2!=FinalPartMap.end()){
			if(ppos1->second->actionmother!=ppos2->second->actionmother) printf("| ");
		}
		++ppos1;
	}
	printf("\n");
}

void CB3D::ListFutureCollisions(){
	CActionMap::iterator epos=ActionMap.begin();
	CAction *action;
	CPartMap::iterator p1,p2;
	printf("------------------- LIST OF FUTURE COLLISIONS ---------------------\n");
	while(epos!=ActionMap.end()){
		action=epos->second;
		if(action->type==2){
			p1=action->partmap.begin();
			p2=p1; ++p2;
			printf("%d  %d  will collide at %g\n",p1->second->listid,p2->second->listid,double(action->tau));
		}
		epos++;
	}
}

double CB3D::GetPiBsquared(CPart *part1,CPart *part2){
	int alpha;
	double r[4],P[4],q[4],Pdotq=0.0,Pdotr=0.0,P2=0.0,qdotr=0.0,q2=0.0,rsquared=0.0,g[4]={1,-1,-1,-1},y1,mt,eta1;
	double *p1=part1->p,*p2=part2->p,*r1=part1->r,*r2=part2->r;
	bool flip=false;
	if(BJORKEN && ((part1->cell->ieta==0 && part2->cell->ieta==2*NETA-1) || (part2->cell->ieta==0 && part1->cell->ieta==2*NETA-1))){
		flip=true;
		p1=new double[4];
		r1=new double[4];
		for(alpha=0;alpha<4;alpha++){
			p1[alpha]=part1->p[alpha];
			r1[alpha]=part1->r[alpha];
		}
		mt=sqrt(p1[0]*p1[0]-p1[3]*p1[3]);
		y1=part1->y;
		eta1=part1->eta;
		if(part1->cell->ieta==0){
			y1+=2.0*ETAMAX;
			eta1+=2.0*ETAMAX;
		}
		else{
			y1-=2.0*ETAMAX;
			eta1-=2.0*ETAMAX;
		}
		p1[0]=mt*cosh(y1);
		p1[3]=mt*sinh(y1);
		r1[0]=part1->tau0*cosh(eta1);
		r1[3]=part1->tau0*sinh(eta1);
		//printf("p1=(%g,%g,%g,%g), r1=(%g,%g,%g,%g)\n",p1[0],p1[1],p1[2],p1[3],r1[0],r1[1],r1[2],r1[3]);
		//printf("p2=(%g,%g,%g,%g), r2=(%g,%g,%g,%g)\n",p2[0],p2[1],p2[2],p2[3],r2[0],r2[1],r2[2],r2[3]);
	}
	for(alpha=0;alpha<4;alpha++){
		P[alpha]=p1[alpha]+p2[alpha];
		q[alpha]=p1[alpha]-p2[alpha];
		r[alpha]=r1[alpha]-r2[alpha];
		P2+=g[alpha]*P[alpha]*P[alpha];
		Pdotq+=g[alpha]*P[alpha]*q[alpha];
		Pdotr+=g[alpha]*P[alpha]*r[alpha];
	}
	for(alpha=0;alpha<4;alpha++){
		q[alpha]-=Pdotq*P[alpha]/P2;
		r[alpha]-=Pdotr*P[alpha]/P2;
		q2+=g[alpha]*q[alpha]*q[alpha];
		qdotr+=g[alpha]*q[alpha]*r[alpha];
	}
	for(alpha=0;alpha<4;alpha++){
		r[alpha]-=qdotr*q[alpha]/q2;
		rsquared-=g[alpha]*r[alpha]*r[alpha];
	}
	if(flip){
		//printf("pi*b^2=%g\n",PI*rsquared);
		delete [] p1;
		delete [] r1;
	}
	//else printf("pi*b^2=%g\n",PI*rsquared);
	return PI*rsquared;
}

void CB3D::freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt){
	const double prefactor=1.0/(2.0*PI*PI*pow(HBARC,3));
	double k0,k1,z,k0prime,k1prime,m2,m3,m4,t2,t3,I1,I2,Iomega;
	m2=m*m;
	m3=m2*m;
	m4=m2*m2;
	t2=t*t;
	t3=t2*t;
	z=m/t;
	if(z>1000.0){
		p=e=dens=dedt=0.0;
		printf("z is huge=%g, m=%g, t=%g\n",z,m,t);
	}
	else{
		if(z<0.0){
			printf("___z=%g,m=%g,T=%g ___\n",z,m,t);
			exit(1);
		}
		k0=Bessel::K0(z);
		k1=Bessel::K1(z);
		p=prefactor*(m2*t2*k0+2.0*m*t3*k1);
		e=prefactor*(3.0*m2*t2*k0+(m3*t+6.0*m*t3)*k1);
		dens=p/t;
		k0prime=-k1;
		k1prime=-k0-k1/z;
		dedt=prefactor*(6.0*m2*t*k0+(m3+18.0*m*t2)*k1-3.0*m3*k0prime-((m4/t)+6.0*m2*t)*k1prime);
		Iomega=exp(-m/t)/(30.0*PI*PI*HBARC*HBARC*HBARC);
		I1=pow(m,1.5)*pow(t,3.5)*7.5*sqrt(2.0*PI);
		I2=24.0*pow(t,5);
		sigma2=Iomega*(I1+I2+0.5*sqrt(I1*I2));  // this is an approximation (+/-2%) to messy integral
	}
}

void CB3D::Reset(){
	tau=0.0;
	KillAllActions();
	nactions=0;
	KillAllParts();
}

CB3D::~CB3D(){
	delete h5outfile;
	if(VIZWRITE){
		H5Fflush(viz_file_id, H5F_SCOPE_LOCAL);//redundant with H5Fclose()
		H5Fclose(viz_file_id);
		//delete h5vizfile;
	}
}

#endif
