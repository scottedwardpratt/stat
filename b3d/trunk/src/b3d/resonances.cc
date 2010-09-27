#ifndef __RESONANCES_CC__
#define __RESONANCES_CC__
#include "b3d.h"

using namespace std;

CB3D *CResList::b3d=NULL;

CResList::CResList(){
	parmap=&(b3d->parmap);
	ReadResInfo();
}

CRandom *CResInfo::ranptr=new CRandom(-1234);

CResInfo::CResInfo(){
	count=0;
	minmass=0.0;
	nextResInfoptr=NULL;
	firstbptr=NULL;
}

CMerge::CMerge(CResInfo *resinfo_in,double branching_in){
	resinfo=resinfo_in;
	branching=branching_in;
	next=NULL;
}

void CResInfo::DecayGetResInfoptr(int &nbodies,CResInfo **&daughterresinfoptr){
	double randy,bsum;
	int ibody;
	CBranchInfo *bptr;
	bptr=NULL;
	bsum=0.0;
	randy=ranptr->ran();

	do{
		if(bptr==NULL){
			bptr=firstbptr;
		}
		else{
			bptr=bptr->nextbptr;
		}
		bsum+=bptr->branching;
		if(bsum-1.0>1.0E-6){
			cout << "FATAL: In DecayGetResInfo: bsum too large, = " << bsum << endl;
			exit(1);
		}

	}while(bsum<randy);

	nbodies=bptr->nbodies;
	for(ibody=0;ibody<nbodies;ibody++){
		daughterresinfoptr[ibody]=bptr->resinfoptr[ibody];
	}

}
void CResInfo::DecayGetResInfoptr_minmass(int &nbodies,CResInfo **&daughterresinfoptr){
	nbodies=bptr_minmass->nbodies;
	for(int ibody=0;ibody<nbodies;ibody++){
		daughterresinfoptr[ibody]=bptr_minmass->resinfoptr[ibody];
	}

}

CBranchInfo::CBranchInfo(){
	nextbptr=NULL;
}

void CResList::ReadResInfo(){
	CMerge *merge;
	int mothercode,code,decay,strange,charge,baryon;
	double mass,mothermass,spin,width,bsum,netm;
	int ires,ires1,ires2,iresflip,ichannel,nchannels,ibody,nbodies,length;
	int netq,netb,nets;
	char dummy[120];
	string name;
	CResInfo *resinfoptr=NULL,*oldresinfoptr=NULL;
	CBranchInfo *bptr=NULL,*oldbptr=NULL;
	ifstream resinfofile,decayinfofile;
	stringstream sst;
	string filename;
	filename=parameter::getS(*parmap,"B3D_RESONANCES_INFO_FILE",string("resinfo/resonances_standardhadrons.dat"));
	cout << "will read res info from " << filename << endl;
	resinfofile.open(filename.c_str());
	//resinfofile.getline(dummy,100);
	resinfofile >> NResonances;
	printf("NResonances=%d\n",NResonances);
	MergeArray=new CMerge **[NResonances];
	for(ires=0;ires<NResonances;ires++){
		MergeArray[ires]=new CMerge *[NResonances];
		for(ires2=0;ires2<NResonances;ires2++){
			MergeArray[ires][ires2]=NULL;
		}
	}

	for(ires=0;ires<NResonances;ires++){
		resinfofile >> code >> mass >> charge >> baryon >> strange >> spin 
			>> decay >> width >> name;
		resinfofile.getline(dummy,100);
		resinfoptr=new CResInfo();
		if(oldresinfoptr==NULL){
			GfirstResInfoptr=resinfoptr;
		}
		else{
			oldresinfoptr->nextResInfoptr=resinfoptr;
		}
		resinfoptr->ires=ires;
		resinfoptr->code=code;
		resinfoptr->mass=mass;
		resinfoptr->charge=charge;
		resinfoptr->baryon=baryon;
		resinfoptr->strange=strange;
		resinfoptr->spin=spin;
		resinfoptr->decay=bool(decay);
		resinfoptr->width=width;    
		resinfoptr->name=name;
		/*
		cout << "name=" << name << ":\n";
		cout << "code=" << resinfoptr->code << ", M=" << resinfoptr->mass
			 << ", Q=" << resinfoptr->charge << ", B=" << resinfoptr->baryon
		 << ", S=" << resinfoptr->strange << ", s=" << resinfoptr->spin
		 << ", decay=" << resinfoptr->decay 
		 << ", W=" << resinfoptr->width << endl;*/

		oldresinfoptr=resinfoptr;
	}
	resinfofile.close();

	//cout << "RESONANCES READ, WILL BEGIN READING DECAY INFO\n";

	filename=parameter::getS(*parmap,"B3D_RESONANCES_DECAYS_FILE",string("resinfo/decays_weak.dat"));
	cout << "will read decay info from " << filename << endl;
	decayinfofile.open(filename.c_str());

	while(decayinfofile >> mothercode >> mothermass){
		decayinfofile.getline(dummy,80);
		decayinfofile >> mothercode >> nchannels;
		GetResInfoptr(mothercode,resinfoptr);
		resinfoptr->minmass=1.0E10;
		bsum=0.0;

		for(ichannel=0;ichannel<nchannels;ichannel++){
			bptr=new CBranchInfo();
			if(ichannel==0){
				resinfoptr->firstbptr=bptr;
			}
			else{
				oldbptr->nextbptr=bptr;
			}
			decayinfofile >> bptr->nbodies;
			netq=-resinfoptr->charge;
			netb=-resinfoptr->baryon;
			nets=-resinfoptr->strange;
			netm=0.0;

			for(ibody=0;ibody<bptr->nbodies;ibody++){
				decayinfofile >> code;
				GetResInfoptr(code,bptr->resinfoptr[ibody]);
				netq+=bptr->resinfoptr[ibody]->charge;
				netb+=bptr->resinfoptr[ibody]->baryon;
				nets+=bptr->resinfoptr[ibody]->strange;
				netm+=bptr->resinfoptr[ibody]->mass;
			}
			if(netm<resinfoptr->minmass){
				resinfoptr->minmass=netm;
				resinfoptr->bptr_minmass=bptr;
			}
			if(netq!=0 || netb!=0 || abs(nets)>1){
				cout << "Charge conservation failure while reading decay info,\nnetq=" << netq 
					<< ", netb=" << netb <<", nets=" << nets <<endl;
				cout << "mother=" << resinfoptr->name 
					<< "(" << resinfoptr->code << "), "
					<< "strangeness=" << resinfoptr->strange << ", ichannel="
					<< ichannel << endl;
				cout << "DAUGHTERS: ";
				for(ibody=0;ibody<bptr->nbodies;ibody++){
					cout << bptr->resinfoptr[ibody]->name << ", s=" 
						<< bptr->resinfoptr[ibody]->strange << endl;
				}
				cout << endl;
				if(netq!=0 || netb!=0) exit(1);
			}
			decayinfofile >> bptr->branching;
			if(bptr->nbodies==2){
				ires1=bptr->resinfoptr[0]->ires;
				ires2=bptr->resinfoptr[1]->ires;
				if(ires1>ires2){
					iresflip=ires1; ires1=ires2; ires2=iresflip;
				}
				merge=MergeArray[ires1][ires2];
				if(merge==NULL){
					MergeArray[ires1][ires2]=new CMerge(resinfoptr,bptr->branching);
				}
				else{
					while(merge->next!=NULL){
						merge=merge->next;
					}
					merge->next=new CMerge(resinfoptr,bptr->branching);
				}
			}
			
			bsum+=bptr->branching;
			oldbptr=bptr;
		}
	}
	decayinfofile.close();
}

void CResList::GetResInfoptr(int code,CResInfo *&resinfoptr){
	resinfoptr=GfirstResInfoptr;
	while(resinfoptr->code!=code){
		resinfoptr=resinfoptr->nextResInfoptr;
		if(resinfoptr==NULL){
			cout << "FATAL: failure to identify code=" << code << endl;
			exit(1);
		}
	}
}

void CResList::PrintYields(){
	CResInfo *resinfoptr;
	cout << "_____________________\n     YIELDS\n";
	resinfoptr=GfirstResInfoptr;
	while(resinfoptr!=NULL){
		cout << resinfoptr->name << ": " << resinfoptr->count << endl;
		resinfoptr=resinfoptr->nextResInfoptr;
	}
}

#endif
