#ifndef __ACTION_CC__
#define __ACTION_CC__

#include "b3d.h"

CB3D *CAction::b3d=NULL;
CAction::CAction(CB3D *b3dset){
	b3d=b3dset;
	currentmap=&(b3d->DeadActionMap);
	tau=0.0;
}

// type=0(creation) 1(decay) 2(collision)

CActionMap::iterator CAction::GetPos(CActionMap *emap){
	pair<CActionMap::iterator,CActionMap::iterator> epospair;
	CActionMap::iterator epos;
	epospair=emap->equal_range(key);
	epos=epospair.first;
	while(epos->second!=this && epos!=epospair.second){
		++epos;
	}
	if(epos->second!=this){
		return emap->end();
	}
	else return epos;
}

void CAction::SetKey(){
	key=tau;
}

void CAction::ChangeMap(CActionMap *newmap){
	DeleteFromCurrentMap();
	AddToMap(newmap);
}

CActionMap::iterator CAction::DeleteFromCurrentMap(){
	CActionMap::iterator neighbor,epos=GetPos(currentmap);
	neighbor=epos;
	neighbor++;
	if(epos->second!=this){
		printf("not getting correct actionptr\n");
		exit(1);
	}
	if(epos!=currentmap->end()){
		currentmap->erase(epos);
	}
	else{
		printf("In CAction::DeleteFromCurrenMap, can't find epos!!!\n");
		exit(1);
	}
	currentmap=NULL;
	return neighbor;
}

void CAction::AddToMap(CActionMap *newmap){
	SetKey();
	newmap->insert(CActionPair(key,this));
	currentmap=newmap;
}

void CAction::AddToMap(CActionMap::iterator guess,CActionMap *newmap){
	SetKey();
	newmap->insert(guess,CActionPair(key,this));
	currentmap=newmap;
}

void CAction::Kill(){
	CleanPartMap();
	DeleteFromCurrentMap();
	tau=0.0;
	AddToMap(b3d->DeadActionMap.begin(),&(b3d->DeadActionMap));
}

void CAction::CleanPartMap(){
	CPartMap::iterator ppos;
	CActionMap::iterator epos;
	CActionMap *emap;
	CPart *part;
	for(ppos=partmap.begin();ppos!=partmap.end();++ppos){
		part=ppos->second;
		emap=&(part->actionmap);
		epos=GetPos(emap);
		if(epos!=emap->end())	emap->erase(epos);
	}
	partmap.clear();
}

void CAction::AddPart(CPart *part){
	partmap.insert(CPartPair(part->key,part));
}

void CAction::Print(){
	printf("___________ type=%d, tau=%g, nparts=%d ___________\n",type,tau,int(partmap.size()));
	CPartMap::iterator ppos;
	CPart *part;
	for(ppos=partmap.begin();ppos!=partmap.end();++ppos){
		part=ppos->second;
		part->Print();
	}
	printf("_________________________________________________\n");
}

void CAction::CheckPartList(){
	CPart *part;
	CPartMap::iterator ppos,ppos2;
	CPartMap *pmap;
	ppos=partmap.begin();
	while(ppos!=partmap.end()){
		part=ppos->second;
		ppos2=part->GetPos(&(b3d->PartMap));
		if(ppos2==b3d->PartMap.end()){
			printf("____________ CAction::CheckPartList FATAL, action type=%d ________________\n",type);
			printf("iterator not in expected pmap\n");
			part->Print();
			//exit(1);
		}
		++ppos;
	}
}

#endif
