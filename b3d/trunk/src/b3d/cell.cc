#ifndef __PART_CC__
#define __PART_CC__

#include "b3d.h"
using namespace std;

CB3D *CCell::b3d=NULL;

CCell::CCell(double xminset,double xmaxset,double yminset,double ymaxset,double etaminset,double etamaxset){
	xmin=xminset; xmax=xmaxset; ymin=yminset; ymax=ymaxset; etamin=etaminset; etamax=etamaxset;	
	ireflection=0;
	creflection=NULL;
}

void CCell::Print(){
	printf("___ CELL INFO _____\n");
	printf("ix=%d, iy=%d, ieta=%d, xmin=%g, xmax=%g, ymin=%g, ymax=%g, etamin=%g, etamax=%g", ix,iy,ieta,xmin,xmax,ymin,ymax,etamin,etamax);
	printf("---------------------\n");
}

#endif
