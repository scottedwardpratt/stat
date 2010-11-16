#ifndef __B3DVIZ_CC__
#define __B3DVIZ_CC__

#include"b3d.h"

#ifdef WIN32
#define snprintf sprintf_s
#endif//WIN32

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
	hid_t file_id=b3d->viz_file_id;
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


	//----Write to HDF5 format----
	int eventIdx = b3d->ievent_write+1;
	//char caseString[256];
	char groupString[256];
	//----Create Case level group----
	hid_t group_case_id;
	if((group_case_id = H5Gopen2(file_id, b3d->qualifier.c_str(), H5P_DEFAULT))<0){
		group_case_id = H5Gcreate(file_id, b3d->qualifier.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		//---- Case level attributes ----
		//int attr_int;
		//float attr_float;
		hid_t aid, atype, attr;
		herr_t err;

		//1. code_name - Code Name (string)
		char codeName[256] = "Scott Pratt B3H";
		aid	= H5Screate(H5S_SCALAR);
		atype = H5Tcopy(H5T_C_S1);
		err	= H5Tset_size(atype, strlen(codeName));
		attr	= H5Acreate(group_case_id, codeName, atype, aid, H5P_DEFAULT, H5P_DEFAULT);
		err	= H5Awrite(attr, atype, "Scott Pratt b3h"); 
		err	= H5Aclose(attr);
		err	= H5Sclose(aid);

		////2. version - Version (float)
		//attr_float = qgp_case->GetVersion();
		//aid	= H5Screate(H5S_SCALAR);
		//attr	= H5Acreate(group_case_id, "Version", H5T_NATIVE_FLOAT, aid, H5P_DEFAULT, H5P_DEFAULT);
		//err	= H5Awrite(attr, H5T_NATIVE_FLOAT, &attr_float); 
		//err	= H5Aclose(attr);
		//err	= H5Sclose(aid);

		////3. aproj - ??? (int)
		//attr_int = qgp_case->GetAproj();
		//aid	= H5Screate(H5S_SCALAR);
		//attr	= H5Acreate(group_case_id, "aproj", H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
		//err	= H5Awrite(attr, H5T_NATIVE_FLOAT, &attr_int); 
		//err	= H5Aclose(attr);
		//err	= H5Sclose(aid);

		////4. zproj - ??? (int)
		//attr_int = qgp_case->GetZproj();
		//aid	= H5Screate(H5S_SCALAR);
		//attr	= H5Acreate(group_case_id, "zproj", H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
		//err	= H5Awrite(attr, H5T_NATIVE_FLOAT, &attr_int); 
		//err	= H5Aclose(attr);
		//err	= H5Sclose(aid);

		////5. atarg - ??? (int)
		//attr_int = qgp_case->GetAtarg();
		//aid	= H5Screate(H5S_SCALAR);
		//attr	= H5Acreate(group_case_id, "atarg", H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
		//err	= H5Awrite(attr, H5T_NATIVE_FLOAT, &attr_int); 
		//err	= H5Aclose(attr);
		//err	= H5Sclose(aid);

		////6. ztarg - ??? (int)
		//attr_int = qgp_case->GetZtarg();
		//aid	= H5Screate(H5S_SCALAR);
		//attr	= H5Acreate(group_case_id, "ztarg", H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
		//err	= H5Awrite(attr, H5T_NATIVE_FLOAT, &attr_int); 
		//err	= H5Aclose(attr);
		//err	= H5Sclose(aid);

		////7. refframe - Reference Frame (string)
		//aid	= H5Screate(H5S_SCALAR);
		//atype = H5Tcopy(H5T_C_S1);
		//err	= H5Tset_size(atype, qgp_case->SizeofRefFrame());
		//attr	= H5Acreate(group_case_id, "Reference Frame", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
		//err	= H5Awrite(attr, atype, qgp_case->GetRefFrame()); 
		//err	= H5Aclose(attr);
		//err	= H5Sclose(aid);

		////8. ebeam - ??? (float)
		//attr_float = qgp_case->GetEbeam();
		//aid	= H5Screate(H5S_SCALAR);
		//attr	= H5Acreate(group_case_id, "ebeam", H5T_NATIVE_FLOAT, aid, H5P_DEFAULT, H5P_DEFAULT);
		//err	= H5Awrite(attr, H5T_NATIVE_FLOAT, &attr_float); 
		//err	= H5Aclose(attr);
		//err	= H5Sclose(aid);

		////9. ntestpart - ??? (int)
		//attr_int = qgp_case->GetNtestpart();
		//aid	= H5Screate(H5S_SCALAR);
		//attr	= H5Acreate(group_case_id, "ntestpart", H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
		//err	= H5Awrite(attr, H5T_NATIVE_FLOAT, &attr_int); 
		//err	= H5Aclose(attr);
		//err	= H5Sclose(aid);

		////Add comments
		//err	= H5Gset_comment(group_case_id, ".", "Unit: momentum(GeV/c), energy(GeV), mass(GeV/(c^2))\
		//														 position(fm), time(fm/c), angle(deg), radius(fm), charge(GeV), pressure(???)");	
		//---- Case level attributes ----	
	}
	//----Create Case level group----

	//----Create Event level group----
	hid_t group_event_id;
	hid_t f;
	herr_t status;
	snprintf(groupString, sizeof(groupString), "Event %d", eventIdx);
	if((group_event_id = H5Gopen(group_case_id, groupString, H5P_DEFAULT))<0){
		group_event_id = H5Gcreate(group_case_id,  groupString, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		//---- Event level attributes ----
		hid_t aid, attr;
		herr_t err;

		//1. Event Index
		aid	= H5Screate(H5S_SCALAR);
		attr	= H5Acreate(group_event_id, "Index", H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
		err	= H5Awrite(attr, H5T_NATIVE_INT, &eventIdx); 
		err	= H5Aclose(attr);
		err	= H5Sclose(aid);

		////2. bimp - impact parameter
		//aid	= H5Screate(H5S_SCALAR);
		//attr	= H5Acreate(group_event_id, "Impact Parameter", H5T_NATIVE_FLOAT, aid, H5P_DEFAULT, H5P_DEFAULT);
		//err	= H5Awrite(attr, H5T_NATIVE_FLOAT, &fas[i].bimp); 
		//err	= H5Aclose(attr);
		//err	= H5Sclose(aid);

		////3. phi - azimuthal angle
		//aid	= H5Screate(H5S_SCALAR);
		//attr	= H5Acreate(group_event_id, "Azimuthal Angle", H5T_NATIVE_FLOAT, aid, H5P_DEFAULT, H5P_DEFAULT);
		//err	= H5Awrite(attr, H5T_NATIVE_FLOAT, &fas[i].phi); 
		//err	= H5Aclose(attr);
		//err	= H5Sclose(aid);

		////4. Number of Frames
		//aid	= H5Screate(H5S_SCALAR);
		//attr	= H5Acreate(group_event_id, "Number of Frames", H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
		//err	= H5Awrite(attr, H5T_NATIVE_INT, &fas[i].numFrames); 
		//err	= H5Aclose(attr);
		//err	= H5Sclose(aid);
		//---- Event level attributes ----
	}
	//----Create Event level group----

	//----Create Frame level group----
	hid_t group_frame_id;
	snprintf(groupString, sizeof(groupString), "frame, tau=%g", tauwrite);	
	if((group_frame_id = H5Gopen(group_event_id, groupString, H5P_DEFAULT))<0){
		group_frame_id = H5Gcreate(group_event_id, groupString, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		//---- Frame level attributes ----
		hid_t aid, atype, attr;
		herr_t err;

		//1. Index
		aid	= H5Screate(H5S_SCALAR);
		attr	= H5Acreate(group_frame_id, "Index", H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
		err	= H5Awrite(attr, H5T_NATIVE_INT, &f); 
		err	= H5Aclose(attr);
		err	= H5Sclose(aid);

		//2. Number of Particles
		aid	= H5Screate(H5S_SCALAR);
		attr	= H5Acreate(group_frame_id, "Number of Particles", H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
		err	= H5Awrite(attr, H5T_NATIVE_INT, &nparts); 
		err	= H5Aclose(attr);
		err	= H5Sclose(aid);
		//---- Frame level attributes ----

		hsize_t dim[2];
		dim[0] = nparts;

		//----Write out ID----
		{
			dim[1] = 1;
			hid_t space_id = H5Screate_simple(1, dim, NULL);
			hid_t dataset_id = H5Dcreate(group_frame_id, "ID", H5T_NATIVE_INT, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			if(dataset_id < 0){
				printf("(X) Error in H5Dcreate(ID). Abort...\n");
				exit(-1);
			}
			//----attribute----
			aid	= H5Screate(H5S_SCALAR);
			atype = H5Tcopy(H5T_C_S1);
			err	= H5Tset_size(atype, sizeof(int));
			attr	= H5Acreate(dataset_id, "Unit", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
			err	= H5Awrite(attr, atype, "none"); 
			err	= H5Aclose(attr);
			err	= H5Sclose(aid);
			//----attribute----

			status = H5Dwrite(dataset_id, H5T_NATIVE_INT, space_id, H5S_ALL, H5P_DEFAULT, ID);
			if(status < 0){
				printf("(X) Error in H5Dwrite(ID). Abort...\n");
				exit(-1);
			} 

			H5Dclose(dataset_id);
			H5Sclose(space_id);
		}
		//----Write out ID----

		//----Write out mass----
		{
			dim[1] = 1;
			hid_t space_id = H5Screate_simple(1, dim, NULL);
			hid_t dataset_id = H5Dcreate(group_frame_id, "mass", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			if(dataset_id < 0){
				printf("(X) Error in H5Dcreate(mass). Abort...\n");
				exit(-1);
			}
			//----attribute----
			aid	= H5Screate(H5S_SCALAR);
			atype = H5Tcopy(H5T_C_S1);
			err	= H5Tset_size(atype, sizeof(double));
			attr	= H5Acreate(dataset_id, "Unit", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
			err	= H5Awrite(attr, atype, "GeV/(c^2)"); 
			err	= H5Aclose(attr);
			err	= H5Sclose(aid);
			//----attribute----

			status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, space_id, H5S_ALL, H5P_DEFAULT, mass);
			if(status < 0){
				printf("(X) Error in H5Dwrite(mass). Abort...\n");
				exit(-1);
			} 

			H5Dclose(dataset_id);
			H5Sclose(space_id);
		}
		//----Write out mass----

		//----Write out px----
		{
			dim[1] = 1;
			hid_t space_id = H5Screate_simple(1, dim, NULL);
			hid_t dataset_id = H5Dcreate(group_frame_id, "px", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			if(dataset_id < 0){
				printf("(X) Error in H5Dcreate(px). Abort...\n");
				exit(-1);
			}
			//----attribute----
			aid	= H5Screate(H5S_SCALAR);
			atype = H5Tcopy(H5T_C_S1);
			err	= H5Tset_size(atype, sizeof(double));
			attr	= H5Acreate(dataset_id, "Unit", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
			err	= H5Awrite(attr, atype, "GeV/(fm^3)"); 
			err	= H5Aclose(attr);
			err	= H5Sclose(aid);
			//----attribute----

			status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, space_id, H5S_ALL, H5P_DEFAULT, px);
			if(status < 0){
				printf("(X) Error in H5Dwrite(px). Abort...\n");
				exit(-1);
			} 

			H5Dclose(dataset_id);
			H5Sclose(space_id);
		}
		//----Write out px----

		//----Write out py----
		{
			dim[1] = 1;
			hid_t space_id = H5Screate_simple(1, dim, NULL);
			hid_t dataset_id = H5Dcreate(group_frame_id, "py", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			if(dataset_id < 0){
				printf("(X) Error in H5Dcreate(py). Abort...\n");
				exit(-1);
			}
			//----attribute----
			aid	= H5Screate(H5S_SCALAR);
			atype = H5Tcopy(H5T_C_S1);
			err	= H5Tset_size(atype, sizeof(double));
			attr	= H5Acreate(dataset_id, "Unit", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
			err	= H5Awrite(attr, atype, "GeV/(fm^3)"); 
			err	= H5Aclose(attr);
			err	= H5Sclose(aid);
			//----attribute----

			status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, space_id, H5S_ALL, H5P_DEFAULT, py);
			if(status < 0){
				printf("(X) Error in H5Dwrite(py). Abort...\n");
				exit(-1);
			} 

			H5Dclose(dataset_id);
			H5Sclose(space_id);
		}
		//----Write out py----

		//----Write out rapidity----
		{
			dim[1] = 1;
			hid_t space_id = H5Screate_simple(1, dim, NULL);
			hid_t dataset_id = H5Dcreate(group_frame_id, "rapidity", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			if(dataset_id < 0){
				printf("(X) Error in H5Dcreate(rapidity). Abort...\n");
				exit(-1);
			}
			//----attribute----
			aid	= H5Screate(H5S_SCALAR);
			atype = H5Tcopy(H5T_C_S1);
			err	= H5Tset_size(atype, sizeof(double));
			attr	= H5Acreate(dataset_id, "Unit", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
			err	= H5Awrite(attr, atype, "GeV/(fm^3)"); 
			err	= H5Aclose(attr);
			err	= H5Sclose(aid);
			//----attribute----

			status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, space_id, H5S_ALL, H5P_DEFAULT, rapidity);
			if(status < 0){
				printf("(X) Error in H5Dwrite(rapidity). Abort...\n");
				exit(-1);
			} 

			H5Dclose(dataset_id);
			H5Sclose(space_id);
		}
		//----Write out rapidity----

		//----Write out xyz----
		{
			dim[1] = 3;
			hid_t space_id = H5Screate_simple(2, dim, NULL);
			hid_t dataset_id = H5Dcreate(group_frame_id, "xyz", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			if(dataset_id < 0){
				printf("(X) Error in H5Dcreate(xyz). Abort...\n");
				exit(-1);
			}
			//----attribute----
			aid	= H5Screate(H5S_SCALAR);
			atype = H5Tcopy(H5T_C_S1);
			err	= H5Tset_size(atype, sizeof(double));
			attr	= H5Acreate(dataset_id, "Unit", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
			err	= H5Awrite(attr, atype, "fm"); 
			err	= H5Aclose(attr);
			err	= H5Sclose(aid);
			//----attribute----

			status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, space_id, H5S_ALL, H5P_DEFAULT, xyz);
			if(status < 0){
				printf("(X) Error in H5Dwrite(xyz). Abort...\n");
				exit(-1);
			} 

			H5Dclose(dataset_id);
			H5Sclose(space_id);
		}
		//----Write out xyz----
	}
	H5Gclose(group_frame_id);
	//----Create Frame level group----

	H5Gclose(group_event_id);
	//----Create Event level group----

	H5Gclose(group_case_id);
	//----Create Case level group----
	//----Write to HDF5 format----

	delete [] xyz;
	delete [] px;
	delete [] py;
	delete [] ID;
	delete [] listid;
	delete [] mass;
	delete [] rapidity;	
}

#endif//__B3DVIZ_CC__

