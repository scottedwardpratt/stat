#ifndef __ANALYZE_BALANCE_CC__
#define __ANALYZE_BALANCE_CC__
#include "analyze.h"
// Note, works only for pions-only with CBjMaker::balance=true
void CAnalyze::CalcBalance(){
        long long int naccept=0,dnaccept=0,npairs=0,NBphi=0,NBy=0;
        const int Ny=20,Nphi=24;
        long long int Bphi[Nphi]={0},By[Ny]={0};
        int ievent=1,ipart,nparts,nevents,iy,iphi;
        double delphi=PI/double(Nphi),dely=2.0/double(Ny),dphi;
        double yi,yj,phii,phij;
        CPartH5 *parti,*partj;
        string infilename=input_dataroot+"/"+qualifier+"/"+h5_infilename;
        printf("CAnalyze::ReadData, opening %s\n",infilename.c_str());
        h5file = new H5File(infilename,H5F_ACC_RDONLY);
        nevents=int(h5file->getNumObjs());
        if(nevents>neventsmax) nevents=neventsmax;
        for(ievent=1;ievent<=nevents;ievent++){
                nparts=ReadDataH5(ievent);
                //printf("nparts=%d\n",nparts);
                for(ipart=0;ipart<nparts;ipart+=2){
                        parti=&partH5[ipart];
                        yi=parti->rapidity;
                        phii=atan2(parti->py,parti->px);
                        //printf("yi=%g, ID=%d\n",parti->rapidity,parti->ID);
                        if(abs(parti->ID)==211){
                                NBphi+=1; NBy+=1;
                                partj=&partH5[ipart+1];
                                if(partj->ID!=-parti->ID){
                                        printf("IDs not right, parti->ID=%d, partj->ID=%d\n",parti->ID,partj->ID);
                                        exit(1);
                                }
                                yj=partj->rapidity;
                                iy=lrint(floor(fabs((yi-yj)/dely)));
                                if(iy<Ny){
                                        By[iy]+=1;
                                }
                                phij=atan2(partj->py,partj->px);
                                dphi=fabs(phii-phij);
                                if(dphi>PI) dphi=2.0*PI-dphi;
                                //printf("phii=%g, phij=%g, dphi=%g\n",phii*180.0/PI,phij*180.0/PI,dphi*180.0/PI);
                                iphi=lrint(floor(dphi/delphi));
                                Bphi[iphi]+=1;
                                if(iphi>Nphi){
                                        printf("iphi=%d, out of range\n",iphi);
                                        exit(1);
                                }
                        }
                }
        }
        delete h5file;
        printf("------------- B(y1-y2) ---------------------\n");
        for(iy=0;iy<Ny;iy++){
                printf("%5.3f %g\n",(iy+0.5)*dely,double(By[iy])/double(NBy));
        }
        printf("----------- B(phi1-phi2) -------------------\n");
        for(iphi=0;iphi<Nphi;iphi++){
                printf("%5.3f %g\n",(iphi+0.5)*delphi*180.0/PI,double(Bphi[iphi])/double(NBphi));
        }
}
#endif
