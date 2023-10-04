#include "msu_commonutils/arrays.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/log.h"

using namespace std;

Crandy *C3DArray::randy=NULL;

C3DArray::C3DArray(string arrayparsfilename){
	// Default Pars
	NXMAX=NYMAX=NZMAX=20;
	DELX=DELY=DELZ=1.0;
	XSYM=YSYM=ZSYM=true;
	ReadPars(arrayparsfilename);
	CreateArray();
}

void C3DArray::ReadPars(string arrayparsfilename){
	CparameterMap apars;
	apars.ReadParsFromFile(arrayparsfilename);

	if(apars.getB("IDENTICAL",true)){
		apars.set("XSYM",true);
		apars.set("YSYM",true);
		apars.set("ZSYM",true);
	}
	if(apars.getD("DELXYZ",0) > 1.0E-10){
		DELX=DELY=DELZ=apars.getD("DELXYZ",0);
	}

	NXMAX=apars.getI("NXMAX",NXMAX);
	NYMAX=apars.getI("NYMAX",NYMAX);
	NZMAX=apars.getI("NZMAX",NZMAX);
	DELX=apars.getD("DELX",DELX);
	DELY=apars.getD("DELY",DELY);
	DELZ=apars.getD("DELZ",DELZ);
	XSYM=bool(apars.getB("XSYM",XSYM));
	YSYM=bool(apars.getB("YSYM",YSYM));
	ZSYM=bool(apars.getB("ZSYM",ZSYM));
	if(apars.getI("NXYZMAX",0) != 0){
		NXMAX=NYMAX=NZMAX=apars.getI("NXYZMAX",25);
	}
}

C3DArray::C3DArray(int NXMAXset,double DELXset,int NYMAXset,double DELYset,int NZMAXset,double DELZset,bool XSYMset,bool YSYMset,bool ZSYMset){
	XSYM=XSYMset;
	YSYM=YSYMset;
	ZSYM=ZSYMset;
	NXMAX=NXMAXset;
	NYMAX=NYMAXset;
	NZMAX=NZMAXset;
	DELX=DELXset;
	DELY=DELYset;
	DELZ=DELZset;

	CreateArray();
}

C3DArray::C3DArray(int NXYZMAXset,double DELXYZset,bool XSYMset,bool YSYMset,bool ZSYMset){
	XSYM=XSYMset;
	YSYM=YSYMset;
	ZSYM=ZSYMset;
	NXMAX=NXYZMAXset;
	NYMAX=NXYZMAXset;
	NZMAX=NXYZMAXset;
	DELX=DELXYZset;
	DELY=DELXYZset;
	DELZ=DELXYZset;

	CreateArray();
}

void C3DArray::CreateArray(){
	int ix,iy,iz,isx,isy,isz;
	int nsx,nsy,nsz;
	nsx=nsy=nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	F=new double *****[nsx];
	for(isx=0;isx<nsx;isx++){
		F[isx]=new double ****[NXMAX];
		for(ix=0;ix<NXMAX;ix++){
			F[isx][ix]=new double ***[nsy];
			for(isy=0;isy<nsy;isy++){
				F[isx][ix][isy]=new double **[NYMAX];
				for(iy=0;iy<NYMAX;iy++){
					F[isx][ix][isy][iy]=new double *[nsz];
					for(isz=0;isz<nsz;isz++){
						F[isx][ix][isy][iy][isz]=new double[NZMAX];
						for(iz=0;iz<NZMAX;iz++)
							F[isx][ix][isy][iy][isz][iz]=0.0;
					}
				}
			}
		}
	}
}

C3DArray::~C3DArray(){
	DeleteArray();
}

void C3DArray::DeleteArray(){
	int isx,isy,isz,ix,iy,nsx,nsy,nsz;
	nsx=nsy=nsz=1;
	if(!XSYM) nsx=2;
	if(!YSYM) nsy=2;
	if(!ZSYM) nsz=2;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<NXMAX;ix++){
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<NYMAX;iy++){
					for(isz=0;isz<nsz;isz++){
						delete [] F[isx][ix][isy][iy][isz];
					}
					delete [] F[isx][ix][isy][iy];
				}
				delete [] F[isx][ix][isy];
			}
			delete [] F[isx][ix];
		}
		delete [] F[isx];
	}
	delete [] F;    
}

bool C3DArray::GetXSYM(){
	return XSYM;
}

bool C3DArray::GetYSYM(){
	return YSYM;
}

bool C3DArray::GetZSYM(){
	return ZSYM;
}

void C3DArray::ZeroArray(){
	int ix,iy,iz,isx,isy,isz,nsx,nsy,nsz;
	nsx=nsy=nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<NXMAX;ix++){
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<NYMAX;iy++){
					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<NZMAX;iz++){
							F[isx][ix][isy][iy][isz][iz]=0.0;
						}
					}
				}
			}
		}
	}
}

void C3DArray::ScaleArray(double scalefactor){
	int ix,iy,iz,isx,isy,isz,nsx,nsy,nsz;
	nsx=nsy=nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<NXMAX;ix++){
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<NYMAX;iy++){
					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<NZMAX;iz++)
							F[isx][ix][isy][iy][isz][iz]*=scalefactor;
					}
				}
			}
		}
	}
}

double C3DArray::GetElement(int isx,int ix,int isy,int iy,int isz,int iz){
	return F[isx][ix][isy][iy][isz][iz];
}

void C3DArray::SetElement(int isx,int ix,int isy,int iy,int isz,int iz,double value){
	F[isx][ix][isy][iy][isz][iz]=value;
}

void C3DArray::IncrementElement(int isx,int ix,int isy,int iy,int isz,int iz,double value){
	F[isx][ix][isy][iy][isz][iz]+=value;
}

double C3DArray::GetElement_NoInterpolation(double x,double y,double z){
	int ix,iy,iz,isx,isy,isz;
	double answer;
	answer=0.0;
	ix=int(floor(fabs(x)/DELX));
	if(ix<NXMAX){
		iy=int(floor(fabs(y)/DELY));
		if(iy<NYMAX){
			iz=int(floor(fabs(z)/DELY));
			if(iz<NZMAX){
				isx=isy=isz=0;
				if(!XSYM && x<-1.0E-10) isx=1;
				if(!YSYM && y<-1.0E-10) isy=1;
				if(!ZSYM && z<-1.0E-10) isz=1;
				answer=F[isx][ix][isy][iy][isz][iz];
			}
		}
	}
	return answer;
}

double C3DArray::GetElement(double x,double y,double z){
	int ixa,ixb,iya,iyb,iza,izb,isxa,isxb,isya,isyb,isza,iszb;
	double fx,fy,fz,interpolate;
	ixb=int(floor(0.5+fabs(x)/DELX));
	iyb=int(floor(0.5+fabs(y)/DELY));
	izb=int(floor(0.5+fabs(z)/DELZ));
	ixa=ixb-1;
	iya=iyb-1;
	iza=izb-1;
	isxb=isyb=iszb=0;
	if(x<1.0E-10 && !XSYM) isxb=1;
	if(y<1.0E-10 && !YSYM) isyb=1;
	if(z<1.0E-10 && !ZSYM) iszb=1;
	isxa=isxb;
	isya=isyb;
	isza=iszb;

	if(ixa<0){
		ixa=0;
		if(!XSYM){
			if(isxb==1) isxa=0;
			if(isxb==0) isxa=1;
		}
	}
	if(iya<0){
		iya=0;
		if(!YSYM){
			if(isyb==1) isya=0;
			if(isyb==0) isya=1;
		}
	}
	if(iza<0){
		iza=0;
		if(!ZSYM){
			if(iszb==1) isza=0;
			if(iszb==0) isza=1;
		}
	}

	fx=((0.5+ixb)*DELX-fabs(x))/DELX;
	fy=((0.5+iyb)*DELY-fabs(y))/DELY;
	fz=((0.5+izb)*DELZ-fabs(z))/DELZ;

	if(ixb==NXMAX) ixb=ixa;
	if(iyb==NYMAX) iyb=iya;
	if(izb==NZMAX) izb=iza;

	interpolate=(1.0-fx)*(1.0-fy)*(1.0-fz)
		*GetElement(isxb,ixb,isyb,iyb,iszb,izb);
	interpolate+=fx*(1.0-fy)*(1.0-fz)
		*GetElement(isxa,ixa,isyb,iyb,iszb,izb);
	interpolate+=(1.0-fx)*fy*(1.0-fz)
		*GetElement(isxb,ixb,isya,iya,iszb,izb);
	interpolate+=(1.0-fx)*(1.0-fy)*fz
		*GetElement(isxb,ixb,isyb,iyb,isza,iza);
	interpolate+=fx*fy*(1.0-fz)
		*GetElement(isxa,ixa,isya,iya,iszb,izb);
	interpolate+=fx*(1.0-fy)*fz
		*GetElement(isxa,ixa,isyb,iyb,isza,iza);
	interpolate+=(1.0-fx)*fy*fz
		*GetElement(isxb,ixb,isya,iya,isza,iza);
	interpolate+=fx*fy*fz
		*GetElement(isxa,ixa,isya,iya,isza,iza);
	return interpolate;
}

void C3DArray::SetElement(double x,double y,double z,double value){
	int ix,iy,iz,isx,isy,isz;
	ix=int(floor(fabs(x)/DELX));
	if(ix<NXMAX){
		iy=int(floor(fabs(y)/DELY));
		if(iy<NYMAX){
			iz=int(floor(fabs(z)/DELY));
			if(iz<NZMAX){
				isx=isy=isz=0;
				if(!XSYM && x<-1.0E-15) isx=1;
				if(!YSYM && y<-1.0E-15) isy=1;
				if(!ZSYM && z<-1.0E-15) isz=1;
				F[isx][ix][isy][iy][isz][iz]=value;
			}
		}
	}
}

void C3DArray::IncrementElement(double x,double y,double z,double value){
	int ix,iy,iz,isx,isy,isz;
	ix=int(floor(fabs(x)/DELX));
	if(ix<NXMAX){
		iy=int(floor(fabs(y)/DELY));
		if(iy<NYMAX){
			iz=int(floor(fabs(z)/DELZ));
			if(iz<NZMAX){
				isx=isy=isz=0;
				if(!XSYM && x<-1.0E-10)
					isx=1;
				if(!YSYM && y<-1.0E-10)
					isy=1;
				if(!ZSYM && z<-1.0E-10)
					isz=1;
				F[isx][ix][isy][iy][isz][iz]+=value;
			}
		}
	}
}

int C3DArray::GetNXMAX(){
	return NXMAX;
}

int C3DArray::GetNYMAX(){
	return NYMAX;
}

int C3DArray::GetNZMAX(){
	return NZMAX;
}

double C3DArray::GetDELX(){
	return DELX;
}

double C3DArray::GetDELY(){
	return DELY;
}

double C3DArray::GetDELZ(){
	return DELZ;
}

void C3DArray::PrintPars(){
	printf("_______ 3DArray Pars ________\n");
	printf("XYZSYM=(%d,%d,%d)\n",XSYM,YSYM,ZSYM);
	printf("NXYZMAX=(%d,%d,%d)\n",NXMAX,NYMAX,NZMAX);
	printf("DELXYZ=(%g,%g,%g)\n",DELX,DELY,DELZ);
	printf("_____________________________\n");
}

void C3DArray::WriteArray(string dirname){
	int isx,isy,isz,nsx,nsy,nsz,ix,iy,iz;
	FILE *fptr;
	char shellcommand[300];

	snprintf(shellcommand,300,"mkdir -p %s",dirname.c_str());
	system(shellcommand);
	snprintf(shellcommand,300,"rm -f %s/*.tmp",dirname.c_str());
	system(shellcommand);

	char filename[160];
	snprintf(filename,160,"%s/3Darraypars.dat",dirname.c_str());
	fptr=fopen(filename,"w");
	fprintf(fptr,"bool XSYM %d\n",XSYM);
	fprintf(fptr,"bool YSYM %d\n",YSYM);
	fprintf(fptr,"bool ZSYM %d\n",ZSYM);
	fprintf(fptr,"int NXMAX %d\n",NXMAX);
	fprintf(fptr,"int NYMAX %d\n",NYMAX);
	fprintf(fptr,"int NZMAX %d\n",NZMAX);
	fprintf(fptr,"double DELX %g\n",DELX);
	fprintf(fptr,"double DELY %g\n",DELY);
	fprintf(fptr,"double DELZ %g\n",DELZ);
	fclose(fptr);
	nsx=nsy=nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	for(ix=0;ix<NXMAX;ix++){
		for(iy=0;iy<NYMAX;iy++){
			snprintf(filename,160,"%s/ix%d_iy%d.tmp",
				dirname.c_str(),ix,iy);
			fptr=fopen(filename,"w");
			for(iz=0;iz<NZMAX;iz++){
				for(isz=0;isz<nsz;isz++){
					for(isy=0;isy<nsy;isy++){
						for(isx=0;isx<nsx;isx++){
							fprintf(fptr,"%17.10e ",F[isx][ix][isy][iy][isz][iz]);
						}
					}
				}
				fprintf(fptr,"\n");
			}
			fclose(fptr);

		}
	}
}

void C3DArray::WriteArray_PHENIX(string filename){
	//This writes entire array to one file, and also writes error array as all 1.0s
	int ix0,iy0,iz0,ix,iy,iz,ixx,iyy,izz,isx,isy,isz;
	double x,y,z;
	FILE *fptr;
	fptr=fopen(filename.c_str(),"w");
	ix0=-NXMAX; iy0=-NYMAX; iz0=-NZMAX;
	if(XSYM) ix0=0;
	if(YSYM) iy0=0;
	if(ZSYM) iz0=0;
	for(ixx=ix0;ixx<NXMAX;ixx++){
		x=(ixx+0.5)*DELX;
		isx=0;
		ix=ixx;
		if(ixx<0){
			isx=1;
			ix=-ixx-1;
		}
		for(iyy=iy0;iyy<NYMAX;iyy++){
			isy=0;
			iy=iyy;
			if(iyy<0){
				isy=1;
				iy=-iyy-1;
			}
			y=(iyy+0.5)*DELY;
			for(izz=iz0;izz<NZMAX;izz++){
				isz=0;
				iz=izz;
				if(izz<0){
					isz=1;
					iz=-izz-1;
				}
				z=(izz+0.5)*DELZ;
				fprintf(fptr,"%g %g %g %13.7e 1.0\n",x,y,z,F[isx][ix][isy][iy][isz][iz]);
			}
		}
	}
	fclose(fptr);
}


void C3DArray::ReadArray(string dirname){
	int isx,isy,isz,nsx,nsy,nsz,ix,iy,iz;
	bool XSYMb,YSYMb,ZSYMb;
	int NXMAXb,NYMAXb,NZMAXb;
	double DELXb,DELYb,DELZb;
	CparameterMap apars;
	NXMAXb=NXMAX; NYMAXb=NYMAX; NZMAXb=NZMAX;
	DELXb=DELX; DELYb=DELY; DELZb=DELZ;
	XSYMb=XSYM; YSYMb=YSYM; ZSYMb=ZSYM;
	FILE *fptr;
	char filename[160];
	snprintf(filename,160,"%s/3Darraypars.dat",dirname.c_str());
	ReadPars(filename);
	if(NXMAXb!=NXMAX || NYMAXb!=NYMAX|| NZMAXb!=NZMAX
		|| fabs(DELXb-DELX)>1.0E-10 || fabs(DELYb-DELY)>1.0E-10
		|| fabs(DELZb-DELZ)>1.0E-10 || XSYMb!=XSYM
	|| YSYMb!=YSYM || ZSYMb!=ZSYM){
		printf("Warning: In ReadArray parameters have changed, new array being created\n");
		DeleteArray();
		CreateArray();
	}
	nsx=nsy=nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	for(ix=0;ix<NXMAX;ix++){
		for(iy=0;iy<NYMAX;iy++){
			snprintf(filename,160,"%s/ix%d_iy%d.tmp",dirname.c_str(),ix,iy);
			fptr=fopen(filename,"r");
			for(iz=0;iz<NZMAX;iz++){
				for(isz=0;isz<nsz;isz++){
					for(isy=0;isy<nsy;isy++){
						for(isx=0;isx<nsx;isx++){
							fscanf(fptr,"%lf ",&F[isx][ix][isy][iy][isz][iz]);
						}
					}
				}
			}
			fclose(fptr);

		}
	}
}

double C3DArray::GetBiggest(){
	int isx,isy,isz,ix,iy,iz,nsx,nsy,nsz;
	double biggest=0.0;
	nsx=nsy=nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<NXMAX;ix++){
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<NYMAX;iy++){
					for(iz=0;iz<NZMAX;iz++){
						for(isz=0;isz<nsz;isz++){
							if(fabs(F[isx][ix][isy][iy][isz][iz])>fabs(biggest))
								biggest=F[isx][ix][isy][iy][isz][iz];
						}
					}
				}
			}
		}
	}
	return biggest;
}

void C3DArray::Randomize(double c){
	int isx,isy,isz,ix,iy,iz,nsx,nsy,nsz;
	if(randy==NULL) randy=new Crandy(-1234);
	nsx=nsy=nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<NXMAX;ix++){
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<NYMAX;iy++){
					for(iz=0;iz<NZMAX;iz++){
						for(isz=0;isz<nsz;isz++){
							F[isx][ix][isy][iy][isz][iz]=c*(1.0-2.0*randy->ran());
						}
					}
				}
			}
		}
	}

}

void C3DArray::RandomizeGaussian(double c){
	int isx,isy,isz,ix,iy,iz,nsx,nsy,nsz;
	if(randy==NULL) randy=new Crandy(-1234);
	nsx=nsy=nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<NXMAX;ix++){
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<NYMAX;iy++){
					for(iz=0;iz<NZMAX;iz++){
						for(isz=0;isz<nsz;isz++){
							F[isx][ix][isy][iy][isz][iz]=c*randy->ran_gauss();
						}
					}
				}
			}
		}
	}

}

void C3DArray::MakeConstant(double c){
	int isx,isy,isz,ix,iy,iz,nsx,nsy,nsz;
	nsx=nsy=nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<NXMAX;ix++){
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<NYMAX;iy++){
					for(iz=0;iz<NZMAX;iz++){
						for(isz=0;isz<nsz;isz++){
							F[isx][ix][isy][iy][isz][iz]=c;
						}
					}
				}
			}
		}
	}
}

void C3DArray::PrintProjections(){
	int isx,isy,isz,nsx,nsy,nsz,NMAX,i,nav;
	double fbar;
	nsx=nsy=nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	NMAX=NXMAX;
	if(NYMAX>NMAX) NMAX=NXMAX;
	if(NZMAX>NMAX) NMAX=NZMAX;
	printf("bool XSYM %d ",XSYM);
	printf("bool YSYM %d ",YSYM);
	printf("bool ZSYM %d\n",ZSYM);
	printf("int NXMAX %d ",NXMAX);
	printf("int NYMAX %d ",NYMAX);
	printf("int NZMAX %d\n",NZMAX);
	printf("double DELX %g ",DELX);
	printf("double DELY %g ",DELY);
	printf("double DELZ %g\n",DELZ);
	for(i=0;i<NMAX;i++){
		printf("%3d | ",i);
		for(isx=0;isx<nsx;isx++){
			fbar=0.0;
			if(i<NXMAX){
				nav=0;
				for(isy=0;isy<nsy;isy++){
					for(isz=0;isz<nsz;isz++){
						nav+=1;
						fbar+=F[isx][i][isy][0][isz][0];
					}
				}
				fbar=fbar/double(nav);
			}
			printf("%9.3e ",fbar);
		}
		printf("| ");
		for(isy=0;isy<nsy;isy++){
			fbar=0.0;
			if(i<NYMAX){
				nav=0;
				for(isx=0;isx<nsx;isx++){
					for(isz=0;isz<nsz;isz++){
						nav+=1;
						fbar+=F[isx][0][isy][i][isz][0];
					}
				}
				fbar=fbar/double(nav);
			}
			printf("%9.3e ",fbar);
		}
		printf("| ");
		for(isz=0;isz<nsz;isz++){
			fbar=0.0;
			if(i<NZMAX){
				nav=0;
				for(isy=0;isy<nsy;isy++){
					for(isx=0;isx<nsx;isx++){
						nav+=1;
						fbar+=F[isx][0][isy][0][isz][i];
					}
				}
				fbar=fbar/double(nav);
			}
			printf("%9.3e ",fbar);
		}
		printf("\n");
	}
}

void C3DArray::WriteProjections(string filename){
	int isx,isy,isz,nsx,nsy,nsz,NMAX,i,nav;
	double fbar;
	FILE *fptr=fopen(filename.c_str(),"w");
	nsx=nsy=nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	NMAX=NXMAX;
	if(NYMAX>NMAX) NMAX=NXMAX;
	if(NZMAX>NMAX) NMAX=NZMAX;
	fprintf(fptr,"bool XSYM %d ",XSYM);
	fprintf(fptr,"bool YSYM %d ",YSYM);
	fprintf(fptr,"bool ZSYM %d\n",ZSYM);
	fprintf(fptr,"int NXMAX %d ",NXMAX);
	fprintf(fptr,"int NYMAX %d ",NYMAX);
	fprintf(fptr,"int NZMAX %d\n",NZMAX);
	fprintf(fptr,"double DELX %g ",DELX);
	fprintf(fptr,"double DELY %g ",DELY);
	fprintf(fptr,"double DELZ %g\n",DELZ);
	for(i=0;i<NMAX;i++){
		printf("%3d | ",i);
		for(isx=0;isx<nsx;isx++){
			fbar=0.0;
			if(i<NXMAX){
				nav=0;
				for(isy=0;isy<nsy;isy++){
					for(isz=0;isz<nsz;isz++){
						nav+=1;
						fbar+=F[isx][i][isy][0][isz][0];
					}
				}
				fbar=fbar/double(nav);
			}
			printf("%9.3e ",fbar);
		}
		printf("| ");
		for(isy=0;isy<nsy;isy++){
			fbar=0.0;
			if(i<NYMAX){
				nav=0;
				for(isx=0;isx<nsx;isx++){
					for(isz=0;isz<nsz;isz++){
						nav+=1;
						fbar+=F[isx][0][isy][i][isz][0];
					}
				}
				fbar=fbar/double(nav);
			}
			printf("%9.3e ",fbar);
		}
		printf("| ");
		for(isz=0;isz<nsz;isz++){
			fbar=0.0;
			if(i<NZMAX){
				nav=0;
				for(isy=0;isy<nsy;isy++){
					for(isx=0;isx<nsx;isx++){
						nav+=1;
						fbar+=F[isx][0][isy][0][isz][i];
					}
				}
				fbar=fbar/double(nav);
			}
			printf("%9.3e ",fbar);
		}
		printf("\n");
	}
	fclose(fptr);
}

void C3DArray::CalcMoments(double roff[3],double r2[3][3]){
	int ix,iy,iz,isx,isy,isz,nsx,nsy,nsz;
	double norm;
	int i,j;
	double r[3];
	nsx=nsy=nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	for(i=0;i<3;i++){
		roff[i]=0.0;
		for(j=0;j<3;j++) r2[i][j]=0.0;
	}
	norm=0.0;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<NXMAX;ix++){
			r[0]=(ix+0.5)*DELX;
			if(isx==1) r[0]=-r[0];
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<NYMAX;iy++){
					r[1]=(iy+0.5)*DELY;
					if(isy==1) r[1]=-r[1];
					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<NZMAX;iz++){
							r[2]=(iz+0.5)*DELZ;
							if(isz==0) r[2]=-r[2];
							norm+=F[isx][ix][isy][iy][isz][iz];
							for(i=0;i<3;i++){
								roff[i]+=r[i]*F[isx][ix][isy][iy][isz][iz];
								for(j=i;j<3;j++) r2[i][j]+=r[i]*r[j]*F[isx][ix][isy][iy][isz][iz];
							}
						}
					}
				}
			}
		}
	}
	if(XSYM) roff[0]=r2[0][1]=r2[0][2]=0.0;
	if(YSYM) roff[1]=r2[0][1]=r2[1][2]=0.0;
	if(ZSYM) roff[2]=r2[0][2]=r2[1][2]=0.0;
	//printf("Offsets: ");
	for(i=0;i<3;i++){
		roff[i]=roff[i]/norm;
		//printf("%11.3e ",roff[i]);
		for(j=0;j<i;j++) r2[i][j]=r2[j][i];
	}
	printf("\n<r_ir_j>:\n");
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			r2[i][j]=r2[i][j]/norm;
			r2[i][j]=r2[i][j]-roff[i]*roff[j];
			//printf("%11.3e ",r2[i][j]);
		}
		//printf("\n");
	}
	//printf("xyz_rms = (%g,%g,%g)\n",sqrt(r2[0][0]),sqrt(r2[1][1]),sqrt(r2[2][2]));
}

void C3DArray::PrintMoments(){
	double roff[3],r2[3][3];
	int ix,iy,iz,isx,isy,isz,nsx,nsy,nsz;
	double norm;
	int i,j;
	double r[3];
	nsx=nsy=nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	for(i=0;i<3;i++){
		roff[i]=0.0;
		for(j=0;j<3;j++) r2[i][j]=0.0;
	}
	norm=0.0;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<NXMAX;ix++){
			r[0]=(ix+0.5)*DELX;
			if(isx==1) r[0]=-r[0];
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<NYMAX;iy++){
					r[1]=(iy+0.5)*DELY;
					if(isy==1) r[1]=-r[1];
					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<NZMAX;iz++){
							r[2]=(iz+0.5)*DELZ;
							if(isz==0) r[2]=-r[2];
							norm+=F[isx][ix][isy][iy][isz][iz];
							for(i=0;i<3;i++){
								roff[i]+=r[i]*F[isx][ix][isy][iy][isz][iz];
								for(j=i;j<3;j++) r2[i][j]+=r[i]*r[j]*F[isx][ix][isy][iy][isz][iz];
							}
						}
					}
				}
			}
		}
	}
	if(XSYM) roff[0]=r2[0][1]=r2[0][2]=0.0;
	if(YSYM) roff[1]=r2[0][1]=r2[1][2]=0.0;
	if(ZSYM) roff[2]=r2[0][2]=r2[1][2]=0.0;
	printf("Offsets: ");
	for(i=0;i<3;i++){
		roff[i]=roff[i]/norm;
		printf("%11.3e ",roff[i]);
		for(j=0;j<i;j++) r2[i][j]=r2[j][i];
	}
	printf("\n<r_ir_j>:\n");
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			r2[i][j]=r2[i][j]/norm;
			r2[i][j]=r2[i][j]-roff[i]*roff[j];
			printf("%11.3e ",r2[i][j]);
		}
		printf("\n");
	}
	printf("xyz_rms = (%g,%g,%g)\n",sqrt(r2[0][0]),sqrt(r2[1][1]),sqrt(r2[2][2]));
}

void C3DArray::DivideByArray(C3DArray *threed_denom){
	int isx,ix,isy,iy,isz,iz;
	int nsx=2,nsy=2,nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	double r[4];
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<NXMAX;ix++){
			r[0]=(ix+0.5)*DELX;
			if(isx==1)
				r[0]=-r[0];
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<NYMAX;iy++){
					r[1]=(iy+0.5)*DELY;
					if(isy==1)
						r[1]=-r[1];
					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<NZMAX;iz++){
							r[2]=(iz+0.5)*DELZ;
							if(isz==0)
								r[2]=-r[2];
							if(fabs(threed_denom->F[isx][ix][isy][iy][isz][iz])>1.0E-8)
								F[isx][ix][isy][iy][isz][iz]=F[isx][ix][isy][iy][isz][iz]/threed_denom->F[isx][ix][isy][iy][isz][iz];
							else{
								//CLog::Info("Denom zero for C3DArray::DivideByArray\n");
								F[isx][ix][isy][iy][isz][iz]=1.0;
							}
						}
					}
				}
			}
		}
	}
}
