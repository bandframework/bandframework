#include "msu_commonutils/arrays.h"
#include "msu_commonutils/sf.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/randy.h"

using namespace std;

Crandy *CCHArray::randy=NULL;
CCHCalc *CCHArray::chcalc=NULL;

CCHArray::CCHArray(int LMAXset,int NRADIALset,double RADSTEPset){
  if(chcalc==NULL) chcalc=new CCHCalc();
  NRADIAL=NRADIALset;
  LMAX=LMAXset;
  RADSTEP=RADSTEPset;
  XSYM=YSYM=ZSYM=false;
  dlx=dly=dlz=1;
  CreateArray();
}

CCHArray::CCHArray(int LMAXset,int NRADIALset,double RADSTEPset,bool XSYMset,bool YSYMset,bool ZSYMset){
	if(chcalc==NULL) chcalc=new CCHCalc();
	NRADIAL=NRADIALset;
	RADSTEP=RADSTEPset;
	LMAX=LMAXset;
	XSYM=XSYMset;
	YSYM=YSYMset;
	ZSYM=ZSYMset;
	dlx=dly=dlz=1;
	if(XSYM) dlx=2;
	if(YSYM) dly=2;
	if(ZSYM) dlz=2;
	CreateArray();
}

CCHArray::CCHArray(string arrayparsfilename){
	dlx=dly=dlz=1;
	if(chcalc==NULL) chcalc=new CCHCalc();
	CparameterMap apars;
	
	apars.ReadParsFromFile(arrayparsfilename);
	if(apars.getB("IDENTICAL",false)){
		apars.set("XSYM",true);
		apars.set("YSYM",true);
		apars.set("ZSYM",true);
	}
	
	NRADIAL=apars.getI("NRADIAL",-999);
	LMAX=apars.getI("LMAX",-999);
	RADSTEP=apars.getD("RADSTEP",-999);
	XSYM=apars.getB("XSYM",false);
	YSYM=apars.getB("YSYM",false);
	ZSYM=apars.getB("ZSYM",false);
	if(XSYM) dlx=2;
	if(YSYM) dly=2;
	if(ZSYM) dlz=2;
	PrintPars();
	
	CreateArray();
}

void CCHArray::CreateArray(){
	int lx,ly,lz,ir;
	//PrintPars();
	A=new double ***[LMAX+1];
	for(lx=0;lx<=LMAX;lx+=dlx){
		A[lx]=new double **[LMAX+1-lx];
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			A[lx][ly]=new double *[LMAX+1-lx-ly];
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				A[lx][ly][lz]=new double[NRADIAL];
				for(ir=0;ir<NRADIAL;ir++) A[lx][ly][lz][ir]=0.0;
			}
		}
	}
	
}

CCHArray::~CCHArray(){
	int lx,ly,lz;
	for(lx=0;lx<=LMAX;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz) delete [] A[lx][ly][lz];
			delete [] A[lx][ly];
		}
		delete [] A[lx];
	}
	delete A;
}

int CCHArray::GetLMAX(){
	return LMAX;
}

void CCHArray::SetLMAX(int LMAXset){
	if(LMAXset<LMAX) LMAX=LMAXset;
}

int CCHArray::GetNRADIAL(){
	return NRADIAL;
}

double CCHArray::GetRADSTEP(){
	return RADSTEP;
}

void CCHArray::SetRADSTEP(double RADSTEPset){
	RADSTEP=RADSTEPset;
}

double CCHArray::GetElement(int lx,int ly,int lz,int ir){
	int L=lx+ly+lz;
	if(L<=LMAX && lx%dlx==0 && ly%dly==0 && lz%dlz==0 && ir<NRADIAL) 
		return A[lx][ly][lz][ir];
		else{
			printf("WARNING: Tried to get non-existent element in CCHArray\n");
			printf("LMAX=%d, XSYM=%d, YSYM=%d, ZSYM=%d, ir=%d\n",
				LMAX,int(XSYM),int(YSYM),int(ZSYM),ir);
			printf("l=(%d,%d,%d), ir=%d\n",lx,ly,lz,ir);
			return 0.0;
		}
}

double CCHArray::GetElement(int lx,int ly,int lz,double r){
	int ir=int(floor(r/RADSTEP));
	return GetElement(lx,ly,lz,ir);
}

void CCHArray::SetElement(int lx,int ly,int lz,int ir,double Element){
	int L=lx+ly+lz;
	if(L<=LMAX && L>=0 && lx%dlx==0 && ly%dly==0 && lz%dlz==0 && ir<NRADIAL) 
		A[lx][ly][lz][ir]=Element;
		else{
			printf("Tried to set element out of bounds in CCHArray\n");
			printf("LMAX=%d, XSYM=%d, YSYM=%d, ZSYM=%d, ir=%d\n",
				LMAX,int(XSYM),int(YSYM),int(ZSYM),ir);
			printf("l=(%d,%d,%d), ir=%d\n",lx,ly,lz,ir);
			exit(1);
			
		}
}

void CCHArray::SetElement(int lx,int ly,int lz,double r,double Element){
	int ir=int(floor(r/RADSTEP));
	SetElement(lx,ly,lz,ir,Element);
}

void CCHArray::IncrementElement(int lx,int ly,int lz,int ir,double increment){
	int L=lx+ly+lz;
	if(L<=LMAX && L>=0 && lx%dlx==0 && ly%dly==0 && lz%dlz==0 && ir<NRADIAL) {
		A[lx][ly][lz][ir]+=increment;
	}
	else{
		printf("OUT OF BOUNDS IN INCREMENTELEMENT\n");
		printf("ell=(%d,%d,%d), ir=%d\n",lx,ly,lz,ir);
		exit(1);
	}
	
}

void CCHArray::IncrementElement(int lx,int ly,int lz,double r,double increment){
	int ir=int(floor(r/RADSTEP));
	IncrementElement(lx,ly,lz,ir,increment);
}

void CCHArray::ScaleArray(double scalefactor){
	int lx,ly,lz,ir;
	for(lx=0;lx<=LMAX;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				for(ir=0;ir<NRADIAL;ir++) A[lx][ly][lz][ir]*=scalefactor;
			}
		}
	}
}

void CCHArray::ScaleArray(double scalefactor, int ir){
	int lx,ly,lz;
	for(lx=0;lx<=LMAX;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				A[lx][ly][lz][ir]*=scalefactor;
			}
		}
	}
}

void CCHArray::ZeroArray_Partial(int LMAX_partial){
	int lx,ly,lz,ir;
	for(lx=0;lx<=LMAX_partial;lx+=dlx){
		for(ly=0;ly<=LMAX_partial-lx;ly+=dly){
			for(lz=0;lz<=LMAX_partial-lx-ly;lz+=dlz){
				for(ir=0;ir<NRADIAL;ir++) A[lx][ly][lz][ir]=0.0;
			}
		}
	}
}

void CCHArray::ZeroArray_Partial(int LMAX_partial,int ir){
	int lx,ly,lz;
	for(lx=0;lx<=LMAX_partial;lx+=dlx){
		for(ly=0;ly<=LMAX_partial-lx;ly+=dly){
			for(lz=0;lz<=LMAX_partial-lx-ly;lz+=dlz){
				A[lx][ly][lz][ir]=0.0;
			}
		}
	}
}

void CCHArray::ZeroArray(){
	int lx,ly,lz,ir;
	for(lx=0;lx<=LMAX;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				for(ir=0;ir<NRADIAL;ir++) A[lx][ly][lz][ir]=0.0;
			}
		}
	}
}

void CCHArray::ZeroArray(int lx,int ly,int lz){
	int ir;
	for(ir=0;ir<NRADIAL;ir++) A[lx][ly][lz][ir]=0.0;
}

void CCHArray::ZeroArray(int ir){
	int lx,ly,lz;
	for(lx=0;lx<=LMAX;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz)
				A[lx][ly][lz][ir]=0.0;
		}
	}
}

void CCHArray::Print(int lx,int ly,int lz){
	for(int ir=0;ir<NRADIAL;ir++){
		printf("%6.3f: %10.3e\n",(0.5+ir)*RADSTEP,A[lx][ly][lz][ir]);
	}
}

void CCHArray::PrintProjections(){
	int lx,ly,lz;
	double a0,ax,ay,az;
	printf("_____   lx=ly=lz=0     x-proj      y-proj      z-proj  _____\n");
	for(int ir=0;ir<NRADIAL;ir++){
		a0=A[0][0][0][ir];
		ax=ay=az=0.0;
		ly=lz=0;
		for(lx=0;lx<=LMAX;lx+=dlx) ax+=A[lx][ly][lz][ir];
		lx=lz=0;
		for(ly=0;ly<=LMAX;ly+=dly) ay+=A[lx][ly][lz][ir];
		ly=lx=0;
		for(lz=0;lz<=LMAX;lz+=dlz) az+=A[lx][ly][lz][ir];
		printf("%7.3f: %10.3e  %10.3e  %10.3e  %10.3e\n",(ir+0.5)*RADSTEP,a0,ax,ay,az);
	}
}

void CCHArray::WriteProjections(string filename){
	int lx,ly,lz;
	double a0,ax,ay,az;
	FILE *fptr=fopen(filename.c_str(),"w");
	printf("_____   lx=ly=lz=0     x-proj      y-proj      z-proj  _____\n");
	for(int ir=0;ir<NRADIAL;ir++){
		a0=A[0][0][0][ir];
		ax=ay=az=0.0;
		ly=lz=0;
		for(lx=0;lx<=LMAX;lx+=dlx) ax+=A[lx][ly][lz][ir];
		lx=lz=0;
		for(ly=0;ly<=LMAX;ly+=dly) ay+=A[lx][ly][lz][ir];
		ly=lx=0;
		for(lz=0;lz<=LMAX;lz+=dlz) az+=A[lx][ly][lz][ir];
		fprintf(fptr,"%7.3f: %10.3e  %10.3e  %10.3e  %10.3e\n",(ir+0.5)*RADSTEP,a0,ax,ay,az);
	}
	fclose(fptr);
}

void CCHArray::GetProjections(double **B){
	int lx,ly,lz;
	double a0,ax,ay,az;
	printf("_____   lx=ly=lz=0     x-proj      y-proj      z-proj  _____\n");
	for(int ir=0;ir<NRADIAL;ir++){
		ax=ay=az=0.0;
		ly=lz=0;
		a0=A[0][0][0][ir];
		for(lx=0;lx<=LMAX;lx+=dlx) ax+=A[lx][ly][lz][ir];
		lx=lz=0;
		for(ly=0;ly<=LMAX;ly+=dly) ay+=A[lx][ly][lz][ir];
		ly=lz=0;
		for(lz=0;lz<=LMAX;lz+=dlz) az+=A[lx][ly][lz][ir];
		B[0][ir]=a0; B[1][ir]=ax; B[2][ir]=ay; B[3][ir]=az;
		//printf("%6.3f: %10.3e  %10.3e  %10.3e  %10.3e\n",(ir+0.5)*RADSTEP,a0,ax,ay,az);
	}
}

void CCHArray::PrintArrayFixedIR(int ir){
	PrintArrayFixedIR(LMAX,ir);
}

void CCHArray::PrintArrayFixedIR(int LMAXPrint,int ir){
	int L,lx,ly,lz,dL=1;
	if(XSYM && YSYM && ZSYM) dL=2;
	if(LMAXPrint>LMAX) LMAXPrint=LMAX;
	
	printf("\n______________________________________________________\n");
	for(L=0;L<=LMAXPrint;L+=dL){
		printf("     L=%d\n",L);
		printf(" lx\\ly:");
		for(ly=0;ly<=L;ly+=dly) printf(" %4d       ",ly);
		printf("\n");
		
		for(lx=0;lx<=L;lx+=dlx){
			printf(" %3d ",lx);
			for(ly=0;ly<=L-lx;ly+=dly){
				lz=L-lx-ly;
				if(ZSYM==0 || lz%2==0)
					printf(" %10.3e ",A[lx][ly][lz][ir]);
				else printf("%10.3e ",0.0);
			}
			printf("\n");
		}
		printf("_________________________________________\n");
	}
}

double CCHArray::GetBiggest(int ir){
	int lx,ly,lz;
	double biggy=0.0;
	for(lx=0;lx<=LMAX;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				if(fabs(A[lx][ly][lz][ir])>fabs(biggy))
					biggy=fabs(A[lx][ly][lz][ir]);
			}
		} 
	}
	return biggy;
}

bool CCHArray::GetXSYM(){
	if(XSYM) return true;
	else return false;
}

bool CCHArray::GetYSYM(){
	if(YSYM) return true;
	else return false;
}

bool CCHArray::GetZSYM(){
	if(ZSYM) return true;
	else return false;
}

void CCHArray::WriteAX(string dirname){
	char filename[160],shellcommand[200];
	int ir,lx,ly,lz;
	FILE *fptr;
	
	snprintf(shellcommand,200,"mkdir -p %s",dirname.c_str());
	system(shellcommand);
	
	for(lx=0;lx<=1;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				snprintf(filename,160,"%s/lx%d_ly%d_lz%d.tmp",dirname.c_str(),lx,ly,lz);
				fptr=fopen(filename,"w");
				fprintf(fptr,"%d %g\n",NRADIAL,RADSTEP);
				for(ir=0;ir<NRADIAL;ir++){
					fprintf(fptr,"%g\n",A[lx][ly][lz][ir]);
				}
				fclose(fptr);
			}
		}
	}
}

void CCHArray::ReadAX(string dirname){
	char filename[160],shellcommand[320];
	int ir,lx,ly,lz,NRADIALread;
	double aa,RADSTEPread;
	FILE *fptr;
	for(lx=0;lx<=1;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				snprintf(filename,160,"%s/lx%d_ly%d_lz%d.tmp",dirname.c_str(),lx,ly,lz);
				//printf("READING: L=(%d,%d,%d), filename=%s\n",lx,ly,lz,filename);
				snprintf(shellcommand,320,
					"if [ ! -e %s ]; then echo Reading Error: %s does not exist; fi",filename,filename);
				system(shellcommand);
				fptr=fopen(filename,"r");
				fscanf(fptr,"%d %lf",&NRADIALread,&RADSTEPread);
				if(fabs(RADSTEPread-RADSTEP)>1.0E-6 || NRADIALread<NRADIAL){
					printf("Mesh in %s out of whack: RADSTEP=%g =? %g, NRADIAL=%d=?%d\n",
						filename,RADSTEP,RADSTEPread,NRADIAL,NRADIALread);
					exit(1);
				}
				for(ir=0;ir<NRADIAL;ir++){
					fscanf(fptr,"%lf",&aa);
					A[lx][ly][lz][ir]=aa;
				}
				fclose(fptr);
			}
		}
	}
	FillRemainderX();
	//printf("Reading finished\n");
}

void CCHArray::ReadAllA(string dirname){
	char filename[160],shellcommand[320];
	int ir,lx,ly,lz,NRADIALread;
	double aa,RADSTEPread;
	FILE *fptr;
	
	for(lx=0;lx<=LMAX;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				snprintf(filename,160,"%s/lx%d_ly%d_lz%d.tmp",dirname.c_str(),lx,ly,lz);
				printf("READING: L=(%d,%d,%d), filename=%s\n",lx,ly,lz,filename);
				snprintf(shellcommand,320,
					"if [ ! -e %s ]; then echo Reading Error: %s, does not exist; fi",filename,filename);
				system(shellcommand);
				fptr=fopen(filename,"r");
				fscanf(fptr,"%d %lf",&NRADIALread,&RADSTEPread);
				if(fabs(RADSTEPread-RADSTEP)>1.0E-6 || NRADIALread<NRADIAL){
					printf("Mesh in %s out of whack: RADSTEP=%g =? %g, NRADIAL=%d=?%d\n",
						filename,RADSTEP,RADSTEPread,NRADIAL,NRADIALread);
					exit(1);
				}
				for(ir=0;ir<NRADIAL;ir++){
					fscanf(fptr,"%lf\n",&aa);
					A[lx][ly][lz][ir]=aa;
				}
				fclose(fptr);
			}
		}
	}
}

void CCHArray::PrintPars(){
	printf("XSYM=%d, YSYM=%d, ZSYM=%d\n",int(XSYM),int(YSYM),int(ZSYM));
	printf("LMAX=%d, NRADIAL=%d, RADSTEP=%g\n",LMAX,NRADIAL,RADSTEP);
}

void CCHArray::IncrementAExpArray(double x,double y,double z,double weight){
	double ex,ey,ez;
	double r=sqrt(x*x+y*y+z*z);
	int ir=int(floor(r/RADSTEP));
	if(ir<NRADIAL){
		ex=x/r; ey=y/r; ez=z/r;
		IncrementAExpArrayFromE(ex,ey,ez,weight,ir);
	}
}

void CCHArray::IncrementAExpArrayFromE(double ex,double ey,double ez,double weight,int ir){
	int L,lx,ly,lz;
	double lfact;
	for(lx=0;lx<=1;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				L=lx+ly+lz;
				lfact=chcalc->DoubleFactorial(2*L+1)/chcalc->Factorial(L);
				IncrementElement(lx,ly,lz,ir,
					weight*lfact*chcalc->GetAFromE(lx,ly,lz,ex,ey,ez));
			}
		}
	}
	FillRemainderX(ir);
}

void CCHArray::AltIncrementAExpArrayFromE(double ex,double ey,double ez,double weight,int ir){
	int L,lx,ly,lz;
	double cL,cLL,delB,lfact;
	double ***B;
	B=new double**[2];
	for(lx=0;lx<2;lx++){
		B[lx]=new double *[LMAX+1-lx];
		for(ly=0;ly<=LMAX-lx;ly++){
			B[lx][ly]=new double[LMAX+1-lx-ly];
		}
	}
	
	B[0][0][0]=1.0;
	IncrementElement(0,0,0,ir,weight);
	for(L=1;L<=LMAX;L++){
		lfact=chcalc->DoubleFactorial(2*L+1)/chcalc->Factorial(L);
		cL=1.0/double(L);
		cLL=1.0/(L*(2*L-1));
		for(lx=0;lx<2;lx++){
			for(ly=0;ly<=L-lx;ly++){
				lz=L-lx-ly;
				delB=0.0;
				if(lx>0){
					delB+=ex*(cL*lx-cLL*lx*(lx-1))*B[lx-1][ly][lz];
				}
				if(ly>0)
					delB+=ey*(cL*ly-cLL*ly*(ly-1))*B[lx][ly-1][lz];
				if(lz>0)
					delB+=ez*(cL*lz-cLL*lz*(lz-1))*B[lx][ly][lz-1];
				if(lx>1){
					delB-=ey*cLL*lx*(lx-1)*B[lx-2][ly+1][lz];
					delB-=ez*cLL*lx*(lx-1)*B[lx-2][ly][lz+1];
				}
				if(ly>1){
					if(lx>0){
						delB+=ex*cLL*ly*(ly-1)*(B[lx-1][ly-2][lz+2]+B[lx-1][ly][lz]);
					}
					else delB-=ex*cLL*ly*(ly-1)*B[lx+1][ly-2][lz];
					delB-=ez*cLL*ly*(ly-1)*B[lx][ly-2][lz+1];
				}
				if(lz>1){
					if(lx>0){
						delB+=ex*cLL*lz*(lz-1)*(B[lx-1][ly+2][lz-2]+B[lx-1][ly][lz]);
					}
					else delB-=ex*cLL*lz*(lz-1)*B[lx+1][ly][lz-2];
					delB-=ey*cLL*lz*(lz-1)*B[lx][ly+1][lz-2];
				}
				B[lx][ly][lz]=delB;
				if(lx%dlx==0 && ly%dly==0 && lz%dlz==0){
					IncrementElement(lx,ly,lz,ir,delB*lfact*weight);
				}
			}
		}
		//FillRemainderX(B,0);
	}
	for(lx=0;lx<2;lx++){
		for(ly=0;ly<=LMAX-lx;ly++){
			delete [] B[lx][ly];
		}
		delete [] B[lx];
	}
	delete [] B;
}

double CCHArray::GetAExpElementFromMArray(int lx,int ly,int lz,int ir){
	int L,m,mx,my,mz;
	double factor,answer,lfact;
	L=lx+ly+lz;
	lfact=chcalc->DoubleFactorial(2*L+1)/chcalc->Factorial(L);
	answer=0.0;
	for(mx=0;mx<=lx/2;mx++){
		for(my=0;my<=ly/2;my++){
			for(mz=0;mz<=lz/2;mz++){
				m=mx+my+mz;
				if(m>0)	factor=pow(-0.5,m);
				else factor=1.0;
				if(L>m) factor*=chcalc->DoubleFactorial(2*L-2*m-1);
				if(L>0) factor=factor/chcalc->DoubleFactorial(2*L-1);
				factor*=chcalc->Factorial(lx)/(chcalc->Factorial(lx-2*mx)*chcalc->Factorial(mx));
				factor*=chcalc->Factorial(ly)/(chcalc->Factorial(ly-2*my)*chcalc->Factorial(my));
				factor*=chcalc->Factorial(lz)/(chcalc->Factorial(lz-2*mz)*chcalc->Factorial(mz));
				answer+=factor*lfact*GetElement(lx-2*mx,ly-2*my,lz-2*mz,ir);
			}
		}
	}
	return answer;
}

double CCHArray::GetMElementFromAExpArray(int lx,int ly,int lz,int ir){
	int mx,my,mz,m,L;
	double factor,answer,*lfact;
	L=lx+ly+lz;
	lfact=new double[L+1];
	lfact[0]=1.0;
	for(m=1;m<=L;m++) lfact[m]=lfact[m-1]*double(m)/double(2*m+1);
	answer=0.0;
	for(mx=0;mx<=lx/2;mx++){
		for(my=0;my<=ly/2;my++){
			for(mz=0;mz<=lz/2;mz++){
				m=mx+my+mz;
				factor=pow(0.5,m);
				factor=factor*chcalc->DoubleFactorial(2*L-4*m+1)/chcalc->DoubleFactorial(2*L-2*m+1);
				factor*=chcalc->Factorial(lx)/(chcalc->Factorial(lx-2*mx)*chcalc->Factorial(mx));
				factor*=chcalc->Factorial(ly)/(chcalc->Factorial(ly-2*my)*chcalc->Factorial(my));
				factor*=chcalc->Factorial(lz)/(chcalc->Factorial(lz-2*mz)*chcalc->Factorial(mz));
				answer+=factor*lfact[L-2*m]
				*GetElement(lx-2*mx,ly-2*my,lz-2*mz,ir);
			}
		}
	}
	delete lfact;
	return answer;
}	     

void CCHArray::IncrementAExpArrayFromThetaPhi(double theta,double phi,double weight,int ir){
	double stheta,ex,ey,ez;
	stheta=sin(theta);
	ex=stheta*cos(phi);
	ey=stheta*sin(phi);
	ez=cos(theta);
	IncrementAExpArrayFromE(ex,ey,ez,weight,ir);
}

void CCHArray::IncrementMArrayFromE(double ex,double ey,double ez,double weight,int ir){
	int lx,ly,lz;
	for(lx=0;lx<=LMAX;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				IncrementElement(lx,ly,lz,ir,pow(ex,lx)*pow(ey,ly)*pow(ez,lz)*weight);
			}
		}
	}
}

void CCHArray::IncrementMArrayFromThetaPhi(double theta,double phi,double weight,int ir){
	double stheta,ex,ey,ez;
	stheta=sin(theta);
	ex=stheta*cos(phi);
	ey=stheta*sin(phi);
	ez=cos(theta);
	IncrementMArrayFromE(ex,ey,ez,weight,ir);
}

void CCHArray::FillRemainderX(){
	for(int ir=0;ir<NRADIAL;ir++) FillRemainderX(ir);
}

void CCHArray::FillRemainderY(){
	for(int ir=0;ir<NRADIAL;ir++) FillRemainderY(ir);
}

void CCHArray::FillRemainderZ(){
	for(int ir=0;ir<NRADIAL;ir++) FillRemainderZ(ir);
}

void CCHArray::FillRemainderX(int ir){
	int lx,ly,lz;
	if(LMAX>1){
		for(lx=2;lx<=LMAX;lx+=dlx){
			for(ly=0;ly<=LMAX-lx;ly+=dly){
				for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
					SetElement(lx,ly,lz,ir,-GetElement(lx-2,ly+2,lz,ir)-GetElement(lx-2,ly,lz+2,ir));
				}
			}
		}
	}
}

void CCHArray::FillRemainderY(int ir){
	int lx,ly,lz;
	if(LMAX>1){
		for(ly=2;ly<=LMAX;ly+=dly){
			for(lx=0;lx<=LMAX-ly;lx+=dlx){
				for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
					SetElement(lx,ly,lz,ir,-GetElement(lx,ly-2,lz+2,ir)
						-GetElement(lx+2,ly-2,lz,ir));
				}
			}
		}
	}
}

void CCHArray::FillRemainderZ(int ir){
	int lx,ly,lz;
	if(LMAX>1){
		for(lz=2;lz<=LMAX;lz+=dlz){
			for(lx=0;lx<=LMAX-lz;lx+=dlx){
				for(ly=0;ly<=LMAX-lx-lz;ly+=dly){
					SetElement(lx,ly,lz,ir,-GetElement(lx,ly+2,lz-2,ir)
						-GetElement(lx+2,ly,lz-2,ir));
					
				}
			}
		}
	}
}

double CCHArray::AExpand(double theta,double phi,int ir){
	PrintPars();
	double ex,ey,ez,sthet;
	sthet=sin(theta);
	ex=sthet*cos(phi);
	ey=sthet*sin(phi);
	ez=cos(theta);
	return AExpand(ex,ey,ez,ir);
}

double CCHArray::AExpand(double ex,double ey,double ez,int ir){
	double dela,answer=0.0;
	int lx,ly,lz;
	
	for(lx=0;lx<=LMAX;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				dela=GetElement(lx,ly,lz,ir)*pow(ex,lx)*pow(ey,ly)*pow(ez,lz);
				answer+=dela*chcalc->Trinomial(lx,ly,lz);
			}
		}
	}
	return answer;
}

double CCHArray::AExpand(double x,double y,double z){
	double r,ex,ey,ez;
	int ira,irb;
	double w=0.0,interpolate=0.0,valuea=1.0,valueb=0.0;
	/*
	r=sqrt(x*x+y*y+z*z);
	ir=int(floor(r/RADSTEP));
	if(ir<NRADIAL){
		ex=x/r; ey=y/r; ez=z/r;
		answer=AExpand(ex,ey,ez,ir);
	}
	*/
	
	r=sqrt(x*x+y*y+z*z);
	if(r<NRADIAL*RADSTEP){
		irb=int(floor(0.5+r/RADSTEP));
		ira=irb-1;
		ex=x/r; ey=y/r; ez=z/r;
		if(ira>=0 && irb<NRADIAL){
			valuea=AExpand(ex,ey,ez,ira);
			valueb=AExpand(ex,ey,ez,irb);
			w=((0.5+irb)*RADSTEP-r)/RADSTEP;
		}
		else if(ira<0){
			valuea=AExpand(-ex,-ey,-ez,0);
			valueb=AExpand(ex,ey,ez,irb);
			w=((0.5+irb)*RADSTEP-r)/RADSTEP;
		}
		else if(ira==NRADIAL){
			valuea=AExpand(ex,ey,ez,ira);
			valueb=0.0;
			w=((0.5+irb)*RADSTEP-r)/RADSTEP;
		}
		interpolate=w*valuea+(1.0-w)*valueb;
	}
	return interpolate;
}

void CCHArray::RandomInit(int iseed){
	randy=new Crandy(iseed);
}

void CCHArray::Randomize(double mag,int ir){
	int lx,ly,lz;
	double value;
	if(randy==NULL) RandomInit(-12345);
	for(lx=0;lx<=LMAX;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				value=(1.0-2.0*randy->ran())*mag;
				SetElement(lx,ly,lz,ir,value);
			}
		}
	}  
}

void CCHArray::Randomize(double mag){
	int ir;
	if(randy==NULL) RandomInit(-12345);
	for(ir=0;ir<NRADIAL;ir++){
		Randomize(mag,ir);
	}  
}

void CCHArray::RandomizeA(double mag,int ir){
	int lx,ly,lz;
	double value;
	if(randy==NULL) RandomInit(-12345);
	for(lx=0;lx<=1;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				value=(1.0-2.0*randy->ran())*mag;
				SetElement(lx,ly,lz,ir,value);
			}
		}
		FillRemainderX(ir);
	}
}

void CCHArray::RandomizeA(double mag){
	int ir;
	if(randy==NULL) RandomInit(-12345);
	for(ir=0;ir<NRADIAL;ir++){
		RandomizeA(mag,ir);
	}
}

void CCHArray::RandomizeA_Gaussian(double mag,int ir){
	int lx,ly,lz;
	double value;
	if(randy==NULL) RandomInit(-12345);
	for(lx=0;lx<=1;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				value=randy->ran_gauss()*mag;
				SetElement(lx,ly,lz,ir,value);
			}
		}
		FillRemainderX(ir);
	}
}

void CCHArray::RandomizeA_Gaussian(double mag){
	int ir;
	if(randy==NULL) RandomInit(-12345);
	for(ir=0;ir<NRADIAL;ir++){
		RandomizeA_Gaussian(mag,ir);
	}
}

void CCHArray::Detrace(int ir){
	CCHArray *B;
	B=new CCHArray(LMAX,1,RADSTEP,XSYM,YSYM,ZSYM);
	ArrayCalc::Detrace(this,ir,B,0);
	ArrayCalc::CopyArray(B,0,this,ir);
	delete(B);
}

void CCHArray::Detrace(){
	int ir;
	for(ir=0;ir<NRADIAL;ir++) Detrace(ir);
}

void CCHArray::WriteShort(string filename,int LGMAX){
	int ir,lx,ly,lz,L;
	FILE *fptr;
	if(LGMAX>LMAX) LGMAX=LMAX;
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"!       ");
	for(L=0;L<=LGMAX;L++){
		for(lx=0;lx<=1;lx+=dlx){
			for(ly=0;ly<=L-lx;ly+=dly){
				lz=L-lx-ly;
				if(lz%dlz==0)
					fprintf(fptr,"  (%d,%d,%d)  ",lx,ly,lz);
			}
		}
	}
	fprintf(fptr,"\n");
	for(ir=0;ir<NRADIAL;ir++){
		fprintf(fptr,"%7.2f ",RADSTEP*(ir+0.5));
		for(L=0;L<=LGMAX;L++){
			for(lx=0;lx<=1;lx+=dlx){
				for(ly=0;ly<=L-lx;ly+=dly){
					lz=L-lx-ly;
					if(lz%dlz==0)
						fprintf(fptr," %9.2e ",GetElement(lx,ly,lz,ir));
				}
			}
		}
		fprintf(fptr,"\n");
	}
	fclose(fptr);
	
}

void CCHArray::PrintMoments(){
	int ir;
	double r,r2,r3,r4,norm=0.0,m1[3]={0.0},m2[3][3]={{0.0}};
	for(ir=0;ir<NRADIAL;ir++){
		r=(ir+0.5)*RADSTEP;
		r2=r*r;
		r3=r2*r;
		r4=r2*r2;
		norm+=r2*GetMElementFromAExpArray(0,0,0,ir);;
		if(!XSYM) m1[0]+=r3*GetMElementFromAExpArray(1,0,0,ir);
		if(!YSYM) m1[1]+=r3*GetMElementFromAExpArray(0,1,0,ir);
		if(!ZSYM) m1[2]+=r3*GetMElementFromAExpArray(0,0,1,ir);
		m2[0][0]+=r4*GetMElementFromAExpArray(2,0,0,ir);
		if(!XSYM && !YSYM) m2[1][0]+=r4*GetMElementFromAExpArray(1,1,0,ir);
		if(!XSYM && !ZSYM) m2[2][0]+=r4*GetMElementFromAExpArray(1,0,1,ir);
		m2[1][1]+=r4*GetMElementFromAExpArray(0,2,0,ir);
		if(!YSYM && !ZSYM) m2[2][1]+=r4*GetMElementFromAExpArray(0,1,1,ir);
		m2[2][2]+=r4*GetMElementFromAExpArray(0,0,2,ir);
	}
	m1[0]=m1[0]/norm;
	m1[1]=m1[1]/norm;
	m1[2]=m1[2]/norm;
	m2[0][0]=m2[0][0]/norm;
	m2[1][0]=m2[1][0]/norm;
	m2[2][0]=m2[2][0]/norm;
	m2[1][1]=m2[1][1]/norm;
	m2[2][1]=m2[2][1]/norm;
	m2[2][2]=m2[2][2]/norm;
	m2[1][2]=m2[2][1];
	m2[0][2]=m2[2][0];
	m2[0][1]=m2[1][0];
	printf("%10.3e %10.3e %10.3e\n",m2[0][0],m2[0][1],m2[0][2]);
	printf("%10.3e %10.3e %10.3e\n",m2[1][0],m2[1][1],m2[1][2]);
	printf("%10.3e %10.3e %10.3e\n",m2[2][0],m2[2][1],m2[2][2]);
}
