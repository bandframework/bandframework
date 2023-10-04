#include "msu_commonutils/sf.h"

using namespace std;

int CCHCalc::INITIALIZED=0;
int CCHCalc::LMAXFACT=56;
double *CCHCalc::fact=NULL;
double *CCHCalc::doublefact=NULL;
double **CCHCalc::binomial=NULL;
double *****CCHCalc::overlap=NULL;

CCHCalc::CCHCalc(){
  if(INITIALIZED==0){
    InitStaticData();
  }
  INITIALIZED+=1;
};

CCHCalc::~CCHCalc(){
  INITIALIZED-=1;
  if(INITIALIZED==0){
    ClearStaticData();
    printf("CCHCalc static data cleared\n");
  }
}

void CCHCalc::InitStaticData(){
  int lx,ly,lz,lymax,lzmax,L;
  int m;
  if(fact==NULL){
    fact=new double[LMAXFACT+1];
    fact[0]=fact[1]=1.0;
    for(m=2;m<=LMAXFACT;m++) fact[m]=fact[m-1]*double(m);
  }
  if(doublefact==NULL){
    doublefact=new double[2*LMAXFACT+2];
    doublefact[0]=doublefact[1]=1.0;
    for(m=2;m<=2*LMAXFACT+1;m++) doublefact[m]=doublefact[m-2]*double(m);  
  }

  if(binomial==NULL){
    binomial=new double *[LMAXFACT+1];
    for(L=0;L<=LMAXFACT;L++){
      binomial[L]=new double [L+1];
      binomial[L][0]=1.0;
      for(lx=1;lx<=L;lx++)
	binomial[L][lx]=binomial[L][lx-1]*double(L-lx+1)/double(lx);
    }
  }
  
  if(overlap==NULL){
    overlap=new double ****[LMAXFACT+1];
    for(lx=0;lx<=LMAXFACT;lx++){
      overlap[lx]=new double ***[LMAXFACT-lx+1];
      for(ly=0;ly<=LMAXFACT-lx;ly++){
	for(lz=0;lz<=LMAXFACT-lx-ly;lz++){
	  L=lx+ly+lz;
	}
      }
      lymax=lx;
      if(lymax+lx>LMAXFACT) lymax=LMAXFACT-lx;
      for(ly=0;ly<=lymax;ly++){
	overlap[lx][ly]=new double **[LMAXFACT-lx-ly+1];
	lzmax=ly;
	if(lx+ly+lzmax>LMAXFACT) lzmax=LMAXFACT-lx-ly;
	for(lz=0;lz<=lzmax;lz++){
	  overlap[lx][ly][lz]=NULL;
	  //overlapinit(lx,ly,lz);
	}
      }
    }
  }
}

void CCHCalc::ClearStaticData(){
  int lx,ly,lz,lymax,lzmax,mx,L;
  printf("Clearing static data used for charray calc.s. \nDon't do this unless you are finished with charray objects\n");
  delete [] fact;
  delete [] doublefact;
  for(L=0;L<=LMAXFACT;L++) delete [] binomial[L];
  delete [] binomial;
  for(lx=0;lx<=LMAXFACT;lx++){
    lymax=LMAXFACT-lx;
    if(lymax>lx) lymax=lx;
    for(ly=0;ly<=lymax;ly++){
      lzmax=LMAXFACT-lx-ly;
      if(lzmax>ly) lzmax=ly;
      for(lz=0;lz<=lzmax;lz++){
	if(overlap[lx][ly][lz]!=NULL){
	  L=lx+ly+lz;
	  for(mx=lx%2;mx<=L/2;mx++){
	    delete [] overlap[lx][ly][lz][mx];
	  }
	}
	delete [] overlap[lx][ly][lz];
      }
      delete [] overlap[lx][ly];
    }
    delete [] overlap[lx];
  }
  delete overlap;
  overlap=NULL;
  doublefact=NULL;
  fact=NULL;
  binomial=NULL;
}

double CCHCalc::Factorial(int n){
  return fact[n];
}

double CCHCalc::DoubleFactorial(int n){
  return doublefact[n];
}

double CCHCalc::Binomial(int lx,int ly){
  return binomial[lx+ly][lx];
}

double CCHCalc::Trinomial(int lx,int ly,int lz){
  int L;
  L=lx+ly+lz;
  if(lx<ly) iswitch(lx,ly);
  if(ly<lz) iswitch(ly,lz);
  if(lx<ly) iswitch(lx,ly);
  return binomial[L][lz]*binomial[lx+ly][lx];
}

void CCHCalc::iswitch(int &i,int &j){
  int k;
  k=i;
  i=j;
  j=k;
}

double CCHCalc::GetOverlap(int lx,int ly,int lz,
			     int lxprime,int lyprime,int lzprime){
  int L,Lprime,mx,my;
  L=lx+ly+lz;
  Lprime=lxprime+lyprime+lzprime;
 
  if(L!=Lprime){
    return 0.0;
  }
  else if(((lx-lxprime)%2)!=0){
    return 0.0;
  }
  else if(((ly-lyprime)%2)!=0){
    return 0.0;
  }
  else if(((lz-lzprime)%2)!=0){
    return 0.0;
  }
  else if(L<=LMAXFACT){
    if(ly>lx){
      iswitch(lx,ly);
      iswitch(lxprime,lyprime);
    }
    if(lz>ly){
      iswitch(ly,lz);
      iswitch(lyprime,lzprime);
      if(ly>lx){
	iswitch(lx,ly);
	iswitch(lxprime,lyprime);
      }
    }
    if(lx<ly || ly<lz) {
      printf("l out of order\n");
      exit(1);
    }
    if(overlap[lx][ly][lz]==NULL) overlapinit(lx,ly,lz);
    mx=lxprime/2;
    my=lyprime/2;
    //printf("check: l=(%d,%d,%d), lprime=(%d,%d,%d), mx=%d, my=%d\n",
    //   lx,ly,lz,lxprime,lyprime,lzprime,mx,my);
    //printf("overlap=%g\n",overlap[lx][ly][lz][mx][my]);
    //printf("________________\n");
    return overlap[lx][ly][lz][mx][my];
  }
  else{
    return GetOverlap0(lx,ly,lz,lxprime,lyprime,lzprime);
  }
}

void CCHCalc::overlapinit(int lx,int ly,int lz){
  int lxprime,lyprime,lzprime,L,mx,my;
  L=lx+ly+lz;
  if(overlap[lx][ly][lz]==NULL){
    overlap[lx][ly][lz]=new double *[1+L/2];
    for(mx=0;mx<=L/2;mx++){
      lxprime=lx%2+2*mx;
      overlap[lx][ly][lz][mx]=new double[1+(L-lxprime)/2];
      for(my=0;my<=(L-lxprime)/2;my++){
	lyprime=ly%2+2*my;
	lzprime=L-lxprime-lyprime;
	//printf("___l=(%d,%d,%d), lprime=(%d,%d,%d), ",
	//   lx,ly,lz,lxprime,lyprime,lzprime);
	overlap[lx][ly][lz][mx][my]
	  =GetOverlap0(lx,ly,lz,lxprime,lyprime,lzprime);
	//printf("overlap=%g\n",overlap[lx][ly][lz][lxprime][lyprime]);
	
      }
       }
  }
}


double CCHCalc::GetMFromE(int lx,int ly,int lz,
			    double ex,double ey,double ez){
  return pow(ex,lx)*pow(ey,ly)*pow(ez,lz);
}

double CCHCalc::GetMFromThetaPhi(int lx,int ly,int lz,
				   double theta,double phi){
  double stheta,ex,ey,ez;
  stheta=sin(theta);
  ex=stheta*cos(phi);
  ey=stheta*sin(phi);
  ez=cos(theta);
  return GetMFromE(lx,ly,lz,ex,ey,ez);
}

double CCHCalc::GetAFromE(int lx,int ly,int lz,
			       double ex,double ey,double ez){
  int L,m,mx,my,mz;
  double answer,factor,fx,fy,fz,fl;
  L=lx+ly+lz;
  if(L==0) return 1.0;
  else{
    answer=0.0;
    fl=fact[lx]*fact[ly]*fact[lz]/doublefact[2*L-1];
    for(mx=0;mx<=lx/2;mx++){
      fx=pow(-0.5,mx)*pow(ex,lx-2*mx)/(fact[mx]*fact[lx-2*mx]);
      for(my=0;my<=ly/2;my++){
	fy=pow(-0.5,my)*pow(ey,ly-2*my)/(fact[my]*fact[ly-2*my]);
	for(mz=0;mz<=lz/2;mz++){
	  fz=pow(-0.5,mz)*pow(ez,lz-2*mz)/(fact[mz]*fact[lz-2*mz]);
	  m=mx+my+mz;
	  factor=fl*fx*fy*fz;
	  if(L>m) factor*=doublefact[2*L-2*m-1];
	  answer+=factor;
	}
      }
    }
    return answer;
  }
}

double CCHCalc::GetAFromThetaPhi(int lx,int ly,int lz,double theta,double phi){
  double stheta,ex,ey,ez;
  stheta=sin(theta);
  ex=stheta*cos(phi);
  ey=stheta*sin(phi);
  ez=cos(theta);
  return GetAFromE(lx,ly,lz,ex,ey,ez);
}


double CCHCalc::GetOverlap0(int lx,int ly,int lz,
				int lxprime,int lyprime,int lzprime){
  const double PI=4.0*atan(1.0);
  int L,Lprime;
  double prefactor,factor,sum;
  int mxmin,mymin,mzmin,m,mx,my,mz,kx,ky,kz;
  L=lx+ly+lz;
  Lprime=lxprime+lyprime+lzprime;
  if(L!=Lprime){
    return 0.0;
  }
  else if(((lx-lxprime)%2)!=0){
    return 0.0;
  }
  else if(((ly-lyprime)%2)!=0){
    return 0.0;
  }
  else if(((lz-lzprime)%2)!=0){
    return 0.0;
  }
  else{
    prefactor=4.0*PI*fact[lx]*fact[ly]*fact[lz]*fact[lxprime]*fact[lyprime]
      *fact[lzprime]/((2.0*double(L)+1.0));
    if(L>0)
      prefactor=prefactor/pow(doublefact[2*L-1],2);
    sum=0.0;
    mxmin=mymin=mzmin=0;
    if(lx>lxprime) mxmin=(lx-lxprime)/2;
    if(ly>lyprime) mymin=(ly-lyprime)/2;
    if(lz>lzprime) mzmin=(lz-lzprime)/2;
    for(mx=mxmin;mx<=lx/2;mx++){
      for(my=mymin;my<=ly/2;my++){
	for(mz=mzmin;mz<=lz/2;mz++){
	  m=mx+my+mz;
	  if(2*L-2*m-1>0) factor=doublefact[2*L-2*m-1]*fact[m];
	  else factor=1.0;
	  if(m>0) factor*=pow(-0.5,m);
	  kx=mx+(lxprime-lx)/2;
	  ky=my+(lyprime-ly)/2;
	  kz=mz+(lzprime-lz)/2;
	  factor=factor/double(fact[mx]*fact[kx]*fact[lx-2*mx]);
	  factor=factor/double(fact[my]*fact[ky]*fact[ly-2*my]);
	  factor=factor/double(fact[mz]*fact[kz]*fact[lz-2*mz]);
	  sum+=factor;
	}
      }
    }
    return sum*prefactor;
  }
}
