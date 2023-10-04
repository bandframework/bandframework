#include "msu_commonutils/arrays.h"

CMCList::CMCList(int nmcset){
  int i;
  nmc=nmcset;
  r=new double *[nmc];
  for(i=0;i<nmc;i++) r[i]=new double[4];
	norm=1.0;
}

CMCList::~CMCList(){
  int i;
  for(i=0;i<nmc;i++) delete [] r[i];
  delete [] r;
}

void CMCList::Resize(int nmcset){
  int i;
  for(i=0;i<nmc;i++) delete [] r[i];
  delete [] r;
  nmc=nmcset;
  r=new double *[nmc];
  for(i=0;i<nmc;i++) r[i]=new double[4];  
}

int CMCList::GetNMC(){
  return nmc;
}

void CMCList::SetR(int imc,double *rr){
  int alpha;
  if(imc>=nmc){
    printf("trying to set CMCList.nmc too big, imc=%d, nmc=%d\n",imc,nmc);
    exit(1);
  }
  for(alpha=0;alpha<4;alpha++){
    r[imc][alpha]=rr[alpha];
  }
}

void CMCList::SetR(int imc,double t,double x,double y,double z){
  if(imc>=nmc){
    printf("trying to set CMCList.nmc too big, imc=%d, nmc=%d\n",imc,nmc);
    exit(1);
  }
  r[imc][0]=t;
  r[imc][1]=x;
  r[imc][2]=y;
  r[imc][3]=z;
}

double *CMCList::GetR(int imc){
  if(imc>=nmc){
    printf("trying to set CMCList.nmc too big, imc=%d, nmc=%d\n",imc,nmc);
    exit(1);
  }
  return r[imc];
}

void CMCList::SetNorm(double normset){
  norm=normset;
}

double CMCList::GetNorm(){
  return norm;
}

void CMCList::PrintMoments(double Rmax){
	double r2,xbar[4]={0.0},x2bar[4][4]={{0.0}};
	int imc,alpha,beta;
	norm=0;
	for(imc=0;imc<nmc;imc++){
		r2=r[imc][1]*r[imc][1]+r[imc][2]*r[imc][2]+r[imc][3]*r[imc][3];
		if(r2<Rmax*Rmax){
			for(alpha=0;alpha<4;alpha++){
				xbar[alpha]+=r[imc][alpha];
				for(beta=0;beta<4;beta++){
					x2bar[alpha][beta]+=r[imc][alpha]*r[imc][beta];
				}
				norm+=1;
			}
		}
	}
	printf("--------- CMCList Moments ----------\n");
	printf("<x[alpha]>= ");
	for(alpha=0;alpha<4;alpha++){
		xbar[alpha]=xbar[alpha]/norm;
		printf("%10.3e ",xbar[alpha]);
	}
	printf("\n");
	printf("sigma^2_alpha,beta= \n");
	for(alpha=0;alpha<4;alpha++){
		for(beta=0;beta<4;beta++){
			x2bar[alpha][beta]=x2bar[alpha][beta]/norm;
			x2bar[alpha][beta]-=xbar[alpha]*xbar[beta];
			printf("%10.3e ",x2bar[alpha][beta]);
		}
		printf("\n");
	}
}
