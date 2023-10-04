#include "msu_commonutils/misc.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/sf.h"

void Misc::lorentz(FourVector &u,FourVector &p,FourVector &pprime){
	int alpha;
	double pdotu=p[0]*u[0]-p[1]*u[1]-p[2]*u[2]-p[3]*u[3];
	double A=-(pdotu+p[0])/(1.0+u[0]);
	for(alpha=0;alpha<4;alpha++)
		pprime[alpha]=p[alpha];
	pprime[0]+=A;
	for(alpha=0;alpha<4;alpha++)
		pprime[alpha]+=(A+2*p[0])*u[alpha];
}

void Misc::Boost(FourVector &u,FourVector &p,FourVector &pprime){
	int alpha;
	double pdotu=p[0]*u[0]-p[1]*u[1]-p[2]*u[2]-p[3]*u[3];
	double A=-(pdotu+p[0])/(1.0+u[0]);
	for(alpha=0;alpha<4;alpha++)
		pprime[alpha]=p[alpha];
	pprime[0]+=A;
	for(alpha=0;alpha<4;alpha++)
		pprime[alpha]+=(A+2*p[0])*u[alpha];
}

void Misc::BoostToCM(FourVector &u,FourVector &p,FourVector &ptilde){
	int alpha;
	double pdotu=p[0]*u[0]-p[1]*u[1]-p[2]*u[2]-p[3]*u[3];
	double A=-(pdotu+p[0])/(1.0+u[0]);
	for(alpha=0;alpha<4;alpha++)
		ptilde[alpha]=p[alpha];
	ptilde[0]+=A+2*pdotu;
	for(alpha=0;alpha<4;alpha++)
		ptilde[alpha]+=A*u[alpha];
}

void Misc::BoostToCM(FourVector &u,double **Pi,double **PiTilde){
	int alpha,beta,gamma,delta;
	double Linv[4][4],L[4][4];  // Lorentz Boost Matrices
	double ucontra[4]={u[0],-u[1],-u[2],-u[3]},n[4]={1.0,0,0,0};
	for(alpha=0;alpha<4;alpha++){
		for(beta=0;beta<4;beta++){
			Linv[alpha][beta]=2.0*n[alpha]*ucontra[beta]-(u[alpha]+n[alpha])*(ucontra[beta]+n[beta])/(1.0+u[0]);
			if(alpha==beta) Linv[alpha][beta]+=1.0;
			L[beta][alpha]=Linv[alpha][beta];
		}
	}
	for(alpha=0;alpha<4;alpha++){
		for(delta=0;delta<4;delta++){
			PiTilde[alpha][delta]=0.0;
			for(beta=0;beta<4;beta++){
				for(gamma=0;gamma<4;gamma++){
					PiTilde[alpha][delta]+=Linv[alpha][beta]*Pi[beta][gamma]*L[gamma][delta];
				}
			}
		}
	}
}

void Misc::Boost(FourVector &u,double **PiTilde,double **Pi){
	int alpha,beta,gamma,delta;
	double Linv[4][4],L[4][4];  // Lorentz Boost Matrices
	double ucontra[4]={u[0],-u[1],-u[2],-u[3]},n[4]={1.0,0,0,0};
	for(alpha=0;alpha<4;alpha++){
		for(beta=0;beta<4;beta++){
			Linv[alpha][beta]=2.0*n[alpha]*u[beta]-(ucontra[alpha]+n[alpha])*(u[beta]+n[beta])/(1.0+u[0]);
			if(alpha==beta) Linv[alpha][beta]+=1.0;
			L[beta][alpha]=Linv[alpha][beta];
		}
	}
	for(alpha=0;alpha<4;alpha++){
		for(delta=0;delta<4;delta++){
			Pi[alpha][delta]=0.0;
			for(beta=0;beta<4;beta++){
				for(gamma=0;gamma<4;gamma++){
					Pi[alpha][delta]+=L[alpha][beta]*PiTilde[beta][gamma]*Linv[gamma][delta];
				}
			}
		}
	}	
}

void Misc::Boost(FourVector &u,FourVector &p){
	int mu;
	double ndotp,udotn,udotp;
	ndotp=p[0];
	udotn=u[0];
	udotp=u[0]*p[0]-u[1]*p[1]-u[2]*p[2]-u[3]*p[3];
	p[0]=p[0]-((udotp+ndotp)/(1.0+udotn));
	for(mu=0;mu<4;mu++){
		p[mu]=-((udotp+ndotp)/(1.0+udotn))*u[mu]+2*ndotp*u[mu]+p[mu];
	}
}
