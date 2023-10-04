#include <cstdlib>
#include <algorithm>
#include "msu_commonutils/decay_nbody.h"
#include "msu_commonutils/log.h"
using namespace std;

void CDecay_NBody::SetMasses2(double Mset,double m1set,double m2set){
	M=Mset; m1=m1set; m2=m2set;
	nbodies=2;
}

void CDecay_NBody::GenerateMomenta2(FourVector &p1,FourVector &p2){
	Qmag=sqrt(Misc::triangle(M,m1,m2));
	double ctheta,stheta,phi;
	ctheta=1.0-2.0*randy->ran();
	stheta=sqrt(1.0-ctheta*ctheta);
	phi=2.0*PI*randy->ran();
	p1[3]=Qmag*ctheta;
	p1[1]=Qmag*stheta*cos(phi);
	p1[2]=Qmag*stheta*sin(phi);
	p1[0]=sqrt(Qmag*Qmag+m1*m1);
	p2[1]=-p1[1];
	p2[2]=-p1[2];
	p2[3]=-p1[3];
	p2[0]=sqrt(Qmag*Qmag+m2*m2);
}


void CDecay_NBody::SetMasses3(double Mset,double m1set,double m2set,double m3set){
	const double fracguess=0.5;
	maxfactor3=1.15;  // empircally determined
	double w;
	nbodies=3;
	M=Mset; m1=m1set; m2=m2set; m3=m3set;
	wmax3=1.0;
	KEtot=M-m1-m2-m3;
	M12=m1+m2+fracguess*KEtot;
	w=GetW3();
	wmax3=w*maxfactor3;
}

void CDecay_NBody::SetMasses4(double Mset,double m1set,double m2set,double m3set,double m4set){
	const double fracguess=0.35;
	maxfactor4=1.27;  // empircally determined
	double w;
	nbodies=4;
	M=Mset; m1=m1set; m2=m2set; m3=m3set; m4=m4set;
	wmax4=1.0;
	KEtot=(M-m1-m2-m3-m4);
	M12=m1+m2+fracguess*KEtot;
	M34=m3+m4+fracguess*KEtot;
	w=GetW4();
	wmax4=w*maxfactor4;
}

void CDecay_NBody::SetMasses(int nbodies_set,vector<double> &masses_set){
	nbodies=nbodies_set;
	if(nbodies==2){
		SetMasses2(masses_set[0],masses_set[1],masses_set[2]);
	}
	else if(nbodies==3){
		SetMasses3(masses_set[0],masses_set[1],masses_set[2],masses_set[3]);
	}
	else if(nbodies==4){
		SetMasses4(masses_set[0],masses_set[1],masses_set[2],masses_set[3],masses_set[4]);
	}
	else{
		double w;
		int i;
		allmassless=true;
		if(int(masses.size())<=nbodies){
			masses.resize(nbodies+1);
		}
		for(i=0;i<=nbodies;i++){
			masses[i]=masses_set[i];
			if(i>0 && masses[i]>1.0E-10)
				allmassless=false;
		}
		M=masses[0];
		//masses=masses_set;
		//maxfactor=1.0;
		maxfactor=1.0+0.028*(nbodies-2.0)+0.41*tanh((nbodies-2.0)/2.8);
		if(nbodies>int(Msum.size())){
			Msum.resize(nbodies);
		}
		if(nbodies>int(qsum.size())){
			qsum.resize(nbodies);
		}
	
		KEtot=masses[0];
		for(i=1;i<=nbodies;i++){
			KEtot-=masses[i];
		}
		if(KEtot<0.0){
			for(i=0;i<=nbodies;i++)
				CLog::Info(to_string(masses[i])+" ");
			CLog::Info("\n");
			CLog::Fatal("masses don't add up in CDecay_NBody::SetMasses, KEtot="+to_string(KEtot)+"\n");
		}
		Msum[0]=masses[1];
		for(i=1;i<nbodies;i++){
			Msum[i]=Msum[i-1]+masses[i+1]+KEtot/double(nbodies-1);
		}
		wmax=1.0;
		qbar=sqrt(KEtot*masses[0])/double(nbodies);
		w=GetW();
		wmax=w*maxfactor;
	}
}

// this is experimental
void CDecay_NBody::SetMasses_Trial(int nbodies_set,vector<double> &masses_set){
	double w;
	int i;
	//vector<double> maxfactor={1.0,1.0,1.0,1.03,1.04,1.05,1.06,1.06,}
	nbodies=nbodies_set;
	maxfactor=1.0;
	qsum.resize(nbodies);
	if(int(masses.size())<=nbodies){
		masses.resize(nbodies+1);
	}
	for(i=0;i<=nbodies;i++)
		masses[i]=masses_set[i];
	Msum.resize(nbodies);
	KEtot=masses[0];
	for(i=1;i<=nbodies;i++){
		KEtot-=masses[i];
	}
	if(KEtot<0.0){
		for(i=0;i<=nbodies;i++)
			CLog::Info(to_string(masses[i])+" ");
		CLog::Info("\n");
		CLog::Fatal("masses don't add up in CDecay_NBody::SetMasses, KEtot="+to_string(KEtot)+"\n");
	}
	
	double mu,msum,KEsum,epsilon;
	vector<double> KE(nbodies);
	KE.resize(nbodies);
	Msum[0]=masses[1];
	KEsum=0.0;
	epsilon=KEtot/(nbodies-1);
	msum=masses[1];
	for(i=1;i<nbodies;i++){
		if(msum>1.0E-9 || masses[i+1]>1.0E-9){
			mu=msum*masses[i+1]/(msum+masses[i+1]);
			//KE[i]=1.0+epsilon*epsilon/(mu*mu+epsilon*epsilon);
			KE[i]=1.0+1.5*epsilon/(mu+epsilon);
			msum+=masses[i+1]+epsilon;
		}
		else{
			KE[i]=2.0;
		}
		KEsum+=KE[i];
	}
	for(i=1;i<nbodies;i++){
		KE[i]=KEtot*KE[i]/KEsum;
		Msum[i]=Msum[i-1]+masses[i+1]+KE[i];
	}
	
	wmax=1.0;
	w=GetW();
	qbar=sqrt(KEtot*masses[0])/double(nbodies);
	wmax=w*maxfactor/qbar;
	KE.clear();
}

void CDecay_NBody::ChooseM12(){
	double KEtotmax,KE12,w;
	KEtotmax=M-m1-m2-m3;
	do{
		KE12=KEtotmax*randy->ran();
		M12=m1+m2+KE12;
		w=GetW3();
		if(w>wmaxmax3)
			wmaxmax3=w;
		if(w>1.0){
			CLog::Info("YIKES!!! w="+to_string(w)+" exceeded 1.0, may need to increase wmaxfactor3\n");
		}
		Ntry+=1;
	}while(w<randy->ran());
	Nsuccess+=1;
}

void CDecay_NBody::ChooseM12M34(){
	double KEtotmax,KE12,KE34,w;
	KEtotmax=M-m1-m2-m3-m4;
	do{
		do{
			KE12=KEtotmax*randy->ran();
			KE34=KEtotmax*randy->ran();
		}while(KE12+KE34>KEtotmax);
		M12=m1+m2+KE12;
		M34=m3+m4+KE34;
		w=GetW4();
		if(w>wmaxmax4)
			wmaxmax4=w;
		if(w>1.0){
			CLog::Info("YIKES!!! w="+to_string(w)+" exceeded 1.0, may need to increase wmaxfactor4\n");
		}
		Ntry+=1;
	}while(w<randy->ran());
	Nsuccess+=1;
}

void CDecay_NBody::ChooseMsum(){
	int i;
	double w;
	vector<double> KE;
	KE.resize(nbodies);
	vector<double> x;
	x.resize(nbodies);
	x[0]=0.0;
	do{
		for(i=1;i<nbodies-1;i++)
			x[i]=KEtot*randy->ran();
		x[nbodies-1]=KEtot;
		sort(x.begin(),x.end());
		for(i=0;i<nbodies-1;i++){
			KE[i]=x[i+1]-x[i];
		}
	
		Msum[0]=masses[1];
		for(i=1;i<nbodies;i++){
			Msum[i]=Msum[i-1]+masses[i+1]+KE[i-1];
		}
		w=GetW();
		if(w>wmaxmax)
			wmaxmax=w;
		Ntry+=1;
	}while(w<randy->ran());
	Nsuccess+=1;
	KE.clear();
	x.clear();
}

void CDecay_NBody::CompleteMomenta3(FourVector &p1,FourVector &p2,FourVector &p3){
	FourVector u,Q;
	double ctheta,stheta,phi;
	int alpha;
	
	ctheta=1.0-2.0*randy->ran();
	stheta=sqrt(1.0-ctheta*ctheta);
	phi=2.0*PI*randy->ran();
	p1[3]=q12*ctheta;
	p1[1]=q12*stheta*cos(phi);
	p1[2]=q12*stheta*sin(phi);
	p1[0]=sqrt(p1[1]*p1[1]+p1[2]*p1[2]+p1[3]*p1[3]+m1*m1);
	p2[3]=-p1[3];
	p2[1]=-p1[1];
	p2[2]=-p1[2];
	p2[0]=sqrt(p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3]+m2*m2);
	
	
	ctheta=1.0-2.0*randy->ran();
	stheta=sqrt(1.0-ctheta*ctheta);
	phi=2.0*PI*randy->ran();
	Q[3]=Qmag*ctheta;
	Q[1]=Qmag*stheta*cos(phi);
	Q[2]=Qmag*stheta*sin(phi);
	Q[0]=sqrt(Qmag*Qmag+M12*M12);
	for(alpha=0;alpha<4;alpha++)
		u[alpha]=Q[alpha]/M12;
	Boost(u,p1);
	Boost(u,p2);
	
	p3[0]=sqrt(Qmag*Qmag+m3*m3);
	for(alpha=1;alpha<4;alpha++){
		p3[alpha]=-Q[alpha];
	}
}

void CDecay_NBody::CompleteMomenta4(FourVector &p1,FourVector &p2,FourVector &p3,FourVector &p4){
	FourVector u,Q;
	double ctheta,stheta,phi;
	int alpha;
	
	ctheta=1.0-2.0*randy->ran();
	stheta=sqrt(1.0-ctheta*ctheta);
	phi=2.0*PI*randy->ran();
	p1[3]=q12*ctheta;
	p1[1]=q12*stheta*cos(phi);
	p1[2]=q12*stheta*sin(phi);
	p1[0]=sqrt(p1[1]*p1[1]+p1[2]*p1[2]+p1[3]*p1[3]+m1*m1);
	p2[3]=-p1[3];
	p2[1]=-p1[1];
	p2[2]=-p1[2];
	p2[0]=sqrt(p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3]+m2*m2);
	
	ctheta=1.0-2.0*randy->ran();
	stheta=sqrt(1.0-ctheta*ctheta);
	phi=2.0*PI*randy->ran();
	p3[3]=q34*ctheta;
	p3[1]=q34*stheta*cos(phi);
	p3[2]=q34*stheta*sin(phi);
	p3[0]=sqrt(p3[1]*p3[1]+p3[2]*p3[2]+p3[3]*p3[3]+m3*m3);
	p4[3]=-p3[3];
	p4[1]=-p3[1];
	p4[2]=-p3[2];
	p4[0]=sqrt(p4[1]*p4[1]+p4[2]*p4[2]+p4[3]*p4[3]+m4*m4);
	
	ctheta=1.0-2.0*randy->ran();
	stheta=sqrt(1.0-ctheta*ctheta);
	phi=2.0*PI*randy->ran();
	Q[3]=Qmag*ctheta;
	Q[1]=Qmag*stheta*cos(phi);
	Q[2]=Qmag*stheta*sin(phi);
	Q[0]=sqrt(Qmag*Qmag+M12*M12);
	for(alpha=0;alpha<4;alpha++)
		u[alpha]=Q[alpha]/M12;
	Boost(u,p1);
	Boost(u,p2);
	
	Q[3]=-Q[3];
	Q[1]=-Q[1];
	Q[2]=-Q[2];
	Q[0]=sqrt(Qmag*Qmag+M34*M34);
	for(alpha=0;alpha<4;alpha++)
		u[alpha]=Q[alpha]/M34;
	Boost(u,p3);
	Boost(u,p4);
}

double CDecay_NBody::GetW3(){
	double w;
	q12=sqrt(Misc::triangle(M12,m1,m2));
	Qmag=sqrt(Misc::triangle(M,M12,m3));
	w=q12*Qmag/wmax3;
	return w;
}

double CDecay_NBody::GetW4(){
	double w;
	q12=sqrt(Misc::triangle(M12,m1,m2));
	q34=sqrt(Misc::triangle(M34,m3,m4));
	Qmag=sqrt(Misc::triangle(M,M12,M34));
	w=q12*q34*Qmag/wmax4;
	return w;
}

double CDecay_NBody::GetW(){
	int i;
	double w=1.0;
	for(i=1;i<nbodies;i++){
		qsum[i]=sqrt(Misc::triangle(Msum[i],Msum[i-1],masses[i+1]));
		w*=qsum[i]/qbar;
	}
	w=w/wmax;
	return w;
}

void CDecay_NBody::GenerateMomenta3(FourVector &p1,FourVector &p2,FourVector &p3){
	ChooseM12();
	CompleteMomenta3(p1,p2,p3);
}

void CDecay_NBody::GenerateMomenta4(FourVector &p1,FourVector &p2,FourVector &p3,FourVector &p4){
	ChooseM12M34();
	CompleteMomenta4(p1,p2,p3,p4);
}

void CDecay_NBody::GenerateMomenta(vector<FourVector> &p){
	if(int(p.size())<nbodies){
		CLog::Fatal("danger: not enough fourvectors in array for CDecay_NBody\n");
	}
	if(nbodies==2){
		GenerateMomenta2(p[0],p[1]);
	}
	else if(nbodies==3){
		GenerateMomenta3(p[0],p[1],p[2]);
	}
	else if(nbodies==4){
		GenerateMomenta4(p[0],p[1],p[2],p[3]);
	}
	else{
		if(allmassless){
			GenerateMomentaAllMassless(p);
		}
		else{
			ChooseMsum();
			CompleteMomenta(p);
		}
	}
}

void CDecay_NBody::CompleteMomenta(vector<FourVector> &p){
	FourVector u;
	int i,j;
	double ctheta,stheta,phi;
	int alpha;
	p[0][0]=masses[1];
	p[0][1]=p[0][2]=p[0][3]=0.0;
	for(i=1;i<nbodies;i++){
		ctheta=1.0-2.0*randy->ran();
		stheta=sqrt(1.0-ctheta*ctheta);
		phi=2.0*PI*randy->ran();
		
		p[i][3]=qsum[i]*ctheta;
		p[i][1]=qsum[i]*stheta*cos(phi);
		p[i][2]=qsum[i]*stheta*sin(phi);
		p[i][0]=sqrt(p[i][1]*p[i][1]+p[i][2]*p[i][2]+p[i][3]*p[i][3]+masses[i+1]*masses[i+1]);
		
		u[0]=1.0;
		for(alpha=1;alpha<4;alpha++){
			u[alpha]=-p[i][alpha]/Msum[i-1];
			u[0]+=u[alpha]*u[alpha];
		}
		u[0]=sqrt(u[0]);
		for(j=0;j<i;j++){
			Misc::Boost(u,p[j]);
		}

	}
}

void CDecay_NBody::GenerateMomentaAllMassless(vector<FourVector> &p){
	int ibody,alpha;
	double T=100.0;
	double M=masses[0],Mtot,factor;
	FourVector ptot,u;
	for(alpha=0;alpha<4;alpha++)
		ptot[alpha]=0.0;
	for(ibody=0;ibody<nbodies;ibody++){
		randy->generate_boltzmann(0.0,T,p[ibody]);
		for(alpha=0;alpha<4;alpha++)
			ptot[alpha]+=p[ibody][alpha];
	}
	Mtot=sqrt(ptot[0]*ptot[0]-ptot[1]*ptot[1]-ptot[2]*ptot[2]-ptot[3]*ptot[3]);
	for(alpha=0;alpha<4;alpha++)
		u[alpha]=ptot[alpha]/Mtot;
	for(alpha=0;alpha<4;alpha++)
		ptot[alpha]=0.0;
	for(ibody=0;ibody<nbodies;ibody++){
		Misc::BoostToCM(u,p[ibody],p[ibody]);
		for(alpha=0;alpha<4;alpha++)
			ptot[alpha]+=p[ibody][alpha];
	}
	factor=M/ptot[0];
	for(ibody=0;ibody<nbodies;ibody++){
		for(alpha=0;alpha<4;alpha++)
			p[ibody][alpha]=p[ibody][alpha]*factor;
	}
}
	
	
	
