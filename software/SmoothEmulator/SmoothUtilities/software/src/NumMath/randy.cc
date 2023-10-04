#include "msu_commonutils/randy.h"
#include "msu_commonutils/constants.h"

using namespace std;

Crandy::Crandy(int iseed){
	seed=iseed;
	mt.seed(seed);
	netprob=0.0;
	threshold=-log(ran());
}

void Crandy::reset(int iseed){
	seed=iseed;
	mt.seed(seed);
	netprob=0.0;
	threshold=-log(ran());
}

double Crandy::ran(){
	return ranu(mt);
}

double Crandy::ran_exp(){
	return -log(ran());
}

double Crandy::ran_gauss(void){
	return rang(mt);
}

void Crandy::ran_gauss2(double &ra,double &rb){
	double x,y,r2,r,c,s;
	TRY_AGAIN:
	x=1.0-2.0*ran();
	y=1.0-2.0*ran();
	r2=x*x+y*y;
	if(r2>1.0) goto TRY_AGAIN;
	r=sqrt(r2);
	c=x/r;
	s=y/r;
	ra=c*sqrt(-2.0*log(r2));
	rb=(s/c)*ra;
}

void Crandy::generate_boltzmann_alt(double mass,double T,FourVector &p){
	//const double PI=4.0*atan(1.0);
	double r1,r2,r3,r0,I1,I2,I3,Itot;
	double pmag,E,ctheta,stheta,phi,K;
	array<double,4> pp;
	pp[0]=pp[1]=pp[2]=pp[3]=1.2345;
	GB_TRYAGAIN:
	r0=ran();
	I1=mass*mass;
	I2=2.0*mass*T;
	I3=2.0*T*T;
	Itot=I1+I2+I3;
	if(r0<I1/Itot){
		r1=ran();
		K=-T*log(r1);
	}
	else if(r0<(I1+I2)/Itot){
		r1=ran();
		r2=ran();
		K=-T*log(r1*r2);
	}
	else{
		r1=ran();
		r2=ran();
		r3=ran();
		K=-T*log(r1*r2*r3);
	}
	E=K+mass;
	pmag=sqrt(E*E-mass*mass);
	r0=ran();
	if(r0>pmag/E) goto GB_TRYAGAIN;
	phi=2.0*PI*ran();
	ctheta=1.0-2.0*ran();
	stheta=sqrt(1.0-ctheta*ctheta);
	p[3]=pmag*ctheta;
	p[1]=pmag*stheta*cos(phi);
	p[2]=pmag*stheta*sin(phi);
	p[0]=E;
}

void Crandy::generate_boltzmann(double mass,double T,FourVector &p){
	//const double PI=4.0*atan(1.0);
	double r1,r2,r3,a,b,c;
	double pmag,ctheta,stheta,phi;
	if(T/mass>0.6){
		GB_TRYAGAIN:
		r1=ran();
		r2=ran();
		r3=ran();
		a=-log(r1); b=-log(r2); c=-log(r3);
		pmag=T*(a+b+c);
		p[0]=sqrt(pmag*pmag+mass*mass);
		if(ran()>exp((pmag-p[0])/T)) goto GB_TRYAGAIN;
		ctheta=(a-b)/(a+b);
		stheta=sqrt(1.0-ctheta*ctheta);
		phi=T*T*pow(a+b,2)/(pmag*pmag);
		phi=2.0*PI*phi;
		p[3]=pmag*ctheta;
		p[1]=pmag*stheta*cos(phi);
		p[2]=pmag*stheta*sin(phi);
	}
	else generate_boltzmann_alt(mass,T,p);
}

int Crandy::poisson(){
	return ranp(mt);
}

void Crandy::set_mean(double mu){
	ranp=std::poisson_distribution<int>(mu);

	//ranp.reset(mu);
	//using param_t = std::poisson_distribution<int>::param_type;
	//ranp.param(param_t{mu});
}

double Crandy::ran_lorentzian(){  //
	double r=ran();
	return 0.5*(tan(PI*(r-0.5)));
}

double Crandy::ran_invcosh(){ // return propto 1/cosh(x)
	double r=ran();
	return asinh(tan(PI*(r-0.5)));
}

void Crandy::increase_threshold(){
	threshold-=log(ran());
}

void Crandy::increment_netprob(double delN){
	netprob+=delN;
}

bool Crandy::test_threshold(double delprob){
	if((delprob+netprob)<threshold){
		return false;
	}
	else
		return true;
}
