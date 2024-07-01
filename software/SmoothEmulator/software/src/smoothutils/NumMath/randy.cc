#include "msu_smoothutils/randy.h"
using namespace NMSUUtils;

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
	const double PI=4.0*atan(1.0);
	double r=ran();
	return 0.5*(tan(PI*(r-0.5)));
}

double Crandy::ran_invcosh(){ // return propto 1/cosh(x)
	const double PI=4.0*atan(1.0);
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
