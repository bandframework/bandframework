#include "msu_commonutils/sf.h"

using namespace std;

double Bessel::J0(double x){
	return gsl_sf_bessel_J0(x);
}

double Bessel::J1(double x){
	return gsl_sf_bessel_J1(x);
}

double Bessel::Jn(int n,double x){
	return gsl_sf_bessel_Jn(n,x);
}

double Bessel::Y0(double x){
	return gsl_sf_bessel_Y0(x);
}

double Bessel::Y1(double x){
	return gsl_sf_bessel_Y1(x);
}

double Bessel::Yn(int n,double x){
	return gsl_sf_bessel_Yn(n,x);
}

double Bessel::K0(double x){
	return gsl_sf_bessel_K0(x);
}

double Bessel::K1(double x){
	return gsl_sf_bessel_K1(x);
}

double Bessel::Kn(int n,double x){
	return gsl_sf_bessel_Kn(n,x);
}

double Bessel::I0(double x){
	return gsl_sf_bessel_K0(x);
}

double Bessel::I1(double x){
	return gsl_sf_bessel_K1(x);
}

double Bessel::In(int n,double x){
	return gsl_sf_bessel_Kn(n,x);
}

double Bessel::j0(double x){
	return gsl_sf_bessel_j0(x);
}

double Bessel::j1(double x){
	return gsl_sf_bessel_j1(x);
}

double Bessel::jn(int n,double x){
	double answer,oldanswer,tempanswer;
	int L;
	if(n==0) answer=j0(x);
	else{
		answer=j1(x);
		if(n>1){
			oldanswer=j0(x);
			for(L=2;L<=n;L++){
				tempanswer=answer;
				answer=-oldanswer+(2.0*L-1)*answer/x;
				oldanswer=tempanswer;
			}
		}
	}
	return answer;
}

double Bessel::y0(double x){
	return gsl_sf_bessel_y0(x);
}

double Bessel::y1(double x){
	return gsl_sf_bessel_y1(x);
}

double Bessel::yn(int n,double x){
	double answer,oldanswer,tempanswer;
	int L;
	// Check this out! I'm not sure yn uses same recursion relation as jn
	if(n==0)
		answer=y0(x);
	else{
		answer=y1(x);
		if(n>1){
			oldanswer=y0(x);
			for(L=2;L<=n;L++){
				tempanswer=answer;
				answer=-oldanswer+(2.0*L-1)*answer/x;
				oldanswer=tempanswer;
			}
		}
	}
	return answer;
}

complex<double> Bessel::h0(double x){
	complex<double> ci(0.0,1.0);
	return (gsl_sf_bessel_j0(x)+ci*gsl_sf_bessel_y0(x));
}

complex<double> Bessel::h1(double x){
	complex<double> ci(0.0,1.0);
	return (gsl_sf_bessel_j1(x)+ci*gsl_sf_bessel_y1(x));
}

complex<double> Bessel::hn(int n,double x){
	complex<double> ci(0.0,1.0);
	return (jn(n,x)+ci*yn(n,x));
}

complex<double> Bessel::hstar0(double x){
	complex<double> ci(0.0,1.0);
	return (gsl_sf_bessel_j0(x)-ci*gsl_sf_bessel_y0(x));
}

complex<double> Bessel::hstar1(double x){
	complex<double> ci(0.0,1.0);
	return (gsl_sf_bessel_j1(x)-ci*gsl_sf_bessel_y1(x));
}

complex<double> Bessel::hstarn(int n,double x){
	complex<double> ci(0.0,1.0);
	return (jn(n,x)-ci*yn(n,x));
}

void Bessel::CalcJN_real(int ell,double x,double &jl,double &nl,double &jlprime,double &nlprime){
	int l=0;
	double j[4],n[4],jprime[4],nprime[4];
	double s=sin(x);
	double c=cos(x);
	if(ell>3){
		printf("increase array size in CalcJN_real\n");
		exit(1);
	}
	j[0]=s;
	n[0]=c;
	jprime[0]=c;
	nprime[0]=-s;
	while(l<ell){
		l+=1;
		j[l]=(double(l)/x)*j[l-1]-jprime[l-1];
		n[l]=(double(l)/x)*n[l-1]-nprime[l-1];
		jprime[l]=j[l-1]-(double(l)/x)*j[l];
		nprime[l]=n[l-1]-(double(l)/x)*n[l];
	};
	jl=j[ell];
	nl=n[ell];
	jlprime=jprime[ell];
	nlprime=nprime[ell];
}

void Bessel::CalcJN_imag(int ell,double x,double &jl,double &nl,double &jlprime,double &nlprime){
	int l=0;
	double j[4],n[4],jprime[4],nprime[4];
	double s=sinh(x);
	double c=cosh(x);
	if(ell>3){
		printf("increase array size in CalcJN_imag\n");
		exit(1);
	}
	j[0]=s;
	n[0]=c;
	jprime[0]=c;
	nprime[0]=s;
	while(l<ell){
		l+=1;
		j[l]=(double(l)/x)*j[l-1]-jprime[l-1];
		n[l]=(double(l)/x)*n[l-1]-nprime[l-1];
		jprime[l]=-j[l-1]-(double(l)/x)*j[l];
		nprime[l]=-n[l-1]-(double(l)/x)*n[l];
	};
	jl=j[ell];
	nl=n[ell];
	jlprime=jprime[ell];
	nlprime=nprime[ell];
}
