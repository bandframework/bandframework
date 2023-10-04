#include "msu_commonutils/sf.h"
#include "msu_commonutils/misc.h"

using namespace std;

double SpherHarmonics::legendre(int ell,double ctheta){
	if(ctheta>1.0) ctheta=1.0;
	if(ctheta<-1.0) ctheta=-1.0;
  return gsl_sf_legendre_Pl(ell,ctheta);
}

complex<double> SpherHarmonics::Ylm(int ell, int m, double theta, double phi){
  double ctheta;
  complex<double> answer;
  ctheta=cos(theta);
  answer=gsl_sf_legendre_sphPlm(ell,abs(m),ctheta)*Misc::ceiphi(static_cast<double>(m)*phi);
  if(m<0) answer *= pow(-1.0,(abs(m))); //pow(-1.0,abs(m));
  return answer;
}

complex<double> SpherHarmonics::Ylm(int ell, int m, double x, double y, double z){
    complex<double> answer;
    double ctheta,phi;
    double r = sqrt(x*x+y*y+z*z);
    if ( r < 1e-10 || fabs(z) < 1e-10 ) ctheta = 0.0;
    else ctheta=z/r;
    phi=atan2(y,x);
    answer=gsl_sf_legendre_sphPlm(ell,abs(m),ctheta)*Misc::ceiphi(static_cast<double>(m)*phi);
    if(m<0) answer *= pow(-1.0,(abs(m)));//*pow(-1.0,abs(m));
    return answer;	
}

double SpherHarmonics::ReYlm(int ell, int m, double theta, double phi){
	return real(SpherHarmonics::Ylm(ell,m,theta,phi));
}

double SpherHarmonics::ImYlm(int ell, int m, double theta, double phi){
	return imag(SpherHarmonics::Ylm(ell,m,theta,phi));
}

double SpherHarmonics::ReYlm(int ell, int m, double x,double y,double z){
	return real(SpherHarmonics::Ylm(ell,m,x,y,z));
}

double SpherHarmonics::ImYlm(int ell, int m, double x,double y,double z){
	return imag(SpherHarmonics::Ylm(ell,m,x,y,z));
}
