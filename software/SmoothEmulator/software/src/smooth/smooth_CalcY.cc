#include "msu_smooth/smooth.h"
using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

double CSmooth::CalcY(vector<double> &A,double LAMBDA,vector<double> &theta){
	unsigned int ic,ir;
	double answer=0.0,term,rfactor;
	rfactor=GetRFactor(LAMBDA,theta);
	answer=0.0;
	for(ic=0;ic<NCoefficients;ic++){
		term=A[ic]*sqrt(double(dupfactor[ic])/double(factorial[rank[ic]]));
		for(ir=0;ir<rank[ic];ir++){
			unsigned int ipar=IPar[ic][ir];
			term*=theta[ipar]/LAMBDA;
		}
		answer+=term;
	}
	answer*=rfactor;
	return answer;
}

double CSmooth::CalcY_Remainder(vector<double> &A,double LAMBDA,vector<double> &theta,unsigned int NTrainingPts){
	unsigned int ic,ir;
	double answer=0.0,term,rfactor;
	rfactor=GetRFactor(LAMBDA,theta);
	answer=0.0;
	for(ic=NTrainingPts;ic<NCoefficients;ic++){
		term=A[ic]*sqrt(double(dupfactor[ic])/double(factorial[rank[ic]]));
		for(ir=0;ir<rank[ic];ir++){
			term*=theta[IPar[ic][ir]]/LAMBDA;
		}
		answer+=term;
	}
	answer*=rfactor;
	return answer;
}

double CSmooth::CalcY_Remainder_FromT(vector<double> &A,unsigned int NTrainingPts,vector<double> &T){
	double answer=0.0;
	unsigned int ic;
	for(ic=NTrainingPts;ic<NCoefficients;ic++){
		answer+=T[ic]*A[ic];
	}
	return answer;
}

double CSmooth::CalcY_FromT(vector<double> &A,vector<double> &T){
	double answer=0.0;
	unsigned int ic;
	for(ic=0;ic<NCoefficients;ic++){
		answer+=T[ic]*A[ic];
	}
	return answer;
}

void CSmooth::CalcYDYDTheta(vector<double> &A,double LAMBDA,vector<double> &theta,double &Y,vector<double> &dYdTheta){
	unsigned int ir,ic,ipar,n,r;
	dYdTheta.resize(theta.size(),0.0);
	double rfactor=GetRFactor(LAMBDA, theta);
	double term,prefactor;

	vector<unsigned int> nparvec;
	vector<unsigned int> iparvec;
	dYdTheta.resize(NPars);
	Y=0.0;
	for(ipar=0;ipar<NPars;ipar++)
		dYdTheta[ipar]=0.0;
	
	vector<double> remainder(NPars+1);
	remainder[0]=1.0;
	for(ic=0;ic<NCoefficients;ic++){
		r=rank[ic];
		prefactor=A[ic]*sqrt(double(dupfactor[ic])/double(factorial[r]));
		for(n=1;n<=r;n++){
			ir=r-n;
			ipar=IPar[ic][ir];
			remainder[n]=remainder[n-1]*theta[ipar]/LAMBDA;
		}
		Y+=remainder[r]*prefactor;
		
		term=prefactor;
		for(ir=0;ir<r;ir++){
			ipar=IPar[ic][ir];
			dYdTheta[ipar]+=(term/LAMBDA)*remainder[r-ir-1];
			term*=theta[ipar]/LAMBDA;
		}
	}	
	for(ipar=0;ipar<dYdTheta.size();ipar++) {
		dYdTheta[ipar]*=rfactor;
	}
	Y*=rfactor;
}

/*
void CSmooth::CalcYDYDTheta(vector<double> &A,double LAMBDA,vector<double> &theta,double &Y,vector<double> &dYdTheta){
	unsigned int ir,ic,ipar,n,ndiff,oldipar,npar;
	dYdTheta.resize(theta.size(),0.0);
	double rfactor=GetRFactor(LAMBDA, theta);
	double term,prefactor;

	vector<unsigned int> nparvec;
	vector<unsigned int> iparvec;
	nparvec.resize(0);
	iparvec.resize(0);
	dYdTheta.resize(NPars);
	Y=0.0;
	for(ipar=0;ipar<NPars;ipar++)
		dYdTheta[ipar]=0.0;
	vector<double> remainder(NPars+1);
	remainder[0]=1.0;	

	for(ic=0;ic<NCoefficients;ic++){
		prefactor= A[ic]*sqrt(double(dupfactor[ic])/double(factorial[rank[ic]]));
		iparvec.clear();
		nparvec.clear();
		ndiff=0;
		oldipar=99999;
		for(ir=0;ir<rank[ic];ir++){
			ipar=IPar[ic][ir];
			if(ipar!=oldipar){
				iparvec.push_back(ipar);
				nparvec.push_back(1);
				ndiff+=1;
				oldipar=ipar;
			}
			else{
				nparvec[ndiff-1]+=1;
			}
		}

			
		for(n=1;n<=ndiff;n++){
			ipar=iparvec[ndiff-n];
			npar=nparvec[ndiff-n];
			remainder[n]=remainder[n-1]*pow(theta[ipar]/LAMBDA,npar);
		}
		Y+=remainder[ndiff]*prefactor;
		
		term=prefactor;
		for(n=0;n<ndiff;n++){
			ipar=iparvec[n];
			npar=nparvec[n];
			dYdTheta[ipar]+=term*(double(npar)/LAMBDA)*pow(theta[ipar]/LAMBDA,npar-1)*remainder[ndiff-n-1];
			term*=pow(theta[ipar]/LAMBDA,npar);
		}
		
	}	
	for(ipar=0;ipar<dYdTheta.size();ipar++) {
		dYdTheta[ipar]*=rfactor;
	}
	Y*=rfactor;
}
*/
