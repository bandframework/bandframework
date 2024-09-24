#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUUtils;

double  CSmoothEmulator::GetYOnly(CModelParameters *modpars){
	return GetYOnly(modpars->Theta);
}

double CSmoothEmulator::GetYOnly(vector<double> &Theta){
	double Y;
	if(pca_ignore){
		Y=0.0;
	}
	else{
		Y=smooth->CalcY(ABest,LAMBDA,Theta);
	}
	return Y;
}

void CSmoothEmulator::CalcYDYDTheta(CModelParameters *modpars,double &Y,vector<double> &dYdTheta,double &SigmaY_emulator){
	if(pca_ignore){
		Y=0.0;
		SigmaY_emulator=0.0;
		for(unsigned int ipar=0;ipar<NPars;ipar++){
			dYdTheta[ipar]=0.0;
		}
	}
	else{
		CalcYDYDTheta(modpars->Theta,Y,dYdTheta,SigmaY_emulator);
	}
}

void CSmoothEmulator::CalcYDYDTheta(vector<double> &Theta,double &Y,vector<double> &dYdTheta,double &SigmaY_emulator){
	if(pca_ignore){
		Y=0.0;
		dYdTheta.resize(NPars);
		for(unsigned int ipar=0;ipar<NPars;ipar++){
			dYdTheta[ipar]=0.0;
		}
	}
	else{
		smooth->CalcYDYDTheta(ABest,LAMBDA,Theta,Y,dYdTheta);
		SigmaY_emulator=GetUncertainty(Theta);
	}
}

double CSmoothEmulator::GetUncertainty(CModelParameters *modpars){
	return GetUncertainty(modpars->Theta);
}

double CSmoothEmulator::GetUncertainty(vector<double> &Theta_s){
	double unc2; // squared uncertainty
	unsigned int i,a,b,NCoefficients=smooth->NCoefficients;
	vector<double> T,S;
	T.resize(NCoefficients);
	S.resize(NTrainingPts);

	for(i=0;i<NCoefficients;i++){
		T[i]=smooth->GetT(i,LAMBDA,Theta_s);
	}
	
	for(a=0;a<NTrainingPts;a++){
		S[a]=0.0;
		for(i=0;i<NCoefficients;i++){
			S[a]+=TTrain[a][i]*T[i];
		}
	}

	unc2=0.0;
	for(i=0;i<NCoefficients;i++){
		unc2+=T[i]*T[i];
	}
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			unc2-=S[a]*Binv(a,b)*S[b];
		}
	}
	return SigmaA*sqrt(fabs(unc2));

}
	
void CSmoothEmulator::CalcYAndUncertainty(vector<double> &Theta_s,double &Y,double &uncertainty){
	Y=GetYOnly(Theta_s);
	uncertainty=GetUncertainty(Theta_s);
}
