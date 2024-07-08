#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUUtils;

void CSmoothEmulator::CalcY(CModelParameters *modpars,double &Y,double &SigmaY_emulator){
	CalcY(modpars->Theta,Y,SigmaY_emulator);
}

void CSmoothEmulator::CalcY(vector<double> &Theta,double &Y,double &SigmaY_emulator){
	if(pca_ignore){
		Y=0.0;
		SigmaY_emulator=0.0;
	}
	else{
		if(TuneChooseExact){
			//Y=smooth->CalcY(ABest,LAMBDA,Theta);
			//GetExactUncertainty(Theta,SigmaY_emulator);
			CalcYAndExactUncertainty(Theta,Y,SigmaY_emulator);
		}
		else{
			double y;
			Y=SigmaY_emulator=0.0;
			for(unsigned int isample=0;isample<NASample;isample++){
				y=smooth->CalcY(ASample[isample],LAMBDA,Theta);
				Y+=y;
				SigmaY_emulator+=y*y;
			}
			SigmaY_emulator=SigmaY_emulator/double(NASample);
			Y=Y/double(NASample);
			SigmaY_emulator=sqrt(fabs(SigmaY_emulator-Y*Y));
		}
	}
}

void CSmoothEmulator::CalcYOnly(CModelParameters *modpars,double &Y){
	CalcYOnly(modpars->Theta,Y);
}

void CSmoothEmulator::CalcYOnly(vector<double> &Theta,double &Y){
	if(pca_ignore){
		Y=0.0;
	}
	else{
		if(TuneChooseExact){
			Y=smooth->CalcY(ABest,LAMBDA,Theta);
		}
		else{
			double y;
			Y=0.0;
			for(unsigned int isample=0;isample<NASample;isample++){
				y=smooth->CalcY(ASample[isample],LAMBDA,Theta);
				Y+=y;
			}
			Y=Y/double(NASample);
		}
	}
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
		if(TuneChooseExact){
			smooth->CalcYDYDTheta(ABest,LAMBDA,Theta,Y,dYdTheta);
			GetExactUncertainty(Theta,SigmaY_emulator);
		}
		else{
			double y;
			unsigned int ipar;
			char dummy[200];
			vector<double> dydtheta;
			dYdTheta.resize(NPars);
			dydtheta.resize(NPars);
			for(ipar=0;ipar<NPars;ipar++){
				dYdTheta[ipar]=0.0;
			}
			Y=SigmaY_emulator=0.0;
			for(unsigned int isample=0;isample<NASample;isample++){
				smooth->CalcYDYDTheta(ASample[isample],LAMBDA,Theta,y,dydtheta);
				Y+=y/double(NASample);
				SigmaY_emulator+=y*y/double(NASample);
				for(ipar=0;ipar<NPars;ipar++){
					dYdTheta[ipar]+=dydtheta[ipar]/double(NASample);
					if(dYdTheta[ipar]!=dYdTheta[ipar]){
						snprintf(dummy,200,"ipar=%u, y=%g, dydtheta=%g, LAMBDA=%g, Theta=%g\n",ipar,y,dydtheta[ipar],LAMBDA,Theta[ipar]);
						CLog::Info(dummy);
						CLog::Fatal("disaster in CalcYDYDTheta\n");
					}
				}
			}
			SigmaY_emulator=sqrt(fabs(SigmaY_emulator-Y*Y));
		}
	}
}


void CSmoothEmulator::GetExactUncertainty(vector<double> &Theta_s,double &uncertainty){
	double unc2; // squared uncertainty
	unsigned int i,a,b,NCoefficients=smooth->NCoefficients;
	Eigen::VectorXd T,S;
	T.resize(NCoefficients);
	S.resize(NTrainingPts);

	for(i=0;i<NCoefficients;i++){
		T(i)=smooth->GetT(i,LAMBDA,Theta_s);
	}
	
	for(a=0;a<NTrainingPts;a++){
		S(a)=0.0;
		for(i=NTrainingPts;i<NCoefficients;i++){
			S(a)+=beta(a,i)*T(i);
		}
	}

	unc2=0.0;
	// First term
	for(i=NTrainingPts;i<NCoefficients;i++){
		unc2+=T(i)*T(i);
	}
	// Second  & third terms
	for(a=0;a<NTrainingPts;a++){
		unc2-=2*T(a)*S(a);
	}
	// Fourth term
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			unc2+=T(a)*B[a][b]*T(b);
		}
	}
	// Fifth term
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			unc2+=S(a)*Psi(a,b)*S(b);
		}
	}
	// Sixth & seventh terms
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			unc2-=2*T(a)*H6[a][b]*S(b);
		}
	}
	/// 8th term
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			unc2+=T(a)*H8[a][b]*T(b);
		}
	}
	if(unc2<-1.0E-8){
		CLog::Info("Inside CSmoothEmulator::GetExactUncertainty, sigma^2 is less than zero = "+to_string(unc2)+"\n");
	}
	uncertainty=SigmaA*sqrt(fabs(unc2));

}

void CSmoothEmulator::CalcYAndExactUncertainty(vector<double> &Theta_s,double &Y,double &uncertainty){
	double unc2; // squared uncertainty
	unsigned int i,a,b,NCoefficients=smooth->NCoefficients;
	//Eigen::VectorXd T,S;
	vector<double> T,S;
	T.resize(NCoefficients);
	S.resize(NTrainingPts);

	for(i=0;i<NCoefficients;i++){
		T[i]=smooth->GetT(i,LAMBDA,Theta_s);
	}
	Y=smooth->CalcY_FromT(ABest,T);
	
	for(a=0;a<NTrainingPts;a++){
		S[a]=0.0;
		for(i=NTrainingPts;i<NCoefficients;i++){
			S[a]+=beta(a,i)*T[i];
		}
	}

	unc2=0.0;
	// First term
	for(i=NTrainingPts;i<NCoefficients;i++){
		unc2+=T[i]*T[i];
	}
	// Second  & third terms
	for(a=0;a<NTrainingPts;a++){
		unc2-=2*T[a]*S[a];
	}
	// Fourth term
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			unc2+=T[a]*B[a][b]*T[b];
		}
	}
	// Fifth term
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			unc2+=S[a]*Psi(a,b)*S[b];
		}
	}
	// Sixth & seventh terms
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			unc2-=2*T[a]*H6[a][b]*S[b];
		}
	}
	/// 8th term
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			unc2+=T[a]*H8[a][b]*T[b];
		}
	}
	if(unc2<-1.0E-8){
		CLog::Info("Inside CSmoothEmulator::GetExactUncertainty, sigma^2 is less than zero = "+to_string(unc2)+"\n");
	}
	uncertainty=SigmaA*sqrt(fabs(unc2));

}
