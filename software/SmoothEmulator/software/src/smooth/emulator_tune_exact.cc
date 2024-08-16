#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUUtils;

void CSmoothEmulator::GetExactQuantities(){
	unsigned int NCoefficients=smooth->NCoefficients;
	unsigned int i,a,b,aprime,bprime;
	
	B.resize(NTrainingPts,NTrainingPts);	
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			B(a,b)=0.0;
			for(i=0;i<NCoefficients;i++){
				B(a,b)+=T[a][i]*T[b][i];
			}
		}
	}
	Binv=B.inverse();
	
	chi.resize(NTrainingPts);
	for(a=0;a<NTrainingPts;a++){
		chi(a)=0.0;
		for(b=0;b<NTrainingPts;b++)
			chi(a)+=Binv[a][b]*smoothmaster->traininginfo->YTrain[iY][b];
	}
}

void CSmoothEmulator::TuneExact(){
	unsigned int NCoefficients=smooth->NCoefficients;
	unsigned int itrain,ic,a,b;
	ABest.resize(NCoefficients);
	GetExactQuantities();
	GetExactSigmaA();	
	for(ic=0;ic<NCoefficients;ic++){
		ABest[ic]=0.0;
		for(a=0;a<NTrainingPts;a++){
			for(b=0;b<NTrainingPts;b++){
				ABest[ic]+=T[a][ic]*Binv(a,b)*smoothmaster->traininginfo->YTrain[iY][b];
			}
		}
	}
}

void CSmoothEmulator::GetExactSigmaA(){
	unsigned int a,b;
	double sigmaA2=0.0;
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			sigmaA2+=smoothmaster->traininginfo->YTrain[iY][a]*Binv(a,b)*smoothmaster->traininginfo->YTrain[iY][b];
		}
	}
	sigmaA2=sigmaA2/double(NTrainingPts);
	SigmaA=sqrt(sigmaA2);
}

void CSmoothEmulator::CalcExactLogP(){
	double P=pow(SigmaA,-NTrainingPts)/sqrt(B.determinant());
	logP=log(P);
}


