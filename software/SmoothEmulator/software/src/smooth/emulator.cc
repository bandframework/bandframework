#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUUtils;

unsigned int CSmoothEmulator::NPars=0;
CSmoothMaster *CSmoothEmulator::smoothmaster=NULL;
CparameterMap *CSmoothEmulator::parmap=NULL;
unsigned int CSmoothEmulator::NTrainingPts=0;
unsigned int CSmoothEmulator::NTestingPts=0;

CSmoothEmulator::CSmoothEmulator(string observable_name_set){
	observable_name=observable_name_set;
	NTrainingPts=smoothmaster->traininginfo->NTrainingPts;
   FIXLAMBDA=parmap->getB("SmoothEmulator_FixLambda",false);
	LAMBDA=parmap->getD("SmoothEmulator_LAMBDA",2.5);
	INCLUDE_LAMBDA_UNCERTAINTY=parmap->getB("SmoothEmulator_INCLUDE_LAMBDA_UNCERTAINTY",true);
	iY=smoothmaster->observableinfo->GetIPosition(observable_name);
	ALPHA=smoothmaster->observableinfo->ALPHA[iY];
	ThetaTrain.clear();
}

void CSmoothEmulator::CalcB(){
	unsigned int a,b;
	ThetaTrain.clear();
	ThetaTrain.resize(NTrainingPts);

	for(unsigned int itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		ThetaTrain[itrain]=smoothmaster->traininginfo->modelpars[itrain]->Theta;
	}
	B.resize(NTrainingPts,NTrainingPts);
	Binv.resize(NTrainingPts,NTrainingPts);
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			B(a,b)=GetCorrelation(ThetaTrain[a],ThetaTrain[b]);
		}
	}
	for(a=0;a<NTrainingPts;a++)
		B(a,a)+=ALPHA*ALPHA;
	Binv=B.inverse();
	chi.resize(NTrainingPts);
	for(a=0;a<NTrainingPts;a++){
		chi[a]=0.0;
		for(b=0;b<NTrainingPts;b++){
			chi[a]+=Binv(a,b)*smoothmaster->traininginfo->YTrain[iY][b];
		}
	}
}

double CSmoothEmulator::GetCorrelation(vector<double> &Theta1,vector<double> &Theta2){
	unsigned int ipar;
	double delTheta,delThetaSquared=0.0;
	for(ipar=0;ipar<NPars;ipar++){
		delTheta=Theta1[ipar]-Theta2[ipar];
		delThetaSquared+=delTheta*delTheta;
	}
	return exp(-0.5*delThetaSquared/(LAMBDA*LAMBDA));
}

void CSmoothEmulator::Tune(){
   if(FIXLAMBDA)
      Tune(LAMBDA);
   else{
      CalcSigmaALambda();
   }
}

void CSmoothEmulator::Tune(double LambdaSet){
	LAMBDA=LambdaSet;
	CalcB();
	CalcWBprimeChi();
	CalcSigmaA();
	CalcLogP();
}

