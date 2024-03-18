#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUUtils;

unsigned int CSmoothEmulator::NPars=0;
CSmooth *CSmoothEmulator::smooth=NULL;
CSmoothMaster *CSmoothEmulator::smoothmaster=NULL;
CparameterMap *CSmoothEmulator::parmap=NULL;
Crandy *CSmoothEmulator::randy=NULL;
unsigned int CSmoothEmulator::NTrainingPts=0;

CSmoothEmulator::CSmoothEmulator(string observable_name_set,bool pca_ignore_set){
	observable_name=observable_name_set;
	NTrainingPts=smoothmaster->traininginfo->NTrainingPts;

	LAMBDA=parmap->getD("SmoothEmulator_LAMBDA",2.0);
	NMC=parmap->getI("SmoothEmulator_MCMC_NMC",10000);
	NASample=parmap->getI("SmoothEmulator_MCMC_NASample",8);
	MCStepSize=parmap->getD("SmoothEmulator_MCMC_StepSize",0.01);
	MCSigmaAStepSize=parmap->getD("SmoothEmulator_MCMC_SigmaAStepSize",0.01);
	TuneChooseMCMC=parmap->getB("SmoothEmulator_TuneChooseMCMC",false);
	TuneChooseMCMCPerfect=parmap->getB("SmoothEmulator_TuneChooseMCMCPerfect",false);
	TuneChooseExact=parmap->getB("SmoothEmulator_TuneChooseExact",true);
	UseSigmaY=parmap->getB("SmoothEmulator_MCMC_UseSigmaY",false);
	ConstrainA0=parmap->getB("SmoothEmulator_ConstrainA0",false);
	CutOffA=parmap->getB("SmoothEmulator_MCMC_CutoffA",false);
	iY=smoothmaster->observableinfo->GetIPosition(observable_name);
	SigmaA0=smoothmaster->observableinfo->SigmaA0[iY];
	SigmaA=SigmaA0;
	SigmaAMin=parmap->getD("SmoothEmulator_SigmaAMin",0.1*SigmaA0);
	pca_ignore=pca_ignore_set;
	ThetaTrain.clear();
	Init();
}

void CSmoothEmulator::SetThetaTrain(){
	ThetaTrain.clear();
	ThetaTrain.resize(NTrainingPts);
	for(unsigned int itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(unsigned int ipar=0;ipar<NPars;ipar++){
			ThetaTrain[itrain][ipar]=smoothmaster->traininginfo->modelpars[itrain]->Theta[ipar];
		}
	}
}

void CSmoothEmulator::Init(){
	SigmaA=SigmaA0;
	if(!pca_ignore){
		NSigmaA=0;
		SigmaAbar=0.0;
		FirstTune=true;
		MCStepSize=MCStepSize/double(NPars*NPars);
		MCSigmaAStepSize=MCSigmaAStepSize/double(NPars*NPars);

		if(!TuneChooseExact){
			ASample.resize(NASample);
			for(unsigned int isample=0;isample<NASample;isample++){
				ASample[isample].resize(smooth->NCoefficients);
				//SetA_RanGauss(SigmaA,ASample[isample]);
				SetA_Zero(ASample[isample]);
			}
		}

		A.resize(smooth->NCoefficients);
		SetA_Zero(A);
		ATrial.resize(smooth->NCoefficients);
		SetA_Zero(ATrial);
		Ttilde.resize(NTrainingPts,NTrainingPts);
		TtildeInv.resize(NTrainingPts,NTrainingPts);
		T.resize(NTrainingPts);
		for(unsigned int it=0;it<NTrainingPts;it++){
			T[it].resize(smooth->NCoefficients);
			for(unsigned int ic=0;ic<smooth->NCoefficients;ic++){
				T[it][ic]=0.0;
			}
		}
	}
	ThetaTrain.clear();
}

void CSmoothEmulator::Tune(){
	//FirstTune=true;
	if(!pca_ignore){
		CalcTForTraining();
		if(TuneChooseMCMC==true){
			if(UseSigmaY){
				if(FirstTune){
					TuneMCMC();
					FirstTune=false;
				}
				TuneMCMC_withSigma();
			}
			else{
				if(FirstTune){
					TuneMCMC();
					FirstTune=false;
				}
				TuneMCMC();
			}
		}
		else if(TuneChooseMCMCPerfect==true){
			TunePerfectMCMC();
		}
		else if(TuneChooseExact){
			TuneExact();
			GetExactSigmaA();
			CalcExactLogP();
		}
		else{
			CLog::Fatal("In CSmoothEmulator::Tune(), no tuning method specified\n");
		}
	}
}

void CSmoothEmulator::CalcTForTraining(){
	if(!pca_ignore){
		unsigned int itrain,ic;
		for(itrain=0;itrain<NTrainingPts;itrain++){
			for(ic=0;ic<smooth->NCoefficients;ic++){
				T[itrain][ic]=smooth->GetT(ic,LAMBDA,ThetaTrain[itrain]);
			}
			for(ic=0;ic<NTrainingPts;ic++){
				Ttilde(itrain,ic)=T[itrain][ic];
			}
		}
		TtildeInv=Ttilde.inverse();	
	
		/* cout << "---------\n";
		cout << Ttilde << endl;
		double detTtilde=Ttilde.determinant();
		CLog::Info("ln(determinant(Ttilde))="+to_string(log(detTtilde))+"\n");
		cout << TtildeInv << endl; */
	
		for(itrain=0;itrain<NTrainingPts;itrain++){
			for(ic=0;ic<NTrainingPts;ic++){
				if(TtildeInv(itrain,ic)!=TtildeInv(itrain,ic)){
					CLog::Fatal("TtildeInv != TtildeInv\n");
				}
			}
		}
	}
}

// This adjust first NTrainingPts coefficients to reproduce Y training values

void CSmoothEmulator::CalcAFromTraining(vector<double> &AA){
	if(!pca_ignore){
		vector<double> YTrain=smoothmaster->traininginfo->YTrain[iY];
		unsigned int itrain;
		//vector<double> Ashort;
		//Ashort.resize(smooth->NCoefficients);
		Eigen::VectorXd Ashort,YTarget;
		AA.resize(NTrainingPts);
		YTarget.resize(NTrainingPts);
		if(ThetaTrain.size()!=NTrainingPts){
			CLog::Info("CSmoothEmulator:: array size mismatch!!\n");
			CLog::Fatal("ThetaTrain.size="+to_string(ThetaTrain.size())+", NTrainingPts="+to_string(NTrainingPts)+"\n");
		}
		YTarget.resize(NTrainingPts);

		for(itrain=0;itrain<NTrainingPts;itrain++){
			//YTarget(itrain)=YTrain[itrain]-smooth->CalcY_Remainder(AA,LAMBDA,ThetaTrain[itrain],NTrainingPts);
			YTarget(itrain)=YTrain[itrain]-smooth->CalcY_Remainder_FromT(AA,NTrainingPts,T[itrain]);
		}
		Ashort=TtildeInv*YTarget;
		for(itrain=0;itrain<NTrainingPts;itrain++)
			AA[itrain]=Ashort(itrain);
	}
}

void CSmoothEmulator::OldCalcAFromTraining(vector<double> &AA){
	vector<double> YTrain=smoothmaster->traininginfo->YTrain[iY];
	unsigned int itrain,ic;
	//vector<double> Ashort;
	//Ashort.resize(smooth->NCoefficients);
	Eigen::VectorXd Ashort,YTarget;
	AA.resize(NTrainingPts);
	YTarget.resize(NTrainingPts);
	if(ThetaTrain.size()!=NTrainingPts){
		CLog::Info("CSmoothEmulator:: array size mismatch!!\n");
		CLog::Fatal("ThetaTrain.size="+to_string(ThetaTrain.size())+", NTrainingPts="+to_string(NTrainingPts)+"\n");
	}
	YTarget.resize(NTrainingPts);

	for(itrain=0;itrain<NTrainingPts;itrain++){
		YTarget(itrain)=YTrain[itrain]-smooth->CalcY_Remainder(AA,LAMBDA,ThetaTrain[itrain],NTrainingPts);
		for(ic=0;ic<NTrainingPts;ic++){
			Ttilde(itrain,ic)=smooth->GetT(ic,LAMBDA,ThetaTrain[itrain]);
		}
	}
	//Ashort=Ttilde.colPivHouseholderQr().solve(YTarget);
	Ashort=Ttilde.partialPivLu().solve(YTarget);
	for(itrain=0;itrain<NTrainingPts;itrain++)
		AA[itrain]=Ashort(itrain);
}

double CSmoothEmulator::GetLog_AProb(vector<double> &AA,double ASigmaA){
	double answer=0.0;
	if(!pca_ignore){
		unsigned int ic0=1;
		if(ConstrainA0)
			ic0=0;
		//answer=0.000*smooth->NCoefficients*log(ASigmaA/SigmaA0);
		// Don't prefer small A[0], so start at ic=1
		for(unsigned int ic=ic0;ic<smooth->NCoefficients;ic++){
			answer-=0.5*AA[ic]*AA[ic]/(ASigmaA*ASigmaA);
		}
		if(!UseSigmaY){
			if(ConstrainA0)
				answer-=NTrainingPts*log(ASigmaA);
			else
				answer-=(NTrainingPts-1)*log(ASigmaA);
		}
		// next line keeps A from drifting out to infinity
		if(CutOffA)
			answer-=log(1.0+0.25*(ASigmaA*ASigmaA)/(SigmaA0*SigmaA0));
	}
	return answer;
}

void CSmoothEmulator::SetA_Zero(vector<double> &A){
	if(!pca_ignore){
		for(unsigned int ic=0;ic<A.size();ic++)
			A[ic]=0.0;
	}
}

void CSmoothEmulator::SetA_RanGauss(double SigmaA,vector<double> &A){
	if(!pca_ignore){
		for(unsigned int ic=0;ic<A.size();ic++)
			A[ic]=SigmaA*randy->ran_gauss();
	}
}

void CSmoothEmulator::SetA_Constant(double SigmaA,vector<double> &A){
	if(!pca_ignore){
		for(unsigned int ic=0;ic<A.size();ic++)
			A[ic]=SigmaA;
	}
}

void CSmoothEmulator::SetA_RanSech(double SigmaA,vector<double> &A){
	if(!pca_ignore){
		const double PI=4.0*atan(1.0);
		double r=1.0-2.0*randy->ran();
		for(unsigned int ic=0;ic<A.size();ic++)
			A[ic]=SigmaA*2.0*atanh(tan(0.5*PI*(r-0.5)));
	}
}
