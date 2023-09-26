#include "msu_smooth/emulator.h"
using namespace std;

int CSmoothEmulator::NPars=0;
CSmooth *CSmoothEmulator::smooth=NULL;
CSmoothMaster *CSmoothEmulator::smoothmaster=NULL;
CparameterMap *CSmoothEmulator::parmap=NULL;
Crandy *CSmoothEmulator::randy=NULL;
unsigned int CSmoothEmulator::NTrainingPts=0;

CSmoothEmulator::CSmoothEmulator(string observable_name_set){
	observable_name=observable_name_set;
	NTrainingPts=smoothmaster->traininginfo->NTrainingPts;

	LAMBDA=parmap->getD("SmoothEmulator_LAMBDA",3.0);
	NMC=parmap->getI("SmoothEmulator_NMC",10000);
	NASample=parmap->getI("SmoothEmulator_NASample",8);
	MCStepSize=parmap->getD("SmoothEmulator_MCStepSize",0.01);
	MCSigmaAStepSize=parmap->getD("SmoothEmulator_MCSigmaAStepSize",0.01);
	TuneChooseMCMC=parmap->getB("SmoothEmulator_TuneChooseMCMC",true);
	UseSigmaYReal=parmap->getB("SmoothEmulator_UseSigmaYRreal",false);
	ConstrainA0=parmap->getB("SmoothEmulator_ConstrainA0",false);
	CutOffA=parmap->getB("SmoothEmulator_CutoffA",false);
	iY=smoothmaster->observableinfo->GetIPosition(observable_name);
	SigmaA0=smoothmaster->observableinfo->SigmaA0[iY];
	SigmaA=SigmaA0;
	SigmaAMin=parmap->getD("SmoothEmulator_SigmaAMin",0.1*SigmaA0);

	Init();

}

void CSmoothEmulator::SetThetaTrain(){
	ThetaTrain.resize(NTrainingPts);
	for(int itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(int ipar=0;ipar<NPars;ipar++){
			ThetaTrain[itrain][ipar]=smoothmaster->traininginfo->modelpars[itrain]->Theta[ipar];
		}
	}
}

void CSmoothEmulator::Init(){
	SigmaA=SigmaA0;
	NSigmaA=0;
	SigmaAbar=0.0;
	FirstTune=true;
	MCStepSize=MCStepSize/double(NPars*NPars);
	MCSigmaAStepSize=MCSigmaAStepSize/double(NPars*NPars);

	ASample.resize(NASample);
	for(unsigned int isample=0;isample<NASample;isample++){
		ASample[isample].resize(smooth->NCoefficients);
		//SetA_RanGauss(SigmaA,ASample[isample]);
		SetA_Zero(ASample[isample]);
	}

	A.resize(smooth->NCoefficients);
	SetA_Zero(A);
	ATrial.resize(smooth->NCoefficients);
	SetA_Zero(ATrial);
	M.resize(NTrainingPts,NTrainingPts);
}

void CSmoothEmulator::Tune(){
	FirstTune=true;
	if(TuneChooseMCMC==true){
		if(UseSigmaYReal){
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
	else{
		TunePerfect();
	}
}

void CSmoothEmulator::TuneMCMC(){
	vector<double> *Aswitch,*Aptr,*ATrialptr;
	double dlp,r,SigmaAswitch;
	unsigned int success=0,ic,imc;
	double BestLogP,stepsize;
	CalcAFromTraining(A);
	for(ic=0;ic<smooth->NCoefficients;ic++){
		ATrial[ic]=A[ic];
	}
	Aptr=&A;
	ATrialptr=&ATrial;
	SigmaATrial=SigmaA0;
	double logP,logPTrial;
	logP=GetLog_AProb(*Aptr,SigmaA);
	BestLogP=-1000000.0;
	for(imc=0;imc<NMC;imc++){
		for(ic=NTrainingPts;ic<smooth->NCoefficients;ic++){
			stepsize=SigmaA*MCStepSize*pow(LAMBDA,smooth->rank[ic]);
			(*ATrialptr)[ic]=(*Aptr)[ic]+stepsize*randy->ran_gauss();
		}
		do{
			SigmaATrial=fabs(SigmaA+SigmaA0*MCSigmaAStepSize*randy->ran_gauss());
		}while(SigmaATrial<SigmaAMin);
		for(ic=NTrainingPts;ic<smooth->NCoefficients;ic++){
			(*ATrialptr)[ic]*=(SigmaATrial/SigmaA);
		}
		CalcAFromTraining(*ATrialptr);
		logPTrial=GetLog_AProb(*ATrialptr,SigmaATrial);
		dlp=logPTrial-logP;

		if(dlp>0.0){
			Aswitch=Aptr;
			Aptr=ATrialptr;
			ATrialptr=Aswitch;
			SigmaAswitch=SigmaA;
			SigmaA=SigmaATrial;
			SigmaATrial=SigmaAswitch;
			logP=logPTrial;
			success+=1;
		}
		else{
			if(dlp>-100){
				dlp=exp(dlp);
				r=randy->ran();
				if(dlp>r){
					Aswitch=Aptr;
					Aptr=ATrialptr;
					ATrialptr=Aswitch;
					SigmaAswitch=SigmaA;
					SigmaA=SigmaATrial;
					SigmaATrial=SigmaAswitch;
					logP=logPTrial;
					success+=1;
				}
			}
		}
		if(logP>BestLogP)
			BestLogP=logP;
	}

	if(Aptr!=&A){
		for(ic=0;ic<smooth->NCoefficients;ic++){
			A[ic]=(*Aptr)[ic];
		}
	}

	int Ndof=smooth->NCoefficients-NTrainingPts;
	if(!FirstTune)
		CLog::Info("success percentage="+to_string(double(success)*100.0/double(NMC))+", SigmaA="+to_string(SigmaA)+", logP/Ndof="+to_string(logP/double(Ndof))+",BestLogP/Ndof="+to_string(BestLogP/double(Ndof))+"\n");
}

void CSmoothEmulator::TuneMCMC_withSigma(){
	vector<double> *Aswitch,*Aptr,*ATrialptr;
	double dlp,r,SigmaAswitch,Y,stepsize;
	unsigned int success=0,ic,imc,itrain;
	double BestLogP;
	vector<double> YTrain=smoothmaster->traininginfo->YTrain[iY];
	vector<double> SigmaYTrain=smoothmaster->traininginfo->SigmaYTrain[iY];
	//CalcAFromTraining(A);
	for(ic=0;ic<smooth->NCoefficients;ic++){
		ATrial[ic]=A[ic];
	}
	Aptr=&A;
	ATrialptr=&ATrial;
	SigmaATrial=SigmaA;
	double logP,logPTrial;
	logP=GetLog_AProb(*Aptr,SigmaA);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		Y=smooth->CalcY(*Aptr,LAMBDA,ThetaTrain[itrain]);
		logP-=0.5*(YTrain[itrain]-Y)*(YTrain[itrain]-Y)/(SigmaYTrain[itrain]*SigmaYTrain[itrain]);
	}

	BestLogP=-1000000000.0;
	for(imc=0;imc<NMC;imc++){
		for(ic=0;ic<smooth->NCoefficients;ic++){
			stepsize=SigmaA*MCStepSize*pow(LAMBDA,smooth->rank[ic]);
			(*ATrialptr)[ic]=(*Aptr)[ic]+stepsize*randy->ran_gauss();
		}
		SigmaATrial=fabs(SigmaA+SigmaA0*MCSigmaAStepSize*randy->ran_gauss());
		for(ic=0;ic<smooth->NCoefficients;ic++){
			(*ATrialptr)[ic]*=(SigmaATrial/SigmaA);
		}
		logPTrial=GetLog_AProb(*ATrialptr,SigmaATrial);


		for(itrain=0;itrain<NTrainingPts;itrain++){
			Y=smooth->CalcY(*ATrialptr,LAMBDA,ThetaTrain[itrain]);
			logPTrial-=0.5*(YTrain[itrain]-Y)*(YTrain[itrain]-Y)/(SigmaYTrain[itrain]*SigmaYTrain[itrain]);
		}

		dlp=logPTrial-logP;
		if(dlp>0.0){
			Aswitch=Aptr;
			Aptr=ATrialptr;
			ATrialptr=Aswitch;
			SigmaAswitch=SigmaA;
			SigmaA=SigmaATrial;
			SigmaATrial=SigmaAswitch;
			logP=logPTrial;
			success+=1;
			//			cout << "dlp is:" << dlp << endl;
		}
		else{
			if(dlp>-100){
				dlp=exp(dlp);
				r=randy->ran();
				//				cout << "dlp is:" << dlp << endl;
				//				cout << "randy is: " << r <<endl;
				if(dlp>r){
					Aswitch=Aptr;
					Aptr=ATrialptr;
					ATrialptr=Aswitch;
					SigmaAswitch=SigmaA;
					SigmaA=SigmaATrial;
					SigmaATrial=SigmaAswitch;
					logP=logPTrial;
					success+=1;
				}
			}
		}
		if(logP>BestLogP)
			BestLogP=logP;
	}

	if(Aptr!=&A){
		for(ic=0;ic<smooth->NCoefficients;ic++){
			A[ic]=(*Aptr)[ic];
		}
	}

	int Ndof=smooth->NCoefficients-NTrainingPts;
	if(!FirstTune)
		CLog::Info("success percentage="+to_string(double(success)*100.0/double(NMC))+", SigmaA="+to_string(SigmaA)+", logP/Ndof="+to_string(logP/double(Ndof))+",BestLogP/Ndof="+to_string(BestLogP/double(Ndof))+"\n");
}

void CSmoothEmulator::TunePerfect(){
	unsigned int ic,ic0,ntry=0,ntrymax=100000;
	bool success=false;
	double weight,warg;//sigmafact=1.0;
	if(ConstrainA0){
		ic0=0;
	}
	else
		ic0=1;

	while(ntry<ntrymax && success==false){
		SigmaA=SigmaA0;

		for(ic=NTrainingPts;ic<smooth->NCoefficients;ic++){
			//SigmaA=0.5*SigmaA0*tan((PI/2.0)*(1.0-2.0*randy->ran()));
			ATrial[ic]=SigmaA*randy->ran_gauss();
		}
		warg=0.0;
		CalcAFromTraining(ATrial);
		for(ic=ic0;ic<NTrainingPts;ic++){
			warg-=0.5*ATrial[ic]*ATrial[ic]/(SigmaA*SigmaA);
		}
		//warg-=(NTrainingPts-1)*log(SigmaA/(0.5*SigmaA0));
		if(warg>0.0){
			CLog::Fatal("Disaster, warg="+to_string(warg)+"\n");
		}
		if(warg>-100){
			weight=exp(warg);
			if(weight>randy->ran()){
				success=true;
			}
		}
		ntry+=1;
	}
	if(ntry>=ntrymax){
		CLog::Fatal("TunePerfect Failed, SigmaA0="+to_string(SigmaA0)+"\n");
	}
	else{
		for(ic=0;ic<smooth->NCoefficients;ic++){
			A[ic]=ATrial[ic];
		}
	}
}

// This adjust first NTrainingPts coefficients to reproduce Y training values
void CSmoothEmulator::CalcAFromTraining(vector<double> &AA){
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
		exit(1);
	}
	YTarget.resize(NTrainingPts);

	for(itrain=0;itrain<NTrainingPts;itrain++){
		YTarget(itrain)=YTrain[itrain]-smooth->CalcY_Remainder(AA,LAMBDA,ThetaTrain[itrain],NTrainingPts);
		for(ic=0;ic<NTrainingPts;ic++){
			M(itrain,ic)=smooth->GetM(ic,LAMBDA,ThetaTrain[itrain]);
		}
	}
	//Ashort=M.colPivHouseholderQr().solve(YTarget);
	Ashort=M.partialPivLu().solve(YTarget);
	for(itrain=0;itrain<NTrainingPts;itrain++)
		AA[itrain]=Ashort(itrain);
}

double CSmoothEmulator::GetLog_AProb(vector<double> &AA,double ASigmaA){
	double answer=0.0;
	int ic0=1;
	if(ConstrainA0)
		ic0=0;
	//answer=0.000*smooth->NCoefficients*log(ASigmaA/SigmaA0);
	// Don't prefer small A[0], so start at ic=1
	for(unsigned int ic=ic0;ic<smooth->NCoefficients;ic++){
		answer-=0.5*AA[ic]*AA[ic]/(ASigmaA*ASigmaA);
	}
	if(!UseSigmaYReal){
		if(ConstrainA0)
			answer-=NTrainingPts*log(ASigmaA);
		else
			answer-=(NTrainingPts-1)*log(ASigmaA);
	}
	// next line keeps A from drifting out to infinity
	if(CutOffA)
			answer-=log(1.0+0.25*(ASigmaA*ASigmaA)/(SigmaA0*SigmaA0));
	return answer;
}

void CSmoothEmulator::SetA_Zero(vector<double> &A){
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=0.0;
}

void CSmoothEmulator::SetA_RanGauss(double SigmaA,vector<double> &A){
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=SigmaA*randy->ran_gauss();
}

void CSmoothEmulator::SetA_Constant(double SigmaA,vector<double> &A){
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=SigmaA;
}

void CSmoothEmulator::SetA_RanSech(double SigmaA,vector<double> &A){
	double r=1.0-2.0*randy->ran();
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=SigmaA*2.0*atanh(tan(0.5*PI*(r-0.5)));
}

void CSmoothEmulator::GenerateASamples(){
	unsigned int isample;
	//	cout << "NASample is:" << NASample << endl;
	//NASample = 10
	FirstTune=true;
	for(isample=0;isample<NASample;isample++){
		Tune();
		for(unsigned int ic=0;ic<smooth->NCoefficients;ic++){
			ASample[isample][ic]=A[ic];
		}
	}
	SigmaAbar+=SigmaA;
	NSigmaA+=1;
}

void CSmoothEmulator::PrintA(vector<double> &Aprint){
	for(unsigned int ic=0;ic<smooth->NCoefficients;ic++){
		CLog::Info(to_string(ic)+": "+to_string(Aprint[ic])+"\n");
	}
}

void CSmoothEmulator::CalcY(CModelParameters *modpars,double &Y,double &SigmaY_emulator){
	double y;
	Y=SigmaY_emulator=0.0;
	for(int isample=0;isample<NASample;isample++){
		y=smooth->CalcY(ASample[isample],LAMBDA,modpars->Theta);
		Y+=y;
		SigmaY_emulator+=y*y;
	}
	SigmaY_emulator=SigmaY_emulator/double(NASample);
	Y=Y/double(NASample);
	SigmaY_emulator=sqrt(fabs(SigmaY_emulator-Y*Y));
}

void CSmoothEmulator::CalcY(vector<double> Theta,double &Y,double &SigmaY_emulator){
	double y;
	Y=SigmaY_emulator=0.0;
	for(int isample=0;isample<NASample;isample++){
		y=smooth->CalcY(ASample[isample],LAMBDA,Theta);
		Y+=y;
		SigmaY_emulator+=y*y;
	}
	SigmaY_emulator=SigmaY_emulator/double(NASample);
	Y=Y/double(NASample);
	SigmaY_emulator=sqrt(fabs(SigmaY_emulator-Y*Y));
}

void CSmoothEmulator::WriteCoefficients(){
	int isample,ic;
	FILE *fptr;
	string filename;
	string dirname=smoothmaster->CoefficientsDirName+"/"+observable_name;
	string command="mkdir -p "+dirname;
	system(command.c_str());
	filename=dirname+"/meta.txt";
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"# NPars  MaxRank NCoefficients\n");
	fprintf(fptr,"%d  %d  %d\n",NPars,smooth->MaxRank,smooth->NCoefficients);
	fclose(fptr);
	for(isample=0;isample<NASample;isample++){
		filename=dirname+"/sample"+to_string(isample)+".txt";
		fptr=fopen(filename.c_str(),"w");
		for(ic=0;ic<smooth->NCoefficients;ic++){
			fprintf(fptr,"%15.8e\n",ASample[isample][ic]);
		}
		fclose(fptr);
	}
}

void CSmoothEmulator::ReadCoefficients(){
	int isample,ic,NPars_test,MaxRank_test,NC_test;
	FILE *fptr;
	string filename;
	string dirname=smoothmaster->CoefficientsDirName+"/"+observable_name;
	string command="mkdir -p "+dirname;
	char dummy[100];
	system(command.c_str());
	filename=dirname+"/meta.txt";
	fptr=fopen(filename.c_str(),"r");
	fgets(dummy,100,fptr);
	fscanf(fptr,"%d  %d  %d\n",&NPars_test,&MaxRank_test,&NC_test);
	if(NPars_test!=NPars || MaxRank_test!=smooth->MaxRank || NC_test!=smooth->NCoefficients){
		CLog::Fatal("Mismatch in array sizes in ReadCoefficients");
	}
	fclose(fptr);
	for(isample=0;isample<NASample;isample++){
		filename=dirname+"/sample"+to_string(isample)+".txt";
		fptr=fopen(filename.c_str(),"r");
		for(ic=0;ic<smooth->NCoefficients;ic++){
			fscanf(fptr,"%lf\n",&ASample[isample][ic]);
		}
		fclose(fptr);
	}


}
