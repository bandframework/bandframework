#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUUtils;



void CSmoothEmulator::TuneMCMC(){
	if(!pca_ignore){
		vector<double> *Aswitch,*Aptr,*ATrialptr;
		double dlp,r,SigmaAswitch;
		unsigned int success=0,ic,imc;
		double BestLogP,stepsize;
		for(ic=0;ic<smooth->NCoefficients;ic++){
			ATrial[ic]=A[ic];
		}
		CalcAFromTraining(A);
		Aptr=&A;
		ATrialptr=&ATrial;
		SigmaATrial=SigmaA0;
		double logP,logPTrial;
		logP=GetLog_AProb(*Aptr,SigmaA);
		BestLogP=-1000000.0;
		for(imc=0;imc<NMC;imc++){

			SigmaATrial=fabs(SigmaA+SigmaA0*MCSigmaAStepSize*randy->ran_gauss());
			if(SigmaATrial<SigmaAMin){
				SigmaATrial=SigmaA;
			}
			for(ic=NTrainingPts;ic<smooth->NCoefficients;ic++){
				stepsize=MCStepSize*pow(LAMBDA,smooth->rank[ic]);
				(*ATrialptr)[ic]=SigmaATrial*( ((*Aptr)[ic]/SigmaA)+stepsize*randy->ran_gauss() );
			}
			//
			CalcAFromTraining(*ATrialptr);
			//Aptr=ATrialptr;
		
			//
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

		unsigned int Ndof=smooth->NCoefficients-NTrainingPts;
		if(!FirstTune)
			CLog::Info("success percentage="+to_string(double(success)*100.0/double(NMC))+",SigmaA="+to_string(SigmaA)+",  logP/Ndof="+to_string(logP/double(Ndof))+", BestLogP/Ndof="+to_string(BestLogP/double(Ndof))+"\n");
	}
}

void CSmoothEmulator::TuneMCMC_withSigma(){
	if(!pca_ignore){
		vector<double> *Aswitch,*Aptr,*ATrialptr;
		double dlp,r,SigmaAswitch,Y,stepsize;
		unsigned int success=0,ic,imc,itrain;
		double BestLogP;
		vector<double> YTrain=smoothmaster->traininginfo->YTrain[iY];
		vector<double> SigmaYTrain=smoothmaster->traininginfo->SigmaYTrain[iY];
		for(ic=0;ic<smooth->NCoefficients;ic++){
			ATrial[ic]=A[ic];
		}
		Aptr=&A;
		ATrialptr=&ATrial;
		SigmaATrial=SigmaA;
		double logP,logPTrial;
		logP=GetLog_AProb(*Aptr,SigmaA);
		for(itrain=0;itrain<NTrainingPts;itrain++){
			Y=smooth->CalcY_FromT(*Aptr,T[itrain]);
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
				Y=smooth->CalcY_FromT(*ATrialptr,T[itrain]);
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


		unsigned int Ndof=smooth->NCoefficients-NTrainingPts;
		if(!FirstTune)
			CLog::Info("success percentage="+to_string(double(success)*100.0/double(NMC))+", SigmaA="+to_string(SigmaA)+", logP/Ndof="+to_string(logP/double(Ndof))+",BestLogP/Ndof="+to_string(BestLogP/double(Ndof))+"\n");
	}
}

void CSmoothEmulator::TunePerfectMCMC(){
	if(!pca_ignore){
		unsigned int ic,ic0,ntry=0,ntrymax=100000;
		bool success=false;
		double weight,warg;//sigmafact=1.0;
		CalcTForTraining();
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
			CLog::Fatal("TunePerfectMCMC Failed, SigmaA0="+to_string(SigmaA0)+"\n");
		}
		else{
			for(ic=0;ic<smooth->NCoefficients;ic++){
				A[ic]=ATrial[ic];
			}
		}
	}
}

void CSmoothEmulator::GenerateASamples(){
	if(!pca_ignore){
		unsigned int isample,ic;
		FirstTune=true;
		for(ic=0;ic<smooth->NCoefficients;ic++)
			ABest[ic]=0.0;
		for(isample=0;isample<NASample;isample++){
			Tune();
			for(ic=0;ic<smooth->NCoefficients;ic++){
				ASample[isample][ic]=A[ic];
				ABest[ic]+=A[ic]/double(NASample);
			}
		}
		SigmaAbar+=SigmaA;
		NSigmaA+=1;
	}
}
