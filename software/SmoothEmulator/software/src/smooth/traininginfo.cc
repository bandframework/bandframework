#include "msu_smooth/traininginfo.h"

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;
CSmoothMaster* CTrainingInfo::smoothmaster=NULL;


CTrainingInfo::CTrainingInfo(CObservableInfo *observableinfo_set,CPriorInfo *priorinfo_set){
	observableinfo=observableinfo_set;
	priorinfo=priorinfo_set;
	CModelParameters::priorinfo=priorinfo;
	NObservables=observableinfo->NObservables;
	printf("howdy a\n");
	printf("hmmm --- %s\n",smoothmaster->SmoothEmulator_TrainingFormat.c_str());
	
	if(smoothmaster->SmoothEmulator_TrainingFormat == "training_format_smooth"){
		printf("howdy howdy a\n");
		ReadTrainingInfoSmoothFormat();
				printf("howdy howdy b\n");
	}
	if(smoothmaster->SmoothEmulator_TrainingFormat == "training_format_surmise"){
		string TrainingInfoFileName=smoothmaster->parmap->getS("SmoothEmulator_TrainingInfoFileName","traininginfo.txt");
		printf("howdy howdy c\n");
		ReadTrainingInfoSurmiseFormat();
		printf("howdy howdy d\n");
	}
	else{
		CLog::Fatal("SmoothEmulator_TrainingFormat not recognized,\n should be training_format_smooth or training_format_surmise\n");
	}

	printf("howdy b\n");
	
	
	unsigned int iy,ntrain;
	YTrain.resize(NObservables);
	SigmaYTrain.resize(NObservables);
	for(iy=0;iy<NObservables;iy++){
		YTrain[iy].resize(NTrainingPts);
		SigmaYTrain[iy].resize(NTrainingPts);
		for(ntrain=0;ntrain<NTrainingPts;ntrain++)
			YTrain[iy][ntrain]=SigmaYTrain[iy][ntrain]=0.0;
	}

	CSmoothEmulator::NTrainingPts=NTrainingPts;

}


void CTrainingInfo::ReadTrainingInfoSmoothFormat(){
	
	if(smoothmaster->SmoothEmulator_TrainingFormat != "training_format_smooth"){
		CLog::Fatal("SmoothEmulator_TrainingFormat should be set to training_format_smooth\n if ReadTrainingInfo() is to be used\n");
	}
	//
	string NTrainingStr = smoothmaster->parmap->getS("SmoothEmulator_TrainingPts","1");
	vector<unsigned int> NTrainingList;
	stringstream ss(NTrainingStr);
	string token;
	string rundirname=smoothmaster->ModelRunDirName;

	while(getline(ss, token, ',')) {
		size_t pos = token.find("-");
		if (pos != string::npos) {

			unsigned int start = stoi(token.substr(0, pos));
			unsigned int end = stoi(token.substr(pos+1));

			for (unsigned int i = start; i <= end; i++)
				NTrainingList.push_back(i);
		}
		else {
			NTrainingList.push_back(stoi(token));
		}
	}
	//
	printf("howdy howdy howdy b\n");
	
	NTrainingPts = NTrainingList.size();
	printf("howdy c\n");
	modelpars.resize(NTrainingPts);
	for(ntrain=0;ntrain<NTrainingPts;ntrain++){
		modelpars[ntrain]=new CModelParameters();
	}
	printf("howdy d\n");
	
	unsigned int itrain,ilist,ifile,iy,nsuccess=0,ipar,nread;
	char filename[300],obs_charname[300],mod_par_name[300];
	string obs_name;
	double y,sigmay,x;
	FILE *fptr;
	
	YTrain.resize(NTrainingPts);
	SigmaYTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		YTrain[itrain].resize(smoothmaster->observableinfo->NObservables);
		SigmaYTrain[itrain].resize(smoothmaster->observableinfo->NObservables);
	}
		
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ifile=NTrainingList[itrain];
		if(smoothmaster->UsePCA){
			snprintf(filename,300,"smooth_data/%s/run%u/obs_pca.txt",rundirname.c_str(),ifile);
		}
		else{
			snprintf(filename,300,"smooth_data/%s/run%u/obs.txt",rundirname.c_str(),ifile);
		}
		printf("howdy howdy howdy filename=%s\n",filename);
		fptr=fopen(filename,"r");
		printf("howdy howdy howdy bb\n");
		nsuccess=0;
		do{
			printf("check z\n");
			fscanf(fptr,"%s",obs_charname);
			printf("obs_charname=%s\n",obs_charname);
			if(!feof(fptr)){
				fscanf(fptr,"%lf %lf",&y,&sigmay);
				obs_name=string(obs_charname);
				printf("y=%g, sigmay=%g\n",y,sigmay);
				iy=smoothmaster->observableinfo->GetIPosition(obs_name);
				printf("iy=%d, itrain=%d, YTrain.size=%lu\n",iy,itrain,YTrain.size());
				YTrain[iy][itrain]=y;
				printf("Ytrain=%g\n",YTrain[iy][itrain]);
				SigmaYTrain[iy][itrain]=sigmay;
				nsuccess+=1;
				printf("check xz\n");
			}
		}while(!feof(fptr));
		fclose(fptr);
	}
	printf("howdy howdy howdy c\n");
	
	if(nsuccess!=smoothmaster->observableinfo->NObservables)
		CLog::Fatal("In CTrainingInfo::ReadTrainInfo, only read in "+to_string(nsuccess)+" observables from file "+string(filename)+"\n");

	for(itrain=0;itrain<NTrainingPts;itrain++){
		ilist=NTrainingList[itrain];
		snprintf(filename,300,"smooth_data/%s/run%u/mod_parameters.txt",rundirname.c_str(),ilist);
		printf("---- filename=%s\n",filename);
		fptr=fopen(filename,"r");
		printf("modelpars.size=%lu\n",modelpars.size());
		nread=0;
		do{
			fscanf(fptr,"%s",mod_par_name);
			if(!feof(fptr)){
				fscanf(fptr,"%lf",&x);
				ipar=priorinfo->GetIPosition(mod_par_name);
				modelpars[itrain]->X[ipar]=x;
				nread+=1;
			}
		}while(!feof(fptr));
		fclose(fptr);
	}
	printf("howdy howdy howdy d\n");
	if(nread!=priorinfo->NModelPars){
		CLog::Fatal("Only read in "+to_string(nread)+" parameter values from "+string(filename)+". But there are "+to_string(priorinfo->NModelPars)+" parameters needed.\n");
	}
	printf("howdy howdy howdy e\n");

	for(itrain=0;itrain<NTrainingPts;itrain++){
		modelpars[itrain]->TranslateX_to_Theta();
	}
	printf("howdy howdy howdy f\n");

}



void CTrainingInfo::ReadTrainingInfoSurmiseFormat(){
	if(smoothmaster->SmoothEmulator_TrainingFormat != "training_format_smooth"){
		CLog::Fatal("SmoothEmulator_TrainingFormat should be set to training_format_smooth\n if ReadTrainingInfo(string rundirname) is to be used\n");
	}
	unsigned int itrain,ipar,nread1,nread2,iobs;
	unsigned int NModelPars=CModelParameters::NModelPars;
	unsigned int NObs=smoothmaster->observableinfo->NObservables;
	char dummy[100];
	string obs_name,filename;
	double y,x;
	
	
	filename="smooth_data/"+smoothmaster->TrainingThetasFileName;
	FILE *fptr=fopen(filename.c_str(),"r");
	itrain=0;
	nread1=0;
	do{
		for(ipar=0;ipar<NModelPars;ipar++){
			fscanf(fptr,"%lf",&x);
			if(!feof(fptr)){
				if(ipar==0){
					modelpars.push_back(NULL);
					modelpars[itrain]=new CModelParameters();
					nread1+=1;
				}
				modelpars[itrain]->X[ipar]=x;
			}
		}
		fgets(dummy,100,fptr);
		if(!feof(fptr))
			itrain+=1;
	}while(!feof(fptr));
	fclose(fptr);
	
	filename="smooth_data/"+smoothmaster->TrainingObsFileName;
	fptr=fopen(filename.c_str(),"r");
	itrain=0;
	nread2=0;
	do{
		for(iobs=0;iobs<NObs;iobs++){
			fscanf(fptr,"%lf",&y);
			if(!feof(fptr)){
				if(iobs==0){
					YTrain.resize(YTrain.size()+1);
					YTrain[itrain].resize(NObs);
					nread2+=1;
				}
				YTrain[itrain][itrain]=y;
			}
		}
		fgets(dummy,100,fptr);
		if(!feof(fptr))
			itrain+=1;
		
	}while(!feof(fptr));
	NTrainingPts=itrain;
	fclose(fptr);
	
	
	
	if(nread1!=nread2){
		CLog::Fatal("Read in "+to_string(nread1)+" training thetas but read in "+to_string(nread2)+" training observables\n");
	}

	for(itrain=0;itrain<NTrainingPts;itrain++){
		modelpars[itrain]->TranslateX_to_Theta();
	}

}
