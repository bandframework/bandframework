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
	if(smoothmaster->SmoothEmulator_TrainingFormat == "training_format_smooth"){
		ReadTrainingInfoSmoothFormat();
	}
	else if(smoothmaster->SmoothEmulator_TrainingFormat == "training_format_surmise"){
		string TrainingInfoFileName=smoothmaster->parmap->getS("SmoothEmulator_TrainingInfoFileName","traininginfo.txt");
		ReadTrainingInfoSurmiseFormat();
	}
	else{
		CLog::Fatal("SmoothEmulator_TrainingFormat not recognized,\n should be training_format_smooth or training_format_surmise\n");
	}

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
	unsigned int itrain,ilist,ifile,iy,nsuccess=0,ipar,nread;
	char filename[300],obs_charname[300],mod_par_name[300];
	string obs_name;
	double y,sigmay,x;
	FILE *fptr;
	
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

	NTrainingPts = NTrainingList.size();
	modelpars.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		modelpars[itrain]=new CModelParameters();
	}
	
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
		fptr=fopen(filename,"r");
		nsuccess=0;
		do{
			fscanf(fptr,"%s",obs_charname);
			if(!feof(fptr)){
				fscanf(fptr,"%lf %lf",&y,&sigmay);
				obs_name=string(obs_charname);
				iy=smoothmaster->observableinfo->GetIPosition(obs_name);
				YTrain[iy][itrain]=y;
				SigmaYTrain[iy][itrain]=sigmay;
				nsuccess+=1;
			}
		}while(!feof(fptr));
		fclose(fptr);
	}

	if(nsuccess!=smoothmaster->observableinfo->NObservables)
		CLog::Fatal("In CTrainingInfo::ReadTrainInfo, only read in "+to_string(nsuccess)+" observables from file "+string(filename)+"\n");

	for(itrain=0;itrain<NTrainingPts;itrain++){
		ilist=NTrainingList[itrain];
		snprintf(filename,300,"smooth_data/%s/run%u/mod_parameters.txt",rundirname.c_str(),ilist);
		fptr=fopen(filename,"r");
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
	if(nread!=priorinfo->NModelPars){
		CLog::Fatal("Only read in "+to_string(nread)+" parameter values from "+string(filename)+". But there are "+to_string(priorinfo->NModelPars)+" parameters needed.\n");
	}
	for(itrain=0;itrain<NTrainingPts;itrain++){
		modelpars[itrain]->TranslateX_to_Theta();
		modelpars[itrain]->Print();
	}

}

void CTrainingInfo::ReadTrainingInfoSurmiseFormat(){
	if(smoothmaster->SmoothEmulator_TrainingFormat != "training_format_surmise"){
		CLog::Fatal("SmoothEmulator_TrainingFormat should be set to training_format_surmise\n if ReadTrainingInfoSurmiseFormat() is to be used\n");
	}
	printf("howdy\n");
	unsigned int itrain,ipar,iobs;
	unsigned int NModelPars=CModelParameters::NModelPars;
	unsigned int NObs=smoothmaster->observableinfo->NObservables;
	char dummy[10000];
	string obs_name,filename;
	double y,x;
	
	filename="smooth_data/"+smoothmaster->TrainingThetasFileName;
	FILE *fptr=fopen(filename.c_str(),"r");
	itrain=0;
	do{
		for(ipar=0;ipar<NModelPars;ipar++){
			fscanf(fptr,"%lf",&x);
			if(!feof(fptr)){
				if(ipar==0){
					modelpars.push_back(NULL);
					modelpars[itrain]=new CModelParameters();
				}
				modelpars[itrain]->X[ipar]=x;
				//printf("X[%d][%d]=%g\n",itrain,ipar,modelpars[itrain]->X[ipar]);
			}
		}
		fgets(dummy,10000,fptr);
		if(!feof(fptr))
			itrain+=1;
	}while(!feof(fptr));
	fclose(fptr);
	
	NTrainingPts=itrain;
	filename="smooth_data/"+smoothmaster->TrainingObsFileName;
	printf("ntrain=%d, NObs=%d, filename=%s\n",NTrainingPts,NObs,filename.c_str());
	fptr=fopen(filename.c_str(),"r");
	
	YTrain.resize(NObs);
	SigmaYTrain.resize(NObs);
	for(iobs=0;iobs<NObs;iobs++){
		YTrain[iobs].resize(NTrainingPts);
		SigmaYTrain[iobs].resize(NTrainingPts);
	}
	
	printf("check a\n");
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(iobs=0;iobs<NObs;iobs++){
			fscanf(fptr,"%lf",&y);
			if(feof(fptr)){
				CLog::Fatal("reading training info: not enough lines in "+smoothmaster->TrainingObsFileName+"\n");
			}
			YTrain[iobs][itrain]=y;
			SigmaYTrain[iobs][itrain]=0.0;
		}
		fgets(dummy,10000,fptr);
	}
	
	printf("NTrainingPts=%d\n",NTrainingPts);
	fclose(fptr);
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		modelpars[itrain]->TranslateX_to_Theta();
		modelpars[itrain]->Print();
	}

}
