#include "msu_smooth/traininginfo.h"

CSmoothMaster* CTrainingInfo::smoothmaster=NULL;

using namespace std;

CTrainingInfo::CTrainingInfo(vector<int> NTrainingList_set,CObservableInfo *observableinfo_set,CPriorInfo *priorinfo_set){
	observableinfo=observableinfo_set;
	priorinfo=priorinfo_set;
	NObservables=observableinfo->NObservables;
	NTrainingList = NTrainingList_set;
	NTrainingPts = NTrainingList.size();
	int iy,ntrain;
	YTrain.resize(NObservables);
	SigmaYTrain.resize(NObservables);
	for(iy=0;iy<NObservables;iy++){
		YTrain[iy].resize(NTrainingPts);
		SigmaYTrain[iy].resize(NTrainingPts);
		for(ntrain=0;ntrain<NTrainingPts;ntrain++)
			YTrain[iy][ntrain]=SigmaYTrain[iy][ntrain]=0.0;
	}

	modelpars.resize(NTrainingPts);
	for(ntrain=0;ntrain<NTrainingPts;ntrain++){
		modelpars[ntrain]=new CModelParameters(smoothmaster->priorinfo);
	}

}

void CTrainingInfo::ReadTrainingInfo(string rundirname){
	int itrain,ifile,iy,nsuccess=0,ipar,nread;
	char filename[300],obs_charname[300],mod_par_name[300];
	string obs_name;
	double y,sigmay,x;
	FILE *fptr;
	for(itrain=0;itrain<NTrainingList.size();itrain++){
		ifile=NTrainingList[itrain];
		if(smoothmaster->UsePCA){
			snprintf(filename,300,"%s/run%d/obs_PCA.txt",rundirname.c_str(),ifile);
		}
		else{
			snprintf(filename,300,"%s/run%d/obs.txt",rundirname.c_str(),ifile);
		}
		fptr=fopen(filename,"r");
		nsuccess=0;
		do{
			fscanf(fptr,"%s %lf %lf",obs_charname,&y,&sigmay);
			if(!feof(fptr)){
				obs_name=string(obs_charname);
				iy=smoothmaster->observableinfo->GetIPosition(obs_name);
				YTrain[iy][itrain]=y;
				SigmaYTrain[iy][itrain]=sigmay;
				nsuccess+=1;
			}
		}while(!feof(fptr));
	}
	fclose(fptr);
	if(nsuccess!=smoothmaster->observableinfo->NObservables)
		CLog::Fatal("In CTrainingInfo::ReadTrainInfo, only read in "+to_string(nsuccess)+" observables from file "+string(filename)+"\n");

	for(itrain=0;itrain<NTrainingPts;itrain++){
		snprintf(filename,300,"%s/run%d/mod_parameters.txt",rundirname.c_str(),itrain);
		fptr=fopen(filename,"r");
		nread=0;
		do{
			fscanf(fptr,"%s %lf",mod_par_name,&x);
			if(!feof(fptr)){
				ipar=priorinfo->GetIPosition(mod_par_name);
				modelpars[itrain]->X[ipar]=x;
				nread+=1;
			}
		}while(!feof(fptr));
	}
	if(nread!=priorinfo->NModelPars){
		CLog::Fatal("Only read in "+to_string(nread)+" parameter values from "+string(filename)+". But there are "+to_string(priorinfo->NModelPars)+" parameters needed.\n");
	}

	for(itrain=0;itrain<NTrainingPts;itrain++){
		modelpars[itrain]->TranslateX_to_Theta();
	}

}
