#include "msu_smooth/testinginfo.h"

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;
CSmoothMaster* CTestingInfo::smoothmaster=NULL;

CTestingInfo::CTestingInfo(CObservableInfo *observableinfo_set,CPriorInfo *priorinfo_set){
	observableinfo=observableinfo_set;
	priorinfo=priorinfo_set;
	CModelParameters::priorinfo=priorinfo;
	NObservables=observableinfo->NObservables;
   ReadTestingInfo();
	unsigned int iy,ntest;
	YTest.resize(NObservables);
	for(iy=0;iy<NObservables;iy++){
		YTest[iy].resize(NTestingPts);
		for(ntest=0;ntest<NTestingPts;ntest++)
			YTest[iy][ntest]=0.0;
	}

	CSmoothEmulator::NTestingPts=NTestingPts;

}

void CTestingInfo::ReadTestingInfo(){
	unsigned int itest,ilist,ifile,iy,nsuccess=0,ipar,nread=0;
	unsigned int NObs=smoothmaster->observableinfo->NObservables;
	char filename[300],obs_charname[300],mod_par_name[300];
	string obs_name;
	double y,x;
	FILE *fptr;
	//
	string NTestingStr = smoothmaster->parmap->getS("SmoothEmulator_TestingPts","all");
	vector<unsigned int> NTestingList;
	NTestingList.clear();
	stringstream ss(NTestingStr);
	string token;
	string runfilename;
   
	int irun;
	bool exists;
   
	if(NTestingStr=="all" || NTestingStr=="All" || NTestingStr=="ALL"){
		irun=0;
		exists=false;
		do{
			runfilename=smoothmaster->FullModelTestingRunsDirName+"/run"+to_string(irun);
			if(filesystem::exists(runfilename)){
				NTestingList.push_back(irun);
				exists=true;
			}
			else
				exists=false;
			irun+=1;
		}while(exists);
	}
	else{
		while(getline(ss, token, ',')) {
			size_t pos = token.find("-");
			if (pos != string::npos){

				unsigned int start = stoi(token.substr(0, pos));
				unsigned int end = stoi(token.substr(pos+1));

				for (unsigned int i = start; i <= end; i++)
					NTestingList.push_back(i);
			}
			else {
				NTestingList.push_back(stoi(token));
			}
		}
	}
	NTestingPts=NTestingList.size();
   	
	modelpars.resize(NTestingPts);
	for(itest=0;itest<NTestingPts;itest++){
		modelpars[itest]=new CModelParameters();
	}
	YTest.resize(NObs);
	for(iy=0;iy<NObs;iy++){
		YTest[iy].resize(NTestingPts);
	}
   
	for(itest=0;itest<NTestingPts;itest++){
		ifile=NTestingList[itest];
		snprintf(filename,300,"%s/run%u/obs.txt",smoothmaster->FullModelTestingRunsDirName.c_str(),ifile);
		fptr=fopen(filename,"r");
		nsuccess=0;
		do{
			fscanf(fptr,"%s",obs_charname);
			if(!feof(fptr)){
				fscanf(fptr,"%lf",&y);
				obs_name=string(obs_charname);
				iy=smoothmaster->observableinfo->GetIPosition(obs_name);
				YTest[iy][itest]=y;
				nsuccess+=1;
			}
		}while(!feof(fptr));
		fclose(fptr);
	}

	if(nsuccess!=smoothmaster->observableinfo->NObservables)
		CLog::Fatal("In CTestingInfo::ReadTestInfo, only read in "+to_string(nsuccess)+" observables from file "+string(filename)+"\n");

	for(itest=0;itest<NTestingPts;itest++){
		ilist=NTestingList[itest];
		snprintf(filename,300,"%s/run%u/model_parameters.txt",smoothmaster->FullModelTestingRunsDirName.c_str(),ilist);
		fptr=fopen(filename,"r");
		nread=0;
		do{
			fscanf(fptr,"%s",mod_par_name);
			if(!feof(fptr)){
				fscanf(fptr,"%lf",&x);
				ipar=priorinfo->GetIPosition(mod_par_name);
				modelpars[itest]->X[ipar]=x;
				nread+=1;
			}
		}while(!feof(fptr));
		fclose(fptr);
	}
	if(nread!=priorinfo->NModelPars){
		CLog::Fatal("Only read in "+to_string(nread)+" parameter values from "+string(filename)+". But there are "+to_string(priorinfo->NModelPars)+" parameters needed.\n");
	}
	
	for(itest=0;itest<NTestingPts;itest++){
		modelpars[itest]->TranslateX_to_Theta();
	}

}
