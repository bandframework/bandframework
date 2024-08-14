#include "msu_smooth/master.h"
using namespace std;
using namespace NBandSmooth;

CSmoothMaster::CSmoothMaster(){
	unsigned int NObs;
	unsigned int iZ;
	CPCA *pca=NULL;
	parmap=new CparameterMap;
	parmap->ReadParsFromFile("smooth_data/smooth_parameters/emulator_parameters.txt");
	int ranseed=parmap->getI("RANDY_SEED",time(NULL));
	randy=new Crandy(ranseed);
	
	string logfilename=parmap->getS("SmoothEmulator_LogFileName","Screen");
	if(logfilename!="Screen"){
		CLog::Init(logfilename);
	}
	SmoothEmulator_TrainingFormat=parmap->getS("SmoothEmulator_TrainingFormat","training_format_smooth");
	string filename;
	UsePCA=parmap->getB("SmoothEmulator_UsePCA",false);
	if(UsePCA){
		filename="smooth_data/PCA_Info/observable_info.txt";
		CoefficientsDirName="coefficients_pca";
		pca=new CPCA();
	}
	else{
		filename="smooth_data/Info/observable_info.txt";
		CoefficientsDirName="coefficients";
	}
	observableinfo=new CObservableInfo(filename);
	NObs=observableinfo->NObservables;
	
	ModelRunDirName=parmap->getS("SmoothEmulator_ModelRunDirName","modelruns");
	TrainingThetasFileName=parmap->getS("SmoothEmulator_TrainingThetasFilename","TrainingThetas.txt");
	TrainingObsFileName=parmap->getS("SmoothEmulator_TrainingObsFilename","TrainingObs.txt");
	
	filename="smooth_data/Info/modelpar_info.txt";
	priorinfo=new CPriorInfo(filename);
	NPars=priorinfo->NModelPars;
	parmap->set("SmoothEmulator_NPars",NPars);
	parmap->set("Smooth_NPars",NPars);
	CTrainingInfo::smoothmaster=this;
	traininginfo = new CTrainingInfo(observableinfo,priorinfo);
	unsigned int maxrank=parmap->getI("SmoothEmulator_MAXRANK",4);
	smooth=new CSmooth(NPars,maxrank);

	CSmoothEmulator::NPars=NPars;
	CSmoothEmulator::smooth=smooth;
	CSmoothEmulator::smoothmaster=this;
	CSmoothEmulator::parmap=parmap;
	CSmoothEmulator::randy=randy;
	CSmoothEmulator::smooth=smooth;
	emulator.resize(NObs);
	
	pca_ignore.resize(NObs);
	if(UsePCA){
		pca->ReadTransformationInfo();
		pca_minvariance=parmap->getD("SmoothEmulator_PCAMinVariance",0.0);
		for(iZ=0;iZ<NObs;iZ++)
			pca_ignore[iZ]=false;
		for(iZ=0;iZ<NObs;iZ++){
			if(pca->eigvals(iZ)<pca_minvariance)
				pca_ignore[iZ]=true;
			else
				pca_ignore[iZ]=false;
		}
		delete pca;
	}
	else{
		for(iZ=0;iZ<NObs;iZ++)
			pca_ignore[iZ]=false;
	}
	for(unsigned int iy=0;iy<NObs;iy++){
		if(!UsePCA || !pca_ignore[iy]){
			emulator[iy]=new CSmoothEmulator(observableinfo->observable_name[iy],pca_ignore[iy]);
		}
	}
	ReadTrainingInfo();
	SetThetaTrain();
}

void CSmoothMaster::TuneAllY(){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		if(emulator[iY]->TuneChooseExact){
			if((UsePCA && !pca_ignore[iY]) || !UsePCA){
				CLog::Info("---- Tuning for "+observableinfo->observable_name[iY]+" ----\n");
				emulator[iY]->Tune();
			}
		}
		else{
			CLog::Info("Tuning Emulator for "+observableinfo->GetName(iY)+"\n");
			emulator[iY]->GenerateASamples();
		}
	}
}
	

void CSmoothMaster::TuneY(string obsname){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	emulator[iY]->Tune();
}

void CSmoothMaster::TuneY(unsigned int iY){
	emulator[iY]->Tune();
}

void CSmoothMaster::GenerateCoefficientSamples(){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		CLog::Info("Tuning Emulator for "+observableinfo->GetName(iY)+"\n");
		emulator[iY]->GenerateASamples();
	}
}

void CSmoothMaster::SetThetaTrain(){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->SetThetaTrain();
	}
}

void CSmoothMaster::CalcY(unsigned int iY,CModelParameters *modelpars,double &Y,double &SigmaY_emulator){
	emulator[iY]->CalcY(modelpars,Y,SigmaY_emulator);
}

void CSmoothMaster::CalcY(unsigned int iY,vector<double> &theta,double &Y,double &SigmaY_emulator){
	emulator[iY]->CalcY(theta,Y,SigmaY_emulator);
}

void CSmoothMaster::CalcYdYdTheta(string obsname,CModelParameters *modelpars,double &Y,
double &SigmaY_emulator,vector<double> &dYdTheta){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	emulator[iY]->CalcYDYDTheta(modelpars,Y,dYdTheta,SigmaY_emulator);
}

void CSmoothMaster::CalcYdYdTheta(unsigned int iY,CModelParameters *modelpars,double &Y,
double &SigmaY_emulator,vector<double> &dYdTheta){
	emulator[iY]->CalcYDYDTheta(modelpars,Y,dYdTheta,SigmaY_emulator);
}

void CSmoothMaster::CalcYdYdTheta(unsigned int iY,vector<double> &theta,double &Y,
double &SigmaY_emulator,vector<double> &dYdTheta){
	emulator[iY]->CalcYDYDTheta(theta,Y,dYdTheta,SigmaY_emulator);
}

void CSmoothMaster::CalcY(string obsname,CModelParameters *modelpars,double &Y,double &SigmaY_emulator){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	emulator[iY]->CalcY(modelpars,Y,SigmaY_emulator);
}

void CSmoothMaster::CalcY(string obsname,vector<double> &theta,double &Y,double &SigmaY_emulator){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	emulator[iY]->CalcY(theta,Y,SigmaY_emulator);
}

void CSmoothMaster::CalcAllY(CModelParameters *modelpars,vector<double> &Y,vector<double> &SigmaY_emulator){
	unsigned int NObservables=observableinfo->NObservables;
	Y.resize(NObservables);
	SigmaY_emulator.resize(NObservables);
	for(unsigned int iY=0;iY<NObservables;iY++){
		CalcY(iY,modelpars,Y[iY],SigmaY_emulator[iY]);
	}
}

void CSmoothMaster::CalcAllY(vector<double> &theta,vector<double> &Y,vector<double> &SigmaY_emulator){
	unsigned int NObservables=observableinfo->NObservables;
	Y.resize(NObservables);
	SigmaY_emulator.resize(NObservables);
	for(unsigned int iY=0;iY<NObservables;iY++){
		CalcY(iY,theta,Y[iY],SigmaY_emulator[iY]);
	}
}

void CSmoothMaster::CalcYOnly(string obsname,CModelParameters *modelpars,double &Y){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	emulator[iY]->CalcYOnly(modelpars,Y);
}

void CSmoothMaster::CalcYOnly(unsigned int iY,CModelParameters *modelpars,double &Y){
	emulator[iY]->CalcYOnly(modelpars,Y);
}

void CSmoothMaster::CalcYOnly(string obsname,vector<double> &theta,double &Y){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	emulator[iY]->CalcYOnly(theta,Y);
}

void CSmoothMaster::CalcYOnly(unsigned int iY,vector<double> &theta,double &Y){
	emulator[iY]->CalcYOnly(theta,Y);
}

void CSmoothMaster::CalcAllYOnly(CModelParameters *modelpars,vector<double> &Y){
	unsigned int NObservables=observableinfo->NObservables;
	Y.resize(NObservables);
	for(unsigned int iY=0;iY<NObservables;iY++){
		CalcYOnly(iY,modelpars,Y[iY]);
	}
}

void CSmoothMaster::CalcAllYOnly(vector<double> &theta,vector<double> &Y){
	unsigned int NObservables=observableinfo->NObservables;
	Y.resize(NObservables);
	for(unsigned int iY=0;iY<NObservables;iY++){
		CalcYOnly(iY,theta,Y[iY]);
	}
}

void CSmoothMaster::CalcAllYdYdTheta(CModelParameters *modelpars,vector<double> &Y,
vector<double> &SigmaY_emulator,vector<vector<double>> &dYdTheta){
	unsigned int NObservables=observableinfo->NObservables;
	Y.resize(NObservables);
	dYdTheta.resize(NObservables);
	SigmaY_emulator.resize(NObservables);
	for(unsigned int iY=0;iY<NObservables;iY++){
		CalcYdYdTheta(iY,modelpars,Y[iY],SigmaY_emulator[iY],dYdTheta[iY]);
		dYdTheta.resize(NPars);
	}
}

void CSmoothMaster::CalcAllYdYdTheta(vector<double> &theta,vector<double> &Y,
vector<double> &SigmaY_emulator,vector<vector<double>> &dYdTheta){
	unsigned int NObservables=observableinfo->NObservables;
	Y.resize(NObservables);
	dYdTheta.resize(NObservables);
	SigmaY_emulator.resize(NObservables);
	for(unsigned int iY=0;iY<NObservables;iY++){
		CalcYdYdTheta(iY,theta,Y[iY],SigmaY_emulator[iY],dYdTheta[iY]);
		dYdTheta.resize(NPars);
	}
}

void CSmoothMaster::CalcAllLogP(){
	unsigned int NObservables=observableinfo->NObservables;
	double logP,logPbar=0.0,A2barRatioBar=0.0,SigmaAbar=0.0;
	for(unsigned int iY=0;iY<NObservables;iY++){
		emulator[iY]->CalcExactLogP();
		logP=emulator[iY]->logP;
		CLog::Info("iY="+to_string(iY)+": logP="+to_string(logP)+", A2Ratio="+to_string(emulator[iY]->A2barRatio)+"\n");
		logPbar+=logP;
		A2barRatioBar+=emulator[iY]->A2barRatio;
		SigmaAbar+=emulator[iY]->SigmaA;
	}
	logPbar=logPbar/double(NObservables);
	A2barRatioBar=A2barRatioBar/double(NObservables);
	SigmaAbar=SigmaAbar/double(NObservables);
	//CLog::Info("logPbar="+to_string(logPbar)+", A2barRatio="+to_string(A2barRatioBar)+"\n");
	//CLog::Info(to_string(emulator[0]->LAMBDA)+" "+to_string(logPbar)+" "+to_string(A2barRatioBar)
	//	+" "+to_string(SigmaAbar)+"\n");
}

void CSmoothMaster::TestAtTrainingPts(){
	char pchars[CLog::CHARLENGTH];
	unsigned int itrain,iY;
	unsigned int NObservables=observableinfo->NObservables;
	double Y,SigmaY_emulator;
	CLog::Info("--- Y_train     Y_emulator    Sigma_emulator ----\n");
	for(itrain=0;itrain<traininginfo->NTrainingPts;itrain++){
		CLog::Info("------ itrain="+to_string(itrain)+" --------\n");
		for(iY=0;iY<NObservables;iY++){
			CalcY(iY,traininginfo->modelpars[itrain],Y,SigmaY_emulator);
			snprintf(pchars,CLog::CHARLENGTH,
			"Y[%u]=%10.3e =? %10.3e  +/- %12.5e\n",iY,traininginfo->YTrain[iY][itrain],Y,SigmaY_emulator);
			CLog::Info(pchars);
		}
	}
}

void CSmoothMaster::TestAtTrainingPts(unsigned int iY){
	char pchars[CLog::CHARLENGTH];
	unsigned int itrain;
	double Y,SigmaY_emulator;
	CLog::Info("--- Y_train     Y_emulator    Sigma_emulator ----\n");
	for(itrain=0;itrain<traininginfo->NTrainingPts;itrain++){
		CLog::Info("------ itrain="+to_string(itrain)+" --------\n");
		CalcY(iY,traininginfo->modelpars[itrain],Y,SigmaY_emulator);
		snprintf(pchars,CLog::CHARLENGTH,
		"Y[%u]=%10.3e =? %10.3e  +/- %12.5e\n",iY,traininginfo->YTrain[iY][itrain],Y,SigmaY_emulator);
		CLog::Info(pchars);
	}
}

void CSmoothMaster::TestAtTrainingPts(string obsname){
	char pchars[CLog::CHARLENGTH];
	unsigned int itrain,iY;
	double Y,SigmaY_emulator;
	iY=observableinfo->GetIPosition(obsname);
	CLog::Info("--- TESTING AT TRAINING POINTS ----\n");
	for(itrain=0;itrain<traininginfo->NTrainingPts;itrain++){
		CLog::Info("------ itrain="+to_string(itrain)+" --------\n");
		CalcY(iY,traininginfo->modelpars[itrain],Y,SigmaY_emulator);
		snprintf(pchars,CLog::CHARLENGTH,
		"Y[%u]=%10.3e =? %10.3e,    SigmaY=%12.5e\n",iY,traininginfo->YTrain[iY][itrain],Y,SigmaY_emulator);
		CLog::Info(pchars);
	}
}

void CSmoothMaster::TestVsFullModelAlt(){
	char pchars[CLog::CHARLENGTH];
	unsigned int iY,ipar,nfit=0,ntest=0;
	unsigned int NObservables=observableinfo->NObservables;
	double Y,SigmaY_emulator,realY;
	vector<double> testtheta;
	FILE *fptr,*fptr_out;
	string filename;
	for(iY=0;iY<NObservables;iY++){
		nfit=ntest=0;
		filename="smooth_data/fullmodel_testdata/"+observableinfo->observable_name[iY]+".txt";
		fptr=fopen(filename.c_str(),"r");
		filename="smooth_data/fullmodel_testdata/YvsY_"+observableinfo->observable_name[iY]+".txt";
		fptr_out=fopen(filename.c_str(),"w");
		
		testtheta.resize(NPars);
		do{
			for(ipar=0;ipar<NPars;ipar++){
				fscanf(fptr,"%lf",&testtheta[ipar]);
			}
			fscanf(fptr,"%lf",&realY);
			if(!feof(fptr)){
				ntest+=1;
				CalcY(iY,testtheta,Y,SigmaY_emulator);
				snprintf(pchars,CLog::CHARLENGTH,
				"Y[%u]=%10.3e =? %10.3e,    SigmaY_emulator=%12.5e\n",
				iY,Y,realY,SigmaY_emulator);
				fprintf(fptr_out,"%12.5e  %12.5e %12.5e\n",realY,Y,SigmaY_emulator);
				if(fabs(Y-realY)<SigmaY_emulator)
					nfit+=1;
			}	
		}while(!feof(fptr));
		fclose(fptr);		
		fclose(fptr_out);
		CLog::Info(observableinfo->observable_name[iY]+": "+to_string(nfit)+" out of "+to_string(ntest)+" points within 1 sigma\n");
	}
}

void CSmoothMaster::TestVsFullModel(){
	printf("check aa\n");
	string TestListStr = parmap->getS("SmoothEmulator_TestPts","1");
	
	vector<unsigned int> TestList;
	stringstream ss(TestListStr);
	string token;

	while(getline(ss, token, ',')) {
		size_t pos = token.find("-");
		if (pos != string::npos) {

			unsigned int start = stoi(token.substr(0, pos));
			unsigned int end = stoi(token.substr(pos+1));

			for (unsigned int i = start; i <= end; i++)
			TestList.push_back(i);
		}
		else {
			TestList.push_back(stoi(token));
		}
	}
	
	unsigned int ntestpts=TestList.size();
	char obsnamechars[200],modparnamechars[200];
	string obsname,modparname;
	unsigned int iY,iread,itest,ipar,nfit;
	unsigned int NObservables=observableinfo->NObservables;
	double Y,SigmaY_emulator,realY,realSigmaY;
	double Xread,SigmaXRead;
	CModelParameters testpars;
	FILE *fptr,*fptr_out;
	string filename;
	for(iY=0;iY<NObservables;iY++){
		filename="smooth_data/fullmodel_testdata/YvsY_"+observableinfo->observable_name[iY]+".txt";
		fptr_out=fopen(filename.c_str(),"w");
		nfit=0;
		
		//CLog::Info("Writing test_vs_full_model results to "+filename+"\n");
		
		for(itest=0;itest<ntestpts;itest++){
			printf("check aaa\n");
			filename="smooth_data/modelruns/run"+to_string(TestList[itest])+"/mod_parameters.txt";
			printf("filename=%s\n",filename.c_str());
			fptr=fopen(filename.c_str(),"r");
			for(iread=0;iread<NPars;iread++){
				fscanf(fptr,"%s %lf %lf",modparnamechars,&Xread,&SigmaXRead);
				modparname=string(modparnamechars);
				ipar=priorinfo->GetIPosition(modparname);
				testpars.X[ipar]=Xread;
			}
			printf("check bbb\n");
			fclose(fptr);
			testpars.TranslateX_to_Theta();
			CalcY(iY,testpars.Theta,Y,SigmaY_emulator);
			
			filename="smooth_data/modelruns/run"+to_string(TestList[itest])+"/obs.txt";
			printf("check ccc filename=%s\n",filename.c_str());
			fptr=fopen(filename.c_str(),"r");
			printf("check cccc\n");
			iread=-1;
			do{
				iread+=1;
				fscanf(fptr,"%s %lf %lf",obsnamechars,&realY,&realSigmaY);
				obsname=string(obsnamechars);
				
			}while(iread<observableinfo->NObservables && obsname!=observableinfo->observable_name[iY]);
			printf("check ddd\n");
			fclose(fptr);
			if(fabs(Y-realY)<SigmaY_emulator)
				nfit+=1;
			if(obsname!=observableinfo->observable_name[iY])
				CLog::Fatal("cannot find obsname amongst observales, obsname="+obsname+"\n");
			
			fprintf(fptr_out,"%lf %lf %lf\n",realY,Y,SigmaY_emulator);
			
		}
		fclose(fptr_out);
		CLog::Info(observableinfo->observable_name[iY]+": "+to_string(nfit)+" out of "+to_string(ntestpts)+" points within 1 sigma\n");
	}
}

