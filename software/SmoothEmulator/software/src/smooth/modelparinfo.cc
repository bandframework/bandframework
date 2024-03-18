#include "msu_smooth/modelparinfo.h"
NBandSmooth::CPriorInfo* NBandSmooth::CModelParameters::priorinfo=NULL;
double NBandSmooth::CModelParameters::GSCALE=sqrt(3.0);
unsigned int NBandSmooth::CModelParameters::NModelPars=0;

using namespace std;

using namespace NBandSmooth;

CModelParameters::CModelParameters(){
	X.resize(NModelPars);
	Theta.resize(NModelPars);
};

void CModelParameters::TranslateX_to_Theta(){
	//for min, max range
	double sigmax,xbar;
	unsigned int ipar;

	for(ipar=0;ipar<NModelPars;ipar++){
		if(priorinfo->type[ipar]=="uniform"){
			Theta[ipar]=-1+2*((X[ipar]-priorinfo->xmin[ipar])/(priorinfo->xmax[ipar]-priorinfo->xmin[ipar]));
		}
		else if(priorinfo->type[ipar]=="gaussian"){
			xbar=priorinfo->xmin[ipar];
			sigmax=priorinfo->xmax[ipar];
			Theta[ipar]=(X[ipar]-xbar)/(sigmax*GSCALE);
		}
		else{
			CLog::Fatal("Cannot translate X to Theta because type = "+priorinfo->type[ipar]+" is not recognized\n");
		}
	}
}

void CModelParameters::TranslateTheta_to_X(){
	double sigmax,xbar;
	unsigned int ipar;

	for(ipar=0;ipar<NModelPars;ipar++){
		if(priorinfo->type[ipar]=="uniform"){
			X[ipar]=priorinfo->xmin[ipar]+0.5*(1.0+Theta[ipar])*(priorinfo->xmax[ipar]-priorinfo->xmin[ipar]);
		}
		else if(priorinfo->type[ipar]=="gaussian"){
			xbar=priorinfo->xmin[ipar];
			sigmax=priorinfo->xmax[ipar];
			X[ipar]=xbar+GSCALE*sigmax*Theta[ipar];
		}
		else{
			CLog::Fatal("Cannot translate Theta to X because type = "+priorinfo->type[ipar]+" is not recognized\n");
		}
	}

}

void CModelParameters::Print(){
	char message[CLog::CHARLENGTH];
	unsigned int ipar;
	CLog::Info("--------------------------------------\n");
	for(ipar=0;ipar<NModelPars;ipar++){
		snprintf(message,CLog::CHARLENGTH," %24s (%8s): x=%11.4e, theta=%11.4e\n",
		priorinfo->parname[ipar].c_str(),priorinfo->type[ipar].c_str(),X[ipar],Theta[ipar]);
		CLog::Info(message);
	}
}

void CModelParameters::Write(string filename){
	unsigned int ipar;
	string command;
	command="mkdir -p mcmc_trace";
	system(command.c_str());
	filename="mcmc_trace/XbarThetabar.txt";
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#   parname        X      Theta\n");
	for(ipar=0;ipar<NModelPars;ipar++){
		fprintf(fptr,"%24s %12.5e, theta=%12.5e\n",
		priorinfo->parname[ipar].c_str(),X[ipar],Theta[ipar]);
	}
	fclose(fptr);
}

void CModelParameters::SetX(vector<double> &x){
	for(unsigned int ipar=0;ipar<NModelPars;ipar++){
		X[ipar]=x[ipar];
	}
	TranslateX_to_Theta();
}

void CModelParameters::SetTheta(vector<double> &theta){
	for(unsigned int ipar=0;ipar<NModelPars;ipar++){
		Theta[ipar]=theta[ipar];
	}
	TranslateTheta_to_X();
}

void CModelParameters::Copy(CModelParameters *mp){
	X=mp->X;
	Theta=mp->Theta;
}
