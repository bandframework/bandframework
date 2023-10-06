#include "msu_smooth/modelparinfo.h"

using namespace std;
using namespace NBandSmooth;

CModelParameters::CModelParameters(CPriorInfo *priorinfo_set){
	priorinfo=priorinfo_set;
	NModelPars=priorinfo->NModelPars;
	X.resize(NModelPars);
	Theta.resize(NModelPars);
};

void CModelParameters::TranslateX_to_Theta(){
	//for min, max range
	double sigmax,xbar;
	int ipar;

	for(ipar=0;ipar<NModelPars;ipar++){
		if(priorinfo->type[ipar]=="uniform"){
			Theta[ipar]=-1+2*((X[ipar]-priorinfo->xmin[ipar])/(priorinfo->xmax[ipar]-priorinfo->xmin[ipar]));
		}
		else if(priorinfo->type[ipar]=="gaussian"){
			xbar=priorinfo->xmin[ipar];
			sigmax=priorinfo->xmax[ipar];
			Theta[ipar]=(X[ipar]-xbar)/sigmax;
		}
		else{
			CLog::Fatal("Cannot translate X to Theta because type = "+priorinfo->type[ipar]+" is not recognized\n");
		}
	}
}

void CModelParameters::TranslateTheta_to_X(){
	double sigmax,xbar;
	int ipar;

	for(ipar=0;ipar<NModelPars;ipar++){
		if(priorinfo->type[ipar]=="uniform"){
			X[ipar]=priorinfo->xmin[ipar]+0.5*(1.0+Theta[ipar])*(priorinfo->xmax[ipar]-priorinfo->xmin[ipar]);
		}
		else if(priorinfo->type[ipar]=="gaussian"){
			xbar=priorinfo->xmin[ipar];
			sigmax=priorinfo->xmax[ipar];
			X[ipar]=xbar+sigmax*Theta[ipar];
		}
		else{
			CLog::Fatal("Cannot translate Theta to X because type = "+priorinfo->type[ipar]+" is not recognized\n");
		}
	}

}

void CModelParameters::Print(){
	char message[200];
	int ipar;
	for(ipar=0;ipar<NModelPars;ipar++){
		snprintf(message,200,"   %.24s (%.8s): x=%11.4e, theta=%11.4e\n",
		priorinfo->parname[ipar].c_str(),priorinfo->type[ipar].c_str(),X[ipar],Theta[ipar]);
		CLog::Info(message);
	}
}

void CModelParameters::SetX(vector<double> &x){
	for(int ipar=0;ipar<NModelPars;ipar++){
		X[ipar]=x[ipar];
	}
	TranslateX_to_Theta();
}
