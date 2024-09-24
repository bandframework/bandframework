#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smoothutils/log.h"
using namespace std;
int main(int argc,char *argv[]){
	NBandSmooth::CSmoothMaster master;
	master.TuneAllY();
	NBandSmooth::CModelParameters *modpars=new NBandSmooth::CModelParameters(); // contains info about single point
	modpars->priorinfo=master.priorinfo;
	//master.priorinfo->PrintInfo();
	unsigned int ipar,narg=argc-1,NModelPars=modpars->NModelPars;
	if(NModelPars!=narg){
		printf("NModelPars=%d, but you have %u arguments to the command line\n",NModelPars,argc-1);
		printf("Your input needs to have %d lines\n",NModelPars);
		exit(1);
	}
	
	
	// Set model parameter values
	vector<double> X(modpars->NModelPars);
	for(ipar=0;ipar<NModelPars;ipar++){
		X[ipar]=atof(argv[ipar+1]);
	}
	
	//for(unsigned int ipar=0;ipar<NModelPars;ipar++){
	//	printf("X[%d]=%g\n",ipar,X[ipar]);
	//}
	modpars->SetX(X);
	
	//  Calc Observables
	NBandSmooth::CObservableInfo *obsinfo=master.observableinfo;
	vector<double> Y(obsinfo->NObservables);
	vector<double> SigmaY(obsinfo->NObservables);
	master.CalcAllY(modpars,Y,SigmaY);
	//cout << "---- EMULATED OBSERVABLES ------\n";
	//for(unsigned int iY=0;iY<obsinfo->NObservables;iY++){
	//	cout << obsinfo->GetName(iY) << " = " << Y[iY] << " +/- " << SigmaY[iY] << endl;
	//}
	
	for(unsigned int iY=0;iY<obsinfo->NObservables;iY++)
		cout << Y[iY] << " ";
	cout << endl;
	for(unsigned int iY=0;iY<obsinfo->NObservables;iY++)
		cout << SigmaY[iY] << " ";
	cout << endl;
	

	return 0;
}
