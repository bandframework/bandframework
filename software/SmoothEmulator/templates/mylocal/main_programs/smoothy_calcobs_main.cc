#include "msu_commonutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_commonutils/log.h"

using namespace std;
int main(int argc,char *argv[]){
	if(argc!=2){
		CLog::Info("Usage smoothy_calcobs emulator parameter filename");
		exit(1);
	}
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(string(argv[1]));
	CSmoothMaster master(parmap);
	master.ReadCoefficientsAllY();
	
	CModelParameters *modpars=new CModelParameters(master.priorinfo); // contains info about single point
	
	master.priorinfo->PrintInfo();
	// Prompt user for model parameter values
	vector<double> X(modpars->NModelPars);
	for(int ipar=0;ipar<modpars->NModelPars;ipar++){
		cout << "Enter value for " << master.priorinfo->GetName(ipar) << ":\n";
		cin >> X[ipar];
	}
	modpars->SetX(X);
	
	//  Calc Observables
	CObservableInfo *obsinfo=master.observableinfo;
	vector<double> Y(obsinfo->NObservables);
	vector<double> SigmaY(obsinfo->NObservables);
	master.CalcAllY(modpars,Y,SigmaY);
	for(int iY=0;iY<obsinfo->NObservables;iY++){
		cout << obsinfo->GetName(iY) << " = " << Y[iY] << " +/- " << SigmaY[iY] << endl;
	}

	return 0;
}
