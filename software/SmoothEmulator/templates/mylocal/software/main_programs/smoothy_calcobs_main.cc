#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smoothutils/log.h"
using namespace std;
int main(){
	NBandSmooth::CSmoothMaster master;
	//master.ReadCoefficients();
	master.TuneAllY();
	
	//Create model parameter object to store information about a single set of model parameters
	NBandSmooth::CModelParameters *modpars=new NBandSmooth::CModelParameters(); // contains info about single point
	modpars->priorinfo=master.priorinfo;
	// Print out the prior information
	master.priorinfo->PrintInfo();
	
	// Prompt user for model parameter values and enter them into the modpars object
	vector<double> X(modpars->NModelPars);
	for(unsigned int ipar=0;ipar<modpars->NModelPars;ipar++){
		cout << "Enter value for " << master.priorinfo->GetName(ipar) << ":\n";
		cin >> X[ipar];
	}
	modpars->SetX(X);
	
	//  Calc Observables
	NBandSmooth::CObservableInfo *obsinfo=master.observableinfo;
	vector<double> Y(obsinfo->NObservables);
	vector<double> SigmaY(obsinfo->NObservables);
	master.CalcAllY(modpars,Y,SigmaY);
	cout << "---- EMULATED OBSERVABLES ------\n";
	for(unsigned int iY=0;iY<obsinfo->NObservables;iY++){
		cout << obsinfo->GetName(iY) << " = " << Y[iY] << " +/- " << SigmaY[iY] << endl;
	}

	return 0;
}