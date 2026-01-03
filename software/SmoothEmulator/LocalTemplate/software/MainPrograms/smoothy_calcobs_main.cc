#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smoothutils/log.h"
using namespace std;
int main(){
   unsigned int iY,ipar;
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
	for(ipar=0;ipar<modpars->NModelPars;ipar++){
      printf("Enter value for %s: ",(master.priorinfo->GetName(ipar)).c_str());
      scanf("%lf",&X[ipar]);
	}
	modpars->SetX(X);
	
	//  Calc Observables
	NBandSmooth::CObservableInfo *obsinfo=master.observableinfo;
	vector<double> Y(obsinfo->NObservables);
	vector<double> SigmaY(obsinfo->NObservables);
	master.GetAllY(modpars,Y,SigmaY);
	printf("---- EMULATED OBSERVABLES ------\n");
	for(iY=0;iY<obsinfo->NObservables;iY++){
      printf("%10s = %11.4e +/- %11.4e\n",
             (obsinfo->GetName(iY)).c_str(),Y[iY],SigmaY[iY]);
	}

	return 0;
}
