#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smoothutils/log.h"
#include "msu_smooth/mcmc.h"
using namespace std;
int main(){
	NMSUUtils::CparameterMap *parmap=new CparameterMap();
	NBandSmooth::CSmoothMaster master(parmap);	
	NBandSmooth::CMCMC mcmc(&master);
	master.ReadCoefficientsAllY();
	master.ReadTrainingInfo();
	//master.TestAtTrainingPts();
	
	// Prompt user for model parameter values
	vector<double> theta(master.priorinfo->NModelPars);
	for(unsigned int ipar=0;ipar<master.priorinfo->NModelPars;ipar++){
		theta[ipar]=0.2;
	}
	
	//  Calc Observables
	NBandSmooth::CObservableInfo *obsinfo=master.observableinfo;
	vector<double> Y(obsinfo->NObservables);
	vector<double> SigmaY(obsinfo->NObservables);
	master.CalcAllY(theta,Y,SigmaY);
	cout << "---- EMULATED OBSERVABLES ------\n";
	for(unsigned int iY=0;iY<obsinfo->NObservables;iY++){
		cout << obsinfo->GetName(iY) << " = " << Y[iY] << " +/- " << SigmaY[iY] << endl;
	}
	
	unsigned int Nburn=parmap->getI("MCMC_NBURN",1000);  // Steps for burn in
	unsigned int Ntrace=parmap->getI("MCMC_NTRACE",1000); // Record this many points
	unsigned int Nskip=parmap->getI("MCMC_NSKIP",5); // Only record every Nskip^th point
		
	mcmc.PerformTrace(1,Nburn);	
	CLog::Info("finished burn in\n");
	
	mcmc.PruneTrace(); // Throws away all but last point
	mcmc.PerformTrace(Ntrace,Nskip);
	mcmc.WriteTrace();
	mcmc.EvaluateTrace();

	return 0;
}
