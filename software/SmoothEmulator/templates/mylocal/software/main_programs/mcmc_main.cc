#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smoothutils/log.h"
#include "msu_smooth/mcmc.h"
using namespace std;
int main(){
	NBandSmooth::CSmoothMaster master;
	master.TuneAllY();
	//master.ReadTrainingInfo();
	NMSUUtils::CparameterMap *parmap=master.parmap;
	NBandSmooth::CMCMC mcmc(&master);
	
	unsigned int Nburn=parmap->getI("MCMC_NBURN",1000);  // Steps for burn in
	unsigned int Ntrace=parmap->getI("MCMC_NTRACE",1000); // Record this many points
	unsigned int Nskip=parmap->getI("MCMC_NSKIP",5); // Only record every Nskip^th point
		
	mcmc.PerformTrace(1,Nburn);	
	CLog::Info("finished burn in\n");
	
	mcmc.PruneTrace(); // Throws away all but last point
	mcmc.PerformTrace(Ntrace,Nskip);
	mcmc.WriteTrace();
	mcmc.WriteXTrace();
	mcmc.EvaluateTrace();

	return 0;
}
