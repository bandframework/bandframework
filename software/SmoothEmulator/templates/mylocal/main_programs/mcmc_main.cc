#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smoothutils/log.h"
#include "msu_smooth/mcmc.h"

using namespace std;
using namespace NMSUUtils;
using namespace NBandSmooth;

int main(){
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(string("parameters/emulator_parameters.txt"));
	parmap->ReadParsFromFile(string("parameters/mcmc_parameters.txt"));
	CSmoothMaster master(parmap);	
	CMCMC mcmc(&master);
	master.ReadCoefficientsAllY();
	master.ReadTrainingInfo();
	//master.TestAtTrainingPts();
	
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
