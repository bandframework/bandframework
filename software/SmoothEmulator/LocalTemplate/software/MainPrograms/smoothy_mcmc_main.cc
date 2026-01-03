#include "msu_smoothutils/commonutils.h"
#include "msu_smooth/master.h"
#include "msu_smooth/mcmc.h"
using namespace std;
int main(){
	NBandSmooth::CSmoothMaster master;
	master.TuneAllY();
   //Next two lines can take replace of TuneAllY() if you have previously tuned and wish to speed calculations
   //master.ReadSigmaLambda();
   //master.TuneAllYFixedLambda();
	NBandSmooth::CMCMC mcmc(&master);
	
	unsigned int Nburn=mcmc.parmap->getI("MCMC_NBURN",1000);  // Steps for burn in
	unsigned int Ntrace=mcmc.parmap->getI("MCMC_NTRACE",1000); // Record this many points
	unsigned int Nskip=mcmc.parmap->getI("MCMC_NSKIP",5); // Only record every Nskip^th point
		
	mcmc.PerformTrace(1,Nburn);	
	CLog::Info("finished burn in\n");
	
	mcmc.PruneTrace(); // Throws away all but last point
   printf("Nburn=%u, Nskip=%u, Ntrace=%u\n",Nburn,Nskip,Ntrace);
	mcmc.PerformTrace(Ntrace,Nskip);
	mcmc.WriteTrace(); // Writes trace
	mcmc.EvaluateTrace();

	return 0;
}
