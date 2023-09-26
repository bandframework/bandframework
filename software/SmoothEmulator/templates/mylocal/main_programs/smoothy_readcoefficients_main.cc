#include "msu_commonutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_commonutils/log.h"

using namespace std;
int main(int argc,char *argv[]){
	if(argc!=2){
		printf("Usage smoothy emulator parameter filename");
		exit(1);
	}
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(string(argv[1]));
	
	CSmoothMaster master(parmap);
	master.randy->reset(-time(NULL));
	
	master.ReadCoefficientsAllY();
	
	master.ReadTrainingInfo();

	master.TestAtTrainingPts();

	return 0;
}
