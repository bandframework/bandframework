#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smoothutils/log.h"

using namespace std;
int main(){
	NMSUUtils::CparameterMap *parmap=new NMSUUtils::CparameterMap();
	NBandSmooth::CSmoothMaster master(parmap);
	master.ReadCoefficientsAllY();
	master.ReadTrainingInfo();
	master.TestAtTrainingPts();
	return 0;
}
