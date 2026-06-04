#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smoothutils/log.h"

using namespace std;
int main(){
	NBandSmooth::CSmoothMaster master;
	master.TuneAllY();
	//master.WriteCoefficientsAllY();
	master.TestAtTrainingPts();
	return 0;
}
