#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smoothutils/log.h"

using namespace std;
int main(){
	NBandSmooth::CSmoothMaster master;
   printf("master created\n");
	master.TuneAllY();
   printf("tuning finished\n");
	master.TestAtTrainingPts();
	return 0;
}
