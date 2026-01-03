#include "msu_smoothutils/commonutils.h"
#include "msu_smooth/master.h"

using namespace std;
int main(){
	NBandSmooth::CSmoothMaster master;
	master.TuneAllY();
	master.TestVsFullModel();
	return 0;
}
