#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smoothutils/log.h"

using namespace std;
using namespace  NBandSmooth;
using namespace NMSUUtils;

int main(){
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile("parameters/emulator_parameters.txt");

	CSmoothMaster master(parmap);
	
	master.ReadTrainingInfo();
	
	master.TuneAllY();
	
	//master.CalcAllLogP();
	
	//master.TestAtTrainingPts();
	
	master.WriteCoefficientsAllY();

	return 0;
}
