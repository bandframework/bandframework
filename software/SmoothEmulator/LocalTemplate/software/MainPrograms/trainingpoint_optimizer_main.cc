#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/trainingpoint_optimizer.h"
#include "msu_smoothutils/log.h"
#include "msu_smoothutils/misc.h"
using namespace std;
using namespace NMSUUtils;
int main(){
   NBandSmooth::CTPO *tpo=new NBandSmooth::CTPO();
   tpo->Optimize();
   tpo->WriteModelPars();
   return 0;
}
