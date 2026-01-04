#include "msu_smoothutils/commonutils.h"
#include "msu_smooth/trainingpoint_optimizer.h"
using namespace std;
using namespace NMSUUtils;
int main(){
   NBandSmooth::CTPO *tpo=new NBandSmooth::CTPO();
   tpo->Optimize();
   tpo->WriteModelPars();
   return 0;
}
