#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/simplex.h"
#include "msu_smoothutils/log.h"
#include "msu_smoothutils/misc.h"
using namespace std;
int main(){
	NMSUUtils::CparameterMap *parmap=new NMSUUtils::CparameterMap();
	NBandSmooth::CSimplexSampler *simplex=new NBandSmooth::CSimplexSampler(parmap);
	simplex->SetThetaSimplex();
	simplex->WriteModelPars();
	return 0;
}
