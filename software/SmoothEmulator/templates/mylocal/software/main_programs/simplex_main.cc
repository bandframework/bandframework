#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/simplex.h"
#include "msu_smoothutils/log.h"
#include "msu_smoothutils/misc.h"
using namespace std;
int main(){
	NBandSmooth::CSimplexSampler *simplex=new NBandSmooth::CSimplexSampler();
	simplex->SetThetaSimplex();
	simplex->WriteModelPars();
	return 0;
}
