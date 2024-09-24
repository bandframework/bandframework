#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smoothutils/log.h"

using namespace std;
int main(){
	NBandSmooth::CSmoothMaster master;
	double Lambda;
	unsigned int iY,NObs;
	NObs=master.observableinfo->NObservables;

	FILE *fptr1=fopen("LogPvslambda.txt","w");
	FILE *fptr2=fopen("SigmaAvslambda.txt","w");
	for(Lambda=1.5;Lambda<8.01;Lambda+=0.25){
		
		for(iY=0;iY<NObs;iY++){
			master.emulator[iY]->LAMBDA=Lambda;
		}
		master.TuneAllY();
		fprintf(fptr1,"%5.2f ",Lambda);
		fprintf(fptr2,"%5.2f ",Lambda);
		for(iY=0;iY<NObs;iY++){
			fprintf(fptr1,"%10.3e ",master.emulator[iY]->logP);
			fprintf(fptr2,"%10.3e ",master.emulator[iY]->SigmaA);
		}
		fprintf(fptr1,"\n");
		fprintf(fptr2,"\n");
	}
	fclose(fptr1);
	fclose(fptr2);
	
	return 0;
}
