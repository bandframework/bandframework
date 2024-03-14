#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/smooth.h"
#include "msu_smoothutils/log.h"
#include <filesystem>

using namespace std;
using namespace NMSUUtils;
int main(){

	const unsigned int NObs=10,NPars=10;
	unsigned int NTrain,itrain,iobs,ic,ipar,maxrank=5,it=time(NULL),itest;
	double LAMBDA;
	double y,x;
	vector<double> A;
	vector<double> theta;
	vector<vector<double>> thetatest;
	theta.resize(NPars);
	NMSUUtils::Crandy randy(0);
	
	printf("Enter LAMBDA: ");
	scanf("%lf",&LAMBDA);
	//printf("Enter NTrain: ");
	//scanf("%u",&NTrain);
	NTrain=0;
	bool existence;
	do{
		string filename="modelruns/run"+to_string(NTrain);
		filesystem::path f{filename};
		existence=filesystem::exists(f);
		if(existence){
			NTrain+=1;
		}
	}while(existence);
	CLog::Info("NTraining Pts="+to_string(NTrain)+"\n");
	CLog::Info("NPars="+to_string(NPars)+"\n");
	
	NBandSmooth::CSmooth smooth(NPars,maxrank);
	A.resize(smooth.NCoefficients);
	
	string obsname[NObs];
	for(iobs=0;iobs<NObs;iobs++){
		obsname[iobs]="obs"+to_string(iobs+1);
	}
	string parname[NPars];
	char parname_c[200];
	string filename;
	FILE *fptr,*fptr_out;
	
	
	// Now make some data for later testing, not for tuning
	
	thetatest.resize(100);
	for(itest=0;itest<100;itest++){
		thetatest[itest].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++){
			thetatest[itest][ipar]=-1.0+2.0*randy.ran();
		}
	}

	
	string command="mkdir -p realdata";
	system(command.c_str());
	
	
	for(iobs=0;iobs<NObs;iobs++){
		randy.reset(iobs+it);
		for(ic=0;ic<smooth.NCoefficients;ic++){
			A[ic]=100.0*randy.ran_gauss();
		}
		for(itrain=0;itrain<NTrain;itrain++){
			filename="modelruns/run"+to_string(itrain)+"/mod_parameters.txt";
			fptr=fopen(filename.c_str(),"r");
			for(ipar=0;ipar<NPars;ipar++){
				fscanf(fptr,"%s %lf",parname_c,&x);
				parname[ipar]=parname_c;
				theta[ipar]=-1.0+x/50.0;
			}
			fclose(fptr);
			y=smooth.CalcY(A,LAMBDA,theta);
			filename="modelruns/run"+to_string(itrain)+"/obs.txt";
			if(iobs==0)
				fptr_out=fopen(filename.c_str(),"w");
			else
				fptr_out=fopen(filename.c_str(),"a");
			fprintf(fptr_out,"%s %g 0.0\n",obsname[iobs].c_str(),y);
			fclose(fptr_out);
		}
		
		filename="realdata/"+obsname[iobs]+".txt";
		fptr=fopen(filename.c_str(),"w");
		for(itest=0;itest<100;itest++){
			y=smooth.CalcY(A,LAMBDA,thetatest[itest]);
			for(ipar=0;ipar<NPars;ipar++){
				fprintf(fptr,"%17.10e ",thetatest[itest][ipar]);
			}
			fprintf(fptr,"%17.10e\n",y);
		}
		fclose(fptr);
	}
	
	
	
	return 0;
}
