#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/smooth.h"
#include "msu_smoothutils/log.h"
#include <filesystem>

using namespace std;
using namespace NMSUUtils;
int main(){

	unsigned int NObs=0,NPars=0;
	unsigned int NTrain,itrain,iy,ic,ipar,maxrank=4,it=time(NULL),itest,Ntest=50;
	double LAMBDA; 
	double y,x,xminread,xmaxread;
	vector<double> A;
	vector<double> theta,exptheta;
	vector<string> priortype;
	vector<vector<double>> thetatest;
	char parname_c[200],dummy1[200],dummy2[100];
	string expfilename,filename;
	vector<string> obsname;
	vector<string> parname;
	vector<double> xmin,xmax;
	FILE *fptr,*fptr_out,*expfptr;
	NMSUUtils::Crandy randy(123);
	
	printf("Enter LAMBDA: ");
	scanf("%lf",&LAMBDA);
	
	filename="smooth_data/Info/modelpar_info.txt";
	fptr=fopen(filename.c_str(),"r");
	NPars=0;
	while(!feof(fptr)){
		fscanf(fptr,"%s\n",dummy1);
		if(!feof(fptr)){
			fscanf(fptr,"%s %lf %lf\n",dummy2,&xminread,&xmaxread);
			NPars+=1;
			xmin.push_back(xminread);
			xmax.push_back(xmaxread);
			priortype.push_back(string(dummy2));
			parname.push_back(string(dummy1));
		}
	}
	fclose(fptr);
	
	filename="figs/modelpar_info.txt";
	fptr=fopen(filename.c_str(),"w");
	for(ipar=0;ipar<NPars;ipar++){
		fprintf(fptr,"%s %s %s\n",parname[ipar].c_str(),priortype[ipar].c_str(),parname[ipar].c_str());
	}
	fclose(fptr);
	
	filename="smooth_data/Info/observable_info.txt";
	fptr=fopen(filename.c_str(),"r");
	NObs=0;
	while(!feof(fptr)){
		fscanf(fptr,"%s\n",dummy1);
		//if(!feof(fptr)){
			NObs+=1;
			obsname.push_back(string(dummy1));
			//}
		//fgets(dummy1,200,fptr);
	}
	fclose(fptr);
	CLog::Info("NPars="+to_string(NPars)+", NObs="+to_string(NObs)+"\n");
	
	filename="figs/observable_info.txt";
	fptr=fopen(filename.c_str(),"w");
	for(iy=0;iy<NObs;iy++){
		fprintf(fptr,"%s  %s\n",obsname[iy].c_str(),obsname[iy].c_str());
	}
	fclose(fptr);
		
	NTrain=0;
	bool existence;
	do{
		string filename="smooth_data/modelruns/run"+to_string(NTrain);
		filesystem::path f{filename};
		existence=filesystem::exists(f);
		if(existence){
			NTrain+=1;
		}
	}while(existence);
	CLog::Info("NTraining Pts="+to_string(NTrain)+"\n");
	
	NBandSmooth::CSmooth smooth(NPars,maxrank);
	A.resize(smooth.NCoefficients);
	
	// Now make some data for later testing, not for tuning
	
	thetatest.resize(Ntest);
	for(itest=0;itest<Ntest;itest++){
		thetatest[itest].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++){
			thetatest[itest][ipar]=-1.0+2.0*randy.ran();
		}
	}
	
	string command="mkdir -p smooth_data/fullmodel_testdata";
	system(command.c_str());
	theta.resize(NPars);
	expfilename="smooth_data/Info/experimental_info.txt";
	expfptr=fopen(expfilename.c_str(),"w");
	exptheta.resize(NPars);
	for(ipar=0;ipar<NPars;ipar++){
		exptheta[ipar]=0.2;
	}
	
	for(iy=0;iy<NObs;iy++){
		randy.reset(iy+it);
		for(ic=0;ic<smooth.NCoefficients;ic++){
			A[ic]=100.0*randy.ran_gauss();
		}
		for(itrain=0;itrain<NTrain;itrain++){
			filename="smooth_data/modelruns/run"+to_string(itrain)+"/mod_parameters.txt";
			fptr=fopen(filename.c_str(),"r");
			for(ipar=0;ipar<NPars;ipar++){
				fscanf(fptr,"%s %lf",parname_c,&x);
				parname[ipar]=parname_c;
				theta[ipar]=-1.0+2.0*(x-xmin[ipar])/(xmax[ipar]-xmin[ipar]);
			}
			fclose(fptr);
			
			y=smooth.CalcY(A,LAMBDA,theta);
			filename="smooth_data/modelruns/run"+to_string(itrain)+"/obs.txt";
			if(iy==0)
				fptr_out=fopen(filename.c_str(),"w");
			else
				fptr_out=fopen(filename.c_str(),"a");
			fprintf(fptr_out,"%s %g\n",obsname[iy].c_str(),y);
			fclose(fptr_out);
		}

		filename="smooth_data/fullmodel_testdata/"+obsname[iy]+".txt";
		fptr=fopen(filename.c_str(),"w");
		for(itest=0;itest<Ntest;itest++){
			y=smooth.CalcY(A,LAMBDA,thetatest[itest]);
			for(ipar=0;ipar<NPars;ipar++){
				fprintf(fptr,"%17.10e ",thetatest[itest][ipar]);
			}
			fprintf(fptr,"%17.10e\n",y);
		}
		fclose(fptr);
		
		y=smooth.CalcY(A,LAMBDA,exptheta);
		fprintf(expfptr,"%s  %lf  3.0 0.0\n",obsname[iy].c_str(),y);
	}
	
	fclose(expfptr);
	
	
	
	return 0;
}
