#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cstring>
#include <filesystem>
#include "msu_smoothutils/randy.h"
#include "msu_smoothutils/log.h"

using namespace std;
using namespace NMSUUtils;
const unsigned int NObs=6,NPars=6;

void CalcY(vector<double> &xmin,vector<double> &xmax,vector<double> &x,vector<double> &Y,Crandy *randy){
	
	vector<double> xbar,theta,t;
	xbar.resize(NPars);
	theta.resize(NPars);
	t.resize(NPars);
	for(unsigned int ipar=0;ipar<NPars;ipar++){
		xbar[ipar]=0.5*(xmin[ipar]+xmax[ipar]);
		theta[ipar]=2*(x[ipar]-xbar[ipar])/(xmax[ipar]-xmin[ipar]);
	}
	/*
	t[0]=0.5*theta[0]+0.5*theta[1]+0.5*theta[2];
	t[1]=0.5*theta[0]-0.5*theta[1]+0.5*theta[3];
	t[2]=0.1*theta[0]+0.2*theta[2]+0.3*theta[3]+0.4*theta[4]+0.5*theta[5];
	t[3]=-theta[3]-0.2*theta[5];
	t[4]=theta[4]+0.4*theta[3]-0.5*theta[6];
	t[5]=-0.5*theta[5]-0.3*theta[1]+0.1*theta[3];*/
	
	t[0]=0.8*theta[0]+0.6*theta[1];
	t[1]=0.71*theta[0]-0.5*theta[1]+0.5*theta[4];
	t[2]=-0.6*theta[2]+0.8*theta[3];
	t[3]=0.71*theta[3]+0.71*theta[2];
	t[4]=0.4*theta[0]-0.4*theta[1]+0.4*theta[2]-0.4*theta[3]+0.4*theta[4]-0.4*theta[5];
	t[5]=-0.8*theta[5]-0.4*theta[5]+0.4*theta[0];
	
	randy->reset(123);

	double Lambda=2.0,rg;
	for(unsigned int I=0;I<NPars;I++){
		for(unsigned int i=0;i<NPars;i++){
			for(unsigned int j=0;j<NPars;j++){
				rg=randy->ran_gauss();
				t[I]+=0.5*rg*theta[i]*theta[j]/Lambda;
				for(unsigned int k=0;k<NPars;k++){
					rg=randy->ran_gauss();
					t[I]+=rg*theta[i]*theta[j]*theta[k]/(6.0*Lambda*Lambda);
					for(unsigned int ell=0;ell<NPars;ell++){
						rg=randy->ran_gauss();
						t[I]+=rg*theta[i]*theta[j]*theta[k]*theta[ell]/(24.0*Lambda*Lambda*Lambda);
					}
				}
			}
		}
	}
	
	Y[0]=450.0+100.0*t[0];
	Y[1]=725.0+150.0*t[1];
	Y[2]=1100.0+200.0*t[2];
	Y[3]=5.5+2.5*t[3];
	Y[4]=0.4+0.25*t[4];
	Y[5]=0.6+0.4*t[5];
}

int main(){
	string filename;
	FILE *fptr;
	char dummy[200];
	unsigned int itrain,iobs,ipar,NTrain;
	vector<double> theta,xtrue,Ytrain,Ytrue,xtrain,SigmaY,Y,X;
	vector<double> xmin(NPars),xmax(NPars);
	theta.resize(NPars);
	NMSUUtils::Crandy randy(123);
	string obsname[NObs]={"meanpt_pion","meanpt_kaon","meanpt_proton","Rinv","v2","RAA"};
	string parname[NPars]={"compressibility","etaovers","initial_flow","initial_screening","quenching_length","initial_epsilon"};
	char parname_c[200],type[100];
	
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
	xtrue.resize(NPars);
	Ytrue.resize(NObs);
	Y.resize(NObs);
	SigmaY.resize(NObs);
	
	xtrain.resize(NPars);
	Ytrain.resize(NObs);
	
	// Observable uncertainties
	SigmaY[0]=100.0;
	SigmaY[1]=150.0;
	SigmaY[2]=200.0;
	SigmaY[3]=1.5;
	SigmaY[4]=0.5;
	SigmaY[5]=0.5;

	// read in modelpar_info and set experimental value to theta=0.2
	fptr=fopen("Info/modelpar_info.txt","r");
	fgets(dummy,200,fptr);
	for(ipar=0;ipar<NPars;ipar++){
		fscanf(fptr,"%s %s %lf %lf",parname_c,type,&xmin[ipar],&xmax[ipar]);
		xtrue[ipar]=0.4*xmin[ipar]+0.6*xmax[ipar];
	}
	fclose(fptr);
	CalcY(xmin,xmax,xtrue,Ytrue,&randy);
	fptr=fopen("Info/experimental_info.txt","w");
	for(iobs=0;iobs<NObs;iobs++){
		fprintf(fptr,"%s\t%g\t%g 0.0\n",
		obsname[iobs].c_str(),Ytrue[iobs],SigmaY[iobs]/5.0);
	}
	fclose(fptr);
	
	// Write observable info for every training point
	for(itrain=0;itrain<NTrain;itrain++){
		filename="modelruns/run"+to_string(itrain)+"/mod_parameters.txt";
		fptr=fopen(filename.c_str(),"r");
		for(ipar=0;ipar<NPars;ipar++){
			fscanf(fptr,"%s %lf",parname_c,&xtrain[ipar]);
		}
		fclose(fptr);
		
		CalcY(xmin,xmax,xtrain,Ytrain,&randy);
		
		filename="modelruns/run"+to_string(itrain)+"/obs.txt";
		fptr=fopen(filename.c_str(),"w");
		for(iobs=0;iobs<NObs;iobs++){
			fprintf(fptr,"%s %lf %lf\n",obsname[iobs].c_str(),Ytrain[iobs],SigmaY[iobs]);
		}
		fclose(fptr);
	}
	
	// Write fullmodel test data for random points
	X.resize(NPars);
	string command="mkdir -p fullmodel_testdata";
	system(command.c_str());
	command="rm -f fullmodel_testdata/*.txt";
	system(command.c_str());
	unsigned int itest,Ntest=50;
	for(itest=0;itest<Ntest;itest++){
		randy.reset(itest);
		for(ipar=0;ipar<NPars;ipar++){
			theta[ipar]=-1.0+2.0*randy.ran();
			if(itest==0){
				theta[ipar]=0.2;
			}
			X[ipar]=xmin[ipar]+0.5*(theta[ipar]+1.0)*(xmax[ipar]-xmin[ipar]);
		}
		CalcY(xmin,xmax,X,Y,&randy);
		for(iobs=0;iobs<NObs;iobs++){
			filename="fullmodel_testdata/"+obsname[iobs]+".txt";
			fptr=fopen(filename.c_str(),"a");
			for(ipar=0;ipar<NPars;ipar++){
				fprintf(fptr,"%12.5e ",theta[ipar]);
			}
			fprintf(fptr,"%12.5e\n",Y[iobs]);
			fclose(fptr);
		}
	}
	return 0;
}
