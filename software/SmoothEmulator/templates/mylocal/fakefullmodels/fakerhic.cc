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

void CalcY(vector<double> &xmin,vector<double> &xmax,vector<double> &x,vector<double> &Y){
	
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
	t[0]=theta[0]+theta[1];
	t[1]=theta[1]-0.5*theta[0];
	t[2]=theta[2]+0.5*theta[3];
	t[3]=theta[3]-theta[2];
	t[4]=theta[4]+0.75*theta[5];
	t[5]=theta[5]-0.375*theta[4];
	
	t[0]+=0.3*theta[0]*theta[1]+0.2*theta[0]*theta[0];
	t[1]-=0.2*theta[0]*theta[1];
	t[2]+=0.15*theta[2]*theta[3]+0.15*theta[2]*theta[2];
	t[3]-=0.05*theta[2]*theta[3];
	t[4]+=0.25*theta[4]*theta[5]-0.25*theta[5]*theta[5];
	t[5]-=0.1*theta[4]*theta[5];
	
	/*
	Y[0]=450+75*(t[0]+0.6*t[4]+0.5*t[1]*t[3]+0.1*t[1]*t[1]*t[3]);
	Y[1]=725+100*(t[1]-0.9*t[2]+0.4*t[6]-0.5*t[2]*t[1]+0.1*t[2]*t[3]*t[4]);
	Y[2]=1100+180*(t[2]-0.1*t[0]-0.06*t[3]*t[4]*t[2]);
	Y[3]=5.5+2.5*(t[3]+t[5])-0.1*(t[3]*t[5]+0.1*t[1]*t[2]*t[4]*t[5]);
	Y[4]=0.19+1.2*(t[4]-0.1*t[2]*t[4]);
	Y[5]=0.5-0.6*t[5]-0.5*t[3]+0.1*t[1]-0.1*t[0]*t[1]*t[2]*t[3]*t[4];*/
	Y[0]=450.0+100.0*t[0]+20*t[5];
	Y[1]=725.0+100.0*t[1]+30*t[5];
	Y[2]=1100.0+100.0*t[2]+50*t[5];
	Y[3]=5.5-2.5*t[1]+0.5*t[3]+1.5*t[5];
	Y[4]=0.2+1.25*t[5];
	Y[5]=0.5-0.5*t[5]-t[4]+0.2*t[3];
}

int main(){
	string filename;
	char dummy[200];
	unsigned int itrain,iobs,ipar,NTrain;
	vector<double> theta,xtrue,Ytrain,Ytrue,xtrain,SigmaY;
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
	
	FILE *fptr;
	xtrue.resize(NPars);
	Ytrue.resize(NObs);
	SigmaY.resize(NObs);
	
	xtrain.resize(NPars);
	Ytrain.resize(NObs);
	CLog::Info("NPars="+to_string(NPars)+"\n");
	SigmaY[0]=150.0;
	SigmaY[1]=200.0;
	SigmaY[2]=250.0;
	SigmaY[3]=1.0;
	SigmaY[4]=0.5;
	SigmaY[5]=0.5;

	fptr=fopen("Info/modelpar_info.txt","r");
	fgets(dummy,200,fptr);
	for(ipar=0;ipar<NPars;ipar++){
		fscanf(fptr,"%s %s %lf %lf",parname_c,type,&xmin[ipar],&xmax[ipar]);
		xtrue[ipar]=0.4*xmin[ipar]+0.6*xmax[ipar];
	}
	fclose(fptr);
	
	CalcY(xmin,xmax,xtrue,Ytrue);
	
	fptr=fopen("Info/experimental_info.txt","w");
	for(iobs=0;iobs<NObs;iobs++){
		SigmaY[iobs]=SigmaY[iobs]/5.0;
		fprintf(fptr,"%s\t%g\t%g 0.0\n",
		obsname[iobs].c_str(),Ytrue[iobs],SigmaY[iobs]);
	}
	fclose(fptr);
	

	for(itrain=0;itrain<NTrain;itrain++){
		filename="modelruns/run"+to_string(itrain)+"/mod_parameters.txt";
		fptr=fopen(filename.c_str(),"r");
		for(ipar=0;ipar<NPars;ipar++){
			fscanf(fptr,"%s %lf",parname_c,&xtrain[ipar]);
		}
		fclose(fptr);
		
		CalcY(xmin,xmax,xtrain,Ytrain);
		
		filename="modelruns/run"+to_string(itrain)+"/obs.txt";
		fptr=fopen(filename.c_str(),"w");
		for(iobs=0;iobs<NObs;iobs++){
			fprintf(fptr,"%s %lf %lf\n",obsname[iobs].c_str(),Ytrain[iobs],SigmaY[iobs]);
		}
		fclose(fptr);
	}

	return 0;
}
