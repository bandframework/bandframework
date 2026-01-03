#include "msu_smooth/trainingpoint_optimizer.h"
#include "msu_smooth/modelparinfo.h"
#include <cstdlib>
#include <algorithm>
#include "msu_smoothutils/randy.h"
double my_erfinv (double a);

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

void CTPO::SetThetaLatinHyperCube(vector<vector<double>> &theta){
	double thetamax,dtheta,root2=sqrt(2.0);
	unsigned int ipar,itrain,is;
	theta.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++)
		theta[itrain].resize(NPars);
	
	vector<int> ishuffle(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ishuffle[itrain]=itrain;
	}
	
	for(ipar=0;ipar<NPars;ipar++){
		std::shuffle(std::begin(ishuffle), std::end(ishuffle), randy->mt);
		thetamax=priorinfo->ThetaPrior[ipar];
		dtheta=2.0*thetamax/double(NTrainingPts);
		for(itrain=0;itrain<NTrainingPts;itrain++){
			is=ishuffle[itrain];
			if(priorinfo->type[ipar]=="uniform"){
				theta[itrain][ipar]=-thetamax+dtheta*(is+randy->ran());
			}
			else if(priorinfo->type[ipar]=="gaussian"){
				dtheta=2.0/double(NTrainingPts);
				theta[itrain][ipar]=-1.0+dtheta*(is+randy->ran());
				theta[itrain][ipar]=my_erfinv(theta[itrain][ipar]);
				theta[itrain][ipar]*=thetamax*root2;
			}
		}
	}
}

void CTPO::SetThetaSimplex(double RSimplexSet){
	unsigned int ipar,itrain,jtrain;
	double z,R;
	RSimplex=fabs(RSimplexSet);
	R=1.0;
	NTrainingPts=NPars+1;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=0.0;
	}
	
	ThetaTrain[0][0]=-R;
	ThetaTrain[1][0]=R;
	for(itrain=2;itrain<NTrainingPts;itrain++){
		z=R*itrain/sqrt(double(itrain*itrain)-1.0);
		for(jtrain=0;jtrain<itrain;jtrain++){
			ThetaTrain[jtrain][itrain-1]=-z/double(itrain);
		}
		ThetaTrain[itrain][itrain-1]=z;
		R=z;
	}
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(priorinfo->type[ipar]=="gaussian")
				ThetaTrain[itrain][ipar]*=(RSimplex/R);
			else{
				ThetaTrain[itrain][ipar]*=(RSimplex/R);
			}
		}
	}

	double Rmax=0.95;
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(priorinfo->type[ipar]=="uniform"){
				if(fabs(ThetaTrain[itrain][ipar])>Rmax){
					ThetaTrain[itrain][ipar]*=(Rmax/fabs(ThetaTrain[itrain][ipar]));
				}
			}
		}
	}
}

void CTPO::SetThetaSimplexPlus1(double RSimplexSet){
	unsigned int ipar,itrain,jtrain;
	double z,R;
	RSimplex=fabs(RSimplexSet);
	R=1.0;
	NTrainingPts=NPars+1;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=0.0;
	}
	
	ThetaTrain[0][0]=-R;
	ThetaTrain[1][0]=R;
	for(itrain=2;itrain<NTrainingPts;itrain++){
		z=R*itrain/sqrt(double(itrain*itrain)-1.0);
		for(jtrain=0;jtrain<itrain;jtrain++){
			ThetaTrain[jtrain][itrain-1]=-z/double(itrain);
		}
		ThetaTrain[itrain][itrain-1]=z;
		R=z;
	}
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(priorinfo->type[ipar]=="gaussian")
				ThetaTrain[itrain][ipar]*=(RSimplex/R);
			else{
				ThetaTrain[itrain][ipar]*=(RSimplex/R);
			}
		}
	}
	
	// Add point at origin (differentiates this from type 1
	NTrainingPts+=1;
	ThetaTrain.resize(NTrainingPts);
	ThetaTrain[NTrainingPts-1].resize(NPars);
	for(ipar=0;ipar<NPars;ipar++)
		ThetaTrain[NTrainingPts-1][ipar]=0.0;
	
	
	double Rmax=0.95;
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(priorinfo->type[ipar]=="uniform"){
				if(fabs(ThetaTrain[itrain][ipar])>Rmax){
					ThetaTrain[itrain][ipar]*=(Rmax/fabs(ThetaTrain[itrain][ipar]));
				}
			}
		}
	}
}

void CTPO::SetThetaTrain(vector<vector<double>> &theta){
	unsigned int itrain,ipar;
	for(itrain=0;itrain<ThetaTrain.size();itrain++)
		ThetaTrain[itrain].clear();
	ThetaTrain.clear();
	NTrainingPts=theta.size();
	ThetaTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=theta[itrain][ipar];
	}
}

void CTPO::WriteModelPars(){
	FILE *fptr;
	string filename,dirname,command;
	unsigned int itrain,ipar;
	vector<CModelParameters *> modelparameters(NTrainingPts);

	for(itrain=0;itrain<NTrainingPts;itrain++){
		modelparameters[itrain]=new CModelParameters();
		for(ipar=0;ipar<NPars;ipar++){
			modelparameters[itrain]->Theta[ipar]=ThetaTrain[itrain][ipar];
		}
		modelparameters[itrain]->TranslateTheta_to_X();
	}
	for(itrain=0;itrain<NTrainingPts;itrain++){
		dirname=FullModelRunsDirName+"/run"+to_string(itrain);
		command="mkdir -p "+dirname;
		system(command.c_str());
		filename=dirname+"/model_parameters.txt";
		fptr=fopen(filename.c_str(),"w");
		for(ipar=0;ipar<NPars;ipar++){
			fprintf(fptr,"%s %g\n",
			priorinfo->parname[ipar].c_str(),modelparameters[itrain]->X[ipar]);
		}
		fclose(fptr);
	}
}

void CTPO::FreezeTrainingPts(){
	string token;
	string runfilename;
	unsigned int itrain,i;
	
	TrainingPtsFreeze.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		TrainingPtsFreeze[itrain]=false;
	}
	string TPFreezeStr=parmap->getS("TPO_FreezePoints","None");
	
	if(TPFreezeStr!="None" && TPFreezeStr!="NONE" && TPFreezeStr!="none"){
		stringstream ss(TPFreezeStr);
		while(getline(ss, token, ',')) {
			size_t pos = token.find("-");
			if (pos != string::npos) {
				
				unsigned int start = stoi(token.substr(0, pos));
				unsigned int end = stoi(token.substr(pos+1));

				for (unsigned int i = start; i <= end; i++){
					TrainingPtsFreeze[i]=true;
					if(TrainingPtsRead[i]==false){
						CLog::Info("Warning: For training point "+to_string(i)+": freezing point which was set randomly\n");
					}
				}
			}
			else {
				i=stoi(token);
				TrainingPtsFreeze[i]=true;
			}
		}
	}

}

void CTPO::ReadTrainingPts(){
	string token;
	char filename[300],mod_par_name[300];
	string rundirname;
	unsigned int itrain,ipar;
	bool exists;
	FILE *fptr;
	double x;
	//string TPReadStr=parmap->getS("TrainingPointsRead","None");
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		TrainingPtsRead[itrain]=false;
	}
	string TPReadStr=parmap->getS("TPO_ReadPoints","None");
	
	if(TPReadStr!="None" && TPReadStr!="NONE" && TPReadStr!="none"){
		if(TPReadStr=="All" || TPReadStr=="ALL" || TPReadStr=="all"){
			for(itrain=0;itrain<NTrainingPts;itrain++){
				TrainingPtsRead[itrain]=true;
			}
		}
		else{
			stringstream ss(TPReadStr);
			while(getline(ss, token, ',')) {
				size_t pos = token.find("-");
				if (pos != string::npos) {
					unsigned int start = stoi(token.substr(0, pos));
					unsigned int end = stoi(token.substr(pos+1));

					for (unsigned int itrain = start; itrain <= end; itrain++){
						TrainingPtsRead[itrain]=true;
					}
				}
				else {
					itrain=stoi(token);
					TrainingPtsRead[itrain]=true;
				}
			}
		}
	}
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		if(TrainingPtsRead[itrain]){
			snprintf(filename,300,"%s/run%u/model_parameters.txt",FullModelRunsDirName.c_str(),itrain);
			if(filesystem::exists(filename)){
				exists=true;
			}
			else
				exists=false;
			if(exists)
				TrainingPtsRead[itrain]=false;
		
			if(exists){
				fptr=fopen(filename,"r");
				do{
					fscanf(fptr,"%s",mod_par_name);
					if(!feof(fptr)){
						fscanf(fptr,"%lf",&x);
						ipar=priorinfo->GetIPosition(mod_par_name);
						ThetaTrain[itrain][ipar]=x;
					}
				}while(!feof(fptr));
				fclose(fptr);
			}
		}
	}

}

void CTPO::SetTrainingPts(){
	unsigned int itrain,ipar;
	double thetamax,gausswidth;
	string TPReadStr=parmap->getS("TPO_ReadPoints","None");
	if(TPReadStr=="None" || TPReadStr=="NONE" || TPReadStr=="none"){
		SetThetaLatinHyperCube(ThetaTrain);
	}
	else{
		for(itrain=0;itrain<NTrainingPts;itrain++){
			if(TrainingPtsRead[itrain]==false){
				for(ipar=0;ipar<NPars;ipar++){
					if(priorinfo->type[ipar]=="uniform"){
						thetamax=priorinfo->ThetaPrior[ipar];
						ThetaTrain[itrain][ipar]=-thetamax+2.0*thetamax*randy->ran();
					}
					else if(priorinfo->type[ipar]=="gaussian"){
						gausswidth=priorinfo->ThetaPrior[ipar];
						ThetaTrain[itrain][ipar]=gausswidth*randy->ran_gauss();
					}
				}
			}
		}
	}
	
}

