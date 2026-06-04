#include "msu_smooth/pca.h"

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

CPCA::CPCA(){
	parmap=new CparameterMap;
	parmap->ReadParsFromFile("smooth_parameters/emulator_parameters.txt");
	string command="mkdir -p PCA_Info";
	system(command.c_str());
	
	observable_info=new CObservableInfo("Info/observable_info.txt");
	NObs=observable_info->NObservables;
	observable_info->ReadExperimentalInfo("Info/experimental_info.txt");
	
	modelruns_dirname=parmap->getS("SmoothEmulator_ModelRunDirName","modelruns");
	string NTrainingStr = parmap->getS("SmoothEmulator_TrainingPts","1");
	
	stringstream ss(NTrainingStr);
	string token;

	while(getline(ss, token, ',')) {
		size_t pos = token.find("-");
		if (pos != string::npos){

			unsigned int start = stoi(token.substr(0, pos));
			unsigned int end = stoi(token.substr(pos+1));

			for (unsigned int i = start; i <= end; i++)
				NTrainingList.push_back(i);
		}
		else {
			NTrainingList.push_back(stoi(token));
		}
	}
	nruns=NTrainingList.size();
	
	Y.resize(nruns);
	SigmaY.resize(nruns);
	
}

void CPCA::CalcTransformationInfo(){
	vector<vector<double>> Ytilde,Y;
	vector<double> Ytildebar;
	char filename[300],obs_name[300];
	double y,sigmay;
	unsigned int iy,jy,irun,nruns=NTrainingList.size();
	Eigen::MatrixXd *A;
	FILE *fptr;
	Y.resize(nruns);
	Ytilde.resize(nruns);
	for(irun=0;irun<nruns;irun++){
		Y[irun].resize(NObs);
		Ytilde[irun].resize(NObs);
	}
	Ybar.resize(NObs);
	SigmaY.resize(NObs);
	for(iy=0;iy<NObs;iy++){
		Ybar[iy]=0.0;
		SigmaY[iy]=observable_info->SigmaExp[iy];
	}
	
	for(irun=0;irun<nruns;irun++){
		snprintf(filename,300,"smooth_data/%s/run%u/obs.txt",modelruns_dirname.c_str(),NTrainingList[irun]);
		fptr=fopen(filename,"r");
		while(!feof(fptr)){
			fscanf(fptr,"%s",obs_name);
			if(!feof(fptr)){
				fscanf(fptr,"%lf %lf",&y, &sigmay);
				iy=observable_info->GetIPosition(obs_name);
				Y[irun][iy]=y;
				Ybar[iy]+=y/double(nruns);
			}
		}
		fclose(fptr);
	}
	
	for(iy=0;iy<NObs;iy++){
		for(irun=0;irun<nruns;irun++){
			Ytilde[irun][iy]=(Y[irun][iy]-Ybar[iy])/SigmaY[iy];
		}
	}
	
	// Calculate Covariance Matrix
	
	A = new Eigen::MatrixXd(NObs,NObs);
	A->setZero();
	for(iy=0;iy<NObs;iy++){
		for(jy=0;jy<NObs;jy++){
			for(irun=0;irun<nruns;irun++){
				(*A)(iy,jy)+=(Ytilde[irun][iy]*Ytilde[irun][jy])/double(nruns);
			}
		}
	}
				
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(*A);
	eigvals = es.eigenvalues();
	eigvecs = es.eigenvectors();
	delete A;
	
	// Write Transformation Info
	ofstream PCAInfoFile("PCA_Info/transformation_info.txt");
	
	if (PCAInfoFile.is_open()){
		PCAInfoFile << "NObs  Nruns" << endl;
		PCAInfoFile << NObs <<  " " << nruns << endl;
		PCAInfoFile << "<Y>" << endl;
		for(iy=0;iy<NObs;iy++)
			PCAInfoFile << Ybar[iy] << " ";
		PCAInfoFile << endl;
		PCAInfoFile << "SigmaY" << endl;
		for(iy=0;iy<NObs;iy++)
			PCAInfoFile << SigmaY[iy] << " ";
		PCAInfoFile << endl;
		PCAInfoFile << "EigenValues" << endl;
		for(iy=0;iy<NObs;iy++)
			PCAInfoFile << eigvals(iy) << " ";
		PCAInfoFile << endl;
		PCAInfoFile << "EigenVectors" << endl;
		for(iy=0;iy<NObs;iy++){
			for(jy=0;jy<NObs;jy++){
				PCAInfoFile << eigvecs(iy,jy) << " ";
			}
			PCAInfoFile << endl;
		}
	}
	PCAInfoFile.close();

	// Write PCA Observables for Training Pts
	
	vector<vector<double>> Z;
	unsigned int iz;
	string pcaname;
	
	Z.resize(nruns);
	for(irun=0;irun<nruns;irun++){
		Z[irun].resize(nruns);
	}
	
	for(irun=0;irun<nruns;irun++){
		for(iz=0;iz<NObs;iz++){
			Z[irun][iz]=0.0;
			for(iy=0;iy<NObs;iy++){
				Z[irun][iz]+=eigvecs(iz,iy)*Ytilde[irun][iy];
			}
		}
		snprintf(filename,300,"smooth_data/%s/run%u/obs_pca.txt",modelruns_dirname.c_str(),NTrainingList[irun]);
		fptr=fopen(filename,"w");
		for(iz=0;iz<NObs;iz++){
			pcaname="z"+to_string(iz);
			fprintf(fptr,"%s  %g  0.0\n",pcaname.c_str(),Z[irun][iz]);	
		}
		fclose(fptr);
	}
	
	// Write Info File for PCA
	
	fptr=fopen("PCA_Info/observable_info.txt","w");
	for(iz=0;iz<NObs;iz++){
		pcaname="z"+to_string(iz);
		fprintf(fptr,"%s\n",pcaname.c_str());
	}
	fclose(fptr);
	
	// write Experimental File for PCA
	vector<double> YExpTilde,ZExp;
	YExpTilde.resize(NObs);
	ZExp.resize(NObs);
	fptr=fopen("PCA_Info/experimental_info.txt","w");
	for(iy=0;iy<NObs;iy++){
		YExpTilde[iy]=(observable_info->YExp[iy]-Ybar[iy])/SigmaY[iy];
	}
	for(iz=0;iz<NObs;iz++){
		ZExp[iz]=0.0;
		for(iy=0;iy<NObs;iy++){
			ZExp[iz]+=eigvecs(iz,iy)*YExpTilde[iy];
		}
		pcaname="z"+to_string(iz);
		fprintf(fptr,"%s %g %g\n",pcaname.c_str(),ZExp[iz],1.0);
	}
	
	fclose(fptr);
	
}

void CPCA::ReadTransformationInfo(){
	char dummy[2000];
	unsigned int iy,jy;
	FILE *fptr=fopen("PCA_Info/transformation_info.txt","r");
	fgets(dummy,2000,fptr);
	fscanf(fptr,"%u %u",&NObs,&nruns);
	Ybar.resize(NObs);
	SigmaY.resize(NObs);
	eigvals.resize(NObs);
	eigvecs.resize(NObs,NObs);
	fscanf(fptr,"%s",dummy);
	for(iy=0;iy<NObs;iy++){
		fscanf(fptr,"%lf ",&Ybar[iy]);
	}
	fscanf(fptr,"%s",dummy);
	for(iy=0;iy<NObs;iy++){
		fscanf(fptr,"%lf ",&SigmaY[iy]);
	}
	fscanf(fptr,"%s",dummy);
	for(iy=0;iy<NObs;iy++){
		fscanf(fptr,"%lf",&eigvals(iy));
	}
	fscanf(fptr,"%s",dummy);
	for(iy=0;iy<NObs;iy++){
		for(jy=0;jy<NObs;jy++){
			fscanf(fptr,"%lf ",&eigvecs(iy,jy));
		}
	}
	fclose(fptr);
}

void CPCA::TransformZtoY(vector<double> &Z,vector<double> &SigmaZ_emulator,
vector<double> &Y,vector<double> &SigmaY_emulator){
	unsigned int iy,ky;
	SigmaY_emulator.resize(NObs);
	for(iy=0;iy<NObs;iy++){
		SigmaY_emulator.resize(NObs);
	}
	for(iy=0;iy<NObs;iy++){
		Y[iy]=0.0;
		for(ky=0;ky<NObs;ky++){
			Y[iy]+=Z[ky]*eigvecs(ky,iy);
		}
		SigmaY_emulator[iy]=0.0;
		for(ky=0;ky<NObs;ky++){
			SigmaY_emulator[iy]+=SigmaZ_emulator[ky]*SigmaZ_emulator[ky]*eigvecs(ky,iy)*eigvecs(ky,iy);
		}
		SigmaY_emulator[iy]=SigmaY[iy]*sqrt(SigmaY_emulator[iy]);
	}
}


void CPCA::TransformYtoZ(vector<double> &Z,vector<double> &SigmaZ_emulator,
vector<double> &Y,vector<double> &SigmaY_emulator){
	unsigned int iy,ky;
	SigmaZ_emulator.resize(NObs);
	for(iy=0;iy<NObs;iy++){
		SigmaZ_emulator.resize(NObs);
	}
	for(iy=0;iy<NObs;iy++){
		Z[iy]=0.0;
		for(ky=0;ky<NObs;ky++){
			Z[iy]+=Y[ky]*eigvecs(iy,ky);
		}
		SigmaZ_emulator[iy]=0.0;
		for(ky=0;ky<NObs;ky++){
			SigmaZ_emulator[iy]+=(SigmaY_emulator[ky]*SigmaY_emulator[ky]/(SigmaY[iy]*SigmaY[iy]))*eigvecs(iy,ky)*eigvecs(iy,ky);
		}
		SigmaZ_emulator[iy]=sqrt(SigmaZ_emulator[iy]);
	}
}


