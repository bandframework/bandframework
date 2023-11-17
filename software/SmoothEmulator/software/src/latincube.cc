#include "msu_smooth/simplex.h"
using namespace std;
using namespace NBandSmooth;
using namespace NMSUPratt;

void NAlternativeParameterSampling::GetParsLHC(int NRuns,int NPars,Crandy *randy,vector<vector<double>> &Theta){
	int irun,ipar;
	vector<double> ThetaOrdered;
	ThetaOrdered.resize(NRuns);
	vector<int> ishuffle;
	ishuffle.resize(NRuns);
	double dTheta;
	std::random_device rd;
	std::mt19937 rangen(rd());
	
	if(int(Theta.size())!=NRuns){
		Theta.resize(NRuns);
		for(irun=0;irun<NRuns;irun++){
			Theta[irun].resize(NPars);
		}
	}
	else{
		for(irun=0;irun<NRuns;irun++){
			if(int(Theta[irun].size())!=NPars){
				Theta[irun].resize(NPars);
			}
		}
	}
	
	FILE *fptr=fopen("figs/hypercubedata.txt","w");
	dTheta=2.0/NRuns;
	for(ipar=0;ipar<NPars;ipar++){
		for(irun=0;irun<NRuns;irun++){
			//ThetaOrdered[irun]=-1.0+dTheta*(irun+randy->ran());
			ThetaOrdered[irun]=-1.0+dTheta*(irun+0.5);
			ishuffle[irun]=irun;
		}
		std::shuffle(ishuffle.begin(),ishuffle.end(),rangen);
		for(irun=0;irun<NRuns;irun++){
			Theta[irun][ipar]=ThetaOrdered[ishuffle[irun]];
			fprintf(fptr,"%10.6f ",Theta[irun][ipar]);
		}
		fprintf(fptr,"\n");
		
		
		
	}
	fclose(fptr);
	
}

void NAlternativeParameterSampling::GetParsLHC_Modified(int NRuns,int NPars,Crandy *randy,vector<vector<double>> &Theta){
	int irun,ipar;
	vector<double> ThetaOrdered;
	vector<vector<double>> BestTheta;
	ThetaOrdered.resize(NRuns);
	vector<int> ishuffle;
	double PE,BestPE=1.0E99;
	ishuffle.resize(NRuns);
	double dTheta;
	std::random_device rd;
	std::mt19937_64 rangen(rd());
	int Ntry=1000,itry ;
	
	Theta.resize(NRuns);
	BestTheta.resize(NRuns);
	for(irun=0;irun<NRuns;irun++){
		Theta[irun].resize(NPars);
		BestTheta[irun].resize(NPars);
	}
	ThetaOrdered.resize(NRuns);
	ishuffle.resize(NRuns);
	
	
	dTheta=2.0/NRuns;
	for(ipar=0;ipar<NPars;ipar++){
		for(irun=0;irun<NRuns;irun++){
			//ThetaOrdered[irun]=-1.0+dTheta*(irun+randy->ran());
			ThetaOrdered[irun]=-1.0+dTheta*(irun+0.5);
			ishuffle[irun]=irun;
		}
		for(irun=0;irun<NRuns;irun++){
			BestTheta[irun][ipar]=ThetaOrdered[irun];
		}
	}
	BestPE=GetPEShuffle(BestTheta);
	
	printf("BestPE=%g\n",BestPE);
	
	for(itry=0;itry<Ntry;itry++){
		for(ipar=0;ipar<NPars;ipar++){
			std::shuffle(ishuffle.begin(),ishuffle.end(),randy->mt);
			for(irun=0;irun<NRuns;irun++){
				Theta[irun][ipar]=ThetaOrdered[ishuffle[irun]];
			}
		}
		PE=GetPEShuffle(Theta);
		if(PE<BestPE){
			BestPE=PE;
			BestTheta=Theta;
			printf("BestPE=%g\n",BestPE);
		}
	}
		
	
	FILE *fptr=fopen("figs/hypercubedata.txt","w");
	for(ipar=0;ipar<NPars;ipar++){
		for(irun=0;irun<NRuns;irun++){
			Theta[irun][ipar]=BestTheta[ishuffle[irun]][ipar];
			fprintf(fptr,"%10.6f ",Theta[irun][ipar]);
		}
		fprintf(fptr,"\n");
	}
	fclose(fptr);
	
}

double NAlternativeParameterSampling::GetPEShuffle(vector<vector<double>> x){
	double PE=0.0,r2;
	int idim,ipart,jpart;
	int NParts=x.size();
	int NDim=x[0].size();
	for(ipart=0;ipart<NParts-1;ipart++){
		for(jpart=ipart+1;jpart<NParts;jpart++){
			r2=0.0;
			for(idim=0.0;idim<NDim;idim++){
				r2+=pow(x[ipart][idim]-x[jpart][idim],2);
			}
			PE+=1.0/sqrt(r2);
		}
	}
	return PE;
}

void NAlternativeParameterSampling::GetParsCoulomb(int NParts,int NDim,Crandy *randy,vector<vector<double>> &x){
	int ipart,idim;
	vector<vector<double>> xx,v,vv;
	double t,dt=0.001,PE,Etot,KE;
	if(x.size()!=NParts){
		x.resize(NParts);
	}
	xx.resize(NParts);
	v.resize(NParts);
	vv.resize(NParts);
	for(ipart=0;ipart<NParts;ipart++){
		x[ipart].resize(NDim);
		xx[ipart].resize(NDim);
		v[ipart].resize(NDim);
		vv[ipart].resize(NDim);
		for(idim=0;idim<NDim;idim++){
			xx[ipart][idim]=x[ipart][idim];
			v[ipart][idim]=vv[ipart][idim]=0.0;
		}
	}
	NAlternativeParameterSampling::CalcEnergy(x,vv,v,PE,KE,Etot);
	printf("Etot=%g\n",Etot);
	printf("ready to propagate\n");
	for(idim=0;idim<NDim;idim++){
		for(ipart=0;ipart<NParts;ipart++)
			printf("%8.5f ",x[ipart][idim]);
		printf("\n");
	}
	for(int itest=0;itest<10;itest++){
		for(t=0.0;t<10.0;t+=dt){
			NAlternativeParameterSampling::Propagate(dt,x,xx,v,vv,PE);
			NAlternativeParameterSampling::Propagate(dt,xx,x,vv,v,PE);
		}
		NAlternativeParameterSampling::CalcEnergy(x,vv,v,PE,KE,Etot);
		printf("Etot=%g\n",Etot);
		
		/*
		for(ipart=0;ipart<NParts;ipart++){
			for(idim=0;idim<NDim;idim++){
				v[ipart][idim]*=0.75;
				vv[ipart][idim]*=0.75;
			}
		}
		*/
		NAlternativeParameterSampling::CalcEnergy(x,vv,v,PE,KE,Etot);
		printf("After Cooling Etot=%g\n",Etot);
	}

}

void NAlternativeParameterSampling::CalcEnergy(vector<vector<double>> &x,vector<vector<double>> &vv,vector<vector<double>> &v,double &PE,double &KE,double &Etot){
	int ipart,idim;
	int NParts=x.size(), NDim=x[0].size();
	double vi;
	KE=0.0;
	for(ipart=0;ipart<NParts;ipart++){
		for(idim=0;idim<NDim;idim++){
			vi=0.5*(v[ipart][idim]+vv[ipart][idim]);
			KE+=0.5*vi*vi;
		}
	}
	Etot=PE+KE;
	printf("PE=%g, KE=%g\n",PE,KE);
}

void NAlternativeParameterSampling::GetForcePotential(vector<double> &x,vector<double> &xx,vector<double> &Frel,double &potential){
	int idim,NDim=x.size();
	double r2=0.0,r,delx,alpha=1.0;
	potential=0.0;
	for(idim=0;idim<NDim;idim++){
		delx=x[idim]-xx[idim];
		if(delx>1.0)
			delx-=2.0;
		if(delx<-1.0)
			delx+=2.0;
		r2+=pow(delx,2);
	}
	if(r2<0.0001){
		printf("r2=%g\n",r2);
	}
	if(r2<1.0){
		r=sqrt(r2);
		potential=alpha*(1.0/r+r-2.0);
		for(idim=0.0;idim<NDim;idim++){
			delx=x[idim]-xx[idim];
			if(delx>1.0)
				delx-=2.0;
			if(delx<-1.0)
				delx+=2.0;
			Frel[idim]=alpha*(delx/pow(r,3))-alpha*delx;
		}
	}
	else{
		potential=0.0;
		for(idim=0.0;idim<NDim;idim++){
			Frel[idim]=0.0;
		}
	}
}

void NAlternativeParameterSampling::Propagate(double dt,vector<vector<double>> &x,vector<vector<double>> &xx,vector<vector<double>> &v,vector<vector<double>> &vv,double &PE){
	int idim,ipart,jpart,NParts=x.size(),NDim=x[0].size();
	double potential;
	PE=0.0;
	vector<double> Frel;
	vector<vector<double>> F;
	
	for(ipart=0;ipart<NParts;ipart++){
		for(idim=0;idim<NDim;idim++){
			xx[ipart][idim]=x[ipart][idim]+v[ipart][idim]*dt;
			while(xx[ipart][idim]>1.0)
				xx[ipart][idim]-=2.0;
			while(xx[ipart][idim]<-1.0)
				xx[ipart][idim]+=2.0;
		}
	}
	
	Frel.resize(NDim);
	for(idim=0;idim<NDim;idim++)
		Frel[idim]=0.0;
	F.resize(NParts);
	for(ipart=0;ipart<NParts;ipart++){
		F[ipart].resize(NDim);
		for(idim=0;idim<NDim;idim++){
			F[ipart][idim]=0.0;
		}
	}
	
	
	for(ipart=0;ipart<NParts-1;ipart++){
		for(jpart=ipart+1;jpart<NParts;jpart++){
			GetForcePotential(xx[ipart],xx[jpart],Frel,potential);
			PE+=potential;
			for(idim=0;idim<NDim;idim++){
				F[ipart][idim]+=Frel[idim];
				F[jpart][idim]-=Frel[idim];
			}
		}
	}
	for(ipart=0;ipart<NParts;ipart++){
		for(idim=0;idim<NDim;idim++){
			vv[ipart][idim]=v[ipart][idim]+dt*F[ipart][idim];
		}
	}
}

//-------------------------------------------------

void NAlternativeParameterSampling::GetParsCoulombHO(Crandy *randy,vector<vector<double>> &Theta){
	printf("----- entering GetParsCoulombHO\n");
	int ipart,idim;
	int NParts=Theta.size();
	int NDim=Theta[0].size();
	vector<vector<double>> x,xx,v,vv;
	double t,dt=0.001,PE,Etot,KE;
	
	x.resize(NParts);
	xx.resize(NParts);
	v.resize(NParts);
	vv.resize(NParts);
	for(ipart=0;ipart<NParts;ipart++){
		x[ipart].resize(NDim);
		xx[ipart].resize(NDim);
		v[ipart].resize(NDim);
		vv[ipart].resize(NDim);
		for(idim=0;idim<NDim;idim++){
			x[ipart][idim]=atanh(Theta[ipart][idim]);
			xx[ipart][idim]=x[ipart][idim];
			v[ipart][idim]=vv[ipart][idim]=0.0;
		}
	}
	NAlternativeParameterSampling::CalcEnergyHO(x,vv,v,PE,KE,Etot);
	printf("ready to propagate\n");
	
	for(int itest=0;itest<20;itest++){
		NAlternativeParameterSampling::CalcEnergyHO(x,vv,v,PE,KE,Etot);
		printf("Before: Etot=%g, KE=%g, PE=%g\n",Etot,KE,PE);
		for(t=0.0;t<100.0;t+=dt){
			NAlternativeParameterSampling::PropagateHO(dt,x,xx,v,vv,PE);
			NAlternativeParameterSampling::PropagateHO(dt,xx,x,vv,v,PE);
		}
		NAlternativeParameterSampling::CalcEnergyHO(x,vv,v,PE,KE,Etot);
		printf("After: Etot=%g, KE=%g, PE=%g\n",Etot,KE,PE);
		for(idim=0;idim<NDim;idim++){
			for(ipart=0;ipart<NParts;ipart++){
				printf("%8.5f ",tanh(x[ipart][idim]));
				v[ipart][idim]*=0.75;
				vv[ipart][idim]*=0.75;
			}
			printf("\n");
		}
		
		/*
		for(ipart=0;ipart<NParts;ipart++){
			for(idim=0;idim<NDim;idim++){
				v[ipart][idim]*=0.75;
				vv[ipart][idim]*=0.75;
			}
		}
		NAlternativeParameterSampling::CalcEnergy(x,vv,v,PE,KE,Etot);
		*/
	}
	
	FILE *fptr=fopen("figs/coulomb.txt","w");
	for(idim=0;idim<NDim;idim++){
		for(ipart=0;ipart<NParts;ipart++){
			Theta[ipart][idim]=tanh(x[ipart][idim]);
			fprintf(fptr,"%10.6f ",Theta[ipart][idim]);
		}
		fprintf(fptr,"\n");
	}
	fclose(fptr);

}

void NAlternativeParameterSampling::GetFRelVRelHO(vector<double> &x,vector<double> &xx,vector<double> &Frel,double &Vrel){
	int idim,NDim=x.size();
	double r2=0.0,r,delx,alpha=0.5;
	Vrel=0.0;
	for(idim=0;idim<NDim;idim++){
		delx=x[idim]-xx[idim];
		r2+=pow(delx,2);
	}
	if(r2<0.0001){
		printf("r2=(");
		for(idim=0;idim<NDim;idim++){
			printf("%g,",x[idim]-xx[idim]);
		}
		printf(")\n");
	}
	r=sqrt(r2);
	Vrel=alpha*(1.0/r);
	for(idim=0.0;idim<NDim;idim++){
		delx=x[idim]-xx[idim];
		Frel[idim]=alpha*(delx/pow(r,3));
	}
}

void NAlternativeParameterSampling::PropagateHO(double dt,vector<vector<double>> &x,vector<vector<double>> &xx,vector<vector<double>> &v,vector<vector<double>> &vv,double &PE){
	int idim,ipart,jpart,NParts=x.size(),NDim=x[0].size();
	double Vrel,Vext;
	PE=0.0;
	vector<double> Frel,Fext;
	vector<vector<double>> F;
	
	for(ipart=0;ipart<NParts;ipart++){
		for(idim=0;idim<NDim;idim++){
			xx[ipart][idim]=x[ipart][idim]+v[ipart][idim]*dt;
		}
	}
	
	Frel.resize(NDim);
	for(idim=0;idim<NDim;idim++)
		Frel[idim]=0.0;
	F.resize(NParts);
	for(ipart=0;ipart<NParts;ipart++){
		F[ipart].resize(NDim);
		for(idim=0;idim<NDim;idim++){
			F[ipart][idim]=0.0;
		}
	}
	Fext.resize(NDim);
	
	for(ipart=0;ipart<NParts;ipart++){
		GetFextVextHO(xx[ipart],Fext,Vext);
		PE+=Vext;
		for(idim=0;idim<NDim;idim++)
			F[ipart][idim]+=Fext[idim];
	}
	
	for(ipart=0;ipart<NParts-1;ipart++){
		for(jpart=ipart+1;jpart<NParts;jpart++){
			GetFRelVRelHO(xx[ipart],xx[jpart],Frel,Vrel);
			PE+=Vrel;
			for(idim=0;idim<NDim;idim++){
				F[ipart][idim]+=Frel[idim];
				F[jpart][idim]-=Frel[idim];
			}
		}
	}
	for(ipart=0;ipart<NParts;ipart++){
		for(idim=0;idim<NDim;idim++){
			vv[ipart][idim]=v[ipart][idim]+dt*F[ipart][idim];
		}
	}
}

void NAlternativeParameterSampling::CalcEnergyHO(vector<vector<double>> &x,vector<vector<double>> &vv,vector<vector<double>> &v,double &PE,double &KE,double &Etot){
	int ipart,idim;
	int NParts=x.size(), NDim=x[0].size();
	double vi;
	KE=0.0;
	for(ipart=0;ipart<NParts;ipart++){
		for(idim=0;idim<NDim;idim++){
			vi=0.5*(v[ipart][idim]+vv[ipart][idim]);
			KE+=0.5*vi*vi;
		}
	}
	Etot=PE+KE;
}

void NAlternativeParameterSampling::GetFextVextHO(vector<double> &x,vector<double> &Fext,double &Vext){
	int idim,NDim=x.size();
	double k=0.5,r2=0.0;
	for(idim=0;idim<NDim;idim++){
		r2+=x[idim]*x[idim];
		Fext[idim]=-k*x[idim];
	}
	Vext=0.5*k*r2;
}