#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUUtils;

void CSmoothEmulator::TuneExact(){
	unsigned int NCoefficients=smooth->NCoefficients;
	unsigned int itrain,ic,a,b;
	Psi.resize(NTrainingPts,NTrainingPts);
	Eigen::VectorXd alpha,gamma,YTrain;
	Eigen::MatrixXd C;
	alpha.resize(NCoefficients);
	gamma.resize(NCoefficients);
	beta.resize(NTrainingPts,NCoefficients);
	C.resize(NTrainingPts,NTrainingPts);
	beta.setZero();
	gamma.setZero();
	C.setZero();
	ABest.resize(NCoefficients);
	YTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		YTrain[itrain]=smoothmaster->traininginfo->YTrain[iY][itrain];
	}
	alpha=TtildeInv*YTrain;
	for(a=0;a<NTrainingPts;a++){
		for(ic=NTrainingPts;ic<NCoefficients;ic++){
			for(b=0;b<NTrainingPts;b++){
				beta(a,ic)+=TtildeInv(a,b)*T[b][ic];

			}
		}
	}
	
	GetExactQuantities();
	
	// Will solve for gamma, but first find C
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			C(a,b)=B[a][b];
			if(a==b)
				C(a,b)+=1.0;
		}
	}
	gamma=C.colPivHouseholderQr().solve(alpha);
	
	for(ic=NTrainingPts;ic<NCoefficients;ic++){
		ABest[ic]=0.0;
		for(a=0;a<NTrainingPts;a++){
			ABest[ic]+=gamma(a)*beta(a,ic);
		}
	}
	
	for(a=0;a<NTrainingPts;a++){
		ABest[a]=alpha(a);
		for(ic=NTrainingPts;ic<NCoefficients;ic++){
			ABest[a]-=beta(a,ic)*ABest[ic];
		}
	}
	
}

void CSmoothEmulator::GetExactQuantities(){
	unsigned int NCoefficients=smooth->NCoefficients;
	unsigned int i,a,b,aprime,bprime;
	
	B.resize(NTrainingPts);
	H6.resize(NTrainingPts);
	H8.resize(NTrainingPts);
	
	
	for(a=0;a<NTrainingPts;a++){
		B[a].resize(NTrainingPts);
		H6[a].resize(NTrainingPts);
		H8[a].resize(NTrainingPts);
	}
	
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			B[a][b]=0.0;
			for(i=NTrainingPts;i<NCoefficients;i++){
				B[a][b]+=beta(a,i)*beta(b,i);
			}
		}
	}
	
	// Now calculate arrays used for calculating uncertainty
	
	Eigen::MatrixXd D;
	Psi.resize(NTrainingPts,NTrainingPts);
	D.resize(NTrainingPts,NTrainingPts);
	D.setZero();
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			if(a==b)
				D(a,b)=1.0;
			D(a,b)+=B[a][b];
		}
	}

	Psi=-D.inverse();
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			H6[a][b]=0.0;
			H8[a][b]=0.0;
			for(aprime=0;aprime<NTrainingPts;aprime++){
				H6[a][b]+=B[a][aprime]*Psi(aprime,b);
				for(bprime=0;bprime<NTrainingPts;bprime++){
					H8[a][b]+=B[a][aprime]*Psi(aprime,bprime)*B[bprime][b];
				}
			}
		}
	}	
}

void CSmoothEmulator::GetExactSigmaA(){
	unsigned int ir,ic,ic0,NCoefficients=smooth->NCoefficients,MaxRank=smooth->MaxRank;
	double A2sum=0.0;
	vector<double> A2barByRank;
	vector<int> DenByRank;
	A2barByRank.resize(MaxRank+1);
	DenByRank.resize(MaxRank+1);
	for(ir=0;ir<=MaxRank;ir++){
		A2barByRank[ir]=0.0;
		DenByRank[ir]=0;
	}
	ic0=1;
	if(ConstrainA0)
		ic0=0;
	for(ic=ic0;ic<NCoefficients;ic++){
		A2sum+=ABest[ic]*ABest[ic];
		ir=smooth->rank[ic];
		A2barByRank[ir]+=ABest[ic]*ABest[ic];
		DenByRank[ir]+=1;
	}
	
	//CLog::Info("----"+observable_name+"----LAMBDA="+to_string(LAMBDA)+" ----- \n");

	if(ConstrainA0){
		SigmaA=sqrt(A2sum/double(NTrainingPts));
	}
	else{
		SigmaA=sqrt(A2sum/double(NTrainingPts-1));
	}
	//CLog::Info("SigmaA should be:"+to_string(SigmaA)+"\n");
	for(ir=0;ir<=MaxRank;ir++){
		A2barByRank[ir]=A2barByRank[ir]/double(DenByRank[ir]);
		//CLog::Info("A2barByRank[rank="+to_string(ir)+"] = "+to_string(A2barByRank[ir])+", DenRank="+to_string(DenByRank[ir])+"\n");
	}
	A2barRatio=A2barByRank[2]/A2barByRank[1];
	//CLog::Info("A2bar[2]/A2bar[1]="+to_string(A2barRatio)+"\n");

	
}

void CSmoothEmulator::CalcExactLogP(){
	double Jacobian;
	Jacobian=Ttilde.determinant();
	logP=-log(Jacobian)-NTrainingPts*log(SigmaA);
	//CLog::Info("logP="+to_string(logP)+"\n");
}


