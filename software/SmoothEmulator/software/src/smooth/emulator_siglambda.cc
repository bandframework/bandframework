#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUUtils;

void CSmoothEmulator::CalcSigmaA(){
	CalcB();
	int a,b;
	double sigmaA2=0.0;
	for(a=0;a<int(NTrainingPts);a++){
		for(b=0;b<int(NTrainingPts);b++){
			sigmaA2+=smoothmaster->traininginfo->YTrain[iY][a]
				*Binv(a,b)*smoothmaster->traininginfo->YTrain[iY][b];
		}
	}
	sigmaA2=sigmaA2/double(NTrainingPts);
	SigmaA=sqrt(fabs(sigmaA2));

}

void CSmoothEmulator::CalcLogP(){
	Eigen::MatrixXd BB=B;
	double exparg=0.0;
	double detfactor=4*sqrt(NTrainingPts);
	for(int a=0;a<int(NTrainingPts);a++){
		for(int b=0;b<int(NTrainingPts);b++){
			exparg-=0.5*smoothmaster->traininginfo->YTrain[iY][a]*Binv(a,b)*smoothmaster->traininginfo->YTrain[iY][b];
		}
   }
   exparg=exparg/(SigmaA*SigmaA);
   BB=B*detfactor;
   detBB=BB.determinant();
   logP=-0.5*log(fabs(detBB))-NTrainingPts*log(SigmaA)+exparg;
   if(isinf(logP)){
      CLog::Fatal("XlogP in CalcSigmaA is infinite: detBB="+to_string(detBB)+", Lambda="+to_string(LAMBDA)
                  +", SigmaA="+to_string(SigmaA)+", exparg="+to_string(exparg)+"\n");
   }
}

void CSmoothEmulator::CalcLambdaVariance(){
   double L,Lbar,dL=0.1,L2bar,norm,w;
   norm=Lbar=L2bar=0.0;
   for(L=1.0;L<8;L+=dL){
		LAMBDA=L;
		CalcB();
		CalcSigmaA();
		CalcLogP();
		//w=exp(logP);
		w=sqrt(fabs(detBB));
		double factor=pow(SigmaA/100.0,NTrainingPts);
		w=w/factor;
		norm+=w;
		Lbar+=w*L;
		L2bar+=w*L*L;
	}
	Lbar=Lbar/norm;
	L2bar=L2bar/norm;
	LambdaVariance=L2bar-Lbar*Lbar;
	LAMBDA=Lbar;
}

void CSmoothEmulator::CalcSigmaALambda(){
	double LambdaMin=sqrt(double(NPars)/3.0),LambdaMax=8;
	vector<double> bestLAMBDA(3),bestlogP(3),bestSigmaA(3);
	double dLambda=1.0;
	Eigen::Matrix3d A123;
	Eigen::Vector3d X123,Y123;
	int ntry=0,nsuccess=0,il;
	if(LambdaMin<1.0)
		LambdaMin=1.0;

	LAMBDA=LambdaMin;
	CalcB();
	if(INCLUDE_LAMBDA_UNCERTAINTY){
		CalcWBprimeChi();
	}
	CalcSigmaA();
	CalcLogP();
	for(il=0;il<3;il++){
		bestLAMBDA[il]=-1000-il;
		bestSigmaA[il]=SigmaA;
		bestlogP[il]=-100000.0-il;
	}

	LAMBDA=LambdaMin;
	while(nsuccess<3){
		ntry+=1;
		CalcB();
		if(INCLUDE_LAMBDA_UNCERTAINTY){
			CalcWBprimeChi();
		}
		CalcSigmaA();
		CalcLogP();
		if(!isfinite(detBB)){
			cout << "|B|=" << detBB << "B=\n" << B << endl;
			cout << "Binv=\n" << Binv << endl;
			exit(1);
		}
		if(logP>bestlogP[0]){
			bestLAMBDA[2]=bestLAMBDA[1];
			bestSigmaA[2]=bestSigmaA[1];
			bestlogP[2]=bestlogP[1];
			
			bestLAMBDA[1]=bestLAMBDA[0];
			bestSigmaA[1]=bestSigmaA[0];
			bestlogP[1]=bestlogP[0];
			
			bestLAMBDA[0]=LAMBDA;
			bestSigmaA[0]=SigmaA;
			bestlogP[0]=logP;
		}
		else if(logP>bestlogP[1]){
			bestLAMBDA[2]=bestLAMBDA[1];
			bestSigmaA[2]=bestSigmaA[1];
			bestlogP[2]=bestlogP[1];
			
			bestLAMBDA[1]=LAMBDA;
			bestSigmaA[1]=SigmaA;
			bestlogP[1]=logP;
		}
		else if(logP>bestlogP[2]){
			bestLAMBDA[2]=LAMBDA;
			bestSigmaA[2]=SigmaA;
			bestlogP[2]=logP;
		}

		if(ntry<3){
			LAMBDA=LAMBDA+dLambda;
		}
		else{
			if(bestLAMBDA[0]>bestLAMBDA[1] && bestLAMBDA[0]>bestLAMBDA[2] && nsuccess==0){
				LAMBDA=bestLAMBDA[0]+dLambda;
			}
			else if(bestLAMBDA[0]<bestLAMBDA[1] && bestLAMBDA[0]<bestLAMBDA[2] && nsuccess==0){
				LAMBDA=bestLAMBDA[0]-dLambda;
			}
			else{
				for(il=0;il<3;il++){
					A123(il,0)=bestLAMBDA[il]*bestLAMBDA[il];
					A123(il,1)=bestLAMBDA[il];
					A123(il,2)=1.0;
					Y123(il)=bestlogP[il];
				}
				X123=A123.colPivHouseholderQr().solve(Y123);
				LAMBDA=-0.5*X123(1)/X123(0);
				nsuccess+=1;
				dLambda*=0.5;
			}
		}
		if(LAMBDA<LambdaMin){
			LAMBDA=LambdaMin;
			nsuccess=100000;
		}
		if(LAMBDA>LambdaMax){
			LAMBDA=LambdaMax;
			nsuccess=100000;
		}
	}

	CalcB();
	CalcSigmaA();
	CalcLogP();
	if(INCLUDE_LAMBDA_UNCERTAINTY){
		CalcWBprimeChi();
		d2lndetBBdLambda2=pow(LAMBDA,-6)*((Binv*Bprimeprime)-(Binv*Bprime*Binv*Bprime)).trace()
			-3.0*pow(LAMBDA,-4)*(Binv*Bprime).trace();
	}


}

void CSmoothEmulator::CalcWBprimeChi(){
	unsigned int a,b,ipar;
	double dtheta,dtheta2,C0,sa2=SigmaA*SigmaA,YBinvY;
	Eigen::VectorXd ytrain(NTrainingPts);
	Eigen::MatrixXd BB(NTrainingPts,NTrainingPts);
	
	Bprime.resize(NTrainingPts,NTrainingPts);
	Bprimeprime.resize(NTrainingPts,NTrainingPts);
	chiprime.resize(NTrainingPts);
	for(a=0;a<NTrainingPts;a++){
		ytrain(a)=smoothmaster->traininginfo->YTrain[iY][a];
		for(b=0;b<NTrainingPts;b++){
			dtheta2=0.0;
			for(ipar=0;ipar<NPars;ipar++){
				dtheta=ThetaTrain[a][ipar]-ThetaTrain[b][ipar];
				dtheta2+=dtheta*dtheta;
			}
			C0=GetCorrelation(ThetaTrain[a],ThetaTrain[b]);
			Bprime(a,b)=C0*dtheta2;
			Bprimeprime(a,b)=Bprime(a,b)*dtheta2;
		}
	}
	YBinvY=ytrain.transpose()*Binv*ytrain;
	Winv(0,0)=-double(NTrainingPts)/sa2+(3.0/(sa2*sa2))*YBinvY;
	
	BB=Binv*Bprime*Binv;
	Winv(1,0)=(-1.0/pow(SigmaA*LAMBDA,3))*ytrain.transpose()*BB*ytrain;
	Winv(0,1)=Winv(1,0);
	
	double term1,term2,term3;
	term1=2.0*ytrain.transpose()*BB*Bprime*Binv*ytrain;
	term2=ytrain.transpose()*Binv*Bprimeprime*Binv*ytrain;
	term3=ytrain.transpose()*BB*ytrain;
	Winv(1,1)=-0.5*d2lndetBBdLambda2+(0.5/(sa2*pow(LAMBDA,6)))*(term1-term2)
		+(3.0/(2.0*sa2*pow(LAMBDA,4)))*term3;
	W=Winv.inverse();

	chiprime=Binv*Bprime*Binv*ytrain;
	
}

double CSmoothEmulator::GetSigma2_Lambda(vector<double> &theta){
	unsigned int a,ipar;
	double dtheta2,answer,dEdLambda;
	Eigen::VectorXd c0(NTrainingPts);
	Eigen::VectorXd cprime0(NTrainingPts);
	for(a=0;a<NTrainingPts;a++){
		dtheta2=0.0;
		for(ipar=0;ipar<NPars;ipar++){
			dtheta2+=pow(theta[ipar]-ThetaTrain[a][ipar],2);
		}
		c0(a)=GetCorrelation(theta,ThetaTrain[a]);
		cprime0(a)=dtheta2*c0(a);
	}
	double ABC=(cprime0.transpose()*chi);
	ABC=ABC-(c0.transpose()*chiprime);
	dEdLambda=ABC/(LAMBDA*LAMBDA*LAMBDA);
	answer=dEdLambda*dEdLambda*W(1,1);
	return answer;
}

	
	
