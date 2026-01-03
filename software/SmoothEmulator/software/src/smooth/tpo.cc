#include "msu_smooth/trainingpoint_optimizer.h"
#include "msu_smooth/modelparinfo.h"
#include <cstdlib>
#include <algorithm>
#include "msu_smoothutils/randy.h"

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

CTPO::CTPO(){
	FIRSTCALL=true;
   parmap=new CparameterMap();
	parmap->ReadParsFromFile("smooth_data/Options/tpo_options.txt");
   randy=new Crandy(parmap->getD("TPO_RANSEED",12345));
	string logfilename=parmap->getS("TPO_LogFileName","Screen");
	if(logfilename!="Screen"){
		CLog::Init(logfilename);
	}
	TPO_Method=parmap->getS("TPO_Method","MC");
   FullModelRunsDirName=parmap->getS("TPO_FullModelRunsDirName","smooth_data/FullModelRuns");
	string prior_info_filename="smooth_data/Info/prior_info.txt";
	priorinfo=new CPriorInfo(prior_info_filename);
	CModelParameters::priorinfo=priorinfo;
	NPars=priorinfo->NModelPars;
	INCLUDE_LAMBDA_UNCERTAINTY=parmap->getB("TPO_Include_LAMBDA_Uncertainty",true);
	CreateTrainingPts();
   LAMBDA=parmap->getD("TPO_LAMBDA",2.5);
}

void CTPO::CreateTrainingPts(){
	unsigned int itrain;
	if(TPO_Method=="MC" || TPO_Method=="MCSphere"){
		NTrainingPts=parmap->getI("TPO_NTrainingPts",0);
		if(NTrainingPts==0){
			CLog::Info("Enter NTrainingPts: ");
			scanf("%u",&NTrainingPts);
		}
	}	 
	else if(TPO_Method=="MCSimplex")
		NTrainingPts=NPars+1;
	else if(TPO_Method=="MCSimplexPlus1")
		NTrainingPts=NPars+2;
	else if(TPO_Method=="MCQuadratic" || TPO_Method=="MCSphereQuadratic")
		NTrainingPts=(NPars+1)*(NPars+2)/2;
   
	ThetaTrain.clear();
	ThetaTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
	}
	TrainingPtsRead.resize(NTrainingPts);
	TrainingPtsFreeze.resize(NTrainingPts);
	
	ReadTrainingPts();
	FreezeTrainingPts();
	SetTrainingPts();
	
}

void CTPO::Optimize(){
	ALPHA=parmap->getD("TPO_ALPHA",0.01);
	PLUS1=false;
   NMC=parmap->getI("TPO_NMC",0);
	if(NMC==0){
		CLog::Info("Enter NMC: ");
		scanf("%u",&NMC);
	}
	if(TPO_Method=="MC"){
		PLUS1=false;
		Optimize_MC();
	}
	else if(TPO_Method=="MCSphere"){
		PLUS1=true;
		OptimizeSphere_MC();
	}
	else if(TPO_Method=="MCSimplex"){
		PLUS1=false;
		OptimizeSimplex_MC();
	}
	else if(TPO_Method=="MCSimplexPlus1"){
		PLUS1=true;
		OptimizeSimplex_MC();
	}
	else if(TPO_Method=="MCQuadratic"){
		PLUS1=true;
		NTrainingPts=(NPars+1)*(NPars+2)/2;
		Optimize_MC();
	}
	else if(TPO_Method=="MCSphereQuadratic"){
		PLUS1=true;
		NTrainingPts=(NPars+1)*(NPars+2)/2;
		OptimizeSphere_MC();
	}
}

void CTPO::GetC0DDprime(double LAMBDA,vector<double> &theta1,vector<double> &theta2,double &C0,double &D,double &Dprime){
	unsigned int ipar;
	double delThetaSquared,delTheta;
	delThetaSquared=0.0;
	for(ipar=0;ipar<NPars;ipar++){
		delTheta=theta1[ipar]-theta2[ipar];
		delThetaSquared+=delTheta*delTheta;
	}
	C0=exp(-0.5*delThetaSquared/(LAMBDA*LAMBDA));
	D=C0*delThetaSquared;
	Dprime=D*delThetaSquared;
}

double my_erfinv (double a){
	double p, r, t;
	t = fmaf (a, 0.0f - a, 1.0f);
	t = log(t);
	if (fabs(t) > 6.125f) { // maximum ulp error = 2.35793
		p =              3.03697567e-10; //  0x1.4deb44p-32 
		p = fma (p, t,  2.93243101e-8); //  0x1.f7c9aep-26 
		p = fma (p, t,  1.22150334e-6); //  0x1.47e512p-20 
		p = fma (p, t,  2.84108955e-5); //  0x1.dca7dep-16 
		p = fma (p, t,  3.93552968e-4); //  0x1.9cab92p-12 
		p = fma (p, t,  3.02698812e-3); //  0x1.8cc0dep-9 
		p = fma (p, t,  4.83185798e-3); //  0x1.3ca920p-8 
		p = fma (p, t, -2.64646143e-1); // -0x1.0eff66p-2 
		p = fma (p, t,  8.40016484e-1); //  0x1.ae16a4p-1 
	} else { // maximum ulp error = 2.35002
		p =              5.43877832e-9;  //  0x1.75c000p-28 
		p = fma (p, t,  1.43285448e-7); //  0x1.33b402p-23 
		p = fma (p, t,  1.22774793e-6); //  0x1.499232p-20 
		p = fma (p, t,  1.12963626e-7); //  0x1.e52cd2p-24 
		p = fma (p, t, -5.61530760e-5); // -0x1.d70bd0p-15 
		p = fma (p, t, -1.47697632e-4); // -0x1.35be90p-13 
		p = fma (p, t,  2.31468678e-3); //  0x1.2f6400p-9 
		p = fma (p, t,  1.15392581e-2); //  0x1.7a1e50p-7 
		p = fma (p, t, -2.32015476e-1); // -0x1.db2aeep-3 
		p = fma (p, t,  8.86226892e-1); //  0x1.c5bf88p-1 
	}
	r = a * p;
	return r;
}

void CTPO::Optimize_MC(){
	double Sigma2Bar,bestSigma2,dtheta,W11,successrate;
	unsigned int imc,itrain,ipar,nfail=0,nsuccess=0;
   string command;
	FILE *fptr,*fptr_vsNMC;
	vector<vector<double>> besttheta;
	
	CLog::Info("NTrainingpts="+to_string(NTrainingPts)+", NMC="+to_string(NMC)+", LAMBDA="+to_string(LAMBDA)+", ALPHA="+to_string(ALPHA)+"\n");
	dtheta=0.05/sqrt(double(NTrainingPts));
   command="mkdir -p smooth_data/trainingpoint_data";
   system(command.c_str());
   command="mkdir -p smooth_data/trainingpoint_data/Sigma2vsNMC";
   system(command.c_str());
   command="mkdir -p smooth_data/trainingpoint_data/SigmaVsLambda";
   system(command.c_str());
	if(FIRSTCALL==true){
		SetThetaLatinHyperCube(besttheta);
		SetThetaTrain(besttheta);
		FIRSTCALL=false;
	}
	else{
		besttheta=ThetaTrain;
	}
	
	Sigma2Bar=GetSigma2Bar(LAMBDA,ALPHA,W11);
	bestSigma2=Sigma2Bar;
	CLog::Info("Optimize_MC: at beginning: bestSigma2="+to_string(bestSigma2)+"\n");
	
	string filename="smooth_data/trainingpoint_data/Sigma2vsNMC/NTrain"+to_string(NTrainingPts)+".txt";
	fptr_vsNMC=fopen(filename.c_str(),"w");
	fprintf(fptr_vsNMC,"0 %g\n",bestSigma2);
	
	fptr=fopen("smooth_data/trainingpoint_data/Sigma2vsNTrain_Latin.txt","a");
	fprintf(fptr,"%u %g\n",NTrainingPts,bestSigma2);
	fclose(fptr);
	
	for(imc=0;imc<NMC;imc++){
		for(itrain=0;itrain<NTrainingPts;itrain++){
			for(ipar=0;ipar<NPars;ipar++){
				ThetaTrain[itrain][ipar]=besttheta[itrain][ipar]+dtheta*randy->ran_gauss();
			}
		}
		Sigma2Bar=GetSigma2Bar(LAMBDA,ALPHA,W11);
		if(Sigma2Bar<bestSigma2){
			nsuccess+=1;
			bestSigma2=Sigma2Bar;
			for(itrain=0;itrain<NTrainingPts;itrain++){
				besttheta[itrain]=ThetaTrain[itrain];
			}
		}
		else
			nfail+=1;
		for(int n=0;n<100;n++){
			if(imc==pow(2,n)){
				fprintf(fptr_vsNMC,"%u %g\n",imc,bestSigma2);
				n+=1000;
			}
			if(imc<pow(2,n))
				n+=1000;
		}
		if((100*(imc+1)%NMC)==0){
			successrate=double(nsuccess)/double(nsuccess+nfail);
         stringstream fp,bs,sr,dt,sp;
         fp << setprecision(3) << 100.0*(imc+1.0)/double(NMC);
         bs << setprecision(5) << bestSigma2;
         sp << setprecision(5) << 100.0*successrate;
         dt << setprecision(5) << dtheta;
			CLog::Info("+++ finished "+fp.str()+"%, expected accuracy="+bs.str()+", success %="+sp.str()+", step size="+dt.str()+"\n");
			dtheta*=0.05+4.0*successrate;
			nfail=nsuccess=0;
		}
		
	}
	fprintf(fptr_vsNMC,"%u %g\n",imc,bestSigma2);
	fclose(fptr_vsNMC);
   
   /*
    string rthetastring;
    char rthetachar[11];
    double r;
    for(itrain=0;itrain<NTrainingPts;itrain++){
    r=0.0;
    rthetastring="-- itrain="+to_string(itrain)+" --\ntheta=:\n";
    CLog::Info(rthetastring);
    rthetastring="";
    for(ipar=0;ipar<NPars;ipar++){
    r+=pow(besttheta[itrain][ipar],2);
    snprintf(rthetachar,11,"%9.3f ",besttheta[itrain][ipar]);
    rthetastring+=string(rthetachar);
    
    }
    CLog::Info(rthetastring);
    CLog::Info("r="+to_string(sqrt(r))+"\n");
    
    }*/
   fptr=fopen("smooth_data/trainingpoint_data/Sigma2vsNTrain.txt","a");
   fprintf(fptr,"%3d %8.6f\n",NTrainingPts,bestSigma2);
	fclose(fptr);
	
	fptr=fopen("smooth_data/trainingpoint_data/Sigma2vsLambda.txt","a");
	fprintf(fptr,"%8.5f %8.6f\n",LAMBDA,bestSigma2);
	fclose(fptr);
	
	fptr=fopen("smooth_data/trainingpoint_data/Sigma2vsNPars.txt","a");
	fprintf(fptr,"%8.5u %8.6f %u\n",NPars,bestSigma2,NTrainingPts);
	fclose(fptr);
	
	SetThetaTrain(besttheta);

}

void CTPO::OptimizeSimplex_MC(){
	double R=1.2,bestSigma2,Sigma2bar,dR=0.1,W11,bestR,successrate;
	unsigned int imc,nsuccess,nfail,itrain;
	vector<vector<double>> besttheta;
	bestSigma2=1.0E99;
	bestR=R;
	dR=0.2/sqrt(double(NTrainingPts));
	if(PLUS1){
		SetThetaSimplexPlus1(R);
	}
	else
		SetThetaSimplex(R);

	Sigma2bar=GetSigma2Bar(LAMBDA,ALPHA,W11);
	bestSigma2=Sigma2bar;
	CLog::Info("Beginning simplex optimization:\nR="+to_string(bestR)+", bestSigma2="+to_string(bestSigma2)+"\n");
	
	nsuccess=nfail=0;
	for(imc=0;imc<NMC;imc++){
		R=bestR+dR*randy->ran_gauss();
		if(PLUS1)
			SetThetaSimplexPlus1(R);
		else
			SetThetaSimplex(R);
		Sigma2bar=GetSigma2Bar(LAMBDA,ALPHA,W11);
		if(Sigma2bar<bestSigma2){
			bestSigma2=Sigma2bar;
			bestR=R;
			nsuccess+=1;
			for(itrain=0;itrain<NTrainingPts;itrain++)
				besttheta[itrain]=ThetaTrain[itrain];
		}
		else
			nfail+=1;
		if((100*(imc+1)%NMC)==0){
         successrate=double(nsuccess)/double(nsuccess+nfail);
         stringstream bs,sr,sp;
         sp << setprecision(3) << 100.0*(imc+1.0)/double(NMC);
         bs << setprecision(5) << bestSigma2;
         sr << setprecision(5) << 100.0*successrate;
         CLog::Info("+++ finished "+sp.str()+"%, expected accuracy="+bs.str()+", success %="+sr.str()+"\n");
			dR*=0.05+4.0*successrate;
			nfail=nsuccess=0;
		}
	}

	CLog::Info("Best R for simplex="+to_string(bestR)+"\n");
	SetThetaSimplex(bestR);
	SetThetaTrain(besttheta); 
	Sigma2bar=GetSigma2Bar(LAMBDA,ALPHA,W11);
	bestSigma2=Sigma2bar;	
}

void CTPO::OptimizeSphere_MC(){
	double Sigma2Bar,bestSigma2=1.0E99,dtheta,W11,r,R0=1.0,successrate;
	unsigned int imc,itrain,ipar,nfail=0,nsuccess=0;
	vector<vector<double>> besttheta;

	dtheta=0.2/sqrt(double(NTrainingPts));
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			besttheta[itrain][ipar]=R0*randy->ran_gauss();
			if(itrain==0)
				besttheta[itrain][ipar]=0.0;
			if(itrain==1 && ipar!=0){
				besttheta[itrain][ipar]=0.0;
			}
		}
	}
	for(itrain=1;itrain<NTrainingPts;itrain++){
		r=0.0;
		for(ipar=0;ipar<NPars;ipar++){
			r+=besttheta[itrain][ipar]*besttheta[itrain][ipar];
		}
		r=sqrt(r);
		for(ipar=0;ipar<NPars;ipar++){
			besttheta[itrain][ipar]*=R0/r;
		}
	}
	SetThetaTrain(besttheta); 
	Sigma2Bar=GetSigma2Bar(LAMBDA,ALPHA,W11);
	bestSigma2=Sigma2Bar;
	CLog::Info("beginning: bestSigma2="+to_string(bestSigma2)+"\n");
	
	for(imc=0;imc<NMC;imc++){
		for(itrain=1;itrain<NTrainingPts;itrain++){
			r=0.0;
			for(ipar=0;ipar<NPars;ipar++){
				ThetaTrain[itrain][ipar]=besttheta[itrain][ipar]+dtheta*randy->ran_gauss();
				if(itrain==1 && ipar!=0){
					ThetaTrain[itrain][ipar]=0.0;
				}
				r+=ThetaTrain[itrain][ipar]*ThetaTrain[itrain][ipar];	
			}
			r=sqrt(r);
			for(ipar=0;ipar<NPars;ipar++){
				ThetaTrain[itrain][ipar]*=ThetaTrain[1][0]/r;
			}
		}
		Sigma2Bar=GetSigma2Bar(LAMBDA,ALPHA,W11);
		
		if(Sigma2Bar<bestSigma2){
			nsuccess+=1;
			bestSigma2=Sigma2Bar;
			for(itrain=0;itrain<NTrainingPts;itrain++){
				besttheta[itrain]=ThetaTrain[itrain];
			}
		}
		else
			nfail+=1;
		if((100*(imc+1)%NMC)==0){
         successrate=double(nsuccess)/double(nsuccess+nfail);
         stringstream fp,bs,dt,sp,br;
         fp << setprecision(3) << 100.0*(imc+1.0)/double(NMC);
         bs << setprecision(5) << bestSigma2;
         sp << setprecision(5) << 100.0*successrate;
         dt << setprecision(5) << dtheta;
         br << setprecision(5) << fabs(besttheta[1][0]);
			CLog::Info("+++ finished "+fp.str()+", expected accuracy="+bs.str()+", bestR="+br.str()+"\n success %%="+sp.str()+", step size="+dt.str()+"\n");
			dtheta*=0.05+4.0*successrate;
			nfail=nsuccess=0;
		}
		
	}
	SetThetaTrain(besttheta);
}

double CTPO::GetSigma2Bar(double LAMBDA,double ALPHA,double &W11){
	unsigned int a,b;
	Eigen::MatrixXd B,Binv,D,Dprime,BB0,BB1,BB2,bdd,bddi;
	double SigmaA=1.0; // SigmaA shouldn't matter
	double dEdLambda2,Sigma2Bar;
	double bb,dd,ddprime,S20;
	double detfactor=4*sqrt(NTrainingPts);
	
	I.resize(NTrainingPts,NTrainingPts);
	J.resize(NTrainingPts,NTrainingPts);
	K.resize(NTrainingPts,NTrainingPts);
	bdd.resize(NTrainingPts,NTrainingPts);
	bddi.resize(NTrainingPts,NTrainingPts);
	CalcIJK(LAMBDA,priorinfo->ThetaPrior);
	
	B.resize(NTrainingPts,NTrainingPts);
	Binv.resize(NTrainingPts,NTrainingPts);
	D.resize(NTrainingPts,NTrainingPts);
	Dprime.resize(NTrainingPts,NTrainingPts);
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			GetC0DDprime(LAMBDA,ThetaTrain[a],ThetaTrain[b],bb,dd,ddprime);
			B(a,b)=bb; D(a,b)=dd; Dprime(a,b)=ddprime;
		}
	}
	for(a=0;a<NTrainingPts;a++)
		B(a,a)+=ALPHA*ALPHA;
	Binv=B.inverse();
	S20=1.0-(I*Binv).trace();
	if(INCLUDE_LAMBDA_UNCERTAINTY){

		double Iterm,Jterm,Kterm;
		bdd=D*Binv;
		bddi=Binv*bdd;
		Iterm=(I*bddi*bdd).trace();
		Jterm=-(J*bddi).trace();
		Kterm=(K*Binv).trace();
		dEdLambda2=Iterm+Jterm+Kterm;
		dEdLambda2=dEdLambda2/pow(LAMBDA,6);

		if(dEdLambda2<-0.002){
			CLog::Info("dEdLambda2<0, ="+to_string(dEdLambda2)+"\n");
		}
	
		// Calculate d2|B|/dLambda^2
		BB1=B*detfactor;
	

		Eigen::Matrix2d W,Winv;

		W.resize(2,2);
		Winv.resize(2,2);
		Winv(0,0)=double(2*NTrainingPts)/(SigmaA*SigmaA);
		Winv(1,0)=Winv(0,1)=(-1.0/(SigmaA*pow(LAMBDA,3)))*((Binv*D).trace());
		Winv(1,1)=(0.5/pow(LAMBDA,6))*(D*Binv*D*Binv).trace();

		W=Winv.inverse();
		W11=W(1,1);
		if(W(1,1)<0.0){
			CLog::Info("W(1,1)<0, ="+to_string(W(1,1))+"\n");
			Sigma2Bar=1.0E99;
		}
		double S2_duetoLambda=dEdLambda2*W(1,1);
		Sigma2Bar=S20+S2_duetoLambda;
		return Sigma2Bar;
	}
	else
		return S20;
}
