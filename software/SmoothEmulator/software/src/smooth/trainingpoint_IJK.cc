#include "msu_smooth/trainingpoint_optimizer.h"
#include "msu_smooth/priorinfo.h"
#include <cstdlib>

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

double PI=4.0*atan(1.0);

void CTPO::CalcIJK(double Lambda,vector<double> &ThetaPrior){
	unsigned int a,b,ipar;
	double Ii=0,Jafact=0,Jbfact=0,Kabfact=0;
	double Jasum,Jbsum,Kabsum;
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			Jasum=0.0,Jbsum=0.0,Kabsum=0.0;
			I(a,b)=1.0;
			for(ipar=0;ipar<NPars;ipar++){
				if(priorinfo->type[ipar]=="uniform"){
					GetIiJiKiUniform(ThetaPrior[ipar],Lambda,ThetaTrain[a][ipar],ThetaTrain[b][ipar],
					Ii,Jafact,Jbfact,Kabfact);
				}
				else if(priorinfo->type[ipar]=="gaussian"){
					GetIiJiKiGaussian(ThetaPrior[ipar],Lambda,ThetaTrain[a][ipar],ThetaTrain[b][ipar],
					Ii,Jafact,Jbfact,Kabfact);
				}
				else{
					CLog::Fatal("priorinfo->type not gaussian or uniform\n");
				}
				I(a,b)*=Ii;
				Jasum+=Jafact;
				Jbsum+=Jbfact;
				Kabsum+=Kabfact-Jafact*Jbfact;
			}
			J(a,b)=I(a,b)*(Jasum+Jbsum);
			K(a,b)=I(a,b)*(Jasum*Jbsum+Kabsum);
		}
	}
}

void CTPO::CalcIJK_Gaussian(double LAMBDA,double ThetaPriorGauss){ // only works when all have gaussian priors with same ThetaPrior=ThetaPriorDefault
	unsigned int a,b,ipar;
	double lambda,gamma,dT2,tatb,ta2,tb2,Jab,Jba,Iab,Kab,sumt2,Jfacta,Jfactb,Kfact;
	double X,dXdgamma_a,dXdgamma_b,d2Xdgamma_adgamma_b,beta=1.0/(ThetaPriorGauss*ThetaPriorGauss);
	gamma=1.0/(LAMBDA*LAMBDA);
	lambda=2.0*gamma+beta;
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			dT2=tatb=ta2=tb2=sumt2=0.0;
			for(ipar=0;ipar<NPars;ipar++){
				ta2+=pow(ThetaTrain[a][ipar],2);
				tb2+=pow(ThetaTrain[b][ipar],2);
				dT2+=pow(ThetaTrain[a][ipar]-ThetaTrain[b][ipar],2);
				tatb+=ThetaTrain[a][ipar]*ThetaTrain[b][ipar];
				sumt2+=pow(ThetaTrain[a][ipar]+ThetaTrain[b][ipar],2);
			}
			X=gamma*gamma*dT2+beta*gamma*(ta2+tb2);
			Iab=sqrt(pow(beta/lambda,NPars))*exp(-0.5*X/lambda);
			dXdgamma_a=gamma*dT2+beta*ta2;
			dXdgamma_b=gamma*dT2+beta*tb2;
			Jfacta=-(0.5*double(NPars)/lambda)+(0.5*X/(lambda*lambda))-(0.5*dXdgamma_a/lambda);
			Jab=-Jfacta*Iab;
			Jfactb=-(0.5*double(NPars)/lambda)+(0.5*X/(lambda*lambda))-(0.5*dXdgamma_b/lambda);
			Jba=-Jfactb*Iab;
			Jfacta=-0.5*ta2-(0.5*double(NPars)/lambda) -(gamma*gamma*sumt2/(lambda*lambda))+(gamma/lambda)*(ta2+tatb);
			Jab=Iab*Jfacta;
			Jfactb=-0.5*tb2-(0.5*double(NPars)/lambda) -(gamma*gamma*sumt2/(lambda*lambda))+(gamma/lambda)*(tb2+tatb);
			Jba=Iab*Jfactb;
			d2Xdgamma_adgamma_b=dT2;

			Kfact=+(0.5/(lambda*lambda)) -(X/(lambda*lambda*lambda))
				+(0.5*(dXdgamma_a+dXdgamma_b)/(lambda*lambda)) -(0.5*d2Xdgamma_adgamma_b/lambda);
			Kab=Jab*Jba/Iab +(0.5/(lambda*lambda))*Iab -(X/(lambda*lambda*lambda))*Iab
				+(0.5*(dXdgamma_a+dXdgamma_b)/(lambda*lambda))*Iab -(0.5*d2Xdgamma_adgamma_b/lambda)*Iab;
			Kfact=(0.5*double(NPars)/(lambda*lambda)) +(tatb/lambda)-sumt2*(gamma*lambda+gamma*gamma)/pow(lambda,3);
			Kab=(Jab*Jba/Iab)+Iab*Kfact;
			
			Jab*=-2.0;
			Jba*=-2.0;
			Kab*=4.0;
			I(a,b)=Iab;
			J(a,b)=Jab+Jba;
			K(a,b)=Kab;
		}
	}
}

void CTPO::GetIiJiKiGaussian(double ThetaPrior,double Lambda,double theta_a,double theta_b,
double &I,double &Jaterm,double &Jbterm,double &Kabterm){
	double X,gamma,alpha,lambda,deltheta2,sumt2,Jterm;
	gamma=1.0/(Lambda*Lambda);
	alpha=1.0/(ThetaPrior*ThetaPrior);
	lambda=2.0*gamma+alpha;
	X=gamma*gamma*pow(theta_a-theta_b,2)+alpha*gamma*(theta_a*theta_a+theta_b*theta_b);
	deltheta2=(theta_a-theta_b)*(theta_a-theta_b);
	sumt2=theta_a*theta_a+theta_b*theta_b;
	
	I=sqrt(alpha/lambda)*exp(-0.5*X/lambda);
	Jterm=(-1.0/lambda)
		-deltheta2*((gamma*gamma+alpha*gamma)/(lambda*lambda))
			+sumt2*(-0.5*alpha*alpha/(lambda*lambda));	
	Jaterm=0.5*Jterm-0.25*(alpha/lambda)*(theta_a*theta_a-theta_b*theta_b);
	Jbterm=Jterm-Jaterm;
	//Jab=I*Jterm;
	Kabterm=Jaterm*Jbterm+0.5/(lambda*lambda)
		+(0.5*alpha*alpha*sumt2-0.5*(2*gamma*gamma+alpha*lambda)*deltheta2)/pow(lambda,3);
	Kabterm*=4;
	//Kab=4*I*Kterm;
	//Jab*=-2;
	Jaterm*=-2;
	Jbterm*=-2;
}

void CTPO::GetIiJiKiUniform(double ThetaPrior,double Lambda,double theta_a,double theta_b,
                            double &I,double &Jaterm,double &Jbterm,double &Kabterm){
   double Ja,Jb,J,Kab;
   double deltheta,thetabar,rootgamma,Xplus,Xminus,P,Y,W,bplus2,bminus2,deltheta2,bplus,bminus;
   double gamma=1.0/(Lambda*Lambda);
   rootgamma=sqrt(gamma);
   thetabar=0.5*(theta_a+theta_b);
   deltheta=0.5*(theta_a-theta_b);
   deltheta2=deltheta*deltheta;
   bplus=ThetaPrior-thetabar;
   bminus=-ThetaPrior-thetabar;
   bplus2=bplus*bplus;
   bminus2=bminus*bminus;
   Xplus=exp(-gamma*bplus2);
   Xminus=exp(-gamma*bminus2);
   P=(0.5/ThetaPrior)*exp(-gamma*deltheta*deltheta);
   W=(1.0/rootgamma)*(sqrt(PI)/2.0)*(erf(rootgamma*(ThetaPrior-thetabar))-erf(rootgamma*(-ThetaPrior-thetabar)));
   I=P*W;
   
   J=-(0.5/gamma)*I-deltheta2*I;
   Y=bplus*Xplus-bminus*Xminus;
   J+=(0.5/gamma)*P*Y;
   
   Ja=J-P*Xplus*deltheta/gamma+P*Xminus*deltheta/gamma;
   Jb=2.0*J-Ja;
   
   Kab=J*((-0.5/gamma)-deltheta2);
   Kab+=0.5*I/(gamma*gamma);
   Kab=Kab-P*Xplus*bplus*((0.5/(gamma*gamma))+0.5*(deltheta2/gamma)+0.5*bplus*bplus/gamma);
   Kab=Kab+P*Xminus*bminus*((0.5/(gamma*gamma))+0.5*(deltheta2/gamma)+0.5*bminus*bminus/gamma);
   
   Kab-=I*(2.0*deltheta2/gamma);
   Kab=Kab+P*Xplus*bplus*2.0*deltheta2/gamma;
   Kab=Kab-P*Xminus*bminus*2.0*deltheta2/gamma;
   
   J=-2.0*J;
   Jaterm=-Ja/I;
   Jbterm=-Jb/I;
   Kabterm=Kab/I;
   
}
