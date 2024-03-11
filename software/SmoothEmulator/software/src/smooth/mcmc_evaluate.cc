#include <complex>
#include "msu_smooth/master.h"
#include "msu_smooth/mcmc.h"
#include "msu_smoothutils/log.h"
using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

void CMCMC::EvaluateTrace(){
	unsigned int itrace,ipar,jpar,ntrace=trace.size();
	bool UsePCA=master->UsePCA;
	CPCA *pca=NULL;
	if(UsePCA){
			pca=new CPCA(master->parmap);
			pca->ReadTransformationInfo();
	}
	NObs=master->observableinfo->NObservables;
	vector<double> thetabar,YRMS;
	FILE *fptr;
	string SigmaString;
	char cc[CLog::CHARLENGTH];
	Eigen::MatrixXd CovThetaTheta,CovThetaY,CovYY,CovThetaThetaInv,dYdTheta,RP;
	vector<double> Y,Ybar,Z,Zbar,SigmaYEmulator,SigmaZEmulator;
	double sigmatot2;
	Eigen::VectorXcd evals;
	Eigen::MatrixXcd evecs;
	unsigned int iobs,jobs;
	thetabar.resize(NPars);
	CovThetaTheta.resize(NPars,NPars);
	CovThetaThetaInv.resize(NPars,NPars);
	CovThetaY.resize(NPars,NObs);
	dYdTheta.resize(NPars,NObs);
	evecs.resize(NPars,NPars);
	SigmaYEmulator.resize(NObs);
	CovYY.resize(NObs,NObs);
	RP.resize(NObs,NPars);
	CovThetaTheta.setZero();
	CovThetaY.setZero();
	CovYY.setZero();
	Y.resize(NObs);
	Ybar.resize(NObs);
	YRMS.resize(NObs);
	SigmaYEmulator.resize(NObs);

	if(UsePCA){
		Z.resize(NObs);
		Zbar.resize(NObs);
		SigmaZEmulator.resize(NObs);
	}
	
	for(ipar=0;ipar<NPars;ipar++){
		thetabar[ipar]=0.0;
	}
	for(iobs=0;iobs<NObs;iobs++){
		Ybar[iobs]=0.0;
		if(UsePCA)	
			Zbar[iobs]=0.0;
	}

	for(itrace=0;itrace<ntrace;itrace++){
		master->CalcAllY(trace[itrace],Y,SigmaYEmulator);
		if(UsePCA){
			for(iobs=0;iobs<NObs;iobs++){
				SigmaZEmulator[iobs]=SigmaYEmulator[iobs];
				Z[iobs]=Y[iobs];
			}
			pca->TransformZtoY(Z,SigmaZEmulator,Y,SigmaYEmulator);
		}
		for(ipar=0;ipar<NPars;ipar++){
			thetabar[ipar]+=trace[itrace][ipar];
		}
		for(iobs=0;iobs<NObs;iobs++){
			Ybar[iobs]+=Y[iobs];
			if(UsePCA)
				Zbar[iobs]+=Z[iobs];
		}
		for(ipar=0;ipar<NPars;ipar++){
			for(jpar=0;jpar<NPars;jpar++){
				CovThetaTheta(ipar,jpar)+=trace[itrace][ipar]*trace[itrace][jpar];
			}
		}
		for(ipar=0;ipar<NPars;ipar++){
			for(iobs=0;iobs<NObs;iobs++){
				CovThetaY(ipar,iobs)+=trace[itrace][ipar]*Y[iobs];
			}
		}
		for(iobs=0;iobs<NObs;iobs++){
			for(jobs=0;jobs<NObs;jobs++){
				CovYY(iobs,jobs)+=Y[iobs]*Y[jobs];
			}
		}
	}
	
	for(ipar=0;ipar<NPars;ipar++){
		thetabar[ipar]=thetabar[ipar]/double(ntrace);
		for(jpar=0;jpar<NPars;jpar++){
			CovThetaTheta(ipar,jpar)=CovThetaTheta(ipar,jpar)/double(ntrace);
		}
		for(iobs=0;iobs<NObs;iobs++){
			CovThetaY(ipar,iobs)=CovThetaY(ipar,iobs)/double(ntrace);
		}
	}
	CovThetaThetaInv=CovThetaTheta.inverse();
	
	for(iobs=0;iobs<NObs;iobs++){
		Ybar[iobs]=Ybar[iobs]/double(ntrace);
		if(UsePCA)
			Zbar[iobs]=Zbar[iobs]/double(ntrace);
		for(jobs=0;jobs<NObs;jobs++){
			CovYY(iobs,jobs)=CovYY(iobs,jobs)/double(ntrace);
		}
	}
	
	for(ipar=0;ipar<NPars;ipar++){
		for(jpar=0;jpar<NPars;jpar++){
			CovThetaTheta(ipar,jpar)=CovThetaTheta(ipar,jpar)-thetabar[ipar]*thetabar[jpar];
		}
		for(iobs=0;iobs<NObs;iobs++){
			CovThetaY(ipar,iobs)=CovThetaY(ipar,iobs)-thetabar[ipar]*Ybar[iobs];
		}
	}
	
	for(iobs=0;iobs<NObs;iobs++){
		for(jobs=0;jobs<NObs;jobs++){
			CovYY(iobs,jobs)=CovYY(iobs,jobs)-Ybar[iobs]*Ybar[jobs];
		}
	}
	/*
	cout << "<<ThetaTheta>>\n";
	cout << CovThetaTheta << endl;
	cout << "<<ThetaTheta>>^-1\n";
	cout << CovThetaThetaInv << endl;
	cout << "<<ThetaY>>\n";
	cout << CovThetaY << endl;
	cout << "<<YY>>\n";
	cout << CovYY << endl;
	*/
	
	for(iobs=0;iobs<NObs;iobs++){
		YRMS[iobs]=0.0;
		for(ipar=0;ipar<NPars;ipar++){
			dYdTheta(iobs,ipar)=0.0;
			for(jpar=0;jpar<NPars;jpar++){
				dYdTheta(iobs,ipar)+=CovThetaY(jpar,iobs)*CovThetaThetaInv(ipar,jpar);
			}
		}
		for(ipar=0;ipar<NPars;ipar++){
			YRMS[iobs]+=dYdTheta(iobs,ipar)*dYdTheta(iobs,ipar);
		}
		YRMS[iobs]=sqrt(YRMS[iobs]);
	}
	for(iobs=0;iobs<NObs;iobs++){
		sigmatot2=master->observableinfo->SigmaExp[iobs]*master->observableinfo->SigmaExp[iobs]
			+SigmaYEmulator[iobs]*SigmaYEmulator[iobs];
		for(ipar=0;ipar<NPars;ipar++){
			RP(ipar,iobs)=CovThetaY(ipar,iobs)*YRMS[iobs]/sigmatot2;
		}
	}

	
	// Write output
	
	CModelParameters modpars;
	modpars.priorinfo=master->priorinfo;
	modpars.SetTheta(thetabar);
	modpars.TranslateTheta_to_X();
	//modpars.Print();
	string command="mkdir -p mcmc_trace";
	system(command.c_str());
	modpars.Write("mcmc_trace/xbar_thetabar.txt");
	
	fptr=fopen("mcmc_trace/CovThetaTheta.txt","w");
	for(ipar=0;ipar<NPars;ipar++){
		SigmaString.clear();
		for(jpar=0;jpar<NPars;jpar++){
			CovThetaTheta(ipar,jpar)=CovThetaTheta(ipar,jpar)-thetabar[ipar]*thetabar[jpar];
			snprintf(cc,CLog::CHARLENGTH,"%12.5e ",CovThetaTheta(ipar,jpar));
			SigmaString=SigmaString+cc;
		}
		SigmaString=SigmaString+"\n";
		fprintf(fptr,"%s",SigmaString.c_str());
	}
	fclose(fptr);
	
	fptr=fopen("mcmc_trace/ResolvingPower.txt","w");
	for(ipar=0;ipar<NPars;ipar++){
		for(iobs=0;iobs<NObs;iobs++){
			fprintf(fptr,"%12.5e ",RP(ipar,iobs));
		}
		fprintf(fptr,"\n");
	}
	fclose(fptr);
	
	Eigen::EigenSolver<Eigen::MatrixXd> esolver(CovThetaTheta);
	evals=esolver.eigenvalues();
	evecs=esolver.eigenvectors();
	vector<double> evalnorm;
	evalnorm.resize(NPars);
	fptr=fopen("mcmc_trace/CovThetaTheta_eigenvals.txt","w");
	for(ipar=0;ipar<NPars;ipar++){
		evalnorm[ipar]=sqrt(fabs(real(evals(ipar))));
		fprintf(fptr,"%15.8e\n",evalnorm[ipar]);
	}
	fclose(fptr);
	
	fptr=fopen("mcmc_trace/CovThetaTheta_eigenvecs.txt","w");
	for(ipar=0;ipar<NPars;ipar++){
		SigmaString.clear();
		for(jpar=0;jpar<NPars;jpar++){
			snprintf(cc,CLog::CHARLENGTH,"%12.5e ",real(evecs(ipar,jpar)));
			SigmaString+=cc;
		}
		SigmaString+="\n";
		fprintf(fptr,"%s",SigmaString.c_str());
	}
	fclose(fptr);
		
	
}