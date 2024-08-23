#include "msu_smooth/master.h"
#include "msu_smooth/mcmc.h"
#include "msu_smoothutils/misc.h"

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

CPriorInfo *CMCMC::priorinfo=NULL;

CMCMC::CMCMC(){
	//
}

CMCMC::CMCMC(CSmoothMaster *master_set){
	master=master_set;
	priorinfo=master->priorinfo;
	CLLCalc::priorinfo=priorinfo;
	parmap=master->parmap;
	parmap->ReadParsFromFile("smooth_data/smooth_parameters/mcmc_parameters.txt");
	string logfilename=parmap->getS("MCMC_LogFileName","Screen");
	if(logfilename!="Screen"){
		CLog::Init(logfilename);
	}
	randy=master->randy;
	NPars=master->NPars;
	trace_filename=parmap->getS("MCMC_TRACE_FILENAME","mcmc_trace/trace.txt");
	Xtrace_filename=parmap->getS("MCMC_TRACE_FILENAME","mcmc_trace/Xtrace.txt");
	string command="mkdir -p smooth_data/mcmc_trace";
	system(command.c_str());
	OPTIMIZESTEPS=parmap->getB("MCMC_OPTIMIZESTEPS",false);
	IGNORE_EMULATOR_ERROR=parmap->getB("MCMC_IGNORE_EMULATOR_ERROR",false);
	CLLCalc::IGNORE_EMULATOR_ERROR=IGNORE_EMULATOR_ERROR;
	langevin=parmap->getB("MCMC_LANGEVIN",false);
	if(langevin)
		stepsize=parmap->getD("MCMC_LANGEVIN_STEPSIZE",0.01);
	else
		stepsize=parmap->getD("MCMC_METROPOLIS_STEPSIZE",0.05);
	ClearTrace();
	llcalc=new CLLCalcSmooth(master);
	
	OPTIMIZESTEPS=false;
	dThetadTheta.resize(NPars,NPars);
	dTdTEigenVecs.resize(NPars,NPars);
	dTdTEigenVals.resize(NPars);
	stepvec.resize(NPars);
	stepvecprime.resize(NPars);
	dTdTEigenVecs.setZero();
	for(unsigned int ipar=0;ipar<NPars;ipar++){
		dTdTEigenVecs(ipar,ipar)=1.0;
		dTdTEigenVals(ipar)=stepsize;
	}
}

void CMCMC::ClearTrace(){ // clears trace but adds back first point with theta=0
	for(unsigned int itrace=0;itrace<trace.size();itrace++){
		trace[itrace].clear();
	}
	trace.clear();
	trace.resize(1);
	trace[0].resize(NPars);
	for(unsigned int ipar=0;ipar<NPars;ipar++){
		trace[0][ipar]=0.0;
	}
}

// erases trace, but keeps last point
void CMCMC::PruneTrace(){
	vector<double> theta0;
	theta0=trace[trace.size()-1];
	ClearTrace();
	for(unsigned int ipar=0;ipar<NPars;ipar++){
		trace[0][ipar]=theta0[ipar];
	}
}

void CMCMC::PerformTrace(unsigned int Ntrace,unsigned int Nskip){
	//if(langevin)
	//	PerformLangevinTrace(Ntrace,Nskip);
	//else
		PerformMetropolisTrace(Ntrace,Nskip);
}

void CMCMC::PerformMetropolisTrace(unsigned int Ntrace,unsigned int Nskip){
	unsigned long long int nsuccess=0;
	unsigned int itrace,iskip,ipar,it0;
	double oldLL,newLL,X,bestLL;
	vector<double> oldtheta,newtheta;
	oldtheta.resize(NPars);
	newtheta.resize(NPars);
	OPTIMIZESTEPS=false;
	vector<double> besttheta;
	besttheta.resize(NPars);
	
	it0=trace.size();
	if(trace.size()==0){
		CLog::Fatal("Inside MCMC::PerformMetropolis, no initial point!\n");
	}
	trace.resize(trace.size()+Ntrace);
	if(it0==0){
		for(ipar=0;ipar<NPars;ipar++){
			oldtheta[ipar]=0.2;
		}
	}
	else{
		oldtheta=trace[it0-1];
	}
	llcalc->CalcLL(oldtheta,oldLL);
	besttheta=oldtheta;
	
	bestLL=oldLL;
	CLog::Info("At beginning of Trace, LL="+to_string(oldLL)+"\n");
	
	for(itrace=it0;itrace<it0+Ntrace;itrace++){
		for(iskip=0;iskip<Nskip;iskip++){
			if(OPTIMIZESTEPS){
				for(ipar=0;ipar<NPars;ipar++){
					stepvecprime[ipar]=stepsize*randy->ran_gauss();
				}
				stepvec=dTdTEigenVecs*stepvecprime;
				for(ipar=0;ipar<NPars;ipar++){
					newtheta[ipar]=oldtheta[ipar]+real(stepvec(ipar));
				}
			}
			else{
				for(ipar=0;ipar<NPars;ipar++){
					newtheta[ipar]=oldtheta[ipar]+stepsize*randy->ran_gauss();
				}
			}
			
			llcalc->CalcLL(newtheta,newLL);
			if(newLL>bestLL){
				bestLL=newLL;
				for(ipar=0;ipar<NPars;ipar++)
					besttheta[ipar]=newtheta[ipar];
			}
			
			if(newLL>oldLL){
				for(ipar=0;ipar<NPars;ipar++)
					oldtheta[ipar]=newtheta[ipar];
				oldLL=newLL;
				nsuccess+=1;
			}
			else{
				X=newLL-oldLL;
				if(X>-30.0 && X<0.0){
					if(exp(X)>randy->ran()){
						for(ipar=0;ipar<NPars;ipar++)
							oldtheta[ipar]=newtheta[ipar];
						oldLL=newLL;
						nsuccess+=1;
					}
				}
			}
		}
		trace[itrace].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++){
			trace[itrace][ipar]=oldtheta[ipar];
		}
		if(Ntrace>10 && ((itrace+1)*10)%Ntrace==0){
			CLog::Info("finished "+to_string(lrint(100*double(itrace+1)/double(Ntrace)))+"%\n");
		}
	}
	CLog::Info("At end of trace, best LL="+to_string(bestLL)+"\nBest Theta=\n");
	for(ipar=0;ipar<NPars;ipar++){
		CLog::Info(to_string(besttheta[ipar])+"  ");
	}
	CLog::Info("\nMetropolis success percentage="+to_string(100.0*double(nsuccess)/(double(Ntrace)*double(Nskip)))+"\n");
}

/*(void CMCMC::PerformLangevinTrace(unsigned int Ntrace,unsigned int Nskip){
	unsigned int itrace,iskip,ipar;
	bool inside;
	double LL,ss,sqstep,dLLdTheta2,bestLL=-1.0E99;
	vector<double> dLLdTheta,dTheta;
	vector<double> *oldptr,*newmodpars,*modpars;
	if(trace.size()==0){
		CLog::Fatal("Inside MCMC::PerformMetropolis, no initial point!\n");
	}
	dLLdTheta.resize(NPars);
	dTheta.resize(NPars);
	llcalc->CalcLL(&trace[trace.size()-1],bestLL);
	CLog::Info("At beginning of trace LL="+to_string(bestLL)+"\n");
	
	for(itrace=0;itrace<Ntrace;itrace++){
		oldptr=&trace[trace.size()-1];
		newmodpars=new CModelParameters();
		for(iskip=0;iskip<Nskip;iskip++){
			llcalc->CalcLLPlusDerivatives(oldptr,LL,dLLdTheta);
			if(LL>bestLL){
				bestLL=LL;
			}
			inside=true;
			dLLdTheta2=0.0;
			for(ipar=0;ipar<NPars;ipar++){
				dLLdTheta2+=dLLdTheta[ipar]*dLLdTheta[ipar];
			}
			ss=stepsize/sqrt(dLLdTheta2);
			sqstep=sqrt(2.0*ss);
			for(ipar=0;ipar<NPars;ipar++){
				dTheta[ipar]=ss*dLLdTheta[ipar]+sqstep*randy->ran_gauss();
				if(dTheta[ipar]!=dTheta[ipar]){
					oldptr->Print();
					CLog::Fatal("disaster in PerformLangevinTrace\n");
				}
				newmodpars->Theta[ipar]=oldptr->Theta[ipar]+dTheta[ipar];
				if(fabs(newmodpars->Theta[ipar])>1.0){
					inside=false;
					//ipar=NPars;
				}
			}
			
			if(inside){
				*oldptr=*newmodpars;
			}
		}
		modpars=new CModelParameters();
		*modpars=*oldptr;
		modpars->TranslateTheta_to_X();
		trace.push_back(*modpars);
		//CLog::Info("finished for itrace"+to_string(itrace)+"\n");
		delete newmodpars;
		if(Ntrace>10 && ((itrace+1)*10)%Ntrace==0){
			CLog::Info("finished "+to_string(lrint(100*double(itrace+1)/double(Ntrace)))+"%\n");
		}
	}
	CLog::Info("At end of trace: best LL="+to_string(bestLL)+"\n");
	
}
*/

void CMCMC::WriteTrace(){
	unsigned int itrace,ipar,ntrace=trace.size();
	FILE *fptr;
	CLog::Info("writing Theta values, ntrace = "+to_string(ntrace)+"\n");
	fptr=fopen(("smooth_data/"+trace_filename).c_str(),"w");
	for(itrace=0;itrace<ntrace;itrace++){
		for(ipar=0;ipar<NPars;ipar++){
			fprintf(fptr,"%8.5f ",trace[itrace][ipar]);
		}
		fprintf(fptr,"\n");
	}
	fclose(fptr);
}

void CMCMC::WriteXTrace(){
	CModelParameters modpars;
	
	unsigned int itrace,ipar,ntrace=trace.size();
	FILE *fptr;
	CLog::Info("writing X values, ntrace = "+to_string(ntrace)+"\n");
	fptr=fopen(("smooth_data/"+Xtrace_filename).c_str(),"w");
	for(itrace=0;itrace<ntrace;itrace++){
		for(ipar=0;ipar<NPars;ipar++){
			modpars.Theta[ipar]=trace[itrace][ipar];
		}
		modpars.TranslateTheta_to_X();
		for(ipar=0;ipar<NPars;ipar++){
			fprintf(fptr,"%8.5f ",modpars.X[ipar]);
		}
		fprintf(fptr,"\n");
	}
	fclose(fptr);
}

void CMCMC::ReadTrace(){
	ClearTrace();
	unsigned int ipar,itrace=0;
	vector<double> thetaread;
	thetaread.resize(NPars);
	FILE *fptr;
	fptr=fopen(trace_filename.c_str(),"r");
	while(!feof(fptr)){
		for(ipar=0;ipar<NPars;ipar++){
			fscanf(fptr,"%lf",&thetaread[ipar]);
		}
		if(!feof(fptr)){
			if(itrace!=0){
				trace.resize(trace.size()+1);
			}
			for(ipar=0;ipar<NPars;ipar++)
				trace[itrace][ipar]=thetaread[ipar];
		}
		itrace+=1;
	}
	fclose(fptr);
	CLog::Info("read in "+to_string(trace.size())+" trace elements\n");
}

