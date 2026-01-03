#ifndef __SMOOTH_MCMC_H__
#define __SMOOTH_MCMC_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include <Eigen/Dense>
#include "msu_smoothutils/parametermap.h"
#include "msu_smoothutils/misc.h"
#include "msu_smoothutils/randy.h"
#include "msu_smoothutils/log.h"
#include "msu_smooth/master.h"


using namespace std;
using namespace NMSUUtils;
using namespace NBandSmooth;

namespace NBandSmooth{
	
	class CLLCalc;
	class CLLCalcSmooth;

	class CMCMC{
	public:
		//bool OPTIMIZESTEPS;
		bool IGNORE_EMULATOR_ERROR;
		CparameterMap *parmap;
		CSmoothMaster *master;
		Crandy *randy;
		
		CMCMC();
		CMCMC(CSmoothMaster *master);
		unsigned int NPars,NObs;
		vector<vector<double>> trace;
		string trace_filename,Xtrace_filename;
		double stepsize;
		
		void ClearTrace(); // erases trace info so one can start over, resets at theta=0.
		void PruneTrace(); // erases trace, except for last point
		
		void PerformTrace(unsigned int Ntrace,unsigned int Nskip);
		void PerformMetropolisTrace(unsigned int Ntrace,unsigned int Nskip);
		void WriteTrace();
		void ReadTrace();
		void EvaluateTrace();
		
		//void OptimizeSteps();
		Eigen::VectorXcd stepvec,stepvecprime,dTdTEigenVals;
		Eigen::MatrixXd dThetadTheta;
		Eigen::MatrixXcd dTdTEigenVecs;

		void CalcLL(vector<double> &theta,double &LL);
		//void CalcLLPlusDerivatives(vector<double> &theta,double &LL,vector<double> &dLL_dtheta);
		CLLCalcSmooth *llcalc;
		static CPriorInfo *priorinfo;
	};
	
	class CLLCalc{
	public:
		static bool IGNORE_EMULATOR_ERROR;
		CLLCalc();
		CLLCalc(CSmoothMaster *master);
		unsigned int NPars,NObs;
		vector<double> Y,SigmaY,SigmaY_emulator;
		vector<vector<double>> dYdTheta;
		CObservableInfo *obsinfo;
		CSmoothMaster *master;
		virtual void CalcLL(vector<double> &theta,double &LL);
		//virtual void CalcLLPlusDerivatives(vector<double> &theta,double &LL,vector<double> &dLL_dtheta);
		static CPriorInfo *priorinfo;
	};
	
	class CLLCalcSmooth : public CLLCalc{
	public:
		CLLCalcSmooth(CSmoothMaster *master);
		void CalcLL(vector<double> &theta,double &LL);
		//void CalcLLPlusDerivatives(vector<double> &theta,double &LL,vector<double> &dLL_dtheta);
	};

};

#endif
