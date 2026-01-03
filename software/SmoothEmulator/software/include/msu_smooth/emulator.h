#ifndef __EMULATOR_H__
#define __EMULATOR_H__
#include <cstdlib>
#include <cmath>

#include <cstdio>
#include <vector>
#include <Eigen/Dense>
#include "msu_smoothutils/parametermap.h"
#include "msu_smoothutils/misc.h"
#include "msu_smoothutils/randy.h"
#include "msu_smooth/smooth.h"
#include "msu_smoothutils/log.h"
#include "msu_smooth/master.h"
#include "msu_smooth/modelparinfo.h"
#include "msu_smooth/observableinfo.h"
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/traininginfo.h"
using namespace NMSUUtils;

namespace NBandSmooth{
	class CSmoothMaster;
	class CTrainingInfo;
	class CPriorInfo;
	class CObservableInfo;
	class CModelParInfo;

	class CSmoothEmulator{
	public:
      CSmoothEmulator(string observable_name_set);
      friend class CSmoothMaster;
   private:
      bool INCLUDE_LAMBDA_UNCERTAINTY,FIXLAMBDA;
		unsigned int iY; // labels observable from observable info
		string observable_name;

		double LAMBDA,SigmaA,ALPHA,LambdaVariance;
		double logP,detBB,d2lndetBBdLambda2;
		vector<vector<double>> ThetaTrain;
		Eigen::MatrixXd B,Binv;
		Eigen::MatrixXd Bprime,Bprimeprime;
		Eigen::VectorXd chi,chiprime;
		Eigen::Matrix2d W,Winv;
      void Init();
		double GetCorrelation(vector<double> &theta1,vector<double> &theta2);
		void CalcWBprimeChi();
		double GetSigma2_Lambda(vector<double> &theta);
		void SetThetaTrain();
		void Tune();
		void Tune(double LambdaSet); // fix Lambda
		void CalcSigmaA();
		void CalcSigmaALambda();
		void CalcLambdaVariance();
		void CalcLogP();
		void CalcB();
		void GetYAndUncertaintyFromTheta(vector<double> &Theta,double &Y,double &uncertainty);
		static CSmoothMaster *smoothmaster;
		static unsigned int NPars;
		static CparameterMap *parmap;
   public:
		static unsigned int NTrainingPts,NTestingPts;

	};

};

#endif
