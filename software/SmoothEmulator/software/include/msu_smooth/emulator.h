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
		unsigned int iY; // labels observable from observable info
		string observable_name;

		double LAMBDA,SigmaA;
		double A2barRatio,logP;

		unsigned int NMC;   // NMC is for generating independent samplings of A in Tune
		unsigned int NASample;
		bool ConstrainA0,UseSigmaY;
		bool pca_ignore;
		vector<double> ABest;
		vector<vector<double>> ThetaTrain,TTrain;
		Eigen::MatrixXd B,Binv;
		//vector<vector<double>> H6,H8;
		//Eigen::MatrixXd beta,Psi;

		CSmoothEmulator(string observable_name_set,bool pca_ignore_set);

		void CalcTForTraining();
		void PrintA(vector<double> &Aprint);

		void SetThetaTrain();
		void Tune();
		void GetSigmaA();
		void CalcExactLogP();
		
		double GetLog_AProb(vector<double> &AA);

		void SetA_Zero(vector<double> &A);
		void SetA_RanGauss(double ASigmaA,vector<double> &AA);
		void SetA_Constant(double ASigmaA,vector<double> &AA);
		void SetA_RanSech(double ASigmaA,vector<double> &AA);

		//void GenerateASamples();
		double GetYOnly(CModelParameters *modpars);
		double GetYOnly(vector<double> &Theta);
		double GetUncertainty(CModelParameters *modpars);
		double GetUncertainty(vector<double> &Theta_s);
		void CalcYAndUncertainty(vector<double> &Theta_s,double &Y,double &uncertainty);
		
		
		
		//void CalcYDYDTheta(CModelParameters *modpars,double &Y,vector<double> &dYdTheta,double &SigmaY);
		void CalcYDYDTheta(CModelParameters *modpars,double &Y,vector<double> &dYdTheta,double &SigmaY);
		void CalcYDYDTheta(vector<double> &Theta,double &Y,vector<double> &dYdTheta,double &SigmaY);
		void WriteCoefficients();
		void ReadCoefficients();
		
		void Init();

		static CSmoothMaster *smoothmaster;
		static unsigned int NPars;
		static CSmooth *smooth;
		static CparameterMap *parmap;
		static Crandy *randy;
		static unsigned int NTrainingPts;

	};

};

#endif
