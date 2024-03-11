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
		Eigen::MatrixXd Ttilde,TtildeInv,ExactVariance;
		vector<vector<double>> T;

		double SigmaA0,SigmaAMin,SigmaA,SigmaATrial,MCStepSize,MCSigmaAStepSize,LAMBDA;
		double A2barRatio,logP;

		unsigned int NMC;   // NMC is for generating independent samplings of A in Tune
		unsigned int NASample;
		bool TuneChooseMCMC,ConstrainA0,CutOffA,UseSigmaY,FirstTune,TuneChooseExact,TuneChooseMCMCPerfect;
		bool pca_ignore;
		vector<vector<double>> ASample;
		vector<double> SigmaASample;
		vector<double> A,ATrial,ABest;
		vector<vector<double>> ThetaTrain;
		vector<vector<double>> B;
		vector<vector<double>> H6,H8;
		Eigen::MatrixXd beta,Psi;

		CSmoothEmulator(string observable_name_set,bool pca_ignore_set);

		void CalcTForTraining();
		void CalcAFromTraining(vector<double> &AA);
		void OldCalcAFromTraining(vector<double> &AA);
		void PrintA(vector<double> &Aprint);

		void SetThetaTrain();
		void Tune();
		void TuneMCMC();
		void TuneMCMC_withSigma();
		void TunePerfectMCMC();
		
		void TuneExact();
		void GenerateUncertaintyMatrices();
		void GetExactQuantities();
		void GetExactSigmaA();
		void CalcExactLogP();
		
		double GetLog_AProb(vector<double> &AA,double SigmaA);

		void SetA_Zero(vector<double> &A);
		void SetA_RanGauss(double ASigmaA,vector<double> &AA);
		void SetA_Constant(double ASigmaA,vector<double> &AA);
		void SetA_RanSech(double ASigmaA,vector<double> &AA);

		double SigmaAbar;
		unsigned int NSigmaA;

		void GenerateASamples();
		void CalcY(CModelParameters *modpars,double &Y,double &SigmaY);
		void CalcY(vector<double> &Theta,double &Y,double &SigmaY);
		//void CalcYDYDTheta(CModelParameters *modpars,double &Y,vector<double> &dYdTheta,double &SigmaY);
		void CalcYDYDTheta(CModelParameters *modpars,double &Y,vector<double> &dYdTheta,double &SigmaY);
		void CalcYDYDTheta(vector<double> &Theta,double &Y,vector<double> &dYdTheta,double &SigmaY);
		void WriteCoefficients();
		void ReadCoefficients();
		void GetExactUncertainty(vector<double> &Theta_s,double &sigma);

		void Init();

		static 
		static CSmoothMaster *smoothmaster;
		static unsigned int NPars;
		static CSmooth *smooth;
		static CparameterMap *parmap;
		static Crandy *randy;
		static unsigned int NTrainingPts;

	};

};

#endif
