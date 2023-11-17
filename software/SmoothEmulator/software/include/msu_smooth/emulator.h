#ifndef __EMULATOR_H__
#define __EMULATOR_H__
#include <cstdlib>
#include <cmath>

#include <cstdio>
#include <vector>
#include <Eigen/Dense>
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/constants.h"
#include "msu_smooth/smooth.h"
#include "msu_commonutils/log.h"
#include "msu_smooth/master.h"
#include "msu_smooth/modelparinfo.h"
#include "msu_smooth/observableinfo.h"
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/traininginfo.h"
using namespace NMSUPratt;

namespace NBandSmooth{
	class CSmoothMaster;
	class CTrainingInfo;
	class CPriorInfo;
	class CObservableInfo;
	class CModelParInfo;

	class CSmoothEmulator{
	public:
		int iY; // labels observable from observable info
		string observable_name;
		Eigen::MatrixXd M,Minv;
		vector<vector<double>> Mtot;

		double SigmaA0,SigmaAMin,SigmaA,SigmaATrial,MCStepSize,MCSigmaAStepSize,LAMBDA;
		unsigned int NMC;   // NMC is for generating independent samplings of A in Tune
		unsigned int NASample;
		bool TuneChooseMCMC,ConstrainA0,CutOffA,UseSigmaYReal,FirstTune;
		vector<vector<double>> ASample;
		vector<double> SigmaASample;
		vector<double> A,ATrial;
		vector<vector<double>> ThetaTrain;

		CSmoothEmulator(string observable_name_set);

		void CalcMForTraining();
		void CalcAFromTraining(vector<double> &AA);
		void OldCalcAFromTraining(vector<double> &AA);
		void PrintA(vector<double> &Aprint);

		void SetThetaTrain();
		void Tune();
		void TuneMCMC();
		void TuneMCMC_withSigma();
		void TunePerfect();
		double GetLog_AProb(vector<double> &AA,double SigmaA);

		void SetA_Zero(vector<double> &A);
		void SetA_RanGauss(double ASigmaA,vector<double> &AA);
		void SetA_Constant(double ASigmaA,vector<double> &AA);
		void SetA_RanSech(double ASigmaA,vector<double> &AA);

		double SigmaAbar;
		int NSigmaA;

		void GenerateASamples();
		void CalcY(CModelParameters *modpars,double &Y,double &SigmaY);
		void CalcY(vector<double> Theta,double &Y,double &SigmaY);
		void WriteCoefficients();
		void ReadCoefficients();

		void Init();

		static CSmoothMaster *smoothmaster;
		static int NPars;
		static CSmooth *smooth;
		static CparameterMap *parmap;
		static Crandy *randy;
		static unsigned int NTrainingPts;


	};

}

#endif
