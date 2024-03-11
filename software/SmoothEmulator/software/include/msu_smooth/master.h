#ifndef __SMOOTH_MASTER_H__
#define __SMOOTH_MASTER_H__
#include <cstdlib>
#include <cmath>

#include <cstdio>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <sstream>
#include "msu_smoothutils/parametermap.h"
#include "msu_smoothutils/misc.h"
#include "msu_smoothutils/randy.h"
#include "msu_smooth/pca.h"
#include "msu_smooth/emulator.h"
#include "msu_smooth/modelparinfo.h"
#include "msu_smooth/smooth.h"
#include "msu_smoothutils/log.h"
#include "msu_smooth/observableinfo.h"
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/traininginfo.h"

using namespace NMSUUtils;

namespace NBandSmooth{
	class CSmoothEmulator;
	class CTrainingInfo;
	class CPriorInfo;
	class CObservableInfo;
	class CModelParInfo;
	
	class CSmoothMaster{
	public:
		unsigned int TrainType;
		CSmoothMaster(CparameterMap *parmap_set);
		CparameterMap *parmap;
		unsigned int NPars;
		vector<CSmoothEmulator *> emulator;
		CTrainingInfo *traininginfo;
		CObservableInfo *observableinfo;
		CPriorInfo *priorinfo;
		Crandy *randy;
		CSmooth *smooth;
		string ModelRunDirName,CoefficientsDirName;
		bool UsePCA;
		vector<bool> pca_ignore;
		double pca_minvariance;

		void ReadTrainingInfo();
		void GenerateCoefficientSamples();
		void TuneAllY(); // tune all observables
		void TuneY(string obsname); // tune one observable
		void TuneY(unsigned int iY); // tune one observable
		void SetThetaTrain();
		
		void CalcAllY(CModelParameters *modelpars,vector<double> &Y,vector<double> &SigmaY_emulator);
		void CalcAllY(vector<double> &theta,vector<double> &Y,vector<double> &SigmaY_emulator);
		void CalcY(unsigned int iY,CModelParameters *modelpars,double &Y,double &SigmaY_emulator);
		void CalcY(unsigned int iY,vector<double> &theta,double &Y,double &SigmaY_emulator);
		void CalcY(string obsname,CModelParameters *modelpars,double &Y,double &SigmaY_emulator);
		void CalcY(string obsname,vector<double> &theta,double &Y,double &SigmaY_emulator);
		
		void CalcAllYdYdTheta(CModelParameters *modelpars,vector<double> &Y,
		vector<double> &SigmaY_emulator,vector<vector<double>> &dYdTheta);
		void CalcAllYdYdTheta(vector<double> &theta,vector<double> &Y,
		vector<double> &SigmaY_emulator,vector<vector<double>> &dYdTheta);
		void CalcYdYdTheta(string obsname,CModelParameters *modelpars,double &Y,
		double &SigmaY_emulator,vector<double> &dYdTheta);
		void CalcYdYdTheta(string obsname,vector<double> &theta,double &Y,
		double &SigmaY_emulator,vector<double> &dYdTheta);
		void CalcYdYdTheta(unsigned int iY,CModelParameters *modelpars,double &Y,
		double &SigmaY_emulator,vector<double> &dYdTheta);
		void CalcYdYdTheta(unsigned int iY,vector<double> &theta,double &Y,
		double &SigmaY_emulator,vector<double> &dYdTheta);
		
		void CalcAllLogP();
		void TestAtTrainingPts();
		void TestAtTrainingPts(string obsname);
		void TestAtTrainingPts(unsigned int iY);
		void TestVsFullModel();
		void WriteCoefficientsAllY();
		void WriteCoefficients(string obsname);
		void WriteCoefficients(unsigned int iY);
		void ReadCoefficientsAllY();
		void ReadCoefficients(string obsname);
		void ReadCoefficients(unsigned int iY);

	};

};

#endif
