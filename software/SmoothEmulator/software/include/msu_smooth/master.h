#ifndef __SMOOTH_MASTER_H__
#define __SMOOTH_MASTER_H__
#include <cstdlib>
#include <cmath>

#include <cstdio>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <sstream>
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/constants.h"
#include "msu_smooth/emulator.h"
#include "msu_smooth/modelparinfo.h"
#include "msu_smooth/smooth.h"
#include "msu_commonutils/log.h"
#include "msu_smooth/observableinfo.h"
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/traininginfo.h"

class CSmoothEmulator;
class CTrainingInfo;
class CPriorInfo;
class CObservableInfo;
class CModelParInfo;

class CSmoothMaster{
public:
	int TrainType;
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

	void ReadTrainingInfo();
	void GenerateCoefficientSamples();
	void TuneAllY(); // tune all observables
	void TuneY(string obsname); // tune one observable
	void TuneY(int iY); // tune one observable
	void SetThetaTrain();
	void CalcY(int iY,CModelParameters *modelpars,double &Y,double &SigmaY_emulator);
	void CalcY(string obsname,CModelParameters *modelpars,double &Y,double &SigmaY_emulator);
	void CalcAllY(CModelParameters *modelpars,vector<double> &Y,vector<double> &SigmaY_emulator);
	void TestAtTrainingPts();
	void TestAtTrainingPts(string obsname);
	void TestAtTrainingPts(int iY);
	void WriteCoefficientsAllY();
	void WriteCoefficients(string obsname);
	void WriteCoefficients(int iY);
	void ReadCoefficientsAllY();
	void ReadCoefficients(string obsname);
	void ReadCoefficients(int iY);

};

#endif
