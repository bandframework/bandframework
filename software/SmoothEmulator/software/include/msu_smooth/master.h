#ifndef __SMOOTH_MASTER_H__
#define __SMOOTH_MASTER_H__
#include <cstdlib>
#include <cmath>

#include <cstdio>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <sstream>
#include <iomanip>
#include "msu_smoothutils/parametermap.h"
#include "msu_smoothutils/misc.h"
#include "msu_smoothutils/randy.h"
#include "msu_smooth/emulator.h"
#include "msu_smooth/modelparinfo.h"
#include "msu_smooth/smooth.h"
#include "msu_smoothutils/log.h"
#include "msu_smooth/observableinfo.h"
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/traininginfo.h"
#include "msu_smooth/testinginfo.h"

using namespace NMSUUtils;

namespace NBandSmooth{
   class CSmoothEmulator;
   class CTrainingInfo;
   class CTestingInfo;
   class CPriorInfo;
   class CObservableInfo;
   class CModelParInfo;
   
   class CSmoothMaster{
   public:
      CSmoothMaster();
      CparameterMap *parmap;
      unsigned int NPars;
      vector<CSmoothEmulator *> emulator;
      CTrainingInfo *traininginfo;
      CTestingInfo *testinginfo;
      CObservableInfo *observableinfo;
      CPriorInfo *priorinfo;
      string FullModelRunsDirName,FullModelTestingRunsDirName;
      Crandy *randy;
      string SmoothEmulator_TrainingFormat,SmoothEmulator_TestingFormat;
      double fitpercentage;
      CModelParameters *modelpars;

      int GetNPars(){
         return NPars;
      }
      int GetNObs(){
         return observableinfo->NObservables;
      }
      
      void ReadTrainingInfo();
      void ReadTestingInfo();
      //void GenerateCoefficientSamples();
      void CalcAllSigmaALambda();
      void TuneAllY(); // tune all observables
      void TuneY(string obsname); // tune one observable
      void TuneY(unsigned int iY); // tune one observable
      void TuneAllYFixedLambda(); // tune all observables with fixed Lambda
      void TuneY(string obsname,double LAMBDA); // tune one observable with fixed Lambda
      void TuneY(unsigned int iY,double LAMBDA); // tune one observable with fixed Lambda
      
      void GetAllY(CModelParameters *modelpars,vector<double> &Y,vector<double> &SigmaY_emulator);
      void GetAllYFromTheta(vector<double> &theta,vector<double> &Y,vector<double> &SigmaY_emulator);
      void GetAllYFromX(vector<double> &X,vector<double> &Y,vector<double> &SigmaY_emulator);
      void GetAllYOnly(CModelParameters *modelpars,vector<double> &Y);
      void GetAllYOnlyFromTheta(vector<double> &theta,vector<double> &Y);
      void GetAllYOnlyFromX(vector<double> &theta,vector<double> &X);
      
      
      void GetY(unsigned int iY,CModelParameters *modelpars,double &Y,double &SigmaY_emulator);
      void GetYFromTheta(unsigned int iY,vector<double> &theta,double &Y,double &SigmaY_emulator);
      void GetYFromX(unsigned int iY,vector<double> &X,double &Y,double &SigmaY_emulator);
      void GetY(string obsname,CModelParameters *modelpars,double &Y,double &SigmaY_emulator);
      void GetYFromTheta(string obsname,vector<double> &theta,double &Y,double &SigmaY_emulator);
      void GetYFromX(string obsname,vector<double> &X,double &Y,double &SigmaY_emulator);
      
      double GetYOnly(unsigned int iY,CModelParameters *modelpars);
      double GetYOnlyFromTheta(unsigned int iY,vector<double> &theta);
      double GetYOnlyFromX(unsigned int iY,vector<double> &X);
      double GetYOnly(string obsname,CModelParameters *modelpars);
      double GetYOnlyFromTheta(string obsname,vector<double> &theta);
      double GetYOnlyFromX(string obsname,vector<double> &theta);
      double GetYOnlyFromTheta(int iY,vector<double> theta);
      double GetYOnlyFromX(int iY,vector<double> theta);
      
      double GetYOnlyFromThetaPython(int DiY,vector<double> theta);
      double GetYOnlyFromXPython(int DiY,vector<double> X);
      vector<double> GetYSigmaFromThetaPython(int DiY,vector<double> theta);
      vector<double> GetYSigmaFromXPython(int DiY,vector<double> theta);
      
      double GetUncertaintyFromTheta(string obsname,vector<double> &Theta);
      double GetUncertaintyFromX(string obsname,vector<double> &X);
      double GetUncertaintyFromTheta(unsigned int iY,vector<double> &Theta);
      double GetUncertaintyFromX(unsigned int iY,vector<double> &X);
      
      void TestAtTrainingPts();
      void TestAtTrainingPts(string obsname);
      void TestAtTrainingPts(unsigned int iY);
      void TestVsFullModel();

      vector<double> GetXFromTheta(vector<double> Theta);
      vector<double> GetThetaFromX(vector<double> X);
      
      void ReadSigmaLambda();
      void ReadSigmaLambda(string filename);
      void WriteSigmaLambda();
      void WriteSigmaLambda(string filename);
      
   };
   
};

#endif
