#include "msu_smooth/master.h"
using namespace std;
using namespace NBandSmooth;

CSmoothMaster::CSmoothMaster(){
   unsigned int NObs;
   parmap=new CparameterMap;
   parmap->ReadParsFromFile("smooth_data/Options/emulator_options.txt");
   
   string logfilename=parmap->getS("SmoothEmulator_LogFileName","Screen");
   if(logfilename!="Screen"){
      CLog::Init(logfilename);
   }
   string filename;
   filename="smooth_data/Info/observable_info.txt";
   
   observableinfo=new CObservableInfo(filename);
   NObs=observableinfo->NObservables;
   
   FullModelRunsDirName=parmap->getS("Smooth_FullModelRunsDirName","smooth_data/FullModelRuns");
   FullModelTestingRunsDirName=parmap->getS("Smooth_FullModelTestingRunsDirName","smooth_data/FullModelTestingRuns");
   
   filename="smooth_data/Info/prior_info.txt";
   priorinfo=new CPriorInfo(filename);
   NPars=priorinfo->NModelPars;
   parmap->set("SmoothEmulator_NPars",NPars);
   parmap->set("Smooth_NPars",NPars);
   CTrainingInfo::smoothmaster=this;
   CTestingInfo::smoothmaster=this;
   traininginfo = new CTrainingInfo(observableinfo,priorinfo);
   testinginfo = new CTestingInfo(observableinfo,priorinfo);
   CSmoothEmulator::NPars=NPars;
   CSmoothEmulator::smoothmaster=this;
   CSmoothEmulator::parmap=parmap;
   emulator.resize(NObs);
   
   for(unsigned int iy=0;iy<NObs;iy++){
      emulator[iy]=new CSmoothEmulator(observableinfo->observable_name[iy]);
   }
   ReadTrainingInfo();
   modelpars=new CModelParameters();
   
}

void CSmoothMaster::CalcAllSigmaALambda(){
   for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
      emulator[iY]->CalcSigmaALambda();
   }
}

void CSmoothMaster::TuneAllY(){
   for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
      emulator[iY]->Tune();
   }
   WriteSigmaLambda();
}

void CSmoothMaster::TuneAllYFixedLambda(){
   double Lambda;
   for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
      Lambda=emulator[iY]->LAMBDA;
      emulator[iY]->Tune(Lambda);
   }
}

void CSmoothMaster::TuneY(string obsname){
   unsigned int iY=observableinfo->GetIPosition(obsname);
   emulator[iY]->Tune();
}

void CSmoothMaster::TuneY(string obsname,double LambdaSet){
   unsigned int iY=observableinfo->GetIPosition(obsname);
   emulator[iY]->Tune(LambdaSet);
}

void CSmoothMaster::TuneY(unsigned int iY){
   emulator[iY]->Tune();
}

void CSmoothMaster::TuneY(unsigned int iY,double LambdaSet){
   emulator[iY]->Tune(LambdaSet);
}

void CSmoothMaster::GetY(unsigned int iY,CModelParameters *modelpars,double &Y,double &SigmaY_emulator){
   emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,SigmaY_emulator);
}

void CSmoothMaster::GetYFromTheta(unsigned int iY,vector<double> &theta,double &Y,double &SigmaY_emulator){
   emulator[iY]->GetYAndUncertaintyFromTheta(theta,Y,SigmaY_emulator);
}

void CSmoothMaster::GetYFromX(unsigned int iY,vector<double> &X,double &Y,double &SigmaY_emulator){
   modelpars->SetX(X);
   emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,SigmaY_emulator);
}

void CSmoothMaster::GetY(string obsname,CModelParameters *modelpars,double &Y,double &SigmaY_emulator){
   unsigned int iY=observableinfo->GetIPosition(obsname);
   emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,SigmaY_emulator);
}

void CSmoothMaster::GetYFromTheta(string obsname,vector<double> &theta,double &Y,double &SigmaY_emulator){
   unsigned int iY=observableinfo->GetIPosition(obsname);
   emulator[iY]->GetYAndUncertaintyFromTheta(theta,Y,SigmaY_emulator);
}

void CSmoothMaster::GetYFromX(string obsname,vector<double> &X,double &Y,double &SigmaY_emulator){
   unsigned int iY=observableinfo->GetIPosition(obsname);
   modelpars->SetX(X);
   emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,SigmaY_emulator);
}

void CSmoothMaster::GetAllY(CModelParameters *modelpars,vector<double> &Y,vector<double> &SigmaY_emulator){
   unsigned int NObservables=observableinfo->NObservables;
   Y.resize(NObservables);
   SigmaY_emulator.resize(NObservables);
   for(unsigned int iY=0;iY<NObservables;iY++){
      GetY(iY,modelpars,Y[iY],SigmaY_emulator[iY]);
   }
}

void CSmoothMaster::GetAllYFromTheta(vector<double> &Theta,vector<double> &Y,vector<double> &SigmaY_emulator){
   unsigned int NObservables=observableinfo->NObservables;
   Y.resize(NObservables);
   SigmaY_emulator.resize(NObservables);
   for(unsigned int iY=0;iY<NObservables;iY++){
      emulator[iY]->GetYAndUncertaintyFromTheta(Theta,Y[iY],SigmaY_emulator[iY]);
   }
}

void CSmoothMaster::GetAllYFromX(vector<double> &X,vector<double> &Y,vector<double> &SigmaY_emulator){
   unsigned int NObservables=observableinfo->NObservables;
   modelpars->SetX(X);
   Y.resize(NObservables);
   SigmaY_emulator.resize(NObservables);
   for(unsigned int iY=0;iY<NObservables;iY++){
      emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y[iY],SigmaY_emulator[iY]);
   }
}

double CSmoothMaster::GetYOnly(string obsname,CModelParameters *modelpars){
   unsigned int iY=observableinfo->GetIPosition(obsname);
   double SigmaY_emulator,Y;
   emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,SigmaY_emulator);
   return Y;
   
}

double CSmoothMaster::GetYOnly(unsigned int iY,CModelParameters *modelpars){
   double SigmaY_emulator,Y;
   emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,SigmaY_emulator);
   return Y;
}

double CSmoothMaster::GetYOnlyFromTheta(string obsname,vector<double> &Theta){
   unsigned int iY=observableinfo->GetIPosition(obsname);
   double SigmaY_emulator,Y;
   emulator[iY]->GetYAndUncertaintyFromTheta(Theta,Y,SigmaY_emulator);
   return Y;
}

double CSmoothMaster::GetYOnlyFromX(string obsname,vector<double> &X){
   unsigned int iY=observableinfo->GetIPosition(obsname);
   double SigmaY_emulator,Y;
   modelpars->SetX(X);
   emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,SigmaY_emulator);
   return Y;
}

double CSmoothMaster::GetYOnlyFromTheta(unsigned int iY,vector<double> &Theta){
   double SigmaY_emulator,Y;
   emulator[iY]->GetYAndUncertaintyFromTheta(Theta,Y,SigmaY_emulator);
   return Y;
}

double CSmoothMaster::GetYOnlyFromX(unsigned int iY,vector<double> &X){
   double SigmaY_emulator,Y;
   modelpars->SetX(X);
   emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,SigmaY_emulator);
   return Y;
}

double CSmoothMaster::GetYOnlyFromTheta(int iY,vector<double> Theta){
   double SigmaY_emulator,Y;
   emulator[iY]->GetYAndUncertaintyFromTheta(Theta,Y,SigmaY_emulator);
   return Y;
}

double CSmoothMaster::GetYOnlyFromX(int iY,vector<double> X){
   double SigmaY_emulator,Y;
   modelpars->SetX(X);
   emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,SigmaY_emulator);
   return Y;
}

double CSmoothMaster::GetYOnlyFromThetaPython(int DiY,vector<double> Theta){
   double SigmaY_emulator,Y;
   unsigned int iY=DiY;
   if(iY<observableinfo->NObservables){
      emulator[iY]->GetYAndUncertaintyFromTheta(Theta,Y,SigmaY_emulator);
      return Y;
   }
   else
      return 0.0;
}

double CSmoothMaster::GetYOnlyFromXPython(int DiY,vector<double> X){
   double SigmaY_emulator,Y;
   unsigned int iY=DiY;
   modelpars->SetX(X);
   if(iY<observableinfo->NObservables){
      emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,SigmaY_emulator);
      return Y;
   }
   else
      return 0.0;
}

double CSmoothMaster::GetUncertaintyFromTheta(string obsname,vector<double> &Theta){
   double Y,Uncertainty;
   unsigned int iY=observableinfo->GetIPosition(obsname);
   emulator[iY]->GetYAndUncertaintyFromTheta(Theta,Y,Uncertainty);
   return Uncertainty;
}

double CSmoothMaster::GetUncertaintyFromX(string obsname,vector<double> &X){
   double Y,Uncertainty;
   modelpars->SetX(X);
   unsigned int iY=observableinfo->GetIPosition(obsname);
   emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,Uncertainty);
   return Uncertainty;
}

double CSmoothMaster::GetUncertaintyFromTheta(unsigned int iY,vector<double> &Theta){
   double Y,Uncertainty;
   emulator[iY]->GetYAndUncertaintyFromTheta(Theta,Y,Uncertainty);
   return Uncertainty;
}

double CSmoothMaster::GetUncertaintyFromX(unsigned int iY,vector<double> &X){
   double Y,Uncertainty;
   modelpars->SetX(X);
   emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,Uncertainty);
   return Uncertainty;
}

void CSmoothMaster::GetAllYOnly(CModelParameters *modelpars,vector<double> &Yvec){
   unsigned int NObservables=observableinfo->NObservables;
   Yvec.resize(NObservables);
   for(unsigned int iY=0;iY<NObservables;iY++){
      double Y,Uncertainty;
      emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,Uncertainty);
      Yvec[iY]=Y;
   }
}

void CSmoothMaster::GetAllYOnlyFromTheta(vector<double> &Theta,vector<double> &Yvec){
   unsigned int NObservables=observableinfo->NObservables;
   Yvec.resize(NObservables);
   for(unsigned int iY=0;iY<NObservables;iY++){
      double Y,Uncertainty;
      emulator[iY]->GetYAndUncertaintyFromTheta(Theta,Y,Uncertainty);
      Yvec[iY]=Y;
   }
}

void CSmoothMaster::GetAllYOnlyFromX(vector<double> &X,vector<double> &Yvec){
   unsigned int NObservables=observableinfo->NObservables;
   Yvec.resize(NObservables);
   modelpars->SetX(X);
   for(unsigned int iY=0;iY<NObservables;iY++){
      double Y,Uncertainty;
      emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,Uncertainty);
      Yvec[iY]=Y;
   }
}

vector<double> CSmoothMaster::GetYSigmaFromThetaPython(int DiY,vector<double> theta){
   unsigned int iY=DiY;
   double Y,SigmaY_emulator;
   modelpars->SetTheta(theta);
   if(iY<observableinfo->NObservables)
      emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,SigmaY_emulator);
   else{
      Y=SigmaY_emulator=0.0;
   }
   vector<double> YSigma;
   YSigma.resize(2);
   YSigma[0]=Y;
   YSigma[1]=SigmaY_emulator;
   return YSigma;
}

vector<double> CSmoothMaster::GetYSigmaFromXPython(int DiY,vector<double> X){
   unsigned int iY=DiY;
   double Y,SigmaY_emulator;
   modelpars->SetX(X);
   if(iY<observableinfo->NObservables)
      emulator[iY]->GetYAndUncertaintyFromTheta(modelpars->Theta,Y,SigmaY_emulator);
   else{
      Y=SigmaY_emulator=0.0;
   }
   vector<double> YSigma;
   YSigma.resize(2);
   YSigma[0]=Y;
   YSigma[1]=SigmaY_emulator;
   return YSigma;
}

/*
vector<double> CSmoothMaster::GetYSigmaPython(int DiY,vector<double> theta){
   unsigned int iY=DiY;
   double Y,SigmaY_emulator;
   if(iY>=0 && iY<observableinfo->NObservables)
      emulator[iY]->CalcYAndUncertainty(theta,Y,SigmaY_emulator);
   else{
      Y=SigmaY_emulator=0.0;
   }
   vector<double> YSigma;
   YSigma.resize(2);
   YSigma[0]=Y;
   YSigma[1]=SigmaY_emulator;
   return YSigma;
}
 */


vector<double> CSmoothMaster::GetXFromTheta(vector<double> Theta){
   vector<double> X;
   if(Theta.size()!=CModelParameters::NModelPars){
      CLog::Fatal("In GetXFromTheta, mismatch in vector size, NModelPars="+to_string(CModelParameters::NModelPars)+", Theta.size="+to_string(Theta.size()));
   }
   X.resize(Theta.size());
   modelpars->SetTheta(Theta);
   modelpars->TranslateTheta_to_X();
   X=modelpars->X;
   return X;
}


vector<double> CSmoothMaster::GetThetaFromX(vector<double> X){
   vector<double> Theta;
   if(X.size()!=CModelParameters::NModelPars){
      CLog::Fatal("In GetXFromTheta, mismatch in vector size, NModelPars="+to_string(CModelParameters::NModelPars)+", Theta.size="+to_string(X.size()));
   }
   Theta.resize(X.size());
   modelpars->SetX(X);
   modelpars->TranslateX_to_Theta();
   Theta=modelpars->Theta;
   return Theta;
   
}
