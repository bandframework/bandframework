#include "msu_smooth/master.h"
using namespace std;
using namespace NBandSmooth;

void CSmoothMaster::ReadTrainingInfo(){
	if(SmoothEmulator_TrainingFormat == "training_format_smooth"){
		traininginfo->ReadTrainingInfoSmoothFormat();
	}
	else if(SmoothEmulator_TrainingFormat == "training_format_surmise"){
		traininginfo->ReadTrainingInfoSurmiseFormat();
	}
}

void CSmoothMaster::WriteCoefficientsAllY(){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->WriteCoefficients();
	}
}

void CSmoothMaster::WriteCoefficients(string obsname){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	WriteCoefficients(iY);
}

void CSmoothMaster::WriteCoefficients(unsigned int iY){
	emulator[iY]->WriteCoefficients();
}

void CSmoothMaster::ReadCoefficientsAllY(){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->ReadCoefficients();
	}
}

void CSmoothMaster::ReadCoefficients(string obsname){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	ReadCoefficients(iY);
}

void CSmoothMaster::ReadCoefficients(unsigned int iY){
	emulator[iY]->ReadCoefficients();
}
