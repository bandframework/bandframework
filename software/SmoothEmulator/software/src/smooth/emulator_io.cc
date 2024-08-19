#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUUtils;

void CSmoothEmulator::PrintA(vector<double> &Aprint){
	if(!pca_ignore){
		for(unsigned int ic=0;ic<smooth->NCoefficients;ic++){
			CLog::Info(to_string(ic)+": "+to_string(Aprint[ic])+"\n");
		}
	}
}

void CSmoothEmulator::WriteCoefficients(){
	if(!pca_ignore){
		unsigned int NCoefficients=smooth->NCoefficients;
		unsigned int isample,ic,a;
		FILE *fptr;
		string filename;
		string dirname="smooth_data/"+smoothmaster->CoefficientsDirName+"/"+observable_name;
		string command="mkdir -p "+dirname;
		system(command.c_str());
		filename=dirname+"/meta.txt";
		fptr=fopen(filename.c_str(),"w");
		fprintf(fptr,"# NPars  MaxRank NCoefficients\n");
		fprintf(fptr,"%u  %u  %u\n",NPars,smooth->MaxRank,smooth->NCoefficients);
		fclose(fptr);
	
		filename=dirname+"/ABest.bin";
		fptr=fopen(filename.c_str(),"wb");
		for(ic=0;ic<smooth->NCoefficients;ic++){
			fwrite(&(ABest[ic]),sizeof(double),1,fptr);
		}
		fclose(fptr);
	
	}
}

void CSmoothEmulator::ReadCoefficients(){
	if(!pca_ignore){
		unsigned int isample,ic,a,NPars_test,MaxRank_test,NC_test;
		unsigned int NCoefficients=smooth->NCoefficients;
		FILE *fptr;
		string filename;
		string dirname="smooth_data/"+smoothmaster->CoefficientsDirName+"/"+observable_name;
		char dummy[100];
		filename=dirname+"/meta.txt";
		fptr=fopen(filename.c_str(),"r");
		fgets(dummy,100,fptr);
		fscanf(fptr,"%u  %u  %u\n",&NPars_test,&MaxRank_test,&NC_test);
		if(NPars_test!=NPars || MaxRank_test!=smooth->MaxRank || NC_test!=smooth->NCoefficients){
			CLog::Fatal("Mismatch in array sizes in ReadCoefficients");
		}
		fclose(fptr);
	
		ABest.resize(NCoefficients);
		filename = dirname+"/ABest.bin";
		fptr     = fopen(filename.c_str(),"rb");
		for(ic   = 0;ic<NCoefficients;ic++){
			fread(&(ABest[ic]),sizeof(double),1,fptr);
		}
		fclose(fptr);
		
		GetSigmaA();
		
	}

}
