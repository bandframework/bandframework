#include "msu_smooth/master.h"
using namespace std;
using namespace NBandSmooth;

void CSmoothMaster::ReadTrainingInfo(){
   traininginfo->ReadTrainingInfo();
}

void CSmoothMaster::ReadTestingInfo(){
   testinginfo->ReadTestingInfo();
}

void CSmoothMaster::WriteSigmaLambda(string filename){
   string command="mkdir -p smooth_data/output_stuff";
   system(command.c_str());
   FILE *fptr=fopen(filename.c_str(),"w");
   for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
      fprintf(fptr,"%24s %10.3f %10.5f\n",
              observableinfo->observable_name[iY].c_str(),emulator[iY]->SigmaA,emulator[iY]->LAMBDA);
   }
   fclose(fptr);
}

void CSmoothMaster::WriteSigmaLambda(){
   string filename="smooth_data/output_stuff/sigmalambda.txt";
   WriteSigmaLambda(filename);
}

void CSmoothMaster::ReadSigmaLambda(){
   string filename="smooth_data/output_stuff/sigmalambda.txt";
   ReadSigmaLambda(filename);
}

void CSmoothMaster::ReadSigmaLambda(string filename){
   FILE *fptr=fopen(filename.c_str(),"r");
   char obsname[200];
   string obsstring;
   for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
      fscanf(fptr,"%s %lf %lf",obsname,&emulator[iY]->SigmaA,&emulator[iY]->LAMBDA);
      obsstring=obsname;
      if(obsstring!=observableinfo->observable_name[iY]){
         CLog::Fatal("Reading in obs name ("+obsstring+") from "+filename+" - does not match name in observable info ("+observableinfo->observable_name[iY]+")");
      }
   }
   fclose(fptr);
}
