#include "msu_smooth/priorinfo.h"
#include "msu_smooth/modelparinfo.h"

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

CPriorInfo::CPriorInfo(string parinfo_filename_set){
	CModelParameters::priorinfo=this;
	parinfo_filename=parinfo_filename_set;
	double minval,maxval;
	char dummy1[200],dummy2[40],dummy3[200];

	FILE *fptr;
	fptr=fopen(parinfo_filename.c_str(), "r");
	NModelPars=0;
	do{
		fscanf(fptr, "%s",dummy1);
		if(dummy1[0]!='#' && !feof(fptr)){
			fscanf(fptr, "%s %lf %lf",dummy2,&minval,&maxval);
			if(!feof(fptr)){
				if(string(dummy2)!="uniform" && string(dummy2)!="gaussian"){
					CLog::Fatal("reading priorinfo: type="+string(dummy2)+". Must be uniform or gaussian.\n");
				}
				parname.push_back(string(dummy1));
				type.push_back(string(dummy2));
				xmin.push_back(minval);  // Gaussian type xmin refers to <x> and <xmax> is sigma
				xmax.push_back(maxval);
				name_map.insert(pair<string,unsigned int>(parname[NModelPars],NModelPars));
				NModelPars+=1;
			}
		}
		else{
			fgets(dummy3,200,fptr);
		}
	}while(!feof(fptr));
	CModelParameters::NModelPars=NModelPars;
	fclose(fptr);
}


unsigned int CPriorInfo::GetIPosition(string par_name){
	map<string,unsigned int>::iterator iter;
	pair<string,unsigned int> mpair;
	iter=name_map.find(par_name);
	if(iter==name_map.end()){
		CLog::Fatal("In CPriorInfo::GetIposition, cannot find parameter "+par_name+"\n");
	}
	return iter->second;
}

string CPriorInfo::GetName(unsigned int i){
	return parname[i];
}

void CPriorInfo::PrintInfo(){
	char cstring[CLog::CHARLENGTH];
	CLog::Info("Prior Info\n");
	CLog::Info("#         ParameterName Type   Xmin_or_Xbar  Xmax_or_SigmaX\n");
	for(unsigned int ipar=0;ipar<NModelPars;ipar++){
		snprintf(cstring,CLog::CHARLENGTH,"%2u: %15s %9s %10g %10g\n",
		ipar,GetName(ipar).c_str(),type[ipar].c_str(),xmin[ipar],xmax[ipar]);
		CLog::Info(cstring);
	}
}
