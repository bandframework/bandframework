#include "msu_commonutils/qualifier.h"
#include "msu_commonutils/parametermap.h"

using namespace std;

void CQualifiers::Read(string qfilename){
	CQualifier *qptr;
	char cread[120],dummy[120];
	int npars=0,iqual=-1;
	string sread;
	FILE *fptr=fopen(qfilename.c_str(),"r");
	while(!feof(fptr)){
		fscanf(fptr,"%s",cread);
		sread=cread;
		if(sread[0]=='#'){
			fgets(dummy,120,fptr);
		}
		else if(sread=="qualifier"){
			npars=0;
			iqual+=1;
			qptr=new CQualifier();
			qualifier.push_back(qptr);
			fscanf(fptr,"%s",cread);
			sread=string(cread);
			qualifier[iqual]->qualname=sread;
			fgets(dummy,120,fptr);
		}
		else{
			if(sread=="int" || sread=="double" || sread=="bool"){
				qualifier[iqual]->type.push_back(sread);
				fscanf(fptr,"%s",cread);
				sread=string(cread);
				qualifier[iqual]->parname.push_back(sread);
				fscanf(fptr,"%s",cread);
				sread=string(cread);
				qualifier[iqual]->value.push_back(sread);
			}
			else{	
				qualifier[iqual]->type.push_back("unknown");
				qualifier[iqual]->parname.push_back(sread);
				fscanf(fptr,"%s",cread);
				sread=string(cread);
				qualifier[iqual]->value.push_back(sread);
			}
			fgets(dummy,120,fptr);
			npars+=1;
		}
		if(iqual>=0)
			qualifier[iqual]->npars=npars;
	}
	nqualifiers=qualifier.size();
	fclose(fptr);
}

void CQualifiers::SetPars(CparameterMap *pmap,int iqual){
	for(int ipar=0;ipar<qualifier[iqual]->npars;ipar++){
		pmap->set(qualifier[iqual]->parname[ipar],qualifier[iqual]->value[ipar]);
	}
}

void CQualifiers::Print(){
	printf("------ QUALIFIERS -------\n");
	for(int iqual=0;iqual<nqualifiers;iqual++){
		for(int ipar=0;ipar<qualifier[iqual]->npars;ipar++){
			printf("%s %s %s\n",qualifier[iqual]->type[ipar].c_str(),qualifier[iqual]->parname[ipar].c_str(),qualifier[iqual]->value[ipar].c_str());
		}
	}
}
