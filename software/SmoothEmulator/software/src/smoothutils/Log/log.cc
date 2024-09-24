#include "msu_smoothutils/log.h"
using namespace NMSUUtils;

string CLog::logfilename="log.txt";
bool CLog::INTERACTIVE=true;
FILE *CLog::fptr=NULL;

using namespace std;

void CLog::Init(string &logfilename_in){
	INTERACTIVE=true;
	logfilename=logfilename_in;
	if(logfilename!="INTERACTIVE"){
		INTERACTIVE=false;
		fptr=fopen(logfilename.c_str(),"w");
	}
}

void CLog::Init(char *logfilename_in){
	INTERACTIVE=false;
	logfilename=logfilename_in;
	fptr=fopen(logfilename.c_str(),"w");
}

void CLog::Fatal(string message){
	if(INTERACTIVE){
		printf("FATAL error:\n");
		printf("%s",message.c_str());
	}
	else{
		fprintf(fptr,"FATAL error:\n");
		fprintf(fptr,"%s",message.c_str());
		fflush(fptr);
	}
	exit(1);
}

void CLog::Info(string message){
	if(INTERACTIVE){
		printf("%s",message.c_str());
	}
	else{
		fprintf(fptr,"%s",message.c_str());
		fflush(fptr);
	}
}

