#ifndef __QUALIFIER_H__
#define __QUALIFIER_H__
#include "commondefs.h"

using namespace std;

class CQualifier{
public:
	int npars;
	string qualname;
	vector<string> type;
	vector<string> parname;
	vector<string> value;	
};

class CQualifiers{
public:
	vector<CQualifier *> qualifier;
	int nqualifiers;
	void Read(string qfilename);
	void SetPars(CparameterMap *pmap,int iqualifier);
	void Print();
};


#endif
