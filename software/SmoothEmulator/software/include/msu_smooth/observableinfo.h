#ifndef __OBSERVABLE_INFO_H__
#define __OBSERVABLE_INFO_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "msu_commonutils/misc.h"
#include <vector>
#include <map>

namespace NBandSmooth{

	class CObservableInfo{
	public:
		CObservableInfo(string filename);
		int NObservables;
		vector<string> observable_name;
		vector<double> SigmaA0; // representative spread of coefficients
		map<string,int> name_map;
		int GetIPosition(string obsname);  // finds position given name of observable
		string GetName(int iposition);  // finds name give position
		void ReadObservableInfo(string filename);
		void ReadExperimentalInfo(string filename);
		vector<double> YExp,SigmaExp;
		void PrintInfo();
	};

}

#endif
