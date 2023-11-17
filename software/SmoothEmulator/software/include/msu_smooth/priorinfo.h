#ifndef __PRIORINFO_H__
#define __PRIORINFO_H__

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/log.h"
//#include "msu_commonutils/constants.h"
//#include <list>
//#include <iostream>
//#include <Eigen/Dense>
//using namespace NMSUPratt;

namespace NBandSmooth{

	class CPriorInfo{
	public:
		CPriorInfo(string parinfo_filename);
		int NModelPars;
		string parinfo_filename;
		vector<string> parname,type; // type is gaussian or linear
		vector<double> xmin, xmax;
		map<string,int> name_map;
		int GetIPosition(string par_name);  // finds position given name of observable
		string GetName(int iposition);
		void PrintInfo();
	};

}

#endif
