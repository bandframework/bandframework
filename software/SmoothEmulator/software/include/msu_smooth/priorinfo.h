#ifndef __PRIORINFO_H__
#define __PRIORINFO_H__

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "msu_smoothutils/parametermap.h"
#include "msu_smoothutils/misc.h"
#include "msu_smoothutils/randy.h"
#include "msu_smoothutils/log.h"
//#include "msu_smoothutils/constants.h"
//#include <list>
//#include <iostream>
//#include <Eigen/Dense>
//using namespace NMSUUtils;

namespace NBandSmooth{

	class CPriorInfo{
	public:
		CPriorInfo(string parinfo_filename);
		unsigned int NModelPars;
		string parinfo_filename;
		vector<string> parname,type; // type is gaussian or linear
		vector<double> xmin, xmax;
		map<string,unsigned int> name_map;
		unsigned int GetIPosition(string par_name);  // finds position given name of parameter
		string GetName(unsigned int iposition);
		void PrintInfo();
	};

};

#endif
