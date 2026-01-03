#ifndef __testingINFO_H__
#define __testingINFO_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include <iomanip>
#include <filesystem>

#include "msu_smoothutils/parametermap.h"
#include "msu_smoothutils/misc.h"
#include "msu_smoothutils/randy.h"
#include <list>
#include "msu_smooth/smooth.h"
#include <iostream>
#include <Eigen/Dense>
#include "msu_smoothutils/log.h"
#include "msu_smooth/master.h"
#include "msu_smooth/modelparinfo.h"
#include "msu_smooth/observableinfo.h"
#include "msu_smooth/priorinfo.h"
using namespace NMSUUtils;

namespace NBandSmooth{
	class CSmoothMaster;

	class CTestingInfo{
	public:
		CObservableInfo *observableinfo;
		CPriorInfo *priorinfo;
		CTestingInfo(CObservableInfo *observableinfo,CPriorInfo *priorinfo);
		unsigned int NTestingPts,NObservables;
		vector<unsigned int> NtestingList;
		vector<vector<double>> YTest,SigmaYTest;
		vector<CModelParameters *> modelpars;
		void ReadTestingInfo();
		static CSmoothMaster *smoothmaster;
	};

};

#endif
