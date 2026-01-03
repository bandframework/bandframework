#ifndef __TRAININGINFO_H__
#define __TRAININGINFO_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include <filesystem>

#include "msu_smoothutils/commonutils.h"
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

	class CTrainingInfo{
	public:
      CTrainingInfo(CObservableInfo *observableinfo,CPriorInfo *priorinfo);
      friend class CSmoothMaster;
      friend class CSmoothEmulator;
   private:
		CPriorInfo *priorinfo;
      CObservableInfo *observableinfo;
		unsigned int NTrainingPts,NObservables;
		vector<unsigned int> NTrainingList;
		vector<vector<double>> YTrain,SigmaYTrain;
		vector<CModelParameters *> modelpars;
		void ReadTrainingInfo();
		static CSmoothMaster *smoothmaster;
	};

};

#endif
