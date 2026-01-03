#ifndef __MODPARINFO_H__
#define __MODPARINFO_H__

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include "msu_smoothutils/parametermap.h"
#include "msu_smoothutils/misc.h"
#include "msu_smoothutils/randy.h"
#include "msu_smoothutils/log.h"
#include "msu_smooth/priorinfo.h"
#include <Eigen/Dense>
using namespace NMSUUtils;

namespace NBandSmooth{

	class CModelParameters{
	public:
		//has actual values of the parameter
		CModelParameters();
		void TranslateTheta_to_X();
		void TranslateX_to_Theta();
		void Print();
		void Write(string filename);
		void Copy(CModelParameters *mp);
		void SetX(vector<double> &xset);
		void SetTheta(vector<double> &thetaset);
      vector<double> X;
      vector<double> Theta;
      static CPriorInfo *priorinfo;
      static unsigned int NModelPars;
	};

};

#endif
