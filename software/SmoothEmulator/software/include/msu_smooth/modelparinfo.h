#ifndef __MODPARINFO_H__
#define __MODPARINFO_H__

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/log.h"
#include "msu_commonutils/constants.h"
#include "msu_smooth/priorinfo.h"
#include <Eigen/Dense>

class CModelParameters{
public:
//has actual values of the parameter
	vector<double> X;
	vector<double> Theta;
	CModelParameters(CPriorInfo *priorinfo_set);
	void TranslateTheta_to_X();
	void TranslateX_to_Theta();
	void Print();
	int NModelPars;
	CPriorInfo *priorinfo;
	void SetX(vector<double> &xset);
};

#endif
