#ifndef __SCORECARD_H__
#define __SCORECARD_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include "msu_smoothutils/misc.h"
#include "msu_smoothutils/randy.h"
#include "msu_smooth/emulator.h"
#include "msu_smooth/smooth.h"
using namespace NMSUUtils;

namespace NBandSmooth{

	class CScoreCard{
	public:
		double score,YExp,SigmaYExp,SigmaYReal;
		unsigned int NTest;
		double yi,Pi,Pibar,Pi2bar;
		void CalcScore(CSmoothEmulator *emulator,vector<vector<double>> &ThetaTest,double YExpSet,double SigmaYExpSet);
	};

};

#endif
