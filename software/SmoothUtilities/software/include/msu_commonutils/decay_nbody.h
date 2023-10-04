#ifndef __MSU_DECAY_NBODY_H__
#define __MSU_DECAY_NBODY_H__

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <string>
#include <cstring>
#include "constants.h"
#include "commondefs.h"
#include "misc.h"
#include "randy.h"

using namespace std;
using namespace Misc;

class CDecay_NBody{
public:
	CDecay_NBody(Crandy *randyset){
		randy=randyset;
		wmaxmax3=wmaxmax4=0.0;
		Ntry=Nsuccess=0;
	}
	double M,m1,m2,m3,m4,wmax3,wmax4,M12,M34,q12,q34,Qmag,maxfactor3,maxfactor4,wmaxmax3,wmaxmax4,wmaxmax;
	double wmax,maxfactor,KEtot,qbar;
	vector<double> Msum,masses,qsum;
	int nbodies;
	bool allmassless;
	long long int Ntry,Nsuccess;
	Crandy *randy;
	// These functions also set weights
	void SetMasses2(double M,double m1,double m2);
	void SetMasses3(double M,double m1,double m2,double m3);
	void SetMasses4(double M,double m1,double m2,double m3,double m4);
	void SetMasses(int nbodies,vector<double> &masses_set);
	void SetMasses_Trial(int nbodies,vector<double> &masses_set);
	//
	void GenerateMomenta2(FourVector &p1,FourVector &p2);
	void GenerateMomenta3(FourVector &p1,FourVector &p2,FourVector &p3);
	void GenerateMomenta4(FourVector &p1,FourVector &p2,FourVector &p3,FourVector &p4);
	void GenerateMomenta(vector<FourVector> &p);
	void GenerateMomentaAllMassless(vector<FourVector> &p);

private:
	// These calculate weights for specific M12 (M34) and divide by wmax. They also return q12 (q13) and Qmag
	// The need wmax3 or wmax4 to give correct weights.
	double GetW3();
	double GetW4();
	double GetW();
	//
	void ChooseM12();     // for 3-body decays
	void ChooseM12M34();  // for 4-body decays
	void ChooseMsum();   // for n-body decays

	// These functions generate directions of momenta in c.o.m. given q12 (and q13) and Qmag
	void CompleteMomenta3(FourVector &p1,FourVector &p2,FourVector &p3);
	void CompleteMomenta4(FourVector &p1,FourVector &p2,FourVector &p3,FourVector &p4);
	void CompleteMomenta(vector<FourVector> &p);

};

#endif