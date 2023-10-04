#ifndef __COMMON_DEFS_H__
#define __COMMON_DEFS_H__

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <complex>
#include <sys/stat.h>
#include <ctime>
#include <vector>
#include <array>
#include <unordered_map>
#include <map>
#include <cmath>
#include <random>
#include <vector>
#include <map>
#include <list>
#include <unordered_map>
#include <Eigen/Dense>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

using namespace std;

class CCHCalc;
class CCHArray;
class CHydroBalance;
class CHydroMesh;
class CHyperMesh;
class CHBEoS;
class CHBCharge;
class CdecayInfo;
class CresList;
class CresInfo;
class CbranchInfo;
class CHyperElement;
class CBalance;
class CAcceptance;
class CAction;
class CB3D;
class CB3DCell;
class Crandy;
class CparameterMap;
class CBalanceArrays;
class Csampler;
class CLocalInfo;
class CAction;
class CRegenerate;
class CSEInfo;
class CHYDROtoB3D;
class CMCList;
class CKernel;
class C3DArray;
class CWaveFunction;
class CSourceCalc;
class CKernelWF;
class CMSUPart;
class CdecayInfo;

typedef multimap<long int,CresInfo *> CresInfoMap;
typedef multimap<double,CresInfo *> CresMassMap;
typedef unordered_map<int,CdecayInfo *> CdecayInfoMap;
typedef pair<long int, CresInfo*> CresInfoPair;
typedef pair<long int, CdecayInfo*> CdecayInfoPair;
typedef pair<double,CresInfo*> CresMassPair;
typedef vector<CbranchInfo *> CbranchList; //gives branchlist name
typedef multimap<int,CHBCharge* > CHBChargeMap;
typedef pair<int,CHBCharge* > CHBChargePair;
//typedef multimap<int,CPart* > CPartMap;
//typedef pair<int,CPart* > CPartPair;
typedef multimap<int,CHBCharge* > mapic;
//typedef multimap<int,CPart* > mapip;
typedef array<double,4> FourVector;
//typedef double FourVector[4];
typedef double FourTensor[4][4];
typedef multimap<double,CAction *> CActionMap;
typedef pair<double,CAction*> CActionPair;
typedef pair<int,CHBCharge* > pairic;
typedef pair<int,CHBCharge* > pairip;

#endif
