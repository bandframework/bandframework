#ifndef __INCLUDE_MISC_H_
#define __INCLUDE_MISC_H_
#include "msu_smoothutils/commondefs.h"

using namespace std;
namespace NMSUUtils{

	namespace Misc{

		bool comparestrings(char *s1,char *s2);
		
		void Cubic(double a0,double a1,double a2,double a3,
		complex<double>& z1,complex<double>& z2,complex<double>& z3);

		// Quartic routines provided by Andrew Steiner
		double signswitch(double a, double b);
		void Quartic(double a0,double a1,double a2,double a3,double a4,complex<double> &z1, complex<double> &z2,complex<double> &z3,complex<double> &z4);
		void Quartic(double a0,double a1,double a2,double a3,double a4,complex<double> *z);
		// this corresponds to z0=x[0], z1=x[1]+ix[2], z2=x[1]-ix[2]
		void CubicResolvant(double r,double s,double t,double x[],double &d);

		int CubicReal(double a0,double a1,double a2,double a3,double *x);
		void CubicComplex(double a0,double a1,double a2,double a3,complex<double> &z1,complex<double> &z2,complex<double> &z3);
		void Pause();
		void Pause(int seconds);
		bool file_exists( const char* file_name );
		bool file_exists( const std::string& file_name );
		double mod(double x, double y);

	};
}

#endif
