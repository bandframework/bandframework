#ifndef __RANDY_H__
#define __RANDY_H__

#include "commondefs.h"

using namespace std;

class Crandy{
public:
	Crandy(int iseed);
	int seed;
	double threshold,netprob;
	double ran();
	double ran_gauss(); // return x propto exp(-x^2/2), <x^2> = 1
	void ran_gauss2(double &g1,double &g2);
	double ran_exp(); // return x propto exp(-x)
	double ran_lorentzian(); // return x propto 1/(1+4x^2), i.e. multiply by full width
	double ran_invcosh();    // return x propto 1/cosh(x), <x^2> = PI^2/4.
	void reset(int iseed);
	void generate_boltzmann(double mass,double T,FourVector &p);
	void generate_boltzmann_alt(double mass,double T,FourVector &p);
	int poisson();
	void set_mean(double mu);  // For Poisson Dist

	void increase_threshold();
	void increment_netprob(double delN);
	bool test_threshold(double delprob);

	
	//private:
	//std::mt19937 mt;   // 32 bit
	std::mt19937_64 mt;   // 64 bit
	std::uniform_real_distribution<double> ranu;
	std::normal_distribution<double> rang;
	std::poisson_distribution<int> ranp;
};

#endif