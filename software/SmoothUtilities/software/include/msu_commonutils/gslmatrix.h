#ifndef __GSLMATRIX_H__
#define __GSLMATRIX_H__

#include "commondefs.h"

using namespace std;

// REAL MATRICES

class CGSLMatrix_Real{
public:
	int dim;
	// solve y=A*x for x
	void SolveLinearEqs(double *y,double **A,double *x);
	void SolveLinearEqs(vector<double> &y,vector<vector<double>> &A,vector<double> &x);
	void EigenFind(double **A,double **eigenvec,double *eigenval); // A must be symmetric
	void Invert(double **A,double **Ainv); // A must be symmetric
	void Cholesky_Invert(double **A,double **Ainv);
	void Invert_NonSymm(double **A,double **Ainv);
	void Print(double **A);
	double Determinant(vector<vector<double>> &A);
	CGSLMatrix_Real(int dimset);
	~CGSLMatrix_Real();
private:
	gsl_eigen_symmv_workspace *w;
	gsl_vector *eval;
	gsl_matrix *evec;
	gsl_vector *g;
	gsl_permutation *p;
	gsl_matrix *m,*minv;
	gsl_vector *v;
	double **U;
};

// COMPLEX MATRICES
  
class CGSLMatrix_Complex{
public:
  int dim;
  void SolveLinearEqs(complex<double> *y,complex<double> **A,complex<double> *x);
  void EigenFind(complex<double> **A,complex<double> **eigenvec,double *eigenval); // A must be Hermittian
  void Invert(complex<double> **A,complex<double> **Ainv); // A must  be Herm.
  CGSLMatrix_Complex(int dimset);
  ~CGSLMatrix_Complex();
private:
  gsl_vector *eval;
  gsl_matrix_complex *evec;
  gsl_vector_complex *g;
  gsl_eigen_hermv_workspace *w;
  gsl_permutation *p;
  complex<double> **U;
  gsl_matrix_complex *m;
  gsl_vector_complex *v;
};

#endif
