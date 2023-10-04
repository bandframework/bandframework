#include "msu_commonutils/gslmatrix.h"


CGSLMatrix_Real::CGSLMatrix_Real(int dimset){
	dim=dimset;
	w=NULL;
	eval=NULL;
	evec=NULL;
	g=NULL;
	p=NULL;
	U=NULL;
	m=NULL;
	minv=NULL;
	v=NULL;
}

CGSLMatrix_Real::~CGSLMatrix_Real(){
	if(w!=NULL) gsl_eigen_symmv_free(w);
	if(eval!=NULL) gsl_vector_free(eval);
	if(evec!=NULL) gsl_matrix_free(evec);
	if(g!=NULL) gsl_vector_free(g);
	if(p!=NULL) gsl_permutation_free(p);
	if(m!=NULL) gsl_matrix_free(m);
	if(minv!=NULL) gsl_matrix_free(minv);
	if(v!=NULL) gsl_vector_free(v);
	if(U!=NULL){
		for(int i=0;i<dim;i++)
			delete [] U[i];
		delete [] U;
	}
}

void CGSLMatrix_Real::SolveLinearEqs(double *y,double **A,double *x){
	int i,j,s;
	if(v==NULL)
		v=gsl_vector_alloc(dim);

	if(m==NULL)
		m=gsl_matrix_alloc(dim,dim);

	for(i=0;i<dim;i++){
		gsl_vector_set(v,i,y[i]);
		for(j=0;j<dim;j++)
			gsl_matrix_set(m,i,j,A[i][j]);
	}

	if(g==NULL) g = gsl_vector_alloc (dim);
	if(p==NULL) p = gsl_permutation_alloc (dim);

	gsl_linalg_LU_decomp (m,p,&s);
	gsl_linalg_LU_solve (m, p,v,g);

	for(i=0;i<dim;i++)
		x[i]=gsl_vector_get(g,i);
}

void CGSLMatrix_Real::SolveLinearEqs(vector<double> &y,vector<vector<double>> &A,vector<double> &x){
	int i,j,s;
	if(v==NULL)
		v=gsl_vector_alloc(dim);

	if(m==NULL)
		m=gsl_matrix_alloc(dim,dim);

	for(i=0;i<dim;i++){
		gsl_vector_set(v,i,y[i]);
		for(j=0;j<dim;j++)
			gsl_matrix_set(m,i,j,A[i][j]);
	}

	if(g==NULL) g = gsl_vector_alloc (dim);
	if(p==NULL) p = gsl_permutation_alloc (dim);

	gsl_linalg_LU_decomp (m,p,&s);
	gsl_linalg_LU_solve (m, p,v,g);

	for(i=0;i<dim;i++)
		x[i]=gsl_vector_get(g,i);
}

double CGSLMatrix_Real::Determinant(vector<vector<double>> &A){
	int i,j,s=0;
	if(m==NULL)
		m=gsl_matrix_alloc(dim,dim);
	if(p==NULL) p = gsl_permutation_alloc (dim);
	
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			gsl_matrix_set(m,i,j,A[i][j]);
			printf("m[%d][%d]=%g\n",i,j,gsl_matrix_get(m,i,j));
		}
	}
	
	printf("howdy, dim=%d\n",dim);
	gsl_linalg_LU_decomp (m, p, &s);
	printf("hmmmm\n");
	return gsl_linalg_LU_det (m, s);
}

void CGSLMatrix_Real::EigenFind(double **A,double **eigenvec,double *eigenval){
	int i,j;
	if(m==NULL) m=gsl_matrix_alloc(dim,dim);
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++)
			gsl_matrix_set(m,i,j,A[i][j]);
	}

	if(eval==NULL) eval=gsl_vector_alloc(dim);
	if(evec==NULL) evec=gsl_matrix_alloc(dim,dim);
	if(w==NULL) w=gsl_eigen_symmv_alloc(dim);

	gsl_eigen_symmv(m,eval,evec,w);
	//gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);

	for(i=0;i<dim;i++){
		eigenval[i]=gsl_vector_get(eval,i);
		for(j=0;j<dim;j++)
			eigenvec[i][j]=gsl_matrix_get(evec,i,j);
	}

}

void CGSLMatrix_Real::Invert(double **A,double **Ainv){
	int i,j;
	if(m==NULL) m=gsl_matrix_alloc(dim,dim);
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++) gsl_matrix_set(m,i,j,A[i][j]);
	}

if(eval==NULL) eval=gsl_vector_alloc(dim);
if(evec==NULL) evec=gsl_matrix_alloc(dim,dim);
if(w==NULL) w=gsl_eigen_symmv_alloc(dim);

gsl_eigen_symmv(m,eval,evec,w);
	//gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);

if(U==NULL){
	U=new double*[dim];
	for(i=0;i<dim;i++) U[i]=new double[dim];
}

for(i=0;i<dim;i++){
	for(j=0;j<dim;j++) {
		U[i][j]=gsl_matrix_get(evec,j,i);
		Ainv[i][j]=0.0;
	}
}  
int k;
for(i=0;i<dim;i++){
	for(j=0;j<dim;j++){
		for(k=0;k<dim;k++) Ainv[i][j]+=U[k][i]*U[k][j]/gsl_vector_get(eval,k);
	}
}
}

void CGSLMatrix_Real::Invert_NonSymm(double **A,double **Ainv){
	int i,j,signum;
	if(p==NULL) p = gsl_permutation_alloc (dim);
	if(m==NULL) m=gsl_matrix_alloc(dim,dim);
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++) gsl_matrix_set(m,i,j,A[i][j]);
	}

	if(minv==NULL) minv=gsl_matrix_alloc(dim,dim);
	gsl_linalg_LU_decomp(m, p, &signum);
	gsl_linalg_LU_invert(m, p, minv);

	int k;
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			for(k=0;k<dim;k++)
				Ainv[i][j]=gsl_matrix_get(minv,i,j);
		}
	}
	double **testm=new double *[dim];
	for(i=0;i<dim;i++)
		testm[i]=new double[dim];
	/**
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			for(k=0;k<dim;k++) testm[i][j]+=A[i][k]*Ainv[k][j];
			printf("%10.3e ",testm[i][j]);
		}
		printf("\n");
	}
	
	for(i=0;i<dim;i++) delete [] testm[i];
	delete [] testm;
	*/
}


void CGSLMatrix_Real::Print(double **A){
	int i,j;
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			printf("%11.4e ",A[i][j]);
		}
		printf("\n");
	}
}

// HERMITIAN COMPLEX MATRICES

CGSLMatrix_Complex::CGSLMatrix_Complex(int dimset){
	dim=dimset;
	eval=NULL;
	evec=NULL;
	w=NULL;
	g=NULL;
	p=NULL;
	U=NULL;
	m=NULL;
	v=NULL;
}

CGSLMatrix_Complex::~CGSLMatrix_Complex(){
	if(eval!=NULL) gsl_vector_free(eval);
	if(evec!=NULL) gsl_matrix_complex_free(evec);
	if(w!=NULL) gsl_eigen_hermv_free(w);
	if(g!=NULL) gsl_vector_complex_free(g);
	if(p!=NULL) gsl_permutation_free(p);
	if(m!=NULL) gsl_matrix_complex_free(m);
	if(v!=NULL) gsl_vector_complex_free(v);
	if(U!=NULL){
		for(int i=0;i<dim;i++)
			delete [] U[i];
		delete [] U;
	}
}

void CGSLMatrix_Complex::EigenFind(complex<double> **A,complex<double> **eigenvec,double *eigenval){
	complex<double> ci(0.0,1.0);
	gsl_complex z;
	//gsl_matrix_complex *m=gsl_matrix_complex_alloc(dim,dim);
	if(m==NULL) m=gsl_matrix_complex_alloc(dim,dim);
	int i,j;
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			GSL_SET_COMPLEX(&z,real(A[i][j]),imag(A[i][j]));
			gsl_matrix_complex_set(m,i,j,z);
		}
	}

	if(eval==NULL) eval=gsl_vector_alloc(dim);
	if(evec==NULL) evec=gsl_matrix_complex_alloc(dim,dim);
	if(w==NULL) w=gsl_eigen_hermv_alloc(dim);

	gsl_eigen_hermv(m,eval,evec,w);

	//gsl_eigen_hermv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);

	for(i=0;i<dim;i++){
		eigenval[i]=gsl_vector_get(eval,i);
		for(j=0;j<dim;j++){
			z=gsl_matrix_complex_get(evec,i,j);
			eigenvec[i][j]=GSL_REAL(z)+ci*GSL_IMAG(z);
		}
	}

}

void CGSLMatrix_Complex::Invert(complex<double> **A,complex<double> **Ainv){
	complex<double> ci(0.0,1.0);
	gsl_complex z;
	//gsl_matrix_complex *m=gsl_matrix_complex_alloc(dim,dim);
	if(m==NULL) m=gsl_matrix_complex_alloc(dim,dim);
	int i,j;
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			GSL_SET_COMPLEX(&z,real(A[i][j]),imag(A[i][j]));
			gsl_matrix_complex_set(m,i,j,z);
		}
	}

	if(eval==NULL) eval=gsl_vector_alloc(dim);
	if(evec==NULL) evec=gsl_matrix_complex_alloc(dim,dim);
	if(w==NULL) w=gsl_eigen_hermv_alloc(dim);

	gsl_eigen_hermv(m,eval,evec,w);
	//gsl_eigen_hermv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);

	if(U==NULL){
		U=new complex<double> *[dim];
		for(i=0;i<dim;i++) U[i]=new complex<double>[dim];
	}

for(i=0;i<dim;i++){
	for(j=0;j<dim;j++) {
		z=gsl_matrix_complex_get(evec,j,i);
		U[i][j]=GSL_REAL(z)+ci*GSL_IMAG(z);
		Ainv[i][j]=0.0;
	}
}
int k;
for(i=0;i<dim;i++){
	for(j=0;j<dim;j++){
		for(k=0;k<dim;k++)
			Ainv[i][j]+=U[k][i]*conj(U[k][j])/gsl_vector_get(eval,k);
	}
}
}

void CGSLMatrix_Complex::SolveLinearEqs(complex<double> *y,complex<double> **A,complex<double> *x){
	complex<double> ci(0.0,1.0);
	int i,j,s;
	gsl_complex z;

	if(m==NULL) m=gsl_matrix_complex_alloc(dim,dim);
	if(v==NULL) v=gsl_vector_complex_alloc(dim);
	for(i=0;i<dim;i++){
		GSL_SET_COMPLEX(&z,real(y[i]),imag(y[i]));
		gsl_vector_complex_set(v,i,z);
		for(j=0;j<dim;j++){
			GSL_SET_COMPLEX(&z,real(A[i][j]),imag(A[i][j]));
			gsl_matrix_complex_set(m,i,j,z);
		}
	}

	if(g==NULL) g = gsl_vector_complex_alloc (dim);
	if(p==NULL) p = gsl_permutation_alloc (dim);


	gsl_linalg_complex_LU_decomp (m, p, &s);
	gsl_linalg_complex_LU_solve (m, p, v, g);

	for(i=0;i<dim;i++){
		z=gsl_vector_complex_get(g,i);
		x[i]=GSL_REAL(z)+ci*GSL_IMAG(z);
	}
}


void CGSLMatrix_Real::Cholesky_Invert(double **A,double **Ainv){
	int cholesky_test,i,j;
	double determinant_c=0.0;
	if(m==NULL) m=gsl_matrix_alloc(dim,dim);
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			gsl_matrix_set(m,i,j,A[i][j]);
		}
	}
	gsl_error_handler_t *temp_handler; // disable the default handler
	// do a cholesky decomp of the cov matrix, LU is not stable for ill conditioned matrices
	temp_handler = gsl_set_error_handler_off();
	cholesky_test = gsl_linalg_cholesky_decomp(m);
	if(cholesky_test == GSL_EDOM){
		fprintf(stderr, "trying to cholesky a non postive def matrix, in emulate-fns.c sorry...\n");
		determinant_c = 1.0;
		for(i = 0; i < dim; i++)
			determinant_c *= gsl_matrix_get(m, i, i);
		printf("det CHOL:%g\n", determinant_c);	
		determinant_c = determinant_c * determinant_c;
		exit(1);
	}
	gsl_set_error_handler(temp_handler);
	// find the determinant and then invert 
	// the determinant is just the trace squared
	gsl_linalg_cholesky_invert(m);
	//gsl_matrix_memcpy(result_matrix,m);
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			Ainv[i][j]=gsl_matrix_get(m,i,j);
		}
	}
	/**
	double **mcheck=new double *[dim];
	for(i=0;i<dim;i++){
		mcheck[i]=new double[dim];
		for(j=0;j<dim;j++){
			mcheck[i][j]=0.0;
			for(k=0;k<dim;k++){
				mcheck[i][j]+=A[i][k]*Ainv[k][j];
			}
			printf("%6.4f ",mcheck[i][j]);
		}
		printf("\n");
		delete [] mcheck[i];
	}
	delete [] mcheck;
	**/
}

