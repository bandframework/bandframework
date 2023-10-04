#include "msu_commonutils/sf.h"

using namespace std;

complex<double> CoulWave::cgamma(complex<double> z){
	complex<double> ci(0.0,1.0);
	gsl_sf_result gsl_lnr,gsl_arg;
	gsl_sf_lngamma_complex_e(real(z),imag(z),&gsl_lnr,&gsl_arg);
	return exp(gsl_lnr.val)*(cos(gsl_arg.val)+ci*sin(gsl_arg.val));
}

/* ************************************* */
void CoulWave::GetFG(int L,double x,double eta,double *FL,double *GL){
	double expF,expG;
	complex<double> cg,ci(0.0,1.0);
	double *fc,*gc;
	fc=new double[L+1];
	gc=new double[L+1];
	cg=CoulWave::cgamma(double(L+1)+ci*eta);
	//double sigma=atan2(imag(cg),real(cg));

	// This calculates fc and gc arrays for indices L to L+k  
	gsl_sf_coulomb_wave_FG_array(0,L,eta,x,fc,gc,&expF,&expG);
	*FL=fc[L]*exp(expF);
	*GL=gc[L]*exp(expG);
	/*
#ifndef NO_GSLCOULWAVE_BUG
	if(x>2.0){
		phase=atan2(*GL,*FL);
		phase=phase+x-0.5*PI*(L+1)+sigma-eta*log(2.0*x);
		phase+=0.5*double(L*(L+1))/x;
		phase=phase-2.0*PI*floor(0.5*phase/PI);
		if(phase>PI) phase=phase-2.0*PI;
		if(fabs(phase)>2.0){
			*FL=-*FL;
			*GL=-*GL;
			//printf("L=%d, x=%6.3f, phase=%10.6f, sigma=%g\n",L,x,phase,sigma);
		}
	}
#endif */
	delete [] fc;
	delete [] gc;
}

// This calcules dCW/dx
void CoulWave::GetFGprime(int L,double x,double eta,double *FL,double *GL,double *FLprime,double *GLprime){
	double expF,expG;
	double *fc,*gc,*fcp,*gcp;
	int k=0;
	fc=new double[k+1];
	gc=new double[k+1];
	fcp=new double[k+1];
	gcp=new double[k+1];
	// This calculates fc and gc arrays for indices L to L+k  
	gsl_sf_coulomb_wave_FGp_array(L,k,eta,x,fc,fcp,gc,gcp,
		&expF,&expG);
	*FL=fc[0]*exp(expF);
	*GL=gc[0]*exp(expG);
	*FLprime=fcp[0]*exp(expF);
	*GLprime=gcp[0]*exp(expG);
	delete [] fc;
	delete [] gc;
	delete [] fcp;
	delete [] gcp;
}

void CoulWave::GetFGPrimeArray(int Lmax,double x,double eta,
vector<double> &FL,vector<double> &GL,vector<double> &FLprime,vector<double> &GLprime){
	double expF,expG;
	double *fc,*gc,*fcp,*gcp;
	int ell;
	fc=new double[Lmax+1];
	gc=new double[Lmax+1];
	fcp=new double[Lmax+1];
	gcp=new double[Lmax+1];
	// This calculates fc and gc arrays for indices L to L+k  
	gsl_sf_coulomb_wave_FGp_array(0,Lmax,eta,x,fc,fcp,gc,gcp,&expF,&expG);
	
	if(int(FL.size())!=Lmax+1){
		FL.resize(Lmax+1);
		GL.resize(Lmax+1);
		FLprime.resize(Lmax+1);
		GLprime.resize(Lmax+1);			
	}
	for(ell=0;ell<=Lmax;ell++){
		FL[ell]=fc[ell]*exp(expF);
		GL[ell]=gc[ell]*exp(expG);
		FLprime[ell]=fcp[ell]*exp(expF);
		GLprime[ell]=gcp[ell]*exp(expG);
	}
	delete [] fc;
	delete [] gc;
	delete [] fcp;
	delete [] gcp;
}

complex<double> CoulWave::CWoutgoing(int ell,double x,double eta){
	double FL,GL;
	complex<double> ci(0.0,1.0);
	GetFG(ell,x,eta,&FL,&GL);
	return FL-ci*GL;
} 

complex<double> CoulWave::CWincoming(int ell,double x,double eta){
	double FL,GL;
	complex<double> ci(0.0,1.0);
	GetFG(ell,x,eta,&FL,&GL);
	return FL+ci*GL;
}

/* ************************************* */

void CoulWave::phaseshift_CoulombCorrect(int ell,double q,double eta,
double &delta,double &ddeltadq){
	const double PI=4.0*atan(1.0);
	double *gamow,*dgamowdq;
	double tandelta0,tandelta,dtandelta0dq,dtandeltadq;
	int i;

	if(fabs(eta)>1.0E-10){
		gamow=new double[ell+1];
		dgamowdq=new double[ell+1];
		gamow[0]=2.0*PI*eta/(exp(2.0*PI*eta)-1.0);
		dgamowdq[0]=(gamow[0]/q)*(gamow[0]*exp(2.0*PI*eta)-1.0);
		for(i=0;i<ell;i++){
			gamow[i+1]=gamow[i]*(double((i+1)*(i+1))+eta*eta);
			dgamowdq[i+1]=dgamowdq[i]*(double((i+1)*(i+1))+eta*eta)
				-(2.0*eta*eta/q)*gamow[i];
		}

		tandelta0=tan(delta);
		dtandelta0dq=ddeltadq*(1.0+tandelta0*tandelta0);

		tandelta=tandelta0*gamow[ell];
		if(cos(delta)>0.0) delta=atan(tandelta);
		else delta=atan(tandelta)+PI;
		if(delta>PI) delta=delta-2.0*PI;

		dtandeltadq=gamow[ell]*dtandelta0dq+tandelta0*dgamowdq[ell];
		ddeltadq=dtandeltadq/(1.0+tandelta*tandelta);
		delete [] gamow;
		delete [] dgamowdq;
	}
}

void CoulWave::GetFGprime_ImagQ(int ellmax,double x,double eta,double *FL,double *GL,double *FLprime,double *GLprime){
	double *F,*G,*Fprime,*Gprime;
	F=new double[ellmax+1]; G=new double[ellmax+1]; Fprime=new double[ellmax+1]; Gprime=new double[ellmax+1];
	double ff,root,sign;
	int n,l,nmax;
	double *A,*B;
	nmax=24+lrint(fabs(x));
	A=new double[nmax+1];
	B=new double[nmax+1];
	// Calc F and Fprime
	A[0]=0.0; A[1]=1.0;
	for(n=0;n<=nmax;n++) B[n]=0.0;
	for(n=0;n<=nmax-2;n++) A[n+2]=(2.0*eta*A[n+1]+A[n])/double((n+1)*(n+2));
	F[0]=Fprime[0]=0.0;	
	for(n=1;n<=nmax;n++){
		F[0]+=A[n]*pow(x,n);
		Fprime[0]+=double(n)*A[n]*pow(x,n-1);
	}
	for(l=0;l<ellmax;l++){
		root=(l+1.0)*(l+1.0)-eta*eta;
		if(root>=0){
			root=sqrt(root);
			sign=1;
		}
		else{
			root=sqrt(-root);
			sign=-1;
		}
		//printf("root=%g\n",root);
		ff=(((l+1.0)*(l+1.0)/x)+eta);
		F[l+1]=(ff*F[l]-(l+1.0)*Fprime[l])/root;
		Fprime[l+1]=-(sign*root*F[l]+ff*F[l+1])/(l+1.0);
	}
	// CALC G and Gprime
	A[0]=1.0; A[1]=0.0; B[1]=2.0*eta*A[0]; B[0]=0.0;
	for(n=0;n<=nmax-2;n++){
		B[n+2]=(2.0*eta*B[n+1]+B[n])/((n+1.0)*(n+2.0));
		A[n+2]=(2.0*eta*A[n+1]+A[n]+(1.0-2.0*(n+2.0))*B[n+2])/((n+1.0)*(n+2.0));
	}
	G[0]=1; Gprime[0]=0.0;
	for(n=1;n<=nmax;n++){
		G[0]+=(A[n]+B[n]*log(fabs(x)))*pow(x,n);
		Gprime[0]+=double(n)*(A[n]+B[n]*log(fabs(x)))*pow(x,n-1)+B[n]*pow(x,n-1);
	}
	for(l=0;l<ellmax;l++){
		root=(l+1.0)*(l+1.0)-eta*eta;
		if(root>=0){
			root=sqrt(root);
			sign=1;
		}
		else{
			root=sqrt(-root);
			sign=-1;
		}
		ff=((l+1.0)*(l+1.0)/x)+eta;
		G[l+1]=(ff*G[l]-(l+1.0)*Gprime[l])/root;
		Gprime[l+1]=-(sign*root*G[l]+ff*G[l+1])/(l+1.0);
	}

	*FL=F[ellmax];
	*GL=G[ellmax];
	*FLprime=Fprime[ellmax];
	*GLprime=Gprime[ellmax];

	delete[] A;
	delete[] B;
	delete[] F;
	delete [] G;
	delete [] Fprime;
	delete [] Gprime;
}

void CoulWave::GetFGprime_ComplexQ(int ell,complex<double> cx,complex<double> ceta,double *FL,double *GL,double *FLprime,double *GLprime){
	// This works for q either real or imaginary, not really for arbitrarily complex q
	if(fabs(imag(cx)) > fabs(real(cx))){
		// Using the recurrence relation to get the derivatives of F and G
		// (L+1)du/dx = ((L+1)^2/x + eta)u - sqrt((L+1)^2+eta^2)uL+1
		if(abs(ceta)>1.0E-8) CoulWave::GetFGprime_ImagQ(ell,imag(cx),-imag(ceta),FL,GL,FLprime,GLprime);
		else Bessel::CalcJN_imag(ell,imag(cx),*FL,*GL,*FLprime,*GLprime);
	}
	else{
		if(abs(ceta)>1.0E-8) CoulWave::GetFGprime(ell,real(cx),real(ceta),FL,GL,FLprime,GLprime);
		else Bessel::CalcJN_real(ell,real(cx),*FL,*GL,*FLprime,*GLprime);
	}
	return;
}

complex<double> CoulWave::cw2_small_r(int l,complex<double> r,
complex<double>  eta){
	const double PI=4.0*atan(1.0);
	const double EULER=0.5772156649015328606;

	/* The notation is like Gr. + R, page 1063.
	The Coulomb wave function is the same as W(i*eta,l+1/2,2*i*rho) */
		double psi1,psi2,sum2;
	complex<double> delsum1,sum1,lterm,answer,ci(0.0,1.0);
	double cdcon1,cdcon2,factor,cx,delp,psi3,fact1,fact2,rr,ceta;
	int k;
	rr=real(-ci*r);
	ceta=real(-ci*eta);

	//printf("eta=(%g,%g)\n",real(eta),imag(eta));
	cdcon1=real(CoulWave::cgamma(-double(l)+ceta));
	cdcon2=real(CoulWave::cgamma(double(l+1)+ceta));

	factor=pow(-1,double(2*l+1))*pow(2.0*rr,double(l+1))
		*exp(rr)/(cdcon1*cdcon2);
	psi1=-EULER;
	psi2=-EULER;
	for(k=1;k<=2*l+1;k++){
		psi2=psi2+1.0/double(k);
	}
	cx=l+1+ceta;
	psi3=-EULER-(1.0/cx)+cx*(PI*PI/6.0);
	for(k=1;k<100000;k++){
		delp=-cx*cx/(double(k*k)*(cx+double(k)));
		psi3=psi3+delp;
		//     if(abs(delp)<1.0E-12) goto CONVERGE1;
		if(fabs(delp)<1.0E-12) goto CONVERGE1;
	}
	printf("never escaped loop1 in cw2_small_r!\n");
	CONVERGE1:
	lterm=ci*PI+log(2.0*rr);
	fact1=cdcon2/dgamma(2*l+2);
	sum1=fact1*(psi1+psi2-psi3-lterm);
	for(k=1;k<=10000;k++){
		fact1=fact1*(-2.0*rr)*(double(l+k)+ceta)
			/(double(k)*double(2*l+1+k));
		psi1=psi1+1.0/double(k);
		psi2=psi2+1.0/double(2*l+1+k);
		psi3=psi3+1.0/(double(k-1)+cx);
		delsum1=fact1*(psi1+psi2-psi3-lterm);
		sum1=sum1+delsum1;
		if(abs(delsum1)<1.0E-15) goto CONVERGE2;
	}
	printf("never escaped loop2 in cw2_small_r!\n");
	CONVERGE2:
	fact2=dgamma(2.0*l+1)*cdcon1/pow(2.0*rr,2*l+1);
	sum2=fact2;
	//printf("r=(%g,%g), cdcon1=(%g,%g)\n",real(r),imag(r),real(cdcon1),imag(cdcon1));
	//printf("k=0, fact2=(%g,%g)\n",real(fact2),imag(fact2));
	for(k=1;k<=2*l;k++){
		fact2=fact2*(double(k-l-1)+ceta)*rr/(double(k)*double(2*l-k+1));
		//printf("k=%d, fact2=(%g,%g)\n",k,real(fact2),imag(fact2));
		sum2=sum2+fact2;
	}
	sum1=factor*sum1;
	sum2=factor*sum2;
	answer=(sum1+sum2);//*exp(0.5*PI*eta);
	//printf("factor=(%g,%g)\n",real(factor),imag(factor));
	/*
	printf("r=(%g,%g), eta=(%g,%g), factor=(%g,%g), sum1=(%g,%g), sum2=(%g,%g), answer=(%g,%g)\n",real(r),imag(r),
		real(eta),imag(eta),real(factor),imag(factor),real(sum1),imag(sum1),real(sum2),imag(sum2),real(answer),imag(answer));
	*/
	return answer;
}
/* ****************************************** */

double CoulWave::dgamma(int mm){
	/* This calc.s gamma functions which are in the form gamma(n)
	where n is an int > 0. */
		double dg;
	int j;
	dg=1.0;
	if(mm<1) {
		for(j=1;j<=-mm+1;j++){
			dg=dg/(1.0+double(-j));
		}
	}
	if(mm>1){
		for(j=1;j<=mm-1;j++){
			dg=dg*double(j);
		}
	}
	return dg;
}

void CoulWave::SphericalCW(int L,double x,double eta,double *FL,double *GL){
	double expF,expG;
	double *fc,*gc;
	fc=new double[L+1];
	gc=new double[L+1];

	// This calculates fc and gc arrays for indices L to L+k
	gsl_sf_coulomb_wave_FG_array(0,L,eta,x,fc,gc,&expF,&expG);
	*FL=fc[L]*exp(expF);
	*GL=gc[L]*exp(expG);
	delete [] fc;
	delete [] gc;
}

// This calcules dCW/dx
void CoulWave::SphericalCWprime(int L,double x,double eta,double *FL,double *GL,double *FLprime,double *GLprime){
	double expF,expG;
	double *fc,*gc,*fcp,*gcp;
	int k=0;
	fc=new double[k+1];
	gc=new double[k+1];
	fcp=new double[k+1];
	gcp=new double[k+1];
	// This calculates fc and gc arrays for indices L to L+k
	gsl_sf_coulomb_wave_FGp_array(L,k,eta,x,fc,fcp,gc,gcp,&expF,&expG);
	*FL=fc[0]*exp(expF);
	*GL=gc[0]*exp(expG);
	*FLprime=fcp[0]*exp(expF);
	*GLprime=gcp[0]*exp(expG);
	delete [] fc;
	delete [] gc;
	delete [] fcp;
	delete [] gcp;
}

