#include "msu_commonutils/misc.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/sf.h"
#include "msu_commonutils/log.h"

bool Misc::comparestrings(char *s1,char *s2){
	int length1,length2,ic;
	bool answer;
	answer=0;
	length1=strlen(s1);
	length2=strlen(s2);
	if(length1==length2){
		answer=1;
		for(ic=0;ic<length1;ic++){
			if(s1[ic]!=s2[ic]){
				answer=0;
				return answer;
			}
		}
	}
	return answer;
}

double Misc::DotProduct(FourVector &p,FourVector &q){
	return p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3];
}

double Misc::triangle(double m0,double m1,double m2){
	double answer,m0sq,m1sq,m2sq;
	char message[CLog::CHARLENGTH];
	if(m0<m1+m2) {
		snprintf(message,CLog::CHARLENGTH,"Disaster with triangle, M=%g, m1=%g, m2=%g\n",m0,m1,m2);
		CLog::Fatal(message);
	}
	m0sq=m0*m0;m1sq=m1*m1;m2sq=m2*m2;
	answer=(m0sq*m0sq+m1sq*m1sq+m2sq*m2sq-2.0*(m0sq*m1sq+m0sq*m2sq+m1sq*m2sq))/(4.0*m0sq);
	return answer;
}

double Misc::triangle2(double m0sq,double m1sq,double m2sq){
	double answer;
	answer=(m0sq*m0sq+m1sq*m1sq+m2sq*m2sq-2.0*(m0sq*m1sq+m0sq*m2sq+m1sq*m2sq))/(4.0*m0sq);
	return answer;
	return (m0sq*m0sq+m1sq*m1sq+m2sq*m2sq
		-2.0*(m0sq*m1sq+m0sq*m2sq+m1sq*m2sq))/(4.0*m0sq);
}

double Misc::GetRapidity(FourVector &pa){
	return 0.5*log((pa[0]+pa[3])/(pa[0]-pa[3]));
}

double Misc::GetDely(FourVector &pa,FourVector &pb){
	return GetRapidity(pa)-GetRapidity(pb);
}

double Misc::GetQinv(FourVector &pa,FourVector &pb){
	double answer;
	answer=0;
	answer=pow(pa[1]-pb[1],2)+pow(pa[2]-pb[2],2)+pow(pa[3]-pb[3],2)-pow(pa[0]-pb[0],2);
	return sqrt(answer);
}

complex<double> Misc::cexp(complex<double> z){
	return exp(real(z))*ceiphi(imag(z));
}

complex<double> Misc::ceiphi(double phi){
	return complex<double>(cos(phi),sin(phi));
}

complex<double> Misc::cpow(complex<double> z,complex<double> a){
	double zr=real(z);
	double zi=imag(z);
	double phi=atan2(zi,zr);
	double zmag=sqrt(zr*zr+zi*zi);
	complex<double> alnz=a*(log(zmag)+ci*phi);
	return cexp(alnz);
}

int Misc::iround(double x){
	return int(floor(x+0.5));
}

// Using GSL
double Misc::cgc(double j1,double m1,double j2,double m2,double J,double M){
	int sign=-1;
	if((lrint((M+j1-j2)))%2==0) sign=1;
	return sign*sqrt(2.0*J+1)
		*gsl_sf_coupling_3j(lrint(2*j1),lrint(2*j2),lrint(2*J),
		lrint(2*m1),lrint(2*m2),-lrint(2*M));
}

// from Alexander Volya, Declan Mulhall
//page 44 edmonds
double Misc::oldcgc(double j1,double m1,double j2,double m2,double j,double m){
	//rule check
	if (fabs(m1+m2-m)>0.01) return 0.0; //sum of m's must be zero
	if (j1+j2<j) return 0.0; //triangle rule
	if (j+j1<j2) return 0.0; 
	double p1,p2,p3,p4,p5,thesum = 0,ans;
	double z=0;
	if(m1 + m2 == m && m<=j && -j <=m){
		p1=(2*j+1)*cgc_fractorial((j1+j2-j),(j1+j2+j+1));
		p2=cgc_fractorial((j1-m1),(j1-j2+j));
		p3=cgc_fractorial((j2-m2),(j2-j1+j));
		p4=cgc_fractorial((j+m),(j1+m1));
		p5=cgc_fractorial((j-m),(j2+m2));
		for(z=0;z<=2*(j1+j2)+j;z++){
			if(j1-m1-z >=0 && j-m-z >=0 && j2-j+m1+z >=0 ){
				thesum = thesum +
					pow(-1,z+j1-m1)*cgc_fractorial(j1+m1+z,j1-m1-z)
					*cgc_fractorial(j2+j-m1-z,j-m-z)
					/(cgc_factorial(z)*cgc_factorial(j2-j+m1+z));
			}
		}
		ans = sqrt(p1*p2*p3*p4*p5) * thesum;
		return ans;
	}
	else return 0;
}


double Misc::cgc_factorial(double n){
	double i,j=1;
	for (i=1;i<=n;i++) j=j*i;
	return j;
}

double Misc::cgc_fractorial(double n,double m){
	double  i,j=1;
	if(n>m){
		for(i=m+1;i<n+1;i++){j=j*i;}
		return j;
	}
	if(n<m){
		for(i=n+1;i<m+1;i++){j=j*i;}
		return 1/j;
	}
	return 1;
}

int Misc::cgc_delta (int x, int y) {if (x==y) return 1; else return 0;}

void Misc::outsidelong(FourVector &pa,FourVector &pb, double &qinv, double &qout, double &qside, double &qlong){
	double vs,gamma,ptot[4],q[4],qtemp,ptot_perp;
	int alpha;
	for(alpha=0;alpha<4;alpha++){
		q[alpha]=0.5*(pa[alpha]-pb[alpha]);
		ptot[alpha]=pa[alpha]+pb[alpha];
	}

	// Perform long. comoving boost
	vs=ptot[3]/ptot[0];
	gamma=1.0/sqrt(1.0-vs*vs);
	ptot[0]=gamma*(ptot[0]-vs*ptot[3]);
	ptot[3]=0.0;
	qtemp=q[3];
	q[3]=gamma*(q[3]-vs*q[0]);
	q[0]=gamma*(q[0]-vs*qtemp);

	ptot_perp=sqrt(ptot[1]*ptot[1]+ptot[2]*ptot[2]);
	qout=(q[1]*ptot[1]+q[2]*ptot[2])/ptot_perp;
	qside=fabs(q[2]*ptot[1]-q[1]*ptot[2])/ptot_perp;
	qlong=q[3];
	qinv=sqrt(qlong*qlong+qside*qside+qout*qout);

	vs=ptot_perp/ptot[0];
	gamma=1.0/sqrt(1.0-vs*vs);
	qout=gamma*(qout-vs*q[0]);

}

void Misc::Pause(){
	double pausedummy;
	printf("PAUSED, enter anything to continue: ");
	scanf("%lf",&pausedummy);
}

void Misc::Pause(int seconds){
	//const double CLOCKS_PER_SEC=1.0;
	clock_t endwait;
	endwait = clock () + seconds * CLOCKS_PER_SEC ;
	while (clock() < endwait) {}
}

double Misc::CalcDelta_FromSqWells(int ell,double mu,int nwells,double q,double *V0,double *r){
	int iwell;
	double norm,dd;
	double *qq,*A,*B;
	bool *realcheck;
	A=new double[nwells+1];
	B=new double[nwells+1];
	qq=new double[nwells+1];
	realcheck=new bool[nwells+1];
	double c,cprime,qr,jl,jlprime,nl,nlprime,ek;

	A[0]=1.0; B[0]=0.0;

	qq[nwells]=q;
	realcheck[nwells]=1;
	for(iwell=0;iwell<nwells;iwell++){
		realcheck[iwell]=1;
		ek=q*q-2.0*mu*V0[iwell];
		qq[iwell]=sqrt(fabs(ek));
		if(ek<0.0) realcheck[iwell]=0;
		//printf("q[%d]=%g, realcheck=%d\n",iwell,qq[iwell],realcheck[iwell]);
	}
	for(iwell=1;iwell<=nwells;iwell++){

		qr=qq[iwell-1]*r[iwell-1]/HBARC;
		if(realcheck[iwell-1]==1)	Bessel::CalcJN_real(ell,qr,jl,nl,jlprime,nlprime);
		else	Bessel::CalcJN_imag(ell,qr,jl,nl,jlprime,nlprime);
		c=A[iwell-1]*jl+B[iwell-1]*nl;
		cprime=qq[iwell-1]*(A[iwell-1]*jlprime+B[iwell-1]*nlprime);

		cprime=cprime/qq[iwell];
		qr=qq[iwell]*r[iwell-1]/HBARC;
		if(realcheck[iwell]==1) Bessel::CalcJN_real(ell,qr,jl,nl,jlprime,nlprime);
		else Bessel::CalcJN_imag(ell,qr,jl,nl,jlprime,nlprime);
		A[iwell]=(cprime*nl-c*nlprime)/(jlprime*nl-jl*nlprime);
		B[iwell]=(c*jlprime-cprime*jl)/(jlprime*nl-jl*nlprime);
	}

	norm=1.0/sqrt(A[nwells]*A[nwells]+B[nwells]*B[nwells]);
	for(iwell=0;iwell<=nwells;iwell++){
		A[iwell]*=norm;
		B[iwell]*=norm;
		//printf("iwell=%d, A=%g, B=%g\n",iwell,A[iwell],B[iwell]);
	}
	dd=atan2(B[nwells],A[nwells]);
	if(q<20) printf("A=%g, B=%g\n",A[nwells],B[nwells]);
	dd-=2.0*PI*rint(dd/(2.0*PI));
	if(dd<0.0) dd+=PI;

	delete [] A;
	delete [] B;
	delete [] qq;
	delete [] realcheck;
	return dd;

}

void Misc::Cubic(double a0,double a1,double a2,double a3,
complex<double>& z1, complex<double>& z2,complex<double>& z3){
	complex<double> x,w,check,dZdz,dz,z[3];
	double Q,R,p,q,arg;
	int i,ncheck=0;
	a0=a0/a3;
	a1=a1/a3;
	a2=a2/a3;
	p=a1-(a2*a2/3.0);
	q=(a1*a2/3.0)-a0-2.0*a2*a2*a2/27.0;
	Q=p/3.0;
	R=0.5*q;
	arg=R*R+Q*Q*Q;
		//printf("p=%g, q=%g, Q=%g,R=%g, arg=%g\n",p,q,Q,R,arg);
	if(arg>0.0){
		w=pow(fabs(R+sqrt(arg)),1.0/3.0);
		if(R+sqrt(arg)<0.0) w=-w;
		x=w-p/(3.0*w);
			//printf("X CHECK : 0 =? (%g,%g)\n",real(x*x*x+p*x-q),imag(x*x*x+p*x-q));
		z1=x-a2/3.0;
		w=w*exp(2.0*ci*PI/3.0);
		x=w-p/(3.0*w);
		z2=x-a2/3.0;
		w=w*exp(2.0*ci*PI/3.0);
		x=w-p/(3.0*w);
		z3=x-a2/3.0;
	}
	else{
		w=pow(R+ci*sqrt(-arg),1.0/3.0);
		x=w-p/(3.0*w);
		z1=x-a2/3.0;
		w=w*exp(2.0*ci*PI/3.0);
		x=w-p/(3.0*w);
		z2=x-a2/3.0;
		w=w*exp(2.0*ci*PI/3.0);
		x=w-p/(3.0*w);
		z3=x-a2/3.0;
	}
	z[0]=z1;
	z[1]=z2;
	z[2]=z3;

	for(i=0;i<3;i++){
		check=(pow(z[i],3)+a2*pow(z[i],2)+a1*z[i]+a0)
		/(abs(pow(z[i],3))+abs(a2*pow(z[i],2))+abs(a1*z[i]+a0));
		ncheck=0;
		while(abs(check)>1.0E-12){
			if(ncheck>10){
				printf("FAILURE IN CUBIC SOLVER: check=(%g,%g) =? 0\n",real(check),imag(check));
				exit(1);
			}
			dZdz=3.0*pow(z[i],2)+2.0*a2*z[i]+a1;
			dz=-(pow(z[i],3)+a2*pow(z[i],2)+a1*z[i]+a0)/dZdz;
			z[i]+=dz;
			check=(pow(z[i],3)+a2*pow(z[i],2)+a1*z[i]+a0)
			/(abs(pow(z[i],3))+abs(a2*pow(z[i],2))+abs(a1*z[i]+a0));
			ncheck+=1;
		}
	}
	z1=z[0]; z2=z[1]; z3=z[2];
}


double Misc::signswitch(double a, double b) {
	if (b>=0.0) return fabs(a);
	return -fabs(a);
}

void Misc::Quartic(double a0,double a1,double a2,double a3,double a4,complex<double> &z1, complex<double> &z2,complex<double> &z3,complex<double> &z4){
	complex<double> z[4];
	if(a4==0.0){
		printf("Leading coefficient zero in cern_quartic_real_coeff::solve_rc().\n");
		exit(-1);
	}
	Quartic(a0,a1,a2,a3,a4,z);
	z1=z[0]; z2=z[1]; z3=z[2]; z4=z[3];	
}

// There are a couple differences with the original routine.
// The arrays z[] and u[] are now zero-indexed.
//int cern_quartic_real_coeff::rrteq4(double a, double b, double c, double d, 
	//			    complex<double> z[], double &dc, 
		//		    int &mt) {
void Misc::Quartic(double a0,double a1,double a2,double a3,double a4,complex<double> *z){
	double dc;
	int mt;
	double a=a3/a4,b=a2/a4,c=a1/a4,d=a0/a4;

	complex<double> i(0.0,1.0), z0[5];
	complex<double> w1(0.0,0.0), w2(0.0,0.0), w3;
	double r4=1.0/4.0, r12=1.0/12.0;
	double q2=1.0/2.0, q4=1.0/4.0, q8=1.0/8.0;
	double q1=3.0/8.0, q3=3.0/16.0;
	double u[3], v[4], v1, v2;
	int j, k1=0, k2=0;

// degenerate cases
	if (b==0 && c==0) {
		if (d==0) {
			mt=1;
			z[0]=-a;
			z[1]=0;
			z[2]=0;
			z[3]=0;
			dc=0;
			return;
		}
		else if (a==0) {
			if (d>0) {
				mt=2;
				z[0]=sqrt(i*sqrt(d));
				z[1]=-z[0];
				z[3]=sqrt(-z[0]*z[0]);
				z[2]=-z[3];
			}
			else {
				mt=3;
				z[0]=sqrt(sqrt(-d));
				z[1]=-z[0];
				z[2]=sqrt(-z[0]*z[0]);
				z[3]=-z[2];
			}
			dc=-r12*d*r12*d*r12*d;
			return;
		}
	}

// Solve the resolvant cubic
	double aa=a*a;
	double pp=b-q1*aa;
	double qq=c-q2*a*(b-q4*aa);
	double rr=d-q4*(a*c-q4*aa*(b-q3*aa));
	double rc=q2*pp;
	double sc=q4*(q4*pp*pp-rr);
	double tc=-(q8*qq*q8*qq);

	//cub_obj.rrteq3(rc,sc,tc,u,dc);
	Misc::CubicResolvant(rc,sc,tc,u,dc);

	double q=qq;
	double h=r4*a;
	if (dc==0) u[2]=u[1];
	if (dc<=0) {
		mt=2;
		v[1]=fabs(u[0]);
		v[2]=fabs(u[1]);
		v[3]=fabs(u[2]);
		v1=max(max(v[1],v[2]),v[3]);
		if (v1==v[1]) {
			k1=0;
			v2=max(v[2],v[3]);
		} else if (v1==v[2]) {
			k1=1;
			v2=max(v[1],v[3]);
		} else {
			k1=2;
			v2=max(v[1],v[2]);
		}
		if (v2==v[1]) {
			k2=0;
		} else if (v2==v[2]) {
			k2=1;
		} else {
			k2=2;
		}
		w1=sqrt(static_cast<complex<double>>(u[k1]));
		w2=sqrt(static_cast<complex<double>>(u[k2]));
	} else {
		mt=3;
		w1=sqrt(u[1]+i*u[2]);
		w2=sqrt(u[1]-i*u[2]);
	}
	w3=0;
	if (w1*w2!=0.0) w3=-q/(8.0*w1*w2);
	z0[1]=w1+w2+w3-h;
	z0[2]=-w1-w2+w3-h;
	z0[3]=-w1+w2-w3-h;
	z0[4]=w1-w2-w3-h;
	if (mt==2) {
		if (u[k1]>=0 && u[k2]>=0) {
			mt=1;
			for(j=1;j<=4;j++) {
				z[j-1]=z0[j].real();
			}
		} else if (u[k1]>=0 && u[k2]<0) {
			z[0]=z0[1];
			z[1]=z0[4];
			z[2]=z0[3];
			z[3]=z0[2];
		} else if (u[k1]<0 && u[k2]>=0) {
			z[0]=z0[1];
			z[1]=z0[3];
			z[2]=z0[4];
			z[3]=z0[2];
		} else if (u[k1]<0 && u[k2]<0) {
			z[0]=z0[1];
			z[1]=z0[2];
			z[2]=z0[4];
			z[3]=z0[3];
		}
	} else if (mt==3) {
		for(j=1;j<=2;j++) {
			z[j-1]=z0[j].real();
		}
		z[2]=z0[4];
		z[3]=z0[3];
	}
}


// This essentially follows the original Fortran code
// except the roots are returned in x[0], x[1], and x[2]
// instead of x[1], x[2], and x[3]. (The arrays y and z
// were already zero-indexed in the CERN routine.)
// Similar to the original, in the case of complex roots, x[0]
// is the real root and x[1] and x[2] contain the real and
// imaginary parts of the complex roots.
//int cern_cubic_real_coeff::rrteq3(double r, double s, double t, 
//double x[], double &d) {
void Misc::CubicResolvant(double r,double s,double t,double x[],double &d){
	double eps=1.0e-6, delta=1.0e-15;
	double r1=2.0/27.0, r2=0.5, r3=1.0/3.0;
	double w3=sqrt(3.0), r4=w3/2.0;
	double q1=2.0/27.0, q2=0.5, q3=1.0/3.0;
	double y[3];
	complex<double> z[3], i(0.0,1.0);
	double h2, h3;
	int j,k;

	if (s==0.0 && t==0.0) {
		x[0]=-r;
		x[1]=0.0;
		x[2]=0.0;
		d=0;
		return;
	}
	double p=s-r3*r*r;
	double q=(r1*r*r-r3*s)*r+t;
	d=r2*r2*q*q+r3*p*r3*p*r3*p;
	if (fabs(d)<=eps) {
		double pp=s-q3*r*r;
		double qq=(q1*r*r-q3*s)*r+t;
		d=q2*q2*qq*qq+q3*pp*q3*pp*q3*pp;
		p=pp;
		q=qq;
	}
	double h=r3*r;
	double h1=r2*q;
	double u,v,d_new;
// AWS hack 
// The discriminant in 'd' has units of [x]^6 so it is very
// sensitive to the absolute magnitude of the roots.  We attempt to
// fix this by using the ratio instead of the sum.
	if (true) {
		double da=r2*r2*q*q;
		double db=r3*p*r3*p*r3*p;
		if (db==0.0) {
			delta=0.0;
			d_new=da;
		} else if (db>0.0) {
			d_new=da/db+1.0;
		} else {
			d_new=-da/db-1.0;
		}
	} else {
		d_new=d;
	}
	if (d_new>delta) {
		h2=sqrt(d);
		double u0=-h1+h2;
		double v0=-h1-h2;
		if (fabs(u0)==0.0) u=Misc::signswitch(0.0,u0);
		else u=Misc::signswitch(pow(fabs(u0),r3),u0);
		if (fabs(v0)==0.0) v=Misc::signswitch(0.0,v0);
		else v=Misc::signswitch(pow(fabs(v0),r3),v0);
		x[0]=u+v-h;
		x[1]=-r2*(u+v)-h;
		x[2]=r4*fabs(u-v);
		if (fabs(u0)<=eps || fabs(v0)<=eps) {
			y[0]=x[0];
			for(k=0;k<=1;k++) {
				y[k+1]=y[k]-(((y[k]+r)*y[k]+s)*y[k]+t)/((3.0*y[k]+2.0*r)*y[k]+s);
			}
			x[0]=y[2];
			z[0]=x[1]+i*x[2];
			for(k=0;k<=1;k++) {
				z[k+1]=z[k]-(((z[k]+r)*z[k]+s)*z[k]+t)/((3.0*z[k]+2.0*r)*z[k]+s);
			}
			x[1]=z[2].real();
			x[2]=z[2].imag();
		}
	} else if (fabs(d_new)<=delta) {
		d=0.0;
		if (fabs(h1)==0.0) u=Misc::signswitch(0.0,-h1);
		else u=Misc::signswitch(pow(fabs(h1),r3),-h1);
		x[0]=u+u-h;
		x[1]=-u-h;
		x[2]=x[1];
		if (fabs(h1)<=eps) {
			y[0]=x[0];
			for(k=0;k<=1;k++) {
				h1=(3.0*y[k]+2.0*r)*y[k]+s;
				if (fabs(h1)>delta) {
					y[k+1]=y[k]-(((y[k]+r)*y[k]+s)*y[k]+t)/h1;
				} else {
					x[0]=-r3*r;
					x[1]=x[0];
					x[2]=x[0];
					return;
				}
			}
			x[0]=y[2];
			x[1]=-r2*(r+x[0]);
			x[2]=x[1];
		}
	}
	else {
		h3=fabs(r3*p);
		h3=sqrt(h3*h3*h3);
		h2=r3*acos(-h1/h3);
		if (h3==0.0) h1=0.0;
		else h1=pow(h3,r3);
		u=h1*cos(h2);
		v=w3*h1*sin(h2);
		x[0]=u+u-h;
		x[1]=-u-v-h;
		x[2]=-u+v-h;
		if (h3<=eps || x[0]<=eps || x[1]<=eps || x[2]<=eps) {
			for(j=0;j<3;j++) {
				y[0]=x[j];
				for(k=0;k<=1;k++) {
					y[k+1]=y[k]-(((y[k]+r)*y[k]+s)*y[k]+t)/((3.0*y[k]+2.0*r)*y[k]+s);
				}
				x[j]=y[2];
			}
		}
	}
}


/* int cern_cubic_real_coeff::solve_rc
(const double a3, const double b3, const double c3, const double d3, 
double &r1, complex<double> &r2, complex<double> &r3) {

	if (a3==0.0) {
		cout << 
			"Leading coefficient zero in cern_cubic_real_coeff::solve_rc()."
			<< endl;
		exit(-1);
	}

	double x[3],d;
	complex<double> i(0.0,1.0);

	rrteq3(b3/a3,c3/a3,d3/a3,x,d);
	if (d>0.0) {
		r1=x[0];
		r2=x[1]+i*x[2];
		r3=x[1]-i*x[2];
	} else {
		r1=x[0];
		r2=x[1];
		r3=x[2];
	}

	return 0;
}

int main(void) {

	cern_quartic_real_coeff cq;
	std::complex<double> x1,x2,x3,x4;

	cq.solve_rc(1.0,-10.0,35.0,-50.0,24.0,x1,x2,x3,x4);
	cout << x1 << " " << x2 << " " << x3 << " " << x4 << endl;
	return 0;
} */

int Misc::CubicReal(double a0,double a1,double a2,double a3,double *x){
	a0=a0/a3; a1=a1/a3; a2=a2/a3;
	int nrealroots=gsl_poly_solve_cubic(a2,a1,a0,&x[0],&x[1],&x[2]);
	return nrealroots;
}

void Misc::CubicComplex(double a0,double a1,double a2,double a3,complex<double> &z1,complex<double> &z2,complex<double> &z3){
	a0=a0/a3; a1=a1/a3; a2=a2/a3;
	gsl_complex gslz1,gslz2,gslz3;
	gsl_poly_complex_solve_cubic(a2,a1,a0,&gslz1,&gslz2,&gslz3);
	z1=GSL_REAL(gslz1)+ci*GSL_IMAG(gslz1);
	z2=GSL_REAL(gslz2)+ci*GSL_IMAG(gslz2);
	z3=GSL_REAL(gslz3)+ci*GSL_IMAG(gslz3);
}
			
bool Misc::file_exists(const char* file_name){
  std::ifstream file; 
  file.open(file_name,std::ios::in);
  bool result;
  if (file.fail()) result = false;
  else result = true;
  file.close();
  return result;
}

bool Misc::file_exists( const std::string& file_name ){
    return file_exists( file_name.c_str() );
}

double Misc::mod(double x, double y){
    return x-int(x/y)*y;
}

int Misc::Sign(int a){
	if (a < 0){
		return 0;
	}else{
		return 1;
	}
}

void Misc::outsidelong(FourVector &pa,FourVector &pb, double &qinv, double &qout,
double &qside,double &qlong,double &deleta,double &dely,double &delphi){ // qout is in pair frame
	// q.. refer to half relative momenta in pair CM frame
	double vs,gamma,ptot[4],q[4],qtemp,ptot_perp,pmaga,pmagb,ya,yb,etaa,etab,phia,phib;
	int alpha;
	for(alpha=0;alpha<4;alpha++){
		q[alpha]=0.5*(pa[alpha]-pb[alpha]);
		ptot[alpha]=pa[alpha]+pb[alpha];
	}
	// Perform long. comoving boost
	vs=ptot[3]/ptot[0];
	gamma=1.0/sqrt(1.0-vs*vs);
	ptot[0]=gamma*(ptot[0]-vs*ptot[3]);
	ptot[3]=0.0;
	qtemp=q[3];
	q[3]=gamma*(q[3]-vs*q[0]);
	q[0]=gamma*(q[0]-vs*qtemp);

	ptot_perp=sqrt(ptot[1]*ptot[1]+ptot[2]*ptot[2]);
	qout=(q[1]*ptot[1]+q[2]*ptot[2])/ptot_perp;
	qside=(q[2]*ptot[1]-q[1]*ptot[2])/ptot_perp;
	qlong=(q[3]);

	vs=ptot_perp/ptot[0];
	gamma=1.0/sqrt(1.0-vs*vs);
	qout=gamma*(qout-vs*q[0]);
	qinv=sqrt(qout*qout+qlong*qlong+qside*qside);
	
	pmaga=sqrt(pa[1]*pa[1]+pa[2]*pa[2]+pa[3]*pa[3]);
	pmagb=sqrt(pb[1]*pb[1]+pb[2]*pb[2]+pb[3]*pb[3]);
	etaa=0.5*log((pmaga+pa[3])/(pmaga-pa[3]));
	etab=0.5*log((pmagb+pb[3])/(pmagb-pb[3]));
	ya=0.5*log((pa[0]+pa[3])/(pa[0]-pa[3]));
	yb=0.5*log((pb[0]+pb[3])/(pb[0]-pb[3]));
	phia=atan2(pa[2],pa[1]);
	phib=atan2(pb[2],pb[1]);
	deleta=fabs(etaa-etab);
	dely=fabs(ya-yb);
	delphi=fabs(phib-phia)*180.0/PI;
	if(delphi>180.0)
		delphi=360.0-delphi;
}

void Misc::outsidelong_lcms(FourVector &pa,FourVector &pb, double &qinv, double &qout,double &qout_lcms,
double &qside,double &qlong,double &deleta,double &dely,double &delphi){
	// q.. refer to half relative momenta in pair CM frame, qout is in pair frame, qout_lcms is in LCMS frame
	double vs,gamma,ptot[4],q[4],qtemp,ptot_perp,pmaga,pmagb,ya,yb,etaa,etab,phia,phib;
	int alpha;
	for(alpha=0;alpha<4;alpha++){
		q[alpha]=0.5*(pa[alpha]-pb[alpha]);
		ptot[alpha]=pa[alpha]+pb[alpha];
	}
	// Perform long. comoving boost
	vs=ptot[3]/ptot[0];
	gamma=1.0/sqrt(1.0-vs*vs);
	ptot[0]=gamma*(ptot[0]-vs*ptot[3]);
	ptot[3]=0.0;
	qtemp=q[3];
	q[3]=gamma*(q[3]-vs*q[0]);
	q[0]=gamma*(q[0]-vs*qtemp);

	ptot_perp=sqrt(ptot[1]*ptot[1]+ptot[2]*ptot[2]);
	qout_lcms=qout=(q[1]*ptot[1]+q[2]*ptot[2])/ptot_perp;
	qside=(q[2]*ptot[1]-q[1]*ptot[2])/ptot_perp;
	qlong=q[3];

	vs=ptot_perp/ptot[0];
	gamma=1.0/sqrt(1.0-vs*vs);
	qout=gamma*(qout_lcms-vs*q[0]);
	qinv=sqrt(qout*qout+qlong*qlong+qside*qside);
	
	pmaga=sqrt(pa[1]*pa[1]+pa[2]*pa[2]+pa[3]*pa[3]);
	pmagb=sqrt(pb[1]*pb[1]+pb[2]*pb[2]+pb[3]*pb[3]);
	etaa=0.5*log((pmaga+pa[3])/(pmaga-pa[3]));
	etab=0.5*log((pmagb+pb[3])/(pmagb-pb[3]));
	ya=0.5*log((pa[0]+pa[3])/(pa[0]-pa[3]));
	yb=0.5*log((pb[0]+pb[3])/(pb[0]-pb[3]));
	phia=atan2(pa[2],pa[1]);
	phib=atan2(pb[2],pb[1]);
	deleta=fabs(etaa-etab);
	dely=fabs(ya-yb);
	delphi=fabs(phib-phia)*180.0/PI;
	if(delphi>180.0)
		delphi=360.0-delphi;
}

void Misc::CalcDCA(FourVector &p,FourVector &r,FourVector &dca){
	int alpha;
	double pdotr,p2;
	p2=p[1]*p[1]+p[2]*p[2]+p[3]*p[3];
	pdotr=(p[1]*r[1]+p[2]*r[2]+p[3]*r[3])/p2;
	for(alpha=1;alpha<4;alpha++)
		dca[alpha]=(r[alpha]-pdotr*p[alpha])/1.0E13;
	dca[0]=sqrt(dca[1]*dca[1]+dca[2]*dca[2]+dca[3]*dca[3]);
}
