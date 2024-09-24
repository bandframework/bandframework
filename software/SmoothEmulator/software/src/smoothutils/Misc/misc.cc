#include "msu_smoothutils/misc.h"
#include "msu_smoothutils/log.h"
using namespace NMSUUtils;

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

void Misc::Cubic(double a0,double a1,double a2,double a3,
complex<double>& z1, complex<double>& z2,complex<double>& z3){
	complex<double> ci(0.0,1.0);
	complex<double> x,w,check,dZdz,dz,z[3];
	const double PI=4.0*atan(1.0);
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

double Misc::signswitch(double a, double b) {
	if (b>=0.0) return fabs(a);
	return -fabs(a);
}

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
