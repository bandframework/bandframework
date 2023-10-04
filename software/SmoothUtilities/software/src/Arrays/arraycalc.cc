#include "msu_commonutils/arrays.h"
#include "msu_commonutils/sf.h"
using namespace std;

void ArrayCalc::CalcAExpArrayFromMArray(CCHArray *M,int irm,CCHArray *A,int ira){
	int dlx=1,dly=1,dlz=1;
	if(A->GetXSYM()) dlx=2;
	if(A->GetYSYM()) dly=2;
	if(A->GetZSYM()) dlz=2;

	int LMAX=A->GetLMAX();
	if(M->GetLMAX()!=LMAX){
		printf("LMAX from M does not match LMAX from A\n");
		exit(1);
	}
	int L,lx,ly,lz;
	A->SetElement(0,0,0,ira,M->GetElement(0,0,0,ira));

	for(lx=0;lx<=1;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				L=lx+ly+lz;
				if(L>0) A->SetElement(lx,ly,lz,ira,
					M->GetAExpElementFromMArray(lx,ly,lz,irm));
			}
		}
	}
	A->FillRemainderX(ira);
}

void ArrayCalc::CalcMArrayFromAExpArray(CCHArray *A,int ira,CCHArray *M,int irm){
	int LMAX=A->GetLMAX();
	if(M->GetLMAX()!=LMAX){
		printf("LMAX from M does not match LMAX from A\n");
		exit(1);
	}
	int lx,ly,lz,L;
	int dlx=1,dly=1,dlz=1;
	if(M->GetXSYM()) dlx=2;
	if(M->GetYSYM()) dly=2;
	if(M->GetZSYM()) dlz=2;

	M->SetElement(0,0,0,irm,A->GetElement(0,0,0,ira));
	for(lx=0;lx<=LMAX;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				L=lx+ly+lz;
				if(L>0) M->SetElement(lx,ly,lz,irm,
					A->GetMElementFromAExpArray(lx,ly,lz,ira));
			}
		}	     
	}
}

void ArrayCalc::AddArrays(CCHArray *A,int ira,CCHArray *B,int irb,CCHArray *C,int irc){
	int LMAX,lx,ly,lz;
	int dlx=1,dly=1,dlz=1;
	if(C->GetXSYM()) dlx=2;
	if(C->GetYSYM()) dly=2;
	if(C->GetZSYM()) dlz=2;
	LMAX=A->GetLMAX();
	if(LMAX>B->GetLMAX()) LMAX=B->GetLMAX();
	for(lx=0;lx<=LMAX;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				C->SetElement(lx,ly,lz,irc,A->GetElement(lx,ly,lz,ira)
					+B->GetElement(lx,ly,lz,irb));
			}
		}
	}
}

void ArrayCalc::SubtractArrays(CCHArray *A,int ira,CCHArray *B,int irb,CCHArray *C,int irc){
	int LMAX,lx,ly,lz;
	int dlx=1,dly=1,dlz=1;
	if(C->GetXSYM()) dlx=2;
	if(C->GetYSYM()) dly=2;
	if(C->GetZSYM()) dlz=2;
	LMAX=A->GetLMAX();
	if(LMAX>B->GetLMAX()) LMAX=B->GetLMAX();
	for(lx=0;lx<=LMAX;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				C->SetElement(lx,ly,lz,irc,A->GetElement(lx,ly,lz,ira)
					-B->GetElement(lx,ly,lz,irb));
			}
		}
	}
}

void ArrayCalc::AddArrays(CCHArray *A,CCHArray *B,CCHArray *C){
	int ir;
	if((!CompareArrayParameters(A,B)) || (!CompareArrayParameters(B,C))){
		printf("fatal error, parameters mismatch in ArrayCalc::AddArrays\n");
		exit(1);
	}
	for(ir=0;ir<A->GetNRADIAL();ir++) AddArrays(A,ir,B,ir,C,ir);
}

void ArrayCalc::SubtractArrays(CCHArray *A,CCHArray *B,CCHArray *C){
	int ir;
	if((!CompareArrayParameters(A,B)) || (!CompareArrayParameters(B,C))){
		printf("fatal error, parameters mismatch in ArrayCalc::AddArrays\n");
		exit(1);
	}
	for(ir=0;ir<A->GetNRADIAL();ir++) SubtractArrays(A,ir,B,ir,C,ir);
}

void ArrayCalc::DivideArrays(CCHArray *A,CCHArray *B,CCHArray *C){
	int ir;
	if((A->GetNRADIAL()!=B->GetNRADIAL()) || (B->GetNRADIAL()!=C->GetNRADIAL())){
		printf("FATAL: NRADIAL mismatch in ArrayCalc::DivideArrays\n");
		exit(1);
	}
	for(ir=0;ir<A->GetNRADIAL();ir++) DivideArrays(A,ir,B,ir,C,ir);
}

void ArrayCalc::MultiplyArrays(CCHArray *A,CCHArray *B,CCHArray *C){
	int ir;
	if((!CompareArrayParameters(A,B)) || (!CompareArrayParameters(B,C))){
		printf("fatal error, parameters mismatch in ArrayCalc::MultiplyArrays\n");
		exit(1);
	}
	for(ir=0;ir<A->GetNRADIAL();ir++) MultiplyArrays(A,ir,B,ir,C,ir);
}

void ArrayCalc::CopyArray(CCHArray *A,CCHArray *B){
	int ir;
	if(!CompareArrayParameters(A,B)){
		printf("fatal error, parameters mismatch in ArrayCalc::CopyArray\n");
		exit(1);
	}
	for(ir=0;ir<A->GetNRADIAL();ir++) CopyArray(A,ir,B,ir);
}

void ArrayCalc::CalcMArrayFromAExpArray(CCHArray *A,CCHArray *M){
	int ir;
	if(!CompareArrayParameters(A,M)){
		printf("fatal error, parameters mismatch in ArrayCalc::CalcMArrayFromAExpArray\n");
		exit(1);
	}
	for(ir=0;ir<A->GetNRADIAL();ir++) CalcMArrayFromAExpArray(A,ir,M,ir);
}

void ArrayCalc::CalcAExpArrayFromMArray(CCHArray *M,CCHArray *A){
	int ir;
	if(!CompareArrayParameters(A,M)){
		printf("fatal error, parameters mismatch in ArrayCalc::CalcAExpArrayFromMArray\n");
		exit(1);
	}
	for(ir=0;ir<A->GetNRADIAL();ir++) CalcAExpArrayFromMArray(M,ir,A,ir);
}

void ArrayCalc::CalcAExpArrayFromXExpArray(CCHArray *X,CCHArray *A){
	int ir;
	if(!CompareArrayParameters(A,X)){
		printf("fatal error, parameters mismatch in  ArrayCalc::CalcAExpArrayFromXExpArray\n");
		exit(1);
	}
	for(ir=0;ir<A->GetNRADIAL();ir++) CalcAExpArrayFromXExpArray(X,ir,A,ir);
}

void ArrayCalc::CalcXExpArrayFromAExpArray(CCHArray *A,CCHArray *X){
	int ir;
	if(!CompareArrayParameters(A,X)){
		printf("fatal error, parameters mismatch in ArrayCalc::CalcXExpArrayFromAExpArray\n");
		exit(1);
	}
	for(ir=0;ir<A->GetNRADIAL();ir++) CalcXExpArrayFromAExpArray(X,ir,A,ir);
}

void ArrayCalc::Detrace(CCHArray *M,CCHArray *A){
	int ir;
	if(!CompareArrayParameters(M,A)){
		printf("fatal error, parameters mismatch in ArrayCalc::Detrace\n");
		exit(1);
	}
	for(ir=0;ir<A->GetNRADIAL();ir++) Detrace(M,ir,A,ir);
}

void ArrayCalc::CopyArray(CCHArray *A,int ira,CCHArray *B,int irb){
	int lx,ly,lz,LMAX,LMAXA,L;
	int dlx=1,dly=1,dlz=1;
	if(A->GetXSYM()) dlx=2;
	if(A->GetYSYM()) dly=2;
	if(A->GetZSYM()) dlz=2;
	LMAX=B->GetLMAX();
	LMAXA=A->GetLMAX();
	for(lx=0;lx<=LMAX;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				L=lx+ly+lz;
				if(L<=LMAXA)
					B->SetElement(lx,ly,lz,irb,A->GetElement(lx,ly,lz,ira));
				else
					B->SetElement(lx,ly,lz,irb,0);
			}
		}
	}
}

void ArrayCalc::CalcYlmExpArrayFromAExpArray(CCHArray *A,int ira,CYlmArray *YlmArray,int irlm){
	CCHCalc chcalc;
	const double PI=4.0*atan(1.0);
	int LMAX,L,lx,ly,lz,M;
	double flm;
	complex<double> factor,ci(0.0,1.0);

	int dlx=1,dly=1,dlm=1,delL=1,Mstart;
	if(A->GetXSYM()) dlx=2;
	if(A->GetZSYM()) dlm=2;
	if(A->GetXSYM() && A->GetYSYM() && A->GetZSYM()) delL=2;
	if(A->GetYSYM()) dly=2;

	LMAX=YlmArray->GetLMAX();

	YlmArray->SetElement(0,0,irlm,A->GetElement(0,0,0,ira));
	for(L=delL;L<=LMAX;L+=delL){
		Mstart=0;
		if(A->GetZSYM()) Mstart=L%2;
		for(M=Mstart;M<=L;M+=dlm){
			YlmArray->SetElement(L,M,irlm,0.0);
			lz=L-M;
			if(dlm==1 || lz%2==0){
				flm=pow(-1.0,M)/sqrt((2.0*L+1)*4*PI*chcalc.Factorial(L+M)*chcalc.Factorial(L-M));
				flm=flm*chcalc.Factorial(M)*sqrt(4.0*PI)*chcalc.Factorial(L);

				for(lx=0;lx<=M;lx+=dlx){
					ly=M-lx;
					if(dly==1 || ly%2==0){
						factor=flm*pow(ci,M-lx)/(chcalc.Factorial(lx)*chcalc.Factorial(M-lx));
						YlmArray->IncrementElement(L,M,irlm,
							factor*A->GetElement(lx,ly,lz,ira));

					}
				}
			}
		}
	}  
}

void ArrayCalc::CalcAExpArrayFromYlmExpArray(CYlmArray *YlmArray,int irlm,CCHArray *A,int ira){
	CCHCalc chcalc;
	const double PI=4.0*atan(1.0);
	int LMAX,L,lx,ly,lz,M,k;
	complex<double> factor,ci(0.0,1.0);
	int dlx=1,dly=1,dlz=1;
	if(A->GetXSYM()) dlx=2;
	if(A->GetYSYM()) dly=2;
	if(A->GetZSYM()) dlz=2;

	LMAX=A->GetLMAX();

	A->SetElement(0,0,0,ira,real(YlmArray->GetElement(0,0,irlm)));
	for(lx=0;lx<=1;lx+=dlx){
		for(ly=0;ly<=LMAX-lx;ly+=dly){
			for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
				L=lx+ly+lz;
				if(L>0){
					A->SetElement(lx,ly,lz,ira,0);
					for(M=0;M<=L;M++){
						factor=0.0;
						for(k=0;k<=M;k++){
							factor+=pow(-ci,M-k)
								*(chcalc.Factorial(M)/(chcalc.Factorial(k)*chcalc.Factorial(M-k)))*chcalc.GetOverlap(lx,ly,lz,k,M-k,L-M);
						}
						factor*=pow(-1.0,M)*double(chcalc.DoubleFactorial(2*L-1));
						factor*=sqrt(double(2*L+1)/(4*PI*chcalc.Factorial(L+M)*chcalc.Factorial(L-M)));
						factor*=chcalc.DoubleFactorial(2*L+1)/(chcalc.Factorial(L)*sqrt(4*PI));

						if(M==0)
							A->IncrementElement(lx,ly,lz,ira,
							real(factor
							*YlmArray->GetElement(L,M,irlm)));
						else
							A->IncrementElement(lx,ly,lz,ira,
							real(2.0*factor
							*YlmArray->GetElement(L,M,irlm)));

					}
				}  
			}
		}
	}
	A->FillRemainderX(ira);
}

// This takes a function exapanded as:
// \sum_{\vec\ell} M_{\vec\ell} e_x^{\ell_x}e_y^{\ell_y}e_z^{\ell_z}
// and returns the array of exp coeff.s of Cartesian Harmonics
// If M satisfies traceless condition, A=M
void ArrayCalc::Detrace(CCHArray *M,int irm,CCHArray *A,int ira){

	CCHCalc chcalc;
	int dlx=1,dly=1,dlz=1;
	double factor;
	const double PI=4.0*atan(1.0);

	//CompareArrayParameters(M,A);
	if(M->GetXSYM()) dlx=2;
	if(M->GetYSYM()) dly=2;
	if(M->GetZSYM()) dlz=2;

	A->ZeroArray(ira);
	int L,LMAX,lx,ly,lz,m,mx,my,mz,lxprime,lyprime,lzprime;
	int LMAXM=M->GetLMAX();
	int LMAXA=A->GetLMAX();
	LMAX=LMAXM;
	if(LMAXA>LMAX) LMAX=LMAXA;

	CCHArray *Atilde;
	Atilde=new CCHArray(LMAXA,1,1.0,M->GetXSYM(),
		M->GetYSYM(),M->GetZSYM());
	for(lx=0;lx<=LMAXA;lx+=dlx){
		for(ly=0;ly<=LMAXA-lx;ly+=dly){
			for(lz=0;lz<=LMAXA-lx-ly;lz+=dlz){
				L=lx+ly+lz;
				for(mx=0;mx<=(LMAXM-L)/2;mx++){
					for(my=0;my<=(LMAXM-L-2*mx)/2;my++){
						for(mz=0;mz<=(LMAXM-L-2*mx-2*my)/2;mz++){
							m=mx+my+mz;
							if(L+2*m>LMAXM){
								printf("L+2m is tooooo big\n");
								exit(1);
							}
							factor=pow(0.5,m);
							factor*=chcalc.Factorial(L+2*m)*chcalc.DoubleFactorial(2*L+1)
								/(chcalc.Factorial(L)*chcalc.Factorial(mx)*chcalc.Factorial(my)*chcalc.Factorial(mz)*chcalc.DoubleFactorial(2*L+2*m+1));
							Atilde->IncrementElement(lx,ly,lz,0,
								factor*M->GetElement(lx+2*mx,ly+2*my,
								lz+2*mz,irm));
						}
					}
				}
			}
		}
	}
	for(lx=0;lx<=1;lx+=dlx){
		for(ly=0;ly<=LMAXA-lx;ly+=dly){
			for(lz=0;lz<=LMAXA-lx-ly;lz+=dlz){
				L=lx+ly+lz;
				for(lxprime=lx;lxprime<=L;lxprime+=2){
					for(lyprime=ly%2;lyprime<=L-lxprime;lyprime+=2){
						lzprime=L-lxprime-lyprime;
						factor=chcalc.DoubleFactorial(2*L+1)/(chcalc.Factorial(L)*4.0*PI);
						factor*=chcalc.GetOverlap(lx,ly,lz,lxprime,lyprime,lzprime)
							*Atilde->GetElement(lxprime,lyprime,lzprime,0);
						factor*=chcalc.Factorial(L)/(chcalc.Factorial(lxprime)*chcalc.Factorial(lyprime)*chcalc.Factorial(lzprime));
						A->IncrementElement(lx,ly,lz,ira,factor);
					}
				}
			}
		}
	}
	A->FillRemainderX(ira);
	delete Atilde;
}

void ArrayCalc::MultiplyArrays(CCHArray *A,int ira,CCHArray *B,int irb,CCHArray *C,int irc){
	int LMAXA=A->GetLMAX();
	int LMAXB=B->GetLMAX();
	int LMAXC=C->GetLMAX();
	bool XSYMd=0,YSYMd=0,ZSYMd=0;

	if(A->GetXSYM() && B->GetXSYM()) XSYMd=1; 
	if(A->GetYSYM() && B->GetYSYM()) YSYMd=1; 
	if(A->GetZSYM() && B->GetZSYM()) ZSYMd=1; 

	CCHArray *D;
	D=new CCHArray(LMAXC,1,1.0,XSYMd,YSYMd,ZSYMd);
	MultiplyArrays_Partial(LMAXA,A,ira,LMAXB,B,irb,LMAXC,D,0);
	Detrace(D,0,C,irc);
	delete D;

}

void ArrayCalc::MultiplyArrays_Partial(int LMAXA,CCHArray *A,int ira,int LMAXB,CCHArray *B,int irb,int LMAXC,CCHArray *C,int irc){
	CCHCalc chcalc;
	int dlxa=1,dlya=1,dlza=1;
	if(A->GetXSYM()) dlxa=2;
	if(A->GetYSYM()) dlya=2;
	if(A->GetZSYM()) dlza=2;
	int dlxb=1,dlyb=1,dlzb=1;
	if(B->GetXSYM()) dlxb=2;
	if(B->GetYSYM()) dlyb=2;
	if(B->GetZSYM()) dlzb=2;

	if(LMAXC>LMAXA+LMAXB) LMAXC=LMAXA+LMAXB;
	if(LMAXC>C->GetLMAX()) LMAXC=C->GetLMAX();
	if(LMAXA>LMAXC) LMAXA=LMAXC;
	if(LMAXB>LMAXC) LMAXB=LMAXC;

	int La,Lb,lxa,lya,lza,lxb,lyb,lzb,lxc,lyc,lzc,lmaxbprime;
	double gammaa,gammab,gammac,delC;

	C->ZeroArray(irc);

	for(lxa=0;lxa<=LMAXA;lxa+=dlxa){
		for(lya=0;lya<=LMAXA-lxa;lya+=dlya){
			for(lza=0;lza<=LMAXA-lxa-lya;lza+=dlza){
				La=lxa+lya+lza;
				gammaa=chcalc.Trinomial(lxa,lya,lza);
				lmaxbprime=LMAXC-La;
				if(lmaxbprime>LMAXB) lmaxbprime=LMAXB;
				for(lxb=0;lxb<=lmaxbprime;lxb+=dlxb){
					for(lyb=0;lyb<=lmaxbprime-lxb;lyb+=dlyb){
						for(lzb=0;lzb<=lmaxbprime-lxb-lyb;lzb+=dlzb){
							Lb=lxb+lyb+lzb;
							if(La+Lb<=LMAXC){
								lxc=lxa+lxb;
								lyc=lya+lyb;
								lzc=lza+lzb;
								gammab=chcalc.Trinomial(lxb,lyb,lzb);
								gammac=chcalc.Trinomial(lxc,lyc,lzc);
								delC=gammaa*gammab*A->GetElement(lxa,lya,lza,ira)
									*B->GetElement(lxb,lyb,lzb,irb)/gammac;
								C->IncrementElement(lxc,lyc,lzc,irc,delC);
							}
						}
					}
				}
			}
		}
	}
}

void ArrayCalc::DivideArrays(CCHArray *A,int ira,CCHArray *B,int irb,CCHArray *C,int irc){
	if( ((!A->GetXSYM() || !B->GetXSYM()) && C->GetXSYM())
		|| ((!A->GetYSYM() || !B->GetYSYM()) && C->GetYSYM())
	|| ((!A->GetZSYM() || !B->GetZSYM()) && C->GetZSYM())){
		printf("FATAL: Symmetry mismatch in ArrayCalc::DivdeArrays\n");
		exit(1);
	}
	CCHCalc chcalc;
	CCHArray *D;
	int LMAXA=A->GetLMAX();
	int LMAXB=B->GetLMAX();
	int LMAXD=55;
	bool XSYMd=0,YSYMd=0,ZSYMd=0;
	if(A->GetXSYM() && B->GetXSYM()) XSYMd=1; 
	if(A->GetYSYM() && B->GetYSYM()) YSYMd=1; 
	if(A->GetZSYM() && B->GetZSYM()) ZSYMd=1; 
	D=new CCHArray(LMAXD,1,1.0,XSYMd,YSYMd,ZSYMd);
	int Lb,Ld,lxb,lyb,lzb,lxd,lyd,lzd,m,mx,my,mz;
	double gammad,delD,delDsumtest=1.0,olddelDsumtest=1.0,delDsum;
	double qthresh,Asum=0.0,B0=B->GetElement(0,0,0,irb);
	bool qtest=0;
	const double EPSILON=1.0E-6;
	int dlxa=1,dlya=1,dlza=1;
	if(A->GetXSYM()) dlxa=2;
	if(A->GetYSYM()) dlya=2;
	if(A->GetZSYM()) dlza=2;
	int dlxb=1,dlyb=1,dlzb=1;
	if(B->GetXSYM()) dlxb=2;
	if(B->GetYSYM()) dlyb=2;
	if(B->GetZSYM()) dlzb=2;
	int dlxd=1,dlyd=1,dlzd=1;
	if(A->GetXSYM() && B->GetXSYM()) dlxd=2;
	if(A->GetYSYM() && B->GetYSYM()) dlyd=2;
	if(A->GetZSYM() && B->GetZSYM()) dlzd=2;
	int dLd=1;
	if(A->GetXSYM() && A->GetYSYM() && A->GetZSYM()
		&& B->GetXSYM() && B->GetYSYM() && B->GetZSYM())
		dLd=2;

	Ld=0;
	do{
		delDsumtest=0.0;
		for(lxd=0;lxd<=Ld;lxd+=dlxd){
			for(lyd=0;lyd<=Ld-lxd;lyd+=dlyd){
				lzd=Ld-lxd-lyd;
				if(lzd%dlzd==0){
					gammad=chcalc.Trinomial(lxd,lyd,lzd);
					if(Ld<=LMAXA && lxd%dlxa==0 && lyd%dlya==0 && lzd%dlza==0){
						delD=A->GetElement(lxd,lyd,lzd,ira)/B0;
						Asum+=fabs(A->GetElement(lxd,lyd,lzd,ira));
					}
					else delD=0.0;
					D->SetElement(lxd,lyd,lzd,0,delD);
					delDsum=0.0;
					for(m=0;m<Ld;m+=dLd){
						for(mx=0;mx<=m;mx+=dlxd){
							for(my=0;my<=m-mx;my+=dlyd){
								mz=m-mx-my;
								if(mz%dlzd==0){
									lxb=lxd-mx;
									lyb=lyd-my;
									lzb=lzd-mz;
									Lb=lxb+lyb+lzb;
									if(Lb<=LMAXB && lxb>=0 && lyb>=0 && lzb>=0
									&& lxb%dlxb==0 && lyb%dlyb==0 && lzb%dlzb==0){
										delD=-B->GetElement(lxb,lyb,lzb,irb)
											*D->GetElement(mx,my,mz,0)
											*(chcalc.Trinomial(lxb,lyb,lzb)*chcalc.Trinomial(mx,my,mz))
											/(gammad*B0);
										delDsum+=delD;
									}
								}
							}
						}
					}
					D->IncrementElement(lxd,lyd,lzd,0,delDsum);
					delDsumtest+=fabs(delDsum);
				}
			}
		}
		qtest=0;
		qthresh=EPSILON*Asum/B0;
		if(olddelDsumtest<qthresh && delDsumtest<qthresh) qtest=1;
		olddelDsumtest=delDsumtest;
		Ld+=dLd;

	} while((qtest==0 || Ld<=LMAXA) && Ld<=LMAXD);
	printf("INFO: Sum in DivideAExpArrays satisfied convergence criteria at L=%d\n",
		Ld-1);
	D->SetLMAX(Ld-1);
	//printf("expD=%g\n",AExpand(0.5,0.5,sqrt(0.5),D));

	Detrace(D,0,C,irc);
	delete D;
}

void ArrayCalc::CalcAExpArrayFromXExpArray(CCHArray *X,int irx,CCHArray *A,int ira){

	double snorm,biggy,btest;
	CCHArray *bb,*C,*AA,*oldC;
	bool XSYM,YSYM,ZSYM;
	int n,LMAX,XLMAX,GNMAX;
	XSYM=YSYM=ZSYM=0;
	if(X->GetXSYM()) XSYM=1;
	if(X->GetYSYM()) YSYM=1;
	if(X->GetZSYM()) ZSYM=1;

	GNMAX=48;
	XLMAX=X->GetLMAX();
	LMAX=A->GetLMAX();
	if(LMAX<XLMAX) LMAX=XLMAX;

	C=new CCHArray(LMAX,1,1.0,XSYM,YSYM,ZSYM);
	AA=new CCHArray(LMAX,1,1.0,XSYM,YSYM,ZSYM);
	oldC=new CCHArray(LMAX,1,1.0,XSYM,YSYM,ZSYM);
	bb=new CCHArray(XLMAX,1,1.0,XSYM,YSYM,ZSYM);
	CopyArray(X,irx,bb,0);

	snorm=exp(bb->GetElement(0,0,0,0));
	bb->SetElement(0,0,0,0,0.0);

	AA->ZeroArray(0);
	AA->SetElement(0,0,0,0,1);

	oldC->ZeroArray(0);
	oldC->SetElement(0,0,0,0,1.0);
	n=1;
	biggy=1.0;
	while(n<GNMAX && biggy>1.0E-8){
		if(n%6==0) biggy=0.0;
		MultiplyArrays_Partial(XLMAX,bb,0,XLMAX*(n-1),oldC,0,XLMAX*n,C,0);
		C->ScaleArray(1.0/double(n));
		btest=fabs(C->GetBiggest(0));
		if(btest>biggy) biggy=btest;
		AddArrays(C,0,AA,0,AA,0);
		CopyArray(C,0,oldC,0);
		n+=1;
	}

	//CopyArray(AA,0,A,ira);
	Detrace(AA,0,A,ira);
	A->ScaleArray(snorm,ira);

	delete bb;
	delete C;
	delete oldC;
	delete AA;
}

void ArrayCalc::CalcXExpArrayFromAExpArray(CCHArray *A,int ira,CCHArray *X,int irx){

	double snorm,biggy,btest;
	CCHArray *bb,*C,*XX,*oldC;
	bool XSYM,YSYM,ZSYM;
	int n,LMAX,ALMAX,GNMAX;
	XSYM=YSYM=ZSYM=0;
	if(A->GetXSYM()) XSYM=1;
	if(A->GetYSYM()) YSYM=1;
	if(A->GetZSYM()) ZSYM=1;

	GNMAX=50;
	ALMAX=A->GetLMAX();
	LMAX=X->GetLMAX();
	if(LMAX<ALMAX) LMAX=ALMAX;

	C=new CCHArray(LMAX,1,1.0,XSYM,YSYM,ZSYM);
	XX=new CCHArray(LMAX,1,1.0,XSYM,YSYM,ZSYM);
	oldC=new CCHArray(LMAX,1,1.0,XSYM,YSYM,ZSYM);
	bb=new CCHArray(LMAX,1,1.0,XSYM,YSYM,ZSYM);
	CopyArray(A,ira,bb,0);

	snorm=bb->GetElement(0,0,0,0);
	bb->ScaleArray(1.0/snorm,0);
	bb->SetElement(0,0,0,0,0);

	XX->ZeroArray(0);
	XX->SetElement(0,0,0,0,log(snorm));

	oldC->ZeroArray(0);
	oldC->SetElement(0,0,0,0,1.0);
	n=1;
	biggy=1.0;
	while(n<GNMAX && biggy>1.0E-8){
		if(n%6==0) biggy=0.0;
		//printf("n=%d, biggest bb=%g, biggestoldC=%g\n",n,bb->GetBiggest(0),
		//  oldC->GetBiggest(0));
		MultiplyArrays_Partial(LMAX,bb,0,LMAX*(n-1),oldC,0,LMAX*n,C,0);
		if(n>1) C->ScaleArray(-double(n-1)/double(n),0);
		btest=fabs(C->GetBiggest(0));
		if(btest>biggy) biggy=btest;
		AddArrays(C,0,XX,0,XX,0);
		CopyArray(C,0,oldC,0);
		n+=1;
	}
	//printf("nmax was %d\n",n);

	Detrace(XX,0,X,irx);

	delete bb;
	delete C;
	delete oldC;
	delete XX;
}

void ArrayCalc::CalcAExpArrayFrom3DArray(C3DArray *threed,CCHArray *A){
	double x,y,z,r,cthetamin,phimin,phimax,delx,dely,delz,ex,ey,ez;
	double ctheta,stheta,phi,interpolate;
	const double PI=4.0*atan(1.0);
	int ictheta,iphi,nctheta=60,nphi=60,ir;
	bool XSYM,YSYM,ZSYM;
	CompareArrayParameters(threed,A);
	XSYM=A->GetXSYM();
	YSYM=A->GetYSYM();
	ZSYM=A->GetZSYM();
	if(ZSYM) cthetamin=0.0;
	else cthetamin=-1.0;
	if(XSYM && ZSYM){
		phimin=0.0;
		phimax=0.5*PI;
	}
	else if(XSYM && !YSYM){
		phimin=-0.5*PI;
		phimax=0.5*PI;
	}
	else if(YSYM && !XSYM){
		phimin=-0.5*PI;
		phimax=0.5*PI;
	}
	else{
		phimin=-PI;
		phimax=PI;
	}
	delx=threed->GetDELX();
	dely=threed->GetDELY();
	delz=threed->GetDELZ();
	for(ir=0;ir<A->GetNRADIAL();ir++){
		r=(ir+0.5)*A->GetRADSTEP();
		for(ictheta=0;ictheta<nctheta;ictheta++){
			ctheta=cthetamin+(ictheta+0.5)*(1.0-cthetamin)/double(nctheta);
			stheta=sqrt(1.0-ctheta*ctheta);
			for(iphi=0;iphi<nphi;iphi++){
				phi=phimin+(0.5+iphi)*(phimax-phimin)/double(nphi);
				x=r*stheta*cos(phi);
				y=r*stheta*sin(phi);
				z=r*ctheta;
				if(fabs(x)<delx*threed->GetNXMAX() && fabs(y)<dely*threed->GetNYMAX()
				&& fabs(z)<delz*threed->GetNZMAX()){
					interpolate=threed->GetElement(x,y,z);
					ex=x/r;
					ey=y/r;
					ez=z/r;
					interpolate=interpolate/(double(nctheta*nphi));
					A->IncrementAExpArrayFromE(ex,ey,ez,interpolate,ir);
				}
			}
		}
	}
}

void ArrayCalc::Calc3DArrayFromAExpArray(CCHArray *A,C3DArray *threed){
	int ix,iy,iz,isx,isy,isz,nsx,nsy,nsz;
	double interpolate,x,y,z;
	if(!CompareArrayParameters(threed,A)) exit(1);

	nsx=nsy=nsz=2;
	if(A->GetXSYM()) nsx=1;
	if(A->GetYSYM()) nsy=1;
	if(A->GetZSYM()) nsz=1;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<threed->GetNXMAX();ix++){
			x=(0.5+ix)*threed->GetDELX();
			if(isx==1) x=-x;
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<threed->GetNYMAX();iy++){
					y=(0.5+iy)*threed->GetDELY();
					if(isy==1) y=-y;
					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<threed->GetNZMAX();iz++){
							z=(0.5+iz)*threed->GetDELZ();
							if(isz==1) z=-z;
							interpolate=A->AExpand(x,y,z);
							/*
							r=sqrt(x*x+y*y+z*z);
							if(r<A->GetNRADIAL()*A->GetRADSTEP()){
								irb=int(floor(0.5+r/A->GetRADSTEP()));
								ira=irb-1;
								ex=x/r; ey=y/r; ez=z/r;
								if(ira>=0 && irb<A->GetNRADIAL()){
									valuea=A->AExpand(ex,ey,ez,ira);
									valueb=A->AExpand(ex,ey,ez,irb);
									f=((0.5+irb)*A->GetRADSTEP()-r)/A->GetRADSTEP();
								}
								else if(ira<0){
									valuea=A->AExpand(-ex,-ey,-ez,0);
									valueb=A->AExpand(ex,ey,ez,irb);
									f=((0.5+irb)*A->GetRADSTEP()-r)/A->GetRADSTEP();
								}
								else if(ira==A->GetNRADIAL()){
									valuea=A->AExpand(ex,ey,ez,ira);
									valueb=valuea;
									f=0.0;
								}
								interpolate=f*valuea+(1.0-f)*valueb;
							}*/
							threed->SetElement(isx,ix,isy,iy,isz,iz,interpolate);
						}
					}
				}
			}
		}
	}
}

bool ArrayCalc::CompareArrayParameters(C3DArray *threed,CCHArray *A){
	if(A->GetXSYM()!=threed->GetXSYM() || A->GetYSYM()!=threed->GetYSYM()
	|| A->GetZSYM()!=threed->GetZSYM()){
		printf("X: %d=?%d, Y: %d=?%d Z: %d=?%d\n",
			A->GetXSYM(),threed->GetXSYM(),
			A->GetYSYM(),threed->GetYSYM(),
			A->GetZSYM(),threed->GetZSYM());
		printf("Symmetry mismatch between CHArray and 3DArray!!!\n");
		return 0;
	}
	return 1;
}

bool ArrayCalc::CompareArrayParameters(CCHArray *A,C3DArray*threed){
	if(A->GetXSYM()!=threed->GetXSYM() || A->GetYSYM()!=threed->GetYSYM()
	|| A->GetZSYM()!=threed->GetZSYM()){
		printf("Symmetry mismatch between CHArray and 3DArray!!!\n");
		return 0;
	}
	return 1;
}

bool ArrayCalc::CompareArrayParameters(CCHArray *A,CCHArray *B){
	if(A->GetXSYM()!=B->GetXSYM() || A->GetYSYM()!=B->GetYSYM()
	|| A->GetZSYM()!=B->GetZSYM()){
		printf("Symmetry mismatch between CHArrays!!!\n");
		return 0;
	}
	if(A->GetNRADIAL()!=B->GetNRADIAL()){
		printf("NRADIAL mismatch between CCHArrays!!!\n");
		return 0;
	}
	if(fabs(A->GetRADSTEP())-fabs(B->GetRADSTEP())>1.0E-10){
		printf("RADSTEP mismatch between CCHArrays!!!\n");
		return 0;
	}
	if(A->GetLMAX()!=B->GetLMAX()){
		printf("LMAX mismatch metween CCHArrays!!!\n");
		return 0;
	}
	return 1;
}

bool ArrayCalc::CompareArrayParameters(C3DArray *threeda,C3DArray *threedb){
	bool XSYM,YSYM,ZSYM;
	XSYM=threeda->GetXSYM();
	YSYM=threeda->GetYSYM();
	ZSYM=threeda->GetZSYM();
	if(XSYM!=threedb->GetXSYM() || YSYM!=threedb->GetYSYM()
	|| ZSYM!=threedb->GetZSYM()){
		printf("Symmetry mismatch between 3DArrays!!!\n");
		threeda->PrintPars();
		threedb->PrintPars();
		return 0;
	}
	if(threeda->GetNXMAX()!=threedb->GetNXMAX()
		|| threeda->GetNYMAX()!=threedb->GetNYMAX()
	|| threedb->GetNZMAX()!=threedb->GetNZMAX()){
		printf("NXYZMAX mismatch between 3DArrays!!!\n");
		return 0;
	}
	if(fabs(threeda->GetDELX()-threedb->GetDELX())>1.0E-10
		|| fabs(threeda->GetDELY()-threedb->GetDELY())>1.0E-10
	||	fabs(threeda->GetDELZ()-threedb->GetDELZ())>1.0E-10){
		printf("DELXYZ mismatch between 3DArrays!!!\n");
		return 0;
	}
	return 1;
}

void ArrayCalc::MultiplyArrays(C3DArray *threeda,C3DArray *threedb,C3DArray *threedc){
	int ix,iy,iz,isx,isy,isz,nsx,nsy,nsz;
	if(!CompareArrayParameters(threeda,threedb)) exit(1);
	if(!CompareArrayParameters(threeda,threedc)) exit(1);
	nsx=nsy=nsz=2;
	if(threeda->GetXSYM()) nsx=1;
	if(threeda->GetYSYM()) nsy=1;
	if(threeda->GetZSYM()) nsz=1;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<threeda->GetNXMAX();ix++){
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<threeda->GetNYMAX();iy++){
					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<threeda->GetNZMAX();iz++){
							threedc->SetElement(isx,ix,isy,iy,isz,iz,
								threeda->GetElement(isx,ix,isy,iy,isz,iz)
								*threedb->GetElement(isx,ix,isy,iy,isz,iz));
						}
					}
				}
			}
		}
	}

}

void ArrayCalc::DivideArrays(C3DArray *threeda,C3DArray *threedb,C3DArray *threedc){
	int ix,iy,iz,isx,isy,isz,nsx,nsy,nsz;
	if(!CompareArrayParameters(threeda,threedb)) exit(1);
	if(!CompareArrayParameters(threeda,threedc)) exit(1);
	nsx=nsy=nsz=2;
	if(threeda->GetXSYM()) nsx=1;
	if(threeda->GetYSYM()) nsy=1;
	if(threeda->GetZSYM()) nsz=1;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<threeda->GetNXMAX();ix++){
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<threeda->GetNYMAX();iy++){
					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<threeda->GetNZMAX();iz++){
							threedc->SetElement(isx,ix,isy,iy,isz,iz,
								threeda->GetElement(isx,ix,isy,iy,isz,iz)
								/threedb->GetElement(isx,ix,isy,iy,isz,iz));
						}
					}
				}
			}
		}
	}

}

void ArrayCalc::AddArrays(C3DArray *threeda,C3DArray *threedb,C3DArray *threedc){
	int ix,iy,iz,isx,isy,isz,nsx,nsy,nsz;
	if(!CompareArrayParameters(threeda,threedb)) exit(1);
	if(!CompareArrayParameters(threeda,threedc)) exit(1);
	nsx=nsy=nsz=2;
	if(threeda->GetXSYM()) nsx=1;
	if(threeda->GetYSYM()) nsy=1;
	if(threeda->GetZSYM()) nsz=1;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<threeda->GetNXMAX();ix++){
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<threeda->GetNYMAX();iy++){
					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<threeda->GetNZMAX();iz++){
							threedc->SetElement(isx,ix,isy,iy,isz,iz,
								threeda->GetElement(isx,ix,isy,iy,isz,iz)
								+threedb->GetElement(isx,ix,isy,iy,isz,iz));
						}
					}
				}
			}
		}
	}

}

void ArrayCalc::SubtractArrays(C3DArray *threeda,C3DArray *threedb,C3DArray *threedc){
	int ix,iy,iz,isx,isy,isz,nsx,nsy,nsz;
	if(!CompareArrayParameters(threeda,threedb)) exit(1);
	if(!CompareArrayParameters(threeda,threedc)) exit(1);
	nsx=nsy=nsz=2;
	if(threeda->GetXSYM()) nsx=1;
	if(threeda->GetYSYM()) nsy=1;
	if(threeda->GetZSYM()) nsz=1;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<threeda->GetNXMAX();ix++){
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<threeda->GetNYMAX();iy++){
					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<threeda->GetNZMAX();iz++){
							threedc->SetElement(isx,ix,isy,iy,isz,iz,
								threeda->GetElement(isx,ix,isy,iy,isz,iz)
								-threedb->GetElement(isx,ix,isy,iy,isz,iz));
						}
					}
				}
			}
		}
	}
}

void ArrayCalc::CopyArray(C3DArray *threeda,C3DArray *threedb){
	int ix,iy,iz,isx,isy,isz,nsx,nsy,nsz;
	if(!CompareArrayParameters(threeda,threedb)) exit(1);
	nsx=nsy=nsz=2;
	if(threeda->GetXSYM()) nsx=1;
	if(threeda->GetYSYM()) nsy=1;
	if(threeda->GetZSYM()) nsz=1;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<threeda->GetNXMAX();ix++){
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<threeda->GetNYMAX();iy++){
					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<threeda->GetNZMAX();iz++){
							threedb->SetElement(isx,ix,isy,iy,isz,iz,
								threeda->GetElement(isx,ix,isy,iy,isz,iz));
						}
					}
				}
			}
		}
	}

}

void ArrayCalc::InvertArray(C3DArray *A,C3DArray *B){
	int ix,iy,iz,isx,isy,isz,nsx,nsy,nsz;
	if(!CompareArrayParameters(A,B)) exit(1);
	nsx=nsy=nsz=2;
	if(A->GetXSYM()) nsx=1;
	if(A->GetYSYM()) nsy=1;
	if(A->GetZSYM()) nsz=1;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<A->GetNXMAX();ix++){
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<A->GetNYMAX();iy++){
					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<A->GetNZMAX();iz++){
							B->SetElement(isx,ix,isy,iy,isz,iz,
								1.0/A->GetElement(isx,ix,isy,iy,isz,iz));
						}
					}
				}
			}
		}
	}  
}

void ArrayCalc::InvertArray(CCHArray *A,CCHArray *B){
	int ir;
	if(A->GetNRADIAL()!=B->GetNRADIAL()){
		printf("FATAL: NRADIAL mismatch in ArrayCalc::InvertArrays\n");
		exit(1);
	}
	for(ir=0;ir<A->GetNRADIAL();ir++) InvertArray(A,ir,B,ir);
}

void ArrayCalc::InvertArray(CCHArray *A,int ira,CCHArray *B,int irb){
	if( (A->GetXSYM()!=B->GetXSYM())
		|| (!A->GetYSYM()!=B->GetYSYM())
	|| (!A->GetZSYM()!=B->GetZSYM()) ){
		printf("FATAL: Symmetry mismatch in ArrayCalc::InvertArray\n");
		exit(1);
	}
	CCHArray *C;
	C=new CCHArray(A->GetLMAX(),1,1.0,A->GetXSYM(),
		A->GetYSYM(),A->GetZSYM());
	C->ZeroArray();
	C->IncrementElement(0,0,0,0,1.0);
	DivideArrays(A,ira,C,0,B,irb);
	delete (C);
}
