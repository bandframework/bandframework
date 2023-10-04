#ifndef __INCLUDE_ARRAYS_H__
#define __INCLUDE_ARRAYS_H__

#include "commondefs.h"

using namespace std;

class CCHArray{
public:
	CCHArray(string arrayparsfilename);
	CCHArray(int LMAXset,int NRADIALset,double RADSTEPset);
	CCHArray(int LMAXset,int NRADIALset,double RADSTEPset,
	bool XSYMset,bool YSYMset,bool ZSYMset);
	~CCHArray();

	int GetLMAX();
	int GetNRADIAL();
	double GetRADSTEP();
	void SetLMAX(int LMAXset);
	void SetRADSTEP(double RADSTEPset);
	bool GetXSYM();
	bool GetYSYM();
	bool GetZSYM();
	void PrintPars();

	void ZeroArray();
	void ZeroArray(int lx,int ly,int lz);
	void ZeroArray(int ir);
	void ZeroArray_Partial(int LMAX_Partial);
	void ZeroArray_Partial(int LMAX_Partial,int ir);
	void ScaleArray(double scalefactor);
	void ScaleArray(double scalefactor,int ir);

	double GetElement(int lx,int ly,int lz,int ir);
	double GetElement(int lx,int ly,int lz,double r);
	void SetElement(int lx,int ly,int lz,int ir,double Element);
	void SetElement(int lx,int ly,int lz,double r,double Element);
	void IncrementElement(int lx,int ly,int lz,int ir,double increment);
	void IncrementElement(int lx,int ly,int lz,double r,double increment);

	void PrintArrayFixedIR(int ir);
	void PrintArrayFixedIR(int LMAXPrint,int ir);
	void Print(int lx,int ly,int lz);
	void PrintProjections();
	void WriteProjections(string filename);
	void GetProjections(double **A); // A[4][NRADIAL], A[0] is angle-averaged...
	//FillRemainder functions use identity
	//  A_{lx+2,ly,lz}+A_{lx,ly+2,lz}+A_{lx,ly,lz+2}=0
	// to find entire array if (2L+1) values for lx=0,1 are known
	void PrintMoments();

	double GetBiggest(int ir);

	// Read and Write to a given directory, "AX" suffix means only lx=0,1 elements
	// will be read/written. To write only up to WLMAX, use WriteShort()
	void ReadAX(string dirname);
	void WriteAX(string dirname);
	void ReadAllA(string dirname);
	void WriteAllA(string dirname);
	void WriteShort(string filename,int WLMAX);

	void IncrementAExpArray(double x,double y,double z,double weight);
	void IncrementAExpArrayFromE(double ex,double ey,double ez,
	double weight,int ir);
	void AltIncrementAExpArrayFromE(double ex,double ey,double ez,
	double weight,int ir);
	void AltAltIncrementAExpArrayFromE(double ex,double ey,double ez,
	double weight,int ir);
	void IncrementMArrayFromE(double ex,double ey,double ez,double weight,int ir);
	void IncrementAExpArrayFromThetaPhi(double theta,double phi,
	double weight,int ir);
	void IncrementMArrayFromThetaPhi(double theta,double phi,
	double weight,int ir);

	// Returns < ex^lx ey^ly ez^lz > where (ex,ey,ez) is unit vector
	double GetMElementFromAExpArray(int lx,int ly,int lz,int ir);
	// Given array, < ex^lx ey^ly ez^lz >, this returns the CH expansion 
	// element 
	double GetAExpElementFromMArray(int lx,int ly,int lz,int ir);

	void FillRemainderX(int ir);
	void FillRemainderY(int ir);
	void FillRemainderZ(int ir);
	void FillRemainderX();
	void FillRemainderY();
	void FillRemainderZ();

	// For CH Expansion array, this returns
	// F(Omega)=\sum_{lx,ly,lz} (L!/(lx!ly!lz!))
	// F_{lx,ly,lz} e_x^lx e_y^ly e_z^lz
	double AExpand(double ex,double ey,double ez,int ir);
	double AExpand(double x,double y,double z);
	double AExpand(double theta,double phi,int ir);

	void Detrace(int ir);
	void Detrace();

	void Randomize(double mag,int ir);
	void RandomizeA(double mag,int ir);
	void Randomize(double mag);
	void RandomizeA(double mag); 
	void RandomizeA_Gaussian(double mag,int ir);
	void RandomizeA_Gaussian(double mag);

	static CCHCalc *chcalc;

private:
	static Crandy *randy;
	bool XSYM,YSYM,ZSYM;
	double RADSTEP;
	int NRADIAL;
	int LMAX;
	double ****A;
	int dlx,dly,dlz;
	void CreateArray();
	void RandomInit(int iseed);

};

class C3DArray{
public:
	C3DArray(string arrayparsfilename);
	C3DArray(int NXYZMAX,double DELXYZ,bool XSYM,bool YSYM,bool ZSYM);
	C3DArray(int NXMAX,double DELX,int NYMAX,double DELY,int NZMAX,double DELZ,
	bool XSYM,bool YSYM,bool ZSYM);
	~C3DArray();

	int GetNXMAX();
	int GetNYMAX();
	int GetNZMAX();
	double GetDELX();
	double GetDELY();
	double GetDELZ();
	bool GetXSYM();
	bool GetYSYM();
	bool GetZSYM();
	void PrintPars();

	double GetElement(double x,double y,double z);
	double GetElement_NoInterpolation(double x,double y,double z);
	double GetElement(int isx,int ix,int isy,int iy,int isz,int iz);
	void SetElement(int isx,int ix,int isy,int iy,int isz,int iz,double value);
	void IncrementElement(int isx,int ix,int isy,int iy,int isz,int iz,
	double value);
	void SetElement(double x,double y,double z,double value);
	void IncrementElement(double x,double y,double z,double increment);
	void CalcMoments(double roff[3],double r2[3][3]);
	void PrintMoments();

	void ZeroArray();
	void MakeConstant(double c);
	void ScaleArray(double scalefactor);

	void Randomize(double c);
	void RandomizeGaussian(double c);
	void DivideByArray(C3DArray *threed_denom);

	void PrintArray();
	void PrintProjections();
	void WriteProjections(string filename);
	double GetBiggest();

	void ReadArray(string dirname);
	void WriteArray(string dirname);
	void WriteArray_PHENIX(string filename);

private:
	bool XSYM,YSYM,ZSYM;
	double DELX,DELY,DELZ;
	int NXMAX,NYMAX,NZMAX;
	double ******F;
	void ReadPars(string arrayparsfilename);
	void CreateArray();
	void DeleteArray();
	static Crandy *randy;
};

class CYlmArray{
public:
	CYlmArray(int LMAXset,int NRADIALset);
	~CYlmArray();
	int GetLMAX();
	complex<double> GetElement(int L,int M,int ir);
	void SetElement(int L,int M,int ir,complex<double>);
	void IncrementElement(int L,int M,int ir,complex<double> increment);
	void ScaleArray(double scalefactor);
	void ScaleArray(double scalefactor,int ir);
	void ZeroArray();
	void ZeroArray(int ir);
	void PrintArrayFixedIR(int ir);
	void PrintArrayFixedIR(int LMAXPrint,int ir);
private:
	int NRADIAL;
	int LMAX;
	complex<double> ***ylm;
};

class CMCList{
public:
	CMCList(int nmcset);
	~CMCList();
	void Resize(int nmcset);
	int GetNMC();
	void SetR(int imc,double t,double x,double y,double z);
	void SetR(int imc,double *r1);
	void SetNorm(double normset);
	double GetNorm();
	double *GetR(int imc);
	void PrintMoments(double Rmax);
private:
	int nmc;
	double norm;
	double **r;
};

class CMCPRList{
public:
	CMCPRList(int nmcset);
	~CMCPRList();
	void Resize(int nmcset);
	int GetNMC();
	void SetPR(int imc,double p0,double px,double py,double pz,double t,double x,double y,double z);
	void SetPR(int imc,double *p1,double *r1);
	void SetNorm(double normset);
	double GetNorm();
	double *GetR(int imc);
	double *GetP(int imc);
	void PrintMoments(double Rmax);
private:
	int nmc;
	double norm;
	double **r;
	double **p;
};

namespace ArrayCalc{

	void CalcMArrayFromAExpArray(CCHArray *A,CCHArray *M);
	void CalcMArrayFromAExpArray(CCHArray *A,int ira,CCHArray *M,int irm);
	void CalcAExpArrayFromMArray(CCHArray *M,CCHArray *A);
	void CalcAExpArrayFromMArray(CCHArray *M,int irm,CCHArray *A,int ira);

	void AddArrays(CCHArray *A,CCHArray *B,CCHArray *C);
	void SubtractArrays(CCHArray *A,CCHArray *B,CCHArray *C);
	void AddArrays(CCHArray *A,int ira,CCHArray *B,int irb,CCHArray *C,int irc);
	void SubtractArrays(CCHArray *A,int ira,CCHArray *B,int irb,CCHArray *C,int irc);
	// If A(Omega)=B(Omega)*C(Omega), this finds C in terms of A and B
	void DivideArrays(CCHArray *A,CCHArray *B,CCHArray *C);
	void DivideArrays(CCHArray *A,int ira,CCHArray *B,int irb,
	CCHArray *C,int irc);
	// If C(Omega)=A(Omega)*B(Omega), this finds A in terms of A and B
	void MultiplyArrays(CCHArray *A,CCHArray *B,CCHArray *C);
	void MultiplyArrays(CCHArray *A,int ira,CCHArray *B,
	int irb,CCHArray *C,int irc);
	// If you know array is zero up to given Ls, or don't care to detrace,
	// this can save time
	void MultiplyArrays_Partial(int LMAXA,CCHArray *A,int ira,
	int LMAXB,CCHArray *B,int irb,
	int LMAXC,CCHArray *C,int irc);

	void CopyArray(CCHArray *A,CCHArray *B);
	void CopyArray(CCHArray *A,int ira,CCHArray *B,int irb);

	void CalcYlmExpArrayFromAExpArray(CCHArray *A,int ir,
	CYlmArray *YlmArray,int irlm);
	void CalcAExpArrayFromYlmExpArray(CYlmArray *YlmArray,int irlm,
	CCHArray *A,int ira);

	// If one has an angular function
	// F=\sum_{\vec\ell} M_{\vec\ell} e_x^{\ell_x}e_y^{\ell_y}e_z^{\ell_z}
	// and wants to find expansion coefficients A which give the same answer
	// but satisfy traceless condition, this finds the array A
	void Detrace(CCHArray *M,CCHArray *A);
	void Detrace(CCHArray *M,int irm,CCHArray *A,int ira);

	// XExpArray express functions in form exp( X_ell * nhat^ell ) 
	void CalcAExpArrayFromXExpArray(CCHArray *X,CCHArray *A);
	void CalcAExpArrayFromXExpArray(CCHArray *X,int irx,CCHArray *A,int ira);
	void CalcXExpArrayFromAExpArray(CCHArray *A,CCHArray *X);
	void CalcXExpArrayFromAExpArray(CCHArray *A,int ira,CCHArray *X,int irx);


	void CalcAExpArrayFrom3DArray(C3DArray *threedarray,CCHArray *A);
	void Calc3DArrayFromAExpArray(CCHArray *A,C3DArray *threed);

	void MultiplyArrays(C3DArray *A,C3DArray *B,C3DArray *C);
	void DivideArrays(C3DArray *A,C3DArray *B,C3DArray *C);
	void AddArrays(C3DArray *A,C3DArray *B,C3DArray *C);
	void SubtractArrays(C3DArray *A,C3DArray *B,C3DArray *C);
	void CopyArray(C3DArray *A,C3DArray *B);

	void InvertArray(C3DArray *A,C3DArray *B);
	void InvertArray(CCHArray *A,int ira,CCHArray *B,int irb);
	void InvertArray(CCHArray *A,CCHArray *B);

	/* These check that two arrays have the same parameters
	For arrays of different types (e.g., 3D and Cart.Harmonic) it checks that
	they have the same reflection symmetries */
	bool CompareArrayParameters(C3DArray *threed,CCHArray *A);
	bool CompareArrayParameters(CCHArray *A,C3DArray *threed);
	bool CompareArrayParameters(CCHArray *A,CCHArray *B);
	bool CompareArrayParameters(C3DArray *threeda,C3DArray *threedb);
};

#endif

