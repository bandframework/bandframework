#ifndef __PCA_H__
#define __PCA_H__
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <filesystem>
#include "msu_commonutils/parametermap.h"
#include "msu_smooth/master.h"
using namespace NMSUPratt;
namespace NBandSmooth{
	
	class PCA{
	public:

		PCA(string parameter_filename);
		int nruns,Nobs;
		Eigen::MatrixXd eigvecs;
		Eigen::VectorXd eigvals;
		vector<vector<double>> Y;
		vector<double> SigmaY,Ybar;
		vector<vector<double>> SigmaY_emulator;
		vector<int> NTrainingList;
		string modelruns_dirname;
		CObservableInfo *observable_info;

		void CalcTransformationInfo();
		void WriteTransformationInfo();
		void ReadTransformationInfo();
	
		void TransformZtoY(vector<double> &Z,vector<double> &SigmaZ_emulator,
		vector<double> &Y,vector<double> &SigmaY_emulator);
	
		void TransformYtoZ(vector<double> &Y,vector<double> &SigmaY_emulator,
		vector<double> &Z,vector<double> &SigmaZ_emulator);
	
	};

}

#endif
