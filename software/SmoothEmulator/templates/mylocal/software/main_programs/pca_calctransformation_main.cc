#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smooth/pca.h"

using namespace std;

int main() {
	CparameterMap parmap;
  NBandSmooth::CPCA *pca = new NBandSmooth::CPCA(&parmap);

  pca->CalcTransformationInfo();

  //pca->ReadPCATransformationInfo();


}
