#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smooth/pca.h"

using namespace std;

int main(){
  NBandSmooth::CPCA *pca = new NBandSmooth::CPCA();
  pca->CalcTransformationInfo();
}
