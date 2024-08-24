#include <pybind11/pybind11.h>
namespace py = pybind11;
using namespace NBandSmooth;

#include "msu_smooth/master.h"

PYBIND11_MODULE(msu_smooth,m){
		py::class_<NBandSmooth::CSmoothMaster>(m,"smoothmaster")
			.def(py::init<>())
			.def("AddXY",&NBandSmooth::CSmoothMaster::addxy);
}


//PYBIND11_MODULE(msu_smooth,m){
//	py::class_<NBandSmooth::CSmoothMaster>(m,"smoothmaster")
//	.def("GetY",[](NBandSmooth::CSmoothMaster *smoothmaster,int iY,vector<double> theta){
//		double Y=smoothmaster->GetYOnly(iY,theta);
//		return Y;
//	});
//};