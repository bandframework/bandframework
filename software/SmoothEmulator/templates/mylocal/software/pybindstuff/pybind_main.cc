#include "msu_smooth/master.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using namespace std;
PYBIND11_MODULE(emulator_smooth, m) {

    py::class_<NBandSmooth::CSmoothMaster>(m, "emulator_smooth")

        .def(py::init<>())
					
				.def("GetNPars", static_cast<int (NBandSmooth::CSmoothMaster::*)()>(&NBandSmooth::CSmoothMaster::GetNPars), "Get NPars")
				
				.def("GetNObs", static_cast<int (NBandSmooth::CSmoothMaster::*)()>(&NBandSmooth::CSmoothMaster::GetNObs), "Get NObs")

        .def("GetYOnlyPython", static_cast<double (NBandSmooth::CSmoothMaster::*)(int,vector<double> V)>(&NBandSmooth::CSmoothMaster::GetYOnlyPython), "Get Y Only, I guess")
					
				.def("GetYSigmaPython", static_cast<vector<double> (NBandSmooth::CSmoothMaster::*)(int,vector<double> V)>(&NBandSmooth::CSmoothMaster::GetYSigmaPython), "Get YSigma, I guess")
		
		    .def("GetThetaFromX", static_cast<vector<double> (NBandSmooth::CSmoothMaster::*)(vector<double> V)>(&NBandSmooth::CSmoothMaster::GetThetaFromX), "Get ThetaFromX, I guess")
					
			  .def("GetXhetaFromTheta", static_cast<vector<double> (NBandSmooth::CSmoothMaster::*)(vector<double> V)>(&NBandSmooth::CSmoothMaster::GetXFromTheta), "Get XFromTheta, I guess")
		
				.def("TuneAllY",&NBandSmooth::CSmoothMaster::TuneAllY);

}
