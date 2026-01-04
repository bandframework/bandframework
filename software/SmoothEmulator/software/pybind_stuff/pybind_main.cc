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
   
      .def("GetYOnlyFromThetaPython", static_cast<double (NBandSmooth::CSmoothMaster::*)(int,vector<double> V)>(&NBandSmooth::CSmoothMaster::GetYOnlyFromThetaPython), "Get Y Only From Theta")
   
      .def("GetYSigmaFromThetaPython", static_cast<vector<double> (NBandSmooth::CSmoothMaster::*)(int,vector<double> V)>(&NBandSmooth::CSmoothMaster::GetYSigmaFromThetaPython), "Get YSigma From Theta")
   
      .def("GetYOnlyFromXPython", static_cast<double (NBandSmooth::CSmoothMaster::*)(int,vector<double> V)>(&NBandSmooth::CSmoothMaster::GetYOnlyFromXPython), "Get Y Only From X")
   
      .def("GetYSigmaFromXPython", static_cast<vector<double> (NBandSmooth::CSmoothMaster::*)(int,vector<double> V)>(&NBandSmooth::CSmoothMaster::GetYSigmaFromXPython), "Get YSigma From X")
   
      .def("GetThetaFromX", static_cast<vector<double> (NBandSmooth::CSmoothMaster::*)(vector<double> V)>(&NBandSmooth::CSmoothMaster::GetThetaFromX), "Get ThetaFromX")
   
      .def("GetXFromTheta", static_cast<vector<double> (NBandSmooth::CSmoothMaster::*)(vector<double> V)>(&NBandSmooth::CSmoothMaster::GetXFromTheta), "Get XFromTheta")
   
      .def("TuneAllY",&NBandSmooth::CSmoothMaster::TuneAllY);
   
}
