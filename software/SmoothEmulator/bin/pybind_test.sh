#! /bin/sh
echo "This script requires pybind11 to be installed"
echo "It also assumes that smooth_test.sh ran successfully"
#export PYTHONPATH=/opt/homebrew/bin
echo SMOOTH_HOME = ${SMOOTH_HOME}
echo -------------------------------------
thisdir=${PWD}
cd ${SMOOTH_HOME}/software/pybind_stuff
cmake .
make
echo ------- made pybind11 libraries ---------
cd ${thisdir}
/opt/homebrew/bin/python3 ${SMOOTH_HOME}/software/pybind_stuff/smoothy_emulate.py
echo --------- ran script smoothy_emulate.py ---------

