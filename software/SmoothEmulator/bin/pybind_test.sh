#! /bin/sh
echo "This script requires pybind11 to be installed"
echo "It also assumes that smoothy_test.sh ran successfully"
#export PYTHONPATH=/opt/homebrew/bin
case $# in
0)
	echo "Usage: pybind_test.sh  SMOOTH_HOME_PATH";
   echo "PATH is to home directory of Smooth Installation, e.g. ../bandframework/software/SmoothEmulator"
	exit 1 ;;
*)
   SMOOTH_HOME=$1
   echo SMOOTH_HOME = ${SMOOTH_HOME}
   PATH=${PATH}:${SMOOTH_HOME}/bin
   echo -------------------------------------
   /opt/homebrew/bin/python3 ${SMOOTH_HOME}/software/pybind_stuff/smoothy_emulate.py
   echo === FINISHED TEST OF PYBIND11 FUNCTIONALITY ===
esac
