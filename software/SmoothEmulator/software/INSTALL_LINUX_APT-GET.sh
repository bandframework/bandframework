#!/bin/bash
SMOOTHHOME=`echo ${PWD} | awk '{print substr($1,1,length($1)-9)}'`
MYPYTHON=/usr/bin/python3

oldstring="REPLACEMEWITHSMOOTHHOME"
sed -e "s|${oldstring}|${SMOOTHHOME}|g" scripts/smoothy_emulate.py > scripts/smoothy_emulate.py.tmp
sed -e "s|${oldstring}|${SMOOTHHOME}|g" scripts/smoothy_software_test.sh > scripts/smoothy_software_test.sh.tmp

oldstring="REPLACEMEWITHPYTHON"
sed -e "s|${oldstring}|${MYPYTHON}|g" scripts/smoothy_emulate.py.tmp > ../bin/smoothy_emulate.py
sed -e "s|${oldstring}|${MYPYTHON}|g" scripts/smoothy_software_test.sh.tmp > ../bin/smoothy_software_test.sh

chmod +x ../bin/smoothy_emulate.py
chmod +x ../bin/smoothy_software_test.sh

sudo apt-get update
if dpkg -s g++-11 >/dev/null 2>&1; then
   sudo apt-get upgrade g++-11
else
   sudo apt-get install g++-11
fi
if dpkg -s cmake >/dev/null 2>&1; then
   sudo apt-get upgrade cmake
else
   sudo apt-get install cmake
fi
if dpkg -s libeigen3-dev  >/dev/null 2>&1; then
   sudo apt-get upgrade libeigen3-dev
else
   sudo apt-get install libeigen3-dev
fi
if dpkg -s python3  >/dev/null 2>&1; then
   sudo apt-get upgrade python3
else
   sudo apt-get install python3
fi
if dpkg -s python3-numpy  >/dev/null 2>&1; then
   sudo apt-get upgrade python3-numpy
else
   sudo apt-get install python3-numpy
fi
if dpkg -s python-matplotlib  >/dev/null 2>&1; then
   sudo apt-get upgrade python-matplotlib
else
   sudo apt-get install python-matplotlib
fi
if dpkg -s python3-pybind11  >/dev/null 2>&1; then
   sudo apt-get upgrade python3-pybind11
else
   sudo apt-get install python3-pybind11
fi
if dpkg -s pybind11-dev  >/dev/null 2>&1; then
   sudo apt-get upgrade pybind11-dev
else
   sudo apt-get install pybind11-dev
fi
cmake . -D CMAKE_VERSION=3.2 -D EIGEN3_INCLUDE_DIR=/usr/include/eigen3  -D CMAKE_CXX_COMPILER=/usr/bin/g++-11
make

thisdir=${PWD}
cd pybind_stuff
cmake . -D CMAKE_VERSION=3.2 -D EIGEN3_INCLUDE_DIR=/usr/include/eigen3  -D CMAKE_CXX_COMPILER=/usr/bin/g++-11
make
cd ${thisdir}
exit;
