#!/bin/bash
echo  "=== NOTE: This script will install (if not already installed) homebrew, cmake (version 3.2x or higher), Eigen, python3, numpy, matplotlib, pybind11, and g++-15 ==="
echo "All packages are installed via homebrew and will be placed in /opt/homebrew/... so you may wish to review this script in case you think there might be a conflict with updating any existing software"
echo "=== You can easily change the g++ version (from 15) by editing this script (line 5) ==="
SMOOTH_GCC_VERSION=15
SMOOTHHOME=`echo ${PWD} | awk '{print substr($1,1,length($1)-9)}'`
MYPYTHON=/opt/homebrew/bin/python3
mkdir -p ../bin
nthreads=`nproc`
nthreads=`expr ${nthreads} - 2`

oldstring="REPLACEMEWITHSMOOTHHOME"
sed -e "s|${oldstring}|${SMOOTHHOME}|g" scripts/smoothy_emulate.py > scripts/smoothy_emulate.py.tmp
sed -e "s|${oldstring}|${SMOOTHHOME}|g" scripts/smoothy_software_test.sh > scripts/smoothy_software_test.sh.tmp

oldstring="REPLACEMEWITHMYPYTHON"
sed -e "s|${oldstring}|${MYPYTHON}|g" scripts/smoothy_emulate.py.tmp > ../bin/smoothy_emulate.py
sed -e "s|${oldstring}|${MYPYTHON}|g" scripts/smoothy_software_test.sh.tmp > ../bin/smoothy_software_test.sh

chmod +x ../bin/smoothy_emulate.py
chmod +x ../bin/smoothy_software_test.sh

if ! `command -v /opt/homebrew/bin/brew >/dev/null 2>&1`; then
   echo Installing homebrew
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
fi
if `command -v /opt/homebrew/bin/brew >/dev/null 2>&1`; then
   /opt/homebrew/bin/brew update
   if /opt/homebrew/bin/brew list gcc@${SMOOTH_GCC_VERSION} >/dev/null 2>&1; then
      /opt/homebrew/bin/brew upgrade gcc@${SMOOTH_GCC_VERSION}
   else
      /opt/homebrew/bin/brew install gcc@${SMOOTH_GCC_VERSION}
   fi
   if /opt/homebrew/bin/brew list cmake >/dev/null 2>&1; then
      /opt/homebrew/bin/brew upgrade cmake
   else
      /opt/homebrew/bin/brew install cmake
   fi
   if /opt/homebrew/bin/brew list eigen@5  >/dev/null 2>&1; then
      /opt/homebrew/bin/brew upgrade eigen@5
   else
      /opt/homebrew/bin/brew install eigen@5
   fi
   if /opt/homebrew/bin/brew list python@3  >/dev/null 2>&1; then
      /opt/homebrew/bin/brew upgrade python@3
   else
      /opt/homebrew/bin/brew install python@3
   fi
   if /opt/homebrew/bin/brew list numpy  >/dev/null 2>&1; then
      /opt/homebrew/bin/brew upgrade numpy
   else
      /opt/homebrew/bin/brew install numpy
   fi
   if /opt/homebrew/bin/brew list python-matplotlib  >/dev/null 2>&1; then
      /opt/homebrew/bin/brew upgrade python-matplotlib
   else
      /opt/homebrew/bin/brew install python-matplotlib
   fi
   if /opt/homebrew/bin/brew list pybind11  >/dev/null 2>&1; then
      /opt/homebrew/bin/brew upgrade pybind11
   else
      /opt/homebrew/bin/brew install pybind11
   fi
   cmake . -D CMAKE_VERSION=4.0 -D EIGEN3_INCLUDE_DIR=/opt/homebrew/include/eigen3  -D CMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-${SMOOTH_GCC_VERSION}
   make -j ${nthreads}
fi
   thisdir=${PWD}
   cd pybind_stuff
   cmake . -D CMAKE_VERSION_MAC=4.0 -D EIGEN3_INCLUDE_DIR=/opt/homebrew/include/eigen3  -D CMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-${SMOOTH_GCC_VERSION}
   make -j ${nthreads}
   cd ${thisdir}   
exit;
