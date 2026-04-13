#!/bin/bash
SMOOTHHOME=`echo ${PWD} | awk '{print substr($1,1,length($1)-9)}'`
MYPYTHON=/opt/homebrew/bin/python3

oldstring="REPLACEMEWITHSMOOTHHOME"
sed -e "s|${oldstring}|${SMOOTHHOME}|g" scripts/smoothy_emulate.py > scripts/smoothy_emulate.py.tmp
sed -e "s|${oldstring}|${SMOOTHHOME}|g" scripts/smoothy_software_test.sh > scripts/smoothy_software_test.sh.tmp

oldstring="REPLACEMEWITHPYTHON"
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
   if /opt/homebrew/bin/brew list gcc@15 >/dev/null 2>&1; then
      /opt/homebrew/bin/brew upgrade gcc@15
   else
      /opt/homebrew/bin/brew install gcc@15
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
   cmake . -D EIGEN3_INCLUDE_DIR=/opt/homebrew/include/eigen3  -D CMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15
   make
fi
   thisdir=${PWD}
   cd pybind_stuff
   cmake . -D EIGEN3_INCLUDE_DIR=/opt/homebrew/include/eigen3  -D CMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15
   make
   cd ${thisdir}   
exit;
