#!/bin/bash
echo  "=== NOTE: This script will install (if not already installed) cmake (version 4.x or higher), Eigen, python3, numpy, matplotlib, pybind11, okular and g++-11 ==="
echo "=== You can easily change the g++ version by editing this script (line 5) or the pdf previewer from okular to evince (line 8) ==="
echo "=== The apt-get package installer will request su priveleges to install packages in /usr/... so you may wish to review this script in case you think there might be a conflict with updating any existing software"
SMOOTH_GCC_VERSION=11
SMOOTHHOME=`echo ${PWD} | awk '{print substr($1,1,length($1)-9)}'`
MYPYTHON=/usr/bin/python3
MYPDFPREVIEWER=okular
#PDFPREVIEWER=evince
mkdir -p ../bin
nthreads=`nproc`

oldstring="REPLACEMEWITHSMOOTHHOME"
sed -e "s|${oldstring}|${SMOOTHHOME}|g" scripts/smoothy_emulate.py > scripts/smoothy_emulate.py.tmp
sed -e "s|${oldstring}|${SMOOTHHOME}|g" scripts/smoothy_software_test.sh > scripts/smoothy_software_test.sh.tmp

oldstring="REPLACEMEWITHMYPDFPREVIEWER"
sed -e "s|${oldstring}|${MYPDFPREVIEWER}|g" scripts/smoothy_software_test.sh.tmp > scripts/smoothy_software_test.sh.tmp2

oldstring="REPLACEMEWITHMYPYTHON"
sed -e "s|${oldstring}|${MYPYTHON}|g" scripts/smoothy_emulate.py.tmp > ../bin/smoothy_emulate.py
sed -e "s|${oldstring}|${MYPYTHON}|g" scripts/smoothy_software_test.sh.tmp2 > ../bin/smoothy_software_test.sh

chmod +x ../bin/smoothy_emulate.py
chmod +x ../bin/smoothy_software_test.sh

sudo apt-get update
if dpkg -s g++-${SMOOTH_GCC_VERSION} >/dev/null 2>&1; then
   sudo apt-get upgrade g++-${SMOOTH_GCC_VERSION}
else
   sudo apt-get install g++-${SMOOTH_GCC_VERSION}
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
if dpkg -s python3-matplotlib  >/dev/null 2>&1; then
   sudo apt-get upgrade python3-matplotlib
else
   sudo apt-get install python3-matplotlib
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
if dpkg -s ${MYPDFPREVIEWER}  >/dev/null 2>&1; then
   sudo ${MYPDFPREVIEWER} upgrade ${MYPDFPREVIEWER}
else
   sudo apt-get install ${MYPDFPREVIEWER}
fi
cmake . -D CMAKE_VERSION=3.2 -D EIGEN3_INCLUDE_DIR=/usr/include/eigen3  -D CMAKE_CXX_COMPILER=/usr/bin/g++-${SMOOTH_GCC_VERSION}
make -j ${nthreads}

thisdir=${PWD}
cd pybind_stuff
cmake . -D CMAKE_VERSION=3.2 -D EIGEN3_INCLUDE_DIR=/usr/include/eigen3  -D CMAKE_CXX_COMPILER=/usr/bin/g++-${SMOOTH_GCC_VERSION}
make -j ${nthreads}
cd ${thisdir}
exit;
