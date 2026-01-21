#! /bin/bash
echo "Running script to test Smooth Emulator Software. If successful, pdf images will be generated which should match those in figs/testfigs"
echo "If you are running linux, the comparison figs will automatically be generated if you have okular or evince installed"
echo "If error messages arise and script fails, you need to see where failures occured by viewing output of this script"
echo "If you don't have pybind working, the rest of the scripts and the comparisons should function"
exit
echo SMOOTH_HOME = ${SMOOTH_HOME}
echo -------------------------------------
thisdir=${PWD}
cd ${SMOOTH_HOME}/software
cmake .
cmake
echo ------- made C++ programs ----------
cd pybind_stuff
cmake .
make
echo ------- made pybind11 libraries ---------
cd ${thisdir}
rm -f -r smooth_data/FullModelRuns/run*
rm -f -r smooth_data/FullModelTestingRuns/run*
rm -f -r smooth_data/fullmodel_testdata/*
fakeinfo;
echo --------- ran fakeinfo ---------
trainingpoint_optimizer;
echo --------- ran trainingpoint_optimizer ---------
fakefullmodel;
echo --------- ran fakefullmodel ---------
smoothy_testattrainingpts;
echo --------- ran smooth_testattrainingpts ---------
smoothy_testvsfullmodel;
echo --------- ran smoothy_testvsfullmodel ---------
python3 ${SMOOTH_HOME}/software/pybind_stuff/smoothy_emulate.py
echo --------- ran smoothy_emulate.py ---------
smoothy_mcmc;
echo --------- ran smoothy_mcmc ---------
mkdir -p figs/figdata
cd figs/
\cp -f ../smooth_data/MCMC/trace_theta.txt figdata/
\cp -f ../smooth_data/MCMC/ResolvingPower.txt figdata/
\cp -f -r ../smooth_data/output_stuff/fullmodel_testdata figdata/
cd YvsY
echo 3 | python3 YvsY.py
cd ../posterior
python3 posterior.py
cd ../ResolvingPower
python3 RP.py
cd ${thisdir}
cd figs/
osname=`uname -s`
echo --- osname=${osname}  ---
if [ ${osname} = "Darwin" ]
then
   open YvsY/YvsY_obs3.pdf
   open testfigs/YvsY_obs3_test.pdf
   open posterior/posterior.pdf
   open testfigs/posterior_test.pdf
   open resolvingpower/RP.pdf
   open testfigs/RP_test.pdf
elif [ osname = "Linux" ]
then
   if command -v okular &> /dev/null
   then
      okular YvsY/YvsY_obs3.pdf &
      okular testfigs/YvsY_obs3_test.pdf &
      okular posterior/posterior.pdf &
      okular testfigs/posterior_test.pdf &
      okular resolvingpower/RP.pdf &
      okular testfigs/RP_test.pdf &
   elif command -v evince &> /dev/null
   then
      evince YvsY/YvsY_obs3.pdf &
      evince testfigs/YvsY_obs3t.pdf &
      evince posterior/posterior.pdf &
      evince testfigs/posterior.pdf &
      evince resolvingpower/RP.pdf &
      evince testfigs/RP.pdf &
   else
      echo "Need to install okular or evince pdf viewers for this script to work"
      echo "You can compare figures named figs/testfigs/*_test.pdf to new figures by hand"
   fi
else
   echo "Script written for Linux or Mac"
   echo "You can compare figures named figs/testfigs/*_test.pdf to new figures by hand"
fi
cd ${thisdir}
echo FINISHED SOFTWARE TEST COMMANDS

