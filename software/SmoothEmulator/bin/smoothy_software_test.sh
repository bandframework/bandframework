#! /bin/sh
echo "Running script to test Smooth Emulator Software. If successful, pdf images will be generated which should match those in figs/testfigs"
echo "If you are running linux, the comparison figs will automatically be generated if you have okular or evince installed"
echo "If error messages arise and script fails, you need to see where failures occured by viewing output of this script"

   SMOOTH_HOME="/Users/scottpratt/SmoothEmulator"
   MYPYTHON=/opt/homebrew/bin/python3
   
   thisdir=${PWD}
   echo -------------------------------------
   rm -f -r smooth_data/FullModelRuns/run*
   rm -f -r smooth_data/FullModelTestingRuns/run*
   rm -f -r smooth_data/fullmodel_testdata/*
   ${SMOOTH_HOME}/bin/smoothy_fakeinfo;
   echo --------- ran smoothy_fakeinfo ---------
   ${SMOOTH_HOME}/bin/smoothy_trainingpoint_optimizer;
   echo --------- ran smoothy_trainingpoint_optimizer ---------
   ${SMOOTH_HOME}/bin/smoothy_fakefullmodel;
   echo --------- ran smoothy_fakefullmodel ---------
   ${SMOOTH_HOME}/bin/smoothy_testattrainingpts;
   echo --------- ran smoothy_testattrainingpts ---------
   ${SMOOTH_HOME}/bin/smoothy_testvsfullmodel;
   echo --------- ran smoothy_testvsfullmodel ---------
   ${SMOOTH_HOME}/bin/smoothy_mcmc;
   echo --------- ran smoothy_mcmc ---------
   mkdir -p figs/figdata
   cd figs/
   \cp -f ../smooth_data/MCMC/trace_theta.txt figdata/
   \cp -f ../smooth_data/MCMC/ResolvingPower.txt figdata/
   \cp -f -r ../smooth_data/output_stuff/fullmodel_testdata figdata/
   cd YvsY
   echo 3 | MYPYTHON YvsY.py
   cd ../posterior
   MYPYTHON posterior.py
   cd ../ResolvingPower
   MYPYTHON RP.py
   cd ${thisdir}
   cd figs/
   osname=`uname -s`
   echo --- osname=${osname}  ---
   if [ ${osname} = "Darwin" ]
   then
      open YvsY/YvsY_obs3.pdf
      open testfigs/YvsY_obs3.pdf
      open posterior/posterior.pdf
      open testfigs/posterior.pdf
      open resolvingpower/RP.pdf
      open testfigs/RP.pdf
   elif [ osname = "Linux" ]
   then
   if command -v okular &> /dev/null
   then
         okular YvsY/YvsY_obs3.pdf &
         okular testfigs/YvsY_obs3.pdf &
         okular posterior/posterior.pdf &
         okular testfigs/posterior.pdf &
         okular resolvingpower/RP.pdf &
         okular testfigs/RP.pdf &
      elif command -v evince &> /dev/null
      then
         evince YvsY/YvsY_obs3.pdf &
         evince testfigs/YvsY_obs3.pdf &
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
   echo === FINISHED SOFTWARE TEST COMMANDS ===

