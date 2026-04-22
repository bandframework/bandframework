#! /bin/sh
echo "Running script to test Smooth Emulator Software. If successful, pdf images will be generated which should match those in figs/testfigs"
echo "If you are running linux, the comparison figs will automatically be generated if you have okular or evince installed"
echo "If error messages arise and script fails, you need to see where failures occured by viewing output of this script"

   SMOOTH_HOME="REPLACEMEWITHSMOOTHHOME"
   MYPYTHON="REPLACEMEWITHMYPYTHON"
   MYPDFPREVIEWER="REPLACEWITHMYPDFPREVIEWER"
   
   analdir=${PWD}
   mkdir -p figs/figdata
   figsdir=${analdir}/figs
   
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

   cd ${figsdir}
   \cp -f ../smooth_data/MCMC/trace_theta.txt figdata/
   \cp -f ../smooth_data/MCMC/ResolvingPower.txt figdata/
   \cp -f -r ../smooth_data/output_stuff/fullmodel_testdata figdata/
   cd YvsY
   echo 3 | ${MYPYTHON} YvsY.py
   cd ${figsdir}/posterior
   ${MYPYTHON} posterior.py
   cd ${figsdir}/resolvingpower
   ${MYPYTHON} RP.py
   
   cd ${figsdir}
   osname=`uname -s`
   echo --- osname=${osname}  ---
   if [ ${osname} = "Darwin" ]
   then
      open testfigs/YvsY_obs3.pdf
      open testfigs/posterior.pdf
      open testfigs/RP.pdf
   elif [ ${osname} = "Linux" ]
   then
      ${MYPDFVIEWER} testfigs/YvsY_obs3.pdf &
      ${MYPDFVIEWER} testfigs/posterior.pdf &
      ${MYPDFVIEWER} testfigs/RP.pdf &
   else
      echo "Script written for Linux or Mac"
      echo "You can compare figures named figs/testfigs/*_test.pdf to new figures by hand"
   fi
   cd ${analdir}
   echo === FINISHED SOFTWARE TEST COMMANDS ===

