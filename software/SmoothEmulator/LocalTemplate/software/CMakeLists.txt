cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(smoothproject)
set(CMAKE_VERBOSE_MAKEFILE on)
project(smooth_mainprograms CXX)
set(SMOOTH_HOME $ENV{SMOOTH_HOME})
set (MAIN_PROGRAM_DIR "MainPrograms")
set (TUTORIAL_PROGRAM_DIR "TutorialPrograms")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "$ENV{PWD}/../bin")
SET( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()

add_custom_target(extern_smooth COMMAND make WORKING_DIRECTORY ${SMOOTH_HOME}/software)

message("Will Compile with " ${CMAKE_CXX_COMPILER})
message("For different compiler enter: cmake -D YOUR_PREFERRED_COMPILER .")
set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++20 -O2")
set(CMAKE_MODULE_PATH "${MAIN_PROGRAM_DIR}/cmake")
find_package(Eigen3 5.0.1 REQUIRED NO_MODULE)
message("^^^^^ eigen3_include_dir is ${EIGEN3_INCLUDE_DIR}")
include_directories(SYSTEM ${SMOOTH_HOME}/software/include $ENV{EIGEN3_INCLUDE_DIR})

link_directories(${SMOOTH_HOME}/software/lib)

find_library(${SMOOTH_HOME}/software/lib msu_smoothutils)
find_library(${SMOOTH_HOME}/software/lib msu_smooth)


#################################################
# trainingpoint_optimizer
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/trainingpoint_optimizer_main.cc)
add_executable(trainingpoint_optimizer ${MAIN_PROGRAM_DIR}/trainingpoint_optimizer_main.cc)

get_target_property(TARGET_INCLUDES trainingpoint_optimizer INCLUDE_DIRECTORIES)
message(include dirs are ${INCLUDE_DIRECTORIES})

add_dependencies(trainingpoint_optimizer extern_smooth)
target_link_libraries(trainingpoint_optimizer ${SMOOTH_HOME}/software/lib/libsmooth.a  ${SMOOTH_HOME}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY})


#################################################
# smoothy_testattrainingpts
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/smoothy_testattrainingpts_main.cc)
add_executable(smoothy_testattrainingpts ${MAIN_PROGRAM_DIR}/smoothy_testattrainingpts_main.cc)
add_dependencies(smoothy_testattrainingpts extern_smooth)
target_link_libraries(smoothy_testattrainingpts ${SMOOTH_HOME}/software/lib/libsmooth.a  ${SMOOTH_HOME}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY})

#################################################
# smoothy_testvsfullmodel
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/smoothy_testvsfullmodel_main.cc)
add_executable(smoothy_testvsfullmodel ${MAIN_PROGRAM_DIR}/smoothy_testvsfullmodel_main.cc)
add_dependencies(smoothy_testvsfullmodel extern_smooth)
target_link_libraries(smoothy_testvsfullmodel ${SMOOTH_HOME}/software/lib/libsmooth.a  ${SMOOTH_HOME}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY})

#################################################
# smoothy_calcobs
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/smoothy_calcobs_main.cc)
add_executable(smoothy_calcobs ${MAIN_PROGRAM_DIR}/smoothy_calcobs_main.cc)
add_dependencies(smoothy_calcobs extern_smooth)
target_link_libraries(smoothy_calcobs ${SMOOTH_HOME}/software/lib/libsmooth.a ${SMOOTH_HOME}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY})

#################################################
# mcmc
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/smoothy_mcmc_main.cc)
add_executable(smoothy_mcmc ${MAIN_PROGRAM_DIR}/smoothy_mcmc_main.cc)
add_dependencies(smoothy_mcmc extern_smooth)
target_link_libraries(smoothy_mcmc ${SMOOTH_HOME}/software/lib/libsmooth.a ${SMOOTH_HOME}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY})

#################################################
# fake full model
#################################################
set(SOURCE_FILES ${TUTORIAL_PROGRAM_DIR}/fakefullmodel.cc)
add_executable(fakefullmodel ${TUTORIAL_PROGRAM_DIR}/fakefullmodel.cc)
add_dependencies(fakefullmodel extern_smooth)
target_link_libraries(fakefullmodel
${SMOOTH_HOME}/software/lib/libsmooth.a ${SMOOTH_HOME}/software/lib/libmsu_smoothutils.a ${EIGEN3_LIBRARY})

#################################################
# create fake info files
#################################################
set(SOURCE_FILES ${TUTORIAL_PROGRAM_DIR}/fakeinfo.cc)
add_executable(fakeinfo ${TUTORIAL_PROGRAM_DIR}/fakeinfo.cc)
add_dependencies(fakeinfo extern_smooth)
target_link_libraries(fakeinfo
${SMOOTH_HOME}/software/lib/libsmooth.a ${SMOOTH_HOME}/software/lib/libmsu_smoothutils.a ${EIGEN3_LIBRARY})
#################################################
