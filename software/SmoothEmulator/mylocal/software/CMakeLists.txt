cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
set(CMAKE_VERBOSE_MAKEFILE on)
project(smooth_mainprograms CXX)
#set(GITHOME_BAND_SMOOTH $ENV{GITHOME_BAND_SMOOTH})
set(GITHOME_BAND_SMOOTH "$ENV{PWD}/../..")
set (MAIN_PROGRAM_DIR "main_programs")
set (FAKEFULLMODEL_DIR "fakefullmodels")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "$ENV{PWD}/../bin")

add_custom_target(extern_smooth COMMAND make WORKING_DIRECTORY ${GITHOME_BAND_SMOOTH}/software)

message("Will Compile with " ${CMAKE_CXX_COMPILER})
message("For different compiler enter: cmake -D YOUR_PREFERRED_COMPILER .")
set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")
set(CMAKE_MODULE_PATH "${MAIN_PROGRAM_DIR}/cmake")
find_package(Eigen3 REQUIRED)
message("eigen3_include_dir is ${EIGEN3_INCLUDE_DIR}")
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR} ${GITHOME_BAND_SMOOTH}/software/include $ENV{EIGEN3_INCLUDE_DIR})      
link_directories(${GITHOME_BAND_SMOOTH}/software/lib)

find_library(${GITHOME_BAND_SMOOTH}/software/lib msu_smoothutils)
find_library(${GITHOME_BAND_SMOOTH}/software/lib msu_smooth)

#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/simplex_main.cc)
add_executable(simplex ${MAIN_PROGRAM_DIR}/simplex_main.cc)
add_dependencies(simplex extern_smooth)
target_link_libraries(simplex ${GITHOME_BAND_SMOOTH}/software/lib/libsmooth.a  ${GITHOME_BAND_SMOOTH}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY})
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/smoothy_tune_main.cc)
add_executable(smoothy_tune ${MAIN_PROGRAM_DIR}/smoothy_tune_main.cc)
add_dependencies(smoothy_tune extern_smooth)
target_link_libraries(smoothy_tune ${GITHOME_BAND_SMOOTH}/software/lib/libsmooth.a  ${GITHOME_BAND_SMOOTH}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY})
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/smoothy_testattrainingpts_main.cc)
add_executable(smoothy_testattrainingpts ${MAIN_PROGRAM_DIR}/smoothy_testattrainingpts_main.cc)
add_dependencies(smoothy_testattrainingpts extern_smooth)
target_link_libraries(smoothy_testattrainingpts ${GITHOME_BAND_SMOOTH}/software/lib/libsmooth.a  ${GITHOME_BAND_SMOOTH}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY})
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/smoothy_testvsfullmodel_main.cc)
add_executable(smoothy_testvsfullmodel ${MAIN_PROGRAM_DIR}/smoothy_testvsfullmodel_main.cc)
add_dependencies(smoothy_testvsfullmodel extern_smooth)
target_link_libraries(smoothy_testvsfullmodel ${GITHOME_BAND_SMOOTH}/software/lib/libsmooth.a  ${GITHOME_BAND_SMOOTH}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY})
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/pca_calctransformation_main.cc)
add_executable(pca_calctransformation ${MAIN_PROGRAM_DIR}/pca_calctransformation_main.cc)
add_dependencies(pca_calctransformation extern_smooth)
target_link_libraries(pca_calctransformation ${GITHOME_BAND_SMOOTH}/software/lib/libsmooth.a ${GITHOME_BAND_SMOOTH}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY})
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/pca_readinfo_calcy_main.cc)
add_executable(pca_readinfo_calcy ${MAIN_PROGRAM_DIR}/pca_readinfo_calcy_main.cc)
add_dependencies(pca_readinfo_calcy extern_smooth)
target_link_libraries(pca_readinfo_calcy ${GITHOME_BAND_SMOOTH}/software/lib/libsmooth.a ${GITHOME_BAND_SMOOTH}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY})
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/smoothy_calcobs_main.cc)
add_executable(smoothy_calcobs ${MAIN_PROGRAM_DIR}/smoothy_calcobs_main.cc)
add_dependencies(smoothy_calcobs extern_smooth)
target_link_libraries(smoothy_calcobs ${GITHOME_BAND_SMOOTH}/software/lib/libsmooth.a ${GITHOME_BAND_SMOOTH}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY})
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/mcmc_main.cc)
add_executable(mcmc ${MAIN_PROGRAM_DIR}/mcmc_main.cc)
add_dependencies(mcmc extern_smooth)
target_link_libraries(mcmc ${GITHOME_BAND_SMOOTH}/software/lib/libsmooth.a ${GITHOME_BAND_SMOOTH}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY})
#################################################
#
# Some fake full models (for tutorial or emulator-testing purposes
#################################################
# fakesmooth
#################################################
set(SOURCE_FILES ${FAKEFULLMODEL_DIR}/fakesmooth.cc)
add_executable(fakesmooth ${FAKEFULLMODEL_DIR}/fakesmooth.cc)
add_dependencies(fakesmooth extern_smooth)
target_link_libraries(fakesmooth 
${GITHOME_BAND_SMOOTH}/software/lib/libsmooth.a ${GITHOME_BAND_SMOOTH}/software/lib/libmsu_smoothutils.a ${EIGEN3_LIBRARY})
#################################################
# fakerhic
#################################################
set(SOURCE_FILES ${FAKEFULLMODEL_DIR}/fakerhic.cc)
add_executable(fakerhic ${FAKEFULLMODEL_DIR}/fakerhic.cc)
add_dependencies(fakerhic extern_smooth)
target_link_libraries(fakerhic ${GITHOME_BAND_SMOOTH}/software/lib/libmsu_smoothutils.a
 ${EIGEN3_LIBRARY})
#################################################
