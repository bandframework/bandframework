cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
set(CMAKE_VERBOSE_MAKEFILE on)
project(smooth)
set(GITHOME_MSU $ENV{GITHOME_MSU})
set (MAIN_PROGRAM_DIR "${GITHOME_MSU}/smooth/local/main_programs")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY $ENV{GITHOME_MSU}/smooth/local/bin)

add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_smooth COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/smooth/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -std=c++17")
set(CMAKE_MODULE_PATH "${MAIN_PROGRAM_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/smooth/software/include)      
link_directories(${GITHOME_MSU}/smooth/software/lib ${GITHOME_MSU}/commonutils/software/lib)

find_library(${GITHOME_MSU}/commonutils/software/lib msu_commonutils)
find_library(${GITHOME_MSU}/smooth/software/lib smooth)

#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/simplex_main.cc)
add_executable(simplex ${MAIN_PROGRAM_DIR}/simplex_main.cc)
add_dependencies(simplex extern_smooth extern_commonutils)
target_link_libraries(simplex ${GITHOME_MSU}/smooth/software/lib/libsmooth.a  ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/smoothy_tune_main.cc)
add_executable(smoothy_tune ${MAIN_PROGRAM_DIR}/smoothy_tune_main.cc)
add_dependencies(smoothy_tune extern_smooth extern_commonutils)
target_link_libraries(smoothy_tune ${GITHOME_MSU}/smooth/software/lib/libsmooth.a  ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/pca_calctransformation_main.cc)
add_executable(pca_calctransformation ${MAIN_PROGRAM_DIR}/pca_calctransformation_main.cc)
add_dependencies(pca_calctransformation extern_smooth extern_commonutils)
target_link_libraries(pca_calctransformation ${GITHOME_MSU}/smooth/software/lib/libsmooth.a  ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/pca_readinfo_calcy_main.cc)
add_executable(pca_readinfo_calcy ${MAIN_PROGRAM_DIR}/pca_readinfo_calcy_main.cc)
add_dependencies(pca_readinfo_calcy extern_smooth extern_commonutils)
target_link_libraries(pca_readinfo_calcy ${GITHOME_MSU}/smooth/software/lib/libsmooth.a  ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
#################################################
set(SOURCE_FILES ${MAIN_PROGRAM_DIR}/smoothy_calcobs_main.cc)
add_executable(smoothy_calcobs ${MAIN_PROGRAM_DIR}/smoothy_calcobs_main.cc)
add_dependencies(smoothy_calcobs extern_smooth extern_commonutils)
target_link_libraries(smoothy_calcobs ${GITHOME_MSU}/smooth/software/lib/libsmooth.a  ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
#################################################

