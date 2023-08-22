#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(RadiaCorrecHallC)

# Setup GNUInstallDirs for installation:
#----------------------------------------------------------------------------
include(GNUInstallDirs)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# test to link fortran function (check)

enable_language(Fortran)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()
    set(CMAKE_Fortran_FLAGS "-ffree-line-length-none  -lowercase -extend_source -g -m64 -Wno-tabs -fno-f2c -DLinux -fno-automatic -finit-local-zero -ffixed-line-length-none -fno-second-underscore -funroll-loops")

#----------------------------------------------------------------------------
# The following lines are required so that the executable can dynamically locate the library dependencies independent of 
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
ENDIF("${isSystemDir}" STREQUAL "-1")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------

#
# Force use of c++11 standard
#
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#


#ROOT path
# commented out the following line so that this procedure will evolve with ROOT:
# set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
#include(${ROOT_USE_FILE})

include_directories(${ROOT_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/src)

include_directories(${PROJECT_SOURCE_DIR}/include)


file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.f  ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
add_executable(a1nd2n radcorrHC.cpp ${sources} ${headers})

target_link_libraries(a1nd2n ${ROOT_LIBRARIES} )


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build N03. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(CLASMODEL_DATAFILES
  f00.ful
  f11.ful
  f12.ful
  f21.ful
  f22.ful
  f90.ful
  f99.ful
  Proton_Data.in
  disa1n_fit_error.dat
  disa1n_fit_param.dat
  disa1p_fit_error.dat
  disa1p_fit_param.dat
  resa1n_fit_final.dat
  resa1p_fit_final.dat
  resa2n_fit_final.dat
  resa2p_fit_final.dat
  resaR_fit_final.dat
  options.ini 
  nhd.tbl1
  nhd.tbl2
  F1F207D2emat.dat
  F2ALLM07.inp1
  F2ALLM07.inp2
  yelas.dat
)

foreach(_script ${CLASMODEL_DATAFILES})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
