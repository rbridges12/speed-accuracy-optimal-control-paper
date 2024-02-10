###########################################
## For building with the rest of OpenSim ##
###########################################

cmake_minimum_required(VERSION 2.6)

set( CMAKE_C_COMPILER "/share/software/user/open/gcc/6.3.0/bin/gcc")
set( CMAKE_CXX_COMPILER "/share/software/user/open/gcc/6.3.0/bin/g++")

#set(CMAKE_INCLUDE_PATH=$CPATH)
#set(CMAKE_LIBRARY_PATH=$LIBRARY_PATH)

# Define project
PROJECT (OptimizationExample)

message(status "right cmake file")

#list(APPEND OpenSim_SOURCE_DIR "/XXX/opensim/opensim-3.3.0/source")
SET(OpenSim_SOURCE_DIR "XXX/OpenSim33-source")

INCLUDE_DIRECTORIES(${OpenSim_SOURCE_DIR} ${OpenSim_SOURCE_DIR}/Vendors)

SET(OPENSIM_INSTALL_DIR "/home/users/XXX/opensim33" 
		CACHE PATH "Top-level directory of OpenSim install")

find_package(MPI REQUIRED)
find_package(OpenCV REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${OpenCV_INCLUDE_DIRS})

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#elseif(COMPILER_SUPPORTS_CXX0X)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


SET(BUILD_USING_OTHER_LAPACK "/share/software/user/open/openblas/0.2.19/lib/liblapack.so;/share/software/user/open/openblas/0.2.19/lib/libblas.so" CACHE STRING
 "If you have your own Lapack, put its library basename here. Default is to use lapack.")

IF(BUILD_USING_OTHER_LAPACK)
    SET(LAPACK_BEING_USED ${BUILD_USING_OTHER_LAPACK} CACHE STRING
"Basename of the actual Lapack library we're depending on; can't change here; see variable BUILD_USING_OTHER_LAPACK." FORCE)
ENDIF(BUILD_USING_OTHER_LAPACK)


# Change name of build target
SET(TARGET optimizationExample CACHE TYPE STRING)
FILE(GLOB SOURCE_FILES *.cpp *.c)
FILE(GLOB INCLUDE_FILES *.h)
SET(SOURCE ${SOURCE_FILES})
SET_SOURCE_FILES_PROPERTIES( ${CMAKE_SOURCE_DIR}/cmaes.c PROPERTIES LANGUAGE CXX )


# To add Debug feature add ";Debug" after Release on the line below
SET(CMAKE_CONFIGURATION_TYPES "RelWithDebInfo;Release;Debug" 
      CACHE STRING "Semicolon separated list of supported configuration types, only supports Debug, Release, MinSizeRel, and RelWithDebInfo, anything else will be ignored." FORCE )
	
# Location of headers
#SET(SIMTK_HEADERS_DIR ${OPENSIM_INSTALL_DIR}/sdk/include/SimTK/include)
SET(SIMTK_HEADERS_DIR ${OPENSIM_INSTALL_DIR}/sdk/include/SimTK/simbody)
SET(OPENSIM_HEADERS_DIR ${OPENSIM_INSTALL_DIR}/sdk/include)
INCLUDE_DIRECTORIES(${SIMTK_HEADERS_DIR} ${OPENSIM_HEADERS_DIR})
# Libraries and dlls
SET(OPENSIM_LIBS_DIR ${OPENSIM_INSTALL_DIR}/sdk/lib ${OPENSIM_INSTALL_DIR}/lib)
SET(OPENSIM_DLLS_DIR ${OPENSIM_INSTALL_DIR}/bin)
LINK_DIRECTORIES(${OPENSIM_LIBS_DIR} ${OPENSIM_DLLS_DIR})

# Namespace
#SET(NameSpace "OpenSim_" CACHE STRING "Prefix for simtk lib names, includes trailing '_'. Leave empty to use stock SimTK libraries.")
#MARK_AS_ADVANCED(NameSpace)

ADD_EXECUTABLE(${TARGET} ${SOURCE})

TARGET_LINK_LIBRARIES(${TARGET}
	debug osimSimulation_d	optimized osimSimulation
	debug osimActuators_d	optimized osimActuators
	debug osimCommon_d		optimized osimCommon
	debug osimAnalyses_d	optimized osimAnalyses
	debug osimTools_d		optimized osimTools
	debug ${NameSpace}SimTKcommon_d optimized ${NameSpace}SimTKcommon
	debug ${NameSpace}SimTKmath_d optimized  ${NameSpace}SimTKmath
	debug ${NameSpace}SimTKsimbody_d optimized ${NameSpace}SimTKsimbody
	${PLATFORM_LIBS}
	${MPI_LIBRARIES}
	${OpenCV_LIBS}
)

IF(WIN32)
	SET(PLATFORM_LIBS  pthreadVC2)
ELSE (WIN32)
  SET(NameSpace "")
  IF(APPLE)
	SET(PLATFORM_LIBS  SimTKAtlas)
  ELSE(APPLE)
	SET(PLTAFORM_LIBS SimTKAtlas_Lin_generic)
  ENDIF(APPLE)
ENDIF (WIN32)

# This block copies the additional files into the running directory
# For example vtp, obj files. Add to the end for more extentions
FILE(GLOB DATA_FILES *.vtp *.obj *.osim)
FOREACH (dataFile ${DATA_FILES})
 ADD_CUSTOM_COMMAND(
    TARGET ${TARGET}
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy
    ${dataFile}
    ${OptimizationExample_BINARY_DIR})
 ENDFOREACH (dataFile) 

MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)
MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH)
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
