###########################################
## For building with the rest of OpenSim ##
###########################################

cmake_minimum_required(VERSION 3.2)
project(OptimizationExample)

set(OpenSim_SOURCE_DIR "~/opensim-workspace/opensim-core-source")
set(OPENSIM_INSTALL_DIR "~/opensim-core")
# set(SIMTK_HEADERS_DIR ${OPENSIM_INSTALL_DIR}/sdk/Simbody)
# set(OPENSIM_HEADERS_DIR ${OPENSIM_INSTALL_DIR}/sdk/include)
# set(OPENSIM_LIBS_DIR ${OPENSIM_INSTALL_DIR}/sdk/lib ${OPENSIM_INSTALL_DIR}/lib)
# set(OPENSIM_DLLS_DIR ${OPENSIM_INSTALL_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(EXTRA_COMPILE_FLAGS "-w")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Generate compile_commands.json for clangd

find_package(OpenSim REQUIRED PATHS ${OPENSIM_INSTALL_DIR})
find_package(MPI REQUIRED)
find_package(OpenCV REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include(CheckCXXCompilerFlag)

include_directories(${MPI_INCLUDE_PATH})
# include_directories(${OpenCV_INCLUDE_DIRS})
# include_directories(${OpenSim_SOURCE_DIR} ${OpenSim_SOURCE_DIR}/Vendors)
# include_directories(${SIMTK_HEADERS_DIR} ${OPENSIM_HEADERS_DIR})
# link_directories(${OPENSIM_LIBS_DIR} ${OPENSIM_DLLS_DIR})

SET(TARGET optimizationExample CACHE TYPE STRING)
FILE(GLOB SOURCE_FILES *.cpp *.c)
FILE(GLOB INCLUDE_FILES *.h)
SET_SOURCE_FILES_PROPERTIES( ${CMAKE_SOURCE_DIR}/cmaes.c PROPERTIES LANGUAGE CXX )

add_executable(${TARGET} ${SOURCE_FILES} ${INCLUDE_FILES})
target_link_libraries(${TARGET}
	${OpenSim_LIBRARIES}
	${MPI_LIBRARIES}
	${OpenCV_LIBS}
)

file(GLOB DATA_FILES *.vtp *.obj)
foreach(dataFile ${DATA_FILES})
    add_custom_command(
        TARGET ${TARGET}
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy
        ${dataFile}
        ${CMAKE_BINARY_DIR})
endforeach(dataFile)
