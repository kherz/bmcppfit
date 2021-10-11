# Simple CMake file for bmcppfit executable 
# Kai Herz, 2021

cmake_minimum_required(VERSION 3.16.0)

project(bmcppfit)

# O2 is faster than O3, but CMAKE makes O3 standard for GCC
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   string(REPLACE "-O3 " "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
   set(CMAKE_CXX_FLAGS_RELEASE "-O2 ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

# we need at least c++ 14 standard for experimental filesystem
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-std=c++17 HAVE_FLAG_STD_CXX17)
if(HAVE_FLAG_STD_CXX17)
    set (CMAKE_CXX_STANDARD 17)
else()
    set (CMAKE_CXX_STANDARD 14)
endif()


include(FetchContent)

#fetch pulseq
FetchContent_Declare(
    pulseq
    GIT_REPOSITORY https://github.com/pulseq/pulseq.git
    GIT_SHALLOW ON
    GIT_TAG v1.3.1
)  
if(NOT pulseq_POPULATED)
	message(STATUS "Populating pulseq...")
	FetchContent_Populate(pulseq)
	message(STATUS "Done.")
endif()

#fetch pulseq-cest
FetchContent_Declare(
    pulseq_cest
    GIT_REPOSITORY https://github.com/kherz/pulseq-cest/
    GIT_TAG d0600a7bd3b8bf874fc000e1fe1b3f2a6298b635
)  
if(NOT pulseq_cest_POPULATED)
	message(STATUS "Populating pulseq-cest...")
	FetchContent_Populate(pulseq_cest)
	message(STATUS "Done.")
endif()


set(SOURCE_FILES 
    ${pulseq_cest_SOURCE_DIR}/pulseq-cest-sim/src/SimulationParameters.cpp
    ${pulseq_cest_SOURCE_DIR}/pulseq-cest-sim/src/BMCSim.cpp
    ${pulseq_SOURCE_DIR}/src/ExternalSequence.h
    ${pulseq_SOURCE_DIR}/src/ExternalSequence.cpp
    src/BMSimFit.h
    src/BMSimFit.cpp
    src/YamlParser.h
    src/YamlParser.cpp
    src/main.cpp)
			 		
add_executable(bmcppfit ${SOURCE_FILES})

include_directories(${pulseq_SOURCE_DIR}/src)
include_directories(${pulseq_cest_SOURCE_DIR}/pulseq-cest-sim/src/)

find_package(Eigen3 CONFIG REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR}/Eigen)
find_package(Ceres CONFIG REQUIRED)
target_link_libraries(bmcppfit PRIVATE ceres Ceres::ceres)
find_package(yaml-cpp CONFIG REQUIRED)
target_link_libraries(bmcppfit PRIVATE yaml-cpp)
find_package(argh CONFIG REQUIRED)
target_link_libraries(bmcppfit PRIVATE argh)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    add_compile_options(-DEIGEN_DONT_PARALLELIZE)
    target_link_libraries(bmcppfit PRIVATE OpenMP::OpenMP_CXX)
endif()

install(TARGETS bmcppfit)

# todo: we assume dynamic linking on windows here, that should be checked
IF (WIN32)
install(CODE [[
file(COPY ${CMAKE_CURRENT_BINARY_DIR}/Release/yaml-cpp.dll DESTINATION ${CMAKE_INSTALL_PREFIX}/bin) 
file(COPY ${CMAKE_CURRENT_BINARY_DIR}/Release/ceres.dll DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
file(COPY ${CMAKE_CURRENT_BINARY_DIR}/Release/gflags.dll DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
file(COPY ${CMAKE_CURRENT_BINARY_DIR}/Release/glog.dll DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
]])
ENDIF()
