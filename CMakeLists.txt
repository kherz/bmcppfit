# Simple CMake file for Sim_pulseqSBB mex-file creation
# Kai Herz, 2020

cmake_minimum_required(VERSION 3.18.0)

project(bmcppfit)


# O2 is faster than O3, but CMAKE makes O3 standard for GCC
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   string(REPLACE "-O3 " "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
   set(CMAKE_CXX_FLAGS_RELEASE "-O2 ${CMAKE_CXX_FLAGS_RELEASE}")
endif()


include(FetchContent)

#fetch pulseq
FetchContent_Declare(
    pulseq
    GIT_REPOSITORY https://github.com/pulseq/pulseq.git
	GIT_SHALLOW	ON
	GIT_TAG v1.3.1
)  
if(NOT pulseq_POPULATED)
	message(STATUS "Populating pulseq...")
	FetchContent_Populate(pulseq)
	message(STATUS "Done.")
endif()


set(SOURCE_FILES BlochMcConnellSolver.h
                 NLSFit.h
                 SimulationParameters.h
				 SimulationParameters.cpp
                 BMSim_T.h
				 BMSimFit.h
				 YamlIO.h
				 ${pulseq_SOURCE_DIR}/src/ExternalSequence.cpp
                 main.cpp)
			 		
add_executable(bmcppfit ${SOURCE_FILES})

find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(bmcppfit PRIVATE Eigen3::Eigen)
find_package(Ceres CONFIG REQUIRED)
target_link_libraries(bmcppfit PRIVATE ceres Ceres::ceres)
find_package(yaml-cpp CONFIG REQUIRED)
target_link_libraries(bmcppfit PRIVATE yaml-cpp)
find_package(argh CONFIG REQUIRED)
target_link_libraries(bmcppfit PRIVATE argh)
	
include_directories(${pulseq_SOURCE_DIR}/src)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    add_compile_options(-DEIGEN_DONT_PARALLELIZE)
    target_link_libraries(bmcppfit PRIVATE OpenMP::OpenMP_CXX)
endif()



