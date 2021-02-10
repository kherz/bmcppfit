# Simple CMake file for Sim_pulseqSBB mex-file creation
# Kai Herz, 2020

cmake_minimum_required(VERSION 3.18.0)

project(bmcppfit)

# need c++11 standards
set (CMAKE_CXX_STANDARD 11)

# O2 is faster than O3, but CMAKE makes O3 standard for GCC
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   string(REPLACE "-O3 " "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
   set(CMAKE_CXX_FLAGS_RELEASE "-O2 ${CMAKE_CXX_FLAGS_RELEASE}")
endif()


include(FetchContent)

#fetch gsl
FetchContent_Declare(
		gsl
		GIT_REPOSITORY https://github.com/ampl/gsl
		GIT_SHALLOW	ON
		GIT_TAG v2.5.0
	)
FetchContent_GetProperties(gsl)
if(NOT gsl_POPULATED)
	message(STATUS "Populating gsl...")
	FetchContent_Populate(gsl)
	set(GSL_DISABLE_TESTS TRUE)
        set(NO_AMPL_BINDINGS TRUE)
	add_subdirectory(${gsl_SOURCE_DIR} ${gsl_BINARY_DIR})
	message(STATUS "Done.")
endif()


# fetch yaml cpp
FetchContent_Declare(
		yaml-cpp
		GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
		GIT_SHALLOW	ON
		GIT_TAG yaml-cpp-0.6.3
	)
FetchContent_GetProperties(yaml-cpp)
if(NOT yaml-cpp_POPULATED)
	message(STATUS "Populating yaml-cpp...")
	FetchContent_Populate(yaml-cpp)
	#do not build yaml-cpp tests
	set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "overwrite" FORCE)
	add_subdirectory(${yaml-cpp_SOURCE_DIR} ${yaml-cpp_BINARY_DIR})
	message(STATUS "Done.")
endif()

# fetch eigen
FetchContent_Declare(
		eigen
		GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
		GIT_SHALLOW	ON
		GIT_TAG 3.3.8
	)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
	message(STATUS "Populating eigen...")
	FetchContent_Populate(eigen)
	message(STATUS "Done.")
endif()

#fetch argh
FetchContent_Declare(
    argh
    GIT_REPOSITORY https://github.com/adishavit/argh
	GIT_SHALLOW	ON
	GIT_TAG v1.3.1
)  
if(NOT argh_POPULATED)
	message(STATUS "Populating argh...")
	FetchContent_Populate(argh)
	message(STATUS "Done.")
endif()

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


set(SOURCE_FILES BMS.h
                 BlochMcConnellSolver.h
                 BlochMatrix.h
                 SimulationParameters.h
                 BMSim_T.h
				 BMSimFit.h
				 YamlIO.h
				 ${pulseq_SOURCE_DIR}/src/ExternalSequence.cpp
                 main.cpp)
			 		
add_executable(bmcppfit ${SOURCE_FILES})

include_directories(BEFORE SYSTEM ${gsl_SOURCE_DIR} ${gsl_BINARY_DIR}/include)
include_directories(BEFORE SYSTEM ${yaml-cpp_SOURCE_DIR} ${yaml-cpp_BINARY_DIR}/include)
include_directories(${eigen_SOURCE_DIR})
include_directories(${argh_SOURCE_DIR})
include_directories(${pulseq_SOURCE_DIR}/src)


target_link_libraries(bmcppfit gsl)
target_link_libraries(bmcppfit yaml-cpp)



