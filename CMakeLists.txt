cmake_minimum_required(VERSION 3.14)
project(SQUARE_LATTICE_TJ_MODEL)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
set(CMAKE_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src/")

option(RealCode "tensor element type double or complex" OFF)
if (RealCode)
    add_definitions(-DUSE_REAL=1)
else ()
    add_definitions(-DUSE_COMPLEX=1)
endif ()

option(HIGHER_SYMMETRY "particle number conservation." ON)
if (HIGHER_SYMMETRY)
    add_definitions(-DSYM_LEVEL=1)
else ()
    add_definitions(-DSYM_LEVEL=0)
endif ()


# Build type.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES Debug)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif ()

option(QLTEN_TIMING_MODE "Timing mode." OFF)
if (QLTEN_TIMING_MODE)
    add_definitions(-DQLTEN_TIMING_MODE)
endif ()

option(QLMPS_TIMING_MODE "Timing mode." OFF)
if (QLMPS_TIMING_MODE)
    add_definitions(-DQLMPS_TIMING_MODE)
endif ()

option(QLTEN_MPI_TIMING_MODE "Timing mode." OFF)
if (QLTEN_MPI_TIMING_MODE)
    add_definitions(-DQLTEN_MPI_TIMING_MODE)
endif ()

option(QLMPS_MPI_TIMING_MODE "Timing mode." OFF)
if (QLMPS_MPI_TIMING_MODE)
    add_definitions(-DQLMPS_MPI_TIMING_MODE)
endif ()

option(COMPILE_FIX_CODE "compile the fixing code" OFF)

# Compilation and linking control.
set(CMAKE_CXX_STANDARD 17)
add_definitions(-Wall -g)

option(QLMPS_USE_GPU "Use GPU" OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/") # HPTT/cutensor

if (QLMPS_USE_GPU)
    add_definitions(-DUSE_GPU=1)

    # Enable CUDA language support
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)  # Modern CUDA Toolkit detection

    include_directories(${CUDAToolkit_INCLUDE_DIRS})
    message(STATUS "CUDA_INCLUDE_DIRS: '${CUDAToolkit_INCLUDE_DIRS}'")

    # Set CUDA architectures
    set(CMAKE_CUDA_ARCHITECTURES 70 80) # For V100 & A100, adjust if needed
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -lineinfo -Xcompiler=-fPIC")

    # -----------------------------------------------------------
    # Attempt to find cuTENSOR in the same base directory as cuBLAS
    # -----------------------------------------------------------
    # Extract the CUDA Toolkit root path (e.g., /usr/local/cuda or HPC SDK's CUDA dir)
    get_filename_component(CUDA_TOOLKIT_ROOT "${CUDAToolkit_INCLUDE_DIRS}" DIRECTORY)
    get_filename_component(CUDA_TOOLKIT_ROOT "${CUDA_TOOLKIT_ROOT}" DIRECTORY)  # Strips off "include"

    # Check if CUDA is part of the NVIDIA HPC SDK (path contains "hpc_sdk")
    if (CUDA_TOOLKIT_ROOT MATCHES "hpc_sdk")
        # Deduce HPC SDK root by removing "/cuda/..." from the path
        string(REGEX REPLACE "/cuda/.*" "" HPC_SDK_ROOT "${CUDA_TOOLKIT_ROOT}")
        set(CUTENSOR_HPC_INC_DIR "${HPC_SDK_ROOT}/math_libs/include")
        set(CUTENSOR_HPC_LIB_DIR "${HPC_SDK_ROOT}/math_libs/lib64")
    else ()
        # For non-HPC SDK installations, check sibling directories
        set(CUTENSOR_HPC_INC_DIR "${CUDA_TOOLKIT_ROOT}/include")
        set(CUTENSOR_HPC_LIB_DIR "${CUDA_TOOLKIT_ROOT}/lib64")
    endif ()

    # Search for cuTENSOR header (account for subdirectory)
    find_path(CUTENSOR_INCLUDE_DIR
            NAMES cutensor.h
            HINTS ${CUTENSOR_HPC_INC_DIR}
            PATH_SUFFIXES cutensor  # Header is in a subdirectory
            DOC "cuTENSOR include directory"
    )

    # Search for cuTENSOR library
    find_library(CUTENSOR_LIBRARY
            NAMES cutensor
            HINTS ${CUTENSOR_HPC_LIB_DIR}
            DOC "cuTENSOR library"
    )

    if (CUTENSOR_INCLUDE_DIR AND CUTENSOR_LIBRARY)
        message(STATUS "Found cuTENSOR in CUDA Toolkit/HPC SDK paths:")
        message(STATUS "  Include dir: ${CUTENSOR_INCLUDE_DIR}")
        message(STATUS "  Library: ${CUTENSOR_LIBRARY}")
    else ()
        message(STATUS "cuTENSOR not found in system. Searching in user directory...")
        # Use custom FindCUTENSOR.cmake from project's cmake directory
        list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")
        find_package(CUTENSOR REQUIRED)
    endif ()

    file(GLOB_RECURSE TEST_SOURCES src/*.cpp)
    message(STATUS "TEST_SOURCES" ${TEST_SOURCES})
    # Treat all .cc files in tests/ as CUDA files
    set_source_files_properties(${TEST_SOURCES} PROPERTIES LANGUAGE CUDA)
else ()
    find_package(hptt)
    message(STATUS "HPTT Lib : " "${hptt_LIBRARY}")
endif ()

include(MathBackend)

FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

find_path(TENSOR_HEADER_PATH "qlten")
find_path(MPS_HEADER_PATH "qlmps")


# Intel MKL multi-thread may conflict with other compiler's OpenMP
# Close OpenMP if using MKL with other compiler
if (NOT BLA_VENDOR STREQUAL "Intel10_64lp" OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    find_package(OpenMP REQUIRED)
    link_libraries(OpenMP::OpenMP_CXX)
endif ()

include_directories(
        ${BLAS_INCLUDE_DIRS}
        ${MPI_CXX_HEADER_DIR}
        ${TENSOR_HEADER_PATH}
        ${MPS_HEADER_PATH})

link_libraries(BLAS::BLAS
        LAPACK::LAPACK
        MPI::MPI_CXX)

if (NOT QLMPS_USE_GPU)
    include_directories(${hptt_INCLUDE_DIR})
    link_libraries(${hptt_LIBRARY})
else ()
    include_directories(${CUTENSOR_INCLUDE_DIR})
    link_libraries(
            CUDA::cublas_static
            CUDA::cusolver_static
            ${CUTENSOR_LIBRARY}
    )
endif ()

add_library(common_objects OBJECT src/myutil.cpp)

add_executable(dmrg src/dmrg.cpp $<TARGET_OBJECTS:common_objects>)
add_executable(mpogen src/mpogen.cpp $<TARGET_OBJECTS:common_objects>)
add_executable(vmps1FM src/vmps1FM.cpp $<TARGET_OBJECTS:common_objects>)
add_executable(vmps src/vmps.cpp $<TARGET_OBJECTS:common_objects>)
add_executable(noised_vmps src/noised_vmps.cpp $<TARGET_OBJECTS:common_objects>)
add_executable(measure1 src/measure1.cpp $<TARGET_OBJECTS:common_objects>)
add_executable(measureSC src/measureSC_fix_ref.cpp $<TARGET_OBJECTS:common_objects>)
# Extended version: user-defined reference bond
add_executable(measureSC_fix_ref_ext
    src/measureSC_fix_ref_extension.cpp
    $<TARGET_OBJECTS:common_objects>
)

add_executable(plot_lattice src/plot_lattice.cpp)



add_executable(move_center src/move_center.cpp $<TARGET_OBJECTS:common_objects>)
if (NOT QLMPS_USE_GPU)
    add_executable(measureS src/measure_total_spin.cpp $<TARGET_OBJECTS:common_objects>)
    #add_executable(measure2 src/measure2.cpp src/myutil.cpp src/operators.cpp)
    #add_executable(measure3 src/measure3.cpp src/myutil.cpp src/operators.cpp)
    #
    #add_executable(measureSCxx src/measureSCxx.cpp src/myutil.cpp src/operators.cpp)
    #add_executable(measurePlasc src/measurePlaquetteDSC.cpp src/myutil.cpp src/operators.cpp)
    #add_executable(test1 src/test1_parallel.cpp src/myutil.cpp)
endif ()
if (COMPILE_FIX_CODE)

    add_executable(fix_update_renv src/fix_update_renv.cpp src/myutil.cpp)
    add_executable(fix_update_lenv src/fix_update_lenv.cpp src/myutil.cpp)
    add_executable(fix_mps src/fix_mps.cpp src/myutil.cpp)
    add_executable(fix_env src/fix_env.cpp src/myutil.cpp)
    add_executable(continue src/continue.cpp src/myutil.cpp)
endif ()
