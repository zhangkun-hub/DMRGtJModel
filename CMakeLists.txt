cmake_minimum_required(VERSION 3.14)
project(SQUARE_LATTICE_TJ_MODEL)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
set(CMAKE_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src/")

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


if (QLMPS_USE_GPU)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules/") # find cutensor

    add_definitions(-DUSE_GPU=1)

    # Enable CUDA language support
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)  # Modern CUDA Toolkit detection

    include_directories(${CUDAToolkit_INCLUDE_DIRS})
    message(STATUS "CUDA_INCLUDE_DIRS: '${CUDAToolkit_INCLUDE_DIRS}'")

    # Set CUDA architectures
    set(CMAKE_CUDA_ARCHITECTURES 70 80) # For V100 & A100, adjust if needed
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -lineinfo -Xcompiler=-fPIC")


    find_package(CUTENSOR REQUIRED)

    file(GLOB_RECURSE TEST_SOURCES src/*.cpp)
    message(STATUS "TEST_SOURCES" ${TEST_SOURCES})
    # Treat all .cc files in tests/ as CUDA files
    set_source_files_properties(${TEST_SOURCES} PROPERTIES LANGUAGE CUDA)
else ()
    find_path(hptt_INCLUDE_DIR "hptt.h")
    find_library(hptt_LIBRARY "libhptt.a")
    message(STATUS "HPTT Lib : " "${hptt_LIBRARY}")
endif ()


if (APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(MATH_LIB_COMPILE_FLAGS "-I$ENV{MKLROOT}/include")
        #Need test
        set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl)
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        # May not work
        set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -L$ENV{MKLROOT}/lib -L$ENV{CMPLR_ROOT}/mac/compiler/lib/ -liomp5 -lpthread -lm -ldl)
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        # Note as of Intel oneAPI 2021.2, "source /opt/intel/oneapi/compiler/latest/env/vars.sh"
        set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -L$ENV{MKLROOT}/lib -L$ENV{CMPLR_ROOT}/mac/compiler/lib/ -Wl, -rpath $ENV{CMPLR_ROOT}/mac/compiler/lib/libiomp5.dylib -liomp5 -lpthread -lm -ldl)
    endif ()
elseif (APPLE) #ARM
    add_definitions(-DUSE_OPENBLAS)
    set(BLA_VENDOR OpenBLAS)
    #        FIND_PACKAGE(BLAS REQUIRED)
    #        FIND_PACKAGE(LAPACK REQUIRED)
    set(OpenBLAS_ROOT "/opt/homebrew/opt/openblas/")
    set(Lapack_ROOT "/opt/homebrew/opt/lapack")
    message(${OpenBLAS_ROOT})
    set(OpenBLAS_INCLUDE_DIRS "${OpenBLAS_ROOT}/include")
    set(OpenBLAS_LIBRARIES "${OpenBLAS_ROOT}/lib/libblas.dylib")
    message(${OpenBLAS_LIBRARIES})
    set(MATH_LIB_COMPILE_FLAGS -I${OpenBLAS_INCLUDE_DIRS} -pthread)
    set(MATH_LIB_LINK_FLAGS ${OpenBLAS_LIBRARIES} ${OpenBLAS_ROOT}/lib/liblapack.dylib -lm -lpthread -ldl -fopenmp -lclapack)

elseif (UNIX)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "~/.local/")
    message("CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")

    #UNIX include APPLE, but we except it here
    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(MATH_LIB_COMPILE_FLAGS "-I$ENV{MKLROOT}/include")
        set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl)
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        # Link the Intel's OpenMP library to avoid performance issue when the library calls the MKL's gesdd function.
        set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -L$ENV{MKLROOT}/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -lpthread -Wl,--no-as-needed -ldl)
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -L$ENV{MKLROOT}/lib/intel64 -liomp5 -lpthread -lm -ldl)
    endif ()
endif ()


find_path(TENSOR_HEADER_PATH "qlten")
find_path(MPS_HEADER_PATH "qlmps")
find_package(hptt)
find_package(MPI REQUIRED)
#set omp flag
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(OMP_FLAGS -qopenmp)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(OMP_FLAGS -fopenmp)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(OMP_FLAGS -fopenmp)
endif ()


add_compile_options(${MATH_LIB_COMPILE_FLAGS} ${OMP_FLAGS})

if (NOT QLMPS_USE_GPU)
    include_directories(${MPS_HEADER_PATH}
            ${TENSOR_HEADER_PATH}
            ${hptt_INCLUDE_DIR}
            ${MPI_CXX_HEADER_DIR}
    )
    link_libraries(
            ${hptt_LIBRARY}
            ${MATH_LIB_LINK_FLAGS}
            MPI::MPI_CXX
    )
else ()
    include_directories(${MPS_HEADER_PATH}
            ${TENSOR_HEADER_PATH}
            ${CUTENSOR_INCLUDE_DIR}
            ${MPI_CXX_HEADER_DIR}
    )
    link_libraries(
            CUDA::cublas_static
            CUDA::cusolver_static
            ${CUTENSOR_LIBRARY}
            ${MATH_LIB_LINK_FLAGS}
            MPI::MPI_CXX
    )
endif ()

add_library(common_objects OBJECT src/myutil.cpp src/operators.cpp)

add_executable(dmrg src/dmrg.cpp $<TARGET_OBJECTS:common_objects>)
add_executable(mpogen src/mpogen.cpp $<TARGET_OBJECTS:common_objects>)
add_executable(vmps1FM src/vmps1FM.cpp $<TARGET_OBJECTS:common_objects>)
add_executable(vmps src/vmps.cpp $<TARGET_OBJECTS:common_objects>)
add_executable(measure1 src/measure1.cpp $<TARGET_OBJECTS:common_objects>)
add_executable(measureSC src/measureSC_fix_ref.cpp $<TARGET_OBJECTS:common_objects>)

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
    add_executable(move_center src/move_center.cpp src/myutil.cpp)
    add_executable(fix_update_renv src/fix_update_renv.cpp src/myutil.cpp)
    add_executable(fix_update_lenv src/fix_update_lenv.cpp src/myutil.cpp)
    add_executable(fix_mps src/fix_mps.cpp src/myutil.cpp)
    add_executable(fix_env src/fix_env.cpp src/myutil.cpp)
    add_executable(continue src/continue.cpp src/myutil.cpp)
endif ()
