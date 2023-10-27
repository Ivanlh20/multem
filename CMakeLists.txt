cmake_minimum_required(VERSION 3.27)
project(MULTEM LANGUAGES CXX CUDA)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CUDA_STANDARD 20)

# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)

# Find FFTW
# set(FFTW_DIR "/usr/lib/cmake/fftw3")
set(FFTW_DIR "${CMAKE_SOURCE_DIR}/include")

find_package(FFTW REQUIRED)
find_package(CUDAToolkit REQUIRED)

# Include directories
include_directories(
    ${CUDAToolkit_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/morph/tomlplusplus/include
)

# set(CMAKE_CUDA_COMPILER /usr/local/cuda-11.2/bin/nvcc)

# CUDA compiler flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -allow-unsupported-compiler")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --compiler-bindir=/usr/bin") 
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --default-stream per-thread --expt-extended-lambda --expt-relaxed-constexpr")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --compiler-options=-fexceptions,-fpermissive,-fPIC,-fno-omit-frame-pointer,-pthread,-fwrapv,-march=native,-D_GNU_SOURCE")
set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -g -G -lineinfo -Xcompiler -rdynamic")
set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -O3 -DNDEBUG -Xcompiler -rdynamic")
set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "${CMAKE_CUDA_FLAGS_RELWITHDEBINFO} -g -O2 -DNDEBUG -Xcompiler -rdynamic")
set(CMAKE_CUDA_FLAGS_MINSIZEREL "${CMAKE_CUDA_FLAGS_MINSIZEREL} -O1 -DNDEBUG -Xcompiler -rdynamic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")

# Source file(s)
add_executable(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/morph/morp_op.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# Linked libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
    ${HDF5_CXX_LIBRARIES}
    fftw3f
    fftw3
    fftw3f_threads
    fftw3_threads
)
