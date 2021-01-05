cmake_minimum_required(VERSION 3.11.0)

# Set the project name
project(multem CXX CUDA)

# Set the cmake module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Need CUDA and FFTW
find_package(CUDA REQUIRED)
find_package(FFTW REQUIRED)

# Turn off LTO (has problems with cuda)
set(HAS_FLTO False)
set(HAS_FLTO_THIN False)

# Set the coverage compile and link flags for gcc or clang
set(COVERAGE_COMPILE_FLAGS -coverage)
set(COVERAGE_LINK_FLAGS    -coverage)

# Add the automatically determined parts of the RPATH which point to directories
# outside the build tree to the install RPATH. Required for submission to
# clusters which may not allow export of LD_LIBRARY_PATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)

# Add a C/C++ extension
add_library(multem SHARED src/multem.cpp)

# Ensure we are using C++11
target_compile_features(multem PUBLIC cxx_std_11)

# Not specifying CUDA architecture throws and error.
# Setting this option does not pass arch flag to compiler
set_property(TARGET multem PROPERTY CUDA_ARCHITECTURES OFF)

# Set the include directory
target_include_directories(multem PUBLIC src)

# Link to the CUDA libraries
target_link_libraries(multem PUBLIC
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES}
  ${FFTW_LIBRARIES})

# Set the coverage compile and link flags in the debug build
#target_compile_options(multem 
#  PUBLIC "$<$<CONFIG:DEBUG>:${COVERAGE_COMPILE_FLAGS}>")
#target_link_options(multem 
#  PUBLIC "$<$<CONFIG:DEBUG>:${COVERAGE_LINK_FLAGS}>")

# Install the python extension
install(TARGETS multem LIBRARY DESTINATION lib)
install(FILES include/multem.h DESTINATION include)
