##############################################################################
# Build file for multem project
##############################################################################

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

# Add the automatically determined parts of the RPATH which point to directories
# outside the build tree to the install RPATH. Required for submission to
# clusters which may not allow export of LD_LIBRARY_PATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)

##############################################################################
# Build the multem library
##############################################################################

# Add a C/C++ library
add_library(multem SHARED src/multem.cu)

# Ensure we are using C++11
target_compile_features(multem PUBLIC cxx_std_11)

# Not specifying CUDA architecture throws and error.
# Setting this option does not pass arch flag to compiler
set_property(TARGET multem PROPERTY CUDA_ARCHITECTURES OFF)

# Set the include directory
target_include_directories(multem PRIVATE src)
target_include_directories(multem PUBLIC include)

# Link to the CUDA libraries
target_link_libraries(multem PUBLIC
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES}
  ${FFTW_LIBRARIES})

# Install the python extension
install(TARGETS multem LIBRARY DESTINATION lib)
install(FILES include/multem.h DESTINATION include)
install(FILES include/lin_alg_def.cuh DESTINATION include)
install(FILES include/math.cuh DESTINATION include)
install(FILES include/safe_types.cuh DESTINATION include)
install(FILES include/multem/constants.h DESTINATION include/multem)

##############################################################################
# Build the test executable
##############################################################################

# Enable testing so we can just call "ctest" if we want
enable_testing()

# Add an executable
add_executable(multem-test test/tests.cpp)

# Ensure we are using C++11
target_compile_features(multem-test PUBLIC cxx_std_11)

# Set the include directory
target_include_directories(multem-test PUBLIC include)

# Link to the multem shared library
target_link_libraries(multem-test multem)

# Add the test to ctest runner
add_test(NAME test COMMAND multem-test)
