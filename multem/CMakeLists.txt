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
# set(HAS_FLTO False)
# set(HAS_FLTO_THIN False)

# Add the automatically determined parts of the RPATH which point to directories
# outside the build tree to the install RPATH. Required for submission to
# clusters which may not allow export of LD_LIBRARY_PATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)

# Set the cmake build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Set some compile flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

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

# Define all symbols as hidden by default
set_target_properties(multem PROPERTIES CXX_VISIBILITY_PRESET "hidden")
set_target_properties(multem PROPERTIES CUDA_VISIBILITY_PRESET "hidden")

# Turn on LTO
set_property(TARGET multem PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

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
install(DIRECTORY include/multem DESTINATION include)

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
