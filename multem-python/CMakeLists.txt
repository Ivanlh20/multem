##############################################################################
# Build file for multem-python project
##############################################################################

cmake_minimum_required(VERSION 3.12.0)

# Set the project name
project(multem-python CXX)

# Set the cmake module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Need pybind11 for Python C/C++ extensions
find_package(MULTEM REQUIRED)

# Add pybind sub directory
add_subdirectory(pybind11)

# Add the automatically determined parts of the RPATH which point to directories
# outside the build tree to the install RPATH. Required for submission to
# clusters which may not allow export of LD_LIBRARY_PATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)

##############################################################################
# Build the multem-python library
##############################################################################

# Add a C/C++ extension
pybind11_add_module(multem_ext src/multem/multem_ext.cpp)

# Ensure we are using C++11
target_compile_features(multem_ext PUBLIC cxx_std_14)

# Enable LTO
set_property(TARGET multem_ext PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

# Set the include directory
target_include_directories(multem_ext PUBLIC 
  src
  ${MULTEM_INCLUDES})

# Link to the CUDA libraries
target_link_libraries(multem_ext PUBLIC ${MULTEM_LIB})

# Install the python extension
install(TARGETS multem_ext LIBRARY DESTINATION src/multem)
