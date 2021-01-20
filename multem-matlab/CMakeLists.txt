##############################################################################
# Build file for multem-matlab project
##############################################################################

cmake_minimum_required(VERSION 3.12.0)

# Set the project name
project(multem-matlab CXX)

# Set the cmake module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Need pybind11 for Python C/C++ extensions
find_package(Matlab)
find_package(MULTEM REQUIRED)

# Add the automatically determined parts of the RPATH which point to directories
# outside the build tree to the install RPATH. Required for submission to
# clusters which may not allow export of LD_LIBRARY_PATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)

# Add to the include path
include_directories(mex_files_multem ${MULTEM_INCLUDES})
link_libraries(${MULTEM_LIB})

##############################################################################
# Build the multem-matlab mex files
##############################################################################

matlab_add_mex(NAME ilc_add_amorp_lay SRC mex_files_multem/ilc_add_amorp_lay.cpp)
matlab_add_mex(NAME ilc_amorp_spec SRC mex_files_multem/ilc_amorp_spec.cpp)
matlab_add_mex(NAME ilc_crystal_by_lays SRC mex_files_multem/ilc_crystal_by_lays.cpp)
matlab_add_mex(NAME ilc_feg SRC mex_files_multem/ilc_feg.cpp)
matlab_add_mex(NAME ilc_fwhm_2_sigma SRC mex_files_multem/ilc_fwhm_2_sigma.cpp)
matlab_add_mex(NAME ilc_fxeg_data SRC mex_files_multem/ilc_fxeg_data.cpp)
matlab_add_mex(NAME ilc_fxg SRC mex_files_multem/ilc_fxg.cpp)
matlab_add_mex(NAME ilc_gamma SRC mex_files_multem/ilc_gamma.cpp)
matlab_add_mex(NAME ilc_gmax SRC mex_files_multem/ilc_gmax.cpp)
matlab_add_mex(NAME ilc_hwhm_2_sigma SRC mex_files_multem/ilc_hwhm_2_sigma.cpp)
matlab_add_mex(NAME ilc_iehwgd_2_sigma SRC mex_files_multem/ilc_iehwgd_2_sigma.cpp)
matlab_add_mex(NAME ilc_incident_wave SRC mex_files_multem/ilc_incident_wave.cpp)
matlab_add_mex(NAME ilc_lambda SRC mex_files_multem/ilc_lambda.cpp)
matlab_add_mex(NAME ilc_microscope_aberrations SRC mex_files_multem/ilc_microscope_aberrations.cpp)
matlab_add_mex(NAME ilc_min_spl SRC mex_files_multem/ilc_min_spl.cpp)
matlab_add_mex(NAME ilc_mrad_2_rAng SRC mex_files_multem/ilc_mrad_2_rAng.cpp)
matlab_add_mex(NAME ilc_mrad_2_sigma SRC mex_files_multem/ilc_mrad_2_sigma.cpp)
matlab_add_mex(NAME ilc_multem SRC mex_files_multem/ilc_multem.cpp)
matlab_add_mex(NAME ilc_pr SRC mex_files_multem/ilc_pr.cpp)
matlab_add_mex(NAME ilc_projected_potential SRC mex_files_multem/ilc_projected_potential.cpp)
matlab_add_mex(NAME ilc_propagate SRC mex_files_multem/ilc_propagate.cpp)
matlab_add_mex(NAME ilc_rdf_3d SRC mex_files_multem/ilc_rdf_3d.cpp)
matlab_add_mex(NAME ilc_scherzer_aperture SRC mex_files_multem/ilc_scherzer_aperture.cpp)
matlab_add_mex(NAME ilc_scherzer_defocus SRC mex_files_multem/ilc_scherzer_defocus.cpp)
matlab_add_mex(NAME ilc_sigma SRC mex_files_multem/ilc_sigma.cpp)
matlab_add_mex(NAME ilc_spec_planes SRC mex_files_multem/ilc_spec_planes.cpp)
matlab_add_mex(NAME ilc_spec_rot SRC mex_files_multem/ilc_spec_rot.cpp)
matlab_add_mex(NAME ilc_spec_slicing SRC mex_files_multem/ilc_spec_slicing.cpp)
matlab_add_mex(NAME ilc_transmission_function SRC mex_files_multem/ilc_transmission_function.cpp)
matlab_add_mex(NAME ilc_vp SRC mex_files_multem/ilc_vp.cpp)
matlab_add_mex(NAME ilc_vr SRC mex_files_multem/ilc_vr.cpp)
matlab_add_mex(NAME ilc_vz SRC mex_files_multem/ilc_vz.cpp)
matlab_add_mex(NAME ilc_wave_function SRC mex_files_multem/ilc_wave_function.cpp)
