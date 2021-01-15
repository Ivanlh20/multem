/*
 *  misc.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of
 *  which is included in the root directory of this package.
 */
#ifndef MULTEM_PYTHON_MISC
#define MULTEM_PYTHON_MISC

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {

template <typename T>
T mrad_2_rAngs(T E0, T theta) {
  return rad_2_rAngs(E0, (T)(theta * mt::c_mrad_2_rad));
}

template <typename T>
T mrad_2_sigma(T E0, T theta) {
  return rad_2_sigma(E0, (T)(theta * mt::c_mrad_2_rad));
}

template <typename T>
T scherzer_aperture(T E0, T c_30) {
  return mt::get_Scherzer_aperture(E0, c_30 * mt::c_mm_2_Angs) / mt::c_mrad_2_rad;
}

template <typename T>
T scherzer_defocus(T E0, T c_30) {
  return mt::get_Scherzer_defocus(E0, c_30 * mt::c_mm_2_Angs);
}

}  // namespace mt

void export_misc(py::module_ m) {
  m.def("fwhm_2_sigma", &mt::fwhm_2_sigma<float>);
  m.def("fwhm_2_sigma", &mt::fwhm_2_sigma<double>);

  m.def("hwhm_2_sigma", &mt::hwhm_2_sigma<float>);
  m.def("hwhm_2_sigma", &mt::hwhm_2_sigma<double>);

  m.def("iehwgd_2_sigma", &mt::iehwgd_2_sigma<float>);
  m.def("iehwgd_2_sigma", &mt::iehwgd_2_sigma<double>);

  m.def("get_gamma", &mt::get_gamma<float>);
  m.def("get_gamma", &mt::get_gamma<double>);

  m.def("get_lambda", &mt::get_lambda<float>);
  m.def("get_lambda", &mt::get_lambda<double>);

  m.def("get_sigma", &mt::get_sigma<float>);
  m.def("get_sigma", &mt::get_sigma<double>);

  m.def("mrad_2_rAng", &mt::mrad_2_rAngs<float>);
  m.def("mrad_2_rAng", &mt::mrad_2_rAngs<double>);

  m.def("mrad_2_sigma", &mt::mrad_2_sigma<float>);
  m.def("mrad_2_sigma", &mt::mrad_2_sigma<double>);

  m.def("rad_2_rAng", &mt::rad_2_rAngs<float>);
  m.def("rad_2_rAng", &mt::rad_2_rAngs<double>);

  m.def("rad_2_sigma", &mt::rad_2_sigma<float>);
  m.def("rad_2_sigma", &mt::rad_2_sigma<double>);

  m.def("scherzer_aperture", mt::scherzer_aperture<double>);
  m.def("scherzer_defocus", mt::scherzer_defocus<double>);
}

#endif
