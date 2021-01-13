/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * MULTEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef INPUT_MULTISLICE_H
#define INPUT_MULTISLICE_H

#include <input_multislice_api.h>
#include "types.cuh"
#include "atomic_data_mt.hpp"
#include "slicing.hpp"
#include "memory_info.cuh"

namespace mt {

  template <typename T>
  void Input_Multislice<T>::validate_parameters()
  {
    pn_seed = max(1, pn_seed);

    pn_nconf = (!is_frozen_phonon()) ? 1 : max(1, pn_nconf);

    fp_iconf_0 = (!is_frozen_phonon()) ? 1 : (pn_single_conf) ? pn_nconf : 1;

    islice = max(0, islice);

    if (isZero(Vrl))
    {
      Vrl = c_Vrl;
    }

    if (isZero(nR))
    {
      nR = c_nR;
    }

    dp_Shift = (is_PED()) ? true : false;

    if (!is_scanning())
    {
      scanning.set_default();
    }
    scanning.set_grid();

    lambda = get_lambda(E_0);

    // tem_simulation sign
    mul_sign = (reverse_multislice) ? -1 : 1;

    theta = set_incident_angle(theta);
    nrot = max(1, nrot);
    if (!is_PED_HCTEM())
    {
      nrot = 1;
    }

    // set beam position
    set_beam_position(iw_x, iw_y);

    if (is_EELS_EFTEM())
    {
      pn_coh_contrib = false;
      interaction_model = mt::eESIM_Multislice;
      illumination_model = eIM_Coherent;
    }

    if (is_EWFS_EWRS())
    {
      pn_coh_contrib = true;
      illumination_model = eIM_Coherent;
    }

    // It will be changed when you include spatial incoherences (beta)
    if (is_CBED_CBEI() && !is_illu_mod_full_integration())
    {
      illumination_model = eIM_Coherent;
    }

    slice_storage = is_slice_storage();

    if (!is_multislice())
    {
      islice = 0;
      potential_slicing = ePS_Planes;
      if (is_through_slices())
      {
        thick_type = eTT_Through_Thick;
      }
      slice_storage = slice_storage || !is_whole_spec();
    }

    if (is_subslicing() && is_through_thick())
    {
      thick_type = eTT_Whole_Spec;
    }

    if (!is_subslicing())
    {
      pn_dim.z = false;
    }

    if (is_spec_rot_active())
    {
      thick_type = eTT_Whole_Spec;
      spec_rot_u0.normalized();
      // get geometric center
      if (spec_rot_center_type == eRPT_geometric_center)
      {
        spec_rot_center_p = r3d<T>(atoms.x_mean, atoms.y_mean, atoms.z_mean);
      }
      // rotate atoms
      rotate_atoms(atoms, spec_rot_theta, spec_rot_u0, spec_rot_center_p);
      // get statistic
      atoms.get_statistic();
      // reset theta
      spec_rot_theta = 0;
    }
      
    // temporal fix for phase object
    if (!is_multislice() && is_whole_spec())
    {
      potential_slicing = ePS_dz_Proj;
    }

    // match slicing with the require thickness
    Slicing<T> slicing;
    slicing.match_thickness(potential_slicing, atoms, thick_type, thick);

    // remove atoms outside the thickness range
    T ee_z = 0.1;
    T z_e = thick.back() + ee_z;

    if (atoms.z_max > z_e)
    {
      T z_0 = atoms.z_min - ee_z;
      remove_atoms_outside_z_range(atoms, z_0, z_e);
      // get statistic
      atoms.get_statistic();
    }

    /************* verify lenses parameters **************/

    if (isZero(cond_lens.si_sigma))
    {
      temporal_spatial_incoh = eTSI_Temporal;
    }

    if (!is_illu_mod_full_integration())
    {
      cond_lens.si_rad_npts = 1;
      cond_lens.si_azm_npts = 1;
      cond_lens.ti_npts = 1;

      obj_lens.si_rad_npts = 1;
      obj_lens.si_azm_npts = 1;
      obj_lens.ti_npts = 1;
    }

    if (is_incoh_temporal())
    {
      cond_lens.si_rad_npts = 1;
      cond_lens.si_azm_npts = 1;

      obj_lens.si_rad_npts = 1;
      obj_lens.si_azm_npts = 1;
    }
    else if (is_incoh_spatial())
    {
      cond_lens.ti_npts = 1;
      obj_lens.ti_npts = 1;
    }

    // set condenser lens parameters
    cond_lens.set_input_data(E_0, grid_2d);
    if (atoms.empty())
    {
      cond_lens.zero_defocus_plane = 0;
    }
    else
    {
      cond_lens.zero_defocus_plane = cond_lens.get_zero_defocus_plane(atoms.z_min, atoms.z_max);
    }

    cond_lens.zero_defocus_type = eZDT_User_Define;

    // set objetive lens parameters
    obj_lens.set_input_data(E_0, grid_2d);

    // validate output area
    validate_output_area();
  }

}

#endif
