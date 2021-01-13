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

#ifndef OUTPUT_MULTISLICE_H
#define OUTPUT_MULTISLICE_H

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "matlab_types.cuh"
#include "atomic_data_mt.hpp"
#include "input_multislice.cuh"
#include "cpu_fcns.hpp"
#include "gpu_fcns.cuh"
#include <output_multislice_api.h>

namespace mt {
  
  template <typename T>
  struct Output_Multislice<T>::Implementation {

    Output_Multislice& parent;
		
    using T_r = T;
		using T_c = complex<T>;
		
    using TVector_hr = host_vector<T>;
		using TVector_hc = host_vector<complex<T>>;
    
		using TVector_dr = device_vector<T>;
		using TVector_dc = device_vector<complex<T>>;
		
    eTEM_Output_Type output_type;

		int ndetector;
		int nx;
		int ny;
		T_r dx;
		T_r dy;
		T_r dr;
    
    host_vector<T_r> x;
    host_vector<T_r> y;
    host_vector<T_r> r;

    host_vector<Det_Int<TVector_hr>> image_tot;
    host_vector<Det_Int<TVector_hr>> image_coh;

    host_vector<TVector_hr> m2psi_tot;
    host_vector<TVector_hr> m2psi_coh;
    host_vector<TVector_hc> psi_coh;
    host_vector<TVector_hr> V;
    host_vector<TVector_hc> trans;
    host_vector<TVector_hc> psi_0;

    std::vector<bool> thk_gpu;
		host_vector<TVector_dr> m2psi_tot_d;
		host_vector<TVector_dr> m2psi_coh_d;
		host_vector<TVector_dc> psi_coh_d;

    Vector<T_c, e_host> psi_zh;
    Vector<T_r, e_host> m2psi_zh;

		int n_thk;
		int n_thk_d;

		Stream<e_host> stream;
		
    Implementation(Output_Multislice<T> &parent_in) 
      : parent(parent_in),
        output_type(eTEMOT_m2psi_tot), 
		    ndetector(0), 
        nx(0), 
        ny(0), 
        dx(0), 
        dy(0), 
        dr(0), 
        n_thk(0), 
        n_thk_d(0) {}
  
    template <class TImplementation>
    void assign(TImplementation &other)
    {
      output_type = other.output_type;
      ndetector = other.ndetector;
      nx = other.nx;
      ny = other.ny;
      dx = other.dx;
      dy = other.dy;
      dr = other.dr;

      x = other.x;
      y = other.y;
      r = other.r;

      image_tot.resize(other.image_tot.size());
      for (auto ithk = 0; ithk < other.image_tot.size(); ithk++)
      {
        image_tot[ithk].image.resize(other.image_tot[ithk].image.size());
        for (auto idet = 0; idet < other.image_tot[ithk].image.size(); idet++)
        {
          image_tot[ithk].image[idet] = other.image_tot[ithk].image[idet];
        }
      }

      image_coh.resize(other.image_coh.size());
      for (auto ithk = 0; ithk < other.image_coh.size(); ithk++)
      {
        image_coh[ithk].image.resize(other.image_coh[ithk].image.size());
        for (auto idet = 0; idet < other.image_coh[ithk].image.size(); idet++)
        {
          image_coh[ithk].image[idet] = other.image_coh[ithk].image[idet];
        }
      }

      m2psi_tot.resize(other.m2psi_tot.size());
      for (auto ithk = 0; ithk < other.m2psi_tot.size(); ithk++)
      {
        m2psi_tot[ithk] = other.m2psi_tot[ithk];
      }

      m2psi_coh.resize(other.m2psi_coh.size());
      for (auto ithk = 0; ithk < other.m2psi_coh.size(); ithk++)
      {
        m2psi_coh[ithk] = other.m2psi_coh[ithk];
      }

      psi_coh.resize(other.psi_coh.size());
      for (auto ithk = 0; ithk < other.psi_coh.size(); ithk++)
      {
        psi_coh[ithk] = other.psi_coh[ithk];
      }

      V.resize(other.V.size());
      for (auto ithk = 0; ithk < other.V.size(); ithk++)
      {
        V[ithk] = other.V[ithk];
      }

      trans.resize(other.trans.size());
      for (auto ithk = 0; ithk < other.trans.size(); ithk++)
      {
        trans[ithk] = other.trans[ithk];
      }

      psi_0.resize(other.psi_0.size());
      for (auto ithk = 0; ithk < other.psi_0.size(); ithk++)
      {
        psi_0[ithk] = other.psi_0[ithk];
      }
    }

    template <class TOutput_Multislice>
    Output_Multislice<T>& operator=(TOutput_Multislice &output_multislice)
    {
      assign(output_multislice);
      return *this;
    }

    void clear()
    {
      output_type = eTEMOT_m2psi_tot;
      ndetector = 0;
      nx = 0;
      ny = 0;
      dx = 0;
      dy = 0;
      dr = 0;

      n_thk = 0;
      n_thk_d = 0;

      x.clear();
      x.shrink_to_fit();

      y.clear();
      y.shrink_to_fit();

      r.clear();
      r.shrink_to_fit();

      image_tot.clear();
      image_tot.shrink_to_fit();

      image_coh.clear();
      image_coh.shrink_to_fit();

      m2psi_tot.clear();
      m2psi_tot.shrink_to_fit();

      m2psi_coh.clear();
      m2psi_coh.shrink_to_fit();

      psi_coh.clear();
      psi_coh.shrink_to_fit();

      /*******************************************/
      psi_zh.clear();
      psi_zh.shrink_to_fit();

      m2psi_zh.clear();
      m2psi_zh.shrink_to_fit();

      thk_gpu.clear();
      thk_gpu.shrink_to_fit();

      m2psi_tot_d.clear();
      m2psi_tot_d.shrink_to_fit();

      m2psi_coh_d.clear();
      m2psi_coh_d.shrink_to_fit();

      psi_coh_d.clear();
      psi_coh_d.shrink_to_fit();

      /*******************************************/
      V.clear();
      V.shrink_to_fit();

      trans.clear();
      trans.shrink_to_fit();

      psi_0.clear();
      psi_0.shrink_to_fit();
    }

    void clean_temporal()
    {
      psi_zh.clear();
      //psi_zh.shrink_to_fit();

      m2psi_zh.clear();
      //m2psi_zh.shrink_to_fit();

      switch (output_type)
      {
        case eTEMOT_image_tot_coh:
        {
          for (auto ithk = 0; ithk < n_thk; ithk++)
          {
            psi_coh[ithk].clear();
            //psi_coh[ithk].shrink_to_fit();

            thk_gpu[ithk] = false;

            if(ithk<n_thk_d)
            {
              psi_coh_d[ithk].clear();
              //psi_coh_d[ithk].shrink_to_fit();
            }
          }
          psi_coh.clear();
          //psi_coh.shrink_to_fit();

          thk_gpu.clear();
          //thk_gpu.shrink_to_fit();

          psi_coh_d.clear();
          //psi_coh_d.shrink_to_fit();
        }
        break;
        case eTEMOT_image_tot:
        {

        }
        break;
        case eTEMOT_m2psi_tot_coh:
        {
          for (auto ithk = 0; ithk < n_thk; ithk++)
          {
            psi_coh[ithk].clear();
            //psi_coh[ithk].shrink_to_fit();

            thk_gpu[ithk] = false;

            if(ithk<n_thk_d)
            {
              m2psi_tot_d[ithk].clear();
              //m2psi_tot_d[ithk].shrink_to_fit();

              m2psi_coh_d[ithk].clear();
              //m2psi_coh_d[ithk].shrink_to_fit();

              psi_coh_d[ithk].clear();
              //psi_coh_d[ithk].shrink_to_fit();
            }
          }
          psi_coh.clear();
          //psi_coh.shrink_to_fit();

          thk_gpu.clear();
          //thk_gpu.shrink_to_fit();

          m2psi_tot_d.clear();
          //m2psi_tot_d.shrink_to_fit();

          m2psi_coh_d.clear();
          //m2psi_coh_d.shrink_to_fit();

          psi_coh_d.clear();
          //psi_coh_d.shrink_to_fit();
        }
        break;
        case eTEMOT_m2psi_tot:
        {
          for (auto ithk = 0; ithk < n_thk; ithk++)
          {
            thk_gpu[ithk] = false;

            if(ithk<n_thk_d)
            {
              m2psi_tot_d[ithk].clear();
              //m2psi_tot_d[ithk].shrink_to_fit();
            }
          }

          thk_gpu.clear();
          //thk_gpu.shrink_to_fit();

          m2psi_tot_d.clear();
          //m2psi_tot_d.shrink_to_fit();
        }
        break;
        case eTEMOT_m2psi_tot_psi_coh:
        {
          for (auto ithk = 0; ithk < n_thk; ithk++)
          {
            thk_gpu[ithk] = false;
            if(ithk<n_thk_d)
            {
              m2psi_tot_d[ithk].clear();
              //m2psi_tot_d[ithk].shrink_to_fit();

              psi_coh_d[ithk].clear();
              //psi_coh_d[ithk].shrink_to_fit();
            }
          }

          thk_gpu.clear();
          //thk_gpu.shrink_to_fit();

          m2psi_tot_d.clear();
          //m2psi_tot_d.shrink_to_fit();

          psi_coh_d.clear();
          //psi_coh_d.shrink_to_fit();
        }
        break;
        case eTEMOT_psi_coh:
        {
          for (auto ithk = 0; ithk < n_thk; ithk++)
          {
            thk_gpu[ithk] = false;
            if(ithk<n_thk_d)
            {
              psi_coh_d[ithk].clear();
              //psi_coh_d[ithk].shrink_to_fit();
            }
          }

          thk_gpu.clear();
          //thk_gpu.shrink_to_fit();

          psi_coh_d.clear();
          psi_coh_d.shrink_to_fit();
        }
        break;
        case eTEMOT_psi_0:
        {

        }
        break;
        case eTEMOT_V:
        {

        }
        break;
        case eTEMOT_trans:
        {

        }
        break;
      }
    }

    template <class TInput_Multislice>
    void set_input_data(TInput_Multislice *input_multislice)
    {
      clear();

      stream.resize(1);

      parent.assign_input_multislice(*input_multislice);

      // set required number of thickness
      n_thk = parent.thick.size();
      n_thk_d = 0;

      thk_gpu.resize(n_thk);
      std::fill(thk_gpu.begin(), thk_gpu.end(), false);

      // check selected device
      auto bb_is_device = parent.system_conf.is_device();

      // get available gpu free memory
      double free_memory_mb = get_free_memory<e_device>() - 10;
      int nxy_r = parent.output_area.nxy();
      int nxy_g = parent.grid_2d.nxy();

      if(bb_is_device)
      {
        psi_zh.resize(nxy_g);
        m2psi_zh.resize(nxy_g);
      }

      ndetector = (parent.is_EELS()) ? 1 : parent.detector.size();

      set_output_grid();

      parent.set_output_type();

      switch (output_type)
      {
        case eTEMOT_image_tot_coh:
        {
          image_tot.resize(n_thk);
          image_coh.resize(n_thk);
          psi_coh.resize(n_thk);

          n_thk_d = (bb_is_device)?cal_n_thk_a<T_c>(free_memory_mb, nxy_g):0;
          n_thk_d = min(n_thk_d, n_thk);

          psi_coh_d.resize(n_thk_d);

          for (auto ithk = 0; ithk < n_thk; ithk++)
          {
            image_tot[ithk].image.resize(ndetector);
            image_coh[ithk].image.resize(ndetector);

            for (auto idet = 0; idet < ndetector; idet++)
            {
              image_tot[ithk].image[idet].resize(nxy());
              image_coh[ithk].image[idet].resize(nxy());
            }

            psi_coh[ithk].resize(nxy_g);

            if(ithk<n_thk_d)
            {
              thk_gpu[ithk] = true;

              psi_coh_d[ithk].resize(nxy_g);
            }
          }
        }
        break;
        case eTEMOT_image_tot:
        {
          image_tot.resize(n_thk);

          for (auto ithk = 0; ithk < n_thk; ithk++)
          {
            image_tot[ithk].image.resize(ndetector);

            for (auto idet = 0; idet < ndetector; idet++)
            {
              image_tot[ithk].image[idet].resize(nxy());
            }
          }
        }
        break;
        case eTEMOT_m2psi_tot_coh:
        {
          m2psi_tot.resize(n_thk);
          m2psi_coh.resize(n_thk);
          psi_coh.resize(n_thk);

          n_thk_d = (bb_is_device)?cal_n_thk_a<T_c>(free_memory_mb, nxy_r+nxy_g):0;
          n_thk_d = min(n_thk_d, n_thk);

          m2psi_tot_d.resize(n_thk_d);
          m2psi_coh_d.resize(n_thk_d);
          psi_coh_d.resize(n_thk_d);

          for (auto ithk = 0; ithk < n_thk; ithk++)
          {
            m2psi_tot[ithk].resize(nxy_r);
            m2psi_coh[ithk].resize(nxy_r);
            psi_coh[ithk].resize(nxy_g);

            if(ithk<n_thk_d)
            {
              thk_gpu[ithk] = true;

              m2psi_tot_d[ithk].resize(nxy_r);
              m2psi_coh_d[ithk].resize(nxy_r);
              psi_coh_d[ithk].resize(nxy_g);
            }
          }
        }
        break;
        case eTEMOT_m2psi_tot:
        {
          m2psi_tot.resize(n_thk);

          n_thk_d = (bb_is_device)?cal_n_thk_a<T_r>(free_memory_mb, nxy_r):0;
          n_thk_d = min(n_thk_d, n_thk);

          m2psi_tot_d.resize(n_thk_d);

          for (auto ithk = 0; ithk < n_thk; ithk++)
          {
            m2psi_tot[ithk].resize(nxy_r);

            if(ithk<n_thk_d)
            {
              thk_gpu[ithk] = true;

              m2psi_tot_d[ithk].resize(nxy_r);
            }
          }
        }
        break;
        case eTEMOT_m2psi_tot_psi_coh:
        {
          m2psi_tot.resize(n_thk);
          psi_coh.resize(n_thk);

          n_thk_d = (bb_is_device)?cal_n_thk_a<T_r>(free_memory_mb, nxy_r+2*nxy_g):0;
          n_thk_d = min(n_thk_d, n_thk);

          m2psi_tot_d.resize(n_thk_d);
          psi_coh_d.resize(n_thk_d);

          for (auto ithk = 0; ithk < n_thk; ithk++)
          {
            m2psi_tot[ithk].resize(nxy_r);
            psi_coh[ithk].resize(nxy_g);

            if(ithk<n_thk_d)
            {
              thk_gpu[ithk] = true;

              m2psi_tot_d[ithk].resize(nxy_r);
              psi_coh_d[ithk].resize(nxy_g);
            }				
          }
        }
        break;
        case eTEMOT_psi_coh:
        {
          psi_coh.resize(n_thk);

          n_thk_d = (bb_is_device)?cal_n_thk_a<T_c>(free_memory_mb, nxy_g):0;
          n_thk_d = min(n_thk_d, n_thk);

          psi_coh_d.resize(n_thk_d);

          for (auto ithk = 0; ithk < n_thk; ithk++)
          {
            psi_coh[ithk].resize(nxy_g);

            if(ithk<n_thk_d)
            {
              thk_gpu[ithk] = true;

              psi_coh_d[ithk].resize(nxy_g);
            }
          }
        }
        break;
        case eTEMOT_psi_0:
        {
          psi_0.resize(n_thk);

          for (auto ithk = 0; ithk < n_thk; ithk++)
          {
            psi_0[ithk].resize(nxy_r);
          }
        }
        break;
        case eTEMOT_V:
        {
          V.resize(n_thk);

          for (auto ithk = 0; ithk < n_thk; ithk++)
          {
            V[ithk].resize(nxy_r);
          }
        }
        break;
        case eTEMOT_trans:
        {
          trans.resize(n_thk);

          for (auto ithk = 0; ithk < n_thk; ithk++)
          {
            trans[ithk].resize(nxy_r);
          }
        }
        break;
      }
    }
      

    void init()
    {
      switch (output_type)
      {
      case eTEMOT_image_tot_coh:
      {
        for (auto ithk = 0; ithk < parent.thick.size(); ithk++)
        {
          for (auto idet = 0; idet < image_tot[ithk].image.size(); idet++)
          {
            mt::fill(stream, image_tot[ithk].image[idet], T_r(0));
            mt::fill(stream, image_coh[ithk].image[idet], T_r(0));
          }
          mt::fill(stream, psi_coh[ithk], T_c(0));

          if(ithk<n_thk_d)
          {
            thrust::fill(psi_coh_d[ithk].begin(), psi_coh_d[ithk].end(), T_c(0));
          }
        }
      }
      break;
      case eTEMOT_image_tot:
      {
        for (auto ithk = 0; ithk < parent.thick.size(); ithk++)
        {
          for (auto idet = 0; idet < image_tot[ithk].image.size(); idet++)
          {
            mt::fill(stream, image_tot[ithk].image[idet], T_r(0));
          }
        }
      }
      break;
      case eTEMOT_m2psi_tot_coh:
      {
        for (auto ithk = 0; ithk < parent.thick.size(); ithk++)
        {
          mt::fill(stream, m2psi_tot[ithk], T_r(0));
          mt::fill(stream, m2psi_coh[ithk], T_r(0));
          mt::fill(stream, psi_coh[ithk], T_c(0));

          if(ithk<n_thk_d)
          {
            thrust::fill(m2psi_tot_d[ithk].begin(), m2psi_tot_d[ithk].end(), T_r(0));
            thrust::fill(m2psi_coh_d[ithk].begin(), m2psi_coh_d[ithk].end(), T_r(0));
            thrust::fill(psi_coh_d[ithk].begin(), psi_coh_d[ithk].end(), T_c(0));
          }
        }

      }
      break;
      case eTEMOT_m2psi_tot:
      {
        for (auto ithk = 0; ithk < parent.thick.size(); ithk++)
        {
          mt::fill(stream, m2psi_tot[ithk], T_r(0));

          if(ithk<n_thk_d)
          {
            thrust::fill(m2psi_tot_d[ithk].begin(), m2psi_tot_d[ithk].end(), T_r(0));
          }
        }
      }
      break;
      case eTEMOT_m2psi_tot_psi_coh:
      {
        for (auto ithk = 0; ithk < parent.thick.size(); ithk++)
        {
          mt::fill(stream, m2psi_tot[ithk], T_r(0));
          mt::fill(stream, psi_coh[ithk], T_c(0));

          if(ithk<n_thk_d)
          {
            thrust::fill(m2psi_tot_d[ithk].begin(), m2psi_tot_d[ithk].end(), T_r(0));
            thrust::fill(psi_coh_d[ithk].begin(), psi_coh_d[ithk].end(), T_c(0));
          }
        }
      }
      break;
      case eTEMOT_psi_coh:
      {
        for (auto ithk = 0; ithk < parent.thick.size(); ithk++)
        {
          mt::fill(stream, psi_coh[ithk], T_c(0));

          if(ithk<n_thk_d)
          {
            thrust::fill(psi_coh_d[ithk].begin(), psi_coh_d[ithk].end(), T_c(0));
          }
        }
      }
      break;
      case eTEMOT_psi_0:
      {
        for (auto ithk = 0; ithk < parent.thick.size(); ithk++)
        {
          mt::fill(stream, psi_0[ithk], T_c(0));
        }
      }
      break;
      case eTEMOT_V:
      {
        for (auto ithk = 0; ithk < parent.thick.size(); ithk++)
        {
          mt::fill(stream, V[ithk], T_r(0));
        }
      }
      break;
      case eTEMOT_trans:
      {
        for (auto ithk = 0; ithk < parent.thick.size(); ithk++)
        {
          mt::fill(stream, trans[ithk], T_c(0));
        }
      }
      break;
      }
    }

    void init_psi_coh()
    {
      for (auto ithk = 0; ithk < parent.thick.size(); ithk++)
      {
        if(ithk<n_thk_d)
        {
          thrust::fill(psi_coh_d[ithk].begin(), psi_coh_d[ithk].end(), T_c(0));
        }
        else
        {
          mt::fill(stream, psi_coh[ithk], T_c(0));
        }
      }
    }

    /***************************************************************************/
    template<class TVector>
    void add_scale_psi_coh(int ithk, T_c w, TVector &phi)
    {
      if(thk_gpu[ithk])
      {
        thrust::transform(thrust::device, phi.begin(), phi.end(), psi_coh_d[ithk].begin(), psi_coh_d[ithk].begin(), functor::add_scale<T_c>(w));
      }
      else
      {
         mt::add_scale_to_host(stream, w, phi, psi_coh[ithk], &psi_zh);
      }
    }

    template<class TVector>
    void add_scale_shift_psi_coh(int ithk, T_c w, TVector &phi)
    {
      if(thk_gpu[ithk])
      {
        mt::add_scale_shift_2d(parent.grid_2d, w, phi, psi_coh_d[ithk], &psi_zh);
      }
      else
      {
         mt::add_scale_shift_2d(parent.grid_2d, w, phi, psi_coh[ithk], &psi_zh);
      }
    }

    template<class TVector>
    void add_scale_crop_shift_psi_coh(int ithk, T_c w, TVector &phi)
    {
      if(thk_gpu[ithk])
      {
         mt::add_scale_crop_shift_2d(parent.grid_2d, w, phi, parent.output_area, psi_coh_d[ithk], &psi_zh);
      }
      else
      {
         mt::add_scale_crop_shift_2d(parent.grid_2d, w, phi, parent.output_area, psi_coh[ithk], &psi_zh);
      }
    }

    template<class TVector>
    void add_scale_crop_shift_m2psi_coh_from_psi(int ithk, T_r w, TVector &phi)
    {
      if(thk_gpu[ithk])
      {
         mt::add_scale_square_crop_shift_2d(parent.grid_2d, w, phi, parent.output_area, m2psi_coh_d[ithk], &psi_zh);
      }
      else
      {
         mt::add_scale_square_crop_shift_2d(parent.grid_2d, w, phi, parent.output_area, m2psi_coh[ithk], &psi_zh);
      }
    }

    template<class TVector>
    void add_scale_crop_shift_m2psi_coh_from_m2psi(int ithk, T_r w, TVector &m2phi)
    {
      if(thk_gpu[ithk])
      {
         mt::add_scale_crop_shift_2d(parent.grid_2d, w, m2phi, parent.output_area, m2psi_coh_d[ithk], &m2psi_zh);
      }
      else
      {
         mt::add_scale_crop_shift_2d(parent.grid_2d, w, m2phi, parent.output_area, m2psi_coh[ithk], &m2psi_zh);
      }
    }

    template<class TVector>
    void add_scale_crop_shift_m2psi_tot_from_psi(int ithk, T_r w, TVector &phi)
    {
      if(thk_gpu[ithk])
      {
         mt::add_scale_square_crop_shift_2d(parent.grid_2d, w, phi, parent.output_area, m2psi_tot_d[ithk], &psi_zh);
      }
      else
      {
         mt::add_scale_square_crop_shift_2d(parent.grid_2d, w, phi, parent.output_area, m2psi_tot[ithk], &psi_zh);
      }
    }

    template<class TVector>
    void add_scale_crop_shift_m2psi_tot_from_m2psi(int ithk, T_r w, TVector &m2phi)
    {
      if(thk_gpu[ithk])
      {
         mt::add_scale_crop_shift_2d(parent.grid_2d, w, m2phi, parent.output_area, m2psi_tot_d[ithk], &m2psi_zh);
      }
      else
      {
         mt::add_scale_crop_shift_2d(parent.grid_2d, w, m2phi, parent.output_area, m2psi_tot[ithk], &m2psi_zh);
      }
    }

    /***************************************************************************/
    template<class TVector>
    void set_psi_coh(int ithk, TVector &phi)
    {
      if(thk_gpu[ithk])
      {
        thrust::copy(phi.begin(), phi.end(), psi_coh_d[ithk].begin());
      }
      else
      {
        thrust::copy(phi.begin(), phi.end(), psi_coh[ithk].begin());
      }
    }

    template<class TVector>
    void set_shift_psi_coh(int ithk, TVector &phi)
    {
      if(thk_gpu[ithk])
      {
         mt::assign_shift_2d(parent.grid_2d, phi, psi_coh_d[ithk], &psi_zh);
      }
      else
      {
         mt::assign_shift_2d(parent.grid_2d, phi, psi_coh[ithk], &psi_zh);
      }
    }

      template<class TVector>
    void set_crop_shift_psi_coh(int ithk, TVector &phi)
    {
      if(thk_gpu[ithk])
      {
         mt::assign_crop_shift_2d(parent.grid_2d, phi, parent.output_area, psi_coh_d[ithk], &psi_zh);
      }
      else
      {
         mt::assign_crop_shift_2d(parent.grid_2d, phi, parent.output_area, psi_coh[ithk], &psi_zh);
      }
    }

      template<class TVector>
    void set_crop_shift_m2psi_coh(int ithk, TVector &m2phi)
    {
      if(thk_gpu[ithk])
      {
         mt::assign_crop_shift_2d(parent.grid_2d, m2phi, parent.output_area, m2psi_coh_d[ithk], &m2psi_zh);
      }
      else
      {
         mt::assign_crop_shift_2d(parent.grid_2d, m2phi, parent.output_area, m2psi_coh[ithk], &m2psi_zh);
      }
    }

      template<class TVector>
    void set_crop_shift_m2psi_tot(int ithk, TVector &m2phi)
    {
      if(thk_gpu[ithk])
      {
         mt::assign_crop_shift_2d(parent.grid_2d, m2phi, parent.output_area, m2psi_tot_d[ithk], &m2psi_zh);
      }
      else
      {
         mt::assign_crop_shift_2d(parent.grid_2d, m2phi, parent.output_area, m2psi_tot[ithk], &m2psi_zh);
      }
    }

    /***************************************************************************/
      template<class TVector>
    void from_psi_coh_2_phi(int ithk, TVector &phi)
    {
      if(thk_gpu[ithk])
      {
         thrust::copy(psi_coh_d[ithk].begin(), psi_coh_d[ithk].end(), phi.begin());
      }
      else
      {
         thrust::copy(psi_coh[ithk].begin(), psi_coh[ithk].end(), phi.begin());
      }
    }

      template<class TVector>
    void from_m2psi_coh_2_m2phi(int ithk, TVector &m2phi)
    {
      if(thk_gpu[ithk])
      {
         thrust::copy(m2psi_coh_d[ithk].begin(), m2psi_coh_d[ithk].end(), m2phi.begin());
      }
      else
      {
         thrust::copy(m2psi_coh[ithk].begin(), m2psi_coh[ithk].end(), m2phi.begin());
      }		
    }

      template<class TVector>
    void from_m2psi_tot_2_m2phi(int ithk, TVector &m2phi)
    {
      if(thk_gpu[ithk])
      {
        thrust::copy(m2psi_tot_d[ithk].begin(), m2psi_tot_d[ithk].end(), m2phi.begin());
      }
      else
      {
        thrust::copy(m2psi_tot[ithk].begin(), m2psi_tot[ithk].end(), m2phi.begin());
      }
    }

      void gather()
    {
      const int n_thk = parent.thick.size();

      switch (output_type)
      {
      case eTEMOT_m2psi_tot_coh:
      {
        for (auto ithk = 0; ithk < n_thk; ithk++)
        {
          if(ithk<n_thk_d)
          {
            thrust::copy(m2psi_tot_d[ithk].begin(), m2psi_tot_d[ithk].end(), m2psi_tot[ithk].begin());
            thrust::copy(m2psi_coh_d[ithk].begin(), m2psi_coh_d[ithk].end(), m2psi_coh[ithk].begin());
          }
        }
      }
      break;
      case eTEMOT_m2psi_tot:
      {
        for (auto ithk = 0; ithk < n_thk; ithk++)
        {
          if(ithk<n_thk_d)
          {
            thrust::copy(m2psi_tot_d[ithk].begin(), m2psi_tot_d[ithk].end(), m2psi_tot[ithk].begin());
          }
        }
      }
      break;
      case eTEMOT_m2psi_tot_psi_coh:
      {
        for (auto ithk = 0; ithk < n_thk; ithk++)
        {
          if(ithk<n_thk_d)
          {
            thrust::copy(m2psi_tot_d[ithk].begin(), m2psi_tot_d[ithk].end(), m2psi_tot[ithk].begin());
            thrust::copy(psi_coh_d[ithk].begin(), psi_coh_d[ithk].end(), psi_coh[ithk].begin());
          }
        }
      }
      break;
      case eTEMOT_psi_coh:
      {
        for (auto ithk = 0; ithk < n_thk; ithk++)
        {
          if(ithk<n_thk_d)
          {
            thrust::copy(psi_coh_d[ithk].begin(), psi_coh_d[ithk].end(), psi_coh[ithk].begin());
          }
        }
      }
      break;
      case eTEMOT_psi_0:
      {
        for (auto ithk = 0; ithk < n_thk; ithk++)
        {
          mt::fft2_shift(stream, parent.grid_2d, psi_0[ithk]);
        }
      }
      break;
      case eTEMOT_V:
      {
        for (auto ithk = 0; ithk < n_thk; ithk++)
        {
          mt::fft2_shift(stream, parent.grid_2d, V[ithk]);
        }
      }
      break;
      case eTEMOT_trans:
      {
        for (auto ithk = 0; ithk < n_thk; ithk++)
        {
          mt::fft2_shift(stream, parent.grid_2d, trans[ithk]);
        }
      }
      break;
      }
    }
		
    void set_output_grid()
		{
			if (parent.is_STEM() || parent.is_EELS())
			{
 				nx = parent.scanning.nx;
				ny = parent.scanning.ny;

				dx = parent.scanning.dRx;
				dy = parent.scanning.dRy;

				x = parent.scanning.x;
				y = parent.scanning.y;
				r = parent.scanning.r;
			}
			else
			{
				nx = parent.output_area.nx();
				ny = parent.output_area.ny();

				const bool is_RS = parent.is_grid_RS();
				dx = (is_RS)?parent.grid_2d.dRx:parent.grid_2d.dgx;
				dy = (is_RS)?parent.grid_2d.dRy:parent.grid_2d.dgy;

				x.resize(nx);
				y.resize(ny);

				for (auto ix = parent.output_area.ix_0; ix < parent.output_area.ix_e; ix++)
				{
					int ix_s = ix-parent.output_area.ix_0;
					x[ix_s] = (is_RS)?parent.grid_2d.Rx(ix):parent.grid_2d.gx(ix);
				}

				for (auto iy = parent.output_area.iy_0; iy < parent.output_area.iy_e; iy++)
				{
					int iy_s = iy-parent.output_area.iy_0;
					y[iy_s] = (is_RS)?parent.grid_2d.Ry(iy):parent.grid_2d.gy(iy);
				}
			}
		}
		
		int nxy() const { return nx*ny; }
      
    std::vector<T> extract_data(ePhonon_Model_Output fp_ctr, eShow_CData show_data, int ithk, int idet)
    {
      std::vector<T> data(nxy());

      switch (output_type)
      {
      case eTEMOT_image_tot_coh:
      {
        if (fp_ctr == eFMO_Total) {
          data.assign(
            image_tot[ithk].image[idet].begin(),
            image_tot[ithk].image[idet].end());
        } else {
          data.assign(
            image_coh[ithk].image[idet].begin(),
            image_coh[ithk].image[idet].end());
        }
      }
      break;
      case eTEMOT_image_tot:
      {
        data.assign(
            image_tot[ithk].image[idet].begin(),
            image_tot[ithk].image[idet].end());
      }
      break;
      case eTEMOT_m2psi_tot_coh:
      {
        if (fp_ctr == eFMO_Total) {
          data.assign(m2psi_tot[ithk].begin(), m2psi_tot[ithk].end());
        } else {
          data.assign(m2psi_coh[ithk].begin(), m2psi_coh[ithk].end());
        }
      }
      break;
      case eTEMOT_m2psi_tot:
      {
        data.assign(m2psi_tot[ithk].begin(), m2psi_tot[ithk].end());
      }
      break;
      case eTEMOT_m2psi_tot_psi_coh:
      {
        if (fp_ctr == eFMO_Total)
        {
          data.assign(m2psi_tot[ithk].begin(), m2psi_tot[ithk].end());
        }
        else
        {
          from_complex_to_real(show_data, psi_coh[ithk], data);
        }
      }
      break;
      case eTEMOT_psi_coh:
      {
        from_complex_to_real(show_data, psi_coh[ithk], data);
      }
      break;
      case eTEMOT_psi_0:
      {
        from_complex_to_real(show_data, psi_0[ithk], data);
      }
      break;
      case eTEMOT_V:
      {
        data.assign(V[ithk].begin(), V[ithk].end());
      }
      break;
      case eTEMOT_trans:
      {
        from_complex_to_real(show_data, trans[ithk], data);
      }
      break;
      }

      return data;
    }
 		
    template <class U>
		int cal_n_thk_a(double memory, int nxy)
		{
			return static_cast<int>(floor(memory/mt::sizeMb<U>(nxy)));
		}


  };
		
  template <typename T>
  Output_Multislice<T>::Output_Multislice()
    : Input_Multislice<T_r>(),
      pimpl(std::make_unique<Output_Multislice<T>::Implementation>(*this)) {}

  template <typename T>
  Output_Multislice<T>::Output_Multislice(const Output_Multislice &other)
    : Input_Multislice<T_r>(), 
      pimpl(std::make_unique<Output_Multislice<T>::Implementation>(*this)) {
    assign(other);
  }
  
  template <typename T>
  Output_Multislice<T>::Output_Multislice(Output_Multislice<T> &&other) = default;
  
  template <typename T>
  Output_Multislice<T>::~Output_Multislice() = default;

  template <typename T>
  Output_Multislice<T>::Implementation& Output_Multislice<T>::data() {
    return *pimpl;
  }
  
  template <typename T>
  const Output_Multislice<T>::Implementation& Output_Multislice<T>::data() const {
    return *pimpl;
  }

  template <typename T>
  eTEM_Output_Type& Output_Multislice<T>::output_type() {
    return pimpl->output_type;
  }

  template <typename T>
  const eTEM_Output_Type& Output_Multislice<T>::output_type() const {
    return pimpl->output_type;
  }

  template <typename T>
  int& Output_Multislice<T>::ndetector() {
    return pimpl->ndetector;
  }

  template <typename T>
  const int& Output_Multislice<T>::ndetector() const {
    return pimpl->ndetector;
  }

  template <typename T>
  int& Output_Multislice<T>::nx() {
    return pimpl->nx;
  }

  template <typename T>
  const int& Output_Multislice<T>::nx() const {
    return pimpl->nx;
  }

  template <typename T>
  int& Output_Multislice<T>::ny() {
    return pimpl->ny;
  }

  template <typename T>
  const int& Output_Multislice<T>::ny() const {
    return pimpl->ny;
  }

  template <typename T>
  Output_Multislice<T>::T_r& Output_Multislice<T>::dx() {
    return pimpl->dx;
  }

  template <typename T>
  const Output_Multislice<T>::T_r& Output_Multislice<T>::dx() const {
    return pimpl->dx;
  }

  template <typename T>
  Output_Multislice<T>::T_r& Output_Multislice<T>::dy() {
    return pimpl->dy;
  }

  template <typename T>
  const Output_Multislice<T>::T_r& Output_Multislice<T>::dy() const {
    return pimpl->dy;
  }

  template <typename T>
  Output_Multislice<T>::T_r& Output_Multislice<T>::dr() {
    return pimpl->dr;
  }

  template <typename T>
  const Output_Multislice<T>::T_r& Output_Multislice<T>::dr() const {
    return pimpl->dr;
  }

  template <typename T>
  template <class TOutput_Multislice>
  void Output_Multislice<T>::assign(TOutput_Multislice &output_multislice)
  {
    assign_input_multislice(output_multislice);

    pimpl->assign(*output_multislice.pimpl);
  }

  template <typename T>
  template <class TOutput_Multislice>
  Output_Multislice<T>& Output_Multislice<T>::operator=(TOutput_Multislice &output_multislice)
  {
    assign(output_multislice);
    return *this;
  }

  template <typename T>
  void Output_Multislice<T>::clear()
  {
    pimpl->clear();
  }

  template <typename T>
  void Output_Multislice<T>::clean_temporal()
  {
    pimpl->clean_temporal();
  }

  template <typename T>
  template <class TInput_Multislice>
  void Output_Multislice<T>::set_input_data(TInput_Multislice *input_multislice)
  {
    pimpl->set_input_data(input_multislice);
  }
		

  template <typename T>
  void Output_Multislice<T>::init()
  {
    pimpl->init();
  }

  template <typename T>
  void Output_Multislice<T>::init_psi_coh()
  {
    pimpl->init_psi_coh();
  }

  /***************************************************************************/
  template <typename T>
  template<class TVector>
  void Output_Multislice<T>::add_scale_psi_coh(int ithk, T_c w, TVector &phi)
  {
    pimpl->add_scale_psi_coh(ithk, w, phi);
  }

  template <typename T>
  template<class TVector>
  void Output_Multislice<T>::add_scale_shift_psi_coh(int ithk, T_c w, TVector &phi)
  {
    pimpl->add_scale_shift_psi_coh(ithk, w, phi);
  }

  template <typename T>
  template<class TVector>
  void Output_Multislice<T>::add_scale_crop_shift_psi_coh(int ithk, T_c w, TVector &phi)
  {
    pimpl->add_scale_crop_shift_psi_coh(ithk, w, phi);
  }

  template <typename T>
  template<class TVector>
  void Output_Multislice<T>::add_scale_crop_shift_m2psi_coh_from_psi(int ithk, T_r w, TVector &phi)
  {
    pimpl->add_scale_crop_shift_m2psi_coh_from_psi(ithk, w, phi);
  }

  template <typename T>
  template<class TVector>
  void Output_Multislice<T>::add_scale_crop_shift_m2psi_coh_from_m2psi(int ithk, T_r w, TVector &m2phi)
  {
    pimpl->add_scale_crop_shift_m2psi_coh_from_m2psi(ithk, w, m2phi);
  }

  template <typename T>
  template<class TVector>
  void Output_Multislice<T>::add_scale_crop_shift_m2psi_tot_from_psi(int ithk, T_r w, TVector &phi)
  {
    pimpl->add_scale_crop_shift_m2psi_tot_from_psi(ithk, w, phi);
  }

  template <typename T>
  template<class TVector>
  void Output_Multislice<T>::add_scale_crop_shift_m2psi_tot_from_m2psi(int ithk, T_r w, TVector &m2phi)
  {
    pimpl->add_scale_crop_shift_m2psi_tot_from_m2psi(ithk, w, m2phi);
  }

  /***************************************************************************/
  template <typename T>
  template<class TVector>
  void Output_Multislice<T>::set_psi_coh(int ithk, TVector &phi)
  {
    pimpl->set_psi_coh(ithk, phi);
  }

  template <typename T>
  template<class TVector>
  void Output_Multislice<T>::set_shift_psi_coh(int ithk, TVector &phi)
  {
    pimpl->set_shift_psi_coh(ithk, phi);
  }

  template <typename T>
    template<class TVector>
  void Output_Multislice<T>::set_crop_shift_psi_coh(int ithk, TVector &phi)
  {
    pimpl->set_crop_shift_psi_coh(ithk, phi);
  }

  template <typename T>
  template<class TVector>
  void Output_Multislice<T>::set_crop_shift_m2psi_coh(int ithk, TVector &m2phi)
  {
    pimpl->set_crop_shift_m2psi_coh(ithk, m2phi);
  }

  template <typename T>
  template<class TVector>
  void Output_Multislice<T>::set_crop_shift_m2psi_tot(int ithk, TVector &m2phi)
  {
    pimpl->set_crop_shift_m2psi_tot(ithk, m2phi);
  }

  /***************************************************************************/
  template <typename T>
  template<class TVector>
  void Output_Multislice<T>::from_psi_coh_2_phi(int ithk, TVector &phi)
  {
    pimpl->from_psi_coh_2_phi(ithk, phi);
  }

  template <typename T>
  template<class TVector>
  void Output_Multislice<T>::from_m2psi_coh_2_m2phi(int ithk, TVector &m2phi)
  {
    pimpl->from_m2psi_coh_2_m2phi(ithk, m2phi);
  }

  template <typename T>
  template<class TVector>
  void Output_Multislice<T>::from_m2psi_tot_2_m2phi(int ithk, TVector &m2phi)
  {
    pimpl->from_m2psi_tot_2_m2phi(ithk, m2phi);
  }

  template <typename T>
  void Output_Multislice<T>::gather()
  {
    pimpl->gather();
  }
    
  template <typename T>
  std::vector<T> Output_Multislice<T>::extract_data(ePhonon_Model_Output fp_ctr, eShow_CData show_data, int ithk, int idet)
  {
    return pimpl->extract_data(fp_ctr, show_data, ithk, idet);
  }
    
  template <typename T>
  void Output_Multislice<T>::set_output_grid() {
    pimpl->set_output_grid();
  }
		
  template <typename T>
  template <class TInput_Multislice>
  void Output_Multislice<T>::assign_input_multislice(TInput_Multislice &input_multislice)
  {
    this->system_conf = input_multislice.system_conf;

    this->interaction_model = input_multislice.interaction_model;
    this->potential_type = input_multislice.potential_type;

    this->operation_mode = input_multislice.operation_mode;
    this->slice_storage = input_multislice.slice_storage;
    this->reverse_multislice = input_multislice.reverse_multislice;
    this->mul_sign = input_multislice.mul_sign;
    this->Vrl = input_multislice.Vrl;
    this->nR = input_multislice.nR;

    this->pn_model = input_multislice.pn_model;
    this->pn_coh_contrib = input_multislice.pn_coh_contrib;
    this->pn_dim = input_multislice.pn_dim;
    this->fp_dist = input_multislice.fp_dist;
    this->pn_seed = input_multislice.pn_seed;
    this->pn_single_conf = input_multislice.pn_single_conf;
    this->pn_nconf = input_multislice.pn_nconf;

    this->atoms = input_multislice.atoms;
    this->is_crystal = input_multislice.is_crystal;

    this->spec_rot_theta = input_multislice.spec_rot_theta;
    this->spec_rot_u0 = input_multislice.spec_rot_u0;
    this->spec_rot_center_type = input_multislice.spec_rot_center_type;
    this->spec_rot_center_p = input_multislice.spec_rot_center_p;

    this->thick_type = input_multislice.thick_type;
    this->thick = input_multislice.thick;

    this->potential_slicing = input_multislice.potential_slicing;

    this->grid_2d = input_multislice.grid_2d;
    this->output_area = input_multislice.output_area;

    this->simulation_type = input_multislice.simulation_type;

    this->iw_type = input_multislice.iw_type;
    this->iw_psi = input_multislice.iw_psi;
    this->iw_x = input_multislice.iw_x;
    this->iw_y = input_multislice.iw_y;

    this->E_0 = input_multislice.E_0;
    this->theta = input_multislice.theta;
    this->phi = input_multislice.phi;
    this->nrot = input_multislice.nrot;

    this->illumination_model = input_multislice.illumination_model;
    this->temporal_spatial_incoh = input_multislice.temporal_spatial_incoh;

    this->cond_lens = input_multislice.cond_lens;
    this->obj_lens = input_multislice.obj_lens;

    this->scanning = input_multislice.scanning;
    this->detector = input_multislice.detector;

    this->eels_fr = input_multislice.eels_fr;

    this->cdl_var_type = input_multislice.cdl_var_type;
    this->cdl_var = input_multislice.cdl_var;

    this->iscan = input_multislice.iscan;
    this->beam_x = input_multislice.beam_x;
    this->beam_y = input_multislice.beam_y;

    this->islice = input_multislice.islice;
    this->dp_Shift = input_multislice.dp_Shift;
  }
		
  // 1:(image_tot, image_coh); 2:(image_tot); 3:(m2psi_tot, m2psi_coh); 4:(m2psi_tot); 
  // 5:(m2psi_tot, psi_coh); 6:(psi_coh); 7:(psi_0); 8:(V); 9:(trans)
  template <typename T>
  void Output_Multislice<T>::set_output_type()
  {
    if (this->is_STEM() || this->is_EELS())
    {
      pimpl->output_type = (this->pn_coh_contrib) ? eTEMOT_image_tot_coh : eTEMOT_image_tot;
    }
    else if (this->is_ISTEM() || this->is_CBED_CBEI() || this->is_PED_HCTEM() ||
      this->is_ED_HRTEM() || this->is_EFTEM())
    {
      pimpl->output_type = (this->pn_coh_contrib) ? eTEMOT_m2psi_tot_coh : eTEMOT_m2psi_tot;
    }
    else if (this->is_EWFS_EWRS())
    {
      pimpl->output_type = (this->is_EWFS_EWRS_SC()) ? eTEMOT_psi_coh : eTEMOT_m2psi_tot_psi_coh;
    }
    else if (this->is_PropFS_PropRS())
    {
      pimpl->output_type = eTEMOT_psi_coh;
    }
    else if (this->is_IWFS_IWRS())
    {
      pimpl->output_type = eTEMOT_psi_0;
    }
    else if (this->is_PPFS_PPRS())
    {
      pimpl->output_type = eTEMOT_V;
    }
    else if (this->is_TFFS_TFRS())
    {
      pimpl->output_type = eTEMOT_trans;
    }
  }

}

#endif
