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

#include <vector>
#include <algorithm>

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "matlab_types.cuh"
#include "atomic_data_mt.hpp"
#include "input_multislice.cuh"
#include "cpu_fcns.hpp"
#include "gpu_fcns.cuh"

namespace mt
{
	inline
		bool is_ot_image_tot_coh(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_image_tot_coh;
	}

	inline
		bool is_ot_image_tot(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_image_tot;
	}

	inline
		bool is_ot_m2psi_tot_coh(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_m2psi_tot_coh;
	}

	inline
		bool is_ot_m2psi_tot(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_m2psi_tot;
	}

	inline
		bool is_ot_m2psi_tot_psi_coh(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_m2psi_tot_psi_coh;
	}

	inline
		bool is_ot_psi_coh(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_psi_coh;
	}

	inline
		bool is_ot_psi_0(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_psi_0;
	}

	inline
		bool is_ot_V(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_V;
	}

	inline
		bool is_ot_trans(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_trans;
	}

	/**************************************************************************************/
	template <class T>
	class Output_Multislice : public Input_Multislice<T>
	{
	public:
		using T_r = T;
		using T_c = complex<T>;

		using TVector_hr = host_vector<T>;
		using TVector_hc = host_vector<complex<T>>;

		using TVector_dr = device_vector<T>;
		using TVector_dc = device_vector<complex<T>>;

		Output_Multislice() : Input_Multislice<T_r>(), output_type(eTEMOT_m2psi_tot), 
		ndetector(0), nx(0), ny(0), dx(0), dy(0), dr(0), n_thk(0), n_thk_d(0) {}

		template <class TOutput_Multislice>
		void assign(TOutput_Multislice &output_multislice)
		{
			assign_input_multislice(output_multislice);

			output_type = output_multislice.output_type;
			ndetector = output_multislice.ndetector;
			nx = output_multislice.nx;
			ny = output_multislice.ny;
			dx = output_multislice.dx;
			dy = output_multislice.dy;
			dr = output_multislice.dr;

			x = output_multislice.x;
			y = output_multislice.y;
			r = output_multislice.r;

			image_tot.resize(output_multislice.image_tot.size());
			for (auto ithk = 0; ithk < output_multislice.image_tot.size(); ithk++)
			{
				image_tot[ithk].image.resize(output_multislice.image_tot[ithk].image.size());
				for (auto idet = 0; idet < output_multislice.image_tot[ithk].image.size(); idet++)
				{
					image_tot[ithk].image[idet] = output_multislice.image_tot[ithk].image[idet];
				}
			}

			image_coh.resize(output_multislice.image_coh.size());
			for (auto ithk = 0; ithk < output_multislice.image_coh.size(); ithk++)
			{
				image_coh[ithk].image.resize(output_multislice.image_coh[ithk].image.size());
				for (auto idet = 0; idet < output_multislice.image_coh[ithk].image.size(); idet++)
				{
					image_coh[ithk].image[idet] = output_multislice.image_coh[ithk].image[idet];
				}
			}

			m2psi_tot.resize(output_multislice.m2psi_tot.size());
			for (auto ithk = 0; ithk < output_multislice.m2psi_tot.size(); ithk++)
			{
				m2psi_tot[ithk] = output_multislice.m2psi_tot[ithk];
			}

			m2psi_coh.resize(output_multislice.m2psi_coh.size());
			for (auto ithk = 0; ithk < output_multislice.m2psi_coh.size(); ithk++)
			{
				m2psi_coh[ithk] = output_multislice.m2psi_coh[ithk];
			}

			psi_coh.resize(output_multislice.psi_coh.size());
			for (auto ithk = 0; ithk < output_multislice.psi_coh.size(); ithk++)
			{
				psi_coh[ithk] = output_multislice.psi_coh[ithk];
			}

			V.resize(output_multislice.V.size());
			for (auto ithk = 0; ithk < output_multislice.V.size(); ithk++)
			{
				V[ithk] = output_multislice.V[ithk];
			}

			trans.resize(output_multislice.trans.size());
			for (auto ithk = 0; ithk < output_multislice.trans.size(); ithk++)
			{
				trans[ithk] = output_multislice.trans[ithk];
			}

			psi_0.resize(output_multislice.psi_0.size());
			for (auto ithk = 0; ithk < output_multislice.psi_0.size(); ithk++)
			{
				psi_0[ithk] = output_multislice.psi_0[ithk];
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

			assign_input_multislice(*input_multislice);

			// set required number of thickness
			n_thk = this->thick.size();
			n_thk_d = 0;

 			thk_gpu.resize(n_thk);
			std::fill(thk_gpu.begin(), thk_gpu.end(), false);

			// check selected device
			auto bb_is_device = this->system_conf.is_device();

			// get available gpu free memory
			double free_memory_mb = get_free_memory<e_device>() - 10;
			int nxy_r = this->output_area.nxy();
			int nxy_g = this->grid_2d.nxy();

			if(bb_is_device)
			{
				psi_zh.resize(nxy_g);
				m2psi_zh.resize(nxy_g);
			}

			ndetector = (this->is_EELS()) ? 1 : this->detector.size();

			set_output_grid();

			set_output_type();

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
				for (auto ithk = 0; ithk < this->thick.size(); ithk++)
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
				for (auto ithk = 0; ithk < this->thick.size(); ithk++)
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
				for (auto ithk = 0; ithk < this->thick.size(); ithk++)
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
				for (auto ithk = 0; ithk < this->thick.size(); ithk++)
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
				for (auto ithk = 0; ithk < this->thick.size(); ithk++)
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
				for (auto ithk = 0; ithk < this->thick.size(); ithk++)
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
				for (auto ithk = 0; ithk < this->thick.size(); ithk++)
				{
					mt::fill(stream, psi_0[ithk], T_c(0));
				}
			}
			break;
			case eTEMOT_V:
			{
				for (auto ithk = 0; ithk < this->thick.size(); ithk++)
				{
					mt::fill(stream, V[ithk], T_r(0));
				}
			}
			break;
			case eTEMOT_trans:
			{
				for (auto ithk = 0; ithk < this->thick.size(); ithk++)
				{
					mt::fill(stream, trans[ithk], T_c(0));
				}
			}
			break;
			}
		}

		void init_psi_coh()
		{
			for (auto ithk = 0; ithk < this->thick.size(); ithk++)
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
				mt::add_scale_shift_2d(this->grid_2d, w, phi, psi_coh_d[ithk], &psi_zh);
			}
			else
			{
				 mt::add_scale_shift_2d(this->grid_2d, w, phi, psi_coh[ithk], &psi_zh);
			}
		}

 		template<class TVector>
		void add_scale_crop_shift_psi_coh(int ithk, T_c w, TVector &phi)
		{
			if(thk_gpu[ithk])
			{
				 mt::add_scale_crop_shift_2d(this->grid_2d, w, phi, this->output_area, psi_coh_d[ithk], &psi_zh);
			}
			else
			{
				 mt::add_scale_crop_shift_2d(this->grid_2d, w, phi, this->output_area, psi_coh[ithk], &psi_zh);
			}
		}

   		template<class TVector>
		void add_scale_crop_shift_m2psi_coh_from_psi(int ithk, T_r w, TVector &phi)
		{
			if(thk_gpu[ithk])
			{
				 mt::add_scale_square_crop_shift_2d(this->grid_2d, w, phi, this->output_area, m2psi_coh_d[ithk], &psi_zh);
			}
			else
			{
				 mt::add_scale_square_crop_shift_2d(this->grid_2d, w, phi, this->output_area, m2psi_coh[ithk], &psi_zh);
			}
		}

 		template<class TVector>
		void add_scale_crop_shift_m2psi_coh_from_m2psi(int ithk, T_r w, TVector &m2phi)
		{
			if(thk_gpu[ithk])
			{
				 mt::add_scale_crop_shift_2d(this->grid_2d, w, m2phi, this->output_area, m2psi_coh_d[ithk], &m2psi_zh);
			}
			else
			{
				 mt::add_scale_crop_shift_2d(this->grid_2d, w, m2phi, this->output_area, m2psi_coh[ithk], &m2psi_zh);
			}
		}

 		template<class TVector>
		void add_scale_crop_shift_m2psi_tot_from_psi(int ithk, T_r w, TVector phi)
		{
			if(thk_gpu[ithk])
			{
				 mt::add_scale_square_crop_shift_2d(this->grid_2d, w, phi, this->output_area, m2psi_tot_d[ithk], &psi_zh);
			}
			else
			{
				 mt::add_scale_square_crop_shift_2d(this->grid_2d, w, phi, this->output_area, m2psi_tot[ithk], &psi_zh);
			}
		}

 		template<class TVector>
		void add_scale_crop_shift_m2psi_tot_from_m2psi(int ithk, T_r w, TVector &m2phi)
		{
			if(thk_gpu[ithk])
			{
				 mt::add_scale_crop_shift_2d(this->grid_2d, w, m2phi, this->output_area, m2psi_tot_d[ithk], &m2psi_zh);
			}
			else
			{
				 mt::add_scale_crop_shift_2d(this->grid_2d, w, m2phi, this->output_area, m2psi_tot[ithk], &m2psi_zh);
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
				 mt::assign_shift_2d(this->grid_2d, phi, psi_coh_d[ithk], &psi_zh);
			}
			else
			{
				 mt::assign_shift_2d(this->grid_2d, phi, psi_coh[ithk], &psi_zh);
			}
		}

  		template<class TVector>
		void set_crop_shift_psi_coh(int ithk, TVector &phi)
		{
			if(thk_gpu[ithk])
			{
				 mt::assign_crop_shift_2d(this->grid_2d, phi, this->output_area, psi_coh_d[ithk], &psi_zh);
			}
			else
			{
				 mt::assign_crop_shift_2d(this->grid_2d, phi, this->output_area, psi_coh[ithk], &psi_zh);
			}
		}

 		template<class TVector>
		void set_crop_shift_m2psi_coh(int ithk, TVector &m2phi)
		{
			if(thk_gpu[ithk])
			{
				 mt::assign_crop_shift_2d(this->grid_2d, m2phi, this->output_area, m2psi_coh_d[ithk], &m2psi_zh);
			}
			else
			{
				 mt::assign_crop_shift_2d(this->grid_2d, m2phi, this->output_area, m2psi_coh[ithk], &m2psi_zh);
			}
		}

 		template<class TVector>
		void set_crop_shift_m2psi_tot(int ithk, TVector &m2phi)
		{
			if(thk_gpu[ithk])
			{
				 mt::assign_crop_shift_2d(this->grid_2d, m2phi, this->output_area, m2psi_tot_d[ithk], &m2psi_zh);
			}
			else
			{
				 mt::assign_crop_shift_2d(this->grid_2d, m2phi, this->output_area, m2psi_tot[ithk], &m2psi_zh);
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
			const int n_thk = this->thick.size();

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
					mt::fft2_shift(stream, this->grid_2d, psi_0[ithk]);
				}
			}
			break;
			case eTEMOT_V:
			{
				for (auto ithk = 0; ithk < n_thk; ithk++)
				{
					mt::fft2_shift(stream, this->grid_2d, V[ithk]);
				}
			}
			break;
			case eTEMOT_trans:
			{
				for (auto ithk = 0; ithk < n_thk; ithk++)
				{
					mt::fft2_shift(stream, this->grid_2d, trans[ithk]);
				}
			}
			break;
			}
		}

		int nxy() const { return nx*ny; }

		/**************************************************************************************/
		inline
		bool is_ot_image_tot_coh() const
		{
			return mt::is_ot_image_tot_coh(output_type);
		}

		inline
		bool is_ot_image_tot() const
		{
			return mt::is_ot_image_tot(output_type);
		}

		inline
		bool is_ot_m2psi_tot_coh() const
		{
			return mt::is_ot_m2psi_tot_coh(output_type);
		}

		inline
		bool is_ot_m2psi_tot() const
		{
			return mt::is_ot_m2psi_tot(output_type);
		}

		inline
		bool is_ot_m2psi_tot_psi_coh() const
		{
			return mt::is_ot_m2psi_tot_psi_coh(output_type);
		}

		inline
		bool is_ot_psi_coh() const
		{
			return mt::is_ot_psi_coh(output_type);
		}

		inline
		bool is_ot_psi_0() const
		{
			return mt::is_ot_psi_0(output_type);
		}

		inline
		bool is_ot_V() const
		{
			return mt::is_ot_V(output_type);
		}

		inline
		bool is_ot_trans() const
		{
			return mt::is_ot_trans(output_type);
		}

		host_vector<float> extract_data(ePhonon_Model_Output fp_ctr, eShow_CData show_data, int ithk, int idet = 0)
		{
			host_vector<float> data(nxy());

			switch (output_type)
			{
			case eTEMOT_image_tot_coh:
			{
				data = (fp_ctr == eFMO_Total) ? image_tot[ithk].image[idet] : image_coh[ithk].image[idet];
			}
			break;
			case eTEMOT_image_tot:
			{
				data = image_tot[ithk].image[idet];
			}
			break;
			case eTEMOT_m2psi_tot_coh:
			{
				data = (fp_ctr == eFMO_Total) ? m2psi_tot[ithk] : m2psi_coh[ithk];
			}
			break;
			case eTEMOT_m2psi_tot:
			{
				data = m2psi_tot[ithk];
			}
			break;
			case eTEMOT_m2psi_tot_psi_coh:
			{
				if (fp_ctr == eFMO_Total)
				{
					data = m2psi_tot[ithk];
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
				data = V[ithk];
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

		Stream<e_host> stream;
	private:
		Vector<T_c, e_host> psi_zh;
		Vector<T_r, e_host> m2psi_zh;

		int n_thk;
		int n_thk_d;

 		template <class U>
		int cal_n_thk_a(double memory, int nxy)
		{
			return static_cast<int>(floor(memory/mt::sizeMb<U>(nxy)));
		}

		void set_output_grid()
		{
			if (this->is_STEM() || this->is_EELS())
			{
 				nx = this->scanning.nx;
				ny = this->scanning.ny;

				dx = this->scanning.dRx;
				dy = this->scanning.dRy;

				x = this->scanning.x;
				y = this->scanning.y;
				r = this->scanning.r;
			}
			else
			{
				nx = this->output_area.nx();
				ny = this->output_area.ny();

				const bool is_RS = this->is_grid_RS();
				dx = (is_RS)?this->grid_2d.dRx:this->grid_2d.dgx;
				dy = (is_RS)?this->grid_2d.dRy:this->grid_2d.dgy;

				x.resize(nx);
				y.resize(ny);

				for (auto ix = this->output_area.ix_0; ix < this->output_area.ix_e; ix++)
				{
					int ix_s = ix-this->output_area.ix_0;
					x[ix_s] = (is_RS)?this->grid_2d.Rx(ix):this->grid_2d.gx(ix);
				}

				for (auto iy = this->output_area.iy_0; iy < this->output_area.iy_e; iy++)
				{
					int iy_s = iy-this->output_area.iy_0;
					y[iy_s] = (is_RS)?this->grid_2d.Ry(iy):this->grid_2d.gy(iy);
				}
			}
		}

		// 1:(image_tot, image_coh); 2:(image_tot); 3:(m2psi_tot, m2psi_coh); 4:(m2psi_tot); 
		// 5:(m2psi_tot, psi_coh); 6:(psi_coh); 7:(psi_0); 8:(V); 9:(trans)
		void set_output_type()
		{
			if (this->is_STEM() || this->is_EELS())
			{
				output_type = (this->pn_coh_contrib) ? eTEMOT_image_tot_coh : eTEMOT_image_tot;
			}
			else if (this->is_ISTEM() || this->is_CBED_CBEI() || this->is_PED_HCTEM() ||
				this->is_ED_HRTEM() || this->is_EFTEM())
			{
				output_type = (this->pn_coh_contrib) ? eTEMOT_m2psi_tot_coh : eTEMOT_m2psi_tot;
			}
			else if (this->is_EWFS_EWRS())
			{
				output_type = (this->is_EWFS_EWRS_SC()) ? eTEMOT_psi_coh : eTEMOT_m2psi_tot_psi_coh;
			}
			else if (this->is_PropFS_PropRS())
			{
				output_type = eTEMOT_psi_coh;
			}
			else if (this->is_IWFS_IWRS())
			{
				output_type = eTEMOT_psi_0;
			}
			else if (this->is_PPFS_PPRS())
			{
				output_type = eTEMOT_V;
			}
			else if (this->is_TFFS_TFRS())
			{
				output_type = eTEMOT_trans;
			}
		}

		template <class TInput_Multislice>
		void assign_input_multislice(TInput_Multislice &input_multislice)
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
	};

} // namespace mt

#endif