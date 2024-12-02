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

#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H

#include "math.cuh"
#include "types.cuh"
#include "fft.cuh"
#include "input_multislice.cuh"
#include "cpu_fcns.hpp"
#include "gpu_fcns.cuh"
#include "cgpu_fcns.cuh"
#include "transmission_function.cuh"
#include "incident_wave.cuh"
#include "propagator.cuh"
#include "microscope_effects.cuh"

namespace mt
{
	template <class T, eDevice dev>
	class Wave_Function: public Transmission_Function<T, dev>
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T_r, dev>;
			using TVector_c = Vector<T_c, dev>;
			using size_type = std::size_t;

			Wave_Function(): Transmission_Function<T, dev>(){}

			void set_input_data(Input_Multislice<T_r> *input_multislice_i, Stream<dev> *stream_i, FFT<T_r, dev> *fft2_i)
			{
				psi_z.resize(input_multislice_i->grid_2d.nxy());
				m2psi_z.resize(input_multislice_i->grid_2d.nxy());

				if(input_multislice_i->is_STEM())
				{
					detector.assign(input_multislice_i->detector);
					if(input_multislice_i->is_detector_matrix())
					{
						for(auto i = 0; i<detector.size(); i++)
						{
							mt::fft2_shift(*stream_i, input_multislice_i->grid_2d, detector.fR[i]);
						}
					}
				}

				incident_wave.set_input_data(input_multislice_i, stream_i, fft2_i);

				propagator.set_input_data(input_multislice_i, stream_i, fft2_i);

				if(input_multislice_i->is_ISTEM_CBEI_HRTEM_HCTEM_EFTEM())
				{
					microscope_effects.set_input_data(input_multislice_i, stream_i, fft2_i);
				} 

				Transmission_Function<T, dev>::set_input_data(input_multislice_i, stream_i, fft2_i);
			}

			void phase_multiplication(const T_r &gxu, const T_r &gyu, TVector_c &psi_i, TVector_c &psi_o, bool scaling=true)
			{
				if(this->input_multislice->dp_Shift || isZero(gxu, gyu))
				{
					if (psi_i.data() != psi_o.data())
					{
						psi_o.assign(psi_i.begin(), psi_i.end());
					}
					return;
				}

				mt::exp_r_factor_2d(*(this->stream), this->input_multislice->grid_2d, c_2Pi*gxu, c_2Pi*gyu, psi_i, psi_o, scaling);
			}

			void phase_multiplication(const T_r &gxu, const T_r &gyu, TVector_c &psi_io, bool scaling=true)
			{
				phase_multiplication(gxu, gyu, psi_io, psi_io, scaling);
			}

			TVector_c* get_psi(const eSpace &space, const T_r &gxu, const T_r &gyu, 
			T_r z, TVector_c &psi_i)
			{
				TVector_c *psi_o = &(this->trans_0);
				// real space not need to include scaling and phase multiplication oposite sign as the propagation
				phase_multiplication(-gxu, -gyu, psi_i, *psi_o, false); 
				propagator(space, gxu, gyu, z, *psi_o);

				return psi_o;
			}

			T_r integrated_intensity_over_det(T_r w_i, const int &iDet, TVector_c &psi_z)
			{
				T_r int_val = 0;
				switch (detector.type)
				{
					case mt::eDT_Circular:
					{
						auto g_inner = detector.g_inner[iDet];
						auto g_outer = detector.g_outer[iDet];
							
						int_val = w_i*mt::sum_square_over_Det(*(this->stream), this->input_multislice->grid_2d, g_inner, g_outer, psi_z);
					}
					break;
					case mt::eDT_Radial:
					{
						int_val = 0;
					}
					break;
					case mt::eDT_Matrix:
					{
						int_val = w_i*mt::sum_square_over_Det(*(this->stream), this->input_multislice->grid_2d, detector.fR[iDet], psi_z);
					}
					break;
				}

				return int_val;
			}

			template <class TOutput_multislice>
			void set_m2psi_tot_psi_coh(TVector_c &psi_z_i, const T_r &gxu, const T_r &gyu, 
			const int &islice, const T_r &w_i, TOutput_multislice &output_multislice)
			{
				int ithk = this->slicing.slice[islice].ithk;
				if(0 <= ithk)
				{
					auto *psi_z_o = get_psi(this->input_multislice->get_simulation_space(), gxu, gyu, this->slicing.thick[ithk].z_back_prop, psi_z_i);

					if(this->input_multislice->is_STEM())
					{
						for(auto iDet = 0; iDet<detector.size(); iDet++)
						{
							int iscan = this->input_multislice->iscan[0];
							output_multislice.image_tot[ithk].image[iDet][iscan] += integrated_intensity_over_det(w_i, iDet, *psi_z_o);
						}

						if(this->input_multislice->pn_coh_contrib)
						{
							output_multislice.add_scale_psi_coh(ithk, w_i, *psi_z_o);
						}
					}
					else if(this->input_multislice->is_EWFS_EWRS_SC())
					{
						output_multislice.set_crop_shift_psi_coh(ithk, *psi_z_o);
					}
					else if(this->input_multislice->is_EWFS_EWRS())
					{
						output_multislice.add_scale_crop_shift_m2psi_tot_from_psi(ithk, w_i, *psi_z_o);
						output_multislice.add_scale_crop_shift_psi_coh(ithk, w_i, *psi_z_o);
					}
					else if(this->input_multislice->is_ISTEM_CBEI_HRTEM_HCTEM_EFTEM())
					{
						microscope_effects(*psi_z_o, m2psi_z);
						output_multislice.add_scale_crop_shift_m2psi_tot_from_m2psi(ithk, w_i, m2psi_z);

						if(this->input_multislice->pn_coh_contrib)
						{
							output_multislice.add_scale_psi_coh(ithk, w_i, *psi_z_o);
						}
					}
					else
					{
						output_multislice.add_scale_crop_shift_m2psi_tot_from_psi(ithk, w_i, *psi_z_o);

						if(this->input_multislice->pn_coh_contrib)
						{
							output_multislice.add_scale_psi_coh(ithk, w_i, *psi_z_o);
						}
					}
				}
			}

			template <class TOutput_multislice>
			void set_m2psi_coh(TOutput_multislice &output_multislice)
			{
				if(!this->input_multislice->pn_coh_contrib || this->input_multislice->is_EWFS_EWRS())
				{
					return;
				}

				int n_thk = this->input_multislice->thick.size();

				if(this->input_multislice->is_STEM())
				{
					for(auto ithk = 0; ithk < n_thk; ithk++)
					{
						output_multislice.from_psi_coh_2_phi(ithk, psi_z);
						for(auto iDet = 0; iDet<detector.size(); iDet++)
						{
							int iscan = this->input_multislice->iscan[0];
							output_multislice.image_coh[ithk].image[iDet][iscan] = integrated_intensity_over_det(1, iDet, psi_z);
						}
					}
				}
				else if(this->input_multislice->is_ISTEM_CBEI_HRTEM_HCTEM_EFTEM())
				{
					for(auto ithk = 0; ithk < n_thk; ithk++)
					{
						output_multislice.from_psi_coh_2_phi(ithk, psi_z);
						microscope_effects(psi_z, m2psi_z);
						output_multislice.set_crop_shift_m2psi_coh(ithk, m2psi_z);
					}
				}
				else
				{
					for(auto ithk = 0; ithk < n_thk; ithk++)
					{
						output_multislice.from_psi_coh_2_phi(ithk, psi_z);
						output_multislice.add_scale_crop_shift_m2psi_coh_from_psi(ithk, 1.0, psi_z);
					}
				}
			}

			template <class TVector_c>
			void psi_slice(const T_r &gxu, const T_r &gyu, const int &islice, TVector_c &psi_z)
			{
				this->transmit(islice, psi_z);

				if(this->input_multislice->is_multislice())
				{
					propagator(eS_Real, gxu, gyu, this->dz(islice), psi_z);
				}
			}

			template <class TOutput_multislice>
			void psi(T_r w_i, TVector_c &psi_z, TOutput_multislice &output_multislice)
			{
				T_r gx_0 = this->input_multislice->gx_0();
				T_r gy_0 = this->input_multislice->gy_0();

				for(auto islice = 0; islice<this->slicing.slice.size(); islice++)
				{
					psi_slice(gx_0, gy_0, islice, psi_z);

					set_m2psi_tot_psi_coh(psi_z, gx_0, gy_0, islice, w_i, output_multislice);
				}
			}

			template <class TOutput_multislice>
			void psi(int islice_0, int islice_e, T_r w_i, TVector_c &trans, TOutput_multislice &output_multislice)
			{
				int ithk = this->slicing.slice[islice_e].ithk;
				if(0 <= ithk)
				{
					T_r gx_0 = this->input_multislice->gx_0();
					T_r gy_0 = this->input_multislice->gy_0();

					if(this->input_multislice->eels_fr.is_Single_Channelling())
					{
						T_r dz = this->dz_m(islice_0, islice_e);
						propagator(eS_Real, gx_0, gy_0, dz, psi_z);
					}
					else if(this->input_multislice->eels_fr.is_Mixed_Channelling())
					{
						T_r dz = 0.5*this->dz_m(islice_0, islice_e);
						propagator(eS_Real, gx_0, gy_0, dz, psi_z);
						mt::multiply(*(this->stream), trans, psi_z);
						propagator(eS_Real, gx_0, gy_0, dz, psi_z);
					}
					else if(this->input_multislice->eels_fr.is_Double_Channelling())
					{
						for(auto islice = islice_0; islice<= islice_e; islice++)
						{
							psi_slice(gx_0, gy_0, islice, psi_z);
						}
					}

					phase_multiplication(gx_0, gy_0, psi_z);
					propagator(eS_Reciprocal, gx_0, gy_0, this->slicing.thick[ithk].z_back_prop, psi_z);

					if(this->input_multislice->is_EELS())
					{
						int iscan = this->input_multislice->iscan[0];
						output_multislice.image_tot[ithk].image[0][iscan] += w_i*mt::sum_square_over_Det(*(this->stream), this->input_multislice->grid_2d, 0, this->input_multislice->eels_fr.g_collection, psi_z);
					}
					else
					{
						mt::hard_aperture(*(this->stream), this->input_multislice->grid_2d, this->input_multislice->eels_fr.g_collection, 1.0, psi_z);
						microscope_effects(psi_z, m2psi_z);
						output_multislice.add_scale_crop_shift_m2psi_tot_from_m2psi(ithk, w_i, m2psi_z);
					}
				}
			}

			void set_incident_wave(TVector_c &psi, Vector<T_r, e_host> &beam_x, Vector<T_r, e_host> &beam_y)
			{
				T_r gxu = 0;
				T_r gyu = 0;
				auto z_init = this->slicing.z_m(0);

				incident_wave(psi, gxu, gyu, beam_x, beam_y, z_init);
			}

			void set_incident_wave(TVector_c &psi)
			{
				auto &beam_x = this->input_multislice->beam_x;
				auto &beam_y = this->input_multislice->beam_y;

				set_incident_wave(psi, beam_x, beam_y);
			}

			Propagator<T_r, dev> propagator;

			TVector_c psi_z;
			TVector_r m2psi_z;

			Detector<T_r, dev> detector; 	
			Microscope_Effects<T_r, dev> microscope_effects;
			Incident_Wave<T_r, dev> incident_wave;		

			//mt::Timing<mt::e_device> time;
	};

} // namespace mt

#endif