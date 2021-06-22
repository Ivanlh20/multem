/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * Multem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version of the License, or
 * (at your option) any later version.
 *
 * Multem is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Multem. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H

#include "math.cuh"
#include "types.cuh"
#include "cgpu_fft.cuh"
#include "in_classes.cuh"
#include "cpu_fcns.hpp"
#include "gpu_fcns.cuh"
#include "cgpu_fcns.cuh"
#include "transmission_fcn.cuh"
#include "incident_wave.cuh"
#include "propagator.cuh"
#include "microscope_effects.cuh"

#include "timing.cuh"

#include <mex.h>

namespace mt
{
	template <class T, eDev Dev>
	class Wave_Function: public Transmission_Fcn<T, Dev>
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVctr_r = Vctr<T_r, Dev>;
			using TVctr_c = Vctr<T_c, Dev>;
			using size_type = dt_uint64;

			Wave_Function(): Transmission_Fcn<T, Dev>() {}

			void set_in_data(In_Multem<T_r> *in_multem_i, Stream<Dev> *stream_i, FFT<T_r, Dev> *fft2_i)
			{
				psi_z.resize(in_multem_i->grid_2d.size());
				m2psi_z.resize(in_multem_i->grid_2d.size());

				if (in_multem_i->is_STEM())
				{
					detector.assign(in_multem_i->detector);
					if (in_multem_i->is_detector_matrix())
					{
						for(auto i = 0; i<detector.size(); i++)
						{
							mt::fcn_fftsft_2d(*stream_i, in_multem_i->grid_2d, detector.fR[i]);
						}
					}
				}

				incident_wave.set_in_data(in_multem_i, stream_i, fft2_i);

				propagator.set_in_data(in_multem_i, stream_i, fft2_i);

				if (in_multem_i->is_ISTEM_CBEI_HRTEM_HCTEM_EFTEMRS())
				{
					microscope_effects.set_in_data(in_multem_i, stream_i, fft2_i);
				} 

				Transmission_Fcn<T, Dev>::set_in_data(in_multem_i, stream_i, fft2_i);
			}

			void phase_multiplication(const T_r &gxu, const T_r &gyu, TVctr_c& psi_i, TVctr_c& psi_o)
			{
				if (this->in_multem->dp_Shift || fcn_is_zero(gxu, gyu))
				{
					if (psi_i.data() != psi_o.data())
					{
						psi_o.assign(psi_i.begin(), psi_i.end());
					}
					return;
				}

				mt::fcn_rs_exp_factor_2d(*(this->stream), this->in_multem->grid_2d, c_2pi<T>*gxu, c_2pi<T>*gyu, psi_i, psi_o);
			}

			void phase_multiplication(const T_r &gxu, const T_r &gyu, TVctr_c& psi_io)
			{
				phase_multiplication(gxu, gyu, psi_io, psi_io);
			}

			TVctr_c* get_psi(const eSpace &space, const T_r &gxu, const T_r &gyu, 
			T_r z, TVctr_c& psi_i)
			{
				TVctr_c *psi_o = &(this->trans_0);
				phase_multiplication(gxu, gyu, psi_i, *psi_o);
				propagator(space, gxu, gyu, z, *psi_o);

				return psi_o;
			}

			T_r integrated_intensity_over_det(T_r w_i, const dt_int32& iDet, TVctr_c& psi_z)
			{
				T_r int_val = 0;
				switch (detector.type)
				{
					case mt::edt_circular:
					{
						auto g_inner = detector.g_inner[iDet];
						auto g_outer = detector.g_outer[iDet];
							
						int_val = w_i*mt::fcn_int_det_ring_norm_2(*(this->stream), this->in_multem->grid_2d, g_inner, g_outer, psi_z);
					}
					break;
					case mt::edt_radial:
					{
						int_val = 0;
					}
					break;
					case mt::edt_matrix:
					{
						int_val = w_i*mt::fcn_int_det_ring_norm_2(*(this->stream), this->in_multem->grid_2d, detector.fR[iDet], psi_z);
					}
					break;
				}

				return int_val;
			}

			template <class TOutput_multislice>
			void set_m2psi_tot_psi_coh(TVctr_c& psi_z_i, const T_r &gxu, const T_r &gyu, 
			const dt_int32& islice, const T_r &w_i, TOutput_multislice &output_multem)
			{
				dt_int32 ithk = this->slicing.slice[islice].ithk;
				if (0 <= ithk)
				{
					auto *psi_z_o = get_psi(this->in_multem->get_simulation_space(), gxu, gyu, this->slicing.thick[ithk].z_back_prop, psi_z_i);

					if (this->in_multem->is_STEM())
					{
						for(auto iDet = 0; iDet<detector.size(); iDet++)
						{
							dt_int32 ibeam = this->in_multem->ibeam[0];
							output_multem.image_tot(ithk, iDet, ibeam) += integrated_intensity_over_det(w_i, iDet, *psi_z_o);
						}

						if (this->in_multem->atomic_vib.coh_contrib)
						{
							output_multem.add_sc_psi_coh(ithk, 0, w_i, *psi_z_o);
						}
					}
					else if (this->in_multem->is_EWFS_EWRS_SC())
					{
						output_multem.set_crop_sft_psi_coh(ithk, *psi_z_o);
					}
					else if (this->in_multem->is_EWFS_EWRS())
					{
						output_multem.add_sc_crop_sft_m2psi_tot_from_psi(ithk, w_i, *psi_z_o);
						output_multem.add_sc_crop_sft_psi_coh(ithk, w_i, *psi_z_o);
					}
					else if (this->in_multem->is_CBED())
					{
						dt_int32 ithk_beam = (this->in_multem->thick.size()==1)?this->in_multem->ibeam[0]:ithk;

						output_multem.add_sc_crop_sft_m2psi_tot_from_psi(ithk_beam, w_i, *psi_z_o);

 						if (this->in_multem->atomic_vib.coh_contrib)
						{
							output_multem.add_sc_psi_coh(ithk_beam, w_i, *psi_z_o);
						}
					}
					else if (this->in_multem->is_CBEI())
					{
						dt_int32 ithk_beam = (this->in_multem->thick.size()==1)?this->in_multem->ibeam[0]:ithk;

						microscope_effects(*psi_z_o, m2psi_z);
						output_multem.add_sc_crop_sft_m2psi_tot_from_m2psi(ithk_beam, w_i, m2psi_z);

 						if (this->in_multem->atomic_vib.coh_contrib)
						{
							output_multem.add_sc_psi_coh(ithk_beam, w_i, *psi_z_o);
						}
					}
					else if (this->in_multem->is_ISTEM_CBEI_HRTEM_HCTEM_EFTEMRS())
					{
						microscope_effects(*psi_z_o, m2psi_z);
						output_multem.add_sc_crop_sft_m2psi_tot_from_m2psi(ithk, w_i, m2psi_z);

						if (this->in_multem->atomic_vib.coh_contrib)
						{
							output_multem.add_sc_psi_coh(ithk, w_i, *psi_z_o);
						}
					}
					else
					{
						output_multem.add_sc_crop_sft_m2psi_tot_from_psi(ithk, w_i, *psi_z_o);

						if (this->in_multem->atomic_vib.coh_contrib)
						{
							output_multem.add_sc_psi_coh(ithk, w_i, *psi_z_o);
						}
					}

					// this->stream->synchronize();
				}
			}

			template <class TOutput_multislice>
			void set_m2psi_coh(TOutput_multislice &output_multem)
			{
				if (!this->in_multem->atomic_vib.coh_contrib || this->in_multem->is_EWFS_EWRS())
				{
					return;
				}

				dt_int32 n_thk = this->in_multem->thick.size();

				if (this->in_multem->is_STEM())
				{
					for(auto ithk = 0; ithk < n_thk; ithk++)
					{
						output_multem.from_psi_coh_2_phi(ithk, psi_z);
						for(auto iDet = 0; iDet<detector.size(); iDet++)
						{
							dt_int32 ibeam = this->in_multem->ibeam[0];
							output_multem.image_coh(ithk, iDet, ibeam) = integrated_intensity_over_det(1, iDet, psi_z);
						}
					}
				}
				else if (this->in_multem->is_CBED())
				{
					n_thk = (n_thk==1)?this->in_multem->number_of_beams():n_thk;

					for(auto ithk = 0; ithk < n_thk; ithk++)
					{
 						output_multem.from_psi_coh_2_phi(ithk, psi_z);
						output_multem.add_sc_crop_sft_m2psi_coh_from_psi(ithk, 1.0, psi_z);
					}
				}
 				else if (this->in_multem->is_CBEI())
				{
					n_thk = (n_thk==1)?this->in_multem->number_of_beams():n_thk;

					for(auto ithk = 0; ithk < n_thk; ithk++)
					{
						output_multem.from_psi_coh_2_phi(ithk, psi_z);
						microscope_effects(psi_z, m2psi_z);
						output_multem.set_crop_sft_m2psi_coh(ithk, m2psi_z);
					}
				}
				else if (this->in_multem->is_ISTEM_CBEI_HRTEM_HCTEM_EFTEMRS())
				{
					for(auto ithk = 0; ithk < n_thk; ithk++)
					{
						output_multem.from_psi_coh_2_phi(ithk, psi_z);
						microscope_effects(psi_z, m2psi_z);
						output_multem.set_crop_sft_m2psi_coh(ithk, m2psi_z);
					}
				}
				else
				{
					for(auto ithk = 0; ithk < n_thk; ithk++)
					{
						output_multem.from_psi_coh_2_phi(ithk, psi_z);
						output_multem.add_sc_crop_sft_m2psi_coh_from_psi(ithk, 1.0, psi_z);
					}
				}

				// this->stream->synchronize();
			}

			template <class TVctr_c>
			void psi_slice(const T_r &gxu, const T_r &gyu, const dt_int32& islice, TVctr_c& psi_z)
			{
				this->transmit(islice, psi_z);

				if (this->in_multem->is_multislice())
				{
					// time.tic();
					propagator(esp_real, gxu, gyu, this->sli_thk(islice), psi_z);
					// time.toc();
					// mexPrintf("time = %7.5f\n", time.elapsed_ms());
				}
			}

			template <class TOutput_multislice>
			void psi(T_r w_i, TVctr_c& psi_z, TOutput_multislice &output_multem)
			{
				T_r gx_0 = this->in_multem->gx_0();
				T_r gy_0 = this->in_multem->gy_0();

				for(auto islice = 0; islice<this->slicing.slice.size(); islice++)
				{
					psi_slice(gx_0, gy_0, islice, psi_z);

					set_m2psi_tot_psi_coh(psi_z, gx_0, gy_0, islice, w_i, output_multem);
				}
			}

			template <class TOutput_multislice>
			void psi(dt_int32 islice_0, dt_int32 islice_e, T_r w_i, TVctr_c& trans, TOutput_multislice &output_multem)
			{
				dt_int32 ithk = this->slicing.slice[islice_e].ithk;
				if (0 <= ithk)
				{
					T_r gx_0 = this->in_multem->gx_0();
					T_r gy_0 = this->in_multem->gy_0();

					if (this->in_multem->eels_fr.is_Single_Chan())
					{
						T_r sli_thk = this->dz_m(islice_0, islice_e);
						propagator(esp_real, gx_0, gy_0, sli_thk, psi_z);
					}
					else if (this->in_multem->eels_fr.is_Mixed_Chan())
					{
						T_r sli_thk = 0.5*this->dz_m(islice_0, islice_e);
						propagator(esp_real, gx_0, gy_0, sli_thk, psi_z);
						mt::ew_mult(*(this->stream), trans, psi_z);
						propagator(esp_real, gx_0, gy_0, sli_thk, psi_z);
					}
					else if (this->in_multem->eels_fr.is_Double_Chan())
					{
						for(auto islice = islice_0; islice<= islice_e; islice++)
						{
							psi_slice(gx_0, gy_0, islice, psi_z);
						}
					}

					phase_multiplication(gx_0, gy_0, psi_z);
					propagator(esp_fourier, gx_0, gy_0, this->slicing.thick[ithk].z_back_prop, psi_z);

					if (this->in_multem->is_STEM_ISTEM_EELS())
					{
						dt_int32 ibeam = this->in_multem->ibeam[0];
						output_multem.image_tot(ithk, 0, ibeam) += w_i*mt::fcn_int_det_ring_norm_2(*(this->stream), this->in_multem->grid_2d, 0, this->in_multem->eels_fr.g_coll, psi_z);
					}
					if (this->in_multem->is_EFTEMRS())
					{
						microscope_effects(psi_z, m2psi_z);
						output_multem.add_sc_crop_sft_m2psi_tot_from_m2psi(ithk, w_i, m2psi_z);
					}
					else
					{
		 output_multem.add_sc_crop_sft_m2psi_tot_from_psi(ithk, w_i, psi_z);
					}
				}
			}

			void set_incident_wave(TVctr_c& psi, Beam_Pos_2d<T_r>& beam_pos)
			{
				T_r gxu = 0;
				T_r gyu = 0;
				auto z_init = this->slicing.z_m(0);

				incident_wave(psi, gxu, gyu, beam_pos, z_init);
			}

			void set_incident_wave(TVctr_c& psi)
			{
				auto &beam_pos = this->in_multem->beam_pos;

				set_incident_wave(psi, beam_pos);
			}

			Propagator<T_r, Dev> propagator;

			TVctr_c psi_z;
			TVctr_r m2psi_z;

			Detector<T_r, Dev> detector;
			Microscope_Effects<T_r, Dev> microscope_effects;
			Incident_Wave<T_r, Dev> incident_wave;

			mt::Timing<mt::edev_gpu> time;
	};

}

#endif