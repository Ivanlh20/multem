/*
 * This file is part of Multem.
 * Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef TEM_SIMULATION_H
#define TEM_SIMULATION_H

#include <fftw3.h>
#include "math.cuh"
#include "types.cuh"
#include "type_traits_gen.cuh"
#include "in_classes.cuh"
#include "output_multem.hpp"
#include "cpu_fcns.hpp"
#include "gpu_fcns.cuh"
#include "cgpu_fcns.cuh"
#include "energy_loss.cuh"
#include "wave_function.cuh"
#include "timing.cuh"

namespace mt
{
	template <class T, eDev Dev>
	class Tem_Simulation
	{
		public:
			using T_r = T;
			using T_c = complex<T>;

			static const eDev device = Dev;

			static dt_bool ext_stop_sim;
			static dt_int32 ext_niter;
			static dt_int32 ext_iter;

			Tem_Simulation(): in_multem(nullptr) {}

			Tem_Simulation(In_Multem<T> *in_multem_i)
			{
				set_in_data(in_multem_i);
			}

			void set_in_data(In_Multem<T> *in_multem_i)
			{
				in_multem = in_multem_i;

				stream.resize(in_multem->system_config.n_stream);

				fft_2d.create_plan_2d(in_multem->grid_2d.ny, in_multem->grid_2d.nx, in_multem->system_config.n_stream);

				if (in_multem->is_EELS_EFTEM())
				{
					energy_loss.set_in_data(in_multem, &stream, &fft_2d);
					psi_thk.resize(in_multem->grid_2d.size());
					if (in_multem->eels_fr.is_Mixed_Chan())
					{
						trans_thk.resize(in_multem->grid_2d.size());
					}
				}

				wave_function.set_in_data(in_multem, &stream, &fft_2d);
			}

			template <class TOutput_multislice>
			void operator()(TOutput_multislice &output_multem)
			{
				if (in_multem->is_STEM_ISTEM())
				{
					STEM_ISTEM(output_multem);
				}
				else if (in_multem->is_CBED_CBEI())
				{
					CBED_CBEI(output_multem);
				}
				else if (in_multem->is_ED_HRTEM())
				{
					ED_HRTEM(output_multem);
				}
				else if (in_multem->is_PED_HCTEM())
				{
					PED_HCTEM(output_multem);
				}
				else if (in_multem->is_EWFS_EWRS())
				{
					EWFS_EWRS(output_multem);
				}
				else if (in_multem->is_EELS_EFTEM())
				{
					EELS_EFTEM(output_multem);
				}

				stream.synchronize();

				output_multem.gather();
				output_multem.clean_temporal();
			}

			void cleanup()
			{
				psi_thk.clear();
				trans_thk.clear();

				fft_2d.cleanup();
				stream.cleanup();
			}

		private:
			template <class TOutput_multislice>
			void STEM_ISTEM(TOutput_multislice &output_multem)
			{
				ext_niter = in_multem->scanning.size()*in_multem->number_pn_conf();
				ext_iter = 0;
				/***************************************************************************************/

				T_r w = in_multem->get_phonon_rot_weight();

				output_multem.init();

				in_multem->ibeam.resize(1);
				in_multem->beam_x.resize(1);
				in_multem->beam_y.resize(1);

				
				if (in_multem->is_STEM() && in_multem->phonon_par.coh_contrib)
				{
					for(auto ibeam = 0; ibeam < in_multem->scanning.size(); ibeam++)
					{	
						output_multem.init_psi_coh();
						in_multem->ibeam[0] = ibeam;
						in_multem->set_iscan_beam_position();
						for(auto iconf = in_multem->phonon_par.iconf_0; iconf <= in_multem->phonon_par.nconf; iconf++)
						{
							wave_function.move_atoms(iconf);
							wave_function.set_incident_wave(wave_function.psi_z);
							wave_function.psi(w, wave_function.psi_z, output_multem);

							ext_iter++;
							if (ext_stop_sim) break;
						}
						wave_function.set_m2psi_coh(output_multem);

						if (ext_stop_sim) break;
					}
				}
				else
				{
					for(auto iconf = in_multem->phonon_par.iconf_0; iconf <= in_multem->phonon_par.nconf; iconf++)
					{
						wave_function.move_atoms(iconf);
						for(auto ibeam = 0; ibeam < in_multem->scanning.size(); ibeam++)
						{
							in_multem->ibeam[0] = ibeam;
							in_multem->set_iscan_beam_position();
							wave_function.set_incident_wave(wave_function.psi_z);
							wave_function.psi(w, wave_function.psi_z, output_multem);

							ext_iter++;
							if (ext_stop_sim) break;
						}
						if (ext_stop_sim) break;
					}

					wave_function.set_m2psi_coh(output_multem);
				}
			}

			template <class TOutput_multislice>
			void CBED_CBEI(TOutput_multislice &output_multem)
			{
				output_multem.init();

				in_multem->ibeam.resize(1);
				in_multem->beam_x.resize(1);
				in_multem->beam_y.resize(1);

				Quad_Coef_1d<dt_float64, edev_cpu> qt;
				Quad_Coef_2d<dt_float64, edev_cpu> qs;

				// Load quadratures
				cond_lens_temporal_spatial_quadratures(in_multem->cond_lens, qt, qs);

				/***************************************************************************************/
				dt_float64 w_pr_0 = in_multem->get_phonon_rot_weight();
				dt_float64 c_10_0 = in_multem->cond_lens.c_10;
				const dt_int32 n_beams = in_multem->number_of_beams();

				ext_niter = qs.size()*qt.size()*in_multem->number_pn_conf();
				ext_iter = 0;

				for(auto iconf = in_multem->phonon_par.iconf_0; iconf <= in_multem->phonon_par.nconf; iconf++)
				{
					wave_function.move_atoms(iconf);

					for(auto ibeam=0; ibeam<n_beams; ibeam++)
					{
						in_multem->ibeam[0] = ibeam;
						in_multem->beam_x[0] = in_multem->beam_x[ibeam];
						in_multem->beam_y[0] = in_multem->beam_y[ibeam];

						// spatial incoherence
						for(auto ispat = 0; ispat<qs.size(); ispat++)
						{
							auto beam_x = qs.x[ispat] + in_multem->beam_x[0];
							auto beam_y = qs.y[ispat] + in_multem->beam_y[0];

							// temporal incoherence
							for(auto itemp = 0; itemp<qt.size(); itemp++)
							{
								auto c_10 = c_10_0 + qt.x[itemp];
								auto w = w_pr_0*qs.w[ispat]*qt.w[itemp];

								in_multem->cond_lens.set_defocus(c_10);
								wave_function.set_incident_wave(wave_function.psi_z, beam_x, beam_y);
								wave_function.psi(w, wave_function.psi_z, output_multem);

								ext_iter++;
								if (ext_stop_sim) break;
							}
						}
					}

					if (ext_stop_sim) break;
				}
				wave_function.set_m2psi_coh(output_multem);

				in_multem->cond_lens.set_defocus(c_10_0);
				in_multem->set_beam_position(in_multem->beam_x, in_multem->beam_y);
			}

			template <class TOutput_multislice>
			void ED_HRTEM(TOutput_multislice &output_multem)
			{
				EWFS_EWRS(output_multem);
			}

			template <class TOutput_multislice>
			void PED_HCTEM(TOutput_multislice &output_multem)
			{
				ext_niter = in_multem->nrot*in_multem->number_pn_conf();
				ext_iter = 0;
				/***************************************************************************************/

				T_r w = in_multem->get_phonon_rot_weight();

				output_multem.init();

				for(auto iconf = in_multem->phonon_par.iconf_0; iconf <= in_multem->phonon_par.nconf; iconf++)
				{
					wave_function.move_atoms(iconf);
					for(auto irot = 0; irot < in_multem->nrot; irot++)
					{
						in_multem->set_phi(irot);
						wave_function.set_incident_wave(wave_function.psi_z);
						wave_function.psi(w, wave_function.psi_z, output_multem);

						ext_iter++;
						if (ext_stop_sim) break;
					}
					if (ext_stop_sim) break;
				}

				wave_function.set_m2psi_coh(output_multem);
			}

			template <class TOutput_multislice>
			void EWFS_EWRS(TOutput_multislice &output_multem)
			{
				ext_niter = in_multem->number_pn_conf();
				ext_iter = 0;
				/***************************************************************************************/

				T_r w = in_multem->get_phonon_rot_weight();

				output_multem.init();

				for(auto iconf = in_multem->phonon_par.iconf_0; iconf <= in_multem->phonon_par.nconf; iconf++)
				{
					wave_function.move_atoms(iconf);

					wave_function.set_incident_wave(wave_function.psi_z);

					wave_function.psi(w, wave_function.psi_z, output_multem);

					ext_iter++;
					if (ext_stop_sim) break;
				}

				wave_function.set_m2psi_coh(output_multem);
			}

			template <class TOutput_multislice>
			void EELS_EFTEM(TOutput_multislice &output_multem)
			{
				ext_niter = wave_function.slicing.slice.size()*in_multem->number_pn_conf();
				ext_iter = 0;

				if (in_multem->is_STEM_ISTEM_EELS())
				{
					ext_niter *= in_multem->scanning.size();
				}
				/***************************************************************************************/

				T_r w = in_multem->get_phonon_rot_weight();

				auto psi = [&](T_r w, Vctr<T_c, Dev>& psi_z, TOutput_multislice &output_multem)
				{
					T_r gx_0 = in_multem->gx_0();
					T_r gy_0 = in_multem->gy_0();

					for(auto islice = 0; islice < wave_function.slicing.slice.size(); islice++)
					{
						if (in_multem->eels_fr.is_Mixed_Chan())
						{
							wave_function.trans(islice, wave_function.slicing.slice.size()-1, trans_thk);
						}			

						for(auto iatoms = wave_function.slicing.slice[islice].iatom_0; iatoms <= wave_function.slicing.slice[islice].iatom_e; iatoms++)
						{
							if (wave_function.atoms.Z[iatoms] == in_multem->eels_fr.Z)
							{
								in_multem->set_eels_fr_atom(iatoms, wave_function.atoms);
								energy_loss.set_atom_type(in_multem->eels_fr);

								for(auto ikn = 0; ikn < energy_loss.kernel.size(); ikn++)
								{
									mt::ew_mult(stream, energy_loss.kernel[ikn], psi_z, wave_function.psi_z);
									wave_function.psi(islice, wave_function.slicing.slice.size()-1, w, trans_thk, output_multem);
								}
							}

							if (ext_stop_sim) break;
						}
						wave_function.psi_slice(gx_0, gy_0, islice, psi_z);

						ext_iter++;
						if (ext_stop_sim) break;
					}
				};

				output_multem.init();

				in_multem->ibeam.resize(1);
				in_multem->beam_x.resize(1);
				in_multem->beam_y.resize(1);

				if (in_multem->is_STEM_ISTEM_EELS())
				{
					for(auto iconf = in_multem->phonon_par.iconf_0; iconf <= in_multem->phonon_par.nconf; iconf++)
					{
						wave_function.move_atoms(iconf);
						for(auto ibeam = 0; ibeam < in_multem->scanning.size(); ibeam++)
						{
							in_multem->ibeam[0] = ibeam;
							in_multem->set_iscan_beam_position();
							wave_function.set_incident_wave(psi_thk);
							psi(w, psi_thk, output_multem);

							if (ext_stop_sim) break;
						}

						if (ext_stop_sim) break;
					}
				}
				else
				{
					for(auto iconf = in_multem->phonon_par.iconf_0; iconf <= in_multem->phonon_par.nconf; iconf++)
					{
						wave_function.move_atoms(iconf);
						wave_function.set_incident_wave(psi_thk);
						psi(w, psi_thk, output_multem);

						if (ext_stop_sim) break;
					}
				}
			}

			In_Multem<T_r> *in_multem;
			Stream<Dev> stream;
			FFT<T_r, Dev> fft_2d;

			Wave_Function<T_r, Dev> wave_function;
			Energy_Loss<T_r, Dev> energy_loss;

			Vctr<T_c, Dev> psi_thk;
			Vctr<T_c, Dev> trans_thk;
	};

	template <class T, eDev Dev>
	dt_bool Tem_Simulation<T, Dev>::ext_stop_sim = false;

	template <class T, eDev Dev>
	dt_int32 Tem_Simulation<T, Dev>::ext_niter = 0;

	template <class T, eDev Dev>
	dt_int32 Tem_Simulation<T, Dev>::ext_iter = 0;

	template <class T, eDev Dev>
	class Multem
	{
		public:
			Multem(): n_devices(1), in_multem(nullptr) {}

			Multem(In_Multem<T> *in_multem_i)
			{
				set_in_data(in_multem_i);
			}

			void set_in_data(In_Multem<T> *in_multem_i)
			{
				in_multem = in_multem_i;
				n_devices = 1;
				if (in_multem->is_STEM_ISTEM()||in_multem->is_CBED_CBEI())
				{
					n_devices = in_multem->system_config.get_n_sel_gpu();
				}
				output_multem_v.resize(n_devices);
			}

			template <class TOutput_Multem>
			void operator()(TOutput_Multem &output_multem)
			{
				if (in_multem->is_STEM_ISTEM())
				{
					STEM_ISTEM(output_multem);
				}
				else if (in_multem->is_CBED_CBEI())
				{
					CBED_CBEI(output_multem);
				}
				else if (in_multem->is_ED_HRTEM())
				{
					ED_HRTEM(output_multem);
				}
				else if (in_multem->is_PED_HCTEM())
				{
					PED_HCTEM(output_multem);
				}
				else if (in_multem->is_EWFS_EWRS())
				{
					EWFS_EWRS(output_multem);
				}
				else if (in_multem->is_EELS_EFTEM())
				{
					EELS_EFTEM(output_multem);
				}
			}
		private:
			template <class TOutput_Multem>
			void STEM_ISTEM(TOutput_Multem &output_multem)
			{
 				vector<std::thread> threads;
				threads.reserve(n_devices);

				auto stem_istem_thr =[&](dt_int32 ithr)
				{
					In_Multem<T> in_multem_thr = *in_multem;
					in_multem_thr.scanning.type = eST_user_def;
					in_multem_thr.scanning.R = in_multem->extract_beam_pos(ithr, n_devices);
					in_multem_thr.validate_parameters();

 					in_multem_thr.system_config.set_gpu_by_ind(ithr);

					Tem_Simulation<T, Dev> tem_simulation(&in_multem_thr);
					output_multem_v[ithr].set_in_data(&in_multem_thr);

					tem_simulation(output_multem_v[ithr]);
					tem_simulation.cleanup();
				};

				for(auto ithr=0; ithr<n_devices; ithr++)
				{
					threads.emplace_back(stem_istem_thr, ithr);
				}

				for(auto ithr=0; ithr<n_devices; ithr++)
				{
					if (threads[ithr].joinable())
					{
						threads[ithr].join();
					}
				}

				output_multem.joint_data(output_multem_v);
			}

			template <class TOutput_Multem>
			void CBED_CBEI(TOutput_Multem &output_multem)
			{
 				vector<std::thread> threads;
				threads.reserve(n_devices);

				auto cbed_cbei_thr =[&](dt_int32 ithr)
				{
					In_Multem<T> in_multem_thr = *in_multem;
					in_multem_thr.beam_x = in_multem->extract_probe_pos_x(ithr, n_devices);
					in_multem_thr.beam_y = in_multem->extract_probe_pos_y(ithr, n_devices);
					in_multem_thr.validate_parameters();

 					in_multem_thr.system_config.set_gpu_by_ind(ithr);

					Tem_Simulation<T, Dev> tem_simulation;
					tem_simulation.set_in_data(&in_multem_thr);
					output_multem_v[ithr].set_in_data(&in_multem_thr);

					tem_simulation(output_multem_v[ithr]);
					tem_simulation.cleanup();
				};

				for(auto ithr=0; ithr<n_devices; ithr++)
				{
					threads.emplace_back(cbed_cbei_thr, ithr);
				}

				for(auto ithr=0; ithr<n_devices; ithr++)
				{
					if (threads[ithr].joinable())
					{
						threads[ithr].join();
					}
				}

				output_multem.joint_data(output_multem_v);
			}

			template <class TOutput_Multem>
			void ED_HRTEM(TOutput_Multem &output_multem)
			{
				In_Multem<T> in_multem_thr = *in_multem;
				in_multem_thr.system_config.set_gpu_by_ind(0);

				Tem_Simulation<T, Dev> tem_simulation;
				tem_simulation.set_in_data(&in_multem_thr);
				output_multem.set_in_data(&in_multem_thr);

				tem_simulation(output_multem);
				tem_simulation.cleanup();
			}

			template <class TOutput_Multem>
			void PED_HCTEM(TOutput_Multem &output_multem)
			{
				In_Multem<T> in_multem_thr = *in_multem;
				in_multem_thr.system_config.set_gpu_by_ind(0);

				Tem_Simulation<T, Dev> tem_simulation;
				tem_simulation.set_in_data(&in_multem_thr);
				output_multem.set_in_data(&in_multem_thr);

				tem_simulation(output_multem);
				tem_simulation.cleanup();
			}

			template <class TOutput_Multem>
			void EWFS_EWRS(TOutput_Multem &output_multem)
			{
				In_Multem<T> in_multem_thr = *in_multem;
				in_multem_thr.system_config.set_gpu_by_ind(0);

				Tem_Simulation<T, Dev> tem_simulation;
				tem_simulation.set_in_data(&in_multem_thr);
				output_multem.set_in_data(&in_multem_thr);

				tem_simulation(output_multem);
				tem_simulation.cleanup();
			}

			template <class TOutput_Multem>
			void EELS_EFTEM(TOutput_Multem &output_multem)
			{
				In_Multem<T> in_multem_thr = *in_multem;
				in_multem_thr.system_config.set_gpu_by_ind(0);

				Tem_Simulation<T, Dev> tem_simulation;
				tem_simulation.set_in_data(&in_multem_thr);
				output_multem.set_in_data(&in_multem_thr);

				tem_simulation(output_multem);
				tem_simulation.cleanup();
			}

			dt_int32 n_devices;
			In_Multem<T> *in_multem;
			vector<Output_Multem<T>> output_multem_v;
	};

}

#endif
