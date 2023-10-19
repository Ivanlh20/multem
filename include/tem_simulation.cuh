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
	#include "math_mt.h"
	#include "types.cuh"
	#include "type_traits_gen.h"
	#include "in_classes.cuh"
	#include "output_multem.hpp"
	#include "fcns_cpu.h"
	#include "fcns_gpu.h"
	#include "fcns_gpu.h"
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

				Tem_Simulation(): multem_in_parm(nullptr) {}

				Tem_Simulation(Multem_In_Parm<T> *multem_in_parm_i)
				{
					set_in_data(multem_in_parm_i);
				}

				void set_in_data(Multem_In_Parm<T> *multem_in_parm_i)
				{
					multem_in_parm = multem_in_parm_i;

					stream.resize(multem_in_parm->system_config.n_stream);

					fft_2d.create_plan_2d(multem_in_parm->grid_2d.ny, multem_in_parm->grid_2d.nx, multem_in_parm->system_config.n_stream);

					if (multem_in_parm->is_EELS_EFTEM())
					{
						energy_loss.set_in_data(multem_in_parm, &stream, &fft_2d);
						psi_thk.resize(multem_in_parm->grid_2d.size());
						if (multem_in_parm->eels_fr.is_Mixed_Chan())
						{
							trans_thk.resize(multem_in_parm->grid_2d.size());
						}
					}

					wave_function.set_in_data(multem_in_parm, &stream, &fft_2d);
				}

				template <class TOutput_multislice>
				void operator()(TOutput_multislice &output_multem)
				{
					if (multem_in_parm->is_STEM_ISTEM())
					{
						STEM_ISTEM(output_multem);
					}
					else if (multem_in_parm->is_CBED_CBEI())
					{
						CBED_CBEI(output_multem);
					}
					else if (multem_in_parm->is_ED_HRTEM())
					{
						ED_HRTEM(output_multem);
					}
					else if (multem_in_parm->is_PED_HCTEM())
					{
						PED_HCTEM(output_multem);
					}
					else if (multem_in_parm->is_EWFS_EWRS())
					{
						EWFS_EWRS(output_multem);
					}
					else if (multem_in_parm->is_EELS_EFTEM())
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
					ext_niter = multem_in_parm->scanning.size()*multem_in_parm->number_pn_conf();
					ext_iter = 0;
					/***************************************************************************************/

					T_r w = multem_in_parm->get_phonon_rot_weight();

					output_multem.init();

					multem_in_parm->ibeam.resize(1);
					multem_in_parm->beam_x.resize(1);
					multem_in_parm->beam_y.resize(1);

				
					if (multem_in_parm->is_STEM() && multem_in_parm->atomic_vib.coh_contrib)
					{
						for(auto ibeam = 0; ibeam < multem_in_parm->scanning.size(); ibeam++)
						{	
							output_multem.init_psi_coh();
							multem_in_parm->ibeam[0] = ibeam;
							multem_in_parm->set_iscan_beam_position();
							for(auto iconf = multem_in_parm->atomic_vib.iconf_0; iconf <= multem_in_parm->atomic_vib.nconf; iconf++)
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
						for(auto iconf = multem_in_parm->atomic_vib.iconf_0; iconf <= multem_in_parm->atomic_vib.nconf; iconf++)
						{
							wave_function.move_atoms(iconf);
							for(auto ibeam = 0; ibeam < multem_in_parm->scanning.size(); ibeam++)
							{
								multem_in_parm->ibeam[0] = ibeam;
								multem_in_parm->set_iscan_beam_position();
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

					multem_in_parm->ibeam.resize(1);
					multem_in_parm->beam_x.resize(1);
					multem_in_parm->beam_y.resize(1);

					Quad_Coef_1d<dt_float64, edev_cpu> qt;
					Quad_Coef_2d<dt_float64, edev_cpu> qs;

					// Load quadratures
					cond_lens_temporal_spatial_quadratures(multem_in_parm->cond_lens, qt, qs);

					/***************************************************************************************/
					dt_float64 w_pr_0 = multem_in_parm->get_phonon_rot_weight();
					dt_float64 c_10_0 = multem_in_parm->cond_lens.c_10;
					const dt_int32 n_beams = multem_in_parm->number_of_beams();

					ext_niter = qs.size()*qt.size()*multem_in_parm->number_pn_conf();
					ext_iter = 0;

					for(auto iconf = multem_in_parm->atomic_vib.iconf_0; iconf <= multem_in_parm->atomic_vib.nconf; iconf++)
					{
						wave_function.move_atoms(iconf);

						for(auto ibeam=0; ibeam<n_beams; ibeam++)
						{
							multem_in_parm->ibeam[0] = ibeam;
							multem_in_parm->beam_x[0] = multem_in_parm->beam_x[ibeam];
							multem_in_parm->beam_y[0] = multem_in_parm->beam_y[ibeam];

							// spatial incoherence
							for(auto ispat = 0; ispat<qs.size(); ispat++)
							{
								auto beam_x = qs.x[ispat] + multem_in_parm->beam_x[0];
								auto beam_y = qs.y[ispat] + multem_in_parm->beam_y[0];

								// temporal incoherence
								for(auto itemp = 0; itemp<qt.size(); itemp++)
								{
									auto c_10 = c_10_0 + qt.x[itemp];
									auto w = w_pr_0*qs.w[ispat]*qt.w[itemp];

									multem_in_parm->cond_lens.set_defocus(c_10);
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

					multem_in_parm->cond_lens.set_defocus(c_10_0);
					multem_in_parm->set_beam_position(multem_in_parm->beam_x, multem_in_parm->beam_y);
				}

				template <class TOutput_multislice>
				void ED_HRTEM(TOutput_multislice &output_multem)
				{
					EWFS_EWRS(output_multem);
				}

				template <class TOutput_multislice>
				void PED_HCTEM(TOutput_multislice &output_multem)
				{
					ext_niter = multem_in_parm->nrot*multem_in_parm->number_pn_conf();
					ext_iter = 0;
					/***************************************************************************************/

					T_r w = multem_in_parm->get_phonon_rot_weight();

					output_multem.init();

					for(auto iconf = multem_in_parm->atomic_vib.iconf_0; iconf <= multem_in_parm->atomic_vib.nconf; iconf++)
					{
						wave_function.move_atoms(iconf);
						for(auto irot = 0; irot < multem_in_parm->nrot; irot++)
						{
							multem_in_parm->set_phi(irot);
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
					ext_niter = multem_in_parm->number_pn_conf();
					ext_iter = 0;
					/***************************************************************************************/

					T_r w = multem_in_parm->get_phonon_rot_weight();

					output_multem.init();

					for(auto iconf = multem_in_parm->atomic_vib.iconf_0; iconf <= multem_in_parm->atomic_vib.nconf; iconf++)
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
					ext_niter = wave_function.slicing.slice.size()*multem_in_parm->number_pn_conf();
					ext_iter = 0;

					if (multem_in_parm->is_STEM_ISTEM_EELS())
					{
						ext_niter *= multem_in_parm->scanning.size();
					}
					/***************************************************************************************/

					T_r w = multem_in_parm->get_phonon_rot_weight();

					auto psi = [&](T_r w, Vctr<T_c, Dev>& psi_z, TOutput_multislice &output_multem)
					{
						T_r gx_0 = multem_in_parm->gx_0();
						T_r gy_0 = multem_in_parm->gy_0();

						for(auto islice = 0; islice < wave_function.slicing.slice.size(); islice++)
						{
							if (multem_in_parm->eels_fr.is_Mixed_Chan())
							{
								wave_function.trans(islice, wave_function.slicing.slice.size()-1, trans_thk);
							}			

							for(auto iatoms = wave_function.slicing.slice[islice].iatom_0; iatoms <= wave_function.slicing.slice[islice].iatom_e; iatoms++)
							{
								if (wave_function.atoms.Z[iatoms] == multem_in_parm->eels_fr.Z)
								{
									multem_in_parm->set_eels_fr_atom(iatoms, wave_function.atoms);
									energy_loss.set_atom_type(multem_in_parm->eels_fr);

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

					multem_in_parm->ibeam.resize(1);
					multem_in_parm->beam_x.resize(1);
					multem_in_parm->beam_y.resize(1);

					if (multem_in_parm->is_STEM_ISTEM_EELS())
					{
						for(auto iconf = multem_in_parm->atomic_vib.iconf_0; iconf <= multem_in_parm->atomic_vib.nconf; iconf++)
						{
							wave_function.move_atoms(iconf);
							for(auto ibeam = 0; ibeam < multem_in_parm->scanning.size(); ibeam++)
							{
								multem_in_parm->ibeam[0] = ibeam;
								multem_in_parm->set_iscan_beam_position();
								wave_function.set_incident_wave(psi_thk);
								psi(w, psi_thk, output_multem);

								if (ext_stop_sim) break;
							}

							if (ext_stop_sim) break;
						}
					}
					else
					{
						for(auto iconf = multem_in_parm->atomic_vib.iconf_0; iconf <= multem_in_parm->atomic_vib.nconf; iconf++)
						{
							wave_function.move_atoms(iconf);
							wave_function.set_incident_wave(psi_thk);
							psi(w, psi_thk, output_multem);

							if (ext_stop_sim) break;
						}
					}
				}

				Multem_In_Parm<T_r> *multem_in_parm;
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
				Multem(): n_devices(1), multem_in_parm(nullptr) {}

				Multem(Multem_In_Parm<T> *multem_in_parm_i)
				{
					set_in_data(multem_in_parm_i);
				}

				void set_in_data(Multem_In_Parm<T> *multem_in_parm_i)
				{
					multem_in_parm = multem_in_parm_i;
					n_devices = 1;
					if (multem_in_parm->is_STEM_ISTEM()||multem_in_parm->is_CBED_CBEI())
					{
						n_devices = multem_in_parm->system_config.get_n_sel_gpu();
					}
					output_multem_v.resize(n_devices);
				}

				template <class TOutput_Multem>
				void operator()(TOutput_Multem &output_multem)
				{
					if (multem_in_parm->is_STEM_ISTEM())
					{
						STEM_ISTEM(output_multem);
					}
					else if (multem_in_parm->is_CBED_CBEI())
					{
						CBED_CBEI(output_multem);
					}
					else if (multem_in_parm->is_ED_HRTEM())
					{
						ED_HRTEM(output_multem);
					}
					else if (multem_in_parm->is_PED_HCTEM())
					{
						PED_HCTEM(output_multem);
					}
					else if (multem_in_parm->is_EWFS_EWRS())
					{
						EWFS_EWRS(output_multem);
					}
					else if (multem_in_parm->is_EELS_EFTEM())
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
						Multem_In_Parm<T> multem_in_parm_thr = *multem_in_parm;
						multem_in_parm_thr.scanning.type = espt_user_def;
						multem_in_parm_thr.scanning.R = multem_in_parm->extract_beam_pos(ithr, n_devices);
						multem_in_parm_thr.set_dep_var();

 						multem_in_parm_thr.system_config.set_gpu_by_ind(ithr);

						Tem_Simulation<T, Dev> tem_simulation(&multem_in_parm_thr);
						output_multem_v[ithr].set_in_data(&multem_in_parm_thr);

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
						Multem_In_Parm<T> multem_in_parm_thr = *multem_in_parm;
						multem_in_parm_thr.beam_x = multem_in_parm->extract_probe_pos_x(ithr, n_devices);
						multem_in_parm_thr.beam_y = multem_in_parm->extract_probe_pos_y(ithr, n_devices);
						multem_in_parm_thr.set_dep_var();

 						multem_in_parm_thr.system_config.set_gpu_by_ind(ithr);

						Tem_Simulation<T, Dev> tem_simulation;
						tem_simulation.set_in_data(&multem_in_parm_thr);
						output_multem_v[ithr].set_in_data(&multem_in_parm_thr);

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
					Multem_In_Parm<T> multem_in_parm_thr = *multem_in_parm;
					multem_in_parm_thr.system_config.set_gpu_by_ind(0);

					Tem_Simulation<T, Dev> tem_simulation;
					tem_simulation.set_in_data(&multem_in_parm_thr);
					output_multem.set_in_data(&multem_in_parm_thr);

					tem_simulation(output_multem);
					tem_simulation.cleanup();
				}

				template <class TOutput_Multem>
				void PED_HCTEM(TOutput_Multem &output_multem)
				{
					Multem_In_Parm<T> multem_in_parm_thr = *multem_in_parm;
					multem_in_parm_thr.system_config.set_gpu_by_ind(0);

					Tem_Simulation<T, Dev> tem_simulation;
					tem_simulation.set_in_data(&multem_in_parm_thr);
					output_multem.set_in_data(&multem_in_parm_thr);

					tem_simulation(output_multem);
					tem_simulation.cleanup();
				}

				template <class TOutput_Multem>
				void EWFS_EWRS(TOutput_Multem &output_multem)
				{
					Multem_In_Parm<T> multem_in_parm_thr = *multem_in_parm;
					multem_in_parm_thr.system_config.set_gpu_by_ind(0);

					Tem_Simulation<T, Dev> tem_simulation;
					tem_simulation.set_in_data(&multem_in_parm_thr);
					output_multem.set_in_data(&multem_in_parm_thr);

					tem_simulation(output_multem);
					tem_simulation.cleanup();
				}

				template <class TOutput_Multem>
				void EELS_EFTEM(TOutput_Multem &output_multem)
				{
					Multem_In_Parm<T> multem_in_parm_thr = *multem_in_parm;
					multem_in_parm_thr.system_config.set_gpu_by_ind(0);

					Tem_Simulation<T, Dev> tem_simulation;
					tem_simulation.set_in_data(&multem_in_parm_thr);
					output_multem.set_in_data(&multem_in_parm_thr);

					tem_simulation(output_multem);
					tem_simulation.cleanup();
				}

				dt_int32 n_devices;
				Multem_In_Parm<T> *multem_in_parm;
				vector<Output_Multem<T>> output_multem_v;
		};

	}

#endif
