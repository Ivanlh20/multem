/*
 * This file is part of MULTEM.
 * Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
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
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MULTISLICE_H
#define MULTISLICE_H

#include "fftw3.h"
#include "types.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "input_multislice.hpp"
#include "microscope_effects.cuh"
#include "energy_loss.cuh"
#include "wave_function.cuh"

namespace multem
{
	template<class T, eDevice dev>
	class Multislice
	{
		public:
			using value_type_r = typename T;
			using value_type_c = typename complex<T>;

			void set_input_data(Input_Multislice<value_type_r, dev> *input_multislice_io)
			{
				input_multislice = input_multislice_io;
				if(dev == e_Device)
				{
					cudaSetDevice(input_multislice->gpu_device);
				}
				stream.resize(input_multislice->nstream);
				fft2.create_plan(input_multislice->grid.ny, input_multislice->grid.nx, input_multislice->nstream);

				psi_coh.resize(input_multislice->grid.nxy());
				m2psi_coh.resize(input_multislice->grid.nxy());
				m2psi_tot.resize(input_multislice->grid.nxy());
				
				if(input_multislice->is_EELS_EFTEM())
				{
					energy_loss.set_input_data(input_multislice, &fft2);
					psi_thk.resize(input_multislice->grid.nxy());
					psi_kn.resize(input_multislice->grid.nxy());
				}

				if(input_multislice->is_HRTEM() || input_multislice->is_HCI())
				{
					microscope_effects.set_input_data(input_multislice, &stream, &fft2);
				}

				wave_function.set_input_data(input_multislice, &stream, &fft2);
			}

			template<class TVector_Host_r>
			void STEM(TVector_Host_r &det_int_tot, TVector_Host_r &det_int_coh)
			{
				value_type_r w = input_multislice->ifp_nconf();

				for(auto is=0; is < input_multislice->scanning.size(); is++)
				{
					input_multislice->set_stem_beam_position(is);

					reset_m2psi_tot_and_psi_coh(m2psi_tot, psi_coh);

					for(auto iconf=1; iconf <= input_multislice->fp_nconf; iconf++)
					{
						wave_function.move_atoms(iconf);
						get_total_intensity_and_coherent_wave(w, m2psi_tot, psi_coh);
					}
					get_coherent_intensity(psi_coh, m2psi_coh);

					for(auto iDet=0; iDet<input_multislice->det_cir.size(); iDet++)
					{
						value_type_r g_inner = input_multislice->det_cir.g_inner(iDet);
						value_type_r g_outer = input_multislice->det_cir.g_outer(iDet);
						det_int_tot.image[iDet][is] = w*sum_over_Det(input_multislice->grid, g_inner, g_outer, m2psi_tot);
						det_int_coh.image[iDet][is] = w*sum_over_Det(input_multislice->grid, g_inner, g_outer, m2psi_coh);
					}
				}
			}

			template<class TVector_Host_r>
			void ISTEM(TVector_Host_r &host_m2psi_tot, TVector_Host_r &host_m2psi_coh)
			{
				value_type_r w = input_multislice->ifp_nconf();

				reset_m2psi_tot_and_psi_coh(m2psi_tot, psi_coh);

				for(auto iconf=1; iconf <= input_multislice->fp_nconf; iconf++)
				{
					wave_function.move_atoms(iconf);
					for(auto is=0; is < input_multislice->scanning.size(); is++)
					{
						input_multislice->set_stem_beam_position(is);
						get_total_intensity_and_coherent_wave(w, m2psi_tot, psi_coh);
					}
				}
				get_coherent_intensity(psi_coh, m2psi_coh);

				multem::to_host_shift(input_multislice->grid, m2psi_tot, host_m2psi_tot);
				if(input_multislice->coherent_contribution)
				{
					multem::to_host_shift(input_multislice->grid, m2psi_coh, host_m2psi_coh);
				}
			}

			template<class TVector_Host_r>
			void CBED_CBEI(TVector_Host_r &host_m2psi_tot, TVector_Host_r &host_m2psi_coh)
			{
				tot_coh_FS_RS(m2psi_tot, psi_coh);
				get_coherent_intensity(psi_coh, m2psi_coh);

				multem::to_host_shift(input_multislice->grid, m2psi_tot, host_m2psi_tot);
				if(input_multislice->coherent_contribution)
				{
					multem::to_host_shift(input_multislice->grid, m2psi_coh, host_m2psi_coh);
				}
			}

			template<class TVector_Host_r>
			void ED_HRTEM(TVector_Host_r &host_m2psi_tot, TVector_Host_r &host_m2psi_coh)
			{
				tot_coh_FS_RS(m2psi_tot, psi_coh);
				get_coherent_intensity(psi_coh, m2psi_coh);

				multem::to_host_shift(input_multislice->grid, m2psi_tot, host_m2psi_tot);
				if(input_multislice->coherent_contribution)
				{
					multem::to_host_shift(input_multislice->grid, m2psi_coh, host_m2psi_coh);
				}
			}

			template<class TVector_Host_r>
			void PED_HCI(TVector_Host_r &host_m2psi_tot, TVector_Host_r &host_m2psi_coh)
			{
				input_multislice->theta = input_multislice->pe_fr.theta;
				value_type_r w = input_multislice->pe_fr.weight(input_multislice->fp_nconf);

				reset_m2psi_tot_and_psi_coh(m2psi_tot, psi_coh);

				for(auto iconf=1; iconf <= input_multislice->fp_nconf; iconf++)
				{
					wave_function.move_atoms(iconf);
					for(auto irot=0; irot < input_multislice->pe_fr.nrot; irot++)
					{
						input_multislice->phi = input_multislice->pe_fr.phi(irot);
						get_total_intensity_and_coherent_wave(w, m2psi_tot, psi_coh);
					}
				}
				get_coherent_intensity(psi_coh, m2psi_coh);
				
				multem::to_host_shift(input_multislice->grid, m2psi_tot, host_m2psi_tot);
				if(input_multislice->coherent_contribution)
				{
					multem::to_host_shift(input_multislice->grid, m2psi_coh, host_m2psi_coh);
				}
			}

			template<class TVector_Host_r, class TVector_Host_c>
			void EWFS_EWRS(TVector_Host_r &host_m2psi_tot, TVector_Host_c &host_psi_coh)
			{
				tot_coh_FS_RS(m2psi_tot, psi_coh);

				multem::to_host_shift(input_multislice->grid, m2psi_tot, host_m2psi_tot);
				if(input_multislice->coherent_contribution)
				{
					multem::to_host_shift(input_multislice->grid, psi_coh, host_psi_coh);
				}
			}

			template<class TVector_Host_r>
			void STEM(TVector_Host_r &det_int_tot)
			{
				value_type_r w = input_multislice->ifp_nconf();

				for(auto iDet=0; iDet<input_multislice->det_cir.size(); iDet++)
				{
					multem::fill(det_int_tot.image[iDet], 0.0);
				}

				for(auto iconf=1; iconf <= input_multislice->fp_nconf; iconf++)
				{
					wave_function.move_atoms(iconf);
					for(auto is=0; is < input_multislice->scanning.size(); is++)
					{
						input_multislice->set_stem_beam_position(is);

						wave_function.psi_0();
						wave_function.psi(eS_Reciprocal);

						for(auto iDet=0; iDet<input_multislice->det_cir.size(); iDet++)
						{
							value_type_r g_inner = input_multislice->det_cir.g_inner(iDet);
							value_type_r g_outer = input_multislice->det_cir.g_outer(iDet);
							det_int_tot.image[iDet][is] += w*sum_square_over_Det(input_multislice->grid, g_inner, g_outer, wave_function.psi_z);
						}
					}
				}
			}

			template<class TVector_Host_r>
			void EFTEM(TVector_Host_r &host_m2psi_tot)
			{
				value_type_r w = input_multislice->ifp_nconf();

				multem::fill(m2psi_tot, 0.0);

				for(auto iconf=1; iconf <= input_multislice->fp_nconf; iconf++)
				{
					wave_function.move_atoms(iconf);
					for(auto ithk=0; ithk < wave_function.slice.size(); ithk++)
					{
						wave_function.psi_0();
						wave_function.psi(eS_Real, ithk);
						multem::assign(wave_function.psi_z, psi_thk);

						for(auto iatom=wave_function.slice.iatom_0[ithk]; iatom <= wave_function.slice.iatom_e[ithk]; iatom++)
						{
							if(wave_function.atoms.Z[iatom] == input_multislice->eels_fr.Z)
							{
								input_multislice->eels_fr.x = -c_2Pi*(wave_function.atoms.x[iatom]-input_multislice->grid.lxh());
								input_multislice->eels_fr.y = -c_2Pi*(wave_function.atoms.y[iatom]-input_multislice->grid.lyh());
								input_multislice->eels_fr.occ = wave_function.atoms.occ[iatom];

								energy_loss.set_atom_type(input_multislice->eels_fr);
								for(auto ikn=0; ikn < energy_loss.kernel.size(); ikn++)
								{
									multem::multiply(energy_loss.kernel[ikn], psi_thk, psi_kn);
									wave_function.psi(eS_Real, ithk, wave_function.slice.size(), &psi_kn);
									multem::add_square(wave_function.psi_z, m2psi_tot);
								}
							}
						}
					}
				}

				multem::to_host_shift(input_multislice->grid, m2psi_tot, host_m2psi_tot);
			}

			template<class TVector_Host_r>
			void EELS(TVector_Host_r &host_stem_eels)
			{
				value_type_r w = input_multislice->ifp_nconf();

				multem::fill(host_stem_eels, 0.0);

				for(auto iconf=1; iconf <= input_multislice->fp_nconf; iconf++)
				{
					wave_function.move_atoms(iconf);
					for(auto is=0; is < input_multislice->scanning.size(); is++)
					{
						input_multislice->set_stem_beam_position(is);
						for(auto ithk=0; ithk < wave_function.slice.size(); ithk++)
						{
							wave_function.psi_0();
							wave_function.psi(eS_Real, ithk);
							multem::assign(wave_function.psi_z, psi_thk);

							for(auto iatom=wave_function.slice.iatom_0[ithk]; iatom <= wave_function.slice.iatom_e[ithk]; iatom++)
							{
								if(wave_function.atoms.Z[iatom] == input_multislice->eels_fr.Z)
								{
									input_multislice->eels_fr.x = -c_2Pi*(wave_function.atoms.x[iatom]-input_multislice->grid.lxh());
									input_multislice->eels_fr.y = -c_2Pi*(wave_function.atoms.y[iatom]-input_multislice->grid.lyh());
									input_multislice->eels_fr.occ = wave_function.atoms.occ[iatom];

									energy_loss.set_atom_type(input_multislice->eels_fr);
									for(auto ikn=0; ikn < energy_loss.kernel.size(); ikn++)
									{
										multem::multiply(energy_loss.kernel[ikn], psi_thk, psi_kn);
										wave_function.psi(eS_Reciprocal, ithk, wave_function.slice.size(), &psi_kn);
										host_stem_eels[is] += sum_square_over_Det(input_multislice->grid, 0, input_multislice->eels_fr.g_collection, wave_function.psi_z);
									}
								}
							}
						}
					}
				}
			}

			template<class TMatlab_1>
			void output_matlab(TMatlab_1 *host_tot)
			{
				using value_type_1 = multem::traits::Value_type<TMatlab_1>;

				if(input_multislice->is_STEM())
				{
					auto *tot = reinterpret_cast<Det_Int<value_type_1, e_Host>*>(host_tot);

					STEM(*tot);
				}
				else if(input_multislice->is_ISTEM())
				{
					auto *tot = reinterpret_cast<m_matrix_r*>(host_tot);
					m_matrix_r coh;

					ISTEM(*tot, coh);
				}
				else if(input_multislice->is_CBED_CBEI())
				{
					auto *tot = reinterpret_cast<m_matrix_r*>(host_tot);
					m_matrix_r coh;

					CBED_CBEI(*tot, coh);
				}
				else if(input_multislice->is_ED_HRTEM())
				{
					auto *tot = reinterpret_cast<m_matrix_r*>(host_tot);
					m_matrix_r coh;

					ED_HRTEM(*tot, coh);
				}
				else if(input_multislice->is_PED_HCI())
				{
					auto *tot = reinterpret_cast<m_matrix_r*>(host_tot);
					m_matrix_r coh;

					PED_HCI(*tot, coh);
				}
				else if(input_multislice->is_EWFS_EWRS())
				{
					auto *tot = reinterpret_cast<m_matrix_r*>(host_tot);
					m_matrix_c coh;

					EWFS_EWRS(*tot, coh);
				}
				else if(input_multislice->is_EFTEM())
				{
					auto *tot = reinterpret_cast<m_matrix_r*>(host_tot);

					EFTEM(*tot);
				}
				else if(input_multislice->is_EELS())
				{
					auto *tot = reinterpret_cast<Vector<value_type_1, e_Host>*>(host_tot);

					EELS(*tot);
				}
			}

			template<class TMatlab_1, class TMatlab_2>
			void output_matlab(TMatlab_1 *host_tot, TMatlab_2 *host_coh)
			{
				using value_type_1 = multem::traits::Value_type<TMatlab_1>;
				using value_type_2 = multem::traits::Value_type<TMatlab_2>;

				if(input_multislice->is_STEM())
				{
					auto *tot = reinterpret_cast<Det_Int<value_type_1, e_Host>*>(host_tot);
					auto *coh = reinterpret_cast<Det_Int<value_type_2, e_Host>*>(host_coh);

					STEM(*tot, *coh);
				}
				else if(input_multislice->is_ISTEM())
				{
					auto *tot = reinterpret_cast<m_matrix_r*>(host_tot);
					auto *coh = reinterpret_cast<m_matrix_r*>(host_coh);

					ISTEM(*tot, *coh);
				}
				else if(input_multislice->is_CBED_CBEI())
				{
					auto *tot = reinterpret_cast<m_matrix_r*>(host_tot);
					auto *coh = reinterpret_cast<m_matrix_r*>(host_coh);

					CBED_CBEI(*tot, *coh);
				}
				else if(input_multislice->is_ED_HRTEM())
				{
					auto *tot = reinterpret_cast<m_matrix_r*>(host_tot);
					auto *coh = reinterpret_cast<m_matrix_r*>(host_coh);

					ED_HRTEM(*tot, *coh);
				}
				else if(input_multislice->is_PED_HCI())
				{
					auto *tot = reinterpret_cast<m_matrix_r*>(host_tot);
					auto *coh = reinterpret_cast<m_matrix_r*>(host_coh);

					PED_HCI(*tot, *coh);
				}
				if(input_multislice->is_EWFS_EWRS())
				{
					auto *tot = reinterpret_cast<m_matrix_r*>(host_tot);
					auto *coh = reinterpret_cast<m_matrix_c*>(host_coh);

					EWFS_EWRS(*tot, *coh);
				}
			}

			void cleanup()
			{
				fft2.cleanup();
			}

		private:
			void reset_m2psi_tot_and_psi_coh(Vector<value_type_r, dev> &m2psi_tot, Vector<value_type_c, dev> &psi_coh)
			{
				multem::fill(m2psi_tot, 0.0);
				if(input_multislice->coherent_contribution)
				{
					multem::fill(psi_coh, 0.0);
				}
			}

			void get_total_intensity_and_coherent_wave(const value_type_r &w, Vector<value_type_r, dev> &m2psi_tot, Vector<value_type_c, dev> &psi_coh)
			{
				wave_function.psi_0();
							
				if(input_multislice->is_simulation_type_FS())
				{
					wave_function.psi(eS_Reciprocal);
					multem::add_square_scale(w, wave_function.psi_z, m2psi_tot);
				}
				else if(input_multislice->is_HRTEM_HCI_EFTEM())
				{
					wave_function.psi(eS_Reciprocal);
					microscope_effects.apply(wave_function.psi_z, m2psi_coh);
					multem::add_scale(w, m2psi_coh, m2psi_tot);
				}
				else
				{
					wave_function.psi(eS_Real);
					multem::add_square_scale(w, wave_function.psi_z, m2psi_tot);
				}

				if(input_multislice->coherent_contribution)
				{
					multem::add_scale(w, wave_function.psi_z, psi_coh);
				}
			}

			void get_coherent_intensity(Vector<value_type_c, dev> &psi_coh, Vector<value_type_r, dev> &m2psi_coh)
			{
				if(!input_multislice->coherent_contribution)
				{
					return;
				}

				if(input_multislice->is_HRTEM_HCI_EFTEM())
				{
					microscope_effects.apply(psi_coh, m2psi_coh);
				}
				else
				{
					multem::assign_square(psi_coh, m2psi_coh);
				}
			}

			void tot_coh_FS_RS(Vector<value_type_r, dev> &m2psi_tot, Vector<value_type_c, dev> &psi_coh)
			{
				value_type_r w = input_multislice->ifp_nconf();

				reset_m2psi_tot_and_psi_coh(m2psi_tot, psi_coh);

				for(auto iconf=1; iconf <= input_multislice->fp_nconf; iconf++)
				{
					wave_function.move_atoms(iconf);
					get_total_intensity_and_coherent_wave(w, m2psi_tot, psi_coh);
				}
			}

			Input_Multislice<value_type_r, dev> *input_multislice;
			Stream<value_type_r, dev> stream;
			FFT2<value_type_r, dev> fft2;

			Microscope_Effects<value_type_r, dev> microscope_effects;
			Wave_Function<value_type_r, dev> wave_function;
			Energy_Loss<value_type_r, dev> energy_loss;

			Vector<value_type_c, dev> psi_coh;
			Vector<value_type_r, dev> m2psi_coh;
			Vector<value_type_r, dev> m2psi_tot;

			Vector<value_type_c, dev> psi_thk;
			Vector<value_type_c, dev> psi_kn;
		};

} // namespace multem

#endif