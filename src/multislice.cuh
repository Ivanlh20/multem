/*
 * This file iscan part of MULTEM.
 * Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * MULTEM iscan free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM iscan distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MULTISLICE_H
#define MULTISLICE_H

#include <fftw3.h>
#include "types.cuh"
#include "traits.cuh"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "energy_loss.cuh"
#include "wave_function.cuh"

namespace multem
{
	template<class T, eDevice dev>
	class Multislice
	{
		public:
			using value_type_r = T;
			using value_type_c = complex<T>;

			void set_input_data(Input_Multislice<value_type_r, dev> *input_multislice_i)
			{
				input_multislice = input_multislice_i;
				if(dev == e_device)
				{
					cudaSetDevice(input_multislice->gpu_device);
				}
				stream.resize(input_multislice->nstream);
				fft2.create_plan(input_multislice->grid.ny, input_multislice->grid.nx, input_multislice->nstream);

				if(input_multislice->is_EELS_EFTEM())
				{
					energy_loss.set_input_data(input_multislice, &stream, &fft2);
					psi_thk.resize(input_multislice->grid.nxy());
					if(input_multislice->eels_fr.is_Double_Channelling_FOMS_SOMS())
					{
						trans_thk.resize(input_multislice->grid.nxy());
					}
				}

				wave_function.set_input_data(input_multislice, &stream, &fft2);
			}

			template<class TOutput_multislice>
			void run(TOutput_multislice &output_multislice)
			{
				if(input_multislice->is_STEM_ISTEM())
				{
					STEM_ISTEM(output_multislice);
				}
				else if(input_multislice->is_CBED_CBEI())
				{
					CBED_CBEI(output_multislice);
				}
				else if(input_multislice->is_ED_HRTEM())
				{
					ED_HRTEM(output_multislice);
				}
				else if(input_multislice->is_PED_HCI())
				{
					PED_HCI(output_multislice);
				}
				else if(input_multislice->is_EWFS_EWRS())
				{
					EWFS_EWRS(output_multislice);
				}
				else if(input_multislice->is_EELS_EFTEM())
				{
					EELS_EFTEM(output_multislice);
				}
			}

			void cleanup()
			{
				fft2.cleanup();
			}

		private:
			template<class TOutput_multislice>
			void STEM_ISTEM(TOutput_multislice &output_multislice)
			{
				value_type_r w = input_multislice->get_weight();

				output_multislice.init();

				if(input_multislice->is_STEM() && input_multislice->coherent_contribution)
				{
					for(auto iscan=0; iscan < input_multislice->scanning.size(); iscan++)
					{	
						input_multislice->set_beam_position(iscan);					
						for(auto iconf=input_multislice->fp_iconf_0; iconf <= input_multislice->fp_nconf; iconf++)
						{
							wave_function.move_atoms(iconf, input_multislice->tm_irot);		
							wave_function.psi_0(wave_function.psi_z);
							wave_function.psi(w, wave_function.psi_z, output_multislice);
						}

						wave_function.set_m2psi_coh(output_multislice);
					}
				}
				else
				{
					for(auto iconf=input_multislice->fp_iconf_0; iconf <= input_multislice->fp_nconf; iconf++)
					{
						wave_function.move_atoms(iconf, input_multislice->tm_irot);		
						for(auto iscan=0; iscan < input_multislice->scanning.size(); iscan++)
						{
							input_multislice->set_beam_position(iscan);
							wave_function.psi_0(wave_function.psi_z);
							wave_function.psi(w, wave_function.psi_z, output_multislice);
						}
					}

					wave_function.set_m2psi_coh(output_multislice);
				}

				output_multislice.shift();
				output_multislice.clear_temporal_data();
			}

			template<class TOutput_multislice>
			void CBED_CBEI(TOutput_multislice &output_multislice)
			{
				EWFS_EWRS(output_multislice);
			}

			template<class TOutput_multislice>
			void ED_HRTEM(TOutput_multislice &output_multislice)
			{
				EWFS_EWRS(output_multislice);
			}

			template<class TOutput_multislice>
			void PED_HCI(TOutput_multislice &output_multislice)
			{
				value_type_r w = input_multislice->get_weight();

				output_multislice.init();

				for(auto iconf=input_multislice->fp_iconf_0; iconf <= input_multislice->fp_nconf; iconf++)
				{
					wave_function.move_atoms(iconf, input_multislice->tm_irot);		
					for(auto irot=0; irot < input_multislice->nrot; irot++)
					{
						input_multislice->set_phi(irot);
						wave_function.psi_0(wave_function.psi_z);
						wave_function.psi(w, wave_function.psi_z, output_multislice);
					}
				}

				wave_function.set_m2psi_coh(output_multislice);

				output_multislice.shift();
				output_multislice.clear_temporal_data();
			}

			template<class TOutput_multislice>
			void EWFS_EWRS(TOutput_multislice &output_multislice)
			{
				value_type_r w = input_multislice->get_weight();

				output_multislice.init();

				for(auto iconf=input_multislice->fp_iconf_0; iconf <= input_multislice->fp_nconf; iconf++)
				{
					wave_function.move_atoms(iconf, input_multislice->tm_irot);		
					wave_function.psi_0(wave_function.psi_z);
					wave_function.psi(w, wave_function.psi_z, output_multislice);
				}

				wave_function.set_m2psi_coh(output_multislice);

				output_multislice.shift();
				output_multislice.clear_temporal_data();
			}

			template<class TOutput_multislice>
			void EELS_EFTEM(TOutput_multislice &output_multislice)
			{
				value_type_r w = input_multislice->get_weight();

				auto psi = [&](value_type_r w, Vector<value_type_c, dev> &psi_z, TOutput_multislice &output_multislice)
				{
					value_type_r gx_0 = input_multislice->gx_0();
					value_type_r gy_0 = input_multislice->gy_0();

					for(auto islice=0; islice < wave_function.slice.size(); islice++)
					{
						if(input_multislice->eels_fr.is_Double_Channelling_FOMS_SOMS())
						{
							wave_function.trans(islice, wave_function.slice.size()-1, trans_thk);
						}			

						for(auto iatom=wave_function.slice.iatom_0[islice]; iatom <= wave_function.slice.iatom_e[islice]; iatom++)
						{
							if(wave_function.atoms.Z[iatom] == input_multislice->eels_fr.Z)
							{
								input_multislice->set_eels_fr_atom(iatom, wave_function.atoms);
								energy_loss.set_atom_type(input_multislice->eels_fr);

								for(auto ikn=0; ikn < energy_loss.kernel.size(); ikn++)
								{
									multem::multiply(energy_loss.kernel[ikn], psi_z, wave_function.psi_z);
									wave_function.psi(islice, wave_function.slice.size()-1, w, trans_thk, output_multislice);
								}
							}
						}
						wave_function.psi_slice(gx_0, gy_0, islice, psi_z);
					}
				};

				output_multislice.init();

				if(input_multislice->is_EELS())
				{
					for(auto iconf=input_multislice->fp_iconf_0; iconf <= input_multislice->fp_nconf; iconf++)
					{
						wave_function.move_atoms(iconf, input_multislice->tm_irot);		
						for(auto iscan=0; iscan < input_multislice->scanning.size(); iscan++)
						{
							input_multislice->set_beam_position(iscan);
							wave_function.psi_0(psi_thk);
							psi(w, psi_thk, output_multislice);
						}
					}
				}
				else
				{
					for(auto iconf=input_multislice->fp_iconf_0; iconf <= input_multislice->fp_nconf; iconf++)
					{
						wave_function.move_atoms(iconf, input_multislice->tm_irot);		
						wave_function.psi_0(psi_thk);
						psi(w, psi_thk, output_multislice);
					}
				}

				output_multislice.shift();
				output_multislice.clear_temporal_data();
			}

			Input_Multislice<value_type_r, dev> *input_multislice;
			Stream<dev> stream;
			FFT2<value_type_r, dev> fft2;

			Wave_Function<value_type_r, dev> wave_function;
			Energy_Loss<value_type_r, dev> energy_loss;

			Vector<value_type_c, dev> psi_thk;
			Vector<value_type_c, dev> trans_thk;
		};

} // namespace multem

#endif