/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef OUTPUT_MULTISLICE_H
#define OUTPUT_MULTISLICE_H

#include <algorithm>
#include <vector>

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "atom_data.hpp"
#include "input_multislice.cuh"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"

namespace multem
{
	template<class TVector_r, class TVector_c>
	class Output_Multislice: public Input_Multislice<Value_type<TVector_r>, e_host>
	{
		public:
			using value_type_r = Value_type<TVector_r>;
			using value_type_c = Value_type<TVector_c>;
			static const bool is_vector = is_host_device_vector<TVector_c>::value;

			Output_Multislice(): Input_Multislice<value_type_r, e_host>(), output_type(0), nx(0), ny(0), dx(0), dy(0){}

			template<class TInput_Multislice>
			void set_input_data(TInput_Multislice *input_multislice_i)
			{ 
				set_input_multislice(*input_multislice_i);

				stream.resize(this->cpu_nthread);

				// 1:(image_tot, image_coh); 2:(image_tot); 3:(m2psi_tot, m2psi_coh); 4:(m2psi_tot); 5:(m2psi_tot, psi_coh); 6:(psi_coh); 7:(probe); 8:(V); 9:(trans)
				switch(output_type)
				{
					case 1:
					{
						image_tot.resize(this->thickness.size());
						image_coh.resize(this->thickness.size());
						int ndet = (this->is_EELS())?1:this->det_cir.size();
						for(auto ithk=0; ithk < this->thickness.size(); ithk++)
						{
							image_tot[ithk].image.resize(ndet);
							image_coh[ithk].image.resize(ndet);
						}
						psi_coh.resize(this->thickness.size());
					}
					break;
					case 2:
					{
						image_tot.resize(this->thickness.size());
						int ndet = (this->is_EELS())?1:this->det_cir.size();
						for(auto ithk=0; ithk < this->thickness.size(); ithk++)
						{
							image_tot[ithk].image.resize(ndet);
						}
					}
					break;
					case 3:
					{
						m2psi_tot.resize(this->thickness.size());
						m2psi_coh.resize(this->thickness.size());
						psi_coh.resize(this->thickness.size());
					}
					break;
					case 4:
					{
						m2psi_tot.resize(this->thickness.size());
					}
					break;
					case 5:
					{
						m2psi_tot.resize(this->thickness.size());
						psi_coh.resize(this->thickness.size());										
					}
					break;
					case 6:
					{
						psi_coh.resize(this->thickness.size());					
					}
					break;
					case 7:
					{
						probe.resize(this->thickness.size());
					}
					break;
					case 8:
					{
						V.resize(this->thickness.size());
					}
					break;
					case 9:
					{
						trans.resize(this->thickness.size());
					}
					break;
				}

				if(is_vector)
				{
					for(auto ithk=0; ithk < this->thickness.size(); ithk++)
					{
						switch(output_type)
						{
							case 1:
							{
								int ndet = (this->is_EELS())?1:this->det_cir.size();
								for(auto idet=0; idet < ndet; idet++)
								{
									image_tot[ithk].image[idet].resize(nxy());
									image_coh[ithk].image[idet].resize(nxy());
								}
								psi_coh[ithk].resize(nxy());
							}
							break;
							case 2:
							{
								int ndet = (this->is_EELS())?1:this->det_cir.size();
								for(auto idet=0; idet < ndet; idet++)
								{
									image_tot[ithk].image[idet].resize(nxy());
								}
							}
							break;
							case 3:
							{
								m2psi_tot[ithk].resize(nxy());
								m2psi_coh[ithk].resize(nxy());
								psi_coh[ithk].resize(nxy());
							}
							break;
							case 4:
							{
								m2psi_tot[ithk].resize(nxy());
							}
							break;
							case 5:
							{
								m2psi_tot[ithk].resize(nxy());
								psi_coh[ithk].resize(nxy());							
							}
							break;
							case 6:
							{
								psi_coh[ithk].resize(nxy());							
							}
							break;
							case 7:
							{
								probe[ithk].resize(nxy());
							}
							break;
							case 8:
							{
								V[ithk].resize(nxy());
							}
							break;
							case 9:
							{
								trans[ithk].resize(nxy());
							}
							break;
						}
					}
				}
				else
				{
					if((output_type==1)||(output_type==3))
					{
						for(auto ithk=0; ithk < this->thickness.size(); ithk++)
						{
							psi_coh[ithk].resize(nxy());
						}
					}
				}
			}

			void init()
			{ 
				for(auto ithk=0; ithk < this->thickness.size(); ithk++)
				{
					switch(output_type)
					{
						case 1:
						{
							for(auto idet=0; idet < image_tot[ithk].image.size(); idet++)
							{
								multem::fill(image_tot[ithk].image[idet], 0);
								multem::fill(image_coh[ithk].image[idet], 0);
							}
							multem::fill(psi_coh[ithk], value_type_c(0));
						}
						break;
						case 2:
						{
							for(auto idet=0; idet < image_tot[ithk].image.size(); idet++)
							{
								multem::fill(image_tot[ithk].image[idet], 0);
							}
						}
						break;
						case 3:
						{
							multem::fill(m2psi_tot[ithk], 0);
							multem::fill(psi_coh[ithk], value_type_c(0));
						}
						break;
						case 4:
						{
							multem::fill(m2psi_tot[ithk], 0);
						}
						break;
						case 5:
						{
							multem::fill(m2psi_tot[ithk], 0);
							multem::fill(psi_coh[ithk], value_type_c(0));
						}
						break;
						case 6:
						{
							multem::fill(psi_coh[ithk], value_type_c(0));							
						}
						break;
						case 7:
						{
							multem::fill(probe[ithk], value_type_c(0));
						}
						break;
						case 8:
						{
							multem::fill(V[ithk], 0);
						}
						break;
						case 9:
						{
							multem::fill(trans[ithk], value_type_c(0));
						}
						break;
					}
				}
			}

			void shift()
			{ 
				for(auto ithk=0; ithk < this->thickness.size(); ithk++)
				{
					switch(output_type)
					{
						case 3:
						{
							multem::fft2_shift(stream, this->grid, m2psi_tot[ithk]);
							multem::fft2_shift(stream, this->grid, m2psi_coh[ithk]);
						}
						break;
						case 4:
						{
							multem::fft2_shift(stream, this->grid, m2psi_tot[ithk]);
						}
						break;
						case 5:
						{
							multem::fft2_shift(stream, this->grid, m2psi_tot[ithk]);
							multem::fft2_shift(stream, this->grid, psi_coh[ithk]);
						}
						break;
						case 6:
						{
							multem::fft2_shift(stream, this->grid, psi_coh[ithk]);							
						}
						break;
						case 7:
						{
							multem::fft2_shift(stream, this->grid, probe[ithk]);
						}
						break;
						case 8:
						{
							multem::fft2_shift(stream, this->grid, V[ithk]);
						}
						break;
						case 9:
						{
							multem::fft2_shift(stream, this->grid, trans[ithk]);
						}
						break;
					}
				}
			}

			void clear_temporal_data()
			{
				if((output_type==1)||(output_type==3))
				{
					for(auto ithk=0; ithk < this->thickness.size(); ithk++)
					{
						psi_coh[ithk].clear();
					}

					psi_coh.clear();
					psi_coh.shrink_to_fit();
				}
			}

			int nxy() const { return nx*ny; }

			bool is_grid_FS() const
			{
				return this->is_CBED() || this->is_ED() || this->is_PED() || this->is_EWFS() ||
					this->is_ProbeFS() || this->is_PPFS() || this->is_TFFS() || this->is_EWSFS(); 
			}

			bool is_grid_RS() const
			{
				return !is_grid_FS();
			}

			// 1:(image_tot, image_coh); 2:(image_tot); 3:(m2psi_tot, m2psi_coh); 4:(m2psi_tot); 5:(m2psi_tot, psi_coh); 6:(psi_coh); 7:(probe); 8:(V); 9:(trans)
			int output_type;	
			int nx;
			int ny;
			value_type_r dx;
			value_type_r dy;

			Vector<value_type_r, e_host> x;
			Vector<value_type_r, e_host> y;

			Vector<Det_Int<TVector_r>, e_host> image_tot;
			Vector<Det_Int<TVector_r>, e_host> image_coh;
			Vector<TVector_r, e_host> m2psi_tot;
			Vector<TVector_r, e_host> m2psi_coh;
			Vector<TVector_c, e_host> psi_coh;
			Vector<TVector_r, e_host> V;
			Vector<TVector_c, e_host> trans;
			Vector<TVector_c, e_host> probe;

			Stream<e_host> stream;
		private:
			template<class TInput_Multislice>
			void set_input_multislice(TInput_Multislice &input_multislice)
			{ 
				this->precision = input_multislice.precision;
				this->device = input_multislice.device;
				this->cpu_ncores = input_multislice.cpu_ncores;
				this->cpu_nthread = input_multislice.cpu_nthread;
				this->gpu_device = input_multislice.gpu_nstream;
				this->simulation_type = input_multislice.simulation_type;
				this->phonon_model = input_multislice.phonon_model;
				this->interaction_model = input_multislice.interaction_model;
				this->potential_slicing = input_multislice.potential_slicing;
				this->potential_type = input_multislice.potential_type;

				this->fp_dim = input_multislice.fp_dim;
				this->fp_dist = input_multislice.fp_dist;
				this->fp_seed = input_multislice.fp_seed;
				this->fp_single_conf = input_multislice.fp_single_conf;
				this->fp_nconf = input_multislice.fp_nconf;

				this->tm_active = input_multislice.tm_active;
				this->tm_nrot = input_multislice.tm_nrot;
				this->tm_irot = input_multislice.tm_irot;
				this->tm_theta_0 = input_multislice.tm_theta_0;
				this->tm_theta_e = input_multislice.tm_theta_e;
				this->tm_u0 = input_multislice.tm_u0;
				this->tm_rot_point_type = input_multislice.tm_rot_point_type;
				this->tm_p0 = input_multislice.tm_p0;

				this->microscope_effect = input_multislice.microscope_effect;
				this->spatial_temporal_effect = input_multislice.spatial_temporal_effect;
				this->zero_defocus_type = input_multislice.zero_defocus_type;
				this->zero_defocus_plane = input_multislice.zero_defocus_plane;
				this->thickness_type = input_multislice.thickness_type;
				this->thickness = input_multislice.thickness;
				this->input_wave_type = input_multislice.input_wave_type;
				//this->psi_0 = input_multislice.psi_0;
				this->operation_mode = input_multislice.operation_mode;
				this->coherent_contribution = input_multislice.coherent_contribution;
				this->slice_storage = input_multislice.slice_storage;
				this->E_0 = input_multislice.E_0;
				this->theta = input_multislice.theta;
				this->phi = input_multislice.phi;
				this->grid = input_multislice.grid;
				this->Vrl = input_multislice.Vrl;
				this->nR = input_multislice.nR;
				this->lens = input_multislice.lens;
				this->is_crystal = input_multislice.is_crystal;
				//this->atoms = input_multislice.atoms;
				this->eels_fr = input_multislice.eels_fr;
				this->scanning = input_multislice.scanning;
				this->det_cir = input_multislice.det_cir;
				this->hrtem = input_multislice.hrtem;
				this->cbe_fr = input_multislice.cbe_fr;
				this->ew_fr = input_multislice.ew_fr;
				this->beam_type = input_multislice.beam_type;
				this->conv_beam_wave_x = input_multislice.conv_beam_wave_x;
				this->conv_beam_wave_y = input_multislice.conv_beam_wave_y;
				this->islice = input_multislice.islice;
				this->nstream = input_multislice.nstream;
				this->dp_Shift = input_multislice.dp_Shift;

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

				V.clear();
				V.shrink_to_fit();

				trans.clear();
				trans.shrink_to_fit();

				probe.clear();
				probe.shrink_to_fit();
				
				if(this->is_STEM() || this->is_EELS())
				{
					nx = this->scanning.nx;
					ny = this->scanning.ny;
					dx = this->scanning.dRx;
					dy = this->scanning.dRy;

					x.resize(nx);
					y.resize(ny);

					for(auto ix=0; ix<nx; ix++)
					{
						x[ix] = this->scanning.Rx(ix);
					}
					for(auto iy=0; iy<ny; iy++)
					{
						y[iy] = this->scanning.Ry(iy);
					}
				}
				else
				{
					nx = this->grid.nx;
					ny = this->grid.ny;
					dx = (is_grid_RS())?this->grid.dRx:this->grid.dgx;
					dy = (is_grid_RS())?this->grid.dRy:this->grid.dgy;

					x.resize(nx);
					y.resize(ny);

					for(auto ix=0; ix<nx; ix++)
					{
						x[ix] = (is_grid_RS())?this->grid.Rx(ix):this->grid.gx(ix);
					}
					for(auto iy=0; iy<ny; iy++)
					{
						y[iy] = (is_grid_RS())?this->grid.Ry(iy):this->grid.gy(iy);
					}
				}

				// 1:(image_tot, image_coh); 2:(image_tot); 3:(m2psi_tot, m2psi_coh); 4:(m2psi_tot); 5:(m2psi_tot, psi_coh); 6:(psi_coh); 7:(probe); 8:(V); 9:(trans)
				if(this->is_STEM() || this->is_EELS())
				{
					output_type = (this->coherent_contribution)?1:2;
				}
				else if(this->is_ISTEM() || this->is_CBED_CBEI() ||this->is_PED_HCI() || 
					this->is_ED_HRTEM() || this->is_EFTEM())
				{
					output_type = (this->coherent_contribution)?3:4;
				}
				else if(this->is_EWFS_EWRS())
				{
					output_type = 5;
				}
				else if(this->is_EWSFS_EWSRS() || this->is_PropFS_PropRS())
				{
					output_type = 6;
				}
				else if(this->is_ProbeFS_ProbeRS())
				{
					output_type = 7;
				}
				else if(this->is_PPFS_PPRS())
				{
					output_type = 8;
				}
				else if(this->is_TFFS_TFRS())
				{
					output_type = 9;
				}
			}

	};

} //namespace multem

#endif