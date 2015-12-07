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
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef OUTPUT_MULTISLICE_H
#define OUTPUT_MULTISLICE_H

#include <algorithm>
#include <vector>

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "atom_data.hpp"
#include "input_multislice.cuh"
#include "host_functions.hpp"
#include "device_functions.cuh"

namespace multem
{
	template<class TVector_r, class TVector_c>
	class Output_Multislice: public Input_Multislice<Value_type<TVector_r>>
	{
		public:
			using value_type_r = Value_type<TVector_r>;
			using value_type_c = Value_type<TVector_c>;
			static const bool is_vector = is_host_device_vector<TVector_c>::value;

			Output_Multislice(): Input_Multislice<value_type_r>(), output_type(0), ndetector(0), nx(0), ny(0), dx(0), dy(0), dr(0){}

			template<class TInput_Multislice>
			void set_input_data(TInput_Multislice *input_multislice_i)
			{ 
				set_input_multislice(*input_multislice_i);

				stream.resize(this->cpu_nthread);

				// 1:(image_tot, image_coh); 2:(image_tot); 3:(m2psi_tot, m2psi_coh); 4:(m2psi_tot); 5:(m2psi_tot, psi_coh); 6:(psi_coh); 7:(psi_0); 8:(V); 9:(trans)
				switch(output_type)
				{
					case 1:
					{
						image_tot.resize(this->thickness.size());
						image_coh.resize(this->thickness.size());
						for(auto ithk =0; ithk < this->thickness.size(); ithk++)
						{
							image_tot[ithk].image.resize(ndetector);
							image_coh[ithk].image.resize(ndetector);
						}
						psi_coh.resize(this->thickness.size());
					}
					break;
					case 2:
					{
						image_tot.resize(this->thickness.size());
						for(auto ithk =0; ithk < this->thickness.size(); ithk++)
						{
							image_tot[ithk].image.resize(ndetector);
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
						psi_0.resize(this->thickness.size());
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
					for(auto ithk =0; ithk < this->thickness.size(); ithk++)
					{
						switch(output_type)
						{
							case 1:
							{
								for(auto idet =0; idet < ndetector; idet++)
								{
									image_tot[ithk].image[idet].resize(nxy());
									image_coh[ithk].image[idet].resize(nxy());
								}
								psi_coh[ithk].resize(nxy());
							}
							break;
							case 2:
							{
								for(auto idet =0; idet < ndetector; idet++)
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
								psi_0[ithk].resize(nxy());
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
					if((output_type == 1)||(output_type == 3))
					{
						for(auto ithk =0; ithk < this->thickness.size(); ithk++)
						{
							psi_coh[ithk].resize(nxy());
						}
					}
				}
			}

			void init()
			{ 
				for(auto ithk =0; ithk < this->thickness.size(); ithk++)
				{
					switch(output_type)
					{
						case 1:
						{
							for(auto idet =0; idet < image_tot[ithk].image.size(); idet++)
							{
								multem::fill(stream, image_tot[ithk].image[idet], 0);
								multem::fill(stream, image_coh[ithk].image[idet], 0);
							}
							multem::fill(stream, psi_coh[ithk], 0);
						}
						break;
						case 2:
						{
							for(auto idet =0; idet < image_tot[ithk].image.size(); idet++)
							{
								multem::fill(stream, image_tot[ithk].image[idet], 0);
							}
						}
						break;
						case 3:
						{
							multem::fill(stream, m2psi_tot[ithk], 0);
							multem::fill(stream, psi_coh[ithk], 0);
						}
						break;
						case 4:
						{
							multem::fill(stream, m2psi_tot[ithk], 0);
						}
						break;
						case 5:
						{
							multem::fill(stream, m2psi_tot[ithk], 0);
							multem::fill(stream, psi_coh[ithk], 0);
						}
						break;
						case 6:
						{
							multem::fill(stream, psi_coh[ithk], 0);							
						}
						break;
						case 7:
						{
							multem::fill(stream, psi_0[ithk], 0);
						}
						break;
						case 8:
						{
							multem::fill(stream, V[ithk], 0);
						}
						break;
						case 9:
						{
							multem::fill(stream, trans[ithk], 0);
						}
						break;
					}
				}
			}

			void shift()
			{ 
				for(auto ithk =0; ithk < this->thickness.size(); ithk++)
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
							multem::fft2_shift(stream, this->grid, psi_0[ithk]);
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
				if((output_type == 1)||(output_type == 3))
				{
					for(auto ithk =0; ithk < this->thickness.size(); ithk++)
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
					this->is_IWFS() || this->is_PPFS() || this->is_TFFS(); 
			}

			bool is_grid_RS() const
			{
				return !is_grid_FS();
			}

			// 1:(image_tot, image_coh); 2:(image_tot); 3:(m2psi_tot, m2psi_coh); 4:(m2psi_tot); 5:(m2psi_tot, psi_coh); 6:(psi_coh); 7:(psi_0); 8:(V); 9:(trans)
			int output_type;	
			int ndetector;
			int nx;
			int ny;
			value_type_r dx;
			value_type_r dy;
			value_type_r dr;

			Vector<value_type_r, e_host> x;
			Vector<value_type_r, e_host> y;
			Vector<value_type_r, e_host> r;

			Vector<Det_Int<TVector_r>, e_host> image_tot;
			Vector<Det_Int<TVector_r>, e_host> image_coh;
			Vector<TVector_r, e_host> m2psi_tot;
			Vector<TVector_r, e_host> m2psi_coh;
			Vector<TVector_c, e_host> psi_coh;
			Vector<TVector_r, e_host> V;
			Vector<TVector_c, e_host> trans;
			Vector<TVector_c, e_host> psi_0;

			Stream<e_host> stream;
		private:
			template<class TInput_Multislice>
			void set_input_multislice(TInput_Multislice &input_multislice)
			{ 
				this->precision = input_multislice.precision;
				this->device = input_multislice.device;
				this->cpu_ncores = input_multislice.cpu_ncores;
				this->cpu_nthread = input_multislice.cpu_nthread;
				this->gpu_device = input_multislice.gpu_device;
				this->gpu_nstream = input_multislice.gpu_nstream;
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
				this->tm_theta = input_multislice.tm_theta;
				this->tm_u0 = input_multislice.tm_u0;
				this->tm_rot_point_type = input_multislice.tm_rot_point_type;
				this->tm_p0 = input_multislice.tm_p0;

				this->microscope_effect = input_multislice.microscope_effect;
				this->spatial_temporal_effect = input_multislice.spatial_temporal_effect;
				this->zero_defocus_type = input_multislice.zero_defocus_type;
				this->zero_defocus_plane = input_multislice.zero_defocus_plane;
				this->thickness_type = input_multislice.thickness_type;
				this->thickness = input_multislice.thickness;
				this->operation_mode = input_multislice.operation_mode;
				this->coherent_contribution = input_multislice.coherent_contribution;
				this->slice_storage = input_multislice.slice_storage;
				this->E_0 = input_multislice.E_0;
				this->theta = input_multislice.theta;
				this->phi = input_multislice.phi;
				this->grid = input_multislice.grid;
				this->Vrl = input_multislice.Vrl;
				this->nR = input_multislice.nR;
				this->iw_type = input_multislice.iw_type;
				// this->iw_psi = input_multislice.iw_psi;
				this->iw_x = input_multislice.iw_x;
				this->iw_y = input_multislice.iw_y;
				this->lens = input_multislice.lens;
				this->is_crystal = input_multislice.is_crystal;
				// this->atoms = input_multislice.atoms;
				this->eels_fr = input_multislice.eels_fr;
				this->scanning = input_multislice.scanning;
				// this->detector = input_multislice.detector;
				ndetector = (input_multislice.is_EELS())?1:input_multislice.detector.size();
				this->beam_x = input_multislice.beam_x;
				this->beam_y = input_multislice.beam_y;
				this->iscan = input_multislice.iscan;
				this->islice = input_multislice.islice;
				this->dp_Shift = input_multislice.dp_Shift;
				this->nstream = input_multislice.nstream;

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

				psi_0.clear();
				psi_0.shrink_to_fit();
				
				if(this->is_STEM() || this->is_EELS())
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
					nx = this->grid.nx;
					ny = this->grid.ny;
					dx = (is_grid_RS())?this->grid.dRx:this->grid.dgx;
					dy = (is_grid_RS())?this->grid.dRy:this->grid.dgy;

					x.resize(nx);
					y.resize(ny);

					for(auto ix =0; ix<nx; ix++)
					{
						x[ix] = (is_grid_RS())?this->grid.Rx(ix):this->grid.gx(ix);
					}
					for(auto iy =0; iy<ny; iy++)
					{
						y[iy] = (is_grid_RS())?this->grid.Ry(iy):this->grid.gy(iy);
					}
				}

				// 1:(image_tot, image_coh); 2:(image_tot); 3:(m2psi_tot, m2psi_coh); 4:(m2psi_tot); 5:(m2psi_tot, psi_coh); 6:(psi_coh); 7:(psi_0); 8:(V); 9:(trans)
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
					output_type = (this->is_EWFS_EWRS_SC())?6:5;
				}
				else if(this->is_PropFS_PropRS())
				{
					output_type = 6;
				}
				else if(this->is_IWFS_IWRS())
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

	using Output_Multislice_Matlab = Output_Multislice<rmatrix_r, rmatrix_c>;

	template<class T>
	using Output_Multislice_Vector = Output_Multislice<Vector<T, e_host>, Vector<complex<T>, e_host>>;
} // namespace multem

#endif