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
#include "types.hpp"
#include "traits.cuh"
#include "atom_data.hpp"
#include "input_multislice.hpp"

namespace multem
{
	template<class TVector_r, class TVector_c>
	class Output_Multislice: public Input_Multislice<Value_type<TVector_r>, e_Host>
	{
		public:
			using value_type_r = Value_type<TVector_r>;

			Output_Multislice(): Input_Multislice<value_type_r, e_Host>(), output_type(0), nx(0), ny(0), dx(0), dy(0){}

			//template<class TInput_Multislice>
			//enable_if_host_vector<TVector_r, void>
			//set_input_data(TInput_Multislice *input_multislice_i)
			//{ 
			//	set_input_multislice(*input_multislice_i);

			//	switch(output_type)
			//	{
			//		case 1:
			//		{
			//			m2psi.resize(this->thickness.size());
			//			psi.resize(this->thickness.size());
			//		}
			//		case 2:
			//		{
			//			m2psi.resize(this->thickness.size());
			//		}
			//		case 3:
			//		{
			//			probe.resize(this->thickness.size());
			//		}
			//		case 4:
			//		{
			//			V.resize(this->thickness.size());
			//		}
			//		case 5:
			//		{
			//			trans.resize(this->thickness.size());
			//		}
			//		case 6:
			//		{
			//			psi.resize(this->thickness.size());					
			//		}
			//		break;
			//	}

			//	for(auto i=0; i < m2psi.size(); i++)
			//	{
			//		switch(output_type)
			//		{
			//			case 1:
			//			{
			//				m2psi[i].resize(nxy());
			//				psi[i].resize(nxy());
			//			}
			//			case 2:
			//			{
			//				m2psi[i].resize(nxy());
			//			}
			//			case 3:
			//			{
			//				probe[i].resize(nxy());
			//			}
			//			case 4:
			//			{
			//				V[i].resize(nxy());
			//			}
			//			case 5:
			//			{
			//				trans[i].resize(nxy());
			//			}
			//			case 6:
			//			{
			//				psi[i].resize(nxy());							
			//			}
			//			break;
			//		}
			//	}
			//}

			//enable_if_rmatrix<TVector_r, void>
			template<class TInput_Multislice>
			void set_input_data(TInput_Multislice *input_multislice_i)
			{ 
				set_input_multislice(*input_multislice_i);

				switch(output_type)
				{
					case 1:
					{
						m2psi.resize(this->thickness.size());
						psi.resize(this->thickness.size());
					}
					case 2:
					{
						m2psi.resize(this->thickness.size());
					}
					case 3:
					{
						probe.resize(this->thickness.size());
					}
					case 4:
					{
						V.resize(this->thickness.size());
					}
					case 5:
					{
						trans.resize(this->thickness.size());
					}
					case 6:
					{
						psi.resize(this->thickness.size());					
					}
					break;
				}
			}

			int nxy() const { return nx*ny; }

			bool is_grid_FS() const
			{
				return is_CBED() || is_ED() || is_PED() || is_EWFS() ||
					is_ProbeFS() || is_PPFS() || is_TFFS() || is_EWSFS(); 
			}

			bool is_grid_RS() const
			{
				return !is_grid_FS();
			}

			int output_type;	// 1:(m2psi, psi); 2:(m2psi); 3:(probe); 4:(V); 5:(trans); 6:(psi)
			int nx;
			int ny;
			value_type_r dx;
			value_type_r dy;

			Vector<value_type_r, e_Host> x;
			Vector<value_type_r, e_Host> y;

			Vector<TVector_r, e_Host> m2psi;
			Vector<TVector_c, e_Host> psi;
			Vector<TVector_r, e_Host> V;
			Vector<TVector_c, e_Host> trans;
			Vector<TVector_c, e_Host> probe;

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
				this->pe_fr = input_multislice.pe_fr;
				this->ew_fr = input_multislice.ew_fr;
				this->beam_type = input_multislice.beam_type;
				this->conv_beam_wave_x = input_multislice.conv_beam_wave_x;
				this->conv_beam_wave_y = input_multislice.conv_beam_wave_y;
				this->islice = input_multislice.islice;
				this->nstream = input_multislice.nstream;
				this->dp_Shift = input_multislice.dp_Shift;

				m2psi.clear();
				m2psi.shrink_to_fit();

				psi.clear();
				psi.shrink_to_fit();

				V.clear();
				V.shrink_to_fit();

				trans.clear();
				trans.shrink_to_fit();

				probe.clear();
				probe.shrink_to_fit();
				
				if(this->is_scanning())
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

				// 1:(m2psi, psi); 2:(m2psi); 3:(probe); 4:(V); 5:(trans); 6:(psi)
				if(this->is_STEM_ISTEM()||this->is_CBED_CBEI() || this->is_ED_HRTEM()||
					this->is_PED_HCI() || this->is_EWFS_EWRS() || this->is_EELS_EFTEM())
				{
					output_type = (this->coherent_contribution)?1:2;
				}
				else if(this->is_ProbeFS_ProbeRS())
				{
					output_type = 3;
				}
				else if(this->is_PPFS_PPRS())
				{
					output_type = 4;
				}
				else if(this->is_TFFS_TFRS())
				{
					output_type = 5;
				}
				else if(this->is_EWSFS_EWSRS())
				{
					output_type = 6;
				}
			}

	};

} //namespace multem

#endif