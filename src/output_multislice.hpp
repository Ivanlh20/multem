/*
 * This file is part of MULTEM.
 * Copyright 2017 Ivan Lobato <Ivanlh20@gmail.com>
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

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "matlab_types.cuh"
#include "atom_data.hpp"
#include "input_multislice.cuh"
#include "host_functions.hpp"
#include "device_functions.cuh"

namespace mt
{
	inline
	bool is_ot_image_tot_coh(const eTEM_Output_Type &output_type)
	{
		return output_type==mt::eTEMOT_image_tot_coh;
	}

	inline
	bool is_ot_image_tot(const eTEM_Output_Type &output_type)
	{
		return output_type==mt::eTEMOT_image_tot;
	}

	inline
	bool is_ot_m2psi_tot_coh(const eTEM_Output_Type &output_type)
	{
		return output_type==mt::eTEMOT_m2psi_tot_coh;
	}

	inline
	bool is_ot_m2psi_tot(const eTEM_Output_Type &output_type)
	{
		return output_type==mt::eTEMOT_m2psi_tot;
	}

	inline
	bool is_ot_m2psi_tot_psi_coh(const eTEM_Output_Type &output_type)
	{
		return output_type==mt::eTEMOT_m2psi_tot_psi_coh;
	}

	inline
	bool is_ot_psi_coh(const eTEM_Output_Type &output_type)
	{
		return output_type==mt::eTEMOT_psi_coh;
	}

	inline
	bool is_ot_psi_0(const eTEM_Output_Type &output_type)
	{
		return output_type==mt::eTEMOT_psi_0;
	}

	inline
	bool is_ot_V(const eTEM_Output_Type &output_type)
	{
		return output_type==mt::eTEMOT_V;
	}

	inline
	bool is_ot_trans(const eTEM_Output_Type &output_type)
	{
		return output_type==mt::eTEMOT_trans;
	}

	/**************************************************************************************/
	template <class T>
	class Output_Multislice: public Input_Multislice<T>
	{
		public:
			using T_r = T;
			using T_c = complex<T>;

			using TVector_r = host_vector<T>;
			using TVector_c = host_vector<complex<T>>;

			Output_Multislice(): Input_Multislice<T_r>(), output_type(eTEMOT_m2psi_tot), ndetector(0), nx(0), ny(0), dx(0), dy(0), dr(0){}

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
				for(auto ithk = 0; ithk < output_multislice.image_tot.size(); ithk++)
				{
					 image_tot[ithk].image.resize(output_multislice.image_tot[ithk].image.size());
					 for(auto idet = 0; idet < output_multislice.image_tot[ithk].image.size(); idet++)
					 {
						image_tot[ithk].image[idet] = output_multislice.image_tot[ithk].image[idet];
					 }
				}

				image_coh.resize(output_multislice.image_coh.size());
				for(auto ithk = 0; ithk < output_multislice.image_coh.size(); ithk++)
				{
					image_coh[ithk].image.resize(output_multislice.image_coh[ithk].image.size());
					for(auto idet = 0; idet < output_multislice.image_coh[ithk].image.size(); idet++)
					{
						image_coh[ithk].image[idet] = output_multislice.image_coh[ithk].image[idet];
					}
				}

				m2psi_tot.resize(output_multislice.m2psi_tot.size());
				for(auto ithk = 0; ithk < output_multislice.m2psi_tot.size(); ithk++)
				{
					m2psi_tot[ithk] = output_multislice.m2psi_tot[ithk];
				}

				m2psi_coh.resize(output_multislice.m2psi_coh.size());
				for(auto ithk = 0; ithk < output_multislice.m2psi_coh.size(); ithk++)
				{
					m2psi_coh[ithk] = output_multislice.m2psi_coh[ithk];
				}

				psi_coh.resize(output_multislice.psi_coh.size());
				for(auto ithk = 0; ithk < output_multislice.psi_coh.size(); ithk++)
				{
					psi_coh[ithk] = output_multislice.psi_coh[ithk];
				}

				V.resize(output_multislice.V.size());
				for(auto ithk = 0; ithk < output_multislice.V.size(); ithk++)
				{
					V[ithk] = output_multislice.V[ithk];
				}

				trans.resize(output_multislice.trans.size());
				for(auto ithk = 0; ithk < output_multislice.trans.size(); ithk++)
				{
					trans[ithk] = output_multislice.trans[ithk];
				}

				psi_0.resize(output_multislice.psi_0.size());
				for(auto ithk = 0; ithk < output_multislice.psi_0.size(); ithk++)
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

				V.clear();
				V.shrink_to_fit();

				trans.clear();
				trans.shrink_to_fit();

				psi_0.clear();
				psi_0.shrink_to_fit();
			}

			template <class TInput_Multislice>
			void set_input_data(TInput_Multislice *input_multislice)
			{ 
				clear();

				stream.resize(1);

				assign_input_multislice(*input_multislice);

				ndetector = (this->is_EELS())?1:this->detector.size();
				
				set_output_grid();

				set_output_type();

				switch(output_type)
				{
					case eTEMOT_image_tot_coh:
					{
						image_tot.resize(this->thick.size());
						image_coh.resize(this->thick.size());
						psi_coh.resize(this->thick.size());

						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							image_tot[ithk].image.resize(ndetector);
							image_coh[ithk].image.resize(ndetector);

							for(auto idet = 0; idet < ndetector; idet++)
							{
								image_tot[ithk].image[idet].resize(nxy());
								image_coh[ithk].image[idet].resize(nxy());

								mt::fill(stream, image_tot[ithk].image[idet], T_r(0));
								mt::fill(stream, image_coh[ithk].image[idet], T_r(0));
							}

							psi_coh[ithk].resize(this->grid_2d.nxy());
							mt::fill(stream, psi_coh[ithk], T_c(0));
						}
					}
					break;
					case eTEMOT_image_tot:
					{
						image_tot.resize(this->thick.size());

						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							image_tot[ithk].image.resize(ndetector);

							for(auto idet = 0; idet < ndetector; idet++)
							{
								image_tot[ithk].image[idet].resize(nxy());

								mt::fill(stream, image_tot[ithk].image[idet], T_r(0));
							}
						}
					}
					break;
					case eTEMOT_m2psi_tot_coh:
					{
						m2psi_tot.resize(this->thick.size());
						m2psi_coh.resize(this->thick.size());
						psi_coh.resize(this->thick.size());

						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							m2psi_tot[ithk].resize(nxy());
							m2psi_coh[ithk].resize(nxy());
							psi_coh[ithk].resize(nxy());

							mt::fill(stream, m2psi_tot[ithk], T_r(0));
							mt::fill(stream, m2psi_coh[ithk], T_r(0));
							mt::fill(stream, psi_coh[ithk], T_c(0));
						}
					}
					break;
					case eTEMOT_m2psi_tot:
					{
						m2psi_tot.resize(this->thick.size());

						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							m2psi_tot[ithk].resize(nxy());

							mt::fill(stream, m2psi_tot[ithk], T_r(0));
						}
					}
					break;
					case eTEMOT_m2psi_tot_psi_coh:
					{
						m2psi_tot.resize(this->thick.size());
						psi_coh.resize(this->thick.size());		

						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							m2psi_tot[ithk].resize(nxy());
							psi_coh[ithk].resize(nxy());

							mt::fill(stream, m2psi_tot[ithk], T_r(0));
							mt::fill(stream, psi_coh[ithk], T_c(0));
						}
					}
					break;
					case eTEMOT_psi_coh:
					{
						psi_coh.resize(this->thick.size());		

						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							psi_coh[ithk].resize(nxy());

							mt::fill(stream, psi_coh[ithk], T_c(0));
						}
					}
					break;
					case eTEMOT_psi_0:
					{
						psi_0.resize(this->thick.size());

						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							psi_0[ithk].resize(nxy());

							mt::fill(stream, psi_0[ithk], T_c(0));
						}
					}
					break;
					case eTEMOT_V:
					{
						V.resize(this->thick.size());

						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							V[ithk].resize(nxy());

							mt::fill(stream, V[ithk], T_r(0));
						}
					}
					break;
					case eTEMOT_trans:
					{
						trans.resize(this->thick.size());

						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							trans[ithk].resize(nxy());

							mt::fill(stream, trans[ithk], T_c(0));
						}
					}
					break;
				}
			}

			void init()
			{ 
				switch(output_type)
				{
					case eTEMOT_image_tot_coh:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							for(auto idet = 0; idet < image_tot[ithk].image.size(); idet++)
							{
								mt::fill(stream, image_tot[ithk].image[idet], T_r(0));
								mt::fill(stream, image_coh[ithk].image[idet], T_r(0));
							}
							mt::fill(stream, psi_coh[ithk], 0);
						}
					}
					break;
					case eTEMOT_image_tot:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							for(auto idet = 0; idet < image_tot[ithk].image.size(); idet++)
							{
								mt::fill(stream, image_tot[ithk].image[idet], T_r(0));
							}
						}
					}
					break;
					case eTEMOT_m2psi_tot_coh:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							mt::fill(stream, m2psi_tot[ithk], T_r(0));
							mt::fill(stream, m2psi_coh[ithk], T_r(0));
							mt::fill(stream, psi_coh[ithk], T_c(0));
						}
					}
					break;
					case eTEMOT_m2psi_tot:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							mt::fill(stream, m2psi_tot[ithk], T_r(0));
						}
					}
					break;
					case eTEMOT_m2psi_tot_psi_coh:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							mt::fill(stream, m2psi_tot[ithk], T_r(0));
							mt::fill(stream, psi_coh[ithk], T_c(0));
						}
					}
					break;
					case eTEMOT_psi_coh:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							mt::fill(stream, psi_coh[ithk], T_c(0));		
						}
					}
					break;
					case eTEMOT_psi_0:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							mt::fill(stream, psi_0[ithk], T_c(0));
						}
					}
					break;
					case eTEMOT_V:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							mt::fill(stream, V[ithk], T_r(0));
						}
					}
					break;
					case eTEMOT_trans:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							mt::fill(stream, trans[ithk], T_c(0));
						}
					}
					break;
				}
			}

			void init_psi_coh()
			{
				for(auto ithk = 0; ithk < this->thick.size(); ithk++)
				{
					mt::fill(stream, psi_coh[ithk], 0);
				}
			}

			void shift()
			{ 
				switch(output_type)
				{
					case eTEMOT_m2psi_tot_coh:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							mt::fft2_shift(stream, this->grid_2d, m2psi_tot[ithk]);
							mt::fft2_shift(stream, this->grid_2d, m2psi_coh[ithk]);
						}
					}
					break;
					case eTEMOT_m2psi_tot:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							mt::fft2_shift(stream, this->grid_2d, m2psi_tot[ithk]);
						}
					}
					break;
					case eTEMOT_m2psi_tot_psi_coh:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							mt::fft2_shift(stream, this->grid_2d, m2psi_tot[ithk]);
							mt::fft2_shift(stream, this->grid_2d, psi_coh[ithk]);
						}
					}
					break;
					case eTEMOT_psi_coh:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							mt::fft2_shift(stream, this->grid_2d, psi_coh[ithk]);	
						}
					}
					break;
					case eTEMOT_psi_0:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							mt::fft2_shift(stream, this->grid_2d, psi_0[ithk]);
						}
					}
					break;
					case eTEMOT_V:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							mt::fft2_shift(stream, this->grid_2d, V[ithk]);
						}
					}
					break;
					case eTEMOT_trans:
					{
						for(auto ithk = 0; ithk < this->thick.size(); ithk++)
						{
							mt::fft2_shift(stream, this->grid_2d, trans[ithk]);
						}
					}
					break;
				}
			}

			void clear_temporal_data()
			{
				if(is_ot_image_tot_coh()||is_ot_m2psi_tot_coh())
				{
					for(auto ithk = 0; ithk < this->thick.size(); ithk++)
					{
						psi_coh[ithk].clear();
					}
				}
			}

			int nxy() const { return nx*ny; }

			/**************************************************************************************/
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

			host_vector<float> extract_data(ePhonon_Model_Output fp_ctr, eShow_CData show_data, int ithk, int idet=0)
			{
				host_vector<float> data(nxy());

				switch(output_type)
				{
					case eTEMOT_image_tot_coh:
					{
						data = (fp_ctr==eFMO_Total)?image_tot[ithk].image[idet]:image_coh[ithk].image[idet];
					}
					break;
					case eTEMOT_image_tot:
					{
						data = image_tot[ithk].image[idet];
					}
					break;
					case eTEMOT_m2psi_tot_coh:
					{
						data = (fp_ctr==eFMO_Total)?m2psi_tot[ithk]:m2psi_coh[ithk];
					}
					break;
					case eTEMOT_m2psi_tot:
					{
						 data = m2psi_tot[ithk];
					}
					break;
					case eTEMOT_m2psi_tot_psi_coh:
					{
						if(fp_ctr==eFMO_Total)
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

			host_vector<Det_Int<TVector_r>> image_tot;
			host_vector<Det_Int<TVector_r>> image_coh;
			host_vector<TVector_r> m2psi_tot;
			host_vector<TVector_r> m2psi_coh;
			host_vector<TVector_c> psi_coh;
			host_vector<TVector_r> V;
			host_vector<TVector_c> trans;
			host_vector<TVector_c> psi_0;

			Stream<e_host> stream;
		private:

			void set_output_grid()
			{
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
					nx = this->grid_2d.nx;
					ny = this->grid_2d.ny;
					const bool is_RS = this->is_grid_RS();
					dx = (is_RS)?this->grid_2d.dRx:this->grid_2d.dgx;
					dy = (is_RS)?this->grid_2d.dRy:this->grid_2d.dgy;

					x.resize(nx);
					y.resize(ny);

					for(auto ix = 0; ix<nx; ix++)
					{
						x[ix] = (is_RS)?this->grid_2d.Rx(ix):this->grid_2d.gx(ix);
					}

					for(auto iy = 0; iy<ny; iy++)
					{
						y[iy] = (is_RS)?this->grid_2d.Ry(iy):this->grid_2d.gy(iy);
					}
				}
			}

			// 1:(image_tot, image_coh); 2:(image_tot); 3:(m2psi_tot, m2psi_coh); 4:(m2psi_tot); 
			// 5:(m2psi_tot, psi_coh); 6:(psi_coh); 7:(psi_0); 8:(V); 9:(trans)
			void set_output_type()
			{
				if(this->is_STEM() || this->is_EELS())
				{
					output_type = (this->pn_coh_contrib)?eTEMOT_image_tot_coh:eTEMOT_image_tot;
				}
				else if(this->is_ISTEM() || this->is_CBED_CBEI() ||this->is_PED_HCTEM() || 
				this->is_ED_HRTEM() || this->is_EFTEM())
				{
					output_type = (this->pn_coh_contrib)?eTEMOT_m2psi_tot_coh:eTEMOT_m2psi_tot;
				}
				else if(this->is_EWFS_EWRS())
				{
					output_type = (this->is_EWFS_EWRS_SC())?eTEMOT_psi_coh:eTEMOT_m2psi_tot_psi_coh;
				}
				else if(this->is_PropFS_PropRS())
				{
					output_type = eTEMOT_psi_coh;
				}
				else if(this->is_IWFS_IWRS())
				{
					output_type = eTEMOT_psi_0;
				}
				else if(this->is_PPFS_PPRS())
				{
					output_type = eTEMOT_V;
				}
				else if(this->is_TFFS_TFRS())
				{
					output_type = eTEMOT_trans;
				}
			}

			template <class TInput_Multislice>
			void assign_input_multislice(TInput_Multislice &input_multislice)
			{ 
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
