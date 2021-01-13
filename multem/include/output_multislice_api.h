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

#ifndef OUTPUT_MULTISLICE_API_H
#define OUTPUT_MULTISLICE_API_H

#include <vector>
#include <algorithm>

#include "math.cuh"
#include "safe_types.cuh"
#include "input_multislice_api.h"

namespace mt
{
	inline
		bool is_ot_image_tot_coh(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_image_tot_coh;
	}

	inline
		bool is_ot_image_tot(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_image_tot;
	}

	inline
		bool is_ot_m2psi_tot_coh(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_m2psi_tot_coh;
	}

	inline
		bool is_ot_m2psi_tot(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_m2psi_tot;
	}

	inline
		bool is_ot_m2psi_tot_psi_coh(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_m2psi_tot_psi_coh;
	}

	inline
		bool is_ot_psi_coh(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_psi_coh;
	}

	inline
		bool is_ot_psi_0(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_psi_0;
	}

	inline
		bool is_ot_V(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_V;
	}

	inline
		bool is_ot_trans(const eTEM_Output_Type &output_type)
	{
		return output_type == mt::eTEMOT_trans;
	}


	/**************************************************************************************/
	template <class T>
	class Output_Multislice : public Input_Multislice<T>
	{
	public:
		using T_r = T;
		using T_c = complex<T>;

		DLL_PUBLIC Output_Multislice();
    DLL_PUBLIC Output_Multislice(const Output_Multislice &other);
    DLL_PUBLIC Output_Multislice(Output_Multislice &&other);
    DLL_PUBLIC ~Output_Multislice();

    DLL_PUBLIC eTEM_Output_Type& output_type();
    DLL_PUBLIC const eTEM_Output_Type& output_type() const;

		DLL_PUBLIC int& ndetector();
		DLL_PUBLIC const int& ndetector() const;

    DLL_PUBLIC int& nx();
		DLL_PUBLIC const int& nx() const;

    DLL_PUBLIC int& ny();
		DLL_PUBLIC const int& ny() const;

    DLL_PUBLIC T_r& dx();
		DLL_PUBLIC const T_r& dx() const;

    DLL_PUBLIC T_r& dy();
		DLL_PUBLIC const T_r& dy() const;

		DLL_PUBLIC T_r& dr();
		DLL_PUBLIC const T_r& dr() const;

    DLL_PUBLIC std::vector<T_r> get_x() const;
    DLL_PUBLIC void set_x(const std::vector<T_r>& x);

    DLL_PUBLIC std::vector<T_r> get_y() const;
    DLL_PUBLIC void set_y(const std::vector<T_r>& y);

    DLL_PUBLIC std::vector<T_r> get_r() const;
    DLL_PUBLIC void set_r(const std::vector<T_r>& r);

    DLL_PUBLIC std::vector<std::vector<T_r>> get_m2psi_tot() const;
    DLL_PUBLIC void set_m2psi_tot(const std::vector<std::vector<T_r>>& m2psi_totx);

    DLL_PUBLIC std::vector<std::vector<T_r>> get_m2psi_coh() const;
    DLL_PUBLIC void set_m2psi_coh(const std::vector<std::vector<T_r>>& m2psi_coh);

    DLL_PUBLIC std::vector<std::vector<T_c>> get_psi_coh() const;
    DLL_PUBLIC void set_psi_coh(const std::vector<std::vector<T_c>>& psi_coh);

    DLL_PUBLIC std::vector<std::vector<T_r>> get_V() const;
    DLL_PUBLIC void set_V(const std::vector<std::vector<T_r>>& V);

    DLL_PUBLIC std::vector<std::vector<T_c>> get_trans() const;
    DLL_PUBLIC void set_trans(const std::vector<std::vector<T_c>>& trans);

    DLL_PUBLIC std::vector<std::vector<T_c>> get_psi_0() const;
    DLL_PUBLIC void set_psi_0(const std::vector<std::vector<T_c>>& psi_0);

		template <class TOutput_Multislice>
		DLL_PUBLIC void assign(TOutput_Multislice &output_multislice);

		template <class TOutput_Multislice>
		DLL_PUBLIC Output_Multislice<T>& operator=(TOutput_Multislice &output_multislice);

		DLL_PUBLIC void clear();
		DLL_PUBLIC void clean_temporal();

    template <class TInput_Multislice>
    DLL_PUBLIC void set_input_data(TInput_Multislice *input_multislice);

		DLL_PUBLIC void init();
		DLL_PUBLIC void init_psi_coh();

		/***************************************************************************/
		template<class TVector>
		DLL_PUBLIC void add_scale_psi_coh(int ithk, T_c w, TVector &phi);

		template<class TVector>
		DLL_PUBLIC void add_scale_shift_psi_coh(int ithk, T_c w, TVector &phi);

 		template<class TVector>
		DLL_PUBLIC void add_scale_crop_shift_psi_coh(int ithk, T_c w, TVector &phi);

   	template<class TVector>
		DLL_PUBLIC void add_scale_crop_shift_m2psi_coh_from_psi(int ithk, T_r w, TVector &phi);

 		template<class TVector>
		DLL_PUBLIC void add_scale_crop_shift_m2psi_coh_from_m2psi(int ithk, T_r w, TVector &m2phi);

 		template<class TVector>
		DLL_PUBLIC void add_scale_crop_shift_m2psi_tot_from_psi(int ithk, T_r w, TVector &phi);

 		template<class TVector>
		DLL_PUBLIC void add_scale_crop_shift_m2psi_tot_from_m2psi(int ithk, T_r w, TVector &m2phi);

		/***************************************************************************/
    template<class TVector>
		DLL_PUBLIC void set_psi_coh(int ithk, TVector &phi);

 		template<class TVector>
		DLL_PUBLIC void set_shift_psi_coh(int ithk, TVector &phi);

    template<class TVector>
		DLL_PUBLIC void set_crop_shift_psi_coh(int ithk, TVector &phi);

 		template<class TVector>
		DLL_PUBLIC void set_crop_shift_m2psi_coh(int ithk, TVector &m2phi);

 		template<class TVector>
		DLL_PUBLIC void set_crop_shift_m2psi_tot(int ithk, TVector &m2phi);

 		/***************************************************************************/
 		template<class TVector>
		DLL_PUBLIC void from_psi_coh_2_phi(int ithk, TVector &phi);

 		template<class TVector>
		DLL_PUBLIC void from_m2psi_coh_2_m2phi(int ithk, TVector &m2phi);

 		template<class TVector>
		DLL_PUBLIC void from_m2psi_tot_2_m2phi(int ithk, TVector &m2phi);

		DLL_PUBLIC void gather();
  
    DLL_PUBLIC std::vector<T> extract_data(ePhonon_Model_Output fp_ctr, eShow_CData show_data, int ithk, int idet = 0);

		/**************************************************************************************/
		inline
		bool is_ot_image_tot_coh() const
		{
			return mt::is_ot_image_tot_coh(output_type());
		}

		inline
		bool is_ot_image_tot() const
		{
			return mt::is_ot_image_tot(output_type());
		}

		inline
		bool is_ot_m2psi_tot_coh() const
		{
			return mt::is_ot_m2psi_tot_coh(output_type());
		}

		inline
		bool is_ot_m2psi_tot() const
		{
			return mt::is_ot_m2psi_tot(output_type());
		}

		inline
		bool is_ot_m2psi_tot_psi_coh() const
		{
			return mt::is_ot_m2psi_tot_psi_coh(output_type());
		}

		inline
		bool is_ot_psi_coh() const
		{
			return mt::is_ot_psi_coh(output_type());
		}

		inline
		bool is_ot_psi_0() const
		{
			return mt::is_ot_psi_0(output_type());
		}

		inline
		bool is_ot_V() const
		{
			return mt::is_ot_V(output_type());
		}

		inline
		bool is_ot_trans() const
		{
			return mt::is_ot_trans(output_type());
		}

    struct Implementation;

    Implementation& data();
    const Implementation& data() const;

	private:

    std::unique_ptr<Implementation> pimpl;

    void set_output_grid();

		void set_output_type();
  
    template <class TInput_Multislice>
    void assign_input_multislice(TInput_Multislice &input_multislice);

	};

} // namespace mt

#endif
