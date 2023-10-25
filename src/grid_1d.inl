/*
 * This file is part of Multem.
 * Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>
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

#include "grid_1d.h"

/* template specialization 1d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T, class ST>
	CGPU_EXEC
	Grid_sxd<T, ST, edim_1>::Grid_sxd(): iGrid_sxd<ST, edim_1>(), bs_x(0), rx_0(0), pbc_x(true), 
		bwl(false), sli_thick(0), drx(0), dgx(0) {}

	template <class T, class ST>
	Grid_sxd<T, ST, edim_1>::Grid_sxd(const ST& nx)
	{
		set_in_data(nx);
	}

	template <class T, class ST>
	template <class U, class SU>
	Grid_sxd<T, ST, edim_1>::Grid_sxd(const U& bs_x, const SU& nx)
	{
		set_in_data(bs_x, nx);
	}

	template <class T, class ST>
	template <class U, class SU>
	Grid_sxd<T, ST, edim_1>::Grid_sxd(const U& bs_x, const SU& nx, const U& rx_0, 
	dt_bool pbc_x, dt_bool bwl, U sli_thick)
	{
		set_in_data(bs_x, nx, rx_0, pbc_x, bwl, sli_thick);
	}

	/* copy constructor */
	template <class T, class ST>
	CGPU_EXEC
	Grid_sxd<T, ST, edim_1>::Grid_sxd(const Grid_sxd<T, ST, edim_1>& grid)
	{
		*this = grid;
	}

	/* converting constructor */
	template <class T, class ST>
	template <class U, class SU>
	CGPU_EXEC
	Grid_sxd<T, ST, edim_1>::Grid_sxd(const Grid_sxd<U, SU, edim_1>& grid)
	{
		*this = grid;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, class ST>
	CGPU_EXEC
	Grid_sxd<T, ST, edim_1>& Grid_sxd<T, ST, edim_1>::operator=(const Grid_sxd<T, ST, edim_1>& grid)
	{
		if (this != &grid)
		{
			iGrid_sxd<ST, edim_1>::operator=(grid);

			bs_x = grid.bs_x;
			rx_0 = grid.rx_0;
			pbc_x = grid.pbc_x;
			bwl = grid.bwl;
			sli_thick = grid.sli_thick;

			drx = grid.drx;
			dgx = grid.dgx;
		}

		return *this;
	}

	/* converting assignment operator */
	template <class T, class ST>
	template <class U, class SU>
	CGPU_EXEC
	Grid_sxd<T, ST, edim_1>& Grid_sxd<T, ST, edim_1>::operator=(const Grid_sxd<U, SU, edim_1>& grid)
	{
		iGrid_sxd<ST, edim_1>::operator=(grid);

		bs_x = T(grid.bs_x);
		rx_0 = T(grid.rx_0);
		pbc_x = grid.pbc_x;
		bwl = grid.bwl;
		sli_thick = T(grid.sli_thick);

		drx = T(grid.drx);
		dgx = T(grid.dgx);

		return *this;
	}

	template <class T, class ST>
	template <class U, class SU> 
	CGPU_EXEC
	void Grid_sxd<T, ST, edim_1>::assign(const Grid_sxd<U, SU, edim_1>& grid)
	{
		*this = grid;
	}

	/************************ user define conversion operators ****************************/
	template <class T, class ST>
	Grid_sxd<T, ST, edim_1>::operator iRegion_Rect_xd<edim_1>() const
	{
		return {ST(0), this->nx};
	}

	/***************************************************************************************/
	template <class T, class ST>
	void Grid_sxd<T, ST, edim_1>::set_in_data(const ST& nx)
	{
		set_in_data(T(nx), nx, T(0));
	}

	template <class T, class ST>
	template <class U, class SU>
	void Grid_sxd<T, ST, edim_1>::set_in_data(const U& bs_x, const SU& nx)
	{
		set_in_data(T(bs_x), ST(nx), T(0));
	}

	template <class T, class ST>
	template <class U, class SU>
	void Grid_sxd<T, ST, edim_1>::set_in_data(const U& bs_x, const SU& nx, 
	const U& rx_0, dt_bool pbc_x, dt_bool bwl, U sli_thick)
	{
		this->set_size(nx);

		this->bs_x = T(bs_x);
		this->rx_0 = T(rx_0);
		this->pbc_x = pbc_x;
		this->bwl = bwl;
		this->sli_thick = T(sli_thick);

		set_dep_var();
	}

	template <class T, class ST>
	void Grid_sxd<T, ST, edim_1>::set_dep_var()
	{
		drx = mt::fcn_div(bs_x, this->nx);
		dgx = mt::fcn_div(T(1), bs_x);
	}

	template <class T, class ST>
	CGPU_EXEC
	void Grid_sxd<T, ST, edim_1>::set_r_0(const T& rx_0) 
	{	
		this->rx_0 = rx_0;
	}

	/***************************************************************************************/
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::nx_r() const
	{ 
		return T(this->nx);
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::size_r() const 
	{ 
		return T(this->size());
	}

	template <class T, class ST>
	T Grid_sxd<T, ST, edim_1>::isize_r() const
	{ 
		return T(1)/size_r();
	}

	/***************************************************************************************/
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::bs_x_h() const 
	{ 
		return T(0.5)*bs_x;
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::bs_h() const 
	{ 
		return bs_x_h();
	}

	/***************************************************************************************/
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::rx_c() const 
	{ 
		return rx_0 + bs_x_h();
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::rv_c() const 
	{ 
		return rx_c();
	}

	/***************************************************************************************/
	// maximum frequency
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::g_max() const 
	{ 
		return gx_back();
	}

	// maximum square frequency
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::g2_max() const 
	{ 
		return ::square(g_max());
	}

	// maximum allowed frequency
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::gl_max() const
	{
		return g_max()*T(2.0/3.0);
	}

	// maximum square allowed frequency
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::gl2_max() const
	{
		return ::square(gl_max());
	}

	/***************************************************************************************/
	/****************************** Fourier space positions ********************************/
	/***************************************************************************************/
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::gx(const ST& ix) const 
	{ 
		return T(this->igx(ix))*dgx;
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::gx2(const ST& ix) const 
	{ 
		return ::square(gx(ix));
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::g2(const ST& ix) const 
	{ 
		return gx2(ix);
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::g(const ST& ix) const 
	{ 
		return ::fabs(gx(ix));
	}

	/***************************************************************************************/
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::gx(const ST& ix, const T& gx_0) const 
	{ 
		return gx(ix) - gx_0;
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::gx2(const ST& ix, const T& gx_0) const 
	{ 
		return ::square(gx(ix, gx_0));
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::g2(const ST& ix, const T& gx_0) const 
	{ 
		return gx2(ix, gx_0);
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::g(const ST& ix, const T& gx_0) const 
	{ 
		return ::fabs(gx(ix, gx_0));
	}

	/************************************* shift *******************************************/
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::gx_sft(const ST& ix) const 
	{ 
		return T(this->igx_sft(ix))*dgx;
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::gx2_sft(const ST& ix) const 
	{ 
		return ::square(gx_sft(ix));
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::g2_sft(const ST& ix) const 
	{ 
		return gx2_sft(ix);
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::g_sft(const ST& ix) const 
	{ 
		return ::fabs(gx_sft(ix));
	}

	/***************************************************************************************/
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::gx_sft(const ST& ix, const T& gx_0) const 
	{ 
		return gx_sft(ix) - gx_0;
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::gx2_sft(const ST& ix, const T& gx_0) const 
	{ 
		return ::square(gx_sft(ix, gx_0));
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::g2_sft(const ST& ix, const T& gx_0) const 
	{ 
		return gx2_sft(ix, gx_0);
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::g_sft(const ST& ix, const T& gx_0) const 
	{ 
		return ::fabs(gx_sft(ix, gx_0));
	}

	/***************************************************************************************/
	/******************************** real space positions *********************************/
	/***************************************************************************************/
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::rx(const ST& ix) const 
	{ 
		return T(this->irx(ix))*drx + rx_0;
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::rx2(const ST& ix) const 
	{ 
		return ::square(rx(ix));
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::r2(const ST& ix) const 
	{ 
		return rx2(ix);
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::r(const ST& ix) const 
	{ 
		return ::fabs(rx(ix));
	}

	/***************************************************************************************/
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::rx(const ST& ix, const T& x0) const 
	{ 
		return rx(ix) - x0;
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::rx2(const ST& ix, const T& x0) const 
	{ 
		return ::square(rx(ix, x0));
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::r2(const ST& ix, const T& x0) const 
	{ 
		return rx2(ix, x0);
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::r(const ST& ix, const T& x0) const 
	{ 
		return ::fabs(rx(ix, x0));
	}

	/***************************************************************************************/
	/***************************************************************************************/
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::rx_sft(const ST& ix) const 
	{ 
		return T(this->irx_sft(ix))*drx + rx_0;
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::rx2_sft(const ST& ix) const 
	{ 
		return ::square(rx_sft(ix));
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::r2_sft(const ST& ix) const 
	{ 
		return rx2_sft(ix);
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::r_sft(const ST& ix) const 
	{ 
		return ::fabs(rx_sft(ix));
	}

	/***************************************************************************************/
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::rx_sft(const ST& ix, const T& x0) const 
	{ 
		return rx_sft(ix) - x0;
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::rx2_sft(const ST& ix, const T& x0) const 
	{ 
		return ::square(rx_sft(ix, x0));
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::r2_sft(const ST& ix, const T& x0) const 
	{ 
		return rx2_sft(ix, x0);
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::r_sft(const ST& ix, const T& x0) const 
	{ 
		return ::fabs(rx_sft(ix, x0));
	}

	/***************************************************************************************/
	/******************************* from position to index ********************************/
	/***************************************************************************************/
	template <class T, class ST>
	template <class SU>
	void Grid_sxd<T, ST, edim_1>::ix_0_ix_n(const T& x, const T& x_max, SU& ix_0, SU& ix_n) const 
	{
		fcn_get_idx_0_idx_n(x, x_max, drx, pbc_x, this->nx-1, ix_0, ix_n);
	}	
			
	template <class T, class ST>
	template <class SU>
	void Grid_sxd<T, ST, edim_1>::ix_0_ix_e(const T& x, const T& x_max, SU& ix_0, SU& ix_e) const 
	{
		fcn_get_idx_0_idx_n(x, x_max, drx, pbc_x, this->nx-1, ix_0, ix_e);
		ix_e += ix_0;
	}

	/************ fds = floor/division by pixel size ************/
	// locate x -> irx
	template <class T, class ST>
	CGPU_EXEC
	ST Grid_sxd<T, ST, edim_1>::rx_2_irx_fd(const T& x) const 
	{ 
		return fcn_cfloor<ST>(x/drx);
	}

	/********* bfds = bound/floor/division by pixel size ********/
	// locate x -> irx
	template <class T, class ST>
	CGPU_EXEC
	ST Grid_sxd<T, ST, edim_1>::rx_2_irx_bfd(const T& x) const 
	{ 
		return fcn_set_bound(rx_2_irx_fd(x), ST(0), this->nx-1);
	}

	/********* cds = ceil/division by pixel size/shift **********/
	// locate x -> irx
	template <class T, class ST>
	CGPU_EXEC
	ST Grid_sxd<T, ST, edim_1>::rx_2_irx_cd(const T& x) const 
	{ 
		return fcn_cceil<ST>(x/drx);
	}

	/****** bcds = bound/ceil/division by pixel size/shift ******/
	// locate x -> irx
	template <class T, class ST>
	CGPU_EXEC
	ST Grid_sxd<T, ST, edim_1>::rx_2_irx_bcd(const T& x) const 
	{ 
		return fcn_set_bound(rx_2_irx_cd(x), ST(0), this->nx-1);
	}

	/********* fds = floor/division by pixel size/shift *********/
	// locate x -> irx
	template <class T, class ST>
	CGPU_EXEC
	ST Grid_sxd<T, ST, edim_1>::rx_2_irx_fds(const T& x) const 
	{ 
		return fcn_cfloor<ST>((x - rx_0)/drx);
	}

	/****** bfds = bound/floor/division by pixel size/shift ******/
	// locate x -> irx
	template <class T, class ST>
	CGPU_EXEC
	ST Grid_sxd<T, ST, edim_1>::rx_2_irx_bfds(const T& x) const 
	{ 
		return fcn_set_bound(rx_2_irx_fds(x), ST(0), this->nx-1);
	}

	/********* cds = ceil/division by pixel size/shift **********/
	// locate x -> irx
	template <class T, class ST>
	CGPU_EXEC
	ST Grid_sxd<T, ST, edim_1>::rx_2_irx_cds(const T& x) const 
	{ 
		return fcn_cceil<ST>((x - rx_0)/drx);
	}

	/****** bcds = bound/ceil/division by pixel size/shift *******/
	// locate x -> irx
	template <class T, class ST>
	CGPU_EXEC
	ST Grid_sxd<T, ST, edim_1>::rx_2_irx_bcds(const T& x) const 
	{ 
		return fcn_set_bound(rx_2_irx_cds(x), ST(0), this->nx-1);
	}

	/************ From position to index by searching ***********/
	// locate x -> irx
	template <class T, class ST>
	CGPU_EXEC
	ST Grid_sxd<T, ST, edim_1>::rx_2_irx_b(const T& x, ST ix_min, ST ix_max) const 
	{
		if (ix_min == ix_max)
		{
			ix_min = ST(0);
			ix_max = this->nx-1;
		}

		return fcn_r_2_ir_b_by_fcn(rx, ix_min, ix_max);
	}

	/***************************************************************************************/
	/********************************** check bounds ***************************************/
	/***************************************************************************************/
	template <class T, class ST>
	CGPU_EXEC
	dt_bool Grid_sxd<T, ST, edim_1>::chk_bound_x(const T& x) const 
	{ 
		return fcn_chk_bound(x, rx_front(), rx_back());
	}

	template <class T, class ST>
	CGPU_EXEC
	dt_bool Grid_sxd<T, ST, edim_1>::chk_bound_x_eps(const T& x) const 
	{ 
		return fcn_chk_bound_eps(x, rx_front(), rx_back());
	}

	/***************************************************************************************/
	/************************************ set bounds ***************************************/
	/***************************************************************************************/
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::set_bound_x(const T& x) const 
	{ 
		return fcn_set_bound(x, rx_front(), rx_back());
	}

	/***************************************************************************************/
	/*********************************** front/back ****************************************/
	/***************************************************************************************/
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::gx_front() const 
	{ 
		return gx(ST(0));
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::gx_back() const 
	{ 
		return gx(this->nx-1);
	}

	/***************************************************************************************/
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::rx_front() const
	{ 
		return rx(ST(0));
	}

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::rx_back() const
	{ 
		return rx(this->nx-1);
	}

	/***************************************************************************************/
	/************************************** factors ****************************************/
	/***************************************************************************************/
	/* calculate fermi low-pass filter alpha parameter */
	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::fermi_lpf_alpha() const
	{ 
		return fcn_fermi_lpf_alpha(gl_max(), T(0.25), T(1e-02));
	}	

	template <class T, class ST>
	CGPU_EXEC
	T Grid_sxd<T, ST, edim_1>::factor_2pi_rx_ctr(const T& x) const 
	{ 
		return fcn_n_2pi_sft(x, bs_x_h());
	}

	template <class T, class ST>
	CPU_EXEC
	Vctr<T, edev_cpu> Grid_sxd<T, ST, edim_1>::factor_2pi_rv_ctr(const Vctr<T, edev_cpu>& rv) const
	{
		Vctr<T, edev_cpu> rv_o(rv.size());

		for(auto ik = 0; ik<rv.size(); ik++)
		{
			rv_o[ik] = factor_2pi_rx_ctr(rv[ik]);
		}

		return rv_o;
	}

	/***************************************************************************************/
	template <class T, class ST>
	iRegion_Rect_xd<edim_1> Grid_sxd<T, ST, edim_1>::iregion_rect(const T& r, const T& radius) const 
	{ 
		return {rx_2_irx_bfds(r - radius), rx_2_irx_bcds(r + radius)};
	}
}

/* traits */
namespace mt
{
	template <class T>
	struct is_grid_1d: std::integral_constant<dt_bool, std::is_same<T, Grid_1d_st<typename T::value_type, typename T::size_type>>::value> {};	

	/***************************************************************************************/
	template <class T, class U>
	struct is_grid_1d_and_vctr_cpu: std::integral_constant<dt_bool, is_grid_1d<T>::value && is_vctr_cpu<U>::value> {};

	template <class T, class U>
	struct is_grid_1d_and_vctr_gpu: std::integral_constant<dt_bool, is_grid_1d<T>::value && is_vctr_gpu<U>::value> {};

	/***************************************************************************************/
	template <class T, class U>
	struct is_grid_1d_and_cvctr_cpu: std::integral_constant<dt_bool, is_grid_1d<T>::value && is_cvctr_cpu<U>::value> {};

	template <class T, class U>
	struct is_grid_1d_and_cvctr_gpu: std::integral_constant<dt_bool, is_grid_1d<T>::value && is_cvctr_gpu<U>::value> {};

	/***************************************************************************************/
	template <class T, class U, class V=void>
	using enable_if_grid_1d_and_vctr_cpu = typename std::enable_if<is_grid_1d_and_vctr_cpu<T, U>::value, V>::type;

	template <class T, class U, class V=void>
	using enable_if_grid_1d_and_vctr_gpu = typename std::enable_if<is_grid_1d_and_vctr_gpu<T, U>::value, V>::type;

	/***************************************************************************************/
	template <class T, class U, class V=void>
	using enable_if_grid_1d_and_cvctr_cpu = typename std::enable_if<is_grid_1d_and_cvctr_cpu<T, U>::value, V>::type;

	template <class T, class U, class V=void>
	using enable_if_grid_1d_and_cvctr_gpu = typename std::enable_if<is_grid_1d_and_cvctr_gpu<T, U>::value, V>::type;
}