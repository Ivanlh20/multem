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

#ifndef CGPU_CLASSES_H
	#define CGPU_CLASSES_H

	#include "math_mt.h"
	#include "types.cuh"
	#include "type_traits_gen.h"
	#include "cgpu_stream.cuh"
	#include "cgpu_fft.cuh"
	//#include "eval_fit_gaussians.hpp"
	#include "cgpu_rand.cuh"

	#ifdef __CUDACC__
		#include <cuda.h>
		#include <cuda_runtime.h>
		#include <cufft.h>
	#endif

	#include "in_classes.cuh"
	#include "fcns_gpu.h"
	#include "fcns_cpu.h"
	#include "particles.cuh"

	#ifdef __CUDACC__
		#include "fcns_gpu.h"
	#endif

	namespace mt
	{	
		/********************************* Gaussian Conv ***************************************/
		template <class T, eDev Dev>
		class Gauss_Cv_1d
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				static const eDev device = Dev;

				Gauss_Cv_1d():fft_1d(nullptr) {}

				Gauss_Cv_1d(FFT<T, Dev> *fft_1d_i, Grid_1d<T>& grid_1d_i)
				{
					set_in_data(fft_1d_i, grid_1d_i);
				}

				inline
				void set_in_data(FFT<T, Dev> *fft_1d_i, Grid_1d<T>& grid_1d_i)
				{
					fft_1d = fft_1d_i;
					grid_1d = grid_1d_i;
				}

				void set_fft_plan()
				{
					fft_1d->create_plan_1d(grid_1d.nx, 1);
				}

				void operator()(T sigma_r, TVctr_c& Im)
				{
					fcn_fftsft_1d(grid_1d, Im);
					fft_1d->forward(Im);

					gauss_cv_1d(sigma_r, Im);

					fft_1d->inverse(Im);
					fcn_fftsft_1d(grid_1d, Im);
				}

				void operator()(T sigma_r, TVctr_r &Im)
				{
					TVctr_c Im_c(Im.begin(), Im.end());

					this->operator()(sigma_r, Im_c);

					fcn_assign_real(Im_c, Im);
				}

				void cleanup()
				{
					fft_1d->cleanup();
				}
		protected:
				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				gauss_cv_1d(T sigma_r, TVctr_c& Im)
				{
					auto alpha = 2*c_pi2*sigma_r*sigma_r;

					for(auto ix = 0; ix < grid_1d.nx; ix++)
					{
						cgpu_detail::gauss_cv_1d<Grid_1d<T>, TVctr_c>(ix, grid_1d, alpha, Im);
					}
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				gauss_cv_1d(T sigma_r, TVctr_c& Im)
				{
					auto alpha = 2*c_pi2*sigma_r*sigma_r;

					auto d_grid_blk = grid_1d.d_grid_blk();
					gpu_detail::gauss_cv_1d<Grid_1d<T>, typename TVctr_c::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_1d, alpha, Im);
				}
			#endif

				FFT<T, Dev> *fft_1d;
				Grid_1d<T> grid_1d;
		};

		template <class T, eDev Dev>
		class Gauss_Cv_2d_BC
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				static const eDev device = Dev;

				Gauss_Cv_2d_BC():stream(nullptr), fft_2d(nullptr) {}

				Gauss_Cv_2d_BC(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					set_in_data(stream_i, fft_2d_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					stream = stream_i;
					fft_2d = fft_2d_i;
					grid_2d = grid_2d_i;
				}

				void set_fft_plan()
				{
					fft_2d->create_plan_1d_batch(grid_2d.ny, grid_2d.nx, stream->size());
				}

				void operator()(T sigma_r, TVctr_c& Im)
				{
					fcn_fftsft_bc_2d(*stream, grid_2d, Im);
					fft_2d->forward(Im);

					gauss_cv_2d_bc(sigma_r, Im);

					fft_2d->inverse(Im);
					fcn_fftsft_bc_2d(*stream, grid_2d, Im);
				}

				void operator()(T sigma_r, TVctr_r &Im)
				{
					TVctr_c Im_c(Im.begin(), Im.end());

					this->operator()(sigma_r, Im_c);

					fcn_assign_real(Im_c, Im);
				}

				void cleanup()
				{
					fft_2d->cleanup();
				}

		protected:
				Vctr<T, edev_cpu> gauss_vector_1d(T alpha)
				{
					TVctr_r fg;
					fg.reserve(grid_2d.ny);

					for(auto iy = 0; iy < grid_2d.ny; iy++)
					{
						auto v = exp(-alpha*grid_2d.gy2_sft(iy))/grid_2d.ny_r();
						fg.push_back(v);
					}
					return fg;
				}

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				gauss_cv_2d_bc(T sigma_r, TVctr_c& M_g)
				{
					auto alpha = 2*c_pi2*sigma_r*sigma_r;
					auto fg = gauss_vector_1d(alpha);

					stream->set_n_stream_act(grid_2d.nx);
					stream->set_grid(grid_2d.nx, grid_2d.ny);
					stream->exec_2d(cgpu_detail::fcn_ew_mult_mx_vctr_col<Grid_2d<T>, TVctr_r, TVctr_c>, grid_2d, fg, M_g);
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				gauss_cv_2d_bc(T sigma_r, TVctr_c& M_g)
				{
					auto alpha = 2*c_pi2*sigma_r*sigma_r;
					auto fg_h = gauss_vector_1d(alpha);
					TVctr_r fg = fg_h;

					auto d_grid_blk = grid_2d.d_grid_blk();
					gpu_detail::fcn_ew_mult_mx_vctr_col<Grid_2d<T>, typename TVctr_c::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, fg, M_g);
				}
			#endif

				Stream<Dev> *stream;
				FFT<T, Dev> *fft_2d;
				Grid_2d<T> grid_2d;
		};

		template <class T, eDev Dev>
		class Gauss_Cv_2d
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				static const eDev device = Dev;

				Gauss_Cv_2d():stream(nullptr), fft_2d(nullptr) {}

				Gauss_Cv_2d(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					set_in_data(stream_i, fft_2d_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					stream = stream_i;
					fft_2d = fft_2d_i;
					grid_2d = grid_2d_i;
				}

				void set_fft_plan()
				{
					fft_2d->create_plan_2d(grid_2d.ny, grid_2d.nx, stream->size());
				}

				void operator()(T sigma_r, TVctr_c& Im)
				{
					if (fcn_is_zero(sigma_r))
					{
						return;
					}

					fcn_fftsft_2d(*stream, grid_2d, Im);
					fft_2d->forward(Im);

					gauss_cv_2d(sigma_r, Im);

					fft_2d->inverse(Im);
					fcn_fftsft_2d(*stream, grid_2d, Im);
				}

				void operator()(T sigma_r, TVctr_r &Im)
				{
					if (fcn_is_zero(sigma_r))
					{
						return;
					}

					Im_c.assign(Im.begin(), Im.end());

					this->operator()(sigma_r, Im_c);

					fcn_assign_real(Im_c, Im);
				}

				void cleanup()
				{
					fft_2d->cleanup();
				}
		protected:
				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				gauss_cv_2d(T sigma_r, TVctr_c& Im)
				{
					auto alpha = 2*c_pi2*sigma_r*sigma_r;

					stream->set_n_stream_act(grid_2d.nx);
					stream->set_grid(grid_2d.nx, grid_2d.ny);
					stream->exec_2d(cgpu_detail::gauss_cv_2d<Grid_2d<T>, TVctr_c>, grid_2d, alpha, Im);
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				gauss_cv_2d(T sigma_r, TVctr_c& Im)
				{
					auto alpha = 2*c_pi2*sigma_r*sigma_r;

					auto d_grid_blk = grid_2d.d_grid_blk();
					gpu_detail::gauss_cv_2d<Grid_2d<T>, typename TVctr_c::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, alpha, Im);
				}
			#endif

				Stream<Dev> *stream;
				FFT<T, Dev> *fft_2d;
				Grid_2d<T> grid_2d;

				TVctr_c Im_c;
		};

		/********************************* Gaussian Deconv *************************************/
		template <class T, eDev Dev>
		class Gauss_Dcv_1d
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				static const eDev device = Dev;

				Gauss_Dcv_1d():fft_1d(nullptr) {}

				Gauss_Dcv_1d(FFT<T, Dev> *fft_1d_i, Grid_1d<T>& grid_1d_i)
				{
					set_in_data(fft_1d_i, grid_1d_i);
				}

				inline
				void set_in_data(FFT<T, Dev> *fft_1d_i, Grid_1d<T>& grid_1d_i)
				{
					fft_1d = fft_1d_i;
					grid_1d = grid_1d_i;
				}

				void set_fft_plan()
				{
					fft_1d->create_plan_1d(grid_1d.nx, 1);
				}

				void operator()(T sigma_r, T PSNR, TVctr_c& Im)
				{
					fcn_fftsft_1d(grid_1d, Im);
					fft_1d->forward(Im);

					gauss_dcv_1d(sigma_r, PSNR, Im);

					fft_1d->inverse(Im);
					fcn_fftsft_1d(grid_1d, Im);
				}

				void operator()(T sigma_r, T PSNR, TVctr_r &Im)
				{
					TVctr_c Im_c(Im.begin(), Im.end());

					this->operator()(sigma_r, PSNR, Im_c);

					fcn_assign_real(Im_c, Im);
				}

				void cleanup()
				{
					fft_1d->cleanup();
				}
		protected:
				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				gauss_dcv_1d(T sigma_r, T PSNR, TVctr_c& Im)
				{
					auto alpha = 2*c_pi2*sigma_r*sigma_r;

					for(auto ix = 0; ix < grid_1d.nx; ix++)
					{
						cgpu_detail::gauss_dcv_1d<Grid_1d<T>, TVctr_c>(ix, grid_1d, alpha, PSNR, Im);
					}
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				gauss_dcv_1d(T sigma_r, T PSNR, TVctr_c& Im)
				{
					auto alpha = 2*c_pi2*sigma_r*sigma_r;

					auto d_grid_blk = grid_1d.d_grid_blk();
					gpu_detail::gauss_dcv_1d<Grid_1d<T>, typename TVctr_c::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_1d, alpha, PSNR, Im);
				}
			#endif

				FFT<T, Dev> *fft_1d;
				Grid_1d<T> grid_1d;
		};

		template <class T, eDev Dev>
		class Gauss_Dcv_2d_BC
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				static const eDev device = Dev;

				Gauss_Dcv_2d_BC():stream(nullptr), fft_2d(nullptr) {}

				Gauss_Dcv_2d_BC(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					set_in_data(stream_i, fft_2d_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					stream = stream_i;
					fft_2d = fft_2d_i;
					grid_2d = grid_2d_i;
				}

				void set_fft_plan()
				{
					fft_2d->create_plan_1d_batch(grid_2d.ny, grid_2d.nx, stream->size());
				}

				void operator()(T sigma_r, T PSNR, TVctr_c& Im)
				{
					fcn_fftsft_bc_2d(*stream, grid_2d, Im);
					fft_2d->forward(Im);

					gauss_cv_2d_bc(sigma_r, PSNR, Im);

					fft_2d->inverse(Im);
					fcn_fftsft_bc_2d(*stream, grid_2d, Im);
				}

				void operator()(T sigma_r, T PSNR, TVctr_r &Im)
				{
					TVctr_c Im_c(Im.begin(), Im.end());

					this->operator()(sigma_r, PSNR, Im_c);

					fcn_assign_real(Im_c, Im);
				}

				void cleanup()
				{
					fft_2d->cleanup();
				}

		protected:
				Vctr<T, edev_cpu> gauss_vector_1d(T alpha, T PSNR)
				{
					TVctr_r fg;
					fg.reserve(grid_2d.ny);

					for(auto iy = 0; iy < grid_2d.ny; iy++)
					{
						auto v = exp(-alpha*grid_2d.gy2_sft(iy));
						fg.push_back(v/((v*v+PSNR)*grid_2d.ny_r()));
					}
					return fg;
				}

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				gauss_cv_2d_bc(T sigma_r, T PSNR, TVctr_c& M_g)
				{
					auto alpha = 2*c_pi2*sigma_r*sigma_r;
					auto fg = gauss_vector_1d(alpha, PSNR);

					stream->set_n_stream_act(grid_2d.nx);
					stream->set_grid(grid_2d.nx, grid_2d.ny);
					stream->exec_2d(cgpu_detail::fcn_ew_mult_mx_vctr_col<Grid_2d<T>, TVctr_r, TVctr_c>, grid_2d, fg, M_g);
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				gauss_cv_2d_bc(T sigma_r, T PSNR, TVctr_c& M_g)
				{
					auto alpha = 2*c_pi2*sigma_r*sigma_r;
					auto fg_h = gauss_vector_1d(alpha, PSNR);
					TVctr_r fg = fg_h;

					auto d_grid_blk = grid_2d.d_grid_blk();
					gpu_detail::fcn_ew_mult_mx_vctr_col<Grid_2d<T>, typename TVctr_c::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, fg, M_g);
				}
			#endif

				Stream<Dev> *stream;
				FFT<T, Dev> *fft_2d;
				Grid_2d<T> grid_2d;
		};

		template <class T, eDev Dev>
		class Gauss_Dcv_2d
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				static const eDev device = Dev;

				Gauss_Dcv_2d():stream(nullptr), fft_2d(nullptr) {}

				Gauss_Dcv_2d(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					set_in_data(stream_i, fft_2d_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					stream = stream_i;
					fft_2d = fft_2d_i;
					grid_2d = grid_2d_i;
				}

				void set_fft_plan()
				{
					fft_2d->create_plan_2d(grid_2d.ny, grid_2d.nx, stream->size());
				}

				void operator()(T sigma_r, T PSNR, TVctr_c& Im)
				{
					fcn_fftsft_2d(*stream, grid_2d, Im);
					fft_2d->forward(Im);

					gauss_dcv_2d(sigma_r, PSNR, Im);

					fft_2d->inverse(Im);
					fcn_fftsft_2d(*stream, grid_2d, Im);
				}

				void operator()(T sigma_r, T PSNR, TVctr_r &Im)
				{
					TVctr_c Im_c(Im.begin(), Im.end());

					this->operator()(sigma_r, PSNR, Im_c);

					fcn_assign_real(Im_c, Im);
				}

				void cleanup()
				{
					fft_2d->cleanup();
				}
		protected:
				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				gauss_dcv_2d(T sigma_r, T PSNR, TVctr_c& Im)
				{
					auto alpha = 2*c_pi2*sigma_r*sigma_r;

					stream->set_n_stream_act(grid_2d.nx);
					stream->set_grid(grid_2d.nx, grid_2d.ny);
					stream->exec_2d(cgpu_detail::gauss_dcv_2d<Grid_2d<T>, TVctr_c>, grid_2d, alpha, PSNR, Im);
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				gauss_dcv_2d(T sigma_r, T PSNR, TVctr_c& Im)
				{
					auto alpha = 2*c_pi2*sigma_r*sigma_r;

					auto d_grid_blk = grid_2d.d_grid_blk();
					gpu_detail::gauss_dcv_2d<Grid_2d<T>, typename TVctr_c::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, alpha, PSNR, Im);
				}
			#endif

				Stream<Dev> *stream;
				FFT<T, Dev> *fft_2d;
				Grid_2d<T> grid_2d;
		};

		/******************************* affine transformations ********************************/	
		template <class T, eDev Dev>
		class Sd_2d
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Sd_2d(): stream(nullptr) {}

				Sd_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					set_in_data(stream_i, grid_2d_i, bg_opt_i, bg_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					stream = stream_i;
					grid_2d = grid_2d_i;
					bg_opt = bg_opt_i;
					bg = bg_i;
				}

				void generate_dx_dy(T ds_x, T phi_x, T ds_y, T phi_y, dt_int32 seed, TVctr& dx_o, TVctr& dy_o)
				{
					if (seed<1)
					{
						std::random_device rd;
						seed = rd();
					}

					gen.seed(seed);
					rnd_n.reset();

 					dx_o = generate_ds(ds_x, phi_x);
					dy_o = generate_ds(ds_y, phi_y);
				}

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, T>
				operator()(TVctr& mx_i, TVctr& dx, TVctr& dy, TVctr& mx_o)
				{
					// calculate background
					T bg = get_bg(mx_i);

					stream->set_n_stream_act(grid_2d.nx);
					stream->set_grid(grid_2d.nx, grid_2d.ny);
					stream->exec_2d(cgpu_detail::sd_2d<Grid_2d<T>, TVctr>, grid_2d, mx_i, dx, dy, bg, mx_o);

					return bg;
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, T>
				operator()(TVctr& mx_i, TVctr& dx, TVctr& dy, TVctr& mx_o)
				{
					// calculate background
					T bg = get_bg(mx_i);

					auto d_grid_blk = grid_2d.d_grid_blk();
					gpu_detail::sd_2d<Grid_2d<T>, typename TVctr::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, mx_i, dx, dy, bg, mx_o);

					return bg;
				}
			#endif

			protected:
				std::mt19937_64 gen;
				std::normal_distribution<T> rnd_n;

				Grid_2d<T> grid_2d;
				Stream<Dev> *stream;
				eFil_Sel_Typ bg_opt;
				T bg;

				Vctr<T, edev_cpu> generate_ds(T ds, T phi)
				{
					Vctr<T, edev_cpu> ds_v;
					ds_v.reserve(grid_2d.ny);

					T a = ds*rnd_n(gen);
					ds_v.push_back(a/::sqrt(1.0-phi*phi));
					for(auto iy = 1; iy < grid_2d.ny; iy++)
					{
			 T a = ds*rnd_n(gen);
			 ds_v.push_back(phi*ds_v.back() + a);
					}
					return ds_v;
				}

				T get_bg(TVctr& M)
				{
					T bg_r = 0;
					switch (bg_opt) 
					{
						case efst_min:
							bg_r = fcn_min_element(M);
						break;
						case efst_max:
							bg_r = fcn_max_element(M);
						break;
						case efst_mean:
							bg_r = fcn_mean(M);
						break;
						case efst_min_mean:
							bg_r = 0.5*(fcn_mean(M) + fcn_min_element(M));
						break;
						case efst_max_mean:
							bg_r = 0.5*(fcn_mean(M) + fcn_max_element(M));
						break;
						case efst_user_def:
							bg_r = bg;
						break;
					}

					return bg_r;
				}
		};

		template <class T, eDev Dev>
		class Sd_nr_2d
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Sd_nr_2d(): stream(nullptr) {}

				Sd_nr_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					set_in_data(stream_i, grid_2d_i, bg_opt_i, bg_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					stream = stream_i;
					grid_2d = grid_2d_i;
					bg_opt = bg_opt_i;
					bg = bg_i;
				}

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, T>
				operator()(TVctr& mx_i, TVctr& ds_x_i, TVctr& ds_y_i, TVctr& mx_o)
				{
					preprocessing(mx_i, ds_x_i, ds_y_i, parm);

					stream->set_n_stream_act(grid_2d.nx);
					stream->set_grid(grid_2d.nx, grid_2d.ny);
					stream->exec_2d(cgpu_detail::sd_nr_2d<Grid_2d<T>, TVctr, SPar>, grid_2d, mx_i, parm, mx_o);

					return parm.bg;
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, T>
				operator()(TVctr& mx_i, TVctr& ds_x_i, TVctr& ds_y_i, TVctr& mx_o)
				{
					// preprocessing(ds_x_i, ds_y_i, parm);
					// auto d_grid_blk = grid_2d.d_grid_blk();
					// gpu_detail::sd_nr_2d<Grid_2d<T>, typename TVctr::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, mx_i, ds_x_i, ds_y_i, mx_o);
				
					return parm.bg;
				}

			#endif

			protected:
				struct SPar
				{
					eFil_Sel_Typ bg_opt;
					T bg;
					TVctr dx;
					TVctr dy;
					TVctr ry_s;
					Vctr<dt_int32, Dev> iy;
				} parm;

				struct sort_by_first
				{
					template <class Ttuple1, class Ttuple2>
					CGPU_EXEC
					dt_bool operator()(const Ttuple1 &t1, const Ttuple2 &t2)
					{
						return thrust::get<0>(t1) < thrust::get<0>(t2);
					}
				};

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				preprocessing(TVctr& mx_i, TVctr& ds_x_i, TVctr& ds_y_i, SPar &parm)
				{
					// calculate background
					parm.bg_opt = bg_opt;
					parm.bg = get_bg(mx_i);

					parm.dx = ds_x_i;
					parm.dy = ds_y_i;
					parm.ry_s.resize(parm.dy.size());
					parm.iy.resize(parm.dy.size());

					for(auto ik = 0; ik < parm.dy.size(); ik++)
					{
						parm.ry_s[ik] = grid_2d.ry(ik) + parm.dy[ik];
						parm.iy[ik] = ik;
					}

					auto first = thrust::make_zip_iterator(thrust::make_tuple(parm.ry_s.begin(), parm.iy.begin()));
					auto last = thrust::make_zip_iterator(thrust::make_tuple(parm.ry_s.end(), parm.iy.end()));

					thrust::sort(first, last, sort_by_first());
				}

				/**********************Device**********************/
				// not working yet
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				preprocessing(TVctr& mx_i, TVctr& ds_x_i, TVctr& ds_y_i, SPar &parm)
				{
					// calculate background
					parm.bg_opt = bg_opt;
					parm.bg = get_bg(mx_i);

					parm.dx = ds_x_i;
					parm.dy = ds_y_i;
					parm.ry_s.resize(parm.dy.size());
					parm.iy.resize(parm.dy.size());

					for(auto ik = 0; ik < parm.dy.size(); ik++)
					{
						parm.ry_s[ik] = grid_2d.ry(ik) + parm.dy[ik];
						parm.iy[ik] = ik;
					}

					auto first = thrust::make_zip_iterator(thrust::make_tuple(parm.ry_s.begin(), parm.iy.begin()));
					auto last = thrust::make_zip_iterator(thrust::make_tuple(parm.ry_s.end(), parm.iy.end()));

					thrust::sort(first, last, sort_by_first());
				}
			#endif

			protected:
				Grid_2d<T> grid_2d;
				Stream<Dev> *stream;
				eFil_Sel_Typ bg_opt;
				T bg;

				T get_bg(TVctr& M)
				{
					T bg_r = 0;
					switch (bg_opt) 
					{
						case efst_min:
							bg_r = fcn_min_element(M);
						break;
						case efst_max:
							bg_r = fcn_max_element(M);
						break;
						case efst_mean:
							bg_r = fcn_mean(M);
						break;
						case efst_min_mean:
							bg_r = 0.5*(fcn_mean(M) + fcn_min_element(M));
						break;
						case efst_max_mean:
							bg_r = 0.5*(fcn_mean(M) + fcn_max_element(M));
						break;
						case efst_user_def:
							bg_r = bg;
						break;
					}

					return bg_r;
				}
		};

		template <class T, eDev Dev>
		class Data_Aug_2d
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Data_Aug_2d():stream(nullptr) {}

				Data_Aug_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i)
				{
					set_in_data(stream_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i)
				{
					stream = stream_i;
					grid_2d = grid_2d_i;
				}

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				operator()(TVctr& mx_i, T theta, R_2d<T> p0, T sx, T sy, R_2d<T> ps, T sim, 
				dt_int32 nx_o, dt_int32 ny_o, TVctr& mx_o)
				{
					// calculate maximum
					auto M_max = *thrust::max_element(mx_i.begin(), mx_i.end());

					// calculate background
					auto g_max = grid_2d.gx_back();
					auto g_min = 0.9*g_max;
					// T bg = fcn_int_det_ring(*stream, grid_2d, g_min, g_max, mx_i);

					T bg = 0.6*fcn_mean(*stream, mx_i);

					T rx_0 = p0.x*sx-(nx_o/2)*grid_2d.drx;
					T ry_0 = p0.y*sy-(ny_o/2)*grid_2d.dry;

					Grid_2d<T> grid_2d_o(nx_o, ny_o, nx_o*grid_2d.drx, ny_o*grid_2d.dry);
					grid_2d_o.set_r_0(rx_0, ry_0);

					stream->set_n_stream_act(grid_2d_o.nx);
					stream->set_grid(grid_2d_o.nx, grid_2d_o.ny);
					stream->exec_2d(cgpu_detail::rot_sca_sft_2d<Grid_2d<T>, TVctr>, grid_2d, mx_i, theta, p0, sx, sy, ps, bg, grid_2d_o, mx_o);

					normalized_data(mx_o, M_max);

					// add Poisson noise
					mx_o = fcn_add_poiss_nois(*stream, mx_o, sim, -1);

					M_max = *thrust::max_element(mx_o.begin(), mx_o.end());
					normalized_data(mx_o, M_max);
				}

			protected:
				void normalized_data(TVctr& M, T M_max)
				{
					std::for_each(M.begin(), M.end(), [M_max](T& v){ v = (v>0)?v/M_max:0; });
				}

				Grid_2d<T> grid_2d;
				Stream<Dev> *stream;
		};

		template <class T, eDev Dev>
		class Interp_rn_2d
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Interp_rn_2d(): stream(nullptr) {}

				Interp_rn_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, Grid_2d<T>& grid_2d_o, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					set_in_data(stream_i, grid_2d_i, grid_2d_o, bg_opt_i, bg_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, Grid_2d<T>& grid_2d_o, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					stream = stream_i;
					grid_2d = grid_2d_i;
					grid_2d_mo = grid_2d_o;
					bg_opt = bg_opt_i;
					bg = bg_i;
				}

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, T>
				operator()(TVctr& mx_i, TVctr& Rx_i, TVctr& Ry_i, TVctr& mx_o)
				{
					// calculate background
					T bg = get_bg(mx_i);

					stream->set_n_stream_act(grid_2d_mo.nx);
					stream->set_grid(grid_2d_mo.nx, grid_2d_mo.ny);
					stream->exec_2d(cgpu_detail::intrpl_rg_2d<Grid_2d<T>, TVctr>, grid_2d, mx_i, Rx_i, Ry_i, grid_2d_mo, bg, mx_o);
				
					return bg;
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, T>
				operator()(TVctr& mx_i, TVctr& Rx_i, TVctr& Ry_i, TVctr& mx_o)
				{
					// calculate background
					T bg = get_bg(mx_i);

					auto d_grid_blk = grid_2d_mo.d_grid_blk();
					gpu_detail::intrpl_rg_2d<Grid_2d<T>, typename TVctr::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, mx_i, Rx_i, Ry_i, grid_2d_mo, bg, mx_o);
				
					return bg;
				}
			#endif

			protected:
				Grid_2d<T> grid_2d;
				Grid_2d<T> grid_2d_mo;
				Stream<Dev> *stream;
				eFil_Sel_Typ bg_opt;
				T bg;

				T get_bg(TVctr& M)
				{
					T bg_r = 0;
					switch (bg_opt) 
					{
						case efst_min:
							bg_r = fcn_min_element(M);
						break;
						case efst_max:
							bg_r = fcn_max_element(M);
						break;
						case efst_mean:
							bg_r = fcn_mean(M);
						break;
						case efst_min_mean:
							bg_r = 0.5*(fcn_mean(M) + fcn_min_element(M));
						break;
						case efst_max_mean:
							bg_r = 0.5*(fcn_mean(M) + fcn_max_element(M));
						break;
						case efst_user_def:
							bg_r = bg;
						break;
					}

					return bg_r;
				}
		};

		template <class T, eDev Dev>
		class Rs_2d
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Rs_2d(): stream(nullptr) {}

				Rs_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i)
				{
					set_in_data(stream_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i)
				{
					stream = stream_i;
					grid_2d = grid_2d_i;
				}

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				operator()(TVctr& mx_i, T sxy, TVctr& mx_o)
				{
					if (fcn_is_equal<T>(sxy, T(1)))
					{
						mx_o = mx_i;
					}

					dt_int32 nx_o = fcn_sc_size_c(grid_2d.nx, sxy);
					dt_int32 ny_o = fcn_sc_size_c(grid_2d.ny, sxy);

					Grid_2d<T> grid_2d_o(nx_o, ny_o, nx_o*grid_2d.drx, ny_o*grid_2d.dry);

					stream->set_n_stream_act(grid_2d_o.nx);
					stream->set_grid(grid_2d_o.nx, grid_2d_o.ny);
					stream->exec_2d(cgpu_detail::sc_2d<Grid_2d<T>, TVctr>, grid_2d, mx_i, sxy, grid_2d_o, mx_o);
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				operator()(TVctr& mx_i, T sxy, TVctr& mx_o)
				{
					if (fcn_is_equal<T>(sxy, T(1)))
					{
						mx_o = mx_i;
					}

					dt_int32 nx_o = max(dt_int32(::floor(grid_2d.nx*sxy)), 1);
					dt_int32 ny_o = max(dt_int32(::floor(grid_2d.ny*sxy)), 1);

					Grid_2d<T> grid_2d_o(nx_o, ny_o, nx_o*grid_2d.drx, ny_o*grid_2d.dry);

					auto d_grid_blk = grid_2d_o.d_grid_blk();
					gpu_detail::sc_2d<Grid_2d<T>, typename TVctr::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, mx_i, sxy, grid_2d_o, mx_o);
				}
			#endif

			protected:
				Grid_2d<T> grid_2d;
				Stream<Dev> *stream;
		};

		template <class T, eDev Dev>
		class Rot_Sca_sft_2d
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Rot_Sca_sft_2d():stream(nullptr) {}

				Rot_Sca_sft_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i)
				{
					set_in_data(stream_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i)
				{
					stream = stream_i;
					grid_2d = grid_2d_i;
				}

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				operator()(TVctr& mx_i, T theta, R_2d<T> p0, T sx, T sy, R_2d<T> ps, TVctr& mx_o)
				{
					// calculate background
					T bg = fcn_mean(*stream, mx_i);

					dt_int32 nx_o = fcn_sc_size_c(grid_2d.nx, sx);
					dt_int32 ny_o = fcn_sc_size_c(grid_2d.ny, sy);

					Grid_2d<T> grid_2d_o(nx_o, ny_o, nx_o*grid_2d.drx, ny_o*grid_2d.dry);

					stream->set_n_stream_act(grid_2d_o.nx);
					stream->set_grid(grid_2d_o.nx, grid_2d_o.ny);
					stream->exec_2d(cgpu_detail::rot_sca_sft_2d<Grid_2d<T>, TVctr>, grid_2d, mx_i, theta, p0, sx, sy, ps, bg, grid_2d_o, mx_o);
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				operator()(TVctr& mx_i, T theta, R_2d<T> p0, T sx, T sy, R_2d<T> ps, TVctr& mx_o)
				{
					// calculate background
					T bg = fcn_mean(*stream, mx_i);

					dt_int32 nx_o = max(dt_int32(::floor(grid_2d.nx*sx)), 1);
					dt_int32 ny_o = max(dt_int32(::floor(grid_2d.ny*sy)), 1);

					Grid_2d<T> grid_2d_o(nx_o, ny_o, nx_o*grid_2d.drx, ny_o*grid_2d.dry);

					auto d_grid_blk = grid_2d_o.d_grid_blk();
					gpu_detail::rot_sca_sft_2d<Grid_2d<T>, typename TVctr::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, mx_i, theta, p0, sx, sy, ps, bg, grid_2d_o, mx_o);
				}
			#endif

			protected:
				Grid_2d<T> grid_2d;
				Stream<Dev> *stream;
		};

		/************************************ Gradient ******************************************/
		template <class T, eDev Dev>
		class Gradient
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Gradient(): stream(nullptr) {}

				Gradient(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i)
				{
					set_in_data(stream_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i)
				{
					stream = stream_i;
					grid_2d = grid_2d_i;
				}

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				operator()(TVctr& mx_i, TVctr& dM_x, TVctr& dM_y)
				{
					stream->set_n_stream_act(grid_2d.nx);
					stream->set_grid(grid_2d.nx, grid_2d.ny);
					stream->exec_2d(cgpu_detail::gradient<Grid_2d<T>, TVctr>, grid_2d, mx_i, dM_x, dM_y);
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				operator()(TVctr& mx_i, TVctr& dM_x, TVctr& dM_y)
				{
					auto d_grid_blk = grid_2d.d_grid_blk();
					gpu_detail::gradient<Grid_2d<T>, typename TVctr::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, mx_i, dM_x, dM_y);
				}
			#endif

			protected:
				Grid_2d<T> grid_2d;
				Stream<Dev> *stream;
		};

		/***************************** 2d affine transformation *********************************/
		template <class T, eDev Dev>
		class AT_2d
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				AT_2d(): stream(nullptr) {}

				AT_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					set_in_data(stream_i, grid_2d_i, bg_opt_i, bg_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					stream = stream_i;
					grid_2d = grid_2d_i;
					bg_opt = bg_opt_i;
					bg = bg_i;
				}

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, T>
				operator()(TVctr& mx_i, Mx_2x2<T> A, R_2d<T> txy, TVctr& mx_o)
				{
					if (fcn_is_zero<T>(txy) && is_I2x2<T>(A))
					{
						mx_o = mx_i;
					}

					// calculate background
					T bg = get_bg(mx_i);

					A = inv(A);
					txy = T(-1)*A*txy;

					stream->set_n_stream_act(grid_2d.nx);
					stream->set_grid(grid_2d.nx, grid_2d.ny);
					stream->exec_2d(cgpu_detail::at_2d<Grid_2d<T>, TVctr>, grid_2d, mx_i, A, txy, bg, mx_o);

					return bg;
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, T>
				operator()(TVctr& mx_i, Mx_2x2<T> A, R_2d<T> txy, TVctr& mx_o)
				{
					if (fcn_is_zero<T>(txy) && is_I2x2<T>(A))
					{
						mx_o = mx_i;
					}

					// calculate background
					T bg = get_bg(mx_i);

					A = inv(A);
					txy = T(-1)*A*txy;

					auto d_grid_blk = grid_2d.d_grid_blk();
					gpu_detail::at_2d<Grid_2d<T>, typename TVctr::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, mx_i, A, txy, bg, mx_o);

					return bg;
				}
			#endif

			protected:
				Grid_2d<T> grid_2d;
				Stream<Dev> *stream;
				eFil_Sel_Typ bg_opt;
				T bg;

				T get_bg(TVctr& M)
				{
					T bg_r = 0;
					switch (bg_opt) 
					{
						case efst_min:
							bg_r = fcn_min_element(M);
						break;
						case efst_max:
							bg_r = fcn_max_element(M);
						break;
						case efst_mean:
							bg_r = fcn_mean(M);
						break;
						case efst_min_mean:
							bg_r = 0.5*(fcn_mean(M) + fcn_min_element(M));
						break;
						case efst_max_mean:
							bg_r = 0.5*(fcn_mean(M) + fcn_max_element(M));
						break;
						case efst_user_def:
							bg_r = bg;
						break;
					}

					return bg_r;
				}
		};

 		/************************************ translation **************************************/
		template <class T, eDev Dev>
		class Tr_2d: public AT_2d<T, Dev>
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Tr_2d(): AT_2d<T, Dev>() {}

				Tr_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					this->set_in_data(stream_i, grid_2d_i, bg_opt_i, bg_i);
				}

				T operator()(TVctr& mx_i, R_2d<T> txy, TVctr& mx_o)
				{
					Mx_2x2<T> A;
					A(1, 1) = 1.0;
					A(2, 1) = 0;
					A(1, 2) = 0;
					A(2, 2) = 1.0;

					T bg = AT_2d<T, Dev>::operator()(mx_i, A, txy, mx_o);

					return bg;
				}

		};

		/***************************************************************************************/
		template <class T, eDev Dev>
		class Rot_Tr_2d: public AT_2d<T, Dev>
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Rot_Tr_2d(): AT_2d<T, Dev>() {}

				Rot_Tr_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					this->set_in_data(stream_i, grid_2d_i, bg_opt_i, bg_i);
				}

				T operator()(TVctr& mx_i, T theta, R_2d<T> pr, R_2d<T> txy, TVctr& mx_o)
				{
					auto A = fcn_rot_mx_2d(theta);
					txy = pr - A*(pr-txy);

					T bg = AT_2d<T, Dev>::operator()(mx_i, A, txy, mx_o);

					return bg;
				}
		};

		template <class T, eDev Dev>
		class Tr_Rot_2d: public AT_2d<T, Dev>
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Tr_Rot_2d(): AT_2d<T, Dev>() {}

				Tr_Rot_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					this->set_in_data(stream_i, grid_2d_i, bg_opt_i, bg_i);
				}

				T operator()(TVctr& mx_i, T theta, R_2d<T> pr, R_2d<T> txy, TVctr& mx_o)
				{
					auto A = fcn_rot_mx_2d(theta);
					txy = pr - A*pr + txy;

					T bg = AT_2d<T, Dev>::operator()(mx_i, A, txy, mx_o);

					return bg;
				}
		};

 		/************************************** scy_shx_tr *************************************/
 		template <class T, eDev Dev>
		class Scy_Shx_Tr_2d: public AT_2d<T, Dev>
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Scy_Shx_Tr_2d(): AT_2d<T, Dev>() {}

				Scy_Shx_Tr_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					this->set_in_data(stream_i, grid_2d_i, bg_opt_i, bg_i);
				}

				T operator()(TVctr& mx_i, T scy, T shx, R_2d<T> txy, TVctr& mx_o)
				{
					Mx_2x2<T> A;
					A(1, 1) = 1;
					A(2, 1) = 0;
					A(1, 2) = shx;
					A(2, 2) = scy;

					txy = A*txy;

					T bg = AT_2d<T, Dev>::operator()(mx_i, A, txy, mx_o);

					return bg;
				}
		};

 		template <class T, eDev Dev>
		class Tr_Scy_Shx_2d: public AT_2d<T, Dev>
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Tr_Scy_Shx_2d(): AT_2d<T, Dev>() {}

				Tr_Scy_Shx_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					this->set_in_data(stream_i, grid_2d_i, bg_opt_i, bg_i);
				}

				T operator()(TVctr& mx_i, T scy, T shx, R_2d<T> txy, TVctr& mx_o)
				{
					Mx_2x2<T> A;
					A(1, 1) = 1;
					A(2, 1) = 0;
					A(1, 2) = shx;
					A(2, 2) = scy;

					T bg = AT_2d<T, Dev>::operator()(mx_i, A, txy, mx_o);

					return bg;
				}
		};

 		/************************************ scy_scx_tr ***************************************/
 		template <class T, eDev Dev>
		class Scy_Scx_Tr_2d: public AT_2d<T, Dev>
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Scy_Scx_Tr_2d(): AT_2d<T, Dev>() {}

				Scy_Scx_Tr_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					this->set_in_data(stream_i, grid_2d_i, bg_opt_i, bg_i);
				}

				T operator()(TVctr& mx_i, T scy, T scx, R_2d<T> txy, TVctr& mx_o)
				{
					Mx_2x2<T> A;
					A(1, 1) = scx;
					A(2, 1) = 0;
					A(1, 2) = 0;
					A(2, 2) = scy;

					txy = A*txy;

					T bg = AT_2d<T, Dev>::operator()(mx_i, A, txy, mx_o);

					return bg;
				}
		};

 		template <class T, eDev Dev>
		class Tr_Scy_Scx_2d: public AT_2d<T, Dev>
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Tr_Scy_Scx_2d(): AT_2d<T, Dev>() {}

				Tr_Scy_Scx_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					this->set_in_data(stream_i, grid_2d_i, bg_opt_i, bg_i);
				}

				T operator()(TVctr& mx_i, T scy, T scx, R_2d<T> txy, TVctr& mx_o)
				{
					Mx_2x2<T> A;
					A(1, 1) = scx;
					A(2, 1) = 0;
					A(1, 2) = 0;
					A(2, 2) = scy;

					T bg = AT_2d<T, Dev>::operator()(mx_i, A, txy, mx_o);

					return bg;
				}
		};

		/************************************* shift *******************************************/
		template <class T, eDev Dev>
		class Sft_1d
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				static const eDev device = Dev;

				Sft_1d():fft_1d(nullptr) {}

				Sft_1d(FFT<T, Dev> *fft_1d_i, Grid_1d<T>& grid_1d_i)
				{
					set_in_data(fft_1d_i, grid_1d_i);
				}

				inline
				void set_in_data(FFT<T, Dev> *fft_1d_i, Grid_1d<T>& grid_1d_i)
				{
					fft_1d = fft_1d_i;
					grid_1d = grid_1d_i;
				}

				void set_fft_plan()
				{
					fft_1d->create_plan_1d(grid_1d.nx);
				}

				void operator()(T xs, TVctr_c& Im)
				{
					fcn_fftsft_1d(grid_1d, Im);
					fft_1d->forward(Im);
					fcn_exp_g_factor_1d(grid_1d, -c_2pi<T>*xs, Im, Im);
					fft_1d->inverse(Im);
					fcn_fftsft_1d(grid_1d, Im);
				}

				void operator()(T xs, TVctr_r &Im)
				{
					TVctr_c Im_c(Im.begin(), Im.end());
					T Im_min = fcn_min_element(Im);

					this->operator()(xs, Im_c);
					fcn_assign_real(Im_c, Im, Im_min);
				}

				void cleanup()
				{
					fft_1d->cleanup();
				}

			protected:
				FFT<T, Dev> *fft_1d;
				Grid_1d<T> grid_1d;
		};

		template <class T, eDev Dev>
		class Sft_2d_BC
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				static const eDev device = Dev;

				Sft_2d_BC():stream(nullptr), fft_2d(nullptr) {}

				Sft_2d_BC(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					set_in_data(stream_i, fft_2d_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					stream = stream_i;
					fft_2d = fft_2d_i;
					grid_2d = grid_2d_i;
				}

				void set_fft_plan()
				{
					fft_2d->create_plan_1d_batch(grid_2d.ny, grid_2d.nx, stream->size());
				}

				void operator()(T alpha, TVctr_r ys, TVctr_c& Im)
				{
					fcn_fftsft_bc_2d(*stream, grid_2d, Im);
					fft_2d->forward(Im);
					fcn_exp_g_factor_2d_bc(*stream, grid_2d, -c_2pi<T>*alpha, ys, Im, Im);
					fft_2d->inverse(Im);
					fcn_fftsft_bc_2d(*stream, grid_2d, Im);
				}

				void operator()(T alpha, TVctr_r ys, TVctr_r &Im)
				{
					TVctr_c Im_c(Im.begin(), Im.end());
					T Im_min = fcn_min_element(Im);

					this->operator()(alpha, ys, Im_c);
					fcn_assign_real(Im_c, Im, Im_min);
				}

				void cleanup()
				{
					fft_2d->cleanup();
				}

			protected:
				Stream<Dev> *stream;
				FFT<T, Dev> *fft_2d;
				Grid_2d<T> grid_2d;
		};

		template <class T, eDev Dev>
		class Sft_2d
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				static const eDev device = Dev;

				Sft_2d():stream(nullptr), fft_2d(nullptr) {}

				Sft_2d(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					set_in_data(stream_i, fft_2d_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					stream = stream_i;
					fft_2d = fft_2d_i;
					grid_2d = grid_2d_i;
				}

				void set_fft_plan()
				{
					fft_2d->create_plan_2d(grid_2d.ny, grid_2d.nx, stream->size());
				}

				void operator()(R_2d<T> p, TVctr_c& Im)
				{
					fcn_fftsft_2d(*stream, grid_2d, Im);
					fft_2d->forward(Im);
					fcn_exp_g_factor_2d(*stream, grid_2d, T(-c_2pi)*p, Im, Im);
					fft_2d->inverse(Im);
					fcn_fftsft_2d(*stream, grid_2d, Im);
				}

				void operator()(R_2d<T> p, TVctr_r &Im)
				{
					TVctr_c Im_c(Im.begin(), Im.end());
					T Im_min = fcn_min_element(Im);

					this->operator()(p, Im_c);
					fcn_assign_real(Im_c, Im, Im_min);
				}

				void cleanup()
				{
					fft_2d->cleanup();
				}
			protected:
				Stream<Dev> *stream;
				FFT<T, Dev> *fft_2d;
				Grid_2d<T> grid_2d;
		};

		/******************************** phase correlation ************************************/
		template <class T, eDev Dev>
		class Pcf_1d
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				static const eDev device = Dev;

				Pcf_1d(): fft_1d(nullptr) {}

				Pcf_1d(FFT<T, Dev> *fft_1d_i, Grid_1d<T>& grid_1d_i)
				{
					set_in_data(fft_1d_i, grid_1d_i);
				}

				inline
				void set_in_data(FFT<T, Dev> *fft_1d_i, Grid_1d<T>& grid_1d_i)
				{
					fft_1d = fft_1d_i;
					grid_1d = grid_1d_i;

					PN_Fact pn;

					dt_int32 nx = pn(2*grid_1d.nx-1, edst_greater_than);

					grid_1d_e.set_in_data(nx, grid_1d.drx*nx);

					M_r_c.resize(grid_1d_e.nx);
					M_s_c.resize(grid_1d_e.nx);
				}

				void set_fft_plan()
				{
					fft_1d->create_plan_1d(grid_1d.nx, 1);
					fft_1d_e.create_plan_1d(grid_1d_e.nx, 1);
				}

				void operator()(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, 
				Region _Rect_1d<T> bd, dt_bool b_pv, TVctr_r& mx_o)
				{
					thrust::fill(M_r_c.begin(), M_r_c.end(), T(0));
					thrust::fill(M_s_c.begin(), M_s_c.end(), T(0));

					preprocessing(M_r, p, bd, M_r_c);
					preprocessing(M_s, p, bd, M_s_c);

					TVctr_c& pcf = M_s_c;

					// shift matrix
					fcn_fftsft_1d(grid_1d_e, M_r_c);
					fcn_fftsft_1d(grid_1d_e, M_s_c);

					// fft_1d_e
					fft_1d_e.forward(M_r_c);
					fft_1d_e.forward(M_s_c);

					pcf_g(M_r_c, M_s_c, sigma_g, pcf);

					fft_1d_e.inverse(pcf);

					// shift pcf
					fcn_fftsft_1d(grid_1d_e, pcf);

					dt_int32 ix_s = (grid_1d_e.nx-grid_1d.nx)/2;
					this->fcn_assign_real(pcf, ix_s, mx_o, b_pv);
				}

				void cleanup()
				{
					fft_1d->destroy_plan();
					fft_1d_e.cleanup();
				}

			protected:

				void fcn_assign_real(TVctr_c& mx_i, dt_int32 ix_s, TVctr_r& mx_o, dt_bool b_pos = false)
				{
					auto first = mx_i.begin() + ix_s;
					auto last = first + grid_1d.nx;

					if (b_pos)
					{
						thrust::transform(first, last, mx_o.begin(), cgpu_fctr::assign_max_real<T>(0));
					}
					else
					{
						thrust::transform(first, last, mx_o.begin(), cgpu_fctr::assign_real<T>());
					}
				}

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				preprocessing(TVctr_r& mx_i, T p, Region _Rect_1d<T>& bd, TVctr_c& mx_o)
				{
					dt_int32 ix_s = (grid_1d_e.nx-grid_1d.nx)/2;

					Wd_Butwth_1d<T> bw_1d(bd, bd.radius_p(p), 32);

					for(auto ix = 0; ix < grid_1d.nx; ix++)
					{
						cgpu_detail::fcn_rs_pcf_1d_dp(ix, ix_s, grid_1d, bw_1d, mx_i, mx_o);
					}
				}

				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				pcf_g(TVctr_c& M_r_c, TVctr_c& M_s_c, T sigma_g, TVctr_c& pcf)
				{
					Wd_Gauss_1d<T> gs_1d(0, sigma_g);

					for(auto ix = 0; ix < grid_1d_e.nx; ix++)
					{
						cgpu_detail::fcn_fs_pcf_1d_dp(ix, grid_1d_e, gs_1d, M_r_c, M_s_c, pcf);
					}
				}

			/*************************************** device ****************************************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				preprocessing(TVctr_r& mx_i, T p, Region _Rect_1d<T>& bd, TVctr_c& mx_o)
				{
					dt_int32 ix_s = (grid_1d_e.nx-grid_1d.nx)/2;

					Wd_Butwth_1d<T> bw_1d(bd, bd.radius_p(p), 32);

					auto d_grid_blk = grid_1d.d_grid_blk();
					gpu_detail::fcn_rs_pcf_1d_dp<Grid_1d<T>, typename TVctr_c::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(ix_s, grid_1d, bw_1d, mx_i, mx_o);
				}

				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				pcf_g(TVctr_c& M_r_c, TVctr_c& M_s_c, T sigma_g, TVctr_c& pcf)
				{
					Wd_Gauss_1d<T> gs_1d(0, sigma_g);

					auto d_grid_blk = grid_1d_e.d_grid_blk();
					gpu_detail::fcn_fs_pcf_1d_dp<Grid_1d<T>, typename TVctr_c::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_1d_e, gs_1d, M_r_c, M_s_c, pcf);
				}
			#endif

				TVctr_c M_r_c;
				TVctr_c M_s_c;

				FFT<T, Dev> *fft_1d;
				Grid_1d<T> grid_1d;

			private:
				FFT<T, Dev> fft_1d_e;
				Grid_1d<T> grid_1d_e;
		};

		template <class T, eDev Dev>
		class Pcf_2d_BC
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				static const eDev device = Dev;

				Pcf_2d_BC():stream(nullptr), fft_2d(nullptr) {}

				Pcf_2d_BC(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					set_in_data(stream_i, fft_2d_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					stream = stream_i;
					fft_2d = fft_2d_i;
					grid_2d = grid_2d_i;
					grid_1d.set_in_data(grid_2d.ny, grid_2d.bs_y);

					PN_Fact pn;

					dt_int32 ny = pn(2*grid_2d.ny-1, edst_greater_than);

					grid_2d_e.set_in_data(grid_2d.nx, ny, grid_2d.bs_x, grid_2d.dry*ny);
					grid_1d_e.set_in_data(grid_2d_e.ny, grid_2d_e.bs_y);

					M_r_c.resize(grid_2d_e.size());
					M_s_c.resize(grid_2d_e.size());
				}

				void set_fft_plan()
				{
					fft_2d->create_plan_1d_batch(grid_2d.ny, grid_2d.nx, stream->size());
					fft_2d_e.create_plan_1d_batch(grid_2d_e.ny, grid_2d_e.nx, stream->size());
				}

				void operator()(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, 
				Region _Rect_1d<T> bd, dt_bool b_pv, TVctr_r& mx_o)
				{
					thrust::fill(M_r_c.begin(), M_r_c.end(), T(0));
					thrust::fill(M_s_c.begin(), M_s_c.end(), T(0));

					preprocessing(M_r, p, bd, M_r_c);
					preprocessing(M_s, p, bd, M_s_c);

					TVctr_c& pcf = M_s_c;

					// shift matrix
					fcn_fftsft_bc_2d(*stream, grid_2d_e, M_r_c);
					fcn_fftsft_bc_2d(*stream, grid_2d_e, M_s_c);

					// fft_2d
					fft_2d_e.forward(M_r_c);
					fft_2d_e.forward(M_s_c);

					pcf_g(M_r_c, M_s_c, sigma_g, pcf);
					fft_2d_e.inverse(pcf);

					// shift pcf
					fcn_fftsft_bc_2d(*stream, grid_2d_e, pcf);

					dt_int32 iy_s = (grid_2d_e.ny-grid_2d.ny)/2;
					this->fcn_assign_real(pcf, iy_s, mx_o, b_pv);
				}

				void cleanup()
				{
					fft_2d->destroy_plan();
					fft_2d_e.cleanup();
				}

			protected:

				TVctr_r gauss_vector_1d(Grid_1d<T>& grid_1d, T sigma_g, dt_bool b_norm = false)
				{
					Wd_Gauss_1d<T> gs_1d(0, sigma_g);

					TVctr_r fg;
					fg.reserve(grid_1d.nx);

					for(auto ix = 0; ix < grid_1d.nx; ix++)
					{
						T g2 = grid_1d.g2_sft(ix);
						T v = gs_1d(g2);
						fg.push_back(v);
					}
					return fg;
				}

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				preprocessing(TVctr_r& mx_i, T p, Region _Rect_1d<T>& bd, TVctr_c& mx_o)
				{
					dt_int32 iy_s = (grid_2d_e.ny-grid_2d.ny)/2;

					auto fh = fcn_wd_butwth_1d<TVctr_r>(grid_1d, bd.radius_p(p), 32, false, bd);

					stream->set_n_stream_act(grid_2d.nx);
					stream->set_grid(grid_2d.nx, grid_2d.ny);
					stream->exec_2d(cgpu_detail::pcf_2d_bc_pp<Grid_2d<T>, TVctr_r, TVctr_c>, iy_s, grid_2d, grid_2d_e, mx_i, fh, mx_o);
				}

				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				pcf_g(TVctr_c& M_r_c, TVctr_c& M_s_c, T sigma_g, TVctr_c& pcf)
				{
					auto fg = gauss_vector_1d(grid_1d_e, sigma_g);

					stream->set_n_stream_act(grid_2d_e.nx);
					stream->set_grid(grid_2d_e.nx, grid_2d_e.ny);
					stream->exec_2d(cgpu_detail::pcf_2d_bc_gaussian<Grid_2d<T>, TVctr_r, TVctr_c>, grid_2d_e, M_r_c, M_s_c, fg, pcf);
				}

				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				fcn_assign_real(TVctr_c& mx_i, dt_int32 iy_s, TVctr_r& mx_o, dt_bool b_pos = false)
				{
					for(auto ix = 0; ix < grid_2d.nx; ix++)
					{
						for(auto iy = 0; iy < grid_2d.ny; iy++)
						{
							dt_int32 ixy_i = grid_2d_e.sub_2_ind(ix, iy+iy_s);
							dt_int32 ixy_o = grid_2d.sub_2_ind(ix, iy);
							auto v = mx_i[ixy_i].real();
							mx_o[ixy_o] = (!b_pos || (v>0))?v:0;
						}
					}
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				preprocessing(TVctr_r& mx_i, T p, Region _Rect_1d<T>& bd, TVctr_c& mx_o)
				{
					dt_int32 iy_s = (grid_2d_e.ny-grid_2d.ny)/2;

					auto fh_h = fcn_wd_butwth_1d<Vctr<T, edev_cpu>>(grid_1d, bd.radius_p(p), 32, false, bd);
					TVctr_r fh = fh_h;

					auto d_grid_blk = grid_2d.d_grid_blk();
					gpu_detail::pcf_2d_bc_pp<Grid_2d<T>, typename TVctr_c::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(iy_s, grid_2d, grid_2d_e, mx_i, fh, mx_o);
				}

				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				pcf_g(TVctr_c& M_r_c, TVctr_c& M_s_c, T sigma_g, TVctr_c& pcf)
				{
					auto fg_h = gauss_vector_1d(grid_1d_e, sigma_g);
					TVctr_r fg = fg_h;

					auto d_grid_blk = grid_2d_e.d_grid_blk();
					gpu_detail::pcf_2d_bc_gaussian<Grid_2d<T>, typename TVctr_c::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d_e, M_r_c, M_s_c, fg, pcf);
				}

				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				fcn_assign_real(TVctr_c& mx_i, dt_int32 iy_s, TVctr_r& mx_o, dt_bool b_pos = false)
				{
					auto d_grid_blk = grid_2d.d_grid_blk();
					gpu_detail::pcf_2d_bc_assign_real<Grid_2d<T>, typename TVctr_c::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(iy_s, grid_2d, grid_2d_e, mx_i, mx_o, b_pos);
				}
			#endif

				Vctr<T_c, Dev> M_r_c;
				Vctr<T_c, Dev> M_s_c;

				Stream<Dev> *stream;
				FFT<T, Dev> *fft_2d;
				Grid_2d<T> grid_2d;
				Grid_1d<T> grid_1d;

			private:
				FFT<T, Dev> fft_2d_e;
				Grid_2d<T> grid_2d_e;
				Grid_1d<T> grid_1d_e;
		};

		template <class T, eDev Dev>
		class Pcf_2d
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				static const eDev device = Dev;

				Pcf_2d():stream(nullptr), fft_2d(nullptr) {}

				Pcf_2d(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					set_in_data(stream_i, fft_2d_i, grid_2d_i, bg_opt_i, bg_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0)
				{
					stream = stream_i;
					fft_2d = fft_2d_i;
					grid_2d = grid_2d_i;

					af_tr.set_in_data(stream_i, grid_2d_i, bg_opt_i, bg_i);

					M_r_c.resize(grid_2d.size());
					M_s_c.resize(grid_2d.size());
					pcf_r.resize(grid_2d.size());

					Fit_Ellipt_Gauss_2d.init_variables(grid_2d, 1.0);
				}

				void set_fft_plan()
				{
					fft_2d->create_plan_2d(grid_2d.ny, grid_2d.nx, stream->size());
				}

				void operator()(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, 
				Region _Rect_2d<T> bd, TVctr_r &pcf_o)
				{
					preprocessing(M_r, p, bd, M_r_c);
					preprocessing(M_s, p, bd, M_s_c);

					// shift matrix
					fcn_fftsft_2d(*stream, grid_2d, M_r_c);
					fcn_fftsft_2d(*stream, grid_2d, M_s_c);

					// fft_2d
					fft_2d->forward(M_r_c);
					fft_2d->forward(M_s_c);

					TVctr_c& pcf_c = M_s_c;
					pcf_g(M_r_c, M_s_c, sigma_g, pcf_c);
					fft_2d->inverse(pcf_c);

					// shift pcf
					fcn_fftsft_2d(*stream, grid_2d, pcf_c);

					fcn_assign_real(pcf_c, pcf_o);
				}

				void operator()(TVctr_r &M_r, TVctr_r &M_s, Mx_2x2<T> A, R_2d<T> txy, 
				T p, T sigma_g, Region _Rect_2d<T> bd, TVctr_r &pcf_o)
				{
					auto &M_s_t = pcf_r;
					af_tr(M_s, A, txy, M_s_t);
					this->operator()(M_r, M_s_t, p, sigma_g, bd, pcf_o);
				}

				T chi_2(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, 
				Region _Rect_2d<T> bd)
				{
					preprocessing(M_r, p, bd, M_r_c);
					preprocessing(M_s, p, bd, M_s_c);

					// shift matrix
					fcn_fftsft_2d(*stream, grid_2d, M_r_c);
					fcn_fftsft_2d(*stream, grid_2d, M_s_c);

					// fft_2d
					fft_2d->forward(M_r_c);
					fft_2d->forward(M_s_c);

					TVctr_c& pcf_c = M_s_c;
					pcf_g(M_r_c, M_s_c, sigma_g, pcf_c);
					fft_2d->inverse(pcf_c);

					return mt::fcn_mean_abs_real(pcf_c);
				}

				T chi_2(TVctr_r &M_r, TVctr_r &M_s, Mx_2x2<T> A, R_2d<T> txy, 
				T p, T sigma_g, Region _Rect_2d<T> bd)
				{
					auto &M_s_t = pcf_r;
					af_tr(M_s, A, txy, M_s_t);
					return this->chi_2(M_r, M_s_t, p, sigma_g, bd);
				}

				T chi_2(TVctr_r &M_r, TVctr_r &M_s, Mx_2x2<T> A, R_2d<T> txy, 
				T p, T sigma_g, T radius, Region _Rect_2d<T> bd, Vctr<T, edev_cpu>& coef)
				{
					// pcf
					this->operator()(M_r, M_s, A, txy, p, sigma_g, bd, this->pcf_r);

					// get coefficients
					coef = fit_coef(this->pcf_r, txy, sigma_g, radius);

					return mt::fcn_mean_abs(this->pcf_r);
				}

				void cleanup()
				{
					fft_2d->cleanup();
				}

			protected:

				Vctr<T, edev_cpu> fit_coef(TVctr_r &pcf_r, R_2d<T> txy, T sigma_g, T radius)
				{
					// get maximum position
					R_2d<T> r_c = Fit_Ellipt_Gauss_2d.fd_max_peak_pos(pcf_r);

					// fitting
					T sigma_r = 1.0/(c_2pi<T>*sigma_g);
					auto coef = Fit_Ellipt_Gauss_2d.fit(pcf_r, r_c, sigma_r, radius);
					coef[0] = txy.x + this->grid_2d.bs_x_h();
					coef[1] = txy.y + this->grid_2d.bs_y_h();

					return coef;
				}

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				preprocessing(TVctr_r& mx_i, T p, Region _Rect_2d<T>& bd, TVctr_c& mx_o)
				{
					Wd_Butwth_2d<T> bw_2d(bd, bd.radius_p(p), 16);

					stream->set_n_stream_act(grid_2d.nx);
					stream->set_grid(grid_2d.nx, grid_2d.ny);
					stream->exec_2d(cgpu_detail::fcn_rs_pcf_2d_dp<Grid_2d<T>, TVctr_r, TVctr_c>, grid_2d, bw_2d, mx_i, mx_o);
				}

				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				pcf_g(TVctr_c& M_r_c, TVctr_c& M_s_c, T sigma_g, TVctr_c& pcf)
				{
					Wd_Gauss_2d<T> gs_2d(R_2d<T>(), sigma_g);

					stream->set_n_stream_act(grid_2d.nx);
					stream->set_grid(grid_2d.nx, grid_2d.ny);
					stream->exec_2d(cgpu_detail::fcn_fs_pcf_2d_dp<Grid_2d<T>, TVctr_c>, grid_2d, gs_2d, M_r_c, M_s_c, pcf);
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				preprocessing(TVctr_r& mx_i, T p, Region _Rect_2d<T>& bd, TVctr_c& mx_o)
				{
					Wd_Butwth_2d<T> bw_2d(bd, bd.radius_p(p), 16);

					auto d_grid_blk = grid_2d.d_grid_blk();
					gpu_detail::fcn_rs_pcf_2d_dp<Grid_2d<T>, typename TVctr_c::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, bw_2d, mx_i, mx_o);
				}

				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				pcf_g(TVctr_c& M_r_c, TVctr_c& M_s_c, T sigma_g, TVctr_c& pcf)
				{
					Wd_Gauss_2d<T> gs_2d(R_2d<T>(), sigma_g);

					auto d_grid_blk = grid_2d.d_grid_blk();
					gpu_detail::fcn_fs_pcf_2d_dp<Grid_2d<T>, typename TVctr_c::value_type><<< d_grid_blk.grid, d_grid_blk.blk >>>(grid_2d, gs_2d, M_r_c, M_s_c, pcf);
				}
			#endif

				Fit_Ellipt_Gauss_2d<T, Dev> Fit_Ellipt_Gauss_2d;

				TVctr_c M_r_c;
				TVctr_c M_s_c;
				TVctr_r pcf_r;

				Stream<Dev> *stream;
				FFT<T, Dev> *fft_2d;
				Grid_2d<T> grid_2d;

				AT_2d<T, Dev> af_tr;
		};

		/*********************************** find shift ****************************************/
		template <class T, eDev Dev>
		class Fd_Sft_1d: public Pcf_1d<T, Dev>
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				Fd_Sft_1d(): Pcf_1d<T, Dev>() {}

				Fd_Sft_1d(FFT<T, Dev> *fft_1d_i, Grid_1d<T>& grid_1d_i): Pcf_1d<T, Dev>()
				{
					set_in_data(fft_1d_i, grid_1d_i);
				}

				inline
				void set_in_data(FFT<T, Dev> *fft_1d_i, Grid_1d<T>& grid_1d_i)
				{
					Pcf_1d<T, Dev>::set_in_data(fft_1d_i, grid_1d_i);
					sft_1d.set_in_data(fft_1d_i, grid_1d_i);
				}

				T operator()(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, T dx, 
				Region _Rect_1d<T> bd, dt_int32 nit_pcf)
				{
					T sigma_r = 1.0/(c_2pi<T>*sigma_g);
					T radius = ::fmax(3*this->grid_1d.drx, 0.9*sigma_r);

					TVctr_r M = M_s;
					TVctr_r pcf(M.size());

					if (fcn_is_nzero(dx))
					{
						sft_1d(dx, M);
					}

					for(auto it=0; it<nit_pcf; it++)
					{
						Pcf_1d<T, Dev>::operator()(M_r, M, p, sigma_g, bd, true, pcf);

						// get maximum position
						T x_c = cgpu_detail::fd_max_peak_pos(this->grid_1d, pcf);

						if (it==nit_pcf-1)
						{
							x_c = fit_max_pos_1d(this->grid_1d, pcf, x_c, sigma_r, radius);
						}

						T dx_t = -(x_c - this->grid_1d.bs_x_h());

						if (it<nit_pcf-1)
						{
							sft_1d(dx_t, M);
						}
						dx += dx_t;

						// calculated new borders
						bd.fcn_repl_bdr(dx);
					}

					return dx;
				}

			protected:
				Sft_1d<T, Dev> sft_1d;
		};

		template <class T, eDev Dev>
		class Fd_Sft_2d_BC: public Pcf_2d_BC<T, Dev>
		{
			public:
				using TB = std::pair<T, dt_bool>;
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;
				using TVctr_rh = Vctr<T, edev_cpu>;
				using TVctr_ch = Vctr<T_c, edev_cpu>;

				Fd_Sft_2d_BC(): Pcf_2d_BC<T, Dev>() {}

				Fd_Sft_2d_BC(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i): Pcf_2d_BC<T, Dev>()
				{
					set_in_data(stream_i, fft_2d_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					Pcf_2d_BC<T, Dev>::set_in_data(stream_i, fft_2d_i, grid_2d_i);
					gauss_cv_2d_bc.set_in_data(stream_i, fft_2d_i, grid_2d_i);
					sft_2d_bc.set_in_data(stream_i, fft_2d_i, grid_2d_i);
				}

				TVctr_rh operator()(TVctr_r &M_r_i, TVctr_r &M_s_i, T p, T sigma_g, TVctr_rh dx, 
				Region _Rect_2d<T> bd_2d, dt_int32 nit_pcf)
				{
					T sigma_r = 1.0/(c_2pi<T>*sigma_g);
					T radius = ::fmax(3*this->grid_1d.drx, 0.9*sigma_r);
					Region _Rect_1d<T> bd_1d(bd_2d.bs_y, bd_2d.yb_0, bd_2d.yb_e);
					Peaks peaks(this->grid_1d, &bd_1d);

					TVctr_r M_r = M_r_i;
					TVctr_r M_s = M_s_i;
					TVctr_r pcf(M_r.size());

					// Gaussian filter
					gauss_cv_2d_bc(sigma_r, M_r);
					gauss_cv_2d_bc(sigma_r, M_s);

					const dt_int32 ix_0 = this->grid_2d.rx_2_irx_cds(bd_2d.x_0());
					const dt_int32 ix_e = this->grid_2d.rx_2_irx_fds(bd_2d.x_e());

					dt_bool dx_s = false;
					for(auto ix = 0; ix<dx.size(); ix++)
					{
						dx_s = dx_s || fcn_is_nzero(dx[ix]);
					}
					if (dx_s)
					{
						sft_2d_bc(1, dx, M_s);
					}

					TVctr_rh dx_t(this->grid_1d.nx, T(0));

					for(auto it=0; it<nit_pcf; it++)
					{
						Pcf_2d_BC<T, Dev>::operator()(M_r, M_s, p, sigma_g, bd_1d, false, pcf);

						// get maximum distance
						T d_max = (it==0)?peaks.get_d_max_0(ix_0, ix_e, pcf):peaks.get_d_max(ix_0, ix_e, dx_t);
						d_max = ::fmax(sigma_r, d_max);

						for(auto ix=ix_0; ix<ix_e; ix++)
						{
							TVctr_rh y(pcf.begin()+ix*this->grid_1d.nx, pcf.begin()+(ix+1)*this->grid_1d.nx);

							auto py = peaks.find(ix, M_r, M_s, y, d_max);
							if (it==nit_pcf-1)
							{
								py = fit_max_pos_1d(this->grid_1d, y, py, sigma_r, radius);
							}
							dx_t[ix] = -(py - this->grid_1d.bs_x_h());
							dx[ix] += dx_t[ix];
						}

						// shift by column
						if (it<nit_pcf-1)
						{
							sft_2d_bc(1, dx_t, M_s);
						}

						// get new border
						T d_b = fabs(*thrust::max_element(dx.begin(), dx.end(), [](T a, T b){ return fabs(a)<fabs(b); }));

						// calculated new borders
						bd_1d.xb_0 = max(bd_1d.xb_0, d_b);
						bd_1d.xb_e = max(bd_1d.xb_e, d_b);
					}

					peaks.destroy_plan();

					return dx;
				}

			protected:
				Gauss_Cv_2d_BC<T, Dev> gauss_cv_2d_bc;
				Sft_2d_BC<T, Dev> sft_2d_bc;

			private:
				struct Peaks
				{
					public:
						using TVctr_ih = Vctr<dt_int32, edev_cpu>;

						Peaks(): y_c(0), bd_1d(nullptr) {}

						Peaks(Grid_1d<T>& grid_1d_i, Region _Rect_1d<T> *bd_1d_i)
						{
							set_in_data(grid_1d_i, bd_1d_i);
						}

						inline
						void set_in_data(Grid_1d<T>& grid_1d_i, Region _Rect_1d<T> *bd_1d_i)
						{
							grid_1d = grid_1d_i;
							bd_1d = bd_1d_i;
							fft_1d.create_plan_1d(grid_1d.nx, 1);
							sft_1d.set_in_data(&fft_1d, grid_1d);
						
							ix_pk.resize(grid_1d.nx);
							x_pk.resize(grid_1d.nx);
							y_pk.resize(grid_1d.nx);
						}

						T get_d_max_0(dt_int32 ix_0, dt_int32 ix_e, TVctr_r &pcf)
						{
							T x_c = grid_1d.bs_x_h();
							TVctr_rh dx;
							dx.reserve(grid_1d.nx);

							for(auto ix=ix_0; ix<ix_e; ix++)
							{
								TVctr_rh y(pcf.begin()+ix*grid_1d.nx, pcf.begin()+(ix+1)*grid_1d.nx);

								// find peaks
								fd_peaks(y);

								// add point
								dx.push_back(fabs(x_pk[idx_x_clt()]-x_c));
							}
							T dx_min, dx_max;
							get_limits(dx, dx_min, dx_max);

							// remove far away points
							auto it = std::remove_if (dx.begin(), dx.end(), [dx_min, dx_max](T a){return (a<dx_min)||(a>dx_max); });
							dx.resize(std::distance(dx.begin(), it));

							T dx_std = ::sqrt(fcn_variance(dx));
							T d_max = fcn_max_element(dx) + 2*dx_std;

							return ::fmax(10*grid_1d.drx, d_max);
						}

						T get_d_max(dt_int32 ix_0, dt_int32 ix_e, TVctr_rh &dx)
						{
							TVctr_rh dx_h(dx.begin()+ix_0, dx.begin()+ix_e);
							return ::fmax(10*grid_1d.drx, 3*::sqrt(fcn_variance(dx_h)));
						}

						T find(dt_int32 icol, TVctr_r &M_r, TVctr_r &M_s, TVctr_rh &y, T d_max)
						{
							const T x_c = grid_1d.bs_x_h();
							y_c = y[grid_1d.nx_h];

							// find peaks
							fd_peaks(y);

							// remove peaks
							remove_peaks(x_c-d_max, x_c+d_max);

							// return if x_pk is empty
							if (empty())
							{
								return x_c;
							}

							// return if there is only one element
							if (size() == 1)
							{
								return x_pk.front();
							}

							const dt_int32 ix_max = idx_y_max();
							const dt_int32 ix_clt = idx_x_clt();

							// return if closest position is equal to maximum intensity position
							// and the maximum intensity is greater than the second maximum intensity
							if (fabs(x_pk[ix_max]-x_pk[ix_clt])<1)
							{
								T y_thr = 0.61*y_pk[ix_max];
								for(auto ix = 0; ix<x_pk.size(); ix++)
								{
									if ((ix!=ix_max) && (y_thr<y_pk[ix]))
									{
										return fd_sft(icol, M_r, M_s);
									}
								}
								return x_pk[ix_max];
							}

							return fd_sft(icol, M_r, M_s);
						}

						void destroy_plan()
						{
							fft_1d.destroy_plan();
						}

					private:
						TVctr_ih ix_pk;
						TVctr_rh x_pk;
						TVctr_rh y_pk;

						Grid_1d<T> grid_1d;
						Region _Rect_1d<T> *bd_1d;

						FFT<T, edev_cpu> fft_1d;
						Sft_1d<T, edev_cpu> sft_1d;

						T y_c;

						dt_bool empty() const
						{
							return ix_pk.empty();
						}

						dt_int32 size() const
						{
							return ix_pk.size();
						}

						dt_int32 idx_y_max() const
						{
							dt_int32 idx_max = (std::max_element(y_pk.begin(), y_pk.end())-y_pk.begin());

							return idx_max;
						}

						dt_int32 idx_x_clt() const
						{
							T x_c = grid_1d.bs_x_h();

							dt_int32 idx_clt = 0;
							if (x_pk.size()>1)
							{
								idx_clt = std::min_element(x_pk.begin(), x_pk.end(), [x_c](T a, T b){return fabs(a-x_c)<fabs(b-x_c); })-x_pk.begin();
							}

							return idx_clt;
						}

						dt_int32 idx_x_fht() const
						{
							T x_c = grid_1d.bs_x_h();

							dt_int32 idx_clt = 0;
							if (x_pk.size()>1)
							{
								idx_clt = std::max_element(x_pk.begin(), x_pk.end(), [x_c](T a, T b){return fabs(a-x_c)<fabs(b-x_c); })-x_pk.begin();
							}

							return idx_clt;
						}

						void fd_peaks(TVctr_rh &y)
						{
							ix_pk.resize(y.size());
							x_pk.resize(y.size());
							y_pk.resize(y.size());

							const T y_thr = 0;
							dt_int32 ic = 0;
							for(auto ix=1; ix<y.size()-1; ix++)
							{
								if (y[ix]>y_thr)
								{
									if ((y[ix-1]<y[ix]) && (y[ix+1]<y[ix]))
									{
										ix_pk[ic] = ix;
										x_pk[ic] = grid_1d.rx(ix);
										y_pk[ic] = y[ix];
										ic++;
									}
								}
							}

							ix_pk.resize(ic);
							x_pk.resize(ic);
							y_pk.resize(ic);
						}

						void remove_peaks(T x_min, T x_max)
						{
							dt_int32 ic = 0;
							for(auto ix = 0; ix<size(); ix++)
							{
								if ((x_min<=x_pk[ix]) && (x_pk[ix]<=x_max))
								{
									ix_pk[ic] = ix_pk[ix];
									x_pk[ic] = x_pk[ix];
									y_pk[ic] = y_pk[ix];
									ic++;
								}
							}

							ix_pk.resize(ic);
							x_pk.resize(ic);
							y_pk.resize(ic);
						}

						void add_central_point()
						{
							const dt_int32 ix_c = grid_1d.nx_h;
							for(auto ix = 0; ix<ix_pk.size(); ix++)
							{
								if (ix_pk[ix]==ix_c)
								{
									return;
								}
							}

							ix_pk.push_back(ix_c);
							x_pk.push_back(grid_1d.rx(ix_c));
							y_pk.push_back(y_c);
						}

						void get_limits(TVctr_rh &y, T& y_min, T& y_max, T f = 4.0)
						{				
							// get fcn_mean and standard deviation
							T y_mean, y_std;
							fcn_mean_var(y, y_mean, y_std);
							y_std = ::sqrt(y_std);

							// set limits
							y_min = y_mean - f*y_std;
							y_max = y_mean + f*y_std;
						};

						T fd_sft(dt_int32 icol, TVctr_r &M_r, TVctr_r &M_s)
						{
							const T x_c = grid_1d.bs_x_h();

							if ((icol==0)||(icol==grid_1d.nx-1))
							{
								return x_c;
							}

							const dt_int32 ix_max = idx_y_max();
							const dt_int32 ix_clt = idx_x_clt();

							// select range
							const T x_min = ::fmin(x_pk[ix_clt], x_pk[ix_max])- grid_1d.drx;
							const T x_max = ::fmax(x_pk[ix_clt], x_pk[ix_max]) + grid_1d.drx;

							// remove far away peaks
							remove_peaks(x_min, x_max);

							// add central point
							add_central_point();

							// find maximum diplacement
							const T d_fht = fabs(x_pk[idx_x_fht()]-x_c);

							// set borders
							Region _Rect_1d<T> bd_l = *bd_1d;
							bd_l.xb_0 = max(bd_l.xb_0, d_fht);
							bd_l.xb_e = max(bd_l.xb_e, d_fht);

							// get indexes
							const dt_int32 ix_0 = grid_1d.rx_2_irx_cds(bd_l.x_0());
							const dt_int32 ix_e = grid_1d.rx_2_irx_fds(bd_l.x_e());

							// get reference data
							TVctr_rh yr_b(M_r.begin()+(icol-1)*grid_1d.nx, M_r.begin()+icol*grid_1d.nx);
							TVctr_rh yr(M_r.begin()+icol*grid_1d.nx, M_r.begin()+(icol+1)*grid_1d.nx);
							TVctr_rh yr_n(M_r.begin()+(icol+1)*grid_1d.nx, M_r.begin()+(icol+2)*grid_1d.nx);

							// get shifted data
							TVctr_rh ys_0(M_s.begin()+icol*grid_1d.nx, M_s.begin()+(icol+1)*grid_1d.nx);

							TVctr_rh chi2_pk(size());

							for(auto ix_pk=0; ix_pk<size(); ix_pk++)
							{
								TVctr_rh ys = ys_0;

								// shift 1d
								T dx = x_pk[ix_pk] - x_c;
								sft_1d(-dx, ys);

								// cost function
								const T alpha = 0.5;
								const T beta = 1-alpha;
								T chi2 = 0;
								T chi2_ee = 0;
								for(auto ix=ix_0; ix<ix_e; ix++)
								{
									T ee_f = pow(yr[ix]-ys[ix], 2);
									T ee_df = pow(yr_n[ix]+yr_b[ix]-2*ys[ix], 2);
									T chi2_p = alpha*ee_f + beta*ee_df;
									fcn_kh_sum(chi2, chi2_p, chi2_ee);
								}
								chi2_pk[ix_pk] = chi2;
							}
							dt_int32 ix_min = std::min_element(chi2_pk.begin(), chi2_pk.end())-chi2_pk.begin();

							return x_pk[ix_min];
						}
				};
		};

		template <class T, eDev Dev>
		class Fd_Sft_2d: public Pcf_2d<T, Dev>
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				Fd_Sft_2d(): Pcf_2d<T, Dev>() {}

				Fd_Sft_2d(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i): Pcf_2d<T, Dev>()
				{
					set_in_data(stream_i, fft_2d_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					Pcf_2d<T, Dev>::set_in_data(stream_i, fft_2d_i, grid_2d_i);
					sft_2d.set_in_data(stream_i, fft_2d_i, grid_2d_i);
				}

				R_2d<T> operator()(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, R_2d<T> dr, 
				Region _Rect_2d<T> bd, dt_int32 nit_pcf)
				{
					T sigma_r = 1.0/(c_2pi<T>*sigma_g);
					T radius = ::fmax(3*this->grid_2d.dR_min(), 0.9*sigma_r);

					TVctr_r M = M_s;
					TVctr_r pcf(M.size());

					if (fcn_is_nzero(dr))
					{
						sft_2d(dr, M);
					}

					for(auto it=0; it<nit_pcf; it++)
					{
						Pcf_2d<T, Dev>::operator()(M_r, M, p, sigma_g, bd, pcf);

						// get maximum position
						R_2d<T> r_c = this->Fit_Ellipt_Gauss_2d.fd_max_peak_pos(pcf);

						if (it==nit_pcf-1)
						{
							r_c = this->Fit_Ellipt_Gauss_2d.fit_peak_pos(pcf, r_c, sigma_r, radius);
						}

						R_2d<T> dr_t = -(r_c - this->grid_2d.bs_h());

						if (it<nit_pcf-1)
						{
							sft_2d(dr_t, M);
						}
						dr += dr_t;

						// calculated new borders
						// bd.fcn_repl_bdr(dr);
					}

					return dr;
				}

			protected:
				Sft_2d<T, Dev> sft_2d;
		};

		template <class T, eDev Dev>
		class Fd_Tr_2d: public Pcf_2d<T, Dev>
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				Fd_Tr_2d(): Pcf_2d<T, Dev>(), b_fit(true) {}

				Fd_Tr_2d(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i, 
				eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0, dt_bool b_fit_i = true): 
				Pcf_2d<T, Dev>(stream_i, fft_2d_i, grid_2d_i, bg_opt_i, bg_i), b_fit(b_fit_i) {}

				R_2d<T> operator()(TVctr_r &M_r, TVctr_r &M_s, Mx_2x2<T> A, R_2d<T> txy, 
				T p, T sigma_g, T radius_f, Region _Rect_2d<T> bd, dt_int32 nit_pcf)
				{
					T sigma_r = fcn_sigma_r_2_sigma_g(sigma_g);
					radius_f = fcn_set_bound(radius_f, 3*this->grid_2d.dR_min(), sigma_r);

					for(auto it=0; it<nit_pcf; it++)
					{
						Pcf_2d<T, Dev>::operator()(M_r, M_s, A, txy, p, sigma_g, bd, this->pcf_r);

						// get maximum position
						R_2d<T> r_c = this->Fit_Ellipt_Gauss_2d.fd_max_peak_pos(this->pcf_r);

						if ((it==nit_pcf-1) & b_fit)
						{
							r_c = this->Fit_Ellipt_Gauss_2d.fit_peak_pos(this->pcf_r, r_c, sigma_r, radius_f);
						}

						txy += -(r_c - this->grid_2d.bs_h());
					}

					return txy;
				}
			protected:
				dt_bool b_fit;
		};

		/******************************** Corrrect shift ***************************************/
		template <class T, eDev Dev>
		class Crt_Sft_1d: public Fd_Sft_1d<T, Dev>
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				Crt_Sft_1d(): Fd_Sft_1d<T, Dev>() {}

				Crt_Sft_1d(FFT<T, Dev> *fft_1d_i, Grid_1d<T>& grid_1d_i): 
				Fd_Sft_1d<T, Dev>(fft_1d_i, grid_1d_i) {}

				T operator()(TVctr_r &M_r_i, TVctr_r &M_s_io, T p, T sigma_g, 
				Region _Rect_1d<T> bd, dt_int32 nit_pcf)
				{
					T dx = 0;
					dx = Fd_Sft_1d<T, Dev>::operator()(M_r_i, M_s_io, p, sigma_g, dx, bd, nit_pcf);
					this->sft_1d(dx, M_s_io);

					return dx;
				}
		};

		template <class T, eDev Dev>
		class Crt_Sft_2d_BC: public Fd_Sft_2d_BC<T, Dev>
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;
				using TVctr_rh = Vctr<T, edev_cpu>;
				using TVctr_ch = Vctr<T_c, edev_cpu>;

				Crt_Sft_2d_BC(): Fd_Sft_2d_BC<T, Dev>() {}

				Crt_Sft_2d_BC(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i): 
				Fd_Sft_2d_BC<T, Dev>(stream_i, fft_2d_i, grid_2d_i) {}

				TVctr_rh operator()(TVctr_r &M_r_i, TVctr_r &M_s_io, T p, T sigma_g, 
				Region _Rect_2d<T> bd, dt_int32 nit_pcf)
				{
					TVctr_rh dr(this->grid_2d.ny, T(0));
					dr = Fd_Sft_2d_BC<T, Dev>::operator()(M_r_i, M_s_io, p, sigma_g, dr, bd, nit_pcf);
					this->sft_2d_bc(1, dr, M_s_io);

					return dr;
				}
		};

		template <class T, eDev Dev>
		class Crt_Sft_2d: public Fd_Sft_2d<T, Dev>
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				Crt_Sft_2d(): Fd_Sft_2d<T, Dev>() {}

				Crt_Sft_2d(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i): 
					Fd_Sft_2d<T, Dev>(stream_i, fft_2d_i, grid_2d_i) {}

				R_2d<T> operator()(TVctr_r &M_r_i, TVctr_r &M_s_io, T p, T sigma_g, 
				Region _Rect_2d<T> bd, dt_int32 nit_pcf)
				{
					R_2d<T> dr(0, 0);
					dr = Fd_Sft_2d<T, Dev>::operator()(M_r_i, M_s_io, p, sigma_g, dr, bd, nit_pcf);
					this->sft_2d(dr, M_s_io);

					return dr;
				}
		};

		template <class T, eDev Dev>
		class Crt_Tr_2d: public Fd_Tr_2d<T, Dev>
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				Crt_Tr_2d(): Fd_Tr_2d<T, Dev>() {}

				Crt_Tr_2d(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i, 
				eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0, dt_bool b_fit_i = true): 
				Fd_Tr_2d<T, Dev>(stream_i, fft_2d_i, grid_2d_i, bg_opt_i, bg_i, b_fit_i) {}

				R_2d<T> operator()(TVctr_r &M_r, TVctr_r &M_s, Mx_2x2<T> A, R_2d<T> txy, 
				T p, T sigma_g, T radius, Region _Rect_2d<T> bd, dt_int32 nit_pcf, TVctr_r& mx_o)
				{
					txy = Fd_Tr_2d<T, Dev>::operator()(M_r, M_s, A, txy, p, sigma_g, radius, bd, nit_pcf);
					this->af_tr(M_s, A, txy, mx_o);

					return txy;
				}
		};

		/*********************************** calculate Chi^2 ***********************************/
 		template <class T, eDev Dev>
		class Scy_Shx
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Scy_Shx():stream(nullptr) {}

				Scy_Shx(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i)
				{
					set_in_data(stream_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i)
				{
					stream = stream_i;
					grid_2d = grid_2d_i;
				}

				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				operator()(TVctr& mx_i, R_2d<T> af, TVctr& mx_o)
				{
					if (fcn_is_zero<T>(af.x) && fcn_is_equal<T>(af.y, T(1)))
					{
						mx_o = mx_i;
					}

					TVctr M_o_t;
					TVctr *pM_o = &mx_o;

					if (mx_i.data() == mx_o.data())
					{
						M_o_t.resize(mx_i.size());
						pM_o = &M_o_t;
					}

					// calculate background
					T bg = fcn_mean(*stream, mx_i);

					Grid_2d<T> grid_2d_o = grid_2d;

					stream->set_n_stream_act(grid_2d_o.nx);
					stream->set_grid(grid_2d_o.nx, grid_2d_o.ny);
					stream->exec_2d(cgpu_detail::shx_scy<Grid_2d<T>, TVctr>, grid_2d, mx_i, af.x, af.y, bg, grid_2d_o, *pM_o);
				
					if (mx_i.data() == mx_o.data())
					{
						mx_o.assign(pM_o->begin(), pM_o->end());
					}
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				operator()(TVctr& mx_i, R_2d<T> af, TVctr& mx_o)
				{
					if (fcn_is_zero<T>(af.x) && fcn_is_equal<T>(af.y, T(1)))
					{
						mx_o = mx_i;
					}

					TVctr M_o_t;
					TVctr *pM_o = &mx_o;

					if (mx_i.data() == mx_o.data())
					{
						M_o_t.resize(mx_i.size());
						pM_o = &M_o_t;
					}

					// calculate background
					T bg = fcn_mean(*stream, mx_i);

					Grid_2d<T> grid_2d_o = grid_2d;

					auto d_grid_blk = grid_2d.d_grid_blk();
					gpu_detail::shx_scy<Grid_2d<T>, typename TVctr::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, mx_i, af.x, af.y, bg, grid_2d_o, *pM_o);
				
					if (mx_i.data() == mx_o.data())
					{
						mx_o.assign(pM_o->begin(), pM_o->end());
					}
				}
			#endif

			protected:
				Grid_2d<T> grid_2d;
				Stream<Dev> *stream;
		};

		template <class T, eDev Dev>
		class Chi2_Pcf_Sft_Scy_Shx_2d: public Crt_Sft_2d<T, Dev>
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				Chi2_Pcf_Sft_Scy_Shx_2d(): Crt_Sft_2d<T, Dev>() {}

				Chi2_Pcf_Sft_Scy_Shx_2d(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i): Crt_Sft_2d<T, Dev>()
				{
					set_in_data(stream_i, fft_2d_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					Crt_Sft_2d<T, Dev>::set_in_data(stream_i, fft_2d_i, grid_2d_i);
					shx_scy.set_in_data(this->stream, this->grid_2d);
				}

				T operator()(TVctr_r &M_r_i, TVctr_r &M_s_i, T p, T sigma_g, 
				R_2d<T> af, Region _Rect_2d<T> bd, dt_int32 nit_pcf, R_2d<T>& ds)
				{
					TVctr_r M(M_r_i.size());
					shx_scy(M_s_i, af, M);

					// correct shift and set borders
					ds = Crt_Sft_2d<T, Dev>::operator()(M_r_i, M, p, sigma_g, bd, nit_pcf);

					// pcf
					TVctr_r &pcf = M;
					Pcf_2d<T, Dev>::operator()(M_r_i, M, p, sigma_g, bd, pcf);

					// cost function
					return fcn_mean(*(this->stream), pcf);
				}

				T operator()(TVctr_r &M_r_i, TVctr_r &M_s_i, T p, T sigma_g, 
				R_2d<T> af, Region _Rect_2d<T> bd, dt_int32 nit_pcf, R_2d<T>& ds, Vctr<T, edev_cpu>& coef)
				{
					TVctr_r M(M_r_i.size());
					shx_scy(M_s_i, af, M);

					// correct shift and set borders
					ds = Crt_Sft_2d<T, Dev>::operator()(M_r_i, M, p, sigma_g, bd, nit_pcf);

					// pcf
					TVctr_r &pcf = M;
					Pcf_2d<T, Dev>::operator()(M_r_i, M, p, sigma_g, bd, pcf);

					// cost function
					auto chi2 = fcn_mean(*(this->stream), pcf);

					// get maximum position
					R_2d<T> r_c = this->Fit_Ellipt_Gauss_2d.fd_max_peak_pos(pcf);

					// fitting
					T sigma_r = mt::fcn_sigma_r_2_sigma_g(sigma_g);
					T radius = ::fmax(3*this->grid_2d.dR_min(), 1.5*sigma_r);
					coef = this->Fit_Ellipt_Gauss_2d.fit(pcf, r_c, sigma_r, radius);
					coef[0] = ds.x + this->grid_2d.bs_x_h();
					coef[1] = ds.y + this->grid_2d.bs_y_h();

					return chi2;
				}

			protected:
				Scy_Shx<T, Dev> shx_scy;
		};
 	
		template <class T, eDev Dev>
		class Chi2_Pcf_2d: public Fd_Tr_2d<T, Dev>
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				Chi2_Pcf_2d(): Fd_Tr_2d<T, Dev>() {}

				Chi2_Pcf_2d(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, 
				Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0):
				Fd_Tr_2d<T, Dev>(stream_i, fft_2d_i, grid_2d_i, bg_opt_i, bg_i) {}

				T operator()(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, Region _Rect_2d<T> bd)
				{
					return this->chi_2(M_r, M_s, p, sigma_g, bd);
				}

				T operator()(TVctr_r &M_r, TVctr_r &M_s, Mx_2x2<T> A, R_2d<T> txy, 
				T p, T sigma_g, Region _Rect_2d<T> bd)
				{
					return this->chi_2(M_r, M_s, A, txy, p, sigma_g, bd);
				}

				T operator()(TVctr_r &M_r, TVctr_r &M_s, Mx_2x2<T> A, R_2d<T> txy, 
				T p, T sigma_g, T radius, Region _Rect_2d<T> bd, Vctr<T, edev_cpu>& coef)
				{
					return this->chi_2(M_r, M_s, A, txy, p, sigma_g, radius, bd, coef);
				}

				R_2d<T> find_txy(TVctr_r &M_r, TVctr_r &M_s, Mx_2x2<T> A, R_2d<T> txy, 
				T p, T sigma_g, T radius, Region _Rect_2d<T> bd, dt_int32 nit_pcf)
				{
					return Fd_Tr_2d<T, Dev>::operator()(M_r, M_s, A, txy, p, sigma_g, radius, bd, nit_pcf);
				}
		};

		template <class T, eDev Dev>
		class Fd_Sft_Scy_Shx_2d: public Chi2_Pcf_Sft_Scy_Shx_2d<T, Dev>
		{
			public:
				using T_r = T;
				using T_c = complex<T>;
				using TVctr_r = Vctr<T, Dev>;
				using TVctr_c = Vctr<T_c, Dev>;

				Fd_Sft_Scy_Shx_2d():Chi2_Pcf_Sft_Scy_Shx_2d<T, Dev>(), nit_pcf(2) {}

				Fd_Sft_Scy_Shx_2d(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i): 
				Chi2_Pcf_Sft_Scy_Shx_2d<T, Dev>(stream_i, fft_2d_i, grid_2d_i), nit_pcf(2) {}

				void operator()(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, 
				Vctr<Atp_1<T>, edev_cpu>& spx, Region _Rect_2d<T> bd, dt_int32 nit_nm)
				{
					dt_int32 ic = 1;

					T d_pix_min = this->grid_2d.dg_min();
					T d_min = 0.05*d_pix_min;
					T alpha = 1.0;
					T beta = 0.5;
					T gamma = 2.0;

					for(auto it=0; it<nit_nm; it++)
					{
						sort(spx);
						T d_l = max_length(spx);
						if (d_l<d_min)
						{
							break;
						}

						if (d_l<d_pix_min)
						{
							if (ic==1)
							{
								orthogonal_simplex(M_r, M_s, p, sigma_g, spx, bd, nit_pcf);
								ic++;
							}
						}

						auto x_m = best_centroid(spx);

						auto x_b = spx[0].f;
						auto chi2_b = spx[0].chi2;

						auto x_w = spx[2].f;
						auto chi2_w = spx[2].chi2;

						auto x_r = x_m+alpha*(x_m-x_w);
						R_2d<T> ds_r(0, 0);
						auto chi2_r = chi2_pcf(M_r, M_s, p, sigma_g, x_r, bd, nit_pcf, ds_r);

						if (chi2_r<chi2_b)
						{
							auto x_e = x_m + gamma*(x_r-x_m);
							R_2d<T> ds_e(0, 0);
							auto chi2_e = chi2_pcf(M_r, M_s, p, sigma_g, x_e, bd, nit_pcf, ds_e);
							if (chi2_e<chi2_r)
							{
								spx[2] = Atp_1<T>(x_e, ds_e, chi2_e);
							}
							else
							{
								spx[2] = Atp_1<T>(x_r, ds_r, chi2_r);
							}
						}
						else if (chi2_r<spx[1].chi2)
						{
							spx[2] = Atp_1<T>(x_r, ds_r, chi2_r);
						}
						else if (chi2_r<chi2_w)
						{
							auto x_rc = x_m + beta*(x_r-x_m);
							R_2d<T> ds_rc(0, 0);
							auto chi2_rc = chi2_pcf(M_r, M_s, p, sigma_g, x_rc, bd, nit_pcf, ds_rc);
							if (chi2_rc<chi2_r)
							{
								spx[2] = Atp_1<T>(x_rc, ds_rc, chi2_rc);
							}
							else
							{
								contract_simplex(M_r, M_s, p, sigma_g, spx, bd, nit_pcf);
							}
						}
						else
						{
							auto x_wc = x_m + beta*(x_w-x_m);
							R_2d<T> ds_wc(0, 0);
							auto chi2_wc = chi2_pcf(M_r, M_s, p, sigma_g, x_wc, bd, nit_pcf, ds_wc);
							if (chi2_wc<chi2_w)
							{
								spx[2] = Atp_1<T>(x_wc, ds_wc, chi2_wc);
							}
							else
							{
								contract_simplex(M_r, M_s, p, sigma_g, spx, bd, nit_pcf);
							}
						}
					}

				} 

				void operator()(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, 
				R_2d<T>& af, Region _Rect_2d<T> bd, R_2d<T>& tr, dt_int32 nit_nm)
				{
					auto spx = set_simplex(M_r, M_s, p, sigma_g, af, tr, bd, nit_nm);
					this->operator()(M_r, M_s, p, sigma_g, spx, bd, nit_nm);

					af = spx[0].f;
					tr = spx[0].tr;

					// // shear and scaling
					// shx_scy(M_s, af, M_s);
					// // correct shift
					// tr = Fd_Sft_2d<T, Dev>::operator()(M_r, M_s, p, sigma_g, tr, bd, nit_pcf);
				}

			protected:
				R_2d<T> x_af(const R_2d<T>& af, const R_2d<T>& r)
				{
					return R_2d<T>(r.x+af.x*r.y, af.y*r.y);
				}

				R_2d<T> x_iaf(const R_2d<T>& af, const R_2d<T>& r)
				{
					return R_2d<T>(r.x-af.x*r.y/af.y, r.y/af.y);
				}

				Vctr<Atp_1<T>, edev_cpu> set_simplex(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, 
				R_2d<T> x_0, R_2d<T> ds_0, Region _Rect_2d<T>& bd, dt_int32 nit_nm)
				{
					// global shift
					if (fcn_is_zero(ds_0))
					{
						ds_0 = Fd_Sft_2d<T, Dev>::operator()(M_r, M_s, p, sigma_g, ds_0, bd, nit_pcf);
						// bd.fcn_repl_bdr(ds_0);
					}

					TVctr_r M_s_t = M_s;
					this->sft_2d(x_iaf(x_0, ds_0), M_s_t);

					// determine coefficients
					Vctr<T, edev_cpu> coef(6);
					chi2_pcf(M_r, M_s_t, p, sigma_g, x_0, bd, nit_pcf, ds_0, coef);

					// calculate dd0
					T ff = coef[4]/coef[3];
					ff = (ff<1)?1/ff:ff;
					ff = ::fmax(1.05, ff);

					T dd0 = ::sqrt(1.354e-04*pow(ff-1, 2)+6.622e-4*(ff-1));
					dd0 = ::fmax(dd0, this->grid_2d.dg_min());

					T sin_t = sin(coef[5]);
					T cos_t = cos(coef[5]);

					T sigma_x = pow(cos_t/coef[3], 2)+pow(sin_t/coef[4], 2);
					sigma_x = ::sqrt(1/sigma_x);

					T sigma_y = pow(sin_t/coef[3], 2)+pow(cos_t/coef[4], 2);
					sigma_y = ::sqrt(1/sigma_y);

					T theta = atan2(sigma_y, sigma_x);

					R_2d<T> u(dd0*cos(theta), dd0*sin(theta));
					R_2d<T> v(-u.y, u.x);

					// set simplex
					Vctr<Atp_1<T>, edev_cpu> spx(3);

					spx[0].f = x_0;
					spx[0].chi2 = chi2_pcf(M_r, M_s, p, sigma_g, spx[0].f, bd, nit_pcf, spx[0].tr);

					spx[1].f = x_0+u;
					spx[1].chi2 = chi2_pcf(M_r, M_s, p, sigma_g, spx[1].f, bd, nit_pcf, spx[1].tr);

					spx[2].f = x_0+v;
					spx[2].chi2 = chi2_pcf(M_r, M_s, p, sigma_g, spx[2].f, bd, nit_pcf, spx[2].tr);

					// sort simplex
					sort(spx);

					return spx;
				}

				T chi2_pcf(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, 
				R_2d<T> af, Region _Rect_2d<T>& bd, dt_int32 nit_pcf, R_2d<T>& tr)
				{
					return Chi2_Pcf_Sft_Scy_Shx_2d<T, Dev>::operator()(M_r, M_s, p, sigma_g, af, bd, nit_pcf, tr);
				}

				T chi2_pcf(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, 
				R_2d<T> af, Region _Rect_2d<T>& bd, dt_int32 nit_pcf, R_2d<T>& tr, Vctr<T, edev_cpu>& coef)
				{
					return Chi2_Pcf_Sft_Scy_Shx_2d<T, Dev>::operator()(M_r, M_s, p, sigma_g, af, bd, nit_pcf, tr, coef);
				}

				void sort(Vctr<Atp_1<T>, edev_cpu>& spx)
				{
					std::sort(spx.begin(), spx.end(), [](const Atp_1<T>& x, const Atp_1<T>& y){ return x.chi2<y.chi2; });
				}

				T max_length(Vctr<Atp_1<T>, edev_cpu>& spx)
				{
					auto r0 = spx[0].f;

					R_2d<T> dr = spx[1].f - r0;
					T d_max = dr.norm();
					for(auto i=2; i<spx.size(); i++)
					{
						R_2d<T> dr = spx[i].f - r0;
						d_max = ::fmax(d_max, dr.norm());
					}
					return d_max;
				}

				T min_length(Vctr<Atp_1<T>, edev_cpu>& spx)
				{
					auto r0 = spx[0].f;
					R_2d<T> dr = spx[1].f - r0;
					T d_min = dr.norm();
					for(auto i=2; i<spx.size(); i++)
					{
						R_2d<T> dr = spx[i].f - r0;
						d_min = ::fmin(d_min, dr.norm());
					}
					return d_min;
				}

				R_2d<T> best_centroid(Vctr<Atp_1<T>, edev_cpu>& spx)
				{
					R_2d<T> r_c(0, 0);
					for(auto i=0; i<spx.size()-1; i++)
					{
						r_c += spx[i].f;
					}
					return r_c/T(spx.size()-1);
				}

				void orthogonal_simplex(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, 
				Vctr<Atp_1<T>, edev_cpu>& spx, Region _Rect_2d<T>& bd, dt_int32 nit_pcf)
				{
					auto p12 = spx[1].f-spx[0].f;
					auto p13 = spx[2].f-spx[0].f;
					auto mp12 = p12.norm();
					auto mp13 = p13.norm();
					auto mp_max = ::fmax(mp12, mp13);

					T theta = angle(p12, p13);

					T theta_min = 10;
					T theta_0 = theta_min*c_pi/180;
					T theta_e = c_pi-theta_min*c_pi/180;
	
					dt_bool b_m = (mp12<mp_max/4)||(mp13<mp_max/4);

					if ((theta<theta_0)||(theta>theta_e)||b_m)
					{
			 if (b_m && (mp12<mp13))
			 {
							auto u = p13/norm(p13);
							auto x_o = spx[0].f + mp_max*R_2d<T>(-u.y, u.x);
							R_2d<T> ds_o(0, 0);
							auto chi2_o = chi2_pcf(M_r, M_s, p, sigma_g, x_o, bd, nit_pcf, ds_o);

							spx[1] = Atp_1<T>(x_o, ds_o, chi2_o);
			 }
			 else
			 {
							T m = max_length(spx);
							auto u = p12/norm(p12);
							auto x_o = spx[0].f + mp_max*R_2d<T>(-u.y, u.x);
							R_2d<T> ds_o(0, 0);
							auto chi2_o = chi2_pcf(M_r, M_s, p, sigma_g, x_o, bd, nit_pcf, ds_o);

							spx[2] = Atp_1<T>(x_o, ds_o, chi2_o);
			 }

						sort(spx);
					}
				}

				void contract_simplex(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, 
				Vctr<Atp_1<T>, edev_cpu>& spx, Region _Rect_2d<T>& bd, dt_int32 nit_pcf)
				{
					T ff = 0.5;

					auto p12 = spx[1].f-spx[0].f;
					auto p13 = spx[2].f-spx[0].f;
					auto mp12 = p12.norm();
					auto mp13 = p13.norm();

					if (mp12<mp13)
					{
						auto u = p12/mp12;
						auto x_c = spx[0].f + ff*mp13*R_2d<T>(-u.y, u.x);
						R_2d<T> ds_c(0, 0);
						auto chi2_c = chi2_pcf(M_r, M_s, p, sigma_g, x_c, bd, nit_pcf, ds_c);

						spx[2] = Atp_1<T>(x_c, ds_c, chi2_c);
					}
					else
					{
						auto u = p13/mp13;
						auto x_c = spx[0].f + ff*mp12*R_2d<T>(-u.y, u.x);
						R_2d<T> ds_c(0, 0);
						auto chi2_c = chi2_pcf(M_r, M_s, p, sigma_g, x_c, bd, nit_pcf, ds_c);

						spx[1] = Atp_1<T>(x_c, ds_c, chi2_c);
					}

					sort(spx);
				}

				dt_int32 nit_pcf;
		};

		template <class T, eDev Dev>
		class Fd_Tr_Scy_Shx_2d:public Chi2_Pcf_2d<T, Dev>
		{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVctr_r = Vctr<T, Dev>;
			using TVctr_c = Vctr<T_c, Dev>;

			Fd_Tr_Scy_Shx_2d(): Chi2_Pcf_2d<T, Dev>(), 
			nit_pcf(2), c_r(1.0), c_e(2.0), c_c(0.5), c_s(0.5) {}

			Fd_Tr_Scy_Shx_2d(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0):
			Chi2_Pcf_2d<T, Dev>(stream_i, fft_2d_i, grid_2d_i, bg_opt_i, bg_i), 
			nit_pcf(2), c_r(1.0), c_e(2.0), c_c(0.5), c_s(0.5) {}

			void operator()(TVctr_r &M_r, TVctr_r &M_s, T& scy, T& shx, R_2d<T>& txy, 
			T p, T sigma_g, T radius, Region _Rect_2d<T> bd, dt_int32 nit_nm)
			{
				Simplex spx = spx_0(M_r, M_s, scy, shx, txy, p, sigma_g, radius, bd);

				dt_int32 ic = 1;

				T d_pix_min = this->grid_2d.dg_min();
				T d_min = 0.05*d_pix_min;

				for(auto it = 0; it<nit_nm; it++)
				{
					T d_l = spx.max_length();
					if (d_l<d_min)
					{
						break;
					}

					if ((d_l<d_pix_min)&(ic == 1))
					{
						orthogonal_simplex(M_r, M_s, p, sigma_g, radius, bd, spx);
						ic++;
					}

					auto x_m = spx.centroid();

					auto x_b = spx[0].f;
					auto t_b = spx[0].tr;
					auto chi2_b = spx[0].chi2;

					auto x_w = spx[2].f;
					auto chi2_w = spx[2].chi2;

					auto x_r = x_m + c_r*(x_m - x_w);
					auto t_r = find_txy(M_r, M_s, x_r, t_b, p, sigma_g, radius, bd);
					auto chi2_r = chi2_pcf(M_r, M_s, x_r, t_r, p, sigma_g, bd);

					if (chi2_r<spx[1].chi2)
					{
						if (chi2_b<chi2_r)
						{
							spx[2] = Atp_1<T>(x_r, t_r, chi2_r);
						}
						else
						{
							auto x_e = x_m + c_e*(x_r - x_m);
							auto t_e = find_txy(M_r, M_s, x_e, t_b, p, sigma_g, radius, bd);
							auto chi2_e = chi2_pcf(M_r, M_s, x_e, t_e, p, sigma_g, bd);

							spx[2] = (chi2_e<chi2_b)?Atp_1<T>(x_e, t_e, chi2_e):Atp_1<T>(x_r, t_r, chi2_r);
						}
					}
					else
					{
						if (chi2_r<chi2_w)
						{
							spx[2] = Atp_1<T>(x_r, t_r, chi2_r);
						}

						auto x_c = x_m + c_c*(x_r - x_m);
						auto t_c = find_txy(M_r, M_s, x_c, t_b, p, sigma_g, radius, bd);
						auto chi2_c = chi2_pcf(M_r, M_s, x_c, t_c, p, sigma_g, bd);

						if (chi2_c<chi2_r)
						{
							spx[2] = Atp_1<T>(x_c, t_c, chi2_c);
						}
						else
						{
							shrink_simplex(M_r, M_s, p, sigma_g, radius, bd, spx);
						}
					}

					spx.sort();
				}

				shx = spx[0].f.x;
				scy = spx[0].f.y;
				txy = spx[0].tr;
			}

		protected:
			const dt_int32 nit_pcf;

			const T c_r;
			const T c_e;
			const T c_c;
			const T c_s;

			struct Simplex
			{
				Simplex() {}

				Simplex(dt_int32 new_size)
				{
					resize(new_size);
				}

				dt_int32 size() const
				{
					return m_spx.size();
				}

				void resize(dt_int32 new_size)
				{
					m_spx.resize(new_size);
				}

				Atp_1<T>& operator[](const dt_int32 idx) { return m_spx[idx]; }

				const Atp_1<T>& operator[](const dt_int32 idx) const { return m_spx[idx]; }

				void sort()
				{
					std::sort(m_spx.begin(), m_spx.end(), [](const Atp_1<T>& x, const Atp_1<T>& y) { return x.chi2 < y.chi2; });
				}

				T max_length()
				{
					auto f0 = m_spx[0].f;
					R_2d<T> dr = m_spx[1].f - f0;
					T d_max = dr.norm();
					for(auto ik = 2; ik < m_spx.size(); ik++)
					{
						R_2d<T> dr = m_spx[ik].f - f0;
						d_max = ::fmax(d_max, dr.norm());
					}
					return d_max;
				}

				T min_length()
				{
					auto f0 = m_spx[0].f;
					R_2d<T> dr = m_spx[1].f - f0;
					T d_min = dr.norm();
					for(auto ik = 2; ik < m_spx.size(); ik++)
					{
						R_2d<T> dr = m_spx[ik].f - f0;
						d_min = ::fmin(d_min, dr.norm());
					}
					return d_min;
				}

				R_2d<T> centroid()
				{
					R_2d<T> r_c(0, 0);
					dt_int32 n_spx = m_spx.size()-1;
					for(auto ik = 0; ik<n_spx; ik++)
					{
						r_c += m_spx[ik].f;
					}
					return r_c/T(n_spx);
				}

				private:
					Vctr<Atp_1<T>, edev_cpu> m_spx;
			};

			R_2d<T> find_txy(TVctr_r &M_r, TVctr_r &M_s, R_2d<T> f, R_2d<T> txy, T p, T sigma_g, T radius, Region _Rect_2d<T>& bd)
			{
				Mx_2x2<T> A(1, 0, f.x, f.y);
				return Chi2_Pcf_2d<T, Dev>::find_txy(M_r, M_s, A, txy, p, sigma_g, radius, bd, nit_pcf);
			}

			T chi2_pcf(TVctr_r &M_r, TVctr_r &M_s, R_2d<T> f, R_2d<T> txy, T p, T sigma_g, Region _Rect_2d<T>& bd)
			{
				Mx_2x2<T> A(1, 0, f.x, f.y);
				return Chi2_Pcf_2d<T, Dev>::operator()(M_r, M_s, A, txy, p, sigma_g, bd);
			}

			T chi2_pcf(TVctr_r &M_r, TVctr_r &M_s, R_2d<T> f, R_2d<T> txy, T p, T sigma_g, T radius, 
			Region _Rect_2d<T>& bd, Vctr<T, edev_cpu>& coef)
			{
				Mx_2x2<T> A(1, 0, f.x, f.y);
				return Chi2_Pcf_2d<T, Dev>::operator()(M_r, M_s, A, txy, p, sigma_g, radius, bd, coef);
			}

			Simplex spx_0(TVctr_r &M_r, TVctr_r &M_s, T scy, T shx, R_2d<T> txy, 
			T p, T sigma_g, T radius, Region _Rect_2d<T> bd)
			{
				// global shift
				Vctr<T, edev_cpu> coef(6);

				auto x_0 = R_2d<T>(shx, scy);
				auto t_0 = find_txy(M_r, M_s, x_0, txy, p, sigma_g, radius, bd);
				auto chi2_0 = chi2_pcf(M_r, M_s, x_0, t_0, p, sigma_g, radius, bd, coef);

				// bd.fcn_repl_bdr(t_0);

				// calculate dd0
				T ff = coef[4]/coef[3];
				ff = (ff<1)?1/ff:ff;
				ff = ::fmax(1.05, ff);

				T dd0 = ::sqrt(1.354e-04*pow(ff - 1, 2) + 6.622e-4*(ff - 1));
				dd0 = ::fmax(dd0, this->grid_2d.dg_min());

				T sin_t = sin(coef[5]);
				T cos_t = cos(coef[5]);

				T sigma_x = pow(cos_t/coef[3], 2) + pow(sin_t/coef[4], 2);
				sigma_x = ::sqrt(1/sigma_x);

				T sigma_y = pow(sin_t/coef[3], 2) + pow(cos_t/coef[4], 2);
				sigma_y = ::sqrt(1/sigma_y);

				T theta = atan2(sigma_y, sigma_x);

				R_2d<T> u(dd0*cos(theta), dd0*sin(theta));
				R_2d<T> v(-u.y, u.x);

				// set simplex
				Simplex spx(3);

				spx[0].f = x_0;
				spx[0].tr = t_0;
				spx[0].chi2 = chi2_0;

				spx[1].f = x_0 + u;
				spx[1].tr = find_txy(M_r, M_s, spx[1].f, t_0, p, sigma_g, radius, bd);
				spx[1].chi2 = chi2_pcf(M_r, M_s, spx[1].f, spx[1].tr, p, sigma_g, bd);

				spx[2].f = x_0 + v;
				spx[2].tr = find_txy(M_r, M_s, spx[2].f, t_0, p, sigma_g, radius, bd);
				spx[2].chi2 = chi2_pcf(M_r, M_s, spx[2].f, spx[2].tr, p, sigma_g, bd);

				spx.sort();

				return spx;
			}

			void orthogonal_simplex(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, T radius, 
			Region _Rect_2d<T>& bd, Simplex &spx)
			{
				auto x_b = spx[0].f;
				auto t_b = spx[0].tr;

				auto p12 = spx[1].f - x_b;
				auto p13 = spx[2].f - x_b;
				auto mp12 = p12.norm();
				auto mp13 = p13.norm();
				auto mp_max = ::fmax(mp12, mp13);

				T theta = angle(p12, p13);

				T theta_min = 10;
				T theta_0 = theta_min*c_pi/180;
				T theta_e = c_pi - theta_min*c_pi/180;

				dt_bool b_m = (mp12<mp_max/4) || (mp13<mp_max/4);

				if ((theta<theta_0) || (theta>theta_e) || b_m)
				{
					if (b_m && (mp12<mp13))
					{
						auto u = normalize(p13);
						auto x_o = x_b + mp_max*R_2d<T>(-u.y, u.x);
						auto t_o = find_txy(M_r, M_s, x_o, t_b, p, sigma_g, radius, bd);
						auto chi2_o = chi2_pcf(M_r, M_s, x_o, t_o, p, sigma_g, bd);

						spx[1] = Atp_1<T>(x_o, t_o, chi2_o);
					}
					else
					{
						auto u = normalize(p12);
						auto x_o = x_b + mp_max*R_2d<T>(-u.y, u.x);
						auto t_o = find_txy(M_r, M_s, x_o, t_b, p, sigma_g, radius, bd);
						auto chi2_o = chi2_pcf(M_r, M_s, x_o, t_o, p, sigma_g, bd);

						spx[2] = Atp_1<T>(x_o, t_o, chi2_o);
					}

					spx.sort();
				}
			}

			void shrink_simplex(TVctr_r &M_r, TVctr_r &M_s, T p, T sigma_g, T radius, 
			Region _Rect_2d<T>& bd, Simplex &spx)
			{
				auto x_b = spx[0].f;
				auto t_b = spx[0].tr;
				for(auto ik = 1; ik < spx.size(); ik++)
				{
					auto x_s = x_b + c_s*(spx[ik].f-x_b);
					auto t_s = find_txy(M_r, M_s, x_s, t_b, p, sigma_g, radius, bd);
					auto chi2_s = chi2_pcf(M_r, M_s, x_s, t_s, p, sigma_g, bd);

					spx[ik] = Atp_1<T>(x_s, t_s, chi2_s);
				}
			}
		};

		template <class T, eDev Dev>
		class Fd_Tr_Rot_2d:public Chi2_Pcf_2d<T, Dev>
		{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVctr_r = Vctr<T, Dev>;
			using TVctr_c = Vctr<T_c, Dev>;

			Fd_Tr_Rot_2d(): Chi2_Pcf_2d<T, Dev>(), 
			nit_pcf(2), c_r(1.0), c_e(2.0), c_c(0.5), c_s(0.5), p_c(0, 0) {}

			Fd_Tr_Rot_2d(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0):
			Chi2_Pcf_2d<T, Dev>(stream_i, fft_2d_i, grid_2d_i, bg_opt_i, bg_i), 
			nit_pcf(2), c_r(1.0), c_e(2.0), c_c(0.5), c_s(0.5), p_c(grid_2d_i.rx_c(), grid_2d_i.ry_c()) {}

			void operator()(TVctr_r &M_r, TVctr_r &M_s, T& theta, R_2d<T>& txy, 
			T p, T sigma_g, T radius, Region _Rect_2d<T> bd, dt_int32 nit_nm)
			{
				Simplex spx = spx_0(M_r, M_s, theta, txy, p, sigma_g, radius, bd);

				T d_min = 0.025*c_deg_2_rad;

				for(auto it = 0; it<nit_nm; it++)
				{
					if (spx.max_length()<d_min)
					{
						break;
					}

					auto x_m = spx.centroid();

					auto x_b = spx[0].theta;
					auto t_b = spx[0].tr;
					auto chi2_b = spx[0].chi2;

					auto x_w = spx[1].theta;
					auto chi2_w = spx[1].chi2;

					auto x_r = x_m + c_r*(x_m - x_w);
					auto t_r = find_txy(M_r, M_s, x_r, t_b, p, sigma_g, radius, bd);
					auto chi2_r = chi2_pcf(M_r, M_s, x_r, t_r, p, sigma_g, bd);

					if (chi2_r<spx[1].chi2)
					{
						if (chi2_b<chi2_r)
						{
							spx[1] = Atp_2<T>(x_r, t_r, chi2_r);
						}
						else
						{
							auto x_e = x_m + c_e*(x_r - x_m);
							auto t_e = find_txy(M_r, M_s, x_e, t_b, p, sigma_g, radius, bd);
							auto chi2_e = chi2_pcf(M_r, M_s, x_e, t_e, p, sigma_g, bd);

							spx[1] = (chi2_e<chi2_b)?Atp_2<T>(x_e, t_e, chi2_e):Atp_2<T>(x_r, t_r, chi2_r);
						}
					}
					else
					{
						if (chi2_r<chi2_w)
						{
							spx[1] = Atp_2<T>(x_r, t_r, chi2_r);
						}

						auto x_c = x_m + c_c*(x_r - x_m);
						auto t_c = find_txy(M_r, M_s, x_c, t_b, p, sigma_g, radius, bd);
						auto chi2_c = chi2_pcf(M_r, M_s, x_c, t_c, p, sigma_g, bd);

						if (chi2_c<chi2_r)
						{
							spx[1] = Atp_2<T>(x_c, t_c, chi2_c);
						}
					}

					spx.sort();
				}

				theta = spx[0].theta;
				txy = spx[0].tr;
			}

		protected:
			const dt_int32 nit_pcf;

			const T c_r;
			const T c_e;
			const T c_c;
			const T c_s;
			const R_2d<T> p_c;

			struct Simplex
			{
				Simplex() {}

				Simplex(dt_int32 new_size)
				{
					resize(new_size);
				}

				dt_int32 size() const
				{
					return m_spx.size();
				}

				void resize(dt_int32 new_size)
				{
					m_spx.resize(new_size);
				}

				Atp_2<T>& operator[](const dt_int32 idx) { return m_spx[idx]; }

				const Atp_2<T>& operator[](const dt_int32 idx) const { return m_spx[idx]; }

				void sort()
				{
					std::sort(m_spx.begin(), m_spx.end(), [](const Atp_2<T>& x, const Atp_2<T>& y) { return x.chi2 < y.chi2; });
				}

				T max_length()
				{
					return ::fabs(m_spx[1].theta - m_spx[0].theta);
				}

				T min_length()
				{
					return ::fabs(m_spx[1].theta - m_spx[0].theta);
				}

				T centroid()
				{
					return m_spx[0].theta;
				}

				private:
				Vctr<Atp_2<T>, edev_cpu> m_spx;
			};

			R_2d<T> find_txy(TVctr_r &M_r, TVctr_r &M_s, T theta, R_2d<T> txy, T p, T sigma_g, T radius, Region _Rect_2d<T>& bd)
			{
				auto A = fcn_rot_mx_2d(theta);
				txy = p_c - A*p_c + txy;

				return Chi2_Pcf_2d<T, Dev>::find_txy(M_r, M_s, A, txy, p, sigma_g, radius, bd, nit_pcf);
			}

			T chi2_pcf(TVctr_r &M_r, TVctr_r &M_s, T theta, R_2d<T> txy, T p, T sigma_g, Region _Rect_2d<T>& bd)
			{
				auto A = mt::fcn_rot_mx_2d(theta);
				txy = p_c - A*p_c + txy;

				return Chi2_Pcf_2d<T, Dev>::operator()(M_r, M_s, A, txy, p, sigma_g, bd);
			}

			Simplex spx_0(TVctr_r &M_r, TVctr_r &M_s, T theta, R_2d<T> txy, 
			T p, T sigma_g, T radius, Region _Rect_2d<T> bd)
			{
				auto x_0 = theta;
				auto t_0 = find_txy(M_r, M_s, x_0, txy, p, sigma_g, radius, bd);
				auto chi2_0 = chi2_pcf(M_r, M_s, x_0, t_0, p, sigma_g, bd);

				// set simplex
				Simplex spx(2);

				spx[0].theta = x_0;
				spx[0].tr = t_0;
				spx[0].chi2 = chi2_0;

				spx[1].theta = x_0 + 1.0*c_deg_2_rad;
				spx[1].tr = find_txy(M_r, M_s, spx[1].theta, t_0, p, sigma_g, radius, bd);
				spx[1].chi2 = chi2_pcf(M_r, M_s, spx[1].theta, spx[1].tr, p, sigma_g, bd);
				spx.sort();

				return spx;
			}
		};

		/********************************* calculate d_phi *************************************/
		template <class T, eDev Dev>
		class Opt_Flow
		{
			public:
				using T_r = T;
				using TVctr = Vctr<T, Dev>;

				static const eDev device = Dev;

				Opt_Flow(): stream(nullptr) {}

				Opt_Flow(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					set_in_data(stream_i, fft_2d_i, grid_2d_i);
				}

				inline
				void set_in_data(Stream<Dev> *stream_i, FFT<T, Dev> *fft_2d_i, Grid_2d<T>& grid_2d_i)
				{
					stream = stream_i;
					grid_2d = grid_2d_i;

					intrpl_rg_2d.set_in_data(stream, grid_2d, grid_2d);

					gauss_cv_2d.set_in_data(stream, fft_2d_i, grid_2d);

					v_x.resize(grid_2d.size());
					v_y.resize(grid_2d.size());

					Rx.resize(grid_2d.size());
					Ry.resize(grid_2d.size());
					M.resize(grid_2d.size());

				}

				void operator()(TVctr& M_s, TVctr& M_m, T alpha, T sigma, dt_int32 n_iter, TVctr& v_xt, TVctr& v_yt)
				{
					// create rectangular grid
					set_regular_grid(Rx, Ry);

					// set initial optical flow
					v_x = v_xt;
					v_y = v_yt;

					for(auto iter = 0; iter < n_iter; iter++)
					{
						// create new grid
						mt::add(*stream, v_x, Rx);
						mt::add(*stream, v_y, Ry);

						// resample distored image in a new grid
						intrpl_rg_2d(M_m, Rx, Ry, M);

						// calculate optical flow
						fcn_opt_flow(M_s, M, alpha, v_x, v_y);

						// regularization based on convolution
						if (fcn_is_nzero(sigma))
						{
							gauss_cv_2d(sigma, v_x);
							gauss_cv_2d(sigma, v_y);
						}

						// add optical flow
						mt::add(*stream, v_x, v_xt);
						mt::add(*stream, v_y, v_yt);

					}
				}

				void set_fft_plan()
				{
					gauss_cv_2d.set_fft_plan();
				}

				void cleanup()
				{
					gauss_cv_2d.cleanup();
				}

			protected:

				void set_regular_grid(TVctr& Rx, TVctr& Ry)
				{
					Vctr<T, edev_cpu> Rx_h;
					Vctr<T, edev_cpu> Ry_h;

					Rx_h.reserve(grid_2d.size());
					Ry_h.reserve(grid_2d.size());

					for(auto ix = 0; ix < grid_2d.nx; ix++)
					{
						for(auto iy = 0; iy < grid_2d.ny; iy++)
						{
							Rx_h.push_back(grid_2d.rx(ix));
							Ry_h.push_back(grid_2d.ry(iy));
						}
					}

					thrust::copy(Rx_h.begin(), Rx_h.end(), Rx.begin());
					thrust::copy(Ry_h.begin(), Ry_h.end(), Ry.begin());
				}

				/***************************************** cpu *****************************************/
				template <eDev devn = Dev>
				enable_if_edev_cpu<devn, void>
				fcn_opt_flow(TVctr& M_s, TVctr& M_m, T alpha, TVctr& v_x, TVctr& v_y)
				{
					stream->set_n_stream_act(grid_2d.nx);
					stream->set_grid(grid_2d.nx, grid_2d.ny);
					stream->exec_2d(cgpu_detail::fcn_opt_flow<Grid_2d<T>, TVctr>, grid_2d, M_s, M_m, alpha, v_x, v_y);
				}

				/**********************Device**********************/
			#ifdef __CUDACC__
				template <eDev devn = Dev>
				enable_if_edev_gpu<devn, void>
				fcn_opt_flow(TVctr& M_s, TVctr& M_m, T alpha, TVctr& v_x, TVctr& v_y)
				{
					auto d_grid_blk = grid_2d.d_grid_blk();
					gpu_detail::fcn_opt_flow<Grid_2d<T>, typename TVctr::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, M_s, M_m, alpha, v_x, v_y);
				}
			#endif

				Stream<Dev> *stream;
				Grid_2d<T> grid_2d;

				Interp_rn_2d<T, Dev> intrpl_rg_2d;
				Gauss_Cv_2d<T, Dev> gauss_cv_2d;

				TVctr v_x;
				TVctr v_y;

				TVctr Rx;
				TVctr Ry;
				TVctr M;
		};

	}

#endif