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

#ifndef CGPU_CLASSES_H
#define CGPU_CLASSES_H

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "fft.cuh"
#include "cgpu_rand.cuh"

#ifdef __CUDACC__
	#include <cuda.h>
	#include <cuda_runtime.h>
	#include <cufft.h>
#endif

#include "cgpu_fcns.cuh"
#include "cpu_fcns.hpp"

#ifdef __CUDACC__
	#include "gpu_fcns.cuh"
#endif
namespace mt
{
	/**************** Gaussian Conv ***************/
	template <class T, eDevice dev>
	class Gauss_Cv_1d
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			static const eDevice device = dev;

			Gauss_Cv_1d():fft_1d(nullptr){}

			Gauss_Cv_1d(FFT<T, dev> *fft_1d_i, Grid_1d<T> &grid_1d_i)
			{
				set_input_data(fft_1d_i, grid_1d_i);
			}

			inline
			void set_input_data(FFT<T, dev> *fft_1d_i, Grid_1d<T> &grid_1d_i)
			{
				fft_1d = fft_1d_i;
				grid_1d = grid_1d_i;
			}

			void set_fft_plan()
			{
				fft_1d->create_plan_1d(grid_1d.nx, 1);
			}

			void operator()(T sigma_r, TVector_c &Im)
			{
				fft1_shift(grid_1d, Im);
				fft_1d->forward(Im);

				gauss_cv_1d(sigma_r, Im);

				fft_1d->inverse(Im);
				fft1_shift(grid_1d, Im);
			}

			void operator()(T sigma_r, TVector_r &Im)
			{
				TVector_c Im_c(Im.begin(), Im.end());

				this->operator()(sigma_r, Im_c);

				assign_real(Im_c, Im);
			}

			void cleanup()
			{
				fft_1d->cleanup();
			}
	protected:
			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			gauss_cv_1d(T sigma_r, TVector_c &Im)
			{
				auto alpha = 2*c_Pi2*sigma_r*sigma_r;

				for (auto ix = 0; ix < grid_1d.nx; ix++)
				{
					host_device_detail::gauss_cv_1d<Grid_1d<T>, TVector_c>(ix, grid_1d, alpha, Im);
				}
			}

			/**********************Device**********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			gauss_cv_1d(T sigma_r, TVector_c &Im)
			{
				auto alpha = 2*c_Pi2*sigma_r*sigma_r;

				auto grid_bt = grid_1d.cuda_grid();
				device_detail::gauss_cv_1d<Grid_1d<T>, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_1d, alpha, Im);
			}
		#endif

			FFT<T, dev> *fft_1d;
			Grid_1d<T> grid_1d;
	};

	template <class T, eDevice dev>
	class Gauss_Cv_2d_BC
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			static const eDevice device = dev;

			Gauss_Cv_2d_BC():stream(nullptr), fft_2d(nullptr){}

			Gauss_Cv_2d_BC(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				set_input_data(stream_i, fft_2d_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				stream = stream_i;
				fft_2d = fft_2d_i;
				grid_2d = grid_2d_i;
			}

			void set_fft_plan()
			{
				fft_2d->create_plan_1d_batch(grid_2d.ny, grid_2d.nx, stream->size());
			}

			void operator()(T sigma_r, TVector_c &Im)
			{
				fft2_sft_bc(*stream, grid_2d, Im);
				fft_2d->forward(Im);

				gauss_cv_2d_bc(sigma_r, Im);

				fft_2d->inverse(Im);
				fft2_sft_bc(*stream, grid_2d, Im);
			}

			void operator()(T sigma_r, TVector_r &Im)
			{
				TVector_c Im_c(Im.begin(), Im.end());

				this->operator()(sigma_r, Im_c);

				assign_real(Im_c, Im);
			}

			void cleanup()
			{
				fft_2d->cleanup();
			}

	protected:
			Vector<T, e_host> gauss_vector_1d(T alpha)
			{
				TVector_r fg;
				fg.reserve(grid_2d.ny);

				for (auto iy=0; iy < grid_2d.ny; iy++)
				{
					auto v = exp(-alpha*grid_2d.gy2_shift(iy))/grid_2d.ny_r();
					fg.push_back(v);
				}
				return fg;
			}

			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			gauss_cv_2d_bc(T sigma_r, TVector_c &M_g)
			{
				auto alpha = 2*c_Pi2*sigma_r*sigma_r;
				auto fg = gauss_vector_1d(alpha);

				stream->set_n_act_stream(grid_2d.nx);
				stream->set_grid(grid_2d.nx, grid_2d.ny);
				stream->exec_matrix(host_device_detail::vector_col_x_matrix<Grid_2d<T>, TVector_r, TVector_c>, grid_2d, fg, M_g);	
			}

			/**********************Device**********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			gauss_cv_2d_bc(T sigma_r, TVector_c &M_g)
			{
				auto alpha = 2*c_Pi2*sigma_r*sigma_r;
				auto fg_h = gauss_vector_1d(alpha);
				TVector_r fg = fg_h;

				auto grid_bt = grid_2d.cuda_grid();
				device_detail::vector_col_x_matrix<Grid_2d<T>, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, fg, M_g);
			}
		#endif

			Stream<dev> *stream;
			FFT<T, dev> *fft_2d;
			Grid_2d<T> grid_2d;
	};

	template <class T, eDevice dev>
	class Gauss_Cv_2d
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			static const eDevice device = dev;

			Gauss_Cv_2d():stream(nullptr), fft_2d(nullptr){}

			Gauss_Cv_2d(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				set_input_data(stream_i, fft_2d_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				stream = stream_i;
				fft_2d = fft_2d_i;
				grid_2d = grid_2d_i;
			}

			void set_fft_plan()
			{
				fft_2d->create_plan_2d(grid_2d.ny, grid_2d.nx, stream->size());
			}

			void operator()(T sigma_r, TVector_c &Im)
			{
				fft2_shift(*stream, grid_2d, Im);
				fft_2d->forward(Im);

				gauss_cv_2d(sigma_r, Im);

				fft_2d->inverse(Im);
				fft2_shift(*stream, grid_2d, Im);
			}

			void operator()(T sigma_r, TVector_r &Im)
			{
				TVector_c Im_c(Im.begin(), Im.end());

				this->operator()(sigma_r, Im_c);

				assign_real(Im_c, Im);
			}

			void cleanup()
			{
				fft_2d->cleanup();
			}
	protected:
			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			gauss_cv_2d(T sigma_r, TVector_c &Im)
			{
				auto alpha = 2*c_Pi2*sigma_r*sigma_r;

				stream->set_n_act_stream(grid_2d.nx);
				stream->set_grid(grid_2d.nx, grid_2d.ny);
				stream->exec_matrix(host_device_detail::gauss_cv_2d<Grid_2d<T>, TVector_c>, grid_2d, alpha, Im);
			}

			/**********************Device**********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			gauss_cv_2d(T sigma_r, TVector_c &Im)
			{
				auto alpha = 2*c_Pi2*sigma_r*sigma_r;

				auto grid_bt = grid_2d.cuda_grid();
				device_detail::gauss_cv_2d<Grid_2d<T>, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, alpha, Im);
			}
		#endif

			Stream<dev> *stream;
			FFT<T, dev> *fft_2d;
			Grid_2d<T> grid_2d;
	};

	/**************** Gaussian Deconv ***************/
	template <class T, eDevice dev>
	class Gauss_Dcv_1d
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			static const eDevice device = dev;

			Gauss_Dcv_1d():fft_1d(nullptr){}

			Gauss_Dcv_1d(FFT<T, dev> *fft_1d_i, Grid_1d<T> &grid_1d_i)
			{
				set_input_data(fft_1d_i, grid_1d_i);
			}

			inline
			void set_input_data(FFT<T, dev> *fft_1d_i, Grid_1d<T> &grid_1d_i)
			{
				fft_1d = fft_1d_i;
				grid_1d = grid_1d_i;
			}

			void set_fft_plan()
			{
				fft_1d->create_plan_1d(grid_1d.nx, 1);
			}

			void operator()(T sigma_r, T PSNR, TVector_c &Im)
			{
				fft1_shift(grid_1d, Im);
				fft_1d->forward(Im);

				gauss_dcv_1d(sigma_r, PSNR, Im);

				fft_1d->inverse(Im);
				fft1_shift(grid_1d, Im);
			}

			void operator()(T sigma_r, T PSNR, TVector_r &Im)
			{
				TVector_c Im_c(Im.begin(), Im.end());

				this->operator()(sigma_r, PSNR, Im_c);

				assign_real(Im_c, Im);
			}

			void cleanup()
			{
				fft_1d->cleanup();
			}
	protected:
			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			gauss_dcv_1d(T sigma_r, T PSNR, TVector_c &Im)
			{
				auto alpha = 2*c_Pi2*sigma_r*sigma_r;

				for (auto ix = 0; ix < grid_1d.nx; ix++)
				{
					host_device_detail::gauss_dcv_1d<Grid_1d<T>, TVector_c>(ix, grid_1d, alpha, PSNR, Im);
				}
			}

			/**********************Device**********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			gauss_dcv_1d(T sigma_r, T PSNR, TVector_c &Im)
			{
				auto alpha = 2*c_Pi2*sigma_r*sigma_r;

				auto grid_bt = grid_1d.cuda_grid();
				device_detail::gauss_dcv_1d<Grid_1d<T>, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_1d, alpha, PSNR, Im);
			}
		#endif

			FFT<T, dev> *fft_1d;
			Grid_1d<T> grid_1d;
	};

	template <class T, eDevice dev>
	class Gauss_Dcv_2d_BC
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			static const eDevice device = dev;

			Gauss_Dcv_2d_BC():stream(nullptr), fft_2d(nullptr){}

			Gauss_Dcv_2d_BC(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				set_input_data(stream_i, fft_2d_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				stream = stream_i;
				fft_2d = fft_2d_i;
				grid_2d = grid_2d_i;
			}

			void set_fft_plan()
			{
				fft_2d->create_plan_1d_batch(grid_2d.ny, grid_2d.nx, stream->size());
			}

			void operator()(T sigma_r, T PSNR, TVector_c &Im)
			{
				fft2_sft_bc(*stream, grid_2d, Im);
				fft_2d->forward(Im);

				gauss_cv_2d_bc(sigma_r, PSNR, Im);

				fft_2d->inverse(Im);
				fft2_sft_bc(*stream, grid_2d, Im);
			}

			void operator()(T sigma_r, T PSNR, TVector_r &Im)
			{
				TVector_c Im_c(Im.begin(), Im.end());

				this->operator()(sigma_r, PSNR, Im_c);

				assign_real(Im_c, Im);
			}

			void cleanup()
			{
				fft_2d->cleanup();
			}

	protected:
			Vector<T, e_host> gauss_vector_1d(T alpha, T PSNR)
			{
				TVector_r fg;
				fg.reserve(grid_2d.ny);

				for (auto iy=0; iy < grid_2d.ny; iy++)
				{
					auto v = exp(-alpha*grid_2d.gy2_shift(iy));
					fg.push_back(v/((v*v+PSNR)*grid_2d.ny_r()));
				}
				return fg;
			}

			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			gauss_cv_2d_bc(T sigma_r, T PSNR, TVector_c &M_g)
			{
				auto alpha = 2*c_Pi2*sigma_r*sigma_r;
				auto fg = gauss_vector_1d(alpha, PSNR);

				stream->set_n_act_stream(grid_2d.nx);
				stream->set_grid(grid_2d.nx, grid_2d.ny);
				stream->exec_matrix(host_device_detail::vector_col_x_matrix<Grid_2d<T>, TVector_r, TVector_c>, grid_2d, fg, M_g);	
			}

			/**********************Device**********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			gauss_cv_2d_bc(T sigma_r, T PSNR, TVector_c &M_g)
			{
				auto alpha = 2*c_Pi2*sigma_r*sigma_r;
				auto fg_h = gauss_vector_1d(alpha, PSNR);
				TVector_r fg = fg_h;

				auto grid_bt = grid_2d.cuda_grid();
				device_detail::vector_col_x_matrix<Grid_2d<T>, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, fg, M_g);
			}
		#endif

			Stream<dev> *stream;
			FFT<T, dev> *fft_2d;
			Grid_2d<T> grid_2d;
	};

	template <class T, eDevice dev>
	class Gauss_Dcv_2d
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			static const eDevice device = dev;

			Gauss_Dcv_2d():stream(nullptr), fft_2d(nullptr){}

			Gauss_Dcv_2d(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				set_input_data(stream_i, fft_2d_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				stream = stream_i;
				fft_2d = fft_2d_i;
				grid_2d = grid_2d_i;
			}

			void set_fft_plan()
			{
				fft_2d->create_plan_2d(grid_2d.ny, grid_2d.nx, stream->size());
			}

			void operator()(T sigma_r, T PSNR, TVector_c &Im)
			{
				fft2_shift(*stream, grid_2d, Im);
				fft_2d->forward(Im);

				gauss_dcv_2d(sigma_r, PSNR, Im);

				fft_2d->inverse(Im);
				fft2_shift(*stream, grid_2d, Im);
			}

			void operator()(T sigma_r, T PSNR, TVector_r &Im)
			{
				TVector_c Im_c(Im.begin(), Im.end());

				this->operator()(sigma_r, PSNR, Im_c);

				assign_real(Im_c, Im);
			}

			void cleanup()
			{
				fft_2d->cleanup();
			}
	protected:
			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			gauss_dcv_2d(T sigma_r, T PSNR, TVector_c &Im)
			{
				auto alpha = 2*c_Pi2*sigma_r*sigma_r;

				stream->set_n_act_stream(grid_2d.nx);
				stream->set_grid(grid_2d.nx, grid_2d.ny);
				stream->exec_matrix(host_device_detail::gauss_dcv_2d<Grid_2d<T>, TVector_c>, grid_2d, alpha, PSNR, Im);
			}

			/**********************Device**********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			gauss_dcv_2d(T sigma_r, T PSNR, TVector_c &Im)
			{
				auto alpha = 2*c_Pi2*sigma_r*sigma_r;

				auto grid_bt = grid_2d.cuda_grid();
				device_detail::gauss_dcv_2d<Grid_2d<T>, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, alpha, PSNR, Im);
			}
		#endif

			Stream<dev> *stream;
			FFT<T, dev> *fft_2d;
			Grid_2d<T> grid_2d;
	};

	/****************Gauss_Spt***************/
	template <class T, eDevice dev>
	class Gauss_Spt
	{
		public:
			using T_r = T;
			using TVector_r = Vector<T, dev>;

			using size_type = std::size_t;

			static const eDevice device = dev;

			Gauss_Spt(): input_gauss_spt(nullptr), stream(nullptr), n_atoms_p(512){}

			void set_input_data(Input_Gauss_Spt<T> *input_gauss_spt_i, Stream<dev> *stream_i)
			{	
				input_gauss_spt = input_gauss_spt_i;
				stream = stream_i;

				n_atoms_p = (device==e_host)?(stream->size()):512;

				if(device==e_host)
				{
					int nv = input_gauss_spt->get_nv();
					stream_data.resize(n_atoms_p);
					for(auto i = 0; i<stream_data.size(); i++)
					{
						stream_data.iv[i].resize(nv);
						stream_data.v[i].resize(nv);
					}
				}

				gauss_sp_h.resize(n_atoms_p);
				gauss_sp.resize(n_atoms_p);
			}
			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			operator()(TVector_r &Im)
			{
				auto gauss_eval = [](Stream<e_host> &stream, Grid_2d<T> &grid_2d, 
				Vector<Gauss_Sp<T>, e_host> &gauss, TVector_r &M_o)
				{
					if(stream.n_act_stream<= 0)
					{
						return;
					}

					for(auto istream = 0; istream < stream.n_act_stream-1; istream++)
					{
						stream[istream] = std::thread(std::bind(host_detail::gauss_eval<T>, std::ref(stream), std::ref(grid_2d), std::ref(gauss[istream]), std::ref(M_o)));
					}

					host_detail::gauss_eval<T>(stream, grid_2d, gauss[stream.n_act_stream-1], M_o);

					stream.synchronize();
				};

				mt::fill(*stream, Im, 0.0);

				int iatom_0 = 0;
				int iatom_e = input_gauss_spt->atoms.size()-1;

				int iatoms = iatom_0;
				while (iatoms <= iatom_e)
				{
					stream->set_n_act_stream(iatom_e-iatoms+1);
					set_gauss_sp(iatoms, stream->n_act_stream, gauss_sp);

					gauss_eval(*stream, input_gauss_spt->grid_2d, gauss_sp, Im);
					iatoms += stream->n_act_stream;
				}

				stream->synchronize();
			}

			/***********************Device***********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			operator()(TVector_r &Im)
			{
				auto get_eval_cubic_poly_gridBT = [](int natoms)->Grid_BT
				{
					Grid_BT grid_bt;
					grid_bt.Blk = dim3(natoms, 1, 1);
					grid_bt.Thr = dim3(c_thrnxny, c_thrnxny, 1);

					return grid_bt;
				};

				mt::fill(*stream, Im, 0.0);

				int iatom_0 = 0;
				int iatom_e = input_gauss_spt->atoms.size()-1;

				int iatoms = iatom_0;
				while (iatoms <= iatom_e)
				{
					int n_atoms = min(n_atoms_p, iatom_e-iatoms+1);
					set_gauss_sp(iatoms, n_atoms, gauss_sp);

					auto grid_bt = get_eval_cubic_poly_gridBT(n_atoms);
					device_detail::gauss_eval<T><<<grid_bt.Blk, grid_bt.Thr>>>(input_gauss_spt->grid_2d, gauss_sp, Im);

					iatoms += n_atoms;
				}
			}
		#endif

			Input_Gauss_Spt<T> *input_gauss_spt;
			Stream<dev> *stream;
		private:
			int n_atoms_p;

			struct Stream_Data
			{
				using value_type = T;
				using size_type = std::size_t;

				static const eDevice device = e_host;

				size_type size() const
				{
					return iv.size();
				}

				void resize(const size_type &new_size)
				{
					iv.resize(new_size);
					v.resize(new_size);
				}

				Vector<Vector<int, dev>, e_host> iv;
				Vector<Vector<T, dev>, e_host> v;
			};

			void set_gauss_sp(int iatoms, int n_atoms_p, Vector<Gauss_Sp<T>, dev> &gauss_sp)
			{
				for(auto istream = 0; istream < n_atoms_p; istream++)
				{
					gauss_sp_h[istream].x = input_gauss_spt->atoms.x[iatoms];
					gauss_sp_h[istream].y = input_gauss_spt->atoms.y[iatoms];
					gauss_sp_h[istream].a = input_gauss_spt->atoms.a[iatoms];
					gauss_sp_h[istream].alpha = input_gauss_spt->alpha(iatoms);
					auto R_max = input_gauss_spt->R_max(iatoms);
					auto R2_max = R_max*R_max;
					gauss_sp_h[istream].R2_tap = pow(0.85*R_max, 2);
					gauss_sp_h[istream].tap_cf = c_i2Pi/(R2_max-gauss_sp_h[istream].R2_tap);
					gauss_sp_h[istream].R2_max = R2_max;
					gauss_sp_h[istream].set_ix0_ixn(input_gauss_spt->grid_2d, R_max);
					gauss_sp_h[istream].set_iy0_iyn(input_gauss_spt->grid_2d, R_max);
					if(device==e_host)
					{
						gauss_sp_h[istream].iv = raw_pointer_cast(stream_data.iv[istream].data());
						gauss_sp_h[istream].v = raw_pointer_cast(stream_data.v[istream].data());
					}
					iatoms++;
				}
				thrust::copy(gauss_sp_h.begin(), gauss_sp_h.end(), gauss_sp.begin());
			}

			Stream_Data stream_data;
			Vector<Gauss_Sp<T>, e_host> gauss_sp_h;
			Vector<Gauss_Sp<T>, dev> gauss_sp;
	};

	/****************affine transformations***************/
	template <class T, eDevice dev>
	class Sc_2d
	{
		public:
			using T_r = T;
			using TVector = Vector<T, dev>;

			static const eDevice device = dev;

			Sc_2d(): stream(nullptr){}

			Sc_2d(Stream<dev> *stream_i, Grid_2d<T> &grid_2d_i)
			{
				set_input_data(stream_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, Grid_2d<T> &grid_2d_i)
			{
				stream = stream_i;
				grid_2d = grid_2d_i;
			}

			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			operator()(TVector &M_i, T sxy, TVector &M_o)
			{
				if(isEqual<T>(sxy, T(1)))
				{
					M_o = M_i;
				}

				int nx_o = get_new_size(grid_2d.nx, sxy);
				int ny_o = get_new_size(grid_2d.ny, sxy);

				Grid_2d<T> grid_2d_o(nx_o, ny_o, nx_o*grid_2d.dRx, ny_o*grid_2d.dRy);

				stream->set_n_act_stream(grid_2d_o.nx);
				stream->set_grid(grid_2d_o.nx, grid_2d_o.ny);
				stream->exec_matrix(host_device_detail::sc_2d<Grid_2d<T>, TVector>, grid_2d, M_i, sxy, grid_2d_o, M_o);
			}

			/**********************Device**********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			operator()(TVector &M_i, T sxy, TVector &M_o)
			{
				if(isEqual<T>(sxy, T(1)))
				{
					M_o = M_i;
				}

				int nx_o = max(int(floor(grid_2d.nx*sxy)), 1);
				int ny_o = max(int(floor(grid_2d.ny*sxy)), 1);

				Grid_2d<T> grid_2d_o(nx_o, ny_o, nx_o*grid_2d.dRx, ny_o*grid_2d.dRy);

				auto grid_bt = grid_2d_o.cuda_grid();
				device_detail::sc_2d<Grid_2d<T>, typename TVector::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, M_i, sxy, grid_2d_o, M_o);
			}
		#endif

		protected:
			Grid_2d<T> grid_2d;
			Stream<dev> *stream;
	};

	template <class T, eDevice dev>
	class Rot_2d
	{
		public:
			using T_r = T;
			using TVector = Vector<T, dev>;

			static const eDevice device = dev;

			Rot_2d(): stream(nullptr){}

			Rot_2d(Stream<dev> *stream_i, Grid_2d<T> &grid_2d_i)
			{
				set_input_data(stream_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, Grid_2d<T> &grid_2d_i)
			{
				stream = stream_i;
				grid_2d = grid_2d_i;
			}

			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			operator()(TVector &M_i, T theta, r2d<T> p0, TVector &M_o)
			{
				if(isZero<T>(theta))
				{
					M_o = M_i;
				}

				// calculate background
				T bg = mean(*stream, M_i);

				Grid_2d<T> grid_2d_o = grid_2d;

				stream->set_n_act_stream(grid_2d_o.nx);
				stream->set_grid(grid_2d_o.nx, grid_2d_o.ny);
				stream->exec_matrix(host_device_detail::rot_2d<Grid_2d<T>, TVector>, grid_2d, M_i, theta, p0, bg, grid_2d_o, M_o);
			}

			/**********************Device**********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			operator()(TVector &M_i, T theta, r2d<T> p0, TVector &M_o)
			{
				if(isZero<T>(theta))
				{
					M_o = M_i;
				}

				// calculate background
				T bg = mean(*stream, M_i);

				Grid_2d<T> grid_2d_o = grid_2d;

				auto grid_bt = grid_2d_o.cuda_grid();
				device_detail::rot_2d<Grid_2d<T>, typename TVector::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, M_i, theta, p0, bg, grid_2d_o, M_o);
			}
		#endif

		protected:
			Grid_2d<T> grid_2d;
			Stream<dev> *stream;
	};

	template <class T, eDevice dev>
	class Rot_Sca_sft_2d
	{
		public:
			using T_r = T;
			using TVector = Vector<T, dev>;

			static const eDevice device = dev;

			Rot_Sca_sft_2d():stream(nullptr){}

			Rot_Sca_sft_2d(Stream<dev> *stream_i, Grid_2d<T> &grid_2d_i)
			{
				set_input_data(stream_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, Grid_2d<T> &grid_2d_i)
			{
				stream = stream_i;
				grid_2d = grid_2d_i;
			}

			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			operator()(TVector &M_i, T theta, r2d<T> p0, T sx, T sy, r2d<T> ps, TVector &M_o)
			{
				// calculate background
				T bg = mean(*stream, M_i);

				int nx_o = get_new_size(grid_2d.nx, sx);
				int ny_o = get_new_size(grid_2d.ny, sy);

				Grid_2d<T> grid_2d_o(nx_o, ny_o, nx_o*grid_2d.dRx, ny_o*grid_2d.dRy);

				stream->set_n_act_stream(grid_2d_o.nx);
				stream->set_grid(grid_2d_o.nx, grid_2d_o.ny);
				stream->exec_matrix(host_device_detail::rot_sca_sft_2d<Grid_2d<T>, TVector>, grid_2d, M_i, theta, p0, sx, sy, ps, bg, grid_2d_o, M_o);
			}

			/**********************Device**********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			operator()(TVector &M_i, T theta, r2d<T> p0, T sx, T sy, r2d<T> ps, TVector &M_o)
			{
				// calculate background
				T bg = mean(*stream, M_i);

				int nx_o = max(int(floor(grid_2d.nx*sx)), 1);
				int ny_o = max(int(floor(grid_2d.ny*sy)), 1);

				Grid_2d<T> grid_2d_o(nx_o, ny_o, nx_o*grid_2d.dRx, ny_o*grid_2d.dRy);

				auto grid_bt = grid_2d_o.cuda_grid();
				device_detail::rot_sca_sft_2d<Grid_2d<T>, typename TVector::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, M_i, theta, p0, sx, sy, ps, bg, grid_2d_o, M_o);
			}
		#endif

		protected:
			Grid_2d<T> grid_2d;
			Stream<dev> *stream;
	};

	template <class T, eDevice dev>
	class Data_Aug_2d
	{
		public:
			using T_r = T;
			using TVector = Vector<T, dev>;

			static const eDevice device = dev;

			Data_Aug_2d():stream(nullptr){}

			Data_Aug_2d(Stream<dev> *stream_i, Grid_2d<T> &grid_2d_i)
			{
				set_input_data(stream_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, Grid_2d<T> &grid_2d_i)
			{
				stream = stream_i;
				grid_2d = grid_2d_i;
			}

			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			operator()(TVector &M_i, T theta, r2d<T> p0, T sx, T sy, r2d<T> ps, T sim, 
			int nx_o, int ny_o, TVector &M_o)
			{
				// calculate maximum
				auto M_max = *thrust::max_element(M_i.begin(), M_i.end());

				// calculate background
				auto g_max = grid_2d.gx_last();
				auto g_min = 0.9*g_max;
				//T bg = sum_over_Det(*stream, grid_2d, g_min, g_max, M_i);

				T bg = 0.6*mean(*stream, M_i);

				T Rx_0 = p0.x*sx-(nx_o/2)*grid_2d.dRx;
				T Ry_0 = p0.y*sy-(ny_o/2)*grid_2d.dRy;

				Grid_2d<T> grid_2d_o(nx_o, ny_o, nx_o*grid_2d.dRx, ny_o*grid_2d.dRy);
				grid_2d_o.set_R_0(Rx_0, Ry_0);

				stream->set_n_act_stream(grid_2d_o.nx);
				stream->set_grid(grid_2d_o.nx, grid_2d_o.ny);
				stream->exec_matrix(host_device_detail::rot_sca_sft_2d<Grid_2d<T>, TVector>, grid_2d, M_i, theta, p0, sx, sy, ps, bg, grid_2d_o, M_o);

				normalized_data(M_o, M_max);

				// add Poisson noise
				M_o = add_poiss_nois(*stream, M_o, sim);

				M_max = *thrust::max_element(M_o.begin(), M_o.end());
				normalized_data(M_o, M_max);
			}

		protected:
			void normalized_data(TVector &M, T M_max)
			{
				std::for_each(M.begin(), M.end(), [M_max](T &v){ v = (v>0)?v/M_max:0;});
			}

			Grid_2d<T> grid_2d;
			Stream<dev> *stream;
	};

	template <class T, eDevice dev>
	class Shx_Scy
	{
		public:
			using T_r = T;
			using TVector = Vector<T, dev>;

			static const eDevice device = dev;

			Shx_Scy():stream(nullptr){}

			Shx_Scy(Stream<dev> *stream_i, Grid_2d<T> &grid_2d_i)
			{
				set_input_data(stream_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, Grid_2d<T> &grid_2d_i)
			{
				stream = stream_i;
				grid_2d = grid_2d_i;
			}

			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			operator()(TVector &M_i, r2d<T> af, TVector &M_o)
			{
				if(isZero<T>(af.x) && isEqual<T>(af.y, T(1)))
				{
					M_o = M_i;
				}

				TVector M_o_t;
				TVector *pM_o = &M_o;

				if (M_i.data() == M_o.data())
				{
					M_o_t.resize(M_i.size());
					pM_o = &M_o_t;
				}

				// calculate background
				T bg = mean(*stream, M_i);

				Grid_2d<T> grid_2d_o = grid_2d;

				stream->set_n_act_stream(grid_2d_o.nx);
				stream->set_grid(grid_2d_o.nx, grid_2d_o.ny);
				stream->exec_matrix(host_device_detail::shx_scy<Grid_2d<T>, TVector>, grid_2d, M_i, af.x, af.y, bg, grid_2d_o, *pM_o);
				
				if (M_i.data() == M_o.data())
				{
					M_o.assign(pM_o->begin(), pM_o->end());
				}
			}

			/**********************Device**********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			operator()(TVector &M_i, r2d<T> af, TVector &M_o)
			{
				if(isZero<T>(af.x) && isEqual<T>(af.y, T(1)))
				{
					M_o = M_i;
				}

				TVector M_o_t;
				TVector *pM_o = &M_o;

				if (M_i.data() == M_o.data())
				{
					M_o_t.resize(M_i.size());
					pM_o = &M_o_t;
				}

				// calculate background
				T bg = mean(*stream, M_i);

				Grid_2d<T> grid_2d_o = grid_2d;

				auto grid_bt = grid_2d.cuda_grid();
				device_detail::shx_scy<Grid_2d<T>, typename TVector::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, M_i, af.x, af.y, bg, grid_2d_o, *pM_o);
				
				if (M_i.data() == M_o.data())
				{
					M_o.assign(pM_o->begin(), pM_o->end());
				}
			}
		#endif

		protected:
			Grid_2d<T> grid_2d;
			Stream<dev> *stream;
	};

	/**********************shift*************************/
	template <class T, eDevice dev>
	class Sft_1d
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			static const eDevice device = dev;

			Sft_1d():fft_1d(nullptr){}

			Sft_1d(FFT<T, dev> *fft_1d_i, Grid_1d<T> &grid_1d_i)
			{
				set_input_data(fft_1d_i, grid_1d_i);
			}

			inline
			void set_input_data(FFT<T, dev> *fft_1d_i, Grid_1d<T> &grid_1d_i)
			{
				fft_1d = fft_1d_i;
				grid_1d = grid_1d_i;
			}

			void set_fft_plan()
			{
				fft_1d->create_plan_1d(grid_1d.nx);
			}

			void operator()(T xs, TVector_c &Im)
			{
				fft1_shift(grid_1d, Im);
				fft_1d->forward(Im);
				exp_g_factor_1d(grid_1d, -c_2Pi*xs, Im, Im);
				fft_1d->inverse(Im);
				fft1_shift(grid_1d, Im);
			}

			void operator()(T xs, TVector_r &Im)
			{
				TVector_c Im_c(Im.begin(), Im.end());
				T Im_min = min_element(Im);

				this->operator()(xs, Im_c);
				assign_real(Im_c, Im, Im_min);
			}

			void cleanup()
			{
				fft_1d->cleanup();
			}

		protected:
			FFT<T, dev> *fft_1d;
			Grid_1d<T> grid_1d;
	};

	template <class T, eDevice dev>
	class Sft_2d_BC
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			static const eDevice device = dev;

			Sft_2d_BC():stream(nullptr), fft_2d(nullptr){}

			Sft_2d_BC(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				set_input_data(stream_i, fft_2d_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				stream = stream_i;
				fft_2d = fft_2d_i;
				grid_2d = grid_2d_i;
			}

			void set_fft_plan()
			{
				fft_2d->create_plan_1d_batch(grid_2d.ny, grid_2d.nx, stream->size());
			}

			void operator()(T alpha, TVector_r ys, TVector_c &Im)
			{
				fft2_sft_bc(*stream, grid_2d, Im);
				fft_2d->forward(Im);
				exp_g_factor_2d_bc(*stream, grid_2d, -c_2Pi*alpha, ys, Im, Im);
				fft_2d->inverse(Im);
				fft2_sft_bc(*stream, grid_2d, Im);
			}

			void operator()(T alpha, TVector_r ys, TVector_r &Im)
			{
				TVector_c Im_c(Im.begin(), Im.end());
				T Im_min = min_element(Im);

				this->operator()(alpha, ys, Im_c);
				assign_real(Im_c, Im, Im_min);
			}

			void cleanup()
			{
				fft_2d->cleanup();
			}

		protected:
			Stream<dev> *stream;
			FFT<T, dev> *fft_2d;
			Grid_2d<T> grid_2d;
	};

	template <class T, eDevice dev>
	class Sft_2d
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			static const eDevice device = dev;

			Sft_2d():stream(nullptr), fft_2d(nullptr){}

			Sft_2d(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				set_input_data(stream_i, fft_2d_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				stream = stream_i;
				fft_2d = fft_2d_i;
				grid_2d = grid_2d_i;
			}

			void set_fft_plan()
			{
				fft_2d->create_plan_2d(grid_2d.ny, grid_2d.nx, stream->size());
			}

			void operator()(r2d<T> p, TVector_c &Im)
			{
				fft2_shift(*stream, grid_2d, Im);
				fft_2d->forward(Im);
				exp_g_factor_2d(*stream, grid_2d, -c_2Pi*p.x, -c_2Pi*p.y, Im, Im);
				fft_2d->inverse(Im);
				fft2_shift(*stream, grid_2d, Im);
			}

			void operator()(r2d<T> p, TVector_r &Im)
			{
				TVector_c Im_c(Im.begin(), Im.end());
				T Im_min = min_element(Im);

				this->operator()(p, Im_c);
				assign_real(Im_c, Im, Im_min);
			}

			void cleanup()
			{
				fft_2d->cleanup();
			}
		protected:
			Stream<dev> *stream;
			FFT<T, dev> *fft_2d;
			Grid_2d<T> grid_2d;
	};

	/******************phase correlation*****************/
	template <class T, eDevice dev>
	class Pcf_1d
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			static const eDevice device = dev;

			Pcf_1d(): fft_1d(nullptr){}

			Pcf_1d(FFT<T, dev> *fft_1d_i, Grid_1d<T> &grid_1d_i)
			{
				set_input_data(fft_1d_i, grid_1d_i);
			}

			inline
			void set_input_data(FFT<T, dev> *fft_1d_i, Grid_1d<T> &grid_1d_i)
			{
				fft_1d = fft_1d_i;
				grid_1d = grid_1d_i;

				Prime_Num pn;

				int nx = pn(2*grid_1d.nx-1, eDST_Greater_Than);

				grid_1d_e.set_input_data(nx, grid_1d.dRx*nx);

				M_r_c.resize(grid_1d_e.nx);
				M_s_c.resize(grid_1d_e.nx);
			}

			void set_fft_plan()
			{
				fft_1d->create_plan_1d(grid_1d.nx, 1);
				fft_1d_e.create_plan_1d(grid_1d_e.nx, 1);
			}

			void operator()(TVector_r &M_r, TVector_r &M_s, T p, T sigma_g, 
			Border_1d<T> bd, bool b_pv, TVector_r &M_o)
			{
				thrust::fill(M_r_c.begin(), M_r_c.end(), T(0));
				thrust::fill(M_s_c.begin(), M_s_c.end(), T(0));

				preprocessing(M_r, p, bd, M_r_c);
				preprocessing(M_s, p, bd, M_s_c);

				TVector_c &pcf = M_s_c;

				// shift matrix
				fft1_shift(grid_1d_e, M_r_c);
				fft1_shift(grid_1d_e, M_s_c);

				// fft_1d_e
				fft_1d_e.forward(M_r_c);
				fft_1d_e.forward(M_s_c);

				pcf_g(M_r_c, M_s_c, sigma_g, pcf);

				fft_1d_e.inverse(pcf);

				// shift pcf
				fft1_shift(grid_1d_e, pcf);

				int ix_s = (grid_1d_e.nx-grid_1d.nx)/2;
				this->assign_real(pcf, ix_s, M_o, b_pv);
			}

			void cleanup()
			{
				fft_1d->destroy_plan();
				fft_1d_e.cleanup();
			}

		protected:

			void assign_real(TVector_c &M_i, int ix_s, TVector_r &M_o, bool b_pos = false)
			{
				auto first = M_i.begin() + ix_s;
				auto last = first + grid_1d.nx;

				if(b_pos)
				{
					thrust::transform(first, last, M_o.begin(), functor::assign_max_real<T>(0));
				}
				else
				{
					thrust::transform(first, last, M_o.begin(), functor::assign_real<T>());
				}
			}

			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			preprocessing(TVector_r &M_i, T p, Border_1d<T> &bd, TVector_c &M_o)
			{
				int ix_s = (grid_1d_e.nx-grid_1d.nx)/2;

				Butterworth_1d<T> bw_1d(bd, bd.radius_ptl(p), 32);

				for (auto ix = 0; ix < grid_1d.nx; ix++)
				{
					host_device_detail::pcf_1d_pp(ix, ix_s, grid_1d, bw_1d, M_i, M_o);
				}
			}

			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			pcf_g(TVector_c &M_r_c, TVector_c &M_s_c, T sigma_g, TVector_c &pcf)
			{
				Gauss_1d<T> gs_1d(0, sigma_g);

				for (auto ix = 0; ix < grid_1d_e.nx; ix++)
				{
					host_device_detail::pcf_1d_gaussian(ix, grid_1d_e, gs_1d, M_r_c, M_s_c, pcf);
				}
			}

		/***********************Device***********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			preprocessing(TVector_r &M_i, T p, Border_1d<T> &bd, TVector_c &M_o)
			{
				int ix_s = (grid_1d_e.nx-grid_1d.nx)/2;

				Butterworth_1d<T> bw_1d(bd, bd.radius_ptl(p), 32);

				auto grid_bt = grid_1d.cuda_grid();
				device_detail::pcf_1d_pp<Grid_1d<T>, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(ix_s, grid_1d, bw_1d, M_i, M_o);
			}

			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			pcf_g(TVector_c &M_r_c, TVector_c &M_s_c, T sigma_g, TVector_c &pcf)
			{
				Gauss_1d<T> gs_1d(0, sigma_g);

				auto grid_bt = grid_1d_e.cuda_grid();
				device_detail::pcf_1d_gaussian<Grid_1d<T>, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_1d_e, gs_1d, M_r_c, M_s_c, pcf);
			}
		#endif

			TVector_c M_r_c;
			TVector_c M_s_c;

			FFT<T, dev> *fft_1d;
			Grid_1d<T> grid_1d;

		private:
			FFT<T, dev> fft_1d_e;
			Grid_1d<T> grid_1d_e;
	};

	template <class T, eDevice dev>
	class Pcf_2d_BC
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			static const eDevice device = dev;

			Pcf_2d_BC():stream(nullptr), fft_2d(nullptr){}

			Pcf_2d_BC(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				set_input_data(stream_i, fft_2d_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				stream = stream_i;
				fft_2d = fft_2d_i;
				grid_2d = grid_2d_i;
				grid_1d.set_input_data(grid_2d.ny, grid_2d.ly);

				Prime_Num pn;

				int ny = pn(2*grid_2d.ny-1, eDST_Greater_Than);

				grid_2d_e.set_input_data(grid_2d.nx, ny, grid_2d.lx, grid_2d.dRy*ny);
				grid_1d_e.set_input_data(grid_2d_e.ny, grid_2d_e.ly);

				M_r_c.resize(grid_2d_e.nxy());
				M_s_c.resize(grid_2d_e.nxy());
			}

			void set_fft_plan()
			{
				fft_2d->create_plan_1d_batch(grid_2d.ny, grid_2d.nx, stream->size());
				fft_2d_e.create_plan_1d_batch(grid_2d_e.ny, grid_2d_e.nx, stream->size());
			}

			void operator()(TVector_r &M_r, TVector_r &M_s, T p, T sigma_g, 
			Border_1d<T> bd, bool b_pv, TVector_r &M_o)
			{
				thrust::fill(M_r_c.begin(), M_r_c.end(), T(0));
				thrust::fill(M_s_c.begin(), M_s_c.end(), T(0));

				preprocessing(M_r, p, bd, M_r_c);
				preprocessing(M_s, p, bd, M_s_c);

				TVector_c &pcf = M_s_c;

				// shift matrix
				fft2_sft_bc(*stream, grid_2d_e, M_r_c);
				fft2_sft_bc(*stream, grid_2d_e, M_s_c);

				// fft_2d
				fft_2d_e.forward(M_r_c);
				fft_2d_e.forward(M_s_c);

				pcf_g(M_r_c, M_s_c, sigma_g, pcf);
				fft_2d_e.inverse(pcf);

				// shift pcf
				fft2_sft_bc(*stream, grid_2d_e, pcf);

				int iy_s = (grid_2d_e.ny-grid_2d.ny)/2;
				this->assign_real(pcf, iy_s, M_o, b_pv);
			}

			void cleanup()
			{
				fft_2d->destroy_plan();
				fft_2d_e.cleanup();
			}

		protected:

			TVector_r gauss_vector_1d(Grid_1d<T> &grid_1d, T sigma_g, bool b_norm = false)
			{
				Gauss_1d<T> gs_1d(0, sigma_g);

				TVector_r fg;
				fg.reserve(grid_1d.nx);

				for (auto ix=0; ix < grid_1d.nx; ix++)
				{
					T g2 = grid_1d.g2_shift(ix);
					T v = gs_1d(g2);
					fg.push_back(v);
				}
				return fg;
			}

			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			preprocessing(TVector_r &M_i, T p, Border_1d<T> &bd, TVector_c &M_o)
			{
				int iy_s = (grid_2d_e.ny-grid_2d.ny)/2;

				auto fh = func_butterworth_1d<TVector_r>(grid_1d, bd.radius_ptl(p), 32, false, bd);

				stream->set_n_act_stream(grid_2d.nx);
				stream->set_grid(grid_2d.nx, grid_2d.ny);
				stream->exec_matrix(host_device_detail::pcf_2d_bc_pp<Grid_2d<T>, TVector_r, TVector_c>, iy_s, grid_2d, grid_2d_e, M_i, fh, M_o);	
			}

			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			pcf_g(TVector_c &M_r_c, TVector_c &M_s_c, T sigma_g, TVector_c &pcf)
			{
				auto fg = gauss_vector_1d(grid_1d_e, sigma_g);

				stream->set_n_act_stream(grid_2d_e.nx);
				stream->set_grid(grid_2d_e.nx, grid_2d_e.ny);
				stream->exec_matrix(host_device_detail::pcf_2d_bc_gaussian<Grid_2d<T>, TVector_r, TVector_c>, grid_2d_e, M_r_c, M_s_c, fg, pcf);	
			}

			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			assign_real(TVector_c &M_i, int iy_s, TVector_r &M_o, bool b_pos = false)
			{
				for(auto ix = 0; ix < grid_2d.nx; ix++)
				{
					for(auto iy = 0; iy < grid_2d.ny; iy++)
					{
						int ixy_i = grid_2d_e.ind_col(ix, iy+iy_s);
						int ixy_o = grid_2d.ind_col(ix, iy);
						auto v = M_i[ixy_i].real();
						M_o[ixy_o] = (!b_pos || (v>0))?v:0;
					}
				}
			}

			/**********************Device**********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			preprocessing(TVector_r &M_i, T p, Border_1d<T> &bd, TVector_c &M_o)
			{
				int iy_s = (grid_2d_e.ny-grid_2d.ny)/2;

				auto fh_h = func_butterworth_1d<Vector<T, e_host>>(grid_1d, bd.radius_ptl(p), 32, false, bd);
				TVector_r fh = fh_h;

				auto grid_bt = grid_2d.cuda_grid();
				device_detail::pcf_2d_bc_pp<Grid_2d<T>, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(iy_s, grid_2d, grid_2d_e, M_i, fh, M_o);	
			}

			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			pcf_g(TVector_c &M_r_c, TVector_c &M_s_c, T sigma_g, TVector_c &pcf)
			{
				auto fg_h = gauss_vector_1d(grid_1d_e, sigma_g);
				TVector_r fg = fg_h;

				auto grid_bt = grid_2d_e.cuda_grid();
				device_detail::pcf_2d_bc_gaussian<Grid_2d<T>, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d_e, M_r_c, M_s_c, fg, pcf);	
			}

			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			assign_real(TVector_c &M_i, int iy_s, TVector_r &M_o, bool b_pos = false)
			{
				auto grid_bt = grid_2d.cuda_grid();
				device_detail::pcf_2d_bc_assign_real<Grid_2d<T>, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(iy_s, grid_2d, grid_2d_e, M_i, M_o, b_pos);	
			}
		#endif

			Vector<T_c, dev> M_r_c;
			Vector<T_c, dev> M_s_c;

			Stream<dev> *stream;
			FFT<T, dev> *fft_2d;
			Grid_2d<T> grid_2d;
			Grid_1d<T> grid_1d;

		private:
			FFT<T, dev> fft_2d_e;
			Grid_2d<T> grid_2d_e;
			Grid_1d<T> grid_1d_e;
	};

	template <class T, eDevice dev>
	class Pcf_2d
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			static const eDevice device = dev;

			Pcf_2d():stream(nullptr), fft_2d(nullptr){}

			Pcf_2d(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				set_input_data(stream_i, fft_2d_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				stream = stream_i;
				fft_2d = fft_2d_i;
				grid_2d = grid_2d_i;

				M_r_c.resize(grid_2d.nxy());
				M_s_c.resize(grid_2d.nxy());
			}

			void set_fft_plan()
			{
				fft_2d->create_plan_2d(grid_2d.ny, grid_2d.nx, stream->size());
			}

			void operator()(TVector_r &M_r, TVector_r &M_s, T p, T sigma_g, 
			Border_2d<T> bd, bool b_pv, TVector_r &M_o)
			{
				preprocessing(M_r, p, bd, M_r_c);
				preprocessing(M_s, p, bd, M_s_c);

				TVector_c &pcf = M_s_c;

				// shift matrix
				fft2_shift(*stream, grid_2d, M_r_c);
				fft2_shift(*stream, grid_2d, M_s_c);

				// fft_2d
				fft_2d->forward(M_r_c);
				fft_2d->forward(M_s_c);

				pcf_g(M_r_c, M_s_c, sigma_g, pcf);
				fft_2d->inverse(pcf);

				// shift pcf
				fft2_shift(*stream, grid_2d, pcf);

				T pv = (b_pv)?1:-1;
				assign_real(pcf, M_o, pv);
			}

			void cleanup()
			{
				fft_2d->cleanup();
			}

		protected:
			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			preprocessing(TVector_r &M_i, T p, Border_2d<T> &bd, TVector_c &M_o)
			{
				Butterworth_2d<T> bw_2d(bd, bd.radius_ptl(p), 32);

				stream->set_n_act_stream(grid_2d.nx);
				stream->set_grid(grid_2d.nx, grid_2d.ny);
				stream->exec_matrix(host_device_detail::pcf_2d_pp<Grid_2d<T>, TVector_r, TVector_c>, grid_2d, bw_2d, M_i, M_o);
			}

			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			pcf_g(TVector_c &M_r_c, TVector_c &M_s_c, T sigma_g, TVector_c &pcf)
			{
				Gauss_2d<T> gs_2d(0, 0, sigma_g);

				stream->set_n_act_stream(grid_2d.nx);
				stream->set_grid(grid_2d.nx, grid_2d.ny);
				stream->exec_matrix(host_device_detail::pcf_2d_gaussian<Grid_2d<T>, TVector_c>, grid_2d, gs_2d, M_r_c, M_s_c, pcf);
			}

			/**********************Device**********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			preprocessing(TVector_r &M_i, T p, Border_2d<T> &bd, TVector_c &M_o)
			{
				Butterworth_2d<T> bw_2d(bd, bd.radius_ptl(p), 32);

				auto grid_bt = grid_2d.cuda_grid();
				device_detail::pcf_2d_pp<Grid_2d<T>, typename TVector_c::value_type><<<grid_bt.Blk, grid_bt.Thr>>>(grid_2d, bw_2d, M_i, M_o);
			}

			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			pcf_g(TVector_c &M_r_c, TVector_c &M_s_c, T sigma_g, TVector_c &pcf)
			{
				Gauss_2d<T> gs_2d(0, 0, sigma_g);

				auto grid_bt = grid_2d.cuda_grid();
				device_detail::pcf_2d_gaussian<Grid_2d<T>, typename TVector_c::value_type><<< grid_bt.Blk, grid_bt.Thr >>>(grid_2d, gs_2d, M_r_c, M_s_c, pcf);
			}
		#endif

			Vector<T_c, dev> M_r_c;
			Vector<T_c, dev> M_s_c;

			Stream<dev> *stream;
			FFT<T, dev> *fft_2d;
			Grid_2d<T> grid_2d;
	};

	/********************find shift**********************/
	template <class T, eDevice dev>
	class Fd_Sft_1d: public Pcf_1d<T, dev>
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			Fd_Sft_1d(): Pcf_1d<T, dev>(){}

			Fd_Sft_1d(FFT<T, dev> *fft_1d_i, Grid_1d<T> &grid_1d_i): Pcf_1d<T, dev>()
			{
				set_input_data(fft_1d_i, grid_1d_i);
			}

			inline
			void set_input_data(FFT<T, dev> *fft_1d_i, Grid_1d<T> &grid_1d_i)
			{
				Pcf_1d<T, dev>::set_input_data(fft_1d_i, grid_1d_i);
				sft_1d.set_input_data(fft_1d_i, grid_1d_i);
			}

			T operator()(TVector_r &M_r, TVector_r &M_s, T p, T sigma_g, T dx, 
			Border_1d<T> bd, int nit_pcf)
			{
				T sigma_r = 1.0/(c_2Pi*sigma_g);
				T radius = ::fmax(3*this->grid_1d.dRx, 0.9*sigma_r);

				TVector_r M = M_s;
				TVector_r pcf(M.size());

				if(nonZero(dx))
				{
					sft_1d(dx, M);
				}

				for (auto it=0; it<nit_pcf; it++)
				{
					Pcf_1d<T, dev>::operator()(M_r, M, p, sigma_g, bd, true, pcf);

					// get maximun position
					T x_c = host_device_detail::max_pos_1d(this->grid_1d, pcf);

					if(it==nit_pcf-1)
					{
						x_c = fit_max_pos_1d(this->grid_1d, pcf, x_c, sigma_r, radius);
					}

					T dx_t = -(x_c - this->grid_1d.lxh());

					if(it<nit_pcf-1)
					{
						sft_1d(dx_t, M);
					}
					dx += dx_t;

					// calculated new borders
					bd.set_bd(dx);
				}

				return dx;
			}

		protected:
			Sft_1d<T, dev> sft_1d;
	};

	template <class T, eDevice dev>
	class Fd_Sft_2d_BC: public Pcf_2d_BC<T, dev>
	{
		public:
			using TB = std::pair<T, bool>;
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;
			using TVector_rh = Vector<T, e_host>;
			using TVector_ch = Vector<T_c, e_host>;

			Fd_Sft_2d_BC(): Pcf_2d_BC<T, dev>() {}

			Fd_Sft_2d_BC(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i): Pcf_2d_BC<T, dev>()
			{
				set_input_data(stream_i, fft_2d_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				Pcf_2d_BC<T, dev>::set_input_data(stream_i, fft_2d_i, grid_2d_i);
				gauss_cv_2d_bc.set_input_data(stream_i, fft_2d_i, grid_2d_i);
				sft_2d_bc.set_input_data(stream_i, fft_2d_i, grid_2d_i);
			}

			TVector_rh operator()(TVector_r &M_r_i, TVector_r &M_s_i, T p, T sigma_g, TVector_rh dx, 
			Border_2d<T> bd_2d, int nit_pcf)
			{
				T sigma_r = 1.0/(c_2Pi*sigma_g);
				T radius = ::fmax(3*this->grid_1d.dRx, 0.9*sigma_r);
				Border_1d<T> bd_1d(bd_2d.ly, bd_2d.yb_0, bd_2d.yb_e);
				Peaks peaks(this->grid_1d, &bd_1d);

				TVector_r M_r = M_r_i;
				TVector_r M_s = M_s_i;
				TVector_r pcf(M_r.size());

				// Gaussian filter
				gauss_cv_2d_bc(sigma_r, M_r);
				gauss_cv_2d_bc(sigma_r, M_s);

				const int ix_0 = this->grid_2d.ceil_dRx(bd_2d.x_0());
				const int ix_e = this->grid_2d.floor_dRx(bd_2d.x_e());

				bool dx_s = false;
				for(auto ix=0; ix<dx.size(); ix++)
				{
					dx_s = dx_s || nonZero(dx[ix]);
				}
				if(dx_s)
				{
					sft_2d_bc(1, dx, M_s);
				}

				TVector_rh dx_t(this->grid_1d.nx, T(0));

				for (auto it=0; it<nit_pcf; it++)
				{
					Pcf_2d_BC<T, dev>::operator()(M_r, M_s, p, sigma_g, bd_1d, false, pcf);

					// get maximum distance
					T d_max = (it==0)?peaks.get_d_max_0(ix_0, ix_e, pcf):peaks.get_d_max(ix_0, ix_e, dx_t);
					d_max = ::fmax(sigma_r, d_max);

					for(auto ix=ix_0; ix<ix_e; ix++)
					{
						TVector_rh y(pcf.begin()+ix*this->grid_1d.nx, pcf.begin()+(ix+1)*this->grid_1d.nx);

						auto py = peaks.find(ix, M_r, M_s, y, d_max);
						if(it==nit_pcf-1)
						{
							py = fit_max_pos_1d(this->grid_1d, y, py, sigma_r, radius);
						}
						dx_t[ix] = -(py - this->grid_1d.lxh());
						dx[ix] += dx_t[ix];
					}

					// shift by column
					if(it<nit_pcf-1)
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
			Gauss_Cv_2d_BC<T, dev> gauss_cv_2d_bc;
			Sft_2d_BC<T, dev> sft_2d_bc;

		private:
			struct Peaks
			{
				public:
					using TVector_ih = Vector<int, e_host>;

					Peaks(): y_c(0), bd_1d(nullptr){}

					Peaks(Grid_1d<T> &grid_1d_i, Border_1d<T> *bd_1d_i)
					{
						set_input_data(grid_1d_i, bd_1d_i);
					}

					inline
					void set_input_data(Grid_1d<T> &grid_1d_i, Border_1d<T> *bd_1d_i)
					{
						grid_1d = grid_1d_i;
						bd_1d = bd_1d_i;
						fft_1d.create_plan_1d(grid_1d.nx, 1);
						sft_1d.set_input_data(&fft_1d, grid_1d);
						
						ix_pk.resize(grid_1d.nx);
						x_pk.resize(grid_1d.nx);
						y_pk.resize(grid_1d.nx);
					}

					T get_d_max_0(int ix_0, int ix_e, TVector_r &pcf)
					{
						T x_c = grid_1d.lxh();
						TVector_rh dx;
						dx.reserve(grid_1d.nx);

						for(auto ix=ix_0; ix<ix_e; ix++)
						{
							TVector_rh y(pcf.begin()+ix*grid_1d.nx, pcf.begin()+(ix+1)*grid_1d.nx);

							// find peaks
							fd_peaks(y);

							// add point
							dx.push_back(fabs(x_pk[idx_x_clt()]-x_c));
						}
						T dx_min, dx_max;
						get_limits(dx, dx_min, dx_max);

						// remove far away points
						auto it = std::remove_if(dx.begin(), dx.end(), [dx_min, dx_max](T a){return (a<dx_min)||(a>dx_max); });
						dx.resize(std::distance(dx.begin(), it));

						T dx_std = sqrt(variance(dx));
						T d_max = max_element(dx) + 2*dx_std;

						return ::fmax(10*grid_1d.dRx, d_max);
					}

					T get_d_max(int ix_0, int ix_e, TVector_rh &dx)
					{
						TVector_rh dx_h(dx.begin()+ix_0, dx.begin()+ix_e);
						return ::fmax(10*grid_1d.dRx, 3*sqrt(variance(dx_h)));
					}

					T find(int icol, TVector_r &M_r, TVector_r &M_s, TVector_rh &y, T d_max)
					{
						const T x_c = grid_1d.lxh();
						y_c = y[grid_1d.nxh];

						// find peaks
						fd_peaks(y);

						// remove peaks
						remove_peaks(x_c-d_max, x_c+d_max);

						// return if x_pk is empty
						if(empty())
						{
							return x_c;
						}

						// return if there is only one element
						if(size() == 1)
						{
							return x_pk.front();
						}

						const int ix_max = idx_y_max();
						const int ix_clt = idx_x_clt();

						// return if closest position is equal to maximum intensity position
						// and the maximum intensity is greater than the second maximum intensity
						if(fabs(x_pk[ix_max]-x_pk[ix_clt])<1)
						{
							T y_thr = 0.61*y_pk[ix_max];
							for(auto ix=0; ix<x_pk.size(); ix++)
							{
								if((ix!=ix_max) && (y_thr<y_pk[ix]))
								{
									return fd_shift(icol, M_r, M_s);
								}
							}
							return x_pk[ix_max];
						}

						return fd_shift(icol, M_r, M_s);
					}

					void destroy_plan()
					{
						fft_1d.destroy_plan();
					}

				private:
					TVector_ih ix_pk;
					TVector_rh x_pk;
					TVector_rh y_pk;

					Grid_1d<T> grid_1d;
					Border_1d<T> *bd_1d;

					FFT<T, e_host> fft_1d;
					Sft_1d<T, e_host> sft_1d;

					T y_c;

					bool empty() const
					{
						return ix_pk.empty();
					}

					int size() const
					{
						return ix_pk.size();
					}

					int idx_y_max() const
					{
						int idx_max = (std::max_element(y_pk.begin(), y_pk.end())-y_pk.begin());

						return idx_max;
					}

					int idx_x_clt() const
					{
						T x_c = grid_1d.lxh();

						int idx_clt = 0;
						if(x_pk.size()>1)
						{
							idx_clt = std::min_element(x_pk.begin(), x_pk.end(), [x_c](T a, T b){return fabs(a-x_c)<fabs(b-x_c); })-x_pk.begin();
						}

						return idx_clt;
					}

					int idx_x_fht() const
					{
						T x_c = grid_1d.lxh();

						int idx_clt = 0;
						if(x_pk.size()>1)
						{
							idx_clt = std::max_element(x_pk.begin(), x_pk.end(), [x_c](T a, T b){return fabs(a-x_c)<fabs(b-x_c); })-x_pk.begin();
						}

						return idx_clt;
					}

					void fd_peaks(TVector_rh &y)
					{
						ix_pk.resize(y.size());
						x_pk.resize(y.size());
						y_pk.resize(y.size());

						const T y_thr = 0;
						int ic = 0;
						for(auto ix=1; ix<y.size()-1; ix++)
						{
							if(y[ix]>y_thr)
							{
								if((y[ix-1]<y[ix]) && (y[ix+1]<y[ix]))
								{
									ix_pk[ic] = ix;
									x_pk[ic] = grid_1d.Rx(ix);
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
						int ic = 0;
						for(auto ix=0; ix<size(); ix++)
						{
							if((x_min<=x_pk[ix]) && (x_pk[ix]<=x_max))
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
						const int ix_c = grid_1d.nxh;
						for(auto ix=0; ix<ix_pk.size(); ix++)
						{
							if(ix_pk[ix]==ix_c)
							{
								return;
							}
						}

						ix_pk.push_back(ix_c);
						x_pk.push_back(grid_1d.Rx(ix_c));
						y_pk.push_back(y_c);
					}

					void get_limits(TVector_rh &y, T &y_min, T &y_max, T f = 4.0)
					{				
						// get mean and standard deviation
						T y_mean, y_std;
						mean_var(y, y_mean, y_std);
						y_std = sqrt(y_std);

						// set limits
						y_min = y_mean - f*y_std;
						y_max = y_mean + f*y_std;
					};

					T fd_shift(int icol, TVector_r &M_r, TVector_r &M_s)
					{
						const T x_c = grid_1d.lxh();

						if((icol==0)||(icol==grid_1d.nx-1))
						{
							return x_c;
						}

						const int ix_max = idx_y_max();
						const int ix_clt = idx_x_clt();

						// select range
						const T x_min = ::fmin(x_pk[ix_clt], x_pk[ix_max])- grid_1d.dRx;
						const T x_max = ::fmax(x_pk[ix_clt], x_pk[ix_max]) + grid_1d.dRx;

						// remove far away peaks
						remove_peaks(x_min, x_max);

						// add central point
						add_central_point();

						// find maximun diplacement
						const T d_fht = fabs(x_pk[idx_x_fht()]-x_c);

						// set borders
						Border_1d<T> bd_l = *bd_1d;
						bd_l.xb_0 = max(bd_l.xb_0, d_fht);
						bd_l.xb_e = max(bd_l.xb_e, d_fht);

						// get indexes
						const int ix_0 = grid_1d.ceil_dRx(bd_l.x_0());
						const int ix_e = grid_1d.floor_dRx(bd_l.x_e());

						// get reference data
						TVector_rh yr_b(M_r.begin()+(icol-1)*grid_1d.nx, M_r.begin()+icol*grid_1d.nx);
						TVector_rh yr(M_r.begin()+icol*grid_1d.nx, M_r.begin()+(icol+1)*grid_1d.nx);
						TVector_rh yr_n(M_r.begin()+(icol+1)*grid_1d.nx, M_r.begin()+(icol+2)*grid_1d.nx);

						// get shifted data
						TVector_rh ys_0(M_s.begin()+icol*grid_1d.nx, M_s.begin()+(icol+1)*grid_1d.nx);

						TVector_rh chi2_pk(size());

						for(auto ix_pk=0; ix_pk<size(); ix_pk++)
						{
							TVector_rh ys = ys_0;

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
								host_device_detail::kh_sum(chi2, chi2_p, chi2_ee);
							}
							chi2_pk[ix_pk] = chi2;
						}
						int ix_min = std::min_element(chi2_pk.begin(), chi2_pk.end())-chi2_pk.begin();

						return x_pk[ix_min];
					}
			};
	};

	template <class T, eDevice dev>
	class Fd_Sft_2d: public Pcf_2d<T, dev>
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			Fd_Sft_2d(): Pcf_2d<T, dev>(){}

			Fd_Sft_2d(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i): Pcf_2d<T, dev>()
			{
				set_input_data(stream_i, fft_2d_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				Pcf_2d<T, dev>::set_input_data(stream_i, fft_2d_i, grid_2d_i);
				sft_2d.set_input_data(stream_i, fft_2d_i, grid_2d_i);
			}

			r2d<T> operator()(TVector_r &M_r, TVector_r &M_s, T p, T sigma_g, r2d<T> dr, 
			Border_2d<T> bd, int nit_pcf)
			{
				T sigma_r = 1.0/(c_2Pi*sigma_g);
				T radius = ::fmax(3*this->grid_2d.dR_min(), 0.9*sigma_r);

				TVector_r M = M_s;
				TVector_r pcf(M.size());

				if(nonZero(dr))
				{
					sft_2d(dr, M);
				}

				for (auto it=0; it<nit_pcf; it++)
				{
					Pcf_2d<T, dev>::operator()(M_r, M, p, sigma_g, bd, true, pcf);

					// get maximun position
					r2d<T> r_c = host_device_detail::max_pos_2d(this->grid_2d, pcf);

					if(it==nit_pcf-1)
					{
						r_c = fit_max_pos_2d(this->grid_2d, pcf, r_c, sigma_r, radius);
					}

					r2d<T> dr_t = -(r_c - r2d<T>(this->grid_2d.lxh(), this->grid_2d.lyh()));

					if(it<nit_pcf-1)
					{
						sft_2d(dr_t, M);
					}
					dr += dr_t;

					// calculated new borders
					bd.set_bd(dr);
				}

				return dr;
			}

		protected:
			Sft_2d<T, dev> sft_2d;
	};

	/*******************Corrrect shift*******************/
	template <class T, eDevice dev>
	class Crt_Sft_1d: public Fd_Sft_1d<T, dev>
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			Crt_Sft_1d(): Fd_Sft_1d<T, dev>(){}

			Crt_Sft_1d(FFT<T, dev> *fft_1d_i, Grid_1d<T> &grid_1d_i): 
			Fd_Sft_1d<T, dev>(fft_1d_i, grid_1d_i){}

			T operator()(TVector_r &M_r_i, TVector_r &M_s_io, T p, T sigma_g, 
			Border_1d<T> bd, int nit_pcf)
			{
				T dx = 0;
				dx = Fd_Sft_1d<T, dev>::operator()(M_r_i, M_s_io, p, sigma_g, dx, bd, nit_pcf);
				this->sft_1d(dx, M_s_io);

				return dx;
			}
	};

	template <class T, eDevice dev>
	class Crt_Sft_2d_BC: public Fd_Sft_2d_BC<T, dev>
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;
			using TVector_rh = Vector<T, e_host>;
			using TVector_ch = Vector<T_c, e_host>;

			Crt_Sft_2d_BC(): Fd_Sft_2d_BC<T, dev>(){}

			Crt_Sft_2d_BC(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i): 
			Fd_Sft_2d_BC<T, dev>(stream_i, fft_2d_i, grid_2d_i){}

			TVector_rh operator()(TVector_r &M_r_i, TVector_r &M_s_io, T p, T sigma_g, 
			Border_2d<T> bd, int nit_pcf)
			{
				TVector_rh dr(this->grid_2d.ny, T(0));
				dr = Fd_Sft_2d_BC<T, dev>::operator()(M_r_i, M_s_io, p, sigma_g, dr, bd, nit_pcf);
				this->sft_2d_bc(1, dr, M_s_io);

				return dr;
			}
	};

	template <class T, eDevice dev>
	class Crt_Sft_2d: public Fd_Sft_2d<T, dev>
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			Crt_Sft_2d(): Fd_Sft_2d<T, dev>(){}

			Crt_Sft_2d(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i): 
				Fd_Sft_2d<T, dev>(stream_i, fft_2d_i, grid_2d_i){}

			r2d<T> operator()(TVector_r &M_r_i, TVector_r &M_s_io, T p, T sigma_g, 
			Border_2d<T> bd, int nit_pcf)
			{
				r2d<T> dr(0, 0);
				dr = Fd_Sft_2d<T, dev>::operator()(M_r_i, M_s_io, p, sigma_g, dr, bd, nit_pcf);
				this->sft_2d(dr, M_s_io);

				return dr;
			}
	};

	/****************calculate Chi^2*****************/
	template <class T, eDevice dev>
	class Chi2_Pcf_2d: public Crt_Sft_2d<T, dev>
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			Chi2_Pcf_2d(): Crt_Sft_2d<T, dev>(){}

			Chi2_Pcf_2d(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i): Crt_Sft_2d<T, dev>()
			{
				set_input_data(stream_i, fft_2d_i, grid_2d_i);
			}

			inline
			void set_input_data(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i)
			{
				Crt_Sft_2d<T, dev>::set_input_data(stream_i, fft_2d_i, grid_2d_i);
				shx_scy.set_input_data(this->stream, this->grid_2d);
			}

			T operator()(TVector_r &M_r_i, TVector_r &M_s_i, T p, T sigma_g, 
			r2d<T> af, Border_2d<T> bd, int nit_pcf, r2d<T> &ds)
			{
				TVector_r M(M_r_i.size());
				shx_scy(M_s_i, af, M);

				// correct shift and set borders
				ds = Crt_Sft_2d<T, dev>::operator()(M_r_i, M, p, sigma_g, bd, nit_pcf);

				// pcf
				TVector_r &pcf = M;
				Pcf_2d<T, dev>::operator()(M_r_i, M, p, sigma_g, bd, true, pcf);

				// cost function
				return mean(*(this->stream), pcf);
			}

			T operator()(TVector_r &M_r_i, TVector_r &M_s_i, T p, T sigma_g, 
			r2d<T> af, Border_2d<T> bd, int nit_pcf, r2d<T> &ds, Vector<T, e_host> &coef)
			{
				TVector_r M(M_r_i.size());
				shx_scy(M_s_i, af, M);

				// correct shift and set borders
				ds = Crt_Sft_2d<T, dev>::operator()(M_r_i, M, p, sigma_g, bd, nit_pcf);

				// pcf
				TVector_r &pcf = M;
				Pcf_2d<T, dev>::operator()(M_r_i, M, p, sigma_g, bd, true, pcf);

				// cost function
				auto chi2 = mean(*(this->stream), pcf);

				// get maximun position
				r2d<T> r_c = host_device_detail::max_pos_2d(this->grid_2d, pcf);

				// fitting
				T sigma_r = 1.0/(c_2Pi*sigma_g);
				T radius = ::fmax(3*this->grid_2d.dR_min(), 1.5*sigma_r);
				Vector<T, e_host> pcf_h = pcf;
				coef = fit_ellipt_gauss_2d(this->grid_2d, pcf_h, r_c, sigma_r, radius);
				coef[0] = ds.x + this->grid_2d.lxh();
				coef[1] = ds.y + this->grid_2d.lyh();

				return chi2;
			}

		protected:
			Shx_Scy<T, dev> shx_scy;
	};

	template <class T, eDevice dev>
	class NM_Shx_Scy: public Chi2_Pcf_2d<T, dev>
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using TVector_r = Vector<T, dev>;
			using TVector_c = Vector<T_c, dev>;

			NM_Shx_Scy():Chi2_Pcf_2d<T, dev>(), nit_pcf(2){}

			NM_Shx_Scy(Stream<dev> *stream_i, FFT<T, dev> *fft_2d_i, Grid_2d<T> &grid_2d_i): 
			Chi2_Pcf_2d<T, dev>(stream_i, fft_2d_i, grid_2d_i), nit_pcf(2){}

			void operator()(TVector_r &M_r, TVector_r &M_s, T p, T sigma_g, 
			Vector<Afp_2<T>, e_host> &spx, Border_2d<T> bd, int nit_nm)
			{
				int ic = 1;

				T d_pix_min = this->grid_2d.dg_min();
				T d_min = 0.05*d_pix_min;
				T alpha = 1.0;
				T beta = 0.5;
				T gamma = 2.0;

				for(auto it=0; it<nit_nm; it++)
				{
					sort(spx);
					T d_l = max_length(spx);
					if(d_l<d_min)
					{
						break;
					}

					if(d_l<d_pix_min)
					{
						if(ic==1)
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
					r2d<T> ds_r(0, 0);
					auto chi2_r = chi2_pcf(M_r, M_s, p, sigma_g, x_r, bd, nit_pcf, ds_r);

					if(chi2_r<chi2_b)
					{
						auto x_e = x_m + gamma*(x_r-x_m);
						r2d<T> ds_e(0, 0);
						auto chi2_e = chi2_pcf(M_r, M_s, p, sigma_g, x_e, bd, nit_pcf, ds_e); 
						if(chi2_e<chi2_r)
						{
							spx[2] = Afp_2<T>(x_e, ds_e, chi2_e); 
						}
						else
						{
							spx[2] = Afp_2<T>(x_r, ds_r, chi2_r); 
						}
					}
					else if(chi2_r<spx[1].chi2)
					{
						spx[2] = Afp_2<T>(x_r, ds_r, chi2_r);
					}
					else if(chi2_r<chi2_w)
					{
						auto x_rc = x_m + beta*(x_r-x_m);
						r2d<T> ds_rc(0, 0);
						auto chi2_rc = chi2_pcf(M_r, M_s, p, sigma_g, x_rc, bd, nit_pcf, ds_rc); 
						if(chi2_rc<chi2_r)
						{
							spx[2] = Afp_2<T>(x_rc, ds_rc, chi2_rc); 
						}
						else
						{
							contract_simplex(M_r, M_s, p, sigma_g, spx, bd, nit_pcf);
						}
					}
					else
					{
						auto x_wc = x_m + beta*(x_w-x_m);
						r2d<T> ds_wc(0, 0);
						auto chi2_wc = chi2_pcf(M_r, M_s, p, sigma_g, x_wc, bd, nit_pcf, ds_wc); 
						if(chi2_wc<chi2_w)
						{
							spx[2] = Afp_2<T>(x_wc, ds_wc, chi2_wc); 
						}
						else
						{
							contract_simplex(M_r, M_s, p, sigma_g, spx, bd, nit_pcf);
						}
					}
				}

			} 

			void operator()(TVector_r &M_r, TVector_r &M_s, T p, T sigma_g, 
			r2d<T> &af, Border_2d<T> bd, r2d<T> &ds, int nit_nm)
			{
				auto spx = set_simplex(M_r, M_s, p, sigma_g, af, ds, bd, nit_nm);
				this->operator()(M_r, M_s, p, sigma_g, spx, bd, nit_nm);

				af = spx[0].f;
				ds = spx[0].ds;

				// // shear and scaling
				// shx_scy(M_s, af, M_s);
				// // correct shift
				// ds = Fd_Sft_2d<T, dev>::operator()(M_r, M_s, p, sigma_g, ds, bd, nit_pcf);
			}

		protected:
			r2d<T> x_af(const r2d<T> &af, const r2d<T> &r)
			{
				return r2d<T>(r.x+af.x*r.y, af.y*r.y);
			}

			r2d<T> x_iaf(const r2d<T> &af, const r2d<T> &r)
			{
				return r2d<T>(r.x-af.x*r.y/af.y, r.y/af.y);
			}

			Vector<Afp_2<T>, e_host> set_simplex(TVector_r &M_r, TVector_r &M_s, T p, T sigma_g, 
			r2d<T> x_0, r2d<T> ds_0, Border_2d<T> &bd, int nit_nm)
			{
				// global shift
				if(isZero(ds_0))
				{
					ds_0 = Fd_Sft_2d<T, dev>::operator()(M_r, M_s, p, sigma_g, ds_0, bd, nit_pcf);
					bd.set_bd(ds_0);
				}

				TVector_r M_s_t = M_s;
				this->sft_2d(x_iaf(x_0, ds_0), M_s_t);

				// determine coefficients
				Vector<T, e_host> coef(6);
				chi2_pcf(M_r, M_s_t, p, sigma_g, x_0, bd, nit_pcf, ds_0, coef);

				// calculate dd0
				T ff = coef[4]/coef[3];
				ff = (ff<1)?1/ff:ff;
				ff = ::fmax(1.05, ff);

				T dd0 = sqrt(1.354e-04*pow(ff-1, 2)+6.622e-4*(ff-1));
				dd0 = ::fmax(dd0, this->grid_2d.dg_min());

				T sin_t = sin(coef[5]);
				T cos_t = cos(coef[5]);

				T sigma_x = pow(cos_t/coef[3], 2)+pow(sin_t/coef[4], 2);
				sigma_x = sqrt(1/sigma_x);

				T sigma_y = pow(sin_t/coef[3], 2)+pow(cos_t/coef[4], 2);
				sigma_y = sqrt(1/sigma_y);

				T theta = atan2(sigma_y, sigma_x); 

				r2d<T> u(dd0*cos(theta), dd0*sin(theta));
				r2d<T> v(-u.y, u.x);

				// set simplex
				Vector<Afp_2<T>, e_host> spx(3);

				spx[0].f = x_0;
				spx[0].chi2 = chi2_pcf(M_r, M_s, p, sigma_g, spx[0].f, bd, nit_pcf, spx[0].ds);

				spx[1].f = x_0+u;
				spx[1].chi2 = chi2_pcf(M_r, M_s, p, sigma_g, spx[1].f, bd, nit_pcf, spx[1].ds);

				spx[2].f = x_0+v;
				spx[2].chi2 = chi2_pcf(M_r, M_s, p, sigma_g, spx[2].f, bd, nit_pcf, spx[2].ds);

				// sort simplex
				sort(spx);

				return spx;
			}

			T chi2_pcf(TVector_r &M_r, TVector_r &M_s, T p, T sigma_g, 
			r2d<T> af, Border_2d<T> &bd, int nit_pcf, r2d<T> &ds)
			{
				return Chi2_Pcf_2d<T, dev>::operator()(M_r, M_s, p, sigma_g, af, bd, nit_pcf, ds);
			}

			T chi2_pcf(TVector_r &M_r, TVector_r &M_s, T p, T sigma_g, 
			r2d<T> af, Border_2d<T> &bd, int nit_pcf, r2d<T> &ds, Vector<T, e_host> &coef)
			{
				return Chi2_Pcf_2d<T, dev>::operator()(M_r, M_s, p, sigma_g, af, bd, nit_pcf, ds, coef);
			}

			void sort(Vector<Afp_2<T>, e_host> &spx)
			{
				std::sort(spx.begin(), spx.end(), [](const Afp_2<T> &x, const Afp_2<T> &y){ return x.chi2<y.chi2; });
			}

			T max_length(Vector<Afp_2<T>, e_host> &spx)
			{
				auto r0 = spx[0].f;

				r2d<T> dr = spx[1].f - r0;
				T d_max = dr.module();
				for(auto i=2; i<spx.size(); i++)
				{
					r2d<T> dr = spx[i].f - r0;
					d_max = ::fmax(d_max, dr.module());
				}
				return d_max;
			}

			T min_length(Vector<Afp_2<T>, e_host> &spx)
			{
				auto r0 = spx[0].f;
				r2d<T> dr = spx[1].f - r0;
				T d_min = dr.module();
				for(auto i=2; i<spx.size(); i++)
				{
					r2d<T> dr = spx[i].f - r0;
					d_min = ::fmin(d_min, dr.module());
				}
				return d_min;
			}

			r2d<T> best_centroid(Vector<Afp_2<T>, e_host> &spx)
			{
				r2d<T> r_c(0, 0);
				for(auto i=0; i<spx.size()-1; i++)
				{
					r_c += spx[i].f;
				}
				return r_c/T(spx.size()-1);
			}

			void orthogonal_simplex(TVector_r &M_r, TVector_r &M_s, T p, T sigma_g, 
			Vector<Afp_2<T>, e_host> &spx, Border_2d<T> &bd, int nit_pcf)
			{
				auto p12 = spx[1].f-spx[0].f;
				auto p13 = spx[2].f-spx[0].f;
				auto mp12 = p12.module();
				auto mp13 = p13.module();
				auto mp_max = ::fmax(mp12, mp13);

				T theta = angle(p12, p13);

				T theta_min = 10;
				T theta_0 = theta_min*c_Pi/180;
				T theta_e = c_Pi-theta_min*c_Pi/180;
	
				bool b_m = (mp12<mp_max/4)||(mp13<mp_max/4);

				if((theta<theta_0)||(theta>theta_e)||b_m)
				{
					 if(b_m && (mp12<mp13))
					 {
						auto u = p13/module(p13);
						auto x_o = spx[0].f + mp_max*r2d<T>(-u.y, u.x);
						r2d<T> ds_o(0, 0);
						auto chi2_o = chi2_pcf(M_r, M_s, p, sigma_g, x_o, bd, nit_pcf, ds_o);

						spx[1] = Afp_2<T>(x_o, ds_o, chi2_o);
					 }
					 else
					 {
						T m = max_length(spx);
						auto u = p12/module(p12);
						auto x_o = spx[0].f + mp_max*r2d<T>(-u.y, u.x);
						r2d<T> ds_o(0, 0);
						auto chi2_o = chi2_pcf(M_r, M_s, p, sigma_g, x_o, bd, nit_pcf, ds_o);

						spx[2] = Afp_2<T>(x_o, ds_o, chi2_o);
					 }

					sort(spx);
				}
			}

			void contract_simplex(TVector_r &M_r, TVector_r &M_s, T p, T sigma_g, 
			Vector<Afp_2<T>, e_host> &spx, Border_2d<T> &bd, int nit_pcf)
			{
				T ff = 0.5;

				auto p12 = spx[1].f-spx[0].f;
				auto p13 = spx[2].f-spx[0].f;
				auto mp12 = p12.module();
				auto mp13 = p13.module();

				if(mp12<mp13)
				{
					auto u = p12/mp12;
					auto x_c = spx[0].f + ff*mp13*r2d<T>(-u.y, u.x);
					r2d<T> ds_c(0, 0);
					auto chi2_c = chi2_pcf(M_r, M_s, p, sigma_g, x_c, bd, nit_pcf, ds_c);

					spx[2] = Afp_2<T>(x_c, ds_c, chi2_c);
				}
				else
				{
					auto u = p13/mp13;
					auto x_c = spx[0].f + ff*mp12*r2d<T>(-u.y, u.x);
					r2d<T> ds_c(0, 0);
					auto chi2_c = chi2_pcf(M_r, M_s, p, sigma_g, x_c, bd, nit_pcf, ds_c);

					spx[1] = Afp_2<T>(x_c, ds_c, chi2_c);
				}

				sort(spx);
			}

			int nit_pcf;
	};

} // namespace mt

#endif