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

#include "cpu_fcns_image.h"

namespace mt
{
	/* threshold */
	namespace cpu_fcn_image
	{
		/* otsu threshold*/
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type<TVctr>>
		fcn_otsu_thr(const TVctr& vctr, const dt_uint32& n_bins)
		{
			using T = dt_float64;
			using size_type = Size_type<TVctr>;

			Vctr_cpu<T> vctr_rang(n_bins);
			Vctr_cpu<T> vctr_hist(n_bins);

			/* get histogram */
			fcn_hist(vctr, n_bins, vctr_hist, vctr_rang);

			T sum_i = T();
			for(size_type ibins = 0; ibins < n_bins; ibins++)
			{
				sum_i += vctr_rang[ibins]*vctr_hist[ibins];
			}

			T sum_b = T();

			dt_int64 w_b = dt_int64();
			dt_int64 w_f = dt_int64();

			T var_max = T();
			T thr = T();

			dt_int64 vctr_size = dt_int64(vctr.size());

			for(size_type ibins = 0; ibins < n_bins; ibins++)
			{
				/* weight background */
				w_b += dt_int64(vctr_hist[ibins]);

				if (w_b == 0)
				{
					continue;
				}

				/* weight foreground */
				w_f = vctr_size - w_b;
				if (w_f == 0)
				{
					break;
				}

				sum_b += vctr_rang[ibins]*vctr_hist[ibins];

				/* fcn_mean background */
				T m_B = sum_b/T(w_b);

				/* fcn_mean Foreground */
				T m_F = (sum_i - sum_b)/T(w_f);

				/* calculate Between class Variance */
				auto var_between = T(w_b)*T(w_f)*::square(m_B-m_F);

				/* check if new maximum found */
				if (var_between > var_max) 
				{
					var_max = var_between;
					thr = vctr_rang[ibins];
				}
			}

			return Value_type<TVctr>(thr);
		}

		/* matrix binarization */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_binarize_mx(const TVctr& mx_i, const Value_type<TVctr>& thr, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			Vctr_cpu<T> mx_o(mx_i.shape());

			auto thr_binarize_mx = [&](const iThread_Rect_1d &ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, mx_o.begin() + ithread.ind_0, cgpu_fctr::binarize<T>(thr));
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_binarize_mx);

			return mx_o;
		}

		/* maximum threshold */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_threshold_max(const TVctr& mx_i, const Value_type<TVctr>& thr, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			Vctr_cpu<T> mx_o(mx_i.shape());

			auto thr_threshold_max = [&](const iThread_Rect_1d &ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, mx_o.begin() + ithread.ind_0, cgpu_fctr::threshold_max<T>(thr));
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_threshold_max);

			return mx_o;
		}

		/* maximum threshold */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_threshold_min(const TVctr& mx_i, const Value_type<TVctr>& thr, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			Vctr_cpu<T> mx_o(mx_i.shape());

			auto thr_threshold_min = [&](const iThread_Rect_1d &ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, mx_o.begin() + ithread.ind_0, cgpu_fctr::threshold_min<T>(thr));
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_threshold_min);

			return mx_o;
		}
	}

	/* morphological operations */
	namespace cpu_fcn_image
	{
		/* gray scale dilation */
		template <class TVctr, class TFcn>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_morp_op(TVctr& mx_i, dt_int32 nkr, TFcn fcn, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			Vctr_cpu<T> mx_o(mx_i.shape());
			
			/* column dilation and transpose*/
			auto thr_morp_op = [nkr, fcn](const iThread_Rect_2d& ithread, pVctr_cpu_32<T>& mx_i, pVctr_cpu_32<T>& mx_o)
			{
				std::deque<dt_int32> wd;

				const dt_int32 ny = mx_i.m_s0;
				const dt_int32 nx = mx_i.m_s1;
				const dt_int32 wkr = 2*nkr+1;

				for(auto ix = ithread.ix_0; ix < ithread.ix_e; ix++)
				{
					dt_int32 iy = 0;
					dt_int32 iy_tk = 0;
					dt_int32 iy_o = 0;

					for(; iy < nkr; iy++)
					{
						fcn(ix, iy, mx_i, wd);
					}

					for(; iy < ny; iy++, iy_o++)
					{
						fcn(ix, iy, mx_i, wd);

						mx_o[iy_o*nx + ix] = mx_i(wd.front(), ix);
						if (iy + 1 >= wkr)
						{
							if (wd.front() == iy_tk)
							{
								wd.pop_front();
							}
							iy_tk++;
						}
					}

					for(; iy_o < ny; iy_o++)
					{
						mx_o[iy_o*nx + ix] = mx_i(wd.front(), ix);
						if (wd.front() == iy_tk)
						{
							wd.pop_front();
						}
						iy_tk++;
					}
					wd.clear();
				}
			};
			
			/* run dilation along y */
			fcn_stream_exec_xd_fcn<edim_2>(pstream, mx_i.s1_32(), mx_i.s0_32(), thr_morp_op, mx_i.ptr_32(), mx_o.ptr_32());

			auto mx_t = mx_o;

			/* transpose shape */
			mx_t.trs_shape_2d(); 

			/* run dilation along y */
			fcn_stream_exec_xd_fcn<edim_2>(pstream, mx_t.s1_32(), mx_t.s0_32(), thr_morp_op, mx_t.ptr_32(), mx_o.ptr_32());

			return mx_o;
		}

		/* gray value data dilation */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_morp_op_dilate(TVctr& mx_i, dt_int32 nkr, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto fcn_wedge_push = [](const dt_int32& ix, const dt_int32& iy, pVctr_cpu_32<T>& mx_i, std::deque<dt_int32>& wd)
			{
				while ((!wd.empty()) && (mx_i(wd.back(), ix) <= mx_i(iy, ix)))
				{
					wd.pop_back();
				}
				wd.push_back(iy);
			};

			return fcn_morp_op(mx_i, nkr, fcn_wedge_push, pstream);
		}
		
		/* gray value data erosion */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_morp_op_erode(TVctr& mx_i, dt_int32 nkr, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto fcn_wedge_push = [](const dt_int32& ix, const dt_int32& iy, pVctr_cpu_32<T>& mx_i, std::deque<dt_int32>& wd)
			{
				while ((!wd.empty()) && (mx_i(wd.back(), ix) >= mx_i(iy, ix)))
				{
					wd.pop_back();
				}
				wd.push_back(iy);
			};

			return fcn_morp_op(mx_i, nkr, fcn_wedge_push, pstream);
		}

		/* gray value data opening */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_morp_op_open(TVctr& mx_i, dt_int32 nkr, Stream_cpu* pstream = nullptr)
		{
			auto mx_o = fcn_morp_op_erode(mx_i, nkr, pstream);
			return fcn_morp_op_dilate(mx_o, nkr, pstream);
		}

		/* gray value data closing */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_morp_op_close(TVctr& mx_i, dt_int32 nkr, Stream_cpu* pstream = nullptr)
		{
			auto mx_o = fcn_morp_op_dilate(mx_i, nkr, pstream);
			return fcn_morp_op_erode(mx_o, nkr, pstream);
		}

		/* gray value datatophat */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_morp_op_tophat(TVctr& mx_i, dt_int32 nkr, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			auto mx_o = fcn_morp_op_open(mx_i, nkr, pstream);
			fcn_sub(mx_i.ptr_32(), mx_o.ptr_32(), mx_o.ptr_32(), pstream);

			return mx_o;
		}
	}

	/* filter wiener */
	namespace cpu_fcn_image
	{
		/* 1d */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_wiener_1d(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			if (n_kr<=0)
			{
				memcpy_cpu_cpu(mx_o.m_data, mx_i.m_data, mx_i.size_64());
				return;
			}

			const dt_int32 n_k0 = -n_kr;
			const dt_int32 n_ke = n_kr+1;
			const dt_int32 nx_i = mx_i.size_32();

			Vctr_cpu<T> mx_mean(nx_i, T(0));
			Vctr_cpu<T> mx_var(nx_i, T(0));

			/* get local mean and var */
			auto krn_mean_var = [n_k0, n_ke, nx_i, &mx_i](const dt_int32& ix_i, Vctr_cpu<T>& mx_mean, Vctr_cpu<T>& mx_var, KS<T>& v2)
			{
				const dt_int32 ix_0 = fcn_max(ix_i+n_k0, 0);
				const dt_int32 ix_e = fcn_min(ix_i+n_ke, nx_i);

				T v_mean = 0;
				T v_var = 0;
				for(auto ix = ix_0; ix < ix_e; ix++)
				{
					const auto x = mx_i[ix];
					v_mean += x;
					v_var += x*x;
				}

				v_mean = v_mean/T(ix_e-ix_0);
				v_var = v_var/T(ix_e-ix_0) - v_mean*v_mean;

				mx_mean[ix_i] = v_mean;
				mx_var[ix_i] = v_var;

				v2 += v_var;
			};

			auto v2 = fcn_stream_exec_xd_krn_reduce<edim_1>(pstream, nx_i, std::plus<KS<T>>(), T(0), krn_mean_var, mx_mean, mx_var);
			v2 /= T(nx_i);

			/* apply wiener filter */
			auto krn_fltr_wiener = [v2, &mx_i, &mx_o](const dt_int32& ix, Vctr_cpu<T>& mx_mean, Vctr_cpu<T>& mx_var)
			{
				mx_o[ix] = mx_mean[ix] + ::fmax(T(0), mx_var[ix]-v2)*(mx_i[ix]-mx_mean[ix])/::fmax(mx_var[ix], v2);
			};

			fcn_stream_exec_xd_krn<edim_1>(pstream, nx_i, krn_fltr_wiener, mx_mean, mx_var);
		}
		
		/* 2d by col */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_wiener_2d_bc(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			if (n_kr<=0)
			{
				memcpy_cpu_cpu(mx_o.m_data, mx_i.m_data, mx_i.size_64());
				return;
			}

			const dt_int32 n_k0 = -n_kr;
			const dt_int32 n_ke = n_kr+1;
			const dt_int32 nx_i = mx_i.s1_32();
			const dt_int32 ny_i = mx_i.s0_32();

			Vctr_cpu<T> mx_mean(mx_i.shape(), T(0));
			Vctr_cpu<T> mx_var(mx_i.shape(), T(0));

			/* get local mean and var */
			auto krn_mean_var = [n_k0, n_ke, ny_i, &mx_i](const dt_int32& ix_i, const dt_int32& iy_i, Vctr_cpu<T>& mx_mean, Vctr_cpu<T>& mx_var)
			{
				const dt_int32 iy_0 = max(iy_i+n_k0, 0);
				const dt_int32 iy_e = min(iy_i+n_ke, ny_i);

				T v_mean = 0;
				T v_var = 0;
				for(auto iy = iy_0; iy < iy_e; iy++)
				{
					const auto y = mx_i(iy, ix_i);
					v_mean += y;
					v_var += y*y;
				}

				v_mean = v_mean/T(iy_e-iy_0);
				v_var = v_var/T(iy_e-iy_0) - v_mean*v_mean;

				const auto ind = mx_i.sub_2_ind(iy_i, ix_i);
				mx_mean[ind] = v_mean;
				mx_var[ind] = v_var;
			};

			fcn_stream_exec_xd_krn<edim_2>(pstream, nx_i, ny_i, krn_mean_var, mx_mean, mx_var);

			/* apply wiener filter */
			auto krn_fltr_wiener = [ny_i, &mx_i, &mx_o](const dt_int32& ix_i, Vctr_cpu<T>& mx_mean, Vctr_cpu<T>& mx_var)
			{
				KS<T> v2 = T(0);
				for(auto iy = 0; iy < ny_i; iy++)
				{
					v2 += mx_var(iy, ix_i);
				}
				v2 /= T(ny_i);

				for(auto iy = 0; iy < ny_i; iy++)
				{
					auto ind = mx_i.sub_2_ind(iy, ix_i);
					mx_o[ind] = mx_mean[ind] + ::fmax(T(0), mx_var[ind]-v2.m_sum)*(mx_i[ind]-mx_mean[ind])/::fmax(mx_var[ind], v2.m_sum);
				}
			};

			fcn_stream_exec_xd_krn<edim_1>(pstream, nx_i, krn_fltr_wiener, mx_mean, mx_var);
		}

		/* 2d */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_wiener_2d(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			if (n_kr<=0)
			{
				memcpy_cpu_cpu(mx_o.m_data, mx_i.m_data, mx_i.size_64());
				return;
			}

			const dt_int32 n_k0 = -n_kr;
			const dt_int32 n_ke = n_kr+1;
			const dt_int32 nx_i = mx_i.s1_32();
			const dt_int32 ny_i = mx_i.s0_32();

			const T radius_2 = ::square(T(n_kr)) + epsilon_rel<T>();

			Vctr_cpu<T> mx_mean(mx_i.shape(), T(0));
			Vctr_cpu<T> mx_var(mx_i.shape(), T(0));

			/* get local mean and var */
			auto krn_mean_var = [n_k0, n_ke, radius_2, nx_i, ny_i, &mx_i](const dt_int32& ix_i, const dt_int32& iy_i, Vctr_cpu<T>& mx_mean, Vctr_cpu<T>& mx_var, KS<T>& v2)
			{
				const dt_int32 ix_0 = max(ix_i+n_k0, 0);
				const dt_int32 ix_e = min(ix_i+n_ke, nx_i);

				const dt_int32 iy_0 = max(iy_i+n_k0, 0);
				const dt_int32 iy_e = min(iy_i+n_ke, ny_i);

				T v_mean = 0;
				T v_var = 0;
				for(auto ix = ix_0; ix < ix_e; ix++)
				{					
					for(auto iy = iy_0; iy < iy_e; iy++)
					{
						if (R_2d<T>(ix-ix_i, iy-iy_i).norm_2()<=radius_2)
						{
							const auto y = mx_i(iy, ix);
							v_mean += y;
							v_var += y*y;
						}
					}
				}

				v_mean = v_mean/T((iy_e-iy_0)*(ix_e-ix_0));
				v_var = v_var/T((iy_e-iy_0)*(ix_e-ix_0)) - v_mean*v_mean;

				const auto ind = mx_i.sub_2_ind(iy_i, ix_i);
				mx_mean[ind] = v_mean;
				mx_var[ind] = v_var;

				v2 += v_var;
			};

			auto v2 = fcn_stream_exec_xd_krn_reduce<edim_2>(pstream, nx_i, ny_i, std::plus<KS<T>>(), T(0), krn_mean_var, mx_mean, mx_var);
			v2 /= T(nx_i*ny_i);

			/* apply wiener filter */
			auto krn_fltr_wiener = [nx_i, ny_i, v2, &mx_i, &mx_o](const dt_int32& ix, const dt_int32& iy, Vctr_cpu<T>& mx_mean, Vctr_cpu<T>& mx_var)
			{
				auto ind = mx_i.sub_2_ind(iy, ix);
				mx_o[ind] = mx_mean[ind] + ::fmax(T(0), mx_var[ind]-v2)*(mx_i[ind]-mx_mean[ind])/::fmax(mx_var[ind], v2);
			};

			fcn_stream_exec_xd_krn<edim_2>(pstream, nx_i, ny_i, krn_fltr_wiener, mx_mean, mx_var);
		}
	}

	/* filter mean */
	namespace cpu_fcn_image
	{
		/* 1d */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_mean_1d(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream= nullptr)
		{
			using T = Value_type<TVctr>;

			if (n_kr<=0)
			{
				memcpy_cpu_cpu(mx_o.m_data, mx_i.m_data, mx_i.size_64());
				return;
			}

			const dt_int32 n_k0 = -n_kr;
			const dt_int32 n_ke = n_kr+1;
			const dt_int32 nx_i = mx_i.size_32();

			auto thr_fltr_mean = [n_k0, n_ke, nx_i, &mx_i, &mx_o](const iThread_Rect_1d& ithread)
			{
				std::vector<T> vctr_t(n_ke-n_k0);

				for(auto ind = ithread.ind_0; ind < ithread.ind_e; ind++)
				{
					const dt_int32 ix_0 = fcn_max(ind+n_k0, 0);
					const dt_int32 ix_e = fcn_min(ind+n_ke, nx_i);

					T mv = T(0);
					dt_int32 ic = 0;
					for(auto ix = ix_0; ix < ix_e; ix++)
					{
						mv += mx_i[ix];
						ic++;
					}

					mx_o[ind] = mv/T(ic);
				}
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, nx_i, thr_fltr_mean);
		}

		/* 2d */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_mean_2d(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream= nullptr)
		{
			using T = Value_type<TVctr>;

			if (n_kr<=0)
			{
				memcpy_cpu_cpu(mx_o.m_data, mx_i.m_data, mx_i.size_64());
				return;
			}

			const dt_int32 n_k0 = -n_kr;
			const dt_int32 n_ke = n_kr+1;
			const dt_int32 nx_i = mx_i.s1_32();
			const dt_int32 ny_i = mx_i.s0_32();

			const T radius_2 = ::square(T(n_kr)) + epsilon_rel<T>();

			auto thr_fltr_mean = [n_k0, n_ke, radius_2, nx_i, ny_i, &mx_i, &mx_o](const iThread_Rect_2d& ithread)
			{
				std::vector<T> vctr_t(::square(n_ke-n_k0));

				for(auto ix = ithread.ix_0; ix < ithread.ix_e; ix++)
				{
					for(auto iy = ithread.iy_0; iy < ithread.iy_e; iy++)
					{
						const dt_int32 ix_0 = max(ix+n_k0, 0);
						const dt_int32 ix_e = min(ix+n_ke, nx_i);

						const dt_int32 iy_0 = max(iy+n_k0, 0);
						const dt_int32 iy_e = min(iy+n_ke, ny_i);

						T mv = T(0);
						dt_int32 ic = 0;
						for(auto ix_t = ix_0; ix_t < ix_e; ix_t++)
						{
							for(auto iy_t = iy_0; iy_t < iy_e; iy_t++)
							{
								if (R_2d<T>(ix_t-ix, iy_t-iy).norm_2()<=radius_2)
								{
									mv += mx_i(iy_t, ix_t);
									ic++;
								}
							}
						}

						mx_o(iy, ix) = mv/T(ic);
					}
				}
			};

			fcn_stream_exec_xd_fcn<edim_2>(pstream, nx_i, ny_i, thr_fltr_mean);
		}
	}

	/* filter median */
	namespace cpu_fcn_image
	{
		/* 1d */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_median_1d(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream= nullptr)
		{
			using T = Value_type<TVctr>;

			if (n_kr<=0)
			{
				memcpy_cpu_cpu(mx_o.m_data, mx_i.m_data, mx_i.size_64());
				return;
			}

			const dt_int32 n_k0 = -n_kr;
			const dt_int32 n_ke = n_kr+1;
			const dt_int32 nx_i = mx_i.size_32();

			auto thr_fltr_median = [n_k0, n_ke, nx_i, &mx_i, &mx_o](const iThread_Rect_1d& ithread)
			{
				std::vector<T> vctr_t(n_ke-n_k0);

				for(auto ind = ithread.ind_0; ind < ithread.ind_e; ind++)
				{
					const dt_int32 ix_0 = fcn_max(ind+n_k0, 0);
					const dt_int32 ix_e = fcn_min(ind+n_ke, nx_i);

					dt_int32 ic = 0;
					for(auto ix = ix_0; ix < ix_e; ix++)
					{
						vctr_t[ic++] = mx_i[ix];
					}

					auto median = vctr_t.begin() + ic/2;
					std::nth_element(vctr_t.begin(), median, vctr_t.begin() + ic);
					mx_o[ind] = *median;
				}
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, nx_i, thr_fltr_median);
		}

		/* 2d by col */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_median_2d_bc(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream= nullptr)
		{
			using T = Value_type<TVctr>;

			if (n_kr<=0)
			{
				memcpy_cpu_cpu(mx_o.m_data, mx_i.m_data, mx_i.size_64());
				return;
			}

			const dt_int32 n_k0 = -n_kr;
			const dt_int32 n_ke = n_kr+1;
			const dt_int32 nx_i = mx_i.s1_32();
			const dt_int32 ny_i = mx_i.s0_32();

			auto thr_fltr_median = [n_k0, n_ke, ny_i, &mx_i, &mx_o](const iThread_Rect_2d& ithread)
			{
				std::vector<T> vctr_t(n_ke-n_k0);

				for(auto ix = ithread.ix_0; ix < ithread.ix_e; ix++)
				{
					for(auto iy = ithread.iy_0; iy < ithread.iy_e; iy++)
					{
						const dt_int32 iy_0 = fcn_max(iy+n_k0, 0);
						const dt_int32 iy_e = fcn_min(iy+n_ke, ny_i);

						dt_int32 ic = 0;
						for(auto iy_t = iy_0; iy_t < iy_e; iy_t++)
						{
							vctr_t[ic++] = mx_i(iy_t, ix);
						}

						auto median = vctr_t.begin() + ic/2;
						std::nth_element(vctr_t.begin(), median, vctr_t.begin() + ic);
						mx_o(iy, ix) = *median;
					}
				}
			};

			fcn_stream_exec_xd_fcn<edim_2>(pstream, nx_i, ny_i, thr_fltr_median);
		}

		/* 2d */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_median_2d(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream= nullptr)
		{
			using T = Value_type<TVctr>;

			if (n_kr<=0)
			{
				memcpy_cpu_cpu(mx_o.m_data, mx_i.m_data, mx_i.size_64());
				return;
			}

			const dt_int32 n_k0 = -n_kr;
			const dt_int32 n_ke = n_kr+1;
			const dt_int32 nx_i = mx_i.s1_32();
			const dt_int32 ny_i = mx_i.s0_32();

			const T radius_2 = ::square(T(n_kr)) + epsilon_rel<T>();

			auto thr_fltr_median = [n_k0, n_ke, radius_2, nx_i, ny_i, &mx_i, &mx_o](const iThread_Rect_2d& ithread)
			{
				std::vector<T> vctr_t(::square(n_ke-n_k0));

				for(auto ix = ithread.ix_0; ix < ithread.ix_e; ix++)
				{
					for(auto iy = ithread.iy_0; iy < ithread.iy_e; iy++)
					{
						const dt_int32 ix_0 = max(ix+n_k0, 0);
						const dt_int32 ix_e = min(ix+n_ke, nx_i);

						const dt_int32 iy_0 = max(iy+n_k0, 0);
						const dt_int32 iy_e = min(iy+n_ke, ny_i);

						dt_int32 ic = 0;
						for(auto ix_t = ix_0; ix_t < ix_e; ix_t++)
						{
							for(auto iy_t = iy_0; iy_t < iy_e; iy_t++)
							{
								if (R_2d<T>(ix_t-ix, iy_t-iy).norm_2()<=radius_2)
								{
									vctr_t[ic++] = mx_i(iy_t, ix_t);
								}
							}
						}

						auto median = vctr_t.begin() + ic/2;
						std::nth_element(vctr_t.begin(), median, vctr_t.begin() + ic);
						mx_o(iy, ix) = *median;
					}
				}
			};

			fcn_stream_exec_xd_fcn<edim_2>(pstream, nx_i, ny_i, thr_fltr_median);
		}
	}

	/* filter median by position */
	namespace cpu_fcn_image
	{
		/* 1d */
		template <class TVctr, class TVctr_idx>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr, TVctr_idx, void>
		fltr_median_1d_pos(TVctr& mx_i, dt_int32 n_kr, TVctr_idx& vctr_idx, TVctr& mx_o, Stream_cpu* pstream= nullptr)
		{
			using T = Value_type<TVctr>;
			memcpy_cpu_cpu(mx_o.m_data, mx_i.m_data, mx_i.size_64());

			if (n_kr<=0)
			{
				return;
			}

			const dt_int32 n_k0 = -n_kr;
			const dt_int32 n_ke = n_kr+1;
			const dt_int32 nx_i = mx_i.size_32();

			auto thr_fltr_median = [n_k0, n_ke, nx_i, &mx_i, &vctr_idx, &mx_o](const iThread_Rect_1d& ithread)
			{
				std::vector<T> vctr_t(n_ke-n_k0);

				for(auto ind = ithread.ind_0; ind < ithread.ind_e; ind++)
				{
					const dt_int32 idx = vctr_idx[ind];

					const dt_int32 ix_0 = fcn_max(idx+n_k0, 0);
					const dt_int32 ix_e = fcn_min(idx+n_ke, nx_i);

					dt_int32 ic = 0;
					for(auto ix = ix_0; ix < ix_e; ix++)
					{
						vctr_t[ic++] = mx_i[ix];
					}

					auto median = vctr_t.begin() + ic/2;
					std::nth_element(vctr_t.begin(), median, vctr_t.begin() + ic);
					mx_o[idx] = *median;
				}
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, vctr_idx.size_32(), thr_fltr_median);
		}
	
		// 2d for specific points
		template <class TVctr, class TVctr_idx>
		enable_if_vctr_cpu_r_2d_and_vctr_cpu<TVctr_idx, TVctr, void>
		fltr_median_2d_pos(TVctr& mx_i, dt_int32 n_kr, TVctr_idx& vctr_idx, TVctr& mx_o, Stream_cpu* pstream)
		{
			using T = Value_type<TVctr>;
			memcpy_cpu_cpu(mx_o.m_data, mx_i.m_data, mx_i.size_64());

			if (n_kr<=0)
			{
				return;
			}

			const dt_int32 n_k0 = -n_kr;
			const dt_int32 n_ke = n_kr+1;
			const dt_int32 nx_i = mx_i.s1_32();
			const dt_int32 ny_i = mx_i.s0_32();

			const T radius_2 = ::square(T(n_kr)) + epsilon_rel<T>();

			auto thr_fltr_median = [n_k0, n_ke, radius_2, nx_i, ny_i, &mx_i, &vctr_idx, &mx_o](const iThread_Rect_1d& ithread)
			{
				std::vector<T> vctr_t(::square(n_ke-n_k0));

				for(auto ind = ithread.ind_0; ind < ithread.ind_e; ind++)
				{
					const R_2d<dt_int32> r_2d = vctr_idx[ind];

					const dt_int32 ix_0 = max(r_2d.x+n_k0, 0);
					const dt_int32 ix_e = min(r_2d.x+n_ke, nx_i);

					const dt_int32 iy_0 = max(r_2d.y+n_k0, 0);
					const dt_int32 iy_e = min(r_2d.y+n_ke, ny_i);

					dt_int32 ic = 0;
					for(auto ix_t = ix_0; ix_t < ix_e; ix_t++)
					{
						for(auto iy_t = iy_0; iy_t < iy_e; iy_t++)
						{
							if (R_2d<T>(ix_t-r_2d.x, iy_t-r_2d.y).norm_2()<=radius_2)
							{
								vctr_t[ic++] = mx_i(iy_t, ix_t);
							}
						}
					}

					auto median = vctr_t.begin() + ic/2;
					std::nth_element(vctr_t.begin(), median, vctr_t.begin() + ic);
					mx_o(r_2d.y, r_2d.x) = *median;
				}
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, vctr_idx.size_32(), thr_fltr_median);
		}
	}

	/* filter poisson noise */
	namespace cpu_fcn_image
	{
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_poiss_nois_1d(TVctr& mx_i, dt_int32 nkr_w, dt_int32 nkr_m, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			if ((nkr_w == 0) && (nkr_m == 0))
			{
				memcpy_cpu_cpu(mx_o.m_data, mx_i.m_data, mx_i.size_64());
				return;
			}

			Vctr_cpu<T> mx_tmp(mx_i.shape());

			fcn_anscombe(mx_i, mx_tmp, pstream);
			fltr_wiener_1d(mx_tmp, nkr_w, mx_o, pstream);
			fltr_median_1d(mx_o, nkr_m, mx_tmp, pstream);
			fcn_anscombe_inv(mx_tmp, mx_o, pstream);
		}

		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_poiss_nois_2d_bc(TVctr& mx_i, dt_int32 nkr_w, dt_int32 nkr_m, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			if ((nkr_w == 0) && (nkr_m == 0))
			{
				memcpy_cpu_cpu(mx_o.m_data, mx_i.m_data, mx_i.size_64());
				return;
			}

			Vctr_cpu<T> mx_tmp(mx_i.shape());

			fcn_anscombe(mx_i, mx_tmp, pstream);
			fltr_wiener_2d_bc(mx_tmp, nkr_w, mx_o, pstream);
			fltr_median_2d_bc(mx_o, nkr_m, mx_tmp, pstream);
			fcn_anscombe_inv(mx_tmp, mx_o, pstream);
		}

		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fltr_poiss_nois_2d(TVctr& mx_i, dt_int32 nkr_w, dt_int32 nkr_m, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			if ((nkr_w == 0) && (nkr_m == 0))
			{
				memcpy_cpu_cpu(mx_o.m_data, mx_i.m_data, mx_i.size_64());
				return;
			}

			Vctr_cpu<T> mx_tmp(mx_i.shape());

			fcn_anscombe(mx_i, mx_tmp, pstream);
			fltr_wiener_2d(mx_tmp, nkr_w, mx_o, pstream);
			fltr_median_2d(mx_o, nkr_m, mx_tmp, pstream);
			fcn_anscombe_inv(mx_tmp, mx_o, pstream);
		}
	}


	/* peak signal to noise ratio */
	namespace cpu_fcn_image
	{
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type<TVctr>>
		fcn_get_psnr(TVctr& mx_i, TVctr& mx_r, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto var_r = fcn_variance(mx_r, pstream);
			Vctr_cpu<T> mx_n(mx_r.size());
			fcn_sub(mx_i, mx_r, mx_n, pstream);

			auto var_n = fcn_variance(mx_n, pstream);

			return var_n/var_r;
		}

		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type<TVctr>>
		fcn_get_psnr(TVctr& mx_i, dt_int32 nkr_w, dt_int32 nkr_m, Stream_cpu* pstream = nullptr)	
		{
			auto mx_r = fltr_poiss_nois_2d(mx_i, nkr_w, nkr_m, pstream);

			auto var_r = fcn_variance(mx_r, pstream);
			auto& mx_n = mx_r;
			fcn_sub(mx_i, mx_r, mx_n, pstream);

			auto var_n = fcn_variance(mx_n, pstream);

			return var_n/var_r;
		}
	}

	/* others */
	namespace cpu_fcn_image
	{
		//// scale_image
		//template <class TVctr>
		//TVctr scale_image_mean(Stream<edev_cpu>& stream, dt_int32 ny_i, dt_int32 nx_i, TVctr& Im_i, Value_type<TVctr> shrink_factor, dt_int32& ny_o, dt_int32& nx_o)
		//{
		//	using T = Value_type<TVctr>;
		//	TVctr Im_o;

		//	dt_int32 nkr = max(dt_int32(1.0/(2.0*shrink_factor)), 1);
		//	dt_int32 nk0 = -nkr;
		//	dt_int32 nke = nkr+1;

		//	nx_o = max(dt_int32(::floor(nx_i*shrink_factor)), 1);
		//	ny_o = max(dt_int32(::floor(ny_i*shrink_factor)), 1);

		//	if ((nx_i == nx_o) && (ny_i == ny_o))
		//	{
		//		Im_o.assign(Im_i.begin(), Im_i.end());
		//		return Im_o;
		//	}
		//	Im_o.resize(nx_o*ny_o);

		//	auto krn_sc_image = [&](const dt_int32& ix_i, const dt_int32& iy_i, TVctr& Im_i, TVctr& Im_o)
		//	{
		//		auto ix_t = static_cast<dt_int32>(ix_i/shrink_factor);
		//		auto iy_t = static_cast<dt_int32>(iy_i/shrink_factor);

		//		dt_int32 ix_0 = max(ix_t+nk0, 0);
		//		dt_int32 ix_e = min(ix_t+nke, nx_i);

		//		dt_int32 iy_0 = max(iy_t+nk0, 0);
		//		dt_int32 iy_e = min(iy_t+nke, ny_i);

		//		T sum = 0;
		//		for(auto ix = ix_0; ix < ix_e; ix++)
		//		{
		//				for(auto iy = iy_0; iy < iy_e; iy++)
		//				{
		//					sum += Im_i[ix*ny_i+iy];
		//				}
		//		}

		//		Im_o[ix_i*ny_o+iy_i] = sum/((ix_e-ix_0)*(iy_e-iy_0));
		//	};

		//	auto thr_sc_image = [&](const iThread_Rect_2d& ithread)
		//	{
		//		cpu_detail::matrix_iter(ithread, krn_sc_image, Im_i, Im_o);
		//	};

		//	stream.set_n_stream_act(nx_o);
		//	stream.set_grid(nx_o, ny_o);
		//	stream.exec(thr_sc_image);

		//	return Im_o;
		//}

		//// copy image
		//template <class TVctr>
		//TVctr cpy_image(Stream<edev_cpu>& stream, dt_int32 ny_src, dt_int32 nx_src, TVctr& Im_src, dt_int32 iy_0, dt_int32 ix_0, dt_int32 iy_e, dt_int32 ix_e, dt_int32 ny_dst, dt_int32 nx_dst)
		//{
		//	TVctr Im_dst(nx_dst*ny_dst);

		//	dt_int32 nx = min(nx_src, nx_dst);
		//	ix_e = (ix_e<nx)?ix_e+1:nx;
		//	nx = ix_e-ix_0;

		//	dt_int32 ny = min(ny_src, ny_dst);
		//	iy_e = (iy_e<ny)?iy_e+1:ny;
		//	ny = iy_e-iy_0;

		//	auto thr_cpy_image = [&](const iThread_Rect_2d& ithread)
		//	{
		//		for(auto ix = ithread.ix_0; ix < ithread.ix_e; ix++)
		//		{
		//				for(auto iy = ithread.iy_0; iy < ithread.iy_e; iy++)
		//				{
		//					Im_dst[ix*ny_dst+iy] = Im_src[(ix+ix_0)*ny_src+(iy+iy_0)];
		//				}
		//		}
		//	};

		//	stream.set_n_stream_act(nx);
		//	stream.set_grid(nx, ny);
		//	stream.exec(thr_cpy_image);

		//	return Im_dst;
		//}
	}
}