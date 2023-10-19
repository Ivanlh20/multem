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

#ifndef PARTICLE_FCNS_H
	#define PARTICLE_FCNS_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include <numeric>
	#include <vector>
	#include <tuple>
	#include <functional>
	#include <algorithm>

	#include "macros.h"
	#include "math_mt.h"
	#include "type_traits_gen.h"
	#include "r_2d.h"
	#include "r_3d.h"
	#include "vctr_cpu.h"

	namespace mt
	{
		// calculate the mean position using xy cordinates and it frequency
		template <class TVctr>
		enable_if_vctr_cpu_r_2d<TVctr, Vctr_r_3d_cpu<Value_type_r<TVctr>>>
		fcn_xy_2_xyc(TVctr& vr_2d, Value_type_r<TVctr> radius)
		{
			using T = Value_type_r<TVctr>;
			using ST = Size_type<TVctr>;

			const ST n_vr_2d = vr_2d.size();

			const T r2_max = ::square(radius);

			Vctr_r_3d_cpu<T> vr_3d;
			vr_3d.reserve(n_vr_2d);

			std::vector<dt_bool> bb(n_vr_2d, true);

			for(ST ip = 0; ip < n_vr_2d; ip++)
			{
				if (bb[ip])
				{
					const auto r_i = vr_2d[ip];

					ST ic = 0;
					R_2d<T> p_m = 0;

					for(ST ip_s = 0; ip_s < n_vr_2d; ip_s++)
					{
						if (bb[ip_s])
						{
							if (norm_2(vr_2d[ip_s]-r_i) < r2_max)
							{
								p_m += vr_2d[ip_s];
								bb[ip_s] = false;
								ic++;
							}
						}
					}
					p_m /= T(ic);

					vr_3d.push_back(R_3d<T>(p_m.x, p_m.y, T(ic)));
				}
			}

			return vr_3d;
		}

		// replace xy positions by its mean positions using a xy reference
		template <class TVctr>
		enable_if_vctr_cpu_r_2d<TVctr, void>
		fcn_xy_2_xym(const TVctr& vr_2d_r, Value_type_r<TVctr> radius, TVctr& vr_2d_io)
		{
			using T = Value_type_r<TVctr>;
			using ST = Size_type<TVctr>;

			const ST n_vr_2d = vr_2d_io.size();
			const ST n_vr_2d_r = vr_2d_r.size();
			//const ST n_vr_2d_r = min(vr_2d_r.size(), n_vr_2d);

			const T r2_max = ::square(radius);

			std::vector<dt_bool> bb(n_vr_2d, true);
			std::vector<ST> ind_s(n_vr_2d, ST(0));

			for(ST ip = 0; ip < n_vr_2d_r; ip++)
			{
				const auto r_c = vr_2d_r[ip];

				ST ic = 0;
				R_2d<T> p_m = 0;

				for(ST ip_s = 0; ip_s < n_vr_2d; ip_s++)
				{
					if (bb[ip_s])
					{
						if (norm_2(vr_2d_io[ip_s]-r_c) < r2_max)
						{
							p_m += vr_2d_io[ip_s];
							bb[ip_s] = false;
							ind_s[ic] = ip_s;
							ic++;
						}
					}
				}
				p_m /= T(ic);

				for(ST ip_s = 0; ip_s < ic; ip_s++)
				{
					vr_2d_io[ind_s[ip_s]] = p_m;
				}
			}
		}

		// match one to one xy positions: vr_2d_r: array of R_2d and vr_2d: array of R_2d
		template <class TVctr>
		enable_if_vctr_cpu_r_2d<TVctr, void>
		fcn_xy_match_one_2_one(const TVctr& vr_2d_r, TVctr& vr_2d, Vctr_cpu<dt_uint64>& ind_match)
		{
			if (vr_2d_r.size() != vr_2d.size())
				return;

			using T = Value_type<TVctr>;
			using U = Value_type_r<TVctr>;
			using ST = Size_type<TVctr>;

			const ST n_vr_2d = vr_2d.size();

			using typ_ind_d2 = std::tuple<ST, ST, U, T>;

			std::vector<typ_ind_d2> data(n_vr_2d, std::make_tuple(ST(), ST(), U(), T()));

			// find minimum distance for each reference position
			for(ST ip = 0; ip < n_vr_2d; ip++)
			{
				const auto r_c = vr_2d_r[ip];

				ST ip_s_min = 0;
				U d2_s_min = norm_2(vr_2d[ip_s_min]-r_c);

				for(ST ip_s = 1; ip_s < n_vr_2d; ip_s++)
				{
					auto d2_ip = norm_2(vr_2d[ip_s]-r_c);

					if (d2_ip < d2_s_min)
					{
						d2_s_min = d2_ip;
						ip_s_min =ip_s;
					}
				}

				data[ip] = std::make_tuple(ip, ip_s_min, d2_s_min, vr_2d[ip_s_min]);
			}

			// sort by distance
			std::sort(data.begin(), data.end(), [](const typ_ind_d2 &a, const typ_ind_d2 &b)
			{ return (std::get<1>(a) < std::get<1>(b)) || 
				((std::get<1>(a) == std::get<1>(b)) && (std::get<2>(a) < std::get<2>(b)));
			});

			// std::sort(data.begin(), data.end(), [](const typ_ind_d2 &a, const typ_ind_d2 &b){ return std::get<2>(a) < std::get<2>(b); });

			// get unique indices
			auto last = std::unique(data.begin(), data.end(), [](const typ_ind_d2 &a, const typ_ind_d2 &b){ return std::get<1>(a) == std::get<1>(b); });
			ST n_uniq = std::distance(data.begin(), last);
			ST n_left = n_vr_2d - n_uniq;

			// return if there is not duplicate indices
			if (n_left==0)
			{
				for(ST ip = 0; ip < n_vr_2d; ip++)
				{
					ind_match[std::get<0>(data[ip])] = std::get<1>(data[ip]);
					vr_2d[std::get<0>(data[ip])] = std::get<3>(data[ip]);
				}

				return;
			}
		
			// select unique matches
			std::vector<dt_bool> bb_r(n_vr_2d, true);
			std::vector<dt_bool> bb(n_vr_2d, true);

			for(ST ip = 0; ip < n_uniq; ip++)
			{
				bb_r[std::get<0>(data[ip])] = false;
				bb[std::get<1>(data[ip])] = false;
			}

			// get duplicated indices
			std::vector<ST> ind_r;
			ind_r.reserve(n_left);

			std::vector<ST> ind;
			ind.reserve(n_left);

			for(ST ip = 0; ip < n_vr_2d; ip++)
			{
				if (bb_r[ip])
					ind_r.push_back(ip);

				if (bb[ip])
					ind.push_back(ip);
			}

			// match duplicate indices
			for(ST ip = 0; ip < n_left; ip++)
			{
				const auto r_c = vr_2d_r[ind_r[ip]];

				// initial index have to be calculated for each iteration
				ST ip_s_g_min = 0;

				for(ST ip_s = 0; ip_s < n_left; ip_s++)
				{
					if (bb[ind[ip_s]])
					{
						ip_s_g_min = ind[ip_s];
						break;
					}
				}

				U d2 = norm_2(vr_2d[ip_s_g_min]-r_c);

				for(ST ip_s = 1; ip_s < n_left; ip_s++)
				{
					auto ip_s_g = ind[ip_s];

					if (bb[ip_s_g])
					{
						auto d2_ip = norm_2(vr_2d[ip_s_g]-r_c);

						if (d2_ip < d2)
						{
							d2 = d2_ip;
							ip_s_g_min =ip_s_g;
						}
					}
				}

				bb[ip_s_g_min] = false;
				data[n_uniq + ip] = std::make_tuple(ind_r[ip], ip_s_g_min, d2, vr_2d[ip_s_g_min]);
			}

			for(ST ip = 0; ip < n_vr_2d; ip++)
			{
				ind_match[std::get<0>(data[ip])] = std::get<1>(data[ip]);
				vr_2d[std::get<0>(data[ip])] = std::get<3>(data[ip]);
			}
		}

		/***************************************************************************************/
		// radial distribution function
		template <class TVctr_r_nd, class TVctr>
		enable_if_vctr_cpu_r_nd_and_vctr_cpu<TVctr_r_nd, TVctr, void>
		fcn_rdf(TVctr_r_nd& vr_nd, Value_type<TVctr> r_max, dt_int32 nr, TVctr& r_o, TVctr& rdf_o)
		{
			using T = Value_type<TVctr>;
			using ST = Size_type<TVctr>;

			const T dr = r_max/T(nr);		

			for(ST ir = 0; ir < nr; ir++)
			{
				r_o[ir] = ir*dr;
				rdf_o[ir] = T(0);
			}

			const auto r2_max = ::square(r_max);
				
			const ST n_vr_nd = vr_nd.size();
			for(ST ixy = 0; ixy < n_vr_nd; ixy++)
			{
				const auto r_i = vr_nd[ixy];
				for(ST ixy_s = 0; ixy_s < n_vr_nd; ixy_s++)
				{
					auto d2 = norm_2(vr_nd[ixy_s]-r_i);
					if ((d2 < r2_max) && (ixy != ixy_s))
					{
						const auto ir = static_cast<ST>(::floor(::sqrt(d2)/dr));
						rdf_o[ir] += 1;
					}
				}
			}
			if (typeid(Value_type<TVctr_r_nd>) == typeid(R_2d<T>))
			{
				for(auto ir = 1; ir < nr; ir++)
				{
					rdf_o[ir] = rdf_o[ir]/(c_pi<T>*r_o[ir]*dr);
				}
			}
			else if (typeid(Value_type<TVctr_r_nd>) == typeid(R_3d<T>))
			{
				for(auto ir = 1; ir < nr; ir++)
				{
					rdf_o[ir] = rdf_o[ir]/(T(4)*c_pi<T>*::square(r_o[ir])*dr);
				}
			}
		}

		/***************************************************************************************/
		template <class TVctr, class TFcn>
		enable_if_vctr_cpu_r_nd<TVctr, std::tuple<Value_type_r<TVctr>, Size_type<TVctr>, Size_type<TVctr>>>
		fcn_vctr_r_nd_dist(const TVctr& vr_nd, TFcn fcn)
		{
			using T = Value_type_r<TVctr>;
			using ST = Size_type<TVctr>;

			const ST n_vr_nd = vr_nd.size();

			if (n_vr_nd == 1)
			{
				return std::make_tuple(T(0), ST(0), ST(0));
			}

			ST idx_0 = 0;
			ST idx_e = 1;
			T d2_op = norm_2(vr_nd[0]-vr_nd[1]);

			for(ST ixy = 0; ixy < n_vr_nd; ixy++)
			{
				const auto r_i = vr_nd[ixy];
				//for (ST ixy_s = 0; ixy_s < n_vr_nd; ixy_s++) 	// check out this line because it can be optmize by 
				for(ST ixy_s = ixy+1; ixy_s < n_vr_nd; ixy_s++)
				{
					auto d2 = norm_2(vr_nd[ixy_s]-r_i);
					if ((ixy != ixy_s) && fcn(d2, d2_op))
					{
						idx_0 = ixy;
						idx_e = ixy_s;
						d2_op = d2;
					}
				}
			}

			return std::make_tuple(::sqrt(d2_op), idx_0, idx_e);
		}

		template <class TVctr>
		enable_if_vctr_cpu_r_nd<TVctr, Value_type_r<TVctr>>
		fcn_min_dist(const TVctr& vr_nd)
		{
			using T = Value_type_r<TVctr>;

			auto dist_t = fcn_vctr_r_nd_dist(vr_nd, std::less<T>());

			return std::get<0>(dist_t);
		}

		template <class TVctr>
		enable_if_vctr_cpu_r_nd<TVctr, Value_type_r<TVctr>>
		fcn_max_dist(const TVctr& vr_nd)
		{
			using T = Value_type_r<TVctr>;

			auto dist_t = fcn_vctr_r_nd_dist(vr_nd, std::greater<T>());

			return std::get<0>(dist_t);
		}
	}

#endif