/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef EVAL_FIT_GAUSSIANS_H
#define EVAL_FIT_GAUSSIANS_H

#include <random>
#include <vector>
#include <deque>
#include <map>

#include "math.cuh"
#include "type_traits_gen.cuh"
#include "types.cuh"
#include "lapack.hpp"

 namespace mt
{
	template <class T, eDev Dev>
	class Fit_one_rad_fcn_2d
	{
		public:
			using T_r = T;
			using TVctr_r = vector<T>;
			using TRegion_r = Region_Rad_2d<TVctr_r>;
			using size_type = dt_uint64;

			static const eDev device = Dev;

			Fit_one_rad_fcn_2d():n_coef(6), n_coef_nl(4), n_coef_l(2), lambda(1), 
			lambda_fu(8), lambda_du(5), ee_ee(1e-5), ee_G_max(1e-5), ee_coef(1e-7) {};

			template <class TVctr>
			TVctr_r fit(TVctr& Ixy_i, R_2d<T> pos_i, T sigma_i, T radius_i)
			{
				thrust::copy(Ixy_i.begin(), Ixy_i.end(), Ixy.begin());

				// select circular region
				select_cir_reg(Ixy, pos_i, radius_i, region);

				// set initial coefficients and limits
				set_init_min_max_coef_0(pos_i, region.Ixy_sc, sigma_i, radius_i, coef, coef_min, coef_max);

				const dt_int32 n_it_o = 4;
				for(auto it = 0; it < n_it_o; it++)
				{
					const dt_int32 n_it_i = (it == n_it_o-1)?35:15;

					T ee_b = fit_region(region, coef, coef_min, coef_max, lambda, n_it_i);

					if (it < n_it_o - 1)
					{
						// select elliptical region
						select_ellip_reg(Ixy, coef, 0.5, region);

						// set coefficients and limits
						set_min_max_coef(coef, coef_min, coef_max);
					}
					lambda = lambda/2;
				}

				return coef;
			}

		protected:
			T lambda;
			T lambda_fu;
			T lambda_du;

			dt_int32 n_iter;
			T ee_ee;
			T ee_G_max;
			T ee_coef;

			dt_int32 n_coef;
			dt_int32 n_coef_nl;
			dt_int32 n_coef_l;

			std::mt19937_64 gen_c;
			std::uniform_real_distribution<T> rand_c;

			Grid_2d<T> grid_2d;

			TVctr_r Ixy;
			TVctr_r dIxy;
			TVctr_r J;
			TVctr_r A;
			TVctr_r A_lambda;
			TVctr_r b;

			TVctr_r coef_min;
			TVctr_r coef_max;
			TVctr_r coef;
			TVctr_r d_coef;
			TVctr_r coef_t;

			lapack::MNS_SVD<T> mns_svd;

			TRegion_r region;

			/***************************************************************************************/
			// set initial coefficients and limits
			virtual void set_init_min_max_coef_0(R_2d<T> p, T Ip_max, T sigma, T radius, TVctr_r &coef, TVctr_r &coef_min, TVctr_r &coef_max) = 0;

			// set coefficients and limits
			virtual void set_min_max_coef(TVctr_r &coef, TVctr_r &coef_min, TVctr_r &coef_max) = 0;

			// forward transformation coefficients
			virtual void fw_transformation_coef(TRegion_r &region, TVctr_r &coef) = 0;

			// backward transformation coefficients
			virtual void bk_transformation_coef(TRegion_r &region, TVctr_r &coef) =0;

			// set order coef
			virtual void set_order_coef(TVctr_r &coef) = 0;

			// calculate error, jacobian, and dIxy
			virtual T cal_ee_dIxy_J(TRegion_r &region, TVctr_r &coef, TVctr_r &dIxy, TVctr_r &J) = 0;

			/***************************************************************************************/
			T fit_region(TRegion_r &region, TVctr_r &coef, TVctr_r coef_min, TVctr_r coef_max, T lambda, dt_int32 n_iter)
			{
				const dt_int32 nxy = region.size();

				dt_bool bb_A_b = true;
				dt_bool bb_shrink_lambda = false;

				T G_max;
				T ee_coef_max;
				T ee_t;

				// forward transformation coefficients
				fw_transformation_coef(region, coef);
				fw_transformation_coef(region, coef_min);
				fw_transformation_coef(region, coef_max);

				// set order coef
				set_order_coef(coef);

				// call ee
				T ee = cal_ee_dIxy_J(region, coef, dIxy, J);

				for(auto iter = 0; iter < n_iter; iter++)
				{
					if (bb_A_b)
					{
						// calculate matrix A and b
						cal_A_b(nxy, J, dIxy, A, b);

						// get the maximum gradient
						G_max = maximum_gradient(b);
					}

					// add lamba x diagonal to matrix A_lambda
					add_lambda(A, lambda, A_lambda);

					// solve the system of equations for d_coef
					mns_svd(A_lambda.data(), b.data(), d_coef.data());

					// add d_coef
					add_d_coef(coef, d_coef, coef_t);

					// set order coef
					set_order_coef(coef_t);

					// set constraints
					set_constraints(coef, coef_min, coef_max, coef_t);

					// calculate maximum error coefficients
					ee_coef_max = diff_coef(coef, coef_t);

					// calculate error
					ee_t = cal_ee_dIxy_J(region, coef_t, dIxy, J);

					if ((G_max < ee_G_max)||(fabs(ee-ee_t) < ee_ee)||(ee_coef_max < ee_coef))
					{
						break;
					}

					if (ee > ee_t)
					{
						coef = coef_t;
						ee = ee_t;
						lambda = (bb_shrink_lambda)?::fmax(lambda/lambda_du, 1e-8):lambda;
						bb_A_b = true;
						bb_shrink_lambda = true;
					}
					else
					{
						lambda = ::fmin(lambda*lambda_fu, 1e+8);
 						bb_A_b = false;
						bb_shrink_lambda = false;
					}
				}

				// backward transformation coefficients
				bk_transformation_coef(region, coef);

				return ee;
			}

			/***************************************************************************************/
			// set random seed
			void set_random_seed()
			{
				std::random_device rd;
				gen_c.seed(rd());
			}

			// allocate memory
			void set_size_local_variables(dt_int32 nxy, dt_int32 n_coef)
			{
				Ixy.resize(nxy);
				dIxy.resize(nxy);
				J.resize(nxy*n_coef);
				A.resize(n_coef*n_coef);
				A_lambda.resize(n_coef*n_coef);
				b.resize(n_coef);

				coef_min.resize(n_coef);
				coef_max.resize(n_coef);
				coef.resize(n_coef);
				d_coef.resize(n_coef);
				coef_t.resize(n_coef);
			}

			void init_variables(Grid_2d<T>& grid_2d_i, T lambda_i, 
			dt_int32 n_coef_i, dt_int32 n_coef_nl_i)
			{
				// set random seed
				set_random_seed();

				grid_2d = grid_2d_i;

				n_coef = n_coef_i;
				n_coef_nl = n_coef_nl_i;
				n_coef_l = n_coef - n_coef_nl;

				// set default values
				ee_ee = 1e-5;
				ee_G_max = 1e-5;
				ee_coef = 1e-7;

				lambda = lambda_i;
				lambda_fu = 8;
				lambda_du = 5;

				// allocate memory
				set_size_local_variables(grid_2d.size(), n_coef);

				// SVD minimum-norm solution
				mns_svd.init(n_coef, n_coef);
			}

			/***************************************************************************************/
			// add lamba x diagonal to matrix A_lambda
			void add_lambda(TVctr_r &A, T lambda, TVctr_r &B)
			{
				std::copy(A.begin(), A.end(), B.begin());
				for(auto ic=0; ic<n_coef; ic++)
				{
					const dt_int32 id = ic*n_coef+ic;
					B[id] += lambda*B[id];
				}
			}

			// get the maximum gradient
			T maximum_gradient(TVctr_r &y)
			{
				T G_max = 0;
				for(auto ic=0; ic<n_coef; ic++)
				{
					G_max = ::fmax(G_max, fabs(y[ic]));
				}

				return G_max;
			}

			// calculate maximum error coefficients
			T diff_coef(TVctr_r &coef, TVctr_r &coef_t)
			{
				T d_coef_sum = 0;
				for(dt_int32 ic = 0; ic < n_coef; ic++)
				{
					d_coef_sum += ::fabs(coef_t[ic]-coef[ic]);
				}

				return d_coef_sum;
			}

			// add d_coef
			void add_d_coef(TVctr_r &coef_i, TVctr_r &d_coef, TVctr_r &coef_o)
			{
				std::transform(coef_i.begin(), coef_i.end(), d_coef.begin(), coef_o.begin(), std::plus<T>());
			}

			// set constraints
			void set_constraints(TVctr_r &coef_i, TVctr_r &coef_min, 
			TVctr_r &coef_max, TVctr_r &coef_o)
			{
				dt_int32 is = 0;
				for(auto ic = 0; ic < n_coef; ic++)
				{
					if ((coef_min[ic]>coef_o[ic])||(coef_o[ic]>coef_max[ic]))
					{
						coef_o[ic] = coef_i[ic];
						is++;
					}
				}

				if (is==n_coef)
				{
					for(auto ic = 0; ic < n_coef; ic++)
					{
						T f = 0.9 + 0.0999*rand_c(gen_c);
						coef_o[ic] = f*coef_i[ic];
					}
				}
			}

			// calculate matrix A and b
			void cal_A_b(dt_int32 nxy, TVctr_r &J, TVctr_r &dIxy, TVctr_r &A, TVctr_r &b)
			{
				/***************************************************************************************/
				// get ATxA
				KS<T> ks;
				for(auto ic=0; ic<n_coef; ic++)
				{
					for(auto ir=0; ir<n_coef; ir++)
					{
						ks = 0;
						for(auto k=0;k<nxy;k++)
						{
							ks += J[ir*nxy+k]*J[ic*nxy+k];
						}

						A[ic*n_coef+ir] = ks;
					}
				}

				/***************************************************************************************/
				// get b
				for(auto ir=0; ir<n_coef; ir++)
				{
					ks = 0;
					for(auto k=0;k<nxy;k++)
					{
						ks += J[ir*nxy+k]*dIxy[k];
					}
					b[ir] = ks;
				}
			}

			void ellipse_var(T sx, T sy, T theta, T& a, T& b, T& c)
			{
				const T cos_1t = cos(theta);
				const T sin_1t = sin(theta);
				const T cos_2t = cos(2*theta);
				const T sin_2t = sin(2*theta);

				const T sx2 = pow(sx, 2);
				const T sy2 = pow(sy, 2);

				a = 0.5*(cos_1t*cos_1t/sx2 + sin_1t*sin_1t/sy2);
				b = 0.5*(sin_1t*sin_1t/sx2 + cos_1t*cos_1t/sy2);
				c = 0.5*(-sin_2t/sx2 + sin_2t/sy2);
			}

			void select_cir_reg(TVctr_r &Ixy, R_2d<T>& p, T radius, TRegion_r &region)
			{
				const T R_max = radius;
				const T R2_max = pow(R_max, 2);
				const T Rl2_max = pow(2.5*grid_2d.dR_min(), 2);

				auto range = grid_2d.region_ind(p, R_max);

				region.clear();
				region.reserve(range.ind_e);

				// select circular region
				KS<T> I_xy_m = 0;
				dt_int32 I_c = 0;
				for(auto ix = range.ix_0; ix < range.ix_e; ix++)
				{
					for(auto iy = range.iy_0; iy < range.iy_e; iy++)
					{
						const T R2 = grid_2d.r2(ix, iy, p.x, p.y);
						if (R2 < R2_max)
						{
							auto v = Ixy[grid_2d.sub_2_ind(ix, iy)];
							region.Rx.push_back(grid_2d.rx(ix));
							region.Ry.push_back(grid_2d.ry(iy));
							region.Ixy.push_back(v);
							if (R2 < Rl2_max)
							{
								I_xy_m += v;
								I_c++;
							}
						}

					}
				}
				region.shrink_to_fit();

				I_xy_m = I_xy_m/I_c;
				region.R_max = R_max;
				region.Rx_sf = p.x;
				region.Ry_sf = p.y;
				region.Rxy_sc = region.R_max;
				region.Ixy_sf = 0.0;
				region.Ixy_sc = 0.5*(I_xy_m + fcn_max_element(region.Ixy));

				// scaling and shifting data
				region.sf_sc();
			};

			void select_ellip_reg(TVctr_r &Ixy, TVctr_r &coef, T f0, TRegion_r &region)
			{
				const R_2d<T> p(coef[0], coef[1]);

				T a, b, c;
				ellipse_var(coef[3], coef[4], coef[5], a, b, c);

				const T d = log(f0);
				const T dd = c*c - 4*a*b;

				const T radius_y = sqrt(fabs(4*a*d/dd));
				const T radius_x = sqrt(fabs(4*b*d/dd));

				const T x_0 = p.x - radius_x;
				const T x_e = p.x + radius_x;

				const T y_0 = p.y - radius_y;
				const T y_e = p.y + radius_y;

				const dt_int32 ix_0 = grid_2d.rx_2_irx_cds(x_0);
				const dt_int32 ix_e = grid_2d.rx_2_irx_fds(x_e) + 1;

				const dt_int32 iy_0 = grid_2d.ry_2_iry_cds(y_0);
				const dt_int32 iy_e = grid_2d.ry_2_iry_fds(y_e) + 1;

				region.clear();
				region.reserve((ix_e - ix_0)*(iy_e - iy_0));

				// select elliptical region
				for(auto ix = ix_0; ix < ix_e; ix++)
				{
					T dx = grid_2d.rx(ix) - p.x;
					T ddy = sqrt(fabs(dd*dx*dx - 4*b*d));

					T y_0 = (-c*dx - ddy)/(2*b) + p.y;
					T y_e = (-c*dx + ddy)/(2*b) + p.y;

					dt_int32 iy_0 = grid_2d.ry_2_iry_bfds(y_0);
					dt_int32 iy_e = grid_2d.ry_2_iry_bcds(y_e) + 1;

					for(auto iy = iy_0; iy < iy_e; iy++)
					{
						region.Rx.push_back(grid_2d.rx(ix));
						region.Ry.push_back(grid_2d.ry(iy));
						region.Ixy.push_back(Ixy[grid_2d.sub_2_ind(ix, iy)]);
					}
				}
				region.shrink_to_fit();

				region.R_max = max(coef[3], coef[4]);
				region.Rx_sf = p.x;
				region.Ry_sf = p.y;
				region.Rxy_sc = region.R_max;
				region.Ixy_sf = 0.0;
				region.Ixy_sc = fcn_max_element(region.Ixy);

				// scaling and shifting data
				region.sf_sc();
			};
	};

	 /***************************************************************************************/
	 /***************************************************************************************/

	template <class T, eDev Dev>
	class Eval_Ellipt_Gauss_2d
	{
		public:
			using T_r = T;
			using TVctr_r = vector<T>;
			using TRegion_r = Region_Rad_2d<TVctr_r>;
			using size_type = dt_uint64;

			static const eDev device = Dev;

			Eval_Ellipt_Gauss_2d() {};

			void operator()(TVctr_r &coef, TVctr_r &Rx, TVctr_r &Ry, TVctr_r &Ixy)
			{
				const T x_0 = coef[0];
				const T y_0 = coef[1];
				const T A = coef[2];
				const T theta = coef[5];
				const T cos_1t = cos(theta);
				const T sin_1t = sin(theta);
				const T cos_2t = cos(2*theta);
				const T sin_2t = sin(2*theta);

				const T sx2 = pow(coef[3], 2);
				const T sy2 = pow(coef[4], 2);

				const T a = 0.5*(cos_1t*cos_1t/sx2 + sin_1t*sin_1t/sy2);
				const T b = 0.5*(sin_1t*sin_1t/sx2 + cos_1t*cos_1t/sy2);
				const T c = 0.5*(-sin_2t/sx2 + sin_2t/sy2);

				const dt_int32 nxy = Ixy.size();
				for(auto ixy = 0; ixy < nxy; ixy++)
				{
					T x = Rx[ixy] - x_0;
					T y = Ry[ixy] - y_0;

					T v = exp(-a*x*x - b*y*y - c*x*y);
					Ixy[ixy] = A*v;
				}
			}
	};

	template <class T, eDev Dev>
	class Fit_Ellipt_Gauss_2d: public Fit_one_rad_fcn_2d<T, Dev>{
		public:
			using T_r = T;
			using TVctr_r = vector<T>;
			using TRegion_r = Region_Rad_2d<TVctr_r>;
			using size_type = dt_uint64;

			static const eDev device = Dev;

			Fit_Ellipt_Gauss_2d():Fit_one_rad_fcn_2d<T, Dev>() {};

			void init_variables(Grid_2d<T>& grid_2d_i, T lambda_i)
			{
				Fit_one_rad_fcn_2d<T, Dev>::init_variables(grid_2d_i, lambda_i, 6, 4);
			}

			template <class TVctr>
			R_2d<T> fd_max_peak_pos(TVctr& Ixy_i)
			{
				dt_int32 ixy_max = thrust::max_element(Ixy_i.begin(), Ixy_i.end())-Ixy_i.begin();
				dt_int32 ix_max, iy_max;
				grid_2d.ind_2_sub(ixy_max, ix_max, iy_max);

				return R_2d<T>(grid_2d.rx(ix_max), grid_2d.ry(iy_max));
			}

			template <class TVctr>
			R_2d<T> fit_peak_pos(TVctr& Ixy_i, R_2d<T> pos_i, T sigma_i, T radius_i)
			{
				auto coef = this->fit(Ixy_i, pos_i, sigma_i, radius_i);
				return R_2d<T>(coef[0], coef[1]);
			}

		protected:
			// set initial coefficients and limits
			void set_init_min_max_coef_0(R_2d<T> p, T Ip_max, T sigma, T radius, TVctr_r &coef, TVctr_r &coef_min, TVctr_r &coef_max) 
			{
				const T sigma_x = sigma;
				const T sigma_y = sigma;

				/***************************************************************************************/
				// 0.8 = exp(-0.5*f^2)
				const T dx = ::fmax(0.668*sigma_x, 1.5*grid_2d.drx);
				const T x_min = grid_2d.set_bound_x(p.x-dx);
				const T x_max = grid_2d.set_bound_x(p.x+dx);

				// 0.8 = exp(-0.5*f^2)
				const T dy = ::fmax(0.668*sigma_y, 1.5*grid_2d.dry);
				const T y_min = grid_2d.set_bound_y(p.y-dy);
				const T y_max = grid_2d.set_bound_y(p.y+dy);

				const T sigma_x_min = ::fmax(sigma_x/5, 0.5*grid_2d.drx);
				const T sigma_x_max = ::fmin(10*radius, grid_2d.bs_x_h());

				const T sigma_y_min = ::fmax(sigma_y/5, 0.5*grid_2d.dry);
				const T sigma_y_max = ::fmin(10*radius, grid_2d.bs_y_h());

				/***************************************************************************************/
				coef = { p.x, p.y, Ip_max, sigma, sigma, T(c_pi)};
				coef_min = { x_min, y_min, T(0.25*Ip_max), sigma_x_min, sigma_y_min, T(0) };
				coef_max = { x_max, y_max, T(1.5*Ip_max), sigma_x_max, sigma_y_max, T(c_2pi) };
			}

			// set coefficients and limits
			void set_min_max_coef(TVctr_r &coef, TVctr_r &coef_min, TVctr_r &coef_max)
			{
				T Ip_max = coef[2];

				R_2d<T> p(coef[0], coef[1]);

				const T sigma_x = coef[3];
				const T sigma_y = coef[4];

				/***************************************************************************************/
				const T dx = ::fmax(sigma_x/4.0, 1.5*grid_2d.drx);
				const T x_min = grid_2d.set_bound_x(p.x-dx);
				const T x_max = grid_2d.set_bound_x(p.x+dx);

				const T dy = ::fmax(sigma_y/4.0, 1.5*grid_2d.dry);
				const T y_min = grid_2d.set_bound_y(p.y-dy);
				const T y_max = grid_2d.set_bound_y(p.y+dy);

				const T sigma_x_min = ::fmax(sigma_x/5.0, 0.5*grid_2d.drx);
				const T sigma_x_max = ::fmin(3.0*sigma_x, grid_2d.bs_x_h());

				const T sigma_y_min = ::fmax(sigma_y/5.0, 0.5*grid_2d.dry);
				const T sigma_y_max = ::fmin(3.0*sigma_y, grid_2d.bs_y_h());

				/***************************************************************************************/
				coef_min = { x_min, y_min, T(0.5*Ip_max), sigma_x_min, sigma_y_min, T(0) };
				coef_max = { x_max, y_max, T(1.50*Ip_max), sigma_x_max, sigma_y_max, T(c_2pi) };
			}

			// forward transformation coefficients
			void fw_transformation_coef(TRegion_r &region, TVctr_r &coef) 
			{
				coef[0] = (coef[0] - region.Rx_sf)/region.Rxy_sc;
				coef[1] = (coef[1] - region.Ry_sf)/region.Rxy_sc;
				coef[2] = coef[2]/region.Ixy_sc;
				coef[3] = coef[3]/region.Rxy_sc;
				coef[4] = coef[4]/region.Rxy_sc;
			}

			// backward transformation coefficients
			void bk_transformation_coef(TRegion_r &region, TVctr_r &coef) 
			{
				coef[0] = coef[0]*region.Rxy_sc + region.Rx_sf;
				coef[1] = coef[1]*region.Rxy_sc + region.Ry_sf;
				coef[2] = coef[2]*region.Ixy_sc;
				coef[3] = coef[3]*region.Rxy_sc;
				coef[4] = coef[4]*region.Rxy_sc;
			}

			// set order
			void set_order_coef(TVctr_r &coef)
			{
				T theta = coef[5];
				theta = theta - ::floor(theta/c_2pi)*c_2pi;
				coef[5] = theta;
			}

			T cal_ee_dIxy_J(TRegion_r &region, TVctr_r &coef, TVctr_r &d_Ixy, TVctr_r &J)
			{
				const T x_0 = coef[0];
				const T y_0 = coef[1];
				const T A = coef[2];
				const T theta = coef[5];
				const T cos_1t = cos(theta);
				const T sin_1t = sin(theta);
				const T cos_2t = cos(2*theta);
				const T sin_2t = sin(2*theta);

				const T sx2 = pow(coef[3], 2);
				const T sy2 = pow(coef[4], 2);

				const T sx3 = pow(coef[3], 3);
				const T sy3 = pow(coef[4], 3);

				const T a = 0.5*(cos_1t*cos_1t/sx2 + sin_1t*sin_1t/sy2);
				const T b = 0.5*(sin_1t*sin_1t/sx2 + cos_1t*cos_1t/sy2);
				const T c = 0.5*(-sin_2t/sx2 + sin_2t/sy2);

				KS<T> ee = 0;

				const dt_int32 nxy = region.size();
				for(auto ixy = 0; ixy < nxy; ixy++)
				{
					T x = region.Rx[ixy] - x_0;
					T y = region.Ry[ixy] - y_0;

					T v = exp(-a*x*x - b*y*y - c*x*y);
					T Ixy_p = A*v;

					J[0*nxy + ixy] = (2*a*x + c*y)*Ixy_p;
					J[1*nxy + ixy] = (2*b*y + c*x)*Ixy_p;
					J[2*nxy + ixy] = v;
					J[3*nxy + ixy] = (pow(cos_1t*x, 2) + pow(sin_1t*y, 2) - sin_2t*x*y)*Ixy_p/sx3;
					J[4*nxy + ixy] = (pow(sin_1t*x, 2) + pow(cos_1t*y, 2) + sin_2t*x*y)*Ixy_p/sy3;
					J[5*nxy + ixy] = (sin_2t*x*x - sin_2t*y*y + 2*cos_2t*x*y)*Ixy_p*(0.5/sx2 - 0.5/sy2);

					T dIxy_p = region.Ixy[ixy] - Ixy_p;
					d_Ixy[ixy] = dIxy_p;
					ee += ::fabs(dIxy_p);
				}
				ee = 100*ee/region.Ixy_sum;

				return ee;
			}
	};

	 /***************************************************************************************/
	 /***************************************************************************************/
	template <class T, eDev Dev>
	class Eval_2d_PAC
	{
		public:
			using T_r = T;
			using TVctr_r = vector<T>;
			using TVctr_r2d = vector<R_2d<T>>;
			using size_type = dt_uint64;

			Eval_2d_PAC() {};

			virtual void operator()(TVctr_r &Rx_i, TVctr_r &Ry_i, TVctr_r &pos_i, TVctr_r &coef, TVctr_r &Ixy)=0;

		protected:
			struct PTypes
			{
				public:
					TVctr_r2d p;

					dt_int32 size() const
					{
						return p.size();
					}

					void reserve(size_type n_pos_typ)
					{
						p.reserve(n_pos_typ);
					}	

					void shrink_to_fit()
					{
						p.shrink_to_fit();
					}
			};

			using VPTypes = vector<PTypes>;

			VPTypes split_pos_by_types(TVctr_r &pos_i)
			{
				const auto n_pos_typ = pos_i.size()/3;

				// get unique types
				TVctr_r ptypes_u(pos_i.begin()+2*n_pos_typ, pos_i.end());
				std::sort(ptypes_u.begin(), ptypes_u.end());
				ptypes_u.erase(std::unique(ptypes_u.begin(), ptypes_u.end(), mt::fcn_is_equal<T>), ptypes_u.end());
				ptypes_u.shrink_to_fit();

				dt_int32 n_ptypes_u = ptypes_u.size();
				VPTypes pos_typ(n_ptypes_u);

				// create map
				std::map<dt_int32, dt_int32> mptypes_u;
				for(auto idx=0; idx<n_ptypes_u; idx++)
				{
					mptypes_u.emplace(ptypes_u[idx], idx);
					pos_typ[idx].reserve(n_pos_typ);
				}

				// assign types of positions
				for(auto ip=0; ip<n_pos_typ; ip++)
				{
					auto idx = static_cast<dt_int32>(pos_i[n_pos_typ*2+ip]);
					idx = mptypes_u[idx];
					T x = pos_i[n_pos_typ*0+ip];
					T y = pos_i[n_pos_typ*1+ip];
					pos_typ[idx].p.push_back(R_2d<T>(x, y));
				}

				// shrink to fit
				for(dt_int32 idx=0; idx<n_ptypes_u; idx++)
				{
					pos_typ[idx].shrink_to_fit();
				}

				return pos_typ;
			}
	};

	template <class T, eDev Dev>
	class Fit_2d_PAC
	{
		public:
			using T_r = T;
			using TVctr_r = vector<T>;
			using TVctr_r2d = vector<R_2d<T>>;
			using size_type = dt_uint64;

			static const eDev device = Dev;

			Fit_2d_PAC():n_pos_typ(1), n_coef_pa(2), n_coef_nl_pa(1), n_coef_l_pa(1), n_coef(0), n_coef_nl(0), n_cstr(0), 
			lambda(1), lambda_fu(8), lambda_du(5), r2_min(1e-10), Ixy_i_sc(1), Ixy_0_sum(1), 
			ee_ee(1e-5), ee_G_max(1e-5), ee_coef(1e-7), n_iter(150) {};

			T fit(TVctr_r &Ixy_i, TVctr_r &cstr_i, TVctr_r &coef)
			{
				dt_bool bb_A_b = true;
				dt_bool bb_shrink_lambda = false;

				T G_max;
				T ee_coef_max;
				T ee_t;

				// set Ixy and constraints
				set_Ixy_cstr(Ixy_i, cstr_i);

				// forward transformation coefficients
				fw_transformation_coef(coef);

				// set order coef
				set_order_coef(coef);

				// set initial coefficients and limits
				set_min_max_coef(coef_min, coef_max);

				// call ee
				T ee = cal_ee_Ixy_dIxy_J(pos_typ, Rxy_0, Ixy_0, coef, Ixy, dIxy, J);

				for(auto iter = 0; iter < n_iter; iter++)
				{
					if (bb_A_b)
					{
						// calculate matrix A and b
						cal_A_b(J, dIxy, A, b);

						// get the maximum gradient
						G_max = maximum_gradient(b);
					}

					// add lamba x diagonal to matrix A_lambda
					add_lambda(A, lambda, A_lambda);

					// solve the system of equations for d_coef
					mns_svd(A_lambda.data(), b.data(), d_coef.data());

					// add d_coef
					add_d_coef(coef, d_coef, coef_t);

					// set constraints
					set_constraints(coef, coef_min, coef_max, coef_t);

					// calculate maximum error coefficients
					ee_coef_max = diff_coef(coef, coef_t);

					// calculate error
					ee_t = cal_ee_Ixy_dIxy_J(pos_typ, Rxy_0, Ixy_0, coef_t, Ixy, dIxy, J);

					if ((G_max < ee_G_max)||(fabs(ee-ee_t) < ee_ee)||(ee_coef_max < ee_coef))
					{
						break;
					}

					if (ee > ee_t)
					{
						set_order_coef(coef_t);

						coef = coef_t;
						ee = ee_t;
						lambda = (bb_shrink_lambda)?::fmax(lambda/lambda_du, 1e-8):lambda;
						bb_A_b = true;
						bb_shrink_lambda = true;
					}
					else
					{
						lambda = ::fmin(lambda*lambda_fu, 1e+8);
 						bb_A_b = false;
						bb_shrink_lambda = false;
					}
				}

				// backward transformation coefficients
				bk_transformation_coef(coef);

				return ee;
			}

		protected:
			struct PTypes
			{
				public:
					TVctr_r2d p;

					dt_int32 size() const
					{
						return p.size();
					}

					void reserve(size_type n_pos)
					{
						p.reserve(n_pos);
					}	

					void shrink_to_fit()
					{
						p.shrink_to_fit();
					}
			};

			using VPTypes = vector<PTypes>;

			T lambda;
			T lambda_fu;
			T lambda_du;

			dt_int32 n_iter;
			T ee_ee;
			T ee_G_max;
			T ee_coef;
			T r2_min;

			T Ixy_i_sc;
			T Ixy_0_sum;

			dt_int32 n_pos_typ;
			dt_int32 n_coef_pa;
			dt_int32 n_coef_nl_pa;
			dt_int32 n_coef_l_pa;
			dt_int32 n_coef;
			dt_int32 n_coef_nl;
			dt_int32 nxy;
			dt_int32 n_cstr;

			std::mt19937_64 gen_c;
			std::uniform_real_distribution<T> rand_c;

			VPTypes pos_typ;
			TVctr_r2d Rxy_0;
			TVctr_r Ixy_0;
			TVctr_r Ixy;
			TVctr_r dIxy;
			TVctr_r J;
			TVctr_r A;
			TVctr_r A_lambda;
			TVctr_r b;
			TVctr_r cstr;

			TVctr_r coef_min;
			TVctr_r coef_max;
			TVctr_r coef;
			TVctr_r d_coef;
			TVctr_r coef_t;

			lapack::MNS_SVD<T> mns_svd;

			/***************************************************************************************/
			// set initial coefficients and limits
			virtual void set_min_max_coef(TVctr_r &coef_min, TVctr_r &coef_max) = 0;

			// set order
			virtual void set_order_coef(TVctr_r &coef) = 0;

			// calculate error, jacobian, and dIxy
			virtual T cal_ee_Ixy_dIxy_J(VPTypes &pos_typ, TVctr_r2d &Rxy_0, TVctr_r &Ixy_0, 
			TVctr_r &coef, TVctr_r &Ixy, TVctr_r &dIxy, TVctr_r &J) = 0;

			/***************************************************************************************/
			// set Ixy and constraints
			void set_Ixy_cstr(TVctr_r &Ixy_i, TVctr_r &cstr_i)
			{
				// calculate intensity scaling factor
				Ixy_i_sc = fcn_max_element(Ixy_i);

				// assign data
				KS<T> sum = 0;
				for(auto ixy = 0; ixy < nxy; ixy++)
				{
					T Ixy = Ixy_i[ixy]/Ixy_i_sc;
					sum += Ixy;

					Ixy_0[ixy] = Ixy;
				}
				Ixy_0_sum = sum;

				cstr = cstr_i;
				n_cstr = cstr_i.size();
				for(auto ik = 0; ik < n_cstr; ik++)
				{
					cstr[ik] /= Ixy_i_sc;
				}
			}

			// set random seed
			void set_random_seed()
			{
				std::random_device rd;
				gen_c.seed(rd());
			}

			// allocate memory
			void set_size_local_variables(dt_int32 nxy, dt_int32 n_coef)
			{
				Rxy_0.resize(nxy);
				Ixy_0.resize(nxy);
				Ixy.resize(nxy);
				dIxy.resize(nxy);
				J.resize(nxy*n_coef);
				A.resize(n_coef*n_coef);
				A_lambda.resize(n_coef*n_coef);
				b.resize(n_coef);

				coef_min.resize(n_coef);
				coef_max.resize(n_coef);
				coef.resize(n_coef);
				d_coef.resize(n_coef);
				coef_t.resize(n_coef);
			}

			// set Rxy_0
			void set_rv(TVctr_r &Rx_i, TVctr_r &Ry_i)
			{
				for(auto ixy = 0; ixy < nxy; ixy++)
				{
					Rxy_0[ixy] = R_2d<T>(Rx_i[ixy], Ry_i[ixy]);
				}
			}

			VPTypes split_pos_by_types(TVctr_r &pos_i)
			{
				const auto n_pos = pos_i.size()/3;

				// get unique types
				TVctr_r ptypes_u(pos_i.begin()+2*n_pos, pos_i.end());
				std::sort(ptypes_u.begin(), ptypes_u.end());
				ptypes_u.erase(std::unique(ptypes_u.begin(), ptypes_u.end(), mt::fcn_is_equal<T>), ptypes_u.end());
				ptypes_u.shrink_to_fit();

				dt_int32 n_ptypes_u = ptypes_u.size();
				VPTypes pos_typ(n_ptypes_u);

				// create map
				std::map<dt_int32, dt_int32> mptypes_u;
				for(auto idx=0; idx<n_ptypes_u; idx++)
				{
					mptypes_u.emplace(ptypes_u[idx], idx);
					pos_typ[idx].reserve(n_pos);
				}

				// assign types of positions
				for(auto ip=0; ip<n_pos; ip++)
				{
					auto idx = static_cast<dt_int32>(pos_i[n_pos*2+ip]);
					idx = mptypes_u[idx];
					T x = pos_i[n_pos*0+ip];
					T y = pos_i[n_pos*1+ip];
					pos_typ[idx].p.push_back(R_2d<T>(x, y));
				}

				// shrink to fit
				for(dt_int32 idx=0; idx<n_ptypes_u; idx++)
				{
					pos_typ[idx].shrink_to_fit();
				}

				return pos_typ;
			}

			void init_variables(TVctr_r &Rx_i, TVctr_r &Ry_i, TVctr_r &pos_i, T lambda_i, 
			dt_int32 n_coef_pa_i, dt_int32 n_coef_nl_pa_i)
			{
				// set random seed
				set_random_seed();

				pos_typ = split_pos_by_types(pos_i);
				n_pos_typ = pos_typ.size();

				n_coef_pa = n_coef_pa_i;
				n_coef_nl_pa = n_coef_nl_pa_i;
				n_coef_l_pa = n_coef_pa - n_coef_nl_pa;

				n_coef = n_coef_pa*n_pos_typ;
				n_coef_nl = n_coef_nl_pa*n_pos_typ;
				nxy = Rx_i.size();

				// set default values
				ee_ee = 1e-5;
				ee_G_max = 1e-5;
				ee_coef = 1e-7;
				n_iter = 150;

				lambda = lambda_i;
				lambda_fu = 8;
				lambda_du = 5;

				// allocate memory
				set_size_local_variables(nxy, n_coef);

				// set Rxy_0
				set_rv(Rx_i, Ry_i);

				// SVD minimum-norm solution
				mns_svd.init(n_coef, n_coef);
			}

			/***************************************************************************************/
			// add lamba x diagonal to matrix A_lambda
			void add_lambda(TVctr_r &A, T lambda, TVctr_r &B)
			{
				std::copy(A.begin(), A.end(), B.begin());
				for(auto ic=0; ic<n_coef; ic++)
				{
					const dt_int32 id = ic*n_coef+ic;
					B[id] += lambda*B[id];
				}
			}

			// get the maximum gradient
			T maximum_gradient(TVctr_r &y)
			{
				T G_max = 0;
				for(auto ic=0; ic<n_coef; ic++)
				{
					G_max = ::fmax(G_max, fabs(y[ic]));
				}

				return G_max;
			}

			// calculate maximum error coefficients
			T diff_coef(TVctr_r &coef, TVctr_r &coef_t)
			{
				T d_coef_sum = 0;
				for(dt_int32 ic = 0; ic < n_coef; ic++)
				{
					d_coef_sum += ::fabs(coef_t[ic]-coef[ic]);
				}

				return d_coef_sum;
			}

			// add d_coef
			void add_d_coef(TVctr_r &coef_i, TVctr_r &d_coef, TVctr_r &coef_o)
			{
				std::transform(coef_i.begin(), coef_i.end(), d_coef.begin(), coef_o.begin(), std::plus<T>());
			}

			// check lambda
			dt_bool increase_lambda(TVctr_r &coef_min, TVctr_r &coef_max, TVctr_r &coef)
			{
				dt_bool bb_lambda = false;
				for(auto ic = 0; ic < n_coef; ic++)
				{
					if ((coef_min[ic]>coef[ic])||(coef[ic]>coef_max[ic]))
					{
						bb_lambda = true;
						break;
					}
				}

				return bb_lambda;
			}

			// set constraints
			void set_constraints(TVctr_r &coef_i, TVctr_r &coef_min, 
			TVctr_r &coef_max, TVctr_r &coef_o)
			{
				dt_int32 is = 0;
				for(auto ic = 0; ic < n_coef; ic++)
				{
					if ((coef_min[ic]>coef_o[ic])||(coef_o[ic]>coef_max[ic]))
					{
						coef_o[ic] = coef_i[ic];
						is++;
					}
				}

				if (is==n_coef)
				{
					for(auto ic = 0; ic < n_coef; ic++)
					{
						T f = 0.9 + 0.0999*rand_c(gen_c);
						coef_o[ic] = f*coef_i[ic];
					}
				}
			}

			// from mixed tononlinear -linear coeff
			void mixed_2_nonlinear_linear(TVctr_r &coef)
			{
				auto coef_t = coef;

				for(auto ipt = 0; ipt < n_pos_typ; ipt++)
				{
					const dt_int32 idx = n_coef_pa*ipt;
					const dt_int32 idx_nl = n_coef_nl_pa*ipt;
					const dt_int32 idx_l = n_coef_nl + idx_nl;

					for(auto ic = 0; ic < n_coef_pa; ic++)
					{
						auto idx_s = (ic<n_coef_l_pa)?(idx_l+ic):(idx_nl+ic-n_coef_l_pa);
						coef[idx_s]= coef_t[idx+ic];
					}
				}
			}

			// from nonlinear - linear coeff to mixed
			void nonlinear_linear_2_mixed(TVctr_r &coef)
			{
				auto coef_t = coef;

				for(auto ipt = 0; ipt < n_pos_typ; ipt++)
				{
					const dt_int32 idx = n_coef_pa*ipt;
					const dt_int32 idx_nl = n_coef_nl_pa*ipt;
					const dt_int32 idx_l = n_coef_nl + idx_nl;

					for(auto ic = 0; ic < n_coef_pa; ic++)
					{
						auto idx_s = (ic<n_coef_l_pa)?(idx_l+ic):(idx_nl+ic-n_coef_l_pa);
						coef[idx+ic]= coef_t[idx_s];
					}
				}
			}

			// forward transformation coefficients
			void fw_transformation_coef(TVctr_r &coef)
			{
				mixed_2_nonlinear_linear(coef);

				for(auto ic = n_coef_nl; ic < n_coef; ic++)
				{
					coef[ic] = log(coef[ic]/Ixy_i_sc);
				}
			}

			// backward transformation coefficients
			void bk_transformation_coef(TVctr_r &coef)
			{
				for(auto ic = n_coef_nl; ic < n_coef; ic++)
				{
					coef[ic] = Ixy_i_sc*exp(coef[ic]);
				}

				nonlinear_linear_2_mixed(coef);
			}

			// calculate matrix A and b
			void cal_A_b(TVctr_r &J, TVctr_r &dIxy, TVctr_r &A, TVctr_r &b)
			{
				const dt_int32 nxy = dIxy.size();
				const dt_int32 n_coef = b.size();

				/***************************************************************************************/
				// get ATxA
				KS<T> ks;
				for(auto ic=0; ic<n_coef; ic++)
				{
					for(auto ir=0; ir<n_coef; ir++)
					{
						ks = 0;
						for(auto k=0;k<nxy;k++)
						{
							ks += J[ir*nxy+k]*J[ic*nxy+k];
						}

						A[ic*n_coef+ir] = ks;
					}
				}

				/***************************************************************************************/
				// get b
				for(auto ir=0; ir<n_coef; ir++)
				{
					ks = 0;
					for(auto k=0;k<nxy;k++)
					{
						ks += J[ir*nxy+k]*dIxy[k];
					}
					b[ir] = ks;
				}
			}

			// calculate dIxy and ee
			T cal_dIxy_ee(TVctr_r &Ixy_0, TVctr_r &Ixy, TVctr_r &dIxy)
			{
				KS<T> ee = 0;
				for(auto ixy = 0; ixy < nxy; ixy++)
				{
					const T dIxy_p = Ixy_0[ixy] - Ixy[ixy];
					dIxy[ixy] = dIxy_p;
					ee += ::fabs(dIxy_p);
				}

				ee = 100*ee/Ixy_0_sum;

				return ee;
			}

			// set vector to zero
			void fill_vector(TVctr_r &Ixy, T val)
			{
				std::fill(Ixy.begin(), Ixy.end(), val);
			}
	};

	 /***************************************************************************************/
	template <class T, eDev Dev>
	class Eval_2_Gauss_PAC: public Eval_2d_PAC<T, Dev>{
		public:
			using T_r = T;
			using TVctr_r = vector<T>;
			using TVctr_r2d = vector<R_2d<T>>;
			using size_type = dt_uint64;

			Eval_2_Gauss_PAC(): Eval_2d_PAC<T, Dev>() {};

			void operator()(TVctr_r &Rx_i, TVctr_r &Ry_i, TVctr_r &pos_i, TVctr_r &coef, TVctr_r &Ixy)
			{
				auto pos_typ = split_pos_by_types(pos_i);

				const dt_int32 n_pos_typ = pos_typ.size();
				const dt_int32 n_coef_pa = coef.size()/n_pos_typ;
				const dt_int32 nxy = Rx_i.size();

				KS<T> Ixy_p;

				for(auto ixy = 0; ixy < nxy; ixy++)
				{
					Ixy_p = 0;

					const R_2d<T> Rxy = R_2d<T>(Rx_i[ixy], Ry_i[ixy]);
					for(auto ipt = 0; ipt < n_pos_typ; ipt++)
					{
						const dt_int32 idx_c = n_coef_pa*ipt;

						const T a1 = coef[idx_c + 0];
						const T a2 = coef[idx_c + 1];

						const T sigma1 = coef[idx_c + 2];
						const T f_sigma1 = 1.0/(2*sigma1*sigma1);

						const T sigma2 = coef[idx_c + 3];
						const T f_sigma2 = 1.0/(2*sigma2*sigma2);

						const auto &p = pos_typ[ipt].p;
						const dt_int32 np = p.size();

						// sum over positions
						for(auto ip = 0; ip < np; ip++)
						{
							const T r2 = norm_2(Rxy-p[ip]);
							const T ep_1 = a1*exp(-f_sigma1*r2);
							const T ep_2 = a2*exp(-f_sigma2*r2);

							Ixy_p += ep_1 + ep_2;
						}
					}

					Ixy[ixy] = Ixy_p;
				}
			}
	};

	// fit two gaussians per atomic column
	template <class T, eDev Dev>
	class Fit_2_Gauss_PAC: public Fit_2d_PAC<T, Dev>{
		public:
			using T_r = T;
			using TVctr_r = vector<T>;
			using TVctr_r2d = vector<R_2d<T>>;
			using size_type = dt_uint64;

			static const eDev device = Dev;

			Fit_2_Gauss_PAC(): Fit_2d_PAC<T, Dev>() {};

			void init_variables(TVctr_r &Rx_i, TVctr_r &Ry_i, TVctr_r &pos_i, T lambda_i)
			{
				Fit_2d_PAC<T, Dev>::init_variables(Rx_i, Ry_i, pos_i, lambda_i, 4, 2);
			}

		protected:
			// set initial coefficients and limits
			void set_min_max_coef(TVctr_r &coef_min, TVctr_r &coef_max)
			{
				for(auto ipt = 0; ipt < n_pos_typ; ipt++)
				{
					const dt_int32 idx = n_coef_pa*ipt;

					coef_min[idx + 0] = 1e-2*Ixy_i_sc;
					coef_min[idx + 1] = 1e-2*Ixy_i_sc;
					coef_min[idx + 2] = 0.10;
					coef_min[idx + 3] = 0.10;

					coef_max[idx + 0] = 1.10*Ixy_i_sc;
					coef_max[idx + 1] = 1.10*Ixy_i_sc;
					coef_max[idx + 2] = 5.0;
					coef_max[idx + 3] = 5.0;
				}

				this->fw_transformation_coef(coef_min);
				this->fw_transformation_coef(coef_max);
			}

			void set_order_coef(TVctr_r &coef)
			{

			}

			// calculate error, jacobian, and dIxy
			T cal_ee_Ixy_dIxy_J(VPTypes &pos_typ, TVctr_r2d &Rxy_0, TVctr_r &Ixy_0, 
			TVctr_r &coef, TVctr_r &Ixy, TVctr_r &dIxy, TVctr_r &J)
			{
				KS<T> Ixy_p;

				KS<T> J0;
				KS<T> J1;

				KS<T> J2;
				KS<T> J3;

				// set zero Ixy
				this->fill_vector(Ixy, 0);

				for(auto ipt = 0; ipt < n_pos_typ; ipt++)
				{
					const dt_int32 idx_nl = n_coef_nl_pa*ipt;
					const dt_int32 idx_l = n_coef_nl + idx_nl;

					const T sigma1 = coef[idx_nl + 0];
					const T sigma2 = coef[idx_nl + 1];

					const T a1 = coef[idx_l + 0];
					const T a2 = coef[idx_l + 1];

					const T sigma1_2 = sigma1*sigma1;
					const T sigma1_3 = sigma1_2*sigma1;
					const T f_sigma1 = 1.0/(2*sigma1_2);

					const T sigma2_2 = sigma2*sigma2;
					const T sigma2_3 = sigma2_2*sigma2;
					const T f_sigma2 = 1.0/(2*sigma2_2);

					const auto &p = pos_typ[ipt].p;
					const dt_int32 np = p.size();

					for(auto ixy = 0; ixy < nxy; ixy++)
					{
						Ixy_p = 0;
						J0 = J1 = J2 = J3 = 0;

						const auto Rxy = Rxy_0[ixy];
						for(auto ip = 0; ip < np; ip++)
						{
							const T r2 = norm_2(Rxy-p[ip]);
							const T ep_1 = exp(a1-f_sigma1*r2);
							const T ep_2 = exp(a2-f_sigma2*r2);

							Ixy_p += ep_1 + ep_2;

							J0 += ep_1*r2;
							J1 += ep_2*r2;
							J2 += ep_1;
							J3 += ep_2;
						}

						// assign values
						Ixy[ixy] += Ixy_p;

						J[(idx_nl + 0)*nxy + ixy] = J0/sigma1_3;
						J[(idx_nl + 1)*nxy + ixy] = J1/sigma2_3;
						J[(idx_l + 0)*nxy + ixy] = J2;
						J[(idx_l + 1)*nxy + ixy] = J3;
					}
				}

				T ee = this->cal_dIxy_ee(Ixy_0, Ixy, dIxy);

				return ee;
			}
	};

	 /***************************************************************************************/
	template <class T, eDev Dev>
	class Eval_Gauss_1P_PAC: public Eval_2d_PAC<T, Dev>{
		public:
			using T_r = T;
			using TVctr_r = vector<T>;
			using TVctr_r2d = vector<R_2d<T>>;
			using size_type = dt_uint64;

			Eval_Gauss_1P_PAC(): Eval_2d_PAC<T, Dev>() {};

			void operator()(TVctr_r &Rx_i, TVctr_r &Ry_i, TVctr_r &pos_i, TVctr_r &coef, TVctr_r &Ixy)
			{
				auto pos_typ = split_pos_by_types(pos_i);

				const dt_int32 n_pos_typ = pos_typ.size();
				const dt_int32 n_coef_pa = coef.size()/n_pos_typ;
				const dt_int32 nxy = Rx_i.size();

				KS<T> Ixy_p;

				for(auto ixy = 0; ixy < nxy; ixy++)
				{
					Ixy_p = 0;
					const R_2d<T> Rxy = R_2d<T>(Rx_i[ixy], Ry_i[ixy]);

					for(auto ipt = 0; ipt < n_pos_typ; ipt++)
					{
						const dt_int32 idx_c = n_coef_pa*ipt;

						const T a1 = coef[idx_c + 0];
						const T a2 = coef[idx_c + 1];
						const T alpha = coef[idx_c + 2];
						const T sigma = coef[idx_c + 3];
						const T f_sigma = 1.0/(2*sigma*sigma);

						const auto &p = pos_typ[ipt].p;
						const dt_int32 np = p.size();

						// sum over positions
						for(auto ip = 0; ip < np; ip++)
						{
							const T r2 = norm_2(Rxy-p[ip]);
							const T ep_p = pow(r2, alpha);
							const T ep_g = exp(-f_sigma*r2);

							Ixy_p += (a1 + a2*ep_p)*ep_g;
						}
					}

					Ixy[ixy] = Ixy_p;
				}
			}
	};

	template <class T, eDev Dev>
	class Fit_Gauss_1P_PAC: public Fit_2d_PAC<T, Dev>{
		public:
			using T_r = T;
			using TVctr_r = vector<T>;
			using TVctr_r2d = vector<R_2d<T>>;
			using size_type = dt_uint64;

			static const eDev device = Dev;

			Fit_Gauss_1P_PAC(): Fit_2d_PAC<T, Dev>() {};

			void init_variables(TVctr_r &Rx_i, TVctr_r &Ry_i, TVctr_r &pos_i, T lambda_i)
			{
				Fit_2d_PAC<T, Dev>::init_variables(Rx_i, Ry_i, pos_i, lambda_i, 4, 2);
			}

		protected:
			void set_min_max_coef(TVctr_r &coef_min, TVctr_r &coef_max)
			{
				for(auto ipt = 0; ipt < n_pos_typ; ipt++)
				{
					const dt_int32 idx = n_coef_pa*ipt;

					coef_min[idx + 0] = 1e-2*Ixy_i_sc;
					coef_min[idx + 1] = 1e-3*Ixy_i_sc;
					coef_min[idx + 2] = 0.01;
					coef_min[idx + 3] = 0.10;

					coef_max[idx + 0] = 1.10*Ixy_i_sc;
					coef_max[idx + 1] = 1.0*Ixy_i_sc;
					coef_max[idx + 2] = 8.0;
					coef_max[idx + 3] = 5.0;
				}

				this->fw_transformation_coef(coef_min);
				this->fw_transformation_coef(coef_max);
			}

			void set_order_coef(TVctr_r &coef)
			{

			}

			// calculate error, jacobian, and dIxy
			T cal_ee_Ixy_dIxy_J(VPTypes &pos_typ, TVctr_r2d &Rxy_0, TVctr_r &Ixy_0, 
			TVctr_r &coef, TVctr_r &Ixy, TVctr_r &dIxy, TVctr_r &J)
			{
				KS<T> Ixy_p;

				KS<T> J0;
				KS<T> J1;

				KS<T> J2;
				KS<T> J3;

				// set zero Ixy
				this->fill_vector(Ixy, 0);

				for(auto ipt = 0; ipt < n_pos_typ; ipt++)
				{
					const dt_int32 idx_nl = n_coef_nl_pa*ipt;
					const dt_int32 idx_l = n_coef_nl + idx_nl;

					const T alpha = coef[idx_nl + 0];
					const T sigma1 = coef[idx_nl + 1];

					const T a1 = coef[idx_l + 0];
					const T a2 = coef[idx_l + 1];

					const T sigma1_2 = sigma1*sigma1;
					const T sigma1_3 = sigma1_2*sigma1;
					const T f_sigma1 = 1.0/(2*sigma1_2);

					const auto &p = pos_typ[ipt].p;
					const dt_int32 np = p.size();

					for(auto ixy = 0; ixy < nxy; ixy++)
					{
						Ixy_p = 0;
						J0 = J1 = J2 = J3 = 0;

						const auto Rxy = Rxy_0[ixy];
						for(auto ip = 0; ip < np; ip++)
						{
							const T r2 = norm_2(Rxy-p[ip]);
							const T ln_r2 = (r2>r2_min)?log(r2):0;

							const T ep_1 = exp(a1-f_sigma1*r2);
							const T ep_2 = (r2>r2_min)?exp(a2+ln_r2*alpha-f_sigma1*r2):0;

							Ixy_p += ep_1 + ep_2;

							J0 += ln_r2*ep_2;
							J1 += ep_1*r2;
							J2 += ep_1;
							J3 += ep_2;
						}

						// assign values
						Ixy[ixy] += Ixy_p;

						J[(idx_nl + 0)*nxy + ixy] = J0;
						J[(idx_nl + 1)*nxy + ixy] = J1/sigma1_3;
						J[(idx_l + 0)*nxy + ixy] = J2;
						J[(idx_l + 1)*nxy + ixy] = J3;
					}
				}

				T ee = this->cal_dIxy_ee(Ixy_0, Ixy, dIxy);

				return ee;
			}
	};

}
#endif