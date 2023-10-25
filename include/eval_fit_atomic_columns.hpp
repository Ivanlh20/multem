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

#ifndef EVAL_FIT_ATOMIC_COLUMN_H
#define EVAL_FIT_ATOMIC_COLUMN_H

#include <random>
#include <vector>
#include <deque>
#include <map>

#include <fftw3.h>
#include "math_mt.h"
#include "type_traits_gen.h"
#include "types.cuh"
#include "lapack.hpp"
 namespace mt
{
	template <class T, eDev Dev>
	class Eval_Atomic_Columns
	{
		public:
			using T_r = T;
			using TVctr_r = vector<T>;
			using TVctr_r2d = vector<R_2d<T>>;
			using size_type = dt_uint64;

			Eval_Atomic_Columns() {};

			void operator()(TVctr_r &Rx_i, TVctr_r &Ry_i, TVctr_r &pos_i, TVctr_r &coef, TVctr_r &Ixy)
			{
				auto pos = split_pos_by_types(pos_i);

				const dt_int32 n_pos = pos.size();
				const dt_int32 n_coef_pa = coef.size()/n_pos;
				const dt_int32 nxy = Rx_i.size();

				KS<T> Ixy_p;

				for(auto ixy = 0; ixy < nxy; ixy++)
				{
					Ixy_p = 0;

					for(auto ipt = 0; ipt < n_pos; ipt++)
					{
						const dt_int32 idx_c = n_coef_pa*ipt;

						const T a1 = coef[idx_c + 0];
						const T a2 = coef[idx_c + 1];
						const T alpha = coef[idx_c + 2];
						const T sigma1 = coef[idx_c + 3];
						const T sigma2 = coef[idx_c + 4];

						const T f_sigma1 = 1.0/(2*sigma1*sigma1);
						const T f_sigma2 = 1.0/(2*sigma2*sigma2);

						const auto &p = pos[ipt].p;
						const dt_int32 np = p.size();

						const R_2d<T> Rxy = R_2d<T>(Rx_i[ixy], Ry_i[ixy]);

						// sum over positions
						for(auto ip = 0; ip < np; ip++)
						{
							const T r2 = norm_2(Rxy-p[ip]);
							const T ep_1 = a1*exp(-f_sigma1*r2);
							const T ep_2 = a2*pow(r2, alpha)*exp(-f_sigma2*r2);

							Ixy_p += ep_1 + ep_2;
						}
					}

					Ixy[ixy] = Ixy_p;
				}
			}

		private:
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
	};

	template <class T, eDev Dev>
	class Fit_Atomic_Columns
	{
		public:
			using T_r = T;
			using TVctr_r = vector<T>;
			using TVctr_r2d = vector<R_2d<T>>;
			using size_type = dt_uint64;

			static const eDev device = Dev;

			Fit_Atomic_Columns():n_pos_typ(0), n_coef_pa(5), n_coef_nl_pa(3), n_coef(0), n_coef_nl(0), n_cstr(0), r2_min(1e-10) {};

			T fit(TVctr_r &Ixy_i, TVctr_r &Rx_i, TVctr_r &Ry_i, TVctr_r &pos_i, T lambda_0, TVctr_r &coef, TVctr_r &cstr_i)
			{
				// set random seed
				std::random_device rd;
				gen_c.seed(rd());

				auto pos_typ = split_pos_by_types(pos_i);

				n_pos_typ = pos_typ.size();
				n_coef = n_coef_pa*n_pos_typ;
				n_coef_nl = n_coef_nl_pa*n_pos_typ;
				nxy = Ixy_i.size();
				n_cstr = cstr_i.size();

				// allocate memory
				set_size_local_variables(nxy, n_coef);

				// assign in data
				assign_in_data(Ixy_i, Rx_i, Ry_i, cstr_i);

				TVctr_r coef_min(n_coef);
				TVctr_r coef_max(n_coef);

				// set initial coefficients and limits
				set_coef_min_max(coef_min, coef_max);

				// SVD minimum-norm solution
				mns_svd.init(n_coef, n_coef);

				const T ee_ee = 1e-5;
				const T ee_G_max = 1e-5;
				const T ee_coef = 1e-7;
				const dt_int32 n_iter = 150;
				dt_bool bb_A_b = true;
				dt_bool bb_shrink_lambda = false;

				T G_max;
				T ee_coef_max;
				T ee_t;
				TVctr_r d_coef(n_coef);
				TVctr_r coef_t(n_coef);

				T lambda = lambda_0;
				const T lambda_fu = 8;
				const T lambda_du = 5;

				// from mixed to nonlinear -linear coeff
				mixed_2_nonlinear_linear(coef);

				// forward transformation coefficients
				fw_transformation_coef(coef);

				// set dependent coefficients
				// set_dependent_coef(cstr, coef);

				// call ee
				T ee_b = cal_ee_Ixy_dIxy_J(pos_typ, Rxy_0, Ixy_0, coef, Ixy, dIxy, J);

				for(auto iter = 0; iter < n_iter; iter++)
				{
					if (bb_A_b)
					{
						// calculate matrix A and b
						get_A_b(J, dIxy, A, b);

						// get the maximum gradient
						G_max = maximum_gradient(n_coef, b);
					}

					// add lamba x diagonal to matrix A_lambda
					add_lambda(n_coef, A, lambda, A_lambda);

					// solve the system of equations for d_coef
					mns_svd(A_lambda.data(), b.data(), d_coef.data());

					// add d_coef
					add_d_coef(coef, d_coef, coef_t);

					// increase lambda
					// if (increase_lambda(coef_min, coef_max, coef_t))
					// {
					// 	lambda = ::fmin(lambda*lambda_fu, 1e+8);
					// 	continue;
					// }

					// set constraints
					set_constraints(coef, coef_min, coef_max, coef_t);

					// set dependent coefficients
					// set_dependent_coef(cstr, coef_t);

					// calculate maximum error coefficients
					ee_coef_max = diff_coef(coef, coef_t);

					// calculate error
					ee_t = cal_ee_Ixy_dIxy_J(pos_typ, Rxy_0, Ixy_0, coef_t, Ixy, dIxy, J);

					if ((G_max < ee_G_max)||(fabs(ee_b-ee_t) < ee_ee)||(ee_coef_max < ee_coef))
					{
						break;
					}

					if (ee_b > ee_t)
					{
						coef = coef_t;
						ee_b = ee_t;
						lambda = (bb_shrink_lambda)?::fmax(lambda/lambda_du, 1e-10):lambda;
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

				// from nonlinear - linear coeff to mixed
				nonlinear_linear_2_mixed(coef);

				return ee_b;
			}

		private:
			const T r2_min;
			T Ixy_i_sc;
			T Ixy_0_sum;
			dt_int32 n_pos_typ;
			dt_int32 n_coef_pa;
			dt_int32 n_coef_nl_pa;
			dt_int32 n_coef;
			dt_int32 n_coef_nl;
			dt_int32 nxy;
			dt_int32 n_cstr;

			std::mt19937_64 gen_c;
			std::uniform_real_distribution<T> rnd_c;

			TVctr_r2d Rxy_0;
			TVctr_r Ixy_0;
			TVctr_r Ixy;
			TVctr_r dIxy;
			TVctr_r J;
			TVctr_r A;
			TVctr_r A_lambda;
			TVctr_r b;
			TVctr_r cstr;

			lapack::MNS_SVD<T> mns_svd;

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
			}

			// assign in data
			void assign_in_data(TVctr_r &Ixy_i, TVctr_r &Rx_i, TVctr_r &Ry_i, TVctr_r &cstr_i)
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
					Rxy_0[ixy] = R_2d<T>(Rx_i[ixy], Ry_i[ixy]);
				}
				Ixy_0_sum = sum;

				cstr = cstr_i;
				for(auto ik = 2; ik < n_cstr; ik++)
				{
					cstr[ik] /= Ixy_i_sc;
				}
			}

			// add lamba x diagonal to matrix A_lambda
			void add_lambda(dt_int32 n, TVctr_r &A, T lambda, TVctr_r &B)
			{
				std::copy(A.begin(), A.end(), B.begin());
				for(auto ic=0; ic<n; ic++)
				{
					B[ic*n+ic] += lambda*B[ic*n+ic];
				}
			}

			// get the maximum gradient
			T maximum_gradient(dt_int32 n, TVctr_r &y)
			{
				T G_max = 0;
				for(auto ic=0; ic<n; ic++)
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

			// set initial coefficients and limits
			void set_coef_min_max(TVctr_r &coef_min, TVctr_r &coef_max)
			{
				for(auto ipt = 0; ipt < n_pos_typ; ipt++)
				{
					const dt_int32 idx_nl = n_coef_nl_pa*ipt;
					const dt_int32 idx_l = n_coef_nl + idx_nl;

					coef_min[idx_nl + 0] = 0.01;
					coef_min[idx_nl + 1] = 0.10;
					coef_min[idx_nl + 2] = 0.10;

					coef_max[idx_nl + 0] = 10.0;
					coef_max[idx_nl + 1] = 6.0;
					coef_max[idx_nl + 2] = 6.0;

					coef_min[idx_l + 0] = log(0.01);
					coef_min[idx_l + 1] = log(0.000001);

					coef_max[idx_l + 0] = log(1.25);
					coef_max[idx_l + 1] = log(3.0);
				}
			}

			// set dependent coefficients
			void set_dependent_coef(TVctr_r &cstr, TVctr_r &coef)
			{
				// T ss = 0;
				// for(auto ic = 0; ic < n_coef_nl; ic++)
				// {
				// 	const T a = coef[n_coef_nl+ic];
				// 	const T sigma = coef[ic];
				// 	ss += a*sigma*sigma;
				// }

				// const T f = sqrt(cstr[2]/(c_2pi<T>*ss));

				// for(auto ic = 0; ic < n_coef_nl; ic++)
				// {
				// 	coef[ic] *= f;
				// }
			}

			// add d_coef
			void add_d_coef(TVctr_r &coef_i, TVctr_r &d_coef, TVctr_r &coef_o)
			{
				std::transform(coef_i.begin(), coef_i.end(), d_coef.begin(), coef_o.begin(), std::plus<T>());
			}

			// forward transformation coefficients
			void fw_transformation_coef(TVctr_r &coef)
			{
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
						T f = 0.9 + 0.0999*rnd_c(gen_c);
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

					coef[idx_l+0]= coef_t[idx+0];
					coef[idx_l+1] = coef_t[idx+1];

					coef[idx_nl+0] = coef_t[idx+2];
					coef[idx_nl+1] = coef_t[idx+3];
					coef[idx_nl+2] = coef_t[idx+4];
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

					coef[idx+0] = coef_t[idx_l+0];
					coef[idx+1] = coef_t[idx_l+1];

					coef[idx+2] = coef_t[idx_nl+0];
					coef[idx+3]= coef_t[idx_nl+1];
					coef[idx+4] = coef_t[idx_nl+2];
				}
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
				KS<T> J4;

				// set zero Ixy
				std::fill(Ixy.begin(), Ixy.end(), T(0));

				for(auto ipt = 0; ipt < n_pos_typ; ipt++)
				{
					const dt_int32 idx_nl = n_coef_nl_pa*ipt;
					const dt_int32 idx_l = n_coef_nl + idx_nl;

					const T alpha = coef[idx_nl + 0];
					const T sigma1 = coef[idx_nl + 1];
					const T sigma2 = coef[idx_nl + 2];

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
						J0 = J1 = J2 = J3 = J4 = 0;

						const auto Rxy = Rxy_0[ixy];
						// sum over positions
						for(auto ip = 0; ip < np; ip++)
						{
							const T r2 = norm_2(Rxy-p[ip]);
							const T ln_r2 = (r2>r2_min)?log(r2):0;

							const T ep_1 = exp(a1-f_sigma1*r2);
							const T ep_2 = (r2>r2_min)?exp(a2+ln_r2*alpha-f_sigma2*r2):0;

							Ixy_p += ep_1 + ep_2;

							J0 += ln_r2*ep_2;
							J1 += ep_1*r2;
							J2 += ep_2*r2;
							J3 += ep_1;
							J4 += ep_2;
						}

						// assign values
						Ixy[ixy] += Ixy_p;

						J[(idx_nl + 0)*nxy + ixy] = J0;
						J[(idx_nl + 1)*nxy + ixy] = J1/sigma1_3;
						J[(idx_nl + 2)*nxy + ixy] = J2/sigma2_3;
						J[(idx_l + 0)*nxy + ixy] = J3;
						J[(idx_l + 1)*nxy + ixy] = J4;
					}
				}

				/***************************************************************************************/
				KS<T> ee = 0;
				for(auto ixy = 0; ixy < nxy; ixy++)
				{
					T dIxy_p = Ixy_0[ixy] - Ixy[ixy];
					dIxy[ixy] = dIxy_p;
					ee += ::fabs(dIxy_p);
				}

				ee = 100*ee/Ixy_0_sum;

				return ee;
			}
			
			// calculate matrix A and b
			void get_A_b(TVctr_r &J, TVctr_r &dIxy, TVctr_r &A, TVctr_r &b)
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
	};

}
#endif