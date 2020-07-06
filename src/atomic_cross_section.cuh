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

#ifndef ATOMIC_CROSS_SECTION_H
#define ATOMIC_CROSS_SECTION_H

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "fft.cuh"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "tem_simulation.cuh"
#include "cubic_spline.hpp"
#include "cgpu_classes.cuh"
#include <thrust/transform.h>

namespace mt
{
	template <class T, eDevice dev>
	class Atomic_Cross_Section{
		public:
			using T_r = T;
			using T_c = complex<T>;

			Atomic_Cross_Section(): input_multislice(nullptr){}

			void set_input_data(Input_Multislice<T_r> *input_multislice_i)
			{
				input_multislice = input_multislice_i;
			}

			template <class TVector>
			void get(TVector &r_o, TVector &fr_o)
			{
				/**************************tem_simulation calculation*******************************/
				mt::Output_Multislice_Vector<T> output_multislice;
				output_multislice.set_input_data(input_multislice);

				mt::Multislice<T, dev> tem_simulation;
				tem_simulation.set_input_data(input_multislice);

				tem_simulation.run(output_multislice);

				tem_simulation.cleanup();

				/******************************* Cross section ********************************/
				mt::Grid_2d<double> grid_2d;
				grid_2d.assign(input_multislice->grid_2d);

				mt::Cubic_Spline<double> spline;
				auto sigma = input_multislice->atoms.sigma[0];

				auto pr = &(output_multislice.scanning.r);
				std::vector<double> r2(pr->size());
				thrust::transform(pr->begin(), pr->end(), r2.begin(), [](T &x)->double{return x*x; });

				auto pfr2 = &(output_multislice.image_tot[0].image[0]);
				std::vector<double> fr2(pr->size());
				fr2.assign(pfr2->begin(), pfr2->begin()+fr2.size());

				mt::Vector<complex<double>, e_host> M(grid_2d.nxy());
				T x = 0.5*grid_2d.lx;
				T y = 0.5*grid_2d.ly;

				// spline interpolation
				spline.set_points(r2, fr2);
				spline.eval_radial_function(grid_2d, x, y, M);

				// Convolution
				mt::Stream<e_host> stream(input_multislice->cpu_nthread);
				mt::FFT<double, e_host> fft_2d;

				mt::Gauss_Cv_2d<double, e_host> gauss_cv_2d(&stream, &fft_2d, grid_2d);
				gauss_cv_2d.set_fft_plan();
				gauss_cv_2d(sigma, M);
				gauss_cv_2d.cleanup();

				// extract radial values
				std::vector<double> r2_c(grid_2d.nyh);
				std::vector<double> fr_c(grid_2d.nyh);
				int ix = grid_2d.nxh;
				int ir = 0;
				for(auto iy =grid_2d.nyh; iy<grid_2d.ny; iy++)
				{
					r2_c[ir] = grid_2d.R2(ix, iy, x, y);
					fr_c[ir] = M[grid_2d.ind_col(ix, iy)].real();
					ir++;
				}

				// spline interpolation
				spline.set_points(r2_c, fr_c);
				spline.eval_function(r2, fr2);

				r_o.assign(pr->begin(), pr->end());
				fr_o.assign(fr2.begin(), fr2.end());
			}

		private:
			Input_Multislice<T_r> *input_multislice;
	};

} // namespace mt

#endif