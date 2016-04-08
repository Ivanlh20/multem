/*
 * This file is part of MULTEM.
 * Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>
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
#include "fft2.cuh"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "multislice.cuh"
#include "cubic_spline.hpp"
#include <thrust/transform.h>

namespace multem
{
	template<class T, eDevice dev>
	class Atomic_Cross_Section{
		public:
			using value_type_r = T;
			using value_type_c = complex<T>;

			Atomic_Cross_Section(): input_multislice(nullptr){}

			void set_input_data(Input_Multislice<value_type_r> *input_multislice_i)
			{
				input_multislice = input_multislice_i;
			}

			template<class TVector>
			void get(TVector &r_o, TVector &fr_o)
			{
				/**************************multislice calculation*******************************/
				multem::Output_Multislice_Vector<T> output_multislice;
				output_multislice.set_input_data(input_multislice);

				multem::Multislice<T, dev> multislice;
				multislice.set_input_data(input_multislice);

				multislice.run(output_multislice);

				multislice.cleanup();

				/******************************* Cross section ********************************/
				multem::Grid<double> grid;
				grid.assign(input_multislice->grid);

				multem::Cubic_Spline<double> spline;
				auto sigma = input_multislice->atoms.sigma[0];

				auto pr = &(output_multislice.scanning.r);
				std::vector<double> r2(pr->size());
				thrust::transform(pr->begin(), pr->end(), r2.begin(), [](T &x)->double{return x*x; });

				auto pfr2 = &(output_multislice.image_tot[0].image[0]);
				std::vector<double> fr2(pr->size());
				fr2.assign(pfr2->begin(), pfr2->begin()+fr2.size());

				multem::Vector<complex<double>, e_host> M(grid.nxy());
				T x = 0.5*grid.lx;
				T y = 0.5*grid.ly;

				// spline interpolation
				spline.set_points(r2, fr2);
				spline.eval_radial_function(grid, x, y, M);

				// Convolution
				multem::Stream<e_host> stream(input_multislice->cpu_nthread);

				multem::FFT2<double, e_host> fft2;
				fft2.create_plan_2d(grid.ny, grid.nx, input_multislice->cpu_nthread);

				multem::fft2_shift(stream, grid, M);
				multem::gaussian_convolution(stream, fft2, grid, sigma, M);
				multem::fft2_shift(stream, grid, M);
				fft2.cleanup();

				// extract radial values
				std::vector<double> r2_c(grid.nyh);
				std::vector<double> fr_c(grid.nyh);
				int ix = grid.nxh;
				int ir = 0;
				for(auto iy =grid.nyh; iy<grid.ny; iy++)
				{
					r2_c[ir] = grid.R2(ix, iy, x, y);
					fr_c[ir] = M[grid.ind_col(ix, iy)].real();
					ir++;
				}

				// spline interpolation
				spline.set_points(r2_c, fr_c);
				spline.eval_function(r2, fr2);

				r_o.assign(pr->begin(), pr->end());
				fr_o.assign(fr2.begin(), fr2.end());
			}

		private:
			Input_Multislice<value_type_r> *input_multislice;
	};

} // namespace multem

#endif