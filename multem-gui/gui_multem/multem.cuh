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

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "fft.cuh"
#include "stream.cuh"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "tem_simulation.cuh"

template <class T, mt::eDevice dev>
void run_multislice(mt::System_Configuration &system_conf,
	mt::Input_Multislice<T> &input_multislice, 
	mt::Output_Multislice<T> &output_multislice)
{
	mt::Input_Multislice<T> input_multislice;
	read_input_multislice(input_multislice);

	output_multislice.set_input_data(&input_multislice);

	mt::Stream<dev> stream(system_conf.nstream);
	mt::FFT<T, dev> fft_2d;
	fft_2d.create_plan_2d(input_multislice.grid_2d.ny, input_multislice.grid_2d.nx, system_conf.nstream);

	mt::Multislice<T, dev> tem_simulation;
	tem_simulation.set_input_data(&input_multislice, &stream, &fft_2d);
	tem_simulation(output_multislice);

	stream.synchronize();
	fft_2d.cleanup();
}