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
#include "fft.cuh"
#include "stream.cuh"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "tem_simulation.cuh"

template <class T, mt::eDevice dev>
void mt_run_multislice(mt::System_Configuration &system_conf,
 mt::Input_Multislice<T> &input_multislice, mt::Output_Multislice<T> &output_multislice)
{
	mt::Stream<dev> stream(system_conf.nstream);
	mt::FFT<T, dev> fft_2d;
	fft_2d.create_plan_2d(input_multislice.grid_2d.ny, input_multislice.grid_2d.nx, system_conf.nstream);

	if(input_multislice.is_IWFS_IWRS())
	{
		mt::Incident_Wave<T, dev> incident_wave;
		incident_wave.set_input_data(&input_multislice, &stream, &fft_2d);
	auto space = (input_multislice.is_IWRS())?mt::eS_Real:mt::eS_Reciprocal;
		incident_wave(space, output_multislice);
	}
	else
	{
		mt::Multislice<T, dev>::ext_stop_sim = false;

		mt::Multislice<T, dev> tem_simulation;
		tem_simulation.set_input_data(&input_multislice, &stream, &fft_2d);
		tem_simulation(output_multislice);
	}

	 stream.synchronize();
	 fft_2d.cleanup();
}

inline
void mt_stop_multislice(mt::System_Configuration &system_conf)
{
	if (system_conf.is_float_host())
	{
		mt::Multislice<float, mt::e_host>::ext_stop_sim = true;
	}
	else if (system_conf.is_double_host())
	{
		mt::Multislice<double, mt::e_host>::ext_stop_sim = true;
	}
	else if (system_conf.is_float_device())
	{
		mt::Multislice<float, mt::e_device>::ext_stop_sim = true;
	}
	else if (system_conf.is_double_device())
	{
		mt::Multislice<double, mt::e_device>::ext_stop_sim = true;
	}
}

inline
bool mt_success_multislice(mt::System_Configuration &system_conf)
{
	if (system_conf.is_float_host())
	{
		return !mt::Multislice<float, mt::e_host>::ext_stop_sim;
	}
	else if (system_conf.is_double_host())
	{
		return !mt::Multislice<double, mt::e_host>::ext_stop_sim;
	}
	else if (system_conf.is_float_device())
	{
		return !mt::Multislice<float, mt::e_device>::ext_stop_sim;
	}
	else if (system_conf.is_double_device())
	{
		return !mt::Multislice<double, mt::e_device>::ext_stop_sim;
	}
}

inline
int mt_niter(mt::System_Configuration &system_conf)
{
	if (system_conf.is_float_host())
	{
		return mt::Multislice<float, mt::e_host>::ext_niter;
	}
	else if (system_conf.is_double_host())
	{
		return mt::Multislice<double, mt::e_host>::ext_niter;
	}
	else if (system_conf.is_float_device())
	{
		return mt::Multislice<float, mt::e_device>::ext_niter;
	}
	else if (system_conf.is_double_device())
	{
		return mt::Multislice<double, mt::e_device>::ext_niter;
	}
}

inline
int mt_iter(mt::System_Configuration &system_conf)
{
	if (system_conf.is_float_host())
	{
		return mt::Multislice<float, mt::e_host>::ext_iter;
	}
	else if (system_conf.is_double_host())
	{
		return mt::Multislice<double, mt::e_host>::ext_iter;
	}
	else if (system_conf.is_float_device())
	{
		return mt::Multislice<float, mt::e_device>::ext_iter;
	}
	else if (system_conf.is_double_device())
	{
		return mt::Multislice<double, mt::e_device>::ext_iter;
	}
}

inline
void mt_init()
{
	mt::Stream<mt::e_device> stream;
	mt::FFT<float, mt::e_device> fft_2d;
}