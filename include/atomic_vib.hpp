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

#ifndef ATOMIC_VIB_H
	#define ATOMIC_VIB_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include "const_enum_mt.cuh"
	#include "math_mt.h"
	#include "r_3d.h"

	namespace mt
	{
		/**************** atomic vibration model parameters ******************/
		class Atomic_Vib
		{
			public:
				eAtomic_Vib_Mod model;	// 1: Still atom model, 2: Absorptive potential model, 3: Frozen phonon model
				dt_bool coh_contrib;	// true, false		
				dt_int32 dist;			// 1: Gaussian (Phonon distribution)
				dt_int32 seed;			// random seed (frozen phonon)
				dt_bool sgl_conf;		// single configuration: true, false	
				dt_int32 nconf;			// if single phonon configuration == true then phonon configuration if not number of frozen phonon configurations
				R_3d<dt_bool> dim;		// atomic vibration dimension (x, y, z)

				dt_int32 iconf_0;		// initial configuration
				dt_int32 iconf_e;		// final configuration

				Atomic_Vib(): model(eavm_still_atom), coh_contrib(false), sgl_conf(false), dist(1), seed(300183), 
				nconf(1), iconf_0(1), iconf_e(1), dim{true, true, false}{}

				void assign(Atomic_Vib& atomic_vib)
				{
					if (this !=  &atomic_vib)
					{
						model = atomic_vib.model;
						coh_contrib = atomic_vib.coh_contrib;
						sgl_conf = atomic_vib.sgl_conf;
						dist = atomic_vib.dist;
						seed = atomic_vib.seed;
						nconf = atomic_vib.nconf;
						iconf_0 = atomic_vib.iconf_0;
						dim = atomic_vib.dim;
					}
				}

				void set_in_data(const eAtomic_Vib_Mod& model, const dt_bool& coh_contrib, const dt_bool& sgl_conf, const dt_int32& dist, 
					const dt_int32& seed, const dt_int32& nconf, const dt_int32& iconf_0, const R_3d<dt_bool>& dim)
				{
					this->model = model;
					this->coh_contrib = coh_contrib;
					this->sgl_conf = sgl_conf;
					this->dist = dist;
					this->seed = seed;
					this->nconf = nconf;
					this->iconf_0 = iconf_0;
					this->dim = dim;

					set_dep_var();
				}

				void set_dep_var()
				{
					seed = max(1, seed);

					if (is_avm_frozen_phonon())
					{
						nconf = max(1, nconf);
						iconf_0 = (sgl_conf)?nconf:1;
						iconf_e = nconf;
					}
					else
					{
						iconf_0 = iconf_e = nconf = 1;
					}
				}

				Atomic_Vib& operator=(Atomic_Vib& atomic_vib)
				{
					assign(atomic_vib);
					return *this;
				}

				dt_bool is_avm_still_atom() const
				{
					return mt::is_avm_still_atom(model);
				}

				dt_bool is_avm_absorptive_pot() const
				{
					return  mt::is_avm_absorptive_pot(model);
				}

				dt_bool is_avm_frozen_phonon() const
				{
					return  mt::is_avm_frozen_phonon(model);
				}

				dt_bool is_avm_frozen_phonon_sgl_conf() const
				{
					return  mt::is_avm_frozen_phonon_sgl_conf(model, sgl_conf);
				}

				dt_int32 number_conf()
				{
					return (nconf-iconf_0+1);
				}
		};
	}

#endif