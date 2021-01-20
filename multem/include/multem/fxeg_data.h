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

#ifndef MULTEM_FXEG_DATA_H
#define MULTEM_FXEG_DATA_H

#include <vector>
#include <multem/config.h>

namespace mt {

  class fxeg_Tabulated_Data{
    public:

      DLL_PUBLIC void ReadTabData(int Z_i, int Type_i, int dng_i, int &n_o, double *g_o, double *g2_o, double *fx_o, double *fe_o);

    private:

      void fxegActaCrys(int Z);
      void fxegRez(int Z);
      void fxegKirkland(int Z);
      void fxegLobato(int Z);

      std::vector<double> g;
      std::vector<double> g2;
      std::vector<double> feg;
      std::vector<double> fxg;
  };

}

#endif

