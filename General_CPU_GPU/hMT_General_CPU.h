/*
 * This file is part of MULTEM.
 * Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
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
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef hMT_GeneralCPU_H
#define hMT_GeneralCPU_H

#include "fftw3.h"
#include "hConstTypes.h"
#include "hMT_MGP_CPU.h"
#include "hMT_InMulSli_CPU.h"

// Grid's parameter initialization
void f_sGP_Init(sGP &GP);

// Grid's parameter calculation
void f_sGP_Cal(int nx, int ny, double lx, double ly, double dz, int BWL, int PBC_xy, sGP &GP);

// Grid's parameter calculation
void f_sGP_SetInputData(cMT_MGP_CPU *MT_MGP_CPU, sGP &GP);

/***************************************************************************/
/***************************************************************************/

// Block and Thread parameter initialization
void f_sLens_Init(sLens &Lens);

// Lens' parameter calculation
void f_sLens_Cal(double E0, sGP &GP, sLens &Lens);

// Set input data Lens' parameter
void f_sLens_SetInputData(cMT_InMulSli_CPU &MT_InMulSli_CPU, sGP &GP, sLens &Lens);

/***************************************************************************/
/***************************************************************************/

void f_ReadQuadratureCPU(int typ, int nQCPU, sQ1 &QCPU);

/***************************************************************************/
/***************************************************************************/

void f_sCoefPar_Free(sCoefPar &CoefPar);

void f_sCoefPar_Init(sCoefPar &CoefPar);

void f_sCoefPar_Malloc(int nCoefPar, sCoefPar &CoefPar);

/***************************************************************************/
/***************************************************************************/

void f_sciVn_Free(sciVn &ciVn);

void f_sciVn_Init(sciVn &ciVn);

void f_sciVn_Malloc(int nciVn, sciVn &ciVn);

/***************************************************************************/
/***************************************************************************/

void f_sDetCir_Free(sDetCir &DetCir);

void f_sDetCir_Init(sDetCir &DetCir);

void f_sDetCir_Malloc(int nDetCir, sDetCir &DetCir);

/***************************************************************************/
/***************************************************************************/
void f_scVp_Init(scVp &cVp);

void f_scVp_Init(int ncVp, scVp *cVp);

/***************************************************************************/
/***************************************************************************/

void f_Hadamard_Product_CPU(sGP &GP, double * __restrict A_i, double * __restrict B_io);

void f_Hadamard_Product_CPU(sGP &GP, double * __restrict A_i, double * __restrict B_i, double * __restrict C_o);

void f_Hadamard_Product_CPU(sGP &GP, fftw_complex * __restrict A_i, fftw_complex * __restrict B_io);

void f_Hadamard_Product_CPU(sGP &GP, fftw_complex * __restrict A_i, fftw_complex * __restrict B_i, fftw_complex * __restrict C_o);

/***************************************************************************/
/***************************************************************************/

// From double to float
void f_Double_2_Float_CPU(sGP &GP, double * __restrict V_i, float * __restrict V_o);

// From double to float
void f_Double_2_Float_CPU(sGP &GP, double * __restrict V1_i, double * __restrict V2_i, float * __restrict V1_o, float * __restrict V2_o);

// From float to double
void f_Float_2_Double_CPU(sGP &GP, float * __restrict V_i, double * __restrict V_o);

// From float to double
void f_Float_2_Double_CPU(sGP &GP, float * __restrict V1_i, float * __restrict V2_i, double * __restrict V1_o, double * __restrict V2_o);

/***************************************************************************/
/***************************************************************************/

void f_Scale_MD_CPU(sGP &GP, double w, double * __restrict MD_io);

void f_Scale_MD_CPU(sGP &GP, double w, double * __restrict MD1_io, double * __restrict MD2_io);

void f_Scale_MC_CPU(sGP &GP, double w, fftw_complex * __restrict MC_io);

void f_Scale_MC_CPU(sGP &GP, double w, fftw_complex * __restrict MC1_io, fftw_complex * __restrict MC2_io);

/***************************************************************************/
/***************************************************************************/

// Set value to Double vector:
void f_Set_MD_CPU(sGP &GP, double M, double * __restrict MD_o);

// Set value to 2 Double vector:
void f_Set_MD_CPU(sGP &GP, double M, double * __restrict MD1_o, double * __restrict MD2_o);

// Set value to Complex vector
void f_Set_MC_CPU(sGP &GP, double Mr, double Mi, fftw_complex * __restrict MC_o);

// Set value to 2 Double complex vector:
void f_Set_MC_CPU(sGP &GP, double Mr, double Mi, fftw_complex * __restrict MC1_o, fftw_complex * __restrict MC2_o);

// Set value to Complex vector and Double vector
void f_Set_MC_MD_CPU(sGP &GP, double Mr, double Mi, fftw_complex * __restrict MC_o, double M, double * __restrict MD_o);

// Set Real and Imaginary part of a Complex vector
void f_Set_MC_CPU(sGP &GP, sComplex &MC_i, fftw_complex * __restrict MC_o);

// Get Real and Imaginary part of a Complex vector
void f_Get_MC_CPU(sGP &GP, fftw_complex * __restrict MC_i, sComplex &MC_o);

/***************************************************************************/
/***************************************************************************/

void f_Set_wMD_CPU(sGP &GP, double w, double * __restrict MD_i, double * __restrict MD_io);

void f_Set_wMD_CPU(sGP &GP, double w, double * __restrict MD1_i, double * __restrict MD2_i, double * __restrict MD1_io, double * __restrict MD2_io);

void f_Set_wMC_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, fftw_complex * __restrict MC_io);

void f_Set_wMC_CPU(sGP &GP, double w, fftw_complex * __restrict MC1_i, fftw_complex * __restrict MC2_i, fftw_complex * __restrict MC1_io, fftw_complex * __restrict MC2_io);

void f_Set_wMC2_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, double * __restrict MD_io);

void f_Set_wMC_wMD_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, fftw_complex * __restrict MC_io, double * __restrict MD_io);

void f_Set_wMC_wMD_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, double * __restrict MD_i, fftw_complex * __restrict MC_io, double * __restrict MD_io);


void f_Add_wMD_CPU(sGP &GP, double w, double * __restrict MD_i, double * __restrict MD_io);

void f_Add_wMD_CPU(sGP &GP, double w, double * __restrict MD1_i, double * __restrict MD2_i, double * __restrict MD1_io, double * __restrict MD2_io);

void f_Add_wMC_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, fftw_complex * __restrict MC_io);

void f_Add_wMC_CPU(sGP &GP, double w, fftw_complex * __restrict MC1_i, fftw_complex * __restrict MC2_i, fftw_complex * __restrict MC1_io, fftw_complex * __restrict MC2_io);

void f_Add_wMC2_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, double * __restrict MD_io);

void f_Add_wMC_wMD_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, fftw_complex * __restrict MC_io, double * __restrict MD_io);

void f_Add_wMC_wMD_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, double * __restrict MD_i, fftw_complex * __restrict MC_io, double * __restrict MD_io);

/***************************************************************************/
/***************************************************************************/

// Shift Double matrix respect to (nxh, nyh)
void f_fft2Shift_MD_CPU(sGP &GP, double * __restrict MD_io);

// Shift Double matrix respect to (nxh, nyh)
void f_fft2Shift_MD_CPU(sGP &GP, double * __restrict MD1_io, double * __restrict MD2_io);

// Shift Complex matrix respect to (nxh, nyh)
void f_fft2Shift_MC_CPU(sGP &GP, fftw_complex * __restrict MC_io);

// Shift Complex matrix respect to (nxh, nyh)
void f_fft2Shift_MC_CPU(sGP &GP, fftw_complex * __restrict MC1_io, fftw_complex * __restrict MC2_io);

// Shift Complex matrix respect to (nxh, nyh)
void f_fft2Shift_MC_MD_CPU(sGP& GP, fftw_complex * __restrict MC_io, double * __restrict MD_io);

// Shift Complex matrix respect to (nxh, nyh)
void f_fft2Shift_MC_MD_CPU(sGP &GP, fftw_complex * __restrict MC_io, double * __restrict MD1_io, double * __restrict MD2_io);

// Shift Complex matrix respect to (nxh, nyh)
void f_fft2Shift_MC_CPU(sGP &GP, sComplex &MC_io);

/***************************************************************************/
/***************************************************************************/

double f_Sum_MD_CPU(sGP &GP, double w_i, const double * __restrict MD_i);

double f_Sum_MC2_CPU(sGP &GP, double w_i, const fftw_complex * __restrict MC_i);

void f_Sum_MD_Det_CPU(sGP &GP, double w_i, const double * __restrict MD_i, double gmin2, double gmax2, int iSum, double * __restrict Sum_MD_o);

void f_Sum_MD_Det_CPU(sGP &GP, double w_i, const double * __restrict MD1_i, const double * __restrict MD2_i, double gmin2, double gmax2, int iSum, double * __restrict Sum_MD1_o, double * __restrict Sum_MD2_o);

void f_Sum_MC_Det_CPU(sGP &GP, double w_i, const fftw_complex * __restrict MC_i, double gmin2, double gmax2, double * __restrict MDp_i, int iSum, double * __restrict Sum_MD_o);

/***************************************************************************/
/***************************************************************************/

// Phase
void f_Phase_CPU(sGP &GP, double gxu, double gyu, sACD &ExpR_x_o, sACD &ExpR_y_o);

// Phase multiplication
void f_PhaseMul_CPU(sGP &GP, double gxu, double gyu, sACD &ExpR_x_o, sACD &ExpR_y_o, fftw_complex *Psi_io);

// BandWidth Limit function
void f_BandwidthLimit2D_CPU(fftw_plan &PlanForward, fftw_plan &PlanBackward, sGP &GP, fftw_complex *MC_io);

// Build propagator function
void f_Propagator_CPU(sGP &GP, double gxu, double gyu, double scale, sACD &Prop_x_o, sACD &Prop_y_o);

// Propagate, scale with cut-off (2/3)gmax
void f_Propagate_CPU(fftw_plan &PlanForward, fftw_plan &PlanBackward, sGP &GP, eSpace Space, double gxu, double gyu, double lambda, double z, sACD &Prop_x_o, sACD &Prop_y_o, fftw_complex *Psi_io);

// Probe in Fourier space
void f_Probe_FS_CPU(sGP &GP, sLens Lens, double x, double y, fftw_complex *fPsi_o);

// Apply Coherent transfer function
void f_Apply_CTF_CPU(sGP &GP, sLens &Lens, double gxu, double gyu, fftw_complex *fPsi_io);

// Apply Coherent transfer function
void f_Apply_CTF_CPU(sGP &GP, sLens &Lens, double gxu, double gyu, fftw_complex *fPsi_i, fftw_complex *fPsi_o);

// Partially coherent transfer function, linear image model and weak phase object
void f_Apply_PCTF_CPU(sGP &GP, sLens &Lens, fftw_complex *fPsi_io);

// Partially coherent transfer function, linear image model and weak phase object
void f_Apply_PCTF_CPU(sGP &GP, sLens &Lens, fftw_complex *fPsi_i, fftw_complex *fPsi_o);

/***************************************************************************/
/***************************************************************************/

void f_Row_2_Column_Format_CPU(sGP &GP, double *&MC_i, double *&MC_o);

void f_Column_2_Row_Format_CPU(sGP &GP, double *&MC_i, double *&MC_o);

void f_fftw_complex_2_sComplex_CPU(sGP &GP, fftw_complex *&MC_i, sComplex &MC_o);

void f_sComplex_2_fftw_complex_CPU(sGP &GP, sComplex &MC_i, fftw_complex *&MC_o);

#endif