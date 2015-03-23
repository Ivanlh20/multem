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

#ifndef hMT_General_GPU_H
#define hMT_General_GPU_H

#include "hConstTypes.h"
#include <cuda.h>
#include <cufft.h>

/***************************************************************************/
/***************************************************************************/

void f_GPU_Sync_CPU(int &iSynCPU, int &cSynCPU);

/***************************************************************************/
/***************************************************************************/

void f_get_BTnxny(sGP &GP, dim3 &B, dim3 &T);

void f_get_BTnxhnyh(sGP &GP, dim3 &B, dim3 &T);

void f_get_BTmnxny(sGP &GP, dim3 &B, dim3 &T);

void f_get_BTnxy(sGP &GP, dim3 &B, dim3 &T);

/***************************************************************************/
/***************************************************************************/

void f_ReadQuadratureGPU(int typ, int nQGPU, sQ1 &QGPU);

void f_ReadQuadratureCPUGPU(int typ, int nQ, sQ1 &QCPU, sQ1 &QGPU);

/***************************************************************************/
/***************************************************************************/

void f_sCoefPar_Free_GPU(sCoefPar &CoefPar);

void f_sCoefPar_Init_GPU(sCoefPar &CoefPar);

void f_sCoefPar_Malloc_GPU(int nCoefPar, sCoefPar &CoefPar);

/***************************************************************************/
/***************************************************************************/

void f_sciVn_Free_GPU(sciVn &ciVn);

void f_sciVn_Init_GPU(sciVn &ciVn);

void f_sciVn_Malloc_GPU(int nciVn, sciVn &ciVn);

/***************************************************************************/
/***************************************************************************/

void f_sDetCir_Free_GPU(sDetCir &DetCir);

void f_sDetCir_Init_GPU(sDetCir &DetCir);

void f_sDetCir_Malloc_GPU(int nDetCir, sDetCir &DetCir);

/***************************************************************************/
/***************************************************************************/

void f_sACD_Free_GPU(sACD &ACD);

void f_sACD_Init_GPU(sACD &ACD);

void f_sACD_Malloc_GPU(int nACD, sACD &ACD);

/***************************************************************************/
/***************************************************************************/

// Hadamard product
void f_Hadamard_Product_GPU(sGP &GP, double *&A_i, double *&B_io);

// Hadamard product
void f_Hadamard_Product_GPU(sGP &GP, double *&A_i, double *&B_i, double *&C_o);

// Hadamard product
void f_Hadamard_Product_GPU(sGP &GP, double2 *&A_i, double2 *&B_io);

// Hadamard product
void f_Hadamard_Product_GPU(sGP &GP, double2 *&A_i, double2 *&B_i, double2 *&C_o);

/***************************************************************************/
/***************************************************************************/

void f_Double_2_Float_GPU(sGP &GP, double *&V_i, float *&V_o);

void f_Double_2_Float_GPU(sGP &GP, double *&V1_i, double *&V2_i, float *&V1_o, float *&V2_o);

void f_Float_2_Double_GPU(sGP &GP, float *&V_i, double *&V_o);

void f_Float_2_Double_GPU(sGP &GP, float *&V1_i, float *&V2_i, double *&V1_o, double *&V2_o);

/***************************************************************************/
/***************************************************************************/

void f_Scale_MD_GPU(sGP &GP, double f, double *&MD_io);

void f_Scale_MD_GPU(sGP &GP, double f, double *&MD1_io, double *&MD2_io);

void f_Scale_MC_GPU(sGP &GP, double f, double2 *&MC_io);

void f_Scale_MC_GPU(sGP &GP, double f, double2 *&MC1_io, double2 *&MC2_io);

/***************************************************************************/
/***************************************************************************/

void f_Set_MD_GPU(sGP &GP, double M, double *&MD_o);

void f_Set_MD_GPU(sGP &GP, double M, double *&MD1_o, double *&MD2_o);

void f_Set_MC_GPU(sGP &GP, double Mr, double Mi, double2 *&MC_o);

void f_Set_MC_GPU(sGP &GP, double Mr, double Mi, double2 *&MC1_o, double2 *&MC2_o);

void f_Set_MC_MD_GPU(sGP &GP, double Mr, double Mi, double2 *&MC_o, double M, double *&MD_o);

/***************************************************************************/
/***************************************************************************/

void f_Set_wMD_GPU(sGP &GP, double w, double *&MD_i, double *&MD_io);

void f_Set_wMD_GPU(sGP &GP, double w, double *&MD1_i, double *&MD2_i, double *&MD1_io, double *&MD2_io);

void f_Set_wMC_GPU(sGP &GP, double w, double2 *&MC_i, double2 *&MC_io);

void f_Set_wMC_GPU(sGP &GP, double w, double2 *&MC1_i, double2 *&MC2_i, double2 *&MC1_io, double2 *&MC2_io);

void f_Set_wMC2_GPU(sGP &GP, double w, double2 *&MC_i, double *&MD_io);

void f_Set_wMC_wMD_GPU(sGP &GP, double w, double2 *&MC_i, double2 *&MC_io, double *&MD_io);

void f_Set_wMC_wMD_GPU(sGP &GP, double w, double2 *&MC_i, double *&MD_i, double2 *&MC_io, double *&MD_io);


void f_Add_wMD_GPU(sGP &GP, double w, double *&MD_i, double *&MD_io);

void f_Add_wMD_GPU(sGP &GP, double w, double *&MD1_i, double *&MD2_i, double *&MD1_io, double *&MD2_io);

void f_Add_wMC_GPU(sGP &GP, double w, double2 *&MC_i, double2 *&MC_io);

void f_Add_wMC_GPU(sGP &GP, double w, double2 *&MC1_i, double2 *&MC2_i, double2 *&MC1_io, double2 *&MC2_io);

void f_Add_wMC2_GPU(sGP &GP, double w, double2 *&MC_i, double *&MD_io);

void f_Add_wMC_wMD_GPU(sGP &GP, double w, double2 *&MC_i, double2 *&MC_io, double *&MD_io);

void f_Add_wMC_wMD_GPU(sGP &GP, double w, double2 *&MC_i, double *&MD_i, double2 *&MC_io, double *&MD_io);

/***************************************************************************/
/***************************************************************************/

void f_fft2Shift_MD_GPU(sGP &GP, double *&MD_io);

void f_fft2Shift_MD_GPU(sGP &GP, double *&MD1_io, double *&MD2_io);

void f_fft2Shift_MC_GPU(sGP &GP, double2 *&MC_io);

void f_fft2Shift_MC_GPU(sGP &GP, double2 *&MC1_io, double2 *&MC2_io);

void f_fft2Shift_MC_MD_GPU(sGP &GP, double2 *&MC_io, double *&MD_io);

void f_fft2Shift_MC_MD_GPU(sGP &GP, double2 *&MC_io, double *&MD1_io, double *&MD2_io);

void f_fft2Shift_MC_GPU(sGP &GP, sComplex &MC_io);

/***************************************************************************/
/***************************************************************************/

double f_Sum_MD_GPU(sGP &GP, double w_i, double *&MD_i, double *&MDp_i);

double f_Sum_MC2_GPU(sGP &GP, double w_i, double2 *&MC_i, double *&MDp_i);

void f_Sum_MD_Det_GPU(sGP &GP, double w_i, double *&MD_i, double gmin2, double gmax2, double *&MDp_i, int iSum, double *&Sum_MD_o);

void f_Sum_MD_Det_GPU(sGP &GP, double w_i, double *&MD1_i, double *&MD2_i, double gmin2, double gmax2, double *&MD1p_i, double *&MD2p_i, int iSum, double *&Sum_MD1_o, double *&Sum_MD2_o);

void f_Sum_MC_Det_GPU(sGP &GP, double w_i, double2 *&MC_i, double gmin2, double gmax2, double *&MDp_i, int iSum, double *&Sum_MD_o);

/***************************************************************************/
/***************************************************************************/

void f_PhaseMul_GPU(sGP &GP, double gxu, double gyu, sACD &ExpRg_x_o, sACD &ExpRg_y_o, double2 *&Psi_io);

void f_BandwidthLimit2D_GPU(cufftHandle &PlanPsi, sGP &GP, double2 *&MC_io);

void f_Propagate_GPU(cufftHandle &PlanPsi, sGP &GP, eSpace Space, double gxu, double gyu, double lambda, double z, sACD &Prop_x_o, sACD &Prop_y_o, double2 *&Psi_io);

void f_Probe_FS_GPU(sGP &GP, sLens Lens, double x, double y, double2 *&fPsi_o);

void f_Apply_CTF_GPU(sGP &GP, sLens &Lens, double gxu, double gyu, double2 *&fPsi_io);

void f_Apply_CTF_GPU(sGP &GP, sLens &Lens, double gxu, double gyu, double2 *&fPsi_i, double2 *&fPsi_o);

void f_Apply_PCTF_GPU(sGP &GP, sLens &Lens, double2 *&fPsi_io);

void f_Apply_PCTF_GPU(sGP &GP, sLens &Lens, double2 *&fPsi_i, double2 *&fPsi_o);

/***************************************************************************/
/***************************************************************************/

void f_cuDouble_2_double_GPU(sGP &GP, double *&MD_d_i, double *&MD_h_o);

void f_double_2_cuDouble_GPU(sGP &GP, double *&MD_h_i, double *&MD_d_o);

void f_cuDoubleComplex_2_sComplex_GPU(sGP &GP, double2 *&MC_d_i, double2 *&MC_h_t, sComplex &MC_h_o);

void f_sComplex_2_cuDoubleComplex_GPU(sGP &GP, sComplex &MC_h_i, double2 *&MC_h_t, double2 *&MC_d_o);

#endif