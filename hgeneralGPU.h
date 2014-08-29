#ifndef hgeneralGPU_H
#define hgeneralGPU_H

#include "C:\cuda\include\vector_types.h"

__device__ unsigned int stRetCount = 0;
__device__ double stTotSumd = 0;
__device__ double2 stTotSumc;

void ScaleVectorD(dim3 grid, dim3 threads, int n, double f, double *Vg);

void ScaleVectorC(dim3 grid, dim3 threads, int n, double f, double2 *Vg);

void SetValueVectorD(dim3 grid, dim3 threads, int n, double v, double *Vg1, double *Vg2);

void SetValueVectorD(dim3 grid, dim3 threads, int n, double v, double *Vg);

void SetValueVectorC(dim3 grid, dim3 threads, int n, double vr, double vi, double2 *Vg);

void SetValueVectorCD(dim3 grid, dim3 threads, int n, double vr, double vi, double2 *Vcg, double v, double *Vdg);

void GetRealImagVectorC(dim3 grid, dim3 threads, int n, double2 *Vg, double *Vrg, double *Vig);

void SetRealImagVectorC(dim3 grid, dim3 threads, int n, double *Vrg, double *Vig, double2 *Vg);

void AddwDtoD(dim3 grid, dim3 threads, bool add, int n, double w, double *Vig, double *Vog);

void AddwM2CtoD(dim3 grid, dim3 threads, bool add, int n, double w, double2 *Vg, double *Vg2);

void AddwCtoC(dim3 grid, dim3 threads, bool add, int n, double w, double2 *Vig, double2 *Vog);

void AddwCwM2CtoCD(dim3 grid, dim3 threads, bool add, int n, double w, double2 *Vg, double2 *Vga, double *Vg2a);

void AddwCwDtoCD(dim3 grid, dim3 threads, bool add, int n, double w, double2 *Vg, double *Vg2, double2 *Vga, double *Vg2a);

void fft2ShiftD(dim3 grid, dim3 threads, int nxh, int nyh, double *Mg);

void fft2ShiftC(dim3 grid, dim3 threads, int nxh, int nyh, double2 *Mg);

void SetRetCount(unsigned int ih=0);

double SumAd(dim3 grid, dim3 threads, unsigned int nxy, const double *A, double *TB);

double SumAc2(dim3 grid, dim3 threads, unsigned int nxy, const double2 *A, double *TB);

double SumAc2Ind(dim3 grid, dim3 threads, int nind, int *ind, const double2 *A, double *TB);

void ReadQuadratureGPU(int typ, int nQGPU, sQ1 &QGPU);

void ReadQuadratureCPUGPU(int typ, int nQ, sQ1 &QCPU, sQ1 &QGPU);

#endif