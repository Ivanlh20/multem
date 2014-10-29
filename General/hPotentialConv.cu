#include "hConstTypes.h"
#include "hPotentialGPU.h"
#include "hQuadrature.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>

void GetQuadratureParameters(int nr, double *r, int nz, double sigma, sQp Qp){
	double xmin, xmax, x0n, xl;
	double sigman, a, b;
	double tmin, tmax;

	xmin = -1+1e-14; xmax = 1.0-1e-05;
	for (int i = 0; i<nr; i++){
		Qp.r[i] = r[i];
		Qp.rlim[i] = xl = r[i] + 7.5*sigma;
		x0n = -1.0 + 2.0*r[i]/xl;
		sigman = 2.0*sigma/xl;
		Qp.b[i] = b = 0.08096*sinh(8.545*sigman);
		Qp.a[i] = a = atanh(x0n)-(1.0-sqrt(abs(1.0-SQR(4.0*x0n*b))))/(4.0*x0n); 
		Qp.tmin[i] = tmin = asinh((atanh(xmin)-a)/b);
		tmax = asinh((atanh(xmax)-a)/b);
		Qp.h[i] = (tmax-tmin)/(nz-1);
	}
}

void PotentialConv(sAtomTypesCPU AtomTypesCPU, int nr, double *rh, int n, double *Vrh){

	/**********************Quadratures R0************************/
	sQ1 QRh, QR;
	QRh.x = new double [stnR0c]; 
	QRh.w = new double [stnR0c];
	cQuadrature Quadrature;
	Quadrature.ReadQuadrature(2, stnR0c, QRh); // ExpExp Quadrature

	cudaMalloc((void**)&QR.x, stnR0c*cSizeofRD);
	cudaMalloc((void**)&QR.w, stnR0c*cSizeofRD);

	cudaMemcpy(QR.x, QRh.x, stnR0c*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(QR.w, QRh.w, stnR0c*cSizeofRD, cudaMemcpyHostToDevice);

	delete [] QRh.x;
	delete [] QRh.w;

	/******************* Quadratures Parameters z0 ******************/
	double sigma = 0; //AtomTypesCPU.sigma;//////////////////////
	double alpha = 1.0/(2.0*SQR(sigma));

	sQp Qph, Qp;
	Qph.r = new double [nr];
	Qph.rlim = new double [nr];
	Qph.a = new double [nr];
	Qph.b = new double [nr];
	Qph.tmin = new double [nr];
	Qph.h = new double [nr];

	GetQuadratureParameters(nr, rh, stnz0c, sigma, Qph);

	cudaMalloc((void**)&Qp.r, nr*cSizeofRD);
	cudaMalloc((void**)&Qp.rlim, nr*cSizeofRD); 
	cudaMalloc((void**)&Qp.a, nr*cSizeofRD);
	cudaMalloc((void**)&Qp.b, nr*cSizeofRD);
	cudaMalloc((void**)&Qp.tmin, nr*cSizeofRD);
	cudaMalloc((void**)&Qp.h, nr*cSizeofRD);

	cudaMemcpy(Qp.r, Qph.r, nr*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(Qp.rlim, Qph.rlim, nr*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(Qp.a, Qph.a, nr*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(Qp.b, Qph.b, nr*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(Qp.tmin, Qph.tmin, nr*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(Qp.h, Qph.h, nr*cSizeofRD, cudaMemcpyHostToDevice);

	delete [] Qph.r;
	delete [] Qph.rlim;
	delete [] Qph.a;
	delete [] Qph.b;
	delete [] Qph.tmin;
	delete [] Qph.h;

	/***************************** r ******************************/
	double *cl, *cnl; 
	cudaMalloc((void**)&cl, 6*cSizeofRD);
	cudaMalloc((void**)&cnl, 6*cSizeofRD);

	cudaMemcpy(cl, AtomTypesCPU.cVr.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(cnl, AtomTypesCPU.cVr.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	double *Vr, f = CUB(1.0/(sigma*sqrt(c2Pi)));

	cudaMalloc((void**)&Vr, nr*cSizeofRD);
	//BlurredPotential(nr, stnz0c, AtomTypesCPU.PotPar, alpha, Qp, stnz0c, stnR0c, QR, cl, cnl, n, f, Vr);

	cudaMemcpy(Vrh, Vr, nr*cSizeofRD, cudaMemcpyDeviceToHost);

	cudaFreen(QR.x);
	cudaFreen(QR.w);
	cudaFreen(cl);
	cudaFreen(cnl);
	cudaFreen(Qp.r);
	cudaFreen(Qp.rlim);
	cudaFreen(Qp.a);
	cudaFreen(Qp.b);
	cudaFreen(Qp.tmin);
	cudaFreen(Qp.h);
	cudaFreen(Vr);
}