#include "math.h"
#include "hConstTypes.h"
#include "hgeneralGPU.h"
#include "hSTEMGPU.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include <device_functions.h>

/****************************************************************************/
/****************************************************************************/
void cSTEMGPU::PowerGrid(double x0, double xe, int np, int n, double *x){
	int i, nh = (n-1)/2+1;
	double xmax = 0.5*(xe-x0);
	double t, tt, tmin = 0, tmax = pow(xmax, 1.0/np), dt = (tmax-tmin)/(nh-1);

	for (i=0; i<nh; i++){
		t = tmin + dt*i; tt = pow(t, np);
		x[nh-1+i] = x0 + xmax + tt;
		if (i>0)
			x[nh-1-i] = x0 + xmax - tt;
	}
}

void cSTEMGPU::BuildGrid(int typ, int line, double x0s, double y0s, double xes, double yes){
	int i, np = 2;
	double lxs = xes-x0s, lys = yes-y0s, ls = sqrt(lxs*lxs+lys*lys);
	double ds, theta = atan(lys/lxs);
	double costheta = cos(theta), sintheta = sin(theta);
	if (ns%2==0)
		ns = ns + 1;

	if (line){
			nxs = ns;
			nys = ns;
	}else{
		if (abs(lxs)>abs(lys)){
			nxs = ns;
			nys = (int)floor(ns*abs(lys/lxs));
		}
		else{
			nxs = (int)floor(ns*abs(lxs/lys));
			nys = ns;
		}
	}

	delete [] ax;
	delete [] ay;
	delete [] als;

	ax = new double [nxs];
	ay = new double [nys];
	als = new double [ns];

	switch(typ){
		case 0:
			if (line){
				ds = ls/(ns-1);

				for (i=0; i<ns; i++){
					als[i] = i*ds;
					ax[i] = x0s + als[i]*costheta;
					ay[i] = y0s + als[i]*sintheta;
				}
			}else{
				dxs = lxs/(nxs-1);
				for (i=0; i<nxs; i++)
					ax[i] = x0s + i*dxs;

				dys = lys/(nys-1);
				for (i=0; i<nys; i++)
					ay[i] = y0s + i*dys;
			}
			break;
		case 1:
			if (line){
				PowerGrid(0, ls, np, ns, als);
				for (i=0; i<ns; i++){
					ax[i] = x0s + als[i]*costheta;
					ay[i] = y0s + als[i]*sintheta;
				}
			}else{
				PowerGrid(x0s, xes, np, nxs, ax);
				PowerGrid(y0s, yes, np, nys, ay);
			}
			break;
	}
}

void cSTEMGPU::GetIndDet(){
	int i, j, *indD;
	int nxh = nx/2, nyh = ny/2;
	double gx, gy, g2;
	double g2min = Det.g2min, g2max = Det.g2max;

	indD = new int [nx*ny];
	nIndDet = 0;
	for (i=0; i<nx; i++){	
		gx = ((i<nxh)?i:(i-nx))*dgx;
		for (j=0; j<ny; j++){	
			gy = ((j<nyh)?j:(j-ny))*dgy;
			g2 = gx*gx + gy*gy;
			if ((g2min <= g2)&&(g2 <= g2max))
				indD[nIndDet++] = i*ny+j;
		}
	}
	cudaFreen(IndDet);
	cudaMalloc((void**)&IndDet, nIndDet*cSizeofI);
	cudaMemcpy(IndDet, indD, nIndDet*cSizeofI, cudaMemcpyHostToDevice);
	delete [] indD;
}

void cSTEMGPU::DeriveParameters(){
	double a = lx/nucx, b = ly/nucy;
	double x0 = (x1u-0.5)*a, y0 = (y1u-0.5)*b, xe = (x2u-0.5)*a, ye = (y2u-0.5)*b;

	BuildGrid(GridType, line, x0, y0, xe, ye);

	GetIndDet();

	nxy = nx*ny;
}

/******************************************************************************/
/******************************************************************************/

cSTEMGPU::cSTEMGPU(){
	Psig = 0;
	nIndDet = 0;
	IndDet = 0;
	Sd = 0;
	ax = 0;
	ay = 0;
	als = 0;
	GridType = 0;	// regular grid
	line = 0;		// area
}

void cSTEMGPU::freeMemory(){
	Psig = 0;
	cudaFreen(IndDet);
	cudaFreen(Sd);
	delete [] ax; ax = 0;
	delete [] ay; ay = 0;
	delete [] als; als = 0;
}

void cSTEMGPU::GenerateParameters(int nx_i, int ny_i, double dgx_i, double dgy_i, double2 *&Psig_i){
	nx = nx_i;
	ny = ny_i;
	dgx = dgx_i;
	dgy = dgy_i;
	Psig = Psig_i;

	DeriveParameters();

	BnDet.x = MIN(64, (nIndDet+thrnxy-1)/thrnxy); BnDet.y = 1; BnDet.z = 1;
	TnDet.x = thrnxy; TnDet.y = 1; TnDet.z = 1;

	cudaFreen(Sd);
	cudaMalloc((void**)&Sd, 2*thrnxy*cSizeofRD);
}

double cSTEMGPU::DetIntensity(){
	return SumAc2Ind(BnDet, TnDet, nIndDet, IndDet, Psig, Sd);
}