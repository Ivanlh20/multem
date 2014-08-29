#include "hmathCPU.h"
#include "hGridCPU.h"

// Build logarithm grid
void cGridCPU::RegularGrid(double rmin, double rmax, int nr, double *r){
	double dr = (rmax-rmin)/(nr-1);
	for (int i=0; i<nr; i++)
		r[i] = rmin + i*dr;
}

// Build logarithm grid
void cGridCPU::LogarithmicGrid(double rmin, double rmax, int nr, double *r){
	double dlnr = log(rmax/rmin)/double(nr-1);
	for (int i=0; i<nr; i++)
		r[i] = rmin*exp(i*dlnr);
}

// Build logarithm grid
void cGridCPU::LogarithmicGridShifted(double rmin, double rmax, int nr, double *r){
	double dlnr = log(rmax/rmin + 1.0)/double(nr-1);
	for (int i=0; i<nr; i++)
		r[i] = rmin*(exp(i*dlnr)-1.0);
}

// Read logarithm grid
void cGridCPU::ReadGrid(double rmin, double rmax, int nr, int gridTyp, double *r){
	switch (gridTyp){
		case 0:
			RegularGrid(rmin, rmax, nr, r);
			break;
		case 1:
			LogarithmicGrid(rmin, rmax, nr, r);
			break;
		case 2:
			LogarithmicGridShifted(rmin, rmax, nr, r);
			break;
	}
}