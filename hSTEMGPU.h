#ifndef hSTEMGPU_H
#define hSTEMGPU_H

#include "hConstTypes.h"
#include "cuda.h"

/**************************STEM*************************/
class cSTEMGPU{
	private:
		int nx;
		int ny;
		double dgx; 
		double dgy;
		double2 *Psig;

		/********************Calculated data******************/
		int nxy;

		int nIndDet;
		int *IndDet;
		double *Sd;

		double dxs;			// x-pixel size
		double dys;			// y-pixel size
		double x0s, y0s;	// Initial scan coordinates respect to the total area (x0, y0)

		dim3 BnDet;
		dim3 TnDet;

		void PowerGrid(double x0, double xe, int np, int n, double *x);
		void BuildGrid(int typ, int line, double x0, double y0, double xe, double ye);
		void GetIndDet();
		void DeriveParameters();
	public:
		int ns;		// Number of scan points
		double lx;	// a distance
		double ly;	// b distance
		int nucx;	// x-number of unit cell per supercells
		int nucy;	// y-number of unit cell per supercells
		double x1u, y1u;	// scan coordinates (x1, y1) units of unit cell
		double x2u, y2u;	// scan coordinates (x2, y2) units of unit cell
		bool line;			// True: Line; False: area
		sDetCir Det;		// Detector
		int GridType;		// Grid type

		int nxs;			// Number of scan points x
		int nys;			// Number of scan points x

		double *ax;			// array x-direction
		double *ay;			// array y-direction
		double *als;		// array y-direction

		void GenerateParameters(int nx_i, int ny_i, double dgx_i, double dgy_i, double2 *&Psig_i);

		double DetIntensity();

		cSTEMGPU();
		void freeMemory();
};

#endif