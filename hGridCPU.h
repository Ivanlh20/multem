#ifndef hGridCPU_H
#define hGridCPU_H

class cGridCPU{
	private:
		void RegularGrid(double rmin, double rmax, int nr, double *r);
		void LogarithmicGrid(double rmin, double rmax, int nr, double *r);
		void LogarithmicGridShifted(double rmin, double rmax, int nr, double *r);
	public:
		void ReadGrid(double rmin, double rmax, int nr, int gridTyp, double *r);
};

#endif
