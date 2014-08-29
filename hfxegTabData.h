#ifndef hfxegTabData_H
#define hfxegTabData_H

class cfxegTabData{
	private:
		int nmax;
		int ng;
		double *g;
		double *g2;
		double *feg;
		double *fxg;
		void fxegActaCrys(int Z);
		void fxegRez(int Z);
		void fxegKirkland(int Z);
		void fxegLobato(int Z);
	public:
		cfxegTabData();
		~cfxegTabData();
		void ReadTabData(int Zi, int Typi, int dni, int &no, double *go, double *g2o, double *fxo, double *feo);
};

#endif