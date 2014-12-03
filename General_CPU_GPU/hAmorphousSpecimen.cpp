#include "math.h"
#include <cstring>
#include "hConstTypes.h"
#include "hRandGen.h"

#include "hAmorphousSpecimen.h"

cAmorphousSpecimen::cAmorphousSpecimen(){
	c3d = 0;
	x = 0;
	y = 0;
	z = 0;
}

void cAmorphousSpecimen::freeMemory(){
	delete [] c3d; c3d = 0;
	delete [] x; x = 0;
	delete [] y; y = 0;
	delete [] z; z = 0;
}

cAmorphousSpecimen::~cAmorphousSpecimen(){
	freeMemory();
}

void cAmorphousSpecimen::SetInputData(double lxi, double lyi, double lzi, double dmini, int Zi){
	RandGen.reset();

	lx = lxi;
	ly = lyi;
	lz = lzi;
	Z = Zi;

	dmin = dmini;
	dmin2 = dmin*dmin;
	rmin = 0.5*dmin;
	Vsmin = 4.0*cPi*rmin*rmin*rmin/3.0;
	na = (int)floor((lx+dmin)*(ly+dmin)*(lz+dmin)/Vsmin);
	delete [] x; x = new double[na];
	delete [] y; y = new double[na];
	delete [] z; z = new double[na];

	lc = dmin/c3i2;
	ilc = 1.0/lc;
	ncx = (int)floor(lx*ilc+1); 
	ncy = (int)floor(ly*ilc+1);
	ncz = (int)floor(lz*ilc+1),
	ncxy = ncx*ncy;
	ncxyz = ncxy*ncz;

	delete [] c3d; c3d = new int[ncxyz];
	for (int i=0; i<ncxyz; i++)
		c3d[i] = -1;

}

void cAmorphousSpecimen::getLimits(int k, int kd, int nk, int &k0, int &ke){
	k0 = k-kd; if (k0<0) k0 = 0;
	ke = k+kd+1; if (ke>nk) ke = nk;
}

int cAmorphousSpecimen::CheckDistance(double &xn, double &yn, double &zn){
	int i, j, k;
	int i0, j0, k0, ie, je, ke;
	int p, pt;
	double xd, yd, zd, d2;

	i = (int)floor(xn*ilc);
	j = (int)floor(yn*ilc);
	k = (int)floor(zn*ilc);
	p = i+j*ncx+k*ncxy;

	if (c3d[p]>-1)
		return -1;

	getLimits(i, 2, ncx, i0, ie);
	getLimits(j, 2, ncy, j0, je);
	getLimits(k, 2, ncz, k0, ke);

	for (i=i0; i<ie; i++)
		for (j=j0; j<je; j++)
			for (k=k0; k<ke; k++){
				pt = c3d[i+j*ncx+k*ncxy];
				if (pt>-1){	
					xd = x[pt]-xn; yd = y[pt]-yn; zd = z[pt]-zn;
					d2 = xd*xd + yd*yd + zd*zd;
					if (d2<dmin2)
						return -1;
				}
			}

	return p;
}

int cAmorphousSpecimen::GeneratePoint(double &xn, double &yn, double &zn){
	int i, kc;
	for(i=0; i<ng; i++){
		xn = lx*RandGen.randu(); yn = ly*RandGen.randu(); zn = lz*RandGen.randu();
		kc = CheckDistance(xn, yn, zn);
		if (kc>-1) break;		
	}
	return kc;
}

int cAmorphousSpecimen::CreateAmorphous(){
int i, k, kc;
	double xn, yn, zn;

k = 0;
for (i=0; i<na; i++){	
		kc = GeneratePoint(xn, yn, zn);
		if(kc>-1){
			x[k] = xn; y[k] = yn; z[k] = zn;
			c3d[kc] = k;
			k++;
		}
	}
	na = k;
	return na;
}

void cAmorphousSpecimen::Amorphous2Atoms(double *Atoms){
	double lxh = 0.5*lx, lyh = 0.5*ly, lzh = 0.5*lz;
for (int i=0; i<na; i++){
		Atoms[0*na+i]= x[i]-lxh;
		Atoms[1*na+i]= y[i]-lyh;
		Atoms[2*na+i]= z[i]-lzh;
		Atoms[3*na+i]= Z;
		Atoms[5*na+i]= 1.0;
	}
}