#ifndef hRandGen_H
#define hRandGen_H

#include "hmathCPU.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>

class cRandGen { 
 private:
	static const unsigned long N = 624;
	static const unsigned long M = 397;
	unsigned long mt[N];
	int mti;	
	unsigned long mag01[2];	
	unsigned long y;
	int j, kk;
	bool bn;
	double rnt;
 public:
	inline cRandGen();
	inline void reset();
 inline void seed(unsigned long s);
	inline int randiu(void);
 inline double randu(void);
	inline double randn(void);
};

inline cRandGen::cRandGen(){
	unsigned int t;
	mag01[0] = 0x0UL;
	mag01[1] = 0x9908b0dfUL;
	mti = N + 1;
	t = (unsigned)time(NULL);
	srand (t);
	seed((unsigned)rand());
	bn = true;
}

inline void cRandGen::reset(){
	unsigned int t;
	mag01[0] = 0x0UL;
	mag01[1] = 0x9908b0dfUL;
	mti = N + 1;
	t = (unsigned)time(NULL);
	srand (t);
	seed((unsigned)rand());
	bn = true;
}

inline void cRandGen::seed(unsigned long s){
	mt[0]= s & 0xffffffffUL;
	for (mti=1; mti<N; mti++){
		mt[mti] = (1812433253UL*(mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
		mt[mti] &= 0xffffffffUL;
	}
	bn = true;
}

inline double cRandGen::randu(void){
	if (mti >= N){ /* generate N words at one time */
		if (mti == N+1) /* if init_genrand() has not been called, */
		seed(5489UL); /* a default initial seed is used */

		for (kk=0;kk<N-M;kk++){
			y = (mt[kk]&0x80000000UL)|(mt[kk+1]&0x7fffffffUL);
			mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (;kk<N-1;kk++){
			y = (mt[kk]&0x80000000UL)|(mt[kk+1]&0x7fffffffUL);
			mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[N-1]&0x80000000UL)|(mt[0]&0x7fffffffUL);
		mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
	}
 
	y = mt[mti++];

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return ((double(y) + 0.5)/4294967296.0); 
}

inline int cRandGen::randiu(void){
	 if (mti >= N){ /* generate N words at one time */
		 if (mti == N+1) /* if init_genrand() has not been called, */
		 seed(5489UL); /* a default initial seed is used */

		 for (kk=0;kk<N-M;kk++){
			 y = (mt[kk]&0x80000000UL)|(mt[kk+1]&0x7fffffffUL);
			 mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		 }
		 for (;kk<N-1;kk++){
			 y = (mt[kk]&0x80000000UL)|(mt[kk+1]&0x7fffffffUL);
			 mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		 }
		 y = (mt[N-1]&0x80000000UL)|(mt[0]&0x7fffffffUL);
		 mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		 mti = 0;
	 }
 
	 y = mt[mti++];

	 /* Tempering */
	 y ^= (y >> 11);
	 y ^= (y << 7) & 0x9d2c5680UL;
	 y ^= (y << 15) & 0xefc60000UL;
	 y ^= (y >> 18);

	return y%2147483647; 
}

inline double cRandGen::randn(void){
	double u, v, s, r;
	if (bn){
		do{
			u = 2.0*randu()-1.0;
			v = 2.0*randu()-1.0;
			s = u*u+v*v;
		}while ((s==0)||(s>=1));
		s = sqrt(-2.0*log(s)/s);
		r = u*s; rnt = v*s;
		bn = false;
	}else{
		r = rnt;
		bn = true;
	}

	return r;
}

#endif