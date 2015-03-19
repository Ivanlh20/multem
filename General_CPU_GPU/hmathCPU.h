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

#ifndef hmathCPU_H
#define hmathCPU_H

#include "hConstTypes.h"
#include <cmath>

//extern double __cdecl abs (double);
//extern double __cdecl fabs (double);
//extern double __cdecl ceil (double);
//extern double __cdecl floor (double);
//extern double __cdecl sqrt (double);
//extern double __cdecl pow (double, double);
//extern double __cdecl log (double);
//extern double __cdecl log10 (double);
//extern double __cdecl fmod (double, double);
//extern double __cdecl modf (double, double*);
//extern double __cdecl exp (double);
//extern double __cdecl frexp (double, int*);
//extern double __cdecl ldexp (double, int);
//extern double __cdecl sin (double);
//extern double __cdecl cos (double);
//extern double __cdecl tan (double);
//extern double __cdecl sinh (double);
//extern double __cdecl cosh (double);
//extern double __cdecl tanh (double);
//extern double __cdecl asin (double);
//extern double __cdecl acos (double);
//extern double __cdecl atan (double);

inline double rsqrtCPU(double a);

inline double rcbrtCPU(double a);

inline double erfinvCPU(double a);

inline double erfcinvCPU(double a);

inline double erfcxCPU(double a);

inline double log1pCPU(double a);

inline double log2CPU(double a);

inline double exp2CPU(double a);

inline double acoshCPU(double a);

inline double asinhCPU(double a);

inline double atanhCPU(double a);

inline double erfCPU(double a);

inline double erfcCPU(double a);

inline double gammaCPU(double a);

inline double lgammaCPU(double a);

inline double bessjCPU(int n, double x);

inline double bessyCPU(int n, double x);

inline double bessiCPU(int n, double x);

inline double besskCPU(int n, double x);

/*********************************************************************/

inline int __isnanCPU(double a)
{
	volatile union
	{
		double d;
		unsigned long long int l;
	} cvt;

	cvt.d = a;

	return cvt.l << 1 > 0xffe0000000000000ull;
}

inline int __isinfCPU(double a)
{
	volatile union 
	{
		double d;
		unsigned long long int l;
	} cvt;

	cvt.d = a;

	return cvt.l << 1 == 0xffe0000000000000ull;
}

inline int __signbitCPU(double a)
{
	volatile union {
		double d;
		signed long long int l;
	} cvt;

	cvt.d = a;
	return cvt.l < 0ll;
}

inline double copysignCPU(double a, double b)
{
	volatile union
	{
		double d;
		unsigned long long int l;
	} cvta, cvtb;

	cvta.d = a;
	cvtb.d = b;
	cvta.l = (cvta.l & 0x7fffffffffffffffULL) | (cvtb.l & 0x8000000000000000ULL);
	return cvta.d;
}

inline double rintCPU(double a)
{
	double fa = fabs(a);
	double CUDART_TWO_TO_52 = 4503599627370496.0;
	double u = CUDART_TWO_TO_52 + fa;
	if(fa >= CUDART_TWO_TO_52)
	{
	
	}
	else
	{
		u = u - CUDART_TWO_TO_52;
		u = copysignCPU (u, a);
	}
	return u; 
}

inline double rsqrtCPU(double a)
{
	return 1.0/sqrt(a);
}

inline double rcbrtCPU(double a)
{
	double s, t;

	if(__isnanCPU(a))
	{
		return a + a;
	}

	if(a == 0.0 || __isinfCPU(a))
	{
		return 1.0 / a;
	}

	s = fabs(a);
	t = exp2CPU(-3.3333333333333333e-1*log2CPU(s));			/* initial approximation */
	t = ((t*t)*(-s*t)+1.0)*(3.3333333333333333e-1*t) + t;		/* refine approximation */
	if(__signbitCPU(a))
	{
		t = -t;
	}

	return t;
}

inline double erfinvCPU(double a)
{
	double p, q, t, fa;
	volatile union {
		double d;
		unsigned long long int l;
	} cvt;

	fa = fabs(a);
	if(fa >= 1.0)
	{
		cvt.l = 0xfff8000000000000ull;
		t = cvt.d; /* INDEFINITE */
		if(fa == 1.0)
		{
			t = a * exp(1000.0); /* Infinity */
		}
	}
	else if(fa >= 0.9375)
	{
		/* Based on: J.M. Blair, C.A. Edwards, J.H. Johnson: Rational Chebyshev
		Approximations for the Inverse of the Error Function. Mathematics of
		Computation, Vol. 30, No. 136 (Oct. 1976), PotPar. 827-830. Table 59
		*/
		t = log1pCPU(-fa);
		t = 1.0/sqrt(-t);
		p = 2.7834010353747001060e-3;
		p = p * t + 8.6030097526280260580e-1;
		p = p * t + 2.1371214997265515515e+0;
		p = p * t + 3.1598519601132090206e+0;
		p = p * t + 3.5780402569085996758e+0;
		p = p * t + 1.5335297523989890804e+0;
		p = p * t + 3.4839207139657522572e-1;
		p = p * t + 5.3644861147153648366e-2;
		p = p * t + 4.3836709877126095665e-3;
		p = p * t + 1.3858518113496718808e-4;
		p = p * t + 1.1738352509991666680e-6;
		q = t + 2.2859981272422905412e+0;
		q = q * t + 4.3859045256449554654e+0;
		q = q * t + 4.6632960348736635331e+0;
		q = q * t + 3.9846608184671757296e+0;
		q = q * t + 1.6068377709719017609e+0;
		q = q * t + 3.5609087305900265560e-1;
		q = q * t + 5.3963550303200816744e-2;
		q = q * t + 4.3873424022706935023e-3;
		q = q * t + 1.3858762165532246059e-4;
		q = q * t + 1.1738313872397777529e-6;
		t = p / (q * t);
		if(a < 0.0) t = -t;
	}
	else if(fa >= 0.75)
	{
		/* Based on: J.M. Blair, C.A. Edwards, J.H. Johnson: Rational Chebyshev
		Approximations for the Inverse of the Error Function. Mathematics of
		Computation, Vol. 30, No. 136 (Oct. 1976), PotPar. 827-830. Table 39
		*/
		t = a*a - .87890625;
		p = .21489185007307062000e+0;
		p = p * t - .64200071507209448655e+1;
		p = p * t + .29631331505876308123e+2;
		p = p * t - .47644367129787181803e+2;
		p = p * t + .34810057749357500873e+2;
		p = p * t - .12954198980646771502e+2;
		p = p * t + .25349389220714893917e+1;
		p = p * t - .24758242362823355486e+0;
		p = p * t + .94897362808681080020e-2;
		q = t - .12831383833953226499e+2;
		q = q * t + .41409991778428888716e+2;
		q = q * t - .53715373448862143349e+2;
		q = q * t + .33880176779595142685e+2;
		q = q * t - .11315360624238054876e+2;
		q = q * t + .20369295047216351160e+1;
		q = q * t - .18611650627372178511e+0;
		q = q * t + .67544512778850945940e-2;
		p = p / q;
		t = a*p;
	}
	else
	{
		/* Based on: J.M. Blair, C.A. Edwards, J.H. Johnson: Rational Chebyshev
		Approximations for the Inverse of the Error Function. Mathematics of
		Computation, Vol. 30, No. 136 (Oct. 1976), PotPar. 827-830. Table 18
		*/
		t = a * a - .5625;
		p = - .23886240104308755900e+2;
		p = p * t + .45560204272689128170e+3;
		p = p * t - .22977467176607144887e+4;
		p = p * t + .46631433533434331287e+4;
		p = p * t - .43799652308386926161e+4;
		p = p * t + .19007153590528134753e+4;
		p = p * t - .30786872642313695280e+3;
		q = t - .83288327901936570000e+2;
		q = q * t + .92741319160935318800e+3;
		q = q * t - .35088976383877264098e+4;
		q = q * t + .59039348134843665626e+4;
		q = q * t - .48481635430048872102e+4;
		q = q * t + .18997769186453057810e+4;
		q = q * t - .28386514725366621129e+3;
		p = p / q;
		t = a * p;
	}
	return t;
}

inline double erfcinvCPU(double a)
{
	double t;
	volatile union{
		double d;
		unsigned long long int l;
	} cvt;

	if(__isnanCPU(a))
	{
		return a + a;
	}
	if(a <= 0.0)
	{
		cvt.l = 0xfff8000000000000ull;
		t = cvt.d; /* INDEFINITE */
		if(a == 0.0)
		{
			t = (1.0 - a) * exp(1000.0); /* Infinity */
		}
	}
	else if(a >= 0.0625)
	{
		t = erfinvCPU (1.0 - a);
	}
	else if(a >= 1e-100)
	{
		/* Based on: J.M. Blair, C.A. Edwards, J.H. Johnson: Rational Chebyshev
		Approximations for the Inverse of the Error Function. Mathematics of
		Computation, Vol. 30, No. 136 (Oct. 1976), PotPar. 827-830. Table 59
		*/
		double p, q;
		t = log(a);
		t = 1.0 / sqrt(-t);
		p = 2.7834010353747001060e-3;
		p = p * t + 8.6030097526280260580e-1;
		p = p * t + 2.1371214997265515515e+0;
		p = p * t + 3.1598519601132090206e+0;
		p = p * t + 3.5780402569085996758e+0;
		p = p * t + 1.5335297523989890804e+0;
		p = p * t + 3.4839207139657522572e-1;
		p = p * t + 5.3644861147153648366e-2;
		p = p * t + 4.3836709877126095665e-3;
		p = p * t + 1.3858518113496718808e-4;
		p = p * t + 1.1738352509991666680e-6;
		q = t + 2.2859981272422905412e+0;
		q = q * t + 4.3859045256449554654e+0;
		q = q * t + 4.6632960348736635331e+0;
		q = q * t + 3.9846608184671757296e+0;
		q = q * t + 1.6068377709719017609e+0;
		q = q * t + 3.5609087305900265560e-1;
		q = q * t + 5.3963550303200816744e-2;
		q = q * t + 4.3873424022706935023e-3;
		q = q * t + 1.3858762165532246059e-4;
		q = q * t + 1.1738313872397777529e-6;
		t = p / (q * t);
	}
	else 
	{
		/* Based on: J.M. Blair, C.A. Edwards, J.H. Johnson: Rational Chebyshev
		Approximations for the Inverse of the Error Function. Mathematics of
		Computation, Vol. 30, No. 136 (Oct. 1976), PotPar. 827-830. Table 82
		*/
		double p, q;
		t = log(a);
		t = 1.0 / sqrt(-t);
		p = 6.9952990607058154858e-1;
		p = p * t + 1.9507620287580568829e+0;
		p = p * t + 8.2810030904462690216e-1;
		p = p * t + 1.1279046353630280005e-1;
		p = p * t + 6.0537914739162189689e-3;
		p = p * t + 1.3714329569665128933e-4;
		p = p * t + 1.2964481560643197452e-6;
		p = p * t + 4.6156006321345332510e-9;
		p = p * t + 4.5344689563209398450e-12;
		q = t + 1.5771922386662040546e+0;
		q = q * t + 2.1238242087454993542e+0;
		q = q * t + 8.4001814918178042919e-1;
		q = q * t + 1.1311889334355782065e-1;
		q = q * t + 6.0574830550097140404e-3;
		q = q * t + 1.3715891988350205065e-4;
		q = q * t + 1.2964671850944981713e-6;
		q = q * t + 4.6156017600933592558e-9;
		q = q * t + 4.5344687377088206783e-12;
		t = p / (q * t);
	}
	return t;
}

inline double erfcxCPU(double a)
{
	double x, t1, t2, t3;

	if(__isnanCPU(a))
	{
		return a + a;
	}
	x = fabs(a); 
	if(x < 32.0)
	{
		/* 
		* This implementation of erfcxCPU() is based on the algorithm in: M. M. 
		* Shepherd and J. G. Laframboise, "Chebyshev Approximation of (1 + 2x)
		* exp(x^2)erfcCPU x in 0 <= x < INF", Mathematics of Computation, Vol. 
		* 36, No. 153, January 1981, PotPar. 249-253. For the core approximation,
		* the input domain [0,INF] is transformed via (x-k) / (x+k) where k is
		* a precision-dependent constant. Here, we choose k = 4.0, so the input 
		* domain [0, 27.3] is transformed into the core approximation domain 
		* [-1, 0.744409]. 
		*/ 
		/* (1+2*x)*exp(x*x)*erfcCPU(x) */ 
		/* t2 = (x-4.0)/(x+4.0), transforming [0,INF] to [-1,+1] */ 
		t1 = x - 4.0; 
		t2 = x + 4.0; 
		t2 = t1 / t2;
		/* approximate on [-1, 0.744409] */ 
		t1 = - 3.5602694826817400E-010; 
		t1 = t1 * t2 - 9.7239122591447274E-009; 
		t1 = t1 * t2 - 8.9350224851649119E-009; 
		t1 = t1 * t2 + 1.0404430921625484E-007; 
		t1 = t1 * t2 + 5.8806698585341259E-008; 
		t1 = t1 * t2 - 8.2147414929116908E-007; 
		t1 = t1 * t2 + 3.0956409853306241E-007; 
		t1 = t1 * t2 + 5.7087871844325649E-006; 
		t1 = t1 * t2 - 1.1231787437600085E-005; 
		t1 = t1 * t2 - 2.4399558857200190E-005; 
		t1 = t1 * t2 + 1.5062557169571788E-004; 
		t1 = t1 * t2 - 1.9925637684786154E-004; 
		t1 = t1 * t2 - 7.5777429182785833E-004; 
		t1 = t1 * t2 + 5.0319698792599572E-003; 
		t1 = t1 * t2 - 1.6197733895953217E-002; 
		t1 = t1 * t2 + 3.7167515553018733E-002; 
		t1 = t1 * t2 - 6.6330365827532434E-002; 
		t1 = t1 * t2 + 9.3732834997115544E-002; 
		t1 = t1 * t2 - 1.0103906603555676E-001; 
		t1 = t1 * t2 + 6.8097054254735140E-002; 
		t1 = t1 * t2 + 1.5379652102605428E-002; 
		t1 = t1 * t2 - 1.3962111684056291E-001; 
		t1 = t1 * t2 + 1.2329951186255526E+000; 
		/* (1+2*x)*exp(x*x)*erfcCPU(x) / (1+2*x) = exp(x*x)*erfcCPU(x) */ 
		t2 = 2.0 * x + 1.0; 
		t1 = t1 / t2;
	}
	else
	{
		/* asymptotic expansion for large aguments */
		t2 = 1.0 / x;
		t3 = t2 * t2;
		t1 = -29.53125;
		t1 = t1 * t3 + 6.5625;
		t1 = t1 * t3 - 1.875;
		t1 = t1 * t3 + 0.75;
		t1 = t1 * t3 - 0.5;
		t1 = t1 * t3 + 1.0;
		t2 = t2 * 5.6418958354775628e-001;
		t1 = t1 * t2;
	}

	if(a < 0.0)
	{
			/* erfcxCPU(x) = 2*exp(x^2) - erfcxCPU(|x|) */
			t2 = ((int)(x * 16.0)) * 0.0625;
			t3 = (x - t2) * (x + t2);
			t3 = exp(t2 * t2) * exp(t3);
			t3 = t3 + t3;
			t1 = t3 - t1;
	}
	return t1;
}

inline double log1pCPU(double a)
{
	volatile double u, m;

	u = 1.0 + a;
	if(u == 1.0)
	{
		/* a very close to zero */
		u = a;
	}
	else
	{
		m = u - 1.0;
		u = log(u);
		if(a < 1.0)
		{
			/* a somewhat close to zero */
			u = a * u;
			u = u / m;
		}
	}
	return u;
}

inline double log2CPU(double a)
{
	return log(a)/log(2.0);
}

inline double exp2CPU(double a)
{
	return pow(2.0, a);
}

inline double acoshCPU(double a)
{
	double s, t;

	t = a - 1.0;
	if(t == a)
	{
		return log(2.0) + log(a);
	}
	else
	{
		s = a + 1.0;
		t = t + sqrt(s * t);
		return log1pCPU(t);
	}
}

inline double asinhCPU(double a)
{
	double fa, oofa, t;

	fa = fabs(a);
	if(fa > 1e18)
	{
		t = log(2.0) + log(fa);
	}
	else
	{
		oofa = 1.0 / fa;
		t = fa + fa / (oofa + sqrt(1.0 + oofa * oofa));
		t = log1pCPU(t);
	}
	t = copysignCPU(t, a);
	return t;
}

inline double atanhCPU(double a)
{
	double fa, t;

	fa = fabs(a);
	t = (2.0*fa)/(1.0 - fa);
	t = 0.5 * log1pCPU(t);
	if(__isnanCPU(t) || !__signbitCPU(a))
	{
		return t;
	}
		return -t;
}

inline double erfCPU(double a)
{
	double t, r, q;

	t = fabs(a);
	if(t >= 1.0)
	{
		r = -1.28836351230756500E-019;
		r = r * t + 1.30597472161093370E-017;
		r = r * t - 6.33924401259620500E-016;
		r = r * t + 1.96231865908940140E-014;
		r = r * t - 4.35272243559990750E-013;
		r = r * t + 7.37083927929352150E-012;
		r = r * t - 9.91402142550461630E-011;
		r = r * t + 1.08817017167760820E-009;
		r = r * t - 9.93918713097634620E-009;
		r = r * t + 7.66739923255145500E-008;
		r = r * t - 5.05440278302806720E-007;
		r = r * t + 2.87474157099000620E-006;
		r = r * t - 1.42246725399722510E-005;
		r = r * t + 6.16994555079419460E-005;
		r = r * t - 2.36305221938908790E-004;
		r = r * t + 8.05032844055371070E-004;
		r = r * t - 2.45833366629108140E-003;
		r = r * t + 6.78340988296706120E-003;
		r = r * t - 1.70509103597554640E-002;
		r = r * t + 3.93322852515666300E-002;
		r = r * t - 8.37271292613764040E-002;
		r = r * t + 1.64870423707623280E-001;
		r = r * t - 2.99729521787681470E-001;
		r = r * t + 4.99394435612628580E-001;
		r = r * t - 7.52014596480123030E-001;
		r = r * t + 9.99933138314926250E-001;
		r = r * t - 1.12836725321102670E+000;
		r = r * t + 9.99998988715182450E-001;
		q = exp (-t * t);
		r = 1.0 - r * q;
		if(t >= 6.5)
		{
			r = 1.0;
		}
		a = copysignCPU(r, a);
	}
	else
	{
		q = a * a;
		r = -7.77946848895991420E-010;
		r = r * q + 1.37109803980285950E-008;
		r = r * q - 1.62063137584932240E-007;
		r = r * q + 1.64471315712790040E-006;
		r = r * q - 1.49247123020098620E-005;
		r = r * q + 1.20552935769006260E-004;
		r = r * q - 8.54832592931448980E-004;
		r = r * q + 5.22397760611847340E-003;
		r = r * q - 2.68661706431114690E-002;
		r = r * q + 1.12837916709441850E-001;
		r = r * q - 3.76126389031835210E-001;
		r = r * q + 1.12837916709551260E+000;
		a = r * a;
	}
	return a;
}

inline double erfcCPU(double a)
{
	double p, q, h, l;

	if(a < 0.75)
	{
		return 1.0 - erfCPU(a);
	}
	if(a > 27.3)
	{
		return 0.0;
	}
	if(a < 5.0)
	{
		double t;
		t = 1.0/a;
		p = 1.9759923722227928E-008;
		p = p * t - 1.0000002670474897E+000;
		p = p * t - 7.4935303236347828E-001;
		p = p * t - 1.5648136328071860E-001;
		p = p * t + 1.2871196242447239E-001;
		p = p * t + 1.1126459974811195E-001;
		p = p * t + 4.0678642255914332E-002;
		p = p * t + 7.9915414156678296E-003;
		p = p * t + 7.1458332107840234E-004;
		q = t + 2.7493547525030619E+000;
		q = q * t + 3.3984254815725423E+000;
		q = q * t + 2.4635304979947761E+000;
		q = q * t + 1.1405284734691286E+000;
		q = q * t + 3.4130157606195649E-001;
		q = q * t + 6.2250967676044953E-002;
		q = q * t + 5.5661370941268700E-003;
		q = q * t + 1.0575248365468671E-009;
		p = p / q;
		p = p * t;
		h = ((int)(a * 16.0)) * 0.0625;
		l = (a - h) * (a + h);
		q = exp(-h * h) * exp(-l);
		q = q * 0.5;
		p = p * q + q;
		p = p * t;
	}
	else
	{
		double ooa, ooasq;

		ooa = 1.0/a;
		ooasq = ooa * ooa;
		p = -4.0025406686930527E+005;
		p = p * ooasq + 1.4420582543942123E+005;
		p = p * ooasq - 2.7664185780951841E+004;
		p = p * ooasq + 4.1144611644767283E+003;
		p = p * ooasq - 5.8706000519209351E+002;
		p = p * ooasq + 9.1490086446323375E+001;
		p = p * ooasq - 1.6659491387740221E+001;
		p = p * ooasq + 3.7024804085481784E+000;
		p = p * ooasq - 1.0578553994424316E+000;
		p = p * ooasq + 4.2314218745087778E-001;
		p = p * ooasq - 2.8209479177354962E-001;
		p = p * ooasq + 5.6418958354775606E-001;
		h = a * a;
		h = ((int)(a * 16.0)) * 0.0625;
		l = (a - h) * (a + h);
		q = exp(-h * h) * exp(-l);
		p = p * ooa;
		p = p * q;
	}
	return p;
}

inline double gammaCPU(double a)
{
	// Split the function domain into three intervals:
	// (0, 0.001), [0.001, 12), and (12, infinity)

	///
	// First interval: (0, 0.001)
	//
	// For small a, 1/gamma(a) has power series a + gamma a^2 - ...
	// So in this range, 1/gamma(a) = a + gamma a^2 with error on the order of a^3.
	// The relative error over this interval is less than 6e-7.

	const double gamma = 0.577215664901532860606512090; // Euler's gamma constant

	if(a < 0.001)
	{
		return 1.0/(a*(1.0 + gamma*a));
	}

	///
	// Second interval: [0.001, 12)
 
	if(a < 12.0)
	{
		// The algorithm directly approximates gamma over (1,2) and uses
		// reduction identities to reduce other arguments to this interval.
		
		double y = a;
		int n = 0;
		bool arg_was_less_than_one = (y < 1.0);

		// Add or subtract integers as necessary to bring y into (1,2)
		// Will correct for this below
		if(arg_was_less_than_one)
		{
			y += 1.0;
		}
		else
		{
			n = static_cast<int> (floor(y)) - 1; // will use n later
			y -= n;
		}

		// numerator coefficients for approximation over the interval (1,2)
		static const double p[] =
		{
			-1.71618513886549492533811E+0,
				2.47656508055759199108314E+1,
			-3.79804256470945635097577E+2,
				6.29331155312818442661052E+2,
				8.66966202790413211295064E+2,
			-3.14512729688483675254357E+4,
			-3.61444134186911729807069E+4,
				6.64561438202405440627855E+4
		};

		// denominator coefficients for approximation over the interval (1,2)
		static const double q[] =
		{
			-3.08402300119738975254353E+1,
				3.15350626979604161529144E+2,
			-1.01515636749021914166146E+3,
			-3.10777167157231109440444E+3,
				2.25381184209801510330112E+4,
				4.75584627752788110767815E+3,
			-1.34659959864969306392456E+5,
			-1.15132259675553483497211E+5
		};

		double num = 0.0;
		double den = 1.0;
		int i;

		double z = y - 1;
		for(i = 0; i < 8; i++)
		{
			num = (num + p[i])*z;
			den = den*z + q[i];
		}
		double result = num/den + 1.0;

		// Apply correction if argument was not initially in (1,2)
		if(arg_was_less_than_one)
		{
			// Use identity gamma(z) = gamma(z+1)/z
			// The variable "result" now holds gamma of the original y + 1
			// Thus we use y-1 to get back the orginal y.
			result /= (y-1.0);
		}
		else
		{
			// Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
			for(i = 0; i < n; i++)
				result *= y++;
		}

		return result;
	}

	///
	// Third interval: [12, infinity)

	if(a > 171.624)
	{
		// Correct answer too large to display. Force +infinity.
		double temp = 1.7976931348623157e+308;
		return temp*2.0;
	}

	return exp(lgammaCPU(a));
}

inline double lgammaCPU(double a)
{
	double s;
	double t;
	double i;
	double fa;
	double sum;
	long long int quot;
	if(__isnanCPU(a) || __isinfCPU(a))
	{
		return a * a;
	}
	fa = fabs(a);
	if(fa >= 3.0)
	{
		if(fa >= 8.0)
		{
			/* Stirling approximation; coefficients from Hart et al, "Computer 
			* Approximations", Wiley 1968. Approximation 5404. 
			*/
			s = 1.0 / fa;
			t = s * s;
			sum = -0.1633436431e-2;
			sum = sum * t + 0.83645878922e-3;
			sum = sum * t - 0.5951896861197e-3;
			sum = sum * t + 0.793650576493454e-3;
			sum = sum * t - 0.277777777735865004e-2;
			sum = sum * t + 0.833333333333331018375e-1;
			sum = sum * s + 0.918938533204672;
			s = 0.5 * log (fa);
			t = fa - 0.5;
			s = s * t;
			t = s - fa;
			s = s + sum;
			t = t + s;
		}
		else
		{
			i = fa - 3.0;
			s = -4.02412642744125560E+003;
			s = s * i - 2.97693796998962000E+005;
			s = s * i - 6.38367087682528790E+006;
			s = s * i - 5.57807214576539320E+007;
			s = s * i - 2.24585140671479230E+008;
			s = s * i - 4.70690608529125090E+008;
			s = s * i - 7.62587065363263010E+008;
			s = s * i - 9.71405112477113250E+008;
			t = i - 1.02277248359873170E+003;
			t = t * i - 1.34815350617954480E+005;
			t = t * i - 4.64321188814343610E+006;
			t = t * i - 6.48011106025542540E+007;
			t = t * i - 4.19763847787431360E+008;
			t = t * i - 1.25629926018000720E+009;
			t = t * i - 1.40144133846491690E+009;
			t = s / t;
			t = t + i;
		}
	}
	else if(fa >= 1.5)
	{
		i = fa - 2.0;
		t = 9.84839283076310610E-009;
		t = t * i - 6.69743850483466500E-008;
		t = t * i + 2.16565148880011450E-007;
		t = t * i - 4.86170275781575260E-007;
		t = t * i + 9.77962097401114400E-007;
		t = t * i - 2.03041287574791810E-006;
		t = t * i + 4.36119725805364580E-006;
		t = t * i - 9.43829310866446590E-006;
		t = t * i + 2.05106878496644220E-005;
		t = t * i - 4.49271383742108440E-005;
		t = t * i + 9.94570466342226000E-005;
		t = t * i - 2.23154589559238440E-004;
		t = t * i + 5.09669559149637430E-004;
		t = t * i - 1.19275392649162300E-003;
		t = t * i + 2.89051032936815490E-003;
		t = t * i - 7.38555102806811700E-003;
		t = t * i + 2.05808084278121250E-002;
		t = t * i - 6.73523010532073720E-002;
		t = t * i + 3.22467033424113040E-001;
		t = t * i + 4.22784335098467190E-001;
		t = t * i;
	}
	else if(fa >= 0.7)
	{
		i = 1.0 - fa;
		t = 1.17786911519331130E-002; 
		t = t * i + 3.89046747413522300E-002;
		t = t * i + 5.90045711362049900E-002;
		t = t * i + 6.02143305254344420E-002;
		t = t * i + 5.61652708964839180E-002;
		t = t * i + 5.75052755193461370E-002;
		t = t * i + 6.21061973447320710E-002;
		t = t * i + 6.67614724532521880E-002;
		t = t * i + 7.14856037245421020E-002;
		t = t * i + 7.69311251313347100E-002;
		t = t * i + 8.33503129714946310E-002;
		t = t * i + 9.09538288991182800E-002;
		t = t * i + 1.00099591546322310E-001;
		t = t * i + 1.11334278141734510E-001;
		t = t * i + 1.25509666613462880E-001;
		t = t * i + 1.44049896457704160E-001;
		t = t * i + 1.69557177031481600E-001;
		t = t * i + 2.07385551032182120E-001;
		t = t * i + 2.70580808427600350E-001;
		t = t * i + 4.00685634386517050E-001;
		t = t * i + 8.22467033424113540E-001;
		t = t * i + 5.77215664901532870E-001;
		t = t * i;
	}
	else
	{
		t = -9.04051686831357990E-008;
		t = t * fa + 7.06814224969349250E-007;
		t = t * fa - 3.80702154637902830E-007;
		t = t * fa - 2.12880892189316100E-005;
		t = t * fa + 1.29108470307156190E-004;
		t = t * fa - 2.15932815215386580E-004;
		t = t * fa - 1.16484324388538480E-003;
		t = t * fa + 7.21883433044470670E-003;
		t = t * fa - 9.62194579514229560E-003;
		t = t * fa - 4.21977386992884450E-002;
		t = t * fa + 1.66538611813682460E-001;
		t = t * fa - 4.20026350606819980E-002;
		t = t * fa - 6.55878071519427450E-001;
		t = t * fa + 5.77215664901523870E-001;
		t = t * fa;
		t = t * fa + fa;
		t = -log (t);
	}
	if(a >= 0.0)
	{
		return t;
	}
	if(fa < 1e-19)
	{
		return -log(fa);
	}
	i = floor(fa); 
	if(fa == i)
	{
		return 1.0 / (fa - i); /* a is an integer: return infinity */
	}
	i = rintCPU (2.0 * fa);
	quot = (long long int)i;
	i = fa - 0.5 * i;
	i = i * cPi;
	if(quot & 1)
	{
		i = cos(i);
	}
	else
	{
		i = sin(i);
	}
	i = fabs(i);
	t = log(cPi / (i * fa)) - t;
	return t;
}
 
#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order */
/* 0 at input x */
/*------------------------------------------------------------*/
inline double bessj0CPU(double x)
{
	double ax,z;
	double xx,y,ans,ans1,ans2;

	if((ax=fabs(x)) < 8.0)
	{
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
		+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
		+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	}
	else
	{
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
		+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
		+y*(-0.6911147651e-5+y*(0.7621095161e-6
		-y*0.934935152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order */
/* 1 at input x */
/*------------------------------------------------------------*/
inline double bessj1CPU(double x)
{
	double ax,z;
	double xx,y,ans,ans1,ans2;

	if((ax=fabs(x)) < 8.0)
	{
		y=x*x;
		ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
		+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
		+y*(99447.43394+y*(376.9991397+y*1.0))));
		ans=ans1/ans2;
	}
	else
	{
		z=8.0/ax;
		y=z*z;
		xx=ax-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
		+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
		+y*(0.8449199096e-5+y*(-0.88228987e-6
		+y*0.105787412e-6)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		if(x < 0.0)
		{
			ans = -ans;
		}
	}
	return ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order */
/* n at input x */
/* The function can also be called for n = 0 and n = 1. */
/*------------------------------------------------------------*/
inline double bessjCPU(int n, double x)
{
	int j, jsum, m;
	double ax, bj, bjm, bjp, sum, tox, ans;

	//if(n < 0)
	//{
	// double dblank;
	// setdblank_c( &dblank );
	// return( dblank );
	//}
	ax=fabs(x);
	if(n == 0)
	{
		return( bessj0CPU(ax) );
	}
	if(n == 1)
	{
		return( bessj1CPU(ax) );
	}

	if(ax == 0.0)
	{
		return 0.0;
	}
	else if(ax > (double) n)
	{
		tox=2.0/ax;
		bjm=bessj0CPU(ax);
		bj=bessj1CPU(ax);
		for(j=1;j<n;j++)
		{
			bjp=j*tox*bj-bjm;
			bjm=bj;
			bj=bjp;
		}
		ans=bj;
	}
	else
	{
		tox=2.0/ax;
		m=2*((n+(int)sqrt(ACC*n))/2);
		jsum=0;
		bjp=ans=sum=0.0;
		bj=1.0;
		for(j=m;j>0;j--)
		{
			bjm=j*tox*bj-bjp;
			bjp=bj;
			bj=bjm;
			if(fabs(bj) > BIGNO)
			{
				bj *= BIGNI;
				bjp *= BIGNI;
				ans *= BIGNI;
				sum *= BIGNI;
			}
			if(jsum)
			{
				sum += bj;
			}
			jsum=!jsum;
			if(j == n)
			{
				ans=bjp;
			}
		}
		sum=2.0*sum-bj;
		ans /= sum;
	}
	return x < 0.0 && n%2 == 1 ? -ans : ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of second kind and order */
/* 0 at input x. */
/*------------------------------------------------------------*/
inline double bessy0CPU(double x)
{
	double z;
	double xx,y,ans,ans1,ans2;

	if(x < 8.0)
	{
		y=x*x;
		ans1 = -2957821389.0+y*(7062834065.0+y*(-512359803.6
		+y*(10879881.29+y*(-86327.92757+y*228.4622733))));
		ans2=40076544269.0+y*(745249964.8+y*(7189466.438
		+y*(47447.26470+y*(226.1030244+y*1.0))));
		ans= (ans1/ans2)+0.636619772*bessj0CPU(x)*log(x);
	}
	else
	{
		z=8.0/x;
		y=z*z;
		xx=x-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
		+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
		+y*(-0.6911147651e-5+y*(0.7621095161e-6
		+y*(-0.934945152e-7))));
		ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
	}
	return ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of second kind and order */
/* 1 at input x. */
/*------------------------------------------------------------*/
inline double bessy1CPU(double x)
{
	double z;
	double xx,y,ans,ans1,ans2;

	if(x < 8.0)
	{
		y=x*x;
		ans1=x*(-0.4900604943e13+y*(0.1275274390e13
		+y*(-0.5153438139e11+y*(0.7349264551e9
		+y*(-0.4237922726e7+y*0.8511937935e4)))));
		ans2=0.2499580570e14+y*(0.4244419664e12
		+y*(0.3733650367e10+y*(0.2245904002e8
		+y*(0.1020426050e6+y*(0.3549632885e3+y)))));
		ans= (ans1/ans2)+0.636619772*(bessj1CPU(x)*log(x)-1.0/x);
	}
	else
	{
		z=8.0/x;
		y=z*z;
		xx=x-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
		+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
		+y*(0.8449199096e-5+y*(-0.88228987e-6
		+y*0.105787412e-6)));
		ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
	}
	return ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of second kind and order */
/* n for input x. (n >= 0) */
/* Note that for x == 0 the functions bessyCPU and besskCPU are not */
/* defined and a blank is returned. */
/*------------------------------------------------------------*/
inline double bessyCPU(int n, double x)
{
	int j;
	double by,bym,byp,tox;

	//if(n < 0 || x == 0.0)
	//{
	// double dblank;
	// setdblank_c( &dblank );
	// return( dblank );
	//}
	if(n == 0)
	{
		return( bessy0CPU(x) );
	}
	if(n == 1)
	{
		return( bessy1CPU(x) );
	}

	tox=2.0/x;
	by=bessy1CPU(x);
	bym=bessy0CPU(x);
	for(j=1;j<n;j++)
	{
		byp=j*tox*by-bym;
		bym=by;
		by=byp;
	}
	return by;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=0. */
/*------------------------------------------------------------*/
inline double bessi0CPU(double x)
{
	double ax,ans;
	double y;

	if((ax=fabs(x)) < 3.75)
	{
		y=x/3.75,y=y*y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
		+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	}
	else
	{
		y=3.75/ax;
		ans= (exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
		+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
		+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
		+y*0.392377e-2))))))));
	}
	return ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=1. */
/*------------------------------------------------------------*/
inline double bessi1CPU(double x)
{
	double ax,ans;
	double y;

	if((ax=fabs(x)) < 3.75)
	{
		y=x/3.75,y=y*y;
		ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
		+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
	}
	else
	{
		y=3.75/ax;
		ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
		-y*0.420059e-2));
		ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
		+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
		ans *= (exp(ax)/sqrt(ax));
	}
	return x < 0.0 ? -ans : ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) for n >= 0*/
/*------------------------------------------------------------*/
inline double bessiCPU(int n, double x)
{
	int j;
	double bi,bim,bip,tox,ans;

	//if(n < 0)
	//{
	// double dblank;
	// setdblank_c( &dblank );
	// return( dblank );
	//}
	if(n == 0)
	{
		return( bessi0CPU(x) );
	}
	if(n == 1)
	{
		return( bessi1CPU(x) );
	}

	if(x == 0.0)
	{
		return 0.0;
	}
	else
	{
		tox=2.0/fabs(x);
		bip=ans=0.0;
		bi=1.0;
		for(j=2*(n+(int) sqrt(ACC*n));j>0;j--)
		{
			bim=bip+j*tox*bi;
			bip=bi;
			bi=bim;
			if(fabs(bi) > BIGNO)
			{
				ans *= BIGNI;
				bi *= BIGNI;
				bip *= BIGNI;
			}
			if(j == n)
			{
				ans=bip;
			}
		}
		ans *= bessi0CPU(x)/bi;
		return x < 0.0 && n%2 == 1 ? -ans : ans;
	}
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=0. */
/*------------------------------------------------------------*/
inline double bessk0CPU(double x)
{
	double y,ans;

	if(x <= 2.0)
	{
		y=x*x/4.0;
		ans= (-log(x/2.0)*bessi0CPU(x))+(-0.57721566+y*(0.42278420
		+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
		+y*(0.10750e-3+y*0.74e-5))))));
	}
	else
	{
		y=2.0/x;
		ans= (exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
		+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
		+y*(-0.251540e-2+y*0.53208e-3))))));
	}
	return ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=1. */
/*------------------------------------------------------------*/
inline double bessk1CPU(double x)
{
	double y,ans;

	if(x <= 2.0)
	{
		y=x*x/4.0;
		ans= (log(x/2.0)*bessi1CPU(x))+(1.0/x)*(1.0+y*(0.15443144
		+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
		+y*(-0.110404e-2+y*(-0.4686e-4)))))));
	}
	else
	{
		y=2.0/x;
		ans= (exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
		+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
		+y*(0.325614e-2+y*(-0.68245e-3)))))));
	}
	return ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n >= 0*/
/* Note that for x == 0 the functions bessyCPU and besskCPU are not */
/* defined and a blank is returned. */
/*------------------------------------------------------------*/
inline double besskCPU(int n, double x)
{
	int j;
	double bk,bkm,bkp,tox;

	//if(n < 0 || x == 0.0)
	//{
	// double dblank;
	// setdblank_c( &dblank );
	// return( dblank );
	//}
	if(n == 0)
	{
		return( bessk0CPU(x) );
	}
	if(n == 1)
	{
		return( bessk1CPU(x) );
	}

	tox=2.0/x;
	bkm=bessk0CPU(x);
	bk=bessk1CPU(x);
	for(j=1;j<n;j++)
	{
		bkp=bkm+j*tox*bk;
		bkm=bk;
		bk=bkp;
	}
	return bk;
}

#undef ACC
#undef BIGNO
#undef BIGNI

#endif