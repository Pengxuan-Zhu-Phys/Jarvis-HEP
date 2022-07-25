* util.h
* prototypes for the util functions
* this file is part of FormCalc
* last modified 20 Feb 13 th

#include "types.h"

#ifndef LEGS
#define LEGS 1
#endif

	integer nvec
	parameter (nvec = 10)

* vec(*,*,0) is q1 in num.h
	ComplexType vec(2,2,0:nvec*LEGS), vec_end(0)
	common /vectors/ vec, vec_end

	RealType momspec(8*nvec,LEGS)
	equivalence (vec(1,1,1), momspec)

* encoding base for momenta
	integer*8 JK
	parameter (JK = 256)


#ifndef SPEC_M

#define SPEC_M 65
#define SPEC_K 66
#define SPEC_E 67
#define SPEC_KT 68
#define SPEC_ET 69
#define SPEC_PRAP 70
#define SPEC_RAP 71
#define SPEC_DELTAK 72
#define SPEC_PHI 73
#define SPEC_EX 74
#define SPEC_EY 75
#define SPEC_EZ 76

#define k(i) (nvec*(i-1)+1)
#define s(i) (nvec*(i-1)+3)
#define e(i) (nvec*(i-1)+3+Hel(i))
#define ec(i) (nvec*(i-1)+3-Hel(i))
#define Spinor(i,s,d) (s*Hel(i)+nvec*(i-1)+d+5)

#define MomEncoding(f,i) iand(f,JK-1)*JK**(i-1)

#define signbit(i) ibits(i,31,1)
#define IndexDelta(i,j) signbit(ieor(i,j)-1)
#define Digit(i) char(i+48)
#define Polar(r,theta) r*exp(cI*degree*theta)

#define Error(err,msg) call m_(err, __LINE__, __FILE__, msg)
#define Warning(msg) call m_(0, 0, __FILE__, msg)
#define INFO print *,
#define DEB(a,x) print *, a, x
#define LOOP(var,from,to,step) do var = from, to, step
#define ENDLOOP(var) enddo
#define TEST(i,b) if( btest(i,b) ) then
#define ENDTEST(i,b) endif

#define BIT_RESET 0
#define BIT_LOOP 1
#define BIT_HEL(i) (5*(LEGS-i)+Hel(i)+2)
#define LOOP_HEL(h) do h = -2, 2
#define ENDLOOP_HEL(h) enddo

#define INI_S(seq) call clearcache
#define INI_ANGLE(seq) call markcache
#define DEINI(seq) call restorecache

#ifdef PARALLEL
#define PREP(h,he, v,ve, a,ae, s,se) call sqmeprep(h,he, v,ve, a,ae, s,se)
#define EXEC(f, res, flags) call sqmeexec(f, res, flags)
#define SYNC(res) call sqmesync(res)
#else
#define PREP(h,he, v,ve, a,ae, s,se)
#define EXEC(f, res, flags) call f(res, flags)
#define SYNC(res)
#endif

#define Cut(c,m) (m)*(c)

#define CUT_MIN 1
#define CUT_MAX 2

#define CUT_COSTH 4
#define CUT_COSTHCMS 16
#define CUT_COSTH_E 64
#define CUT_COSTH_K 65
#define CUT_MREM 256
#define CUT_MREM_E 1024
#define CUT_MREM_K 1025
#define CUT_MREM_ET 4096
#define CUT_MREM_KT 4097
#define CUT_MREM_RAP 16384
#define CUT_MREM_PRAP 16385

#endif

