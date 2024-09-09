* Split.F
* splits a parent momentum into two child momenta
* loosely based on code by T. Fritzsche
* this file is part of FormCalc
* last modified 30 Nov 11 th

#include "types.h"


* The first child momentum is stored in vector i, the second
* replaces the input variables such that repeated invocation
* builds up the entire phase-space.

* on entry:
*   minv = available invariant mass
*   msum = sum of particle masses 1..n
*   p (ex,ey,ez) = parent momentum
*   xmrem = the fraction of invariant mass that goes into
*     child particles 2..n
*   xcosth and xphi, the fractions of costh and phi which
*     specify the direction relative to (ex,ey,ez) in the CMS
*     frame of the child particles

* on exit:
*   vector i is assigned child particle 1
*   minv = invariant mass remaining for child particles 2..n
*   msum = sum of particle masses 2..n
*   p (ex,ey,ez) = sum of momenta of child particles 2..n


* NB: the following enumeration must match the one in xsection.h
#define Var(v) v(1)
#define Show(v) v(2)
#define Lower(v) v(3)
#define Upper(v) v(4)
#define Step(v) v(5)
#define CutMin(v) v(6)
#define CutMax(v) v(7)

#define Test(c, m) iand(key, (m)*(c)) .ne. 0
#define Zero(c, m) key = iand(key, not((m)*(c)))
#define MomQ(m) ibits(key, m-1, 1)

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

#define CUT_ANY CUT_MIN + CUT_MAX
#define CUT_COSTH_INDIRECT CUT_COSTHCMS + CUT_COSTH_E
#define CUT_MREM_INDIRECT CUT_MREM_E + CUT_MREM_ET + CUT_MREM_RAP

	subroutine Split(i, m, p, ex, ey, ez, minv, msum,
     &    fac, cutkey, xmrem, xcosth, xphi)
	implicit none
	integer i, cutkey
	RealType m, p, ex, ey, ez, minv, msum
	RealType fac, xmrem(*), xcosth(*), xphi(*)

	RealType costh, sinth, cosphi, sinphi
	RealType bmrem(2), bcosth(2)
	RealType fxz, fyz, fxy
	RealType fx, fy, fz
	RealType nx, ny, nz
	RealType ref, cut, c0, c1, c2, a, b
	RealType r, r0, r0costh, rz, r1
	RealType m2, mrem, k, k10, q, qe, qn
	integer key

	RealType rpart(0:1)
	equivalence (rpart(0), r0), (rpart(1), r)

	RealType pi
	parameter (pi = 3.1415926535897932384626433832795029D0)

	RealType ksq, mremsq
	ksq() = (2*c1*(sqrt(c0*(c0 + c2*m2) + (c1*m2)**2) - c1*m2) -
     &    c0*c2)/((2*c1 - c2)*(2*c1 + c2))
	mremsq() = m2 + minv*(minv - 2*sqrt(ksq() + m2))

	key = cutkey

* simultaneous indirect cuts on mrem and costh are not allowed:
	if( Test(CUT_COSTH_INDIRECT, CUT_ANY) .and.
     &      Test(CUT_MREM_INDIRECT, CUT_ANY) )
     &    stop "Illegal key"

* rotate the unit vector e = (ex,ey,ez) into the desired direction:
* rotate by theta in the e-n-plane and then by phi around e.

	if( Step(xphi) .eq. 0 ) fac = 2*pi*fac
	Show(xphi) = 2*pi*Var(xphi)
	cosphi = cos(Show(xphi))
	sinphi = sin(Show(xphi))

	fx = ex
	fy = ey
	fz = ez
	fxz = fx**2 + fz**2
	fyz = fy**2 + fz**2
	if( fyz .lt. fxz ) then
	  fxz = sqrt(fxz)
	  nx = (-cosphi*fx*fy - sinphi*fz)/fxz
	  ny = cosphi*fxz
	  nz = (-cosphi*fy*fz + sinphi*fx)/fxz
	else
	  fyz = sqrt(fyz)
	  nx = cosphi*fyz
	  ny = (-cosphi*fx*fy + sinphi*fz)/fyz
	  nz = (-cosphi*fx*fz - sinphi*fy)/fyz
	endif

	r = p/minv
	r0 = sqrt(r**2 + 1)
	m2 = m**2

	msum = msum - m
	bmrem(1) = msum
	bmrem(2) = minv - m

* direct cut on mrem
	if( Test(CUT_MREM, CUT_MIN) ) then
	  bmrem(1) = max(bmrem(1), CutMin(xmrem))
	  Zero(CUT_MREM, CUT_MIN)
	endif
	if( Test(CUT_MREM, CUT_MAX) ) then
	  bmrem(2) = min(bmrem(2), CutMax(xmrem))
	  Zero(CUT_MREM, CUT_MAX)
	endif

	bcosth(1) = -1
	bcosth(2) = 1

* direct cut on costh
	if( Test(CUT_COSTH, CUT_MIN) ) then
	  bcosth(1) = max(bcosth(1), CutMin(xcosth))
	  Zero(CUT_COSTH, CUT_MIN)
	endif
	if( Test(CUT_COSTH, CUT_MAX) ) then
	  bcosth(2) = min(bcosth(2), CutMax(xcosth))
	  Zero(CUT_COSTH, CUT_MAX)
	endif

* indirect cuts on costh:
* - costh(CMS) cut
* - energy cut
* - momentum cut

	if( Test(CUT_COSTH_INDIRECT, CUT_ANY) ) then
	  mrem = bmrem(1) + Var(xmrem)*(bmrem(2) - bmrem(1))
	  k10 = .5D0*(minv - (mrem - m)*(mrem + m)/minv)
	  k = sqrt((k10 - m)*(k10 + m))

	  if( Test(CUT_COSTHCMS, CUT_MIN) ) then
	    if( k .ne. 0 ) then
	      costh = CutMin(xcosth)
	      sinth = (1 - costh)*(1 + costh)
	      if( sinth .ge. 0 ) then
	        sinth = sqrt(sinth)
	        cut = max((k - m*r*sinth)*(k + m*r*sinth), 0D0)
	        cut = (costh*sqrt(cut) - k10*r0*r*sinth**2)/
     &            (k*((r*sinth)**2 + 1))
	        bcosth(1) = max(bcosth(1), cut)
	      endif
	    endif
	    Zero(CUT_COSTHCMS, CUT_MIN)
	  else if( Test(CUT_COSTH_E, CUT_MIN) ) then
	    if( k*r .ne. 0 ) then
	      cut = CutMin(xcosth)
	      if( MomQ(CUT_MIN) .ne. 0 ) cut = sqrt(cut**2 + m2)
	      bcosth(1) = max(bcosth(1), (cut - k10*r0)/(k*r))
	    endif
	    Zero(CUT_COSTH_K, CUT_MIN)
	  endif

	  if( Test(CUT_COSTHCMS, CUT_MAX) ) then
	    if( k .ne. 0 ) then
	      costh = CutMin(xcosth)
	      sinth = (1 - costh)*(1 + costh)
	      if( sinth .ge. 0 ) then
	        sinth = sqrt(sinth)
	        cut = max((k - m*r*sinth)*(k + m*r*sinth), 0D0)
	        cut = (costh*sqrt(cut) - k10*r0*r*sinth**2)/
     &            (k*((r*sinth)**2 + 1))
	        bcosth(2) = min(bcosth(2), cut)
	      endif
	    endif
	    Zero(CUT_COSTHCMS, CUT_MAX)
	  else if( Test(CUT_COSTH_E, CUT_MAX) ) then
	    if( k*r .ne. 0 ) then
	      cut = CutMax(xcosth)
	      if( MomQ(CUT_MAX) .ne. 0 ) cut = sqrt(cut**2 + m2)
	      bcosth(2) = min(bcosth(2), (cut - k10*r0)/(k*r))
	    endif
	    Zero(CUT_COSTH_K, CUT_MAX)
	  endif
	endif

	bcosth(2) = bcosth(2) - bcosth(1)
	if( bcosth(2) .lt. 0 ) then
	  fac = 0
	  return
	endif
	if( Step(xcosth) .eq. 0 ) fac = fac*bcosth(2)
	costh = bcosth(1) + Var(xcosth)*bcosth(2)
	Show(xcosth) = costh

	r0costh = r0*costh
	sinth = sqrt((1 - costh)*(1 + costh))

* indirect cuts on mrem:
* - energy cut
* - momentum cut
* - transverse-energy cut
* - transverse-momentum cut
* - rapidity cut
* - pseudo-rapidity cut

	if( Test(CUT_MREM_INDIRECT, CUT_ANY) ) then
	  rz = r*fz
	  r1 = sqrt(r**2*(1 + costh**2) + 1)

	  if( Test(CUT_MREM_E, CUT_MIN) ) then
	    c1 = r*r0costh
	    c2 = r1**2
	    ref = m*rpart(MomQ(CUT_MIN))
	    if( CutMin(xmrem) .gt. ref ) then
	      c0 = (CutMin(xmrem) - ref)*(CutMin(xmrem) + ref)
	      cut = mremsq()
	      if( cut .gt. 0 ) bmrem(2) = min(bmrem(2), sqrt(cut))
	    endif
	    Zero(CUT_MREM_K, CUT_MIN)
	  else if( Test(CUT_MREM_ET, CUT_MIN) ) then
	    fxy = fx**2 + fy**2
	    c1 = r*(r0costh*fxy - sinth*fz*nz)
	    a = sqrt((r1 - rz)*(r1 + rz))
	    b = r0costh*fz + sinth*nz
	    c2 = (a - b)*(a + b)
	    ref = m*sqrt(r**2*fxy + (1 - MomQ(CUT_MIN)))
	    if( CutMin(xmrem) .gt. ref ) then
	      c0 = (CutMin(xmrem) - ref)*(CutMin(xmrem) + ref)
	      cut = mremsq()
	      if( cut .gt. 0 ) bmrem(2) = min(bmrem(2), sqrt(cut))
	    endif
	    Zero(CUT_MREM_KT, CUT_MIN)
	  else if( Test(CUT_MREM_RAP, CUT_MIN) ) then
	    cut = tanh(CutMin(xmrem))
	    ref = rpart(MomQ(CUT_MIN))*cut
	    if( ref .gt. rz ) then
	      c0 = m2*(ref - rz)*(ref + rz)
	      c1 = r*(sinth*fz*nz - r0costh*(cut - fz)*(cut + fz))
	      c2 = (r0costh*fz + sinth*nz)**2 -
     &          (cut*r1 - rz)*(cut*r1 + rz)
	      cut = mremsq()
	      if( cut .gt. 0 ) then
	        cut = sqrt(cut)
	        if( c2 .gt. 0 ) then
	          bmrem(1) = max(bmrem(1), cut)
	        else
	          bmrem(2) = min(bmrem(2), cut)
	        endif
	      endif
	    endif
	    Zero(CUT_MREM_PRAP, CUT_MIN)
	  endif

	  if( Test(CUT_MREM_E, CUT_MAX) ) then
	    c1 = r*r0costh
	    c2 = r1**2
	    ref = m*rpart(MomQ(CUT_MAX))
	    if( CutMax(xmrem) .gt. ref ) then
	      c0 = (CutMax(xmrem) - ref)*(CutMax(xmrem) + ref)
	      cut = mremsq()
	      if( cut .gt. 0 ) bmrem(1) = max(bmrem(1), sqrt(cut))
	    endif
	    Zero(CUT_MREM_K, CUT_MAX)
	  else if( Test(CUT_MREM_ET, CUT_MAX) ) then
	    fxy = fx**2 + fy**2
	    c1 = r*(r0costh*fxy - sinth*fz*nz)
	    a = sqrt((r1 - rz)*(r1 + rz))
	    b = r0costh*fz + sinth*nz
	    c2 = (a - b)*(a + b)
	    ref = m*sqrt(r**2*fxy + (1 - MomQ(CUT_MAX)))
	    if( CutMax(xmrem) .gt. ref ) then
	      c0 = (CutMax(xmrem) - ref)*(CutMax(xmrem) + ref)
	      cut = mremsq()
	      if( cut .gt. 0 ) bmrem(1) = max(bmrem(1), sqrt(cut))
	    endif
	    Zero(CUT_MREM_KT, CUT_MAX)
	  else if( Test(CUT_MREM_RAP, CUT_MAX) ) then
	    cut = tanh(CutMax(xmrem))
	    ref = rpart(MomQ(CUT_MAX))*cut
	    if( ref .gt. rz ) then
	      c0 = m2*(ref - rz)*(ref + rz)
	      c1 = r*(sinth*fz*nz - r0costh*(cut - fz)*(cut + fz))
	      c2 = (r0costh*fz + sinth*nz)**2 -
     &          (cut*r1 - rz)*(cut*r1 + rz)
	      cut = mremsq()
	      if( cut .gt. 0 ) then
	        cut = sqrt(cut)
	        if( c2 .gt. 0 ) then
	          bmrem(2) = min(bmrem(2), cut)
	        else
	          bmrem(1) = max(bmrem(1), cut)
	        endif
	      endif
	    endif
	    Zero(CUT_MREM_PRAP, CUT_MAX)
	  endif
	endif

	bmrem(2) = bmrem(2) - bmrem(1)
	if( bmrem(2) .lt. 0 ) then
	  fac = 0
	  return
	endif
	if( Step(xmrem) .eq. 0 ) fac = fac*bmrem(2)
	mrem = bmrem(1) + Var(xmrem)*bmrem(2)
	Show(xmrem) = mrem

	k10 = .5D0*(minv - (mrem - m)*(mrem + m)/minv)
	k = sqrt((k10 - m)*(k10 + m))
	fac = .5D0*k*fac

* k is in the CMS of the decay minv -> m + mrem and must finally
* be boosted into the CMS of the entire process

	qe = k10*r + k*r0costh
	qn = k*sinth
	q = sqrt(qe**2 + qn**2)
	if( q .ne. 0 ) then
	  fx = (qe*fx + qn*nx)/q
	  fy = (qe*fy + qn*ny)/q
	  fz = (qe*fz + qn*nz)/q
	endif
	call VecSet(i, m, q, fx, fy, fz)

	fx = p*ex - q*fx
	fy = p*ey - q*fy
	fz = p*ez - q*fz
	p = sqrt(fx**2 + fy**2 + fz**2)
	if( p .ne. 0 ) then
	  ex = fx/p
	  ey = fy/p
	  ez = fz/p
	endif

	minv = mrem
	end

