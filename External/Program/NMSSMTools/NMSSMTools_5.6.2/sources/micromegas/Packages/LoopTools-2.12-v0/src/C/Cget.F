* Cget.F
* the three-point tensor coefficients
* this file is part of LoopTools
* improvements by M. Rauch
* last modified 23 Mar 15 th

#include "externals.h"
#include "types.h"

#define npoint 3
#include "defs.h"


	subroutine XCpara(para, p1, p2, p1p2, m1, m2, m3)
	implicit none
	DVAR para(1,*)
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	P(1) = p1
	P(2) = p2
	P(3) = p1p2

	M(1) = m1
	if( abs(M(1)) .lt. minmass ) M(1) = 0
	M(2) = m2
	if( abs(M(2)) .lt. minmass ) M(2) = 0
	M(3) = m3
	if( abs(M(3)) .lt. minmass ) M(3) = 0
	end

************************************************************************

	memindex function XCget(p1, p2, p1p2, m1, m2, m3)
	implicit none
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	memindex cacheindex
	external cacheindex, XCcoefx

#ifdef COMPLEXPARA
	memindex Cget
	external Cget
#endif

	DVAR para(1,Pcc)

#ifdef COMPLEXPARA
	if( abs(Im(p1)) + abs(Im(p2)) + abs(Im(p1p2)) .gt. 0 )
     &    print *, "Complex momenta not implemented"
	if( abs(Im(m1)) + abs(Im(m2)) + abs(Im(m3)) .eq. 0 ) then
	  XCget = Cget(p1, p2, p1p2, m1, m2, m3) - offsetC
	  return
	endif
#endif

	call XCpara(para, p1, p2, p1p2, m1, m2, m3)
	XCget = cacheindex(para, Cval(1,0), XCcoefx, RC*Pcc, Ncc, Cno)
	end

************************************************************************

	subroutine XCput(res, p1, p2, p1p2, m1, m2, m3)
	implicit none
	ComplexType res(*)
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	external XCcoefx

	DVAR para(1,Pcc)

#ifdef COMPLEXPARA
	if( abs(Im(p1)) + abs(Im(p2)) + abs(Im(p1p2)) .gt. 0 )
     &    print *, "Complex momenta not implemented"
	if( abs(Im(m1)) + abs(Im(m2)) + abs(Im(m3)) .eq. 0 ) then
	  call Cput(res, p1, p2, p1p2, m1, m2, m3)
	  return
	endif
#endif

	call XCpara(para, p1, p2, p1p2, m1, m2, m3)
	call cachecopy(res, para, Cval(1,0), XCcoefx, RC*Pcc, Ncc, Cno)
	end

************************************************************************

	subroutine XC0nocache(res, p1, p2, p1p2, m1, m2, m3)
	implicit none
	ComplexType res(*)
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	DVAR para(1,Pcc)

#ifdef COMPLEXPARA
	if( abs(Im(p1)) + abs(Im(p2)) + abs(Im(p1p2)) .gt. 0 )
     &    print *, "Complex momenta not implemented"
	if( abs(Im(m1)) + abs(Im(m2)) + abs(Im(m3)) .eq. 0 ) then
	  call C0nocache(res, p1, p2, p1p2, m1, m2, m3)
	  return
	endif
#endif

	call XCpara(para, p1, p2, p1p2, m1, m2, m3)
	call C0func(res, para)
	end

************************************************************************

	ComplexType function XC0i(i, p1, p2, p1p2, m1, m2, m3)
	implicit none
	integer i
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	memindex XCget
	external XCget

	memindex b

	b = XCget(p1, p2, p1p2, m1, m2, m3)
	XC0i = Cval(i+epsi,b)
	end

************************************************************************

	ComplexType function XC0(p1, p2, p1p2, m1, m2, m3)
	implicit none
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	ComplexType XC0i
	external XC0i

	XC0 = XC0i(cc0, p1, p2, p1p2, m1, m2, m3)
	end

************************************************************************

	subroutine XCcoefx(C, para)
	implicit none
	ComplexType C(*)
	DVAR para(1,*)

#include "lt.h"

	memindex b12, b23, b13

	memindex XBget
	external XBget

	b12 = XBget(P(1), M(1), M(2))
	b23 = XBget(P(2), M(2), M(3))
	b13 = XBget(P(3), M(1), M(3))

	call XCcoeff(C, para, Bval(1,b12), Bval(1,b23), Bval(1,b13))
	end

************************************************************************

	subroutine XCcoeff(C, para, B12, B23, B13)
	implicit none
	ComplexType C(*), B12(*), B23(*), B13(*)
	DVAR para(1,*)

#include "lt.h"

	memindex XBget
	external XBget

	DVAR p1, p2, p1p2, m1, m2, m3
	DVAR f1, f2
	QVAR G(2,2)
	ComplexType bsum(0:2), b1sum(0:2), b00sum(0:2), b11sum(0:2)
	ComplexType in(0:2,2)
	logical dump

#ifdef SOLVE_EIGEN
	QVAR Ginv(2,2)
#define IN(i) in(:,i)
#define SOLVE_SETUP XInverse(2, G,2, Ginv,2)
#define SOLVE(b) XSolve(2, G,2, Ginv,2, b)
#else
	integer perm(2)
#define IN(i) in(:,perm(i))
#define SOLVE_SETUP XDecomp(2, G,2, perm)
#define SOLVE(b) XSolve(2, G,2, b)
#endif
#define OUT(i) in(:,i)

	integer finite(0:2)
	data finite /1, 0, 0/

	QVAR aa,aa2,aa3,aa4,aa5,aa6
	double precision pmax,repsT2
	DVAR p(3),m(3)
	QVAR detG,epsT2,alpha(0:2)
	integer i

	serial = serial + 1
	dump = ibits(debugkey, DebugC, 1) .ne. 0 .and.
     &    serial .ge. debugfrom .and. serial .le. debugto
	if( dump ) call XDumpPara(3, para, "Ccoeff")

	m1 = M(1)
	m2 = M(2)
	m3 = M(3)
	p1 = P(1)
	p2 = P(2)
	p1p2 = P(3)

	f1 = m2
	  f1 = f1 - m1
	  f1 = f1 - p1
	f2 = m3
	  f2 = f2 - m1
	  f2 = f2 - p1p2

	G(1,1) = 2*p1
	G(2,2) = 2*p1p2
	G(1,2) = p1
	  G(1,2) = G(1,2) + p1p2
	  G(1,2) = G(1,2) - p2
	G(2,1) = G(1,2)


	detG=G(1,1)*G(2,2)-G(1,2)*G(1,2)
	pmax=Max(abs(p1),abs(p2),abs(p1p2), diffeps)
	if (pmax.eq.diffeps) then
	print*,'warning in Cget.F: All external momenta too small'
c	stop
	endif
	repsT2=abs(detG)/pmax**2/4D0
	if (repsT2 .ge. 1D-8) then


	call SOLVE_SETUP

	bsum = B23(BB0) + B23(BB1)
	b1sum = B23(BB1) + B23(BB11)
	b00sum = B23(BB00) + B23(BB001)
	b11sum = B23(BB11) + B23(BB111)

	call XC0func(C(CC0), para)

	IN(1) = f1*C(CC0) - B23(BB0) + B13(BB0)
	IN(2) = f2*C(CC0) - B23(BB0) + B12(BB0)
	call SOLVE(in)
	C(CC1) = OUT(1)
	C(CC2) = OUT(2)

	C(CC00) = .5D0*(m1*C(CC0) -
     &    .5D0*(f1*C(CC1) + f2*C(CC2) - B23(BB0) - finite))

	IN(1) = f1*C(CC1) + bsum - 2*C(CC00)
	IN(2) = f2*C(CC1) + bsum + B12(BB1)
	call SOLVE(in)
	C(CC11) = OUT(1)
	C(CC12) = OUT(2)

	IN(1) = f1*C(CC2) - B23(BB1) + B13(BB1)
	IN(2) = f2*C(CC2) - B23(BB1) - 2*C(CC00)
	call SOLVE(in)
	C(CC12) = .5D0*(C(CC12) + OUT(1))
	C(CC22) = OUT(2)

	C(CC001) = 1/3D0*(m1*C(CC1) -
     &    .5D0*(f1*C(CC11) + f2*C(CC12) + bsum + finite/3D0))
	C(CC002) = 1/3D0*(m1*C(CC2) -
     &    .5D0*(f1*C(CC12) + f2*C(CC22) - B23(BB1) + finite/3D0))

	bsum = bsum + b1sum
	IN(1) = f1*C(CC11) - bsum - 4*C(CC001)
	IN(2) = f2*C(CC11) - bsum + B12(BB11)
	call SOLVE(in)
	C(CC111) = OUT(1)
	C(CC112) = OUT(2)

	IN(1) = f1*C(CC22) - B23(BB11) + B13(BB11)
	IN(2) = f2*C(CC22) - B23(BB11) - 4*C(CC002)
	call SOLVE(in)
	C(CC122) = OUT(1)
	C(CC222) = OUT(2)

	C(CC0000) = 1/4D0*(m1*C(CC00) -
     &    .5D0*(f1*C(CC001) + f2*C(CC002) - B23(BB00) -
     &      finite*(m1 + m2 + m3 - .25D0*(p1 + p2 + p1p2))/6D0))

	IN(1) = f1*C(CC001) + b00sum - 2*C(CC0000)
	IN(2) = f2*C(CC001) + b00sum + B12(BB001)
	call SOLVE(in)
	C(CC0011) = OUT(1)
	C(CC0012) = OUT(2)

	IN(1) = f1*C(CC002) - B23(BB001) + B13(BB001)
	IN(2) = f2*C(CC002) - B23(BB001) - 2*C(CC0000)
	call SOLVE(in)
	C(CC0012) = .5D0*(C(CC0012) + OUT(1))
	C(CC0022) = OUT(2)

	bsum = bsum + b1sum + b11sum
	IN(1) = f1*C(CC111) + bsum - 6*C(CC0011)
	IN(2) = f2*C(CC111) + bsum + B12(BB111)
	call SOLVE(in)
	C(CC1111) = OUT(1)
	C(CC1112) = OUT(2)

	IN(1) = f1*C(CC222) - B23(BB111) + B13(BB111)
	IN(2) = f2*C(CC222) - B23(BB111) - 6*C(CC0022)
	call SOLVE(in)
	C(CC1222) = OUT(1)
	C(CC2222) = OUT(2)

	IN(1) = f1*C(CC122) + b11sum - 2*C(CC0022)
	IN(2) = f2*C(CC122) + b11sum - 4*C(CC0012)
	call SOLVE(in)
	C(CC1122) = OUT(1)
	C(CC1222) = .5D0*(C(CC1222) + OUT(2))


	else

	do i=1,3
	 p(i)=P(i)
	 m(i)=M(i)
	enddo

	aa = 0d0
        
	if (DBLE(p1) .gt. 1D-12) then
        alpha(2)=1D0
        alpha(1)=(p(2)-p(1)-p(3))/p(1)/2D0
        alpha(0)=-alpha(1)-alpha(2)
        aa=-alpha(0)*m(1)-alpha(1)*m(2)-alpha(2)*m(3)
     >   +alpha(1)*p(1)+alpha(2)*p(3)                           
        epsT2=detG/p(1)/4D0
        if ( abs(aa) .gt. pmax*1D-12 ) goto 1
        endif

        if (DBLE(p1p2) .gt. 1D-12) then
        alpha(1)=1D0
        alpha(2)=(p2-p1p2-p1)/p1p2/2D0
        alpha(0)=-alpha(1)-alpha(2)
        aa=-alpha(0)*m(1)-alpha(1)*m(2)-alpha(2)*m(3)
     >   +alpha(1)*p(1)+alpha(2)*p(3)
        epsT2=detG/p1p2/4D0
        if ( abs(aa) .gt. pmax*1D-12 ) goto 1
        endif

        if (DBLE(p1).le.0 .and. DBLE(p1p2).le.0) then
        print*,'error in Cget.F (segment): All momenta NOT time-like'
        print*,'p1 =',p1
        print*,'p2 =',p2
        print*,'p1p2 =',p1p2
        print*,'m1 =',m1
        print*,'m2 =',m2
        print*,'m3 =',m3
        stop
        endif

        if (abs(aa) .le. 1D-12*pmax ) then
        print*,'Error in Cget.F: aa too small'
        print*,'aa =',aa
        print*,'p1 =',p1
        print*,'p2 =',p2
        print*,'p1p2 =',p1p2
        print*,'m1 =',m1
        print*,'m2 =',m2
        print*,'m3 =',m3
        stop
        endif

 1      Continue

        aa=1D0/aa
        aa2=aa**2
        aa3=aa2*aa
        aa4=aa3*aa
        aa5=aa4*aa
        aa6=aa5*aa

*** rank 0 ***  
       C(cc0) =  + aa * (
     &     + alpha(0)*B23(bb0)
     &     + alpha(1)*B13(bb0)
     &     + alpha(2)*B12(bb0)
     &     )

*** rank 1 ***  
      C(cc1) =  + aa * (
     &     - alpha(0)*B23(bb1)
     &     - alpha(0)*B23(bb0)
     &     + alpha(2)*B12(bb1)
     &     )

      C(cc2) =  + aa * (
     &     + alpha(0)*B23(bb1)
     &     + alpha(1)*B13(bb1)
     &     )

*** rank 2 ***  
      C(cc00) =  + aa * (
     &     + alpha(0)*B23(bb00)
     &     + alpha(1)*B13(bb00)
     &     + alpha(2)*B12(bb00)
     &     )

      C(cc11) =  + aa * (
     &     + 2.D0*alpha(0)*B23(bb1)
     &     + alpha(0)*B23(bb11)
     &     + alpha(0)*B23(bb0)
     &     + alpha(2)*B12(bb11)
     &     )

      C(cc12) =  + aa * (
     &     - alpha(0)*B23(bb1)
     &     - alpha(0)*B23(bb11)
     &     )

      C(cc22) =  + aa * (
     &     + alpha(0)*B23(bb11)
     &     + alpha(1)*B13(bb11)
     &     )

*** rank 3 ***
      C(cc001) =  + aa * (
     &     - alpha(0)*B23(bb00)
     &     - alpha(0)*B23(bb001)
     &     + alpha(2)*B12(bb001)
     &     )

      C(cc111) =  + aa * (
     &     - 3.D0*alpha(0)*B23(bb1)
     &     - 3.D0*alpha(0)*B23(bb11)
     &     - alpha(0)*B23(bb111)
     &     - alpha(0)*B23(bb0)
     &     + alpha(2)*B12(bb111)
     &     )

      C(cc112) =  + aa * (
     &     + alpha(0)*B23(bb1)
     &     + 2.D0*alpha(0)*B23(bb11)
     &     + alpha(0)*B23(bb111)
     &     )

      C(cc122) =  + aa * (
     &     - alpha(0)*B23(bb11)
     &     - alpha(0)*B23(bb111)
     &     )

      C(cc002) =  + aa * (
     &     + alpha(0)*B23(bb001)
     &     + alpha(1)*B13(bb001)
     &     )

      C(cc222) =  + aa * (
     &     + alpha(0)*B23(bb111)
     &     + alpha(1)*B13(bb111)
     &     )

*** rank 4 ***
      C(cc0000) =  + aa * (
     &     + alpha(0)*B23(bb0000)
     &     + alpha(1)*B13(bb0000)
     &     + alpha(2)*B12(bb0000)
     &     )

      C(cc0011) =  + aa * (
     &     + alpha(0)*B23(bb00)
     &     + 2.D0*alpha(0)*B23(bb001)
     &     + alpha(0)*B23(bb0011)
     &     + alpha(2)*B12(bb0011)
     &     )

      C(cc1111) =  + aa * (
     &     + 4.D0*alpha(0)*B23(bb1)
     &     + 6.D0*alpha(0)*B23(bb11)
     &     + 4.D0*alpha(0)*B23(bb111)
     &     + alpha(0)*B23(bb1111)
     &     + alpha(0)*B23(bb0)
     &     + alpha(2)*B12(bb1111)
     &     )

      C(cc1112) =  + aa * (
     &     - alpha(0)*B23(bb1)
     &     - 3.D0*alpha(0)*B23(bb11)
     &     - 3.D0*alpha(0)*B23(bb111)
     &     - alpha(0)*B23(bb1111)
     &     )

      C(cc1122) =  + aa * (
     &     + alpha(0)*B23(bb11)
     &     + 2.D0*alpha(0)*B23(bb111)
     &     + alpha(0)*B23(bb1111)
     &     )

      C(cc0012) =  + aa * (
     &     - alpha(0)*B23(bb001)
     &     - alpha(0)*B23(bb0011)
     &     )

      C(cc1222) =  + aa * (
     &     - alpha(0)*B23(bb111)
     &     - alpha(0)*B23(bb1111)
     &     )

      C(cc0022) =  + aa * (
     &     + alpha(0)*B23(bb0011)
     &     + alpha(1)*B13(bb0011)
     &     )

      C(cc2222) =  + aa * (
     &     + alpha(0)*B23(bb1111)
     &     + alpha(1)*B13(bb1111)
     &     )
	endif


	if( dump ) call XDumpCoeff(3, C)
	end

