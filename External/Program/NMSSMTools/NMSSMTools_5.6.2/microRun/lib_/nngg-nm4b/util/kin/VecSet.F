* VecSet.F
* construct the CMS representation of a momentum and the
* three corresponding polarization vectors
* this file is part of FormCalc
* last modified 9 Nov 12 th

#include "types.h"


* i: the index of the momentum
* m, p: mass and three-momentum of the particle
* ex,ey,ez: the unit three-vector of the momentum

	subroutine VecSet(i, m, p, ex,ey,ez)
	implicit none
	integer i
	RealType m, p, ex,ey,ez

#include "util.h"

	ComplexType vec_(2,2, nvec, 1)

#define K(x,y)     vec_(x,y, 1, i)
#define EpsTR(x,y) vec_(x,y, 2, i)
#define EpsL(x,y)  vec_(x,y, 3, i)
#define EpsTL(x,y) vec_(x,y, 4, i)

	ComplexType spi_(2, 6:7, nvec, 1)

#define Spi(hel,om,d,x) spi_(x, om, hel+d+5, i)

* i.e., spi(x, 6, 5, i) = left-handed,  undotted, upper
*       spi(x, 7, 5, i) = left-handed,  undotted, lower
*       spi(x, 6, 6, i) = left-handed,  dotted,   upper
*       spi(x, 7, 6, i) = left-handed,  dotted,   lower
*       spi(x, 6, 7, i) = right-handed, undotted, upper
*       spi(x, 7, 7, i) = right-handed, undotted, lower
*       spi(x, 6, 8, i) = right-handed, dotted,   upper
*       spi(x, 7, 8, i) = right-handed, dotted,   lower

	equivalence (momspec, vec_, spi_)

	RealType m2, p2, p0, deltap, sump
	RealType sinth, onePez, oneMez
	ComplexType expIphi

	RealType sqrt2
	parameter (sqrt2 = 1.4142135623730950488016887242096981D0)

c      write(*,*) i,m,p,ex,ey,ez

	momspec(SPEC_M,i) = m
	momspec(SPEC_K,i) = p

	m2 = m**2
	if( m2 .lt. 1D-14 ) then
	  p0 = p
	  deltap = 0
	else
	  p2 = p**2
	  p0 = sqrt(p2 + m2)
	  deltap = m2/(p0 + p)
	endif

	momspec(SPEC_E,i) = p0
	momspec(SPEC_DELTAK,i) = deltap
	momspec(SPEC_EX,i) = ex
	momspec(SPEC_EY,i) = ey
	momspec(SPEC_EZ,i) = ez

	sinth = ex**2 + ey**2
	onePez = 1 + ez
	if( onePez .lt. .5D0 ) onePez = sinth/(1 - ez)
	oneMez = 1 - ez
	if( oneMez .lt. .5D0 ) oneMez = sinth/(1 + ez)

	if( sinth .lt. 1D-14 ) then
* phi is irrelevant when theta = 0
	  expIphi = 1
	  momspec(SPEC_PHI,i) = 0
	else
	  sinth = sqrt(sinth)
	  expIphi = ToComplex(ex, ey)/sinth
	  momspec(SPEC_PHI,i) = atan2(ey, ex)
	endif

	K(1,1) = p0*onePez - deltap*ez
	K(2,2) = p0*oneMez + deltap*ez
	K(2,1) = p*ToComplex(ex, ey)
	K(1,2) = Conjugate(K(2,1))

	if( m .ne. 0 ) then
	  EpsL(1,1) = (p*onePez + deltap*ez)/m
	  EpsL(2,2) = (p*oneMez - deltap*ez)/m
	  EpsL(2,1) = p0/m*ToComplex(ex, ey)
	  EpsL(1,2) = Conjugate(EpsL(2,1))
	endif

	EpsTR(2,2) = 1/sqrt2*sinth
	EpsTL(2,2) = EpsTR(2,2)
	EpsTR(1,1) = -EpsTR(2,2)
	EpsTL(1,1) = EpsTR(1,1)
	EpsTR(2,1) = -1/sqrt2*oneMez*expIphi
	EpsTL(1,2) = Conjugate(EpsTR(2,1))
	EpsTL(2,1) = 1/sqrt2*onePez*expIphi
	EpsTR(1,2) = Conjugate(EpsTL(2,1))

#define safelog(x) log(max(x, 1D-300))

	momspec(SPEC_KT,i) = p*sinth
	momspec(SPEC_ET,i) = sqrt(momspec(SPEC_KT,i)**2 + m2)
	momspec(SPEC_PRAP,i) = safelog(onePez) -
     &    safelog(sinth)
	momspec(SPEC_RAP,i) = safelog(Re(K(1,1))) -
     &    safelog(momspec(SPEC_ET,i))

* this is E^(I phi) cos(th/2) = 1/sqrt2*sqrt(1 + ez)*expIphi:
	expIphi = 1/sqrt2*sqrt(onePez)*expIphi
* this is sin(th/2):
	sinth = 1/sqrt2*sqrt(oneMez)

* 'big' components go as
	sump = sqrt(p0 + p)

	Spi(-1,7,1, 1) = sump*sinth
	Spi(-1,7,2, 1) = Spi(-1,7,1, 1)
	Spi(-1,7,1, 2) = -sump*expIphi
	Spi(-1,7,2, 2) = Conjugate(Spi(-1,7,1, 2))

	Spi(+1,6,2, 1) = sump*expIphi
	Spi(+1,6,1, 1) = Conjugate(Spi(+1,6,2, 1))
	Spi(+1,6,2, 2) = sump*sinth
	Spi(+1,6,1, 2) = Spi(+1,6,2, 2)

* 'small' components go as
	deltap = sqrt(deltap)

	Spi(-1,6,1, 1) = deltap*sinth
	Spi(-1,6,2, 1) = Spi(-1,6,1, 1)
	Spi(-1,6,1, 2) = -deltap*expIphi
	Spi(-1,6,2, 2) = Conjugate(Spi(-1,6,1, 2))

	Spi(+1,7,2, 1) = deltap*expIphi
	Spi(+1,7,1, 1) = Conjugate(Spi(+1,7,2, 1))
	Spi(+1,7,2, 2) = deltap*sinth
	Spi(+1,7,1, 2) = Spi(+1,7,2, 2)
	end

