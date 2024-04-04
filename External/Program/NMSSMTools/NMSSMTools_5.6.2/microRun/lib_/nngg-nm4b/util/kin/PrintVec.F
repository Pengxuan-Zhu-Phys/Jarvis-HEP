* PrintVec.F
* print out the Minkowskian representation of a vector
* this file is part of FormCalc
* last modified 17 Jul 12 th

#include "types.h"


* Note: the arguments a, b specify the vector's position in
* the array vec, NOT the vector itself.  The macros k, e, ec, s
* (defined in util.h) convert from 1...LEGS to a, b.

	subroutine PrintVec(a, s)
	implicit none
	integer a
	character*(*) s

	ComplexType vec(2,2, 0:1)
	common /vectors/ vec

#define A(i,j) vec(i,j,a)

	print *, s, " =", .5D0*(A(1,1) + A(2,2))
	print *, "    ", .5D0*(A(1,2) + A(2,1))
	print *, "    ", (0D0,.5D0)*(A(1,2) - A(2,1))
	print *, "    ", .5D0*(A(1,1) - A(2,2))
	end

