* fortranflush.F
* C-callable function to flush the Fortran stdout (unit 6)
* this file is part of FormCalc
* last modified 3 Dec 08 th


	subroutine fortranflush()
	implicit none

	call flush(6)
	end

