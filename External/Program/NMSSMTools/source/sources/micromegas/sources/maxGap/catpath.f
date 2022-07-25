      integer function strlen(st)

*     Logical length of a string (omitting blanks)
*
*     | F | O | R | T | R | A | N |  |  |  |  |  |  |  |
*     <------------ full length returned by len() ----->
*     <---- logical length ------> <- trailing blank -->
*
      implicit none
      character      st*(*)
      strlen = len(st)
      do while (st(strlen:strlen) .eq. ' ')
       strlen = strlen - 1
      enddo
      return
      end



      character*1000 function catpath(st1,st2)

*     st1 is a string containing the path.
*     st2 is the file name or an other directory name
*     return value is the string "st1/st2"

      implicit none
      character*(*) st1
      character*(*) st2
      integer strlen
      catpath = st1(1:strlen(st1)) // '/' // st2
      return
      end
