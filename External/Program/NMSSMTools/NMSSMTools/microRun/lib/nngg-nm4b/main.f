      real*8 function findValW(name)
      implicit none
      character*(*) name
      findValW=0
      end

      implicit none

      character*100  fInput 
      character*100  fOutput
      character*20 rdBuff
      real *8 v,vcsnngg,vcsnngz,vcsgg,vcsgz
      integer nArgs,err ,slharead,mode
	  real*8 Wh(5), slhaWidth
      integer ModelConstIni, Nch    
      external slhaWidth,ModelConstIni  
      
      OPEN(UNIT=78,FILE='nngg.out',STATUS='UNKNOWN')
      write(78,fmt='("BLOCK lGamma  # AZ and AA cross sections")')
              
      if(slharead('spectr',0).ne.0.or.slharead('decay',1).ne.0) goto 3
      nArgs=iargc()
      if(nArgs.eq.0) then
           mode=0
      else
          mode=1  
      endif
      

      Wh(1)=slhaWidth(25) 
      Wh(2)=slhaWidth(35)
      Wh(3)=slhaWidth(45)
      Wh(4)=slhaWidth(36)
      Wh(5)=slhaWidth(46)
                                                     
	  
	                            
      Nch= ModelConstIni(mode,Wh,err)
      if(err.ne.0) goto 3
      v=0.0
      vcsgg=vcsnngg(v)
      if(Nch.gt.1) then  
        vcsgz=vcsnngz(v)      
      else
        vcsgz=0
      endif    
      write(78,fmt='(A6, 1PE10.4,A20)') '  1   ', vcsgz,'# ~o1,~o1->A,Z [pb]'
      write(78,fmt='(A6, 1PE10.4,A20)') '  2   ', vcsgg,'# ~o1,~o1->A,A [pb]'
      close(78)
      return 
3     close(78)
      return   
      end 
