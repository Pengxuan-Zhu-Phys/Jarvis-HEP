      subroutine setCh(n)
      include 'dswacom.h'

      integer n,k
      do  k=1, 29 
        wabr(k)=0
        wabr(n)=1        
      enddo
      return
      end

      subroutine printCh
      include 'dswacom.h'
      integer n,k
      if(wabr(1).gt.0) write(*,*) "  S1 S1            ",wabr(1)      
      if(wabr(2).gt.0) write(*,*) "  S1 S2            ",wabr(2)
      if(wabr(3).gt.0) write(*,*) "  S2 S2            ",wabr(3)
      if(wabr(4).gt.0) write(*,*) "  S3 S3            ",wabr(4)
      if(wabr(5).gt.0) write(*,*) "  S1 S3            ",wabr(5)
      if(wabr(6).gt.0) write(*,*) "  S2 S3            ",wabr(6)
      if(wabr(7).gt.0) write(*,*) "  S- S+            ",wabr(7)
      if(wabr(8).gt.0) write(*,*) "  Z S1             ",wabr(8)
      if(wabr(9).gt.0) write(*,*) "  Z S2             ",wabr(9)
      if(wabr(10).gt.0) write(*,*) "  Z S3             ",wabr(10)
      if(wabr(11).gt.0) write(*,*) "  W- S+ and W+ S-  ",wabr(11)
      if(wabr(12).gt.0) write(*,*) "  Z0 Z0            ",wabr(12)
      if(wabr(13).gt.0) write(*,*) "  W+ W-            ",wabr(13)
      if(wabr(14).gt.0) write(*,*) "  nu_e nu_e-bar    ",wabr(14)
      if(wabr(15).gt.0) write(*,*) "  e+ e-            ",wabr(15)
      if(wabr(16).gt.0) write(*,*) "  nu_mu nu_mu-bar  ",wabr(16)
      if(wabr(17).gt.0) write(*,*) "  mu+ mu-          ",wabr(17)
      if(wabr(18).gt.0) write(*,*) "  nu_tau nu_tau-bar",wabr(18)
      if(wabr(19).gt.0) write(*,*) "  tau+ tau-        ",wabr(19)
      if(wabr(20).gt.0) write(*,*) "  u u-bar          ",wabr(20)
      if(wabr(21).gt.0) write(*,*) "  d d-bar          ",wabr(21)
      if(wabr(22).gt.0) write(*,*) "  c c-bar          ",wabr(22)
      if(wabr(23).gt.0) write(*,*) "  s s-bar          ",wabr(23)
      if(wabr(24).gt.0) write(*,*) "  t t-bar          ",wabr(24)
      if(wabr(25).gt.0) write(*,*) "  b b-bar          ",wabr(25)
      if(wabr(26).gt.0) write(*,*) "  gluon gluon      ",wabr(26)
      if(wabr(27).gt.0) write(*,*) "  q q gluon (not im",wabr(27)
      if(wabr(28).gt.0) write(*,*) "  gamma gamma (1-lo",wabr(28)
      if(wabr(29).gt.0) write(*,*) "  Z0 gamma (1-loop)",wabr(29)
      end
      
      subroutine setMcdmSdSI(Mcdm,csSD,csSI)
      include 'dswacom.h'
      real*4 Mcdm,csSI,cdSD
       write(*,*) "test",  Mcdm,csSD,csSI

       wamwimp=Mcdm
       wasigsip=csSI
       wasigsdp=csSD
      end

      subroutine printMcdmSdSI
      include 'dswacom.h'
      write(*,*) '!!!', wamwimp,wasigsip,wasigsdp
      end
