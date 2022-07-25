      SUBROUTINE MNEU_CPV(IFAIL)

c      Neutralino masses: diag(MNEU(I)) = NEU*.MNEU.NEU+
c      The squared neutralino masses MNEU(i) (i=1..5) as well as the
c      rotation matrix NEU(i,j,k) (i,j=1..5; k=1..2 for real/imaginary part)
c      are stored in the common NEUSPEC.
c      In case of a negative MNEU(i), IFAIL=7

      IMPLICIT NONE

      INTEGER I,J,IFAIL

      DOUBLE PRECISION MNEU2(10,10),VALP(10),VECP(10,10),MNEUP(5,2)
      DOUBLE PRECISION aux,aux1,aux2

      DOUBLE PRECISION l,k,Alcos1,Akcos2,muq,nuq
      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION mur,M1r,M2r,msi
      DOUBLE PRECISION MNEU(5),NEU(5,5,2)
      DOUBLE PRECISION phiF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      DOUBLE PRECISION DDCOS,DDSIN,Pi

      COMMON/QPAR/l,k,Alcos1,Akcos2,muq,NUQ

      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/GAUGINOPAR/mur,M1r,M2r,msi
      COMMON/NEUSPEC/MNEU,NEU
      COMMON/Z3VAUX/phiF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si


      PI=4d0*DATAN(1d0)

c       1) Squared mass matrix in the format [[Re,Im],[-Im,Re]]
      MNEU2(1,1)=M1r**2+g1/2d0*(vu**2+vd**2)
      MNEU2(6,6)=M1r**2+g1/2d0*(vu**2+vd**2)
      MNEU2(6,1)=0d0
      MNEU2(1,6)=0d0

      MNEU2(1,2)=-dsqrt(g1*g2)/2d0*(vu**2+vd**2)
      MNEU2(2,1)=-dsqrt(g1*g2)/2d0*(vu**2+vd**2)
      MNEU2(6,7)=-dsqrt(g1*g2)/2d0*(vu**2+vd**2)
      MNEU2(7,6)=-dsqrt(g1*g2)/2d0*(vu**2+vd**2)
      MNEU2(1,7)=0d0
      MNEU2(2,6)=0d0
      MNEU2(7,1)=0d0
      MNEU2(6,2)=0d0
      
      MNEU2(1,3)=dsqrt(g1/2d0)*(M1r*vu*DDCOS(PhiM1)
     .                                 +mur*vd*DDCOS(Phi01))
      MNEU2(3,1)=dsqrt(g1/2d0)*(M1r*vu*DDCOS(PhiM1)
     .                                 +mur*vd*DDCOS(Phi01))
      MNEU2(6,8)=dsqrt(g1/2d0)*(M1r*vu*DDCOS(PhiM1)
     .                                 +mur*vd*DDCOS(Phi01))
      MNEU2(8,6)=dsqrt(g1/2d0)*(M1r*vu*DDCOS(PhiM1)
     .                                 +mur*vd*DDCOS(Phi01))
      MNEU2(1,8)=dsqrt(g1/2d0)*(-M1r*vu*DDSIN(PhiM1)
     .                                 +mur*vd*DDSIN(Phi01))
      MNEU2(3,6)=dsqrt(g1/2d0)*(M1r*vu*DDSIN(PhiM1)
     .                                 -mur*vd*DDSIN(Phi01))
      MNEU2(8,1)=dsqrt(g1/2d0)*(-M1r*vu*DDSIN(PhiM1)
     .                                 +mur*vd*DDSIN(Phi01))
      MNEU2(6,3)=dsqrt(g1/2d0)*(M1r*vu*DDSIN(PhiM1)
     .                                 -mur*vd*DDSIN(Phi01))

      MNEU2(1,4)=-dsqrt(g1/2d0)*(M1r*vd*DDCOS(PhiM1)
     .                                 +mur*vu*DDCOS(Phi01))
      MNEU2(4,1)=-dsqrt(g1/2d0)*(M1r*vd*DDCOS(PhiM1)
     .                                 +mur*vu*DDCOS(Phi01))
      MNEU2(6,9)=-dsqrt(g1/2d0)*(M1r*vd*DDCOS(PhiM1)
     .                                 +mur*vu*DDCOS(Phi01))
      MNEU2(9,6)=-dsqrt(g1/2d0)*(M1r*vd*DDCOS(PhiM1)
     .                                 +mur*vu*DDCOS(Phi01))
      MNEU2(1,9)=dsqrt(g1/2d0)*(M1r*vd*DDSIN(PhiM1)
     .                                 -mur*vu*DDSIN(Phi01))
      MNEU2(9,1)=dsqrt(g1/2d0)*(M1r*vd*DDSIN(PhiM1)
     .                                 -mur*vu*DDSIN(Phi01))
      MNEU2(4,6)=dsqrt(g1/2d0)*(-M1r*vd*DDSIN(PhiM1)
     .                                 +mur*vu*DDSIN(Phi01))
      MNEU2(6,4)=dsqrt(g1/2d0)*(-M1r*vd*DDSIN(PhiM1)
     .                                 +mur*vu*DDSIN(Phi01))

      MNEU2(5,1)=0d0
      MNEU2(1,5)=0d0
      MNEU2(6,10)=0d0
      MNEU2(10,6)=0d0
      MNEU2(10,1)=0d0
      MNEU2(1,10)=0d0
      MNEU2(10,1)=0d0
      MNEU2(1,10)=0d0
      MNEU2(6,5)=0d0
      MNEU2(5,6)=0d0

      MNEU2(2,2)=M2r**2+g2/2d0*(vu**2+vd**2)
      MNEU2(7,7)=M2r**2+g2/2d0*(vu**2+vd**2)
      MNEU2(2,7)=0d0
      MNEU2(7,2)=0d0

      MNEU2(2,3)=-dsqrt(g2/2d0)*(M2r*vu*DDCOS(PhiM2)
     .                                 +mur*vd*DDCOS(Phi01))
      MNEU2(3,2)=-dsqrt(g2/2d0)*(M2r*vu*DDCOS(PhiM2)
     .                                 +mur*vd*DDCOS(Phi01))
      MNEU2(7,8)=-dsqrt(g2/2d0)*(M2r*vu*DDCOS(PhiM2)
     .                                 +mur*vd*DDCOS(Phi01))
      MNEU2(8,7)=-dsqrt(g2/2d0)*(M2r*vu*DDCOS(PhiM2)
     .                                 +mur*vd*DDCOS(Phi01))
      MNEU2(2,8)=dsqrt(g2/2d0)*(M2r*vu*DDSIN(PhiM2)
     .                                 -mur*vd*DDSIN(Phi01))
      MNEU2(8,2)=dsqrt(g2/2d0)*(M2r*vu*DDSIN(PhiM2)
     .                                 -mur*vd*DDSIN(Phi01))
      MNEU2(3,7)=dsqrt(g2/2d0)*(-M2r*vu*DDSIN(PhiM2)
     .                                 +mur*vd*DDSIN(Phi01))
      MNEU2(7,3)=dsqrt(g2/2d0)*(-M2r*vu*DDSIN(PhiM2)
     .                                 +mur*vd*DDSIN(Phi01))

      MNEU2(2,4)=dsqrt(g2/2d0)*(M2r*vd*DDCOS(PhiM2)
     .                                 +mur*vu*DDCOS(Phi01))
      MNEU2(4,2)=dsqrt(g2/2d0)*(M2r*vd*DDCOS(PhiM2)
     .                                 +mur*vu*DDCOS(Phi01))
      MNEU2(7,9)=dsqrt(g2/2d0)*(M2r*vd*DDCOS(PhiM2)
     .                                 +mur*vu*DDCOS(Phi01))
      MNEU2(9,7)=dsqrt(g2/2d0)*(M2r*vd*DDCOS(PhiM2)
     .                                 +mur*vu*DDCOS(Phi01))
      MNEU2(2,9)=dsqrt(g2/2d0)*(-M2r*vd*DDSIN(PhiM2)
     .                                 +mur*vu*DDSIN(Phi01))
      MNEU2(9,2)=dsqrt(g2/2d0)*(-M2r*vd*DDSIN(PhiM2)
     .                                 +mur*vu*DDSIN(Phi01))
      MNEU2(7,4)=dsqrt(g2/2d0)*(M2r*vd*DDSIN(PhiM2)
     .                                 -mur*vu*DDSIN(Phi01))
      MNEU2(4,7)=dsqrt(g2/2d0)*(M2r*vd*DDSIN(PhiM2)
     .                                 -mur*vu*DDSIN(Phi01))

      MNEU2(2,5)=0d0
      MNEU2(5,2)=0d0
      MNEU2(7,10)=0d0
      MNEU2(10,7)=0d0
      MNEU2(2,10)=0d0
      MNEU2(10,2)=0d0
      MNEU2(5,7)=0d0
      MNEU2(7,5)=0d0

      MNEU2(3,3)=mur**2+l**2*vd**2+(g1+g2)/2d0*vu**2
      MNEU2(8,8)=mur**2+l**2*vd**2+(g1+g2)/2d0*vu**2
      MNEU2(3,8)=0d0
      MNEU2(8,3)=0d0

      MNEU2(3,4)=(l**2-(g1+g2)/2d0)*vu*vd
      MNEU2(4,3)=(l**2-(g1+g2)/2d0)*vu*vd
      MNEU2(8,9)=(l**2-(g1+g2)/2d0)*vu*vd
      MNEU2(9,8)=(l**2-(g1+g2)/2d0)*vu*vd
      MNEU2(3,9)=0d0
      MNEU2(9,3)=0d0
      MNEU2(4,8)=0d0
      MNEU2(8,4)=0d0

      MNEU2(3,5)=l*mur*vu-l*vd*(mupsi*DDCOS(Phi01-phip)
     .                              +ks2si*DDCOS(Phi01-Phi02))
      MNEU2(5,3)=l*mur*vu-l*vd*(mupsi*DDCOS(Phi01-phip)
     .                              +ks2si*DDCOS(Phi01-Phi02))
      MNEU2(8,10)=l*mur*vu-l*vd*(mupsi*DDCOS(Phi01-phip)
     .                              +ks2si*DDCOS(Phi01-Phi02))
      MNEU2(10,8)=l*mur*vu-l*vd*(mupsi*DDCOS(Phi01-phip)
     .                              +ks2si*DDCOS(Phi01-Phi02))
      MNEU2(3,10)=l*vd*(mupsi*DDSIN(Phi01-phip)
     .                           +ks2si*DDSIN(Phi01-Phi02))
      MNEU2(10,3)=l*vd*(mupsi*DDSIN(Phi01-phip)
     .                           +ks2si*DDSIN(Phi01-Phi02))
      MNEU2(5,8)=-l*vd*(mupsi*DDSIN(Phi01-phip)
     .                           +ks2si*DDSIN(Phi01-Phi02))
      MNEU2(8,5)=-l*vd*(mupsi*DDSIN(Phi01-phip)
     .                           +ks2si*DDSIN(Phi01-Phi02))

      MNEU2(4,4)=mur**2+l**2*vu**2+(g1+g2)/2d0*vd**2
      MNEU2(9,9)=mur**2+l**2*vu**2+(g1+g2)/2d0*vd**2
      MNEU2(4,9)=0d0
      MNEU2(9,4)=0d0

      MNEU2(4,5)=l*mur*vd-l*vu*(mupsi*DDCOS(Phi01-phip)
     .                              +ks2si*DDCOS(Phi01-Phi02))
      MNEU2(5,4)=l*mur*vd-l*vu*(mupsi*DDCOS(Phi01-phip)
     .                              +ks2si*DDCOS(Phi01-Phi02))
      MNEU2(9,10)=l*mur*vd-l*vu*(mupsi*DDCOS(Phi01-phip)
     .                              +ks2si*DDCOS(Phi01-Phi02))
      MNEU2(10,9)=l*mur*vd-l*vu*(mupsi*DDCOS(Phi01-phip)
     .                              +ks2si*DDCOS(Phi01-Phi02))
      MNEU2(4,10)=l*vu*(mupsi*DDSIN(Phi01-phip)
     .                           +ks2si*DDSIN(Phi01-Phi02))
      MNEU2(10,4)=l*vu*(mupsi*DDSIN(Phi01-phip)
     .                           +ks2si*DDSIN(Phi01-Phi02))
      MNEU2(5,9)=-l*vu*(mupsi*DDSIN(Phi01-phip)
     .                           +ks2si*DDSIN(Phi01-Phi02))
      MNEU2(9,5)=-l*vu*(mupsi*DDSIN(Phi01-phip)
     .                           +ks2si*DDSIN(Phi01-Phi02))

      MNEU2(5,5)=msi**2+l**2*(vu**2+vd**2)
      MNEU2(10,10)=msi**2+l**2*(vu**2+vd**2)
      MNEU2(5,10)=0d0
      MNEU2(10,5)=0d0

c       2) Diagonalization - Eigenvalues
      CALL DIAGN(10,MNEU2,VALP,VECP,1d-10)
      CALL SORTNA(10,VALP,VECP)
      DO I=1,5
      IF(VALP(2*I-1).le.0d0)THEN
c       print*,'neupro'
       IFAIL=7
       GOTO 617
      ENDIF
       MNEU(I)=VALP(2*I-1)

c      aux1=VECP(1,2*I-1)**2+VECP(2,2*I-1)**2+VECP(3,2*I-1)**2
c     c +VECP(4,2*I-1)**2+VECP(5,2*I-1)**2
c      aux2=VECP(6,2*I-1)**2+VECP(7,2*I-1)**2+VECP(8,2*I-1)**2
c     c +VECP(9,2*I-1)**2+VECP(10,2*I-1)**2
c      IF(aux1.ge.aux2)THEN
       DO J=1,5
        NEU(I,J,1)=VECP(J,2*I-1)
        NEU(I,J,2)=VECP(J+5,2*I-1)
       ENDDO
c      ELSE
c       DO J=1,5
c        NEU(I,J,1)=VECP(J,2*I)
c        NEU(I,J,2)=-VECP(J+5,2*I)
c       ENDDO
c      ENDIF
      ENDDO

c       3) Phase shift
c      Pseudomasses
        DO I=1,5
      MNEUP(I,1)=0d0
      MNEUP(I,2)=0d0
      MNEUP(I,1)=MNEUP(I,1)
     . +M1r*DDCOS(PhiM1)*(NEU(I,1,1)*NEU(I,1,1)-NEU(I,1,2)*NEU(I,1,2))
     . +M1r*DDSIN(PhiM1)*(NEU(I,1,2)*NEU(I,1,1)+NEU(I,1,1)*NEU(I,1,2))
     . +M2r*DDCOS(PhiM2)*(NEU(I,2,1)*NEU(I,2,1)-NEU(I,2,2)*NEU(I,2,2))
     . +M2r*DDSIN(PhiM2)*(NEU(I,2,2)*NEU(I,2,1)+NEU(I,2,1)*NEU(I,2,2))
     . -2d0*mur*DDCOS(Phi01)
     .                 *(NEU(I,3,1)*NEU(I,4,1)-NEU(I,3,2)*NEU(I,4,2))
     . -2d0*mur*DDSIN(Phi01)
     .                 *(NEU(I,3,2)*NEU(I,4,1)+NEU(I,3,1)*NEU(I,4,2))
     . +2d0*dsqrt(g1/2d0)*vu
     .                 *(NEU(I,3,1)*NEU(I,1,1)-NEU(I,3,2)*NEU(I,1,2))
     . -2d0*dsqrt(g1/2d0)*vd
     .                 *(NEU(I,4,1)*NEU(I,1,1)-NEU(I,4,2)*NEU(I,1,2))
     . -2d0*dsqrt(g2/2d0)*vu
     .                 *(NEU(I,3,1)*NEU(I,2,1)-NEU(I,3,2)*NEU(I,2,2))
     . +2d0*dsqrt(g2/2d0)*vd
     .                 *(NEU(I,4,1)*NEU(I,2,1)-NEU(I,4,2)*NEU(I,2,2))
     . +(mupsi*DDCOS(phip)+ks2si*DDCOS(Phi02))
     .                 *(NEU(I,5,1)*NEU(I,5,1)-NEU(I,5,2)*NEU(I,5,2))
     . +(mupsi*DDSIN(phip)+ks2si*DDSIN(Phi02))
     .                 *(NEU(I,5,1)*NEU(I,5,2)+NEU(I,5,2)*NEU(I,5,1))
     . -2d0*l*vd*DDCOS(Phi01)
     .                 *(NEU(I,3,1)*NEU(I,5,1)-NEU(I,3,2)*NEU(I,5,2))
     . -2d0*l*vd*DDSIN(Phi01)
     .                 *(NEU(I,3,1)*NEU(I,5,2)+NEU(I,3,2)*NEU(I,5,1))
     . -2d0*l*vu*DDCOS(Phi01)
     .                 *(NEU(I,4,1)*NEU(I,5,1)-NEU(I,4,2)*NEU(I,5,2))
     . -2d0*l*vu*DDSIN(Phi01)
     .                 *(NEU(I,4,1)*NEU(I,5,2)+NEU(I,4,2)*NEU(I,5,1))
      MNEUP(I,2)=MNEUP(I,2)
     . -M1r*DDCOS(PhiM1)*(NEU(I,1,2)*NEU(I,1,1)+NEU(I,1,1)*NEU(I,1,2))
     . +M1r*DDSIN(PhiM1)*(NEU(I,1,1)*NEU(I,1,1)-NEU(I,1,2)*NEU(I,1,2))
     . -M2r*DDCOS(PhiM2)*(NEU(I,2,2)*NEU(I,2,1)+NEU(I,2,1)*NEU(I,2,2))
     . +M2r*DDSIN(PhiM2)*(NEU(I,2,1)*NEU(I,2,1)-NEU(I,2,2)*NEU(I,2,2))
     . +2d0*mur*DDCOS(Phi01)
     .                 *(NEU(I,3,2)*NEU(I,4,1)+NEU(I,3,1)*NEU(I,4,2))
     . -2d0*mur*DDSIN(Phi01)
     .                 *(NEU(I,3,1)*NEU(I,4,1)-NEU(I,3,2)*NEU(I,4,2))
     . -2d0*dsqrt(g1/2d0)*vu
     .                 *(NEU(I,3,2)*NEU(I,1,1)+NEU(I,3,1)*NEU(I,1,2))
     . +2d0*dsqrt(g1/2d0)*vd
     .                 *(NEU(I,4,2)*NEU(I,1,1)+NEU(I,4,1)*NEU(I,1,2))
     . +2d0*dsqrt(g2/2d0)*vu
     .                 *(NEU(I,3,2)*NEU(I,2,1)+NEU(I,3,1)*NEU(I,2,2))
     . -2d0*dsqrt(g2/2d0)*vd
     .                 *(NEU(I,4,1)*NEU(I,2,2)+NEU(I,4,2)*NEU(I,2,1))
     . -(mupsi*DDCOS(phip)+ks2si*DDCOS(Phi02))
     .                 *(NEU(I,5,1)*NEU(I,5,2)+NEU(I,5,2)*NEU(I,5,1))
     . +(mupsi*DDSIN(phip)+ks2si*DDSIN(Phi02))
     .                 *(NEU(I,5,1)*NEU(I,5,1)-NEU(I,5,2)*NEU(I,5,2))
     . +2d0*l*vd*DDCOS(Phi01)
     .                 *(NEU(I,3,1)*NEU(I,5,2)+NEU(I,3,2)*NEU(I,5,1))
     . -2d0*l*vd*DDSIN(Phi01)
     .                 *(NEU(I,3,1)*NEU(I,5,1)-NEU(I,3,2)*NEU(I,5,2))
     . +2d0*l*vu*DDCOS(Phi01)
     .                 *(NEU(I,4,1)*NEU(I,5,2)+NEU(I,4,2)*NEU(I,5,1))
     . -2d0*l*vu*DDSIN(Phi01)
     .                 *(NEU(I,4,1)*NEU(I,5,1)-NEU(I,4,2)*NEU(I,5,2))
        ENDDO

      DO I=1,5
      If(dabs(MNEUP(I,1)).ge.dabs(MNEUP(I,2))*1.d-10)THEN
       aux=datan(MNEUP(I,2)/MNEUP(I,1))
       IF(MNEUP(I,1).lt.0d0)THEN
        IF(MNEUP(I,2).ge.0d0)aux=aux+Pi
        IF(MNEUP(I,2).lt.0d0)aux=aux-Pi
       ENDIF
      ELSE
      if(MNEUP(I,1).ge.0d0)THEN
       aux=Pi/2d0
      else
       aux=-Pi/2d0
      endif
      ENDIF
c      aux=-aux

      DO J=1,5
       aux1=NEU(I,J,1)
       aux2=NEU(I,J,2)
       NEU(I,J,1)=aux1*DDCOS(aux/2d0)-aux2*DDSIN(aux/2d0)
       NEU(I,J,2)=aux1*DDSIN(aux/2d0)+aux2*DDCOS(aux/2d0)
c       NEUCOMP(I,J)=DCMPLX(NEU(I,J,1),NEU(I,J,2))
      ENDDO
      ENDDO

        DO I=1,5
      MNEUP(I,1)=0d0
      MNEUP(I,2)=0d0
      MNEUP(I,1)=MNEUP(I,1)
     . +M1r*DDCOS(PhiM1)*(NEU(I,1,1)*NEU(I,1,1)-NEU(I,1,2)*NEU(I,1,2))
     . +M1r*DDSIN(PhiM1)*(NEU(I,1,2)*NEU(I,1,1)+NEU(I,1,1)*NEU(I,1,2))
     . +M2r*DDCOS(PhiM2)*(NEU(I,2,1)*NEU(I,2,1)-NEU(I,2,2)*NEU(I,2,2))
     . +M2r*DDSIN(PhiM2)*(NEU(I,2,2)*NEU(I,2,1)+NEU(I,2,1)*NEU(I,2,2))
     . -2d0*mur*DDCOS(Phi01)
     .                 *(NEU(I,3,1)*NEU(I,4,1)-NEU(I,3,2)*NEU(I,4,2))
     . -2d0*mur*DDSIN(Phi01)
     .                 *(NEU(I,3,2)*NEU(I,4,1)+NEU(I,3,1)*NEU(I,4,2))
     . +2d0*dsqrt(g1/2d0)*vu
     .                 *(NEU(I,3,1)*NEU(I,1,1)-NEU(I,3,2)*NEU(I,1,2))
     . -2d0*dsqrt(g1/2d0)*vd
     .                 *(NEU(I,4,1)*NEU(I,1,1)-NEU(I,4,2)*NEU(I,1,2))
     . -2d0*dsqrt(g2/2d0)*vu
     .                 *(NEU(I,3,1)*NEU(I,2,1)-NEU(I,3,2)*NEU(I,2,2))
     . +2d0*dsqrt(g2/2d0)*vd
     .                 *(NEU(I,4,1)*NEU(I,2,1)-NEU(I,4,2)*NEU(I,2,2))
     . +(mupsi*DDCOS(phip)+ks2si*DDCOS(Phi02))
     .                 *(NEU(I,5,1)*NEU(I,5,1)-NEU(I,5,2)*NEU(I,5,2))
     . +(mupsi*DDSIN(phip)+ks2si*DDSIN(Phi02))
     .                 *(NEU(I,5,1)*NEU(I,5,2)+NEU(I,5,2)*NEU(I,5,1))
     . -2d0*l*vd*DDCOS(Phi01)
     .                 *(NEU(I,3,1)*NEU(I,5,1)-NEU(I,3,2)*NEU(I,5,2))
     . -2d0*l*vd*DDSIN(Phi01)
     .                 *(NEU(I,3,1)*NEU(I,5,2)+NEU(I,3,2)*NEU(I,5,1))
     . -2d0*l*vu*DDCOS(Phi01)
     .                 *(NEU(I,4,1)*NEU(I,5,1)-NEU(I,4,2)*NEU(I,5,2))
     . -2d0*l*vu*DDSIN(Phi01)
     .                 *(NEU(I,4,1)*NEU(I,5,2)+NEU(I,4,2)*NEU(I,5,1))
      MNEUP(I,2)=MNEUP(I,2)
     . -M1r*DDCOS(PhiM1)*(NEU(I,1,2)*NEU(I,1,1)+NEU(I,1,1)*NEU(I,1,2))
     . +M1r*DDSIN(PhiM1)*(NEU(I,1,1)*NEU(I,1,1)-NEU(I,1,2)*NEU(I,1,2))
     . -M2r*DDCOS(PhiM2)*(NEU(I,2,2)*NEU(I,2,1)+NEU(I,2,1)*NEU(I,2,2))
     . +M2r*DDSIN(PhiM2)*(NEU(I,2,1)*NEU(I,2,1)-NEU(I,2,2)*NEU(I,2,2))
     . +2d0*mur*DDCOS(Phi01)
     .                 *(NEU(I,3,2)*NEU(I,4,1)+NEU(I,3,1)*NEU(I,4,2))
     . -2d0*mur*DDSIN(Phi01)
     .                 *(NEU(I,3,1)*NEU(I,4,1)-NEU(I,3,2)*NEU(I,4,2))
     . -2d0*dsqrt(g1/2d0)*vu
     .                 *(NEU(I,3,2)*NEU(I,1,1)+NEU(I,3,1)*NEU(I,1,2))
     . +2d0*dsqrt(g1/2d0)*vd
     .                 *(NEU(I,4,2)*NEU(I,1,1)+NEU(I,4,1)*NEU(I,1,2))
     . +2d0*dsqrt(g2/2d0)*vu
     .                 *(NEU(I,3,2)*NEU(I,2,1)+NEU(I,3,1)*NEU(I,2,2))
     . -2d0*dsqrt(g2/2d0)*vd
     .                 *(NEU(I,4,1)*NEU(I,2,2)+NEU(I,4,2)*NEU(I,2,1))
     . -(mupsi*DDCOS(phip)+ks2si*DDCOS(Phi02))
     .                 *(NEU(I,5,1)*NEU(I,5,2)+NEU(I,5,2)*NEU(I,5,1))
     . +(mupsi*DDSIN(phip)+ks2si*DDSIN(Phi02))
     .                 *(NEU(I,5,1)*NEU(I,5,1)-NEU(I,5,2)*NEU(I,5,2))
     . +2d0*l*vd*DDCOS(Phi01)
     .                 *(NEU(I,3,1)*NEU(I,5,2)+NEU(I,3,2)*NEU(I,5,1))
     . -2d0*l*vd*DDSIN(Phi01)
     .                 *(NEU(I,3,1)*NEU(I,5,1)-NEU(I,3,2)*NEU(I,5,2))
     . +2d0*l*vu*DDCOS(Phi01)
     .                 *(NEU(I,4,1)*NEU(I,5,2)+NEU(I,4,2)*NEU(I,5,1))
     . -2d0*l*vu*DDSIN(Phi01)
     .                 *(NEU(I,4,1)*NEU(I,5,1)-NEU(I,4,2)*NEU(I,5,2))
        ENDDO

c      print*,'M11',M1r*DDCOS(PhiM1),M1r*DDSIN(PhiM1)
c      print*,'M12',0d0
c      print*,'M13',dsqrt(g1/2d0)*vu
c      print*,'M14',-dsqrt(g1/2d0)*vd
c      print*,'M15',0d0
c      print*,'M21',0d0
c      print*,'M22',M2r*DDCOS(PhiM2),M2r*DDSIN(PhiM2)
c      print*,'M23',-dsqrt(g2/2d0)*vu
c      print*,'M24',dsqrt(g2/2d0)*vd
c      print*,'M25',0d0
c      print*,'M31',dsqrt(g1/2d0)*vu
c      print*,'M32',-dsqrt(g2/2d0)*vu
c      print*,'M33',0d0
c      print*,'M34',-mur*DDCOS(Phi01),-mur*DDSIN(Phi01)
c      print*,'M35',-l*vd*DDCOS(Phi01),-l*vd*DDSIN(Phi01)
c      print*,'M41',-dsqrt(g1/2d0)*vd
c      print*,'M42',dsqrt(g2/2d0)*vd
c      print*,'M43',-mur*DDCOS(Phi01),-mur*DDSIN(Phi01)
c      print*,'M44',0d0
c      print*,'M45',-l*vu*DDCOS(Phi01),-l*vu*DDSIN(Phi01)
c      print*,'M51',0d0
c      print*,'M52',0d0
c      print*,'M53',-l*vd*DDCOS(Phi01),-l*vd*DDSIN(Phi01)
c      print*,'M54',-l*vu*DDCOS(Phi01),-l*vu*DDSIN(Phi01)
c      print*,'M55',msi*DDCOS(Phi02),msi*DDSIN(Phi02)

c      print*,'mneu',dsqrt(MNEU(1)),dsqrt(MNEU(2)),dsqrt(MNEU(3)),
c     . dsqrt(MNEU(4)),dsqrt(MNEU(5))
c      print*,'N11',NEU(1,1,1),NEU(1,1,2),MNEU2(1,1),MNEU2(1,6)
c      print*,'N12',NEU(1,2,1),NEU(1,2,2),MNEU2(1,2),MNEU2(1,7)
c      print*,'N13',NEU(1,3,1),NEU(1,3,2),MNEU2(1,3),MNEU2(1,8)
c      print*,'N14',NEU(1,4,1),NEU(1,4,2),MNEU2(1,4),MNEU2(1,9)
c      print*,'N15',NEU(1,5,1),NEU(1,5,2),MNEU2(1,5),MNEU2(1,10)
c      print*,'N21',NEU(2,1,1),NEU(2,1,2),MNEU2(2,1),MNEU2(2,6)
c      print*,'N22',NEU(2,2,1),NEU(2,2,2),MNEU2(2,2),MNEU2(2,7)
c      print*,'N23',NEU(2,3,1),NEU(2,3,2),MNEU2(2,3),MNEU2(2,8)
c      print*,'N24',NEU(2,4,1),NEU(2,4,2),MNEU2(2,4),MNEU2(2,9)
c      print*,'N25',NEU(2,5,1),NEU(2,5,2),MNEU2(2,5),MNEU2(2,10)
c      print*,'N31',NEU(3,1,1),NEU(3,1,2),MNEU2(3,1),MNEU2(3,6)
c      print*,'N32',NEU(3,2,1),NEU(3,2,2),MNEU2(3,2),MNEU2(3,7)
c      print*,'N33',NEU(3,3,1),NEU(3,3,2),MNEU2(3,3),MNEU2(3,8)
c      print*,'N34',NEU(3,4,1),NEU(3,4,2),MNEU2(3,4),MNEU2(3,9)
c      print*,'N35',NEU(3,5,1),NEU(3,5,2),MNEU2(3,5),MNEU2(3,10)
c      print*,'N41',NEU(4,1,1),NEU(4,1,2),MNEU2(4,1),MNEU2(4,6)
c      print*,'N42',NEU(4,2,1),NEU(4,2,2),MNEU2(4,2),MNEU2(4,7)
c      print*,'N43',NEU(4,3,1),NEU(4,3,2),MNEU2(4,3),MNEU2(4,8)
c      print*,'N44',NEU(4,4,1),NEU(4,4,2),MNEU2(4,4),MNEU2(4,9)
c      print*,'N45',NEU(4,5,1),NEU(4,5,2),MNEU2(4,5),MNEU2(4,10)
c      print*,'N51',NEU(5,1,1),NEU(5,1,2),MNEU2(5,1),MNEU2(5,6)
c      print*,'N52',NEU(5,2,1),NEU(5,2,2),MNEU2(5,2),MNEU2(5,7)
c      print*,'N53',NEU(5,3,1),NEU(5,3,2),MNEU2(5,3),MNEU2(5,8)
c      print*,'N54',NEU(5,4,1),NEU(5,4,2),MNEU2(5,4),MNEU2(5,9)
c      print*,'N55',NEU(5,5,1),NEU(5,5,2),MNEU2(5,5),MNEU2(5,10)

 617  RETURN
      END
