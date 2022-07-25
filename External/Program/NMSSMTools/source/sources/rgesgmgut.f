      SUBROUTINE RGESGMGUT(PROB,IFAIL,GUTEST)

***********************************************************************
*   Subroutine to integrate the RG equations for the gauge and Yukawa
*   couplings from MMESS up to the GUT scale (which is determined here)
*   through a CALL of the subroutine ODEINTGMGUT that is part of
*   the file integgmgut.f
*
*   It checks whether there is a Landau Pole below M_GUT
*   for the couplings lambda, kappa, HT and HB
*   If yes: PROB(27) =/= 0
*
***********************************************************************

      IMPLICIT NONE

      INTEGER I,IFAIL,JM,JL,GMFLAG,NN,CFLAG(5)
      PARAMETER (NN=15)

      DOUBLE PRECISION PROB(*),EPS,X1,X2,Y(NN)
      DOUBLE PRECISION PI,COEF,YMAX,GUTEST
      DOUBLE PRECISION G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTGUT,HBGUT,HLGUT
      DOUBLE PRECISION MSUSYEFF,MMESS,N5,MGUT,MSREF,D,DMIN
      DOUBLE PRECISION LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU
      DOUBLE PRECISION LPPGUT,LTTGUT,LUGUT,LDGUT,LTGUT,LBGUT,LLGUT
      DOUBLE PRECISION MSM,MST,LM,LT,KM,KT,HTM,HTT,LPPM,LPPT,LTTM,LTTT

      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG
      COMMON/MESCOUP/G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTGUT,HBGUT,HLGUT
      COMMON/MESGUT/LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU
      COMMON/GUTMES/LPPGUT,LTTGUT,LUGUT,LDGUT,LTGUT,LBGUT,LLGUT
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/MGUT/MGUT
      COMMON/GMSAVE/MSM,MST,LM,LT,KM,KT,HTM,HTT,
     . LPPM,LPPT,LTTM,LTTT,JM,JL
      COMMON/CFLAG/CFLAG

      EXTERNAL DERIVSGMGUT,RKQSGMGUT

      EPS=1d-8
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)

* Definition of the couplings squared Y(I) at MMESS

      Y(1)=G1MES
      Y(2)=G2MES
      Y(3)=G3MES
      Y(4)=DSQRT(LM)
      Y(5)=KM
      Y(6)=DSQRT(HTM)
      Y(7)=DSQRT(HBMES)
      Y(8)=DSQRT(HLMES)
      Y(9)=LPPM
      Y(10)=LTTM
      Y(11)=LUMES
      Y(12)=LDMES
      Y(13)=LTMES
      Y(14)=LBMES
      Y(15)=LLMES

      X1=0d0
      X2=(3d0/G1MES-5d0/G2MES)/28d0

!      WRITE(0,*)"CALL RGESGMGUT"
!      WRITE(0,*)""
!      WRITE(0,*)"MMESS =",MMESS
!      WRITE(0,*)"G1 =",5d0/3d0*Y(1)
!      WRITE(0,*)"G2 =",Y(2)
!      WRITE(0,*)"G3 =",Y(3)
!      WRITE(0,*)"L =",Y(4)
!      WRITE(0,*)"K =",Y(5)
!      WRITE(0,*)"HT =",Y(6)
!      WRITE(0,*)"HB =",Y(7)
!      WRITE(0,*)"HL =",Y(8)
!      WRITE(0,*)"LPP =",Y(9)
!      WRITE(0,*)"LTT =",Y(10)
!      WRITE(0,*)"LU =",Y(11)
!      WRITE(0,*)"LD =",Y(12)
!      WRITE(0,*)"LT =",Y(13)
!      WRITE(0,*)"LB =",Y(14)
!      WRITE(0,*)"LL =",Y(15)
!      WRITE(0,*)""

      CALL ODEINTGMGUT(Y,NN,X1,X2,EPS,DERIVSGMGUT,RKQSGMGUT,IFAIL)

      IF(IFAIL.GT.0)THEN
!       WRITE(0,*)"IFAIL =",IFAIL
!       WRITE(0,*)""
!       WRITE(0,*)""
       IFAIL=14
      ENDIF
!      WRITE(0,*)""

      YMAX=0d0
      DO I=1,3
       YMAX=MAX(YMAX,Y(I))
      ENDDO
      DO I=4,NN
       YMAX=MAX(YMAX,Y(I)**2)
      ENDDO

      IF(CFLAG(1).EQ.1)PROB(27)=DDIM(YMAX/(4d0*PI),1d0)
      
* The GUT scale in GeV:

      MGUT=MMESS*DEXP(8d0*PI**2*X2)

* Couplings at the GUT scale

      G1GUT=Y(1)
      G2GUT=Y(2)
      G3GUT=Y(3)
      LGUT=Y(4)
      KGUT=Y(5)
      HTGUT=Y(6)
      HBGUT=Y(7)
      HLGUT=Y(8)
      LPPGUT=Y(9)
      LTTGUT=Y(10)
      LUGUT=Y(11)
      LDGUT=Y(12)
      LTGUT=Y(13)
      LBGUT=Y(14)
      LLGUT=Y(15)

      GUTEST=(LPPGUT-XIU)**2+(LTTGUT-XIU)**2

!      WRITE(0,*)"MGUT =",MGUT
!      WRITE(0,*)"G1 =",5d0/3d0*Y(1)
!      WRITE(0,*)"G2 =",Y(2)
!      WRITE(0,*)"G3 =",Y(3)
!      WRITE(0,*)"L =",Y(4)
!      WRITE(0,*)"K =",Y(5)
!      WRITE(0,*)"HT =",Y(6)
!      WRITE(0,*)"HB =",Y(7)
!      WRITE(0,*)"HL =",Y(8)
!      WRITE(0,*)"LPP =",Y(9)
!      WRITE(0,*)"LTT =",Y(10)
!      WRITE(0,*)"LU =",Y(11)
!      WRITE(0,*)"LD =",Y(12)
!      WRITE(0,*)"LT =",Y(13)
!      WRITE(0,*)"LB =",Y(14)
!      WRITE(0,*)"LL =",Y(15)
!      WRITE(0,*)""
!      WRITE(0,*)""

      END
