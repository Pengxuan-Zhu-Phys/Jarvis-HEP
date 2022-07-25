      SUBROUTINE RGES(PAR,PROB,IFAIL,FLAG)

***********************************************************************
*   Subroutine to integrate the RG equations for the gauge and Yukawa
*   couplings from Q2 up to the GUT scale (which is determined here)
*   through a CALL of the subroutine ODEINT that is part of
*   the file integ.f
*
*   It checks whether there is a Landau Pole below M_GUT
*   for the couplings lambda, kappa, htop, hbot and htau
*   If yes: PROB(27) =/= 0
*
*   The gauge/Yukawa couplings at Q2 are taken from COMMON/SUSYCOUP/
*   defined in the subroutine GETSUSYCOUP(PAR)
*   where now DELMB is included for HBOTS
*   Above the Susy scale Q2 the two loop beta functions are used.
*   (Note: The sparticle thresholds are consistent even if a
*   sparticle mass is above Q2: then the threshold effect between
*   MT and Q2 "anticipates" the threshold effect between Q2 and MGUT)
*
***********************************************************************

      IMPLICIT NONE

      INTEGER I,IFAIL,NN,FLAG,CFLAG(5)
      PARAMETER (NN=8)

      DOUBLE PRECISION PAR(*),PROB(*),EPS,X1,X2,Y(NN)
      DOUBLE PRECISION PI,QSTSB,Q2,YMAX
      DOUBLE PRECISION MGUT,g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     .      HBOTGUT,HTAUGUT
      DOUBLE PRECISION SIGNKAPPA

      COMMON/RENSCALE/Q2
      COMMON/MGUT/MGUT
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     .      HBOTGUT,HTAUGUT      
      COMMON/STSBSCALE/QSTSB      
      COMMON/CFLAG/CFLAG

      EXTERNAL DERIVS,RKQS

      EPS=1d-8
      PI=4d0*DATAN(1d0)

* CALL GETSUSYCOUP(PAR) in order to determine the gauge and
* Yukawa couplings in "COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS" at Q

      IF(FLAG.NE.1)THEN
       CALL GETSUSYCOUP(PAR)
      ENDIF

      IF(PAR(2).GE.0d0)THEN
       SIGNKAPPA=1d0
      ELSE
       SIGNKAPPA=-1d0
      ENDIF


* Definition of the couplings squared Y(I) at M_SUSY

      Y(1)=g1s
      Y(2)=g2s
      Y(3)=g3s
      Y(4)=PAR(1)**2
      Y(5)=PAR(2)**2
      Y(6)=HTOPS**2
      Y(7)=HBOTS**2
      Y(8)=HTAUS**2

      X1=0d0
      X2=(3d0/g1s-5d0/g2s)/28d0

!      WRITE(0,*)"CALL RGES"
!      WRITE(0,*)""
!      WRITE(0,*)"MSUSY =",DSQRT(Q2)
!      WRITE(0,*)"G1 =",5d0/3d0*Y(1)
!      WRITE(0,*)"G2 =",Y(2)
!      WRITE(0,*)"G3 =",Y(3)
!      WRITE(0,*)"L2 =",Y(4)
!      WRITE(0,*)"K2 =",Y(5)
!      WRITE(0,*)"HT2 =",Y(6)
!      WRITE(0,*)"HB2 =",Y(7)
!      WRITE(0,*)"HL2 =",Y(8)
!      WRITE(0,*)""

      CALL ODEINT(Y,NN,X1,X2,EPS,DERIVS,RKQS,IFAIL)

* The GUT scale in GeV:

      MGUT=DSQRT(Q2)*DEXP(8d0*PI**2*X2)

!      WRITE(0,*)"MGUT =",MGUT
!      WRITE(0,*)"G1 =",5d0/3d0*Y(1)
!      WRITE(0,*)"G2 =",Y(2)
!      WRITE(0,*)"G3 =",Y(3)
!      WRITE(0,*)"L2 =",Y(4)
!      WRITE(0,*)"K2 =",Y(5)
!      WRITE(0,*)"HT2 =",Y(6)
!      WRITE(0,*)"HB2 =",Y(7)
!      WRITE(0,*)"HL2 =",Y(8)
!      WRITE(0,*)""

      IF(IFAIL.GT.0)THEN
!       WRITE(0,*)"IFAIL =",IFAIL
!       WRITE(0,*)""
!       WRITE(0,*)""
       IFAIL=11
       RETURN
      ELSE
!       WRITE(0,*)""
       IFAIL=0
      ENDIF

      YMAX=0d0
      DO I=1,NN
       YMAX=MAX(YMAX,Y(I))
      ENDDO

      IF(CFLAG(1).EQ.1)PROB(27)=DDIM(YMAX/(4d0*PI),1d0)
      
* Couplings at the GUT scale

      G1GUT=Y(1)
      G2GUT=Y(2)
      G3GUT=Y(3)
      LGUT=Y(4)
      KGUT=DSQRT(Y(5))*SIGNKAPPA
      HTOPGUT=Y(6)
      HBOTGUT=Y(7)
      HTAUGUT=Y(8)

      END
