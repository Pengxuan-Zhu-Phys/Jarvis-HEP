      SUBROUTINE RGESGM(PAR,IFAIL)

***********************************************************************
*   Subroutine to integrate the RG equations for the gauge and Yukawa
*   couplings up to the messenger scale MMESS, through a CALL of the
*   subroutine ODEINTGM that is part of the file integgm.f
*
*   It checks whether there is a Landau Pole below MMES
*   for the couplings lambda, kappa, htop, hbot and htau
*   If yes: IFAIL = 11
*
*   Below Q2 all sparticle/heavy Higgs thresholds are taken into
*   account in the naive step function approximation.
*   Above the Susy scale Q2 the two loop beta functions are used.
*   (Note: The sparticle thresholds are consistent even if a
*   sparticle mass is above Q2: then the threshold effect between
*   MT and Q2 "anticipates" the threshold effect between Q2 and MMESS)
*
***********************************************************************

      IMPLICIT NONE

      INTEGER I,IFAIL,GMFLAG,JM,JL,NN
      PARAMETER (NN=8)

      DOUBLE PRECISION PAR(*),EPS,X1,X2,Y(NN)
      DOUBLE PRECISION PI,Q2,QSTSB,YMAX
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      DOUBLE PRECISION SIGNKAPPA,MSUSYEFF,MMESS,N5
      DOUBLE PRECISION MSM,MST,LM,LT,KM,KT,HTM,HTT,LPPM,LPPT,LTTM,LTTT
      DOUBLE PRECISION LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES
      DOUBLE PRECISION LLMES,XIU,MSREF,D,DMIN

      COMMON/RENSCALE/Q2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/MESCOUP/G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES      
      COMMON/STSBSCALE/QSTSB      
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/GMSAVE/MSM,MST,LM,LT,KM,KT,HTM,HTT,
     . LPPM,LPPT,LTTM,LTTT,JM,JL
      COMMON/MESGUT/LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU
      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG

      EXTERNAL DERIVS,RKQSGM

      EPS=1d-8
      PI=4d0*DATAN(1d0)

* CALL GETSUSYCOUP(PAR) in order to determine the gauge and
* Yukawa couplings in "COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS" at Q

      CALL GETSUSYCOUP(PAR)

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
      X2=DLOG(MMESS**2/Q2)/(4d0*PI)**2
      
!      WRITE(0,*)"CALL RGESGM"
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

      CALL ODEINTGM(Y,NN,X1,X2,EPS,DERIVS,RKQSGM,IFAIL)
      
!      WRITE(0,*)"MMESS =",MMESS
!      WRITE(0,*)"G1 =",5d0/3d0*Y(1)
!      WRITE(0,*)"G2 =",Y(2)
!      WRITE(0,*)"G3 =",Y(3)
!      WRITE(0,*)"L2 =",Y(4)
!      WRITE(0,*)"K2 =",Y(5)
!      WRITE(0,*)"HT2 =",Y(6)
!      WRITE(0,*)"HB2 =",Y(7)
!      WRITE(0,*)"HL2 =",Y(8)
!      WRITE(0,*)""

      YMAX=0d0
      DO I=1,NN
       YMAX=MAX(YMAX,Y(I))
      ENDDO

      IF(YMAX.GT.4d0*PI)THEN
       IFAIL=5
      ENDIF

      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"IFAIL =",IFAIL
!       WRITE(0,*)""
!       WRITE(0,*)""
       IFAIL=11
       RETURN
      ENDIF
!      WRITE(0,*)""

* Couplings at the messenger scale

      G1MES=Y(1)
      G2MES=Y(2)
      G3MES=Y(3)
      LM=Y(4)
      LMES=LM
      KM=DSQRT(Y(5))*SIGNKAPPA
      KMES=KM
      HTM=Y(6)
      HTMES=HTM
      HBMES=Y(7)
      HLMES=Y(8)
      IF(GMFLAG.EQ.3 .OR. GMFLAG.EQ.4)THEN
       LDMES=LTMES*DSQRT(LM/HTM)
       LUMES=LBMES*DSQRT(LM/HBMES)
       LLMES=LBMES*DSQRT(HLMES/HBMES)
      ENDIF

      END
