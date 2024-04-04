      SUBROUTINE RGESINVGM(PAR,IFAIL)

***********************************************************************
*   Subroutine to integrate the 2-loop RGEs for all soft terms
*   from MMESS down to Q2 = M_SUSY, through a CALL of the
*   subroutine ODEINTS that is part of the file integs.f
*
*   At MMESS:
*   MS^2 = MSMES (MAFLAG=-1,-3), MSINP (MAFLAG=-2,-4)
*   XIS = XISMES (MAFLAG=-2,-4), XISINP (MAFLAG=-1,-3)
*   XIF = XIFMES (MAFLAG=-3,-4), XIFINP (MAFLAG=-1,-2)
*
*   It uses COMMON/INPPAR and /MESEXT for the soft terms,
*   COMMON/MESCOUP for the gauge/Yukawa couplings at MMES.
*
*   The output (soft terms at Q2) is written into PAR(*) and
*   COMMON/SUSYMH and /SUSYEXT
*
***********************************************************************

      IMPLICIT NONE

      INTEGER IFAIL,OMGFLAG,MAFLAG,MOFLAG,GMFLAG,JM,JL,NN
      PARAMETER (NN=35)

      DOUBLE PRECISION PAR(*),EPS,X1,X2,Y(NN),PI,ALP1,ALP2,ALP3
      DOUBLE PRECISION Q2,COEF
      DOUBLE PRECISION G1,G2,G3,LMES,KMES,HTMES,HB,HL
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MSM,MST,LM,LT,KM,KT,HTM,HTT,LPPM,LPPT,LTTM,LTTT
      DOUBLE PRECISION MH1S,MH2S,MSS,MSUSY,MMESS,N5
      DOUBLE PRECISION ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      DOUBLE PRECISION XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      DOUBLE PRECISION LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION M1INP,M2INP,M3INP,AKINP,ATINP,ABINP,
     .      ATAUINP,AMUINP,MH1INP,MH2INP,MQ3INP,MU3INP,MD3INP,
     .      MQINP,MUINP,MDINP,ML3INP,ME3INP,MLINP,MEINP
      DOUBLE PRECISION X,F1,F2,SP2,MSREF,D,DMIN

      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/RENSCALE/Q2
      COMMON/INPPAR/ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      COMMON/SOFTINP/M1INP,M2INP,M3INP,AKINP,ATINP,ABINP,
     .      ATAUINP,AMUINP,MH1INP,MH2INP,MQ3INP,MU3INP,MD3INP,
     .      MQINP,MUINP,MDINP,ML3INP,ME3INP,MLINP,MEINP
      COMMON/MESCAL/MSUSY,MMESS,N5
      COMMON/MESCOUP/G1,G2,G3,LMES,KMES,HTMES,HB,HL      
      COMMON/MESEXT/XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      COMMON/MESGUT/LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/SUSYMH/MH1S,MH2S,MSS
      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG
      COMMON/GMSAVE/MSM,MST,LM,LT,KM,KT,HTM,HTT,
     . LPPM,LPPT,LTTM,LTTT,JM,JL

      EXTERNAL DERIVSS,RKQSS

      EPS=1d-8
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)

* Definition of the couplings squared Y(I) at MMESS

      Y(1)=G1
      Y(2)=G2
      Y(3)=G3
      Y(4)=LM
* NOTE: KM = KAPPA, NOT KAPPA**2
      Y(5)=KM
      Y(6)=HTM
      Y(7)=HB
      Y(8)=HL

* Corrections higher order in X=MSUSY/MMESS:

      X=MSUSY/MMESS
      IF(X.GE.1d-7) THEN
        IF(X.LT.1d0) THEN
          F1=((1d0+X)*DLOG(1d0+X)+(1d0-X)*DLOG(1d0-X))/X**2
          F2=(1d0+X)/X**2*(DLOG(1d0+X)-2d0*SP2(X/(1d0+X))
     .       +SP2(2d0*X/(1d0+X))/2d0)
     .       +(1d0-X)/X**2*(DLOG(1d0-X)-2d0*SP2(-X/(1d0-X))
     .       +SP2(-2d0*X/(1d0-X))/2d0)
        ELSE
          F1=0d0
          F2=0d0
        ENDIF
      ELSE
        F1=1d0
        F2=1d0
      ENDIF

* (SP2 = Li_2  is defined in the subroutine bsg.f)
      
* Input values for the soft terms at MMESS:

      ALP1=G1/(4d0*PI)
      ALP2=G2/(4d0*PI)
      ALP3=G3/(4d0*PI)
      
      Y(9)=5d0/3d0*ALP1*N5*MSUSY*F1/(4d0*PI)
      Y(10)=ALP2*N5*MSUSY*F1/(4d0*PI)
      Y(11)=ALP3*N5*MSUSY*F1/(4d0*PI)
      IF(GMFLAG.NE.0)THEN
       ALINP=-COEF*MSUSY*
     .       (2d0*LPPM**2+3d0*LTTM**2+3d0*LUMES**2+3d0*LDMES**2)
      ENDIF
      Y(12)=ALINP
      IF(PAR(2).NE.0d0)THEN
       IF(GMFLAG.NE.0)THEN
        Y(13)=-3d0*COEF*MSUSY*
     .        (2d0*LPPM**2+3d0*LTTM**2+2d0*LUMES**2+2d0*LDMES**2)
       ELSE
        Y(13)=3d0*ALINP
       ENDIF
      ELSE
       Y(13)=0d0
      ENDIF
      Y(14)=-COEF*MSUSY*(LUMES**2+3d0*LTMES**2+LBMES**2)
      Y(15)=-COEF*MSUSY*(LDMES**2+LTMES**2+3d0*LBMES**2)
      Y(16)=-COEF*MSUSY*(LDMES**2+3d0*LLMES**2)
      Y(17)=COEF*(5d0/6d0*ALP1**2+3d0/2d0*ALP2**2)*N5*F2*MSUSY**2
     .     +(COEF*MSUSY)**2*(
     .     - 3d0*HTM*LBMES**2 - 3d0*LDMES**2*LM - 2d0*LM*LPPM**2
     .     - 9d0*HTM*LTMES**2 - 3d0*LM*LTTM**2
     .     + 3d0*LPPM**2*LUMES**2 + 3d0*LTTM**2*LUMES**2
     .     + 6d0*DSQRT(HB*LM)*LBMES*LUMES + 2d0*DSQRT(HL*LM)*LLMES*LUMES
     .     + 2d0*LDMES*DSQRT(LM)*LPPM*LUMES - G1*LUMES**2
     .     - 3d0*G2*LUMES**2 + 2d0*KM**2*LUMES**2
     .     + 3d0*LBMES**2*LUMES**2 + 2d0*LDMES**2*LUMES**2
     .     + LLMES**2*LUMES**2 + 2d0*LM*LUMES**2 + 4d0*LUMES**4)
     .     - MSUSY**4/MMESS**2*COEF/6d0*LUMES**2
      Y(18)=COEF*(5d0/6d0*ALP1**2+3d0/2d0*ALP2**2)*N5*F2*MSUSY**2
     .     +(COEF*MSUSY)**2*(
     .     - 9d0*HB*LBMES**2 - G1*LDMES**2 - 3d0*G2*LDMES**2
     .     + 2d0*KM**2*LDMES**2 + 4d0*LDMES**4 - 3d0*HL*LLMES**2
     .     + 2d0*LDMES**2*LM + 3d0*LDMES**2*LPPM**2
     .     - 2d0*LM*LPPM**2 + 6d0*DSQRT(HTM*LM)*LDMES*LTMES
     .     - 3d0*HB*LTMES**2 + 3d0*LDMES**2*LTMES**2
     .     + 3d0*LDMES**2*LTTM**2 - 3d0*LM*LTTM**2
     .     + 2d0*LDMES*DSQRT(LM)*LPPM*LUMES + 2d0*LDMES**2*LUMES**2
     .     - 3d0*LM*LUMES**2) - MSUSY**4/MMESS**2*COEF/6d0*LDMES**2
      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
       Y(19)=MSM
      ELSE
       Y(19)=MSINP
      ENDIF
      Y(20)=COEF*(5d0/54d0*ALP1**2+3d0/2d0*ALP2**2
     .     +8d0/3d0*ALP3**2)*F2*N5*MSUSY**2
     .     +(COEF*MSUSY)**2*(
     .     - 7d0/9d0*G1*LBMES**2 - 3d0*G2*LBMES**2
     .     - 16d0/3d0*G3*LBMES**2 + 6d0*HB*LBMES**2 + 6d0*LBMES**4
     .     - HB*LDMES**2 + LBMES**2*LLMES**2
     .     + 2d0*DSQRT(HB*HL)*LBMES*LLMES
     .     + 2d0*DSQRT(HB)*LBMES*LDMES*LPPM + LBMES**2*LPPM**2
     .     + 2d0*DSQRT(HTM*LM)*LDMES*LTMES - 13d0/9d0*G1*LTMES**2
     .     - 3d0*G2*LTMES**2 - 16d0/3d0*G3*LTMES**2 + 6d0*HTM*LTMES**2
     .     + 2d0*LBMES**2*LTMES**2 + LDMES**2*LTMES**2
     .     + LPPM**2*LTMES**2 + 6d0*LTMES**4
     .     + 2d0*DSQRT(HB*LM)*LBMES*LUMES
     .     + 2d0*DSQRT(HTM)*LPPM*LTMES*LUMES
     .     - HTM*LUMES**2 + LBMES**2*LUMES**2)
     .     - MSUSY**4/MMESS**2*COEF/6d0*(LTMES**2+LBMES**2)
      Y(21)=COEF*(40d0/27d0*ALP1**2
     .     +8d0/3d0*ALP3**2)*F2*N5*MSUSY**2
     .     +(COEF*MSUSY)**2*(
     .     - 2d0*HTM*LBMES**2 + 4d0*DSQRT(HTM*LM)*LDMES*LTMES
     .     - 26d0/9d0*G1*LTMES**2 - 6d0*G2*LTMES**2
     .     - 32d0/3d0*G3*LTMES**2 + 2d0*HB*LTMES**2 - 2d0*HTM*LUMES**2
     .     + 12d0*HTM*LTMES**2 + 2d0*LBMES**2*LTMES**2
     .     + 2d0*LDMES**2*LTMES**2 + 2d0*LPPM**2*LTMES**2
     .     + 12d0*LTMES**4 + 4d0*DSQRT(HTM)*LPPM*LTMES*LUMES)
     .     - MSUSY**4/MMESS**2*COEF/3d0*LTMES**2
      Y(22)=COEF*(10d0/27d0*ALP1**2
     .     +8d0/3d0*ALP3**2)*F2*N5*MSUSY**2
     .     +(COEF*MSUSY)**2*(
     .     - 14d0/9d0*G1*LBMES**2 - 6d0*G2*LBMES**2
     .     - 32d0/3d0*G3*LBMES**2 + 12d0*HB*LBMES**2 + 2d0*HTM*LBMES**2
     .     + 12d0*LBMES**4 - 2d0*HB*LDMES**2
     .     + 4d0*DSQRT(HB*HL)*LBMES*LLMES
     .     + 4d0*DSQRT(HB)*LBMES*LDMES*LPPM
     .     + 2d0*LBMES**2*LLMES**2 + 2d0*LBMES**2*LPPM**2
     .     - 2d0*HB*LTMES**2 + 2d0*LBMES**2*LTMES**2
     .     + 4d0*DSQRT(HB*LM)*LBMES*LUMES + 2d0*LBMES**2*LUMES**2)
     .     - MSUSY**4/MMESS**2*COEF/3d0*LBMES**2
      Y(23)=COEF*(5d0/54d0*ALP1**2+3d0/2d0*ALP2**2
     .     +8d0/3d0*ALP3**2)*F2*N5*MSUSY**2
      Y(24)=COEF*(40d0/27d0*ALP1**2
     .     +8d0/3d0*ALP3**2)*F2*N5*MSUSY**2
      Y(25)=COEF*(10d0/27d0*ALP1**2
     .     +8d0/3d0*ALP3**2)*F2*N5*MSUSY**2
      Y(26)=COEF*(5d0/6d0*ALP1**2
     .     +3d0/2d0*ALP2**2)*F2*N5*MSUSY**2
     .     +(COEF*MSUSY)**2*(
     .     - HL*LDMES**2 + 6d0*DSQRT(HB*HL)*LBMES*LLMES
     .     - 3d0*G1*LLMES**2 - 3d0*G2*LLMES**2 + 2d0*HL*LLMES**2
     .     + 3d0*LBMES**2*LLMES**2 + 4d0*LLMES**4
     .     + 2d0*DSQRT(HL)*LDMES*LLMES*LPPM + LLMES**2*LPPM**2
     .     + 2d0*DSQRT(HL*LM)*LLMES*LUMES + LLMES**2*LUMES**2)
     .     - MSUSY**4/MMESS**2*COEF/6d0*LLMES**2
      Y(27)=COEF*(10d0/3d0*ALP1**2)*F2*N5*MSUSY**2
     .     +(COEF*MSUSY)**2*(
     .     - 2*HL*LDMES**2 + 12*DSQRT(HB*HL)*LBMES*LLMES
     .     - 6*G1*LLMES**2 - 6*G2*LLMES**2 + 4*HL*LLMES**2
     .     + 6*LBMES**2*LLMES**2 + 8*LLMES**4
     .     + 4*DSQRT(HL)*LDMES*LLMES*LPPM + 2*LLMES**2*LPPM**2
     .     + 4*DSQRT(HL*LM)*LLMES*LUMES + 2*LLMES**2*LUMES**2)
     .     - MSUSY**4/MMESS**2*COEF/3d0*LLMES**2
      Y(28)=COEF*(5d0/6d0*ALP1**2
     .     +3d0/2d0*ALP2**2)*F2*N5*MSUSY**2
      Y(29)=COEF*(10d0/3d0*ALP1**2)*F2*N5*MSUSY**2
      IF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4)THEN
       Y(30)=XIFMES
      ELSE
       Y(30)=XIFINP
      ENDIF
      IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)THEN
       Y(31)=XISMES
      ELSE
       Y(31)=XISINP
      ENDIF
      Y(32)=MUPINP
      Y(33)=MSPINP
      Y(34)=0d0
      Y(35)=0d0

* Save the required input values for soft masses and couplings in
* COMMON/SOFTINP

      M1INP=Y(9)
      M2INP=Y(10)
      M3INP=Y(11)
      IF(GMFLAG.NE.0)THEN
       ALINP=Y(12)
      ENDIF
      AKINP=Y(13)
      ATINP=Y(14)
      ABINP=Y(15)
      ATAUINP=Y(16)
      MH1INP=Y(17)
      MH2INP=Y(18)
      MQ3INP=Y(20)
      MU3INP=Y(21)
      MD3INP=Y(22)
      MQINP=Y(23)
      MUINP=Y(24)
      MDINP=Y(25)
      ML3INP=Y(26)
      ME3INP=Y(27)
      MLINP=Y(28)
      MEINP=Y(29)
      AMUINP=Y(35)

      X1=COEF*DLOG(MMESS**2/Q2)
      X2=0d0

!      WRITE(0,*)"CALL RGESINVGM"
!      WRITE(0,*)""
!      WRITE(0,*)"MMESS =",MMESS
!      WRITE(0,*)"G1MES =",DSQRT(Y(1))
!      WRITE(0,*)"G2MES =",DSQRT(Y(2))
!      WRITE(0,*)"G3MES =",DSQRT(Y(3))
!      WRITE(0,*)"LMES =",DSQRT(Y(4))
!      WRITE(0,*)"KMES =",Y(5)
!      WRITE(0,*)"HTMES =",DSQRT(Y(6))
!      WRITE(0,*)"HBMES =",DSQRT(Y(7))
!      WRITE(0,*)"HLMES =",DSQRT(Y(8))
!      WRITE(0,*)"M1MES =",Y(9)
!      WRITE(0,*)"M2MES =",Y(10)
!      WRITE(0,*)"M3MES =",Y(11)
!      WRITE(0,*)"ALMES =",Y(12)
!      WRITE(0,*)"AKMES =",Y(13)
!      WRITE(0,*)"ATOPMES =",Y(14)
!      WRITE(0,*)"ABOTMES =",Y(15)
!      WRITE(0,*)"ATAUMES =",Y(16)
!      WRITE(0,*)"AMUMES =",Y(35)
!      WRITE(0,*)"MH1MES =",Y(17)
!      WRITE(0,*)"MH2MES =",Y(18)
!      WRITE(0,*)"MSMES =",Y(19)
!      WRITE(0,*)"MQ3MES =",Y(20)
!      WRITE(0,*)"MU3MES =",Y(21)
!      WRITE(0,*)"MD3MES =",Y(22)
!      WRITE(0,*)"MQMES =",Y(23)
!      WRITE(0,*)"MUMES =",Y(24)
!      WRITE(0,*)"MDMES =",Y(25)
!      WRITE(0,*)"ML3MES =",Y(26)
!      WRITE(0,*)"ME3MES =",Y(27)
!      WRITE(0,*)"MLMES =",Y(28)
!      WRITE(0,*)"MEMES =",Y(29)
!      WRITE(0,*)"XIFMES =",Y(30)
!      WRITE(0,*)"XISMES =",Y(31)
!      WRITE(0,*)"MUPMES =",Y(32)
!      WRITE(0,*)"MSPMES =",Y(33)
!      WRITE(0,*)"M3HMES =",Y(34)
!      WRITE(0,*)""

      CALL ODEINTS(Y,NN,X1,X2,EPS,DERIVSS,RKQSS,IFAIL)      

!      WRITE(0,*)"MSUSY =",DSQRT(Q2)
!      WRITE(0,*)"G1 =",DSQRT(Y(1))
!      WRITE(0,*)"G2 =",DSQRT(Y(2))
!      WRITE(0,*)"G3 =",DSQRT(Y(3))
!      WRITE(0,*)"L =",DSQRT(Y(4))
!      WRITE(0,*)"K =",Y(5)
!      WRITE(0,*)"HT =",DSQRT(Y(6))
!      WRITE(0,*)"HB =",DSQRT(Y(7))
!      WRITE(0,*)"HL =",DSQRT(Y(8))
!      WRITE(0,*)"M1 =",Y(9)
!      WRITE(0,*)"M2 =",Y(10)
!      WRITE(0,*)"M3 =",Y(11)
!      WRITE(0,*)"AL =",Y(12)
!      WRITE(0,*)"AK =",Y(13)
!      WRITE(0,*)"ATOP =",Y(14)
!      WRITE(0,*)"ABOT =",Y(15)
!      WRITE(0,*)"ATAU =",Y(16)
!      WRITE(0,*)"AMUON =",Y(35)
!      WRITE(0,*)"MH1 =",Y(17)
!      WRITE(0,*)"MH2 =",Y(18)
!      WRITE(0,*)"MS =",Y(19)
!      WRITE(0,*)"MQ3 =",Y(20)
!      WRITE(0,*)"MU3 =",Y(21)
!      WRITE(0,*)"MD3 =",Y(22)
!      WRITE(0,*)"MQ =",Y(23)
!      WRITE(0,*)"MU =",Y(24)
!      WRITE(0,*)"MD =",Y(25)
!      WRITE(0,*)"ML3 =",Y(26)
!      WRITE(0,*)"ME3 =",Y(27)
!      WRITE(0,*)"ML =",Y(28)
!      WRITE(0,*)"ME =",Y(29)
!      WRITE(0,*)"XIF =",Y(30)
!      WRITE(0,*)"XIS =",Y(31)
!      WRITE(0,*)"MUP =",Y(32)
!      WRITE(0,*)"MSP =",Y(33)
!      WRITE(0,*)"M3H =",Y(34)
!      WRITE(0,*)""

      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"IFAIL =",IFAIL
!       WRITE(0,*)""
!       WRITE(0,*)""
       IFAIL=13
       RETURN
      ENDIF
!      WRITE(0,*)""

* Y(5) = KAPPA, NOT KAPPA**2

      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)PAR(2)=Y(5)

* SOFT TERMS AT THE SUSY SCALE

      PAR(5)=Y(12)
      PAR(6)=Y(13)
      PAR(7)=Y(20)
      PAR(8)=Y(21)
      PAR(9)=Y(22)
      PAR(10)=Y(26)
      PAR(11)=Y(27)
      PAR(12)=Y(14)
      PAR(13)=Y(15)
      PAR(14)=Y(16)
      PAR(15)=Y(23)
      PAR(16)=Y(24)
      PAR(17)=Y(25)
      PAR(18)=Y(28)
      PAR(19)=Y(29)
      PAR(20)=Y(9)
      PAR(21)=Y(10)
      PAR(22)=Y(11)
      PAR(25)=Y(35)

* MH1S, MH2S AND MSS at Q2 are stored in COMMON/SUSYMH:
      
      MH1S=Y(17)
      MH2S=Y(18)
      MSS=Y(19)
      
* EXT parameters at Q2, stored in COMMON/SUSYEXT:

      XIFSUSY=Y(30)
      XISSUSY=Y(31)
      MUPSUSY=Y(32)
      MSPSUSY=Y(33)
      M3HSUSY=Y(34)
      
      END
