      SUBROUTINE MSFERM(PAR,IFAIL,FLAG)
      
*******************************************************************
* Subroutine to compute the squark pole masses/mixing angles and
* slepton masses/mixing angles (with thanks to S. Kraml).
*
* The (large tan(beta)) SUSY correction DELMB is computed here.
*
* Negative squark masses squared give IFAIL = 8.
*
* The running masses of the 3rd generation of squarks are first
* evaluated at a scale QSTSB taking all possible thresholds
* (i.e. possibly large logs) into account. Pole masses are computed
* including 1 loop (S)QCD rad. corrs.
* No electroweak rad. corrs. are included at present.
*
* Concerning the 1st and 2nd generation squarks, it is assumed that
* the input parameters PAR(15)=MQ, PAR(16)=MU, PAR(17)=MD are the
* running squark masses at the SUSY scale Q2 ~ MQ^2 ~ MU^2 ~ MD^2.
*
* The conventions for the cosines/sines CSF/SSF (F=T, B, L=Stau) of the
* 3rd generation sfermion mixing angle are
* F1 = CSF*FL + SSF*FR, F2 = CSF*FR - SSF*FL
* where FR, FL are the left/right handed weak eigenstates, and
* F1, F2 are the physical states ordered in mass (M_F1 < M_F2)
*******************************************************************

      IMPLICIT NONE

      INTEGER IFAIL,OMGFLAG,MAFLAG,MOFLAG,FLAG,PFLAG

      DOUBLE PRECISION PAR(*)
      DOUBLE PRECISION cor,Q2,QSTSB,pi
      DOUBLE PRECISION Wt,At,mstL,mstR
      DOUBLE PRECISION Wb,Ab,msbL,msbR
      DOUBLE PRECISION Wl,Atau,AMUON,Xl,mslL,mslR
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW,X
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      DOUBLE PRECISION MGL,nen,fac,DMSQUARK
      DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      DOUBLE PRECISION M1,M2,M3,MQ3,MU3,MD3,ML3,ME3,MQ,MU,MD,ML,ME
      DOUBLE PRECISION MQ3P,MU3P,MD3P,MH1S,MH2S,MSS,B,GM
      DOUBLE PRECISION LQSTSB,SB1,SB2,ST1,ST2
      DOUBLE PRECISION COEF,ATP,ABP
      DOUBLE PRECISION INTEG,DELMB,DELML,DEL1
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION HTQ,HBQ,MTQ,MBQ
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      DOUBLE PRECISION XI,HT2,HB2,MH1,MH2,L2,K2,HTAU2,AL,AK
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS,SP,SIG1,SIG2,SIG3
      DOUBLE PRECISION F20,F21,F22,F14,F15
      DOUBLE PRECISION DPI11ST,DPI22ST,DPI11SB,DPI22SB

      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      COMMON/RADCOR2/MQ3P,MU3P,MD3P,ATP,ABP
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QQUARK/HTQ,HBQ,MTQ,MBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/QEXT/XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      COMMON/SUSYMH/MH1S,MH2S,MSS
      COMMON/PFLAG/PFLAG
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS

      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)

!      WRITE(0,*)"CALL MSFERM"
!      WRITE(0,*)""

      AT=PAR(12)
      AB=PAR(13)
      ATAU=PAR(14)
      AMUON=PAR(25)

*   In order to run the sfermion masses from Q2 to QSTSB
*   tree level expressions for the soft
*   Higgs masses and various logarithms are needed:

      IF(MAFLAG.GE.0)THEN
* Otherwise: taken from COMMON/SUSYMH
       B=ALQ+NUQ
       GM= LQ*XIFQ+MUPQ*MUQ+M3HQ
       MH1S= -LQ**2*H2Q**2 - MUQ**2 + (MUQ*B+GM)/TANBQ
     .      + GQ/2d0*(H2Q**2-H1Q**2)
       MH2S= -LQ**2*H1Q**2 - MUQ**2 + (MUQ*B+GM)*TANBQ
     .      + GQ/2d0*(H1Q**2-H2Q**2)
       MSS=-LQ**2*(H1Q**2+H2Q**2) - 2d0*NUQ**2
     .     +LQ**2*H1Q*H2Q/MUQ*(ALQ+2d0*NUQ+MUPQ) - NUQ*AKQ
     .     -XIFQ*(2d0*KQ+LQ*MUPQ/MUQ)-MUPQ**2-3d0*MUPQ*NUQ
     .     -MSPQ-LQ*XISQ/MUQ
      ENDIF

      M1=PAR(20)
      M2=PAR(21)
      M3=PAR(22)
      MGL=PAR(22)
      
      MQ3=PAR(7)
      MU3=PAR(8)
      MD3=PAR(9)
      ML3=PAR(10)
      ME3=PAR(11)
      MQ=PAR(15)
      MU=PAR(16)
      MD=PAR(17)
      ML=PAR(18)
      ME=PAR(19)

      LQSTSB=DLOG(Q2/QSTSB)

*   Running stop/sbottom masses squared and mixings:
*   On input, they are defined at the scale Q2.
*   First, they have to be evaluated at the scale QSTSB
*   (as used in the calculation of the pole masses)
*   Use 2-loop RGEs from Q2 to QSTSB in order to comply
*   with SLHA (input parameters in DR_bar without thresholds)

* CALL GETSUSYCOUP(PAR) in order to determine the gauge and
* Yukawa couplings in "COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS" at Q

      IF(FLAG.NE.1)THEN
       CALL GETSUSYCOUP(PAR)
      ENDIF

      XI= g1S*(MH1S-MH2S+MQ3+2d0*MQ-2d0*(MU3+2d0*MU)+MD3+2d0*MD
     .  +ME3+2d0*ME-ML3+2d0*ML)
      HT2=HTOPS**2
      HB2=HBOTS**2
      HTAU2=HTAUS**2
      MH1=MH1S
      MH2=MH2S
      L2=PAR(1)**2
      K2=PAR(2)**2
      AL=PAR(5)
      AK=PAR(6)

      SP= HT2*(-3d0*MH1 - MQ3 + 4d0*MU3)
     .  + HB2*(3d0*MH2 - MQ3 - 2d0*MD3)
     .  + HTAU2*(MH2 + ML3 - 2d0*ME3) + L2*(MH2 - MH1)
     .  + (G1S/18d0 + 3d0/2d0*G2S + 8d0/3d0*G3S)*(MQ3+2d0*MQ)
     .  - (16d0/9d0*G1S + 16d0/3d0*G3S)*(MU3+2d0*MU)
     .  + (2d0/9d0*G1S + 8d0/3d0*G3S)*(MD3+2d0*MD)
     .  + (G1S/2d0 + 3d0/2d0*G2S)*(MH1-MH2-(ML3+2d0*ML)) 
     .  + 2d0*G1S*(ME3+2d0*ME)
      SIG1= G1S*(MH1 + MH2 + (MQ3+2d0*MQ)/3d0 + 8d0/3d0*(MU3+2d0*MU) 
     .    + 2d0/3d0*(MD3+2d0*MD)+ ML3+2d0*ML + 2d0*(ME3+2d0*ME))
      SIG2=G2S*(MH1 + MH2 + 3d0*(MQ3+2d0*MQ) + ML3+2d0*ML)
      SIG3=G3S*(2d0*(MQ3+2d0*MQ) + MU3+2d0*MU + MD3+2d0*MD)

* Copy F20, F21, F22 from integs.f, and S -> XI, 
* G1,G2,G3 -> G1S,G2S,G3S, MS -> MSS

      F20= HT2*(MH1+MQ3+MU3+AT**2) + HB2*(MH2+MQ3+MD3+AB**2)
     .       - g1S*M1**2/9d0 - 3d0*g2S*M2**2
     .       - 16d0/3d0*g3S*M3**2 + XI/6d0
     .       + COEF*(-10d0*HT2**2*(MH1+MQ3+MU3+2d0*AT**2)
     .       - 10d0*HB2**2*(MH2+MQ3+MD3+2d0*AB**2)
     .       - 10d0*HB2**4*(MH2+MQ3+MD3+2d0*AB**2)
     .       - HB2*HTAU2*(2d0*MH2+MQ3+MD3+ML3+ME3+(AB+ATAU)**2)
     .       - L2*HT2*(2d0*MH1+MH2+MSS+MQ3+MU3+(AL+AT)**2)
     .       - L2*HB2*(MH1+2d0*MH2+MSS+MQ3+MD3+(AL+AB)**2)
     .       + 4d0/3d0*G1S*HT2*(MH1+MQ3+MU3+AT**2+2d0*M1*(M1-AT))
     .       + 2d0/3d0*G1S*HB2*(MH2+MQ3+MD3+AB**2+2d0*M1*(M1-AB))
     .       + 199d0/54d0*G1S**2*M1**2 + 33d0/2d0*G2S**2*M2**2
     .     - 64d0/3d0*G3S**2*M3**2
     .     + 1d0/3d0*G1S*G2S*(M1**2+M2**2+M1*M2)
     .       + 16d0/27d0*G1S*G3S*(M1**2+M3**2+M1*M3)
     .     + 16d0*G2S*G3S*(M2**2+M3**2+M2*M3)
     .       + 1d0/3d0*G1S*SP + 1d0/18d0*G1S*SIG1
     .       + 3d0/2d0*G2S*SIG2 + 8d0/3d0*G3S*SIG3)

      F21= 2d0*HT2*(MH1+MQ3+MU3+AT**2)
     .       - 16d0/9d0*g1S*M1**2 - 16d0/3d0*g3S*M3**2 - 2d0*XI/3d0
     .       + COEF*(-16d0*HT2**2*(MH1+MQ3+MU3+2d0*AT**2)
     .       - 2d0*HT2*HB2*(MH1+MH2+2d0*MQ3+MU3+MD3+(AT+AB)**2)
     .       - 2d0*HT2*L2*(2d0*MH1+MH2+MSS+MQ3+MU3+(AT+AL)**2)
     .       - 2d0/3d0*G1S*HT2*(MH1+MQ3+MU3+AT**2+2d0*M1*(M1-AT))
     .       + 6d0*G2S*HT2*(MH1+MQ3+MU3+AT**2+2d0*M2*(M2-AT))
     .       + 1712d0/27d0*G1S**2*M1**2 - 64d0/3d0*G3S**2*M3**2
     .       + 256d0/27d0*G1S*G3S*(M1**2+M3**2+M1*M3)
     .     - 4d0/3d0*G1S*SP + 8d0/9d0*G1S*SIG1 + 8d0/3d0*G3S*SIG3)

      F22= 2d0*HB2*(MH2+MQ3+MD3+AB**2)
     .       - 4d0/9d0*g1S*M1**2 - 16d0/3d0*g3S*M3**2 + XI/3d0
     .       + COEF*(-16d0*HB2**2*(MH2+MQ3+MD3+2d0*AB**2)
     .       - 2d0*HB2*HT2*(MH1+MH2+2d0*MQ3+MU3+MD3+(AB+AT)**2)
     .       - 2d0*HB2*HTAU2*(2d0*MH2+MQ3+MD3+ML3+ME3+(AB+ATAU)**2)
     .       - 2d0*HB2*L2*(MH1+2d0*MH2+MSS+MQ3+MD3+(AB+AL)**2)
     .       + 2d0/3d0*G1S*HB2*(MH2+MQ3+MD3+AB**2+2d0*M1*(M1-AB))
     .       + 6d0*G2S*HB2*(MH2+MQ3+MD3+AB**2+2d0*M2*(M2-AB))
     .       + 404d0/27d0*G1S**2*M1**2 - 64d0/3d0*G3S**2*M3**2
     .       + 64d0/27d0*G1S*G3S*(M1**2+M3**2+M1*M3)
     .     + 2d0/3d0*G1S*SP + 2d0/9d0*G1S*SIG1 + 8d0/3d0*G3S*SIG3)

      MQ3P=MQ3-COEF*LQSTSB*F20
      MU3P=MU3-COEF*LQSTSB*F21
      MD3P=MD3-COEF*LQSTSB*F22

      F14= 6d0*HT2*AT + HB2*AB + L2*AL
     .     + 13d0/9d0*g1S*M1 + 3d0*g2S*M2 + 16d0/3d0*g3S*M3
     .     + COEF*(-44d0*HT2**2*AT - 5d0*HT2*HB2*(AT+AB)
     .     - 3d0*HT2*L2*(AT+AL) - 10d0*HB2**2*AB
     .     - HB2*HTAU2*(AB+ATAU) - 4d0*HB2*L2*(AB+AL)
     .     - 6d0*L2**2*AL - L2*HTAU2*(AL+ATAU)
     .     - 2d0*L2*K2*(AL+AK)
     .     + 2d0*g1S*HT2*(AT-M1) + 6d0*g2S*HT2*(AT-M2)
     .     + 16d0*g3S*HT2*(AT-M3) + 2d0/3d0*g1S*HB2*(AB-M1)
     .     - 2743d0/81d0*g1S**2*M1 - 5d0/3d0*g1S*g2S*(M1+M2)
     .     - 136d0/27d0*g1S*g3S*(M1+M3) - 15d0*g2S**2*M2
     .     - 8d0*g2S*g3S*(M2+M3) + 32d0/9d0*g3S**2*M3)

      F15= 6d0*HB2*AB + HT2*AT + HTAU2*ATAU + L2*AL
     .     + 7d0/9d0*g1S*M1 + 3d0*g2S*M2 + 16d0/3d0*g3S*M3
     .     + COEF*(-44d0*HB2**2*AB - 5d0*HT2*HB2*(AT+AB)
     .     - 3d0*HB2*HTAU2*(AB+ATAU) - 3d0*HB2*L2*(AB+AL)
     .     - 10d0*HT2**2*AT - 4d0*HT2*L2*(AT+AL)
     .     - 6d0*HTAU2**2*ATAU - 6d0*L2**2*AL - 2d0*L2*K2*(AL+AK)
     .     + 2d0/3d0*g1S*HB2*(AB-M1) + 6d0*g2S*HB2*(AB-M2)
     .     + 16d0*g3S*HB2*(AB-M3) + 4d0/3d0*g1S*HT2*(AT-M1)
     .     + 2d0*HTAU2*g1S*(ATAU-M1)
     .     - 1435d0/81d0*g1S**2*M1 - 5d0/3d0*g1S*g2S*(M1+M2)
     .     - 40d0/27d0*g1S*g3S*(M1+M3) - 15d0*g2S**2*M2
     .     - 8d0*g2S*g3S*(M2+M3) + 32d0/9d0*g3S**2*M3)


      ATP=AT-COEF*LQSTSB*F14
      ABP=AB-COEF*LQSTSB*F15

*   Running stop masses squared and mixings

      mstL= MQ3P + mtq**2 + (gQ/2d0-g1Q/3d0)*(h2q**2-h1q**2)
      mstR= MU3P + mtq**2 + g1Q/3d0*(h2q**2-h1q**2)
      Xt= ATP-muq/tanbq
      Wt= DSQRT( (mstL-mstR)**2 + 4d0*(Xt*mtq)**2)
      RMST1= 0.5d0*(mstL+mstR-Wt)
      RMST2= 0.5d0*(mstL+mstR+Wt)

      IF(RMST1.LE.1d1)THEN
       IFAIL=8
       nen=0d0
       IF(FLAG.EQ.0)RMST1=1d1
!      WRITE(0,*)"MST1^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ELSE
       nen= DSQRT((mstL-RMST1)**2 + (Xt*mtq)**2)
      ENDIF

      IF(RMST2.LE.1d1)THEN
       IFAIL=8
       IF(FLAG.EQ.0)RMST2=1d1
!      WRITE(0,*)"MST2^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(nen.EQ.0d0)THEN
       IF(mstL.LE.mstR)THEN
        CST= 1d0
       ELSE
        CST= 0d0
       ENDIF
      ELSE
       CST= -mtq*Xt/nen
      ENDIF
      S2T=2d0*CST*DSQRT(1d0-CST**2)

*   Running sbottom masses squared and mixings

      msbL= MQ3P + mbq**2 - (gQ/2d0-g1Q/6d0)*(h2q**2-h1q**2)
      msbR= MD3P + mbq**2 - g1Q/6d0*(h2q**2-h1q**2)
      Xb= ABP-muq*tanbq

      Wb= DSQRT((msbL-msbR)**2 + 4d0*(Xb*mbq)**2)
      RMSB1= 0.5d0*(msbL+msbR-Wb)
      RMSB2= 0.5d0*(msbL+msbR+Wb)

      IF(RMSB1.LE.1d1)THEN
       IFAIL=8
       nen=0d0
       IF(FLAG.EQ.0)RMSB1=1d1
!      WRITE(0,*)"MSB1^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ELSE
       nen= DSQRT((msbL-RMSB1)**2 + (Xb*mbq)**2)
      ENDIF

      IF(RMSB2.LE.1d1)THEN
       IFAIL=8
       IF(FLAG.EQ.0)RMSB2=1d1
!      WRITE(0,*)"MSB2^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(nen.EQ.0d0)THEN
       IF(msbL.LE.msbR)THEN
        CSB= 1d0
       ELSE
        CSB= 0d0
       ENDIF
      ELSE
       CSB= -mbq*Xb/nen
      ENDIF
      S2B=2*CSB*DSQRT(1d0-CSB**2)

*   Pole Squark Masses squared

      fac=ALSQ/(3d0*PI)
      MST1=RMST1
      MST2=RMST2
      MSB1=RMSB1
      MSB2=RMSB2

      CALL DMSQUARK_YUK(DPI11ST,DPI22ST,DPI11SB,DPI22SB)

      IF(MST1.GT.1d1)THEN
       MST1=RMST1-fac*DMSQUARK(1,RMST1,RMST2,S2T,mtq,MGL,QSTSB)
       MST1=MST1-DPI11ST
       MST2=RMST2-fac*DMSQUARK(2,RMST1,RMST2,S2T,mtq,MGL,QSTSB)
       MST2=MST2-DPI22ST
      ENDIF
      IF(MSB1.GT.1d1)THEN
       MSB1=RMSB1-fac*DMSQUARK(1,RMSB1,RMSB2,S2B,mbq,MGL,QSTSB)
       MSB1=MSB1-DPI11SB
       MSB2=RMSB2-fac*DMSQUARK(2,RMSB1,RMSB2,S2B,mbq,MGL,QSTSB)
       MSB2=MSB2-DPI22SB
      ENDIF

      IF(MST1.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MST1=1d2
!      WRITE(0,*)"MST1^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MST2.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MST2=1d2
!      WRITE(0,*)"MST2^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSB1.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSB1=1d2
!      WRITE(0,*)"MSB1^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSB2.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSB2=1d2
!      WRITE(0,*)"MSB2^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

*   First generation squark masses squared

      MUR=MU+g1Q/3d0*(h2q**2-h1q**2)
      MUL=MQ+(gQ/2d0-g1Q/3d0)*(h2q**2-h1q**2)
      MDR=MD-g1Q/6d0*(h2q**2-h1q**2)
      MDL=MQ+(-gQ/2d0+g1Q/6d0)*(h2q**2-h1q**2)

      IF(MUR.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MUR=1d2
!      WRITE(0,*)"MUR^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MUL.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MUL=1d2
!      WRITE(0,*)"MUL^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MDR.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MDR=1d2
!      WRITE(0,*)"MDR^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MDL.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MDL=1d2
!      WRITE(0,*)"MDL^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

*   Pole Squark Masses squared
*   On input they are defined at the scale Q2 that is assumed to not
*   too far from their actual masses
*   (Threshold effects are neglected)

      IF(MUR.GT.1d2)THEN
       X=M3**2/MUR
       IF(X.NE.1d0)THEN
        cor=1d0+3d0*X+1d0/2d0*(1d0-X)**2*DLOG((1d0-X)**2)
     .     -X**2*DLOG(X)+2d0*X*DLOG(Q2/MUR)
       ELSE
        cor=4d0+2d0*DLOG(Q2/MUR)
       ENDIF
       MUR=MUR*(1d0+2d0*fac*cor)
      ENDIF

      IF(MUL.GT.1d2)THEN
       X=M3**2/MUL
       IF(X.NE.1d0)THEN
        cor=1d0+3d0*X+1d0/2d0*(1d0-X)**2*DLOG((1d0-X)**2)
     .     -X**2*DLOG(X)+2d0*X*DLOG(Q2/MUL)
       ELSE
        cor=4d0+2d0*DLOG(Q2/MUL)
       ENDIF
       MUL=MUL*(1d0+2d0*fac*cor)
      ENDIF

      IF(MUL.GT.1d2)THEN
       X=M3**2/MDR
       IF(X.NE.1d0)THEN
        cor=1d0+3d0*X+1d0/2d0*(1d0-X)**2*DLOG((1d0-X)**2)
     .     -X**2*DLOG(X)+2d0*X*DLOG(Q2/MDR)
       ELSE
        cor=4d0+2d0*DLOG(Q2/MDR)
       ENDIF
       MDR=MDR*(1d0+2d0*fac*cor)
      ENDIF

      IF(MUL.GT.1d2)THEN
       X=M3**2/MDL
       IF(X.NE.1d0)THEN
        cor=1d0+3d0*X+1d0/2d0*(1d0-X)**2*DLOG((1d0-X)**2)
     .     -X**2*DLOG(X)+2d0*X*DLOG(Q2/MDL)
       ELSE
        cor=4d0+2d0*DLOG(Q2/MDL)
       ENDIF
       MDL=MDL*(1d0+2d0*fac*cor)
      ENDIF

      IF(MUR.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MUR=1d2
!      WRITE(0,*)"MUR^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MUL.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MUL=1d2
!      WRITE(0,*)"MUL^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MDR.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MDR=1d2
!      WRITE(0,*)"MDR^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MDL.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MDL=1d2
!      WRITE(0,*)"MDL^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MST1.GE.0d0)MST1=DSQRT(MST1)
      IF(MST2.GE.0d0)MST2=DSQRT(MST2)
      IF(MSB1.GE.0d0)MSB1=DSQRT(MSB1)
      IF(MSB2.GE.0d0)MSB2=DSQRT(MSB2)
      IF(MUL.GE.0d0)MUL=DSQRT(MUL)
      IF(MUR.GE.0d0)MUR=DSQRT(MUR)
      IF(MDL.GE.0d0)MDL=DSQRT(MDL)
      IF(MDR.GE.0d0)MDR=DSQRT(MDR)

* Calculation of the SUSY corrections to h_bot, DELMB, as in
* Carena et al., hep-ph/9912516
      
      IF(FLAG.NE.1)THEN
      SB1=DSQRT(MAX(RMSB1,1d0))
      SB2=DSQRT(MAX(RMSB2,1d0))
      ST1=DSQRT(MAX(RMST1,1d0))
      ST2=DSQRT(MAX(RMST2,1d0))

      DEL1=-2d0/(3d0*PI)*ALSQ*M3*ABP*INTEG(SB1,SB2,M3)

      DELMB=MUQ*TANBQ*(2d0/(3d0*PI)*ALSQ*M3*INTEG(SB1,SB2,M3)
     .   +COEF*HTQ**2*AT*INTEG(ST1,ST2,MUQ)
     .   -COEF*G2Q*M2*(CST**2*INTEG(ST1,M2,MUQ)
     .   +(1d0-CST**2)*INTEG(ST2,M2,MUQ)
     .   +1d0/2d0*(CSB**2*INTEG(SB1,M2,MUQ)
     .   +(1d0-CSB**2)*INTEG(SB2,M2,MUQ)))
     .   -COEF*G1Q/3d0*M1*(1d0/3d0*INTEG(SB1,SB2,M1)
     .   +1d0/2d0*(CSB**2*INTEG(SB1,M1,MUQ)
     .   +(1d0-CSB**2)*INTEG(SB2,M1,MUQ))
     .   +(1d0-CSB**2)*INTEG(SB1,M1,MUQ)
     .   +CSB**2*INTEG(SB2,M1,MUQ)))
     .   /(1d0+DEL1)

!      WRITE(0,*)"DELMB =",DELMB
!      WRITE(0,*)"DEL1 =",DEL1
!      WRITE(0,*)""
      ENDIF

*   Stau masses squared and mixings

      mslL= ML3 + MTAU**2 - (gQ/2d0-g1Q/2d0)*(h2q**2-h1q**2)
      mslR= ME3 + MTAU**2 - g1Q/2d0*(h2q**2-h1q**2)
      Xl= Atau-MUQ*tanbq
      Wl= DSQRT((mslL-mslR)**2 + 4d0*(Xl*MTAU)**2)
      msl1= 0.5d0*(mslL+mslR-Wl)
      msl2= 0.5d0*(mslL+mslR+Wl)

      nen= DSQRT((mslL-msl1)**2 + (Xl*MTAU)**2)
      IF(nen.EQ.0d0)THEN
       IF(mslL.LE.mslR)THEN
        CSL= 1d0
       ELSE
        CSL= 0d0
       ENDIF
      ELSE
       CSL= -MTAU*Xl/nen
      ENDIF

      MSNT= ML3 + gQ/2d0*(h2q**2-h1q**2)

      IF(MSL1.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSL1=1d2
!      WRITE(0,*)"MSL1^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSL2.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSL2=1d2
!      WRITE(0,*)"MSL2^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSNT.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSNT=1d2
!      WRITE(0,*)"MSNT^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSL1.GE.0d0)MSL1=DSQRT(MSL1)
      IF(MSL2.GE.0d0)MSL2=DSQRT(MSL2)
      IF(MSNT.GE.0d0)MSNT=DSQRT(MSNT)

* Calculation of the SUSY corrections to h_tau, DELML
      
      IF(FLAG.NE.1)THEN
      DELML=MUQ*TANBQ*(
     .   -COEF*G2Q*M2*(INTEG(MSNT,M2,MUQ)
     .   +1d0/2d0*(CSL**2*INTEG(MSL1,M2,MUQ)
     .   +(1d0-CSL**2)*INTEG(MSL2,M2,MUQ)))
     .   +COEF*G1Q*M1*(INTEG(MSL1,MSL2,M1)
     .   +1d0/2d0*(CSL**2*INTEG(MSL1,M1,MUQ)
     .   +(1d0-CSL**2)*INTEG(MSL2,M1,MUQ))
     .   -(1d0-CSL**2)*INTEG(MSL1,M1,MUQ)
     .   -CSL**2*INTEG(MSL2,M1,MUQ)))

!      WRITE(0,*)"DELML =",DELML
!      WRITE(0,*)""
      ENDIF

*   Smuon masses squared and mixings

      mslL= ML + MMUON**2 - (gQ/2d0-g1Q/2d0)*(h2q**2-h1q**2)
      mslR= ME + MMUON**2 - g1Q/2d0*(h2q**2-h1q**2)
      Xl= AMUON-MUQ*tanbq
      Wl= DSQRT((mslL-mslR)**2 + 4d0*(Xl*MMUON)**2)
      MSMU1= 0.5d0*(mslL+mslR-Wl)
      MSMU2= 0.5d0*(mslL+mslR+Wl)

      nen= DSQRT((mslL-MSMU1)**2 + (Xl*MMUON)**2)
      IF(nen.EQ.0d0)THEN
       IF(mslL.LE.mslR)THEN
        CSMU= 1d0
       ELSE
        CSMU= 0d0
       ENDIF
      ELSE
       CSMU= -MMUON*Xl/nen
      ENDIF

      MSNM= ML + gQ/2d0*(h2q**2-h1q**2)

      IF(MSMU1.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSMU1=1d2
!      WRITE(0,*)"MSMU1^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSMU2.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSMU2=1d2
!      WRITE(0,*)"MSMU2^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSNM.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSNM=1d2
!      WRITE(0,*)"MSNM^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSMU1.GE.0d0)MSMU1=DSQRT(MSMU1)
      IF(MSMU1.GE.0d0)MSMU2=DSQRT(MSMU2)
      IF(MSNM.GE.0d0)MSNM=DSQRT(MSNM)

*   Selectron masses squared

      MLR=ME - g1Q/2d0*(h2q**2-h1q**2)
      MLL=ML + (-gQ/2d0+g1Q/2d0)*(h2q**2-h1q**2)
      MNL=ML + gQ/2d0*(h2q**2-h1q**2)

      IF(MLR.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MLR=1d2
!      WRITE(0,*)"MLR^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MLL.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MLL=1d2
!      WRITE(0,*)"MLL^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MNL.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MNL=1d2
!      WRITE(0,*)"MNL^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MLR.GE.0d0)MLR=DSQRT(MLR)
      IF(MLL.GE.0d0)MLL=DSQRT(MLL)
      IF(MNL.GE.0d0)MNL=DSQRT(MNL)

      END


      DOUBLE PRECISION FUNCTION DMSQUARK(i,msq1,msq2,s2t,mq,mglu,q2)

*    with thanks to S. Kraml
*    msq1, msq2 are the squark masses squared

      IMPLICIT NONE
      DOUBLE PRECISION msq1,msq2,s2t,mq,mglu,msq,msqp
      DOUBLE PRECISION mq2,mglu2,sgn,dmg,dmsg,dmsq,q2
      DOUBLE PRECISION NMA0,NMB0,NMB1
      INTEGER i

      mglu=dabs(mglu)
      mq2 = mq*mq
      mglu2 = mglu*mglu
      IF (i.eq.1) THEN
        msq=msq1
        msqp=msq2
        sgn = -1d0
      ELSE
        msq=msq2
        msqp=msq1
        sgn = 1d0
      ENDIF

      dmg = -2d0*msq*( 2d0*NMB0(msq,0d0,msq,q2)-NMB1(msq,0d0,msq,q2) )
      dmsg = -4d0*( NMA0(mq2,q2) - msq*NMB1(msq,mglu2,mq2,q2) +
     .    (mglu2 + sgn*mglu*mq*s2t)*NMB0(msq,mglu2,mq2,q2) )
      dmsq = (1d0-s2t**2)*NMA0(msq,q2) + s2t**2*NMA0(msqp,q2)
      DMSQUARK = dmg + dmsg + dmsq

      END


      DOUBLE PRECISION FUNCTION INTEG(X,Y,Z)

* Function for DELMB:

      IMPLICIT NONE

      DOUBLE PRECISION X,Y,Z

       IF(DABS(X).EQ.DABS(Y) .AND. DABS(X).EQ.DABS(Z))THEN
        INTEG=.5d0/X
      ELSEIF(DABS(X).EQ.DABS(Y))THEN
        INTEG=(X**2-Z**2+Z**2*DLOG(Z**2/X**2))/(X**2-Z**2)**2
      ELSEIF(DABS(Y).EQ.DABS(Z))THEN
        INTEG=(Y**2-X**2+X**2*DLOG(X**2/Y**2))/(Y**2-X**2)**2
      ELSEIF(DABS(X).EQ.DABS(Z))THEN
        INTEG=(X**2-Y**2+Y**2*DLOG(Y**2/X**2))/(X**2-Y**2)**2
      ELSE
        INTEG=(X**2*Y**2*DLOG(X**2/Y**2)
     .   +Y**2*Z**2*DLOG(Y**2/Z**2)+Z**2*X**2*DLOG(Z**2/X**2))/
     .   ((X**2-Y**2)*(Y**2-Z**2)*(X**2-Z**2))
      ENDIF

      END


      SUBROUTINE DMSQUARK_YUK(DPI11ST,DPI22ST,DPI11SB,DPI22SB)
      
*******************************************************************
* Subroutine to compute the Yukawa contributions to the squark pole masses
* based on Pierce/Bagger et al., hep-ph/9606211 eqs. (D.47-49)
*
* NOTE: mu(Pierce/Bagger) = -MUQ(this code)!!!
*
* The conventions for the cosines/sines CSF/SSF (F=T, B) of the
* 3rd generation sfermion mixing angle are
* F1 = CSF*FL + SSF*FR, F2 = CSF*FR - SSF*FL
* where FR, FL are the left/right handed weak eigenstates, and
* F1, F2 are the physical states ordered in mass (M_F1 < M_F2)
*
* Comparing to (D.54) in hep-ph/9606211 one also finds
* F1 = CSF*FL + SSF*FR, F2 = CSF*FR - SSF*FL although F1, F2 are
* not ordered on mass, which is irrelevant
*
* Used trilin. couplings A_top, A_bottom: 
* ATP, ABP at QSTSB from COMMON/RADCOR2/
* Used running stop/sbottom masses at QSTSB:
* RMST1/2, RMSB1/2 from COMMON/RADCOR/
*
*      SCOMP(1-3,1-3): Mixing angles: IF HB(I) are the weak eigenstates,
*        HB(I) = Re(H1), Re(H2), Re(S) where H1 == H_u, H2 == H_d
*        and HM(I) are the mass eigenstates,
*        the convention is HB(I) = SUM_(J=1,3) SCOMP(J,I)*HM(J)
*        which is equivalent to HM(I) = SUM_(J=1,3) SCOMP(I,J)*HB(J)
* -> SCOMP(K,1)) is the H_u-component of the mass eigenstate K
* -> SCOMP(K,2)) is the H_d-component of the mass eigenstate K
*
*      PCOMP(1-2,1-2): Mixing angles: IF AB(I) are the weak eigenstates,
*        AB(I) = Im(H1), Im(H2), Im(S) where A1 == A_u, A2 == A_d
*        and AM(I) are the mass eigenstates,
*        the convention is
*        AM(I) = PCOMP(I,1)*(COSBETA*AB(1)+SINBETA*AB(2))
*              + PCOMP(I,2)*AB(3)
* or, using
*       P(I,1)= PCOMP(I,1)*CB
*       P(I,2)= PCOMP(I,1)*SB
*       P(I,3)= PCOMP(I,2)
*       AM(I) = P(I,1)*AB(1)+P(I,2)*AB(2)) + P(I,3)*AB(3)

*******************************************************************

      IMPLICIT NONE

      INTEGER I,J,K,N

      DOUBLE PRECISION QSTSB,pi,coef
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      DOUBLE PRECISION RMST(2),RMSB(2)
      DOUBLE PRECISION MQ3P,MU3P,MD3P,ATP,ABP
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION HTQ,HBQ,MTQ,MBQ
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION SCOMP(3,3),PCOMP(2,2),CMASS,P(2,3)
      DOUBLE PRECISION H(3,3),A(2,2),MH(3),MA(2),VEC3(3,3),VEC2(2,2),EPS
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION RT,RB,FMT,FMB,AlshIFt,B,GMCOMB,SQR2

      DOUBLE PRECISION DPI11ST,DPI22ST,DPI11SB,DPI22SB
      DOUBLE PRECISION NMA0,NMB0
      DOUBLE PRECISION A0T1,A0T2,A0B1,A0B2,A0H0(6),A0HC(2)
      DOUBLE PRECISION DNUST(6),DNDST(2),DNUSB(6),DNDSB(2)
      DOUBLE PRECISION B011HNTI(6,2),B022HNTI(6,2)
      DOUBLE PRECISION B011HCBI(2,2),B022HCBI(2,2)
      DOUBLE PRECISION B011HNBI(6,2),B022HNBI(6,2)
      DOUBLE PRECISION B011HCTI(2,2),B022HCTI(2,2)

      DOUBLE PRECISION SST,SSB,CB,SB
      DOUBLE PRECISION HRULUL,HRDLDL,HRURUR,HRDRDR,HRULUR
      DOUBLE PRECISION HRDLDR,HIULUR,HIDLDR,HGBULUR,HGBDLDR
      DOUBLE PRECISION LHNTLT(6,2),LHNTRT(6,2),LHNBLB(6,2),LHNBRB(6,2)
      DOUBLE PRECISION LHCTLBL(2),LHCTRBL(2),LHCTLBR(2),LHCTRBR(2)
      DOUBLE PRECISION LHCTLB(2,2),LHCTRB(2,2)
      DOUBLE PRECISION LHCBLTL(2),LHCBRTL(2),LHCBLTR(2),LHCBRTR(2)
      DOUBLE PRECISION LHCBLT(2,2),LHCBRT(2,2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      COMMON/RADCOR2/MQ3P,MU3P,MD3P,ATP,ABP
      COMMON/STSBSCALE/QSTSB
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QQUARK/HTQ,HBQ,MTQ,MBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/QEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

* DEFINITION OF COUPLINGS AS FUNCTIONS:

      HRULUL(K)= SQR2*(HTQ**2*H1Q*SCOMP(K,1))
      HRDLDL(K)= SQR2*(HBQ**2*H2Q*SCOMP(K,2))
      HRURUR(K)= SQR2*(HTQ**2*H1Q*SCOMP(K,1))
      HRDRDR(K)= SQR2*(HBQ**2*H2Q*SCOMP(K,2))
      HRULUR(K)= HTQ/SQR2*(ATP*SCOMP(K,1)-MUQ*SCOMP(K,2)
     . - LQ*H2Q*SCOMP(K,3))
      HRDLDR(K)= HBQ/SQR2*(-MUQ*SCOMP(K,1)+ABP*SCOMP(K,2)
     . - LQ*H1Q*SCOMP(K,3))
      HIULUR(K)= HTQ/SQR2*(ATP*P(K,1)+MUQ*P(K,2)
     . + LQ*H2Q*P(K,3))
      HIDLDR(K)= HBQ/SQR2*(MUQ*P(K,1)+ABP*P(K,2)
     . + LQ*H1Q*P(K,3))
* END: DEFINITION OF COUPLINGS AS FUNCTIONS

      EPS=1d-8
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)
      SQR2= DSQRT(2d0)

      SST= DSQRT(1d0-CST**2)
      SSB= DSQRT(1d0-CSB**2)

*   Trig. functions of beta
      cb=1d0/DSQRT(1d0+tanbq**2)
      sb=tanbq*cb

********************************************
* First: Need an estimate of Higgs masses and mixings,
* tree level + top/stop-contr.; copied from mhiggs.f:
      rt= 3d0/2d0*COEF*htq**2
      rb= 3d0/2d0*COEF*hbq**2
       fmt= (rmst2*DLOG(rmst2/QSTSB)-rmst1*DLOG(rmst1/QSTSB))/
     .     (rmst2-rmst1)-1d0

       fmb= (rmsb2*DLOG(rmsb2/QSTSB)-rmsb1*DLOG(rmsb1/QSTSB))/
     .    (rmsb2-rmsb1)-1d0
      Alshift= ALQ+2d0*rt*Atp*fmt+2d0*rb*Abp*fmb
      B= Alshift+NUQ
      GMCOMB= LQ*XIF+MUP*MUQ+M3H

*   Tree level CP-even Higgs mass matrix
      H(1,1) = GQ*h1q**2 + (MUQ*B+GMCOMB)/tanbq
      H(1,1)= H(1,1) + rt*(4d0*mtq**2*DLOG(rmst1*rmst2/mtq**4))
      H(2,2) = GQ*h2q**2 + (MUQ*B+GMCOMB)*tanbq
      H(3,3) = LQ**2*(AlshIFt+MUP)*h1q*h2q/MUQ
     .   + NUQ*(AKQ+4d0*NUQ+3d0*MUP)-LQ/MUQ*(XIS+XIF*MUP)
      H(1,2) = (2d0*LQ**2-GQ)*h1q*h2q - MUQ*B-GMCOMB
      H(1,3) = LQ*(2d0*MUQ*h1q - (B+NUQ+MUP)*h2q)
      H(2,3) = LQ*(2d0*MUQ*h2q - (B+NUQ+MUP)*h1q)
*   Diagonalization
      CALL DIAGN(3,H,MH,VEC3,EPS)
      CALL SORTN(3,MH,VEC3)
      DO I= 1,3
       DO J= 1,3
        SCOMP(I,J)= VEC3(J,I)
       ENDDO
      ENDDO

*   CP-odd Higgs mass matrix including
      A(1,1)= (MUQ*B+GMCOMB)*(TANBQ+1D0/TANBQ)
      A(2,2)= (LQ**2*h1q*h2q/MUQ*(B+3d0*NUQ+MUP)-3d0*AKQ*NUQ
     .       -XIF*(4d0*KQ+LQ*MUP/MUQ)-2d0*MSP-MUP*NUQ-LQ*XIS/MUQ)
      A(1,2)= LQ*(ALQ+2d0*rt*Atp*fmt+2d0*rb*Abp*fmb)
     .       *DSQRT(h1q**2+h2q**2)
*   Diagonalization
      CALL DIAGN(2,A,MA,VEC2,EPS)
      CALL SORTN(2,MA,VEC2)
      DO I= 1,2
       DO J= 1,2
        PCOMP(I,J)= VEC2(J,I)
       ENDDO
      ENDDO
*   CP odd mixing angles
      DO I=1,2
       P(I,1)= PCOMP(I,1)*CB
       P(I,2)= PCOMP(I,1)*SB
       P(I,3)= PCOMP(I,2)
      ENDDO

*   Charged Higgs mass
      CMASS=((G2Q/2d0-LQ**2)*h1q*h2q+MUQ*B+GMCOMB)*(TANBQ+1D0/TANBQ)

      DO I=1,3
       MH(I)=MAX(MH(I),10d0)
      ENDDO
      DO I=1,2
       MA(I)=MAX(MA(I),10d0)
      ENDDO
      CMASS=MAX(CMASS,10d0)

* END estimate of Higgs masses
********************************************

      RMST(1)=MAX(RMST1,10d0)
      RMST(2)=MAX(RMST2,10d0)
      RMSB(1)=MAX(RMSB1,10d0)
      RMSB(2)=MAX(RMSB2,10d0)

      A0T1=NMA0(RMST(1),QSTSB)
      A0T2=NMA0(RMST(2),QSTSB)
      A0B1=NMA0(RMSB(1),QSTSB)
      A0B2=NMA0(RMSB(2),QSTSB)

* DEFINITIONS OF COUPLINGS AND PASSARINO-VELTMAN-FUNCTIONS
* defined in subfun.f
* CP-even neutral Higgses: J: MH(1-3)
      DO J=1,3
          A0H0(J)=NMA0(MH(J),QSTSB)
          DNUST(J)=SCOMP(J,1)**2
          DNUSB(J)=SCOMP(J,2)**2
* Couplings lambda_Higgs(J)_STOP(L,R)_STOP(1,2):
          LHNTLT(J,1)=CST*HRULUL(J)+SST*HRULUR(J)
          LHNTLT(J,2)=-SST*HRULUL(J)+CST*HRULUR(J)
          LHNTRT(J,1)=CST*HRULUR(J)+SST*HRURUR(J)
          LHNTRT(J,2)=-SST*HRULUR(J)+CST*HRURUR(J)
* Couplings lambda_Higgs(J)_SBOT(L,R)_SBOT(1,2):
          LHNBLB(J,1)=CSB*HRDLDL(J)+SSB*HRDLDR(J)
          LHNBLB(J,2)=-SSB*HRDLDL(J)+CSB*HRDLDR(J)
          LHNBRB(J,1)=CSB*HRDLDR(J)+SSB*HRDRDR(J)
          LHNBRB(J,2)=-SSB*HRDLDR(J)+CSB*HRDRDR(J)
* I: stop1/stop2 or sbot1/sbot2
        DO I=1,2
          B011HNTI(J,I)=NMB0(RMST(1),MH(J),RMST(I),QSTSB)
          B022HNTI(J,I)=NMB0(RMST(2),MH(J),RMST(I),QSTSB)
          B011HNBI(J,I)=NMB0(RMSB(1),MH(J),RMSB(I),QSTSB)
          B022HNBI(J,I)=NMB0(RMSB(2),MH(J),RMSB(I),QSTSB)
        ENDDO
      ENDDO
* END CP-even neutral Higgses
* CP-odd neutral Higgses:  J: MA(1-2)
      DO J=1,2
          A0H0(3+J)=NMA0(MA(J),QSTSB)
          DNUST(3+J)=P(J,1)**2
          DNUSB(3+J)=P(J,2)**2
* Couplings lambda_Higgs(J)_STOP(L,R)_STOP(1,2):
          LHNTLT(3+J,1)=SST*HIULUR(J)
          LHNTLT(3+J,2)=CST*HIULUR(J)
          LHNTRT(3+J,1)=CST*HIULUR(J)
          LHNTRT(3+J,2)=-SST*HIULUR(J)
* Couplings lambda_Higgs(J)_SBOT(L,R)_SBOT(1,2):
          LHNBLB(3+J,1)=SSB*HIDLDR(J)
          LHNBLB(3+J,2)=CSB*HIDLDR(J)
          LHNBRB(3+J,1)=CSB*HIDLDR(J)
          LHNBRB(3+J,2)=-SSB*HIDLDR(J)
* I: stop1/stop2 or sbot1/sbot2
        DO I=1,2
          B011HNTI(3+J,I)=NMB0(RMST(1),MA(J),RMST(I),QSTSB)
          B022HNTI(3+J,I)=NMB0(RMST(2),MA(J),RMST(I),QSTSB)
          B011HNBI(3+J,I)=NMB0(RMSB(1),MA(J),RMSB(I),QSTSB)
          B022HNBI(3+J,I)=NMB0(RMSB(2),MA(J),RMSB(I),QSTSB)
        ENDDO
      ENDDO
* END CP-odd neutral Higgses
* Couplings of the neutral Goldstone:
      HGBULUR=HTQ/SQR2*(ATP*SB-MUQ*CB)
      HGBDLDR= HBQ/SQR2*(MUQ*SB-ABP*CB)
* J=3: Goldstone Boson
          A0H0(3+3)=NMA0(MZ**2,QSTSB)
          DNUST(3+3)=SB**2
          DNUSB(3+3)=CB**2
* Couplings lambda_Higgs(J)_STOP(L,R)_STOP(1,2):
          LHNTLT(3+3,1)=SST*HGBULUR
          LHNTLT(3+3,2)=CST*HGBULUR
          LHNTRT(3+3,1)=CST*HGBULUR
          LHNTRT(3+3,2)=-SST*HGBULUR
* Couplings lambda_Higgs(J)_SBOT(L,R)_SBOT(1,2):
          LHNBLB(3+3,1)=SSB*HGBDLDR
          LHNBLB(3+3,2)=CSB*HGBDLDR
          LHNBRB(3+3,1)=CSB*HGBDLDR
          LHNBRB(3+3,2)=-SSB*HGBDLDR
* I: stop1/stop2 or sbot1/sbot2
        DO I=1,2
          B011HNTI(3+3,I)=NMB0(RMST(1),MZ**2,RMST(I),QSTSB)
          B022HNTI(3+3,I)=NMB0(RMST(2),MZ**2,RMST(I),QSTSB)
          B011HNBI(3+3,I)=NMB0(RMSB(1),MZ**2,RMSB(I),QSTSB)
          B022HNBI(3+3,I)=NMB0(RMSB(2),MZ**2,RMSB(I),QSTSB)
        ENDDO
* END Goldstone Boson
* Charged Higgs:
        A0HC(1) = NMA0(CMASS,QSTSB)
        DNDST(1)=CB**2
        DNDSB(1)=SB**2
* Couplings lambda_CHiggs_STOP(L,R)_SBOT(L/R):
        LHCTLBL(1)=-HTQ**2*H1Q*CB-HBQ**2*H2Q*SB
        LHCTRBR(1)=-HTQ*HBQ*(H2Q*CB+H1Q*SB)
        LHCTLBR(1)=-HBQ*(MUQ*CB+ABP*SB)
        LHCTRBL(1)=-HTQ*(MUQ*SB+ATP*CB)
* Couplings lambda_CHiggs_SBOT(L,R)_STOP(L/R):
        LHCBLTL(1)=-HTQ**2*H1Q*CB-HBQ**2*H2Q*SB
        LHCBRTR(1)=-HTQ*HBQ*(H2Q*CB+H1Q*SB)
        LHCBLTR(1)=-HTQ*(MUQ*SB+ATP*CB)
        LHCBRTL(1)=-HBQ*(MUQ*CB+ABP*SB)
* LHCTLB(i,j): i=1: charged Higgs; j: sbottom1/2
        LHCTLB(1,1)= CSB*LHCTLBL(1)+SSB*LHCTLBR(1)
        LHCTLB(1,2)=-SSB*LHCTLBL(1)+CSB*LHCTLBR(1)
        LHCTRB(1,1)= CSB*LHCTRBL(1)+SSB*LHCTRBR(1)
        LHCTRB(1,2)=-SSB*LHCTRBL(1)+CSB*LHCTRBR(1)
* LHCBLT(i,j): i=1: charged Higgs; j: stop1/2
        LHCBLT(1,1)= CST*LHCBLTL(1)+SST*LHCBLTR(1)
        LHCBLT(1,2)=-SST*LHCBLTL(1)+CST*LHCBLTR(1)
        LHCBRT(1,1)= CST*LHCBRTL(1)+SST*LHCBRTR(1)
        LHCBRT(1,2)=-SST*LHCBRTL(1)+CST*LHCBRTR(1)
* I: sbot1/sbot2 or stop1/stop2
      DO I=1,2
        B011HCBI(1,I)=NMB0(RMST(1),CMASS,RMSB(I),QSTSB)
        B022HCBI(1,I)=NMB0(RMST(2),CMASS,RMSB(I),QSTSB)
        B011HCTI(1,I)=NMB0(RMSB(1),CMASS,RMST(I),QSTSB)
        B022HCTI(1,I)=NMB0(RMSB(2),CMASS,RMST(I),QSTSB)
      ENDDO
* Charged Goldstone:
        A0HC(2) = NMA0(MW**2,QSTSB)
        DNDST(2)=SB**2
        DNDSB(2)=CB**2
        LHCTLBL(2)=-HTQ**2*H1Q*SB+HBQ**2*H2Q*CB
        LHCTRBR(2)=0D0
        LHCTLBR(2)=HBQ*(-MUQ*SB+ABP*CB)
        LHCTRBL(2)=HTQ*(MUQ*CB-ATP*SB)
        LHCBLTL(2)=-HTQ**2*H1Q*SB+HBQ**2*H2Q*CB
        LHCBRTR(2)=0D0
        LHCBLTR(2)=HTQ*( MUQ*CB-ATP*SB)
        LHCBRTL(2)=HBQ*(-MUQ*SB+ABP*CB)
* LHCTLB(i,j): i=2: charged Goldstone; j: sbottom1/2
        LHCTLB(2,1)=CSB*LHCTLBL(2)+SSB*LHCTLBR(2)
        LHCTLB(2,2)=-SSB*LHCTLBL(2)+CSB*LHCTLBR(2)
        LHCTRB(2,1)=CSB*LHCTRBL(2)+SSB*LHCTRBR(2)
        LHCTRB(2,2)=-SSB*LHCTRBL(2)+CSB*LHCTRBR(2)
* LHCBLT(i,j): i=2: charged Goldstone; j: stop1/2
        LHCBLT(2,1)= CST*LHCBLTL(2)+SST*LHCBLTR(2)
        LHCBLT(2,2)=-SST*LHCBLTL(2)+CST*LHCBLTR(2)
        LHCBRT(2,1)= CST*LHCBRTL(2)+SST*LHCBRTR(2)
        LHCBRT(2,2)=-SST*LHCBRTL(2)+CST*LHCBRTR(2)
* I: sbot1/sbot2 or stop1/stop2
      DO I=1,2
        B011HCBI(2,I)=NMB0(RMST(1),MW**2,RMSB(I),QSTSB)
        B022HCBI(2,I)=NMB0(RMST(2),MW**2,RMSB(I),QSTSB)
        B011HCTI(2,I)=NMB0(RMSB(1),MW**2,RMST(I),QSTSB)
        B022HCTI(2,I)=NMB0(RMSB(2),MW**2,RMST(I),QSTSB)
      ENDDO
* END Charged Goldstone

* CONTRIBUTIONS TO DPI_11ST:
      DPI11ST=HTQ**2*(2D0*(CST*SST)**2*A0T1+(CST**4+SST**4)*A0T2)
     .       +((HBQ**2*(CST*SSB)**2+HTQ**2*(SST*CSB)**2)*A0B1
     .        +(HBQ**2*(CST*CSB)**2+HTQ**2*(SST*SSB)**2)*A0B2)
     .       +6d0*HTQ**2*CST**2*SST**2*(A0T1-A0T2)
*  N=1,2,3: CP-EVEN; I=4,5: CP-ODD; I=6: GOLDSTONE
      DO N=1,6
        DPI11ST=DPI11ST+HTQ**2/2D0*DNUST(N)*A0H0(N)

*   I: STOP1/STOP2
        DO I=1,2
        If(N.le.3)then
          DPI11ST=DPI11ST+(CST*LHNTLT(N,I)+SST*LHNTRT(N,I))**2
     .                   *B011HNTI(N,I)
        else
          DPI11ST=DPI11ST+(CST*LHNTLT(N,I)-SST*LHNTRT(N,I))**2
     .                   *B011HNTI(N,I)
        endif
        ENDDO
      ENDDO
*  N=1: CHARGED HIGGS; N=2: CHARGED GOLDSTONE
      DO N=1,2
        DPI11ST=DPI11ST
     . +(HBQ**2*CST**2*DNDST(N)+HTQ**2*SST**2*DNUST(N))*A0HC(N)
*   I: SBOT1/SBOT2
        DO I=1,2
          DPI11ST=DPI11ST+(CST*LHCTLB(N,I)+SST*LHCTRB(N,I))**2
     .                   *B011HCBI(N,I)
        ENDDO
      ENDDO

        DPI11ST=DPI11ST+HTQ**2
     . *((RMST(1)-muq**2-mtq**2)*NMB0(RMST(1),muq**2,mtq**2,QSTSB)
     .                       -NMA0(muq**2,QSTSB)-NMA0(mtq**2,QSTSB))
     . +(HBQ**2*CST**2+HTQ**2*SST**2)
     . *((RMST(1)-muq**2-mbq**2)*NMB0(RMST(1),muq**2,mbq**2,QSTSB)
     .                       -NMA0(muq**2,QSTSB)-NMA0(mbq**2,QSTSB))
     . -4d0*CST*SST
     .     *HTQ*HBQ*muq*mbq*NMB0(RMST(1),muq**2,mbq**2,QSTSB)
      DPI11ST=DPI11ST*COEF

* CONTRIBUTIONS TO DPI_11SB:
      DPI11SB=HBQ**2*(2D0*(CSB*SSB)**2*A0B1+(CSB**4+SSB**4)*A0B2)
     .       +((HTQ**2*(CSB*SST)**2+HBQ**2*(SSB*CST)**2)*A0T1
     .        +(HTQ**2*(CSB*CST)**2+HBQ**2*(SSB*SST)**2)*A0T2)
     .       +6d0*HBQ**2*CSB**2*SSB**2*(A0B1-A0B2)
*  N=1,2,3: CP-EVEN; I=4,5: CP-ODD; I=6: GOLDSTONE
      DO N=1,6
        DPI11SB=DPI11SB+HBQ**2/2D0*DNUSB(N)*A0H0(N)
*   I: SBOT1/SBOT2
        DO I=1,2
        If(N.le.3)then
          DPI11SB=DPI11SB+(CSB*LHNBLB(N,I)+SSB*LHNBRB(N,I))**2
     .                   *B011HNBI(N,I)
        else
          DPI11SB=DPI11SB+(CSB*LHNBLB(N,I)-SSB*LHNBRB(N,I))**2
     .                   *B011HNBI(N,I)
        endif
        ENDDO
      ENDDO
*  N=1: CHARGED HIGGS; N=2: CHARGED GOLDSTONE
      DO N=1,2
        DPI11SB=DPI11SB
     .      +(HTQ**2*CSB**2*DNUSB(N)+HBQ**2*SSB**2*DNDSB(N))*A0HC(N)
*   I: STOP1/STOP2
        DO I=1,2
          DPI11SB=DPI11SB+(CSB*LHCBLT(N,I)+SSB*LHCBRT(N,I))**2
     .                   *B011HCTI(N,I)
        ENDDO
      ENDDO

      DPI11SB=DPI11SB+HBQ**2
     . *((RMSB(1)-muq**2-mbq**2)*NMB0(RMSB(1),muq**2,mbq**2,QSTSB)
     .                       -NMA0(muq**2,QSTSB)-NMA0(mbq**2,QSTSB))
     . +(HTQ**2*CSB**2+HBQ**2*SSB**2)
     . *((RMSB(1)-muq**2-mtq**2)*NMB0(RMSB(1),muq**2,mtq**2,QSTSB)
     .                       -NMA0(muq**2,QSTSB)-NMA0(mtq**2,QSTSB))
     . -4d0*CSB*SSB
     .     *HTQ*HBQ*muq*mtq*NMB0(RMSB(1),muq**2,mtq**2,QSTSB)
      DPI11SB=DPI11SB*COEF

* CONTRIBUTIONS TO DPI_22ST:
      DPI22ST=HTQ**2*(2D0*(CST*SST)**2*A0T2+(CST**4+SST**4)*A0T1)
     .       +((HTQ**2*(CST*SSB)**2+HBQ**2*(SST*CSB)**2)*A0B2
     .        +(HTQ**2*(CST*CSB)**2+HBQ**2*(SST*SSB)**2)*A0B1)
     .       +6d0*HTQ**2*CST**2*SST**2*(A0T2-A0T1)
*  N=1,2,3: CP-EVEN; I=4,5: CP-ODD; I=6: GOLDSTONE
      DO N=1,6
        DPI22ST=DPI22ST+HTQ**2/2D0*DNUST(N)*A0H0(N)
*   I: STOP1/STOP2
        DO I=1,2
        if(N.le.3)then
          DPI22ST=DPI22ST+(-SST*LHNTLT(N,I)+CST*LHNTRT(N,I))**2
     .                   *B022HNTI(N,I)
        else
          DPI22ST=DPI22ST+(-SST*LHNTLT(N,I)-CST*LHNTRT(N,I))**2
     .                   *B022HNTI(N,I)
        endif
        ENDDO
      ENDDO
*  N=1: CHARGED HIGGS; N=2: CHARGED GOLDSTONE
      DO N=1,2
        DPI22ST=DPI22ST
     .   +(HBQ**2*SST**2*DNDST(N)+HTQ**2*CST**2*DNUST(N))*A0HC(N)
*   I: SBOT1/SBOT2
        DO I=1,2
          DPI22ST=DPI22ST+(-SST*LHCTLB(N,I)+CST*LHCTRB(N,I))**2
     .                   *B022HCBI(N,I)
        ENDDO
      ENDDO

        DPI22ST=DPI22ST+HTQ**2
     . *((RMST(2)-muq**2-mtq**2)*NMB0(RMST(2),muq**2,mtq**2,QSTSB)
     .                       -NMA0(muq**2,QSTSB)-NMA0(mtq**2,QSTSB))
     . +(HBQ**2*SST**2+HTQ**2*CST**2)
     . *((RMST(2)-muq**2-mbq**2)*NMB0(RMST(2),muq**2,mbq**2,QSTSB)
     .                       -NMA0(muq**2,QSTSB)-NMA0(mbq**2,QSTSB))
     . +4d0*CST*SST
     .     *HTQ*HBQ*muq*mbq*NMB0(RMST(2),muq**2,mbq**2,QSTSB)
      DPI22ST=DPI22ST*COEF

* CONTRIBUTIONS TO DPI_22SB:
      DPI22SB=HBQ**2*(2D0*(CSB*SSB)**2*A0B2+(CSB**4+SSB**4)*A0B1)
     .       +((HBQ**2*(CSB*SST)**2+HTQ**2*(SSB*CST)**2)*A0T2
     .        +(HBQ**2*(CSB*CST)**2+HTQ**2*(SSB*SST)**2)*A0T1)
     .       +6d0*HBQ**2*CSB**2*SSB**2*(A0B2-A0B1)
*  N=1,2,3: CP-EVEN; I=4,5: CP-ODD; I=6: GOLDSTONE
      DO N=1,6
        DPI22SB=DPI22SB+HBQ**2/2D0*DNUSB(N)*A0H0(N)
*   I: SBOT1/SBOT2
        DO I=1,2
        if(N.le.3)then
          DPI22SB=DPI22SB+(-SSB*LHNBLB(N,I)+CSB*LHNBRB(N,I))**2
     .                   *B022HNBI(N,I)
        else
          DPI22SB=DPI22SB+(-SSB*LHNBLB(N,I)-CSB*LHNBRB(N,I))**2
     .                   *B022HNBI(N,I)
        endif
        ENDDO
      ENDDO
*  N=1: CHARGED HIGGS; N=2: CHARGED GOLDSTONE
      DO N=1,2
        DPI22SB=DPI22SB
     .    +(HTQ**2*SSB**2*DNUSB(N)+HBQ**2*CSB**2*DNDSB(N))*A0HC(N)
*   I: STOP1/STOP2
        DO I=1,2
          DPI22SB=DPI22SB+(-SSB*LHCBLT(N,I)+CSB*LHCBRT(N,I))**2
     .                   *B022HCTI(N,I)
        ENDDO
      ENDDO

      DPI22SB=DPI22SB+HBQ**2
     . *((RMSB(2)-muq**2-mbq**2)*NMB0(RMSB(2),muq**2,mbq**2,QSTSB)
     .                       -NMA0(muq**2,QSTSB)-NMA0(mbq**2,QSTSB))
     . +(HTQ**2*SSB**2+HBQ**2*CSB**2)
     . *((RMSB(2)-muq**2-mtq**2)*NMB0(RMSB(2),muq**2,mtq**2,QSTSB)
     .                       -NMA0(muq**2,QSTSB)-NMA0(mtq**2,QSTSB))
     . +4d0*CSB*SSB
     .     *HTQ*HBQ*muq*mtq*NMB0(RMSB(2),muq**2,mtq**2,QSTSB)
      DPI22SB=DPI22SB*COEF

      END

