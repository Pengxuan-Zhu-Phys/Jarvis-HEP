      SUBROUTINE HDECAY_lightA(PAR,PROB)

      IMPLICIT NONE

      INTEGER K,L,INDA,INDPI,INDETA,INDETAP,INDEB1,MIXPROB,N0,NF

      DOUBLE PRECISION PAR(*),PROB(*),MA,GACHACHA(2,2),EPS,PI,aux
      DOUBLE PRECISION FPI,THETAE,MPI,MPIC,MPI8,MPI9,META,METAP,MKc,MK0
      DOUBLE PRECISION DM3,DM8,DM9,MMIX(4,4),VALP(4),VECP(4,4)
     .                 ,OMIX(4,4)
      DOUBLE PRECISION WidthPiinv,WidthEinv,WidthEPinv,PROBpiinv
      DOUBLE PRECISION CGPI,CGETA,CGETAP,WidthPigaga,WidthEgaga,
     .                WidthEPgaga,WidthECgaga,PROBpigaga
      DOUBLE PRECISION YUE,YDE,YSE
      DOUBLE PRECISION X,SP,BETA,QCd0,AQCDM,AQCD,QCDA,TQCDA
      DOUBLE PRECISION ARPill,AIPill,Rchirho,Ichirho
      DOUBLE PRECISION RMS,ASH,AS3,AS4,RATCOUP,HIGTOP,ALPHAS,RUNM
      DOUBLE PRECISION AAee,APIee,AEee,AEPee,AIPIee,AIEee,AIEPee
      DOUBLE PRECISION AAmumu,APImumu,AEmumu,AEPmumu
      DOUBLE PRECISION AAtata,AIPImumu,AIEmumu,AIEPmumu
      DOUBLE PRECISION PROBpill,WidthPiee,WidthEee,WidthEPee,WidthEmumu
      DOUBLE PRECISION WidthEB1mumu,WidthEB1tata
      DOUBLE PRECISION AA3Pi,API3Pi,AE3Pi,AEP3Pi,GamA3Pi,KinA3Pi
      DOUBLE PRECISION AAPi3PiC,APIPi3PiC,AEPi3PiC,AEPPi3PiC,GamAPi3PiC
     .                 ,KinAPi3PiC
      DOUBLE PRECISION AAEPi3,APIEPi3,AEEPi3,AEPEPi3,GamAEPi3
     .                 ,KinAEPi3
      DOUBLE PRECISION AAEPiC,APIEPiC,AEEPiC,AEPEPiC,GamAEPiC
     .                 ,KinAEPiC
      DOUBLE PRECISION AAEPPi3,APIEPPi3,AEEPPi3,AEPEPPi3,GamAEPPi3
     .                 ,KinAEPPi3
      DOUBLE PRECISION AAEPPiC,APIEPPiC,AEEPPiC,AEPEPPiC,GamAEPPiC
     .                 ,KinAEPPiC
      DOUBLE PRECISION AAPiEE,APIPiEE,AEPiEE,AEPPiEE,GamAPiEE,KinAPiEE
      DOUBLE PRECISION AAPiEEP,APIPiEEP,AEPiEEP,AEPPiEEP,GamAPiEEP
     .                 ,KinAPiEEP
      DOUBLE PRECISION AAPiEPEP,APIPiEPEP,AEPiEPEP,AEPPiEPEP,GamAPiEPEP
     .                 ,KinAPiEPEP
      DOUBLE PRECISION AA3E,API3E,AE3E,AEP3E,GamA3E,KinA3E
      DOUBLE PRECISION AAE2EP,APIE2EP,AEE2EP,AEPE2EP,GamAE2EP,KinAE2EP
      DOUBLE PRECISION AAEEp2,APIEEp2,AEEEp2,AEPEEp2,GamAEEp2,KinAEEp2
      DOUBLE PRECISION AA3EP,API3EP,AE3EP,AEP3EP,GamA3EP,KinA3EP
      DOUBLE PRECISION AAPiKC,APIPiKC,AEPiKC,AEPPiKC,GamAPiKC,KinAPiKC
      DOUBLE PRECISION AAPiK0,APIPiK0,AEPiK0,AEPPiK0,GamAPiK0,KinAPiK0
      DOUBLE PRECISION AAPiKCK0,APIPiKCK0,AEPiKCK0,AEPPiKCK0,GamAPiKCK0
     .                 ,KinAPiKCK0
      DOUBLE PRECISION AAEKC,APIEKC,AEEKC,AEPEKC,GamAEKC,KinAEKC
      DOUBLE PRECISION AAEK0,APIEK0,AEEK0,AEPEK0,GamAEK0,KinAEK0
      DOUBLE PRECISION AAEPKC,APIEPKC,AEEPKC,AEPEPKC,GamAEPKC,KinAEPKC
      DOUBLE PRECISION AAEPK0,APIEPK0,AEEPK0,AEPEPK0,GamAEPK0,KinAEPK0
      DOUBLE PRECISION AARhogam,APIRhogam,AERhogam,AEPRhogam,GamARhogam
     .                 ,KinARhogam
      DOUBLE PRECISION AAuu,AAdd,AAss,GamAuu,GamAdd,GamAss,GamAjj
      DOUBLE PRECISION METAC,DMC,CMIX,MD,RMC,DCC,KinAPiDD
      DOUBLE PRECISION METAB1,METAB2,METAB3,DMB1,DMB2,DMB3
      DOUBLE PRECISION MBm,RMB,runmb,DBB,KinAPiBB
      DOUBLE PRECISION ALEM0
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION LAMBDA,KAPPA,ALAMBDA,AKAPPA,MUEFF,NUEFF
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MCC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MPI0,MEL,MSTRANGE
      DOUBLE PRECISION XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION GamAGAGA,GamAee,GamAmumu,GamAtata,GamAhadr,
     . GamAcc,GamAbb,GamAinv(3)
      DOUBLE PRECISION DDCOS,DDSIN
* New May 2019:
      DOUBLE PRECISION ZETA2,ZETA3,ASG,AGGQCD2
* End New

      DOUBLE COMPLEX XC,TC,BC,CC,LC,MC,EC,CH1C,CH2C,CJA,CGA,FGG,CM,LI2

      COMMON/ALEM0/ALEM0
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL
      COMMON/QPAR/LAMBDA,KAPPA,ALAMBDA,AKAPPA,MUEFF,NUEFF
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MCC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SMEXT/MPI0,MEL,MSTRANGE
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/LIGHTADECAYS/GamAGAGA,GamAee,GamAmumu,GamAtata,GamAhadr,
     . GamAcc,GamAbb,GamAinv

      CM(X)= DCMPLX(MIN(1d3,X)**2,-EPS/4d0)
      FGG(XC)=2d0*XC
     .        *CDLOG((CDSQRT(1d0-4d0*XC)-1d0)
     .               /(CDSQRT(1d0-4d0*XC)+1d0))**2
      BETA(X)= DSQRT(1d0-4d0*X)
      QCd0(X)= (1d0+X**2)*(4d0*SP((1d0-X)/(1d0+X))
     . +2d0*SP((X-1d0)/(X+1d0))
     . - 3d0*DLOG((1d0+X)/(1d0-X))*DLOG(2d0/(1d0+X))
     . - 2d0*DLOG((1d0+X)/(1d0-X))*DLOG(X))
     . - 3d0*X*DLOG(4d0/(1d0-X**2))-4d0*X*DLOG(X)
      AQCDM(X)= QCd0(X)/X+(19d0+2d0*X**2+3d0*X**4)/16d0/X
     . * DLOG((1d0+X)/(1d0-X))+3d0/8d0*(7d0-X**2)
* July 2010:
c      AQCD(X)=(4d0/3d0*AQCDM(BETA(X))
c     .   +2d0*(4d0/3d0-DLOG(X))*(1d0-6d0*X)/(1d0-4d0*X))*ASH/PI
c     .   + (29.14671d0 + RATCOUP*(23d0/6d0 - DLOG(HIGTOP)
c     .   + DLOG(X)**2/6d0))*(ASH/PI)**2
c     .   + (164.14d0 - 25.77d0*5 + 0.259d0*5**2)*(ASH/PI)**3
* New May 2019:
      AQCD(X)=(4d0/3d0*AQCDM(BETA(X))
     .       + 2d0*(4d0/3d0-DLOG(X))*(1d0-6d0*X)/(1d0-4d0*X))*ASH/PI
     .       + (29.14671d0 + RATCOUP*(23d0/6d0 - DLOG(HIGTOP)
     .         + DLOG(X)**2/6d0))*(ASH/PI)**2
     .       + (164.14d0 - 25.77d0*5d0 + 0.259d0*5d0**2)*(ASH/PI)**3
     .       + (39.34d0-220.9d0*5d0+9.685d0*5d0**2
     .         - 0.0205d0*5d0**3)*(ASH/PI)**4
* End New
      QCDA(X)= 1d0+AQCD(X)
      TQCDA(X)= 1d0+4d0/3d0*AQCDM(BETA(X))*ASH/PI
      ARPill(X)=(dlog((1d0-X)/(1d0+X))**2+4d0*SP(-(1d0-X)/(1d0+X))
     .   +Pi**2/3d0)/(4d0*X)-5d0/2d0+Rchirho
      AIPill(X)=dimag(LI2(DCMPLX(-(1d0-X)/(1d0+X),0d0)))/X
     .   +Ichirho+Pi/(2d0*X)*dlog((1d0-X)/(1d0+X))
* New May 2019:
      AGGQCD2(ASG,NF,MA,MT)=1d0+ASG/PI*(97d0/4d0-NF*7d0/6d0)
     . +(ASG/PI)**2*(237311d0/864d0-529d0/24d0*ZETA2-445d0/8d0*ZETA3
     . +5d0*DLOG(MA**2/MT**2))
* End New

      EPS=1d-8
      PI=4d0*DATAN(1d0)
      PROBpiinv=0d0
      PROBpigaga=0d0
      PROBpill=0d0
* New May 2019:
      ZETA2=PI**2/6d0
      ZETA3=1.202056903159594d0
* End New

      MA=PMASS(1)

      FPI=0.093d0
      THETAE=-13d0*Pi/180d0
      MPI=0.135d0
      MPIC=0.1396d0
      META=0.548d0
      METAP=0.958d0
      MKc=0.494d0
      MK0=0.498d0
      MPI8=dsqrt(META**2*DDCOS(THETAE)**2+METAP**2*DDSIN(THETAE)**2)
      MPI9=dsqrt(META**2*DDSIN(THETAE)**2+METAP**2*DDCOS(THETAE)**2)

C       Coupling to photons / gluons

      TC=CM(MT/MA)
      BC=CM(MBP/MA)
      CC=CM(MCC/MA)
      LC=CM(MTAU/MA)
      MC=CM(MMUON/MA)
      EC=CM(MEL/MA)
      CH1C=CM(MCHA(1)/MA)
      CH2C=CM(MCHA(2)/MA)

      DO K=1,2
      DO L=1,2
       GACHACHA(K,L)=LAMBDA/dsqrt(2d0)*P2(1,2)*U(K,2)*V(L,2)
     .    -dsqrt(g2/2d0)*P2(1,1)/dsqrt(1d0+PAR(3)**2)
     .           *(U(K,1)*V(L,2)+PAR(3)*U(K,2)*V(L,1))
      ENDDO
      ENDDO

      CJA=-DSQRT(dsqrt(2d0)*GF)/4d0
     .         *(CU(4)*(FGG(TC)+FGG(CC))+CB(4)*FGG(BC))

      CGA=-DSQRT(dsqrt(2d0)*GF)/2d0*(4d0/3d0*CU(4)*(FGG(TC)+FGG(CC))
     .            +CB(4)/3d0*FGG(BC)+CD(4)*(FGG(LC)+FGG(MC)+FGG(EC)))
     . -GACHACHA(1,1)/MCHA(1)*FGG(CH1C)/2d0
     . -GACHACHA(2,2)/MCHA(2)*FGG(CH2C)/2d0

C       Initializing diphoton and leptonic decays

C  * A -> 2gamma
      GamAGAGA=8d0*CDABS(CGA)**2*MA**3
     .          /(32d0*Pi)*(ALEM0/4d0/Pi)**2

C  * A -> ee
      AAee=MEL*CD(4)*DSQRT(dsqrt(2d0)*GF)
      GamAee=AAee**2/(8d0*Pi)*MA*dsqrt(max(0d0,1d0-4d0*MEL**2/MA**2))

C  * A -> mumu
      AAmumu=MMUON*CD(4)*DSQRT(dsqrt(2d0)*GF)
      GamAmumu=AAmumu**2/(8d0*Pi)*MA
     .               *dsqrt(max(0d0,1d0-4d0*MMUON**2/MA**2))

C  * A -> tautau
      AAtata=MTAU*CL(4)*DSQRT(dsqrt(2d0)*GF)
      GamAtata=AAtata**2/(8d0*Pi)*MA
     .               *dsqrt(max(0d0,1d0-4d0*MTAU**2/MA**2))

C       Decay to light neutralinos

      IF(MA.LE.2d0*DABS(MNEU(1)))THEN
       GamAinv(1)=0d0
      ELSE
       aux=dsqrt(2d0)*LAMBDA*(P2(1,1)/dsqrt(1d0+PAR(3)**2)
     .                        *(N(1,4)*N(1,5)+PAR(3)*N(1,3)*N(1,5))
     .                        +P2(1,2)*N(1,3)*N(1,4))
     .    -dsqrt(2d0)*KAPPA*P2(1,2)*N(1,5)*N(1,5)
     .    -dsqrt(g1)*P2(1,1)/dsqrt(1d0+PAR(3)**2)
     .              *(-N(1,1)*N(1,3)+PAR(3)*N(1,1)*N(1,4))
     .    -dsqrt(g2)*P2(1,1)/dsqrt(1d0+PAR(3)**2)
     .              *(N(1,2)*N(1,3)-PAR(3)*N(1,2)*N(1,4))

       GamAinv(1)=aux**2/(16d0*PI)*MA*dsqrt(1d0-4d0*(MNEU(1)/MA)**2)
      ENDIF

      IF(MA.LE.DABS(MNEU(1))+DABS(MNEU(2)))THEN
       GamAinv(2)=0d0
      ELSE
       aux=LAMBDA/dsqrt(2d0)*(P2(1,1)/dsqrt(1d0+PAR(3)**2)
     .        *(N(1,4)*N(2,5)+N(2,4)*N(1,5)
     .         +PAR(3)*(N(1,3)*N(2,5)+N(2,3)*N(1,5)))
     .   +P2(1,2)*(N(1,3)*N(2,4)+N(2,3)*N(1,4)))
     .    -dsqrt(2d0)*KAPPA*P2(1,2)*N(1,5)*N(2,5)
     .    -dsqrt(g1)/2d0*P2(1,1)/dsqrt(1d0+PAR(3)**2)
     .            *(-N(1,1)*N(2,3)-N(2,1)*N(1,3)
     .              +PAR(3)*(N(1,1)*N(2,4)+N(2,1)*N(1,4)))
     .    -dsqrt(g2)/2d0*P2(1,1)/dsqrt(1d0+PAR(3)**2)
     .            *(N(1,2)*N(2,3)+N(2,2)*N(1,3)
     .              -PAR(3)*(N(1,2)*N(2,4)+N(2,2)*N(3,4)))

       GamAinv(2)=aux**2/(8d0*PI)*MA
     .             *dsqrt(1d0-((DABS(MNEU(1))+DABS(MNEU(2)))/MA)**2)
      ENDIF

      IF(MA.LE.2d0*DABS(MNEU(2)))THEN
       GamAinv(3)=0d0
      ELSE
       aux=dsqrt(2d0)*LAMBDA*(P2(1,1)/dsqrt(1d0+PAR(3)**2)
     .                        *(N(2,4)*N(2,5)+PAR(3)*N(2,3)*N(2,5))
     .                        +P2(1,2)*N(2,3)*N(2,4))
     .    -dsqrt(2d0)*KAPPA*P2(1,2)*N(2,5)*N(2,5)
     .    -dsqrt(g1)*P2(1,1)/dsqrt(1d0+PAR(3)**2)
     .              *(-N(2,1)*N(2,3)+PAR(3)*N(2,1)*N(2,4))
     .    -dsqrt(g2)*P2(1,1)/dsqrt(1d0+PAR(3)**2)
     .              *(N(2,2)*N(2,3)-PAR(3)*N(2,2)*N(2,4))

       GamAinv(3)=aux**2/(16d0*PI)*MA*dsqrt(1d0-4d0*(MNEU(2)/MA)**2)
      ENDIF

C       Mixing with mesons [1612.06538[hep-ph]]

      GamAhadr=0d0

      IF(MA.lt.4d0)then

      DM3=-DSQRT(dsqrt(2d0)*GF)/2d0*FPI*MPI**2*(CU(4)-CD(4))
      DM8=-DSQRT(dsqrt(2d0)*GF)/2d0*FPI*(-dsqrt(3d0)*MPI8**2*CD(4)
     .                           +MPI**2/dsqrt(3d0)*(CU(4)+2d0*CD(4)))
      DM9=-DSQRT(dsqrt(2d0)*GF)/2d0*FPI*(dsqrt(3d0/2d0)*MPI8**2*CD(4)
     .                           +MPI**2/dsqrt(6d0)*(2d0*CU(4)+CD(4)))
     .    -dsqrt(2d0/3d0)*DREAL(CJA)*FPI*(MPI9**2-(MPI8**2+MPI**2)/2d0)

      MMIX(1,1)=MPI**2
      MMIX(1,2)=0d0
      MMIX(1,3)=0d0
      MMIX(1,4)=DM3
      MMIX(2,1)=0d0
      MMIX(2,2)=META**2
      MMIX(2,3)=0d0
      MMIX(2,4)=DDCOS(THETAE)*DM8-DDSIN(THETAE)*DM9
      MMIX(3,1)=0d0
      MMIX(3,2)=0d0
      MMIX(3,3)=METAP**2
      MMIX(3,4)=DDSIN(THETAE)*DM8+DDCOS(THETAE)*DM9
      MMIX(4,1)=DM3
      MMIX(4,2)=DDCOS(THETAE)*DM8-DDSIN(THETAE)*DM9
      MMIX(4,3)=DDSIN(THETAE)*DM8+DDCOS(THETAE)*DM9
      MMIX(4,4)=MA**2
     .    +2d0/3d0*FPI**2*(MPI9**2-(MPI8**2+MPI**2)/3d0)*CDABS(CJA)**2

      CALL DIAGN(4,MMIX,VALP,VECP,EPS)
      CALL SORTN(4,VALP,VECP)
      DO K=1,4
       DO L=1,4
        OMIX(K,L)=VECP(L,K)
       ENDDO
      ENDDO

      MIXPROB=0
      INDA=1
      DO K=1,4
       IF(dabs(OMIX(INDA,4))**2.lt.dabs(OMIX(K,4))**2)INDA=K
      ENDDO
      INDPI=Max(1,3-INDA)
      DO K=INDPI,4
       IF(dabs(OMIX(INDPI,1))**2.lt.dabs(OMIX(K,1))**2.and.K.ne.INDA)
     .   INDPI=K
      ENDDO
      IF(dabs(OMIX(INDA,1))**2.gt.dabs(OMIX(INDPI,1))**2)MIXPROB=1
      INDETA=Max(1,3-INDA)
      DO K=INDETA,4
       IF(dabs(OMIX(INDETA,2))**2.lt.dabs(OMIX(K,2))**2.and.K.ne.INDA)
     .   INDETA=K
      ENDDO
      IF(dabs(OMIX(INDA,2))**2.gt.dabs(OMIX(INDETA,2))**2)MIXPROB=1
      IF(INDETA.eq.INDPI)MIXPROB=1
      INDETAP=Max(1,3-INDA)
      DO K=INDETAP,4
       IF(dabs(OMIX(INDETAP,3))**2.lt.dabs(OMIX(K,3))**2.and.K.ne.INDA)
     .   INDETAP=K
      ENDDO
      IF(dabs(OMIX(INDA,3))**2.gt.dabs(OMIX(INDETAP,3))**2)MIXPROB=1
      IF(INDETAP.eq.INDPI)MIXPROB=1
      IF(INDETAP.eq.INDETA)MIXPROB=1


C       Invisible decay A -> 2 neu_{1,2}

      GamAinv(1)=GamAinv(1)*dabs(OMIX(INDA,4))**2
      GamAinv(2)=GamAinv(2)*dabs(OMIX(INDA,4))**2
      GamAinv(3)=GamAinv(3)*dabs(OMIX(INDA,4))**2

c Test invisible meson decays
       aux=dsqrt(2d0)*LAMBDA*(P2(1,1)/dsqrt(1d0+PAR(3)**2)
     .                        *(N(1,4)*N(1,5)+PAR(3)*N(1,3)*N(1,5))
     .                        +P2(1,2)*N(1,3)*N(1,4))
     .    -dsqrt(2d0)*KAPPA*P2(1,2)*N(1,5)*N(1,5)
     .    -dsqrt(g1)*P2(1,1)/dsqrt(1d0+PAR(3)**2)
     .              *(-N(1,1)*N(1,3)+PAR(3)*N(1,1)*N(1,4))
     .    -dsqrt(g2)*P2(1,1)/dsqrt(1d0+PAR(3)**2)
     .              *(N(1,2)*N(1,3)-PAR(3)*N(1,2)*N(1,4))
      IF(MPI.LE.2d0*DABS(MNEU(1)))THEN
       WidthPiinv=0d0
      ELSE
       WidthPiinv=aux**2/(16d0*PI)*MPI*dabs(OMIX(INDPI,4))**2
     .            *dsqrt(1d0-4d0*(MNEU(1)/MPI)**2)
      ENDIF
      PROBpiinv=Max(0d0,WidthPiinv*0.7d0/2.1d-15-1d0)
      IF(META.LE.2d0*DABS(MNEU(1)))THEN
       WidthEinv=0d0
      ELSE
       WidthEinv=aux**2/(16d0*PI)*META*dabs(OMIX(INDETA,4))**2
     .            *dsqrt(1d0-4d0*(MNEU(1)/META)**2)
      ENDIF
      PROBpiinv=Max(PROBpiinv,WidthEinv*0.7d0/3.4d-10-1d0)
      IF(METAP.LE.2d0*DABS(MNEU(1)))THEN
       WidthEPinv=0d0
      ELSE
       WidthEPinv=aux**2/(16d0*PI)*METAP*dabs(OMIX(INDETAP,4))**2
     .            *dsqrt(1d0-4d0*(MNEU(1)/METAP)**2)
      ENDIF
      PROBpiinv=Max(PROBpiinv,WidthEPinv*0.7d0/9.8d-8-1d0)


C       Diphoton decay A -> 2gamma

      CGPI=-dsqrt((6.582d-25/8.52d-17*0.98823d0)
     .            /(8d0*(ALEM0/(4d0*Pi))**2*MPI**3/(32d0*Pi)))
      CGETA=dsqrt((1.31d-6*0.3941d0)
     .            /(8d0*(ALEM0/(4d0*Pi))**2*META**3/(32d0*Pi)))
     . *(-(DDCOS(THETAE)-2d0*dsqrt(2d0)*DDSIN(THETAE))
     .    /dabs(DDCOS(THETAE)-2d0*dsqrt(2d0)*DDSIN(THETAE)))
      CGETAP=dsqrt((0.000197d0*0.0221d0)
     .            /(8d0*(ALEM0/(4d0*Pi))**2*METAP**3/(32d0*Pi)))
     . *(-(DDSIN(THETAE)+2d0*dsqrt(2d0)*DDCOS(THETAE))
     .    /dabs(DDSIN(THETAE)+2d0*dsqrt(2d0)*DDCOS(THETAE)))

      aux=OMIX(INDA,1)*CGPI+OMIX(INDA,2)*CGETA+OMIX(INDA,3)*CGETAP
     .          +OMIX(INDA,4)*DREAL(CGA)
      GamAGAGA=8d0*(aux**2+(OMIX(INDA,4)*DIMAG(CGA))**2)*MA**3
     .          /(32d0*Pi)*(ALEM0/4d0/Pi)**2

c Test meson decays into photons
      aux=(OMIX(INDPI,1)*CGPI+OMIX(INDPI,2)*CGETA+OMIX(INDPI,3)*CGETAP
     .     +OMIX(INDPI,4)*DREAL(CGA))**2+(OMIX(INDPI,4)*DIMAG(CGA))**2
      WidthPigaga=8d0*aux*MPI**3/(32d0*Pi)*(ALEM0/4d0/Pi)**2   ! Decay width of the Pi0 state
      aux=(dabs(aux-CGPI**2)*0.3d0+0.07d0*CGPI**2)*MPI**3
     .          /(32d0*Pi)*(ALEM0/4d0/Pi)**2                   ! Uncertainty
      PROBpigaga=Max(0d0,dabs(WidthPigaga-6.582d-25/8.52d-17*0.98823d0)
     .                                        /aux-1d0)

      aux=(OMIX(INDETA,1)*CGPI+OMIX(INDETA,2)*CGETA
     .     +OMIX(INDETA,3)*CGETAP+OMIX(INDETA,4)*DREAL(CGA))**2
     .    +(OMIX(INDETA,4)*DIMAG(CGA))**2
      WidthEgaga=8d0*aux*META**3/(32d0*Pi)*(ALEM0/4d0/Pi)**2   ! Decay width of the Eta state
      aux=(dabs(aux-CGETA**2)*0.3d0+0.07d0*CGETA**2)*META**3
     .          /(32d0*Pi)*(ALEM0/4d0/Pi)**2                   ! Uncertainty
      PROBpigaga=Max(PROBpigaga,dabs(WidthEgaga-1.31d-6*0.3941d0)
     .                                        /aux-1d0)

      aux=(OMIX(INDETAP,1)*CGPI+OMIX(INDETAP,2)*CGETA
     .     +OMIX(INDETAP,3)*CGETAP+OMIX(INDETAP,4)*DREAL(CGA))**2
     .    +(OMIX(INDETAP,4)*DIMAG(CGA))**2
      WidthEPgaga=8d0*aux*METAP**3/(32d0*Pi)*(ALEM0/4d0/Pi)**2 ! Decay width of the Eta state
      aux=(dabs(aux-CGETAP**2)*0.3d0+0.07d0*CGETAP**2)*METAP**3
     .          /(32d0*Pi)*(ALEM0/4d0/Pi)**2                   ! Uncertainty
      PROBpigaga=Max(PROBpigaga,dabs(WidthEPgaga-0.000197d0*0.0221d0)
     .                                        /aux-1d0)


C       Leptonic decays

C  * A -> ee
      AAee=MEL*CD(4)*DSQRT(dsqrt(2d0)*GF)
      APiee=0d0
      AIPiee=0d0
      AEee=0d0
      AIEee=0d0
      AEPee=0d0
      AIEPee=0d0
      IF(MA.ge.2d0*MEL)then
       Rchirho=3d0/2d0*dlog(MEL**2/0.77d0**2)+2.99d0
       Ichirho=0d0
       APiee=(ALEM0/Pi)**2*MEL*CGPI/2d0
     .                *ARPill(dsqrt(1d0-4d0*MEL**2/MA**2))
! dsqrt(6.582d-25/8.52d-17*6.46d-8
!     .            *8d0*Pi/(MPI*dsqrt(1d0-4d0*MEL**2/MPI**2)))
       AIPiee=(ALEM0/Pi)**2*MEL*CGPI/2d0
     .                *AIPill(dsqrt(1d0-4d0*MEL**2/MA**2))
       Rchirho=3d0/2d0*dlog(MEL**2/0.77d0**2)+6.46d0
       Ichirho=0d0
       AEee=(ALEM0/Pi)**2*MEL*CGETA/2d0
     .                *ARPill(dsqrt(1d0-4d0*MEL**2/MA**2))
       AIEee=(ALEM0/Pi)**2*MEL*CGETA/2d0
     .                *AIPill(dsqrt(1d0-4d0*MEL**2/MA**2))
       Rchirho=3d0/2d0*dlog(MEL**2/0.77d0**2)+14.9d0
       Ichirho=2.52d0
       AEPee=(ALEM0/Pi)**2*MEL*CGETAP/2d0
     .                *ARPill(dsqrt(1d0-4d0*MEL**2/MA**2))
       AIEPee=(ALEM0/Pi)**2*MEL*CGETAP/2d0
     .                *AIPill(dsqrt(1d0-4d0*MEL**2/MA**2))
      ENDIF

      aux=(OMIX(INDA,1)*APiee+OMIX(INDA,2)*AEee+OMIX(INDA,3)*AEPee
     .          +OMIX(INDA,4)*AAee)**2
     . +(OMIX(INDA,1)*AIPiee+OMIX(INDA,2)*AIEee+OMIX(INDA,3)*AIEPee)**2
      GamAee=aux/(8d0*Pi)*MA*dsqrt(max(0d0,1d0-4d0*MEL**2/MA**2))

c Test Pi0 -> e+e-:
      Rchirho=3d0/2d0*dlog(MEL**2/0.77d0**2)+2.99d0
      Ichirho=0d0
      APiee=(ALEM0/Pi)**2*MEL*CGPI/2d0
     .                *ARPill(dsqrt(1d0-4d0*MEL**2/MPI**2))
      AIPiee=(ALEM0/Pi)**2*MEL*CGPI/2d0
     .                *AIPill(dsqrt(1d0-4d0*MEL**2/MPI**2))
      Rchirho=3d0/2d0*dlog(MEL**2/0.77d0**2)+6.46d0
      Ichirho=0d0
      AEee=(ALEM0/Pi)**2*MEL*CGETA/2d0
     .                *ARPill(dsqrt(1d0-4d0*MEL**2/MPI**2))
      AIEee=(ALEM0/Pi)**2*MEL*CGETA/2d0
     .                *AIPill(dsqrt(1d0-4d0*MEL**2/MPI**2))
      Rchirho=3d0/2d0*dlog(MEL**2/0.77d0**2)+14.9d0
      Ichirho=2.52d0
      AEPee=(ALEM0/Pi)**2*MEL*CGETAP/2d0
     .                *ARPill(dsqrt(1d0-4d0*MEL**2/MPI**2))
      AIEPee=(ALEM0/Pi)**2*MEL*CGETAP/2d0
     .                *AIPill(dsqrt(1d0-4d0*MEL**2/MPI**2))
      aux=(OMIX(INDPI,1)*APiee+OMIX(INDPI,2)*AEee+OMIX(INDPI,3)*AEPee
     .          +OMIX(INDPI,4)*AAee)**2
     . +(OMIX(INDPI,1)*AIPiee+OMIX(INDPI,2)*AIEee
     .          +OMIX(INDPI,3)*AIEPee)**2
      WidthPiee=aux/(8d0*Pi)*MPI*dsqrt(max(0d0,1d0-4d0*MEL**2/MPI**2)) ! Width Pi0 -> e+e-
      aux=0.3d0*dabs(aux-APiee**2-AIPiee**2)/(8d0*Pi)*MPI              ! Mixing uncertainty
     .     *dsqrt(max(0d0,1d0-4d0*MEL**2/MPI**2))
     .   +0.05d0*(APiee**2+AIPiee**2)*MPI/(8d0*Pi)                     ! SM uncertainty             
     .          *dsqrt(max(0d0,1d0-4d0*MEL**2/MPI**2))
     .   +2d0*0.055d0*6.582d-25/8.52d-17*6.46d-8                       ! Exp uncertainty
      aux=Max(aux,1.1d0*dabs((APiee**2+AIPiee**2)*MPI/(8d0*Pi)         ! SM-Exp. Error
     .          *dsqrt(max(0d0,1d0-4d0*MEL**2/MPI**2))
     .                     -6.582d-25/8.52d-17*6.46d-8))
      PROBpill=Max(0d0,dabs(WidthPiee-6.582d-25/8.52d-17*6.46d-8)
     .                                                  /aux-1d0)

c Test Eta -> e+e-:
      Rchirho=3d0/2d0*dlog(MEL**2/0.77d0**2)+2.99d0
      Ichirho=0d0
      APiee=(ALEM0/Pi)**2*MEL*CGPI/2d0
     .                *ARPill(dsqrt(1d0-4d0*MEL**2/META**2))
      AIPiee=(ALEM0/Pi)**2*MEL*CGPI/2d0
     .                *AIPill(dsqrt(1d0-4d0*MEL**2/META**2))
      Rchirho=3d0/2d0*dlog(MEL**2/0.77d0**2)+6.46d0
      Ichirho=0d0
      AEee=(ALEM0/Pi)**2*MEL*CGETA/2d0
     .                *ARPill(dsqrt(1d0-4d0*MEL**2/META**2))
      AIEee=(ALEM0/Pi)**2*MEL*CGETA/2d0
     .                *AIPill(dsqrt(1d0-4d0*MEL**2/META**2))
      Rchirho=3d0/2d0*dlog(MEL**2/0.77d0**2)+14.9d0
      Ichirho=2.52d0
      AEPee=(ALEM0/Pi)**2*MEL*CGETAP/2d0
     .                *ARPill(dsqrt(1d0-4d0*MEL**2/META**2))
      AIEPee=(ALEM0/Pi)**2*MEL*CGETAP/2d0
     .                *AIPill(dsqrt(1d0-4d0*MEL**2/META**2))
      aux=(OMIX(INDETA,1)*APiee+OMIX(INDETA,2)*AEee
     .          +OMIX(INDETA,3)*AEPee+OMIX(INDETA,4)*AAee)**2
     . +(OMIX(INDETA,1)*AIPiee+OMIX(INDETA,2)*AIEee
     .          +OMIX(INDETA,3)*AIEPee)**2
      WidthEee=aux/(8d0*Pi)*META*dsqrt(max(0d0,1d0-4d0*MEL**2/META**2))! Width Eta -> e+e-
      aux=0.3d0*dabs(aux-AEee**2-AIEee**2)/(8d0*Pi)*META               ! Mixing uncertainty
     .     *dsqrt(max(0d0,1d0-4d0*MEL**2/META**2))
     .   +0.1d0*(AEee**2+AIEee**2)*META/(8d0*Pi)                       ! SM uncertainty             
     .          *dsqrt(max(0d0,1d0-4d0*MEL**2/META**2))
      PROBpill=Max(PROBpill,(WidthPiee-aux)/(3.013d-12)-1d0)

c Test Eta' -> e+e-:
      Rchirho=3d0/2d0*dlog(MEL**2/0.77d0**2)+2.99d0
      Ichirho=0d0
      APiee=(ALEM0/Pi)**2*MEL*CGPI/2d0
     .                *ARPill(dsqrt(1d0-4d0*MEL**2/METAP**2))
      AIPiee=(ALEM0/Pi)**2*MEL*CGPI/2d0
     .                *AIPill(dsqrt(1d0-4d0*MEL**2/METAP**2))
      Rchirho=3d0/2d0*dlog(MEL**2/0.77d0**2)+6.46d0
      Ichirho=0d0
      AEee=(ALEM0/Pi)**2*MEL*CGETA/2d0
     .                *ARPill(dsqrt(1d0-4d0*MEL**2/METAP**2))
      AIEee=(ALEM0/Pi)**2*MEL*CGETA/2d0
     .                *AIPill(dsqrt(1d0-4d0*MEL**2/METAP**2))
      Rchirho=3d0/2d0*dlog(MEL**2/0.77d0**2)+14.9d0
      Ichirho=2.52d0
      AEPee=(ALEM0/Pi)**2*MEL*CGETAP/2d0
     .                *ARPill(dsqrt(1d0-4d0*MEL**2/METAP**2))
      AIEPee=(ALEM0/Pi)**2*MEL*CGETAP/2d0
     .                *AIPill(dsqrt(1d0-4d0*MEL**2/METAP**2))
      aux=(OMIX(INDETAP,1)*APiee+OMIX(INDETAP,2)*AEee
     .          +OMIX(INDETAP,3)*AEPee+OMIX(INDETAP,4)*AAee)**2
     . +(OMIX(INDETAP,1)*AIPiee+OMIX(INDETAP,2)*AIEee
     .          +OMIX(INDETAP,3)*AIEPee)**2
      WidthEPee=aux/(8d0*Pi)*METAP
     .             *dsqrt(max(0d0,1d0-4d0*MEL**2/METAP**2))          ! Width Eta' -> e+e-
      aux=0.3d0*dabs(aux-AEPee**2-AIEPee**2)/(8d0*Pi)*METAP          ! Mixing uncertainty
     .     *dsqrt(max(0d0,1d0-4d0*MEL**2/METAP**2))
     .   +0.1d0*(AEPee**2+AIEPee**2)*METAP/(8d0*Pi)                  ! SM uncertainty             
     .          *dsqrt(max(0d0,1d0-4d0*MEL**2/METAP**2))
      PROBpill=Max(PROBpill,(WidthPiee-aux)/(1.1d-12)-1d0)

C  * A -> mumu
      AAmumu=MMUON*CD(4)*DSQRT(dsqrt(2d0)*GF)
      APimumu=0d0
      AIPimumu=0d0
      AEmumu=0d0
      AIEmumu=0d0
      AEPmumu=0d0
      AIEPmumu=0d0
      IF(MA.ge.2d0*MMUON)then
       Rchirho=3d0/2d0*dlog(MMUON**2/0.77d0**2)+3.82d0
       Ichirho=0d0
       AEmumu=(ALEM0/Pi)**2*MMUON*CGETA/2d0
     .                *ARPill(dsqrt(1d0-4d0*MMUON**2/MA**2))
! dsqrt(5.8d-6*1.31d-6
!     .            *8d0*Pi/(META*dsqrt(1d0-4d0*MMUON**2/META**2)))
       AIEmumu=(ALEM0/Pi)**2*MMUON*CGETA/2d0
     .                *AIPill(dsqrt(1d0-4d0*MMUON**2/MA**2))
       Rchirho=3d0/2d0*dlog(MMUON**2/0.77d0**2)+6.31d0
       Ichirho=0.75d0
       AEPmumu=(ALEM0/Pi)**2*MMUON*CGETAP/2d0
     .                *ARPill(dsqrt(1d0-4d0*MMUON**2/MA**2))
      AIEPmumu=(ALEM0/Pi)**2*MMUON*CGETAP/2d0
     .                *AIPill(dsqrt(1d0-4d0*MMUON**2/MA**2))
      ENDIF

      aux=(OMIX(INDA,1)*APimumu+OMIX(INDA,2)*AEmumu+OMIX(INDA,3)*AEPmumu
     .          +OMIX(INDA,4)*AAmumu)**2
     . +(OMIX(INDA,1)*AIPimumu+OMIX(INDA,2)*AIEmumu
     .                          +OMIX(INDA,3)*AIEPmumu)**2
      GamAmumu=aux/(8d0*Pi)*MA
     .               *dsqrt(max(0d0,1d0-4d0*MMUON**2/MA**2))

       Rchirho=3d0/2d0*dlog(MMUON**2/0.77d0**2)+3.82d0
       Ichirho=0d0
       aux=(ALEM0/Pi)**2*MMUON*CGETA/2d0*dsqrt(
     . ARPill(dsqrt(1d0-4d0*MMUON**2/META**2))**2+
     . AIPill(dsqrt(1d0-4d0*MMUON**2/META**2))**2)

c Test Eta -> mu+mu-:
      Rchirho=3d0/2d0*dlog(MMUON**2/0.77d0**2)+3.82d0
      Ichirho=0d0
      AEmumu=(ALEM0/Pi)**2*MMUON*CGETA/2d0
     .                *ARPill(dsqrt(1d0-4d0*MMUON**2/META**2))
      AIEmumu=(ALEM0/Pi)**2*MMUON*CGETA/2d0
     .                *AIPill(dsqrt(1d0-4d0*MMUON**2/META**2))
      Rchirho=3d0/2d0*dlog(MMUON**2/0.77d0**2)+6.31d0
      Ichirho=0.75d0
      AEPmumu=(ALEM0/Pi)**2*MMUON*CGETAP/2d0
     .                *ARPill(dsqrt(1d0-4d0*MMUON**2/META**2))
      AIEPmumu=(ALEM0/Pi)**2*MMUON*CGETAP/2d0
     .                *AIPill(dsqrt(1d0-4d0*MMUON**2/META**2))
      aux=(OMIX(INDETA,1)*APimumu+OMIX(INDETA,2)*AEmumu
     .          +OMIX(INDETA,3)*AEPmumu+OMIX(INDETA,4)*AAmumu)**2
     . +(OMIX(INDETA,1)*AIPimumu+OMIX(INDETA,2)*AIEmumu
     .          +OMIX(INDETA,3)*AIEPmumu)**2
      WidthEmumu=aux/(8d0*Pi)*META
     .             *dsqrt(max(0d0,1d0-4d0*MMUON**2/META**2))        ! Width Eta -> mu+mu-
      aux=0.3d0*dabs(aux-AEmumu**2-AIEmumu**2)/(8d0*Pi)*META        ! Mixing uncertainty
     .     *dsqrt(max(0d0,1d0-4d0*MMUON**2/META**2))
     .   +0.1d0*(AEmumu**2+AIEmumu**2)*META/(8d0*Pi)                ! SM uncertainty             
     .          *dsqrt(max(0d0,1d0-4d0*MMUON**2/META**2))
     .   +2d0*0.15d0*5.8d-6*1.31d-6                                 ! Exp uncertainty
      aux=Max(aux,1.1d0*dabs((AEmumu**2+AIEmumu**2)*META/(8d0*Pi)   ! SM-Exp. Error
     .          *dsqrt(max(0d0,1d0-4d0*MMUON**2/META**2))
     .                     -5.8d-6*1.31d-6))
      PROBpill=Max(0d0,dabs(WidthEmumu-5.8d-6*1.31d-6)/aux-1d0)

C  * A -> tautau
      AAtata=MTAU*CL(4)*DSQRT(dsqrt(2d0)*GF)
      GamAtata=(OMIX(INDA,4)*AAtata)**2/(8d0*Pi)*MA
     .               *dsqrt(max(0d0,1d0-4d0*MTAU**2/MA**2))


C       3-meson decays

C   Effective yukawas
      YUE=MPI**2*CU(4)*DSQRT(2d0*dsqrt(2d0)*GF/3d0)/FPI
      YDE=MPI**2*CD(4)*DSQRT(2d0*dsqrt(2d0)*GF/3d0)/FPI
      YSE=(MKC**2+MK0**2-MPI**2)*CD(4)*DSQRT(2d0*dsqrt(2d0)*GF/3d0)/FPI

C  * A -> 3Pi0
      IF(MA.lt.1.1d0)then
       AA3Pi=-DSQRT(dsqrt(2d0)*GF)*MPI**2/(2d0*FPI)*(CU(4)-CD(4))
      ELSE
       IF(MA.lt.1.3d0)then
        AA3Pi=-DSQRT(dsqrt(2d0)*GF)*MPI**2/(2d0*FPI)*(CU(4)-CD(4))
     .        *(1d0-(MA-1.1d0)/.2d0)
     . +dsqrt(6d0*(YUE**2+YDE**2)/16d0)*SIGN(1d0,CD(4)-CU(4))
     .        *(MA-1.1d0)/.2d0
       ELSE
        AA3Pi=dsqrt(6d0*(YUE**2+YDE**2)/16d0)
     .   *SIGN(1d0,CD(4)-CU(4))
       ENDIF
      ENDIF
      APi3Pi=(MPI/FPI)**2
      AE3Pi=0.722665d0
      AEP3Pi=0.278538d0

      aux=OMIX(INDA,1)*APi3Pi+OMIX(INDA,2)*AE3Pi+OMIX(INDA,3)*AEP3Pi
     .          +OMIX(INDA,4)*AA3Pi
      GamA3Pi=aux**2/(6d0*2d0**8*Pi**3*MA)*KinA3Pi(MA)

C  * A -> Pi0Pi+Pi-
      IF(MA.lt.1.1d0)then
       AAPi3PiC=-DSQRT(dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)*(CU(4)-CD(4))
      ELSE
       IF(MA.lt.1.3d0)then
        AAPi3PiC=-DSQRT(dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)*(CU(4)-CD(4))
     .        *(1d0-(MA-1.1d0)/.2d0)
     . +dsqrt((YUE**2+YDE**2)/24d0)*SIGN(1d0,CD(4)-CU(4))
     .       *(MA-1.1d0)/.2d0
       ELSE
        AAPi3PiC=dsqrt((YUE**2+YDE**2)/24d0)*SIGN(1d0,CD(4)-CU(4))
       ENDIF
      ENDIF
      APiPi3PiC=(MPI/FPI)**2/3d0
      AEPi3PiC=0.262726d0
      AEPPi3PiC=0.151766d0

      aux=OMIX(INDA,1)*APiPi3PiC+OMIX(INDA,2)*AEPi3PiC+OMIX(INDA,3)
     .          *AEPPi3PiC+OMIX(INDA,4)*AAPi3PiC
      GamAPi3PiC=aux**2/(2d0**8*Pi**3*MA)*KinAPi3PiC(MA)

C  * A -> Eta2Pi0
      IF(MA.lt.1.1d0)then
       AAEPi3=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)*(CU(4)+CD(4))
      ELSE
       IF(MA.lt.1.3d0)then
        AAEPi3=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)*(CU(4)+CD(4))
     .        *(1d0-(MA-1.1d0)/.2d0)
     . -dsqrt(2d0*(YUE**2+YDE**2)/16d0)*SIGN(1d0,CD(4)+CU(4))
     .               *(MA-1.1d0)/.2d0
       ELSE
        AAEPi3=-dsqrt(2d0*(YUE**2+YDE**2)/16d0)
     .               *SIGN(1d0,CD(4)+CU(4))
       ENDIF
      ENDIF
      AAEPi3=AAEPi3*(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
      APIEPi3=(MKC**2-MK0**2-MPIC**2+MPI**2)/(dsqrt(3d0)*FPI**2)
     .        *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
      AEEPi3=(MPI/FPI)**2/3d0*(DDCOS(THETAE)
     .       -dsqrt(2d0)*DDSIN(THETAE))**2
      AEPEPi3=6.33755d0

      aux=OMIX(INDA,1)*APIEPi3+OMIX(INDA,2)*AEEPi3+OMIX(INDA,3)
     .          *AEPEPi3+OMIX(INDA,4)*AAEPi3
      GamAEPi3=aux**2/(2d0*2d0**8*Pi**3*MA)*KinAEPi3(MA)

C  * A -> EtaPi+Pi-
      IF(MA.lt.1.1d0)then
       AAEPiC=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)*(CU(4)+CD(4))
      ELSE
       IF(MA.lt.1.3d0)then
        AAEPiC=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)*(CU(4)+CD(4))
     .        *(1d0-(MA-1.1d0)/.2d0)
     . -dsqrt((YUE**2+YDE**2)/8d0)*SIGN(1d0,CD(4)+CU(4))
     .               *(MA-1.1d0)/.2d0
       ELSE
        AAEPiC=-dsqrt((YUE**2+YDE**2)/8d0)
     .               *SIGN(1d0,CD(4)+CU(4))
       ENDIF
      ENDIF
      AAEPiC=AAEPiC*(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
      APIEPiC=(MKC**2-MK0**2-MPIC**2+MPI**2)/(3d0*dsqrt(3d0)*FPI**2)
     .        *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
      AEEPiC=(MPI/FPI)**2/3d0*(DDCOS(THETAE)
     .       -dsqrt(2d0)*DDSIN(THETAE))**2
      AEPEPiC=6.60772d0

      aux=OMIX(INDA,1)*APIEPiC+OMIX(INDA,2)*AEEPiC+OMIX(INDA,3)
     .          *AEPEPiC+OMIX(INDA,4)*AAEPiC
      GamAEPiC=aux**2/(2d0**8*Pi**3*MA)*KinAEPiC(MA)

C  * A -> Eta'2Pi0
      IF(MA.lt.1.1d0)then
       AAEPPi3=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)*(CU(4)+CD(4))
      ELSE
       IF(MA.lt.1.3d0)then
        AAEPPi3=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .        *(CU(4)+CD(4))*(1d0-(MA-1.1d0)/.2d0)
     . -dsqrt(2d0*(YUE**2+YDE**2)/16d0)*SIGN(1d0,CD(4)+CU(4))
     .               *(MA-1.1d0)/.2d0
       ELSE
        AAEPPi3=-dsqrt(2d0*(YUE**2+YDE**2)/16d0)
     .               *SIGN(1d0,CD(4)+CU(4))
       ENDIF
      ENDIF
      AAEPPi3=AAEPPi3*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
      APIEPPi3=(MKC**2-MK0**2-MPIC**2+MPI**2)/(dsqrt(3d0)*FPI**2)
     .        *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
      AEEPPi3=(MPI/FPI)**2/3d0*(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
     .                        *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
      AEPEPPi3=5.65362d0

      aux=OMIX(INDA,1)*APIEPPi3+OMIX(INDA,2)*AEEPPi3+OMIX(INDA,3)
     .          *AEPEPPi3+OMIX(INDA,4)*AAEPPi3
      GamAEPPi3=aux**2/(2d0*2d0**8*Pi**3*MA)*KinAEPPi3(MA)

C  * A -> Eta'Pi+Pi-
      IF(MA.lt.1.1d0)then
       AAEPPiC=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)*(CU(4)+CD(4))
      ELSE
       IF(MA.lt.1.3d0)then
        AAEPPiC=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .        *(CU(4)+CD(4))*(1d0-(MA-1.1d0)/.2d0)
     . -dsqrt((YUE**2+YDE**2)/8d0)*SIGN(1d0,CD(4)+CU(4))
     .               *(MA-1.1d0)/.2d0
       ELSE
        AAEPPiC=-dsqrt((YUE**2+YDE**2)/8d0)*SIGN(1d0,CD(4)+CU(4))
       ENDIF
      ENDIF
      AAEPPiC=AAEPPiC*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
      APIEPPiC=(MKC**2-MK0**2-MPIC**2+MPI**2)/(3d0*dsqrt(3d0)*FPI**2)
     .        *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
      AEEPPiC=(MPI/FPI)**2/3d0*(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
     .        *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
      AEPEPPiC=5.89464d0

      aux=OMIX(INDA,1)*APIEPPiC+OMIX(INDA,2)*AEEPPiC+OMIX(INDA,3)
     .          *AEPEPPiC+OMIX(INDA,4)*AAEPPiC
      GamAEPPiC=aux**2/(2d0**8*Pi**3*MA)*KinAEPPiC(MA)

C  * A -> Pi02Eta
      IF(MA.lt.1.1d0)then
       AAPiEE=-DSQRT(dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)*(CU(4)-CD(4))
      ELSE
       IF(MA.lt.1.3d0)then
        AAPiEE=-DSQRT(dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .        *(CU(4)-CD(4))*(1d0-(MA-1.1d0)/.2d0)
     . +dsqrt(2d0*(YUE**2+YDE**2)/48d0)*SIGN(1d0,CD(4)-CU(4))
     .               *(MA-1.1d0)/.2d0
       ELSE
        AAPiEE=dsqrt(2d0*(YUE**2+YDE**2)/48d0)*SIGN(1d0,CD(4)-CU(4))
       ENDIF
      ENDIF
      AAPiEE=AAPiEE*(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))**2
      APIPiEE=(MPI/FPI)**2/3d0
     .        *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))**2
      AEPiEE=(MKC**2-MK0**2-MPIC**2+MPI**2)/(3d0*dsqrt(3d0)*FPI**2)
     .        *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))**3
      AEPPiEE=(MKC**2-MK0**2-MPIC**2+MPI**2)/(3d0*dsqrt(3d0)*FPI**2)
     .        *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))**2
     .        *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))

      aux=OMIX(INDA,1)*APIPiEE+OMIX(INDA,2)*AEPiEE+OMIX(INDA,3)
     .          *AEPPiEE+OMIX(INDA,4)*AAPiEE
      GamAPiEE=aux**2/(2d0*2d0**8*Pi**3*MA)*KinAPiEE(MA)

C  * A -> Pi0EtaEta'
      AAPiEEP=dsqrt((YUE**2+YDE**2)/24d0)*SIGN(1d0,CD(4)-CU(4))
      AAPiEEP=AAPiEEP*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
     .               *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
      APIPiEEP=(MPI/FPI)**2/3d0
     .               *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
     .               *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
      AEPiEEP=(MKC**2-MK0**2-MPIC**2+MPI**2)/(3d0*dsqrt(3d0)*FPI**2)
     .        *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))**2
     .        *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
      AEPPiEEP=(MKC**2-MK0**2-MPIC**2+MPI**2)/(3d0*dsqrt(3d0)*FPI**2)
     .        *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
     .        *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))**2

      aux=OMIX(INDA,1)*APIPiEEP+OMIX(INDA,2)*AEPiEEP+OMIX(INDA,3)
     .          *AEPPiEEP+OMIX(INDA,4)*AAPiEEP
      GamAPiEEP=aux**2/(2d0**8*Pi**3*MA)*KinAPiEEP(MA)

C  * A -> Pi02Eta'
      AAPiEPEP=dsqrt(2d0*(YUE**2+YDE**2)/48d0)*SIGN(1d0,CD(4)-CU(4))
      AAPiEPEP=AAPiEPEP*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))**2
      APIPiEPEP=(MPI/FPI)**2/3d0
     .        *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))**2
      AEPiEPEP=(MKC**2-MK0**2-MPIC**2+MPI**2)/(3d0*dsqrt(3d0)*FPI**2)
     .        *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))**2
     .        *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
      AEPPiEPEP=(MKC**2-MK0**2-MPIC**2+MPI**2)/(3d0*dsqrt(3d0)*FPI**2)
     .        *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))**3

      aux=OMIX(INDA,1)*APIPiEPEP+OMIX(INDA,2)*AEPiEPEP+OMIX(INDA,3)
     .          *AEPPiEPEP+OMIX(INDA,4)*AAPiEPEP
      GamAPiEPEP=aux**2/(2d0*2d0**8*Pi**3*MA)*KinAPiEPEP(MA)

C  * A -> 3Eta
      AA3E=-dsqrt(6d0*(YUE**2+YDE**2+64d0*YSE**2)/432d0)
     . *SIGN(1d0,CU(4)*(MKC**2+MPI**2-MK0**2)+CD(4)
     .     *(MK0**2+MPI**2-MKC**2)-8d0*CD(4)*(MKC**2+MK0**2-MPI**2))
      API3E=(MKC**2-MK0**2-MPIC**2+MPI**2)/(3d0*dsqrt(3d0)*FPI**2)
!     .        *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))**3
      AE3E=(2d0*MPI**2+16d0*(MKC**2+MK0**2-MPI**2))/(9d0*FPI**2)
      AEP3E=dsqrt(2d0)*(2d0*MPI**2-8d0*(MKC**2+MK0**2-MPI**2))
     .             /(9d0*FPI**2)

      aux=OMIX(INDA,1)*API3E+OMIX(INDA,2)*AE3E+OMIX(INDA,3)
     .          *AEP3E+OMIX(INDA,4)*AA3E
      GamA3E=aux**2/(6d0*2d0**8*Pi**3*MA)*KinA3E(MA)

C  * A -> 2EtaEta'
      AAE2EP=-dsqrt(2d0*(YUE**2+YDE**2+16d0*YSE**2)/72d0)
     . *SIGN(1d0,CU(4)*(MKC**2+MPI**2-MK0**2)+CD(4)
     .     *(MK0**2+MPI**2-MKC**2)+4d0*CD(4)*(MKC**2+MK0**2-MPI**2))
      APIE2EP=(MKC**2-MK0**2-MPIC**2+MPI**2)/(3d0*dsqrt(3d0)*FPI**2)
     .        *dsqrt(2d0)
!     .        *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))**2
!     .        *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
      AEE2EP=dsqrt(2d0)*(2d0*MPI**2-8d0*(MKC**2+MK0**2-MPI**2))
     .             /(9d0*FPI**2)
      AEPE2EP=2d0*(2d0*MPI**2+4d0*(MKC**2+MK0**2-MPI**2))
     .             /(9d0*FPI**2)

      aux=OMIX(INDA,1)*APIE2EP+OMIX(INDA,2)*AEE2EP+OMIX(INDA,3)
     .          *AEPE2EP+OMIX(INDA,4)*AAE2EP
      GamAE2EP=aux**2/(2d0*2d0**8*Pi**3*MA)*KinAE2EP(MA)

C  * A -> Eta2Eta'
      AAEEP2=-dsqrt(2d0*(YUE**2+YDE**2+4d0*YSE**2)/36d0)
     . *SIGN(1d0,CU(4)*(MKC**2+MPI**2-MK0**2)+CD(4)
     .     *(MK0**2+MPI**2-MKC**2)-2d0*CD(4)*(MKC**2+MK0**2-MPI**2))
      APIEEP2=(MKC**2-MK0**2-MPIC**2+MPI**2)/(3d0*dsqrt(3d0)*FPI**2)
     .        *dsqrt(2d0)**2
!     .        *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
!     .        *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))**2
      AEEEP2=2d0*(2d0*MPI**2+4d0*(MKC**2+MK0**2-MPI**2))
     .             /(9d0*FPI**2)
      AEPEEP2=2d0*dsqrt(2d0)*(2d0*MPI**2-2d0*(MKC**2+MK0**2-MPI**2))
     .             /(9d0*FPI**2)

      aux=OMIX(INDA,1)*APIEEP2+OMIX(INDA,2)*AEEEP2+OMIX(INDA,3)
     .          *AEPEEP2+OMIX(INDA,4)*AAEEP2
      GamAEEP2=aux**2/(2d0*2d0**8*Pi**3*MA)*KinAEEP2(MA)

C  * A -> 3Eta'
      AA3EP=-dsqrt(6d0*(YUE**2+YDE**2+YSE**2)/54d0)
     . *SIGN(1d0,CU(4)*(MKC**2+MPI**2-MK0**2)+CD(4)
     .     *(MK0**2+MPI**2-MKC**2)+CD(4)*(MKC**2+MK0**2-MPI**2))
      API3EP=(MKC**2-MK0**2-MPIC**2+MPI**2)/(3d0*dsqrt(3d0)*FPI**2)
     .        *dsqrt(2d0)**3
!     .        *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))**3
      AE3EP=2d0*dsqrt(2d0)*(2d0*MPI**2-2d0*(MKC**2+MK0**2-MPI**2))
     .             /(9d0*FPI**2)
      AEP3EP=4d0*(2d0*MPI**2+(MKC**2+MK0**2-MPI**2))/(9d0*FPI**2)

      aux=OMIX(INDA,1)*API3EP+OMIX(INDA,2)*AE3EP+OMIX(INDA,3)
     .          *AEP3EP+OMIX(INDA,4)*AA3EP
      GamA3EP=aux**2/(6d0*2d0**8*Pi**3*MA)*KinA3EP(MA)

C  * A -> Pi0K+K-
      IF(MA.lt.1.1d0)then
       AAPiKC=-DSQRT(dsqrt(2d0)*GF)/(6d0*FPI)*(MKC**2*(2d0*CU(4)+CD(4))
     .                             +(MPI**2-MK0**2)*(2d0*CU(4)-CD(4)))
      ELSE
       IF(MA.lt.1.3d0)then
        AAPiKC=-DSQRT(dsqrt(2d0)*GF)/(6d0*FPI)
     .    *(MKC**2*(2d0*CU(4)+CD(4))+(MPI**2-MK0**2)*(2d0*CU(4)-CD(4)))
     .               *(1d0-(MA-1.1d0)/.2d0)
     . -dsqrt((4d0*YUE**2+YSE**2)/24d0)*SIGN(1d0,MKC**2
     .   *(2d0*CU(4)+CD(4))+(MPI**2-MK0**2)*(2d0*CU(4)-CD(4)))
     .               *(MA-1.1d0)/.2d0
       ELSE
        AAPiKC=-dsqrt((4d0*YUE**2+YSE**2)/24d0)*SIGN(1d0,MKC**2
     .    *(2d0*CU(4)+CD(4))+(MPI**2-MK0**2)*(2d0*CU(4)-CD(4)))
       ENDIF
      ENDIF
      APIPiKC=(2d0*MKC**2-MK0**2+MPI**2)/(6d0*FPI**2)
      AEPiKC=(-3d0*dsqrt(2d0)*DDSIN(THETAE)*MKC**2
     .        -(MK0**2-MPI**2)*(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE)))
     . /(6d0*dsqrt(3d0)*FPI**2)
      AEPPiKC=(3d0*dsqrt(2d0)*DDCOS(THETAE)*MKC**2
     .        -(MK0**2-MPI**2)*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE)))
     . /(6d0*dsqrt(3d0)*FPI**2)

      aux=OMIX(INDA,1)*APIPiKC+OMIX(INDA,2)*AEPiKC+OMIX(INDA,3)
     .          *AEPPiKC+OMIX(INDA,4)*AAPiKC
      GamAPiKC=aux**2/(2d0**8*Pi**3*MA)*KinAPiKC(MA)

C  * A -> Pi0K0K0b
      IF(MA.lt.1.1d0)then
       AAPiK0=-DSQRT(dsqrt(2d0)*GF)/(6d0*FPI)*CD(4)
     .                              *(MKC**2-MPI**2-3d0*MK0**2)
      ELSE
       IF(MA.lt.1.3d0)then
        AAPiK0=-DSQRT(dsqrt(2d0)*GF)/(6d0*FPI)*CD(4)
     .                              *(MKC**2-MPI**2-3d0*MK0**2)
     .               *(1d0-(MA-1.1d0)/.2d0)
     . -dsqrt((4d0*YDE**2+YSE**2)/24d0)
     .    *SIGN(1d0,CD(4)*(MKC**2-MPI**2-3d0*MK0**2))
     .               *(MA-1.1d0)/.2d0
       ELSE
        AAPiK0=-dsqrt((4d0*YDE**2+YSE**2)/24d0)
     .    *SIGN(1d0,CD(4)*(MKC**2-MPI**2-3d0*MK0**2))
       ENDIF
      ENDIF
      APIPiK0=(2d0*MK0**2-MKC**2+MPI**2)/(6d0*FPI**2)
      AEPiK0=(3d0*dsqrt(2d0)*DDSIN(THETAE)*MK0**2
     .        +(MKC**2-MPI**2)*(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE)))
     . /(6d0*dsqrt(3d0)*FPI**2)
      AEPPiK0=(-3d0*dsqrt(2d0)*DDCOS(THETAE)*MK0**2
     .        -(MKC**2-MPI**2)*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE)))
     . /(6d0*dsqrt(3d0)*FPI**2)

      aux=OMIX(INDA,1)*APIPiK0+OMIX(INDA,2)*AEPiK0+OMIX(INDA,3)
     .          *AEPPiK0+OMIX(INDA,4)*AAPiK0
      GamAPiK0=aux**2/(2d0**8*Pi**3*MA)*KinAPiK0(MA)

C  * A -> PiK+K0b + PiK-K0
      IF(MA.lt.1.1d0)then
       AAPiKCK0=-DSQRT(2d0*dsqrt(2d0)*GF)/(6d0*FPI)
     .     *(CU(4)*(MKC**2+MPI**2-MK0**2)+2d0*CD(4)*MK0**2)
      ELSE
       IF(MA.lt.1.3d0)then
        AAPiKCK0=-DSQRT(2d0*dsqrt(2d0)*GF)/(6d0*FPI)
     .     *(CU(4)*(MKC**2+MPI**2-MK0**2)+2d0*CD(4)*MK0**2)
     .               *(1d0-(MA-1.1d0)/.2d0)
     . -dsqrt((YUE**2+YDE**2+YSE**2)/12d0)
     .    *SIGN(1d0,CU(4)*(MKC**2+MPI**2-MK0**2)+2d0*CD(4)*MK0**2)
     .               *(MA-1.1d0)/.2d0
       ELSE
        AAPiKCK0=-dsqrt((YUE**2+YDE**2+YSE**2)/12d0)
     .    *SIGN(1d0,CU(4)*(MKC**2+MPI**2-MK0**2)+2d0*CD(4)*MK0**2)
       ENDIF
      ENDIF
      APIPiKCK0=(MKC**2-MK0**2+MPI**2-MPIC**2)/(6d0*dsqrt(2d0)*FPI**2)
      AEPiKCK0=((2d0*MPI**2-MKC**2-MK0**2)*DDCOS(THETAE)
     .    -2d0*dsqrt(2d0)*(MPI**2+MKC**2+MK0**2)*DDSIN(THETAE))
     . /(6d0*dsqrt(6d0)*FPI**2)
      AEPPiKCK0=((2d0*MPI**2-MKC**2-MK0**2)*DDSIN(THETAE)
     .    +2d0*dsqrt(2d0)*(MPI**2+MKC**2+MK0**2)*DDCOS(THETAE))
     . /(6d0*dsqrt(6d0)*FPI**2)

      aux=OMIX(INDA,1)*APIPiKCK0+OMIX(INDA,2)*AEPiKCK0+OMIX(INDA,3)
     .          *AEPPiKCK0+OMIX(INDA,4)*AAPiKCK0
      GamAPiKCK0=2d0*aux**2/(2d0**8*Pi**3*MA)*KinAPiKCK0(MA)

C  * A -> EtaK+K-
      AAEKC=dsqrt((YSE**2)/8d0)
      APIEKC=-(MK0**2-MPI**2)/(6d0*dsqrt(3d0)*FPI**2)
c     (-3d0*dsqrt(2d0)*DDSIN(THETAE)*MKC**2
c     .        -(MK0**2-MPI**2)*(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE)))
c     . /(6d0*dsqrt(3d0)*FPI**2)
      AEEKC=(MKC**2+(MK0**2-MPI**2)/2d0)/3d0/FPI**2
c     ((MKC**2+(MK0**2-MPI**2)/2d0)*DDCOS(THETAE)**2
c     .  +dsqrt(2d0)*(MKC**2+MK0**2-MPI**2)*DDCOS(THETAE)*DDSIN(THETAE)
c     .  +MKC**2*DDSIN(THETAE)**2)/3d0/FPI**2
      AEPEKC=-dsqrt(2d0)*(MKC**2+MK0**2-MPI**2)/6d0/FPI**2
c    ((2d0*MKC**2+MK0**2-MPI**2)*DDCOS(THETAE)*DDSIN(THETAE)
c     .  -dsqrt(2d0)*(MKC**2+MK0**2-MPI**2)*DDCOS(2d0*THETAE)
c     .  -2d0*MKC**2*DDSIN(THETAE)*DDCOS(THETAE))/6d0/FPI**2

      aux=OMIX(INDA,1)*APIEKC+OMIX(INDA,2)*AEEKC+OMIX(INDA,3)
     .          *AEPEKC+OMIX(INDA,4)*AAEKC
      GamAEKC=aux**2/(2d0**8*Pi**3*MA)*KinAEKC(MA)

C  * A -> EtaK0K0b
      AAEK0=dsqrt((YSE**2)/8d0)*SIGN(1d0,CD(4)*(MKC**2+MK0**2-MPI**2))
      APIEK0=(MKC**2-MPI**2)/(6d0*dsqrt(3d0)*FPI**2)
c     (3d0*dsqrt(2d0)*DDSIN(THETAE)*MK0**2
c     .        +(MKC**2-MPI**2)*(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE)))
c     . /(6d0*dsqrt(3d0)*FPI**2)
      AEEK0=(MK0**2+(MKC**2-MPI**2)/2d0)/3d0/FPI**2
c     ((MK0**2+(MKC**2-MPI**2)/2d0)*DDCOS(THETAE)**2
c     .  +dsqrt(2d0)*(MK0**2+MKC**2-MPI**2)*DDCOS(THETAE)*DDSIN(THETAE)
c     .  +MK0**2*DDSIN(THETAE)**2)/3d0/FPI**2
      AEPEK0=-dsqrt(2d0)*(MK0**2+MKC**2-MPI**2)/6d0/FPI**2
c     ((2d0*MK0**2+MKC**2-MPI**2)*DDCOS(THETAE)*DDSIN(THETAE)
c     .  -dsqrt(2d0)*(MK0**2+MKC**2-MPI**2)*DDCOS(2d0*THETAE)
c     .  -2d0*MK0**2*DDSIN(THETAE)*DDCOS(THETAE))/6d0/FPI**2

      aux=OMIX(INDA,1)*APIEK0+OMIX(INDA,2)*AEEK0+OMIX(INDA,3)
     .          *AEPEK0+OMIX(INDA,4)*AAEK0
      GamAEK0=aux**2/(2d0**8*Pi**3*MA)*KinAEK0(MA)

C  * A -> Eta'K+K-
      AAEPKC=-dsqrt((YUE**2+YSE**2)/4d0)*SIGN(1d0,
     . CU(4)*(MPI**2+MKC**2-MK0**2)+CD(4)*(MK0**2+MPI**2-MKC**2))
      APIEPKC=dsqrt(2d0/3d0)*(3d0*MKC**2-(MK0**2-MPI**2))/(6d0*FPI**2)
c     (3d0*dsqrt(2d0)*DDCOS(THETAE)*MKC**2
c     .        -(MK0**2-MPI**2)*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE)))
c     . /(6d0*dsqrt(3d0)*FPI**2)
      AEEPKC=-dsqrt(2d0)*(MKC**2+MK0**2-MPI**2)/6d0/FPI**2
c     ((2d0*MKC**2+MK0**2-MPI**2)*DDCOS(THETAE)*DDSIN(THETAE)
c     .  -dsqrt(2d0)*(MKC**2+MK0**2-MPI**2)*DDCOS(2d0*THETAE)
c     .  -2d0*MKC**2*DDSIN(THETAE)*DDCOS(THETAE))/6d0/FPI**2
      AEPEPKC=MKC**2/3d0/FPI**2
c     ((MKC**2+(MK0**2-MPI**2)/2d0)*DDSIN(THETAE)**2
c     .  -dsqrt(2d0)*(MKC**2+MK0**2-MPI**2)*DDCOS(THETAE)*DDSIN(THETAE)
c     .  +MKC**2*DDCOS(THETAE)**2)/3d0/FPI**2

      aux=OMIX(INDA,1)*APIEPKC+OMIX(INDA,2)*AEEPKC+OMIX(INDA,3)
     .          *AEPEPKC+OMIX(INDA,4)*AAEPKC
      GamAEPKC=aux**2/(2d0**8*Pi**3*MA)*KinAEPKC(MA)

C  * A -> Eta'K0K0b
      AAEPK0=dsqrt((YDE**2+YSE**2)/4d0)
      APIEPK0=dsqrt(2d0/3d0)*(-3d0*MK0**2-MKC**2+MPI**2)/(6d0*FPI**2)
c     (-3d0*dsqrt(2d0)*DDCOS(THETAE)*MK0**2
c     .        -(MKC**2-MPI**2)*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE)))
c     . /(6d0*dsqrt(3d0)*FPI**2)
      AEEPK0=-dsqrt(2d0)*(MK0**2+MKC**2-MPI**2)/6d0/FPI**2
c     ((2d0*MK0**2+MKC**2-MPI**2)*DDCOS(THETAE)*DDSIN(THETAE)
c     .  -dsqrt(2d0)*(MK0**2+MKC**2-MPI**2)*DDCOS(2d0*THETAE)
c     .  -2d0*MK0**2*DDSIN(THETAE)*DDCOS(THETAE))/6d0/FPI**2
      AEPEPK0=MK0**2/3d0/FPI**2
c     ((MK0**2+(MKC**2-MPI**2)/2d0)*DDSIN(THETAE)**2
c     .  -dsqrt(2d0)*(MK0**2+MKC**2-MPI**2)*DDCOS(THETAE)*DDSIN(THETAE)
c     .  +MK0**2*DDCOS(THETAE)**2)/3d0/FPI**2

      aux=OMIX(INDA,1)*APIEPK0+OMIX(INDA,2)*AEEPK0+OMIX(INDA,3)
     .          *AEPEPK0+OMIX(INDA,4)*AAEPK0
      GamAEPK0=aux**2/(2d0**8*Pi**3*MA)*KinAEPK0(MA)

C       Radiative hadronic decays

C  * A -> gamma (Rho,omega->)Pi+Pi-
      AARhogam=0d0
      APIRhogam=1d0
      AERhogam=0.745672d0 !1.22374d0
      AEPRhogam=0.9274d0  !1.39445d0

      aux=OMIX(INDA,1)*APIRhogam+OMIX(INDA,2)*AERhogam+OMIX(INDA,3)
     .          *AEPRhogam+OMIX(INDA,4)*AARhogam
      GamARhogam=4d0*Pi*ALEM0/(4d0*Pi**2*FPI**3)**2
     .            *aux**2/(2d0**8*Pi**3*MA)*KinARhogam(MA)

C       Total hadronic decays

      GamAhadr=GamA3Pi+GamAPi3PiC+GamAEPi3+GamAEPiC+GamAEPPi3+GamAEPPiC
     .   +GamAPiEE+GamAPiEEP+GamAPiEPEP+GamA3E+GamAE2EP+GamAEEP2
     .   +GamA3EP+GamAPiKC+GamAPiK0+GamAPiKCK0+GamAEKC+GamAEK0
     .   +GamAEPKC+GamAEPK0+GamARhogam

      ENDIF

C       Light-quark and gluon decays

      AAuu=2d-3*CU(4)*DSQRT(dsqrt(2d0)*GF)
      AAdd=4d-3*CD(4)*DSQRT(dsqrt(2d0)*GF)
      AAss=MS*CD(4)*DSQRT(dsqrt(2d0)*GF)

C   Initializing the strong coupling constant and running masses
      IF(MA.LT.3d0)THEN
       GamAuu=0d0
       GamAdd=0d0
       GamAss=0d0
       GamAjj=0d0
      ELSE
       HIGTOP=(MA/MT)**2
       MT0=3d8
* New May 2019:
       ASH=ALPHAS(MA,3)
       MC0=1d8
       MB0=2d8
       AS3=ALPHAS(MA,3)
       MC0=MCC
       AS4=ALPHAS(MA,3)
* End New
       MB0=MBP
       MT0=MT
       RMS=RUNM(MA,3)

       RATCOUP=1d0
       GamAuu=3d0*dabs(AAuu)**2/(8d0*Pi)*MA
     .        *(4d0*(2d-3/MA)**2*TQCDA((2d-3/MA)**2)
     .   +(1d0-4d0*(2d-3/MA)**2)*max(0d0,QCDA((2d-3/MA)**2)))
     .               *dsqrt(max(0d0,1d0-4d0*(2d-3/MA)**2))

       RATCOUP=0d0
       IF(CD(4).NE.0d0)RATCOUP=CU(4)/CD(4)
       GamAdd=3d0*dabs(AAdd)**2/(8d0*Pi)*MA
     .         *(4d0*(4d-3/MA)**2*TQCDA((4d-3/MA)**2)
     .   +(1d0-4d0*(4d-3/MA)**2)*max(0d0,QCDA((4d-3/MA)**2)))
     .               *dsqrt(max(0d0,1d0-4d0*(4d-3/MA)**2))
       GamAss=3d0*dabs(AAss)**2/(8d0*Pi)*MA
     .         *(4d0*(MS/MA)**2*TQCDA((MS/MA)**2)
     .   +(1d0-4d0*(MS/MA)**2)*max(0d0,(RMS/MS)**2*QCDA((RMS/MA)**2)))
     .               *dsqrt(max(0d0,1d0-4d0*(MS/MA)**2))
       GamAjj=CDABS(CJA)**2/(8d0*PI**3)*MA**3
     .     *AS3**2*(1d0+AS3/Pi*(97d0/4d0-3d0*7d0/6d0))
      ENDIF

      aux=GamAhadr
      IF(MA.lt.3d0)then
       GamAhadr=aux
      ELSE
       IF(MA.lt.4d0)THEN
        GamAhadr=aux*(4d0-MA)+(MA-3d0)*(0d0*GamAuu+0d0*GamAdd
     .                                         +GamAss+GamAjj)
       ELSE
        GamAhadr=0d0*GamAuu+0d0*GamAdd+GamAss+GamAjj
       ENDIF
      ENDIF

C       Mixing with the Eta_c(1S)

      METAC=2.9834d0
! eta_c(1S) wave-function from J/Psi(1S)->e+e-
      DMC=dsqrt(3d0/4d0*METAC**3*0.05971d0*92.9d-6*2d0*dsqrt(2d0)*GF) 
     .    *3.0969d0*CU(4)/(4d0*Pi*ALEM0*2d0/3d0)

c      IF(MA.gt.2d0.and.MA.lt.2d0*MTAU)THEN
       CMIX=0d0
       aux=Max(0d0,min(1d0,MA-1d0))
       IF(DMC.ne.0d0)THEN
        CMIX=datan(2d0*DMC
     .        /(MA**2-METAC**2+dsqrt((MA**2-METAC**2)**2+4d0*DMC**2)))
       ENDIF
       IF(DDCOS(CMIX)**2.gt.DDSIN(CMIX)**2)then
        GamAhadr=GamAhadr*DDCOS(CMIX)**2+31.8d-3*aux*DDSIN(CMIX)**2
        GamAee=GamAee*DDCOS(CMIX)**2
        GamAmumu=GamAmumu*DDCOS(CMIX)**2
        GamAtata=GamAtata*DDCOS(CMIX)**2
        GamAgaga=GamAgaga*DDCOS(CMIX)**2+31.8d-3*1.59d-4*aux
     .                                                 *DDSIN(CMIX)**2
        GamAinv=GamAinv*DDCOS(CMIX)**2
       ELSE
        GamAhadr=GamAhadr*DDSIN(CMIX)**2+31.8d-3*aux*DDCOS(CMIX)**2
        GamAee=GamAee*DDSIN(CMIX)**2
        GamAmumu=GamAmumu*DDSIN(CMIX)**2
        GamAtata=GamAtata*DDSIN(CMIX)**2
        GamAgaga=GamAgaga*DDSIN(CMIX)**2+31.8d-3*1.59d-4*aux
     .                                                 *DDCOS(CMIX)**2
        GamAinv=GamAinv*DDSIN(CMIX)**2
       ENDIF
c      ENDIF
! We neglect the mixing with the eta_c(2,3S) because they matter after the opening of the tautau channel

c Test eta_c(1S)->2gamma remains within 15% of the experimental value:
       IF(DDCOS(CMIX)**2.gt.DDSIN(CMIX)**2)then
      WidthECgaga=8d0*CDABS(CGA)**2*METAC**3
     .          /(32d0*Pi)*(ALEM0/4d0/Pi)**2*DDSIN(CMIX)**2*aux
     .           +31.8d-3*1.59d-4*DDCOS(CMIX)**2
       ELSE
      WidthECgaga=8d0*CDABS(CGA)**2*METAC**3
     .          /(32d0*Pi)*(ALEM0/4d0/Pi)**2*DDCOS(CMIX)**2*aux
     .           +31.8d-3*1.59d-4*DDSIN(CMIX)**2
       ENDIF
      PROBpigaga=Max(PROBpigaga,
     .           dabs(WidthECgaga/(31.8d-3*1.59d-4)-1d0)/0.15d0-1d0)


C       A -> cc decays

      MD=1.865d0

      IF(MA.LE.2d0*MD+MPI)THEN
       GamAcc= 0d0
      ELSE
       RMC=RUNM(MA,4)
       RATCOUP= 1d0

       DCC=CDABS(CJA)**2/(8d0*PI**3)*MA**3
     .     *(AS4**2*(1d0+AS4/Pi*(97d0/4d0-4d0*7d0/6d0))
     .      -AS3**2*(1d0+AS3/Pi*(97d0/4d0-3d0*7d0/6d0)))

* New July 2019:
       GamAcc=4d0*(MCC/MA)**2*
     .        3d0*GF*(MCC*CU(4))**2*MA/(4d0*dsqrt(2d0)*Pi)
     .        *TQCDA((MCC/MA)**2)*dsqrt(1d0-4d0*MCC**2/MA**2)
     .         +(1d0-4d0*(MCC/MA)**2)*max(0d0,
     .         3d0*GF*(RMC*CU(4))**2*MA/(4d0*dsqrt(2d0)*Pi)
     .         *QCDA((RMC/MA)**2)*dsqrt(1d0-4d0*RMC**2/MA**2)!+DCC
     .         )
* End new
! Modifying the shape of the phase space function, so that it looks "3-body" at low mass
       IF(MA.le.50d0)THEN
        GamAcc=GamAcc*dsqrt(1d0-(2d0*MD+MPI)**2/MA**2)
     .                /dsqrt(1d0-4d0*MCC**2/MA**2)
     . *(KinAPiDD(MA)*((50d0-MA)/(50d0-(2d0*MD+MPI)))**2
     .                 +1d0-((50d0-MA)/(50d0-(2d0*MD+MPI)))**2)

* New July 2019:
        GamAhadr=GamAhadr+DCC*dsqrt(1d0-4d0*MCC**2/MA**2)
     . *(dsqrt(1d0-(2d0*MD+MPI)**2/MA**2)
     . *(KinAPiDD(MA)-1d0)*((50d0-MA)/(50d0-(2d0*MD+MPI)))**2+1d0)
* End new
       ENDIF
      ENDIF

C       Mixing with the Eta_b(1,2,3S) [1105.1722[hep-ph]]

      METAB1=9.399d0
      METAB2=9.999d0
      METAB3=10.343d0

      DMB1=0.14d0*CB(4)
      DMB2=0.11d0*CB(4)
      DMB3=0.10d0*CB(4)

c      IF(MA.gt.4d0.and.MA.lt.2d0*5.2795d0+MPI)THEN

       MMIX(1,1)=METAB1**2
       MMIX(1,2)=0d0
       MMIX(1,3)=0d0
       MMIX(1,4)=DMB1
       MMIX(2,1)=0d0
       MMIX(2,2)=METAB2**2
       MMIX(2,3)=0d0
       MMIX(2,4)=DMB2
       MMIX(3,1)=0d0
       MMIX(3,2)=0d0
       MMIX(3,3)=METAB3**2
       MMIX(3,4)=DMB3
       MMIX(4,1)=DMB1
       MMIX(4,2)=DMB2
       MMIX(4,3)=DMB3
       MMIX(4,4)=MA**2

       CALL DIAGN(4,MMIX,VALP,VECP,EPS)
       CALL SORTN(4,VALP,VECP)
       DO K=1,4
        DO L=1,4
         OMIX(K,L)=VECP(L,K)
        ENDDO
       ENDDO

       INDA=1
       DO K=1,4
        IF(dabs(OMIX(INDA,4))**2.lt.dabs(OMIX(K,4))**2)INDA=K
       ENDDO
       INDEB1=Max(1,3-INDA)
       DO K=INDEB1,4
        IF(dabs(OMIX(INDEB1,1))**2.lt.dabs(OMIX(K,1))**2)INDEB1=K
       ENDDO

       aux=Max(0d0,min(1d0,MA/2d0-1.5d0))
       GamAhadr=GamAhadr*OMIX(INDA,4)**2+aux*(11.8d-3*OMIX(INDA,1)**2
     .                 +5.4d-3*OMIX(INDA,2)**2+3.9d-3*OMIX(INDA,3)**2)
       GamAee=GamAee*OMIX(INDA,4)**2
       GamAmumu=GamAmumu*OMIX(INDA,4)**2
       GamAtata=GamAtata*OMIX(INDA,4)**2
c Mesonic BRs to 2photons from 1601.05093
       GamAgaga=GamAgaga*OMIX(INDA,4)**2
     .          +aux*(11.8d-3*3.42d-3*OMIX(INDA,1)**2
     .                +5.4d-3*3.38d-3*OMIX(INDA,2)**2
     .                +3.9d-3*3.40d-3*OMIX(INDA,3)**2)
       GamAcc=GamAcc*OMIX(INDA,4)**2
       GamAinv=GamAinv*OMIX(INDA,4)**2
! GamAhadr also contains contributions to the cc final state
c      ENDIF

c Test leptonic eta_b(1S) BR
      WidthEB1mumu=AAmumu**2/(8d0*Pi)*METAB1
     .    *dsqrt(max(0d0,1d0-4d0*MMUON**2/METAB1**2))*OMIX(INDEB1,4)**2
      PROBpill=Max(PROBpill,0.7d0*WidthEB1mumu/(9d-3*11.8d-3)-1d0)
      WidthEB1tata=AAtata**2/(8d0*Pi)*METAB1
     .    *dsqrt(max(0d0,1d0-4d0*MMUON**2/METAB1**2))*OMIX(INDEB1,4)**2
      PROBpill=Max(PROBpill,0.7d0*WidthEB1tata/(8d-2*11.8d-3)-1d0)


C       A -> bb decays

      MBm=5.2795d0

      IF(MA.LE.2d0*MBm+MPI)THEN
       GamAbb= 0d0
      ELSE
       RMB=RUNMB(MA)
       IF(CB(4).ne.0d0)THEN
        RATCOUP= CU(4)/CB(4)
       ELSE
        RATCOUP=0d0
       ENDIF

* New July 2019:
       DBB=CDABS(CJA)**2/(8d0*PI**3)*MA**3
     .     *(ASH**2*AGGQCD2(ASH,5,MA,MT)
     .      -AS4**2*(1d0+AS4/Pi*(97d0/4d0-4d0*7d0/6d0)))

       GamAbb=4d0*(MBP/MA)**2*
     .        3d0*GF*(MBP*CB(4))**2*MA/(4d0*dsqrt(2d0)*Pi)
     .        *TQCDA((MBP/MA)**2)*dsqrt(1d0-4d0*MBP**2/MA**2)
     .         +(1d0-4d0*(MBP/MA)**2)*max(0d0,
     .         3d0*GF*(RMB*CB(4))**2*MA/(4d0*dsqrt(2d0)*Pi)
     .         *QCDA((RMB/MA)**2)*dsqrt(1d0-4d0*RMB**2/MA**2)!+DBB
     .         )
* End new
! Modifying the shape of the phase space function, so that it looks "3-body" at low mass
       IF(MA.le.50d0)THEN
        GamAbb=GamAbb*dsqrt(1d0-(2d0*MBm+MPI)**2/MA**2)
     .                /dsqrt(1d0-4d0*MBP**2/MA**2)
     . *(KinAPiBB(MA)*((50d0-MA)/(50d0-(2d0*MBm+MPI)))**2
     .                 +1d0-((50d0-MA)/(50d0-(2d0*MBm+MPI)))**2)

* New July 2019:
        GamAhadr=GamAhadr+DBB*dsqrt(1d0-4d0*(MBP/MA)**2)
     . *(dsqrt(1d0-(2d0*MBm+MPI)**2/MA**2)
     . *(KinAPiBB(MA)-1d0)*((50d0-MA)/(50d0-(2d0*MBm+MPI)))**2+1d0)
* End new
       ENDIF
      ENDIF

      PROB(65)=PROBpiinv+PROBpigaga+PROBpill

      END

*********************************************************************

      double precision function KinA3Pi(MA1)
      
      implicit none
      double precision MA1,MPI0,aux,f0,f1,f2,f3

      MPI0=0.135d0

      f0=0.024263679d0-0.40366984d0*MA1+0.74504339d0*MA1**2
     . -0.092152445d0*MA1**3+0.018634676d0*MA1**4-0.0015536639d0*MA1**5
      f1=0.24515663d0-1.4645441d0*MA1+2.872038d0*MA1**2
     . -2.3062847d0*MA1**3+1.2069459d0*MA1**4-0.26269506d0*MA1**5
      f2=0.41180543d0-2.9425554d0*MA1+8.1783005d0*MA1**2
     . -11.948466d0*MA1**3+10.074616d0*MA1**4-3.5634090d0*MA1**5
      f3=-0.11370280d0-0.13568551d0*MA1+0.53018965d0*MA1**2
     .   -0.0039151428d0*MA1**3+0.00021232667d0*MA1**4

      if(MA1.lt.3d0*MPI0)then
      aux=0d0
      else
       if(MA1.lt.0.519192d0)then
       aux=f2
       else
        if(MA1.lt.1d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinA3Pi=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAPi3PiC(MA1)
      
      implicit none
      double precision MA1,MPI0,MPIC,aux,f0,f1,f2,f3

      MPI0=0.135d0
      MPIC=0.1396d0

      f0=0.029893595d0-0.42550451d0*MA1+0.75883575d0*MA1**2
     . -0.097429755d0*MA1**3+0.019711936d0*MA1**4-0.0016440058d0*MA1**5
      f1=0.26174772d0-1.5344088d0*MA1+2.9716865d0*MA1**2
     . -2.3894606d0*MA1**3+1.2436537d0*MA1**4-0.26930112d0*MA1**5
      f2=0.43211403d0-3.0248990d0*MA1+8.2469923d0*MA1**2
     . -11.834045d0*MA1**3+9.7967702d0*MA1**4-3.4027148d0*MA1**5
      f3=-0.11582827d0-0.14232105d0*MA1+0.53169789d0*MA1**2
     .   -0.0041128955d0*MA1**3+0.00022312216d0*MA1**4

      if(MA1.lt.MPI0+2d0*MPIC)then
      aux=0d0
      else
       if(MA1.lt.0.526626d0)then
       aux=f2
       else
        if(MA1.lt.1d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAPi3PiC=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAEPi3(MA1)
      
      implicit none
      double precision MA1,MPI0,META,aux,f0,f1,f2,f3

      MPI0=0.135d0
      META=0.548d0

      f0=0.64990732d0-1.7452121d0*MA1+1.3467387d0*MA1**2
     . -0.26834882d0*MA1**3+0.047205187d0*MA1**4-0.0035008773d0*MA1**5
      f1=0.57359739d0-1.2774367d0*MA1+0.38476266d0*MA1**2
     . +0.64516169d0*MA1**3-0.36821687d0*MA1**4+0.070153785d0*MA1**5
      f2=0.34720588d0-0.34301093d0*MA1-1.0057699d0*MA1**2
     . +1.459601d0*MA1**3-0.4300414d0*MA1**4
      f3=0.12068896d0-0.8165721d0*MA1+0.67919395d0*MA1**2
     . -0.023030227d0*MA1**3+0.0012411159d0*MA1**4

      if(MA1.lt.META+2d0*MPI0)then
      aux=0d0
      else
       if(MA1.lt.0.910505d0)then
       aux=f2
       else
        if(MA1.lt.1.3d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAEPi3=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAEPiC(MA1)
      
      implicit none
      double precision MA1,MPIC,META,aux,f0,f1,f2,f3

      MPIC=0.1396d0
      META=0.548d0

      f0=0.67059089d0-1.78635406d0*MA1+1.3721348d0*MA1**2
     . -0.2773546d0*MA1**3+0.048897900d0*MA1**4-0.0036320543d0*MA1**5
      f1=0.60906261d0-1.3723764d0*MA1+0.49128146d0*MA1**2
     . +0.57292304d0*MA1**3-0.34115244d0*MA1**4+0.065852841d0*MA1**5
      f2=0.39620094d0-0.50276489d0*MA1-0.78693086d0*MA1**2
     . +1.3070126d0*MA1**3-0.38795188d0*MA1**4
      f3=0.12498622d0-0.82768627d0*MA1+0.68207647d0*MA1**2
     . -0.023431339d0*MA1**3+0.0012637590d0*MA1**4

      if(MA1.lt.META+2d0*MPIC)then
      aux=0d0
      else
       if(MA1.lt.0.917939d0)then
       aux=f2
       else
        if(MA1.lt.1.3d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAEPiC=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAEPPi3(MA1)
      
      implicit none
      double precision MA1,MPI0,METAP,aux,f0,f1,f2,f3

      MPI0=0.135d0
      METAP=0.958d0

      f0=1.2618863d0-2.0147285d0*MA1+0.78837207d0*MA1**2
     . +0.037033316d0*MA1**3-0.020833739d0*MA1**4+0.0022454179d0*MA1**5
      f1=-1.5349096d0+6.6546959d0*MA1-10.046096d0*MA1**2
     . +6.8581882d0*MA1**3-2.1832176d0*MA1**4+0.27822999d0*MA1**5
      f2=0.34678251d0-0.093908891d0*MA1-0.53691542d0*MA1**2
     . +0.3122347d0*MA1**3
      f3=1.1465540d0-1.9975303d0*MA1+0.90494812d0*MA1**2
     . -0.04948079d0*MA1**3+0.0025742677d0*MA1**4

      if(MA1.lt.METAP+2d0*MPI0)then
      aux=0d0
      else
       if(MA1.lt.1.28909d0)then
       aux=f2
       else
        if(MA1.lt.1.56061d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAEPPi3=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAEPPiC(MA1)
      
      implicit none
      double precision MA1,MPIC,METAP,aux,f0,f1,f2,f3

      MPIC=0.1396d0
      METAP=0.958d0

      f0=1.2950736d0-2.0645192d0*MA1+0.81569259d0*MA1**2
     . +0.028266044d0*MA1**3-0.019309250d0*MA1**4+0.0021341552d0*MA1**5
      f1=-1.4877804d0+6.5312593d0*MA1-9.8872647d0*MA1**2
     . +6.7409858d0*MA1**3-2.1389842d0*MA1**4+0.27158353d0*MA1**5
      f2=0.40142672d0-2.0147285d0*MA1-0.467092288d0*MA1**2
     . +0.29474635d0*MA1**3
      f3=1.1612658d0-2.0155320d0*MA1+0.90987116d0*MA1**2
     . -0.050176372d0*MA1**3+0.0026137474d0*MA1**4

      if(MA1.lt.METAP+2d0*MPIC)then
      aux=0d0
      else
       if(MA1.lt.1.29745d0)then
       aux=f2
       else
        if(MA1.lt.1.56525d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAEPPiC=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAPiEE(MA1)
      
      implicit none
      double precision MA1,MPI0,META,aux,f0,f1,f2,f3

      MPI0=0.135d0
      META=0.548d0

      f0=1.9698720d0-3.6574445d0*MA1+2.1476427d0*MA1**2
     . -0.47073681d0*MA1**3+0.074502976d0*MA1**4-0.0049971571d0*MA1**5
      f1=0.098088398d0+2.4886458d0*MA1-5.9545178d0*MA1**2
     . +4.8885726d0*MA1**3-1.7039829d0*MA1**4+0.23183219d0*MA1**5
      f2=0.34933172d0+1.0673304d0*MA1-3.0579117d0*MA1**2
     . +2.1060809d0*MA1**3-0.41722356d0*MA1**4
      f3=0.86284269d0-1.8343802d0*MA1+0.92840463d0*MA1**2
     . -0.056657617d0*MA1**3+0.0031021849d0*MA1**4

      if(MA1.lt.MPI0+2d0*META)then
      aux=0d0
      else
       if(MA1.lt.1.28273d0)then
       aux=f2
       else
        if(MA1.lt.1.51263d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAPiEE=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAPiEPEP(MA1)
      
      implicit none
      double precision MA1,MPI0,METAP,aux,f0,f1,f2

      MPI0=0.135d0
      METAP=0.958d0

      f0=-0.21328657d0+3.896706d0*MA1-4.8196189d0*MA1**2
     . +2.1435040d0*MA1**3-0.39813691d0*MA1**4+0.028863984d0*MA1**5
      f1=-4.5986145d0+10.970899d0*MA1-8.7327382d0*MA1**2
     . +2.8235570d0*MA1**3-0.31242227d0*MA1**4
      f2=4.8520850d0-5.0859282d0*MA1+1.5625481d0*MA1**2
     . -0.12724282d0*MA1**3+0.0063956343d0*MA1**4

      if(MA1.lt.MPI0+2d0*METAP)then
      aux=0d0
      else
       if(MA1.lt.2.23313d0)then
       aux=f1
       else
        if(MA1.lt.3d0)then
        aux=f0
        else
         if(MA1.lt.4d0)then
          aux=f2
         else
          aux=0d0
         endif
        endif
       endif
      endif

      KinAPiEPEP=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAPiEEP(MA1)
      
      implicit none
      double precision MA1,MPI0,META,METAP,aux,f0,f1,f2,f3

      MPI0=0.135d0
      META=0.548d0
      METAP=0.958d0

      f0=2.8642497d0-3.6569669d0*MA1+1.2735147d0*MA1**2
     . -0.052944707d0*MA1**3-0.010224751d0*MA1**4+0.0015652681d0*MA1**5
      f1=-1.9004084d0+7.8196861d0*MA1-9.8235874d0*MA1**2
     . +5.3303425d0*MA1**3-1.3201466d0*MA1**4+0.12944595d0*MA1**5
      f2=1.4776110d0-1.3890111d0*MA1+0.046790986d0*MA1**2
     . +0.15291947d0*MA1**3
      f3=2.5842822d0-3.4182999d0*MA1+1.2620314d0*MA1**2
     . -0.09746226d0*MA1**3+0.0052029761d0*MA1**4

      if(MA1.lt.MPI0+META+METAP)then
      aux=0d0
      else
       if(MA1.lt.1.72309d0)then
       aux=f2
       else
        if(MA1.lt.2.08792d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAPiEEP=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinA3E(MA1)
      
      implicit none
      double precision MA1,META,aux,f0,f1,f2,f3

      META=0.548d0

      f0=4.931444d0-7.6241437d0*MA1+4.1306376d0*MA1**2
     . -1.0377254d0*MA1**3+0.16048075d0*MA1**4-0.010385007d0*MA1**5
      f1=6.6576464d0-11.610145d0*MA1+7.8289996d0*MA1**2
     . -2.761359d0*MA1**3+0.56398663d0*MA1**4-0.048342641d0*MA1**5
      f2=6.1829128d0-10.077599d0*MA1+5.8730418d0*MA1**2
     . -1.5248800d0*MA1**3+0.17616880d0*MA1**4
      f3=2.2935015d0-3.3511520d0*MA1+1.3349613d0*MA1**2
     . -0.11403480d0*MA1**3+0.0063677382d0*MA1**4

      if(MA1.lt.3d0*META)then
      aux=0d0
      else
       if(MA1.lt.1.72591d0)then
       aux=f2
       else
        if(MA1.lt.2.08993d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinA3E=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAE2EP(MA1)
      
      implicit none
      double precision MA1,META,METAP,aux,f0,f1,f2,f3

      META=0.548d0
      METAP=0.958d0

      f0=8.1577023d0-10.451988d0*MA1+4.8833683d0*MA1**2
     . -1.1287598d0*MA1**3+0.15939030d0*MA1**4-0.0095062710d0*MA1**5
      f1=8.8981257d0-11.956894d0*MA1+6.1087606d0*MA1**2
     . -1.6284664d0*MA1**3+0.26145002d0*MA1**4-0.017858556d0*MA1**5
      f2=6.780322d0-7.5522603d0*MA1+2.5323358d0*MA1**2
     . -0.22522148d0*MA1**3
      f3=6.4360426d0-7.6224031d0*MA1+3.0155995d0*MA1**2
     . -0.50979654d0*MA1**3+0.056413735d0*MA1**4-0.002625951*MA1**5

      if(MA1.lt.2d0*META+METAP)then
      aux=0d0
      else
       if(MA1.lt.2.11114d0)then
       aux=f2
       else
        if(MA1.lt.2.3651d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAE2EP=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAEEP2(MA1)
      
      implicit none
      double precision MA1,META,METAP,aux,f0,f1,f2,f3

      META=0.548d0
      METAP=0.958d0

      f0=12.629192d0-13.959622d0*MA1+5.7978647d0*MA1**2
     . -1.2380980d0*MA1**3+0.15955828d0*MA1**4-0.0087323961d0*MA1**5
      f1=12.762432d0-14.201008d0*MA1+5.9727066d0*MA1**2
     . -1.3013874d0*MA1**3+0.17100697d0*MA1**4-0.0095603275d0*MA1**5
      f2=9.880472158d0-9.191292858d0*MA1+2.578242545d0*MA1**2
     . -0.192945403d0*MA1**3
      f3=11.152317d0-11.612590d0*MA1+4.3026519d0*MA1**2
     . -0.76073205d0*MA1**3+0.083173891d0*MA1**4-0.0038313683*MA1**5

      if(MA1.lt.2d0*METAP+META)then
      aux=0d0
      else
       if(MA1.lt.2.49638d0)then
       aux=f2
       else
        if(MA1.lt.2.64027d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAEEP2=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinA3EP(MA1)
      
      implicit none
      double precision MA1,METAP,aux,f0,f1,f2

      METAP=0.958d0

      f0=22.5681998d0-23.854110d0*MA1+10.109227d0*MA1**2
     . -2.3108416d0*MA1**3+0.30188926d0*MA1**4-0.016586127d0*MA1**5
      f1=19.217308058d0-18.058918599d0*MA1+6.1003309745d0*MA1**2
     . -0.92427222990d0*MA1**3+0.062106313474d0*MA1**4
      f2=22.756613d0-24.667569d0*MA1+11.129349d0*MA1**2
     . -2.9074938d0*MA1**3+0.48483710d0*MA1**4-0.045235221d0*MA1**5
     . +0.0018154043d0*MA1**6

      if(MA1.lt.3d0*METAP)then
      aux=0d0
      else
       if(MA1.lt.2.91544d0)then
       aux=f1
       else
        if(MA1.lt.3d0)then
        aux=f0
        else
         if(MA1.lt.4d0)then
          aux=f2
         else
          aux=0d0
         endif
        endif
       endif
      endif

      KinA3EP=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAPiKC(MA1)
      
      implicit none
      double precision MA1,MPI0,MKC,aux,f0,f1,f2,f3

      MPI0=0.135d0
      MKC=0.494d0

      f0=1.6184070d0-3.2969036d0*MA1+2.1105929d0*MA1**2
     . -0.49440392d0*MA1**3+0.083445992d0*MA1**4-0.0059305000d0*MA1**5
      f1=0.80123212d0-0.349357371d0*MA1-2.1228457d0*MA1**2
     . +2.5366896d0*MA1**3-0.99969151d0*MA1**4+0.14874120d0*MA1**5
      f2=0.6984555428d0-0.3068607074d0*MA1-1.544771723d0*MA1**2
     . +1.508426459d0*MA1**3-0.3407850274d0*MA1**4
      f3=0.55238932d0-1.4779014d0*MA1+0.84479081d0*MA1**2
     . -0.045640920d0*MA1**3+0.0025019861d0*MA1**4

      if(MA1.lt.MPI0+2d0*MKC)then
      aux=0d0
      else
       if(MA1.lt.1.18455d0)then
       aux=f2
       else
        if(MA1.lt.1.45808d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAPiKC=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAPiK0(MA1)
      
      implicit none
      double precision MA1,MPI0,MK0,aux,f0,f1,f2,f3

      MPI0=0.135d0
      MK0=0.498d0

      f0=1.6454207d0-3.3285868d0*MA1+2.1180662d0*MA1**2
     . -0.49460758d0*MA1**3+0.083176117d0*MA1**4-0.0058925040d0*MA1**5
      f1=0.76465544d0-0.18166273d0*MA1-2.3648306d0*MA1**2
     . +2.6919416d0*MA1**3-1.0481136d0*MA1**4+0.15469361d0*MA1**5
      f2=1.6184070d0-3.2969036d0*MA1+2.1105929d0*MA1**2
     . -0.49440392d0*MA1**3+0.083445992d0*MA1**4-0.0059305000d0*MA1**5
      f3=0.67816973d0-0.21518244d0*MA1-1.6542307d0*MA1**2
     . +1.5548786d0*MA1**3-0.34729486d0*MA1**4

      if(MA1.lt.MPI0+2d0*MK0)then
      aux=0d0
      else
       if(MA1.lt.1.19182d0)then
       aux=f2
       else
        if(MA1.lt.1.46212d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAPiK0=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAPiKcK0(MA1)
      
      implicit none
      double precision MA1,MPIC,MKC,MK0,aux,f0,f1,f2,f3

      MPIC=0.1396d0
      MKC=0.494d0
      MK0=0.498d0

      f0=1.6542043d0-3.3524966d0*MA1+2.1405766d0*MA1**2
     . -0.50400710d0*MA1**3+0.085100830d0*MA1**4-0.0060492292d0*MA1**5
      f1=0.88259317d0-0.56911161d0*MA1-1.8534525d0*MA1**2
     . +2.3509667d0*MA1**3-0.93284682d0*MA1**4+0.13892869d0*MA1**5
      f2=0.79820142d0-0.57758071d0*MA1-1.2413712d0*MA1**2
     . +1.3437664d0*MA1**3-0.30625602d0*MA1**4
      f3=0.56765794d0-1.4980334d0*MA1+0.84983877d0*MA1**2
     . -0.046331854d0*MA1**3+0.0025406034d0*MA1**4

      if(MA1.lt.MPIC+MKC+MK0)then
      aux=0d0
      else
       if(MA1.lt.1.19236d0)then
       aux=f2
       else
        if(MA1.lt.1.46242d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAPiKcK0=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAEKC(MA1)
      
      implicit none
      double precision MA1,META,MKC,aux,f0,f1,f2,f3

      META=0.548d0
      MKC=0.494d0

      f0=4.0378221d0-6.5656473d0*MA1+3.6672310d0*MA1**2
     . -0.91894486d0*MA1**3+0.14426331d0*MA1**4-0.0094731738d0*MA1**5
      f1=5.6695386d0-10.463516d0*MA1+7.4126237d0*MA1**2
     . -2.7285844d0*MA1**3+0.58391277d0*MA1**4-0.052436215d0*MA1**5
      f2=5.348836142d0-9.303973989d0*MA1+5.775979320d0*MA1**2
     . -1.594182292d0*MA1**3+0.1961717071d0*MA1**4
      f3=1.7516996d0-2.8309136d0*MA1+1.2004068d0*MA1**2
     . -0.095396951d0*MA1**3+0.0053196349d0*MA1**4

      if(MA1.lt.META+2d0*MKC)then
      aux=0d0
      else
       if(MA1.lt.1.62443d0)then
       aux=f2
       else
        if(MA1.lt.2.01745d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAEKC=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAEK0(MA1)
      
      implicit none
      double precision MA1,META,MK0,aux,f0,f1,f2,f3

      META=0.548d0
      MK0=0.498d0

      f0=4.1009729d0-6.6424338d0*MA1+3.7015327d0*MA1**2
     . -0.92791299d0*MA1**3+0.14551367d0*MA1**4-0.0095451364d0*MA1**5
      f1=5.7419986d0-10.552553d0*MA1+7.4488077d0*MA1**2
     . -2.7335440d0*MA1**3+0.58296173d0*MA1**4-0.052170109d0*MA1**5
      f2=5.410569593d0-9.365543879d0*MA1+5.787008605d0*MA1**2
     . -1.590089333d0*MA1**3+0.1947455514d0*MA1**4
      f3=1.7887610d0-2.8674416d0*MA1+1.2098626d0*MA1**2
     . -0.096707931d0*MA1**3+0.0053934134d0*MA1**4

      if(MA1.lt.META+2d0*MK0)then
      aux=0d0
      else
       if(MA1.lt.1.63195d0)then
       aux=f2
       else
        if(MA1.lt.2.02282d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAEK0=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAEPKC(MA1)
      
      implicit none
      double precision MA1,METAP,MKC,aux,f0,f1,f2,f3

      METAP=0.958d0
      MKC=0.494d0

      f0=6.7874942d0-8.9513768d0*MA1+4.2149725d0*MA1**2
     . -0.95613039d0*MA1**3+0.13546848d0*MA1**4-0.0081235587d0*MA1**5
      f1=7.2607734d0-9.9228439d0*MA1+5.0134778d0*MA1**2
     . -1.2846903d0*MA1**3+0.20315081d0*MA1**4-0.013708078d0*MA1**5
      f2=6.82817926d0-8.838795311d0*MA1+3.926626571d0*MA1**2
     . -0.7397461805d0*MA1**3+0.06650646595d0*MA1**4
      f3=5.3444409d0-6.5686216d0*MA1+2.6341031d0*MA1**2
     . -0.42933370d0*MA1**3+0.047302253d0*MA1**4-0.0021953075d0*MA1**5

      if(MA1.lt.METAP+2d0*MKC)then
      aux=0d0
      else
       if(MA1.lt.2.00966d0)then
       aux=f2
       else
        if(MA1.lt.2.29262d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAEPKC=max(0d0,aux)

      return

      end

*********************************************************************

      double precision function KinAEPK0(MA1)
      
      implicit none
      double precision MA1,METAP,MK0,aux,f0,f1,f2,f3

      METAP=0.958d0
      MK0=0.498d0

      f0=6.8838000d0-9.0584799d0*MA1+4.2631525d0*MA1**2
     . -0.96870538d0*MA1**3+0.13723064d0*MA1**4-0.0082266343d0*MA1**5
      f1=7.3788298d0-10.074798d0*MA1+5.0987964d0*MA1**2
     . -1.3127033d0*MA1**3+0.20813483d0*MA1**4-0.014081224d0*MA1**5
      f2=6.929587159d0-8.951366199d0*MA1+3.974857529d0*MA1**2
     . -0.7503920953d0*MA1**3+0.06745018995d0*MA1**4
      f3=5.4202542d0-6.6426440d0*MA1+2.6608901d0*MA1**2
     . -0.43498190d0*MA1**3+0.047941703d0*MA1**4-0.0022255247d0*MA1**5

      if(MA1.lt.METAP+2d0*MK0)then
      aux=0d0
      else
       if(MA1.lt.2.01718d0)then
       aux=f2
       else
        if(MA1.lt.2.29799d0)then
        aux=f1
        else
         if(MA1.lt.3d0)then
          aux=f0
         else
          if(MA1.lt.4d0)then
           aux=f3
          else
           aux=0d0
          endif
         endif
        endif
       endif
      endif

      KinAEPK0=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinARhogam(MA1)
      
      implicit none
      double precision MA1,MPIC,aux,f0,f1,f2,f3,f4,f5

      MPIC=0.1396d0

      f0=-0.18242916d0+0.68164651d0*MA1-0.85152503d0*MA1**2
     . +0.38260017d0*MA1**3-8.6998880d-2*MA1**4+7.5335087d-2*MA1**5
     . -1.0076012d-2*MA1**6+9.8901154d-4*MA1**7+1.9707307d-4*MA1**8
      f1=-5.590472138d0+54.85645674d0*MA1-233.2654325d0*MA1**2
     . +561.5186568d0*MA1**3-837.2576233d0*MA1**4+792.4264267d0*MA1**5
     . -465.4373957d0*MA1**6+155.3287077d0*MA1**7-22.56958544d0*MA1**8
      f2=7.700949045d-2-1.100798389d0*MA1+6.886613070d0*MA1**2
     . -24.63449002d0*MA1**3+55.13018382d0*MA1**4-79.07367195d0*MA1**5
     . +71.02603065d0*MA1**6-36.55738269d0*MA1**7+8.265684609d0*MA1**8
      f3=1.108073418d-5-2.601297711d-4*MA1+2.687574905d-3*MA1**2
     . -1.593264879d-2*MA1**3+5.917005055d-2*MA1**4
     . -1.408698209d-1*MA1**5+2.104861084d-1*MA1**6
     . -1.820515767d-1*MA1**7+7.109517316d-2*MA1**8
      f4=-2.457776218207d-5+6.026288437041d-4*MA1
     . -6.452738235925d-3*MA1**2+3.945732240042d-2*MA1**3
     . -1.508583185681d-1*MA1**4+3.694614492014d-1*MA1**5
     . -5.655714822692d-1*MA1**6+4.932756176607d-1*MA1**7
     . -1.863868994640d-1*MA1**8
      f5=41.968136d0-82.383631d0*MA1+68.635980d0*MA1**2
     . -31.455760d0*MA1**3+8.4604253d0*MA1**4-1.2374164d0*MA1**5
     . +0.087875503d0*MA1**6

      if(MA1.lt.2d0*MPIC)then
      aux=0d0
      else
       if(MA1.lt.0.344727d0)then
       aux=f4
       else
        if(MA1.lt.0.490343d0)then
        aux=f3
        else
         if(MA1.lt.0.708768d0)then
         aux=f2
         else
          if(MA1.lt.1d0)then
          aux=f1
          else
           if(MA1.lt.3d0)then
           aux=f0
           else
            if(MA1.lt.4d0)then
            aux=f5
            else
            aux=0d0
            endif
           endif
          endif
         endif
        endif
       endif
      endif

      KinARhogam=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAPIDD(MA1)
      
      implicit none
      double precision MA1,MPI,MD,aux,f0,f1,f2,f3

      MPI=0.135d0
      MD=1.86483d0

      f0=-1.3326448d0+0.42254184d0*MA1-0.038101705d0*MA1**2
     . +0.0019325144d0*MA1**3-5.2692088d-5*MA1**4+6.0301076d-7*MA1**5
      f1=7.3769868d0-6.7131360d0*MA1+2.5139293d0*MA1**2
     . -0.51302330d0*MA1**3+0.063079400d0*MA1**4-0.0046893174d0*MA1**5
     . +1.9501589d-4*MA1**6-3.4945281d-6*MA1**7
      f2=-116.793482669d0+165.515661567d0*MA1-97.2345923407d0*MA1**2
     . +30.3600771681d0*MA1**3-5.32241364118d0*MA1**4
     . +0.497391376090d0*MA1**5-0.0193725847592d0*MA1**6
      f3=-0.44076166d0+0.17647968d0*MA1-0.010473319d0*MA1**2
     . +0.00035401583d0*MA1**3-6.9392524d-6*MA1**4
     . +7.3503518d-8*MA1**5-3.2564926d-10*MA1**6

      if(MA1.lt.MPI+2d0*MD)then
      aux=0d0
      else
       if(MA1.lt.4.36827d0)then
       aux=f2
       else
        if(MA1.lt.9.40435d0)then
        aux=f1
        else
         if(MA1.lt.15d0)then
          aux=f0
         else
          if(MA1.lt.50d0)then
           aux=f3
          else
           aux=1d0
          endif
         endif
        endif
       endif
      endif

      KinAPiDD=max(0d0,aux)

      return
      end

*********************************************************************

      double precision function KinAPIBB(MA1)
      
      implicit none
      double precision MA1,MPI,MB,aux,f0,f1,f2,f3

      MPI=0.135d0
      MB=5.2795d0

      f0=4.5191789d0-1.2840862d0*MA1+0.13868018d0*MA1**2
     . -7.2155875d-3*MA1**3+1.8761613d-4*MA1**4-1.9645580d-6*MA1**5
      f1=-1.3461087d0+0.68217702d0*MA1-0.11195317d0*MA1**2
     . +0.0074532096d0*MA1**3-1.7287542d-4*MA1**4
      f2=0.86212930d0-0.16107566d0*MA1+0.0075236223d0*MA1**2
      f3=-1.0071683d0+0.095212536d0*MA1-1.0000702d-3*MA1**2
     . -5.2877368d-5*MA1**3+2.0062512d-6*MA1**4
     . -2.7213842d-8*MA1**5+1.3590193d-10*MA1**6

      if(MA1.lt.MPI+2d0*MB)then
      aux=0d0
      else
       if(MA1.lt.11.089d0)then
       aux=f2
       else
        if(MA1.lt.12.4717d0)then
        aux=f1
        else
         if(MA1.lt.20.3724d0)then
          aux=f0
         else
          if(MA1.lt.50d0)then
           aux=f3
          else
           aux=1d0
          endif
         endif
        endif
       endif
      endif

      KinAPiBB=max(0d0,aux)

      return
      end



