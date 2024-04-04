      SUBROUTINE HDECAY_lightH(PAR)

      IMPLICIT NONE

      INTEGER K,L,N0,VFLAG,INDH,NF

      DOUBLE PRECISION PAR(*),MH,EPS,PI,aux,HF,H1,H2,SINB,COSB,GHCC
      DOUBLE PRECISION BETA,X,SP,QCD0,HQCDM,HQCD,QCDH,TQCDH,RATCOUP
      DOUBLE PRECISION FPI,THETAE,MPI,MPIC,MPI8,MPI9,META,METAP,MKc,MK0
      DOUBLE PRECISION CRUL,CRDL,CRUR,CRDR,CRLL,CRLR,CRULR,CRDLR,CRLLR
      DOUBLE PRECISION GHCHACHA(2,2),MS2,GamS2,MS3,GamS3,MSIG,GamSIG
      DOUBLE PRECISION AHee,AHmumu,AHtata,MG0,GamG0
      DOUBLE PRECISION RMS,ASH,ASH2,AS3,AS4,HIGTOP,ALPHAS,RUNM
      DOUBLE PRECISION AH2Pi,AHPiE,GamH2Pi,GamH2PiC,GamHPiE,GamHPiEP
      DOUBLE PRECISION AH2K,GamH2KC,GamH2K0,AH2E,GamH2E,GamHEEP,GamH2EP
      DOUBLE PRECISION AHuu,AHdd,AHss,GamHuu,GamHdd,GamHss,GamHjj
      DOUBLE PRECISION MD,DCC,RMC,MBm,DBB,RMB,runmb,HTWW,HTZZ
      DOUBLE PRECISION MCHIC,DMC,CMIX,MCHIB1,MCHIB2,DMB1,DMB2
      DOUBLE PRECISION MMIX(3,3),VALP(3),VECP(3,3),OMIX(3,3)
      DOUBLE PRECISION ALEM0
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION LAMBDA,KAPPA,ALAMBDA,AKAPPA,MUEFF,NUEFF
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MCC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MPI0,MEL,MSTRANGE
      DOUBLE PRECISION XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION HUQ,HDQ,MTQ,MBQ
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION GamHGAGA,GamHee,GamHmumu,GamHtata,GamHhadr,
     . GamHcc,GamHbb,GamHinv(3),GamHWW,GamHZZ,GamHAA
      DOUBLE PRECISION DDCOS,DDSIN
* New May 2019:
      DOUBLE PRECISION ZETA2,ZETA3,ASG,HGGQCD2,CIH
* End New

      DOUBLE COMPLEX XC,TC,BC,CC,LC,MC,EC,CH1C,CH2C,WC,HC,PropI,PropII
      DOUBLE COMPLEX ULC,URC,DLC,DRC,T1C,T2C,B1C,B2C,LLC,LRC,L1C,L2C
      DOUBLE COMPLEX CJH,CGH,F0,FF,FS,FV,CM,CJHF

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
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QQUARK/HUQ,HDQ,MTQ,MBQ
      COMMON/QEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/STSBSCALE/QSTSB
      COMMON/VFLAG/VFLAG
      COMMON/LIGHTHDECAYS/GamHGAGA,GamHee,GamHmumu,GamHtata,GamHhadr,
     . GamHcc,GamHbb,GamHinv,GamHWW,GamHZZ,GamHAA

      CM(X)= DCMPLX(MIN(1d3,X)**2,-EPS/4d0)
      F0(XC)=-CDLOG((CDSQRT(1d0-4d0*XC)-1d0)
     .               /(CDSQRT(1d0-4d0*XC)+1d0))**2/4d0
      FF(XC)=8d0*XC*(1d0+(1d0-4d0*XC)*F0(XC))
      FS(XC)=2d0*XC*(4d0*XC*F0(XC)-1d0)
      FV(XC)=-(2d0+12d0*XC+24d0*XC*(1d0-2d0*XC)*F0(XC))
      CRUL(HF)= dsqrt(2d0)*(HF**2*H1Q*S(1,1)
     . + (g1/12d0-g2/4d0)*(H1Q*S(1,1)-H2Q*S(1,2)))
      CRDL(HF)= dsqrt(2d0)*(HF**2*H2Q*S(1,2)
     . + (g1/12d0+g2/4d0)*(H1Q*S(1,1)-H2Q*S(1,2)))
      CRUR(HF)= dsqrt(2d0)*(HF**2*H1Q*S(1,1)
     . - g1/3d0*(H1Q*S(1,1)-H2Q*S(1,2)))
      CRDR(HF)= dsqrt(2d0)*(HF**2*H2Q*S(1,2)
     . + g1/6d0*(H1Q*S(1,1)-H2Q*S(1,2)))
      CRLL(HF)= dsqrt(2d0)*(HF**2*H2Q*S(1,2)
     .      + (-g1/4d0+g2/4d0)*(H1Q*S(1,1)-H2Q*S(1,2)))
      CRLR(HF)= dsqrt(2d0)*(HF**2*H2Q*S(1,2)
     .                + g1/2d0*(H1Q*S(1,1)-H2Q*S(1,2)))
      BETA(X)= DSQRT(1d0-4d0*X)
      QCd0(X)= (1d0+X**2)*(4d0*SP((1d0-X)/(1d0+X))
     . +2d0*SP((X-1d0)/(X+1d0))
     . - 3d0*DLOG((1d0+X)/(1d0-X))*DLOG(2d0/(1d0+X))
     . - 2d0*DLOG((1d0+X)/(1d0-X))*DLOG(X))
     . - 3d0*X*DLOG(4d0/(1d0-X**2))-4d0*X*DLOG(X)
      HQCDM(X)= QCd0(X)/X+(3d0+34d0*X**2-13d0*X**4)/16d0/X**3
     . * DLOG((1d0+X)/(1d0-X))+3d0/8d0/X**2*(7d0*X**2-1d0)
* July 2010:
c      HQCD(X)=(4d0/3d0*HQCDM(BETA(X))
c     .   +2d0*(4d0/3d0-DLOG(X))*(1d0-10d0*X)/(1d0-4d0*X))*ASH/PI
c     .   + (29.14671d0 + RATCOUP*(1.570d0 - 2d0*DLOG(HIGTOP)/3d0
c     .   + DLOG(X)**2/9d0))*(ASH/PI)**2
c     .   +(164.14d0 - 25.77d0*5 + 0.259d0*5**2)*(ASH/PI)**3
* New May 2019:
      HQCD(X)=(4d0/3d0*HQCDM(BETA(X))
     .       + 2d0*(4d0/3d0-DLOG(X))*(1d0-10d0*X)/(1d0-4d0*X))*ASH/PI
     .       + (29.14671d0
     .         + RATCOUP*(1.570d0 - 2d0*DLOG(HIGTOP)/3d0
     .            + DLOG(X)**2/9d0))*(ASH/PI)**2
     .       + (164.14d0 - 25.77d0*5d0 + 0.259d0*5d0**2)*(ASH/PI)**3
     .       + (39.34d0-220.9d0*5d0+9.685d0*5d0**2
     .         - 0.0205d0*5d0**3)*(ASH/PI)**4
* End New
      QCDH(X)= 1d0+HQCD(X)
      TQCDH(X)= 1d0+4d0/3d0*HQCDM(BETA(X))*ASH/PI
* New May 2019:
      HGGQCD2(ASG,NF,MH,MT)= 1d0+ASG/PI*(95d0/4d0-NF*7d0/6d0)
     . +(ASG/PI)**2*(149533d0/288d0-363d0/8d0*ZETA2-495d0/8d0*ZETA3
     .              +19d0/8d0*DLOG(MH**2/MT**2)
     . +NF*(-4157d0/72d0+11d0/2d0*ZETA2+5d0/4d0*ZETA3
     . +2d0/3d0*DLOG(MH**2/MT**2))
     . +NF**2*(127d0/108d0-1d0/6d0*ZETA2))+(ASG/PI)**3
     . *(467.683620788d0+122.440972222d0*DLOG(MH**2/MT**2)
     .              +10.9409722222d0*DLOG(MH**2/MT**2)**2)
* End New

      EPS=1d-8
      PI=4d0*DATAN(1d0)
* New May 2019:
      ZETA2=PI**2/6d0
      ZETA3=1.202056903159594d0
* End New

      MH=SMASS(1)
      COSB=1d0/DSQRT(1d0+PAR(3)**2)
      SINB=PAR(3)*COSB
      H1=SINB/DSQRT(2d0*dsqrt(2d0)*GF)
      H2=COSB/DSQRT(2d0*dsqrt(2d0)*GF)

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

      TC=CM(MT/MH)
      BC=CM(MBP/MH)
      CC=CM(MCC/MH)
      LC=CM(MTAU/MH)
      MC=CM(MMUON/MH)
      EC=CM(MEL/MH)
      CH1C=CM(MCHA(1)/MH)
      CH2C=CM(MCHA(2)/MH)
      WC=CM(MW/MH)
      HC=CM(CMASS/MH)
      ULC=CM(MUL/MH)
      URC=CM(MUR/MH)
      DLC=CM(MDL/MH)
      DRC=CM(MDR/MH)
      T1C=CM(MST1/MH)
      T2C=CM(MST2/MH)
      B1C=CM(MSB1/MH)
      B2C=CM(MSB2/MH)
      LLC=CM(MLL/MH)
      LRC=CM(MLR/MH)
      L1C=CM(MSL1/MH)
      L2C=CM(MSL2/MH)

      CRULR=HUQ/dsqrt(2d0)*(PAR(12)*S(1,1)-MUEFF*S(1,2)
     .                            -LAMBDA*H2Q*S(1,3))
      CRDLR=HDQ/dsqrt(2d0)*(-MUEFF*S(1,1)+PAR(13)*S(1,2)
     .                            -LAMBDA*H1Q*S(1,3))
      CRLLR= MTAU/H2Q/dsqrt(2d0)*(-MUEFF*S(1,1)+PAR(14)*S(1,2)
     .                                      - LAMBDA*H1Q*S(1,3))

      GHCC=LAMBDA*dsqrt(2d0)*(MUEFF*S(1,3)
     .                 -LAMBDA*SINB*COSB*(H1*S(1,2)+H2*S(1,1)))
     . +MUEFF*KAPPA*2d0*dsqrt(2d0)*S(1,3)*SINB*COSB
     . +LAMBDA*ALAMBDA*dsqrt(2d0)*S(1,3)*SINB*COSB
     . +g1/(2d0*dsqrt(2d0))*(H1*S(1,1)*(COSB**2-SINB**2)
     .                      +H2*S(1,2)*(SINB**2-COSB**2))
     . +g2/(2d0*dsqrt(2d0))*(H1*(S(1,1)+2d0*S(1,2)*SINB*COSB)
     .                      +H2*(S(1,2)+2d0*S(1,1)*SINB*COSB))
     . +LAMBDA*MUP*S(1,3)*SINB*COSB*dsqrt(2d0)
     . +6d0*dsqrt(2d0)*DLOG(MAX(QSTSB,MH**2)/(MAX(MT,MH)**2))
     .     *(MT**4/H1**3*S(1,1)*COSB**2+MB**4/H2**3*S(1,2)*SINB**2
     .      +MT**2*MB**2/(H1**2*H2**2)
     .        *(H1*S(1,1)*SINB**2+H2*S(1,1)*SINB*COSB
     .         +H2*S(1,2)*COSB**2+H1*S(1,2)*SINB*COSB))/(16d0*PI**2)

      DO K=1,2
      DO L=1,2
       GHCHACHA(K,L)=LAMBDA/dsqrt(2d0)*S(1,3)*U(K,2)*V(L,2)
     .    +dsqrt(g2/2d0)*(S(1,1)*U(K,1)*V(L,2)+S(1,2)*U(K,2)*V(L,1))
      ENDDO
      ENDDO

      CJHF=DSQRT(dsqrt(2d0)*GF)/4d0*(CU(1)*(FF(TC)+FF(CC))+CB(1)*FF(BC))

      CJH=CJHF
     .  +(CRUL(0d0)*FS(ULC)/MUL**2+CRUR(0d0)*FS(URC)/MUR**2
     .   +CRDL(0d0)*FS(DLC)/MDL**2+CRDR(0d0)*FS(DRC)/MDR**2)/2d0
     .  +(CST**2*CRUL(HUQ)+(1d0-CST**2)*CRUR(HUQ)+2d0*CST
     .                  *DSQRT(1d0-CST**2)*CRULR)*FS(T1C)/(4d0*MST1**2)
     .  +((1d0-CST**2)*CRUL(HUQ)+CST**2*CRUR(HUQ)-2d0*CST
     .                  *DSQRT(1d0-CST**2)*CRULR)*FS(T2C)/(4d0*MST2**2)
     .  +(CSB**2*CRDL(HDQ)+(1d0-CSB**2)*CRDR(HDQ)+2d0*CSB
     .                  *DSQRT(1d0-CSB**2)*CRDLR)*FS(B1C)/(4d0*MSB1**2)
     .  +((1d0-CSB**2)*CRDL(HDQ)+CSB**2*CRDR(HDQ)-2d0*CSB
     .                  *DSQRT(1d0-CSB**2)*CRDLR)*FS(B2C)/(4d0*MSB2**2)

* New May 2019:
      CIH=DREAL(DCONJG(CJH)*(CJH-CJHF))
* End New

      CGH=DSQRT(dsqrt(2d0)*GF)/2d0*(4d0/3d0*CU(1)*(FF(TC)+FF(CC))
     .    +CB(1)*FF(BC)/3d0+CD(1)*(FF(LC)+FF(MC)+FF(EC))+CV(1)*FV(WC))
     .  +GHCC/CMASS**2*FS(HC)/2d0
     .  +GHCHACHA(1,1)/MCHA(1)*FF(CH1C)/2d0
     .  +GHCHACHA(2,2)/MCHA(2)*FF(CH2C)/2d0
     .  +2d0/3d0*((CST**2*CRUL(HUQ)+(1d0-CST**2)*CRUR(HUQ)+2d0*CST
     .                     *DSQRT(1d0-CST**2)*CRULR)*FS(T1C)/MST1**2
     .           +((1d0-CST**2)*CRUL(HUQ)+CST**2*CRUR(HUQ)-2d0*CST
     .                     *DSQRT(1d0-CST**2)*CRULR)*FS(T2C)/MST2**2
     .    +2d0*CRUL(0d0)*FS(ULC)/MUL**2+2d0*CRUR(0d0)*FS(URC)/MUR**2)
     .  +1d0/6d0*((CSB**2*CRDL(HDQ)+(1d0-CSB**2)*CRDR(HDQ)+2d0*CSB
     .                     *DSQRT(1d0-CSB**2)*CRDLR)*FS(B1C)/MSB1**2
     .           +((1d0-CSB**2)*CRDL(HDQ)+CSB**2*CRDR(HDQ)-2d0*CSB
     .                     *DSQRT(1d0-CSB**2)*CRDLR)*FS(B2C)/MSB2**2
     .    +2d0*CRDL(0d0)*FS(DLC)/MDL**2+2d0*CRDR(0d0)*FS(DRC)/MDR**2)
     .  +CRLL(0d0)*FS(LLC)/MLL**2+CRLR(0d0)*FS(LRC)/MLR**2
     .  +(CSL**2*CRLL(MTAU/H2Q)+(1d0-CSL**2)*CRLR(MTAU/H2Q)+2d0*CSL
     .                   *DSQRT(1d0-CSL**2)*CRLLR)*FS(L1C)/MSL1**2/2d0
     .  +((1d0-CSL**2)*CRLL(MTAU/H2Q)+CSL**2*CRLR(MTAU/H2Q)-2d0*CSL
     .                   *DSQRT(1d0-CSL**2)*CRLLR)*FS(L2C)/MSL2**2/2d0

C       Leptonic decays

C  * H -> ee
      AHee=MEL*CD(1)*DSQRT(dsqrt(2d0)*GF)
      GamHee=AHee**2/(8d0*Pi)*MH
     .               *dsqrt(max(0d0,1d0-4d0*MEL**2/MH**2))**3

C  * H -> mumu
      AHmumu=MMUON*CD(1)*DSQRT(dsqrt(2d0)*GF)
      GamHmumu=AHmumu**2/(8d0*Pi)*MH
     .               *dsqrt(max(0d0,1d0-4d0*MMUON**2/MH**2))**3

C  * H -> tautau
      AHtata=MTAU*CL(1)*DSQRT(dsqrt(2d0)*GF)
      GamHtata=AHtata**2/(8d0*Pi)*MH
     .               *dsqrt(max(0d0,1d0-4d0*MTAU**2/MH**2))**3

c       Decay to light neutralinos

      IF(MH.LE.2d0*DABS(MNEU(1)))THEN
       GamHinv(1)=0d0
      ELSE
       aux=dsqrt(2d0)*LAMBDA*(S(1,1)*N(1,4)*N(1,5)
     .                       +S(1,2)*N(1,3)*N(1,5)
     .                       +S(1,3)*N(1,3)*N(1,4))
     .    -dsqrt(2d0)*KAPPA*S(1,3)*N(1,5)*N(1,5)
     .    +dsqrt(g1)*(-S(1,1)*N(1,1)*N(1,3)+S(1,2)*N(1,1)*N(1,4))
     .    +dsqrt(g2)*(S(1,1)*N(1,2)*N(1,3)-S(1,2)*N(1,2)*N(1,4))

       GamHinv(1)=aux**2/(16d0*PI)*MH*dsqrt(1d0-4d0*(MNEU(1)/MH)**2)**3
      ENDIF

      IF(MH.LE.DABS(MNEU(1))+DABS(MNEU(2)))THEN
       GamHinv(2)=0d0
      ELSE
       aux=LAMBDA*(S(1,1)*(N(1,4)*N(2,5)+N(2,4)*N(1,5))
     .            +S(1,2)*(N(1,3)*N(2,5)+N(2,3)*N(1,5))
     .            +S(1,3)*(N(1,3)*N(2,4)+N(2,3)*N(1,4)))/dsqrt(2d0)
     .    -dsqrt(2d0)*KAPPA*S(1,3)*N(1,5)*N(2,5)
     .    +dsqrt(g1)/2d0*(-S(1,1)*(N(1,1)*N(2,3)+N(2,1)*N(1,3))
     .                    +S(1,2)*(N(1,1)*N(2,4)+N(2,1)*N(1,4)))
     .    +dsqrt(g2)/2d0*(S(1,1)*(N(1,2)*N(2,3)+N(2,2)*N(1,3))
     .                   -S(1,2)*(N(1,2)*N(2,4)+N(2,2)*N(1,4)))

       GamHinv(2)=aux**2/(8d0*PI)*MH
     .               *dsqrt(1d0-((MNEU(1)+MNEU(2))/MH)**2)**3
     .               *dsqrt(1d0-((MNEU(1)-MNEU(2))/MH)**2)
      ENDIF

      IF(MH.LE.2d0*DABS(MNEU(2)))THEN
       GamHinv(3)=0d0
      ELSE
       aux=dsqrt(2d0)*LAMBDA*(S(1,1)*N(2,4)*N(2,5)
     .                       +S(1,2)*N(2,3)*N(2,5)
     .                       +S(1,3)*N(2,3)*N(2,4))
     .    -dsqrt(2d0)*KAPPA*S(1,3)*N(2,5)*N(2,5)
     .    +dsqrt(g1)*(-S(1,1)*N(2,1)*N(2,3)+S(1,2)*N(2,1)*N(2,4))
     .    +dsqrt(g2)*(S(1,1)*N(2,2)*N(2,3)-S(1,2)*N(2,2)*N(2,4))

       GamHinv(3)=aux**2/(16d0*PI)*MH*dsqrt(1d0-4d0*(MNEU(2)/MH)**2)**3
      ENDIF

c       Decay to light pseudoscalars

       IF(MH.LE.2d0*PMASS(1))THEN
        GamHAA=0d0
       ELSE
      aux=(g1+g2)/(2d0*dsqrt(2d0))
     . *P2(1,1)**2*(H1*S(1,1)*(COSB**2-SINB**2)
     .             +H2*S(1,2)*(SINB**2-COSB**2))
     . +LAMBDA*ALAMBDA*dsqrt(2d0)*(S(1,1)*P2(1,1)*P2(1,2)*SINB
     .   +S(1,2)*P2(1,1)*P2(1,2)*COSB+S(1,3)*P2(1,1)**2*SINB*COSB)
     . -KAPPA*AKAPPA*dsqrt(2d0)*S(1,3)*P2(1,2)**2
     . +LAMBDA**2*dsqrt(2d0)*(H1*S(1,1)*(P2(1,1)**2*SINB**2+P2(1,2)**2)
     .    +H2*S(1,2)*(P2(1,1)**2*COSB**2+P2(1,2)**2)
     .    +MUEFF/LAMBDA*S(1,3)*P2(1,1)**2)
     . +KAPPA**2*2d0*dsqrt(2d0)*MUEFF/LAMBDA*S(1,3)*P2(1,2)**2
     . +LAMBDA*KAPPA*dsqrt(2d0)
     .   *(H1*(S(1,2)*P2(1,2)**2-2d0*S(1,3)*P2(1,1)*P2(1,2)*SINB)
     .    +H2*(S(1,1)*P2(1,2)**2-2d0*S(1,3)*P2(1,1)*P2(1,2)*COSB)
     .    +2d0*MUEFF/LAMBDA*(S(1,3)*P2(1,1)**2*SINB*COSB
     .     -S(1,1)*P2(1,1)*P2(1,2)*SINB-S(1,2)*P2(1,1)*P2(1,2)*COSB))
     . +LAMBDA*MUP*dsqrt(2d0)*(-S(1,1)*P2(1,1)*P2(1,2)*SINB
     .      -S(1,2)*P2(1,1)*P2(1,2)*COSB+S(1,3)*P2(1,1)**2*SINB*COSB)
     . +KAPPA*MUP*S(1,3)*P2(1,2)**2*dsqrt(2d0)
     . +6d0*dsqrt(2d0)*DLOG(MAX(QSTSB,MH**2)/(MAX(MT,MH)**2))
     .    *(MT**4/H1**3*S(1,1)*P2(1,1)**2*COSB**2
     .     +MB**4/H2**3*S(1,2)*P2(1,1)**2*SINB**2)/(16d0*PI**2)
        GamHAA=aux**2/(32d0*PI*MH)*dsqrt(1d0-4d0*(PMASS(1)/MH)**2)
       ENDIF

c       Decay to W*W*/Z*Z*

      IF(VFLAG.ne.0)then

        CALL HTOVV(MW,2.08856d0,MH,HTWW)
        GamHWW= 3d0/2d0*GF*MW**4/DSQRT(2d0)/PI/MH**3*HTWW*CV(1)**2

        CALL HTOVV(MZ,2.49581d0,MH,HTZZ)
        GamHZZ= 3d0/4d0*GF*MZ**4/DSQRT(2d0)/PI/MH**3*HTZZ*CV(1)**2

      ELSE

       GamHWW=0d0
       GamHZZ=0d0

      ENDIF

C       Hadronic decays

      GamHhadr=0d0

C   Initializing the strong coupling constant and running masses
      IF(MH.ge.2d0)then
       HIGTOP=(MH/MT)**2
       MT0=3d8
* New May 2019:
       ASH=ALPHAS(MH,3)
       MC0=1d8
       MB0=2d8
       AS3=ALPHAS(MH,3)
       ASH2=AS3*dsqrt(1d0+AS3/Pi*(95d0/4d0-3d0*7d0/6d0))
       MC0=MCC
       AS4=ALPHAS(MH,3)
* End New
       MB0=MBP
       MT0=MT
       RMS=RUNM(MH,3)
      ELSE
       IF(MH.lt.0.8d0)ASH2=Pi/3d0
       IF(MH.ge.0.8d0.and.MH.lt.2d0)ASH2=
     . (Pi/3d0-0.458d0)*dexp(-((MH-0.8d0)/0.732d0)**2)+0.458d0
       IF(MH.gt.1d0)RMS=RUNM(MH,3)
      ENDIF

      IF(MH.lt.4d0)then

      MS2=1d0
      GamS2=0.18d0

      MS3=1.5d0
      GamS3=0.1d0

      MSIG=1d0
      GamSIG=0.7d0

      MG0=1.6d0
      GamG0=0.06d0

C  * H -> 2Pi0
      AH2Pi=CDABS(-DSQRT(dsqrt(2d0)*GF)/2d0*PropII(MH,MS2,GamS2)
     .   *((MPI**2+MKC**2-MK0**2)*CU(1)+(MPI**2+MK0**2-MKC**2)*CD(1))
     .  -2d0*ASH2/PI/3d0*CJH*(MH**2+MPI**2)*PropII(MH,MSIG,GamSIG))**2

      GamH2Pi=AH2Pi/(32d0*PI*MH)*dsqrt(max(0d0,1d0-4d0*(MPI/MH)**2))

C  * H -> Pi+Pi-
      AH2Pi=CDABS(-DSQRT(dsqrt(2d0)*GF)/2d0*PropII(MH,MS2,GamS2)
     .   *((MPI**2+MKC**2-MK0**2)*CU(1)+(MPI**2+MK0**2-MKC**2)*CD(1))
     . -2d0*ASH2/PI/3d0*CJH*(MH**2+MPIC**2)*PropII(MH,MSIG,GamSIG))**2

      GamH2PiC=AH2Pi/(16d0*PI*MH)*dsqrt(max(0d0,1d0-4d0*(MPIC/MH)**2))

C  * H -> EtaPi0
      AHPiE=CDABS(-DSQRT(dsqrt(2d0)*GF/3d0)/2d0*PropII(MH,MS2,GamS2)
     .     *((MPI**2+MKC**2-MK0**2)*CU(1)-(MPI**2+MK0**2-MKC**2)*CD(1))
     .        -2d0*ASH2/PI/dsqrt(3d0)*CJH*PropII(MH,MSIG,GamSIG)
     .                   *(MKC**2-MK0**2-MPIC**2+MPI**2))**2

      GamHPiE=AHPiE/(16d0*PI*MH)
     .       *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))**2
     .       *dsqrt(max(0d0,1d0-(MPI+META)**2/MH**2))
     .       *dsqrt(max(0d0,1d0-(MPI-META)**2/MH**2))

C  * H -> Eta'Pi0
      AHPiE=CDABS(-DSQRT(dsqrt(2d0)*GF/3d0)/2d0*PropII(MH,MS2,GamS2)
     .     *((MPI**2+MKC**2-MK0**2)*CU(1)-(MPI**2+MK0**2-MKC**2)*CD(1))
     .        -2d0*ASH2/PI/dsqrt(3d0)*CJH*PropII(MH,MSIG,GamSIG)
     .                   *(MKC**2-MK0**2-MPIC**2+MPI**2))**2

      GamHPiEP=AHPiE/(16d0*PI*MH)
     .       *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))**2
     .       *dsqrt(max(0d0,1d0-(MPI+METAP)**2/MH**2))
     .       *dsqrt(max(0d0,1d0-(MPI-METAP)**2/MH**2))

C  * H -> K+K-
      AH2K=CDABS(-DSQRT(dsqrt(2d0)*GF)/2d0
     .     *((MPI**2+MKC**2-MK0**2)*CU(1)*PropII(MH,MS2,GamS2)
     .      +(MKC**2+MK0**2-MPI**2)*CD(1)*PropII(MH,MS3,GamS3))
     . -2d0*ASH2/PI/3d0*CJH*(MH**2+MKC**2)*PropII(MH,MSIG,GamSIG))**2

      GamH2KC=AH2K/(16d0*PI*MH)
     .       *dsqrt(max(0d0,1d0-4d0*MKC**2/MH**2))

C  * H -> K0K0b
      AH2K=CDABS(-DSQRT(dsqrt(2d0)*GF)/2d0
     .     *((MPI**2+MK0**2-MKC**2)*CD(1)*PropII(MH,MS2,GamS2)
     .      +(MKC**2+MK0**2-MPI**2)*CD(1)*PropII(MH,MS3,GamS3))
     . -2d0*ASH2/PI/3d0*CJH*(MH**2+MK0**2)*PropII(MH,MSIG,GamSIG))**2

      GamH2K0=AH2K/(16d0*PI*MH)
     .       *dsqrt(max(0d0,1d0-4d0*MK0**2/MH**2))

C  * H -> 2Eta
      AH2E=CDABS(-DSQRT(dsqrt(2d0)*GF)/6d0*
     . (((MPI**2+MKC**2-MK0**2)*CU(1)+(MPI**2+MK0**2-MKC**2)*CD(1))
     . *PropII(MH,MS2,GamS2)*(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))**2
     .   +4d0*(MKC**2+MK0**2-MPI**2)*CD(1)*PropII(MH,MS3,GamS3)
     .              *(DDCOS(THETAE)+DDSIN(THETAE)/dsqrt(2d0))**2)
     .        -2d0*ASH2/PI/3d0*CJH*(PropII(MH,MSIG,GamSIG)
     .     *(MH**2-2d0*META**2+MPI**2*(DDCOS(THETAE)
     .       -dsqrt(2d0)*DDSIN(THETAE))**2+2d0*(MKC**2+MK0**2-MPI**2)
     .          *(DDCOS(THETAE)+DDSIN(THETAE)/dsqrt(2d0))**2)
     .     +PropI(MH,MG0,GamG0)*DDSIN(THETAE)**2))**2

      GamH2E=AH2E/(32d0*PI*MH)*dsqrt(max(0d0,1d0-4d0*(META/MH)**2))

C  * H -> 2Eta'
      AH2E=CDABS(-DSQRT(dsqrt(2d0)*GF)/6d0*
     . (((MPI**2+MKC**2-MK0**2)*CU(1)+(MPI**2+MK0**2-MKC**2)*CD(1))
     . *PropII(MH,MS2,GamS2)*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))**2
     .   +4d0*(MKC**2+MK0**2-MPI**2)*CD(1)*PropII(MH,MS3,GamS3)
     .              *(DDSIN(THETAE)-DDCOS(THETAE)/dsqrt(2d0))**2)
     .        -2d0*ASH2/PI/3d0*CJH*(PropII(MH,MSIG,GamSIG)
     .    *(MH**2-2d0*METAP**2+MPI**2*(DDSIN(THETAE)
     .       +dsqrt(2d0)*DDCOS(THETAE))**2+2d0*(MKC**2+MK0**2-MPI**2)
     .          *(DDSIN(THETAE)-DDCOS(THETAE)/dsqrt(2d0))**2)
     .     +PropI(MH,MG0,GamG0)*DDCOS(THETAE)**2))**2

      GamH2EP=AH2E/(32d0*PI*MH)*dsqrt(max(0d0,1d0-4d0*(METAP/MH)**2))

C  * H -> EtaEta'
      AH2E=CDABS(-DSQRT(dsqrt(2d0)*GF)/6d0*
     . (((MPI**2+MKC**2-MK0**2)*CU(1)+(MPI**2+MK0**2-MKC**2)*CD(1))
     .  *PropII(MH,MS2,GamS2)*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
     .              *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
     .   +4d0*(MKC**2+MK0**2-MPI**2)*CD(1)*PropII(MH,MS3,GamS3)
     .              *(DDSIN(THETAE)-DDCOS(THETAE)/dsqrt(2d0))
     .              *(DDCOS(THETAE)+DDSIN(THETAE)/dsqrt(2d0)))
     .        -2d0*ASH2/PI/3d0*CJH*(PropII(MH,MSIG,GamSIG)*(
     .       MPI**2*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
     .             *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
     .       +2d0*(MKC**2+MK0**2-MPI**2)
     .          *(DDSIN(THETAE)-DDCOS(THETAE)/dsqrt(2d0))
     .          *(DDCOS(THETAE)+DDSIN(THETAE)/dsqrt(2d0)))
     .     +PropI(MH,MG0,GamG0)*DDSIN(THETAE)*DDCOS(THETAE)))**2

      GamHEEP=AH2E/(16d0*PI*MH)
     .                  *dsqrt(max(0d0,1d0-((META+METAP)/MH)**2))
     .                  *dsqrt(max(0d0,1d0-((META-METAP)/MH)**2))

      GamHhadr=GamH2Pi+GamH2PiC+GamHPiE+GamHPiEP+GamH2KC+GamH2K0
     .         +GamH2E+GamH2EP+GamHEEP

C       Diphoton decays

      GamHGAGA=CDABS(ALEM0*CGH/Sqrt(2d0)/Pi
     .         +DSQRT(dsqrt(2d0)*GF)/2d0*PropI(MH,MS2,GamS2)
     .   *((MPI**2+MKC**2-MK0**2)*CU(1)+(MPI**2+MK0**2-MKC**2)*CD(1))
     .   *dsqrt(16d0*Pi*6.7d-7/MS2**3)/2.5d0
     .         +DSQRT(dsqrt(2d0)*GF)/2d0*PropI(MH,MS3,GamS3)
     .   *(MKC**2+MK0**2-MPI**2)*CD(1)*dsqrt(16d0*Pi*4d-7/MS3**3)/2.7d0
     . +2d0*ASH2/PI/3d0*CJH*
     .  (MSIG**2*PropI(MH,MSIG,GamSIG)*dsqrt(16d0*Pi*2d-7/MSIG**3)/5d0
     .   +MG0**2*PropI(MH,MG0,GamG0)*dsqrt(16d0*Pi*1d-6/MG0**3)/2d0)
     .          )**2*MH**3/(32d0*Pi)

      ENDIF

      aux=8d0*CDABS(CGH)**2*MH**3
     .          /(32d0*Pi)*(ALEM0/4d0/Pi)**2
      IF(MH.ge.3d0)then
       IF(MH.lt.4d0)THEN
        GamHGAGA=GamHGAGA*(1d0-((MH-3d0)/1d0)**2)+((MH-3d0)/1d0)**2*aux
       ELSE
        GamHGAGA=aux
       ENDIF
      ENDIF

C       Light-quark and gluon decays

      AHuu=2d-3*CU(1)*DSQRT(dsqrt(2d0)*GF)
      AHdd=4d-3*CD(1)*DSQRT(dsqrt(2d0)*GF)
      AHss=MS*CD(1)*DSQRT(dsqrt(2d0)*GF)

      IF(MH.LT.2d0)THEN
       GamHuu=0d0
       GamHdd=0d0
       GamHss=0d0
       GamHjj=0d0
      ELSE

       RATCOUP=1d0
       GamHuu=3d0*dabs(AHuu)**2/(8d0*Pi)*MH
     .        *(4d0*(2d-3/MH)**2*TQCDH((2d-3/MH)**2)
     .   +(1d0-4d0*(2d-3/MH)**2)*max(0d0,QCDH((2d-3/MH)**2)))
     .               *dsqrt(max(0d0,1d0-4d0*(2d-3/MH)**2))**3

       RATCOUP=0d0
       IF(CD(1).NE.0d0)RATCOUP=CU(1)/CD(1)
       GamHdd=3d0*dabs(AHdd)**2/(8d0*Pi)*MH
     .         *(4d0*(4d-3/MH)**2*TQCDH((4d-3/MH)**2)
     .   +(1d0-4d0*(4d-3/MH)**2)*max(0d0,QCDH((4d-3/MH)**2)))
     .               *dsqrt(max(0d0,1d0-4d0*(4d-3/MH)**2))**3
       GamHss=3d0*dabs(AHss)**2/(8d0*Pi)*MH
     .         *(4d0*(MS/MH)**2*TQCDH((MS/MH)**2)
     .   +(1d0-4d0*(MS/MH)**2)*max(0d0,(RMS/MS)**2*QCDH((RMS/MH)**2)))
     .               *dsqrt(max(0d0,1d0-4d0*(MS/MH)**2))**3
* New May 2019:
       GamHjj=AS3**2/(8d0*PI**3)*MH**3
     .  *Max(0d0,CDABS(CJH)**2*(1d0+AS3/Pi*(95d0/4d0-3d0*7d0/6d0))
     .           +CIH*AS3/PI*7d0/2d0)
* End New
      ENDIF

      aux=GamHhadr
      IF(MH.lt.3d0)then
       GamHhadr=aux
      ELSE
       IF(MH.lt.4d0)THEN
        GamHhadr=aux*(1d0-((MH-3d0)/1d0))
     .          +((MH-3d0)/1d0)*(0d0*GamHuu+0d0*GamHdd+GamHss+GamHjj)
       ELSE
        GamHhadr=0d0*GamHuu+0d0*GamHdd+GamHss+GamHjj
       ENDIF
      ENDIF

C       Mixing with the chi_c(1P)  - hep-ph/9503356

      MCHIC=3.41475d0
      DMC=dsqrt(27d0*dsqrt(2d0)/Pi*GF*MCHIC*0.1d0)*CU(1)    !0.007d0*CU(1)

c      IF(MH.gt.2d0.and.MH.lt.2d0*MTAU)THEN
       CMIX=0d0
       aux=Max(0d0,min(1d0,MH-1d0))
       IF(DMC.ne.0d0)THEN
        CMIX=datan(2d0*DMC
     .        /(MH**2-MCHIC**2+dsqrt((MH**2-MCHIC**2)**2+4d0*DMC**2)))
       ENDIF
       IF(DDCOS(CMIX)**2.gt.DDSIN(CMIX)**2)then
        GamHhadr=GamHhadr*DDCOS(CMIX)**2+10.5d-3*aux*DDSIN(CMIX)**2
        GamHee=GamHee*DDCOS(CMIX)**2
        GamHmumu=GamHmumu*DDCOS(CMIX)**2
        GamHtata=GamHtata*DDCOS(CMIX)**2
        GamHgaga=GamHgaga*DDCOS(CMIX)**2+10.5d-3*2.23d-4*aux
     .                                                 *DDSIN(CMIX)**2
        GamHinv=GamHinv*DDCOS(CMIX)**2
       ELSE
        GamHhadr=GamHhadr*DDSIN(CMIX)**2+10.5d-3*aux*DDCOS(CMIX)**2
        GamHee=GamHee*DDSIN(CMIX)**2
        GamHmumu=GamHmumu*DDSIN(CMIX)**2
        GamHtata=GamHtata*DDSIN(CMIX)**2
        GamHgaga=GamHgaga*DDSIN(CMIX)**2+10.5d-3*2.23d-4*aux
     .                                                 *DDCOS(CMIX)**2
        GamHinv=GamHinv*DDSIN(CMIX)**2
       ENDIF
c      ENDIF

C       H -> cc decays

      MD=1.865d0

      IF(MH.LE.2d0*MD)THEN
       GamHcc= 0d0
      ELSE
       RMC=RUNM(MH,4)
       RATCOUP= 1d0

* New July 2019:
       DCC=MH**3/(8d0*PI**3)*(AS4**2*Max(0d0,CDABS(CJH)**2*
     .     (1d0+AS4/Pi*(95d0/4d0-4d0*7d0/6d0))+CIH*AS4/PI*7d0/2d0)
     .   -AS3**2*Max(0d0,CDABS(CJH)**2*
     .    (1d0+AS3/Pi*(95d0/4d0-3d0*7d0/6d0))+CIH*AS3/PI*7d0/2d0))

       GamHcc=4d0*(MCC/MH)**2*
     .        3d0*GF*(MCC*CU(1))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .        *TQCDH((MCC/MH)**2)*dsqrt(1d0-4d0*MCC**2/MH**2)**3
     .         +(1d0-4d0*(MCC/MH)**2)*max(0d0,
     .         3d0*GF*(RMC*CU(1))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .         *QCDH((RMC/MH)**2)*dsqrt(1d0-4d0*RMC**2/MH**2)**3!+DCC
     .         )
* End new

       IF(MH.le.50d0)then
        GamHcc=GamHcc*((1d0-(MH-2d0*MD)/(50d0-2d0*MD))*
     .  (dsqrt(1d0-4d0*MD**2/MH**2)/dsqrt(1d0-4d0*MCC**2/MH**2))**3
     .  +(MH-2d0*MD)/(50d0-2d0*MD))
* New July 2019:
        GamHhadr=GamHhadr+DCC*
     .  ((1d0-4d0*(MCC/MH)**2)*(1d0-(MH-2d0*MD)/(50d0-2d0*MD))*
     .  (dsqrt(1d0-4d0*MD**2/MH**2)/dsqrt(1d0-4d0*MCC**2/MH**2))**3
     .  +(MH-2d0*MD)/(50d0-2d0*MD))
* End new
       ENDIF
      ENDIF

C       Mixing with the chi_b(1,2P)   - Drees, Hikasa (1990)

      MCHIB1=9.85944d0
      MCHIB2=10.2325d0

      DMB1=dsqrt(27d0*dsqrt(2d0)/Pi*GF*MCHIB1*1.7d0)*CB(1) ! 0.049d0*CB(1)
      DMB2=dsqrt(27d0*dsqrt(2d0)/Pi*GF*MCHIB2*2.0d0)*CB(1) ! 0.054d0*CB(1)

c      IF(MH.gt.4d0.and.MH.lt.2d0*5.2795d0)THEN

       MMIX(1,1)=MCHIB1**2
       MMIX(1,2)=0d0
       MMIX(1,3)=DMB1
       MMIX(2,1)=0d0
       MMIX(2,2)=MCHIB2**2
       MMIX(2,3)=DMB2
       MMIX(3,1)=DMB1
       MMIX(3,2)=DMB2
       MMIX(3,3)=MH**2

       CALL DIAGN(3,MMIX,VALP,VECP,EPS)
       CALL SORTN(3,VALP,VECP)
       DO K=1,3
        DO L=1,3
         OMIX(K,L)=VECP(L,K)
        ENDDO
       ENDDO

       INDH=1
       DO K=1,3
        IF(dabs(OMIX(INDH,3))**2.lt.dabs(OMIX(K,3))**2)INDH=K
       ENDDO
c chi_b widths from 1601.05093
       aux=Max(0d0,min(1d0,MH/2d0-1.5d0))
       GamHhadr=GamHhadr*OMIX(INDH,3)**2+aux*(2.03d-3*OMIX(INDH,1)**2
     .                                       +2.39d-3*OMIX(INDH,2)**2)
       GamHee=GamHee*OMIX(INDH,3)**2
       GamHmumu=GamHmumu*OMIX(INDH,3)**2
       GamHtata=GamHtata*OMIX(INDH,3)**2
       GamHgaga=GamHgaga*OMIX(INDH,3)**2+aux*(0.12d-3*OMIX(INDH,1)**2
     .                                       +0.14d-3*OMIX(INDH,2)**2)
       GamHcc=GamHcc*OMIX(INDH,3)**2
       GamHinv=GamHinv*OMIX(INDH,3)**2
! GamAhadr also contains contributions to the cc final state
c      ENDIF

C       H -> bb decays

      MBm=5.2795d0

      IF(MH.LE.2d0*MBm)THEN
       GamHbb= 0d0
      ELSE
       RMB=RUNMB(MH)
       IF(CB(1).ne.0d0)THEN
        RATCOUP= CU(1)/CB(1)
       ELSE
        RATCOUP=0d0
       ENDIF

* New July 2019:
       DBB=MH**3/(8d0*PI**3)*(ASH**2*Max(0d0,CDABS(CJH)**2*
     .    (1d0+ASH/Pi*(95d0/4d0-5d0*7d0/6d0))+CIH*ASH/PI*7d0/2d0)
     .   -AS4**2*Max(0d0,CDABS(CJH)**2*
     .    (1d0+AS4/Pi*(95d0/4d0-4d0*7d0/6d0))+CIH*AS4/PI*7d0/2d0))

       DBB=DBB+MH**3/(8d0*PI**3)*CDABS(CJHF)**2*ASH**2*(
     .      HGGQCD2(ASH,5,MH,MT)
     .   -(1d0+ASH/Pi*(95d0/4d0-5d0*7d0/6d0)))

       GamHbb=4d0*(MBP/MH)**2*
     .        3d0*GF*(MBP*CB(1))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .        *TQCDH((MBP/MH)**2)*dsqrt(1d0-4d0*MBP**2/MH**2)**3
     .         +(1d0-4d0*(MBP/MH)**2)*max(0d0,
     .         3d0*GF*(RMB*CB(1))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .      *QCDH((RMB/MH)**2)*dsqrt(1d0-4d0*RMB**2/MH**2)**3!+DBB
     .       )
* End new
       IF(MH.le.50d0)THEN
        GamHbb=GamHbb*((1d0-(MH-2d0*MBm)/(50d0-2d0*MBm))*
     .  (dsqrt(1d0-4d0*MBm**2/MH**2)/dsqrt(1d0-4d0*MBP**2/MH**2))**3
     .  +(MH-2d0*MBm)/(50d0-2d0*MBm))
* New July 2019:
        GamHhadr=GamHhadr+DBB*
     .  ((1d0-4d0*(MBP/MH)**2)*(1d0-(MH-2d0*MBm)/(50d0-2d0*MBm))*
     .  (dsqrt(1d0-4d0*MBm**2/MH**2)/dsqrt(1d0-4d0*MBP**2/MH**2))**3
     .  +(MH-2d0*MBm)/(50d0-2d0*MBm))
* End new
       ENDIF
      ENDIF

      RETURN
      END

*********************************************************************

      DOUBLE COMPLEX FUNCTION PropI(MH,M,G)

      DOUBLE PRECISION MH,M,G

       PropI=M**2/DCMPLX(MH**2-M**2,M*G)  

      END

*********************************************************************

      DOUBLE COMPLEX FUNCTION PropII(MH,M,G)

      DOUBLE PRECISION MH,M,G

      IF(MH.le.M)THEN
       PropII=M**2/DCMPLX(MH**2-M**2,M*G)    
      ELSE
       PropII=((MH**2-M**2)**3+M**6)**(1d0/3d0)/DCMPLX(MH**2-M**2,M*G)
      ENDIF 

      END

