      SUBROUTINE HDECAY_lightCPV(IND)

      IMPLICIT NONE

      INTEGER IND,INDA,K,L,N0,VFLAG,MIXPROB,NF

      DOUBLE PRECISION MH,EPS,PI,aux
      DOUBLE PRECISION BETA,X,SP,QCD0,HQCDM,HQCD,QCDH,TQCDH,RATCOUP
      DOUBLE PRECISION AQCDM,AQCD,QCDA,TQCDA
      DOUBLE PRECISION FPI,THETAE,MPI,MPIC,MPI8,MPI9,META,METAP,MKc,MK0
      DOUBLE PRECISION MS2,GamS2,MS3,GamS3,MSIG,GamSIG
      DOUBLE PRECISION DM3,DM8,DM9,MMIX(4,4),VALP(4),VECP(4,4),
     .                 OMIX(4,4)
      DOUBLE PRECISION AHee,AHmumu,AHtata,MG0,GamG0
      DOUBLE PRECISION RMS,ASH,ASH2,AS3,AS4,HIGTOP,ALPHAS,RUNM
      DOUBLE PRECISION AH2Pi,AHPiE,GamH2Pi,GamH2PiC,GamHPiE,GamHPiEP
      DOUBLE PRECISION AH2K,GamH2KC,GamH2K0,AH2E,GamH2E,GamHEEP,GamH2EP
      DOUBLE PRECISION AHuu,AHdd,AHss,GamHuu,GamHdd,GamHss,GamHjj
      DOUBLE PRECISION MD,HCC,ACC,RMC,MBm,HBB,ABB,RMB,runmb,HTWW,HTZZ
      DOUBLE PRECISION GamAGAGA,GamAee,GamAmumu,GamAtata,GamAhadr,
     . GamAcc,GamAbb
      DOUBLE PRECISION CGPI,CGETA,CGETAP
      DOUBLE PRECISION APIee,AEee,AEPee
      DOUBLE PRECISION APImumu,AEmumu,AEPmumu
      DOUBLE PRECISION YUE,YDE,YSE
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
      DOUBLE PRECISION METAC,DMC,KinAPiDD,MCHIC,DMCP,MMIXC(3,3),
     .                 VALPC(3),VECPC(3,3),OMIXC(3,3)
      DOUBLE PRECISION METAB1,METAB2,METAB3,DMB1,DMB2,DMB3,KinAPiBB
      DOUBLE PRECISION MMIX2(6,6),VALP2(6),VECP2(6,6),OMIX2(6,6),
     .                 DMB1P,DMB2P,MCHIB1,MCHIB2
* New May 2019:
      DOUBLE PRECISION ZETA2,ZETA3,ASG,HGGQCD2,AGGQCD2,CIH

      DOUBLE COMPLEX XXC,TC,BC,CC,LC,MC,EC,CH1C,CH2C,WC,HC,PropI,PropII
      DOUBLE COMPLEX ULC,URC,DLC,DRC,T1C,T2C,B1C,B2C,LLC,LRC,L1C,L2C
      DOUBLE COMPLEX CJH,CGH,CJA,CGA,F0,FF,FS,FV,CM,FGG,CJHF

      DOUBLE PRECISION ALEM0
      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION MCH2(2),U(2,2,2),V(2,2,2)
      DOUBLE PRECISION MNEU(5),NEU(5,5,2)
      DOUBLE PRECISION CU(5),CUP(5),CD(5),CDP(5),CB(5),CBP(5),CJ(5),
     . CJP(5),CI(5),CG(5),CGP(5),CV(5),CZG(5),CZGP(5),CL(5),CLP(5)
      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION mt,mb,mtau,mmuon,mel,MS,MCC,MBP,MPI0,MSTRANGE
      DOUBLE PRECISION XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION MST2(2),UT(2,2,2),MSB2(2),UB(2,2,2),MSL2(2),
     . UTAU(2,2,2),MSNT2
      DOUBLE PRECISION MSU2(2),MSD2(2),MSE2(2),MSNE2,MSMU2(2),
     . UMU(2,2,2)
      DOUBLE PRECISION MST2P(2),MSB2P(2),MSU2P(2),MSD2P(2)
      DOUBLE PRECISION GRHSTST(5,2,2),GRHSBSB(5,2,2),GRHSLSL(5,2,2),
     . GRHSUSU(5,2,2),GRHSDSD(5,2,2),GRHSESE(5,2,2),GRHSNSN(5)
      DOUBLE PRECISION GIHSTST(5,2,2),GIHSBSB(5,2,2),GIHSLSL(5,2,2)
      DOUBLE PRECISION GRHCSTSB(2,2),GRHCSNSL(2),GRHCSUSD(2,2),
     . GRHCSNSE(2),GIHCSTSB(2,2),GIHCSNSL(2)
      DOUBLE PRECISION COH0CH(5,2,2,2),COH0NEU(5,5,5,2),
     . COHPNEUCHM(2,5,2,2),COHMNEUCHP(2,5,2,2)
      DOUBLE PRECISION gH0H0H0(5,5,5),gRH0HPHM(5,2,2),gIH0HPHM(5,2,2)
      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION GamHGAGA,GamHee,GamHmumu,GamHtata,GamHhadr,
     . GamHcc,GamHbb,GamHinv(3),GamHWW,GamHZZ,GamHAA
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/ALEM0/ALEM0
      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/CHASPEC/MCH2,U,V
      COMMON/NEUSPEC/MNEU,NEU
      COMMON/HNSMCOUP/CU,CUP,CD,CDP,CB,CBP,CJ,CJP,CI,CG,CGP,CV,CZG,CZGP,
     . CL,CLP
      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/SMFERM/mt,mb,mtau,mmuon,mel,MS,MCC,MBP,MPI0,MSTRANGE
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/SFERM3SPEC/MST2,UT,MSB2,UB,MSL2,UTAU,MSNT2
      COMMON/SFERM1SPEC/MSU2,MSD2,MSE2,MSNE2,MSMU2,UMU
      COMMON/SFERMPSPEC/MST2P,MSB2P,MSU2P,MSD2P
      COMMON/HISFCOUP/GRHSTST,GRHSBSB,GRHSLSL,GRHSUSU,GRHSDSD,
     . GRHSESE,GRHSNSN,GIHSTST,GIHSBSB,GIHSLSL,GRHCSTSB,GRHCSNSL,
     . GRHCSUSD,GRHCSNSE,GIHCSTSB,GIHCSNSL
      COMMON/HINOCOUP/COH0CH,COH0NEU,COHPNEUCHM,COHMNEUCHP
      COMMON/HICOUP/gH0H0H0,gRH0HPHM,gIH0HPHM
      COMMON/STSBSCALE/QSTSB
      COMMON/VFLAG/VFLAG
      COMMON/LIGHTHDECAYS/GamHGAGA,GamHee,GamHmumu,GamHtata,GamHhadr,
     . GamHcc,GamHbb,GamHinv,GamHWW,GamHZZ,GamHAA

      CM(X)= DCMPLX(MIN(1d3,X)**2,-EPS/4d0)
      F0(XXC)=-CDLOG((CDSQRT(1d0-4d0*XXC)-1d0)
     .               /(CDSQRT(1d0-4d0*XXC)+1d0))**2/4d0
      FF(XXC)=8d0*XXC*(1d0+(1d0-4d0*XXC)*F0(XXC))
      FS(XXC)=2d0*XXC*(4d0*XXC*F0(XXC)-1d0)
      FV(XXC)=-(2d0+12d0*XXC+24d0*XXC*(1d0-2d0*XXC)*F0(XXC))
      FGG(XXC)=2d0*XXC
     .        *CDLOG((CDSQRT(1d0-4d0*XXC)-1d0)
     .               /(CDSQRT(1d0-4d0*XXC)+1d0))**2
      BETA(X)= DSQRT(1d0-4d0*X)
      QCd0(X)= (1d0+X**2)*(4d0*SP((1d0-X)/(1d0+X))
     . +2d0*SP((X-1d0)/(X+1d0))
     . - 3d0*DLOG((1d0+X)/(1d0-X))*DLOG(2d0/(1d0+X))
     . - 2d0*DLOG((1d0+X)/(1d0-X))*DLOG(X))
     . - 3d0*X*DLOG(4d0/(1d0-X**2))-4d0*X*DLOG(X)
      HQCDM(X)= QCd0(X)/X+(3d0+34d0*X**2-13d0*X**4)/16d0/X**3
     . * DLOG((1d0+X)/(1d0-X))+3d0/8d0/X**2*(7d0*X**2-1d0)
      AQCDM(X)= QCd0(X)/X+(19d0+2d0*X**2+3d0*X**4)/16d0/X
     . * DLOG((1d0+X)/(1d0-X))+3d0/8d0*(7d0-X**2)
* July 2010:
c      HQCD(X)=(4d0/3d0*HQCDM(BETA(X))
c     .   +2d0*(4d0/3d0-DLOG(X))*(1d0-10d0*X)/(1d0-4d0*X))*ASH/PI
c     .   + (29.14671d0 + RATCOUP*(1.570d0 - 2d0*DLOG(HIGTOP)/3d0
c     .   + DLOG(X)**2/9d0))*(ASH/PI)**2
c     .   +(164.14d0 - 25.77d0*5 + 0.259d0*5**2)*(ASH/PI)**3
c      AQCD(X)=(4d0/3d0*AQCDM(BETA(X))
c     .   +2d0*(4d0/3d0-DLOG(X))*(1d0-6d0*X)/(1d0-4d0*X))*ASH/PI
c     .   + (29.14671d0 + RATCOUP*(23d0/6d0 - DLOG(HIGTOP)
c     .   + DLOG(X)**2/6d0))*(ASH/PI)**2
c     .   + (164.14d0 - 25.77d0*5 + 0.259d0*5**2)*(ASH/PI)**3
* New May 2019:
      HQCD(X)=(4d0/3d0*HQCDM(BETA(X))
     .       + 2d0*(4d0/3d0-DLOG(X))*(1d0-10d0*X)/(1d0-4d0*X))*ASH/PI
     .       + (29.14671d0
     .         + RATCOUP*(1.570d0 - 2d0*DLOG(HIGTOP)/3d0
     .            + DLOG(X)**2/9d0))*(ASH/PI)**2
     .       + (164.14d0 - 25.77d0*5d0 + 0.259d0*5d0**2)*(ASH/PI)**3
     .       + (39.34d0-220.9d0*5d0+9.685d0*5d0**2
     .         - 0.0205d0*5d0**3)*(ASH/PI)**4
      AQCD(X)=(4d0/3d0*AQCDM(BETA(X))
     .       + 2d0*(4d0/3d0-DLOG(X))*(1d0-6d0*X)/(1d0-4d0*X))*ASH/PI
     .       + (29.14671d0 + RATCOUP*(23d0/6d0 - DLOG(HIGTOP)
     .         + DLOG(X)**2/6d0))*(ASH/PI)**2
     .       + (164.14d0 - 25.77d0*5d0 + 0.259d0*5d0**2)*(ASH/PI)**3
     .       + (39.34d0-220.9d0*5d0+9.685d0*5d0**2
     .         - 0.0205d0*5d0**3)*(ASH/PI)**4
* End New
      QCDH(X)= 1d0+HQCD(X)
      TQCDH(X)= 1d0+4d0/3d0*HQCDM(BETA(X))*ASH/PI
      QCDA(X)= 1d0+AQCD(X)
      TQCDA(X)= 1d0+4d0/3d0*AQCDM(BETA(X))*ASH/PI
* New May 2019:
      HGGQCD2(ASG,NF,MH,MT)= 1d0+ASG/PI*(95d0/4d0-NF*7d0/6d0)
     . +(ASG/PI)**2*(149533d0/288d0-363d0/8d0*ZETA2-495d0/8d0*ZETA3
     .              +19d0/8d0*DLOG(MH**2/MT**2)
     . +NF*(-4157d0/72d0+11d0/2d0*ZETA2+5d0/4d0*ZETA3
     . +2d0/3d0*DLOG(MH**2/MT**2))
     . +NF**2*(127d0/108d0-1d0/6d0*ZETA2))+(ASG/PI)**3
     . *(467.683620788d0+122.440972222d0*DLOG(MH**2/MT**2)
     .              +10.9409722222d0*DLOG(MH**2/MT**2)**2)
      AGGQCD2(ASG,NF,MH,MT)=1d0+ASG/PI*(97d0/4d0-NF*7d0/6d0)
     . +(ASG/PI)**2*(237311d0/864d0-529d0/24d0*ZETA2-445d0/8d0*ZETA3
     . +5d0*DLOG(MH**2/MT**2))
* End New

      EPS=1d-8
      PI=4d0*DATAN(1d0)
* New May 2019:
      ZETA2=PI**2/6d0
      ZETA3=1.202056903159594d0
* End New

      MH=dsqrt(MH0(IND))

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
      CH1C=CM(dsqrt(MCH2(1))/MH)
      CH2C=CM(dsqrt(MCH2(2))/MH)
      WC=CM(MW/MH)
      HC=CM(dsqrt(MHC)/MH)
      ULC=CM(dsqrt(MSU2(1))/MH)
      URC=CM(dsqrt(MSU2(2))/MH)
      DLC=CM(dsqrt(MSD2(1))/MH)
      DRC=CM(dsqrt(MSD2(2))/MH)
      T1C=CM(dsqrt(MST2(1))/MH)
      T2C=CM(dsqrt(MST2(2))/MH)
      B1C=CM(dsqrt(MSB2(1))/MH)
      B2C=CM(dsqrt(MSB2(2))/MH)
      LLC=CM(dsqrt(MSE2(1))/MH)
      LRC=CM(dsqrt(MSE2(2))/MH)
      L1C=CM(dsqrt(MSL2(1))/MH)
      L2C=CM(dsqrt(MSL2(2))/MH)

      CJHF=DSQRT(dsqrt(2d0)*GF)/4d0*(CU(IND)*(FF(TC)+FF(CC))
     .                                                 +CB(IND)*FF(BC))

      CJH=CJHF
     .  +(gRHSUSU(IND,1,1)*FS(ULC)/MSU2(1)
     .   +gRHSUSU(IND,2,2)*FS(URC)/MSU2(2)
     .   +gRHSDSD(IND,1,1)*FS(DLC)/MSD2(1)
     .   +gRHSDSD(IND,2,2)*FS(DRC)/MSD2(2))/2d0
     .  +DCMPLX(GRHSTST(IND,1,1),GIHSTST(IND,1,1))*FS(T1C)/4d0/MST2(1)
     .  +DCMPLX(GRHSTST(IND,2,2),GIHSTST(IND,2,2))*FS(T2C)/4d0/MST2(2)
     .  +DCMPLX(GRHSBSB(IND,1,1),GIHSBSB(IND,1,1))*FS(B1C)/4d0/MSB2(1)
     .  +DCMPLX(GRHSBSB(IND,2,2),GIHSBSB(IND,2,2))*FS(B2C)/4d0/MSB2(2)

* New July 2019:
      CIH=DREAL(DCONJG(CJH)*(CJH-CJHF))
* End New

      CGH=DSQRT(dsqrt(2d0)*GF)/2d0*(4d0/3d0*CU(IND)*(FF(TC)+FF(CC))
     .             +CB(IND)*FF(BC)/3d0+CD(IND)*(FF(MC)+FF(EC))
     .             +CL(IND)*FF(LC)+CV(IND)*FV(WC))
     .  +GRH0HPHM(IND,2,2)*FS(HC)/2d0/MHC
     .  +COH0CH(IND,1,1,1)/dsqrt(MCH2(1))*FF(CH1C)/2d0
     .  +COH0CH(IND,2,2,1)/dsqrt(MCH2(2))*FF(CH2C)/2d0
     .  +2d0/3d0*(
     .    DCMPLX(GRHSTST(IND,1,1),GIHSTST(IND,1,1))*FS(T1C)/MST2(1)
     .   +DCMPLX(GRHSTST(IND,2,2),GIHSTST(IND,2,2))*FS(T2C)/MST2(2)
     .   +2d0*gRHSUSU(IND,1,1)*FS(ULC)/MSU2(1)
     .   +2d0*gRHSUSU(IND,2,2)*FS(URC)/MSU2(2))
     .  +1d0/6d0*(
     .    DCMPLX(GRHSBSB(IND,1,1),GIHSBSB(IND,1,1))*FS(B1C)/MSB2(1)
     .   +DCMPLX(GRHSBSB(IND,2,2),GIHSBSB(IND,2,2))*FS(B2C)/MSB2(2)
     .   +2d0*gRHSDSD(IND,1,1)*FS(DLC)/MSD2(1)
     .   +2d0*gRHSDSD(IND,2,2)*FS(DRC)/MSD2(2))
     .  +GRHSESE(IND,1,1)*FS(LLC)/MSE2(1)
     .  +GRHSESE(IND,1,1)*FS(LRC)/MSE2(2)
     .  +DCMPLX(GRHSLSL(IND,1,1),GIHSLSL(IND,1,1))*FS(L1C)/MSL2(1)/2d0
     .  +DCMPLX(GRHSLSL(IND,2,2),GIHSLSL(IND,2,2))*FS(L2C)/MSL2(2)/2d0

      CJA=-DSQRT(dsqrt(2d0)*GF)/4d0
     .         *(CUP(IND)*(FGG(TC)+FGG(CC))+CBP(IND)*FGG(BC))

      CGA=-DSQRT(dsqrt(2d0)*GF)/2d0*(4d0/3d0*CUP(IND)*(FGG(TC)+FGG(CC))
     .        +CBP(IND)/3d0*FGG(BC)+CDP(IND)*(FGG(MC)+FGG(EC))
     .        +CLP(IND)*FGG(LC))
     . -COH0CH(IND,1,1,2)/dsqrt(MCH2(1))*FGG(CH1C)/2d0
     . -COH0CH(IND,2,2,2)/dsqrt(MCH2(2))*FGG(CH2C)/2d0

C       Leptonic decays

C  * H -> ee
      AHee=MEL*CD(IND)*DSQRT(dsqrt(2d0)*GF)
      GamHee=AHee**2/(8d0*Pi)*MH
     .               *dsqrt(max(0d0,1d0-4d0*MEL**2/MH**2))**3

      AHee=MEL*CDP(IND)*DSQRT(dsqrt(2d0)*GF)
      GamAee=AHee**2/(8d0*Pi)*MH*dsqrt(max(0d0,1d0-4d0*MEL**2/MH**2))

C  * H -> mumu
      AHmumu=MMUON*CD(IND)*DSQRT(dsqrt(2d0)*GF)
      GamHmumu=AHmumu**2/(8d0*Pi)*MH
     .               *dsqrt(max(0d0,1d0-4d0*MMUON**2/MH**2))**3

      AHmumu=MMUON*CDP(IND)*DSQRT(dsqrt(2d0)*GF)
      GamAmumu=AHmumu**2/(8d0*Pi)*MH
     .               *dsqrt(max(0d0,1d0-4d0*MMUON**2/MH**2))

C  * H -> tautau
      AHtata=MTAU*CL(IND)*DSQRT(dsqrt(2d0)*GF)
      GamHtata=AHtata**2/(8d0*Pi)*MH
     .               *dsqrt(max(0d0,1d0-4d0*MTAU**2/MH**2))**3

      AHtata=MTAU*CLP(IND)*DSQRT(dsqrt(2d0)*GF)
      GamAtata=AHtata**2/(8d0*Pi)*MH
     .               *dsqrt(max(0d0,1d0-4d0*MTAU**2/MH**2))

C  * A -> 2gamma
      GamAGAGA=8d0*CDABS(CGA)**2*MH**3
     .          /(32d0*Pi)*(ALEM0/4d0/Pi)**2

c       Decay to light neutralinos

      IF(MH.LE.2d0*DSQRT(MNEU(1)))THEN
       GamHinv(1)=0d0
      ELSE
       aux=COH0NEU(IND,1,1,1)**2
       GamHinv(1)=aux**2/(16d0*PI)*MH*dsqrt(1d0-4d0*MNEU(1)/MH**2)**3
       aux=COH0NEU(IND,1,1,2)**2
       GamHinv(1)=GamHinv(1)+aux**2/(16d0*PI)*MH
     .                               *dsqrt(1d0-4d0*MNEU(1)/MH**2)
      ENDIF

      IF(MH.LE.DSQRT(MNEU(1))+DSQRT(MNEU(2)))THEN
       GamHinv(2)=0d0
      ELSE
       aux=COH0NEU(IND,1,2,1)**2
       GamHinv(2)=aux**2/(8d0*PI)*MH
     .   *dsqrt(1d0-((dsqrt(MNEU(1))+dsqrt(MNEU(2)))/MH)**2)**3
     .   *dsqrt(1d0-((dsqrt(MNEU(1))-dsqrt(MNEU(2)))/MH)**2)

       aux=COH0NEU(IND,1,2,2)**2
       GamHinv(2)=GamHinv(2)+aux**2/(8d0*PI)*MH
     .   *dsqrt(1d0-((dsqrt(MNEU(1))+dsqrt(MNEU(2)))/MH)**2)
     .   *dsqrt(1d0-((dsqrt(MNEU(1))-dsqrt(MNEU(2)))/MH)**2)**3
      ENDIF

      IF(MH.LE.2d0*DSQRT(MNEU(2)))THEN
       GamHinv(3)=0d0
      ELSE
       aux=COH0NEU(IND,2,2,1)**2
       GamHinv(3)=aux**2/(16d0*PI)*MH*dsqrt(1d0-4d0*MNEU(2)/MH**2)**3

       aux=COH0NEU(IND,2,2,2)**2
       GamHinv(3)=GamHinv(3)+aux**2/(16d0*PI)*MH
     .                       *dsqrt(1d0-4d0*MNEU(2)/MH**2)
      ENDIF

c       Decay to light Higgs states

       IF(MH.LE.2d0*dsqrt(MH0(1)))THEN
        GamHAA=0d0
       ELSE
        aux=GH0H0H0(IND,1,1)

        GamHAA=aux**2/(32d0*PI*MH)*dsqrt(1d0-4d0*MH0(1)/MH**2)
       ENDIF

c       Decay to W*W*/Z*Z*

      IF(VFLAG.ne.0)then

        CALL HTOVV(MW,2.08856d0,MH,HTWW)
        GamHWW= 3d0/2d0*GF*MW**4/DSQRT(2d0)/PI/MH**3*HTWW*CV(IND)**2

        CALL HTOVV(MZ,2.49581d0,MH,HTZZ)
        GamHZZ= 3d0/4d0*GF*MZ**4/DSQRT(2d0)/PI/MH**3*HTZZ*CV(IND)**2

      ELSE

       GamHWW=0d0
       GamHZZ=0d0

      ENDIF

C       Initializing the strong coupling constant and running masses
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

C       Mixing with mesons [1612.06538[hep-ph]]

      GamHhadr=0d0
      GamAhadr=0d0

      IF(MH.lt.4d0)then

      DM3=-DSQRT(dsqrt(2d0)*GF)/2d0*FPI*MPI**2*(CUP(IND)-CDP(IND))
      DM8=-DSQRT(dsqrt(2d0)*GF)/2d0*FPI*(-dsqrt(3d0)*MPI8**2*CDP(IND)
     .                      +MPI**2/dsqrt(3d0)*(CUP(IND)+2d0*CDP(IND)))
      DM9=-DSQRT(dsqrt(2d0)*GF)/2d0*FPI*(dsqrt(3d0/2d0)*MPI8**2
     .             *CDP(IND)+MPI**2/dsqrt(6d0)*(2d0*CUP(IND)+CDP(IND)))
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
      MMIX(4,4)=MH**2
     .    +2d0/3d0*FPI**2*(MPI9**2-(MPI8**2+MPI**2)/3d0)*CDABS(CJA)**2

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
      IF(dabs(OMIX(INDA,4))**2.lt.0.6d0)MIXPROB=1

C       Diphoton decay A -> 2gamma - CP-odd modes

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
      GamAGAGA=8d0*(aux**2+(OMIX(INDA,4)*DIMAG(CGA))**2)*MH**3
     .          /(32d0*Pi)*(ALEM0/4d0/Pi)**2

C       Leptonic decays - CP-odd modes

C  * A -> ee
      AHee=MEL*CDP(IND)*DSQRT(dsqrt(2d0)*GF)
      APiee=dsqrt(6.582d-25/8.52d-17*6.46d-8
     .            *8d0*Pi/(MPI*dsqrt(1d0-4d0*MEL**2/MPI**2)))
      AEee=0d0
      AEPee=0d0

      aux=OMIX(INDA,1)*APiee+OMIX(INDA,2)*AEee+OMIX(INDA,3)*AEPee
     .          +OMIX(INDA,4)*AHee
      GamAee=aux**2/(8d0*Pi)*MH*dsqrt(max(0d0,1d0-4d0*MEL**2/MH**2))

C  * A -> mumu
      AHmumu=MMUON*CDP(IND)*DSQRT(dsqrt(2d0)*GF)
      APimumu=0d0
      AEmumu=dsqrt(5.8d-6*1.31d-6
     .            *8d0*Pi/(META*dsqrt(1d0-4d0*MMUON**2/META**2)))
      AEPmumu=0d0

      aux=OMIX(INDA,1)*APimumu+OMIX(INDA,2)*AEmumu+OMIX(INDA,3)*AEPmumu
     .          +OMIX(INDA,4)*AHmumu
      GamAmumu=aux**2/(8d0*Pi)*MH
     .               *dsqrt(max(0d0,1d0-4d0*MMUON**2/MH**2))

C  * A -> tautau
      AHtata=MTAU*CLP(IND)*DSQRT(dsqrt(2d0)*GF)
      GamAtata=(OMIX(INDA,4)*AHtata)**2/(8d0*Pi)*MH
     .               *dsqrt(max(0d0,1d0-4d0*MTAU**2/MH**2))

C       Hadronic decays

c              - 2 meson channels (CP-even modes)

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
     . *((MPI**2+MKC**2-MK0**2)*CU(IND)+(MPI**2+MK0**2-MKC**2)*CD(IND))
     .  -2d0*ASH2/PI/3d0*CJH*(MH**2+MPI**2)*PropII(MH,MSIG,GamSIG))**2

      GamH2Pi=AH2Pi/(32d0*PI*MH)*dsqrt(max(0d0,1d0-4d0*(MPI/MH)**2))

C  * H -> Pi+Pi-
      AH2Pi=CDABS(-DSQRT(dsqrt(2d0)*GF)/2d0*PropII(MH,MS2,GamS2)
     . *((MPI**2+MKC**2-MK0**2)*CU(IND)+(MPI**2+MK0**2-MKC**2)*CD(IND))
     . -2d0*ASH2/PI/3d0*CJH*(MH**2+MPIC**2)*PropII(MH,MSIG,GamSIG))**2

      GamH2PiC=AH2Pi/(16d0*PI*MH)*dsqrt(max(0d0,1d0-4d0*(MPIC/MH)**2))

C  * H -> EtaPi0
      AHPiE=CDABS(-DSQRT(dsqrt(2d0)*GF/3d0)/2d0*PropII(MH,MS2,GamS2)
     . *((MPI**2+MKC**2-MK0**2)*CU(IND)-(MPI**2+MK0**2-MKC**2)*CD(IND))
     .        -2d0*ASH2/PI/dsqrt(3d0)*CJH*PropII(MH,MSIG,GamSIG)
     .                   *(MKC**2-MK0**2-MPIC**2+MPI**2))**2

      GamHPiE=AHPiE/(16d0*PI*MH)
     .       *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))**2
     .       *dsqrt(max(0d0,1d0-(MPI+META)**2/MH**2))
     .       *dsqrt(max(0d0,1d0-(MPI-META)**2/MH**2))

C  * H -> Eta'Pi0
      AHPiE=CDABS(-DSQRT(dsqrt(2d0)*GF/3d0)/2d0*PropII(MH,MS2,GamS2)
     . *((MPI**2+MKC**2-MK0**2)*CU(IND)-(MPI**2+MK0**2-MKC**2)*CD(IND))
     .        -2d0*ASH2/PI/dsqrt(3d0)*CJH*PropII(MH,MSIG,GamSIG)
     .                   *(MKC**2-MK0**2-MPIC**2+MPI**2))**2

      GamHPiEP=AHPiE/(16d0*PI*MH)
     .       *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))**2
     .       *dsqrt(max(0d0,1d0-(MPI+METAP)**2/MH**2))
     .       *dsqrt(max(0d0,1d0-(MPI-METAP)**2/MH**2))

C  * H -> K+K-
      AH2K=CDABS(-DSQRT(dsqrt(2d0)*GF)/2d0
     .     *((MPI**2+MKC**2-MK0**2)*CU(IND)*PropII(MH,MS2,GamS2)
     .      +(MKC**2+MK0**2-MPI**2)*CD(IND)*PropII(MH,MS3,GamS3))
     . -2d0*ASH2/PI/3d0*CJH*(MH**2+MKC**2)*PropII(MH,MSIG,GamSIG))**2

      GamH2KC=AH2K/(16d0*PI*MH)
     .       *dsqrt(max(0d0,1d0-4d0*MKC**2/MH**2))

C  * H -> K0K0b
      AH2K=CDABS(-DSQRT(dsqrt(2d0)*GF)/2d0
     .     *((MPI**2+MK0**2-MKC**2)*CD(IND)*PropII(MH,MS2,GamS2)
     .      +(MKC**2+MK0**2-MPI**2)*CD(IND)*PropII(MH,MS3,GamS3))
     . -2d0*ASH2/PI/3d0*CJH*(MH**2+MK0**2)*PropII(MH,MSIG,GamSIG))**2

      GamH2K0=AH2K/(16d0*PI*MH)
     .       *dsqrt(max(0d0,1d0-4d0*MK0**2/MH**2))

C  * H -> 2Eta
      AH2E=CDABS(-DSQRT(dsqrt(2d0)*GF)/6d0*
     . (((MPI**2+MKC**2-MK0**2)*CU(IND)+(MPI**2+MK0**2-MKC**2)*CD(IND))
     . *PropII(MH,MS2,GamS2)*(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))**2
     .   +4d0*(MKC**2+MK0**2-MPI**2)*CD(IND)*PropII(MH,MS3,GamS3)
     .              *(DDCOS(THETAE)+DDSIN(THETAE)/dsqrt(2d0))**2)
     .        -2d0*ASH2/PI/3d0*CJH*(PropII(MH,MSIG,GamSIG)
     .     *(MH**2-2d0*META**2+MPI**2*(DDCOS(THETAE)
     .       -dsqrt(2d0)*DDSIN(THETAE))**2+2d0*(MKC**2+MK0**2-MPI**2)
     .          *(DDCOS(THETAE)+DDSIN(THETAE)/dsqrt(2d0))**2)
     .     +PropI(MH,MG0,GamG0)*DDSIN(THETAE)**2))**2

      GamH2E=AH2E/(32d0*PI*MH)*dsqrt(max(0d0,1d0-4d0*(META/MH)**2))

C  * H -> 2Eta'
      AH2E=CDABS(-DSQRT(dsqrt(2d0)*GF)/6d0*
     . (((MPI**2+MKC**2-MK0**2)*CU(IND)+(MPI**2+MK0**2-MKC**2)*CD(IND))
     . *PropII(MH,MS2,GamS2)*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))**2
     .   +4d0*(MKC**2+MK0**2-MPI**2)*CD(IND)*PropII(MH,MS3,GamS3)
     .              *(DDSIN(THETAE)-DDCOS(THETAE)/dsqrt(2d0))**2)
     .        -2d0*ASH2/PI/3d0*CJH*(PropII(MH,MSIG,GamSIG)
     .    *(MH**2-2d0*METAP**2+MPI**2*(DDSIN(THETAE)
     .       +dsqrt(2d0)*DDCOS(THETAE))**2+2d0*(MKC**2+MK0**2-MPI**2)
     .          *(DDSIN(THETAE)-DDCOS(THETAE)/dsqrt(2d0))**2)
     .     +PropI(MH,MG0,GamG0)*DDCOS(THETAE)**2))**2

      GamH2EP=AH2E/(32d0*PI*MH)*dsqrt(max(0d0,1d0-4d0*(METAP/MH)**2))

C  * H -> EtaEta'
      AH2E=CDABS(-DSQRT(dsqrt(2d0)*GF)/6d0*
     . (((MPI**2+MKC**2-MK0**2)*CU(IND)+(MPI**2+MK0**2-MKC**2)*CD(IND))
     .  *PropII(MH,MS2,GamS2)*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
     .              *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
     .   +4d0*(MKC**2+MK0**2-MPI**2)*CD(IND)*PropII(MH,MS3,GamS3)
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

c              - 3 meson channels (CP-odd modes)

C   Effective yukawas
      YUE=MPI**2*CUP(IND)*DSQRT(2d0*dsqrt(2d0)*GF/3d0)/FPI
      YDE=MPI**2*CDP(IND)*DSQRT(2d0*dsqrt(2d0)*GF/3d0)/FPI
      YSE=(MKC**2+MK0**2-MPI**2)*CDP(IND)*DSQRT(2d0*dsqrt(2d0)*GF/3d0)
     .                                  /FPI

C  * A -> 3Pi0
      IF(MH.lt.1.1d0)then
       AA3Pi=-DSQRT(dsqrt(2d0)*GF)*MPI**2/(2d0*FPI)*(CUP(IND)-CDP(IND))
      ELSE
       IF(MH.lt.1.3d0)then
        AA3Pi=-DSQRT(dsqrt(2d0)*GF)*MPI**2/(2d0*FPI)
     .        *(CUP(IND)-CDP(IND))*(1d0-(MH-1.1d0)/.2d0)
     . +dsqrt(6d0*(YUE**2+YDE**2)/16d0)*SIGN(1d0,CDP(IND)-CUP(IND))
     .        *(MH-1.1d0)/.2d0
       ELSE
        AA3Pi=dsqrt(6d0*(YUE**2+YDE**2)/16d0)
     .      *SIGN(1d0,CDP(IND)-CUP(IND))
       ENDIF
      ENDIF
      APi3Pi=(MPI/FPI)**2
      AE3Pi=0.722665d0
      AEP3Pi=0.278538d0

      aux=OMIX(INDA,1)*APi3Pi+OMIX(INDA,2)*AE3Pi+OMIX(INDA,3)*AEP3Pi
     .          +OMIX(INDA,4)*AA3Pi
      GamA3Pi=aux**2/(6d0*2d0**8*Pi**3*MH)*KinA3Pi(MH)

C  * A -> Pi0Pi+Pi-
      IF(MH.lt.1.1d0)then
       AAPi3PiC=-DSQRT(dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .                                 *(CUP(IND)-CDP(IND))
      ELSE
       IF(MH.lt.1.3d0)then
        AAPi3PiC=-DSQRT(dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .        *(CUP(IND)-CDP(IND))*(1d0-(MH-1.1d0)/.2d0)
     . +dsqrt((YUE**2+YDE**2)/24d0)*SIGN(1d0,CDP(IND)-CUP(IND))
     .        *(MH-1.1d0)/.2d0
       ELSE
        AAPi3PiC=dsqrt((YUE**2+YDE**2)/24d0)
     .       *SIGN(1d0,CDP(IND)-CUP(IND))
       ENDIF
      ENDIF
      APiPi3PiC=(MPI/FPI)**2/3d0
      AEPi3PiC=0.262726d0
      AEPPi3PiC=0.151766d0

      aux=OMIX(INDA,1)*APiPi3PiC+OMIX(INDA,2)*AEPi3PiC+OMIX(INDA,3)
     .          *AEPPi3PiC+OMIX(INDA,4)*AAPi3PiC
      GamAPi3PiC=aux**2/(2d0**8*Pi**3*MH)*KinAPi3PiC(MH)

C  * A -> Eta2Pi0
      IF(MH.lt.1.1d0)then
       AAEPi3=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .                                  *(CUP(IND)+CDP(IND))
      ELSE
       IF(MH.lt.1.3d0)then
        AAEPi3=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .        *(CUP(IND)+CDP(IND))*(1d0-(MH-1.1d0)/.2d0)
     . -dsqrt(2d0*(YUE**2+YDE**2)/16d0)*SIGN(1d0,CDP(IND)+CUP(IND))
     .           *(MH-1.1d0)/.2d0
       ELSE
        AAEPi3=-dsqrt(2d0*(YUE**2+YDE**2)/16d0)
     .        *SIGN(1d0,CDP(IND)+CUP(IND))
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
      GamAEPi3=aux**2/(2d0*2d0**8*Pi**3*MH)*KinAEPi3(MH)

C  * A -> EtaPi+Pi-
      IF(MH.lt.1.1d0)then
       AAEPiC=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .                                  *(CUP(IND)+CDP(IND))
      ELSE
       IF(MH.lt.1.3d0)then
        AAEPiC=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .        *(CUP(IND)+CDP(IND))*(1d0-(MH-1.1d0)/.2d0)
     . -dsqrt((YUE**2+YDE**2)/8d0)*SIGN(1d0,CDP(IND)+CUP(IND))
     .             *(MH-1.1d0)/.2d0
       ELSE
        AAEPiC=-dsqrt((YUE**2+YDE**2)/8d0)*SIGN(1d0,CDP(IND)+CUP(IND))
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
      GamAEPiC=aux**2/(2d0**8*Pi**3*MH)*KinAEPiC(MH)

C  * A -> Eta'2Pi0
      IF(MH.lt.1.1d0)then
       AAEPPi3=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .                                  *(CUP(IND)+CDP(IND))
      ELSE
       IF(MH.lt.1.3d0)then
        AAEPPi3=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .        *(CUP(IND)+CDP(IND))*(1d0-(MH-1.1d0)/.2d0)
     . -dsqrt(2d0*(YUE**2+YDE**2)/16d0)*SIGN(1d0,CDP(IND)+CUP(IND))
     .             *(MH-1.1d0)/.2d0
       ELSE
        AAEPPi3=-dsqrt(2d0*(YUE**2+YDE**2)/16d0)
     .            *SIGN(1d0,CDP(IND)+CUP(IND))
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
      GamAEPPi3=aux**2/(2d0*2d0**8*Pi**3*MH)*KinAEPPi3(MH)

C  * A -> Eta'Pi+Pi-
      IF(MH.lt.1.1d0)then
       AAEPPiC=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .                                      *(CUP(IND)+CDP(IND))
      ELSE
       IF(MH.lt.1.3d0)then
        AAEPPiC=-DSQRT(3d0*dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .        *(CUP(IND)+CDP(IND))*(1d0-(MH-1.1d0)/.2d0)
     . -dsqrt((YUE**2+YDE**2)/8d0)*SIGN(1d0,CDP(IND)+CUP(IND))
     .           *(MH-1.1d0)/.2d0
       ELSE
        AAEPPiC=-dsqrt((YUE**2+YDE**2)/8d0)*SIGN(1d0,CDP(IND)+CUP(IND))
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
      GamAEPPiC=aux**2/(2d0**8*Pi**3*MH)*KinAEPPiC(MH)

C  * A -> Pi02Eta
      IF(MH.lt.1.1d0)then
       AAPiEE=-DSQRT(dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .                                 *(CUP(IND)-CDP(IND))
      ELSE
       IF(MH.lt.1.3d0)then
        AAPiEE=-DSQRT(dsqrt(2d0)*GF)*MPI**2/(6d0*FPI)
     .        *(CUP(IND)-CDP(IND))*(1d0-(MH-1.1d0)/.2d0)
     . +dsqrt(2d0*(YUE**2+YDE**2)/48d0)*SIGN(1d0,CDP(IND)-CUP(IND))
     .           *(MH-1.1d0)/.2d0
       ELSE
        AAPiEE=dsqrt(2d0*(YUE**2+YDE**2)/48d0)
     .     *SIGN(1d0,CDP(IND)+CUP(IND))
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
      GamAPiEE=aux**2/(2d0*2d0**8*Pi**3*MH)*KinAPiEE(MH)

C  * A -> Pi0EtaEta'
      AAPiEEP=dsqrt((YUE**2+YDE**2)/24d0)*SIGN(1d0,CDP(IND)+CUP(IND))
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
      GamAPiEEP=aux**2/(2d0**8*Pi**3*MH)*KinAPiEEP(MH)

C  * A -> Pi02Eta'
      AAPiEPEP=dsqrt(2d0*(YUE**2+YDE**2)/48d0)
     .        *SIGN(1d0,CDP(IND)+CUP(IND))
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
      GamAPiEPEP=aux**2/(2d0*2d0**8*Pi**3*MH)*KinAPiEPEP(MH)

C  * A -> 3Eta
      AA3E=-dsqrt(6d0*(YUE**2+YDE**2+64d0*YSE**2)/432d0)*SIGN(1d0,
     . CUP(IND)*(MKC**2+MPI**2-MK0**2)+CDP(IND)
     .    *(MK0**2+MPI**2-MKC**2)-8d0*CDP(IND)*(MKC**2+MK0**2-MPI**2))
      API3E=(MKC**2-MK0**2-MPIC**2+MPI**2)/(3d0*dsqrt(3d0)*FPI**2)
!     .        *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))**3
      AE3E=(2d0*MPI**2+16d0*(MKC**2+MK0**2-MPI**2))/(9d0*FPI**2)
      AEP3E=dsqrt(2d0)*(2d0*MPI**2-8d0*(MKC**2+MK0**2-MPI**2))
     .             /(9d0*FPI**2)

      aux=OMIX(INDA,1)*API3E+OMIX(INDA,2)*AE3E+OMIX(INDA,3)
     .          *AEP3E+OMIX(INDA,4)*AA3E
      GamA3E=aux**2/(6d0*2d0**8*Pi**3*MH)*KinA3E(MH)

C  * A -> 2EtaEta'
      AAE2EP=-dsqrt(2d0*(YUE**2+YDE**2+16d0*YSE**2)/72d0)*SIGN(1d0,
     . CUP(IND)*(MKC**2+MPI**2-MK0**2)+CDP(IND)
     .     *(MK0**2+MPI**2-MKC**2)+4d0*CDP(IND)*(MKC**2+MK0**2-MPI**2))
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
      GamAE2EP=aux**2/(2d0*2d0**8*Pi**3*MH)*KinAE2EP(MH)

C  * A -> Eta2Eta'
      AAEEP2=-dsqrt(2d0*(YUE**2+YDE**2+4d0*YSE**2)/36d0)*SIGN(1d0,
     . CUP(IND)*(MKC**2+MPI**2-MK0**2)+CDP(IND)
     .    *(MK0**2+MPI**2-MKC**2)-2d0*CDP(IND)*(MKC**2+MK0**2-MPI**2))
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
      GamAEEP2=aux**2/(2d0*2d0**8*Pi**3*MH)*KinAEEP2(MH)

C  * A -> 3Eta'
      AA3EP=-dsqrt(6d0*(YUE**2+YDE**2+YSE**2)/54d0)*SIGN(1d0,
     . CUP(IND)*(MKC**2+MPI**2-MK0**2)+CDP(IND)
     .    *(MK0**2+MPI**2-MKC**2)+CDP(IND)*(MKC**2+MK0**2-MPI**2))
      API3EP=(MKC**2-MK0**2-MPIC**2+MPI**2)/(3d0*dsqrt(3d0)*FPI**2)
     .        *dsqrt(2d0)**3
!     .        *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))**3
      AE3EP=2d0*dsqrt(2d0)*(2d0*MPI**2-2d0*(MKC**2+MK0**2-MPI**2))
     .             /(9d0*FPI**2)
      AEP3EP=4d0*(2d0*MPI**2+(MKC**2+MK0**2-MPI**2))/(9d0*FPI**2)

      aux=OMIX(INDA,1)*API3EP+OMIX(INDA,2)*AE3EP+OMIX(INDA,3)
     .          *AEP3EP+OMIX(INDA,4)*AA3EP
      GamA3EP=aux**2/(6d0*2d0**8*Pi**3*MH)*KinA3EP(MH)

C  * A -> Pi0K+K-
      IF(MH.lt.1.1d0)then
       AAPiKC=-DSQRT(dsqrt(2d0)*GF)/(6d0*FPI)
     .   *(MKC**2*(2d0*CUP(IND)+CDP(IND))
     .    +(MPI**2-MK0**2)*(2d0*CUP(IND)-CDP(IND)))
      ELSE
       IF(MH.lt.1.3d0)then
        AAPiKC=-DSQRT(dsqrt(2d0)*GF)/(6d0*FPI)
     .    *(MKC**2*(2d0*CUP(IND)+CDP(IND))+(MPI**2-MK0**2)
     .               *(2d0*CUP(IND)-CDP(IND)))*(1d0-(MH-1.1d0)/.2d0)
     . -dsqrt((4d0*YUE**2+YSE**2)/24d0)*SIGN(1d0,MKC**2*(2d0*CUP(IND)
     .    +CDP(IND))+(MPI**2-MK0**2)*(2d0*CUP(IND)-CDP(IND)))
     .               *(MH-1.1d0)/.2d0
       ELSE
        AAPiKC=-dsqrt((4d0*YUE**2+YSE**2)/24d0)*SIGN(1d0,
     .    MKC**2*(2d0*CUP(IND)+CDP(IND))+(MPI**2-MK0**2)
     .    *(2d0*CUP(IND)-CDP(IND)))
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
      GamAPiKC=aux**2/(2d0**8*Pi**3*MH)*KinAPiKC(MH)

C  * A -> Pi0K0K0b
      IF(MH.lt.1.1d0)then
       AAPiK0=-DSQRT(dsqrt(2d0)*GF)/(6d0*FPI)*CDP(IND)
     .                              *(MKC**2-MPI**2-3d0*MK0**2)
      ELSE
       IF(MH.lt.1.3d0)then
        AAPiK0=-DSQRT(dsqrt(2d0)*GF)/(6d0*FPI)*CDP(IND)
     .                              *(MKC**2-MPI**2-3d0*MK0**2)
     .               *(1d0-(MH-1.1d0)/.2d0)
     . -dsqrt((4d0*YDE**2+YSE**2)/24d0)*SIGN(1d0,
     .    CDP(IND)*(MKC**2-MPI**2-3d0*MK0**2))
     .    
     .               *(MH-1.1d0)/.2d0
       ELSE
        AAPiK0=-dsqrt((4d0*YDE**2+YSE**2)/24d0)
     .    *SIGN(1d0,CDP(IND)*(MKC**2-MPI**2-3d0*MK0**2))
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
      GamAPiK0=aux**2/(2d0**8*Pi**3*MH)*KinAPiK0(MH)

C  * A -> PiK+K0b + PiK-K0
      IF(MH.lt.1.1d0)then
       AAPiKCK0=-DSQRT(2d0*dsqrt(2d0)*GF)/(6d0*FPI)
     .     *(CUP(IND)*(MKC**2+MPI**2-MK0**2)+2d0*CDP(IND)*MK0**2)
      ELSE
       IF(MH.lt.1.3d0)then
        AAPiKCK0=-DSQRT(2d0*dsqrt(2d0)*GF)/(6d0*FPI)
     .     *(CUP(IND)*(MKC**2+MPI**2-MK0**2)+2d0*CDP(IND)*MK0**2)
     .               *(1d0-(MH-1.1d0)/.2d0)
     . -dsqrt((YUE**2+YDE**2+YSE**2)/12d0)
     .  *SIGN(1d0,CUP(IND)*(MKC**2+MPI**2-MK0**2)+2d0*CDP(IND)*MK0**2)
     .    *(MH-1.1d0)/.2d0
       ELSE
        AAPiKCK0=-dsqrt((YUE**2+YDE**2+YSE**2)/12d0)
     .  *SIGN(1d0,CUP(IND)*(MKC**2+MPI**2-MK0**2)+2d0*CDP(IND)*MK0**2)
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
      GamAPiKCK0=2d0*aux**2/(2d0**8*Pi**3*MH)*KinAPiKCK0(MH)

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
      GamAEKC=aux**2/(2d0**8*Pi**3*MH)*KinAEKC(MH)

C  * A -> EtaK0K0b
      AAEK0=dsqrt((YSE**2)/8d0)
     .     *SIGN(1d0,CDP(IND)*(MKC**2+MK0**2-MPI**2))
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
      GamAEK0=aux**2/(2d0**8*Pi**3*MH)*KinAEK0(MH)

C  * A -> Eta'K+K-
      AAEPKC=-dsqrt((YUE**2+YSE**2)/4d0)*SIGN(1d0,CUP(IND)*(MPI**2+
     . MKC**2-MK0**2)+CDP(IND)*(MK0**2+MPI**2-MKC**2))
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
      GamAEPKC=aux**2/(2d0**8*Pi**3*MH)*KinAEPKC(MH)

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
      GamAEPK0=aux**2/(2d0**8*Pi**3*MH)*KinAEPK0(MH)

C       Radiative hadronic decays

C  * A -> gamma (Rho,omega->)Pi+Pi-
      AARhogam=0d0
      APIRhogam=1d0
      AERhogam=0.745672d0 !1.22374d0
      AEPRhogam=0.9274d0  !1.39445d0

      aux=OMIX(INDA,1)*APIRhogam+OMIX(INDA,2)*AERhogam+OMIX(INDA,3)
     .          *AEPRhogam+OMIX(INDA,4)*AARhogam
      GamARhogam=4d0*Pi*ALEM0/(4d0*Pi**2*FPI**3)**2
     .            *aux**2/(2d0**8*Pi**3*MH)*KinARhogam(MH)

C       Total hadronic decays

      GamHhadr=GamH2Pi+GamH2PiC+GamHPiE+GamHPiEP+GamH2KC+GamH2K0
     .         +GamH2E+GamH2EP+GamHEEP

      GamAhadr=GamA3Pi+GamAPi3PiC+GamAEPi3+GamAEPiC+GamAEPPi3+GamAEPPiC
     .   +GamAPiEE+GamAPiEEP+GamAPiEPEP+GamA3E+GamAE2EP+GamAEEP2
     .   +GamA3EP+GamAPiKC+GamAPiK0+GamAPiKCK0+GamAEKC+GamAEK0
     .   +GamAEPKC+GamAEPK0+GamARhogam

C       Diphoton decays - CP-even mode

      GamHGAGA=CDABS(ALEM0*CGH/Sqrt(2d0)/Pi
     .         +DSQRT(dsqrt(2d0)*GF)/2d0*PropI(MH,MS2,GamS2)
     . *((MPI**2+MKC**2-MK0**2)*CU(IND)+(MPI**2+MK0**2-MKC**2)*CD(IND))
     .   *dsqrt(16d0*Pi*6.7d-7/MS2**3)/2.5d0
     .         +DSQRT(dsqrt(2d0)*GF)/2d0*PropI(MH,MS3,GamS3)
     . *(MKC**2+MK0**2-MPI**2)*CD(IND)*dsqrt(16d0*Pi*4d-7/MS3**3)/2.7d0
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

      AHuu=2d-3*CU(IND)*DSQRT(dsqrt(2d0)*GF)
      AHdd=4d-3*CD(IND)*DSQRT(dsqrt(2d0)*GF)
      AHss=MS*CD(IND)*DSQRT(dsqrt(2d0)*GF)

      AAuu=2d-3*CUP(IND)*DSQRT(dsqrt(2d0)*GF)
      AAdd=4d-3*CDP(IND)*DSQRT(dsqrt(2d0)*GF)
      AAss=MS*CDP(IND)*DSQRT(dsqrt(2d0)*GF)

      IF(MH.LT.2d0)THEN
       GamHuu=0d0
       GamHdd=0d0
       GamHss=0d0
       GamHjj=0d0
       GamAuu=0d0
       GamAdd=0d0
       GamAss=0d0
       GamAjj=0d0
      ELSE

       RATCOUP=1d0
       GamHuu=3d0*dabs(AHuu)**2/(8d0*Pi)*MH
     .        *(4d0*(2d-3/MH)**2*TQCDH((2d-3/MH)**2)
     .   +(1d0-4d0*(2d-3/MH)**2)*max(0d0,QCDH((2d-3/MH)**2)))
     .               *dsqrt(max(0d0,1d0-4d0*(2d-3/MH)**2))**3

       GamAuu=3d0*dabs(AAuu)**2/(8d0*Pi)*MH
     .        *(4d0*(2d-3/MH)**2*TQCDA((2d-3/MH)**2)
     .   +(1d0-4d0*(2d-3/MH)**2)*max(0d0,QCDA((2d-3/MH)**2)))
     .               *dsqrt(max(0d0,1d0-4d0*(2d-3/MH)**2))

       RATCOUP=0d0
       IF(CD(IND).NE.0d0)RATCOUP=CU(IND)/CD(IND)
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

       RATCOUP=0d0
       IF(CDP(IND).NE.0d0)RATCOUP=CUP(IND)/CDP(IND)
       GamAdd=3d0*dabs(AAdd)**2/(8d0*Pi)*MH
     .         *(4d0*(4d-3/MH)**2*TQCDA((4d-3/MH)**2)
     .   +(1d0-4d0*(4d-3/MH)**2)*max(0d0,QCDA((4d-3/MH)**2)))
     .               *dsqrt(max(0d0,1d0-4d0*(4d-3/MH)**2))
       GamAss=3d0*dabs(AAss)**2/(8d0*Pi)*MH
     .         *(4d0*(MS/MH)**2*TQCDA((MS/MH)**2)
     .   +(1d0-4d0*(MS/MH)**2)*max(0d0,(RMS/MS)**2*QCDA((RMS/MH)**2)))
     .               *dsqrt(max(0d0,1d0-4d0*(MS/MH)**2))
       GamAjj=CDABS(CJA)**2/(8d0*PI**3)*MH**3
     .     *AS3**2*(1d0+AS3/Pi*(97d0/4d0-3d0*7d0/6d0))
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

      aux=GamAhadr
      IF(MH.lt.3d0)then
       GamAhadr=aux
      ELSE
       IF(MH.lt.4d0)THEN
        GamAhadr=aux*(4d0-MH)+(MH-3d0)*(0d0*GamAuu+0d0*GamAdd
     .                                         +GamAss+GamAjj)
       ELSE
        GamAhadr=0d0*GamAuu+0d0*GamAdd+GamAss+GamAjj
       ENDIF
      ENDIF

C       Mixing with the Eta_c(1S), chi_c(1P)

      METAC=2.9834d0
      MCHIC=3.41475d0
! eta_c(1S) wave-function from J/Psi(1S)->e+e-
      DMC=dsqrt(3d0/4d0*METAC**3*0.05971d0*92.9d-6*2d0*dsqrt(2d0)*GF) 
     .    *3.0969d0*CUP(IND)/(4d0*Pi*ALEM0*2d0/3d0)
      DMCP=dsqrt(27d0*dsqrt(2d0)/Pi*GF*MCHIC*0.1d0)*CU(IND)    !0.007d0*CU(1)

c      IF(MH.gt.2d0.and.MH.lt.2d0*MTAU)THEN
       MMIXC(1,1)=METAC**2
       MMIXC(1,2)=0d0
       MMIXC(1,3)=DMC
       MMIXC(2,1)=0d0
       MMIXC(2,2)=MCHIC**2
       MMIXC(2,3)=DMCP
       MMIXC(3,1)=DMC
       MMIXC(3,2)=DMCP
       MMIXC(3,3)=MH**2

       CALL DIAGN(3,MMIXC,VALPC,VECPC,EPS)
       CALL SORTN(3,VALPC,VECPC)
       DO K=1,3
        DO L=1,3
         OMIXC(K,L)=VECPC(L,K)
        ENDDO
       ENDDO

       INDA=1
       DO K=1,3
        IF(dabs(OMIXC(INDA,3))**2.lt.dabs(OMIXC(K,3))**2)INDA=K
       ENDDO

       aux=Max(0d0,min(1d0,MH-1d0))
       GamAhadr=GamAhadr*OMIXC(INDA,3)**2+aux*31.8d-3*OMIXC(INDA,1)**2
       GamHhadr=GamHhadr*OMIXC(INDA,3)**2+aux*10.5d-3*OMIXC(INDA,2)**2
       GamAee=GamAee*OMIXC(INDA,3)**2
       GamHee=GamHee*OMIXC(INDA,3)**2
       GamAmumu=GamAmumu*OMIXC(INDA,3)**2
       GamHmumu=GamHmumu*OMIXC(INDA,3)**2
       GamAtata=GamAtata*OMIXC(INDA,3)**2
       GamHtata=GamHtata*OMIXC(INDA,3)**2
       GamAgaga=GamAgaga*OMIXC(INDA,3)**2
     .          +aux*31.8d-3*1.59d-4*OMIXC(INDA,1)**2
       GamHgaga=GamHgaga*OMIXC(INDA,3)**2
     .          +aux*10.5d-3*2.23d-4*OMIXC(INDA,2)**2
       GamHinv=GamHinv*OMIXC(INDA,3)**2
! GamAhadr also contains contributions to the cc final state
c      ENDIF

c      ENDIF
! We neglect the mixing with the eta_c(2,3S) because they matter after the opening of the tautau channel

C       H -> cc decays

      MD=1.865d0

      IF(MH.LE.2d0*MD)THEN
       GamHcc= 0d0
      ELSE
       RMC=RUNM(MH,4)
       RATCOUP= 1d0

* New July 2019:
       HCC=MH**3/(8d0*PI**3)*(AS4**2*Max(0d0,CDABS(CJH)**2*
     .     (1d0+AS4/Pi*(95d0/4d0-4d0*7d0/6d0))+CIH*AS4/PI*7d0/2d0)
     .   -AS3**2*Max(0d0,CDABS(CJH)**2*
     .    (1d0+AS3/Pi*(95d0/4d0-3d0*7d0/6d0))+CIH*AS3/PI*7d0/2d0))

       GamHcc=4d0*(MCC/MH)**2*
     .        3d0*GF*(MCC*CU(IND))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .        *TQCDH((MCC/MH)**2)*dsqrt(1d0-4d0*MCC**2/MH**2)**3
     .         +(1d0-4d0*(MCC/MH)**2)*max(0d0,
     .         3d0*GF*(RMC*CU(IND))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .         *QCDH((RMC/MH)**2)*dsqrt(1d0-4d0*RMC**2/MH**2)**3 !+HCC)
     .        )
* End new

      IF(MH.le.50d0)then
       GamHcc=GamHcc*((1d0-(MH-2d0*MD)/(50d0-2d0*MD))*
     . (dsqrt(1d0-4d0*MD**2/MH**2)/dsqrt(1d0-4d0*MCC**2/MH**2))**3
     . +(MH-2d0*MD)/(50d0-2d0*MD))
* New July 2019:
        GamHhadr=GamHhadr+HCC*
     .  ((1d0-4d0*(MCC/MH)**2)*(1d0-(MH-2d0*MD)/(50d0-2d0*MD))*
     .  (dsqrt(1d0-4d0*MD**2/MH**2)/dsqrt(1d0-4d0*MCC**2/MH**2))**3
     .  +(MH-2d0*MD)/(50d0-2d0*MD))
* End new
      ENDIF
      ENDIF

      IF(MH.LE.2d0*MD+MPI)THEN
       GamAcc= 0d0
      ELSE
       RMC=RUNM(MH,4)
       RATCOUP= 1d0
       ACC=CDABS(CJA)**2/(8d0*PI**3)*MH**3
     .     *(AS4**2*(1d0+AS4/Pi*(97d0/4d0-4d0*7d0/6d0))
     .      -AS3**2*(1d0+AS3/Pi*(97d0/4d0-3d0*7d0/6d0)))

* New July 2019:
       GamAcc=4d0*(MCC/MH)**2*
     .        3d0*GF*(MCC*CUP(IND))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .        *TQCDA((MCC/MH)**2)*dsqrt(1d0-4d0*MCC**2/MH**2)
     .         +(1d0-4d0*(MCC/MH)**2)*max(0d0,
     .         3d0*GF*(RMC*CUP(IND))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .        *QCDA((RMC/MH)**2)*dsqrt(1d0-4d0*RMC**2/MH**2) !+ACC)
     .       )
* End new
! Modifying the shape of the phase space function, so that it looks "3-body" at low mass
       IF(MH.le.50d0)THEN
        GamAcc=GamAcc*dsqrt(1d0-(2d0*MD+MPI)**2/MH**2)
     .                /dsqrt(1d0-4d0*MCC**2/MH**2)
     . *(KinAPiDD(MH)*((50d0-MH)/(50d0-(2d0*MD+MPI)))**2
     .                 +1d0-((50d0-MH)/(50d0-(2d0*MD+MPI)))**2)
* New July 2019:
        GamAhadr=GamAhadr+ACC*dsqrt(1d0-4d0*MCC**2/MH**2)
     . *(dsqrt(1d0-(2d0*MD+MPI)**2/MH**2)
     . *(KinAPiDD(MH)-1d0)*((50d0-MH)/(50d0-(2d0*MD+MPI)))**2+1d0)
* End new
       ENDIF
      ENDIF

C       Mixing with the Eta_b(1,2,3S), chi_b(1,2P)

      METAB1=9.399d0
      METAB2=9.999d0
      METAB3=10.343d0
      MCHIB1=9.85944d0
      MCHIB2=10.2325d0

      DMB1=0.14d0*CBP(IND)
      DMB2=0.11d0*CBP(IND)
      DMB3=0.10d0*CBP(IND)
      DMB1P=dsqrt(27d0*dsqrt(2d0)/Pi*GF*MCHIB1*1.7d0)*CB(IND) ! 0.049d0*CB(1)
      DMB2P=dsqrt(27d0*dsqrt(2d0)/Pi*GF*MCHIB2*2.0d0)*CB(IND) ! 0.054d0*CB(1)

c      IF(MA.gt.4d0.and.MA.lt.2d0*5.2795d0+MPI)THEN

       MMIX2(1,1)=METAB1**2
       MMIX2(1,2)=0d0
       MMIX2(1,3)=0d0
       MMIX2(1,4)=0d0
       MMIX2(1,5)=0d0
       MMIX2(1,6)=DMB1
       MMIX2(2,1)=0d0
       MMIX2(2,2)=METAB2**2
       MMIX2(2,3)=0d0
       MMIX2(2,4)=0d0
       MMIX2(2,5)=0d0
       MMIX2(2,6)=DMB2
       MMIX2(3,1)=0d0
       MMIX2(3,2)=0d0
       MMIX2(3,3)=METAB3**2
       MMIX2(3,4)=0d0
       MMIX2(3,5)=0d0
       MMIX2(3,6)=DMB3
       MMIX2(4,1)=0d0
       MMIX2(4,2)=0d0
       MMIX2(4,3)=0d0
       MMIX2(4,4)=MCHIB1**2
       MMIX2(4,5)=0d0
       MMIX2(4,6)=DMB1P
       MMIX2(5,1)=0d0
       MMIX2(5,2)=0d0
       MMIX2(5,3)=0d0
       MMIX2(5,4)=0d0
       MMIX2(5,5)=MCHIB2**2
       MMIX2(5,6)=DMB2P
       MMIX2(6,1)=DMB1
       MMIX2(6,2)=DMB2
       MMIX2(6,3)=DMB3
       MMIX2(6,4)=DMB1P
       MMIX2(6,5)=DMB2P
       MMIX2(6,6)=MH**2

       CALL DIAGN(6,MMIX2,VALP2,VECP2,EPS)
       CALL SORTN(6,VALP2,VECP2)
       DO K=1,6
        DO L=1,6
         OMIX2(K,L)=VECP2(L,K)
        ENDDO
       ENDDO

       INDA=1
       DO K=1,6
        IF(dabs(OMIX2(INDA,6))**2.lt.dabs(OMIX2(K,6))**2)INDA=K
       ENDDO

       aux=Max(0d0,min(1d0,MH/2d0-1.5d0))
       GamAhadr=GamAhadr*OMIX2(INDA,6)**2+aux*(11.8d-3*OMIX2(INDA,1)**2
     .                +5.4d-3*OMIX2(INDA,2)**2+3.9d-3*OMIX2(INDA,3)**2)
       GamHhadr=GamHhadr*OMIX2(INDA,6)**2+aux*(
     .               2.03d-3*OMIX2(INDA,4)**2+2.39d-3*OMIX2(INDA,5)**2)
       GamAee=GamAee*OMIX2(INDA,6)**2
       GamHee=GamHee*OMIX2(INDA,6)**2
       GamAmumu=GamAmumu*OMIX2(INDA,6)**2
       GamHmumu=GamHmumu*OMIX2(INDA,6)**2
       GamAtata=GamAtata*OMIX2(INDA,6)**2
       GamHtata=GamHtata*OMIX2(INDA,6)**2
       GamAgaga=GamAgaga*OMIX2(INDA,6)**2
     .          +aux*(11.8d-3*3.42d-3*OMIX2(INDA,1)**2
     .                +5.4d-3*3.38d-3*OMIX2(INDA,2)**2
     .                +3.9d-3*3.40d-3*OMIX2(INDA,3)**2)
       GamHgaga=GamHgaga*OMIX2(INDA,6)**2+aux*(
     .                 0.12d-3*OMIX2(INDA,4)**2
     .                +0.14d-3*OMIX2(INDA,5)**2)
       GamAcc=GamAcc*OMIX2(INDA,6)**2
       GamHcc=GamHcc*OMIX2(INDA,6)**2
       GamHinv=GamHinv*OMIX2(INDA,6)**2
! GamAhadr also contains contributions to the cc final state
c      ENDIF

C       H -> bb decays

      MBm=5.2795d0

      IF(MH.LE.2d0*MBm)THEN
       GamHbb= 0d0
      ELSE
       RMB=RUNMB(MH)
       IF(CB(IND).ne.0d0)THEN
        RATCOUP= CU(IND)/CB(IND)
       ELSE
        RATCOUP=0d0
       ENDIF

* New July 2019:
       HBB=MH**3/(8d0*PI**3)*(ASH**2*Max(0d0,CDABS(CJH)**2*
     .    (1d0+ASH/Pi*(95d0/4d0-5d0*7d0/6d0))+CIH*ASH/PI*7d0/2d0)
     .   -AS4**2*Max(0d0,CDABS(CJH)**2*
     .    (1d0+AS4/Pi*(95d0/4d0-4d0*7d0/6d0))+CIH*AS4/PI*7d0/2d0))

       HBB=HBB+MH**3/(8d0*PI**3)*ASH**2*CDABS(CJHF)**2*(
     .     HGGQCD2(ASH,5,MH,MT)-(1d0+ASH/Pi*(95d0/4d0-5d0*7d0/6d0)))

       GamHbb=4d0*(MBP/MH)**2*
     .        3d0*GF*(MBP*CB(IND))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .        *TQCDH((MBP/MH)**2)*dsqrt(1d0-4d0*MBP**2/MH**2)**3
     .         +(1d0-4d0*(MBP/MH)**2)*max(0d0,
     .         3d0*GF*(RMB*CB(IND))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .         *QCDH((RMB/MH)**2)*dsqrt(1d0-4d0*RMB**2/MH**2)**3 !+HBB)
     .         )
* End new
       aux=GamHbb
       IF(MH.lt.50d0)THEN
        GamHbb=GamHbb*((1d0-(MH-2d0*MBm)/(50d0-2d0*MBm))*
     . (dsqrt(1d0-4d0*MBm**2/MH**2)/dsqrt(1d0-4d0*MBP**2/MH**2))**3
     . +(MH-2d0*MBm)/(50d0-2d0*MBm))
* New July 2019:
        GamHhadr=GamHhadr+HBB*
     .  ((1d0-4d0*(MBP/MH)**2)*(1d0-(MH-2d0*MBm)/(50d0-2d0*MBm))*
     .  (dsqrt(1d0-4d0*MBm**2/MH**2)/dsqrt(1d0-4d0*MBP**2/MH**2))**3
     .  +(MH-2d0*MBm)/(50d0-2d0*MBm))
* End new
       ENDIF
      ENDIF

      IF(MH.LE.2d0*MBm+MPI)THEN
       GamAbb= 0d0
      ELSE
       RMB=RUNMB(MH)
       IF(CBP(IND).ne.0d0)THEN
        RATCOUP= CUP(IND)/CBP(IND)
       ELSE
        RATCOUP=0d0
       ENDIF

* New July 2019:
       ABB=CDABS(CJA)**2/(8d0*PI**3)*MH**3
     .     *(ASH**2*AGGQCD2(ASH,5,MH,MT)
     .      -AS4**2*(1d0+AS4/Pi*(97d0/4d0-4d0*7d0/6d0)))

       GamAbb=4d0*(MBP/MH)**2*
     .        3d0*GF*(MBP*CBP(IND))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .        *TQCDA((MBP/MH)**2)*dsqrt(1d0-4d0*MBP**2/MH**2)
     .         +(1d0-4d0*(MBP/MH)**2)*max(0d0,
     .         3d0*GF*(RMB*CBP(IND))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .         *QCDA((RMB/MH)**2)*dsqrt(1d0-4d0*RMB**2/MH**2)!+ABB
     .         )
* End new
! Modifying the shape of the phase space function, so that it looks "3-body" at low mass
       IF(MH.lt.50d0)THEN
        GamAbb=GamAbb*dsqrt(1d0-(2d0*MBm+MPI)**2/MH**2)
     .                /dsqrt(1d0-4d0*MBP**2/MH**2)
     . *(KinAPiBB(MH)*((50d0-MH)/(50d0-(2d0*MBm+MPI)))**2
     .                 +1d0-((50d0-MH)/(50d0-(2d0*MBm+MPI)))**2)
* New July 2019:
        GamAhadr=GamAhadr+ABB*dsqrt(1d0-4d0*(MBP/MH)**2)
     . *(dsqrt(1d0-(2d0*MBm+MPI)**2/MH**2)
     . *(KinAPiBB(MH)-1d0)*((50d0-MH)/(50d0-(2d0*MBm+MPI)))**2+1d0)
* End new
       ENDIF
      ENDIF

       GamHGAGA=GamHGAGA+GamAGAGA
       GamHee=GamHee+GamAee
       GamHmumu=GamHmumu+GamAmumu
       GamHtata=GamHtata+GamAtata
       GamHhadr=GamHhadr+GamAhadr
       GamHcc=GamHcc+GamAcc
       GamHbb=GamHbb+GamAbb
       GamHjj=GamHjj+GamAjj
       HCC=HCC+ACC
       HBB=HBB+ABB

      RETURN
      END
