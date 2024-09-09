      SUBROUTINE EDM_CPV(PAR,PROB)

c     Estimates the EDMs

      IMPLICIT NONE

      INTEGER I,J,M

      DOUBLE PRECISION PAR(*),PROB(*),Pi
      DOUBLE PRECISION errdEe,dEu,errdEu,dCu,errdCu,dEd,errdEd,
     . dCd,errdCd,dEs,errdEs
      DOUBLE PRECISION Cee,Cdd,Cuu,Ced,Cde,Cue,Ceu,Cdu,Cud,Ces,Cse
      DOUBLE PRECISION Cec,Cet,Ceb,Cbd,Cdb,Csd,CSG,dGW
      DOUBLE PRECISION errCSG,errdGW
      DOUBLE PRECISION AEDM,BEDM,CEDM,fSF,fS,fPS,fWW,HEDM,HQEDM
      DOUBLE PRECISION aux,aux1,aux2,DELT(2,2),me,mqu,mqd,mqs
      DOUBLE PRECISION USU(2,2,2),USD(2,2,2),USC(2,2,2),USS(2,2,2)
      DOUBLE PRECISION UE(2,2,2),UMU(2,2,2)
      DOUBLE PRECISION MSU2(2),MSD2(2),MSC2(2),MSS2(2),MSE2(2),MSMU2(2)
      DOUBLE PRECISION CONESSsL(5,2,2),CONESSsR(5,2,2),COCHSCsL(2,2,2),
     . COCHSCsR(2,2,2)
      DOUBLE PRECISION gRHSCSC(5,2,2),gRHSSSS(5,2,2),gRHSMSM(5,2,2)
      DOUBLE PRECISION MSF2_11,MSF2_22,MSF4_12,phiSF,thetaSF
      DOUBLE PRECISION AU,AD,AC,AS,AE,phiAU,phiAD,phiAE,XHG(5,6)

      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION l,k,Alcos1,Akcos2,muq,nuq
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION MST2(2),UT(2,2,2),MSB2(2),UB(2,2,2),MSL2(2),
     . UTAU(2,2,2),MSNT2
      DOUBLE PRECISION MSU2_CPV(2),MSD2_CPV(2),MSE2_CPV(2),MSNE2,
     . MSMU2_CPV(2),UMU_CPV(2,2,2)
      DOUBLE PRECISION MGL
      DOUBLE PRECISION MSQ1,MSU1,MSD1
      DOUBLE PRECISION MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      DOUBLE PRECISION MCH2(2),U(2,2,2),V(2,2,2)
      DOUBLE PRECISION MNEU(5),NEU(5,5,2)
      DOUBLE PRECISION COCHSTbL(2,2,2),COCHSTbR(2,2,2),
     . COCHSBtL(2,2,2),COCHSBtR(2,2,2),COCHSLnL(2,2,2),
     . COCHSNlL(2,2),COCHSNlR(2,2),COCHSUdL(2,2,2),COCHSUdR(2,2,2),
     . COCHSDuL(2,2,2),COCHSDuR(2,2,2),COCHSEnL(2,2,2),
     . COCHSNeL(2,2),COCHSNeR(2,2)
      DOUBLE PRECISION CONESTtL(5,2,2),CONESTtR(5,2,2),
     . CONESBbL(5,2,2),CONESBbR(5,2,2),CONESLlL(5,2,2),CONESLlR(5,2,2),
     . CONESNnL(5,2),CONESUuL(5,2,2),CONESUuR(5,2,2),
     . CONESDdL(5,2,2),CONESDdR(5,2,2),CONESEeL(5,2,2),
     . CONESEeR(5,2,2)
      DOUBLE PRECISION GRHSTST(5,2,2),GRHSBSB(5,2,2),GRHSLSL(5,2,2),
     . GRHSUSU(5,2,2),GRHSDSD(5,2,2),GRHSESE(5,2,2),GRHSNSN(5)
      DOUBLE PRECISION GIHSTST(5,2,2),GIHSBSB(5,2,2),GIHSLSL(5,2,2)
      DOUBLE PRECISION GRHCSTSB(2,2),GRHCSNSL(2),GRHCSUSD(2,2),
     . GRHCSNSE(2),GIHCSTSB(2,2),GIHCSNSL(2)
      DOUBLE PRECISION COH0CH(5,2,2,2),COH0NEU(5,5,5,2),
     . COHPNEUCHM(2,5,2,2),COHMNEUCHP(2,5,2,2)
      DOUBLE PRECISION CU(5),CUP(5),CD(5),CDP(5),CB(5),CBP(5),CJ(5),
     . CJP(5),CI(5),CG(5),CGP(5),CV(5),CZG(5),CZGP(5),CL(5),CLP(5)
      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION dEe,dETl,dEnCQM,dEnPQM,dEnQSR,dEHg
      DOUBLE PRECISION dEemin,dETlmin,dEnCQMmin,dEnPQMmin,dEnQSRmin,
     . dEHgmin
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/SMFERM/mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/QPAR/l,k,Alcos1,Akcos2,muq,NUQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/SFERM3SPEC/MST2,UT,MSB2,UB,MSL2,UTAU,MSNT2
      COMMON/SFERM1SPEC/MSU2_CPV,MSD2_CPV,MSE2_CPV,MSNE2,
     . MSMU2_CPV,UMU_CPV
      COMMON/GLUSPEC/MGL
      COMMON/SQUPAR/MSQ1,MSU1,MSD1
      COMMON/SLEPPAR/MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      COMMON/CHASPEC/MCH2,U,V
      COMMON/NEUSPEC/MNEU,NEU
      COMMON/CHSFfCOUP/COCHSTbL,COCHSTbR,COCHSBtL,COCHSBtR,COCHSLnL,
     . COCHSNlL,COCHSNlR,COCHSUdL,COCHSUdR,COCHSDuL,COCHSDuR,COCHSEnL,
     . COCHSNeL,COCHSNeR
      COMMON/NEUSFfCOUP/CONESTtL,CONESTtR,CONESBbL,CONESBbR,CONESLlL,
     . CONESLlR,CONESNnL,CONESUuL,CONESUuR,CONESDdL,CONESDdR,CONESEeL,
     . CONESEeR
      COMMON/HISFCOUP/GRHSTST,GRHSBSB,GRHSLSL,GRHSUSU,GRHSDSD,
     . GRHSESE,GRHSNSN,GIHSTST,GIHSBSB,GIHSLSL,GRHCSTSB,GRHCSNSL,
     . GRHCSUSD,GRHCSNSE,GIHCSTSB,GIHCSNSL
      COMMON/HINOCOUP/COH0CH,COH0NEU,COHPNEUCHM,COHMNEUCHP
      COMMON/HNSMCOUP/CU,CUP,CD,CDP,CB,CBP,CJ,CJP,CI,CG,CGP,CV,CZG,CZGP,
     . CL,CLP
      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/EDM/dEe,dETl,dEnCQM,dEnPQM,dEnQSR,dEHg,
     . dEemin,dETlmin,dEnCQMmin,dEnPQMmin,dEnQSRmin,dEHgmin

      Pi=4d0*DATAN(1d0)
      me=mel
      mqu=2.5d-3
      mqd=5.0d-3
      mqs=0.101d0

c      DO I=1,5
c       CU(I)= XH(I,1)/SINB
c       CD(I)= XH(I,2)/COSB
c       CUP(I)= XH(I,4)/TANB
c       CDP(I)= XH(I,4)*TANB
c      ENDDO

      DELT(1,1)=1.d0
      DELT(1,2)=0.d0
      DELT(2,1)=0.d0
      DELT(2,2)=1.d0

      DO I=1,2
       DO J=1,2
        USU(I,J,1)=DELT(I,J)
        USU(I,J,2)=0d0
        USD(I,J,1)=DELT(I,J)
        USD(I,J,2)=0d0
        USS(I,J,1)=DELT(I,J)
        USS(I,J,2)=0d0
        COCHSCsL(I,J,1)=COCHSUdL(I,J,1)
        COCHSCsL(I,J,2)=COCHSUdL(I,J,2)
        COCHSCsR(I,J,1)=COCHSUdR(I,J,1)
        COCHSCsR(I,J,2)=COCHSUdR(I,J,2)
       ENDDO
       MSC2(I)=MSU2_CPV(I)
       MSS2(I)=MSD2_CPV(I)
       MSMU2(I)=MSE2_CPV(I)
      ENDDO

      DO I=1,5
       DO J=1,2
        CONESSsL(I,J,1)=CONESDdL(I,J,1)
        CONESSsL(I,J,2)=CONESDdL(I,J,2)
        CONESSsR(I,J,1)=CONESDdR(I,J,1)
        CONESSsR(I,J,2)=CONESDdR(I,J,2)
        DO M=1,2
         gRHSCSC(I,J,M)=gRHSUSU(I,J,M)
         gRHSSSS(I,J,M)=gRHSDSD(I,J,M)
         gRHSMSM(I,J,M)=gRHSESE(I,J,M)
        ENDDO
       ENDDO
      ENDDO

c        0: Masses/couplings for the squarks of 1st/2nd generation

      M=1                           ! IF M=1, will compute new masses
      IF(M.eq.1)THEN
      AU=PAR(12)
      AD=PAR(13)
      AC=PAR(12)
      AS=PAR(13)
      AMU=PAR(25)
      AE=PAR(25)
      phiAU=phiAC
      phiAD=phiAS
      phiAE=PhiAMU

c      I- Sup / Scharm

      MSF2_11=MSQ1+mqu**2+(g1q/3.d0-g2q)/4.d0
     .                                 *(vuq**2-vdq**2)
      MSF2_22=MSU1+mqu**2+(-4.d0/3.d0*g1q)/4.d0
     .                                 *(vuq**2-vdq**2)
      MSF4_12=mqu**2*(AU**2+(muq*vdq/vuq)**2
     .             -2.d0*AU*muq*vdq/vuq*DDCOS(Phi01+phiAU))
     
      aux1=MSF2_11+MSF2_22
      aux2=(MSF2_11-MSF2_22)**2+4.d0*MSF4_12
      MSU2(1)=(aux1-dsqrt(aux2))/2.d0
      MSU2(2)=(aux1+dsqrt(aux2))/2.d0
      
      IF(MSU2(1).le.0.d0)THEN
       MSU2(1)=MSQ1+(g1q/3.d0-g2q)/4.d0*(vuq**2-vdq**2)
       MSU2(2)=MSU1+(-4.d0/3.d0*g1q)/4.d0*(vuq**2-vdq**2)
      ELSE
       aux1=mqu*(-AU*DDSIN(phiAU)-muq*vdq/vuq*DDSIN(Phi01))
       aux2=mqu*(AU*DDCOS(phiAU)-muq*vdq/vuq*DDCOS(Phi01))
       IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
        phiSF=datan(aux1/aux2)
        IF(aux2.le.0.d0)phiSF=phiSF+Pi
       ELSEIF(aux1.ge.0.d0)THEN
        phiSF=Pi/2.d0
       ELSE
        phiSF=-Pi/2.d0
       ENDIF

       aux1=MSQ1-MSU1+(g1q*(1.d0/3.d0-(-4.d0/3.d0))-g2q)
     .           *(vuq**2-vdq**2)/4.d0
     .        +dsqrt((MSF2_11-MSF2_22)**2+4.d0*MSF4_12)
       aux2=dsqrt(4.d0*MSF4_12)
       IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
        thetaSF=datan(aux1/aux2)
        IF(aux2.le.0.d0)thetaSF=thetaSF+Pi
       ELSEIF(aux1.ge.0.d0)THEN
        thetaSF=Pi/2.d0
       ELSE
        thetaSF=-Pi/2.d0
       ENDIF

       USU(1,1,1)=DDCOS(thetaSF)
       USU(1,1,2)=0.d0
       USU(1,2,1)=-DDSIN(thetaSF)*DDCOS(phiSF)
       USU(1,2,2)=-DDSIN(thetaSF)*DDSIN(phiSF)
       USU(2,1,1)=DDSIN(thetaSF)*DDCOS(phiSF)
       USU(2,1,2)=-DDSIN(thetaSF)*DDSIN(phiSF)
       USU(2,2,1)=DDCOS(thetaSF)
       USU(2,2,2)=0.d0
      ENDIF
      
      
      MSF2_11=MSQ1+mc**2+(g1q/3.d0-g2q)/4.d0
     .                                 *(vuq**2-vdq**2)
      MSF2_22=MSU1+mc**2+(-4.d0/3.d0*g1q)/4.d0
     .                                 *(vuq**2-vdq**2)
      MSF4_12=mc**2*(AC**2+(muq*vdq/vuq)**2
     .             -2.d0*AC*muq*vdq/vuq*DDCOS(Phi01+PhiAC))
     
      aux1=MSF2_11+MSF2_22
      aux2=(MSF2_11-MSF2_22)**2+4.d0*MSF4_12
      MSU2(1)=(aux1-dsqrt(aux2))/2.d0
      MSU2(2)=(aux1+dsqrt(aux2))/2.d0
      
      IF(MSC2(1).le.0.d0)THEN
       MSC2(1)=MSQ1+(g1q/3.d0-g2q)/4.d0*(vuq**2-vdq**2)
       MSC2(2)=MSU1+(-4.d0/3.d0*g1q)/4.d0*(vuq**2-vdq**2)
      ELSE
       aux1=mc*(-AC*DDSIN(PhiAC)-muq*vdq/vuq*DDSIN(Phi01))
       aux2=mc*(AC*DDCOS(PhiAC)-muq*vdq/vuq*DDCOS(Phi01))
       IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
        phiSF=datan(aux1/aux2)
        IF(aux2.le.0.d0)phiSF=phiSF+Pi
       ELSEIF(aux1.ge.0.d0)THEN
        phiSF=Pi/2.d0
       ELSE
        phiSF=-Pi/2.d0
       ENDIF

       aux1=MSQ1-MSU1+(g1q*(1.d0/3.d0-(-4.d0/3.d0))-g2q)
     .           *(vuq**2-vdq**2)/4.d0
     .        +dsqrt((MSF2_11-MSF2_22)**2+4.d0*MSF4_12)
       aux2=dsqrt(4.d0*MSF4_12)
       IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
        thetaSF=datan(aux1/aux2)
        IF(aux2.le.0.d0)thetaSF=thetaSF+Pi
       ELSEIF(aux1.ge.0.d0)THEN
        thetaSF=Pi/2.d0
       ELSE
        thetaSF=-Pi/2.d0
       ENDIF

       USC(1,1,1)=DDCOS(thetaSF)
       USC(1,1,2)=0.d0
       USC(1,2,1)=-DDSIN(thetaSF)*DDCOS(phiSF)
       USC(1,2,2)=-DDSIN(thetaSF)*DDSIN(phiSF)
       USC(2,1,1)=DDSIN(thetaSF)*DDCOS(phiSF)
       USC(2,1,2)=-DDSIN(thetaSF)*DDSIN(phiSF)
       USC(2,2,1)=DDCOS(thetaSF)
       USC(2,2,2)=0.d0
      ENDIF
      
c      II- Sdown / Sstrange

      MSF2_11=MSQ1+mqd**2+(g1q/3.d0+g2q)/4.d0
     .                                 *(vuq**2-vdq**2)
      MSF2_22=MSD1+mqd**2+(2.d0/3.d0*g1q)/4.d0
     .                                 *(vuq**2-vdq**2)
      MSF4_12=mqd**2*(AD**2+(muq*vuq/vdq)**2
     .               -2.d0*AD*muq*vuq/vdq*DDCOS(Phi01+PhiAD))
     
      aux1=MSF2_11+MSF2_22
      aux2=(MSF2_11-MSF2_22)**2+4.d0*MSF4_12
      MSD2(1)=(aux1-dsqrt(aux2))/2.d0
      MSD2(2)=(aux1+dsqrt(aux2))/2.d0
      
      IF(MSD2(1).le.0.d0)THEN
      MSD2(1)=MSQ1+(g1q/3.d0+g2q)/4.d0*(vuq**2-vdq**2)
      MSD2(2)=MSD1+(2.d0/3.d0*g1q)/4.d0*(vuq**2-vdq**2)
      ELSE
      aux1=mqd*(-AD*DDSIN(PhiAD)-muq*vuq/vdq*DDSIN(Phi01))
      aux2=mqd*(AD*DDCOS(PhiAD)-muq*vuq/vdq*DDCOS(Phi01))
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       phiSF=datan(aux1/aux2)
       IF(aux2.le.0.d0)PhiSF=PhiSF+Pi
      ELSEIF(aux1.ge.0.d0)THEN
       phiSF=Pi/2.d0
      ELSE
       phiSF=-Pi/2.d0
      ENDIF

      aux1=MSQ1-MSD1+(g1q*(1.d0/3.d0-(2.d0/3.d0))+g2q)
     .           *(vuq**2-vdq**2)/4.d0
     .   +dsqrt((MSF2_11-MSF2_22)**2+4.d0*MSF4_12)
      aux2=dsqrt(4.d0*MSF4_12)
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       thetaSF=datan(aux1/aux2)
       IF(aux2.le.0.d0)thetaSF=thetaSF+Pi
      ELSEIF(aux1.ge.0.d0)THEN
       thetaSF=Pi/2.d0
      ELSE
       thetaSF=-Pi/2.d0
      ENDIF

      USD(1,1,1)=DDCOS(thetaSF)
      USD(1,1,2)=0.d0
      USD(1,2,1)=-DDSIN(thetaSF)*DDCOS(PhiSF)
      USD(1,2,2)=-DDSIN(thetaSF)*DDSIN(PhiSF)
      USD(2,1,1)=DDSIN(thetaSF)*DDCOS(PhiSF)
      USD(2,1,2)=-DDSIN(thetaSF)*DDSIN(PhiSF)
      USD(2,2,1)=DDCOS(thetaSF)
      USD(2,2,2)=0.d0
      ENDIF

      MSF2_11=MSQ1+mqs**2+(g1q/3.d0+g2q)/4.d0
     .                                 *(vuq**2-vdq**2)
      MSF2_22=MSD1+mqs**2+(2.d0/3.d0*g1q)/4.d0
     .                                 *(vuq**2-vdq**2)
      MSF4_12=mqs**2*(AS**2+(muq*vuq/vdq)**2
     .               -2.d0*AS*muq*vuq/vdq*DDCOS(Phi01+PhiAS))
     
      aux1=MSF2_11+MSF2_22
      aux2=(MSF2_11-MSF2_22)**2+4.d0*MSF4_12
      MSS2(1)=(aux1-dsqrt(aux2))/2.d0
      MSS2(2)=(aux1+dsqrt(aux2))/2.d0
      
      IF(MSS2(1).le.0.d0)THEN
       MSS2(1)=MSQ1+(g1q/3.d0+g2q)/4.d0*(vuq**2-vdq**2)
       MSS2(2)=MSD1+(2.d0/3.d0*g1q)/4.d0*(vuq**2-vdq**2)
      ELSE
      aux1=mqs*(-AS*DDSIN(PhiAS)-muq*vuq/vdq*DDSIN(Phi01))
      aux2=mqs*(AS*DDCOS(PhiAS)-muq*vuq/vdq*DDCOS(Phi01))
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       phiSF=datan(aux1/aux2)
       IF(aux2.le.0.d0)PhiSF=PhiSF+Pi
      ELSEIF(aux1.ge.0.d0)THEN
       phiSF=Pi/2.d0
      ELSE
       phiSF=-Pi/2.d0
      ENDIF

      aux1=MSQ1-MSD1+(g1q*(1.d0/3.d0-(2.d0/3.d0))+g2q)
     .           *(vuq**2-vdq**2)/4.d0
     .   +dsqrt((MSF2_11-MSF2_22)**2+4.d0*MSF4_12)
      aux2=dsqrt(4.d0*MSF4_12)
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       thetaSF=datan(aux1/aux2)
       IF(aux2.le.0.d0)thetaSF=thetaSF+Pi
      ELSEIF(aux1.ge.0.d0)THEN
       thetaSF=Pi/2.d0
      ELSE
       thetaSF=-Pi/2.d0
      ENDIF

      USS(1,1,1)=DDCOS(thetaSF)
      USS(1,1,2)=0.d0
      USS(1,2,1)=-DDSIN(thetaSF)*DDCOS(PhiSF)
      USS(1,2,2)=-DDSIN(thetaSF)*DDSIN(PhiSF)
      USS(2,1,1)=DDSIN(thetaSF)*DDCOS(PhiSF)
      USS(2,1,2)=-DDSIN(thetaSF)*DDSIN(PhiSF)
      USS(2,2,1)=DDCOS(thetaSF)
      USS(2,2,2)=0.d0
      ENDIF

c      III- Selectron / Smuon

      MSF2_11=MSL1+me**2+(-g1q+g2q)/4.d0*(vuq**2-vdq**2)
      MSF2_22=MSE1+me**2+(2.d0*g1q)/4.d0*(vuq**2-vdq**2)
      MSF4_12=me**2*(AE**2+muq**2*(vuq/vdq)**2
     .                -2.d0*AE*muq*vuq/vdq*DDCOS(Phi01+PhiAE))

      aux1=MSF2_11+MSF2_22
      aux2=(MSF2_11-MSF2_22)**2+4.d0*MSF4_12
      MSE2(1)=(aux1-dsqrt(aux2))/2.d0
      MSE2(2)=(aux1+dsqrt(aux2))/2.d0
      IF(MSE2(1).le.0.d0)THEN
       MSE2(1)=MSL1+(-g1q+g2q)/4.d0*(vuq**2-vdq**2)
       MSE2(2)=MSE1+(2.d0*g1q)/4.d0*(vuq**2-vdq**2)
      ELSE

      aux1=me*(-AE*DDSIN(PhiAE)-muq*vuq/vdq*DDSIN(Phi01))
      aux2=me*(AE*DDCOS(PhiAE)-muq*vuq/vdq*DDCOS(Phi01))
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       phiSF=datan(aux1/aux2)
       IF(aux2.le.0.d0)PhiSF=PhiSF+Pi
      ELSEIF(aux1.ge.0.d0)THEN
       phiSF=Pi/2.d0
      ELSE
       phiSF=-Pi/2.d0
      ENDIF

      aux1=MSL1-MSE1+(g1q*(-1.d0-2.d0)+g2q)
     .           *(vuq**2-vdq**2)/4.d0
     .      +dsqrt((MSF2_11-MSF2_22)**2+4.d0*MSF4_12)
      aux2=dsqrt(4.d0*MSF4_12)
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       thetaSF=datan(aux1/aux2)
       IF(aux2.le.0.d0)thetaSF=thetaSF+Pi
      ELSEIF(aux1.ge.0.d0)THEN
       thetaSF=Pi/2.d0
      ELSE
       thetaSF=-Pi/2.d0
      ENDIF

      UE(1,1,1)=DDCOS(thetaSF)
      UE(1,1,2)=0.d0
      UE(1,2,1)=-DDSIN(thetaSF)*DDCOS(PhiSF)
      UE(1,2,2)=-DDSIN(thetaSF)*DDSIN(PhiSF)
      UE(2,1,1)=DDSIN(thetaSF)*DDCOS(PhiSF)
      UE(2,1,2)=-DDSIN(thetaSF)*DDSIN(PhiSF)
      UE(2,2,1)=DDCOS(thetaSF)
      UE(2,2,2)=0.d0
      ENDIF

      MSF2_11=MSL1+mmu**2+(-g1q+g2q)/4.d0*(vuq**2-vdq**2)
      MSF2_22=MSE1+mmu**2+(2.d0*g1q)/4.d0*(vuq**2-vdq**2)
      MSF4_12=mmu**2*(Amu**2+muq**2*(vuq/vdq)**2
     .                -2.d0*Amu*muq*vuq/vdq*DDCOS(Phi01+PhiAmu))

      aux1=MSF2_11+MSF2_22
      aux2=(MSF2_11-MSF2_22)**2+4.d0*MSF4_12
      MSMU2(1)=(aux1-dsqrt(aux2))/2.d0
      MSMU2(2)=(aux1+dsqrt(aux2))/2.d0
      IF(MSMU2(1).le.0.d0)THEN
       MSMU2(1)=MSL1+(-g1q+g2q)/4.d0*(vuq**2-vdq**2)
       MSMU2(2)=MSE1+(2.d0*g1q)/4.d0*(vuq**2-vdq**2)
      ELSE

      aux1=mmu*(-Amu*DDSIN(PhiAmu)-muq*vuq/vdq*DDSIN(Phi01))
      aux2=mmu*(Amu*DDCOS(PhiAmu)-muq*vuq/vdq*DDCOS(Phi01))
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       phiSF=datan(aux1/aux2)
       IF(aux2.le.0.d0)PhiSF=PhiSF+Pi
      ELSEIF(aux1.ge.0.d0)THEN
       phiSF=Pi/2.d0
      ELSE
       phiSF=-Pi/2.d0
      ENDIF

      aux1=MSL1-MSE1+(g1q*(-1.d0-2.d0)+g2q)
     .           *(vuq**2-vdq**2)/4.d0
     .      +dsqrt((MSF2_11-MSF2_22)**2+4.d0*MSF4_12)
      aux2=dsqrt(4.d0*MSF4_12)
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       thetaSF=datan(aux1/aux2)
       IF(aux2.le.0.d0)thetaSF=thetaSF+Pi
      ELSEIF(aux1.ge.0.d0)THEN
       thetaSF=Pi/2.d0
      ELSE
       thetaSF=-Pi/2.d0
      ENDIF

      UMU(1,1,1)=DDCOS(thetaSF)
      UMU(1,1,2)=0.d0
      UMU(1,2,1)=-DDSIN(thetaSF)*DDCOS(PhiSF)
      UMU(1,2,2)=-DDSIN(thetaSF)*DDSIN(PhiSF)
      UMU(2,1,1)=DDSIN(thetaSF)*DDCOS(PhiSF)
      UMU(2,1,2)=-DDSIN(thetaSF)*DDSIN(PhiSF)
      UMU(2,2,1)=DDCOS(thetaSF)
      UMU(2,2,2)=0.d0
      ENDIF

c      IV- Couplings to charginos/neutralinos-SM fermions

        DO I=1,2
        DO J=1,2

      COCHSUdL(I,J,1)=mqu/vuq
     .                 *(V(I,2,1)*USU(J,2,1)+V(I,2,2)*USU(J,2,2))
     .     -dsqrt(g2q)*(V(I,1,1)*USU(J,1,1)+V(I,1,2)*USU(J,1,2))
      COCHSUdL(I,J,2)=mqu/vuq
     .                 *(V(I,2,1)*USU(J,2,2)-V(I,2,2)*USU(J,2,1))
     .     -dsqrt(g2q)*(V(I,1,1)*USU(J,1,2)-V(I,1,2)*USU(J,1,1))
      COCHSUdR(I,J,1)=mqd/vdq
     .                 *(U(I,2,1)*USU(J,1,1)-U(I,2,2)*USU(J,1,2))
      COCHSUdR(I,J,2)=mqd/vdq
     .                 *(U(I,2,2)*USU(J,1,1)+U(I,2,1)*USU(J,1,2))

      COCHSDuL(I,J,1)=mqd/vdq
     .                 *(U(I,2,1)*USD(J,2,1)+U(I,2,2)*USD(J,2,2))
     .     -dsqrt(g2q)*(U(I,1,1)*USD(J,1,1)+U(I,1,2)*USD(J,1,2))
      COCHSDuL(I,J,2)=mqd/vdq
     .                 *(U(I,2,1)*USD(J,2,2)-U(I,2,2)*USD(J,2,1))
     .     -dsqrt(g2q)*(U(I,1,1)*USD(J,1,2)-U(I,1,2)*USD(J,1,1))
      COCHSDuR(I,J,1)=mqu/vuq
     .                 *(V(I,2,1)*USD(J,1,1)-V(I,2,2)*USD(J,1,2))
      COCHSDuR(I,J,2)=mqu/vuq
     .                 *(V(I,2,2)*USD(J,1,1)+V(I,2,1)*USD(J,1,2))

      COCHSCsL(I,J,1)=mc/vuq
     .                 *(V(I,2,1)*USC(J,2,1)+V(I,2,2)*USC(J,2,2))
     .     -dsqrt(g2q)*(V(I,1,1)*USC(J,1,1)+V(I,1,2)*USC(J,1,2))
      COCHSCsL(I,J,2)=mc/vuq
     .                 *(V(I,2,1)*USC(J,2,2)-V(I,2,2)*USC(J,2,1))
     .     -dsqrt(g2q)*(V(I,1,1)*USC(J,1,2)-V(I,1,2)*USC(J,1,1))
      COCHSCsR(I,J,1)=mqs/vdq
     .                 *(U(I,2,1)*USC(J,1,1)-U(I,2,2)*USC(J,1,2))
      COCHSCsR(I,J,2)=mqs/vdq
     .                 *(U(I,2,2)*USC(J,1,1)+U(I,2,1)*USC(J,1,2))

      COCHSNeR(I,1)=me/vdq*U(I,2,1)
      COCHSNeR(I,2)=me/vdq*U(I,2,2)
      
        ENDDO
        ENDDO

        DO I=1,5
        DO J=1,2

      CONESUuL(I,J,1)=-mqu/vuq
     .          *(NEU(I,3,1)*USU(J,2,1)+NEU(I,3,2)*USU(J,2,2))
     . -(dsqrt(g1q)/3.d0*NEU(I,1,1)+dsqrt(g2q)*NEU(I,2,1))
     .                   *USU(J,1,1)/dsqrt(2.d0)
     . -(dsqrt(g1q)/3.d0*NEU(I,1,2)+dsqrt(g2q)*NEU(I,2,2))
     .                   *USU(J,1,2)/dsqrt(2.d0)
      CONESUuL(I,J,2)=-mqu/vuq
     .          *(NEU(I,3,1)*USU(J,2,2)-NEU(I,3,2)*USU(J,2,1))
     . -(dsqrt(g1q)/3.d0*NEU(I,1,1)+dsqrt(g2q)*NEU(I,2,1))
     .                   *USU(J,1,2)/dsqrt(2.d0)
     . +(dsqrt(g1q)/3.d0*NEU(I,1,2)+dsqrt(g2q)*NEU(I,2,2))
     .                   *USU(J,1,1)/dsqrt(2.d0)
      CONESUuR(I,J,1)=-mqu/vuq
     .          *(NEU(I,3,1)*USU(J,1,1)-NEU(I,3,2)*USU(J,1,2))
     . +2.d0*dsqrt(2.d0*g1q)
     .          *(NEU(I,1,1)*USU(J,2,1)-NEU(I,1,2)*USU(J,2,2))
      CONESUuR(I,J,2)=-mqu/vuq
     .          *(NEU(I,3,2)*USU(J,1,1)+NEU(I,3,1)*USU(J,1,2))
     . +2.d0*dsqrt(2.d0*g1q)
     .          *(NEU(I,1,1)*USU(J,2,2)+NEU(I,1,2)*USU(J,2,1))

      CONESDdL(I,J,1)=-mqd/vdq
     .          *(NEU(I,4,1)*USD(J,2,1)+NEU(I,4,2)*USD(J,2,2))
     . -(dsqrt(g1q)/3.d0*NEU(I,1,1)-dsqrt(g2q)*NEU(I,2,1))
     .                  *USD(J,1,1)/dsqrt(2.d0)
     . -(dsqrt(g1q)/3.d0*NEU(I,1,2)-dsqrt(g2q)*NEU(I,2,2))
     .                  *USD(J,1,2)/dsqrt(2.d0)
      CONESDdL(I,J,2)=-mqd/vdq
     .         *(NEU(I,4,1)*USD(J,2,2)-NEU(I,4,2)*USD(J,2,1))
     . -(dsqrt(g1q)/3.d0*NEU(I,1,1)-dsqrt(g2q)*NEU(I,2,1))
     .                  *USD(J,1,2)/dsqrt(2.d0)
     . +(dsqrt(g1q)/3.d0*NEU(I,1,2)-dsqrt(g2q)*NEU(I,2,2))
     .                  *USD(J,1,1)/dsqrt(2.d0)
      CONESDdR(I,J,1)=-mqd/vdq
     .        *(NEU(I,4,1)*USD(J,1,1)-NEU(I,4,2)*USD(J,1,2))
     .-dsqrt(2.d0*g1q)/3.d0
     .        *(NEU(I,1,1)*USD(J,2,1)-NEU(I,1,2)*USD(J,2,2))
      CONESDdR(I,J,2)=-mqd/vdq
     .        *(NEU(I,4,2)*USD(J,1,1)+NEU(I,4,1)*USD(J,1,2))
     .-dsqrt(2.d0*g1q)/3.d0
     .        *(NEU(I,1,1)*USD(J,2,2)+NEU(I,1,2)*USD(J,2,1))

      CONESSsL(I,J,1)=-mqs/vdq
     .        *(NEU(I,4,1)*USS(J,2,1)+NEU(I,4,2)*USS(J,2,2))
     . -(dsqrt(g1q)/3.d0*NEU(I,1,1)-dsqrt(g2q)*NEU(I,2,1))
     .                  *USS(J,1,1)/dsqrt(2.d0)
     . -(dsqrt(g1q)/3.d0*NEU(I,1,2)-dsqrt(g2q)*NEU(I,2,2))
     .                  *USS(J,1,2)/dsqrt(2.d0)
      CONESSsL(I,J,2)=-mqs/vdq
     .        *(NEU(I,4,1)*USS(J,2,2)-NEU(I,4,2)*USS(J,2,1))
     . -(dsqrt(g1q)/3.d0*NEU(I,1,1)-dsqrt(g2q)*NEU(I,2,1))
     .        *USS(J,1,2)/dsqrt(2.d0)
     . +(dsqrt(g1q)/3.d0*NEU(I,1,2)-dsqrt(g2q)*NEU(I,2,2))
     .                  *USS(J,1,1)/dsqrt(2.d0)
      CONESSsR(I,J,1)=-mqs/vdq
     .        *(NEU(I,4,1)*USS(J,1,1)-NEU(I,4,2)*USS(J,1,2))
     .-dsqrt(2.d0*g1q)/3.d0
     .        *(NEU(I,1,1)*USS(J,2,1)-NEU(I,1,2)*USS(J,2,2))
      CONESSsR(I,J,2)=-mqs/vdq
     .        *(NEU(I,4,2)*USS(J,1,1)+NEU(I,4,1)*USS(J,1,2))
     .-dsqrt(2.d0*g1q)/3.d0
     .        *(NEU(I,1,1)*USS(J,2,2)+NEU(I,1,2)*USS(J,2,1))

      CONESEeL(I,J,1)=-me/vdq*(NEU(I,4,1)*UE(J,2,1)+NEU(I,4,2)
     . *UE(J,2,2))+(dsqrt(g1q)*NEU(I,1,1)+dsqrt(g2q/2.d0)*NEU(I,2,1))
     . *UE(J,1,1)
     . +(dsqrt(g1q)*NEU(I,1,2)+dsqrt(g2q/2.d0)*NEU(I,2,2))*UE(J,1,2)
      CONESEeL(I,J,2)=-me/vdq*(NEU(I,4,1)*UE(J,2,2)-NEU(I,4,2)
     . *UE(J,2,1))+(dsqrt(g1q)*NEU(I,1,1)+dsqrt(g2q/2.d0)*NEU(I,2,1))
     . *UE(J,1,2)
     . -(dsqrt(g1q)*NEU(I,1,2)+dsqrt(g2q/2.d0)*NEU(I,2,2))*UE(J,1,1)
      CONESEeR(I,J,1)=-me/vdq*(NEU(I,4,1)*UE(J,1,1)-NEU(I,4,2)
     . *UE(J,1,2))-dsqrt(2.d0)*dsqrt(g1q)
     . *(NEU(I,1,1)*UE(J,2,1)-NEU(I,1,2)*UE(J,2,2))
      CONESEeR(I,J,2)=-me/vdq*(NEU(I,4,2)*UE(J,1,1)+NEU(I,4,1)
     . *UE(J,1,2))-dsqrt(2.d0)*dsqrt(g1q)
     . *(NEU(I,1,1)*UE(J,2,2)+NEU(I,1,2)*UE(J,2,1))
     
        ENDDO
        ENDDO

c      V- Higgs-Sfermion couplings

      DO I=1,5
       XHG(I,1)=XH(I,1)/dsqrt(ZHU)
       XHG(I,2)=XH(I,2)/dsqrt(ZHD)
       XHG(I,3)=XH(I,3)/dsqrt(ZS)
       XHG(I,4)=XH(I,4)*cosb/dsqrt(ZHU)
       XHG(I,5)=XH(I,4)*sinb/dsqrt(ZHD)
       XHG(I,6)=XH(I,5)/dsqrt(ZS)

      DO M=1,2

       gRHSUSU(I,M,M)=dsqrt(2.d0)
     .               *(USU(M,1,1)*USU(M,1,1)+USU(M,1,2)*USU(M,1,2))
     . *(mqu**2/vuq*XHG(I,1)+(g1q/3.d0-g2q)/4.d0
     .                               *(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +dsqrt(2.d0)*(USU(M,2,1)*USU(M,2,1)+USU(M,2,2)*USU(M,2,2))
     . *(mqu**2/vuq*XHG(I,1)-g1q/3.d0*(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +mqu/vuq/dsqrt(2.d0)
     .   *(AU*(DDCOS(PhiAU)*XHG(I,1)-DDSIN(PhiAT)*XHG(I,4))
     .     -DDCOS(Phi01)*(muq*XHG(I,2)+l*vdq*XHG(I,3))
     .     +DDSIN(Phi01)*(muq*XHG(I,5)+l*vdq*XHG(I,6)))
     .  *(USU(M,1,1)*USU(M,2,1)+USU(M,2,1)*USU(M,1,1)
     .      +USU(M,1,2)*USU(M,2,2)+USU(M,2,2)*USU(M,1,2))
     . +mqu/vuq/dsqrt(2.d0)
     .   *(AU*(DDSIN(PhiAT)*XHG(I,1)+DDCOS(PhiAU)*XHG(I,4))
     .    +DDSIN(Phi01)*(muq*XHG(I,2)+l*vdq*XHG(I,3))
     .    +DDCOS(Phi01)*(muq*XHG(I,5)+l*vdq*XHG(I,6)))
     .  *(USU(M,2,1)*USU(M,1,2)-USU(M,2,2)*USU(M,1,1)
     .      +USU(M,1,2)*USU(M,2,1)-USU(M,1,1)*USU(M,2,2))

       gRHSCSC(I,M,M)=dsqrt(2.d0)
     .               *(USC(M,1,1)*USC(M,1,1)+USC(M,1,2)*USC(M,1,2))
     . *(mc**2/vuq*XHG(I,1)+(g1q/3.d0-g2q)/4.d0
     .                               *(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +dsqrt(2.d0)*(USC(M,2,1)*USC(M,2,1)+USC(M,2,2)*USC(M,2,2))
     . *(mc**2/vuq*XHG(I,1)-g1q/3.d0*(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +mc/vuq/dsqrt(2.d0)
     .   *(AC*(DDCOS(PhiAC)*XHG(I,1)-DDSIN(PhiAT)*XHG(I,4))
     .     -DDCOS(Phi01)*(muq*XHG(I,2)+l*vdq*XHG(I,3))
     .     +DDSIN(Phi01)*(muq*XHG(I,5)+l*vdq*XHG(I,6)))
     .  *(USC(M,1,1)*USC(M,2,1)+USC(M,2,1)*USC(M,1,1)
     .      +USC(M,1,2)*USC(M,2,2)+USC(M,2,2)*USC(M,1,2))
     . +mc/vuq/dsqrt(2.d0)
     .   *(AC*(DDSIN(PhiAT)*XHG(I,1)+DDCOS(PhiAC)*XHG(I,4))
     .    +DDSIN(Phi01)*(muq*XHG(I,2)+l*vdq*XHG(I,3))
     .    +DDCOS(Phi01)*(muq*XHG(I,5)+l*vdq*XHG(I,6)))
     .  *(USC(M,2,1)*USC(M,1,2)-USC(M,2,2)*USC(M,1,1)
     .      +USC(M,1,2)*USC(M,2,1)-USC(M,1,1)*USC(M,2,2))

       gRHSDSD(I,M,M)=dsqrt(2.d0)
     .               *(USD(M,1,1)*USD(M,1,1)+USD(M,1,2)*USD(M,1,2))
     . *(mqd**2/vdq*XHG(I,2)+(g1q/3.d0+g2q)/4.d0
     .                               *(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +dsqrt(2.d0)*(USD(M,2,1)*USD(M,2,1)+USD(M,2,2)*USD(M,2,2))
     . *(mqd**2/vdq*XHG(I,2)+g1q/6.d0*(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +mqd/vdq/dsqrt(2.d0)
     .    *(AD*(DDCOS(PhiAD)*XHG(I,2)-DDSIN(PhiAD)*XHG(I,5))
     .     -DDCOS(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .      +DDSIN(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(USD(M,1,1)*USD(M,2,1)+USD(M,2,1)*USD(M,1,1)
     .      +USD(M,1,2)*USD(M,2,2)+USD(M,2,2)*USD(M,1,2))
     . +mqd/vdq/dsqrt(2.d0)
     .    *(AD*(DDSIN(PhiAD)*XHG(I,2)+DDCOS(PhiAD)*XHG(I,5))
     .   +DDSIN(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .   +DDCOS(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(USD(M,2,1)*USD(M,1,2)-USD(M,2,2)*USD(M,1,1)
     .      +USD(M,1,2)*USD(M,2,1)-USD(M,1,1)*USD(M,2,2))

       gRHSSSS(I,M,M)=dsqrt(2.d0)
     .               *(USS(M,1,1)*USS(M,1,1)+USS(M,1,2)*USS(M,1,2))
     . *(mqs**2/vdq*XHG(I,2)+(g1q/3.d0+g2q)/4.d0
     .                               *(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +dsqrt(2.d0)*(USS(M,2,1)*USS(M,2,1)+USS(M,2,2)*USS(M,2,2))
     . *(mqs**2/vdq*XHG(I,2)+g1q/6.d0*(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +mqs/vdq/dsqrt(2.d0)
     .    *(AS*(DDCOS(PhiAS)*XHG(I,2)-DDSIN(PhiAS)*XHG(I,5))
     .     -DDCOS(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .      +DDSIN(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(USS(M,1,1)*USS(M,2,1)+USS(M,2,1)*USS(M,1,1)
     .      +USS(M,1,2)*USS(M,2,2)+USS(M,2,2)*USS(M,1,2))
     . +mqs/vdq/dsqrt(2.d0)
     .    *(AS*(DDSIN(PhiAS)*XHG(I,2)+DDCOS(PhiAS)*XHG(I,5))
     .   +DDSIN(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .   +DDCOS(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(USS(M,2,1)*USS(M,1,2)-USS(M,2,2)*USS(M,1,1)
     .      +USS(M,1,2)*USS(M,2,1)-USS(M,1,1)*USS(M,2,2))

       gRHSESE(I,M,M)=dsqrt(2.d0)*
     .        (UE(M,1,1)*UE(M,1,1)+UE(M,1,2)*UE(M,1,2))
     . *(me**2/vdq*XHG(I,2)+(-g1q+g2q)/4.d0
     .                               *(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +dsqrt(2.d0)*
     .        (UE(M,2,1)*UE(M,2,1)+UE(M,2,2)*UE(M,2,2))
     . *(me**2/vdq*XHG(I,2)+g1q/2.d0*(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +me/vdq/dsqrt(2.d0)*
     .            (AE*(DDCOS(PhiAE)*XHG(I,2)-DDSIN(PhiAE)*XHG(I,5))
     .   -DDCOS(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .   +DDSIN(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(UE(M,1,1)*UE(M,2,1)+UE(M,2,1)*UE(M,1,1)
     .      +UE(M,1,2)*UE(M,2,2)+UE(M,2,2)*UE(M,1,2))
     . +me/vdq/dsqrt(2.d0)*
     .            (AE*(DDSIN(PhiAE)*XHG(I,2)+DDCOS(PhiAE)*XHG(I,5))
     .   +DDSIN(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .   +DDCOS(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(UE(M,2,1)*UE(M,1,2)-UE(M,2,2)*UE(M,1,1)
     .      +UE(M,1,2)*UE(M,2,1)-UE(M,1,1)*UE(M,2,2))

       gRHSMSM(I,M,M)=dsqrt(2.d0)*
     .        (UMU(M,1,1)*UMU(M,1,1)+UMU(M,1,2)*UMU(M,1,2))
     . *(mmu**2/vdq*XHG(I,2)+(-g1q+g2q)/4.d0
     .                               *(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +dsqrt(2.d0)*
     .        (UMU(M,2,1)*UMU(M,2,1)+UMU(M,2,2)*UMU(M,2,2))
     . *(mmu**2/vdq*XHG(I,2)+g1q/2.d0*(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +mmu/vdq/dsqrt(2.d0)*
     .            (AMU*(DDCOS(PhiAMU)*XHG(I,2)-DDSIN(PhiAMU)*XHG(I,5))
     .   -DDCOS(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .   +DDSIN(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(UMU(M,1,1)*UMU(M,2,1)+UMU(M,2,1)*UMU(M,1,1)
     .      +UMU(M,1,2)*UMU(M,2,2)+UMU(M,2,2)*UMU(M,1,2))
     . +mmu/vdq/dsqrt(2.d0)*
     .            (AMU*(DDSIN(PhiAMU)*XHG(I,2)+DDCOS(PhiAMU)*XHG(I,5))
     .   +DDSIN(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .   +DDCOS(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(UMU(M,2,1)*UMU(M,1,2)-UMU(M,2,2)*UMU(M,1,1)
     .      +UMU(M,1,2)*UMU(M,2,1)-UMU(M,1,1)*UMU(M,2,2))
     
      ENDDO
      ENDDO
      
      ENDIF

c        A: Electron / quark (C)EDMs

c      I- 1-loop SUSY contributions

c       1) Electron EDM

      dEe=0.d0
      errdEe=0.d0

* Chargino/sneutrino
      aux1=0.d0
      DO I=1,2
       aux=dsqrt(MCH2(I))/16.d0/Pi**2/MSNE2*
     .   (COCHSNeR(I,1)*COCHSNeL(I,2)-COCHSNeR(I,2)*COCHSNeL(I,1))
       aux1=aux1+aux*(-AEDM(MCH2(I)/MSNE2))
      ENDDO
      dEe=dEe+aux1
      errdEe=errdEe+0.1d0*dabs(aux1)

* Neutralino/Selectron
      aux1=0.d0
      DO I=1,5
      DO J=1,2
       aux=dsqrt(MNEU(I))/16.d0/Pi**2/MSE2(J)*
     . (CONESEeR(I,J,1)*CONESEeL(I,J,2)-CONESEeR(I,J,2)*CONESEeL(I,J,1))
       aux1=aux1+aux*(-BEDM(MNEU(I)/MSE2(J)))
      ENDDO
      ENDDO
      dEe=dEe+aux1
      errdEe=errdEe+0.1d0*dabs(aux1)

c       2) Up (C)EDM

      dEu=0.d0
      dCu=0.d0
      errdEu=0.d0
      errdCu=0.d0

* Chargino/sdown
      aux1=0.d0
      aux2=0.d0
      DO I=1,2
      DO J=1,2
       aux=dsqrt(MCH2(I))/16.d0/Pi**2/MSD2(J)*
     . (COCHSDuR(I,J,1)*COCHSDuL(I,J,2)-COCHSDuR(I,J,2)*COCHSDuL(I,J,1))
       aux1=aux1+aux*(AEDM(MCH2(I)/MSD2(J))-BEDM(MCH2(I)/MSD2(J))/3.d0)
       aux2=aux2+aux*BEDM(MCH2(I)/MSD2(J))
      ENDDO
      ENDDO
      dEu=dEu+aux1
      dCu=dCu+aux2
      errdEu=0.3d0*dabs(aux1)
      errdCu=0.3d0*dabs(aux2)

* Neutralino/sup
      aux1=0.d0
      aux2=0.d0
      DO I=1,5
      DO J=1,2
       aux=dsqrt(MNEU(I))/16.d0/Pi**2/MSU2(J)*
     . (CONESUuR(I,J,1)*CONESUuL(I,J,2)-CONESUuR(I,J,2)*CONESUuL(I,J,1))
       aux1=aux1+aux*(2.d0*BEDM(MNEU(I)/MSU2(J))/3.d0)
       aux2=aux2+aux*BEDM(MNEU(I)/MSU2(J))
      ENDDO
      ENDDO
       dEu=dEu+aux1
       dCu=dCu+aux2
      errdEu=errdEu+0.3d0*dabs(aux1)
      errdCu=errdCu+0.3d0*dabs(aux2)

* Gluino/sup
      aux1=0.d0
      aux2=0.d0
      DO J=1,2
       aux=MGL/4.d0/Pi/MSU2(J)*ALSQ*
     . (DDCOS(PhiM3)*(USU(J,1,2)*USU(J,2,1)-USU(J,1,1)*USU(J,2,2))
     . +DDSIN(PhiM3)*(USU(J,1,1)*USU(J,2,1)+USU(J,1,2)*USU(J,2,2)))
       aux1=aux1+8.d0/3.d0*aux*(2.d0*BEDM(MGL**2/MSU2(J))/3.d0)
       aux2=aux2+aux*(-CEDM(MGL**2/MSU2(J)))
      ENDDO
       dEu=dEu+aux1
       dCu=dCu+aux2
      errdEu=errdEu+0.3d0*dabs(aux1)
      errdCu=errdCu+0.3d0*dabs(aux2)

c       3) Down (C)EDM

      dEd=0.d0
      dCd=0.d0
      errdEd=0.d0
      errdCd=0.d0

* Chargino/sup
      aux1=0.d0
      aux2=0.d0
      DO I=1,2
      DO J=1,2
       aux=dsqrt(MCH2(I))/16.d0/Pi**2/MSU2(J)*
     . (COCHSUdR(I,J,1)*COCHSUdL(I,J,2)-COCHSUdR(I,J,2)*COCHSUdL(I,J,1))
       aux1=aux1+aux*(-AEDM(MCH2(I)/MSU2(J))
     .                                +2.d0*BEDM(MCH2(I)/MSU2(J))/3.d0)
       aux2=aux2+aux*BEDM(MCH2(I)/MSU2(J))
      ENDDO
      ENDDO
      dEd=dEd+aux1
      dCd=dCd+aux2
      errdEd=errdEd+0.3d0*dabs(aux1)
      errdCd=errdCd+0.3d0*dabs(aux2)

* Neutralino/sdown
      aux1=0.d0
      aux2=0.d0
      DO I=1,5
      DO J=1,2
       aux=dsqrt(MNEU(I))/16.d0/Pi**2/MSD2(J)*
     . (CONESDdR(I,J,1)*CONESDdL(I,J,2)-CONESDdR(I,J,2)*CONESDdL(I,J,1))
       aux1=aux1+aux*(-BEDM(MNEU(I)/MSD2(J))/3.d0)
       aux2=aux2+aux*BEDM(MNEU(I)/MSD2(J))
      ENDDO
      ENDDO
       dEd=dEd+aux1
       dCd=dCd+aux2
      errdEd=errdEd+0.3d0*dabs(aux1)
      errdCd=errdCd+0.3d0*dabs(aux2)
* Gluino/sdown
      aux1=0.d0
      aux2=0.d0
      DO J=1,2
       aux=MGL/4.d0/Pi/MSD2(J)*ALSQ*
     . (DDCOS(PhiM3)*(USD(J,1,2)*USD(J,2,1)-USD(J,1,1)*USD(J,2,2))
     . +DDSIN(PhiM3)*(USD(J,1,1)*USD(J,2,1)+USD(J,1,2)*USD(J,2,2)))
       aux1=aux1+aux*8.d0/3.d0*(-BEDM(MGL**2/MSD2(J))/3.d0)
       aux2=aux2+aux*(-CEDM(MGL**2/MSD2(J)))
      ENDDO
       dEd=dEd+aux1
       dCd=dCd+aux2
      errdEd=errdEd+0.3d0*dabs(aux1)
      errdCd=errdCd+0.3d0*dabs(aux2)

c       4) Strange EDM

      dEs=0.d0
      errdEs=0.d0

* Chargino/scharm
      aux1=0.d0
      DO I=1,2
      DO J=1,2
       aux=dsqrt(MCH2(I))/16.d0/Pi**2/MSC2(J)*
     . (COCHSCsR(I,J,1)*COCHSCsL(I,J,2)-COCHSCsR(I,J,2)*COCHSCsL(I,J,1))
       aux1=aux1+aux*(-AEDM(MCH2(I)/MSC2(J))
     .                                +2.d0*BEDM(MCH2(I)/MSC2(J))/3.d0)
      ENDDO
      ENDDO
      dEs=dEs+aux1
      errdEs=0.3d0*dabs(aux1)

* Neutralino/sstrange
      aux1=0.d0
      DO I=1,5
      DO J=1,2
       aux=dsqrt(MNEU(I))/16.d0/Pi**2/MSS2(J)*
     . (CONESSsR(I,J,1)*CONESSsL(I,J,2)-CONESSsR(I,J,2)*CONESSsL(I,J,1))
       aux1=aux1+aux*(-BEDM(MNEU(I)/MSS2(J))/3.d0)
      ENDDO
      ENDDO
       dEs=dEs+aux1
      errdEs=errdEs+0.3d0*dabs(aux1)

* Gluino/sstrange
      aux1=0.d0
      DO J=1,2
       aux=MGL/4.d0/Pi/MSS2(J)*ALSQ*
     . (DDCOS(PhiM3)*(USS(J,1,2)*USS(J,2,1)-USS(J,1,1)*USS(J,2,2))
     . +DDSIN(PhiM3)*(USS(J,1,1)*USS(J,2,1)+USS(J,1,2)*USS(J,2,2)))
       aux1=aux1+aux*8.d0/3.d0*(-BEDM(MGL**2/MSS2(J))/3.d0)
      ENDDO
      dEs=dEs+aux1
      errdEs=errdEs+0.3d0*dabs(aux1)

c      II- 2-loop Bar-Zee diagrams

c       1) Electron

      aux=0.d0

      DO I=1,5
      DO J=1,2
* Stop
       aux1=CDP(I)*(gRHSTST(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                         *2.d0*fSF(MST2(J)/MH0(I))/MST2(J)
* Sup/Scharm
       aux2=CDP(I)*(gRHSUSU(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                         *2.d0*fSF(MSU2(J)/MH0(I))/MSU2(J)
     .       +CDP(I)*(gRHSCSC(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                         *2.d0*fSF(MSC2(J)/MH0(I))/MSC2(J)
       aux=aux+3.d0*(2.d0/3.d0)**2*(aux1+aux2)
* Sbottom
       aux1=CDP(I)*(gRHSBSB(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                          *2.d0*fSF(MSB2(J)/MH0(I))/MSB2(J)
* Sdown/Sstrange
       aux2=CDP(I)*(gRHSDSD(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                          *2.d0*fSF(MSD2(J)/MH0(I))/MSD2(J)
     .       +CDP(I)*(gRHSSSS(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                          *2.d0*fSF(MSS2(J)/MH0(I))/MSS2(J)
       aux=aux+3.d0*(-1.d0/3.d0)**2*(aux1+aux2)
* Stau
       aux1=CDP(I)*(gRHSLSL(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                          *2.d0*fSF(MSL2(J)/MH0(I))/MSL2(J)
* Selectron/Smuon
       aux2=CDP(I)*(gRHSESE(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                          *2.d0*fSF(MSE2(J)/MH0(I))/MSE2(J)
     .       +CDP(I)*(gRHSMSM(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                          *2.d0*fSF(MSMU2(J)/MH0(I))/MSMU2(J)

       aux=aux+(-1.d0)**2*(aux1+aux2)
      ENDDO
* top
       aux1=CDP(I)*CU(I)*(-2.d0*fS(mt**2/MH0(I)))
       aux2=CD(I)*CUP(I)*(2.d0*fPS(mt**2/MH0(I)))
       aux=aux+3.d0*(2.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* charm
       aux1=CDP(I)*CU(I)*(-2.d0*fS(mc**2/MH0(I)))
       aux2=CD(I)*CUP(I)*(2.d0*fPS(mc**2/MH0(I)))
       aux=aux+3.d0*(2.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* up
       aux1=CDP(I)*CU(I)*(-2.d0*fS(mqu**2/MH0(I)))
       aux2=CD(I)*CUP(I)*(2.d0*fPS(mqu**2/MH0(I)))
       aux=aux+3.d0*(2.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* bottom
       aux1=CDP(I)*CB(I)*(-2.d0*fS(mb**2/MH0(I)))
       aux2=CD(I)*CBP(I)*(2.d0*fPS(mb**2/MH0(I)))
       aux=aux+3.d0*(-1.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* strange
       aux1=CDP(I)*CD(I)*(-2.d0*fS(mqs**2/MH0(I)))
       aux2=CD(I)*CDP(I)*(2.d0*fPS(mqs**2/MH0(I)))
       aux=aux+3.d0*(-1.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* down
       aux1=CDP(I)*CD(I)*(-2.d0*fS(mqd**2/MH0(I)))
       aux2=CD(I)*CDP(I)*(2.d0*fPS(mqd**2/MH0(I)))
       aux=aux+3.d0*(-1.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* tau
       aux1=CDP(I)*CL(I)*(-2.d0*fS(mtau**2/MH0(I)))
       aux2=CD(I)*CLP(I)*(2.d0*fPS(mtau**2/MH0(I)))
       aux=aux+(-1.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* mu
       aux1=CDP(I)*CD(I)*(-2.d0*fS(mmu**2/MH0(I)))
       aux2=CD(I)*CDP(I)*(2.d0*fPS(mmu**2/MH0(I)))
       aux=aux+(-1.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* e
       aux1=CDP(I)*CD(I)*(-2.d0*fS(me**2/MH0(I)))
       aux2=CD(I)*CDP(I)*(2.d0*fPS(me**2/MH0(I)))
       aux=aux+(-1.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* chargino
      DO J=1,2
       aux1=CDP(I)*COH0CH(I,J,J,1)/dsqrt(2.d0*g2)
     .                                *(-2.d0*fS(MCH2(J)/MH0(I)))
       aux2=CD(I)*COH0CH(I,J,J,2)/dsqrt(2.d0*g2)
     .                                *(2.d0*fPS(MCH2(J)/MH0(I)))
       aux=aux+4.d0*dsqrt(2.d0)*Pi*ALEMMZ/S2TW/MW/dsqrt(MCH2(J))
     .                                                  *(aux1+aux2)
      ENDDO
      ENDDO

      dEe=dEe-ALEMMZ*(-1.d0)*me/32.d0/Pi**3*aux
      errdEe=errdEe+0.3d0*dabs(ALEMMZ*me/32.d0/Pi**3*aux)

c       2) Up (C)EDM

      aux=0.d0

      DO I=1,5
      DO J=1,2
* Stop
       aux1=CUP(I)*(gRHSTST(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MST2(J)/MH0(I))/MST2(J)
* Sup/Scharm
       aux2=CUP(I)*(gRHSUSU(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSU2(J)/MH0(I))/MSU2(J)
     .       +CUP(I)*(gRHSCSC(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSC2(J)/MH0(I))/MSC2(J)
       aux=aux+3.d0*(2.d0/3.d0)**2*(aux1+aux2)
* Sbottom
       aux1=CUP(I)*(gRHSBSB(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSB2(J)/MH0(I))/MSB2(J)
* Sdown/Sstrange
       aux2=CUP(I)*(gRHSDSD(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSD2(J)/MH0(I))/MSD2(J)
     .       +CUP(I)*(gRHSSSS(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSS2(J)/MH0(I))/MSS2(J)
       aux=aux+3.d0*(-1.d0/3.d0)**2*(aux1+aux2)
* Stau
       aux1=CUP(I)*(gRHSLSL(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSL2(J)/MH0(I))/MSL2(J)
* Selectron/Smuon
       aux2=CUP(I)*(gRHSESE(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSE2(J)/MH0(I))/MSE2(J)
     .       +CUP(I)*(gRHSMSM(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSMU2(J)/MH0(I))/MSMU2(J)
       aux=aux+(-1.d0)**2*(aux1+aux2)
      ENDDO
* top
       aux1=CUP(I)*CU(I)*(-2.d0*fS(mt**2/MH0(I)))
       aux2=CU(I)*CUP(I)*(2.d0*fPS(mt**2/MH0(I)))
       aux=aux+3.d0*(2.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* charm
       aux1=CUP(I)*CU(I)*(-2.d0*fS(mc**2/MH0(I)))
       aux2=CU(I)*CUP(I)*(2.d0*fPS(mc**2/MH0(I)))
       aux=aux+3.d0*(2.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* up
       aux1=CUP(I)*CU(I)*(-2.d0*fS(mqu**2/MH0(I)))
       aux2=CU(I)*CUP(I)*(2.d0*fPS(mqu**2/MH0(I)))
       aux=aux+3.d0*(2.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* bottom
       aux1=CUP(I)*CB(I)*(-2.d0*fS(mb**2/MH0(I)))
       aux2=CU(I)*CBP(I)*(2.d0*fPS(mb**2/MH0(I)))
       aux=aux+3.d0*(-1.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* strange
       aux1=CUP(I)*CD(I)*(-2.d0*fS(mqs**2/MH0(I)))
       aux2=CU(I)*CDP(I)*(2.d0*fPS(mqs**2/MH0(I)))
       aux=aux+3.d0*(-1.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* down
       aux1=CUP(I)*CD(I)*(-2.d0*fS(mqd**2/MH0(I)))
       aux2=CU(I)*CDP(I)*(2.d0*fPS(mqd**2/MH0(I)))
       aux=aux+3.d0*(-1.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* tau
       aux1=CUP(I)*CL(I)*(-2.d0*fS(mtau**2/MH0(I)))
       aux2=CU(I)*CLP(I)*(2.d0*fPS(mtau**2/MH0(I)))
       aux=aux+(-1.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* mu
       aux1=CUP(I)*CD(I)*(-2.d0*fS(mmu**2/MH0(I)))
       aux2=CU(I)*CDP(I)*(2.d0*fPS(mmu**2/MH0(I)))
       aux=aux+(-1.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* e
       aux1=CUP(I)*CD(I)*(-2.d0*fS(me**2/MH0(I)))
       aux2=CU(I)*CDP(I)*(2.d0*fPS(me**2/MH0(I)))
       aux=aux+(-1.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* chargino
      DO J=1,2
       aux1=CUP(I)*COH0CH(I,J,J,1)/dsqrt(2.d0*g2)
     .                                *(-2.d0*fS(MCH2(J)/MH0(I)))
       aux2=CU(I)*COH0CH(I,J,J,2)/dsqrt(2.d0*g2)
     .                                *(2.d0*fPS(MCH2(J)/MH0(I)))
       aux=aux+4.d0*dsqrt(2.d0)*Pi*ALEMMZ/S2TW/MW/dsqrt(MCH2(J))
     .                                                  *(aux1+aux2)
      ENDDO
      ENDDO

      dEu=dEu-ALEMMZ*(2.d0/3.d0)*mqu/32.d0/Pi**3*aux
      errdEu=errdEu
     .            +0.3d0*dabs(ALEMMZ*(2.d0/3.d0)*mqu/32.d0/Pi**3*aux)

      aux=0.d0

      DO I=1,5
      DO J=1,2
* Stop
       aux1=CUP(I)*(gRHSTST(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                          *2.d0*fSF(MST2(J)/MH0(I))/MST2(J)
* Sup/Scharm
       aux2=CUP(I)*(gRHSUSU(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                          *2.d0*fSF(MSU2(J)/MH0(I))/MSU2(J)
     .       +CUP(I)*(gRHSCSC(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                          *2.d0*fSF(MSC2(J)/MH0(I))/MSC2(J)
       aux=aux+(aux1+aux2)
* Sbottom
       aux1=CUP(I)*(gRHSBSB(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                          *2.d0*fSF(MSB2(J)/MH0(I))/MSB2(J)
* Sdown/Sstrange
       aux2=CUP(I)*(gRHSDSD(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                          *2.d0*fSF(MSD2(J)/MH0(I))/MSD2(J)
     .       +CUP(I)*(gRHSSSS(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                          *2.d0*fSF(MSS2(J)/MH0(I))/MSS2(J)
       aux=aux+(aux1+aux2)
      ENDDO
* top
       aux1=CUP(I)*CU(I)*(-2.d0*fS(mt**2/MH0(I)))
       aux2=CU(I)*CUP(I)*(2.d0*fPS(mt**2/MH0(I)))
       aux=aux+4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* charm
       aux1=CUP(I)*CU(I)*(-2.d0*fS(mc**2/MH0(I)))
       aux2=CU(I)*CUP(I)*(2.d0*fPS(mc**2/MH0(I)))
       aux=aux+4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* up
       aux1=CUP(I)*CU(I)*(-2.d0*fS(mqu**2/MH0(I)))
       aux2=CU(I)*CUP(I)*(2.d0*fPS(mqu**2/MH0(I)))
       aux=aux+4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* bottom
       aux1=CUP(I)*CB(I)*(-2.d0*fS(mb**2/MH0(I)))
       aux2=CU(I)*CBP(I)*(2.d0*fPS(mb**2/MH0(I)))
       aux=aux+4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* strange
       aux1=CUP(I)*CD(I)*(-2.d0*fS(mqs**2/MH0(I)))
       aux2=CU(I)*CDP(I)*(2.d0*fPS(mqs**2/MH0(I)))
       aux=aux+4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* down
       aux1=CUP(I)*CD(I)*(-2.d0*fS(mqd**2/MH0(I)))
       aux2=CU(I)*CDP(I)*(2.d0*fPS(mqd**2/MH0(I)))
       aux=aux+4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
      ENDDO

      dCu=dCu-alSMZ*mqu/64.d0/Pi**3*aux
      errdCu=errdCu+0.3d0*dabs(alSMZ*mqu/64.d0/Pi**3*aux)

c       3) Down / strange (C)EDM

      aux=0.d0

      DO I=1,5
      DO J=1,2
* Stop
       aux1=CDP(I)*(gRHSTST(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MST2(J)/MH0(I))/MST2(J)
* Sup/Scharm
       aux2=CDP(I)*(gRHSUSU(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSU2(J)/MH0(I))/MSU2(J)
     .       +CDP(I)*(gRHSCSC(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSC2(J)/MH0(I))/MSC2(J)
       aux=aux+3.d0*(2.d0/3.d0)**2*(aux1+aux2)
* Sbottom
       aux1=CDP(I)*(gRHSBSB(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSB2(J)/MH0(I))/MSB2(J)
* Sdown/Sstrange
       aux2=CDP(I)*(gRHSDSD(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSD2(J)/MH0(I))/MSD2(J)
     .       +CDP(I)*(gRHSSSS(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSS2(J)/MH0(I))/MSS2(J)
       aux=aux+3.d0*(-1.d0/3.d0)**2*(aux1+aux2)
* Stau
       aux1=CDP(I)*(gRHSLSL(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSL2(J)/MH0(I))/MSL2(J)
* Selectron/Smuon
       aux2=CDP(I)*(gRHSESE(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSE2(J)/MH0(I))/MSE2(J)
     .       +CDP(I)*(gRHSMSM(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSMU2(J)/MH0(I))/MSMU2(J)
       aux=aux+(-1.d0)**2*(aux1+aux2)
      ENDDO
* top
       aux1=CDP(I)*CU(I)*(-2.d0*fS(mt**2/MH0(I)))
       aux2=CD(I)*CUP(I)*(2.d0*fPS(mt**2/MH0(I)))
       aux=aux+3.d0*(2.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* charm
       aux1=CDP(I)*CU(I)*(-2.d0*fS(mc**2/MH0(I)))
       aux2=CD(I)*CUP(I)*(2.d0*fPS(mc**2/MH0(I)))
       aux=aux+3.d0*(2.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* up
       aux1=CDP(I)*CU(I)*(-2.d0*fS(mqu**2/MH0(I)))
       aux2=CD(I)*CUP(I)*(2.d0*fPS(mqu**2/MH0(I)))
       aux=aux+3.d0*(2.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* bottom
       aux1=CDP(I)*CB(I)*(-2.d0*fS(mb**2/MH0(I)))
       aux2=CD(I)*CBP(I)*(2.d0*fPS(mb**2/MH0(I)))
       aux=aux+3.d0*(-1.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* strange
       aux1=CDP(I)*CD(I)*(-2.d0*fS(mqs**2/MH0(I)))
       aux2=CD(I)*CDP(I)*(2.d0*fPS(mqs**2/MH0(I)))
       aux=aux+3.d0*(-1.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* down
       aux1=CDP(I)*CD(I)*(-2.d0*fS(mqd**2/MH0(I)))
       aux2=CD(I)*CDP(I)*(2.d0*fPS(mqd**2/MH0(I)))
       aux=aux+3.d0*(-1.d0/3.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2
     .                                                  *(aux1+aux2)
* tau
       aux1=CDP(I)*CL(I)*(-2.d0*fS(mtau**2/MH0(I)))
       aux2=CD(I)*CLP(I)*(2.d0*fPS(mtau**2/MH0(I)))
       aux=aux+(-1.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* mu
       aux1=CDP(I)*CD(I)*(-2.d0*fS(mmu**2/MH0(I)))
       aux2=CD(I)*CDP(I)*(2.d0*fPS(mmu**2/MH0(I)))
       aux=aux+(-1.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* e
       aux1=CDP(I)*CD(I)*(-2.d0*fS(me**2/MH0(I)))
       aux2=CD(I)*CDP(I)*(2.d0*fPS(me**2/MH0(I)))
       aux=aux+(-1.d0)**2*4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)

* chargino
      DO J=1,2
       aux1=CDP(I)*COH0CH(I,J,J,1)/dsqrt(2.d0*g2)
     .                                *(-2.d0*fS(MCH2(J)/MH0(I)))
       aux2=CD(I)*COH0CH(I,J,J,2)/dsqrt(2.d0*g2)
     .                                *(2.d0*fPS(MCH2(J)/MH0(I)))
       aux=aux+4.d0*dsqrt(2.d0)*Pi*ALEMMZ/S2TW/MW/dsqrt(MCH2(J))
     .                                                  *(aux1+aux2)
      ENDDO
      ENDDO

      dEd=dEd-ALEMMZ*(-1.d0/3.d0)*mqd/32.d0/Pi**3*aux
      dEs=dEs-ALEMMZ*(-1.d0/3.d0)*mqs/32.d0/Pi**3*aux
      errdEd=errdEd+0.3d0*dabs(ALEMMZ/3.d0*mqd/32.d0/Pi**3*aux)
      errdEs=errdEs+0.3d0*dabs(ALEMMZ/3.d0*mqs/32.d0/Pi**3*aux)

      aux=0.d0

      DO I=1,5
      DO J=1,2
* Stop
       aux1=CDP(I)*(gRHSTST(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MST2(J)/MH0(I))/MST2(J)
* Sup/Scharm
       aux2=CDP(I)*(gRHSUSU(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSU2(J)/MH0(I))/MSU2(J)
     .       +CDP(I)*(gRHSCSC(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSC2(J)/MH0(I))/MSC2(J)
       aux=aux+(aux1+aux2)
* Sbottom
       aux1=CDP(I)*(gRHSBSB(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSB2(J)/MH0(I))/MSB2(J)
* Sdown/Sstrange
       aux2=CDP(I)*(gRHSDSD(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSD2(J)/MH0(I))/MSD2(J)
     .       +CDP(I)*(gRHSSSS(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))
     .                            *2.d0*fSF(MSS2(J)/MH0(I))/MSS2(J)
       aux=aux+(aux1+aux2)
      ENDDO
* top
       aux1=CDP(I)*CU(I)*(-2.d0*fS(mt**2/MH0(I)))
       aux2=CD(I)*CUP(I)*(2.d0*fPS(mt**2/MH0(I)))
       aux=aux+4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* charm
       aux1=CDP(I)*CU(I)*(-2.d0*fS(mc**2/MH0(I)))
       aux2=CD(I)*CUP(I)*(2.d0*fPS(mc**2/MH0(I)))
       aux=aux+4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* up
       aux1=CDP(I)*CU(I)*(-2.d0*fS(mqu**2/MH0(I)))
       aux2=CD(I)*CUP(I)*(2.d0*fPS(mqu**2/MH0(I)))
       aux=aux+4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* bottom
       aux1=CDP(I)*CB(I)*(-2.d0*fS(mb**2/MH0(I)))
       aux2=CD(I)*CBP(I)*(2.d0*fPS(mb**2/MH0(I)))
       aux=aux+4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* strange
       aux1=CDP(I)*CD(I)*(-2.d0*fS(mqs**2/MH0(I)))
       aux2=CD(I)*CDP(I)*(2.d0*fPS(mqs**2/MH0(I)))
       aux=aux+4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
* down
       aux1=CDP(I)*CD(I)*(-2.d0*fS(mqd**2/MH0(I)))
       aux2=CD(I)*CDP(I)*(2.d0*fPS(mqd**2/MH0(I)))
       aux=aux+4.d0*Pi*ALEMMZ/S2TW/MW**2*(aux1+aux2)
      ENDDO

      dCd=dCd-alSMZ*mqd/64.d0/Pi**3*aux
      errdCd=errdCd+0.3d0*dabs(alSMZ*mqd/64.d0/Pi**3*aux)

c       4) WW contribution

      aux=0.d0
      DO I=1,5
      DO J=1,2
       aux1=NEU(I,4,2)*U(J,2,1)-NEU(I,4,1)*U(J,2,2)
     .      +dsqrt(2.d0)*(NEU(I,2,2)*U(J,1,1)-NEU(I,2,1)*U(J,1,2))
       aux2=-(NEU(I,3,1)*V(J,2,1)+NEU(I,3,2)*V(J,2,2))
     .      +dsqrt(2.d0)*(NEU(I,2,1)*V(J,1,1)+NEU(I,2,2)*V(J,1,2))
       aux=aux+dsqrt(MCH2(J)*MNEU(I))*aux1*aux2
     .               *fWW(MNEU(I)/MW**2,MCH2(J)/MW**2)
       aux1=NEU(I,4,1)*U(J,2,1)+NEU(I,4,2)*U(J,2,2)
     .      +dsqrt(2.d0)*(NEU(I,2,1)*U(J,1,1)+NEU(I,2,2)*U(J,1,2))
       aux2=-(NEU(I,3,2)*V(J,2,1)-NEU(I,3,1)*V(J,2,2))
     .      +dsqrt(2.d0)*(NEU(I,2,2)*V(J,1,1)-NEU(I,2,1)*V(J,1,2))
       aux=aux+dsqrt(MCH2(J)*MNEU(I))*aux1*aux2
     .               *fWW(MNEU(I)/MW**2,MCH2(J)/MW**2)
      ENDDO
      ENDDO

c      dEe=dEe+ALEMMZ**2*me/(32.d0*Pi**2*S2TW**2*MW**4)*aux
c      dEu=dEu+ALEMMZ**2*mqu/(32.d0*Pi**2*S2TW**2*MW**4)*aux
c      dEd=dEd+ALEMMZ**2*mqd/(32.d0*Pi**2*S2TW**2*MW**4)*aux
c      dEs=dEs+ALEMMZ**2*mqs/(32.d0*Pi**2*S2TW**2*MW**4)*aux
c      errdEe=errdEe
c     c     +0.1d0*dabs(ALEMMZ**2*me/(32.d0*Pi**2*S2TW**2*MW**4)*aux)
c      errdEu=errdEu
c     c     +0.1d0*dabs(ALEMMZ**2*mqu/(32.d0*Pi**2*S2TW**2*MW**4)*aux)
c      errdEd=errdEd
c     c     +0.1d0*dabs(ALEMMZ**2*mqd/(32.d0*Pi**2*S2TW**2*MW**4)*aux)
c      errdEs=errdEs
c     c     +0.1d0*dabs(ALEMMZ**2*mqs/(32.d0*Pi**2*S2TW**2*MW**4)*aux)

c        B: Other operators

c      I- 4-fermion interactions

      aux=0.d0
      DO I=1,5
      aux=aux+CD(I)*CDP(I)/MH0(I)
      ENDDO
      Cee=g2/4d0*(me/MW)**2*aux
      Cdd=g2/4d0*(mqd/MW)**2*aux
      Ced=g2/4d0*mqd*me/MW**2*aux
      Cde=Ced
      Ces=g2/4d0*mqs*me/MW**2*aux
      Cse=Ces
      Ceb=g2/4d0*me*mb/MW**2*aux
      Cbd=g2/4d0*mqd*mb/MW**2*aux
      Cdb=Cbd
      Csd=g2/4d0*mqd*mqs/MW**2*aux

      aux=0.d0
      DO I=1,5
      aux=aux+CU(I)*CUP(I)/MH0(I)
      ENDDO
      Cuu=g2/4d0*(mqu/MW)**2*aux

      aux=0.d0
      DO I=1,5
      aux=aux+CU(I)*CDP(I)/MH0(I)
      ENDDO
      Cue=g2/4d0*mqu*me/MW**2*aux
      Cud=g2/4d0*mqu*mqd/MW**2*aux

      aux=0.d0
      DO I=1,5
      aux=aux+CD(I)*CUP(I)/MH0(I)
      ENDDO
      Ceu=g2/4d0*mqu*me/MW**2*aux
      Cdu=g2/4d0*mqu*mqd/MW**2*aux
      Cec=g2/4d0*mc*me/MW**2*aux
      Cet=g2/4d0*mt*me/MW**2*aux

c      II- Weinberg operator

* Higgs/quark contribution
      aux=0.d0
      DO I=1,5
       aux=aux+CU(I)*CUP(I)*HEDM(MH0(I)/mt**2)
     .          +CD(I)*CDP(I)*HEDM(MH0(I)/mb**2)
      ENDDO
      dGW=4d0*dsqrt(2d0)*GF*alsMZ*dsqrt(alsMZ/4d0/Pi)/(4d0*Pi)**2
     .              *aux
      errdGW=0.3d0*dabs(dGW)

* gluino/quark/squark contribution
      aux=0.d0
      DO I=1,2
      aux=aux+2d0*Pi*alsMZ*mt*MST2(I)/MGL**2*(
     .  DDSIN(PhiM3)*(UT(I,1,1)*UT(I,2,1)+UT(I,1,2)*UT(I,2,2))
     . -DDCOS(PhiM3)*(UT(I,1,1)*UT(I,2,2)-UT(I,1,2)*UT(I,2,1)))
      ENDDO
      aux1=(MST2(1)+MST2(2))/2d0/MGL**2
      aux2=(MST2(2)-MST2(1))/(MST2(1)+MST2(2))
      dGW=dGW-3d0/(2d0*Pi)*dsqrt(alsMZ/4d0/Pi)**3/MGL**3*aux
     .    *HQEDM((mt/MGL)**2,aux1,aux2)
      errdGW=errdGW+0.3d0*dabs(3d0/(2d0*Pi)*dsqrt(alsMZ/4d0/Pi)**3
     .      /MGL**3*aux*HQEDM((mt/MGL)**2,aux1,aux2))

      aux=0.d0
      DO I=1,2
      aux=aux+2d0*Pi*alsMZ*mb*MSB2(I)/MGL**2*(
     .  DDSIN(PhiM3)*(UB(I,1,1)*UB(I,2,1)+UB(I,1,2)*UB(I,2,2))
     . -DDCOS(PhiM3)*(UB(I,1,1)*UB(I,2,2)-UB(I,1,2)*UB(I,2,1)))
      ENDDO
      aux1=(MSB2(1)+MSB2(2))/2d0/MGL**2
      aux2=(MSB2(2)-MSB2(1))/(MSB2(1)+MSB2(2))
      dGW=dGW-3d0/(2d0*Pi)*dsqrt(alsMZ/4d0/Pi)**3/MGL**3*aux
     .    *HQEDM((mb/MGL)**2,aux1,aux2)
      errdGW=errdGW+0.3d0*dabs(3d0/(2d0*Pi)*dsqrt(alsMZ/4d0/Pi)**3
     .      /MGL**3*aux*HQEDM((mb/MGL)**2,aux1,aux2))

c      III- CS Operator

* Gluon-gluon-Higgs contribution
      aux=0d0
      DO I=1,5
       aux1=2d0/3d0*(CU(I)+(1-0.25d0*0.5d0)*CD(I))
       DO J=1,2
        aux1=aux1-(vu**2+vd**2)/6d0
     .  *(gRHSTST(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))/MST2(J)
     .             -(vu**2+vd**2)/6d0
     .  *(gRHSBSB(I,J,J)/dsqrt(2.d0*(vu**2+vd**2)))/MSB2(J)
       ENDDO
       aux=aux+aux1*CDP(I)/MH0(I)
      ENDDO
      CSG=0.1d0*me/2.d0/(vu**2+vd**2)*aux
      errCSG=0.5d0*dabs(CSG)

* 4-fermion contribution
      CSG=CSG+0.029d0/mqd*Cde+0.5d0*0.22d0/mqs*Cse
      errCSG=errCSG+0.5d0*dabs(0.029d0/mqd*Cde+0.5d0*0.22d0/mqs*Cse)

c       C: Observable EDMs

c      I- Thalium EDM

      dETl=-585d0*dEe-43d0*CSG
      aux=585d0*errdEe+43d0*errCSG
      
      dETlmin=min(dabs(dETl-aux),dabs(dETl+aux))
      IF(dabs(dETl).le.dabs(aux))dETlmin=0d0

      dETl=dETl*1.97d-14                                ! conversion to e.cm
      dETlmin=dETlmin*1.97d-14                          ! conversion to e.cm

*      Comparison with experiment: |dETl| < 1.3d-24 e.cm ; cf. Phys. Rev. Lett. 88 (2002) 071805.
      IF(dETlmin.ge.1.3d-24)PROB(54)=dETlmin/1.3d-24-1d0

c      II- Neutron EDM

c       1) Chiral model

      dEnCQM=4d0/3d0*(1.53d0*dEd+3.4d0*(dCd*dsqrt(ALSMZ/4d0/Pi)
     .                               +1.19d0/4d0/Pi*dGW))
     .        -1d0/3d0*(1.53d0*dEu+3.4d0*(dCu*dsqrt(ALSMZ/4d0/Pi)
     .                               +1.19d0/4d0/Pi*dGW))
      aux=4d0/3d0*(1.53d0*errdEd+3.4d0*(errdCd*dsqrt(ALSMZ/4d0/Pi)
     .                               +1.19d0/4d0/Pi*errdGW))
     .     +1d0/3d0*(1.53d0*errdEu+3.4d0*(errdCu*dsqrt(ALSMZ/4d0/Pi)
     .                               +1.19d0/4d0/Pi*errdGW))
     .     +0.2d0*dabs(dEnCQM)

c       2) Parton model

      dEnPQM=1.53d0*(0.746d0*dEd-0.508d0*dEu-0.226d0*dEs)
      aux1=1.53d0*(0.746d0*errdEd+0.508d0*errdEu+0.226d0*errdEs)
     . +0.2d0*dabs(dEnPQM)

c       3) QCD sum rules

      dEnQSR=1.4d0*(dEd-0.25d0*dEu)
     .        +1.1d0*(dCd+0.5d0*dCu)
     .        +0.02d0*8.5d0*dGW
     .        +2.6d-3*(Cbd/mb+0.75d0*Cdb/mb)
      aux2=1.4d0*(errdEd+0.25d0*errdEu)+1.1d0*(errdCd+0.5d0*errdCu)
     . +0.02d0*8.5d0*errdGW+2.6d-3*0.3d0*dabs(Cbd/mb+0.75d0*Cdb
     . /mb)+0.3d0*dabs(dEnQSR)

c       4) Summary
      dEnCQMmin=min(dabs(dEnCQM-aux),dabs(dEnCQM+aux))
      IF(dabs(dEnCQM).le.dabs(aux))dEnCQMmin=0d0
      dEnPQMmin=min(dabs(dEnPQM-aux1),dabs(dEnPQM+aux1))
      IF(dabs(dEnPQM).le.dabs(aux1))dEnPQMmin=0d0
      dEnQSRmin=min(dabs(dEnQSR-aux2),dabs(dEnQSR+aux2))
      IF(dabs(dEnQSR).le.dabs(aux2))dEnQSRmin=0d0
      
      dEnCQM=dEnCQM*1.97d-14                                ! conversion to e.cm
      dEnCQMmin=dEnCQMmin*1.97d-14                          ! conversion to e.cm
      
      dEnPQM=dEnPQM*1.97d-14                                ! conversion to e.cm
      dEnPQMmin=dEnPQMmin*1.97d-14                          ! conversion to e.cm
      
      dEnQSR=dEnQSR*1.97d-14                                ! conversion to e.cm
      dEnQSRmin=dEnQSRmin*1.97d-14                          ! conversion to e.cm
      
*      Comparison with experiment: |dEn| < 3d-26 e.cm ; cf. arXiv:hep-ex/0602020
      aux=min(dEnCQMmin,dEnPQMmin,dEnQSRmin)
      IF(aux.ge.3d-26)PROB(54)=PROB(54)+aux/3d-26-1d0


c      III- Mercury EDM (QCD sum rules)

      aux=-0.375d0*(Ces/mqs+Cec/mc+Ceb/mb+Cet/mt)-(-0.2d0)*
     . (0.806d0*Ced/mqd+0.181d0*(Cec/mc+Ces/mqs+Ceb/mb+Cet/mt))
      dEHg=0.007d0*(dCu-dCd)+0.01d0*dEe
     . -1.4d-5*(0.5d0*Cdd/mqd+3.3d0*0.5d0*Csd/mqs
     .                                   +(1d0-0.25d0*0.5d0)*Cbd/mb)
     . +0.0035d0*CSG
     . +4d-4*aux
      aux1=0.007d0*(errdCu+errdCd)+0.01d0*errdEe+1.4d-5*0.3d0*dabs(
     . 0.5d0*Cdd/mqd+3.3d0*0.5d0*Csd/mqs+(1d0-0.25d0*0.5d0)*Cbd/mb)
     . +0.0035d0*errCSG
     . +4d-4*0.5d0*dabs(aux)

      dEHgmin=min(dabs(dEHg-aux1),dabs(dEHg+aux1))
      IF(dabs(dEHg).le.dabs(aux1))dEHgmin=0d0
      dEHg=dEHg*1.97d-14                                ! conversion to e.cm
      dEHgmin=dEHgmin*1.97d-14                          ! conversion to e.cm

*      Comparison with experiment: |dEHg| < 3.1d-29 e.cm ; cf. Phys. Rev. Lett. 102 (2009), 101601
      IF(dEHgmin.ge.3.1d-29)PROB(54)=PROB(54)+dEHgmin/3.1d-29-1d0


c      IV- Electron EDM
      dEemin=min(dabs(dEe-errdEe),dabs(dEe+errdEe))
      IF(dabs(dEe).le.dabs(errdEe))dEemin=0d0
      
      dEe=dEe*1.97d-14                                ! conversion to e.cm
      dEemin=dEemin*1.97d-14                          ! conversion to e.cm
      
*      Comparison with experiment: |dEe| < 1d-28 e.cm ; cf. arXiv:1310.7534
      IF(dEemin.ge.1d-28)PROB(54)=PROB(54)+dEemin/1d-28-1d0

c      print*,'here',dEe,dEu,dEd,dEs,dCu,dCd,CSG,dGW
c      print*,errdEe,errdEu,errdEd,errdEs,errdCu,errdCd,errCSG,errdGW
c      print*,PROB(54),dETl,dEnCQM,dEnPQM,dEnQSR,dEHg
c      print*,dETlmin,dEnCQMmin,dEnPQMmin,dEnQSRmin,dEHgmin,dEemin

      RETURN
      END

*********************************************************************

      DOUBLE PRECISION FUNCTION AEDM(x)

      IMPLICIT NONE
      DOUBLE PRECISION x
      IF(dabs(x-1d0).gt.1.d-8)THEN
      AEDM=1d0/(1d0-x)**3*(3d0-4d0*x+x**2+2d0*dlog(x))/2d0
      ELSE
      AEDM=-1d0/3.d0
      ENDIF
      RETURN
      END

*********************************************************************

      DOUBLE PRECISION FUNCTION BEDM(x)

      IMPLICIT NONE
      DOUBLE PRECISION x
      IF(dabs(x).gt.1.d-8)THEN
      IF(dabs(x-1d0).gt.1.d-3)THEN
      BEDM=1d0/(1d0-x)**3*(1d0-x**2+2d0*x*dlog(x))/2d0
      ELSE
      BEDM=1d0/6d0
      ENDIF
      ELSE
      BEDM=1d0/2d0      
      ENDIF
      RETURN
      END

*********************************************************************

      DOUBLE PRECISION FUNCTION CEDM(x)

      IMPLICIT NONE
      DOUBLE PRECISION x
      IF(dabs(x-1d0).gt.1.d-8)THEN
      CEDM=1d0/(1d0-x)**3*(-26d0+36*x-10*x**2+(2d0*x-18d0)*dlog(x))
     .                                                /6d0
      ELSE
      CEDM=19d0/18d0
      ENDIF
      RETURN
      END


*********************************************************************

      DOUBLE COMPLEX FUNCTION Cdilog(z)

      IMPLICIT NONE
      INTEGER I
      DOUBLE PRECISION Pi,A(8),W(8)
      DOUBLE COMPLEX z,x,aux,auxi

      Pi=4.d0*datan(1.d0)
      A(1)=0.095012509837637440185d0
      A(2)=0.281603550779258913230d0
      A(3)=0.458016777657227386342d0
      A(4)=0.617876244402643748447d0
      A(5)=0.755404408355003033895d0
      A(6)=0.865631202387831743880d0
      A(7)=0.944575023073232576078d0
      A(8)=0.989400934991649932595d0
      W(1)=0.189450610455068496285d0
      W(2)=0.182603415044923588867d0
      W(3)=0.169156519395002538189d0
      W(4)=0.149595988816576732081d0
      W(5)=0.124628971255533872052d0
      W(6)=0.095158511682492784810d0
      W(7)=0.062253523938647892863d0
      W(8)=0.027152459411754094852d0

      x=z
      IF(CDABS(x).ge.1.d0.and.CDABS(x).ne.1.d0)x=1.d0/x
      IF(CDABS(x).le.1.d0.and.CDABS(1.d0-x).le.0.5d0)x=1.d0-x

      aux=DCMPLX(0.d0,0.d0)
      IF(CDABS(x).le.1.d0)THEN
      IF(CDABS(x).le.0.5d0)THEN
       IF(CDabs(x-DCMPLX(-2.d0+dsqrt(3.d0),0.d0)).le.1d-6)THEN
        aux=DCMPLX(-.2518620182d0,0.d0)
       ELSE
        auxi=DCMPLX(0.d0,0.d0)
        DO I=1,15
        auxi=auxi+x**I/(I*(I+1.d0)*(I+2.d0))**2
        ENDDO
        aux=(4.d0*x+23.d0/4.d0*x**2+3.d0*(1-x**2)*CDlog(1.d0-x)
     .    +4.d0*x**2*auxi)/(1+4.d0*x+x**2)
       ENDIF
      ENDIF
      IF(CDABS(x).ge.0.5d0.and.CDABS(x).ne.0.5d0)THEN
       auxi=DCMPLX(0.d0,0.d0)
       DO I=1,8
        auxi=auxi-W(I)*CDlog(1.d0-x*(A(I)+1.d0)/2.d0)/(A(I)+1.d0)
     .           -W(I)*CDlog(1.d0-x*(-A(I)+1.d0)/2.d0)/(-A(I)+1.d0)
       ENDDO
       aux=auxi
      ENDIF
      ENDIF

      IF(CDABS(z).ge.1.d0.and.CDABS(z).ne.1.d0)THEN
       IF(CDABS(1.d0-1.d0/z).le.0.5d0)THEN
        aux=-aux-CDlog(1.d0/z)*CDlog(1.d0-1.d0/z)+Pi**2/6.d0
       ENDIF
       Cdilog=-aux-CDlog(-z)**2/2.d0-Pi**2/6.d0
      ENDIF
      IF(CDABS(z).le.1.d0)THEN
       IF(CDABS(1.d0-z).le.0.5d0)THEN
        Cdilog=-aux-CDlog(z)*CDlog(1.d0-z)+Pi**2/6.d0
       ELSE
        Cdilog=aux
       ENDIF
      ENDIF
      IF(CDABS(z-1.d0).le.1d-6)Cdilog=DCMPLX(Pi**2/6.d0,0.d0)

      RETURN
      END

*********************************************************************

      DOUBLE PRECISION FUNCTION fWW(a,b)

      IMPLICIT NONE
      DOUBLE PRECISION a,b,Delt
      DOUBLE COMPLEX xp,xm,aux,Cdilog,eps

      Delt=(1.d0+a-b)**2-4.d0*a
      IF(Delt.ge.0.d0)THEN
       xp=DCMPLX((1.d0+a-b+dsqrt(Delt))/2.d0,-1d-10)
       xm=DCMPLX((1.d0+a-b-dsqrt(Delt))/2.d0,-1d-10)
       eps=DCMPLX(1.d0,0.d0)
      ELSE
       xp=DCMPLX((1.d0+a-b)/2.d0,dsqrt(-Delt)/2.d0)
       xm=DCMPLX((1.d0+a-b)/2.d0,-dsqrt(-Delt)/2.d0)
       eps=DCMPLX(0.d0,-1.d0)
      ENDIF
      aux=xp*(CDlog((a*(1.d0-xp)+b*xp)/(xp*(1.d0-xp)))
     .                                        *CDlog(-(1-xp)/xp)
     .         -Cdilog((b-a)*(xp-1.d0)/(a*(1.d0-xp)+b*xp))
     .         +Cdilog((b-a)*xp/(a*(1.d0-xp)+b*xp))
     .         -Cdilog(-xp/(1.d0-xp))+Cdilog(-(1.d0-xp)/xp))
     .     -xm*(CDlog((a*(1.d0-xm)+b*xm)/(xm*(1.d0-xm)))
     .                                        *CDlog(-(1-xm)/xm)
     .         -Cdilog((b-a)*(xm-1.d0)/(a*(1.d0-xm)+b*xm))
     .         +Cdilog((b-a)*xm/(a*(1.d0-xm)+b*xm))
     .         -Cdilog(-xm/(1.d0-xm))+Cdilog(-(1.d0-xm)/xm))
      aux=aux*eps/dsqrt(dabs(Delt))
      fWW=DREAL(aux)

      RETURN
      END

*********************************************************************

      DOUBLE PRECISION FUNCTION HEDM(x)

      IMPLICIT NONE
      DOUBLE PRECISION x,aux
      IF(dabs(x-1d0).gt.1.d-4)THEN
      aux=x**2/(x-1)**3*(dlog(x)-3d0/2d0+2/x-1/2d0/x**2)
      ELSE
      aux=1d0/3d0
      ENDIF

      HEDM=aux*(1d0-0.813301412284338499d0+5.46546428719477811d-2*
     . datan(0.405192104968222777d0*dlog(x)-0.696827935267380294d0))
      RETURN
      END

*********************************************************************

      DOUBLE PRECISION FUNCTION HQEDM(zg,zs,delt)

      IMPLICIT NONE
      DOUBLE PRECISION zg,zs,delt,aux1,aux2

      IF(dabs(zs-1d0).ge.1d-4)THEN
      aux1=(2d0*(1d0-zs)*(1d0+11d0*zs)-(1d0-16d0*zs-9d0*zs**2)*dlog(zs))
     .          /18d0/(1d0-zs)**4
      aux2=((1d0-zs)*(1d0+7d0*zs+295d0*zs**2+177d0*zs**3)+6d0*zs**2*
     .    (21d0+50d0*zs+9d0*zs**2)*dlog(zs))/108d0/(1d0-zs)**6
      ELSE
      aux1=5d0/108d0
      aux2=11d0/1080d0
      ENDIF
      HQEDM=(aux1+delt**2*aux2)/zg

      RETURN
      END
