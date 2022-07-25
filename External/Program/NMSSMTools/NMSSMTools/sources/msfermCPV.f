      SUBROUTINE MSFERM_CPV(PAR,IFAIL)

c       Sfermion masses: diag(MSF2(I)) = UF.MSF2.UF+ (L,R) base
c      - The squared masses of the sfermions of 3rd generation MST2(i) (stops),
c        MSB2(i) (sbottoms), MSL2(i) (staus), (i=1,2), MSNT2 (sneutrino) as 
c        well as the rotation matrices UT(i,j,k) (stops), UB(i,j,k) (sbottoms),
c        UTAU(i,j,k) (staus), (i,j=1..2; k=1..2 for real/imaginary part),
c        are stored in the common SFERM3SPEC.
c      - The squared masses of the sfermion of 1st/2nd generation, assumed 
c        unmixed, are stored as MSU2(i) (sups), MSD2(i) (sdowns), MSE2(i)
c        (selectrons), MSNE2 (sneutrino), (i=1..2), in the common SFERM1SPEC.
c      - The squared squark masses including alpha_S corrections, MST2P(i)
c        (stops), MSB2P(i) (sbottoms), MSU2P(i) (sups), MSD2P(i) (sdowns),
c        (i=1..2) are stored in SFERMPSPEC.
c      IF a negative squared mass appears, IFAIL=8.
c
c      The corrections to the bottom Yukawa DELMB: Yb = mb/[vd*(1+DELMB)] are
c      also computed (generalizing the results of hep-ph/9912516) and stored
c      in the common DELMB.

      IMPLICIT NONE

      INTEGER I,IFAIL

      DOUBLE PRECISION PAR(*)
      DOUBLE PRECISION Pi,M3,DMSQUARK_CPV,PhiST,thetaST,PhiSB,thetaSB,
     . phiSL,thetaSL,Phismu,thetasmu,integ
      DOUBLE PRECISION MSF2_11,MSF2_22,MSF4_12,Tr,det,aux1,aux2
      DOUBLE PRECISION DPIST(2),DPISB(2)

      DOUBLE PRECISION Q2
      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      DOUBLE PRECISION Ytq,Ybq,MTOPQ,MBOTQ
      DOUBLE PRECISION MSQ3,MSU3,MSD3,AT,AB
      DOUBLE PRECISION MSQ1,MSU1,MSD1
      DOUBLE PRECISION l,k,Alcos1,Akcos2,muq,nuq
      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION DELMB,DELML,DEL1
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION MST2(2),UT(2,2,2),MSB2(2),UB(2,2,2),MSL2(2),
     . UTAU(2,2,2),MSNT2
      DOUBLE PRECISION MSU2(2),MSD2(2),MSE2(2),MSNE2,MSMU2(2),
     . UMU(2,2,2)
      DOUBLE PRECISION MST2P(2),MSB2P(2),MSU2P(2),MSD2P(2)
      DOUBLE PRECISION MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/SMFERM/mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      COMMON/QQUARK/Ytq,Ybq,MTOPQ,MBOTQ
      COMMON/RADCOR2/MSQ3,MSU3,MSD3,AT,AB
      COMMON/SQUPAR/MSQ1,MSU1,MSD1
      COMMON/QPAR/l,k,Alcos1,Akcos2,muq,NUQ
      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/SFERM3SPEC/MST2,UT,MSB2,UB,MSL2,UTAU,MSNT2
      COMMON/SFERM1SPEC/MSU2,MSD2,MSE2,MSNE2,MSMU2,UMU
      COMMON/SFERMPSPEC/MST2P,MSB2P,MSU2P,MSD2P
      COMMON/SLEPPAR/MSL3,MSE3,MSL1,MSE1,ATAU,AMU

************************************************************************************************

      PI=4d0*DATAN(1d0)
      M3=PAR(22)

c                 A: 3rd generation
c        I- STOPS

c      1) Squared-mass matrix
      MSF2_11=MSQ3+(Ytq*vuq)**2+(g1q/3.d0-g2q)/4.d0
     .                                 *(vuq**2-vdq**2)
      MSF2_22=MSU3+(Ytq*vuq)**2+(-4.d0/3.d0*g1q)/4.d0
     .                                 *(vuq**2-vdq**2)
      MSF4_12=(Ytq*vuq)**2*(At**2+(muq*vdq/vuq)**2
     .             -2d0*At*muq*vdq/vuq*DDCOS(Phi01+PhiAT))

c      2) Eigenvalues
      Tr=MSF2_11+MSF2_22
      det=(MSF2_11-MSF2_22)**2+4.d0*MSF4_12
      MST2(1)=(Tr-dsqrt(det))/2d0
      MST2(2)=(Tr+dsqrt(det))/2d0
      IF(MST2(1).le.0d0)THEN
       IFAIL=8
c       print*,'stoppro'
       GOTO 618
      ENDIF

c      3) Mixing angles
      aux1=Ytq*vuq*(-At*DDSIN(PhiAT)-muq*vdq/vuq*DDSIN(Phi01))
      aux2=Ytq*vuq*(At*DDCOS(PhiAT)-muq*vdq/vuq*DDCOS(Phi01))
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       phiST=datan(aux1/aux2)
       IF(aux2.le.0d0)PhiST=PhiST+Pi
      ELSEIF(aux1.ge.0d0)THEN
       phiST=Pi/2d0
      ELSE
       phiST=-Pi/2d0
      ENDIF

      aux1=MSQ3-MSU3+(g1q*(1d0/3.d0-(-4.d0/3.d0))-g2q)
     .           *(vuq**2-vdq**2)/4.d0+dsqrt(det)
      aux2=dsqrt(4.d0*MSF4_12)
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       thetaST=datan(aux1/aux2)
       IF(aux2.le.0d0)thetaST=thetaST+Pi
      ELSEIF(aux1.ge.0d0)THEN
       thetaST=Pi/2d0
      ELSE
       thetaST=-Pi/2d0
      ENDIF

c      4) Rotation matrix
      UT(1,1,1)=DDCOS(thetaST)
      UT(1,1,2)=0d0
      UT(1,2,1)=-DDSIN(thetaST)*DDCOS(PhiST)
      UT(1,2,2)=-DDSIN(thetaST)*DDSIN(PhiST)
      UT(2,1,1)=DDSIN(thetaST)*DDCOS(PhiST)
      UT(2,1,2)=-DDSIN(thetaST)*DDSIN(PhiST)
      UT(2,2,1)=DDCOS(thetaST)
      UT(2,2,2)=0d0
c      UTCOMP(1,1)=DCMPLX(UT(1,1,1),UT(1,1,2))
c      UTCOMP(1,2)=DCMPLX(UT(1,2,1),UT(1,2,2))
c      UTCOMP(2,1)=DCMPLX(UT(2,1,1),UT(2,1,2))
c      UTCOMP(2,2)=DCMPLX(UT(2,2,1),UT(2,2,2))

c      5) Radiative corrections
      aux1=MST2(1)
      aux2=MST2(2)
      Tr=2d0*(UT(1,1,1)*(UT(1,2,1)*DDCOS(PhiM3)
     .                                -UT(1,2,2)*DDSIN(PhiM3)))
      MST2P(1)=MST2(1)-alsq/3.d0/Pi
     . *DMSQUARK_CPV(1,aux1,aux2,Tr,DDCOS(2d0*thetaST)**2,
     .  DDSIN(2d0*thetaST)**2*(1d0+DDCOS(2d0*PhiST))/2d0,
     .  mtopq,M3,QSTSB)
      MST2P(2)=MST2(2)-alsq/3.d0/Pi
     . *DMSQUARK_CPV(2,aux1,aux2,Tr,DDCOS(2d0*thetaST)**2,
     .  DDSIN(2d0*thetaST)**2*(1d0+DDCOS(2d0*PhiST))/2d0,
     .  mtopq,M3,QSTSB)

c      print*,'MST(1)',dsqrt(MST2P(1))
c      print*,'MST(2)',dsqrt(MST2P(2))
c      print*,'UST11',UT(1,1,1),UT(1,1,2)
c      print*,'UST12',UT(1,2,1),UT(1,2,2)
c      print*,'UST21',UT(2,1,1),UT(2,1,2)
c      print*,'UST22',UT(2,2,1),UT(2,2,2)


c        II- SBOTTOMS

c      1) Squared-mass matrix
      MSF2_11=MSQ3+(Ybq*vdq)**2+(g1q/3.d0+g2q)/4.d0
     .                                 *(vuq**2-vdq**2)
      MSF2_22=MSD3+(Ybq*vdq)**2+(2d0/3.d0*g1q)/4.d0
     .                                 *(vuq**2-vdq**2)
      MSF4_12=(Ybq*vdq)**2*(Ab**2+(muq*vuq/vdq)**2
     .               -2d0*Ab*muq*vuq/vdq*DDCOS(Phi01+PhiAb))

c      2) Eigenvalues
      Tr=MSF2_11+MSF2_22
      det=(MSF2_11-MSF2_22)**2+4.d0*MSF4_12
      MSB2(1)=(Tr-dsqrt(det))/2d0
      MSB2(2)=(Tr+dsqrt(det))/2d0
      IF(MSB2(1).le.0d0)THEN
c       print*,'sbotpro'
       IFAIL=8
       GOTO 618
      ENDIF

c      3) Mixing angles
      aux1=Ybq*vdq*(-Ab*DDSIN(PhiAb)-muq*vuq/vdq*DDSIN(Phi01))
      aux2=Ybq*vdq*(Ab*DDCOS(PhiAb)-muq*vuq/vdq*DDCOS(Phi01))
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       phiSB=datan(aux1/aux2)
       IF(aux2.le.0d0)PhiSB=PhiSB+Pi
      ELSEIF(aux1.ge.0d0)THEN
       phiSB=Pi/2d0
      ELSE
       phiSB=-Pi/2d0
      ENDIF

      aux1=MSQ3-MSD3+(g1q*(1d0/3.d0-(2d0/3.d0))+g2q)
     .           *(vuq**2-vdq**2)/4.d0+dsqrt(det)
      aux2=dsqrt(4.d0*MSF4_12)
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       thetaSB=datan(aux1/aux2)
       IF(aux2.le.0d0)thetaSB=thetaSB+Pi
      ELSEIF(aux1.ge.0d0)THEN
       thetaSB=Pi/2d0
      ELSE
       thetaSB=-Pi/2d0
      ENDIF

c      4) Rotation matrix
      UB(1,1,1)=DDCOS(thetaSB)
      UB(1,1,2)=0d0
      UB(1,2,1)=-DDSIN(thetaSB)*DDCOS(PhiSB)
      UB(1,2,2)=-DDSIN(thetaSB)*DDSIN(PhiSB)
      UB(2,1,1)=DDSIN(thetaSB)*DDCOS(PhiSB)
      UB(2,1,2)=-DDSIN(thetaSB)*DDSIN(PhiSB)
      UB(2,2,1)=DDCOS(thetaSB)
      UB(2,2,2)=0d0
c      UBCOMP(1,1)=DCMPLX(UB(1,1,1),UB(1,1,2))
c      UBCOMP(1,2)=DCMPLX(UB(1,2,1),UB(1,2,2))
c      UBCOMP(2,1)=DCMPLX(UB(2,1,1),UB(2,1,2))
c      UBCOMP(2,2)=DCMPLX(UB(2,2,1),UB(2,2,2))

c      5) Radiative corrections
      aux1=MSB2(1)
      aux2=MSB2(2)
      Tr=2d0*(UB(1,1,1)*(UB(1,2,1)*DDCOS(PhiM3)
     .                                -UB(1,2,2)*DDSIN(PhiM3)))
      MSB2P(1)=MSB2(1)-alsq/3.d0/Pi
     . *DMSQUARK_CPV(1,aux1,aux2,Tr,DDCOS(2d0*thetaSB)**2,
     .  DDSIN(2d0*thetaSB)**2*(1d0+DDCOS(2d0*PhiSB))/2d0,
     .  mbotq,M3,QSTSB)
      MSB2P(2)=MSB2(2)-alsq/3.d0/Pi
     . *DMSQUARK_CPV(2,aux1,aux2,Tr,DDCOS(2d0*thetaSB)**2,
     .  DDSIN(2d0*thetaSB)**2*(1d0+DDCOS(2d0*PhiSB))/2d0,
     .  mbotq,M3,QSTSB)

c      print*,'MSB(1)',dsqrt(MSB2P(1))
c      print*,'MSB(2)',dsqrt(MSB2P(2))
c      print*,'USB11',UB(1,1,1),UB(1,1,2)
c      print*,'USB12',UB(1,2,1),UB(1,2,2)
c      print*,'USB21',UB(2,1,1),UB(2,1,2)
c      print*,'USB22',UB(2,2,1),UB(2,2,2)

c      6) Yukawa corrections to the stop and sbottom masses
      CALL DMSQUARK_YUK_CPV(DPIST,DPISB)
      DO I=1,2
       MST2P(I)=MST2P(I)-DPIST(I)
       MSB2P(I)=MSB2P(I)-DPISB(I)
      ENDDO

c        III- STAUS

      MSF2_11=MSL3+mtau**2+(-g1q+g2q)/4.d0*(vuq**2-vdq**2)
      MSF2_22=MSE3+mtau**2+(2d0*g1q)/4.d0*(vuq**2-vdq**2)
      MSF4_12=mtau**2*(Atau**2+muq**2*(vuq/vdq)**2
     .                -2d0*Atau*muq*vuq/vdq*DDCOS(Phi01+PhiAtau))

      Tr=MSF2_11+MSF2_22
      det=(MSF2_11-MSF2_22)**2+4.d0*MSF4_12
      MSL2(1)=(Tr-dsqrt(det))/2d0
      MSL2(2)=(Tr+dsqrt(det))/2d0
      IF(MSL2(1).le.0d0)THEN
c       print*,'staupro'
       IFAIL=8
       GOTO 618
      ENDIF

      aux1=mtau*(-Atau*DDSIN(PhiAtau)-muq*vuq/vdq*DDSIN(Phi01))
      aux2=mtau*(Atau*DDCOS(PhiAtau)-muq*vuq/vdq*DDCOS(Phi01))
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       phiSL=datan(aux1/aux2)
       IF(aux2.le.0d0)PhiSL=PhiSL+Pi
      ELSEIF(aux1.ge.0d0)THEN
       phiSL=Pi/2d0
      ELSE
       phiSL=-Pi/2d0
      ENDIF

      aux1=MSL3-MSE3+(g1q*(-1d0-2d0)+g2q)
     .           *(vuq**2-vdq**2)/4.d0+dsqrt(det)
      aux2=dsqrt(4.d0*MSF4_12)
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       thetaSL=datan(aux1/aux2)
       IF(aux2.le.0d0)thetaSL=thetaSL+Pi
      ELSEIF(aux1.ge.0d0)THEN
       thetaSL=Pi/2d0
      ELSE
       thetaSL=-Pi/2d0
      ENDIF

      UTAU(1,1,1)=DDCOS(thetaSL)
      UTAU(1,1,2)=0d0
      UTAU(1,2,1)=-DDSIN(thetaSL)*DDCOS(PhiSL)
      UTAU(1,2,2)=-DDSIN(thetaSL)*DDSIN(PhiSL)
      UTAU(2,1,1)=DDSIN(thetaSL)*DDCOS(PhiSL)
      UTAU(2,1,2)=-DDSIN(thetaSL)*DDSIN(PhiSL)
      UTAU(2,2,1)=DDCOS(thetaSL)
      UTAU(2,2,2)=0d0
c      UTAUCOMP(1,1)=DCMPLX(UTAU(1,1,1),UTAU(1,1,2))
c      UTAUCOMP(1,2)=DCMPLX(UTAU(1,2,1),UTAU(1,2,2))
c      UTAUCOMP(2,1)=DCMPLX(UTAU(2,1,1),UTAU(2,1,2))
c      UTAUCOMP(2,2)=DCMPLX(UTAU(2,2,1),UTAU(2,2,2))

c      print*,'MSL(1)',dsqrt(MSL2(1))
c      print*,'MSL(2)',dsqrt(MSL2(2))
c      print*,'USTAU11',UTAU(1,1,1),UTAU(1,1,2)
c      print*,'USTAU12',UTAU(1,2,1),UTAU(1,2,2)
c      print*,'USTAU21',UTAU(2,1,1),UTAU(2,1,2)
c      print*,'USTAU22',UTAU(2,2,1),UTAU(2,2,2)


c        IV- SNEUTRINO

      MSNT2=MSL3+(-g1q-g2q)/4.d0*(vuq**2-vdq**2)
      IF(MSNT2.le.0d0)THEN
c       print*,'sneupro'
       IFAIL=8
       GOTO 618
      ENDIF

c      print*,'MSNT',dsqrt(MSNT2)

c                 B: Light generations (assumed degenerate)

c        I- SUP/SCHARM

      MSU2(1)=MSQ1+(g1q/3.d0-g2q)/4.d0*(vuq**2-vdq**2) ! L
      MSU2(2)=MSU1-4.d0*g1q/3.d0/4.d0*(vuq**2-vdq**2)  ! R
      IF(min(MSU2(1),MSU2(2)).le.0d0)THEN
c       print*,'suppro'
       IFAIL=8
       GOTO 618
      ENDIF

c      Radiative corrections
      DO I=1,2
      aux2=M3**2/MSU2(I)
      IF(aux2.NE.1d0)THEN
       aux1=1d0+3d0*aux2+1d0/2d0*(1d0-aux2)**2*DLOG((1d0-aux2)**2)
     .     -aux2**2*DLOG(aux2)+2d0*aux2*DLOG(Q2/MSU2(I))
      ELSE
       aux1=4d0+2d0*DLOG(Q2/MSU2(I))
      ENDIF
      MSU2P(I)=MSU2(I)*(1d0+2d0*alsq*aux1/3.d0/Pi)
      ENDDO

c      print*,'MSU',dsqrt(MSU2(1)),dsqrt(MSU2(2))

c        II- SDOWN/SSTRANGE

      MSD2(1)=MSQ1+(g1q/3.d0+g2q)/4.d0*(vuq**2-vdq**2) ! L
      MSD2(2)=MSD1+2d0*g1q/3.d0/4.d0*(vuq**2-vdq**2)  ! R
      IF(min(MSD2(1),MSD2(2)).le.0d0)THEN
c       print*,'sdownpro'
       IFAIL=8
       GOTO 618
      ENDIF

c      Radiative corrections
      DO I=1,2
      aux2=M3**2/MSD2(I)
      IF(aux2.NE.1d0)THEN
       aux1=1d0+3d0*aux2+1d0/2d0*(1d0-aux2)**2*DLOG((1d0-aux2)**2)
     .     -aux2**2*DLOG(aux2)+2d0*aux2*DLOG(Q2/MSD2(I))
      ELSE
       aux1=4d0+2d0*DLOG(Q2/MSD2(I))
      ENDIF
      MSD2P(I)=MSD2(I)*(1d0+2d0*alsq*aux1/3.d0/Pi)
      ENDDO

c      print*,'MSD',dsqrt(MSD2(1)),dsqrt(MSD2(2))

c        III- SELECTRON/SMUON

      MSE2(1)=MSL1+(-g1q+g2q)/4.d0*(vuq**2-vdq**2)
      MSE2(2)=MSE1+2d0*g1q/4.d0*(vuq**2-vdq**2)
      IF(min(MSE2(1),MSE2(2)).le.0d0)THEN
c       print*,'selpro'
       IFAIL=8
       GOTO 618
      ENDIF

c      print*,'MSE',dsqrt(MSE2(1)),dsqrt(MSE2(2))

c      smuons: for g-2

      MSF2_11=MSL1+mmu**2+(-g1q+g2q)/4.d0*(vuq**2-vdq**2)
      MSF2_22=MSE1+mmu**2+(2d0*g1q)/4.d0*(vuq**2-vdq**2)
      MSF4_12=mmu**2*(Amu**2+muq**2*(vuq/vdq)**2
     .                -2d0*Amu*muq*vuq/vdq*DDCOS(Phi01+PhiAmu))

      Tr=MSF2_11+MSF2_22
      det=(MSF2_11-MSF2_22)**2+4.d0*MSF4_12
      MSMU2(1)=(Tr-dsqrt(det))/2d0
      MSMU2(2)=(Tr+dsqrt(det))/2d0
      IF(MSMU2(1).le.0d0)THEN
c       print*,'smupro'
       IFAIL=8
       GOTO 618
      ENDIF

      aux1=mmu*(-Amu*DDSIN(PhiAmu)-muq*vuq/vdq*DDSIN(Phi01))
      aux2=mmu*(Amu*DDCOS(PhiAmu)-muq*vuq/vdq*DDCOS(Phi01))
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       phiSMU=datan(aux1/aux2)
       IF(aux2.le.0d0)PhiSMU=PhiSMU+Pi
      ELSEIF(aux1.ge.0d0)THEN
       phiSMU=Pi/2d0
      ELSE
       phiSMU=-Pi/2d0
      ENDIF

      aux1=MSL1-MSE1+(g1q*(-1d0-2d0)+g2q)
     .           *(vuq**2-vdq**2)/4.d0+dsqrt(det)
      aux2=dsqrt(4.d0*MSF4_12)
      IF(dabs(aux2).ge.dabs(aux1)*1.d-10)THEN
       thetaSMU=datan(aux1/aux2)
       IF(aux2.le.0d0)thetaSmu=thetaSmu+Pi
      ELSEIF(aux1.ge.0d0)THEN
       thetaSmu=Pi/2d0
      ELSE
       thetaSmu=-Pi/2d0
      ENDIF

      UMU(1,1,1)=DDCOS(thetaSmu)
      UMU(1,1,2)=0d0
      UMU(1,2,1)=-DDSIN(thetaSmu)*DDCOS(PhiSmu)
      UMU(1,2,2)=-DDSIN(thetaSmu)*DDSIN(PhiSmu)
      UMU(2,1,1)=DDSIN(thetaSmu)*DDCOS(PhiSmu)
      UMU(2,1,2)=-DDSIN(thetaSmu)*DDSIN(PhiSmu)
      UMU(2,2,1)=DDCOS(thetaSmu)
      UMU(2,2,2)=0d0

c      print*,'MSL(1)',dsqrt(MSmu2(1))
c      print*,'MSL(2)',dsqrt(MSmu2(2))
c      print*,'USMU11',UmU(1,1,1),UmU(1,1,2)
c      print*,'USMU12',UmU(1,2,1),UmU(1,2,2)
c      print*,'USMU21',UmU(2,1,1),UmU(2,1,2)
c      print*,'USMU22',UmU(2,2,1),UmU(2,2,2)

c        IV- SNEUTRINOS

      MSNE2=MSL1+(-g1q-g2q)/4.d0*(vuq**2-vdq**2)
      IF(MSNE2.le.0d0)THEN
c       print*,'sneupro'
       IFAIL=8
       GOTO 618
      ENDIF

c      print*,'MSNE',dsqrt(MSNE2)

c      print*,'MSL11',MSL3+mtau**2+(-g1q+g2q)/4.d0*(vuq**2-vdq**2)
c      print*,'MSL12',mtau*(Atau*DDCOS(PhiATAU)-muq*vuq/vdq*DDCOS(Phi01)),
c     c -mtau*(Atau*DDSIN(PhiATAU)+muq*vuq/vdq*DDSIN(Phi01))
c      print*,'MSL21',mtau*(Atau*DDCOS(PhiATAU)-muq*vuq/vdq*DDCOS(Phi01)),
c     c mtau*(Atau*DDSIN(PhiATAU)+muq*vuq/vdq*DDSIN(Phi01))
c      print*,'MSL22',MSE3+mtau**2+(2d0*g1q)/4.d0*(vuq**2-vdq**2)


c                 C: Corrections to the b and tau Yukawas

      DEL1=-2d0/(3d0*PI)*ALSQ*M3*AB*DDCOS(PhiM3+PhiAB)
     .                     *INTEG(dsqrt(MSB2P(1)),dsqrt(MSB2P(2)),M3)

      DELMB=MUQ*TANBQ*(2d0/(3d0*PI)*ALSMZ*DDCOS(PhiM3+Phi01)*M3
     .                     *INTEG(dsqrt(MSB2P(1)),dsqrt(MSB2P(2)),M3)
     . +Ytq**2/16.d0/Pi**2*AT*DDCOS(PhiAT+Phi01)
     .                     *INTEG(dsqrt(MST2P(1)),dsqrt(MST2P(2)),MUQ)
     . -G2Q/16.d0/Pi**2*PAR(21)*DDCOS(PhiM2+Phi01)*(
     .   (UT(1,1,1)**2+UT(1,1,2)**2)*INTEG(dsqrt(MST2P(1)),PAR(21),MUQ)
     .  +(UT(2,1,1)**2+UT(2,1,2)**2)*INTEG(dsqrt(MST2P(2)),PAR(21),MUQ)
     .  +((UB(1,1,1)**2+UB(1,1,2)**2)*INTEG(dsqrt(MSB2P(1)),PAR(21),MUQ)
     .  +(UB(1,2,1)**2+UB(1,2,2)**2)*INTEG(dsqrt(MSB2P(2)),PAR(21),MUQ))
     .                                                           /2d0)
     . -G1Q/16.d0/Pi**2/3d0*PAR(20)*DDCOS(PhiM1+Phi01)
     .   *(1d0/3d0*INTEG(dsqrt(MSB2P(1)),dsqrt(MSB2P(2)),PAR(20))
     .   +1d0/2d0*((UB(1,1,1)**2+UB(1,1,2)**2)
     .    *INTEG(dsqrt(MSB2P(1)),PAR(20),MUQ)
     .   +(UB(1,2,1)**2+UB(1,2,2)**2)
     .    *INTEG(dsqrt(MSB2P(2)),PAR(20),MUQ))
     .   +(UB(1,2,1)**2+UB(1,2,2)**2)
     .    *INTEG(dsqrt(MSB2P(1)),PAR(20),MUQ)
     .   +(UB(1,1,1)**2+UB(1,1,2)**2)
     .    *INTEG(dsqrt(MSB2P(2)),PAR(20),MUQ)))
     .   /(1d0+DEL1)

      DELML=MUQ*TANBQ*(
     . -G2Q/16.d0/Pi**2*PAR(21)*DDCOS(PhiM2+Phi01)
     .   *(INTEG(dsqrt(MSNT2),PAR(21),MUQ)
     .   +1d0/2d0*((UTAU(1,1,1)**2+UTAU(1,1,2)**2)
     .    *INTEG(dsqrt(MSL2(1)),PAR(21),MUQ)
     .   +(UTAU(1,2,1)**2+UTAU(1,2,2)**2)
     .    *INTEG(dsqrt(MSL2(2)),PAR(21),MUQ)))
     . +G1Q/16.d0/Pi**2*PAR(20)*DDCOS(PhiM1+Phi01)
     .   *(INTEG(dsqrt(MSL2(1)),dsqrt(MSL2(2)),PAR(20))
     .   +1d0/2d0*((UTAU(1,1,1)**2+UTAU(1,1,2)**2)
     .    *INTEG(dsqrt(MSL2(1)),PAR(20),MUQ)
     .   +(UTAU(1,2,1)**2+UTAU(1,2,2)**2)
     .    *INTEG(dsqrt(MSL2(2)),PAR(20),MUQ))
     .   -(UTAU(1,2,1)**2+UTAU(1,2,2)**2)
     .    *INTEG(dsqrt(MSL2(1)),PAR(20),MUQ)
     .   -(UTAU(1,1,1)**2+UTAU(1,1,2)**2)
     .    *INTEG(dsqrt(MSL2(2)),PAR(20),MUQ)))

c      print*,'DEL1',DEL1
c      print*,'DELMB',DELMB
c      print*,'DELML',DELML

 618  RETURN
      END


      DOUBLE PRECISION FUNCTION DMSQUARK_CPV(i,msq1,msq2,s2t,c2t2,s2t2,
     . mq,mglu,q2)

*    with thanks to S. Kraml
*    msq1, msq2 are the squark masses squared

      IMPLICIT NONE
      DOUBLE PRECISION msq1,msq2,s2t,c2t2,s2t2,mq,mglu,msq,msqp
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
      dmsq = c2t2*NMA0(msq,q2) + s2t2*NMA0(msqp,q2)
      DMSQUARK_CPV = dmg + dmsg + dmsq

      END


      SUBROUTINE DMSQUARK_YUK_CPV(DPIST,DPISB)
      
*******************************************************************
* Subroutine to compute the Yukawa contributions to the squark pole masses
* based on Pierce/Bagger et al., hep-ph/9606211 eqs. (D.47-49)
*
* NOTE: mu(Pierce/Bagger) = -MUQ(this code)!!!

*******************************************************************

      IMPLICIT NONE
      
      INTEGER I,J,K
      
      DOUBLE PRECISION DPIST(2),DPISB(2)
      DOUBLE PRECISION Pi,MW2,MZ2,sinbq,cosbq,NMA0,NMB0
      DOUBLE PRECISION MHTC(2),MHT2(6),XHT(6,6)
      DOUBLE PRECISION COH0SFL(6,2,2),COH0SFR(6,2,2),COHCSFL(6,2,2),
     .      COHCSFR(6,2,2)

      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION Ytq,Ybq,MTOPQ,MBOTQ
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION l,ka,Alcos1,Akcos2,muq,NUQ
      DOUBLE PRECISION MSQ3,MSU3,MSD3,AT,AB
      DOUBLE PRECISION mur,M1r,M2r,msi
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION MST2(2),UT(2,2,2),MSB2(2),UB(2,2,2),MSL2(2),
     . UTAU(2,2,2),MSNT2
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/STSBSCALE/QSTSB
      COMMON/QQUARK/Ytq,Ybq,MTOPQ,MBOTQ
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/QPAR/l,ka,Alcos1,Akcos2,muq,NUQ
      COMMON/RADCOR2/MSQ3,MSU3,MSD3,AT,AB
      COMMON/GAUGINOPAR/mur,M1r,M2r,msi
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/SFERM3SPEC/MST2,UT,MSB2,UB,MSL2,UTAU,MSNT2

      PI=4d0*DATAN(1d0)

*      I- Adding the Goldstone to the Higgs matrices

      MW2=g2q/2d0*(vuq**2+vdq**2)
      MZ2=(g1q+g2q)/2d0*(vuq**2+vdq**2)
      sinbq=vuq/dsqrt(vuq**2+vdq**2)
      cosbq=vdq/dsqrt(vuq**2+vdq**2)

      MHTC(1)=MW2
      MHTC(2)=MHC

      DO I=1,5
       MHT2(I)=MH0(I)
       XHT(I,1)=XH(I,1)
       XHT(I,2)=XH(I,2)
       XHT(I,3)=XH(I,3)
       XHT(I,4)=XH(I,4)*cosbq
       XHT(I,5)=XH(I,4)*sinbq
       XHT(I,6)=XH(I,5)
      ENDDO
      MHT2(6)=MZ2
      DO I=1,6
       XHT(6,I)=0d0
      ENDDO
      XHT(6,4)=-sinbq
      XHT(6,5)=cosbq


*      II- Corrections to the stop masses

      DO I=1,2
       DPIST(I)=0d0

*       1) Sfermion Quartic couplings
      DO k=1,2
       DPIST(I)=DPIST(I)+Ytq**2*(
     .  (UT(I,1,1)**2+UT(I,1,2)**2)*(UT(k,2,1)**2+UT(k,2,2)**2)
     . +(UT(I,2,1)**2+UT(I,2,2)**2)*(UT(k,1,1)**2+UT(k,1,2)**2)
     . +6d0*((UT(I,1,1)*UT(I,2,1)+UT(I,1,2)*UT(I,2,2))
     .       *(UT(k,1,1)*UT(k,2,1)+UT(k,1,2)*UT(k,2,2))
     .      +(UT(I,1,1)*UT(I,2,2)-UT(I,1,2)*UT(I,2,1))
     .       *(UT(k,1,1)*UT(k,2,2)-UT(k,1,2)*UT(k,2,1))))
     .                *NMA0(MST2(k),QSTSB)
     .       +(Ybq**2*(UT(I,1,1)**2+UT(I,1,2)**2)
     .                             *(UB(k,2,1)**2+UB(k,2,2)**2)
     .        +Ytq**2*(UT(I,2,1)**2+UT(I,2,2)**2)
     .                             *(UB(k,1,1)**2+UB(k,1,2)**2))
     .                *NMA0(MSB2(k),QSTSB)
      ENDDO
      
*       2) Higgs/sfermion Quartic couplings
*  - neutral
      DO k=1,6
       DPIST(I)=DPIST(I)+Ytq**2/2d0*
     .    (UT(I,1,1)**2+UT(I,1,2)**2+UT(I,2,1)**2+UT(I,2,2)**2)
     .     *(XHT(k,1)**2+XHT(k,4)**2)*NMA0(MHT2(k),QSTSB)
      ENDDO

*  - charged
      DO k=1,2
       DPIST(I)=DPIST(I)+
     .   (Ybq**2*(UT(I,1,1)**2+UT(I,1,2)**2)*XC(k,2)**2
     .   +Ytq**2*(UT(I,2,1)**2+UT(I,2,2)**2)*XC(k,1)**2)
     .                *NMA0(MHTC(k),QSTSB)
      ENDDO
      
*       3) Higgs/sfermion loop
*  - neutral
      DO j=1,6
      DO k=1,2
c       COH0SFL(j,k,1)=dsqrt(2d0)*mtopq*XHT(j,1)*UT(k,1,1)
c     .        +(AT*(DDCOS(phiAT)*XHT(j,1)-DDSIN(phiAT)*XHT(j,4))
c     .          -muq*(DDCOS(phi01)*XHT(j,2)-DDSIN(phi01)*XHT(j,5))
c     .          -l*vdq*(DDCOS(phi01)*XHT(j,3)-DDSIN(phi01)*XHT(j,6)))
c     .                             *UT(k,2,1)/dsqrt(2d0)
c     .        -(AT*(DDSIN(phiAT)*XHT(j,1)+DDCOS(phiAT)*XHT(j,4))
c     .          +muq*(DDSIN(phi01)*XHT(j,2)+DDCOS(phi01)*XHT(j,5))
c     .          +l*vdq*(DDSIN(phi01)*XHT(j,3)+DDCOS(phi01)*XHT(j,6)))
c     .                             *UT(k,2,2)/dsqrt(2d0)
c       COH0SFL(j,k,2)=dsqrt(2d0)*mtopq*XHT(j,1)*UT(k,1,2)
c     .        +(AT*(DDCOS(phiAT)*XHT(j,1)-DDSIN(phiAT)*XHT(j,4))
c     .          -muq*(DDCOS(phi01)*XHT(j,2)-DDSIN(phi01)*XHT(j,5))
c     .          -l*vdq*(DDCOS(phi01)*XHT(j,3)-DDSIN(phi01)*XHT(j,6)))
c     .                             *UT(k,2,2)/dsqrt(2d0)
c     .        +(AT*(DDSIN(phiAT)*XHT(j,1)+DDCOS(phiAT)*XHT(j,4))
c     .          +muq*(DDSIN(phi01)*XHT(j,2)+DDCOS(phi01)*XHT(j,5))
c     .          +l*vdq*(DDSIN(phi01)*XHT(j,3)+DDCOS(phi01)*XHT(j,6)))
c     .                             *UT(k,2,1)/dsqrt(2d0)

c       COH0SFR(j,k,1)=dsqrt(2d0)*mtopq*XHT(j,1)*UT(k,2,1)
c     .        +(AT*(DDCOS(phiAT)*XHT(j,1)-DDSIN(phiAT)*XHT(j,4))
c     .          -muq*(DDCOS(phi01)*XHT(j,2)-DDSIN(phi01)*XHT(j,5))
c     .          -l*vdq*(DDCOS(phi01)*XHT(j,3)-DDSIN(phi01)*XHT(j,6)))
c     .                             *UT(k,1,1)/dsqrt(2d0)
c     .        +(AT*(DDSIN(phiAT)*XHT(j,1)+DDCOS(phiAT)*XHT(j,4))
c     .          +muq*(DDSIN(phi01)*XHT(j,2)+DDCOS(phi01)*XHT(j,5))
c     .          +l*vdq*(DDSIN(phi01)*XHT(j,3)+DDCOS(phi01)*XHT(j,6)))
c     .                             *UT(k,1,2)/dsqrt(2d0)
c       COH0SFR(j,k,2)=-dsqrt(2d0)*mtopq*XHT(j,1)*UT(k,2,2)
c     .        -(AT*(DDCOS(phiAT)*XHT(j,1)-DDSIN(phiAT)*XHT(j,4))
c     .          -muq*(DDCOS(phi01)*XHT(j,2)-DDSIN(phi01)*XHT(j,5))
c     .          -l*vdq*(DDCOS(phi01)*XHT(j,3)-DDSIN(phi01)*XHT(j,6)))
c     .                             *UT(k,1,2)/dsqrt(2d0)
c     .        +(AT*(DDSIN(phiAT)*XHT(j,1)+DDCOS(phiAT)*XHT(j,4))
c     .          +muq*(DDSIN(phi01)*XHT(j,2)+DDCOS(phi01)*XHT(j,5))
c     .          +l*vdq*(DDSIN(phi01)*XHT(j,3)+DDCOS(phi01)*XHT(j,6)))
c     .                             *UT(k,2,1)/dsqrt(2d0)
c       DPIST(I)=DPIST(I)+Ytq**2*(
c     .           (UT(I,1,1)**2+UT(I,1,2)**2)
c     .                   *(COH0SFL(j,k,1)**2+COH0SFL(j,k,2)**2)
c     .          +(UT(I,2,1)**2+UT(I,2,2)**2)
c     .                   *(COH0SFR(j,k,1)**2+COH0SFR(j,k,2)**2)
c     .          +2d0*(UT(I,1,1)*UT(I,2,1)+UT(I,1,2)*UT(I,2,2))
c     .                   *(COH0SFL(j,k,1)*COH0SFR(j,k,1)
c     .                      -COH0SFL(j,k,2)*COH0SFR(j,k,2))
c     .          +2d0*(UT(I,1,2)*UT(I,2,1)-UT(I,1,1)*UT(I,2,2))
c     .                   *(COH0SFL(j,k,2)*COH0SFR(j,k,1)
c     .                      +COH0SFL(j,k,1)*COH0SFR(j,k,2)))
c     .                       *NMB0(MST2(i),MHT2(j),MST2(k),QSTSB)

      COH0SFL(j,k,1)=dsqrt(2d0)*Ytq**2*vuq*XHT(j,1)
     .               *(UT(I,1,1)*UT(k,1,1)+UT(I,1,2)*UT(k,1,2))
     . +dsqrt(2d0)*Ytq**2*vuq*XHT(j,1)
     .               *(UT(I,2,1)*UT(k,2,1)+UT(I,2,2)*UT(k,2,2))
     . +Ytq/dsqrt(2d0)*(AT*(DDCOS(PhiAT)*XHT(j,1)
     .                  -DDSIN(PhiAT)*XHT(j,4))
     .                  -DDCOS(Phi01)*(muq*XHT(j,2)+l*vdq*XHT(j,3))
     .                  +DDSIN(Phi01)*(muq*XHT(j,5)+l*vdq*XHT(j,6)))
     .               *(UT(I,1,1)*UT(k,2,1)+UT(I,2,1)*UT(k,1,1)
     .                +UT(I,1,2)*UT(k,2,2)+UT(I,2,2)*UT(k,1,2))
     . +Ytq/dsqrt(2d0)*(AT*(DDSIN(PhiAT)*XHT(j,1)
     .                  +DDCOS(PhiAT)*XHT(j,4))
     .                  +DDSIN(Phi01)*(muq*XHT(j,2)+l*vdq*XHT(j,3))
     .                  +DDCOS(Phi01)*(muq*XHT(j,5)+l*vdq*XHT(j,6)))
     .               *(UT(I,2,1)*UT(k,1,2)-UT(I,2,2)*UT(k,1,1)
     .                +UT(I,1,2)*UT(k,2,1)-UT(I,1,1)*UT(k,2,2))
      COH0SFL(j,k,2)=dsqrt(2d0)*Ytq**2*vuq*XHT(j,1)
     .               *(UT(i,1,2)*UT(k,1,1)-UT(i,1,1)*UT(k,1,2))
     . +dsqrt(2d0)*Ytq**2*vuq*XHT(j,1)
     .               *(UT(i,2,2)*UT(k,2,1)-UT(i,2,1)*UT(k,2,2))
     . +Ytq/dsqrt(2d0)*(AT*(DDCOS(PhiAT)*XHT(j,1)
     .                  -DDSIN(PhiAT)*XHT(j,4))
     .                  -DDCOS(Phi01)*(muq*XHT(j,2)+l*vdq*XHT(j,3))
     .                  +DDSIN(Phi01)*(muq*XHT(j,5)+l*vdq*XHT(j,6)))
     .               *(UT(i,2,2)*UT(k,1,1)-UT(i,2,1)*UT(k,1,2)
     .                +UT(i,1,2)*UT(k,2,1)-UT(i,1,1)*UT(k,2,2))
     . +Ytq/dsqrt(2d0)*(AT*(DDSIN(PhiAT)*XHT(j,1)
     .                  +DDCOS(PhiAT)*XHT(j,4))
     .                  +DDSIN(Phi01)*(muq*XHT(j,2)+l*vdq*XHT(j,3))
     .                  +DDCOS(Phi01)*(muq*XHT(j,5)+l*vdq*XHT(j,6)))
     .               *(UT(i,2,1)*UT(k,1,1)-UT(i,1,1)*UT(k,2,1)
     .                +UT(i,2,2)*UT(k,1,2)-UT(i,1,2)*UT(k,2,2))
       DPIST(I)=DPIST(I)+(COH0SFL(j,k,1)**2+COH0SFL(j,k,2)**2)
     .                       *NMB0(MST2(i),MHT2(j),MST2(k),QSTSB)

      ENDDO
      ENDDO

*  - charged
      DO j=1,2
      DO k=1,2
       COHCSFL(j,k,1)=-(Ytq*mtopq*XC(j,1)+Ybq*mbotq*XC(j,2))
     .                                                *UB(k,1,1)
     .  -Ybq*(Ab*DDCOS(phiAb)*XC(j,2)+muq*DDCOS(phi01)*XC(j,1))
     .                                                *UB(k,2,1)
     .  +Ybq*(Ab*DDSIN(phiAb)*XC(j,2)-muq*DDSIN(phi01)*XC(j,1))
     .                                                *UB(k,2,2)
       COHCSFL(j,k,2)=(Ytq*mtopq*XC(j,1)+Ybq*mbotq*XC(j,2))
     .                                                *UB(k,1,2)
     .  +Ybq*(Ab*DDCOS(phiAb)*XC(j,2)+muq*DDCOS(phi01)*XC(j,1))
     .                                                *UB(k,2,2)
     .  +Ybq*(Ab*DDSIN(phiAb)*XC(j,2)-muq*DDSIN(phi01)*XC(j,1))
     .                                                *UB(k,2,1)

       COHCSFR(j,k,1)=-Ytq*Ybq*(vuq*XC(j,2)+vdq*XC(j,1))*UB(k,2,1)
     .  -Ytq*(At*DDCOS(phiAt)*XC(j,1)+muq*DDCOS(phi01)*XC(j,2))
     .                                                *UB(k,1,1)
     .  -Ytq*(At*DDSIN(phiAt)*XC(j,1)-muq*DDSIN(phi01)*XC(j,2))
     .                                                *UB(k,1,2)
       COHCSFR(j,k,2)=Ytq*Ybq*(vuq*XC(j,2)+vdq*XC(j,1))*UB(k,2,2)
     .  +Ytq*(At*DDCOS(phiAt)*XC(j,1)+muq*DDCOS(phi01)*XC(j,2))
     .                                                *UB(k,1,2)
     .  -Ytq*(At*DDSIN(phiAt)*XC(j,1)-muq*DDSIN(phi01)*XC(j,2))
     .                                                *UB(k,1,1)
     
       DPIST(I)=DPIST(I)+((UT(I,1,1)**2+UT(I,1,2)**2)
     .                   *(COHCSFL(j,k,1)**2+COHCSFL(j,k,2)**2)
     .          +(UT(I,2,1)**2+UT(I,2,2)**2)
     .                   *(COHCSFR(j,k,1)**2+COHCSFR(j,k,2)**2)
     .          +2d0*(UT(I,1,1)*UT(I,2,1)+UT(I,1,2)*UT(I,2,2))
     .                   *(COHCSFL(j,k,1)*COHCSFR(j,k,1)
     .                      +COHCSFL(j,k,2)*COHCSFR(j,k,2))
     .          +2d0*(UT(I,1,2)*UT(I,2,1)-UT(I,1,1)*UT(I,2,2))
     .                   *(COHCSFL(j,k,2)*COHCSFR(j,k,1)
     .                      +COHCSFL(j,k,1)*COHCSFR(j,k,2)))
     .                       *NMB0(MST2(i),MHTC(j),MSB2(k),QSTSB)

c       COHCSFL(j,k,1)=-(Ytq*mtopq*XC(j,1)+Ybq*mbotq*XC(j,2))
c     .                   *(UT(i,1,1)*UB(k,1,1)+UT(i,1,2)*UB(k,1,2))
c     .  -Ytq*Ybq*(vuq*XC(j,2)+vdq*XC(j,1))
c     .                   *(UT(i,2,1)*UB(k,2,1)+UT(i,2,2)*UB(k,2,2))
c     .  -Ytq*(At*DDCOS(phiAt)*XC(j,1)+muq*DDCOS(phi01)*XC(j,2))
c     .                   *(UT(i,2,1)*UB(k,1,1)+UT(i,2,2)*UB(k,1,2))
c     .  -Ytq*(At*DDSIN(phiAt)*XC(j,1)-muq*DDSIN(phi01)*XC(j,2))
c     .                   *(UT(i,2,1)*UB(k,1,2)-UT(i,2,2)*UB(k,1,1))
c     .  -Ybq*(Ab*DDCOS(phiAb)*XC(j,2)+muq*DDCOS(phi01)*XC(j,1))
c     .                   *(UT(i,1,1)*UB(k,2,1)+UT(i,1,2)*UB(k,2,2))
c     .  +Ybq*(Ab*DDSIN(phiAb)*XC(j,2)-muq*DDSIN(phi01)*XC(j,1))
c     .                   *(UT(i,1,1)*UB(k,2,2)-UT(i,1,2)*UB(k,2,1))
c       COHCSFL(j,k,2)=-(Ytq*mtopq*XC(j,1)+Ybq*mbotq*XC(j,2))
c     .                   *(UT(i,1,2)*UB(k,1,1)-UT(i,1,1)*UB(k,1,2))
c     .  -Ytq*Ybq*(vuq*XC(j,2)+vdq*XC(j,1))
c     .                   *(UT(i,2,2)*UB(k,2,1)-UT(i,2,1)*UB(k,2,2))
c     .  +Ytq*(At*DDCOS(phiAt)*XC(j,1)+muq*DDCOS(phi01)*XC(j,2))
c     .                   *(UT(i,2,1)*UB(k,1,2)-UT(i,2,2)*UB(k,1,1))
c     .  -Ytq*(At*DDSIN(phiAt)*XC(j,1)-muq*DDSIN(phi01)*XC(j,2))
c     .                   *(UT(i,2,1)*UB(k,1,1)+UT(i,2,2)*UB(k,1,2))
c     .  +Ybq*(Ab*DDCOS(phiAb)*XC(j,2)+muq*DDCOS(phi01)*XC(j,1))
c     .                   *(UT(i,1,1)*UB(k,2,2)-UT(i,1,2)*UB(k,2,1))
c     .  +Ybq*(Ab*DDSIN(phiAb)*XC(j,2)-muq*DDSIN(phi01)*XC(j,1))
c     .                   *(UT(i,1,1)*UB(k,2,1)+UT(i,1,2)*UB(k,2,2))
c       DPIST(I)=DPIST(I)+(COHCSFL(j,k,1)**2+COHCSFL(j,k,2)**2)
c     .                       *NMB0(MST2(i),MHTC(j),MSB2(k),QSTSB)

      ENDDO
      ENDDO
      
*       4) SM fermion / higgsino loop

       DPIST(I)=DPIST(I)+Ytq**2*
     .  (UT(I,1,1)**2+UT(I,1,2)**2+UT(I,2,1)**2+UT(I,2,2)**2)
     . *((MST2(i)-mur**2-mtopq**2)*NMB0(MST2(i),mur**2,mtopq**2,QSTSB)
     .                       -NMA0(mur**2,QSTSB)-NMA0(mtopq**2,QSTSB))
     . +(Ybq**2*(UT(I,1,1)**2+UT(I,1,2)**2)
     .                            +Ytq**2*(UT(I,2,1)**2+UT(I,2,2)**2))
     . *((MST2(i)-mur**2-mbotq**2)*NMB0(MST2(i),mur**2,mbotq**2,QSTSB)
     .                       -NMA0(mur**2,QSTSB)-NMA0(mbotq**2,QSTSB))
     . -4d0*(UT(I,1,1)*UT(I,2,1)+UT(I,1,2)*UT(I,2,2))
     .     *Ytq*Ybq*mur*mbotq*NMB0(MST2(i),mur**2,mbotq**2,QSTSB)
     
       DPIST(I)=DPIST(I)/16d0/Pi**2
      ENDDO


*      III- Corrections to the sbottom masses

      DO I=1,2
       DPISB(I)=0d0

*       1) Sfermion Quartic couplings
      DO k=1,2
       DPISB(I)=DPISB(I)+Ybq**2*(
     .  (UB(I,1,1)**2+UB(I,1,2)**2)*(UB(k,2,1)**2+UB(k,2,2)**2)
     . +(UB(I,2,1)**2+UB(I,2,2)**2)*(UB(k,1,1)**2+UB(k,1,2)**2)
     . +6d0*((UB(I,1,1)*UB(I,2,1)+UB(I,1,2)*UB(I,2,2))
     .       *(UB(k,1,1)*UB(k,2,1)+UB(k,1,2)*UB(k,2,2))
     .      +(UB(I,1,1)*UB(I,2,2)+UB(I,1,2)*UB(I,2,1))
     .       *(UB(k,1,1)*UB(k,2,2)+UB(k,1,2)*UB(k,2,1))))
     .                *NMA0(MSB2(k),QSTSB)
     .       +(Ytq**2*(UB(I,1,1)**2+UB(I,1,2)**2)
     .                             *(UT(k,2,1)**2+UT(k,2,2)**2)
     .        +Ybq**2*(UB(I,2,1)**2+UB(I,2,2)**2)
     .                             *(UT(k,1,1)**2+UT(k,1,2)**2))
     .                *NMA0(MST2(k),QSTSB)
      ENDDO
      
*       2) Higgs/sfermion Quartic couplings
*  - neutral
      DO k=1,6
       DPISB(I)=DPISB(I)+Ybq**2/2d0
     .  *(UB(I,1,1)**2+UB(I,1,2)**2+UB(I,2,1)**2+UB(I,2,2)**2)
     .           *(XHT(k,2)**2+XHT(k,5)**2)*NMA0(MHT2(k),QSTSB)
      ENDDO

*  - charged
      DO k=1,2
       DPISB(I)=DPISB(I)+
     .   (Ytq**2*(UB(I,1,1)**2+UB(I,1,2)**2)*XC(k,1)**2
     .   +Ybq**2*(UB(I,2,1)**2+UB(I,2,2)**2)*XC(k,2)**2)
     .                *NMA0(MHTC(k),QSTSB)
      ENDDO
      
*       3) Higgs/sfermion loop
*  - neutral
      DO j=1,6
      DO k=1,2
       COH0SFL(j,k,1)=dsqrt(2d0)*mbotq*XHT(j,2)*UB(k,1,1)
     .        +(AB*(DDCOS(phiAB)*XHT(j,2)-DDSIN(phiAB)*XHT(i,5))
     .          -muq*(DDCOS(phi01)*XHT(j,1)-DDSIN(phi01)*XHT(j,4))
     .          -l*vuq*(DDCOS(phi01)*XHT(j,3)-DDSIN(phi01)*XHT(j,6)))
     .                             *UB(k,2,1)/dsqrt(2d0)
     .        -(AB*(DDSIN(phiAB)*XHT(j,2)+DDCOS(phiAB)*XHT(i,5))
     .          +muq*(DDSIN(phi01)*XHT(j,1)+DDCOS(phi01)*XHT(j,4))
     .          +l*vuq*(DDSIN(phi01)*XHT(j,3)+DDCOS(phi01)*XHT(j,6)))
     .                             *UB(k,2,2)/dsqrt(2d0)
       COH0SFL(j,k,2)=dsqrt(2d0)*mbotq*XHT(j,2)*UB(k,1,2)
     .        +(AB*(DDCOS(phiAB)*XHT(j,2)-DDSIN(phiAB)*XHT(i,5))
     .          -muq*(DDCOS(phi01)*XHT(j,1)-DDSIN(phi01)*XHT(j,4))
     .          -l*vuq*(DDCOS(phi01)*XHT(j,3)-DDSIN(phi01)*XHT(j,6)))
     .                             *UB(k,2,2)/dsqrt(2d0)
     .        +(AB*(DDSIN(phiAB)*XHT(j,2)+DDCOS(phiAB)*XHT(i,5))
     .          +muq*(DDSIN(phi01)*XHT(j,1)+DDCOS(phi01)*XHT(j,4))
     .          +l*vuq*(DDSIN(phi01)*XHT(j,3)+DDCOS(phi01)*XHT(j,6)))
     .                             *UB(k,2,1)/dsqrt(2d0)

       COH0SFR(j,k,1)=dsqrt(2d0)*mbotq*XHT(j,2)*UB(k,2,1)
     .        +(AB*(DDCOS(phiAB)*XHT(j,2)-DDSIN(phiAB)*XHT(i,5))
     .          -muq*(DDCOS(phi01)*XHT(j,1)-DDSIN(phi01)*XHT(j,4))
     .          -l*vuq*(DDCOS(phi01)*XHT(j,3)-DDSIN(phi01)*XHT(j,6)))
     .                             *UB(k,1,1)/dsqrt(2d0)
     .        +(AB*(DDSIN(phiAB)*XHT(j,2)+DDCOS(phiAB)*XHT(i,5))
     .          +muq*(DDSIN(phi01)*XHT(j,1)+DDCOS(phi01)*XHT(j,4))
     .          +l*vuq*(DDSIN(phi01)*XHT(j,3)+DDCOS(phi01)*XHT(j,6)))
     .                             *UB(k,1,2)/dsqrt(2d0)
       COH0SFR(j,k,2)=-dsqrt(2d0)*mbotq*XHT(j,2)*UB(k,2,2)
     .        -(AB*(DDCOS(phiAB)*XHT(j,2)-DDSIN(phiAB)*XHT(i,5))
     .          -muq*(DDCOS(phi01)*XHT(j,1)-DDSIN(phi01)*XHT(j,4))
     .          -l*vuq*(DDCOS(phi01)*XHT(j,3)-DDSIN(phi01)*XHT(j,6)))
     .                             *UB(k,1,2)/dsqrt(2d0)
     .        +(AB*(DDSIN(phiAB)*XHT(j,2)+DDCOS(phiAB)*XHT(i,5))
     .          +muq*(DDSIN(phi01)*XHT(j,1)+DDCOS(phi01)*XHT(j,4))
     .          +l*vuq*(DDSIN(phi01)*XHT(j,3)+DDCOS(phi01)*XHT(j,6)))
     .                             *UB(k,2,1)/dsqrt(2d0)
     
c       DPISB(I)=DPISB(I)+Ybq**2*(
c     .           (UB(I,1,1)**2+UB(I,1,2)**2)
c     .                   *(COH0SFL(j,k,1)**2+COH0SFL(j,k,2)**2)
c     .          +(UB(I,2,1)**2+UB(I,2,2)**2)
c     .                   *(COH0SFR(j,k,1)**2+COH0SFR(j,k,2)**2)
c     .          +2d0*(UB(I,1,1)*UB(I,2,1)+UB(I,1,2)*UB(I,2,2))
c     .                   *(COH0SFL(j,k,1)*COH0SFR(j,k,1)
c     .                      -COH0SFL(j,k,2)*COH0SFR(j,k,2))
c     .          +2d0*(UB(I,1,2)*UB(I,2,1)-UB(I,1,1)*UB(I,2,2))
c     .                   *(COH0SFL(j,k,2)*COH0SFR(j,k,1)
c     .                      +COH0SFL(j,k,1)*COH0SFR(j,k,2)))
c     .                       *NMB0(MSB2(i),MHT2(j),MSB2(k),QSTSB)

       COH0SFL(j,k,1)=dsqrt(2d0)*Ybq**2*vdq*XHT(j,2)
     .          *(UB(i,1,1)*UB(k,1,1)+UB(i,1,2)*UB(k,1,2)
     .           +UB(i,2,1)*UB(k,2,1)+UB(i,2,2)*UB(k,2,2))
     . +Ybq/dsqrt(2d0)*(AB*(DDCOS(phiAB)*XHT(j,2)-DDSIN(phiAB)*XHT(i,5))
     .             -muq*(DDCOS(phi01)*XHT(j,1)-DDSIN(phi01)*XHT(j,4))
     .             -l*vuq*(DDCOS(phi01)*XHT(j,3)-DDSIN(phi01)*XHT(j,6)))
     .           *(UB(i,2,1)*UB(k,1,1)+UB(i,2,2)*UB(k,1,2)
     .            +UB(i,1,1)*UB(k,2,1)+UB(i,1,2)*UB(k,2,2))
     . +Ybq/dsqrt(2d0)*(AB*(DDSIN(phiAB)*XHT(j,2)+DDCOS(phiAB)*XHT(i,5))
     .             +muq*(DDSIN(phi01)*XHT(j,1)+DDCOS(phi01)*XHT(j,4))
     .             +l*vuq*(DDSIN(phi01)*XHT(j,3)+DDCOS(phi01)*XHT(j,6)))
     .           *(UB(i,2,1)*UB(k,1,2)-UB(i,2,2)*UB(k,1,1)
     .            +UB(i,1,2)*UB(k,2,1)-UB(i,1,1)*UB(k,2,2))
       COH0SFL(j,k,2)=dsqrt(2d0)*Ybq**2*vdq*XHT(j,2)
     .          *(UB(i,1,2)*UB(k,1,1)-UB(i,1,1)*UB(k,1,2)
     .           +UB(i,2,1)*UB(k,2,1)+UB(i,2,2)*UB(k,2,2))
     . +Ybq/dsqrt(2d0)*(AB*(DDCOS(phiAB)*XHT(j,2)-DDSIN(phiAB)*XHT(i,5))
     .             -muq*(DDCOS(phi01)*XHT(j,1)-DDSIN(phi01)*XHT(j,4))
     .             -l*vuq*(DDCOS(phi01)*XHT(j,3)-DDSIN(phi01)*XHT(j,6)))
     .           *(UB(i,2,2)*UB(k,1,1)-UB(i,2,1)*UB(k,1,2)
     .            +UB(i,1,2)*UB(k,2,1)+UB(i,1,1)*UB(k,2,2))
     . +Ybq/dsqrt(2d0)*(AB*(DDSIN(phiAB)*XHT(j,2)+DDCOS(phiAB)*XHT(i,5))
     .             +muq*(DDSIN(phi01)*XHT(j,1)+DDCOS(phi01)*XHT(j,4))
     .             +l*vuq*(DDSIN(phi01)*XHT(j,3)+DDCOS(phi01)*XHT(j,6)))
     .           *(UB(i,2,1)*UB(k,1,1)+UB(i,2,2)*UB(k,1,2)
     .            -UB(i,1,1)*UB(k,2,1)-UB(i,1,2)*UB(k,2,2))
       DPISB(I)=DPISB(I)+(COH0SFL(j,k,1)**2+COH0SFL(j,k,2)**2)
     .                       *NMB0(MSB2(i),MHT2(j),MSB2(k),QSTSB)
      ENDDO
      ENDDO

*  - charged
      DO j=1,2
      DO k=1,2
       COHCSFL(j,k,1)=-(Ytq*mtopq*XC(j,1)+Ybq*mbotq*XC(j,2))
     .                                                *UT(k,1,1)
     .  -Ytq*(At*DDCOS(phiAt)*XC(j,1)+muq*DDCOS(phi01)*XC(j,2))
     .                                                *UT(k,2,1)
     .  +Ytq*(At*DDSIN(phiAt)*XC(j,1)-muq*DDSIN(phi01)*XC(j,2))
     .                                                *UT(k,2,2)
       COHCSFL(j,k,2)=-(Ytq*mtopq*XC(j,1)+Ybq*mbotq*XC(j,2))
     .                                                *UT(k,1,2)
     .  -Ytq*(At*DDCOS(phiAt)*XC(j,1)+muq*DDCOS(phi01)*XC(j,2))
     .                                                *UT(k,2,2)
     .  -Ytq*(At*DDSIN(phiAt)*XC(j,1)-muq*DDSIN(phi01)*XC(j,2))
     .                                                *UT(k,2,1)

       COHCSFR(j,k,1)=-Ytq*Ybq*(vuq*XC(j,2)+vdq*XC(j,1))*UT(k,2,1)
     .  -Ybq*(Ab*DDCOS(phiAb)*XC(j,2)+muq*DDCOS(phi01)*XC(j,1))
     .                                                *UT(k,1,1)
     .  -Ybq*(Ab*DDSIN(phiAb)*XC(j,2)-muq*DDSIN(phi01)*XC(j,1))
     .                                                *UT(k,1,2)
       COHCSFR(j,k,2)=-Ytq*Ybq*(vuq*XC(j,2)+vdq*XC(j,1))*UT(k,2,2)
     .  -Ybq*(Ab*DDCOS(phiAb)*XC(j,2)+muq*DDCOS(phi01)*XC(j,1))
     .                                                *UT(k,1,2)
     .  +Ybq*(Ab*DDSIN(phiAb)*XC(j,2)-muq*DDSIN(phi01)*XC(j,1))
     .                                                *UT(k,1,1)
     
       DPISB(I)=DPISB(I)+((UB(I,1,1)**2+UB(I,1,2)**2)
     .                   *(COHCSFL(j,k,1)**2+COHCSFL(j,k,2)**2)
     .          +(UB(I,2,1)**2+UB(I,2,2)**2)
     .                   *(COHCSFR(j,k,1)**2+COHCSFR(j,k,2)**2)
     .          +2d0*(UB(I,1,1)*UB(I,2,1)+UB(I,1,2)*UB(I,2,2))
     .                   *(COHCSFL(j,k,1)*COHCSFR(j,k,1)
     .                      +COHCSFL(j,k,2)*COHCSFR(j,k,2))
     .          +2d0*(UB(I,1,2)*UB(I,2,1)-UB(I,1,1)*UB(I,2,2))
     .                   *(COHCSFL(j,k,2)*COHCSFR(j,k,1)
     .                      -COHCSFL(j,k,1)*COHCSFR(j,k,2)))
     .                       *NMB0(MSB2(i),MHTC(j),MST2(k),QSTSB)
      ENDDO
      ENDDO
      
*       4) SM fermion / higgsino loop

       DPISB(I)=DPISB(I)+Ybq**2*
     .  (UB(I,1,1)**2+UB(I,1,2)**2+UB(I,2,1)**2+UB(I,2,2)**2)
     . *((MSB2(i)-mur**2-mbotq**2)*NMB0(MSB2(i),mur**2,mbotq**2,QSTSB)
     .                       -NMA0(mur**2,QSTSB)-NMA0(mbotq**2,QSTSB))
     . +(Ytq**2*(UB(I,1,1)**2+UB(I,1,2)**2)
     .                            +Ybq**2*(UB(I,2,1)**2+UB(I,2,2)**2))
     . *((MSB2(i)-mur**2-mtopq**2)*NMB0(MSB2(i),mur**2,mtopq**2,QSTSB)
     .                       -NMA0(mur**2,QSTSB)-NMA0(mtopq**2,QSTSB))
     . -4d0*(UB(I,1,1)*UB(I,2,1)+UB(I,1,2)*UB(I,2,2))
     .     *Ytq*Ybq*mur*mtopq*NMB0(MSB2(i),mur**2,mtopq**2,QSTSB)
     
       DPISB(I)=DPISB(I)/16d0/Pi**2
      ENDDO

      END
