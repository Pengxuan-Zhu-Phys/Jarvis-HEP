      SUBROUTINE CHECKMIN_CPV(PAR,PROB)

**********************************************************************
c The effective potential (up to a constant term) is evaluated at the 
c electroweak minimum, VEW. It is then compared to the value of the
c potential at several other points:
c  * VZ3: for <S>=<Hu>=<Hd>=0
c  * Vgauge: for <Hu>=<Hd>=0, <S>=mu/l
c  * Vvu: for <S>=<Hd>=0, <Hu>=vu
c  * Vvd: for <S>=<Hu>=0, <Hd>=vd
c  * Vphi: for <S>=mu/l.exp(i phiS), <Hu.Hd>=vu.vd.exp(i phiH)
c where phiS and phiH are scanned over on a Nphi*Nphi grid.
c Should one of these minima prove deeper than VEW, PROB(28) is. 
c switched on
c The quadratic terms of the potential mHu2 and mHd2 are computed and
c it is checked whether they lie too far from the SUSY scale, in which 
c case PROB(29) neq 0.
**********************************************************************

      IMPLICIT NONE

      INTEGER I,J,Nphi,CFLAG(5)

      DOUBLE PRECISION Pi,aux,Taux,Daux,Ytau
      DOUBLE PRECISION MW2,MZ2
      DOUBLE PRECISION VEW,VZ3,Vgauge,Vvu,Vvd,Vphi
      DOUBLE PRECISION phiS,phiH
      DOUBLE PRECISION MSTaux1,MSTaux2,MSBaux1,MSBaux2,MSLaux1,MSLaux2,
     . MSNTaux,MSUaux1,MSUaux2,MSDaux1,MSDaux2,MSEaux1,MSEaux2,MSNEaux
      DOUBLE PRECISION MHD1aux,MHD2aux,MS1aux,MS2aux
      DOUBLE PRECISION MH2(6,6),VALPH(6),VECPH(6,6),Mch1aux,Mch2aux
      DOUBLE PRECISION MN2(10,10),VALPN(10),VECPN(10,10)
      DOUBLE PRECISION Aroot,Broot,Croot,Droot,Proot,Qroot,Delt,root1,
     . root2,root3

      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION Q2
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION l,k,Alcos1,Akcos2,muq,nuq
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      DOUBLE PRECISION mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      DOUBLE PRECISION Ytq,Ybq,MTOPQ,MBOTQ
      DOUBLE PRECISION MST2(2),UT(2,2,2),MSB2(2),UB(2,2,2),MSL2(2),
     . UTAU(2,2,2),MSNT2
      DOUBLE PRECISION MSU2(2),MSD2(2),MSE2(2),MSNE2,MSMU2(2),
     . UMU(2,2,2)
      DOUBLE PRECISION MSQ3,MSU3,MSD3,ATP,ABP
      DOUBLE PRECISION MSQ1,MSU1,MSD1
      DOUBLE PRECISION MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      DOUBLE PRECISION MCH2(2),U(2,2,2),V(2,2,2)
      DOUBLE PRECISION MNEU(5),NEU(5,5,2)
      DOUBLE PRECISION mur,M1r,M2r,msi
      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION PAR(*),PROB(*)
      DOUBLE PRECISION MHuS,MHdS,MSS
      DOUBLE PRECISION IAL,IAK,IXIS
      DOUBLE PRECISION phiF,phiP,phi3,phiSs,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION RXIS,RXIF,IAlQ2,IAkQ2,IXISQ2,MHuQ2,MHdQ2,MSQ2
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/STSBSCALE/QSTSB
      COMMON/RENSCALE/Q2
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/QPAR/l,k,Alcos1,Akcos2,muq,NUQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/QEXT/XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      COMMON/SMFERM/mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      COMMON/QQUARK/Ytq,Ybq,MTOPQ,MBOTQ
      COMMON/SFERM3SPEC/MST2,UT,MSB2,UB,MSL2,UTAU,MSNT2
      COMMON/SFERM1SPEC/MSU2,MSD2,MSE2,MSNE2,MSMU2,UMU
      COMMON/RADCOR2/MSQ3,MSU3,MSD3,ATP,ABP
      COMMON/SQUPAR/MSQ1,MSU1,MSD1
      COMMON/SLEPPAR/MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      COMMON/CHASPEC/MCH2,U,V
      COMMON/NEUSPEC/MNEU,NEU
      COMMON/GAUGINOPAR/mur,M1r,M2r,msi
      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/MH2TREE/MHuS,MHdS,MSS
      COMMON/IMALAK/IAL,IAK,IXIS
      COMMON/Z3VAUX/phiF,phiP,phi3,phiSs,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/REXI/RXIS,RXIF,IAlQ2,IAkQ2,IXISQ2,MHuQ2,MHdQ2,MSQ2
      COMMON/CFLAG/CFLAG

**********************************************************************

      Pi=4d0*datan(1d0)

c      Running to Q2

      Ytau=mtau/vdq

      IAlQ2=IAl-(4d0*l**2*IAl
     .    +2d0*k**2*(IAk*DDCOS(phi0)+Akcos2*DDSIN(phi0))
     .    +3d0*Ytq**2*PAR(12)*DDSIN(phiAt+phi01)
     .    +3d0*Ybq**2*PAR(13)*DDSIN(phiAb+phi01)
     .    +Ytau**2*PAR(14)*DDSIN(phiAtau+phi01)
     .    +g1q*PAR(20)*DDSIN(phiM1+phi01)
     .    +3d0*g2q*PAR(21)*DDSIN(phiM2+phi01))*DLOG(QSTSB/Q2)/16d0/Pi**2

       If(k.ne.0d0)then
      IAkQ2=IAk-6d0*(l**2*(IAl*DDCOS(phi0)-PAR(5)*DDSIN(phi0))
     .     +k**2*IAk)*DLOG(QSTSB/Q2)/16d0/Pi**2

      IXISQ2=XIS*DDSIN(PhiSs)
      else
      IAkQ2=0d0

      IXISQ2=IXIS
     .        *(1d0-(L**2+K**2)/16d0/Pi**2*DLOG(QSTSB/Q2))
     .  -(2d0*XIF*(L**2*(IAl*DDCOS(phiF-phi01)+PAR(5)*DDSIN(phiF-phi01))
     .           +K**2*(IAk*DDCOS(phiF-phi02)+PAR(6)*DDSIN(phiF-phi02)))
     .   +2d0*L*M3H*(IAl*DDCOS(phi3)-PAR(5)*DDSIN(phi3)
     .                             +MUP*DDSIN(phi3+phiP-phi01))
     .   +K*MSP*(IAk*DDCOS(phiSP)-PAR(6)*DDSIN(phiSP)
     .                             +MUP*DDSIN(phip+phiSP-phi02))
     .   +2d0*k*mup*MSS*DDSIN(phi02-phiP))/16d0/Pi**2*DLOG(QSTSB/Q2)
      endif

      MHuQ2=MHuS-(3d0*Ytq**2*(MSQ3+MSU3+MHuS+ATP**2)
     .           +L**2*(MHuS+MHdS+MSS+Alcos1**2+IAl**2)
     .           -g1q*PAR(20)**2-3*g2q*PAR(21)**2
     .           +g1q/2d0*(2d0*MSQ1+MSQ3-4d0*MSU1-2d0*MSU3+2d0*MSD1+MSD3
     .                    -2d0*MSL1-MSL3+2d0*MSE1+MSE3+MHuS-MHdS))
     .           /16d0/Pi**2*DLOG(QSTSB/Q2)

      MHdQ2=MHdS-(3d0*Ybq**2*(MSQ3+MSD3+MHuS+ABP**2)
     .           +Ytau**2*(MSL3+MSE3+MHdS+ATAU**2)
     .           +L**2*(MHuS+MHdS+MSS+Alcos1**2+IAl**2)
     .           -g1q*PAR(20)**2-3*g2q*PAR(21)**2
     .           +g1q/2d0*(2d0*MSQ1+MSQ3-4d0*MSU1-2d0*MSU3+2d0*MSD1+MSD3
     .                    -2d0*MSL1-MSL3+2d0*MSE1+MSE3+MHuS-MHdS))
     .           /16d0/Pi**2*DLOG(QSTSB/Q2)

      MSQ2=MSS-(2d0*L**2*(MHuS+MHdS+MSS+Alcos1**2+IAl**2)
     .         +2d0*K**2*(3d0*MSS+Akcos2**2+IAk**2))
     .           /16d0/Pi**2*DLOG(QSTSB/Q2)

      IF(CFLAG(1).EQ.0)RETURN

c      Preliminaries

      Nphi=10
      MW2=g2q*(vuq**2+vdq**2)/2d0
      MZ2=(g1q+g2q)*(vuq**2+vdq**2)/2d0

c      Electroweak minimum

      VEW=MHuS*vuq**2+MHdS*vdq**2+MSS*(muq/l)**2-2d0*muq*vuq*vdq*Alcos1
     .    +2d0/3d0*k*Akcos2*(muq/l)**3+muq**2*(vuq**2+vdq**2)
     .    +(l*vuq*vdq)**2-2d0*k/l*muq**2*DDCOS(phi0)*vuq*vdq
     .    +k**2*(muq/l)**4+(G1Q+G2Q)/8d0*(vuq**2-vdq**2)**2
     .    -2d0*(M3HQ*DDCOS(phi3q)+l*XIFQ*DDCOS(phi01-phiF)
     .                        +MUPQ*muq*DDCOS(phi01-phiP))*vuq*vdq
     .    +2d0*(XISQ*DDCOS(phiSq)+XIFQ*MUPQ*DDCOS(phiP-phiF))*muq/l
     .    +2d0*(MSPQ/2d0*DDCOS(phiSPq)+k*XIFQ*DDCOS(phi02-phiF))
     .     *(muq/l)**2+2d0*k*MUPQ*DDCOS(phi02-phiP)*(muq/l)**3
     .    +(MUPQ*muq/l)**2
     .    -(3d0*mtopq**4*(dlog(mtopq**2/QSTSB)-1.5d0)
     .     +3d0*mbotq**4*(dlog(mbotq**2/QSTSB)-1.5d0)
     .     +(mtau*vdq/vd)**4*(dlog(mtau**2/QSTSB)-1.5d0))/16d0/Pi**2
     .    +3d0*(2d0*MW2**2*(dlog(MW2/QSTSB)-1.5d0)
     .                +MZ2*(dlog(MZ2/QSTSB)-1.5d0))/64d0/Pi**2
     .    +(3d0*MST2(1)**2*(dlog(MST2(1)/QSTSB)-1.5d0)
     .     +3d0*MST2(2)**2*(dlog(MST2(2)/QSTSB)-1.5d0)
     .     +3d0*MSB2(1)**2*(dlog(MSB2(1)/QSTSB)-1.5d0)
     .     +3d0*MSB2(2)**2*(dlog(MSB2(2)/QSTSB)-1.5d0)
     .     +MSL2(1)**2*(dlog(MSL2(1)/QSTSB)-1.5d0)
     .     +MSL2(2)**2*(dlog(MSL2(2)/QSTSB)-1.5d0)
     .     +MSNT2**2*(dlog(MSNT2/QSTSB)-1.5d0)
     .     +6d0*MSU2(1)**2*(dlog(MSU2(1)/QSTSB)-1.5d0)
     .     +6d0*MSU2(2)**2*(dlog(MSU2(2)/QSTSB)-1.5d0)
     .     +6d0*MSD2(1)**2*(dlog(MSD2(1)/QSTSB)-1.5d0)
     .     +6d0*MSD2(2)**2*(dlog(MSD2(2)/QSTSB)-1.5d0)
     .     +2d0*MSE2(1)**2*(dlog(MSE2(1)/QSTSB)-1.5d0)
     .     +2d0*MSE2(2)**2*(dlog(MSE2(2)/QSTSB)-1.5d0)
     .     +2d0*MSNE2**2*(dlog(MSNE2/QSTSB)-1.5d0))/32d0/Pi**2
     .    +3d0/512d0/Pi**4*Ytq**4*((dlog(QSTSB/mtopq**2))**2
     .      *(64d0*Pi*ALSQ+4d0/3d0*g1q-3d0*sinb**2*Ytq**2
     .       +3d0*cosb**2*Ybq**2)+
     .       ((dlog(MA2/mtopq**2))**2-(dlog(QSTSB/mtopq**2))**2)
     .       *(3d0*cosb**2*Ytq**2+(3d0*cosb**2+1d0)*Ybq**2))*vuq**4
     .    +3d0/512d0/Pi**4*Ybq**4*(dlog(QSTSB/mtopq**2)**2
     .      *(64d0*Pi*ALSQ-2d0/3d0*g1q+3d0*sinb**2*Ytq**2
     .       -3d0*cosb**2*Ybq**2)+
     .       (dlog(MA2/mtopq**2)**2-dlog(QSTSB/mtopq**2)**2)
     .       *(3d0*sinb**2*Ybq**2+(3d0*sinb**2+1d0)*Ytq**2))*vdq**4
     .    -(MCH2(1)**2*(dlog(MCH2(1)/QSTSB)-1.5d0)
     .     +MCH2(2)**2*(dlog(MCH2(2)/QSTSB)-1.5d0))/16d0/Pi**2
     .    -(MNEU(1)**2*(dlog(MNEU(1)/QSTSB)-1.5d0)
     .     +MNEU(2)**2*(dlog(MNEU(2)/QSTSB)-1.5d0)
     .     +MNEU(3)**2*(dlog(MNEU(3)/QSTSB)-1.5d0)
     .     +MNEU(4)**2*(dlog(MNEU(4)/QSTSB)-1.5d0)
     .     +MNEU(5)**2*(dlog(MNEU(5)/QSTSB)-1.5d0))/32d0/Pi**2
     .    +MHC**2*(dlog(MHC/QSTSB)-1.5d0)/32d0/Pi**2
     .    +(MH0(1)**2*(dlog(MH0(1)/QSTSB)-1.5d0)
     .     +MH0(2)**2*(dlog(MH0(2)/QSTSB)-1.5d0)
     .     +MH0(3)**2*(dlog(MH0(3)/QSTSB)-1.5d0)
     .     +MH0(4)**2*(dlog(MH0(4)/QSTSB)-1.5d0)
     .     +MH0(5)**2*(dlog(MH0(5)/QSTSB)-1.5d0))/64d0/Pi**2


c      Z3-symmetric minimum (or with vanishing vev's for the Z3-violating case)

      Taux=MHuS+MHdS
      Daux=MHuS-MHdS
      aux=M3HQ**2+(l*XIFQ)**2+2d0*M3HQ*l*XIFQ*DDCOS(phi3q+phiF-phi01)
      MHD1aux=(Taux-dsqrt(Daux**2+4d0*aux))/2d0
      MHD2aux=(Taux+dsqrt(Daux**2+4d0*aux))/2d0

      Taux=2d0*(MSS+MUPQ**2)
      Daux=MSPQ**2+(2d0*k*XIFQ)**2+4d0*k*XIFQ*MSPQ*DDCOS(phi02-phiF)
      MS1aux=(Taux-dsqrt(Daux))/2d0
      MS2aux=(Taux+dsqrt(Daux))/2d0

      VZ3=(3d0*MSQ3**2*(dlog(MSQ3/QSTSB)-1.5d0)
     .     +3d0*MSU3**2*(dlog(MSU3/QSTSB)-1.5d0)
     .     +3d0*MSQ3**2*(dlog(MSQ3/QSTSB)-1.5d0)
     .     +3d0*MSD3**2*(dlog(MSD3/QSTSB)-1.5d0)
     .     +MSL3**2*(dlog(MSL3/QSTSB)-1.5d0)
     .     +MSE3**2*(dlog(MSE3/QSTSB)-1.5d0)
     .     +MSL3**2*(dlog(MSL3/QSTSB)-1.5d0)
     .     +6d0*MSQ1**2*(dlog(MSQ1/QSTSB)-1.5d0)
     .     +6d0*MSU1**2*(dlog(MSU1/QSTSB)-1.5d0)
     .     +6d0*MSQ1**2*(dlog(MSQ1/QSTSB)-1.5d0)
     .     +6d0*MSD1**2*(dlog(MSD1/QSTSB)-1.5d0)
     .     +2d0*MSL1**2*(dlog(MSL1/QSTSB)-1.5d0)
     .     +2d0*MSE1**2*(dlog(MSE1/QSTSB)-1.5d0)
     .     +2d0*MSL1**2*(dlog(MSL1/QSTSB)-1.5d0))/32d0/Pi**2
     .    -(M2r**4*(dlog(M2r**2/QSTSB)-1.5d0))/16d0/Pi**2
     .    -(M1r**4*(dlog(M1r**2/QSTSB)-1.5d0)
     .     +M2r**4*(dlog(M2r**2/QSTSB)-1.5d0)
     .     +mupsi**4*(dlog(Max(mupsi**2,1d0)/QSTSB)-1.5d0))/32d0/Pi**2
     .     +(4d0*MHD1aux**2*(dlog(Max(MHD1aux,MZ2)/QSTSB)-1.5d0)
     .      +4d0*MHD2aux**2*(dlog(Max(MHD2aux,MZ2)/QSTSB)-1.5d0)
     .      +2d0*MS1aux**2*(dlog(Max(MS1aux,MZ2)/QSTSB)-1.5d0)
     .      +2d0*MS2aux**2*(dlog(Max(MS2aux,MZ2)/QSTSB)-1.5d0))
     .                                             /64d0/Pi**2

      IF(VZ3.lt.VEW)PROB(28)=DDIM(-1d-2,(VZ3-VEW)/Max(dabs(VEW),1d-10))


c      EW-symmetric minimum

      Taux=MHuS+MHdS+2d0*muq**2
      Daux=MHuS-MHdS
      aux=M3HQ**2+(l*XIFQ)**2+2d0*M3HQ*l*XIFQ*DDCOS(phi3q+phiF-phi01)
     . +(MUPQ*muq)**2+2d0*M3HQ*MUPQ*muq*DDCOS(PhiP-Phi01-Phi3q)
     . +2d0*l*XIFQ*MUPQ*muq*DDCOS(PhiP-PhiF)
     . +muq**2*(Alcos1**2+IAL**2)+(k*muq**2/l)**2
     . +2d0*muq**3*k/l*(Alcos1*DDCOS(phi0)+IAL*DDSIN(phi0))
     . +2d0*muq*Alcos1*(M3HQ*DDCOS(phi3q)+l*XIFQ*DDCOS(phi01-PhiF)
     .                  +muq*MUPQ*DDCOS(phi01-phiP))
     . +2d0*muq*IAL*(M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-PhiF)
     .                  +muq*MUPQ*DDSIN(phi01-phiP))
     . +2d0*k/l*muq**2*(M3HQ*DDCOS(phi3q-phi0)+l*XIFQ*DDCOS(phi02-PhiF)
     .                  +muq*MUPQ*DDCOS(phi02-phiP))
      MHD1aux=(Taux-dsqrt(Daux**2+4d0*aux))/2d0
      MHD2aux=(Taux+dsqrt(Daux**2+4d0*aux))/2d0

      Taux=2d0*(MSS+MUPQ**2)
     .     +4d0*((k*muq/l)**2+k/l*MUPQ*muq*DDCOS(phi02-phiP))
      Daux=2d0*MSPQ*DDCOS(phiSPq)+4d0*k*XIFQ*DDCOS(phi02-phiF)
     . +8d0*(k/l*muq)**2+4d0*k/l*muq*(Akcos2+2d0*MUPQ*DDCOS(phi02-phiP))
      aux=MSPQ*DDSIN(PhiSPq)+2d0*k*XIFQ*DDSIN(phi02-phiF)
     . +2d0*k/l*muq*(IAk+MUPQ*DDSIN(phi02-phiP))
      MS1aux=(Taux-dsqrt(Daux**2+4d0*aux**2))/2d0
      MS2aux=(Taux+dsqrt(Daux**2+4d0*aux**2))/2d0

      Vgauge=MSS*(muq/l)**2+2d0/3d0*k*Akcos2*(muq/l)**3+k**2*(muq/l)**4
     .    +2d0*(XISQ*DDCOS(phiSq)+XIFQ*MUPQ*DDCOS(phiP-phiF))*muq/l
     .    +2d0*(MSPQ/2d0*DDCOS(phiSPq)+k*XIFQ*DDCOS(phi02-phiF))
     .     *(muq/l)**2+2d0*k*MUPQ*DDCOS(phi02-phiP)*(muq/l)**3
     .    +(MUPQ*muq/l)**2
     .    +(3d0*MSQ3**2*(dlog(MSQ3/QSTSB)-1.5d0)
     .     +3d0*MSU3**2*(dlog(MSU3/QSTSB)-1.5d0)
     .     +3d0*MSQ3**2*(dlog(MSQ3/QSTSB)-1.5d0)
     .     +3d0*MSD3**2*(dlog(MSD3/QSTSB)-1.5d0)
     .     +MSL3**2*(dlog(MSL3/QSTSB)-1.5d0)
     .     +MSE3**2*(dlog(MSE3/QSTSB)-1.5d0)
     .     +MSL3**2*(dlog(MSL3/QSTSB)-1.5d0)
     .     +6d0*MSQ1**2*(dlog(MSQ1/QSTSB)-1.5d0)
     .     +6d0*MSU1**2*(dlog(MSU1/QSTSB)-1.5d0)
     .     +6d0*MSQ1**2*(dlog(MSQ1/QSTSB)-1.5d0)
     .     +6d0*MSD1**2*(dlog(MSD1/QSTSB)-1.5d0)
     .     +2d0*MSL1**2*(dlog(MSL1/QSTSB)-1.5d0)
     .     +2d0*MSE1**2*(dlog(MSE1/QSTSB)-1.5d0)
     .     +2d0*MSL1**2*(dlog(MSL1/QSTSB)-1.5d0))/32d0/Pi**2
     .    -(M2r**4*(dlog(M2r**2/QSTSB)-1.5d0)
     .     +mur**4*(dlog(mur**2/QSTSB)-1.5d0))/16d0/Pi**2
     .    -(M1r**4*(dlog(M1r**2/QSTSB)-1.5d0)
     .     +M2r**4*(dlog(M2r**2/QSTSB)-1.5d0)
     .     +2d0*mur**4*(dlog(mur**2/QSTSB)-1.5d0)
     .     +mupsi**4*(dlog(Max(mupsi**2,1d0)/QSTSB)-1.5d0))/32d0/Pi**2
     .     +(4d0*MHD1aux**2*(dlog(Max(MHD1aux,MZ2)/QSTSB)-1.5d0)
     .      +4d0*MHD2aux**2*(dlog(Max(MHD2aux,MZ2)/QSTSB)-1.5d0)
     .      +2d0*MS1aux**2*(dlog(Max(MS1aux,MZ2)/QSTSB)-1.5d0)
     .      +2d0*MS2aux**2*(dlog(Max(MS2aux,MZ2)/QSTSB)-1.5d0))
     .                                             /64d0/Pi**2

      Aroot=4d0*k**2
      Broot=2d0*k*Akcos2+6d0*k*MUPQ*DDCOS(phi02-phiP)
      Croot=2d0*MSS+4d0*(MSPQ/2d0*DDCOS(phiSPq)
     .     +k*XIFQ*DDCOS(phi02-phiF))+2d0*MUPQ**2
      Droot=2d0*(XISQ*DDCOS(phiSq)+XIFQ*MUPQ*DDCOS(phiP-phiF))
      root1=0d0
      root2=0d0
      root3=0d0

        IF(Aroot.ne.0d0)THEN
      Proot=-(Broot/Aroot)**2/3d0+Croot/Aroot
      Qroot=Broot/27d0/Aroot*(2d0*(Broot/Aroot)**2-9d0*Croot/Aroot)
     .     +Droot/Aroot
      Delt=-4d0*Proot**3-27d0*Qroot**2
         IF(Delt.gt.0d0)THEN
      root1=-Broot/3d0/Aroot+2d0*dsqrt(-Proot/3d0)
     .      *DDCOS(DACOS(-Qroot/2d0*dsqrt(27d0/(-Proot**3))))
      root2=-Broot/3d0/Aroot+2d0*dsqrt(-Proot/3d0)
     .      *DDCOS(DACOS(-Qroot/2d0*dsqrt(27d0/(-Proot**3)))+2d0/3d0*Pi)
      root3=-Broot/3d0/Aroot+2d0*dsqrt(-Proot/3d0)
     .      *DDCOS(DACOS(-Qroot/2d0*dsqrt(27d0/(-Proot**3)))+4d0/3d0*Pi)
         ELSEIF(Delt.eq.0d0)THEN
      root1=-Broot/3d0/Aroot
         ELSEIF(Delt.lt.0d0)THEN
      root1=-Broot/3d0/Aroot
     .     +((-Qroot+dsqrt(-Delt/27d0))/2d0)**(1d0/3d0)
     .     +((-Qroot-dsqrt(-Delt/27d0))/2d0)**(1d0/3d0)
         ENDIF
        ELSE
         IF(Broot.ne.0d0)THEN
      Proot=Croot/Broot
      Qroot=Droot/Broot
      Delt=Proot**2-4d0*Qroot
          IF(Delt.gt.0d0)THEN
      root1=(-Proot-dsqrt(Delt))/2d0
      root2=(-Proot+dsqrt(Delt))/2d0
          ELSEIF(Delt.eq.0d0)THEN
      root1=-Proot/2d0
          ENDIF
         ELSE
          IF(Croot.ne.0d0)root1=-Droot/Croot
         ENDIF
        ENDIF

        If(root1.ne.0d0)then
      Taux=MHuS+MHdS+2d0*(l*root1)**2
      Daux=MHuS-MHdS
      aux=M3HQ**2+(l*XIFQ)**2+2d0*M3HQ*l*XIFQ*DDCOS(phi3q+phiF-phi01)
     . +(MUPQ*l*root1)**2+2d0*M3HQ*MUPQ*l*root1*DDCOS(PhiP-Phi01-Phi3q)
     . +2d0*l*XIFQ*MUPQ*l*root1*DDCOS(PhiP-PhiF)
     . +(l*root1)**2*(Alcos1**2+IAL**2)+(k*l*root1**2)**2
     . +2d0*k*l**2*root1**3*(Alcos1*DDCOS(phi0)+IAL*DDSIN(phi0))
     . +2d0*l*root1*Alcos1*(M3HQ*DDCOS(phi3q)+l*XIFQ*DDCOS(phi01-PhiF)
     .                  +l*root1*MUPQ*DDCOS(phi01-phiP))
     . +2d0*l*root1*IAL*(M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-PhiF)
     .                  +l*root1*MUPQ*DDSIN(phi01-phiP))
     . +2d0*k*l*root1**2*(M3HQ*DDCOS(phi3q-phi0)
     .                  +l*XIFQ*DDCOS(phi02-PhiF)
     .                  +l*root1*MUPQ*DDCOS(phi02-phiP))
      MHD1aux=(Taux-dsqrt(Daux**2+4d0*aux))/2d0
      MHD2aux=(Taux+dsqrt(Daux**2+4d0*aux))/2d0

      Taux=2d0*(MSS+MUPQ**2)
     .     +4d0*((k*root1)**2+k*MUPQ*root1*DDCOS(phi02-phiP))
      Daux=2d0*MSPQ*DDCOS(phiSPq)+4d0*k*XIFQ*DDCOS(phi02-phiF)
     . +8d0*(k*root1)**2+4d0*k*root1*(Akcos2+2d0*MUPQ*DDCOS(phi02-phiP))
      aux=MSPQ*DDSIN(PhiSPq)+2d0*k*XIFQ*DDSIN(phi02-phiF)
     . +2d0*k*root1*(IAk+MUPQ*DDSIN(phi02-phiP))
      MS1aux=(Taux-dsqrt(Daux**2+4d0*aux**2))/2d0
      MS2aux=(Taux+dsqrt(Daux**2+4d0*aux**2))/2d0

      Vgauge=min(Vgauge,
     .     MSS*(root1)**2+2d0/3d0*k*Akcos2*(root1)**3+k**2*(root1)**4
     .    +2d0*(XISQ*DDCOS(phiSq)+XIFQ*MUPQ*DDCOS(phiP-phiF))*root1
     .    +2d0*(MSPQ/2d0*DDCOS(phiSPq)+k*XIFQ*DDCOS(phi02-phiF))
     .     *(root1)**2+2d0*k*MUPQ*DDCOS(phi02-phiP)*(root1)**3
     .    +(MUPQ*root1)**2
     .    +(3d0*MSQ3**2*(dlog(MSQ3/QSTSB)-1.5d0)
     .     +3d0*MSU3**2*(dlog(MSU3/QSTSB)-1.5d0)
     .     +3d0*MSQ3**2*(dlog(MSQ3/QSTSB)-1.5d0)
     .     +3d0*MSD3**2*(dlog(MSD3/QSTSB)-1.5d0)
     .     +MSL3**2*(dlog(MSL3/QSTSB)-1.5d0)
     .     +MSE3**2*(dlog(MSE3/QSTSB)-1.5d0)
     .     +MSL3**2*(dlog(MSL3/QSTSB)-1.5d0)
     .     +6d0*MSQ1**2*(dlog(MSQ1/QSTSB)-1.5d0)
     .     +6d0*MSU1**2*(dlog(MSU1/QSTSB)-1.5d0)
     .     +6d0*MSQ1**2*(dlog(MSQ1/QSTSB)-1.5d0)
     .     +6d0*MSD1**2*(dlog(MSD1/QSTSB)-1.5d0)
     .     +2d0*MSL1**2*(dlog(MSL1/QSTSB)-1.5d0)
     .     +2d0*MSE1**2*(dlog(MSE1/QSTSB)-1.5d0)
     .     +2d0*MSL1**2*(dlog(MSL1/QSTSB)-1.5d0))/32d0/Pi**2
     .    -(M2r**4*(dlog(M2r**2/QSTSB)-1.5d0)
     .     +(l*root1)**4*(dlog((l*root1)**2/QSTSB)-1.5d0))/16d0/Pi**2
     .    -(M1r**4*(dlog(M1r**2/QSTSB)-1.5d0)
     .     +M2r**4*(dlog(M2r**2/QSTSB)-1.5d0)
     .     +2d0*(l*root1)**4*(dlog((l*root1)**2/QSTSB)-1.5d0)
     .     +(mupsi+2d0*k/l*root1)**4
     .         *(dlog(Max((mupsi+2d0*k/l*root1)**2,1d0)/QSTSB)-1.5d0)
     .                                                   )/32d0/Pi**2
     .     +(4d0*MHD1aux**2*(dlog(Max(MHD1aux,MZ2)/QSTSB)-1.5d0)
     .      +4d0*MHD2aux**2*(dlog(Max(MHD2aux,MZ2)/QSTSB)-1.5d0)
     .      +2d0*MS1aux**2*(dlog(Max(MS1aux,MZ2)/QSTSB)-1.5d0)
     .      +2d0*MS2aux**2*(dlog(Max(MS2aux,MZ2)/QSTSB)-1.5d0))
     .                                             /64d0/Pi**2)
        Endif
        If(root2.ne.0d0)then
      Taux=MHuS+MHdS+2d0*(l*root2)**2
      Daux=MHuS-MHdS
      aux=M3HQ**2+(l*XIFQ)**2+2d0*M3HQ*l*XIFQ*DDCOS(phi3q+phiF-phi01)
     . +(MUPQ*l*root2)**2+2d0*M3HQ*MUPQ*l*root2*DDCOS(PhiP-Phi01-Phi3q)
     . +2d0*l*XIFQ*MUPQ*l*root2*DDCOS(PhiP-PhiF)
     . +(l*root2)**2*(Alcos1**2+IAL**2)+(k*l*root2**2)**2
     . +2d0*k*l**2*root2**3*(Alcos1*DDCOS(phi0)+IAL*DDSIN(phi0))
     . +2d0*l*root2*Alcos1*(M3HQ*DDCOS(phi3q)+l*XIFQ*DDCOS(phi01-PhiF)
     .                  +l*root2*MUPQ*DDCOS(phi01-phiP))
     . +2d0*l*root2*IAL*(M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-PhiF)
     .                  +l*root2*MUPQ*DDSIN(phi01-phiP))
     . +2d0*k*l*root2**2*(M3HQ*DDCOS(phi3q-phi0)
     .                  +l*XIFQ*DDCOS(phi02-PhiF)
     .                  +l*root2*MUPQ*DDCOS(phi02-phiP))
      MHD1aux=(Taux-dsqrt(Daux**2+4d0*aux))/2d0
      MHD2aux=(Taux+dsqrt(Daux**2+4d0*aux))/2d0

      Taux=2d0*(MSS+MUPQ**2)
     .     +4d0*((k*root2)**2+k*MUPQ*root2*DDCOS(phi02-phiP))
      Daux=2d0*MSPQ*DDCOS(phiSPq)+4d0*k*XIFQ*DDCOS(phi02-phiF)
     . +8d0*(k*root2)**2+4d0*k*root2*(Akcos2+2d0*MUPQ*DDCOS(phi02-phiP))
      aux=MSPQ*DDSIN(PhiSPq)+2d0*k*XIFQ*DDSIN(phi02-phiF)
     . +2d0*k*root2*(IAk+MUPQ*DDSIN(phi02-phiP))
      MS1aux=(Taux-dsqrt(Daux**2+4d0*aux**2))/2d0
      MS2aux=(Taux+dsqrt(Daux**2+4d0*aux**2))/2d0

      Vgauge=min(Vgauge,
     .     MSS*(root2)**2+2d0/3d0*k*Akcos2*(root2)**3+k**2*(root2)**4
     .    +2d0*(XISQ*DDCOS(phiSq)+XIFQ*MUPQ*DDCOS(phiP-phiF))*root2
     .    +2d0*(MSPQ/2d0*DDCOS(phiSPq)+k*XIFQ*DDCOS(phi02-phiF))
     .     *(root2)**2+2d0*k*MUPQ*DDCOS(phi02-phiP)*(root2)**3
     .    +(MUPQ*root2)**2
     .    +(3d0*MSQ3**2*(dlog(MSQ3/QSTSB)-1.5d0)
     .     +3d0*MSU3**2*(dlog(MSU3/QSTSB)-1.5d0)
     .     +3d0*MSQ3**2*(dlog(MSQ3/QSTSB)-1.5d0)
     .     +3d0*MSD3**2*(dlog(MSD3/QSTSB)-1.5d0)
     .     +MSL3**2*(dlog(MSL3/QSTSB)-1.5d0)
     .     +MSE3**2*(dlog(MSE3/QSTSB)-1.5d0)
     .     +MSL3**2*(dlog(MSL3/QSTSB)-1.5d0)
     .     +6d0*MSQ1**2*(dlog(MSQ1/QSTSB)-1.5d0)
     .     +6d0*MSU1**2*(dlog(MSU1/QSTSB)-1.5d0)
     .     +6d0*MSQ1**2*(dlog(MSQ1/QSTSB)-1.5d0)
     .     +6d0*MSD1**2*(dlog(MSD1/QSTSB)-1.5d0)
     .     +2d0*MSL1**2*(dlog(MSL1/QSTSB)-1.5d0)
     .     +2d0*MSE1**2*(dlog(MSE1/QSTSB)-1.5d0)
     .     +2d0*MSL1**2*(dlog(MSL1/QSTSB)-1.5d0))/32d0/Pi**2
     .    -(M2r**4*(dlog(M2r**2/QSTSB)-1.5d0)
     .     +(l*root2)**4*(dlog((l*root2)**2/QSTSB)-1.5d0))/16d0/Pi**2
     .    -(M1r**4*(dlog(M1r**2/QSTSB)-1.5d0)
     .     +M2r**4*(dlog(M2r**2/QSTSB)-1.5d0)
     .     +2d0*(l*root2)**4*(dlog((l*root2)**2/QSTSB)-1.5d0)
     .     +(mupsi+2d0*k/l*root2)**4
     .         *(dlog(Max((mupsi+2d0*k/l*root2)**2,1d0)/QSTSB)-1.5d0)
     .                                                   )/32d0/Pi**2
     .     +(4d0*MHD1aux**2*(dlog(Max(MHD1aux,MZ2)/QSTSB)-1.5d0)
     .      +4d0*MHD2aux**2*(dlog(Max(MHD2aux,MZ2)/QSTSB)-1.5d0)
     .      +2d0*MS1aux**2*(dlog(Max(MS1aux,MZ2)/QSTSB)-1.5d0)
     .      +2d0*MS2aux**2*(dlog(Max(MS2aux,MZ2)/QSTSB)-1.5d0))
     .                                             /64d0/Pi**2)
        Endif
        If(root3.ne.0d0)then
      Taux=MHuS+MHdS+2d0*(l*root3)**2
      Daux=MHuS-MHdS
      aux=M3HQ**2+(l*XIFQ)**2+2d0*M3HQ*l*XIFQ*DDCOS(phi3q+phiF-phi01)
     . +(MUPQ*l*root3)**2+2d0*M3HQ*MUPQ*l*root3*DDCOS(PhiP-Phi01-Phi3q)
     . +2d0*l*XIFQ*MUPQ*l*root3*DDCOS(PhiP-PhiF)
     . +(l*root3)**2*(Alcos1**2+IAL**2)+(k*l*root3**2)**2
     . +2d0*k*l**2*root3**3*(Alcos1*DDCOS(phi0)+IAL*DDSIN(phi0))
     . +2d0*l*root3*Alcos1*(M3HQ*DDCOS(phi3q)+l*XIFQ*DDCOS(phi01-PhiF)
     .                  +l*root3*MUPQ*DDCOS(phi01-phiP))
     . +2d0*l*root3*IAL*(M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-PhiF)
     .                  +l*root3*MUPQ*DDSIN(phi01-phiP))
     . +2d0*k*l*root3**2*(M3HQ*DDCOS(phi3q-phi0)
     .                  +l*XIFQ*DDCOS(phi02-PhiF)
     .                  +l*root3*MUPQ*DDCOS(phi02-phiP))
      MHD1aux=(Taux-dsqrt(Daux**2+4d0*aux))/2d0
      MHD2aux=(Taux+dsqrt(Daux**2+4d0*aux))/2d0

      Taux=2d0*(MSS+MUPQ**2)
     .     +4d0*((k*root3)**2+k*MUPQ*root3*DDCOS(phi02-phiP))
      Daux=2d0*MSPQ*DDCOS(phiSPq)+4d0*k*XIFQ*DDCOS(phi02-phiF)
     . +8d0*(k*root3)**2+4d0*k*root3*(Akcos2+2d0*MUPQ*DDCOS(phi02-phiP))
      aux=MSPQ*DDSIN(PhiSPq)+2d0*k*XIFQ*DDSIN(phi02-phiF)
     . +2d0*k*root3*(IAk+MUPQ*DDSIN(phi02-phiP))
      MS1aux=(Taux-dsqrt(Daux**2+4d0*aux**2))/2d0
      MS2aux=(Taux+dsqrt(Daux**2+4d0*aux**2))/2d0
      Vgauge=min(Vgauge,
     .     MSS*(root3)**2+2d0/3d0*k*Akcos2*(root3)**3+k**2*(root3)**4
     .    +2d0*(XISQ*DDCOS(phiSq)+XIFQ*MUPQ*DDCOS(phiP-phiF))*root3
     .    +2d0*(MSPQ/2d0*DDCOS(phiSPq)+k*XIFQ*DDCOS(phi02-phiF))
     .     *(root3)**2+2d0*k*MUPQ*DDCOS(phi02-phiP)*(root3)**3
     .    +(MUPQ*root3)**2
     .    +(3d0*MSQ3**2*(dlog(MSQ3/QSTSB)-1.5d0)
     .     +3d0*MSU3**2*(dlog(MSU3/QSTSB)-1.5d0)
     .     +3d0*MSQ3**2*(dlog(MSQ3/QSTSB)-1.5d0)
     .     +3d0*MSD3**2*(dlog(MSD3/QSTSB)-1.5d0)
     .     +MSL3**2*(dlog(MSL3/QSTSB)-1.5d0)
     .     +MSE3**2*(dlog(MSE3/QSTSB)-1.5d0)
     .     +MSL3**2*(dlog(MSL3/QSTSB)-1.5d0)
     .     +6d0*MSQ1**2*(dlog(MSQ1/QSTSB)-1.5d0)
     .     +6d0*MSU1**2*(dlog(MSU1/QSTSB)-1.5d0)
     .     +6d0*MSQ1**2*(dlog(MSQ1/QSTSB)-1.5d0)
     .     +6d0*MSD1**2*(dlog(MSD1/QSTSB)-1.5d0)
     .     +2d0*MSL1**2*(dlog(MSL1/QSTSB)-1.5d0)
     .     +2d0*MSE1**2*(dlog(MSE1/QSTSB)-1.5d0)
     .     +2d0*MSL1**2*(dlog(MSL1/QSTSB)-1.5d0))/32d0/Pi**2
     .    -(M2r**4*(dlog(M2r**2/QSTSB)-1.5d0)
     .     +(l*root3)**4*(dlog((l*root3)**2/QSTSB)-1.5d0))/16d0/Pi**2
     .    -(M1r**4*(dlog(M1r**2/QSTSB)-1.5d0)
     .     +M2r**4*(dlog(M2r**2/QSTSB)-1.5d0)
     .     +2d0*(l*root3)**4*(dlog((l*root3)**2/QSTSB)-1.5d0)
     .     +(mupsi+2d0*k/l*root3)**4
     .         *(dlog(Max((mupsi+2d0*k/l*root3)**2,1d0)/QSTSB)-1.5d0)
     .                                                   )/32d0/Pi**2
     .     +(4d0*MHD1aux**2*(dlog(Max(MHD1aux,MZ2)/QSTSB)-1.5d0)
     .      +4d0*MHD2aux**2*(dlog(Max(MHD2aux,MZ2)/QSTSB)-1.5d0)
     .      +2d0*MS1aux**2*(dlog(Max(MS1aux,MZ2)/QSTSB)-1.5d0)
     .      +2d0*MS2aux**2*(dlog(Max(MS2aux,MZ2)/QSTSB)-1.5d0))
     .                                             /64d0/Pi**2)
        Endif

      IF(Vgauge.lt.VEW)PROB(28)=Max(PROB(28),
     .                   DDIM(-1d-2,(Vgauge-VEW)/Max(dabs(VEW),1d-10)))


c      Minimum with <S>=<Hd>=0

      Taux=MHuS+MHdS+g2q/2d0*vuq**2
      Daux=MHuS-MHdS+g1q/2d0*vuq**2
      aux=M3HQ**2+(l*XIFQ)**2+2d0*M3HQ*l*XIFQ*DDCOS(phi3q+phiF-phi01)
      MHD1aux=(Taux-dsqrt(Daux**2+4d0*aux))/2d0
      MHD2aux=(Taux+dsqrt(Daux**2+4d0*aux))/2d0

      MH2(1,1)=MHuS+3d0/4d0*(g1q+g2q)*vuq**2
      MH2(1,2)=-M3HQ*DDCOS(phi3q)-l*XIFQ*DDCOS(phi01-phiF)
      MH2(2,1)=MH2(1,2)
      MH2(1,3)=0d0
      MH2(3,1)=MH2(1,3)
      MH2(2,2)=MHdS+(l**2-(g1q+g2q)/4d0)*vuq**2
      MH2(2,3)=-l*(Alcos1+MUPQ*DDCOS(phi01-phiP))*vuq
      MH2(3,2)=MH2(2,3)
      MH2(3,3)=MSS+MUPQ**2+(l*vuq)**2+MSPQ*DDCOS(phiSPq)
     .                               +2d0*k*XIFQ*DDCOS(phi02-phiF)
      MH2(1,4)=0d0
      MH2(4,1)=MH2(1,4)
      MH2(2,4)=M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-phiF)
      MH2(4,2)=MH2(2,4)
      MH2(3,4)=0d0
      MH2(4,3)=MH2(3,4)
      MH2(4,4)=MHuS+(g1q+g2q)/4d0*vuq**2
      MH2(1,5)=M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-phiF)
      MH2(5,1)=MH2(1,5)
      MH2(2,5)=0d0
      MH2(5,2)=MH2(2,5)
      MH2(3,5)=l*(IAl+MUPQ*DDSIN(phi01-phiP))*vuq
      MH2(5,3)=MH2(3,5)
      MH2(4,5)=M3HQ*DDCOS(phi3q)+l*XIFQ*DDCOS(phi01-phiF)
      MH2(5,4)=MH2(4,5)
      MH2(5,5)=MHdS+(l**2-(g1q+g2q)/4d0)*vuq**2
      MH2(1,6)=0d0
      MH2(6,1)=MH2(1,6)
      MH2(2,6)=l*(IAl+MUPQ*DDSIN(phi01-phiP))*vuq
      MH2(6,2)=MH2(2,6)
      MH2(3,6)=-MSPQ*DDSIN(phiSPq)-2d0*k*XIFQ*DDSIN(phi02-phiF)
      MH2(6,3)=MH2(3,6)
      MH2(4,6)=0d0
      MH2(6,4)=MH2(4,6)
      MH2(5,6)=l*(Alcos1+MUPQ*DDCOS(phi01-phiP))*vuq
      MH2(6,5)=MH2(4,5)
      MH2(6,6)=MSS+MUPQ**2+l**2*vuq**2-MSPQ*DDCOS(phiSPq)
     .                               -2d0*k*XIFQ*DDCOS(phi02-phiF)

      CALL DIAGN(6,MH2,VALPH,VECPH,1d-10)

        MSTaux1=(MSQ3+MSU3+2d0*mtopq**2-(g1q+g2q)/4d0*vuq**2
     . -dsqrt((MSQ3-MSU3+(5d0/3d0*g1q-g2q)/4d0*vuq**2)**2
     .        +4d0*(Ytq*vuq*ATP)**2))/2d0
        MSTaux2=(MSQ3+MSU3+2d0*mtopq**2-(g1q+g2q)/4d0*vuq**2
     . +dsqrt((MSQ3-MSU3+(5d0/3d0*g1q-g2q)/4d0*vuq**2)**2
     .        +4d0*(Ytq*vuq*ATP)**2))/2d0
        MSBaux1=MSQ3+(g1q/3d0+g2q)/4d0*vuq**2
        MSBaux2=MSD3+g1q/6d0*vuq**2
        MSLaux1=MSL3+(-g1q+g2q)/4d0*vuq**2
        MSLaux2=MSE3+g1q/2d0*vuq**2
        MSNTaux=MSL3-(g1q+g2q)/4d0*vuq**2
        MSUaux1=MSQ1+(g1q/3d0-g2q)/4d0*vuq**2
        MSUaux2=MSU1-g1q/3d0*vuq**2
        MSDaux1=MSQ1+(g1q/3d0+g2q)/4d0*vuq**2
        MSDaux2=MSD1+g1q/6d0*vuq**2
        MSEaux1=MSL1+(-g1q+g2q)/4d0*vuq**2
        MSEaux2=MSE1+g1q/2d0*vuq**2
        MSNEaux=MSL1-(g1q+g2q)/4d0*vuq**2

      Vvu=MHuS*vuq**2+(G1Q+G2Q)/8d0*vuq**4
     .    -3d0*mtopq**4*(dlog(mtopq**2/QSTSB)-1.5d0)/16d0/Pi**2
     .    +3d0*(2d0*(g2q/2d0*vuq**2)**2
     .           *(dlog(g2q*vuq**2/2d0/QSTSB)-1.5d0)
     .         +((g1q+g2q)*vuq**2/2d0)**2
     .           *(dlog((g1q+g2q)*vuq**2/2d0/QSTSB)-1.5d0))/64d0/Pi**2
     .    +(3d0*MSTaux1**2*(dlog(Max(MSTaux1,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSTaux2**2*(dlog(Max(MSTaux2,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSBaux1**2*(dlog(Max(MSBaux1,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSBaux2**2*(dlog(Max(MSBaux2,MZ2)/QSTSB)-1.5d0)
     .     +MSLaux1**2*(dlog(Max(MSLaux1,MZ2)/QSTSB)-1.5d0)
     .     +MSLaux2**2*(dlog(Max(MSLaux2,MZ2)/QSTSB)-1.5d0)
     .     +MSNTaux**2*(dlog(Max(MSNTaux,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSUaux1**2*(dlog(Max(MSUaux1,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSUaux2**2*(dlog(Max(MSUaux2,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSDaux1**2*(dlog(Max(MSDaux1,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSDaux2**2*(dlog(Max(MSDaux2,MZ2)/QSTSB)-1.5d0)
     .     +2d0*MSEaux1**2*(dlog(Max(MSEaux1,MZ2)/QSTSB)-1.5d0)
     .     +2d0*MSEaux2**2*(dlog(Max(MSEaux2,MZ2)/QSTSB)-1.5d0)
     . +2d0*MSNEaux**2*(dlog(Max(MSNEaux,MZ2)/QSTSB)-1.5d0))/32d0/Pi**2
     .    +3d0/512d0/Pi**4*Ytq**4*((dlog(QSTSB/mtopq**2))**2
     .      *(64d0*Pi*ALSQ+4d0/3d0*g1q-3d0*sinb**2*Ytq**2
     .       +3d0*cosb**2*Ybq**2)+
     .       ((dlog(MA2/mtopq**2))**2-(dlog(QSTSB/mtopq**2))**2)
     .       *(3d0*cosb**2*Ytq**2+(3d0*cosb**2+1d0)*Ybq**2))*vuq**4
     .    -((M2r**2+g2q*vuq**2)**2
     .              *(dlog((M2r**2+g2q*vuq**2)/QSTSB)-1.5d0))/16d0/Pi**2
     .    -((M1r**2+g1q*vuq**2)**2
     .                   *(dlog((M1r**2+g1q*vuq**2)/QSTSB)-1.5d0)
     .     +(M2r**2+g2q*vuq**2)**2
     .                   *(dlog((M2r**2+g2q*vuq**2)/QSTSB)-1.5d0)
     .     +(mupsi**2+2d0*(l*vuq)**2)**2*
     . (dlog(Max(mupsi**2+2d0*(l*vuq)**2,1d0)/QSTSB)-1.5d0))/32d0/Pi**2
     .     +(MHD1aux**2*(dlog(Max(MHD1aux,MW2)/QSTSB)-1.5d0)
     .     +MHD2aux**2*(dlog(Max(MHD2aux,MW2)/QSTSB)-1.5d0))/32d0/Pi**2
     .    +(VALPH(1)**2*(dlog(Max(MZ2,VALPH(1))/QSTSB)-1.5d0)
     .     +VALPH(2)**2*(dlog(Max(MZ2,VALPH(2))/QSTSB)-1.5d0)
     .     +VALPH(3)**2*(dlog(Max(MZ2,VALPH(3))/QSTSB)-1.5d0)
     .     +VALPH(4)**2*(dlog(Max(MZ2,VALPH(4))/QSTSB)-1.5d0)
     .     +VALPH(5)**2*(dlog(Max(MZ2,VALPH(5))/QSTSB)-1.5d0)
     .   +VALPH(6)**2*(dlog(Max(MZ2,VALPH(6))/QSTSB)-1.5d0))/64d0/Pi**2

        If(MHuS.lt.0d0)then
      aux=dsqrt(-4d0*MHuS/(G1Q+G2Q))

      Taux=MHuS-MHdS+g1q/2d0*aux**2
      Daux=M3HQ**2+(l*XIFQ)**2+2d0*M3HQ*l*XIFQ*DDCOS(phi3q+phiF-phi01)
      MHD1aux=(MHuS+MHdS+g2q/2d0*aux**2-dsqrt(Taux**2+4d0*Daux))/2d0
      MHD2aux=(MHuS+MHdS+g2q/2d0*aux**2+dsqrt(Taux**2+4d0*Daux))/2d0

      MH2(1,1)=MHuS+3d0/4d0*(g1q+g2q)*aux**2
      MH2(1,2)=-M3HQ*DDCOS(phi3q)-l*XIFQ*DDCOS(phi01-phiF)
      MH2(2,1)=MH2(1,2)
      MH2(1,3)=0d0
      MH2(3,1)=MH2(1,3)
      MH2(2,2)=MHdS+(l**2-(g1q+g2q)/4d0)*aux**2
      MH2(2,3)=-l*(Alcos1+MUPQ*DDCOS(phi01-phiP))*aux
      MH2(3,2)=MH2(2,3)
      MH2(3,3)=MSS+MUPQ**2+(l*aux)**2+MSPQ*DDCOS(phiSPq)
     .                               +2d0*k*XIFQ*DDCOS(phi02-phiF)
      MH2(1,4)=0d0
      MH2(4,1)=MH2(1,4)
      MH2(2,4)=M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-phiF)
      MH2(4,2)=MH2(2,4)
      MH2(3,4)=0d0
      MH2(4,3)=MH2(3,4)
      MH2(4,4)=MHuS+(g1q+g2q)/4d0*aux**2
      MH2(1,5)=M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-phiF)
      MH2(5,1)=MH2(1,5)
      MH2(2,5)=0d0
      MH2(5,2)=MH2(2,5)
      MH2(3,5)=l*(IAl+MUPQ*DDSIN(phi01-phiP))*aux
      MH2(5,3)=MH2(3,5)
      MH2(4,5)=M3HQ*DDCOS(phi3q)+l*XIFQ*DDCOS(phi01-phiF)
      MH2(5,4)=MH2(4,5)
      MH2(5,5)=MHdS+(l**2-(g1q+g2q)/4d0)*aux**2
      MH2(1,6)=0d0
      MH2(6,1)=MH2(1,6)
      MH2(2,6)=l*(IAl+MUPQ*DDSIN(phi01-phiP))*aux
      MH2(6,2)=MH2(2,6)
      MH2(3,6)=-MSPQ*DDSIN(phiSPq)-2d0*k*XIFQ*DDSIN(phi02-phiF)
      MH2(6,3)=MH2(3,6)
      MH2(4,6)=0d0
      MH2(6,4)=MH2(4,6)
      MH2(5,6)=l*(Alcos1+MUPQ*DDCOS(phi01-phiP))*aux
      MH2(6,5)=MH2(4,5)
      MH2(6,6)=MSS+MUPQ**2+l**2*aux**2-MSPQ*DDCOS(phiSPq)
     .                               -2d0*k*XIFQ*DDCOS(phi02-phiF)

      CALL DIAGN(6,MH2,VALPH,VECPH,1d-10)

      MSTaux1=(MSQ3+MSU3+2d0*(Ytq*aux)**2-(g1q+g2q)/4d0*aux**2
     . -dsqrt((MSQ3-MSU3+(5d0/3d0*g1q-g2q)/4d0*aux**2)**2
     .        +4d0*(Ytq*aux*ATP)**2))/2d0
      MSTaux2=(MSQ3+MSU3+2d0*(Ytq*aux)**2-(g1q+g2q)/4d0*aux**2
     . +dsqrt((MSQ3-MSU3+(5d0/3d0*g1q-g2q)/4d0*aux**2)**2
     .        +4d0*(Ytq*aux*ATP)**2))/2d0
      MSBaux1=MSQ3+(g1q/3d0+g2q)/4d0*aux**2
      MSBaux2=MSD3+g1q/6d0*aux**2
      MSLaux1=MSL3+(-g1q+g2q)/4d0*aux**2
      MSLaux2=MSE3+g1q/2d0*aux**2
      MSNTaux=MSL3-(g1q+g2q)/4d0*aux**2
      MSNTaux=MSL3-(g1q+g2q)/4d0*aux**2
      MSUaux1=MSQ1+(g1q/3d0-g2q)/4d0*aux**2
      MSUaux2=MSU1-g1q/3d0*aux**2
      MSDaux1=MSQ1+(g1q/3d0+g2q)/4d0*aux**2
      MSDaux2=MSD1+g1q/6d0*aux**2
      MSEaux1=MSL1+(-g1q+g2q)/4d0*aux**2
      MSEaux2=MSE1+g1q/2d0*aux**2
      MSNEaux=MSL1-(g1q+g2q)/4d0*aux**2

      Vvu=min(Vvu,MHuS*aux**2+(G1Q+G2Q)/8d0*aux**4
     .    -3d0*(Ytq*aux)**4*(dlog((Ytq*aux)**2/QSTSB)-1.5d0)/16d0/Pi**2
     .    +3d0*(2d0*(g2q/2d0*aux**2)**2
     .           *(dlog(g2q*aux**2/2d0/QSTSB)-1.5d0)
     .         +((g1q+g2q)*aux**2/2d0)**2
     .           *(dlog((g1q+g2q)*aux**2/2d0/QSTSB)-1.5d0))/64d0/Pi**2
     .    +(3d0*MSTaux1**2*(dlog(Max(MSTaux1,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSTaux2**2*(dlog(Max(MSTaux2,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSBaux1**2*(dlog(Max(MSBaux1,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSBaux2**2*(dlog(Max(MSBaux2,MZ2)/QSTSB)-1.5d0)
     .     +MSLaux1**2*(dlog(Max(MSLaux1,MZ2)/QSTSB)-1.5d0)
     .     +MSLaux2**2*(dlog(Max(MSLaux2,MZ2)/QSTSB)-1.5d0)
     .     +MSNTaux**2*(dlog(Max(MSNTaux,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSUaux1**2*(dlog(Max(MSUaux1,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSUaux2**2*(dlog(Max(MSUaux2,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSDaux1**2*(dlog(Max(MSDaux1,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSDaux2**2*(dlog(Max(MSDaux2,MZ2)/QSTSB)-1.5d0)
     .     +2d0*MSEaux1**2*(dlog(Max(MSEaux1,MZ2)/QSTSB)-1.5d0)
     .     +2d0*MSEaux2**2*(dlog(Max(MSEaux2,MZ2)/QSTSB)-1.5d0)
     . +2d0*MSNEaux**2*(dlog(Max(MSNEaux,MZ2)/QSTSB)-1.5d0))/32d0/Pi**2
     .    +3d0/512d0/Pi**4*Ytq**4*((dlog(QSTSB/mtopq**2))**2
     .      *(64d0*Pi*ALSQ+4d0/3d0*g1q-3d0*sinb**2*Ytq**2
     .       +3d0*cosb**2*Ybq**2)+
     .       ((dlog(MA2/mtopq**2))**2-(dlog(QSTSB/mtopq**2))**2)
     .       *(3d0*cosb**2*Ytq**2+(3d0*cosb**2+1d0)*Ybq**2))*aux**4
     .    -((M2r**2+g2q*aux**2)**2
     .            *(dlog((M2r**2+g2q*aux**2)/QSTSB)-1.5d0))/16d0/Pi**2
     .    -((M1r**2+g1q*aux**2)**2
     .                   *(dlog((M1r**2+g1q*aux**2)/QSTSB)-1.5d0)
     .     +(M2r**2+g2q*aux**2)**2
     .                   *(dlog((M2r**2+g2q*aux**2)/QSTSB)-1.5d0)
     .     +(mupsi**2+2d0*(l*aux)**2)**2*
     . (dlog(Max(mupsi**2+2d0*(l*aux)**2,1d0)/QSTSB)-1.5d0))/32d0/Pi**2
     .     +(MHD1aux**2*(dlog(Max(MHD1aux,MW2)/QSTSB)-1.5d0)
     .     +MHD2aux**2*(dlog(Max(MHD2aux,MW2)/QSTSB)-1.5d0))/32d0/Pi**2
     .    +(VALPH(1)**2*(dlog(Max(MZ2,VALPH(1))/QSTSB)-1.5d0)
     .     +VALPH(2)**2*(dlog(Max(MZ2,VALPH(2))/QSTSB)-1.5d0)
     .     +VALPH(3)**2*(dlog(Max(MZ2,VALPH(3))/QSTSB)-1.5d0)
     .     +VALPH(4)**2*(dlog(Max(MZ2,VALPH(4))/QSTSB)-1.5d0)
     .     +VALPH(5)**2*(dlog(Max(MZ2,VALPH(5))/QSTSB)-1.5d0)
     .   +VALPH(6)**2*(dlog(Max(MZ2,VALPH(6))/QSTSB)-1.5d0))/64d0/Pi**2)
        endif

      IF(Vvu.lt.VEW)PROB(28)=Max(PROB(28),
     .                      DDIM(-1d-2,(VEW-Vvu)/Max(dabs(VEW),1d-10)))


c      Minimum with <S>=<Hu>=0

      Taux=MHuS+MHdS+g2q/2d0*vdq**2
      Daux=MHuS-MHdS-g1q/2d0*vdq**2
      aux=M3HQ**2+(l*XIFQ)**2+2d0*M3HQ*l*XIFQ*DDCOS(phi3q+phiF-phi01)
      MHD1aux=(Taux-dsqrt(Daux**2+4d0*aux))/2d0
      MHD2aux=(Taux+dsqrt(Daux**2+4d0*aux))/2d0

      MH2(1,1)=MHuS+(l**2-(g1q+g2q)/4d0)*vdq**2
      MH2(1,2)=-M3HQ*DDCOS(phi3q)-l*XIFQ*DDCOS(phi01-phiF)
      MH2(2,1)=MH2(1,2)
      MH2(1,3)=-l*(Alcos1+MUPQ*DDCOS(phi01-phiP))*vdq
      MH2(3,1)=MH2(1,3)
      MH2(2,2)=MHdS+3d0*(g1q+g2q)/4d0*vdq**2
      MH2(2,3)=0d0
      MH2(3,2)=MH2(2,3)
      MH2(3,3)=MSS+MUPQ**2+(l*vdq)**2+MSPQ*DDCOS(phiSPq)
     .                               +2d0*k*XIFQ*DDCOS(phi02-phiF)
      MH2(1,4)=0d0
      MH2(4,1)=MH2(1,4)
      MH2(2,4)=M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-phiF)
      MH2(4,2)=MH2(2,4)
      MH2(3,4)=l*(IAl+MUPQ*DDSIN(phi01-phiP))*vdq
      MH2(4,3)=MH2(3,4)
      MH2(4,4)=MHuS+(l**2-(g1q+g2q)/4d0)*vdq**2
      MH2(1,5)=M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-phiF)
      MH2(5,1)=MH2(1,5)
      MH2(2,5)=0d0
      MH2(5,2)=MH2(2,5)
      MH2(3,5)=0d0
      MH2(5,3)=MH2(3,5)
      MH2(4,5)=M3HQ*DDCOS(phi3q)+l*XIFQ*DDCOS(phi01-phiF)
      MH2(5,4)=MH2(4,5)
      MH2(5,5)=MHdS+(g1q+g2q)/4d0*vdq**2
      MH2(1,6)=l*(IAl+MUPQ*DDSIN(phi01-phiP))*vdq
      MH2(6,1)=MH2(1,6)
      MH2(2,6)=0d0
      MH2(6,2)=MH2(2,6)
      MH2(3,6)=-MSPQ*DDSIN(phiSPq)-2d0*k*XIFQ*DDSIN(phi02-phiF)
      MH2(6,3)=MH2(3,6)
      MH2(4,6)=l*(Alcos1+MUPQ*DDCOS(phi01-phiP))*vdq
      MH2(6,4)=MH2(4,6)
      MH2(5,6)=0d0
      MH2(6,5)=MH2(4,5)
      MH2(6,6)=MSS+MUPQ**2+l**2*vdq**2-MSPQ*DDCOS(phiSPq)
     .                               -2d0*k*XIFQ*DDCOS(phi02-phiF)

      CALL DIAGN(6,MH2,VALPH,VECPH,1d-10)

        MSTaux1=MSQ3-(g1q/3d0-g2q)/4d0*vdq**2
        MSTaux2=MSU3+g1q/3d0*vdq**2
        MSBaux1=(MSQ3+MSD3+2d0*mbotq**2-(g1q+g2q)/4d0*vdq**2
     . -dsqrt((MSQ3-MSD3+(-g1q/3d0+g2q)/4d0*vdq**2)**2
     .        +4d0*(Ybq*vdq*ABP)**2))/2d0
        MSBaux2=(MSQ3+MSD3+2d0*mbotq**2-(g1q+g2q)/4d0*vdq**2
     . +dsqrt((MSQ3-MSD3+(-g1q/3d0+g2q)/4d0*vdq**2)**2
     .        +4d0*(Ybq*vdq*ABP)**2))/2d0
        MSLaux1=(MSL3+MSE3+2d0*(mtau*vdq/vd)**2-(g1q+g2q)/4d0*vdq**2
     . -dsqrt((MSL3-MSE3-(-3d0/2d0*g1q+g2q)/4d0*vdq**2)**2
     .        +4d0*(mtau/vd*vdq*ATAU)**2))/2d0
        MSLaux2=(MSL3+MSE3+2d0*(mtau*vdq/vd)**2-(g1q+g2q)/4d0*vdq**2
     . +dsqrt((MSL3-MSE3-(-3d0/2d0*g1q+g2q)/4d0*vdq**2)**2
     .        +4d0*(mtau/vd*vdq*ATAU)**2))/2d0
        MSNTaux=MSL3+(g1q+g2q)/4d0*vdq**2
        MSUaux1=MSQ1-(g1q/3d0-g2q)/4d0*vdq**2
        MSUaux2=MSU1+g1q/3d0*vdq**2
        MSDaux1=MSQ1-(g1q/3d0+g2q)/4d0*vdq**2
        MSDaux2=MSD1-g1q/6d0*vdq**2
        MSEaux1=MSL1-(-g1q+g2q)/4d0*vdq**2
        MSEaux2=MSE1-g1q/2d0*vdq**2
        MSNEaux=MSL1+(g1q+g2q)/4d0*vdq**2

      Vvd=MHdS*vdq**2+(G1Q+G2Q)/8d0*vdq**4
     .    -(3d0*mbotq**4*(dlog(mbotq**2/QSTSB)-1.5d0)
     .     +(mtau*vdq/vd)**4*(dlog(mtau**2/QSTSB)-1.5d0))/16d0/Pi**2
     .    +3d0*(2d0*(g2q*vdq**2/2d0)**2
     .           *(dlog(g2q*vdq**2/2d0/QSTSB)-1.5d0)
     .         +((g1q+g2q)/2d0*vdq**2)**2
     .           *(dlog((g1q+g2q)/2d0*vdq**2/QSTSB)-1.5d0))/64d0/Pi**2
     .    +(3d0*MSTaux1**2*(dlog(Max(MSTaux1,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSTaux2**2*(dlog(Max(MSTaux2,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSBaux1**2*(dlog(Max(MSBaux1,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSBaux2**2*(dlog(Max(MSBaux2,MZ2)/QSTSB)-1.5d0)
     .     +MSLaux1**2*(dlog(Max(MSLaux1,MZ2)/QSTSB)-1.5d0)
     .     +MSLaux2**2*(dlog(Max(MSLaux2,MZ2)/QSTSB)-1.5d0)
     .     +MSNTaux**2*(dlog(Max(MSNTaux,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSUaux1**2*(dlog(Max(MSUaux1,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSUaux2**2*(dlog(Max(MSUaux2,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSDaux1**2*(dlog(Max(MSDaux1,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSDaux2**2*(dlog(Max(MSDaux2,MZ2)/QSTSB)-1.5d0)
     .     +2d0*MSEaux1**2*(dlog(Max(MSEaux1,MZ2)/QSTSB)-1.5d0)
     .     +2d0*MSEaux2**2*(dlog(Max(MSEaux2,MZ2)/QSTSB)-1.5d0)
     . +2d0*MSNEaux**2*(dlog(Max(MSNEaux,MZ2)/QSTSB)-1.5d0))/32d0/Pi**2
     .    +3d0/512d0/Pi**4*Ybq**4*(dlog(QSTSB/mtopq**2)**2
     .      *(64d0*Pi*ALSQ-2d0/3d0*g1q+3d0*sinb**2*Ytq**2
     .       -3d0*cosb**2*Ybq**2)+
     .       (dlog(MA2/mtopq**2)**2-dlog(QSTSB/mtopq**2)**2)
     .       *(3d0*sinb**2*Ybq**2+(3d0*sinb**2+1d0)*Ytq**2))*vdq**4
     .    -((M2r**2+g2q*vdq**2)**2
     .             *(dlog((M2r**2+g2q*vdq**2)/QSTSB)-1.5d0))/16d0/Pi**2
     .    -((M1r**2+g1q*vdq**2)**2
     .                   *(dlog((M1r**2+g1q*vdq**2)/QSTSB)-1.5d0)
     .     +(M2r**2+g2q*vdq**2)**2
     .                   *(dlog((M2r**2+g2q*vdq**2)/QSTSB)-1.5d0)
     .     +(mupsi**2+2d0*(l*vdq)**2)**2*
     . (dlog(Max(mupsi**2+2d0*(l*vuq)**2,1d0)/QSTSB)-1.5d0))/32d0/Pi**2
     .     +(MHD1aux**2*(dlog(Max(MHD1aux,MW2)/QSTSB)-1.5d0)
     .     +MHD2aux**2*(dlog(Max(MHD2aux,MW2)/QSTSB)-1.5d0))/32d0/Pi**2
     .    +(VALPH(1)**2*(dlog(Max(MZ2,VALPH(1))/QSTSB)-1.5d0)
     .     +VALPH(2)**2*(dlog(Max(MZ2,VALPH(2))/QSTSB)-1.5d0)
     .     +VALPH(3)**2*(dlog(Max(MZ2,VALPH(3))/QSTSB)-1.5d0)
     .     +VALPH(4)**2*(dlog(Max(MZ2,VALPH(4))/QSTSB)-1.5d0)
     .     +VALPH(5)**2*(dlog(Max(MZ2,VALPH(5))/QSTSB)-1.5d0)
     .   +VALPH(6)**2*(dlog(Max(MZ2,VALPH(6))/QSTSB)-1.5d0))/64d0/Pi**2

        If(MHdS.lt.0d0)then
      aux=dsqrt(-4d0*MHdS/(G1Q+G2Q))

      Taux=MHuS-MHdS-g1q/2d0*aux**2
      Daux=M3HQ**2+(l*XIFQ)**2+2d0*M3HQ*l*XIFQ*DDCOS(phi3q+phiF-phi01)
      MHD1aux=(MHuS+MHdS+g2q/2d0*aux**2-dsqrt(Taux**2+4d0*Daux))/2d0
      MHD2aux=(MHuS+MHdS+g2q/2d0*aux**2+dsqrt(Taux**2+4d0*Daux))/2d0

      MH2(1,1)=MHuS+(l**2-(g1q+g2q)/4d0)*aux**2
      MH2(1,2)=-M3HQ*DDCOS(phi3q)-l*XIFQ*DDCOS(phi01-phiF)
      MH2(2,1)=MH2(1,2)
      MH2(1,3)=-l*(Alcos1+MUPQ*DDCOS(phi01-phiP))*aux
      MH2(3,1)=MH2(1,3)
      MH2(2,2)=MHdS+3d0*(g1q+g2q)/4d0*aux**2
      MH2(2,3)=0d0
      MH2(3,2)=MH2(2,3)
      MH2(3,3)=MSS+MUPQ**2+(l*aux)**2+MSPQ*DDCOS(phiSPq)
     .                               +2d0*k*XIFQ*DDCOS(phi02-phiF)
      MH2(1,4)=0d0
      MH2(4,1)=MH2(1,4)
      MH2(2,4)=M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-phiF)
      MH2(4,2)=MH2(2,4)
      MH2(3,4)=l*(IAl+MUPQ*DDSIN(phi01-phiP))*aux
      MH2(4,3)=MH2(3,4)
      MH2(4,4)=MHuS+(l**2-(g1q+g2q)/4d0)*aux**2
      MH2(1,5)=M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-phiF)
      MH2(5,1)=MH2(1,5)
      MH2(2,5)=0d0
      MH2(5,2)=MH2(2,5)
      MH2(3,5)=0d0
      MH2(5,3)=MH2(3,5)
      MH2(4,5)=M3HQ*DDCOS(phi3q)+l*XIFQ*DDCOS(phi01-phiF)
      MH2(5,4)=MH2(4,5)
      MH2(5,5)=MHdS+(g1q+g2q)/4d0*aux**2
      MH2(1,6)=l*(IAl+MUPQ*DDSIN(phi01-phiP))*aux
      MH2(6,1)=MH2(1,6)
      MH2(2,6)=0d0
      MH2(6,2)=MH2(2,6)
      MH2(3,6)=-MSPQ*DDSIN(phiSPq)-2d0*k*XIFQ*DDSIN(phi02-phiF)
      MH2(6,3)=MH2(3,6)
      MH2(4,6)=l*(Alcos1+MUPQ*DDCOS(phi01-phiP))*aux
      MH2(6,4)=MH2(4,6)
      MH2(5,6)=0d0
      MH2(6,5)=MH2(4,5)
      MH2(6,6)=MSS+MUPQ**2+l**2*aux**2-MSPQ*DDCOS(phiSPq)
     .                               -2d0*k*XIFQ*DDCOS(phi02-phiF)

      CALL DIAGN(6,MH2,VALPH,VECPH,1d-10)

      MSTaux1=MSQ3-(g1q/3d0-g2q)/4d0*aux**2
        MSTaux2=MSU3+g1q/3d0*aux**2
        MSBaux1=(MSQ3+MSD3+2d0*(Ybq*aux)**2-(g1q+g2q)/4d0*aux**2
     . -dsqrt((MSQ3-MSD3+(-g1q/3d0+g2q)/4d0*aux**2)**2
     .        +4d0*(Ybq*aux*ABP)**2))/2d0
        MSBaux2=(MSQ3+MSD3+2d0*(Ybq*aux)**2-(g1q+g2q)/4d0*aux**2
     . +dsqrt((MSQ3-MSD3+(-g1q/3d0+g2q)/4d0*aux**2)**2
     .        +4d0*(Ybq*aux*ABP)**2))/2d0
        MSLaux1=(MSL3+MSE3+2d0*(mtau*aux/vd)**2-(g1q+g2q)/4d0*aux**2
     . -dsqrt((MSL3-MSE3-(-3d0/2d0*g1q+g2q)/4d0*aux**2)**2
     .        +4d0*(mtau/vd*aux*ATAU)**2))/2d0
        MSLaux2=(MSL3+MSE3+2d0*(mtau*aux/vd)**2-(g1q+g2q)/4d0*aux**2
     . +dsqrt((MSL3-MSE3-(-3d0/2d0*g1q+g2q)/4d0*aux**2)**2
     .        +4d0*(mtau/vd*aux*ATAU)**2))/2d0
        MSNTaux=MSL3+(g1q+g2q)/4d0*aux**2
        MSUaux1=MSQ1-(g1q/3d0-g2q)/4d0*aux**2
        MSUaux2=MSU1+g1q/3d0*aux**2
        MSDaux1=MSQ1-(g1q/3d0+g2q)/4d0*aux**2
        MSDaux2=MSD1-g1q/6d0*aux**2
        MSEaux1=MSL1-(-g1q+g2q)/4d0*aux**2
        MSEaux2=MSE1-g1q/2d0*aux**2
        MSNEaux=MSL1+(g1q+g2q)/4d0*aux**2

      Vvd=min(Vvd,MHdS*aux**2+(G1Q+G2Q)/8d0*aux**4
     .    -(3d0*(Ybq*aux)**4*(dlog((Ybq*aux)**2/QSTSB)-1.5d0)
     .     +(mtau*aux/vd)**4*(dlog((mtau*aux/vd)**2/QSTSB)-1.5d0)
     .                                                    )/16d0/Pi**2
     .    +3d0*(2d0*(g2q*aux**2/2d0)**2
     .           *(dlog(g2q*aux**2/2d0/QSTSB)-1.5d0)
     .         +((g1q+g2q)/2d0*aux**2)**2
     .           *(dlog((g1q+g2q)/2d0*aux**2/QSTSB)-1.5d0))/64d0/Pi**2
     .    +(3d0*MSTaux1**2*(dlog(Max(MSTaux1,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSTaux2**2*(dlog(Max(MSTaux2,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSBaux1**2*(dlog(Max(MSBaux1,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSBaux2**2*(dlog(Max(MSBaux2,MZ2)/QSTSB)-1.5d0)
     .     +MSLaux1**2*(dlog(Max(MSLaux1,MZ2)/QSTSB)-1.5d0)
     .     +MSLaux2**2*(dlog(Max(MSLaux2,MZ2)/QSTSB)-1.5d0)
     .     +MSNTaux**2*(dlog(Max(MSNTaux,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSUaux1**2*(dlog(Max(MSUaux1,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSUaux2**2*(dlog(Max(MSUaux2,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSDaux1**2*(dlog(Max(MSDaux1,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSDaux2**2*(dlog(Max(MSDaux2,MZ2)/QSTSB)-1.5d0)
     .     +2d0*MSEaux1**2*(dlog(Max(MSEaux1,MZ2)/QSTSB)-1.5d0)
     .     +2d0*MSEaux2**2*(dlog(Max(MSEaux2,MZ2)/QSTSB)-1.5d0)
     . +2d0*MSNEaux**2*(dlog(Max(MSNEaux,MZ2)/QSTSB)-1.5d0))/32d0/Pi**2
     .    +3d0/512d0/Pi**4*Ybq**4*(dlog(QSTSB/mtopq**2)**2
     .      *(64d0*Pi*ALSQ-2d0/3d0*g1q+3d0*sinb**2*Ytq**2
     .       -3d0*cosb**2*Ybq**2)+
     .       (dlog(MA2/mtopq**2)**2-dlog(QSTSB/mtopq**2)**2)
     .       *(3d0*sinb**2*Ybq**2+(3d0*sinb**2+1d0)*Ytq**2))*aux**4
     .    -((M2r**2+g2q*aux**2)**2
     .             *(dlog((M2r**2+g2q*aux**2)/QSTSB)-1.5d0))/16d0/Pi**2
     .    -((M1r**2+g1q*aux**2)**2
     .                   *(dlog((M1r**2+g1q*aux**2)/QSTSB)-1.5d0)
     .     +(M2r**2+g2q*aux**2)**2
     .                   *(dlog((M2r**2+g2q*aux**2)/QSTSB)-1.5d0)
     .     +(mupsi**2+2d0*(l*aux)**2)**2*
     . (dlog(Max(mupsi**2+2d0*(l*vuq)**2,1d0)/QSTSB)-1.5d0))/32d0/Pi**2
     .     +(MHD1aux**2*(dlog(Max(MHD1aux,MW2)/QSTSB)-1.5d0)
     .     +MHD2aux**2*(dlog(Max(MHD2aux,MW2)/QSTSB)-1.5d0))/32d0/Pi**2
     .    +(VALPH(1)**2*(dlog(Max(MZ2,VALPH(1))/QSTSB)-1.5d0)
     .     +VALPH(2)**2*(dlog(Max(MZ2,VALPH(2))/QSTSB)-1.5d0)
     .     +VALPH(3)**2*(dlog(Max(MZ2,VALPH(3))/QSTSB)-1.5d0)
     .     +VALPH(4)**2*(dlog(Max(MZ2,VALPH(4))/QSTSB)-1.5d0)
     .     +VALPH(5)**2*(dlog(Max(MZ2,VALPH(5))/QSTSB)-1.5d0)
     .   +VALPH(6)**2*(dlog(Max(MZ2,VALPH(6))/QSTSB)-1.5d0))/64d0/Pi**2)
        endif

      IF(Vvd.lt.VEW)PROB(28)=Max(PROB(28),
     .                      DDIM(-1d-2,(Vvd-VEW)/Max(dabs(VEW),1d-10)))


c      Phase variations

      DO I=1,Nphi
      DO J=1,Nphi

      phiS=2d0*Pi*I/DFLOAT(Nphi+1)
      phiH=2d0*Pi*J/DFLOAT(Nphi+1)

        MSTaux1=(MSQ3+MSU3+2d0*mtopq**2-(g1q+g2q)/4d0*(vuq**2-vdq**2)
     . -dsqrt((MSQ3-MSU3+(5d0/3d0*g1q-g2q)/4d0*(vuq**2-vdq**2))**2
     .        +4d0*Ytq**2*(ATP**2*vuq**2+muq**2*vdq**2
     .              -2d0*ATP*muq*vuq*vdq*DDCOS(phiAT+phiS+phiH))))/2d0
        MSTaux2=(MSQ3+MSU3+2d0*mtopq**2-(g1q+g2q)/4d0*(vuq**2-vdq**2)
     . +dsqrt((MSQ3-MSU3+(5d0/3d0*g1q-g2q)/4d0*(vuq**2-vdq**2))**2
     .        +4d0*Ytq**2*(ATP**2*vuq**2+muq**2*vdq**2
     .              -2d0*ATP*muq*vuq*vdq*DDCOS(phiAT+phiS+phiH))))/2d0
        MSBaux1=(MSQ3+MSD3+2d0*mbotq**2-(g1q+g2q)/4d0*(vdq**2-vuq**2)
     . -dsqrt((MSQ3-MSD3+(-g1q/3d0+g2q)/4d0*(vdq**2-vuq**2))**2
     .        +4d0*Ybq**2*(ABP**2*vdq**2+muq**2*vuq**2
     .              -2d0*ABP*muq*vuq*vdq*DDCOS(phiAB+phiS+phiH))))/2d0
        MSBaux2=(MSQ3+MSD3+2d0*mbotq**2-(g1q+g2q)/4d0*(vdq**2-vuq**2)
     . +dsqrt((MSQ3-MSD3+(-g1q/3d0+g2q)/4d0*(vdq**2-vuq**2))**2
     .        +4d0*Ybq**2*(ABP**2*vdq**2+muq**2*vuq**2
     .              -2d0*ABP*muq*vuq*vdq*DDCOS(phiAB+phiS+phiH))))/2d0
        MSLaux1=(MSL3+MSE3+2d0*(mtau*vdq/vd)**2
     .                                 -(g1q+g2q)/4d0*(vdq**2-vuq**2)
     . -dsqrt((MSL3-MSE3-(-3d0/2d0*g1q+g2q)/4d0*(vdq**2-vuq**2))**2
     .        +4d0*(mtau/vd)**2*(ATAU**2*vdq**2+muq**2*vuq**2
     .           -2d0*ATAU*muq*vuq*vdq*DDCOS(PhiATAU+PhiS+PhiH))))/2d0
        MSLaux2=(MSL3+MSE3+2d0*(mtau*vdq/vd)**2
     .                                 -(g1q+g2q)/4d0*(vdq**2-vuq**2)
     . +dsqrt((MSL3-MSE3-(-3d0/2d0*g1q+g2q)/4d0*(vdq**2-vuq**2))**2
     .        +4d0*(mtau/vd)**2*(ATAU**2*vdq**2+muq**2*vuq**2
     .           -2d0*ATAU*muq*vuq*vdq*DDCOS(PhiATAU+PhiS+PhiH))))/2d0
        MSNTaux=MSNT2

      Taux=M2r**2+mur**2+g2q*(vuq**2+vdq**2)
      Daux=M2r**2-mur**2+g2q*(vuq**2-vdq**2)
      aux=g2q*((M2r*vuq)**2+(mur*vdq)**2
     .         +2d0*mur*M2r*vuq*vdq*DDCOS(phiM2-phi01-phiS-phiH))
        MCH1aux=(Taux-dsqrt(Daux**2+4d0*aux))/2d0
        MCH2aux=(Taux+dsqrt(Daux**2+4d0*aux))/2d0

      MN2(1,1)=M1r**2+g1q/2d0*(vuq**2+vdq**2)
      MN2(6,6)=M1r**2+g1q/2d0*(vuq**2+vdq**2)
      MN2(6,1)=0d0
      MN2(1,6)=0d0
      MN2(1,2)=-dsqrt(g1q*g2q)/2d0*(vuq**2+vdq**2)
      MN2(2,1)=-dsqrt(g1q*g2q)/2d0*(vuq**2+vdq**2)
      MN2(6,7)=-dsqrt(g1q*g2q)/2d0*(vuq**2+vdq**2)
      MN2(7,6)=-dsqrt(g1q*g2q)/2d0*(vuq**2+vdq**2)
      MN2(1,7)=0d0
      MN2(2,6)=0d0
      MN2(7,1)=0d0
      MN2(6,2)=0d0
      MN2(1,3)=dsqrt(g1q/2d0)*(M1r*vuq*DDCOS(PhiM1+phiH/2d0)
     .                           +mur*vdq*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(3,1)=dsqrt(g1q/2d0)*(M1r*vuq*DDCOS(PhiM1)
     .                           +mur*vdq*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(6,8)=dsqrt(g1q/2d0)*(M1r*vuq*DDCOS(PhiM1)
     .                           +mur*vdq*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(8,6)=dsqrt(g1q/2d0)*(M1r*vuq*DDCOS(PhiM1)
     .                           +mur*vdq*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(1,8)=dsqrt(g1q/2d0)*(-M1r*vuq*DDSIN(PhiM1+phiH/2d0)
     .                           +mur*vdq*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(3,6)=dsqrt(g1q/2d0)*(M1r*vuq*DDSIN(PhiM1+phiH/2d0)
     .                           -mur*vdq*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(8,1)=dsqrt(g1q/2d0)*(-M1r*vuq*DDSIN(PhiM1+phiH/2d0)
     .                           +mur*vdq*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(6,3)=dsqrt(g1q/2d0)*(M1r*vuq*DDSIN(PhiM1+phiH/2d0)
     .                           -mur*vdq*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(1,4)=-dsqrt(g1q/2d0)*(M1r*vdq*DDCOS(PhiM1+phiH/2d0)
     .                           +mur*vuq*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(4,1)=-dsqrt(g1q/2d0)*(M1r*vdq*DDCOS(PhiM1+phiH/2d0)
     .                           +mur*vuq*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(6,9)=-dsqrt(g1q/2d0)*(M1r*vdq*DDCOS(PhiM1+phiH/2d0)
     .                           +mur*vuq*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(9,6)=-dsqrt(g1q/2d0)*(M1r*vdq*DDCOS(PhiM1+phiH/2d0)
     .                           +mur*vuq*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(1,9)=dsqrt(g1q/2d0)*(M1r*vdq*DDSIN(PhiM1+phiH/2d0)
     .                           -mur*vuq*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(9,1)=dsqrt(g1q/2d0)*(M1r*vdq*DDSIN(PhiM1+phiH/2d0)
     .                           -mur*vuq*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(4,6)=dsqrt(g1q/2d0)*(-M1r*vdq*DDSIN(PhiM1+phiH/2d0)
     .                           +mur*vuq*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(6,4)=dsqrt(g1q/2d0)*(-M1r*vdq*DDSIN(PhiM1+phiH/2d0)
     .                           +mur*vuq*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(5,1)=0d0
      MN2(1,5)=0d0
      MN2(6,10)=0d0
      MN2(10,6)=0d0
      MN2(10,1)=0d0
      MN2(1,10)=0d0
      MN2(10,1)=0d0
      MN2(1,10)=0d0
      MN2(6,5)=0d0
      MN2(5,6)=0d0
      MN2(2,2)=M2r**2+g2q/2d0*(vuq**2+vdq**2)
      MN2(7,7)=M2r**2+g2q/2d0*(vuq**2+vdq**2)
      MN2(2,7)=0d0
      MN2(7,2)=0d0
      MN2(2,3)=-dsqrt(g2q/2d0)*(M2r*vuq*DDCOS(PhiM2+phiH/2d0)
     .                           +mur*vd*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(3,2)=-dsqrt(g2q/2d0)*(M2r*vuq*DDCOS(PhiM2+phiH/2d0)
     .                           +mur*vd*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(7,8)=-dsqrt(g2q/2d0)*(M2r*vuq*DDCOS(PhiM2+phiH/2d0)
     .                           +mur*vd*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(8,7)=-dsqrt(g2q/2d0)*(M2r*vuq*DDCOS(PhiM2+phiH/2d0)
     .                           +mur*vd*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(2,8)=dsqrt(g2q/2d0)*(M2r*vuq*DDSIN(PhiM2+phiH/2d0)
     .                           -mur*vd*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(8,2)=dsqrt(g2q/2d0)*(M2r*vuq*DDSIN(PhiM2+phiH/2d0)
     .                           -mur*vd*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(3,7)=dsqrt(g2q/2d0)*(-M2r*vuq*DDSIN(PhiM2+phiH/2d0)
     .                           +mur*vd*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(7,3)=dsqrt(g2q/2d0)*(-M2r*vuq*DDSIN(PhiM2+phiH/2d0)
     .                           +mur*vd*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(2,4)=dsqrt(g2q/2d0)*(M2r*vdq*DDCOS(PhiM2+phiH/2d0)
     .                           +mur*vu*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(4,2)=dsqrt(g2q/2d0)*(M2r*vdq*DDCOS(PhiM2+phiH/2d0)
     .                           +mur*vu*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(7,9)=dsqrt(g2q/2d0)*(M2r*vdq*DDCOS(PhiM2+phiH/2d0)
     .                           +mur*vu*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(9,7)=dsqrt(g2q/2d0)*(M2r*vdq*DDCOS(PhiM2+phiH/2d0)
     .                           +mur*vu*DDCOS(Phi01+phiS+phiH/2d0))
      MN2(2,9)=dsqrt(g2q/2d0)*(-M2r*vdq*DDSIN(PhiM2+phiH/2d0)
     .                           +mur*vu*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(9,2)=dsqrt(g2q/2d0)*(-M2r*vdq*DDSIN(PhiM2+phiH/2d0)
     .                           +mur*vu*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(7,4)=dsqrt(g2q/2d0)*(M2r*vdq*DDSIN(PhiM2+phiH/2d0)
     .                           -mur*vu*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(4,7)=dsqrt(g2q/2d0)*(M2r*vdq*DDSIN(PhiM2+phiH/2d0)
     .                           -mur*vu*DDSIN(Phi01+phiS+phiH/2d0))
      MN2(2,5)=0d0
      MN2(5,2)=0d0
      MN2(7,10)=0d0
      MN2(10,7)=0d0
      MN2(2,10)=0d0
      MN2(10,2)=0d0
      MN2(5,7)=0d0
      MN2(7,5)=0d0
      MN2(3,3)=mur**2+l**2*vdq**2+(g1q+g2q)/2d0*vuq**2
      MN2(8,8)=mur**2+l**2*vdq**2+(g1q+g2q)/2d0*vuq**2
      MN2(3,8)=0d0
      MN2(8,3)=0d0
      MN2(3,4)=(l**2-(g1q+g2q)/2d0)*vuq*vdq
      MN2(4,3)=(l**2-(g1q+g2q)/2d0)*vuq*vdq
      MN2(8,9)=(l**2-(g1q+g2q)/2d0)*vuq*vdq
      MN2(9,8)=(l**2-(g1q+g2q)/2d0)*vuq*vdq
      MN2(3,9)=0d0
      MN2(9,3)=0d0
      MN2(4,8)=0d0
      MN2(8,4)=0d0
      MN2(3,5)=l*mur*vuq*DDCOS(phiS-phiH/2d0)
     .        -l*vdq*(mupsi*DDCOS(Phi01-phip+phiH/2d0)
     .               +ks2si*DDCOS(Phi01-Phi02+phiH/2d0-phiS))
      MN2(5,3)=l*mur*vuq*DDCOS(phiS-phiH/2d0)
     .        -l*vdq*(mupsi*DDCOS(Phi01-phip+phiH/2d0)
     .               +ks2si*DDCOS(Phi01-Phi02+phiH/2d0-phiS))
      MN2(8,10)=l*mur*vuq*DDCOS(phiS-phiH/2d0)
     .        -l*vdq*(mupsi*DDCOS(Phi01-phip+phiH/2d0)
     .               +ks2si*DDCOS(Phi01-Phi02+phiH/2d0-phiS))
      MN2(10,8)=l*mur*vuq*DDCOS(phiS-phiH/2d0)
     .        -l*vdq*(mupsi*DDCOS(Phi01-phip+phiH/2d0)
     .               +ks2si*DDCOS(Phi01-Phi02+phiH/2d0-phiS))
      MN2(3,10)=-l*mur*vuq*DDSIN(phiS-phiH/2d0)
     .          +l*vdq*(mupsi*DDSIN(Phi01-phip+phiH/2d0)
     .                 +ks2si*DDSIN(Phi01-Phi02+phiH/2d0-phiS))
      MN2(10,3)=-l*mur*vuq*DDSIN(phiS-phiH/2d0)
     .          +l*vdq*(mupsi*DDSIN(Phi01-phip+phiH/2d0)
     .                 +ks2si*DDSIN(Phi01-Phi02+phiH/2d0-phiS))
      MN2(5,8)=l*mur*vuq*DDSIN(phiS-phiH/2d0)
     .         -l*vdq*(mupsi*DDSIN(Phi01-phip+phiH/2d0)
     .                 +ks2si*DDSIN(Phi01-Phi02+phiH/2d0-phiS))
      MN2(8,5)=l*mur*vuq*DDSIN(phiS-phiH/2d0)
     .         -l*vdq*(mupsi*DDSIN(Phi01-phip+phiH/2d0)
     .                 +ks2si*DDSIN(Phi01-Phi02+phiH/2d0-phiS))
      MN2(4,4)=mur**2+l**2*vuq**2+(g1q+g2q)/2d0*vdq**2
      MN2(9,9)=mur**2+l**2*vuq**2+(g1q+g2q)/2d0*vdq**2
      MN2(4,9)=0d0
      MN2(9,4)=0d0
      MN2(4,5)=l*mur*vdq*DDCOS(phiS-phiH/2d0)
     .        -l*vuq*(mupsi*DDCOS(Phi01-phip+phiH/2d0)
     .               +ks2si*DDCOS(Phi01-Phi02+phiH/2d0-phiS))
      MN2(5,4)=l*mur*vdq*DDCOS(phiS-phiH/2d0)
     .        -l*vuq*(mupsi*DDCOS(Phi01-phip+phiH/2d0)
     .               +ks2si*DDCOS(Phi01-Phi02+phiH/2d0-phiS))
      MN2(9,10)=l*mur*vdq*DDCOS(phiS-phiH/2d0)
     .        -l*vuq*(mupsi*DDCOS(Phi01-phip+phiH/2d0)
     .               +ks2si*DDCOS(Phi01-Phi02+phiH/2d0-phiS))
      MN2(10,9)=l*mur*vdq*DDCOS(phiS-phiH/2d0)
     .        -l*vuq*(mupsi*DDCOS(Phi01-phip+phiH/2d0)
     .               +ks2si*DDCOS(Phi01-Phi02+phiH/2d0-phiS))
      MN2(4,10)=-l*mur*vdq*DDSIN(phiS-phiH/2d0)
     .          +l*vuq*(mupsi*DDSIN(Phi01-phip+phiH/2d0)
     .               +ks2si*DDSIN(Phi01-Phi02+phiH/2d0-phiS))
      MN2(10,4)=-l*mur*vdq*DDSIN(phiS-phiH/2d0)
     .          +l*vuq*(mupsi*DDSIN(Phi01-phip+phiH/2d0)
     .               +ks2si*DDSIN(Phi01-Phi02+phiH/2d0-phiS))
      MN2(5,9)=l*mur*vdq*DDSIN(phiS-phiH/2d0)
     .         -l*vuq*(mupsi*DDSIN(Phi01-phip+phiH/2d0)
     .               +ks2si*DDSIN(Phi01-Phi02+phiH/2d0-phiS))
      MN2(9,5)=l*mur*vdq*DDSIN(phiS-phiH/2d0)
     .         -l*vuq*(mupsi*DDSIN(Phi01-phip+phiH/2d0)
     .               +ks2si*DDSIN(Phi01-Phi02+phiH/2d0-phiS))
      MN2(5,5)=mupsi**2+ks2si**2+2d0*mupsi*ks2si*DDCOS(phip-phi02-phiS)
     .        +l**2*(vuq**2+vdq**2)
      MN2(10,10)=mupsi**2+ks2si**2+2d0*mupsi*ks2si
     .        *DDCOS(phip-phi02-phiS)+l**2*(vuq**2+vdq**2)
      MN2(5,10)=0d0
      MN2(10,5)=0d0
      CALL DIAGN(10,MN2,VALPN,VECPN,1d-10)
      CALL SORTNA(10,VALPN,VECPN)

      Taux=MHuS+MHdS+2d0*muq**2+g2q/2d0*(vuq**2+vdq**2)
      Daux=MHuS-MHdS+g1q/2d0*(vuq**2-vdq**2)
      aux=M3HQ**2+(l*XIFQ)**2+2d0*M3HQ*l*XIFQ*DDCOS(phi3q+phiF-phi01)
     .   +(MUPQ*muq)**2+2d0*MUPQ*muq*(M3HQ*DDCOS(phi3q+phiP-phi01+phiS)
     .                                    +l*XIFQ*DDCOS(phiP-phiF+phiS))
     .   +muq**2*(Alcos1**2+IAl**2+(k/l*muq)**2+2d0*k/l*muq
     .       *(Alcos1*DDCOS(phi0-3d0*phiS)+IAl*DDSIN(phi0-3d0*phiS)))
     .   +(l**2-g2q/2d0)**2*(vuq*vdq)**2
     . -2d0*(l**2-g2q/2d0)*muq*vuq*vdq*(Alcos1*DDCOS(phiH+phiS)
     .       -IAl*DDSIN(phiH+phiS)+k/l*muq*DDCOS(phi0+phiH-2d0*phiS))
     . +2d0*M3HQ*(muq*(Alcos1*DDCOS(phi3q-phiS)+IAL*DDSIN(phi3q-phiS))
     .           +k/l*muq**2*DDCOS(phi3q-phi0+2d0*phiS)
     .           -(l**2-g2q/2d0)*vuq*vdq*DDCOS(phi3q+phiH))
     . +2d0*l*XIFQ*(muq*(Alcos1*DDCOS(phi01-phiF-phiS)
     .                       +IAL*DDSIN(phi01-phiF-phiS))
     .           +k/l*muq**2*DDCOS(phi02-phiF+2d0*phiS)
     .           -(l**2-g2q/2d0)*vuq*vdq*DDCOS(phi01-phiF+phiH))
     . +2d0*MUPQ*muq*(muq*(Alcos1*DDCOS(phi01-phiP-2d0*phiS)
     .                       +IAL*DDSIN(phi01-phiP-2d0*phiS))
     .           +k/l*muq**2*DDCOS(phi02-phiP+3d0*phiS)
     .           -(l**2-g2q/2d0)*vuq*vdq*DDCOS(phi01-phiP+phiS+phiH))
      MHD1aux=(Taux-dsqrt(Daux**2+4d0*aux))/2d0
      MHD2aux=(Taux+dsqrt(Daux**2+4d0*aux))/2d0

      MH2(1,1)=MHuS+muq**2+(l**2-(g1q+g2q)/4d0)*vdq**2
     .        +(g1q+g2q)/4d0*vuq**2*(1d0+DDCOS(phiH/2d0)**2)
      MH2(1,2)=-muq*(Alcos1*DDCOS(phiS)-IAl*DDSIN(phiS))
     .        -k/l*muq**2*DDCOS(phi0-2d0*phiS)
     .        +2d0*(l**2-(g1q+g2q)/4d0)*vuq*vdq*DDCOS(phiH/2d0)**2
     .        -M3HQ*DDCOS(phi3q)-l*XIFQ*DDCOS(phi01-phiF)
     .        -muq*MUPQ*DDCOS(phi01-phiP-phiS)
      MH2(2,1)=MH2(1,2)
      MH2(1,3)=-l*vdq*(Alcos1*DDCOS(phiH/2d0)-IAl*DDSIN(phiH/2d0)
     .                +2d0*k/l*muq*DDCOS(phi0-phiS+phiH/2d0)
     .                +MUPQ*DDCOS(phi01-phiP+phiH/2))
     .         +2d0*l*muq*vuq*DDCOS(phiS)*DDSIN(phiH/2d0)
      MH2(3,1)=MH2(1,3)
      MH2(2,2)=MHdS+muq**2+(l**2-(g1q+g2q)/4d0)*vuq**2
     .        +(g1q+g2q)/4d0*vdq**2*(1d0+DDCOS(phiH/2d0)**2)
      MH2(2,3)=-l*vuq*(Alcos1*DDCOS(phiH/2d0)-IAl*DDSIN(phiH/2d0)
     .                +2d0*k/l*muq*DDCOS(phi0-phiS+phiH/2d0)
     .                +MUPQ*DDCOS(phi01-phiP+phiH/2d0))
     .         +2d0*l*muq*vdq*DDCOS(phiS)*DDSIN(phiH/2d0)
      MH2(3,2)=MH2(2,3)
      MH2(3,3)=MSS+l**2*(vuq**2+vdq**2)+MUPQ**2
     .        +2d0*k/l*muq*(Akcos2*DDCOS(phiS)-IAk*DDSIN(phiS))
     .        +2d0*(k/l*muq)**2*(1d0+1d0*DDCOS(phiS)**2)
     .        -2d0*l*k*vuq*vdq*DDCOS(phi0+phiH)
     .        +MSPQ*DDCOS(phiSPq)+2d0*k*XIFQ*DDCOS(phi02-phiF)
     .        +2d0*k/l*muq*MUPQ*(DDCOS(phi02-phiP+phiS)
     .                          +2d0*DDCOS(phi02-phiP)*DDCOS(phiS))
      MH2(1,4)=(g1q+g2q)/4d0*vuq**2*DDSIN(phiH)
      MH2(4,1)=MH2(1,4)
      MH2(2,4)=muq*(IAl*DDCOS(phiS)+Alcos1*DDSIN(phiS)
     .              +k/l*muq*DDSIN(phi0-2d0*phiS))
     .        +(l**2-(g1q+g2q)/4d0)*vuq*vdq*DDSIN(phiH)
     .        +M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-phiF)
     .        +muq*MUPQ*DDSIN(phi01-phiP-phiS)
      MH2(4,2)=MH2(2,4)
      MH2(3,4)=l*vdq*(IAl*DDCOS(phiH/2d0)+Alcos1*DDSIN(phiH/2d0)
     .               +2d0*k/l*muq*DDSIN(phi0-phiS+phiH/2d0)
     .               +MUPQ*DDSIN(phi01-phiP+phiH/2d0))
     .        +2d0*l*muq*DDCOS(phiS)*DDSIN(phiH/2d0)
      MH2(4,3)=MH2(3,4)
      MH2(4,4)=MHuS+muq**2+(l**2-(g1q+g2q)/4d0)*vdq**2
     .        +(g1q+g2q)/4d0*vuq**2*(1d0+2d0*DDSIN(phiH/2d0)**2)
      MH2(1,5)=muq*(IAl*DDCOS(phiS)+Alcos1*DDSIN(phiS)
     .             +k/l*muq*DDSIN(phi0-2d0*phiS))
     .       +(l**2-(g1q+g2q)/4d0)*vuq*vdq*DDSIN(phiH)
     .       +M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-phiF)
     .       +muq*MUPQ*DDSIN(phi01-phiP-phiS)
      MH2(5,1)=MH2(1,5)
      MH2(2,5)=(g1q+g2q)/4d0*vdq**2*DDSIN(phiH)
      MH2(5,2)=MH2(2,5)
      MH2(3,5)=l*vuq*(IAl*DDCOS(phiH/2d0)+Alcos1*DDSIN(phiH/2d0)
     .               +2d0*k/l*muq*DDSIN(phi0-phiS+phiH/2d0)
     .               +MUPQ*DDSIN(phi01-phiP+phiH/2d0))
     .         +2d0*l*muq*vdq*DDCOS(phiS)*DDSIN(phiH/2d0)
      MH2(5,3)=MH2(3,5)
      MH2(4,5)=muq*(Alcos1*DDCOS(phiS)-IAl*DDSIN(phiS)
     .             +k/l*muq*DDCOS(phi0-2d0*phiS))
     .        +2d0*(l**2-(g1q+g2q)/4d0)*vuq*vdq*DDSIN(phiH/2d0)**2
     .        +M3HQ*DDCOS(phi3q)+l*XIFQ*DDCOS(phi01-phiF)
     .        +muq*MUPQ*DDCOS(phi01-phiP-phiS)
      MH2(5,4)=MH2(4,5)
      MH2(5,5)=MHdS+muq**2+(l**2-(g1q+g2q)/4d0)*vuq**2
     .        +(g1q+g2q)/4d0*vdq**2*(1d0+2d0*DDSIN(phiH/2d0)**2)
      MH2(1,6)=l*vdq*(IAl*DDCOS(phiH/2d0)+Alcos1*DDSIN(phiH/2d0)
     .                -2d0*k/l*muq*DDSIN(phi0-phiS+phiH/2d0)
     .                +MUPQ*DDSIN(phi01-phiP+phiH/2d0))
     .         +2d0*l*muq*vuq*DDSIN(phiS)*DDCOS(phiH/2d0)
      MH2(6,1)=MH2(1,6)
      MH2(2,6)=l*vuq*(IAl*DDCOS(phiH/2d0)+Alcos1*DDSIN(phiH/2d0)
     .               -2d0*k/l*muq*DDSIN(phi0-phiS+phiH/2d0)
     .               +MUPQ*DDSIN(phi01-phiP+phiH/2d0))
     .         +2d0*l*muq*vdq*DDSIN(phiS)*DDCOS(phiH/2d0)
      MH2(6,2)=MH2(2,6)
      MH2(3,6)=-2d0*k/l*muq*(IAk*DDCOS(phiS)+Akcos2*DDSIN(phiS)
     .                      -k/l*muq*DDSIN(2d0*phiS))
     .   -2d0*k*l*vuq*vdq*DDSIN(phi0+phiH)
     .   -MSPQ*DDSIN(phiSPq)-2d0*k*XIFQ*DDSIN(phi02-phiF)
     .   +2d0*k/l*muq*MUPQ*(-DDSIN(phi02-phiP)*DDCOS(phiS)
     .                        +DDCOS(phi02-phiP)*DDSIN(phiS))
      MH2(6,3)=MH2(3,6)
      MH2(4,6)=l*vdq*(Alcos1*DDCOS(phiH/2d0)-IAl*DDSIN(phiH/2d0)
     .               -2d0*k/l*muq*DDCOS(phi0-phiS+phiH/2d0)
     .               +MUPQ*DDCOS(phi01-phiP+phiH/2d0))
     .        +2d0*l*muq*vuq*DDSIN(phiS)*DDCOS(phiH/2d0)
      MH2(6,4)=MH2(4,6)
      MH2(5,6)=l*vuq*(Alcos1*DDCOS(phiH/2d0)-IAl*DDSIN(phiH/2d0)
     .               -2d0*k/l*muq*DDCOS(phi0-phiS+phiH/2d0)
     .               +MUPQ*DDCOS(phi01-phiP+phiH/2d0))
     .        +2d0*l*muq*vdq*DDSIN(phiS)*DDCOS(phiH/2d0)
      MH2(6,5)=MH2(4,5)
      MH2(6,6)=MSS-2d0*k/l*muq*(Akcos2*DDCOS(phiS)-IAk*DDSIN(phiS)
     .                         +k/l*muq*(1d0+2d0*DDSIN(phiS)**2))
     .      +l**2*(vuq**2+vdq**2)+2d0*l*k*vuq*vdq*DDCOS(phi0+phiH)
     .      +MUPQ**2-MSPQ*DDCOS(phiSPq)-2d0*k*XIFQ*DDCOS(phi02-phiF)
     .      +2d0*k/l*muq*MUPQ*(DDCOS(phi02-phiP+phiS)
     .                        -2d0*DDSIN(phi02-phiP)*DDSIN(phiS))
      CALL DIAGN(6,MH2,VALPH,VECPH,1d-10)

      Vphi=MHuS*vuq**2+MHdS*vdq**2+MSS*(muq/l)**2
     . -2d0*muq*vuq*vdq*(Alcos1*DDCOS(phiS+phiH)-IAl*DDSIN(phiS+phiH))
     . +2d0/3d0*k*(muq/l)**3*(Akcos2*DDCOS(3*PhiS)-IAk*DDSIN(3d0*phiS))
     . +muq**2*(vuq**2+vdq**2)+(l*vuq*vdq)**2
     . -2d0*k/l*muq**2*vuq*vdq*DDCOS(phi0+phiH-2d0*phiS)
     . +k**2*(muq/l)**4+(G1Q+G2Q)/8d0*(vuq**2-vdq**2)**2
     . -2d0*(M3HQ*DDCOS(phi3q+phiH)+l*XIFQ*DDCOS(phi01-phiF+phiH)
     .                 +MUPQ*muq*DDCOS(phi01-phiP+phiH-phiS))*vuq*vdq
     . +2d0*(XISQ*DDCOS(phiSq)*DDCOS(phiS)-IXIS*DDSIN(PhiS)
     .                   +XIFQ*MUPQ*DDCOS(phiP-phiF+phiS))*muq/l
     . +2d0*(MSPQ/2d0*DDCOS(phiSPq+2d0*phiS)
     .                    +k*XIFQ*DDCOS(phi02-phiF+2d0*phiS))*(muq/l)**2
     . +2d0*k*MUPQ*DDCOS(phi02-phiP+phiS)*(muq/l)**3
     . +(MUPQ*muq/l)**2
     . -(3d0*mtopq**4*(dlog(mtopq**2/QSTSB)-1.5d0)
     .  +3d0*mbotq**4*(dlog(mbotq**2/QSTSB)-1.5d0)
     .  +(mtau*vdq/vd)**4*(dlog(mtau**2/QSTSB)-1.5d0))/16d0/Pi**2
     . +3d0*(2d0*MW2**2*(dlog(MW2/QSTSB)-1.5d0)
     .             +MZ2*(dlog(MZ2/QSTSB)-1.5d0))/64d0/Pi**2
     . +(3d0*MSTaux1**2*(dlog(Max(MSTaux1,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSTaux2**2*(dlog(Max(MSTaux2,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSBaux1**2*(dlog(Max(MSBaux1,MZ2)/QSTSB)-1.5d0)
     .     +3d0*MSBaux2**2*(dlog(Max(MSBaux2,MZ2)/QSTSB)-1.5d0)
     .     +MSLaux1**2*(dlog(Max(MSLaux1,MZ2)/QSTSB)-1.5d0)
     .     +MSLaux2**2*(dlog(Max(MSLaux2,MZ2)/QSTSB)-1.5d0)
     .     +MSNTaux**2*(dlog(Max(MSNTaux,MZ2)/QSTSB)-1.5d0)
     .     +6d0*MSU2(1)**2*(dlog(MSU2(1)/QSTSB)-1.5d0)
     .     +6d0*MSU2(2)**2*(dlog(MSU2(2)/QSTSB)-1.5d0)
     .     +6d0*MSD2(1)**2*(dlog(MSD2(1)/QSTSB)-1.5d0)
     .     +6d0*MSD2(2)**2*(dlog(MSD2(2)/QSTSB)-1.5d0)
     .     +2d0*MSE2(1)**2*(dlog(MSE2(1)/QSTSB)-1.5d0)
     .     +2d0*MSE2(2)**2*(dlog(MSE2(2)/QSTSB)-1.5d0)
     .     +2d0*MSNE2**2*(dlog(MSNE2/QSTSB)-1.5d0))/32d0/Pi**2
     .    +3d0/512d0/Pi**4*Ytq**4*((dlog(QSTSB/mtopq**2))**2
     .      *(64d0*Pi*ALSQ+4d0/3d0*g1q-3d0*sinb**2*Ytq**2
     .       +3d0*cosb**2*Ybq**2)+
     .       ((dlog(MA2/mtopq**2))**2-(dlog(QSTSB/mtopq**2))**2)
     .       *(3d0*cosb**2*Ytq**2+(3d0*cosb**2+1d0)*Ybq**2))*vuq**4
     .    +3d0/512d0/Pi**4*Ybq**4*(dlog(QSTSB/mtopq**2)**2
     .      *(64d0*Pi*ALSQ-2d0/3d0*g1q+3d0*sinb**2*Ytq**2
     .       -3d0*cosb**2*Ybq**2)+
     .       (dlog(MA2/mtopq**2)**2-dlog(QSTSB/mtopq**2)**2)
     .       *(3d0*sinb**2*Ybq**2+(3d0*sinb**2+1d0)*Ytq**2))*vdq**4
     .    -(MCH1aux**2*(dlog(Max(MCH1aux,MZ2)/QSTSB)-1.5d0)
     .     +MCH2aux**2*(dlog(Max(MCH2aux,MZ2)/QSTSB)-1.5d0))/16d0/Pi**2
     .    -(VALPN(1)**2*(dlog(Max(VALPN(1),1d0)/QSTSB)-1.5d0)
     .     +VALPN(3)**2*(dlog(Max(VALPN(3),1d0)/QSTSB)-1.5d0)
     .     +VALPN(5)**2*(dlog(Max(VALPN(5),1d0)/QSTSB)-1.5d0)
     .     +VALPN(7)**2*(dlog(Max(VALPN(7),1d0)/QSTSB)-1.5d0)
     .   +VALPN(9)**2*(dlog(Max(VALPN(9),1d0)/QSTSB)-1.5d0))/32d0/Pi**2
     .     +(MHD1aux**2*(dlog(Max(MHD1aux,MW2)/QSTSB)-1.5d0)
     .     +MHD2aux**2*(dlog(Max(MHD2aux,MW2)/QSTSB)-1.5d0))/32d0/Pi**2
     .    +(VALPH(1)**2*(dlog(Max(MZ2,VALPH(1))/QSTSB)-1.5d0)
     .     +VALPH(2)**2*(dlog(Max(MZ2,VALPH(2))/QSTSB)-1.5d0)
     .     +VALPH(3)**2*(dlog(Max(MZ2,VALPH(3))/QSTSB)-1.5d0)
     .     +VALPH(4)**2*(dlog(Max(MZ2,VALPH(4))/QSTSB)-1.5d0)
     .     +VALPH(5)**2*(dlog(Max(MZ2,VALPH(5))/QSTSB)-1.5d0)
     .   +VALPH(6)**2*(dlog(Max(MZ2,VALPH(6))/QSTSB)-1.5d0))/64d0/Pi**2

      ENDDO
      ENDDO

      IF(Vphi.lt.VEW)PROB(28)=Max(PROB(28),
     .                     DDIM(-1d-2,(Vphi-VEW)/Max(dabs(VEW),1d-10)))


c      Naturalness of the doublet squared masses

      PROB(29)=DDIM(MAX(DABS(MHuS),DABS(MHdS))/Q2,10d0)

      RETURN
      END
