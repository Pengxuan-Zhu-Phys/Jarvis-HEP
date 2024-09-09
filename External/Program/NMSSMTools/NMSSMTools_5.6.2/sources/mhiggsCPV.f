      SUBROUTINE MHIGGSTREE_CPV()

c       Tree-level Higgs masses.
c      - The tree-level neutral Higgs squared-mass matrix MH02(i,j), (i,j=1..5;
c        1 -> hu, 2 -> hd, 3 -> hs, 4 -> cosb*au+sinb*ad, 5 -> as) and the
c        tree-level charged-higgs squared mass MHC2 are stored in the common
c        SQUHIMASSM, for further processing at the squark squared-scale QSTSB.
c      - The tree-level (at the weak-scale) squared-masses for the charged
c        Higgs, MHC, and neutral Higgses, MH0(i), (i=1..5), are stored 
c        (temporarily) in the common HISPEC, together with the rotation
c        matrices XC(i,j), (i,j=1..2) for the charged-Higgs: beta angle,
c        and XH(i,j), (i,j=1..5) for the neutral Higgses, as well as the 
c        heavy-doublet squared-mass MA2. These tree-level quantities shall
c        be used for the implementation of the Higgs-loop contributions to
c        the Higgs masses (mhiggsloop_gaugehiggs_CPV.f and
c        mhiggsloop_pole_CPV.f), before being replaced by the loop-corrected
c        quantities.
c      - The tree-level contributions to the effective Z3-conserving parameters
c        of the Higgs potential are stored in the common EFFPOTPAR.

      IMPLICIT NONE

      INTEGER I,J

      DOUBLE PRECISION VALPH(5),VECPH(5,5),MHT2(5,5)
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION l,k,Alcos1,Akcos2,muq,nuq
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION MH02(5,5),MHC2
      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION lu,ld,l3,l4,Rel5,Iml5,Rel6,Iml6,Rel7,Iml7,RAud,
     . RAS,K2,lPu,lPd,RlPM,IlPM
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      DOUBLE PRECISION phIF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      DOUBLE PRECISION Rm3,Im3,RAudt,IAudt,RxS,IxS,Rmsp,Imsp,RAst,IAst
      DOUBLE PRECISION RlPMt,IlPMt,RlM,IlM,RAqs,IAqs,Rltqs,Iltqs
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/QPAR/l,k,Alcos1,Akcos2,muq,NUQ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/SQUHIMASSM/MH02,MHC2
      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/EFFPOTPAR/lu,ld,l3,l4,Rel5,Iml5,Rel6,Iml6,Rel7,Iml7,RAud,
     . RAS,K2,lPu,lPd,RlPM,IlPM
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/QEXT/XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      COMMON/Z3VAUX/phIF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      COMMON/Z3VPOT/Rm3,Im3,RAudt,IAudt,RxS,IxS,Rmsp,Imsp,RAst,IAst,
     . RlPMt,IlPMt,RlM,IlM,RAqs,IAqs,Rltqs,Iltqs

c         A: Charged Higgs

c      I- Squared mass

      MHC2=((muq/vuq/vdq*(Alcos1+k/l*muq*DDCOS(phi0))+g2q/2d0-l**2)
     .       +(m3Hq*DDCOS(phi3q)+l*xIFQ*DDCOS(phi01-phIF)
     .        +muq*mupq*DDCOS(phi01-phip))/vuq/vdq)*(vu**2+vd**2)          !/(ZHu*ZHd)

c      II- Tree level input for radiative corrections

      MHC=((muq/vuq/vdq*(Alcos1+k/l*muq*DDCOS(phi0))+g2q/2d0-l**2)
     .       +(m3Hq*DDCOS(phi3q)+l*xIFQ*DDCOS(phi01-phIF)
     .        +muq*mupq*DDCOS(phi01-phip))/vuq/vdq)*(vuq**2+vdq**2)
      XC(1,1)=-vuq/dsqrt(vuq**2+vdq**2)
      XC(1,2)=vdq/dsqrt(vuq**2+vdq**2)
      XC(2,1)=vdq/dsqrt(vuq**2+vdq**2)
      XC(2,2)=vuq/dsqrt(vuq**2+vdq**2)
      IF(MHC.le.0d0)then
c       print*,higgspro
       MHC=1d0
      ENDIF


c         B: Neutral states: diag(MH0(I)^2) = XH.MH02.XHt

c      I- Squared mass-matrix

      MH02(1,1)=muq*(Alcos1+k/l*muq*DDCOS(phi0))*vdq/vuq
     .             +(g1q+g2q)/2d0*vuq**2
     .            +(m3Hq*DDCOS(phi3q)+l*xIFQ*DDCOS(phi01-phIF)
     .                +muq*mupq*DDCOS(phi01-phip))*vdq/vuq
      MH02(1,2)=-muq*(Alcos1+k/l*muq*DDCOS(phi0))
     .             +(2d0*l**2-(g1q+g2q)/2d0)*vuq*vdq
     .            -(m3Hq*DDCOS(phi3q)+l*xIFQ*DDCOS(phi01-phIF)
     .                +muq*mupq*DDCOS(phi01-phip))
      MH02(2,1)=MH02(1,2)
      MH02(1,3)=-l*vdq*(Alcos1+2d0*k/l*muq*DDCOS(phi0)
     .                    +mupq*DDCOS(Phi01-phip))
     .            +2d0*l*muq*vuq
      MH02(3,1)=MH02(1,3)
      MH02(1,4)=0d0
      MH02(4,1)=MH02(1,4)
      MH02(1,5)=-3.d0*k*muq*vdq*DDSIN(phi0)
     .    -l*vdq/muq*(m3Hq*DDSIN(phi3q)+l*xIFq*DDSIN(phi01-phIF)
     .                    +2d0*muq*mupq*DDSIN(phi01-phip))
      MH02(5,1)=MH02(1,5)
      MH02(2,2)=muq*(Alcos1+k/l*muq*DDCOS(phi0))*vuq/vdq
     .             +(g1q+g2q)/2d0*vdq**2
     .            +(m3Hq*DDCOS(phi3q)+l*xIFQ*DDCOS(phi01-phIF)
     .             +muq*mupq*DDCOS(phi01-phip))*vuq/vdq
      MH02(2,3)=-l*vuq*(Alcos1+2d0*k/l*muq*DDCOS(phi0)
     .                    +mupq*DDCOS(Phi01-phip))
     .            +2d0*l*muq*vdq
      MH02(3,2)=MH02(2,3)
      MH02(2,4)=0d0
      MH02(4,2)=MH02(2,4)
      MH02(2,5)=-3.d0*k*muq*vuq*DDSIN(phi0)
     .    -l*vuq/muq*(m3Hq*DDSIN(phi3q)+l*xIFq*DDSIN(phi01-phIF)
     .                    +2d0*muq*mupq*DDSIN(phi01-phip))
      MH02(5,2)=MH02(2,5)
      MH02(3,3)=l**2*vuq*vdq/muq*(Alcos1+mupq*DDCOS(phi01-phip))
     .   +k/l*muq*(Akcos2+4.d0*k/l*muq+3d0*mupq*DDCOS(phi02-phip))
     .   -l/muq*(xiSq*DDCOS(phiSq)+xIFq*mupq*DDCOS(phip-phIF))
      MH02(3,4)=dsqrt(vu**2+vd**2)*(muq*k*DDSIN(phi0)
     .    -l/muq*(m3Hq*DDSIN(phi3q)+l*xIFq*DDSIN(phi01-phIF)))  ! /dsqrt(Zs*ZHu*ZHd)
      MH02(4,3)=MH02(3,4)
       IF(k.ne.0d0)then
      MH02(3,5)=4*l*k*vuq*vdq*DDSIN(phi0)
     .   +2d0*l/muq*(xiSq*DDSIN(phiSq)+xIFq*mupq*DDSIN(phip-phIF))
     .   +MSPq*DDSIN(phiSPq)+2d0*k*xIFq*DDSIN(phi02-phIF)
     .   +2d0*l**2*vuq*vdq/muq**2*(m3hq*DDSIN(phi3q)
     .      +l*xIFq*DDSIN(phi01-phIF)+2d0*muq*mupq*DDSIN(phi01-phip))
       ELSE
      MH02(3,5)=-MSPq*DDSIN(phiSPq)
       ENDIF
      MH02(5,3)=MH02(3,5)
      MH02(4,4)=(muq*(Alcos1+k/l*muq*DDCOS(phi0))
     .               +m3Hq*DDCOS(phi3q)+l*xIFQ*DDCOS(phi01-phIF)
     .                 +muq*mupq*DDCOS(phi01-phip))
     .                 /vuq/vdq*(vu**2+vd**2)       ! /(ZHu*ZHd)
      MH02(4,5)=l*dsqrt(vu**2+vd**2)*(Alcos1
     .       -2d0*k/l*muq*DDCOS(phi0)-mupq*DDCOS(Phi01-phip)) ! /dsqrt(Zs*ZHu*ZHd)
      MH02(5,4)=MH02(4,5)
      MH02(5,5)=-k/l*muq*(3.d0*Akcos2+mupq*DDCOS(phi02-phip))
     .   +l**2*vuq*vdq/muq*(Alcos1+4.d0*k/l*muq*DDCOS(phi0)
     .                        +mupq*DDCOS(phi01-phip))
     .   -l/muq*(xiSq*DDCOS(phiSq)+xIFq*mupq*DDCOS(phip-phIF))
     .    -2d0*(MSPq*DDCOS(phiSPq)+2d0*k*xIFq*DDCOS(phi02-phIF))

c      print*,'MH02_1*',MH02(1,1),MH02(1,2),MH02(1,3),MH02(1,4),MH02(1,5)
c      print*,'MH02_2*',MH02(2,1),MH02(2,2),MH02(2,3),MH02(2,4),MH02(2,5)
c      print*,'MH02_3*',MH02(3,1),MH02(3,2),MH02(3,3),MH02(3,4),MH02(3,5)
c      print*,'MH02_4*',MH02(4,1),MH02(4,2),MH02(4,3),MH02(4,4),MH02(4,5)
c      print*,'MH02_5*',MH02(5,1),MH02(5,2),MH02(5,3),MH02(5,4),MH02(5,5)

cUE:
c      WRITE(0,*)"MH02 =",
c     .   MH02(1,1),MH02(2,2),MH02(3,3),MH02(1,2),MH02(1,3),MH02(2,3)


c      II- Tree level input for radiative corrections

      MA2=(muq*(Alcos1+k/l*muq*DDCOS(phi0))
     .            +m3Hq*DDCOS(phi3q)+l*xIFQ*DDCOS(phi01-phIF)
     .                +muq*mupq*DDCOS(phi01-phip))
     .          /vuq/vdq*(vu**2+vd**2)/(ZHu*ZHd)

      MHT2(1,1)=muq*(Alcos1+k/l*muq*DDCOS(phi0))*vdq/vuq
     .             +(g1q+g2q)/2d0*vuq**2
     .            +(m3Hq*DDCOS(phi3q)+l*xIFQ*DDCOS(phi01-phIF)
     .                +muq*mupq*DDCOS(phi01-phip))*vdq/vuq
      MHT2(1,2)=-muq*(Alcos1+k/l*muq*DDCOS(phi0))
     .             +(2d0*l**2-(g1q+g2q)/2d0)*vuq*vdq
     .            -(m3Hq*DDCOS(phi3q)+l*xIFQ*DDCOS(phi01-phIF)
     .                +muq*mupq*DDCOS(phi01-phip))
      MHT2(2,1)=MHT2(1,2)
      MHT2(1,3)=-l*vdq*(Alcos1+2d0*k/l*muq*DDCOS(phi0)
     .                    +mupq*DDCOS(Phi01-phip))
     .            +2d0*l*muq*vuq
      MHT2(3,1)=MHT2(1,3)
      MHT2(1,4)=0d0
      MHT2(4,1)=MHT2(1,4)
      MHT2(1,5)=-3.d0*k*muq*vdq*DDSIN(phi0)
     .    -l*vdq/muq*(m3Hq*DDSIN(phi3q)+l*xIFq*DDSIN(phi01-phIF)
     .                    +2d0*muq*mupq*DDSIN(phi01-phip))
      MHT2(5,1)=MHT2(1,5)
      MHT2(2,2)=muq*(Alcos1+k/l*muq*DDCOS(phi0))*vuq/vdq
     .             +(g1q+g2q)/2d0*vdq**2
     .            +(m3Hq*DDCOS(phi3q)+l*xIFQ*DDCOS(phi01-phIF)
     .             +muq*mupq*DDCOS(phi01-phip))*vuq/vdq
      MHT2(2,3)=-l*vuq*(Alcos1+2d0*k/l*muq*DDCOS(phi0)
     .                    +mupq*DDCOS(Phi01-phip))
     .            +2d0*l*muq*vdq
      MHT2(3,2)=MHT2(2,3)
      MHT2(2,4)=0d0
      MHT2(4,2)=MHT2(2,4)
      MHT2(2,5)=-3.d0*k*muq*vuq*DDSIN(phi0)
     .    -l*vuq/muq*(m3Hq*DDSIN(phi3q)+l*xIFq*DDSIN(phi01-phIF)
     .                    +2d0*muq*mupq*DDSIN(phi01-phip))
      MHT2(5,2)=MHT2(2,5)
      MHT2(3,3)=l**2*vuq*vdq/muq*(Alcos1+mupq*DDCOS(phi01-phip))
     .   +k/l*muq*(Akcos2+4.d0*k/l*muq+3d0*mupq*DDCOS(phi02-phip))
     .   -l/muq*(xiSq*DDCOS(phiSq)+xIFq*mupq*DDCOS(phip-phIF))
      MHT2(3,4)=dsqrt(vuq**2+vdq**2)*(muq*k*DDSIN(phi0)
     .    -l/muq*(m3Hq*DDSIN(phi3q)+l*xIFq*DDSIN(phi01-phIF)))
      MHT2(4,3)=MHT2(3,4)
       IF(k.ne.0d0)then
      MHT2(3,5)=4*l*k*vuq*vdq*DDSIN(phi0)
     .   +2d0*l/muq*(xiSq*DDSIN(phiSq)+xIFq*mupq*DDSIN(phip-phIF))
     .   +MSPq*DDSIN(phiSPq)+2d0*k*xIFq*DDSIN(phi02-phIF)
     .   +2d0*l**2*vuq*vdq/muq**2*(m3hq*DDSIN(phi3q)
     .      +l*xIFq*DDSIN(phi01-phIF)+2d0*muq*mupq*DDSIN(phi01-phip))
       ELSE
      MHT2(3,5)=-MSPq*DDSIN(phiSPq)
       ENDIF
      MHT2(5,3)=MHT2(3,5)
      MHT2(4,4)=(muq*(Alcos1+k/l*muq*DDCOS(phi0))
     .               +m3Hq*DDCOS(phi3q)+l*xIFQ*DDCOS(phi01-phIF)
     .                 +muq*mupq*DDCOS(phi01-phip))
     .                 /vuq/vdq*(vuq**2+vdq**2)
      MHT2(4,5)=l*(Alcos1-2d0*k/l*muq*DDCOS(phi0)
     .                -mupq*DDCOS(Phi01-phip))*dsqrt(vuq**2+vdq**2)
      MHT2(5,4)=MHT2(4,5)
      MHT2(5,5)=-k/l*muq*(3.d0*Akcos2+mupq*DDCOS(phi02-phip))
     .   +l**2*vuq*vdq/muq*(Alcos1+4.d0*k/l*muq*DDCOS(phi0)
     .                        +mupq*DDCOS(phi01-phip))
     .   -l/muq*(xiSq*DDCOS(phiSq)+xIFq*mupq*DDCOS(phip-phIF))
     .    -2d0*(MSPq*DDCOS(phiSPq)+2d0*k*xIFq*DDCOS(phi02-phIF))

      CALL DIAGN(5,MHT2,VALPH,VECPH,1.d-10)
      CALL SORTNA(5,VALPH,VECPH)

cUE:
c      write(*,*) "In mhiggstree_CPV: MH02 off-diag:",
c     .   MH02(1,4),MH02(2,4),MH02(3,4)
c     .   ,MH02(1,5),MH02(2,5),MH02(3,5)
c      write(*,*) "In mhiggstree_CPV: MH02 4-5:",
c     .   MH02(4,4),MH02(5,5),MH02(4,5)
c      write(*,*) "In mhiggstree_CPV: VALPH:",VALPH

      DO I=1,5
      IF(VALPH(I).le.0d0)then
       VALPH(I)=1d0
      ENDIF
      MH0(I)=VALPH(I)
       DO J=1,5
      XH(I,J)=VECPH(J,I)
       ENDDO
      ENDDO

c        C: Parameters of the Effective Potential
      lu=(g1q+g2q)/4.d0
      ld=(g1q+g2q)/4.d0
      l3=(-g1q+g2q)/4.d0
      l4=l**2-g2q/2d0
      Rel5=0d0
      Iml5=0d0
      Rel6=0d0
      Iml6=0d0
      Rel7=0d0
      Iml7=0d0
      RAud=l*Alcos1
      RAS=k*Akcos2
      K2=k**2
      lPu=l**2
      lPd=l**2
      RlPM=k*l*DDCOS(Phi0)
      IlPM=k*l*DDSIN(Phi0)

      Rm3=m3hq*DDCOS(phi3q)+l*xIFq*DDCOS(phi01-phIF)
      Im3=m3hq*DDSIN(phi3q)+l*xIFq*DDSIN(phi01-phIF)
      RAudt=l*mupq*DDCOS(phi01-phip)
      IAudt=l*mupq*DDSIN(phi01-phip)
      RxS=xiSq*DDCOS(phiSq)+mupq*xIFq*DDCOS(phip-phIF)
      IxS=xiSq*DDSIN(phiSq)+mupq*xIFq*DDSIN(phip-phIF)
      Rmsp=mspq*DDCOS(phiSPq)+2d0*k*xIFq*DDCOS(phi02-phIF)
      Imsp=mspq*DDSIN(phiSPq)+2d0*k*xIFq*DDSIN(phi02-phIF)
      RAst=k*mupq*DDCOS(phi02-phip)
      IAst=k*mupq*DDSIN(phi02-phip)
      RlPMt=0d0
      IlPMt=0d0
      RlM=0d0
      IlM=0d0
      RAqs=0d0
      IAqs=0d0
      Rltqs=0d0
      Iltqs=0d0

c      print*,'mHc',dsqrt(mHc2)
c      print*,'mH01',MH0(1)
c      print*,'mH02',MH0(2)
c      print*,'mH03',MH0(3)
c      print*,'mH04',MH0(4)
c      print*,'mH05',MH0(5)
c      print*,'XH1*',XH(1,1),XH(1,2),XH(1,3),XH(1,4),XH(1,5)
c      print*,'XH2*',XH(2,1),XH(2,2),XH(2,3),XH(2,4),XH(2,5)
c      print*,'XH3*',XH(3,1),XH(3,2),XH(3,3),XH(3,4),XH(3,5)
c      print*,'XH4*',XH(4,1),XH(4,2),XH(4,3),XH(4,4),XH(4,5)
c      print*,'XH5*',XH(5,1),XH(5,2),XH(5,3),XH(5,4),XH(5,5)

      RETURN
      END

************************************************************************************************

      SUBROUTINE MHIGGSLOOP_SFERM_CPV()

c         One-loop corrections to the Higgs potential + leading 2-loop
c                 - SM-fermions + Sfermions contribution
c      - The 1-loop + leading 2-loop corrections from SM-fermions and sfermions
c        to the parameters of the effective Higgs potential are added and stored
c        in the common EFFPOTPAR.
c      - The corresponding corrections to the squared-mass matrices of the 
c        neutral-Higgs states, as well as the charged one, are added to MH0(i,j)
c        or MHC2 and stored within the common SQUHIMASSM.

      IMPLICIT NONE

      INTEGER I,J

      DOUBLE PRECISION Pi,aux,Ytau,fsf1,fsf2,fsf3,fsf5,fsf6,fsf7,s
      DOUBLE PRECISION dTdh(5),dRdh(5),R2,Rhh(5,5)
      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION Ytq,Ybq,MTOPQ,MBOTQ
      DOUBLE PRECISION mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION MSQ3,MSU3,MSD3,AT,AB
      DOUBLE PRECISION MSQ1,MSU1,MSD1
      DOUBLE PRECISION l,k,Alcos1,Akcos2,muq,nuq
      DOUBLE PRECISION MST2(2),UT(2,2,2),MSB2(2),UB(2,2,2),MSL2(2),
     . UTAU(2,2,2),MSNT2
      DOUBLE PRECISION MSU2(2),MSD2(2),MSE2(2),MSNE2,MSMU2(2),
     . UMU(2,2,2)
      DOUBLE PRECISION MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      DOUBLE PRECISION MH02(5,5),MHC2
      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION lu,ld,l3,l4,Rel5,Iml5,Rel6,Iml6,Rel7,Iml7,RAud,
     . RAS,K2,lPu,lPd,RlPM,IlPM
      DOUBLE PRECISION Rm3,Im3,RAudt,IAudt,RxS,IxS,Rmsp,Imsp,RAst,IAst,
     . RlPMt,IlPMt,RlM,IlM,RAqs,IAqs,Rltqs,Iltqs
      DOUBLE PRECISION MHuS,MHdS,MSS
      DOUBLE PRECISION IAL,IAK,IXIS
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/STSBSCALE/QSTSB
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/QQUARK/Ytq,Ybq,MTOPQ,MBOTQ
      COMMON/SMFERM/mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/RADCOR2/MSQ3,MSU3,MSD3,AT,AB
      COMMON/SQUPAR/MSQ1,MSU1,MSD1
      COMMON/QPAR/l,k,Alcos1,Akcos2,muq,NUQ
      COMMON/SFERM3SPEC/MST2,UT,MSB2,UB,MSL2,UTAU,MSNT2
      COMMON/SFERM1SPEC/MSU2,MSD2,MSE2,MSNE2,MSMU2,UMU
      COMMON/SLEPPAR/MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      COMMON/SQUHIMASSM/MH02,MHC2
      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/EFFPOTPAR/lu,ld,l3,l4,Rel5,Iml5,Rel6,Iml6,Rel7,Iml7,RAud,
     . RAS,K2,lPu,lPd,RlPM,IlPM
      COMMON/Z3VPOT/Rm3,Im3,RAudt,IAudt,RxS,IxS,Rmsp,Imsp,RAst,IAst,
     . RlPMt,IlPMt,RlM,IlM,RAqs,IAqs,Rltqs,Iltqs
      COMMON/MH2TREE/MHuS,MHdS,MSS
      COMMON/IMALAK/IAL,IAK,IXIS

      PI=4d0*DATAN(1d0)

c            A: Corrections to the neutral mass-matrix / 3rd generation

c        I - Tops/Stops

c      a) Tops

      aux=4.d0*Ytq**4*vuq**2*dlog(mtopq**2/QSTSB)
      MH02(1,1)=MH02(1,1)-3.d0/16.d0/Pi**2*aux

c      b) Stops

      DO I=1,5
       dTdh(I)=0d0
       dRdh(I)=0d0
       DO J=1,5
        Rhh(I,J)=0d0
       ENDDO
      ENDDO

      R2=(MSQ3-MSU3+(5.d0/3.d0*g1q-g2q)/4.d0*
     .  (vuq**2-vdq**2))**2+4.d0*Ytq**2*(At**2*vuq**2
     .     +muq**2*vdq**2-2d0*muq*At*vuq*vdq*DDCOS(PhiAt+Phi01))

      dTdh(1)=4.d0*(Ytq**2-(g1q+g2q)/8.d0)*vuq
      dTdh(2)=(g1q+g2q)/2d0*vdq

      dRdh(1)=8.d0*Ytq**2*At*(vuq*At-muq*vdq*DDCOS(PhiAt+Phi01))
     .          +(5.d0/3.d0*g1q-g2q)*vuq*
     .   (MSQ3-MSU3+(5.d0/3.d0*g1q-g2q)/4.d0*(vuq**2-vdq**2))
      dRdh(2)=8.d0*Ytq**2*muq*(vdq*muq-At*vuq*DDCOS(PhiAt+Phi01))
     .          -(5.d0/3.d0*g1q-g2q)*vdq*
     .   (MSQ3-MSU3+(5.d0/3.d0*g1q-g2q)/4.d0*(vuq**2-vdq**2))
      dRdh(3)=-8.d0*Ytq**2*l*vdq*(At*vuq*DDCOS(PhiAt+Phi01)-muq*vdq)
      dRdh(4)=8.d0*Ytq**2*muq*dsqrt(vu**2+vd**2)*At*DDSIN(PhiAt+Phi01)
      dRdh(5)=8.d0*Ytq**2*l*vuq*vdq*At*DDSIN(PhiAt+Phi01)

      Rhh(1,1)=4.d0*(((5.d0/3.d0*g1q-g2q)/4.d0*vuq)**2
     .     +Ytq**2*muq*At*vdq/vuq*DDCOS(PhiAt+Phi01))
      Rhh(2,2)=4.d0*((5.d0/3.d0*g1q-g2q)/4.d0*vdq)**2
     .     +4.d0*Ytq**2*muq*At*vuq/vdq*DDCOS(PhiAt+Phi01)
      Rhh(1,2)=-4.d0*(((5.d0/3.d0*g1q-g2q)/4.d0)**2*vuq*vdq
     .     +Ytq**2*muq*At*DDCOS(PhiAt+Phi01))
      Rhh(2,1)=Rhh(1,2)
      Rhh(3,3)=4.d0*Ytq**2*l**2/muq*vuq*vdq*At*DDCOS(PhiAt+Phi01)
      Rhh(1,3)=-4.d0*Ytq**2*l*vdq*At*DDCOS(PhiAt+Phi01)
      Rhh(3,1)=Rhh(1,3)
      Rhh(2,3)=4.d0*Ytq**2*l*(2d0*muq*vdq-vuq*At*DDCOS(PhiAt+Phi01))
      Rhh(3,2)=Rhh(2,3)
      Rhh(4,4)=4.d0*Ytq**2*muq*At*DDCOS(PhiAt+Phi01)
     .                                    /vdq/vuq*(vu**2+vd**2)
      Rhh(5,5)=4.d0*Ytq**2*l**2/muq*vdq*vuq*At*DDCOS(PhiAt+Phi01)
      Rhh(4,5)=4.d0*Ytq**2*l*At*DDCOS(PhiAt+Phi01)*dsqrt(vu**2+vd**2)
      Rhh(5,4)=Rhh(4,5)

      DO I=1,5
      DO J=I,5
      aux=2d0*Rhh(I,J)*Fsf1(MST2(1),MST2(2),QSTSB)
     .      +dTdh(I)*dTdh(J)*dlog(MST2(1)*MST2(2)/QSTSB**2)
     .      +(dTdh(I)*dRdh(J)+dTdh(J)*dRdh(I))/2d0
     .                        *dlog(MST2(2)/MST2(1))/dsqrt(R2)
     .      +dRdh(I)*dRdh(J)*Fsf2(MST2(1),MST2(2))/(4.d0*R2)
      MH02(I,J)=MH02(I,J)+3.d0/128.d0/Pi**2*aux
      MH02(J,I)=MH02(I,J)
      ENDDO
      ENDDO


c      II- Bottoms/Sbottoms

c      a) Bottoms

      aux=4.d0*Ybq**4*vdq**2*dlog(mbotq**2/QSTSB)
      MH02(2,2)=MH02(2,2)-3.d0/16.d0/Pi**2*aux

c      b) Sbottoms

      DO I=1,5
       dTdh(I)=0d0
       dRdh(I)=0d0
       DO J=1,5
        Rhh(I,J)=0d0
       ENDDO
      ENDDO

      R2=(MSQ3-MSD3+(-1d0/3.d0*g1q+g2q)/4.d0*
     .  (vuq**2-vdq**2))**2+4.d0*Ybq**2*(Ab**2*vdq**2+muq**2*vuq**2
     .                  -2d0*muq*Ab*vuq*vdq*DDCOS(PhiAb+Phi01))

      dTdh(1)=(g1q+g2q)/2d0*vuq
      dTdh(2)=4.d0*(Ybq**2-(g1q+g2q)/8.d0)*vdq

      dRdh(1)=8.d0*Ybq**2*(vuq*muq**2-muq*Ab*vdq*DDCOS(PhiAb+Phi01))
     .   +(-1d0/3.d0*g1q+g2q)*vuq*
     .   (MSQ3-MSD3+(-1d0/3.d0*g1q+g2q)/4.d0*(vuq**2-vdq**2))
      dRdh(2)=8.d0*Ybq**2*(vdq*Ab**2-muq*Ab*vuq*DDCOS(PhiAb+Phi01))
     .   -(-1d0/3.d0*g1q+g2q)*vdq*
     .   (MSQ3-MSD3+(-1d0/3.d0*g1q+g2q)/4.d0*(vuq**2-vdq**2))
      dRdh(3)=8.d0*Ybq**2*l
     .         *(muq*vuq**2-Ab*vuq*vdq*DDCOS(PhiAb+Phi01))
      dRdh(4)=8.d0*Ybq**2*muq*Ab*DDSIN(PhiAb+Phi01)*dsqrt(vu**2+vd**2)
      dRdh(5)=8.d0*Ybq**2*l*vuq*vdq*Ab*DDSIN(PhiAb+Phi01)

      Rhh(1,1)=4.d0*(((-1d0/3.d0*g1q+g2q)/4.d0*vuq)**2
     .       +Ybq**2*muq*Ab*vdq/vuq*DDCOS(PhiAb+Phi01))
      Rhh(2,2)=4.d0*(((-1d0/3.d0*g1q+g2q)/4.d0*vdq)**2
     .     +Ybq**2*muq*Ab*vuq/vdq*DDCOS(PhiAb+Phi01))
      Rhh(1,2)=-4.d0*(((-1d0/3.d0*g1q+g2q)/4.d0)**2*vuq*vdq
     .     +Ybq**2*muq*Ab*DDCOS(PhiAb+Phi01))
      Rhh(2,1)=Rhh(1,2)
      Rhh(3,3)=4.d0*Ybq**2*l**2/muq*vuq*vdq*Ab*DDCOS(PhiAb+Phi01)
      Rhh(1,3)=4.d0*Ybq**2*l*(2d0*muq*vuq-vdq*Ab*DDCOS(PhiAb+Phi01))
      Rhh(3,1)=Rhh(1,3)
      Rhh(2,3)=-4.d0*Ybq**2*l*vuq*Ab*DDCOS(PhiAb+Phi01)
      Rhh(3,2)=Rhh(2,3)
      Rhh(4,4)=4.d0*Ybq**2*muq*Ab*DDCOS(PhiAb+Phi01)
     .                       /vdq/vuq*(vu**2+vd**2)
      Rhh(5,5)=4.d0*Ybq**2*l**2/muq*vdq*vuq*Ab*DDCOS(PhiAb+Phi01)
      Rhh(4,5)=4.d0*Ybq**2*l*Ab*DDCOS(PhiAb+Phi01)*dsqrt(vu**2+vd**2)
      Rhh(5,4)=Rhh(4,5)

      DO I=1,5
      DO J=I,5
      aux=2d0*Rhh(I,J)*Fsf1(MSB2(1),MSB2(2),QSTSB)
     .      +dTdh(I)*dTdh(J)*dlog(MSB2(1)*MSB2(2)/QSTSB**2)
     .      +(dTdh(I)*dRdh(J)+dTdh(J)*dRdh(I))/2d0
     .                        *dlog(MSB2(2)/MSB2(1))/dsqrt(R2)
     .      +dRdh(I)*dRdh(J)*Fsf2(MSB2(1),MSB2(2))/(4.d0*R2)
      MH02(I,J)=MH02(I,J)+3.d0/128.d0/Pi**2*aux
      MH02(J,I)=MH02(I,J)
      ENDDO
      ENDDO


c      III- Taus/Staus

      Ytau=mtau/vdq

c      a) Taus

      aux=4.d0*Ytau**4*vdq**2*dlog(mtau**2/QSTSB)
      MH02(2,2)=MH02(2,2)-1d0/16.d0/Pi**2*aux

c      b) Staus

      DO I=1,5
       dTdh(I)=0d0
       dRdh(I)=0d0
       DO J=1,5
        Rhh(I,J)=0d0
       ENDDO
      ENDDO

      R2=(MSL3-MSE3+(-3.d0*g1q+g2q)/4.d0*(vuq**2-vdq**2)
     .     )**2+4.d0*Ytau**2*(Atau**2*vdq**2+muq**2*vuq**2
     .               -2d0*muq*Atau*vuq*vdq*DDCOS(PhiAtau+Phi01))

      dTdh(1)=(g1q+g2q)/2d0*vuq
      dTdh(2)=4.d0*(Ytau**2-(g1q+g2q)/8.d0)*vdq

      dRdh(1)=8.d0*Ytau**2*(vuq*muq**2-muq*Atau*vdq*
     .   DDCOS(PhiAtau+Phi01))+vuq*(-3.d0*g1q+g2q)*
     .   (MSL3-MSE3+(-3.d0*g1q+g2q)/4.d0*(vuq**2-vdq**2))
      dRdh(2)=8.d0*Ytau**2*(Atau**2*vdq-muq*Atau*vuq*
     .   DDCOS(PhiAtau+Phi01))-vdq*(-3.d0*g1q+g2q)*
     .   (MSL3-MSE3+(-3.d0*g1q+g2q)/4.d0*(vuq**2-vdq**2))
      dRdh(3)=8.d0*Ytau**2*l*(muq*vuq**2-Atau*vuq*vdq*
     .                              DDCOS(PhiAtau+Phi01))
      dRdh(4)=8.d0*Ytau**2*muq*Atau*DDSIN(PhiAtau+Phi01)
     .                                            *dsqrt(vu**2+vd**2)
      dRdh(5)=8.d0*Ytau**2*l*vuq*vdq*Atau*DDSIN(PhiAtau+Phi01)

      Rhh(1,1)=4.d0*(((-3.d0*g1q+g2q)/4.d0*vuq)**2
     .     +Ytau**2*muq*Atau*vdq/vuq*DDCOS(PhiAtau+Phi01))
      Rhh(2,2)=4.d0*(((-3.d0*g1q+g2q)/4.d0*vdq)**2
     .     +Ytau**2*muq*Atau*vuq/vdq*DDCOS(PhiAtau+Phi01))
      Rhh(1,2)=-4.d0*(((-3.d0*g1q+g2q)/4.d0)**2*vuq*vdq
     .     +Ytau**2*muq*Atau*DDCOS(PhiAtau+Phi01))
      Rhh(2,1)=Rhh(1,2)
      Rhh(3,3)=4.d0*Ytau**2*l**2/muq*vuq*vdq*Atau*DDCOS(PhiAtau+Phi01)
      Rhh(1,3)=4.d0*Ytau**2*l*(2d0*muq*vuq-vdq*Atau*
     .                                        DDCOS(PhiAtau+Phi01))
      Rhh(3,1)=Rhh(1,3)
      Rhh(2,3)=-4.d0*Ytau**2*l*vuq*Atau*DDCOS(PhiAtau+Phi01)
      Rhh(3,2)=Rhh(2,3)
      Rhh(4,4)=4.d0*Ytau**2*muq*Atau*DDCOS(PhiAtau+Phi01)
     .                       /vdq/vuq*(vu**2+vd**2)
      Rhh(5,5)=4.d0*Ytau**2*l**2/muq*vdq*vuq*Atau
     .                                        *DDCOS(PhiAtau+Phi01)
      Rhh(4,5)=4.d0*Ytau**2*l*Atau*DDCOS(PhiAtau+Phi01)
     .                                              *dsqrt(vu**2+vd**2)
      Rhh(5,4)=Rhh(4,5)

      DO I=1,5
      DO J=I,5
      aux=2d0*Rhh(I,J)*Fsf1(MSL2(1),MSL2(2),QSTSB)
     .      +dTdh(I)*dTdh(J)*dlog(MSL2(1)*MSL2(2)/QSTSB**2)
     .      +(dTdh(I)*dRdh(J)+dTdh(J)*dRdh(I))/2d0
     .                        *dlog(MSL2(2)/MSL2(1))/dsqrt(R2)
     .      +dRdh(I)*dRdh(J)*Fsf2(MSL2(1),MSL2(2))/(4.d0*R2)
      MH02(I,J)=MH02(I,J)+1d0/128.d0/Pi**2*aux
      MH02(J,I)=MH02(I,J)
      ENDDO
      ENDDO


c      IV- Sneutrinos

      aux=1d0/32d0/Pi**2*((-g1q-g2q)/2d0)**2*dlog(MSNT2/QSTSB)

      MH02(1,1)=MH02(1,1)+aux*vuq**2

      MH02(2,2)=MH02(2,2)+aux*vdq**2

      MH02(1,2)=MH02(1,2)-aux*vuq*vdq

      MH02(2,1)=MH02(1,2)


c            B: Corrections to the eff. parameters /3rd generation

c      I- Stops/Sbottoms

      lu=lu+6.d0/32d0/Pi**2*(
     . (Ytq**2-g1q/3.d0)**2*dlog(MSU3/QSTSB)
     .  +(g1q/6.d0)**2*dlog(MSD3/QSTSB)
     .  +(Ytq**4+Ytq**2/2d0*(g1q/3.d0-g2q)+2d0*(g1q/12d0)**2
     .   +2d0*(g2q/4.d0)**2)*dlog(MSQ3/QSTSB)
     .  +2d0*Ytq**2*AT**2*(Ytq**2-g1q/3.d0)*Fsf3(MSQ3,MSU3)
     .  +2d0*Ybq**2*g1q/6.d0*muq**2*Fsf3(MSQ3,MSD3)
     .  +2d0*Ytq**2*AT**2*(Ytq**2+g1q/12d0-g2q/4.d0)
     .                                       *Fsf3(MSU3,MSQ3)
     .  +2d0*Ybq**2*(g1q/3.d0+g2q)/4.d0*muq**2*Fsf3(MSD3,MSQ3)
     .  +(Ytq*AT)**4*Fsf7(MSQ3,MSU3)+(Ybq*muq)**4*Fsf7(MSQ3,MSD3))

      ld=ld+6.d0/32d0/Pi**2*
     . ((g1q/3.d0)**2*dlog(MSU3/QSTSB)
     .  +(Ybq**2-g1q/6.d0)**2*dlog(MSD3/QSTSB)
     .  +(Ybq**4-Ybq**2/2d0*(g1q/3.d0+g2q)+2d0*(g1q/12d0)**2
     .   +2d0*(g2q/4.d0)**2)*dlog(MSQ3/QSTSB)
     .  +2d0*Ytq**2*g1q/3.d0*muq**2*Fsf3(MSQ3,MSU3)
     .  +2d0*Ybq**2*(Ybq**2-g1q/6.d0)*Ab**2*Fsf3(MSQ3,MSD3)
     .  -2d0*Ytq**2*(g1q/3.d0-g2q)/4.d0*muq**2*Fsf3(MSU3,MSQ3)
     .  +2d0*Ybq**2*(Ybq**2-g1q/12d0-g2q/4.d0)*Ab**2
     .                                         *Fsf3(MSD3,MSQ3)
     .  +(Ytq*muq)**4*Fsf7(MSQ3,MSU3)+(Ybq*AB)**4*Fsf7(MSQ3,MSD3))

      l3=l3+6.d0/32d0/Pi**2*
     . ((Ytq**2-g1q/3.d0)*g1q/3.d0*dlog(MSU3/QSTSB)
     .  +(Ybq**2-g1q/6.d0)*g1q/6.d0*dlog(MSD3/QSTSB)
     .  +(Ytq**2*Ybq**2-Ytq**2/4.d0*(g1q/3.d0+g2q)+Ybq**2/4.d0
     .   *(g1q/3.d0-g2q)-2d0*(g1q/12d0)**2
     .   +2d0*(g2q/4.d0)**2)*dlog(MSQ3/QSTSB)
     .  +Ytq**2*((Ytq**2-g1q/3.d0)*muq**2+g1q/3.d0*At**2)
     .                                      *Fsf3(MSQ3,MSU3)
     .  +Ybq**2*((Ybq**2-g1q/6.d0)*muq**2+g1q/6.d0*Ab**2)
     .                                      *Fsf3(MSQ3,MSD3)
     .  +Ytq**2*((Ybq**2-g1q/12d0-g2q/4.d0)*At**2
     .          +(g1q/3.d0+g2q)/4.d0*muq**2)*Fsf3(MSU3,MSQ3)
     .  +Ybq**2*((Ytq**2+g1q/12d0-g2q/4.d0)*Ab**2-(g1q/3.d0-g2q)
     .                            *muq**2/4.d0)*Fsf3(MSD3,MSQ3)
     .  +Ytq**2*Ybq**2*Fsf1(MSU3,MSD3,QSTSB)
     .  +2d0*Ytq**2*Ybq**2*(At*Ab*DDCOS(PhiAt-PhiAb)-muq**2)
     .    *Fsf5(MSU3,MSD3,MSQ3)+Ytq**4*At**2*muq**2*Fsf7(MSQ3,MSU3)
     .  +Ybq**4*Ab**2*muq**2*Fsf7(MSQ3,MSD3)+Ytq**2*Ybq**2*
     .   (At**2*Ab**2-2d0*At*Ab*muq**2*DDCOS(PhiAt-PhiAb)+muq**4)
     .   *Fsf6(MSU3,MSD3,MSQ3))

      l4=l4+6.d0/32d0/Pi**2*
     . (-(Ytq**2-g2q/2d0)*(Ybq**2-g2q/2d0)*dlog(MSQ3/QSTSB)
     .  +Ytq**2*((Ytq**2-g2q/2d0)*muq**2-(Ybq**2-g2q/2d0)*At**2)
     .                                      *Fsf3(MSU3,MSQ3)
     .  +Ybq**2*((Ybq**2-g2q/2d0)*muq**2-(Ytq**2-g2q/2d0)*Ab**2)
     .                                      *Fsf3(MSD3,MSQ3)
     .  -Ytq**2*Ybq**2*Fsf1(MSU3,MSD3,QSTSB)
     .  -2d0*Ytq**2*Ybq**2*(At*Ab*DDCOS(PhiAt-PhiAb)-muq**2)
     .    *Fsf5(MSU3,MSD3,MSQ3)+Ytq**4*At**2*muq**2*Fsf7(MSQ3,MSU3)
     .  +Ybq**4*Ab**2*muq**2*Fsf7(MSQ3,MSD3)-Ytq**2*Ybq**2*
     .   (At**2*Ab**2-2d0*At*Ab*muq**2*DDCOS(PhiAt-PhiAb)+muq**4)
     .   *Fsf6(MSU3,MSD3,MSQ3))

      Rel5=Rel5+6.d0/32d0/Pi**2*muq**2*
     . (Ytq**4*At**2*DDCOS(2d0*(PhiAt+Phi01))*Fsf7(MSQ3,MSU3)
     . +Ybq**4*Ab**2*DDCOS(2d0*(PhiAb+Phi01))*Fsf7(MSQ3,MSD3))

      Iml5=Iml5+6/32d0/Pi**2*muq**2*
     . (Ytq**4*At**2*DDSIN(2d0*(PhiAt+Phi01))*Fsf7(MSQ3,MSU3)
     . +Ybq**4*Ab**2*DDSIN(2d0*(PhiAb+Phi01))*Fsf7(MSQ3,MSD3))

      Rel6=Rel6+6.d0/32d0/Pi**2*muq*
     . (Ytq**2*At*(Ytq**2-g1q/3.d0)*DDCOS(PhiAt+Phi01)
     .                                      *Fsf3(MSQ3,MSU3)
     .  +Ybq**2*Ab*g1q/6.d0*DDCOS(PhiAb+Phi01)*Fsf3(MSQ3,MSD3)
     .  +Ytq**2*At*(Ytq**2+g1q/12d0-g2q/4.d0)*DDCOS(PhiAt+Phi01)
     .                                      *Fsf3(MSU3,MSQ3)
     .  +Ybq**2*Ab*(g1q/3.d0+g2q)/4.d0*DDCOS(PhiAb+Phi01)
     .                                      *Fsf3(MSD3,MSQ3)
     .  +Ytq**4*At**3*DDCOS(PhiAt+Phi01)*Fsf7(MSQ3,MSU3)
     .  +Ybq**4*Ab*muq**2*DDCOS(PhiAb+Phi01)*Fsf7(MSQ3,MSD3))

      Iml6=Iml6+6.d0/32d0/Pi**2*muq*
     . (Ytq**2*At*(Ytq**2-g1q/3.d0)*DDSIN(PhiAt+Phi01)
     .                                      *Fsf3(MSQ3,MSU3)
     .  +Ybq**2*Ab*g1q/6.d0*DDSIN(PhiAb+Phi01)*Fsf3(MSQ3,MSD3)
     .  +Ytq**2*At*(Ytq**2+g1q/12d0-g2q/4.d0)*DDSIN(PhiAt+Phi01)
     .                                      *Fsf3(MSU3,MSQ3)
     .  +Ybq**2*Ab*(g1q/3.d0+g2q)/4.d0*DDSIN(PhiAb+Phi01)
     .                                      *Fsf3(MSD3,MSQ3)
     .  +Ytq**4*At**3*DDSIN(PhiAt+Phi01)*Fsf7(MSQ3,MSU3)
     .  +Ybq**4*Ab*muq**2*DDSIN(PhiAb+Phi01)*Fsf7(MSQ3,MSD3))

      Rel7=Rel7+6.d0/32d0/Pi**2*muq*
     . (Ytq**2*At*g1q/3.d0*DDCOS(PhiAt+Phi01)*Fsf3(MSQ3,MSU3)
     .  +Ybq**2*Ab*(Ybq**2-g1q/6.d0)*DDCOS(PhiAb+Phi01)
     .                                      *Fsf3(MSQ3,MSD3)
     .  -Ytq**2*At*(g1q/3.d0-g2q)/4.d0*DDCOS(PhiAt+Phi01)
     .                                      *Fsf3(MSU3,MSQ3)
     .  +Ybq**2*Ab*(Ybq**2-g1q/12d0-g2q/4.d0)*DDCOS(PhiAb+Phi01)
     .                                      *Fsf3(MSD3,MSQ3)
     .  +Ytq**4*At*muq**2*DDCOS(PhiAt+Phi01)*Fsf7(MSQ3,MSU3)
     .  +Ybq**4*Ab**3*DDCOS(PhiAb+Phi01)*Fsf7(MSQ3,MSD3))

      Iml7=Iml7+6.d0/32d0/Pi**2*muq*
     . (Ytq**2*At*g1q/3.d0*DDSIN(PhiAt+Phi01)*Fsf3(MSQ3,MSU3)
     .  +Ybq**2*Ab*(Ybq**2-g1q/6.d0)*DDSIN(PhiAb+Phi01)
     .                                      *Fsf3(MSQ3,MSD3)
     .  -Ytq**2*At*(g1q/3.d0-g2q)/4.d0*DDSIN(PhiAt+Phi01)
     .                                      *Fsf3(MSU3,MSQ3)
     .  +Ybq**2*Ab*(Ybq**2-g1q/12d0-g2q/4.d0)*DDSIN(PhiAb+Phi01)
     .                                      *Fsf3(MSD3,MSQ3)
     .  +Ytq**4*At*muq**2*DDSIN(PhiAt+Phi01)*Fsf7(MSQ3,MSU3)
     .  +Ybq**4*Ab**3*DDSIN(PhiAb+Phi01)*Fsf7(MSQ3,MSD3))

      RAud=RAud+6.d0*l/32d0/Pi**2*
     .  (Ytq**2*At*DDCOS(PhiAt+Phi01)*Fsf1(MSQ3,MSU3,QSTSB)
     .  +Ybq**2*Ab*DDCOS(PhiAb+Phi01)*Fsf1(MSQ3,MSD3,QSTSB))

      lPu=lPu+6.d0/32d0/Pi**2*Ybq**2*l**2*Fsf1(MSQ3,MSU3,QSTSB)

      lPd=lPd+6.d0/32d0/Pi**2*Ytq**2*l**2*Fsf1(MSQ3,MSU3,QSTSB)

      s=muq/l


c      II- Staus/Sneutrinos

      lu=lu+2d0/32d0/Pi**2*(
     .  (g1q/2d0)**2*dlog(MSE3/QSTSB)
     .  +2d0*((g1q/4.d0)**2+(g2q/4.d0)**2)*dlog(MSL3/QSTSB)
     .  +2d0*Ytau**2*g1q/2d0*muq**2*Fsf3(MSL3,MSE3)
     .  +2d0*Ytau**2*(-g1q+g2q)/4.d0*muq**2*Fsf3(MSE3,MSL3)
     .  +(Ytau*muq)**4*Fsf7(MSL3,MSE3))

      ld=ld+2d0/32d0/Pi**2*(
     .  (Ytau**2-g1q/2d0)**2*dlog(MSE3/QSTSB)
     .  +(Ytau**4+Ytau**2/2d0*(g1q-g2q)+2d0*(g1q/4.d0)**2
     .                     +2d0*(g2q/4.d0)**2)*dlog(MSL3/QSTSB)
     .  +2d0*Ytau**2*(Ytau**2-g1q/2d0)*Atau**2*Fsf3(MSL3,MSE3)
     .  +2d0*Ytau**2*(Ytau**2+g1q/4.d0-g2q/4.d0)*Atau**2
     .                                         *Fsf3(MSE3,MSL3)
     .  +(Ytau*Atau)**4*Fsf7(MSL3,MSE3))

      l3=l3+2d0/32d0/Pi**2*(
     .  (Ytau**2-g1q/2d0)*g1q/2d0*dlog(MSE3/QSTSB)
     .  -(Ytau**2+g1q/2d0-g2q/2d0)/4.d0*(g1q+g2q)*dlog(MSL3/QSTSB)
     .  +Ytau**2*((Ytau**2-g1q/2d0)*muq**2+g1q/2d0*Atau**2)
     .                                      *Fsf3(MSL3,MSE3)
     .  +Ytau**2*(g1q+g2q)/4.d0*(muq**2-Atau**2)
     .                                      *Fsf3(MSE3,MSL3)
     .  +Ytau**4*Atau**2*muq**2*Fsf7(MSL3,MSE3))

      l4=l4+2d0/32d0/Pi**2*(
     .  g2q/2d0*(Ytau**2-g2q/2d0)*dlog(MSL3/QSTSB)
     .  +Ytau**2*((Ytau**2-g2q/2d0)*muq**2+g2q/2d0*Atau**2)
     .                                      *Fsf3(MSE3,MSL3)
     .  +Ytau**4*Atau**2*muq**2*Fsf7(MSL3,MSE3))

      Rel5=Rel5+2d0/32d0/Pi**2*muq**2*Ytau**4*Atau**2
     . *DDCOS(2d0*(PhiAtau+Phi01))*Fsf7(MSL3,MSE3)

      Iml5=Iml5+2d0/32d0/Pi**2*muq**2*Ytau**4*Atau**2
     . *DDSIN(2d0*(PhiAtau+Phi01))*Fsf7(MSL3,MSE3)

      Rel6=Rel6+2d0/32d0/Pi**2*muq*Ytau**2*Atau*
     .  DDCOS(PhiAtau+Phi01)*(g1q/2d0*Fsf3(MSL3,MSE3)
     .  +(-g1q+g2q)/4.d0*Fsf3(MSE3,MSL3)
     .  +Ytau**2*muq**2*Fsf7(MSL3,MSE3))

      Iml6=Iml6+2d0/32d0/Pi**2*muq*Ytau**2*Atau*
     .  DDSIN(PhiAtau+Phi01)*(g1q/2d0*Fsf3(MSL3,MSE3)
     .  +(-g1q+g2q)/4.d0*Fsf3(MSE3,MSL3)
     .  +Ytau**2*muq**2*Fsf7(MSL3,MSE3))

      Rel7=Rel7+2d0/32d0/Pi**2*muq*Ytau**2*Atau*
     .  DDCOS(PhiAtau+Phi01)*((Ytau**2-g1q/2d0)*Fsf3(MSL3,MSE3)
     .  +(Ytau**2+g1q/4.d0-g2q/4.d0)*Fsf3(MSE3,MSL3)
     .  +Ytau**2*Atau**2*Fsf7(MSL3,MSE3))

      Iml7=Iml7+2d0/32d0/Pi**2*muq*Ytau**2*Atau*
     .  DDSIN(PhiAtau+Phi01)*((Ytau**2-g1q/2d0)*Fsf3(MSL3,MSE3)
     .  +(Ytau**2+g1q/4.d0-g2q/4.d0)*Fsf3(MSE3,MSL3)
     .  +Ytau**2*Atau**2*Fsf7(MSL3,MSE3))

      RAud=RAud+2d0*l*Ytau**2*Atau/32d0/Pi**2*
     .  DDCOS(PhiAtau+Phi01)*Fsf1(MSL3,MSE3,QSTSB)

      lPu=lPu+2d0/32d0/Pi**2*l**2*Ytau**2*Fsf1(MSL3,MSE3,QSTSB)

c      Minimization conditions

        MHuS=MHuS
     .    +3d0*mtopq**4/vuq**2*(dlog(mtopq**2/QSTSB)-1d0)/8d0/Pi**2
     .    -(3d0*(Ytq**2-(g1q+g2q)/8d0)
     .      *(MST2(1)*(dlog(MST2(1)/QSTSB)-1d0)
     .       +MST2(2)*(dlog(MST2(2)/QSTSB)-1d0))
     .     +3d0*Fsf1(MST2(1),MST2(2),QSTSB)
     .      *(Ytq**2*AT*(AT-muq*vdq/vuq*DDCOS(PhiAT+Phi01))
     .        +(5d0/3d0*g1q-g2q)/8d0
     .         *(MSQ3-MSU3+(5d0/3d0*g1q-g2q)/4d0*(vuq**2-vdq**2)))
     .     +3d0*(g1q+g2q)/8d0
     .      *(MSB2(1)*(dlog(MSB2(1)/QSTSB)-1d0)
     .       +MSB2(2)*(dlog(MSB2(2)/QSTSB)-1d0))
     .     +3d0*Fsf1(MSB2(1),MSB2(2),QSTSB)
     .      *(Ybq**2*muq*(muq-AB*vdq/vuq*DDCOS(PhiAB+Phi01))
     .        +(g1q/3+g2q)/8d0
     .         *(MSQ3-MSD3-(g1q/3d0+g2q)/4d0*(vuq**2-vdq**2)))
     .     +(g1q+g2q)/8d0
     .      *(MSL2(1)*(dlog(MSL2(1)/QSTSB)-1d0)
     .       +MSL2(2)*(dlog(MSL2(2)/QSTSB)-1d0))
     .     +Fsf1(MSL2(1),MSL2(2),QSTSB)
     .      *((mtau/vd)**2*muq*(muq-ATAU*vdq/vuq*DDCOS(PhiATAU+Phi01))
     .        +(2d0*g1q-g2q)/8d0
     .         *(MSL3-MSE3-(2d0*g1q-g2q)/4d0*(vuq**2-vdq**2)))
     .     -(g1q+g2q)/8d0*MSNT2*(dlog(MSNT2/QSTSB)-1d0))/16d0/Pi**2

        MHdS=MHdS
     .    +(3d0*mbotq**4/vdq**2*(dlog(mbotq**2/QSTSB)-1d0)
     .     +(mtau/vd)**4*vdq**2*(dlog(mtau**2/QSTSB)-1d0))/8d0/Pi**2
     .    -(3d0*(g1q+g2q)/8d0
     .      *(MST2(1)*(dlog(MST2(1)/QSTSB)-1d0)
     .       +MST2(2)*(dlog(MST2(2)/QSTSB)-1d0))
     .     +3d0*Fsf1(MST2(1),MST2(2),QSTSB)
     .      *(Ytq**2*muq*(muq-AT*vuq/vdq*DDCOS(PhiAT+Phi01))
     .        -(5d0/3d0*g1q-g2q)/8d0
     .         *(MSQ3-MSU3+(5d0/3d0*g1q-g2q)/4d0*(vuq**2-vdq**2)))
     .     +3d0*(Ybq**2-(g1q+g2q)/8d0)
     .      *(MSB2(1)*(dlog(MSB2(1)/QSTSB)-1d0)
     .       +MSB2(2)*(dlog(MSB2(2)/QSTSB)-1d0))
     .     +3d0*Fsf1(MSB2(1),MSB2(2),QSTSB)
     .      *(Ybq**2*AB*(AB-muq*vuq/vdq*DDCOS(PhiAB+Phi01))
     .        -(g1q/3+g2q)/8d0
     .         *(MSQ3-MSD3-(g1q/3d0+g2q)/4d0*(vuq**2-vdq**2)))
     .     +((mtau/vd)**2-(g1q+g2q)/8d0)
     .      *(MSL2(1)*(dlog(MSL2(1)/QSTSB)-1d0)
     .       +MSL2(2)*(dlog(MSL2(2)/QSTSB)-1d0))
     .     +Fsf1(MSL2(1),MSL2(2),QSTSB)
     .      *((mtau/vd)**2*ATAU*(ATAU-muq*vuq/vdq*DDCOS(PhiATAU+Phi01))
     .        -(2d0*g1q-g2q)/8d0
     .         *(MSL3-MSE3-(2d0*g1q-g2q)/4d0*(vuq**2-vdq**2)))
     .     +(g1q+g2q)/8d0*MSNT2*(dlog(MSNT2/QSTSB)-1d0))/16d0/Pi**2

        MSS=MSS
     .    -(3d0*Fsf1(MST2(1),MST2(2),QSTSB)
     .      *Ytq**2*l**2*vdq*(vdq-AT*vuq/muq)
     .     +3d0*Fsf1(MSB2(1),MSB2(2),QSTSB)
     .      *Ybq**2*l**2*vuq*(vuq-AB*vdq/muq*DDCOS(PhiAB+Phi01))
     .     +Fsf1(MSL2(1),MSL2(2),QSTSB)
     .      *(mtau/vd)**2*l**2*vuq
     .              *(vuq-ATAU*vdq/muq*DDCOS(PhiATAU+Phi01)))/16d0/Pi**2

        IAl=IAl
     .    -(3d0*Fsf1(MST2(1),MST2(2),QSTSB)
     .      *Ytq**2*AT*DDSIN(PhiAT+Phi01)
     .     +3d0*Fsf1(MSB2(1),MSB2(2),QSTSB)
     .      *Ybq**2*AB*DDSIN(PhiAB+Phi01)
     .     +Fsf1(MSL2(1),MSL2(2),QSTSB)
     .      *(mtau/vd)**2*ATAU*DDSIN(PhiATAU+Phi01))/16d0/Pi**2

c            C: Corrections to the neutral mass-matrix / 1st & 2nd generations

      aux=2d0/32d0/Pi**2*(
     .   3.d0*((g1q/3.d0-g2q)/2d0)**2*dlog(MSU2(1)/QSTSB)
     .  +3.d0*((-4.d0*g1q/3.d0)/2d0)**2*dlog(MSU2(2)/QSTSB)
     .  +3.d0*((g1q/3.d0+g2q)/2d0)**2*dlog(MSD2(1)/QSTSB)
     .  +3.d0*((2d0*g1q/3.d0)/2d0)**2*dlog(MSD2(2)/QSTSB)
     .  +((-g1q+g2q)/2d0)**2*dlog(MSE2(1)/QSTSB)
     .  +((2d0*g1q)/2d0)**2*dlog(MSE2(2)/QSTSB)
     .  +((-g1q-g2q)/2d0)**2*dlog(MSNE2/QSTSB))

      MH02(1,1)=MH02(1,1)+aux*vuq**2

      MH02(2,2)=MH02(2,2)+aux*vdq**2

      MH02(1,2)=MH02(1,2)-aux*vuq*vdq

      MH02(2,1)=MH02(1,2)


c            D: Corrections to the eff. parameters / 1st & 2nd generations

c      I- Squarks

      lu=lu+12d0/32d0/Pi**2*(
     . (g1q/3.d0)**2*dlog(MSU1/QSTSB)
     . +(g1q/6.d0)**2*dlog(MSD1/QSTSB)
     . +2d0*((g1q/12d0)**2+(g2q/4.d0)**2)*dlog(MSQ1/QSTSB))

      ld=ld+12d0/32d0/Pi**2*(
     . (g1q/3.d0)**2*dlog(MSU1/QSTSB)
     . +(g1q/6.d0)**2*dlog(MSD1/QSTSB)
     . +2d0*((g1q/12d0)**2+(g2q/4.d0)**2)*dlog(MSQ1/QSTSB))

      l3=l3+12d0/32d0/Pi**2*(
     . -(g1q/3.d0)**2*dlog(MSU1/QSTSB)
     . -(g1q/6.d0)**2*dlog(MSD1/QSTSB)
     . +2d0*(-(g1q/12d0)**2+(g2q/4.d0)**2)*dlog(MSQ1/QSTSB))

      l4=l4-12d0/32d0/Pi**2*(g2q/2d0)**2*dlog(MSQ1/QSTSB)

!      Rel5=Rel5+0d0
!      Iml5=Iml5+0d0
!      Rel6=Rel6+0d0
!      Iml6=Iml6+0d0
!      Rel7=Rel7+0d0
!      Iml7=Iml7+0d0
!      RAud=RAud+0d0

c      II- 2HDM parameters 1st/2nd gen sleptons

      lu=lu+4.d0/32d0/Pi**2*(
     . (g1q/2d0)**2*dlog(MSE1/QSTSB)
     . +2d0*((g1q/4.d0)**2+(g2q/4.d0)**2)*dlog(MSL1/QSTSB))

      ld=ld+4.d0/32d0/Pi**2*(
     . (g1q/2d0)**2*dlog(MSE1/QSTSB)
     . +2d0*((g1q/4.d0)**2+(g2q/4.d0)**2)*dlog(MSL1/QSTSB))

      l3=l3+4.d0/32d0/Pi**2*(
     . -(g1q/2d0)**2*dlog(MSE1/QSTSB)
     . -2d0*((g1q/4.d0)**2-(g2q/4.d0)**2)*dlog(MSL1/QSTSB))

      l4=l4-4.d0/32d0/Pi**2*(g2q/2d0)**2*dlog(MSL1/QSTSB)

!      Rel5=Rel5+0d0
!      Iml5=Iml5+0d0
!      Rel6=Rel6+0d0
!      Iml6=Iml6+0d0
!      Rel7=Rel7+0d0
!      Iml7=Iml7+0d0
!      RAud=RAud+0d0

c      Minimization conditions

        MHuS=MHuS+2d0*
     .     (-3d0*(g1q+g2q)/8d0
     .      *(MSU2(1)*(dlog(MSU2(1)/QSTSB)-1d0)
     .       +MSU2(2)*(dlog(MSU2(2)/QSTSB)-1d0))
     .     +3d0*Fsf1(MSU2(1),MSU2(2),QSTSB)
     .      *((5d0/3d0*g1q-g2q)/8d0
     .         *(MSQ1-MSU1+(5d0/3d0*g1q-g2q)/4d0*(vuq**2-vdq**2)))
     .     +3d0*(g1q+g2q)/8d0
     .      *(MSD2(1)*(dlog(MSD2(1)/QSTSB)-1d0)
     .       +MSD2(2)*(dlog(MSD2(2)/QSTSB)-1d0))
     .     +3d0*Fsf1(MSD2(1),MSD2(2),QSTSB)
     .      *((g1q/3+g2q)/8d0
     .         *(MSQ1-MSD1-(g1q/3d0+g2q)/4d0*(vuq**2-vdq**2)))
     .     +(g1q+g2q)/8d0
     .      *(MSE2(1)*(dlog(MSE2(1)/QSTSB)-1d0)
     .       +MSE2(2)*(dlog(MSE2(2)/QSTSB)-1d0))
     .     +Fsf1(MSE2(1),MSE2(2),QSTSB)
     .      *((2d0*g1q-g2q)/8d0
     .         *(MSL1-MSE1-(2d0*g1q-g2q)/4d0*(vuq**2-vdq**2)))
     .     -(g1q+g2q)/8d0*MSNE2*(dlog(MSNE2/QSTSB)-1d0))/16d0/Pi**2

        MHdS=MHdS+2d0*
     .     (-3d0*(g1q+g2q)/8d0
     .      *(MSU2(1)*(dlog(MSU2(1)/QSTSB)-1d0)
     .       +MSU2(2)*(dlog(MSU2(2)/QSTSB)-1d0))
     .     +3d0*Fsf1(MSU2(1),MSU2(2),QSTSB)
     .      *(-(5d0/3d0*g1q-g2q)/8d0
     .         *(MSQ1-MSU1+(5d0/3d0*g1q-g2q)/4d0*(vuq**2-vdq**2)))
     .     +3d0*(Ybq**2-(g1q+g2q)/8d0)
     .      *(MSD2(1)*(dlog(MSD2(1)/QSTSB)-1d0)
     .       +MSD2(2)*(dlog(MSD2(2)/QSTSB)-1d0))
     .     +3d0*Fsf1(MSD2(1),MSD2(2),QSTSB)
     .      *(-(g1q/3+g2q)/8d0
     .         *(MSQ1-MSD1-(g1q/3d0+g2q)/4d0*(vuq**2-vdq**2)))
     .     -(g1q+g2q)/8d0
     .      *(MSE2(1)*(dlog(MSE2(1)/QSTSB)-1d0)
     .       +MSE2(2)*(dlog(MSE2(2)/QSTSB)-1d0))
     .     +Fsf1(MSE2(1),MSE2(2),QSTSB)
     .      *(-(2d0*g1q-g2q)/8d0
     .         *(MSL1-MSE1-(2d0*g1q-g2q)/4d0*(vuq**2-vdq**2)))
     .     +(g1q+g2q)/8d0*MSNE2*(dlog(MSNE2/QSTSB)-1d0))/16d0/Pi**2


c            E: Leading 2-loop corrections (Lead. Log Approx.)

      aux=3.d0/256.d0/Pi**4*
     . Ytq**4*((dlog(QSTSB/mtopq**2))**2
     . *(64.d0*Pi*ALSQ+4.d0/3.d0*g1q-3.d0*sinb**2*Ytq**2
     . +3.d0*cosb**2*Ybq**2)+
     . ((dlog(MA2/mtopq**2))**2-(dlog(QSTSB/mtopq**2))**2)
     . *(3.d0*cosb**2*Ytq**2+(3.d0*cosb**2+1d0)*Ybq**2))

      R2=3.d0/256.d0/Pi**4*
     . Ybq**4*(dlog(QSTSB/mtopq**2)**2
     . *(64.d0*Pi*ALSQ-2d0/3.d0*g1q+3.d0*sinb**2*Ytq**2
     . -3.d0*cosb**2*Ybq**2)+
     . (dlog(MA2/mtopq**2)**2-dlog(QSTSB/mtopq**2)**2)
     . *(3.d0*sinb**2*Ybq**2+(3.d0*sinb**2+1d0)*Ytq**2))

      MH02(1,1)=MH02(1,1)+2d0*aux*vuq**2

      MH02(2,2)=MH02(2,2)+2d0*R2*vdq**2

      lu=lu+aux

      ld=ld+R2

c      Minimization conditions

        MHuS=MHuS-aux*vuq**2

        MHdS=MHdS-R2*vdq**2

c      s=muq/l

c      print*,'MH02_1*',(RAud+RlPM*s)*s*vdq/vuq+2d0*lu*vuq**2
c     c +Rel7*vdq**3/vuq-3.d0*Rel6*vuq*vdq,
c     c -(RAud+RlPM*s)*s+2d0*(l3+l4+Rel5)*vuq*vdq
c     c -3.d0*(Rel6*vuq**2+Rel7*vdq**2),
c     c -(RAud+2d0*RlPM*s)*vdq+2d0*lPU*s*vuq,
c     c (2d0*Iml6*vuq-Iml5*vdq)*dsqrt(vu**2+vd**2),
c     c (-Iml5*vuq*vdq+2d0*Iml6*vuq**2-3.d0*IlPM*s**2)*vdq/s
c      print*,'MH02_2*',-(RAud+RlPM*s)*s+2d0*(l3+l4+Rel5)*vuq*vdq
c     c -3.d0*(Rel6*vuq**2+Rel7*vdq**2),
c     c (RAud+RlPM*s)*s*vuq/vdq+2d0*ld*vdq**2+Rel6*vuq**3/vdq
c     c -3.d0*Rel7*vdq*vuq,
c     c -(RAud+2d0*RlPM*s)*vuq+2d0*lPD*s*vdq,
c     c (2d0*Iml7*vdq-Iml5*vuq)*dsqrt(vu**2+vd**2),
c     c (-Iml5*vuq*vdq+2d0*Iml7*vdq**2-3.d0*IlPM*s**2)*vuq/s
c      print*,'MH02_3*',-(RAud+2d0*RlPM*s)*vdq+2d0*lPU*s*vuq,
c     c -(RAud+2d0*RlPM*s)*vuq+2d0*lPD*s*vdq,
c     c (RAs+4*K2*s)*s+RAud*vuq*vdq/s,
c     c (IlPM*s**2-Iml5*vuq*vdq)*dsqrt(vu**2+vd**2)/s,
c     c (4.d0*IlPM*s**2-Iml5*vuq*vdq)*vuq*vdq/s**2
c      print*,'MH02_4*',(2d0*Iml6*vuq-Iml5*vdq)*dsqrt(vu**2+vd**2),
c     c (2d0*Iml7*vdq-Iml5*vuq)*dsqrt(vu**2+vd**2),
c     c (IlPM*s**2-Iml5*vuq*vdq)*dsqrt(vu**2+vd**2)/s,
c     c ((RAud+RlPM*s)*s-2d0*Rel5*vuq*vdq+Rel6*vuq**2+Rel7*vdq**2)
c     c *(vu**2+vd**2)/vuq/vdq/(ZHu*ZHd),
c     c (RAud-2d0*RlPM*s)*dsqrt(vu**2+vd**2)/dsqrt(Zs*ZHu*ZHd)
c      print*,'MH02_5*',(-Iml5*vuq*vdq+2d0*Iml6*vuq**2
c     c -3.d0*IlPM*s**2)*vdq/s,(-Iml5*vuq*vdq+2d0*Iml7*vdq**2
c     c -3.d0*IlPM*s**2)*vuq/s,(4.d0*IlPM*s**2-Iml5*vuq*vdq)
c     c *vuq*vdq/s**2,(RAud-2d0*RlPM*s)
c     c *dsqrt(vu**2+vd**2),-3.d0*RAs*s+(RAud+4.d0*RlPM*s)*vu*vd/s


c            F: Corrections to the charged-Higgs mass

      aux=-3.d0*Ytq**2*Ybq**2*(vu**2+vd**2)/(8.d0*Pi**2)
     .                        *Fsf1(mtopq**2,mbotq**2,QSTSB)

      MHC2=(Rm3+(RAud+RAudt+(RlPM+RlPMt+RlM)*muq/l)*muq/l
     .     -(l4+Rel5)*vuq*vdq+Rel6*vuq**2+Rel7*vdq**2+aux)
     .                 *(vu**2+vd**2)/vuq/vdq              !/(ZHu*ZHd)

c      print*,'MH02_1*',MH02(1,1),MH02(1,2),MH02(1,3),MH02(1,4),MH02(1,5)
c      print*,'MH02_2*',MH02(2,1),MH02(2,2),MH02(2,3),MH02(2,4),MH02(2,5)
c      print*,'MH02_3*',MH02(3,1),MH02(3,2),MH02(3,3),MH02(3,4),MH02(3,5)
c      print*,'MH02_4*',MH02(4,1),MH02(4,2),MH02(4,3),MH02(4,4)
c     c /(ZHu*ZHd),MH02(4,5)/dsqrt(Zs*ZHu*ZHd)
c      print*,'MH02_5*',MH02(5,1),MH02(5,2),MH02(5,3),MH02(5,4),MH02(5,5)
c      print*,'MHC2',MHC2                                   !/(ZHu*ZHd)

      RETURN
      END


      DOUBLE PRECISION function Fsf1(x,y,z)

c            ->Fsf1(m1^2,m2^2,Q^2)
      
      IMPLICIT NONE
      DOUBLE PRECISION x,y,z,aux

      IF(min(x,y).ge.1.d-10)THEN
       IF(dabs(x-y).ge.1.d-10)THEN
        aux=(y*dlog(y/z)-x*dlog(x/z))/(y-x)-1d0
       ELSE
        aux=dlog(x/z)
       ENDIF
      ELSEIF(min(x,y).le.1.d-10)THEN
       aux=dlog(Max(x,y)/z)-1d0
      ENDIF

      Fsf1=aux

      RETURN
      END

************************************************************************************************

      SUBROUTINE MHIGGSLOOP_INOS_CPV()

c         One-loop corrections to the Higgs potential
c                 - chargino + neutralino contribution
c      - The one-loop corrections from charginos/neutralinos to the parameters
c        of the effective Higgs potential are added and stored in the common
c        EFFPOTPAR.
c      - The corresponding corrections to the squared-mass matrices of the Higgs 
c        states are added and recorded in the common SQUHIMASSM.

      IMPLICIT NONE

      DOUBLE PRECISION Pi,aux,Mbi2,Mwi2,Mhi2,Msi2
      DOUBLE PRECISION fsf0,fsf1,fsf3,fsf4,fsf5,fsf6,fsf7
      DOUBLE PRECISION ludchi,l34chi,l4chi,Rel5chi,Iml5chi,Rel67chi,
     . Iml67chi
      DOUBLE PRECISION Rm3chi,Im3chi,RAudchi,RAudtchi,IAudtchi,RlPMchi,
     . IlPMchi,RlPMtchi,IlPMtchi,RlMchi,IlMchi,RAqschi,IAqschi,Rltqschi,
     . Iltqschi,lPqschi
      DOUBLE PRECISION dM2dHs,RdAuddHs,IdAuddHs,RdlPMdHs,IdlPMdHs
      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION mur,M1r,M2r,msi
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION l,k,Alcos1,Akcos2,muq,nuq
      DOUBLE PRECISION MH02(5,5),MHC2
      DOUBLE PRECISION lu,ld,l3,l4,Rel5,Iml5,Rel6,Iml6,Rel7,Iml7,RAud,
     . RAS,K2,lPu,lPd,RlPM,IlPM
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION phIF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      DOUBLE PRECISION Rm3,Im3,RAudt,IAudt,RxS,IxS,Rmsp,Imsp,Ast,IAst
      DOUBLE PRECISION RlPMt,IlPMt,RlM,IlM,RAqs,IAqs,Rltqs,Iltqs
      DOUBLE PRECISION MHuS,MHdS,MSS
      DOUBLE PRECISION IAL,IAK,IXIS
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/STSBSCALE/QSTSB
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/GAUGINOPAR/mur,M1r,M2r,msi
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/QPAR/l,k,Alcos1,Akcos2,muq,NUQ
      COMMON/SQUHIMASSM/MH02,MHC2
      COMMON/EFFPOTPAR/lu,ld,l3,l4,Rel5,Iml5,Rel6,Iml6,Rel7,Iml7,RAud,
     . RAS,K2,lPu,lPd,RlPM,IlPM
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/Z3VAUX/phIF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      COMMON/Z3VPOT/Rm3,Im3,RAudt,IAudt,RxS,IxS,Rmsp,Imsp,Ast,IAst,
     . RlPMt,IlPMt,RlM,IlM,RAqs,IAqs,Rltqs,Iltqs
      COMMON/MH2TREE/MHuS,MHdS,MSS
      COMMON/IMALAK/IAL,IAK,IXIS

      PI=4d0*DATAN(1d0)

      Mbi2=M1r**2
      Mwi2=M2r**2
      Mhi2=mur**2
      Msi2=msi**2


c      I- Pure singlet corrections

      aux=4.d0*(2d0*l*mur)**2*dlog(Mhi2/QSTSB)
     . +(4d0*k*(mupsi*DDCOS(phi02-phip)+ks2si))**2*dlog(Msi2/QSTSB)
     . -4d0*k*mupsi*DDCOS(Phi02-phiP)*l/muq*Msi2*(dlog(Msi2/QSTSB)-1d0)

      MH02(3,3)=MH02(3,3)-1d0/32d0/Pi**2*aux

      aux=(4d0*k*mupsi*DDSIN(phi02-phip))**2*dlog(Msi2/QSTSB)
     . -4d0*k*mupsi*DDCOS(Phi02-phiP)*l/muq*Msi2*(dlog(Msi2/QSTSB)-1d0)

      MH02(5,5)=MH02(5,5)-1d0/32d0/Pi**2*aux

       IF(k.ne.0d0)then
      aux=16d0*k**2*mupsi*DDSIN(Phi02-phiP)
     .   *(mupsi*DDCOS(phi02-phip)+ks2si)*dlog(Msi2/QSTSB)
     . +4d0*k*mupsi*DDSIN(Phi02-phiP)*l/muq*Msi2*(dlog(Msi2/QSTSB)-1d0)
       else
      aux=0d0
       endif

      MH02(3,5)=MH02(3,5)-1d0/32d0/Pi**2*aux
      MH02(5,3)=MH02(5,3)-1d0/32d0/Pi**2*aux

c        Minimization conditions

        MSS=MSS+(l*mur)**2*(dlog(Mhi2/QSTSB)-1d0)/4d0/Pi**2
     .    +k*(2d0*k+l*mupsi/muq)*Msi2*(dlog(Msi2/QSTSB)-1d0)/8d0/Pi**2

        aux=k*mupsi*DDSIN(phi02-phiP)/8d0/Pi**2
     .                              *Msi2*(dlog(Msi2/QSTSB)-1d0)
        IF(k.ne.0d0)THEN
         IAk=IAk+(l/muq)**2/k*aux
        ELSE
         IXIS=IXIS+aux*(muq/l)**2
        ENDIF


c      II- Corrections to the O(v^2) parameters of the potential

      RAudchi=2d0*l/32d0/Pi**2*(
     . g1q*M1r*DDCOS(PhiM1+Phi01)*Fsf1(Mbi2,Mhi2,QSTSB)
     . +3.d0*g2q*M2r*DDCOS(PhiM2+Phi01)*Fsf1(Mwi2,Mhi2,QSTSB))

      RlPMchi=-4.d0*ks2si*l**4/32d0/Pi**2/muq
     .         *DDCOS(Phi0)*Fsf1(Mhi2,Msi2,QSTSB)


      IlPMchi=-4.d0*ks2si*l**4/32d0/Pi**2/muq
     .         *DDSIN(Phi0)*Fsf1(Mhi2,Msi2,QSTSB)

c      differentiating coefficients wrt singlet fields

      aux=(2d0*l**2+g1q+3.d0*g2q)*(2d0*l*mur)*dlog(Mhi2/QSTSB)
     . +2d0*l**2*(2d0*ks2si**2)*dlog(Msi2/QSTSB)*l/mur
     . +g1q*(2d0*l*mur)*(fSF1(Mbi2,Mhi2,QSTSB)
     .                              +(Mbi2+Mhi2)*fSF3(Mbi2,Mhi2))
     . +3.d0*g2q*(2d0*l*mur)*(fSF1(Mwi2,Mhi2,QSTSB)
     .                              +(Mwi2+Mhi2)*fSF3(Mwi2,Mhi2))
     . +2d0*l**2*(2d0*l*mur)*(fSF1(Msi2,Mhi2,QSTSB)
     .                              +(Msi2+Mhi2)*fSF3(Msi2,Mhi2))
     . +2d0*l**2*(2d0*ks2si**2)*l/mur*(fSF1(Msi2,Mhi2,QSTSB)
     .                              +(Msi2+Mhi2)*fSF3(Mhi2,Msi2))

      dM2dHs=-1d0/32d0/Pi**2*aux*mur/muq

      aux=g1q*M1r*DDCOS(PhiM1+Phi01)*(2d0*l*mur)*fSF3(Mbi2,Mhi2)
     . +3.d0*g2q*M2r*DDCOS(PhiM2+Phi01)*(2d0*l*mur)*fSF3(Mwi2,Mhi2)

      RdAuddHs=2d0/32d0/Pi**2*l*aux

      aux=g1q*M1r*DDSIN(PhiM1+Phi01)*(2d0*l*mur)*fSF3(Mbi2,Mhi2)
     . +3.d0*g2q*M2r*DDSIN(PhiM2+Phi01)*(2d0*l*mur)*fSF3(Mwi2,Mhi2)

      IdAuddHs=2d0/32d0/Pi**2*l*aux

      aux=(2d0*l*mur)*fSF3(Msi2,Mhi2)
     . +(4.d0*k*msi)*fSF3(Mhi2,Msi2)

      RdlPMdHs=-8.d0/32d0/Pi**2*k*l**3*DDCOS(Phi01-Phi02)*aux
      IdlPMdHs=-8.d0/32d0/Pi**2*k*l**3*DDSIN(Phi01-Phi02)*aux

      aux=Max(dabs(XIS),dabs(MSP),dabs(XIF),dabs(MUP),dabs(M3H))
      IF(aux.ge.1d-4)THEN

      RlPMchi=g1q*M1r*DDCOS(PhiM1+Phi01)*Mhi2**2*FSF4(Mbi2,Mhi2)
     . +3d0*g2q*M2r*DDCOS(PhiM2+Phi01)*Mhi2**2*FSF4(Mwi2,Mhi2)
     . -4d0*l**2*(ks2si*DDCOS(Phi0)*FSF1(Msi2,Mhi2,QSTSB)
     .   +Mhi2*(mupsi*DDCOS(Phi01-phiP)+2d0*ks2si*DDCOS(Phi0))
     .                                        *FSF3(Msi2,Mhi2)
     .   +ks2si*(DDCOS(Phi0)*(2d0*Msi2-mupsi**2)-ks2si*mupsi
     .                      *DDCOS(Phi01-phiP))*FSF3(Mhi2,Msi2)
     .   +Mhi2**2*(mupsi*DDCOS(Phi01-phiP)+ks2si*DDCOS(Phi0))/2d0
     .                                        *FSF4(Msi2,Mhi2)
     .   +ks2si**2*(mupsi*DDCOS(Phi0+phiP-Phi02)+ks2si*DDCOS(Phi0))
     .                               *Msi2/2d0*FSF4(Mhi2,Msi2)
     .   +ks2si*Mhi2*Msi2*DDCOS(phi0)*FSF7(Msi2,Mhi2))
      RlPMchi=RlPMchi*l**2/muq/32d0/Pi**2

      IlPMchi=g1q*M1r*DDSIN(PhiM1+Phi01)*Mhi2**2*FSF4(Mbi2,Mhi2)
     . +3d0*g2q*M2r*DDSIN(PhiM2+Phi01)*Mhi2**2*FSF4(Mwi2,Mhi2)
     . -4d0*l**2*(ks2si*DDSIN(Phi0)*FSF1(Msi2,Mhi2,QSTSB)
     .   +Mhi2*(mupsi*DDSIN(Phi01-phiP)+2d0*ks2si*DDSIN(Phi0))
     .                                        *FSF3(Msi2,Mhi2)
     .   +ks2si*(DDSIN(Phi0)*(2d0*Msi2-mupsi**2)-ks2si*mupsi
     .                      *DDSIN(Phi01-phiP))*FSF3(Mhi2,Msi2)
     .   +Mhi2**2*(mupsi*DDSIN(Phi01-phiP)+ks2si*DDSIN(Phi0))/2d0
     .                                        *FSF4(Msi2,Mhi2)
     .   +ks2si**2*(mupsi*DDSIN(Phi0+phiP-Phi02)+ks2si*DDSIN(Phi0))
     .                               *Msi2/2d0*FSF4(Mhi2,Msi2)
     .   +ks2si*Mhi2*Msi2*DDSIN(phi0)*FSF7(Msi2,Mhi2))
      IlPMchi=IlPMchi*l**2/muq/32d0/Pi**2

      RlPMtchi=g1q*M1r*DDCOS(PhiM1+Phi01)*Mhi2
     .                *(2d0*FSF3(Mbi2,Mhi2)+Mhi2*FSF4(Mbi2,Mhi2))
     . +3d0*g2q*M2r*DDCOS(PhiM2+Phi01)*Mhi2
     .                *(2d0*FSF3(Mwi2,Mhi2)+Mhi2*FSF4(Mwi2,Mhi2))
     . -2d0*l**2*(Mhi2**2*(mupsi*DDCOS(Phi01-phiP)+ks2si*DDCOS(Phi0))
     .                                        *FSF4(Msi2,Mhi2)
     .   +ks2si**2*(mupsi**3*DDCOS(Phi0+3d0*(Phi02-phiP))
     .    +3d0*ks2si*mupsi**2*DDCOS(Phi0+2d0*(Phi02-phiP))
     .    +3d0*ks2si**2*mupsi*DDCOS(Phi01-phiP)+ks2si**3*DDCOS(Phi0))
     .                                        *FSF4(Mhi2,Msi2)
     .   +2d0*ks2si*Mhi2*(mupsi**2*DDCOS(Phi0+2d0*(Phi02-phiP))
     .                    +2d0*mupsi*ks2si*DDCOS(Phi01-phiP)
     .                    +ks2si**2*DDCOS(phi0))*FSF7(Msi2,Mhi2))
      RlPMtchi=RlPMtchi*l**2/muq/32d0/Pi**2

      IlPMtchi=g1q*M1r*DDSIN(PhiM1+Phi01)*Mhi2
     .                *(2d0*FSF3(Mbi2,Mhi2)+Mhi2*FSF4(Mbi2,Mhi2))
     . +3d0*g2q*M2r*DDSIN(PhiM2+Phi01)*Mhi2
     .                *(2d0*FSF3(Mwi2,Mhi2)+Mhi2*FSF4(Mwi2,Mhi2))
     . -2d0*l**2*(Mhi2**2*(mupsi*DDSIN(Phi01-phiP)+ks2si*DDSIN(Phi0))
     .                                        *FSF4(Msi2,Mhi2)
     .   +ks2si**2*(mupsi**3*DDSIN(Phi0+3d0*(Phi02-phiP))
     .    +3d0*ks2si*mupsi**2*DDSIN(Phi0+2d0*(Phi02-phiP))
     .    +3d0*ks2si**2*mupsi*DDSIN(Phi01-phiP)+ks2si**3*DDSIN(Phi0))
     .                                        *FSF4(Mhi2,Msi2)
     .   +2d0*ks2si*Mhi2*(mupsi**2*DDSIN(Phi0+2d0*(Phi02-phiP))
     .                    +2d0*mupsi*ks2si*DDSIN(Phi01-phiP)
     .                    +ks2si**2*DDSIN(phi0))*FSF7(Msi2,Mhi2))
      IlPMtchi=IlPMtchi*l**2/muq/32d0/Pi**2

      RlMchi=g1q*M1r*DDCOS(PhiM1+Phi01)*Mhi2
     .                *(2d0*FSF3(Mbi2,Mhi2)+Mhi2*FSF4(Mbi2,Mhi2))
     . +3d0*g2q*M2r*DDCOS(PhiM2+Phi01)*Mhi2
     .                *(2d0*FSF3(Mwi2,Mhi2)+Mhi2*FSF4(Mwi2,Mhi2))
     . -2d0*l**2*(Mhi2*(2d0*mupsi*DDCOS(Phi01-phiP)+3d0*ks2si
     .                            *DDCOS(Phi0))*FSF3(Msi2,Mhi2)
     .   +ks2si*(mupsi**2*DDCOS(Phi0+2d0*(Phi02-phiP))
     .    +4d0*ks2si*mupsi*DDCOS(Phi01-phiP)+3d0*ks2si**2*DDCOS(Phi0))
     .                                        *FSF3(Mhi2,Msi2)
     .   +(mupsi*DDCOS(Phi01-phiP)+ks2si*DDCOS(Phi0))
     .      *(Mhi2**2*FSF4(Msi2,Mhi2)+ks2si**2*Msi2*FSF4(Mhi2,Msi2)
     .       +2d0*ks2si*Mhi2*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .                                             *FSF7(Msi2,Mhi2)))
      RlMchi=2d0*RlMchi*l**2/muq/32d0/Pi**2

      IlMchi=g1q*M1r*DDSIN(PhiM1+Phi01)*Mhi2
     .                *(2d0*FSF3(Mbi2,Mhi2)+Mhi2*FSF4(Mbi2,Mhi2))
     . +3d0*g2q*M2r*DDSIN(PhiM2+Phi01)*Mhi2
     .                *(2d0*FSF3(Mwi2,Mhi2)+Mhi2*FSF4(Mwi2,Mhi2))
     . -2d0*l**2*(Mhi2*(2d0*mupsi*DDSIN(Phi01-phiP)+3d0*ks2si
     .                            *DDSIN(Phi0))*FSF3(Msi2,Mhi2)
     .   +ks2si*(mupsi**2*DDSIN(Phi0+2d0*(Phi02-phiP))
     .    +4d0*ks2si*mupsi*DDSIN(Phi01-phiP)+3d0*ks2si**2*DDSIN(Phi0))
     .                                        *FSF3(Mhi2,Msi2)
     .   +(mupsi*DDSIN(Phi01-phiP)+ks2si*DDSIN(Phi0))
     .      *(Mhi2**2*FSF4(Msi2,Mhi2)+ks2si**2*Msi2*FSF4(Mhi2,Msi2)
     .       +2d0*ks2si*Mhi2*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .                                             *FSF7(Msi2,Mhi2)))
      IlMchi=2d0*IlMchi*l**2/muq/32d0/Pi**2

      RAudchi=g1q*M1r*DDCOS(PhiM1+Phi01)*(FSF1(Mbi2,Mhi2,QSTSB)
     .     -3d0*Mhi2*FSF3(Mbi2,Mhi2)-2d0*Mhi2**2*FSF4(Mbi2,Mhi2))
     . +3d0*g2q*M2r*DDCOS(PhiM2+Phi01)*(FSF1(Mwi2,Mhi2,QSTSB)
     .     -3d0*Mhi2*FSF3(Mwi2,Mhi2)-2d0*Mhi2**2*FSF4(Mwi2,Mhi2))
     . +2d0*l**2*(Mhi2*(mupsi*DDCOS(Phi01-phiP)+2d0*ks2si
     .                            *DDCOS(Phi0))*FSF3(Msi2,Mhi2)
     .   +2d0*ks2si**2*(mupsi*DDCOS(Phi01-phiP)+ks2si*DDCOS(Phi0))
     .                                        *FSF3(Mhi2,Msi2)
     .   +2d0*Mhi2**2*(mupsi*DDCOS(Phi01-phiP)+ks2si*DDCOS(Phi0))
     .                                        *FSF4(Msi2,Mhi2)
     .   +2d0*ks2si**2*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .     *(mupsi**2*DDCOS(Phi0+2d0*(Phi02-phiP))
     .      +2d0*mupsi*ks2si*DDCOS(Phi01-phiP)+ks2si**2*DDCOS(Phi0))
     .                                        *FSF4(Mhi2,Msi2)
     .   +ks2si*Mhi2*(3d0*(mupsi**2*DDCOS(Phi0+2d0*(Phi02-phiP))
     .      +2d0*mupsi*ks2si*DDCOS(Phi01-phiP)+ks2si**2*DDCOS(Phi0))
     .                        +Msi2*DDCOS(Phi0))*FSF7(Msi2,Mhi2))
      RAudchi=2d0*RAudchi*l/32d0/Pi**2

      RAudtchi=-g1q*M1r*DDCOS(PhiM1+Phi01)*Mhi2
     .                *(FSF3(Mbi2,Mhi2)+2d0*Mhi2*FSF4(Mbi2,Mhi2))
     . -3d0*g2q*M2r*DDCOS(PhiM2+Phi01)*Mhi2
     .                *(FSF3(Mwi2,Mhi2)+2d0*Mhi2*FSF4(Mwi2,Mhi2))
     . -2d0*l**2*(mupsi*DDCOS(Phi01-phiP)*FSF1(Msi2,Mhi2,QSTSB)
     .   -3d0*Mhi2*(mupsi*DDCOS(Phi01-phiP)+2d0*ks2si*DDCOS(Phi0))
     .                                        *FSF3(Msi2,Mhi2)
     .   -2d0*ks2si*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .         *(mupsi*DDCOS(Phi01-phiP)+3d0*ks2si*DDCOS(Phi0))
     .                                        *FSF3(Mhi2,Msi2)
     .   -2d0*Mhi2**2*(mupsi*DDCOS(Phi01-phiP)+ks2si*DDCOS(Phi0))
     .                                        *FSF4(Msi2,Mhi2)
     .   -2d0*ks2si**2*(mupsi*DDCOS(Phi02-phiP)+ks2si)*Msi2*DDCOS(Phi0)
     .                                        *FSF4(Mhi2,Msi2)
     .   -ks2si*Mhi2*(mupsi**2*DDCOS(Phi0+2d0*(Phi02-phiP))
     .      +2d0*mupsi*ks2si*DDCOS(Phi01-phiP)+ks2si**2*DDCOS(Phi0)
     .                    +3d0*Msi2*DDCOS(Phi0))*FSF7(Msi2,Mhi2))
      RAudtchi=2d0*RAudtchi*l/32d0/Pi**2

      IAudtchi=-g1q*M1r*DDSIN(PhiM1+Phi01)*Mhi2
     .                *(FSF3(Mbi2,Mhi2)+2d0*Mhi2*FSF4(Mbi2,Mhi2))
     . -3d0*g2q*M2r*DDSIN(PhiM2+Phi01)*Mhi2
     .                *(FSF3(Mwi2,Mhi2)+2d0*Mhi2*FSF4(Mwi2,Mhi2))
     . -2d0*l**2*(mupsi*DDSIN(Phi01-phiP)*FSF1(Msi2,Mhi2,QSTSB)
     .   -3d0*Mhi2*(mupsi*DDSIN(Phi01-phiP)+2d0*ks2si*DDSIN(Phi0))
     .                                        *FSF3(Msi2,Mhi2)
     .   -2d0*ks2si*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .         *(mupsi*DDSIN(Phi01-phiP)+3d0*ks2si*DDSIN(Phi0))
     .                                        *FSF3(Mhi2,Msi2)
     .   -2d0*Mhi2**2*(mupsi*DDSIN(Phi01-phiP)+ks2si*DDSIN(Phi0))
     .                                        *FSF4(Msi2,Mhi2)
     .   -2d0*ks2si**2*(mupsi*DDCOS(Phi02-phiP)+ks2si)*Msi2*DDSIN(Phi0)
     .                                        *FSF4(Mhi2,Msi2)
     .   -ks2si*Mhi2*(mupsi**2*DDSIN(Phi0+2d0*(Phi02-phiP))
     .      +2d0*mupsi*ks2si*DDSIN(Phi01-phiP)+ks2si**2*DDSIN(Phi0)
     .                    +3d0*Msi2*DDSIN(Phi0))*FSF7(Msi2,Mhi2))
      IAudtchi=2d0*IAudtchi*l/32d0/Pi**2

      Rm3chi=g1q*M1r*DDCOS(PhiM1+Phi01)*Mhi2
     .                *(FSF3(Mbi2,Mhi2)+2d0*Mhi2*FSF4(Mbi2,Mhi2))
     . +3d0*g2q*M2r*DDCOS(PhiM2+Phi01)*Mhi2
     .                *(FSF3(Mwi2,Mhi2)+2d0*Mhi2*FSF4(Mwi2,Mhi2))
     . -2d0*l**2*(Mhi2*(mupsi*DDCOS(Phi01-phiP)+3d0*ks2si*DDCOS(Phi0))
     .                                          *FSF3(Msi2,Mhi2)
     .   +ks2si**2*(2d0*mupsi*DDCOS(Phi01-phiP)
     .    +mupsi*DDCOS(Phi0+phiP-Phi02)+3d0*ks2si*DDCOS(Phi0))
     .                                          *FSF3(Mhi2,Msi2)
     .   +2d0*(mupsi*DDCOS(Phi01-phiP)+ks2si*DDCOS(Phi0))
     .    *(Mhi2**2*FSF4(Msi2,Mhi2)
     .     +ks2si**2*(mupsi*DDCOS(Phi02-phiP)+ks2si)**2*FSF4(Mhi2,Msi2)
     .       +2d0*ks2si*Mhi2*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .                                             *FSF7(Msi2,Mhi2)))
      Rm3chi=2d0*Rm3chi*muq/32d0/Pi**2

      Im3chi=g1q*M1r*DDSIN(PhiM1+Phi01)*Mhi2
     .                *(FSF3(Mbi2,Mhi2)+2d0*Mhi2*FSF4(Mbi2,Mhi2))
     . +3d0*g2q*M2r*DDSIN(PhiM2+Phi01)*Mhi2
     .                *(FSF3(Mwi2,Mhi2)+2d0*Mhi2*FSF4(Mwi2,Mhi2))
     . -2d0*l**2*(Mhi2*(mupsi*DDSIN(Phi01-phiP)+3d0*ks2si*DDSIN(Phi0))
     .                                          *FSF3(Msi2,Mhi2)
     .   +ks2si**2*(2d0*mupsi*DDSIN(Phi01-phiP)
     .    +mupsi*DDSIN(Phi0+phiP-Phi02)+3d0*ks2si*DDSIN(Phi0))
     .                                          *FSF3(Mhi2,Msi2)
     .   +2d0*(mupsi*DDSIN(Phi01-phiP)+ks2si*DDSIN(Phi0))
     .    *(Mhi2**2*FSF4(Msi2,Mhi2)
     .     +ks2si**2*(mupsi*DDCOS(Phi02-phiP)+ks2si)**2*FSF4(Mhi2,Msi2)
     .       +2d0*ks2si*Mhi2*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .                                             *FSF7(Msi2,Mhi2)))
      Im3chi=2d0*Im3chi*muq/32d0/Pi**2

      RAqschi=-2d0*(2d0*l**2+g1q+3d0*g2q)*Mhi2
     .  +2d0*l**2*ks2si*(mupsi*DDCOS(Phi02-phiP)*dlog(Msi2/QSTSB)
     .    -2d0*ks2si*(mupsi*DDCOS(Phi02-phiP)+ks2si)**2/Msi2)
     .  -2d0*g1q*Mhi2**2*(2d0*FSF3(Mbi2,Mhi2)
     .                             +(Mbi2+Mhi2)*FSF4(Mbi2,Mhi2))
     .  -6d0*g2q*Mhi2**2*(2d0*FSF3(Mwi2,Mhi2)
     .                             +(Mwi2+Mhi2)*FSF4(Mwi2,Mhi2))
     .  +2d0*l**2*(ks2si*mupsi*DDCOS(Phi02-phiP)
     .           *(FSF1(Msi2,Mhi2,QSTSB)+(Msi2+Mhi2)*FSF3(Mhi2,Msi2))
     .   -2d0*ks2si**2*(mupsi*DDCOS(Phi02-phiP)+ks2si)**2
     .           *(2d0*FSF3(Mhi2,Msi2)
     . +(Msi2+Mhi2)*FSF4(Mhi2,Msi2))
     .   -2d0*Mhi2**2*(2d0*FSF3(Msi2,Mhi2)+(Msi2+Mhi2)*FSF4(Msi2,Mhi2))
     .   -4d0*ks2si*Mhi2*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .           *(FSF3(Mhi2,Msi2)+FSF3(Msi2,Mhi2)
     .                            +(Msi2+Mhi2)*FSF7(Msi2,Mhi2)))
      RAqschi=-RAqschi*l/muq/32d0/Pi**2

      IAqschi=2d0*l**2*ks2si*mupsi*DDSIN(Phi02-phiP)*(dlog(Msi2/QSTSB)
     .    -2d0*ks2si*(mupsi*DDCOS(Phi02-phiP)+ks2si)/Msi2)
     .  +2d0*l**2*mupsi*DDSIN(Phi02-phiP)*(ks2si
     .           *(FSF1(Msi2,Mhi2,QSTSB)+(Msi2+Mhi2)*FSF3(Mhi2,Msi2))
     .   -2d0*ks2si**2*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .           *(2d0*FSF3(Mhi2,Msi2)+(Msi2+Mhi2)*FSF4(Mhi2,Msi2))
     .   -2d0*ks2si*Mhi2*(FSF3(Mhi2,Msi2)+FSF3(Msi2,Mhi2)
     .                            +(Msi2+Mhi2)*FSF7(Msi2,Mhi2)))
      IAqschi=-IAqschi*l/muq/32d0/Pi**2

      Rltqschi=(2d0*l**2+g1q+3d0*g2q)*Mhi2/2d0
     .  +l**2*ks2si**2*(mupsi**2*DDCOS(2d0*(Phi02-phiP))
     .          +2d0*ks2si*mupsi*DDCOS(Phi02-phiP)+ks2si**2)/Msi2
     .  +g1q*Mhi2**2/2d0*(2d0*FSF3(Mbi2,Mhi2)
     .                             +(Mbi2+Mhi2)*FSF4(Mbi2,Mhi2))
     .  +3d0*g2q*Mhi2**2/2d0*(2d0*FSF3(Mwi2,Mhi2)
     .                             +(Mwi2+Mhi2)*FSF4(Mwi2,Mhi2))
     .  +l**2*(ks2si**2*(mupsi**2*DDCOS(2d0*(Phi02-phiP))
     .         +2d0*ks2si*mupsi*DDCOS(Phi02-phiP)+ks2si**2)
     .           *(2d0*FSF3(Mhi2,Msi2)+(Msi2+Mhi2)*FSF4(Mhi2,Msi2))
     .   +Mhi2**2*(2d0*FSF3(Msi2,Mhi2)+(Msi2+Mhi2)*FSF4(Msi2,Mhi2))
     .   +2d0*ks2si*Mhi2*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .           *(FSF3(Mhi2,Msi2)+FSF3(Msi2,Mhi2)
     .                            +(Msi2+Mhi2)*FSF7(Msi2,Mhi2)))
      Rltqschi=-Rltqschi*(l/muq)**2/32d0/Pi**2

      Iltqschi=2d0*l**2*ks2si**2/Msi2*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .  +l**2*(2d0*ks2si**2*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .           *(2d0*FSF3(Mhi2,Msi2)+(Msi2+Mhi2)*FSF4(Mhi2,Msi2))
     .   +2d0*ks2si*Mhi2*(FSF3(Mhi2,Msi2)+FSF3(Msi2,Mhi2)
     .                            +(Msi2+Mhi2)*FSF7(Msi2,Mhi2)))
      Iltqschi=-Iltqschi*(l/muq)**2/32d0/Pi**2
     .                          *mupsi*DDSIN(Phi02-phiP)

      lPqschi=(2d0*l**2+g1q+3d0*g2q)*Mhi2*(dlog(Mhi2/QSTSB)+1d0)
     .  +2d0*l**2*ks2si**2*(dlog(Msi2/QSTSB)+1d0)
     .  +g1q*Mhi2*(FSF1(Mbi2,Mhi2,QSTSB)
     .              +(Mbi2+3d0*Mhi2)*FSF3(Mbi2,Mhi2)
     .              +Mhi2*(Mbi2+Mhi2)*FSF4(Mbi2,Mhi2))
     .  +3d0*g2q*Mhi2*(FSF1(Mwi2,Mhi2,QSTSB)
     .              +(Mwi2+3d0*Mhi2)*FSF3(Mwi2,Mhi2)
     .              +Mhi2*(Mwi2+Mhi2)*FSF4(Mwi2,Mhi2))
     .  +2d0*l**2*((Mhi2+ks2si**2)*FSF1(Msi2,Mhi2,QSTSB)
     .    +Mhi2*(Msi2+3d0*Mhi2+2d0*ks2si
     .          *(mupsi*DDCOS(Phi02-phiP)+ks2si))*FSF3(Msi2,Mhi2)
     .    +(ks2si**2*(3d0*Msi2+Mhi2)+2d0*Mhi2*ks2si
     .          *(mupsi*DDCOS(Phi02-phiP)+ks2si))*FSF3(Mhi2,Msi2)
     .   +(Msi2+Mhi2)*(Mhi2**2*FSF4(Msi2,Mhi2)
     .     +ks2si**2*Msi2*FSF4(Mhi2,Msi2)
     .     +2d0*ks2si*Mhi2*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .                                           *FSF7(Msi2,Mhi2)))
      lPqschi=-lPqschi*(l/muq)**2/32d0/Pi**2

      ENDIF


c      III- Corrections to the O(v^4) parameters of the potential

      ludchi=-1d0/32d0/Pi**2*(
     . g1q**2/2d0*dlog(Mbi2/QSTSB)
     .  +5.d0*g2q**2/2d0*dlog(Mwi2/QSTSB)
     .  +(2d0*l**4+(5.d0*g2q**2+g1q**2+2d0*g1q*g2q)/2d0)
     .                              *dlog(Mhi2/QSTSB)
     .  +2d0*l**4*dlog(Msi2/QSTSB)
     . +g1q*g2q*Fsf1(Mbi2,Mwi2,QSTSB)
     . +g1q**2*(Mbi2+Mhi2)*Fsf3(Mhi2,Mbi2)
     . +g1q*(2d0*l**2*Mhi2+(g1q+g2q)*Mbi2)*Fsf3(Mbi2,Mhi2)
     . +g2q**2*(Mwi2+5.d0*Mhi2)*Fsf3(Mhi2,Mwi2)
     . +g2q*(2d0*l**2*Mhi2+(g1q+5.d0*g2q)*Mwi2)*Fsf3(Mwi2,Mhi2)
     . +4.d0*l**4*(Mhi2+Msi2)*Fsf3(Mhi2,Msi2)
     . +2d0*l**2*(8.d0*k**2+g1q+g2q)*Mhi2*Fsf3(Msi2,Mhi2)
     . +2d0*g1q*g2q*(M1r*M2r*DDCOS(PhiM1-PhiM2)+Mhi2)
     .                                   *Fsf5(Mhi2,Mbi2,Mwi2)
     . +g1q**2/2d0*(Mbi2+Mhi2)**2*Fsf7(Mbi2,Mhi2)
     . +g2q**2/2d0*(5.d0*Mwi2**2+5.d0*Mhi2**2+2d0*Mwi2*Mhi2)
     .                                   *Fsf7(Mwi2,Mhi2)
     . +2d0*l**4*(Mhi2+Msi2)**2*Fsf7(Msi2,Mhi2)
     . +g1q*g2q*(Mbi2*Mwi2+Mhi2**2+2d0*M1r*M2r*Mhi2
     .                *DDCOS(PhiM1-PhiM2))*Fsf6(Mbi2,Mwi2,Mhi2)
     . +2d0*g1q*l**2*Mhi2*(Mbi2+Msi2+2d0*M1r*(mupsi*DDCOS(PhiM1+phiP)
     .        +ks2si*DDCOS(PhiM1+Phi02)))*Fsf6(Mbi2,Msi2,Mhi2)
     . +2d0*g2q*l**2*Mhi2*(Mwi2+Msi2+2d0*M2r*(mupsi*DDCOS(PhiM2+phiP)
     .        +ks2si*DDCOS(PhiM2+Phi02)))*Fsf6(Mwi2,Msi2,Mhi2))

      l34chi=-1d0/32d0/Pi**2*(
     . g1q**2/2d0*dlog(Mbi2/QSTSB)
     .  +g2q**2/2d0*dlog(Mwi2/QSTSB)
     .  +(2d0*l**4+(g2q+g1q)**2/2d0)*dlog(Mhi2/QSTSB)
     .  +2d0*l**4*dlog(Msi2/QSTSB)
     . +g1q*g2q*Fsf1(Mbi2,Mwi2,QSTSB)
     . +g1q**2*(Mbi2+Mhi2)*Fsf3(Mhi2,Mbi2)
     . +g1q*((g1q+g2q)*Mbi2-2d0*(l**2-g1q-g2q)*Mhi2)*Fsf3(Mbi2,Mhi2)
     . +g2q**2*(5.d0*Mwi2+Mhi2)*Fsf3(Mhi2,Mwi2)
     . +g2q*((g1q+g2q)*Mwi2-2d0*(l**2-g1q-3.d0*g2q)*Mhi2)
     .                                        *Fsf3(Mwi2,Mhi2)
     . +4.d0*l**4*(Mhi2+Msi2)*Fsf3(Mhi2,Msi2)
     . +2d0*l**2*(l**2*(2d0*Msi2+4.d0*Mhi2)-(g1q+g2q)*Mhi2)
     .                                             *Fsf3(Msi2,Mhi2)
     . +2d0*g1q*g2q*(M1r*M2r*DDCOS(PhiM1-PhiM2)+Mhi2)
     .                                        *Fsf5(Mhi2,Mbi2,Mwi2)
     . +g1q**2/2d0*(Mbi2**2+Mhi2**2+6.d0*Mbi2*Mhi2)*Fsf7(Mbi2,Mhi2)
     . +g2q**2/2d0*(Mwi2**2+Mhi2**2+22d0*Mwi2*Mhi2)
     .                                        *Fsf7(Mwi2,Mhi2)
     . +2d0*l**4*(Mhi2**2+Msi2**2+6.d0*Mhi2*Msi2)*Fsf7(Msi2,Mhi2)
     . +g1q*g2q*(Mbi2*Mwi2+Mhi2**2+2d0*Mhi2*
     .  (Mbi2+Mwi2+M1r*M2r*DDCOS(PhiM1-PhiM2)))*Fsf6(Mbi2,Mwi2,Mhi2)
     . -2d0*g1q*l**2*Mhi2*(Mbi2+Msi2+2d0*M1r*(mupsi*DDCOS(PhiM1+phiP)
     .        +ks2si*DDCOS(PhiM1+Phi02)))*Fsf6(Mbi2,Msi2,Mhi2)
     . -2d0*g2q*l**2*Mhi2*(Mwi2+Msi2+2d0*M2r*(mupsi*DDCOS(PhiM2+phiP)
     .        +ks2si*DDCOS(PhiM2+Phi02)))*Fsf6(Mwi2,Msi2,Mhi2))

      l4chi=-1d0/32d0/Pi**2*(
     . -2d0*g2q**2*(dlog(Mwi2/QSTSB)+dlog(Mhi2/QSTSB))
     . +2d0*g1q*g2q*(Fsf1(Mbi2,Mwi2,QSTSB)+dlog(Mhi2/QSTSB))
     . +2d0*g1q*(g2q*Mbi2+(g1q-g2q)*Mhi2)*Fsf3(Mbi2,Mhi2)
     . +4.d0*g2q**2*(Mwi2-Mhi2)*Fsf3(Mhi2,Mwi2)
     . +2d0*g2q*((g1q-2d0*g2q)*Mwi2-(2d0*l**2+g1q-3.d0*g2q)
     .                                     *Mhi2)*Fsf3(Mwi2,Mhi2)
     . +4.d0*l**2*(2d0*l**2-g2q)*Mhi2*Fsf3(Msi2,Mhi2)
     . +4.d0*g1q*g2q*(M1r*M2r*DDCOS(PhiM1-PhiM2)+Mhi2)
     .                                      *Fsf5(Mhi2,Mbi2,Mwi2)
     . +2d0*g1q**2*Mbi2*Mhi2*Fsf7(Mbi2,Mhi2)
     . +2d0*g2q**2*(-Mwi2**2-Mhi2**2+5.d0*Mwi2*Mhi2)*Fsf7(Mwi2,Mhi2)
     . +8.d0*l**4*Mhi2*Msi2*Fsf7(Msi2,Mhi2)
     . +2d0*g1q*g2q*(Mbi2*Mwi2+Mhi2**2-Mhi2*(Mbi2+Mwi2
     .      -2d0*M1r*M2r*DDCOS(PhiM1-PhiM2)))*Fsf6(Mbi2,Mwi2,Mhi2)
     . -4.d0*g2q*l**2*Mhi2*((Mwi2+Msi2)+2d0*M2r*(mupsi*
     . DDCOS(PhiM2+phiP)+ks2si*DDCOS(PhiM2+Phi02)))
     .*Fsf6(Mwi2,Msi2,Mhi2))

      Rel5chi=-2d0/32d0/Pi**2*Mhi2*(
     . g1q**2*Mbi2*DDCOS(2d0*(PhiM1+Phi01))*Fsf7(Mbi2,Mhi2)
     . +3.d0*g2q**2*Mwi2*DDCOS(2d0*(PhiM2+Phi01))*Fsf7(Mwi2,Mhi2)
     . +2d0*g1q*g2q*M1r*M2r*DDCOS(PhiM1+PhiM2+2d0*Phi01)
     .                                      *Fsf6(Mbi2,Mwi2,Mhi2)
     . +4.d0*l**4*(mupsi**2*DDCOS(2d0*(Phi01-Phip))
     .             +2d0*mupsi*ks2si*DDCOS(2d0*Phi01-phiP-Phi02)
     .             +ks2si**2*DDCOS(2d0*(Phi01-Phi02)))
     .*Fsf7(Msi2,Mhi2))

      Iml5chi=-2d0/32d0/Pi**2*Mhi2*(
     . g1q**2*Mbi2*DDSIN(2d0*(PhiM1+Phi01))*Fsf7(Mbi2,Mhi2)
     . +3.d0*g2q**2*Mwi2*DDSIN(2d0*(PhiM2+Phi01))*Fsf7(Mwi2,Mhi2)
     . +2d0*g1q*g2q*M1r*M2r*DDSIN(PhiM1+PhiM2+2d0*Phi01)
     .                            *Fsf6(Mbi2,Mwi2,Mhi2)
     . +4.d0*l**4*(mupsi**2*DDSIN(2d0*(Phi01-Phip))
     .             +2d0*mupsi*ks2si*DDSIN(2d0*Phi01-phiP-Phi02)
     .             +ks2si**2*DDSIN(2d0*(Phi01-Phi02)))
     .*Fsf7(Msi2,Mhi2))

      Rel67chi=-1d0/32d0/Pi**2*(
     . -g1q**2*M1r*mur*DDCOS(PhiM1+Phi01)*Fsf3(Mhi2,Mbi2)
     . -g1q*(g1q+g2q)*M1r*mur*DDCOS(PhiM1+Phi01)*Fsf3(Mbi2,Mhi2)
     . -3.d0*g2q**2*M2r*mur*DDCOS(PhiM2+Phi01)*Fsf3(Mhi2,Mwi2)
     . -g2q*(g1q+3.d0*g2q)*M2r*mur*DDCOS(PhiM2+Phi01)
     .                                         *Fsf3(Mwi2,Mhi2)
     . +4.d0*l**4*mur*(mupsi*DDCOS(Phi01-phip)+ks2si*DDCOS(Phi01-Phi02))
     .                       *(Fsf3(Mhi2,Msi2)+Fsf3(Msi2,Mhi2))
     . -g1q*g2q*mur*(M1r*DDCOS(PhiM1+Phi01)+M2r*DDCOS(PhiM2+Phi01))
     .                                    *Fsf5(Mhi2,Mbi2,Mwi2)
     . -g1q**2*M1r*mur*DDCOS(PhiM1+Phi01)*(Mbi2+Mhi2)
     .                                         *Fsf7(Mbi2,Mhi2)
     . -3.d0*g2q**2*M2r*mur*DDCOS(PhiM2+Phi01)*(Mwi2+Mhi2)
     .                                         *Fsf7(Mwi2,Mhi2)
     . +4.d0*l**4*mur*(mupsi*DDCOS(Phi01-phip)+ks2si*DDCOS(Phi01-Phi02))
     .                             *(Mhi2+Msi2)*Fsf7(Msi2,Mhi2)
     . -g1q*g2q*mur*(M1r*(Mwi2+Mhi2)*DDCOS(PhiM1+Phi01)
     .               +M2r*(Mbi2+Mhi2)*DDCOS(PhiM2+Phi01))
     .                                    *Fsf6(Mbi2,Mwi2,Mhi2))

      Iml67chi=-1d0/32d0/Pi**2*(
     . -g1q**2*M1r*mur*DDSIN(PhiM1+Phi01)*Fsf3(Mhi2,Mbi2)
     . -g1q*(g1q+g2q)*M1r*mur*DDSIN(PhiM1+Phi01)*Fsf3(Mbi2,Mhi2)
     . -3.d0*g2q**2*M2r*mur*DDSIN(PhiM2+Phi01)*Fsf3(Mhi2,Mwi2)
     . -g2q*(g1q+3.d0*g2q)*M2r*mur*DDSIN(PhiM2+Phi01)
     .                                         *Fsf3(Mwi2,Mhi2)
     . +4.d0*l**4*mur*(mupsi*DDSIN(Phi01-phip)+ks2si*DDSIN(Phi01-Phi02))
     .                       *(Fsf3(Mhi2,Msi2)+Fsf3(Msi2,Mhi2))
     . -g1q*g2q*mur*(M1r*DDSIN(PhiM1+Phi01)+M2r*DDSIN(PhiM2+Phi01))
     .                                    *Fsf5(Mhi2,Mbi2,Mwi2)
     . -g1q**2*M1r*mur*DDSIN(PhiM1+Phi01)*(Mbi2+Mhi2)
     .                                         *Fsf7(Mbi2,Mhi2)
     . -3.d0*g2q**2*M2r*mur*DDSIN(PhiM2+Phi01)*(Mwi2+Mhi2)
     .                                         *Fsf7(Mwi2,Mhi2)
     . +4.d0*l**4*mur*(mupsi*DDSIN(Phi01-phip)+ks2si*DDSIN(Phi01-Phi02))
     .                             *(Mhi2+Msi2)*Fsf7(Msi2,Mhi2)
     . -g1q*g2q*mur*(M1r*(Mwi2+Mhi2)*DDSIN(PhiM1+Phi01)
     .               +M2r*(Mbi2+Mhi2)*DDSIN(PhiM2+Phi01))
     .                                    *Fsf6(Mbi2,Mwi2,Mhi2))


      aux=Max(dabs(XIS),dabs(MSP),dabs(XIF),dabs(MUP),dabs(M3H))
      IF(aux.le.1d-4)THEN

c      IV- Corrections to the neutral Higgs mass-matrix

      aux=2d0*ludchi*vuq**2+((RAudchi+RlPMchi*muq/l)*muq/l
     .                     +Rel67chi*(vdq**2-3.d0*vuq**2))*vdq/vuq

      MH02(1,1)=MH02(1,1)+aux

      aux=2d0*ludchi*vdq**2+((RAudchi+RlPMchi*muq/l)*muq/l
     .                     +Rel67chi*(vuq**2-3.d0*vdq**2))*vuq/vdq

      MH02(2,2)=MH02(2,2)+aux

      aux=-((RAudchi+RlPMchi*muq/l)*muq/l+3.d0*Rel67chi
     .    *(vuq**2+vdq**2))+2d0*(l34chi+Rel5chi)*vuq*vdq

      MH02(1,2)=MH02(1,2)+aux
      MH02(2,1)=MH02(2,1)+aux

      aux=-(RAudchi+2d0*RlPMchi*muq/l)*vdq+dM2dHs*vuq
     . -(RdAuddHs+RdlPMdHs*muq/l)*vdq*muq/l+l*vdq/muq*
     . (2d0*Rel5chi*vuq*vdq-Rel67chi*(3.d0*vuq**2+vdq**2))

      MH02(1,3)=MH02(1,3)+aux
      MH02(3,1)=MH02(3,1)+aux

      aux=-(RAudchi+2d0*RlPMchi*muq/l)*vuq+dM2dHs*vdq
     . -(RdAuddHs+RdlPMdHs*muq/l)*vuq*muq/l+l*vuq/muq*
     . (2d0*Rel5chi*vuq*vdq-Rel67chi*(vuq**2+3.d0*vdq**2))

      MH02(2,3)=MH02(2,3)+aux
      MH02(3,2)=MH02(3,2)+aux

      MH02(3,3)=MH02(3,3)+l/muq*vuq*vdq*(
     . RAudchi-2d0*l/muq*RdAuddHs-4.d0*(l/muq)**2*RdlPMdHs)

      aux=(RAudchi+RlPMchi*muq/l)*muq/l-2d0*Rel5chi*vuq*vdq
     .     +Rel67chi*(vuq**2+vdq**2)

      MH02(4,4)=MH02(4,4)+aux*(vu**2+vd**2)/vuq/vdq

      aux=(RAudchi-2d0*RlPMchi*muq/l)
     . -(2d0*Rel5chi*vuq*vdq-Rel67chi*(vuq**2+vdq**2))*l/muq

      MH02(4,5)=MH02(4,5)+aux*dsqrt(vu**2+vd**2)
      MH02(5,4)=MH02(5,4)+aux*dsqrt(vu**2+vd**2)

      aux=(RAudchi+4.d0*RlPMchi*muq/l)*vuq*vdq*l/muq

      MH02(5,5)=MH02(5,5)+aux

      aux=(2d0*Iml67chi*vuq-Iml5chi*vdq)

      MH02(1,4)=MH02(1,4)+aux*dsqrt(vu**2+vd**2)
      MH02(4,1)=MH02(4,1)+aux*dsqrt(vu**2+vd**2)

      aux=(2d0*Iml67chi*vdq-Iml5chi*vuq)

      MH02(2,4)=MH02(2,4)+aux*dsqrt(vu**2+vd**2)
      MH02(4,2)=MH02(4,2)+aux*dsqrt(vu**2+vd**2)

      aux=-3.d0*IlPMchi*muq/l
     .      -(Iml5chi*vuq*vdq-2d0*Iml67chi*vuq**2)*l/muq

      MH02(1,5)=MH02(1,5)+aux*vdq
      MH02(5,1)=MH02(5,1)+aux*vdq

      aux=-3.d0*IlPMchi*muq/l
     .      -(Iml5chi*vuq*vdq-2d0*Iml67chi*vdq**2)*l/muq

      MH02(2,5)=MH02(2,5)+aux*vuq
      MH02(5,2)=MH02(5,2)+aux*vuq

      aux=IlPMchi*muq/l-Iml5chi*vuq*vdq*l/muq

      MH02(3,4)=MH02(3,4)+aux*dsqrt(vu**2+vd**2)
      MH02(4,3)=MH02(4,3)+aux*dsqrt(vu**2+vd**2)

       IF(k.ne.0d0)then
      aux=(4.d0*IlPMchi-Iml5chi*vuq*vdq/(muq/l)**2)*vuq*vdq
       else
      aux=-2d0*(IlPMchi+Iml5chi*vuq*vdq/(muq/l)**2)*vuq*vdq
       endif

      MH02(3,5)=MH02(3,5)+aux
      MH02(5,3)=MH02(5,3)+aux


c      V- Corrections to the charged-Higgs mass

      aux=(RAudchi+RlPMchi*muq/l)*muq/l-(l4chi+Rel5chi)*vuq*vdq
     . +Rel67chi*(vuq**2+vdq**2)

      MHC2=MHC2+aux*(vu**2+vd**2)/vuq/vdq


c      VI- Inclusion in the corrected parameters

      lu=lu+ludchi
      ld=ld+ludchi
      l3=l3+l34chi-l4chi
      l4=l4+l4chi
      Rel5=Rel5+Rel5chi
      Iml5=Iml5+Iml5chi
      Rel6=Rel6+Rel67chi
      Iml6=Iml6+Iml67chi
      Rel7=Rel7+Rel67chi
      Iml7=Iml7+Iml67chi
      RAud=RAud+RAudchi
      RlPM=RlPM+RlPMchi
      IlPM=IlPM+IlPMchi
      lPu=lPu+l*dM2dHs/(2d0*muq)
      lPd=lPd+l*dM2dHs/(2d0*muq)

c        minimization conditions

      aux=g1q*fSF0(Mbi2,QSTSB)+3d0*fSF0(Mwi2,QSTSB)
     .    +(2d0*l**2+g1q+3d0*g2q)*fSF0(Mhi2,QSTSB)
     .    +2d0*l**2*fSF0(Msi2,QSTSB)
     .    +g1q*(Mbi2+Mhi2)*fSF1(Mbi2,Mhi2,QSTSB)
     .    +3d0*g2q*(Mwi2+Mhi2)*fSF1(Mwi2,Mhi2,QSTSB)
     .    +2d0*l**2*(Msi2+Mhi2)*fSF1(Msi2,Mhi2,QSTSB)

        MHuS=MHuS+aux/32d0/Pi**2
     . +((RAudchi+RlPMchi*muq/l)*muq/l+Rel67chi*(3d0*vuq**2+vdq**2))
     .                                                   *vdq/vuq
     . -(ludchi*vuq**2+(l34chi+Rel5chi)*vdq**2)

        MHdS=MHdS+aux/32d0/Pi**2
     . +((RAudchi+RlPMchi*muq/l)*muq/l+Rel67chi*(3d0*vuq**2+vdq**2))
     .                                                   *vdq/vuq
     . -(ludchi*vuq**2+(l34chi+Rel5chi)*vdq**2)

        MSS=MSS
     . -l*dM2dHs/(2d0*muq)*(vuq**2+vdq**2)
     . +l*vuq*vdq/muq*(RAudchi+RdAuddHs*muq/l+2d0*RlPMchi*muq/l)

      aux=g1q*M1r*DDSIN(PhiM1+phi01)*fSF1(Mbi2,Mhi2,QSTSB)
     .   +3d0*g2q*M2r*DDSIN(phiM2+phi01)*fSF1(Mwi2,Mhi2,QSTSB)

        IAl=IAl-aux/16d0/Pi**2-IlPMchi*muq/l**2
     .     +(Iml5chi*vuq*vdq-Iml67chi*(vuq**2+vdq**2))/muq


      ELSE
c      with violation of Z3

c      IVbis- Corrections to the neutral Higgs mass-matrix

      aux=2d0*ludchi*vuq**2
     . +((RAudchi+RAudtchi+(RlPMchi+RlPMtchi+RlMchi)*muq/l)*muq/l
     .         +Rm3chi+Rel67chi*(vdq**2-3.d0*vuq**2))*vdq/vuq

      MH02(1,1)=MH02(1,1)+aux

      aux=2d0*ludchi*vdq**2
     . +((RAudchi+RAudtchi+(RlPMchi+RlPMtchi+RlMchi)*muq/l)*muq/l
     .          +Rm3chi+Rel67chi*(vuq**2-3.d0*vdq**2))*vuq/vdq

      MH02(2,2)=MH02(2,2)+aux

      aux=2d0*(l34chi+Rel5chi)*vuq*vdq
     . -((RAudchi+RAudtchi+(RlPMchi+RlPMtchi+RlMchi)*muq/l)*muq/l
     .    +Rm3chi+3.d0*Rel67chi*(vuq**2+vdq**2))

      MH02(1,2)=MH02(1,2)+aux
      MH02(2,1)=MH02(2,1)+aux

      aux=2d0*(RAqschi+(lPqschi+2d0*Rltqschi)*muq/l)*vuq
     . -(RAudchi+RAudtchi+2d0*(RlPMchi+RlPMtchi+RlMchi)*muq/l)*vdq

      MH02(1,3)=MH02(1,3)+aux
      MH02(3,1)=MH02(3,1)+aux

      aux=2d0*(RAqschi+(lPqschi+2d0*Rltqschi)*muq/l)*vdq
     . -(RAudchi+RAudtchi+2d0*(RlPMchi+RlPMtchi+RlMchi)*muq/l)*vuq

      MH02(2,3)=MH02(2,3)+aux
      MH02(3,2)=MH02(3,2)+aux

      MH02(3,3)=MH02(3,3)
     . +l/muq*(vuq*vdq*(RAudchi+RAudtchi)-RAqschi*(vuq**2+vdq**2))

      aux=(RAudchi+RAudtchi+(RlPMchi+RlPMtchi+RlMchi)*muq/l)*muq/l
     .    +Rm3chi-2d0*Rel5chi*vuq*vdq+Rel67chi*(vuq**2+vdq**2)

      MH02(4,4)=MH02(4,4)+aux*(vu**2+vd**2)/vuq/vdq

      aux=(RAudchi-RAudtchi+2d0*(RlPMtchi-RlPMchi)*muq/l)

      MH02(4,5)=MH02(4,5)+aux*dsqrt(vu**2+vd**2)
      MH02(5,4)=MH02(5,4)+aux*dsqrt(vu**2+vd**2)

      aux=((RAudchi+RAudtchi+4.d0*(RlPMtchi+RlPMchi)*muq/l)*vuq*vdq
     . -(RAqschi+4d0*Rltqschi*muq/l)*(vuq**2+vdq**2))*l/muq

      MH02(5,5)=MH02(5,5)+aux

      aux=(2d0*Iml67chi*vuq-Iml5chi*vdq)

      MH02(1,4)=MH02(1,4)+aux*dsqrt(vu**2+vd**2)
      MH02(4,1)=MH02(4,1)+aux*dsqrt(vu**2+vd**2)

      aux=(2d0*Iml67chi*vdq-Iml5chi*vuq)

      MH02(2,4)=MH02(2,4)+aux*dsqrt(vu**2+vd**2)
      MH02(4,2)=MH02(4,2)+aux*dsqrt(vu**2+vd**2)

      aux=-(IAqschi+2d0*Iltqschi*muq/l)*vuq
     .   -(Im3chi*l/muq
     .    +2d0*IAudtchi+(3.d0*IlPMchi-IlPMtchi+IlMchi)*muq/l)*vdq

      MH02(1,5)=MH02(1,5)+aux
      MH02(5,1)=MH02(5,1)+aux

      aux=-(IAqschi+2d0*Iltqschi*muq/l)*vdq
     .   -(Im3chi*l/muq
     .    +2d0*IAudtchi+(3.d0*IlPMchi-IlPMtchi+IlMchi)*muq/l)*vuq

      MH02(2,5)=MH02(2,5)+aux
      MH02(5,2)=MH02(5,2)+aux

      aux=(IlPMchi+IlPMtchi+IlMchi)*muq/l-Im3chi*l/muq

      MH02(3,4)=MH02(3,4)+aux*dsqrt(vu**2+vd**2)
      MH02(4,3)=MH02(4,3)+aux*dsqrt(vu**2+vd**2)

       IF(k.ne.0d0)then
      aux=2d0*(Im3chi*(l/muq)**2+2d0*IAudtchi*l/muq
     .                +(2d0*IlPMchi+IlMchi))*vuq*vdq
     . +2d0*(IAqschi*l/muq+Iltqschi)*(vuq**2+vdq**2)
       else
      aux=-2d0*(IlPMchi*vuq*vdq+Iltqschi*(vuq**2+vdq**2))
       endif

      MH02(3,5)=MH02(3,5)+aux
      MH02(5,3)=MH02(5,3)+aux


c      Vbis- Corrections to the charged-Higgs mass

      aux=(RAudchi+RAudtchi+(RlPMchi+RlPMtchi+RlMchi)*muq/l)*muq/l
     . +Rm3chi-(l4chi+Rel5chi)*vuq*vdq+Rel67chi*(vuq**2+vdq**2)

      MHC2=MHC2+aux*(vu**2+vd**2)/vuq/vdq                              !/(ZHu*ZHd)


c      VIbis- Inclusion in the corrected parameters

      lu=lu+ludchi
      ld=ld+ludchi
      l3=l3+l34chi-l4chi
      l4=l4+l4chi
      Rel5=Rel5+Rel5chi
      Iml5=Iml5+Iml5chi
      Rel6=Rel6+Rel67chi
      Iml6=Iml6+Iml67chi
      Rel7=Rel7+Rel67chi
      Iml7=Iml7+Iml67chi
      Rm3=Rm3+Rm3chi
      Im3=Im3+Im3chi
      RAud=RAud+RAudchi
      RAudt=RAudt+RAudtchi
      IAudt=IAudt+IAudtchi
      RlPM=RlPM+RlPMchi
      IlPM=IlPM+IlPMchi
      RlPMt=RlPMt+RlPMtchi
      IlPMt=IlPMt+IlPMtchi
      RlM=RlM+RlMchi
      IlM=IlM+IlMchi
      lPu=lPu+lPqschi
      lPd=lPd+lPqschi
      RAqs=RAqschi
      IAqs=IAqschi
      Rltqs=Rltqschi
      Iltqs=Iltqschi

c        minimization conditions

      aux=g1q*fSF0(Mbi2,QSTSB)+3d0*fSF0(Mwi2,QSTSB)
     .    +(2d0*l**2+g1q+3d0*g2q)*fSF0(Mhi2,QSTSB)
     .    +2d0*l**2*fSF0(Msi2,QSTSB)
     .    +g1q*(Mbi2+Mhi2)*fSF1(Mbi2,Mhi2,QSTSB)
     .    +3d0*g2q*(Mwi2+Mhi2)*fSF1(Mwi2,Mhi2,QSTSB)
     .    +2d0*l**2*(Msi2+Mhi2)*fSF1(Msi2,Mhi2,QSTSB)

        MHuS=MHuS+aux/32d0/Pi**2
     . +((RAudchi+RAudtchi+(RlPMchi+RlPMtchi+RlMchi)*muq/l)*muq/l
     .          +Rm3chi+Rel67chi*(3d0*vuq**2+vdq**2))*vdq/vuq
     . -(ludchi*vuq**2+(l34chi+Rel5chi)*vdq**2)

        MHdS=MHdS+aux/32d0/Pi**2
     . +((RAudchi+RAudtchi+(RlPMchi+RlPMtchi+RlMchi)*muq/l)*muq/l
     .          +Rm3chi+Rel67chi*(3d0*vuq**2+vdq**2))*vdq/vuq
     . -(ludchi*vuq**2+(l34chi+Rel5chi)*vdq**2)

        MSS=MSS
     .    -(RAqschi*l/muq+(lPqschi+2d0*Rltqschi))*(vuq**2+vdq**2)
     .    +l*vuq*vdq/muq
     .     *(RAudchi+RAudtchi+2d0*(RlPMchi+RlPMtchi+RlMchi)*muq/l)

      aux=g1q*M1r*DDSIN(PhiM1+phi01)*fSF1(Mbi2,Mhi2,QSTSB)
     .   +3d0*g2q*M2r*DDSIN(phiM2+phi01)*fSF1(Mwi2,Mhi2,QSTSB)

        IAl=IAl-aux/16d0/Pi**2-(IlPMchi+IlPMtchi+IlMchi)*muq/l**2
     .     -Im3chi/muq-IAudt/l
     .     +(Iml5chi*vuq*vdq-Iml67chi*(vuq**2+vdq**2))/muq

      aux=-l*aux/16d0/Pi**2-(IlPMchi+IlPMtchi+IlMchi)*muq/l
     .     -l*Im3chi/muq-IAudt
     .     +l*(Iml5chi*vuq*vdq-Iml67chi*(vuq**2+vdq**2))/muq
     . -(IAudtchi+2d0*(IlPMchi-IlPMtchi)*muq/l)*vuq*vdq*(l/muq)**2
     c -(IAqschi+2d0*Iltqschi*muq/l)*(vuq**2+vdq**2)*(l/muq)**2

        IF(k.ne.0d0)THEN
      IAk=IAk+aux/k
        ELSE
      IXIS=IXIS+aux*(muq/l)**2
        ENDIF

      ENDIF


c      print*,'MH02_1*',MH02(1,1),MH02(1,2),MH02(1,3),MH02(1,4),MH02(1,5)
c      print*,'MH02_2*',MH02(2,1),MH02(2,2),MH02(2,3),MH02(2,4),MH02(2,5)
c      print*,'MH02_3*',MH02(3,1),MH02(3,2),MH02(3,3),MH02(3,4),MH02(3,5)
c      print*,'MH02_4*',MH02(4,1),MH02(4,2),MH02(4,3),MH02(4,4)
c     c /(ZHu*ZHd),MH02(4,5)/dsqrt(Zs*ZHu*ZHd)
c      print*,'MH02_5*',MH02(5,1),MH02(5,2),MH02(5,3),MH02(5,4),MH02(5,5)
c      print*,'MHC2',MHC2     

      RETURN
      END

************************************************************************************************

      SUBROUTINE MHIGGSLOOP_GAUGEHIGGS_CPV()

c         One-loop corrections to the Higgs potential
c                 - EW gauge + Higgs contribution
c      - The gauge and Higgs contributions to the Higgs squared-mass matrices
c        are computed (from 0-momentum self - tadpole) in the Feynmann gauge
c        and added to the quantities stored in SQUHIMASSP.
c      - The corresponding contributions to the Z3-conserving parameters of the
c        effective Higgs potential are deduced from `inverting the system'
c        and added to the quantities stored in EFFPOTPAR.

      IMPLICIT NONE

      INTEGER I,J,M,N

      DOUBLE PRECISION fSF0,fSF1
      DOUBLE PRECISION Pi,aux,MW2,MZ2,MHT2(6),MHTC(2),XHT(6,6)
      DOUBLE PRECISION dMH0(5,5),dMHC,sinbq,cosbq,gRH0HpHm(5,2,2),
     . gIH0HpHm(5,2,2),gHHH(6,6,6),gHHHH(6,6,6),XHG(6,6)
      DOUBLE PRECISION PIS111(6,6,6),PIS122(6,6,6),PIS133(6,6,6),
     . PIS144(6,6,6),PIS155(6,6,6),PIS166(6,6,6),PIS211(6,6,6),
     . PIS222(6,6,6),PIS233(6,6,6),PIS244(6,6,6),PIS255(6,6,6),
     . PIS266(6,6,6),PIS311(6,6,6),PIS322(6,6,6),PIS333(6,6,6),
     . PIS344(6,6,6),PIS355(6,6,6),PIS366(6,6,6),PIS312(6,6,6),
     . PIS345(6,6,6),PIS156(6,6,6),PIS256(6,6,6),PIS346(6,6,6),
     . PIS356(6,6,6),PIS513(6,6,6),PIS423(6,6,6),PIS612(6,6,6),
     . PIS456(6,6,6),PIS433(6,6,6),PIS613(6,6,6),PIS466(6,6,6),
     . PIS533(6,6,6),PIS623(6,6,6),PIS566(6,6,6),PIS666(6,6,6),
     . PIS633(6,6,6),PIS246(6,6,6)
      DOUBLE PRECISION PISS1111(6,6,6),PISS1122(6,6,6),
     . PISS1133(6,6,6),PISS1144(6,6,6),PISS1155(6,6,6),
     . PISS1166(6,6,6),PISS2222(6,6,6),PISS2233(6,6,6),
     . PISS2244(6,6,6),PISS2255(6,6,6),PISS2266(6,6,6),
     . PISS3333(6,6,6),PISS3344(6,6,6),PISS3355(6,6,6),
     . PISS3366(6,6,6),PISS4444(6,6,6),PISS4455(6,6,6),
     . PISS4466(6,6,6),PISS5555(6,6,6),PISS5566(6,6,6),
     . PISS6666(6,6,6),PISS1233(6,6,6),PISS1266(6,6,6),
     . PISS1356(6,6,6),PISS1335(6,6,6),PISS1236(6,6,6),
     . PISS3345(6,6,6),PISS2346(6,6,6),PISS2334(6,6,6),
     . PISS4566(6,6,6),PISS1566(6,6,6),PISS2466(6,6,6),
     . PISS3456(6,6,6)
      DOUBLE PRECISION RASH,K2H,RAudH,RlPMH,IlPMH,lPuH,lPdH,
     .                             luH,ldH,l3H,l4H,Il5H,Il6H,Il7H
      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION l,k,Alcos1,Akcos2,muq,nuq
      DOUBLE PRECISION MH02(5,5),MHC2
      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      DOUBLE PRECISION phIF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      DOUBLE PRECISION lu,ld,l3,l4,Rel5,Iml5,Rel6,Iml6,Rel7,Iml7,RAud,
     . RAS,K2,lPu,lPd,RlPM,IlPM
      DOUBLE PRECISION MHuS,MHdS,MSS
      DOUBLE PRECISION IAL,IAK,IXIS
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/STSBSCALE/QSTSB
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/QPAR/l,k,Alcos1,Akcos2,muq,NUQ
      COMMON/SQUHIMASSM/MH02,MHC2
      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/QEXT/XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      COMMON/Z3VAUX/phIF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      COMMON/EFFPOTPAR/lu,ld,l3,l4,Rel5,Iml5,Rel6,Iml6,Rel7,Iml7,RAud,
     . RAS,K2,lPu,lPd,RlPM,IlPM
      COMMON/MH2TREE/MHuS,MHdS,MSS
      COMMON/IMALAK/IAL,IAK,IXIS

      PI=4d0*DATAN(1d0)

c      Masses and rotation matrices

      MW2=g2q/2d0*(vuq**2+vdq**2)
      MZ2=(g1q+g2q)/2d0*(vuq**2+vdq**2)

      MHTC(1)=MW2
      MHTC(2)=MHC

      DO I=1,5
       MHT2(I)=MH0(I)
      ENDDO
      MHT2(6)=MZ2

      sinbq=vuq/dsqrt(vuq**2+vdq**2)
      cosbq=vdq/dsqrt(vuq**2+vdq**2)
      DO I=1,5
        XHT(I,1)=XH(I,1)
        XHT(I,2)=XH(I,2)
        XHT(I,3)=XH(I,3)
        XHT(I,4)=XH(I,4)*cosbq
        XHT(I,5)=XH(I,4)*sinbq
        XHT(I,6)=XH(I,5)
      ENDDO
      DO J=1,6
       XHT(6,J)=0d0
      ENDDO
      XHT(6,4)=-sinbq
      XHT(6,5)=cosbq

      DO I=1,6
      DO J=1,6
       XHG(I,J)=0d0
       IF(I.eq.J)XHG(I,J)=1d0
      ENDDO
      ENDDO

c      Couplings

      aux=0d0
      DO M=1,2
      DO N=1,2

      gRH0HpHm(1,M,N)=((g1q+g2q)/2d0*vuq*XC(M,1)*XC(N,1)
     .                 +(g2q-g1q)/2d0*vuq*XC(M,2)*XC(N,2)
     . -(l**2-g2q/2d0)*vdq*(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1)))
     .                                           /dsqrt(2d0)
      gIH0HpHm(1,M,N)=0d0

      gRH0HpHm(2,M,N)=((g2q-g1q)/2d0*vdq*XC(M,1)*XC(N,1)
     .                 +(g2q+g1q)/2d0*vdq*XC(M,2)*XC(N,2)
     . -(l**2-g2q/2d0)*vuq*(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1)))
     .                                           /dsqrt(2d0)
      gIH0HpHm(2,M,N)=0d0

      gRH0HpHm(3,M,N)=(2d0*l*muq*(XC(M,1)*XC(N,1)+XC(M,2)*XC(N,2))
     .                  +l*(Alcos1+2d0*k/l*muq*DDCOS(Phi0)
     .                              +MUPQ*DDCOS(Phi01-phiP))
     .                *(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1)))/dsqrt(2d0)
      gIH0HpHm(3,M,N)=-(k*muq*DDSIN(Phi0)
     .      +(M3HQ*DDSIN(phi3Q)+l*XIFQ*DDSIN(Phi01-phIF))*l/muq)
     .          *(XC(M,1)*XC(N,2)-XC(M,2)*XC(N,1))/dsqrt(2d0)

      gRH0HpHm(4,M,N)=0d0
      gIH0HpHm(4,M,N)=-(l**2-g2q/2d0)*(vuq*sinbq+vdq*cosbq)
     .          *(XC(M,1)*XC(N,2)-XC(M,2)*XC(N,1))/dsqrt(2d0)

      gRH0HpHm(5,M,N)=(3.d0*k*muq*DDSIN(Phi0)+2d0*l
     .      *MUPQ*DDSIN(Phi01-phiP)
     .      +(M3HQ*DDSIN(phi3Q)+l*XIFQ*DDSIN(Phi01-phIF))*l/muq)
     .                *(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1))/dsqrt(2d0)
      gIH0HpHm(5,M,N)=-l*(Alcos1-2d0*k/l*muq*DDCOS(Phi0)
     .                               +MUPQ*DDCOS(Phi01-PhiP))
     .          *(XC(M,1)*XC(N,2)-XC(M,2)*XC(N,1))/dsqrt(2d0)

      ENDDO
      ENDDO

c      PIS(A,B,C,I,J,M)=XHT(I,A)*XHT(J,B)*XHT(M,C)
c     c +XHT(I,A)*XHT(M,B)*XHT(J,C)+XHT(J,A)*XHT(I,B)*XHT(M,C)
c     c +XHT(J,A)*XHT(M,B)*XHT(I,C)+XHT(M,A)*XHT(J,B)*XHT(I,C)
c     c +XHT(M,A)*XHT(I,B)*XHT(J,C)

      DO I=1,6
      DO J=1,6
      DO M=1,6

      PIS111(I,J,M)=6.d0*XHG(I,1)*XHT(J,1)*XHT(M,1)
      PIS122(I,J,M)=2d0*(XHG(I,1)*XHT(J,2)*XHT(M,2)
     . +XHG(I,2)*XHT(J,1)*XHT(M,2)+XHG(I,2)*XHT(J,2)*XHT(M,1))
      PIS133(I,J,M)=2d0*(XHG(I,1)*XHT(J,3)*XHT(M,3)
     . +XHG(I,3)*XHT(J,1)*XHT(M,3)+XHG(I,3)*XHT(J,3)*XHT(M,1))
      PIS144(I,J,M)=2d0*(XHG(I,1)*XHT(J,4)*XHT(M,4)
     . +XHG(I,4)*XHT(J,1)*XHT(M,4)+XHG(I,4)*XHT(J,4)*XHT(M,1))
      PIS155(I,J,M)=2d0*(XHG(I,1)*XHT(J,5)*XHT(M,5)
     . +XHG(I,5)*XHT(J,1)*XHT(M,5)+XHG(I,5)*XHT(J,5)*XHT(M,1))
      PIS166(I,J,M)=2d0*(XHG(I,1)*XHT(J,6)*XHT(M,6)
     . +XHG(I,6)*XHT(J,1)*XHT(M,6)+XHG(I,6)*XHT(J,6)*XHT(M,1))
      PIS211(I,J,M)=2d0*(XHG(I,2)*XHT(J,1)*XHT(M,1)
     . +XHG(I,1)*XHT(J,2)*XHT(M,1)+XHG(I,1)*XHT(J,1)*XHT(M,2))
      PIS222(I,J,M)=6.d0*XHG(I,2)*XHT(J,2)*XHT(M,2)
      PIS233(I,J,M)=2d0*(XHG(I,2)*XHT(J,3)*XHT(M,3)
     . +XHG(I,3)*XHT(J,2)*XHT(M,3)+XHG(I,3)*XHT(J,3)*XHT(M,2))
      PIS244(I,J,M)=2d0*(XHG(I,2)*XHT(J,4)*XHT(M,4)
     . +XHG(I,4)*XHT(J,2)*XHT(M,4)+XHG(I,4)*XHT(J,4)*XHT(M,2))
      PIS255(I,J,M)=2d0*(XHG(I,2)*XHT(J,5)*XHT(M,5)
     . +XHG(I,5)*XHT(J,2)*XHT(M,5)+XHG(I,5)*XHT(J,5)*XHT(M,2))
      PIS266(I,J,M)=2d0*(XHG(I,2)*XHT(J,6)*XHT(M,6)
     . +XHG(I,6)*XHT(J,2)*XHT(M,6)+XHG(I,6)*XHT(J,6)*XHT(M,2))
      PIS311(I,J,M)=2d0*(XHG(I,3)*XHT(J,1)*XHT(M,1)
     . +XHG(I,1)*XHT(J,3)*XHT(M,1)+XHG(I,1)*XHT(J,1)*XHT(M,3))
      PIS322(I,J,M)=2d0*(XHG(I,3)*XHT(J,2)*XHT(M,2)
     . +XHG(I,2)*XHT(J,3)*XHT(M,2)+XHG(I,2)*XHT(J,2)*XHT(M,3))
      PIS333(I,J,M)=6.d0*XHG(I,3)*XHT(J,3)*XHT(M,3)
      PIS344(I,J,M)=2d0*(XHG(I,3)*XHT(J,4)*XHT(M,4)
     . +XHG(I,4)*XHT(J,3)*XHT(M,4)+XHG(I,4)*XHT(J,4)*XHT(M,3))
      PIS355(I,J,M)=2d0*(XHG(I,3)*XHT(J,5)*XHT(M,5)
     . +XHG(I,5)*XHT(J,3)*XHT(M,5)+XHG(I,5)*XHT(J,5)*XHT(M,3))
      PIS366(I,J,M)=2d0*(XHG(I,3)*XHT(J,6)*XHT(M,6)
     . +XHG(I,6)*XHT(J,3)*XHT(M,6)+XHG(I,6)*XHT(J,6)*XHT(M,3))
      PIS312(I,J,M)=XHG(I,3)*XHT(J,1)*XHT(M,2)
     . +XHG(I,3)*XHT(J,2)*XHT(M,1)+XHG(I,1)*XHT(J,3)*XHT(M,2)
     . +XHG(I,1)*XHT(J,2)*XHT(M,3)+XHG(I,2)*XHT(J,1)*XHT(M,3)
     . +XHG(I,2)*XHT(J,3)*XHT(M,1)
      PIS345(I,J,M)=XHG(I,3)*XHT(J,4)*XHT(M,5)
     . +XHG(I,3)*XHT(J,5)*XHT(M,4)+XHG(I,4)*XHT(J,3)*XHT(M,5)
     . +XHG(I,4)*XHT(J,5)*XHT(M,3)+XHG(I,5)*XHT(J,4)*XHT(M,3)
     . +XHG(I,5)*XHT(J,3)*XHT(M,4)
      PIS156(I,J,M)=XHG(I,1)*XHT(J,5)*XHT(M,6)
     . +XHG(I,1)*XHT(J,6)*XHT(M,5)+XHG(I,5)*XHT(J,1)*XHT(M,6)
     . +XHG(I,5)*XHT(J,6)*XHT(M,1)+XHG(I,6)*XHT(J,5)*XHT(M,1)
     . +XHG(I,6)*XHT(J,1)*XHT(M,5)
      PIS256(I,J,M)=XHG(I,2)*XHT(J,5)*XHT(M,6)
     . +XHG(I,2)*XHT(J,6)*XHT(M,5)+XHG(I,5)*XHT(J,2)*XHT(M,6)
     . +XHG(I,5)*XHT(J,6)*XHT(M,2)+XHG(I,6)*XHT(J,5)*XHT(M,2)
     . +XHG(I,6)*XHT(J,2)*XHT(M,5)
      PIS346(I,J,M)=XHG(I,3)*XHT(J,4)*XHT(M,6)
     . +XHG(I,3)*XHT(J,6)*XHT(M,4)+XHG(I,4)*XHT(J,3)*XHT(M,6)
     . +XHG(I,4)*XHT(J,6)*XHT(M,3)+XHG(I,6)*XHT(J,4)*XHT(M,3)
     . +XHG(I,6)*XHT(J,3)*XHT(M,4)
      PIS356(I,J,M)=XHG(I,3)*XHT(J,5)*XHT(M,6)
     . +XHG(I,3)*XHT(J,6)*XHT(M,5)+XHG(I,5)*XHT(J,3)*XHT(M,6)
     . +XHG(I,5)*XHT(J,6)*XHT(M,3)+XHG(I,6)*XHT(J,5)*XHT(M,3)
     . +XHG(I,6)*XHT(J,3)*XHT(M,5)
      PIS513(I,J,M)=XHG(I,5)*XHT(J,1)*XHT(M,3)
     . +XHG(I,5)*XHT(J,3)*XHT(M,1)+XHG(I,1)*XHT(J,5)*XHT(M,3)
     . +XHG(I,1)*XHT(J,3)*XHT(M,5)+XHG(I,3)*XHT(J,1)*XHT(M,5)
     . +XHG(I,3)*XHT(J,5)*XHT(M,1)
      PIS423(I,J,M)=XHG(I,4)*XHT(J,2)*XHT(M,3)
     . +XHG(I,4)*XHT(J,3)*XHT(M,2)+XHG(I,2)*XHT(J,4)*XHT(M,3)
     . +XHG(I,2)*XHT(J,3)*XHT(M,4)+XHG(I,3)*XHT(J,2)*XHT(M,4)
     . +XHG(I,3)*XHT(J,4)*XHT(M,2)
      PIS612(I,J,M)=XHG(I,6)*XHT(J,1)*XHT(M,2)
     . +XHG(I,6)*XHT(J,2)*XHT(M,1)+XHG(I,1)*XHT(J,6)*XHT(M,2)
     . +XHG(I,1)*XHT(J,2)*XHT(M,6)+XHG(I,6)*XHT(J,1)*XHT(M,2)
     . +XHG(I,6)*XHT(J,2)*XHT(M,1)
      PIS456(I,J,M)=XHG(I,4)*XHT(J,5)*XHT(M,6)
     . +XHG(I,4)*XHT(J,6)*XHT(M,5)+XHG(I,5)*XHT(J,4)*XHT(M,6)
     . +XHG(I,5)*XHT(J,6)*XHT(M,4)+XHG(I,6)*XHT(J,5)*XHT(M,4)
     . +XHG(I,6)*XHT(J,4)*XHT(M,5)
      PIS433(I,J,M)=2d0*(XHG(I,4)*XHT(J,3)*XHT(M,3)
     . +XHG(I,3)*XHT(J,4)*XHT(M,3)+XHG(I,3)*XHT(J,3)*XHT(M,4))
      PIS613(I,J,M)=XHG(I,6)*XHT(J,1)*XHT(M,3)
     . +XHG(I,6)*XHT(J,3)*XHT(M,1)+XHG(I,1)*XHT(J,6)*XHT(M,3)
     . +XHG(I,1)*XHT(J,3)*XHT(M,6)+XHG(I,3)*XHT(J,1)*XHT(M,6)
     . +XHG(I,3)*XHT(J,6)*XHT(M,1)
      PIS466(I,J,M)=2d0*(XHG(I,4)*XHT(J,6)*XHT(M,6)
     . +XHG(I,6)*XHT(J,4)*XHT(M,6)+XHG(I,6)*XHT(J,6)*XHT(M,4))
      PIS533(I,J,M)=2d0*(XHG(I,5)*XHT(J,3)*XHT(M,3)
     . +XHG(I,3)*XHT(J,5)*XHT(M,3)+XHG(I,3)*XHT(J,3)*XHT(M,5))
      PIS623(I,J,M)=XHG(I,6)*XHT(J,2)*XHT(M,3)
     . +XHG(I,6)*XHT(J,3)*XHT(M,2)+XHG(I,2)*XHT(J,6)*XHT(M,3)
     . +XHG(I,2)*XHT(J,3)*XHT(M,6)+XHG(I,3)*XHT(J,2)*XHT(M,6)
     . +XHG(I,3)*XHT(J,6)*XHT(M,2)
      PIS566(I,J,M)=2d0*(XHG(I,5)*XHT(J,6)*XHT(M,6)
     . +XHG(I,6)*XHT(J,5)*XHT(M,6)+XHG(I,6)*XHT(J,6)*XHT(M,5))
      PIS666(I,J,M)=6.d0*XHG(I,6)*XHT(J,6)*XHT(M,6)
      PIS633(I,J,M)=2d0*(XHG(I,6)*XHT(J,3)*XHT(M,3)
     . +XHG(I,3)*XHT(J,6)*XHT(M,3)+XHG(I,3)*XHT(J,3)*XHT(M,6))
      PIS246(I,J,M)=XHG(I,2)*XHT(J,4)*XHT(M,6)
     . +XHG(I,2)*XHT(J,6)*XHT(M,4)+XHG(I,4)*XHT(J,2)*XHT(M,6)
     . +XHG(I,4)*XHT(J,6)*XHT(M,2)+XHG(I,6)*XHT(J,4)*XHT(M,2)
     . +XHG(I,6)*XHT(J,2)*XHT(M,4)

       IF(k.ne.0d0)then
      gHHH(I,J,M)=((g1q+g2q)/4.d0*(vuq*(PIS111(I,J,M)
     . +PIS144(I,J,M)-PIS122(I,J,M)-PIS155(I,J,M))+vdq*
     . (PIS222(I,J,M)+PIS255(I,J,M)-PIS211(I,J,M)-PIS244(I,J,M)))
     . +l*muq*(PIS311(I,J,M)+PIS344(I,J,M)+PIS322(I,J,M)
     . +PIS355(I,J,M))
     . +l**2*vuq*(PIS122(I,J,M)+PIS155(I,J,M)
     . +PIS133(I,J,M)+PIS166(I,J,M))+l**2*vdq*(PIS211(I,J,M)
     . +PIS244(I,J,M)+PIS233(I,J,M)+PIS266(I,J,M))
     . -l*Alcos1*(PIS312(I,J,M)-PIS345(I,J,M)-PIS156(I,J,M)
     . -PIS246(I,J,M))
     . +k/3.d0*Akcos2*(PIS333(I,J,M)-3.d0*PIS366(I,J,M))
     . +2d0*k**2*muq/l*(PIS333(I,J,M)+PIS366(I,J,M))
     . -k*DDCOS(Phi0)*(2d0*muq*(PIS312(I,J,M)-PIS345(I,J,M)
     . +PIS156(I,J,M)+PIS246(I,J,M))+l*vdq*(PIS133(I,J,M)
     . -PIS166(I,J,M)+2d0*PIS346(I,J,M))+l*vuq*(PIS233(I,J,M)
     . -PIS266(I,J,M)+2d0*PIS356(I,J,M)))
     . +k*DDSIN(Phi0)*(muq*(PIS513(I,J,M)+PIS423(I,J,M)
     . -3.d0*PIS612(I,J,M)+3.d0*PIS456(I,J,M))+l*vdq*
     . (PIS433(I,J,M)-2d0*PIS613(I,J,M)-PIS466(I,J,M))
     . +l*vuq*(PIS533(I,J,M)-2d0*PIS623(I,J,M)
     . -PIS566(I,J,M))
     . +l*vuq*vdq/muq*(3.d0*PIS633(I,J,M)-PIS666(I,J,M)))
     . -l*MUPQ*DDCOS(Phi01-phiP)*(PIS312(I,J,M)-PIS345(I,J,M)
     .                    +PIS156(I,J,M)+PIS246(I,J,M))
     . +k*MUPQ*DDCOS(Phi02-phiP)*(PIS333(I,J,M)+PIS366(I,J,M))
     . -l/muq*(M3HQ*DDSIN(phi3Q)+l*XIFQ*DDSIN(Phi01-phIF))
     .  *(PIS612(I,J,M)+PIS423(I,J,M)+PIS513(I,J,M)-PIS456(I,J,M)
     .    +vuq*vdq*l**2/muq**2*(PIS633(I,J,M)-PIS666(I,J,M)/3d0))
     . -2d0*l*MUPQ*DDSIN(Phi01-phiP)*(PIS612(I,J,M)
     .    +vuq*vdq*l**2/muq**2*(PIS666(I,J,M)/3d0-PIS633(I,J,M)))
     . -4d0/3d0*k*MUPQ*DDSIN(Phi02-phiP)*PIS666(I,J,M)
     . -l/muq*(l/muq*(XISQ*DDSIN(phiSQ)+XIFQ*MUPQ*DDSIN(PhiP-phIF))
     .      +MSPQ*DDSIN(PhiSPQ))*(PIS666(I,J,M)/3d0-PIS633(I,J,M))
     .                                         )/dsqrt(2d0)
        else
      gHHH(I,J,M)=((g1q+g2q)/4.d0*(vuq*(PIS111(I,J,M)
     . +PIS144(I,J,M)-PIS122(I,J,M)-PIS155(I,J,M))+vdq*
     . (PIS222(I,J,M)+PIS255(I,J,M)-PIS211(I,J,M)-PIS244(I,J,M)))
     . +l*muq*(PIS311(I,J,M)+PIS344(I,J,M)+PIS322(I,J,M)
     . +PIS355(I,J,M))
     . +l**2*vuq*(PIS122(I,J,M)+PIS155(I,J,M)
     . +PIS133(I,J,M)+PIS166(I,J,M))+l**2*vdq*(PIS211(I,J,M)
     . +PIS244(I,J,M)+PIS233(I,J,M)+PIS266(I,J,M))
     . -l*Alcos1*(PIS312(I,J,M)-PIS345(I,J,M)-PIS156(I,J,M)
     . -PIS246(I,J,M))
     . -l*MUPQ*DDCOS(Phi01-phiP)*(PIS312(I,J,M)-PIS345(I,J,M)
     .                    +PIS156(I,J,M)+PIS246(I,J,M))
     . -l/muq*(M3HQ*DDSIN(phi3Q)+l*XIFQ*DDSIN(Phi01-phIF))
     .  *(PIS612(I,J,M)+PIS423(I,J,M)+PIS513(I,J,M)-PIS456(I,J,M))
     . -2d0*l*MUPQ*DDSIN(Phi01-phiP)*PIS612(I,J,M)
     .                                         )/dsqrt(2d0)
       endif

c      PISS(A,B,C,D,I,J,M,N)=XHT(I,A)*XHT(J,B)*XHT(N,C)*XHT(M,D)
c     c +XHT(I,A)*XHT(J,B)*XHT(M,C)*XHT(N,D)
c     c +XHT(I,A)*XHT(N,B)*XHT(J,C)*XHT(M,D)
c     c +XHT(I,A)*XHT(N,B)*XHT(M,C)*XHT(J,D)
c     c +XHT(I,A)*XHT(M,B)*XHT(N,C)*XHT(J,D)
c     c +XHT(I,A)*XHT(M,B)*XHT(J,C)*XHT(N,D)
c     c +XHT(J,A)*XHT(I,B)*XHT(N,C)*XHT(M,D)
c     c +XHT(J,A)*XHT(I,B)*XHT(M,C)*XHT(N,D)
c     c +XHT(J,A)*XHT(N,B)*XHT(I,C)*XHT(M,D)
c     c +XHT(J,A)*XHT(N,B)*XHT(M,C)*XHT(I,D)
c     c +XHT(J,A)*XHT(M,B)*XHT(I,C)*XHT(N,D)
c     c +XHT(J,A)*XHT(M,B)*XHT(N,C)*XHT(I,D)
c     c +XHT(N,A)*XHT(I,B)*XHT(J,C)*XHT(M,D)
c     c +XHT(N,A)*XHT(I,B)*XHT(M,C)*XHT(J,D)
c     c +XHT(N,A)*XHT(J,B)*XHT(I,C)*XHT(M,D)
c     c +XHT(N,A)*XHT(J,B)*XHT(M,C)*XHT(I,D)
c     c +XHT(N,A)*XHT(M,B)*XHT(I,C)*XHT(J,D)
c     c +XHT(N,A)*XHT(M,B)*XHT(J,C)*XHT(I,D)
c     c +XHT(M,A)*XHT(I,B)*XHT(J,C)*XHT(N,D)
c     c +XHT(M,A)*XHT(I,B)*XHT(N,C)*XHT(J,D)
c     c +XHT(M,A)*XHT(J,B)*XHT(I,C)*XHT(N,D)
c     c +XHT(M,A)*XHT(J,B)*XHT(N,C)*XHT(I,D)
c     c +XHT(M,A)*XHT(N,B)*XHT(I,C)*XHT(J,D)
c     c +XHT(M,A)*XHT(N,B)*XHT(J,C)*XHT(I,D)

      PISS1111(I,J,M)=24d0*XHG(I,1)*XHG(J,1)*XHT(M,1)*XHT(M,1)
      PISS1122(I,J,M)=4.d0*(XHG(I,1)*XHG(J,1)*XHT(M,2)*XHT(M,2)
     . +XHG(I,2)*XHG(J,2)*XHT(M,1)*XHT(M,1)
     . +2d0*(XHG(I,1)*XHG(J,2)+XHG(I,2)*XHG(J,1))
     .                                 *XHT(M,1)*XHT(M,2))
      PISS1133(I,J,M)=4.d0*(XHG(I,1)*XHG(J,1)*XHT(M,3)*XHT(M,3)
     . +XHG(I,3)*XHG(J,3)*XHT(M,1)*XHT(M,1)
     . +2d0*(XHG(I,1)*XHG(J,3)+XHG(I,3)*XHG(J,1))
     .                                 *XHT(M,1)*XHT(M,3))
      PISS1144(I,J,M)=4.d0*(XHG(I,1)*XHG(J,1)*XHT(M,4)*XHT(M,4)
     . +XHG(I,4)*XHG(J,4)*XHT(M,1)*XHT(M,1)
     . +2d0*(XHG(I,1)*XHG(J,4)+XHG(I,4)*XHG(J,1))
     .                                 *XHT(M,1)*XHT(M,4))
      PISS1155(I,J,M)=4.d0*(XHG(I,1)*XHG(J,1)*XHT(M,5)*XHT(M,5)
     . +XHG(I,5)*XHG(J,5)*XHT(M,1)*XHT(M,1)
     . +2d0*(XHG(I,1)*XHG(J,5)+XHG(I,5)*XHG(J,1))
     .                                 *XHT(M,1)*XHT(M,5))
      PISS1166(I,J,M)=4.d0*(XHG(I,1)*XHG(J,1)*XHT(M,6)*XHT(M,6)
     . +XHG(I,6)*XHG(J,6)*XHT(M,1)*XHT(M,1)
     . +2d0*(XHG(I,1)*XHG(J,6)+XHG(I,6)*XHG(J,1))
     .                                 *XHT(M,1)*XHT(M,6))
      PISS2222(I,J,M)=24d0*XHG(I,2)*XHG(J,2)*XHT(M,2)*XHT(M,2)
      PISS2233(I,J,M)=4.d0*(XHG(I,2)*XHG(J,2)*XHT(M,3)*XHT(M,3)
     . +XHG(I,3)*XHG(J,3)*XHT(M,2)*XHT(M,2)
     . +2d0*(XHG(I,2)*XHG(J,3)+XHG(I,3)*XHG(J,2))
     .                                 *XHT(M,2)*XHT(M,3))
      PISS2244(I,J,M)=4.d0*(XHG(I,2)*XHG(J,2)*XHT(M,4)*XHT(M,4)
     . +XHG(I,4)*XHG(J,4)*XHT(M,2)*XHT(M,2)
     . +2d0*(XHG(I,2)*XHG(J,4)+XHG(I,4)*XHG(J,2))
     .                                 *XHT(M,2)*XHT(M,4))
      PISS2255(I,J,M)=4.d0*(XHG(I,2)*XHG(J,2)*XHT(M,5)*XHT(M,5)
     . +XHG(I,5)*XHG(J,5)*XHT(M,2)*XHT(M,2)
     . +2d0*(XHG(I,2)*XHG(J,5)+XHG(I,5)*XHG(J,2))
     .                                 *XHT(M,2)*XHT(M,5))
      PISS2266(I,J,M)=4.d0*(XHG(I,2)*XHG(J,2)*XHT(M,6)*XHT(M,6)
     . +XHG(I,6)*XHG(J,6)*XHT(M,2)*XHT(M,2)
     . +2d0*(XHG(I,2)*XHG(J,6)+XHG(I,6)*XHG(J,2))
     .                                 *XHT(M,2)*XHT(M,6))
      PISS3333(I,J,M)=24d0*XHG(I,3)*XHG(J,3)*XHT(M,3)*XHT(M,3)
      PISS3344(I,J,M)=4.d0*(XHG(I,3)*XHG(J,3)*XHT(M,4)*XHT(M,4)
     . +XHG(I,4)*XHG(J,4)*XHT(M,3)*XHT(M,3)
     . +2d0*(XHG(I,3)*XHG(J,4)+XHG(I,4)*XHG(J,3))
     .                                 *XHT(M,3)*XHT(M,4))
      PISS3355(I,J,M)=4.d0*(XHG(I,3)*XHG(J,3)*XHT(M,5)*XHT(M,5)
     . +XHG(I,5)*XHG(J,5)*XHT(M,3)*XHT(M,3)
     . +2d0*(XHG(I,3)*XHG(J,5)+XHG(I,5)*XHG(J,3))
     .                                 *XHT(M,3)*XHT(M,5))
      PISS3366(I,J,M)=4.d0*(XHG(I,3)*XHG(J,3)*XHT(M,6)*XHT(M,6)
     . +XHG(I,6)*XHG(J,6)*XHT(M,3)*XHT(M,3)
     . +2d0*(XHG(I,3)*XHG(J,6)+XHG(I,6)*XHG(J,3))
     .                                 *XHT(M,3)*XHT(M,6))
      PISS4444(I,J,M)=24d0*XHG(I,4)*XHG(J,4)*XHT(M,4)*XHT(M,4)
      PISS4455(I,J,M)=4.d0*(XHG(I,4)*XHG(J,4)*XHT(M,5)*XHT(M,5)
     . +XHG(I,5)*XHG(J,5)*XHT(M,4)*XHT(M,4)
     . +2d0*(XHG(I,4)*XHG(J,5)+XHG(I,5)*XHG(J,4))
     .                                 *XHT(M,4)*XHT(M,5))
      PISS4466(I,J,M)=4.d0*(XHG(I,4)*XHG(J,4)*XHT(M,6)*XHT(M,6)
     . +XHG(I,6)*XHG(J,6)*XHT(M,4)*XHT(M,4)
     . +2d0*(XHG(I,4)*XHG(J,6)+XHG(I,6)*XHG(J,4))
     .                                 *XHT(M,4)*XHT(M,6))
      PISS5555(I,J,M)=24d0*XHG(I,5)*XHG(J,5)*XHT(M,5)*XHT(M,5)
      PISS5566(I,J,M)=4.d0*(XHG(I,5)*XHG(J,5)*XHT(M,6)*XHT(M,6)
     . +XHG(I,6)*XHG(J,6)*XHT(M,5)*XHT(M,5)
     . +2d0*(XHG(I,5)*XHG(J,6)+XHG(I,6)*XHG(J,5))
     .                                 *XHT(M,5)*XHT(M,6))
      PISS6666(I,J,M)=24d0*XHG(I,6)*XHG(J,6)*XHT(M,6)*XHT(M,6)
      PISS1233(I,J,M)=2d0*((XHG(I,1)*XHG(J,2)+XHG(I,2)*XHG(J,1))
     .                                       *XHT(M,3)*XHT(M,3)
     . +2d0*(XHG(I,1)*XHG(J,3)*XHT(M,2)*XHT(M,3)
     . +XHG(I,3)*XHG(J,1)*XHT(M,2)*XHT(M,3)
     . +XHG(I,2)*XHG(J,3)*XHT(M,1)*XHT(M,3)
     . +XHG(I,3)*XHG(J,2)*XHT(M,1)*XHT(M,3)
     . +XHG(I,3)*XHG(J,3)*XHT(M,1)*XHT(M,2)))
      PISS1266(I,J,M)=2d0*((XHG(I,1)*XHG(J,2)+XHG(I,2)*XHG(J,1))
     .                                       *XHT(M,6)*XHT(M,6)
     . +2d0*(XHG(I,1)*XHG(J,6)*XHT(M,2)*XHT(M,6)
     . +XHG(I,6)*XHG(J,1)*XHT(M,2)*XHT(M,6)
     . +XHG(I,2)*XHG(J,6)*XHT(M,1)*XHT(M,6)
     . +XHG(I,6)*XHG(J,2)*XHT(M,1)*XHT(M,6)
     . +XHG(I,6)*XHG(J,6)*XHT(M,1)*XHT(M,2)))
      PISS3345(I,J,M)=2d0*((XHG(I,4)*XHG(J,5)+XHG(I,5)*XHG(J,4))
     .                                       *XHT(M,3)*XHT(M,3)
     . +2d0*(XHG(I,4)*XHG(J,3)*XHT(M,5)*XHT(M,3)
     . +XHG(I,3)*XHG(J,4)*XHT(M,5)*XHT(M,3)
     . +XHG(I,5)*XHG(J,3)*XHT(M,4)*XHT(M,3)
     . +XHG(I,3)*XHG(J,5)*XHT(M,4)*XHT(M,3)
     . +XHG(I,3)*XHG(J,3)*XHT(M,4)*XHT(M,5)))
      PISS1335(I,J,M)=2d0*((XHG(I,1)*XHG(J,5)+XHG(I,5)*XHG(J,1))
     .                                       *XHT(M,3)*XHT(M,3)
     . +2d0*(XHG(I,1)*XHG(J,3)*XHT(M,5)*XHT(M,3)
     . +XHG(I,3)*XHG(J,1)*XHT(M,5)*XHT(M,3)
     . +XHG(I,5)*XHG(J,3)*XHT(M,1)*XHT(M,3)
     . +XHG(I,3)*XHG(J,5)*XHT(M,1)*XHT(M,3)
     . +XHG(I,3)*XHG(J,3)*XHT(M,1)*XHT(M,5)))
      PISS2334(I,J,M)=2d0*((XHG(I,4)*XHG(J,2)+XHG(I,2)*XHG(J,4))
     .                                       *XHT(M,3)*XHT(M,3)
     . +2d0*(XHG(I,4)*XHG(J,3)*XHT(M,2)*XHT(M,3)
     . +XHG(I,3)*XHG(J,4)*XHT(M,2)*XHT(M,3)
     . +XHG(I,2)*XHG(J,3)*XHT(M,4)*XHT(M,3)
     . +XHG(I,3)*XHG(J,2)*XHT(M,4)*XHT(M,3)
     . +XHG(I,3)*XHG(J,3)*XHT(M,4)*XHT(M,2)))
      PISS4566(I,J,M)=2d0*((XHG(I,4)*XHG(J,5)+XHG(I,5)*XHG(J,4))
     .                                       *XHT(M,6)*XHT(M,6)
     . +2d0*(XHG(I,4)*XHG(J,6)*XHT(M,5)*XHT(M,6)
     . +XHG(I,6)*XHG(J,4)*XHT(M,5)*XHT(M,6)
     . +XHG(I,5)*XHG(J,6)*XHT(M,4)*XHT(M,6)
     . +XHG(I,6)*XHG(J,5)*XHT(M,4)*XHT(M,6)
     . +XHG(I,6)*XHG(J,6)*XHT(M,4)*XHT(M,5)))
      PISS1566(I,J,M)=2d0*((XHG(I,1)*XHG(J,5)+XHG(I,5)*XHG(J,1))
     .                                       *XHT(M,6)*XHT(M,6)
     . +2d0*(XHG(I,1)*XHG(J,6)*XHT(M,5)*XHT(M,6)
     . +XHG(I,6)*XHG(J,1)*XHT(M,5)*XHT(M,6)
     . +XHG(I,5)*XHG(J,6)*XHT(M,1)*XHT(M,6)
     . +XHG(I,6)*XHG(J,5)*XHT(M,1)*XHT(M,6)
     . +XHG(I,6)*XHG(J,6)*XHT(M,1)*XHT(M,5)))
      PISS2466(I,J,M)=2d0*((XHG(I,4)*XHG(J,2)+XHG(I,2)*XHG(J,4))
     .                                       *XHT(M,6)*XHT(M,6)
     . +2d0*(XHG(I,4)*XHG(J,6)*XHT(M,2)*XHT(M,6)
     . +XHG(I,6)*XHG(J,4)*XHT(M,2)*XHT(M,6)
     . +XHG(I,2)*XHG(J,6)*XHT(M,4)*XHT(M,6)
     . +XHG(I,6)*XHG(J,2)*XHT(M,4)*XHT(M,6)
     . +XHG(I,6)*XHG(J,6)*XHT(M,4)*XHT(M,2)))
      PISS1356(I,J,M)=2d0*(XHG(I,1)*XHG(J,3)*XHT(M,5)*XHT(M,6)
     . +XHG(I,1)*XHT(M,3)*XHG(J,5)*XHT(M,6)
     . +XHG(I,1)*XHT(M,3)*XHT(M,5)*XHG(J,6)
     . +XHG(J,1)*XHG(I,3)*XHT(M,5)*XHT(M,6)
     . +XHG(J,1)*XHT(M,3)*XHG(I,5)*XHT(M,6)
     . +XHG(J,1)*XHT(M,3)*XHT(M,5)*XHG(I,6)
     . +XHT(M,1)*XHG(I,3)*XHG(J,5)*XHT(M,6)
     . +XHT(M,1)*XHG(I,3)*XHT(M,5)*XHG(J,6)
     . +XHT(M,1)*XHG(J,3)*XHG(I,5)*XHT(M,6)
     . +XHT(M,1)*XHG(J,3)*XHT(M,5)*XHG(I,6)
     . +XHT(M,1)*XHT(M,3)*XHG(I,5)*XHG(J,6)
     . +XHT(M,1)*XHT(M,3)*XHG(J,5)*XHG(I,6))
      PISS1236(I,J,M)=2d0*(XHG(I,1)*XHG(J,2)*XHT(M,3)*XHT(M,6)
     . +XHG(I,1)*XHT(M,2)*XHG(J,3)*XHT(M,6)
     . +XHG(I,1)*XHT(M,2)*XHT(M,3)*XHG(J,6)
     . +XHG(J,1)*XHG(I,2)*XHT(M,3)*XHT(M,6)
     . +XHG(J,1)*XHT(M,2)*XHG(I,3)*XHT(M,6)
     . +XHG(J,1)*XHT(M,2)*XHT(M,3)*XHG(I,6)
     . +XHT(M,1)*XHG(I,2)*XHG(J,3)*XHT(M,6)
     . +XHT(M,1)*XHG(I,2)*XHT(M,3)*XHG(J,6)
     . +XHT(M,1)*XHG(J,2)*XHG(I,3)*XHT(M,6)
     . +XHT(M,1)*XHG(J,2)*XHT(M,3)*XHG(I,6)
     . +XHT(M,1)*XHT(M,2)*XHG(I,3)*XHG(J,6)
     . +XHT(M,1)*XHT(M,2)*XHG(J,3)*XHG(I,6))
      PISS2346(I,J,M)=2d0*(XHG(I,2)*XHG(J,3)*XHT(M,4)*XHT(M,6)
     . +XHG(I,2)*XHT(M,3)*XHG(J,4)*XHT(M,6)
     . +XHG(I,2)*XHT(M,3)*XHT(M,4)*XHG(J,6)
     . +XHG(J,2)*XHG(I,3)*XHT(M,4)*XHT(M,6)
     . +XHG(J,2)*XHT(M,3)*XHG(I,4)*XHT(M,6)
     . +XHG(J,2)*XHT(M,3)*XHT(M,4)*XHG(I,6)
     . +XHT(M,2)*XHG(I,3)*XHG(J,4)*XHT(M,6)
     . +XHT(M,2)*XHG(I,3)*XHT(M,4)*XHG(J,6)
     . +XHT(M,2)*XHG(J,3)*XHG(I,4)*XHT(M,6)
     . +XHT(M,2)*XHG(J,3)*XHT(M,4)*XHG(I,6)
     . +XHT(M,2)*XHT(M,3)*XHG(I,4)*XHG(J,6)
     . +XHT(M,2)*XHT(M,3)*XHG(J,4)*XHG(I,6))
      PISS3456(I,J,M)=2d0*(XHG(I,3)*XHG(J,4)*XHT(M,5)*XHT(M,6)
     . +XHG(I,3)*XHT(M,4)*XHG(J,5)*XHT(M,6)
     . +XHG(I,3)*XHT(M,4)*XHT(M,5)*XHG(J,6)
     . +XHG(J,3)*XHG(I,4)*XHT(M,5)*XHT(M,6)
     . +XHG(J,3)*XHT(M,4)*XHG(I,5)*XHT(M,6)
     . +XHG(J,3)*XHT(M,4)*XHT(M,5)*XHG(I,6)
     . +XHT(M,3)*XHG(I,4)*XHG(J,5)*XHT(M,6)
     . +XHT(M,3)*XHG(I,4)*XHT(M,5)*XHG(J,6)
     . +XHT(M,3)*XHG(J,4)*XHG(I,5)*XHT(M,6)
     . +XHT(M,3)*XHG(J,4)*XHT(M,5)*XHG(I,6)
     . +XHT(M,3)*XHT(M,4)*XHG(I,5)*XHG(J,6)
     . +XHT(M,3)*XHT(M,4)*XHG(J,5)*XHG(I,6))

      gHHHH(I,J,M)=((g1q+g2q)/8.d0*(PISS1111(I,J,M)
     . +PISS2222(I,J,M)-2d0*PISS1122(I,J,M)
     . +PISS4444(I,J,M)+PISS5555(I,J,M)
     . -2d0*PISS4455(I,J,M)+2d0*PISS1144(I,J,M)
     . +2d0*PISS2255(I,J,M)-2d0*PISS1155(I,J,M)
     . -2d0*PISS2244(I,J,M))
     . +l**2*(PISS1122(I,J,M)+PISS1133(I,J,M)
     . +PISS2233(I,J,M)+PISS4455(I,J,M)
     . +PISS4466(I,J,M)+PISS5566(I,J,M)
     . +PISS1155(I,J,M)+PISS2244(I,J,M)
     . +PISS1166(I,J,M)+PISS2266(I,J,M)
     . +PISS3344(I,J,M)+PISS3355(I,J,M))
     . -2d0*k*l*DDCOS(Phi0)*(PISS1233(I,J,M)
     . +PISS4566(I,J,M)-PISS3345(I,J,M)
     . -PISS1266(I,J,M)+2d0*PISS1356(I,J,M)
     . +2d0*PISS2346(I,J,M))
     . +k**2*(PISS3333(I,J,M)+PISS6666(I,J,M)
     . +2d0*PISS3366(I,J,M))
     . +2d0*k*l*DDSIN(Phi0)*(PISS1335(I,J,M)
     . +PISS2334(I,J,M)-2d0*PISS1236(I,J,M)
     . -PISS1566(I,J,M)-PISS2466(I,J,M)
     . +2d0*PISS3456(I,J,M)))/4.d0

      ENDDO
      ENDDO
      ENDDO


************************************************************************************************

c            A: Gauge contributions in the potential approach

c      I- Neutral mass-matrix

      aux=2d0*g2q**2*dlog(MW2/QSTSB)
     .        +(g1q+g2q)**2*dlog(MZ2/QSTSB)

      MH02(1,1)=MH02(1,1)+3.d0/64.d0/Pi**2*aux*vuq**2
      MH02(2,2)=MH02(2,2)+3.d0/64.d0/Pi**2*aux*vdq**2
      MH02(1,2)=MH02(1,2)+3.d0/64.d0/Pi**2*aux*vuq*vdq
      MH02(2,1)=MH02(2,1)+3.d0/64.d0/Pi**2*aux*vuq*vdq


c      II- Charged Higgs mass

      aux=-2d0*g1q*g2q*(dlog(MZ2/QSTSB)-1d0)

      MHC2=MHC2+3.d0/64.d0/Pi**2*aux*(vu**2+vd**2)

c        minimization equations

        MHuS=MHuS
     .    -3d0*g2q*MW2*(dlog(MW2/QSTSB)-1d0)/32d0/Pi**2
     .    -3d0*(g1q+g2q)*MZ2*(dlog(MZ2/QSTSB)-1d0)/64d0/Pi**2

        MHdS=MHdS
     .    -3d0*g2q*MW2*(dlog(MW2/QSTSB)-1d0)/32d0/Pi**2
     .    -3d0*(g1q+g2q)*MZ2*(dlog(MZ2/QSTSB)-1d0)/64d0/Pi**2


c            B: Higgs contributions

      DO I=1,5
      DO J=1,5
      dMH0(I,J)=0d0
      ENDDO
      ENDDO


c      I- Higgs/Gauge loop contribution
c             + (Remainder from the gauge contribution in the Feynman gauge)

c       a) (W,charged Higgs) loops

      aux=(2d0*MHTC(2)-MW2)*fSF1(MHTC(2),MW2,QSTSB)
     .               -fSF0(MHTC(2),QSTSB)+fSF0(MW2,QSTSB)
      dMH0(1,1)=dMH0(1,1)+g2q*aux/32d0/Pi**2*cosbq**2
      dMH0(2,2)=dMH0(2,2)+g2q*aux/32d0/Pi**2*sinbq**2
      dMH0(1,2)=dMH0(1,2)-g2q*aux/32d0/Pi**2*sinbq*cosbq
      dMH0(2,1)=dMH0(2,1)-g2q*aux/32d0/Pi**2*sinbq*cosbq
      dMH0(4,4)=dMH0(4,4)+g2q*aux/32d0/Pi**2

c       b) (Z,neutral Higgs) loops

      aux=0d0
      DO M=1,5
       aux=aux+((2d0*MHT2(M)-MZ2)*fSF1(MHT2(M),MZ2,QSTSB)
     .        -fSF0(MHT2(M),QSTSB)+fSF0(MZ2,QSTSB))*XHT(M,4)**2
      ENDDO
      dMH0(1,1)=dMH0(1,1)+(g1q+g2q)*aux/64.d0/Pi**2

      aux=0d0
      DO M=1,5
       aux=aux+((2d0*MHT2(M)-MZ2)*fSF1(MHT2(M),MZ2,QSTSB)
     .        -fSF0(MHT2(M),QSTSB)+fSF0(MZ2,QSTSB))*XHT(M,5)**2
      ENDDO
      dMH0(2,2)=dMH0(2,2)+(g1q+g2q)*aux/64.d0/Pi**2

      aux=0d0
      DO M=1,5
       aux=aux-((2d0*MHT2(M)-MZ2)*fSF1(MHT2(M),MZ2,QSTSB)
     .  -fSF0(MHT2(M),QSTSB)+fSF0(MZ2,QSTSB))*XHT(M,4)*XHT(M,5)
      ENDDO
      dMH0(1,2)=dMH0(1,2)+(g1q+g2q)*aux/64.d0/Pi**2
      dMH0(2,1)=dMH0(2,1)+(g1q+g2q)*aux/64.d0/Pi**2

      aux=0d0
      DO M=1,5
       aux=aux+((2d0*MHT2(M)-MZ2)*fSF1(MHT2(M),MZ2,QSTSB)
     .        -fSF0(MHT2(M),QSTSB)+fSF0(MZ2,QSTSB))
     .         *(cosbq*XHT(M,1)-sinbq*XHT(M,2))**2
      ENDDO
      dMH0(4,4)=dMH0(4,4)+(g1q+g2q)*aux/64.d0/Pi**2

      aux=0d0
      DO M=1,5
       aux=aux-((2d0*MHT2(M)-MZ2)*fSF1(MHT2(M),MZ2,QSTSB)
     .        -fSF0(MHT2(M),QSTSB)+fSF0(MZ2,QSTSB))
     .   *XHT(M,4)*(cosbq*XHT(M,1)-sinbq*XHT(M,2))
      ENDDO
      dMH0(1,4)=dMH0(1,4)+(g1q+g2q)*aux/64.d0/Pi**2
      dMH0(4,1)=dMH0(4,1)+(g1q+g2q)*aux/64.d0/Pi**2

      aux=0d0
      DO M=1,5
       aux=aux+((2d0*MHT2(M)-MZ2)*fSF1(MHT2(M),MZ2,QSTSB)
     .        -fSF0(MHT2(M),QSTSB)+fSF0(MZ2,QSTSB))
     .   *XHT(M,5)*(cosbq*XHT(M,1)-sinbq*XHT(M,2))
      ENDDO
      dMH0(2,4)=dMH0(2,4)+(g1q+g2q)*aux/64.d0/Pi**2
      dMH0(4,2)=dMH0(4,2)+(g1q+g2q)*aux/64.d0/Pi**2


c      II- Charged Higgs loop contribution

c       a) A0 loop - Tadpole

      aux=0d0
      DO M=1,2
      aux=aux-fSF0(MHTC(M),QSTSB)
     . *(l**2-g2q/2d0)*vdq/vuq*XC(M,1)*XC(M,2)
      ENDDO
      dMH0(1,1)=dMH0(1,1)+aux/16.d0/Pi**2

      aux=0d0
      DO M=1,2
      aux=aux-fSF0(MHTC(M),QSTSB)
     . *(l**2-g2q/2d0)*vuq/vdq*XC(M,1)*XC(M,2)
      ENDDO
      dMH0(2,2)=dMH0(2,2)+aux/16.d0/Pi**2

      aux=0d0
      DO M=1,2
      aux=aux+fSF0(MHTC(M),QSTSB)
     . *(l**2-g2q/2d0)*XC(M,1)*XC(M,2)
      ENDDO
      dMH0(1,2)=dMH0(1,2)+aux/16.d0/Pi**2
      dMH0(2,1)=dMH0(2,1)+aux/16.d0/Pi**2

      aux=0d0
      DO M=1,2
      aux=aux-fSF0(MHTC(M),QSTSB)
     . *(l**2-g2q/2d0)/sinbq/cosbq*XC(M,1)*XC(M,2)
      ENDDO
      dMH0(4,4)=dMH0(4,4)+aux/16.d0/Pi**2

      aux=0d0
      DO M=1,2
      aux=aux+fSF0(MHTC(M),QSTSB)
     . *l**2*(Alcos1+MUPQ*DDCOS(phi01-phiP))/muq*XC(M,1)*XC(M,2)
      ENDDO
      dMH0(3,3)=dMH0(3,3)+aux/16.d0/Pi**2

      aux=0d0
      DO M=1,2
      aux=aux+fSF0(MHTC(M),QSTSB)
     . *l**2/muq*(Alcos1+4.d0*k/l*muq*DDCOS(Phi0))*XC(M,1)*XC(M,2)
      ENDDO
      dMH0(5,5)=dMH0(5,5)+aux/16.d0/Pi**2

      aux=0d0
      DO M=1,2
      aux=aux+4.d0*fSF0(MHTC(M),QSTSB)
     . *k*l*DDSIN(Phi0)*XC(M,1)*XC(M,2)
      ENDDO
      dMH0(3,5)=dMH0(3,5)+aux/16.d0/Pi**2
      dMH0(5,3)=dMH0(5,3)+aux/16.d0/Pi**2

c        minimization conditions

      aux=0d0
      DO M=1,2
      aux=aux+fSF0(MHTC(M),QSTSB)
     . *((g1q+g2q)/4d0*XC(M,1)*XC(M,1)
     .  +(-g1q+g2q)/4d0*XC(M,2)*XC(M,2)
     .  -(l**2-g2q/2d0)*vdq/vuq*XC(M,1)*XC(M,2))
      ENDDO
      MHuS=MHuS-aux/16.d0/Pi**2

      aux=0d0
      DO M=1,2
      aux=aux+fSF0(MHTC(M),QSTSB)
     . *((g1q+g2q)/4d0*XC(M,2)*XC(M,2)
     .  +(-g1q+g2q)/4d0*XC(M,1)*XC(M,1)
     .  -(l**2-g2q/2d0)*vuq/vdq*XC(M,1)*XC(M,2))
      ENDDO
      MHdS=MHdS-aux/16.d0/Pi**2

      aux=0d0
      DO M=1,2
      aux=aux+fSF0(MHTC(M),QSTSB)
     . *(l**2*(XC(M,1)*XC(M,1)+XC(M,2)*XC(M,2))
     .  +l**2/muq*(Alcos1+MUPQ*DDCOS(phi01-phiP)
     .             +2d0*k/l*muq*DDCOS(phi0))*XC(M,1)*XC(M,2))
      ENDDO
      MSS=MSS-aux/16.d0/Pi**2

      aux=0d0
      DO M=1,2
      aux=aux+fSF0(MHTC(M),QSTSB)*XC(M,1)*XC(M,2)
     . *(l/muq*(M3HQ*DDSIN(phi3q)+l*XIFQ*DDSIN(phi01-phiF))
     .      +3d0*k*muq*DDSIN(phi0)+2d0*l*MUPQ*DDSIN(phi01-phiP))
      ENDDO

      IF(k.ne.0d0)THEN
      IAk=IAk+aux/k*(l/muq)**2/16.d0/Pi**2
      ELSE
      IXIS=IXIS+aux/16.d0/Pi**2
      ENDIF


c       b) B0 loop

      DO I=1,5
      DO J=1,5
      DO M=1,2
      DO N=1,2
       aux=-(gRH0HpHm(I,M,N)*gRH0HpHm(J,M,N)+gIH0HpHm(I,M,N)*
     .        gIH0HpHm(J,M,N))*fSF1(MHTC(M),MHTC(N),QSTSB)
       dMH0(I,J)=dMH0(I,J)+aux/16.d0/Pi**2
      ENDDO
      ENDDO
      ENDDO
      ENDDO


c      III- Neutral Higgs loop contribution

c       a) A0 loop

      DO I=1,3
      DO J=1,3
       aux=0d0
       DO M=1,6
      aux=aux-gHHHH(I,J,M)*fSF0(MHT2(M),QSTSB)
       ENDDO
       dMH0(I,J)=dMH0(I,J)+aux/32d0/Pi**2
      ENDDO
      ENDDO

      aux=0d0
      DO M=1,6
      aux=aux-(cosbq**2*gHHHH(4,4,M)+sinbq**2*gHHHH(5,5,M)
     .           +sinbq*cosbq*(gHHHH(4,5,M)+gHHHH(5,4,M)))
     .                                    *fSF0(MHT2(M),QSTSB)
      ENDDO
      dMH0(4,4)=dMH0(4,4)+aux/32d0/Pi**2

      aux=0d0
      DO M=1,6
      aux=aux-gHHHH(6,6,M)*fSF0(MHT2(M),QSTSB)
      ENDDO
      dMH0(5,5)=dMH0(5,5)+aux/32d0/Pi**2

      aux=0d0
      DO M=1,6
      aux=aux-(cosbq*gHHHH(4,6,M)+sinbq*gHHHH(5,6,M))
     .                       *fSF0(MHT2(M),QSTSB)
      ENDDO
      dMH0(4,5)=dMH0(4,5)+aux/32d0/Pi**2
      dMH0(5,4)=dMH0(5,4)+aux/32d0/Pi**2

      DO I=1,3
       aux=0d0
       DO M=1,6
      aux=aux-(cosbq*gHHHH(4,I,M)+sinbq*gHHHH(5,I,M))
     .                       *fSF0(MHT2(M),QSTSB)
       ENDDO
       dMH0(4,I)=dMH0(4,I)+aux/32d0/Pi**2
       dMH0(I,4)=dMH0(I,4)+aux/32d0/Pi**2
      ENDDO

      DO I=1,3
       aux=0d0
       DO M=1,6
      aux=aux-gHHHH(6,I,M)*fSF0(MHT2(M),QSTSB)
       ENDDO
       dMH0(5,I)=dMH0(5,I)+aux/32d0/Pi**2
       dMH0(I,5)=dMH0(I,5)+aux/32d0/Pi**2
      ENDDO


c       b) Tadpoles

      aux=0d0
      DO M=1,6
       aux=aux-gHHH(1,M,M)*fSF0(MHT2(M),QSTSB)/dsqrt(2d0)/vuq
      ENDDO
      dMH0(1,1)=dMH0(1,1)-aux/32d0/Pi**2
      dMH0(4,4)=dMH0(4,4)-aux/32d0/Pi**2*cosbq**2
      MHuS=MHuS+aux/32d0/Pi**2

      aux=0d0
      DO M=1,6
       aux=aux-gHHH(2,M,M)*fSF0(MHT2(M),QSTSB)/dsqrt(2d0)/vdq
      ENDDO
      dMH0(2,2)=dMH0(2,2)-aux/32d0/Pi**2
      dMH0(4,4)=dMH0(4,4)-aux/32d0/Pi**2*sinbq**2
      MHdS=MHdS+aux/32d0/Pi**2

      aux=0d0
      DO M=1,6
       aux=aux-gHHH(3,M,M)*fSF0(MHT2(M),QSTSB)/dsqrt(2d0)*l/muq
      ENDDO
      dMH0(3,3)=dMH0(3,3)-aux/32d0/Pi**2
      dMH0(5,5)=dMH0(5,5)-aux/32d0/Pi**2
      MSS=MSS+aux/32d0/Pi**2

      aux=0d0
      DO M=1,6
       aux=aux-gHHH(4,M,M)*fSF0(MHT2(M),QSTSB)/dsqrt(2d0)/vdq
      ENDDO
      dMH0(1,4)=dMH0(1,4)-aux/32d0/Pi**2*sinbq
      dMH0(4,1)=dMH0(4,1)-aux/32d0/Pi**2*sinbq
      dMH0(2,4)=dMH0(2,4)-aux/32d0/Pi**2*cosbq
      dMH0(4,2)=dMH0(4,2)-aux/32d0/Pi**2*cosbq
      dMH0(1,5)=dMH0(1,5)-aux/32d0/Pi**2*vdq*l/muq
      dMH0(5,1)=dMH0(5,1)-aux/32d0/Pi**2*vdq*l/muq
      dMH0(2,5)=dMH0(2,5)-aux/32d0/Pi**2*vuq*l/muq
      dMH0(5,2)=dMH0(5,2)-aux/32d0/Pi**2*vuq*l/muq
      dMH0(3,4)=dMH0(3,4)-aux/32d0/Pi**2*l/muq
     .                             *dsqrt(vuq**2+vdq**2)
      dMH0(4,3)=dMH0(4,3)-aux/32d0/Pi**2*l/muq
     .                             *dsqrt(vuq**2+vdq**2)
      dMH0(3,5)=dMH0(3,5)+2d0*aux/32d0/Pi**2*vuq*vdq*(l/muq)**2
      dMH0(5,3)=dMH0(5,3)+2d0*aux/32d0/Pi**2*vuq*vdq*(l/muq)**2
      IAl=IAl+aux/32d0/Pi**2/muq

      aux=0d0
      DO M=1,6
       aux=aux-gHHH(6,M,M)*fSF0(MHT2(M),QSTSB)/dsqrt(2d0)*l/muq
      ENDDO
      dMH0(3,5)=dMH0(3,5)-2d0*aux/32d0/Pi**2
      dMH0(5,3)=dMH0(5,3)-2d0*aux/32d0/Pi**2
      IF(k.ne.0d0)THEN
       IAk=IAk-aux/32d0/Pi**2/k*(l/muq)**2
      ELSE
       IXIS=IXIS-aux/32d0/Pi**2
      ENDIF


c       c) B0 loop

      DO I=1,3
      DO J=1,3
       aux=0d0
       DO M=1,6
       DO N=1,6
        aux=aux-gHHH(I,M,N)*gHHH(J,N,M)
     .                     *fSF1(MHT2(M),MHT2(N),QSTSB)
       ENDDO
       ENDDO
       dMH0(I,J)=dMH0(I,J)+aux/32d0/Pi**2
      ENDDO
      ENDDO

       aux=0d0
       DO M=1,6
       DO N=1,6
        aux=aux-(gHHH(4,M,N)*gHHH(4,N,M)*cosbq**2
     . +gHHH(5,M,N)*gHHH(5,N,M)*sinbq**2+sinbq*cosbq*
     . (gHHH(4,M,N)*gHHH(5,N,M)+gHHH(5,M,N)*gHHH(4,N,M)))
     .                     *fSF1(MHT2(M),MHT2(N),QSTSB)
       ENDDO
       ENDDO
       dMH0(4,4)=dMH0(4,4)+aux/32d0/Pi**2

       aux=0d0
       DO M=1,6
       DO N=1,6
        aux=aux-gHHH(6,M,N)*gHHH(6,N,M)
     .                     *fSF1(MHT2(M),MHT2(N),QSTSB)
       ENDDO
       ENDDO
       dMH0(5,5)=dMH0(5,5)+aux/32d0/Pi**2

       aux=0d0
       DO M=1,6
       DO N=1,6
        aux=aux-(cosbq*gHHH(4,M,N)+sinbq*gHHH(5,M,N))
     .                 *gHHH(6,N,M)*fSF1(MHT2(M),MHT2(N),QSTSB)
       ENDDO
       ENDDO
       dMH0(4,5)=dMH0(4,5)+aux/32d0/Pi**2
       dMH0(5,4)=dMH0(5,4)+aux/32d0/Pi**2

      DO I=1,3
       aux=0d0
       DO M=1,6
       DO N=1,6
        aux=aux-(cosbq*gHHH(4,M,N)+sinbq*gHHH(5,M,N))
     .                *gHHH(I,N,M)*fSF1(MHT2(M),MHT2(N),QSTSB)
       ENDDO
       ENDDO
       dMH0(4,I)=dMH0(4,I)+aux/32d0/Pi**2
       dMH0(I,4)=dMH0(I,4)+aux/32d0/Pi**2

       aux=0d0
       DO M=1,6
       DO N=1,6
        aux=aux-gHHH(6,M,N)*gHHH(I,N,M)
     .                          *fSF1(MHT2(M),MHT2(N),QSTSB)
       ENDDO
       ENDDO
       dMH0(5,I)=dMH0(5,I)+aux/32d0/Pi**2
       dMH0(I,5)=dMH0(I,5)+aux/32d0/Pi**2
      ENDDO


c      IV- Corrections to the charged-Higgs mass

c       a) Higgs/gauge loop
c            (+ remainder from gauge contr. in the Feynman gauge)

      aux=-g2q*fSF0(MW2,QSTSB)
     .    -(g1q-g2q)**2/2d0/(g1q+g2q)*fSF0(MZ2,QSTSB)
      dMHC=aux/32d0/Pi**2

      aux=0d0
      DO M=1,5
       aux=aux+g2q*(cosbq**2*(XHT(M,1)**2+XHT(M,4)**2)
     .          +sinbq**2*(XHT(M,2)**2+XHT(M,5)**2)
     . -2d0*sinbq*cosbq*(XHT(M,1)*XHT(M,2)-XHT(M,4)*XHT(M,5)))
     .   *((2d0*MHT2(M)-MW2)*fSF1(MW2,MHT2(M),QSTSB)
     .      -fSF0(MHT2(M),QSTSB)+2d0*fSF0(MW2,QSTSB))
      ENDDO
      dMHC=dMHC+aux/64.d0/Pi**2

      aux=g1q*g2q/(g1q+g2q)*fSF0(MHTC(2),QSTSB)
      dMHC=dMHC+aux/16.d0/Pi**2

      aux=(g1q-g2q)**2/(g1q+g2q)/4.d0
     .   *((2d0*MHTC(2)-MZ2)*fSF1(MZ2,MHTC(2),QSTSB)
     .      -fSF0(MHTC(2),QSTSB)+2d0*fSF0(MZ2,QSTSB))
      dMHC=dMHC+aux/16.d0/Pi**2


c       b) Charged Higgs A0 loop - Tadpole

      aux=0d0
      DO M=1,2
      aux=aux-(l**2*(cosbq*XC(M,2)+sinbq*XC(M,1))**2
     . +(g1q+g2q)/4.d0*((2d0*cosbq**2-sinbq**2)*XC(M,1)**2
     . +(2d0*sinbq**2-cosbq**2)*XC(M,2)**2-2d0*sinbq*cosbq
     .                              *XC(M,1)*XC(M,2))
     . +(l**2-g2q/2d0)*(cosbq**4+sinbq**4)/sinbq/cosbq
     .                              *XC(M,1)*XC(M,2)
     . -(g2q+g1q*(cosbq**2-sinbq**2)*(XC(M,1)**2-XC(M,2)**2))/4.d0)
     .                           *fSF0(MHTC(M),QSTSB)
      ENDDO
      dMHC=dMHC+aux/16.d0/Pi**2


c       c) Neutral Higgs A0 loop - Tadpole

      aux=0d0
      DO M=1,6
       aux=aux-((g2q+g1q*(cosbq**2-sinbq**2))/4.d0
     .                                  *(XHT(M,1)**2+XHT(M,4)**2)
     . +(g2q-g1q*(cosbq**2-sinbq**2))/4.d0
     .                                  *(XHT(M,2)**2+XHT(M,5)**2)
     . -(2d0*l**2-g2q)*(XHT(M,1)*XHT(M,2)-XHT(M,4)*XHT(M,5))
     .                                                 *sinbq*cosbq
     . +l*(l+2d0*k*DDCOS(Phi0)*sinbq*cosbq)*XHT(M,3)**2
     . +l*(l-2d0*k*DDCOS(Phi0)*sinbq*cosbq)*XHT(M,6)**2
     . +4.d0*k*l*DDSIN(Phi0)*sinbq*cosbq*XHT(M,3)*XHT(M,6))
     .                         *fSF0(MHT2(M),QSTSB)
      ENDDO
      dMHC=dMHC+aux/32d0/Pi**2

      aux=0d0
      DO M=1,6
       aux=aux-(gHHH(1,M,M)/vuq*cosbq**2+gHHH(2,M,M)/vdq*sinbq**2)
     .                             *fSF0(MHT2(M),QSTSB)/dsqrt(2d0)
      ENDDO
      dMHC=dMHC-aux/32d0/Pi**2


c       d) B0 loop

      aux=0d0
      DO M=1,2
      DO N=1,6
      DO I=1,5
      DO J=1,5
       aux=aux-XHT(N,I)*XHT(N,J)*(gRH0HpHm(I,M,2)*gRH0HpHm(J,2,M)
     . +gIH0HpHm(I,M,2)*gIH0HpHm(J,2,M))*fSF1(MHT2(N),MHTC(M),QSTSB)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      dMHC=dMHC+aux/16/Pi**2

c      print*,'dMH0_1*',dMH0(1,1),dMH0(1,2),dMH0(1,3),dMH0(1,4),
c     c dMH0(1,5)
c      print*,'dMH0_2*',dMH0(2,1),dMH0(2,2),dMH0(2,3),dMH0(2,4),
c     c dMH0(2,5)
c      print*,'dMH0_3*',dMH0(3,1),dMH0(3,2),dMH0(3,3),dMH0(3,4),
c     c dMH0(3,5)
c      print*,'dMH0_4*',dMH0(4,1),dMH0(4,2),dMH0(4,3),dMH0(4,4),
c     c dMH0(4,5)
c      print*,'dMH0_5*',dMH0(5,1),dMH0(5,2),dMH0(5,3),dMH0(5,4),
c     c dMH0(5,5)
c      print*,'dMHC',dMHC


      aux=Max(dabs(XIS),dabs(MSP),dabs(XIF),dabs(MUP),dabs(M3H))
      IF(aux.le.1d-4)THEN
c Z3-conserving version
c      V- Reconstruction of the Z3-conserving parameters
c                    (-sign because we computed self-tadpole and dm2=-that)

      RASH=l/muq*(dMH0(5,5)-vuq*vdq*l/muq*(2d0*l/muq*sinbq*cosbq
     . *dMH0(4,4)-dMH0(4,5)/dsqrt(vuq**2+vdq**2)))/3.d0

      K2H=-(l/muq)**2*(dMH0(3,3)+dMH0(5,5)/3.d0
     .    -4.d0/3.d0*vuq*vdq*sinbq*cosbq*(l/muq)**2*dMH0(4,4))/4.d0

      RAudH=-(2d0*l/muq*sinbq*cosbq*dMH0(4,4)
     .                         +dMH0(4,5)/dsqrt(vuq**2+vdq**2))/3.d0

      RlPMH=-l/muq*(l/muq*sinbq*cosbq*dMH0(4,4)
     .                         -dMH0(4,5)/dsqrt(vuq**2+vdq**2))/3.d0

      IlPMH=-(dMH0(3,5)/vuq/vdq-dMH0(3,4)*l/muq/dsqrt(vuq**2+vdq**2))
     .                                                         /3.d0

      lPuH=-l/muq*(dMH0(1,3)+cosbq/3.d0*(-dMH0(4,5)
     .       +4.d0*vuq*cosbq*l/muq*dMH0(4,4)))/2d0/vuq

      lPdH=-l/muq*(dMH0(2,3)+sinbq/3.d0*(-dMH0(4,5)
     .       +4.d0*vuq*cosbq*l/muq*dMH0(4,4)))/2d0/vdq

      luH=-(dMH0(1,1)-cosbq**2*dMH0(4,4))/2d0/vuq**2

      ldH=-(dMH0(2,2)-sinbq**2*dMH0(4,4))/2d0/vdq**2

      l3H=-(dMH0(1,2)+sinbq*cosbq*(2d0*dMHC-dMH0(4,4)))/2d0/vuq/vdq

      l4H=-(dMH0(4,4)-dMHC)/(vuq**2+vdq**2)

      Il5H=-(l/muq)**2/vuq/vdq*(dMH0(3,5)/vuq/vdq
     .                -4.d0*dMH0(3,4)*l/muq/dsqrt(vuq**2+vdq**2))/3.d0

      Il6H=-(dMH0(1,4)/dsqrt(vuq**2+vdq**2)+Il5H*vdq)/2d0/vuq

      Il7H=-(dMH0(2,4)/dsqrt(vuq**2+vdq**2)+Il5H*vuq)/2d0/vdq

      RAS=RAS+RASH
      K2=K2+K2H
      RAud=RAud+RAudH
      RlPM=RlPM+RlPMH
      IlPM=IlPM+IlPMH
      lPu=lPu+lPuH
      lPd=lPd+lPdH
      lu=lu+luH
      ld=ld+ldH
      l3=l3+l3H
      l4=l4+l4H
      Iml5=Iml5+Il5H
      Iml6=Iml6+Il6H
      Iml7=Iml7+Il7H


c      VI- Inclusion of the corrections to the Higgs mass

      aux=2d0*luH*vuq**2+(RAudH+RlPMH*muq/l)*muq/l*vdq/vuq

      MH02(1,1)=MH02(1,1)+aux

      aux=2d0*ldH*vdq**2+(RAudH+RlPMH*muq/l)*muq/l*vuq/vdq

      MH02(2,2)=MH02(2,2)+aux

      aux=-(RAudH+RlPMH*muq/l)*muq/l+2d0*(l3H+l4H)*vuq*vdq

      MH02(1,2)=MH02(1,2)+aux
      MH02(2,1)=MH02(2,1)+aux

      aux=-(RAudH+2d0*RlPMH*muq/l)*vdq+2d0*lPuH*vuq*muq/l

      MH02(1,3)=MH02(1,3)+aux
      MH02(3,1)=MH02(3,1)+aux

      aux=-(RAudH+2d0*RlPMH*muq/l)*vuq+2d0*lPdH*vdq*muq/l

      MH02(2,3)=MH02(2,3)+aux
      MH02(3,2)=MH02(3,2)+aux

      MH02(3,3)=MH02(3,3)+(RASH+4.d0*K2H*muq/l)*muq/l
     .                                      +l/muq*vuq*vdq*RAudH

      aux=(RAudH+RlPMH*muq/l)*muq/l

      MH02(4,4)=MH02(4,4)+aux*(vu**2+vd**2)/vuq/vdq

      aux=(RAudH-2d0*RlPMH*muq/l)

      MH02(4,5)=MH02(4,5)+aux*dsqrt(vu**2+vd**2)
      MH02(5,4)=MH02(5,4)+aux*dsqrt(vu**2+vd**2)

      aux=-3.d0*muq/l*RASH
     .            +(RAudH+4.d0*RlPMH*muq/l)*vuq*vdq*l/muq

      MH02(5,5)=MH02(5,5)+aux

      aux=(2d0*Il6H*vuq-Il5H*vdq)

      MH02(1,4)=MH02(1,4)+aux*dsqrt(vu**2+vd**2)
      MH02(4,1)=MH02(4,1)+aux*dsqrt(vu**2+vd**2)

      aux=(2d0*Il7H*vdq-Il5H*vuq)

      MH02(2,4)=MH02(2,4)+aux*dsqrt(vu**2+vd**2)
      MH02(4,2)=MH02(4,2)+aux*dsqrt(vu**2+vd**2)

      aux=-3.d0*IlPMH*muq/l
     .      -(Il5H*vuq*vdq-2d0*Il6H*vuq**2)*l/muq

      MH02(1,5)=MH02(1,5)+aux*vdq
      MH02(5,1)=MH02(5,1)+aux*vdq

      aux=-3.d0*IlPMH*muq/l
     .      -(Il5H*vuq*vdq-2d0*Il7H*vdq**2)*l/muq

      MH02(2,5)=MH02(2,5)+aux*vuq
      MH02(5,2)=MH02(5,2)+aux*vuq

      aux=IlPMH*muq/l-Il5H*vuq*vdq*l/muq

      MH02(3,4)=MH02(3,4)+aux*dsqrt(vu**2+vd**2)
      MH02(4,3)=MH02(4,3)+aux*dsqrt(vu**2+vd**2)

      aux=(4.d0*IlPMH-Il5H*vuq*vdq/(muq/l)**2)*vuq*vdq

      MH02(3,5)=MH02(3,5)+aux
      MH02(5,3)=MH02(5,3)+aux

      aux=(RAudH+RlPMH*muq/l)*muq/l-l4H*vuq*vdq

      MHC2=MHC2+aux*(vu**2+vd**2)/vuq/vdq


      ELSE
c Z3-violating version
c      Vbis- Inclusion of the corrections to the Higgs mass

      MH02(1,1)=MH02(1,1)-dMH0(1,1)

      MH02(2,2)=MH02(2,2)-dMH0(2,2)

      MH02(1,2)=MH02(1,2)-dMH0(1,2)
      MH02(2,1)=MH02(2,1)-dMH0(2,1)

      MH02(1,3)=MH02(1,3)-dMH0(1,3)
      MH02(3,1)=MH02(3,1)-dMH0(3,1)

      MH02(2,3)=MH02(2,3)-dMH0(2,3)
      MH02(3,2)=MH02(3,2)-dMH0(3,2)

      MH02(3,3)=MH02(3,3)-dMH0(3,3)

      MH02(4,4)=MH02(4,4)-dMH0(4,4)*(vu**2+vd**2)/(vuq**2+vdq**2)

      MH02(4,5)=MH02(4,5)-dMH0(4,5)
     .                       *dsqrt((vu**2+vd**2)/(vuq**2+vdq**2))
      MH02(5,4)=MH02(5,4)-dMH0(5,4)
     .                       *dsqrt((vu**2+vd**2)/(vuq**2+vdq**2))

      MH02(5,5)=MH02(5,5)-dMH0(5,5)

      MH02(1,4)=MH02(1,4)-dMH0(1,4)
     .                       *dsqrt((vu**2+vd**2)/(vuq**2+vdq**2))
      MH02(4,1)=MH02(4,1)-dMH0(4,1)
     .                       *dsqrt((vu**2+vd**2)/(vuq**2+vdq**2))

      MH02(2,4)=MH02(2,4)-dMH0(2,4)
     .                       *dsqrt((vu**2+vd**2)/(vuq**2+vdq**2))
      MH02(4,2)=MH02(4,2)-dMH0(4,2)
     .                       *dsqrt((vu**2+vd**2)/(vuq**2+vdq**2))

      MH02(1,5)=MH02(1,5)-dMH0(1,5)
      MH02(5,1)=MH02(5,1)-dMH0(5,1)

      MH02(2,5)=MH02(2,5)-dMH0(2,5)
      MH02(5,2)=MH02(5,2)-dMH0(5,2)

      MH02(3,4)=MH02(3,4)-dMH0(3,4)
     .                       *dsqrt((vu**2+vd**2)/(vuq**2+vdq**2))
      MH02(4,3)=MH02(4,3)-dMH0(4,3)
     .                       *dsqrt((vu**2+vd**2)/(vuq**2+vdq**2))

      MH02(3,5)=MH02(3,5)-dMH0(3,5)
      MH02(5,3)=MH02(5,3)-dMH0(5,3)

      MHC2=MHC2-dMHC*(vu**2+vd**2)/(vuq**2+vdq**2)


c      VIbis- Effective potential parameters (under construction)
      ENDIF


c      print*,'MH02_1*',MH02(1,1),MH02(1,2),MH02(1,3),MH02(1,4),MH02(1,5)
c      print*,'MH02_2*',MH02(2,1),MH02(2,2),MH02(2,3),MH02(2,4),MH02(2,5)
c      print*,'MH02_3*',MH02(3,1),MH02(3,2),MH02(3,3),MH02(3,4),MH02(3,5)
c      print*,'MH02_4*',MH02(4,1),MH02(4,2),MH02(4,3),MH02(4,4)
c     c /(ZHu*ZHd),MH02(4,5)/dsqrt(Zs*ZHu*ZHd)
c      print*,'MH02_5*',MH02(5,1),MH02(5,2),MH02(5,3),MH02(5,4),MH02(5,5)
c      print*,'MHC2',MHC2                                   !/(ZHu*ZHd)

      RETURN
      END

************************************************************************************************

      SUBROUTINE MHIGGSLOOP_POLE_CPV(IFAIL)

c         One-loop corrections to the Higgs potential
c                 - DRbar-masses and pole corrections 
c      - The corrected squared-mass matrices in the Higgs sector (stored in
c        the common SQUHIMASSM) are rescaled according to the wave-function
c        renormalization factors (see init_CPV.f) and diagonalized.
c      - Momentum-depENDent (pole) corrections are added to the squared-mass 
c        eigenvalues.
c      - The pole-corrected squared-masses for the charged-Higgs, MHC2, and 
c        neutral-Higgs, MH0(i), (i=1..5), replace their tree-level counterparts
c        in the common HISPEC. Same thing for the rotation matrices XC and XH.
c      This concludes the calculation of the radiative corrections to the Higgs
c      squared masses. IF a squared-mass becomes negative, IFAIL.NE.0.

c      The trilinear Higgs-Sfermion couplings have also been computed 
c      meanwhile, and they are stored in the common HISFCOUP.

      IMPLICIT NONE

      INTEGER I,J,M,N,IFAIL

      DOUBLE PRECISION aux,Pi,VALPH(5),VECPH(5,5),NMB0,fSF0,fSF1
      DOUBLE PRECISION muH2,MW2,MZ2,sinbq,cosbq,MHT2(6),XHT(6,6)
      DOUBLE PRECISION MHTC(2),Ytau,gRHSS,gIHSS,XHG(5,6),DELT(2,2)
      DOUBLE PRECISION PIS111,PIS122,PIS133,PIS144,PIS155,PIS166,
     . PIS211,PIS222,PIS233,PIS244,PIS255,PIS266,PIS311,PIS322,
     . PIS333,PIS344,PIS355,PIS366,PIS312,PIS345,PIS156,PIS256,
     . PIS346,PIS356,PIS513,PIS423,PIS612,PIS456,PIS433,PIS613,
     . PIS466,PIS533,PIS623,PIS566,PIS666,PIS633,PIS246,gHSS
      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION l,k,Alcos1,Akcos2,muq,nuq
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      DOUBLE PRECISION Ytq,Ybq,MTOPQ,MBOTQ
      DOUBLE PRECISION mur,M1r,M2r,msi
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION MST2(2),UT(2,2,2),MSB2(2),UB(2,2,2),MSL2(2),
     . UTAU(2,2,2),MSNT2
      DOUBLE PRECISION MSU2(2),MSD2(2),MSE2(2),MSNE2,MSMU2(2),
     . UMU(2,2,2)
      DOUBLE PRECISION MSQ3,MSU3,MSD3,AT,AB
      DOUBLE PRECISION MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      DOUBLE PRECISION MHC,XC(2,2),MH0T(5),XH0(5,5),MA2
      DOUBLE PRECISION MH02(5,5),MHC2
      DOUBLE PRECISION MH0(5),XH(5,5)
      DOUBLE PRECISION GRHSTST(5,2,2),GRHSBSB(5,2,2),GRHSLSL(5,2,2),
     . GRHSUSU(5,2,2),GRHSDSD(5,2,2),GRHSESE(5,2,2),GRHSNSN(5)
      DOUBLE PRECISION GIHSTST(5,2,2),GIHSBSB(5,2,2),GIHSLSL(5,2,2)
      DOUBLE PRECISION GRHCSTSB(2,2),GRHCSNSL(2),GRHCSUSD(2,2),
     . GRHCSNSE(2),GIHCSTSB(2,2),GIHCSNSL(2)
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      DOUBLE PRECISION phIF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/STSBSCALE/QSTSB
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/QPAR/l,k,Alcos1,Akcos2,muq,NUQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/SMFERM/mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      COMMON/QQUARK/Ytq,Ybq,MTOPQ,MBOTQ
      COMMON/GAUGINOPAR/mur,M1r,M2r,msi
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/HISPEC/MHC,XC,MH0T,XH0,MA2
      COMMON/SFERM3SPEC/MST2,UT,MSB2,UB,MSL2,UTAU,MSNT2
      COMMON/SFERM1SPEC/MSU2,MSD2,MSE2,MSNE2,MSMU2,UMU
      COMMON/RADCOR2/MSQ3,MSU3,MSD3,AT,AB
      COMMON/SLEPPAR/MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      COMMON/SQUHIMASSM/MH02,MHC2
      COMMON/HISFCOUP/GRHSTST,GRHSBSB,GRHSLSL,GRHSUSU,GRHSDSD,
     . GRHSESE,GRHSNSN,GIHSTST,GIHSBSB,GIHSLSL,GRHCSTSB,GRHCSNSL,
     . GRHCSUSD,GRHCSNSE,GIHCSTSB,GIHCSNSL
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/QEXT/XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      COMMON/Z3VAUX/phIF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      COMMON/HIGGSMS/muH2


      PI=4d0*DATAN(1d0)
      DELT(1,1)=1d0
      DELT(1,2)=0d0
      DELT(2,1)=0d0
      DELT(2,2)=1d0

c            A: Higgs DRbar masses

c      I- Rescaling (Wave-function renorm.)

      MH02(1,1)=MH02(1,1)/ZHU
      MH02(1,2)=MH02(1,2)/dsqrt(ZHU*ZHD)
      MH02(2,1)=MH02(1,2)
      MH02(1,3)=MH02(1,3)/dsqrt(ZHU*ZS)
      MH02(3,1)=MH02(1,3)
      MH02(1,4)=MH02(1,4)/dsqrt(ZHD)/ZHU
      MH02(4,1)=MH02(1,4)
      MH02(1,5)=MH02(1,5)/dsqrt(ZHU*ZS)
      MH02(5,1)=MH02(1,5)
      MH02(2,2)=MH02(2,2)/ZHD
      MH02(2,3)=MH02(2,3)/dsqrt(ZHD*ZS)
      MH02(3,2)=MH02(2,3)
      MH02(2,4)=MH02(2,4)/dsqrt(ZHU)/ZHD
      MH02(4,2)=MH02(2,4)
      MH02(2,5)=MH02(2,5)/dsqrt(ZHD*ZS)
      MH02(5,2)=MH02(2,5)
      MH02(3,3)=MH02(3,3)/ZS
      MH02(3,4)=MH02(3,4)/dsqrt(ZHU*ZHD*ZS)
      MH02(4,3)=MH02(3,4)
      MH02(3,5)=MH02(3,5)/ZS
      MH02(5,3)=MH02(3,5)
      MH02(4,4)=MH02(4,4)/(ZHU*ZHD)
      MH02(4,5)=MH02(4,5)/dsqrt(ZHU*ZHD*ZS)
      MH02(5,4)=MH02(4,5)
      MH02(5,5)=MH02(5,5)/ZS

      MHC2=MHC2/(ZHU*ZHD)


c      II- Diagonalization of the neutral sector

      CALL DIAGN(5,MH02,VALPH,VECPH,1.d-10)
      CALL SORTNA(5,VALPH,VECPH)
      DO I=1,5
      MH0(I)=VALPH(I)
       DO J=1,5
      XH(I,J)=VECPH(J,I)
       ENDDO
      ENDDO

      DO I=1,5
      IF(MH0(I).le.0d0)THEN
       IFAIL=1
       GOTO 621
      ENDIF
      ENDDO


c            B: Pole corrections to the neutral states

      MW2=g2q/2d0*(vuq**2+vdq**2)
      MZ2=(g1q+g2q)/2d0*(vuq**2+vdq**2)
      Ytau=mtau/vdq
      sinbq=vuq/dsqrt(vuq**2+vdq**2)
      cosbq=vdq/dsqrt(vuq**2+vdq**2)

      MHTC(1)=MW2
      MHTC(2)=MHC
      DO I=1,5
       MHT2(I)=MH0T(I)
       XHT(I,1)=XH0(I,1)
       XHT(I,2)=XH0(I,2)
       XHT(I,3)=XH0(I,3)
       XHT(I,4)=XH0(I,4)*cosbq
       XHT(I,5)=XH0(I,4)*sinbq
       XHT(I,6)=XH0(I,5)
      ENDDO
      MHT2(6)=MZ2
      DO J=1,6
       XHT(6,J)=0d0
      ENDDO
      XHT(6,4)=-sinbq
      XHT(6,5)=cosbq
      DO I=1,5
       XHG(I,1)=XH(I,1)/dsqrt(ZHU)
       XHG(I,2)=XH(I,2)/dsqrt(ZHD)
       XHG(I,3)=XH(I,3)/dsqrt(ZS)
       XHG(I,4)=XH(I,4)*cosb/dsqrt(ZHU)
       XHG(I,5)=XH(I,4)*sinb/dsqrt(ZHD)
       XHG(I,6)=XH(I,5)/dsqrt(ZS)
      ENDDO


      DO I=1,5
       aux=0d0

c      I- SM Fermions

      aux=aux-3.d0*Ytq**2*((XHG(I,1)**2+XHG(I,4)**2)
     . *MH0(I)*(NMB0(MH0(I),mtopq**2,mtopq**2,QSTSB)
     .           -NMB0(muH2,mtopq**2,mtopq**2,QSTSB))
     .  -4.d0*mtopq**2*(NMB0(MH0(I),mtopq**2,mtopq**2,QSTSB)
     .         +fSF1(mtopq**2,mtopq**2,QSTSB))*XHG(I,1)**2)
     . -3.d0*Ybq**2*((XHG(I,2)**2+XHG(I,5)**2)
     . *MH0(I)*(NMB0(MH0(I),mbotq**2,mbotq**2,QSTSB)
     .           -NMB0(muH2,mbotq**2,mbotq**2,QSTSB))
     .  -4.d0*mbotq**2*(NMB0(MH0(I),mbotq**2,mbotq**2,QSTSB)
     .         +fSF1(mbotq**2,mbotq**2,QSTSB))*XHG(I,2)**2)
     . -(mtau/vdq)**2*((XHG(I,2)**2+XHG(I,5)**2)
     . *MH0(I)*(NMB0(MH0(I),mtau**2,mtau**2,QSTSB)
     .           -NMB0(muH2,mtau**2,mtau**2,QSTSB))
     .  -4.d0*mtau**2*(NMB0(MH0(I),mtau**2,mtau**2,QSTSB)
     .         +fSF1(mtau**2,mtau**2,QSTSB))*XHG(I,2)**2)

c      II- Gauginos/Higgsinos

      aux=aux-g1q/2d0*(MH0(I)*(NMB0(MH0(I),mur**2,M1r**2,QSTSB)
     .                          -NMB0(muH2,mur**2,M1r**2,QSTSB))
     .        -(mur**2+M1r**2)*(NMB0(MH0(I),mur**2,M1r**2,QSTSB)
     .                          +fSF1(mur**2,M1r**2,QSTSB)))
     .     *(XHG(I,1)**2+XHG(I,2)**2+XHG(I,4)**2+XHG(I,5)**2)
     . -3.d0*g2q/2d0*(MH0(I)*(NMB0(MH0(I),mur**2,M2r**2,QSTSB)
     .                         -NMB0(muH2,mur**2,M2r**2,QSTSB))
     .        -(mur**2+M2r**2)*(NMB0(MH0(I),mur**2,M2r**2,QSTSB)
     .                          +fSF1(mur**2,M2r**2,QSTSB)))
     .     *(XHG(I,1)**2+XHG(I,2)**2+XHG(I,4)**2+XHG(I,5)**2)
     . -l**2*(MH0(I)*(NMB0(MH0(I),mur**2,msi**2,QSTSB)
     .                     -NMB0(muH2,mur**2,msi**2,QSTSB))
     .        -(mur**2+msi**2)*
     .              (NMB0(MH0(I),mur**2,msi**2,QSTSB)
     .                          +fSF1(mur**2,msi**2,QSTSB)))
     .     *(XHG(I,1)**2+XHG(I,2)**2+XHG(I,4)**2+XHG(I,5)**2)
     . -2d0*(l**2*(MH0(I)*(NMB0(MH0(I),mur**2,mur**2,QSTSB)
     .                      -NMB0(muH2,mur**2,mur**2,QSTSB))
     .        -2d0*mur**2*(NMB0(MH0(I),mur**2,mur**2,QSTSB)
     .                          +fSF1(mur**2,mur**2,QSTSB)))
     .   +k**2*(MH0(I)*
     .       (NMB0(MH0(I),msi**2,msi**2,QSTSB)
     .        -NMB0(muH2,msi**2,msi**2,QSTSB))
     .        -2d0*msi**2*(NMB0(MH0(I),msi**2,msi**2,QSTSB)
     .              +fSF1(msi**2,msi**2,QSTSB)))
     .          )*(XHG(I,3)**2+XHG(I,6)**2)

c      III- Gauge

      aux=aux
     . +MH0(I)*(g2q*(NMB0(MH0(I),MW2,MW2,QSTSB)
     .                       -NMB0(muH2,MW2,MW2,QSTSB))
     .   +(g1q+g2q)/2d0*(NMB0(MH0(I),MZ2,MZ2,QSTSB)
     .                       -NMB0(muH2,MZ2,MZ2,QSTSB)))*sinb**2
     .             *(XHG(I,1)**2+XHG(I,4)**2)!+XHG(I,2)**2+XHG(I,5)**2
     . -MH0(I)*(g2q*NMB0(MH0(I),MW2,MW2,QSTSB)
     .   +(g1q+g2q)/2d0*NMB0(MH0(I),MZ2,MZ2,QSTSB))
     . *(0d0*(XHG(I,1)**2+XHG(I,4))-cosb**2*(XHG(I,2)**2+XHG(I,5)**2)
     .     -2d0*sinb*cosb*(XHG(I,1)*XHG(I,2)-XHG(I,4)*XHG(I,5)))
     . -3.d0*(g2q*MW2*(NMB0(MH0(I),MW2,MW2,QSTSB)+fSF1(MW2,MW2,QSTSB))
     .        +(g1q+g2q)/2d0*MZ2
     .         *(NMB0(MH0(I),MZ2,MZ2,QSTSB)+fSF1(MZ2,MZ2,QSTSB)))
     . *(sinb**2*(XHG(I,1)**2+XHG(I,4)**2)
     .    +cosb**2*(XHG(I,2)**2+XHG(I,5)**2)
     .    +2d0*sinb*cosb*(XHG(I,1)*XHG(I,2)-XHG(I,4)*XHG(I,5)))

      aux=aux+g2q*((MH0(I)+MHC-MW2/2d0)*NMB0(MH0(I),MHC,MW2,QSTSB)
     .                    +(MHC-MW2/2d0)*fSF1(MHC,MW2,QSTSB))
     . *(cosb**2*(XHG(I,1)**2+XHG(I,4)**2)
     .     +sinb**2*(XHG(I,2)**2+XHG(I,5)**2)
     .     -2d0*sinb*cosb*(XHG(I,1)*XHG(I,2)-XHG(I,4)*XHG(I,5)))

      DO J=1,5
       aux=aux+(g1q+g2q)/2d0*((MH0(I)+MH0T(J)-MZ2/2d0)
     .              *NMB0(MH0(I),MH0T(J),MZ2,QSTSB)
     .      +(MH0T(J)-MZ2/2d0)*fSF1(MH0T(J),MZ2,QSTSB))
     .  *(XHG(I,4)*XHT(J,1)-XHG(I,1)*XHT(J,4)
     .    +XHG(I,2)*XHT(J,5)-XHG(I,5)*XHT(J,2))**2
      ENDDO


c      IV- Higgs

      DO M=1,2
      DO N=1,2

      gRHSS=(XHG(I,1)*((g1q+g2q)/2d0*vuq*XC(M,1)*XC(N,1)
     .                 +(g2q-g1q)/2d0*vuq*XC(M,2)*XC(N,2)
     . -(l**2-g2q/2d0)*vdq*(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1)))
     .           +XHG(I,2)*((g2q-g1q)/2d0*vdq*XC(M,1)*XC(N,1)
     .                 +(g2q+g1q)/2d0*vdq*XC(M,2)*XC(N,2)
     . -(l**2-g2q/2d0)*vuq*(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1)))
     .           +XHG(I,3)*(l*(Alcos1+2d0*k/l*muq*DDCOS(Phi0)
     .                              +MUPQ*DDCOS(Phi01-phiP))
     .                *(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1))
     .               +2d0*l*muq*(XC(M,1)*XC(N,1)+XC(M,2)*XC(N,2)))
     .           +XHG(I,6)*(3.d0*k*muq*DDSIN(Phi0)
     .      +2d0*l*MUPQ*DDSIN(Phi01-phiP)
     .      +(M3HQ*DDSIN(phi3Q)+l*XIFQ*DDSIN(Phi01-phIF))*l/muq)
     .                *(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1))
     .                                           )/dsqrt(2d0)
      gIHSS=-(XHG(I,3)*(k*muq*DDSIN(Phi0)
     .      +(M3HQ*DDSIN(phi3Q)+l*XIFQ*DDSIN(Phi01-phIF))*l/muq)
     .          *(XC(M,1)*XC(N,2)-XC(M,2)*XC(N,1))
     .      +(l**2-g2q/2d0)*(vuq*XHG(I,5)+vdq*XHG(I,5))
     .          *(XC(M,1)*XC(N,2)-XC(M,2)*XC(N,1))
     .      +XHG(I,6)*l*(Alcos1-2d0*k/l*muq*DDCOS(Phi0)
     .                               +MUPQ*DDCOS(Phi01-PhiP))
     .          *(XC(M,1)*XC(N,2)-XC(M,2)*XC(N,1))
     .                                )/dsqrt(2d0)


       aux=aux-(gRHSS*gRHSS+gIHSS*gIHSS)*
     . (NMB0(MH0(I),MHTC(M),MHTC(N),QSTSB)+fSF1(MHTC(M),MHTC(N),QSTSB))

      ENDDO
      ENDDO


      DO M=1,6
      DO N=1,6

      PIS111=6.d0*XHG(I,1)*XHT(M,1)*XHT(N,1)
      PIS122=2d0*(XHG(I,1)*XHT(M,2)*XHT(N,2)
     . +XHG(I,2)*XHT(M,1)*XHT(N,2)+XHG(I,2)*XHT(M,2)*XHT(N,1))
      PIS133=2d0*(XHG(I,1)*XHT(M,3)*XHT(N,3)
     . +XHG(I,3)*XHT(M,1)*XHT(N,3)+XHG(I,3)*XHT(M,3)*XHT(N,1))
      PIS144=2d0*(XHG(I,1)*XHT(M,4)*XHT(N,4)
     . +XHG(I,4)*XHT(M,1)*XHT(N,4)+XHG(I,4)*XHT(M,4)*XHT(N,1))
      PIS155=2d0*(XHG(I,1)*XHT(M,5)*XHT(N,5)
     . +XHG(I,5)*XHT(M,1)*XHT(N,5)+XHG(I,5)*XHT(M,5)*XHT(N,1))
      PIS166=2d0*(XHG(I,1)*XHT(M,6)*XHT(N,6)
     . +XHG(I,6)*XHT(M,1)*XHT(N,6)+XHG(I,6)*XHT(M,6)*XHT(N,1))
      PIS211=2d0*(XHG(I,2)*XHT(M,1)*XHT(N,1)
     . +XHG(I,1)*XHT(M,2)*XHT(N,1)+XHG(I,1)*XHT(M,1)*XHT(N,2))
      PIS222=6.d0*XHG(I,2)*XHT(M,2)*XHT(N,2)
      PIS233=2d0*(XHG(I,2)*XHT(M,3)*XHT(N,3)
     . +XHG(I,3)*XHT(M,2)*XHT(N,3)+XHG(I,3)*XHT(M,3)*XHT(N,2))
      PIS244=2d0*(XHG(I,2)*XHT(M,4)*XHT(N,4)
     . +XHG(I,4)*XHT(M,2)*XHT(N,4)+XHG(I,4)*XHT(M,4)*XHT(N,2))
      PIS255=2d0*(XHG(I,2)*XHT(M,5)*XHT(N,5)
     . +XHG(I,5)*XHT(M,2)*XHT(N,5)+XHG(I,5)*XHT(M,5)*XHT(N,2))
      PIS266=2d0*(XHG(I,2)*XHT(M,6)*XHT(N,6)
     . +XHG(I,6)*XHT(M,2)*XHT(N,6)+XHG(I,6)*XHT(M,6)*XHT(N,2))
      PIS311=2d0*(XHG(I,3)*XHT(M,1)*XHT(N,1)
     . +XHG(I,1)*XHT(M,3)*XHT(N,1)+XHG(I,1)*XHT(M,1)*XHT(N,3))
      PIS322=2d0*(XHG(I,3)*XHT(M,2)*XHT(N,2)
     . +XHG(I,2)*XHT(M,3)*XHT(N,2)+XHG(I,2)*XHT(M,2)*XHT(N,3))
      PIS333=6.d0*XHG(I,3)*XHT(M,3)*XHT(N,3)
      PIS344=2d0*(XHG(I,3)*XHT(M,4)*XHT(N,4)
     . +XHG(I,4)*XHT(M,3)*XHT(N,4)+XHG(I,4)*XHT(M,4)*XHT(N,3))
      PIS355=2d0*(XHG(I,3)*XHT(M,5)*XHT(N,5)
     . +XHG(I,5)*XHT(M,3)*XHT(N,5)+XHG(I,5)*XHT(M,5)*XHT(N,3))
      PIS366=2d0*(XHG(I,3)*XHT(M,6)*XHT(N,6)
     . +XHG(I,6)*XHT(M,3)*XHT(N,6)+XHG(I,6)*XHT(M,6)*XHT(N,3))
      PIS312=XHG(I,3)*XHT(M,1)*XHT(N,2)
     . +XHG(I,3)*XHT(M,2)*XHT(N,1)+XHG(I,1)*XHT(M,3)*XHT(N,2)
     . +XHG(I,1)*XHT(M,2)*XHT(N,3)+XHG(I,2)*XHT(M,1)*XHT(N,3)
     . +XHG(I,2)*XHT(M,3)*XHT(N,1)
      PIS345=XHG(I,3)*XHT(M,4)*XHT(N,5)
     . +XHG(I,3)*XHT(M,5)*XHT(N,4)+XHG(I,4)*XHT(M,3)*XHT(N,5)
     . +XHG(I,4)*XHT(M,5)*XHT(N,3)+XHG(I,5)*XHT(M,4)*XHT(N,3)
     . +XHG(I,5)*XHT(M,3)*XHT(N,4)
      PIS156=XHG(I,1)*XHT(M,5)*XHT(N,6)
     . +XHG(I,1)*XHT(M,6)*XHT(N,5)+XHG(I,5)*XHT(M,1)*XHT(N,6)
     . +XHG(I,5)*XHT(M,6)*XHT(N,1)+XHG(I,6)*XHT(M,5)*XHT(N,1)
     . +XHG(I,6)*XHT(M,1)*XHT(N,5)
      PIS256=XHG(I,2)*XHT(M,5)*XHT(N,6)
     . +XHG(I,2)*XHT(M,6)*XHT(N,5)+XHG(I,5)*XHT(M,2)*XHT(N,6)
     . +XHG(I,5)*XHT(M,6)*XHT(N,2)+XHG(I,6)*XHT(M,5)*XHT(N,2)
     . +XHG(I,6)*XHT(M,2)*XHT(N,5)
      PIS346=XHG(I,3)*XHT(M,4)*XHT(N,6)
     . +XHG(I,3)*XHT(M,6)*XHT(N,4)+XHG(I,4)*XHT(M,3)*XHT(N,6)
     . +XHG(I,4)*XHT(M,6)*XHT(N,3)+XHG(I,6)*XHT(M,4)*XHT(N,3)
     . +XHG(I,6)*XHT(M,3)*XHT(N,4)
      PIS356=XHG(I,3)*XHT(M,5)*XHT(N,6)
     . +XHG(I,3)*XHT(M,6)*XHT(N,5)+XHG(I,5)*XHT(M,3)*XHT(N,6)
     . +XHG(I,5)*XHT(M,6)*XHT(N,3)+XHG(I,6)*XHT(M,5)*XHT(N,3)
     . +XHG(I,6)*XHT(M,3)*XHT(N,5)
      PIS513=XHG(I,5)*XHT(M,1)*XHT(N,3)
     . +XHG(I,5)*XHT(M,3)*XHT(N,1)+XHG(I,1)*XHT(M,5)*XHT(N,3)
     . +XHG(I,1)*XHT(M,3)*XHT(N,5)+XHG(I,3)*XHT(M,1)*XHT(N,5)
     . +XHG(I,3)*XHT(M,5)*XHT(N,1)
      PIS423=XHG(I,4)*XHT(M,2)*XHT(N,3)
     . +XHG(I,4)*XHT(M,3)*XHT(N,2)+XHG(I,2)*XHT(M,4)*XHT(N,3)
     . +XHG(I,2)*XHT(M,3)*XHT(N,4)+XHG(I,3)*XHT(M,2)*XHT(N,4)
     . +XHG(I,3)*XHT(M,4)*XHT(N,2)
      PIS612=XHG(I,6)*XHT(M,1)*XHT(N,2)
     . +XHG(I,6)*XHT(M,2)*XHT(N,1)+XHG(I,1)*XHT(M,6)*XHT(N,2)
     . +XHG(I,1)*XHT(M,2)*XHT(N,6)+XHG(I,6)*XHT(M,1)*XHT(N,2)
     . +XHG(I,6)*XHT(M,2)*XHT(N,1)
      PIS456=XHG(I,4)*XHT(M,5)*XHT(N,6)
     . +XHG(I,4)*XHT(M,6)*XHT(N,5)+XHG(I,5)*XHT(M,4)*XHT(N,6)
     . +XHG(I,5)*XHT(M,6)*XHT(N,4)+XHG(I,6)*XHT(M,5)*XHT(N,4)
     . +XHG(I,6)*XHT(M,4)*XHT(N,5)
      PIS433=2d0*(XHG(I,4)*XHT(M,3)*XHT(N,3)
     . +XHG(I,3)*XHT(M,4)*XHT(N,3)+XHG(I,3)*XHT(M,3)*XHT(N,4))
      PIS613=XHG(I,6)*XHT(M,1)*XHT(N,3)
     . +XHG(I,6)*XHT(M,3)*XHT(N,1)+XHG(I,1)*XHT(M,6)*XHT(N,3)
     . +XHG(I,1)*XHT(M,3)*XHT(N,6)+XHG(I,3)*XHT(M,1)*XHT(N,6)
     . +XHG(I,3)*XHT(M,6)*XHT(N,1)
      PIS466=2d0*(XHG(I,4)*XHT(M,6)*XHT(N,6)
     . +XHG(I,6)*XHT(M,4)*XHT(N,6)+XHG(I,6)*XHT(M,6)*XHT(N,4))
      PIS533=2d0*(XHG(I,5)*XHT(M,3)*XHT(N,3)
     . +XHG(I,3)*XHT(M,5)*XHT(N,3)+XHG(I,3)*XHT(M,3)*XHT(N,5))
      PIS623=XHG(I,6)*XHT(M,2)*XHT(N,3)
     . +XHG(I,6)*XHT(M,3)*XHT(N,2)+XHG(I,2)*XHT(M,6)*XHT(N,3)
     . +XHG(I,2)*XHT(M,3)*XHT(N,6)+XHG(I,3)*XHT(M,2)*XHT(N,6)
     . +XHG(I,3)*XHT(M,6)*XHT(N,2)
      PIS566=2d0*(XHG(I,5)*XHT(M,6)*XHT(N,6)
     . +XHG(I,6)*XHT(M,5)*XHT(N,6)+XHG(I,6)*XHT(M,6)*XHT(N,5))
      PIS666=6.d0*XHG(I,6)*XHT(M,6)*XHT(N,6)
      PIS633=2d0*(XHG(I,6)*XHT(M,3)*XHT(N,3)
     . +XHG(I,3)*XHT(M,6)*XHT(N,3)+XHG(I,3)*XHT(M,3)*XHT(N,6))
      PIS246=XHG(I,2)*XHT(M,4)*XHT(N,6)
     . +XHG(I,2)*XHT(M,6)*XHT(N,4)+XHG(I,4)*XHT(M,2)*XHT(N,6)
     . +XHG(I,4)*XHT(M,6)*XHT(N,2)+XHG(I,6)*XHT(M,4)*XHT(N,2)
     . +XHG(I,6)*XHT(M,2)*XHT(N,4)

        IF(k.ne.0d0)then
       gHSS=((g1q+g2q)/4.d0*(vuq*(PIS111+PIS144-PIS122-PIS155)
     .                      +vdq*(PIS222+PIS255-PIS211-PIS244))
     . +l*muq*(PIS311+PIS344+PIS322+PIS355)
     . +l**2*vuq*(PIS122+PIS155+PIS133+PIS166)
     . +l**2*vdq*(PIS211+PIS244+PIS233+PIS266)
     . -l*Alcos1*(PIS312-PIS345-PIS156-PIS246)
     . +k/3.d0*Akcos2*(PIS333-3.d0*PIS366)
     . +2d0*k**2*muq/l*(PIS333+PIS366)
     . -k*DDCOS(Phi0)*(2d0*muq*(PIS312-PIS345+PIS156+PIS246)
     .  +l*vdq*(PIS133-PIS166+2d0*PIS346)
     .  +l*vuq*(PIS233-PIS266+2d0*PIS356))
     . +k*DDSIN(Phi0)*(muq*(PIS513+PIS423-3.d0*PIS612+3.d0*PIS456)
     .  +l*vdq*(PIS433-2d0*PIS613-PIS466)
     .  +l*vuq*(PIS533-2d0*PIS623-PIS566)
     . +l*vuq*vdq/muq*(3.d0*PIS633-PIS666))
     . -l*MUPQ*DDCOS(Phi01-phiP)*(PIS312-PIS345+PIS156+PIS246)
     . +k*MUPQ*DDCOS(Phi02-phiP)*(PIS333+PIS366)
     . -l/muq*(M3HQ*DDSIN(phi3Q)+l*XIFQ*DDSIN(Phi01-phIF))
     .  *(PIS612+PIS423+PIS513-PIS456
     .          +vuq*vdq*l**2/muq**2*(PIS633-PIS666/3d0))
     . -2d0*l*MUPQ*DDSIN(Phi01-phiP)*(PIS612
     .    +vuq*vdq*l**2/muq**2*(PIS666/3d0-PIS633))
     . -4d0/3d0*k*MUPQ*DDSIN(Phi02-phiP)*PIS666
     . -l/muq*(l/muq*(XISQ*DDSIN(phiSQ)+XIFQ*MUPQ*DDSIN(PhiP-phIF))
     .      +MSPQ*DDSIN(PhiSPQ))*(PIS666/3d0-PIS633))/dsqrt(2d0)
        else
       gHSS=((g1q+g2q)/4.d0*(vuq*(PIS111+PIS144-PIS122-PIS155)
     .                      +vdq*(PIS222+PIS255-PIS211-PIS244))
     . +l*muq*(PIS311+PIS344+PIS322+PIS355)
     . +l**2*vuq*(PIS122+PIS155+PIS133+PIS166)
     . +l**2*vdq*(PIS211+PIS244+PIS233+PIS266)
     . -l*Alcos1*(PIS312-PIS345-PIS156-PIS246)
     . -l*MUPQ*DDCOS(Phi01-phiP)*(PIS312-PIS345+PIS156+PIS246)
     . -l/muq*(M3HQ*DDSIN(phi3Q)+l*XIFQ*DDSIN(Phi01-phIF))
     .  *(PIS612+PIS423+PIS513-PIS456)
     . -2d0*l*MUPQ*DDSIN(Phi01-phiP)*PIS612)/dsqrt(2d0)
        endif

       aux=aux-gHSS**2*(NMB0(MH0(I),MHT2(M),MHT2(N),QSTSB)
     .                     +fSF1(MHT2(M),MHT2(N),QSTSB))/2d0

      ENDDO
      ENDDO

c      V- Sfermions

      DO M=1,2
      DO N=1,2

       gRHSS=dsqrt(2d0)*(UT(M,1,1)*UT(N,1,1)+UT(M,1,2)*UT(N,1,2))
     . *(Ytq**2*vuq*XHG(I,1)+(g1q/3.d0-g2q)/4.d0
     .                               *(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +dsqrt(2d0)*(UT(M,2,1)*UT(N,2,1)+UT(M,2,2)*UT(N,2,2))
     . *(Ytq**2*vuq*XHG(I,1)-g1q/3.d0*(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +Ytq/dsqrt(2d0)*(AT*(DDCOS(PhiAT)*XHG(I,1)
     .   -DDSIN(PhiAT)*XHG(I,4))
     .   -DDCOS(Phi01)*(muq*XHG(I,2)+l*vdq*XHG(I,3))
     .   +DDSIN(Phi01)*(muq*XHG(I,5)+l*vdq*XHG(I,6)))
     .  *(UT(M,1,1)*UT(N,2,1)+UT(M,2,1)*UT(N,1,1)
     .      +UT(M,1,2)*UT(N,2,2)+UT(M,2,2)*UT(N,1,2))
     . +Ytq/dsqrt(2d0)*(AT*(DDSIN(PhiAT)*XHG(I,1)
     .   +DDCOS(PhiAT)*XHG(I,4))
     .   +DDSIN(Phi01)*(muq*XHG(I,2)+l*vdq*XHG(I,3))
     .   +DDCOS(Phi01)*(muq*XHG(I,5)+l*vdq*XHG(I,6)))
     .  *(UT(M,2,1)*UT(N,1,2)-UT(M,2,2)*UT(N,1,1)
     .      +UT(M,1,2)*UT(N,2,1)-UT(M,1,1)*UT(N,2,2))

       gIHSS=dsqrt(2d0)*(UT(M,1,2)*UT(N,1,1)-UT(M,1,1)*UT(N,1,2))
     . *(Ytq**2*vuq*XHG(I,1)+(g1q/3.d0-g2q)/4.d0
     .                               *(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +dsqrt(2d0)*(UT(M,2,2)*UT(N,2,1)-UT(M,2,1)*UT(N,2,2))
     . *(Ytq**2*vuq*XHG(I,1)-g1q/3.d0*(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +Ytq/dsqrt(2d0)*(AT*(DDCOS(PhiAT)*XHG(I,1)
     .   -DDSIN(PhiAT)*XHG(I,4))
     .   -DDCOS(Phi01)*(muq*XHG(I,2)+l*vdq*XHG(I,3))
     .   +DDSIN(Phi01)*(muq*XHG(I,5)+l*vdq*XHG(I,6)))
     .  *(UT(M,2,2)*UT(N,1,1)-UT(M,2,1)*UT(N,1,2)
     .      +UT(M,1,2)*UT(N,2,1)-UT(M,1,1)*UT(N,2,2))
     . +Ytq/dsqrt(2d0)*(AT*(DDSIN(PhiAT)*XHG(I,1)
     .   +DDCOS(PhiAT)*XHG(I,4))
     .   +DDSIN(Phi01)*(muq*XHG(I,2)+l*vdq*XHG(I,3))
     .   +DDCOS(Phi01)*(muq*XHG(I,5)+l*vdq*XHG(I,6)))
     .  *(UT(M,2,1)*UT(N,1,1)-UT(M,1,1)*UT(N,2,1)
     .      +UT(M,2,2)*UT(N,1,2)-UT(M,1,2)*UT(N,2,2))

      aux=aux-3.d0*(gRHSS**2+gIHSS**2)
     . *(NMB0(MH0(I),MST2(M),MST2(N),QSTSB)+fSF1(MST2(M),MST2(N),QSTSB))

      gRHSTST(I,M,N)=gRHSS
      gIHSTST(I,M,N)=gIHSS

       gRHSS=dsqrt(2d0)*(UB(M,1,1)*UB(N,1,1)+UB(M,1,2)*UB(N,1,2))
     . *(Ybq**2*vdq*XHG(I,2)+(g1q/3.d0+g2q)/4.d0
     .                               *(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +dsqrt(2d0)*(UB(M,2,1)*UB(N,2,1)+UB(M,2,2)*UB(N,2,2))
     . *(Ybq**2*vdq*XHG(I,2)+g1q/6.d0*(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +Ybq/dsqrt(2d0)*(AB*(DDCOS(PhiAB)*XHG(I,2)
     .   -DDSIN(PhiAB)*XHG(I,5))
     .   -DDCOS(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .   +DDSIN(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(UB(M,1,1)*UB(N,2,1)+UB(M,2,1)*UB(N,1,1)
     .      +UB(M,1,2)*UB(N,2,2)+UB(M,2,2)*UB(N,1,2))
     . +Ybq/dsqrt(2d0)*(AB*(DDSIN(PhiAB)*XHG(I,2)
     .   +DDCOS(PhiAB)*XHG(I,5))
     .   +DDSIN(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .   +DDCOS(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(UB(M,2,1)*UB(N,1,2)-UB(M,2,2)*UB(N,1,1)
     .      +UB(M,1,2)*UB(N,2,1)-UB(M,1,1)*UB(N,2,2))

       gIHSS=dsqrt(2d0)*(UB(M,1,2)*UB(N,1,1)-UB(M,1,1)*UB(N,1,2))
     . *(Ybq**2*vdq*XHG(I,2)+(g1q/3.d0+g2q)/4.d0
     .                               *(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +dsqrt(2d0)*(UB(M,2,2)*UB(N,2,1)-UB(M,2,1)*UB(N,2,2))
     . *(Ybq**2*vdq*XHG(I,2)+g1q/6.d0*(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +Ybq/dsqrt(2d0)*(AB*(DDCOS(PhiAB)*XHG(I,2)
     .   -DDSIN(PhiAB)*XHG(I,5))
     .   -DDCOS(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .   +DDSIN(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(UB(M,2,2)*UB(N,1,1)-UB(M,2,1)*UB(N,1,2)
     .      +UB(M,1,2)*UB(N,2,1)-UB(M,1,1)*UB(N,2,2))
     . +Ybq/dsqrt(2d0)*(AB*(DDSIN(PhiAB)*XHG(I,2)
     .   +DDCOS(PhiAB)*XHG(I,5))
     .   +DDSIN(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .   +DDCOS(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(UB(M,2,1)*UB(N,1,1)-UB(M,1,1)*UB(N,2,1)
     .      +UB(M,2,2)*UB(N,1,2)-UB(M,1,2)*UB(N,2,2))

      aux=aux-3.d0*(gRHSS**2+gIHSS**2)
     . *(NMB0(MH0(I),MSB2(M),MSB2(N),QSTSB)+fSF1(MSB2(M),MSB2(N),QSTSB))

      gRHSBSB(I,M,N)=gRHSS
      gIHSBSB(I,M,N)=gIHSS

       gRHSS=dsqrt(2d0)*
     .        (UTAU(M,1,1)*UTAU(N,1,1)+UTAU(M,1,2)*UTAU(N,1,2))
     . *(Ytau**2*vdq*XHG(I,2)+(-g1q+g2q)/4.d0
     .                               *(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +dsqrt(2d0)*
     .        (UTAU(M,2,1)*UTAU(N,2,1)+UTAU(M,2,2)*UTAU(N,2,2))
     . *(Ytau**2*vdq*XHG(I,2)+g1q/2d0*(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +Ytau/dsqrt(2d0)*
     .            (ATAU*(DDCOS(PhiATAU)*XHG(I,2)
     .   -DDSIN(PhiATAU)*XHG(I,5))
     .   -DDCOS(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .   +DDSIN(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(UTAU(M,1,1)*UTAU(N,2,1)+UTAU(M,2,1)*UTAU(N,1,1)
     .      +UTAU(M,1,2)*UTAU(N,2,2)+UTAU(M,2,2)*UTAU(N,1,2))
     . +Ytau/dsqrt(2d0)*
     .            (ATAU*(DDSIN(PhiATAU)*XHG(I,2)
     .   +DDCOS(PhiATAU)*XHG(I,5))
     .   +DDSIN(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .   +DDCOS(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(UTAU(M,2,1)*UTAU(N,1,2)-UTAU(M,2,2)*UTAU(N,1,1)
     .      +UTAU(M,1,2)*UTAU(N,2,1)-UTAU(M,1,1)*UTAU(N,2,2))

       gIHSS=dsqrt(2d0)*
     .              (UTAU(M,1,2)*UTAU(N,1,1)-UTAU(M,1,1)*UTAU(N,1,2))
     . *(Ytau**2*vdq*XHG(I,2)+(-g1q+g2q)/4.d0
     .                               *(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +dsqrt(2d0)*(UTAU(M,2,2)*UTAU(N,2,1)-UTAU(M,2,1)*UTAU(N,2,2))
     . *(Ytau**2*vdq*XHG(I,2)+g1q/2d0*(vuq*XHG(I,1)-vdq*XHG(I,2)))
     . +Ytau/dsqrt(2d0)*
     .          (ATAU*(DDCOS(PhiATAU)*XHG(I,2)-DDSIN(PhiATAU)*XHG(I,5))
     .   -DDCOS(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .   +DDSIN(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(UTAU(M,2,2)*UTAU(N,1,1)-UTAU(M,2,1)*UTAU(N,1,2)
     .      +UTAU(M,1,2)*UTAU(N,2,1)-UTAU(M,1,1)*UTAU(N,2,2))
     . +Ytau/dsqrt(2d0)*
     .          (ATAU*(DDSIN(PhiATAU)*XHG(I,2)+DDCOS(PhiATAU)*XHG(I,5))
     .   +DDSIN(Phi01)*(muq*XHG(I,1)+l*vuq*XHG(I,3))
     .   +DDCOS(Phi01)*(muq*XHG(I,4)+l*vuq*XHG(I,6)))
     .  *(UTAU(M,2,1)*UTAU(N,1,1)-UTAU(M,1,1)*UTAU(N,2,1)
     .      +UTAU(M,2,2)*UTAU(N,1,2)-UTAU(M,1,2)*UTAU(N,2,2))

      aux=aux-(gRHSS**2+gIHSS**2)
     . *(NMB0(MH0(I),MSL2(M),MSL2(N),QSTSB)+fSF1(MSL2(M),MSL2(N),QSTSB))

      gRHSLSL(I,M,N)=gRHSS
      gIHSLSL(I,M,N)=gIHSS

      gRHSS=dsqrt(2d0)*DELT(M,1)*DELT(N,1)*
     .          (g1q/3.d0-g2q)/4.d0*(vuq*XHG(I,1)-vdq*XHG(I,2))
     .     +dsqrt(2d0)*DELT(M,2)*DELT(N,2)*
     .          (-g1q/3.d0)*(vuq*XHG(I,1)-vdq*XHG(I,2))

      aux=aux-6.d0*gRHSS**2
     . *(NMB0(MH0(I),MSU2(M),MSU2(N),QSTSB)+fSF1(MSU2(M),MSU2(N),QSTSB))

      gRHSUSU(I,M,N)=gRHSS

      gRHSS=dsqrt(2d0)*DELT(M,1)*DELT(N,1)*
     .          (g1q/3.d0+g2q)/4.d0*(vuq*XHG(I,1)-vdq*XHG(I,2))
     .     +dsqrt(2d0)*DELT(M,2)*DELT(N,2)*
     .           g1q/6.d0*(vuq*XHG(I,1)-vdq*XHG(I,2))

      aux=aux-6.d0*gRHSS**2
     . *(NMB0(MH0(I),MSD2(M),MSD2(N),QSTSB)+fSF1(MSD2(M),MSD2(N),QSTSB))

      gRHSDSD(I,M,N)=gRHSS

      gRHSS=dsqrt(2d0)*DELT(M,1)*DELT(N,1)*
     .          (-g1q+g2q)/4.d0*(vuq*XHG(I,1)-vdq*XHG(I,2))
     .     +dsqrt(2d0)*DELT(M,2)*DELT(N,2)*
     .           g1q/2d0*(vuq*XHG(I,1)-vdq*XHG(I,2))

      aux=aux-2d0*gRHSS**2
     . *(NMB0(MH0(I),MSE2(M),MSE2(N),QSTSB)+fSF1(MSE2(M),MSE2(N),QSTSB))

      gRHSESE(I,M,N)=gRHSS

      ENDDO
      ENDDO

      gRHSS=-dsqrt(2d0)*(g1q+g2q)/4.d0*(vuq*XHG(I,1)-vdq*XHG(I,2))

      aux=aux-gRHSS**2
     . *(NMB0(MH0(I),MSNT2,MSNT2,QSTSB)+fSF1(MSNT2,MSNT2,QSTSB))

      aux=aux-2d0*gRHSS**2
     . *(NMB0(MH0(I),MSNE2,MSNE2,QSTSB)+fSF1(MSNE2,MSNE2,QSTSB))

      gRHSNSN(I)=gRHSS

c      VI- Inclusion of the pole corrections

      MH0(I)=MH0(I)
     .    +aux/16.d0/Pi**2

      IF(MH0(I).le.0d0)THEN
       IFAIL=1
       GOTO 621
      ENDIF

      ENDDO


c            C: Pole corrections to the charged state

      aux=0d0

c      I- SM fermions

      aux=aux-3.d0*(Ytq**2*cosb**2+Ybq**2*sinb**2)
     . *(MHC2*NMB0(MHC2,mtopq**2,mbotq**2,QSTSB)
     .  -(mtopq**2+mbotq**2)*(NMB0(MHC2,mtopq**2,mbotq**2,QSTSB)
     .                          +fSF1(mtopq**2,mbotq**2,QSTSB)))
     . +3.d0*MHC2*(Ytq**2*cosb**2*NMB0(muH2,mtopq**2,mtopq**2,QSTSB)
     .            +Ybq**2*sinb**2*NMB0(muH2,mbotq**2,mbotq**2,QSTSB))
     . +4.d0*mtopq*mbotq*Ytq*Ybq*sinb*cosb*
     .                             (NMB0(MHC2,mtopq**2,mbotq**2,QSTSB)
     .                                  +fSF1(mtopq**2,mbotq**2,QSTSB))
     . -Ytau**2*sinb**2*(MHC2*(NMB0(MHC2,0d0,mtau**2,QSTSB)
     .                        -NMB0(muH2,mtau**2,mtau**2,QSTSB))
     .  -mtau**2*(NMB0(MHC2,0d0,mtau**2,QSTSB)+fSF1(0d0,mtau**2,QSTSB)))


c      II- Higgsinos/gauginos

      aux=aux-g1q/2d0*(MHC2*(NMB0(MHC2,mur**2,M1r**2,QSTSB)
     .                       -NMB0(muH2,mur**2,M1r**2,QSTSB))
     .        -(mur**2+M1r**2)*(NMB0(MHC2,mur**2,M1r**2,QSTSB)
     .                          +fSF1(mur**2,M1r**2,QSTSB)))
     .        -3.d0*g2q/2d0*(MHC2*(NMB0(MHC2,mur**2,M2r**2,QSTSB)
     .                       -NMB0(muH2,mur**2,M2r**2,QSTSB))
     .        -(mur**2+M2r**2)*(NMB0(MHC2,mur**2,M2r**2,QSTSB)
     .                          +fSF1(mur**2,M2r**2,QSTSB)))
     .        -l**2*(MHC2*(NMB0(MHC2,mur**2,msi**2,QSTSB)
     .                    -NMB0(muH2,mur**2,msi**2,QSTSB))
     .        -(mur**2+msi**2)
     .                 *(NMB0(MHC2,mur**2,msi**2,QSTSB)
     .                          +fSF1(mur**2,msi**2,QSTSB)))

c      III- Gauge

      aux=aux-MHC2*(g2q*NMB0(muH2,MW2,MW2,QSTSB)
     .   +(g1q+g2q)/2d0*NMB0(muH2,MZ2,MZ2,QSTSB))*sinb**2*cosb**2


c      IV- Higgs

      DO M=1,5
       aux=aux+g2q/4.d0*(cosbq**2*(XHT(M,1)**2+XHT(M,4)**2)
     .          +sinbq**2*(XHT(M,2)**2+XHT(M,5)**2)
     . -2d0*sinbq*cosbq*(XHT(M,1)*XHT(M,2)-XHT(M,4)*XHT(M,5)))
     .   *(2d0*MHC2*NMB0(MHC2,MW2,MH0T(M),QSTSB)+(2d0*MHT2(M)-MW2)
     .      *(NMB0(MHC2,MW2,MH0T(M),QSTSB)+fSF1(MW2,MH0T(M),QSTSB)))
      ENDDO

      aux=aux+g1q*g2q/(g1q+g2q)*(2d0*MHC2*NMB0(MHC2,0d0,MHC,QSTSB)
     . +(fSF0(MHC,QSTSB)-NMB0(MHC2,0d0,MHC,QSTSB)))

      aux=aux+(g1q-g2q)**2/(g1q+g2q)/4.d0*(
     . 2d0*MHC2*NMB0(MHC2,MZ2,MHC,QSTSB)+(2d0*MHC-MZ2)*
     . (NMB0(MHC2,MZ2,MHC,QSTSB)+fSF1(MZ2,MHC,QSTSB)))

      DO M=1,2
      DO N=1,5

      gRHSS=(XHT(N,1)*((g1q+g2q)/2d0*vuq*XC(M,1)*cosb
     .                 +(g2q-g1q)/2d0*vuq*XC(M,2)*sinb
     . -(l**2-g2q/2d0)*vdq*(XC(M,1)*sinb+XC(M,2)*cosb))
     .           +XHT(N,2)*((g2q-g1q)/2d0*vdq*XC(M,1)*cosb
     .                 +(g2q+g1q)/2d0*vdq*XC(M,2)*sinb
     . -(l**2-g2q/2d0)*vuq*(XC(M,1)*sinb+XC(M,2)*cosb))
     .           +XHT(N,3)*(l*(Alcos1+2d0*k/l*muq*DDCOS(Phi0)
     .                              +MUPQ*DDCOS(Phi01-phiP))
     .               +2d0*l*muq*(XC(M,1)*cosb+XC(M,2)*sinb)
     .                *(XC(M,1)*sinb+XC(M,2)*cosb))
     .           +XHT(N,5)*(3.d0*k*muq*DDSIN(Phi0)
     .      +2d0*l*MUPQ*DDSIN(Phi01-phiP)
     .      +(M3HQ*DDSIN(phi3Q)+l*XIFQ*DDSIN(Phi01-phIF))*l/muq)
     .                *(XC(M,1)*sinb+XC(M,2)*cosb)
     .                                           )/dsqrt(2d0)
      gIHSS=-(XHT(N,3)*k*muq*DDSIN(Phi0)
     .          *(XC(M,1)*sinb-XC(M,2)*cosb)
     .      +XHT(N,4)*(l**2-g2q/2d0)*(vuq*sinbq+vdq*cosbq)
     .          *(XC(M,1)*sinb-XC(M,2)*cosb)
     .      +XHT(N,5)*l*(Alcos1-2d0*k/l*muq*DDCOS(Phi0)
     .                               +MUPQ*DDCOS(Phi01-PhiP))
     .          *(XC(M,1)*sinb-XC(M,2)*cosb)
     .                                )/dsqrt(2d0)

      aux=aux-(gRHSS**2+gIHSS**2)*
     . (NMB0(MHC2,MH0(N),MHTC(M),QSTSB)+fSF1(MH0(N),MHTC(M),QSTSB))

      ENDDO
      ENDDO


c      V- Sfermions

      DO M=1,2
      DO N=1,2

       gRHSS=-(Ytq**2+Ybq**2-g2q)*dsqrt(vuq**2+vdq**2)*sinb*cosb
     .               *(UT(M,1,1)*UB(N,1,1)+UT(M,1,2)*UB(N,1,2))
     . -Ytq*Ybq*dsqrt(vuq**2+vdq**2)
     .               *(UT(M,2,1)*UB(N,2,1)+UT(M,2,2)*UB(N,2,2))
     . -Ytq*(AT*DDCOS(PhiAT)*cosb+muq*DDCOS(Phi01)*sinb)
     .    *(UT(M,2,1)*UB(N,1,1)+UT(M,2,2)*UB(N,1,2))
     . +Ytq*(AT*DDSIN(PhiAT)*cosb-muq*DDSIN(Phi01)*sinb)
     .    *(UT(M,2,2)*UB(N,1,1)-UT(M,2,1)*UB(N,1,2))
     . -Ybq*(AB*DDCOS(PhiAB)*sinb+muq*DDCOS(Phi01)*cosb)
     .    *(UT(M,1,1)*UB(N,2,1)+UT(M,1,2)*UB(N,2,2))
     . +Ybq*(AB*DDSIN(PhiAB)*sinb-muq*DDSIN(Phi01)*cosb)
     .    *(UT(M,1,1)*UB(N,2,2)-UT(M,1,2)*UB(N,2,1))

       gIHSS=-(Ytq**2+Ybq**2-g2q)*dsqrt(vuq**2+vdq**2)*sinb*cosb
     .               *(UT(M,1,2)*UB(N,1,1)-UT(M,1,1)*UB(N,1,2))
     . -Ytq*Ybq*dsqrt(vuq**2+vdq**2)
     .               *(UT(M,2,2)*UB(N,2,1)-UT(M,2,1)*UB(N,2,2))
     . -Ytq*(AT*DDCOS(PhiAT)*cosb+muq*DDCOS(Phi01)*sinb)
     .    *(UT(M,2,2)*UB(N,1,1)-UT(M,2,1)*UB(N,1,2))
     . -Ytq*(AT*DDSIN(PhiAT)*cosb-muq*DDSIN(Phi01)*sinb)
     .    *(UT(M,2,1)*UB(N,1,1)+UT(M,2,2)*UB(N,1,2))
     . +Ybq*(AB*DDCOS(PhiAB)*sinb+muq*DDCOS(Phi01)*cosb)
     .    *(UT(M,1,1)*UB(N,2,2)-UT(M,1,2)*UB(N,2,1))
     . +Ybq*(AB*DDSIN(PhiAB)*sinb-muq*DDSIN(Phi01)*cosb)
     .    *(UT(M,1,1)*UB(N,2,1)+UT(M,1,2)*UB(N,2,2))

       aux=aux-3.d0*(gRHSS**2+gIHSS**2)*
     .           (NMB0(MHC2,MST2(M),MSB2(N),QSTSB)
     .                            +fSF1(MST2(M),MSB2(N),QSTSB))

      gRHCSTSB(M,N)=gRHSS
      gIHCSTSB(M,N)=gIHSS

       gRHSS=g2q*dsqrt(vuq**2+vdq**2)*sinb*cosb*DELT(M,1)*DELT(N,1)

       aux=aux-6.d0*gRHSS**2*
     .           (NMB0(MHC2,MSU2(M),MSD2(N),QSTSB)
     .                            +fSF1(MSU2(M),MSD2(N),QSTSB))

      gRHCSUSD(M,N)=gRHSS

      ENDDO

       gRHSS=-(Ytau**2-g2q)*dsqrt(vuq**2+vdq**2)*sinb*cosb
     .                                               *UTAU(M,1,1)
     . -Ytau*(ATAU*DDCOS(PhiATAU)*sinb+muq*DDCOS(Phi01)*cosb)
     .                                               *UTAU(M,2,1)
     . +Ytau*(ATAU*DDSIN(PhiATAU)*sinb-muq*DDSIN(Phi01)*cosb)
     .                                               *UTAU(M,2,2)

       gIHSS=(Ytau**2-g2q)*dsqrt(vuq**2+vdq**2)*sinb*cosb
     .                                               *UTAU(M,1,2)
     . +Ytau*(ATAU*DDCOS(PhiATAU)*sinb+muq*DDCOS(Phi01)*cosb)
     .                                               *UTAU(M,2,2)
     . +Ytau*(ATAU*DDSIN(PhiATAU)*sinb-muq*DDCOS(Phi01)*cosb)
     .                                               *UTAU(M,2,1)

       aux=aux-(gRHSS**2+gIHSS**2)*
     .           (NMB0(MHC2,MSNT2,MSL2(M),QSTSB)
     .                            +fSF1(MSNT2,MSL2(M),QSTSB))

      gRHCSNSL(M)=gRHSS
      gIHCSNSL(M)=gIHSS

       gRHSS=g2q*dsqrt(vuq**2+vdq**2)*sinb*cosb*DELT(M,1)

       aux=aux-2d0*gRHSS**2*
     .           (NMB0(MHC2,MSNE2,MSE2(M),QSTSB)
     .                            +fSF1(MSNE2,MSE2(M),QSTSB))

      gRHCSNSE(M)=gRHSS

      ENDDO


c      VI- Inclusion of the pole corrections

      MHC2=MHC2+aux/16.d0/Pi**2


      IF(MHC2.le.0d0)THEN
       IFAIL=4
       GOTO 621
      ENDIF
      MHC=MHC2
      MH0T=MH0
      XH0=XH

c      print*,dsqrt(MH0(1)),dsqrt(MH0(2)),dsqrt(MH0(3)),dsqrt(MH0(4)),
c     c dsqrt(MH0(5)),dsqrt(MHC2)


 621      RETURN
      END

************************************************************************************************

      DOUBLE PRECISION function Fsf0(x,z)
      
c            ->Fsf0(m^2,Q^2)

      IMPLICIT NONE
      DOUBLE PRECISION x,z,aux

      aux=x*(dlog(x/z)-1d0)

      Fsf0=aux

      RETURN
      END

************************************************************************************************

      DOUBLE PRECISION function Fsf2(x,y)
      
c            ->Fsf2(m1^2,m2^2)

      IMPLICIT NONE
      DOUBLE PRECISION x,y,aux

      IF(dabs(x-y).ge.1.d-10)THEN
       aux=(x+y)/(y-x)*dlog(x/y)+2d0
      ELSE
       aux=0d0
      ENDIF

      Fsf2=aux

      RETURN
      END

************************************************************************************************

      DOUBLE PRECISION function Fsf3(x,z)
      
c            ->Fsf3(m1^2,m2^2)

      IMPLICIT NONE
      DOUBLE PRECISION x,z,aux

      IF(dabs(x-z).ge.1.d-10)THEN
       aux=(x*dlog(x/z)+z-x)/(z-x)**2
      ELSE
       aux=1d0/2d0/x
      ENDIF

      Fsf3=aux

      RETURN
      END

************************************************************************************************

      DOUBLE PRECISION function Fsf4(x,z)
      
c            ->Fsf4(m1^2,m2^2)

      IMPLICIT NONE
      DOUBLE PRECISION x,z,aux

       IF(dabs(x-z).ge.1.d-10)THEN
      aux=(2d0*x*z*dlog(x/z)+z**2-x**2)/(x-z)**3/z
       ELSE
      aux=-1d0/3.d0/x**2
       ENDIF

      Fsf4=aux

      RETURN
      END

************************************************************************************************

      DOUBLE PRECISION function Fsf5(x,y,z)
      
c            ->Fsf5(m1^2,m2^2,m3^2)

      IMPLICIT NONE
      DOUBLE PRECISION x,y,z,aux,Fsf3

      IF(dabs(x-z).ge.1.d-10.and.dabs(y-z).ge.1.d-10)THEN
       IF(dabs(x-y).ge.1.d-10)THEN
        aux=(x*y*dlog(x/y)-x*z*dlog(x/z)+y*z*dlog(y/z))
     .      /((x-y)*(x-z)*(y-z))
       ELSE
        aux=Fsf3(z,y)
       ENDIF
      ELSEIF(dabs(x-z).le.1.d-10.and.dabs(y-z).ge.1.d-10)THEN
       aux=Fsf3(y,z)
      ELSE
       aux=Fsf3(x,z)
      ENDIF

      Fsf5=aux

      RETURN
      END

************************************************************************************************

      DOUBLE PRECISION function Fsf6(x,y,z)
      
c            ->Fsf6(m1^2,m2^2,m3^2)

      IMPLICIT NONE
      DOUBLE PRECISION x,y,z,aux,Fsf7

      IF(dabs(x-z).ge.1.d-10.and.dabs(y-z).ge.1.d-10)THEN
       IF(dabs(x-y).ge.1.d-10)THEN
        aux=x*dlog(x/z)/(x-y)/(x-z)**2
     .     +y*dlog(y/z)/(y-x)/(y-z)**2
     .     +1d0/(x-z)/(y-z)
       ELSE
        aux=Fsf7(x,z)
       ENDIF
      ELSEIF(dabs(x-z).le.1.d-10.and.dabs(y-x).ge.1.d-10)THEN
       aux=(2d0*x*y*dlog(x/y)+y**2-x**2)/(2d0*x*(x-y)**3)
      ELSEIF(dabs(y-z).le.1.d-10.and.dabs(x-y).ge.1.d-10)THEN
       aux=(2d0*x*y*dlog(x/y)+y**2-x**2)/(2d0*y*(x-y)**3)
      ELSE
       aux=-1/6.d0/x**2
      ENDIF

      Fsf6=aux

      RETURN
      END

************************************************************************************************

      DOUBLE PRECISION function Fsf7(x,z)
      
c            ->Fsf7(m1^2,m2^2)

      IMPLICIT NONE
      DOUBLE PRECISION x,z,aux

      IF(dabs(x-z).ge.1.d-10)THEN
       aux=(x*(dlog(z/x)+2d0)-z*(dlog(x/z)+2d0))/(x-z)**3
      ELSE
       aux=-1d0/6.d0/x**2
      ENDIF

      Fsf7=aux

      RETURN
      END
