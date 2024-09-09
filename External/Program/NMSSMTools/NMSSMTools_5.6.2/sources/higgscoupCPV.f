      SUBROUTINE HIGGSCOUP_CPV(PAR)

c         Trilinear (effective) Higgs couplings
c      - The tree-level + loop sfermion/chargino/neutralino/higgs
c        contributions to the trilinear neutral Higgs-couplings gH0H0H0(i,j,k)
c        (i,j,k=1..5) and the charged Higgs-couplings gRH0HPHM(i,j,k)
c        gIH0HPHM(i,j,k), (i=1..5, j,k=1..2, R/I: real/imaginary part) are
c        computed from the parameters of the effective Higgs potential,
c        stored in the common EFFPOTPAR.
c      - The loop-corrections from SM-fermions and gauge-bosons are added
c        before the couplings are stored in the common HICOUP.
c        Here, the Yukawa and gauge couplings are taken at / run to the scale 
c        of the Higgs-mass of 1st index ('i').

      IMPLICIT NONE

      INTEGER I,M,N

      DOUBLE PRECISION PAR(*),ALPHAS,RUNMB,gg1,gg2,aux
      DOUBLE PRECISION XHG(5,6),Pi,Yt,Yb,Ytau,MW2,MZ2,quf1,quf2,Fsf1
      DOUBLE PRECISION PIS111,PIS122,PIS133,PIS144,PIS155,PIS166,
     . PIS211,PIS222,PIS233,PIS244,PIS255,PIS266,PIS311,PIS322,
     . PIS333,PIS344,PIS355,PIS366,PIS312,PIS345,PIS156,PIS256,
     . PIS346,PIS356,PIS513,PIS423,PIS612,PIS456,PIS433,PIS613,
     . PIS466,PIS533,PIS623,PIS566,PIS666,PIS633,PIS246,PIS245,
     . PIS145,PIS455,PIS422,PIS512,PIS544,PIS511,PIS412,PIS411,
     . PIS444,PIS522,PIS555,PIS611,PIS644,PIS622,PIS655
      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION lu,ld,l3,l4,Rel5,Iml5,Rel6,Iml6,Rel7,Iml7,RAud,
     . RAS,K2,lPu,lPd,RlPM,IlPM
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      DOUBLE PRECISION Ytq,Ybq,MTOPQ,MBOTQ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION l,k,Alcos1,Akcos2,muq,nuq
      DOUBLE PRECISION mur,M1r,M2r,msi
      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION gH0H0H0(5,5,5),gRH0HPHM(5,2,2),gIH0HPHM(5,2,2)
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION phiF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      DOUBLE PRECISION Rm3,Im3,RAudt,IAudt,RxS,IxS,Rmsp,Imsp,Ast,IAst
      DOUBLE PRECISION RlPMt,IlPMt,RlM,IlM,RAqs,IAqs,Rltqs,Iltqs
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/STSBSCALE/QSTSB
      COMMON/EFFPOTPAR/lu,ld,l3,l4,Rel5,Iml5,Rel6,Iml6,Rel7,Iml7,RAud,
     . RAS,K2,lPu,lPd,RlPM,IlPM
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/SMFERM/mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      COMMON/QQUARK/Ytq,Ybq,MTOPQ,MBOTQ
      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QPAR/l,k,Alcos1,Akcos2,muq,NUQ
      COMMON/GAUGINOPAR/mur,M1r,M2r,msi
      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/HICOUP/gH0H0H0,gRH0HPHM,gIH0HPHM
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/QEXT/XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/Z3VAUX/phiF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      COMMON/Z3VPOT/Rm3,Im3,RAudt,IAudt,RxS,IxS,Rmsp,Imsp,Ast,IAst,
     . RlPMt,IlPMt,RlM,IlM,RAqs,IAqs,Rltqs,Iltqs


      PI=4d0*DATAN(1d0)
      Ytau=mtau/vd

      DO I=1,5
       XHG(I,1)=XH(I,1)/dsqrt(ZHU)
       XHG(I,2)=XH(I,2)/dsqrt(ZHD)
       XHG(I,3)=XH(I,3)/dsqrt(ZS)
       XHG(I,4)=XH(I,4)*cosb/dsqrt(ZHU)
       XHG(I,5)=XH(I,4)*sinb/dsqrt(ZHD)
       XHG(I,6)=XH(I,5)/dsqrt(ZS)
      ENDDO

c           A: Trilinear couplings of the neutral Higgs

      DO I=1,5

      IF(MH0(I).ge.mt**2)THEN
      Yt=mt/vu
     .  /(1d0+4d0*ALPHAS(MT,2)/(3d0*PI)+11d0*(ALPHAS(MT,2)/PI)**2)
     .      *(1d0+7d0/(4d0*PI)*ALPHAS(MT,2)
     .                        *DLOG(MH0(I)/mt**2))**(-4d0/7d0)
      ELSE
      Yt=mt/vu
     .  /(1d0+4d0*ALPHAS(MT,2)/(3d0*PI)+11d0*(ALPHAS(MT,2)/PI)**2)
     .      *(1d0+23d0/(12d0*PI)*ALPHAS(MT,2)
     .                        *DLOG(MH0(I)/mt**2))**(-12d0/23d0)
      ENDIF

      Yb=RUNMB(dsqrt(MH0(I)))/vd

      MW2=g2/2d0*(vu**2+vd**2)
      MZ2=(g1+g2)/2d0*(vu**2+vd**2)
      gg1=g1/(1d0-g1/16.d0/PI**2*(
     .     DLOG(MH0(I)/MIN(MH0(I),MZ2))*53d0/9d0
     .    +DLOG(MH0(I)/MIN(MH0(I),MT**2))*17d0/18d0
     .    +DLOG(MH0(I)/MIN(MH0(I),MAX(MA2,MZ2)))/6d0
     .    +DLOG(MH0(I)/MIN(MH0(I),MAX(MUQ**2,MZ2)))*2d0/3d0
     .    +DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(7),MZ2)))/18d0
     .    +DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(8),MZ2)))*4d0/9d0
     .    +DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(9),MZ2)))/9d0
     .    +DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(10),MZ2)))/6d0
     .    +DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(11),MZ2)))/3d0
     .    +DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(15),MZ2)))/9d0
     .    +DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(16),MZ2)))*8d0/9d0
     .    +DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(17),MZ2)))*2d0/9d0
     .    +DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(18),MZ2)))/3d0
     .    +DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(19),MZ2)))*2d0/3d0))
      gg2=g2/(1d0+g2/16.d0/PI**2*(
     .     DLOG(MH0(I)/MIN(MH0(I),MZ2))*19d0/6d0
     .    -DLOG(MH0(I)/MIN(MH0(I),MAX(MA2,MZ2)))/6d0
     .    -DLOG(MH0(I)/MIN(MH0(I),MAX(MUQ**2,MZ2)))*2d0/3d0
     .    -DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(7),MZ2)))/2d0
     .    -DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(10),MZ2)))/6d0
     .    -DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(15),MZ2)))
     .    -DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(18),MZ2)))/3d0
     .    -DLOG(MH0(I)/MIN(MH0(I),MAX(PAR(21)**2,MZ2)))*4d0/3d0))
      MW2=gg2/2d0*(vu**2+vd**2)
      MZ2=(gg1+gg2)/2d0*(vu**2+vd**2)

      DO M=1,5
      DO N=1,5

      PIS111=6.d0*XHG(I,1)*XHG(M,1)*XHG(N,1)
      PIS122=2d0*(XHG(I,1)*XHG(M,2)*XHG(N,2)
     . +XHG(I,2)*XHG(M,1)*XHG(N,2)+XHG(I,2)*XHG(M,2)*XHG(N,1))
      PIS133=2d0*(XHG(I,1)*XHG(M,3)*XHG(N,3)
     . +XHG(I,3)*XHG(M,1)*XHG(N,3)+XHG(I,3)*XHG(M,3)*XHG(N,1))
      PIS144=2d0*(XHG(I,1)*XHG(M,4)*XHG(N,4)
     . +XHG(I,4)*XHG(M,1)*XHG(N,4)+XHG(I,4)*XHG(M,4)*XHG(N,1))
      PIS155=2d0*(XHG(I,1)*XHG(M,5)*XHG(N,5)
     . +XHG(I,5)*XHG(M,1)*XHG(N,5)+XHG(I,5)*XHG(M,5)*XHG(N,1))
      PIS166=2d0*(XHG(I,1)*XHG(M,6)*XHG(N,6)
     . +XHG(I,6)*XHG(M,1)*XHG(N,6)+XHG(I,6)*XHG(M,6)*XHG(N,1))
      PIS211=2d0*(XHG(I,2)*XHG(M,1)*XHG(N,1)
     . +XHG(I,1)*XHG(M,2)*XHG(N,1)+XHG(I,1)*XHG(M,1)*XHG(N,2))
      PIS222=6.d0*XHG(I,2)*XHG(M,2)*XHG(N,2)
      PIS233=2d0*(XHG(I,2)*XHG(M,3)*XHG(N,3)
     . +XHG(I,3)*XHG(M,2)*XHG(N,3)+XHG(I,3)*XHG(M,3)*XHG(N,2))
      PIS244=2d0*(XHG(I,2)*XHG(M,4)*XHG(N,4)
     . +XHG(I,4)*XHG(M,2)*XHG(N,4)+XHG(I,4)*XHG(M,4)*XHG(N,2))
      PIS255=2d0*(XHG(I,2)*XHG(M,5)*XHG(N,5)
     . +XHG(I,5)*XHG(M,2)*XHG(N,5)+XHG(I,5)*XHG(M,5)*XHG(N,2))
      PIS266=2d0*(XHG(I,2)*XHG(M,6)*XHG(N,6)
     . +XHG(I,6)*XHG(M,2)*XHG(N,6)+XHG(I,6)*XHG(M,6)*XHG(N,2))
      PIS311=2d0*(XHG(I,3)*XHG(M,1)*XHG(N,1)
     . +XHG(I,1)*XHG(M,3)*XHG(N,1)+XHG(I,1)*XHG(M,1)*XHG(N,3))
      PIS322=2d0*(XHG(I,3)*XHG(M,2)*XHG(N,2)
     . +XHG(I,2)*XHG(M,3)*XHG(N,2)+XHG(I,2)*XHG(M,2)*XHG(N,3))
      PIS333=6.d0*XHG(I,3)*XHG(M,3)*XHG(N,3)
      PIS344=2d0*(XHG(I,3)*XHG(M,4)*XHG(N,4)
     . +XHG(I,4)*XHG(M,3)*XHG(N,4)+XHG(I,4)*XHG(M,4)*XHG(N,3))
      PIS355=2d0*(XHG(I,3)*XHG(M,5)*XHG(N,5)
     . +XHG(I,5)*XHG(M,3)*XHG(N,5)+XHG(I,5)*XHG(M,5)*XHG(N,3))
      PIS366=2d0*(XHG(I,3)*XHG(M,6)*XHG(N,6)
     . +XHG(I,6)*XHG(M,3)*XHG(N,6)+XHG(I,6)*XHG(M,6)*XHG(N,3))
      PIS312=(XHG(I,3)*XHG(M,1)*XHG(N,2)
     . +XHG(I,3)*XHG(M,2)*XHG(N,1)+XHG(I,1)*XHG(M,3)*XHG(N,2)
     . +XHG(I,1)*XHG(M,2)*XHG(N,3)+XHG(I,2)*XHG(M,1)*XHG(N,3)
     . +XHG(I,2)*XHG(M,3)*XHG(N,1))
      PIS345=(XHG(I,3)*XHG(M,4)*XHG(N,5)
     . +XHG(I,3)*XHG(M,5)*XHG(N,4)+XHG(I,4)*XHG(M,3)*XHG(N,5)
     . +XHG(I,4)*XHG(M,5)*XHG(N,3)+XHG(I,5)*XHG(M,4)*XHG(N,3)
     . +XHG(I,5)*XHG(M,3)*XHG(N,4))
      PIS156=(XHG(I,1)*XHG(M,5)*XHG(N,6)
     . +XHG(I,1)*XHG(M,6)*XHG(N,5)+XHG(I,5)*XHG(M,1)*XHG(N,6)
     . +XHG(I,5)*XHG(M,6)*XHG(N,1)+XHG(I,6)*XHG(M,5)*XHG(N,1)
     . +XHG(I,6)*XHG(M,1)*XHG(N,5))
      PIS256=(XHG(I,2)*XHG(M,5)*XHG(N,6)
     . +XHG(I,2)*XHG(M,6)*XHG(N,5)+XHG(I,5)*XHG(M,2)*XHG(N,6)
     . +XHG(I,5)*XHG(M,6)*XHG(N,2)+XHG(I,6)*XHG(M,5)*XHG(N,2)
     . +XHG(I,6)*XHG(M,2)*XHG(N,5))
      PIS346=(XHG(I,3)*XHG(M,4)*XHG(N,6)
     . +XHG(I,3)*XHG(M,6)*XHG(N,4)+XHG(I,4)*XHG(M,3)*XHG(N,6)
     . +XHG(I,4)*XHG(M,6)*XHG(N,3)+XHG(I,6)*XHG(M,4)*XHG(N,3)
     . +XHG(I,6)*XHG(M,3)*XHG(N,4))
      PIS356=(XHG(I,3)*XHG(M,5)*XHG(N,6)
     . +XHG(I,3)*XHG(M,6)*XHG(N,5)+XHG(I,5)*XHG(M,3)*XHG(N,6)
     . +XHG(I,5)*XHG(M,6)*XHG(N,3)+XHG(I,6)*XHG(M,5)*XHG(N,3)
     . +XHG(I,6)*XHG(M,3)*XHG(N,5))
      PIS513=(XHG(I,5)*XHG(M,1)*XHG(N,3)
     . +XHG(I,5)*XHG(M,3)*XHG(N,1)+XHG(I,1)*XHG(M,5)*XHG(N,3)
     . +XHG(I,1)*XHG(M,3)*XHG(N,5)+XHG(I,3)*XHG(M,1)*XHG(N,5)
     . +XHG(I,3)*XHG(M,5)*XHG(N,1))
      PIS423=(XHG(I,4)*XHG(M,2)*XHG(N,3)
     . +XHG(I,4)*XHG(M,3)*XHG(N,2)+XHG(I,2)*XHG(M,4)*XHG(N,3)
     . +XHG(I,2)*XHG(M,3)*XHG(N,4)+XHG(I,3)*XHG(M,2)*XHG(N,4)
     . +XHG(I,3)*XHG(M,4)*XHG(N,2))
      PIS612=(XHG(I,6)*XHG(M,1)*XHG(N,2)
     . +XHG(I,6)*XHG(M,2)*XHG(N,1)+XHG(I,1)*XHG(M,6)*XHG(N,2)
     . +XHG(I,1)*XHG(M,2)*XHG(N,6)+XHG(I,6)*XHG(M,1)*XHG(N,2)
     . +XHG(I,6)*XHG(M,2)*XHG(N,1))
      PIS456=(XHG(I,4)*XHG(M,5)*XHG(N,6)
     . +XHG(I,4)*XHG(M,6)*XHG(N,5)+XHG(I,5)*XHG(M,4)*XHG(N,6)
     . +XHG(I,5)*XHG(M,6)*XHG(N,4)+XHG(I,6)*XHG(M,5)*XHG(N,4)
     . +XHG(I,6)*XHG(M,4)*XHG(N,5))
      PIS433=2d0*(XHG(I,4)*XHG(M,3)*XHG(N,3)
     . +XHG(I,3)*XHG(M,4)*XHG(N,3)+XHG(I,3)*XHG(M,3)*XHG(N,4))
      PIS613=(XHG(I,6)*XHG(M,1)*XHG(N,3)
     . +XHG(I,6)*XHG(M,3)*XHG(N,1)+XHG(I,1)*XHG(M,6)*XHG(N,3)
     . +XHG(I,1)*XHG(M,3)*XHG(N,6)+XHG(I,3)*XHG(M,1)*XHG(N,6)
     . +XHG(I,3)*XHG(M,6)*XHG(N,1))/dsqrt(ZHU)/ZS
      PIS466=2d0*(XHG(I,4)*XHG(M,6)*XHG(N,6)
     . +XHG(I,6)*XHG(M,4)*XHG(N,6)+XHG(I,6)*XHG(M,6)*XHG(N,4))
      PIS533=2d0*(XHG(I,5)*XHG(M,3)*XHG(N,3)
     . +XHG(I,3)*XHG(M,5)*XHG(N,3)+XHG(I,3)*XHG(M,3)*XHG(N,5))
      PIS623=(XHG(I,6)*XHG(M,2)*XHG(N,3)
     . +XHG(I,6)*XHG(M,3)*XHG(N,2)+XHG(I,2)*XHG(M,6)*XHG(N,3)
     . +XHG(I,2)*XHG(M,3)*XHG(N,6)+XHG(I,3)*XHG(M,2)*XHG(N,6)
     . +XHG(I,3)*XHG(M,6)*XHG(N,2))/dsqrt(ZHD)/ZS
      PIS566=2d0*(XHG(I,5)*XHG(M,6)*XHG(N,6)
     . +XHG(I,6)*XHG(M,5)*XHG(N,6)+XHG(I,6)*XHG(M,6)*XHG(N,5))
      PIS666=6.d0*XHG(I,6)*XHG(M,6)*XHG(N,6)
      PIS633=2d0*(XHG(I,6)*XHG(M,3)*XHG(N,3)
     . +XHG(I,3)*XHG(M,6)*XHG(N,3)+XHG(I,3)*XHG(M,3)*XHG(N,6))
      PIS246=(XHG(I,2)*XHG(M,4)*XHG(N,6)
     . +XHG(I,2)*XHG(M,6)*XHG(N,4)+XHG(I,4)*XHG(M,2)*XHG(N,6)
     . +XHG(I,4)*XHG(M,6)*XHG(N,2)+XHG(I,6)*XHG(M,4)*XHG(N,2)
     . +XHG(I,6)*XHG(M,2)*XHG(N,4))
      PIS245=(XHG(I,2)*XHG(M,4)*XHG(N,5)
     . +XHG(I,2)*XHG(M,5)*XHG(N,4)+XHG(I,4)*XHG(M,2)*XHG(N,5)
     . +XHG(I,4)*XHG(M,5)*XHG(N,2)+XHG(I,5)*XHG(M,4)*XHG(N,2)
     . +XHG(I,5)*XHG(M,2)*XHG(N,4))
      PIS145=(XHG(I,1)*XHG(M,4)*XHG(N,5)
     . +XHG(I,1)*XHG(M,5)*XHG(N,4)+XHG(I,4)*XHG(M,1)*XHG(N,5)
     . +XHG(I,4)*XHG(M,5)*XHG(N,1)+XHG(I,5)*XHG(M,4)*XHG(N,1)
     . +XHG(I,5)*XHG(M,1)*XHG(N,4))
      PIS422=2d0*(XHG(I,4)*XHG(M,2)*XHG(N,2)
     . +XHG(I,2)*XHG(M,4)*XHG(N,2)+XHG(I,2)*XHG(M,2)*XHG(N,4))
      PIS455=2d0*(XHG(I,4)*XHG(M,5)*XHG(N,5)
     . +XHG(I,5)*XHG(M,4)*XHG(N,5)+XHG(I,5)*XHG(M,5)*XHG(N,4))
      PIS512=(XHG(I,2)*XHG(M,1)*XHG(N,5)
     . +XHG(I,2)*XHG(M,5)*XHG(N,1)+XHG(I,1)*XHG(M,2)*XHG(N,5)
     . +XHG(I,1)*XHG(M,5)*XHG(N,2)+XHG(I,5)*XHG(M,1)*XHG(N,2)
     . +XHG(I,5)*XHG(M,2)*XHG(N,1))
      PIS511=2d0*(XHG(I,5)*XHG(M,1)*XHG(N,1)
     . +XHG(I,1)*XHG(M,5)*XHG(N,1)+XHG(I,1)*XHG(M,1)*XHG(N,5))
      PIS544=2d0*(XHG(I,5)*XHG(M,4)*XHG(N,4)
     . +XHG(I,4)*XHG(M,5)*XHG(N,4)+XHG(I,4)*XHG(M,4)*XHG(N,5))
      PIS412=(XHG(I,2)*XHG(M,1)*XHG(N,4)
     . +XHG(I,2)*XHG(M,4)*XHG(N,1)+XHG(I,1)*XHG(M,2)*XHG(N,4)
     . +XHG(I,1)*XHG(M,4)*XHG(N,2)+XHG(I,4)*XHG(M,1)*XHG(N,2)
     . +XHG(I,4)*XHG(M,2)*XHG(N,1))
      PIS444=6.d0*XHG(I,4)*XHG(M,4)*XHG(N,4)
      PIS411=2d0*(XHG(I,4)*XHG(M,1)*XHG(N,1)
     . +XHG(I,1)*XHG(M,4)*XHG(N,1)+XHG(I,1)*XHG(M,1)*XHG(N,4))
      PIS555=6.d0*XHG(I,5)*XHG(M,5)*XHG(N,5)
      PIS522=2d0*(XHG(I,5)*XHG(M,2)*XHG(N,2)
     . +XHG(I,2)*XHG(M,5)*XHG(N,2)+XHG(I,2)*XHG(M,2)*XHG(N,5))
      PIS611=2d0*(XHG(I,6)*XHG(M,1)*XHG(N,1)
     . +XHG(I,1)*XHG(M,6)*XHG(N,1)+XHG(I,1)*XHG(M,1)*XHG(N,6))
      PIS644=2d0*(XHG(I,6)*XHG(M,4)*XHG(N,4)
     . +XHG(I,4)*XHG(M,6)*XHG(N,4)+XHG(I,4)*XHG(M,4)*XHG(N,6))
      PIS622=2d0*(XHG(I,6)*XHG(M,2)*XHG(N,2)
     . +XHG(I,2)*XHG(M,6)*XHG(N,2)+XHG(I,2)*XHG(M,2)*XHG(N,6))
      PIS655=2d0*(XHG(I,6)*XHG(M,5)*XHG(N,5)
     . +XHG(I,5)*XHG(M,6)*XHG(N,5)+XHG(I,5)*XHG(M,5)*XHG(N,6))

c      I- Couplings from effective potential parameters (Tree-level+Sfermions+Inos)

        IF(k.ne.0d0)THEN
      gH0H0H0(I,M,N)=(lu*vu*(PIS111+PIS144)+ld*vd*(PIS222+PIS255)
     . +(l3+l4)*(vu*(PIS122+PIS155)+vd*(PIS211+PIS244))
     . +Rel5*(vu*(PIS122-PIS155-2d0*PIS245)
     .                 +vd*(PIS211-PIS244-2d0*PIS145))
     . +Iml5*(vu*(PIS455-PIS422-2d0*PIS512)
     .                 +vd*(PIS544-PIS511-2d0*PIS412))
     . -Rel6*(vu*(3.d0*PIS211+PIS244-2d0*PIS145)+vd*(PIS111+PIS144))
     . +Iml6*(vu*(3.d0*PIS511+2d0*PIS412+PIS544)+vd*(PIS444+PIS411))
     . -Rel7*(vu*(PIS222+PIS255)+vd*(3.d0*PIS122+PIS155-2d0*PIS245))
     . +Iml7*(vu*(PIS522+PIS555)+vd*(3.d0*PIS422+2d0*PIS512+PIS455))
     . -RAud*(PIS312-PIS345-PIS156-PIS246)
     . -RlPM*(2d0*muq/l*(PIS312-PIS345+PIS156+PIS246)
     .             +vd*(PIS133-PIS166+2d0*PIS346)
     .             +vu*(PIS233-PIS266+2d0*PIS356))
     . +IlPM*(muq/l*(PIS423+PIS513-3.d0*PIS612+3.d0*PIS456)
     .            +vd*(PIS433-2d0*PIS613-PIS466)
     .            +vu*(PIS533-2d0*PIS623-PIS566)
     .       -vu*vd*l/muq*(PIS666-3.d0*PIS633))
     . +lPu*(muq/l*(PIS311+PIS344)+vu*(PIS133+PIS166))
     . +lPd*(muq/l*(PIS322+PIS355)+vd*(PIS233+PIS266))
     . +RAS/3.d0*(PIS333-3.d0*PIS366)
     . +2d0*K2*muq/l*(PIS333+PIS366)
     .                                                 )/dsqrt(2d0)

      aux=Max(dabs(XIS),dabs(MSP),dabs(XIF),dabs(MUP),dabs(M3H))
      IF(aux.ge.1d-4)THEN
      gH0H0H0(I,M,N)=gH0H0H0(I,M,N)+1d0/dsqrt(2d0)*(
     . -RAudt*(PIS312-PIS345+PIS156+PIS246)
     . +IAudt*(PIS423+PIS513-PIS612+PIS456)
     . -RlPMt*(2d0*muq/l*(PIS312-PIS345-PIS156-PIS246)
     .             +vd*(PIS133-PIS166-2d0*PIS346)
     .             +vu*(PIS233-PIS266-2d0*PIS356))
     . +IlPMt*(2d0*muq/l*(PIS423+PIS513+PIS612-PIS456)
     .            +vd*(PIS433+2d0*PIS613-PIS466)
     .            +vu*(PIS533+2d0*PIS623-PIS566))
     . -RlM*(2d0*muq/l*(PIS312-PIS345)
     .             +vd*(PIS133+PIS166)+vu*(PIS233+PIS266))
     . +IlM*(2d0*muq/l*(PIS423+PIS513)
     .             +vd*(PIS433+PIS466)+vu*(PIS533+PIS566))
     . -(Im3*l/muq+IAudt+(IlPMt+IlM)*muq/l)
     .                      *(PIS423+PIS513+PIS612-PIS456)
     . +RAqs*(PIS311+PIS344+PIS322+PIS355)
     . +2d0*Rltqs*(muq/l*(PIS311+PIS344+PIS322+PIS355)
     .             +vuq*(PIS133-PIS166)+vdq*(PIS233-PIS266))
     . -2d0*Iltqs*(muq/l*(PIS611+PIS644+PIS622+PIS655)
     .                   +2d0*vuq*PIS613+2d0*vdq*PIS623)
     . -((IAudt-2d0*IlPMt*muq/l)*vuq*vdq-(IAqs+2d0*Iltqs*muq/l)
     .     *(vuq**2+vdq**2)+IxS+Imsp*muq/l)*(l/muq)**2
     .                               *(PIS666/3d0-PIS633)
     . +Ast*(PIS333+PIS366)-4d0/3d0*IASt*PIS666)
      ENDIF
        ELSE
      gH0H0H0(I,M,N)=(lu*vu*(PIS111+PIS144)+ld*vd*(PIS222+PIS255)
     . +(l3+l4)*(vu*(PIS122+PIS155)+vd*(PIS211+PIS244))
     . +Rel5*(vu*(PIS122-PIS155-2d0*PIS245)
     .                 +vd*(PIS211-PIS244-2d0*PIS145))
     . +Iml5*(vu*(PIS455-PIS422-2d0*PIS512)
     .                 +vd*(PIS544-PIS511-2d0*PIS412))
     . -Rel6*(vu*(3.d0*PIS211+PIS244-2d0*PIS145)+vd*(PIS111+PIS144))
     . +Iml6*(vu*(3.d0*PIS511+2d0*PIS412+PIS544)+vd*(PIS444+PIS411))
     . -Rel7*(vu*(PIS222+PIS255)+vd*(3.d0*PIS122+PIS155-2d0*PIS245))
     . +Iml7*(vu*(PIS522+PIS555)+vd*(3.d0*PIS422+2d0*PIS512+PIS455))
     . -RAud*(PIS312-PIS345-PIS156-PIS246)
     . -RlPM*(2d0*muq/l*(PIS312-PIS345+PIS156+PIS246)
     .             +vd*(PIS133-PIS166+2d0*PIS346)
     .             +vu*(PIS233-PIS266+2d0*PIS356))
     . +IlPM*(muq/l*(PIS423+PIS513-3.d0*PIS612+3.d0*PIS456)
     .            +vd*(PIS433-2d0*PIS613-PIS466)
     .            +vu*(PIS533-2d0*PIS623-PIS566))
     . +lPu*(muq/l*(PIS311+PIS344)+vu*(PIS133+PIS166))
     . +lPd*(muq/l*(PIS322+PIS355)+vd*(PIS233+PIS266))
     . -RAudt*(PIS312-PIS345+PIS156+PIS246)
     . +IAudt*(PIS423+PIS513-PIS612+PIS456)
     . -RlPMt*(2d0*muq/l*(PIS312-PIS345-PIS156-PIS246)
     .             +vd*(PIS133-PIS166-2d0*PIS346)
     .             +vu*(PIS233-PIS266-2d0*PIS356))
     . +IlPMt*(2d0*muq/l*(PIS423+PIS513+PIS612-PIS456)
     .            +vd*(PIS433+2d0*PIS613-PIS466)
     .            +vu*(PIS533+2d0*PIS623-PIS566))
     . -RlM*(2d0*muq/l*(PIS312-PIS345)
     .             +vd*(PIS133+PIS166)+vu*(PIS233+PIS266))
     . +IlM*(2d0*muq/l*(PIS423+PIS513)
     .             +vd*(PIS433+PIS466)+vu*(PIS533+PIS566))
     . -(Im3*l/muq+IAudt+(IlPMt+IlM)*muq/l)
     .                      *(PIS423+PIS513+PIS612-PIS456)
     . +RAqs*(PIS311+PIS344+PIS322+PIS355)
     . +2d0*Rltqs*(muq/l*(PIS311+PIS344+PIS322+PIS355)
     .             +vuq*(PIS133-PIS166)+vdq*(PIS233-PIS266))
     . -2d0*Iltqs*(muq/l*(PIS611+PIS644+PIS622+PIS655)
     .                   +2d0*vuq*PIS613+2d0*vdq*PIS623)
     . +Ast*(PIS333+PIS366)-4d0/3d0*IASt*PIS666
     .                                                 )/dsqrt(2d0)
        ENDIF

c      II- Singlet potential from higgsinos

      gH0H0H0(I,M,N)=gH0H0H0(I,M,N)-1d0/16.d0/Pi**2/dsqrt(2d0)
     . *(l/muq)**3*(4.d0*mur**4*(dlog(mur**2/QSTSB)+2d0/3.d0)*PIS333
     . +4.d0*mur**4*dlog(mur**2/QSTSB)*PIS366
     . +2d0*ks2si**3*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .      *(dlog(Max(1d0,msi**2)/QSTSB)
     .      +2d0/3.d0*(mupsi*DDCOS(Phi02-phiP)+ks2si)**2/msi**2)*PIS333
     . +2d0*ks2si**3*(mupsi*DDCOS(Phi02-phiP)+ks2si)
     .      *(dlog(Max(1d0,msi**2)/QSTSB)
     .      +2d0*(mupsi*DDSIN(Phi02-phiP))**2/msi**2)*PIS366
     . -2d0*ks2si**3*mupsi*DDSIN(Phi02-phiP)
     .         *(dlog(Max(1d0,msi**2)/QSTSB)
     .      +2d0*(mupsi*DDCOS(Phi02-phiP)+ks2si)**2/msi**2)*PIS633
     . -2d0*ks2si**3*mupsi*DDSIN(Phi02-phiP)
     .         *(dlog(Max(1d0,msi**2)/QSTSB)
     .      +2d0/3d0*(mupsi*DDSIN(Phi02-phiP))**2/msi**2)*PIS666)

c      III- Doublet potential from SM fermions

      gH0H0H0(I,M,N)=gH0H0H0(I,M,N)-3.d0/8.d0/Pi**2/dsqrt(2d0)*(
     .  Yt**4*vu*(dlog((Yt*vu)**2/QSTSB)+2d0/3.d0)*PIS111
     . +Yt**4*vu*dlog((Yt*vu)**2/QSTSB)*PIS144
     . +Yb**4*vd*(dlog((Yb*vd)**2/QSTSB)+2d0/3.d0)*PIS222
     . +Yb**4*vd*dlog((Yb*vd)**2/QSTSB)*PIS255)
     . -Ytau**4*vd/8.d0/Pi**2/dsqrt(2d0)*(
     .                (dlog(mtau**2/QSTSB)+2d0/3.d0)*PIS222
     .                         +dlog(mtau**2/QSTSB)*PIS255)

c      IV- Doublet potential from SM gauge bosons

      gH0H0H0(I,M,N)=gH0H0H0(I,M,N)+3.d0/64.d0/Pi**2/dsqrt(2d0)*(
     .  (gg2**2*dlog(MW2/QSTSB)+(gg1+gg2)**2/2d0*dlog(MZ2/QSTSB))
     .             *(vu*(PIS111+PIS122+PIS144+PIS155)
     .               +vd*(PIS222+PIS211+PIS255+PIS244))
     .  +(2d0*gg2**2+(gg1+gg2)**2)/3.d0/(vu**2+vd**2)
     .    *(vu**3*PIS111+3.d0*vu**2*vd*PIS211+3.d0*vu*vd**2*PIS122
     .                                                  +vd**3*PIS222))

      ENDDO
      ENDDO
      ENDDO

c           B: Couplings of 1 neutral Higgs to 2 charged ones

      DO I=1,5
      DO M=1,2
      DO N=1,2

c      I- Couplings from effective potential parameters (Tree-level+Sfermions+Inos)

      gRH0HPHM(I,M,N)=(lu*vu*XHG(I,1)*XC(M,1)*XC(N,1)
     . +ld*vd*XHG(I,2)*XC(M,2)*XC(N,2)
     . +l3*(vu*XHG(I,1)*XC(M,2)*XC(N,2)+vd*XHG(I,2)*XC(M,1)*XC(N,1))
     . -(l4+Rel5)*(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1))/2d0
     .                          *(vu*XHG(I,2)+vd*XHG(I,1))
     . +Iml5*(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1))/2d0
     .                          *(vu*XHG(I,5)+vd*XHG(I,4))
     . +Rel6*(vu*(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1))*XHG(I,1)
     .           -XC(M,1)*XC(N,1)*(vu*XHG(I,2)+vd*XHG(I,1)))
     . +Iml6*XC(M,1)*XC(N,1)*(vu*XHG(I,5)+vd*XHG(I,4))
     . +Rel7*(vd*(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1))*XHG(I,2)
     .           -XC(M,2)*XC(N,2)*(vu*XHG(I,2)+vd*XHG(I,1)))
     . +Iml7*XC(M,2)*XC(N,2)*(vu*XHG(I,5)+vd*XHG(I,4))
     . +lPu*muq/l*XC(M,1)*XC(N,1)*XHG(I,3)
     . +lPd*muq/l*XC(M,2)*XC(N,2)*XHG(I,3)
     . +RAud*(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1))/2d0*XHG(I,3)
     . +RlPM*muq/l*(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1))*XHG(I,3)
     . +3.d0/2d0*IlPM*muq/l*(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1))*XHG(I,6)
     .                                        )*dsqrt(2d0)

      gIH0HPHM(I,M,N)=((l4-Rel5)/2d0*(vd*XHG(I,4)+vu*XHG(I,5))
     . -Iml5/2d0*(vd*XHG(I,1)+vu*XHG(I,2))
     . +Iml6*vu*XHG(I,1)+Iml7*vd*XHG(I,2)
     . +(RAud/2d0-RlPM*muq/l)*XHG(I,6)+IlPM*muq/l/2d0*XHG(I,3)
     .                 )*(XC(M,1)*XC(N,2)-XC(M,2)*XC(N,1))*dsqrt(2d0)


c      II- Doublet potential from SM fermions

      gRH0HPHM(I,M,N)=gRH0HPHM(I,M,N)-3.d0/8.d0/Pi**2*dsqrt(2d0)*(
     .  2d0*(Yt**2*vu*XHG(I,1)*quf1((Yt*vu)**2,(Yb*vd)**2,QSTSB)
     .       +Yb**2*vd*XHG(I,2)*quf1((Yb*vd)**2,(Yt*vu)**2,QSTSB))
     .            *(Yt**2*XC(M,1)*XC(N,1)+Yb**2*XC(M,2)*XC(N,2))
     . +Yt**2*Yb**2*(vd*XHG(I,1)*quf2((Yt*vu)**2,(Yb*vd)**2,QSTSB)
     .               +vu*XHG(I,2)*quf2((Yb*vd)**2,(Yt*vu)**2,QSTSB))
     .            *(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1)))
     .       -1d0/8.d0/Pi**2/dsqrt(2d0)*quf1(mtau**2,0d0,QSTSB)
     .        *2d0*Ytau**4*vd*XHG(I,2)*XC(M,2)*XC(N,2)

      gIH0HPHM(I,M,N)=gIH0HPHM(I,M,N)-3.d0/8.d0/Pi**2*dsqrt(2d0)
     .      *(-Yt**2*Yb**2)*Fsf1((Yt*vu)**2,(Yb*vd)**2,QSTSB)
     .       *(vd*XHG(I,1)+vu*XHG(I,2))
     .       *(XC(M,1)*XC(N,2)-XC(M,2)*XC(N,1))

      IF(aux.ge.1d-4)THEN

      gRH0HPHM(I,M,N)=gRH0HPHM(I,M,N)+1d0/dsqrt(2d0)*(
     . ((RAudt+2d0*(RlPMt+RlM)*muq/l)*XHG(I,3)
     .      +(Im3*l/muq+2d0*IAudt+(IlM-IlPMt)*muq/l)*XHG(I,6))
     .         *(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1))/2d0
     . +(RAqs+2d0*Rltqs*muq/l)*XHG(I,3)
     .                   *(XC(M,1)*XC(N,1)+XC(M,2)*XC(N,2)))

      gIH0HPHM(I,M,N)=gIH0HPHM(I,M,N)+1d0/dsqrt(2d0)/2d0*(
     . (-(RAudt-2d0*RlPMt*muq/l)*XHG(I,6)
     .  +(Im3*l/muq+(IlPMt+IlM)*muq/l)*XHG(I,3))
     .                *(XC(M,1)*XC(N,2)-XC(M,2)*XC(N,1)))

      ENDIF

c      III- Doublet potential from SM gauge bosons

      gRH0HPHM(I,M,N)=gRH0HPHM(I,M,N)+3.d0/64.d0*dsqrt(2d0)/Pi**2*(
     . (gg2**2*dlog(MW2/QSTSB)+(gg1+gg2)**2/2d0*dlog(MZ2/QSTSB))
     .   *(vu*XHG(I,1)+vd*XHG(I,2))*(XC(M,1)*XC(N,1)+XC(M,2)*XC(N,2))
     . -gg1*gg2*(dlog(MZ2/QSTSB)-1d0)*(
     .          2d0*vu*XHG(I,1)*XC(M,2)*XC(N,2)
     .         +2d0*vd*XHG(I,2)*XC(M,1)*XC(N,1)
     .         +(vd*XHG(I,1)+vu*XHG(I,2))
     .                   *(XC(M,1)*XC(N,2)+XC(M,2)*XC(N,1)))
     . -2d0*gg1*gg2*(vu*XHG(I,1)+vd*XHG(I,2))/(vu**2+vd**2)
     .         *(vd**2*XC(M,1)*XC(N,1)+vu**2*XC(M,2)*XC(N,2)
     .                   +vu*vd*(XC(M,2)*XC(N,1)+XC(M,1)*XC(N,2))))

      gIH0HPHM(I,M,N)=gIH0HPHM(I,M,N)+3.d0/64.d0*dsqrt(2d0)/Pi**2*(
     . gg1*gg2*(dlog(MZ2/QSTSB)-1d0)*(vd*XHG(I,4)+vu*XHG(I,5))
     .                   *(XC(M,1)*XC(N,2)-XC(M,2)*XC(N,1)))


      ENDDO
      ENDDO
      ENDDO


      RETURN
      END


************************************************************************************************

      DOUBLE PRECISION function quf1(x,y,z)
      
c            ->quf1(m1^2,m2^2,Q^2)

      IMPLICIT NONE
      DOUBLE PRECISION x,y,z,aux
      IF(min(x,y).ge.1d-10)THEN
       IF(dabs(x-y).ge.1d-10)THEN
      aux=(x**2*dlog(x/z)-x*y*(2d0*dlog(x/z)-1d0)
     . +y**2*(dlog(y/z)-1d0))/(x-y)**2
       ELSE
      aux=dlog(x/z)+1d0/2d0
       ENDIF
      ELSEIF(min(x,y).le.1d-10)THEN
      IF(x.ge.y)aux=dlog(x/z)
      IF(y.ge.x)aux=dlog(y/z)-1d0
      ENDIF

      quf1=aux

      RETURN
      END


************************************************************************************************

      DOUBLE PRECISION function quf2(x,y,z)
      
c            ->quf2(m1^2,m2^2,Q^2)

      IMPLICIT NONE
      DOUBLE PRECISION x,y,z,aux
      IF(min(x,y).ge.1d-10)THEN
       IF(dabs(x-y).ge.1d-10)THEN
      aux=(x**2*(dlog(x/z)+1d0)-x*y*(3.d0*dlog(x/z)-dlog(y/z))
     . +y**2*(dlog(y/z)-1d0))/(x-y)**2
       ELSE
      aux=dlog(x/z)+1d0
       ENDIF
      ELSEIF(min(x,y).le.1d-10)THEN
      IF(x.ge.y)aux=dlog(x/z)+1d0
      IF(y.ge.x)aux=dlog(y/z)-1d0
      ENDIF

      quf2=aux

      RETURN
      END
