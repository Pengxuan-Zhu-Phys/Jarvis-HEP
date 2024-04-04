      SUBROUTINE runpar_CPV(PAR)

      IMPLICIT NONE

      INTEGER OMGFLAG,MAFLAG,MOFLAG,Q2FIX

      DOUBLE PRECISION PAR(*),MYPHASES(16)
      DOUBLE PRECISION Yt,Yb,Ytau,alsmt,runmb
      DOUBLE PRECISION mu,M1,M2,M3,mhs2,mas2
      DOUBLE PRECISION QSUSY,QSL,QSQ,l,k,IAL,IAk,IXIS
      DOUBLE PRECISION Pi,NMB0,NMB1,MA2,MP2
      DOUBLE PRECISION ALSMA,DLA,DLQA,F1,F2,HTMA
      DOUBLE PRECISION Q2MIN,Q2,QSTSB
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      DOUBLE PRECISION lq,kq,Alcos1,Akcos2,muq,nuq
      DOUBLE PRECISION Ytq,Ybq,MTOPQ,MBOTQ
      DOUBLE PRECISION PhiAT,PhiAB,aux1,aux2
      DOUBLE PRECISION MQ3,MU3,MD3,ML3,ME3,MQ,MU1,MD,ML,ME
      DOUBLE PRECISION XI,HT2,HB2,L2,K2,HTAU2,MH1,MH2
      DOUBLE PRECISION SP,SIG1,SIG2,SIG3,F20,F21,F22
      DOUBLE PRECISION MHuS,MHdS,MSS,muH2,AT,AB
      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION ALSMZC,ALEMMZC,GFC,g1C,g2C,S2TWC
      DOUBLE PRECISION MSC,MCC,MBC,MBPC,MTC,MTAUC,MMUC,MZC,MWC
      DOUBLE PRECISION MPIC,MELC,MSTRANGEC
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiATQ,phiABQ,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION mur,M1r,M2r,msi
      DOUBLE PRECISION MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      DOUBLE PRECISION MSQ1,MSU1,MSD1
      DOUBLE PRECISION Drv
      DOUBLE PRECISION phiF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      DOUBLE PRECISION MSQ3,MSU3,MSD3,ATP,ABP
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION RXIS,RXIF,IAlQ2,IAkQ2,IXISQ2,MHuQ2,MHdQ2,MSQ2
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/HIGGSMS/muH2
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/QEXT/XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      COMMON/QPAR/lq,kq,Alcos1,Akcos2,muq,NUQ
      COMMON/QQUARK/Ytq,Ybq,MTOPQ,MBOTQ
      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/SMFERM/mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      COMMON/GAUGE/ALSMZC,ALEMMZC,GFC,g1C,g2C,S2TWC
      COMMON/SMSPEC/MSC,MCC,MBC,MBPC,MTC,MTAUC,MMUC,MZC,MWC
      COMMON/SMEXT/MPIC,MELC,MSTRANGEC
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiATQ,phiABQ,phiATAU,phiAC,phiAS,phiAMU
      COMMON/GAUGINOPAR/mur,M1r,M2r,msi
      COMMON/SLEPPAR/MSL3,MSE3,MSL1,MSE1,ATAU,AMU
      COMMON/SQUPAR/MSQ1,MSU1,MSD1
      COMMON/VEVSHIFT/Drv
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/Z3VAUX/phiF,phiP,phi3,phiS,phiSP,phi3q,phiSq,phiSPq,
     .              mupsi,ks2si
      COMMON/MYPHASES/MYPHASES
      COMMON/MH2TREE/MHuS,MHdS,MSS
      COMMON/IMALAK/IAL,IAK,IXIS
      COMMON/RADCOR2/MSQ3,MSU3,MSD3,ATP,ABP
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/REXI/RXIS,RXIF,IAlQ2,IAkQ2,IXISQ2,MHuQ2,MHdQ2,MSQ2

      PI=4d0*DATAN(1d0)

c      1) Gauge (at MZ)
      GF=GFC
      MZ=MZC
      MW=MWC
      g1=g1C
      g2=g2C
      ALSMZ=ALSMZC
      S2TW=S2TWC
      ALEMMZ=ALEMMZC

c      2) tanb
      tanb=PAR(3)
      cosb=1d0/dsqrt(1d0+tanb**2)
      sinb=tanb*cosb
      Drv=0d0
      vu=sinb/dsqrt(2d0*dsqrt(2d0)*GF*(1d0-Drv))
      vd=cosb/dsqrt(2d0*dsqrt(2d0)*GF*(1d0-Drv))

c      3) Phases
      phi01=MYPHASES(1) !phil+phi1
      phi02=MYPHASES(2) !phik+phi2
      phi0=phi01-phi02
      phiM1=MYPHASES(3)
      phiM2=MYPHASES(4)
      phiM3=MYPHASES(5)
      phiAT=MYPHASES(6)
      phiAB=MYPHASES(7)
      phiATAU=MYPHASES(8)
      phiAC=MYPHASES(9)
      phiAS=MYPHASES(10)
      phiAMU=MYPHASES(11)
      IF(PAR(2).ne.0d0)THEN
       phiS=MYPHASES(12)
      ELSE
       phiS=0d0
      ENDIF
      phiSP=MYPHASES(13)
      phiF=MYPHASES(14)
      phiP=MYPHASES(15)
      phi3=MYPHASES(16)
      RXIF=XIF*DDCOS(phiF)
      RXIS=XIS*DDCOS(phiS)
      IXIS=XIS*DDSIN(phiS)

c      4) Yukawa (at mt) and SM fermions
      mt=MTC
      mb=MBC            ! MSbar, scale MB
      mtau=MTAUC
      mmu=MMUC
      mel=MELC

      MS=MSC
      MC=MCC
      MBP=MBPC

      MSTRANGE=MSTRANGEC
      MPI=MPIC
      alsmt=ALSMZ/(1d0+23.d0/12d0/Pi*ALSMZ*dlog((mt/MZ)**2))
      Yt=mt/vu
      Yt=Yt/(1d0+4.d0/3.d0/Pi*alsmt+11d0/Pi**2*alsmt**2)
c      Yb=mb/vd*(1d0-23.d0/12d0/Pi*alsmt
c     c   *dlog(mt**2/MZ**2))**(12d0/23d0)
      Yb=RUNMB(mt)/vd
      Ytau=mtau/vd

c      4) Sfermion parameters
      MSL3=PAR(10)
      MSE3=PAR(11)
      MSL1=PAR(18)
      MSE1=PAR(19)
      ATAU=PAR(14)
      AMU=PAR(25)
      MSQ1=PAR(15)
      MSU1=PAR(16)
      MSD1=PAR(17)

*   Definition of the SUSY scale Q2, unless defined initially
      IF(Q2FIX.EQ.0)THEN
       Q2=MAX((2d0*PAR(15)+PAR(16)+PAR(17))/4d0,Q2MIN)
      ENDIF

*   Definition of the scale QSTSB
      QSTSB=DSQRT(MAX(PAR(7)*PAR(8),Q2MIN**2))
      QSUSY=dsqrt(QSTSB)
      QSL=dsqrt((2d0*MSL3+MSE3+2d0*MSE1+4.d0*MSL1)/9.d0)
      QSQ=dsqrt((2d0*PAR(7)+PAR(8)+PAR(9)
     .       +4d0*PAR(15)+2d0*PAR(16)+2d0*PAR(17))/12d0)

c      5) Higgs parameters
       l=PAR(1)
       k=PAR(2)
       mu=PAR(4)
       M1=PAR(20)
       IF(M1.ge.0d0)THEN
        M1=Max(1d0,M1)
       ELSE
        M1=-Max(1d0,dabs(M1))
       ENDIF
       M2=PAR(21)
       IF(M2.ge.0d0)THEN
        M2=Max(1d0,M2)
       ELSE
        M2=-Max(1d0,dabs(M2))
       ENDIF
       M3=PAR(22)
       IF(M3.ge.0d0)THEN
        M3=Max(1d0,M3)
       ELSE
        M3=-Max(1d0,dabs(M3))
       ENDIF
       msi=Max(1d0,
     . dsqrt(4d0*(k/l*mu)**2+muP**2+4d0*k/l*muP*mu*DDCOS(phi02-phiP)))

      IF(MAFLAG.LT.0)THEN

       Alcos1=PAR(5)
       Akcos2=PAR(6)
       MA2=PAR(23)**2

      ELSE

       IF(MOD(MAFLAG,3).EQ.0)THEN
        Alcos1=PAR(5)
        MA2=((Alcos1+k*mu/l*DDCOS(Phi0))*mu+M3H*DDCOS(phi3)
     .           +mu*mup*DDCOS(phi01-phip)+l*xif*DDCOS(phi01-phiF))
     .          *(tanb+1d0/tanb)
        PAR(23)=DSQRT(MAX(MA2,1d0))
       ELSEIF(MOD(MAFLAG,3).EQ.1)THEN
        MA2=PAR(23)**2
        Alcos1=(MA2*TANB/(1d0+TANB**2)-M3H*DDCOS(phi3)
     .          -l*xif*DDCOS(phi01-phiF))/mu-MUP*DDCOS(phi01-phip)
     .          -k*mu/l*DDCOS(Phi0)
        PAR(5)=Alcos1
       ELSE
        Alcos1=PAR(5)
        MA2=PAR(23)**2
        aux1=XIF
        RXIF=(MA2*TANB/(1d0+TANB**2)-(Alcos1+k*mu/l*DDCOS(Phi0))*mu
     .  -M3H*DDCOS(phi3)-mu*mup*DDCOS(phi01-phip))/(l*DDCOS(phi01))
     .  -datan(phi01)*XIF
        IF(RXIF.ne.0d0)THEN
         phiF=datan(aux1/RXIF)
        ELSEIF(aux1.gt.0d0)THEN
         phiF=Pi/2d0
        ELSEIF(aux1.lt.0d0)THEN
         phiF=-Pi/2d0
        ELSE
         phiF=0d0
        ENDIF
        XIF=DSQRT(RXIF**2+aux1**2)
        IF(RXIF.lt.0d0)XIF=-XIF
       ENDIF

       IF(MAFLAG/3.EQ.0)THEN 
        Akcos2=PAR(6)
        MP2=-k/l*mu*(3.d0*Akcos2+mup*DDCOS(phi02-phip))
     . +l**2*vu*vd/mu*(Alcos1+muP*DDCOS(phiP-phi01)
     .                              +4.d0*k/l*mu*DDCOS(phi01-phi02))
     . -l/mu*(xiS*DDCOS(phiS)+XIF*mup*DDCOS(phiP-phiF))
     . -2d0*MSP*DDCOS(phiSP)-4d0*k*XIF*DDCOS(phi02-phiF)
        PAR(24)=DSQRT(MAX(MP2,1d0))
       ELSEIF(MAFLAG/3.EQ.1)THEN
        MP2=PAR(24)**2
        IF(K.EQ.0d0)THEN
         Akcos2=0d0
         PAR(6)=0d0
        ELSE
         Akcos2=(l**2*(Alcos1+4d0*k*mu/l*DDCOS(Phi0)
     .                        +mup*DDCOS(phi01-phip))*vu*vd/mu
     .  -l*(XIF*mup*DDCOS(phip-phiF)+xiS*DDCOS(phiS))/mu
     .  -mup*k*mu/l*DDCOS(phi02-phip)-4d0*k*XIF*DDCOS(phi02-phiF)
     .  -2d0*MSP*DDCOS(phiSP)-MP2)/(3d0*k*mu/l*DDCOS(Phi0))
         PAR(6)=Akcos2
        ENDIF
       ELSE
        MP2=PAR(24)**2
        Akcos2=PAR(6)
        aux1=XIS
        RXIS=(l*(Alcos1+4d0*k*mu/l*DDCOS(Phi0)+mup*DDCOS(phi01-phip))
     .      *vu*vd-(k*mu/l*(3d0*Akcos2+mup*DDCOS(phi02-phip))
     .       +4d0*k*XIF*DDCOS(phi02-phiF)+2d0*MSP*DDCOS(phiSP)+MP2)*mu/l
     .      -XIF*muP*DDCOS(phip-phiF))
        IF(k.ne.0d0)THEN
         IF(RXIS.ne.0d0)THEN
          phiS=datan(aux1/RXIS)
         ELSEIF(aux1.gt.0d0)THEN
          phiS=Pi/2d0
         ELSEIF(aux1.lt.0d0)THEN
          phiS=-Pi/2d0
         ELSE
          phiS=0d0
         ENDIF
         XIS=DSQRT(RXIS**2+aux1**2)
         IF(RXIS.lt.0d0)XIS=-XIS
        ELSE
         XIS=RXIS
        ENDIF
       ENDIF

      ENDIF

      MHuS=(MU*(Alcos1+MUP*DDCOS(phi01-phiP)+k/l*mu*DDCOS(Phi0))
     .         +M3H*DDCOS(Phi3)+L*XIF*DDCOS(phi01-phiF))/TANB
     .      -L**2*VD**2-MU**2+(g1+g2)/4d0*(VD**2-VU**2)
      MHdS=(MU*(Alcos1+MUP*DDCOS(phi01-phiP)+k/l*mu*DDCOS(Phi0))
     .         +M3H*DDCOS(Phi3)+L*XIF*DDCOS(phi01-phiF))*TANB
     .      -L**2*VU**2-MU**2+(g1+g2)/4d0*(VU**2-VD**2)
      MSS=L**2*vu*vd/mu
     .        *(Alcos1+2d0*k/l*mu*DDCOS(phi0)+mup*DDCOS(Phi01-PhiP))
     . -k/l*mu*(Akcos2+2d0*k/l*mu)-l**2*(vu**2+vd**2)
     .  -XIF*(2d0*K*DDCOS(phi02-PhiF)+L*MUP*DDCOS(PhiP-PhiF)/MU)
     .  -MUP**2-3d0*MUP*K/L*MU*DDCOS(Phi02-PhiP)-MSP*DDCOS(PhiSP)
     .  -L*XIS/MU*DDCOS(PhiS)

      IAl=-k/l*mu*DDSIN(phi0)-(M3H*DDSIN(phi3)+l*XIF*DDSIN(phi01-phiF)
     .                              +mu*mup*DDSIN(phi01-phip))/mu
      IF(k.ne.0d0)THEN
         IAk=l**2/mu**2/k*(l*vu*vd*(IAl-muP*DDSIN(phi01-phiP)
     .                                     -2d0*k/l*mu*DDSIN(phi0))
     .      -xiS*DDSIN(phiS)-mu/l*MSP*DDSIN(phiSP)
     .      -muP*k*(mu/l)**2*DDSIN(phi02-phiP)
     .      -XIF*(muP*DDSIN(phiP-phiF)+2d0*k/l*mu*DDSIN(phi02-phiF)))
      ELSE
         IAk=0d0
         aux1=XIS
         IXIS=l*vu*vd*(IAl-muP*DDSIN(phi01-phiP))-mu/l*MSP*DDSIN(phiSP)
     .          -XIF*muP*DDSIN(phiP-phiF)
         IF(aux1.ne.0d0)THEN
          phiS=datan(IXIS/aux1)
         ELSEIF(IXIS.gt.0d0)THEN
          phiS=Pi/2d0
         ELSEIF(IXIS.lt.0d0)THEN
          phiS=-Pi/2d0
         ELSE
          phiS=0d0
         ENDIF
         XIS=dsqrt(aux1**2+IXIS**2)
         IF(aux1.lt.0d0)XIS=-XIS
      ENDIF

c      6) Running gauge couplings

      ALSQ=ALSMT/(1d0+ALSMT/(4d0*PI)*(7d0*DLOG(MAX(QSTSB,MT**2)/MT**2)
     .         -2d0*DLOG(MAX(QSTSB,PAR(22)**2)/MAX(PAR(22)**2,MT**2))))

*   g_2**2 including the Higgs and sparticle thresholds:
      g2q=g2/(1d0+g2/16d0/Pi**2*(DLOG(QSTSB/MZ**2)*19d0/6d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(MA2,MZ**2)))/6d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(MU**2,MZ**2)))*2d0/3d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(7),MZ**2)))/2d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(10),MZ**2)))/6d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(15),MZ**2)))
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(18),MZ**2)))/3d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(M2**2,MZ**2)))*4d0/3d0))

*   g_1**2 including the top, Higgs and sparticle thresholds:
      g1q=g1/(1d0-g1/16d0/Pi**2*(DLOG(QSTSB/MZ**2)*53d0/9d0
     .    +DLOG(QSTSB/MT**2)*17d0/18d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(MA2,MZ**2)))/6d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(MU**2,MZ**2)))*2d0/3d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(7),MZ**2)))/18d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(8),MZ**2)))*4d0/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(9),MZ**2)))/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(10),MZ**2)))/6d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(11),MZ**2)))/3d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(15),MZ**2)))/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(16),MZ**2)))*8d0/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(17),MZ**2)))*2d0/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(18),MZ**2)))/3d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(19),MZ**2)))*2d0/3d0))

      gq=(g1q+g2q)/2d0

c      7) Running Yukawa couplings

      ALSMA=ALSMT/(1d0+ALSMT/(4d0*PI)
     . *(7d0*DLOG(MIN(MAX(MA2,MT**2),QSTSB)/MT**2)
     .   -2d0*DLOG(MAX(MA2,PAR(22)**2)/MAX(PAR(22)**2,MT**2))))
      DLA=(ALSMA/ALSMT)**(1d0/7d0)
      DLQA=(ALSQ/ALSMA)**(1d0/7d0)
      F1=DABS(1d0-9d0*sinb**2*Yt**2*(1d0-DLA)/(8d0*PI*ALSMT))
      HTMA=YT*DLA**4/DSQRT(F1)
      F2=DABS(1d0-9d0*HTMA**2*(1d0-DLQA)/(8d0*PI*ALSMA))

*   Top/Bottom Yukawas at QSTSB
*   (Note: RUNMB(Q) includes QCD corrections only)
*   including electroweak contributions, Logs ht^2*LQT resummed:
*   + Conversion to DR_bar:
      Ytq=Yt/DSQRT(F1*F2)*(1d0+7d0/(4d0*PI)*ALSMT
     .                    *DLOG(MAX(QSTSB,MT**2)/MT**2))**(-4d0/7d0)
     .   *(1d0+((-17d0/6d0*g1q-9d0/2d0*g2q+Yb**2)
     .                    *DLOG(MAX(QSTSB,MT**2)/MT**2)
     .   +((3d0*cosb**2-1d0)*Yb**2+2d0*Ytau**2*cosb**2)
     .                    *DLOG(MIN(MAX(MA2,MT**2),QSTSB)/MT**2)
     .   -2d0*l**2*DLOG(MIN(MAX(mu**2,msi**2,MZ**2),QSTSB)/QSTSB)
     .   -G1Q*DLOG(MIN(MAX(M1**2,mu**2,MZ**2),QSTSB)/QSTSB)
     .   -3d0*G2Q*DLOG(MIN(MAX(M2**2,mu**2,MZ**2),QSTSB)/QSTSB))
     .                                     /64d0/Pi**2)
     .   *(1d0-ALSQ/(3d0*PI)+g2q*3d0/128d0/Pi**2)

      Ybq=RUNMB(DSQRT(QSTSB))/vd*F1**(-1d0/6d0)
     .   *(1d0-3d0*HTMA**2*(1d0-DLQA)/(8d0*PI*ALSMA))**(-1d0/6d0)
     .   *(1d0+((-5d0/6d0*g1q-9d0/2d0*g2q+9d0*Yb**2+2d0*Ytau**2)
     .                     *DLOG(MAX(QSTSB,MT**2)/MT**2)
     .   +(-9d0*Yb**2-2d0*Ytau**2)*sinb**2
     .                     *DLOG(MIN(MAX(MA2,MT**2),QSTSB)/MT**2)
     .   -2d0*L**2*DLOG(MIN(MAX(mu**2,msi**2,MZ**2),QSTSB)/QSTSB)
     .   -G1Q*DLOG(MIN(MAX(M1**2,mu**2,MZ**2),QSTSB)/QSTSB)
     .   -3d0*G2Q*DLOG(MIN(MAX(M2**2,mu**2,MZ**2),QSTSB)/QSTSB))
     .                                     /64d0/Pi**2)
     .   *(1d0-ALSQ/(3d0*PI)+g2q*3d0/128d0/Pi**2)

c      8) Higgs parameters at the SUSY scale

      LQ=L*(1d0+(-G1Q-3d0*G2Q+4d0*L**2+2d0*K**2
     .           +3d0*Ytq**2+3d0*Ybq**2+Ytau**2)/32d0/Pi**2
     .                                      *DLOG(QSTSB/Q2))
      KQ=K*(1d0+3d0/16d0/Pi**2*(L**2+K**2)*DLOG(QSTSB/Q2))
      MUQ=MU*(1d0+(3d0*Ytq**2+3d0*Ybq**2+Ytau**2+2d0*l**2
     .             -G1Q-3d0*G2Q)/32d0/Pi**2*DLOG(QSTSB/Q2))
      NUQ=MUQ*KQ/LQ
      MUPQ=MUP*(1d0+(L**2+K**2)/8d0/Pi**2*DLOG(QSTSB/Q2))
      XIFQ=XIF*(1d0+(L**2+K**2)/16d0/Pi**2*DLOG(QSTSB/Q2))
      MSPQ=MSP+(2d0*L**2*(MSP+2d0*MUP
     .      *(Alcos1*DDCOS(phip-phi01-phiSP)
     .      -IAl*DDSIN(phip-phi01-phiSP)))
     . +4d0*K**2*(MSP+MUP
     .      *(Akcos2*DDCOS(phip-phi02-phiSP)
     .      -IAk*DDSIN(phip-phi02-phiSP)))
     . +4d0*L*K*M3H*DDCOS(phi3-phiSP-phi0))/16d0/Pi**2*DLOG(QSTSB/Q2)
      M3HQ=M3H+((3d0*Ytq**2+3d0*Ybq**2+Ytau**2+6d0*l**2-G1Q-3d0*G2Q)
     .    *M3H+2d0*L*K*MSP*DDCOS(phiSP-phi3+phi0))/32d0/Pi**2
     .                                *DLOG(QSTSB/Q2)
      IF(XIS.ne.0d0)THEN
       phiSq=phiS+(2d0*XIF*L**2
     .  *(Alcos1*DDSIN(phiF-phi01-phiS)+IAl*DDCOS(phiF-phi01-phiS))
     . +2d0*XIF*K**2
     .  *(Akcos2*DDSIN(phiF-phi02-phiS)+IAk*DDCOS(phiF-phi02-phiS))
     . +2d0*L*M3H*(IAl*DDCOS(phi3+phiS)-Alcos1*DDSIN(phi3+phiS)
     .                        +MUP*DDSIN(phi3+phiP-phi01-phiS))
     . +K*MSP*(IAk*DDCOS(phiSP-phiS)-Akcos2*DDSIN(phiSP-phiS)
     .         +MUP*DDSIN(phip+phiSP-phi02-phiS))
     . +2d0*k*mup*MSS*DDSIN(phi02-phiP-phiS))
     .       /16d0/Pi**2/XIS*DLOG(QSTSB/Q2)
      ELSE
       phiSq=0d0
      ENDIF
      IF(MSP.ne.0d0)THEN
       phiSPq=phiSP+(4d0*L**2*(IAl*DDCOS(phip-phi01-phiSP)
     .                         +Alcos1*DDSIN(phip-phi01-phiSP))
     . +4d0*K**2*MUP*(IAk*DDCOS(phip-phi02-phiSP)
     .                     +Akcos2*DDSIN(phip-phi02-phiSP))
     . +4d0*L*K*M3H*DDSIN(phi3-phiSP-phi0))
     .    /MSP/16d0/Pi**2*DLOG(QSTSB/Q2)
      ELSE
       phiSPq=0d0
      ENDIF
      IF(M3H.ne.0d0)THEN
       phi3q=phi3+l*k*MSP/M3H*DDSIN(phiSP-phi3+phi0)
     .                        /16d0/Pi**2*DLOG(QSTSB/Q2)
      ELSE
       phi3q=0d0
      ENDIF

      IF(MAFLAG.LT.0)THEN

       Alcos1=Alcos1+(4d0*l**2*Alcos1
     .    +2d0*k**2*(Akcos2*DDCOS(phi0)-IAk*DDSIN(phi0))
     .    +3d0*Ytq**2*PAR(12)*DDCOS(phiAt+phi01)
     .    +3d0*Ybq**2*PAR(13)*DDCOS(phiAb+phi01)
     .    +Ytau**2*PAR(14)*DDCOS(phiAtau+phi01)
     .    +g1q*M1*DDCOS(phiM1+phi01)
     .    +3d0*g2q*M2*DDCOS(phiM2+phi01))*DLOG(QSTSB/Q2)/16d0/Pi**2
       Akcos2=Akcos2+6d0*(l**2*(PAR(5)*DDCOS(phi0)+IAl*DDSIN(phi0))
     .     +k**2*Akcos2)*DLOG(QSTSB/Q2)/16d0/Pi**2
       XIFQ=XIF*(1d0+(L**2+K**2)/16d0/Pi**2*DLOG(QSTSB/Q2))
       XISQ=XIS+(L**2*(XIS+2d0*XIF
     .  *(Alcos1*DDCOS(phiF-phi01-phiS)-IAl*DDSIN(phiF-phi01-phiS)))
     . +K**2*(XIS+2d0*XIF
     .  *(Akcos2*DDCOS(phiF-phi02-phiS)-IAk*DDSIN(phiF-phi02-phiS)))
     . +2d0*L*M3H*(Alcos1*DDCOS(phi3+phiS)+IAl*DDSIN(phi3+phiS)
     .                        +MUP*DDCOS(phi3+phiP-phi01-phiS))
     . +K*MSP*(Akcos2*DDCOS(phiSP+phiS)+IAk*DDSIN(phiSP+phiS)
     .                +MUP*DDCOS(phip+phiSP-phi02-phiS))
     . +2d0*k*mup*MSS*DDCOS(phi02-phiP-phiS))
     .                      /16d0/Pi**2*DLOG(QSTSB/Q2)

      ELSE

       IF(MOD(MAFLAG,3).EQ.0)THEN
        Alcos1=Alcos1+(4d0*l**2*Alcos1
     .    +2d0*k**2*(Akcos2*DDCOS(phi0)-IAk*DDSIN(phi0))
     .    +3d0*Ytq**2*PAR(12)*DDCOS(phiAt+phi01)
     .    +3d0*Ybq**2*PAR(13)*DDCOS(phiAb+phi01)
     .    +Ytau**2*PAR(14)*DDCOS(phiAtau+phi01)
     .    +g1q*M1*DDCOS(phiM1+phi01)
     .    +3d0*g2q*M2*DDCOS(phiM2+phi01))*DLOG(QSTSB/Q2)/16d0/Pi**2
        XIFQ=XIF*(1d0+(L**2+K**2)/16d0/Pi**2*DLOG(QSTSB/Q2))
        MA2=((Alcos1+NUQ*DDCOS(Phi0)+MUPQ*DDCOS(PhiP-Phi01))*MUQ
     . +M3HQ*DDCOS(phi3)+LQ*XIFQ*DDCOS(PhIF-Phi01))*(TANB+1d0/TANB)
        PAR(23)=DSQRT(MAX(MA2,1d0))
       ELSEIF(MOD(MAFLAG,3).EQ.1)THEN
        XIFQ=XIF*(1d0+(L**2+K**2)/16d0/Pi**2*DLOG(QSTSB/Q2))
        Alcos1=(MA2*TANB/(1d0+TANB**2)-M3HQ*DDCOS(phi3)
     . -LQ*XIFQ*DDCOS(PhIF-Phi01))/MUQ-MUPQ*DDCOS(PhiP-Phi01)
     . -NUQ*DDCOS(Phi0)
        PAR(5)=Alcos1-(4d0*l**2*Alcos1
     .    +2d0*k**2*(Akcos2*DDCOS(phi0)-IAk*DDSIN(phi0))
     .    +3d0*Ytq**2*PAR(12)*DDCOS(phiAt+phi01)
     .    +3d0*Ybq**2*PAR(13)*DDCOS(phiAb+phi01)
     .    +Ytau**2*PAR(14)*DDCOS(phiAtau+phi01)
     .    +g1q*M1*DDCOS(phiM1+phi01)
     .    +3d0*g2q*M2*DDCOS(phiM2+phi01))*DLOG(QSTSB/Q2)/16d0/Pi**2
       ELSE
        Alcos1=Alcos1+(4d0*l**2*Alcos1
     .    +2d0*k**2*(Akcos2*DDCOS(phi0)-IAk*DDSIN(phi0))
     .    +3d0*Ytq**2*PAR(12)*DDCOS(phiAt+phi01)
     .    +3d0*Ybq**2*PAR(13)*DDCOS(phiAb+phi01)
     .    +Ytau**2*PAR(14)*DDCOS(phiAtau+phi01)
     .    +g1q*M1*DDCOS(phiM1+phi01)
     .    +3d0*g2q*M2*DDCOS(phiM2)+phi01)*DLOG(QSTSB/Q2)/16d0/Pi**2
        aux1=((MA2*TANB/(1d0+TANB**2)
     .   -(Alcos1+NUQ*DDCOS(Phi0))*MUQ-M3HQ*DDCOS(Phi3)
     .   -MUQ*MUPQ*DDCOS(PhiP-Phi01))/LQ)
        aux2=XIF*DDSIN(PhiF)*(1d0+(L**2+K**2)/16d0/Pi**2*DLOG(QSTSB/Q2)) ! Im(XIFQ)
        XIFQ=(aux1/DDCOS(Phi01)-aux2*dtan(Phi01))                        ! Re(XIFQ)
        IF(XIFQ.ne.0d0)THEN
         PhiF=datan(aux2/XIFQ)                                          ! arg(XIFQ)
        ELSEIF(aux2.gt.0d0)THEN
         PhiF=Pi/2d0
        ELSEIF(aux2.lt.0d0)THEN
         PhiF=-Pi/2d0
        ELSE
         PhiF=0d0
        ENDIF
        IF(XIFQ.ge.0d0)THEN
         XIFQ=dsqrt(XIFQ**2+aux2**2)                                    ! |XIFQ|
        ELSE
         XIFQ=-dsqrt(XIFQ**2+aux2**2)
        ENDIF
        XIF=XIFQ*(1d0-(L**2+K**2)/16d0/Pi**2*DLOG(QSTSB/Q2))
        RXIF=XIF*DDCOS(PhiF)
       ENDIF

       IF(MAFLAG/3.EQ.0)THEN
        IF(k.ne.0d0)then
         Akcos2=Akcos2+6d0*(l**2*(PAR(5)*DDCOS(phi0)+IAl*DDSIN(phi0))
     .         +k**2*Akcos2)*DLOG(QSTSB/Q2)/16d0/Pi**2
        ELSE
         Akcos2=0d0
        ENDIF
        XISQ=XIS+(L**2*(XIS+2d0*XIF
     .  *(Alcos1*DDCOS(phiF-phi01-phiS)-IAl*DDSIN(phiF-phi01-phiS)))
     . +K**2*(XIS+2d0*XIF
     .  *(Akcos2*DDCOS(phiF-phi02-phiS)-IAk*DDSIN(phiF-phi02-phiS)))
     . +2d0*L*M3H*(Alcos1*DDCOS(phi3+phiS)+IAl*DDSIN(phi3+phiS)
     .                        +MUP*DDCOS(phi3+phiP-phi01-phiS))
     . +K*MSP*(Akcos2*DDCOS(phiSP+phiS)+IAk*DDSIN(phiSP+phiS)
     .                +MUP*DDCOS(phip+phiSP-phi02-phiS))
     . +2d0*k*mup*MSS*DDCOS(phi02-phiP-phiS))
     .                      /16d0/Pi**2*DLOG(QSTSB/Q2)
        MP2=LQ**2*(Alcos1+4d0*NUQ+MUPQ)*vu*vd/MUQ-3d0*Akcos2*NUQ
     .    -LQ*(XIFQ*MUPQ+XISQ)/MUQ-MUPQ*NUQ-4d0*KQ*XIFQ-2d0*MSPQ
        PAR(24)=DSQRT(MAX(MP2,1d0))
       ELSEIF(MAFLAG/3.EQ.1)THEN
        XISQ=XIS+(L**2*(XIS+2d0*XIF
     .  *(Alcos1*DDCOS(phiF-phi01-phiS)-IAl*DDSIN(phiF-phi01-phiS)))
     . +K**2*(XIS+2d0*XIF
     .  *(Akcos2*DDCOS(phiF-phi02-phiS)-IAk*DDSIN(phiF-phi02-phiS)))
     . +2d0*L*M3H*(Alcos1*DDCOS(phi3+phiS)+IAl*DDSIN(phi3+phiS)
     .                        +MUP*DDCOS(phi3+phiP-phi01-phiS))
     . +K*MSP*(Akcos2*DDCOS(phiSP+phiS)+IAk*DDSIN(phiSP+phiS)
     .                +MUP*DDCOS(phip+phiSP-phi02-phiS))
     . +2d0*k*mup*MSS*DDCOS(phi02-phiP-phiS))
     .                      /16d0/Pi**2*DLOG(QSTSB/Q2)
        IF(K.EQ.0d0)THEN
         Akcos2=0d0
         PAR(6)=0d0
        ELSE
         Akcos2=(LQ**2*(Alcos1+4d0*NUQ+MUPQ)*vu*vd/MUQ
     .      -LQ*(XIFQ*MUPQ+XISQ)/MUQ-MUPQ*NUQ
     .      -4d0*KQ*XIFQ-2d0*MSPQ-MP2)/(3d0*NUQ)
         PAR(6)=Akcos2-6d0*(l**2*(PAR(5)*DDCOS(phi0)+IAl*DDSIN(phi0))
     .     +k**2*Akcos2)*DLOG(QSTSB/Q2)/16d0/Pi**2
        ENDIF
       ELSE
        Akcos2=Akcos2+6d0*(l**2*(PAR(5)*DDCOS(phi0)+IAl*DDSIN(phi0))
     .     +k**2*Akcos2)*DLOG(QSTSB/Q2)/16d0/Pi**2
        XISQ=LQ*(Alcos1+4d0*NUQ*DDCOS(Phi0)
     .       +MUPQ*DDCOS(phi01-phiP))*vu*vd
     .       -MUQ/LQ*((3d0*Akcos2+MUPQ*DDCOS(phi02-phiP))*NUQ
     .       +4d0*KQ*XIFQ*DDCOS(phi02-phiF)+2d0*MSPQ*DDCOS(phiSPq)+MP2)
     .       -XIFQ*MUPQ*DDCOS(phip-phiF)                                  ! Re(XISQ)
        aux2=XIS*DDSIN(PhiS)*(1d0+(L**2+K**2)/16d0/Pi**2*DLOG(QSTSB/Q2))
     .  +(2d0*XIF*(L**2*(PAR(5)*DDSIN(phiF-phi01)+IAl*DDCOS(phiF-phi01))
     .            +K**2*(PAR(6)*DDSIN(phiF-phi02)
     .            +IAk*DDCOS(phiF-phi02)))
     .   +2d0*L*M3H*(IAl*DDCOS(phi3)-PAR(5)*DDSIN(phi3)
     .                             +MUP*DDSIN(phi3+phiP-phi01))
     .   +K*MSP*(IAk*DDCOS(phiSP)-PAR(6)*DDSIN(phiSP)
     .                             +MUP*DDSIN(phip+phiSP-phi02))
     .   +2d0*k*mup*MSS*DDSIN(phi02-phiP))/16d0/Pi**2*DLOG(QSTSB/Q2)      ! Im(XISQ)
        IF(XISQ.ne.0d0)THEN
         phiSq=datan(aux2/XISQ)
        ELSEIF(aux2.gt.0d0)THEN
         phiSq=Pi/2d0
        ELSEIF(aux2.lt.0d0)THEN
         phiSq=-Pi/2d0
        ELSE
         phiSq=0d0
        ENDIF
        IF(XISQ.ge.0d0)THEN
         XISQ=DSQRT(XISQ**2+aux2**2)
        ELSE
         XISQ=-DSQRT(XISQ**2+aux2**2)
        ENDIF
        RXIS=XISQ*DDCOS(PhiSq)
     .        *(1d0-(L**2+K**2)/16d0/Pi**2*DLOG(QSTSB/Q2))
     .  -(2d0*XIF*(L**2*(PAR(5)*DDCOS(phiF-phi01)-IAl*DDSIN(phiF-phi01))
     .            +K**2*(PAR(6)*DDCOS(phiF-phi02)
     .            -IAk*DDSIN(phiF-phi02)))
     .   +2d0*L*M3H*(PAR(5)*DDCOS(phi3)+IAl*DDSIN(phi3)
     .                             +MUP*DDCOS(phi3+phiP-phi01))
     .   +K*MSP*(PAR(6)*DDCOS(phiSP)+IAk*DDSIN(phiSP)
     .                             +MUP*DDCOS(phip+phiSP-phi02))
     .   +2d0*k*mup*MSS*DDCOS(phi02-phiP))/16d0/Pi**2*DLOG(QSTSB/Q2)
        aux2=XIS*DDSIN(PhiS)
        IF(RXIS.ne.0d0)THEN
         phiS=datan(aux2/RXIS)
        ELSEIF(RXIS.eq.0d0.and.aux2.gt.0d0)THEN
         phiS=Pi/2d0
        ELSEIF(RXIS.eq.0d0.and.aux2.lt.0d0)THEN
         phiS=-Pi/2d0
        ELSE
         phiS=0d0
        ENDIF
        XIS=DSQRT(RXIS**2+aux2**2)
        IF(RXIS.lt.0d0)XIS=-XIS
       ENDIF
      ENDIF

c      9) Gaugino parameters

      mhs2=Max(l**2*vu*vd/mu*(PAR(5)+mup*DDCOS(phip-phi01))
     . +k/l*mu*(PAR(6)+3d0*mup*DDCOS(phiP-phi02)+4.d0*k/l*mu)
     . -l/mu*(xiS*DDCOS(phiS)+XIF*mup*DDCOS(phiP-phiF)),1d0)

      mas2=Max(-k/l*mu*(3.d0*PAR(6)+mup*DDCOS(phi02-phip))
     . +l**2*vu*vd/mu*(PAR(5)+muP*DDCOS(phiP-phi01)
     .                              +4.d0*k/l*mu*DDCOS(phi01-phi02))
     . -l/mu*(xiS*DDCOS(phiS)+XIF*mup*DDCOS(phiP-phiF))
     . -2d0*MSP*DDCOS(phiSP)-4d0*k*XIF*DDCOS(phi02-phiF),1d0)

      mur=mu*(1d0-1d0/(64d0*Pi**2)*(
     .       12d0*(Ytq**2+Ybq**2)*NMB1(mu**2,0d0,QSUSY**2,Q2)
     . +3d0*g2*(NMB1(mu**2,M2**2,PAR(23)**2,Q2)
     .    +NMB1(MU**2,M2**2,MZ**2,Q2)+2d0*NMB1(MU**2,MU**2,MZ**2,Q2)
     .    -4d0*NMB0(MU**2,MU**2,MZ**2,Q2)
     .     +2d0*sinb*cosb*mu/M2*(NMB0(MU**2,M2**2,PAR(23)**2,Q2)
     .     -NMB0(MU**2,M2**2,MZ**2,Q2))*DDCOS(PhiM2+Phi01))
     . +g1*(NMB1(mu**2,M1**2,PAR(23)**2,Q2)
     .    +NMB1(MU**2,M1**2,MZ**2,Q2)+2d0*NMB1(MU**2,MU**2,MZ**2,Q2)
     .    -4d0*NMB0(MU**2,MU**2,MZ**2,Q2)
     .   +2d0*sinb*cosb*M1/mu*(NMB0(MU**2,M1**2,PAR(23)**2,Q2)
     .     -NMB0(MU**2,M1**2,MZ**2,Q2))*DDCOS(PhiM1+Phi01))
     . +2d0*l**2*(NMB1(mu**2,msi**2,PAR(23)**2,Q2)
     .    +NMB1(MU**2,msi**2,MZ**2,Q2)
     .    -2d0*sinb*cosb*2d0*k/l*(NMB0(MU**2,msi**2,PAR(23)**2,Q2)
     .     -NMB0(MU**2,msi**2,MZ**2,Q2))*DDCOS(Phi01-Phi02))))

      M1r=M1*(1d0-g1/(16d0*Pi**2)*(11d0*NMB1(M1**2,0d0,QSQ**2,Q2)
     .     +9d0* NMB1(M1**2,0d0,QSL**2,Q2)
     .     +2d0*sinb*cosb*mu/M1*(NMB0(M1**2,MU**2,PAR(23)**2,Q2)
     .     -NMB0(M1**2,MU**2,MZ**2,Q2))*DDCOS(PhiM1+Phi01)
     .     +NMB1(M1**2,MU**2,PAR(23)**2,Q2)+NMB1(M1**2,MU**2,MZ**2,Q2)))

      M2r=M2*(1d0-g2/(16d0*Pi**2)*(9d0*NMB1(M2**2,0d0,QSQ**2,Q2)
     .     +3d0* NMB1(M2**2,0d0,QSL**2,Q2)
     .     +2d0*sinb*cosb*mu/M2*(NMB0(M2**2,MU**2,PAR(23)**2,Q2)
     .     -NMB0(M2**2,MU**2,MZ**2,Q2))*DDCOS(PhiM2+Phi01)
     .     +NMB1(M2**2,MU**2,PAR(23)**2,Q2)+NMB1(M2**2,MU**2,MZ**2,Q2)
     .  -8d0*NMB0(M2**2,M2**2,MW**2,Q2)+4d0*NMB1(M2**2,M2**2,MW**2,Q2)))

      msi=msi*(1d0-2d0*l**2/(16d0*Pi**2)*(
     .   NMB1(msi**2,MU**2,PAR(23)**2,Q2)+NMB1(msi**2,MU**2,MZ**2,Q2)))

      msi=msi-dsqrt(4d0*(k/l*mu)**2+muP**2
     .   +4d0*k/l*muP*mu*DDCOS(phi02-phiP))*2d0*k**2/(16d0*Pi**2)*(
     .   -(NMB0(Max(1d0,msi**2),msi**2,mhs2,Q2)
     .     -NMB0(Max(1d0,msi**2),msi**2,mas2,Q2))
     .   +NMB1(Max(1d0,msi**2),msi**2,mhs2,Q2)
     .   +NMB1(Max(1d0,msi**2),msi**2,mas2,Q2))


      mupsi=mup*msi/Max(dsqrt(4d0*(k/l*mu)**2+muP**2
     .                         +4d0*k/l*muP*mu*DDCOS(phi02-phiP)),1d0)

      ks2si=2d0*k/l*mu*msi/Max(dsqrt(4d0*(k/l*mu)**2+muP**2
     .                         +4d0*k/l*muP*mu*DDCOS(phi02-phiP)),1d0)

      ZHU=1d0+1d0/16.d0/Pi**2*(3.d0*Ytq**2
cUE     c                         *NMB0(muH2,mtopq**2,mtopq**2,QSTSB)
     .                         *NMB0(muH2,mt**2,mt**2,QSTSB)
     .    -(g1q+g2q)/2d0*NMB0(muH2,MZ**2,MZ**2,QSTSB)*sinb**2 !PAR(23)**2
     .    -g2q*NMB0(muH2,MW**2,MW**2,QSTSB)*sinb**2 !PAR(23)**2
     .    +g1q/2d0*NMB0(muH2,mur**2,M1r**2,QSTSB)
     .    +3.d0*g2q/2d0*NMB0(muH2,mur**2,M2r**2,QSTSB)
     .    +l**2*NMB0(muH2,mur**2,msi**2,QSTSB))
      ZHD=1d0+1d0/16.d0/Pi**2*(3.d0*Ybq**2
cUE     c                         *NMB0(muH2,mbotq**2,mbotq**2,QSTSB)
     .                         *NMB0(muH2,mb**2,mb**2,QSTSB)
     .    +Ytau**2*NMB0(muH2,mtau**2,mtau**2,QSTSB)
c     c    -(g1q+g2q)/2d0*NMB0(muH2,MZ**2,MZ**2,QSTSB) !PAR(23)**2
c     c    -g2q*NMB0(muH2,MW**2,MW**2,QSTSB) !PAR(23)**2
     .    +g1q/2d0*NMB0(muH2,mur**2,M1r**2,QSTSB)
     .    +3.d0*g2q/2d0*NMB0(muH2,mur**2,M2r**2,QSTSB)
     .    +l**2*NMB0(muH2,mur**2,msi**2,QSTSB))
      ZS=1d0+1d0/8.d0/Pi**2
     .    *(l**2*NMB0(muH2,mur**2,mur**2,QSTSB)
     .     +k**2*NMB0(muH2,msi**2,msi**2,QSTSB))

      vuq=vu/dsqrt(ZHU)
      vdq=vd/dsqrt(ZHD)
      tanbq=vuq/vdq
      MTOPQ=YTQ*vuQ
      MBOTQ=YBQ*vdQ
      muq=muq/dsqrt(ZS)

      IAl=nuq*DDSIN(phi0)-(M3HQ*DDSIN(phi3q)+lq*XIFQ*DDSIN(phi01-phiF)
     .                              +muq*mupq*DDSIN(phi01-phip))/muq
      IF(kq.ne.0d0)THEN
         IAk=lq**2/muq**2/kq*(lq*vuq*vdq*(IAl-muPQ*DDSIN(phi01-phiP)
     .                                     -2d0*kq/lq*muq*DDSIN(phi0))
     .     -xiSQ*DDSIN(phiSQ)-muq/lq*MSPQ*DDSIN(phiSPq)
     .     -muPQ*kq*(muq/lq)**2*DDSIN(phi02-phiP)
     .     -XIFQ*(muPQ*DDSIN(phiP-phiF)
     .     +2d0*kq/lq*muq*DDSIN(phi02-phiF)))
      ELSE
         IAk=0d0
         aux1=XISQ*DDCOS(PhiSq)
         IXIS=lq*vuq*vdq*(IAl-muPQ*DDSIN(phi01-phiP))
     .       -muq/lq*MSPQ*DDSIN(phiSPQ)-XIFQ*muPQ*DDSIN(phiP-phiF)
      ENDIF

c      10) Squark and slepton
*   In order to run the squark masses from Q2 to QSTSB
*   including all thresholds, tree level expressions for the soft
*   Higgs masses and various logarithms are needed:

      AT=PAR(12)
      AB=PAR(13)

      MHuS=(MUQ*(Alcos1+MUPQ*DDCOS(phi01-phiP)+kq/lq*muq*DDCOS(Phi0))
     .         +M3HQ*DDCOS(Phi3q)+LQ*XIFQ*DDCOS(phi01-phiF))/TANBQ
     .      -LQ**2*VDQ**2-MUQ**2+GQ/2d0*(VDQ**2-VUQ**2)
      MHdS=(MUQ*(Alcos1+MUPQ*DDCOS(phi01-phiP)+kq/lq*muq*DDCOS(Phi0))
     .         +M3HQ*DDCOS(Phi3q)+LQ*XIFQ*DDCOS(phi01-phiF))*TANBQ
     .      -LQ**2*VUQ**2-MUQ**2+GQ/2d0*(VUQ**2-VDQ**2)
      MSS=LQ**2*vuq*vdq/muq
     .        *(Alcos1+2d0*kq/lq*muq*DDCOS(phi0)+mupq*DDCOS(Phi01-PhiP))
     . -kq/lq*muq*(Akcos2+2d0*kq/lq*muq)-lq**2*(vuq**2+vdq**2)
     .  -XIFQ*(2d0*KQ*DDCOS(phi02-PhiF)+LQ*MUPQ*DDCOS(PhiP-PhiF)/MUQ)
     .  -MUPQ**2-3d0*MUPQ*KQ/LQ*MUQ*DDCOS(Phi02-PhiP)-MSPQ*DDCOS(PhiSPq)
     .  -LQ*XISQ/MUQ*DDCOS(PhiSq)

      MQ3=PAR(7)
      MU3=PAR(8)
      MD3=PAR(9)
      ML3=PAR(10)
      ME3=PAR(11)
      MQ=PAR(15)
      MU1=PAR(16)
      MD=PAR(17)
      ML=PAR(18)
      ME=PAR(19)

      CALL GETSUSYCOUP(PAR)

      XI=g1S*(MHuS-MHdS+MQ3+2d0*MQ-2d0*(MU3+2d0*MU1)+MD3+2d0*MD
     .   +ME3+2d0*ME-ML3+2d0*ML)
      HT2=HTOPS**2
      HB2=HBOTS**2
      HTAU2=HTAUS**2
      MH1=MHuS
      MH2=MHdS
      L2=PAR(1)**2
      K2=PAR(2)**2

      SP= HT2*(-3d0*MH1 - MQ3 + 4d0*MU3)
     .  + HB2*(3d0*MH2 - MQ3 - 2d0*MD3)
     .  + HTAU2*(MH2 + ML3 - 2d0*ME3) + L2*(MH2 - MH1)
     .  + (G1S/18d0 + 3d0/2d0*G2S + 8d0/3d0*G3S)*(MQ3+2d0*MQ)
     .  - (16d0/9d0*G1S + 16d0/3d0*G3S)*(MU3+2d0*MU1)
     .  + (2d0/9d0*G1S + 8d0/3d0*G3S)*(MD3+2d0*MD)
     .  + (G1S/2d0 + 3d0/2d0*G2S)*(MH1-MH2-(ML3+2d0*ML)) 
     .  + 2d0*G1S*(ME3+2d0*ME)
      SIG1= G1S*(MH1 + MH2 + (MQ3+2d0*MQ)/3d0 + 8d0/3d0*(MU3+2d0*MU1) 
     .    + 2d0/3d0*(MD3+2d0*MD)+ ML3+2d0*ML + 2d0*(ME3+2d0*ME))
      SIG2=G2S*(MH1 + MH2 + 3d0*(MQ3+2d0*MQ) + ML3+2d0*ML)
      SIG3=G3S*(2d0*(MQ3+2d0*MQ) + MU3+2d0*MU1 + MD3+2d0*MD)

      F20= HT2*(MH1+MQ3+MU3+AT**2) + HB2*(MH2+MQ3+MD3+AB**2)
     .       - g1S*M1**2/9d0 - 3d0*g2S*M2**2
     .       - 16d0/3d0*g3S*M3**2 + XI/6d0
     .       + (-10d0*HT2**2*(MH1+MQ3+MU3+2d0*AT**2)
     .       - 10d0*HB2**2*(MH2+MQ3+MD3+2d0*AB**2)
     .       - 10d0*HB2**4*(MH2+MQ3+MD3+2d0*AB**2)
     .       - HB2*HTAU2*(2d0*MH2+MQ3+MD3+ML3+ME3+AB**2+ATAU**2
     .                         +2d0*AB*ATAU*DDCOS(PhiAb-PhiAtau))
     .       - L2*HT2*(2d0*MH1+MH2+MSS+MQ3+MU3+AT**2+Alcos1**2+IAl**2
     .                    +2d0*AT*(Alcos1*DDCOS(PhiAT)
     .                    +IAl*DDSIN(PhiAT)))
     .       - L2*HB2*(MH1+2d0*MH2+MSS+MQ3+MD3+AB**2+Alcos1**2+IAl**2
     .                    +2d0*AB*(Alcos1*DDCOS(PhiAB)
     .                    +IAl*DDSIN(PhiAB)))
     .       + 4d0/3d0*G1S*HT2*(MH1+MQ3+MU3+AT**2+2d0*M1**2
     .                          -2d0*M1*AT*DDCOS(PhiAT-PhiM1))
     .       + 2d0/3d0*G1S*HB2*(MH2+MQ3+MD3+AB**2+2d0*M1**2
     .                          -2d0*M1*AB*DDCOS(PhiAB-PhiM1))
     .       + 199d0/54d0*G1S**2*M1**2 + 33d0/2d0*G2S**2*M2**2
     .     - 64d0/3d0*G3S**2*M3**2
     .     + 1d0/3d0*G1S*G2S*(M1**2+M2**2+M1*M2*DDCOS(PhiM1-PhiM2))
     .       + 16d0/27d0*G1S*G3S*(M1**2+M3**2+M1*M3*DDCOS(PhiM1-PhiM3))
     .     + 16d0*G2S*G3S*(M2**2+M3**2+M2*M3*DDCOS(PhiM1-PhiM2))
     .       + 1d0/3d0*G1S*SP + 1d0/18d0*G1S*SIG1
     .       + 3d0/2d0*G2S*SIG2 + 8d0/3d0*G3S*SIG3)/16d0/Pi**2

      F21= 2d0*HT2*(MH1+MQ3+MU3+AT**2)
     .       - 16d0/9d0*g1S*M1**2 - 16d0/3d0*g3S*M3**2 - 2d0*XI/3d0
     .       + (-16d0*HT2**2*(MH1+MQ3+MU3+2d0*AT**2)
     .       - 2d0*HT2*HB2*(MH1+MH2+2d0*MQ3+MU3+MD3+AT**2+AB**2
     .                           +2d0*AT*AB*DDCOS(PhiAT-PhiAB))
     .       - 2d0*HT2*L2*(2d0*MH1+MH2+MSS+MQ3+MU3+AT**2+Alcos1**2
     .             +IAL**2+2d0*AT*(Alcos1*DDCOS(PhiAT)
     .             +IAl*DDSIN(PhiAT)))
     .       - 2d0/3d0*G1S*HT2*(MH1+MQ3+MU3+AT**2+2d0*M1**2
     .                          -2d0*AT*M1*DDCOS(PhiAT-PhiM1))
     .       + 6d0*G2S*HT2*(MH1+MQ3+MU3+AT**2+2d0*M2**2
     .                          -2d0*AT*M2*DDCOS(PhiAT-PhiM2))
     .       + 1712d0/27d0*G1S**2*M1**2 - 64d0/3d0*G3S**2*M3**2
     .       + 256d0/27d0*G1S*G3S*(M1**2+M3**2+M1*M3*DDCOS(PhiM1-PhiM3))
     .     - 4d0/3d0*G1S*SP + 8d0/9d0*G1S*SIG1 + 8d0/3d0*G3S*SIG3)
     .                                /16d0/Pi**2

      F22= 2d0*HB2*(MH2+MQ3+MD3+AB**2)
     .       - 4d0/9d0*g1S*M1**2 - 16d0/3d0*g3S*M3**2 + XI/3d0
     .       + (-16d0*HB2**2*(MH2+MQ3+MD3+2d0*AB**2)
     .       - 2d0*HB2*HT2*(MH1+MH2+2d0*MQ3+MU3+MD3+AB**2+AT**2
     .                  +2d0*AT*AB*DDCOS(PhiAB-PhiAT))
     .       - 2d0*HB2*HTAU2*(2d0*MH2+MQ3+MD3+ML3+ME3+AB**2+ATAU**2
     .                  +2d0*AB*ATAU*DDCOS(PhiAB-PhiAtau))
     .       - 2d0*HB2*L2*(MH1+2d0*MH2+MSS+MQ3+MD3+AB**2+Alcos1**2
     .            +IAL**2+2d0*AB*(Alcos1*DDCOS(PhiAB)+IAl*DDSIN(PhiAB)))
     .       + 2d0/3d0*G1S*HB2*(MH2+MQ3+MD3+AB**2+2d0*M1**2
     .                            -2d0*M1*AB*DDCOS(PhiAB-PhiM1))
     .       + 6d0*G2S*HB2*(MH2+MQ3+MD3+AB**2+2d0*M2**2
     .                            -2d0*M2*AB*DDCOS(PhiAB-PhiM2))
     .       + 404d0/27d0*G1S**2*M1**2 - 64d0/3d0*G3S**2*M3**2
     .       + 64d0/27d0*G1S*G3S*(M1**2+M3**2+M1*M3*DDCOS(PhiM1-PhiM3))
     .     + 2d0/3d0*G1S*SP + 2d0/9d0*G1S*SIG1 + 8d0/3d0*G3S*SIG3)
     .                                /16d0/Pi**2

        MSQ3=MQ3-dlog(Q2/QSTSB)*F20/16d0/Pi**2
        MSU3=MU3-dlog(Q2/QSTSB)*F21/16d0/Pi**2
        MSD3=MD3-dlog(Q2/QSTSB)*F22/16d0/Pi**2

      aux1=AT*DDCOS(PhiAT)+
     .   (6d0*HT2*AT*DDCOS(PhiAt)+HB2*AB*DDCOS(PhiAb)
     .     +L2*(Alcos1*DDCOS(phi01)+IAl*DDSIN(Phi01))
     .     +13d0/9d0*g1S*M1*DDCOS(PhiM1)+3d0*g2S*M2*DDCOS(PhiM2)
     .     +16d0/3d0*g3S*M3*DDCOS(PhiM3)
     .     +(-44d0*HT2**2*AT*DDCOS(PhiAt)
     .       -5d0*HT2*HB2*(AT*DDCOS(PhiAt)+AB*DDCOS(PhiAb))
     .       -3d0*HT2*L2
     .         *(AT*DDCOS(PhiAt)+Alcos1*DDCOS(phi01)+IAl*DDSIN(Phi01))
     .       -10d0*HB2**2*AB*DDCOS(PhiAb)
     .       -HB2*HTAU2*(AB*DDCOS(PhiAb)+ATAU*DDCOS(PhiAtau))
     .       -4d0*HB2*L2
     .         *(AB*DDCOS(PhiAb)+Alcos1*DDCOS(phi01)+IAl*DDSIN(Phi01))
     .       -6d0*L2**2*(Alcos1*DDCOS(phi01)+IAl*DDSIN(Phi01))
     .       -L2*HTAU2
     .         *(Alcos1*DDCOS(phi01)+IAl*DDSIN(Phi01)
     .          +ATAU*DDCOS(PhiAtau))
     .       -2d0*L2*K2*(Alcos1*DDCOS(phi01)+IAl*DDSIN(Phi01)
     .                  +Akcos2*DDCOS(phi02)+IAk*DDSIN(Phi02))
     .       +2d0*g1S*HT2*(AT*DDCOS(PhiAt)-M1*DDCOS(PhiM1))
     .       +6d0*g2S*HT2*(AT*DDCOS(PhiAt)-M2*DDCOS(PhiM2))
     .       +16d0*g3S*HT2*(AT*DDCOS(PhiAt)-M3*DDCOS(PhiM3))
     .       +2d0/3d0*g1S*HB2*(AB*DDCOS(PhiAb)-M1*DDCOS(PhiM1))
     .       -2743d0/81d0*g1S**2*M1*DDCOS(PhiM1)
     .       -5d0/3d0*g1S*g2S*(M1*DDCOS(PhiM1)+M2*DDCOS(PhiM2))
     .       -136d0/27d0*g1S*g3S*(M1*DDCOS(PhiM1)+M3*DDCOS(PhiM3))
     .       -15d0*g2S**2*M2*DDCOS(PhiM2)
     .       -8d0*g2S*g3S*(M2*DDCOS(PhiM2)+M3*DDCOS(PhiM3))
     .       +32d0/9d0*g3S**2*M3*DDCOS(PhiM3))/16d0/Pi**2)
     .             /16d0/Pi**2*dlog(QSTSB/Q2)
     
      aux2=AT*DDSIN(PhiAT)+
     .   (6d0*HT2*AT*DDSIN(PhiAt)+HB2*AB*DDSIN(PhiAB)
     .     +L2*(IAl*DDCOS(phi01)-Alcos1*DDSIN(Phi01))                                   
     .     +13d0/9d0*g1S*M1*DDSIN(PhiM1)+3d0*g2S*M2*DDSIN(PhiM2)
     .     +16d0/3d0*g3S*M3*DDSIN(PhiM3)
     .     +(-44d0*HT2**2*AT*DDSIN(PhiAt)
     .       -5d0*HT2*HB2*(AT*DDSIN(PhiAt)+AB*DDSIN(PhiAb))
     .       -3d0*HT2*L2
     .         *(AT*DDSIN(PhiAt)+IAl*DDCOS(phi01)-Alcos1*DDSIN(Phi01))
     .       -10d0*HB2**2*AB*DDSIN(PhiAb)
     .       -HB2*HTAU2*(AB*DDSIN(PhiAb)+ATAU*DDSIN(PhiAtau))
     .       -4d0*HB2*L2
     .         *(AB*DDSIN(PhiAb)+IAl*DDCOS(phi01)-Alcos1*DDSIN(Phi01))
     .       -6d0*L2**2*(IAl*DDCOS(phi01)-Alcos1*DDSIN(Phi01))
     .       -L2*HTAU2
     .         *(IAl*DDCOS(phi01)-Alcos1*DDSIN(Phi01)
     .          +ATAU*DDSIN(PhiAtau))
     .       -2d0*L2*K2*(IAl*DDCOS(phi01)-Alcos1*DDSIN(Phi01)
     .                  +IAk*DDCOS(phi02)-Akcos2*DDSIN(Phi02))
     .       +2d0*g1S*HT2*(AT*DDSIN(PhiAt)-M1*DDSIN(PhiM1))
     .       +6d0*g2S*HT2*(AT*DDSIN(PhiAt)-M2*DDSIN(PhiM2))
     .       +16d0*g3S*HT2*(AT*DDSIN(PhiAt)-M3*DDSIN(PhiM3))
     .       +2d0/3d0*g1S*HB2*(AB*DDSIN(PhiAb)-M1*DDSIN(PhiM1))
     .       -2743d0/81d0*g1S**2*M1*DDSIN(PhiM1)
     .       -5d0/3d0*g1S*g2S*(M1*DDSIN(PhiM1)+M2*DDSIN(PhiM2))
     .       -136d0/27d0*g1S*g3S*(M1*DDSIN(PhiM1)+M3*DDSIN(PhiM3))
     .       -15d0*g2S**2*M2*DDSIN(PhiM2)
     .       -8d0*g2S*g3S*(M2*DDSIN(PhiM2)+M3*DDSIN(PhiM3))
     .       +32d0/9d0*g3S**2*M3*DDSIN(PhiM3))/16d0/Pi**2)
     .             /16d0/Pi**2*dlog(QSTSB/Q2)

      IF(aux1.ne.0d0)THEN
       phiATQ=datan(aux2/aux1)
       IF(aux1.lt.0d0)PhiATQ=PhiATQ+Pi
      ELSEIF(aux2.gt.0d0)THEN
       phiATQ=Pi/2d0
      ELSEIF(aux2.lt.0d0)THEN
       phiATQ=-Pi/2d0
      ELSE
       phiATQ=0d0
      ENDIF

        ATP=dsqrt(aux1**2+aux2**2)

      aux1=AB*DDCOS(PhiAB)+
     .   (6d0*HB2**2*AB*DDCOS(PhiAB)+HT2*AT*DDCOS(PhiAT)
     .      +L2*(Alcos1*DDCOS(phi01)+IAl*DDSIN(Phi01))
     .      +7d0/9d0*G1S*M1*DDCOS(PhiM1)+3d0*G2S*M2*DDCOS(PhiM2)
     .      +16d0/3d0*g3S*M3*DDCOS(PhiM3)
     .      +(-44d0*HB2**2*AB*DDCOS(PhiAB)
     .      -5d0*HT2*HB2*(AT*DDCOS(PhiAt)+AB*DDCOS(PhiAB))
     .      -3d0*HB2*HTAU2*(AB*DDCOS(PhiAB)+ATAU*DDCOS(PhiAtau))
     .      -3d0*HB2*L2
     .           *(AB*DDCOS(PhiAB)+Alcos1*DDCOS(phi01)+IAl*DDSIN(Phi01))
     .      -10d0*HT2**2*AT*DDCOS(PhiAT)
     .      -4d0*HT2*L2
     .           *(AT*DDCOS(PhiAT)+Alcos1*DDCOS(phi01)+IAl*DDSIN(Phi01))
     .      -6d0*HTAU2**2*ATAU*DDCOS(PhiAtau)
     .      -6d0*L2**2*(Alcos1*DDCOS(phi01)+IAl*DDSIN(Phi01))
     .      -2d0*L2*K2*(Alcos1*DDCOS(phi01)+IAl*DDSIN(Phi01)
     .                     +Akcos2*DDCOS(phi02)+IAk*DDSIN(Phi02))
     .      +2d0/3d0*g1S*HB2*(AB*DDCOS(PhiAB)-M1*DDCOS(PhiM1))
     .      +6d0*g2S*HB2*(AB*DDCOS(PhiAB)-M2*DDCOS(PhiM2))
     .      +16d0*g3S*HB2*(AB*DDCOS(PhiAB)-M3*DDCOS(PhiM3))
     .      +4d0/3d0*g1S*HT2*(AT*DDCOS(PhiAT)-M1*DDCOS(PhiM1))
     .      +2d0*HTAU2*g1S*(ATAU*DDCOS(PhiAtau)-M1*DDCOS(PhiM1))
     .      -1435d0/81d0*g1S**2*M1*DDCOS(PhiM1)
     .      -5d0/3d0*g1S*g2S*(M1*DDCOS(PhiM1)+M2*DDCOS(PhiM2))
     .      -40d0/27d0*g1S*g3S*(M1*DDCOS(PhiM1)+M3*DDCOS(PhiM3))
     .      -15d0*g2S**2*M2*DDCOS(PhiM2)
     .      -8d0*g2S*g3S*(M2*DDCOS(PhiM2)+M3*DDCOS(PhiM3))
     .      +32d0/9d0*g3S**2*M3*DDCOS(PhiM3))/16d0/Pi**2)
     .             /16d0/Pi**2*dlog(QSTSB/Q2)
     
      aux2=AB*DDSIN(PhiAB)+
     .   (6d0*HB2*AB*DDSIN(PhiAB)+HT2*AT*DDSIN(PhiAT)
     .      +L2*(IAl*DDCOS(phi01)-Alcos1*DDSIN(Phi01))                                  
     .      +7d0/9d0*g1S*M1*DDSIN(PhiM1)+3d0*g2S*M2*DDSIN(PhiM2)
     .      +16d0/3d0*g3S*M3*DDSIN(PhiM3)
     .      +(-44d0*HB2**2*AB*DDSIN(PhiAB)
     .      -5d0*HT2*HB2*(AT*DDSIN(PhiAt)+AB*DDSIN(PhiAB))
     .      -3d0*HB2*HTAU2*(AB*DDSIN(PhiAB)+ATAU*DDSIN(PhiAtau))
     .      -3d0*HB2*L2
     .           *(AB*DDSIN(PhiAB)+IAl*DDCOS(phi01)-Alcos1*DDSIN(Phi01))
     .      -10d0*HT2**2*AT*DDSIN(PhiAT)
     .      -4d0*HT2*L2
     .           *(AT*DDSIN(PhiAT)+IAl*DDCOS(phi01)-Alcos1*DDSIN(Phi01))
     .      -6d0*HTAU2**2*ATAU*DDSIN(PhiAtau)
     .      -6d0*L2**2*(IAl*DDCOS(phi01)-Alcos1*DDSIN(Phi01))
     .      -2d0*L2*K2*(IAl*DDCOS(phi01)-Alcos1*DDSIN(Phi01)
     .                     +IAk*DDCOS(phi02)-Akcos2*DDSIN(Phi02))
     .      +2d0/3d0*g1S*HB2*(AB*DDSIN(PhiAB)-M1*DDSIN(PhiM1))
     .      +6d0*g2S*HB2*(AB*DDSIN(PhiAB)-M2*DDSIN(PhiM2))
     .      +16d0*g3S*HB2*(AB*DDSIN(PhiAB)-M3*DDSIN(PhiM3))
     .      +4d0/3d0*g1S*HT2*(AT*DDSIN(PhiAT)-M1*DDSIN(PhiM1))
     .      +2d0*HTAU2*g1S*(ATAU*DDSIN(PhiAtau)-M1*DDSIN(PhiM1))
     .      -1435d0/81d0*g1S**2*M1*DDSIN(PhiM1)
     .      -5d0/3d0*g1S*g2S*(M1*DDSIN(PhiM1)+M2*DDSIN(PhiM2))
     .      -40d0/27d0*g1S*g3S*(M1*DDSIN(PhiM1)+M3*DDSIN(PhiM3))
     .      -15d0*g2S**2*M2*DDSIN(PhiM2)
     .      -8d0*g2S*g3S*(M2*DDSIN(PhiM2)+M3*DDSIN(PhiM3))
     .      +32d0/9d0*g3S**2*M3*DDSIN(PhiM3))/16d0/Pi**2)
     .             /16d0/Pi**2*dlog(QSTSB/Q2)

      IF(aux1.ne.0d0)THEN
       phiABQ=datan(aux2/aux1)
       IF(aux1.lt.0d0)PhiABQ=PhiABQ+Pi
      ELSEIF(aux2.gt.0d0)THEN
       phiABQ=Pi/2d0
      ELSEIF(aux2.lt.0d0)THEN
       phiABQ=-Pi/2d0
      ELSE
       phiABQ=0d0
      ENDIF

        ABP=dsqrt(aux1**2+aux2**2)

!      write(0,*) "init_CPV,bottom:"
!      write(0,*) "RENSCALE",Q2,"STSBSCALE",QSTSB
!      write(0,*) "GAUGE",ALSMZ,ALEMMZ,GF,g1,g2,S2TW
!      write(0,*) "QGAUGE",G1Q,G2Q,GQ,ALSQ
!      write(0,*) "QHIGGS",ZHU,ZHD,ZS,vuq,vdq,TANBQ
!      write(0,*) "QQUARK",Ytq,Ybq,MTOPQ,MBOTQ
!      write(0,*) "QPAR",lq,kq,Alcos1,Akcos2,muq,NUQ
!      write(0,*) Ytq**2,muH2,mtopq**2,QSTSB
!      write(0,*) g1q+g2q,MZ**2,sinb**2
!      write(0,*) g2q,MW**2
!      write(0,*) g1q,mur**2,M1r**2,mur**2,M2r**2
!      write(0,*) YTQ,vuQ

      RETURN
      END
