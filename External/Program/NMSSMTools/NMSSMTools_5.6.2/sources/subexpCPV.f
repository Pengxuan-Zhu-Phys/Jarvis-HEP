      SUBROUTINE CONSTSUSYPART_CPV(PROB)

*   Subroutine to check experimental constraints on SUSY particles
*
*   The required data files and numbers (inv. Z width, lower bounds on
*   sparticle masses) are transferred via COMMON/LEP/...

*      PROB(1) =/= 0   chargino too light
*      PROB(2) =/= 0   excluded by Z -> neutralinos
*      PROB(20) =/= 0  excluded by stop -> b l sneutrino
*      PROB(21) =/= 0  excluded by stop -> neutralino c
*      PROB(22) =/= 0  excluded by sbottom -> neutralino b
*      PROB(23) =/= 0  squark/gluino too light
*      PROB(24) =/= 0  selectron/smuon too light
*      PROB(25) =/= 0  stau too light


      IMPLICIT NONE

      INTEGER I,J
      INTEGER NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg
      INTEGER NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b
      INTEGER NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau
      INTEGER Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1
      INTEGER Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1
      INTEGER Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1
      INTEGER Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1
      INTEGER Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1
      INTEGER Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1
      INTEGER Nstblsn,Nstnc,Nsbnb,Nglsq

      DOUBLE PRECISION PROB(*),PI,MMIN
      DOUBLE PRECISION LAMBDA,X,Y,Z,O,Q,DZ,E1,E2,SIG,S,SQRS,CONV

      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION MGL
      DOUBLE PRECISION MCH2(2),U(2,2,2),V(2,2,2)
      DOUBLE PRECISION MNEU(5),NEU(5,5,2)
      DOUBLE PRECISION MST2(2),UT(2,2,2),MSB2(2),UB(2,2,2),MSL2(2),
     . UTAU(2,2,2),MSNT2
      DOUBLE PRECISION MSU2(2),MSD2(2),MSE2(2),MSNE2,MSMU2(2),
     . UMU(2,2,2)
      DOUBLE PRECISION MST2P(2),MSB2P(2),MSU2P(2),MSD2P(2)
      DOUBLE PRECISION GZMAX,GZINV,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1
      DOUBLE PRECISION SIGNEU,MSLMIN,MSTMIN,MSQMIN,MGLMIN
      DOUBLE PRECISION hZind(1000,2),hZbb(1000,2),hZll(1000,2)
      DOUBLE PRECISION hZinv(1000,2),hZjj(1000,2),hZgg(1000,2)
      DOUBLE PRECISION hA4b(10000,3),hA4tau(10000,3)
      DOUBLE PRECISION hA2b2tau(10000,3),hA2tau2b(10000,3)
      DOUBLE PRECISION AAA6b(10000,3),AAA6tau(10000,3)
      DOUBLE PRECISION AAZ4b(10000,3),AAZ4tau(10000,3)
      DOUBLE PRECISION AAZ2b2tau(10000,3)
      DOUBLE PRECISION cccc02(100,2),cccc04(100,2),cccc05(100,2)
      DOUBLE PRECISION cccc06(100,2),cccc08(100,2),cccc1(100,2)
      DOUBLE PRECISION ccgg02(100,2),ccgg04(100,2),ccgg05(100,2)
      DOUBLE PRECISION ccgg06(100,2),ccgg08(100,2),ccgg1(100,2)
      DOUBLE PRECISION cctt02(100,2),cctt04(100,2),cctt05(100,2)
      DOUBLE PRECISION cctt06(100,2),cctt08(100,2),cctt1(100,2)
      DOUBLE PRECISION gggg02(100,2),gggg04(100,2),gggg05(100,2)
      DOUBLE PRECISION gggg06(100,2),gggg08(100,2),gggg1(100,2)
      DOUBLE PRECISION ttgg02(100,2),ttgg04(100,2),ttgg05(100,2)
      DOUBLE PRECISION ttgg06(100,2),ttgg08(100,2),ttgg1(100,2)
      DOUBLE PRECISION tttt02(100,2),tttt04(100,2),tttt05(100,2)
      DOUBLE PRECISION tttt06(100,2),tttt08(100,2),tttt1(100,2)
      DOUBLE PRECISION stblsn(100,2),stnc(100,2),sbnb(100,2)
      DOUBLE PRECISION glsq(100,2)

      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/GLUSPEC/MGL
      COMMON/CHASPEC/MCH2,U,V
      COMMON/NEUSPEC/MNEU,NEU
      COMMON/SFERM3SPEC/MST2,UT,MSB2,UB,MSL2,UTAU,MSNT2
      COMMON/SFERM1SPEC/MSU2,MSD2,MSE2,MSNE2,MSMU2,UMU
      COMMON/SFERMPSPEC/MST2P,MSB2P,MSU2P,MSD2P
      COMMON/LEP/GZMAX,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1,SIGNEU,
     .      MSLMIN,MSTMIN,MSQMIN,MGLMIN,
     .      hZind,hZbb,hZll,hZinv,hZjj,hZgg,
     .      hA4b,hA4tau,hA2b2tau,hA2tau2b,
     .      AAA6b,AAA6tau,AAZ4b,AAZ4tau,AAZ2b2tau,
     .      cccc02,cccc04,cccc05,cccc06,cccc08,cccc1,
     .      ccgg02,ccgg04,ccgg05,ccgg06,ccgg08,ccgg1,
     .      cctt02,cctt04,cctt05,cctt06,cctt08,cctt1,
     .      gggg02,gggg04,gggg05,gggg06,gggg08,gggg1,
     .      ttgg02,ttgg04,ttgg05,ttgg06,ttgg08,ttgg1,
     .      tttt02,tttt04,tttt05,tttt06,tttt08,tttt1,
     .      stblsn,stnc,sbnb,glsq,
     .      NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg,
     .      NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b,
     .      NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau,
     .      Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1,
     .      Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1,
     .      Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1,
     .      Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1,
     .      Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1,
     .      Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1,
     .      Nstblsn,Nstnc,Nsbnb,Nglsq

c********************************************************************************

      LAMBDA(X,Y,Z)= DSQRT(X**2+Y**2+Z**2-2d0*X*Y-2d0*X*Z-2d0*Y*Z)

      PI=4d0*DATAN(1d0)

c      A: LEP searches

* Test on stop -> b l sneutrino

      I=1
      DOWHILE(stblsn(I,1).LE.dsqrt(MST2P(1)).AND. I.LT.Nstblsn)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MST2P(1)).LT.stblsn(Nstblsn,1))THEN
       MMIN=stblsn(I-1,2)+(dsqrt(MST2P(1))-stblsn(I-1,1))
     .  /(stblsn(I,1)-stblsn(I-1,1))*(stblsn(I,2)-stblsn(I-1,2))
       PROB(20)=DDIM(1d0,dsqrt(MSNE2)/MMIN)
      ENDIF

* Test on stop -> neutralino c

      I=1
      DOWHILE(stnc(I,1).LE.dsqrt(MNEU(1)) .AND. I.LT.Nstnc)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MNEU(1)).LT.stnc(Nstnc,1))THEN
       MMIN=stnc(I-1,2)+(dsqrt(MNEU(1))-stnc(I-1,1))
     .  /(stnc(I,1)-stnc(I-1,1))*(stnc(I,2)-stnc(I-1,2))
       PROB(21)=DDIM(1d0,dsqrt(MST2P(1))/MMIN)
      ENDIF

* Test on sbottom -> neutralino b

      I=1
      DOWHILE(sbnb(I,1).LE.dsqrt(MSB2P(1)) .AND. I.LT.Nsbnb)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. dsqrt(MSB2P(1)).LT.sbnb(Nsbnb,1))THEN
       MMIN=sbnb(I-1,2)+(dsqrt(MSB2P(1))-sbnb(I-1,1))
     .  /(sbnb(I,1)-sbnb(I-1,1))*(sbnb(I,2)-sbnb(I-1,2))
       PROB(22)=DDIM(1d0,dsqrt(MNEU(1))/MMIN)
      ENDIF

* Test on gluino/squark masses

      I=1
      DOWHILE(glsq(I,1).LE.DABS(MGL) .AND. I.LT.Nglsq)
       I=I+1
      ENDDO
      MMIN=MSQMIN
      IF(I.GT.1 .AND. DABS(MGL).LT.glsq(Nglsq,1))THEN
       MMIN=glsq(I-1,2)+(DABS(MGL)-glsq(I-1,1))
     .  /(glsq(I,1)-glsq(I-1,1))*(glsq(I,2)-glsq(I-1,2))
      ENDIF
      PROB(23)=DDIM(1d0,
     .          dsqrt(MIN(MSU2P(1),MSU2P(2),MSD2P(1),MSD2P(2)))/MMIN)
     .       +DDIM(1d0,DABS(MGL)/MGLMIN)

* Test on slepton masses

      PROB(24)=DDIM(1d0,dsqrt(MIN(MSE2(1),MSE2(2)))/MSLMIN)
      PROB(25)=DDIM(1d0,dsqrt(MSL2(1))/MSTMIN)

* Test on chargino mass

      PROB(1)=DDIM(1d0,dsqrt(MCH2(1))/MCHAMIN)

* Test on Z width into neutralinos

      IF(dsqrt(MNEU(1)).LT.MZ/2d0)THEN
       GZINV=MZ**3*GF/(12d0*DSQRT(2d0)*PI)
     .   *(NEU(1,3,1)**2-NEU(1,4,1)**2+NEU(1,3,2)**2-NEU(1,4,2)**2)**2
     .   *(1d0-4d0*MNEU(1)/MZ**2)**1.5d0
        PROB(2)=DDIM(GZINV/GZINVMAX,1d0)
      ENDIF

* Test on neutralinos

      SQRS=209d0
      S=SQRS**2
      CONV=.3894D9
      DZ = 1d0/(S-MZ**2)**2

      DO I=1,5
       DO J=I,5
        IF(dsqrt(MNEU(I))+dsqrt(MNEU(J)).LT.SQRS)THEN
         O= (NEU(I,3,1)*NEU(J,3,1)+NEU(I,3,2)*NEU(J,3,2)
     .      -NEU(I,4,1)*NEU(J,4,1)-NEU(I,4,2)*NEU(J,4,2))**2
     .     +(NEU(I,3,1)*NEU(J,3,2)-NEU(I,3,2)*NEU(J,3,1)
     .      -NEU(I,4,1)*NEU(J,4,2)+NEU(I,4,2)*NEU(J,4,1))**2
         E1= (S+MNEU(I)-MNEU(J))/(2d0*SQRS)
         E2= (S+MNEU(J)-MNEU(I))/(2d0*SQRS)
         Q= LAMBDA(S,MNEU(I),MNEU(J))/(2d0*SQRS)
         SIG= 1d0/(4d0*PI)*Q/SQRS*DZ
     .       *(E1*E2+Q**2/3d0-dsqrt(MNEU(I)*MNEU(J)))
     .       *(g1+g2)**2/4d0*O*(.25d0-S2TW+2d0*S2TW**2)
         SIG= CONV*SIG
         IF(I.EQ.J)SIG= SIG/2d0
         IF(I.EQ.1)THEN
          IF(J.GT.1)PROB(2)=PROB(2)+DDIM(SIG/SIGNEU1,1d0)
         ELSE
          PROB(2)=PROB(2)+DDIM(SIG/SIGNEU,1d0)
         ENDIF
        ENDIF
       ENDDO
      ENDDO


c      B: LHC limits -> not implemented

      END


      SUBROUTINE LEP_HIGGS_CPV(PROB)

***********************************************************************
*   Subroutine to check Higgs LEP constraints:
*      PROB(3) =/= 0   charged Higgs too light
*      PROB(4) =/= 0   excluded by ee -> hZ 
*      PROB(5) =/= 0   excluded by ee -> hZ, h -> bb
*      PROB(6) =/= 0   excluded by ee -> hZ, h -> tautau
*      PROB(7) =/= 0   excluded by ee -> hZ, h -> invisible 
*      PROB(8) =/= 0   excluded by ee -> hZ, h -> 2jets
*      PROB(9) =/= 0   excluded by ee -> hZ, h -> 2photons
*      PROB(10) =/= 0  excluded by ee -> hZ, h -> AA -> 4bs
*      PROB(11) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus
*      PROB(12) =/= 0  excluded by ee -> hZ, h -> AA -> 2bs 2taus
*      PROB(13) =/= 0  excluded by Z -> hA (Z width)
*      PROB(14) =/= 0  excluded by ee -> hA -> 4bs
*      PROB(15) =/= 0  excluded by ee -> hA -> 4taus
*      PROB(16) =/= 0  excluded by ee -> hA -> 2bs 2taus
*      PROB(17) =/= 0  excluded by ee -> hA -> AAA -> 6bs
*      PROB(18) =/= 0  excluded by ee -> hA -> AAA -> 6taus
*      PROB(19) =/= 0  excluded by ee -> Zh -> ZAA -> Z + light pairs
*      PROB(41) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus (ALEPH analysis)
***********************************************************************

      IMPLICIT NONE

      INTEGER I,J,K
      INTEGER NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg
      INTEGER NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b
      INTEGER NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau
      INTEGER Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1
      INTEGER Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1
      INTEGER Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1
      INTEGER Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1
      INTEGER Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1
      INTEGER Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1
      INTEGER Nstblsn,Nstnc,Nsbnb,Nglsq

      DOUBLE PRECISION PROB(*),PI,S,SQRS,SQR2,CONV
      DOUBLE PRECISION LAMBDA,X,Y,Z,MA,MH,BRINV,DMA,Rmax
      DOUBLE PRECISION MAtest,ceff,HMIN,HMAX,M1,M2,R

      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd,CMASS
      DOUBLE PRECISION MHC2,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION WIDTH(5),HCWIDTH
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),
     . BRCC(5),BRBB(5),BRTT(5),BRWW(5),BRZZ(5),BRGG(5),BRZG(5)
      DOUBLE PRECISION BRHHH(5,10),BRHCHC(5),BRHAZ(5,4),BRHCW(5),
     . BRHIGGS(5)
      DOUBLE PRECISION BRNEU(5,5,5),BRCHA(5,3),BRHSQ(5,10),BRHSL(5,7),
     . BRSUSY(5)
      DOUBLE PRECISION GZ,GZMAX,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1
      DOUBLE PRECISION SIGNEU,MSLMIN,MSTMIN,MSQMIN,MGLMIN
      DOUBLE PRECISION hZind(1000,2),hZbb(1000,2),hZll(1000,2)
      DOUBLE PRECISION hZinv(1000,2),hZjj(1000,2),hZgg(1000,2)
      DOUBLE PRECISION hA4b(10000,3),hA4tau(10000,3)
      DOUBLE PRECISION hA2b2tau(10000,3),hA2tau2b(10000,3)
      DOUBLE PRECISION AAA6b(10000,3),AAA6tau(10000,3)
      DOUBLE PRECISION AAZ4b(10000,3),AAZ4tau(10000,3)
      DOUBLE PRECISION AAZ2b2tau(10000,3)
      DOUBLE PRECISION cccc02(100,2),cccc04(100,2),cccc05(100,2)
      DOUBLE PRECISION cccc06(100,2),cccc08(100,2),cccc1(100,2)
      DOUBLE PRECISION ccgg02(100,2),ccgg04(100,2),ccgg05(100,2)
      DOUBLE PRECISION ccgg06(100,2),ccgg08(100,2),ccgg1(100,2)
      DOUBLE PRECISION cctt02(100,2),cctt04(100,2),cctt05(100,2)
      DOUBLE PRECISION cctt06(100,2),cctt08(100,2),cctt1(100,2)
      DOUBLE PRECISION gggg02(100,2),gggg04(100,2),gggg05(100,2)
      DOUBLE PRECISION gggg06(100,2),gggg08(100,2),gggg1(100,2)
      DOUBLE PRECISION ttgg02(100,2),ttgg04(100,2),ttgg05(100,2)
      DOUBLE PRECISION ttgg06(100,2),ttgg08(100,2),ttgg1(100,2)
      DOUBLE PRECISION tttt02(100,2),tttt04(100,2),tttt05(100,2)
      DOUBLE PRECISION tttt06(100,2),tttt08(100,2),tttt1(100,2)
      DOUBLE PRECISION stblsn(100,2),stnc(100,2),sbnb(100,2)
      DOUBLE PRECISION glsq(100,2)
      DOUBLE PRECISION XIN(6),HM1(6),HM2(6),HM3(6),HM4(6),HM5(6)

      DATA XIN/.2d0,.25d0,.4d0,.5d0,.7d0,1d0/
      DATA HM1/70d0,94d0,99.7d0,102.3d0,105.2d0,108.4d0/
      DATA HM2/70d0,93d0,98.6d0,101.6d0,105d0,108d0/
      DATA HM3/70d0,92.3d0,98.3d0,100.3d0,104.6d0,107.6d0/
      DATA HM4/70d0,90.8d0,98d0,100.1d0,104d0,107d0/
      DATA HM5/70d0,86.7d0,97d0,99d0,103.5d0,106d0/

      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/HISPEC/MHC2,XC,MH0,XH,MA2
      COMMON/HIWIDTH/WIDTH,HCWIDTH
      COMMON/HNSMBR/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG
      COMMON/HNHIBR/BRHHH,BRHCHC,BRHAZ,BRHCW,BRHIGGS
      COMMON/HNSUSYBR/BRNEU,BRCHA,BRHSQ,BRHSL,BRSUSY
      COMMON/LEP/GZMAX,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1,SIGNEU,
     .      MSLMIN,MSTMIN,MSQMIN,MGLMIN,
     .      hZind,hZbb,hZll,hZinv,hZjj,hZgg,
     .      hA4b,hA4tau,hA2b2tau,hA2tau2b,
     .      AAA6b,AAA6tau,AAZ4b,AAZ4tau,AAZ2b2tau,
     .      cccc02,cccc04,cccc05,cccc06,cccc08,cccc1,
     .      ccgg02,ccgg04,ccgg05,ccgg06,ccgg08,ccgg1,
     .      cctt02,cctt04,cctt05,cctt06,cctt08,cctt1,
     .      gggg02,gggg04,gggg05,gggg06,gggg08,gggg1,
     .      ttgg02,ttgg04,ttgg05,ttgg06,ttgg08,ttgg1,
     .      tttt02,tttt04,tttt05,tttt06,tttt08,tttt1,
     .      stblsn,stnc,sbnb,glsq,
     .      NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg,
     .      NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b,
     .      NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau,
     .      Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1,
     .      Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1,
     .      Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1,
     .      Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1,
     .      Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1,
     .      Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1,
     .      Nstblsn,Nstnc,Nsbnb,Nglsq

      LAMBDA(X,Y,Z)= DSQRT(X**2+Y**2+Z**2-2d0*X*Y-2d0*X*Z-2d0*Y*Z)

      PI=4d0*DATAN(1d0)
      SQR2=DSQRT(2d0)
      SQRS=209d0
      S=SQRS**2
      CONV=.3894D9

*   Test on charged Higgs

      CMASS=dsqrt(MHC2)
      PROB(3)=DDIM(1d0,CMASS/MCMIN)

*   Tests on neutral Higgses

c      I- Associated Z / CP-even Higgs production (Higgsstrahlung)
c          ee -> Z* -> hZ

      DO I=1,5

      IF(I.GT.1)THEN
       IF(DABS(dsqrt(MH0(I))-dsqrt(MH0(I-1))).LT.3d0)THEN
        MH=dsqrt(MH0(I-1))
        R=R+(XH(I,1)*vu+XH(I,2)*vd)**2/(vu**2+vd**2)
       ENDIF
      ELSE
       MH=dsqrt(MH0(I))
       R=(XH(I,1)*vu+XH(I,2)*vd)**2/(vu**2+vd**2)
      ENDIF

      IF(MH+MZ.LT.SQRS)THEN

*  ee -> hZ flavor independent

       Rmax=1d0
       J=1
       DOWHILE(hZind(J,1).LE.MH .AND. J.LT.NhZind)
        J=J+1
       ENDDO
       IF(J.GT.1 .AND. MH.LE.hZind(NhZind,1))
     .    Rmax=hZind(J-1,2)+(MH-hZind(J-1,1))/(hZind(J,1)
     .       -hZind(J-1,1))*(hZind(J,2)-hZind(J-1,2))
       PROB(4)=PROB(4)+DDIM(R/Rmax,1d0)

*  ee -> hZ, h -> bb

       Rmax=1d0
       J=1
       DOWHILE(hZbb(J,1).LE.MH .AND. J.LT.NhZbb)
        J=J+1
       ENDDO
       IF(J.GT.1 .AND. MH.LE.hZbb(NhZbb,1))
     .    Rmax=hZbb(J-1,2)+(MH-hZbb(J-1,1))/(hZbb(J,1)-hZbb(J-1,1))
     .       *(hZbb(J,2)-hZbb(J-1,2))
       PROB(5)=PROB(5)+DDIM(R*BRBB(I)/Rmax,1d0)

*  ee -> hZ, h -> tautau

       Rmax=1d0
       J=1
       DOWHILE(hZll(J,1).LE.MH .AND. J.LT.NhZll)
        J=J+1
       ENDDO
       IF(J.GT.1 .AND. MH.LE.hZll(NhZll,1))
     .    Rmax=hZll(J-1,2)+(MH-hZll(J-1,1))/(hZll(J,1)-hZll(J-1,1))
     .       *(hZll(J,2)-hZll(J-1,2))
       PROB(6)=PROB(6)+DDIM(R*BRLL(I)/Rmax,1d0)

*  ee -> hZ, h -> invisible

       Rmax=1d0
       J=1
       DOWHILE(hZinv(J,1).LE.MH .AND. J.LT.NhZinv)
        J=J+1
       ENDDO
       IF(J.GT.1 .AND. MH.LE.hZinv(NhZinv,1))
     .    Rmax=hZinv(J-1,2)+(MH-hZinv(J-1,1))/(hZinv(J,1)-hZinv(J-1,1))
     .       *(hZinv(J,2)-hZinv(J-1,2))
       BRINV=BRNEU(I,1,1)
       IF(I.GT.1)
     .  BRINV=BRINV+BRHHH(I,1)*BRNEU(1,1,1)**2
       IF(I.GT.2)
     .  BRINV=BRINV+BRHHH(I,2)*BRNEU(1,1,1)*BRNEU(2,1,1)
     .             +BRHHH(I,3)*BRNEU(2,1,1)**2
       IF(I.GT.3)
     .  BRINV=BRINV+BRHHH(I,4)*BRNEU(1,1,1)*BRNEU(3,1,1)
     .             +BRHHH(I,5)*BRNEU(2,1,1)*BRNEU(3,1,1)
     .             +BRHHH(I,6)*BRNEU(3,1,1)**2
       IF(I.GT.4)
     .  BRINV=BRINV+BRHHH(I,7)*BRNEU(1,1,1)*BRNEU(4,1,1)
     .             +BRHHH(I,8)*BRNEU(2,1,1)*BRNEU(4,1,1)
     .             +BRHHH(I,9)*BRNEU(3,1,1)*BRNEU(4,1,1)
     .             +BRHHH(I,10)*BRNEU(4,1,1)**2
       PROB(7)=PROB(7)+DDIM(R*BRINV/Rmax,1d0)

*  ee -> hZ, h -> 2jets

       Rmax=1d0
       J=1
       DOWHILE(hZjj(J,1).LE.MH .AND. J.LT.NhZjj)
        J=J+1
       ENDDO
       IF(J.GT.1 .AND. MH.LE.hZjj(NhZjj,1))
     .    Rmax=hZjj(J-1,2)+(MH-hZjj(J-1,1))/(hZjj(J,1)-hZjj(J-1,1))
     .       *(hZjj(J,2)-hZjj(J-1,2))
       PROB(8)=PROB(8)
     .       +DDIM(R*(BRJJ(I)+BRCC(I)+BRBB(I))/Rmax,1d0)

*  ee -> hZ, h -> 2photons

       Rmax=1d0
       J=1
       DOWHILE(hZgg(J,1).LE.MH .AND. J.LT.NhZgg)
        J=J+1
       ENDDO
       IF(J.GT.1 .AND. MH.LE.hZgg(NhZgg,1))
     .    Rmax=hZgg(J-1,2)+(MH-hZgg(J-1,1))/(hZgg(J,1)-hZgg(J-1,1))
     .       *(hZgg(J,2)-hZgg(J-1,2))
       PROB(9)=PROB(9)+DDIM(R*BRGG(I)/Rmax,1d0)

*  ee -> hZ, h -> AA -> 4bs

       MA=dsqrt(MH0(1))
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ4b
         IF(AINT(MH).EQ.AAZ4b(K,1) .AND. AINT(MA).EQ.AAZ4b(K,2))THEN
          Rmax=AAZ4b(K,3)
          GOTO 1
         ENDIF
        ENDDO
 1      PROB(10)=PROB(10)+DDIM(R*BRHHH(I,1)*BRBB(1)**2/Rmax,1d0)
       ENDIF

       MA=dsqrt(MH0(2))
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ4b
         IF(AINT(MH).EQ.AAZ4b(K,1) .AND. AINT(MA).EQ.AAZ4b(K,2))THEN
          Rmax=AAZ4b(K,3)
          GOTO 2
         ENDIF
        ENDDO
 2      PROB(10)=PROB(10)+DDIM(R*BRHHH(I,3)*BRBB(2)**2/Rmax,1d0)
       ENDIF

*  ee -> hZ, h -> AA -> 4taus

       MA=dsqrt(MH0(1))
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ4tau
         IF(AINT(MH).EQ.AAZ4tau(K,1) .AND.
     .     AINT(MA).EQ.AAZ4tau(K,2))THEN
          Rmax=AAZ4tau(K,3)
          GOTO 11
         ENDIF
        ENDDO
* Apply only for MA < 9.4GeV where A <-> eta_b mixing is absent:
 11     IF(MA.LE.9.4d0) 
     . PROB(11)=PROB(11)+DDIM(R*BRHHH(I,1)*BRLL(1)**2/Rmax,1d0)
       ENDIF

       MA=dsqrt(MH0(2))
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ4tau
         IF(AINT(MH).EQ.AAZ4tau(K,1) .AND.
     .     AINT(MA).EQ.AAZ4tau(K,2))THEN
          Rmax=AAZ4tau(K,3)
          GOTO 12
         ENDIF
        ENDDO
 12     PROB(11)=PROB(11)+DDIM(R*BRHHH(I,3)*BRLL(2)**2/Rmax,1d0)
       ENDIF

*  ee -> hZ, h -> AA -> 4taus (new ALEPH analysis)

       MA=dsqrt(MH0(1))
       RMAX=1d0
       IF(MH.GE.70d0)THEN
        IF(MA.GE.4d0.AND.MA.LT.6d0)THEN
             DMA=(6d0-MA)/2d0
             K=0
 991     K=K+1
             HMIN=HM2(K)+(HM1(K)-HM2(K))*DMA
             HMAX=HM2(K+1)+(HM1(K+1)-HM2(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 991
         ENDIF
        ELSEIF(MA.GE.6d0.AND.MA.LT.8d0)THEN
         DMA=(8d0-MA)/2d0
              K=0
 992     K=K+1
         HMIN=HM3(K)+(HM2(K)-HM3(K))*DMA
             HMAX=HM3(K+1)+(HM2(K+1)-HM3(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 992
         ENDIF
        ELSEIF(MA.GE.8d0.AND.MA.LT.10d0)THEN
              DMA=(10d0-MA)/2d0
              K=0
 993     K=K+1
             HMIN=HM4(K)+(HM3(K)-HM4(K))*DMA
             HMAX=HM4(K+1)+(HM3(K+1)-HM4(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 993
         ENDIF
        ELSEIF(MA.GE.10d0.AND.MA.LE.12d0)THEN
              DMA=(12d0-MA)/2d0
         K=0
 994     K=K+1
             HMIN=HM5(K)+(HM4(K)-HM5(K))*DMA
             HMAX=HM5(K+1)+(HM4(K+1)-HM5(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 994
         ENDIF
        ENDIF
       ENDIF
* Apply only for MA < 9.4GeV where A <-> eta_b mixing is absent:
       IF(MA.LE.9.4d0) 
     .   PROB(41)=PROB(41)+DDIM(R*BRHHH(I,1)*BRLL(1)**2/Rmax,1d0)


       MA=dsqrt(MH0(2))
       RMAX=1d0
       IF(MH.GE.70d0)THEN
        IF(MA.GE.4d0.AND.MA.LT.6d0)THEN
             DMA=(6d0-MA)/2d0
             K=0
 995     K=K+1
             HMIN=HM2(K)+(HM1(K)-HM2(K))*DMA
             HMAX=HM2(K+1)+(HM1(K+1)-HM2(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 995
         ENDIF
        ELSEIF(MA.GE.6d0.AND.MA.LT.8d0)THEN
             DMA=(8d0-MA)/2d0
             K=0
 996     K=K+1
             HMIN=HM3(K)+(HM2(K)-HM3(K))*DMA
             HMAX=HM3(K+1)+(HM2(K+1)-HM3(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 996
         ENDIF
        ELSEIF(MA.GE.8d0.AND.MA.LT.10d0)THEN
             DMA=(10d0-MA)/2d0
             K=0
 997     K=K+1
             HMIN=HM4(K)+(HM3(K)-HM4(K))*DMA
             HMAX=HM4(K+1)+(HM3(K+1)-HM4(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 997
         ENDIF
        ELSEIF(MA.GE.10d0.AND.MA.LE.12d0)THEN
             DMA=(12d0-MA)/2d0
             K=0
 998     K=K+1
             HMIN=HM5(K)+(HM4(K)-HM5(K))*DMA
             HMAX=HM5(K+1)+(HM4(K+1)-HM5(K+1))*DMA
             IF(MH.LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(MH-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 998
         ENDIF
        ENDIF
       ENDIF
       PROB(41)=PROB(41)+DDIM(R*BRHHH(I,3)*BRLL(2)**2/Rmax,1d0)

*  ee -> hZ, h -> AA -> 2bs 2taus

       MA=dsqrt(MH0(1))
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ2b2tau
         IF(AINT(MH).EQ.AAZ2b2tau(K,1) .AND.
     .     AINT(MA).EQ.AAZ2b2tau(K,2))THEN
          Rmax=AAZ2b2tau(K,3)
          GOTO 21
         ENDIF
        ENDDO
* Apply only for MA < 9.4GeV or MA > 10.5GeV without A <-> eta_b mixing:
 21    IF(MA.LE.9.4d0.OR.MA.GE.10.5d0) 
     .   PROB(12)=PROB(12)
     .   +DDIM(R*BRHHH(I,1)*2d0*BRLL(1)*BRBB(1)/Rmax,1d0)
       ENDIF

       MA=dsqrt(MH0(2))
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ2b2tau
         IF(AINT(MH).EQ.AAZ2b2tau(K,1) .AND.
     .     AINT(MA).EQ.AAZ2b2tau(K,2))THEN
          Rmax=AAZ2b2tau(K,3)
          GOTO 22
         ENDIF
        ENDDO
 22     PROB(12)=PROB(12)
     .   +DDIM(R*BRHHH(I,3)*2d0*BRBB(2)*BRLL(2)/Rmax,1d0)
       ENDIF


*  ee -> hZ -> AAZ -> Z + light pairs

       MA=dsqrt(MH0(1))
       IF(MH.LT.2d0*MA .OR. MH.LT.40d0 .OR. MH.GT.90d0
     .   .OR. MA.LT.2d0 .OR. MA.GT.12d0)GOTO 752

*      AA -> cccc

       ceff=R*BRHHH(I,1)*BRCC(1)**2
       IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 102
       IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 202
       IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 302
       IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 402
       IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 502
       IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 602
       GOTO 702

 102   Rmax=.2d0
       J=1
       DOWHILE(cccc02(J,1).LE.MH .AND. J.LT.Ncccc02)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.cccc02(Ncccc02,1))
     .  MAtest=cccc02(J-1,2)+(MH-cccc02(J-1,1))/
     .   (cccc02(J,1)-cccc02(J-1,1))
     .   *(cccc02(J,2)-cccc02(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 702

 202   Rmax=.4d0
       J=1
       DOWHILE(cccc04(J,1).LE.MH .AND. J.LT.Ncccc04)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.cccc04(Ncccc04,1))
     .  MAtest=cccc04(J-1,2)+(MH-cccc04(J-1,1))/
     .   (cccc04(J,1)-cccc04(J-1,1))
     .   *(cccc04(J,2)-cccc04(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 702

 302   Rmax=.5d0
       J=1
       DOWHILE(cccc05(J,1).LE.MH .AND. J.LT.Ncccc05)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.cccc05(Ncccc05,1))
     .  MAtest=cccc05(J-1,2)+(MH-cccc05(J-1,1))/
     .   (cccc05(J,1)-cccc05(J-1,1))
     .   *(cccc05(J,2)-cccc05(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 702

 402   Rmax=.6d0
       J=1
       DOWHILE(cccc06(J,1).LE.MH .AND. J.LT.Ncccc06)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.cccc06(Ncccc06,1))
     .  MAtest=cccc06(J-1,2)+(MH-cccc06(J-1,1))/
     .   (cccc06(J,1)-cccc06(J-1,1))
     .   *(cccc06(J,2)-cccc06(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 702

 502   Rmax=.8d0
       J=1
       DOWHILE(cccc08(J,1).LE.MH .AND. J.LT.Ncccc08)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.cccc08(Ncccc08,1))
     .  MAtest=cccc08(J-1,2)+(MH-cccc08(J-1,1))/
     .   (cccc08(J,1)-cccc08(J-1,1))
     .   *(cccc08(J,2)-cccc08(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 702

 602   Rmax=1d0
       J=1
       DOWHILE(cccc1(J,1).LE.MH .AND. J.LT.Ncccc1)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.cccc1(Ncccc1,1))
     .  MAtest=cccc1(J-1,2)+(MH-cccc1(J-1,1))/
     .   (cccc1(J,1)-cccc1(J-1,1))
     .   *(cccc1(J,2)-cccc1(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)

 702   CONTINUE

*      AA -> ccjj

       ceff=R*BRHHH(I,1)*BRCC(1)*BRJJ(1)*2d0
       IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 112
       IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 212
       IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 312
       IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 412
       IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 512
       IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 612
       GOTO 712

 112   Rmax=.2d0
       J=1
       DOWHILE(ccgg02(J,1).LE.MH .AND. J.LT.Nccgg02)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.ccgg02(Nccgg02,1))
     .  MAtest=ccgg02(J-1,2)+(MH-ccgg02(J-1,1))/
     .   (ccgg02(J,1)-ccgg02(J-1,1))
     .   *(ccgg02(J,2)-ccgg02(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 712

 212   Rmax=.4d0
       J=1
       DOWHILE(ccgg04(J,1).LE.MH .AND. J.LT.Nccgg04)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.ccgg04(Nccgg04,1))
     .  MAtest=ccgg04(J-1,2)+(MH-ccgg04(J-1,1))/
     .   (ccgg04(J,1)-ccgg04(J-1,1))
     .   *(ccgg04(J,2)-ccgg04(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 712

 312   Rmax=.5d0
       J=1
       DOWHILE(ccgg05(J,1).LE.MH .AND. J.LT.Nccgg05)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.ccgg05(Nccgg05,1))
     .  MAtest=ccgg05(J-1,2)+(MH-ccgg05(J-1,1))/
     .   (ccgg05(J,1)-ccgg05(J-1,1))
     .   *(ccgg05(J,2)-ccgg05(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 712

 412   Rmax=.6d0
       J=1
       DOWHILE(ccgg06(J,1).LE.MH .AND. J.LT.Nccgg06)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.ccgg06(Nccgg06,1))
     .    MAtest=ccgg06(J-1,2)+(MH-ccgg06(J-1,1))/
     .   (ccgg06(J,1)-ccgg06(J-1,1))
     .   *(ccgg06(J,2)-ccgg06(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 712

 512   Rmax=.8d0
       J=1
       DOWHILE(ccgg08(J,1).LE.MH .AND. J.LT.Nccgg08)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.ccgg08(Nccgg08,1))
     .  MAtest=ccgg08(J-1,2)+(MH-ccgg08(J-1,1))/
     .   (ccgg08(J,1)-ccgg08(J-1,1))
     .   *(ccgg08(J,2)-ccgg08(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 712

 612   Rmax=1d0
       J=1
       DOWHILE(ccgg1(J,1).LE.MH .AND. J.LT.Nccgg1)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.ccgg1(Nccgg1,1))
     .  MAtest=ccgg1(J-1,2)+(MH-ccgg1(J-1,1))/
     .   (ccgg1(J,1)-ccgg1(J-1,1))
     .   *(ccgg1(J,2)-ccgg1(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
      
 712   CONTINUE

*      AA -> cctautau

       ceff=R*BRHHH(I,1)*BRCC(1)*BRLL(1)*2d0
       IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 122
       IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 222
       IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 322
       IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 422
       IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 522
       IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 622
       GOTO 722

 122   Rmax=.2d0
       J=1
       DOWHILE(cctt02(J,1).LE.MH .AND. J.LT.Ncctt02)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.cctt02(Ncctt02,1))
     .  MAtest=cctt02(J-1,2)+(MH-cctt02(J-1,1))/
     .   (cctt02(J,1)-cctt02(J-1,1))
     .   *(cctt02(J,2)-cctt02(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 722

 222   Rmax=.4d0
       J=1
       DOWHILE(cctt04(J,1).LE.MH .AND. J.LT.Ncctt04)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.cctt04(Ncctt04,1))
     .  MAtest=cctt04(J-1,2)+(MH-cctt04(J-1,1))/
     .   (cctt04(J,1)-cctt04(J-1,1))
     .   *(cctt04(J,2)-cctt04(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 722

 322   Rmax=.5d0
       J=1
       DOWHILE(cctt05(J,1).LE.MH .AND. J.LT.Ncctt05)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.cctt05(Ncctt05,1))
     .  MAtest=cctt05(J-1,2)+(MH-cctt05(J-1,1))/
     .   (cctt05(J,1)-cctt05(J-1,1))
     .   *(cctt05(J,2)-cctt05(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 722

 422   Rmax=.6d0
       J=1
       DOWHILE(cctt06(J,1).LE.MH .AND. J.LT.Ncctt06)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.cctt06(Ncctt06,1))
     .  MAtest=cctt06(J-1,2)+(MH-cctt06(J-1,1))/
     .   (cctt06(J,1)-cctt06(J-1,1))
     .   *(cctt06(J,2)-cctt06(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 722

 522   Rmax=.8d0
       J=1
       DOWHILE(cctt08(J,1).LE.MH .AND. J.LT.Ncctt08)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.cctt08(Ncctt08,1))
     .  MAtest=cctt08(J-1,2)+(MH-cctt08(J-1,1))/
     .   (cctt08(J,1)-cctt08(J-1,1))
     .   *(cctt08(J,2)-cctt08(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 722

 622   Rmax=1d0
       J=1
       DOWHILE(cctt1(J,1).LE.MH .AND. J.LT.Ncctt1)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.cctt1(Ncctt1,1))
     .  MAtest=cctt1(J-1,2)+(MH-cctt1(J-1,1))/
     .   (cctt1(J,1)-cctt1(J-1,1))
     .   *(cctt1(J,2)-cctt1(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
      
 722   CONTINUE

*      AA -> jjjj

       ceff=R*BRHHH(I,1)*BRJJ(1)**2
       IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 132
       IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 232
       IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 332
       IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 432
       IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 532
       IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 632
       GOTO 732

 132   Rmax=.2d0
       J=1
       DOWHILE(gggg02(J,1).LE.MH .AND. J.LT.Ngggg02)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.gggg02(Ngggg02,1))
     .  MAtest=gggg02(J-1,2)+(MH-gggg02(J-1,1))/
     .   (gggg02(J,1)-gggg02(J-1,1))
     .   *(gggg02(J,2)-gggg02(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 732

 232   Rmax=.4d0
       J=1
       DOWHILE(gggg04(J,1).LE.MH .AND. J.LT.Ngggg04)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.gggg04(Ngggg04,1))
     .  MAtest=gggg04(J-1,2)+(MH-gggg04(J-1,1))/
     .   (gggg04(J,1)-gggg04(J-1,1))
     .   *(gggg04(J,2)-gggg04(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 732

 332   Rmax=.5d0
       J=1
       DOWHILE(gggg05(J,1).LE.MH .AND. J.LT.Ngggg05)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.gggg05(Ngggg05,1))
     .  MAtest=gggg05(J-1,2)+(MH-gggg05(J-1,1))/
     .   (gggg05(J,1)-gggg05(J-1,1))
     .   *(gggg05(J,2)-gggg05(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 732

 432   Rmax=.6d0
       J=1
       DOWHILE(gggg06(J,1).LE.MH .AND. J.LT.Ngggg06)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.gggg06(Ngggg06,1))
     .  MAtest=gggg06(J-1,2)+(MH-gggg06(J-1,1))/
     .   (gggg06(J,1)-gggg06(J-1,1))
     .   *(gggg06(J,2)-gggg06(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 732

 532   Rmax=.8d0
       J=1
       DOWHILE(gggg08(J,1).LE.MH .AND. J.LT.Ngggg08)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.gggg08(Ngggg08,1))
     .  MAtest=gggg08(J-1,2)+(MH-gggg08(J-1,1))/
     .   (gggg08(J,1)-gggg08(J-1,1))
     .   *(gggg08(J,2)-gggg08(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 732

 632   Rmax=1d0
       J=1
       DOWHILE(gggg1(J,1).LE.MH .AND. J.LT.Ngggg1)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.gggg1(Ngggg1,1))
     .  MAtest=gggg1(J-1,2)+(MH-gggg1(J-1,1))/
     .   (gggg1(J,1)-gggg1(J-1,1))
     .   *(gggg1(J,2)-gggg1(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
      
 732   CONTINUE

*      AA -> tautaujj

       ceff=R*BRHHH(I,1)*BRLL(1)*BRJJ(1)*2d0
       IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 142
       IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 242
       IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 342
       IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 442
       IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 542
       IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 642
       GOTO 742

 142   Rmax=.2d0
       J=1
       DOWHILE(ttgg02(J,1).LE.MH .AND. J.LT.Nttgg02)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.ttgg02(Nttgg02,1))
     .  MAtest=ttgg02(J-1,2)+(MH-ttgg02(J-1,1))/
     .   (ttgg02(J,1)-ttgg02(J-1,1))
     .   *(ttgg02(J,2)-ttgg02(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 742

 242   Rmax=.4d0
       J=1
       DOWHILE(ttgg04(J,1).LE.MH .AND. J.LT.Nttgg04)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.ttgg04(Nttgg04,1))
     .  MAtest=ttgg04(J-1,2)+(MH-ttgg04(J-1,1))/
     .   (ttgg04(J,1)-ttgg04(J-1,1))
     .   *(ttgg04(J,2)-ttgg04(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 742

 342   Rmax=.5d0
       J=1
       DOWHILE(ttgg05(J,1).LE.MH .AND. J.LT.Nttgg05)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.ttgg05(Nttgg05,1))
     .  MAtest=ttgg05(J-1,2)+(MH-ttgg05(J-1,1))/
     .   (ttgg05(J,1)-ttgg05(J-1,1))
     .   *(ttgg05(J,2)-ttgg05(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 742

 442   Rmax=.6d0
       J=1
       DOWHILE(ttgg06(J,1).LE.MH .AND. J.LT.Nttgg06)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.ttgg06(Nttgg06,1))
     .  MAtest=ttgg06(J-1,2)+(MH-ttgg06(J-1,1))/
     .   (ttgg06(J,1)-ttgg06(J-1,1))
     .   *(ttgg06(J,2)-ttgg06(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 742

 542   Rmax=.8d0
       J=1
       DOWHILE(ttgg08(J,1).LE.MH .AND. J.LT.Nttgg08)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.ttgg08(Nttgg08,1))
     .  MAtest=ttgg08(J-1,2)+(MH-ttgg08(J-1,1))/
     .   (ttgg08(J,1)-ttgg08(J-1,1))
     .   *(ttgg08(J,2)-ttgg08(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 742

 642   Rmax=1d0
       J=1
       DOWHILE(ttgg1(J,1).LE.MH .AND. J.LT.Nttgg1)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.ttgg1(Nttgg1,1))
     .  MAtest=ttgg1(J-1,2)+(MH-ttgg1(J-1,1))/
     .   (ttgg1(J,1)-ttgg1(J-1,1))
     .   *(ttgg1(J,2)-ttgg1(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
      
 742   CONTINUE

*      AA -> tautautautau

       ceff=R*BRHHH(I,1)*BRLL(1)**2
       IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 152
       IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 252
       IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 352
       IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 452
       IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 552
       IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 652
       GOTO 752

 152   Rmax=.2d0
       J=1
       DOWHILE(tttt02(J,1).LE.MH .AND. J.LT.Ntttt02)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.tttt02(Ntttt02,1))
     .  MAtest=tttt02(J-1,2)+(MH-tttt02(J-1,1))/
     .   (tttt02(J,1)-tttt02(J-1,1))
     .   *(tttt02(J,2)-tttt02(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 752

 252   Rmax=.4d0
       J=1
       DOWHILE(tttt04(J,1).LE.MH .AND. J.LT.Ntttt04)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.tttt04(Ntttt04,1))
     .  MAtest=tttt04(J-1,2)+(MH-tttt04(J-1,1))/
     .   (tttt04(J,1)-tttt04(J-1,1))
     .   *(tttt04(J,2)-tttt04(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 752

 352   Rmax=.5d0
       J=1
       DOWHILE(tttt05(J,1).LE.MH .AND. J.LT.Ntttt05)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.tttt05(Ntttt05,1))
     .  MAtest=tttt05(J-1,2)+(MH-tttt05(J-1,1))/
     .   (tttt05(J,1)-tttt05(J-1,1))
     .   *(tttt05(J,2)-tttt05(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 752

 452   Rmax=.6d0
       J=1
       DOWHILE(tttt06(J,1).LE.MH .AND. J.LT.Ntttt06)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.tttt06(Ntttt06,1))
     .  MAtest=tttt06(J-1,2)+(MH-tttt06(J-1,1))/
     .   (tttt06(J,1)-tttt06(J-1,1))
     .   *(tttt06(J,2)-tttt06(J-1,2))
       PROB(19)=PROB(19)+10d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 752

 552   Rmax=.8d0
       J=1
       DOWHILE(tttt08(J,1).LE.MH .AND. J.LT.Ntttt08)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.tttt08(Ntttt08,1))
     .  MAtest=tttt08(J-1,2)+(MH-tttt08(J-1,1))/
     .   (tttt08(J,1)-tttt08(J-1,1))
     .   *(tttt08(J,2)-tttt08(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
       GOTO 752

 652   Rmax=1d0
       J=1
       DOWHILE(tttt1(J,1).LE.MH .AND. J.LT.Ntttt1)
        J=J+1
       ENDDO
       MAtest=0d0
       IF(J.GT.1 .AND. MH.LE.tttt1(Ntttt1,1))
     .  MAtest=tttt1(J-1,2)+(MH-tttt1(J-1,1))/
     .   (tttt1(J,1)-tttt1(J-1,1))
     .   *(tttt1(J,2)-tttt1(J-1,2))
       PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MA,1d0)
      
 752   CONTINUE
      ENDIF

      ENDDO


c      II- Pair production: ee -> Z* -> h_i+h_j

      DO I=1,5
      DO J=1,I

      MH=dsqrt(MH0(I))
      MA=dsqrt(MH0(J))
      R=(XH(I,4)*(XH(J,1)*vd-XH(J,2)*vu)
     .  -XH(J,4)*(XH(I,1)*vd-XH(I,2)*vu))**2/(vu**2+vd**2)
      IF(I.eq.J)R=R/2d0

*  Z width

      IF(MH+MA.LT.MZ)THEN
       GZ=SQR2*GF*R*LAMBDA(MZ**2,MA**2,MH**2)**3/(48d0*PI*MZ**3)
        PROB(13)=PROB(13)+DDIM(GZ/GZMAX,1d0)
      ENDIF

*  ee -> hA -> 4bs

      IF(MH+MA.LT.SQRS)THEN

       M1=MIN(MH,MA)
       M2=MAX(MH,MA)

       Rmax=1d0
       DO K=1,NhA4b
        IF(AINT(M2).EQ.hA4b(K,1) .AND. AINT(M1).EQ.hA4b(K,2))THEN
         Rmax=hA4b(K,3)
         GOTO 3
        ENDIF
       ENDDO
 3     PROB(14)=PROB(14)+DDIM(R*BRBB(I)*BRBB(J)/Rmax,1d0)

*  ee -> hA -> 4taus

       Rmax=1d0
       DO K=1,NhA4tau
        IF(AINT(M2).EQ.hA4tau(K,1) .AND. AINT(M1).EQ.hA4tau(K,2))THEN
         Rmax=hA4tau(K,3)
         GOTO 4
        ENDIF
       ENDDO
 4     PROB(15)=PROB(15)+DDIM(R*BRLL(I)*BRLL(J)/Rmax,1d0)

*  ee -> hA -> 2b 2taus

       Rmax=1d0
       IF(MA.GT.MH)THEN
        DO K=1,NhA2b2tau
         IF(AINT(M2).EQ.hA2b2tau(K,1) .AND.
     .     AINT(M1).EQ.hA2b2tau(K,2))THEN
          Rmax=hA2b2tau(K,3)
          GOTO 5
         ENDIF
        ENDDO
 5      PROB(16)=PROB(16)+DDIM(R*BRLL(I)*BRBB(J)/Rmax,1d0)
       ENDIF

       Rmax=1d0
       IF(MH.GT.MA)THEN
        DO K=1,NhA2b2tau
         IF(AINT(M2).EQ.hA2b2tau(K,1) .AND.
     .     AINT(M1).EQ.hA2b2tau(K,2))THEN
          Rmax=hA2b2tau(K,3)
          GOTO 6
         ENDIF
        ENDDO
 6      PROB(16)=PROB(16)+DDIM(R*BRBB(I)*BRLL(J)/Rmax,1d0)
       ENDIF

       Rmax=1d0
       IF(MA.GT.MH)THEN
        DO K=1,NhA2tau2b
         IF(AINT(M2).EQ.hA2tau2b(K,1) .AND.
     .     AINT(M1).EQ.hA2tau2b(K,2))THEN
          Rmax=hA2tau2b(K,3)
          GOTO 7
         ENDIF
        ENDDO
 7      PROB(16)=PROB(16)+DDIM(R*BRBB(I)*BRLL(J)/Rmax,1d0)
       ENDIF

       Rmax=1d0
       IF(MH.GT.MA)THEN
        DO K=1,NhA2tau2b
         IF(AINT(M2).EQ.hA2tau2b(K,1) .AND.
     .     AINT(M1).EQ.hA2tau2b(K,2))THEN
          Rmax=hA2tau2b(K,3)
          GOTO 8
         ENDIF
        ENDDO
 8      PROB(16)=PROB(16)+DDIM(R*BRLL(I)*BRBB(J)/Rmax,1d0)
       ENDIF

*  ee -> hA -> AAA -> 6bs

       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAA6b
         IF(AINT(MH).EQ.AAA6b(K,1) .AND. AINT(MA).EQ.AAA6b(K,2))THEN
          Rmax=AAA6b(K,3)
          GOTO 9
         ENDIF
        ENDDO
 9      PROB(17)=PROB(17)+DDIM(R*BRHHH(I,J*(J+1)/2)*BRBB(J)**3/Rmax,
     .   1d0)

*  ee -> hA -> AAA -> 6taus

        Rmax=1d0
        DO K=1,NAAA6tau
         IF(AINT(MH).EQ.AAA6tau(K,1) .AND.
     .     AINT(MA).EQ.AAA6tau(K,2))THEN
          Rmax=AAA6tau(K,3)
          GOTO 10
         ENDIF
        ENDDO
 10     PROB(18)=PROB(18)+DDIM(R*BRHHH(I,J*(J+1)/2)*BRLL(J)**3/Rmax,
     .   1d0)
      ENDIF

      ENDIF

      ENDDO
      ENDDO

      END


      SUBROUTINE TEVATRON_CHIGGS_CPV(PROB)

***********************************************************************
*   Subroutine to check TeVatron constraints on a charged Higgs:
*      PROB(42) =/= 0  excluded by top -> b H+, H+ -> c s (CDF, D0)
*      PROB(43) =/= 0  excluded by top -> b H+, H+ -> tau nu_tau (D0)
*      PROB(44) =/= 0  excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus (CDF)
***********************************************************************

      IMPLICIT NONE

      INTEGER I,J

      DOUBLE PRECISION PROB(*),MA,CMASS,brtoplim
      DOUBLE PRECISION minf(9),msup(9),binf(9),bsup(9)
      DOUBLE PRECISION mh4m(4),mh4p(4),br4m(4),br4p(4),brtau(8)
      DOUBLE PRECISION mh7(5),br7(5),br9(5)
      DOUBLE PRECISION brtopcs,brtoptau,brtopa1

      DOUBLE PRECISION MHC2,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,HCBRBT
      DOUBLE PRECISION HCBRWH(5),HCBRWHT
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),
     . BRCC(5),BRBB(5),BRTT(5),BRWW(5),BRZZ(5),BRGG(5),BRZG(5)

      DATA minf/60d0,70d0,80d0,100d0,110d0,120d0,130d0,140d0,150d0/
      DATA msup/70d0,80d0,100d0,110d0,120d0,130d0,140d0,150d0,155d0/
      DATA binf/.09d0,0d0,.21d0,.21d0,.15d0,.12d0,.08d0,.1d0,.20d0/
      DATA bsup/.12d0,0d0,.21d0,.15d0,.12d0,.08d0,.1d0,.13d0,.19d0/
      DATA mh4m/85d0,90d0,100d0,120d0/
      DATA mh4p/90d0,100d0,120d0,140d0/
      DATA br4m/.5d0,.33d0,.27d0,.35d0/
      DATA br4p/.33d0,.27d0,.35d0,.52d0/
      DATA mh7/90d0,100d0,120d0,140d0,160d0/
      DATA br7/.13d0,.1d0,.13d0,.2d0,.3d0/
      DATA br9/.11d0,.08d0,.09d0,.14d0,.21d0/
      DATA brtau/.16d0,.15d0,.16,.17d0,.175,.18d0,.19d0,.18d0/

      COMMON/HISPEC/MHC2,XC,MH0,XH,MA2
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/HCSMBR/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,HCBRBT
      COMMON/HCHIBR/HCBRWH,HCBRWHT
      COMMON/HNSMBR/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG

***********************************************************************

      CMASS=dsqrt(MHC2)

* top -> H+ b, H+ -> c s (from CDF, 0907.1269, and D0, 0908.1811)
       
       brtopcs=brtopbh*HCBRSC
       PROB(42)=0d0
       
       DO I=1,9
        IF((minf(i).LE.CMASS).AND.(CMASS.LE.msup(i))) THEN
        brtoplim=binf(i)+
     .      (CMASS-minf(i))*(bsup(i)-binf(i))/(msup(i)-minf(i))
         PROB(42)=PROB(42)+DDIM(brtopcs/brtoplim,1d0)
        ENDIF
       ENDDO

* top -> H+ b, H+ -> tau nu_tau (from D0, 0908.1811)

       brtoptau=brtopbh*HCBRL
       PROB(43)=0d0
    
       DO I=1,7
        IF((minf(i+2).LE.CMASS).AND.(CMASS.LE.msup(i+2))) THEN
        brtoplim=brtau(i)+
     .    (CMASS-minf(i+2))*(brtau(i+1)-brtau(i))/(msup(i+2)-minf(i+2))
         PROB(43)=PROB(43)+DDIM(brtoptau/brtoplim,1d0)
        ENDIF
       ENDDO

* top -> H+ b, H+ -> W+ A_1, A_1 -> 2taus (from CDF Note 10104)

       PROB(44)=0d0

      DO J=1,4
       MA=dsqrt(MH0(J))
       brtopa1=brtopbh*HCBRWH(J)*BRLL(J)

      IF(MA.LE.4d0) then     
       DO I=1,4
        IF((mh4m(i).LE.CMASS).AND.(CMASS.LE.mh4p(i))) THEN
          brtoplim=br4m(i)+
     .      (CMASS-mh4m(i))*(br4p(i)-br4m(i))/(mh4p(i)-mh4m(i))
         PROB(44)=PROB(44)+DDIM(brtopa1/brtoplim,1d0)
        ENDIF
       ENDDO
      ENDIF

      IF((MA.gt.4d0).and.(MA.LE.7d0)) then     
       DO I=1,3
        IF((mh7(i).LE.CMASS).AND.(CMASS.LE.mh7(i+1))) THEN
          brtoplim=br4p(i)+(MA-4d0)*(br7(i)-br4p(i))/3d0
     .     +(br4p(i+1)-br4p(i)
     .      +(br7(i+1)-br7(i)-br4p(i+1)+br4p(i))*(MA-4d0)/3d0)
     .     *(CMASS-mh7(i))/(mh7(i+1)-mh7(i))
         PROB(44)=PROB(44)+DDIM(brtopa1/brtoplim,1d0)
        ENDIF
       ENDDO
      ENDIF

      IF((MA.gt.7d0).and.(MA.LE.9d0)) then     
       DO I=1,4
        IF((mh7(i).LE.CMASS).AND.(CMASS.LE.mh7(i+1))) THEN
          brtoplim=br7(i)+(MA-7d0)*(br9(i)-br7(i))/2d0
     .     +(br7(i+1)-br7(i)
     .      +(br9(i+1)-br9(i)-br7(i+1)+br7(i))*(MA-7d0)/2d0)
     .     *(CMASS-mh7(i))/(mh7(i+1)-mh7(i))
         PROB(44)=PROB(44)+DDIM(brtopa1/brtoplim,1d0)
        ENDIF
       ENDDO
      ENDIF

      ENDDO

      END
