      SUBROUTINE SUBEXP(PAR,PROB)

***********************************************************************
*   Subroutine to check experimental constraints
*
*   The required data files and numbers (inv. Z width, lower bounds on
*   sparticle masses) are transferred via COMMON/LEP/...
*            
*      PROB(1) =/= 0   chargino too light
*      PROB(2) =/= 0   excluded by Z -> neutralinos
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
*      PROB(20) =/= 0  excluded by stop -> b l sneutrino
*      PROB(21) =/= 0  excluded by stop -> neutralino c
*      PROB(22) =/= 0  excluded by sbottom -> neutralino b
*      PROB(23) =/= 0  squark/gluino too light
*      PROB(24) =/= 0  selectron/smuon too light
*      PROB(25) =/= 0  stau too light
*      PROB(41) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus (new ALEPH analysis)
*      PROB(42) =/= 0  excluded by top -> b H+, H+ -> c s (CDF, D0)
*      PROB(43) =/= 0  excluded by top -> b H+, H+ -> tau nu_tau (D0)
*      PROB(44) =/= 0  excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus (CDF)
*
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

      DOUBLE PRECISION PAR(*),PROB(*),PI,S,SQRS,SQR2,CONV
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2)
      DOUBLE PRECISION brtopcs,brtoplim,brtoptau,brtopa1
      DOUBLE PRECISION minf(9),msup(9),binf(9),bsup(9)
      DOUBLE PRECISION mh4m(4),mh4p(4),br4m(4),br4p(4),brtau(8)
      DOUBLE PRECISION mh7(5),br7(5),br9(5)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH
      DOUBLE PRECISION MH,MA,TANB,h1,h2,R,Rmax,BRINV,M1,M2
      DOUBLE PRECISION GZ,GZMAX,ALSMZ,ALEMMZ,GF,g1,g2,S2TW,g
      DOUBLE PRECISION GZINV,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1,SIGNEU
      DOUBLE PRECISION MSLMIN,MSTMIN,MSQMIN,MGLMIN,MMIN
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION LAMBDA,X,Y,Z,Q,O,DZ,E1,E2,SIG
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU,MAtest,ceff
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
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION XIN(6),HM1(6),HM2(6),HM3(6),HM4(6),HM5(6)
      DOUBLE PRECISION DMA,HMIN,HMAX
   
      DATA XIN/.2d0,.25d0,.4d0,.5d0,.7d0,1d0/
      DATA HM1/70d0,94d0,99.7d0,102.3d0,105.2d0,108.4d0/
      DATA HM2/70d0,93d0,98.6d0,101.6d0,105d0,108d0/
      DATA HM3/70d0,92.3d0,98.3d0,100.3d0,104.6d0,107.6d0/
      DATA HM4/70d0,90.8d0,98d0,100.1d0,104d0,107d0/
      DATA HM5/70d0,86.7d0,97d0,99d0,103.5d0,106d0/
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

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     . HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     . HCBRSUSY,HCWIDTH
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
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
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL

      LAMBDA(X,Y,Z)= DSQRT(X**2+Y**2+Z**2-2d0*X*Y-2d0*X*Z-2d0*Y*Z)

      PI=4d0*DATAN(1d0)
      SQR2=DSQRT(2d0)
      SQRS=209d0
      S=SQRS**2
      CONV=.3894D9
      g=(g1+g2)/2d0
      TANB=PAR(3)
      h2=1d0/DSQRT(2d0*DSQRT(2d0)*(1d0+TANB**2)*GF)
      h1=h2*TANB

* Test on stop -> b l sneutrino

      I=1
      DOWHILE(stblsn(I,1).LE.MST1 .AND. I.LT.Nstblsn)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MST1.LE.stblsn(Nstblsn,1))THEN
       MMIN=stblsn(I-1,2)+(MST1-stblsn(I-1,1))
     .  /(stblsn(I,1)-stblsn(I-1,1))*(stblsn(I,2)-stblsn(I-1,2))
       PROB(20)=DDIM(1d0,MNL/MMIN)
      ENDIF

* Test on stop -> neutralino c

      I=1
      DOWHILE(stnc(I,1).LE.DABS(MNEU(1)) .AND. I.LT.Nstnc)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. DABS(MNEU(1)).LE.stnc(Nstnc,1))THEN
       MMIN=stnc(I-1,2)+(DABS(MNEU(1))-stnc(I-1,1))
     .  /(stnc(I,1)-stnc(I-1,1))*(stnc(I,2)-stnc(I-1,2))
       PROB(21)=DDIM(1d0,MST1/MMIN)
      ENDIF

* Test on sbottom -> neutralino b

      I=1
      DOWHILE(sbnb(I,1).LE.MSB1 .AND. I.LT.Nsbnb)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MSB1.LE.sbnb(Nsbnb,1))THEN
       MMIN=sbnb(I-1,2)+(MSB1-sbnb(I-1,1))/(sbnb(I,1)-sbnb(I-1,1))
     .  *(sbnb(I,2)-sbnb(I-1,2))
       PROB(22)=DDIM(1d0,DABS(MNEU(1))/MMIN)
      ENDIF

* Test on gluino/squark masses

      I=1
      DOWHILE(glsq(I,1).LE.DABS(MGL) .AND. I.LT.Nglsq)
       I=I+1
      ENDDO
      MMIN=MSQMIN
      IF(I.GT.1 .AND. DABS(MGL).LE.glsq(Nglsq,1))THEN
       MMIN=glsq(I-1,2)+(DABS(MGL)-glsq(I-1,1))
     .  /(glsq(I,1)-glsq(I-1,1))*(glsq(I,2)-glsq(I-1,2))
      ENDIF
      PROB(23)=DDIM(1d0,MIN(MUR,MUL,MDR,MDL)/MMIN)
     .       +DDIM(1d0,DABS(MGL)/MGLMIN)

* Test on slepton masses

      PROB(24)=DDIM(1d0,MIN(MLR,MLL)/MSLMIN)
      PROB(25)=DDIM(1d0,MSL1/MSTMIN)

* Test on chargino mass

      PROB(1)=DDIM(1d0,DABS(MCH(1))/MCHAMIN)

* Test on Z width into neutralinos

      IF(DABS(MNEU(1)).LT.MZ/2d0)THEN
       GZINV=MZ**3*GF/(12d0*DSQRT(2d0)*PI)
     .     *(NEU(1,3)**2-NEU(1,4)**2)**2*(1d0-4*MNEU(1)**2/MZ**2)**1.5d0
        PROB(2)=DDIM(GZINV/GZINVMAX,1d0)
      ENDIF

* Test on neutralinos

      DZ = 1d0/(S-MZ**2)**2
      DO I=1,5
       DO J=I,5
        IF(DABS(MNEU(I))+DABS(MNEU(J)).LT.SQRS)THEN
         O= (NEU(I,3)*NEU(J,3)-NEU(I,4)*NEU(J,4))
         E1= (S+MNEU(I)**2-MNEU(J)**2)/(2d0*SQRS)
         E2= (S+MNEU(J)**2-MNEU(I)**2)/(2d0*SQRS)
         Q= LAMBDA(S,MNEU(I)**2,MNEU(J)**2)/(2d0*SQRS)
         SIG= 1d0/(4d0*PI)*Q/SQRS*DZ*(E1*E2+Q**2/3d0
     .       -MNEU(I)*MNEU(J))* g**2*O**2*(.25d0-S2TW+2d0*S2TW**2)
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

* Test on charged Higgs

      PROB(3)=DDIM(1d0,CMASS/MCMIN)
      
* Higgs Strahlung

      DO I=1,3

      IF(I.GT.1.AND.DABS(SMASS(I)-SMASS(I-1)).LT.3d0)THEN
       MH=SMASS(I-1)
       R=R+(SCOMP(I,1)*h1+SCOMP(I,2)*h2)**2/(h1**2+h2**2)
      ELSE
       MH=SMASS(I)
       R=(SCOMP(I,1)*h1+SCOMP(I,2)*h2)**2/(h1**2+h2**2)
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
     .  +BRHAA(I,1)*BRNEU(4,1,1)**2
     .  +BRHAA(I,2)*BRNEU(4,1,1)*BRNEU(5,1,1)
     .  +BRHAA(I,3)*BRNEU(5,1,1)**2
       IF(I.EQ.2)
     .  BRINV=BRINV+BRHHH(1)*BRNEU(1,1,1)**2
       IF(I.EQ.3)
     .  BRINV=BRINV+BRHHH(2)*BRNEU(1,1,1)**2
     .   +BRHHH(3)*BRNEU(1,1,1)*BRNEU(2,1,1)
     .   +BRHHH(4)*BRNEU(2,1,1)**2
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

       MA=PMASS(1)
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ4b
         IF(AINT(MH).EQ.AAZ4b(K,1) .AND. AINT(MA).EQ.AAZ4b(K,2))THEN
          Rmax=AAZ4b(K,3)
          GOTO 1
         ENDIF
        ENDDO
 1      PROB(10)=PROB(10)+DDIM(R*BRHAA(I,1)*BRBB(4)**2/Rmax,1d0)
       ENDIF

       MA=SMASS(1)
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ4b
         IF(AINT(MH).EQ.AAZ4b(K,1) .AND. AINT(MA).EQ.AAZ4b(K,2))THEN
          Rmax=AAZ4b(K,3)
          GOTO 2
         ENDIF
        ENDDO
 2      PROB(10)=PROB(10)+DDIM(R*BRHHH(I-1)*BRBB(1)**2/Rmax,1d0)
       ENDIF

*  ee -> hZ, h -> AA -> 4taus

       MA=PMASS(1)
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
     . PROB(11)=PROB(11)+DDIM(R*BRHAA(I,1)*BRLL(4)**2/Rmax,1d0)
       ENDIF

       MA=SMASS(1)
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ4tau
         IF(AINT(MH).EQ.AAZ4tau(K,1) .AND.
     .     AINT(MA).EQ.AAZ4tau(K,2))THEN
          Rmax=AAZ4tau(K,3)
          GOTO 12
         ENDIF
        ENDDO
 12     PROB(11)=PROB(11)+DDIM(R*BRHHH(I-1)*BRLL(1)**2/Rmax,1d0)
       ENDIF

*  ee -> hZ, h -> AA -> 4taus (new ALEPH analysis)

       MA=PMASS(1)
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
     .   PROB(41)=PROB(41)+DDIM(R*BRHAA(I,1)*BRLL(4)**2/Rmax,1d0)

       MA=SMASS(1)
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
       IF(I.GT.1)
     . PROB(41)=PROB(41)+DDIM(R*BRHHH(I-1)*BRLL(1)**2/Rmax,1d0)

*  ee -> hZ, h -> AA -> 2bs 2taus

       MA=PMASS(1)
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
     .   +DDIM(R*BRHAA(I,1)*2d0*BRLL(4)*BRBB(4)/Rmax,1d0)
       ENDIF

       MA=SMASS(1)
       IF(MH.GT.2d0*MA)THEN
        Rmax=1d0
        DO K=1,NAAZ2b2tau
         IF(AINT(MH).EQ.AAZ2b2tau(K,1) .AND.
     .     AINT(MA).EQ.AAZ2b2tau(K,2))THEN
          Rmax=AAZ2b2tau(K,3)
          GOTO 22
         ENDIF
        ENDDO
 22     IF(I.GT.1)
     . PROB(12)=PROB(12)
     .   +DDIM(R*BRHHH(I-1)*2d0*BRBB(1)*BRLL(1)/Rmax,1d0)
       ENDIF

*  ee -> hZ -> AAZ -> Z + light pairs

       MA=PMASS(1)
       IF(MH.LT.2d0*MA .OR. MH.LT.40d0 .OR. MH.GT.90d0
     .   .OR. MA.LT.2d0 .OR. MA.GT.12d0)GOTO 752

*      AA -> cccc

       ceff=R*BRHAA(I,1)*BRCC(4)**2
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

       ceff=R*BRHAA(I,1)*BRCC(4)*BRJJ(4)*2d0
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

       ceff=R*BRHAA(I,1)*BRCC(4)*BRLL(4)*2d0
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

       ceff=R*BRHAA(I,1)*BRJJ(4)**2
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

       ceff=R*BRHAA(I,1)*BRLL(4)*BRJJ(4)*2d0
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

       ceff=R*BRHAA(I,1)*BRLL(4)**2
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

* Associated production

      DO I=1,3
      DO J=1,2
      MH=SMASS(I)
      MA=PMASS(J)
      R=PCOMP(J,1)**2*(SCOMP(I,1)*h2-SCOMP(I,2)*h1)**2
     .  /(h1**2+h2**2)

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
 3     PROB(14)=PROB(14)+DDIM(R*BRBB(I)*BRBB(J+3)/Rmax,1d0)

*  ee -> hA -> 4taus

       Rmax=1d0
       DO K=1,NhA4tau
        IF(AINT(M2).EQ.hA4tau(K,1) .AND. AINT(M1).EQ.hA4tau(K,2))THEN
         Rmax=hA4tau(K,3)
         GOTO 4
        ENDIF
       ENDDO
 4     PROB(15)=PROB(15)+DDIM(R*BRLL(I)*BRLL(J+3)/Rmax,1d0)

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
 5      PROB(16)=PROB(16)+DDIM(R*BRLL(I)*BRBB(J+3)/Rmax,1d0)
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
 6      PROB(16)=PROB(16)+DDIM(R*BRBB(I)*BRLL(J+3)/Rmax,1d0)
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
 7      PROB(16)=PROB(16)+DDIM(R*BRBB(I)*BRLL(J+3)/Rmax,1d0)
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
 8      PROB(16)=PROB(16)+DDIM(R*BRLL(I)*BRBB(J+3)/Rmax,1d0)
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
 9      PROB(17)=PROB(17)+DDIM(R*BRHAA(I,2*J-1)*BRBB(J+3)**3/Rmax,
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
 10     PROB(18)=PROB(18)+DDIM(R*BRHAA(I,2*J-1)*BRLL(J+3)**3/Rmax,
     .   1d0)
      ENDIF

      ENDIF
      ENDDO
      ENDDO

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

       MA=PMASS(1)
       brtopa1=brtopbh*HCBRWH(4)*BRLL(4)
       PROB(44)=0d0

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

      END


      SUBROUTINE Higgsino_CMS_Trilep(PROB)

      IMPLICIT NONE

      INTEGER I,J,ICHARG,INEUT(2),IChargMAX,ICH,INE(2),IEFF
      INTEGER IX,IY,NX,NY,ND,NM
      PARAMETER (NX=120,NY=80,ND=17,NM=50)

      DOUBLE PRECISION PROB(*)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION brcharsel(2),brcharser(2),brcharsmu1(2),
     .         brcharsmu2(2),brcharstau1(2),brcharstau2(2),
     .         brcharsne1(2),brcharsne2(2),brcharsnm1(2),brcharsnm2(2),
     .         brcharsnt1(2),brcharsnt2(2),brcharsupl(2),brcharsupr(2),
     .         brcharsdownl(2),brcharsdownr(2),brcharst1(2),
     .         brcharst2(2),brcharsb1(2),brcharsb2(2),brcharwneut(2,5),
     .         brcharhcneut(2,5),brcharzchic,brcharHchic(3),
     .         brcharAchic(2),brntaunut(2,5),brnelnue(2,5),
     .         brnmunumu(2,5),brnupdb(2,5),brnchsb(2,5),brntopbb(2,5),
     .         brglupdb(2),brglchsb(2),brgltopbb(2),brchee,brchmumu,
     .         brchtautau,brchnene,brchnmunmu,brchntauntau,brchupup,
     .         brchdodo,brchchch,brchstst,brchtoptop,brchbotbot
      DOUBLE PRECISION brneutst1(5),brneutst2(5),brneutsb1(5),
     .         brneutsb2(5),
     .         brneutsupl(5),brneutsupr(5),brneutsdownl(5),
     .         brneutsdownr(5),brneutsnel(5),brneutsn1(5),
     .         brneutsn2(5),brneutsell(5),brneutselr(5),
     .         brneutsnmu(5),brneutsmu1(5),brneutsmu2(5),
     .         brneutstau1(5),brneutstau2(5),brneutwchar(5,2),
     .         brneuthcchar(5,2),brneutzneut(5,5),
     .         brneutHneut(5,5,3),brneutAneut(5,5,2),brnraddec(5,5)
      DOUBLE PRECISION brneutup(5,5),brneutdow(5,5),brneutch(5,5),
     .         brneutst(5,5),brneutbot(5,5),brneuttop(5,5),
     .         brneutel(5,5),brneutmu(5,5),brneuttau(5,5),
     .         brneutnue(5,5),brneutnumu(5,5),brneutnutau(5,5),
     .         brchubd(5,2),brchcbs(5,2),brchtbb(5,2),brchelne(5,2),
     .         brchmunmu(5,2),brchtauntau(5,2),brglup(5),brgldo(5),
     .         brglch(5),brglst(5),brgltop(5),brglbot(5)
      DOUBLE PRECISION neuttot2(5),neuttot(5),neuttot3(5),neuttotrad(5)
      DOUBLE PRECISION chartot2(2),chartot(2),chartot3(2)
      DOUBLE PRECISION REDCOUP(2),REDCBR(2),mchargino,mneutralino(2)
      DOUBLE PRECISION HCMS(NX,NY),XMIN,XMAX,YMIN,YMAX
      DOUBLE PRECISION DmNeutmCharg(20),mCharg(NM),prosp(ND,NM)
      DOUBLE PRECISION prosp11(2),prosp12(2),prosp21(2),prosp22(2)
      DOUBLE PRECISION xsect(2),aux1,aux2,xsectaver,mchargeff
      DOUBLE PRECISION xsectot,limtrilep
      DOUBLE PRECISION brneuz(5,5), brcharw(2,5)

      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/CHARGINO_BR_2BD/brcharsel,brcharser,brcharsmu1,
     .         brcharsmu2,brcharstau1,brcharstau2,
     .         brcharsne1,brcharsne2,brcharsnm1,brcharsnm2,
     .         brcharsnt1,brcharsnt2,brcharsupl,brcharsupr,
     .         brcharsdownl,brcharsdownr,brcharst1,
     .         brcharst2,brcharsb1,brcharsb2,brcharwneut,
     .         brcharhcneut,brcharzchic,brcharHchic,
     .         brcharAchic
      COMMON/CHARGINO_BR_3BD/brntaunut,brnelnue,brnmunumu,
     .         brnupdb,brnchsb,brntopbb,
     .         brglupdb,brglchsb,brgltopbb,
     .         brchee,brchmumu,brchtautau,brchnene,
     .         brchnmunmu,brchntauntau,brchupup,brchdodo,
     .         brchchch,brchstst,brchtoptop,brchbotbot
      COMMON/NEUTRALINO_BR_2BD/brneutst1,brneutst2,brneutsb1,brneutsb2,
     .         brneutsupl,brneutsupr,brneutsdownl,brneutsdownr,
     .         brneutsnel,brneutsn1,brneutsn2,brneutsell,brneutselr,
     .         brneutsnmu,brneutsmu1,brneutsmu2,
     .         brneutstau1,brneutstau2,brneutwchar,brneuthcchar,
     .         brneutzneut,brneutHneut,brneutAneut,brnraddec
      COMMON/NEUTRALINO_BR_3BD/brneutup,brneutdow,brneutch,brneutst,
     .         brneutbot,brneuttop,brneutel,brneutmu,brneuttau,
     .         brneutnue,brneutnumu,brneutnutau,brchubd,brchcbs,
     .         brchtbb,brchelne,brchmunmu,brchtauntau,brglup,brgldo,
     .         brglch,brglst,brgltop,brglbot
      COMMON/NEUTRALINO_WIDTH/neuttot2,neuttot,neuttot3,neuttotrad
      COMMON/CHARGINO_WIDTH/chartot2,chartot,chartot3
      COMMON/CMSINOS/HCMS,XMIN,XMAX,YMIN,YMAX,DmNeutmCharg,
     .      mCharg,prosp,IChargMAX
      COMMON/XSECTRILEP/xsectot,limtrilep

* Search for the most higgsino-like chargino:
        icharg=1
        IF(u(2,2)**2+v(2,2)**2.ge.u(1,2)**2+v(1,2)**2)
     .    icharg=2
* Search for the most higgsino-like neutralino:
        INEUT(1)=1
        DO I=2,5
          IF(NEU(I,3)**2+NEU(I,4)**2.ge.
     .      NEU(INEUT(1),3)**2+NEU(INEUT(1),4)**2) INEUT(1)=I
        ENDDO
* Search for the next-to-most higgsino-like neutralino:
        INEUT(2)=1
        DO I=2,5
          IF(I.NE.INEUT(1).AND.NEU(I,3)**2+NEU(I,4)**2.ge.
     .      NEU(INEUT(2),3)**2+NEU(INEUT(2),4)**2) INEUT(2)=I
        ENDDO
        mchargino=mcha(icharg)
        mneutralino(1)=dabs(MNEU(INEUT(1)))
        mneutralino(2)=dabs(MNEU(INEUT(2)))

* Reduced couplings squared W-chargino-neutralino(i)
        DO I=1,2
            REDCOUP(I)=(NEU(INEUT(I),4)*U(ICHARG,2)
     .      +NEU(INEUT(I),2)*U(ICHARG,1)*dsqrt(2.D0))**2 
     .      + (NEU(INEUT(I),3)*V(ICHARG,2)
     .      -NEU(INEUT(I),2)*V(ICHARG,1)*dsqrt(2.D0))**2 
* Add possible decays via off-shell Z to brneuz(INEUT(I),1):
          brneuz(INEUT(I),1)=brneutzneut(INEUT(I),1)
     .       +brneutel(INEUT(I),1)*29.7d0
* Add possible decays via off-shell W to brcharw(ICHARG,1):
          brcharw(ICHARG,1)=brcharwneut(ICHARG,1)
     .       +brnelnue(ICHARG,1)/.1071d0
        ENDDO

* Reduced couplings weighted by the BRs, as they contribute to the signal
* relative to pure higgsinos:
        DO I=1,2
          IF(neuttot(INEUT(I)).ne.0d0.AND.chartot(ICHARG).NE.0D0) then
            REDCBR(I)=REDCOUP(I)*brcharw(ICHARG,1)*brneuz(INEUT(I),1)
          ELSE
            REDCBR(I)=0D0
          ENDIF
        ENDDO

        IF(REDCBR(1)+REDCBR(2).EQ.0d0)THEN
         RETURN
        ENDIF

* Find the values of ICH such that Mcharg(ich) is just below mchargino:
       DO I=1,ICHARGMAX-1
         IF(MCHARG(I).LE.MCHARGINO) ICH=I
       ENDDO
       IF(MCHARG(ICHARGMAX).LE.MCHARGINO) THEN
         RETURN
       ELSEIF(MCHARG(1).GE.MCHARGINO) THEN
         RETURN
       ENDIF

* Find the values of INE(I) such that Mneut(INE(I)) is just below mneut(I):
* Use Mneut(J)=MCHARGINO+dmeutmcharg(J)
       DO I=1,2
         DO J=1,ND
           IF(MCHARGINO+DmNeutmCharg(J).LE.mneutralino(I)) INE(I)=J
         ENDDO
         IF(MCHARG(ICH)+DmNeutmCharg(ND).LE.mneutralino(I)) THEN
           RETURN
         ELSEIF(MCHARG(ICH)+DmNeutmCharg(1).GE.mneutralino(I)) THEN
           RETURN
         ENDIF
       ENDDO

* Find the four values of prosp(i,j) just below and above the charg. and neutr. masses:
* Recall prosp(INE,ICH)
       DO I=1,2
         prosp11(I)=prosp(INE(I),ICH)
         prosp12(I)=prosp(INE(I),ICH+1)
         prosp21(I)=prosp(INE(I)+1,ICH)
         prosp22(I)=prosp(INE(I)+1,ICH+1)
* Use the 2dim interpolation formula for xsect(I):
* Along Mchargino:
         aux1=(prosp11(I)*(MCHARG(ICH+1)-mchargino)
     .         +prosp12(I)*(mchargino-MCHARG(ICH)))
     .         /(MCHARG(ICH+1)-MCHARG(ICH))
         aux2=(prosp21(I)*(MCHARG(ICH+1)-mchargino)
     .         +prosp22(I)*(mchargino-MCHARG(ICH)))
     .         /(MCHARG(ICH+1)-MCHARG(ICH))
* Along Mneutrino:
         xsect(I)=
     .     (aux1*(MCHARGINO+DmNeutmCharg(INE(I)+1)-mneutralino(I))
     .      +aux2*(mneutralino(I)-MCHARGINO-DmNeutmCharg(INE(I))))
     .       /(10d0)
      ENDDO

* Find the average cross section, weighted by the reduced coupl.*BRs:
      xsectaver=(REDCBR(1)*xsect(1)+REDCBR(2)*xsect(2))
     .          /(REDCBR(1)+REDCBR(2))

* Find mchargeff which reproduces this cross section
* for degenerate charginos/neutralinos:
      IEFF=0
      DO I=1,NM
        IF(PROSP(9,I).GE.xsectaver) IEFF=I
      ENDDO
      mchargeff=MCHARG(IEFF)+(MCHARG(IEFF+1)-MCHARG(IEFF))*
     .   (xsectaver-PROSP(9,IEFF))/(PROSP(9,IEFF+1)-PROSP(9,IEFF))
      xsectot=REDCBR(1)*xsect(1)+REDCBR(2)*xsect(2)
* Rescaling twiki/prospino:
      xsectot=xsectot*(0.936d0 + 3.87d-4*mchargeff)

* Excluded cross section
      IX=1+INT(DFLOAT(NX)*(mchargeff-XMIN)/(XMAX-XMIN))
      IY=1+INT(DFLOAT(NY)*(DABS(MNEU(1))-YMIN)/(YMAX-YMIN))

* Use the 2dim interpolation formula for limtrilep:
* Along Mchargino:
        aux1=(HCMS(IX,IY)*
     .   (DFLOAT(IX)-DFLOAT(NX)*(mchargeff-XMIN)/(XMAX-XMIN))
     .         +HCMS(IX+1,IY)*
     .   (DFLOAT(NX)*(mchargeff-XMIN)/(XMAX-XMIN)-DFLOAT(IX)+1d0))
        aux2=(HCMS(IX,IY+1)*
     .   (DFLOAT(IX)-DFLOAT(NX)*(mchargeff-XMIN)/(XMAX-XMIN))
     .         +HCMS(IX+1,IY+1)*
     .   (DFLOAT(NX)*(mchargeff-XMIN)/(XMAX-XMIN)-DFLOAT(IX)+1d0))
* Along Mneutralino:
        limtrilep=
     .    aux1*
     .     (DFLOAT(IY)-DFLOAT(NY)*(DABS(MNEU(1))-YMIN)/(YMAX-YMIN))
     .    +aux2*
     .     (DFLOAT(NY)*(DABS(MNEU(1))-YMIN)/(YMAX-YMIN)-DFLOAT(IY)+1d0)

      IF(IX.GE.1 .AND. IX.LE.NX .AND. IY.GE.1 .AND. IY.LE.NY)THEN
       IF(xsectot.GT.limtrilep)THEN
        PROB(64)=xsectot/limtrilep-1d0
       ENDIF
      ENDIF

      RETURN

      END
