      SUBROUTINE HIDECAY_CPV()

*         Higgs decays
*        calculates the Higgs BRs
*
*      CU,CD,CV,CJ,CG(i) Reduced scalar couplings of hi (i=1..5) to up type 
*                        fermions, down type fermions, gauge bosons, gluons 
*                        and photons
*      CUP,CDP,CJP,CGP(i) idem for pseudoscalar component
*      CB,CBP(I)         Reduced couplings of hi (i=1..5) to b-quarks 
*                        including DELMB corrections
*
*      WIDTH(i) Total decay width of hi (i=1..5)with the following 
*               branching ratios:
*      BRJJ(i)   hi (i=1..5)  -> hadrons
*      BRMM(i)        "       -> mu mu
*      BRLL(i)        "       -> tau tau
*      BRCC(i)        "       -> cc
*      BRBB(i)        "       -> bb
*      BRTT(i)        "       -> tt
*      BRWW(i)        "       -> WW 
*      BRZZ(i)        "       -> ZZ 
*      BRGG(i)        "       -> gamma gamma
*      BRZG(i)        "       -> Z gamma
*      BRHIGGS(i)     "       -> other Higgses, including:
*        BRHCHC(i)    "       -> h+h-
*        BRHAZ(i,j)   "       -> Zhj  (j=1..4)
*        BRHCW(i)     "       -> h+W- 
*        BRHHH(i,j)   "       -> h1h1 (1), h1h2 (2), h2h2 (3), h1h3 (4), 
*                                h2h3 (5), h3h3 (6), h1h4 (7), h2h4 (8),
*                                h3h4 (9), h4h4 (10)
*      BRSUSY(i)      "       -> susy particles, including:
*        BRNEU(i,j,k)         -> neutralinos j,k (j,k=1..5)
*        BRCHA(i,j)           -> charginos 11, 12, 22 (j=1..3)
*        BRHSQ(i,j)   "       -> uLuL, uRuR, dLdL, dRdR, t1t1, t2t2,
*                                t1t2, b1b1, b2b2, b1b2 (j=1..10)
*        BRHSL(i,j)   "       -> lLlL, lRlR, nLnL, l1l1, l2l2, l1l2,
*                                ntnt (i=1..3, j=1..7)
*
*      HCWIDTH  Total decay width of the charged Higgs
*               with the following branching ratios:
*      HCBRM         h+ -> mu nu_mu
*      HCBRL         "  -> tau nu_tau
*      HCBRSU        "  -> s u
*      HCBRBU        "  -> b u
*      HCBRSC        "  -> s c
*      HCBRBC        "  -> b c
*      HCBRBT        "  -> b t
*      HCBRWHT       "  -> neutral Higgs W+, including:
*        HCBRWH(i)   "  -> hiW+ (i=1..5)
*      HCBRSUSY      "  -> susy particles,including
*        HCBRNC(i,j) "  -> neutralino i chargino j (i=1..5, j=1..2)
*        HCBRSQ(i)   "  -> uLdL, t1b1, t1b2, t2b1, t2b2 (i=1..5)
*        HCBRSL(i)   "  -> lLnL, t1nt, t2nt (i=1..3)
*

      IMPLICIT NONE

      INTEGER I,J,M,NF,N0,NFGG,NFEXT,VFLAG

      DOUBLE PRECISION SQR2,PI,C2TW,T2TW,HIGTOP,ASG,ASH,AS3,AS4,ASMT
      DOUBLE PRECISION EPS,FQCD,SQCD,XFAC,X,Y,RATCOUP,RAT,RMTTOP
      DOUBLE PRECISION HVV,HV,HFF,QCd0,HQCDM,HQCD,QCDH,TQCDH,HGGQCD
      DOUBLE PRECISION AFF,AQCDM,AQCD,QCDA,TQCDA,AGGQCD,SGGQCD
      DOUBLE PRECISION QCDC,QCDCM,CQCD,CQCDM,QCDCI,QCDCMI
      DOUBLE PRECISION BETA,LAMB,SP,ALPHAS,RUNM,QQINT,FINT
      DOUBLE PRECISION T,Z,XI,BIJ,CFF,LQ,MH,Ytau
      DOUBLE PRECISION RMS,RMC,RMB,RMT,RUNMB
      DOUBLE PRECISION HJJ,HEE,HMM,HLL,HSS,HCC,HBB,HTT,HWW,HZZ,HGG
      DOUBLE PRECISION HZG,HS1,HS2,HC1,HC2,HB1,HB2,HT1,HT2,DCC,DBB
      DOUBLE PRECISION DLU,DLD,XM1,XM2,CWW,CZZ,XX(4),YY(4)
      DOUBLE PRECISION FCH1,FCH2,HTWW,HTZZ,FT,FB,ACOUP
      DOUBLE PRECISION HHH(10),HHCHC,HTOT,CH,RH,HAZ(4),HHCW
      DOUBLE PRECISION HNEU(5,5),HCHA(3),HSQ(10),HSL(7),STOT
      DOUBLE PRECISION HMN,HLN,HSU,HSU1,HSU2,HSC,HSC1,HSC2
      DOUBLE PRECISION HBC,HBC1,HBC2,HBU,HBU1,HBU2,HBT,HBT1,HBT2
      DOUBLE PRECISION HCWH(5),HCNC(5,2),HCSQ(5),HCSL(3)

      DOUBLE COMPLEX CTT,CTB,CTC,CTL,CTW,CTHC,CTCH1,CTCH2,CTM,CTE
      DOUBLE COMPLEX CXT,CXB,CXC,CXL,CXW,CXHC,CXCH1,CXCH2,CXM,CXE
      DOUBLE COMPLEX CXTP,CXBP,CXCP,CXLP,CXCH1P,CXCH2P,CXMP,CXEP
      DOUBLE COMPLEX CTUL,CTUR,CTDL,CTDR,CTST1,CTST2,CTSB1,CTSB2
      DOUBLE COMPLEX CXUL,CXUR,CXDL,CXDR,CXST1,CXST2,CXSB1,CXSB2
      DOUBLE COMPLEX CTLL,CTLR,CTSL1,CTSL2,CXLL,CXLR,CXSL1,CXSL2
      DOUBLE COMPLEX CLT,CLB,CLC,CLW,CLH,CXTZ,CXBZ,CXWZ,CXHZ,CXCZ
      DOUBLE COMPLEX CI1,CI2,CGZ,CF,CA,CBC,CM
      DOUBLE COMPLEX CLCH1,CLCH2,CXCH1Z,CXCH2Z

      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION ALEM0
      DOUBLE PRECISION XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION VUS,VCB,VUB
      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      DOUBLE PRECISION l,k,Alcos1,Akcos2,muq,nuq
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION DELMB,DELML,DEL1
      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION MCH2(2),U(2,2,2),V(2,2,2)
      DOUBLE PRECISION MNEU(5),NEU(5,5,2)
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
      DOUBLE PRECISION CU(5),CUP(5),CD(5),CDP(5),CB(5),CBP(5),CJ(5),
     . CJP(5),CI(5),CG(5),CGP(5),CV(5),CZG(5),CZGP(5),CL(5),CLP(5)
      DOUBLE PRECISION WIDTH(5),HCWIDTH
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),
     . BRCC(5),BRBB(5),BRTT(5),BRWW(5),BRZZ(5),BRGG(5),BRZG(5)
      DOUBLE PRECISION BRHHH(5,10),BRHCHC(5),BRHAZ(5,4),BRHCW(5),
     . BRHIGGS(5)
      DOUBLE PRECISION BRNEU(5,5,5),BRCHA(5,3),BRHSQ(5,10),BRHSL(5,7),
     . BRSUSY(5)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,HCBRBT
      DOUBLE PRECISION HCBRWH(5),HCBRWHT
      DOUBLE PRECISION HCBRNC(5,2),HCBRSQ(5),HCBRSL(3),HCBRSUSY
      DOUBLE PRECISION GamHGAGA,GamHee,GamHmumu,GamHtata,GamHhadr,
     . GamHcc,GamHbb,GamHinv(3),GamHWW,GamHZZ,GamHAA
* New May 2019:
      DOUBLE PRECISION ZETA2,ZETA3,HGGQCD2,FQCD0,AGGQCD2
* End New

      COMMON/STSBSCALE/QSTSB
      COMMON/ALEM0/ALEM0
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/CKM/VUS,VCB,VUB
      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/SMFERM/mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      COMMON/QPAR/l,k,Alcos1,Akcos2,muq,NUQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/CHASPEC/MCH2,U,V
      COMMON/NEUSPEC/MNEU,NEU
      COMMON/SFERM3SPEC/MST2,UT,MSB2,UB,MSL2,UTAU,MSNT2
      COMMON/SFERM1SPEC/MSU2,MSD2,MSE2,MSNE2,MSMU2,UMU
      COMMON/SFERMPSPEC/MST2P,MSB2P,MSU2P,MSD2P
      COMMON/HISFCOUP/GRHSTST,GRHSBSB,GRHSLSL,GRHSUSU,GRHSDSD,
     . GRHSESE,GRHSNSN,GIHSTST,GIHSBSB,GIHSLSL,GRHCSTSB,GRHCSNSL,
     . GRHCSUSD,GRHCSNSE,GIHCSTSB,GIHCSNSL
      COMMON/HINOCOUP/COH0CH,COH0NEU,COHPNEUCHM,COHMNEUCHP
      COMMON/HICOUP/gH0H0H0,gRH0HPHM,gIH0HPHM
      COMMON/HNSMCOUP/CU,CUP,CD,CDP,CB,CBP,CJ,CJP,CI,CG,CGP,CV,CZG,CZGP,
     . CL,CLP
      COMMON/HIWIDTH/WIDTH,HCWIDTH
      COMMON/HNSMBR/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG
      COMMON/HNHIBR/BRHHH,BRHCHC,BRHAZ,BRHCW,BRHIGGS
      COMMON/HNSUSYBR/BRNEU,BRCHA,BRHSQ,BRHSL,BRSUSY
      COMMON/HCSMBR/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,HCBRBT
      COMMON/HCHIBR/HCBRWH,HCBRWHT
      COMMON/HCSUSYBR/HCBRNC,HCBRSQ,HCBRSL,HCBRSUSY
      COMMON/VFLAG/VFLAG
      COMMON/LIGHTHDECAYS/GamHGAGA,GamHee,GamHmumu,GamHtata,GamHhadr,
     . GamHcc,GamHbb,GamHinv,GamHWW,GamHZZ,GamHAA

      QQINT(RAT,X,Y)= RAT**2*X+(1d0-RAT**2)*Y
      BETA(X)= DSQRT(1d0-4d0*X)
      LAMB(X,Y)= DSQRT((1d0-X-Y)**2-4d0*X*Y)
      CM(X)= DCMPLX(4d0*MIN(1d3,X)**2,-EPS)
      CF(CA)= -CDLOG(-(1d0+CDSQRT(1d0-CA))
     . / (1d0-CDSQRT(1d0-CA)))**2/4d0
      CGZ(CA)= CDSQRT(1d0-CA)/2d0*CDLOG(-(1d0+CDSQRT(1d0-CA))
     . / (1d0-CDSQRT(1d0-CA)))
      CI1(CA,CBC)= CA*CBC/2d0/(CA-CBC)
     . + CA**2*CBC**2/2/(CA-CBC)**2*(CF(CA)-CF(CBC))
     . + CA**2*CBC/(CA-CBC)**2*(CGZ(CA)-CGZ(CBC))
      CI2(CA,CBC)= -CA*CBC/2d0/(CA-CBC)*(CF(CA)-CF(CBC))
      HV(X)= 3d0*(1d0-8d0*X+20d0*X**2)/DSQRT((4d0*X-1d0))
     . * DACOS((3d0*X-1d0)/2d0/DSQRT(X**3))
     . - (1d0-X)*(47d0/2d0*X-13d0/2d0+1d0/X)
     . - 3d0/2d0*(1d0-6d0*X+4d0*X**2)*DLOG(X)
      HVV(X,Y)= GF/(4d0*PI*SQR2)*X**3/2d0*BETA(Y)
     . * (1d0-4d0*Y+12d0*Y**2)
      HFF(X,Y)= GF/(4d0*PI*SQR2)*X**3*Y*(BETA(Y))**3
      AFF(X,Y)= GF/(4d0*PI*SQR2)*X**3*Y*(BETA(Y))
      CFF(Z,T,X,Y)= GF/(4d0*PI*SQR2)*Z**3*LAMB(X,Y)
     . * ((1d0-X-Y)*(X*T**2+Y/T**2)-4d0*X*Y)
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
      HGGQCD(ASG,NF)= 1d0+ASG/PI*(95d0/4d0-NF*7d0/6d0)
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
      AGGQCD(ASG,NF)= 1d0+ASG/PI*(97d0/4d0-NF*7d0/6d0)
* July 2010:
c      SGGQCD(ASG)= ASG/PI*17d0/6d0
* New May 2019:
      AGGQCD2(ASG,NF,MH,MT)=1d0+ASG/PI*(97d0/4d0-NF*7d0/6d0)
     . +(ASG/PI)**2*(237311d0/864d0-529d0/24d0*ZETA2-445d0/8d0*ZETA3
     . +5d0*DLOG(MH**2/MT**2))
      SGGQCD(ASG)=ASG/PI*7d0/2d0
* End New
      XI(X,Y)= 2d0*X/(1d0-X-Y+LAMB(X,Y))
      BIJ(X,Y)= (1d0-X-Y)/LAMB(X,Y)
     . * (4d0*SP(XI(X,Y)*XI(Y,X))
     . - 2d0*SP(-XI(X,Y))-2d0*SP(-XI(Y,X))
     . + 2d0*DLOG(XI(X,Y)*XI(Y,X))*DLOG(1d0-XI(X,Y)*XI(Y,X))
     . - DLOG(XI(X,Y))*DLOG(1d0+XI(X,Y))
     . - DLOG(XI(Y,X))*DLOG(1d0+XI(Y,X)))
     . - 4d0*(DLOG(1d0-XI(X,Y)*XI(Y,X))
     . + XI(X,Y)*XI(Y,X)/(1d0-XI(X,Y)*XI(Y,X))*DLOG(XI(X,Y)*XI(Y,X)))
     . + (LAMB(X,Y)+X-Y)/LAMB(X,Y)*(DLOG(1d0+XI(X,Y))
     . - XI(X,Y)/(1d0+XI(X,Y))*DLOG(XI(X,Y)))
     . + (LAMB(X,Y)-X+Y)/LAMB(X,Y)*(DLOG(1d0+XI(Y,X))
     . - XI(Y,X)/(1d0+XI(Y,X))*DLOG(XI(Y,X)))
      QCDC(X,Y)= 1d0+4d0/3d0*ASH/PI*(9d0/4d0+(3d0-2d0*X+2d0*Y)
     . / 4d0*DLOG(X/Y)+((1.5d0-X-Y)*LAMB(X,Y)**2+5d0*X*Y)/2d0
     . / LAMB(X,Y)/(1d0-X-Y)*DLOG(XI(X,Y)*XI(Y,X))+BIJ(X,Y))
     . + ASH/PI*(2d0*(4d0/3d0-DLOG(X))
     . - (X*2d0*(4d0/3d0-DLOG(X))+Y*2d0*(4d0/3d0-DLOG(Y)))
     . / (1d0-X-Y)-(X*2d0*(4d0/3d0-DLOG(X))*(1d0-X+Y)
     . + Y*2d0*(4d0/3d0-DLOG(Y))*(1d0+X-Y))/LAMB(X,Y)**2)
      QCDCI(X,Y)= 1d0+4d0/3d0*ASH/PI*(3d0+(Y-X)/2d0*DLOG(X/Y)
     . + (2d0*(1d0-X-Y)+LAMB(X,Y)**2)/2d0/LAMB(X,Y)
     . * DLOG(XI(X,Y)*XI(Y,X))+BIJ(X,Y))
     . + ASH/PI*(2d0*(4d0/3d0-DLOG(X))+2d0*(4d0/3d0-DLOG(Y))
     . - (X*2d0*(4d0/3d0-DLOG(X))*(1d0-X+Y)
     . + Y*2d0*(4d0/3d0-DLOG(Y))*(1d0+X-Y))/LAMB(X,Y)**2)
      QCDCM(X,Y)= 1d0+4d0/3d0*ASH/PI*(9d0/4d0
     . + (3d0-2d0*X+2d0*Y)/4d0*DLOG(X/Y)+((1.5d0-X-Y)
     . * LAMB(X,Y)**2+5d0*X*Y)/2d0/LAMB(X,Y)/(1d0-X-Y)
     . * DLOG(4d0*X*Y/(1d0-X-Y+LAMB(X,Y))**2)
     . + BIJ(X,Y))
      QCDCMI(X,Y)= 1d0+4d0/3d0*ASH/PI*(3d0+(Y-X)/2d0*DLOG(X/Y)
     . + (2d0*(1d0-X-Y)*LAMB(X,Y)**2)/2d0/LAMB(X,Y)
     . * DLOG(4d0*X*Y/(1d0-X-Y+LAMB(X,Y))**2)
     . + BIJ(X,Y))
      CQCD(Z,T,X,Y)= GF/(4d0*PI*SQR2)*Z**3*LAMB(X,Y)
     . * ((1d0-X-Y)*(X*T**2*QCDC(Y,X)
     . + Y/T**2*QCDC(X,Y))
     . - 4d0*X*Y*QCDCI(X,Y))
      CQCDM(Z,T,X,Y)= GF/(4d0*PI*SQR2)*Z**3*LAMB(X,Y)
     . * ((1d0-X-Y)*(X*T**2*QCDCM(Y,X)
     . + Y/T**2*QCDCM(X,Y))
     . - 4d0*X*Y*QCDCMI(X,Y))

      EPS= 1d-8
      PI= 4d0*DATAN(1d0)
      SQR2= DSQRT(2d0)
* New May 2019:
      ZETA2=PI**2/6d0
      ZETA3=1.202056903159594d0
* End New

*   Number of light flavours included in the gluonic decays
*   Higgs -> gg* -> gqq (see hdecay): NFGG = 5
* New May 2019:
      NFGG= 5
* End New

*   Weak angle theta_W (S2TW = sin(theta_W)**2):
      C2TW= 1d0-S2TW
      T2TW= S2TW/C2TW

*   Higgs rescaled rotation matrix

c      DO I=1,5
c       XHG(I,1)=XH(I,1)/dsqrt(ZHU)
c       XHG(I,2)=XH(I,2)/dsqrt(ZHD)
c       XHG(I,3)=XH(I,3)/dsqrt(ZS)
c       XHG(I,4)=XH(I,4)*cosb/dsqrt(ZHU)
c       XHG(I,5)=XH(I,4)*sinb/dsqrt(ZHD)
c       XHG(I,6)=XH(I,5)/dsqrt(ZS)
c      ENDDO

*   Alpha_s at the top pole mass scales, used for the running
*   Yukawa coupling ht and running quark masses RMT below
*   NOTE: MT = top pole mass
* New May 2019:
      ASMT= ALPHAS(MT,3)
* End New
      
*   Tau Yukawa coupling for Higgs-stau coupling
      Ytau=mtau/vdq

*   MT = Top pole mass; RMTTOP = running mass at Mtop (MS_bar):
      RMTTOP= MT/(1d0+4d0*ASMT/(3d0*PI)+11d0*(ASMT/PI)**2)

* Loop over the neutral Higgs bosons

      DO I=1,5

       MH= dsqrt(MH0(I))

*  Running quark masses at MH

       RMS= RUNM(MH,3)
       RMC= RUNM(MH,4)
       RMB= RUNMB(MH)
       IF(MH.GE.MT)THEN
        RMT= RMTTOP
     .   *(1d0+7d0/(4d0*PI)*ASMT*DLOG(MH**2/MT**2))
     .   **(-4d0/7d0)
       ELSE
         RMT= RMTTOP
     .   *(1d0+23d0/(12d0*PI)*ASMT*DLOG(MH**2/MT**2))
     .   **(-12d0/23d0)
       ENDIF

*   Log for rad. corrs to Higgs self couplings
*   (QSTSB was initialized in the subroutine RUNPAR):

       LQ=DLOG(MAX(QSTSB,MH**2)/(MAX(MT,MH)**2))

*  Strong coupling constant at MH

       HIGTOP= (MH/MT)**2
       MT0= 3d8
* New May 2019:
       ASH= ALPHAS(MAX(1d0,MH),3)
       MC0= 1d8
       MB0= 2d8
       AS3= ALPHAS(MAX(1d0,MH),3)
       MC0= MC
       AS4= ALPHAS(MAX(1d0,MH),3)
* End New
       MB0= MBP
       MT0= MT

*  Scalar couplings relative to the standard model

       CV(I)= XH(I,1)*SINB+XH(I,2)*COSB
       CU(I)= XH(I,1)/SINB
       CD(I)= XH(I,2)/COSB
       CB(I)=(CD(I)+DELMB*(CU(I)+DSQRT(vuq**2+vdq**2)*l/muq*XH(I,3)))
     .       /(1d0+DELMB)
       CL(I)=(CD(I)+DELML*(CU(I)+DSQRT(vuq**2+vdq**2)*l/muq*XH(I,3)))
     .       /(1d0+DELML)

*  Pseudoscalar couplings relative to the standard model

       CUP(I)= XH(I,4)/TANB
       CDP(I)= XH(I,4)*TANB
       CBP(I)= CDP(I)/(1d0+DELMB)
       CBP(I)=(CDP(I)-DELMB*(CUP(I)+DSQRT(vuq**2+vdq**2)*l/muq*XH(I,5)))
     .         /(1d0+DELMB)
       CLP(I)=(CDP(I)-DELML*(CUP(I)+DSQRT(vuq**2+vdq**2)*l/muq*XH(I,5)))
     .         /(1d0+DELML)

*  Effective coupling to 2 gluons/2 photons - scalar part

       CTT= CM(MT/MH)
       CTB= CM(MBP/MH)
       CTC= CM(MC/MH)
       CTL= CM(MTAU/MH)
       CTM= CM(MMU/MH)
       CTE= CM(MEL/MH)
       CTW= CM(MW/MH)
       CTHC= CM(DSQRT(MHC)/MH)
       CTCH1= CM(DSQRT(MCH2(1))/MH)
       CTCH2= CM(DSQRT(MCH2(2))/MH)
       CTUL= CM(DSQRT(MSU2(1))/MH)
       CTUR= CM(DSQRT(MSU2(2))/MH)
       CTDL= CM(DSQRT(MSD2(1))/MH)
       CTDR= CM(DSQRT(MSD2(2))/MH)
       CTST1= CM(DSQRT(MST2(1))/MH)
       CTST2= CM(DSQRT(MST2(2))/MH)
       CTSB1= CM(DSQRT(MSB2(1))/MH)
       CTSB2= CM(DSQRT(MSB2(2))/MH)
       CTLL= CM(DSQRT(MSE2(1))/MH)
       CTLR= CM(DSQRT(MSE2(2))/MH)
       CTSL1= CM(DSQRT(MSL2(1))/MH)
       CTSL2= CM(DSQRT(MSL2(2))/MH)
       CXT= 2d0*CTT*(1d0+(1d0-CTT)*CF(CTT))
       CXB= 2d0*CTB*(1d0+(1d0-CTB)*CF(CTB))
       CXC= 2d0*CTC*(1d0+(1d0-CTC)*CF(CTC))
       CXL= 2d0*CTL*(1d0+(1d0-CTL)*CF(CTL))
       CXM= 2d0*CTM*(1d0+(1d0-CTM)*CF(CTM))
       CXE= 2d0*CTE*(1d0+(1d0-CTE)*CF(CTE))
       CXW= -(2d0+3d0*CTW+3d0*CTW*(2d0-CTW)*CF(CTW))
       CXHC= CTHC*(CTHC*CF(CTHC)-1d0)
       CXCH1= 2d0*CTCH1*(1d0+(1d0-CTCH1)*CF(CTCH1))
       CXCH2= 2d0*CTCH2*(1d0+(1d0-CTCH2)*CF(CTCH2))
       CXUL= CTUL*(CTUL*CF(CTUL)-1d0)/(DSQRT(SQR2*GF)*MSU2(1))
     .      *GRHSUSU(I,1,1)
       CXUR= CTUR*(CTUR*CF(CTUR)-1d0)/(DSQRT(SQR2*GF)*MSU2(2))
     .      *GRHSUSU(I,2,2)
       CXDL= CTDL*(CTDL*CF(CTDL)-1d0)/(DSQRT(SQR2*GF)*MSD2(1))
     .      *GRHSDSD(I,1,1)
       CXDR= CTDR*(CTDR*CF(CTDR)-1d0)/(DSQRT(SQR2*GF)*MSD2(2))
     .      *GRHSDSD(I,2,2)
       CXST1= CTST1*(CTST1*CF(CTST1)-1d0)
     .       /(2d0*DSQRT(SQR2*GF)*MST2(1))
     .       *DCMPLX(GRHSTST(I,1,1),GIHSTST(I,1,1))
       CXST2= CTST2*(CTST2*CF(CTST2)-1d0)
     .       /(2d0*DSQRT(SQR2*GF)*MST2(2))
     .       *DCMPLX(GRHSTST(I,2,2),GIHSTST(I,2,2))
       CXSB1= CTSB1*(CTSB1*CF(CTSB1)-1d0)
     .       /(2d0*DSQRT(SQR2*GF)*MSB2(1))
     .       *DCMPLX(GRHSBSB(I,1,1),GIHSBSB(I,1,1))
       CXSB2= CTSB2*(CTSB2*CF(CTSB2)-1d0)
     .       /(2d0*DSQRT(SQR2*GF)*MSB2(2))
     .       *DCMPLX(GRHSBSB(I,2,2),GIHSBSB(I,2,2))
       CXLL= CTLL*(CTLL*CF(CTLL)-1d0)/(DSQRT(SQR2*GF)*MSE2(1))
     .      *GRHSESE(I,1,1)
       CXLR= CTLR*(CTLR*CF(CTLR)-1d0)/(DSQRT(SQR2*GF)*MSE2(2))
     .      *GRHSESE(I,2,2)
       CXSL1= CTSL1*(CTSL1*CF(CTSL1)-1d0)
     .       /(2d0*DSQRT(SQR2*GF)*MSL2(1))
     .       *DCMPLX(GRHSLSL(I,1,1),GIHSLSL(I,1,1))
       CXSL2= CTSL2*(CTSL2*CF(CTSL2)-1d0)
     .       /(2d0*DSQRT(SQR2*GF)*MSL2(2))
     .       *DCMPLX(GRHSLSL(I,2,2),GIHSLSL(I,2,2))

*   Here CJ and CG are the actual couplings. Later we divide
*   by the SM couplings in order to obtain reduced couplings

       CJ(I)= CDABS(CU(I)*(CXT+CXC)+CB(I)*CXB
     .       +CXUL+CXUR+CXDL+CXDR+CXST1+CXST2+CXSB1+CXSB2)

       CI(I)= DREAL(DCONJG(CU(I)*(CXT+CXC)+CB(I)*CXB
     .       +CXUL+CXUR+CXDL+CXDR+CXST1+CXST2+CXSB1+CXSB2)
     .      *(CXUL+CXUR+CXDL+CXDR+CXST1+CXST2+CXSB1+CXSB2))

       CG(I)= CDABS(4d0/3d0*CU(I)*(CXT+CXC) + CB(I)*CXB/3d0
     .      + CL(I)*CXL + CD(I)*(CXM+CXE) + CV(I)*CXW
     .      + GRH0HPHM(I,2,2)/(2d0*MHC*DSQRT(SQR2*GF))*CXHC
     .      + COH0CH(I,1,1,1)/DSQRT(SQR2*GF*MCH2(1))*CXCH1
     .      + COH0CH(I,2,2,1)/DSQRT(SQR2*GF*MCH2(2))*CXCH2
     .      + 4d0/3d0*(CXUL+CXUR+CXST1+CXST2)
     .      + 1d0/3d0*(CXDL+CXDR+CXSB1+CXSB2)
     .      + CXLL+CXLR+CXSL1+CXSL2)

*  Effective coupling to 2 gluons/2 photons - pseudoscalar part

       CXTP= 2d0*CTT*CF(CTT)
       CXBP= 2d0*CTB*CF(CTB)
       CXCP= 2d0*CTC*CF(CTC)
       CXLP= 2d0*CTL*CF(CTL)
       CXMP= 2d0*CTM*CF(CTM)
       CXEP= 2d0*CTE*CF(CTE)
       CXCH1P= 2d0*CTCH1*CF(CTCH1)
       CXCH2P= 2d0*CTCH2*CF(CTCH2)

*   Here CJP and CGP are the actual couplings. Later we divide
*   by the SM couplings in order to obtain  reduced couplings

       CJP(I)= CDABS(CUP(I)*(CXTP+CXCP) + CBP(I)*CXBP)

       CGP(I)= CDABS(4d0/3d0*CUP(I)*(CXTP+CXCP)
     .  + CBP(I)*CXBP/3d0 + CLP(I)*CXLP + CDP(I)*(CXMP+CXEP)
     .  + COH0CH(I,1,1,2)/DSQRT(SQR2*GF*MCH2(1))*CXCH1P
     .  + COH0CH(I,2,2,2)/DSQRT(SQR2*GF*MCH2(2))*CXCH2P)

*   LIGHT STATE

         IF(MH.LE.50d0)THEN

*    CJ(P)/CG(P) becomes the REDUCED couplings to 2 gluons/photons
*    as compared to the couplings of a SM Higgs with same mass
*    (CJ(P) does not contain strange quark loops)

       NFEXT= 3
       ASG= AS3
       FQCD= HGGQCD(ASG,NFEXT)
       SQCD= SGGQCD(ASG)
       CJ(I)= DSQRT(MAX(0d0,CJ(I)**2+CI(I)*SQCD/FQCD))/
     .        CDABS(CXT+CXC+CXB)
       CJP(I)= CJP(I)/CDABS(CXT+CXC+CXB)
       CG(I)= CG(I)/CDABS(4d0/3d0*(CXT+CXC)+CXB/3d0+CXL+CXM+CXE+CXW)
       CGP(I)= CGP(I)/CDABS(4d0/3d0*(CXT+CXC)+CXB/3d0+CXL+CXM+CXE+CXW)

       CALL HDECAY_lightCPV(I)

       WIDTH(I)= GamHGAGA+GamHee+GamHmumu+GamHtata+GamHhadr+GamHcc
     . +GamHbb+GamHinv(1)+GamHinv(2)+GamHinv(3)+GamHWW+GamHZZ+GamHAA
       BRJJ(I)= GamHhadr/WIDTH(I)
       BREE(I)= GamHee/WIDTH(I)
       BRMM(I)= GamHmumu/WIDTH(I)
       BRLL(I)= GamHtata/WIDTH(I)
       BRCC(I)= GamHcc/WIDTH(I)
       BRBB(I)= GamHbb/WIDTH(I)
       BRGG(I)= GamHGAGA/WIDTH(I)
       BRNEU(I,1,1)= GamHinv(1)/WIDTH(I)
       BRNEU(I,1,2)= GamHinv(2)/WIDTH(I)
       BRNEU(I,2,1)= BRNEU(I,1,2)
       BRNEU(I,2,2)= GamHinv(3)/WIDTH(I)
       BRSUSY(I)= 0d0
       DO J=1,2
       DO M=J,2
        BRSUSY(I)= BRSUSY(I)+BRNEU(I,J,M)
       ENDDO
       ENDDO
       BRWW(I)= GamHWW/WIDTH(I)
       BRZZ(I)= GamHZZ/WIDTH(I)
       BRTT(I)= 0d0
       BRZG(I)= 0d0
       BRHHH(I,1)= GamHAA/WIDTH(I)
       BRHIGGS(I)= BRHHH(I,1)
       DO J=2,10
        BRHHH(I,J)= 0d0
       ENDDO
       BRHCHC(I)= 0d0
       DO J=1,2
        BRHAZ(I,J)= 0d0
       ENDDO
       BRHCW(I)=  0d0
       DO J=1,5
       DO M=3,5
        BRNEU(I,J,M)= 0d0
        BRNEU(I,M,J)= 0d0
       ENDDO
       ENDDO
       DO J=1,3
        BRCHA(I,J)= 0d0
       ENDDO
       DO J=1,10
        BRHSQ(I,J)= 0d0
       ENDDO
       DO J=1,7
        BRHSL(I,J)= 0d0
       ENDDO

*   HEAVY STATE

         ELSE

*  Partial widths

*   h -> gg

       NFEXT= 3
       ASG= AS3
       FQCD= HGGQCD(ASG,NFEXT)
       SQCD= SGGQCD(ASG)
       XFAC= MAX(0d0,CJ(I)**2*FQCD+CI(I)*SQCD)
       HJJ= GF/(64d0*PI*SQR2)*MH**3*(ASG/PI)**2*XFAC
       FQCD= AGGQCD(ASG,NFEXT)
       XFAC= CJP(I)**2*FQCD
       HJJ= HJJ+GF/(64d0*PI*SQR2)*MH**3*(ASG/PI)**2*XFAC

*   h -> gg* -> gcc to be added to h -> cc

       NFEXT= 4
       ASG= AS4
       FQCD= HGGQCD(ASG,NFEXT)
       SQCD= SGGQCD(ASG)
       XFAC= MAX(0d0,CJ(I)**2*FQCD+CI(I)*SQCD)
       DCC= GF/(64d0*PI*SQR2)*MH**3*(ASG/PI)**2*XFAC
       FQCD= AGGQCD(ASG,NFEXT)
       XFAC= CJP(I)**2*FQCD
       DCC= DCC+GF/(64d0*PI*SQR2)*MH**3*(ASG/PI)**2*XFAC-HJJ

*   h -> gg* -> gbb to be added to h -> bb

       NFEXT= 5
       ASG= ASH
       FQCD= HGGQCD(ASG,NFEXT)
       SQCD= SGGQCD(ASG)
       XFAC= MAX(0d0,CJ(I)**2*FQCD+CI(I)*SQCD)
       DBB= GF/(64d0*PI*SQR2)*MH**3*(ASG/PI)**2*XFAC
* New July 2019:
       HJJ= DBB
* End new
       FQCD= AGGQCD(ASG,NFEXT)
       XFAC= CJP(I)**2*FQCD
       DBB= DBB+GF/(64d0*PI*SQR2)*MH**3*(ASG/PI)**2*XFAC-HJJ-DCC
* New July 2019:
       HJJ= HJJ+GF/(64d0*PI*SQR2)*MH**3*(ASG/PI)**2*XFAC

C  H ---> G G: FULL NNNLO CORRECTIONS TO TOP LOOPS FOR NF=5
       FQCD0= HGGQCD(ASG,5)
       FQCD= HGGQCD2(ASG,5,MH,MT)
       XFAC= CDABS(CU(I)*(CXT+CXC)+CB(I)*CXB)**2*(FQCD-FQCD0)
       HJJ= HJJ+GF/(64d0*PI*SQR2)*MH**3*(ASG/PI)**2*XFAC
       FQCD0= AGGQCD(ASG,5)
       FQCD= AGGQCD2(ASG,5,MH,MT)
       XFAC= CJP(I)**2*(FQCD-FQCD0)
       HJJ= HJJ+GF/(64d0*PI*SQR2)*MH**3*(ASG/PI)**2*XFAC
       IF(NFGG.EQ.3)THEN
        HJJ= HJJ - DBB - DCC
       ELSEIF(NFGG.EQ.4)THEN
        HJJ= HJJ - DBB
        DCC= 0d0
       ELSE
        DCC= 0d0
        DBB= 0d0
       ENDIF
* End new
!       IF(NFGG.EQ.5)THEN
!        HJJ= HJJ+DBB+DCC
!        DBB= 0d0
!        DCC= 0d0
!       ELSEIF(NFGG.EQ.4)THEN
!        HJJ= HJJ+DCC
!        DCC= 0d0
!       ENDIF

*    Below CJ(P) becomes the REDUCED coupling to 2 gluons
*    as compared to the coupling of a SM Higgs with same mass
*    (CJ(P) does not contain strange quark loops)

       NFEXT= 3
       ASG= AS3
       FQCD= HGGQCD(ASG,NFEXT)
       SQCD= SGGQCD(ASG)
       CJ(I)= DSQRT(MAX(0d0,CJ(I)**2+CI(I)*SQCD/FQCD))/
     .        CDABS(CXT+CXC+CXB)
       CJP(I)= CJP(I)/CDABS(CXT+CXC+CXB)

*   h -> ee

       IF(MH.LE.2d0*MEL)THEN
        HEE= 0d0
       ELSE
        HEE= CD(I)**2*HFF(MH,(MEL/MH)**2)
     .      +CDP(I)**2*AFF(MH,(MEL/MH)**2)
       ENDIF

*   h -> mumu

       IF(MH.LE.2d0*MMU)THEN
        HMM= 0d0
       ELSE
        HMM= CD(I)**2*HFF(MH,(MMU/MH)**2)
     .      +CDP(I)**2*AFF(MH,(MMU/MH)**2)
       ENDIF

*   h -> tautau

       IF(MH.LE.2d0*MTAU)THEN
        HLL= 0d0
       ELSE
        HLL= CL(I)**2*HFF(MH,(MTAU/MH)**2)
     .      +CLP(I)**2*AFF(MH,(MTAU/MH)**2)
       ENDIF

*   h -> ss

       IF(MH.LE.2d0*MS)THEN
        HSS= 0d0
       ELSE
        IF(CD(I).NE.0d0)THEN
         RATCOUP= CU(I)/CD(I)
        ELSE
         RATCOUP=0d0
        ENDIF
        HS1= 3d0*HFF(MH,(MS/MH)**2)*CD(I)**2
     .             * TQCDH((MS/MH)**2)
        HS2= 3d0*HFF(MH,(RMS/MH)**2)*CD(I)**2
     .             * QCDH((RMS/MH)**2)
        IF(CDP(I).NE.0d0)THEN
         RATCOUP= CUP(I)/CDP(I)
        ELSE
         RATCOUP=0d0
        ENDIF
        HS1= HS1+3d0*AFF(MH,(MS/MH)**2)*CDP(I)**2
     .             * TQCDA((MS/MH)**2)
        HS2= HS2+3d0*AFF(MH,(RMS/MH)**2)*CDP(I)**2
     .             * QCDA((RMS/MH)**2)
        IF(HS2.LT.0d0) HS2=0d0
        RAT= 2d0*MS/MH
        HSS= QQINT(RAT,HS1,HS2)
       ENDIF

*   h -> cc

       IF(MH.LE.2d0*MC)THEN
        HCC= 0d0
       ELSE
        RATCOUP= 1d0
        HC1= 3d0*HFF(MH,(MC/MH)**2)*CU(I)**2
     .     * TQCDH((MC/MH)**2)
     .      +3d0*AFF(MH,(MC/MH)**2)*CUP(I)**2
     .     * TQCDA((MC/MH)**2)
        HC2= 3d0*HFF(MH,(RMC/MH)**2)*CU(I)**2
     .     * QCDH((RMC/MH)**2)
     .      +3d0*AFF(MH,(RMC/MH)**2)*CUP(I)**2
     .     * QCDA((RMC/MH)**2)
     .     + DCC
        IF(HC2.LT.0d0) HC2=0d0
        RAT= 2d0*MC/MH
        HCC= QQINT(RAT,HC1,HC2)
       ENDIF

*   h -> bb

       IF(MH.LE.10.6d0)THEN
        HBB= 0d0
       ELSE
        IF(CB(I).NE.0d0)THEN
         RATCOUP= CU(I)/CB(I)
        ELSE
         RATCOUP=0d0
        ENDIF
        HB1= 3d0*HFF(MH,(MBP/MH)**2)*CB(I)**2
     .            * TQCDH((MBP/MH)**2)
        HB2= 3d0*HFF(MH,(RMB/MH)**2)*CB(I)**2
     .            * QCDH((RMB/MH)**2)
        IF(CBP(I).NE.0d0)THEN
         RATCOUP= CUP(I)/CBP(I)
        ELSE
         RATCOUP=0d0
        ENDIF
        HB1= HB1+3d0*AFF(MH,(MBP/MH)**2)*CBP(I)**2
     .            * TQCDA((MBP/MH)**2)
        HB2= HB2+3d0*AFF(MH,(RMB/MH)**2)*CBP(I)**2
     .            * QCDA((RMB/MH)**2)
     .       + DBB
        IF(HB2.LT.0d0) HB2=0d0
        RAT= 2d0*MBP/MH
        HBB= QQINT(RAT,HB1,HB2)
       ENDIF

*   h -> tt

       IF (MH.LE.2d0*MT)THEN
        HTT= 0d0
       ELSE
        RATCOUP= 0d0
        HT1= 3d0*HFF(MH,(MT/MH)**2)*CU(I)**2
     .             * TQCDH((MT/MH)**2)
     .      +3d0*AFF(MH,(MT/MH)**2)*CUP(I)**2
     .             * TQCDA((MT/MH)**2)
        IF (MH.LE.2d0*RMT)THEN
         HT2= 0d0
        ELSE
         HT2= 3d0*HFF(MH,(RMT/MH)**2)*CU(I)**2
     .      * QCDH((RMT/MH)**2)
     .       +3d0*AFF(MH,(RMT/MH)**2)*CUP(I)**2
     .      * QCDA((RMT/MH)**2)
        ENDIF
        IF(HT2.LT.0d0) HT2=0d0
        RAT= 2d0*MT/MH
        HTT= QQINT(RAT,HT1,HT2)
       ENDIF

*   h -> WW

       IF(VFLAG.EQ.0)THEN
        DLD= 2d0
        DLU= 2d0
        XM1= 2d0*MW-DLD
        XM2= 2d0*MW+DLU
        IF(MH.LE.MW)THEN
         HWW= 0d0
        ELSEIF(MH.LE.XM1)THEN
         CWW= 3d0*GF**2*MW**4/16d0/PI**3
         HWW= CV(I)**2*HV((MW/MH)**2)*CWW*MH
        ELSEIF(MH.LT.XM2)THEN
         CWW= 3d0*GF**2*MW**4/16d0/PI**3
         XX(1)= XM1-1d0
         XX(2)= XM1
         XX(3)= XM2
         XX(4)= XM2+1d0
         YY(1)= HV((MW/XX(1))**2)*CWW*XX(1)
         YY(2)= HV((MW/XX(2))**2)*CWW*XX(2)
         YY(3)= HVV(XX(3),(MW/XX(3))**2)
         YY(4)= HVV(XX(4),(MW/XX(4))**2)
         HWW= CV(I)**2*FINT(MH,XX,YY)
        ELSE
         HWW= CV(I)**2*HVV(MH,(MW/MH)**2)
        ENDIF
       ELSE
        CALL HTOVV(MW,2.08856d0,MH,HTWW)
        HWW= 3d0/2d0*GF*MW**4/DSQRT(2d0)/PI/MH**3*HTWW*CV(I)**2
       ENDIF

*   h -> ZZ

       IF(VFLAG.EQ.0)THEN
        DLD= 2d0
        DLU= 2d0
        XM1= 2d0*MZ-DLD
        XM2= 2d0*MZ+DLU
        IF(MH.LE.MZ)THEN
         HZZ= 0d0
        ELSEIF(MH.LE.XM1)THEN
         CZZ= 3d0*GF**2*MZ**4/192d0/PI**3
     .      * (7d0-40d0/3d0*S2TW+160d0/9d0*S2TW**2)
         HZZ= CV(I)**2*HV((MZ/MH)**2)*CZZ*MH
        ELSEIF(MH.LT.XM2)THEN
         CZZ= 3d0*GF**2*MZ**4/192d0/PI**3
     .      * (7d0-40d0/3d0*S2TW+160d0/9d0*S2TW**2)
         XX(1)= XM1-1d0
         XX(2)= XM1
         XX(3)= XM2
         XX(4)= XM2+1d0
         YY(1)= HV((MZ/XX(1))**2)*CZZ*XX(1)
         YY(2)= HV((MZ/XX(2))**2)*CZZ*XX(2)
         YY(3)= HVV(XX(3),(MZ/XX(3))**2)/2d0
         YY(4)= HVV(XX(4),(MZ/XX(4))**2)/2d0
         HZZ= CV(I)**2*FINT(MH,XX,YY)
        ELSE
         HZZ= CV(I)**2*HVV(MH,(MZ/MH)**2)/2d0
        ENDIF
       ELSE
        CALL HTOVV(MZ,2.49581d0,MH,HTZZ)
        HZZ= 3d0/4d0*GF*MZ**4/DSQRT(2d0)/PI/MH**3*HTZZ*CV(I)**2
       ENDIF

*   h -> gamma gamma

       XFAC= CG(I)**2
       HGG= GF/(128d0*PI*SQR2)*MH**3*(ALEM0/PI)**2*XFAC
       XFAC= CGP(I)**2
       HGG= HGG+GF/(128d0*PI*SQR2)*MH**3*(ALEM0/PI)**2*XFAC

*    Below CG(P) becomes the REDUCDED coupling to 2 photons
*    as compared to the coupling of a SM Higgs with same mass

       CG(I)= CG(I)/CDABS(4d0/3d0*(CXT+CXC)+CXB/3d0+CXL+CXM+CXE+CXW)
       CGP(I)= CGP(I)/CDABS(4d0/3d0*(CXT+CXC)+CXB/3d0+CXL+CXM+CXE+CXW)
      
*  h -> Z gamma

      IF(MH.LE.MZ)THEN
       HZG= 0d0
      ELSE
       FT= -2d0*(1d0-8d0/3d0*S2TW)/DSQRT(S2TW*C2TW)*CU(I)
       FB= (-1d0+4d0/3d0*S2TW)/DSQRT(S2TW*C2TW)*CB(I)
       FCH1=4d0*MW/DSQRT(G2*S2TW*C2TW*MCH2(1))*COH0CH(I,1,1,1)
     .  *(-V(1,1,1)*V(1,1,1)-0.5d0*V(1,2,1)*V(1,2,1)+S2TW
     .    -V(1,1,2)*V(1,1,2)-0.5d0*V(1,2,2)*V(1,2,2)
     .    -U(1,1,1)*U(1,1,1)-0.5d0*U(1,2,1)*U(1,2,1)+S2TW
     .    -U(1,1,2)*U(1,1,2)-0.5d0*U(1,2,2)*U(1,2,2))
       FCH2=4d0*MW/DSQRT(G2*S2TW*C2TW*MCH2(2))*COH0CH(I,2,2,1)
     .  *(-V(2,1,1)*V(2,1,1)-0.5d0*V(2,2,1)*V(2,2,1)+S2TW
     .    -V(2,1,2)*V(2,1,2)-0.5d0*V(2,2,2)*V(2,2,2)
     .    -U(2,1,1)*U(2,1,1)-0.5d0*U(2,2,1)*U(2,2,1)+S2TW
     .    -U(2,1,2)*U(2,1,2)-0.5d0*U(2,2,2)*U(2,2,2))
       CLT= 4d0*(MT/MZ)**2*DCMPLX(1d0,-EPS)
       CLB= 4d0*(MBP/MZ)**2*DCMPLX(1d0,-EPS)
       CLC= 4d0*(MC/MZ)**2*DCMPLX(1d0,-EPS)
       CLW= 4d0*(MW/MZ)**2*DCMPLX(1d0,-EPS)
       CLH= 4d0*MHC/MZ**2*DCMPLX(1d0,-EPS)
       CLCH1= 4d0*MCH2(1)/MZ**2*DCMPLX(1d0,-EPS)
       CLCH2= 4d0*MCH2(2)/MZ**2*DCMPLX(1d0,-EPS)
       CXTZ= FT*(CI1(CTT,CLT) - CI2(CTT,CLT))
       CXBZ= FB*(CI1(CTB,CLB) - CI2(CTB,CLB))
       CXCZ= FT*(CI1(CTC,CLC) - CI2(CTC,CLC))
       CXWZ= -1d0/DSQRT(T2TW)*(4d0*(3d0-T2TW)*CI2(CTW,CLW)
     .     + ((1d0+2d0/CTW)*T2TW
     .     - (5d0+2d0/CTW))*CI1(CTW,CLW))*CV(I)
       CXHZ= (1d0-2d0*S2TW)/(DSQRT(S2TW*C2TW)*2d0*MHC)
     .     * CI1(CTHC,CLH)*GRH0HPHM(I,2,2)/(DSQRT(SQR2*GF))
       CXCH1Z=FCH1*(CI1(CTCH1,CLCH1) - CI2(CTCH1,CLCH1))
       CXCH2Z=FCH2*(CI1(CTCH2,CLCH2) - CI2(CTCH2,CLCH2))
       XFAC= CDABS(CXTZ+CXBZ+CXCZ+CXWZ+CXHZ+CXCH1Z+CXCH2Z)**2
       ACOUP= SQR2*GF*MZ**2*S2TW*C2TW/PI**2
       HZG= GF/(4d0*PI*SQR2)*MH**3*(ALEM0/PI)*ACOUP/16d0
     .    * XFAC*(1d0-(MZ/MH)**2)**3

       CXWZ= -1d0/DSQRT(T2TW)*(4d0*(3d0-T2TW)*CI2(CTW,CLW)
     .     + ((1d0+2d0/CTW)*T2TW
     .     - (5d0+2d0/CTW))*CI1(CTW,CLW))
       IF(XFAC.ge.EPS)THEN
        CZG(I)=dsqrt(XFAC/MAX(CDABS(CXTZ+CXBZ+CXCZ+CXWZ),EPS))
       ELSE
        CZG(I)=0d0
       ENDIF

       FT= -2d0*(1d0-8d0/3d0*S2TW)/DSQRT(S2TW*C2TW)*CUP(I)
       FB= (-1d0+4d0/3d0*S2TW)/DSQRT(S2TW*C2TW)*CDP(I)
       FCH1=4d0*MW/DSQRT(G2*S2TW*C2TW*MCH2(1))*COH0CH(I,1,1,2)
     .  *(-V(1,1,1)*V(1,1,1)-0.5d0*V(1,2,1)*V(1,2,1)+S2TW
     .    -V(1,1,2)*V(1,1,2)-0.5d0*V(1,2,2)*V(1,2,2)
     .    -U(1,1,1)*U(1,1,1)-0.5d0*U(1,2,1)*U(1,2,1)+S2TW
     .    -U(1,1,2)*U(1,1,2)-0.5d0*U(1,2,2)*U(1,2,2))
       FCH2=4d0*MW/DSQRT(G2*S2TW*C2TW*MCH2(2))*COH0CH(I,2,2,2)
     .  *(-V(2,1,1)*V(2,1,1)-0.5d0*V(2,2,1)*V(2,2,1)+S2TW
     .    -V(2,1,2)*V(2,1,2)-0.5d0*V(2,2,2)*V(2,2,2)
     .    -U(2,1,1)*U(2,1,1)-0.5d0*U(2,2,1)*U(2,2,1)+S2TW
     .    -U(2,1,2)*U(2,1,2)-0.5d0*U(2,2,2)*U(2,2,2))
       CLT= 4d0*(MT/MZ)**2*DCMPLX(1d0,-EPS)
       CLB= 4d0*(MBP/MZ)**2*DCMPLX(1d0,-EPS)
       CLC= 4d0*(MC/MZ)**2*DCMPLX(1d0,-EPS)
       CLCH1= 4d0*MCH2(1)/MZ**2*DCMPLX(1d0,-EPS)
       CLCH2= 4d0*MCH2(2)/MZ**2*DCMPLX(1d0,-EPS)
       CXTZ= FT*(-CI2(CTT,CLT))
       CXBZ= FB*(-CI2(CTB,CLB))
       CXCZ= FT*(-CI2(CTC,CLC))
       CXCH1Z=FCH1*(-CI2(CTCH1,CLCH1))
       CXCH2Z=FCH2*(-CI2(CTCH2,CLCH2))
       XFAC= CDABS(CXTZ+CXBZ+CXCZ+CXCH1Z+CXCH2Z)**2
       ACOUP= SQR2*GF*MZ**2*S2TW*C2TW/PI**2
       HZG= HZG+GF/(4d0*PI*SQR2)*MH**3*(ALEM0/PI)*ACOUP/16d0
     .    * XFAC*(1d0-(MZ/MH)**2)**3
      ENDIF

       CXTZ= FT*(CI1(CTT,CLT) - CI2(CTT,CLT))
       CXBZ= FB*(CI1(CTB,CLB) - CI2(CTB,CLB))
       CXCZ= FT*(CI1(CTC,CLC) - CI2(CTC,CLC))
       CXWZ= -1d0/DSQRT(T2TW)*(4d0*(3d0-T2TW)*CI2(CTW,CLW)
     .     + ((1d0+2d0/CTW)*T2TW
     .     - (5d0+2d0/CTW))*CI1(CTW,CLW))
       IF(XFAC.ge.EPS)THEN
        CZGP(I)=dsqrt(XFAC/MAX(CDABS(CXTZ+CXBZ+CXCZ+CXWZ),EPS))
       ELSE
        CZGP(I)=0d0
       ENDIF

*   h -> hh

       HTOT= 0d0

       DO J=1,10
        HHH(J)=0d0
       ENDDO

       IF(I.EQ.2)THEN

        IF(MH.LE.2d0*dsqrt(MH0(1)))THEN
         HHH(1)= 0d0
        ELSE
         RH= GH0H0H0(2,1,1)
         CH= BETA(MH0(1)/MH**2)
         HHH(1)= CH*RH**2/(32d0*PI*MH)
         HTOT= HTOT+HHH(1)
        ENDIF

       ENDIF

       IF(I.EQ.3)THEN

        IF(MH.LE.2d0*dsqrt(MH0(1)))THEN
         HHH(1)= 0d0
        ELSE
         RH= GH0H0H0(3,1,1)
         CH= BETA(MH0(1)/MH**2)
         HHH(1)= CH*RH**2/(32d0*PI*MH)
         HTOT= HTOT+HHH(1)
        ENDIF

        IF(MH.LE.dsqrt(MH0(1))+dsqrt(MH0(2)))THEN
         HHH(2)= 0d0
        ELSE
         RH= GH0H0H0(3,1,2)
         CH= LAMB(MH0(1)/MH**2,MH0(2)/MH**2)
         HHH(2)= CH*RH**2/(16d0*PI*MH)
         HTOT= HTOT+HHH(2)
        ENDIF

        IF(MH.LE.2d0*dsqrt(MH0(2)))THEN
         HHH(3)= 0d0
        ELSE
         RH= GH0H0H0(3,2,2)
         CH= BETA(MH0(2)/MH**2)
         HHH(3)= CH*RH**2/(32d0*PI*MH)
         HTOT= HTOT+HHH(3)
        ENDIF

       ENDIF

       IF(I.EQ.4)THEN

        IF(MH.LE.2d0*dsqrt(MH0(1)))THEN
         HHH(1)= 0d0
        ELSE
         RH= GH0H0H0(4,1,1)
         CH= BETA(MH0(1)/MH**2)
         HHH(1)= CH*RH**2/(32d0*PI*MH)
         HTOT= HTOT+HHH(1)
        ENDIF

        IF(MH.LE.dsqrt(MH0(1))+dsqrt(MH0(2)))THEN
         HHH(2)= 0d0
        ELSE
         RH= GH0H0H0(4,1,2)
         CH= LAMB(MH0(1)/MH**2,MH0(2)/MH**2)
         HHH(2)= CH*RH**2/(16d0*PI*MH)
         HTOT= HTOT+HHH(2)
        ENDIF

        IF(MH.LE.2d0*dsqrt(MH0(2)))THEN
         HHH(3)= 0d0
        ELSE
         RH= GH0H0H0(4,2,2)
         CH= BETA(MH0(2)/MH**2)
         HHH(3)= CH*RH**2/(32d0*PI*MH)
         HTOT= HTOT+HHH(3)
        ENDIF

        IF(MH.LE.dsqrt(MH0(1))+dsqrt(MH0(3)))THEN
         HHH(4)= 0d0
        ELSE
         RH= GH0H0H0(4,1,3)
         CH= LAMB(MH0(1)/MH**2,MH0(3)/MH**2)
         HHH(4)= CH*RH**2/(16d0*PI*MH)
         HTOT= HTOT+HHH(4)
        ENDIF

        IF(MH.LE.dsqrt(MH0(2))+dsqrt(MH0(3)))THEN
         HHH(5)= 0d0
        ELSE
         RH= GH0H0H0(4,2,3)
         CH= LAMB(MH0(2)/MH**2,MH0(3)/MH**2)
         HHH(5)= CH*RH**2/(16d0*PI*MH)
         HTOT= HTOT+HHH(5)
        ENDIF

        IF(MH.LE.2d0*dsqrt(MH0(3)))THEN
         HHH(6)= 0d0
        ELSE
         RH= GH0H0H0(4,3,3)
         CH= BETA(MH0(3)/MH**2)
         HHH(6)= CH*RH**2/(32d0*PI*MH)
         HTOT= HTOT+HHH(6)
        ENDIF
        
       ENDIF

       IF(I.EQ.5)THEN

        IF(MH.LE.2d0*dsqrt(MH0(1)))THEN
         HHH(1)= 0d0
        ELSE
         RH= GH0H0H0(5,1,1)
         CH= BETA(MH0(1)/MH**2)
         HHH(1)= CH*RH**2/(32d0*PI*MH)
         HTOT= HTOT+HHH(1)
        ENDIF

        IF(MH.LE.dsqrt(MH0(1))+dsqrt(MH0(2)))THEN
         HHH(2)= 0d0
        ELSE
         RH= GH0H0H0(5,1,2)
         CH= LAMB(MH0(1)/MH**2,MH0(2)/MH**2)
         HHH(2)= CH*RH**2/(16d0*PI*MH)
         HTOT= HTOT+HHH(2)
        ENDIF

        IF(MH.LE.2d0*dsqrt(MH0(2)))THEN
         HHH(3)= 0d0
        ELSE
         RH= GH0H0H0(5,2,2)
         CH= BETA(MH0(2)/MH**2)
         HHH(3)= CH*RH**2/(32d0*PI*MH)
         HTOT= HTOT+HHH(3)
        ENDIF

        IF(MH.LE.dsqrt(MH0(1))+dsqrt(MH0(3)))THEN
         HHH(4)= 0d0
        ELSE
         RH= GH0H0H0(5,1,3)
         CH= LAMB(MH0(1)/MH**2,MH0(3)/MH**2)
         HHH(4)= CH*RH**2/(16d0*PI*MH)
         HTOT= HTOT+HHH(4)
        ENDIF

        IF(MH.LE.dsqrt(MH0(2))+dsqrt(MH0(3)))THEN
         HHH(5)= 0d0
        ELSE
         RH= GH0H0H0(5,2,3)
         CH= LAMB(MH0(2)/MH**2,MH0(3)/MH**2)
         HHH(5)= CH*RH**2/(16d0*PI*MH)
         HTOT= HTOT+HHH(5)
        ENDIF

        IF(MH.LE.2d0*dsqrt(MH0(3)))THEN
         HHH(6)= 0d0
        ELSE
         RH= GH0H0H0(5,3,3)
         CH= BETA(MH0(3)/MH**2)
         HHH(6)= CH*RH**2/(32d0*PI*MH)
         HTOT= HTOT+HHH(6)
        ENDIF

        IF(MH.LE.dsqrt(MH0(1))+dsqrt(MH0(4)))THEN
         HHH(7)= 0d0
        ELSE
         RH= GH0H0H0(5,1,4)
         CH= LAMB(MH0(1)/MH**2,MH0(4)/MH**2)
         HHH(7)= CH*RH**2/(16d0*PI*MH)
         HTOT= HTOT+HHH(7)
        ENDIF

        IF(MH.LE.dsqrt(MH0(2))+dsqrt(MH0(4)))THEN
         HHH(8)= 0d0
        ELSE
         RH= GH0H0H0(5,2,4)
         CH= LAMB(MH0(2)/MH**2,MH0(4)/MH**2)
         HHH(8)= CH*RH**2/(16d0*PI*MH)
         HTOT= HTOT+HHH(8)
        ENDIF

        IF(MH.LE.dsqrt(MH0(3))+dsqrt(MH0(4)))THEN
         HHH(9)= 0d0
        ELSE
         RH= GH0H0H0(5,3,4)
         CH= LAMB(MH0(3)/MH**2,MH0(4)/MH**2)
         HHH(9)= CH*RH**2/(16d0*PI*MH)
         HTOT= HTOT+HHH(9)
        ENDIF

        IF(MH.LE.2d0*dsqrt(MH0(4)))THEN
         HHH(10)= 0d0
        ELSE
         RH= GH0H0H0(5,4,4)
         CH= BETA(MH0(4)/MH**2)
         HHH(10)= CH*RH**2/(32d0*PI*MH)
         HTOT= HTOT+HHH(10)
        ENDIF

       ENDIF

*   h -> h+h-

       IF(MH.LE.2d0*dsqrt(MHC))THEN
        HHCHC= 0d0
       ELSE
        RH= GRH0HPHM(I,2,2)
        CH= BETA(MHC/MH**2)
        HHCHC= CH*RH**2/(16d0*PI*MH)
        HTOT= HTOT+HHCHC
       ENDIF

*   h -> h'Z

       DO J=1,4
        IF(MH.LE.dsqrt(MH0(J))+MZ)THEN
         HAZ(J)= 0d0
        ELSE
         RH= XH(I,4)*(XH(J,1)*COSB-XH(J,2)*SINB)
     .      -XH(J,4)*(XH(I,1)*COSB-XH(I,2)*SINB)
         CH= LAMB(MH0(J)/MH**2,(MZ/MH)**2)
     .     * LAMB((MH/MZ)**2,MH0(J)/MZ**2)**2
         HAZ(J)= GF/8d0/SQR2/PI*MZ**4/MH*CH*RH**2
         HTOT= HTOT+HAZ(J)
        ENDIF
       ENDDO

*   h -> h+W-

       IF(MH.LE.dsqrt(MHC)+MW)THEN
        HHCW= 0d0
       ELSE
        RH= XH(I,1)*COSB-XH(I,2)*SINB
        CH= LAMB(MHC/MH**2,(MW/MH)**2)
     .    * LAMB((MH/MW)**2,MHC/MW**2)**2
        HHCW= GF/8d0/SQR2/PI*MW**4/MH*CH*(RH**2+XH(I,4)**2)
        HTOT= HTOT+2d0*HHCW
       ENDIF

*   h -> neutralinos/charginos

       STOT= 0d0

       DO J=1,5
       DO M=J,5
        IF(MH.LE.dsqrt(MNEU(J))+dsqrt(MNEU(M)))THEN
         HNEU(J,M)= 0d0
        ELSE
         HNEU(J,M)= 1d0/(16d0*PI)*LAMB(MNEU(J)/MH**2,MNEU(M)/MH**2)/MH
     .   *((MH**2-(dsqrt(MNEU(J))+dsqrt(MNEU(M)))**2)
     .                               *COH0NEU(I,J,M,1)**2
     .    +(MH**2-(dsqrt(MNEU(J))-dsqrt(MNEU(M)))**2)
     .                               *COH0NEU(I,J,M,2)**2)
         IF(J.NE.M)HNEU(J,M)= 2d0*HNEU(J,M)
         STOT= STOT+HNEU(J,M)
        ENDIF
        HNEU(M,J)= HNEU(J,M)
       ENDDO
       ENDDO

       IF(MH.LE.2d0*dsqrt(MCH2(1)))THEN
        HCHA(1)= 0d0
       ELSE
        HCHA(1)= 1d0/(8d0*PI)*MH*LAMB(MCH2(1)/MH**2,MCH2(1)/MH**2)
     .  *(LAMB(MCH2(1)/MH**2,MCH2(1)/MH**2)**2*COH0CH(I,1,1,1)**2
     .                                        +COH0CH(I,1,1,2)**2)
        STOT= STOT+HCHA(1)
       ENDIF

       IF(MH.LE.dsqrt(MCH2(1))+dsqrt(MCH2(2)))THEN
        HCHA(2)= 0d0
       ELSE
        HCHA(2)= 1d0/(16d0*PI)*MH
     .      * LAMB(MCH2(1)/MH**2,MCH2(2)/MH**2)
     .      * ((COH0CH(I,1,2,1)**2+COH0CH(I,1,2,2)**2
     .         +COH0CH(I,2,1,1)**2+COH0CH(I,2,1,2)**2)
     .      * (1d0-MCH2(1)/MH**2-MCH2(2)/MH**2)
     .      - 4d0*(COH0CH(I,1,2,1)*COH0CH(I,2,1,1)
     .             -COH0CH(I,1,2,2)*COH0CH(I,2,1,2))
     .      * dsqrt(MCH2(1)*MCH2(2))/MH**2)
        STOT= STOT+2d0*HCHA(2)
       ENDIF

       IF(MH.LE.2d0*dsqrt(MCH2(2)))THEN
        HCHA(3)= 0d0
       ELSE
        HCHA(3)= 1d0/(8d0*PI)*MH*LAMB(MCH2(2)/MH**2,MCH2(2)/MH**2)
     .  *(LAMB(MCH2(2)/MH**2,MCH2(2)/MH**2)**2*COH0CH(I,2,2,1)**2
     .                                        +COH0CH(I,2,2,2)**2)
        STOT= STOT+HCHA(3)
       ENDIF

*   h -> squarks
*   UL, UR, DL, DR are the first two families

       IF(MH.LE.2d0*dsqrt(MSU2P(1)))THEN
        HSQ(1)= 0d0
       ELSE
        RH= GRHSUSU(I,1,1)
        CH= BETA(MSU2P(1)/MH**2)
        HSQ(1)= 3d0*CH*RH**2/(16d0*PI*MH)
        STOT= STOT+2d0*HSQ(1)
       ENDIF

       IF(MH.LE.2d0*dsqrt(MSU2P(2)))THEN
        HSQ(2)= 0d0
       ELSE
        RH= GRHSUSU(I,2,2)
        CH= BETA(MSU2P(2)/MH**2)
        HSQ(2)= 3d0*CH*RH**2/(16d0*PI*MH)
        STOT= STOT+2d0*HSQ(2)
       ENDIF

       IF(MH.LE.2d0*dsqrt(MSD2P(1)))THEN
        HSQ(3)= 0d0
       ELSE
        RH= GRHSDSD(I,1,1)
        CH= BETA(MSD2P(1)/MH**2)
        HSQ(3)= 3d0*CH*RH**2/(16d0*PI*MH)
        STOT= STOT+2d0*HSQ(3)
       ENDIF

       IF(MH.LE.2d0*dsqrt(MSD2P(2)))THEN
        HSQ(4)= 0d0
       ELSE
        RH= GRHSDSD(I,2,2)
        CH= BETA(MSD2P(2)/MH**2)
        HSQ(4)= 3d0*CH*RH**2/(16d0*PI*MH)
        STOT= STOT+2d0*HSQ(4)
       ENDIF

       IF(MH.LE.2d0*dsqrt(MST2P(1)))THEN
        HSQ(5)= 0d0
       ELSE
        RH= dsqrt(GRHSTST(I,1,1)**2+GIHSTST(I,1,1)**2)
        CH= BETA(MST2P(1)/MH**2)
        HSQ(5)= 3d0*CH*RH**2/(16d0*PI*MH)
        STOT= STOT+HSQ(5)
       ENDIF

       IF(MH.LE.2d0*dsqrt(MST2P(2)))THEN
        HSQ(6)= 0d0
       ELSE
        RH= dsqrt(GRHSTST(I,2,2)**2+GIHSTST(I,2,2)**2)
        CH= BETA(MST2P(2)/MH**2)
        HSQ(6)= 3d0*CH*RH**2/(16d0*PI*MH)
        STOT= STOT+HSQ(6)
       ENDIF

       IF(MH.LE.dsqrt(MST2P(1))+dsqrt(MST2P(2)))THEN
        HSQ(7)= 0d0
       ELSE
        RH= dsqrt(GRHSTST(I,1,2)**2+GIHSTST(I,1,2)**2)
        CH= LAMB(MST2P(1)/MH**2,MST2P(2)/MH**2)
        HSQ(7)= 3d0*CH*RH**2/(16d0*PI*MH)
        STOT= STOT+2d0*HSQ(7)
       ENDIF

       IF(MH.LE.2d0*dsqrt(MSB2P(1)))THEN
        HSQ(8)= 0d0
       ELSE
        RH= dsqrt(GRHSBSB(I,1,1)**2+GIHSBSB(I,1,1)**2)
        CH= BETA(MSB2P(1)/MH**2)
        HSQ(8)= 3d0*CH*RH**2/(16d0*PI*MH)
        STOT= STOT+HSQ(8)
       ENDIF

       IF(MH.LE.2d0*dsqrt(MSB2P(2)))THEN
        HSQ(9)= 0d0
       ELSE
        RH= dsqrt(GRHSBSB(I,2,2)**2+GIHSBSB(I,2,2)**2)
        CH= BETA(MSB2P(2)/MH**2)
        HSQ(9)= 3d0*CH*RH**2/(16d0*PI*MH)
        STOT= STOT+HSQ(9)
       ENDIF

       IF(MH.LE.dsqrt(MSB2P(1))+dsqrt(MSB2P(2)))THEN
        HSQ(10)= 0d0
       ELSE
        RH= dsqrt(GRHSTST(I,1,2)**2+GIHSTST(I,1,2)**2)
        CH= LAMB(MSB2P(1)/MH**2,MSB2P(2)/MH**2)
        HSQ(10)= 3d0*CH*RH**2/(16d0*PI*MH)
        STOT= STOT+2d0*HSQ(10)
       ENDIF

*   h -> sleptons
*   LL, LR, NL are the first two families

       IF(MH.LE.2d0*dsqrt(MSE2(1)))THEN
        HSL(1)= 0d0
       ELSE
        RH= GRHSESE(I,1,1)
        CH= BETA(MSE2(1)/MH**2)
        HSL(1)= CH*RH**2/(16d0*PI*MH)
        STOT= STOT+2d0*HSL(1)
       ENDIF

       IF(MH.LE.2d0*dsqrt(MSE2(2)))THEN
        HSL(2)= 0d0
       ELSE
        RH= GRHSESE(I,2,2)
        CH= BETA(MSE2(2)/MH**2)
        HSL(2)= CH*RH**2/(16d0*PI*MH)
        STOT= STOT+2d0*HSL(2)
       ENDIF

       IF(MH.LE.2d0*dsqrt(MSNE2))THEN
        HSL(3)= 0d0
       ELSE
        RH= GRHSNSN(I)
        CH= BETA(MSNE2/MH**2)
        HSL(3)= CH*RH**2/(16d0*PI*MH)
        STOT= STOT+2d0*HSL(3)
       ENDIF

       IF(MH.LE.2d0*dsqrt(MSL2(1)))THEN
        HSL(4)= 0d0
       ELSE
        RH= dsqrt(GRHSLSL(I,1,1)**2+GIHSLSL(I,1,1)**2)
        CH= BETA(MSL2(1)/MH**2)
        HSL(4)= CH*RH**2/(16d0*PI*MH)
        STOT= STOT+HSL(4)
       ENDIF

       IF(MH.LE.2d0*dsqrt(MSL2(2)))THEN
        HSL(5)= 0d0
       ELSE
        RH= dsqrt(GRHSLSL(I,2,2)**2+GIHSLSL(I,2,2)**2)
        CH= BETA(MSL2(2)/MH**2)
        HSL(5)= CH*RH**2/(16d0*PI*MH)
        STOT= STOT+HSL(5)
       ENDIF

       IF(MH.LE.dsqrt(MSL2(1))+dsqrt(MSL2(2)))THEN
        HSL(6)= 0d0
       ELSE
        RH= dsqrt(GRHSLSL(I,1,2)**2+GIHSLSL(I,1,2)**2)
        CH= LAMB(MSL2(1)/MH**2,MSL2(2)/MH**2)
        HSL(6)= CH*RH**2/(16d0*PI*MH)
        STOT= STOT+2d0*HSL(6)
       ENDIF

       IF(MH.LE.2d0*dsqrt(MSNT2))THEN
        HSL(7)= 0d0
       ELSE
        RH= GRHSNSN(I)
        CH= BETA(MSNT2/MH**2)
        HSL(7)= CH*RH**2/(16d0*PI*MH)
        STOT= STOT+HSL(7)
       ENDIF

*  Branching ratios

       WIDTH(I)= HJJ+HEE+HMM+HLL+HSS+HCC+HBB+HTT+HWW+HZZ+HGG+HZG
     .          +HTOT+STOT
       BRJJ(I)= (HJJ+HSS)/WIDTH(I)
       BREE(I)= HEE/WIDTH(I)
       BRMM(I)= HMM/WIDTH(I)
       BRLL(I)= HLL/WIDTH(I)
       BRCC(I)= HCC/WIDTH(I)
       BRBB(I)= HBB/WIDTH(I)
       BRTT(I)= HTT/WIDTH(I)
       BRWW(I)= HWW/WIDTH(I)
       BRZZ(I)= HZZ/WIDTH(I)
       BRGG(I)= HGG/WIDTH(I)
       BRZG(I)= HZG/WIDTH(I)
       BRHIGGS(I)= HTOT/WIDTH(I)
       BRSUSY(I)= STOT/WIDTH(I)
       DO J=1,10
        BRHHH(I,J)=HHH(J)/WIDTH(I)
       ENDDO
       BRHCHC(I)= HHCHC/WIDTH(I)
       DO J=1,4
        BRHAZ(I,J)=HAZ(J)/WIDTH(I)
       ENDDO
       BRHCW(I)= HHCW/WIDTH(I)
       DO J=1,5
       DO M=1,5
        BRNEU(I,J,M)=HNEU(J,M)/WIDTH(I)
       ENDDO
       ENDDO
       DO J=1,3
        BRCHA(I,J)=HCHA(J)/WIDTH(I)
       ENDDO
       DO J=1,10
        BRHSQ(I,J)=HSQ(J)/WIDTH(I)
       ENDDO
       DO J=1,7
        BRHSL(I,J)=HSL(J)/WIDTH(I)
       ENDDO

         ENDIF

      ENDDO

* Charged Higgs boson

      MH=dsqrt(MHC)
      ASH= ALPHAS(MH,2)

*  Running quark masses

      RMS= RUNM(MH,3)
      RMC= RUNM(MH,4)

*   Running bottom/top masses at MH:

       IF(MH.GE.MT)THEN
        RMT= RMTTOP
     .   *(1d0+7d0/(4d0*PI)*ASMT*DLOG(MH**2/MT**2))
     .   **(-4d0/7d0)
       ELSE
         RMT= RMTTOP
     .   *(1d0+23d0/(12d0*PI)*ASMT*DLOG(MH**2/MT**2))
     .   **(-12d0/23d0)

       ENDIF

       RMB=RUNMB(MH)

*  Partial widths

*   h+ -> munu

      IF(MH.LE.MMU)THEN
       HMN= 0d0
      ELSE
       HMN= CFF(MH,TANB,(MMU/MH)**2,0d0)
      ENDIF

*   h+ -> taunu

      IF(MH.LE.MTAU)THEN
       HLN= 0d0
      ELSE
       HLN= CFF(MH,TANB,(MTAU/MH)**2,0d0)
      ENDIF
      HLN=HLN/(1d0+DELML)**2

*   h+ -> su

      IF(MH.LE.MS+EPS)THEN
       HSU= 0d0
      ELSE
       HSU1= 3d0*VUS**2*CQCDM(MH,TANB,(MS/MH)**2,EPS)
       HSU2= 3d0*VUS**2*CQCD(MH,TANB,(RMS/MH)**2,EPS)
       IF(HSU2.LT.0d0) HSU2= 0d0
       RAT= MS/MH
       HSU= QQINT(RAT,HSU1,HSU2)
      ENDIF

*   h+ -> cs

      IF(MH.LE.MS+MC)THEN
       HSC= 0d0
      ELSE
       HSC1= 3d0*CQCDM(MH,TANB,(MS/MH)**2,(MC/MH)**2)
       HSC2= 3d0*CQCD(MH,TANB,(RMS/MH)**2,(RMC/MH)**2)
       IF(HSC2.LT.0d0) HSC2= 0d0
       RAT= (MS+MC)/MH
       HSC= QQINT(RAT,HSC1,HSC2)
      ENDIF

*   h+ -> cb

      IF(MH.LE.MBP+MC)THEN
       HBC= 0d0
      ELSE
       HBC1= 3d0*VCB**2
     .  *CQCDM(MH,TANB,(MBP/MH)**2,(MC/MH)**2)
       HBC2= 3d0*VCB**2
     .  *CQCD(MH,TANB,(RMB/MH)**2,(RMC/MH)**2)
       IF(HBC2.LT.0d0) HBC2= 0d0
       RAT= (MBP+MC)/MH
       HBC= QQINT(RAT,HBC1,HBC2)
      ENDIF
      HBC=HBC/(1d0+DELMB)**2

*   h+ -> bu

      IF(MH.LE.MBP+EPS)THEN
       HBU= 0d0
      ELSE
       HBU1= 3d0*VUB**2*CQCDM(MH,TANB,(MBP/MH)**2,EPS)
       HBU2= 3d0*VUB**2*CQCD(MH,TANB,(RMB/MH)**2,EPS)
       IF(HBU2.LT.0d0) HBU2= 0d0
       RAT= MBP/MH
       HBU= QQINT(RAT,HBU1,HBU2)
      ENDIF
      HBU=HBU/(1d0+DELMB)**2

*   h+ -> tb :

      IF(MH.LE.MT+MBP)THEN
       HBT= 0d0
      ELSE
       HBT1= 3d0*CQCDM(MH,TANB,(MBP/MH)**2,(MT/MH)**2)
       IF(MH.LE.RMT+RMB)THEN
        HBT2= 0d0
       ELSE
        HBT2= 3d0*CQCD(MH,TANB,(RMB/MH)**2,(RMT/MH)**2)
       ENDIF
       IF(HBT2.LT.0d0) HBT2= 0d0
       RAT= (MBP+MT)/MH
       HBT= QQINT(RAT,HBT1,HBT2)
      ENDIF
      HBT=HBT/(1d0+DELMB)**2

*   h+ -> hW

      HTOT= 0d0

      DO I=1,5
       IF(MH.LE.MW+dsqrt(MH0(I)))THEN
        HCWH(I)= 0d0
       ELSE
        RH= XH(I,1)*COSB-XH(I,2)*SINB
        CH= LAMB(MH0(I)/MH**2,(MW/MH)**2)
     .   *LAMB((MH/MW)**2,MH0(I)/MW**2)**2
        HCWH(I)= GF/8d0/SQR2/PI*MW**4/MH*CH*RH**2
        RH= XH(I,4)
        HCWH(I)= HCWH(I)+GF/8d0/SQR2/PI*MW**4/MH*CH*RH**2
       HTOT= HTOT+HCWH(I)
       ENDIF
      ENDDO

*   h+ -> charginos+neutralinos

      STOT= 0d0

      DO I=1,5
      DO J=1,2
       IF (MH.LE.dsqrt(MNEU(I))+dsqrt(MCH2(J)))THEN
        HCNC(I,J)= 0d0
       ELSE
        HCNC(I,J)= 1d0/(16d0*PI)*MH
     .       * LAMB(MNEU(I)/MH**2,MCH2(J)/MH**2)
     .       * ((COHPNEUCHM(2,I,J,1)**2+COHMNEUCHP(2,I,J,1)**2
     .          +COHPNEUCHM(2,I,J,2)**2+COHMNEUCHP(2,I,J,2)**2)
     .       * (1d0-MNEU(I)/MH**2-MCH2(J)/MH**2)
     .       - 4d0*(COHPNEUCHM(2,I,J,1)*COHMNEUCHP(2,I,J,1)
     .             -COHPNEUCHM(2,I,J,2)*COHMNEUCHP(2,I,J,2))
     .       * dsqrt(MNEU(I)*MCH2(J))/MH**2)

        STOT= STOT+HCNC(I,J)
       ENDIF
      ENDDO
      ENDDO

*   h+ -> squarks

       IF(MH.LE.dsqrt(MSU2P(1))+dsqrt(MSD2P(1)))THEN
        HCSQ(1)= 0d0
       ELSE
        RH= gRHCSUSD(1,1)
        CH= LAMB(MSU2P(1)/MH**2,MSD2P(1)/MH**2)
        HCSQ(1)= 3d0*CH*RH**2/(16d0*PI*MH)
        STOT= STOT+2d0*HCSQ(1)
       ENDIF

      DO I=1,2
      DO J=1,2
       IF(MH.LE.dsqrt(MST2P(I))+dsqrt(MSB2P(J)))THEN
        HCSQ(2*I+J-1)= 0d0
       ELSE
        RH=dsqrt(gRHCSTSB(I,J)**2+gIHCSTSB(I,J)**2)
        CH= LAMB(MST2P(I)/MH**2,MSB2P(J)/MH**2)
        HCSQ(2*I+J-1)= 3d0*CH*RH**2/(16d0*PI*MH)
        STOT= STOT+HCSQ(2*I+J-1)
       ENDIF
      ENDDO
      ENDDO

*   h+ -> sleptons

       IF(MH.LE.dsqrt(MSE2(1))+dsqrt(MSNE2))THEN
        HCSL(1)= 0d0
       ELSE
        RH= gRHCSNSE(1)
        CH= LAMB(MSE2(1)/MH**2,MSNE2/MH**2)
        HCSL(1)= CH*RH**2/(16d0*PI*MH)
        STOT= STOT+2d0*HCSL(1)
       ENDIF

      DO I=1,2
       IF(MH.LE.dsqrt(MSL2(I))+dsqrt(MSNT2))THEN
        HCSL(I+1)= 0d0
       ELSE
        RH= dsqrt(gRHCSNSL(I)**2+gIHCSNSL(I)**2)
        CH= LAMB(MSL2(1)/MH**2,MSNT2/MH**2)
        HCSL(I+1)= CH*RH**2/(16d0*PI*MH)
        STOT= STOT+HCSL(I+1)
       ENDIF
      ENDDO

*  Branching ratios

      HCWIDTH= HLN+HMN+HSU+HBU+HSC+HBC+HBT+HTOT+STOT

      HCBRM= HMN/HCWIDTH
      HCBRL= HLN/HCWIDTH
      HCBRSU= HSU/HCWIDTH
      HCBRBU= HBU/HCWIDTH
      HCBRSC= HSC/HCWIDTH
      HCBRBC= HBC/HCWIDTH
      HCBRBT= HBT/HCWIDTH
      HCBRWHT= HTOT/HCWIDTH
      HCBRSUSY= STOT/HCWIDTH
      DO I=1,5
       HCBRWH(I)= HCWH(I)/HCWIDTH
      ENDDO
      DO I=1,5
       DO J=1,2
        HCBRNC(I,J)= HCNC(I,J)/HCWIDTH
       ENDDO
      ENDDO
      DO I=1,5
       HCBRSQ(I)= HCSQ(I)/HCWIDTH
      ENDDO
      DO I=1,3
       HCBRSL(I)= HCSL(I)/HCWIDTH
      ENDDO


      RETURN
      END
