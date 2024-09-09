      PROGRAM MAIN

*  Program to compute the NMSSM Higgs masses, couplings and
*  Branching ratios, with experimental and theoretical constraints
*  in the CP-violating version.
*
*   On input:
*
*      REALP(1)  = Re(lambda)
*      REALP(2)  = Re(kappa)
*      REALP(3)  = Re(M1)
*      REALP(4)  = Re(M2)
*      REALP(5)  = Re(M3)
*      REALP(6)  = Re(AU3)
*      REALP(7)  = Re(AD3)
*      REALP(8)  = Re(AE3)
*      REALP(9)  = Re(XIF)
*      REALP(10) = Re(XIS)
*      REALP(11) = Re(MUP)
*      REALP(12) = Re(MSP)
*      REALP(13) = Re(M3H) 
*      REALP(14) = Re(MUEFF)
*
*      IMAGP(1)  = Im(lambda)
*      IMAGP(2)  = Im(kappa)
*      IMAGP(3)  = Im(M1)
*      IMAGP(4)  = Im(M2)
*      IMAGP(5)  = Im(M3)
*      IMAGP(6)  = Im(AU3)
*      IMAGP(7)  = Im(AD3)
*      IMAGP(8)  = Im(AE3)
*      IMAGP(9)  = Im(XIF)
*      IMAGP(10) = Im(XIS)
*      IMAGP(11) = Im(MUP)
*      IMAGP(12) = Im(MSP)
*      IMAGP(13) = Im(M3H)
*    [ IMAGP(14) = Im(MUEFF) -> aligned with lambda
*                = REALP(14)*IMAGP(1)/REALP(1)                      ] 
*
*    [ PAR(1)  = |lambda| 
*                -> Sqrt( REALP(1)**2 + IMAGP(1)**2 )               ]
*    [ PAR(2)  = sgn(Re(kappa)).|kappa| 
*                -> sgn(REALP(2)).Sqrt( REALP(2)**2 + IMAGP(2)**2 ) ]
*      PAR(3)  = tan(beta)
*    [ PAR(4)  = sgn(Re(mueff)).|mueff| -> aligned with lambda
*                = REALP(14)*PAR(1)/REALP(1)                        ] 
*      PAR(5)  = Re(Alambda)     (if (MA,XIF) is not an input)
*      PAR(6)  = Re(Akappa)      (if (MP,XIS) is not an input)
*      PAR(7)  = mQ3**2
*      PAR(8)  = mU3**2
*      PAR(9)  = mD3**2
*      PAR(10) = mL3**2
*      PAR(11) = mE3**2
*      PAR(12) = AU3
*      PAR(13) = AD3
*      PAR(14) = AE3
*      PAR(15) = mQ2**2
*      PAR(16) = mU2**2
*      PAR(17) = mD2**2
*      PAR(18) = mL2**2
*      PAR(19) = mE2**2
*    [ PAR(20) = sgn(Re(M1)).|M1|
*                -> sgn(REALP(3)).Sqrt( REALP(3)**2 + IMAGP(3)**2 ) ]
*    [ PAR(21) = sgn(Re(M2)).|M2|
*                -> sgn(REALP(4)).Sqrt( REALP(4)**2 + IMAGP(4)**2 ) ]
*    [ PAR(22) = sgn(Re(M3)).|M3|
*                -> sgn(REALP(5)).Sqrt( REALP(5)**2 + IMAGP(5)**2 ) ]
*      PAR(23) = MA (diagonal doublet CP-odd mass matrix element at tree-level)
*      PAR(24) = MP (diagonal singlet CP-odd mass matrix element at tree-level)
*      PAR(25) = AE2
*
*   Phases in [-Pi/2,Pi/2]
*    [ MYPHASES(1)  = Phi_lambda 
*                -> arctan( IMAGP(1) / REALP(1) )                   ]
*    [ MYPHASES(2)  = Phi_kappa
*                -> arctan( IMAGP(2) / REALP(2) )                   ]
*    [ MYPHASES(3)  = Phi_M1
*                -> arctan( IMAGP(3) / REALP(3) )                   ]
*    [ MYPHASES(4)  = Phi_M2
*                -> arctan( IMAGP(4) / REALP(4) )                   ]
*    [ MYPHASES(5)  = Phi_M3
*                -> arctan( IMAGP(5) / REALP(5) )                   ]
*    [ MYPHASES(6)  = Phi_AU3
*                -> arctan( IMAGP(6) / REALP(6) )                   ]
*    [ MYPHASES(7)  = Phi_AD3
*                -> arctan( IMAGP(7) / REALP(7) )                   ]
*    [ MYPHASES(8)  = Phi_AE3
*                -> arctan( IMAGP(8) / REALP(8) )                   ]
*    [ MYPHASES(9)  = Phi_AU2 (unused)
*                -> 0                                               ]
*    [ MYPHASES(10) = Phi_AD2 (unused)
*                -> 0                                               ]
*    [ MYPHASES(11) = Phi_AE2 (unused)
*                -> 0                                               ]
*    [ MYPHASES(12) = Phi_XIS
*                -> arctan( IMAGP(10) / REALP(10) )                 ]
*    [ MYPHASES(13) = Phi_MSP
*                -> arctan( IMAGP(12) / REALP(12) )                 ]
*    [ MYPHASES(14) = Phi_XIF
*                -> arctan( IMAGP(9) / REALP(9) )                   ]
*    [ MYPHASES(15) = Phi_MUP
*                -> arctan( IMAGP(11) / REALP(11) )                 ]
*    [ MYPHASES(16) = Phi_M3H
*                -> arctan( IMAGP(13) / REALP(13) )                 ]
*
*      All these parameters are assumed to be defined in DRbar
*      at the scale Q2, except for tan(beta) defined at MZ.
*      Q2 is either defined by the user in the input file or
*      computed as Q2 = (2*mQ2+mU2+mD2)/4
*
*   On output:
*
*      MH0(1-5): Higgs masses (squared, ordered)
*
*      XH0(1-5,1-5): Mixing matrix in the basis H1_R=Hu_R, H2_R=Hd_R, S_R, A, S_I
*                    where the Goldstone mode has been rotated away
*
*      MHC: Charged Higgs mass (squared)
*
*      CU,CD,CV,CJ,CG(i),CZG(i) Reduced scalar couplings of hi (i=1..5) to up type 
*                               fermions, down type fermions, gauge bosons, gluons, 
*                               photons and (Z+photon)
*      CUP,CDP,CJP,CGP(i),CZGP(i) idem for pseudoscalar component
*      CB(i)                    Reduced couplings of hi (i=1..5) to b-quarks 
*                               including DELMB corrections
*      CBP(i)                   idem for pseudoscalar component
*
*      WIDTH(i) Total decay width of hi (i=1..5)
*               with the following branching ratios:
*      BRJJ(i)   hi (i=1..5)  -> hadrons
*      BREE(i)        "       -> e+ e-
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
*      MNEU(i)     Mass of neutralino chi_i (i=1,5, squared, ordered in mass)
*      NEU(i,j,k)  chi_i components of bino, wino, higgsino u&d, singlino
*                  (i,j=1..5), real part (k=1) and imaginary part (k=2)
*
*      MCH2(i)           Chargino masses (squared)
*      U(i,j,k),V(i,j,k) Chargino mixing matrices,
*                        real part (k=1) and imaginary part (k=2)
*
*      Sfermion masses: 
*        MSU2(i)     S-up/charm, masses squared (i=1->L,2->R), DRbar
*        MSU2P(i)    idem, including QCD+Yuk. corrections
*        MSD2(i)     S-down/strange, masses squared (i=1->L,2->R), DRbar
*        MSD2P(i)    idem, including QCD+Yuk. corrections
*        MSE2(i)     S-electron, masses squared (i=1->L,2->R), DRbar
*        MSMU2(i)    S-muon, masses squared ordered (i=1,2), DRbar
*         UMU(i,j,k) S-muon rotation matrix (i->mass, j=1->L, 2->R)*
*        MSNE2       Sneutrino 1st gen., mass squared
*        MST2(i)     S-top, masses squared ordered (i=1->L,2->R), DRbar
*        MST2P(i)    idem, including QCD+Yuk. corrections
*         UT(i,j,k)  S-top rotation matrix (i->mass, j=1->L, 2->R)*
*        MSB2(i)     S-bottom, masses squared ordered (i=1->L,2->R), DRbar
*        MSB2P(i)    idem, including QCD+Yuk. corrections
*         UB(i,j,k)  S-bottom rotation matrix (i->mass, j=1->L, 2->R)*
*        MSL2(i)     S-tau, masses squared ordered (i=1,2), DRbar
*         UL(i,j,k)  S-tau rotation matrix (i->mass, j=1->L, 2->R)*
*        MSNT2       Sneutrino 3rd gen., mass squared
*                       * [real part (k=1) and imaginary part (k=2)]
*
*      Minimisation conditions:
*        MHuS    MH1^2 soft squared mass (minimisation wrt. h1_R)
*        MHdS    MH2^2 soft squared mass (minimisation wrt. h2_R)
*        MSS     MS^2 soft squared mass (minimisation wrt. R_R)
*        IAl     Im(Alambda) (minimisation wrt. h1_I/h2_I)
*        IAk     Im(Akappa) (minimisation wrt. S_I, if kappa =/= 0)
*        IXIS    Im(XIS) (minimisation wrt. h1_I, if kappa = 0)
*
*  ERRORS: IFAIL = 0..10
*
*  IFAIL = 0         OK
*          1         m_h1^2 < 0
*          4         m_h+^2 < 0
*          8         m_sfermion^2 < 0
*          9         mu = 0 | k=0 & (Ak|XiS)=/=0
*          10        Violation of phenomenological constraint(s)
*
*  Phenomenological constraints:
*
*      PROB(I)  = 0, I = 1..82: OK
*
*      PROB(1) =/= 0   chargino too light
*      PROB(2) =/= 0   excluded by Z -> neutralinos
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
*      PROB(28) =/= 0  unphysical global minimum
*      PROB(29) =/= 0  Higgs soft masses >> Msusy
*      PROB(32) =/= 0  b->s gamma more than 2 sigma away
*      PROB(33) =/= 0  Delta M_s more than 2 sigma away
*      PROB(34) =/= 0  Delta M_d more than 2 sigma away
*      PROB(35) =/= 0  B_s->mu+mu- more than 2 sigma away
*      PROB(36) =/= 0  B+-> tau+nu_tau more than 2 sigma away
*      PROB(37) =/= 0  (g-2)_muon more than 2 sigma away
*      PROB(38) =/= 0  excluded by Upsilon(1S) -> H/A gamma
*      PROB(39) =/= 0  excluded by eta_b(1S) mass measurement
*      PROB(41) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus (ALEPH analysis)
*      PROB(42) =/= 0  excluded by top -> b H+, H+ -> c s (CDF, D0)
*      PROB(43) =/= 0  excluded by top -> b H+, H+ -> tau nu_tau (D0)
*      PROB(44) =/= 0  excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus (CDF)
*      PROB(45) =/= 0  excluded by t -> bH+ (ATLAS)
*      PROB(46) =/= 0  No Higgs in the MHmin-MHmax GeV range
*      PROB(51) =/= 0: excluded by H/A->tautau (ATLAS+CMS)
*      PROB(52) =/= 0: excluded by H->AA->4leptons/2lept.+2b (ATLAS+CMS)
*      PROB(53) =/= 0: excluded by ggF->H/A->gamgam (ATLAS)
*      PROB(54) =/= 0: excluded from EDMs
*      PROB(55) =/= 0: b -> d gamma more than 2 sigma away
*      PROB(56) =/= 0: B_d -> mu+ mu- more than 2 sigma away
*      PROB(58) =/= 0: b -> c tau nu more than 2 sigma away (as SM)
*      PROB(63) =/= 0: excluded by H->AA->4gammas (ATLAS+CMS)
*      PROB(67) =/= 0 k_WZ(H_SM) 2 sigma away from LHC measured value
*      PROB(68) =/= 0 k_top(H_SM) 2 sigma away from LHC measured value
*      PROB(69) =/= 0 k_bot(H_SM) 2 sigma away from LHC measured value
*      PROB(70) =/= 0 k_glu(H_SM) 2 sigma away from LHC measured value
*      PROB(71) =/= 0 k_gam(H_SM) 2 sigma away from LHC measured value
*      PROB(72) =/= 0 k_tau(H_SM) 2 sigma away from LHC measured value
*      PROB(73) =/= 0 B_bsm(H_SM) 2 sigma away from LHC measured value
*      PROB(74) =/= 0: excluded by HSM->Z+A (ATLAS)
*      PROB(75) =/= 0: excluded by H/A->toptop (CMS)
*      PROB(76) =/= 0: excluded by A->Z+HSM (CMS)
*      PROB(77) =/= 0: excluded by A->Z+HSM (ATLAS)
*      PROB(78) =/= 0: excluded by H/A->Z+A/H (CMS)
*      PROB(79) =/= 0: excluded by H/A->Z+A/H (ATLAS)
*      PROB(80) =/= 0: excluded by H/A->HSM+H/A->2b2tau (CMS)
*      PROB(81) =/= 0: excluded by H/A->HSM+H/A->4b (CMS)
*      PROB(82) =/= 0: excluded by H->HSM+HSM (ATLAS)
*
************************************************************************

      IMPLICIT NONE

      INTEGER NFL,NPROB,NPAR
      PARAMETER (NFL=10,NPROB=82,NPAR=25)
      INTEGER NFAIL(NFL),IFAIL,I,TOT,ITOT,NTOT,IDUM
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG
      INTEGER MSFLAG,AKFLAG,ALFLAG,XISIFLAG
      INTEGER OMGFLAG,MAFLAG,MOFLAG,GMUFLAG,HFLAG
      INTEGER TOTMIN,TOTMAX,NMAX,IP
      INTEGER CFLAG(5)

      DOUBLE PRECISION PAR(NPAR),PROB(NPROB),GAU,Q2,PI
      DOUBLE PRECISION LCEN,LDEV,KCEN,KDEV,TBCEN,TBDEV,MUCEN,MUDEV,
     . ALCEN,ALDEV,AKCEN,AKDEV,XIFCEN,XIFDEV,XISCEN,XISDEV,MUPCEN,
     . MUPDEV,MSPCEN,MSPDEV,M3HCEN,M3HDEV,MACEN,MADEV,MPCEN,MPDEV,
     . M1CEN,M1DEV,M2CEN,M2DEV,M3CEN,M3DEV,XCEN,XDEV,X,LMIN,KMIN,
     . TBMIN,MUMIN,ALMIN,AKMIN,XIFMIN,XISMIN,MUPMIN,MSPMIN,M3HMIN,
     . MAMIN,MPMIN,M1MIN,M2MIN,M3MIN,
     . AQCEN,AQDEV,AQMIN,MQCEN,MQDEV,MQMIN
      DOUBLE PRECISION M1ICEN,M1IDEV,M2ICEN,M2IDEV,M3ICEN,M3IDEV,
     . AU3ICEN,AU3IDEV,AD3ICEN,AD3IDEV,AE3ICEN,AE3IDEV,LICEN,
     . LIDEV,KICEN,KIDEV,XIFICEN,XIFIDEV,XISICEN,XISIDEV,
     . MUPICEN,MUPIDEV,MSPICEN,MSPIDEV,M3HICEN,M3HIDEV,
     . M1IMIN,M2IMIN,M3IMIN,AU3IMIN,AD3IMIN,AE3IMIN,LIMIN,
     . KIMIN,XIFIMIN,XISIMIN,MUPIMIN,MSPIMIN,M3HIMIN
      DOUBLE PRECISION LN,LNN,KN,KNN,TBN,TBNN,MUN,MUNN
      DOUBLE PRECISION ALN,ALNN,AKN,AKNN,XIFN,XIFNN
      DOUBLE PRECISION XISN,XISNN,MUPN,MUPNN,MSPN,MSPNN
      DOUBLE PRECISION M3HN,M3HNN,MAN,MANN,MPN,MPNN
      DOUBLE PRECISION M1N,M1NN,M2N,M2NN,M3N,M3NN
      DOUBLE PRECISION AQN,AQNN,MQNN,MQN
      DOUBLE PRECISION M1IN,M1INN,M2IN,M2INN,M3IN,M3INN
      DOUBLE PRECISION AU3IN,AU3INN,AD3IN,AD3INN,AE3IN,AE3INN
      DOUBLE PRECISION LIN,LINN,KIN,KINN,XIFIN,XIFINN,XISIN,XISINN
      DOUBLE PRECISION MUPIN,MUPINN,MSPIN,MSPINN,M3HIN,M3HINN
      DOUBLE PRECISION ALIN,ALINN,AKIN,AKINN
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H,DELMB,DELML,DEL1
      DOUBLE PRECISION REALP(14),IMAGP(14),MYPHASES(16)
      DOUBLE PRECISION RXIS,RXIF,IAlQ2,IAkQ2,IXISQ2,MHuQ2,MHdQ2,MSQ2

      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . MSFLAG,AKFLAG,ALFLAG,XISIFLAG
      COMMON/STEPS/NTOT,IDUM,TOTMIN,TOTMAX,NMAX
      COMMON/BOUNDS/LN,LNN,KN,KNN,TBN,TBNN,MUN,MUNN,
     . ALN,ALNN,AKN,AKNN,XIFN,XIFNN,
     . XISN,XISNN,MUPN,MUPNN,MSPN,MSPNN,
     . M3HN,M3HNN,MAN,MANN,MPN,MPNN,
     . M1N,M1NN,M2N,M2NN,M3N,M3NN,
     . AQN,AQNN,MQN,MQNN
      COMMON/IBOUNDS/M1IN,M1INN,M2IN,M2INN,M3IN,M3INN,
     . AU3IN,AU3INN,AD3IN,AD3INN,AE3IN,AE3INN,LIN,LINN,
     . KIN,KINN,XIFIN,XIFINN,XISIN,XISINN,MUPIN,MUPINN,
     . MSPIN,MSPINN,M3HIN,M3HINN,ALIN,ALINN,AKIN,AKINN
      COMMON/MCMCPAR/LCEN,LDEV,KCEN,KDEV,TBCEN,TBDEV,MUCEN,MUDEV,
     . ALCEN,ALDEV,AKCEN,AKDEV,XIFCEN,XIFDEV,XISCEN,XISDEV,MUPCEN,
     . MUPDEV,MSPCEN,MSPDEV,M3HCEN,M3HDEV,MACEN,MADEV,MPCEN,MPDEV,
     . M1CEN,M1DEV,M2CEN,M2DEV,M3CEN,M3DEV,XCEN,XDEV,X,LMIN,KMIN,
     . TBMIN,MUMIN,ALMIN,AKMIN,XIFMIN,XISMIN,MUPMIN,MSPMIN,M3HMIN,
     . MAMIN,MPMIN,M1MIN,M2MIN,M3MIN,
     . AQCEN,AQDEV,AQMIN,MQCEN,MQDEV,MQMIN
      COMMON/MCMCIPAR/M1ICEN,M1IDEV,M2ICEN,M2IDEV,M3ICEN,M3IDEV,
     . AU3ICEN,AU3IDEV,AD3ICEN,AD3IDEV,AE3ICEN,AE3IDEV,LICEN,
     . LIDEV,KICEN,KIDEV,XIFICEN,XIFIDEV,XISICEN,XISIDEV,
     . MUPICEN,MUPIDEV,MSPICEN,MSPIDEV,M3HICEN,M3HIDEV,
     . M1IMIN,M2IMIN,M3IMIN,AU3IMIN,AD3IMIN,AE3IMIN,LIMIN,
     . KIMIN,XIFIMIN,XISIMIN,MUPIMIN,MSPIMIN,M3HIMIN
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/REAL_IMAG/REALP,IMAGP
      COMMON/MYPHASES/MYPHASES
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/RENSCALE/Q2
      COMMON/CFLAG/CFLAG
      COMMON/REXI/RXIS,RXIF,IAlQ2,IAkQ2,IXISQ2,MHuQ2,MHdQ2,MSQ2

      PI=4d0*DATAN(1d0)

*   Initialization

      CALL INITIALIZE()
      DO I=1,NFL
       NFAIL(I)=0
      ENDDO
      TOT=0
      IP=0

*   Reading of the input parameters

      CALL INPUT(PAR,NPAR)

*   Initialization of the range of parameters that has passed all tests

      TBN=1d99
      TBNN=-1d99
      M1N=1d99
      M1NN=-1d99
      M2N=1d99
      M2NN=-1d99
      M3N=1d99
      M3NN=-1d99
      LN=1d99
      LNN=-1d99
      KN=1d99
      KNN=-1d99
      MUN=1d99
      MUNN=-1d99
      MAN=1d99
      MANN=-1d99
      ALN=1d99
      ALNN=-1d99
      XIFN=1d99
      XIFNN=-1d99
      MPN=1d99
      MPNN=-1d99
      AKN=1d99
      AKNN=-1d99
      XISN=1d99
      XISNN=-1d99
      MUPN=1d99
      MUPNN=-1d99
      MSPN=1d99
      MSPNN=-1d99
      M3HN=1d99
      M3HNN=-1d99
      AQN=1d99
      AQNN=-1d99
      MQN=1d99
      MQNN=-1d99
      M1IN=1d99
      M1INN=-1d99
      M2IN=1d99
      M2INN=-1d99
      M3IN=1d99
      M3INN=-1d99
      AU3IN=1d99
      AU3INN=-1d99
      AD3IN=1d99
      AD3INN=-1d99
      AE3IN=1d99
      AE3INN=-1d99
      LIN=1d99
      LINN=-1d99
      KIN=1d99
      KINN=-1d99
      XIFIN=1d99
      XIFINN=-1d99
      XISIN=1d99
      XISINN=-1d99
      MUPIN=1d99
      MUPINN=-1d99
      MSPIN=1d99
      MSPINN=-1d99
      M3HIN=1d99
      M3HINN=-1d99
      ALIN=1d99
      ALINN=-1d99
      AKIN=1d99
      AKINN=-1d99

*   Beginning of the scan

      DO ITOT=1,NTOT

 14   IF(ITOT.EQ.1)THEN

       PAR(3)=TBCEN
       REALP(4)=M2CEN
       IF(M1FLAG.EQ.0)THEN
        REALP(3)=REALP(4)/2d0
       ELSE
        REALP(3)=M1CEN
       ENDIF
       IF(M3FLAG.EQ.0)THEN
        REALP(5)=REALP(4)*3d0
       ELSE
        REALP(5)=M3CEN
       ENDIF
       REALP(1)=LCEN
       REALP(2)=KCEN
       REALP(14)=MUCEN
       IF(MOD(MAFLAG,3).EQ.0)THEN
        PAR(5)=ALCEN
        REALP(9)=XIFCEN
       ELSEIF(MOD(MAFLAG,3).EQ.1)THEN
        PAR(23)=MACEN
        REALP(9)=XIFCEN
       ELSE
        PAR(5)=ALCEN
        PAR(23)=MACEN
        REALP(9)=0d0
       ENDIF
       IF(MAFLAG/3.EQ.0)THEN
        PAR(6)=AKCEN
        REALP(10)=XISCEN
       ELSEIF(MAFLAG/3.EQ.1)THEN
        PAR(24)=MPCEN
        REALP(10)=XISCEN
       ELSE
        PAR(6)=AKCEN
        PAR(24)=MPCEN
        REALP(10)=0d0
       ENDIF
       REALP(11)=MUPCEN
       REALP(12)=MSPCEN
       REALP(13)=M3HCEN
       REALP(6)=AQCEN
       PAR(7)=MQCEN**2

       IMAGP(1)=LICEN
       IMAGP(2)=KICEN
       IMAGP(3)=M1ICEN
       IMAGP(4)=M2ICEN
       IMAGP(5)=M3ICEN
       IMAGP(6)=AU3ICEN
       IMAGP(7)=AD3ICEN
       IMAGP(8)=AE3ICEN
       IMAGP(9)=XIFICEN
       IMAGP(10)=XISICEN
       IMAGP(11)=MUPICEN
       IMAGP(12)=MSPICEN
       IMAGP(13)=M3HICEN

      ELSE

       IF(TBDEV.EQ.0d0)THEN
        PAR(3)=TBCEN
       ELSE
        PAR(3)=TBCEN+MAX(DABS(TBCEN),TBMIN)*TBDEV*GAU(IDUM)
       ENDIF

       IF(M2DEV.EQ.0d0)THEN
        REALP(4)=M2CEN
       ELSE
        REALP(4)=M2CEN+MAX(DABS(M2CEN),M2MIN)*M2DEV*GAU(IDUM)
       ENDIF

       IF(M1FLAG.EQ.0)THEN
        REALP(3)=REALP(4)/2d0
       ELSEIF(M1DEV.EQ.0d0)THEN
        REALP(3)=M1CEN
       ELSE
        REALP(3)=M1CEN+MAX(DABS(M1CEN),M1MIN)*M1DEV*GAU(IDUM)
       ENDIF

       IF(M3FLAG.EQ.0)THEN
        REALP(5)=REALP(4)*3d0
       ELSEIF(M3DEV.EQ.0d0)THEN
        REALP(5)=M3CEN
       ELSE
        REALP(5)=M3CEN+MAX(DABS(M3CEN),M3MIN)*M3DEV*GAU(IDUM)
       ENDIF

       IF(LDEV.EQ.0d0)THEN
        REALP(1)=LCEN
       ELSE
        REALP(1)=LCEN+MAX(DABS(LCEN),LMIN)*LDEV*GAU(IDUM)
       ENDIF

       IF(KDEV.EQ.0d0)THEN
        REALP(2)=KCEN
       ELSE
        REALP(2)=KCEN+MAX(DABS(KCEN),KMIN)*KDEV*GAU(IDUM)
       ENDIF

       IF(MUDEV.EQ.0d0)THEN
        REALP(14)=MUCEN
       ELSE
        REALP(14)=MUCEN+MAX(DABS(MUCEN),MUMIN)*MUDEV*GAU(IDUM)
       ENDIF

       IF(MOD(MAFLAG,3).EQ.0)THEN
        IF(ALDEV.EQ.0d0)THEN
         PAR(5)=ALCEN
        ELSE
         PAR(5)=ALCEN+MAX(DABS(ALCEN),ALMIN)*ALDEV*GAU(IDUM)
        ENDIF
        IF(XIFDEV.EQ.0d0)THEN
         REALP(9)=XIFCEN
        ELSE
         REALP(9)=XIFCEN+MAX(DABS(XIFCEN),XIFMIN)*XIFDEV*GAU(IDUM)
        ENDIF
       ELSEIF(MOD(MAFLAG,3).EQ.1)THEN
        IF(MADEV.EQ.0d0)THEN
         PAR(23)=MACEN
        ELSE
         PAR(23)=MACEN+MAX(DABS(MACEN),MAMIN)*MADEV*GAU(IDUM)
        ENDIF
        IF(XIFDEV.EQ.0d0)THEN
         REALP(9)=XIFCEN
        ELSE
         REALP(9)=XIFCEN+MAX(DABS(XIFCEN),XIFMIN)*XIFDEV*GAU(IDUM)
        ENDIF
       ELSE
        IF(ALDEV.EQ.0d0)THEN
         PAR(5)=ALCEN
        ELSE
         PAR(5)=ALCEN+MAX(DABS(ALCEN),ALMIN)*ALDEV*GAU(IDUM)
        ENDIF
        IF(MADEV.EQ.0d0)THEN
         PAR(23)=MACEN
        ELSE
         PAR(23)=MACEN+MAX(DABS(MACEN),MAMIN)*MADEV*GAU(IDUM)
        ENDIF
        REALP(9)=0d0
       ENDIF

       IF(MAFLAG/3.EQ.0)THEN
        IF(AKDEV.EQ.0d0)THEN
         PAR(6)=AKCEN
        ELSE
         PAR(6)=AKCEN+MAX(DABS(AKCEN),AKMIN)*AKDEV*GAU(IDUM)
        ENDIF
        IF(XISDEV.EQ.0d0)THEN
         REALP(10)=XISCEN
        ELSE
         REALP(10)=XISCEN+MAX(DABS(XISCEN),XISMIN)*XISDEV*GAU(IDUM)
        ENDIF
       ELSEIF(MAFLAG/3.EQ.1)THEN
        IF(MPDEV.EQ.0d0)THEN
         PAR(24)=MPCEN
        ELSE
         PAR(24)=MPCEN+MAX(DABS(MPCEN),MPMIN)*MPDEV*GAU(IDUM)
        ENDIF
        IF(XISDEV.EQ.0d0)THEN
         REALP(10)=XISCEN
        ELSE
         REALP(10)=XISCEN+MAX(DABS(XISCEN),XISMIN)*XISDEV*GAU(IDUM)
        ENDIF
       ELSE
        IF(AKDEV.EQ.0d0)THEN
         PAR(6)=AKCEN
        ELSE
         PAR(6)=AKCEN+MAX(DABS(AKCEN),AKMIN)*AKDEV*GAU(IDUM)
        ENDIF
        IF(MPDEV.EQ.0d0)THEN
         PAR(24)=MPCEN
        ELSE
         PAR(24)=MPCEN+MAX(DABS(MPCEN),MPMIN)*MPDEV*GAU(IDUM)
        ENDIF
        REALP(10)=0d0
       ENDIF

       IF(MUPDEV.EQ.0d0)THEN
        REALP(11)=MUPCEN
       ELSE
        REALP(11)=MUPCEN+MAX(DABS(MUPCEN),MUPMIN)*MUPDEV*GAU(IDUM)
       ENDIF

       IF(MSPDEV.EQ.0d0)THEN
        REALP(12)=MSPCEN
       ELSE
        REALP(12)=MSPCEN+MAX(DABS(MSPCEN),MSPMIN)*MSPDEV*GAU(IDUM)
       ENDIF

       IF(M3HDEV.EQ.0d0)THEN
        REALP(13)=M3HCEN
       ELSE
        REALP(13)=M3HCEN+MAX(DABS(M3HCEN),M3HMIN)*M3HDEV*GAU(IDUM)
       ENDIF

       IF(AQDEV.EQ.0d0)THEN
        REALP(6)=AQCEN
       ELSE
        REALP(6)=AQCEN+MAX(DABS(AQCEN),AQMIN)*AQDEV*GAU(IDUM)
       ENDIF

       IF(MQDEV.EQ.0d0)THEN
        PAR(7)=MQCEN**2
       ELSE
        PAR(7)=(MQCEN+MAX(DABS(MQCEN),MQMIN)*MQDEV*GAU(IDUM))**2
       ENDIF

       IF(LICEN.EQ.0d0)THEN
        IMAGP(1)=LICEN
       ELSE
        IMAGP(1)=LICEN+MAX(DABS(LICEN),LIMIN)*LIDEV*GAU(IDUM)
       ENDIF

       IF(KICEN.EQ.0d0)THEN
        IMAGP(2)=KICEN
       ELSE
        IMAGP(2)=KICEN+MAX(DABS(KICEN),KIMIN)*KIDEV*GAU(IDUM)
       ENDIF

       IF(M1ICEN.EQ.0d0)THEN
        IMAGP(3)=M1ICEN
       ELSE
        IMAGP(3)=M1ICEN+MAX(DABS(M1ICEN),M1IMIN)*M1IDEV*GAU(IDUM)
       ENDIF

       IF(M2ICEN.EQ.0d0)THEN
        IMAGP(4)=M2ICEN
       ELSE
        IMAGP(4)=M2ICEN+MAX(DABS(M2ICEN),M2IMIN)*M2IDEV*GAU(IDUM)
       ENDIF

       IF(M3ICEN.EQ.0d0)THEN
        IMAGP(5)=M3ICEN
       ELSE
        IMAGP(5)=M3ICEN+MAX(DABS(M3ICEN),M3IMIN)*M3IDEV*GAU(IDUM)
       ENDIF

       IF(AU3ICEN.EQ.0d0)THEN
        IMAGP(6)=AU3ICEN
       ELSE
        IMAGP(6)=AU3ICEN+MAX(DABS(AU3ICEN),AU3IMIN)*AU3IDEV*GAU(IDUM)
       ENDIF

       IF(AD3ICEN.EQ.0d0)THEN
        IMAGP(7)=AD3ICEN
       ELSE
        IMAGP(7)=AD3ICEN+MAX(DABS(AD3ICEN),AD3IMIN)*AD3IDEV*GAU(IDUM)
       ENDIF

       IF(AE3ICEN.EQ.0d0)THEN
        IMAGP(8)=AE3ICEN
       ELSE
        IMAGP(8)=AE3ICEN+MAX(DABS(AE3ICEN),AE3IMIN)*AE3IDEV*GAU(IDUM)
       ENDIF

       IF(XIFICEN.EQ.0d0)THEN
        IMAGP(9)=XIFICEN
       ELSE
        IMAGP(9)=XIFICEN+MAX(DABS(XIFICEN),XIFIMIN)*XIFIDEV*GAU(IDUM)
       ENDIF

       IF(XISICEN.EQ.0d0)THEN
        IMAGP(10)=XISICEN
       ELSE
        IMAGP(10)=XISICEN+MAX(DABS(XISICEN),XISIMIN)*XISIDEV*GAU(IDUM)
       ENDIF

       IF(MUPICEN.EQ.0d0)THEN
        IMAGP(11)=MUPICEN
       ELSE
        IMAGP(11)=MUPICEN+MAX(DABS(MUPICEN),MUPIMIN)*MUPIDEV*GAU(IDUM)
       ENDIF

       IF(MSPICEN.EQ.0d0)THEN
        IMAGP(12)=MSPICEN
       ELSE
        IMAGP(12)=MSPICEN+MAX(DABS(MSPICEN),MSPIMIN)*MSPIDEV*GAU(IDUM)
       ENDIF

       IF(M3HICEN.EQ.0d0)THEN
        IMAGP(13)=M3HICEN
       ELSE
        IMAGP(13)=M3HICEN+MAX(DABS(M3HICEN),M3HIMIN)*M3HIDEV*GAU(IDUM)
       ENDIF

      ENDIF

      REALP(7)=REALP(6)
      PAR(8)=PAR(7)
      PAR(9)=PAR(7)
      PAR(15)=PAR(7)
      PAR(16)=PAR(7)
      PAR(17)=PAR(7)

*   Moduli and phases

      PAR(1)=DSQRT(REALP(1)**2+IMAGP(1)**2)
      MYPHASES(1)=DATAN(IMAGP(1)/REALP(1))

      IF(REALP(2).NE.0d0)THEN
       PAR(2)=DSQRT(REALP(2)**2+IMAGP(2)**2)*REALP(2)/DABS(REALP(2))
       MYPHASES(2)=DATAN(IMAGP(2)/REALP(2))
      ELSE
       PAR(2)=IMAGP(2)
       IF(IMAGP(2).NE.0d0)MYPHASES(2)=SIGN(PI/2d0,IMAGP(2))
      ENDIF

      PAR(4)=REALP(14)*PAR(1)/REALP(1)
      IMAGP(14)=REALP(14)*IMAGP(1)/REALP(1)

      IF(REALP(3).NE.0d0)THEN
       PAR(20)=DSQRT(REALP(3)**2+IMAGP(3)**2)*REALP(3)/DABS(REALP(3))
       MYPHASES(3)=DATAN(IMAGP(3)/REALP(3))
      ELSE
       PAR(20)=IMAGP(3)
       IF(IMAGP(3).NE.0d0)MYPHASES(3)=SIGN(PI/2d0,IMAGP(3))
      ENDIF

      IF(REALP(4).NE.0d0)THEN
       PAR(21)=DSQRT(REALP(4)**2+IMAGP(4)**2)*REALP(4)/DABS(REALP(4))
       MYPHASES(4)=DATAN(IMAGP(4)/REALP(4))
      ELSE
       PAR(21)=IMAGP(4)
       IF(IMAGP(4).NE.0d0)MYPHASES(4)=SIGN(PI/2d0,IMAGP(4))
      ENDIF

      IF(REALP(5).NE.0d0)THEN
       PAR(22)=DSQRT(REALP(5)**2+IMAGP(5)**2)*REALP(5)/DABS(REALP(5))
       MYPHASES(5)=DATAN(IMAGP(5)/REALP(5))
      ELSE
       PAR(22)=IMAGP(5)
       IF(IMAGP(5).NE.0d0)MYPHASES(5)=SIGN(PI/2d0,IMAGP(5))
      ENDIF

      IF(REALP(6).NE.0d0)THEN
       PAR(12)=DSQRT(REALP(6)**2+IMAGP(6)**2)*REALP(6)/DABS(REALP(6))
       MYPHASES(6)=DATAN(IMAGP(6)/REALP(6))
      ELSE
       PAR(12)=IMAGP(6)
       IF(IMAGP(6).NE.0d0)MYPHASES(6)=SIGN(PI/2d0,IMAGP(6))
      ENDIF

      IF(REALP(7).NE.0d0)THEN
       PAR(13)=DSQRT(REALP(7)**2+IMAGP(7)**2)*REALP(7)/DABS(REALP(7))
       MYPHASES(7)=DATAN(IMAGP(7)/REALP(7))
      ELSE
       PAR(13)=IMAGP(7)
       IF(IMAGP(7).NE.0d0)MYPHASES(7)=SIGN(PI/2d0,IMAGP(7))
      ENDIF

      IF(REALP(8).NE.0d0)THEN
       PAR(14)=DSQRT(REALP(8)**2+IMAGP(8)**2)*REALP(8)/DABS(REALP(8))
       MYPHASES(8)=DATAN(IMAGP(8)/REALP(8))
      ELSE
       PAR(14)=IMAGP(8)
       IF(IMAGP(8).NE.0d0)MYPHASES(8)=SIGN(PI/2d0,IMAGP(8))
      ENDIF

      IF(REALP(11).NE.0d0)THEN
       MUP=DSQRT(REALP(11)**2+IMAGP(11)**2)*REALP(11)/DABS(REALP(11))
       MYPHASES(15)=DATAN(IMAGP(11)/REALP(11))
      ELSE
       MUP=IMAGP(11)
       IF(IMAGP(11).NE.0d0)MYPHASES(15)=SIGN(PI/2d0,IMAGP(11))
      ENDIF

      IF(REALP(12).NE.0d0)THEN
       MSP=DSQRT(REALP(12)**2+IMAGP(12)**2)*REALP(12)/DABS(REALP(12))
       MYPHASES(13)=DATAN(IMAGP(12)/REALP(12))
      ELSE
       MSP=IMAGP(12)
       IF(IMAGP(12).NE.0d0)MYPHASES(13)=SIGN(PI/2d0,IMAGP(12))
      ENDIF

      IF(REALP(13).NE.0d0)THEN
       M3H=DSQRT(REALP(13)**2+IMAGP(13)**2)*REALP(13)/DABS(REALP(13))
       MYPHASES(16)=DATAN(IMAGP(13)/REALP(13))
      ELSE
       M3H=IMAGP(13)
       IF(IMAGP(13).NE.0d0)MYPHASES(16)=SIGN(PI/2d0,IMAGP(13))
      ENDIF

      IF(REALP(9).NE.0d0)THEN
       XIF=DSQRT(REALP(9)**2+IMAGP(9)**2)*REALP(9)/DABS(REALP(9))
       MYPHASES(14)=DATAN(IMAGP(9)/REALP(9))
      ELSE
       XIF=IMAGP(9)
       IF(IMAGP(9).NE.0d0)MYPHASES(14)=SIGN(PI/2d0,IMAGP(9))
      ENDIF

      IF(REALP(10).NE.0d0)THEN
       XIS=DSQRT(REALP(10)**2+IMAGP(10)**2)*REALP(10)/DABS(REALP(10))
       MYPHASES(12)=DATAN(IMAGP(10)/REALP(10))
      ELSE
       XIS=IMAGP(10)
       IF(IMAGP(10).NE.0d0)MYPHASES(12)=SIGN(PI/2d0,IMAGP(10))
      ENDIF

!      WRITE(0,*)""
!      WRITE(0,*)"------------------------------------------------------"
!      WRITE(0,*)""
!      WRITE(0,*)"Point ",ITOT
!      WRITE(0,*)""
!      WRITE(0,*)"MAFLAG=",MAFLAG
!      WRITE(0,*)""
!      WRITE(0,*)"TANB =",PAR(3)
!      WRITE(0,*)"M1 =",PAR(20)
!      WRITE(0,*)"M2 =",PAR(21)
!      WRITE(0,*)"M3 =",PAR(22)
!      WRITE(0,*)"LAMBDA =",PAR(1)
!      WRITE(0,*)"KAPPA =",PAR(2)
!      WRITE(0,*)"MUEFF =",PAR(4)
!      WRITE(0,*)"ALAMBDA =",PAR(5)
!      WRITE(0,*)"AKAPPA =",PAR(6)
!      WRITE(0,*)"XIF =",XIF
!      WRITE(0,*)"XIS =",XIS
!      WRITE(0,*)"MUP =",MUP
!      WRITE(0,*)"MSP =",MSP
!      WRITE(0,*)"M3H =",M3H
!      WRITE(0,*)"MA =",PAR(23)
!      WRITE(0,*)"MP =",PAR(24)
!      WRITE(0,*)""

*   Initialization of PROB and IFAIL

      DO I=1,NPROB
       PROB(I)=0d0
      ENDDO
      IFAIL=0
      DELMB=.1d0
      DELML=0d0
      DEL1=0d0

*   Check for mu = 0 | k=0 & (Ak|XiS)=/=0

      IF(PAR(4).EQ.0d0 .OR. (PAR(2).EQ.0d0 .AND. (MAFLAG/3.EQ.1
     ..OR.XISIFLAG.EQ.1)))THEN
       IFAIL=9
       GOTO 11
      ENDIF

*   Computation of parameters at QSTSB

      CALL runpar_CPV(PAR)

      CALL mcha_CPV(IFAIL)
      IF(IFAIL.ne.0)THEN
       GOTO 11
      ENDIF

      CALL mneu_CPV(IFAIL)
      IF(IFAIL.ne.0)THEN
       GOTO 11
      ENDIF

      CALL mhiggstree_CPV()

      CALL msferm_CPV(PAR,IFAIL)
      IF(IFAIL.ne.0)THEN
       GOTO 11
      ENDIF

      CALL GETSUSYCOUP(PAR)

      CALL mgluino_CPV(PAR)

      CALL mhiggsloop_sferm_CPV()
      CALL mhiggsloop_inos_CPV()
      CALL mhiggsloop_gaugehiggs_CPV()
      CALL mhiggsloop_pole_CPV(IFAIL)
      IF(IFAIL.ne.0)THEN
       GOTO 11
      ENDIF

      CALL susycoup_CPV()
      CALL higgscoup_CPV(PAR)
      CALL hidecay_CPV()
      CALL tdecay_CPV()

      IF(CFLAG(2).EQ.1)THEN
       CALL constsusypart_CPV(PROB)
       CALL LEP_Higgs_CPV(PROB)
       CALL tevatron_chiggs_CPV(PROB)
      ENDIF

      IF(CFLAG(3).EQ.1)THEN
       CALL LHC_HIGGS_CPV(PROB)
       PROB(69)=0d0
       PROB(72)=0d0
      ENDIF

      IF(CFLAG(4).EQ.1)THEN
       CALL bottomonium_CPV(PROB)
       CALL bsg_CPV(PAR,PROB)
      ENDIF

      IF(GMUFLAG.EQ.1)CALL magnmu_CPV(PROB)

      CALL checkmin_CPV(PAR,PROB)

      IF(CFLAG(5).EQ.1)THEN
       CALL EDM_CPV(PAR,PROB)
      ENDIF

*   Check for problems

      DO I=1,NPROB
       IF(PROB(I).NE.0d0)THEN
!        WRITE(0,*)"PROB",I,PROB(I)
        IFAIL=10
       ENDIF
      ENDDO
!      WRITE(0,*)""

*   Recording of the results

 11   CALL MCMCSTEPCPV(PAR,PROB,NPROB,IFAIL)
      CALL OUTPUT(PAR,PROB,IFAIL)
      IF(IFAIL.EQ.0)THEN
       TOT=TOT+1
       LN=MIN(REALP(1),LN)
       LNN=MAX(REALP(1),LNN)
       KN=MIN(REALP(2),KN)
       KNN=MAX(REALP(2),KNN)
       TBN=MIN(PAR(3),TBN)
       TBNN=MAX(PAR(3),TBNN)
       MUN=MIN(REALP(14),MUN)
       MUNN=MAX(REALP(14),MUNN)
       ALN=MIN(PAR(5),ALN)
       ALNN=MAX(PAR(5),ALNN)
       AKN=MIN(PAR(6),AKN)
       AKNN=MAX(PAR(6),AKNN)
       XIFN=MIN(REALP(9),XIFN)
       XIFNN=MAX(REALP(9),XIFNN)
       XISN=MIN(REALP(10),XISN)
       XISNN=MAX(REALP(10),XISNN)
       MUPN=MIN(REALP(11),MUPN)
       MUPNN=MAX(REALP(11),MUPNN)
       MSPN=MIN(REALP(12),MSPN)
       MSPNN=MAX(REALP(12),MSPNN)
       M3HN=MIN(REALP(13),M3HN)
       M3HNN=MAX(REALP(13),M3HNN)
       M1N=MIN(REALP(3),M1N)
       M1NN=MAX(REALP(3),M1NN)
       M2N=MIN(REALP(4),M2N)
       M2NN=MAX(REALP(4),M2NN)
       M3N=MIN(REALP(5),M3N)
       M3NN=MAX(REALP(5),M3NN)
       MAN=MIN(PAR(23),MAN)
       MANN=MAX(PAR(23),MANN)
       MPN=MIN(PAR(24),MPN)
       MPNN=MAX(PAR(24),MPNN)
       AQN=MIN(REALP(6),AQN)
       AQNN=MAX(REALP(6),AQNN)
       MQN=MIN(PAR(7),MQN)
       MQNN=MAX(PAR(7),MQNN)
       M1IN=MIN(IMAGP(3),M1IN)
       M1INN=MAX(IMAGP(3),M1INN)
       M2IN=MIN(IMAGP(4),M2IN)
       M2INN=MAX(IMAGP(4),M2INN)
       M3IN=MIN(IMAGP(5),M3IN)
       M3INN=MAX(IMAGP(5),M3INN)
       AU3IN=MIN(IMAGP(6),AU3IN)
       AU3INN=MAX(IMAGP(6),AU3INN)
       AD3IN=MIN(IMAGP(7),AD3IN)
       AD3INN=MAX(IMAGP(7),AD3INN)
       AE3IN=MIN(IMAGP(8),AE3IN)
       AE3INN=MAX(IMAGP(8),AE3INN)
       LIN=MIN(IMAGP(1),LIN)
       LINN=MAX(IMAGP(1),LINN)
       KIN=MIN(IMAGP(2),KIN)
       KINN=MAX(IMAGP(2),KINN)
       XIFIN=MIN(IMAGP(9),XIFIN)
       XIFINN=MAX(IMAGP(9),XIFINN)
       XISIN=MIN(IMAGP(10),XISIN)
       XISINN=MAX(IMAGP(10),XISINN)
       MUPIN=MIN(IMAGP(11),MUPIN)
       MUPINN=MAX(IMAGP(11),MUPINN)
       MSPIN=MIN(IMAGP(12),MSPIN)
       MSPINN=MAX(IMAGP(12),MSPINN)
       M3HIN=MIN(IMAGP(13),M3HIN)
       M3HINN=MAX(IMAGP(13),M3HINN)
       ALIN=MIN(IALQ2,ALIN)
       ALINN=MAX(IALQ2,ALINN)
       AKIN=MIN(IAKQ2,AKIN)
       AKINN=MAX(IAKQ2,AKINN)
      ELSE
       NFAIL(IFAIL)=NFAIL(IFAIL)+1
      ENDIF

      IF(TOT.EQ.TOTMAX)THEN
       NTOT=ITOT
       GOTO 12
      ENDIF
      IF(IP.EQ.1)GOTO 13

      ENDDO

      IP=1
 13   IF(TOT.LT.TOTMIN .AND. NTOT.LT.NMAX)THEN
       NTOT=NTOT+1
       ITOT=ITOT+1
       GOTO 14
      ENDIF

*   Summary of the scanning:
*   Number of points that passed/failed the tests
*   and range for scanned parameters

 12   CALL ERROR(TOT,NTOT,NFAIL)

      END


      SUBROUTINE INPUT(PAR,NPAR)

*******************************************************************
*   This subroutine reads SM and NMSSM parameters from input file   .
*******************************************************************

      IMPLICIT NONE

      CHARACTER CHINL*120,CHBLCK*60,CHDUM*120

      INTEGER I,NLINE,INL,ICH,IX,IVAL,Q2FIX,VFLAG,Z3FLAG
      INTEGER N0,NLOOP,NBER,NPAR,ERR,GMUFLAG,HFLAG
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,MSFLAG
      INTEGER AKFLAG,ALFLAG,XISIFLAG,OMGFLAG,MAFLAG,MOFLAG
      INTEGER NTOT,ISEED,TOTMIN,TOTMAX,NMAX,MCFLAG
      INTEGER CFLAG(5)

      DOUBLE PRECISION PAR(*),VAL
      DOUBLE PRECISION ACC,XITLA,XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,Q2,Q2MIN
      DOUBLE PRECISION LCEN,LDEV,KCEN,KDEV,TBCEN,TBDEV,MUCEN,MUDEV,
     . ALCEN,ALDEV,AKCEN,AKDEV,XIFCEN,XIFDEV,XISCEN,XISDEV,MUPCEN,
     . MUPDEV,MSPCEN,MSPDEV,M3HCEN,M3HDEV,MACEN,MADEV,MPCEN,MPDEV,
     . M1CEN,M1DEV,M2CEN,M2DEV,M3CEN,M3DEV,XCEN,XDEV,X,LMIN,KMIN,
     . TBMIN,MUMIN,ALMIN,AKMIN,XIFMIN,XISMIN,MUPMIN,MSPMIN,M3HMIN,
     . MAMIN,MPMIN,M1MIN,M2MIN,M3MIN,
     . AQCEN,AQDEV,AQMIN,MQCEN,MQDEV,MQMIN
      DOUBLE PRECISION M1ICEN,M1IDEV,M2ICEN,M2IDEV,M3ICEN,M3IDEV,
     . AU3ICEN,AU3IDEV,AD3ICEN,AD3IDEV,AE3ICEN,AE3IDEV,LICEN,
     . LIDEV,KICEN,KIDEV,XIFICEN,XIFIDEV,XISICEN,XISIDEV,
     . MUPICEN,MUPIDEV,MSPICEN,MSPIDEV,M3HICEN,M3HIDEV,
     . M1IMIN,M2IMIN,M3IMIN,AU3IMIN,AD3IMIN,AE3IMIN,LIMIN,
     . KIMIN,XIFIMIN,XISIMIN,MUPIMIN,MSPIMIN,M3HIMIN
      DOUBLE PRECISION REALP(14),IMAGP(14)

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/MCMCPAR/LCEN,LDEV,KCEN,KDEV,TBCEN,TBDEV,MUCEN,MUDEV,
     . ALCEN,ALDEV,AKCEN,AKDEV,XIFCEN,XIFDEV,XISCEN,XISDEV,MUPCEN,
     . MUPDEV,MSPCEN,MSPDEV,M3HCEN,M3HDEV,MACEN,MADEV,MPCEN,MPDEV,
     . M1CEN,M1DEV,M2CEN,M2DEV,M3CEN,M3DEV,XCEN,XDEV,X,LMIN,KMIN,
     . TBMIN,MUMIN,ALMIN,AKMIN,XIFMIN,XISMIN,MUPMIN,MSPMIN,M3HMIN,
     . MAMIN,MPMIN,M1MIN,M2MIN,M3MIN,
     . AQCEN,AQDEV,AQMIN,MQCEN,MQDEV,MQMIN
      COMMON/MCMCIPAR/M1ICEN,M1IDEV,M2ICEN,M2IDEV,M3ICEN,M3IDEV,
     . AU3ICEN,AU3IDEV,AD3ICEN,AD3IDEV,AE3ICEN,AE3IDEV,LICEN,
     . LIDEV,KICEN,KIDEV,XIFICEN,XIFIDEV,XISICEN,XISIDEV,
     . MUPICEN,MUPIDEV,MSPICEN,MSPIDEV,M3HICEN,M3HIDEV,
     . M1IMIN,M2IMIN,M3IMIN,AU3IMIN,AD3IMIN,AE3IMIN,LIMIN,
     . KIMIN,XIFIMIN,XISIMIN,MUPIMIN,MSPIMIN,M3HIMIN
      COMMON/STEPS/NTOT,ISEED,TOTMIN,TOTMAX,NMAX
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . MSFLAG,AKFLAG,ALFLAG,XISIFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/VFLAG/VFLAG
      COMMON/MCFLAG/MCFLAG
      COMMON/CFLAG/CFLAG
      COMMON/REAL_IMAG/REALP,IMAGP

*   INITIALIZATION OF THE SUSY PARAMETERS
      DO I=1,NPAR
       PAR(I)=1d99
      ENDDO

*   DEFAULT VALUE FOR REALP,IMAGP
      DO I=1,14
       REALP(I)=1d99
       IMAGP(I)=0d0
      ENDDO

*   INITIALIZATION OF THE SCANNING PARAMETERS
      TBCEN=1d99
      TBDEV=1d99
      TBMIN=0d0
      M1CEN=1d99
      M1DEV=1d99
      M1MIN=0d0
      M2CEN=1d99
      M2DEV=1d99
      M2MIN=0d0
      M3CEN=1d99
      M3DEV=1d99
      M3MIN=0d0
      LCEN=1d99
      LDEV=1d99
      LMIN=0d0
      KCEN=0d0
      KDEV=1d99
      KMIN=0d0
      ALCEN=1d99
      ALDEV=1d99
      ALMIN=0d0
      AKCEN=1d99
      AKDEV=1d99
      AKMIN=0d0
      MUCEN=1d99
      MUDEV=1d99
      MUMIN=0d0
      XIFCEN=1d99
      XIFDEV=1d99
      XIFMIN=0d0
      XISCEN=1d99
      XISDEV=1d99
      XISMIN=0d0
      MUPCEN=0d0
      MUPDEV=1d99
      MUPMIN=0d0
      MSPCEN=0d0
      MSPDEV=1d99
      MSPMIN=0d0
      M3HCEN=0d0
      M3HDEV=1d99
      M3HMIN=0d0
      MACEN=1d99
      MADEV=1d99
      MAMIN=0d0
      MPMIN=0d0
      MPCEN=1d99
      MPDEV=1d99
      AQCEN=1d99
      AQDEV=1d99
      AQMIN=0d0
      MQCEN=1d99
      MQDEV=1d99
      MQMIN=0d0
      XCEN=1d99
      XDEV=1d0
      LICEN=0d0
      LIDEV=1d99
      LIMIN=0d0
      KICEN=0d0
      KIDEV=1d99
      KIMIN=0d0
      M1ICEN=0d0
      M1IDEV=1d99
      M1IMIN=0d0
      M2ICEN=0d0
      M2IDEV=1d99
      M2IMIN=0d0
      M3ICEN=0d0
      M3IDEV=1d99
      M3IMIN=0d0
      AU3ICEN=0d0
      AU3IDEV=1d99
      AU3IMIN=0d0
      AD3ICEN=0d0
      AD3IDEV=1d99
      AD3IMIN=0d0
      AE3ICEN=0d0
      AE3IDEV=1d99
      AE3IMIN=0d0
      XIFICEN=0d0
      XIFIDEV=1d99
      XIFIMIN=0d0
      XISICEN=1d99
      XISIDEV=1d99
      XISIMIN=0d0
      MUPICEN=0d0
      MUPIDEV=1d99
      MUPIMIN=0d0
      MSPICEN=0d0
      MSPIDEV=1d99
      MSPIMIN=0d0
      M3HICEN=0d0
      M3HIDEV=1d99
      M3HIMIN=0d0
      NTOT=0
      TOTMIN=0
      TOTMAX=1000000
      NMAX=1000000


*   DEFAULT VALUES FOR FLAGS
      GMUFLAG=1
      HFLAG=0
      OMGFLAG=0
      MOFLAG=7
      VFLAG=0
      M1FLAG=0
      M3FLAG=0
      MCFLAG=0
      XISIFLAG=1
      DO I=1,5
       CFLAG(I)=1
      ENDDO

*   DEFAULT VALUE FOR THE RANDOM SEED
      ISEED=-1

*   DEFAULT VALUE FOR THE RENSCALE Q2
      Q2=0d0

*   INITIALIZE READ LOOP
      ERR=0
      NLINE=0
      CHBLCK=' '

*   START TO READ NEW LINE INTO CHINL
 21   CHINL=' '

*   LINE NUMBER
      NLINE=NLINE+1

      READ(5,'(A120)',END=29,ERR=999) CHINL

*   CHECK FOR EMPTY OR COMMENT LINES
      IF(CHINL.EQ.' '.OR.CHINL(1:1).EQ.'#'
     .  .OR.CHINL(1:1).EQ.'*') GOTO 21

*   FORCE UPPER CASE LETTERS IN CHINL (AS REQUIRED BELOW)
      INL=0
 22   INL=INL+1
      IF(CHINL(INL:INL).NE.'#')THEN
       DO ICH=97,122
        IF(CHINL(INL:INL).EQ.CHAR(ICH)) CHINL(INL:INL)=CHAR(ICH-32)
       ENDDO
       IF(INL.LT.120) GOTO 22
      ENDIF

*   CHECK FOR BLOCK STATEMENT
      IF(CHINL(1:1).EQ.'B')THEN
       READ(CHINL,'(A6,A)',ERR=999) CHDUM,CHBLCK
       GOTO 21
      ENDIF

*   CHECK FOR NMSSM MODEL IN MODSEL
*   IF THE RELIC DENSITY SHOULD BE COMPUTED
*   THE BLOCK MODSEL MUST CONTAIN THE LINE "  9     1    "
      IF(CHBLCK(1:6).EQ.'MODSEL')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.1) Z3FLAG=IVAL
       IF(IX.EQ.11) GMUFLAG=IVAL
       IF(IX.EQ.14) VFLAG=IVAL
       IF(IX.EQ.17) CFLAG(1)=IVAL
       IF(IX.EQ.18) CFLAG(2)=IVAL
       IF(IX.EQ.19) CFLAG(3)=IVAL
       IF(IX.EQ.20) CFLAG(4)=IVAL
       IF(IX.EQ.21) CFLAG(5)=IVAL

*   READ SMINPUTS
      ELSEIF(CHBLCK(1:8).EQ.'SMINPUTS')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.2) GF=VAL
       IF(IX.EQ.3) ALSMZ=VAL
       IF(IX.EQ.4) MZ=VAL
       IF(IX.EQ.5) MB=VAL
       IF(IX.EQ.6) MT=VAL
       IF(IX.EQ.7) MTAU=VAL

*   READ Q2 AND TANBETA
      ELSEIF(CHBLCK(1:6).EQ.'MINPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0.AND.Q2.EQ.0d0) Q2=VAL**2
       IF(IX.EQ.3) TBCEN=VAL
       IF(IX.EQ.36) TBDEV=VAL
       IF(IX.EQ.37) TBMIN=VAL

*   READ EXTPAR
      ELSEIF(CHBLCK(1:6).EQ.'EXTPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0) Q2=VAL**2
       IF(IX.EQ.-1) XDEV=VAL
       IF(IX.EQ.1) M1CEN=VAL
       IF(IX.EQ.106) M1DEV=VAL
       IF(IX.EQ.107) M1MIN=VAL
       IF(IX.EQ.2) M2CEN=VAL
       IF(IX.EQ.206) M2DEV=VAL
       IF(IX.EQ.207) M2MIN=VAL
       IF(IX.EQ.3) M3CEN=VAL
       IF(IX.EQ.306) M3DEV=VAL
       IF(IX.EQ.307) M3MIN=VAL
       IF(IX.EQ.11) AQCEN=VAL
       IF(IX.EQ.116) AQDEV=VAL
       IF(IX.EQ.117) AQMIN=VAL
       IF(IX.EQ.13) REALP(8)=VAL
       IF(IX.EQ.16) PAR(25)=VAL
       IF(IX.EQ.32) PAR(18)=VAL**2
       IF(IX.EQ.33) PAR(10)=VAL**2
       IF(IX.EQ.35) PAR(19)=VAL**2
       IF(IX.EQ.36) PAR(11)=VAL**2
       IF(IX.EQ.43) MQCEN=VAL
       IF(IX.EQ.436) MQDEV=VAL
       IF(IX.EQ.437) MQMIN=VAL
       IF(IX.EQ.61) LCEN=VAL
       IF(IX.EQ.616) LDEV=VAL
       IF(IX.EQ.617) LMIN=VAL
       IF(IX.EQ.62) KCEN=VAL
       IF(IX.EQ.626) KDEV=VAL
       IF(IX.EQ.627) KMIN=VAL
       IF(IX.EQ.63) ALCEN=VAL
       IF(IX.EQ.636) ALDEV=VAL
       IF(IX.EQ.637) ALMIN=VAL
       IF(IX.EQ.64) AKCEN=VAL
       IF(IX.EQ.646) AKDEV=VAL
       IF(IX.EQ.647) AKMIN=VAL
       IF(IX.EQ.65) MUCEN=VAL
       IF(IX.EQ.656) MUDEV=VAL
       IF(IX.EQ.657) MUMIN=VAL
       IF(IX.EQ.66) XIFCEN=VAL
       IF(IX.EQ.666) XIFDEV=VAL
       IF(IX.EQ.667) XIFMIN=VAL
       IF(IX.EQ.67) XISCEN=VAL
       IF(IX.EQ.676) XISDEV=VAL
       IF(IX.EQ.677) XISMIN=VAL
       IF(IX.EQ.68) MUPCEN=VAL
       IF(IX.EQ.686) MUPDEV=VAL
       IF(IX.EQ.687) MUPMIN=VAL
       IF(IX.EQ.69) MSPCEN=VAL
       IF(IX.EQ.696) MSPDEV=VAL
       IF(IX.EQ.697) MSPMIN=VAL
       IF(IX.EQ.72) M3HCEN=VAL
       IF(IX.EQ.726) M3HDEV=VAL
       IF(IX.EQ.727) M3HMIN=VAL
       IF(IX.EQ.124) MACEN=VAL
       IF(IX.EQ.1246) MADEV=VAL
       IF(IX.EQ.1247) MAMIN=VAL
       IF(IX.EQ.125) MPCEN=VAL
       IF(IX.EQ.1256) MPDEV=VAL
       IF(IX.EQ.1257) MPMIN=VAL

*   READ IMEXTPAR
      ELSEIF(CHBLCK(1:8).EQ.'IMEXTPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.1) M1ICEN=VAL
       IF(IX.EQ.106) M1IDEV=VAL
       IF(IX.EQ.107) M1IMIN=VAL
       IF(IX.EQ.2) M2ICEN=VAL
       IF(IX.EQ.206) M2IDEV=VAL
       IF(IX.EQ.207) M2IMIN=VAL
       IF(IX.EQ.3) M3ICEN=VAL
       IF(IX.EQ.306) M3IDEV=VAL
       IF(IX.EQ.307) M3IMIN=VAL
       IF(IX.EQ.11) AU3ICEN=VAL
       IF(IX.EQ.116) AU3IDEV=VAL
       IF(IX.EQ.117) AU3IMIN=VAL
       IF(IX.EQ.12) AD3ICEN=VAL
       IF(IX.EQ.126) AD3IDEV=VAL
       IF(IX.EQ.127) AD3IMIN=VAL
       IF(IX.EQ.13) AE3ICEN=VAL
       IF(IX.EQ.136) AE3IDEV=VAL
       IF(IX.EQ.137) AE3IMIN=VAL
       IF(IX.EQ.61) LICEN=VAL
       IF(IX.EQ.616) LIDEV=VAL
       IF(IX.EQ.617) LIMIN=VAL
       IF(IX.EQ.62) KICEN=VAL
       IF(IX.EQ.626) KIDEV=VAL
       IF(IX.EQ.627) KIMIN=VAL
       IF(IX.EQ.66) XIFICEN=VAL
       IF(IX.EQ.666) XIFIDEV=VAL
       IF(IX.EQ.667) XIFIMIN=VAL
       IF(IX.EQ.67) XISICEN=VAL
       IF(IX.EQ.676) XISIDEV=VAL
       IF(IX.EQ.677) XISIMIN=VAL
       IF(IX.EQ.68) MUPICEN=VAL
       IF(IX.EQ.686) MUPIDEV=VAL
       IF(IX.EQ.687) MUPIMIN=VAL
       IF(IX.EQ.69) MSPICEN=VAL
       IF(IX.EQ.696) MSPIDEV=VAL
       IF(IX.EQ.697) MSPIMIN=VAL
       IF(IX.EQ.72) M3HICEN=VAL
       IF(IX.EQ.726) M3HIDEV=VAL
       IF(IX.EQ.727) M3HIMIN=VAL
       IF(IX.EQ.63 .OR. IX.EQ.636 .OR. IX.EQ.637)THEN
        WRITE(0,1)"IM(ALAMBDA) IS NOT AN INPUT"
        ERR=1
       ENDIF
       IF(IX.EQ.64 .OR. IX.EQ.646 .OR. IX.EQ.647)THEN
        WRITE(0,1)"IM(AKAPPA) IS NOT AN INPUT"
        ERR=1
       ENDIF
       IF(IX.EQ.65 .OR. IX.EQ.656 .OR. IX.EQ.657)THEN
        WRITE(0,1)"IM(MUEFF) IS NOT AN INPUT"
        ERR=1
       ENDIF

*   READ STEPS
       ELSEIF(CHBLCK(1:6).EQ.'STEPS')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.0) NTOT=IVAL
       IF(IX.EQ.1) ISEED=IVAL
       IF(IX.EQ.3) MCFLAG=IVAL
       IF(IX.EQ.4) TOTMIN=IVAL
       IF(IX.EQ.5) TOTMAX=IVAL
       IF(IX.EQ.6) NMAX=IVAL

      ENDIF

      GOTO 21

*   END OF READING FROM INPUT FILE

*   Check for errors

 29   DO I=1,5
       IF(CFLAG(I).LT.0 .OR. CFLAG(I).GT.1)THEN
        WRITE(0,1)"CONSTRAINT FLAGS MUST BE IN [0-1]"
        ERR=1
       ENDIF
      ENDDO
      IF(GMUFLAG.LT.0 .OR. GMUFLAG.GT.1)THEN
       WRITE(0,1)"GMUFLAG MUST BE IN [0-1]"
       ERR=1
      ENDIF
      IF(VFLAG.LT.0 .OR. VFLAG.GT.1)THEN
       WRITE(0,1)"VFLAG MUST BE IN [0-1]"
       ERR=1
      ENDIF
      IF(XDEV.EQ.1d99)THEN
       WRITE(0,1)"XDEV MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(LCEN.EQ.1d99)THEN
       WRITE(0,1)"RE(LAMBDA)CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ELSEIF(LCEN.LE.0d0)THEN
       WRITE(0,1)"RE(LAMBDA)CEN MUST BE STRICTLY POSITIVE"
       ERR=1
      ENDIF
      IF(TBCEN.EQ.1d99)THEN
       WRITE(0,1)"TBCEN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ELSEIF(TBCEN.LE.0d0)THEN
       WRITE(0,1)"TBCEN MUST BE STRICTLY POSITIVE"
       ERR=1
      ENDIF
      IF(MUCEN.EQ.1d99)THEN
       WRITE(0,1)"RE(MUEFF)CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ELSEIF(MUCEN.EQ.0d0)THEN
       WRITE(0,1)"RE(MUEFF)CEN MUST BE NON ZERO"
       ERR=1
      ENDIF
      IF(M1CEN.EQ.1d99 .AND. M1DEV.NE.1d99)THEN
       WRITE(0,1)"RE(M1)CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M2CEN.EQ.1d99)THEN
       WRITE(0,1)"RE(M2)CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M3CEN.EQ.1d99 .AND. M3DEV.NE.1d99)THEN
       WRITE(0,1)"RE(M3)CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(AQCEN.EQ.1d99)THEN
       WRITE(0,1)"RE(AU3)CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(REALP(8).EQ.1d99)THEN
       WRITE(0,1)"RE(AE3) MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MQCEN.EQ.1d99)THEN
       WRITE(0,1)"MQ3CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(10).EQ.1d99)THEN
       WRITE(0,1)"ML3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(11).EQ.1d99)THEN
       WRITE(0,1)"ME3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MACEN.EQ.1d99 .AND. MADEV.NE.1d99)THEN
       WRITE(0,1)"MACEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MPCEN.EQ.1d99 .AND. MPDEV.NE.1d99)THEN
       WRITE(0,1)"MPCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(ALCEN.EQ.1d99 .AND. ALDEV.NE.1d99)THEN
       WRITE(0,1)"RE(AL)CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(AKCEN.EQ.1d99 .AND. AKDEV.NE.1d99)THEN
       WRITE(0,1)"RE(AK)CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XIFCEN.EQ.1d99 .AND. XIFDEV.NE.1d99)THEN
       WRITE(0,1)"RE(XIF)CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XISCEN.EQ.1d99 .AND. XISDEV.NE.1d99)THEN
       WRITE(0,1)"RE(XIS)CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XISICEN.EQ.1d99 .AND. XISIDEV.NE.1d99)THEN
       WRITE(0,1)"IM(XIS)CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

      IF(PAR(18).EQ.1d99)PAR(18)=PAR(10)
      IF(PAR(19).EQ.1d99)PAR(19)=PAR(11)
      IF(PAR(25).EQ.1d99)PAR(25)=REALP(8)

*   Relations between (RE(ALAMBDA), MA, RE(XIF)) and (RE(AKAPPA), MP, RE(XIS))

      IF(ALCEN.NE.1d99 .AND. XIFCEN.NE.1d99 .AND. MACEN.NE.1d99)THEN
       WRITE(0,1)"AT MOST 2 OF THE 3 PARAMETERS ALAMBDA, MA AND XIF",
     . " CAN BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

      IF(AKCEN.NE.1d99 .AND. XISCEN.NE.1d99 .AND. MPCEN.NE.1d99)THEN
       WRITE(0,1)"AT MOST 2 OF THE 3 PARAMETERS AKAPPA, MP AND XIS",
     . " CAN BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

      IF(KCEN.EQ.0d0)THEN
       IF((AKCEN.NE.0d0 .AND. AKCEN.NE.1d99) .OR. (AKCEN.EQ.0d0
     . .AND. AKDEV.NE.0d0 .AND. AKMIN.NE.0d0))THEN
        WRITE(0,1)"IF KAPPA IS 0, AKAPPA MUST BE 0"
        ERR=1
       ELSE
        AKCEN=0d0
        AKDEV=0d0
        AKMIN=0d0
        IF(XISCEN.NE.1d99 .AND. MPCEN.NE.1d99)THEN
         WRITE(0,1)"IF KAPPA IS 0, EITHER MP OR XIS",
     .   " CAN BE GIVEN IN BLOCK EXTPAR"
         ERR=1
        ENDIF
       ENDIF
       IF(XISICEN.NE.1d99)THEN
        WRITE(0,1)"IF KAPPA IS 0, IM(XIS) IS NOT AN INPUT"
        ERR=1
       ENDIF
       XISIFLAG=0
      ENDIF

*   Set default values

      IF(ALCEN.EQ.1d99.AND.MACEN.EQ.1d99.AND.XIFCEN.EQ.1d99)THEN
       ALCEN=0d0
       XIFCEN=0d0
      ELSEIF(ALCEN.EQ.1d99.AND.MACEN.EQ.1d99)THEN
       ALCEN=0d0
      ELSEIF(ALCEN.EQ.1d99.AND.XIFCEN.EQ.1d99)THEN
       XIFCEN=0d0
      ELSEIF(MACEN.EQ.1d99.AND.XIFCEN.EQ.1d99)THEN
       XIFCEN=0d0
      ENDIF

      IF(AKCEN.EQ.1d99.AND.MPCEN.EQ.1d99.AND.XISCEN.EQ.1d99)THEN
       AKCEN=0d0
       XISCEN=0d0
      ELSEIF(AKCEN.EQ.1d99.AND.MPCEN.EQ.1d99)THEN
       AKCEN=0d0
      ELSEIF(AKCEN.EQ.1d99.AND.XISCEN.EQ.1d99)THEN
       XISCEN=0d0
      ELSEIF(MPCEN.EQ.1d99.AND.XISCEN.EQ.1d99)THEN
       XISCEN=0d0
      ENDIF

*   Set MAFLAG, SCANFLAGS

      IF(MACEN.EQ.1d99)MAFLAG=0
      IF(ALCEN.EQ.1d99)MAFLAG=1
      IF(XIFCEN.EQ.1d99)MAFLAG=2
      IF(AKCEN.EQ.1d99)MAFLAG=MAFLAG+3
      IF(XISCEN.EQ.1d99)MAFLAG=MAFLAG+6
      IF(M1CEN.NE.1d99)M1FLAG=1
      IF(M3CEN.NE.1d99)M3FLAG=1

      IF(XIFCEN.EQ.1d99)XIFCEN=0d0
      IF(XISCEN.EQ.1d99)XISCEN=0d0
      IF(XISICEN.EQ.1d99)XISICEN=0d0

*   Bounds

      IF(TBDEV.EQ.1d99)TBDEV=0d0
      IF(M1DEV.EQ.1d99)M1DEV=0d0
      IF(M2DEV.EQ.1d99)M2DEV=0d0
      IF(M3DEV.EQ.1d99)M3DEV=0d0
      IF(LDEV.EQ.1d99)LDEV=0d0
      IF(KDEV.EQ.1d99)KDEV=0d0
      IF(ALDEV.EQ.1d99)ALDEV=0d0
      IF(AKDEV.EQ.1d99)AKDEV=0d0
      IF(MUDEV.EQ.1d99)MUDEV=0d0
      IF(XIFDEV.EQ.1d99)XIFDEV=0d0
      IF(XISDEV.EQ.1d99)XISDEV=0d0
      IF(MUPDEV.EQ.1d99)MUPDEV=0d0
      IF(MSPDEV.EQ.1d99)MSPDEV=0d0
      IF(M3HDEV.EQ.1d99)M3HDEV=0d0
      IF(MADEV.EQ.1d99)MADEV=0d0
      IF(MPDEV.EQ.1d99)MPDEV=0d0
      IF(AQDEV.EQ.1d99)AQDEV=0d0
      IF(MQDEV.EQ.1d99)MQDEV=0d0
      IF(M1IDEV.EQ.1d99)M1IDEV=0d0
      IF(M2IDEV.EQ.1d99)M2IDEV=0d0
      IF(M3IDEV.EQ.1d99)M3IDEV=0d0
      IF(AU3IDEV.EQ.1d99)AU3IDEV=0d0
      IF(AD3IDEV.EQ.1d99)AD3IDEV=0d0
      IF(AE3IDEV.EQ.1d99)AE3IDEV=0d0
      IF(LIDEV.EQ.1d99)LIDEV=0d0
      IF(KIDEV.EQ.1d99)KIDEV=0d0
      IF(XIFIDEV.EQ.1d99)XIFIDEV=0d0
      IF(XISIDEV.EQ.1d99)XISIDEV=0d0
      IF(MUPIDEV.EQ.1d99)MUPIDEV=0d0
      IF(MSPIDEV.EQ.1d99)MSPIDEV=0d0
      IF(M3HIDEV.EQ.1d99)M3HIDEV=0d0

*   Total number of points

      IF(NTOT.LE.0)THEN
       WRITE(0,1)"WRONG NUMBER OF POINTS IN BLOCK STEPS"
       ERR=1
      ENDIF
      IF(TOTMIN.GT.TOTMAX)THEN
       WRITE(0,1)"TOTMIN must be smaller than TOTMAX"
       ERR=1
      ENDIF

*   Check for Z3 breaking terms

      IF(MOD(MAFLAG,3).EQ.2 .OR. MAFLAG/3.EQ.2
     ..OR. MUPCEN.NE.0d0 .OR. MUPDEV.NE.0d0
     ..OR. MSPCEN.NE.0d0 .OR. MSPDEV.NE.0d0
     ..OR. XIFCEN.NE.0d0 .OR. XIFDEV.NE.0d0
     ..OR. XISCEN.NE.0d0 .OR. XISDEV.NE.0d0
     ..OR. M3HCEN.NE.0d0 .OR. M3HDEV.NE.0d0
     ..OR. MUPICEN.NE.0d0 .OR. MUPIDEV.NE.0d0
     ..OR. MSPICEN.NE.0d0 .OR. MSPIDEV.NE.0d0
     ..OR. XIFICEN.NE.0d0 .OR. XIFIDEV.NE.0d0
     ..OR. XISICEN.NE.0d0 .OR. XISIDEV.NE.0d0
     ..OR. M3HICEN.NE.0d0 .OR. M3HIDEV.NE.0d0
     ..OR. (KCEN.EQ.0d0 .AND. KICEN.EQ.0d0))THEN
       IF(Z3FLAG.GT.2)THEN
        WRITE(0,1)"PRESENCE OF Z3 BREAKING TERMS"
        ERR=1
       ENDIF
      ENDIF

*   Stop if error

      IF(ERR.EQ.1)THEN
       WRITE(0,1)"ERROR IN INPUT FILE"
       STOP 1
      ENDIF

*   Set Q2MIN, Q2FIX:
      Q2MIN=100d0**2
      Q2FIX=1
      IF(Q2.LE.Q2MIN)THEN
       Q2FIX=0
      ENDIF

*   Initialization for ALPHAS and RUNM (as in hdecay)
*   The bottom quark pole mass MBP is set in INIT and can be changed
*   only there (changing its running mass MB above has no effect
*   on MBP, since one would have to compute alpha_s(MB) first)

      MC0=MC
      MB0=MBP
      MT0=MT
      N0=5
      NLOOP=3
      NBER=18
      ACC=1d-10
      XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
      CALL ALSINI(ACC)
      CALL BERNINI(NBER)

*    g1,g2  and sin(theta)^2 in the on-shell scheme in terms of
*    GF, MZ(pole) and MW(pole)

      g2=4d0*DSQRT(2d0)*GF*MW**2
      g1=4d0*DSQRT(2d0)*GF*(MZ**2-MW**2)
      S2TW=1d0-(MW/MZ)**2

      RETURN

 999  WRITE(0,1)"READ ERROR ON LINE:", NLINE
      WRITE(0,*)CHINL(1:80)
      STOP 1

 1    FORMAT(2A)

      END


      SUBROUTINE OUTPUT(PAR,PROB,IFAIL)

*********************************************************************
*   Subroutine writing all the results in the the output files.
*********************************************************************

      IMPLICIT NONE

      INTEGER I,IFAIL,IMAX,IRES,NRES,Q2FIX
      PARAMETER (IMAX=200)

      DOUBLE PRECISION RES(IMAX),PAR(*),PROB(*)
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
      DOUBLE PRECISION CU(5),CUP(5),CD(5),CDP(5),CB(5),CBP(5),CJ(5)
      DOUBLE PRECISION CJP(5),CI(5),CG(5),CGP(5),CV(5),CZG(5),CZGP(5)
      DOUBLE PRECISION CL(5),CLP(5)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,Q2,Q2MIN
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MHUS,MHDS,MSS
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION MHUQ,MHDQ,MSX,LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,QSTSB
      DOUBLE PRECISION BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      DOUBLE PRECISION BRBSll,BRBSllmin,BRBSllMax,
     .             BRBShll,BRBShllmin,BRBShllMax,
     .             BRDG,BRDGmin,BRDGmax,
     .             BRBdMUMU,BRBdMUMUmin,BRBdMUMUmax,
     .             BRBXsnunu,BRBXsnunumin,BRBXsnunumax,
     .             BRBpKpnunu,BRBpKpnunumin,BRBpKpnunuMax,
     .             BRBKsnunu,BRBKsnunumin,BRBKsnunuMax,
     .             RD_taul,RD_taulmin,RD_taulmax,
     .             RDs_taul,RDs_taulmin,RDs_taulmax
      DOUBLE PRECISION delmagmu,damumin,damumax,amuthmax,amuthmin
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2),toptot
      DOUBLE PRECISION MHmin,MHmax
      DOUBLE PRECISION MHC_CPV,XC(2,2),MH0T(5),XH0(5,5),MA2
      DOUBLE PRECISION MCH2(2),U_CPV(2,2,2),V_CPV(2,2,2)
      DOUBLE PRECISION MNEU_CPV(5),NEU_CPV(5,5,2)
      DOUBLE PRECISION MGL_CPV,SIG_CPV(5,8),DELMB,DELML,DEL1
      DOUBLE PRECISION dEe,dETl,dEnCQM,dEnPQM,dEnQSR,dEHg
      DOUBLE PRECISION dEemin,dETlmin,dEnCQMmin,dEnPQMmin,dEnQSRmin,
     . dEHgmin
      DOUBLE PRECISION MST2_CPV(2),UT(2,2,2),MSB2_CPV(2),
     .      UB(2,2,2),MSL2_CPV(2),UTAU(2,2,2),MSNT2_CPV
      DOUBLE PRECISION MSU2_CPV(2),MSD2_CPV(2),MSE2_CPV(2),MSNE2_CPV,
     .      MSMU2_CPV(2),UMU(2,2,2)
      DOUBLE PRECISION MST2P(2),MSB2P(2),MSU2P(2),MSD2P(2)
      DOUBLE PRECISION MSL3_CPV,MSE3_CPV,MSL1_CPV,MSE1_CPV,ATAU_CPV,
     .    AMU_CPV
      DOUBLE PRECISION MYPHASES(16),REALP(14),IMAGP(14)
      DOUBLE PRECISION RXIS,RXIF,IAlQ2,IAkQ2,IXISQ2,MHuQ2,MHdQ2,MSQ2

      COMMON/HIWIDTH/WIDTH,HCWIDTH
      COMMON/HNSMBR/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG
      COMMON/HNHIBR/BRHHH,BRHCHC,BRHAZ,BRHCW,BRHIGGS
      COMMON/HNSUSYBR/BRNEU,BRCHA,BRHSQ,BRHSL,BRSUSY
      COMMON/HCSMBR/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,HCBRBT
      COMMON/HCHIBR/HCBRWH,HCBRWHT
      COMMON/HCSUSYBR/HCBRNC,HCBRSQ,HCBRSL,HCBRSUSY
      COMMON/BRSG/BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      COMMON/FLAV2/BRBSll,BRBSllmin,BRBSllMax,
     .      BRBShll,BRBShllmin,BRBShllMax,
     .      BRDG,BRDGmin,BRDGmax,
     .      BRBdMUMU,BRBdMUMUmin,BRBdMUMUmax,
     .      BRBXsnunu,BRBXsnunumin,BRBXsnunumax,
     .      BRBpKpnunu,BRBpKpnunumin,BRBpKpnunuMax,
     .      BRBKsnunu,BRBKsnunumin,BRBKsnunuMax,
     .      RD_taul,RD_taulmin,RD_taulmax,
     .      RDs_taul,RDs_taulmin,RDs_taulmax
      COMMON/MAGMU/delmagmu,damumin,damumax,amuthmax,amuthmin
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/topwidth/toptot
      COMMON/HNSMCOUP/CU,CUP,CD,CDP,CB,CBP,CJ,CJP,CI,CG,CGP,CV,CZG,CZGP,
     . CL,CLP
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SUSYMH/MHUS,MHDS,MSS
      COMMON/QMHIGGS/MHUQ,MHDQ,MSX
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/HIGGSFIT/MHmin,MHmax
      COMMON/HISPEC/MHC_CPV,XC,MH0T,XH0,MA2
      COMMON/CHASPEC/MCH2,U_CPV,V_CPV
      COMMON/NEUSPEC/MNEU_CPV,NEU_CPV
      COMMON/GLUSPEC/MGL_CPV
      COMMON/LHCSIGCPV/SIG_CPV
      COMMON/EDM/dEe,dETl,dEnCQM,dEnPQM,dEnQSR,dEHg,
     . dEemin,dETlmin,dEnCQMmin,dEnPQMmin,dEnQSRmin,dEHgmin
      COMMON/SFERM3SPEC/MST2_CPV,UT,MSB2_CPV,UB,MSL2_CPV,UTAU,MSNT2_CPV
      COMMON/SFERM1SPEC/MSU2_CPV,MSD2_CPV,MSE2_CPV,MSNE2_CPV,MSMU2_CPV,
     .     UMU
      COMMON/SFERMPSPEC/MST2P,MSB2P,MSU2P,MSD2P
      COMMON/SLEPPAR/MSL3_CPV,MSE3_CPV,MSL1_CPV,MSE1_CPV,ATAU_CPV,
     .     AMU_CPV
      COMMON/MYPHASES/MYPHASES
      COMMON/REAL_IMAG/REALP,IMAGP
      COMMON/REXI/RXIS,RXIF,IAlQ2,IAkQ2,IXISQ2,MHuQ2,MHdQ2,MSQ2

      IF(IFAIL.NE.0)RETURN

      IRES=28
      NRES=IRES+19

      RES(1)=REALP(1)
      RES(2)=REALP(2)
      RES(3)=PAR(3)
      RES(4)=REALP(14)
      RES(5)=REALP(4)
      RES(6)=PAR(5)
      RES(7)=PAR(6)
      RES(8)=REALP(9)
      RES(9)=REALP(10)
      RES(10)=PAR(23)
      RES(11)=PAR(24)
      RES(12)=REALP(6)
      RES(13)=DSQRT(PAR(7))
      DO I=1,13
       IF(I.EQ.10)THEN
        RES(13+I)=IXISQ2
       ELSE
        RES(13+I)=IMAGP(I)
       ENDIF
      ENDDO
      RES(27)=IALQ2
      RES(28)=IAKQ2

      DO I=1,5
       RES(IRES+I)=DSQRT(MH0T(I))
      ENDDO
      DO I=1,5
       RES(IRES+5+I)=DSQRT(MNEU_CPV(I))
      ENDDO
      DO I=1,2
       RES(IRES+10+I)=DSQRT(MCH2(I))
      ENDDO
      DO I=1,2
       RES(IRES+12+I)=DSQRT(MST2P(I))
      ENDDO
      DO I=1,2
       RES(IRES+14+I)=DSQRT(MSB2P(I))
      ENDDO
      DO I=1,2
       RES(IRES+16+I)=DSQRT(MSL2_CPV(I))
      ENDDO
      RES(IRES+19)=DSQRT(MSNT2_CPV)

      WRITE(6,11)(RES(I),I=1,NRES)
 11   FORMAT(200E14.6)

      END


      SUBROUTINE ERROR(TOT,NTOT,NFAIL)

*********************************************************************
*   Subroutine for the error file. It contains a summary of the scan:
*   Number of points that passed/failed the tests
*   and ranges for scanned parameters that passed the tests
*********************************************************************

      IMPLICIT NONE

      INTEGER I,S,TOT,NTOT,NFAIL(*),GMUFLAG,HFLAG,CFLAG(5)
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,MSFLAG
      INTEGER AKFLAG,ALFLAG,XISIFLAG,OMGFLAG,MAFLAG,MOFLAG

      DOUBLE PRECISION LN,LNN,KN,KNN,TBN,TBNN,MUN,MUNN
      DOUBLE PRECISION ALN,ALNN,AKN,AKNN,XIFN,XIFNN
      DOUBLE PRECISION XISN,XISNN,MUPN,MUPNN,MSPN,MSPNN
      DOUBLE PRECISION M3HN,M3HNN,MAN,MANN,MPN,MPNN
      DOUBLE PRECISION M1N,M1NN,M2N,M2NN,M3N,M3NN
      DOUBLE PRECISION AQN,AQNN,MQNN,MQN,DEV
      DOUBLE PRECISION M1IN,M1INN,M2IN,M2INN,M3IN,M3INN
      DOUBLE PRECISION AU3IN,AU3INN,AD3IN,AD3INN,AE3IN,AE3INN
      DOUBLE PRECISION LIN,LINN,KIN,KINN,XIFIN,XIFINN,XISIN,XISINN
      DOUBLE PRECISION MUPIN,MUPINN,MSPIN,MSPINN,M3HIN,M3HINN
      DOUBLE PRECISION ALIN,ALINN,AKIN,AKINN

      COMMON/BOUNDS/LN,LNN,KN,KNN,TBN,TBNN,MUN,MUNN,
     . ALN,ALNN,AKN,AKNN,XIFN,XIFNN,
     . XISN,XISNN,MUPN,MUPNN,MSPN,MSPNN,
     . M3HN,M3HNN,MAN,MANN,MPN,MPNN,
     . M1N,M1NN,M2N,M2NN,M3N,M3NN,
     . AQN,AQNN,MQN,MQNN
      COMMON/IBOUNDS/M1IN,M1INN,M2IN,M2INN,M3IN,M3INN,
     . AU3IN,AU3INN,AD3IN,AD3INN,AE3IN,AE3INN,LIN,LINN,
     . KIN,KINN,XIFIN,XIFINN,XISIN,XISINN,MUPIN,MUPINN,
     . MSPIN,MSPINN,M3HIN,M3HINN,ALIN,ALINN,AKIN,AKINN
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . MSFLAG,AKFLAG,ALFLAG,XISIFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/CFLAG/CFLAG

      WRITE(0,10)"NMSSMTools scan info          "
      WRITE(0,10)"Version number: 5.6.2         "
      WRITE(0,*)
      WRITE(0,*)
      WRITE(0,10)"Number of points:             "
      WRITE(0,*)
      WRITE(0,10)"  scanned                     ",NTOT
      WRITE(0,10)"  mu = 0 | k=0 & (Ak|XiS)=/=0 ",NFAIL(9)
      S=0
      DO I=1,7
       S=S+NFAIL(I)
      ENDDO
      WRITE(0,10)"  with mh1^2 or mhc^2 < 0     ",S
      WRITE(0,10)"  with m_sfermion^2 < 0       ",NFAIL(8)
      WRITE(0,10)"  violating constraints       ",NFAIL(10)
      WRITE(0,*)
      WRITE(0,10)"Remaining good points         ",TOT
      WRITE(0,*)
      WRITE(0,*)

      IF(CFLAG(1).EQ.0 .AND. CFLAG(2).EQ.0 .AND. CFLAG(3).EQ.0 .AND.
     .   CFLAG(4).EQ.0 .AND. CFLAG(5).EQ.0 .AND. GMUFLAG.EQ.0)THEN
       WRITE(0,20)"Contraints taken into account: none"
      ELSE
       WRITE(0,20)"Contraints taken into account:     "
       IF(CFLAG(1).EQ.1)
     .  WRITE(0,20)" - Landau poles and false minima   "
       IF(CFLAG(2).EQ.1)
     .  WRITE(0,20)" - LEP/Tevatron Higgs+sparticle    "
       IF(CFLAG(3).EQ.1)
     .  WRITE(0,20)" - LHC Higgs                       "
       IF(CFLAG(4).EQ.1)
     .  WRITE(0,20)" - Upsilon, B and K decays         "
       IF(CFLAG(5).EQ.1)
     .  WRITE(0,20)" - EDMs                            "
       IF(GMUFLAG.EQ.1)
     .  WRITE(0,20)" - (g-2)_muon                      "
      ENDIF

      IF(TOT.GT.0)THEN

       WRITE(0,*)
       WRITE(0,*)
       WRITE(0,20)"Parameter ranges for good points:  "
       WRITE(0,*)
       WRITE(0,30)" TANB: ",TBN,TBNN,DEV(TBN,TBNN)
       IF(M1FLAG.EQ.1)THEN
        WRITE(0,30)" RE(M1): ",M1N,M1NN,DEV(M1N,M1NN)
       ENDIF
       WRITE(0,30)" RE(M2): ",M2N,M2NN,DEV(M2N,M2NN)
       IF(M3FLAG.EQ.1)THEN
        WRITE(0,30)" RE(M3): ",M3N,M3NN,DEV(M3N,M3NN)
       ENDIF
       WRITE(0,30)" RE(AU3): ",AQN,AQNN,DEV(AQN,AQNN)
       WRITE(0,30)" MQ3: ",MQN,MQNN,DEV(MQN,MQNN)
       WRITE(0,30)" RE(LAMBDA): ",LN,LNN,DEV(LN,LNN)
       WRITE(0,30)" RE(KAPPA): ",KN,KNN,DEV(KN,KNN)
       WRITE(0,30)" RE(MUEFF): ",MUN,MUNN,DEV(MUN,MUNN)
       IF(MOD(MAFLAG,3).EQ.0)THEN
        WRITE(0,30)" RE(ALAMBDA): ",ALN,ALNN,DEV(ALN,ALNN)
        WRITE(0,30)" MA: ",MAN,MANN
        WRITE(0,40)"(MA is not an input parameter)"
        IF(XIFN.NE.0d0 .OR. XIFNN.NE.0d0)
     .   WRITE(0,30)" RE(XIF): ",XIFN,XIFNN,DEV(XIFN,XIFNN)
       ELSEIF(MOD(MAFLAG,3).EQ.1)THEN
        WRITE(0,30)" RE(ALAMBDA): ",ALN,ALNN
        WRITE(0,40)"(RE(ALAMBDA) is not an input parameter)"
        WRITE(0,30)" MA: ",MAN,MANN,DEV(MAN,MANN)
        IF(XIFN.NE.0d0 .OR. XIFNN.NE.0d0)
     .   WRITE(0,30)"RE(XIF): ",XIFN,XIFNN,DEV(MAN,MANN)
       ELSE
        WRITE(0,30)" RE(ALAMBDA): ",ALN,ALNN,DEV(ALN,ALNN)
        WRITE(0,30)" MA: ",MAN,MANN,DEV(MAN,MANN)
        WRITE(0,30)" RE(XIF): ",XIFN,XIFNN
        WRITE(0,40)"(RE(XIF) is not an input parameter)"
       ENDIF
       IF(MAFLAG/3.EQ.0)THEN
        WRITE(0,30)" RE(AKAPPA): ",AKN,AKNN,DEV(AKN,AKNN)
        WRITE(0,30)" MP: ",MPN,MPNN
        WRITE(0,40)"(MP is not an input parameter)"
        IF(XISN.NE.0d0 .OR. XISNN.NE.0d0)
     .   WRITE(0,30)" RE(XIS): ",XISN,XISNN,DEV(XISN,XISNN)
       ELSEIF(MAFLAG/3.EQ.1)THEN
        WRITE(0,30)" RE(AKAPPA): ",AKN,AKNN
        WRITE(0,40)"(RE(AKAPPA) is not an input parameter)"
        WRITE(0,30)" MP: ",MPN,MPNN,DEV(MPN,MPNN)
        IF(XISN.NE.0d0 .OR. XISNN.NE.0d0)
     .   WRITE(0,30)" RE(XIS): ",XISN,XISNN,DEV(XISN,XISNN)
       ELSE
        WRITE(0,30)" RE(AKAPPA): ",AKN,AKNN,DEV(AKN,AKNN)
        WRITE(0,30)" MP: ",MPN,MPNN,DEV(MPN,MPNN)
        WRITE(0,30)" RE(XIS): ",XISN,XISNN
        WRITE(0,40)"(RE(XIS) is not an input parameter)"
       ENDIF
       IF(MUPN.NE.0d0 .OR. MUPNN.NE.0d0)
     .  WRITE(0,30)" RE(MUP): ",MUPN,MUPNN,DEV(MUPN,MUPNN)
       IF(MSPN.NE.0d0 .OR. MSPNN.NE.0d0)
     . WRITE(0,30)" RE(MSP): ",MSPN,MSPNN,DEV(MSPN,MSPNN)
       IF(M3HN.NE.0d0 .OR. M3HNN.NE.0d0)
     .  WRITE(0,30)" RE(M3H): ",M3HN,M3HNN,DEV(M3HN,M3HNN)
       IF(M1IN.NE.0d0 .OR. M1INN.NE.0d0)
     .  WRITE(0,30)" IM(M1): ",M1IN,M1INN,DEV(M1IN,M1INN)
       IF(M2IN.NE.0d0 .OR. M2INN.NE.0d0)
     .  WRITE(0,30)" IM(M2): ",M2IN,M2INN,DEV(M2IN,M2INN)
       IF(M3IN.NE.0d0 .OR. M3INN.NE.0d0)
     .  WRITE(0,30)" IM(M3): ",M3IN,M3INN,DEV(M3IN,M3INN)
       IF(AU3IN.NE.0d0 .OR. AU3INN.NE.0d0)
     .  WRITE(0,30)" IM(AU3): ",AU3IN,AU3INN,DEV(AU3IN,AU3INN)
       IF(AD3IN.NE.0d0 .OR. AD3INN.NE.0d0)
     .  WRITE(0,30)" IM(AD3): ",AD3IN,AD3INN,DEV(AD3IN,AD3INN)
       IF(AE3IN.NE.0d0 .OR. AE3INN.NE.0d0)
     .  WRITE(0,30)" IM(AE3): ",AE3IN,AE3INN,DEV(AE3IN,AE3INN)
       IF(LIN.NE.0d0 .OR.LINN.NE.0d0)
     .  WRITE(0,30)" IM(LAMBDA): ",LIN,LINN,DEV(LIN,LINN)
       IF(KIN.NE.0d0 .OR. KINN.NE.0d0)
     .  WRITE(0,30)" IM(KAPPA): ",KIN,KINN,DEV(KIN,KINN)
       IF(XIFIN.NE.0d0 .OR. XIFINN.NE.0d0)
     .  WRITE(0,30)" IM(XIF): ",XIFIN,XIFINN,DEV(XIFIN,XIFINN)
       IF(XISIN.NE.0d0 .OR. XISINN.NE.0d0)THEN
        WRITE(0,30)" IM(XIS): ",XISIN,XISINN,DEV(XISIN,XISINN)
        IF(XISIFLAG.EQ.0)
     .   WRITE(0,40)"(IM(XIS) is not an input parameter)"
       ENDIF
       IF(MUPIN.NE.0d0 .OR. MUPINN.NE.0d0)
     .  WRITE(0,30)" IM(MUP): ",MUPIN,MUPINN,DEV(MUPIN,MUPINN)
       IF(MSPIN.NE.0d0 .OR. MSPINN.NE.0d0)
     .  WRITE(0,30)" IM(MSP): ",MSPIN,MSPINN,DEV(MSPIN,MSPINN)
       IF(M3HIN.NE.0d0 .OR. M3HINN.NE.0d0)
     .  WRITE(0,30)" IM(M3H): ",M3HIN,M3HINN,DEV(M3HIN,M3HINN)
       IF(ALIN.NE.0d0 .OR. ALINN.NE.0d0)THEN
        WRITE(0,30)" IM(AL): ",ALIN,ALINN,DEV(ALIN,ALINN)
        WRITE(0,40)"(IM(AL) is not an input parameter)"
       ENDIF
       IF(AKIN.NE.0d0 .OR. AKINN.NE.0d0)THEN
        WRITE(0,30)" IM(AK): ",AKIN,AKINN,DEV(AKIN,AKINN)
        WRITE(0,40)"(IM(AK) is not an input parameter)"
       ENDIF

      ENDIF

 10   FORMAT(A30,I10)
 20   FORMAT(A35)
 30   FORMAT(A14,3E15.4)
 40   FORMAT(A39)

      END
