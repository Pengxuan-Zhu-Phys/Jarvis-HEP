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
*      MGl     Gluino mass
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
      INTEGER NFAIL(NFL),IFAIL,I,TOT,ITOT,NTOT
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG
      INTEGER MSFLAG,AKFLAG,ALFLAG,XISIFLAG
      INTEGER OMGFLAG,MAFLAG,MOFLAG,GMUFLAG,HFLAG
      INTEGER NL,NK,NTB,NMU,NMUP,NMSP,NM3H
      INTEGER N1,N2,N3,N4,NM1,NM2,NM3,NAQ,NMQ
      INTEGER NLI,NKI,NM1I,NM2I,NM3I,NAU3I,NAD3I
      INTEGER NAE3I,NXIFI,NXISI,NMUPI,NMSPI,NM3HI
      INTEGER IL,IK,ITB,IMU,IMUP,IMSP,IM3H
      INTEGER I1,I2,I3,I4,IM1,IM2,IM3,IAQ,IMQ
      INTEGER ILI,IKI,IM1I,IM2I,IM3I,IAU3I,IAD3I
      INTEGER IAE3I,IXIFI,IXISI,IMUPI,IMSPI,IM3HI
      INTEGER CFLAG(5)

      DOUBLE PRECISION PAR(NPAR),PROB(NPROB),Q2,PI
      DOUBLE PRECISION LMIN,LMAX,KMIN,KMAX,TBMIN,TBMAX,MUMIN,MUMAX
      DOUBLE PRECISION ALMIN,ALMAX,AKMIN,AKMAX,XIFMIN,XIFMAX
      DOUBLE PRECISION XISMIN,XISMAX,MUPMIN,MUPMAX,MSPMIN,MSPMAX
      DOUBLE PRECISION M3HMIN,M3HMAX,MAMIN,MAMAX,MPMIN,MPMAX
      DOUBLE PRECISION M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX
      DOUBLE PRECISION AQMIN,AQMAX,MQMIN,MQMAX
      DOUBLE PRECISION M1IMIN,M1IMAX,M2IMIN,M2IMAX,M3IMIN,M3IMAX
      DOUBLE PRECISION AU3IMIN,AU3IMAX,AD3IMIN,AD3IMAX,AE3IMIN
      DOUBLE PRECISION AE3IMAX,LIMIN,LIMAX,KIMIN,KIMAX,XIFIMIN
      DOUBLE PRECISION XIFIMAX,XISIMIN,XISIMAX,MUPIMIN,MUPIMAX
      DOUBLE PRECISION MSPIMIN,MSPIMAX,M3HIMIN,M3HIMAX
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
      COMMON/STEPS/NTOT,NL,NK,NTB,NMU,NMUP,NMSP,NM3H,
     . N1,N2,N3,N4,NM1,NM2,NM3,NAQ,NMQ,
     . NLI,NKI,NM1I,NM2I,NM3I,NAU3I,NAD3I,
     . NAE3I,NXIFI,NXISI,NMUPI,NMSPI,NM3HI
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
      COMMON/MINMAX/LMIN,LMAX,KMIN,KMAX,TBMIN,TBMAX,MUMIN,MUMAX,
     . ALMIN,ALMAX,AKMIN,AKMAX,XIFMIN,XIFMAX,XISMIN,XISMAX,
     . MUPMIN,MUPMAX,MSPMIN,MSPMAX,M3HMIN,M3HMAX,MAMIN,MAMAX,
     . MPMIN,MPMAX,M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX,
     . AQMIN,AQMAX,MQMIN,MQMAX
      COMMON/MINIMAX/M1IMIN,M1IMAX,M2IMIN,M2IMAX,M3IMIN,M3IMAX,
     . AU3IMIN,AU3IMAX,AD3IMIN,AD3IMAX,AE3IMIN,AE3IMAX,LIMIN,
     . LIMAX,KIMIN,KIMAX,XIFIMIN,XIFIMAX,XISIMIN,XISIMAX,MUPIMIN,
     . MUPIMAX,MSPIMIN,MSPIMAX,M3HIMIN,M3HIMAX
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
      ITOT=0

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

      DO ITB=1,NTB
      IF(NTB.EQ.1)THEN
       PAR(3)=TBMIN
      ELSE
       PAR(3)=TBMIN+(TBMAX-TBMIN)*DFLOAT(ITB-1)/DFLOAT(NTB-1)
      ENDIF

      DO IM2=1,NM2
      IF(NM2.EQ.1)THEN
       REALP(4)=M2MIN
      ELSE
       REALP(4)=M2MIN+(M2MAX-M2MIN)*DFLOAT(IM2-1)/DFLOAT(NM2-1)
      ENDIF

      DO IM1=1,NM1
      IF(M1FLAG.EQ.0)THEN
       REALP(3)=REALP(4)/2d0
      ELSEIF(NM1.EQ.1)THEN
       REALP(3)=M1MIN
      ELSE
       REALP(3)=M1MIN+(M1MAX-M1MIN)*DFLOAT(IM1-1)/DFLOAT(NM1-1)
      ENDIF

      DO IM3=1,NM3
      IF(M3FLAG.EQ.0)THEN
       REALP(5)=REALP(4)*3d0
      ELSEIF(NM3.EQ.1)THEN
       REALP(5)=M3MIN
      ELSE
       REALP(5)=M3MIN+(M3MAX-M3MIN)*DFLOAT(IM3-1)/DFLOAT(NM3-1)
      ENDIF

      DO IL=1,NL
      IF(NL.EQ.1)THEN
       REALP(1)=LMIN
      ELSE
       REALP(1)=LMIN+(LMAX-LMIN)*DFLOAT(IL-1)/DFLOAT(NL-1)
      ENDIF

      DO IK=1,NK
      IF(NK.EQ.1)THEN
       REALP(2)=KMIN
      ELSE
       REALP(2)=KMIN+(KMAX-KMIN)*DFLOAT(IK-1)/DFLOAT(NK-1)
      ENDIF

      DO IMU=1,NMU
      IF(NMU.EQ.1)THEN
       REALP(14)=MUMIN
      ELSE
       REALP(14)=MUMIN+(MUMAX-MUMIN)*DFLOAT(IMU-1)/DFLOAT(NMU-1)
      ENDIF

      DO I1=1,N1
      IF(N1.EQ.1)THEN
       IF(MOD(MAFLAG,3).NE.1)THEN
        PAR(5)=ALMIN
       ELSE
        PAR(23)=MAMIN
       ENDIF
      ELSE
       IF(MOD(MAFLAG,3).NE.1)THEN
        PAR(5)=ALMIN+(ALMAX-ALMIN)*DFLOAT(I1-1)/DFLOAT(N1-1)
       ELSE
        PAR(23)=MAMIN+(MAMAX-MAMIN)*DFLOAT(I1-1)/DFLOAT(N1-1)
       ENDIF
      ENDIF

      DO I2=1,N2
      IF(N2.EQ.1)THEN
       IF(MOD(MAFLAG,3).NE.2)THEN
        REALP(9)=XIFMIN
       ELSE
        PAR(23)=MAMIN
        REALP(9)=0d0
       ENDIF
      ELSE
       IF(MOD(MAFLAG,3).NE.2)THEN
        REALP(9)=XIFMIN+(XIFMAX-XIFMIN)*DFLOAT(I2-1)/DFLOAT(N2-1)
       ELSE
        PAR(23)=MAMIN+(MAMAX-MAMIN)*DFLOAT(I2-1)/DFLOAT(N2-1)
        REALP(9)=0d0
       ENDIF
      ENDIF

      DO I3=1,N3
      IF(N3.EQ.1)THEN
       IF(MAFLAG/3.NE.1)THEN
        PAR(6)=AKMIN
       ELSE
        PAR(24)=MPMIN
       ENDIF
      ELSE
       IF(MAFLAG/3.NE.1)THEN
        PAR(6)=AKMIN+(AKMAX-AKMIN)*DFLOAT(I3-1)/DFLOAT(N3-1)
       ELSE
        PAR(24)=MPMIN+(MPMAX-MPMIN)*DFLOAT(I3-1)/DFLOAT(N3-1)
       ENDIF
      ENDIF

      DO I4=1,N4
      IF(N4.EQ.1)THEN
       IF(MAFLAG/3.NE.2)THEN
        REALP(10)=XISMIN
       ELSE
        PAR(24)=MPMIN
        REALP(10)=0d0
       ENDIF
      ELSE
       IF(MAFLAG/3.NE.2)THEN
        REALP(10)=XISMIN+(XISMAX-XISMIN)*DFLOAT(I4-1)/DFLOAT(N4-1)
       ELSE
        PAR(24)=MPMIN+(MPMAX-MPMIN)*DFLOAT(I4-1)/DFLOAT(N4-1)
        REALP(10)=0d0
       ENDIF
      ENDIF

      DO IMUP=1,NMUP
      IF(NMUP.EQ.1)THEN
       REALP(11)=MUPMIN
      ELSE
       REALP(11)=MUPMIN+(MUPMAX-MUPMIN)*DFLOAT(IMUP-1)/DFLOAT(NMUP-1)
      ENDIF

      DO IMSP=1,NMSP
      IF(NMSP.EQ.1)THEN
       REALP(12)=MSPMIN
      ELSE
       REALP(12)=MSPMIN+(MSPMAX-MSPMIN)*DFLOAT(IMSP-1)/DFLOAT(NMSP-1)
      ENDIF

      DO IM3H=1,NM3H
      IF(NM3H.EQ.1)THEN
       REALP(13)=M3HMIN
      ELSE
       REALP(13)=M3HMIN+(M3HMAX-M3HMIN)*DFLOAT(IM3H-1)/DFLOAT(NM3H-1)
      ENDIF

      DO IAQ=1,NAQ
      IF(NAQ.EQ.1)THEN
       REALP(6)=AQMIN
      ELSE
       REALP(6)=AQMIN+(AQMAX-AQMIN)*DFLOAT(IAQ-1)/DFLOAT(NAQ-1)
      ENDIF

      DO IMQ=1,NMQ
      IF(NMQ.EQ.1)THEN
       PAR(7)=MQMIN**2
      ELSE
       PAR(7)=(MQMIN+(MQMAX-MQMIN)*DFLOAT(IMQ-1)/DFLOAT(NMQ-1))**2
      ENDIF

      DO ILI=1,NLI
      IF(NLI.EQ.1)THEN
       IMAGP(1)=LIMIN
      ELSE
       IMAGP(1)=LIMIN+(LIMAX-LIMIN)*DFLOAT(ILI-1)/DFLOAT(NLI-1)
      ENDIF

      DO IKI=1,NKI
      IF(NKI.EQ.1)THEN
       IMAGP(2)=KIMIN
      ELSE
       IMAGP(2)=KIMIN+(KIMAX-KIMIN)*DFLOAT(IKI-1)/DFLOAT(NKI-1)
      ENDIF

      DO IM1I=1,NM1I
      IF(NM1I.EQ.1)THEN
       IMAGP(3)=M1IMIN
      ELSE
       IMAGP(3)=M1IMIN+(M1IMAX-M1IMIN)*DFLOAT(IM1I-1)/DFLOAT(NM1I-1)
      ENDIF

      DO IM2I=1,NM2I
      IF(NM2I.EQ.1)THEN
       IMAGP(4)=M2IMIN
      ELSE
       IMAGP(4)=M2IMIN+(M2IMAX-M2IMIN)*DFLOAT(IM2I-1)/DFLOAT(NM2I-1)
      ENDIF

      DO IM3I=1,NM3I
      IF(NM3I.EQ.1)THEN
       IMAGP(5)=M3IMIN
      ELSE
       IMAGP(5)=M3IMIN+(M3IMAX-M3IMIN)*DFLOAT(IM3I-1)/DFLOAT(NM3I-1)
      ENDIF

      DO IAU3I=1,NAU3I
      IF(NAU3I.EQ.1)THEN
       IMAGP(6)=AU3IMIN
      ELSE
       IMAGP(6)=AU3IMIN
     .       +(AU3IMAX-AU3IMIN)*DFLOAT(IAU3I-1)/DFLOAT(NAU3I-1)
      ENDIF

      DO IAD3I=1,NAD3I
      IF(NAD3I.EQ.1)THEN
       IMAGP(7)=AD3IMIN
      ELSE
       IMAGP(7)=AD3IMIN
     .        +(AD3IMAX-AD3IMIN)*DFLOAT(IAD3I-1)/DFLOAT(NAD3I-1)
      ENDIF

      DO IAE3I=1,NAE3I
      IF(NAE3I.EQ.1)THEN
       IMAGP(8)=AE3IMIN
      ELSE
       IMAGP(8)=AE3IMIN
     .        +(AE3IMAX-AE3IMIN)*DFLOAT(IAE3I-1)/DFLOAT(NAE3I-1)
      ENDIF

      DO IXIFI=1,NXIFI
      IF(NXIFI.EQ.1)THEN
       IMAGP(9)=XIFIMIN
      ELSE
       IMAGP(9)=XIFIMIN
     .        +(XIFIMAX-XIFIMIN)*DFLOAT(IXIFI-1)/DFLOAT(NXIFI-1)
      ENDIF

      DO IXISI=1,NXISI
      IF(NXISI.EQ.1)THEN
       IMAGP(10)=XISIMIN
      ELSE
       IMAGP(10)=XISIMIN
     .         +(XISIMAX-XISIMIN)*DFLOAT(IXISI-1)/DFLOAT(NXISI-1)
      ENDIF

      DO IMUPI=1,NMUPI
      IF(NMUPI.EQ.1)THEN
       IMAGP(11)=MUPIMIN
      ELSE
       IMAGP(11)=MUPIMIN
     .         +(MUPIMAX-MUPIMIN)*DFLOAT(IMUPI-1)/DFLOAT(NMUPI-1)
      ENDIF

      DO IMSPI=1,NMSPI
      IF(NMSPI.EQ.1)THEN
       IMAGP(12)=MSPIMIN
      ELSE
       IMAGP(12)=MSPIMIN
     .         +(MSPIMAX-MSPIMIN)*DFLOAT(IMSPI-1)/DFLOAT(NMSPI-1)
      ENDIF

      DO IM3HI=1,NM3HI
      IF(NM3HI.EQ.1)THEN
       IMAGP(13)=M3HIMIN
      ELSE
       IMAGP(13)=M3HIMIN
     .         +(M3HIMAX-M3HIMIN)*DFLOAT(IM3HI-1)/DFLOAT(NM3HI-1)
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

      ITOT=ITOT+1

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

 11   CALL OUTPUT(PAR,PROB,IFAIL)
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
       XIFN=MIN(RXIF,XIFN)
       XIFNN=MAX(RXIF,XIFNN)
       XISN=MIN(RXIS,XISN)
       XISNN=MAX(RXIS,XISNN)
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
       XISIN=MIN(IXISQ2,XISIN)
       XISINN=MAX(IXISQ2,XISINN)
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

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

*   Summary of the scanning:
*   Number of points that passed/failed the tests
*   and range for scanned parameters

      CALL ERROR(TOT,NTOT,NFAIL)

      END


      SUBROUTINE INPUT(PAR,NPAR)

*******************************************************************
*   This subroutine reads SM and NMSSM parameters from input file   .
*******************************************************************

      IMPLICIT NONE

      CHARACTER CHINL*120,CHBLCK*60,CHDUM*120

      INTEGER I,NLINE,INL,ICH,IX,IVAL,Q2FIX,VFLAG,Z3FLAG
      INTEGER N0,NLOOP,NBER,NPAR,ERR,NTOT,GMUFLAG,HFLAG
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,MSFLAG
      INTEGER AKFLAG,ALFLAG,XISIFLAG,OMGFLAG,MAFLAG,MOFLAG
      INTEGER NL,NK,NTB,NMU,NMUP,NMSP,NM3H
      INTEGER N1,N2,N3,N4,NM1,NM2,NM3,NAQ,NMQ
      INTEGER NAL,NMA,NXIF,NAK,NMP,NXIS
      INTEGER NLI,NKI,NM1I,NM2I,NM3I,NAU3I,NAD3I
      INTEGER NAE3I,NXIFI,NXISI,NMUPI,NMSPI,NM3HI
      INTEGER CFLAG(5)

      DOUBLE PRECISION PAR(*),VAL
      DOUBLE PRECISION ACC,XITLA,XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,Q2,Q2MIN
      DOUBLE PRECISION LMIN,LMAX,KMIN,KMAX,TBMIN,TBMAX,MUMIN,MUMAX
      DOUBLE PRECISION ALMIN,ALMAX,AKMIN,AKMAX,XIFMIN,XIFMAX
      DOUBLE PRECISION XISMIN,XISMAX,MUPMIN,MUPMAX,MSPMIN,MSPMAX
      DOUBLE PRECISION M3HMIN,M3HMAX,MAMIN,MAMAX,MPMIN,MPMAX
      DOUBLE PRECISION M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX
      DOUBLE PRECISION AQMIN,AQMAX,MQMIN,MQMAX
      DOUBLE PRECISION M1IMIN,M1IMAX,M2IMIN,M2IMAX,M3IMIN,M3IMAX
      DOUBLE PRECISION AU3IMIN,AU3IMAX,AD3IMIN,AD3IMAX,AE3IMIN
      DOUBLE PRECISION AE3IMAX,LIMIN,LIMAX,KIMIN,KIMAX,XIFIMIN
      DOUBLE PRECISION XIFIMAX,XISIMIN,XISIMAX,MUPIMIN,MUPIMAX
      DOUBLE PRECISION MSPIMIN,MSPIMAX,M3HIMIN,M3HIMAX
      DOUBLE PRECISION REALP(14),IMAGP(14)

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/MINMAX/LMIN,LMAX,KMIN,KMAX,TBMIN,TBMAX,MUMIN,MUMAX,
     . ALMIN,ALMAX,AKMIN,AKMAX,XIFMIN,XIFMAX,XISMIN,XISMAX,
     . MUPMIN,MUPMAX,MSPMIN,MSPMAX,M3HMIN,M3HMAX,MAMIN,MAMAX,
     . MPMIN,MPMAX,M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX,
     . AQMIN,AQMAX,MQMIN,MQMAX
      COMMON/MINIMAX/M1IMIN,M1IMAX,M2IMIN,M2IMAX,M3IMIN,M3IMAX,
     . AU3IMIN,AU3IMAX,AD3IMIN,AD3IMAX,AE3IMIN,AE3IMAX,LIMIN,
     . LIMAX,KIMIN,KIMAX,XIFIMIN,XIFIMAX,XISIMIN,XISIMAX,MUPIMIN,
     . MUPIMAX,MSPIMIN,MSPIMAX,M3HIMIN,M3HIMAX
      COMMON/STEPS/NTOT,NL,NK,NTB,NMU,NMUP,NMSP,NM3H,
     . N1,N2,N3,N4,NM1,NM2,NM3,NAQ,NMQ,
     . NLI,NKI,NM1I,NM2I,NM3I,NAU3I,NAD3I,
     . NAE3I,NXIFI,NXISI,NMUPI,NMSPI,NM3HI
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . MSFLAG,AKFLAG,ALFLAG,XISIFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/VFLAG/VFLAG
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
      TBMIN=1d99
      TBMAX=1d99
      M1MIN=1d99
      M1MAX=1d99
      M2MIN=1d99
      M2MAX=1d99
      M3MIN=1d99
      M3MAX=1d99
      LMIN=1d99
      LMAX=1d99
      KMIN=0d0
      KMAX=1d99
      ALMIN=1d99
      ALMAX=1d99
      AKMIN=1d99
      AKMAX=1d99
      MUMIN=1d99
      MUMAX=1d99
      XIFMIN=1d99
      XIFMAX=1d99
      XISMIN=1d99
      XISMAX=1d99
      MUPMIN=0d0
      MUPMAX=1d99
      MSPMIN=0d0
      MSPMAX=1d99
      M3HMIN=0d0
      M3HMAX=1d99
      MAMIN=1d99
      MAMAX=1d99
      MPMIN=1d99
      MPMAX=1d99
      AQMIN=1d99
      AQMAX=1d99
      MQMIN=1d99
      MQMAX=1d99
      M1IMIN=0d0
      M1IMAX=1d99
      M2IMIN=0d0
      M2IMAX=1d99
      M3IMIN=0d0
      M3IMAX=1d99
      AU3IMIN=0d0
      AU3IMAX=1d99
      AD3IMIN=0d0
      AD3IMAX=1d99
      AE3IMIN=0d0
      AE3IMAX=1d99
      LIMIN=0d0
      LIMAX=1d99
      KIMIN=0d0
      KIMAX=1d99
      XIFIMIN=0d0
      XIFIMAX=1d99
      XISIMIN=1d99
      XISIMAX=1d99
      MUPIMIN=0d0
      MUPIMAX=1d99
      MSPIMIN=0d0
      MSPIMAX=1d99
      M3HIMIN=0d0
      M3HIMAX=1d99
      NM1=0
      NM2=1
      NM3=0
      NL=1
      NK=1
      NTB=1
      NMU=1
      NAL=0
      NAK=0
      NXIF=0
      NXIS=0
      NMA=0
      NMP=0
      NMUP=1
      NMSP=1
      NM3H=1
      NAQ=1
      NMQ=1
      NLI=1
      NKI=1
      NM1I=1
      NM2I=1
      NM3I=1
      NAU3I=1
      NAD3I=1
      NAE3I=1
      NXIFI=1
      NXISI=0
      NMUPI=1
      NMSPI=1
      NM3HI=1

*   DEFAULT VALUES FOR FLAGS
      GMUFLAG=1
      HFLAG=0
      OMGFLAG=0
      MOFLAG=7
      VFLAG=0
      M1FLAG=0
      M3FLAG=0
      XISIFLAG=1
      DO I=1,5
       CFLAG(I)=1
      ENDDO

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
       IF(IX.EQ.37) TBMIN=VAL
       IF(IX.EQ.38) TBMAX=VAL

*   READ EXTPAR
      ELSEIF(CHBLCK(1:6).EQ.'EXTPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0) Q2=VAL**2
       IF(IX.EQ.107) M1MIN=VAL
       IF(IX.EQ.108) M1MAX=VAL
       IF(IX.EQ.207) M2MIN=VAL
       IF(IX.EQ.208) M2MAX=VAL
       IF(IX.EQ.307) M3MIN=VAL
       IF(IX.EQ.308) M3MAX=VAL
       IF(IX.EQ.117) AQMIN=VAL
       IF(IX.EQ.118) AQMAX=VAL
       IF(IX.EQ.13) REALP(8)=VAL
       IF(IX.EQ.16) PAR(25)=VAL
       IF(IX.EQ.32) PAR(18)=VAL**2
       IF(IX.EQ.33) PAR(10)=VAL**2
       IF(IX.EQ.35) PAR(19)=VAL**2
       IF(IX.EQ.36) PAR(11)=VAL**2
       IF(IX.EQ.437) MQMIN=VAL
       IF(IX.EQ.438) MQMAX=VAL
       IF(IX.EQ.617) LMIN=VAL
       IF(IX.EQ.618) LMAX=VAL
       IF(IX.EQ.627) KMIN=VAL
       IF(IX.EQ.628) KMAX=VAL
       IF(IX.EQ.637) ALMIN=VAL
       IF(IX.EQ.638) ALMAX=VAL
       IF(IX.EQ.647) AKMIN=VAL
       IF(IX.EQ.648) AKMAX=VAL
       IF(IX.EQ.657) MUMIN=VAL
       IF(IX.EQ.658) MUMAX=VAL
       IF(IX.EQ.667) XIFMIN=VAL
       IF(IX.EQ.668) XIFMAX=VAL
       IF(IX.EQ.677) XISMIN=VAL
       IF(IX.EQ.678) XISMAX=VAL
       IF(IX.EQ.687) MUPMIN=VAL
       IF(IX.EQ.688) MUPMAX=VAL
       IF(IX.EQ.697) MSPMIN=VAL
       IF(IX.EQ.698) MSPMAX=VAL
       IF(IX.EQ.727) M3HMIN=VAL
       IF(IX.EQ.728) M3HMAX=VAL
       IF(IX.EQ.1247) MAMIN=VAL
       IF(IX.EQ.1248) MAMAX=VAL
       IF(IX.EQ.1257) MPMIN=VAL
       IF(IX.EQ.1258) MPMAX=VAL

*   READ IMEXTPAR
      ELSEIF(CHBLCK(1:8).EQ.'IMEXTPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.107) M1IMIN=VAL
       IF(IX.EQ.108) M1IMAX=VAL
       IF(IX.EQ.207) M2IMIN=VAL
       IF(IX.EQ.208) M2IMAX=VAL
       IF(IX.EQ.307) M3IMIN=VAL
       IF(IX.EQ.308) M3IMAX=VAL
       IF(IX.EQ.117) AU3IMIN=VAL
       IF(IX.EQ.118) AU3IMAX=VAL
       IF(IX.EQ.127) AD3IMIN=VAL
       IF(IX.EQ.128) AD3IMAX=VAL
       IF(IX.EQ.137) AE3IMIN=VAL
       IF(IX.EQ.138) AE3IMAX=VAL
       IF(IX.EQ.617) LIMIN=VAL
       IF(IX.EQ.618) LIMAX=VAL
       IF(IX.EQ.627) KIMIN=VAL
       IF(IX.EQ.628) KIMAX=VAL
       IF(IX.EQ.667) XIFIMIN=VAL
       IF(IX.EQ.668) XIFIMAX=VAL
       IF(IX.EQ.677) XISIMIN=VAL
       IF(IX.EQ.678) XISIMAX=VAL
       IF(IX.EQ.687) MUPIMIN=VAL
       IF(IX.EQ.688) MUPIMAX=VAL
       IF(IX.EQ.697) MSPIMIN=VAL
       IF(IX.EQ.698) MSPIMAX=VAL
       IF(IX.EQ.727) M3HIMIN=VAL
       IF(IX.EQ.728) M3HIMAX=VAL
       IF(IX.EQ.637 .OR. IX.EQ.638)THEN
        WRITE(0,1)"IM(ALAMBDA) IS NOT AN INPUT"
        ERR=1
       ENDIF
       IF(IX.EQ.647 .OR. IX.EQ.648)THEN
        WRITE(0,1)"IM(AKAPPA) IS NOT AN INPUT"
        ERR=1
       ENDIF
       IF(IX.EQ.657 .OR. IX.EQ.658)THEN
        WRITE(0,1)"IM(MUEFF) IS NOT AN INPUT"
        ERR=1
       ENDIF

*   READ STEPS
       ELSEIF(CHBLCK(1:6).EQ.'STEPS')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.109) NM1=IVAL
       IF(IX.EQ.209) NM2=IVAL
       IF(IX.EQ.309) NM3=IVAL
       IF(IX.EQ.39) NTB=IVAL
       IF(IX.EQ.619) NL=IVAL
       IF(IX.EQ.629) NK=IVAL
       IF(IX.EQ.639) NAL=IVAL
       IF(IX.EQ.649) NAK=IVAL
       IF(IX.EQ.659) NMU=IVAL
       IF(IX.EQ.669) NXIF=IVAL
       IF(IX.EQ.679) NXIS=IVAL
       IF(IX.EQ.689) NMUP=IVAL
       IF(IX.EQ.699) NMSP=IVAL
       IF(IX.EQ.729) NM3H=IVAL
       IF(IX.EQ.1249) NMA=IVAL
       IF(IX.EQ.1259) NMP=IVAL
       IF(IX.EQ.119) NAQ=IVAL
       IF(IX.EQ.439) NMQ=IVAL
       IF(IX.EQ.105) NM1I=IVAL
       IF(IX.EQ.205) NM2I=IVAL
       IF(IX.EQ.305) NM3I=IVAL
       IF(IX.EQ.115) NAU3I=IVAL
       IF(IX.EQ.125) NAD3I=IVAL
       IF(IX.EQ.135) NAE3I=IVAL
       IF(IX.EQ.615) NLI=IVAL
       IF(IX.EQ.625) NKI=IVAL
       IF(IX.EQ.665) NXIFI=IVAL
       IF(IX.EQ.675) NXISI=IVAL
       IF(IX.EQ.685) NMUPI=IVAL
       IF(IX.EQ.695) NMSPI=IVAL
       IF(IX.EQ.725) NM3HI=IVAL

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
      IF(LMIN.EQ.1d99)THEN
       WRITE(0,1)"RE(LAMBDA)MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ELSEIF(LMIN.LE.0d0 .OR. LMAX.LE.0d0)THEN
       WRITE(0,1)"LAMBDA MUST BE STRICTLY POSITIVE"
       ERR=1
      ENDIF
      IF(TBMIN.EQ.1d99)THEN
       WRITE(0,1)"TANBMIN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ELSEIF(TBMIN.LE.0d0 .OR. TBMAX.LE.0d0)THEN
       WRITE(0,1)"TANB MUST BE STRICTLY POSITIVE"
       ERR=1
      ENDIF
      IF(MUMIN.EQ.1d99)THEN
       WRITE(0,1)"RE(MUEFF)MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ELSEIF(MUMIN.EQ.0d0 .OR. MUMAX.EQ.0d0)THEN
       WRITE(0,1)"MUEFF MUST BE NON ZERO"
       ERR=1
      ENDIF
      IF(M1MIN.EQ.1d99 .AND. M1MAX.NE.1d99)THEN
       WRITE(0,1)"RE(M1)MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M2MIN.EQ.1d99)THEN
       WRITE(0,1)"RE(M2)MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M3MIN.EQ.1d99 .AND. M3MAX.NE.1d99)THEN
       WRITE(0,1)"RE(M3)MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(AQMIN.EQ.1d99)THEN
       WRITE(0,1)"RE(AU3) MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(REALP(8).EQ.1d99)THEN
       WRITE(0,1)"RE(AE3) MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MQMIN.EQ.1d99)THEN
       WRITE(0,1)"MQ3 MUST BE GIVEN IN BLOCK EXTPAR"
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
      IF(MAMIN.EQ.1d99 .AND. MAMAX.NE.1d99)THEN
       WRITE(0,1)"MAMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MPMIN.EQ.1d99 .AND. MPMAX.NE.1d99)THEN
       WRITE(0,1)"MPMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(ALMIN.EQ.1d99 .AND. ALMAX.NE.1d99)THEN
       WRITE(0,1)"RE(AL)MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(AKMIN.EQ.1d99 .AND. AKMAX.NE.1d99)THEN
       WRITE(0,1)"RE(AK)MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XIFMIN.EQ.1d99 .AND. XIFMAX.NE.1d99)THEN
       WRITE(0,1)"RE(XIF)MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XISMIN.EQ.1d99 .AND. XISMAX.NE.1d99)THEN
       WRITE(0,1)"RE(XIS)MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XISIMIN.EQ.1d99 .AND. XISIMAX.NE.1d99)THEN
       WRITE(0,1)"IM(XIS)MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

      IF(PAR(18).EQ.1d99)PAR(18)=PAR(10)
      IF(PAR(19).EQ.1d99)PAR(19)=PAR(11)
      IF(PAR(25).EQ.1d99)PAR(25)=REALP(8)

*   Relations between (RE(ALAMBDA), MA, RE(XIF)) and (RE(AKAPPA), MP, RE(XIS))

      IF(ALMIN.NE.1d99 .AND. XIFMIN.NE.1d99 .AND. MAMIN.NE.1d99)THEN
       WRITE(0,1)"AT MOST 2 OF THE 3 PARAMETERS ALAMBDA, MA AND XIF",
     . " CAN BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

      IF(AKMIN.NE.1d99 .AND. XISMIN.NE.1d99 .AND. MPMIN.NE.1d99)THEN
       WRITE(0,1)"AT MOST 2 OF THE 3 PARAMETERS AKAPPA, MP AND XIS",
     . " CAN BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

      IF((KIMIN.EQ.0d0 .AND. (KIMAX.EQ.0d0 .OR. KIMAX.EQ.1d99))
     ..AND.(KMIN.EQ.0d0 .AND. (KMAX.EQ.0d0 .OR. KMAX.EQ.1d99)))THEN
       IF((AKMIN.NE.0d0 .AND. AKMIN.NE.1d99) .OR. (AKMIN.EQ.0d0
     ..AND. AKMAX.NE.0d0 .AND. AKMAX.NE.1d99))THEN
        WRITE(0,1)"IF KAPPA IS 0, AKAPPA MUST BE 0"
        ERR=1
       ELSE
        AKMIN=0d0
        AKMAX=0d0
        IF(XISMIN.NE.1d99 .AND. MPMIN.NE.1d99)THEN
         WRITE(0,1)"IF KAPPA IS 0, EITHER MP OR XIS",
     .   " CAN BE GIVEN IN BLOCK EXTPAR"
         ERR=1
        ENDIF
       ENDIF
       IF(XISIMIN.NE.1d99)THEN
        WRITE(0,1)"IF KAPPA IS 0, IM(XIS) IS NOT AN INPUT"
        ERR=1
       ENDIF
       IF(NXISI.NE.0)THEN
        WRITE(0,1)"IF KAPPA IS 0, NXISI CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       XISIFLAG=0
      ENDIF

*   Set default values

      IF(ALMIN.EQ.1d99.AND.MAMIN.EQ.1d99.AND.XIFMIN.EQ.1d99)THEN
       ALMIN=0d0
       XIFMIN=0d0
      ELSEIF(ALMIN.EQ.1d99.AND.MAMIN.EQ.1d99)THEN
       ALMIN=0d0
      ELSEIF(ALMIN.EQ.1d99.AND.XIFMIN.EQ.1d99)THEN
       XIFMIN=0d0
      ELSEIF(MAMIN.EQ.1d99.AND.XIFMIN.EQ.1d99)THEN
       XIFMIN=0d0
      ENDIF

      IF(AKMIN.EQ.1d99.AND.MPMIN.EQ.1d99.AND.XISMIN.EQ.1d99)THEN
       AKMIN=0d0
       XISMIN=0d0
      ELSEIF(AKMIN.EQ.1d99.AND.MPMIN.EQ.1d99)THEN
       AKMIN=0d0
      ELSEIF(AKMIN.EQ.1d99.AND.XISMIN.EQ.1d99)THEN
       XISMIN=0d0
      ELSEIF(MPMIN.EQ.1d99.AND.XISMIN.EQ.1d99)THEN
       XISMIN=0d0
      ENDIF

*   Set MAFLAG, SCANFLAGS

      IF(MAMIN.EQ.1d99)MAFLAG=0
      IF(ALMIN.EQ.1d99)MAFLAG=1
      IF(XIFMIN.EQ.1d99)MAFLAG=2
      IF(AKMIN.EQ.1d99)MAFLAG=MAFLAG+3
      IF(XISMIN.EQ.1d99)MAFLAG=MAFLAG+6
      IF(M1MIN.NE.1d99)M1FLAG=1
      IF(M3MIN.NE.1d99)M3FLAG=1

      IF(XIFMIN.EQ.1d99)XIFMIN=0d0
      IF(XISMIN.EQ.1d99)XISMIN=0d0
      IF(XISIMIN.EQ.1d99)XISIMIN=0d0
      IF(NXISI.EQ.0)NXISI=1

*   Bounds

      IF(TBMAX.EQ.1d99)TBMAX=TBMIN
      IF(M1MAX.EQ.1d99)M1MAX=M1MIN
      IF(M2MAX.EQ.1d99)M2MAX=M2MIN
      IF(M3MAX.EQ.1d99)M3MAX=M3MIN
      IF(LMAX.EQ.1d99)LMAX=LMIN
      IF(KMAX.EQ.1d99)KMAX=KMIN
      IF(ALMAX.EQ.1d99)ALMAX=ALMIN
      IF(AKMAX.EQ.1d99)AKMAX=AKMIN
      IF(MUMAX.EQ.1d99)MUMAX=MUMIN
      IF(XIFMAX.EQ.1d99)XIFMAX=XIFMIN
      IF(XISMAX.EQ.1d99)XISMAX=XISMIN
      IF(MUPMAX.EQ.1d99)MUPMAX=MUPMIN
      IF(MSPMAX.EQ.1d99)MSPMAX=MSPMIN
      IF(M3HMAX.EQ.1d99)M3HMAX=M3HMIN
      IF(MAMAX.EQ.1d99)MAMAX=MAMIN
      IF(MPMAX.EQ.1d99)MPMAX=MPMIN
      IF(AQMAX.EQ.1d99)AQMAX=AQMIN
      IF(MQMAX.EQ.1d99)MQMAX=MQMIN
      IF(M1IMAX.EQ.1d99)M1IMAX=M1IMIN
      IF(M2IMAX.EQ.1d99)M2IMAX=M2IMIN
      IF(M3IMAX.EQ.1d99)M3IMAX=M3IMIN
      IF(AU3IMAX.EQ.1d99)AU3IMAX=AU3IMIN
      IF(AD3IMAX.EQ.1d99)AD3IMAX=AD3IMIN
      IF(AE3IMAX.EQ.1d99)AE3IMAX=AE3IMIN
      IF(LIMAX.EQ.1d99)LIMAX=LIMIN
      IF(KIMAX.EQ.1d99)KIMAX=KIMIN
      IF(XIFIMAX.EQ.1d99)XIFIMAX=XIFIMIN
      IF(XISIMAX.EQ.1d99)XISIMAX=XISIMIN
      IF(MUPIMAX.EQ.1d99)MUPIMAX=MUPIMIN
      IF(MSPIMAX.EQ.1d99)MSPIMAX=MSPIMIN
      IF(M3HIMAX.EQ.1d99)M3HIMAX=M3HIMIN

*   Number of points

      IF(MOD(MAFLAG,3).EQ.0)THEN
       IF(NMA.NE.0)THEN
        WRITE(0,1)"NMA CANNOT BE GIVEN IN BLOCK STEPS"
       ERR=1
       ENDIF
       IF(NAL.EQ.0)NAL=1
       IF(NXIF.EQ.0)NXIF=1
       N1=NAL
       N2=NXIF
      ENDIF

      IF(MOD(MAFLAG,3).EQ.1)THEN
       IF(NAL.NE.0)THEN
        WRITE(0,1)"NAL CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NMA.EQ.0)NMA=1
       IF(NXIF.EQ.0)NXIF=1
       N1=NMA
       N2=NXIF
      ENDIF

      IF(MOD(MAFLAG,3).EQ.2)THEN
       IF(NXIF.NE.0)THEN
        WRITE(0,1)"NXIF CANNOT BE GIVEN IN BLOCK STEPS"
       ERR=1
       ENDIF
       IF(NAL.EQ.0)NAL=1
       IF(NMA.EQ.0)NMA=1
       N1=NAL
       N2=NMA
      ENDIF

      IF(MAFLAG/3.EQ.0)THEN
       IF(NMP.NE.0)THEN
        WRITE(0,1)"NMP CANNOT BE GIVEN IN BLOCK STEPS"
       ERR=1
       ENDIF
       IF(NAK.EQ.0)NAK=1
       IF(NXIS.EQ.0)NXIS=1
       N3=NAK
       N4=NXIS
      ENDIF

      IF(MAFLAG/3.EQ.1)THEN
       IF(NAK.NE.0)THEN
        WRITE(0,1)"NAK CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NMP.EQ.0)NMP=1
       IF(NXIS.EQ.0)NXIS=1
       N3=NMP
       N4=NXIS
      ENDIF

      IF(MAFLAG/3.EQ.2)THEN
       IF(NXIS.NE.0)THEN
        WRITE(0,1)"NXIS CANNOT BE GIVEN IN BLOCK STEPS"
       ERR=1
       ENDIF
       IF(NAK.EQ.0)NAK=1
       IF(NMP.EQ.0)NMP=1
       N3=NAK
       N4=NMP
      ENDIF

      IF(M1FLAG.EQ.0)THEN
       IF(NM1.NE.0)THEN
        WRITE(0,1)"NM1 CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ELSE
        NM1=1
       ENDIF
      ELSE
       IF(NM1.EQ.0)NM1=1
      ENDIF

      IF(M3FLAG.EQ.0)THEN
       IF(NM3.NE.0)THEN
        WRITE(0,1)"NM3 CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ELSE
        NM3=1
       ENDIF
      ELSE
       IF(NM3.EQ.0)NM3=1
      ENDIF

      IF(NM1.LE.0 .OR. NM2.LE.0 .OR. NM3.LE.0 .OR. NTB.LE.0 .OR.
     . NL.LE.0 .OR. NK.LE.0 .OR. N1.LE.0 .OR. N2.LE.0 .OR.
     . N3.LE.0 .OR. N4.LE.0 .OR. NMU.LE.0 .OR. NMUP.LE.0 .OR.
     . NMSP.LE.0 .OR. NM3H.LE.0 .OR. NAQ.LE.0 .OR. NMQ.LE.0 .OR.
     . NLI.LE.0 .OR. NKI.LE.0 .OR. NM1I.LE.0 .OR. NM2I.LE.0 .OR.
     . NM3I.LE.0 .OR. NAU3I.LE.0 .OR. NAD3I.LE.0 .OR.
     . NAE3I.LE.0 .OR. NXIFI.LE.0 .OR. NXISI.LE.0 .OR.
     . NMUPI.LE.0 .OR. NMSPI.LE.0 .OR. NM3HI.LE.0)THEN
       WRITE(0,1)"WRONG NUMBER OF POINTS IN BLOCK STEPS"
       ERR=1
      ENDIF

*   Check for Z3 breaking terms

      IF(MOD(MAFLAG,3).EQ.2 .OR. MAFLAG/3.EQ.2
     ..OR. MUPMIN.NE.0d0 .OR. MUPMAX.NE.0d0
     ..OR. MSPMIN.NE.0d0 .OR. MSPMAX.NE.0d0
     ..OR. XIFMIN.NE.0d0 .OR. XIFMAX.NE.0d0
     ..OR. XISMIN.NE.0d0 .OR. XISMAX.NE.0d0
     ..OR. M3HMIN.NE.0d0 .OR. M3HMAX.NE.0d0
     ..OR. MUPIMIN.NE.0d0 .OR. MUPIMAX.NE.0d0
     ..OR. MSPIMIN.NE.0d0 .OR. MSPIMAX.NE.0d0
     ..OR. XIFIMIN.NE.0d0 .OR. XIFIMAX.NE.0d0
     ..OR. XISIMIN.NE.0d0 .OR. XISIMAX.NE.0d0
     ..OR. M3HIMIN.NE.0d0 .OR. M3HIMAX.NE.0d0
     ..OR. (KIMIN.EQ.0d0 .AND. KIMAX.EQ.0d0
     ..AND. KMIN.EQ.0d0 .AND. KMAX.EQ.0d0))THEN
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

*   Warnings

      IF(TBMIN.EQ.TBMAX .AND. NTB.GT.1)THEN
       WRITE(0,1)"WARNING TANBMIN = TANBMAX => NTB=1"
       NTB=1
      ENDIF
      IF(M1MIN.EQ.M1MAX .AND. NM1.GT.1)THEN
       WRITE(0,1)"WARNING M1MIN = M1MAX => NM1=1"
       NM1=1
      ENDIF
      IF(M2MIN.EQ.M2MAX .AND. NM2.GT.1)THEN
       WRITE(0,1)"WARNING M2MIN = M2MAX => NM2=1"
       NM2=1
      ENDIF
      IF(M3MIN.EQ.M3MAX .AND. NM3.GT.1)THEN
       WRITE(0,1)"WARNING M3MIN = M3MAX => NM3=1"
       NM3=1
      ENDIF
      IF(LMIN.EQ.LMAX .AND. NL.GT.1)THEN
       WRITE(0,1)"WARNING LMIN = LMAX => NL=1"
       NL=1
      ENDIF
      IF(KMIN.EQ.KMAX .AND. NK.GT.1)THEN
       WRITE(0,1)"WARNING KMIN = KMAX => NK=1"
       NK=1
      ENDIF
      IF(MUMIN.EQ.MUMAX .AND. NMU.GT.1)THEN
       WRITE(0,1)"WARNING MUEFFMIN = MUEFFMAX => NMU=1"
       NMU=1
      ENDIF
      IF(MUPMIN.EQ.MUPMAX .AND. NMUP.GT.1)THEN
       WRITE(0,1)"WARNING MUPMIN = MUPMAX => NMUP=1"
       NMUP=1
      ENDIF
      IF(MSPMIN.EQ.MSPMAX .AND. NMSP.GT.1)THEN
       WRITE(0,1)"WARNING MSPMIN = MSPMAX => NMSP=1"
       NMSP=1
      ENDIF
      IF(M3HMIN.EQ.M3HMAX .AND. NM3H.GT.1)THEN
       WRITE(0,1)"WARNING M3HMIN = M3HMAX => N3H=1"
       NM3H=1
      ENDIF
      IF(AQMIN.EQ.AQMAX .AND. NAQ.GT.1)THEN
       WRITE(0,1)"WARNING AU3MIN = AU3MAX => NAU3=1"
       NAQ=1
      ENDIF
      IF(MQMIN.EQ.MQMAX .AND. NMQ.GT.1)THEN
       WRITE(0,1)"WARNING MQ3MIN = MQ3MAX => NMQ3=1"
       NMQ=1
      ENDIF
      IF(LIMIN.EQ.LIMAX .AND. NLI.GT.1)THEN
       WRITE(0,1)"WARNING IM(LAMBDA)MIN = IM(LAMBDA)MAX => NLI=1"
       NLI=1
      ENDIF
      IF(KIMIN.EQ.KIMAX .AND. NKI.GT.1)THEN
       WRITE(0,1)"WARNING IM(KAPPA)MIN = IM(KAPPA)MAX => NKI=1"
       NKI=1
      ENDIF
      IF(M1IMIN.EQ.M1IMAX .AND. NM1I.GT.1)THEN
       WRITE(0,1)"WARNING IM(M1)MIN = IM(M1)MAX => NM1I=1"
       NM1I=1
      ENDIF
      IF(M2IMIN.EQ.M2IMAX .AND. NM2I.GT.1)THEN
       WRITE(0,1)"WARNING IM(M2)MIN = IM(M2)MAX => NM2I=1"
       NM2I=1
      ENDIF
      IF(M3IMIN.EQ.M3IMAX .AND. NM3I.GT.1)THEN
       WRITE(0,1)"WARNING IM(M3)MIN = IM(M3)MAX => NM3I=1"
       NM3I=1
      ENDIF
      IF(AU3IMIN.EQ.AU3IMAX .AND. NAU3I.GT.1)THEN
       WRITE(0,1)"WARNING IM(AU3)MIN = IM(AU3)MAX => NAU3I=1"
       NAU3I=1
      ENDIF
      IF(AD3IMIN.EQ.AD3IMAX .AND. NAD3I.GT.1)THEN
       WRITE(0,1)"WARNING IM(AD3)MIN = IM(AD3)MAX => NAD3I=1"
       NAD3I=1
      ENDIF
      IF(AE3IMIN.EQ.AE3IMAX .AND. NAE3I.GT.1)THEN
       WRITE(0,1)"WARNING IM(AE3)MIN = IM(AE3)MAX => NAE3I=1"
       NAE3I=1
      ENDIF
      IF(XIFIMIN.EQ.XIFIMAX .AND. NXIFI.GT.1)THEN
       WRITE(0,1)"WARNING IM(XIF)MIN = IM(XIF)MAX => NXIFI=1"
       NXIFI=1
      ENDIF
      IF(XISIMIN.EQ.XISIMAX .AND. NXISI.GT.1)THEN
       WRITE(0,1)"WARNING IM(XIS)MIN = IM(XIS)MAX => NXISI=1"
       NXISI=1
      ENDIF
      IF(MUPIMIN.EQ.MUPIMAX .AND. NMUPI.GT.1)THEN
       WRITE(0,1)"WARNING IM(MUP)MIN = IM(MUP)MAX => NMUPI=1"
       NMUPI=1
      ENDIF
      IF(MSPIMIN.EQ.MSPIMAX .AND. NMSPI.GT.1)THEN
       WRITE(0,1)"WARNING IM(MSP)MIN = IM(MSP)MAX => NMSPI=1"
       NMSPI=1
      ENDIF
      IF(M3HIMIN.EQ.M3HIMAX .AND. NM3HI.GT.1)THEN
       WRITE(0,1)"WARNING IM(M3H)MIN = IM(M3H)MAX => NM3HI=1"
       NM3HI=1
      ENDIF

      IF(TBMIN.NE.TBMAX .AND. NTB.EQ.1)THEN
       WRITE(0,10)"WARNING NTB=1 => TANBMAX=TANBMIN=",TBMIN
       TBMAX=TBMIN
      ENDIF
      IF(M1MIN.NE.M1MAX .AND. NM1.EQ.1)THEN
       WRITE(0,10)"WARNING NM1=1 => M1MAX=M1MIN=",M1MIN
       M1MAX=M1MIN
      ENDIF
      IF(M2MIN.NE.M2MAX .AND. NM2.EQ.1)THEN
       WRITE(0,10)"WARNING NM2=1 => M2MAX=M2MIN=",M2MIN
       M2MAX=M2MIN
      ENDIF
      IF(M3MIN.NE.M3MAX .AND. NM3.EQ.1)THEN
       WRITE(0,10)"WARNING NM3=1 => M3MAX=M3MIN=",M3MIN
       M3MAX=M3MIN
      ENDIF
      IF(LMIN.NE.LMAX .AND. NL.EQ.1)THEN
       WRITE(0,10)"WARNING NL=1 => LMAX=LMIN=",LMIN
       LMAX=LMIN
      ENDIF
      IF(KMIN.NE.KMAX .AND. NK.EQ.1)THEN
       WRITE(0,10)"WARNING NK=1 => KMAX=KMIN=",KMIN
       KMAX=KMIN
      ENDIF
      IF(MUMIN.NE.MUMAX .AND. NMU.EQ.1)THEN
       WRITE(0,10)"WARNING NMU=1 => MUEFFMAX=MUEFFMIN=",MUMIN
       MUMAX=MUMIN
      ENDIF
      IF(MUPMIN.NE.MUPMAX .AND. NMUP.EQ.1)THEN
       WRITE(0,10)"WARNING NMUP=1 => MUPMAX=MUPMIN=",MUPMIN
       MUPMAX=MUPMIN
      ENDIF
      IF(MSPMIN.NE.MSPMAX .AND. NMSP.EQ.1)THEN
       WRITE(0,10)"WARNING NMSP=1 =>  MSPMAX= MSPMIN=",MSPMIN
       MSPMAX=MSPMIN
      ENDIF
      IF(M3HMIN.NE.M3HMAX .AND. NM3H.EQ.1)THEN
       WRITE(0,10)"WARNING NM3H=1 => M3HMAX=M3HMIN=",M3HMIN
       M3HMAX=M3HMIN
      ENDIF
      IF(AQMIN.NE.AQMAX .AND. NAQ.EQ.1)THEN
       WRITE(0,10)"WARNING NAU3=1 => AU3MAX=AU3MIN=",AQMIN
       AQMAX=AQMIN
      ENDIF
      IF(MQMIN.NE.MQMAX .AND. NMQ.EQ.1)THEN
       WRITE(0,10)"WARNING NMQ3=1 => MQ3MAX=MQ3MIN=",MQMIN
       MQMAX=MQMIN
      ENDIF
      IF(LIMIN.NE.LIMAX .AND. NLI.EQ.1)THEN
       WRITE(0,10)"WARNING NLI=1 => IM(LAMBDA)MAX=IM(LAMBDA)MIN=",
     . LIMIN
       LIMAX=LIMIN
      ENDIF
      IF(KIMIN.NE.KIMAX .AND. NKI.EQ.1)THEN
       WRITE(0,10)"WARNING NKI=1 => IM(KAPPA)MAX=IM(KAPPA)MIN=",KIMIN
       KIMAX=KIMIN
      ENDIF
      IF(M1IMIN.NE.M1IMAX .AND. NM1I.EQ.1)THEN
       WRITE(0,10)"WARNING NM1I=1 => IM(M1)MAX=IM(M1)MIN=",M1IMIN
       M1IMAX=M1IMIN
      ENDIF
      IF(M2IMIN.NE.M2IMAX .AND. NM2I.EQ.1)THEN
       WRITE(0,10)"WARNING NM2I=1 => IM(M2)MAX=IM(M2)MIN=",M2IMIN
       M2IMAX=M2IMIN
      ENDIF
      IF(M3IMIN.NE.M3IMAX .AND. NM3I.EQ.1)THEN
       WRITE(0,10)"WARNING NM3I=1 => IM(M3)MAX=IM(M3)MIN=",M3IMIN
       M3IMAX=M3IMIN
      ENDIF
      IF(AU3IMIN.NE.AU3IMAX .AND. NAU3I.EQ.1)THEN
       WRITE(0,10)"WARNING NAU3I=1 => IM(AU3)MAX=IM(AU3)MIN=",AU3IMIN
       AU3IMAX=AU3IMIN
      ENDIF
      IF(AD3IMIN.NE.AD3IMAX .AND. NAD3I.EQ.1)THEN
       WRITE(0,10)"WARNING NAD3I=1 => IM(AD3)MAX=IM(AD3)MIN=",AD3IMIN
       AD3IMAX=AD3IMIN
      ENDIF
      IF(AE3IMIN.NE.AE3IMAX .AND. NAE3I.EQ.1)THEN
       WRITE(0,10)"WARNING NAE3I=1 => IM(AE3)MAX=IM(AE3)MIN=",AE3IMIN
       AE3IMAX=AE3IMIN
      ENDIF
      IF(XIFIMIN.NE.XIFIMAX .AND. NXIFI.EQ.1)THEN
       WRITE(0,10)"WARNING NXIFI=1 => IM(XIF)MAX=IM(XIF)MIN=",XIFIMIN
       XIFIMAX=XIFIMIN
      ENDIF
      IF(XISIMIN.NE.XISIMAX .AND. NXISI.EQ.1)THEN
       WRITE(0,10)"WARNING NXISI=1 => IM(XIS)MAX=IM(XIS)MIN=",XISIMIN
       XISIMAX=XISIMIN
      ENDIF
      IF(MUPIMIN.NE.MUPIMAX .AND. NMUPI.EQ.1)THEN
       WRITE(0,10)"WARNING NMUPI=1 => IM(MUP)MAX=IM(MUP)MIN=",MUPIMIN
       MUPIMAX=MUPIMIN
      ENDIF
      IF(MSPIMIN.NE.MSPIMAX .AND. NMSPI.EQ.1)THEN
       WRITE(0,10)"WARNING NMSPI=1 => IM(MSP)MAX=IM(MSP)MIN=",MSPIMIN
       MSPIMAX=MSPIMIN
      ENDIF
      IF(M3HIMIN.NE.M3HIMAX .AND. NM3HI.EQ.1)THEN
       WRITE(0,10)"WARNING NM3HI=1 => IM(M3H)MAX=IM(M3H)MIN=",M3HIMIN
       M3HIMAX=M3HIMIN
      ENDIF

      IF(MOD(MAFLAG,3).EQ.0)THEN
       IF(ALMIN.EQ.ALMAX .AND. NAL.GT.1)THEN
        WRITE(0,1)"WARNING ALMIN=ALMAX => NAL=1"
        NAL=1
        N1=NAL
       ENDIF
       IF(XIFMIN.EQ.XIFMAX .AND. NXIF.GT.1)THEN
        WRITE(0,1)"WARNING XIFMIN=XIFMAX => NXIF=1"
        NXIF=1
        N2=NXIF
       ENDIF
       IF(ALMIN.NE.ALMAX .AND. NAL.EQ.1)THEN
        WRITE(0,10)"WARNING NAL=1 => ALMAX=ALMIN=",ALMIN
        ALMAX=ALMIN
       ENDIF
       IF(XIFMIN.NE.XIFMAX .AND. NXIF.EQ.1)THEN
        WRITE(0,10)"WARNING NXIF=1 => XIFMAX=XIFMIN=",XIFMIN
        XIFMAX=XIFMIN
       ENDIF
      ENDIF

      IF(MOD(MAFLAG,3).EQ.1)THEN
       IF(MAMIN.EQ.MAMAX .AND. NMA.GT.1)THEN
        WRITE(0,1)"WARNING MAMIN=MAMAX => NMA=1"
        NMA=1
        N1=NMA
       ENDIF
       IF(XIFMIN.EQ.XIFMAX .AND. NXIF.GT.1)THEN
        WRITE(0,1)"WARNING XIFMIN=XIFMAX => NXIF=1"
        NXIF=1
        N2=NXIF
       ENDIF
       IF(MAMIN.NE.MAMAX .AND. NMA.EQ.1)THEN
        WRITE(0,10)"WARNING NMA=1 => MAMAX=MAMIN=",MAMIN
        MAMAX=MAMIN
       ENDIF
       IF(XIFMIN.NE.XIFMAX .AND. NXIF.EQ.1)THEN
        WRITE(0,10)"WARNING NXIF=1 => XIFMAX=XIFMIN=",XIFMIN
        XIFMAX=XIFMIN
       ENDIF
      ENDIF

      IF(MOD(MAFLAG,3).EQ.2)THEN
       IF(ALMIN.EQ.ALMAX .AND. NAL.GT.1)THEN
        WRITE(0,1)"WARNING ALMIN=ALMAX => NAL=1"
        NAL=1
        N1=NAL
       ENDIF
       IF(MAMIN.EQ.MAMAX .AND. NMA.GT.1)THEN
        WRITE(0,1)"WARNING MAMIN=MAMAX => NMA=1"
        NMA=1
        N2=NMA
       ENDIF
       IF(ALMIN.NE.ALMAX .AND. NAL.EQ.1)THEN
        WRITE(0,10)"WARNING NAL=1 => ALMAX=ALMIN=",ALMIN
        ALMAX=ALMIN
       ENDIF
       IF(MAMIN.NE.MAMAX .AND. NMA.EQ.1)THEN
        WRITE(0,10)"WARNING NMA=1 => MAMAX=MAMIN=",MAMIN
        MAMAX=MAMIN
       ENDIF
      ENDIF

      IF(MAFLAG/3.EQ.0)THEN
       IF(AKMIN.EQ.AKMAX .AND. NAK.GT.1)THEN
        WRITE(0,1)"WARNING AKMIN=AKMAX => NAK=1"
        NAK=1
        N3=NAK
       ENDIF
       IF(XISMIN.EQ.XISMAX .AND. NXIS.GT.1)THEN
        WRITE(0,1)"WARNING XISMIN=XISMAX => NXIS=1"
        NXIS=1
        N4=NXIS
       ENDIF
       IF(AKMIN.NE.AKMAX .AND. NAK.EQ.1)THEN
        WRITE(0,10)"WARNING NAK=1 => AKMAX=AKMIN=",AKMIN
        AKMAX=AKMIN
       ENDIF
       IF(XISMIN.NE.XISMAX .AND. NXIS.EQ.1)THEN
        WRITE(0,10)"WARNING NXIS=1 => XISMAX=XISMIN=",XISMIN
        XISMAX=XISMIN
       ENDIF
      ENDIF

      IF(MAFLAG/3.EQ.1)THEN
       IF(MPMIN.EQ.MPMAX .AND. NMP.GT.1)THEN
        WRITE(0,1)"WARNING MPMIN=MPMAX => NMP=1"
        NMP=1
        N3=NMP
       ENDIF
       IF(XISMIN.EQ.XISMAX .AND. NXIS.GT.1)THEN
        WRITE(0,1)"WARNING XISMIN=XISMAX => NXIS=1"
        NXIS=1
        N4=NXIS
       ENDIF
       IF(MPMIN.NE.MPMAX .AND. NMP.EQ.1)THEN
        WRITE(0,10)"WARNING NMP=1 => MPMAX=MPMIN=",MPMIN
        MPMAX=MPMIN
       ENDIF
       IF(XISMIN.NE.XISMAX .AND. NXIS.EQ.1)THEN
        WRITE(0,10)"WARNING NXIS=1 => XISMAX=XISMIN=",XISMIN
        XISMAX=XISMIN
       ENDIF
      ENDIF

      IF(MAFLAG/3.EQ.2)THEN
       IF(AKMIN.EQ.AKMAX .AND. NAK.GT.1)THEN
        WRITE(0,1)"WARNING AKMIN=AKMAX => NAK=1"
        NAK=1
        N3=NAK
       ENDIF
       IF(MPMIN.EQ.MPMAX .AND. NMP.GT.1)THEN
        WRITE(0,1)"WARNING MPMIN=MPMAX => NMP=1"
        NMP=1
        N4=NMP
       ENDIF
       IF(AKMIN.NE.AKMAX .AND. NAK.EQ.1)THEN
        WRITE(0,10)"WARNING NAK=1 => AKMAX=AKMIN=",AKMIN
        AKMAX=AKMIN
       ENDIF
       IF(MPMIN.NE.MPMAX .AND. NMP.EQ.1)THEN
        WRITE(0,10)"WARNING NMP=1 => MPMAX=MPMIN=",MPMIN
        MPMAX=MPMIN
       ENDIF
      ENDIF

*   Total number of points

      NTOT=NTB*NM1*NM2*NM3*NL*NMU*NK*NMUP*NMSP*NM3H*N1*N2*N3*N4*NAQ*NMQ
     .    *NLI*NKI*NM1I*NM2I*NM3I*NAU3I*NAD3I*NAE3I
     .    *NXIFI*NXISI*NMUPI*NMSPI*NM3HI

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
 10   FORMAT(A,E10.3)

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
      RES(8)=RXIF
      RES(9)=RXIS
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
