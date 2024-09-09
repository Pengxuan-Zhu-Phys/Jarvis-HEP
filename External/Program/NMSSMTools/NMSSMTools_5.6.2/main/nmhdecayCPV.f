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

      INTEGER NPROB,NPAR
      PARAMETER (NPROB=82,NPAR=25)
      INTEGER IFAIL,I
      INTEGER OMGFLAG,MAFLAG,MOFLAG

      DOUBLE PRECISION PAR(NPAR),PROB(NPROB),PI
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H,DELMB,DELML,DEL1
      DOUBLE PRECISION REALP(14),IMAGP(14),MYPHASES(16)

      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/REAL_IMAG/REALP,IMAGP
      COMMON/MYPHASES/MYPHASES
      COMMON/DELMB/DELMB,DELML,DEL1

      PI=4d0*DATAN(1d0)

*   I/O files

      OPEN(15,FILE='inp',STATUS='UNKNOWN')
      OPEN(17,FILE='spectr', STATUS= 'UNKNOWN')
      OPEN(18,FILE='decay', STATUS= 'UNKNOWN')
      OPEN(19,FILE='omega',STATUS='UNKNOWN')

*   Initialization

      CALL INITIALIZE()

*   Reading of the input parameters

      CALL INPUT(PAR,NPAR)

*   Moduli and phases

      DO I=1,16
       MYPHASES(I)=0d0
      ENDDO

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

      CALL constsusypart_CPV(PROB)
      CALL LEP_Higgs_CPV(PROB)
      CALL tevatron_chiggs_CPV(PROB)

      CALL LHC_HIGGS_CPV(PROB)

      CALL bottomonium_CPV(PROB)
      CALL bsg_CPV(PAR,PROB)

      CALL magnmu_CPV(PROB)

      CALL checkmin_CPV(PAR,PROB)

      CALL EDM_CPV(PAR,PROB)

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
 
      END


      SUBROUTINE INPUT(PAR,NPAR)

*******************************************************************
*   This subroutine reads SM and NMSSM parameters from input file   .
*******************************************************************

      IMPLICIT NONE

      CHARACTER CHINL*120,CHBLCK*60,CHDUM*120

      INTEGER I,NLINE,INL,ICH,IX,IVAL,Q2FIX,OMGFLAG,MAFLAG,MOFLAG
      INTEGER ERR,N0,NLOOP,NBER,NPAR,VFLAG,OUTFLAG,Z3FLAG,CFLAG(5)

      DOUBLE PRECISION PAR(*),VAL,PI
      DOUBLE PRECISION ACC,XITLA,XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,Q2,Q2MIN
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION REALP(14),IMAGP(14)

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/VFLAG/VFLAG
      COMMON/OUTFLAG/OUTFLAG
      COMMON/CFLAG/CFLAG
      COMMON/REAL_IMAG/REALP,IMAGP

      PI=4d0*DATAN(1d0)

*   INITIALIZATION OF THE SUSY PARAMETERS
      DO I=1,NPAR
       PAR(I)=1d99
      ENDDO

*   DEFAULT VALUE FOR REALP,IMAGP
      DO I=1,14
       REALP(I)=1d99
       IMAGP(I)=0d0
      ENDDO
      REALP(2)=0d0
      REALP(11)=0d0
      REALP(12)=0d0
      REALP(13)=0d0
      IMAGP(10)=1d99

*   DEFAULT VALUES FOR FLAGS
      OMGFLAG=0
      MOFLAG=7
      VFLAG=0
      OUTFLAG=0
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

      READ(15,'(A120)',END=29,ERR=999) CHINL

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
       IF(IX.EQ.14) VFLAG=IVAL
       IF(IX.EQ.16) OUTFLAG=IVAL

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
       IF(IX.EQ.3) PAR(3)=VAL

*   READ EXTPAR
      ELSEIF(CHBLCK(1:6).EQ.'EXTPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0) Q2=VAL**2
       IF(IX.EQ.1) REALP(3)=VAL
       IF(IX.EQ.2) REALP(4)=VAL
       IF(IX.EQ.3) REALP(5)=VAL
       IF(IX.EQ.11) REALP(6)=VAL
       IF(IX.EQ.12) REALP(7)=VAL
       IF(IX.EQ.13) REALP(8)=VAL
       IF(IX.EQ.16) PAR(25)=VAL
       IF(IX.EQ.32) PAR(18)=VAL**2
       IF(IX.EQ.33) PAR(10)=VAL**2
       IF(IX.EQ.35) PAR(19)=VAL**2
       IF(IX.EQ.36) PAR(11)=VAL**2
       IF(IX.EQ.42) PAR(15)=VAL**2
       IF(IX.EQ.43) PAR(7)=VAL**2
       IF(IX.EQ.45) PAR(16)=VAL**2
       IF(IX.EQ.46) PAR(8)=VAL**2
       IF(IX.EQ.48) PAR(17)=VAL**2
       IF(IX.EQ.49) PAR(9)=VAL**2
       IF(IX.EQ.61) REALP(1)=VAL
       IF(IX.EQ.62) REALP(2)=VAL
       IF(IX.EQ.63) PAR(5)=VAL
       IF(IX.EQ.64) PAR(6)=VAL
       IF(IX.EQ.65) REALP(14)=VAL
       IF(IX.EQ.66) REALP(9)=VAL
       IF(IX.EQ.67) REALP(10)=VAL
       IF(IX.EQ.68) REALP(11)=VAL
       IF(IX.EQ.69) REALP(12)=VAL
       IF(IX.EQ.72) REALP(13)=VAL
       IF(IX.EQ.124) PAR(23)=VAL
       IF(IX.EQ.125) PAR(24)=VAL

*   READ IMEXTPAR
      ELSEIF(CHBLCK(1:8).EQ.'IMEXTPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.1) IMAGP(3)=VAL
       IF(IX.EQ.2) IMAGP(4)=VAL
       IF(IX.EQ.3) IMAGP(5)=VAL
       IF(IX.EQ.11) IMAGP(6)=VAL
       IF(IX.EQ.12) IMAGP(7)=VAL
       IF(IX.EQ.13) IMAGP(8)=VAL
       IF(IX.EQ.61) IMAGP(1)=VAL
       IF(IX.EQ.62) IMAGP(2)=VAL
       IF(IX.EQ.66) IMAGP(9)=VAL
       IF(IX.EQ.67) IMAGP(10)=VAL
       IF(IX.EQ.68) IMAGP(11)=VAL
       IF(IX.EQ.69) IMAGP(12)=VAL
       IF(IX.EQ.72) IMAGP(13)=VAL
       IF(IX.EQ.63)THEN
        WRITE(0,1)"IM(ALAMBDA) IS NOT AN INPUT"
        ERR=1
       ENDIF
       IF(IX.EQ.64)THEN
        WRITE(0,1)"IM(AKAPPA) IS NOT AN INPUT"
        ERR=1
       ENDIF
       IF(IX.EQ.65)THEN
        WRITE(0,1)"IM(MUEFF) IS NOT AN INPUT"
        ERR=1
       ENDIF

      ENDIF

      GOTO 21

*   END OF READING FROM INPUT FILE

*   Check for errors

 29   IF(VFLAG.LT.0 .OR. VFLAG.GT.1)THEN
       WRITE(0,1)"VFLAG MUST BE IN [0-1]"
       ERR=1
      ENDIF
      IF(OUTFLAG.LT.0 .OR. OUTFLAG.GT.1)THEN
       WRITE(0,1)"OUTFLAG MUST BE IN [0-1]"
       ERR=1
      ENDIF
      IF(REALP(1).EQ.1d99)THEN
       WRITE(0,1)"RE(LAMBDA) MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ELSEIF(REALP(1).LE.0d0)THEN
       WRITE(0,1)"RE(LAMBDA) MUST BE STRICTLY POSITIVE"
       ERR=1
      ENDIF
      IF(PAR(3).EQ.1d99)THEN
       WRITE(0,1)"TANB MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ELSEIF(PAR(3).LE.0d0)THEN
       WRITE(0,1)"TANB MUST BE STRICTLY POSITIVE"
       ERR=1
      ENDIF
      IF(REALP(14).EQ.1d99)THEN
       WRITE(0,1)"RE(MUEFF) MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ELSEIF(REALP(14).EQ.0d0)THEN
       WRITE(0,1)"RE(MUEFF) MUST BE NON ZERO"
       ERR=1
      ENDIF
      IF(REALP(4).EQ.1d99)THEN
       WRITE(0,1)"RE(M2) MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ELSE
       IF(REALP(3).EQ.1d99)REALP(3)=REALP(4)/2d0
       IF(REALP(5).EQ.1d99)REALP(5)=REALP(4)*3d0
      ENDIF
      IF(REALP(6).EQ.1d99)THEN
       WRITE(0,1)"RE(AU3) MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(REALP(7).EQ.1d99)THEN
       WRITE(0,1)"RE(AD3) MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(REALP(8).EQ.1d99)THEN
       WRITE(0,1)"RE(AE3) MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(7).EQ.1d99)THEN
       WRITE(0,1)"MQ3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(8).EQ.1d99)THEN
       WRITE(0,1)"MU3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(9).EQ.1d99)THEN
       WRITE(0,1)"MD3 MUST BE GIVEN IN BLOCK EXTPAR"
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
      DO I=15,19
       IF(PAR(I).EQ.1d99)PAR(I)=PAR(I-8)
      ENDDO
      IF(PAR(25).EQ.1d99)PAR(25)=REALP(8)

*   Relations between (RE(ALAMBDA), MA, RE(XIF)) and (RE(AKAPPA), MP, RE(XIS))

      IF(PAR(5).NE.1d99 .AND. PAR(23).NE.1d99 .AND. REALP(9).NE.1d99)
     . THEN
       WRITE(0,1)"AT MOST 2 OF THE 3 PARAMETERS MA, RE(ALAMBDA)"
     . ," AND RE(XIF) CAN BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

      IF(PAR(6).NE.1d99 .AND. PAR(24).NE.1d99 .AND. REALP(10).NE.1d99)
     . THEN
       WRITE(0,1)"AT MOST 2 OF THE 3 PARAMETERS MP, RE(AKAPPA)"
     . ," AND RE(XIS) CAN BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

      IF(REALP(2).EQ.0d0 .AND. IMAGP(2).EQ.0d0)THEN
       IF(PAR(6).NE.0d0 .AND. PAR(6).NE.1d99)THEN
        WRITE(0,1)"IF KAPPA IS 0, AKAPPA MUST BE 0"
        ERR=1
       ELSE
        PAR(6)=0d0
       ENDIF
       IF(PAR(24).NE.1d99 .AND. REALP(10).NE.1d99)THEN
        WRITE(0,1)"EITHER MP OR RE(XIS) CAN BE GIVEN IN BLOCK EXTPAR"
       ERR=1
       ENDIF
       IF(IMAGP(10).NE.1d99)THEN
        WRITE(0,1)"IF KAPPA IS 0, IM(XIS) IS NOT AN INPUT"
        ERR=1
       ENDIF
      ENDIF

*   Set default values

      IF(PAR(5).EQ.1d99.AND.PAR(23).EQ.1d99.AND.REALP(9).EQ.1d99)THEN
       PAR(5)=0d0
       REALP(9)=0d0
      ELSEIF(PAR(5).EQ.1d99.AND.PAR(23).EQ.1d99)THEN
       PAR(5)=0d0
      ELSEIF(PAR(5).EQ.1d99.AND.REALP(9).EQ.1d99)THEN
       REALP(9)=0d0
      ELSEIF(PAR(23).EQ.1d99.AND.REALP(9).EQ.1d99)THEN
       REALP(9)=0d0
      ENDIF

      IF(PAR(6).EQ.1d99.AND.PAR(24).EQ.1d99.AND.REALP(10).EQ.1d99)THEN
       PAR(6)=0d0
       REALP(10)=0d0
      ELSEIF(PAR(6).EQ.1d99.AND.PAR(24).EQ.1d99)THEN
       PAR(6)=0d0
      ELSEIF(PAR(6).EQ.1d99.AND.REALP(10).EQ.1d99)THEN
       REALP(10)=0d0
      ELSEIF(PAR(24).EQ.1d99.AND.REALP(10).EQ.1d99)THEN
       REALP(10)=0d0
      ENDIF

*   Set MAFLAG

      IF(PAR(23).EQ.1d99)MAFLAG=0
      IF(PAR(5).EQ.1d99)MAFLAG=1
      IF(REALP(9).EQ.1d99)MAFLAG=2
      IF(PAR(6).EQ.1d99)MAFLAG=MAFLAG+3
      IF(REALP(10).EQ.1d99)MAFLAG=MAFLAG+6

      IF(REALP(9).EQ.1d99)REALP(9)=0d0
      IF(REALP(10).EQ.1d99)REALP(10)=0d0
      IF(IMAGP(10).EQ.1d99)IMAGP(10)=0d0

*   Check for Z3 breaking terms

      IF(MOD(MAFLAG,3).EQ.2 .OR. MAFLAG/3.EQ.2
     ..OR. REALP(9).NE.0d0 .OR. IMAGP(9).NE.0d0
     ..OR. REALP(10).NE.0d0 .OR. IMAGP(10).NE.0d0
     ..OR. REALP(11).NE.0d0 .OR. IMAGP(11).NE.0d0
     ..OR. REALP(12).NE.0d0 .OR. IMAGP(12).NE.0d0
     ..OR. REALP(13).NE.0d0 .OR. IMAGP(13).NE.0d0
     ..OR.(REALP(2).EQ.0d0 .AND. IMAGP(2).EQ.0d0))THEN
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

      INTEGER IFAIL,OMGFLAG,MAFLAG,MOFLAG,VFLAG,OUTFLAG,Q2FIX

      DOUBLE PRECISION PAR(*),PROB(*),PI
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
      DOUBLE PRECISION VUS,VCB,VUB,TANB,SINB,COSB,Q2,Q2MIN
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
     .    AMU_CPV,PHASES(16),REALP(14),IMAGP(14)
      DOUBLE PRECISION RXIS,RXIF,IAlQ2,IAkQ2,IXISQ2,MHuQ2,MHdQ2,MSQ2

      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/VFLAG/VFLAG
      COMMON/OUTFLAG/OUTFLAG
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
      COMMON/MYPHASES/PHASES
      COMMON/REAL_IMAG/REALP,IMAGP
      COMMON/REXI/RXIS,RXIF,IAlQ2,IAkQ2,IXISQ2,MHuQ2,MHdQ2,MSQ2

      PI=4d0*DATAN(1d0)
      TANB=PAR(3)
      COSB=1d0/DSQRT(1d0+TANB**2)
      SINB=TANB*COSB

      WRITE(17,899) "# NMSSMTools OUTPUT IN SLHA FORMAT"
      WRITE(17,899) "# Info about spectrum calculator"
      WRITE(17,899) "BLOCK SPINFO   # Program information"
      WRITE(17,900) 1,"NMSSMTools # Spectrum calculator"
      WRITE(17,900) 2,"5.6.2      # Version number"

      IF(DABS(MGL_CPV).LT.1d3 .AND. MGL_CPV.NE.0d0)
     . WRITE(17,900) 3,"# Gluino mass < 1 TeV"
      IF(MIN(dsqrt(MSU2P(1)),dsqrt(MSU2P(2)),
     .       dsqrt(MSD2P(1)),dsqrt(MSD2P(2))).LT.1d3)
     . WRITE(17,900) 3,"# First generation squark masses < 1 TeV"
      IF(PROB(1).NE.0d0)
     . WRITE(17,900) 3,"# Chargino excluded by LEP"
      IF(PROB(2).NE.0d0)
     . WRITE(17,900) 3,"# Neutralinos excluded by LEP"
      IF(PROB(3).NE.0d0)
     . WRITE(17,900) 3,"# Charged Higgs excluded by LEP"
      IF(PROB(4).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, ind. of h decay"
      IF(PROB(5).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> bb"
      IF(PROB(6).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> tautau"
      IF(PROB(7).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> invisible"
      IF(PROB(8).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> 2jets"
      IF(PROB(9).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> 2photons"
      IF(PROB(10).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> h -> 4bs"
      IF(PROB(11).NE.0d0 .OR. PROB(41).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> hh -> 4taus"
      IF(PROB(12).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> hh -> 2bs 2taus"
      IF(PROB(13).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> Z -> hh (Z width)"
      IF(PROB(14).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hh -> 4bs"
      IF(PROB(15).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hh -> 4taus"
      IF(PROB(16).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hh -> 2bs 2taus"
      IF(PROB(17).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hh -> hhh -> 6bs"
      IF(PROB(18).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hh -> hhh -> 6taus"
      IF(PROB(19).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> hh, h -> light
     . pair"
      IF(PROB(20).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by stop -> b l sneutrino"
      IF(PROB(21).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by stop -> neutralino c"
      IF(PROB(22).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by sbottom -> neutralino b"
      IF(PROB(23).NE.0d0)
     . WRITE(17,900) 3,"# Squark/gluino too light"
      IF(PROB(24).NE.0d0)
     . WRITE(17,900) 3,"# Selectron/smuon too light"
      IF(PROB(25).NE.0d0)
     . WRITE(17,900) 3,"# Stau too light"
      IF(PROB(28).NE.0d0)
     . WRITE(17,900) 3,"# Unphysical global minimum"
      IF(PROB(29).NE.0d0)
     . WRITE(17,900) 3,"# Higgs soft masses >> Msusy"
      IF(PROB(32).NE.0d0)
     . WRITE(17,900) 3,"# b -> s gamma more than 2 sigma away"
      IF(PROB(33).NE.0d0)
     . WRITE(17,900) 3,"# Delta M_s more than 2 sigma away"
      IF(PROB(34).NE.0d0)
     . WRITE(17,900) 3,"# Delta M_d more than 2 sigma away"
      IF(PROB(35).NE.0d0)
     . WRITE(17,900) 3,"# B_s -> mu+ mu- more than 2 sigma away"
      IF(PROB(36).NE.0d0)
     . WRITE(17,900) 3,"# B+ -> tau nu_tau more than 2 sigma away"
      IF(PROB(37).NE.0d0)
     . WRITE(17,900) 3,"# Muon magn. mom. more than 2 sigma away"
      IF(PROB(38).NE.0d0)
     . WRITE(17,900) 3,"# excluded by Upsilon(1S) -> H/A gamma"
      IF(PROB(38).LT.0d0)
     . WRITE(17,900) 3,"# (but A width> 10 MeV)"
      IF(PROB(39).NE.0d0)
     . WRITE(17,900) 3,
     . "# Excluded etab(1S) mass difference (BABAR - theory)"
      IF(PROB(42).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by top -> b H+, H+ -> c s"
      IF(PROB(43).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by top -> b H+, H+ -> tau nu_tau"
      IF(PROB(44).NE.0d0)
     . WRITE(17,900) 3,
     . "# Excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus"
      IF(PROB(45).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by t -> bH+ (ATLAS)"
      IF(PROB(46).NE.0d0)
     . WRITE(17,918) 3,"# No Higgs in the",MHMIN,MHMAX," GeV mass range"
      IF(PROB(51).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ggF/bb->H/A->tautau (ATLAS+CMS)"
      IF(PROB(52).NE.0d0)
     . WRITE(17,900) 3,
     ."# Excluded by H_SM->AA->4leptons/2lept.+2b (ATLAS+CMS)"
      IF(PROB(53).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ggF->H/A->gamgam (ATLAS)"
      IF(PROB(54).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by EDMs"
      IF(PROB(55).NE.0d0)
     . WRITE(17,900) 3,"# b -> d gamma more than 2 sigma away"
      IF(PROB(56).NE.0d0)
     . WRITE(17,900) 3,"# B_d -> mu+ mu- more than 2 sigma away"
      IF(PROB(58).NE.0d0)
     . WRITE(17,900) 3,"# b -> c tau nu more than 2 sigma away (as SM)"
      IF(PROB(63).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by H_SM->AA->4gam (ATLAS+CMS)"
      IF(PROB(67).NE.0d0)
     . WRITE(17,900) 3,"# k_W/Z(H_SM) more than 2 sigma away"
      IF(PROB(68).NE.0d0)
     . WRITE(17,900) 3,"# k_top(H_SM) more than 2 sigma away"
      IF(PROB(69).NE.0d0)
     . WRITE(17,900) 3,"# k_bot(H_SM) more than 2 sigma away (as SM)"
      IF(PROB(70).NE.0d0)
     . WRITE(17,900) 3,"# k_glu(H_SM) more than 2 sigma away"
      IF(PROB(71).NE.0d0)
     . WRITE(17,900) 3,"# k_gam(H_SM) more than 2 sigma away"
      IF(PROB(72).NE.0d0)
     . WRITE(17,900) 3,"# k_tau(H_SM) more than 2 sigma away (as SM)"
      IF(PROB(73).NE.0d0)
     . WRITE(17,900) 3,"# B_bsm(H_SM) more than 2 sigma away"
      IF(PROB(74).NE.0d0)
     . WRITE(17,900) 3,"# excluded by H_SM->Z+A->had (ATLAS)"
      IF(PROB(75).NE.0d0)
     . WRITE(17,900) 3,"# excluded by H/A->toptop (CMS)"
      IF(PROB(76).NE.0d0)
     . WRITE(17,900) 3,"# excluded by A->Z+H_SM->bb (CMS)"
      IF(PROB(77).NE.0d0)
     . WRITE(17,900) 3,"# excluded by A->Z+H_SM->bb (ATLAS)"
      IF(PROB(78).NE.0d0)
     . WRITE(17,900) 3,"# excluded by H/A->Z+A/H (CMS)"
      IF(PROB(79).NE.0d0)
     . WRITE(17,900) 3,"# excluded by H/A->Z+A/H (ATLAS)"
      IF(PROB(80).NE.0d0)
     . WRITE(17,900) 3,"# excluded by H/A->HSM+H/A->2b2tau (CMS)"
      IF(PROB(81).NE.0d0)
     . WRITE(17,900) 3,"# excluded by H/A->HSM+H/A->4b (CMS)"
      IF(PROB(82).NE.0d0)
     . WRITE(17,900) 3,"# excluded by H->HSM+HSM (ATLAS)"

      IF(IFAIL.EQ.1.OR.IFAIL.EQ.5)
     . WRITE(17,900) 4,"# M_H1^2<1"
      IF(IFAIL.EQ.4.OR.IFAIL.EQ.5)
     . WRITE(17,900) 4,"# M_HC^2<1"
      IF(IFAIL.EQ.8)
     . WRITE(17,900) 4,"# Negative sfermion mass squared"

      WRITE(17,899) "# Input parameters"
      WRITE(17,899) "BLOCK MODSEL"
      WRITE(17,921) 1,0,"IMOD"
      WRITE(17,921) 3,1,"NMSSM particle content"
      WRITE(17,921) 5,2,"0: without CP-Violation, 2: with CP-Violation"
      WRITE(17,921) 10,0,"ISCAN"
      WRITE(17,921) 14,VFLAG,"H-> VV,VV*,(V*V*)"
      WRITE(17,921) 16,OUTFLAG,"Extra BLOCK's yes/no"

      WRITE(17,899) "BLOCK SMINPUTS"
      WRITE(17,901) 1,1d0/ALEMMZ,"ALPHA_EM^-1(MZ)"
      WRITE(17,901) 2,GF,"GF"
      WRITE(17,901) 3,ALSMZ,"ALPHA_S(MZ)"
      WRITE(17,901) 4,MZ,"MZ"
      WRITE(17,901) 5,MB,"MB(MB)"
      WRITE(17,901) 6,MT,"MTOP (POLE MASS)"
      WRITE(17,901) 7,MTAU,"MTAU"
      WRITE(17,899) "# SMINPUTS Beyond SLHA:"
      WRITE(17,906) "MW:",MW
      WRITE(17,906) "MS:",MS
      WRITE(17,906) "MC:",MC
      WRITE(17,906) "VUS:",VUS
      WRITE(17,906) "VCB:",VCB
      WRITE(17,906) "VUB:",VUB

      WRITE(17,899) "BLOCK MINPAR"
      IF(Q2FIX.EQ.1)WRITE(17,901) 0,DSQRT(Q2),"REN. SCALE"
      WRITE(17,901) 3,TANB,"TANBETA(MZ)"

      WRITE(17,899) "BLOCK EXTPAR"
      WRITE(17,901) 1,REALP(3),"RE(M1)"
      WRITE(17,901) 2,REALP(4),"RE(M2)"
      WRITE(17,901) 3,REALP(5),"RE(M3)"
      WRITE(17,901) 11,REALP(6),"RE(ATOP)"
      WRITE(17,901) 12,REALP(7),"RE(ABOTTOM)"
      WRITE(17,901) 13,REALP(8),"RE(ATAU)"
      WRITE(17,901) 16,PAR(25),"RE(AMUON)"

      WRITE(17,901) 31,DSQRT(PAR(18)),"LEFT SELECTRON"
      WRITE(17,901) 32,DSQRT(PAR(18)),"LEFT SMUON"
      WRITE(17,901) 33,DSQRT(PAR(10)),"LEFT STAU"

      WRITE(17,901) 34,DSQRT(PAR(19)),"RIGHT SELECTRON"
      WRITE(17,901) 35,DSQRT(PAR(19)),"RIGHT SMUON"
      WRITE(17,901) 36,DSQRT(PAR(11)),"RIGHT STAU"

      WRITE(17,901) 41,DSQRT(PAR(15)),"LEFT 1ST GEN. SQUARKS"
      WRITE(17,901) 42,DSQRT(PAR(15)),"LEFT 2ND GEN. SQUARKS"
      WRITE(17,901) 43,DSQRT(PAR(7)),"LEFT 3RD GEN. SQUARKS"

      WRITE(17,901) 44,DSQRT(PAR(16)),"RIGHT U-SQUARKS"
      WRITE(17,901) 45,DSQRT(PAR(16)),"RIGHT C-SQUARKS"
      WRITE(17,901) 46,DSQRT(PAR(8)),"RIGHT T-SQUARKS"

      WRITE(17,901) 47,DSQRT(PAR(17)),"RIGHT D-SQUARKS"
      WRITE(17,901) 48,DSQRT(PAR(17)),"RIGHT S-SQUARKS"
      WRITE(17,901) 49,DSQRT(PAR(9)),"RIGHT B-SQUARKS"

      WRITE(17,901) 61,REALP(1),"RE(LAMBDA)"
      IF(PAR(2).NE.0d0)THEN
       WRITE(17,901) 62,REALP(2),"RE(KAPPA)"
      ENDIF
      IF(MOD(MAFLAG,3).NE.1)THEN
       WRITE(17,901) 63,PAR(5),"RE(ALAMBDA)"
      ELSE
       WRITE(17,920) 63,PAR(5),"RE(ALAMBDA)"
      ENDIF
      IF(PAR(2).NE.0d0)THEN
       IF(MAFLAG/3.NE.1)THEN
        WRITE(17,901) 64,PAR(6),"RE(AKAPPA)"
       ELSE
        WRITE(17,920) 64,PAR(6),"RE(AKAPPA)"
       ENDIF
      ENDIF
      WRITE(17,901) 65,REALP(14),"RE(MUEFF)"
      IF(MOD(MAFLAG,3).NE.2)THEN
       IF(REALP(9).NE.0d0)
     .  WRITE(17,901) 66,RXIF,"RE(XIF)"
      ELSE
       WRITE(17,920) 66,RXIF,"RE(XIF)"
      ENDIF
      IF(MAFLAG/3.NE.2)THEN
       IF(REALP(10).NE.0d0)
     .  WRITE(17,901) 67,RXIS,"RE(XIS)"
      ELSE
       WRITE(17,920) 67,RXIS,"RE(XIS)"
      ENDIF
      IF(REALP(11).NE.0d0)
     . WRITE(17,901) 68,REALP(11),"RE(MUP)"
      IF(REALP(12).NE.0d0)
     . WRITE(17,901) 69,REALP(12),"RE(MSP)"
      IF(REALP(13).NE.0d0)
     .  WRITE(17,901) 72,REALP(13),"RE(M3H)"
      IF(MOD(MAFLAG,3).NE.0)THEN
       WRITE(17,901) 124,PAR(23),"MA AT QSTSB"
      ELSE
       WRITE(17,920) 124,PAR(23),"MA AT QSTSB"
      ENDIF
      IF(MAFLAG/3.NE.0)THEN
       WRITE(17,901) 125,PAR(24),"MP AT QSTSB"
      ELSE
       WRITE(17,920) 125,PAR(24),"MP AT QSTSB"
      ENDIF

      WRITE(17,899) "BLOCK IMEXTPAR"
      WRITE(17,901) 1,IMAGP(3),"IM(M1)"
      WRITE(17,901) 2,IMAGP(4),"IM(M2)"
      WRITE(17,901) 3,IMAGP(5),"IM(M3)"
      WRITE(17,901) 11,IMAGP(6),"IM(ATOP)"
      WRITE(17,901) 12,IMAGP(7),"IM(ABOTTOM)"
      WRITE(17,901) 13,IMAGP(8),"IM(ATAU)"
      WRITE(17,901) 61,IMAGP(1),"IM(LAMBDA)"
      WRITE(17,901) 62,IMAGP(2),"IM(KAPPA)"
      WRITE(17,920) 63,IALQ2,"IM(ALAMBDA)"
      WRITE(17,920) 64,IAKQ2,"IM(AKAPPA)"
      WRITE(17,920) 65,IMAGP(14),"IM(MUEFF)"
      IF(IMAGP(9).NE.0d0)
     . WRITE(17,901) 66,IMAGP(9),"IM(XIF)"
      IF(PAR(2).NE.0d0)THEN
       IF(IXISQ2.NE.0d0)
     .  WRITE(17,901) 67,IXISQ2,"IM(XIS)"
      ELSE
       WRITE(17,920) 67,IXISQ2,"IM(XIS)"
      ENDIF
      IF(IMAGP(11).NE.0d0)
     . WRITE(17,901) 68,IMAGP(11),"IM(MUP)"
      IF(IMAGP(12).NE.0d0)
     . WRITE(17,901) 69,IMAGP(12),"IM(MSP)"
      IF(IMAGP(13).NE.0d0)
     .  WRITE(17,901) 72,IMAGP(13),"IM(M3H)"

      WRITE(17,899) "BLOCK ABS.VALUES mod. SIGN(Real part)"
      WRITE(17,901) 1,PAR(20),"+/-|M1|"
      WRITE(17,901) 2,PAR(21),"+/-|M2|"
      WRITE(17,901) 3,PAR(22),"+/-|M3|"
      WRITE(17,901) 11,PAR(12),"+/-|ATOP|"
      WRITE(17,901) 12,PAR(13),"+/-|ABOTTOM|"
      WRITE(17,901) 13,PAR(14),"+/-|ATAU|"
      WRITE(17,901) 61,PAR(1),"+/-|LAMBDA|"
      WRITE(17,901) 62,PAR(2),"+/-|KAPPA|"
      IF(PAR(5).NE.0d0)THEN
       WRITE(17,901) 63,DSQRT(PAR(5)**2+IALQ2**2)*PAR(5)/DABS(PAR(5)),
     . "+/-|ALAMBDA|"
      ELSE
       WRITE(17,901) 63,DABS(IALQ2),"+/-|ALAMBDA|"
      ENDIF
      IF(PAR(6).NE.0d0)THEN
       WRITE(17,901) 64,DSQRT(PAR(6)**2+IAKQ2**2)*PAR(6)/DABS(PAR(6)),
     . "+/-|AKAPPA|"
      ELSE
       WRITE(17,901) 64,DABS(IAKQ2),"+/-|AKAPPA|"
      ENDIF
      WRITE(17,901) 65,PAR(4),"+/-|MUEFF|"
      IF(RXIF.NE.0d0)THEN
       WRITE(17,901) 66,DSQRT(RXIF**2+IMAGP(9)**2)*RXIF/DABS(RXIF),
     . "+/-|XIF|"
      ELSEIF(IMAGP(9).NE.0d0)THEN
       WRITE(17,901) 66,DABS(IMAGP(9)),"+/-|XIF|"
      ENDIF
      IF(RXIS.NE.0d0)THEN
       WRITE(17,901) 67,DSQRT(RXIS**2+IXISQ2**2)*RXIS/DABS(RXIS),
     . "+/-|XIS|"
      ELSEIF(IXISQ2.NE.0d0)THEN
       WRITE(17,901) 67,DABS(IXISQ2),"+/-|XIS|"
      ENDIF
      IF(MUPSUSY.NE.0d0)
     . WRITE(17,901) 68,MUPSUSY,"+/-|MUP|"
      IF(MSPSUSY.NE.0d0)
     . WRITE(17,901) 69,MSPSUSY,"+/-|MSP|"
      IF(M3HSUSY.NE.0d0)
     . WRITE(17,901) 72,M3HSUSY,"+/-|M3H|"

      WRITE(17,899) "BLOCK PHASES (-Pi/2 < PHI < Pi/2)"
      WRITE(17,901) 1,PHASES(3),"PHI(M1)"
      WRITE(17,901) 2,PHASES(4),"PHI(M2)"
      WRITE(17,901) 3,PHASES(5),"PHI(M3)"
      WRITE(17,901) 11,PHASES(6),"PHI(ATOP)"
      WRITE(17,901) 12,PHASES(7),"PHI(ABOTTOM)"
      WRITE(17,901) 13,PHASES(8),"PHI(ATAU)"
      WRITE(17,901) 61,PHASES(1),"PHI(LAMBDA)"
      WRITE(17,901) 62,PHASES(2),"PHI(KAPPA)"
      IF(PAR(5).NE.0d0)THEN
       WRITE(17,901) 63,DATAN(IALQ2/PAR(5)),"PHI(ALAMBDA)"
      ELSEIF(IALQ2.NE.0d0)THEN
       WRITE(17,901) 63,SIGN(PI/2d0,IALQ2),"PHI(ALAMBDA)"
      ELSE
       WRITE(17,901) 63,0d0,"PHI(ALAMBDA)"
      ENDIF
      IF(PAR(6).NE.0d0)THEN
       WRITE(17,901) 64,DATAN(IAKQ2/PAR(6)),"PHI(AKAPPA)"
      ELSEIF(IAKQ2.NE.0d0)THEN
       WRITE(17,901) 64,SIGN(PI/2d0,IAKQ2),"PHI(AKAPPA)"
      ELSE
       WRITE(17,901) 64,0d0,"PHI(AKAPPA)"
      ENDIF
      WRITE(17,901) 65,DATAN(IMAGP(14)/REALP(14)),"PHI(MUEFF)"
      IF(RXIF.NE.0d0)THEN
       WRITE(17,901) 66,DATAN(IMAGP(9)/RXIF),"PHI(XIF)"
      ELSEIF(IMAGP(9).NE.0d0)THEN
       WRITE(17,901) 66,SIGN(PI/2d0,IMAGP(9)),"PHI(XIF)"
      ENDIF
      IF(RXIS.NE.0d0)THEN
       WRITE(17,901) 67,DATAN(IXISQ2/RXIS),"PHI(XIS)"
      ELSEIF(IXISQ2.NE.0d0)THEN
       WRITE(17,901) 67,SIGN(PI/2d0,IXISQ2),"PHI(XIS)"
      ENDIF
      IF(MUPSUSY.NE.0d0)
     . WRITE(17,901) 68,PHASES(15),"PHI(MUP)"
      IF(MSPSUSY.NE.0d0)
     . WRITE(17,901) 69,PHASES(13),"PHI(MSP)"
      IF(M3HSUSY.NE.0d0)
     . WRITE(17,901) 72,PHASES(16),"PHI(M3H)"

      IF(IFAIL.GT.0.AND.IFAIL.LT.10) RETURN

      WRITE(17,899) "# "
      WRITE(17,899) "BLOCK MASS   # Mass spectrum "
      WRITE(17,899) "#  PDG Code     mass             particle "
      WRITE(17,902) 5,MB,"MB(MB)"
      WRITE(17,902) 6,MT,"MTOP (POLE MASS)"
      WRITE(17,902) 15,MTAU,"MTAU"
      WRITE(17,902) 23,MZ,"MZ"
      WRITE(17,902) 24,MW,"MW"
      WRITE(17,902) 25,DSQRT(MH0T(1)),"lightest neutral CPV_scalar"
      WRITE(17,902) 35,DSQRT(MH0T(2)),"second neutral CPV_scalar"
      WRITE(17,902) 45,DSQRT(MH0T(3)),"third neutral CPV_scalar"
      WRITE(17,902) 36,DSQRT(MH0T(4)),"fourth neutral CPV_scalar"
      WRITE(17,902) 46,DSQRT(MH0T(5)),"fifth neutral CPV_scalar"
      WRITE(17,902) 37,DSQRT(MHC_CPV),"charged Higgs"

      WRITE(17,902) 1000001,dsqrt(MSD2P(1))," ~d_L"
      WRITE(17,902) 2000001,dsqrt(MSD2P(2))," ~d_R"
      WRITE(17,902) 1000002,dsqrt(MSU2P(1))," ~u_L"
      WRITE(17,902) 2000002,dsqrt(MSU2P(2))," ~u_R"
      WRITE(17,902) 1000003,dsqrt(MSD2P(1))," ~s_L"
      WRITE(17,902) 2000003,dsqrt(MSD2P(2))," ~s_R"
      WRITE(17,902) 1000004,dsqrt(MSU2P(1))," ~c_L"
      WRITE(17,902) 2000004,dsqrt(MSU2P(2))," ~c_R"
      WRITE(17,902) 1000006,dsqrt(MSB2P(1))," ~b_1"
      WRITE(17,902) 2000006,dsqrt(MSB2P(2))," ~b_2"
      WRITE(17,902) 1000006,dsqrt(MST2P(1))," ~t_1"
      WRITE(17,902) 2000006,dsqrt(MST2P(2))," ~t_2"

      WRITE(17,902) 1000011,dsqrt(MSE2_CPV(1))," ~e_L"
      WRITE(17,902) 2000011,dsqrt(MSE2_CPV(2))," ~e_R"
      WRITE(17,902) 1000012,dsqrt(MSNE2_CPV)," ~nue_L"
      WRITE(17,902) 1000013,dsqrt(MSE2_CPV(1))," ~mu_L"
      WRITE(17,902) 2000013,dsqrt(MSE2_CPV(2))," ~mu_R"
      WRITE(17,902) 1000014,dsqrt(MSNE2_CPV)," ~numu_L"
      WRITE(17,902) 1000015,DSQRT(MSL2_CPV(1))," ~tau_1"
      WRITE(17,902) 2000015,DSQRT(MSL2_CPV(2))," ~tau_2"
      WRITE(17,902) 1000016,DSQRT(MSNT2_CPV)," ~nutau_L"
      WRITE(17,902) 1000021,MGL_CPV," ~g"
      WRITE(17,902) 1000022,DSQRT(MNEU_CPV(1)),"neutralino(1)"
      WRITE(17,902) 1000023,DSQRT(MNEU_CPV(2)),"neutralino(2)"
      WRITE(17,902) 1000025,DSQRT(MNEU_CPV(3)),"neutralino(3)"
      WRITE(17,902) 1000035,DSQRT(MNEU_CPV(4)),"neutralino(4)"
      WRITE(17,902) 1000045,DSQRT(MNEU_CPV(5)),"neutralino(5)"
      WRITE(17,902) 1000024,DSQRT(MCH2(1)),"chargino(1)"
      WRITE(17,902) 1000037,DSQRT(MCH2(2)),"chargino(2)"
      WRITE(17,899) "# "

       IF(OUTFLAG.NE.1) THEN
      WRITE(17,899) "# Low energy observables"
      WRITE(17,899) "BLOCK LOWEN"
      WRITE(17,899) "# "
      WRITE(17,899) "# EDMs:"
      WRITE(17,899) "# dEe < 1E-28 e.cm (Thorium exp.)"
      WRITE(17,907) "# ",dEe,"   dEe"
      WRITE(17,907) "# ",dEemin,"   |dEemin|: |dEe| - Theor.Err."
      WRITE(17,899) "# dETl < 1.3E-24 e.cm (Thalium exp.)"
      WRITE(17,907) "# ",dETl,"   dETl"
      WRITE(17,907) "# ",dETlmin,"   |dETlmin|: |dETl| - Theor.Err."
      WRITE(17,899) "# dEn < 3E-26 e.cm (neutron exp.)"
      WRITE(17,907) "# ",dEnCQM,"   dEnCQM"
      WRITE(17,907) "# ",dEnCQMmin,
     . "   |dEnCQMmin|: |dEnCQM| - Theor.Err."
      WRITE(17,907) "# ",dEnPQM,"   dEnPQM"
      WRITE(17,907) "# ",dEnPQMmin,
     . "   |dEnPQMmin|: |dEnPQM| - Theor.Err."
      WRITE(17,907) "# ",dEnQSR,"   dEnQSR"
      WRITE(17,907) "# ",dEnQSRmin,
     . "   |dEnQSRmin|: |dEnQSR| - Theor.Err."
      WRITE(17,899) "# dEHg < 3.1E-29 e.cm (Mercury exp.)"
      WRITE(17,907) "# ",dEHg,"   dEHg"
      WRITE(17,907) "# ",dEHgmin,"   |dEHgmin|: |dEHg| - Theor.Err."
      WRITE(17,899) "# "

      WRITE(17,899) "# BSM contr. to the muon anomalous magn. moment:"
      WRITE(17,901) 6,delmagmu,"Del_a_mu"
      WRITE(17,901) 61,amuthmax,"Del_a_mu + Theor.Err."
      WRITE(17,901) 62,amuthmin,"Del_a_mu - Theor.Err."
      WRITE(17,907) "# Minimal Exp.-SM (2 sigma):",damumin
      WRITE(17,907) "# Maximal Exp.-SM (2 sigma):",damumax
      WRITE(17,899) "# "

      WRITE(17,899)
     .   "# Exp. 2 Sigma: 3.02E-4 < BR(b -> s gamma) < 3.62E-4:"
      WRITE(17,901) 1,BRSG,"BR(b -> s gamma)"
      WRITE(17,901) 11,BRSGMAX,"(BR(b -> s gamma)+Theor.Err.)"
      WRITE(17,901) 12,BRSGMIN,"(BR(b -> s gamma)-Theor.Err.)"
      WRITE(17,899) "# Exp. 2 Sigma: 5.027E-1 < Delta M_d < 5.103E-1:"
      WRITE(17,901) 2,DMD,"Delta M_d in ps^-1"
      WRITE(17,901) 21,DMdmax,"Delta M_d +Theor.Err."
      WRITE(17,901) 22,DMdmin,"Delta M_d -Theor.Err."
      WRITE(17,899) 
     .   "# Exp. 2 Sigma: 1.7715E+1 < Delta Ms < 1.7799E+1:"
      WRITE(17,901) 3,DMS,"Delta M_s in ps^-1"
      WRITE(17,901) 31,DMsmax,"Delta M_s +Theor.Err."
      WRITE(17,901) 32,DMsmin,"Delta M_s -Theor.Err."
      WRITE(17,899) "# Exp. 2 Sigma: 1.7E-9 < BR(Bs->mu+mu-) < 4.5E-9:"
      WRITE(17,901) 4,BRBMUMU,"BR(Bs -> mu+mu-)"
      WRITE(17,901) 41,BRBMUMUmax,"BR(Bs -> mu+mu-)+Theor.Err."
      WRITE(17,901) 42,BRBMUMUmin,"BR(Bs -> mu+mu-)-Theor.Err."
      WRITE(17,899) 
     .   "# Exp. 2 Sigma: 0.78E-4 < BR(B+ > tau+ + nu_tau) < 1.44E-4:"
      WRITE(17,901) 5,BRBtaunu,"BR(B+ -> tau+ + nu_tau)"
      WRITE(17,901) 51,BRBtaunumax,
     .   "BR(B+ -> tau+ + nu_tau) + Theor.Err."
      WRITE(17,901) 52,BRBtaunumin,
     .   "BR(B+ -> tau+ + nu_tau) - Theor.Err."
      WRITE(17,899)
     .   "# Exp. 2 Sigma: 2.7E-6 < BR(b -> d gamma) < 2.55E-5:"
      WRITE(17,907) "# ",BRDG,"   BR(b -> d gamma)"
      WRITE(17,907) "# ",BRDGMAX,"   BR(b -> d gamma)+Theor.Err.)"
      WRITE(17,907) "# ",BRDGMIN,"   BR(b -> d gamma)-Theor.Err.)"
      WRITE(17,899) 
     .   "# Exp. 2 Sigma: 1.1E-10 < BR(Bd->mu+mu-) < 7.1E-10:"
      WRITE(17,907) "# ",BRBdMUMU,"   BR(Bd -> mu+mu-)"
      WRITE(17,907) "# ",BRBdMUMUmax,"   BR(Bd -> mu+mu-)+Theor.Err."
      WRITE(17,907) "# ",BRBdMUMUmin,"   BR(Bd -> mu+mu-)-Theor.Err."
      WRITE(17,899) "# RD = BR[B+ -> D tau+ nu_tau]/BR[B+ -> D l+ nu_l]"
      WRITE(17,899) 
     .   "# Exp. 2 Sigma: 2.80E-1 < RD < 4.00E-1:"
      WRITE(17,907) "# ",RD_taul,"   RD"
      WRITE(17,907) "# ",RD_taulmax,"   RD+Theor.Err."
      WRITE(17,907) "# ",RD_taulmin,"   RD-Theor.Err."
      WRITE(17,899)
     .    "# RD* = BR[B+ -> D* tau+ nu_tau]/BR[B+ -> D* l+ nu_l]"
      WRITE(17,899) 
     .   "# Exp. 2 Sigma: 2.67E-1 < RD* < 3.23E-1:"
      WRITE(17,907) "# ",RDs_taul,"   RD*"
      WRITE(17,907) "# ",RDs_taulmax,"   RD*+Theor.Err."
      WRITE(17,907) "# ",RDs_taulmin,"   RD*-Theor.Err."
      WRITE(17,899) 
*  From IF(OUTFLAG.NE.1)
       ENDIF

      WRITE(17,907) "BLOCK HMIX Q=",DSQRT(QSTSB),
     .    " # (STOP/SBOTTOM MASSES)"
      WRITE(17,901) 1,MUQ,"|MUEFF|"
      WRITE(17,901) 2,TANBQ,"TAN(BETA)"
      WRITE(17,901) 3,DSQRT(2d0*(H1Q**2+H2Q**2)),"V(Q)"
      WRITE(17,901) 4,PAR(23)**2,"MA^2"
      WRITE(17,901) 5,PAR(24)**2,"MP^2"

      WRITE(17,899) "# "
      WRITE(17,899) "# 5*5 Higgs mixing"
      WRITE(17,899) "BLOCK CVNMHMIX"
      WRITE(17,903) 1,1,XH0(1,2),"S_(1,1)"
      WRITE(17,903) 1,2,XH0(1,1),"S_(1,2)"
      WRITE(17,903) 1,3,XH0(1,3),"S_(1,3)"
      WRITE(17,903) 1,4,XH0(1,4),"S_(1,4)"
      WRITE(17,903) 1,5,XH0(1,5),"S_(1,5)"

      WRITE(17,903) 2,1,XH0(2,2),"S_(2,1)"
      WRITE(17,903) 2,2,XH0(2,1),"S_(2,2)"
      WRITE(17,903) 2,3,XH0(2,3),"S_(2,3)"
      WRITE(17,903) 2,4,XH0(2,4),"S_(2,4)"
      WRITE(17,903) 2,5,XH0(2,5),"S_(2,5)"

      WRITE(17,903) 3,1,XH0(3,2),"S_(3,1)"
      WRITE(17,903) 3,2,XH0(3,1),"S_(3,2)"
      WRITE(17,903) 3,3,XH0(3,3),"S_(3,3)"
      WRITE(17,903) 3,4,XH0(3,4),"S_(3,4)"
      WRITE(17,903) 3,5,XH0(3,5),"S_(3,5)"

      WRITE(17,903) 4,1,XH0(4,2),"S_(4,1)"
      WRITE(17,903) 4,2,XH0(4,1),"S_(4,2)"
      WRITE(17,903) 4,3,XH0(4,3),"S_(4,3)"
      WRITE(17,903) 4,4,XH0(4,4),"S_(4,4)"
      WRITE(17,903) 4,5,XH0(4,5),"S_(4,5)"

      WRITE(17,903) 5,1,XH0(5,2),"S_(5,1)"
      WRITE(17,903) 5,2,XH0(5,1),"S_(5,2)"
      WRITE(17,903) 5,3,XH0(5,3),"S_(5,3)"
      WRITE(17,903) 5,4,XH0(5,4),"S_(5,4)"
      WRITE(17,903) 5,5,XH0(5,5),"S_(5,5)"

      WRITE(17,899) "# "
      WRITE(17,899) "# 3rd generation sfermion mixing"
      WRITE(17,899) 
     .    "BLOCK STOPMIX  # Real part of the Stop mixing matrix"
      WRITE(17,903) 1,1,UT(1,1,1),"UT(1,1,1)"
      WRITE(17,903) 1,2,UT(1,2,1),"UT(1,2,1)"
      WRITE(17,903) 2,1,UT(2,1,1),"UT(2,1,1)"
      WRITE(17,903) 2,2,UT(2,2,1),"UT(2,2,1)"

      WRITE(17,899) 
     .    "BLOCK IMSTOPMIX  # Imag. part of the Stop mixing matrix"
      WRITE(17,903) 1,1,UT(1,1,2),"UT(1,1,2)"
      WRITE(17,903) 1,2,UT(1,2,2),"UT(1,2,2)"
      WRITE(17,903) 2,1,UT(2,1,2),"UT(2,1,2)"
      WRITE(17,903) 2,2,UT(2,2,2),"UT(2,2,2)"

      WRITE(17,899) 
     .     "BLOCK SBOTMIX  # Real part of the Sbottom mixing matrix"
      WRITE(17,903) 1,1,UB(1,1,1),"UB(1,1,1)"
      WRITE(17,903) 1,2,UB(1,2,1),"UB(1,2,1)"
      WRITE(17,903) 2,1,UB(2,1,1),"UB(2,1,1)"
      WRITE(17,903) 2,2,UB(2,2,1),"UB(2,2,1)"

      WRITE(17,899) 
     .    "BLOCK IMSBOTMIX  # Imag. part of the Sbottom mixing matrix"
      WRITE(17,903) 1,1,UB(1,1,2),"UB(1,1,2)"
      WRITE(17,903) 1,2,UB(1,2,2),"UB(1,2,2)"
      WRITE(17,903) 2,1,UB(2,1,2),"UB(2,1,2)"
      WRITE(17,903) 2,2,UB(2,2,2),"UB(2,2,2)"

      WRITE(17,899) 
     .    "BLOCK STAUMIX  # Real part of the Stau mixing matrix"
      WRITE(17,903) 1,1,UTAU(1,1,1),"UTAU(1,1,1)"
      WRITE(17,903) 1,2,UTAU(1,2,1),"UTAU(1,2,1)"
      WRITE(17,903) 2,1,UTAU(2,1,1),"UTAU(2,1,1)"
      WRITE(17,903) 2,2,UTAU(2,2,1),"UTAU(2,2,1)"

      WRITE(17,899) 
     .    "BLOCK STAUMIX  # Imag. part of the Stau mixing matrix"
      WRITE(17,903) 1,1,UTAU(1,1,2),"UTAU(1,1,2)"
      WRITE(17,903) 1,2,UTAU(1,2,2),"UTAU(1,2,2)"
      WRITE(17,903) 2,1,UTAU(2,1,2),"UTAU(2,1,2)"
      WRITE(17,903) 2,2,UTAU(2,2,2),"UTAU(2,2,2)"

      WRITE(17,899) "# "
      WRITE(17,899) "# Gaugino-Higgsino mixing"
      WRITE(17,899) 
     .   "BLOCK RNMNMIX  # Real part of the Neutralino Mixing Matrix"
      WRITE(17,903) 1,1,NEU_CPV(1,1,1),"RE(N)_(1,1)"
      WRITE(17,903) 1,2,NEU_CPV(1,2,1),"RE(N)_(1,2)"
      WRITE(17,903) 1,3,NEU_CPV(1,4,1),"RE(N)_(1,3)"
      WRITE(17,903) 1,4,NEU_CPV(1,3,1),"RE(N)_(1,4)"
      WRITE(17,903) 1,5,NEU_CPV(1,5,1),"RE(N)_(1,5)"
      WRITE(17,903) 2,1,NEU_CPV(2,1,1),"RE(N)_(2,1)"
      WRITE(17,903) 2,2,NEU_CPV(2,2,1),"RE(N)_(2,2)"
      WRITE(17,903) 2,3,NEU_CPV(2,4,1),"RE(N)_(2,3)"
      WRITE(17,903) 2,4,NEU_CPV(2,3,1),"RE(N)_(2,4)"
      WRITE(17,903) 2,5,NEU_CPV(2,5,1),"RE(N)_(2,5)"
      WRITE(17,903) 3,1,NEU_CPV(3,1,1),"RE(N)_(3,1)"
      WRITE(17,903) 3,2,NEU_CPV(3,2,1),"RE(N)_(3,2)"
      WRITE(17,903) 3,3,NEU_CPV(3,4,1),"RE(N)_(3,3)"
      WRITE(17,903) 3,4,NEU_CPV(3,3,1),"RE(N)_(3,4)"
      WRITE(17,903) 3,5,NEU_CPV(3,5,1),"RE(N)_(3,5)"
      WRITE(17,903) 4,1,NEU_CPV(4,1,1),"RE(N)_(4,1)"
      WRITE(17,903) 4,2,NEU_CPV(4,2,1),"RE(N)_(4,2)"
      WRITE(17,903) 4,3,NEU_CPV(4,4,1),"RE(N)_(4,3)"
      WRITE(17,903) 4,4,NEU_CPV(4,3,1),"RE(N)_(4,4)"
      WRITE(17,903) 4,5,NEU_CPV(4,5,1),"RE(N)_(4,5)"
      WRITE(17,903) 5,1,NEU_CPV(5,1,1),"RE(N)_(5,1)"
      WRITE(17,903) 5,2,NEU_CPV(5,2,1),"RE(N)_(5,2)"
      WRITE(17,903) 5,3,NEU_CPV(5,4,1),"RE(N)_(5,3)"
      WRITE(17,903) 5,4,NEU_CPV(5,3,1),"RE(N)_(5,4)"
      WRITE(17,903) 5,5,NEU_CPV(5,5,1),"RE(N)_(5,5)"

      WRITE(17,899) "# "
      WRITE(17,899) "# Gaugino-Higgsino mixing"
      WRITE(17,899) 
     .     "BLOCK INMNMIX  # Imag. part of the Neutralino Mixing Matrix"
      WRITE(17,903) 1,1,NEU_CPV(1,1,2),"IM(N)_(1,1)"
      WRITE(17,903) 1,2,NEU_CPV(1,2,2),"IM(N)_(1,2)"
      WRITE(17,903) 1,3,NEU_CPV(1,4,2),"IM(N)_(1,3)"
      WRITE(17,903) 1,4,NEU_CPV(1,3,2),"IM(N)_(1,4)"
      WRITE(17,903) 1,5,NEU_CPV(1,5,2),"IM(N)_(1,5)"
      WRITE(17,903) 2,1,NEU_CPV(2,1,2),"IM(N)_(2,1)"
      WRITE(17,903) 2,2,NEU_CPV(2,2,2),"IM(N)_(2,2)"
      WRITE(17,903) 2,3,NEU_CPV(2,4,2),"IM(N)_(2,3)"
      WRITE(17,903) 2,4,NEU_CPV(2,3,2),"IM(N)_(2,4)"
      WRITE(17,903) 2,5,NEU_CPV(2,5,2),"IM(N)_(2,5)"
      WRITE(17,903) 3,1,NEU_CPV(3,1,2),"IM(N)_(3,1)"
      WRITE(17,903) 3,2,NEU_CPV(3,2,2),"IM(N)_(3,2)"
      WRITE(17,903) 3,3,NEU_CPV(3,4,2),"IM(N)_(3,3)"
      WRITE(17,903) 3,4,NEU_CPV(3,3,2),"IM(N)_(3,4)"
      WRITE(17,903) 3,5,NEU_CPV(3,5,2),"IM(N)_(3,5)"
      WRITE(17,903) 4,1,NEU_CPV(4,1,2),"IM(N)_(4,1)"
      WRITE(17,903) 4,2,NEU_CPV(4,2,2),"IM(N)_(4,2)"
      WRITE(17,903) 4,3,NEU_CPV(4,4,2),"IM(N)_(4,3)"
      WRITE(17,903) 4,4,NEU_CPV(4,3,2),"IM(N)_(4,4)"
      WRITE(17,903) 4,5,NEU_CPV(4,5,2),"IM(N)_(4,5)"
      WRITE(17,903) 5,1,NEU_CPV(5,1,2),"IM(N)_(5,1)"
      WRITE(17,903) 5,2,NEU_CPV(5,2,2),"IM(N)_(5,2)"
      WRITE(17,903) 5,3,NEU_CPV(5,4,2),"IM(N)_(5,3)"
      WRITE(17,903) 5,4,NEU_CPV(5,3,2),"IM(N)_(5,4)"
      WRITE(17,903) 5,5,NEU_CPV(5,5,2),"IM(N)_(5,5)"

      WRITE(17,899) "# "
      WRITE(17,899) "BLOCK RUMIX  # Real Chargino U Mixing Matrix"
      WRITE(17,903) 1,1,U_CPV(1,1,1),"RE(U)_(1,1)"
      WRITE(17,903) 1,2,U_CPV(1,2,1),"RE(U)_(1,2)"
      WRITE(17,903) 2,1,U_CPV(2,1,1),"RE(U)_(2,1)"
      WRITE(17,903) 2,2,U_CPV(2,2,1),"RE(U)_(2,2)"

      WRITE(17,899) "# "
      WRITE(17,899) "BLOCK IUMIX  # Imag Chargino U Mixing Matrix"
      WRITE(17,903) 1,1,U_CPV(1,1,2),"IM(U)_(1,1)"
      WRITE(17,903) 1,2,U_CPV(1,2,2),"IM(U)_(1,2)"
      WRITE(17,903) 2,1,U_CPV(2,1,2),"IM(U)_(2,1)"
      WRITE(17,903) 2,2,U_CPV(2,2,2),"IM(U)_(2,2)"

      WRITE(17,899) "# "
      WRITE(17,899) "BLOCK RVMIX  # Real Chargino V Mixing Matrix"
      WRITE(17,903) 1,1,V_CPV(1,1,1),"RE(V)_(1,1)"
      WRITE(17,903) 1,2,V_CPV(1,2,1),"RE(V)_(1,2)"
      WRITE(17,903) 2,1,V_CPV(2,1,1),"RE(V)_(2,1)"
      WRITE(17,903) 2,2,V_CPV(2,2,1),"RE(V)_(2,2)"

      WRITE(17,899) "# "
      WRITE(17,899) "BLOCK IVMIX  # Imag Chargino V Mixing Matrix"
      WRITE(17,903) 1,1,V_CPV(1,1,2),"IM(V)_(1,1)"
      WRITE(17,903) 1,2,V_CPV(1,2,2),"IM(V)_(1,2)"
      WRITE(17,903) 2,1,V_CPV(2,1,2),"IM(V)_(2,1)"
      WRITE(17,903) 2,2,V_CPV(2,2,2),"IM(V)_(2,2)"
      WRITE(17,899) "# "

      WRITE(17,899) "# SM-Higgs reduced couplings"
      WRITE(17,899) "# (as compared to a SM Higgs with same mass)"
      WRITE(17,899) "# Allowed 2sigma intervals combining"
      WRITE(17,'(A,A)') "# 1606.02266 Table 17, 1809.10733 Table 8,",
     .  " ATLAS-CONF-2021-053 Table 8(b)"
      WRITE(17,899) "# 0.841 < CU < 1.161 " ! top quark
      WRITE(17,899) "# 0.655 < CB < 0.917 " ! b quark
      WRITE(17,899) "# 0.788 < CL < 0.984 " ! taus
      WRITE(17,899) "# 0.933 < CV < 1.000 " ! electroweak couplings
      WRITE(17,899) "# 0.810 < CJ < 1.022 " ! gluons
      WRITE(17,899) "# 0.896 < CG < 1.055 " ! photons

      WRITE(17,899) "BLOCK REDCOUP"

      WRITE(17,899) "# H1"
      WRITE(17,903) 1,1,CU(1),"U-type fermions"
      WRITE(17,903) 1,2,CD(1),"D-type fermions"
      WRITE(17,903) 1,3,CB(1),"b-quarks"
      WRITE(17,903) 1,4,CL(1),"taus"
      WRITE(17,903) 1,5,CV(1),"W,Z bosons"
      WRITE(17,903) 1,6,CJ(1),"Gluons"
      WRITE(17,903) 1,7,CG(1),"Photons"
      WRITE(17,899) "# H1, imag. parts:"
      WRITE(17,903) 1,1,CUP(1),"U-type fermions"
      WRITE(17,903) 1,2,CDP(1),"D-type fermions"
      WRITE(17,903) 1,3,CBP(1),"b-quarks"
      WRITE(17,903) 1,4,CLP(1),"taus"
      WRITE(17,903) 1,5,0d0,"W,Z bosons"
      WRITE(17,903) 1,6,CJP(1),"Gluons"
      WRITE(17,903) 1,7,CGP(1),"Photons"

      WRITE(17,899) "# H2"
      WRITE(17,903) 2,1,CU(2),"U-type fermions"
      WRITE(17,903) 2,2,CD(2),"D-type fermions"
      WRITE(17,903) 2,3,CB(2),"b-quarks"
      WRITE(17,903) 1,4,CL(2),"taus"
      WRITE(17,903) 2,5,CV(2),"W,Z bosons"
      WRITE(17,903) 2,6,CJ(2),"Gluons"
      WRITE(17,903) 2,7,CG(2),"Photons"
      WRITE(17,899) "# H2, imag. parts"
      WRITE(17,903) 2,1,CUP(2),"U-type fermions"
      WRITE(17,903) 2,2,CDP(2),"D-type fermions"
      WRITE(17,903) 2,3,CBP(2),"b-quarks"
      WRITE(17,903) 1,4,CLP(2),"taus"
      WRITE(17,903) 1,5,0d0,"W,Z bosons"
      WRITE(17,903) 2,6,CJP(2),"Gluons"
      WRITE(17,903) 2,7,CGP(2),"Photons"

      WRITE(17,899) "# H3"
      WRITE(17,903) 3,1,CU(3),"U-type fermions"
      WRITE(17,903) 3,2,CD(3),"D-type fermions"
      WRITE(17,903) 3,3,CB(3),"b-quarks"
      WRITE(17,903) 1,4,CL(3),"taus"
      WRITE(17,903) 3,5,CV(3),"W,Z bosons"
      WRITE(17,903) 3,6,CJ(3),"Gluons"
      WRITE(17,903) 3,7,CG(3),"Photons"
      WRITE(17,899) "# H3, imag. parts"
      WRITE(17,903) 3,1,CUP(3),"U-type fermions"
      WRITE(17,903) 3,2,CDP(3),"D-type fermions"
      WRITE(17,903) 3,3,CBP(3),"b-quarks"
      WRITE(17,903) 1,4,CLP(3),"taus"
      WRITE(17,903) 1,5,0d0,"W,Z bosons"
      WRITE(17,903) 3,6,CJP(3),"Gluons"
      WRITE(17,903) 3,7,CGP(3),"Photons"

      WRITE(17,899) "# H4"
      WRITE(17,903) 4,1,CU(4),"U-type fermions"
      WRITE(17,903) 4,2,CD(4),"D-type fermions"
      WRITE(17,903) 4,3,CB(4),"b-quarks"
      WRITE(17,903) 1,4,CL(4),"taus"
      WRITE(17,903) 4,5,CV(4),"W,Z bosons"
      WRITE(17,903) 4,6,CJ(4),"Gluons"
      WRITE(17,903) 4,7,CG(4),"Photons"
      WRITE(17,899) "# H4, imag. parts"
      WRITE(17,903) 4,1,CUP(4),"U-type fermions"
      WRITE(17,903) 4,2,CDP(4),"D-type fermions"
      WRITE(17,903) 4,3,CBP(4),"b-quarks"
      WRITE(17,903) 1,4,CLP(4),"taus"
      WRITE(17,903) 1,5,0d0,"W,Z bosons"
      WRITE(17,903) 4,6,CJP(4),"Gluons"
      WRITE(17,903) 4,7,CGP(4),"Photons"

      WRITE(17,899) "# H5"
      WRITE(17,903) 5,1,CU(5),"U-type fermions"
      WRITE(17,903) 5,2,CD(5),"D-type fermions"
      WRITE(17,903) 5,3,CB(5),"b-quarks"
      WRITE(17,903) 1,4,CL(5),"taus"
      WRITE(17,903) 5,5,CV(5),"W,Z bosons"
      WRITE(17,903) 5,6,CJ(5),"Gluons"
      WRITE(17,903) 5,7,CG(5),"Photons"
      WRITE(17,899) "# H5, imag. parts"
      WRITE(17,903) 5,1,CUP(5),"U-type fermions"
      WRITE(17,903) 5,2,CDP(5),"D-type fermions"
      WRITE(17,903) 5,3,CBP(5),"b-quarks"
      WRITE(17,903) 1,4,CLP(5),"taus"
      WRITE(17,903) 1,5,0d0,"W,Z bosons"
      WRITE(17,903) 5,6,CJP(5),"Gluons"
      WRITE(17,903) 5,7,CGP(5),"Photons"

      WRITE(17,899) "# "
      WRITE(17,899) "# GAUGE AND YUKAWA COUPLINGS AT THE SUSY SCALE"
      WRITE(17,907) "BLOCK GAUGE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,901) 1,DSQRT(G1S),"g1(Q,DR_bar)"
      WRITE(17,901) 2,DSQRT(G2S),"g2(Q,DR_bar)"
      WRITE(17,901) 3,DSQRT(G3S),"g3(Q,DR_bar)"

      WRITE(17,907) "BLOCK YU Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 3,3,HTOPS,"HTOP(Q,DR_bar)"
      WRITE(17,907) "BLOCK YD Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 3,3,HBOTS,"HBOT(Q,DR_bar)"
      WRITE(17,907) "BLOCK YE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 3,3,HTAUS,"HTAU(Q,DR_bar)"

      WRITE(17,899) "# "
      WRITE(17,899) "# SOFT TRILINEAR COUPLINGS AT THE SUSY SCALE"
      WRITE(17,899) "# (BOTH SLHA1 AND SLHA2 FORMAT)"
      WRITE(17,907) "BLOCK AU Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 3,3,PAR(12),"ATOP"
      WRITE(17,907) "BLOCK TU Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 3,3,PAR(12),"ATOP"
      WRITE(17,907) "BLOCK AD Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 3,3,PAR(13),"ABOT"
      WRITE(17,907) "BLOCK TD Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 3,3,PAR(13),"ABOT"
      WRITE(17,907) "BLOCK AE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 2,2,PAR(25),"AMUON"
      WRITE(17,903) 3,3,PAR(14),"ATAU"
      WRITE(17,907) "BLOCK TE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 2,2,PAR(25),"AMUON"
      WRITE(17,903) 3,3,PAR(14),"ATAU"

      WRITE(17,899) "# "
      WRITE(17,899) "# SOFT MASSES AT THE SUSY SCALE"
      WRITE(17,907) "BLOCK MSOFT Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,901) 1,PAR(20),"M1"
      WRITE(17,901) 2,PAR(21),"M2"
      WRITE(17,901) 3,PAR(22),"M3"
      WRITE(17,901) 21,MHDS,"M_HD^2"
      WRITE(17,901) 22,MHUS,"M_HU^2"
      WRITE(17,901) 31,PAR(18)/DSQRT(DABS(PAR(18))),"M_eL"
      WRITE(17,901) 32,PAR(18)/DSQRT(DABS(PAR(18))),"M_muL"
      WRITE(17,901) 33,PAR(10)/DSQRT(DABS(PAR(10))),"M_tauL"
      WRITE(17,901) 34,PAR(19)/DSQRT(DABS(PAR(19))),"M_eR"
      WRITE(17,901) 35,PAR(19)/DSQRT(DABS(PAR(19))),"M_muR"
      WRITE(17,901) 36,PAR(11)/DSQRT(DABS(PAR(11))),"M_tauR"
      WRITE(17,901) 41,PAR(15)/DSQRT(DABS(PAR(15))),"M_q1L"
      WRITE(17,901) 42,PAR(15)/DSQRT(DABS(PAR(15))),"M_q2L"
      WRITE(17,901) 43,PAR(7)/DSQRT(DABS(PAR(7))),"M_q3L"
      WRITE(17,901) 44,PAR(16)/DSQRT(DABS(PAR(16))),"M_uR"
      WRITE(17,901) 45,PAR(16)/DSQRT(DABS(PAR(16))),"M_cR"
      WRITE(17,901) 46,PAR(8)/DSQRT(DABS(PAR(8))),"M_tR"
      WRITE(17,901) 47,PAR(17)/DSQRT(DABS(PAR(17))),"M_dR"
      WRITE(17,901) 48,PAR(17)/DSQRT(DABS(PAR(17))),"M_sR"
      WRITE(17,901) 49,PAR(9)/DSQRT(DABS(PAR(9))),"M_bR"

      WRITE(17,899) "# "
      WRITE(17,899) "# REDUCED CROSS SECTIONS AT LHC"
      WRITE(17,899) "BLOCK LHCCROSSSECTIONS"
      WRITE(17,901) 11,SIG_CPV(1,1),"VBF/VH -> H1 -> tautau"
      WRITE(17,901) 12,SIG_CPV(1,2),"ggF -> H1 -> tautau"
      WRITE(17,901) 13,SIG_CPV(1,3),"VBF/VH -> H1 -> bb"
      WRITE(17,901) 14,SIG_CPV(1,4),"ttH -> H1 -> bb"
      WRITE(17,901) 15,SIG_CPV(1,5),"VBF/VH -> H1 -> ZZ/WW"
      WRITE(17,901) 16,SIG_CPV(1,6),"ggF -> H1 -> ZZ/WW"
      WRITE(17,901) 17,SIG_CPV(1,7),"VBF/VH -> H1 -> gammagamma"
      WRITE(17,901) 18,SIG_CPV(1,8),"ggF -> H1 -> gammagamma"
      WRITE(17,901) 21,SIG_CPV(2,1),"VBF/VH -> H2 -> tautau"
      WRITE(17,901) 22,SIG_CPV(2,2),"ggF -> H2 -> tautau"
      WRITE(17,901) 23,SIG_CPV(2,3),"VBF/VH -> H2 -> bb"
      WRITE(17,901) 24,SIG_CPV(2,4),"ttH -> H2 -> bb"
      WRITE(17,901) 25,SIG_CPV(2,5),"VBF/VH -> H2 -> ZZ/WW"
      WRITE(17,901) 26,SIG_CPV(2,6),"ggF -> H2 -> ZZ/WW"
      WRITE(17,901) 27,SIG_CPV(2,7),"VBF/VH -> H2 -> gammagamma"
      WRITE(17,901) 28,SIG_CPV(2,8),"ggF -> H2 -> gammagamma"
      WRITE(17,901) 31,SIG_CPV(3,1),"VBF/VH -> H3 -> tautau"
      WRITE(17,901) 32,SIG_CPV(3,2),"ggF -> H3 -> tautau"
      WRITE(17,901) 33,SIG_CPV(3,3),"VBF/VH -> H3 -> bb"
      WRITE(17,901) 34,SIG_CPV(3,4),"ttH -> H3 -> bb"
      WRITE(17,901) 35,SIG_CPV(3,5),"VBF/VH -> H3 -> ZZ/WW"
      WRITE(17,901) 36,SIG_CPV(3,6),"ggF -> H3 -> ZZ/WW"
      WRITE(17,901) 37,SIG_CPV(3,7),"VBF/VH -> H3 -> gammagamma"
      WRITE(17,901) 38,SIG_CPV(3,8),"ggF -> H3 -> gammagamma"
      WRITE(17,901) 41,SIG_CPV(4,1),"VBF/VH -> H4 -> tautau"
      WRITE(17,901) 42,SIG_CPV(4,2),"ggF -> H4 -> tautau"
      WRITE(17,901) 43,SIG_CPV(4,3),"VBF/VH -> H4 -> bb"
      WRITE(17,901) 44,SIG_CPV(4,4),"ttH -> H4 -> bb"
      WRITE(17,901) 45,SIG_CPV(4,5),"VBF/VH -> H4 -> ZZ/WW"
      WRITE(17,901) 46,SIG_CPV(4,6),"ggF -> H4 -> ZZ/WW"
      WRITE(17,901) 47,SIG_CPV(4,7),"VBF/VH -> H4 -> gammagamma"
      WRITE(17,901) 48,SIG_CPV(4,8),"ggF -> H4 -> gammagamma"
      WRITE(17,901) 51,SIG_CPV(5,1),"VBF/VH -> H5 -> tautau"
      WRITE(17,901) 52,SIG_CPV(5,2),"ggF -> H5 -> tautau"
      WRITE(17,901) 53,SIG_CPV(5,3),"VBF/VH -> H5 -> bb"
      WRITE(17,901) 54,SIG_CPV(5,4),"ttH -> H5 -> bb"
      WRITE(17,901) 55,SIG_CPV(5,5),"VBF/VH -> H5 -> ZZ/WW"
      WRITE(17,901) 56,SIG_CPV(5,6),"ggF -> H5 -> ZZ/WW"
      WRITE(17,901) 57,SIG_CPV(5,7),"VBF/VH -> H5 -> gammagamma"
      WRITE(17,901) 58,SIG_CPV(5,8),"ggF -> H5 -> gammagamma"

      WRITE(18,899) "# "
      WRITE(18,899) "# HIGGS + TOP BRANCHING RATIOS IN SLHA FORMAT"
      WRITE(17,899) "# Allowed BSM branching ratio for the most SM-like"
      WRITE(17,899) "# Higgs state in the 125.1+/-3 GeV window (LHC):"
      WRITE(17,899) "# BR_BSM < 0.194 "
      WRITE(18,899) "# Info about decay package"
      WRITE(18,899) "BLOCK DCINFO   # Program information"
      WRITE(18,900) 1,"NMSSMTools # Decay package"
      WRITE(18,900) 2,"5.6.2      # Version number"

      WRITE(18,899) "#           PDG          Width"
      WRITE(18,904) 25,WIDTH(1),"Lightest neutral Higgs scalar"
      IF(BRJJ(1).GT.0d0)
     .  WRITE(18,905) BRJJ(1),2,21,21,"BR(H_1 -> hadrons)"
      IF(BRMM(1).GT.0d0)
     .  WRITE(18,905) BRMM(1),2,13,-13,"BR(H_1 -> muon muon)"
      IF(BRLL(1).GT.0d0)
     .  WRITE(18,905) BRLL(1),2,15,-15,"BR(H_1 -> tau tau)"
      IF(BRCC(1).GT.0d0)
     .  WRITE(18,905) BRCC(1),2,4,-4,"BR(H_1 -> c cbar)"
      IF(BRBB(1).GT.0d0)
     .  WRITE(18,905) BRBB(1),2,5,-5,"BR(H_1 -> b bbar)"
      IF(BRTT(1).GT.0d0)
     .  WRITE(18,905) BRTT(1),2,6,-6,"BR(H_1 -> t tbar)"
      IF(BRWW(1).GT.0d0)
     .  WRITE(18,905) BRWW(1),2,24,-24,"BR(H_1 -> W+ W-)"
      IF(BRZZ(1).GT.0d0)
     .  WRITE(18,905) BRZZ(1),2,23,23,"BR(H_1 -> Z Z)"
      IF(BRGG(1).GT.0d0)
     .  WRITE(18,905) BRGG(1),2,22,22,"BR(H_1 -> gamma gamma)"
      IF(BRZG(1).GT.0d0)
     .  WRITE(18,905) BRZG(1),2,23,22,"BR(H_1 -> Z gamma)"
      IF(BRHCHC(1).GT.0d0)
     .  WRITE(18,905) BRHCHC(1),2,37,-37,"BR(H_1 -> H+ H-)"
      IF(BRHCW(1).GT.0d0)
     .  WRITE(18,905) BRHCW(1),2,24,-37,"BR(H_1 -> W+ H-)"
      IF(BRHCW(1).GT.0d0)
     .  WRITE(18,905) BRHCW(1),2,-24,37,"BR(H_1 -> W- H+)"
      IF(BRNEU(1,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,1),2,1000022,1000022,
     .    "BR(H_1 -> neu_1 neu_1)"
      IF(BRNEU(1,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,2),2,1000022,1000023,
     .    "BR(H_1 -> neu_1 neu_2)"
      IF(BRNEU(1,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,3),2,1000022,1000025,
     .    "BR(H_1 -> neu_1 neu_3)"
      IF(BRNEU(1,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,4),2,1000022,1000035,
     .    "BR(H_1 -> neu_1 neu_4)"
      IF(BRNEU(1,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,5),2,1000022,1000045,
     .    "BR(H_1 -> neu_1 neu_5)"
      IF(BRNEU(1,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(1,2,2),2,1000023,1000023,
     .    "BR(H_1 -> neu_2 neu_2)"
      IF(BRNEU(1,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(1,2,3),2,1000023,1000025,
     .    "BR(H_1 -> neu_2 neu_3)"
      IF(BRNEU(1,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(1,2,4),2,1000023,1000035,
     .    "BR(H_1 -> neu_2 neu_4)"
      IF(BRNEU(1,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,2,5),2,1000023,1000045,
     .    "BR(H_1 -> neu_2 neu_5)"
      IF(BRNEU(1,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(1,3,3),2,1000025,1000025,
     .    "BR(H_1 -> neu_3 neu_3)"
      IF(BRNEU(1,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(1,3,4),2,1000025,1000035,
     .    "BR(H_1 -> neu_3 neu_4)"
      IF(BRNEU(1,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,3,5),2,1000025,1000045,
     .    "BR(H_1 -> neu_3 neu_5)"
      IF(BRNEU(1,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(1,4,4),2,1000035,1000035,
     .    "BR(H_1 -> neu_4 neu_4)"
      IF(BRNEU(1,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,4,5),2,1000035,1000045,
     .    "BR(H_1 -> neu_4 neu_5)"
      IF(BRNEU(1,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,5,5),2,1000045,1000045,
     .    "BR(H_1 -> neu_5 neu_5)"
      IF(BRCHA(1,1).GT.0d0)
     .  WRITE(18,905) BRCHA(1,1),2,1000024,-1000024,
     .    "BR(H_1 -> cha_1 cha_1bar)"
      IF(BRCHA(1,2).GT.0d0)
     .  WRITE(18,905) BRCHA(1,2),2,1000024,-1000037,
     .    "BR(H_1 -> cha_1 cha_2bar)"
      IF(BRCHA(1,2).GT.0d0)
     .  WRITE(18,905) BRCHA(1,2),2,1000037,-1000024,
     .    "BR(H_1 -> cha_2 cha_1bar)"
      IF(BRCHA(1,3).GT.0d0)
     .  WRITE(18,905) BRCHA(1,3),2,1000037,-1000037,
     .    "BR(H_1 -> cha_2 cha_2bar)"
      IF(BRHSQ(1,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,1),2,1000002,-1000002,
     .    "BR(H_1 -> ~u_L ~ubar_L)"
      IF(BRHSQ(1,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,1),2,1000004,-1000004,
     .    "BR(H_1 -> ~c_L ~cbar_L)"
      IF(BRHSQ(1,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,2),2,2000002,-2000002,
     .    "BR(H_1 -> ~u_R ~ubar_R)"
      IF(BRHSQ(1,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,2),2,2000004,-2000004,
     .    "BR(H_1 -> ~c_R ~cbar_R)"
      IF(BRHSQ(1,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,3),2,1000001,-1000001,
     .    "BR(H_1 -> ~d_L ~dbar_L)"
      IF(BRHSQ(1,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,3),2,1000003,-1000003,
     .    "BR(H_1 -> ~s_L ~sbar_L)"
      IF(BRHSQ(1,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,4),2,2000001,-2000001,
     .    "BR(H_1 -> ~d_R ~dbar_R)"
      IF(BRHSQ(1,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,4),2,2000003,-2000003,
     .    "BR(H_1 -> ~s_R ~sbar_R)"
      IF(BRHSQ(1,5).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,5),2,1000006,-1000006,
     .    "BR(H_1 -> ~t_1 ~tbar_1)"
      IF(BRHSQ(1,6).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,6),2,2000006,-2000006,
     .    "BR(H_1 -> ~t_2 ~tbar_2)"
      IF(BRHSQ(1,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,7),2,1000006,-2000006,
     .    "BR(H_1 -> ~t_1 ~tbar_2)"
      IF(BRHSQ(1,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,7),2,2000006,-1000006,
     .    "BR(H_1 -> ~t_2 ~tbar_1)"
      IF(BRHSQ(1,8).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,8),2,1000005,-1000005,
     .    "BR(H_1 -> ~b_1 ~bbar_1)"
      IF(BRHSQ(1,9).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,9),2,2000005,-2000005,
     .    "BR(H_1 -> ~b_2 ~bbar_2)"
      IF(BRHSQ(1,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,10),2,1000005,-2000005,
     .    "BR(H_1 -> ~b_1 ~bbar_2)"
      IF(BRHSQ(1,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,10),2,2000005,-1000005,
     .    "BR(H_1 -> ~b_2 ~bbar_1)"
      IF(BRHSL(1,1).GT.0d0)
     .  WRITE(18,905) BRHSL(1,1),2,1000011,-1000011,
     .    "BR(H_1 -> ~e_L ~ebar_L)"
      IF(BRHSL(1,1).GT.0d0)
     .  WRITE(18,905) BRHSL(1,1),2,1000013,-1000013,
     .    "BR(H_1 -> ~mu_L ~mubar_L)"
      IF(BRHSL(1,2).GT.0d0)
     .  WRITE(18,905) BRHSL(1,2),2,2000011,-2000011,
     .    "BR(H_1 -> ~e_R ~ebar_R)"
      IF(BRHSL(1,2).GT.0d0)
     .  WRITE(18,905) BRHSL(1,2),2,2000013,-2000013,
     .    "BR(H_1 -> ~mu_R ~mubarRL)"
      IF(BRHSL(1,3).GT.0d0)
     .  WRITE(18,905) BRHSL(1,3),2,1000012,-1000012,
     .    "BR(H_1 -> ~nu_e_L ~nu_ebar_L)"
      IF(BRHSL(1,3).GT.0d0)
     .  WRITE(18,905) BRHSL(1,3),2,1000014,-1000014,
     .    "BR(H_1 -> ~nu_mu_L ~nu_mubar_L)"
      IF(BRHSL(1,4).GT.0d0)
     .  WRITE(18,905) BRHSL(1,4),2,1000015,-1000015,
     .    "BR(H_1 -> ~tau_1 ~taubar_1)"
      IF(BRHSL(1,5).GT.0d0)
     .  WRITE(18,905) BRHSL(1,5),2,2000015,-2000015,
     .    "BR(H_1 -> ~tau_2 ~taubar_2)"
      IF(BRHSL(1,6).GT.0d0)
     .  WRITE(18,905) BRHSL(1,6),2,1000015,-2000015,
     .    "BR(H_1 -> ~tau_1 ~taubar_2)"
      IF(BRHSL(1,6).GT.0d0)
     .  WRITE(18,905) BRHSL(1,6),2,2000015,-1000015,
     .    "BR(H_1 -> ~tau_2 ~taubar_1)"
      IF(BRHSL(1,7).GT.0d0)
     .  WRITE(18,905) BRHSL(1,7),2,1000016,-1000016,
     .    "BR(H_1 -> ~nu_tau_L ~nu_taubar_L)"

      WRITE(18,904) 35,WIDTH(2),"2nd neutral Higgs scalar"
      IF(BRJJ(2).GT.0d0)
     .  WRITE(18,905) BRJJ(2),2,21,21,"BR(H_2 -> hadrons)"
      IF(BRMM(2).GT.0d0)
     .  WRITE(18,905) BRMM(2),2,13,-13,"BR(H_2 -> muon muon)"
      IF(BRLL(2).GT.0d0)
     .  WRITE(18,905) BRLL(2),2,15,-15,"BR(H_2 -> tau tau)"
      IF(BRCC(2).GT.0d0)
     .  WRITE(18,905) BRCC(2),2,4,-4,"BR(H_2 -> c cbar)"
      IF(BRBB(2).GT.0d0)
     .  WRITE(18,905) BRBB(2),2,5,-5,"BR(H_2 -> b bbar)"
      IF(BRTT(2).GT.0d0)
     .  WRITE(18,905) BRTT(2),2,6,-6,"BR(H_2 -> t tbar)"
      IF(BRWW(2).GT.0d0)
     .  WRITE(18,905) BRWW(2),2,24,-24,"BR(H_2 -> W+ W-)"
      IF(BRZZ(2).GT.0d0)
     .  WRITE(18,905) BRZZ(2),2,23,23,"BR(H_2 -> Z Z)"
      IF(BRGG(2).GT.0d0)
     .  WRITE(18,905) BRGG(2),2,22,22,"BR(H_2 -> gamma gamma)"
      IF(BRZG(2).GT.0d0)
     .  WRITE(18,905) BRZG(2),2,23,22,"BR(H_2 -> Z gamma)"
      IF(BRHHH(2,1).GT.0d0)
     .  WRITE(18,905) BRHHH(2,1),2,25,25,"BR(H_2 -> H_1 H_1)"
      IF(BRHAZ(2,1).GT.0d0)
     .  WRITE(18,905) BRHAZ(2,1),2,23,25,"BR(H_2 -> H_1 Z)"
      IF(BRHCHC(2).GT.0d0)
     .  WRITE(18,905) BRHCHC(2),2,37,-37,"BR(H_2 -> H+ H-)"
      IF(BRHCW(2).GT.0d0)
     .  WRITE(18,905) BRHCW(2),2,24,-37,"BR(H_2 -> W+ H-)"
      IF(BRHCW(2).GT.0d0)
     .  WRITE(18,905) BRHCW(2),2,-24,37,"BR(H_2 -> W- H+)"
      IF(BRNEU(2,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,1),2,1000022,1000022,
     .    "BR(H_2 -> neu_1 neu_1)"
      IF(BRNEU(2,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,2),2,1000022,1000023,
     .    "BR(H_2 -> neu_1 neu_2)"
      IF(BRNEU(2,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,3),2,1000022,1000025,
     .    "BR(H_2 -> neu_1 neu_3)"
      IF(BRNEU(2,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,4),2,1000022,1000035,
     .    "BR(H_2 -> neu_1 neu_4)"
      IF(BRNEU(2,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,5),2,1000022,1000045,
     .    "BR(H_2 -> neu_1 neu_5)"
      IF(BRNEU(2,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(2,2,2),2,1000023,1000023,
     .    "BR(H_2 -> neu_2 neu_2)"
      IF(BRNEU(2,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(2,2,3),2,1000023,1000025,
     .    "BR(H_2 -> neu_2 neu_3)"
      IF(BRNEU(2,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(2,2,4),2,1000023,1000035,
     .    "BR(H_2 -> neu_2 neu_4)"
      IF(BRNEU(2,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,2,5),2,1000023,1000045,
     .    "BR(H_2 -> neu_2 neu_5)"
      IF(BRNEU(2,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(2,3,3),2,1000025,1000025,
     .    "BR(H_2 -> neu_3 neu_3)"
      IF(BRNEU(2,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(2,3,4),2,1000025,1000035,
     .    "BR(H_2 -> neu_3 neu_4)"
      IF(BRNEU(2,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,3,5),2,1000025,1000045,
     .    "BR(H_2 -> neu_3 neu_5)"
      IF(BRNEU(2,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(2,4,4),2,1000035,1000035,
     .    "BR(H_2 -> neu_4 neu_4)"
      IF(BRNEU(2,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,4,5),2,1000035,1000045,
     .    "BR(H_2 -> neu_4 neu_5)"
      IF(BRNEU(2,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,5,5),2,1000045,1000045,
     .    "BR(H_2 -> neu_5 neu_5)"
      IF(BRCHA(2,1).GT.0d0)
     .  WRITE(18,905) BRCHA(2,1),2,1000024,-1000024,
     .    "BR(H_2 -> cha_1 cha_1bar)"
      IF(BRCHA(2,2).GT.0d0)
     .  WRITE(18,905) BRCHA(2,2),2,1000024,-1000037,
     .    "BR(H_2 -> cha_1 cha_2bar)"
      IF(BRCHA(2,2).GT.0d0)
     .  WRITE(18,905) BRCHA(2,2),2,1000037,-1000024,
     .    "BR(H_2 -> cha_2 cha_1bar)"
      IF(BRCHA(2,3).GT.0d0)
     .  WRITE(18,905) BRCHA(2,3),2,1000037,-1000037,
     .    "BR(H_2 -> cha_2 cha_2bar)"
      IF(BRHSQ(2,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,1),2,1000002,-1000002,
     .    "BR(H_2 -> ~u_L ~ubar_L)"
      IF(BRHSQ(2,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,1),2,1000004,-1000004,
     .    "BR(H_2 -> ~c_L ~cbar_L)"
      IF(BRHSQ(2,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,2),2,2000002,-2000002,
     .    "BR(H_2 -> ~u_R ~ubar_R)"
      IF(BRHSQ(2,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,2),2,2000004,-2000004,
     .    "BR(H_2 -> ~c_R ~cbar_R)"
      IF(BRHSQ(2,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,3),2,1000001,-1000001,
     .    "BR(H_2 -> ~d_L ~dbar_L)"
      IF(BRHSQ(2,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,3),2,1000003,-1000003,
     .    "BR(H_2 -> ~s_L ~sbar_L)"
      IF(BRHSQ(2,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,4),2,2000001,-2000001,
     .    "BR(H_2 -> ~d_R ~dbar_R)"
      IF(BRHSQ(2,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,4),2,2000003,-2000003,
     .    "BR(H_2 -> ~s_R ~sbar_R)"
      IF(BRHSQ(2,5).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,5),2,1000006,-1000006,
     .    "BR(H_2 -> ~t_1 ~tbar_1)"
      IF(BRHSQ(2,6).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,6),2,2000006,-2000006,
     .    "BR(H_2 -> ~t_2 ~tbar_2)"
      IF(BRHSQ(2,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,7),2,1000006,-2000006,
     .    "BR(H_2 -> ~t_1 ~tbar_2)"
      IF(BRHSQ(2,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,7),2,2000006,-1000006,
     .    "BR(H_2 -> ~t_2 ~tbar_1)"
      IF(BRHSQ(2,8).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,8),2,1000005,-1000005,
     .    "BR(H_2 -> ~b_1 ~bbar_1)"
      IF(BRHSQ(2,9).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,9),2,2000005,-2000005,
     .    "BR(H_2 -> ~b_2 ~bbar_2)"
      IF(BRHSQ(2,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,10),2,1000005,-2000005,
     .    "BR(H_2 -> ~b_1 ~bbar_2)"
      IF(BRHSQ(2,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,10),2,2000005,-1000005,
     .    "BR(H_2 -> ~b_2 ~bbar_1)"
      IF(BRHSL(2,1).GT.0d0)
     .  WRITE(18,905) BRHSL(2,1),2,1000011,-1000011,
     .    "BR(H_2 -> ~e_L ~ebar_L)"
      IF(BRHSL(2,1).GT.0d0)
     .  WRITE(18,905) BRHSL(2,1),2,1000013,-1000013,
     .    "BR(H_2 -> ~mu_L ~mubar_L)"
      IF(BRHSL(2,2).GT.0d0)
     .  WRITE(18,905) BRHSL(2,2),2,2000011,-2000011,
     .    "BR(H_2 -> ~e_R ~ebar_R)"
      IF(BRHSL(2,2).GT.0d0)
     .  WRITE(18,905) BRHSL(2,2),2,2000013,-2000013,
     .    "BR(H_2 -> ~mu_R ~mubarRL)"
      IF(BRHSL(2,3).GT.0d0)
     .  WRITE(18,905) BRHSL(2,3),2,1000012,-1000012,
     .    "BR(H_2 -> ~nu_e_L ~nu_ebar_L)"
      IF(BRHSL(2,3).GT.0d0)
     .  WRITE(18,905) BRHSL(2,3),2,1000014,-1000014,
     .    "BR(H_2 -> ~nu_mu_L ~nu_mubar_L)"
      IF(BRHSL(2,4).GT.0d0)
     .  WRITE(18,905) BRHSL(2,4),2,1000015,-1000015,
     .    "BR(H_2 -> ~tau_1 ~taubar_1)"
      IF(BRHSL(2,5).GT.0d0)
     .  WRITE(18,905) BRHSL(2,5),2,2000015,-2000015,
     .    "BR(H_2 -> ~tau_2 ~taubar_2)"
      IF(BRHSL(2,6).GT.0d0)
     .  WRITE(18,905) BRHSL(2,6),2,1000015,-2000015,
     .    "BR(H_2 -> ~tau_1 ~taubar_2)"
      IF(BRHSL(2,6).GT.0d0)
     .  WRITE(18,905) BRHSL(2,6),2,2000015,-1000015,
     .    "BR(H_2 -> ~tau_2 ~taubar_1)"
      IF(BRHSL(2,7).GT.0d0)
     .  WRITE(18,905) BRHSL(2,7),2,1000016,-1000016,
     .    "BR(H_2 -> ~nu_tau_L ~nu_taubar_L)"

      WRITE(18,904) 45,WIDTH(3),"3rd neutral Higgs scalar"
      IF(BRJJ(3).GT.0d0)
     .  WRITE(18,905) BRJJ(3),2,21,21,"BR(H_3 -> hadrons)"
      IF(BRMM(3).GT.0d0)
     .  WRITE(18,905) BRMM(3),2,13,-13,"BR(H_3 -> muon muon)"
      IF(BRLL(3).GT.0d0)
     .  WRITE(18,905) BRLL(3),2,15,-15,"BR(H_3 -> tau tau)"
      IF(BRCC(3).GT.0d0)
     .  WRITE(18,905) BRCC(3),2,4,-4,"BR(H_3 -> c cbar)"
      IF(BRBB(3).GT.0d0)
     .  WRITE(18,905) BRBB(3),2,5,-5,"BR(H_3 -> b bbar)"
      IF(BRTT(3).GT.0d0)
     .  WRITE(18,905) BRTT(3),2,6,-6,"BR(H_3 -> t tbar)"
      IF(BRWW(3).GT.0d0)
     .  WRITE(18,905) BRWW(3),2,24,-24,"BR(H_3 -> W+ W-)"
      IF(BRZZ(3).GT.0d0)
     .  WRITE(18,905) BRZZ(3),2,23,23,"BR(H_3 -> Z Z)"
      IF(BRGG(3).GT.0d0)
     .  WRITE(18,905) BRGG(3),2,22,22,"BR(H_3 -> gamma gamma)"
      IF(BRZG(3).GT.0d0)
     .  WRITE(18,905) BRZG(3),2,23,22,"BR(H_3 -> Z gamma)"
      IF(BRHHH(3,1).GT.0d0)
     .  WRITE(18,905) BRHHH(3,1),2,25,25,"BR(H_3 -> H_1 H_1)"
      IF(BRHHH(3,2).GT.0d0)
     .  WRITE(18,905) BRHHH(3,2),2,25,35,"BR(H_3 -> H_1 H_2)"
      IF(BRHHH(3,3).GT.0d0)
     .  WRITE(18,905) BRHHH(3,3),2,35,35,"BR(H_3 -> H_2 H_2)"
      IF(BRHAZ(3,1).GT.0d0)
     .  WRITE(18,905) BRHAZ(3,1),2,23,25,"BR(H_3 -> H_1 Z)"
      IF(BRHAZ(3,2).GT.0d0)
     .  WRITE(18,905) BRHAZ(3,2),2,23,35,"BR(H_3 -> H_2 Z)"
      IF(BRHCHC(3).GT.0d0)
     .  WRITE(18,905) BRHCHC(3),2,37,-37,"BR(H_3 -> H+ H-)"
      IF(BRHCW(3).GT.0d0)
     .  WRITE(18,905) BRHCW(3),2,24,-37,"BR(H_3 -> W+ H-)"
      IF(BRHCW(3).GT.0d0)
     .  WRITE(18,905) BRHCW(3),2,-24,37,"BR(H_3 -> W- H+)"
      IF(BRNEU(3,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,1),2,1000022,1000022,
     .    "BR(H_3 -> neu_1 neu_1)"
      IF(BRNEU(3,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,2),2,1000022,1000023,
     .    "BR(H_3 -> neu_1 neu_2)"
      IF(BRNEU(3,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,3),2,1000022,1000025,
     .    "BR(H_3 -> neu_1 neu_3)"
      IF(BRNEU(3,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,4),2,1000022,1000035,
     .    "BR(H_3 -> neu_1 neu_4)"
      IF(BRNEU(3,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,5),2,1000022,1000045,
     .    "BR(H_3 -> neu_1 neu_5)"
      IF(BRNEU(3,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(3,2,2),2,1000023,1000023,
     .    "BR(H_3 -> neu_2 neu_2)"
      IF(BRNEU(3,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(3,2,3),2,1000023,1000025,
     .    "BR(H_3 -> neu_2 neu_3)"
      IF(BRNEU(3,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(3,2,4),2,1000023,1000035,
     .    "BR(H_3 -> neu_2 neu_4)"
      IF(BRNEU(3,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,2,5),2,1000023,1000045,
     .    "BR(H_3 -> neu_2 neu_5)"
      IF(BRNEU(3,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(3,3,3),2,1000025,1000025,
     .    "BR(H_3 -> neu_3 neu_3)"
      IF(BRNEU(3,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(3,3,4),2,1000025,1000035,
     .    "BR(H_3 -> neu_3 neu_4)"
      IF(BRNEU(3,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,3,5),2,1000025,1000045,
     .    "BR(H_3 -> neu_3 neu_5)"
      IF(BRNEU(3,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(3,4,4),2,1000035,1000035,
     .    "BR(H_3 -> neu_4 neu_4)"
      IF(BRNEU(3,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,4,5),2,1000035,1000045,
     .    "BR(H_3 -> neu_4 neu_5)"
      IF(BRNEU(3,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,5,5),2,1000045,1000045,
     .    "BR(H_3 -> neu_5 neu_5)"
      IF(BRCHA(3,1).GT.0d0)
     .  WRITE(18,905) BRCHA(3,1),2,1000024,-1000024,
     .    "BR(H_3 -> cha_1 cha_1bar)"
      IF(BRCHA(3,2).GT.0d0)
     .  WRITE(18,905) BRCHA(3,2),2,1000024,-1000037,
     .    "BR(H_3 -> cha_1 cha_2bar)"
      IF(BRCHA(3,2).GT.0d0)
     .  WRITE(18,905) BRCHA(3,2),2,1000037,-1000024,
     .    "BR(H_3 -> cha_2 cha_1bar)"
      IF(BRCHA(3,3).GT.0d0)
     .  WRITE(18,905) BRCHA(3,3),2,1000037,-1000037,
     .    "BR(H_3 -> cha_2 cha_2bar)"
      IF(BRHSQ(3,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,1),2,1000002,-1000002,
     .    "BR(H_3 -> ~u_L ~ubar_L)"
      IF(BRHSQ(3,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,1),2,1000004,-1000004,
     .    "BR(H_3 -> ~c_L ~cbar_L)"
      IF(BRHSQ(3,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,2),2,2000002,-2000002,
     .    "BR(H_3 -> ~u_R ~ubar_R)"
      IF(BRHSQ(3,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,2),2,2000004,-2000004,
     .    "BR(H_3 -> ~c_R ~cbar_R)"
      IF(BRHSQ(3,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,3),2,1000001,-1000001,
     .    "BR(H_3 -> ~d_L ~dbar_L)"
      IF(BRHSQ(3,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,3),2,1000003,-1000003,
     .    "BR(H_3 -> ~s_L ~sbar_L)"
      IF(BRHSQ(3,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,4),2,2000001,-2000001,
     .    "BR(H_3 -> ~d_R ~dbar_R)"
      IF(BRHSQ(3,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,4),2,2000003,-2000003,
     .    "BR(H_3 -> ~s_R ~sbar_R)"
      IF(BRHSQ(3,5).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,5),2,1000006,-1000006,
     .    "BR(H_3 -> ~t_1 ~tbar_1)"
      IF(BRHSQ(3,6).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,6),2,2000006,-2000006,
     .    "BR(H_3 -> ~t_2 ~tbar_2)"
      IF(BRHSQ(3,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,7),2,1000006,-2000006,
     .    "BR(H_3 -> ~t_1 ~tbar_2)"
      IF(BRHSQ(3,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,7),2,2000006,-1000006,
     .    "BR(H_3 -> ~t_2 ~tbar_1)"
      IF(BRHSQ(3,8).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,8),2,1000005,-1000005,
     .    "BR(H_3 -> ~b_1 ~bbar_1)"
      IF(BRHSQ(3,9).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,9),2,2000005,-2000005,
     .    "BR(H_3 -> ~b_2 ~bbar_2)"
      IF(BRHSQ(3,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,10),2,1000005,-2000005,
     .    "BR(H_3 -> ~b_1 ~bbar_2)"
      IF(BRHSQ(3,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,10),2,2000005,-1000005,
     .    "BR(H_3 -> ~b_2 ~bbar_1)"
      IF(BRHSL(3,1).GT.0d0)
     .  WRITE(18,905) BRHSL(3,1),2,1000011,-1000011,
     .    "BR(H_3 -> ~e_L ~ebar_L)"
      IF(BRHSL(3,1).GT.0d0)
     .  WRITE(18,905) BRHSL(3,1),2,1000013,-1000013,
     .    "BR(H_3 -> ~mu_L ~mubar_L)"
      IF(BRHSL(3,2).GT.0d0)
     .  WRITE(18,905) BRHSL(3,2),2,2000011,-2000011,
     .    "BR(H_3 -> ~e_R ~ebar_R)"
      IF(BRHSL(3,2).GT.0d0)
     .  WRITE(18,905) BRHSL(3,2),2,2000013,-2000013,
     .    "BR(H_3 -> ~mu_R ~mubarRL)"
      IF(BRHSL(3,3).GT.0d0)
     .  WRITE(18,905) BRHSL(3,3),2,1000012,-1000012,
     .    "BR(H_3 -> ~nu_e_L ~nu_ebar_L)"
      IF(BRHSL(3,3).GT.0d0)
     .  WRITE(18,905) BRHSL(3,3),2,1000014,-1000014,
     .    "BR(H_3 -> ~nu_mu_L ~nu_mubar_L)"
      IF(BRHSL(3,4).GT.0d0)
     .  WRITE(18,905) BRHSL(3,4),2,1000015,-1000015,
     .    "BR(H_3 -> ~tau_1 ~taubar_1)"
      IF(BRHSL(3,5).GT.0d0)
     .  WRITE(18,905) BRHSL(3,5),2,2000015,-2000015,
     .    "BR(H_3 -> ~tau_2 ~taubar_2)"
      IF(BRHSL(3,6).GT.0d0)
     .  WRITE(18,905) BRHSL(3,6),2,1000015,-2000015,
     .    "BR(H_3 -> ~tau_1 ~taubar_2)"
      IF(BRHSL(3,6).GT.0d0)
     .  WRITE(18,905) BRHSL(3,6),2,2000015,-1000015,
     .    "BR(H_3 -> ~tau_2 ~taubar_1)"
      IF(BRHSL(3,7).GT.0d0)
     .  WRITE(18,905) BRHSL(3,7),2,1000016,-1000016,
     .    "BR(H_3 -> ~nu_tau_L ~nu_taubar_L)"

      WRITE(18,904) 36,WIDTH(4),"4th neutral Higgs scalar"
      IF(BRJJ(4).GT.0d0)
     .  WRITE(18,905) BRJJ(4),2,21,21,"BR(H_4 -> hadrons)"
      IF(BRMM(4).GT.0d0)
     .  WRITE(18,905) BRMM(4),2,13,-13,"BR(H_4 -> muon muon)"
      IF(BRLL(4).GT.0d0)
     .  WRITE(18,905) BRLL(4),2,15,-15,"BR(H_4 -> tau tau)"
      IF(BRCC(4).GT.0d0)
     .  WRITE(18,905) BRCC(4),2,4,-4,"BR(H_4 -> c cbar)"
      IF(BRBB(4).GT.0d0)
     .  WRITE(18,905) BRBB(4),2,5,-5,"BR(H_4 -> b bbar)"
      IF(BRTT(4).GT.0d0)
     .  WRITE(18,905) BRTT(4),2,6,-6,"BR(H_4 -> t tbar)"
      IF(BRGG(4).GT.0d0)
     .  WRITE(18,905) BRGG(4),2,22,22,"BR(H_4 -> gamma gamma)"
      IF(BRZG(4).GT.0d0)
     .  WRITE(18,905) BRZG(4),2,23,22,"BR(H_4 -> Z gamma)"
      IF(BRHHH(4,1).GT.0d0)
     .  WRITE(18,905) BRHHH(4,1),2,25,25,"BR(H_4 -> H_1 H_1)"
      IF(BRHHH(4,2).GT.0d0)
     .  WRITE(18,905) BRHHH(4,2),2,25,35,"BR(H_4 -> H_1 H_2)"
      IF(BRHHH(4,3).GT.0d0)
     .  WRITE(18,905) BRHHH(4,3),2,35,35,"BR(H_4 -> H_2 H_2)"
      IF(BRHHH(4,4).GT.0d0)
     .  WRITE(18,905) BRHHH(4,4),2,25,45,"BR(H_4 -> H_1 H_3)"
      IF(BRHHH(4,5).GT.0d0)
     .  WRITE(18,905) BRHHH(4,5),2,35,45,"BR(H_4 -> H_2 H_3)"
      IF(BRHHH(4,6).GT.0d0)
     .  WRITE(18,905) BRHHH(4,6),2,45,45,"BR(H_4 -> H_3 H_3)"
      IF(BRHAZ(4,1).GT.0d0)
     .  WRITE(18,905) BRHAZ(4,1),2,23,25,"BR(H_4 -> H_1 Z)"
      IF(BRHAZ(4,2).GT.0d0)
     .  WRITE(18,905) BRHAZ(4,2),2,23,35,"BR(H_4 -> H_2 Z)"
      IF(BRHAZ(4,3).GT.0d0)
     .  WRITE(18,905) BRHAZ(4,3),2,23,45,"BR(H_4 -> H_3 Z)"
      IF(BRHCHC(4).GT.0d0)
     .  WRITE(18,905) BRHCHC(4),2,37,-37,"BR(H_4 -> H+ H-)"
      IF(BRHCW(4).GT.0d0)
     .  WRITE(18,905) BRHCW(4),2,24,-37,"BR(H_4 -> W+ H-)"
      IF(BRHCW(4).GT.0d0)
     .  WRITE(18,905) BRHCW(4),2,-24,37,"BR(H_4 -> W- H+)"
      IF(BRNEU(4,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,1),2,1000022,1000022,
     .    "BR(H_4 -> neu_1 neu_1)"
      IF(BRNEU(4,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,2),2,1000022,1000023,
     .    "BR(H_4 -> neu_1 neu_2)"
      IF(BRNEU(4,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,3),2,1000022,1000025,
     .    "BR(H_4 -> neu_1 neu_3)"
      IF(BRNEU(4,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,4),2,1000022,1000035,
     .    "BR(H_4 -> neu_1 neu_4)"
      IF(BRNEU(4,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,5),2,1000022,1000045,
     .    "BR(H_4 -> neu_1 neu_5)"
      IF(BRNEU(4,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(4,2,2),2,1000023,1000023,
     .    "BR(H_4 -> neu_2 neu_2)"
      IF(BRNEU(4,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(4,2,3),2,1000023,1000025,
     .    "BR(H_4 -> neu_2 neu_3)"
      IF(BRNEU(4,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(4,2,4),2,1000023,1000035,
     .    "BR(H_4 -> neu_2 neu_4)"
      IF(BRNEU(4,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,2,5),2,1000023,1000045,
     .    "BR(H_4 -> neu_2 neu_5)"
      IF(BRNEU(4,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(4,3,3),2,1000025,1000025,
     .    "BR(H_4 -> neu_3 neu_3)"
      IF(BRNEU(4,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(4,3,4),2,1000025,1000035,
     .    "BR(H_4 -> neu_3 neu_4)"
      IF(BRNEU(4,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,3,5),2,1000025,1000045,
     .    "BR(H_4 -> neu_3 neu_5)"
      IF(BRNEU(4,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(4,4,4),2,1000035,1000035,
     .    "BR(H_4 -> neu_4 neu_4)"
      IF(BRNEU(4,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,4,5),2,1000035,1000045,
     .    "BR(H_4 -> neu_4 neu_5)"
      IF(BRNEU(4,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,5,5),2,1000045,1000045,
     .    "BR(H_4 -> neu_5 neu_5)"
      IF(BRCHA(4,1).GT.0d0)
     .  WRITE(18,905) BRCHA(4,1),2,1000024,-1000024,
     .    "BR(H_4 -> cha_1 cha_1bar)"
      IF(BRCHA(4,2).GT.0d0)
     .  WRITE(18,905) BRCHA(4,2),2,1000024,-1000037,
     .    "BR(H_4 -> cha_1 cha_2bar)"
      IF(BRCHA(4,2).GT.0d0)
     .  WRITE(18,905) BRCHA(4,2),2,1000037,-1000024,
     .    "BR(H_4 -> cha_2 cha_1bar)"
      IF(BRCHA(4,3).GT.0d0)
     .  WRITE(18,905) BRCHA(4,3),2,1000037,-1000037,
     .    "BR(H_4 -> cha_2 cha_2bar)"
      IF(BRHSQ(4,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,1),2,1000002,-1000002,
     .    "BR(H_4 -> ~u_L ~ubar_L)"
      IF(BRHSQ(4,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,1),2,1000004,-1000004,
     .    "BR(H_4 -> ~c_L ~cbar_L)"
      IF(BRHSQ(4,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,2),2,2000002,-2000002,
     .    "BR(H_4 -> ~u_R ~ubar_R)"
      IF(BRHSQ(4,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,2),2,2000004,-2000004,
     .    "BR(H_4 -> ~c_R ~cbar_R)"
      IF(BRHSQ(4,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,3),2,1000001,-1000001,
     .    "BR(H_4 -> ~d_L ~dbar_L)"
      IF(BRHSQ(4,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,3),2,1000003,-1000003,
     .    "BR(H_4 -> ~s_L ~sbar_L)"
      IF(BRHSQ(4,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,4),2,2000001,-2000001,
     .    "BR(H_4 -> ~d_R ~dbar_R)"
      IF(BRHSQ(4,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,4),2,2000003,-2000003,
     .    "BR(H_4 -> ~s_R ~sbar_R)"
      IF(BRHSQ(4,5).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,5),2,1000006,-1000006,
     .    "BR(H_4 -> ~t_1 ~tbar_1)"
      IF(BRHSQ(4,6).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,6),2,2000006,-2000006,
     .    "BR(H_4 -> ~t_2 ~tbar_2)"
      IF(BRHSQ(4,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,7),2,1000006,-2000006,
     .    "BR(H_4 -> ~t_1 ~tbar_2)"
      IF(BRHSQ(4,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,7),2,2000006,-1000006,
     .    "BR(H_4 -> ~t_2 ~tbar_1)"
      IF(BRHSQ(4,8).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,8),2,1000005,-1000005,
     .    "BR(H_4 -> ~b_1 ~bbar_1)"
      IF(BRHSQ(4,9).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,9),2,2000005,-2000005,
     .    "BR(H_4 -> ~b_2 ~bbar_2)"
      IF(BRHSQ(4,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,10),2,1000005,-2000005,
     .    "BR(H_4 -> ~b_1 ~bbar_2)"
      IF(BRHSQ(4,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(4,10),2,2000005,-1000005,
     .    "BR(H_4 -> ~b_2 ~bbar_1)"
      IF(BRHSL(4,1).GT.0d0)
     .  WRITE(18,905) BRHSL(4,1),2,1000011,-1000011,
     .    "BR(H_4 -> ~e_L ~ebar_L)"
      IF(BRHSL(4,1).GT.0d0)
     .  WRITE(18,905) BRHSL(4,1),2,1000013,-1000013,
     .    "BR(H_4 -> ~mu_L ~mubar_L)"
      IF(BRHSL(4,2).GT.0d0)
     .  WRITE(18,905) BRHSL(4,2),2,2000011,-2000011,
     .    "BR(H_4 -> ~e_R ~ebar_R)"
      IF(BRHSL(4,2).GT.0d0)
     .  WRITE(18,905) BRHSL(4,2),2,2000013,-2000013,
     .    "BR(H_4 -> ~mu_R ~mubarRL)"
      IF(BRHSL(4,3).GT.0d0)
     .  WRITE(18,905) BRHSL(4,3),2,1000012,-1000012,
     .    "BR(H_4 -> ~nu_e_L ~nu_ebar_L)"
      IF(BRHSL(4,3).GT.0d0)
     .  WRITE(18,905) BRHSL(4,3),2,1000014,-1000014,
     .    "BR(H_4 -> ~nu_mu_L ~nu_mubar_L)"
      IF(BRHSL(4,4).GT.0d0)
     .  WRITE(18,905) BRHSL(4,4),2,1000015,-1000015,
     .    "BR(H_4 -> ~tau_1 ~taubar_1)"
      IF(BRHSL(4,5).GT.0d0)
     .  WRITE(18,905) BRHSL(4,5),2,2000015,-2000015,
     .    "BR(H_4 -> ~tau_2 ~taubar_2)"
      IF(BRHSL(4,6).GT.0d0)
     .  WRITE(18,905) BRHSL(4,6),2,1000015,-2000015,
     .    "BR(H_4 -> ~tau_1 ~taubar_2)"
      IF(BRHSL(4,6).GT.0d0)
     .  WRITE(18,905) BRHSL(4,6),2,2000015,-1000015,
     .    "BR(H_4 -> ~tau_2 ~taubar_1)"
      IF(BRHSL(4,7).GT.0d0)
     .  WRITE(18,905) BRHSL(4,7),2,1000016,-1000016,
     .    "BR(H_4 -> ~nu_tau_L ~nu_taubar_L)"

      WRITE(18,904) 46,WIDTH(5),"5th neutral Higgs scalar"
      IF(BRJJ(5).GT.0d0)
     .  WRITE(18,905) BRJJ(5),2,21,21,"BR(H_5 -> hadrons)"
      IF(BRMM(5).GT.0d0)
     .  WRITE(18,905) BRMM(5),2,13,-13,"BR(H_5 -> muon muon)"
      IF(BRLL(5).GT.0d0)
     .  WRITE(18,905) BRLL(5),2,15,-15,"BR(H_5 -> tau tau)"
      IF(BRCC(5).GT.0d0)
     .  WRITE(18,905) BRCC(5),2,4,-4,"BR(H_5 -> c cbar)"
      IF(BRBB(5).GT.0d0)
     .  WRITE(18,905) BRBB(5),2,5,-5,"BR(H_5 -> b bbar)"
      IF(BRTT(5).GT.0d0)
     .  WRITE(18,905) BRTT(5),2,6,-6,"BR(H_5 -> t tbar)"
      IF(BRGG(5).GT.0d0)
     .  WRITE(18,905) BRGG(5),2,22,22,"BR(H_5 -> gamma gamma)"
      IF(BRZG(5).GT.0d0)
     .  WRITE(18,905) BRZG(5),2,23,22,"BR(H_5 -> Z gamma)"
      IF(BRHHH(5,1).GT.0d0)
     .  WRITE(18,905) BRHHH(5,1),2,25,25,"BR(H_5 -> H_1 H_1)"
      IF(BRHHH(5,2).GT.0d0)
     .  WRITE(18,905) BRHHH(5,2),2,25,35,"BR(H_5 -> H_1 H_2)"
      IF(BRHHH(5,3).GT.0d0)
     .  WRITE(18,905) BRHHH(5,3),2,35,35,"BR(H_5 -> H_2 H_2)"
      IF(BRHHH(5,4).GT.0d0)
     .  WRITE(18,905) BRHHH(5,4),2,25,45,"BR(H_5 -> H_1 H_3)"
      IF(BRHHH(5,5).GT.0d0)
     .  WRITE(18,905) BRHHH(5,5),2,35,45,"BR(H_5 -> H_2 H_3)"
      IF(BRHHH(5,6).GT.0d0)
     .  WRITE(18,905) BRHHH(5,6),2,45,45,"BR(H_5 -> H_3 H_3)"
      IF(BRHHH(5,7).GT.0d0)
     .  WRITE(18,905) BRHHH(5,7),2,25,36,"BR(H_5 -> H_1 H_4)"
      IF(BRHHH(5,8).GT.0d0)
     .  WRITE(18,905) BRHHH(5,8),2,35,36,"BR(H_5 -> H_2 H_4)"
      IF(BRHHH(5,9).GT.0d0)
     .  WRITE(18,905) BRHHH(5,9),2,45,36,"BR(H_5 -> H_3 H_4)"
      IF(BRHHH(5,10).GT.0d0)
     .  WRITE(18,905) BRHHH(5,10),2,36,36,"BR(H_5 -> H_4 H_4)"
      IF(BRHCHC(5).GT.0d0)
     .  WRITE(18,905) BRHCHC(5),2,37,-37,"BR(H_5 -> H+ H-)"
      IF(BRHAZ(5,1).GT.0d0)
     .  WRITE(18,905) BRHAZ(5,1),2,23,25,"BR(H_5 -> H_1 Z)"
      IF(BRHAZ(5,2).GT.0d0)
     .  WRITE(18,905) BRHAZ(5,2),2,23,35,"BR(H_5 -> H_2 Z)"
      IF(BRHAZ(5,3).GT.0d0)
     .  WRITE(18,905) BRHAZ(5,3),2,23,45,"BR(H_5 -> H_3 Z)"
      IF(BRHAZ(5,4).GT.0d0)
     .  WRITE(18,905) BRHAZ(5,4),2,23,36,"BR(H_5 -> H_4 Z)"
      IF(BRHCW(5).GT.0d0)
     .  WRITE(18,905) BRHCW(5),2,24,-37,"BR(H_5 -> W+ H-)"
      IF(BRHCW(5).GT.0d0)
     .  WRITE(18,905) BRHCW(5),2,-24,37,"BR(H_5 -> W- H+)"
      IF(BRNEU(5,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,1),2,1000022,1000022,
     .    "BR(H_5 -> neu_1 neu_1)"
      IF(BRNEU(5,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,2),2,1000022,1000023,
     .    "BR(H_5 -> neu_1 neu_2)"
      IF(BRNEU(5,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,3),2,1000022,1000025,
     .    "BR(H_5 -> neu_1 neu_3)"
      IF(BRNEU(5,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,4),2,1000022,1000035,
     .    "BR(H_5 -> neu_1 neu_4)"
      IF(BRNEU(5,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,5),2,1000022,1000045,
     .    "BR(H_5 -> neu_1 neu_5)"
      IF(BRNEU(5,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(5,2,2),2,1000023,1000023,
     .    "BR(H_5 -> neu_2 neu_2)"
      IF(BRNEU(5,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(5,2,3),2,1000023,1000025,
     .    "BR(H_5 -> neu_2 neu_3)"
      IF(BRNEU(5,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(5,2,4),2,1000023,1000035,
     .    "BR(H_5 -> neu_2 neu_4)"
      IF(BRNEU(5,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,2,5),2,1000023,1000045,
     .    "BR(H_5 -> neu_2 neu_5)"
      IF(BRNEU(5,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(5,3,3),2,1000025,1000025,
     .    "BR(H_5 -> neu_3 neu_3)"
      IF(BRNEU(5,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(5,3,4),2,1000025,1000035,
     .    "BR(H_5 -> neu_3 neu_4)"
      IF(BRNEU(5,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,3,5),2,1000025,1000045,
     .    "BR(H_5 -> neu_3 neu_5)"
      IF(BRNEU(5,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(5,4,4),2,1000035,1000035,
     .    "BR(H_5 -> neu_4 neu_4)"
      IF(BRNEU(5,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,4,5),2,1000035,1000045,
     .    "BR(H_5 -> neu_4 neu_5)"
      IF(BRNEU(5,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,5,5),2,1000045,1000045,
     .    "BR(H_5 -> neu_5 neu_5)"
      IF(BRCHA(5,1).GT.0d0)
     .  WRITE(18,905) BRCHA(5,1),2,1000024,-1000024,
     .    "BR(H_5 -> cha_1 cha_1bar)"
      IF(BRCHA(5,2).GT.0d0)
     .  WRITE(18,905) BRCHA(5,2),2,1000024,-1000037,
     .    "BR(H_5 -> cha_1 cha_2bar)"
      IF(BRCHA(5,2).GT.0d0)
     .  WRITE(18,905) BRCHA(5,2),2,1000037,-1000024,
     .    "BR(H_5 -> cha_2 cha_1bar)"
      IF(BRCHA(5,3).GT.0d0)
     .  WRITE(18,905) BRCHA(5,3),2,1000037,-1000037,
     .    "BR(H_5 -> cha_2 cha_2bar)"
      IF(BRHSQ(5,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,1),2,1000002,-1000002,
     .    "BR(H_5 -> ~u_L ~ubar_L)"
      IF(BRHSQ(5,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,1),2,1000004,-1000004,
     .    "BR(H_5 -> ~c_L ~cbar_L)"
      IF(BRHSQ(5,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,2),2,2000002,-2000002,
     .    "BR(H_5 -> ~u_R ~ubar_R)"
      IF(BRHSQ(5,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,2),2,2000004,-2000004,
     .    "BR(H_5 -> ~c_R ~cbar_R)"
      IF(BRHSQ(5,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,3),2,1000001,-1000001,
     .    "BR(H_5 -> ~d_L ~dbar_L)"
      IF(BRHSQ(5,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,3),2,1000003,-1000003,
     .    "BR(H_5 -> ~s_L ~sbar_L)"
      IF(BRHSQ(5,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,4),2,2000001,-2000001,
     .    "BR(H_5 -> ~d_R ~dbar_R)"
      IF(BRHSQ(5,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,4),2,2000003,-2000003,
     .    "BR(H_5 -> ~s_R ~sbar_R)"
      IF(BRHSQ(5,5).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,5),2,1000006,-1000006,
     .    "BR(H_5 -> ~t_1 ~tbar_1)"
      IF(BRHSQ(5,6).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,6),2,2000006,-2000006,
     .    "BR(H_5 -> ~t_2 ~tbar_2)"
      IF(BRHSQ(5,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,7),2,1000006,-2000006,
     .    "BR(H_5 -> ~t_1 ~tbar_2)"
      IF(BRHSQ(5,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,7),2,2000006,-1000006,
     .    "BR(H_5 -> ~t_2 ~tbar_1)"
      IF(BRHSQ(5,8).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,8),2,1000005,-1000005,
     .    "BR(H_5 -> ~b_1 ~bbar_1)"
      IF(BRHSQ(5,9).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,9),2,2000005,-2000005,
     .    "BR(H_5 -> ~b_2 ~bbar_2)"
      IF(BRHSQ(5,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,10),2,1000005,-2000005,
     .    "BR(H_5 -> ~b_1 ~bbar_2)"
      IF(BRHSQ(5,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(5,10),2,2000005,-1000005,
     .    "BR(H_5 -> ~b_2 ~bbar_1)"
      IF(BRHSL(5,1).GT.0d0)
     .  WRITE(18,905) BRHSL(5,1),2,1000011,-1000011,
     .    "BR(H_5 -> ~e_L ~ebar_L)"
      IF(BRHSL(5,1).GT.0d0)
     .  WRITE(18,905) BRHSL(5,1),2,1000013,-1000013,
     .    "BR(H_5 -> ~mu_L ~mubar_L)"
      IF(BRHSL(5,2).GT.0d0)
     .  WRITE(18,905) BRHSL(5,2),2,2000011,-2000011,
     .    "BR(H_5 -> ~e_R ~ebar_R)"
      IF(BRHSL(5,2).GT.0d0)
     .  WRITE(18,905) BRHSL(5,2),2,2000013,-2000013,
     .    "BR(H_5 -> ~mu_R ~mubarRL)"
      IF(BRHSL(5,3).GT.0d0)
     .  WRITE(18,905) BRHSL(5,3),2,1000012,-1000012,
     .    "BR(H_5 -> ~nu_e_L ~nu_ebar_L)"
      IF(BRHSL(5,3).GT.0d0)
     .  WRITE(18,905) BRHSL(5,3),2,1000014,-1000014,
     .    "BR(H_5 -> ~nu_mu_L ~nu_mubar_L)"
      IF(BRHSL(5,4).GT.0d0)
     .  WRITE(18,905) BRHSL(5,4),2,1000015,-1000015,
     .    "BR(H_5 -> ~tau_1 ~taubar_1)"
      IF(BRHSL(5,5).GT.0d0)
     .  WRITE(18,905) BRHSL(5,5),2,2000015,-2000015,
     .    "BR(H_5 -> ~tau_2 ~taubar_2)"
      IF(BRHSL(5,6).GT.0d0)
     .  WRITE(18,905) BRHSL(5,6),2,1000015,-2000015,
     .    "BR(H_5 -> ~tau_1 ~taubar_2)"
      IF(BRHSL(5,6).GT.0d0)
     .  WRITE(18,905) BRHSL(5,6),2,2000015,-1000015,
     .    "BR(H_5 -> ~tau_2 ~taubar_1)"
      IF(BRHSL(5,7).GT.0d0)
     .  WRITE(18,905) BRHSL(3,7),2,1000016,-1000016,
     .    "BR(H_5 -> ~nu_tau_L ~nu_taubar_L)"

      WRITE(18,904) 37,HCWIDTH,"Charged Higgs"
      IF(HCBRM.GT.0d0)
     .  WRITE(18,905) HCBRM,2,-13,14,"BR(H+ -> muon nu_muon)"
      IF(HCBRL.GT.0d0)
     .  WRITE(18,905) HCBRL,2,-15,16,"BR(H+ -> tau nu_tau)"
      IF(HCBRSU.GT.0d0)
     .  WRITE(18,905) HCBRSU,2,2,-3,"BR(H+ -> u sbar)"
      IF(HCBRSC.GT.0d0)
     .  WRITE(18,905) HCBRSC,2,4,-3,"BR(H+ -> c sbar)"
      IF(HCBRBU.GT.0d0)
     .  WRITE(18,905) HCBRBU,2,2,-5,"BR(H+ -> u bbar)"
      IF(HCBRBC.GT.0d0)
     .  WRITE(18,905) HCBRBC,2,4,-5,"BR(H+ -> c bbar)"
      IF(HCBRBT.GT.0d0)
     .  WRITE(18,905) HCBRBT,2,6,-5,"BR(H+ -> t bbar)"
      IF(HCBRWH(1).GT.0d0)
     .  WRITE(18,905) HCBRWH(1),2,24,25,"BR(H+ -> W+ H_1)"
      IF(HCBRWH(2).GT.0d0)
     .  WRITE(18,905) HCBRWH(2),2,24,35,"BR(H+ -> W+ H_2)"
      IF(HCBRWH(3).GT.0d0)
     .  WRITE(18,905) HCBRWH(3),2,24,45,"BR(H+ -> W+ H_3)"
      IF(HCBRWH(4).GT.0d0)
     .  WRITE(18,905) HCBRWH(4),2,24,36,"BR(H+ -> W+ H_4)"
      IF(HCBRWH(4).GT.0d0)
     .  WRITE(18,905) HCBRWH(5),2,24,46,"BR(H+ -> W+ H_5)"

      IF(HCBRNC(1,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(1,1),2,1000024,1000022,
     .    "BR(H+ -> cha_1 neu_1)"
      IF(HCBRNC(2,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(2,1),2,1000024,1000023,
     .    "BR(H+ -> cha_1 neu_2)"
      IF(HCBRNC(3,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(3,1),2,1000024,1000025,
     .    "BR(H+ -> cha_1 neu_3)"
      IF(HCBRNC(4,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(4,1),2,1000024,1000035,
     .    "BR(H+ -> cha_1 neu_4)"
      IF(HCBRNC(5,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(5,1),2,1000024,1000045,
     .    "BR(H+ -> cha_1 neu_5)"
      IF(HCBRNC(1,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(1,2),2,1000037,1000022,
     .    "BR(H+ -> cha_2 neu_1)"
      IF(HCBRNC(2,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(2,2),2,1000037,1000023,
     .    "BR(H+ -> cha_2 neu_2)"
      IF(HCBRNC(3,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(3,2),2,1000037,1000025,
     .    "BR(H+ -> cha_2 neu_3)"
      IF(HCBRNC(4,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(4,2),2,1000037,1000035,
     .    "BR(H+ -> cha_2 neu_4)"
      IF(HCBRNC(5,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(5,2),2,1000037,1000045,
     .    "BR(H+ -> cha_2 neu_5)"
      IF(HCBRSQ(1).GT.0d0)
     .  WRITE(18,905) HCBRSQ(1),2,1000002,-1000001,
     .    "BR(H+ -> ~u_L ~dbar_L)"
      IF(HCBRSQ(1).GT.0d0)
     .  WRITE(18,905) HCBRSQ(1),2,1000004,-1000003,
     .    "BR(H+ -> ~c_L ~sbar_L)"
      IF(HCBRSQ(2).GT.0d0)
     .  WRITE(18,905) HCBRSQ(2),2,1000006,-1000005,
     .    "BR(H+ -> ~t_1 ~bbar_1)"
      IF(HCBRSQ(3).GT.0d0)
     .  WRITE(18,905) HCBRSQ(3),2,1000006,-2000005,
     .    "BR(H+ -> ~t_1 ~bbar_2)"
      IF(HCBRSQ(4).GT.0d0)
     .  WRITE(18,905) HCBRSQ(4),2,2000006,-1000005,
     .    "BR(H+ -> ~t_2 ~bbar_1)"
      IF(HCBRSQ(5).GT.0d0)
     .  WRITE(18,905) HCBRSQ(5),2,2000006,-2000005,
     .    "BR(H+ -> ~t_2 ~bbar_2)"
      IF(HCBRSL(1).GT.0d0)
     .  WRITE(18,905) HCBRSL(1),2,1000012,-1000011,
     .    "BR(H+ -> ~nu_e_L ~ebar_L)"
      IF(HCBRSL(1).GT.0d0)
     .  WRITE(18,905) HCBRSL(1),2,1000014,-1000013,
     .    "BR(H+ -> ~nu_mu_L ~mubar_L)"
      IF(HCBRSL(2).GT.0d0)
     .  WRITE(18,905) HCBRSL(2),2,1000016,-1000015,
     .    "BR(H+ -> ~nu_tau_L ~taubar_1)"
      IF(HCBRSL(3).GT.0d0)
     .  WRITE(18,905) HCBRSL(3),2,1000016,-2000015,
     .    "BR(H+ -> ~nu_tau_L ~taubar_2)"

      WRITE(18,904) 6,toptot,'Top Quark'
      IF(brtopbw.ne.0.D0)
     .  WRITE(18,905) brtopbw,2,5,24,'BR(t ->  b    W+)'
      IF(brtopbh.ne.0.D0)
     .  WRITE(18,905) brtopbh,2,5,37,'BR(t ->  b    H+)'
      IF(brtopneutrstop(1,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(1,1),2,1000006,1000022,
     . 'BR(t -> ~t_1 ~chi_10)'
      IF(brtopneutrstop(2,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(2,1),2,1000006,1000023,
     . 'BR(t -> ~t_1 ~chi_20)'
      IF(brtopneutrstop(3,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(3,1),2,1000006,1000025,
     . 'BR(t -> ~t_1 ~chi_30)'
      IF(brtopneutrstop(4,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(4,1),2,1000006,1000025,
     . 'BR(t -> ~t_1 ~chi_40)'
      IF(brtopneutrstop(5,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(5,1),2,1000006,1000025,
     . 'BR(t -> ~t_1 ~chi_50)'
      IF(brtopneutrstop(1,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(1,2),2,2000006,1000022,
     . 'BR(t -> ~t_2 ~chi_10)'
      IF(brtopneutrstop(2,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(2,2),2,2000006,1000023,
     . 'BR(t -> ~t_2 ~chi_20)'
      IF(brtopneutrstop(3,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(3,2),2,2000006,1000025,
     .'BR(t -> ~t_2 ~chi_30)'
      IF(brtopneutrstop(4,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(4,2),2,2000006,1000025,
     . 'BR(t -> ~t_2 ~chi_40)'
      IF(brtopneutrstop(5,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(5,2),2,2000006,1000025,
     . 'BR(t -> ~t_2 ~chi_50)'

 899  FORMAT(A)
 900  FORMAT(1X,I5,3X,A)
 901  FORMAT(1X,I5,3X,1P,E16.8,0P,3X,'#',1X,A)
 902  FORMAT(1X,I9,3X,1P,E16.8,0P,3X,'#',1X,A)
 903  FORMAT(1X,I2,1X,I2,3X,1P,E16.8,0P,3X,'#',1X,A)
 904  FORMAT('DECAY',1X,I9,3X,1P,E16.8,0P,3X,'#',1X,A)
 905  FORMAT(3X,1P,E16.8,0P,3X,I2,3X,I9,1X,I9,1X,2X,'#',1X,A)
 906  FORMAT('#',1X,A,3X,E16.8)
 907  FORMAT(A,1P,E16.8,A)
 918  FORMAT(1X,I5,3X,A,F6.1,'-',F5.1,A)
 920  FORMAT('#',0P,I5,3X,1P,E16.8,0P,3X,'#',1X,A)
 921  FORMAT(1X,I2,1X,I2,3X,'#',1X,A)

      END
