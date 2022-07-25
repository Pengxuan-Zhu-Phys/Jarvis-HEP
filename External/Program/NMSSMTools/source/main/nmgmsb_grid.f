      PROGRAM MAIN

*  Program to compute the NMSSM Higgs masses, couplings and
*  Branching ratios, with experimental and theoretical constraints
*
*  Input data corresponds to GMSB model. In addition:
*
*      tan(beta) at the scale MZ, lambda at the scale Q2
*
*      The input file contains lower and upper bounds as well as number
*      of steps for the parameters on which the scan is performed
*
*  On output:
*
*      PAR(1) = lambda
*      PAR(2) = kappa
*      PAR(3) = tan(beta)
*      PAR(4) = mu (effective mu term = lambda*s)
*      PAR(5) = Alambda
*      PAR(6) = Akappa
*      PAR(7) = mQ3**2
*      PAR(8) = mU3**2
*      PAR(9) = mD3**2
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
*      PAR(20) = M1
*      PAR(21) = M2
*      PAR(22) = M3
*      PAR(23) = MA (diagonal doublet CP-odd mass matrix element)
*      PAR(24) = MP (diagonal singlet CP-odd mass matrix element)
*      PAR(25) = AE2
*
*      All these parameters are assumed to be defined in DRbar at the scale
*      Q2 which is either user defined or computed as (2*mQ2+mU2+mD2)/4
*
*      SMASS(1-3): CP-even masses (ordered)
*
*      SCOMP(1-3,1-3): Mixing angles: if HB(I) are the bare states,
*        HB(I) = Re(H1), Re(H2), Re(S), and HM(I) are the mass eigenstates,
*        the convention is HB(I) = SUM_(J=1,3) SCOMP(J,I)*HM(J)
*        which is equivalent to HM(I) = SUM_(J=1,3) SCOMP(I,J)*HB(J)
*
*      PMASS(1-2): CP-odd masses (ordered)
*
*      PCOMP(1-2,1-2): Mixing angles: if AB(I) are the bare states,
*        AB(I) = Im(H1), Im(H2), Im(S), and AM(I) are the mass eigenstates,
*        the convention is
*        AM(I) = PCOMP(I,1)*(COSBETA*AB(1)+SINBETA*AB(2))
*              + PCOMP(I,2)*AB(3)
*
*      CMASS: Charged Higgs mass
*
*      CU,CD,CV,CJ,CG(i) Reduced couplings of h1,h2,h3 (i=1,2,3) or
*                        a1,a2 (i=4,5) to up type fermions, down type
*                        fermions, gauge bosons, gluons and photons
*                        Note: CV(4)=CV(5)=0
*      CB(I)             Reduced couplings of h1,h2,h3 (i=1,2,3) or
*                        a1,a2 (i=4,5) to b-quarks including DELMB corrections
*
*      WIDTH(i) Total decay width of h1,h2,h3,a1,a2 (i=1..5)
*               with the following branching ratios:
*      BRJJ(i) h1,h2,h3,a1,a2 -> hadrons
*      BRMM(i)        "       -> mu mu
*      BRLL(i)        "       -> tau tau
*      BRCC(i)        "       -> cc
*      BRBB(i)        "       -> bb
*      BRTT(i)        "       -> tt
*      BRWW(i)        "       -> WW (BRWW(4)=BRWW(5)=0)
*      BRZZ(i)        "       -> ZZ (BRZZ(4)=BRZZ(5)=0)
*      BRGG(i)        "       -> gamma gamma
*      BRZG(i)        "       -> Z gamma
*      BRHIGGS(i)   (i=1..5)  -> other Higgses, including:
*        BRHAA(i,j)   hi      -> a1a1, a1a2, a2a2 (i=1..3, j=1..3)
*        BRHCHC(i)    hi      -> h+h- (i=1..3)
*        BRHAZ(i,j)   hi      -> Zaj  (i=1..3, j=1..2)
*        BRHCW(i)  h1,h2,h3   -> h+W- (i=1..3), a1,a2 -> h+W- (i=4,5)
*        BRHHH(i)     h2      -> h1h1, h3-> h1h1, h1h2, h2h2 (i=1..4)
*        BRAHA(i)     a2      -> a1hi (i=1..3)
*        BRAHZ(i,j)   ai      -> Zhj  (i=1,2, j=1..3)
*      BRSUSY(i)    (i=1..5)  -> susy particles, including:
*        BRNEU(i,j,k)         -> neutralinos j,k (i=1..5, j,k=1..5)
*        BRCHA(i,j)           -> charginos 11, 12, 22 (i=1..5, j=1..3)
*        BRHSQ(i,j)   hi      -> uLuL, uRuR, dLdL, dRdR, t1t1, t2t2,
*                                t1t2, b1b1, b2b2, b1b2 (i=1..3, j=1..10)
*        BRASQ(i,j)   ai      -> t1t2, b1b2 (i=1,2, j=1,2)
*        BRHSL(i,j)   hi      -> lLlL, lRlR, nLnL, l1l1, l2l2, l1l2,
*                                ntnt (i=1..3, j=1..7)
*        BRASL(i)     ai      -> l1l2 (i=1,2)
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
*        HCBRWH(i)   "  -> H1W+, H2W+, h3W+, a1W+, a2W+ (i=1..5)
*      HCBRSUSY      "  -> susy particles,including
*        HCBRNC(i,j) "  -> neutralino i chargino j (i=1..5, j=1,2)
*        HCBRSQ(i)   "  -> uLdL, t1b1, t1b2, t2b1, t2b2 (i=1..5)
*        HCBRSL(i)   "  -> lLnL, t1nt, t2nt (i=1..3)
*
*      MNEU(i)   Mass of neutralino chi_i (i=1,5, ordered in mass)
*      NEU(i,j)  chi_i components of bino, wino, higgsino u&d, singlino
*                (i,j=1..5)
*
*      MCHA(i)       Chargino masses
*      U(i,j),V(i,j) Chargino mixing matrices
*
*  ERRORS: IFAIL = 0..22
*
*  IFAIL = 0         OK
*          1,3,5,7   m_h1^2 < 0
*          2,3,6,7   m_a1^2 < 0
*          4,5,6,7   m_h+^2 < 0
*          8         m_sfermion^2 < 0
*          10        Violation of phenomenological constraint(s)
*          11-15     Problem in integration of RGEs
*          16-19     Convergence problem
*          20        No electroweak symmetry breaking
*          21        MSMES=/=MSREF (GMSB scenario above MMESS)
*          22        IFAIL = 10 & 21
*
*  Phenomenological constraints:
*
*      PROB(I)  = 0, I = 1..82: OK
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
*      PROB(27) =/= 0  Landau Pole in l, k, ht, hb below MGUT
*      PROB(28) =/= 0  unphysical global minimum
*      PROB(29) =/= 0  Higgs soft masses >> Msusy
*      PROB(30) =/= 0  excluded by DM relic density (checked only if OMGFLAG=/=0)
*      PROB(32) =/= 0  b->s gamma more than 2 sigma away
*      PROB(33) =/= 0  Delta M_s more than 2 sigma away
*      PROB(34) =/= 0  Delta M_d more than 2 sigma away
*      PROB(35) =/= 0  B_s->mu+mu- more than 2 sigma away
*      PROB(36) =/= 0  B+-> tau+nu_tau more than 2 sigma away
*      PROB(37) =/= 0  (g-2)_muon more than 2 sigma away
*      PROB(38) =/= 0  excluded by Upsilon(1S) -> H/A gamma
*      PROB(39) =/= 0  excluded by eta_b(1S) mass measurement
*      PROB(40) =/= 0  BR(B-->X_s mu+ mu-) more than 2 sigma away
*      PROB(41) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus (ALEPH analysis)
*      PROB(42) =/= 0  excluded by top -> b H+, H+ -> c s (CDF, D0)
*      PROB(43) =/= 0  excluded by top -> b H+, H+ -> tau nu_tau (D0)
*      PROB(44) =/= 0  excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus (CDF)
*      PROB(45) =/= 0  excluded by t -> bH+ (ATLAS)
*      PROB(46) =/= 0  No Higgs in the MHmin-MHmax GeV range
*      PROB(51) =/= 0: excluded by H/A->tautau (ATLAS+CMS)
*      PROB(52) =/= 0: excluded by H->AA->4leptons/2lept.+2b (ATLAS+CMS)
*      PROB(53) =/= 0: excluded by ggF->H/A->gamgam (ATLAS)
*      PROB(55) =/= 0: b -> d gamma more than 2 sigma away
*      PROB(56) =/= 0: B_d -> mu+ mu- more than 2 sigma away
*      PROB(57) =/= 0: b -> s nu nubar more than 2 sigma away
*      PROB(58) =/= 0: b -> c tau nu more than 2 sigma away (as SM)
*      PROB(59) =/= 0: K -> pi nu nubar more than 2 sigma away
*      PROB(60) =/= 0: DMK / epsK more than 2 sigma away
*      PROB(63) =/= 0: excluded by H->AA->4gammas (ATLAS+CMS)
*      PROB(65) =/= 0: excluded by light mesons or eta_{c,b} decays
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
      PARAMETER (NFL=22,NPROB=82+1,NPAR=25)
      INTEGER NFAIL(NFL),IFAIL,I,TOT,ITOT,NTOT,JM,JL
      INTEGER IMSUSYEFF,IMMESS,ITB,IL,IAL,IMUP,IMSP
      INTEGER IXIU,ILPP,ILTT,ILU,ILD,ILT,ILB,ILL
      INTEGER NMSUSYEFF,NMMESS,NTB,NL,NAL,NMUP,NMSP
      INTEGER NXIU,NLPP,NLTT,NLU,NLD,NLT,NLB,NLL
      INTEGER N1,N2,I1,I2,ITER,ITERMU,Q2FIX,LOOP,NMSFLAG
      INTEGER OMGFLAG,MAFLAG,MOFLAG,GMUFLAG,HFLAG,GMFLAG
      INTEGER NFLAG,UNCERTFLAG,ITRY,SFFLAG,GRFLAG,CFLAG(5)

      DOUBLE PRECISION PAR(NPAR),PROB(NPROB),CHECK,MESTEST,GUTEST
      DOUBLE PRECISION MSUSYEFFMIN,MSUSYEFFMAX,MMESSMIN,MMESSMAX
      DOUBLE PRECISION TBMIN,TBMAX,LMIN,LMAX,KMIN,KMAX,ALMIN,ALMAX
      DOUBLE PRECISION XIFMIN,XIFMAX,XISMIN,XISMAX,MUPMIN,MUPMAX
      DOUBLE PRECISION MSPMIN,MSPMAX,MSMIN,MSMAX,XIUMIN,XIUMAX
      DOUBLE PRECISION LPPMIN,LPPMAX,LTTMIN,LTTMAX,LUMIN,LUMAX,LDMIN
      DOUBLE PRECISION LDMAX,LTMIN,LTMAX,LBMIN,LBMAX,LLMIN,LLMAX
      DOUBLE PRECISION MSUSYEFFN,MSUSYEFFNN,MMESSN,MMESSNN,TBN,TBNN
      DOUBLE PRECISION LN,LNN,KN,KNN,ALN,ALNN,XIFN,XIFNN,XISN,XISNN
      DOUBLE PRECISION MUPN,MUPNN,MSPN,MSPNN,MSN,MSNN,XIUN,XIUNN
      DOUBLE PRECISION LPPN,LPPNN,LTTN,LTTNN,LUN,LUNN,LDN,LDNN
      DOUBLE PRECISION LTN,LTNN,LBN,LBNN,LLN,LLNN,MUN,MUNN
      DOUBLE PRECISION XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      DOUBLE PRECISION ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      DOUBLE PRECISION MSUSYEFF,MMESS,N5,DETM,SIGMU,Q2,Q2MIN
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW,MSREF,D,DMIN
      DOUBLE PRECISION PI,M32,CGR,MPL
      DOUBLE PRECISION G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      DOUBLE PRECISION LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU
      DOUBLE PRECISION MSM,MST,LM,LT,KM,KT,HTM,HTT,LPPM,LPPT,LTTM,LTTT
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTGUT,HBGUT,HLGUT
      DOUBLE PRECISION LPPGUT,LTTGUT,LUGUT,LDGUT,LTGUT,LBGUT,LLGUT
      DOUBLE PRECISION DELMB,DELML,DEL1,MUFAIL,MUSTEP,MUINIT
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION SMASS(3),PMASS(2),CMASS,SCOMP(3,3),PCOMP(2,2)
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU

      COMMON/NMSFLAG/NMSFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG
      COMMON/SIGMU/SIGMU
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/DETM/DETM
      COMMON/MUFAIL/MUFAIL
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/MESCOUP/G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      COMMON/MESEXT/XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      COMMON/MESGUT/LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTGUT,HBGUT,HLGUT
      COMMON/GUTMES/LPPGUT,LTTGUT,LUGUT,LDGUT,LTGUT,LBGUT,LLGUT
      COMMON/INPPAR/ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/MINMAX/MSUSYEFFMIN,MSUSYEFFMAX,MMESSMIN,MMESSMAX,
     . TBMIN,TBMAX,LMIN,LMAX,KMIN,KMAX,ALMIN,ALMAX,XIFMIN,
     . XIFMAX,XISMIN,XISMAX,MUPMIN,MUPMAX,MSPMIN,MSPMAX,
     . MSMIN,MSMAX,XIUMIN,XIUMAX,LPPMIN,LPPMAX,LTTMIN,LTTMAX,LUMIN,
     . LUMAX,LDMIN,LDMAX,LTMIN,LTMAX,LBMIN,LBMAX,LLMIN,LLMAX
      COMMON/BOUNDS/MSUSYEFFN,MSUSYEFFNN,MMESSN,MMESSNN,TBN,TBNN,
     . LN,LNN,KN,KNN,ALN,ALNN,XIFN,XIFNN,XISN,XISNN,MUPN,MUPNN,
     . MSPN,MSPNN,MSN,MSNN,XIUN,XIUNN,LPPN,LPPNN,LTTN,LTTNN,LUN,LUNN,
     . LDN,LDNN,LTN,LTNN,LBN,LBNN,LLN,LLNN,MUN,MUNN
      COMMON/STEPS/NTOT,NMSUSYEFF,NMMESS,NTB,NL,NAL,NMUP,NMSP,
     . NXIU,NLPP,NLTT,NLU,NLD,NLT,NLB,NLL,N1,N2
      COMMON/GMSAVE/MSM,MST,LM,LT,KM,KT,HTM,HTT,
     . LPPM,LPPT,LTTM,LTTT,JM,JL
      COMMON/M32/M32,CGR,MPL,GRFLAG
      COMMON/UNCERTFLAG/UNCERTFLAG
      COMMON/CFLAG/CFLAG
      COMMON/NFLAG/NFLAG
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     . MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     . CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU

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

      MSUSYEFFN=1d99
      MSUSYEFFNN=-1d99
      MMESSN=1d99
      MMESSNN=-1d99
      TBN=1d99
      TBNN=-1d99
      LN=1d99
      LNN=-1d99
      ALN=1d99
      ALNN=-1d99
      MUN=1d99
      MUNN=-1d99
      XIFN=1d99
      XIFNN=-1d99
      KN=1d99
      KNN=-1d99
      MSN=1d99
      MSNN=-1d99
      XISN=1d99
      XISNN=-1d99
      MUPN=1d99
      MUPNN=-1d99
      MSPN=1d99
      MSPNN=-1d99
      XIUN=1d99
      XIUNN=-1d99
      LPPN=1d99
      LPPNN=-1d99
      LTTN=1d99
      LTTNN=-1d99
      LUN=1d99
      LUNN=-1d99
      LDN=1d99
      LDNN=-1d99
      LTN=1d99
      LTNN=-1d99
      LBN=1d99
      LBNN=-1d99
      LLN=1d99
      LLNN=-1d99

*   Beginning of the scan

      DO IXIU=1,NXIU
      IF(GMFLAG.EQ.1 .OR.GMFLAG.EQ.3)THEN
       IF(NXIU.EQ.1)THEN
        XIU=XIUMIN
       ELSE
        XIU=XIUMIN+(XIUMAX-XIUMIN)*DFLOAT(IXIU-1)/DFLOAT(NXIU-1)
       ENDIF
      ENDIF

      DO ILPP=1,NLPP
      IF(GMFLAG.EQ.2 .OR.GMFLAG.EQ.4)THEN
       IF(NLPP.EQ.1)THEN
        LPPMES=LPPMIN
       ELSE
        LPPMES=LPPMIN+(LPPMAX-LPPMIN)*DFLOAT(ILPP-1)/DFLOAT(NLPP-1)
       ENDIF
      ENDIF

      DO ILTT=1,NLTT
      IF(GMFLAG.EQ.2 .OR.GMFLAG.EQ.4)THEN
       IF(NLTT.EQ.1)THEN
        LTTMES=LTTMIN
       ELSE
        LTTMES=LTTMIN+(LTTMAX-LTTMIN)*DFLOAT(ILTT-1)/DFLOAT(NLTT-1)
       ENDIF
      ENDIF

      DO ILU=1,NLU
      IF(NLU.EQ.1)THEN
       LUMES=LUMIN
      ELSE
       LUMES=LUMIN+(LUMAX-LUMIN)*DFLOAT(ILU-1)/DFLOAT(NLU-1)
      ENDIF

      DO ILD=1,NLD
      IF(NLD.EQ.1)THEN
       LDMES=LDMIN
      ELSE
       LDMES=LDMIN+(LDMAX-LDMIN)*DFLOAT(ILD-1)/DFLOAT(NLD-1)
      ENDIF

      DO ILT=1,NLT
      IF(NLT.EQ.1)THEN
       LTMES=LTMIN
      ELSE
       LTMES=LTMIN+(LTMAX-LTMIN)*DFLOAT(ILT-1)/DFLOAT(NLT-1)
      ENDIF

      DO ILB=1,NLB
      IF(NLB.EQ.1)THEN
       LBMES=LBMIN
      ELSE
       LBMES=LBMIN+(LBMAX-LBMIN)*DFLOAT(ILB-1)/DFLOAT(NLB-1)
      ENDIF

      DO ILL=1,NLL
      IF(NLL.EQ.1)THEN
       LLMES=LLMIN
      ELSE
       LLMES=LLMIN+(LLMAX-LLMIN)*DFLOAT(ILL-1)/DFLOAT(NLL-1)
      ENDIF

      DO IAL=1,NAL
      IF(NAL.EQ.1)THEN
       IF(GMFLAG.EQ.0)THEN
        ALINP=ALMIN
       ENDIF
      ELSE
       ALINP=ALMIN+(ALMAX-ALMIN)*DFLOAT(IAL-1)/DFLOAT(NAL-1)
      ENDIF

      DO IMSUSYEFF=1,NMSUSYEFF
      IF(NMSUSYEFF.EQ.1)THEN
       MSUSYEFF=MSUSYEFFMIN
      ELSE
       MSUSYEFF=MSUSYEFFMIN+(MSUSYEFFMAX-MSUSYEFFMIN)
     .         *DFLOAT(IMSUSYEFF-1)/DFLOAT(NMSUSYEFF-1)
      ENDIF

      DO IMMESS=1,NMMESS
      IF(NMMESS.EQ.1)THEN
       MMESS=MMESSMIN
      ELSE
       MMESS=MMESSMIN+(MMESSMAX-MMESSMIN)
     .      *DFLOAT(IMMESS-1)/DFLOAT(NMMESS-1)
      ENDIF

      DO ITB=1,NTB
      IF(NTB.EQ.1)THEN
       PAR(3)=TBMIN
      ELSE
       PAR(3)=TBMIN+(TBMAX-TBMIN)*DFLOAT(ITB-1)/DFLOAT(NTB-1)
      ENDIF

      DO IL=1,NL
      IF(NL.EQ.1)THEN
       PAR(1)=LMIN
      ELSE
       PAR(1)=LMIN+(LMAX-LMIN)*DFLOAT(IL-1)/DFLOAT(NL-1)
      ENDIF

      DO I1=1,N1
      IF(N1.EQ.1)THEN
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        XIFINP=XIFMIN
       ELSE
        PAR(2)=KMIN
       ENDIF
      ELSE
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        XIFINP=XIFMIN+(XIFMAX-XIFMIN)*DFLOAT(I1-1)/DFLOAT(N1-1)
       ELSE
        PAR(2)=KMIN+(KMAX-KMIN)*DFLOAT(I1-1)/DFLOAT(N1-1)
       ENDIF
      ENDIF

      DO I2=1,N2
      IF(N2.EQ.1)THEN
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        XISINP=XISMIN
       ELSE
        MSINP=MSMIN
       ENDIF
      ELSE
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        XISINP=XISMIN+(XISMAX-XISMIN)*DFLOAT(I2-1)/DFLOAT(N2-1)
       ELSE
        MSINP=MSMIN+(MSMAX-MSMIN)*DFLOAT(I2-1)/DFLOAT(N2-1)
       ENDIF
      ENDIF

      DO IMUP=1,NMUP
      IF(NMUP.EQ.1)THEN
       MUPINP=MUPMIN
      ELSE
       MUPINP=MUPMIN+(MUPMAX-MUPMIN)*DFLOAT(IMUP-1)/DFLOAT(NMUP-1)
      ENDIF

      DO IMSP=1,NMSP
      IF(NMSP.EQ.1)THEN
       MSPINP=MSPMIN
      ELSE
       MSPINP=MSPMIN+(MSPMAX-MSPMIN)*DFLOAT(IMSP-1)/DFLOAT(NMSP-1)
      ENDIF

      ITOT=ITOT+1

!      WRITE(0,*)""
!      WRITE(0,*)"------------------------------------------------------"
!      WRITE(0,*)""
!      WRITE(0,*)"Point ",ITOT
!      WRITE(0,*)""
!      WRITE(0,*)"MAFLAG=",MAFLAG
!      WRITE(0,*)"GMFLAG=",GMFLAG
!      WRITE(0,*)"MSUSYEFF =",MSUSYEFF
!      WRITE(0,*)"MMESS =",MMESS
!      WRITE(0,*)"TANB =",PAR(3)
!      WRITE(0,*)"LAMBDA =",PAR(1)
!      IF(GMFLAG.EQ.0)THEN
!       WRITE(0,*)"ALAMBDA=",ALINP
!       IF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4)WRITE(0,*)"KAPPA =",PAR(2)
!       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)WRITE(0,*)"XIF =",XIFINP
!       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)WRITE(0,*)"XIS =",XISINP
!       IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)WRITE(0,*)"MS =",MSINP
!       IF(MUPINP.NE.0d0)WRITE(0,*)"MUP =",MUPINP
!       IF(MSPINP.NE.0d0)WRITE(0,*)"MSP =",MSPINP
!      ENDIF
!      IF(GMFLAG.EQ.1 .OR. GMFLAG.EQ.3)THEN
!       WRITE(0,*)"XiU =",XIU
!      ENDIF
!      IF(GMFLAG.EQ.2 .OR. GMFLAG.EQ.4)THEN
!       WRITE(0,*)"LPP =",LPPMES
!       WRITE(0,*)"LTT =",LTTMES
!      ENDIF
!      IF(GMFLAG.EQ.1 .OR. GMFLAG.EQ.2)THEN
!       WRITE(0,*)"LU =",LUMES
!       WRITE(0,*)"LD =",LDMES
!       WRITE(0,*)"LT =",LTMES
!       WRITE(0,*)"LB =",LBMES
!       WRITE(0,*)"LL =",LLMES
!      ENDIF
!      IF(GMFLAG.EQ.3 .OR. GMFLAG.EQ.4)THEN
!       WRITE(0,*)"LT =",LTMES
!       WRITE(0,*)"LB =",LBMES
!      ENDIF
!      WRITE(0,*)""
!      WRITE(0,*)""

*   Initialization of PROB

      DO I=1,NPROB
       PROB(I)=0d0
      ENDDO
      UNCERTFLAG=0
      IFAIL=-1

*   Initialization of algorithm parameters

      LPPT=0d0
      LTTT=0d0
      JL=0
      LOOP=0
 2    LOOP=LOOP+1
!      IF(GMFLAG.EQ.1 .OR. GMFLAG.EQ.3)THEN
!       WRITE(0,*)"LOOP =",LOOP
!       WRITE(0,*)""
!       WRITE(0,*)""
!      ENDIF
      SFFLAG=0
      ITRY=1
      MUINIT=.5d0*SIGMU*N5*MSUSYEFF/(4d0*PI)**2
      MUSTEP=MUINIT
 1    MUINIT=MUINIT+MUSTEP
      MUFAIL=MUINIT-MUSTEP
!      WRITE(0,*)"ITRY=",ITRY
!      WRITE(0,*)""
!      WRITE(0,*)"MUINIT=",MUINIT
!      WRITE(0,*)"MUSTEP=",MUSTEP
!      WRITE(0,*)"MUFAIL=",MUFAIL
!      WRITE(0,*)""
!      WRITE(0,*)""
      MST=0d0
      LT=0d0
      KT=0d0
      HTT=0d0
      JM=0

*   Guess parameters at Q2/MMESS

      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"Guesses"
       PAR(4)=MUINIT
       IF(Q2FIX.EQ.0)THEN
        Q2=MAX(N5*MSUSYEFF**2/(4d0*PI)**4,Q2MIN)
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        PAR(2)=PAR(1)/5d0
       ELSE
        XIFMES=MSUSYEFF**2/(4d0*PI)**4
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        MSM=MSUSYEFF**2/(4d0*PI)**4
       ELSE
        XISMES=MSUSYEFF**3/(4d0*PI)**6
       ENDIF
       IF(GMFLAG.EQ.0)THEN
        PAR(5)=ALINP
       ELSE
        PAR(5)=-MSUSYEFF/(4d0*PI)**2
       ENDIF
       IF(PAR(2).NE.0d0)THEN
        PAR(6)=3d0*PAR(5)
       ELSE
        PAR(6)=0d0
       ENDIF
       DO I=7,11
        PAR(I)=N5*MSUSYEFF**2/(4d0*PI)**4
       ENDDO
       DO I=12,14
        PAR(I)=0d0
       ENDDO
       DO I=15,19
        PAR(I)=N5*MSUSYEFF**2/(4d0*PI)**4
       ENDDO
        DO I=20,22
        PAR(I)=N5*MSUSYEFF/(4d0*PI)**2
       ENDDO
       PAR(23)=MSUSYEFF/(4d0*PI)**2
       PAR(24)=MSUSYEFF/(4d0*PI)**2
       PAR(25)=0d0
       DELMB=.1d0
       DELML=0d0
       DEL1=0d0
       IF(GMFLAG.EQ.1 .OR. GMFLAG.EQ.3)THEN
        LPPM=XIU
        LTTM=XIU
       ENDIF
       IF(GMFLAG.EQ.2 .OR. GMFLAG.EQ.4)THEN
        LPPM=LPPMES
        LTTM=LTTMES
       ENDIF
!      ELSE
!       WRITE(0,*)"Current"
      ENDIF
      IFAIL=0

!      WRITE(0,*)""
!      WRITE(0,*)"MU =",PAR(4)
!      WRITE(0,*)"M1 =",PAR(20)
!      WRITE(0,*)"M2 =",PAR(21)
!      WRITE(0,*)"M3 =",PAR(22)
!      WRITE(0,*)"AL =",PAR(5)
!      WRITE(0,*)"AK =",PAR(6)
!      WRITE(0,*)"ATOP =",PAR(12)
!      WRITE(0,*)"ABOT =",PAR(13)
!      WRITE(0,*)"ATAU =",PAR(14)
!      WRITE(0,*)"AMUON =",PAR(25)
!      WRITE(0,*)"MQ3 =",PAR(7)
!      WRITE(0,*)"MU3 =",PAR(8)
!      WRITE(0,*)"MD3 =",PAR(9)
!      WRITE(0,*)"MQ =",PAR(15)
!      WRITE(0,*)"MU =",PAR(16)
!      WRITE(0,*)"MD =",PAR(17)
!      WRITE(0,*)"ML3 =",PAR(10)
!      WRITE(0,*)"ME3 =",PAR(11)
!      WRITE(0,*)"ML =",PAR(18)
!      WRITE(0,*)"ME =",PAR(19)
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
!       WRITE(0,*)"K =",PAR(2)
!      ELSE
!      WRITE(0,*)"XIF =",XIFMES
!       ENDIF
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
!       WRITE(0,*)"MS =",MSM
!      ELSE
!       WRITE(0,*)"XIS =",XISMES
!      ENDIF
!      IF(GMFLAG.EQ.1 .OR. GMFLAG.EQ.3)THEN
!       WRITE(0,*)"LPP =",LPPM
!       WRITE(0,*)"LTT =",LTTM
!      ENDIF
!      WRITE(0,*)""

*   Guess for couplings at MMESS

      CALL RGESGM(PAR,IFAIL)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"RGE integration problem 1"
       IF(DABS(MUINIT).LT.1d2 .AND. DABS(MUSTEP).GT.1d0)THEN
        IF(SFFLAG.EQ.1)MUSTEP=MUSTEP/2d0
        ITRY=ITRY+1
        GOTO 1
       ELSE
!        WRITE(0,*)"Exit (IFAIL=11)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=11
        GOTO 11
       ENDIF
      ENDIF

*   External loop to compute the soft parameters at Q2

      ITER=0
 21   ITER=ITER+1
!      WRITE(0,*)"ITER =",ITER
!      WRITE(0,*)""
!      WRITE(0,*)""

      CALL RGESINVGM(PAR,IFAIL)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"RGE integration problem 2"
       IF(DABS(MUINIT).LT.1d2 .AND. DABS(MUSTEP).GT.1d0)THEN
        IF(SFFLAG.EQ.1)MUSTEP=MUSTEP/2d0
        ITRY=ITRY+1
        GOTO 1
       ELSE
!        WRITE(0,*)"Exit (IFAIL=13)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=13
        GOTO 11
       ENDIF
      ENDIF

*   Internal loop to compute ((k or XiF) and (mS or XiS) and mu)

      ITERMU=0
 22   ITERMU=ITERMU+1
!      WRITE(0,*)"ITERMU =",ITERMU
!      WRITE(0,*)""
!      WRITE(0,*)""

      CALL RUNPAR(PAR)

      CALL MSFERM(PAR,IFAIL,0)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"Negative sfermion mass squared"
       IF(DABS(MUINIT).LT.1d2 .AND. DABS(MUSTEP).GT.1d0)THEN
        MUSTEP=MUSTEP/2d0
        MUINIT=MUSTEP
        SFFLAG=1
        ITRY=ITRY+1
        GOTO 1
       ELSE
!        WRITE(0,*)"Exit (IFAIL=8)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=8
        GOTO 11
       ENDIF
      ENDIF

      CALL MINIMIZE(PAR,CHECK)

*   End of the internal loop

      IF((CHECK.GT.1d-12.AND.ITERMU.LT.10).OR.
     .   (CHECK.GT.1d-8.AND.ITERMU.LT.50).OR.
     .   (CHECK.GT.1d-6.AND.ITERMU.LT.100))GOTO 22
      IF(CHECK.GT.1d-4)THEN
!       WRITE(0,*)"No convergence 1"
!       WRITE(0,*)"Exit (IFAIL=16)"
!       WRITE(0,*)""
!       WRITE(0,*)""
       IFAIL=16
       GOTO 11
      ENDIF

      CALL RGESUNIGM(PAR,IFAIL,MESTEST)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"RGE integration problem 3"
       IF(DABS(MUINIT).LT.1d2 .AND. DABS(MUSTEP).GT.1d0)THEN
        IF(SFFLAG.EQ.1)MUSTEP=MUSTEP/2d0
        ITRY=ITRY+1
        GOTO 1
       ELSE
!        WRITE(0,*)"Exit (IFAIL=12)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=12
        GOTO 11
       ENDIF
      ENDIF

*   End of the external loop

      IF((MESTEST.GT.1d-12.AND.ITER.LT.10).OR.
     .   (MESTEST.GT.1d-8.AND.ITER.LT.50).OR.
     .   (MESTEST.GT.1d-6.AND.ITER.LT.100))GOTO 21
      IF(MESTEST.GT.1d-4)THEN
!       WRITE(0,*)"No convergence 2"
       IF(DABS(MUSTEP).GT.1d0)THEN
        MUSTEP=MUSTEP/2d0
        MUFAIL=MUFAIL+MUSTEP
        MST=0d0
        LT=0d0
        KT=0d0
        HTT=0d0
        JM=0
        ITER=0
        ITRY=ITRY+1
!        WRITE(0,*)"ITRY=",ITRY
!        WRITE(0,*)""
!        WRITE(0,*)"MUINIT=",MUINIT
!        WRITE(0,*)"MUSTEP=",MUSTEP
!        WRITE(0,*)"MUFAIL=",MUFAIL
!        WRITE(0,*)""
!        WRITE(0,*)""
        GOTO 21
       ELSE
!        WRITE(0,*)"Exit (IFAIL=17)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=17
        GOTO 11
       ENDIF
      ENDIF

*   Check if correct EWSB

      IF(DETM.LE.0d0)THEN
!       WRITE(0,*)"Convergence in a false minimum"
       IF(DABS(MUSTEP).GT.1d0)THEN
        MUSTEP=MUSTEP/2d0
        MUFAIL=MUFAIL-MUSTEP
        MST=0d0
        LT=0d0
        KT=0d0
        HTT=0d0
        JM=0
        ITER=0
        ITRY=ITRY+1
!        WRITE(0,*)"ITRY=",ITRY
!        WRITE(0,*)""
!        WRITE(0,*)"MUINIT=",MUINIT
!        WRITE(0,*)"MUSTEP=",MUSTEP
!        WRITE(0,*)"MUFAIL=",MUFAIL
!        WRITE(0,*)""
!        WRITE(0,*)""
        GOTO 21
       ELSEIF(DETM.LT.0d0)THEN
!        WRITE(0,*)"Exit (IFAIL=20)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=20
        GOTO 11
       ENDIF
      ENDIF

*  GUT scale

      CALL RGESGMGUT(PROB,IFAIL,GUTEST)
      IF(IFAIL.GT.0 .AND. GMFLAG.NE.0)THEN
!       WRITE(0,*)"RGE integration problem 4"
!       WRITE(0,*)"Exit (IFAIL=14)"
!       WRITE(0,*)""
!       WRITE(0,*)""
       IFAIL=14
       GOTO 11
      ELSE
       IFAIL=0
      ENDIF
      IF(GMFLAG.EQ.1 .OR. GMFLAG.EQ.3)THEN
!       WRITE(0,*)"GUTEST =",GUTEST
!       WRITE(0,*)""
!       WRITE(0,*)""
       IF((GUTEST.GT.1d-12.AND.LOOP.LT.10).OR.
     .    (GUTEST.GT.1d-8.AND.LOOP.LT.50).OR.
     .    (GUTEST.GT.1d-6.AND.LOOP.LT.100))THEN
        CALL LPTSOLVE(IFAIL)
        IF(IFAIL.NE.0)THEN
!         WRITE(0,*)"RGE integration problem 5"
!         WRITE(0,*)"Exit (IFAIL=15)"
!         WRITE(0,*)""
!         WRITE(0,*)""
         IFAIL=15
         GOTO 11
        ENDIF
        GOTO 2
       ELSEIF(GUTEST.GT.1d-4)THEN
!        WRITE(0,*)"No convergence 3"
!        WRITE(0,*)"Exit (IFAIL=18)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=18
        GOTO 11
       ENDIF
      ENDIF

*   Computation of sfermion masses

      CALL MSFERM(PAR,IFAIL,1)
      IF(IFAIL.NE.0)GOTO 11

*   Computation of Higgs masses

      CALL MHIGGS(PAR,PROB,IFAIL)
      IF(IFAIL.EQ.-1)IFAIL=19
      IF(IFAIL.NE.0)GOTO 11

*   Computation of gluino mass

      CALL GLUINO(PAR)

*   Computation of chargino masses

      CALL CHARGINO(PAR)

*   Computation of neutralino masses

      CALL NEUTRALINO(PAR)
      IF(NFLAG.EQ.1)THEN
       PROB(NPROB)=DDIM(DABS(MNEU(1))/MSL1,1d0)
       PROB(NPROB)=PROB(NPROB)+DDIM(DABS(MNEU(1))/MLR,1d0)
      ELSEIF(NFLAG.EQ.2)THEN
       PROB(NPROB)=DDIM(MSL1/DABS(MNEU(1)),1d0)
       PROB(NPROB)=PROB(NPROB)+DDIM(MSL1/MLR,1d0)
      ELSEIF(NFLAG.EQ.3)THEN
       PROB(NPROB)=DDIM(MLR/DABS(MNEU(1)),1d0)
       PROB(NPROB)=PROB(NPROB)+DDIM(MLR/MSL1,1d0)
      ENDIF

*   Computation of Higgs + top branching ratios

      CALL DECAY(PAR,PROB)
      IF(CFLAG(4).EQ.0)PROB(65)=0d0
      CALL TDECAY(PAR)

*   Exp. constraints on sparticles/Higgses

      IF(CFLAG(2).EQ.1)CALL SUBEXP(PAR,PROB)
      IF(CFLAG(3).EQ.1)THEN
       CALL LHCHIG(PROB)
       PROB(69)=0d0
       PROB(72)=0d0
      ENDIF

*   B + K physics

      IF(CFLAG(4).EQ.1)THEN
       CALL BOTTOMONIUM(PROB)
       CALL BSG(PAR,PROB)
       CALL KPHYS(PAR,PROB)
       PROB(58)=0d0
      ENDIF

*   Anom. magn. moment of the Muon

      IF(GMUFLAG.EQ.1)CALL MAGNMU(PAR,PROB)

*   Global minimum?

      CALL CHECKMIN(PROB)

*   Relic density

      M32=CGR*MSUSYEFF*MMESS/DSQRT(3d0)/MPL
      CALL RELDEN(PAR,PROB)

*   Check for problems

      DO I=1,NPROB
       IF(PROB(I).NE.0d0)THEN
!        WRITE(0,*)"PROB",I,PROB(I)
        IFAIL=10
       ENDIF
      ENDDO
!      WRITE(0,*)""

*   Check MSMES=MSREF if GMFLAG>0

      IF(GMFLAG.NE.0)THEN
       MSREF=MSUSYEFF**2/(4d0*PI)**4*(
     .      - 2d0*G1MES*LDMES**2 - 6d0*G2MES*LDMES**2
     .      + 6d0*HBMES*LDMES**2 + 2d0*HLMES*LDMES**2
     .      - 8d0*KMES**2*LDMES**2 + 8d0*LDMES**4 + 4d0*LDMES**2*LMES
     .      + 12d0*DSQRT(HBMES)*LBMES*LDMES*LPPMES
     .      + 4d0*DSQRT(HLMES)*LDMES*LLMES*LPPMES
     .      - 2d0*G1MES*LPPMES**2 - 6d0*G2MES*LPPMES**2
     .      - 8d0*KMES**2*LPPMES**2 + 6d0*LBMES**2*LPPMES**2
     .      + 16d0*LDMES**2*LPPMES**2 + 2d0*LLMES**2*LPPMES**2
     .      + 8d0*LPPMES**4 + 12d0*DSQRT(HTMES*LMES)*LDMES*LTMES
     .      + 6d0*LDMES**2*LTMES**2 + 6d0*LPPMES**2*LTMES**2
     .      - 4d0/3d0*G1MES*LTTMES**2 - 16d0*G3MES*LTTMES**2
     .      - 12d0*KMES**2*LTTMES**2 + 12d0*LDMES**2*LTTMES**2
     .      + 12d0*LPPMES**2*LTTMES**2 + 15d0*LTTMES**4
     .      + 12d0*DSQRT(HBMES*LMES)*LBMES*LUMES
     .      + 4d0*DSQRT(HLMES*LMES)*LLMES*LUMES
     .      + 16d0*LDMES*DSQRT(LMES)*LPPMES*LUMES
     .      + 12d0*DSQRT(HTMES)*LPPMES*LTMES*LUMES
     .      - 2d0*G1MES*LUMES**2 - 6d0*G2MES*LUMES**2
     .      + 6d0*HTMES*LUMES**2 - 8d0*KMES**2*LUMES**2
     .      + 6d0*LBMES**2*LUMES**2 + 8d0*LDMES**2*LUMES**2
     .      + 2d0*LLMES**2*LUMES**2 + 4d0*LMES*LUMES**2
     .      + 16d0*LPPMES**2*LUMES**2 + 12d0*LTTMES**2*LUMES**2
     .      + 8d0*LUMES**4)
     .      -MSUSYEFF**4/(48d0*PI**2*MMESS**2)
     .       *(2d0*LPPMES**2+3d0*LTTMES**2+LUMES**2+LDMES**2)
       D=(MSMES-MSREF)/(1d0+DABS(MSMES)+DABS(MSREF))
       IF(DABS(D).GT.DMIN)THEN
        IF(IFAIL.EQ.0)IFAIL=21
        IF(IFAIL.EQ.10)IFAIL=22
       ENDIF
!       WRITE(0,*)"MSREF =",MSREF
!       WRITE(0,*)"MSMES =",MSMES
!       WRITE(0,*)"D =",D
!       WRITE(0,*)""
!       WRITE(0,*)""
      ENDIF
      IF(IFAIL.NE.0)GOTO 11

*   Sparticle decays

      IF(NMSFLAG.NE.0)CALL NMSDECAY(PAR)

*   Computation of the fine-tuning

      CALL FTPAR(PAR,2)

*   Recording of the results

11    CALL OUTPUT(PAR,PROB,IFAIL)
      IF(IFAIL.EQ.0)THEN
       TOT=TOT+1
       MSUSYEFFN=MIN(MSUSYEFF,MSUSYEFFN)
       MSUSYEFFNN=MAX(MSUSYEFF,MSUSYEFFNN)
       MMESSN=MIN(MMESS,MMESSN)
       MMESSNN=MAX(MMESS,MMESSNN)
       TBN=MIN(PAR(3),TBN)
       TBNN=MAX(PAR(3),TBNN)
       LN=MIN(PAR(1),LN)
       LNN=MAX(PAR(1),LNN)
       KN=MIN(PAR(2),KN)
       KNN=MAX(PAR(2),KNN)
       ALN=MIN(ALINP,ALN)
       ALNN=MAX(ALINP,ALNN)
       MUN=MIN(PAR(4),MUN)
       MUNN=MAX(PAR(4),MUNN)
       XIFN=MIN(XIFMES,XIFN)
       XIFNN=MAX(XIFMES,XIFNN)
       XISN=MIN(XISMES,XISN)
       XISNN=MAX(XISMES,XISNN)
       MSN=MIN(MSMES,MSN)
       MSNN=MAX(MSMES,MSNN)
       MUPN=MIN(MUPMES,MUPN)
       MUPNN=MAX(MUPMES,MUPNN)
       MSPN=MIN(MSPMES,MSPN)
       MSPNN=MAX(MSPMES,MSPNN)
       XIUN=MIN(XIU,XIUN)
       XIUNN=MAX(XIU,XIUNN)
       LPPN=MIN(LPPMES,LPPN)
       LPPNN=MAX(LPPMES,LPPNN)
       LTTN=MIN(LTTMES,LTTN)
       LTTNN=MAX(LTTMES,LTTNN)
       LUN=MIN(LUMES,LUN)
       LUNN=MAX(LUMES,LUNN)
       LDN=MIN(LDMES,LDN)
       LDNN=MAX(LDMES,LDNN)
       LTN=MIN(LTMES,LTN)
       LTNN=MAX(LTMES,LTNN)
       LBN=MIN(LBMES,LBN)
       LBNN=MAX(LBMES,LBNN)
       LLN=MIN(LLMES,LLN)
       LLNN=MAX(LLMES,LLNN)
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

*   Summary of the scan:
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

      INTEGER I,NLINE,INL,ICH,IX,IVAL,Q2FIX
      INTEGER N0,NLOOP,NBER,NPAR,ERR,VFLAG,NMSFLAG,Z3FLAG
      INTEGER OMGFLAG,MAFLAG,MOFLAG,PFLAG,GMUFLAG,HFLAG
      INTEGER GMFLAG,NFLAG,GRFLAG,CFLAG(5)
      INTEGER NTOT,NMSUSYEFF,NMMESS,NTB,NL,NAL,NMUP,NMSP
      INTEGER NXIU,NLPP,NLTT,NLU,NLD,NLT,NLB,NLL
      INTEGER N1,N2,NK,NXIF,NXIS,NMS

      DOUBLE PRECISION PAR(*),VAL
      DOUBLE PRECISION ACC,XITLA,XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW,M32,CGR,MPL
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,Q2,SIGMU,Q2MIN
      DOUBLE PRECISION MSUSYEFFMIN,MSUSYEFFMAX,MMESSMIN,MMESSMAX
      DOUBLE PRECISION TBMIN,TBMAX,LMIN,LMAX,KMIN,KMAX,ALMIN,ALMAX
      DOUBLE PRECISION XIFMIN,XIFMAX,XISMIN,XISMAX,MUPMIN,MUPMAX
      DOUBLE PRECISION MSPMIN,MSPMAX,MSMIN,MSMAX,XIUMIN,XIUMAX
      DOUBLE PRECISION LPPMIN,LPPMAX,LTTMIN,LTTMAX,LUMIN,LUMAX,LDMIN
      DOUBLE PRECISION LDMAX,LTMIN,LTMAX,LBMIN,LBMAX,LLMIN,LLMAX
      DOUBLE PRECISION MSUSYEFF,MMESS,N5,MSREF,D,DMIN

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/SIGMU/SIGMU
      COMMON/MINMAX/MSUSYEFFMIN,MSUSYEFFMAX,MMESSMIN,MMESSMAX,
     . TBMIN,TBMAX,LMIN,LMAX,KMIN,KMAX,ALMIN,ALMAX,XIFMIN,
     . XIFMAX,XISMIN,XISMAX,MUPMIN,MUPMAX,MSPMIN,MSPMAX,
     . MSMIN,MSMAX,XIUMIN,XIUMAX,LPPMIN,LPPMAX,LTTMIN,LTTMAX,LUMIN,
     . LUMAX,LDMIN,LDMAX,LTMIN,LTMAX,LBMIN,LBMAX,LLMIN,LLMAX
      COMMON/STEPS/NTOT,NMSUSYEFF,NMMESS,NTB,NL,NAL,NMUP,NMSP,
     . NXIU,NLPP,NLTT,NLU,NLD,NLT,NLB,NLL,N1,N2
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/PFLAG/PFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG
      COMMON/VFLAG/VFLAG
      COMMON/NMSFLAG/NMSFLAG
      COMMON/CFLAG/CFLAG
      COMMON/NFLAG/NFLAG
      COMMON/M32/M32,CGR,MPL,GRFLAG

*   INITIALIZATION OF THE SUSY PARAMETERS
      DO I=1,NPAR
       PAR(I)=0d0
      ENDDO

*   INITIALIZATION OF THE SCANNING PARAMETERS
      N5=1d99
      SIGMU=1d99
      MSUSYEFFMIN=1d99
      MSUSYEFFMAX=1d99
      MMESSMIN=1d99
      MMESSMAX=1d99
      TBMIN=1d99
      TBMAX=1d99
      LMIN=1d99
      LMAX=1d99
      KMIN=1d99
      KMAX=1d99
      ALMIN=1d99
      ALMAX=1d99
      XIFMIN=1d99
      XIFMAX=1d99
      XISMIN=1d99
      XISMAX=1d99
      MUPMIN=0d0
      MUPMAX=1d99
      MSPMIN=0d0
      MSPMAX=1d99
      MSMIN=1d99
      MSMAX=1d99
      XIUMIN=0d0
      XIUMAX=1d99
      LPPMIN=0d0
      LPPMAX=1d99
      LTTMIN=0d0
      LTTMAX=1d99
      LUMIN=0d0
      LUMAX=1d99
      LDMIN=0d0
      LDMAX=1d99
      LTMIN=0d0
      LTMAX=1d99
      LBMIN=0d0
      LBMAX=1d99
      LLMIN=0d0
      LLMAX=1d99
      DMIN=1d99
      NMSUSYEFF=1
      NMMESS=1
      NTB=1
      NL=1
      NK=0
      NAL=1
      NXIF=0
      NXIS=0
      NMUP=1
      NMSP=1
      NMS=0
      NXIU=1
      NLPP=1
      NLTT=1
      NLU=1
      NLD=1
      NLT=1
      NLB=1
      NLL=1
      CGR=1d0
      MPL=2.4d18

*   DEFAULT VALUES FOR FLAGS
      GMUFLAG=1
      PFLAG=0
      OMGFLAG=0
      NMSFLAG=0
      HFLAG=0
      VFLAG=0
      MOFLAG=7
      GMFLAG=0
      NFLAG=0
      GRFLAG=1
      DO I=1,4
       CFLAG(I)=1
      ENDDO
      CFLAG(5)=0

*   DEFAULT VALUE FOR THE RENSCALE Q2
      Q2=0d0

*   INITIALIZE READ LOOP
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
       END DO
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
       IF(IX.EQ.8) PFLAG=IVAL
       IF(IX.EQ.9) OMGFLAG=IVAL
       IF(IX.EQ.11) GMUFLAG=IVAL
       IF(IX.EQ.12) GMFLAG=IVAL
       IF(IX.EQ.13) NMSFLAG=IVAL
       IF(IX.EQ.14) VFLAG=IVAL
       IF(IX.EQ.15) MOFLAG=IVAL
       IF(IX.EQ.17) CFLAG(1)=IVAL
       IF(IX.EQ.18) CFLAG(2)=IVAL
       IF(IX.EQ.19) CFLAG(3)=IVAL
       IF(IX.EQ.20) CFLAG(4)=IVAL

*   READ SMINPUTS
      ELSEIF(CHBLCK(1:8).EQ.'SMINPUTS')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.2) GF=VAL
       IF(IX.EQ.3) ALSMZ=VAL
       IF(IX.EQ.4) MZ=VAL
       IF(IX.EQ.5) MB=VAL
       IF(IX.EQ.6) MT=VAL
       IF(IX.EQ.7) MTAU=VAL

*   READ MESS PARAMETERS, SIGMU, Q2 AND TANBETA
      ELSEIF(CHBLCK(1:6).EQ.'MINPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0.AND.Q2.EQ.0d0) Q2=VAL**2
       IF(IX.EQ.17) MSUSYEFFMIN=VAL
       IF(IX.EQ.18) MSUSYEFFMAX=VAL
       IF(IX.EQ.27) MMESSMIN=VAL
       IF(IX.EQ.28) MMESSMAX=VAL
       IF(IX.EQ.37) TBMIN=VAL
       IF(IX.EQ.38) TBMAX=VAL
       IF(IX.EQ.4) SIGMU=VAL
       IF(IX.EQ.5) N5=VAL
       IF(IX.EQ.6) CGR=VAL

*   READ EXTPAR
      ELSEIF(CHBLCK(1:6).EQ.'EXTPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0) Q2=VAL**2
       IF(IX.EQ.-2) DMIN=VAL
       IF(IX.EQ.617) LMIN=VAL
       IF(IX.EQ.618) LMAX=VAL
       IF(IX.EQ.627) KMIN=VAL
       IF(IX.EQ.628) KMAX=VAL
       IF(IX.EQ.637) ALMIN=VAL
       IF(IX.EQ.638) ALMAX=VAL
       IF(IX.EQ.667) XIFMIN=VAL
       IF(IX.EQ.668) XIFMAX=VAL
       IF(IX.EQ.677) XISMIN=VAL
       IF(IX.EQ.678) XISMAX=VAL
       IF(IX.EQ.687) MUPMIN=VAL
       IF(IX.EQ.688) MUPMAX=VAL
       IF(IX.EQ.697) MSPMIN=VAL
       IF(IX.EQ.698) MSPMAX=VAL
       IF(IX.EQ.707) MSMIN=VAL
       IF(IX.EQ.708) MSMAX=VAL
       IF(IX.EQ.737) XIUMIN=VAL
       IF(IX.EQ.738) XIUMAX=VAL
       IF(IX.EQ.747) LPPMIN=VAL
       IF(IX.EQ.748) LPPMAX=VAL
       IF(IX.EQ.757) LTTMIN=VAL
       IF(IX.EQ.758) LTTMAX=VAL
       IF(IX.EQ.767) LUMIN=VAL
       IF(IX.EQ.768) LUMAX=VAL
       IF(IX.EQ.777) LDMIN=VAL
       IF(IX.EQ.778) LDMAX=VAL
       IF(IX.EQ.787) LTMIN=VAL
       IF(IX.EQ.788) LTMAX=VAL
       IF(IX.EQ.797) LBMIN=VAL
       IF(IX.EQ.798) LBMAX=VAL
       IF(IX.EQ.807) LLMIN=VAL
       IF(IX.EQ.808) LLMAX=VAL

*   READ STEPS
       ELSEIF(CHBLCK(1:6).EQ.'STEPS')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.7) NFLAG=IVAL
       IF(IX.EQ.19) NMSUSYEFF=IVAL
       IF(IX.EQ.29) NMMESS=IVAL
       IF(IX.EQ.39) NTB=IVAL
       IF(IX.EQ.619) NL=IVAL
       IF(IX.EQ.629) NK=IVAL
       IF(IX.EQ.639) NAL=IVAL
       IF(IX.EQ.669) NXIF=IVAL
       IF(IX.EQ.679) NXIS=IVAL
       IF(IX.EQ.689) NMUP=IVAL
       IF(IX.EQ.699) NMSP=IVAL
       IF(IX.EQ.709) NMS=IVAL
       IF(IX.EQ.739) NXIU=IVAL
       IF(IX.EQ.749) NLPP=IVAL
       IF(IX.EQ.759) NLTT=IVAL
       IF(IX.EQ.769) NLU=IVAL
       IF(IX.EQ.779) NLD=IVAL
       IF(IX.EQ.789) NLT=IVAL
       IF(IX.EQ.799) NLB=IVAL
       IF(IX.EQ.809) NLL=IVAL

      ENDIF

      GOTO 21

*   END OF READING FROM INPUT FILE

*   Check for errors

 29   ERR=0
      DO I=1,5
       IF(CFLAG(I).LT.0 .OR. CFLAG(I).GT.1)THEN
        WRITE(0,1)"CONSTRAINT FLAGS MUST BE IN [0-1]"
        ERR=1
       ENDIF
      ENDDO
      IF(GMUFLAG.LT.0 .OR. GMUFLAG.GT.1)THEN
       WRITE(0,1)"GMUFLAG MUST BE IN [0-1]"
       ERR=1
      ENDIF
      IF(GMFLAG.LT.0 .OR. GMFLAG.GT.4)THEN
       WRITE(0,1)"GMFLAG MUST BE IN [0-4]"
       ERR=1
      ENDIF
      IF(VFLAG.LT.0 .OR. VFLAG.GT.1)THEN
       WRITE(0,1)"VFLAG MUST BE IN [0-1]"
       ERR=1
      ENDIF
      IF(MOFLAG.LT.0 .OR. MOFLAG.GT.7)THEN
       WRITE(0,1)"MOFLAG MUST BE IN [0-7]"
       ERR=1
      ENDIF
      IF(MSUSYEFFMIN.EQ.1d99)THEN
       WRITE(0,1)"MSUSYEFFMIN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(MMESSMIN.EQ.1d99)THEN
       WRITE(0,1)"MMESSMIN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(TBMIN.EQ.1d99)THEN
       WRITE(0,1)"TANBMIN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ELSEIF(TBMIN.LE.0d0 .OR. TBMAX.LE.0d0)THEN
       WRITE(0,1)"TANB MUST BE STRICTLY POSITIVE"
       ERR=1
      ENDIF
      IF(DABS(SIGMU).NE.1d0)THEN
       WRITE(0,1)"SIGMU IS EITHER 1 OR -1 IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(N5.EQ.1d99)THEN
       WRITE(0,1)"N5 MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(N5.LE.0d0 .OR. ANINT(N5).NE.N5)THEN
       WRITE(0,1)"N5 MUST BE A POSITIVE INTEGER IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(CGR.LE.0d0)THEN
       WRITE(0,1)"CGR MUST BE POSITIVE IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(LMIN.EQ.1d99)THEN
       WRITE(0,1)"LMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ELSEIF(LMIN.LE.0d0 .OR. LMAX.LE.0d0)THEN
       WRITE(0,1)"LAMBDA MUST BE STRICTLY POSITIVE"
       ERR=1
      ENDIF
      IF(ALMIN.EQ.1d99 .AND. ALMAX.NE.1d99)THEN
       WRITE(0,1)"ALMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(KMIN.EQ.1d99 .AND. KMAX.NE.1d99)THEN
       WRITE(0,1)"KMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XIFMIN.EQ.1d99 .AND. XIFMAX.NE.1d99)THEN
       WRITE(0,1)"XIFMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XISMIN.EQ.1d99 .AND. XISMAX.NE.1d99)THEN
       WRITE(0,1)"XISMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MSMIN.EQ.1d99 .AND. MSMAX.NE.1d99)THEN
       WRITE(0,1)"MSMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XIFMIN.NE.1d99 .AND. KMIN.NE.1d99)THEN
       WRITE(0,1)"BOTH XIF AND KAPPA CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XISMIN.NE.1d99 .AND. MSMIN.NE.1d99)THEN
       WRITE(0,1)"BOTH XIS AND MS CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.NE.0 .AND. KMIN.NE.1d99)THEN
       WRITE(0,1)"GMFLAG=/=0 => KAPPA CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.NE.0 .AND. ALMIN.NE.1d99)THEN
       WRITE(0,1)"GMFLAG=/=0 => ALAMBDA CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.NE.0 .AND. MSMIN.NE.1d99)THEN
       WRITE(0,1)"GMFLAG=/=0 => MS CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.NE.0 .AND. XIFMIN.NE.1d99)THEN
       WRITE(0,1)"GMFLAG=/=0 => XIF CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.NE.0 .AND. XISMIN.NE.1d99)THEN
       WRITE(0,1)"GMFLAG=/=0 => XIS CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.EQ.0 .AND. DMIN.NE.1d99)THEN
       WRITE(0,1)"GMFLAG=0 => DMIN CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.NE.0 .AND. DMIN.EQ.1d99)THEN
       WRITE(0,1)"GMFLAG=/=0 => DMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.NE.1 .AND. GMFLAG.NE.3 .AND. (XIUMIN.NE.0d0
     ..OR. XIUMAX.NE.1d99))THEN
       WRITE(0,1)"GMFLAG=/=1,3 => XIU CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.NE.2 .AND. GMFLAG.NE.4 .AND. (LPPMIN.NE.0d0
     ..OR. LPPMAX.NE.1d99))THEN
       WRITE(0,1)"GMFLAG=/=2,4 => LPP CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.NE.2 .AND. GMFLAG.NE.4 .AND. (LTTMIN.NE.0d0
     ..OR. LTTMAX.NE.1d99))THEN
       WRITE(0,1)"GMFLAG=/=2,4 => LTT CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.NE.1 .AND. GMFLAG.NE.2 .AND. (LUMIN.NE.0d0
     ..OR. LUMAX.NE.1d99))THEN
       WRITE(0,1)"GMFLAG=/=1,2 => LU CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.NE.1 .AND. GMFLAG.NE.2 .AND. (LDMIN.NE.0d0
     ..OR. LDMAX.NE.1d99))THEN
       WRITE(0,1)"GMFLAG=/=1,2 => LD CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.EQ.0 .AND. (LTMIN.NE.0d0 .OR. LTMAX.NE.1d99))THEN
       WRITE(0,1)"GMFLAG=0 => LT CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.EQ.0 .AND. (LBMIN.NE.0d0 .OR. LBMAX.NE.1d99))THEN
       WRITE(0,1)"GMFLAG=0 => LB CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(GMFLAG.NE.1 .AND. GMFLAG.NE.2 .AND. (LLMIN.NE.0d0
     ..OR. LLMAX.NE.1d99))THEN
       WRITE(0,1)"GMFLAG=/=1,2 => LL CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

*   Set default values

      IF(KMIN.EQ.1d99 .AND. XIFMIN.EQ.1d99)XIFMIN=0d0
      IF(MSMIN.EQ.1d99 .AND. XISMIN.EQ.1d99)XISMIN=0d0
      IF(ALMIN.EQ.1d99)ALMIN=0d0

      IF(MSUSYEFFMAX.EQ.1d99)MSUSYEFFMAX=MSUSYEFFMIN
      IF(MMESSMAX.EQ.1d99)MMESSMAX=MMESSMIN
      IF(TBMAX.EQ.1d99)TBMAX=TBMIN
      IF(LMAX.EQ.1d99)LMAX=LMIN
      IF(KMAX.EQ.1d99)KMAX=KMIN
      IF(ALMAX.EQ.1d99)ALMAX=ALMIN
      IF(XIFMAX.EQ.1d99)XIFMAX=XIFMIN
      IF(XISMAX.EQ.1d99)XISMAX=XISMIN
      IF(MUPMAX.EQ.1d99)MUPMAX=MUPMIN
      IF(MSPMAX.EQ.1d99)MSPMAX=MSPMIN
      IF(MSMAX.EQ.1d99)MSMAX=MSMIN
      IF(XIUMAX.EQ.1d99)XIUMAX=XIUMIN
      IF(LPPMAX.EQ.1d99)LPPMAX=LPPMIN
      IF(LTTMAX.EQ.1d99)LTTMAX=LTTMIN
      IF(LUMAX.EQ.1d99)LUMAX=LUMIN
      IF(LDMAX.EQ.1d99)LDMAX=LDMIN
      IF(LTMAX.EQ.1d99)LTMAX=LTMIN
      IF(LBMAX.EQ.1d99)LBMAX=LBMIN
      IF(LLMAX.EQ.1d99)LLMAX=LLMIN

*   Set MAFLAG

      IF(KMIN.EQ.1d99 .AND. MSMIN.EQ.1d99)MAFLAG=-1
      IF(KMIN.EQ.1d99 .AND. XISMIN.EQ.1d99)MAFLAG=-2
      IF(XIFMIN.EQ.1d99 .AND. MSMIN.EQ.1d99)MAFLAG=-3
      IF(XIFMIN.EQ.1d99 .AND. XISMIN.EQ.1d99)MAFLAG=-4

*   Number of points

      IF(MAFLAG.EQ.-1)THEN
       IF(NK.NE.0)THEN
        WRITE(0,1)"NK CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NMS.NE.0)THEN
        WRITE(0,1)"NMS CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NXIF.EQ.0)NXIF=1
       IF(NXIS.EQ.0)NXIS=1
       N1=NXIF
       N2=NXIS
      ENDIF

      IF(MAFLAG.EQ.-2)THEN
       IF(NK.NE.0)THEN
        WRITE(0,1)"NK CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NXIS.NE.0)THEN
        WRITE(0,1)"NXIS CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NXIF.EQ.0)NXIF=1
       IF(NMS.EQ.0)NMS=1
       N1=NXIF
       N2=NMS
      ENDIF

      IF(MAFLAG.EQ.-3)THEN
       IF(NXIF.NE.0)THEN
        WRITE(0,1)"NXIF CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NMS.NE.0)THEN
        WRITE(0,1)"NMS CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NK.EQ.0)NK=1
       IF(NXIS.EQ.0)NXIS=1
       N1=NK
       N2=NXIS
      ENDIF

      IF(MAFLAG.EQ.-4)THEN
       IF(NXIF.NE.0)THEN
        WRITE(0,1)"NXIF CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NXIS.NE.0)THEN
        WRITE(0,1)"NXIS CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NK.EQ.0)NK=1
       IF(NMS.EQ.0)NMS=1
       N1=NK
       N2=NMS
      ENDIF

      IF(GMFLAG.NE.0)THEN
       IF(ABS(NAL).NE.1)THEN
        WRITE(0,1)"NAL CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NMS.NE.0)THEN
        WRITE(0,1)"NMS CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
      ENDIF

      IF(GMFLAG.NE.1 .AND. GMFLAG.NE.3)THEN
       IF(ABS(NXIU).NE.1)THEN
        WRITE(0,1)"NXIU CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
      ENDIF

      IF(GMFLAG.NE.2 .AND. GMFLAG.NE.4)THEN
       IF(ABS(NLPP).NE.1)THEN
        WRITE(0,1)"NLPP CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(ABS(NLTT).NE.1)THEN
        WRITE(0,1)"NLTT CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
      ENDIF

      IF(GMFLAG.NE.1 .AND. GMFLAG.NE.2)THEN
       IF(ABS(NLU).NE.1)THEN
        WRITE(0,1)"NLU CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(ABS(NLD).NE.1)THEN
        WRITE(0,1)"NLD CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(ABS(NLL).NE.1)THEN
        WRITE(0,1)"NLL CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
      ENDIF

      IF(GMFLAG.EQ.0)THEN
       IF(ABS(NLT).NE.1)THEN
        WRITE(0,1)"NLT CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(ABS(NLB).NE.1)THEN
        WRITE(0,1)"NLB CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
      ENDIF

      IF(NMSUSYEFF.LE.0 .OR. NMMESS.LE.0 .OR. NTB.LE.0
     . .OR. NL.LE.0 .OR. NAL.LE.0 .OR. NMUP.LE.0 .OR. NMSP.LE.0
     . .OR. N1.LE.0 .OR. N2.LE.0 .OR. NXIU.LE.0 .OR. NLPP.LE.0
     . .OR. NLTT.LE.0 .OR. NLU.LE.0 .OR. NLD.LE.0 .OR. NLT.LE.0
     . .OR. NLB.LE.0 .OR. NLL.LE.0)THEN
       WRITE(0,1)"WRONG NUMBER OF POINTS IN BLOCK STEPS"
       ERR=1
      ENDIF

*   Check for Z3 breaking terms

      IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4 .OR.
     . MUPMIN.NE.0d0 .OR. MUPMAX.NE.0d0 .OR. MSPMIN.NE.0d0 .OR. 
     . MSPMAX.NE.0d0 .OR. XIFMIN.NE.0d0 .OR. XIFMAX.NE.0d0 .OR.
     . XISMIN.NE.0d0 .OR. XISMAX.NE.0d0)THEN
       IF(PFLAG.NE.0)THEN
        WRITE(0,1)"HIGGS MASS PRECISION = 1 OR 2 ONLY FOR Z3-NMSSM"
        ERR=1
       ENDIF
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

      IF(MSUSYEFFMIN.EQ.MSUSYEFFMAX .AND. NMSUSYEFF.GT.1)THEN
       WRITE(0,1)"WARNING MSUSYEFFMIN=MSUSYEFFMAX => NMSUSYEFF=1"
       NMSUSYEFF=1
      ENDIF
      IF(MMESSMIN.EQ.MMESSMAX .AND. NMMESS.GT.1)THEN
       WRITE(0,1)"WARNING MMESSMIN=MMESSMAX => NMMESS=1"
       NMMESS=1
      ENDIF
      IF(TBMIN.EQ.TBMAX .AND. NTB.GT.1)THEN
       WRITE(0,1)"WARNING TANBMIN=TANBMAX => NTB=1"
       NTB=1
      ENDIF
      IF(LMIN.EQ.LMAX .AND. NL.GT.1)THEN
       WRITE(0,1)"WARNING LMIN=LMAX => NL=1"
       NL=1
      ENDIF
      IF(ALMIN.EQ.ALMAX .AND. NAL.GT.1)THEN
       WRITE(0,1)"WARNING ALMIN=ALMAX => NAL=1"
       NAL=1
      ENDIF
      IF(MUPMIN.EQ.MUPMAX .AND. NMUP.GT.1)THEN
       WRITE(0,1)"WARNING MUPMIN=MUPMAX => NMUP=1"
       NMUP=1
      ENDIF
      IF(MSPMIN.EQ.MSPMAX .AND. NMSP.GT.1)THEN
       WRITE(0,1)"WARNING MSPMIN=MSPMAX => NMSP=1"
       NMSP=1
      ENDIF
      IF(XIUMIN.EQ.XIUMAX .AND. NXIU.GT.1)THEN
       WRITE(0,1)"WARNING XIUMIN=XIUMAX => NXIU=1"
       NXIU=1
      ENDIF
      IF(LPPMIN.EQ.LPPMAX .AND. NLPP.GT.1)THEN
       WRITE(0,1)"WARNING LPPMIN=LPPMAX => NLPP=1"
       NLPP=1
      ENDIF
      IF(LTTMIN.EQ.LTTMAX .AND. NLTT.GT.1)THEN
       WRITE(0,1)"WARNING LTTMIN=LTTMAX => NLTT=1"
       NLTT=1
      ENDIF
      IF(LUMIN.EQ.LUMAX .AND. NLU.GT.1)THEN
       WRITE(0,1)"WARNING LUMIN=LUMAX => NLU=1"
       NLU=1
      ENDIF
      IF(LDMIN.EQ.LDMAX .AND. NLD.GT.1)THEN
       WRITE(0,1)"WARNING LDMIN=LDMAX => NLD=1"
       NLD=1
      ENDIF
      IF(LTMIN.EQ.LTMAX .AND. NLT.GT.1)THEN
       WRITE(0,1)"WARNING LTMIN=LTMAX => NLT=1"
       NLT=1
      ENDIF
      IF(LBMIN.EQ.LBMAX .AND. NLB.GT.1)THEN
       WRITE(0,1)"WARNING LBMIN=LBMAX => NLB=1"
       NLB=1
      ENDIF
      IF(LLMIN.EQ.LLMAX .AND. NLL.GT.1)THEN
       WRITE(0,1)"WARNING LLMIN=LLMAX => NLL=1"
       NLL=1
      ENDIF

      IF(MSUSYEFFMIN.NE.MSUSYEFFMAX .AND. NMSUSYEFF.EQ.1)THEN
       WRITE(0,10)"WARNING NMSUSYEFF=1 => MSUSYEFFMAX=MSUSYEFFMIN=",
     . MSUSYEFFMIN
       MSUSYEFFMAX=MSUSYEFFMIN
      ENDIF
      IF(MMESSMIN.NE.MMESSMAX .AND. NMMESS.EQ.1)THEN
       WRITE(0,10)"WARNING NMMESS=1 => MMESSMAX=MMESSMIN=",MMESSMIN
       MMESSMAX=MMESSMIN
      ENDIF
      IF(TBMIN.NE.TBMAX .AND. NTB.EQ.1)THEN
       WRITE(0,10)"WARNING NB=1 => TANBMAX=TANBMIN=",TBMIN
       TBMAX=TBMIN
      ENDIF
      IF(LMIN.NE.LMAX .AND. NL.EQ.1)THEN
       WRITE(0,10)"WARNING NL=1 => LMAX=LMIN=",LMIN
       LMAX=LMIN
      ENDIF
      IF(ALMIN.NE.ALMAX .AND. NAL.EQ.1)THEN
       WRITE(0,10)"WARNING NAL=1 => ALMAX=ALMIN=",ALMIN
       ALMAX=ALMIN
      ENDIF
      IF(MUPMIN.NE.MUPMAX .AND. NMUP.EQ.1)THEN
       WRITE(0,10)"WARNING NMUP=1 => MUPMAX=MUPMIN=",MUPMIN
       MUPMAX=MUPMIN
      ENDIF
      IF(MSPMIN.NE.MSPMAX .AND. NMSP.EQ.1)THEN
       WRITE(0,10)"WARNING NMSP=1 => MSPMAX=MSPMIN=",MSPMIN
       MSPMAX=MSPMIN
      ENDIF
      IF(XIUMIN.NE.XIUMAX .AND. NXIU.EQ.1)THEN
       WRITE(0,10)"WARNING NXIU=1 => XIUMAX=XIUMIN=",XIUMIN
       XIUMAX=XIUMIN
      ENDIF
      IF(LPPMIN.NE.LPPMAX .AND. NLPP.EQ.1)THEN
       WRITE(0,10)"WARNING NLPP=1 => LPPMAX=LPPMIN=",LPPMIN
       LPPMAX=LPPMIN
      ENDIF
      IF(LTTMIN.NE.LTTMAX .AND. NLTT.EQ.1)THEN
       WRITE(0,10)"WARNING NLTT=1 => LTTMAX=LTTMIN=",LTTMIN
       LTTMAX=LTTMIN
      ENDIF
      IF(LUMIN.NE.LUMAX .AND. NLU.EQ.1)THEN
       WRITE(0,10)"WARNING NLU=1 => LUMAX=LUMIN=",LUMIN
       LUMAX=LUMIN
      ENDIF
      IF(LDMIN.NE.LDMAX .AND. NLD.EQ.1)THEN
       WRITE(0,10)"WARNING NLD=1 => LDMAX=LDMIN=",LDMIN
       LDMAX=LDMIN
      ENDIF
      IF(LTMIN.NE.LTMAX .AND. NLT.EQ.1)THEN
       WRITE(0,10)"WARNING NLT=1 => LTMAX=LTMIN=",LTMIN
       LTMAX=LTMIN
      ENDIF
      IF(LBMIN.NE.LBMAX .AND. NLB.EQ.1)THEN
       WRITE(0,10)"WARNING NLB=1 => LBMAX=LBMIN=",LBMIN
       LBMAX=LBMIN
      ENDIF
      IF(LLMIN.NE.LLMAX .AND. NLL.EQ.1)THEN
       WRITE(0,10)"WARNING NLL=1 => LLMAX=LLMIN=",LLMIN
       LLMAX=LLMIN
      ENDIF

      IF(MAFLAG.EQ.-1)THEN
       IF(XIFMIN.EQ.XIFMAX .AND. NXIF.GT.1)THEN
        WRITE(0,1)"WARNING XIFMIN=XIFMAX => NXIF=1"
        NXIF=1
        N1=NXIF
       ENDIF
       IF(XISMIN.EQ.XISMAX .AND. NXIS.GT.1)THEN
        WRITE(0,1)"WARNING XISMIN=XISMAX => NXIS=1"
        NXIS=1
        N2=NXIS
       ENDIF
       IF(XIFMIN.NE.XIFMAX .AND. NXIF.EQ.1)THEN
        WRITE(0,10)"WARNING NXIF=1 => XIFMAX=XIFMIN=", XIFMIN
        XIFMAX=XIFMIN
       ENDIF
       IF(XISMIN.NE.XISMAX .AND. NXIS.EQ.1)THEN
        WRITE(0,10)"WARNING NXIS=1 => XISMAX=XISMIN=", XISMIN
        XISMAX=XISMIN
       ENDIF
      ENDIF

      IF(MAFLAG.EQ.-2)THEN
       IF(XIFMIN.EQ.XIFMAX .AND. NXIF.GT.1)THEN
        WRITE(0,1)"WARNING XIFMIN=XIFMAX => NXIF=1"
        NXIF=1
        N1=NXIF
       ENDIF
       IF(MSMIN.EQ.MSMAX .AND. NMS.GT.1)THEN
        WRITE(0,1)"WARNING MSMIN=MSMAX => NMS=1"
        NMS=1
        N2=NMS
       ENDIF
       IF(XIFMIN.NE.XIFMAX .AND. NXIF.EQ.1)THEN
        WRITE(0,10)"WARNING NXIF=1 => XIFMAX=XIFMIN=", XIFMIN
        XIFMAX=XIFMIN
       ENDIF
       IF(MSMIN.NE.MSMAX .AND. NMS.EQ.1)THEN
        WRITE(0,10)"WARNING NMS=1 => MSMAX=MSMIN=", MSMIN
        MSMAX=MSMIN
       ENDIF
      ENDIF

      IF(MAFLAG.EQ.-3)THEN
       IF(KMIN.EQ.KMAX .AND. NK.GT.1)THEN
        WRITE(0,1)"WARNING KMIN=KMAX => NK=1"
        NK=1
        N1=NK
       ENDIF
       IF(XISMIN.EQ.XISMAX .AND. NXIS.GT.1)THEN
        WRITE(0,1)"WARNING XISMIN=XISMAX => NXIS=1"
        NXIS=1
        N2=NXIS
       ENDIF
       IF(KMIN.NE.KMAX .AND. NK.EQ.1)THEN
        WRITE(0,10)"WARNING NK=1 => KMAX=KMIN=", KMIN
        KMAX=KMIN
       ENDIF
       IF(XISMIN.NE.XISMAX .AND. NXIS.EQ.1)THEN
        WRITE(0,10)"WARNING NXIS=1 => XISMAX=XISMIN=", XISMIN
        XISMAX=XISMIN
       ENDIF
      ENDIF

      IF(MAFLAG.EQ.-4)THEN
       IF(KMIN.EQ.KMAX .AND. NK.GT.1)THEN
        WRITE(0,1)"WARNING KMIN=KMAX => NK=1"
        NK=1
        N1=NK
       ENDIF
       IF(MSMIN.EQ.MSMAX .AND. NMS.GT.1)THEN
        WRITE(0,1)"WARNING MSMIN=MSMAX => NMS=1"
        NMS=1
        N2=NMS
       ENDIF
       IF(KMIN.NE.KMAX .AND. NK.EQ.1)THEN
        WRITE(0,10)"WARNING NK=1 => KMAX=KMIN=", KMIN
        KMAX=KMIN
       ENDIF
       IF(MSMIN.NE.MSMAX .AND. NMS.EQ.1)THEN
        WRITE(0,10)"WARNING NMS=1 => MSMAX=MSMIN=", MSMIN
        MSMAX=MSMIN
       ENDIF
      ENDIF

*   total number of points

      NTOT=NMSUSYEFF*NMMESS*NTB*NL*NAL*NMUP*NMSP*NXIU
     .     *NLPP*NLTT*NLU*NLD*NLT*NLB*NLL*N1*N2

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

 1    FORMAT(A)
 10   FORMAT(A,E10.3)

      END


      SUBROUTINE OUTPUT(PAR,PROB,IFAIL)

*********************************************************************
*   Subroutine writing all the results in the the output file.
*********************************************************************

      IMPLICIT NONE

      INTEGER IFAIL,I,NRES,IRES,NBIN,GMFLAG,GRFLAG
      INTEGER NSUSY,NGUT,NMES,IMAX
      PARAMETER (NSUSY=14,NGUT=21,NMES=21,IMAX=200)

      DOUBLE PRECISION RES(IMAX),PAR(*),PROB(*),SIG(5,8)
      DOUBLE PRECISION SMASS(3),PMASS(2),CMASS,SCOMP(3,3),PCOMP(2,2)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,MSREF,D,DMIN
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU,Q2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MHUS,MHDS,MSS
      DOUBLE PRECISION G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      DOUBLE PRECISION LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU
      DOUBLE PRECISION LPPGUT,LTTGUT,LUGUT,LDGUT,LTGUT,LBGUT,LLGUT
      DOUBLE PRECISION SIGMU,MGUT
      DOUBLE PRECISION MHUQ,MHDQ,MSQ,LQ,KQ,ALQ,AKQ,MUQ,NUQ
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
      DOUBLE PRECISION BRKp_Pipnunub,BRKp_Pipnunubmin,BRKp_PipnunubMax,
     .       BRKL_Pi0nunub,BRKL_Pi0nunubmin,BRKL_Pi0nunubMax,
     .       DMK,DMKmin,DMKmax,epsK,epsKmin,epsKmax
      DOUBLE PRECISION eps0,epst0,epst1,epst2,epst3,epsts,epstb
      DOUBLE PRECISION epsY32,epsY31,epsY23,epsY13,epscs,epscb
      DOUBLE PRECISION delmagmu,damumin,damumax,amuthmax,amuthmin
      DOUBLE PRECISION MSUSYEFF,MMESS,N5
      DOUBLE PRECISION M1MES,M2MES,M3MES,ALMES,AKMES,ATMES,ABMES,
     . ATAUMES,AMUMES,MHUMES,MHDMES,MQ3MES,MU3MES,MD3MES,
     . MQMES,MUMES,MDMES,ML3MES,ME3MES,MLMES,MEMES
      DOUBLE PRECISION XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTGUT,HBGUT,HLGUT
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2),toptot
      DOUBLE PRECISION FTSUSY(NSUSY+2),FTGUT(NGUT+2),FTMES(NMES+2)
      DOUBLE PRECISION PX,PA(6),PB(2),PL(7),PK(8)
      DOUBLE PRECISION MHmin,MHmax
      DOUBLE PRECISION OMG,OMGMIN,OMGMAX
      DOUBLE PRECISION Xf,sigmaV,vcsll,vcsbb,x(100),dNdx(100),EMIN
      DOUBLE PRECISION sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      DOUBLE PRECISION M32,CGR,MPL,DELMB,DELML,DEL1
*
      DOUBLE PRECISION chartot2(2),chartot(2),chartot3(2)
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
*
      DOUBLE PRECISION neuttot2(5),neuttot(5),neuttot3(5),neuttotrad(5)
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
*
      DOUBLE PRECISION selltot,selltot2,selltot3,selrtot,selrtot2,
     .         selrtot3,smu1tot,smu1tot2,smu1tot3,smu2tot,smu2tot2,
     .         smu2tot3,stau1tot2,stau2tot,stau2tot2,stau2tot3,
     .         snelltot,snelltot2,snelltot3,snmu1tot,snmu1tot2,
     .         snmu1tot3,sntautot,sntautot2,sntautot3
      DOUBLE PRECISION brsellneute(5),brselrneute(5),brsellcharnue(2),
     .         brselrcharnue(2),brsnellneut(5),brsnellchar(5),
     .         brsmu1neutmu(5),brsmu2neutmu(5),brsmu1charnumu(2),
     .         brsmu2charnumu(2),brsnmu1neut(5),brsnmu1char(5),
     .         brstau1neut(5),brstau2neut(5),brstau1char(2),
     .         brstau2char(2),brstau1hcsn(2),brstau2hcsn(2),
     .         brstau1wsn(2),brstau2wsn(2),brstau2ztau,brstau2H(3),
     .         brstau2A(2),brsntauneut(5),brsntauchar(2),
     .         brsntau1hcstau(2),brsntau1wstau(2)
      DOUBLE PRECISION brsellstau1star,brsellstau1,
     .         brsellstau1nutau,brselrstau1star,brselrstau1,
     .         brselrstau1nutau,brsnestau1star,brsnestau1,
     .         brsnestau1nutau,brsmu1stau1star,brsmu1stau1,
     .         brsmu1stau1nutau,brsmu2stau1star,brsmu2stau1,
     .         brsmu2stau1nutau,brsnmustau1star,brsnmustau1,
     .         brsnmustau1nutau,brstau2stau1star,brstau2stau1,
     .         brstau2stau1nn,brsntaustau1star,brsntaustau1,
     .         brsntaustau1nutau
*
      DOUBLE PRECISION supltot2,suprtot2,sdowltot2,sdowrtot2
      DOUBLE PRECISION brsuplnup(5),brsuplcdow(2),brsuplglui,
     .         brsuprnup(5),brsuprcdow(2),brsuprglui,
     .         brsdowlndow(5),brsdowlchup(2),brsdowlglui,
     .         brsdowrndow(5),brsdowrchup(2),brsdowrglui
*
      DOUBLE PRECISION stoptot(2),stoptot2(2),stoptot3(2),stoptotrad(2)
      DOUBLE PRECISION brst1neutt(5),brst2neutt(5),brst1charb(2),
     .         brst2charb(2),brst1hcsb(2),brst2hcsb(2),brst1wsb(2),
     .         brst2wsb(2),brst1glui,brst2glui,brst2H(3),brst2A(2),
     .         brst2ztop,brgamma,brgammaup,brgammagluino
      DOUBLE PRECISION brstopw(2,5),brstoph(2,5),brststau(2,2),
     .         brstsntau(2,2),brstsmu(2,2),brstsnmu(2),brstsel(2,2),
     .         brstsnel(2),brstbsbst(2,2),brstbbsbt(2,2),
     .         brsttausbnu(2,2),brstelsbnu(2,2),brstupsbdow(2,2),
     .         brst2st1tt,brst2st1startt,brst2st1bb,brst2st1uu,
     .         brst2st1dd,brst2st1ee,brst2st1nunu,brst2st1tautau
*
      DOUBLE PRECISION sbottot(2),sbottot2(2),sbottot3(2)
      DOUBLE PRECISION brsb1neutt(5),brsb2neutt(5),brsb1chart(2),
     .         brsb2chart(2),brsb1hcst(2),brsb2hcst(2),
     .         brsb1glui,brsb2glui,brsb1wst(2),
     .         brsb2wst(2),brsb2H(3),brsb2A(2),brsb2zbot
      DOUBLE PRECISION  brsbstau(2,2),brsbsntau(2,2),brsbsel(2,2),
     .         brsbtstsb(2,2),brsbtbstb(2,2),brsbtaustnu(2,2),
     .         brsbelstnu(2,2),brsbupstdow(2,2),brsbsnel(2),
     .         brsb2sb1bb,brsb2sb1starbb,brsb2sb1tt,
     .         brsb2sb1uu,brsb2sb1dd,brsb2sb1ee,brsb2sb1nunu,
     .         brsb2sb1tautau,brsbsmu(2,2),brsbsnmu(2)
*
      DOUBLE PRECISION gluitot,gluitot2,gluitot3,gluitotrad
      DOUBLE PRECISION brgst1,brgst2,brgsb1,brgsb2,brgsupl,brgsupr,
     .         brgsdownl,brgsdownr,brglnjgluon(5)
      DOUBLE PRECISION brgoup(5),brgoch(5),brgodn(5),brgost(5),
     .         brgotp(5),brgobt(5),brgoud(2),brgocs(2),brgotb(2),
     .         brhcst1b,brwst1b
*
      DOUBLE PRECISION brcharWgra(2),brcharHCgra(2),brneutGAMgra(5),
     .         brneutZgra(5),brneutHgra(5,3),brneutAgra(5,2),
     .         brgluiGLUgra,brselEgra,brserEgra,brsmu1MUgra,
     .         brsmu2MUgra,brstau1TAUgra,brstau2TAUgra,brsneNEgra,
     .         brsnmNMgra,brsntNTgra,brsulUgra,brsurUgra,brsdlDgra,
     .         brsdrDgra,brst1Tgra,brst2Tgra,brsb1Bgra,brsb2Bgra
*
      COMMON/CHARGINO_WIDTH/chartot2,chartot,chartot3
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
*
      COMMON/NEUTRALINO_WIDTH/neuttot2,neuttot,neuttot3,neuttotrad
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
*
      COMMON/SLEPTON_WIDTH/selltot,selltot2,selltot3,selrtot,
     .         selrtot2,selrtot3,smu1tot,smu1tot2,smu1tot3,smu2tot,
     .         smu2tot2,smu2tot3,stau1tot2,stau2tot,stau2tot2,
     .         stau2tot3,snelltot,snelltot2,snelltot3,snmu1tot,
     .         snmu1tot2,snmu1tot3,sntautot2,sntautot3,sntautot
      COMMON/SLEPTON_BR_2BD/brsellneute,brselrneute,brsellcharnue,
     .         brselrcharnue,brsnellneut,brsnellchar,brsmu1neutmu,
     .         brsmu2neutmu,brsmu1charnumu,brsmu2charnumu,brsnmu1neut,
     .         brsnmu1char,brstau1neut,brstau2neut,brstau1char,
     .         brstau2char,brstau1hcsn,brstau2hcsn,brstau1wsn,
     .         brstau2wsn,brstau2ztau,brstau2H,brstau2A,brsntauneut,
     .         brsntauchar,brsntau1hcstau,brsntau1wstau
      COMMON/SLEPTON_BR_3BD/brsellstau1star,brsellstau1,
     .         brsellstau1nutau,brselrstau1star,brselrstau1,
     .         brselrstau1nutau,brsnestau1star,brsnestau1,
     .         brsnestau1nutau,brsmu1stau1star,brsmu1stau1,
     .         brsmu1stau1nutau,brsmu2stau1star,brsmu2stau1,
     .         brsmu2stau1nutau,brsnmustau1star,brsnmustau1,
     .         brsnmustau1nutau,brstau2stau1star,brstau2stau1,
     .         brstau2stau1nn,brsntaustau1star,brsntaustau1,
     .         brsntaustau1nutau
*
      COMMON/SQUARK_WIDTH/supltot2,suprtot2,sdowltot2,sdowrtot2
      COMMON/SQUARK_BR_2BD/brsuplnup,brsuplcdow,brsuplglui,
     .         brsuprnup,brsuprcdow,brsuprglui,
     .         brsdowlndow,brsdowlchup,brsdowlglui,
     .         brsdowrndow,brsdowrchup,brsdowrglui
*
      COMMON/STOP_WIDTH/stoptot,stoptot2,stoptot3,stoptotrad
      COMMON/STOP_BR_2BD/brst1neutt,brst2neutt,brst1charb,
     .         brst2charb,brst1hcsb,brst2hcsb,brst1wsb,
     .         brst2wsb,brst1glui,brst2glui,brst2H,brst2A,
     .         brst2ztop,brgamma,brgammaup,brgammagluino
      COMMON/STOP_BR_3BD/brstopw,brstoph,brststau,
     .         brstsntau,brstsmu,brstsnmu,brstsel,
     .         brstsnel,brstbsbst,brstbbsbt,
     .         brsttausbnu,brstelsbnu,brstupsbdow,
     .         brst2st1tt,brst2st1startt,brst2st1bb,brst2st1uu,
     .         brst2st1dd,brst2st1ee,brst2st1nunu,brst2st1tautau
*
      COMMON/SBOTTOM_WIDTH/sbottot,sbottot2,sbottot3
      COMMON/SBOTTOM_BR_2BD/brsb1neutt,brsb2neutt,brsb1chart,
     .         brsb2chart,brsb1hcst,brsb2hcst,
     .         brsb1glui,brsb2glui,brsb1wst,
     .         brsb2wst,brsb2H,brsb2A,brsb2zbot
      COMMON/SBOTTOM_BR_3BD/brsbstau,brsbsntau,brsbsel,
     .         brsbtstsb,brsbtbstb,brsbtaustnu,
     .         brsbelstnu,brsbupstdow,brsbsnel,
     .         brsb2sb1bb,brsb2sb1starbb,brsb2sb1tt,
     .         brsb2sb1uu,brsb2sb1dd,brsb2sb1ee,brsb2sb1nunu,
     .         brsb2sb1tautau,brsbsmu,brsbsnmu
*
      COMMON/GLUINO_WIDTH/gluitot,gluitot2,gluitot3,gluitotrad
      COMMON/GLUINO_BR_2BD/brgst1,brgst2,brgsb1,brgsb2,brgsupl,brgsupr,
     .         brgsdownl,brgsdownr,brglnjgluon
      COMMON/GLUINO_BR_3BD/brgoup,brgoch,brgodn,brgost,brgotp,
     .         brgobt,brgoud,brgocs,brgotb,brhcst1b,brwst1b
*
      COMMON/GRAVITINO/brcharWgra,brcharHCgra,brneutGAMgra,
     .         brneutZgra,brneutHgra,brneutAgra,
     .         brgluiGLUgra,brselEgra,brserEgra,brsmu1MUgra,
     .         brsmu2MUgra,brstau1TAUgra,brstau2TAUgra,brsneNEgra,
     .         brsnmNMgra,brsntNTgra,brsulUgra,brsurUgra,brsdlDgra,
     .         brsdrDgra,brst1Tgra,brst2Tgra,brsb1Bgra,brsb2Bgra
*
      COMMON/SIGMU/SIGMU
      COMMON/INPPAR/ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/MESEXT/XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/SOFTMES/M1MES,M2MES,M3MES,ALMES,AKMES,ATMES,ABMES,
     . ATAUMES,AMUMES,MHUMES,MHDMES,MQ3MES,MU3MES,MD3MES,
     . MQMES,MUMES,MDMES,ML3MES,ME3MES,MLMES,MEMES
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     . BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     . BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     . BRSUSY,WIDTH
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     . HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     . HCBRSUSY,HCWIDTH
      COMMON/BRSG/BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     .      DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     .      BRBtaunumax,BRBtaunumin
      COMMON/FLAV2/BRBSll,BRBSllmin,BRBSllMax,
     .      BRBShll,BRBShllmin,BRBShllMax,
     .      BRDG,BRDGmin,BRDGmax,
     .      BRBdMUMU,BRBdMUMUmin,BRBdMUMUmax,
     .      BRBXsnunu,BRBXsnunumin,BRBXsnunumax,
     .      BRBpKpnunu,BRBpKpnunumin,BRBpKpnunuMax,
     .      BRBKsnunu,BRBKsnunumin,BRBKsnunuMax,
     .      RD_taul,RD_taulmin,RD_taulmax,
     .      RDs_taul,RDs_taulmin,RDs_taulmax
      COMMON/FLAV3/BRKp_Pipnunub,BRKp_Pipnunubmin,BRKp_PipnunubMax,
     .             BRKL_Pi0nunub,BRKL_Pi0nunubmin,BRKL_Pi0nunubMax,
     .             DMK,DMKmin,DMKmax,epsK,epsKmin,epsKmax
      COMMON/EPSCOUP/eps0,epst0,epst1,epst2,epst3,epsts,epstb,
     .               epsY32,epsY31,epsY23,epsY13,epscs,epscb
      COMMON/MAGMU/delmagmu,damumin,damumax,amuthmax,amuthmin
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/topwidth/toptot
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     . MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     . CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SUSYMH/MHUS,MHDS,MSS
      COMMON/QMHIGGS/MHUQ,MHDQ,MSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/MESCOUP/G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      COMMON/MESGUT/LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTGUT,HBGUT,HLGUT
      COMMON/GUTMES/LPPGUT,LTTGUT,LUGUT,LDGUT,LTGUT,LBGUT,LLGUT
      COMMON/MGUT/MGUT
      COMMON/FINETUN/FTSUSY,FTGUT,FTMES
      COMMON/EFFCOUP/PX,PA,PB,PL,PK
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/LHCSIG/SIG
      COMMON/HIGGSFIT/MHmin,MHmax
      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,vcsll,vcsbb,
     .      x,dNdx,EMIN,NBIN
      COMMON/MICROMG2/sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      COMMON/M32/M32,CGR,MPL,GRFLAG

      IF(IFAIL.NE.0)RETURN

      IRES=2
      NRES=26+IRES

      RES(1)=PAR(1)
      RES(2)=PAR(3)

      DO I=1,2
       RES(IRES-1+2*I)=SMASS(I)
       RES(IRES+2*I)=SCOMP(I,3)**2
      ENDDO
      RES(IRES+5)=PMASS(1)
      RES(IRES+6)=PCOMP(1,2)**2
      RES(IRES+7)=CMASS
      DO I=1,3
       RES(IRES+4+4*I)=DABS(MNEU(I))
       RES(IRES+5+4*I)=NEU(I,1)**2
       RES(IRES+6+4*I)=NEU(I,3)**2+NEU(I,4)**2
       RES(IRES+7+4*I)=NEU(I,5)**2
      ENDDO
      RES(IRES+20)=DABS(MCHA(1))
      RES(IRES+21)=MGL
      RES(IRES+22)=MIN(MUL,MUR,MDL,MDR)
      RES(IRES+23)=MST1
      RES(IRES+24)=MSL1
      RES(IRES+25)=MSNT
      RES(IRES+26)=FTMES(NMES+1)

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
      INTEGER OMGFLAG,MAFLAG,MOFLAG,GMFLAG

      DOUBLE PRECISION MSUSYEFFN,MSUSYEFFNN,MMESSN,MMESSNN,TBN,TBNN
      DOUBLE PRECISION LN,LNN,KN,KNN,ALN,ALNN,XIFN,XIFNN,XISN,XISNN
      DOUBLE PRECISION MUPN,MUPNN,MSPN,MSPNN,MSN,MSNN,XIUN,XIUNN
      DOUBLE PRECISION LPPN,LPPNN,LTTN,LTTNN,LUN,LUNN,LDN,LDNN
      DOUBLE PRECISION LTN,LTNN,LBN,LBNN,LLN,LLNN,MUN,MUNN,DEV
      DOUBLE PRECISION MSREF,D,DMIN

      COMMON/BOUNDS/MSUSYEFFN,MSUSYEFFNN,MMESSN,MMESSNN,TBN,TBNN,
     . LN,LNN,KN,KNN,ALN,ALNN,XIFN,XIFNN,XISN,XISNN,MUPN,MUPNN,
     . MSPN,MSPNN,MSN,MSNN,XIUN,XIUNN,LPPN,LPPNN,LTTN,LTTNN,LUN,LUNN,
     . LDN,LDNN,LTN,LTNN,LBN,LBNN,LLN,LLNN,MUN,MUNN
      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/CFLAG/CFLAG

      WRITE(0,10)"NMSSMTools scan info                     "
      WRITE(0,10)"Version number: 5.6.2                    "
      WRITE(0,*)
      WRITE(0,*)
      WRITE(0,10)"Number of points:                        "
      WRITE(0,*)
      WRITE(0,10)"  scanned                                ",NTOT
      WRITE(0,10)"  no electroweak symmetry breaking       ",NFAIL(20)
      S=0
      DO I=1,7
       S=S+NFAIL(I)
      ENDDO
      WRITE(0,10)"  with mh1^2 or ma1^2 or mhc^2 < 0      ",S
      WRITE(0,10)"  with m_sfermion^2 < 0                 ",NFAIL(8)
      WRITE(0,10)"  violating constraints                 ",NFAIL(10)
      IF(GMFLAG.NE.0)THEN
       WRITE(0,10)"  MSMES=/=MSREF                         ",NFAIL(21)
       WRITE(0,10)"  MSMES=/=MSREF + violating constraints ",NFAIL(22)
      ENDIF
      S=0
      DO I=11,15
       S=S+NFAIL(I)
      ENDDO
      WRITE(0,10)"  RGE integration problem               ",S
      S=0
      DO I=16,19
       S=S+NFAIL(I)
      ENDDO
      WRITE(0,10)"  convergence problem                   ",S
      WRITE(0,*)
      WRITE(0,10)"Remaining good points                   ",TOT
      WRITE(0,*)
      WRITE(0,*)

      IF(OMGFLAG.EQ.0 .AND. CFLAG(1).EQ.0 .AND. CFLAG(2).EQ.0 .AND.
     .   CFLAG(3).EQ.0 .AND. CFLAG(4).EQ.0 .AND. GMUFLAG.EQ.0)THEN
       WRITE(0,20)"Contraints taken into account: none               "
      ELSE
       WRITE(0,20)"Contraints taken into account:                    "
       IF(OMGFLAG.GT.0)
     .  WRITE(0,20)" - Relic density from Planck +/- 10% [0.107,0.131]"
       IF(OMGFLAG.LT.0)
     .  WRITE(0,20)" - Relic density from Planck upper bound < 0.131  "
       IF(CFLAG(1).EQ.1)
     .  WRITE(0,20)" - Landau poles and false minima                  "
       IF(CFLAG(2).EQ.1)
     .  WRITE(0,20)" - LEP/Tevatron Higgs+sparticle                   "
       IF(CFLAG(3).EQ.1)
     .  WRITE(0,20)" - LHC Higgs                                      "
       IF(CFLAG(4).EQ.1)
     .  WRITE(0,20)" - Upsilon, B and K decays                        "
       IF(GMUFLAG.EQ.1)
     .  WRITE(0,20)" - (g-2)_muon                                     "
      ENDIF

      IF(TOT.GT.0)THEN

       WRITE(0,*)
       WRITE(0,*)
       WRITE(0,20)"Parameter ranges for good points:                 "
       WRITE(0,*)
       WRITE(0,30)" MSUSYEFF: ",MSUSYEFFN,MSUSYEFFNN,
     .                          DEV(MSUSYEFFN,MSUSYEFFNN)
       WRITE(0,30)" MMESS: ",MMESSN,MMESSNN,DEV(MMESSN,MMESSNN)
       WRITE(0,30)" TANB: ",TBN,TBNN,DEV(TBN,TBNN)
       WRITE(0,30)" LAMBDA: ",LN,LNN,DEV(LN,LNN)
       WRITE(0,30)" ALAMBDA: ",ALN,ALNN,DEV(ALN,ALNN)
       IF(GMFLAG.NE.0)
     .  WRITE(0,40)"(ALAMBDA is not an input parameter)"
       WRITE(0,30)"MUEFF: ",MUN,MUNN,DEV(MUN,MUNN)
       WRITE(0,40)"(MUEFF is not an input parameter)"
       IF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4)THEN
        WRITE(0,30)" KAPPA: ",KN,KNN,DEV(KN,KNN)
        WRITE(0,30)" XIF: ",XIFN,XIFNN,DEV(XIFN,XIFNN)
        WRITE(0,40)"(XIF is not an input parameter)"
       ELSEIF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        WRITE(0,30)" KAPPA: ",KN,KNN,DEV(KN,KNN)
        WRITE(0,40)"(KAPPA is not an input parameter)"
        IF(XIFN.NE.0d0 .OR. XIFNN.NE.0d0)
     .   WRITE(0,30)" XIF: ",XIFN,XIFNN,DEV(XIFN,XIFNN)
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        IF(XISN.NE.0d0 .OR. XISNN.NE.0d0)
     .   WRITE(0,30)" XIS: ",XISN,XISNN,DEV(XISN,XISNN)
        WRITE(0,30)" MS: ",MSN,MSNN,DEV(MSN,MSNN)
        WRITE(0,40)"(MS is not an input parameter)"
       ELSEIF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)THEN
        WRITE(0,30)" MS: ",MSN,MSNN,DEV(MSN,MSNN)
        WRITE(0,30)" XIS: ",XISN,XISNN,DEV(XISN,XISNN)
        WRITE(0,40)"(XIS is not an input parameter)"
       ENDIF
       IF(MUPN.NE.0d0 .OR. MUPNN.NE.0d0)
     . WRITE(0,30)" MUP: ",MUPN,MUPNN,DEV(MUPN,MUPNN)
       IF(MSPN.NE.0d0 .OR. MSPNN.NE.0d0)
     . WRITE(0,30)" MSP: ",MSPN,MSPNN,DEV(MSPN,MSPNN)
       IF((GMFLAG.EQ.1 .OR. GMFLAG.EQ.3) .AND. (XIUN.NE.0d0
     . .OR. XIUNN.NE.0d0))
     . WRITE(0,30)" XIU: ",XIUN,XIUNN,DEV(XIUN,XIUNN)
       IF((GMFLAG.EQ.2 .OR. GMFLAG.EQ.4) .AND. (LPPN.NE.0d0
     . .OR. LPPNN.NE.0d0))
     . WRITE(0,30)" LPP: ",LPPN,LPPNN,DEV(LPPN,LPPNN)
       IF((GMFLAG.EQ.2 .OR. GMFLAG.EQ.4) .AND. (LTTN.NE.0d0
     . .OR. LTTNN.NE.0d0))
     . WRITE(0,30)" LTT: ",LTTN,LTTNN,DEV(LTTN,LTTNN)
       IF((GMFLAG.EQ.1 .OR. GMFLAG.EQ.2) .AND. (LUN.NE.0d0
     . .OR. LUNN.NE.0d0))
     . WRITE(0,30)" LU: ",LUN,LUNN,DEV(LUN,LUNN)
       IF((GMFLAG.EQ.1 .OR. GMFLAG.EQ.2) .AND. (LDN.NE.0d0
     . .OR. LDNN.NE.0d0))
     . WRITE(0,30)" LD: ",LDN,LDNN,DEV(LDN,LDNN)
       IF(GMFLAG.NE.0 .AND. (LTN.NE.0d0 .OR. LTNN.NE.0d0))
     . WRITE(0,30)" LT: ",LTN,LTNN,DEV(LTN,LTNN)
       IF(GMFLAG.NE.0 .AND. (LBN.NE.0d0 .OR. LBNN.NE.0d0))
     . WRITE(0,30)" LB: ",LBN,LBNN,DEV(LBN,LBNN)
       IF((GMFLAG.EQ.1 .OR. GMFLAG.EQ.2) .AND. (LLN.NE.0d0
     . .OR. LLNN.NE.0d0))
     . WRITE(0,30)" LL: ",LLN,LLNN,DEV(LLN,LLNN)

      ENDIF

 10   FORMAT(A40,I10)
 20   FORMAT(A50)
 30   FORMAT(A11,3E15.4)
 40   FORMAT(A36)

      END
