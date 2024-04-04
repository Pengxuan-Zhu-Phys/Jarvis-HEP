      PROGRAM MAIN

*  Program to compute the NMSSM Higgs masses, couplings and
*  Branching ratios, with experimental and theoretical constraints
*
*   On input:
*
*      tan(beta) at the scale MZ, lambda at the scale Q2
*      m0, M1/2, A0 at the scale MGUT
*      non-universal parameters are allowed in the Higgs/gaugino sector
*
*      The input file contains the starting point
*      as well as the total number of MCMC steps
*
*   On output:
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
*  ERRORS: IFAIL = 0..17
*
*  IFAIL = 0         OK
*          1,3,5,7   m_h1^2 < 0
*          2,3,6,7   m_a1^2 < 0
*          4,5,6,7   m_h+^2 < 0
*          8         m_sfermion^2 < 0
*          9         mu = 0
*          10        Violation of phenomenological constraint(s)
*          11,12,13  Problem in integration of RGEs
*          14,15,16  Convergence problem
*          17        No electroweak symmetry breaking
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
*      PROB(26) =/= 0  lightest neutralino is not LSP or < 511 keV
*      PROB(27) =/= 0  Landau Pole in l, k, ht, hb below MGUT
*      PROB(28) =/= 0  unphysical global minimum
*      PROB(29) =/= 0  Higgs soft masses >> Msusy
*      PROB(30) =/= 0  excluded by DM relic density (checked only if OMGFLAG=/=0)
*      PROB(31) =/= 0  excluded by DM SI WIMP-nucleon xs (checked if |OMGFLAG|=2 or 4)
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
*      PROB(50) =/= 0  excluded by sparticle searches at the LHC (not called by default)
*      PROB(51) =/= 0: excluded by H/A->tautau (ATLAS+CMS)
*      PROB(52) =/= 0: excluded by H->AA->4leptons/2lept.+2b (ATLAS+CMS)
*      PROB(53) =/= 0: excluded by ggF->H/A->gamgam (ATLAS)
*      PROB(55) =/= 0: b -> d gamma more than 2 sigma away
*      PROB(56) =/= 0: B_d -> mu+ mu- more than 2 sigma away
*      PROB(57) =/= 0: b -> s nu nubar more than 2 sigma away
*      PROB(58) =/= 0: b -> c tau nu more than 2 sigma away (as SM)
*      PROB(59) =/= 0: K -> pi nu nubar more than 2 sigma away
*      PROB(60) =/= 0: DMK / epsK more than 2 sigma away
*      PROB(61) =/= 0  excluded by DM SD WIMP-neutron xs (checked if |OMGFLAG|=2 or 4)
*      PROB(62) =/= 0  excluded by DM SD WIMP-proton xs (checked if |OMGFLAG|=2 or 4)
*      PROB(63) =/= 0: excluded by H->AA->4gammas (ATLAS+CMS)
*      PROB(64) =/= 0: excluded by trilepton searches for charg(neutral)inos (CMS)
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
      PARAMETER (NFL=17,NPROB=82,NPAR=25)
      INTEGER NFAIL(NFL),IFAIL,I,TOT,ITOT,NTOT,IDUM,ISEED
      INTEGER ITER,ITERMU,ITRY,MUFLAG,SFFLAG,Q2FIX,IM
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,MSFLAG
      INTEGER AKFLAG,ALFLAG,OMGFLAG,MAFLAG,MOFLAG,GMUFLAG
      INTEGER HFLAG,NMSFLAG,UNCERTFLAG,GRFLAG,CFLAG(5)
      INTEGER TOTMIN,TOTMAX,NMAX,IP

      DOUBLE PRECISION PAR(NPAR),PROB(NPROB),CHECK,GUTEST
      DOUBLE PRECISION DETM,R,K,MUM,MDM,MSM,MUT,MDT,MST,Q2,Q2MIN
      DOUBLE PRECISION M0,M12,A0,SIGMU,MUFAIL,MUSTEP,MUINIT
      DOUBLE PRECISION M0N,M0NN,M12N,M12NN,TBN,TBNN,A0N,A0NN
      DOUBLE PRECISION M1N,M1NN,M2N,M2NN,M3N,M3NN,MHDN,MHDNN
      DOUBLE PRECISION MHUN,MHUNN,MSN,MSNN,LN,LNN,KN,KNN,ALN,ALNN
      DOUBLE PRECISION AKN,AKNN,MUN,MUNN,XIFN,XIFNN,XISN,XISNN
      DOUBLE PRECISION MUPN,MUPNN,MSPN,MSPNN,M3HN,M3HNN
      DOUBLE PRECISION M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT
      DOUBLE PRECISION ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT
      DOUBLE PRECISION MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      DOUBLE PRECISION MU3GUT,MD3GUT,GAU
      DOUBLE PRECISION M1INP,M2INP,M3INP,MHDINP,MHUINP,ALINP,AKINP
      DOUBLE PRECISION XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      DOUBLE PRECISION XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      DOUBLE PRECISION M0CEN,M0DEV,M12CEN,M12DEV,TBCEN,TBDEV,A0CEN
      DOUBLE PRECISION A0DEV,M1CEN,M1DEV,M2CEN,M2DEV,M3CEN,M3DEV
      DOUBLE PRECISION MHDCEN,MHDDEV,MHUCEN,MHUDEV,MSCEN,MSDEV,LCEN
      DOUBLE PRECISION LDEV,KCEN,KDEV,ALCEN,ALDEV,AKCEN,AKDEV,MUCEN
      DOUBLE PRECISION MUDEV,XIFCEN,XIFDEV,XISCEN,XISDEV,MUPCEN
      DOUBLE PRECISION MUPDEV,MSPCEN,MSPDEV,M3HCEN,M3HDEV,XCEN,XDEV,X
      DOUBLE PRECISION M0MIN,M12MIN,TBMIN,A0MIN,M1MIN,M2MIN,M3MIN
      DOUBLE PRECISION MHDMIN,MHUMIN,MSMIN,LMIN,KMIN,ALMIN,AKMIN
      DOUBLE PRECISION MUMIN,XIFMIN,XISMIN,MUPMIN,MSPMIN,M3HMIN
      DOUBLE PRECISION Q2MEM,DELMEM,XIFMEM,XISMEM,MUMEM,MDMEM,MSMEM
      DOUBLE PRECISION PARMEM(25),M32,CGR,MPL,DELMB,DELML,DEL1

      COMMON/NMSFLAG/NMSFLAG
      COMMON/SIGMU/SIGMU
      COMMON/SOFTGUT/M0,M12,A0
      COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     . ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT,MU3GUT,
     . MD3GUT,MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      COMMON/GUTEXT/XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/DETM/DETM
      COMMON/MUFAIL/MUFAIL
      COMMON/BOUNDS/M0N,M0NN,M12N,M12NN,TBN,TBNN,A0N,A0NN,
     . M1N,M1NN,M2N,M2NN,M3N,M3NN,MHDN,MHDNN,MHUN,MHUNN,
     . MSN,MSNN,LN,LNN,KN,KNN,ALN,ALNN,AKN,AKNN,MUN,MUNN,
     . XIFN,XIFNN,XISN,XISNN,MUPN,MUPNN,MSPN,MSPNN,M3HN,M3HNN
      COMMON/INPPAR/M1INP,M2INP,M3INP,MHDINP,MHUINP,ALINP,AKINP,
     . XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      COMMON/MCMCPAR/M0CEN,M0DEV,M12CEN,M12DEV,TBCEN,TBDEV,A0CEN,
     . A0DEV,M1CEN,M1DEV,M2CEN,M2DEV,M3CEN,M3DEV,MHDCEN,MHDDEV,
     . MHUCEN,MHUDEV,MSCEN,MSDEV,LCEN,LDEV,KCEN,KDEV,ALCEN,ALDEV,
     . AKCEN,AKDEV,MUCEN,MUDEV,XIFCEN,XIFDEV,XISCEN,XISDEV,MUPCEN,
     . MUPDEV,MSPCEN,MSPDEV,M3HCEN,M3HDEV,XCEN,XDEV,X,
     . M0MIN,M12MIN,TBMIN,A0MIN,M1MIN,M2MIN,M3MIN,
     . MHDMIN,MHUMIN,MSMIN,LMIN,KMIN,ALMIN,AKMIN,
     . MUMIN,XIFMIN,XISMIN,MUPMIN,MSPMIN,M3HMIN
      COMMON/STEPS/NTOT,IDUM,TOTMIN,TOTMAX,NMAX
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . MSFLAG,AKFLAG,ALFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/MSAVE/MUM,MDM,MSM,MUT,MDT,MST,IM
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/UNCERTFLAG/UNCERTFLAG
      COMMON/CFLAG/CFLAG
      COMMON/M32/M32,CGR,MPL,GRFLAG
      COMMON/MEM/Q2MEM,DELMEM,XIFMEM,XISMEM,MUMEM,MDMEM,MSMEM,PARMEM

*   Initialization

      CALL INITIALIZE()
      DO I=1,NFL
       NFAIL(I)=0
      ENDDO
      TOT=0
      IP=0

*   Reading of the input parameters

      CALL INPUT(PAR,NPAR)
      ISEED=IDUM

*   Initialization of the range of parameters that has passed all tests

      M0N=1d99
      M0NN=-1d99
      M12N=1d99
      M12NN=-1d99
      A0N=1d99
      A0NN=-1d99
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
      ALN=1d99
      ALNN=-1d99
      AKN=1d99
      AKNN=-1d99
      MUN=1d99
      MUNN=-1d99
      MHDN=1d99
      MHDNN=-1d99
      MHUN=1d99
      MHUNN=-1d99
      MSN=1d99
      MSNN=-1d99
      XIFN=1d99
      XIFNN=-1d99
      XISN=1d99
      XISNN=-1d99
      MUPN=1d99
      MUPNN=-1d99
      MSPN=1d99
      MSPNN=-1d99
      M3HN=1d99
      M3HNN=-1d99

*   Beginning of the scan

      DO ITOT=1,NTOT

 14   IF(ITOT.EQ.1)THEN

       M0=M0CEN
       M12=M12CEN
       A0=A0CEN
       PAR(3)=TBCEN
       PAR(1)=LCEN
       IF(ALFLAG.EQ.0)THEN
        ALINP=A0
       ELSE
        ALINP=ALCEN
       ENDIF
       IF(AKFLAG.EQ.0)THEN
        AKINP=A0
       ELSE
        AKINP=AKCEN
       ENDIF
       IF(M1FLAG.EQ.0)THEN
        M1INP=M12
       ELSE
        M1INP=M1CEN
       ENDIF
       IF(M2FLAG.EQ.0)THEN
        M2INP=M12
       ELSE
        M2INP=M2CEN
       ENDIF
       IF(M3FLAG.EQ.0)THEN
        M3INP=M12
       ELSE
        M3INP=M3CEN
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        XIFINP=XIFCEN
       ELSE
        PAR(2)=KCEN
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        XISINP=XISCEN
       ELSE
        MSINP=MSCEN
       ENDIF
       IF(MAFLAG.EQ.-5)THEN
        PAR(4)=MUCEN
	SIGMU=DABS(PAR(4))/PAR(4)
        PAR(2)=KCEN
        XIFINP=XIFCEN
        XISINP=XISCEN
       ELSE
        IF(MHDFLAG.EQ.0)THEN
         MHDINP=M0**2
        ELSE
         MHDINP=MHDCEN
        ENDIF
        IF(MHUFLAG.EQ.0)THEN
         MHUINP=M0**2
        ELSE
         MHUINP=MHUCEN
        ENDIF
       ENDIF
       MUPINP=MUPCEN
       MSPINP=MSPCEN
       M3HINP=M3HCEN

      ELSE

       IF(M0DEV.EQ.0d0)THEN
        M0=M0CEN
       ELSE
 101    M0=M0CEN+MAX(DABS(M0CEN),M0MIN)*M0DEV*GAU(IDUM)
        IF(M0.LE.0d0)GOTO 101
       ENDIF

       IF(M12DEV.EQ.0d0)THEN
        M12=M12CEN
       ELSE
 102    M12=M12CEN+MAX(DABS(M12CEN),M12MIN)*M12DEV*GAU(IDUM)
        IF(M12.LE.0d0)GOTO 102
       ENDIF

       IF(A0DEV.EQ.0d0)THEN
        A0=A0CEN
       ELSE
        A0=A0CEN+MAX(DABS(A0CEN),A0MIN)*A0DEV*GAU(IDUM)
       ENDIF

       IF(TBDEV.EQ.0d0)THEN
        PAR(3)=TBCEN
       ELSE
 103    PAR(3)=TBCEN+MAX(DABS(TBCEN),TBMIN)*TBDEV*GAU(IDUM)
        IF(PAR(3).LT.1d0)GOTO 103
       ENDIF

       IF(LDEV.EQ.0d0)THEN
        PAR(1)=LCEN
       ELSE
 104    PAR(1)=LCEN+MAX(DABS(LCEN),LMIN)*LDEV*GAU(IDUM)
        IF(PAR(1).LE.0d0)GOTO 104
       ENDIF

       IF(ALFLAG.EQ.0)THEN
        ALINP=A0
       ELSE
        IF(ALDEV.EQ.0d0)THEN
         ALINP=ALCEN
        ELSE
         ALINP=ALCEN+MAX(DABS(ALCEN),ALMIN)*ALDEV*GAU(IDUM)
        ENDIF
       ENDIF

       IF(AKFLAG.EQ.0)THEN
        AKINP=A0
       ELSE
        IF(AKDEV.EQ.0d0)THEN
         AKINP=AKCEN
        ELSE
         AKINP=AKCEN+MAX(DABS(AKCEN),AKMIN)*AKDEV*GAU(IDUM)
        ENDIF
       ENDIF

       IF(M1FLAG.EQ.0)THEN
        M1INP=M12
       ELSE
        IF(M1DEV.EQ.0d0)THEN
         M1INP=M1CEN
        ELSE
         M1INP=M1CEN+MAX(DABS(M1CEN),M1MIN)*M1DEV*GAU(IDUM)
        ENDIF
       ENDIF

       IF(M2FLAG.EQ.0)THEN
        M2INP=M12
       ELSE
        IF(M2DEV.EQ.0d0)THEN
         M2INP=M2CEN
        ELSE
         M2INP=M2CEN+MAX(DABS(M2CEN),M2MIN)*M2DEV*GAU(IDUM)
        ENDIF
       ENDIF

       IF(M3FLAG.EQ.0)THEN
        M3INP=M12
       ELSE
        IF(M3DEV.EQ.0d0)THEN
         M3INP=M3CEN
        ELSE
         M3INP=M3CEN+MAX(DABS(M3CEN),M3MIN)*M3DEV*GAU(IDUM)
        ENDIF
       ENDIF

       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        IF(XIFDEV.EQ.0d0)THEN
         XIFINP=XIFCEN
        ELSE
         XIFINP=XIFCEN+MAX(DABS(XIFCEN),XIFMIN)*XIFDEV*GAU(IDUM)
        ENDIF
       ENDIF

       IF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4)THEN
        IF(KDEV.EQ.0d0)THEN
         PAR(2)=KCEN
        ELSE
         PAR(2)=KCEN+MAX(DABS(KCEN),KMIN)*KDEV*GAU(IDUM)
        ENDIF
       ENDIF

       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        IF(XISDEV.EQ.0d0)THEN
         XISINP=XISCEN
        ELSE
         XISINP=XISCEN+MAX(DABS(XISCEN),XISMIN)*XISDEV*GAU(IDUM)
        ENDIF
       ENDIF

       IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)THEN
        IF(MSFLAG.EQ.0)THEN
         MSINP=M0**2
        ELSE
         IF(MSDEV.EQ.0d0)THEN
          MSINP=MSCEN
         ELSE
          MSINP=MSCEN+MAX(DABS(MSCEN),MSMIN)*MSDEV*GAU(IDUM)
         ENDIF
        ENDIF
       ENDIF

       IF(MAFLAG.EQ.-5)THEN
        IF(MUDEV.EQ.0d0)THEN
         PAR(4)=MUCEN
        ELSE
         PAR(4)=MUCEN+MAX(DABS(MUCEN),MUMIN)*MUDEV*GAU(IDUM)
        ENDIF
	IF(PAR(4).NE.0d0)SIGMU=DABS(PAR(4))/PAR(4)
        IF(KDEV.EQ.0d0)THEN
         PAR(2)=KCEN
        ELSE
         PAR(2)=KCEN+MAX(DABS(KCEN),KMIN)*KDEV*GAU(IDUM)
        ENDIF
        IF(XIFDEV.EQ.0d0)THEN
         XIFINP=XIFCEN
        ELSE
         XIFINP=XIFCEN+MAX(DABS(XIFCEN),XIFMIN)*XIFDEV*GAU(IDUM)
        ENDIF
        IF(XISDEV.EQ.0d0)THEN
         XISINP=XISCEN
        ELSE
         XISINP=XISCEN+MAX(DABS(XISCEN),XISMIN)*XISDEV*GAU(IDUM)
        ENDIF
       ELSE
        IF(MHDFLAG.EQ.0)THEN
         MHDINP=M0**2
        ELSE
         IF(MHDDEV.EQ.0d0)THEN
          MHDINP=MHDCEN
         ELSE
          MHDINP=MHDCEN+MAX(DABS(MHDCEN),MHDMIN)*MHDDEV*GAU(IDUM)
         ENDIF
        ENDIF
        IF(MHUFLAG.EQ.0)THEN
         MHUINP=M0**2
        ELSE
         IF(MHUDEV.EQ.0d0)THEN
          MHUINP=MHUCEN
         ELSE
          MHUINP=MHUCEN+MAX(DABS(MHUCEN),MHUMIN)*MHUDEV*GAU(IDUM)
         ENDIF
        ENDIF
       ENDIF

       IF(MUPDEV.EQ.0d0)THEN
        MUPINP=MUPCEN
       ELSE
        MUPINP=MUPCEN+MAX(DABS(MUPCEN),MUPMIN)*MUPDEV*GAU(IDUM)
       ENDIF

       IF(MSPDEV.EQ.0d0)THEN
        MSPINP=MSPCEN
       ELSE
        MSPINP=MSPCEN+MAX(DABS(MSPCEN),MSPMIN)*MSPDEV*GAU(IDUM)
       ENDIF

       IF(M3HDEV.EQ.0d0)THEN
        M3HINP=M3HCEN
       ELSE
        M3HINP=M3HCEN+MAX(DABS(M3HCEN),M3HMIN)*M3HDEV*GAU(IDUM)
       ENDIF

      ENDIF

!      WRITE(0,*)""
!      WRITE(0,*)"------------------------------------------------------"
!      WRITE(0,*)""
!      WRITE(0,*)"Point ",ITOT
!      WRITE(0,*)""
!      WRITE(0,*)"MAFLAG=",MAFLAG
!      WRITE(0,*)"M0 =",M0
!      IF(M1FLAG*M2FLAG*M3FLAG.EQ.0)WRITE(0,*)"M12 =",M12
!      WRITE(0,*)"TANB =",PAR(3)
!      IF(MAFLAG.NE.-5)WRITE(0,*)"SIGMU =",SIGMU
!      WRITE(0,*)"A0 =",A0
!      IF(M1FLAG.EQ.1)WRITE(0,*)"M1 =",M1INP
!      IF(M2FLAG.EQ.1)WRITE(0,*)"M2 =",M2INP
!      IF(M3FLAG.EQ.1)WRITE(0,*)"M3 =",M3INP
!      IF(MHUFLAG.EQ.1 .AND. MAFLAG.NE.-5)WRITE(0,*)"MHU =",MHUINP
!      IF(MHDFLAG.EQ.1 .AND. MAFLAG.NE.-5)WRITE(0,*)"MHD =",MHDINP
!      IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)WRITE(0,*)"MS =",MSINP
!      WRITE(0,*)"LAMBDA =",PAR(1)
!      IF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4 .OR. MAFLAG.EQ.-5)
!     . WRITE(0,*)"KAPPA =",PAR(2)
!      IF(ALFLAG.EQ.1)WRITE(0,*)"AL =",ALINP
!      IF(AKFLAG.EQ.1)WRITE(0,*)"AK =",AKINP
!      IF(MAFLAG.EQ.-5)WRITE(0,*)"MUEFF =",PAR(4)
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-5)
!     . WRITE(0,*)"XIF =",XIFINP
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-5)
!     . WRITE(0,*)"XIS =",XISINP
!      WRITE(0,*)"MUP =",MUPINP
!      WRITE(0,*)"MSP =",MSPINP
!      WRITE(0,*)"M3H =",M3HINP
!      WRITE(0,*)

*   Initialization of PROB

      DO I=1,NPROB
       PROB(I)=0d0
      ENDDO
      UNCERTFLAG=0

*   Initialization of algorithm parameters

      SFFLAG=0
      ITRY=1
      IF(ITOT.EQ.1)THEN
       MUINIT=.5d0*SIGMU*DSQRT(PAR(3)*Q2MIN)
      ELSE
       MUINIT=.5d0*PARMEM(4)
      ENDIF
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
      MUFLAG=0
      MUT=0d0
      MDT=0d0
      MST=0d0
      IM=0

*   Guess parameters at Q2/MGUT

      IF(ITOT.EQ.1)THEN
!       WRITE(0,*)"Guesses"
       IF(Q2FIX.EQ.0)THEN
        Q2=MAX(M0**2+4d0*M12**2,Q2MIN)
       ENDIF
       R=(1d0+PAR(3)**2)/(1.29d0*PAR(3)**2)
       K=(1d0-R)*(A0-2.24d0*M12)**2+7.84d0*M12**2
       PAR(5)=ALINP-R/2d0*(A0-2.24d0*M12)-.59d0*M12
       PAR(6)=AKINP
       PAR(7)=(1d0-R/2d0)*M0**2+7.02d0*M12**2-R/6d0*K
       PAR(8)=(1d0-R)*M0**2+6.6d0*M12**2-R/3d0*K
       PAR(9)=M0**2+6.55d0*M12**2
       PAR(10)=M0**2+.52d0*M12**2
       PAR(11)=M0**2+.15d0*M12**2
       PAR(12)=A0-R*(A0-2.24d0*M12)-3.97d0*M12
       PAR(13)=A0-R/6d0*(A0-2.24d0*M12)-3.93d0*M12
       PAR(14)=A0-.69d0*M12
       PAR(15)=M0**2+7.02d0*M12**2
       PAR(16)=M0**2+6.6d0*M12**2
       PAR(17)=M0**2+6.55d0*M12**2
       PAR(18)=M0**2+.52d0*M12**2
       PAR(19)=M0**2+.15d0*M12**2
       PAR(20)=.4d0*M1INP
       PAR(21)=.8d0*M2INP
       PAR(22)=2.4d0*M3INP
       PAR(23)=DSQRT(Q2)
       PAR(24)=DSQRT(Q2)
       PAR(25)=A0-.69d0*M12
       DELMB=.1d0
       DELML=0d0
       DEL1=0d0
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        PAR(2)=PAR(1)/5d0
       ELSEIF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4)THEN
        XIFGUT=0d0
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        MSM=M0**2
       ELSEIF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)THEN
        XISGUT=0d0
       ENDIF
       IF(MAFLAG.EQ.-5)THEN
        MUM=5d0*MAX(M0**2,M12**2/5d0)
        MDM=MAX(M0**2,M12**2/5d0)
        MSM=MAX(M0**2,M12**2/5d0)
       ELSE
        PAR(4)=MUINIT
       ENDIF

      ELSE
!	WRITE(0,*)"Memorized"
       IF(Q2FIX.EQ.0)THEN
        Q2=Q2MEM
       ENDIF
       DO I=5,25
        PAR(I)=PARMEM(I)
       ENDDO
       DELMB=DELMEM
       DELML=0d0
       DEL1=0d0
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        PAR(2)=PARMEM(2)
       ELSE
        XIFGUT=XIFMEM
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        MSM=MSMEM
       ELSE
        XISGUT=XISMEM
       ENDIF
       IF(MAFLAG.EQ.-5)THEN
        MUM=MUMEM
        MDM=MDMEM
        MSM=MSMEM
       ELSE
        PAR(4)=MUINIT
       ENDIF

      ENDIF
      IFAIL=0

!      WRITE(0,*)""
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
!      ELSEIF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4)THEN
!      WRITE(0,*)"XIF =",XIFGUT
!       ENDIF
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
!       WRITE(0,*)"MS =",MSM
!      ELSEIF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)THEN
!       WRITE(0,*)"XIS =",XISGUT
!      ENDIF
!      IF(MAFLAG.EQ.-5)THEN
!       WRITE(0,*)"MHU =",MUM
!       WRITE(0,*)"MHD =",MDM
!       WRITE(0,*)"MS =",MSM
!      ELSE
!       WRITE(0,*)"MU =",PAR(4)
!      ENDIF
!      WRITE(0,*)"MA =",PAR(23)
!      WRITE(0,*)"MP =",PAR(24)
!      WRITE(0,*)""

*   Check for mu = 0

      IF(PAR(4).EQ.0d0)THEN
       IFAIL=9
       GOTO 11
      ENDIF

*   Guess for GUT scale and GUT couplings

      CALL RGES(PAR,PROB,IFAIL,0)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"RGE integration problem 1"
       IF(DABS(MUINIT).LT.1d2 .AND. DABS(MUSTEP).GT.1d0
     .    .AND. MAFLAG.NE.-5)THEN
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

      CALL RGESINV(PAR,IFAIL)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"RGE integration problem 2"
       IF(DABS(MUINIT).LT.1d2 .AND. DABS(MUSTEP).GT.1d0
     .    .AND. MAFLAG.NE.-5)THEN
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
*   or (mHu and mHd and mS), no iteration for the latter case

      ITERMU=0
 22   ITERMU=ITERMU+1
!      WRITE(0,*)"ITERMU =",ITERMU
!      WRITE(0,*)""
!      WRITE(0,*)""

      CALL RUNPAR(PAR)

      CALL MSFERM(PAR,IFAIL,0)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"Negative sfermion mass squared"
       IF(MAFLAG.NE.-5)THEN
        IF(DABS(MUINIT).LT.1d2 .AND. DABS(MUSTEP).GT.1d0)THEN
         MUSTEP=MUSTEP/2d0
         MUINIT=MUSTEP
         SFFLAG=1
         ITRY=ITRY+1
         GOTO 1
        ELSE
!         WRITE(0,*)"Exit (IFAIL=8)"
!         WRITE(0,*)""
!         WRITE(0,*)""
         IFAIL=8
         GOTO 11
        ENDIF
       ENDIF
      ENDIF

      CALL MINIMIZE(PAR,CHECK)
      IF(ITER.GT.10.AND.DETM.GE.0d0)MUFLAG=1

*   End of the internal loop

      IF((CHECK.GT.1d-12.AND.ITERMU.LT.10).OR.
     .   (CHECK.GT.1d-8.AND.ITERMU.LT.50).OR.
     .   (CHECK.GT.1d-6.AND.ITERMU.LT.100))GOTO 22
      IF(CHECK.GT.1d-4)THEN
!       WRITE(0,*)"No convergence 1"
!       WRITE(0,*)"Exit (IFAIL=14)"
!       WRITE(0,*)""
!       WRITE(0,*)""
       IFAIL=14
       GOTO 11
      ENDIF

      CALL RGES(PAR,PROB,IFAIL,1)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"RGE integration problem 3"
       IF(DABS(MUINIT).LT.1d2 .AND. DABS(MUSTEP).GT.1d0
     .    .AND. MAFLAG.NE.-5)THEN
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

      CALL RGESUNI(PAR,IFAIL,GUTEST)
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"RGE integration problem 4"
       IF(DABS(MUINIT).LT.1d2 .AND. DABS(MUSTEP).GT.1d0
     .    .AND. MAFLAG.NE.-5)THEN
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

      IF((GUTEST.GT.1d-12.AND.ITER.LT.10).OR.
     .   (GUTEST.GT.1d-8.AND.ITER.LT.50).OR.
     .   (GUTEST.GT.1d-6.AND.ITER.LT.100))GOTO 21
      IF(GUTEST.GT.1d-4)THEN
!       WRITE(0,*)"No convergence 2"
       IF(DABS(MUSTEP).GT.1d0 .AND. MAFLAG.NE.-5)THEN
        MUSTEP=MUSTEP/2d0
        MUFAIL=MUFAIL+MUSTEP
        MUFLAG=0
        MUT=0d0
        MDT=0d0
        MST=0d0
        IM=0
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
!        WRITE(0,*)"Exit (IFAIL=15)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=15
        GOTO 11
       ENDIF
      ENDIF

*   Check if correct EWSB

      IF(DETM.LE.0d0 .AND. MAFLAG.NE.-5)THEN
!       WRITE(0,*)"Convergence in a false minimum"
       IF(MUFLAG.EQ.1 .AND. DABS(MUSTEP).GT.1d0)THEN
        MUSTEP=MUSTEP/2d0
        MUFAIL=MUFAIL+MUSTEP
        MUFLAG=0
        MUT=0d0
        MDT=0d0
        MST=0d0
        IM=0
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
       ELSEIF(MUFLAG.EQ.0 .AND. DABS(MUSTEP).GT.1d0)THEN
        MUSTEP=MUSTEP/2d0
        MUFAIL=MUFAIL-MUSTEP
        MUFLAG=0
        MUT=0d0
        MDT=0d0
        MST=0d0
        IM=0
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
!        WRITE(0,*)"Exit (IFAIL=17)"
!        WRITE(0,*)""
!        WRITE(0,*)""
        IFAIL=17
        GOTO 11
       ENDIF
      ENDIF

*   Computation of sfermion masses

      CALL MSFERM(PAR,IFAIL,1)
      IF(IFAIL.NE.0)GOTO 11

*   Computation of Higgs masses

      CALL MHIGGS(PAR,PROB,IFAIL)
      IF(IFAIL.EQ.-1)IFAIL=16
      IF(IFAIL.NE.0)GOTO 11

*   Computation of gluino mass

      CALL GLUINO(PAR)

*   Computation of chargino masses

      CALL CHARGINO(PAR)

*   Computation of neutralino masses

      CALL NEUTRALINO(PAR)

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

      M32=CGR*DSQRT(Q2/3d0)
      CALL RELDEN(PAR,PROB)

*   Sparticle decays

      IF(NMSFLAG.NE.0)CALL NMSDECAY(PAR)

*   Exp. constraints on sparticles (LHC)

      IF(CFLAG(5).EQ.1)CALL Higgsino_CMS_Trilep(PROB)

*   Computation of the fine-tuning

      CALL FTPAR(PAR,1)

*   Check for problems

      DO I=1,NPROB
       IF(PROB(I).NE.0d0)THEN
!        WRITE(0,*)"PROB",I,PROB(I)
        IFAIL=10
       ENDIF
      ENDDO
!      WRITE(0,*)""

*   Recording of the results

 11   CALL MCMCSTEPSP(PAR,PROB,NPROB,IFAIL)
      IF(ITOT.NE.1 .OR. ISEED.EQ.1)CALL OUTPUT(PAR,PROB,IFAIL)
      IF(IFAIL.EQ.0)THEN
       TOT=TOT+1
       M0N=MIN(M0,M0N)
       M0NN=MAX(M0,M0NN)
       M12N=MIN(M12,M12N)
       M12NN=MAX(M12,M12NN)
       TBN=MIN(PAR(3),TBN)
       TBNN=MAX(PAR(3),TBNN)
       A0N=MIN(A0,A0N)
       A0NN=MAX(A0,A0NN)
       LN=MIN(PAR(1),LN)
       LNN=MAX(PAR(1),LNN)
       ALN=MIN(ALGUT,ALN)
       ALNN=MAX(ALGUT,ALNN)
       AKN=MIN(AKGUT,AKN)
       AKNN=MAX(AKGUT,AKNN)
       KN=MIN(PAR(2),KN)
       KNN=MAX(PAR(2),KNN)
       MUN=MIN(PAR(4),MUN)
       MUNN=MAX(PAR(4),MUNN)
       MHDN=MIN(MHDGUT,MHDN)
       MHDNN=MAX(MHDGUT,MHDNN)
       MHUN=MIN(MHUGUT,MHUN)
       MHUNN=MAX(MHUGUT,MHUNN)
       M1N=MIN(M1GUT,M1N)
       M1NN=MAX(M1GUT,M1NN)
       M2N=MIN(M2GUT,M2N)
       M2NN=MAX(M2GUT,M2NN)
       M3N=MIN(M3GUT,M3N)
       M3NN=MAX(M3GUT,M3NN)
       XIFN=MIN(XIFGUT,XIFN)
       XIFNN=MAX(XIFGUT,XIFNN)
       XISN=MIN(XISGUT,XISN)
       XISNN=MAX(XISGUT,XISNN)
       MSN=MIN(MSGUT,MSN)
       MSNN=MAX(MSGUT,MSNN)
       MUPN=MIN(MUPGUT,MUPN)
       MUPNN=MAX(MUPGUT,MUPNN)
       MSPN=MIN(MSPGUT,MSPN)
       MSPNN=MAX(MSPGUT,MSPNN)
       M3HN=MIN(M3HGUT,M3HN)
       M3HNN=MAX(M3HGUT,M3HNN)
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

      INTEGER I,NLINE,INL,ICH,IX,IVAL,Q2FIX,MCFLAG,TOTMIN,TOTMAX,NMAX
      INTEGER N0,NLOOP,NBER,NPAR,ERR,NTOT,ISEED,GMUFLAG,HFLAG,Z3FLAG
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,MSFLAG,VFLAG
      INTEGER AKFLAG,ALFLAG,OMGFLAG,MAFLAG,MOFLAG,PFLAG,NMSFLAG
      INTEGER GRFLAG,CFLAG(5)

      DOUBLE PRECISION PAR(*),VAL
      DOUBLE PRECISION ACC,XITLA,XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,Q2,SIGMU,Q2MIN
      DOUBLE PRECISION M0CEN,M0DEV,M12CEN,M12DEV,TBCEN,TBDEV,A0CEN
      DOUBLE PRECISION A0DEV,M1CEN,M1DEV,M2CEN,M2DEV,M3CEN,M3DEV
      DOUBLE PRECISION MHDCEN,MHDDEV,MHUCEN,MHUDEV,MSCEN,MSDEV,LCEN
      DOUBLE PRECISION LDEV,KCEN,KDEV,ALCEN,ALDEV,AKCEN,AKDEV,MUCEN
      DOUBLE PRECISION MUDEV,XIFCEN,XIFDEV,XISCEN,XISDEV,MUPCEN
      DOUBLE PRECISION MUPDEV,MSPCEN,MSPDEV,M3HCEN,M3HDEV,XCEN,XDEV,X
      DOUBLE PRECISION M0MIN,M12MIN,TBMIN,A0MIN,M1MIN,M2MIN,M3MIN
      DOUBLE PRECISION MHDMIN,MHUMIN,MSMIN,LMIN,KMIN,ALMIN,AKMIN
      DOUBLE PRECISION MUMIN,XIFMIN,XISMIN,MUPMIN,MSPMIN,M3HMIN
      DOUBLE PRECISION M32,CGR,MPL

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/SIGMU/SIGMU
      COMMON/MCMCPAR/M0CEN,M0DEV,M12CEN,M12DEV,TBCEN,TBDEV,A0CEN,
     . A0DEV,M1CEN,M1DEV,M2CEN,M2DEV,M3CEN,M3DEV,MHDCEN,MHDDEV,
     . MHUCEN,MHUDEV,MSCEN,MSDEV,LCEN,LDEV,KCEN,KDEV,ALCEN,ALDEV,
     . AKCEN,AKDEV,MUCEN,MUDEV,XIFCEN,XIFDEV,XISCEN,XISDEV,MUPCEN,
     . MUPDEV,MSPCEN,MSPDEV,M3HCEN,M3HDEV,XCEN,XDEV,X,
     . M0MIN,M12MIN,TBMIN,A0MIN,M1MIN,M2MIN,M3MIN,
     . MHDMIN,MHUMIN,MSMIN,LMIN,KMIN,ALMIN,AKMIN,
     . MUMIN,XIFMIN,XISMIN,MUPMIN,MSPMIN,M3HMIN
      COMMON/STEPS/NTOT,ISEED,TOTMIN,TOTMAX,NMAX
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . MSFLAG,AKFLAG,ALFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/PFLAG/PFLAG
      COMMON/NMSFLAG/NMSFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/MCFLAG/MCFLAG
      COMMON/VFLAG/VFLAG
      COMMON/CFLAG/CFLAG
      COMMON/M32/M32,CGR,MPL,GRFLAG

*   INITIALIZATION OF THE SUSY PARAMETERS
      DO I=1,NPAR
       PAR(I)=0d0
      ENDDO

*   INITIALIZATION OF THE SCANNING PARAMETERS
      SIGMU=1d99
      M0CEN=1d99
      M0DEV=1d99
      M12CEN=1d99
      M12DEV=1d99
      TBCEN=1d99
      TBDEV=1d99
      A0CEN=1d99
      A0DEV=1d99
      M1CEN=1d99
      M1DEV=1d99
      M2CEN=1d99
      M2DEV=1d99
      M3CEN=1d99
      M3DEV=1d99
      MHDCEN=1d99
      MHDDEV=1d99
      MHUCEN=1d99
      MHUDEV=1d99
      MSCEN=1d99
      MSDEV=1d99
      LCEN=1d99
      LDEV=1d99
      KCEN=1d99
      KDEV=1d99
      ALCEN=1d99
      ALDEV=1d99
      AKCEN=1d99
      AKDEV=1d99
      MUCEN=1d99
      MUDEV=1d99
      XIFCEN=1d99
      XIFDEV=1d99
      XISCEN=1d99
      XISDEV=1d99
      MUPCEN=0d0
      MUPDEV=1d99
      MSPCEN=0d0
      MSPDEV=1d99
      M3HCEN=0d0
      M3HDEV=1d99
      XCEN=1d99
      XDEV=1d0
      M0MIN=0d0
      M12MIN=0d0
      TBMIN=0d0
      A0MIN=0d0
      M1MIN=0d0
      M2MIN=0d0
      M3MIN=0d0
      MHDMIN=0d0
      MHUMIN=0d0
      MSMIN=0d0
      LMIN=0d0
      KMIN=0d0
      ALMIN=0d0
      AKMIN=0d0
      MUMIN=0d0
      XIFMIN=0d0
      XISMIN=0d0
      MUPMIN=0d0
      MSPMIN=0d0
      M3HMIN=0d0
      NTOT=0
      TOTMIN=0
      TOTMAX=1000000
      NMAX=1000000
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
      M1FLAG=0
      M2FLAG=0
      M3FLAG=0
      MHDFLAG=0
      MHUFLAG=0
      MSFLAG=1
      ALFLAG=0
      AKFLAG=0
      MCFLAG=0
      GRFLAG=0
      DO I=1,5
       CFLAG(I)=1
      ENDDO

*   DEFAULT VALUE FOR THE RANDOM SEED
      ISEED=-1

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
       IF(IX.EQ.8) PFLAG=IVAL
       IF(IX.EQ.9) OMGFLAG=IVAL
       IF(IX.EQ.11) GMUFLAG=IVAL
       IF(IX.EQ.13) NMSFLAG=IVAL
       IF(IX.EQ.14) VFLAG=IVAL
       IF(IX.EQ.15) MOFLAG=IVAL
       IF(IX.EQ.17) CFLAG(1)=IVAL
       IF(IX.EQ.18) CFLAG(2)=IVAL
       IF(IX.EQ.19) CFLAG(3)=IVAL
       IF(IX.EQ.20) CFLAG(4)=IVAL
       IF(IX.EQ.22) CFLAG(5)=IVAL

*   READ SMINPUTS
      ELSEIF(CHBLCK(1:8).EQ.'SMINPUTS')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.2) GF=VAL
       IF(IX.EQ.3) ALSMZ=VAL
       IF(IX.EQ.4) MZ=VAL
       IF(IX.EQ.5) MB=VAL
       IF(IX.EQ.6) MT=VAL
       IF(IX.EQ.7) MTAU=VAL

*   READ GUT PARAMETERS, SIGMU, Q2 AND TANBETA
      ELSEIF(CHBLCK(1:6).EQ.'MINPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0.AND.Q2.EQ.0d0) Q2=VAL**2
       IF(IX.EQ.1) M0CEN=VAL
       IF(IX.EQ.16) M0DEV=VAL
       IF(IX.EQ.17) M0MIN=VAL
       IF(IX.EQ.2) M12CEN=VAL
       IF(IX.EQ.26) M12DEV=VAL
       IF(IX.EQ.27) M12MIN=VAL
       IF(IX.EQ.3) TBCEN=VAL
       IF(IX.EQ.36) TBDEV=VAL
       IF(IX.EQ.37) TBMIN=VAL
       IF(IX.EQ.4) SIGMU=VAL
       IF(IX.EQ.5) A0CEN=VAL
       IF(IX.EQ.56) A0DEV=VAL
       IF(IX.EQ.57) A0MIN=VAL
       IF(IX.EQ.6) CGR=VAL

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
       IF(IX.EQ.21) MHDCEN=VAL
       IF(IX.EQ.216) MHDDEV=VAL
       IF(IX.EQ.217) MHDMIN=VAL
       IF(IX.EQ.22) MHUCEN=VAL
       IF(IX.EQ.226) MHUDEV=VAL
       IF(IX.EQ.227) MHUMIN=VAL
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
       IF(IX.EQ.70) MSCEN=VAL
       IF(IX.EQ.706) MSDEV=VAL
       IF(IX.EQ.707) MSMIN=VAL
       IF(IX.EQ.72) M3HCEN=VAL
       IF(IX.EQ.726) M3HDEV=VAL
       IF(IX.EQ.727) M3HMIN=VAL

*   READ STEPS
       ELSEIF(CHBLCK(1:6).EQ.'STEPS')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.0) NTOT=IVAL
       IF(IX.EQ.1) ISEED=IVAL
       IF(IX.EQ.2) HFLAG=IVAL
       IF(IX.EQ.3) MCFLAG=IVAL
       IF(IX.EQ.4) TOTMIN=IVAL
       IF(IX.EQ.5) TOTMAX=IVAL
       IF(IX.EQ.6) NMAX=IVAL

      ENDIF

      GOTO 21

*   END OF READING FROM INPUT FILE

*   Check for errors

 29   ERR=0
      IF(CFLAG(5).EQ.1 .AND. NMSFLAG.EQ.0)THEN
       WRITE(0,2)"CMS CHARG(NEUTRAL)INO CONSTRAINTS CANNOT BE CHECKED ",
     .  "IF NMSDECAY IS NOT CALLED"
       ERR=1
      ENDIF
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
      IF(VFLAG.LT.0 .OR. VFLAG.GT.1)THEN
       WRITE(0,1)"VFLAG MUST BE IN [0-1]"
       ERR=1
      ENDIF
      IF(MOFLAG.LT.0 .OR. MOFLAG.GT.7)THEN
       WRITE(0,1)"MOFLAG MUST BE IN [0-7]"
       ERR=1
      ENDIF
      IF(XDEV.EQ.1d99)THEN
       WRITE(0,1)"XDEV MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M0CEN.EQ.1d99)THEN
       WRITE(0,1)"M0CEN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(A0CEN.EQ.1d99)THEN
       WRITE(0,1)"A0CEN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(LCEN.EQ.1d99)THEN
       WRITE(0,1)"LCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ELSEIF(LCEN.LE.0d0)THEN
       WRITE(0,1)"LCEN MUST BE STRICTLY POSITIVE"
       ERR=1
      ENDIF
      IF(TBCEN.EQ.1d99)THEN
       WRITE(0,1)"TBCEN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ELSEIF(TBCEN.LE.0d0)THEN
       WRITE(0,1)"TBCEN MUST BE STRICTLY POSITIVE"
       ERR=1
      ENDIF
      IF(DABS(SIGMU).NE.1d0 .AND. MUCEN.EQ.1d99)THEN
       WRITE(0,1)"SIGMU IS EITHER 1 OR -1 IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(SIGMU.NE.1d99 .AND. MUCEN.NE.1d99)THEN
       WRITE(0,1)"BOTH MUCEN AND SIGMU CANNOT BE GIVEN"
       ERR=1
      ENDIF
      IF(M12CEN.EQ.1d99 .AND. M12DEV.NE.1d99)THEN
       WRITE(0,1)"M12CEN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(M1CEN.EQ.1d99 .AND. M1DEV.NE.1d99)THEN
       WRITE(0,1)"M1CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M2CEN.EQ.1d99 .AND. M2DEV.NE.1d99)THEN
       WRITE(0,1)"M2CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M3CEN.EQ.1d99 .AND. M3DEV.NE.1d99)THEN
       WRITE(0,1)"M3CEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MHDCEN.EQ.1d99 .AND. MHDDEV.NE.1d99)THEN
       WRITE(0,1)"MHDCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MHUCEN.EQ.1d99 .AND. MHUDEV.NE.1d99)THEN
       WRITE(0,1)"MHUCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MSCEN.EQ.1d99 .AND. MSDEV.NE.1d99)THEN
       WRITE(0,1)"MSCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(KCEN.EQ.1d99 .AND. KDEV.NE.1d99)THEN
       WRITE(0,1)"KCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(ALCEN.EQ.1d99 .AND. ALDEV.NE.1d99)THEN
       WRITE(0,1)"ALCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(AKCEN.EQ.1d99 .AND. AKDEV.NE.1d99)THEN
       WRITE(0,1)"AKCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MUCEN.EQ.1d99 .AND. MUDEV.NE.1d99)THEN
       WRITE(0,1)"MUCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XIFCEN.EQ.1d99 .AND. XIFDEV.NE.1d99)THEN
       WRITE(0,1)"XIFCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XISCEN.EQ.1d99 .AND. XISDEV.NE.1d99)THEN
       WRITE(0,1)"XISCEN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MUCEN.NE.1d99)THEN
       IF(MHDCEN.NE.1d99)THEN
        WRITE(0,1)
     .  "BOTH MUCEN AND MHDCEN CANNOT BE GIVEN IN BLOCK EXTPAR"
        ERR=1
       ENDIF
       IF(MHUCEN.NE.1d99)THEN
        WRITE(0,1)
     .  "BOTH MUCEN AND MHUCEN CANNOT BE GIVEN IN BLOCK EXTPAR"
        ERR=1
       ENDIF
       IF(MSCEN.NE.1d99)THEN
        WRITE(0,1)
     .  "BOTH MUCEN AND MSCEN CANNOT BE GIVEN IN BLOCK EXTPAR"
        ERR=1
       ENDIF
       IF(MUCEN.EQ.0d0)THEN
        WRITE(0,1)"MUCEN MUST BE NON ZERO"
        ERR=1
       ENDIF
      ELSE
       IF(XIFCEN.NE.1d99 .AND. KCEN.NE.1d99)THEN
        WRITE(0,1)
     .  "BOTH XIFCEN AND KCEN CANNOT BE GIVEN IN BLOCK EXTPAR"
        ERR=1
       ENDIF
       IF(XISCEN.NE.1d99 .AND. MSCEN.NE.1d99)THEN
        WRITE(0,1)
     .  "BOTH XISCEN AND MSCEN CANNOT BE GIVEN IN BLOCK EXTPAR"
        ERR=1
       ENDIF
      ENDIF

*   Set default values

      IF(MUCEN.NE.1d99)THEN
       IF(KCEN.EQ.1d99)KCEN=0d0
       IF(XIFCEN.EQ.1d99)XIFCEN=0d0
       IF(XISCEN.EQ.1d99)XISCEN=0d0
      ELSE
       IF(KCEN.EQ.1d99 .AND. XIFCEN.EQ.1d99)XIFCEN=0d0
       IF(MSCEN.EQ.1d99 .AND. XISCEN.EQ.1d99)XISCEN=0d0
      ENDIF

      IF(M0DEV.EQ.1d99)M0DEV=0d0
      IF(M12DEV.EQ.1d99)M12DEV=0d0
      IF(A0DEV.EQ.1d99)A0DEV=0d0
      IF(TBDEV.EQ.1d99)TBDEV=0d0
      IF(M1DEV.EQ.1d99)M1DEV=0d0
      IF(M2DEV.EQ.1d99)M2DEV=0d0
      IF(M3DEV.EQ.1d99)M3DEV=0d0
      IF(MHUDEV.EQ.1d99)MHUDEV=0d0
      IF(MHDDEV.EQ.1d99)MHDDEV=0d0
      IF(MSDEV.EQ.1d99)MSDEV=0d0
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

*   Set MAFLAG, SCANFLAGS

      IF(MUCEN.NE.1d99)THEN
       MAFLAG=-5
       MHDFLAG=1
       MHUFLAG=1
      ELSE
       MAFLAG=-1
       IF(KCEN.NE.1d99)MAFLAG=MAFLAG-2
       IF(MSCEN.NE.1d99)MAFLAG=MAFLAG-1
       IF(MHDCEN.NE.1d99)MHDFLAG=1
       IF(MHUCEN.NE.1d99)MHUFLAG=1
       IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)THEN
        IF(MSCEN.EQ.M0CEN**2 .AND. MSDEV.EQ.M0DEV)MSFLAG=0
       ENDIF
      ENDIF

      IF(M1CEN.NE.1d99)M1FLAG=1
      IF(M2CEN.NE.1d99)M2FLAG=1
      IF(M3CEN.NE.1d99)M3FLAG=1
      IF(ALCEN.NE.1d99)ALFLAG=1
      IF(AKCEN.NE.1d99)AKFLAG=1

      IF(M1FLAG*M2FLAG*M3FLAG.EQ.0)THEN
       IF(M12CEN.EQ.1d99)THEN
        WRITE(0,1)"M12CEN MUST BE GIVEN IN BLOCK MINPAR"
        ERR=1
       ENDIF
      ELSE
       IF(M12CEN.NE.1d99)WRITE(0,1)"WARNING: M12 IS NOT USED"
      ENDIF

      IF(KCEN.EQ.0d0 .AND. KDEV.EQ.0d0 .AND. KMIN.EQ.0d0)THEN
       IF((AKCEN.NE.0d0 .AND. AKCEN.NE.1d99) .OR.
     .  (AKCEN.EQ.0d0 .AND. AKDEV.NE.0d0 .AND. AKMIN.NE.0d0))THEN
        WRITE(0,1)
     .  "WARNING KCEN=KDEV=KMIN0 => AKCEN=AKDEV=AKMIN0"
       ENDIF
       AKCEN=0d0
       AKDEV=0d0
       AKMIN=0d0
       AKFLAG=1
      ENDIF

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

      IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4 .OR.
     . MUPCEN.NE.0d0 .OR. MUPDEV.NE.0d0 .OR. MSPCEN.NE.0d0 .OR. 
     . MSPDEV.NE.0d0 .OR. XIFCEN.NE.0d0 .OR. XIFDEV.NE.0d0 .OR.
     . XISCEN.NE.0d0 .OR. XISDEV.NE.0d0 .OR. M3HCEN.NE.0d0 .OR.
     . M3HDEV.NE.0d0)THEN
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
 2    FORMAT(A,A)

      END


      SUBROUTINE OUTPUT(PAR,PROB,IFAIL)

*********************************************************************
*   Subroutine writing all the results in the the output file.
*********************************************************************

      IMPLICIT NONE

      INTEGER NBIN,I,NRES,IRES,GRFLAG,NSUSY,NGUT,NMES,IMAX,IFAIL
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
      DOUBLE PRECISION BRSUSY(5),WIDTH(5),OMG,OMGMIN,OMGMAX
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,MHUS,MHDS,MSS
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU,Q2
      DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      DOUBLE PRECISION MGUT,g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT
      DOUBLE PRECISION HBOTGUT,HTAUGUT,M1GUT,M2GUT,M3GUT,ALGUT,AKGUT
      DOUBLE PRECISION ATGUT,ABGUT,ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT
      DOUBLE PRECISION MQ3GUT,MU3GUT,MD3GUT,MQGUT,MUGUT,MDGUT,ML3GUT
      DOUBLE PRECISION ME3GUT,MLGUT,MEGUT,M0,M12,A0,SIGMU
      DOUBLE PRECISION M1INP,M2INP,M3INP,MHDINP,MHUINP,ALINP,AKINP
      DOUBLE PRECISION XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      DOUBLE PRECISION XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION Xf,sigmaV,vcsll,vcsbb,x(100),dNdx(100),EMIN
      DOUBLE PRECISION sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
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
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2),toptot
      DOUBLE PRECISION FTSUSY(NSUSY+2),FTGUT(NGUT+2),FTMES(NMES+2)
      DOUBLE PRECISION PX,PA(6),PB(2),PL(7),PK(8)
      DOUBLE PRECISION MHmin,MHmax
      DOUBLE PRECISION M0CEN,M0DEV,M12CEN,M12DEV,TBCEN,TBDEV,A0CEN
      DOUBLE PRECISION A0DEV,M1CEN,M1DEV,M2CEN,M2DEV,M3CEN,M3DEV
      DOUBLE PRECISION MHDCEN,MHDDEV,MHUCEN,MHUDEV,MSCEN,MSDEV,LCEN
      DOUBLE PRECISION LDEV,KCEN,KDEV,ALCEN,ALDEV,AKCEN,AKDEV,MUCEN
      DOUBLE PRECISION MUDEV,XIFCEN,XIFDEV,XISCEN,XISDEV,MUPCEN
      DOUBLE PRECISION MUPDEV,MSPCEN,MSPDEV,M3HCEN,M3HDEV,XCEN,XDEV,XX
      DOUBLE PRECISION M0MIN,M12MIN,TBMIN,A0MIN,M1MIN,M2MIN,M3MIN
      DOUBLE PRECISION MHDMIN,MHUMIN,MSMIN,LMIN,KMIN,ALMIN,AKMIN
      DOUBLE PRECISION MUMIN,XIFMIN,XISMIN,MUPMIN,MSPMIN,M3HMIN
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
      COMMON/SIGMU/SIGMU
      COMMON/SOFTGUT/M0,M12,A0
      COMMON/INPPAR/M1INP,M2INP,M3INP,MHDINP,MHUINP,ALINP,AKINP,
     . XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
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
      COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/MGUT/MGUT
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SUSYMH/MHUS,MHDS,MSS
      COMMON/QMHIGGS/MHUQ,MHDQ,MSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     . HBOTGUT,HTAUGUT
      COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     . ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT,MU3GUT,MD3GUT,
     . MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      COMMON/GUTEXT/XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,vcsll,vcsbb,
     .      x,dNdx,EMIN,NBIN
      COMMON/MICROMG2/sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      COMMON/FINETUN/FTSUSY,FTGUT,FTMES
      COMMON/EFFCOUP/PX,PA,PB,PL,PK
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/LHCSIG/SIG
      COMMON/HIGGSFIT/MHmin,MHmax
      COMMON/MCMCPAR/M0CEN,M0DEV,M12CEN,M12DEV,TBCEN,TBDEV,A0CEN,
     . A0DEV,M1CEN,M1DEV,M2CEN,M2DEV,M3CEN,M3DEV,MHDCEN,MHDDEV,
     . MHUCEN,MHUDEV,MSCEN,MSDEV,LCEN,LDEV,KCEN,KDEV,ALCEN,ALDEV,
     . AKCEN,AKDEV,MUCEN,MUDEV,XIFCEN,XIFDEV,XISCEN,XISDEV,MUPCEN,
     . MUPDEV,MSPCEN,MSPDEV,M3HCEN,M3HDEV,XCEN,XDEV,XX,
     . M0MIN,M12MIN,TBMIN,A0MIN,M1MIN,M2MIN,M3MIN,
     . MHDMIN,MHUMIN,MSMIN,LMIN,KMIN,ALMIN,AKMIN,
     . MUMIN,XIFMIN,XISMIN,MUPMIN,MSPMIN,M3HMIN
      COMMON/M32/M32,CGR,MPL,GRFLAG

      IF(IFAIL.NE.0)RETURN

      IRES=11
      NRES=26+IRES

      RES(1)=M0
      RES(2)=M12
      RES(3)=A0
      RES(4)=ALGUT
      RES(5)=AKGUT
      RES(6)=PAR(1)
      RES(7)=PAR(2)
      RES(8)=PAR(3)
      RES(9)=PAR(4)
      RES(10)=XIFGUT
      RES(11)=XISGUT

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
      RES(IRES+26)=FTGUT(NGUT+1)

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
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG
      INTEGER MSFLAG,AKFLAG,ALFLAG,OMGFLAG,MAFLAG,MOFLAG

      DOUBLE PRECISION M0N,M0NN,M12N,M12NN,TBN,TBNN,A0N,A0NN
      DOUBLE PRECISION M1N,M1NN,M2N,M2NN,M3N,M3NN,MHDN,MHDNN
      DOUBLE PRECISION MHUN,MHUNN,MSN,MSNN,LN,LNN,KN,KNN,ALN,ALNN
      DOUBLE PRECISION AKN,AKNN,MUN,MUNN,XIFN,XIFNN,XISN,XISNN
      DOUBLE PRECISION MUPN,MUPNN,MSPN,MSPNN,M3HN,M3HNN,DEV

      COMMON/BOUNDS/M0N,M0NN,M12N,M12NN,TBN,TBNN,A0N,A0NN,
     . M1N,M1NN,M2N,M2NN,M3N,M3NN,MHDN,MHDNN,MHUN,MHUNN,
     . MSN,MSNN,LN,LNN,KN,KNN,ALN,ALNN,AKN,AKNN,MUN,MUNN,
     . XIFN,XIFNN,XISN,XISNN,MUPN,MUPNN,MSPN,MSPNN,M3HN,M3HNN
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . MSFLAG,AKFLAG,ALFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/CFLAG/CFLAG

      WRITE(0,10)"NMSSMTools scan info               "
      WRITE(0,10)"Version number: 5.6.2              "
      WRITE(0,*)
      WRITE(0,*)
      WRITE(0,10)"Number of points:                  "
      WRITE(0,*)
      WRITE(0,10)"  scanned                          ",NTOT
      WRITE(0,10)"  mu=0                             ",NFAIL(9)
      WRITE(0,10)"  no electroweak symmetry breaking ",NFAIL(17)
      S=0
      DO I=1,7
       S=S+NFAIL(I)
      ENDDO
      WRITE(0,10)"  with mh1^2 or ma1^2 or mhc^2 < 0 ",S
      WRITE(0,10)"  with m_sfermion^2 < 0            ",NFAIL(8)
      WRITE(0,10)"  violating  constraints           ",NFAIL(10)
      S=NFAIL(11)+NFAIL(12)+NFAIL(13)
      WRITE(0,10)"  RGE integration problem          ",S
      S=NFAIL(14)+NFAIL(15)+NFAIL(16)
      WRITE(0,10)"  convergence problem              ",S
      WRITE(0,*)
      WRITE(0,10)"Remaining good points              ",TOT
      WRITE(0,*)
      WRITE(0,*)

      IF(OMGFLAG.EQ.0 .AND. CFLAG(1).EQ.0 .AND. CFLAG(2).EQ.0 .AND.
     .   CFLAG(3).EQ.0 .AND. CFLAG(4).EQ.0 .AND. CFLAG(5).EQ.0
     .  .AND. GMUFLAG.EQ.0)THEN
       WRITE(0,20)"Contraints taken into account: none               "
      ELSE
       WRITE(0,20)"Contraints taken into account:                    "
       IF(OMGFLAG.GT.0)
     .  WRITE(0,20)" - Relic density from Planck +/- 10% [0.107,0.131]"
       IF(OMGFLAG.LT.0)
     .  WRITE(0,20)" - Relic density from Planck upper bound < 0.131  "
       IF(IABS(OMGFLAG).EQ.2 .OR. IABS(OMGFLAG).EQ.4)
     .  WRITE(0,20)" - DM direct detection                            "
       IF(IABS(OMGFLAG).EQ.3 .OR. IABS(OMGFLAG).EQ.4)
     .  WRITE(0,20)" - DM indirect detection                          "
       IF(CFLAG(1).EQ.1)
     .  WRITE(0,20)" - Landau poles and false minima                  "
       IF(CFLAG(2).EQ.1)
     .  WRITE(0,20)" - LEP/Tevatron Higgs+sparticle                   "
       IF(CFLAG(3).EQ.1)
     .  WRITE(0,20)" - LHC Higgs                                      "
       IF(CFLAG(4).EQ.1)
     .  WRITE(0,20)" - Upsilon, B and K decays                        "
       IF(CFLAG(5).EQ.1)
     .  WRITE(0,20)" - CMS charg(neutal)ino                           "
       IF(GMUFLAG.EQ.1)
     .  WRITE(0,20)" - (g-2)_muon                                     "
      ENDIF

      IF(TOT.GT.0)THEN

       WRITE(0,*)
       WRITE(0,*)
       WRITE(0,20)"Parameter ranges for good points:                 "
       WRITE(0,*)
       WRITE(0,30)" LAMBDA: ",LN,LNN,DEV(LN,LNN)
       WRITE(0,30)" TANB: ",TBN,TBNN,DEV(TBN,TBNN)
       WRITE(0,30)" M0: ",M0N,M0NN,DEV(M0N,M0NN)
       IF(M1FLAG*M2FLAG*M3FLAG.EQ.0)
     .  WRITE(0,30)" M12: ",M12N,M12NN,DEV(M12N,M12NN)
       IF(M1FLAG.EQ.1)WRITE(0,30)"M1: ",M1N,M1NN,DEV(M1N,M1NN)
       IF(M2FLAG.EQ.1)WRITE(0,30)"M2: ",M2N,M2NN,DEV(M2N,M2NN)
       IF(M3FLAG.EQ.1)WRITE(0,30)"M3: ",M3N,M3NN,DEV(M3N,M3NN)
       WRITE(0,30)" A0: ",A0N,A0NN,DEV(A0N,A0NN)
       IF(ALFLAG.EQ.1)WRITE(0,30)"ALAMBDA: ",ALN,ALNN,DEV(ALN,ALNN)
       IF(AKFLAG.EQ.1)WRITE(0,30)"AKAPPA: ",AKN,AKNN,DEV(AKN,AKNN)
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
        IF(MSFLAG.EQ.1)
     .   WRITE(0,30)" MS: ",MSN,MSNN,DEV(MSN,MSNN)
        WRITE(0,30)" XIS: ",XISN,XISNN,DEV(XISN,XISNN)
        WRITE(0,40)"(XIS is not an input parameter)"
       ENDIF
       IF(MAFLAG.EQ.-5)THEN
        WRITE(0,30)"KAPPA: ",KN,KNN,DEV(KN,KNN)
        IF(XIFN.NE.0d0 .OR. XIFNN.NE.0d0)
     .   WRITE(0,30)" XIF: ",XIFN,XIFNN,DEV(XIFN,XIFNN)
        IF(XISN.NE.0d0 .OR. XISNN.NE.0d0)
     .   WRITE(0,30)" XIS: ",XISN,XISNN,DEV(XISN,XISNN)
        WRITE(0,30)"MHD: ",MHDN,MHDNN,DEV(MHDN,MHDNN)
        WRITE(0,40)"(MHD is not an input parameter)"
        WRITE(0,30)"MHU: ",MHUN,MHUNN,DEV(MHUN,MHUNN)
        WRITE(0,40)"(MHU is not an input parameter)"
        WRITE(0,30)" MS: ",MSN,MSNN,DEV(MSN,MSNN)
        WRITE(0,40)"(MS is not an input parameter)"
        WRITE(0,30)"MUEFF: ",MUN,MUNN,DEV(MUN,MUNN)
       ELSE
        WRITE(0,30)"MUEFF: ",MUN,MUNN,DEV(MUN,MUNN)
        WRITE(0,40)"(MUEFF is not an input parameter)"
        IF(MHDFLAG.EQ.1)
     .   WRITE(0,30)"MHD: ",MHDN,MHDNN,DEV(MHDN,MHDNN)
        IF(MHUFLAG.EQ.1)
     .   WRITE(0,30)"MHU: ",MHUN,MHUNN,DEV(MHUN,MHUNN)
       ENDIF
       IF(MUPN.NE.0d0 .OR. MUPNN.NE.0d0)
     .  WRITE(0,30)" MUP: ",MUPN,MUPNN,DEV(MUPN,MUPNN)
       IF(MSPN.NE.0d0 .OR. MSPNN.NE.0d0)
     .  WRITE(0,30)" MSP: ",MSPN,MSPNN,DEV(MSPN,MSPNN)
       IF(M3HN.NE.0d0 .OR. M3HNN.NE.0d0)
     .  WRITE(0,30)" M3H: ",M3HN,M3HNN,DEV(M3HN,M3HNN)

      ENDIF

 10   FORMAT(A35,I10)
 20   FORMAT(A50)
 30   FORMAT(A11,3E15.4)
 40   FORMAT(A36)

      END
