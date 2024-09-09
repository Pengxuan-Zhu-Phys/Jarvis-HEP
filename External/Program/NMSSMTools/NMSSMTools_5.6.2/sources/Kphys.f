      SUBROUTINE KPHYS(PAR,PROB)

*      Subroutine for K Physics observables

* Literature:
*
*   - THEORETICAL FORMULAE:
* [1] A.J.Buras, T.Ewerth, S.Jager and J.Rosiek,
*     `K+ ---> pi+ nu anti-nu and K(L) ---> pi0 nu anti-nu decays in the general MSSM,'
*     Nucl.\ Phys.\ B {\bf 714} (2005) 103, ePrint: [hep-ph/0408142].
*
* [2] A.J.Buras and J.Girrbach,
*     `Towards the Identification of New Physics through Quark Flavour Violating Processes,'
*     Rept.\ Prog.\ Phys.\  {\bf 77} (2014) 086201, ePrint: [[arXiv:1306.3775 [hep-ph]].
*  
* [3] A.J.Buras, D.Buttazzo, J.Girrbach-Noe and R.Knegjens,
*     `$K^+\to\pi^+\nu\bar\nu$ and $K_L\to\pi^0\nu\bar\nu$ in the Standard Model: 
*     Status and Perspectives,' ePrint: arXiv:1503.02693 [hep-ph].
*
* [4] A. J. Buras, P. H. Chankowski, J. Rosiek, L. Slawianowska,
*     ' Delta M(d, s), B0(d, s) ---> mu+ mu- and B ---> X(s) gamma
*      in supersymmetry at large tan beta.'
*     Nucl.Phys.B659:3,2003, e-Print: hep-ph/0210145
*
* [5] A.J.Buras, S.Jager and J.Urban,
*    `Master formulae for Delta F=2 NLO QCD factors in the standard model and 
*    beyond,' Nucl.Phys.B 605 (2001) 600
*    e-Print: [hep-ph/0102316].
*
* [6] J.Brod and M.Gorbahn,
*     `Next-to-Next-to-Leading-Order Charm-Quark Contribution to the CP Violation Parameter epsilon_K and Delta M_K,'
*     Phys.Rev.Lett. 108 (2012) 121801, ePrint: [arXiv:1108.2036 [hep-ph]].
*
* [7] J.A.Bailey [MILC s Collaboration],
*     `The $B\to D\ell\nu$ form factors at nonzero recoil and |V$_{cb}$| from 2+1-flavor lattice QCD,'
*     Phys. Rev. D {\bf 92} (2015) 3,  034506, e-Print: [arXiv:1503.07237 [hep-lat]].
*
*   - SOURCES FOR EXPERIMENTAL/LATTICE QCD DATA:
* [1'] K.A. Olive et al. (Particle Data Group), Chin. Phys. C, 38, 090001 (2014).
*
* [2'] A.V.Artamonov et al. [E949 Collaboration],
*     `New measurement of the $K^{+} \to \pi^{+} \nu \bar{\nu}$ branching ratio,'
*     Phys.\ Rev.\ Lett.\  {\bf 101} (2008) 191802, ePrint: arXiv:0808.2459 [hep-ex]].
*
* [3'] J.K.Ahn et al. [E391a Collaboration],
*      `Experimental study of the decay K0(L) ---> pi0 nu nu-bar,'
*      Phys.\ Rev.\ D {\bf 81} (2010) 072004, ePrint: [arXiv:0911.4789 [hep-ex]].
*
* [4'] Y.C.Jang {\it et al.} [SWME Collaboration],
*      `Kaon BSM B-parameters using improved staggered fermions from $N_f=2+1$ unquenched QCD,'
*      ePrint: arXiv:1509.00592 [hep-lat].
*
* [5'] S.Aoki et al.,
*      `Review of lattice results concerning low-energy particle physics,'
*      Eur. Phys. J. C 74 (2014) 2890,  e-Print: [arXiv:1310.8555 [hep-lat]]. 
*      updates at http://itpwiki.unibe.ch/flag/
*

      IMPLICIT NONE

      INTEGER I,J,K,L
      DOUBLE PRECISION PAR(*),PROB(*)
      DOUBLE PRECISION Pi,aux,auxI,scR,VVc,VVcI,VVu,VVuI
      DOUBLE PRECISION MD0,MS0,sc0,scL,asc0,ascL,asmt,asmh,mth,mc_h
      DOUBLE PRECISION MT0,xt,yt,z,TANB,CCT(2,2),CCC(2,2),CCL(2,2)
      DOUBLE PRECISION MC0,AC,SC(2),RSC(2,2),thetaSC,MSF2_11,MSF2_22
      DOUBLE PRECISION SST,ST(2),RST(2,2),SL(2),RSL(2,2),SSL,MSF4_12
      DOUBLE PRECISION fh10,fh20,fh11,fh61,fc30,fc40,fc50,asf,runmass
      DOUBLE PRECISION CLSM,CLH,CLeCHAR,CLtauCHAR,CLeCHARI,CLtauCHARI
      DOUBLE PRECISION CLe,CLtau,lamb,Kp,PcX,DEM,sgn
      DOUBLE PRECISION RVtdVts,IVtdVts,RVcdVcs,IVcdVcs,RVudVus
      DOUBLE PRECISION BRKp_PipnunubexpMax,BRKL_Pi0nunubexpMax,KL
      DOUBLE PRECISION MCHH(2),D0B,D2B,CCD(2),etaH,etaS,S0,S02,xc
      DOUBLE PRECISION PRLH_ts(2),PRLH_td(2),PLRH_ts(2),PLRH_td(2)
      DOUBLE PRECISION PRLH_cs(2),PRLH_cd(2),PLRH_cs(2),PLRH_cd(2)
      DOUBLE PRECISION PLRH_us(2),PLRH_ud(2),sigRLsd,sigLRsd,BB0,BB1
      DOUBLE PRECISION CVLLHIG,C1SLLHIG,C2SLLHIG,C1LRHIG,C2LRHIG
      DOUBLE PRECISION CVLLHIGI,C1SLLHIGI,C2SLLHIGI,C1LRHIGI,C2LRHIGI
      DOUBLE PRECISION CVLLCHAR,C1SLLCHAR,C2SLLCHAR,C1LRCHAR,C2LRCHAR
      DOUBLE PRECISION CVLLCHARI,C1SLLCHARI,C2SLLCHARI,C1LRCHARI
      DOUBLE PRECISION C2LRCHARI,C1SLLDPH,C2LRDPH,C1SLLDPHI,C2LRDPHI
      DOUBLE PRECISION CVLL,C1SLL,C2SLL,C1LR,C2LR,eta5,eta4,MK0,FK0
      DOUBLE PRECISION CVLLI,C1SLLI,C2SLLI,C1LRI,C2LRI
      DOUBLE PRECISION BVLL,B1SLL,B2SLL,B1LR,BK,B2LR,PVLL,P1SLL,P2SLL
      DOUBLE PRECISION P1LR,P2LR,etacc,etact,etatt,DMKSM,epsKSM
      DOUBLE PRECISION DMKexpmin,DMKexpMax,epsKexpmin,epsKexpMax

      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION MS,MC,MBNP,MB,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2)
      DOUBLE PRECISION PCOMP(2,2),CMASS
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      DOUBLE PRECISION HTQ,HBQ,MTOPQ,MBOTQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION eps0,epst0,epst1,epst2,epst3,epsts,epstb
      DOUBLE PRECISION epsY32,epsY31,epsY23,epsY13,epscs,epscb
      DOUBLE PRECISION BRKp_Pipnunub,BRKp_Pipnunubmin,BRKp_PipnunubMax,
     .       BRKL_Pi0nunub,BRKL_Pi0nunubmin,BRKL_Pi0nunubMax,
     .       DMK,DMKmin,DMKmax,epsK,epsKmin,epsKmax
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/STSBSCALE/QSTSB
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/SMSPEC/MS,MC,MBNP,MB,MT,MTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,N
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/FLAV3/BRKp_Pipnunub,BRKp_Pipnunubmin,BRKp_PipnunubMax,
     .             BRKL_Pi0nunub,BRKL_Pi0nunubmin,BRKL_Pi0nunubMax,
     .             DMK,DMKmin,DMKmax,epsK,epsKmin,epsKmax
      COMMON/EPSCOUP/eps0,epst0,epst1,epst2,epst3,epsts,epstb,
     .               epsY32,epsY31,epsY23,epsY13,epscs,epscb
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH

************************************************************************

*	I- Parameters
      pi=4d0*datan(1d0)
      
      TANB=PAR(3)
      
      VVu=-563.043d0                               ! Re[V_ud.V_us^*/V_td.V_ts^*]
      VVuI=-237.432d0                              ! Im[V_ud.V_us^*/V_td.V_ts^*]
      VVc=-VVu-1d0                                 ! Re[V_cd.V_cs^*/V_td.V_ts^*]
      VVcI=-VVuI                                   ! Im[V_cd.V_cs^*/V_td.V_ts^*]

*   m_top(m_top)(MSbar)
      asmt=asf(MT) !0.108002d0
      MT0=MT/(1d0+4d0/(3d0*pi)*asmt+11d0/pi**2*asmt**2) !165d0
      sc0=MT0

*  Diagonalizing the scharm mass-matrices
      scR=dsqrt(QSTSB)
      MC0=runmass(1.25d0,scR)
      AC=PAR(12)

      MSF2_11=MUL**2+MC0**2                           ! scharm
      MSF2_22=MUR**2+MC0**2
      MSF4_12=MC0**2*(AC-muq/tanbq)**2

      SC(1)=(MSF2_11+MSF2_22-dsqrt((MSF2_11-MSF2_22)**2+4d0*MSF4_12))
     .                                                          /2d0
      SC(1)=dsqrt(SC(1))
      SC(2)=(MSF2_11+MSF2_22+dsqrt((MSF2_11-MSF2_22)**2+4d0*MSF4_12))
     .                                                          /2d0
      SC(2)=dsqrt(SC(2))

      aux=MSF2_11-MSF2_22+dsqrt((MSF2_11-MSF2_22)**2+4d0*MSF4_12)
      IF(dsqrt(4d0*MSF4_12).ge.dabs(aux)*1.d-10)THEN
       thetaSC=datan(aux/dsqrt(4d0*MSF4_12))
       IF(dsqrt(4d0*MSF4_12).le.0.d0)thetaSC=thetaSC+Pi
      elseif(aux.ge.0.d0)THEN
       thetaSC=Pi/2.d0
      else
       thetaSC=-Pi/2.d0
      endif

      RSC(1,1)=DDCOS(thetaSC)
      RSC(1,2)=-DDSIN(thetaSC)
      RSC(2,1)=DDSIN(thetaSC)
      RSC(2,2)=DDCOS(thetaSC)
  
*   Stop Masses and Mixing Angles
      ST(1)=MST1
      ST(2)=MST2
      SST=DSQRT(1d0-CST**2)
      RST(1,1)=CST
      RST(2,2)=CST
      RST(1,2)=SST
      RST(2,1)=-SST

*   Stau Masses and Mixing Angles
      SL(1)=MSL1
      SL(2)=MSL2
      SSL=DSQRT(1d0-CSL**2)
      RSL(1,1)=CSL
      RSL(2,2)=CSL
      RSL(1,2)=SSL
      RSL(2,1)=-SSL
      
      MS0=runmass(0.095d0,scR)
      MD0=runmass(0.005d0,scR)
      
*   Experimental limits
*     - BR[K+ -> Pi+ nu nubar]                 [2']
      BRKp_PipnunubexpMax=4.03d-10

*     - BR[KL -> Pi0 nu nubar]                 [3']
      BRKL_Pi0nunubexpMax=2.6d-8
      
*     - DMK=(0.5293 +/- 0.0009).10^-2ps^-1
      DMKexpmin=0.5293d-2-2d0*0.0009d-2
      DMKexpMax=0.5293d-2+2d0*0.0009d-2
      
*     - epsK=(2.228 +/- 0.011).10^-3
      epsKexpmin=2.228d-3-2d0*0.011d-3
      epsKexpMax=2.228d-3+2d0*0.011d-3


*	II- K -> Pi nu nubar

*        1) Wilson coefficients [1]

*  SM contribution
      xt=(MT0/MW)**2
      CLSM=(fh10(xt)+asf(sc0)/4d0/Pi*fh11(xt))/4d0         ! Z-penguin
     . -(fh20(xt)+asf(sc0)/4d0/Pi*fh61(xt))                ! Box
      CLSM=0.994d0*CLSM                                    ! mt(mt) correction
      CLSM=1.481d0                                         ! [3] LO + QCD NLO + EW

*  Charged Higgs contribution
      yt=(MT0/CMASS)**2
      CLH=-xt*fh20(yt)/8d0/tanb**2                  ! Z-penguin

*  SUSY contribution
      DO J=1,2
      DO K=1,2
       CCT(j,k)=V(j,1)*RST(k,1)-MT0/H1Q/dsqrt(g2q)*V(j,2)*RST(k,2)
       CCC(j,k)=V(j,1)*RSC(k,1)-MC0/H1Q/dsqrt(g2q)*V(j,2)*RSC(k,2)
       CCL(j,k)=U(j,1)*RSL(k,1)-MTAU/H2Q/dsqrt(g2q)*U(j,2)*RSL(k,2)
      ENDDO
      ENDDO

      CLeCHAR=0d0
      CLeCHARI=0d0

      DO I=1,2
      DO J=1,2
      DO K=1,2                                            ! Z-penguin
       z=(MCH(J)/MCH(I))**2
       yt=(ST(K)/MCH(I))**2

       CLeCHAR=CLeCHAR+CCT(J,K)*CCT(I,K)        ! 1 stop / 2 charginos
     .     *(2d0*MCH(J)/MCH(I)*U(J,1)*U(I,1)*fc30(z,yt)
     .                        -V(J,1)*V(I,1)*fc40(z,yt))
       z=(ST(J)/MCH(I))**2
       yt=(ST(K)/MCH(I))**2                     ! 2 stops / 1 chargino
       CLeCHAR=CLeCHAR
     .         +CCT(I,J)*CCT(I,K)*RST(K,1)*RST(J,1)*fc40(yt,z)

       z=(MCH(J)/MCH(I))**2
       yt=(SC(K)/MCH(I))**2

       CLeCHAR=CLeCHAR+VVc*CCC(J,K)*CCC(I,K)    ! 1 scharm / 2 charginos
     .     *(2d0*MCH(J)/MCH(I)*U(J,1)*U(I,1)*fc30(z,yt)
     .                        -V(J,1)*V(I,1)*fc40(z,yt))
       CLeCHARI=CLeCHARI+VVcI*CCC(J,K)*CCC(I,K)  ! 1 scharm / 2 charginos
     .     *(2d0*MCH(J)/MCH(I)*U(J,1)*U(I,1)*fc30(z,yt)
     .                        -V(J,1)*V(I,1)*fc40(z,yt))

       z=(SC(J)/MCH(I))**2
       yt=(SC(K)/MCH(I))**2                     ! 2 stops / 1 chargino
       CLeCHAR=CLeCHAR
     .     +VVc*CCC(I,J)*CCC(I,K)*RSC(K,1)*RSC(J,1)*fc40(yt,z)
       CLeCHARI=CLeCHARI
     .     +VVcI*CCC(I,J)*CCC(I,K)*RSC(K,1)*RSC(J,1)*fc40(yt,z)
      ENDDO

       z=(MCH(J)/MCH(I))**2
       yt=(MUL/MCH(I))**2                       ! 1 sup / 2charginos
       CLeCHAR=CLeCHAR+VVu*V(I,1)*V(J,1)
     .                 *(2d0*MCH(J)/MCH(I)*U(J,1)*U(I,1)
     .                       *fc30(z,yt)-V(J,1)*V(I,1)*fc40(z,yt))
       CLeCHARI=CLeCHARI+VVuI*V(I,1)*V(J,1)
     .                 *(2d0*MCH(J)/MCH(I)*U(J,1)*U(I,1)
     .                       *fc30(z,yt)-V(J,1)*V(I,1)*fc40(z,yt))
      ENDDO

       z=(MUL/MCH(I))**2                       ! 2 sups / 1 chargino
       CLeCHAR=CLeCHAR+VVu*V(I,1)**2*fc40(z,z)
       CLeCHARI=CLeCHARI+VVuI*V(I,1)**2*fc40(z,z)
      ENDDO
      CLeCHAR=-CLeCHAR/8d0
      CLeCHARI=-CLeCHARI/8d0
      CLtauCHAR=CLeCHAR
      CLtauCHARI=CLeCHARI

      aux=0d0
      auxI=0d0

      DO I=1,2                                            ! Box
       xt=(MLL/MCH(I))**2
      DO J=1,2
       z=(MCH(J)/MCH(I))**2
      DO K=1,2
       yt=(ST(K)/MCH(I))**2                     ! 2 stops / 2 charginos

       aux=aux+U(J,1)*U(I,1)*CCT(J,K)*CCT(I,K)/MCH(I)**2
     .                                    *fc50(z,yt,xt)
     
       yt=(SC(K)/MCH(I))**2                     ! 2 scharms / 2 charginos

       aux=aux+VVc*U(J,1)*U(I,1)*CCC(J,K)*CCC(I,K)/MCH(I)**2
     .                                        *fc50(z,yt,xt)
       auxI=auxI+VVcI*U(J,1)*U(I,1)*CCC(J,K)*CCC(I,K)/MCH(I)**2
     .                                        *fc50(z,yt,xt)

      ENDDO
       yt=(MUL/MCH(I))**2                       ! 2 sups / 2 charginos

       aux=aux+VVu*U(J,1)*U(I,1)*V(J,1)*V(I,1)/MCH(I)**2
     .                                     *fc50(z,yt,xt)
       auxI=auxI+VVuI*U(J,1)*U(I,1)*V(J,1)*V(I,1)/MCH(I)**2
     .                                     *fc50(z,yt,xt)

      ENDDO
      ENDDO

      CLeCHAR=CLeCHAR+MW**2/4d0*aux
      CLeCHARI=CLeCHARI+MW**2/4d0*auxI

      aux=0d0
      auxI=0d0

      DO I=1,2                                            ! Box
      DO J=1,2
       z=(MCH(J)/MCH(I))**2
      DO L=1,2
       xt=(SL(L)/MCH(I))**2
      DO K=1,2
       yt=(ST(K)/MCH(I))**2                     ! 2 stops / 2 charginos

       aux=aux+CCL(J,L)*CCL(I,L)*CCT(J,K)*CCT(I,K)/MCH(I)**2
     .                                         *fc50(z,yt,xt)

       yt=(SC(K)/MCH(I))**2                     ! 2 scharms / 2 charginos

       aux=aux+VVc*CCL(J,L)*CCL(I,L)*CCC(J,K)*CCC(I,K)/MCH(I)**2
     .                                         *fc50(z,yt,xt)
       auxI=auxI+VVcI*CCL(J,L)*CCL(I,L)*CCC(J,K)*CCC(I,K)/MCH(I)**2
     .                                         *fc50(z,yt,xt)

      ENDDO
       yt=(MUL/MCH(I))**2                       ! 2 sups / 2 charginos

       aux=aux+VVu*V(J,1)*V(I,1)/MCH(I)**2*CCL(J,L)*CCL(I,L)
     .                                      *fc50(z,yt,xt)
       auxI=auxI+VVuI*V(J,1)*V(I,1)/MCH(I)**2*CCL(J,L)*CCL(I,L)
     .                                      *fc50(z,yt,xt)

      ENDDO
      ENDDO
      ENDDO

      CLtauCHAR=CLtauCHAR+MW**2/4d0*aux
      CLtauCHARI=CLtauCHARI+MW**2/4d0*auxI

*   Summary
      CLe=CLSM+CLH+CLeCHAR
      CLtau=CLSM+CLH+CLtauCHAR

*	 2) The branching ratio BR[K+ -> Pi+ nu nu] [2]
* Form factor / charm / EM / CKM parameters [3]
      lamb=0.22537d0
      Kp=5.173d-11*(lamb/0.225d0)**8
      PcX=0.404d0
      DEM=-0.003d0
      RVtdVts=-3.311d-4                 ! [1']
      IVtdVts=1.357d-4
      RVcdVcs=-0.219d0
      IVcdVcs=-IVtdVts
      RVudVus=-RVtdVts-RVcdVcs

* Branching ratio
      aux=((IVtdVts*CLe+RVtdVts*CLeCHARI)/lamb**5)**2
     . +(RVcdVcs*PcX/lamb+(RVtdVts*Cle-IVtdVts*CLeCHARI)/lamb**5)**2
      auxi=((IVtdVts*Cltau+RVtdVts*CltauCHARI)/lamb**5)**2
     . +(RVcdVcs*PcX/lamb+(RVtdVts*Cltau-IVtdVts*CLtauCHARI)/lamb**5)**2
     
      BRKp_Pipnunub=Kp*(1d0+DEM)*(2d0*aux+auxi)/3d0

* Error estimate
      aux=2d0*dabs(IVtdVts*CLSM)/lamb**10
     .  *(dabs(IVtdVts*CLeCHAR+RVtdVts*CLeCHARI)+dabs(IVtdVts*CLH))
     . +2d0*dabs(RVcdVcs*PcX/lamb+RVtdVts*ClSM/lamb**5)
     .  *(dabs(RVtdVts*CleCHAR-IVtdVts*CLeCHARI)+dabs(RVtdVts*CLH))
     .                                                     /lamb**5
      auxi=2d0*dabs(IVtdVts*CLSM)/lamb**10
     .  *(dabs(IVtdVts*CLtauCHAR+RVtdVts*CLtauCHARI)+dabs(IVtdVts*CLH))
     . +2d0*dabs(RVcdVcs*PcX/lamb+RVtdVts*ClSM/lamb**5)
     .  *(dabs(RVtdVts*CltauCHAR-IVtdVts*CLtauCHARI)+dabs(RVtdVts*CLH))
     .                                                         /lamb**5

      BRKp_Pipnunubmin=BRKp_Pipnunub-2d0*1d-11          ! SM [3]
     .  -0.3d0*Kp*(1d0+DEM)*(2d0*aux+auxi)/3d0          ! 30% NP
      BRKp_PipnunubMax=BRKp_Pipnunub+2d0*1d-11          ! SM [3]
     .  +0.3d0*Kp*(1d0+DEM)*(2d0*aux+auxi)/3d0          ! 30% NP
!      print*,BRKp_Pipnunubmin,BRKp_Pipnunub,BRKp_PipnunubMax

*      Comparison with experimental data (source [2'])
      prob(59)=0d0

      IF(BRKp_Pipnunubmin.GE.BRKp_PipnunubexpMax)
     .     PROB(59)=PROB(59)+BRKp_Pipnunubmin/BRKp_PipnunubexpMax-1d0

*	 3) The branching ratio BR[KL -> Pi0 nu nu] [2]
* Form factor [3]
      KL=2.231d-10*(lamb/0.225d0)**8

* Branching ratio
      aux=((IVtdVts*CLe+RVtdVts*CLeCHARI)/lamb**5)**2
      auxi=((IVtdVts*Cltau+RVtdVts*CltauCHARI)/lamb**5)**2
     
      BRKL_Pi0nunub=KL*(2d0*aux+auxi)/3d0

* Error estimate
      aux=2d0*dabs(IVtdVts*CLSM)/lamb**10
     .  *(dabs(IVtdVts*CLeCHAR+RVtdVts*CLeCHARI)+dabs(IVtdVts*CLH))
      auxi=2d0*dabs(IVtdVts*CLSM)/lamb**10
     .  *(dabs(IVtdVts*CLtauCHAR+RVtdVts*CLtauCHARI)+dabs(IVtdVts*CLH))

      BRKL_Pi0nunubmin=BRKL_Pi0nunub-2d0*0.6d-11          ! SM [3]
     .  -0.3d0*KL*(2d0*aux+auxi)/3d0                      ! 30% NP
      BRKL_Pi0nunubMax=BRKL_Pi0nunub+2d0*0.6d-11          ! SM [3]
     .  +0.3d0*KL*(2d0*aux+auxi)/3d0                      ! 30% NP
!      print*,BRKL_Pi0nunubmin,BRKL_Pi0nunub,BRKL_Pi0nunubMax

*      Comparison with experimental data (source [2'])
!      prob(59)=0d0

      IF(BRKL_Pi0nunubmin.GE.BRKL_Pi0nunubexpMax)
     .     PROB(59)=PROB(59)+BRKL_Pi0nunubmin/BRKL_Pi0nunubexpMax-1d0


*	III- K - Kbar mixing

      MK0=0.497611d0                              ! [1']
      FK0=0.1563d0                                ! [5']

*  Matching scale
      sc0=166d0
      asc0=asf(sc0)

*  Low energy scale
      scL=3d0
      ascL=asf(scL)

*	 1) Wilson coefficients for DMs at the matching scale - New-Physics only
*          - Standard Model [1] Eq.(6.7)
      
*          - Charged Higgs Boxes [4] App.(A.4.i)
*  Coefficients at CMASS
      MCHH(1)=MW
      MCHH(2)=CMASS
      
      asmh=asf(CMASS)
      mth=mt*(asmh/asmt)**(4d0/7d0)/(1d0+4d0/(3d0*pi)*asmt)
      mc_h=runmass(1.25d0,CMASS)

      aux=mth/dsqrt(h1q**2+h2q**2)
      PRLH_ts(1)=aux
      PRLH_ts(2)=aux/tanb*(1d0-(1d0/tanb+tanb)
     .            *(epsts-epsY32*(epstb-epsts)/(1d0+eps0)))
      PRLH_td(1)=aux
      PRLH_td(2)=aux/tanb*(1d0-(1d0/tanb+tanb)
     .            *(epsts-epsY31*(epstb-epsts)/(1d0+epst0)))

      aux=mc_h/dsqrt(h1q**2+h2q**2)
      PRLH_cs(1)=aux
      PRLH_cs(2)=aux/tanb*(1d0-(1d0/tanb+tanb)*epscs)
      PRLH_cd(1)=aux
      PRLH_cd(2)=aux/tanb*(1d0-(1d0/tanb+tanb)*epscs)

      aux=runmass(0.095d0,CMASS)/dsqrt(h1q**2+h2q**2)
      PLRH_ts(1)=-aux
      PLRH_ts(2)=aux*((1d0/tanb+tanb)/(1d0+epst2)
     .    *(1d0-epsY23/(1d0+epst2)/(1d0+eps0)
     .         -epsY32*(epst2-epst3)/(1d0+epst3)/(1d0+eps0))-1d0)

      PLRH_cs(1)=-aux
      PLRH_cs(2)=aux*((1d0/tanb+tanb)/(1d0+epst2)-1d0)

      PLRH_us(1)=-aux
      PLRH_us(2)=aux*((1d0/tanb+tanb)/(1d0+epst1)-1d0)
     
      aux=runmass(0.005d0,CMASS)/dsqrt(h1q**2+h2q**2)
      PLRH_td(1)=-aux
      PLRH_td(2)=aux*((1d0/tanb+tanb)/(1d0+epst1)
     .    *(1d0-epsY13/(1d0+epst1)/(1d0+eps0)
     .         -epsY31*(epst1-epst3)/(1d0+epst3)/(1d0+epst0))-1d0)

      PLRH_cd(1)=-aux
      PLRH_cd(2)=aux*((1d0/tanb+tanb)/(1d0+epst1)-1d0)

      PLRH_ud(1)=-aux
      PLRH_ud(2)=aux*((1d0/tanb+tanb)/(1d0+epst1)-1d0)

!      CVLLHIG=-g2q/2d0*PRLH_ts(2)*PRLH_td(2)
!     .                          *mth**2*D0B(MW,MCHH(2),mth,mth)
!     . +(PRLH_ts(2)*PRLH_td(2))**2*D2B(MCHH(2),MCHH(2),mth,mth)/8d0
!     . +PRLH_ts(2)*PRLH_td(2)*PRLH_ts(1)*PRLH_td(1)
!     .                            *D2B(MCHH(1),MCHH(2),mth,mth)/4d0
!      CVLLHIG=CVLLHIG/(GF*MW)**2

!      C1SLLHIG=0d0
!      DO I=1,2
!      DO J=1,2
!       C1SLLHIG=C1SLLHIG+PLRH_ts(J)*PLRH_ts(I)*PRLH_td(I)*PRLH_td(J)
!     .                   *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)/2d0
!      ENDDO
!      ENDDO
!      C1SLLHIG=C1SLLHIG/(GF*MW)**2

!      aux=1d0
!      C1LRHIG=0d0
!      DO I=1,2
!      DO J=1,2
!       C1LRHIG=C1LRHIG+PRLH_ts(J)*PRLH_td(I)*PLRH_ts(I)*PLRH_td(J)
!     .  *(D2B(MCHH(I),MCHH(J),mth,mth)-aux*D2B(MCHH(I),MCHH(J),mth,0d0))
!     .                                               /4d0
!      ENDDO
!      ENDDO
!      C1LRHIG=C1LRHIG/(GF*MW)**2

!      aux=1d0
!      C2LRHIG=0d0
!      DO I=1,2
!       C2LRHIG=C2LRHIG-g2q*PLRH_ts(I)*PLRH_td(I)/2d0
!     .    *(D2B(MCHH(I),MW,mth,mth)
!     .      -2d0*aux*D2B(MCHH(I),MW,mth,0d0)
!     .      +aux**2*D2B(MCHH(I),MW,0d0,0d0))
!      DO J=1,2
!       C2LRHIG=C2LRHIG+PLRH_ts(J)*PRLH_ts(I)*PRLH_td(I)*PLRH_td(J)
!     .                     *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)
!      ENDDO
!      ENDDO
!      C2LRHIG=C2LRHIG/(GF*MW)**2

      aux=-g2q/2d0*PRLH_ts(2)*PRLH_td(2)
     .                          *mth**2*D0B(MW,MCHH(2),mth,mth)
     . +(PRLH_ts(2)*PRLH_td(2))**2*D2B(MCHH(2),MCHH(2),mth,mth)/8d0
     . +PRLH_ts(2)*PRLH_td(2)*PRLH_ts(1)*PRLH_td(1)
     .                            *D2B(MCHH(1),MCHH(2),mth,mth)/4d0
      CVLLHIG=aux*(RVtdVts**2-IVtdVts**2)
      CVLLHIGI=-2d0*aux*RVtdVts*IVtdVts
      
      aux=-g2q/2d0*(PRLH_cs(2)*PRLH_td(2)+PRLH_cd(2)*PRLH_ts(2))
     .                       *mth*mc_h*D0B(MW,MCHH(2),mth,mc_h)
     . +2d0*PRLH_ts(2)*PRLH_td(2)*PRLH_cs(2)*PRLH_cd(2)
     .                            *D2B(MCHH(2),MCHH(2),mth,mc_h)/8d0
     . +(PRLH_ts(2)*PRLH_cd(2)*PRLH_cs(1)*PRLH_td(1)
     .               +PRLH_cs(2)*PRLH_td(2)*PRLH_ts(1)*PRLH_cd(1))
     .                            *D2B(MCHH(1),MCHH(2),mth,mc_h)/4d0
      CVLLHIG=CVLLHIG+aux*(RVtdVts*RVcdVcs-IVtdVts*IVcdVcs)
      CVLLHIGI=CVLLHIGI-aux*(RVtdVts*IVcdVcs+IVtdVts*RVcdVcs)

      aux=-g2q/2d0*PRLH_cs(2)*PRLH_cd(2)
     .                          *mc_h**2*D0B(MW,MCHH(2),mc_h,mc_h)
     . +(PRLH_cs(2)*PRLH_cd(2))**2*D2B(MCHH(2),MCHH(2),mc_h,mc_h)/8d0
     . +PRLH_cs(2)*PRLH_cd(2)*PRLH_cs(1)*PRLH_cd(1)
     .                            *D2B(MCHH(1),MCHH(2),mc_h,mc_h)/4d0
      CVLLHIG=CVLLHIG+aux*(RVcdVcs**2-IVcdVcs**2)
      CVLLHIGI=CVLLHIGI-2d0*aux*RVcdVcs*IVcdVcs
      
      CVLLHIG=CVLLHIG/(GF*MW)**2
      CVLLHIGI=CVLLHIGI/(GF*MW)**2

      aux=0d0
      DO I=1,2
      DO J=1,2
       aux=aux+PLRH_ts(J)*PLRH_ts(I)*PRLH_td(I)*PRLH_td(J)
     .                   *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)/2d0
      ENDDO
      ENDDO
      C1SLLHIG=aux*(RVtdVts**2-IVtdVts**2)
      C1SLLHIGI=-2d0*aux*RVtdVts*IVtdVts

      aux=0d0
      DO I=1,2
      DO J=1,2
       aux=aux+PLRH_ts(J)*PLRH_cs(I)*PRLH_td(I)*PRLH_cd(J)
     .                   *mth*mc_h*D0B(MCHH(I),MCHH(J),mth,mc_h)
      ENDDO
      ENDDO
      C1SLLHIG=C1SLLHIG+aux*(RVtdVts*RVcdVcs-IVtdVts*IVcdVcs)
      C1SLLHIGI=C1SLLHIGI-aux*(RVcdVcs*IVtdVts+IVcdVcs*RVtdVts)

      aux=0d0
      DO I=1,2
      DO J=1,2
       aux=aux+PLRH_cs(J)*PLRH_cs(I)*PRLH_cd(I)*PRLH_cd(J)
     .                *mc_h**2*D0B(MCHH(I),MCHH(J),mc_h,mc_h)/2d0
      ENDDO
      ENDDO
      C1SLLHIG=C1SLLHIG+aux*(RVcdVcs**2-IVcdVcs**2)
      C1SLLHIGI=C1SLLHIGI-2d0*aux*RVcdVcs*IVcdVcs
      
      C1SLLHIG=C1SLLHIG/(GF*MW)**2
      C1SLLHIGI=C1SLLHIGI/(GF*MW)**2

      aux=0d0
      DO I=1,2
      DO J=1,2
       aux=aux+PRLH_ts(J)*PLRH_ts(I)*PRLH_td(I)*PLRH_td(J)
     .                          *D2B(MCHH(I),MCHH(J),mth,mth)/4d0
      ENDDO
      ENDDO
      C1LRHIG=aux*(RVtdVts**2-IVtdVts**2)
      C1LRHIGI=-2d0*aux*RVtdVts*IVtdVts

      aux=0d0
      DO I=1,2
      DO J=1,2
       aux=aux+(PRLH_ts(J)*PLRH_cs(I)*PRLH_td(I)*PLRH_cd(J)
     .         +PRLH_cs(J)*PLRH_ts(I)*PRLH_cd(I)*PLRH_td(J))
     .                         *D2B(MCHH(I),MCHH(J),mth,mc_h)/4d0
      ENDDO
      ENDDO
      C1LRHIG=C1LRHIG+aux*(RVtdVts*RVcdVcs-IVtdVts*IVcdVcs)
      C1LRHIGI=C1LRHIGI-aux*(RVcdVcs*IVtdVts+RVtdVts*IVcdVcs)

      aux=0d0
      DO I=1,2
      DO J=1,2
       aux=aux+PRLH_ts(J)*PLRH_us(I)*PRLH_td(I)*PLRH_ud(J)
     .                         *D2B(MCHH(I),MCHH(J),mth,0d0)/4d0
      ENDDO
      ENDDO
      C1LRHIG=C1LRHIG+aux*RVtdVts*RVudVus
      C1LRHIGI=C1LRHIGI-aux*RVudVus*IVtdVts

      aux=0d0
      DO I=1,2
      DO J=1,2
       aux=aux+PRLH_cs(J)*PLRH_cs(I)*PRLH_cd(I)*PLRH_cd(J)
     .                        *D2B(MCHH(I),MCHH(J),mc_h,mc_h)/4d0
      ENDDO
      ENDDO
      C1LRHIG=C1LRHIG+aux*(RVcdVcs**2-IVcdVcs**2)
      C1LRHIGI=C1LRHIGI-2d0*aux*RVcdVcs*IVcdVcs

      aux=0d0
      DO I=1,2
      DO J=1,2
       aux=aux+PRLH_cs(J)*PLRH_us(I)*PRLH_cd(I)*PLRH_ud(J)
     .                         *D2B(MCHH(I),MCHH(J),mc_h,0d0)/4d0
      ENDDO
      ENDDO
      C1LRHIG=C1LRHIG+aux*RVcdVcs*RVudVus
      C1LRHIGI=C1LRHIGI-aux*RVudVus*IVcdVcs
      
      C1LRHIG=C1LRHIG/(GF*MW)**2
      C1LRHIGI=C1LRHIGI/(GF*MW)**2

      aux=0d0
      DO I=1,2
       aux=aux-g2q*PLRH_ts(I)*PLRH_td(I)/2d0
     .                       *D2B(MCHH(I),MW,mth,mth)
      DO J=1,2
       aux=aux+PLRH_ts(J)*PRLH_ts(I)*PRLH_td(I)*PLRH_td(J)
     .                     *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)
      ENDDO
      ENDDO
      C2LRHIG=aux*(RVtdVts**2-IVtdVts**2)
      C2LRHIGI=-2d0*aux*RVtdVts*IVtdVts

      aux=0d0
      DO I=1,2
       aux=aux-g2q*(PLRH_ts(I)*PLRH_cd(I)+PLRH_cs(I)*PLRH_td(I))/2d0
     .                       *D2B(MCHH(I),MW,mth,mc_h)
      DO J=1,2
       aux=aux+(PLRH_ts(J)*PRLH_cs(I)*PRLH_td(I)*PLRH_cd(J)
     .                 +PLRH_cs(J)*PRLH_ts(I)*PRLH_cd(I)*PLRH_td(J))
     .                   *mth*mc_h*D0B(MCHH(I),MCHH(J),mth,mc_h)
      ENDDO
      ENDDO
      C2LRHIG=C2LRHIG+aux*(RVtdVts*RVcdVcs-IVtdVts*IVcdVcs)
      C2LRHIGI=C2LRHIGI-aux*(RVcdVcs*IVtdVts+IVcdVcs*RVtdVts)

      aux=0d0
      DO I=1,2
       aux=aux-g2q*(PLRH_ts(I)*PLRH_ud(I)+PLRH_us(I)*PLRH_td(I))/2d0
     .                       *D2B(MCHH(I),MW,mth,0d0)
      ENDDO
      C2LRHIG=C2LRHIG+aux*RVtdVts*RVudVus
      C2LRHIGI=C2LRHIGI-aux*RVudVus*IVtdVts

      aux=0d0
      DO I=1,2
       aux=aux-g2q*PLRH_cs(I)*PLRH_cd(I)/2d0
     .                       *D2B(MCHH(I),MW,mc_h,mc_h)
      DO J=1,2
       aux=aux+PLRH_cs(J)*PRLH_cs(I)*PRLH_cd(I)*PLRH_cd(J)
     .                   *mc_h**2*D0B(MCHH(I),MCHH(J),mc_h,mc_h)
      ENDDO
      ENDDO
      C2LRHIG=C2LRHIG+aux*(RVcdVcs**2-IVcdVcs**2)
      C2LRHIGI=C2LRHIGI-2d0*aux*RVcdVcs*IVcdVcs
      
      aux=0d0
      DO I=1,2
       aux=aux-g2q*(PLRH_cs(I)*PLRH_ud(I)+PLRH_us(I)*PLRH_cd(I))/2d0
     .                       *D2B(MCHH(I),MW,mc_h,0d0)
      ENDDO
      C2LRHIG=C2LRHIG+aux*RVcdVcs*RVudVus
      C2LRHIGI=C2LRHIGI-aux*RVudVus*IVcdVcs
      
      aux=0d0
      DO I=1,2
       aux=aux-g2q*PLRH_us(I)*PLRH_ud(I)/2d0
     .                       *D2B(MCHH(I),MW,0d0,0d0)
      ENDDO
      C2LRHIG=C2LRHIG+aux*RVudVus**2

      C2LRHIG=C2LRHIG/(GF*MW)**2
      C2LRHIGI=C2LRHIGI/(GF*MW)**2

*  Running to sc0 [5] App.C
      etaH=asf(CMASS)/asc0

      CVLLHIG=etaH**(6d0/21d0)*(1d0+asc0/4d0/Pi*1.3707d0*(1d0-etaH))
     .                 *CVLLHIG

      aux=C1SLLHIG
      C1SLLHIG=(1.0153d0*etaH**(-0.6916d0)-0.0153d0*etaH**(0.7869d0)
     . +asc0/4d0/Pi*(etaH**(-0.6916d0)*(5.6478d0-6.0350d0*etaH)
     .              +etaH**(0.7869d0)*(0.3272d0+0.0600*etaH)))*aux

      C2SLLHIG=(0.0081d0*(etaH**(0.7869d0)-etaH**(-0.6916d0))
     . +asc0/4d0/Pi*(etaH**(0.7869d0)*(-0.0618d0-0.0315d0*etaH)
     .              +etaH**(-0.6916d0)*(0.0454d0+0.0479d0*etaH)))*aux

      aux=C1LRHIG
      C1LRHIG=(etaH**(3d0/21d0)
     . +asc0/4d0/Pi*(0.9219d0*etaH**(-24d0/21d0)+etaH**(3d0/21d0)
     .   *(-2.2194d0+1.2975d0*etaH)))*C1LRHIG
     . +asc0/4d0/Pi*1.3828d0*(etaH**(24d0/21d0)-etaH**(-24d0/21d0))
     .   *C2LRHIG

      C2LRHIG=(2d0/3d0*(etaH**(3d0/21d0)-etaH**(-24d0/21d0))
     . +asc0/4d0/Pi*(etaH**(-24d0/21d0)*(-6.4603d0+15.7415d0*etaH)
     .   +etaH**(3d0/21d0)*(-10.1463d0+0.8650d0*etaH)))*aux
     . +(etaH**(-24d0/21d0)+asc0/4d0/Pi*(0.9219d0*etaH**(24d0/21d0)
     .   +etaH**(-24d0/21d0)*(9.6904d0-10.6122d0*etaH)))*C2LRHIG

      CVLLHIGI=etaH**(6d0/21d0)*(1d0+asc0/4d0/Pi*1.3707d0*(1d0-etaH))
     .                 *CVLLHIGI

      aux=C1SLLHIGI
      C1SLLHIGI=(1.0153d0*etaH**(-0.6916d0)-0.0153d0*etaH**(0.7869d0)
     . +asc0/4d0/Pi*(etaH**(-0.6916d0)*(5.6478d0-6.0350d0*etaH)
     .              +etaH**(0.7869d0)*(0.3272d0+0.0600*etaH)))*aux

      C2SLLHIGI=(0.0081d0*(etaH**(0.7869d0)-etaH**(-0.6916d0))
     . +asc0/4d0/Pi*(etaH**(0.7869d0)*(-0.0618d0-0.0315d0*etaH)
     .              +etaH**(-0.6916d0)*(0.0454d0+0.0479d0*etaH)))*aux

      aux=C1LRHIGI
      C1LRHIGI=(etaH**(3d0/21d0)
     . +asc0/4d0/Pi*(0.9219d0*etaH**(-24d0/21d0)+etaH**(3d0/21d0)
     .   *(-2.2194d0+1.2975d0*etaH)))*C1LRHIGI
     . +asc0/4d0/Pi*1.3828d0*(etaH**(24d0/21d0)-etaH**(-24d0/21d0))
     .   *C2LRHIGI

      C2LRHIGI=(2d0/3d0*(etaH**(3d0/21d0)-etaH**(-24d0/21d0))
     . +asc0/4d0/Pi*(etaH**(-24d0/21d0)*(-6.4603d0+15.7415d0*etaH)
     .   +etaH**(3d0/21d0)*(-10.1463d0+0.8650d0*etaH)))*aux
     . +(etaH**(-24d0/21d0)+asc0/4d0/Pi*(0.9219d0*etaH**(24d0/21d0)
     .   +etaH**(-24d0/21d0)*(9.6904d0-10.6122d0*etaH)))*C2LRHIGI

*          - Chargino / squark Boxes [4] App.(A.4.ii)
*  Coefficients at scR
      do j=1,2
       CCD(J)=MS0/(1d0+epst2)/H2Q*U(j,2)
      do k=1,2
       CCT(j,k)=dsqrt(g2q)*V(j,1)*RST(k,1)-MT0/H1Q*V(j,2)*RST(k,2)
       CCC(j,k)=dsqrt(g2q)*V(j,1)*RSC(k,1)-MC0/H1Q*V(j,2)*RSC(k,2)
      enddo
      enddo
      
!      aux=0d0                   !3rd family contribution
!      do i=1,2
!      do j=1,2
!      do k=1,2
!      do l=1,2
!       aux=aux+CCT(j,k)*CCT(i,l)*CCT(i,k)*CCT(j,l)
!     .                   *D2B(MCH(i),MCH(j),ST(k),ST(l))
!      enddo                     !3rd family/1st-2nd interference
!       aux=aux-2d0*g2q*V(j,1)*V(i,1)*CCT(i,k)*CCT(j,k)
!     .                   *D2B(MCH(i),MCH(j),ST(k),MUL)
!      enddo                     !1st-2nd family contribution
!       aux=aux+(g2q*V(j,1)*V(i,1))**2*D2B(MCH(i),MCH(j),MUL,MUL)
!      enddo
!      enddo
!      CVLLCHAR=aux/8d0/(GF*MW)**2

!      aux=0d0                   !3rd family contribution
!      do i=1,2
!      do j=1,2
!      do k=1,2
!      do l=1,2
!       aux=aux+CCD(j)*CCD(i)*RST(k,1)*RST(l,1)*CCT(i,k)*CCT(j,l)
!     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),ST(l))
!      enddo                     !3rd family/1st-2nd interference
!       aux=aux-2d0*dsqrt(g2q)*CCD(i)*V(j,1)*CCD(j)*RST(k,1)*CCT(i,k)
!     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),MUL)
!      enddo                     !1st-2nd family contribution
!       aux=aux+g2q*V(j,1)*V(i,1)*CCD(j)*CCD(i)
!     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),MUL,MUL)
!      enddo
!      enddo
!      C1SLLCHAR=-aux/4d0/(GF*MW)**2
!      C2SLLCHAR=aux/16d0/(GF*MW)**2

!      aux=0d0                   !3rd family contribution
!      do i=1,2
!      do j=1,2
!      do k=1,2
!      do l=1,2
!       aux=aux+CCD(j)*CCD(i)*RST(l,1)**2*CCT(j,k)*CCT(i,k)
!     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),ST(l))
!      enddo                     !3rd family/1st-2nd interference
!       aux=aux-CCD(i)*CCD(j)
!     .      *(CCT(i,k)*CCT(j,k)+g2q*V(i,1)*V(j,1)*RST(k,1)**2)
!     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),MUL)
!      enddo                     !1st-2nd family contribution
!       aux=aux+g2q*V(j,1)*V(i,1)*CCD(j)*CCD(i)
!     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),MUL,MUL)
!      enddo
!      enddo
!      C1LRCHAR=-aux/2d0/(GF*MW)**2*MD0/MS0

!      aux=0d0                   !3rd family contribution
!      do i=1,2
!      do j=1,2
!      do k=1,2
!      do l=1,2
!       aux=aux+CCD(j)**2*RST(k,1)*RST(l,1)*CCT(i,k)*CCT(i,l)
!     .                   *D2B(MCH(i),MCH(j),ST(k),ST(l))
!      enddo                     !3rd family/1st-2nd interference
!       aux=aux-2d0*dsqrt(g2q)*V(i,1)*CCD(j)**2*RST(k,1)*CCT(i,k)
!     .                   *D2B(MCH(i),MCH(j),ST(k),MUL)
!      enddo                     !1st-2nd family contribution
!       aux=aux+g2q*V(i,1)**2*CCD(j)**2*D2B(MCH(i),MCH(j),MUL,MUL)
!      enddo
!      enddo
!      C2LRCHAR=-aux/2d0/(GF*MW)**2*MD0/MS0

      aux=0d0
      auxi=0d0
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+(RVtdVts**2-IVtdVts**2)    !3rd family contribution
     .          *CCT(j,k)*CCT(i,l)*CCT(i,k)*CCT(j,l)
     .                       *D2B(MCH(i),MCH(j),ST(k),ST(l))
     . +(RVcdVcs**2-IVcdVcs**2)           !2nd family contribution
     .          *CCC(j,k)*CCC(i,l)*CCC(i,k)*CCC(j,l)
     .                       *D2B(MCH(i),MCH(j),SC(k),SC(l))
     . +(RVtdVts*RVcdVcs-IVtdVts*IVcdVcs) !3rd family/2nd interference
     .          *(CCC(j,k)*CCT(i,l)*CCC(i,k)*CCT(j,l)
     .                       *D2B(MCH(i),MCH(j),SC(k),ST(l))
     .           +CCT(j,k)*CCC(i,l)*CCT(i,k)*CCC(j,l)
     .                       *D2B(MCH(i),MCH(j),ST(k),SC(l)))
       auxi=auxi-2d0*RVtdVts*IVtdVts      !3rd family contribution
     .          *CCT(j,k)*CCT(i,l)*CCT(i,k)*CCT(j,l)
     .                       *D2B(MCH(i),MCH(j),ST(k),ST(l))
     . -2d0*RVcdVcs*IVcdVcs               !2nd family contribution
     .          *CCC(j,k)*CCC(i,l)*CCC(i,k)*CCC(j,l)
     .                       *D2B(MCH(i),MCH(j),SC(k),SC(l))
     . -(IVtdVts*RVcdVcs+RVtdVts*IVcdVcs) !3rd family/2nd interference
     .          *(CCC(j,k)*CCT(i,l)*CCC(i,k)*CCT(j,l)
     .                       *D2B(MCH(i),MCH(j),SC(k),ST(l))
     .           +CCT(j,k)*CCC(i,l)*CCT(i,k)*CCC(j,l)
     .                       *D2B(MCH(i),MCH(j),ST(k),SC(l)))
      enddo                    
       aux=aux+RVtdVts*RVudVus            !3rd family/1st interference
     .  *2d0*g2q*V(j,1)*V(i,1)*CCT(i,k)*CCT(j,k)
     .                   *D2B(MCH(i),MCH(j),ST(k),MUL)
     .        +RVcdVcs*RVudVus            !2nd family/1st interference
     .  *2d0*g2q*V(j,1)*V(i,1)*CCC(i,k)*CCC(j,k)
     .                   *D2B(MCH(i),MCH(j),SC(k),MUL)
       auxi=auxi-IVtdVts*RVudVus          !3rd family/1st interference
     .  *2d0*g2q*V(j,1)*V(i,1)*CCT(i,k)*CCT(j,k)
     .                   *D2B(MCH(i),MCH(j),ST(k),MUL)
     .        -IVcdVcs*RVudVus            !2nd family/1st interference
     .  *2d0*g2q*V(j,1)*V(i,1)*CCC(i,k)*CCC(j,k)
     .                   *D2B(MCH(i),MCH(j),SC(k),MUL)
      enddo                               !1st family contribution
       aux=aux+(RVudVus*g2q*V(j,1)*V(i,1))**2
     .                   *D2B(MCH(i),MCH(j),MUL,MUL)
      enddo
      enddo
      CVLLCHAR=aux/8d0/(GF*MW)**2
      CVLLCHARI=auxi/8d0/(GF*MW)**2

      aux=0d0
      auxi=0d0
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+(RVtdVts**2-IVtdVts**2)    !3rd family contribution
     .          *CCD(j)*CCD(i)*RST(k,1)*RST(l,1)*CCT(i,k)*CCT(j,l)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),ST(l))
     . +(RVcdVcs**2-IVcdVcs**2)           !2nd family contribution
     .          *CCD(j)*CCD(i)*RSC(k,1)*RSC(l,1)*CCC(i,k)*CCC(j,l)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),SC(k),SC(l))
     . +(RVtdVts*RVcdVcs-IVtdVts*IVcdVcs) !3rd family/2nd interference
     .          *CCD(j)*CCD(i)*MCH(i)*MCH(j)
     .             *(RSC(k,1)*RST(l,1)*CCC(i,k)*CCT(j,l)
     .                      *D0B(MCH(i),MCH(j),SC(k),ST(l))
     .              +RST(k,1)*RSC(l,1)*CCT(i,k)*CCC(j,l)
     .                      *D0B(MCH(i),MCH(j),ST(k),SC(l)))
       auxi=auxi-2d0*RVtdVts*IVtdVts      !3rd family contribution
     .          *CCD(j)*CCD(i)*RST(k,1)*RST(l,1)*CCT(i,k)*CCT(j,l)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),ST(l))
     . -2d0*RVcdVcs*IVcdVcs               !2nd family contribution
     .          *CCD(j)*CCD(i)*RSC(k,1)*RSC(l,1)*CCC(i,k)*CCC(j,l)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),SC(k),SC(l))
     . -(IVtdVts*RVcdVcs+RVtdVts*IVcdVcs) !3rd family/2nd interference
     .          *CCD(j)*CCD(i)*MCH(i)*MCH(j)
     .             *(RSC(k,1)*RST(l,1)*CCC(i,k)*CCT(j,l)
     .                      *D0B(MCH(i),MCH(j),SC(k),ST(l))
     .              +RST(k,1)*RSC(l,1)*CCT(i,k)*CCC(j,l)
     .                      *D0B(MCH(i),MCH(j),ST(k),SC(l)))
      enddo
       aux=aux+RVtdVts*RVudVus            !3rd family/1st interference
     .  *2d0*dsqrt(g2q)*CCD(i)*V(j,1)*CCD(j)*RST(k,1)*CCT(i,k)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),MUL)
     .        +RVcdVcs*RVudVus            !2nd family/1st interference
     .  *2d0*dsqrt(g2q)*CCD(i)*V(j,1)*CCD(j)*RSC(k,1)*CCC(i,k)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),SC(k),MUL)
       auxi=auxi-IVtdVts*RVudVus          !3rd family/1st interference
     .  *2d0*dsqrt(g2q)*CCD(i)*V(j,1)*CCD(j)*RST(k,1)*CCT(i,k)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),MUL)
     .        -IVcdVcs*RVudVus            !2nd family/1st interference
     .  *2d0*dsqrt(g2q)*CCD(i)*V(j,1)*CCD(j)*RSC(k,1)*CCC(i,k)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),SC(k),MUL)
      enddo                               !1st family contribution
       aux=aux+RVudVus**2*g2q*V(j,1)*V(i,1)*CCD(j)*CCD(i)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),MUL,MUL)
      enddo
      enddo
      C1SLLCHAR=-aux/4d0/(GF*MW)**2
      C2SLLCHAR=aux/16d0/(GF*MW)**2
      C1SLLCHARI=-auxi/4d0/(GF*MW)**2
      C2SLLCHARI=auxi/16d0/(GF*MW)**2

      aux=0d0
      auxi=0d0
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+(RVtdVts**2-IVtdVts**2)    !3rd family contribution
     .          *CCD(j)*CCD(i)*RST(l,1)**2*CCT(j,k)*CCT(i,k)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),ST(l))
     . +(RVcdVcs**2-IVcdVcs**2)           !2nd family contribution
     .          *CCD(j)*CCD(i)*RSC(l,1)**2*CCC(j,k)*CCC(i,k)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),SC(k),SC(l))
     . +(RVtdVts*RVcdVcs-IVtdVts*IVcdVcs) !3rd family/2nd interference
     .          *CCD(j)*CCD(i)*MCH(i)*MCH(j)
     . *(RST(l,1)**2*CCC(j,k)*CCC(i,k)*D0B(MCH(i),MCH(j),SC(k),ST(l))
     .  +RSC(l,1)**2*CCT(j,k)*CCT(i,k)*D0B(MCH(i),MCH(j),ST(k),SC(l)))
       auxi=auxi-2d0*RVtdVts*IVtdVts      !3rd family contribution
     .          *CCD(j)*CCD(i)*RST(l,1)**2*CCT(j,k)*CCT(i,k)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),ST(l))
     . -2d0*RVcdVcs*IVcdVcs               !2nd family contribution
     .          *CCD(j)*CCD(i)*RSC(l,1)**2*CCC(j,k)*CCC(i,k)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),SC(k),SC(l))
     . -(IVtdVts*RVcdVcs+RVtdVts*IVcdVcs) !3rd family/2nd interference
     .          *CCD(j)*CCD(i)*MCH(i)*MCH(j)
     . *(RST(l,1)**2*CCC(j,k)*CCC(i,k)*D0B(MCH(i),MCH(j),SC(k),ST(l))
     .  +RSC(l,1)**2*CCT(j,k)*CCT(i,k)*D0B(MCH(i),MCH(j),ST(k),SC(l)))
      enddo
       aux=aux+RVtdVts*RVudVus            !3rd family/1st interference
     .  *CCD(i)*CCD(j)*(CCT(i,k)*CCT(j,k)+g2q*V(j,1)*V(i,1)*RST(K,1)**2)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),MUL)
     .        +RVcdVcs*RVudVus            !2nd family/1st interference
     .  *CCD(i)*CCD(j)*(CCC(i,k)*CCC(j,k)+g2q*V(j,1)*V(i,1)*RSC(K,1)**2)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),SC(k),MUL)
       auxi=auxi-IVtdVts*RVudVus          !3rd family/1st interference
     .  *CCD(i)*CCD(j)*(CCT(i,k)*CCT(j,k)+g2q*V(j,1)*V(i,1)*RST(K,1)**2)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),MUL)
     .        -IVcdVcs*RVudVus            !2nd family/1st interference
     .  *CCD(i)*CCD(j)*(CCC(i,k)*CCC(j,k)+g2q*V(j,1)*V(i,1)*RSC(K,1)**2)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),SC(k),MUL)
      enddo                               !1st family contribution
       aux=aux+RVudVus**2*g2q*V(j,1)*V(i,1)*CCD(j)*CCD(i)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),MUL,MUL)
      enddo
      enddo
      C1LRCHAR=-aux/2d0/(GF*MW)**2*MD0/MS0
      C1LRCHARI=-auxi/2d0/(GF*MW)**2*MD0/MS0

      aux=0d0
      auxi=0d0
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+(RVtdVts**2-IVtdVts**2)    !3rd family contribution
     .          *CCD(j)**2*RST(k,1)*RST(l,1)*CCT(i,k)*CCT(i,l)
     .                   *D2B(MCH(i),MCH(j),ST(k),ST(l))
     . +(RVcdVcs**2-IVcdVcs**2)           !2nd family contribution
     .          *CCD(j)**2*RSC(k,1)*RSC(l,1)*CCC(i,k)*CCC(i,l)
     .                   *D2B(MCH(i),MCH(j),SC(k),SC(l))
     . +(RVtdVts*RVcdVcs-IVtdVts*IVcdVcs) !3rd family/2nd interference
     .          *CCD(j)**2*(RST(k,1)*RSC(l,1)*CCT(i,k)*CCC(i,l)
     .                   *D2B(MCH(i),MCH(j),ST(k),SC(l))
     .                     +RSC(k,1)*RST(l,1)*CCC(i,k)*CCT(i,l)
     .                   *D2B(MCH(i),MCH(j),SC(k),ST(l)))
       auxi=auxi-2d0*RVtdVts*IVtdVts      !3rd family contribution
     .          *CCD(j)**2*RST(k,1)*RST(l,1)*CCT(i,k)*CCT(i,l)
     .                   *D2B(MCH(i),MCH(j),ST(k),ST(l))
     . -2d0*RVcdVcs*IVcdVcs               !2nd family contribution
     .          *CCD(j)**2*RSC(k,1)*RSC(l,1)*CCC(i,k)*CCC(i,l)
     .                   *D2B(MCH(i),MCH(j),SC(k),SC(l))
     . -(IVtdVts*RVcdVcs+RVtdVts*IVcdVcs) !3rd family/2nd interference
     .          *CCD(j)**2*(RST(k,1)*RSC(l,1)*CCT(i,k)*CCC(i,l)
     .                   *D2B(MCH(i),MCH(j),ST(k),SC(l))
     .                     +RSC(k,1)*RST(l,1)*CCC(i,k)*CCT(i,l)
     .                   *D2B(MCH(i),MCH(j),SC(k),ST(l)))
      enddo
       aux=aux+RVtdVts*RVudVus            !3rd family/1st interference
     .  *2d0*CCD(j)**2*dsqrt(g2q)*V(i,1)
     .            *RST(k,1)*CCT(i,k)*D2B(MCH(i),MCH(j),ST(k),MUL)
     .        +RVcdVcs*RVudVus            !2nd family/1st interference
     .  *2d0*CCD(j)**2*dsqrt(g2q)*V(i,1)
     .            *RSC(k,1)*CCC(i,k)*D2B(MCH(i),MCH(j),SC(k),MUL)
       auxi=auxi-IVtdVts*RVudVus          !3rd family/1st interference
     .  *2d0*CCD(j)**2*dsqrt(g2q)*V(i,1)
     .            *RST(k,1)*CCT(i,k)*D2B(MCH(i),MCH(j),ST(k),MUL)
     .        -IVcdVcs*RVudVus            !2nd family/1st interference
     .  *2d0*CCD(j)**2*dsqrt(g2q)*V(i,1)
     .            *RSC(k,1)*CCC(i,k)*D2B(MCH(i),MCH(j),SC(k),MUL)
      enddo                               !1st family contribution
       aux=aux+RVudVus**2
     .          *g2q*V(i,1)**2*CCD(j)**2*D2B(MCH(i),MCH(j),MUL,MUL)
      enddo
      enddo
      C2LRCHAR=-aux/2d0/(GF*MW)**2*MD0/MS0
      C2LRCHARI=-auxi/2d0/(GF*MW)**2*MD0/MS0

*  Running to sc0 [5] App.C
      etaS=asf(scR)/asc0

      CVLLCHAR=etaS**(6d0/21d0)*(1d0+asc0/4d0/Pi*1.3707d0*(1d0-etaS))
     .                 *CVLLCHAR

      aux=C1SLLCHAR
      C1SLLCHAR=(1.0153d0*etaS**(-0.6916d0)-0.0153d0*etaS**(0.7869d0)
     . +asc0/4d0/Pi*(etaS**(-0.6916d0)*(5.6478d0-6.0350d0*etaS)
     .              +etaS**(0.7869d0)*(0.3272d0+0.0600*etaS)))*aux
     . +(1.9325d0*(etaS**(-0.6916d0)-etaS**(0.7869d0))
     . +asc0/4d0/Pi*(etaS**(-0.6916d0)*(10.7494d0-37.9209d0*etaS)
     .   +etaS**(0.7869d0)*(41.2556d0-14.0841d0*etaS)))*C2SLLCHAR

      C2SLLCHAR=(0.0081d0*(etaS**(0.7869d0)-etaS**(-0.6916d0))
     . +asc0/4d0/Pi*(etaS**(0.7869d0)*(-0.0618d0-0.0315d0*etaS)
     .              +etaS**(-0.6916d0)*(0.0454d0+0.0479d0*etaS)))*aux
     . +(1.0153d0*etaS**(0.7869d0)-0.0153d0*etaS**(-0.6916)
     . +asc0/4d0/Pi*(etaS**(-0.6916d0)*(0.0865d0+0.3007d0*etaS)
     .    +etaS**(0.7869d0)*(-7.7870d0+7.3999d0*etaS)))*C2SLLCHAR

      aux=C1LRCHAR
      C1LRCHAR=(etaS**(3d0/21d0)
     . +asc0/4d0/Pi*(0.9219d0*etaS**(-24d0/21d0)+etaS**(3d0/21d0)
     .   *(-2.2194d0+1.2975d0*etaS)))*aux
     . +asc0/4d0/Pi*1.3828d0*(etaS**(24d0/21d0)-etaS**(-24d0/21d0))
     .   *C2LRCHAR

      C2LRCHAR=(2d0/3d0*(etaS**(3d0/21d0)-etaS**(-24d0/21d0))
     . +asc0/4d0/Pi*(etaS**(-24d0/21d0)*(-6.4603d0+15.7415d0*etaS)
     .   +etaS**(3d0/21d0)*(-10.1463d0+0.8650d0*etaS)))*aux
     . +(etaS**(-24d0/21d0)+asc0/4d0/Pi*(0.9219d0*etaS**(24d0/21d0)
     .   +etaS**(-24d0/21d0)*(9.6904d0-10.6122d0*etaS)))*C2LRCHAR

      CVLLCHARI=etaS**(6d0/21d0)*(1d0+asc0/4d0/Pi*1.3707d0*(1d0-etaS))
     .                 *CVLLCHARI

      aux=C1SLLCHARI
      C1SLLCHARI=(1.0153d0*etaS**(-0.6916d0)-0.0153d0*etaS**(0.7869d0)
     . +asc0/4d0/Pi*(etaS**(-0.6916d0)*(5.6478d0-6.0350d0*etaS)
     .              +etaS**(0.7869d0)*(0.3272d0+0.0600*etaS)))*aux
     . +(1.9325d0*(etaS**(-0.6916d0)-etaS**(0.7869d0))
     . +asc0/4d0/Pi*(etaS**(-0.6916d0)*(10.7494d0-37.9209d0*etaS)
     .   +etaS**(0.7869d0)*(41.2556d0-14.0841d0*etaS)))*C2SLLCHARI

      C2SLLCHARI=(0.0081d0*(etaS**(0.7869d0)-etaS**(-0.6916d0))
     . +asc0/4d0/Pi*(etaS**(0.7869d0)*(-0.0618d0-0.0315d0*etaS)
     .              +etaS**(-0.6916d0)*(0.0454d0+0.0479d0*etaS)))*aux
     . +(1.0153d0*etaS**(0.7869d0)-0.0153d0*etaS**(-0.6916)
     . +asc0/4d0/Pi*(etaS**(-0.6916d0)*(0.0865d0+0.3007d0*etaS)
     .    +etaS**(0.7869d0)*(-7.7870d0+7.3999d0*etaS)))*C2SLLCHARI

      aux=C1LRCHARI
      C1LRCHARI=(etaS**(3d0/21d0)
     . +asc0/4d0/Pi*(0.9219d0*etaS**(-24d0/21d0)+etaS**(3d0/21d0)
     .   *(-2.2194d0+1.2975d0*etaS)))*aux
     . +asc0/4d0/Pi*1.3828d0*(etaS**(24d0/21d0)-etaS**(-24d0/21d0))
     .   *C2LRCHARI

      C2LRCHARI=(2d0/3d0*(etaS**(3d0/21d0)-etaS**(-24d0/21d0))
     . +asc0/4d0/Pi*(etaS**(-24d0/21d0)*(-6.4603d0+15.7415d0*etaS)
     .   +etaS**(3d0/21d0)*(-10.1463d0+0.8650d0*etaS)))*aux
     . +(etaS**(-24d0/21d0)+asc0/4d0/Pi*(0.9219d0*etaS**(24d0/21d0)
     .   +etaS**(-24d0/21d0)*(9.6904d0-10.6122d0*etaS)))*C2LRCHARI
     
*          - Double Penguin contributions [4] Eqs.(6.12)-(6.22):
      aux=0d0                      ! sd-effective Higgs coupling
      do i=1,2
      do k=1,2
      aux=aux+CCT(i,k)**2*BB1(0d0,MCH(i)**2,ST(k)**2,QSTSB)
     .  -2d0*MCH(i)/MS0*CCD(i)*RST(k,1)*CCT(i,k)
     .                     *BB0(0d0,MCH(i)**2,ST(k)**2,QSTSB)
!     .  +(MD0/MS0)**2*CCD(i)**2*RST(k,1)**2
!     .                     *BB1(0d0,MCH(i)**2,ST(k)**2,QSTSB)
      enddo
      aux=aux-(g2q*V(i,1)**2*BB1(0d0,MCH(i)**2,MUL**2,QSTSB)
     .    -2d0*dsqrt(g2q)*V(i,1)*MCH(i)/MS0*CCD(i)
     .                      *BB0(0d0,MCH(i)**2,MUL**2,QSTSB))
!     .    +(MC/H1Q)**2*V(i,2)**2
!     .                      *BB1(0d0,MCH(i)**2,MUR**2,QSTSB)
!     .    +(MS0/H2Q)**2*U(i,2)**2
!     .                      *BB1(0d0,MCH(i)**2,MUL**2,QSTSB)
      enddo
      aux=aux/(32d0*pi**2)
      sigRLsd=MS0*aux/(H1Q*(1d0+epst2-aux)*(1d0+epst2))

      auxi=0d0                      ! ds-effective Higgs coupling
      do i=1,2
      do k=1,2
      auxi=auxi-2d0*MCH(i)/MS0*CCD(i)*RST(k,1)*CCT(i,k)
     .                        *BB0(0d0,MCH(i)**2,ST(k)**2,QSTSB)
     .            +CCT(i,k)**2*BB1(0d0,MCH(i)**2,ST(k)**2,QSTSB)
     .  +CCD(i)**2*RST(k,1)**2*BB1(0d0,MCH(i)**2,ST(k)**2,QSTSB)
      enddo
      auxi=auxi-(g2q*V(i,1)**2*BB1(0d0,MCH(i)**2,MUL**2,QSTSB)
!     .   +(MC0/H1Q)**2*V(i,2)**2*BB1(0d0,MCH(i)**2,MUR**2,QSTSB)
     .   -2d0*dsqrt(g2q)*V(i,1)*CCD(i)*MCH(i)
     .                        *BB0(0d0,MCH(i)**2,MUL**2,QSTSB)
     .              +CCD(i)**2*BB1(0d0,MCH(i)**2,MUL**2,QSTSB))
      enddo
      auxi=auxi/(32d0*pi**2)
      sigLRsd=MD0*auxi/(H1Q*(1d0+epst1-aux)*(1d0+epst2))
      
      aux=0d0
      do i=1,3
      aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)**2
     . *sgn(SMASS(i)**2-MK0**2)/
     . dsqrt((SMASS(i)**2-MK0**2)**2+(SMASS(i)*WIDTH(i))**2)
      enddo
      do i=1,2
      aux=aux+(PCOMP(i,1)+PCOMP(i,1)*tanb**2)**2
     . /(1d0+tanb**2)*sgn(PMASS(i)**2-MK0**2)/
     . dsqrt((PMASS(i)**2-MK0**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo

      C2LRDPH=-(4d0*pi/(GF*MW))**2*sigRLsd*sigLRsd*aux
     .                   *(RVtdVts**2-IVtdVts**2)
      C2LRDPHI=-(4d0*pi/(GF*MW))**2*sigRLsd*sigLRsd*aux
     .                   *(-2d0*RVtdVts*IVtdVts)

      aux=0d0
      do i=1,3
      aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)**2
     . *sgn(SMASS(i)**2-MK0**2)/
     . dsqrt((SMASS(i)**2-MK0**2)**2+(SMASS(i)*WIDTH(i))**2)
      enddo
      do i=1,2
      aux=aux-(PCOMP(i,1)+PCOMP(i,1)*tanb**2)**2
     . /(1d0+tanb**2)*sgn(PMASS(i)**2-MK0**2)/
     . dsqrt((PMASS(i)**2-MK0**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo

      C1SLLDPH=-(4d0*pi/(GF*MW))**2*sigRLsd**2*aux/2d0
     .                   *(RVtdVts**2-IVtdVts**2)
      C1SLLDPHI=-(4d0*pi/(GF*MW))**2*sigRLsd**2*aux/2d0
     .                   *(-2d0*RVtdVts*IVtdVts)

*          - Summary

      I=1                              ! 0: SM; 1: NMSSM
      IF(I.eq.1)then
      CVLL=CVLLHIG+CVLLCHAR

      C1SLL=C1SLLHIG+C1SLLCHAR+C1SLLDPH
      C2SLL=C2SLLHIG+C2SLLCHAR

      C1LR=C1LRHIG+C1LRCHAR
      C2LR=C2LRHIG+C2LRCHAR+C2LRDPH

      CVLLI=CVLLHIGI+CVLLCHARI

      C1SLLI=C1SLLHIGI+C1SLLCHARI+C1SLLDPHI
      C2SLLI=C2SLLHIGI+C2SLLCHARI

      C1LRI=C1LRHIGI+C1LRCHARI
      C2LRI=C2LRHIGI+C2LRCHARI+C2LRDPHI
      ENDIF

*	 2) Results for DMK
      eta5=asc0/asf(4.6d0)
      eta4=asf(4.6d0)/ascL

*          - `Bag' parameters at scale scL=3 GeV from lattice [4'] Table XIII
* Given the spread of results among lattice collaborations, we choose a somewhat central value
* and larger error bands accordingly.
      BVLL=0.52d0
      B1SLL=0.50d0
      B2SLL=0.35d0
      B1LR=0.65d0
      B2LR=0.90d0

      aux=18.75d0 !(MK0/(runmass(0.005d0,scL)+runmass(0.095d0,scL)))**2
      B1SLL=aux*B1SLL
      B2SLL=aux*B2SLL
      B1LR=aux*B1LR
      B2LR=aux*B2LR

*          - Running between scL and sc0 [5] Eqs.(3.20)-(3.38) and (7.28-7.32)

      PVLL=eta5**(6d0/23d0)*eta4**(6d0/25d0)*(1d0+ascL/4d0/Pi
     .             *(1.7917d0-0.1644d0*eta4-1.6273*eta4*eta5))*BVLL

      P1SLL=-5d0/8d0*(1.0153d0*eta5**(-0.6315d0)*eta4**(-0.5810d0)
     .               -0.0153d0*eta5**(0.7184d0)*eta4**(0.6610d0)
     . +ascL/4d0/pi*(eta5**(-0.6315d0)*(0.0020d0*eta4**(1.6610d0)
     .   +(4.2458d0+0.5700d0*eta4-5.2272d0*eta4*eta5)
     .                                         *eta4**(-0.5810d0))
     .   +eta5**(0.7184d0)*eta4**(0.6610d0)
     .              *(0.3640d0+0.0064d0*eta4+0.0724d0*eta4*eta5)))
     .                                                        *B1SLL
     . -3d0/2d0*(0.0081d0*(eta5**(0.7184d0)*eta4**(0.6610d0)
     .                    -eta5**(-0.6315d0)*eta4**(-0.5810d0))
     . +ascL/4d0/pi*(eta5**(-0.6315d0)*(-0.0011d0*eta4**(1.6610d0)
     .   +(0.0587d0-0.0045d0*eta4+0.0415d0*eta4*eta5)
     .                                          *eta4**(-0.5810d0))
     .              -eta5**(0.7184d0)*(-0.0003d0*eta4**(0.4190d0)
     .   +(0.0534d0+0.0034d0*eta4+0.0380d0*eta4*eta5)
     .                                    *eta4**(0.6610d0))))*B2SLL

      P2SLL=-5d0/8d0*(1.9325d0*(eta5**(-0.6315d0)*eta4**(-0.5810d0)
     .                         -eta5**(0.7184d0)*eta4**(0.6610d0))
     . +ascL/4d0/pi*(eta5**(-0.6315d0)*(0.0038d0*eta4**(1.6610d0)
     .    +(8.0810d0+1.0848d0*eta4-38.8778d0*eta4*eta5)
     .                                         *eta4**(-0.5810d0))
     .            +eta5**(0.7184d0)*(-4.2075d0*eta4**(0.4190d0)
     . +(45.9008d0+0.8087d0*eta4-12.7939d0*eta4*eta5)
     .                                    *eta4**(0.6610d0))))*B1SLL
     . -3d0/2d0*(1.0153d0*eta5**(0.7184d0)*eta4**(0.6610d0)
     .          -0.0153d0*eta5**(-0.6315d0)*eta4**(-0.5810d0)
     . +ascL/4d0/pi*(eta5**(-0.6315d0)*(-0.0020d0*eta4**(1.6610d0)
     .     +(0.1117d0-0.0086d0*eta4+0.3083d0*eta4*eta5)
     .                                       *eta4**(-0.5810d0))
     .   +eta5**(0.7184d0)*(0.0334d0*eta4**(0.4190d0)
     .     +(-6.7398d0-0.4249d0*eta4+6.7219d0*eta4*eta5)
     .                                    *eta4**(0.6610d0))))*B2SLL

      P1LR=-(eta5**(3d0/23d0)*eta4**(3d0/25d0)
     . +ascL/4d0/Pi*(eta5**(-24d0/23d0)*
     .     (0.9279d0*eta4**(-24d0/25d0)-0.0029d0*eta4**(28d0/25d0))
     .   +eta5**(3d0/23d0)*eta4**(3d0/25d0)
     .      *(-2.0241d0-0.0753d0*eta4+1.1744d0*eta4*eta5)))*B1LR/2d0
     . +3d0/4d0*(2d0/3d0*(eta5**(3d0/23d0)*eta4**(3d0/25d0)
     .                   -eta5**(-24d0/23d0)*eta4**(-24d0/25d0))
     . +ascL/4d0/Pi*(eta5**(3d0/23d0)*(5d0*eta4**(1d0/25d0)
     .   +(-16.6828d0-0.0502d0*eta4+0.7829d0*eta4*eta5)
     .                                        *eta4**(3d0/25d0))
     .   +eta5**(-24d0/23d0)*(-0.0019d0*eta4**(28d0/25d0)
     .    +(-4.4701d0-0.8327d0*eta4+16.2548d0*eta4*eta5)
     .                                   *eta4**(-24d0/25d0))))*B2LR

      P2LR=-ascL/4d0/Pi*(1.3875d0*eta5**(26d0/23d0)
     .                                          *eta4**(28d0/25d0)
     . +eta5**(-24d0/23d0)*(-1.3918d0*eta4**(-24d0/25d0)+0.0043d0
     .                                 *eta4**(28d0/25d0)))*B1LR/2d0
     . +3d0/4d0*(eta5**(-24d0/23d0)*eta4**(-24d0/25d0)
     . +ascL/4d0/Pi*(eta5**(-24d0/23d0)*(0.0029d0*eta4**(28d0/25d0)
     .   +(6.7052d0+1.2491d0*eta4-8.8822d0*eta4*eta5)
     .                                        *eta4**(-24d0/25d0))
     .          +0.9250d0*eta5**(26d0/23d0)*eta4**(28d0/25d0)))*B2LR
     

*          - New physics contribution to DMK [5] Eqs.(7.24)-(7.27)
      aux=PVLL*CVLL+P1SLL*C1SLL+P2SLL*C2SLL+P1LR*C1LR+P2LR*C2LR
      
      DMK=GF**2*MW**2/(24d0*pi**2)*MK0*FK0**2*aux

*          - Error estimate
*      First, error bars from uncertainties on lattice Bag parameters: 
*      (2sigma, added quadratically)

      auxi=eta5**(6d0/23d0)*eta4**(6d0/25d0)*(1d0+ascL/4d0/Pi
     .             *(1.7917d0-0.1644d0*eta4-1.6273*eta4*eta5))*CVLL

      DMKMax=4d0*(auxi*0.02d0)**2

      auxi=-5d0/8d0*(1.0153d0*eta5**(-0.6315d0)*eta4**(-0.5810d0)
     .               -0.0153d0*eta5**(0.7184d0)*eta4**(0.6610d0)
     . +ascL/4d0/pi*(eta5**(-0.6315d0)*(0.0020d0*eta4**(1.6610d0)
     .   +(4.2458d0+0.5700d0*eta4-5.2272d0*eta4*eta5)
     .                                         *eta4**(-0.5810d0))
     .   +eta5**(0.7184d0)*eta4**(0.6610d0)
     .              *(0.3640d0+0.0064d0*eta4+0.0724d0*eta4*eta5)))
     .                                                        *C1SLL
     . -5d0/8d0*(1.9325d0*(eta5**(-0.6315d0)*eta4**(-0.5810d0)
     .                         -eta5**(0.7184d0)*eta4**(0.6610d0))
     . +ascL/4d0/pi*(eta5**(-0.6315d0)*(0.0038d0*eta4**(1.6610d0)
     .    +(8.0810d0+1.0848d0*eta4-38.8778d0*eta4*eta5)
     .                                         *eta4**(-0.5810d0))
     .            +eta5**(0.7184d0)*(-4.2075d0*eta4**(0.4190d0)
     . +(45.9008d0+0.8087d0*eta4-12.7939d0*eta4*eta5)
     .                                    *eta4**(0.6610d0))))*C2SLL
     
      DMKMax=DMKMax+4d0*(auxi*0.03d0*18.75d0)**2

      auxi=-3d0/2d0*(0.0081d0*(eta5**(0.7184d0)*eta4**(0.6610d0)
     .                    -eta5**(-0.6315d0)*eta4**(-0.5810d0))
     . +ascL/4d0/pi*(eta5**(-0.6315d0)*(-0.0011d0*eta4**(1.6610d0)
     .   +(0.0587d0-0.0045d0*eta4+0.0415d0*eta4*eta5)
     .                                          *eta4**(-0.5810d0))
     .              -eta5**(0.7184d0)*(-0.0003d0*eta4**(0.4190d0)
     .   +(0.0534d0+0.0034d0*eta4+0.0380d0*eta4*eta5)
     .                                    *eta4**(0.6610d0))))*C1SLL
     . -3d0/2d0*(1.0153d0*eta5**(0.7184d0)*eta4**(0.6610d0)
     .          -0.0153d0*eta5**(-0.6315d0)*eta4**(-0.5810d0)
     . +ascL/4d0/pi*(eta5**(-0.6315d0)*(-0.0020d0*eta4**(1.6610d0)
     .     +(0.1117d0-0.0086d0*eta4+0.3083d0*eta4*eta5)
     .                                       *eta4**(-0.5810d0))
     .   +eta5**(0.7184d0)*(0.0334d0*eta4**(0.4190d0)
     .     +(-6.7398d0-0.4249d0*eta4+6.7219d0*eta4*eta5)
     .                                    *eta4**(0.6610d0))))*C2SLL

      DMKMax=DMKMax+4d0*(auxi*0.06d0*18.75d0)**2

      auxi=-(eta5**(3d0/23d0)*eta4**(3d0/25d0)
     . +ascL/4d0/Pi*(eta5**(-24d0/23d0)*
     .     (0.9279d0*eta4**(-24d0/25d0)-0.0029d0*eta4**(28d0/25d0))
     .   +eta5**(3d0/23d0)*eta4**(3d0/25d0)
     .      *(-2.0241d0-0.0753d0*eta4+1.1744d0*eta4*eta5)))/2d0*C1LR
     . -ascL/4d0/Pi*(1.3875d0*eta5**(26d0/23d0)
     .                                          *eta4**(28d0/25d0)
     . +eta5**(-24d0/23d0)*(-1.3918d0*eta4**(-24d0/25d0)+0.0043d0
     .                                 *eta4**(28d0/25d0)))/2d0*C2LR

      DMKMax=DMKMax+4d0*(auxi*0.07d0*18.75d0)**2
      
      auxi=3d0/4d0*(2d0/3d0*(eta5**(3d0/23d0)*eta4**(3d0/25d0)
     .                   -eta5**(-24d0/23d0)*eta4**(-24d0/25d0))
     . +ascL/4d0/Pi*(eta5**(3d0/23d0)*(5d0*eta4**(1d0/25d0)
     .   +(-16.6828d0-0.0502d0*eta4+0.7829d0*eta4*eta5)
     .                                        *eta4**(3d0/25d0))
     .   +eta5**(-24d0/23d0)*(-0.0019d0*eta4**(28d0/25d0)
     .    +(-4.4701d0-0.8327d0*eta4+16.2548d0*eta4*eta5)
     .                                   *eta4**(-24d0/25d0))))*C1LR
     . +3d0/4d0*(eta5**(-24d0/23d0)*eta4**(-24d0/25d0)
     . +ascL/4d0/Pi*(eta5**(-24d0/23d0)*(0.0029d0*eta4**(28d0/25d0)
     .   +(6.7052d0+1.2491d0*eta4-8.8822d0*eta4*eta5)
     .                                        *eta4**(-24d0/25d0))
     .          +0.9250d0*eta5**(26d0/23d0)*eta4**(28d0/25d0)))*C2LR

      DMKMax=DMKMax+4d0*(auxi*0.06d0*18.75d0)**2

      DMKMax=(GF**2*MW**2/(24d0*pi**2)*MK0*FK0**2)**2*DMKMax
      DMKmin=DMKMax

*      Second: Uncertainty from matching: 30% on each BSM contribution

      aux=0d0

      IF(I.eq.1)then                                             ! 30% BSM
       aux=aux+0.3d0*(dabs(PVLL*CVLLHIG+P1SLL*C1SLLHIG+P2SLL*C2SLLHIG
     . +P1LR*C1LRHIG+P2LR*C2LRHIG)+dabs(PVLL*CVLLCHAR+P1SLL*C1SLLCHAR
     . +P2SLL*C2SLLCHAR+P1LR*C1LRCHAR+P2LR*C2LRCHAR)
     . +dabs(P1SLL*C1SLLDPH+P2LR*C2LRDPH))
      ENDIF

      DMKMax=DMK+dsqrt(DMKMax)+GF**2*MW**2/(24d0*pi**2)*MK0*FK0**2
     .                   *aux
      DMKmin=DMK-dsqrt(DMKmin)-GF**2*MW**2/(24d0*pi**2)*MK0*FK0**2
     .                   *aux

*      Finally: CKM and lattice form factor
      aux=4d0*dsqrt((0.6d-3/8.4d-3)**2+(2.7d-3/40d-3)**2
     .                                      +(0.0009d0/FK0)**2)

*      Total error bars for the BSM contributions
      DMKmax=DMKMax*(1d0+aux)
      DMKmin=DMKmin*(1d0-aux)

*          - SM contribution  [6]
      etacc=1.87d0
      etact=0.496d0
      etatt=0.5765d0
      BK=ascL**(-6d0/25d0)*(1d0+ascL/4d0/Pi*(6719d0/3750d0))*BVLL
      xt=(MT0/MW)**2
      xc=(1.279d0/MW)**2
      
      DMKSM=(GF*MW)**2/6d0/Pi**2*MK0*FK0**2*BK
     .    *((RVcdVcs**2-IVcdVcs**2)*etacc*S0(xc)
     .     +2d0*(RVcdVcs*RVtdVts-IVcdVcs*IVtdVts)*etact*S02(xc,xt)
     .     +(RVtdVts**2-IVtdVts**2)*etatt*S0(xt))

*          - Long-Distance contribution: (20 +/- 10)% 1304.6835
      DMKSM=DMKSM+0.2d0*3.48d-15

      DMK=DMK+DMKSM

*   SM uncertainty
      aux=(FK0**2*BK)**2*(((RVcdVcs**2-IVcdVcs**2)*0.76d0*S0(xc))**2              ! eta factors
     .    +(2d0*(RVcdVcs*RVtdVts-IVcdVcs*IVtdVts)*0.047d0*S02(xc,xt))**2
     .    +((RVtdVts**2-IVtdVts**2)*0.0065d0*S0(xt))**2
     .    +(5d-4*etacc*S0(xc))**2+(9d-6*2d0*etact*S02(xc,xt))**2                      ! CKM
     .    +(3d-8*etatt*S0(xt))**2)
     .    +((2d0*0.0009*FK0*BK)**2                                                ! Form factors
     .      +(FK0**2*ascL**(-6d0/25d0)*(1d0+ascL/4d0/Pi*(6719d0/3750d0))
     .         *0.03d0)**2)
     .    *((RVcdVcs**2-IVcdVcs**2)*etacc*S0(xc)
     .     +2d0*(RVcdVcs*RVtdVts-IVcdVcs*IVtdVts)*etact*S02(xc,xt)
     .     +(RVtdVts**2-IVtdVts**2)*etatt*S0(xt))**2
     .    +(0.1d0*3.48d-15*6d0*Pi**2/(GF*MW)**2/MK0)**2                           ! Long-distance

       DMKMax=DMKMax+DMKSM+(GF*MW)**2/6d0/Pi**2*MK0*2d0*dsqrt(aux)
       DMKmin=DMKmin+DMKSM-(GF*MW)**2/6d0/Pi**2*MK0*2d0*dsqrt(aux)

       IF(dabs(DMKmin).gt.dabs(DMKMax))then
        DMKMax=dabs(DMKmin)
        DMKmin=Max(0d0,DMKmin)
       ENDIF
      
*	- Conversion to ps^-1 and comparison with experiment
      DMK=DMK/(6.58211915d-13)
      DMKmin=DMKmin/(6.58211915d-13)
      DMKMax=DMKMax/(6.58211915d-13)
      
      prob(60)=0d0

      IF(DMKmin.GE.DMKexpMax)
     .     PROB(60)=DMKmin/DMKexpMax-1d0
      IF(DMKmax.LE.DMKexpmin)
     .     PROB(60)=DMKmax/DMKexpMin-1d0
      
!      print*,DMKmin,DMK,DMKmax

*	 2) Results for epsK
*          - New physics contribution to epsK [5] Eqs.(7.26)-(7.27)
      aux=PVLL*CVLLI+P1SLL*C1SLLI+P2SLL*C2SLLI+P1LR*C1LRI+P2LR*C2LRI
      
      epsK=GF**2*MW**2/(48d0*dsqrt(2d0)*pi**2)*MK0*FK0**2*aux/3.48d-15
     .           *dsqrt(2d0)*DDSIN(43.52d0/180d0*PI)    ! Corrections from angle =/= Pi/4

*          - Error estimate
*      First, error bars from uncertainties on lattice Bag parameters: 
*      (2sigma, added quadratically)

      auxi=eta5**(6d0/23d0)*eta4**(6d0/25d0)*(1d0+ascL/4d0/Pi
     .             *(1.7917d0-0.1644d0*eta4-1.6273*eta4*eta5))*CVLLI

      epsKMax=4d0*(auxi*0.02d0)**2

      auxi=-5d0/8d0*(1.0153d0*eta5**(-0.6315d0)*eta4**(-0.5810d0)
     .               -0.0153d0*eta5**(0.7184d0)*eta4**(0.6610d0)
     . +ascL/4d0/pi*(eta5**(-0.6315d0)*(0.0020d0*eta4**(1.6610d0)
     .   +(4.2458d0+0.5700d0*eta4-5.2272d0*eta4*eta5)
     .                                         *eta4**(-0.5810d0))
     .   +eta5**(0.7184d0)*eta4**(0.6610d0)
     .              *(0.3640d0+0.0064d0*eta4+0.0724d0*eta4*eta5)))
     .                                                        *C1SLLI
     . -5d0/8d0*(1.9325d0*(eta5**(-0.6315d0)*eta4**(-0.5810d0)
     .                         -eta5**(0.7184d0)*eta4**(0.6610d0))
     . +ascL/4d0/pi*(eta5**(-0.6315d0)*(0.0038d0*eta4**(1.6610d0)
     .    +(8.0810d0+1.0848d0*eta4-38.8778d0*eta4*eta5)
     .                                         *eta4**(-0.5810d0))
     .            +eta5**(0.7184d0)*(-4.2075d0*eta4**(0.4190d0)
     . +(45.9008d0+0.8087d0*eta4-12.7939d0*eta4*eta5)
     .                                    *eta4**(0.6610d0))))*C2SLLI
     
      epsKMax=epsKMax+4d0*(auxi*0.03d0*18.75d0)**2

      auxi=-3d0/2d0*(0.0081d0*(eta5**(0.7184d0)*eta4**(0.6610d0)
     .                    -eta5**(-0.6315d0)*eta4**(-0.5810d0))
     . +ascL/4d0/pi*(eta5**(-0.6315d0)*(-0.0011d0*eta4**(1.6610d0)
     .   +(0.0587d0-0.0045d0*eta4+0.0415d0*eta4*eta5)
     .                                          *eta4**(-0.5810d0))
     .              -eta5**(0.7184d0)*(-0.0003d0*eta4**(0.4190d0)
     .   +(0.0534d0+0.0034d0*eta4+0.0380d0*eta4*eta5)
     .                                    *eta4**(0.6610d0))))*C1SLLI
     . -3d0/2d0*(1.0153d0*eta5**(0.7184d0)*eta4**(0.6610d0)
     .          -0.0153d0*eta5**(-0.6315d0)*eta4**(-0.5810d0)
     . +ascL/4d0/pi*(eta5**(-0.6315d0)*(-0.0020d0*eta4**(1.6610d0)
     .     +(0.1117d0-0.0086d0*eta4+0.3083d0*eta4*eta5)
     .                                       *eta4**(-0.5810d0))
     .   +eta5**(0.7184d0)*(0.0334d0*eta4**(0.4190d0)
     .     +(-6.7398d0-0.4249d0*eta4+6.7219d0*eta4*eta5)
     .                                    *eta4**(0.6610d0))))*C2SLLI

      epsKMax=epsKMax+4d0*(auxi*0.06d0*18.75d0)**2

      auxi=-(eta5**(3d0/23d0)*eta4**(3d0/25d0)
     . +ascL/4d0/Pi*(eta5**(-24d0/23d0)*
     .     (0.9279d0*eta4**(-24d0/25d0)-0.0029d0*eta4**(28d0/25d0))
     .   +eta5**(3d0/23d0)*eta4**(3d0/25d0)
     .      *(-2.0241d0-0.0753d0*eta4+1.1744d0*eta4*eta5)))/2d0*C1LRI
     . -ascL/4d0/Pi*(1.3875d0*eta5**(26d0/23d0)
     .                                          *eta4**(28d0/25d0)
     . +eta5**(-24d0/23d0)*(-1.3918d0*eta4**(-24d0/25d0)+0.0043d0
     .                                 *eta4**(28d0/25d0)))/2d0*C2LRI

      epsKMax=epsKMax+4d0*(auxi*0.07d0*18.75d0)**2
      
      auxi=3d0/4d0*(2d0/3d0*(eta5**(3d0/23d0)*eta4**(3d0/25d0)
     .                   -eta5**(-24d0/23d0)*eta4**(-24d0/25d0))
     . +ascL/4d0/Pi*(eta5**(3d0/23d0)*(5d0*eta4**(1d0/25d0)
     .   +(-16.6828d0-0.0502d0*eta4+0.7829d0*eta4*eta5)
     .                                        *eta4**(3d0/25d0))
     .   +eta5**(-24d0/23d0)*(-0.0019d0*eta4**(28d0/25d0)
     .    +(-4.4701d0-0.8327d0*eta4+16.2548d0*eta4*eta5)
     .                                   *eta4**(-24d0/25d0))))*C1LRI
     . +3d0/4d0*(eta5**(-24d0/23d0)*eta4**(-24d0/25d0)
     . +ascL/4d0/Pi*(eta5**(-24d0/23d0)*(0.0029d0*eta4**(28d0/25d0)
     .   +(6.7052d0+1.2491d0*eta4-8.8822d0*eta4*eta5)
     .                                        *eta4**(-24d0/25d0))
     .          +0.9250d0*eta5**(26d0/23d0)*eta4**(28d0/25d0)))*C2LRI

      epsKMax=epsKMax+4d0*(auxi*0.06d0*18.75d0)**2

      epsKMax=(GF**2*MW**2/(48d0*dsqrt(2d0)*pi**2)*MK0*FK0**2
     .           *2d0*DDSIN(43.52d0/180d0*PI)/3.48d-15)**2*epsKMax
      epsKmin=epsKMax

*      Second: Uncertainty from matching: 30% on each BSM contribution

      aux=0d0

      IF(I.eq.1)then                                             ! 30% BSM
       aux=aux+0.3d0*(dabs(PVLL*CVLLHIGI+P1SLL*C1SLLHIGI
     .                +P2SLL*C2SLLHIGI+P1LR*C1LRHIGI+P2LR*C2LRHIGI)
     .   +dabs(PVLL*CVLLCHARI+P1SLL*C1SLLCHARI+P2SLL*C2SLLCHARI
     .                              +P1LR*C1LRCHARI+P2LR*C2LRCHARI)
     .   +dabs(P1SLL*C1SLLDPHI+P2LR*C2LRDPHI))
      ENDIF

      epsKMax=epsK+dsqrt(epsKMax)+GF**2*MW**2/(48d0*dsqrt(2d0)*pi**2)
     .      *MK0*FK0**2*aux/3.48d-15*dsqrt(2d0)*DDSIN(43.52d0/180d0*PI)    ! Corrections from angle =/= Pi/4
      epsKmin=epsK-dsqrt(epsKmin)-GF**2*MW**2/(48d0*dsqrt(2d0)*pi**2)
     .      *MK0*FK0**2*aux/3.48d-15*dsqrt(2d0)*DDSIN(43.52d0/180d0*PI)    ! Corrections from angle =/= Pi/4

*          - SM contribution  [6]
      
      epsKSM=(GF*MW)**2/12d0/dsqrt(2d0)/Pi**2*MK0*FK0**2*BK
     .    *(-2d0*RVcdVcs*IVcdVcs*etacc*S0(xc)
     .     -2d0*(IVcdVcs*RVtdVts+RVcdVcs*IVtdVts)*etact*S02(xc,xt)
     .     -2d0*RVtdVts*IVtdVts*etatt*S0(xt))/3.48d-15
     .                 *dsqrt(2d0)*DDSIN(43.52d0/180d0*PI)    ! Corrections from angle =/= Pi/4  
     
!      epsKSM=0.923d0*3.4646d4*0.737d0*0.0406d0**2*0.22537d0**2*0.353d0
!     .*(1d0-0.22537d0**4/8d0)*(0.0406d0**2*(1d0-0.117d0)*etatt*S0(xt)
!     .                  +etact*S02(xc,xt)-etacc*S0(xc))

*          - Long-distance effects [7]

      epsKSM=epsKSM+DDSIN(43.52d0/180d0*PI)*(-1.63d-4)

      epsK=epsK+epsKSM

*   SM uncertainty
       epsKMax=epsKMax+epsKSM+2d0*dsqrt((0.28d-3)**2     ! SM 
     .                 +(0.28d-4)**2+(1.6d-2*epsKSM)**2) ! + long-dist.
       epsKmin=epsKmin+epsKSM-2d0*dsqrt((0.28d-3)**2     ! SM 
     .                 +(0.28d-4)**2+(1.6d-2*epsKSM)**2) ! + long-dist.
      epsKmin=Max(0d0,epsKmin)

*    correction factor
      epsK=0.923d0*epsK
      epsKmin=0.911d0*epsKmin
      epsKMax=0.935d0*epsKMax

!      print*,epsKmin,epsKSM,epsKmax

      IF(epsKmin.GE.epsKexpMax)
     .     PROB(60)=PROB(60)+epsKmin/epsKexpMax-1d0
      IF(epsKmax.LE.epsKexpmin)
     .     PROB(60)=PROB(60)+epsKmax/epsKexpMin-1d0

      END

*********************************************************************

      double precision function S02(x,z)
      
      implicit none
      double precision x,z
      if(dabs(z-1d0).gt.1d-3)then
       S02=x*(dlog(z/x)-3d0*z/4d0/(1d0-z)
     .     -3d0*z**2*dlog(z)/4d0/(1d0-z)**2)
      else
       S02=-x*(dlog(x)+3d0/8d0)
      endif
      return
      end
      
      
      
      
