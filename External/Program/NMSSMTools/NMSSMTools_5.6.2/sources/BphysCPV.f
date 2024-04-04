      SUBROUTINE BOTTOMONIUM_CPV(PROB)

***********************************************************************
*   Subroutine to check bottomonium physics driven by a light CP-odd Higgs:
*      PROB(38) =/= 0  excluded by Upsilon(1S) -> H/A gamma
*      PROB(39) =/= 0  excluded by eta_b(1S) mass measurement
***********************************************************************

      IMPLICIT NONE

      INTEGER I

      DOUBLE PRECISION PROB(*),PI
      DOUBLE PRECISION MY,MBQM,BRYMUMU,C,ZZ,AP,GYEE,Rmax,MA
      DOUBLE PRECISION ALPHAS,ALSMY,DELTA,UPSILONTAU,UPSILONMU
      DOUBLE PRECISION M0,RETA,GAM2,XX,YY,F,YMAX,FMAX,D,MEMAX,UU,VV

      DOUBLE PRECISION ALEM0
      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION WIDTH(5),HCWIDTH
      DOUBLE PRECISION CU(5),CUP(5),CD(5),CDP(5),CB(5),CBP(5),CJ(5),
     . CJP(5),CI(5),CG(5),CGP(5),CV(5),CZG(5),CZGP(5),CL(5),CLP(5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),
     . BRCC(5),BRBB(5),BRTT(5),BRWW(5),BRZZ(5),BRGG(5),BRZG(5)
      DOUBLE PRECISION DDCOS

      COMMON/ALEM0/ALEM0
      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/HIWIDTH/WIDTH,HCWIDTH
      COMMON/HNSMBR/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG
      COMMON/HNSMCOUP/CU,CUP,CD,CDP,CB,CBP,CJ,CJP,CI,CG,CGP,CV,CZG,CZGP,
     . CL,CLP

      PI=4d0*DATAN(1d0)

      PROB(38)=0d0
      PROB(39)=0d0

      DO I=1,5

      MA=dsqrt(MH0(I))

* Test on Upsilon(1S) -> H/A gamma (from CLEO)

      MY=9.46d0 ! Upsilon(1S) mass
      MBQM=4.9d0 ! b quark mass in quark models
      ALSMY=ALPHAS(MY,2) ! alpha_s at MY, 2 loop calculation
      BRYMUMU=2.48d-2 ! BR(Upsilon(1S) -> mu mu)

      IF(MA.LT.MY)THEN

       ZZ=1d0-MA**2/MY**2 ! energy fraction of the photon
       AP=6d0*ZZ**.2d0 ! Nason function for QCD corrections
       C=1d0+4d0*ALSMY/(3d0*PI)*(4d0-AP) ! QCD corrections
       DELTA=1.2d0**2/MBQM**2 ! function for rel. corrections
       C=C* ! relativistic corrections (for MA<~8.8 GeV)
     .  (MY**2-MA**2)**2/(4d0*MBQM**2-MA**2)**2*(1d0-
     .  DELTA/3d0*(36d0*MBQM**2+MA**2)/(4d0*MBQM**2-MA**2))

       C=MAX(C,1d-6)

       RMAX=0d0
       RMAX=dsqrt(2d0)*PI*ALEM0*UPSILONTAU(MA)/(GF*MBQM**2*ZZ*C*BRYMUMU)
       IF(RMAX.NE.0d0)PROB(38)=DDIM(CBP(I)**2*BRLL(I)/RMAX,1d0)

       RMAX=0d0
       RMAX=dsqrt(2d0)*PI*ALEM0*UPSILONMU(MA)/(GF*MBQM**2*ZZ*C*BRYMUMU)
       IF(RMAX.NE.0d0)PROB(38)=PROB(38)+DDIM(CBP(I)**2*BRMM(I)/RMAX,1d0)

       IF(WIDTH(I).GT.1d-2)PROB(38)=-PROB(38)

      ENDIF

* Test on etab(1S) mass difference (BABAR - theory)

      M0=9.389d0 ! etab(1S) mass
      GYEE=1.34d-6 ! Gamma(Upsilon(1S) -> e+ e-)
      RETA=GYEE*9d0*MY**2/(4d0*ALEM0**2)*
     .  (1d0+16d0*ALSMY/(3d0*PI)) ! radial wave fun. at the origin
       ! Resolution of the 3rd degree eq. for the limit on Xd
      GAM2=(M0*2d-2)**2
      XX=MA**2-M0**2

      IF(MA.LT.M0)THEN
       MEMAX=M0-3d-2
      ELSE
       MEMAX=M0+4d-2
      ENDIF
      YMAX=MEMAX**2-M0**2
      FMAX=XX*YMAX*(1d0+GAM2/(XX+YMAX)**2)

      D=XX**2-GAM2/27d0
      IF(D.LT.0d0)THEN
       UU=2d0*DSQRT(GAM2)/DSQRT(3d0)
       VV=-3d0*DSQRT(3d0)*XX/DSQRT(GAM2)
       YY=UU*DDCOS(DACOS(VV)/3d0+4d0*PI/3d0)
       F=XX*YY*(1d0+GAM2/(XX+YY)**2)
       FMAX=MAX(FMAX,F)
      ENDIF

      RMAX=8d0*PI*MZ**2/(3d0*(g1+g2)/2d0*RETA*M0**3)*FMAX
      PROB(39)=PROB(39)+DDIM(CBP(I)**2/RMAX,1d0)

      ENDDO

      END


      SUBROUTINE BSG_CPV(PAR,PROB)

*      Subroutine for B Physics observables

* Literature:
*
*   - THEORETICAL FORMULAE:
* [1] A. J. Buras, P. H. Chankowski, J. Rosiek, L. Slawianowska,
*     ' Delta M(d, s), B0(d, s) ---> mu+ mu- and B ---> X(s) gamma
*      in supersymmetry at large tan beta.'
*     Nucl.Phys.B659:3,2003, e-Print: hep-ph/0210145
*
* [2] K.Chetyrkin, M.Misiak, M.Munz,
*     'Weak Radiative B-Meson Decay Beyond Leading Logarithms'
*     Phys.Lett.B400:206-219,1997,
*     Erratum-ibid.B425:414,1998, e-Print: hep-ph/9612313
*
* [3] M.Ciuchini, G.Degrassi, P.Gambino, G.Giudice
*     'Next to Leading QCD Corrections to B -> Xs gamma: Standard Model
*     and Two-Higgs Doublet Model' Nucl.Phys.B527: (1998), 21
*     e-Print: hep-ph/9710335
*
* [4] M.Ciuchini, G.Degrassi, P.Gambino, G.Giudice,
*     `Next-to-leading QCD corrections to B ---> X(s) gamma in supersymmetry,'
*     Nucl.\ Phys.\ B {\bf 534} (1998) 3, e-Print: [hep-ph/9806308].
*
* [5] C.Bobeth, M.Misiak, J.Urban,
*     `Matching conditions for b ---> s gamma and b ---> s gluon in extensions 
*     of the standard model,'  Nucl. Phys. B bf 567 (2000) 153
*     e-Print: [hep-ph/9904413]
*
* [6] C.Bobeth, A.Buras, T.Ewerth,
*     `Anti-B ---> X(s) l+ l- in the MSSM at NNLO,'
*      Nucl. Phys. B 713 (2005) 522, e-Print: [hep-ph/0409293].
*
* [7] P.Gambino and M.Misiak,
*     `Quark mass effects in anti-B ---> X(s gamma),'
*      Nucl. Phys. B 611 (2001) 338, e-Print: [hep-ph/0104034].
*
* [8] M.Czakon, U.Haisch and M.Misiak,
*     ``Four-Loop Anomalous Dimensions for Radiative Flavour-Changing Decays,''
*     JHEP {\bf 0703} (2007) 008, e-Print: [hep-ph/0612329].
*
* [9] P.Gambino, U.Haisch,
*     'Complete Electroweak Matching for Radiative B Decays'
*     JHEP 0110:020,2001, e-Print: hep-ph/0109058
*
* [10] T.Hurth, E.Lunghi, W.Porod,
*     'Untagged B -> Xs+d gamma CP asymmetry as a probe for New Physics'
*     Nucl.Phys.B704:56-74,2005, e-Print: hep-ph/0312260
*
* [11] M.Czakon, P.Fiedler, T.Huber, M.Misiak, T.Schutzmeier and M.Steinhauser,
*     `The (Q$_{7}$ , Q$_{1,}_{2}$) contribution to $ \overline{B}\to {X}_s\gamma $ 
*     at $ \mathcal{O}\left({\alpha}_{\mathrm{s}}^2\right) $,'
*     JHEP 1504 (2015) 168, e-Print: [arXiv:1503.01791 [hep-ph]].
*
* [12] M.Misiak and M.Steinhauser,
*     `NNLO QCD corrections to the anti-B ---> X(s) gamma matrix elements 
*     using interpolation in m(c),' Nucl. Phys. B 764 (2007) 62
*     e-Print: [hep-ph/0609241].
*
* [13] A.J.Buras, A.Czarnecki, M.Misiak and J.Urban,
*     `Completing the NLO QCD calculation of anti-B ---> X(s gamma),'
*     Nucl. Phys. B 631 (2002) 219, e-Print: [hep-ph/0203135].
*
* [14] M.Misiak et al.,
*     `Updated NNLO QCD predictions for the weak radiative B-meson decays,'
*     Phys. Rev. Lett.  114 (2015) 22,  221801, e-Print: [arXiv:1503.01789 [hep-ph]].
*
* [15] C.W.Bauer,
*      `Corrections to moments of the photon spectrum in the inclusive decay 
*      B ---> X(s) gamma,' Phys. Rev. D 57 (1998) 5611
*      [Phys.\ Rev.\ D {\bf 60} (1999) 099907], e-Print: [hep-ph/9710513].
*
* [16] T.Ewerth, P.Gambino and S.Nandi,
*      `Power suppressed effects in anti-B ---> X(s) gamma at O(alpha(s)),'
*      Nucl. Phys. B 830 (2010) 278, e-Print: [arXiv:0911.2175 [hep-ph]].
*
* [17] M.Benzke, S.J.Lee, M.Neubert and G.Paz,
*      `Factorization at Subleading Power and Irreducible Uncertainties in 
*      $\bar B\to X_s\gamma$ Decay,' JHEP 1008 (2010) 099
*      e-Print: [arXiv:1003.5012 [hep-ph]].
*
* [18] C.Bobeth, M.Gorbahn, T.Hermann, M.Misiak, E.Stamou and M.Steinhauser,
*      `B_{s,d} -> l+ l- in the Standard Model with Reduced Theoretical 
*      Uncertainty,' Phys. Rev. Lett. 112 (2014) 101801
*      e-Print: [arXiv:1311.0903 [hep-ph]].
*
* [19] C.Bobeth, A.J.Buras, F.Kruger and J.Urban,
*      `QCD corrections to $\bar{B} \to X_{d,s} \nu \bar{\nu}$, 
*      $\bar{B}_{d,s} \to \ell^{+} \ell^{-}$, $K \to \pi \nu \bar{\nu}$ 
*      and $K_{L} \to \mu^{+} \mu^{-}$ in the MSSM,'
*      Nucl. Phys. B 630 (2002) 87, e-Print: [hep-ph/0112305].
*
* [20] A.G. Akeroyd, S. Recksiegel,
*     ' The Effect of H+- on B+- ---> tau+- nu(tau) and
*       B+- ---> mu+- muon neutrino'
*     J.Phys.G29:2311-2317,2003, e-Print: hep-ph/0306037
*
* [21] T.Huber, T.Hurth and E.Lunghi,
*     `Inclusive B --> X_s l+l-: complete angular analysis and a thorough study 
*     of collinear photons,' JHEP {\bf 1506} (2015) 176
*     e-Print: [arXiv:1503.04849 [hep-ph]].
*
* [22] A.S.Cornell and N.Gaur,
*     `The Forward backward asymmetries of $B \to X_{s} \tau^{+} \tau^{-}$ in 
*     the MSSM,' JHEP {\bf 0309} (2003) 030
*     e-Print: [hep-ph/0308132].
*
* [23] A.Crivellin, C.Greub and A.Kokulu,
*     `Explaining $B\to D\tau\nu$, $B\to D^*\tau\nu$ and $B\to \tau\nu$ in a 
*     2HDM of type III,'  Phys.\ Rev.\ D {\bf 86} (2012) 054014
*     e-Print: [arXiv:1206.2634 [hep-ph]].
*
* [24] A.J.Buras, S.Jager and J.Urban,
*     `Master formulae for Delta F=2 NLO QCD factors in the standard model and 
*     beyond,' Nucl.Phys.B 605 (2001) 600
*     e-Print: [hep-ph/0102316].
*
* [25] D.Becirevic, V.Gimenez, G.Martinelli, M.Papinutto and J.Reyes,
*     `B parameters of the complete set of matrix elements of delta B = 2 
*     operators from the lattice,' JHEP 0204 (2002) 025
*     e-Print: [hep-lat/0110091].
*
* [26] H.M.Asatrian and C.Greub,
*     `Tree-level contribution to $\bar{B} to X_d \gamma$ using fragmentation 
*     functions,' Phys. Rev. D 88 (2013) 7,  074014
*     e-Print: [arXiv:1305.6464 [hep-ph]].
*
* [27] W.Altmannshofer and D.M.Straub,
*     `New physics in $B \to K^*\mu\mu$?,', Eur.\ Phys.\ J.\ C {\bf 73} (2013) 2646
*     e-Print: [arXiv:1308.1501 [hep-ph]].
*
* [28] A.J.Buras, J.Girrbach-Noe, C.Niehoff and D.M.Straub,
*     `$ B\to {K}^{\left(\ast \right)}\nu \overline{\nu} $ decays in the Standard Model and beyond,'
*     JHEP {\bf 1502} (2015) 184, e-Print: [arXiv:1409.4557 [hep-ph]].
*
*   - SOURCES FOR EXPERIMENTAL/LATTICE QCD DATA:
* [1'] Y.Amhis et al. [Heavy Flavor Averaging Group (HFAG) Collaboration],
*      `Averages of $b$-hadron, $c$-hadron, and $\tau$-lepton properties as of 
*      summer 2014,' e-Print : arXiv:1412.7515 [hep-ex].
*
* [2'] Heavy Flavor Averaging Group,
*     www.slac.stanford.edu/xorg/hfag/osc/PDG_2014/#DMD
*
* [3'] Heavy Flavor Averaging Group,
*     www.slac.stanford.edu/xorg/hfag/osc/PDG_2014/#DMS
*
* [4'] Heavy Flavor Averaging Group,
*      www.slac.stanford.edu/xorg/hfag/rare/2013/radll/btosg.pdf
*
* [5'] S.Aoki et al.,
*      `Review of lattice results concerning low-energy particle physics,'
*      Eur. Phys. J. C 74 (2014) 2890,  e-Print: [arXiv:1310.8555 [hep-lat]]. 
*      updates at http://itpwiki.unibe.ch/flag/
*
* [6'] K.A. Olive et al. (Particle Data Group), Chin. Phys. C, 38, 090001 (2014).
*
* [7'] Heavy Flavor Averaging Group,
*     www.slac.stanford.edu/xorg/hfag/rare/2014/bs/index.html
*     http://www.slac.stanford.edu/xorg/hfag/rare/2014/radll/OUTPUT/HTML/radll_table4.html
*
* [8'] P. Ball and R. Fleischer,
*      Eur.\ Phys.\ J.\  C {\bf 48} (2006) 413 [arXiv:hep-ph/0604249].
*
* [9'] Heavy Flavor Averaging Group,
*      www.slac.stanford.edu/xorg/hfag/rare/2014/radll/OUTPUT/HTML/radll_table4.html
*
* [10'] M. Iwasaki et al. [Belle Collaboration], hep-ex/0503044
*       J. P. Lees et al. [BaBar Collaboration], Phys. Rev. Lett. 112 (2014) 
*       211802 [arXiv:1312.5364[hep-ex]].
*
* [11'] Heavy Flavor Averaging Group,
*      http://www.slac.stanford.edu/xorg/hfag/semi/eps15/eps15_dtaunu.html
*      Formerly: A.Crivellin, J.Heeck and P.Stoffer,
*      `A perturbed lepton-specific two-Higgs-doublet model facing experimental 
*      hints for physics beyond the Standard Model,'
*      e-Print: arXiv:1507.07567 [hep-ph].
*
* [12'] P.del Amo Sanchez et al. [BaBar Collaboration],
*      `Study of $B \to X\gamma$ Decays and Determination of $|V_{td}/V_{ts}|$,'
*      Phys. Rev. D 82 (2010) 051101
*      e-Print: [arXiv:1005.4087 [hep-ex]].
*
* [13'] A.Crivellin and L.Mercolli,
*      `$B -> X_d \gamma$ and constraints on new physics,'
*      Phys. Rev. D 84 (2011) 114005
*      e-Print: [arXiv:1106.5499 [hep-ph]].
*
* [14'] R.Barate et al. [ALEPH Collaboration],
*      `Measurements of BR(b --> tau- anti-nu(tau) X) and BR(b --> tau- anti-nu(tau) D*+- X) 
*      and upper limits on BR (B- ---> tau- anti-nu(tau)) and BR (b---> s nu anti-nu),'
*      Eur. Phys. J. C 19 (2001) 213, e-Print: [hep-ex/0010022].
*
* [15'] J.P.Lees et al. [BaBar Collaboration],
*      `Search for $B → K^{(*)}ν\overlineν$ and invisible quarkonium decays,'
*      Phys. Rev. D 87 (2013) 11,  112005, e-Print: [arXiv:1303.7465 [hep-ex]].
*
* [16'] O.Lutz et al. [Belle Collaboration],
*      `Search for $B \to h^{(*)} \nu \bar{\nu}$ with the full Belle 
*      $\Upsilon(4S)$ data sample,'
*      Phys. Rev. D 87 (2013) 11,  111103, e-Print: [arXiv:1303.3719 [hep-ex]].
*

      IMPLICIT NONE

      INTEGER I,J,K,L,M
      DOUBLE PRECISION PAR(*),PROB(*)
      DOUBLE PRECISION Pi,aux,auxe,asf,H2,BB0,BB1,BB0L,sgn,sp2,ADG
      DOUBLE PRECISION aux1,aux2,au,ad,AAT,AAB,AAS,MTL2,MTR2,MBL2
      DOUBLE PRECISION XT,YT,ZT,mu,Mga2,gg1,gg2,scR,scR2,runmass
      DOUBLE PRECISION asmt,asmh,asmsusy,asc0,asmc,asmb,runmb
      DOUBLE PRECISION scb,sc0,MT0,MTH,mbh,MQU,MD,MS0,MC0,MB0,Vbdg,VVtc
      DOUBLE PRECISION VVc,VVu,Vtb2,Vcs2,Vud2,VtdVtb2,VtsVtb2,Vub2,Vbsg
      DOUBLE PRECISION VVtu,VVtuim,VVtud,VVtudim,dVub,dVtdVtb2,dVtsVtb2
      DOUBLE PRECISION BBs,dBBs,BBd,dBBd,fB,dfB,fBs,dfBs,fBd,dfBd,tauB0
      DOUBLE PRECISION M_Bs,tau_Bs,dtau_Bs,M_Bd,tau_Bd,dtau_Bd,MBu,tauB
      DOUBLE PRECISION DMdexpmin,DMdexpMax,DMsexpmin,DMsexpMax,
     .      BRBMUMUexpMax,BRBMUMUexpMin,BRBTAUNUexpmin,BRBTAUNUexpMax,
     .      BRSGexpmin,BRSGexpMax,BRDGexpmin,BRDGexpMax,
     .      BRBdMUMUexpMax,BRBdMUMUexpMin,RD_taulexpmin,RD_taulexpmax,
     .      RDs_taulexpmin,RDs_taulexpmax
      DOUBLE PRECISION BRBSllminexp,BRBSllMaxexp,BRBShllminexp,
     . BRBShllMaxexp
      DOUBLE PRECISION BRBptoKpllminexp,BRBptoKpllMaxexp
      DOUBLE PRECISION BRBpKpnunuexpMax,BRBKsnunuexpMax,BRBXsnunuexpMax
      DOUBLE PRECISION MCHH(2),etaH,D0B,D2B,etaS,S0
      DOUBLE PRECISION CVLLSM,CVLLHIG,C1SLLHIG,C2SLLHIG,C1LRHIG,C2LRHIG
      DOUBLE PRECISION CVLLHIGI,C1SLLHIGI,C2SLLHIGI,C1LRHIGI,C2LRHIGI
      DOUBLE PRECISION CVLLCHAR,C1SLLCHAR,C2SLLCHAR,C1LRCHAR,C2LRCHAR
      DOUBLE PRECISION CVLLCHARI,C1SLLCHARI,C2SLLCHARI,C1LRCHARI,
     . C2LRCHARI
      DOUBLE PRECISION C2LRDPH,C1SLLDPH,C2LRDPHI,C1SLLDPHI,CVLL,C1SLL
      DOUBLE PRECISION C2SLL,C1LR,C2LR,CVLLI,C1SLLI,C2SLLI,C1LRI,C2LRI
      DOUBLE PRECISION PRLH_tb(2),PRLH_ts(2),PLRH_tb(2),PLRH_ts(2)
      DOUBLE PRECISION PRLHI_tb(2),PRLHI_ts(2),PLRHI_tb(2),PLRHI_ts(2)
      DOUBLE PRECISION BVLL,B1SLL,B2SLL,B1LR,B2LR,PVLL,P1SLL,P2SLL
      DOUBLE PRECISION P1LR,P2LR,rh,ALEMMB
      DOUBLE PRECISION CcR,CcL,CcRI,CcLI,RD_taulSM,RDs_taulSM
      DOUBLE PRECISION sigRLbs,sigRLbd,sigLRbs,sigLRbd
      DOUBLE PRECISION sigRLbsI,sigRLbdI,sigLRbsI,sigLRbdI
      DOUBLE PRECISION C70SM,C80SM,C11SM,C41SM,C71SM,C81SM
      DOUBLE PRECISION C41HIG,C70HIG,C80HIG,dC7HIG,dC8HIG,C7HIGI,C8HIGI
      DOUBLE PRECISION C7HIG,C8HIG,C71HIG,C81HIG
      DOUBLE PRECISION C7CHARS,C8CHARS,C7CHAR,C8CHAR,C41CHAR,C71CHAR
      DOUBLE PRECISION C81CHAR,C70S0,C80S0,eta0
      DOUBLE PRECISION C7CHARSI,C8CHARSI,C7CHARI,C8CHARI
      DOUBLE PRECISION C71CHARI,C81CHARI,C70S0I,C80S0I
      DOUBLE PRECISION C70BSM,C80BSM,C71BSM,C81BSM,C41BSM,DC7BSM,DC8BSM
      DOUBLE PRECISION C70BSMI,C80BSMI,C71BSMI,C81BSMI,DC7BSMI,DC8BSMI
      DOUBLE PRECISION DC7BSM_b,DC8BSM_b,DC7BSM_bI,DC8BSM_bI
      DOUBLE PRECISION C70,C80,C70I,C80I,C11,C41,C71,C81,C71I,C81I,eta
      DOUBLE PRECISION FF1,FF2,FF3,FG1,FG2,FG3,ffh,fgh,esm,eh
      DOUBLE PRECISION echi,esfe,H17,H18,H27,H28,Q11,Q21,Q31,Q41
      DOUBLE PRECISION gg7,gg8,Delt7,Delt8,CCD(2,2),CCT(2,2,2)
      DOUBLE PRECISION aa(8),bb(4),hh(8),h8(4),af,afim,bf,bfim
      DOUBLE PRECISION m0011(8),m0021(8),m0031(8),m0041(8),m0051(8)
      DOUBLE PRECISION m0061(8),m0071(8),m0081(8),m0034(8),m0044(8)
      DOUBLE PRECISION m0054(8),m0064(8),m0074(8),m0084(8),m1012(8)
      DOUBLE PRECISION m1022(8),m1032(8),m1042(8),m1052(8),m1062(8)
      DOUBLE PRECISION m1072(8),m1082(8),m1112(8),m1122(8),m1132(8)
      DOUBLE PRECISION m1142(8),m1152(8),m1162(8),m1172(8),m1182(8)
      DOUBLE PRECISION C10b,C20b,C30b,C40b,C50b,C60b,C70b,C80b,C70bI
      DOUBLE PRECISION C11b,C21b,C31b,C41b,C51b,C61b,C71b,C81b,C80bI
      DOUBLE PRECISION C71bI,C81bI,C7EMb,C7EMbI,EPSew,EPSewI
      DOUBLE PRECISION ff11,ff12,ff17,ff18,ff22,ff27,ff28,ff47,ff48
      DOUBLE PRECISION ff77,ff78,ff88,ff17i,ff27i,ff18i,ff28i
      DOUBLE PRECISION MB_kin,z,MC_scb,MBp,dBPERT,T1,T2,T3,LO4B,NLO4B
      DOUBLE PRECISION K17,K17I,K27,K27I,K37,K37I,K47,K47I,K57,K57I
      DOUBLE PRECISION K67,K67I,K77,K78,K78I
      DOUBLE PRECISION BSGPERT,KC7BSM,KC8BSM,KC7BSMI,KC8BSMI
      DOUBLE PRECISION delt,delt2,lndelt,lndeltp,NP17,NP78,NP88,NP77
      DOUBLE PRECISION lambd1,lambd2,rho1,rho2,MCNP,HQET,CCSL,BRSL
      DOUBLE PRECISION CASM,CAH,CSH,CPH,CACHAR,CSCHAR,CPCHAR,CSHP,CPHP
      DOUBLE PRECISION CACHARI,CSCHARI,CPCHARI,C10eCHAR,C10eCHARI
      DOUBLE PRECISION fh20,fh21,fh20p,fh70,fh71,fh70p,fh30,fh31,fh30p
      DOUBLE PRECISION fc41,fc51,fc50,fc81,fc60,fc121,fc131,fc90,fc100
      DOUBLE PRECISION fc40,fc31,fc30,CSHPI,CPHPI
      DOUBLE PRECISION CA,CS,CP,CAI,CSI,CPI,DCA,DCS,DCP,DCAI,DCSI,DCPI

      DOUBLE PRECISION ALEM0
      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION Ytq,Ybq,MTOPQ,MBOTQ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION mt,mbnp,mtau,mmu,mel,MS,MC,MB,MPI,MSTRANGE
      DOUBLE PRECISION MSQ3,MSU3,MSD3,AT,AB
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION MST2P(2),UT(2,2,2),MSB2P(2),UB(2,2,2),MSL2(2),
     . UTAU(2,2,2),MSNT2
      DOUBLE PRECISION MSU2P(2),MSD2P(2),MSE2(2),MSNE2,MSMU2(2),
     . UMU(2,2,2)
      DOUBLE PRECISION MST2(2),MSB2(2),MSU2(2),MSD2(2)
      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION WIDTH(5),HCWIDTH
      DOUBLE PRECISION MCH2(2),U(2,2,2),V(2,2,2)
      DOUBLE PRECISION MNEU(5),N(5,5,2)
      DOUBLE PRECISION MGL
      DOUBLE PRECISION COCHSTbL(2,2,2),COCHSTbR(2,2,2),
     . COCHSBtL(2,2,2),COCHSBtR(2,2,2),COCHSLnL(2,2,2),
     . COCHSNlL(2,2),COCHSNlR(2,2),COCHSUdL(2,2,2),COCHSUdR(2,2,2),
     . COCHSDuL(2,2,2),COCHSDuR(2,2,2),COCHSEnL(2,2,2),
     . COCHSNeL(2,2),COCHSNeR(2,2)
      DOUBLE PRECISION CONESTtL(5,2,2),CONESTtR(5,2,2),
     . CONESBbL(5,2,2),CONESBbR(5,2,2),CONESLlL(5,2,2),CONESLlR(5,2,2),
     . CONESNnL(5,2),CONESUuL(5,2,2),CONESUuR(5,2,2),
     . CONESDdL(5,2,2),CONESDdR(5,2,2),CONESEeL(5,2,2),
     . CONESEeR(5,2,2)
      DOUBLE PRECISION eps0,epst0,epst1,epst2,epst3,epsts,epstb
      DOUBLE PRECISION epsY32,epsY31,epsY23,epsY13,epscs,epscb
      DOUBLE PRECISION eps0I,epst0I,epst1I,epst2I,epst3I,epstsI,epstbI
      DOUBLE PRECISION epsY32I,epsY31I,epsY23I,epsY13I,epscsI,epscbI
      DOUBLE PRECISION BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     .       DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     .       BRBtaunumax,BRBtaunumin
      DOUBLE PRECISION BRBSll,BRBSllmin,BRBSllMax,
     .             BRBShll,BRBShllmin,BRBShllMax,
     .             BRDG,BRDGmin,BRDGmax,
     .             BRBdMUMU,BRBdMUMUmin,BRBdMUMUmax,
     .             BRBXsnunu,BRBXsnunumin,BRBXsnunumax,
     .             BRBpKpnunu,BRBpKpnunumin,BRBpKpnunuMax,
     .             BRBKsnunu,BRBKsnunumin,BRBKsnunuMax,
     .             RD_taul,RD_taulmin,RD_taulmax,
     .             RDs_taul,RDs_taulmin,RDs_taulmax
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/ALEM0/ALEM0
      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/STSBSCALE/QSTSB
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/QQUARK/Ytq,Ybq,MTOPQ,MBOTQ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/SMFERM/mt,mbnp,mtau,mmu,mel,MS,MC,MB,MPI,MSTRANGE
      COMMON/RADCOR2/MSQ3,MSU3,MSD3,AT,AB
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/SFERM3SPEC/MST2P,UT,MSB2P,UB,MSL2,UTAU,MSNT2
      COMMON/SFERM1SPEC/MSU2P,MSD2P,MSE2,MSNE2,MSMU2,UMU
      COMMON/SFERMPSPEC/MST2,MSB2,MSU2,MSD2
      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/HIWIDTH/WIDTH,HCWIDTH
      COMMON/CHASPEC/MCH2,U,V
      COMMON/NEUSPEC/MNEU,N
      COMMON/GLUSPEC/MGL
      COMMON/CHSFfCOUP/COCHSTbL,COCHSTbR,COCHSBtL,COCHSBtR,COCHSLnL,
     . COCHSNlL,COCHSNlR,COCHSUdL,COCHSUdR,COCHSDuL,COCHSDuR,COCHSEnL,
     . COCHSNeL,COCHSNeR
      COMMON/NEUSFfCOUP/CONESTtL,CONESTtR,CONESBbL,CONESBbR,CONESLlL,
     . CONESLlR,CONESNnL,CONESUuL,CONESUuR,CONESDdL,CONESDdR,CONESEeL,
     . CONESEeR
      COMMON/EPSCOUP/eps0,epst0,epst1,epst2,epst3,epsts,epstb,
     .               epsY32,epsY31,epsY23,epsY13,epscs,epscb
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


***********************************************************************

*	I- PARAMETERS

      pi=4d0*datan(1d0)
!      csqb=0d0

*	 1) Alpha_s and m_top at various scales:
*   Charm Quark Mass:
       asmc=asf(mc) !0.366246d0

*   M_top:
       asmt=asf(MT) !0.108002d0

*   Charged Higgs Mass:
       asmh=asf(dsqrt(MHC))

*   M_top/bot at the Charged Higgs Mass Scale:
      mth=mt*(asmh/asmt)**(4d0/7d0)/(1d0+4d0/(3d0*pi)*asmt)
      mbh=runmb(dsqrt(MHC))

*   Susy scale (squark masses):
       asmsusy=asf(dsqrt(QSTSB))

*   m_top(m_top)(MSbar)
      MT0=MT/(1d0+4d0/(3d0*pi)*asmt+11d0/pi**2*asmt**2) !165d0

*	 2) EW gauge couplings at the SUSY scale:
      gg1=dsqrt(g1q)
      gg2=dsqrt(g2q)

*	 3) SUSY parameters:
*   Trig. Functions of Beta
      au=1d0/tanb
      ad=-tanb

*   Sfermions
      AAT=AT
      AAB=AB    
      AAS=PAR(13)

*   Charginos
      mu=PAR(4)
      Mga2=PAR(21)

*	 4) B-masses and lifetimes

      tauB=1.638d-12/(6.58211915d-25)            ! [1']
      tauB0=1.519d-12/(6.58211915d-25)            ! [1']
      mBu=5.27925d0                              ! [6']

      M_Bs=5.36677d0                       ! [6']
      M_Bd=5.27958d0                       ! [6']

c      tau_Bs=1.607d0/(6.58211915d-13)      ! [1']
c      dtau_Bs=0.020d0                      ! rel. uncertainty 2 sigma
      tau_Bs=1.509d0/(6.58211915d-13)
      dtau_Bs=0.008d0                      ! rel. uncertainty 2 sigma

      tau_Bd=1.520d0/(6.58211915d-13)      ! [1']
      dtau_Bd=0.008d0                      ! rel. uncertainty 2 sigma

*	 5) CKM Matrix:               [6']
      VVc=-0.9879d0                 ! (V_cs.V_cb^*)/(V_ts.V_tb^*)
      VVc=-1d0                      ! (V_cs.V_cb^*+V_us.Vub^*)/(V_ts.V_tb^*)
      VVu=-1d0                      ! (V_ud.V_ub^*+V_cd.V_cb^*)/(V_td.V_tb^*)
      Vtb2=(0.9991d0)**2               ! (V_tb)^2
      Vcs2=(0.97296d0)**2              ! (V_cs)^2
      Vud2=(0.97383d0)**2              ! (V_ud)^2 

*       Since V_ub(incl) differs considerably from V_ub(excl), we allow for 
*       the large range    3.15 10^-3 < V_ub < 5.11 10^-3  [6']
*      (V_ub)^2
      Vub2=(4.13d-3)**2
*      uncertainty on V_ub
      dVub=0.7d-3

*       From [8']: V_tb*V_td=(8.6 +/- 2.8) 10^-3  (2sigma)
*      (V_tb*V_td)^2
      VtdVtb2=(8.6d-3)**2
*      uncertainty on (V_tb*V_td)^2
      dVtdVtb2=2d0*dsqrt(VtdVtb2)*2.8d-3
*       From [8']: V_ts*V_tb=(41.3 +/- 1.4) 10^-3  (2sigma)
*      (Vtb.V_ts)^2
      VtsVtb2=(0.0413d0)**2
      dVtsVtb2=2d0*dsqrt(VtsVtb2)*1.4d-3    ! uncertainty on (Vtb.V_ts)^2

*      V_us.V_ub/V_ts.V_tb
      VVtu=-0.008d0
      VVtuim=0.0180d0
      VVtc=-0.9879d0

*      V_ud.V_ub/V_ts.V_tb
      VVtud=0.007d0
      VVtudim=-0.404d0

*	 6) Quark masses at the SUSY scale:
      scR=dsqrt(QSTSB)
      scR2=QSTSB

      MS0=runmass(0.095d0,scR)
      MC0=runmass(1.25d0,scR)
      MB0=runmb(scR)
      MQU=runmass(0.002d0,scR)
      MD=runmass(0.005d0,scR)


*	 7) Hadronic parameters (form factors):
*   [5'] provides several estimates for fBs from lattice computations:
*   fBs=(0.228+/-0.016) GeV (Nf=2), fBs=(0.2277+/-0.009) GeV (Nf=2+1),
*   fBs=(0.224+/-0.010) GeV (Nf=2+1+1). We shall use the compromise 
*   fBs=(0.226+/-0.012) GeV. (2sigma)
      fBs=0.226d0                         ! 
      dfBs=0.012d0/0.226d0                ! rel. uncertainty 2 sigma
*    BBs from [5']: BBs=1.33+/-0.06
      BBs=1.33d0
      dBBs=0.12d0

*   [5'] provides several estimates for fBs from lattice computations:
*   fBd=(0.189+/-0.016) GeV (Nf=2), fBd=(0.1905+/-0.0084) GeV (Nf=2+1),
*   fBd=(0.186+/-0.008) GeV (Nf=2+1+1). We shall use the compromise 
*   fBd=(0.1885+/-0.0105) GeV. (2sigma)
      fBd=0.1885d0                         ! 
      dfBd=0.0105d0/0.1885d0               ! rel. uncertainty 2 sigma
*    BBd from [5']: BBd=1.27+/-0.10
      BBd=1.270
      dBBd=0.20d0

*      [5'] provides several lattice results with Nf=2,2+1,2+1+1:
*      f_B=(0.189 +/- 0.016) GeV (Nf=2), f_B=(0.190.5 +/- 0.0084) GeV (Nf=2+1),
*      f_B=(0.186 +/- 0.008) GeV (Nf=2+1+1)
*        We decide to use an `averaged' value of f_B=(0.1885 +/- 0.0105) GeV
*        with an ad-hoc uncertainty to cover most of the allowed range.
      fB=0.1885d0             ! GeV       for B+ --> tau+ nu_tau
      dfB=0.0105d0            ! GeV       uncertainty

*       8) Experimental limits on B-processes
*      Delta Md from [1',2'], 2 sigma bounds:
*      0.5027ps-1 < DMd=0.5065ps-1 < 0.5103ps-1
      DMdexpmin=0.5027d0
      DMdexpMax=0.5103d0

*      Delta Ms from [1',3'], 2 sigma bounds:
*      17.715ps-1 < DMs=17.757ps-1 < 17.799ps-1
      DMsexpmin=17.715d0
      DMsexpMax=17.799d0

*      BR(B -> Xs gamma) from [1',4'], 2 sigma bounds:
*      3.02 10^-4 < BR(B -> Xs gamma)=3.32 10^-4 < 3.62 10^-4
      BRSGexpmin=3.02d-4
      BRSGexpMax=3.62d-4

*      BR(B -> Xd gamma) from [1',12',13'], 2 sigma bounds:
*      0.27 10^-5 < BR(B -> Xs gamma)=1.41 10^-5 < 2.55 10^-5
      BRDGexpmin=0.27d-5
      BRDGexpMax=2.55d-5

*      From [1',7']: 1.7 10^-9 < BR(Bs->mu+mu-) < 4.5 10^-9 (95% C.L.)
      BRBMUMUexpmax=4.5d-9
      BRBMUMUexpmin=1.7d-9

*      From [1',7']: 0.11 10^-9 < BR(Bd->mu+mu-) < 0.71 10^-9 (95% C.L.)
      BRBdMUMUexpmax=0.71d-9
      BRBdMUMUexpmin=0.11d-9

*      BR(B+ -> tau+ nu) from [1',9'], 2 sigma bounds:
*      0.78 10^-4 < BR(B+ -> tau+ nu)=1.06 10^-4 < 1.44 10^-4
      BRBTAUNUexpMax=1.44d-4
      BRBTAUNUexpmin=0.78d-4

*      BR(B -> Xsl+l-) from [10'], summarized in [21], 2 sigma bounds:
*   0.84 10^-6 < BR(B -> Xs l+l-)_low=1.58 10^-6 < 2.32 10^-6
*   0.28 10^-6 < BR(B -> Xs l+l-)_high=0.48 10^-6 < 0.68 10^-6
      BRBSllminexp=0.84d-6
      BRBSllMaxexp=2.32d-6
      BRBShllminexp=0.28d-6
      BRBShllMaxexp=0.68d-6

* RD from 1909.12524 Table 92 (UE 23.2.2022):
      RD_taulexpmin=0.340-2d0*0.030d0
      RD_taulexpmax=0.340+2d0*0.030d0

* RD* from 1909.12524 Table 92 (UE 23.2.2022):
      RDs_taulexpmin=0.295-2d0*0.014d0
      RDs_taulexpmax=0.295+2d0*0.014d0

*      BR[B+ -> K+ l+l-] from [27]
      BRBptoKpllminexp=1.04d-7-2d0*0.12d-7
      BRBptoKpllMaxexp=1.04d-7+2d0*0.12d-7

*      BR[B -> Xs nu nubar] from [14']
      BRBXsnunuexpMax=6.4d-4
      
*      BR[B+ -> K+ nu nubar] from [1',15']
      BRBpKpnunuexpMax=1.6d-5
      BRBpKpnunuexpMax=min(1.6d-5,4.9d-5*tauB/tauB0)    ! comparing with bound 
                                                        ! from neutral decay
*      BR[B0 -> K0* nu nubar] from [1',16']
      BRBKsnunuexpMax=5.5d-5
      BRBKsnunuexpMax=min(5.5d-5,4d-5*tauB0/tauB)       ! comparing with bound 
                                                        ! from charged decay


*	II- EFFECTIVE QUARK MASSES AND COUPLINGS TO HIGGS AT LARGE TANB
*       -> EPSILON COEFFICIENTS, notation following [1]:

*        1) Mass corrections
*   epsilontilde_J as in eq. (5.1) in [1] (up to a factor tanb)
*   from Delta m_d as in eq.(2.5), and Sigma as in App. (A.2)
*   epst1 = (epsilontilde as in [1])*tanb, epst2 (* tanb), epst3 (* tanb)

*	 * Gluino / squark contribution
      epst2=asf(scR)/(3d0*pi)*(BB1(0d0,MGL**2,MSD2(1),scR2)
     .                        +BB1(0d0,MGL**2,MSD2(2),scR2)
     . -2d0*MGL*(AAS*DDCOS(PhiAS-PhiM3)-mu*DDCOS(PhiM3+Phi01)*tanbq)
     .                           *BB0L(MGL**2,MSD2(1),MSD2(2),scR2))

      epst2I=asf(scR)/(3d0*pi)*(
     . -2d0*MGL*(AAS*DDSIN(PhiAS-PhiM3)+mu*DDSIN(PhiM3+Phi01)*tanbq)
     .                           *BB0L(MGL**2,MSD2(1),MSD2(2),scR2))

      epst3=0d0
      epst3I=0d0
      do k=1,2
      epst3=epst3+asf(scR)/(3d0*pi)*(BB1(0d0,MGL**2,MSB2(k),scR2)
     . -2d0*MGL/MB0*BB0(0d0,MGL**2,MSB2(k),scR2)
     .   *(DDCOS(PhiM3)*(UB(k,1,1)*UB(k,2,1)+UB(k,1,2)*UB(k,2,2))
     .    +DDSIN(PhiM3)*(UB(k,1,1)*UB(k,2,2)-UB(k,1,2)*UB(k,2,1))))
      epst3I=epst3I+asf(scR)/(3d0*pi)*(
     . -2d0*MGL/MB0*BB0(0d0,MGL**2,MSB2(k),scR2)
     .   *(DDCOS(PhiM3)*(UB(k,1,1)*UB(k,2,2)-UB(k,1,2)*UB(k,2,1))
     .    -DDSIN(PhiM3)*(UB(k,1,1)*UB(k,2,1)+UB(k,1,2)*UB(k,2,2))))
      enddo


*	 * Neutralino / squark contribution
      aux=0d0
      do i=1,5
       aux=aux-2d0/vdq*dsqrt(MNEU(i))                   
     .   *((CONESDdL(i,1,1)*N(i,4,1)+CONESDdL(i,1,2)*N(i,4,2))
     .                     *BB0(0d0,MNEU(i),MSD2(1),scR2)
     .    +(CONESDdR(i,2,1)*N(i,4,1)+CONESDdR(i,2,2)*N(i,4,2))
     .                     *BB0(0d0,MNEU(i),MSD2(2),scR2))
     .  +(CONESDdL(i,1,1)**2+CONESDdL(i,1,2)**2)
     .                     *BB1(0d0,MNEU(i),MSD2(1),scR2)
     .  +(CONESDdR(i,2,1)**2+CONESDdR(i,2,2)**2)
     .                     *BB1(0d0,MNEU(i),MSD2(2),scR2)
      enddo
      epst2=epst2+aux/(32d0*pi**2)

      aux=0d0
      do i=1,5
       aux=aux-2d0/vdq*dsqrt(MNEU(i))                   
     .   *((CONESDdL(i,1,2)*N(i,4,1)-CONESDdL(i,1,1)*N(i,4,2))
     .                     *BB0(0d0,MNEU(i),MSD2(1),scR2)
     .    +(CONESDdR(i,2,2)*N(i,4,1)-CONESDdR(i,2,1)*N(i,4,2))
     .                     *BB0(0d0,MNEU(i),MSD2(2),scR2))
      enddo
      epst2I=epst2I+aux/(32d0*pi**2)

      aux=0d0
      do i=1,5
      do k=1,2
      aux=aux+2d0*dsqrt(MNEU(i))/MB0*(CONESBbL(i,k,1)                    
     .         *CONESBbR(i,k,1)+CONESBbL(i,k,2)*CONESBbR(i,k,2))
     .                    *BB0(0d0,MNEU(i),MSB2(k),scR2)
     . +(CONESBbL(i,k,1)**2+CONESBbL(i,k,2)**2
     .     +CONESBbR(i,k,1)**2+CONESBbR(i,k,2)**2)
     .                    *BB1(0d0,MNEU(i),MSB2(k),scR2)
      enddo
      enddo
      epst3=epst3+1d0/(32d0*pi**2)*aux

      aux=0d0
      do i=1,5
      do k=1,2
      aux=aux+2d0*dsqrt(MNEU(i))/MB0*(CONESBbL(i,k,2)                    
     .         *CONESBbR(i,k,1)-CONESBbL(i,k,1)*CONESBbR(i,k,2))
     .                    *BB0(0d0,MNEU(i),MSB2(k),scR2)
      enddo
      enddo
      epst3I=epst3I+1d0/(32d0*pi**2)*aux

*	 * Chargino / squark contribution
      aux=0d0
      do i=1,2
      aux=aux+2d0*dsqrt(MCH2(i))/vdq
     .     *(COCHSUdL(i,1,1)*U(i,2,1)+COCHSUdL(i,1,2)*U(i,2,2))
     .          *BB0(0d0,MCH2(i),MSU2(1),scR2)
!     .  +(MC/vuq)**2*V(i,2)**2*BB1(0d0,MCH(i)**2,MUR**2,scR2)
     .  +(COCHSUdL(i,1,1)**2+COCHSUdL(i,1,2)**2) !+(MS0/vdq)**2*U(i,2)**2)
     .                        *BB1(0d0,MCH2(i),MSU2(1),scR2)
      enddo
      epst1=epst2+aux*Vud2/(32d0*pi**2)
      epst2=epst2+aux*Vcs2/(32d0*pi**2)

      aux=0d0
      do i=1,2
      aux=aux+2d0*dsqrt(MCH2(i))/vdq
     .     *(COCHSUdL(i,1,2)*U(i,2,1)-COCHSUdL(i,1,1)*U(i,2,2))
     .          *BB0(0d0,MCH2(i),MSU2(1),scR2)
      enddo
      epst1I=epst2I+aux*Vud2/(32d0*pi**2)
      epst2I=epst2I+aux*Vcs2/(32d0*pi**2)

      aux=0d0
      do i=1,2
      do k=1,2
      aux=aux+2d0*dsqrt(MCH2(i))/MB0*(COCHSTbR(i,k,1)*COCHSTbL(i,k,1)
     .                           +COCHSTbR(i,k,2)*COCHSTbL(i,k,2))
     .                            *BB0(0d0,MCH2(i),MST2(k),scR2)
     .  +(COCHSTbR(i,k,1)**2+COCHSTbR(i,k,2)**2+COCHSTbL(i,k,1)**2
     .        +COCHSTbL(i,k,2)**2)*BB1(0d0,MCH2(i),MST2(k),scR2)
      enddo
      enddo
      epst3=epst3+Vtb2/(32d0*pi**2)*aux

      aux=0d0
      do i=1,2
      do k=1,2
      aux=aux+2d0*dsqrt(MCH2(i))/MB0*(COCHSTbR(i,k,1)*COCHSTbL(i,k,2)
     .                           -COCHSTbR(i,k,2)*COCHSTbL(i,k,1))
     .                            *BB0(0d0,MCH2(i),MST2(k),scR2)
      enddo
      enddo
      epst3I=epst3I+Vtb2/(32d0*pi**2)*aux


*        2) Corrections to the neutral Higgs couplings
*     * epsilon_Y^(JI) as in eq. (5.1) in [1] (up to a factor yt^2
*   and a factor tanb), with (3.7) for lambda_0^(JI) and (3.53)
*   for the CKM matrix elements in terms of V^eff
      epsY32=0d0
      do i=1,2
      do k=1,2
      epsY32=epsY32
     .  +(COCHSTbL(i,k,1)**2+COCHSTbL(i,k,2)**2)
     .                     *BB1(0d0,MCH2(i),MST2(k),scR2)
     .  +2d0*dsqrt(MCH2(i))/MB0*(COCHSTbR(i,k,1)
     .     *COCHSTbL(i,k,1)+COCHSTbR(i,k,2)*COCHSTbL(i,k,2))
     .                     *BB0(0d0,MCH2(i),MST2(k),scR2)
!     .  +(MS0/vdq)**2*U(i,2)**2*RST(k,1)**2
!     .                     *BB1(0d0,MCH(i)**2,ST(k)**2,scR2)
      enddo
      epsY32=epsY32
     . +VVc*((COCHSUdL(i,1,1)**2+COCHSUdL(i,1,2)**2)
     .              *BB1(0d0,MCH2(i),MSU2(1),scR2)
     .    +2d0*dsqrt(MCH2(i))/vdq
     .   *(U(i,2,1)*COCHSUdL(i,1,1)+U(i,2,2)*COCHSUdL(i,1,2))
     .                      *BB0(0d0,MCH2(i),MSU2(1),scR2))
!     .    +(MC/vuq)**2*V(i,2)**2
!     .                      *BB1(0d0,MCH(i)**2,MUR**2,scR2)
!     .    +(MS0/vdq)**2*U(i,2)**2
!     .                      *BB1(0d0,MCH(i)**2,MUL**2,scR2)
      enddo
      epsY32=epsY32/(32d0*pi**2)

      epsY32I=0d0
      do i=1,2
      do k=1,2
      epsY32I=epsY32I
     .  +2d0*dsqrt(MCH2(i))/MB0*(COCHSTbR(i,k,1)
     .     *COCHSTbL(i,k,2)-COCHSTbR(i,k,2)*COCHSTbL(i,k,1))
     .                     *BB0(0d0,MCH2(i),MST2(k),scR2)
      enddo
      epsY32I=epsY32I
     .    +VVc*(2d0*dsqrt(MCH2(i))/vdq
     .   *(U(i,2,1)*COCHSUdL(i,1,2)-U(i,2,2)*COCHSUdL(i,1,1))
     .                      *BB0(0d0,MCH2(i),MSU2(1),scR2))
      enddo
      epsY32I=epsY32I/(32d0*pi**2)
      
      epsY31=0d0
      do i=1,2
      do k=1,2
      epsY31=epsY31
     .  +(COCHSTbL(i,k,1)**2+COCHSTbL(i,k,2)**2)
     .                     *BB1(0d0,MCH2(i),MST2(k),scR2)
     .  +2d0*dsqrt(MCH2(i))/MB0*(COCHSTbR(i,k,1)
     .          *COCHSTbL(i,k,1)+COCHSTbR(i,k,2)*COCHSTbL(i,k,2))
     .                     *BB0(0d0,MCH2(i),MST2(k),scR2)
!     .      +(MD/vdq)**2*U(i,2)**2*RST(k,1)**2
!     .                *BB1(0d0,MCH(i)**2,ST(k)**2,scR2)
      enddo
      epsY31=epsY31
     . +VVu*((COCHSUdL(i,1,1)**2+COCHSUdL(i,1,2)**2)
     .              *BB1(0d0,MCH2(i),MSU2(1),scR2)
     .    +2d0*dsqrt(MCH2(i))/vdq
     .   *(U(i,2,1)*COCHSUdL(i,1,1)+U(i,2,2)*COCHSUdL(i,1,2))
     .                      *BB0(0d0,MCH2(i),MSU2(1),scR2))
!     .    +(MQU/vuq)**2*V(i,2)**2
!     .           *BB1(0d0,MCH(i)**2,MUR**2,scR2)
!     .    +(MD/vdq)**2*U(i,2)**2
!     .          *BB1(0d0,MCH(i)**2,MUL**2,scR2))
      enddo
      epsY31=epsY31/(32d0*pi**2)
      
      epsY31I=0d0
      do i=1,2
      do k=1,2
      epsY31I=epsY31I
     .  +2d0*dsqrt(MCH2(i))/MB0*(COCHSTbR(i,k,1)
     .          *COCHSTbL(i,k,2)-COCHSTbR(i,k,2)*COCHSTbL(i,k,1))
     .                     *BB0(0d0,MCH2(i),MST2(k),scR2)
      enddo
      epsY31I=epsY31I
     .    +VVu*(2d0*dsqrt(MCH2(i))/vdq
     .   *(U(i,2,1)*COCHSUdL(i,1,2)-U(i,2,2)*COCHSUdL(i,1,1))
     .                      *BB0(0d0,MCH2(i),MSU2(1),scR2))
      enddo
      epsY31I=epsY31I/(32d0*pi**2)

      epsY13=0d0
      do i=1,2
      do k=1,2
      epsY13=epsY13+2d0*dsqrt(MCH2(i))/MB0*(COCHSTbR(i,k,1)
     .       *COCHSTbL(i,k,1)+COCHSTbR(i,k,2)*COCHSTbL(i,k,2))
     .                       *BB0(0d0,MCH2(i),MST2(k),scR2)
     .  +(COCHSTbL(i,k,1)**2+COCHSTbL(i,k,2)**2
     .           +COCHSTbR(i,k,1)**2+COCHSTbR(i,k,2)**2)
     .                       *BB1(0d0,MCH2(i),MST2(k),scR2)
      enddo
      epsY13=epsY13
     . +VVu*((COCHSUdL(i,1,1)**2+COCHSUdL(i,1,2)**2)
     .                       *BB1(0d0,MCH2(i),MSU2(1),scR2)
!     .    +(MQU/vuq)**2*V(i,2)**2*BB1(0d0,MCH(i)**2,MUR**2,scR2)
     .    +2d0*dsqrt(MCH2(i))/vdq
     .     *(COCHSUdL(i,1,1)*U(i,2,1)+COCHSUdL(i,1,2)*U(i,2,2))
     .                       *BB0(0d0,MCH2(i),MSU2(1),scR2)
     .    +Ybq/vdq*MB0*(U(i,2,1)**2+U(i,2,2)**2)
     .                       *BB1(0d0,MCH2(i),MSU2(1),scR2))
      enddo
      epsY13=epsY13/(32d0*pi**2)

      epsY13I=0d0
      do i=1,2
      do k=1,2
      epsY13I=epsY13I+2d0*dsqrt(MCH2(i))/MB0*(COCHSTbR(i,k,1)
     .       *COCHSTbL(i,k,2)-COCHSTbR(i,k,2)*COCHSTbL(i,k,1))
     .                       *BB0(0d0,MCH2(i),MST2(k),scR2)
      enddo
      epsY13I=epsY13I
     . +VVu*(2d0*dsqrt(MCH2(i))/vdq
     .     *(COCHSUdL(i,1,2)*U(i,2,1)-COCHSUdL(i,1,1)*U(i,2,2))
     .                       *BB0(0d0,MCH2(i),MSU2(1),scR2))
      enddo
      epsY13I=epsY13I/(32d0*pi**2)

      epsY23=0d0
      do i=1,2
      do k=1,2
      epsY23=epsY23+2d0*dsqrt(MCH2(i))/MB0*(COCHSTbR(i,k,1)
     .       *COCHSTbL(i,k,1)+COCHSTbR(i,k,2)*COCHSTbL(i,k,2))
     .                       *BB0(0d0,MCH2(i),MST2(k),scR2)
     .  +(COCHSTbL(i,k,1)**2+COCHSTbL(i,k,2)**2
     .           +COCHSTbR(i,k,1)**2+COCHSTbR(i,k,2)**2)
     .                       *BB1(0d0,MCH2(i),MST2(k),scR2)
      enddo
      epsY23=epsY23
     . +VVc*((COCHSUdL(i,1,1)**2+COCHSUdL(i,1,2)**2)
     .                       *BB1(0d0,MCH2(i),MSU2(1),scR2)
!     .   +(MC0/vuq)**2*V(i,2)**2*BB1(0d0,MCH(i)**2,MUR**2,scR2)
     .   +2d0*dsqrt(MCH2(i))/vdq
     .     *(COCHSUdL(i,1,1)*U(i,2,1)+COCHSUdL(i,1,2)*U(i,2,2))
     .                       *BB0(0d0,MCH2(i),MSU2(1),scR2)
     .    +Ybq/vdq*MB0*(U(i,2,1)**2+U(i,2,2)**2)
     .                       *BB1(0d0,MCH2(i),MSU2(1),scR2))
      enddo
      epsY23=epsY23/(32d0*pi**2)

      epsY23I=0d0
      do i=1,2
      do k=1,2
      epsY23I=epsY23I+2d0*dsqrt(MCH2(i))/MB0*(COCHSTbR(i,k,1)
     .       *COCHSTbL(i,k,2)-COCHSTbR(i,k,2)*COCHSTbL(i,k,1))
     .                       *BB0(0d0,MCH2(i),MST2(k),scR2)
      enddo
      epsY23I=epsY23I
     . +VVc*(2d0*dsqrt(MCH2(i))/vdq
     .     *(COCHSUdL(i,1,2)*U(i,2,1)-COCHSUdL(i,1,2)*U(i,2,1))
     .                       *BB0(0d0,MCH2(i),MSU2(1),scR2))
      enddo
      epsY23I=epsY23I/(32d0*pi**2)

*   * Correction to the CKM matrix: epst0 (*tanb), Eq.(3.53) in [1]
      epst0=epst3-Vtb2*epsY31
      eps0=epst3-Vtb2*epsY32
      epst0I=epst3I-Vtb2*epsY31I
      eps0I=epst3I-Vtb2*epsY32I

*   * Couplings [X^s_RL]^JI as in eqs. (3.55) and (3.56), but
*       WITHOUT the S dependent mixing angles
      sigRLbs=MB0/(vuq*((1d0+eps0)**2+eps0I**2)*(1d0+epst3))
     . *((epsY32*(1d0+eps0)+epsY32I*eps0I)*(1d0-epst3I**2/(1d0+epst3))
     . +(epsY32I*(1d0+eps0)-epsY32*eps0I)*epst3I*(1d0+1d0/(1d0+epst3)))

      sigRLbsI=MB0/(vuq*((1d0+eps0)**2+eps0I**2)*(1d0+epst3))
     . *((epsY32I*(1d0+eps0)-epsY32*eps0I)*(1d0-epst3I**2/(1d0+epst3))
     . +(epsY32*(1d0+eps0)+epsY32I*eps0I)*epst3I*(1d0+1d0/(1d0+epst3)))
      
      sigRLbd=MB0/(vuq*((1d0+epst0)**2+epst0I**2)*(1d0+epst3))
     . *((epsY31*(1d0+epst0)+epsY32I*epst0I)
     .                           *(1d0-epst3I**2/(1d0+epst3))
     .  +(epsY31I*(1d0+epst0)-epsY32*epst0I)
     .                           *epst3I*(1d0+1d0/(1d0+epst3)))
      sigRLbdI=MB0/(vuq*((1d0+epst0)**2+epst0I**2)*(1d0+epst3))
     . *((epsY31I*(1d0+epst0)-epsY31*epst0I)
     .                           *(1d0-epst3I**2/(1d0+epst3))
     .  +(epsY31*(1d0+epst0)+epsY31I*epst0I)
     .                           *epst3I*(1d0+1d0/(1d0+epst3)))

      sigLRbs=MS0/(vuq*((1d0+eps0)**2+eps0I**2)*(1d0+epst2))
     . *((epsY23*(1d0+eps0)-epsY23I*eps0I)*(1d0
     .    +(epsY32*epsY23-epsY32I*epsY23I)/(epsY23**2+epsY23I**2)
     .     *((epst2-epst3)/(1d0+epst3)+epst2I**2/(1d0+epst2))
     .    -(epsY32I*epsY23+epsY23I*epsY32)/(epsY23**2+epsY23I**2)
     .     *(epst3I*(1d0+epst2)/(1d0+epst3)-epst2I/(1d0+epst2)))
     .   +(epsY23I*(1d0+eps0)+epsY23*eps0I)*(epst3I
     .                            +epst2I*(1d0+epst3)/(1d0+epst2)
     .    +(epsY32I*epsY23+epsY23I*epsY32)/(epsY23**2+epsY23I**2)
     .     *((epst2-epst3)/(1d0+epst3)+epst2I**2/(1d0+epst2))
     .    +(epsY32*epsY23-epsY32I*epsY23I)/(epsY23**2+epsY23I**2)
     .     *(epst3I*(1d0+epst2)/(1d0+epst3)-epst2I/(1d0+epst2))))


      sigLRbsI=MS0/(vuq*((1d0+eps0)**2+eps0I**2)*(1d0+epst2))
     . *(-(epsY23I*(1d0+eps0)+epsY23*eps0I)*(1d0
     .    +(epsY32*epsY23-epsY32I*epsY23I)/(epsY23**2+epsY23I**2)
     .     *((epst2-epst3)/(1d0+epst3)+epst2I**2/(1d0+epst2))
     .    -(epsY32I*epsY23+epsY23I*epsY32)/(epsY23**2+epsY23I**2)
     .     *(epst3I*(1d0+epst2)/(1d0+epst3)-epst2I/(1d0+epst2)))
     .   +(epsY23*(1d0+eps0)-epsY23I*eps0I)*(epst3I
     .                            +epst2I*(1d0+epst3)/(1d0+epst2)
     .    +(epsY32I*epsY23+epsY23I*epsY32)/(epsY23**2+epsY23I**2)
     .     *((epst2-epst3)/(1d0+epst3)+epst2I**2/(1d0+epst2))
     .    +(epsY32*epsY23-epsY32I*epsY23I)/(epsY23**2+epsY23I**2)
     .     *(epst3I*(1d0+epst2)/(1d0+epst3)-epst2I/(1d0+epst2))))

      sigLRbd=MD/(vuq*((1d0+epst0)**2+epst0I**2)*(1d0+epst1))
     . *((epsY13*(1d0+epst0)-epsY13I*epst0I)*(1d0
     .    +(epsY31*epsY13-epsY31I*epsY13I)/(epsY13**2+epsY13I**2)
     .     *((epst1-epst3)/(1d0+epst3)+epst1I**2/(1d0+epst1))
     .    -(epsY31I*epsY13+epsY13I*epsY31)/(epsY13**2+epsY13I**2)
     .     *(epst3I*(1d0+epst1)/(1d0+epst3)-epst1I/(1d0+epst1)))
     .   +(epsY13I*(1d0+epst0)+epsY13*epst0I)*(epst3I
     .                            +epst1I*(1d0+epst3)/(1d0+epst1)
     .    +(epsY31I*epsY13+epsY13I*epsY31)/(epsY13**2+epsY13I**2)
     .     *((epst1-epst3)/(1d0+epst3)+epst1I**2/(1d0+epst1))
     .    +(epsY31*epsY13-epsY31I*epsY13I)/(epsY13**2+epsY13I**2)
     .     *(epst3I*(1d0+epst1)/(1d0+epst3)-epst1I/(1d0+epst1))))

      sigLRbdI=MD/(vuq*((1d0+epst0)**2+epst0I**2)*(1d0+epst1))
     . *(-(epsY13I*(1d0+epst0)+epsY13*epst0I)*(1d0
     .    +(epsY31*epsY13-epsY31I*epsY13I)/(epsY13**2+epsY13I**2)
     .     *((epst1-epst3)/(1d0+epst3)+epst1I**2/(1d0+epst1))
     .    -(epsY31I*epsY13+epsY13I*epsY31)/(epsY13**2+epsY13I**2)
     .     *(epst3I*(1d0+epst1)/(1d0+epst3)-epst1I/(1d0+epst1)))
     .   +(epsY13*(1d0+epst0)-epsY13I*epst0I)*(epst3I
     .                            +epst1I*(1d0+epst3)/(1d0+epst1)
     .    +(epsY31I*epsY13+epsY13I*epsY31)/(epsY13**2+epsY13I**2)
     .     *((epst1-epst3)/(1d0+epst3)+epst1I**2/(1d0+epst1))
     .    +(epsY31*epsY13-epsY31I*epsY13I)/(epsY13**2+epsY13I**2)
     .     *(epst3I*(1d0+epst1)/(1d0+epst3)-epst1I/(1d0+epst1))))


*       3) Corrections to the charged Higgs couplings
*          Epsilon_t as in [1] section 5.3
*    * eps_t(s) - gluino/squark contribution
       epsts=0d0
       epstsI=0d0
       do i=1,2
       epsts=epsts-asmsusy*2d0/(3d0*pi)*mu/MGL
     c  *(UT(i,2,1)**2+UT(i,2,2)**2)*DDCOS(PhiM3+Phi01)
     c                      *H2(MST2(i)/MGL**2,MSD2(1)/MGL**2)
       enddo

*               - neutralino/chargino/squark contribution
       MTL2=MST2(1)*(UT(1,1,1)**2+UT(1,1,2)**2)
     c     +MST2(2)*(UT(1,2,1)**2+UT(1,2,2)**2)
       MTR2=MST2(1)*(UT(1,2,1)**2+UT(1,2,2)**2)
     c     +MST2(2)*(UT(2,2,1)**2+UT(2,2,2)**2)
       aux=0d0
       auxe=0d0
       do i=1,2
       do j=1,5
       aux=aux+dsqrt(MNEU(j)/MCH2(i))
     c  *(U(i,2,1)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    -U(i,2,2)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,1)*N(j,4,1)-U(i,1,2)*N(j,4,2)))
     c *(gg2*(V(i,1,1)*N(j,3,1)-V(i,1,2)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MTL2/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c     *gg1*(V(i,1,1)*N(j,1,1)-V(i,1,2)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MTR2/MCH2(i))
     c  -((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,1)
     c    -(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,2))/dsqrt(2d0)
     c                    *H2(MNEU(j)/MCH2(i),MSD2(1)/MCH2(i)))
     c  -dsqrt(MNEU(j)/MCH2(i))
     c   *(U(i,2,2)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    +U(i,2,1)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,2)*N(j,4,1)+U(i,1,1)*N(j,4,2)))
     c *(gg2*(V(i,1,2)*N(j,3,1)+V(i,1,1)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MTL2/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c   *gg1*(V(i,1,2)*N(j,1,1)+V(i,1,1)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MTR2/MCH2(i))
     c  -((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,2)/dsqrt(2d0)
     c     +(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,1)/dsqrt(2d0))
     c                    *H2(MNEU(j)/MCH2(i),MSD2(1)/MCH2(i)))
       auxe=auxe-dsqrt(MNEU(j)/MCH2(i))
     c     *(U(i,2,2)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    +U(i,2,1)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,2)*N(j,4,1)+U(i,1,1)*N(j,4,2)))
     c *(gg2*(V(i,1,1)*N(j,3,1)-V(i,1,2)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MTL2/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c     *gg1*(V(i,1,1)*N(j,1,1)-V(i,1,2)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MTR2/MCH2(i))
     c  +((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,1)
     c    -(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,2))/dsqrt(2d0)
     c                    *H2(MNEU(j)/MCH2(i),MSD2(1)/MCH2(i)))
     c  -dsqrt(MNEU(j)/MCH2(i))
     c   *(U(i,2,1)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    -U(i,2,2)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,1)*N(j,4,1)-U(i,1,2)*N(j,4,2)))
     c *(gg2*(V(i,1,2)*N(j,3,1)+V(i,1,1)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MTL2/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c   *gg1*(V(i,1,2)*N(j,1,1)+V(i,1,1)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MTR2/MCH2(i))
     c  -((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,2)/dsqrt(2d0)
     c     +(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,1)/dsqrt(2d0))
     c                    *H2(MNEU(j)/MCH2(i),MSD2(1)/MCH2(i)))
       enddo
       enddo
       epsts=epsts+aux/16d0/Pi**2
       epstsI=epstsI+auxe/16d0/Pi**2

*    * eps_t(b) - gluino/squark contribution
       epstb=0d0
       epstbI=0d0
       do i=1,2
       do k=1,2
       epstb=epstb-asmsusy*2d0/(3d0*pi)*mu/MGL*DDCOS(PhiM3+Phi01)
     c  *(UT(i,2,1)**2+UT(i,2,2)**2)*(UB(k,1,1)**2+UB(k,1,2)**2)
     c                    *H2(MST2(i)/MGL**2,MSB2(k)/MGL**2)
       epstbI=epstbI+asmsusy*2d0/(3d0*pi)*mu/MGL*DDSIN(PhiM3+Phi01)
     c  *(UT(i,2,1)**2+UT(i,2,2)**2)*(UB(k,1,1)**2+UB(k,1,2)**2)
     c                    *H2(MST2(i)/MGL**2,MSB2(k)/MGL**2)
       enddo
       enddo

*               - neutralino/squark contribution
       aux=0d0
       auxe=0d0
       do i=1,5
       do j=1,2
       do k=1,2
       aux=aux+1d0/dsqrt(MNEU(i))*
     .  ((N(i,4,1)*N(i,3,1)-N(i,4,2)*N(i,3,2))*DDCOS(phiAb)
     .  -(N(i,4,2)*N(i,3,1)+N(i,4,1)*N(i,3,2))*DDSIN(phiAb))
     .  *(UT(j,1,1)**2+UT(j,1,2)**2)*(UB(k,1,1)**2+UB(k,1,2)**2)
     .                 *H2(MST2(j)/MNEU(i),MSB2(k)/MNEU(i))
       auxe=auxe+1d0/dsqrt(MNEU(i))*
     .  ((N(i,4,1)*N(i,3,1)-N(i,4,2)*N(i,3,2))*DDSIN(phiAb)
     .  +(N(i,4,2)*N(i,3,1)+N(i,4,1)*N(i,3,2))*DDCOS(phiAb))
     .  *(UT(j,1,1)**2+UT(j,1,2)**2)*(UB(k,1,1)**2+UB(k,1,2)**2)
     .                 *H2(MST2(j)/MNEU(i),MSB2(k)/MNEU(i))
       enddo
       enddo
       enddo
        epstb=epstb-Ybq**2/(16d0*pi**2)*AAB*aux
        epstbI=epstbI+Ybq**2/(16d0*pi**2)*AAB*auxe

*               - neutralino/chargino/squark contribution
       MBL2=MSB2(1)*(UB(1,1,1)**2+UB(1,1,2)**2)
     .     +MSB2(2)*(UB(2,1,1)**2+UB(2,1,2)**2)
       aux=0d0
       auxe=0d0
       do i=1,2
       do j=1,5
       aux=aux+dsqrt(MNEU(j)/MCH2(i))
     c  *(U(i,2,1)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    -U(i,2,2)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,1)*N(j,4,1)-U(i,1,2)*N(j,4,2)))
     c *(gg2*(V(i,1,1)*N(j,3,1)-V(i,1,2)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MTL2/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c     *gg1*(V(i,1,1)*N(j,1,1)-V(i,1,2)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MTR2/MCH2(i))
     c  -((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,1)
     c    -(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,2))/dsqrt(2d0)
     c                    *H2(MNEU(j)/MCH2(i),MBL2/MCH2(i)))
     c  -dsqrt(MNEU(j)/MCH2(i))
     c   *(U(i,2,2)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    +U(i,2,1)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,2)*N(j,4,1)+U(i,1,1)*N(j,4,2)))
     c *(gg2*(V(i,1,2)*N(j,3,1)+V(i,1,1)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MTL2/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c   *gg1*(V(i,1,2)*N(j,1,1)+V(i,1,1)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MTR2/MCH2(i))
     c  -((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,2)/dsqrt(2d0)
     c     +(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,1)/dsqrt(2d0))
     c                    *H2(MNEU(j)/MCH2(i),MBL2/MCH2(i)))
       auxe=auxe-dsqrt(MNEU(j)/MCH2(i))
     c     *(U(i,2,2)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    +U(i,2,1)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,2)*N(j,4,1)+U(i,1,1)*N(j,4,2)))
     c *(gg2*(V(i,1,1)*N(j,3,1)-V(i,1,2)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MTL2/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c     *gg1*(V(i,1,1)*N(j,1,1)-V(i,1,2)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MTR2/MCH2(i))
     c  +((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,1)
     c    -(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,2))/dsqrt(2d0)
     c                    *H2(MNEU(j)/MCH2(i),MBL2/MCH2(i)))
     c  -dsqrt(MNEU(j)/MCH2(i))
     c   *(U(i,2,1)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    -U(i,2,2)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,1)*N(j,4,1)-U(i,1,2)*N(j,4,2)))
     c *(gg2*(V(i,1,2)*N(j,3,1)+V(i,1,1)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MTL2/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c   *gg1*(V(i,1,2)*N(j,1,1)+V(i,1,1)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MTR2/MCH2(i))
     c  -((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,2)/dsqrt(2d0)
     c     +(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,1)/dsqrt(2d0))
     c                    *H2(MNEU(j)/MCH2(i),MBL2/MCH2(i)))
       enddo
       enddo
       epstb=epstb+aux/16d0/Pi**2
       epstbI=epstbI+auxe/16d0/Pi**2

*    * eps_c(s) - gluino/squark contribution
       epscs=-asmsusy*2d0/(3d0*pi)*mu/MGL*DDCOS(PhiM3+Phi01)
     c                      *H2(MSU2(1)/MGL**2,MSD2(1)/MGL**2)
       epscsI=asmsusy*2d0/(3d0*pi)*mu/MGL*DDSIN(PhiM3+Phi01)
     c                      *H2(MSU2(1)/MGL**2,MSD2(1)/MGL**2)

*               - neutralino/chargino/squark contribution
       aux=0d0
       auxe=0d0
       do i=1,2
       do j=1,5
       aux=aux+dsqrt(MNEU(j)/MCH2(i))
     c  *(U(i,2,1)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    -U(i,2,2)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,1)*N(j,4,1)-U(i,1,2)*N(j,4,2)))
     c *(gg2*(V(i,1,1)*N(j,3,1)-V(i,1,2)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(1)/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c     *gg1*(V(i,1,1)*N(j,1,1)-V(i,1,2)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(2)/MCH2(i))
     c  -((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,1)
     c    -(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,2))/dsqrt(2d0)
     c                    *H2(MNEU(j)/MCH2(i),MSD2(1)/MCH2(i)))
     c  -dsqrt(MNEU(j)/MCH2(i))
     c   *(U(i,2,2)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    +U(i,2,1)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,2)*N(j,4,1)+U(i,1,1)*N(j,4,2)))
     c *(gg2*(V(i,1,2)*N(j,3,1)+V(i,1,1)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(1)/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c   *gg1*(V(i,1,2)*N(j,1,1)+V(i,1,1)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(2)/MCH2(i))
     c  -((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,2)/dsqrt(2d0)
     c     +(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,1)/dsqrt(2d0))
     c                    *H2(MNEU(j)/MCH2(i),MSD2(1)/MCH2(i)))
       auxe=auxe-dsqrt(MNEU(j)/MCH2(i))
     c     *(U(i,2,2)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    +U(i,2,1)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,2)*N(j,4,1)+U(i,1,1)*N(j,4,2)))
     c *(gg2*(V(i,1,1)*N(j,3,1)-V(i,1,2)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(1)/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c     *gg1*(V(i,1,1)*N(j,1,1)-V(i,1,2)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(2)/MCH2(i))
     c  +((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,1)
     c    -(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,2))/dsqrt(2d0)
     c                    *H2(MNEU(j)/MCH2(i),MSD2(1)/MCH2(i)))
     c  -dsqrt(MNEU(j)/MCH2(i))
     c   *(U(i,2,1)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    -U(i,2,2)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,1)*N(j,4,1)-U(i,1,2)*N(j,4,2)))
     c *(gg2*(V(i,1,2)*N(j,3,1)+V(i,1,1)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(1)/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c   *gg1*(V(i,1,2)*N(j,1,1)+V(i,1,1)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(2)/MCH2(i))
     c  -((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,2)/dsqrt(2d0)
     c     +(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,1)/dsqrt(2d0))
     c                    *H2(MNEU(j)/MCH2(i),MSD2(1)/MCH2(i)))
       enddo
       enddo
       epscs=epscs+aux/16d0/Pi**2
       epscsI=epscsI+auxe/16d0/Pi**2

*    * eps_c(b) - gluino/squark contribution
       epscb=0d0
       epscbI=0d0
       do k=1,2
       epscb=epscb-asmsusy*2d0/(3d0*pi)*mu/MGL*DDCOS(PhiM3+Phi01)
     c         *(UB(k,1,1)**2+UB(k,1,2)**2)
     c                      *H2(MSU2(1)/MGL**2,MSB2(k)/MGL**2)
       epscbI=epscbI+asmsusy*2d0/(3d0*pi)*mu/MGL*DDSIN(PhiM3+Phi01)
     c         *(UB(k,1,1)**2+UB(k,1,2)**2)
     c                      *H2(MSU2(1)/MGL**2,MSB2(k)/MGL**2)
       enddo

*               - neutralino/squark contribution
       aux=0d0
       auxe=0d0
       do i=1,5
       do k=1,2
       aux=aux+(UB(k,1,1)**2+UB(k,1,2)**2)/dsqrt(MNEU(i))
     .  *((N(i,4,1)*N(i,3,1)-N(i,4,2)*N(i,3,2))*DDCOS(phiAb)
     .   -(N(i,4,2)*N(i,3,1)+N(i,4,1)*N(i,3,2))*DDSIN(phiAb))
     .                 *H2(MSU2(1)/MNEU(i),MSB2(k)/MNEU(i))
       auxe=auxe+(UB(k,1,1)**2+UB(k,1,2)**2)/dsqrt(MNEU(i))
     .  *((N(i,4,1)*N(i,3,1)-N(i,4,2)*N(i,3,2))*DDSIN(phiAb)
     .   +(N(i,4,2)*N(i,3,1)+N(i,4,1)*N(i,3,2))*DDCOS(phiAb))
     .                 *H2(MSU2(1)/MNEU(i),MSB2(k)/MNEU(i))
       enddo
       enddo
        epscb=epscb-Ybq**2/(16d0*pi**2)*AAB*aux
        epscbI=epscbI+Ybq**2/(16d0*pi**2)*AAB*auxe

*               - neutralino/chargino/squark contribution
       aux=0d0
       auxe=0d0
       do i=1,2
       do j=1,5
       aux=aux+dsqrt(MNEU(j)/MCH2(i))
     c  *(U(i,2,1)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    -U(i,2,2)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,1)*N(j,4,1)-U(i,1,2)*N(j,4,2)))
     c *(gg2*(V(i,1,1)*N(j,3,1)-V(i,1,2)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(1)/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c     *gg1*(V(i,1,1)*N(j,1,1)-V(i,1,2)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(2)/MCH2(i))
     c  -((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,1)
     c    -(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,2))/dsqrt(2d0)
     c                    *H2(MNEU(j)/MCH2(i),MBL2/MCH2(i)))
     c  -dsqrt(MNEU(j)/MCH2(i))
     c   *(U(i,2,2)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    +U(i,2,1)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,2)*N(j,4,1)+U(i,1,1)*N(j,4,2)))
     c *(gg2*(V(i,1,2)*N(j,3,1)+V(i,1,1)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(1)/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c   *gg1*(V(i,1,2)*N(j,1,1)+V(i,1,1)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(2)/MCH2(i))
     c  -((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,2)/dsqrt(2d0)
     c     +(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,1)/dsqrt(2d0))
     c                    *H2(MNEU(j)/MCH2(i),MBL2/MCH2(i)))
       auxe=auxe-dsqrt(MNEU(j)/MCH2(i))
     c     *(U(i,2,2)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    +U(i,2,1)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,2)*N(j,4,1)+U(i,1,1)*N(j,4,2)))
     c *(gg2*(V(i,1,1)*N(j,3,1)-V(i,1,2)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(1)/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c     *gg1*(V(i,1,1)*N(j,1,1)-V(i,1,2)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(2)/MCH2(i))
     c  +((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,1)
     c    -(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,2))/dsqrt(2d0)
     c                    *H2(MNEU(j)/MCH2(i),MBL2/MCH2(i)))
     c  -dsqrt(MNEU(j)/MCH2(i))
     c   *(U(i,2,1)*(gg1*N(j,1,1)+gg2*N(j,2,1))/dsqrt(2d0)
     c    -U(i,2,2)*(gg1*N(j,1,2)+gg2*N(j,2,2))/dsqrt(2d0)
     c    -gg2*(U(i,1,1)*N(j,4,1)-U(i,1,2)*N(j,4,2)))
     c *(gg2*(V(i,1,2)*N(j,3,1)+V(i,1,1)*N(j,3,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(1)/MCH2(i))
     c  +2d0*dsqrt(2d0)/3d0
     c   *gg1*(V(i,1,2)*N(j,1,1)+V(i,1,1)*N(j,1,2))
     c                    *H2(MNEU(j)/MCH2(i),MSU2(2)/MCH2(i))
     c  -((gg1*N(j,1,1)/3d0-gg2*N(j,2,1))*V(i,2,2)/dsqrt(2d0)
     c     +(gg1*N(j,1,2)/3d0-gg2*N(j,2,2))*V(i,2,1)/dsqrt(2d0))
     c                    *H2(MNEU(j)/MCH2(i),MBL2/MCH2(i)))
       enddo
       enddo
       epscb=epscb+aux/16d0/Pi**2
       epscbI=epscbI+auxe/16d0/Pi**2


*	III- Delt B = 2 TRANSITIONS: DMs = m_Bs-m_Bbar_s 
*             and DMd = m_Bd-m_Bbar_d (following [24], [1])

*  Matching scale
      sc0=166d0
      asc0=asf(sc0)

*  Low energy scale
      scb=4.6d0
      asmb=asf(scb)

*	 1) Wilson coefficients for DMs at the matching scale
*          - Standard Model [1] Eq.(6.7)
      xt=(MT0*(asf(sc0)/asf(MT0))**(12d0/23d0)
     .   *(1d0+7462d0/1587d0*(asf(sc0)-asf(MT0))/4d0/Pi)/MW)**2 !(m_t/m_W)^2

      CVLLSM=4d0*S0(xt)*0.985d0
      
*          - Charged Higgs Boxes [1] App.(A.4.i)
*  Coefficients at MHC
      MCHH(1)=MW
      MCHH(2)=dsqrt(MHC)

      aux=mth/dsqrt(vuq**2+vdq**2)
      PRLH_tb(1)=-aux
      PRLHI_tb(1)=0d0
      PRLH_tb(2)=aux/tanb*(1d0-epstb*(1d0/tanb+tanb))
      PRLHI_tb(2)=-aux*epstbI*(1d0/tanb**2+1d0)
      PRLH_ts(1)=-aux
      PRLHI_ts(1)=0d0
      PRLH_ts(2)=aux/tanb*(1d0-(1d0/tanb+tanb)
     . *(epsts-((epsY32*(1d0+eps0)+epsY32I*eps0I)*(epstb-epsts)
     .         -(epsY32I*(1d0+eps0)-epsY32*eps0I)*(epstbI-epstsI))
     .                                    /((1d0+eps0)**2+eps0I**2)))
      PRLHI_ts(2)=-aux*(1d0/tanb**2+1d0)
     . *(epstsI-((epsY32I*(1d0+eps0)-epsY32*eps0I)*(epstb-epsts)
     .         +(epsY32*(1d0+eps0)+epsY32I*eps0I)*(epstbI-epstsI))
     .                                    /((1d0+eps0)**2+eps0I**2))

      aux=mbh/dsqrt(vuq**2+vdq**2)
      PLRH_tb(1)=aux
      PLRHI_tb(1)=0d0
      PLRH_tb(2)=aux*(tanb-(epst3+epst3I**2/(1d0+epst3))/(1d0+epst3)
     .                                        *(tanb+1d0/tanb))
      PLRHI_tb(2)=aux*(tanb+1d0/tanb)*epst3I/(1d0+epst3)**2
      aux=runmass(0.095d0,dsqrt(MHC))/dsqrt(vuq**2+vdq**2)
      PLRH_ts(1)=aux
      PLRHI_ts(1)=0d0
      PLRH_ts(2)=aux*(tanb-(tanb+1d0/tanb)/(1d0+epst2)
     .   *(epst2+epst2I**2/(1d0+epst2)
     . +((1d0+eps0)*(epsY23+epsY23I
     .                       *(epst3+epst2*(1d0+epst3)/(1d0+epst2))
     .               +epsY32*(epst2-epst3)-epsY32I
     .         *(epst3I*(1d0+epst2)-epst2I*(1d0+epst3)/(1d0+epst2)))
     .  -eps0I*(epsY23I-epsY23*(epst3I+epst2I*(1d0+epst3)/(1d0+epst2))
     .         -epsY32I*(epst2-epst3)-epsY32*(epst3I*(1d0+epst2)
     .           -epst2I*(1d0+epst3)/(1d0+epst2))))
     .               /((1d0+eps0)**2+eps0I**2)))
      PLRHI_ts(2)=aux*((tanb+1d0/tanb)/(1d0+epst2)
     .   *(epst2I/(1d0+epst2)
     . -(eps0I*(epsY23+epsY23I*(epst3+epst2*(1d0+epst3)/(1d0+epst2))
     .               +epsY32*(epst2-epst3)-epsY32I
     .         *(epst3I*(1d0+epst2)-epst2I*(1d0+epst3)/(1d0+epst2)))
     .  -(1d0+eps0)*(epsY23I-epsY23*(epst3I
     .                              +epst2I*(1d0+epst3)/(1d0+epst2))
     .         -epsY32I*(epst2-epst3)-epsY32*(epst3I*(1d0+epst2)
     .           -epst2I*(1d0+epst3)/(1d0+epst2))))
     .               /((1d0+eps0)**2+eps0I**2)))

      CVLLHIG=-g2q/2d0*(PRLH_tb(2)*PRLH_ts(2)+PRLHI_tb(2)*PRLHI_ts(2))
     .                          *mth**2*D0B(MW,MCHH(2),mth,mth)
     . +((PRLH_tb(2)*PRLH_tb(2)-PRLHI_tb(2)*PRLHI_tb(2))
     .   *(PRLH_ts(2)*PRLH_ts(2)-PRLHI_ts(2)*PRLHI_ts(2))
     .  +(PRLH_tb(2)*PRLHI_tb(2)+PRLHI_tb(2)*PRLH_tb(2))
     .   *(PRLH_ts(2)*PRLHI_ts(2)+PRLHI_ts(2)*PRLH_ts(2)))
     .              *D2B(MCHH(2),MCHH(2),mth,mth)/8d0
     . +((PRLH_tb(2)*PRLH_tb(1)-PRLHI_tb(2)*PRLHI_tb(1))
     .   *(PRLH_ts(2)*PRLH_ts(1)-PRLHI_ts(2)*PRLHI_ts(1))
     .  +(PRLH_tb(2)*PRLHI_tb(1)+PRLHI_tb(2)*PRLH_tb(1))
     .   *(PRLH_ts(2)*PRLHI_ts(1)+PRLHI_ts(2)*PRLH_ts(1)))
     .              *D2B(MCHH(1),MCHH(2),mth,mth)/4d0
      CVLLHIG=CVLLHIG/(GF*MW)**2

      CVLLHIGI=-g2q/2d0*(PRLH_tb(2)*PRLHI_ts(2)-PRLHI_tb(2)*PRLH_ts(2))
     .                          *mth**2*D0B(MW,MCHH(2),mth,mth)
     . +((PRLH_tb(2)*PRLH_tb(2)-PRLHI_tb(2)*PRLHI_tb(2))
     .   *(PRLH_ts(2)*PRLHI_ts(2)+PRLHI_ts(2)*PRLH_ts(2))
     .  -(PRLH_tb(2)*PRLHI_tb(2)+PRLHI_tb(2)*PRLH_tb(2))
     .   *(PRLH_ts(2)*PRLH_ts(2)-PRLHI_ts(2)*PRLHI_ts(2)))
     .              *D2B(MCHH(2),MCHH(2),mth,mth)/8d0
     . +((PRLH_tb(2)*PRLH_tb(1)-PRLHI_tb(2)*PRLHI_tb(1))
     .   *(PRLH_ts(2)*PRLHI_ts(1)+PRLHI_ts(2)*PRLH_ts(1))
     .  -(PRLH_tb(2)*PRLHI_tb(1)+PRLHI_tb(2)*PRLH_tb(1))
     .   *(PRLH_ts(2)*PRLH_ts(1)-PRLHI_ts(2)*PRLHI_ts(1)))
     .              *D2B(MCHH(1),MCHH(2),mth,mth)/4d0
      CVLLHIGI=CVLLHIGI/(GF*MW)**2

      C1SLLHIG=0d0
      DO I=1,2
      DO J=1,2
       C1SLLHIG=C1SLLHIG
     .  +((PLRH_tb(J)*PLRH_tb(I)-PLRHI_tb(J)*PLRHI_tb(I))
     .   *(PRLH_ts(I)*PRLH_ts(J)-PRLHI_ts(I)*PRLHI_ts(J))
     .  +(PLRH_tb(J)*PLRHI_tb(I)+PLRHI_tb(J)*PLRH_tb(I))
     .   *(PRLH_ts(I)*PRLHI_ts(J)+PRLHI_ts(I)*PRLH_ts(J)))
     .              *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)/2d0
      ENDDO
      ENDDO
      C1SLLHIG=C1SLLHIG/(GF*MW)**2

      C1SLLHIGI=0d0
      DO I=1,2
      DO J=1,2
       C1SLLHIGI=C1SLLHIGI
     .  +((PLRH_tb(J)*PLRH_tb(I)-PLRHI_tb(J)*PLRHI_tb(I))
     .   *(PRLH_ts(I)*PRLHI_ts(J)+PRLHI_ts(I)*PRLH_ts(J))
     .  -(PLRH_tb(J)*PLRHI_tb(I)+PLRHI_tb(J)*PLRH_tb(I))
     .   *(PRLH_ts(I)*PRLH_ts(J)-PRLHI_ts(I)*PRLHI_ts(J)))
     .              *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)/2d0
      ENDDO
      ENDDO
      C1SLLHIGI=C1SLLHIGI/(GF*MW)**2

      C1LRHIG=0d0
      DO I=1,2
      DO J=1,2
       C1LRHIG=C1LRHIG
     .  +((PRLH_tb(J)*PLRH_tb(I)-PRLHI_tb(J)*PLRHI_tb(I))
     .   *(PRLH_ts(I)*PLRH_ts(J)-PRLHI_ts(I)*PLRHI_ts(J))
     .  +(PRLH_tb(J)*PLRHI_tb(I)+PRLHI_tb(J)*PLRH_tb(I))
     .   *(PRLH_ts(I)*PLRHI_ts(J)+PRLHI_ts(I)*PLRH_ts(J)))
     .                          *D2B(MCHH(I),MCHH(J),mth,mth)/4d0
      ENDDO
      ENDDO
      C1LRHIG=C1LRHIG/(GF*MW)**2

      C1LRHIGI=0d0
      DO I=1,2
      DO J=1,2
       C1LRHIGI=C1LRHIGI
     .  +((PRLH_tb(J)*PLRH_tb(I)-PRLHI_tb(J)*PLRHI_tb(I))
     .   *(PRLH_ts(I)*PLRHI_ts(J)+PRLHI_ts(I)*PLRH_ts(J))
     .  -(PRLH_tb(J)*PLRHI_tb(I)+PRLHI_tb(J)*PLRH_tb(I))
     .   *(PRLH_ts(I)*PLRH_ts(J)-PRLHI_ts(I)*PLRHI_ts(J)))
     .                          *D2B(MCHH(I),MCHH(J),mth,mth)/4d0
      ENDDO
      ENDDO
      C1LRHIGI=C1LRHIGI/(GF*MW)**2

      aux=(1d0+epst3)/(1d0+epst2)
      C2LRHIG=0d0
      DO I=1,2
       C2LRHIG=C2LRHIG
     .  -g2q*(PLRH_tb(I)*PLRH_ts(I)+PLRHI_tb(I)*PLRHI_ts(I))/2d0
     .    *(D2B(MCHH(I),MW,mth,mth)
     .      -2d0*aux*D2B(MCHH(I),MW,mth,0d0)
     .      +aux**2*D2B(MCHH(I),MW,0d0,0d0))
      DO J=1,2
       C2LRHIG=C2LRHIG
     .  +((PLRH_tb(J)*PRLH_tb(I)-PLRHI_tb(J)*PRLHI_tb(I))
     .   *(PRLH_ts(I)*PLRH_ts(J)-PRLHI_ts(I)*PLRHI_ts(J))
     .  +(PLRH_tb(J)*PRLHI_tb(I)+PLRHI_tb(J)*PRLH_tb(I))
     .   *(PRLH_ts(I)*PLRHI_ts(J)+PRLHI_ts(I)*PLRH_ts(J)))
     .                     *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)
      ENDDO
      ENDDO
      C2LRHIG=C2LRHIG/(GF*MW)**2

      C2LRHIGI=0d0
      DO I=1,2
       C2LRHIGI=C2LRHIGI
     .  -g2q*(PLRH_tb(I)*PLRHI_ts(I)-PLRHI_tb(I)*PLRH_ts(I))/2d0
     .    *(D2B(MCHH(I),MW,mth,mth)
     .      -2d0*aux*D2B(MCHH(I),MW,mth,0d0)
     .      +aux**2*D2B(MCHH(I),MW,0d0,0d0))
      DO J=1,2
       C2LRHIGI=C2LRHIGI
     .  +((PLRH_tb(J)*PRLH_tb(I)-PLRHI_tb(J)*PRLHI_tb(I))
     .   *(PRLH_ts(I)*PLRHI_ts(J)+PRLHI_ts(I)*PLRH_ts(J))
     .  -(PLRH_tb(J)*PRLHI_tb(I)+PLRHI_tb(J)*PRLH_tb(I))
     .   *(PRLH_ts(I)*PLRH_ts(J)-PRLHI_ts(I)*PLRHI_ts(J)))
     .                     *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)
      ENDDO
      ENDDO
      C2LRHIGI=C2LRHIGI/(GF*MW)**2

*  Running to sc0 [24] App.C
      etaH=asf(dsqrt(MHC))/asc0

      CVLLHIG=etaH**(6d0/21d0)*(1d0+asc0/4d0/Pi*1.3707d0*(1d0-etaH))
     .                 *CVLLHIG
      CVLLHIGI=etaH**(6d0/21d0)*(1d0+asc0/4d0/Pi*1.3707d0*(1d0-etaH))
     .                 *CVLLHIGI

      aux=C1SLLHIG
      C1SLLHIG=(1.0153d0*etaH**(-0.6916d0)-0.0153d0*etaH**(0.7869d0)
     . +asc0/4d0/Pi*(etaH**(-0.6916d0)*(5.6478d0-6.0350d0*etaH)
     .              +etaH**(0.7869d0)*(0.3272d0+0.0600*etaH)))*aux

      C2SLLHIG=(0.0081d0*(etaH**(0.7869d0)-etaH**(-0.6916d0))
     . +asc0/4d0/Pi*(etaH**(0.7869d0)*(-0.0618d0-0.0315d0*etaH)
     .              +etaH**(-0.6916d0)*(0.0454d0+0.0479d0*etaH)))*aux

      aux=C1SLLHIGI
      C1SLLHIGI=(1.0153d0*etaH**(-0.6916d0)-0.0153d0*etaH**(0.7869d0)
     . +asc0/4d0/Pi*(etaH**(-0.6916d0)*(5.6478d0-6.0350d0*etaH)
     .              +etaH**(0.7869d0)*(0.3272d0+0.0600*etaH)))*aux

      C2SLLHIGI=(0.0081d0*(etaH**(0.7869d0)-etaH**(-0.6916d0))
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

*          - Chargino / squark Boxes [1] App.(A.4.ii)
*  Coefficients at scR
      
      aux=0d0                   !3rd family contribution
      auxe=0d0
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+((COCHSTbL(j,k,1)*COCHSTbL(i,l,1)
     .                           -COCHSTbL(j,k,2)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbL(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbL(j,l,2))
     .    +(COCHSTbL(j,k,2)*COCHSTbL(i,l,1)
     .                           +COCHSTbL(j,k,1)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbL(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbL(j,l,2)))
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
       auxe=auxe+((COCHSTbL(j,k,1)*COCHSTbL(i,l,1)
     .                           -COCHSTbL(j,k,2)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbL(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbL(j,l,2))
     .    -(COCHSTbL(j,k,2)*COCHSTbL(i,l,1)
     .                           +COCHSTbL(j,k,1)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbL(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbL(j,l,2)))
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
     .      
      enddo                     !3rd family/1st-2nd interference
       aux=aux-2d0*g2q*((V(j,1,1)*V(i,1,1)+V(j,1,2)*V(i,1,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbL(j,k,1)
     .                           +COCHSTbL(i,k,2)*COCHSTbL(j,k,2))
     .     +(V(j,1,2)*V(i,1,1)-V(j,1,1)*V(i,1,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbL(j,k,2)
     .                           -COCHSTbL(i,k,2)*COCHSTbL(j,k,1)))
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
       auxe=auxe-2d0*g2q*((V(j,1,1)*V(i,1,1)+V(j,1,2)*V(i,1,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbL(j,k,1)
     .                           -COCHSTbL(i,k,1)*COCHSTbL(j,k,2))
     .     +(V(j,1,2)*V(i,1,1)-V(j,1,1)*V(i,1,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbL(j,k,1)
     .                           +COCHSTbL(i,k,2)*COCHSTbL(j,k,2)))
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
      enddo                     !1st-2nd family contribution
       aux=aux+g2q**2*((V(j,1,1)*V(i,1,1)+V(j,1,2)*V(i,1,2))**2
     .                +(V(j,1,2)*V(i,1,1)-V(j,1,1)*V(i,1,2))**2)
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MSU2(1)),dsqrt(MSU2(1)))
      enddo
      enddo
      CVLLCHAR=aux/8d0/(GF*MW)**2
      CVLLCHARI=auxe/8d0/(GF*MW)**2

      aux=0d0                   !3rd family contribution
      auxe=0d0
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+((COCHSTbR(j,k,1)*COCHSTbR(i,l,1)
     .                           -COCHSTbR(j,k,2)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbL(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbL(j,l,2))
     .    +(COCHSTbR(j,k,2)*COCHSTbR(i,l,1)
     .                           +COCHSTbR(j,k,1)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbL(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbL(j,l,2)))
     .                  *dsqrt(MCH2(i)*MCH2(j))/(1d0+epst3)**2
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
       auxe=auxe+((COCHSTbR(j,k,1)*COCHSTbR(i,l,1)
     .                           -COCHSTbR(j,k,2)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbL(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbL(j,l,2))
     .    -(COCHSTbR(j,k,2)*COCHSTbR(i,l,1)
     .                           +COCHSTbR(j,k,1)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbL(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbL(j,l,2)))
     .                  *dsqrt(MCH2(i)*MCH2(j))/(1d0+epst3)**2
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
      enddo                     !3rd family/1st-2nd interference
       aux=aux-2d0*((COCHSTbR(j,k,1)*U(i,2,1)-COCHSTbR(j,k,2)*U(i,2,2))
     .    *(COCHSTbL(i,k,1)*COCHSUdL(j,1,1)
     .                           -COCHSTbL(i,k,2)*COCHSUdL(j,1,2))
     .    +(COCHSTbR(j,k,2)*U(i,2,1)+COCHSTbR(j,k,1)*U(i,2,2))
     .    *(COCHSTbL(i,k,2)*COCHSUdL(j,1,1)
     .                           +COCHSTbL(i,k,1)*COCHSUdL(j,1,2)))
     .         *Ybq/(1+epst3)**2*dsqrt(MCH2(i)*MCH2(j))
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
       auxe=auxe
     .   -2d0*((COCHSTbR(j,k,1)*U(i,2,1)-COCHSTbR(j,k,2)*U(i,2,2))
     .    *(COCHSTbL(i,k,2)*COCHSUdL(j,1,1)
     .                           +COCHSTbL(i,k,1)*COCHSUdL(j,1,2))
     .    -(COCHSTbR(j,k,2)*U(i,2,1)+COCHSTbR(j,k,1)*U(i,2,2))
     .    *(COCHSTbL(i,k,1)*COCHSUdL(j,1,1)
     .                           -COCHSTbL(i,k,2)*COCHSUdL(j,1,2)))
     .         *Ybq/(1+epst3)**2*dsqrt(MCH2(i)*MCH2(j))
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
      enddo                     !1st-2nd family contribution
       aux=aux+((U(j,2,1)*U(i,2,1)-U(j,2,2)*U(i,2,2))
     .    *(COCHSUdL(i,1,1)*COCHSUdL(j,1,1)
     .                           -COCHSUdL(i,1,2)*COCHSUdL(j,1,2))
     .    +(U(j,2,2)*U(i,2,1)+U(j,2,1)*U(i,2,2))
     .    *(COCHSUdL(i,1,2)*COCHSUdL(j,1,1)
     .                           +COCHSUdL(i,1,1)*COCHSUdL(j,1,2)))
     .         *Ybq**2/(1+epst3)**2*dsqrt(MCH2(i)*MCH2(j))
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MSU2(1)),dsqrt(MSU2(1)))
       auxe=auxe+((U(j,2,1)*U(i,2,1)-U(j,2,2)*U(i,2,2))
     .    *(COCHSUdL(i,1,2)*COCHSUdL(j,1,1)
     .                           +COCHSUdL(i,1,1)*COCHSUdL(j,1,2))
     .    -(U(j,2,2)*U(i,2,1)+U(j,2,1)*U(i,2,2))
     .    *(COCHSUdL(i,1,1)*COCHSUdL(j,1,1)
     .                           -COCHSUdL(i,1,2)*COCHSUdL(j,1,2)))
     .         *Ybq**2/(1+epst3)**2*dsqrt(MCH2(i)*MCH2(j))
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MSU2(1)),dsqrt(MSU2(1)))
      enddo
      enddo
      C1SLLCHAR=-aux/4d0/(GF*MW)**2
      C2SLLCHAR=aux/16d0/(GF*MW)**2
      C1SLLCHARI=-auxe/4d0/(GF*MW)**2
      C2SLLCHARI=auxe/16d0/(GF*MW)**2

      aux=0d0                   !3rd family contribution
      auxe=0d0
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+((COCHSTbL(j,k,1)*COCHSTbR(i,l,1)
     .                           -COCHSTbL(j,k,2)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbR(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbR(j,l,2))
     .    +(COCHSTbL(j,k,2)*COCHSTbR(i,l,1)
     .                           +COCHSTbL(j,k,1)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbR(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbR(j,l,2)))
     .                  *dsqrt(MCH2(i)*MCH2(j))/(1d0+epst3)**2
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
       auxe=auxe+((COCHSTbL(j,k,1)*COCHSTbR(i,l,1)
     .                           -COCHSTbL(j,k,2)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbR(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbR(j,l,2))
     .    -(COCHSTbL(j,k,2)*COCHSTbR(i,l,1)
     .                           +COCHSTbL(j,k,1)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbR(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbR(j,l,2)))
     .                  *dsqrt(MCH2(i)*MCH2(j))/(1d0+epst3)**2
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
      enddo                     !3rd family/1st-2nd interference
       aux=aux
     . -(Ybq**2*((COCHSTbL(j,k,1)*U(i,2,1)-COCHSTbL(j,k,2)*U(i,2,2))
     .    *(COCHSTbL(i,k,1)*U(j,2,1)-COCHSTbL(i,k,2)*U(j,2,2))
     .    +(COCHSTbL(j,k,2)*U(i,2,1)+COCHSTbL(j,k,1)*U(i,2,2))
     .    *(COCHSTbL(i,k,2)*U(j,2,1)+COCHSTbL(i,k,1)*U(j,2,2)))
     .    +((COCHSUdL(j,1,1)*COCHSTbR(i,k,1)
     .                           -COCHSUdL(j,1,2)*COCHSTbR(i,k,2))
     .    *(COCHSUdL(i,1,1)*COCHSTbR(j,k,1)
     .                           -COCHSUdL(i,1,2)*COCHSTbR(j,k,2))
     .    +(COCHSUdL(j,1,2)*COCHSTbR(i,k,1)
     .                           +COCHSUdL(j,1,1)*COCHSTbR(i,k,2))
     .    *(COCHSUdL(i,1,2)*COCHSTbR(j,k,1)
     .                           +COCHSUdL(i,1,1)*COCHSTbR(j,k,2))))
     .                  *dsqrt(MCH2(i)*MCH2(j))/(1d0+epst3)**2
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
       auxe=auxe
     . -(Ybq**2*((COCHSTbL(j,k,1)*U(i,2,1)-COCHSTbL(j,k,2)*U(i,2,2))
     .    *(COCHSTbL(i,k,2)*U(j,2,1)+COCHSTbL(i,k,1)*U(j,2,2))
     .    -(COCHSTbL(j,k,2)*U(i,2,1)+COCHSTbL(j,k,1)*U(i,2,2))
     .    *(COCHSTbL(i,k,1)*U(j,2,1)-COCHSTbL(i,k,2)*U(j,2,2)))
     .    +((COCHSUdL(j,1,1)*COCHSTbR(i,k,1)
     .                           -COCHSUdL(j,1,2)*COCHSTbR(i,k,2))
     .    *(COCHSUdL(i,1,2)*COCHSTbR(j,k,1)
     .                           +COCHSUdL(i,1,1)*COCHSTbR(j,k,2))
     .    -(COCHSUdL(j,1,2)*COCHSTbR(i,k,1)
     .                           +COCHSUdL(j,1,1)*COCHSTbR(i,k,2))
     .    *(COCHSUdL(i,1,1)*COCHSTbR(j,k,1)
     .                           -COCHSUdL(i,1,2)*COCHSTbR(j,k,2))))
     .                  *dsqrt(MCH2(i)*MCH2(j))/(1d0+epst3)**2
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
      enddo                     !1st-2nd family contribution
       aux=aux+((COCHSUdL(j,1,1)*U(i,2,1)-COCHSUdL(j,1,2)*U(i,2,2))
     .    *(COCHSUdL(i,1,1)*U(j,2,1)-COCHSUdL(i,1,2)*U(j,2,2))
     .    +(COCHSUdL(j,1,2)*U(i,2,1)+COCHSUdL(j,1,1)*U(i,2,2))
     .    *(COCHSUdL(i,1,2)*U(j,2,1)+COCHSUdL(i,1,1)*U(j,2,2)))
     .      *Ybq**2/(1d0+epst3)**2*dsqrt(MCH2(i)*MCH2(j))
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MSU2(1)),dsqrt(MSU2(1)))
       auxe=auxe+((COCHSUdL(j,1,1)*U(i,2,1)-COCHSUdL(j,1,2)*U(i,2,2))
     .    *(COCHSUdL(i,1,2)*U(j,2,1)+COCHSUdL(i,1,1)*U(j,2,2))
     .    -(COCHSUdL(j,1,2)*U(i,2,1)+COCHSUdL(j,1,1)*U(i,2,2))
     .    *(COCHSUdL(i,1,1)*U(j,2,1)-COCHSUdL(i,1,2)*U(j,2,2)))
     .      *Ybq**2/(1d0+epst3)**2*dsqrt(MCH2(i)*MCH2(j))
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MSU2(1)),dsqrt(MSU2(1)))
      enddo
      enddo
      C1LRCHAR=-aux/2d0/(GF*MW)**2*MS0/MB0
      C1LRCHARI=-auxe/2d0/(GF*MW)**2*MS0/MB0
      
      aux=0d0                   !3rd family contribution
      auxe=0d0
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+((COCHSTbR(j,k,1)*COCHSTbL(i,l,1)
     .                           -COCHSTbR(j,k,2)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbR(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbR(j,l,2))
     .    +(COCHSTbR(j,k,2)*COCHSTbL(i,l,1)
     .                           +COCHSTbR(j,k,1)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbR(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbR(j,l,2)))
     .     /(1d0+epst3)**2
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
       auxe=auxe+((COCHSTbR(j,k,1)*COCHSTbL(i,l,1)
     .                           -COCHSTbR(j,k,2)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbR(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbR(j,l,2))
     .    -(COCHSTbR(j,k,2)*COCHSTbL(i,l,1)
     .                           +COCHSTbR(j,k,1)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbR(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbR(j,l,2)))
     .     /(1d0+epst3)**2
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
      enddo                     !3rd family/1st-2nd interference
       aux=aux-2d0*((COCHSTbR(j,k,1)*COCHSUdL(i,1,1)
     .                           -COCHSTbR(j,k,2)*COCHSUdL(i,1,2))
     .    *(COCHSTbL(i,k,1)*U(j,2,1)-COCHSTbL(i,k,2)*U(j,2,2))
     .    +(COCHSTbR(j,k,2)*COCHSUdL(i,1,1)
     .                           +COCHSTbR(j,k,1)*COCHSUdL(i,1,2))
     .    *(COCHSTbL(i,k,2)*U(j,2,1)+COCHSTbL(i,k,1)*U(j,2,2)))
     .     *Ybq/(1d0+epst3)**2
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
       auxe=auxe-2d0*((COCHSTbR(j,k,1)*COCHSUdL(i,1,1)
     .                           -COCHSTbR(j,k,2)*COCHSUdL(i,1,2))
     .    *(COCHSTbL(i,k,2)*U(j,2,1)+COCHSTbL(i,k,1)*U(j,2,2))
     .    -(COCHSTbR(j,k,2)*COCHSUdL(i,1,1)
     .                           +COCHSTbR(j,k,1)*COCHSUdL(i,1,2))
     .    *(COCHSTbL(i,k,1)*U(j,2,1)-COCHSTbL(i,k,2)*U(j,2,2)))
     .     *Ybq/(1d0+epst3)**2
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
      enddo                     !1st-2nd family contribution
       aux=aux+((U(j,2,1)*COCHSUdL(i,1,1)-U(j,2,2)*COCHSUdL(i,1,2))
     .    *(COCHSUdL(i,1,1)*U(j,2,1)-COCHSUdL(i,1,2)*U(j,2,2))
     .    +(U(j,2,2)*COCHSUdL(i,1,1)+U(j,2,1)*COCHSUdL(i,1,2))
     .    *(COCHSUdL(i,1,2)*U(j,2,1)+COCHSUdL(i,1,1)*U(j,2,2)))
     .     *Ybq**2/(1d0+epst3)**2
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MSU2(1)),dsqrt(MSU2(1)))
       auxe=auxe+((U(j,2,1)*COCHSUdL(i,1,1)-U(j,2,2)*COCHSUdL(i,1,2))
     .    *(COCHSUdL(i,1,2)*U(j,2,1)+COCHSUdL(i,1,1)*U(j,2,2))
     .    -(U(j,2,2)*COCHSUdL(i,1,1)+U(j,2,1)*COCHSUdL(i,1,2))
     .    *(COCHSUdL(i,1,1)*U(j,2,1)-COCHSUdL(i,1,2)*U(j,2,2)))
     .     *Ybq**2/(1d0+epst3)**2
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MSU2(1)),dsqrt(MSU2(1)))
      enddo
      enddo
      C2LRCHAR=-aux/2d0/(GF*MW)**2*MS0/MB0
      C2LRCHARI=-auxe/2d0/(GF*MW)**2*MS0/MB0

*  Running to sc0 [24] App.C
      etaS=asf(scR)/asc0

      CVLLCHAR=etaS**(6d0/21d0)*(1d0+asc0/4d0/Pi*1.3707d0*(1d0-etaS))
     .                 *CVLLCHAR
      CVLLCHARI=etaS**(6d0/21d0)*(1d0+asc0/4d0/Pi*1.3707d0*(1d0-etaS))
     .                 *CVLLCHARI

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

*          - Double Penguin contributions [1] Eqs.(6.12)-(6.22)
      aux=0d0
      do i=1,5
      aux=aux+sgn(MH0(i)-M_Bs**2)/
     .     dsqrt((MH0(i)-M_Bs**2)**2+MH0(i)*WIDTH(i)**2)
     . *((XH(i,1)-XH(i,2)*tanb)**2+XH(i,4)**2*(cosb+sinb*tanb)**2)
      enddo

      C2LRDPH=-(4d0*pi/(GF*MW))**2
     .          *(sigRLbs*sigLRbs-sigRLbsI*sigLRbsI)*aux
      C2LRDPHI=-(4d0*pi/(GF*MW))**2
     .          *(sigRLbs*sigLRbsI+sigRLbsI*sigLRbs)*aux

      aux=0d0
      auxe=0d0
      do i=1,5
      aux=aux+sgn(MH0(i)-M_Bs**2)/
     . dsqrt((MH0(i)-M_Bs**2)**2+MH0(i)*WIDTH(i)**2)
     . *((XH(i,1)-XH(i,2)*tanb)**2-XH(i,4)**2*(cosb+sinb*tanb)**2)
      auxe=auxe+sgn(MH0(i)-M_Bs**2)/
     . dsqrt((MH0(i)-M_Bs**2)**2+MH0(i)*WIDTH(i)**2)
     . *2d0*(XH(i,1)-XH(i,2)*tanb)*XH(i,4)*(cosb+sinb*tanb)
      enddo

      C1SLLDPH=-(4d0*pi/(GF*MW))**2/2d0
     . *((sigRLbs**2-sigRLbsI**2)*aux+2d0*auxe*sigRLbsI*sigRLbs)
      C1SLLDPHI=-(4d0*pi/(GF*MW))**2/2d0
     . *(2d0*aux*sigRLbsI*sigRLbs-(sigRLbs**2-sigRLbsI**2)*auxe)

*          - Summary
      CVLL=CVLLSM 

      I=1                              ! 0: SM; 1: NMSSM

      IF(I.eq.1)then
      CVLL=CVLL+CVLLHIG+CVLLCHAR
      CVLLI=CVLLHIGI+CVLLCHARI

      C1SLL=C1SLLHIG+C1SLLCHAR+C1SLLDPH
      C2SLL=C2SLLHIG+C2SLLCHAR
      C1SLLI=C1SLLHIGI+C1SLLCHARI+C1SLLDPHI
      C2SLLI=C2SLLHIGI+C2SLLCHARI

      C1LR=C1LRHIG+C1LRCHAR
      C2LR=C2LRHIG+C2LRCHAR+C2LRDPH
      C1LRI=C1LRHIGI+C1LRCHARI
      C2LRI=C2LRHIGI+C2LRCHARI+C2LRDPHI
      ENDIF

*	 2) Results for DMs
      eta=asc0/asmb

*          - `Bag' parameters from lattice [25]
      BVLL=BBs*0.551d0/0.985d0/eta**(6d0/23d0)      ! 'SM' Bag parameter from [5']
     .          /(1d0+asmb/4d0/Pi*1.6273d0*(1d0-eta))
      B1SLL=0.80d0
      B2SLL=0.71d0
      B1LR=1.75d0
      B2LR=1.16d0

      aux=1.42d0 !(M_Bs/(runmb(scb)+runmass(0.095d0,scb)))**2
      B1SLL=aux*B1SLL
      B2SLL=aux*B2SLL
      B1LR=aux*B1LR
      B2LR=aux*B2LR

*          - Running between sct and sc0 [24] Eqs.(3.1)-(3.19)

      PVLL=eta**(6d0/23d0)*(1d0+asmb/4d0/Pi*1.6273d0*(1d0-eta))*BVLL

      P1SLL=-5d0/8d0*(1.0153d0*eta**(-0.6315d0)-0.0153d0*eta**(0.7184d0)
     . +asmb/4d0/pi*(eta**(-0.6315d0)*(4.8177d0-5.2272d0*eta)
     .              +eta**(0.7184d0)*(0.3371d0+0.0724d0*eta)))*B1SLL
     . -3d0/2d0*(0.0081d0*(eta**(0.7184d0)-eta**(-0.6315d0))
     . +asmb/4d0/pi*(eta**(-0.6315d0)*(0.0531d0+0.0415d0*eta)
     .              -eta**(0.7184d0)*(0.0566d0+0.0380d0*eta)))*B2SLL

      P2SLL=-5d0/8d0*(1.9325d0*(eta**(-0.6315d0)-eta**(0.7184d0))
     . +asmb/4d0/pi*(eta**(-0.6315d0)*(9.1696d0-38.8778d0*eta)
     .            +eta**(0.7184d0)*(42.5021d0-12.7939d0*eta)))*B1SLL
     . -3d0/2d0*(1.0153d0*eta**(0.7184d0)-0.0153d0*eta**(-0.6315d0)
     . +asmb/4d0/pi*(eta**(-0.6315d0)*(0.1011d0+0.3083d0*eta)
     .      +eta**(0.7184d0)*(-7.1314d0+6.7219d0*eta)))*B2SLL

      P1LR=-(eta**(3d0/23d0)+asmb/4d0/Pi*(0.9250d0*eta**(-24d0/23d0)
     . +eta**(3d0/23d0)*(-2.0994d0+1.1744d0*eta)))*B1LR/2d0
     . +3d0/4d0*(2d0/3d0*(eta**(3d0/23d0)-eta**(-24d0/23d0))
     . +asmb/4d0/Pi*(eta**(3d0/23d0)*(-11.7329d0+0.7829*eta)
     .           +eta**(-24d0/23d0)*(-5.3048d0+16.2548d0*eta)))*B2LR

      P2LR=-(asmb/4d0/Pi*1.3875d0*(eta**(26d0/23d0)-eta**(-24d0/23d0)))
     .                                                      *B1LR/2d0
     . +3d0/4d0*(eta**(-24d0/23d0)+asmb/4d0/Pi*(eta**(-24d0/23d0)
     .      *(7.9572d0-8.8822d0*eta)+0.9250d0*eta**(26d0/23d0)))*B2LR

*          - Bs Mixing parameter DMs
      aux=PVLL*CVLL+P1SLL*C1SLL+P2SLL*C2SLL+P1LR*C1LR+P2LR*C2LR
      auxe=PVLL*CVLLI+P1SLL*C1SLLI+P2SLL*C2SLLI+P1LR*C1LRI+P2LR*C2LRI
      
      DMs=GF**2*MW**2/(24d0*pi**2)*M_Bs*fBs**2*VtsVtb2/(6.58211915d-13)
     .                  *dsqrt(aux**2+auxe**2)

*          - Error estimate
*      First: 2 sigma error bars from lattice Bag parameters: 
*      (2sigma, added quadratically)

      DMsMax=(dBBs*0.551d0*CVLL)**2

      aux1=-5d0/8d0*(1.0153d0*eta**(-0.6315d0)
     .                                   -0.0153d0*eta**(0.7184d0)
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(4.8177d0-5.2272d0*eta)
     .              +eta**(0.7184d0)*(0.3371d0+0.0724d0*eta)))*C1SLL
     . -5d0/8d0*(1.9325d0*(eta**(-0.6315d0)-eta**(0.7184d0))
     .      +asmb/4d0/pi*(eta**(-0.6315d0)*(9.1696d0-38.8778d0*eta)
     .            +eta**(0.7184d0)*(42.5021d0-12.7939d0*eta)))*C2SLL

      aux2=-5d0/8d0*(1.0153d0*eta**(-0.6315d0)
     .                                   -0.0153d0*eta**(0.7184d0)
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(4.8177d0-5.2272d0*eta)
     .              +eta**(0.7184d0)*(0.3371d0+0.0724d0*eta)))*C1SLLI
     . -5d0/8d0*(1.9325d0*(eta**(-0.6315d0)-eta**(0.7184d0))
     .      +asmb/4d0/pi*(eta**(-0.6315d0)*(9.1696d0-38.8778d0*eta)
     .            +eta**(0.7184d0)*(42.5021d0-12.7939d0*eta)))*C2SLLI
     
      DMsMax=DMsMax+4d0*(aux1**2+aux2**2)*(0.041d0*1.42d0)**2

      aux1=-3d0/2d0*(0.0081d0*(eta**(0.7184d0)-eta**(-0.6315d0))
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(0.0531d0+0.0415d0*eta)
     .              -eta**(0.7184d0)*(0.0566d0+0.0380d0*eta)))*C1SLL
     . -3d0/2d0*(1.0153d0*eta**(0.7184d0)-0.0153d0*eta**(-0.6315d0)
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(0.1011d0+0.3083d0*eta)
     .             +eta**(0.7184d0)*(-7.1314d0+6.7219d0*eta)))*C2SLL

      aux2=-3d0/2d0*(0.0081d0*(eta**(0.7184d0)-eta**(-0.6315d0))
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(0.0531d0+0.0415d0*eta)
     .              -eta**(0.7184d0)*(0.0566d0+0.0380d0*eta)))*C1SLLI
     . -3d0/2d0*(1.0153d0*eta**(0.7184d0)-0.0153d0*eta**(-0.6315d0)
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(0.1011d0+0.3083d0*eta)
     .             +eta**(0.7184d0)*(-7.1314d0+6.7219d0*eta)))*C2SLLI

      DMsMax=DMsMax+4d0*(aux1**2+aux2**2)*(0.095d0*1.42d0)**2

      aux1=-(eta**(3d0/23d0)+asmb/4d0/Pi*(0.9250d0*eta**(-24d0/23d0)
     .          +eta**(3d0/23d0)*(-2.0994d0+1.1744d0*eta)))/2d0*C1LR
     . -(asmb/4d0/Pi*1.3875d0*(eta**(26d0/23d0)-eta**(-24d0/23d0)))
     .                                                     *C2LR/2d0

      aux2=-(eta**(3d0/23d0)+asmb/4d0/Pi*(0.9250d0*eta**(-24d0/23d0)
     .          +eta**(3d0/23d0)*(-2.0994d0+1.1744d0*eta)))/2d0*C1LRI
     . -(asmb/4d0/Pi*1.3875d0*(eta**(26d0/23d0)-eta**(-24d0/23d0)))
     .                                                     *C2LRI/2d0

      IF(aux*aux1.gt.0d0)then
       DMsMax=DMsMax+4d0*(aux1*0.212d0*1.42d0)**2
       DMsmin=DMsMax+4d0*(aux1*0.067d0*1.42d0)**2
      ELSE
       DMsMax=DMsMax+4d0*(aux1*0.067d0*1.42d0)**2
       DMsmin=DMsMax+4d0*(aux1*0.212d0*1.42d0)**2
      ENDIF
      IF(auxe*aux2.gt.0d0)then
       DMsMax=DMsMax+4d0*(aux2*0.212d0*1.42d0)**2
       DMsmin=DMsmin+4d0*(aux2*0.067d0*1.42d0)**2
      ELSE
       DMsMax=DMsMax+4d0*(aux2*0.067d0*1.42d0)**2
       DMsmin=DMsmin+4d0*(aux2*0.212d0*1.42d0)**2
      ENDIF
      
      aux1=3d0/4d0*(2d0/3d0*(eta**(3d0/23d0)-eta**(-24d0/23d0))
     .      +asmb/4d0/Pi*(eta**(3d0/23d0)*(-11.7329d0+0.7829*eta)
     .           +eta**(-24d0/23d0)*(-5.3048d0+16.2548d0*eta)))*C1LR
     . +3d0/4d0*(eta**(-24d0/23d0)+asmb/4d0/Pi*(eta**(-24d0/23d0)
     .     *(7.9572d0-8.8822d0*eta)+0.9250d0*eta**(26d0/23d0)))*C2LR

      aux2=3d0/4d0*(2d0/3d0*(eta**(3d0/23d0)-eta**(-24d0/23d0))
     .      +asmb/4d0/Pi*(eta**(3d0/23d0)*(-11.7329d0+0.7829*eta)
     .           +eta**(-24d0/23d0)*(-5.3048d0+16.2548d0*eta)))*C1LRI
     . +3d0/4d0*(eta**(-24d0/23d0)+asmb/4d0/Pi*(eta**(-24d0/23d0)
     .     *(7.9572d0-8.8822d0*eta)+0.9250d0*eta**(26d0/23d0)))*C2LRI

      IF(aux*aux1.gt.0d0)then
       DMsMax=DMsMax+4d0*(aux1*0.054d0*1.42d0)**2
       DMsmin=DMsmin+4d0*(aux1*0.073d0*1.42d0)**2
      ELSE
       DMsMax=DMsMax+4d0*(aux1*0.073d0*1.42d0)**2
       DMsmin=DMsmin+4d0*(aux1*0.054d0*1.42d0)**2
      ENDIF
      IF(aux2*auxe.gt.0d0)then
       DMsMax=DMsMax+4d0*(aux2*0.054d0*1.42d0)**2
       DMsmin=DMsmin+4d0*(aux2*0.073d0*1.42d0)**2
      ELSE
       DMsMax=DMsMax+4d0*(aux2*0.073d0*1.42d0)**2
       DMsmin=DMsmin+4d0*(aux2*0.054d0*1.42d0)**2
      ENDIF

      DMsMax=(GF**2*MW**2/(24d0*pi**2)*M_Bs*fBs**2*VtsVtb2
     .            /6.58211915d-13)**2*DMsMax
      DMsmin=(GF**2*MW**2/(24d0*pi**2)*M_Bs*fBs**2*VtsVtb2
     .            /6.58211915d-13)**2*DMsmin

*      Second: Uncertainty from matching: 1% SM + 30% on each BSM contribution:
*      NB: SM error usually neglected in the literature
      aux=0.01d0*dabs(PVLL*CVLLSM)                               ! 1% SM

      IF(I.eq.1)then                                             ! 30% BSM
       aux=aux+0.3d0*(dsqrt((PVLL*CVLLHIG+P1SLL*C1SLLHIG+P2SLL*C2SLLHIG
     . +P1LR*C1LRHIG+P2LR*C2LRHIG)**2+(PVLL*CVLLHIGI+P1SLL*C1SLLHIGI
     . +P2SLL*C2SLLHIGI+P1LR*C1LRHIGI+P2LR*C2LRHIGI)**2)
     . +dsqrt((PVLL*CVLLCHAR+P1SLL*C1SLLCHAR
     .        +P2SLL*C2SLLCHAR+P1LR*C1LRCHAR+P2LR*C2LRCHAR)**2
     .   +(PVLL*CVLLCHARI+P1SLL*C1SLLCHARI
     .        +P2SLL*C2SLLCHARI+P1LR*C1LRCHARI+P2LR*C2LRCHARI)**2)
     . +dsqrt((P1SLL*C1SLLDPH+P2LR*C2LRDPH)**2
     .       +(P1SLL*C1SLLDPHI+P2LR*C2LRDPHI)**2))
      ENDIF

      DMsmax=DMs+dsqrt(DMsMax)+GF**2*MW**2/(24d0*pi**2)*M_Bs*fBs**2
     .                            *VtsVtb2/(6.58211915d-13)*aux
      DMsmin=DMs-dsqrt(DMsmin)-GF**2*MW**2/(24d0*pi**2)*M_Bs*fBs**2
     .                            *VtsVtb2/(6.58211915d-13)*aux
      DMsmin=Max(0d0,DMsmin)

*      Third: CKM and lattice form factor
      aux=dsqrt((dVtsVtb2/VtsVtb2)**2+(2d0*dfBs)**2)

*      Total error bars
      DMsmax=DMsMax*(1d0+aux)
      DMsmin=DMsmin*(1d0-aux)

!      print*,'DMs',DMsmin,DMs,DMsmax
*      Comparison with experimental data (source [1',3'])
*       (Recall: 2 sigma bounds: 17.715ps-1 < DMs=17.757ps-1 < 17.799ps-1)
      prob(33)=0d0

      IF(DMsmin.GE.DMsexpMax)
     .     PROB(33)=DMsmin/DMsexpMax-1d0
      IF(DMsmax.LE.DMsexpmin)
     .     PROB(33)=DMsmax/DMsexpMin-1d0

!      csqb=csqb+4.d0*(DMs-(DMsexpmin+DMsexpMax)/2d0)**2
!     c /((DMsMax-DMsmin)**2+(DMsexpMax-DMsexpmin)**2)

*	 3) Wilson coefficients for DMd at the matching scale
      
*          - Charged Higgs Boxes [1] App.(A.4.i)
*  Coefficients at MHC
      MCHH(1)=MW
      MCHH(2)=dsqrt(MHC)

      aux=mth/dsqrt(vuq**2+vdq**2)
      PRLH_tb(1)=-aux
      PRLHI_tb(1)=0d0
      PRLH_tb(2)=aux/tanb*(1d0-epstb*(1d0/tanb+tanb))
      PRLHI_tb(2)=-aux*epstbI*(1d0/tanb**2+1d0)
      PRLH_ts(1)=-aux
      PRLHI_ts(1)=0d0
      PRLH_ts(2)=aux/tanb*(1d0-(1d0/tanb+tanb)
     . *(epsts-((epsY31*(1d0+epst0)+epsY31I*epst0I)*(epstb-epsts)
     .         -(epsY31I*(1d0+epst0)-epsY31*epst0I)*(epstbI-epstsI))
     .                                    /((1d0+epst0)**2+epst0I**2)))
      PRLHI_ts(2)=-aux*(1d0/tanb**2+1d0)
     . *(epstsI-((epsY31I*(1d0+epst0)-epsY31*epst0I)*(epstb-epsts)
     .         +(epsY31*(1d0+epst0)+epsY31I*epst0I)*(epstbI-epstsI))
     .                                    /((1d0+epst0)**2+epst0I**2))

      aux=mbh/dsqrt(vuq**2+vdq**2)
      PLRH_tb(1)=aux
      PLRHI_tb(1)=0d0
      PLRH_tb(2)=aux*(tanb-(epst3+epst3I**2/(1d0+epst3))/(1d0+epst3)
     .                                        *(tanb+1d0/tanb))
      PLRHI_tb(2)=aux*(tanb+1d0/tanb)*epst3I/(1d0+epst3)**2
      aux=runmass(0.005d0,dsqrt(MHC))/dsqrt(vuq**2+vdq**2)
      PLRH_ts(1)=aux
      PLRHI_ts(1)=0d0
      PLRH_ts(2)=aux*(tanb-(tanb+1d0/tanb)/(1d0+epst1)
     .   *(epst1+epst1I**2/(1d0+epst1)
     . +((1d0+epst0)*(epsY13+epsY13I
     .                       *(epst3+epst1*(1d0+epst3)/(1d0+epst1))
     .               +epsY31*(epst1-epst3)-epsY31I
     .         *(epst3I*(1d0+epst1)-epst1I*(1d0+epst3)/(1d0+epst1)))
     .  -epst0I*(epsY13I-epsY13*(epst3I+epst1I*(1d0+epst3)/(1d0+epst1))
     .         -epsY31I*(epst1-epst3)-epsY31*(epst3I*(1d0+epst1)
     .           -epst1I*(1d0+epst3)/(1d0+epst1))))
     .               /((1d0+epst0)**2+epst0I**1)))
      PLRHI_ts(2)=aux*((tanb+1d0/tanb)/(1d0+epst1)
     .   *(epst1I/(1d0+epst1)
     . -(epst0I*(epsY13+epsY13I*(epst3+epst1*(1d0+epst3)/(1d0+epst1))
     .               +epsY31*(epst1-epst3)-epsY31I
     .         *(epst3I*(1d0+epst1)-epst1I*(1d0+epst3)/(1d0+epst1)))
     .  -(1d0+epst0)*(epsY13I-epsY13*(epst3I
     .                              +epst1I*(1d0+epst3)/(1d0+epst1))
     .         -epsY31I*(epst1-epst3)-epsY31*(epst3I*(1d0+epst1)
     .           -epst1I*(1d0+epst3)/(1d0+epst1))))
     .               /((1d0+epst0)**2+epst0I**2)))

      CVLLHIG=-g2q/2d0*(PRLH_tb(2)*PRLH_ts(2)+PRLHI_tb(2)*PRLHI_ts(2))
     .                          *mth**2*D0B(MW,MCHH(2),mth,mth)
     . +((PRLH_tb(2)*PRLH_tb(2)-PRLHI_tb(2)*PRLHI_tb(2))
     .   *(PRLH_ts(2)*PRLH_ts(2)-PRLHI_ts(2)*PRLHI_ts(2))
     .  +(PRLH_tb(2)*PRLHI_tb(2)+PRLHI_tb(2)*PRLH_tb(2))
     .   *(PRLH_ts(2)*PRLHI_ts(2)+PRLHI_ts(2)*PRLH_ts(2)))
     .              *D2B(MCHH(2),MCHH(2),mth,mth)/8d0
     . +((PRLH_tb(2)*PRLH_tb(1)-PRLHI_tb(2)*PRLHI_tb(1))
     .   *(PRLH_ts(2)*PRLH_ts(1)-PRLHI_ts(2)*PRLHI_ts(1))
     .  +(PRLH_tb(2)*PRLHI_tb(1)+PRLHI_tb(2)*PRLH_tb(1))
     .   *(PRLH_ts(2)*PRLHI_ts(1)+PRLHI_ts(2)*PRLH_ts(1)))
     .              *D2B(MCHH(1),MCHH(2),mth,mth)/4d0
      CVLLHIG=CVLLHIG/(GF*MW)**2

      CVLLHIGI=-g2q/2d0*(PRLH_tb(2)*PRLHI_ts(2)-PRLHI_tb(2)*PRLH_ts(2))
     .                          *mth**2*D0B(MW,MCHH(2),mth,mth)
     . +((PRLH_tb(2)*PRLH_tb(2)-PRLHI_tb(2)*PRLHI_tb(2))
     .   *(PRLH_ts(2)*PRLHI_ts(2)+PRLHI_ts(2)*PRLH_ts(2))
     .  -(PRLH_tb(2)*PRLHI_tb(2)+PRLHI_tb(2)*PRLH_tb(2))
     .   *(PRLH_ts(2)*PRLH_ts(2)-PRLHI_ts(2)*PRLHI_ts(2)))
     .              *D2B(MCHH(2),MCHH(2),mth,mth)/8d0
     . +((PRLH_tb(2)*PRLH_tb(1)-PRLHI_tb(2)*PRLHI_tb(1))
     .   *(PRLH_ts(2)*PRLHI_ts(1)+PRLHI_ts(2)*PRLH_ts(1))
     .  -(PRLH_tb(2)*PRLHI_tb(1)+PRLHI_tb(2)*PRLH_tb(1))
     .   *(PRLH_ts(2)*PRLH_ts(1)-PRLHI_ts(2)*PRLHI_ts(1)))
     .              *D2B(MCHH(1),MCHH(2),mth,mth)/4d0
      CVLLHIGI=CVLLHIGI/(GF*MW)**2

      C1SLLHIG=0d0
      DO I=1,2
      DO J=1,2
       C1SLLHIG=C1SLLHIG
     .  +((PLRH_tb(J)*PLRH_tb(I)-PLRHI_tb(J)*PLRHI_tb(I))
     .   *(PRLH_ts(I)*PRLH_ts(J)-PRLHI_ts(I)*PRLHI_ts(J))
     .  +(PLRH_tb(J)*PLRHI_tb(I)+PLRHI_tb(J)*PLRH_tb(I))
     .   *(PRLH_ts(I)*PRLHI_ts(J)+PRLHI_ts(I)*PRLH_ts(J)))
     .              *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)/2d0
      ENDDO
      ENDDO
      C1SLLHIG=C1SLLHIG/(GF*MW)**2

      C1SLLHIGI=0d0
      DO I=1,2
      DO J=1,2
       C1SLLHIGI=C1SLLHIGI
     .  +((PLRH_tb(J)*PLRH_tb(I)-PLRHI_tb(J)*PLRHI_tb(I))
     .   *(PRLH_ts(I)*PRLHI_ts(J)+PRLHI_ts(I)*PRLH_ts(J))
     .  -(PLRH_tb(J)*PLRHI_tb(I)+PLRHI_tb(J)*PLRH_tb(I))
     .   *(PRLH_ts(I)*PRLH_ts(J)-PRLHI_ts(I)*PRLHI_ts(J)))
     .              *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)/2d0
      ENDDO
      ENDDO
      C1SLLHIGI=C1SLLHIGI/(GF*MW)**2

      C1LRHIG=0d0
      DO I=1,2
      DO J=1,2
       C1LRHIG=C1LRHIG
     .  +((PRLH_tb(J)*PLRH_tb(I)-PRLHI_tb(J)*PLRHI_tb(I))
     .   *(PRLH_ts(I)*PLRH_ts(J)-PRLHI_ts(I)*PLRHI_ts(J))
     .  +(PRLH_tb(J)*PLRHI_tb(I)+PRLHI_tb(J)*PLRH_tb(I))
     .   *(PRLH_ts(I)*PLRHI_ts(J)+PRLHI_ts(I)*PLRH_ts(J)))
     .                          *D2B(MCHH(I),MCHH(J),mth,mth)/4d0
      ENDDO
      ENDDO
      C1LRHIG=C1LRHIG/(GF*MW)**2

      C1LRHIGI=0d0
      DO I=1,2
      DO J=1,2
       C1LRHIGI=C1LRHIGI
     .  +((PRLH_tb(J)*PLRH_tb(I)-PRLHI_tb(J)*PLRHI_tb(I))
     .   *(PRLH_ts(I)*PLRHI_ts(J)+PRLHI_ts(I)*PLRH_ts(J))
     .  -(PRLH_tb(J)*PLRHI_tb(I)+PRLHI_tb(J)*PLRH_tb(I))
     .   *(PRLH_ts(I)*PLRH_ts(J)-PRLHI_ts(I)*PLRHI_ts(J)))
     .                          *D2B(MCHH(I),MCHH(J),mth,mth)/4d0
      ENDDO
      ENDDO
      C1LRHIGI=C1LRHIGI/(GF*MW)**2

      aux=(1d0+epst3)/(1d0+epst1)
      C2LRHIG=0d0
      DO I=1,2
       C2LRHIG=C2LRHIG
     .  -g2q*(PLRH_tb(I)*PLRH_ts(I)+PLRHI_tb(I)*PLRHI_ts(I))/2d0
     .    *(D2B(MCHH(I),MW,mth,mth)
     .      -2d0*aux*D2B(MCHH(I),MW,mth,0d0)
     .      +aux**2*D2B(MCHH(I),MW,0d0,0d0))
      DO J=1,2
       C2LRHIG=C2LRHIG
     .  +((PLRH_tb(J)*PRLH_tb(I)-PLRHI_tb(J)*PRLHI_tb(I))
     .   *(PRLH_ts(I)*PLRH_ts(J)-PRLHI_ts(I)*PLRHI_ts(J))
     .  +(PLRH_tb(J)*PRLHI_tb(I)+PLRHI_tb(J)*PRLH_tb(I))
     .   *(PRLH_ts(I)*PLRHI_ts(J)+PRLHI_ts(I)*PLRH_ts(J)))
     .                     *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)
      ENDDO
      ENDDO
      C2LRHIG=C2LRHIG/(GF*MW)**2

      C2LRHIGI=0d0
      DO I=1,2
       C2LRHIGI=C2LRHIGI
     .  -g2q*(PLRH_tb(I)*PLRHI_ts(I)-PLRHI_tb(I)*PLRH_ts(I))/2d0
     .    *(D2B(MCHH(I),MW,mth,mth)
     .      -2d0*aux*D2B(MCHH(I),MW,mth,0d0)
     .      +aux**2*D2B(MCHH(I),MW,0d0,0d0))
      DO J=1,2
       C2LRHIGI=C2LRHIGI
     .  +((PLRH_tb(J)*PRLH_tb(I)-PLRHI_tb(J)*PRLHI_tb(I))
     .   *(PRLH_ts(I)*PLRHI_ts(J)+PRLHI_ts(I)*PLRH_ts(J))
     .  -(PLRH_tb(J)*PRLHI_tb(I)+PLRHI_tb(J)*PRLH_tb(I))
     .   *(PRLH_ts(I)*PLRH_ts(J)-PRLHI_ts(I)*PLRHI_ts(J)))
     .                     *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)
      ENDDO
      ENDDO
      C2LRHIGI=C2LRHIGI/(GF*MW)**2

*  Running to sc0 [24] App.C
      etaH=asf(dsqrt(MHC))/asc0

      CVLLHIG=etaH**(6d0/21d0)*(1d0+asc0/4d0/Pi*1.3707d0*(1d0-etaH))
     .                 *CVLLHIG
      CVLLHIGI=etaH**(6d0/21d0)*(1d0+asc0/4d0/Pi*1.3707d0*(1d0-etaH))
     .                 *CVLLHIGI

      aux=C1SLLHIG
      C1SLLHIG=(1.0153d0*etaH**(-0.6916d0)-0.0153d0*etaH**(0.7869d0)
     . +asc0/4d0/Pi*(etaH**(-0.6916d0)*(5.6478d0-6.0350d0*etaH)
     .              +etaH**(0.7869d0)*(0.3272d0+0.0600*etaH)))*aux

      C2SLLHIG=(0.0081d0*(etaH**(0.7869d0)-etaH**(-0.6916d0))
     . +asc0/4d0/Pi*(etaH**(0.7869d0)*(-0.0618d0-0.0315d0*etaH)
     .              +etaH**(-0.6916d0)*(0.0454d0+0.0479d0*etaH)))*aux
      aux=C1SLLHIGI
      C1SLLHIGI=(1.0153d0*etaH**(-0.6916d0)-0.0153d0*etaH**(0.7869d0)
     . +asc0/4d0/Pi*(etaH**(-0.6916d0)*(5.6478d0-6.0350d0*etaH)
     .              +etaH**(0.7869d0)*(0.3272d0+0.0600*etaH)))*aux

      C2SLLHIGI=(0.0081d0*(etaH**(0.7869d0)-etaH**(-0.6916d0))
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

*          - Chargino / squark Boxes [1] App.(A.4.ii)
*  Coefficients at scR
      
      aux=0d0                   !3rd family contribution
      auxe=0d0
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+((COCHSTbL(j,k,1)*COCHSTbL(i,l,1)
     .                           -COCHSTbL(j,k,2)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbL(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbL(j,l,2))
     .    +(COCHSTbL(j,k,2)*COCHSTbL(i,l,1)
     .                           +COCHSTbL(j,k,1)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbL(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbL(j,l,2)))
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
       auxe=auxe+((COCHSTbL(j,k,1)*COCHSTbL(i,l,1)
     .                           -COCHSTbL(j,k,2)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbL(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbL(j,l,2))
     .    -(COCHSTbL(j,k,2)*COCHSTbL(i,l,1)
     .                           +COCHSTbL(j,k,1)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbL(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbL(j,l,2)))
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
     .      
      enddo                     !3rd family/1st-2nd interference
       aux=aux-2d0*g2q*((V(j,1,1)*V(i,1,1)+V(j,1,2)*V(i,1,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbL(j,k,1)
     .                           +COCHSTbL(i,k,2)*COCHSTbL(j,k,2))
     .     +(V(j,1,2)*V(i,1,1)-V(j,1,1)*V(i,1,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbL(j,k,2)
     .                           -COCHSTbL(i,k,2)*COCHSTbL(j,k,1)))
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
       auxe=auxe-2d0*g2q*((V(j,1,1)*V(i,1,1)+V(j,1,2)*V(i,1,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbL(j,k,1)
     .                           -COCHSTbL(i,k,1)*COCHSTbL(j,k,2))
     .     +(V(j,1,2)*V(i,1,1)-V(j,1,1)*V(i,1,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbL(j,k,1)
     .                           +COCHSTbL(i,k,2)*COCHSTbL(j,k,2)))
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
      enddo                     !1st-2nd family contribution
       aux=aux+g2q**2*((V(j,1,1)*V(i,1,1)+V(j,1,2)*V(i,1,2))**2
     .                +(V(j,1,2)*V(i,1,1)-V(j,1,1)*V(i,1,2))**2)
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MSU2(1)),dsqrt(MSU2(1)))
      enddo
      enddo
      CVLLCHAR=aux/8d0/(GF*MW)**2
      CVLLCHARI=auxe/8d0/(GF*MW)**2

      aux=0d0                   !3rd family contribution
      auxe=0d0
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+((COCHSTbR(j,k,1)*COCHSTbR(i,l,1)
     .                           -COCHSTbR(j,k,2)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbL(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbL(j,l,2))
     .    +(COCHSTbR(j,k,2)*COCHSTbR(i,l,1)
     .                           +COCHSTbR(j,k,1)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbL(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbL(j,l,2)))
     .                  *dsqrt(MCH2(i)*MCH2(j))/(1d0+epst3)**2
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
       auxe=auxe+((COCHSTbR(j,k,1)*COCHSTbR(i,l,1)
     .                           -COCHSTbR(j,k,2)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbL(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbL(j,l,2))
     .    -(COCHSTbR(j,k,2)*COCHSTbR(i,l,1)
     .                           +COCHSTbR(j,k,1)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbL(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbL(j,l,2)))
     .                  *dsqrt(MCH2(i)*MCH2(j))/(1d0+epst3)**2
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
      enddo                     !3rd family/1st-2nd interference
       aux=aux-2d0*((COCHSTbR(j,k,1)*U(i,2,1)-COCHSTbR(j,k,2)*U(i,2,2))
     .    *(COCHSTbL(i,k,1)*COCHSUdL(j,1,1)
     .                           -COCHSTbL(i,k,2)*COCHSUdL(j,1,2))
     .    +(COCHSTbR(j,k,2)*U(i,2,1)+COCHSTbR(j,k,1)*U(i,2,2))
     .    *(COCHSTbL(i,k,2)*COCHSUdL(j,1,1)
     .                           +COCHSTbL(i,k,1)*COCHSUdL(j,1,2)))
     .         *Ybq/(1+epst3)**2*dsqrt(MCH2(i)*MCH2(j))
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
       auxe=auxe
     .   -2d0*((COCHSTbR(j,k,1)*U(i,2,1)-COCHSTbR(j,k,2)*U(i,2,2))
     .    *(COCHSTbL(i,k,2)*COCHSUdL(j,1,1)
     .                           +COCHSTbL(i,k,1)*COCHSUdL(j,1,2))
     .    -(COCHSTbR(j,k,2)*U(i,2,1)+COCHSTbR(j,k,1)*U(i,2,2))
     .    *(COCHSTbL(i,k,1)*COCHSUdL(j,1,1)
     .                           -COCHSTbL(i,k,2)*COCHSUdL(j,1,2)))
     .         *Ybq/(1+epst3)**2*dsqrt(MCH2(i)*MCH2(j))
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
      enddo                     !1st-2nd family contribution
       aux=aux+((U(j,2,1)*U(i,2,1)-U(j,2,2)*U(i,2,2))
     .    *(COCHSUdL(i,1,1)*COCHSUdL(j,1,1)
     .                           -COCHSUdL(i,1,2)*COCHSUdL(j,1,2))
     .    +(U(j,2,2)*U(i,2,1)+U(j,2,1)*U(i,2,2))
     .    *(COCHSUdL(i,1,2)*COCHSUdL(j,1,1)
     .                           +COCHSUdL(i,1,1)*COCHSUdL(j,1,2)))
     .         *Ybq**2/(1+epst3)**2*dsqrt(MCH2(i)*MCH2(j))
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MSU2(1)),dsqrt(MSU2(1)))
       auxe=auxe+((U(j,2,1)*U(i,2,1)-U(j,2,2)*U(i,2,2))
     .    *(COCHSUdL(i,1,2)*COCHSUdL(j,1,1)
     .                           +COCHSUdL(i,1,1)*COCHSUdL(j,1,2))
     .    -(U(j,2,2)*U(i,2,1)+U(j,2,1)*U(i,2,2))
     .    *(COCHSUdL(i,1,1)*COCHSUdL(j,1,1)
     .                           -COCHSUdL(i,1,2)*COCHSUdL(j,1,2)))
     .         *Ybq**2/(1+epst3)**2*dsqrt(MCH2(i)*MCH2(j))
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MSU2(1)),dsqrt(MSU2(1)))
      enddo
      enddo
      C1SLLCHAR=-aux/4d0/(GF*MW)**2
      C2SLLCHAR=aux/16d0/(GF*MW)**2
      C1SLLCHARI=-auxe/4d0/(GF*MW)**2
      C2SLLCHARI=auxe/16d0/(GF*MW)**2

      aux=0d0                   !3rd family contribution
      auxe=0d0
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+((COCHSTbL(j,k,1)*COCHSTbR(i,l,1)
     .                           -COCHSTbL(j,k,2)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbR(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbR(j,l,2))
     .    +(COCHSTbL(j,k,2)*COCHSTbR(i,l,1)
     .                           +COCHSTbL(j,k,1)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbR(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbR(j,l,2)))
     .                  *dsqrt(MCH2(i)*MCH2(j))/(1d0+epst3)**2
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
       auxe=auxe+((COCHSTbL(j,k,1)*COCHSTbR(i,l,1)
     .                           -COCHSTbL(j,k,2)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbR(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbR(j,l,2))
     .    -(COCHSTbL(j,k,2)*COCHSTbR(i,l,1)
     .                           +COCHSTbL(j,k,1)*COCHSTbR(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbR(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbR(j,l,2)))
     .                  *dsqrt(MCH2(i)*MCH2(j))/(1d0+epst3)**2
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
      enddo                     !3rd family/1st-2nd interference
       aux=aux
     . -(Ybq**2*((COCHSTbL(j,k,1)*U(i,2,1)-COCHSTbL(j,k,2)*U(i,2,2))
     .    *(COCHSTbL(i,k,1)*U(j,2,1)-COCHSTbL(i,k,2)*U(j,2,2))
     .    +(COCHSTbL(j,k,2)*U(i,2,1)+COCHSTbL(j,k,1)*U(i,2,2))
     .    *(COCHSTbL(i,k,2)*U(j,2,1)+COCHSTbL(i,k,1)*U(j,2,2)))
     .    +((COCHSUdL(j,1,1)*COCHSTbR(i,k,1)
     .                           -COCHSUdL(j,1,2)*COCHSTbR(i,k,2))
     .    *(COCHSUdL(i,1,1)*COCHSTbR(j,k,1)
     .                           -COCHSUdL(i,1,2)*COCHSTbR(j,k,2))
     .    +(COCHSUdL(j,1,2)*COCHSTbR(i,k,1)
     .                           +COCHSUdL(j,1,1)*COCHSTbR(i,k,2))
     .    *(COCHSUdL(i,1,2)*COCHSTbR(j,k,1)
     .                           +COCHSUdL(i,1,1)*COCHSTbR(j,k,2))))
     .                  *dsqrt(MCH2(i)*MCH2(j))/(1d0+epst3)**2
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
       auxe=auxe
     . -(Ybq**2*((COCHSTbL(j,k,1)*U(i,2,1)-COCHSTbL(j,k,2)*U(i,2,2))
     .    *(COCHSTbL(i,k,2)*U(j,2,1)+COCHSTbL(i,k,1)*U(j,2,2))
     .    -(COCHSTbL(j,k,2)*U(i,2,1)+COCHSTbL(j,k,1)*U(i,2,2))
     .    *(COCHSTbL(i,k,1)*U(j,2,1)-COCHSTbL(i,k,2)*U(j,2,2)))
     .    +((COCHSUdL(j,1,1)*COCHSTbR(i,k,1)
     .                           -COCHSUdL(j,1,2)*COCHSTbR(i,k,2))
     .    *(COCHSUdL(i,1,2)*COCHSTbR(j,k,1)
     .                           +COCHSUdL(i,1,1)*COCHSTbR(j,k,2))
     .    -(COCHSUdL(j,1,2)*COCHSTbR(i,k,1)
     .                           +COCHSUdL(j,1,1)*COCHSTbR(i,k,2))
     .    *(COCHSUdL(i,1,1)*COCHSTbR(j,k,1)
     .                           -COCHSUdL(i,1,2)*COCHSTbR(j,k,2))))
     .                  *dsqrt(MCH2(i)*MCH2(j))/(1d0+epst3)**2
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
      enddo                     !1st-2nd family contribution
       aux=aux+((COCHSUdL(j,1,1)*U(i,2,1)-COCHSUdL(j,1,2)*U(i,2,2))
     .    *(COCHSUdL(i,1,1)*U(j,2,1)-COCHSUdL(i,1,2)*U(j,2,2))
     .    +(COCHSUdL(j,1,2)*U(i,2,1)+COCHSUdL(j,1,1)*U(i,2,2))
     .    *(COCHSUdL(i,1,2)*U(j,2,1)+COCHSUdL(i,1,1)*U(j,2,2)))
     .      *Ybq**2/(1d0+epst3)**2*dsqrt(MCH2(i)*MCH2(j))
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MSU2(1)),dsqrt(MSU2(1)))
       auxe=auxe+((COCHSUdL(j,1,1)*U(i,2,1)-COCHSUdL(j,1,2)*U(i,2,2))
     .    *(COCHSUdL(i,1,2)*U(j,2,1)+COCHSUdL(i,1,1)*U(j,2,2))
     .    -(COCHSUdL(j,1,2)*U(i,2,1)+COCHSUdL(j,1,1)*U(i,2,2))
     .    *(COCHSUdL(i,1,1)*U(j,2,1)-COCHSUdL(i,1,2)*U(j,2,2)))
     .      *Ybq**2/(1d0+epst3)**2*dsqrt(MCH2(i)*MCH2(j))
     . *D0B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MSU2(1)),dsqrt(MSU2(1)))
      enddo
      enddo
      C1LRCHAR=-aux/2d0/(GF*MW)**2*MD/MB0
      C1LRCHARI=-auxe/2d0/(GF*MW)**2*MD/MB0
      
      aux=0d0                   !3rd family contribution
      auxe=0d0
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+((COCHSTbR(j,k,1)*COCHSTbL(i,l,1)
     .                           -COCHSTbR(j,k,2)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbR(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbR(j,l,2))
     .    +(COCHSTbR(j,k,2)*COCHSTbL(i,l,1)
     .                           +COCHSTbR(j,k,1)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbR(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbR(j,l,2)))
     .     /(1d0+epst3)**2
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
       auxe=auxe+((COCHSTbR(j,k,1)*COCHSTbL(i,l,1)
     .                           -COCHSTbR(j,k,2)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,2)*COCHSTbR(j,l,1)
     .                           +COCHSTbL(i,k,1)*COCHSTbR(j,l,2))
     .    -(COCHSTbR(j,k,2)*COCHSTbL(i,l,1)
     .                           +COCHSTbR(j,k,1)*COCHSTbL(i,l,2))
     .    *(COCHSTbL(i,k,1)*COCHSTbR(j,l,1)
     .                           -COCHSTbL(i,k,2)*COCHSTbR(j,l,2)))
     .     /(1d0+epst3)**2
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MST2(l)))
      enddo                     !3rd family/1st-2nd interference
       aux=aux-2d0*((COCHSTbR(j,k,1)*COCHSUdL(i,1,1)
     .                           -COCHSTbR(j,k,2)*COCHSUdL(i,1,2))
     .    *(COCHSTbL(i,k,1)*U(j,2,1)-COCHSTbL(i,k,2)*U(j,2,2))
     .    +(COCHSTbR(j,k,2)*COCHSUdL(i,1,1)
     .                           +COCHSTbR(j,k,1)*COCHSUdL(i,1,2))
     .    *(COCHSTbL(i,k,2)*U(j,2,1)+COCHSTbL(i,k,1)*U(j,2,2)))
     .     *Ybq/(1d0+epst3)**2
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
       auxe=auxe-2d0*((COCHSTbR(j,k,1)*COCHSUdL(i,1,1)
     .                           -COCHSTbR(j,k,2)*COCHSUdL(i,1,2))
     .    *(COCHSTbL(i,k,2)*U(j,2,1)+COCHSTbL(i,k,1)*U(j,2,2))
     .    -(COCHSTbR(j,k,2)*COCHSUdL(i,1,1)
     .                           +COCHSTbR(j,k,1)*COCHSUdL(i,1,2))
     .    *(COCHSTbL(i,k,1)*U(j,2,1)-COCHSTbL(i,k,2)*U(j,2,2)))
     .     *Ybq/(1d0+epst3)**2
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MST2(k)),dsqrt(MSU2(1)))
      enddo                     !1st-2nd family contribution
       aux=aux+((U(j,2,1)*COCHSUdL(i,1,1)-U(j,2,2)*COCHSUdL(i,1,2))
     .    *(COCHSUdL(i,1,1)*U(j,2,1)-COCHSUdL(i,1,2)*U(j,2,2))
     .    +(U(j,2,2)*COCHSUdL(i,1,1)+U(j,2,1)*COCHSUdL(i,1,2))
     .    *(COCHSUdL(i,1,2)*U(j,2,1)+COCHSUdL(i,1,1)*U(j,2,2)))
     .     *Ybq**2/(1d0+epst3)**2
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MSU2(1)),dsqrt(MSU2(1)))
       auxe=auxe+((U(j,2,1)*COCHSUdL(i,1,1)-U(j,2,2)*COCHSUdL(i,1,2))
     .    *(COCHSUdL(i,1,2)*U(j,2,1)+COCHSUdL(i,1,1)*U(j,2,2))
     .    -(U(j,2,2)*COCHSUdL(i,1,1)+U(j,2,1)*COCHSUdL(i,1,2))
     .    *(COCHSUdL(i,1,1)*U(j,2,1)-COCHSUdL(i,1,2)*U(j,2,2)))
     .     *Ybq**2/(1d0+epst3)**2
     . *D2B(dsqrt(MCH2(i)),dsqrt(MCH2(j)),dsqrt(MSU2(1)),dsqrt(MSU2(1)))
      enddo
      enddo
      C2LRCHAR=-aux/2d0/(GF*MW)**2*MD/MB0
      C2LRCHARI=-auxe/2d0/(GF*MW)**2*MD/MB0

*  Running to sc0 [24] App.C
      etaS=asf(scR)/asc0

      CVLLCHAR=etaS**(6d0/21d0)*(1d0+asc0/4d0/Pi*1.3707d0*(1d0-etaS))
     .                 *CVLLCHAR
      CVLLCHARI=etaS**(6d0/21d0)*(1d0+asc0/4d0/Pi*1.3707d0*(1d0-etaS))
     .                 *CVLLCHARI

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

*          - Double Penguin contributions [1] Eqs.(6.12)-(6.22)
      aux=0d0
      do i=1,5
      aux=aux+sgn(MH0(i)-M_Bd**2)/
     .     dsqrt((MH0(i)-M_Bd**2)**2+MH0(i)*WIDTH(i)**2)
     . *((XH(i,1)-XH(i,2)*tanb)**2+XH(i,4)**2*(cosb+sinb*tanb)**2)
      enddo

      C2LRDPH=-(4d0*pi/(GF*MW))**2
     .          *(sigRLbd*sigLRbd-sigRLbdI*sigLRbdI)*aux
      C2LRDPHI=-(4d0*pi/(GF*MW))**2
     .          *(sigRLbd*sigLRbdI+sigRLbdI*sigLRbd)*aux

      aux=0d0
      auxe=0d0
      do i=1,5
      aux=aux+sgn(MH0(i)-M_Bd**2)/
     . dsqrt((MH0(i)-M_Bd**2)**2+MH0(i)*WIDTH(i)**2)
     . *((XH(i,1)-XH(i,2)*tanb)**2-XH(i,4)**2*(cosb+sinb*tanb)**2)
      auxe=auxe+sgn(MH0(i)-M_Bs**2)/
     . dsqrt((MH0(i)-M_Bd**2)**2+MH0(i)*WIDTH(i)**2)
     . *2d0*(XH(i,1)-XH(i,2)*tanb)*XH(i,4)*(cosb+sinb*tanb)
      enddo

      C1SLLDPH=-(4d0*pi/(GF*MW))**2/2d0
     . *((sigRLbd**2-sigRLbdI**2)*aux+2d0*auxe*sigRLbdI*sigRLbd)
      C1SLLDPHI=-(4d0*pi/(GF*MW))**2/2d0
     . *(2d0*aux*sigRLbdI*sigRLbd-(sigRLbd**2-sigRLbdI**2)*auxe)

*          - Summary
      CVLL=CVLLSM 

      I=1                              ! 0: SM; 1: NMSSM

      IF(I.eq.1)then
      CVLL=CVLL+CVLLHIG+CVLLCHAR
      CVLLI=CVLLHIGI+CVLLCHARI

      C1SLL=C1SLLHIG+C1SLLCHAR+C1SLLDPH
      C2SLL=C2SLLHIG+C2SLLCHAR
      C1SLLI=C1SLLHIGI+C1SLLCHARI+C1SLLDPHI
      C2SLLI=C2SLLHIGI+C2SLLCHARI

      C1LR=C1LRHIG+C1LRCHAR
      C2LR=C2LRHIG+C2LRCHAR+C2LRDPH
      C1LRI=C1LRHIGI+C1LRCHARI
      C2LRI=C2LRHIGI+C2LRCHARI+C2LRDPHI
      ENDIF

*	 4) Results for DMd

*          - `Bag' parameters from lattice [25]
      BVLL=BBd*0.551d0/0.985d0/eta**(6d0/23d0)      ! 'SM' Bag parameter from [5']
     .          /(1d0+asmb/4d0/Pi*1.6273d0*(1d0-eta))
      B1SLL=0.79d0
      B2SLL=0.70d0
      B1LR=1.72d0
      B2LR=1.15d0

      aux=1.44d0 !(M_Bd/(runmb(scb)+runmass(0.005d0,scb)))**2
      B1SLL=aux*B1SLL
      B2SLL=aux*B2SLL
      B1LR=aux*B1LR
      B2LR=aux*B2LR

*          - Running between sct and sc0 [24] Eqs.(3.1)-(3.19)

      PVLL=eta**(6d0/23d0)*(1d0+asmb/4d0/Pi*1.6273d0*(1d0-eta))*BVLL

      P1SLL=-5d0/8d0*(1.0153d0*eta**(-0.6315d0)-0.0153d0*eta**(0.7184d0)
     . +asmb/4d0/pi*(eta**(-0.6315d0)*(4.8177d0-5.2272d0*eta)
     .              +eta**(0.7184)*(0.3371d0+0.0724d0*eta)))*B1SLL
     . -3d0/2d0*(0.0081d0*(eta**(0.7184d0)-eta**(-0.6315d0))
     . +asmb/4d0/pi*(eta**(-0.6315d0)*(0.0531d0+0.0415d0*eta)
     .              -eta**(0.7184d0)*(0.0566d0+0.0380d0*eta)))*B2SLL

      P2SLL=-5d0/8d0*(1.9325d0*(eta**(-0.6315d0)-eta**(0.7184d0))
     . +asmb/4d0/pi*(eta**(-0.6315d0)*(9.1696d0-38.8778d0*eta)
     .            +eta**(0.7184d0)*(42.5021d0-12.7939d0*eta)))*B1SLL
     . -3d0/2d0*(1.0153d0*eta**(0.7184d0)-0.0153d0*eta**(-0.6315d0)
     . +asmb/4d0/pi*(eta**(-0.6315d0)*(0.1011d0+0.3083d0*eta)
     .      +eta**(0.7184d0)*(-7.1314d0+6.7219d0*eta)))*B2SLL

      P1LR=-(eta**(3d0/23d0)+asmb/4d0/Pi*(0.9250d0*eta**(-24d0/23d0)
     . +eta**(3d0/23d0)*(-2.0994d0+1.1744d0*eta)))*B1LR/2d0
     . +3d0/4d0*(2d0/3d0*(eta**(3d0/23d0)-eta**(-24d0/23d0))
     . +asmb/4d0/Pi*(eta**(3d0/23d0)*(-11.7329d0+0.7829*eta)
     .           +eta**(-24d0/23d0)*(-5.3048d0+16.2548d0*eta)))*B2LR

      P2LR=-(asmb/4d0/Pi*1.3875d0*(eta**(26d0/23d0)-eta**(-24d0/23d0)))
     .                                                      *B1LR/2d0
     . +3d0/4d0*(eta**(-24d0/23d0)+asmb/4d0/Pi*(eta**(-24d0/23d0)
     .      *(7.9572d0-8.8822d0*eta)+0.9250d0*eta**(26d0/23d0)))*B2LR

*          - Bd Mixing parameter DMd
      aux=PVLL*CVLL+P1SLL*C1SLL+P2SLL*C2SLL+P1LR*C1LR+P2LR*C2LR
      auxe=PVLL*CVLLI+P1SLL*C1SLLI+P2SLL*C2SLLI+P1LR*C1LRI+P2LR*C2LRI
      
      DMd=GF**2*MW**2/(24d0*pi**2)*M_Bd*fBd**2*VtdVtb2/(6.58211915d-13)
     .                  *dsqrt(aux**2+auxe**2)

*          - Error estimate
*      First, error bars from uncertainties on lattice Bag parameters: 
*      (2sigma, added quadratically)

      DMdMax=(dBBd*0.551d0*CVLL)**2

      aux1=-5d0/8d0*(1.0153d0*eta**(-0.6315d0)
     .                                   -0.0153d0*eta**(0.7184d0)
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(4.8177d0-5.2272d0*eta)
     .              +eta**(0.7184d0)*(0.3371d0+0.0724d0*eta)))*C1SLL
     . -5d0/8d0*(1.9325d0*(eta**(-0.6315d0)-eta**(0.7184d0))
     .      +asmb/4d0/pi*(eta**(-0.6315d0)*(9.1696d0-38.8778d0*eta)
     .            +eta**(0.7184d0)*(42.5021d0-12.7939d0*eta)))*C2SLL

      aux2=-5d0/8d0*(1.0153d0*eta**(-0.6315d0)
     .                                   -0.0153d0*eta**(0.7184d0)
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(4.8177d0-5.2272d0*eta)
     .              +eta**(0.7184d0)*(0.3371d0+0.0724d0*eta)))*C1SLLI
     . -5d0/8d0*(1.9325d0*(eta**(-0.6315d0)-eta**(0.7184d0))
     .      +asmb/4d0/pi*(eta**(-0.6315d0)*(9.1696d0-38.8778d0*eta)
     .            +eta**(0.7184d0)*(42.5021d0-12.7939d0*eta)))*C2SLLI
     
      DMdMax=DMdMax+4d0*(aux1**2+aux2**2)*(0.045d0*1.44d0)**2

      aux1=-3d0/2d0*(0.0081d0*(eta**(0.7184d0)-eta**(-0.6315d0))
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(0.0531d0+0.0415d0*eta)
     .              -eta**(0.7184d0)*(0.0566d0+0.0380d0*eta)))*C1SLL
     . -3d0/2d0*(1.0153d0*eta**(0.7184d0)-0.0153d0*eta**(-0.6315d0)
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(0.1011d0+0.3083d0*eta)
     .             +eta**(0.7184d0)*(-7.1314d0+6.7219d0*eta)))*C2SLL

      aux2=-3d0/2d0*(0.0081d0*(eta**(0.7184d0)-eta**(-0.6315d0))
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(0.0531d0+0.0415d0*eta)
     .              -eta**(0.7184d0)*(0.0566d0+0.0380d0*eta)))*C1SLLI
     . -3d0/2d0*(1.0153d0*eta**(0.7184d0)-0.0153d0*eta**(-0.6315d0)
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(0.1011d0+0.3083d0*eta)
     .             +eta**(0.7184d0)*(-7.1314d0+6.7219d0*eta)))*C2SLLI

      DMdMax=DMdMax+4d0*(aux1**2+aux2**2)*(0.103d0*1.44d0)**2

      aux1=-(eta**(3d0/23d0)+asmb/4d0/Pi*(0.9250d0*eta**(-24d0/23d0)
     .          +eta**(3d0/23d0)*(-2.0994d0+1.1744d0*eta)))/2d0*C1LR
     . -(asmb/4d0/Pi*1.3875d0*(eta**(26d0/23d0)-eta**(-24d0/23d0)))
     .                                                     *C2LR/2d0

      aux2=-(eta**(3d0/23d0)+asmb/4d0/Pi*(0.9250d0*eta**(-24d0/23d0)
     .          +eta**(3d0/23d0)*(-2.0994d0+1.1744d0*eta)))/2d0*C1LRI
     . -(asmb/4d0/Pi*1.3875d0*(eta**(26d0/23d0)-eta**(-24d0/23d0)))
     .                                                     *C2LRI/2d0

      IF(aux*aux1.gt.0d0)then
       DMdMax=DMdMax+4d0*(aux1*0.204d0*1.44d0)**2
       DMdmin=DMdMax+4d0*(aux1*0.072d0*1.44d0)**2
      ELSE
       DMdMax=DMdMax+4d0*(aux1*0.072d0*1.44d0)**2
       DMdmin=DMdMax+4d0*(aux1*0.204d0*1.44d0)**2
      ENDIF
      IF(auxe*aux2.gt.0d0)then
       DMdMax=DMdMax+4d0*(aux2*0.204d0*1.44d0)**2
       DMdmin=DMdmin+4d0*(aux2*0.072d0*1.44d0)**2
      ELSE
       DMdMax=DMdmin+4d0*(aux2*0.072d0*1.44d0)**2
       DMdmin=DMdMax+4d0*(aux2*0.204d0*1.44d0)**2
      ENDIF

      DMdMax=(GF**2*MW**2/(24d0*pi**2)*M_Bd*fBd**2*VtdVtb2
     .            /6.58211915d-13)**2*DMdMax
      DMdmin=(GF**2*MW**2/(24d0*pi**2)*M_Bd*fBd**2*VtdVtb2
     .            /6.58211915d-13)**2*DMdmin

*      Second, uncertainties from matching: 1% SM + 30% on each BSM contribution:
*      NB: SM error usually neglected in the literature
      aux=0.01d0*dabs(PVLL*CVLLSM)                                ! 1% SM

      IF(I.eq.1)then
       aux=aux+0.3d0*(dsqrt((PVLL*CVLLHIG+P1SLL*C1SLLHIG+P2SLL*C2SLLHIG
     . +P1LR*C1LRHIG+P2LR*C2LRHIG)**2+(PVLL*CVLLHIGI+P1SLL*C1SLLHIGI
     . +P2SLL*C2SLLHIGI+P1LR*C1LRHIGI+P2LR*C2LRHIGI)**2)
     . +dsqrt((PVLL*CVLLCHAR+P1SLL*C1SLLCHAR
     .        +P2SLL*C2SLLCHAR+P1LR*C1LRCHAR+P2LR*C2LRCHAR)**2
     .   +(PVLL*CVLLCHARI+P1SLL*C1SLLCHARI
     .        +P2SLL*C2SLLCHARI+P1LR*C1LRCHARI+P2LR*C2LRCHARI)**2)
     . +dsqrt((P1SLL*C1SLLDPH+P2LR*C2LRDPH)**2
     .       +(P1SLL*C1SLLDPHI+P2LR*C2LRDPHI)**2))
      ENDIF

      DMdmax=DMd+dsqrt(DMdMax)+GF**2*MW**2/(24d0*pi**2)
     .                  *M_Bd*fBd**2*VtdVtb2/(6.58211915d-13)*aux
      DMdmin=DMd-dsqrt(DMdmin)-GF**2*MW**2/(24d0*pi**2)
     .                  *M_Bd*fBd**2*VtdVtb2/(6.58211915d-13)*aux
      DMdmin=Max(0d0,DMdmin)

*      Third: errors from CKM and lattice form factor 2sigma
      aux=dsqrt((dVtdVtb2/VtdVtb2)**2+(2d0*dfBd)**2)

*      Total error bars
      DMdMax=(1d0+aux)*DMdMax
      DMdmin=(1d0-aux)*DMdmin

!      print*,'DMd',DMdmin,DMd,DMdmax
*      Comparison with experimental data (source [1',2']):
*             (0.504ps-1 < DMd=0.510ps-1 < 0.516ps-1)
      prob(34)=0d0

      IF(DMdmin.GE.DMdexpMax)
     .     PROB(34)=DMdmin/DMdexpMax-1d0
      IF(DMdmax.LE.DMdexpmin)
     .     PROB(34)=DMdmax/DMdexpMin-1d0

!      csqb=csqb+4.d0*(DMd-(DMdexpmin+DMdexpMax)/2d0)**2
!     c /((DMdMax-DMdmin)**2+(DMdexpMax-DMdexpmin)**2)
!       DMd=(P1SLL*C1SLLDPH+P2LR*C2LRDPH)/BBd


*	IV- Charged B decays: transitions mediated by a charged Higgs
*   1) BR[B+ -> tau+ nu_tau] following [20]

      rh=(1d0-mBu**2/MHC*tanb**2
     .  *(1d0-(1d0+1d0/tanb**2)/((1d0+epst0)**2+epst0I**2)
     .        *(epst0*(1d0+epst0)+epst0I*(epst0I+epst3I)
     .         +epst3I/(1d0+epst3)*(1d0+epst0)*(epst0I+epst3I)
     .         -epst0I*epst3I/(1d0+epst3))))**2
     .    +mBu**4/MHC**2*tanb**4*(1d0+1d0/tanb**2)**2
     .       /((1d0+epst0)**2+epst0I**2)**2
     .   *((1d0+epst0)*(epst0I+epst3I-epst3I/(1d0+epst3)*(epst0+epst3))
     .      -epst0I*(epst0+epst3I**2/(1d0+epst3)))**2
      aux=GF**2*mBu*MTAU**2/(8d0*pi)*(1d0-(MTAU/mBu)**2)**2
     .       *fB**2*Vub2*tauB
      BRBtaunu=aux*rh

*      Comparison with experimental data:
*  hadronic parameter (source [5']; 2sigma):
*   0.178 GeV < fB=0.1885 GeV <0.199 GeV
*  CKM factor (2sigma): 3.3 10^-3 < Vub=4.13 10^-3 < 4.7 10^-3 [6']

*  2 sigma (absolute) error bars from CKM and
*       lattice uncertainties:

      BRBtaunumax=(1d0+dfB/fB)**2*(1d0+dVub/dsqrt(Vub2))**2*BRBtaunu
      BRBtaunumin=(1d0-dfB/fB)**2*(1d0-dVub/dsqrt(Vub2))**2*BRBtaunu
!      print*,'BRBtaunu',BRBtaunumin,BRBtaunu,BRBtaunumax

      prob(36)=0d0

      IF(BRBtaunumin.GE.BRBTAUNUexpmax)
     .     PROB(36)=BRBtaunumin/BRBTAUNUexpmax-1d0
      IF(BRBtaunumax.LE.BRBTAUNUexpmin)
     .     PROB(36)=BRBtaunumax/BRBTAUNUexpmin-1d0

!      csqb=csqb+4.d0*(BRBtaunu-(BRBtaunuexpmin+BRBtaunuexpMax)/2d0)**2
!     c/((BRBtaunuMax-BRBtaunumin)**2+(BRBtaunuexpMax-BRBtaunuexpmin)**2)

*   2) RD=BR[B+ -> D tau+ nu_tau]/BR[B+ -> D l+ nu_l] following [23]
      RD_taulSM=0.297d0

      CcR=-runmb(mb)*mtau/MHC*tanb**2*(1d0
     .          -(1d0+1d0/tanb**2)/((1d0+epst0)**2+epst0I**2)
     .            *(epst0*(1d0+epst0)+epst0I*(epst0I+epst3I)
     .           +epst3I/(1d0+epst3)*(1d0+epst0)*(epst0I+epst3I)
     .           -epst0I*epst3I/(1d0+epst3)))

      CcRI=runmb(mb)*mtau/MHC*tanb**2*
     .          (1d0+1d0/tanb**2)/((1d0+epst0)**2+epst0I**2)
     .            *((1d0+epst0)*(epst0I+epst3I
     .                           -epst3I*(epst0+epst3)/(1d0+epst3))
     .              -epst0I*(epst0+epst3I**2/(1d0+epst3)))

      CcL=+runmass(1.25d0,mb)*mtau/MHC*(1d0-(1d0/tanb+tanb)
     .            *(epscb-1d0/((1d0+epst0)**2+epst0I**2)
     .           *((epscs-epscb)*(epsY32*(1d0+eps0)+eps0I*epsY32I)
     .            +(epscsI-epscbI)*(epsY32I*(1d0+eps0)-eps0I*epsY32))))

      CcLI=+runmass(1.25d0,mb)*mtau/MHC*(1d0/tanb+tanb)
     .            *(epscbI-1d0/((1d0+epst0)**2+epst0I**2)
     .           *((epscsI-epscbI)*(epsY32*(1d0+eps0)+eps0I*epsY32I)
     .            -(epscs-epscb)*(epsY32I*(1d0+eps0)-eps0I*epsY32)))

      rh=1d0+1.5d0*(CcR+CcL)+(CcR+CcL)**2+(CcRI-CcLI)**2

      RD_taul=RD_taulSM*rh

*  Error estimate
      RD_taulMax=RD_taul+2d0*0.017d0                                ! 2sigma SM
     .  +0.3d0*(dabs(CcR)+dabs(CcL))*RD_taul*dabs(1.5d0+CcR+CcL)    !30% NP
     .  +0.3d0*(dabs(CcRI)+dabs(CcLI))*dabs(CcRI-CcLI)
      RD_taulmin=RD_taul-2d0*0.017d0                                ! 2sigma SM
     .  -0.3d0*(dabs(CcR)+dabs(CcL))*RD_taul*dabs(1.5d0+CcR+CcL)    !30% NP
     .  -0.3d0*(dabs(CcRI)+dabs(CcLI))*dabs(CcRI-CcLI)
      IF(RD_taulmin.le.0d0)RD_taulmin=0d0

!      print*,'RD',RD_taulmin,RD_taul,RD_taulMax
*  Comparison with experiment [11']
      PROB(58)=0d0
      IF(RD_taulmin.GE.RD_taulexpmax)
     .     PROB(58)=RD_taulmin/RD_taulexpmax-1d0
      IF(RD_taulmax.LE.RD_taulexpmin)
     .     PROB(58)=RD_taulmax/RD_taulexpmin-1d0

*   3) RD=BR[B+ -> D* tau+ nu_tau]/BR[B+ -> D* l+ nu_l] following [23]
      RDs_taulSM=0.252d0

      rh=1d0+0.12d0*(CcR-CcL)+0.05d0*((CcR-CcL)**2+(CcRI+CcLI)**2)

      RDs_taul=RDs_taulSM*rh

*  Error estimate
      RDs_taulMax=RDs_taul+2d0*0.003d0                              ! 2sigma SM
     .  +0.3d0*(dabs(CcR)+dabs(CcL))*RD_taul                        !30% NP
     .               *dabs(0.12d0+0.05d0*(CcR-CcL))
     .  +0.3d0*(dabs(CcRI)+dabs(CcLI))*dabs(CcRI+CcLI)
      RDs_taulmin=RDs_taul-2d0*0.003d0                              ! 2sigma SM
     .  -0.3d0*(dabs(CcR)+dabs(CcL))*RD_taul*                       !30% NP
     .               *dabs(0.12d0+0.05d0*(CcR-CcL))
     .  -0.3d0*(dabs(CcRI)+dabs(CcLI))*dabs(CcRI+CcLI)
      IF(RDs_taulmin.le.0d0)RDs_taulmin=0d0

!      print*,'RD*',RDs_taulmin,RDs_taul,RDs_taulMax
*  Comparison with experiment [11']
      IF(RDs_taulmin.GE.RDs_taulexpmax)
     .     PROB(58)=PROB(58)+RDs_taulmin/RDs_taulexpmax-1d0
      IF(RDs_taulmax.LE.RDs_taulexpmin)
     .     PROB(58)=PROB(58)+RDs_taulmax/RDs_taulexpmin-1d0


*	V- Delt B = 1 - BR[B --> Xs gamma]

*   Matching scale
      sc0=160d0 !MT0
      asc0=0.109069d0 ! asf(sc0)

*   Low-energy scale
      scb=2d0 !mb
      asmb=0.293d0 ! asf(scb) ! 0.219d0    

*	 1) Wilson coefficients at sc0
*          - SM [2]
      xt=(MT0*(asf(sc0)/asf(MT0))**(12d0/23d0)
     .   *(1d0+7462d0/1587d0*(asf(sc0)-asf(MT0))/4d0/Pi)/MW)**2 !(m_t/m_W)^2

      C70SM=ff1(xt)                                      ! LO
      C80SM=fg1(xt)

      C11SM=15d0+6d0*dlog(sc0**2/MW**2)                  ! NLO
      C41SM=esm(xt)-2d0/3d0+2d0/3d0*dlog(sc0**2/MW**2)
      C71SM=GG7(xt)+Delt7(xt)*dlog(sc0**2/MW**2)
      C81SM=GG8(xt)+Delt8(xt)*dlog(sc0**2/MW**2)

*          - Charged Higgs [3]
      yt=MTH**2/MHC                                      !(m_t/m_H+)^2

      C70HIG=au**2/3d0*ff1(yt)-au*ad*ff2(yt)             ! LO
      C80HIG=au**2/3d0*fg1(yt)-au*ad*fg2(yt)

      C41HIG=au**2*eh(yt)                                ! NLO
      C71HIG=ffh(yt,tanb)-4d0/9d0*C41HIG
      C81HIG=fgh(yt,tanb)-1d0/6d0*C41HIG

      aux=(1d0-tanb*Vtb2*(epsts                              ! large tanB [1]
     .        -((epsY32*(1d0+eps0)+epsY32I*eps0I)*(epstb-epsts)
     .         -(epsY32I*(1d0+eps0)-epsY32*eps0I)*(epstbI-epstsI))
     .                                    /((1d0+eps0)**2+eps0I**2)))
     . *(1d0-(epst3+epst3I**2/(1d0+epst3))/(1d0+epst3)
     .                                        *(1d0+1d0/tanb**2))            
     . -(1d0+1d0/tanb**2)*epst3I/(1d0+epst3)**2
     .  *(epstsI-((epsY32I*(1d0+eps0)-epsY32*eps0I)*(epstb-epsts)
     .         +(epsY32*(1d0+eps0)+epsY32I*eps0I)*(epstbI-epstsI))
     .                                    /((1d0+eps0)**2+eps0I**2))
      dC7HIG=(aux-1d0)*ff2(yt)
      dC8HIG=(aux-1d0)*fg2(yt)

      aux=-(1d0-tanb*Vtb2*(epsts                            ! large tanB [1]
     .        -((epsY32*(1d0+eps0)+epsY32I*eps0I)*(epstb-epsts)
     .         -(epsY32I*(1d0+eps0)-epsY32*eps0I)*(epstbI-epstsI))
     .                                    /((1d0+eps0)**2+eps0I**2)))
     .         *(1d0+1d0/tanb**2)*epst3I/(1d0+epst3)**2
     .  -(epstsI-((epsY32I*(1d0+eps0)-epsY32*eps0I)*(epstb-epsts)
     .         +(epsY32*(1d0+eps0)+epsY32I*eps0I)*(epstbI-epstsI))
     .                                    /((1d0+eps0)**2+eps0I**2))
     .   *(1d0-(epst3+epst3I**2/(1d0+epst3))/(1d0+epst3)
     .                                        *(1d0+1d0/tanb**2)) 
      C7HIGI=aux*ff2(yt)
      C8HIGI=aux*fg2(yt)

       etaH=asmh/asc0                     ! Evolution from M_Higgs to sc0

       C7HIG=etaH**(16d0/21d0)*(C70HIG+dC7HIG)
     .       +8d0/3d0*(etaH**(2d0/3d0)-etaH**(16d0/21d0))
     .       *(C80HIG+dC8HIG)
       C8HIG=etaH**(2d0/3d0)*(C80HIG+dC8HIG)

       aux=C7HIGI
       C7HIGI=etaH**(16d0/21d0)*aux
     .       +8d0/3d0*(etaH**(2d0/3d0)-etaH**(16d0/21d0))
     .       *C8HIGI
       C8HIGI=etaH**(2d0/3d0)*C8HIGI

*          - Chargino/Squark Contributions [4,5,6]
       do j=1,2
        CCD(J,1)=U(J,2,1)*MW/(dsqrt(2d0)*COSB*dsqrt(MCH2(J)))
     .                      /(1d0+epst3)
        CCD(J,2)=U(J,2,2)*MW/(dsqrt(2d0)*COSB*dsqrt(MCH2(J)))
     .                      /(1d0+epst3)
       do k=1,2
        CCT(j,k,1)=(V(j,1,1)*UT(k,1,1)+V(j,1,2)*UT(k,1,2))
     . *(1d0-asf(dsqrt(MST2(k)))/(4d0*Pi)
     .                       *(7d0/3d0+2d0*dlog(MST2(K)/MGL**2)))
     .     -(V(j,2,1)*UT(k,2,1)+V(j,2,2)*UT(k,2,2))*Ytq/dsqrt(g2)
     . *(1d0+asf(dsqrt(MST2(k)))/(4d0*Pi)
     .                       *(1d0+2d0*dlog(MST2(K)/MGL**2)))
        CCT(j,k,2)=(V(j,1,1)*UT(k,1,2)-V(j,1,2)*UT(k,1,1))
     . *(1d0-asf(dsqrt(MST2(k)))/(4d0*Pi)
     .                       *(7d0/3d0+2d0*dlog(MST2(K)/MGL**2)))
     .     -(V(j,2,1)*UT(k,2,2)-V(j,2,2)*UT(k,2,1))*Ytq/dsqrt(g2)
     . *(1d0+asf(dsqrt(MST2(k)))/(4d0*Pi)
     .                       *(1d0+2d0*dlog(MST2(K)/MGL**2)))
       enddo
       enddo

       C7CHARS=0d0
       C8CHARS=0d0
       C7CHARSI=0d0
       C8CHARSI=0d0

       do j=1,2                            ! LO Scharm/charginos
        C7CHARS=C7CHARS-VVc*(
     .          2d0/3d0*(V(J,1,1)**2+V(J,1,2)**2)
     .                             *MW**2/MSU2(1)*FF1(MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                 *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .         +(CCD(J,1)*V(J,1,1)+CCD(J,2)*V(J,1,2))
     .                             *FF3(MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .               *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))
     .   *(1d0+asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .               *(1d0+2d0*dlog(MSU2(1)/MGL**2))))
        C8CHARS=C8CHARS-VVc*(
     .          2d0/3d0*(V(J,1,1)**2+V(J,1,2)**2)
     .                             *MW**2/MSU2(1)*FG1(MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                 *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .       +(CCD(J,1)*V(J,1,1)+CCD(J,2)*V(J,1,2))
     .                             *FG3(MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .               *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))
     .   *(1d0+asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .               *(1d0+2d0*dlog(MSU2(1)/MGL**2))))
        C7CHARSI=C7CHARSI-VVc*(
     .         -(CCD(J,1)*V(J,1,2)+CCD(J,2)*V(J,1,1))
     .                             *FF3(MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .               *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))
     .   *(1d0+asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .               *(1d0+2d0*dlog(MSU2(1)/MGL**2))))
        C8CHARSI=C8CHARSI-VVc*(
     .          -(CCD(J,1)*V(J,1,2)+CCD(J,2)*V(J,1,1))
     .                             *FG3(MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .               *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))
     .   *(1d0+asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .               *(1d0+2d0*dlog(MSU2(1)/MGL**2))))
       enddo
       etaS=asf(dsqrt(MSU2(1)))/asc0       ! Running between MUL and sc0
       C7CHAR=etaS**(16d0/21d0)*C7CHARS+
     .    8d0/3d0*(etaS**(14d0/21d0)-etaS**(16d0/21d0))*C8CHARS
       C8CHAR=etaS**(14d0/21d0)*C8CHARS
       C7CHARI=etaS**(16d0/21d0)*C7CHARSI+
     .    8d0/3d0*(etaS**(14d0/21d0)-etaS**(16d0/21d0))*C8CHARSI
       C8CHARI=etaS**(14d0/21d0)*C8CHARSI

       do k=1,2                            ! LO Stop/charginos
       C7CHARS=0d0
       C8CHARS=0d0
       C7CHARSI=0d0
       C8CHARSI=0d0
       do j=1,2
         C7CHARS=C7CHARS-2d0/3d0*(CCT(J,K,1)**2+CCT(J,K,2)**2)
     .                           *MW**2/MST2(K)*FF1(MST2(K)/MCH2(J))
     .  -(CCD(J,1)*UT(K,1,1)*CCT(J,K,1)-CCD(J,2)*UT(K,1,2)*CCT(J,K,1)
     .   +CCD(J,1)*UT(K,1,2)*CCT(J,K,2)+CCD(J,2)*UT(K,1,1)*CCT(J,K,2))
     .                                         *FF3(MST2(K)/MCH2(J))
     . *(1d0+asf(dsqrt(MST2(k)))/(4d0*Pi)
     .                              *(1d0+2d0*dlog(MST2(K)/MGL**2)))
         C8CHARS=C8CHARS-2d0/3d0*(CCT(J,K,1)**2+CCT(J,K,2)**2)
     .                           *MW**2/MST2(K)*FG1(MST2(K)/MCH2(J))
     .  -(CCD(J,1)*UT(K,1,1)*CCT(J,K,1)-CCD(J,2)*UT(K,1,2)*CCT(J,K,1)
     .   +CCD(J,1)*UT(K,1,2)*CCT(J,K,2)+CCD(J,2)*UT(K,1,1)*CCT(J,K,2))
     .                                         *FG3(MST2(K)/MCH2(J))
     . *(1d0+asf(dsqrt(MST2(k)))/(4d0*Pi)
     .                              *(1d0+2d0*dlog(MST2(K)/MGL**2)))
         C7CHARSI=C7CHARSI
     .  +(CCD(J,2)*UT(K,1,1)*CCT(J,K,1)+CCD(J,1)*UT(K,1,2)*CCT(J,K,1)
     .   +CCD(J,2)*UT(K,1,2)*CCT(J,K,2)-CCD(J,2)*UT(K,1,1)*CCT(J,K,2))
     .                                         *FF3(MST2(K)/MCH2(J))
     . *(1d0+asf(dsqrt(MST2(k)))/(4d0*Pi)
     .                              *(1d0+2d0*dlog(MST2(K)/MGL**2)))
         C8CHARSI=C8CHARSI
     .  +(CCD(J,2)*UT(K,1,1)*CCT(J,K,1)+CCD(J,2)*UT(K,1,2)*CCT(J,K,1)
     .   +CCD(J,2)*UT(K,1,2)*CCT(J,K,2)-CCD(J,2)*UT(K,1,1)*CCT(J,K,2))
     .                                         *FG3(MST2(K)/MCH2(J))
     . *(1d0+asf(dsqrt(MST2(k)))/(4d0*Pi)
     .                              *(1d0+2d0*dlog(MST2(K)/MGL**2)))
       enddo
       etaS= asf(dsqrt(MST2(k)))/asc0      ! Running between ST(K) and sc0
       C7CHAR=C7CHAR+etaS**(16d0/21d0)*C7CHARS+
     .    8d0/3d0*(etaS**(14d0/21d0)-etaS**(16d0/21d0))*C8CHARS
       C8CHAR=C8CHAR+etaS**(14d0/21d0)*C8CHARS
       C7CHARI=C7CHARI+etaS**(16d0/21d0)*C7CHARSI+
     .    8d0/3d0*(etaS**(14d0/21d0)-etaS**(16d0/21d0))*C8CHARSI
       C8CHARI=C8CHARI+etaS**(14d0/21d0)*C8CHARSI
       enddo

       C41CHAR=0d0
       C71CHAR=0d0
       C81CHAR=0d0
       C71CHARI=0d0
       C81CHARI=0d0

       do j=1,2                            ! NLO Scharm/chargino
       C41CHAR=C41CHAR+VVc*(
     .   (V(J,1,1)**2+V(J,1,2)**2)*MW**2/MSU2(1)*echi(MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))**2)
       C71CHAR=C71CHAR+VVc*(
     .   (V(J,1,1)**2+V(J,1,2)**2)*MW**2/MCH2(J)*H17(MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .   -(CCD(J,1)*V(J,1,1)-CCD(J,2)*V(J,1,2))*H27(MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))
     .   *(1d0+asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(1d0+2d0*dlog(MSU2(1)/MGL**2))))
       C81CHAR=C81CHAR+VVc*(
     .   (V(J,1,1)**2+V(J,1,2)**2)*MW**2/MCH2(J)*H18(MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .   -(CCD(J,1)*V(J,1,1)-CCD(J,2)*V(J,1,2))*H28(MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))
     .   *(1d0+asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(1d0+2d0*dlog(MSU2(1)/MGL**2))))
       C71CHARI=C71CHARI+VVc*
     .   (CCD(J,1)*V(J,1,2)+CCD(J,2)*V(J,1,1))*H27(MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))
     .   *(1d0+asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(1d0+2d0*dlog(MSU2(1)/MGL**2)))
       C81CHARI=C81CHARI+VVc*
     .   (CCD(J,1)*V(J,1,2)+CCD(J,2)*V(J,1,1))*H28(MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))
     .   *(1d0+asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(1d0+2d0*dlog(MSU2(1)/MGL**2)))
       enddo

       do j=1,2                            ! NLO Scharm quartic coupling
       C41CHAR=C41CHAR+VVc*((V(J,1,1)**2+V(J,1,2)**2)*MW**2/MSU2(1)
     .               *esfe(MSU2(1)/MCH2(J),MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))**2)
       C71CHAR=C71CHAR+VVc/2d0*((V(J,1,1)**2+V(J,1,2)**2)*MW**2/MSU2(1)
     .   *(Q11(MSU2(1)/MCH2(J),MSU2(1)/MCH2(J))
     .           -2d0/3d0*Q21(MSU2(1)/MCH2(J),MSU2(1)/MCH2(J)))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .         -(CCD(J,1)*V(J,1,1)-CCD(J,2)*V(J,1,2))*MCH2(J)/MSU2(1)
     .   *(Q31(MSU2(1)/MCH2(J),MSU2(1)/MCH2(J))
     .           -2d0/3d0*Q41(MSU2(1)/MCH2(J),MSU2(1)/MCH2(J)))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))
     .   *(1d0+asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(1d0+2d0*dlog(MSU2(1)/MGL**2))))
       C81CHAR=C81CHAR-VVc/2d0*((V(J,1,1)**2+V(J,1,2)**2)*MW**2/MSU2(1)
     .                        *Q21(MSU2(1)/MCH2(J),MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .         -(CCD(J,1)*V(J,1,1)-CCD(J,2)*V(J,1,2))*MCH2(J)/MSU2(1)
     .                        *Q41(MSU2(1)/MCH2(J),MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))
     .   *(1d0+asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(1d0+2d0*dlog(MSU2(1)/MGL**2))))
       C71CHARI=C71CHARI+VVc/2d0*
     .         (CCD(J,1)*V(J,1,2)+CCD(J,2)*V(J,1,1))*MCH2(J)/MSU2(1)
     .   *(Q31(MSU2(1)/MCH2(J),MSU2(1)/MCH2(J))
     .           -2d0/3d0*Q41(MSU2(1)/MCH2(J),MSU2(1)/MCH2(J)))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))
     .   *(1d0+asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(1d0+2d0*dlog(MSU2(1)/MGL**2)))
       C81CHARI=C81CHARI-VVc/2d0*
     .         (CCD(J,1)*V(J,1,2)+CCD(J,2)*V(J,1,1))*MCH2(J)/MSU2(1)
     .                        *Q41(MSU2(1)/MCH2(J),MSU2(1)/MCH2(J))
     .   *(1d0-asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(7d0/3d0+2d0*dlog(MSU2(1)/MGL**2)))
     .   *(1d0+asf(dsqrt(MSU2(1)))/(4d0*Pi)
     .                    *(1d0+2d0*dlog(MSU2(1)/MGL**2)))
       enddo

       do j=1,2                            ! NLO Stop/chargino
       do k=1,2
       C41CHAR=C41CHAR+(CCT(j,k,1)**2+CCT(j,k,2)**2)
     .              *MW**2/MST2(K)*echi(MST2(K)/MCH2(J))
       C71CHAR=C71CHAR+(CCT(j,k,1)**2+CCT(j,k,2)**2)
     .              *MW**2/MCH2(J)*H17(MST2(K)/MCH2(J))
     .   -(CCD(J,1)*UT(k,1,1)*CCT(J,K,1)+CCD(J,2)*UT(k,1,1)*CCT(J,K,2)
     .    -CCD(J,2)*UT(k,1,2)*CCT(J,K,1)+CCD(J,1)*UT(k,1,2)*CCT(J,K,2))
     .                         *H27(MST2(K)/MCH2(J))
     .   *(1d0+asf(dsqrt(MST2(K)))/(4d0*Pi)
     .                        *(1d0+2d0*dlog(MST2(K)/MGL**2)))
       C81CHAR=C81CHAR+(CCT(j,k,1)**2+CCT(j,k,2)**2)
     .              *MW**2/MCH2(J)*H18(MST2(K)/MCH2(J))
     .   -(CCD(J,1)*UT(k,1,1)*CCT(J,K,1)+CCD(J,2)*UT(k,1,1)*CCT(J,K,2)
     .    -CCD(J,2)*UT(k,1,2)*CCT(J,K,1)+CCD(J,1)*UT(k,1,2)*CCT(J,K,2))
     .                         *H28(MST2(K)/MCH2(J))
     .   *(1d0+asf(dsqrt(MST2(K)))/(4d0*Pi)
     .                        *(1d0+2d0*dlog(MST2(K)/MGL**2)))
       C71CHARI=C71CHARI+
     .   (CCD(J,2)*UT(k,1,1)*CCT(J,K,1)+CCD(J,1)*UT(k,1,2)*CCT(J,K,1)
     .   -CCD(J,1)*UT(k,1,1)*CCT(J,K,2)+CCD(J,2)*UT(k,1,2)*CCT(J,K,2))
     .                         *H27(MST2(K)/MCH2(J))
     .   *(1d0+asf(dsqrt(MST2(K)))/(4d0*Pi)
     .                        *(1d0+2d0*dlog(MST2(K)/MGL**2)))
       C81CHARI=C81CHARI+
     .   (CCD(J,2)*UT(k,1,1)*CCT(J,K,1)+CCD(J,1)*UT(k,2,1)*CCT(J,K,2)
     .   -CCD(J,1)*UT(k,1,1)*CCT(J,K,2)+CCD(J,2)*UT(k,1,2)*CCT(J,K,2))
     .                         *H28(MST2(K)/MCH2(J))
     .   *(1d0+asf(dsqrt(MST2(K)))/(4d0*Pi)
     .                        *(1d0+2d0*dlog(MST2(K)/MGL**2)))
       enddo
       enddo

       do j=1,2                            ! NLO Stop quartic coupling
       do k=1,2
       do l=1,2
       do m=1,2
       C41CHAR=C41CHAR+MW**2/MST2(l)
     . *((CCT(k,J,1)*CCT(m,j,1)+CCT(k,J,2)*CCT(m,j,2))
     .   *((UT(k,1,1)*UT(l,1,1)+UT(k,1,2)*UT(l,1,2)
     .              -UT(k,2,1)*UT(l,2,1)-UT(k,2,2)*UT(l,2,2))
     .     *(UT(l,1,1)*UT(m,1,1)+UT(l,1,2)*UT(m,1,2)
     .              -UT(l,2,1)*UT(m,2,1)-UT(l,2,2)*UT(m,2,2))
     .    +(UT(k,1,1)*UT(l,1,2)-UT(k,1,2)*UT(l,1,1)
     .              -UT(k,2,1)*UT(l,2,2)+UT(k,2,2)*UT(l,2,1))
     .     *(UT(l,1,2)*UT(m,1,1)-UT(l,1,1)*UT(m,1,2)
     .              -UT(l,2,2)*UT(m,2,1)+UT(l,2,1)*UT(m,2,2)))
     .  +(CCT(k,J,2)*CCT(m,j,1)-CCT(k,J,1)*CCT(m,j,2))
     .   *((UT(k,1,1)*UT(l,1,1)+UT(k,1,2)*UT(l,1,2)
     .              -UT(k,2,1)*UT(l,2,1)-UT(k,2,2)*UT(l,2,2))
     .     *(UT(l,1,2)*UT(m,1,1)-UT(l,1,1)*UT(m,1,2)
     .              -UT(l,2,2)*UT(m,2,1)+UT(l,2,1)*UT(m,2,2))
     .    -(UT(k,1,1)*UT(l,1,2)-UT(k,1,2)*UT(l,1,1)
     .              -UT(k,2,1)*UT(l,2,2)+UT(k,2,2)*UT(l,2,1))
     .     *(UT(l,1,1)*UT(m,1,1)+UT(l,1,2)*UT(m,1,2)
     .              -UT(l,2,1)*UT(m,2,1)-UT(l,2,2)*UT(m,2,2))))
     .              *esfe(MST2(k)/MCH2(J),MST2(m)/MCH2(J))
       C71CHAR=C71CHAR+MW**2/MST2(l)/2d0*(
     .    ((UT(k,1,1)*UT(l,1,1)+UT(k,1,2)*UT(l,1,2)
     .              -UT(k,2,1)*UT(l,2,1)-UT(k,2,2)*UT(l,2,2))
     .     *(UT(l,1,1)*UT(m,1,1)+UT(l,1,2)*UT(m,1,2)
     .              -UT(l,2,1)*UT(m,2,1)-UT(l,2,2)*UT(m,2,2))
     .    +(UT(k,1,1)*UT(l,1,2)-UT(k,1,2)*UT(l,1,1)
     .              -UT(k,2,1)*UT(l,2,2)+UT(k,2,2)*UT(l,2,1))
     .     *(UT(l,1,2)*UT(m,1,1)-UT(l,1,1)*UT(m,1,2)
     .              -UT(l,2,2)*UT(m,2,1)+UT(l,2,1)*UT(m,2,2)))
     .   *((CCT(k,J,1)*CCT(m,j,1)+CCT(k,J,2)*CCT(m,j,2))
     .     *(Q11(MST2(k)/MCH2(J),MST2(m)/MCH2(J))
     .               -2d0/3d0*Q21(MST2(k)/MCH2(J),MST2(m)/MCH2(J)))
     .    -(CCT(k,j,1)*CCD(j,1)*UT(m,1,1)+CCT(k,j,2)*CCD(j,2)*UT(m,1,1)
     .     +CCT(k,j,2)*CCD(j,1)*UT(m,1,2)-CCT(k,j,1)*CCD(j,2)*UT(m,1,2))
     .   *(1d0+asf(dsqrt(MST2(K)))/(4d0*Pi)
     .                   *(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .   *MCH2(j)/MST2(l)*(Q31(MST2(k)/MCH2(J),MST2(m)/MCH2(J))
     .               -2d0/3d0*Q41(MST2(k)/MCH2(J),MST2(m)/MCH2(J))))
     . +((UT(k,1,1)*UT(l,1,1)+UT(k,1,2)*UT(l,1,2)
     .              -UT(k,2,1)*UT(l,2,1)-UT(k,2,2)*UT(l,2,2))
     .     *(UT(l,1,2)*UT(m,1,1)-UT(l,1,1)*UT(m,1,2)
     .              -UT(l,2,2)*UT(m,2,1)+UT(l,2,1)*UT(m,2,2))
     .    -(UT(k,1,1)*UT(l,1,2)-UT(k,1,2)*UT(l,1,1)
     .              -UT(k,2,1)*UT(l,2,2)+UT(k,2,2)*UT(l,2,1))
     .     *(UT(l,1,1)*UT(m,1,1)+UT(l,1,2)*UT(m,1,2)
     .              -UT(l,2,1)*UT(m,2,1)-UT(l,2,2)*UT(m,2,2)))
     .   *((CCT(k,J,2)*CCT(m,j,1)-CCT(k,J,1)*CCT(m,j,2))
     .     *(Q11(MST2(k)/MCH2(J),MST2(m)/MCH2(J))
     .               -2d0/3d0*Q21(MST2(k)/MCH2(J),MST2(m)/MCH2(J)))
     .    -(CCT(k,j,2)*CCD(j,1)*UT(m,1,1)-CCT(k,j,2)*CCD(j,2)*UT(m,1,2)
     .     -CCT(k,j,1)*CCD(j,2)*UT(m,1,1)-CCT(k,j,1)*CCD(j,1)*UT(m,1,2))
     .   *(1d0+asf(dsqrt(MST2(K)))/(4d0*Pi)
     .                   *(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .   *MCH2(j)/MST2(l)*(Q31(MST2(k)/MCH2(J),MST2(m)/MCH2(J))
     .               -2d0/3d0*Q41(MST2(k)/MCH2(J),MST2(m)/MCH2(J)))))
       C81CHAR=C81CHAR-MW**2/MST2(l)/2d0*(
     .    ((UT(k,1,1)*UT(l,1,1)+UT(k,1,2)*UT(l,1,2)
     .              -UT(k,2,1)*UT(l,2,1)-UT(k,2,2)*UT(l,2,2))
     .     *(UT(l,1,1)*UT(m,1,1)+UT(l,1,2)*UT(m,1,2)
     .              -UT(l,2,1)*UT(m,2,1)-UT(l,2,2)*UT(m,2,2))
     .    +(UT(k,1,1)*UT(l,1,2)-UT(k,1,2)*UT(l,1,1)
     .              -UT(k,2,1)*UT(l,2,2)+UT(k,2,2)*UT(l,2,1))
     .     *(UT(l,1,2)*UT(m,1,1)-UT(l,1,1)*UT(m,1,2)
     .              -UT(l,2,2)*UT(m,2,1)+UT(l,2,1)*UT(m,2,2)))
     .   *((CCT(k,J,1)*CCT(m,j,1)+CCT(k,J,2)*CCT(m,j,2))
     .                *Q21(MST2(k)/MCH2(J),MST2(m)/MCH2(J))
     .   -(CCT(k,j,1)*CCD(j,1)*UT(m,1,1)+CCT(k,j,2)*CCD(j,2)*UT(m,1,1)
     .    +CCT(k,j,2)*CCD(j,1)*UT(m,1,2)-CCT(k,j,1)*CCD(j,2)*UT(m,1,2))
     .   *(1d0+asf(dsqrt(MST2(K)))/(4d0*Pi)
     .                   *(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .   *MCH2(j)/MST2(l)*Q41(MST2(k)/MCH2(J),MST2(m)/MCH2(J)))
     . +((UT(k,1,1)*UT(l,1,1)+UT(k,1,2)*UT(l,1,2)

     .              -UT(k,2,1)*UT(l,2,1)-UT(k,2,2)*UT(l,2,2))
     .     *(UT(l,1,2)*UT(m,1,1)-UT(l,1,1)*UT(m,1,2)
     .              -UT(l,2,2)*UT(m,2,1)+UT(l,2,1)*UT(m,2,2))
     .    -(UT(k,1,1)*UT(l,1,2)-UT(k,1,2)*UT(l,1,1)
     .              -UT(k,2,1)*UT(l,2,2)+UT(k,2,2)*UT(l,2,1))
     .     *(UT(l,1,1)*UT(m,1,1)+UT(l,1,2)*UT(m,1,2)
     .              -UT(l,2,1)*UT(m,2,1)-UT(l,2,2)*UT(m,2,2)))
     .   *((CCT(k,J,2)*CCT(m,j,1)-CCT(k,J,1)*CCT(m,j,2))
     .                *Q21(MST2(k)/MCH2(J),MST2(m)/MCH2(J))
     .   -(CCT(k,j,2)*CCD(j,1)*UT(m,1,1)-CCT(k,j,2)*CCD(j,2)*UT(m,1,2)
     .    -CCT(k,j,1)*CCD(j,1)*UT(m,1,2)-CCT(k,j,1)*CCD(j,2)*UT(m,1,1))
     .   *(1d0+asf(dsqrt(MST2(K)))/(4d0*Pi)
     .                   *(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .   *MCH2(j)/MST2(l)*Q41(MST2(k)/MCH2(J),MST2(m)/MCH2(J))))
       C71CHARI=C71CHARI+MW**2/MST2(l)/2d0*(
     .    -((UT(k,1,1)*UT(l,1,1)+UT(k,1,2)*UT(l,1,2)
     .              -UT(k,2,1)*UT(l,2,1)-UT(k,2,2)*UT(l,2,2))
     .     *(UT(l,1,2)*UT(m,1,1)-UT(l,1,1)*UT(m,1,2)
     .              -UT(l,2,2)*UT(m,2,1)+UT(l,2,1)*UT(m,2,2))
     .    -(UT(k,1,1)*UT(l,1,2)-UT(k,1,2)*UT(l,1,1)
     .              -UT(k,2,1)*UT(l,2,2)+UT(k,2,2)*UT(l,2,1))
     .     *(UT(l,1,1)*UT(m,1,1)+UT(l,1,2)*UT(m,1,2)
     .              -UT(l,2,1)*UT(m,2,1)-UT(l,2,2)*UT(m,2,2)))
     .   *((CCT(k,J,1)*CCT(m,j,1)+CCT(k,J,2)*CCT(m,j,2))
     .     *(Q11(MST2(k)/MCH2(J),MST2(m)/MCH2(J))
     .               -2d0/3d0*Q21(MST2(k)/MCH2(J),MST2(m)/MCH2(J)))
     .    -(CCT(k,j,1)*CCD(j,1)*UT(m,1,1)+CCT(k,j,2)*CCD(j,2)*UT(m,1,1)
     .     +CCT(k,j,2)*CCD(j,1)*UT(m,1,2)-CCT(k,j,1)*CCD(j,2)*UT(m,1,2))
     .   *(1d0+asf(dsqrt(MST2(K)))/(4d0*Pi)
     .                   *(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .   *MCH2(j)/MST2(l)*(Q31(MST2(k)/MCH2(J),MST2(m)/MCH2(J))
     .               -2d0/3d0*Q41(MST2(k)/MCH2(J),MST2(m)/MCH2(J))))
     . +((UT(k,1,1)*UT(l,1,1)+UT(k,1,2)*UT(l,1,2)
     .              -UT(k,2,1)*UT(l,2,1)-UT(k,2,2)*UT(l,2,2))
     .     *(UT(l,1,1)*UT(m,1,1)+UT(l,1,2)*UT(m,1,2)
     .              -UT(l,2,1)*UT(m,2,1)-UT(l,2,2)*UT(m,2,2))
     .    +(UT(k,1,1)*UT(l,1,2)-UT(k,1,2)*UT(l,1,1)
     .              -UT(k,2,1)*UT(l,2,2)+UT(k,2,2)*UT(l,2,1))
     .     *(UT(l,1,2)*UT(m,1,1)-UT(l,1,1)*UT(m,1,2)
     .              -UT(l,2,2)*UT(m,2,1)+UT(l,2,1)*UT(m,2,2)))
     .   *((CCT(k,J,2)*CCT(m,j,1)-CCT(k,J,1)*CCT(m,j,2))
     .     *(Q11(MST2(k)/MCH2(J),MST2(m)/MCH2(J))
     .               -2d0/3d0*Q21(MST2(k)/MCH2(J),MST2(m)/MCH2(J)))
     .    -(CCT(k,j,2)*CCD(j,1)*UT(m,1,1)-CCT(k,j,2)*CCD(j,2)*UT(m,1,2)
     .     -CCT(k,j,1)*CCD(j,2)*UT(m,1,1)-CCT(k,j,1)*CCD(j,1)*UT(m,1,2))
     .   *(1d0+asf(dsqrt(MST2(K)))/(4d0*Pi)
     .                   *(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .   *MCH2(j)/MST2(l)*(Q31(MST2(k)/MCH2(J),MST2(m)/MCH2(J))
     .               -2d0/3d0*Q41(MST2(k)/MCH2(J),MST2(m)/MCH2(J)))))
       C81CHARI=C81CHARI-MW**2/MST2(l)/2d0*(
     .    -((UT(k,1,1)*UT(l,1,1)+UT(k,1,2)*UT(l,1,2)
     .              -UT(k,2,1)*UT(l,2,1)-UT(k,2,2)*UT(l,2,2))
     .     *(UT(l,1,2)*UT(m,1,1)-UT(l,1,1)*UT(m,1,2)
     .              -UT(l,2,2)*UT(m,2,1)+UT(l,2,1)*UT(m,2,2))
     .    -(UT(k,1,1)*UT(l,1,2)-UT(k,1,2)*UT(l,1,1)
     .              -UT(k,2,1)*UT(l,2,2)+UT(k,2,2)*UT(l,2,1))
     .     *(UT(l,1,1)*UT(m,1,1)+UT(l,1,2)*UT(m,1,2)
     .              -UT(l,2,1)*UT(m,2,1)-UT(l,2,2)*UT(m,2,2)))
     .   *((CCT(k,J,1)*CCT(m,j,1)+CCT(k,J,2)*CCT(m,j,2))
     .                *Q21(MST2(k)/MCH2(J),MST2(m)/MCH2(J))
     .   -(CCT(k,j,1)*CCD(j,1)*UT(m,1,1)+CCT(k,j,2)*CCD(j,2)*UT(m,1,1)
     .    +CCT(k,j,2)*CCD(j,1)*UT(m,1,2)-CCT(k,j,1)*CCD(j,2)*UT(m,1,2))
     .   *(1d0+asf(dsqrt(MST2(K)))/(4d0*Pi)
     .                   *(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .   *MCH2(j)/MST2(l)*Q41(MST2(k)/MCH2(J),MST2(m)/MCH2(J)))
     .    +((UT(k,1,1)*UT(l,1,1)+UT(k,1,2)*UT(l,1,2)
     .              -UT(k,2,1)*UT(l,2,1)-UT(k,2,2)*UT(l,2,2))
     .     *(UT(l,1,1)*UT(m,1,1)+UT(l,1,2)*UT(m,1,2)
     .              -UT(l,2,1)*UT(m,2,1)-UT(l,2,2)*UT(m,2,2))
     .    +(UT(k,1,1)*UT(l,1,2)-UT(k,1,2)*UT(l,1,1)
     .              -UT(k,2,1)*UT(l,2,2)+UT(k,2,2)*UT(l,2,1))
     .     *(UT(l,1,2)*UT(m,1,1)-UT(l,1,1)*UT(m,1,2)
     .              -UT(l,2,2)*UT(m,2,1)+UT(l,2,1)*UT(m,2,2)))
     .   *((CCT(k,J,2)*CCT(m,j,1)-CCT(k,J,1)*CCT(m,j,2))
     .                *Q21(MST2(k)/MCH2(J),MST2(m)/MCH2(J))
     .   -(CCT(k,j,2)*CCD(j,1)*UT(m,1,1)-CCT(k,j,2)*CCD(j,2)*UT(m,1,2)
     .    -CCT(k,j,1)*CCD(j,1)*UT(m,1,2)-CCT(k,j,1)*CCD(j,2)*UT(m,1,1))
     .   *(1d0+asf(dsqrt(MST2(K)))/(4d0*Pi)
     .                   *(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .   *MCH2(j)/MST2(l)*Q41(MST2(k)/MCH2(J),MST2(m)/MCH2(J))))
       enddo
       enddo
       enddo
       enddo

*          - Higgs Penguin Contribution [1]

      C70S0=0d0
      C80S0=0d0
      C70S0I=0d0
      C80S0I=0d0

      do i=1,5
       aux1=(XH(i,1)-XH(i,2)*tanb)
     .   *(XH(i,2)*(1d0-epst3I**2/(1d0+epst3))
     .    +XH(i,1)/tanb*(epst3+epst3I**2/(1d0+epst3))
     .    +XH(i,4)*epst3I/(1d0+epst3)*sinb*(1d0+1d0/tanb**2))
     .  -XH(i,4)*cosb*(1d0+1/tanb**2)
     .   *(XH(i,4)*sinb*(1d0-epst3I**2/(1d0+epst3))
     .    -XH(i,4)*cosb/tanb*(epst3+epst3I**2/(1d0+epst3))
     .    -epst3I/(1d0+epst3)*(XH(i,2)-XH(i,1)/tanb))
       aux2=-XH(i,4)*cosb*(1d0+1/tanb**2)
     .   *(XH(i,2)*(1d0-epst3I**2/(1d0+epst3))
     .    +XH(i,1)/tanb*(epst3+epst3I**2/(1d0+epst3))
     .    +XH(i,4)*epst3I/(1d0+epst3)*sinb*(1d0+1d0/tanb**2))
     .  -(XH(i,1)-XH(i,2)*tanb)
     .   *(XH(i,4)*sinb*(1d0-epst3I**2/(1d0+epst3))
     .    -XH(i,4)*cosb/tanb*(epst3+epst3I**2/(1d0+epst3))
     .    -epst3I/(1d0+epst3)*(XH(i,2)-XH(i,1)/tanb))
       aux=1d0/18d0*MW**2/g2*MB0/(vdq*(1d0+epst3))/MH0(i)
     .                        *(sigRLbs*aux1-sigRLbsI*aux2)
       auxe=1d0/18d0*MW**2/g2*MB0/(vdq*(1d0+epst3))/MH0(i)
     .                        *(sigRLbsI*aux1+sigRLbs*aux2)

       IF(dsqrt(MH0(i)).ge.mb)THEN                    ! Running to sc0
        eta0=asf(dsqrt(MH0(i)))/asc0
       ELSE
        eta0=asmb/asc0
       ENDIF

       IF(dsqrt(MH0(i)).gt.mt0)then
        C80S0=C80S0-3d0*eta0**(14d0/21d0)*aux
        C70S0=C70S0+eta0**(16d0/21d0)*aux
     . -3d0*8d0/3d0*(eta0**(14d0/21d0)-eta0**(16d0/21d0))*aux
        C80S0I=C80S0I-3d0*eta0**(14d0/21d0)*auxe
        C70S0I=C70S0I+eta0**(16d0/21d0)*auxe
     . -3d0*8d0/3d0*(eta0**(14d0/21d0)-eta0**(16d0/21d0))*auxe
       ELSE
        C80S0=C80S0-3d0*eta0**(14d0/23d0)*aux
        C70S0=C70S0+eta0**(16d0/23d0)*aux
     . -3d0*8d0/3d0*(eta0**(14d0/23d0)-eta0**(16d0/23d0))*aux
        C80S0I=C80S0I-3d0*eta0**(14d0/23d0)*auxe
        C70S0I=C70S0I+eta0**(16d0/23d0)*auxe
     . -3d0*8d0/3d0*(eta0**(14d0/23d0)-eta0**(16d0/23d0))*auxe
       ENDIF
      enddo

*          - Summary
      I=1                                             ! I=0: SM; I=1: NMSSM
      IF(I.eq.0)THEN               ! sets all BSM contributions to 0
       C70BSM=0d0
       C80BSM=0d0
       C41BSM=0d0
       C71BSM=0d0
       C81BSM=0d0
       C70BSMI=0d0
       C80BSMI=0d0
       C71BSMI=0d0
       C81BSMI=0d0
      ELSE                         ! collects all BSM contributions
       C70BSM=C7HIG+C7CHAR+C70S0
       C80BSM=C8HIG+C8CHAR+C80S0
       C41BSM=C41HIG+C41CHAR
       C71BSM=C71HIG+C71CHAR
       C81BSM=C81HIG+C81CHAR
       C70BSMI=C7HIGI+C7CHARI+C70S0I
       C80BSMI=C8HIGI+C8CHARI+C80S0I
       C71BSMI=C71CHARI
       C81BSMI=C81CHARI
      ENDIF
                 
*  DC7BSM and DC8BSM are (conservative) error estimates:
*  Dominant errors: HIG: Order (alphas**2) ~ 10%,
*  CHAR: Order(alphas**2) ~ 10%, PENGUIN: Order(alphas) ~ 30%
       DC7BSM=1d-1*DABS(C7HIG)+1d-1*DABS(C7CHAR)+3d-1*DABS(C70S0)
     .      +asc0/4d0/Pi*(3d-1*DABS(C71HIG)+3d-1*DABS(C71CHAR))
       DC8BSM=1d-1*DABS(C8HIG)+1d-1*DABS(C8CHAR)+3d-1*DABS(C80S0)
     .      +asc0/4d0/Pi*(3d-1*DABS(C81HIG)+3d-1*DABS(C81CHAR))
       DC7BSMI=1d-1*DABS(C7HIGI)+1d-1*DABS(C7CHARI)+3d-1*DABS(C70S0I)
     .      +asc0/4d0/Pi*3d-1*DABS(C71CHARI)
       DC8BSMI=1d-1*DABS(C8HIGI)+1d-1*DABS(C8CHARI)+3d-1*DABS(C80S0I)
     .      +asc0/4d0/Pi*3d-1*DABS(C81CHARI)

*  Finally, the total Wilson coefficients at sc0
       C70=C70SM+C70BSM
       C80=C80SM+C80BSM
       C11=C11SM
       C41=C41SM+C41BSM
       C71=C71SM+C71BSM
       C81=C81SM+C81BSM
       C70I=C70BSMI
       C80I=C80BSMI
       C71I=C71BSMI
       C81I=C81BSMI

*	 2) Wilson coefficients at m_b

      eta=asc0/asmb

*          - LO coefficients [2,7]

      aa(1)=14d0/23d0                  ! Array: aa_i
      aa(2)=16d0/23d0
      aa(3)=6d0/23d0
      aa(4)=-12d0/23d0
      aa(5)=0.4086d0
      aa(6)=-0.4230d0
      aa(7)=-0.8994d0
      aa(8)=0.1456d0

      do i=1,4                         ! Array: bb_i
       j=4+i
       bb(i)=aa(j)
      enddo

      hh(1)=626126d0/272277d0          ! Array: hh_i
      hh(2)=-56281d0/51730d0
      hh(3)=-3d0/7d0
      hh(4)=-1d0/14d0
      hh(5)=-0.6494d0
      hh(6)=-0.0380d0
      hh(7)=-0.0186d0
      hh(8)=-0.0057d0

      h8(1)=-0.9135d0                  ! Array: h8_i
      h8(2)=0.0873d0
      h8(3)=-0.0571d0
      h8(4)=0.0209d0

      C10b=eta**(aa(3))-eta**(aa(4))

      C20b=1d0/3d0*(2d0*eta**(aa(3))+eta**(aa(4)))

      C30b=2d0/63d0*eta**(aa(3))-eta**(aa(4))/27d0-0.0659d0*eta**(aa(5))
     .+0.0595d0*eta**(aa(6))-0.0218d0*eta**(aa(7))+0.0335d0*eta**(aa(8))

      C40b=eta**(aa(3))/21d0+eta**(aa(4))/9d0+0.0237d0*eta**(aa(5))
     . -0.0173d0*eta**(aa(6))-0.1336d0*eta**(aa(7))-0.0316*eta**(aa(8))

      C50b=-eta**(aa(3))/126d0+eta**(aa(4))/108d0+0.0094d0*eta**(aa(5))
     . -0.0100d0*eta**(aa(6))+0.0010d0*eta**(aa(7))-0.0017*eta**(aa(8))

      C60b=-eta**(aa(3))/84d0-eta**(aa(4))/36d0+0.0108d0*eta**(aa(5))
     .+0.0163d0*eta**(aa(6))+0.0103d0*eta**(aa(7))+0.0023d0*eta**(aa(8))

      C70b=0d0
      do i=1,8
       C70b=C70b+hh(i)*eta**(aa(i))
      enddo
      C70b=C70b+eta**(aa(2))*C70
     .   +8d0/3d0*(eta**(aa(1))-eta**(aa(2)))*C80

      C80b=0d0
      do i=1,4
       C80b=C80b+h8(i)*eta**(bb(i))
      enddo
      C80b=C80b+eta**(aa(1))*(C80+313063d0/363036d0)

      C70bI=eta**(aa(2))*C70I
     .   +8d0/3d0*(eta**(aa(1))-eta**(aa(2)))*C80I

      C80bI=eta**(aa(1))*C80I

*          - NLO coefficients [2,8]

      m0011(1)=0d0
      m0011(2)=0d0
      m0011(3)=1d0/3d0
      m0011(4)=2d0/3d0
      m0011(5)=0d0
      m0011(6)=0d0
      m0011(7)=0d0
      m0011(8)=0d0

      m0021(1)=0d0
      m0021(2)=0d0
      m0021(3)=0.2222d0
      m0021(4)=-0.2222d0
      m0021(5)=0d0
      m0021(6)=0d0
      m0021(7)=0d0
      m0021(8)=0d0

      m0031(1)=0d0
      m0031(2)=0d0
      m0031(3)=0.0106d0
      m0031(4)=0.0247d0
      m0031(5)=-0.0129d0
      m0031(6)=-0.0497d0
      m0031(7)=0.0092d0
      m0031(8)=0.0182d0

      m0041(1)=0d0
      m0041(2)=0d0
      m0041(3)=0.0159d0
      m0041(4)=-0.0741d0
      m0041(5)=0.0046d0
      m0041(6)=0.0144d0
      m0041(7)=0.0562d0
      m0041(8)=-0.0171d0

      m0051(1)=0d0
      m0051(2)=0d0
      m0051(3)=-0.0026d0
      m0051(4)=-0.0062d0
      m0051(5)=0.0018d0
      m0051(6)=0.0083d0
      m0051(7)=-0.0004d0
      m0051(8)=-0.0009d0

      m0061(1)=0d0
      m0061(2)=0d0
      m0061(3)=-0.0040d0
      m0061(4)=0.0185d0
      m0061(5)=0.0021d0
      m0061(6)=-0.0136d0
      m0061(7)=-0.0043d0
      m0061(8)=0.0012d0

      m0071(1)=0.5784d0
      m0071(2)=-0.3921d0
      m0071(3)=-0.1429d0
      m0071(4)=0.0476d0
      m0071(5)=-0.1275d0
      m0071(6)=0.0317d0
      m0071(7)=0.0078d0
      m0071(8)=-0.0031d0

      m0081(1)=0.2169d0
      m0081(2)=0d0
      m0081(3)=0d0
      m0081(4)=0d0
      m0081(5)=-0.1793d0
      m0081(6)=-0.0730d0
      m0081(7)=0.0240d0
      m0081(8)=0.0113d0

      m0034(1)=0d0
      m0034(2)=0d0
      m0034(3)=0d0
      m0034(4)=0d0
      m0034(5)=-0.1933d0
      m0034(6)=0.1579d0
      m0034(7)=0.1428d0
      m0034(8)=-0.1074d0

      m0044(1)=0d0
      m0044(2)=0d0
      m0044(3)=0d0
      m0044(4)=0d0
      m0044(5)=0.0695d0
      m0044(6)=-0.0459d0
      m0044(7)=0.8752d0
      m0044(8)=0.1012d0

      m0054(1)=0d0
      m0054(2)=0d0
      m0054(3)=0d0
      m0054(4)=0d0
      m0054(5)=0.0274d0
      m0054(6)=-0.0264d0
      m0054(7)=-0.0064d0
      m0054(8)=0.0055d0

      m0064(1)=0d0
      m0064(2)=0d0
      m0064(3)=0d0
      m0064(4)=0d0
      m0064(5)=0.0317d0
      m0064(6)=0.0432d0
      m0064(7)=-0.0675d0
      m0064(8)=-0.0074d0

      m0074(1)=4661194d0/816831d0
      m0074(2)=-8516d0/2217d0
      m0074(3)=0d0
      m0074(4)=0d0
      m0074(5)=-1.9043d0
      m0074(6)=-0.1008d0
      m0074(7)=0.1216d0
      m0074(8)=0.0183d0

      m0084(1)=2.1399d0
      m0084(2)=0d0
      m0084(3)=0d0
      m0084(4)=0d0
      m0084(5)=-2.6788d0
      m0084(6)=0.2318d0
      m0084(7)=0.3741d0
      m0084(8)=-0.0670d0

      m1012(1)=0d0
      m1012(2)=0d0
      m1012(3)=-2.9606d0
      m1012(4)=-4.0951d0
      m1012(5)=0d0
      m1012(6)=0d0
      m1012(7)=0d0
      m1012(8)=0d0

      m1022(1)=0d0
      m1022(2)=0d0
      m1022(3)=-1.9737d0
      m1022(4)=1.3650d0
      m1022(5)=0d0
      m1022(6)=0d0
      m1022(7)=0d0
      m1022(8)=0d0

      m1032(1)=0d0
      m1032(2)=0d0
      m1032(3)=-0.0940d0
      m1032(4)=-0.1517d0
      m1032(5)=-0.2327d0
      m1032(6)=0.2288d0
      m1032(7)=0.1455d0
      m1032(8)=-0.4760d0

      m1042(1)=0d0
      m1042(2)=0d0
      m1042(3)=-0.1410d0
      m1042(4)=0.4550d0
      m1042(5)=0.0836d0
      m1042(6)=-0.0664d0
      m1042(7)=0.8919d0
      m1042(8)=0.4485d0

      m1052(1)=0d0
      m1052(2)=0d0
      m1052(3)=0.0235d0
      m1052(4)=0.0379d0
      m1052(5)=0.0330d0
      m1052(6)=-0.0383d0
      m1052(7)=-0.0066d0
      m1052(8)=0.0242d0

      m1062(1)=0d0
      m1062(2)=0d0
      m1062(3)=0.0352d0
      m1062(4)=-0.1138d0
      m1062(5)=0.0382d0
      m1062(6)=0.0625d0
      m1062(7)=-0.0688d0
      m1062(8)=-0.0327d0

      m1072(1)=9.9372d0
      m1072(2)=-7.4878d0
      m1072(3)=1.2688d0
      m1072(4)=-0.2925d0
      m1072(5)=-2.2923d0
      m1072(6)=-0.1461d0
      m1072(7)=0.1239d0
      m1072(8)=0.0812d0

      m1082(1)=3.7264d0
      m1082(2)=0d0
      m1082(3)=0d0
      m1082(4)=0d0
      m1082(5)=-3.2247d0
      m1082(6)=0.3359d0
      m1082(7)=0.3812d0
      m1082(8)=-0.2968d0

      m1112(1)=0d0
      m1112(2)=0d0
      m1112(3)=5.9606d0
      m1112(4)=1.0951d0
      m1112(5)=0d0
      m1112(6)=0d0
      m1112(7)=0d0
      m1112(8)=0d0

      m1122(1)=0d0
      m1122(2)=0d0
      m1122(3)=1.9737d0
      m1122(4)=-1.3650d0
      m1122(5)=0d0
      m1122(6)=0d0
      m1122(7)=0d0
      m1122(8)=0d0

      m1132(1)=0d0
      m1132(2)=0d0
      m1132(3)=-0.5409d0
      m1132(4)=1.6332d0
      m1132(5)=1.6406d0
      m1132(6)=-1.6702d0
      m1132(7)=-0.2576d0
      m1132(8)=-0.2250d0

      m1142(1)=0d0
      m1142(2)=0d0
      m1142(3)=2.2203d0
      m1142(4)=2.0265d0
      m1142(5)=-4.1830d0
      m1142(6)=-0.7135d0
      m1142(7)=-1.8215d0
      m1142(8)=0.7996d0

      m1152(1)=0d0
      m1152(2)=0d0
      m1152(3)=0.0400d0
      m1152(4)=-0.1861d0
      m1152(5)=-0.1669d0
      m1152(6)=0.1887d0
      m1152(7)=0.0201d0
      m1152(8)=0.0304d0

      m1162(1)=0d0
      m1162(2)=0d0
      m1162(3)=-0.2614d0
      m1162(4)=-0.1918d0
      m1162(5)=0.4197d0
      m1162(6)=0.0295d0
      m1162(7)=0.1474d0
      m1162(8)=-0.0640d0

      m1172(1)=-17.3023d0
      m1172(2)=8.5027d0
      m1172(3)=4.5508d0
      m1172(4)=0.7519d0
      m1172(5)=2.0040d0
      m1172(6)=0.7476d0
      m1172(7)=-0.5385d0
      m1172(8)=0.0914d0

      m1182(1)=-5.8157d0
      m1182(2)=0d0
      m1182(3)=1.4062d0
      m1182(4)=-3.9895d0
      m1182(5)=3.2850d0
      m1182(6)=3.6851d0
      m1182(7)=-0.1424d0
      m1182(8)=0.6492d0

      C11b=0d0
      C21b=0d0
      C31b=0d0
      C41b=0d0
      C51b=0d0
      C61b=0d0
      C71b=0d0
      C81b=0d0
      C71bI=0d0
      C81bI=0d0
      DO I=1,8
       C11b=C11b+eta**(aa(i))*(eta*m0011(i)*C11
     .                                 +eta*m1012(i)+m1112(i))
       C21b=C21b+eta**(aa(i))*(eta*m0021(i)*C11
     .                                 +eta*m1022(i)+m1122(i))
       C31b=C31b+eta**(aa(i))*(eta*m0031(i)*C11+eta*m0034(i)*C41
     .                                 +eta*m1032(i)+m1132(i))
       C41b=C41b+eta**(aa(i))*(eta*m0041(i)*C11+eta*m0044(i)*C41
     .                                 +eta*m1042(i)+m1142(i))
       C51b=C51b+eta**(aa(i))*(eta*m0051(i)*C11+eta*m0054(i)*C41
     .                                 +eta*m1052(i)+m1152(i))
       C61b=C61b+eta**(aa(i))*(eta*m0061(i)*C11+eta*m0064(i)*C41
     .                                 +eta*m1062(i)+m1162(i))
       C71b=C71b+eta**(aa(i))*(eta*m0071(i)*C11+eta*m0074(i)*C41
     .                                 +eta*m1072(i)+m1172(i))
       C81b=C81b+eta**(aa(i))*(eta*m0081(i)*C11+eta*m0084(i)*C41
     .                                 +eta*m1082(i)+m1182(i))
      ENDDO

      C71b=C71b+eta**(aa(2)+1d0)*C71
     .         +8d0/3d0*eta*(eta**(aa(1))-eta**(aa(2)))*C81
     .         +37208d0/4761d0*eta**(aa(2))*(eta-1d0)*C70
     .         +(eta**(aa(1))*(256868d0/14283d0*eta-7164416d0/357075d0)
     .          -eta**(aa(2))*(6698884d0/357075d0*eta-297664d0/14283d0))
     .                                               *C80

      C81b=C81b+eta**(aa(1))*(eta*C81+6.7441d0*(eta-1d0)*C80)

      C71bI=C71bI+eta**(aa(2)+1d0)*C71I
     .         +8d0/3d0*eta*(eta**(aa(1))-eta**(aa(2)))*C81I
     .         +37208d0/4761d0*eta**(aa(2))*(eta-1d0)*C70I
     .         +(eta**(aa(1))*(256868d0/14283d0*eta-7164416d0/357075d0)
     .          -eta**(aa(2))*(6698884d0/357075d0*eta-297664d0/14283d0))
     .                                               *C80I

      C81bI=C81bI+eta**(aa(1))*(eta*C81I+6.7441d0*(eta-1d0)*C80I)

*          - Electroweak corrections [7,9]

      C7EMb=(88d0/575d0*eta**(16d0/23d0)
     .      -40d0/69d0*eta**(-7d0/23d0)
     .      +32d0/75d0*eta**(-9d0/23d0))
     .   *(C70BSM+asc0/(4d0*pi)*C71BSM)
     .     +(640d0/1449d0*eta**(14d0/23d0)
     .     -704d0/1725d0*eta**(16d0/23d0)
     .     +32d0/1449d0*eta**(-7d0/23d0)
     .     -32d0/575d0*eta**(-9d0/23d0))
     .   *(C80BSM+asc0/(4d0*pi)*C81BSM)
      C7EMbI=(88d0/575d0*eta**(16d0/23d0)
     .      -40d0/69d0*eta**(-7d0/23d0)
     .      +32d0/75d0*eta**(-9d0/23d0))
     .   *(C70BSMI+asc0/(4d0*pi)*C71BSMI)
     .     +(640d0/1449d0*eta**(14d0/23d0)
     .     -704d0/1725d0*eta**(16d0/23d0)
     .     +32d0/1449d0*eta**(-7d0/23d0)
     .     -32d0/575d0*eta**(-9d0/23d0))
     .   *(C80BSMI+asc0/(4d0*pi)*C81BSMI)
*      Summation of the Electroweak Corrections, the SM Contribution
*      is set to 0.0071.
      EPSew=0.0071d0
     .     +ALEMMZ*(C7EMb/asmb-(eta**(aa(2))*C70BSM
     .   +8d0/3d0*(eta**(aa(1))-eta**(aa(2)))*C80BSM)*dlog(MZ/scb)/pi)
      EPSewI=ALEMMZ*(C7EMbI/asmb-(eta**(aa(2))*C70BSMI
     .   +8d0/3d0*(eta**(aa(1))-eta**(aa(2)))*C80BSMI)*dlog(MZ/scb)/pi)

*   BSM error estimate from DC7BSM and DC8BSM:
      DC7BSM_b=eta**(aa(2))*DABS(DC7BSM)
     .   +8d0/3d0*DABS(eta**(aa(1))-eta**(aa(2)))*DABS(DC8BSM)
      DC8BSM_b=eta**(aa(1))*DABS(DC8BSM)
      DC7BSM_bI=eta**(aa(2))*DABS(DC7BSMI)
     .   +8d0/3d0*DABS(eta**(aa(1))-eta**(aa(2)))*DABS(DC8BSMI)
      DC8BSM_bI=eta**(aa(1))*DABS(DC8BSMI)

*	 3) Contributions to B --> Xs gamma

*         - Bremsstrahlung Contribution with delta=0.295
*          Fit from gaussian integration of 'G' [7] (precision ~ few permil)
*          CKM-dependence from [10]
      MBp=4.539d0                             ! [11]
      MB_kin=4.564d0                          ! in the 'kinetic scheme'
      MC_scb=1.087d0
      z=(MC_scb/MB_kin)**2                    ! (MC/MB)^2, values from [11], ap. D
      delt=1d0-2d0*1.6d0/MBp
      lndelt=dlog(delt)
      lndeltp=dlog(1d0-delt)
      delt2=delt**2

      ff22=0.0019524d0+0.030238d0*z-0.43189d0*z**2+4.3025d0*z**3
     .    -24.151d0*z**4+57.057d0*z**5+delt*(0.1024d0+0.26005d0*z
     .    -23.169d0*z**2+266.33d0*z**3-1530.0d0*z**4+3664.1d0*z**5)
      ff11=1d0/36d0*ff22
      ff12=-1d0/3d0*ff22

      ff27=-0.020604d0+0.017213d0*z+2.9652d0*z**2-17.755d0*z**3
     .    +47.042d0*z**4+delt*(-0.13382d0+2.0837d0*z+10.895*z**2
     .    -172.57d0*z**3+459.87d0*z**4)
      ff27i=0.0020255d0+0.48376d0*z-3.9246d0*z**2+18.816d0*z**3
     .     -55.611d0*z**4+delt*(0.029659d0+4.0743d0*z-59.714d0*z**2
     .     +310.26d0*z**3-673.18d0*z**4)
      aux=ff27
      ff27=(1d0+VVtu)*ff27-VVtuim*ff27i
      ff27i=aux*VVtuim+ff27i*VVtu
      ff17=-1d0/6d0*ff27
      ff17i=-1d0/6d0*ff27i
      ff28=-1d0/3d0*ff27
      ff28i=-1d0/3d0*ff27i
      ff18=-1d0/6d0*ff28
      ff18i=-1d0/6d0*ff28i

      ff47=Pi/54d0*(3d0*dsqrt(3d0)-Pi)+delt**3/81d0-25d0/108d0*delt2 ! [11], ap. C
     . +5d0/54d0*delt+2d0/9d0*(delt2+2d0*delt+3d0)
     .                *datan(dsqrt((1d0-delt)/(3d0+delt)))**2
     . -(delt2+4d0*delt+3d0)/3d0*dsqrt((1d0-delt)/(3d0+delt))
     .                *datan(dsqrt((1d0-delt)/(3d0+delt)))
     . +(34d0*delt2+59d0*delt-18d0)*delt2*lndelt/486d0/(1d0-delt)
     . +(433d0*delt**3+429d0*delt**2-720d0*delt)/2916d0
!     -delt/54d0*(1d0-delt+delt**2/3d0)          ! [9]
!     .     +(0.00131501d0+ 0.00420692d0*delt)/12d0
      ff48=-ff47/3d0

      ff77=1d0/3d0*(10d0*delt+delt2-2d0/3d0*delt**3
     .     +delt*(delt-4d0)*lndelt)
     .   -1d0/3d0*(2d0*lndelt**2+7d0*lndelt+31d0/3d0)
      ff78=8d0/9d0*(sp2(1d0-delt)-pi**2/6d0-delt*lndelt
     .      +9d0/4d0*delt-delt2/4d0+delt**3/12d0)
      ff88=1d0/27d0*(4d0*sp2(1d0-delt)-2d0/3d0*pi**2
     .       +8d0*lndeltp-delt*(2d0+delt)*lndelt
     .       +7d0*delt+3d0*delt2-2d0/3d0*delt**3
     .       -2d0*(2d0*delt+delt2+4d0*lndeltp)*dlog(MB/0.95d0))
     
*         - LO 4-body decay
      T1=23d0/16d0*delt**4-delt**4/2d0*lndelt-191d0/108d0*delt**3
     . +4d0/9d0*delt**3*lndelt+17d0/18d0*delt**2-delt**2/3d0*lndelt
     . +109d0/18d0*delt-5d0/3d0*delt*lndelt+79d0/18d0*lndeltp
     . -5d0/3d0*lndeltp*lndelt
     . -(delt**4-8d0/9d0*delt**3+2d0/3d0*delt**2+10d0/3d0*delt
     .  +10d0/3d0*lndeltp)*dlog(25d0)-5d0/3d0*Sp2(delt)

      T2=1181d0/2592d0*delt**4-17d0/108d0*delt**4*lndelt
     . -395d0/648d0*delt**3+4d0/27d0*delt**3*lndelt
     . +7d0/18d0*delt**2-delt**2/9d0*lndelt+187d0/108d0*delt
     . -delt*lndelt/2d0+133d0/108d0*lndeltp-lndeltp*lndelt/2d0
     . -(17d0/54d0*delt**4-8d0/27d0*delt**3+2d0/9d0*delt**2+delt
     .  +lndeltp)*dlog(25d0)-Sp2(delt)/2d0

      T3=341d0/7776d0*delt**4-5d0*delt**4/324d0*lndelt
     . -89d0/1944d0*delt**3+delt**3*lndelt/81d0+delt**2/72d0
     . -delt**2*lndelt/108d0+35d0/162d0*delt-delt*lndelt/18d0
     . +13d0/81d0*lndeltp-lndeltp*lndelt/18d0
     . -(5d0/162d0*delt**4-2d0/81d0*delt**3+delt**2/54d0+delt/9d0
     .  +lndeltp/9d0)*dlog(25d0)-Sp2(delt)/18d0

      LO4B=T1*((C30b+asmb/4d0/Pi*C31b)**2
     .    +2d0/9d0*(C40b+asmb/4d0/Pi*C41b)**2
     .    +136d0*(C50b+asmb/4d0/Pi*C51b)**2
     .    +272d0/9d0*(C60b+asmb/4d0/Pi*C61b)**2
     .    +20d0*(C30b+asmb/4d0/Pi*C31b)*(C50b+asmb/4d0/Pi*C51b)
     .    +40d0/9d0*(C40b+asmb/4d0/Pi*C41b)*(C60b+asmb/4d0/Pi*C61b))
     .   +T3*((C30b+asmb/4d0/Pi*C31b)**2
     .    -2d0/9d0*(C40b+asmb/4d0/Pi*C41b)**2
     .    +256d0*(C50b+asmb/4d0/Pi*C51b)**2
     .    -512d0/9d0*(C60b+asmb/4d0/Pi*C61b)**2
     .    +8d0/3d0*(C30b+asmb/4d0/Pi*C31b)*(C40b+asmb/4d0/Pi*C41b)
     .    +32d0*(C30b+asmb/4d0/Pi*C31b)*(C50b+asmb/4d0/Pi*C51b)
     .    +128d0/3d0*(C30b+asmb/4d0/Pi*C31b)*(C60b+asmb/4d0/Pi*C61b)
     .    +128d0/3d0*(C40b+asmb/4d0/Pi*C41b)*(C50b+asmb/4d0/Pi*C51b)
     .    -64d0/9d0*(C40b+asmb/4d0/Pi*C41b)*(C60b+asmb/4d0/Pi*C61b)
     .    +2048d0/3d0*(C50b+asmb/4d0/Pi*C51b)*(C60b+asmb/4d0/Pi*C61b))
     .   +T2*(VVtu**2-VVtuim**2)*(2d0/9d0*(C10b+asmb/4d0/Pi*C11b)**2
     .                             +(C20b+asmb/4d0/Pi*C21b)**2)
     .   -2d0*T2*VVtu*(C10b+asmb/4d0/Pi*C11b)
     .         *(4d0/9d0*(C30b+asmb/4d0/Pi*C31b)
     .          -2d0/27d0*(C40b+asmb/4d0/Pi*C41b)
     .          +64d0/9d0*(C50b+asmb/4d0/Pi*C51b)
     .          -32d0/27d0*(C60b+asmb/4d0/Pi*C61b))
     .   -2d0*T2*VVtu*(C20b+asmb/4d0/Pi*C21b)
     .         *((C30b+asmb/4d0/Pi*C31b)/3d0
     .          +4d0/9d0*(C40b+asmb/4d0/Pi*C41b)
     .          +16d0/3d0*(C50b+asmb/4d0/Pi*C51b)
     .          +64d0/9d0*(C60b+asmb/4d0/Pi*C61b))
     
*	  - NLO 4-body decay
      NLO4B=-0.00165806d0*(VVtu**2-VVtuim**2)
     .                                      *(C10b+asmb/4d0/Pi*C11b)**2
     .      +2d0*0.00994835d0*(VVtu**2-VVtuim**2)
     .                 *(C10b+asmb/4d0/Pi*C11b)*(C20b+asmb/4d0/Pi*C21b)
     .      +2d0*(0.00841465d0*VVtu+0.00607241d0*VVtuim)
     .                 *(C10b+asmb/4d0/Pi*C11b)*(C30b+asmb/4d0/Pi*C31b)
     .      -2d0*(0.0138997d0*VVtu-0.0157189d0*VVtuim)
     .                 *(C10b+asmb/4d0/Pi*C11b)*(C40b+asmb/4d0/Pi*C41b)
     .      +2d0*(0.151345d0*VVtu+0.0971586d0*VVtuim)
     .                 *(C10b+asmb/4d0/Pi*C11b)*(C50b+asmb/4d0/Pi*C51b)
     .      -2d0*(0.0740271d0*VVtu-0.136086d0*VVtuim)
     .                 *(C10b+asmb/4d0/Pi*C11b)*(C60b+asmb/4d0/Pi*C61b)
     .      -2d0*(0.0000758914d0*VVtu-0.00034605d0*VVtuim)*VVtc
     .                 *(C10b+asmb/4d0/Pi*C11b)**2
     .      +2d0*(0.000455348d0*VVtu-0.0020763*VVtuim)*VVtc
     .                 *(C10b+asmb/4d0/Pi*C11b)*(C20b+asmb/4d0/Pi*C21b)
     .      -2d0*VVtu*0.0339475d0
     .                 *(C10b+asmb/4d0/Pi*C11b)*(C70b+asmb/4d0/Pi*C71b)
     .      +2d0*VVtuim*0.0339475d0
     .               *(C10b+asmb/4d0/Pi*C11b)*(C70bi+asmb/4d0/Pi*C71bi)
     .      -2d0*VVtu*0.0136471d0
     .                 *(C10b+asmb/4d0/Pi*C11b)*(C80b+asmb/4d0/Pi*C81b)
     .      +2d0*VVtuim*0.0136471d0
     .               *(C10b+asmb/4d0/Pi*C11b)*(C80bI+asmb/4d0/Pi*C81bI)
     .      -0.0596901d0*(VVtu**2-VVtuim**2)*(C20b+asmb/4d0/Pi*C21b)**2
     .      -2d0*(0.0504879d0*VVtu+0.0364345d0*VVtuim)
     .                 *(C20b+asmb/4d0/Pi*C21b)*(C30b+asmb/4d0/Pi*C31b)
     .      +2d0*(0.0833984d0*VVtu-0.0943134*VVtuim)*VVtc
     .                 *(C20b+asmb/4d0/Pi*C21b)*(C40b+asmb/4d0/Pi*C41b)
     .      -2d0*(0.908071d0*VVtu+0.582951d0*VVtuim)
     .                 *(C20b+asmb/4d0/Pi*C21b)*(C50b+asmb/4d0/Pi*C51b)
     .      +2d0*(0.444163d0*VVtu-0.816517d0*VVtuim)
     .                 *(C20b+asmb/4d0/Pi*C21b)*(C60b+asmb/4d0/Pi*C61b)
     .      +2d0*(0.000455348d0*VVtu-0.0020763d0*VVtuim)*VVtc
     .                 *(C20b+asmb/4d0/Pi*C21b)*(C10b+asmb/4d0/Pi*C11b)
     .      -2d0*(0.00273209d0*VVtu-0.0124578d0*VVtuim)*VVtc
     .                 *(C20b+asmb/4d0/Pi*C21b)**2
     .      +2d0*0.203685d0*VVtu
     .                 *(C20b+asmb/4d0/Pi*C21b)*(C70b+asmb/4d0/Pi*C71b)
     .      -2d0*0.203685d0*VVtuim
     .               *(C20b+asmb/4d0/Pi*C21b)*(C70bI+asmb/4d0/Pi*C71bI)
     .      +2d0*0.0818824d0*VVtu
     .                 *(C20b+asmb/4d0/Pi*C21b)*(C80b+asmb/4d0/Pi*C81b)
     .      -2d0*0.0818824d0*VVtuim
     .               *(C20b+asmb/4d0/Pi*C21b)*(C80bI+asmb/4d0/Pi*C81bI)
     .      -0.00165865d0*(C30b+asmb/4d0/Pi*C31b)**2
     .      +2d0*0.0273824d0
     .                 *(C30b+asmb/4d0/Pi*C31b)*(C40b+asmb/4d0/Pi*C41b)
     .      +2d0*0.00177663d0
     .                 *(C30b+asmb/4d0/Pi*C31b)*(C50b+asmb/4d0/Pi*C51b)
     .      +2d0*0.328096d0
     .                 *(C30b+asmb/4d0/Pi*C31b)*(C60b+asmb/4d0/Pi*C61b)
     .      -2d0*VVtc*0.000877217d0
     .                 *(C30b+asmb/4d0/Pi*C31b)*(C10b+asmb/4d0/Pi*C11b)
     .      +2d0*VVtc*0.0052633d0
     .                 *(C30b+asmb/4d0/Pi*C31b)*(C20b+asmb/4d0/Pi*C21b)
     .      +2d0*0.0489658d0
     .                 *(C30b+asmb/4d0/Pi*C31b)*(C70b+asmb/4d0/Pi*C71b)
     .      -2d0*0.0144535d0
     .                 *(C30b+asmb/4d0/Pi*C31b)*(C80b+asmb/4d0/Pi*C81b)
     .      -0.186252d0*(C40b+asmb/4d0/Pi*C41b)**2
     .      +2d0*0.590294d0
     .                 *(C40b+asmb/4d0/Pi*C41b)*(C50b+asmb/4d0/Pi*C51b)
     .      -2d0*1.35313d0
     .                 *(C40b+asmb/4d0/Pi*C41b)*(C60b+asmb/4d0/Pi*C61b)
     .      -2d0*VVtc*0.00206358d0
     .                 *(C40b+asmb/4d0/Pi*C41b)*(C10b+asmb/4d0/Pi*C11b)
     .      +2d0*VVtc*0.0123815d0
     .                 *(C40b+asmb/4d0/Pi*C41b)*(C20b+asmb/4d0/Pi*C21b)
     .      -2d0*0.113915d0
     .                 *(C40b+asmb/4d0/Pi*C41b)*(C70b+asmb/4d0/Pi*C71b)
     .      -2d0*0.108381d0
     .                 *(C40b+asmb/4d0/Pi*C41b)*(C80b+asmb/4d0/Pi*C81b)
     .      +0.481466d0*(C50b+asmb/4d0/Pi*C51b)**2
     .      +2d0*6.74297d0
     .                 *(C50b+asmb/4d0/Pi*C51b)*(C60b+asmb/4d0/Pi*C61b)
     .      -2d0*0.0140355d0*VVtc
     .                 *(C50b+asmb/4d0/Pi*C51b)*(C10b+asmb/4d0/Pi*C11b)
     .      +2d0*0.0842129d0*VVtc
     .                 *(C50b+asmb/4d0/Pi*C51b)*(C20b+asmb/4d0/Pi*C21b)
     .      +2d0*0.783453d0
     .                 *(C50b+asmb/4d0/Pi*C51b)*(C70b+asmb/4d0/Pi*C71b)
     .      -2d0*0.231256d0
     .                 *(C50b+asmb/4d0/Pi*C51b)*(C80b+asmb/4d0/Pi*C81b)
     .      -8.54439d0*(C60b+asmb/4d0/Pi*C61b)**2
     .      -2d0*0.0197586d0*VVtc
     .                 *(C60b+asmb/4d0/Pi*C61b)*(C10b+asmb/4d0/Pi*C11b)
     .      +2d0*0.118551d0*VVtc
     .                 *(C60b+asmb/4d0/Pi*C61b)*(C20b+asmb/4d0/Pi*C21b)
     .      -2d0*1.18811d0
     .                 *(C60b+asmb/4d0/Pi*C61b)*(C70b+asmb/4d0/Pi*C71b)
     .      -2d0*1.06935d0
     .                 *(C60b+asmb/4d0/Pi*C61b)*(C80b+asmb/4d0/Pi*C81b)

      LO4B=(LO4B+asmb/4d0/Pi*NLO4B)
     .       /(1d0+(50d0/3d0-8d0*Pi**2/3d0)*asmb/4d0/Pi)

*         - NLO perturbative coefficients [12,13]
      K17=833d0/729d0-(af(z)+bf(z))/3d0+208d0/243d0*dlog(scb/mb)
     . +2d0*ff17
      K17I=-40d0/243d0*Pi+(afim(z)+bfim(z))/3d0+2d0*ff17i

      K27=-6d0*K17
      K27I=-6d0*K17I

      K37=2392d0/243d0+8d0*Pi/3d0/dsqrt(3d0)+32d0/9d0*(-0.1684d0)
     . -4.0859d0+2d0*0.0316d0+176d0/81d0*dlog(scb/mb)
      K37I=(4d0/9d0-2d0*4d0/81d0)*Pi-56d0/81d0*Pi

      K47=-761d0/729d0-4d0*Pi/9d0/dsqrt(3d0)-16d0/27d0*(-0.1684d0)
     . +4.0859d0/6d0+5d0/3d0*0.0316d0+2d0*bf(z)
     . +152d0/243d0*dlog(scb/mb)+2d0*ff47
      K47I=-(4d0/9d0/6d0+5d0/3d0*4d0/81d0)*Pi+148d0/243d0*Pi
     . -2d0*bfim(z)

      K57=56680d0/243d0+32d0*Pi/3d0/dsqrt(3d0)+128d0/9d0*(-0.1684d0)
     . -16d0*4.0859d0+32d0*0.0316d0+6272d0/81d0*dlog(scb/mb)
      K57I=(16d0*4d0/9d0-32d0*4d0/81d0)*Pi-896d0/81d0*Pi

      K67=5710d0/729d0-16d0*Pi/9d0/dsqrt(3d0)-64d0/27d0*(-0.1684d0)
     . -10d0/3d0*4.0859d0+44d0/3d0*0.0316d0+12d0*af(z)+20d0*bf(z)
     . -4624d0/243d0*dlog(scb/mb)
      K67I=(10d0/3d0*4d0/9d0-44d0/3d0*4d0/81d0)*Pi+2296d0/243d0*Pi
     . -12d0*afim(z)-20d0*bfim(z)

      K77=-91d0/9d0+4d0/9d0*Pi**2-32d0/3d0*dlog(scb/mb)+2d0*ff77

      K78=44d0/9d0-8d0/27d0*Pi**2+32d0/9d0*dlog(scb/mb)+2d0*ff78
      K78I=8d0/9d0*Pi

      BSGPERT=C70b**2+C70bI**2                          ! NNLO expansion [12]
     . +asmb/2d0/Pi*(C70b*C71b+C70bI*C71bI
     .              +ff11*C10b**2+ff12*C10b*C20b+ff22*C20b**2
     .              +K17*C10b*C70b-K17I*C10b*C70bI
     .              +K27*C20b*C70b-K27I*C20b*C70bI
     .              +K37*C30b*C70b-K37I*C30b*C70bI
     .              +K47*C40b*C70b-K47I*C40b*C70bI
     .              +K57*C50b*C70b-K57I*C50b*C70bI
     .              +K67*C60b*C70b-K67I*C60b*C70bI
     .              +K77*(C70b**2+C70bI**2)+K78*(C70b*C80b+C70bI*C80bI)
     .               +K78I*(C80b*C70bI-C70b*C80bI)
     .              +ff18*C10b*C80b-ff18I*C10b*C80bI
     .              +ff28*C20b*C80b-ff18I*C20b*C80bI
     .              +ff48*C40b*C80b
     .              +ff88*(C80b**2+C80bI**2))
     . +(asmb/4d0/Pi)**2*((C71b**2+C71bI**2)
     .   +2d0*(2d0*ff11*C10b*C11b+ff12*(C10b*C21b+C11b*C20b)
     .        +2d0*ff22*C20b*C21b
     .        +K17*(C10b*C71b+C11b*C70b)-K17I*(C10b*C71bI+C11b*C70bI)
     .        +K27*(C20b*C71b+C21b*C70b)-K27I*(C20b*C71bI+C21b*C70bI)
     .        +K37*(C30b*C71b+C31b*C70b)-K37I*(C30b*C71bI+C31b*C70bI)
     .        +K47*(C40b*C71b+C41b*C70b)-K47I*(C40b*C71bI+C41b*C70bI)
     .        +K57*(C50b*C71b+C51b*C70b)-K57I*(C50b*C71bI+C51b*C70bI)
     .        +K67*(C60b*C71b+C61b*C70b)-K67I*(C60b*C71bI+C61b*C70bI)
     .        +2d0*K77*(C70b*C71b+C70bI*C71bI)
     .        +K78*(C70b*C81b+C71b*C80b)
     .        +K78I*(C80b*C71bI+C81b*C70bI-C80bI*C71b-C81bI*C70b)
     .        +ff18*(C10b*C81b+C11b*C80b)-ff18I*(C10b*C81bI+C11b*C80bI)
     .        +ff28*(C20b*C81b+C21b*C80b)-ff28I*(C20b*C81bI+C21b*C80bI)
     .        +ff48*(C40b*C81b+C41b*C80b)
     .        +2d0*ff88*(C80b*C81b+C80bI*C81bI)))
     . +2d0*(EPSEW*C70b+EPSEWI*C70bI)+0d-4/aux
     . +LO4B
     . +6.256900288d-3 !6.382d-3                       ! Mimicking the missing SM NNLO
     . -1.266116316d-3*C70BSM+1.419489589d-5*C80BSM    ! Mimicking the missing NNLO[SM] BSM
                                                       !  coefficients (from [14])

*         - Non-perturbative effects: Heavy Quark Effective Theory Corrections
*                                                     [15-17]
      lambd1=-0.47d0      ! Numerical input from fit in [11], ap. D
      lambd2=0.309d0/3d0  ! = 1/4(M_B*^2-M_B^2)
      rho1=0.171d0
      rho2=-0.135d0/3d0
      MCNP=1.131d0

      NP17=-(lambd2+(rho1/3d0-13d0*rho2/4d0)/mb)/9d0/MCNP**2 !-2d0/9d0*delt*asmb/Pi

      NP78=0d0 !10d0/9d0*delt*asmb/Pi

      NP88=0d0 !delt*asmb/27d0/Pi*(2d0*dlog(delt*50d0**2)-5d0)

      NP77=(lambd1-9d0*lambd2-(11d0*rho1-27d0*rho2)/3d0/mb)/2d0/mb**2
!     .    +asmb/4d0/Pi*((16d0/9d0*(4d0-Pi**2-3d0*dlog(scb/mb))
!     .         -8d0/3d0*dlog(delt)**2)*(1d0+lambd1/2d0/mb**2)
!     .         -4d0/3d0*(7d0+4d0*delt-delt**2)*dlog(delt)
!     .         -4d0/9d0*(1d0-delt)*(31d0+delt-2d0*delt**2)
!     .      -2d0/9d0*lambd1/mb**2/delt**2*(
!     .         (1d0-delt)*(3d0-14d0*delt+46d0*delt**2-5d0*delt**3)
!     .         +(4d0-12d0*delt+27d0*delt**2+14*delt**3-3d0*delt**4)
!     .                                       *dlog(delt))
!     .      -lambd2/mb**2*((87d0+32d0*Pi**2+27d0*dlog(scb/mb))/9d0
!     .         +16d0/3d0*dlog(delt)**2-(24d0+39d0*delt+92d0*delt**2
!     .                            +7d0*delt**3)*dlog(delt)/3d0/delt
!     .         -(1d0-delt)/delt*(20d0+19d0*delt+15d0*delt**2)))

      HQET=NP17*(C20b-1d0/6d0*C10b)*C70b+NP78*(C70b*C80b+C70bI*C80bI)
     .     +NP77*(C70b**2+C70bI**2)+NP88*(C80b**2+C80bI**2)

*         - The Branching Ratio B -> Xs gamma
      BRSL=0.1067d0            ! [11], ap. D           ! BR[B -> Xc l nub]
      CCSL=(1d0-8d0*z+8*z**3-z**4-12d0*z**2*dlog(z))   ! |V_ub/V_cb|^2 
     .     *(0.903d0-0.588d0*(asf(4.6d0)-0.22d0)       ! *Gam[c e nub]/Gam[u e nub]
     .      +0.0650d0*(4.564d0-4.55d0)-0.1080d0*(1.087d0-1.05)
     .      -0.0122*0.309d0-0.199d0*0.171d0+0.004*(-0.135d0))
      Vbsg=0.9626d0                                    ! (V_ts.V_tb/V_cb)^2
      aux=6d0*ALEM0/(pi*CCSL)*Vbsg*BRSL                !prefactor
      BRSG=aux*(BSGPERT+HQET)

*         - New-Physics coefficients (linearized)
      KC7BSM=(2d0*C70b+asmb/2d0/Pi*(C71b+K17*C10b             ! d/dC70b
     . +K27*C20b+K37*C30b+K47*C40b+K57*C50b+K67*C60b+2d0*K77*C70b
     . +K78*C80b-K78I*C80bI)
     . +2d0*(asmb/4d0/Pi)**2*(K17*C11b+K27*C21b+K37*C31b
     . +K47*C41b+K57*C51b+K67*C61b+2d0*K77*C71b+K78*C81b-K78I*C81bI)
     . +2d0*EPSEW
     . +NP17*(C20b-1d0/6d0*C10b)+NP78*C80b+2d0*NP77*C70b
     . +asmb/4d0/Pi/(1d0+(50d0/3d0-8d0*Pi**2/3d0)*asmb/4d0/Pi)
     .  *(-2d0*VVtu*0.0339475d0*(C10b+asmb/4d0/Pi*C11b)
     .      +2d0*0.203685d0*VVtu*(C20b+asmb/4d0/Pi*C21b)
     .      +2d0*0.0489658d0*(C30b+asmb/4d0/Pi*C31b)
     .      -2d0*0.113915d0*(C40b+asmb/4d0/Pi*C41b)
     .      +2d0*0.783453d0*(C50b+asmb/4d0/Pi*C51b)
     .      -2d0*1.18811d0*(C60b+asmb/4d0/Pi*C61b)))
     .                                     *eta**(aa(2))      ! dC70b/dC70BSM
     . +(asmb/2d0/Pi*C70b+(asmb/4d0/Pi)**2*(2d0*C71b          ! d/dC71b
     . +K17*C10b+K27*C20b+K37*C30b+K47*C40b+K57*C50b+K67*C60b
     . +2d0*K77*C70b+K78*C80b-K78I*C80bI))
     .            *37208d0/4761d0*eta**(aa(2))*(eta-1d0)      ! dC71b/dC70BSM
     . +2d0*C70b                                              ! d/dEPSew
     .         *((88d0/575d0*eta**(16d0/23d0)                 ! dEPSew/dC70BSM 
     .            -40d0/69d0*eta**(-7d0/23d0)
     .            +32d0/75d0*eta**(-9d0/23d0))/asmb
     .            -dlog(MZ/scb)/pi*eta**(aa(2)))*ALEMMZ
     . -1.266116316d-3                                        ! NNLO

      KC7BSMI=(2d0*C70bI+asmb/2d0/Pi*(C71bI-K17I*C10b         ! d/dC70bI
     . -K27I*C20b-K37I*C30b-K47I*C40b-K57I*C50b-K67I*C60b
     . +2d0*K77*C70bI+K78*C80bI+K78I*C80b)
     . +2d0*(asmb/4d0/Pi)**2*(-K17I*C11b-K27I*C21b-K37I*C31b
     . -K47I*C41b-K57I*C51b-K67I*C61b+2d0*K77*C71bI+K78*C81bI
     . +K78I*C81b)+2d0*EPSEWI
     . +NP78*C80bI+2d0*NP77*C70bI
     . +asmb/4d0/Pi/(1d0+(50d0/3d0-8d0*Pi**2/3d0)*asmb/4d0/Pi)
     .  *(2d0*VVtuim*0.0339475d0*(C10b+asmb/4d0/Pi*C11b)
     .      -2d0*0.203685d0*VVtuim*(C20b+asmb/4d0/Pi*C21b)))
     .                                     *eta**(aa(2))      ! dC70bI/dC70BSMI
     . +(asmb/2d0/Pi*C70bI+(asmb/4d0/Pi)**2*(2d0*C71bI        ! d/dC71b
     . -K17I*C10b-K27I*C20b-K37I*C30b-K47I*C40b-K57I*C50b-K67I*C60b
     . +2d0*K77*C70bI+K78*C80bI+K78I*C80b))
     .            *37208d0/4761d0*eta**(aa(2))*(eta-1d0)      ! dC71b/dC70BSM
     . +2d0*C70bI                                             ! d/dEPSewI
     .         *((88d0/575d0*eta**(16d0/23d0)                 ! dEPSewI/dC70BSMI 
     .            -40d0/69d0*eta**(-7d0/23d0)
     .            +32d0/75d0*eta**(-9d0/23d0))/asmb
     .            -dlog(MZ/scb)/pi*eta**(aa(2)))*ALEMMZ

      KC8BSM=(2d0*C70b+asmb/2d0/Pi*(C71b+K17*C10b             ! d/dC70b
     . +K27*C20b+K37*C30b+K47*C40b+K57*C50b+K67*C60b+2d0*K77*C70b
     . +K78*C80b)+2d0*(asmb/4d0/Pi)**2*(K17*C11b+K27*C21b+K37*C31b
     . +K47*C41b+K57*C51b+K67*C61b+2d0*K77*C71b+K78*C81b)+2d0*EPSEW
     . +NP17*(C20b-1d0/6d0*C10b)+NP78*C80b+2d0*NP77*C70b
     . +asmb/4d0/Pi/(1d0+(50d0/3d0-8d0*Pi**2/3d0)*asmb/4d0/Pi)
     .  *(-2d0*VVtu*0.0339475d0*(C10b+asmb/4d0/Pi*C11b)
     .      +2d0*0.203685d0*VVtu*(C20b+asmb/4d0/Pi*C21b)
     .      +2d0*0.0489658d0*(C30b+asmb/4d0/Pi*C31b)
     .      -2d0*0.113915d0*(C40b+asmb/4d0/Pi*C41b)
     .      +2d0*0.783453d0*(C50b+asmb/4d0/Pi*C51b)
     .      -2d0*1.18811d0*(C60b+asmb/4d0/Pi*C61b)))
     .              *8d0/3d0*(eta**(aa(1))-eta**(aa(2)))      ! dC70b/dC80BSM
     . +(asmb/2d0/Pi*C70b+(asmb/4d0/Pi)**2*(2d0*C71b          ! d/dC71b
     . +K17*C10b+K27*C20b+K37*C30b+K47*C40b+K57*C50b+K67*C60b
     . +2d0*K77*C70b+K78*C80b))
     . *(eta**(aa(1))*(256868d0/14283d0*eta-7164416d0/357075d0) ! dC71b/dC80BSM
     .         -eta**(aa(2))*(6698884d0/357075d0*eta-297664d0/14283d0))
     . +2d0*C70b                                              ! d/dEPSew
     .   *ALEMMZ*((640d0/1449d0*eta**(14d0/23d0)              ! dEPSew/dC80BSM 
     .           -704d0/1725d0*eta**(16d0/23d0)
     .           +32d0/1449d0*eta**(-7d0/23d0)
     .           -32d0/575d0*eta**(-9d0/23d0))/asmb
     .    -dlog(MZ/scb)/pi*8d0/3d0*(ETA**(14d0/23d0)-ETA**(12d0/23d0)))
     . +(asmb/2d0/Pi*(K78*C70b+ff18*C10b+ff28*C20b            ! d/dC80b
     . +ff48*C40b+2d0*ff88*C80b)+(asmb/4d0/Pi)**2*(2d0*K78*C71b
     . +ff18*C11b+ff28*C21b+ff48*C41b+2d0*ff88*C81b)
     . +NP78*C70b+2d0*NP88*C80b
     . +asmb/4d0/Pi/(1d0+(50d0/3d0-8d0*Pi**2/3d0)*asmb/4d0/Pi)
     .  *(-2d0*VVtu*0.0136471d0*(C10b+asmb/4d0/Pi*C11b)
     .      +2d0*0.0818824d0*VVtu*(C20b+asmb/4d0/Pi*C21b)
     .      -2d0*0.0144535d0*(C30b+asmb/4d0/Pi*C31b)
     .      -2d0*0.108381d0*(C40b+asmb/4d0/Pi*C41b)
     .      -2d0*0.231256d0*(C50b+asmb/4d0/Pi*C51b)
     .      -2d0*1.06935d0*(C60b+asmb/4d0/Pi*C61b)))
     .                           *eta**(aa(1))                ! dC80b/dC80BSM
     . +((asmb/4d0/Pi)**2*(K78*C70b+ff18*C10b+ff28*C20b       ! d/dC81b
     . +ff48*C40b+2d0*ff88*C80b))
     .                  *eta**(aa(1))*6.7441d0*(eta-1d0)      ! dC81b/dC80BSM
     . +1.419489589d-5                                        ! NNLO

      KC8BSMI=(2d0*C70bI+asmb/2d0/Pi*(C71bI-K17I*C10b         ! d/dC70bI
     . -K27I*C20b-K37I*C30b-K47I*C40b-K57I*C50b-K67I*C60b+2d0*K77*C70bI
     . +K78*C80bI+K78I*C80b)+2d0*(asmb/4d0/Pi)**2*(-K17I*C11b-K27I*C21b
     . -K37I*C31b-K47I*C41b-K57I*C51b-K67I*C61b+2d0*K77*C71bI+K78I*C81b
     . +K78*C81bI)+2d0*EPSEWI+NP78*C80bI+2d0*NP77*C70bI
     . +asmb/4d0/Pi/(1d0+(50d0/3d0-8d0*Pi**2/3d0)*asmb/4d0/Pi)
     .  *(2d0*VVtuim*0.0339475d0*(C10b+asmb/4d0/Pi*C11b)
     .      -2d0*0.203685d0*VVtuim*(C20b+asmb/4d0/Pi*C21b)))
     .              *8d0/3d0*(eta**(aa(1))-eta**(aa(2)))      ! dC70bI/dC80BSMI
     . +(asmb/2d0/Pi*C70b+(asmb/4d0/Pi)**2*(2d0*C71bI         ! d/dC71bI
     . -K17I*C10b-K27I*C20b-K37I*C30b-K47I*C40b-K57I*C50b-K67I*C60b
     . +2d0*K77*C70bI+K78I*C80b+K78*C80bI))
     . *(eta**(aa(1))*(256868d0/14283d0*eta-7164416d0/357075d0) ! dC71bI/dC80BSMI
     .         -eta**(aa(2))*(6698884d0/357075d0*eta-297664d0/14283d0))
     . +2d0*C70bI                                             ! d/dEPSewI
     .   *ALEMMZ*((640d0/1449d0*eta**(14d0/23d0)              ! dEPSewI/dC80BSMI 
     .           -704d0/1725d0*eta**(16d0/23d0)
     .           +32d0/1449d0*eta**(-7d0/23d0)
     .           -32d0/575d0*eta**(-9d0/23d0))/asmb
     .    -dlog(MZ/scb)/pi*8d0/3d0*(ETA**(14d0/23d0)-ETA**(12d0/23d0)))
     . +(asmb/2d0/Pi*(-K78I*C70b-K78*C70bI-ff18I*C10b-ff28I*C20b       ! d/dC80bI
     . +2d0*ff88*C80bI)+(asmb/4d0/Pi)**2*(-2d0*K78I*C71b
     . +2d0*K78*C71bI-ff18I*C11b-ff28I*C21b+2d0*ff88*C81bI)
     . +NP78*C70bI+2d0*NP88*C80bI
     . +asmb/4d0/Pi/(1d0+(50d0/3d0-8d0*Pi**2/3d0)*asmb/4d0/Pi)
     .  *(-2d0*VVtu*0.0136471d0*(C10b+asmb/4d0/Pi*C11b)
     .      +2d0*0.0818824d0*VVtu*(C20b+asmb/4d0/Pi*C21b)
     .      -2d0*0.0144535d0*(C30b+asmb/4d0/Pi*C31b)
     .      -2d0*0.108381d0*(C40b+asmb/4d0/Pi*C41b)
     .      -2d0*0.231256d0*(C50b+asmb/4d0/Pi*C51b)
     .      -2d0*1.06935d0*(C60b+asmb/4d0/Pi*C61b)))
     .                           *eta**(aa(1))                ! dC80bI/dC80BSMI
     . +((asmb/4d0/Pi)**2*(-K78I*C70b+K78*C70bI-ff18I*C10b    ! d/dC81bI
     .         -ff28I*C20b+2d0*ff88*C80bI))
     .                  *eta**(aa(1))*6.7441d0*(eta-1d0)      ! dC81bI/dC80BSMI

*         - Error estimate:
*     First: Variations of BRSG from BSM uncertainties in DC7,8BSM:
      BRSGmax=BRSG+aux*(DABS(KC7BSM)*DC7BSM+DABS(KC8BSM)*DC8BSM
     .                 +DABS(KC7BSMI)*DC7BSMI+DABS(KC8BSMI)*DC8BSMI)
      BRSGmin=BRSG-aux*(DABS(KC7BSM)*DC7BSM+DABS(KC8BSM)*DC8BSM
     .                 +DABS(KC7BSMI)*DC7BSMI+DABS(KC8BSMI)*DC8BSMI)
*     Second: Add the SM + CKM + NP uncertainties +/-0.23d-4 [11] *2sigma
      BRSGmax=BRSGmax+2d0*0.23d-4
      BRSGmin=Max(BRSGmin-2d0*0.23d-4,0d0)
!      print*,'BRBsg',BRSGMIN,BRSG,BRSGMAX,aux*KC7BSM,aux*KC8BSM,
!     . aux*KC7BSMI,aux*KC8BSMI
*     Comparison with experimental data:
      prob(32)=0d0

      IF(BRSGmin.GE.BRSGexpMax)
     .     PROB(32)=BRSGmin/BRSGexpMax-1d0
      IF(BRSGmax.LE.BRSGexpmin)
     .     PROB(32)=BRSGmax/BRSGexpmin-1d0

!      csqb=csqb+4.d0*(BRSG-(BRSGexpmin+BRSGexpMax)/2d0)**2
!     c /((BRSGmin-BRSGMax)**2+(BRSGexpMax-BRSGexpmin)**2)

*	 4) Contributions to B --> Xd gamma

*         - Bremsstrahlung Contribution with delta=0.295
*     Adjusting CKM dependence
      ff27=-0.020604d0+0.017213d0*z+2.9652d0*z**2-17.755d0*z**3
     .    +47.042d0*z**4+delt*(-0.13382d0+2.0837d0*z+10.895*z**2
     .    -172.57d0*z**3+459.87d0*z**4)
      ff27i=0.0020255d0+0.48376d0*z-3.9246d0*z**2+18.816d0*z**3
     .     -55.611d0*z**4+delt*(0.029659d0+4.0743d0*z-59.714d0*z**2
     .     +310.26d0*z**3-673.18d0*z**4)
      aux=ff27
      ff27=(1d0+VVtud)*ff27-VVtudim*ff27i
      ff27i=aux*VVtudim+ff27i*VVtud
      ff17=-1d0/6d0*ff27
      ff17i=-1d0/6d0*ff27i
      ff28=-1d0/3d0*ff27
      ff28i=-1d0/3d0*ff27i
      ff18=-1d0/6d0*ff28
      ff18i=-1d0/6d0*ff28i

*         - Contributions from b -> d uu gamma [26] Eq.(3.1)
      ALEMMB=asmb*23d0/8d0/(460d0*Pi*(1d0-eta)+327d0*asmb*dlog(eta))
     . *(-69d0*Pi+dsqrt(3d0*Pi*(1587d0*Pi+7360d0*Pi*ALEMMZ/asmb
     .   *(1d0-eta)+5232d0*ALEMMZ*dlog(eta))))
      dBPERT=(VVtud**2+VVtudim**2)*(3d0*C20b**2+2d0/3d0*C10b**2)
     . /(1d0-2.41307d0*asmb/Pi-(asmb/Pi)**2*(21.29553d0+4.625050d0
     .  *dlog((scb/mbp)**2))+12d0/23d0*ALEMMB/asmb*(1d0-1d0/eta)
     .  +(lambd1-9d0*lambd2)/2d0/mbp**2)
     . *(0.01454186d0*dlog(29.5d0)-0.02959776d0)    ! quark mass ratio mu/mb=29.5
                                                    ! to recover the central value
                                                    ! of [14] in the SM limit

*         - NLO perturbative coefficients [12,13]
      K17=833d0/729d0-(af(z)+bf(z))/3d0+208d0/243d0*dlog(scb/mb)
     . +2d0*ff17
      K17I=-40d0/243d0*Pi+(afim(z)+bfim(z))/3d0+2d0*ff17i

      K27=-6d0*K17
      K27I=-6d0*K17I

      K37=2392d0/243d0+8d0*Pi/3d0/dsqrt(3d0)+32d0/9d0*(-0.1684d0)
     . -4.0859d0+2d0*0.0316d0+176d0/81d0*dlog(scb/mb)
      K37I=(4d0/9d0-2d0*4d0/81d0)*Pi-56d0/81d0*Pi

      K47=-761d0/729d0-4d0*Pi/9d0/dsqrt(3d0)-16d0/27d0*(-0.1684d0)
     . +4.0859d0/6d0+5d0/3d0*0.0316d0+2d0*bf(z)
     . +152d0/243d0*dlog(scb/mb)+2d0*ff47
      K47I=-(4d0/9d0/6d0+5d0/3d0*4d0/81d0)*Pi+148d0/243d0*Pi
     . -2d0*bfim(z)

      K57=56680d0/243d0+32d0*Pi/3d0/dsqrt(3d0)+128d0/9d0*(-0.1684d0)
     . -16d0*4.0859d0+32d0*0.0316d0+6272d0/81d0*dlog(scb/mb)
      K57I=(16d0*4d0/9d0-32d0*4d0/81d0)*Pi-896d0/81d0*Pi

      K67=5710d0/729d0-16d0*Pi/9d0/dsqrt(3d0)-64d0/27d0*(-0.1684d0)
     . -10d0/3d0*4.0859d0+44d0/3d0*0.0316d0+12d0*af(z)+20d0*bf(z)
     . -4624d0/243d0*dlog(scb/mb)
      K67I=(10d0/3d0*4d0/9d0-44d0/3d0*4d0/81d0)*Pi+2296d0/243d0*Pi
     . -12d0*afim(z)-20d0*bfim(z)

      K77=-91d0/9d0+4d0/9d0*Pi**2-32d0/3d0*dlog(scb/mb)+2d0*ff77

      K78=44d0/9d0-8d0/27d0*Pi**2+32d0/9d0*dlog(scb/mb)+2d0*ff78
      K78I=8d0/9d0*Pi

      BSGPERT=C70b**2+C70bI**2                          ! NNLO expansion [12]
     . +asmb/2d0/Pi*(C70b*C71b+C70bI*C71bI
     .              +ff11*C10b**2+ff12*C10b*C20b+ff22*C20b**2
     .              +K17*C10b*C70b-K17I*C10b*C70bI
     .              +K27*C20b*C70b-K27I*C20b*C70bI
     .              +K37*C30b*C70b-K37I*C30b*C70bI
     .              +K47*C40b*C70b-K47I*C40b*C70bI
     .              +K57*C50b*C70b-K57I*C50b*C70bI
     .              +K67*C60b*C70b-K67I*C60b*C70bI
     .              +K77*(C70b**2+C70bI**2)+K78*(C70b*C80b+C70bI*C80bI)
     .               +K78I*(C80b*C70bI-C70b*C80bI)
     .              +ff18*C10b*C80b-ff18I*C10b*C80bI
     .              +ff28*C20b*C80b-ff18I*C20b*C80bI
     .              +ff48*C40b*C80b
     .              +ff88*(C80b**2+C80bI**2))
     . +(asmb/4d0/Pi)**2*((C71b**2+C71bI**2)
     .   +2d0*(2d0*ff11*C10b*C11b+ff12*(C10b*C21b+C11b*C20b)
     .        +2d0*ff22*C20b*C21b
     .        +K17*(C10b*C71b+C11b*C70b)-K17I*(C10b*C71bI+C11b*C70bI)
     .        +K27*(C20b*C71b+C21b*C70b)-K27I*(C20b*C71bI+C21b*C70bI)
     .        +K37*(C30b*C71b+C31b*C70b)-K37I*(C30b*C71bI+C31b*C70bI)
     .        +K47*(C40b*C71b+C41b*C70b)-K47I*(C40b*C71bI+C41b*C70bI)
     .        +K57*(C50b*C71b+C51b*C70b)-K57I*(C50b*C71bI+C51b*C70bI)
     .        +K67*(C60b*C71b+C61b*C70b)-K67I*(C60b*C71bI+C61b*C70bI)
     .        +2d0*K77*(C70b*C71b+C70bI*C71bI)
     .        +K78*(C70b*C81b+C71b*C80b)
     .        +K78I*(C80b*C71bI+C81b*C70bI-C80bI*C71b-C81bI*C70b)
     .        +ff18*(C10b*C81b+C11b*C80b)-ff18I*(C10b*C81bI+C11b*C80bI)
     .        +ff28*(C20b*C81b+C21b*C80b)-ff28I*(C20b*C81bI+C21b*C80bI)
     .        +ff48*(C40b*C81b+C41b*C80b)
     .        +2d0*ff88*(C80b*C81b+C80bI*C81bI)))
     . -LO4B+dBPERT
     . +2d0*(EPSEW*C70b+EPSEWI*C70bI)+0d-4/aux
     . +6.382d-3                              ! Mimicking the missing SM NNLO
     . -1.266116316d-3*C70BSM+1.419489589d-5*C80BSM         ! Mimicking the missing NNLO[SM] BSM
                                              !  coefficients (from [14])
      Vbdg=(0.2077d0)**2                               ! (V_ts.V_tb/V_cb)^2
      aux=6d0*ALEM0/(pi*CCSL)*Vbdg*BRSL                !prefactor
      BRDG=aux*(BSGPERT+HQET)

*         - New-Physics coefficients (linearized)
      KC7BSM=(2d0*C70b+asmb/2d0/Pi*(C71b+K17*C10b             ! d/dC70b
     . +K27*C20b+K37*C30b+K47*C40b+K57*C50b+K67*C60b+2d0*K77*C70b
     . +K78*C80b-K78I*C80bI)
     . +2d0*(asmb/4d0/Pi)**2*(K17*C11b+K27*C21b+K37*C31b
     . +K47*C41b+K57*C51b+K67*C61b+2d0*K77*C71b+K78*C81b-K78I*C81bI)
     . +2d0*EPSEW
     . +NP17*(C20b-1d0/6d0*C10b)+NP78*C80b+2d0*NP77*C70b)
     .                                     *eta**(aa(2))      ! dC70b/dC70BSM
     . +(asmb/2d0/Pi*C70b+(asmb/4d0/Pi)**2*(2d0*C71b          ! d/dC71b
     . +K17*C10b+K27*C20b+K37*C30b+K47*C40b+K57*C50b+K67*C60b
     . +2d0*K77*C70b+K78*C80b-K78I*C80bI))
     .            *37208d0/4761d0*eta**(aa(2))*(eta-1d0)      ! dC71b/dC70BSM
     . +2d0*C70b                                              ! d/dEPSew
     .         *((88d0/575d0*eta**(16d0/23d0)                 ! dEPSew/dC70BSM 
     .            -40d0/69d0*eta**(-7d0/23d0)
     .            +32d0/75d0*eta**(-9d0/23d0))/asmb
     .            -dlog(MZ/scb)/pi*eta**(aa(2)))*ALEMMZ
     . -1.266116316d-3                                        ! NNLO

      KC7BSMI=(2d0*C70bI+asmb/2d0/Pi*(C71bI-K17I*C10b         ! d/dC70bI
     . -K27I*C20b-K37I*C30b-K47I*C40b-K57I*C50b-K67I*C60b
     . +2d0*K77*C70bI+K78*C80bI+K78I*C80b)
     . +2d0*(asmb/4d0/Pi)**2*(-K17I*C11b-K27I*C21b-K37I*C31b
     . -K47I*C41b-K57I*C51b-K67I*C61b+2d0*K77*C71bI+K78*C81bI
     . +K78I*C81b)+2d0*EPSEWI
     . +NP78*C80bI+2d0*NP77*C70bI)
     .                                     *eta**(aa(2))      ! dC70bI/dC70BSMI
     . +(asmb/2d0/Pi*C70bI+(asmb/4d0/Pi)**2*(2d0*C71bI        ! d/dC71b
     . -K17I*C10b-K27I*C20b-K37I*C30b-K47I*C40b-K57I*C50b-K67I*C60b
     . +2d0*K77*C70bI+K78*C80bI+K78I*C80b))
     .            *37208d0/4761d0*eta**(aa(2))*(eta-1d0)      ! dC71b/dC70BSM
     . +2d0*C70bI                                             ! d/dEPSewI
     .         *((88d0/575d0*eta**(16d0/23d0)                 ! dEPSewI/dC70BSMI 
     .            -40d0/69d0*eta**(-7d0/23d0)
     .            +32d0/75d0*eta**(-9d0/23d0))/asmb
     .            -dlog(MZ/scb)/pi*eta**(aa(2)))*ALEMMZ

      KC8BSM=(2d0*C70b+asmb/2d0/Pi*(C71b+K17*C10b             ! d/dC70b
     . +K27*C20b+K37*C30b+K47*C40b+K57*C50b+K67*C60b+2d0*K77*C70b
     . +K78*C80b)+2d0*(asmb/4d0/Pi)**2*(K17*C11b+K27*C21b+K37*C31b
     . +K47*C41b+K57*C51b+K67*C61b+2d0*K77*C71b+K78*C81b)+2d0*EPSEW
     . +NP17*(C20b-1d0/6d0*C10b)+NP78*C80b+2d0*NP77*C70b)
     .              *8d0/3d0*(eta**(aa(1))-eta**(aa(2)))      ! dC70b/dC80BSM
     . +(asmb/2d0/Pi*C70b+(asmb/4d0/Pi)**2*(2d0*C71b          ! d/dC71b
     . +K17*C10b+K27*C20b+K37*C30b+K47*C40b+K57*C50b+K67*C60b
     . +2d0*K77*C70b+K78*C80b))
     . *(eta**(aa(1))*(256868d0/14283d0*eta-7164416d0/357075d0) ! dC71b/dC80BSM
     .         -eta**(aa(2))*(6698884d0/357075d0*eta-297664d0/14283d0))
     . +2d0*C70b                                              ! d/dEPSew
     .   *ALEMMZ*((640d0/1449d0*eta**(14d0/23d0)              ! dEPSew/dC80BSM 
     .           -704d0/1725d0*eta**(16d0/23d0)
     .           +32d0/1449d0*eta**(-7d0/23d0)
     .           -32d0/575d0*eta**(-9d0/23d0))/asmb
     .    -dlog(MZ/scb)/pi*8d0/3d0*(ETA**(14d0/23d0)-ETA**(12d0/23d0)))
     . +(asmb/2d0/Pi*(K78*C70b+ff18*C10b+ff28*C20b            ! d/dC80b
     . +ff48*C40b+2d0*ff88*C80b)+(asmb/4d0/Pi)**2*(2d0*K78*C71b
     . +ff18*C11b+ff28*C21b+ff48*C41b+2d0*ff88*C81b)
     . +NP78*C70b+2d0*NP88*C80b)
     .                           *eta**(aa(1))                ! dC80b/dC80BSM
     . +((asmb/4d0/Pi)**2*(K78*C70b+ff18*C10b+ff28*C20b       ! d/dC81b
     . +ff48*C40b+2d0*ff88*C80b))
     .                  *eta**(aa(1))*6.7441d0*(eta-1d0)      ! dC81b/dC80BSM
     . +1.419489589d-5                                        ! NNLO

      KC8BSMI=(2d0*C70bI+asmb/2d0/Pi*(C71bI-K17I*C10b         ! d/dC70bI
     . -K27I*C20b-K37I*C30b-K47I*C40b-K57I*C50b-K67I*C60b+2d0*K77*C70bI
     . +K78*C80bI+K78I*C80b)+2d0*(asmb/4d0/Pi)**2*(-K17I*C11b-K27I*C21b
     . -K37I*C31b-K47I*C41b-K57I*C51b-K67I*C61b+2d0*K77*C71bI+K78I*C81b
     . +K78*C81bI)+2d0*EPSEWI+NP78*C80bI+2d0*NP77*C70bI)
     .              *8d0/3d0*(eta**(aa(1))-eta**(aa(2)))      ! dC70bI/dC80BSMI
     . +(asmb/2d0/Pi*C70b+(asmb/4d0/Pi)**2*(2d0*C71bI         ! d/dC71bI
     . -K17I*C10b-K27I*C20b-K37I*C30b-K47I*C40b-K57I*C50b-K67I*C60b
     . +2d0*K77*C70bI+K78I*C80b+K78*C80bI))
     . *(eta**(aa(1))*(256868d0/14283d0*eta-7164416d0/357075d0) ! dC71bI/dC80BSMI
     .         -eta**(aa(2))*(6698884d0/357075d0*eta-297664d0/14283d0))
     . +2d0*C70bI                                             ! d/dEPSewI
     .   *ALEMMZ*((640d0/1449d0*eta**(14d0/23d0)              ! dEPSewI/dC80BSMI 
     .           -704d0/1725d0*eta**(16d0/23d0)
     .           +32d0/1449d0*eta**(-7d0/23d0)
     .           -32d0/575d0*eta**(-9d0/23d0))/asmb
     .    -dlog(MZ/scb)/pi*8d0/3d0*(ETA**(14d0/23d0)-ETA**(12d0/23d0)))
     . +(asmb/2d0/Pi*(-K78I*C70b-K78*C70bI-ff18I*C10b-ff28I*C20b       ! d/dC80bI
     . +2d0*ff88*C80bI)+(asmb/4d0/Pi)**2*(-2d0*K78I*C71b
     . +2d0*K78*C71bI-ff18I*C11b-ff28I*C21b+2d0*ff88*C81bI)
     . +NP78*C70bI+2d0*NP88*C80bI)
     .                           *eta**(aa(1))                ! dC80bI/dC80BSMI
     . +((asmb/4d0/Pi)**2*(-K78I*C70b+K78*C70bI-ff18I*C10b    ! d/dC81bI
     .         -ff28I*C20b+2d0*ff88*C80bI))
     .                  *eta**(aa(1))*6.7441d0*(eta-1d0)      ! dC81bI/dC80BSMI

*         - Error estimate:
*     First: Variations of BRDG from BSM uncertainties in DC7,8BSM:
      BRDGmax=BRDG+aux*(DABS(KC7BSM)*DC7BSM+DABS(KC8BSM)*DC8BSM
     .                 +DABS(KC7BSMI)*DC7BSMI+DABS(KC8BSMI)*DC8BSMI)
      BRDGmin=BRDG-aux*(DABS(KC7BSM)*DC7BSM+DABS(KC8BSM)*DC8BSM
     .                 +DABS(KC7BSMI)*DC7BSMI+DABS(KC8BSMI)*DC8BSMI)
*     Second: Add the SM + CKM + NP uncertainties +0.12/-0.22d-5 [14] *2sigma
      BRDGmax=BRDGmax+2d0*0.12d-5
      BRDGmin=Max(BRDGmin-2d0*0.22d-5,0d0)
!      print*,'BRBdg',BRDGMIN,BRDG,BRDGMAX,aux*KC7BSM,aux*KC8BSM,aux
*     Comparison with experimental data:
      PROB(55)=0d0
      IF(BRDGmin.GE.BRDGexpMax)
     .     PROB(55)=BRDGmin/BRDGexpMax-1d0
      IF(BRDGmax.LE.BRDGexpmin)
     .     PROB(55)=BRDGmax/BRDGexpmin-1d0


*	VI- Delt B = 1 - BR[Bs --> mu+mu-]

*	 1) Wilson coefficients

*         - SM-contribution -> fit of the 3-loop result, Eq.(4) in [18]
      CASM=(0.4802d0*(MT/173.1d0)**(1.52d0)*(alsmz/0.1184d0)**(-0.09)
     .   -0.0112*(MT/173.1d0)**(0.89d0)*(alsmz/0.1184d0)**(-0.09))*2d0

*         - Charged-Higgs [19]
      yt=MT0**2/MHC
      z=MHC/MW**2

      CAH=(MT0/MW/tanb)**2*(fh20(yt)                       ! Z-penguin
     .      +asc0/4d0/Pi*(fh21(yt)+8d0*fh20p(yt)*dlog((sc0/MT0)**2)))
      CSH=mmu*tanb**2/4d0/MW**2*(fh70(xt,z)                ! Box
     .      +asc0/4d0/Pi*(fh71(xt,z)
     .                       +8d0*xt*fh70p(xt,z)*dlog((sc0/MT0)**2)))
      CPH=-CSH
      CSHP=-mmu*tanb**2/4d0/MW**2*(fh30(xt,z)            ! Higgs Penguin
     .      +asc0/4d0/Pi*(fh31(xt,z)
     .                       +8d0*xt*fh30p(xt,z)*dlog((sc0/MT0)**2)))
      CPHP=-CSHP

      aux=0d0
      do i=1,5
      aux=aux+XH(i,4)**2*MH0(i)
     .          *sgn(MH0(i)-M_Bs**2)/
     . dsqrt((MH0(i)-M_Bs**2)**2+MH0(i)*WIDTH(i)**2)
      enddo
      CPHP=CPHP*aux

      aux=0d0
      do i=1,5
      aux=aux+XH(i,2)**2*MH0(i)
     .         *sgn(MH0(i)-M_Bs**2)/
     . dsqrt((MH0(i)-M_Bs**2)**2+MH0(i)*WIDTH(i)**2)
      enddo
      CSHP=CSHP*aux

      CSHPI=0d0
      CPHPI=0d0

*         - Chargino / Squark [19]
      CACHAR=0d0
      CACHARI=0d0

      DO I=1,2
      DO J=1,2
      DO K=1,2                                                      ! Z-penguin
       z=MCH2(J)/MCH2(I)
       yt=MST2(K)/MCH2(I)
       CACHAR=CACHAR+(CCT(J,K,1)*CCT(I,K,1)+CCT(J,K,2)*CCT(I,K,2))  ! 1 stop / 2 charginos
     .     *(2d0*(U(J,1,1)*U(I,1,1)+U(J,1,2)*U(I,1,2))
     .           *dsqrt(MCH2(J)/MCH2(I))*(fc30(z,yt)
     .                        +asf(dsqrt(MST2(K)))/4d0/Pi*fc31(z,yt))
     .      -(V(J,1,1)*V(I,1,1)+V(J,1,2)*V(I,1,2))
     .   *(fc40(z,yt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc41(z,yt)))
     .       +(CCT(J,K,2)*CCT(I,K,1)-CCT(J,K,1)*CCT(I,K,2))
     .     *(2d0*(U(J,1,2)*U(I,1,1)-U(J,1,1)*U(I,1,2))
     .           *dsqrt(MCH2(J)/MCH2(I))*(fc30(z,yt)
     .                        +asf(dsqrt(MST2(K)))/4d0/Pi*fc31(z,yt))
     .      -(V(J,1,2)*V(I,1,1)-V(J,1,1)*V(I,1,2))
     .   *(fc40(z,yt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc41(z,yt)))
       CACHARI=CACHARI+(CCT(J,K,2)*CCT(I,K,1)-CCT(J,K,1)*CCT(I,K,2))
     .     *(2d0*(U(J,1,1)*U(I,1,1)+U(J,1,2)*U(I,1,2))
     .           *dsqrt(MCH2(J)/MCH2(I))*(fc30(z,yt)
     .                        +asf(dsqrt(MST2(K)))/4d0/Pi*fc31(z,yt))
     .      -(V(J,1,1)*V(I,1,1)+V(J,1,2)*V(I,1,2))
     .   *(fc40(z,yt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc41(z,yt)))
     .       -(CCT(J,K,1)*CCT(I,K,1)+CCT(J,K,2)*CCT(I,K,2))
     .     *(2d0*(U(J,1,2)*U(I,1,1)-U(J,1,1)*U(I,1,2))
     .           *dsqrt(MCH2(J)/MCH2(I))*(fc30(z,yt)
     .                        +asf(dsqrt(MST2(K)))/4d0/Pi*fc31(z,yt))
     .      -(V(J,1,2)*V(I,1,1)-V(J,1,1)*V(I,1,2))
     .   *(fc40(z,yt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc41(z,yt)))

       z=MST2(J)/MCH2(I)
       yt=MST2(K)/MCH2(I)                                           ! 2 stops / 1 chargino
       CACHAR=CACHAR+(CCT(I,J,1)*CCT(I,K,1)+CCT(I,J,2)*CCT(I,K,2))
     .      *(UT(K,1,1)*UT(J,1,1)+UT(K,1,2)*UT(J,1,2))
     .              *(fc40(yt,z)+asf(dsqrt(MST2(K)))/4d0/Pi*fc51(yt,z))
     .     +(CCT(I,J,2)*CCT(I,K,1)-CCT(I,J,2)*CCT(I,K,1))
     .      *(UT(K,1,1)*UT(J,1,2)-UT(K,1,2)*UT(J,1,1))
     .              *(fc40(yt,z)+asf(dsqrt(MST2(K)))/4d0/Pi*fc51(yt,z))
       CACHARI=CACHARI+(CCT(I,J,2)*CCT(I,K,1)-CCT(I,J,2)*CCT(I,K,1))
     .      *(UT(K,1,1)*UT(J,1,1)+UT(K,1,2)*UT(J,1,2))
     .              *(fc40(yt,z)+asf(dsqrt(MST2(K)))/4d0/Pi*fc51(yt,z))
     .     -(CCT(I,J,1)*CCT(I,K,1)+CCT(I,J,2)*CCT(I,K,2))
     .      *(UT(K,1,1)*UT(J,1,2)-UT(K,1,2)*UT(J,1,1))
     .              *(fc40(yt,z)+asf(dsqrt(MST2(K)))/4d0/Pi*fc51(yt,z))
      ENDDO
       z=MCH2(J)/MCH2(I)
       yt=MSU2(1)/MCH2(I)
       CACHAR=CACHAR                                                ! 1 scharm / 2charginos
     .  +VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .      *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .      *(2d0*(U(J,1,1)*U(I,1,1)+U(J,1,2)*U(I,1,2))
     .              *dsqrt(MCH2(J)/MCH2(I))*(fc30(z,yt)
     .                          +asf(dsqrt(MSU2(1)))/4d0/Pi*fc31(z,yt))
     .      -(V(J,1,1)*V(I,1,1)+V(J,1,2)*V(I,1,2))
     .     *(fc40(z,yt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc41(z,yt)))
     .  +VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .      *(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))
     .      *(2d0*(U(J,1,2)*U(I,1,1)-U(J,1,1)*U(I,1,2))
     .              *dsqrt(MCH2(J)/MCH2(I))*(fc30(z,yt)
     .                          +asf(dsqrt(MSU2(1)))/4d0/Pi*fc31(z,yt))
     .      -(V(J,1,2)*V(I,1,1)-V(J,1,1)*V(I,1,2))
     .     *(fc40(z,yt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc41(z,yt)))
       CACHARI=CACHARI
     .  +VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .      *(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))
     .      *(2d0*(U(J,1,1)*U(I,1,1)+U(J,1,2)*U(I,1,2))
     .              *dsqrt(MCH2(J)/MCH2(I))*(fc30(z,yt)
     .                          +asf(dsqrt(MSU2(1)))/4d0/Pi*fc31(z,yt))
     .      -(V(J,1,1)*V(I,1,1)+V(J,1,2)*V(I,1,2))
     .     *(fc40(z,yt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc41(z,yt)))
     .  -VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .      *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .      *(2d0*(U(J,1,2)*U(I,1,1)-U(J,1,1)*U(I,1,2))
     .              *dsqrt(MCH2(J)/MCH2(I))*(fc30(z,yt)
     .                          +asf(dsqrt(MSU2(1)))/4d0/Pi*fc31(z,yt))
     .      -(V(J,1,2)*V(I,1,1)-V(J,1,1)*V(I,1,2))
     .     *(fc40(z,yt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc41(z,yt)))
      ENDDO

       z=MSU2(1)/MCH2(I)                                            ! 2 scharms / 1 chargino
       CACHAR=CACHAR
     .  +VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .            *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .    *(V(I,1,1)**2+V(I,1,2)**2)
     .    *(fc40(z,z)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc51(z,z))
      ENDDO
      CACHAR=CACHAR/8d0

      aux=0d0
      auxe=0d0
      aux1=0d0
      aux2=0d0
      CSCHAR=0d0
      CPCHAR=0d0
      CSCHARI=0d0
      CPCHARI=0d0
      DO I=1,2
      DO J=1,2
      DO K=1,2                                            ! Box
       z=MCH2(J)/MCH2(I)
       yt=MST2(K)/MCH2(I)
       xt=MSNE2/MCH2(I)                                      ! 1 stop / 1 sneutrino / 2 charginos
       aux=aux+(CCT(J,K,1)*CCT(I,K,1)+CCT(J,K,2)*CCT(I,K,2))/MCH2(I)*(
     .          (V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .        *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc81(z,yt,xt))
     . +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .    *(U(I,2,1)*U(J,2,1)+U(I,2,2)*U(J,2,2))*fc60(z,yt,xt))
!     .                          +asf(ST(K))/4d0/Pi*fc91(z,yt,xt))
     .    -(CCT(J,K,2)*CCT(I,K,1)-CCT(J,K,1)*CCT(I,K,2))/MCH2(I)*(
     .          (V(I,1,1)*V(J,1,2)-V(I,1,2)*V(J,1,1))
     .        *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc81(z,yt,xt))
     . +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .    *(U(I,2,2)*U(J,2,1)-U(I,2,1)*U(J,2,2))*fc60(z,yt,xt))
!     .                          +asf(ST(K))/4d0/Pi*fc91(z,yt,xt))
       auxe=auxe+(CCT(J,K,1)*CCT(I,K,1)+CCT(J,K,2)*CCT(I,K,2))/MCH2(I)
     .        *((V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .      *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc81(z,yt,xt)))
     .          -(CCT(J,K,2)*CCT(I,K,1)-CCT(J,K,1)*CCT(I,K,2))/MCH2(I)
     .        *((V(I,1,1)*V(J,1,2)-V(I,1,2)*V(J,1,1))
     .      *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc81(z,yt,xt)))
       aux1=aux1+(CCT(J,K,2)*CCT(I,K,1)-CCT(J,K,1)*CCT(I,K,2))/MCH2(I)
     .          *((V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .        *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc81(z,yt,xt))
     . +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .    *(U(I,2,1)*U(J,2,1)+U(I,2,2)*U(J,2,2))*fc60(z,yt,xt))
!     .                          +asf(ST(K))/4d0/Pi*fc91(z,yt,xt))
     .    +(CCT(J,K,1)*CCT(I,K,1)+CCT(J,K,1)*CCT(I,K,2))/MCH2(I)*(
     .          (V(I,1,1)*V(J,1,2)-V(I,1,2)*V(J,1,1))
     .        *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc81(z,yt,xt))
     . +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .    *(U(I,2,2)*U(J,2,1)-U(I,2,1)*U(J,2,2))*fc60(z,yt,xt))
!     .                          +asf(ST(K))/4d0/Pi*fc91(z,yt,xt))
       aux2=aux2+(CCT(J,K,2)*CCT(I,K,1)-CCT(J,K,1)*CCT(I,K,2))/MCH2(I)
     .        *((V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .      *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc81(z,yt,xt)))
     .          +(CCT(J,K,1)*CCT(I,K,1)+CCT(J,K,2)*CCT(I,K,2))/MCH2(I)
     .        *((V(I,1,1)*V(J,1,2)-V(I,1,2)*V(J,1,1))
     .      *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc81(z,yt,xt)))
       CSCHAR=CSCHAR+(CCT(J,K,1)*UT(K,1,1)*U(I,2,1)
     .               -CCT(J,K,1)*UT(K,1,2)*U(I,2,2)
     .               +CCT(J,K,2)*UT(K,1,2)*U(I,2,1)
     .               +CCT(J,K,2)*UT(K,1,1)*U(I,2,2))/MCH2(I)
     . *(1d0+asf(dsqrt(MST2(K)))/4d0/Pi*(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .     *mmu/vdq**2/(1d0+epst3)
     . *((V(J,1,1)*U(I,2,1)-V(J,1,2)*U(I,2,2))
     .   *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc121(z,yt,xt))
     .  +dsqrt(MCH2(J)/MCH2(I))*(U(J,2,1)*V(I,1,1)-U(J,2,2)*V(I,1,2))
     .   *(fc60(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc131(z,yt,xt)))
     .  +(CCT(J,K,2)*UT(K,1,1)*U(I,2,1)
     .               -CCT(J,K,2)*UT(K,1,2)*U(I,2,2)
     .               -CCT(J,K,1)*UT(K,1,2)*U(I,2,1)
     .               -CCT(J,K,1)*UT(K,1,1)*U(I,2,2))/MCH2(I)
     . *(1d0+asf(dsqrt(MST2(K)))/4d0/Pi*(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .     *mmu/vdq**2/(1d0+epst3)
     . *(-(V(J,1,2)*U(I,2,1)+V(J,1,1)*U(I,2,2))
     .   *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc121(z,yt,xt))
     .  +dsqrt(MCH2(J)/MCH2(I))*(U(J,2,2)*V(I,1,1)+U(J,2,1)*V(I,1,2))
     .   *(fc60(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc131(z,yt,xt)))
       CSCHARI=CSCHARI+(CCT(J,K,2)*UT(K,1,1)*U(I,2,1)
     .               -CCT(J,K,2)*UT(K,1,2)*U(I,2,2)
     .               -CCT(J,K,1)*UT(K,1,2)*U(I,2,1)
     .               -CCT(J,K,1)*UT(K,1,1)*U(I,2,2))/MCH2(I)
     . *(1d0+asf(dsqrt(MST2(K)))/4d0/Pi*(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .     *mmu/vdq**2/(1d0+epst3)
     . *((V(J,1,1)*U(I,2,1)-V(J,1,2)*U(I,2,2))
     .   *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc121(z,yt,xt))
     .  +dsqrt(MCH2(J)/MCH2(I))*(U(J,2,1)*V(I,1,1)-U(J,2,2)*V(I,1,2))
     .   *(fc60(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc131(z,yt,xt)))
     .  +(CCT(J,K,1)*UT(K,1,1)*U(I,2,1)
     .               -CCT(J,K,1)*UT(K,1,2)*U(I,2,2)
     .               +CCT(J,K,2)*UT(K,1,2)*U(I,2,1)
     .               +CCT(J,K,2)*UT(K,1,1)*U(I,2,2))/MCH2(I)
     . *(1d0+asf(dsqrt(MST2(K)))/4d0/Pi*(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .     *mmu/vdq**2/(1d0+epst3)
     . *((V(J,1,2)*U(I,2,1)+V(J,1,1)*U(I,2,2))
     .   *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc121(z,yt,xt))
     .  -dsqrt(MCH2(J)/MCH2(I))*(U(J,2,2)*V(I,1,1)+U(J,2,1)*V(I,1,2))
     .   *(fc60(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc131(z,yt,xt)))
       CPCHAR=CPCHAR-(CCT(J,K,1)*UT(K,1,1)*U(I,2,1)
     .               -CCT(J,K,1)*UT(K,1,2)*U(I,2,2)
     .               +CCT(J,K,2)*UT(K,1,2)*U(I,2,1)
     .               +CCT(J,K,2)*UT(K,1,1)*U(I,2,2))/MCH2(I)
     . *(1d0+asf(dsqrt(MST2(K)))/4d0/Pi*(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .     *mmu/vdq**2/(1d0+epst3)
     . *((V(J,1,1)*U(I,2,1)-V(J,1,2)*U(I,2,2))
     .    *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc121(z,yt,xt))
     .  -dsqrt(MCH2(J)/MCH2(I))*(U(J,2,1)*V(I,1,1)-U(J,2,2)*V(I,1,2))
     .      *(fc60(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc131(z,yt,xt)))
     .   -(CCT(J,K,2)*UT(K,1,1)*U(I,2,1)
     .               -CCT(J,K,2)*UT(K,1,2)*U(I,2,2)
     .               -CCT(J,K,1)*UT(K,1,2)*U(I,2,1)
     .               -CCT(J,K,1)*UT(K,1,1)*U(I,2,2))/MCH2(I)
     . *(1d0+asf(dsqrt(MST2(K)))/4d0/Pi*(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .     *mmu/vdq**2/(1d0+epst3)
     . *(-(V(J,1,2)*U(I,2,1)+V(J,1,1)*U(I,2,2))
     .    *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc121(z,yt,xt))
     .  -dsqrt(MCH2(J)/MCH2(I))*(U(J,2,2)*V(I,1,1)+U(J,2,1)*V(I,1,2))
     .      *(fc60(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc131(z,yt,xt)))
       CPCHARI=CPCHARI-(CCT(J,K,2)*UT(K,1,1)*U(I,2,1)
     .               -CCT(J,K,2)*UT(K,1,2)*U(I,2,2)
     .               -CCT(J,K,1)*UT(K,1,2)*U(I,2,1)
     .               -CCT(J,K,1)*UT(K,1,1)*U(I,2,2))/MCH2(I)
     . *(1d0+asf(dsqrt(MST2(K)))/4d0/Pi*(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .     *mmu/vdq**2/(1d0+epst3)
     . *((V(J,1,1)*U(I,2,1)-V(J,1,2)*U(I,2,2))
     .    *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc121(z,yt,xt))
     .  -dsqrt(MCH2(J)/MCH2(I))*(U(J,2,1)*V(I,1,1)-U(J,2,2)*V(I,1,2))
     .      *(fc60(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc131(z,yt,xt)))
     .   -(CCT(J,K,1)*UT(K,1,1)*U(I,2,1)
     .               -CCT(J,K,1)*UT(K,1,2)*U(I,2,2)
     .               +CCT(J,K,2)*UT(K,1,2)*U(I,2,1)
     .               +CCT(J,K,2)*UT(K,1,1)*U(I,2,2))/MCH2(I)
     . *(1d0+asf(dsqrt(MST2(K)))/4d0/Pi*(1d0+2d0*dlog(MST2(K)/MGL**2)))
     .     *mmu/vdq**2/(1d0+epst3)
     . *((V(J,1,2)*U(I,2,1)+V(J,1,1)*U(I,2,2))
     .    *(fc50(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc121(z,yt,xt))
     .  +dsqrt(MCH2(J)/MCH2(I))*(U(J,2,2)*V(I,1,1)+U(J,2,1)*V(I,1,2))
     .      *(fc60(z,yt,xt)+asf(dsqrt(MST2(K)))/4d0/Pi*fc131(z,yt,xt)))
      ENDDO
       z=MCH2(J)/MCH2(I)
       yt=MSU2(1)/MCH2(I)
       xt=MSNE2/MCH2(I)
       aux=aux                                    ! 1 scharm / 1 sneutrino / 2 charginos
     .  +VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .    *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))/MCH2(I)
     .   *((V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .       *(fc50(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc81(z,yt,xt))
     . +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .   *(U(I,2,1)*U(J,2,1)+U(I,2,2)*U(J,2,2))*fc60(z,yt,xt))
!     .                          +asf(MUL)/4d0/Pi*fc91(z,yt,xt)))
     .  +VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .    *(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))/MCH2(I)
     .   *((V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))
     .       *(fc50(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc81(z,yt,xt))
     . +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .   *(U(I,2,1)*U(J,2,2)-U(I,2,2)*U(J,2,1))*fc60(z,yt,xt))
!     .                          +asf(MUL)/4d0/Pi*fc91(z,yt,xt)))
       auxe=auxe
     .  +VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2/MCH2(I)
     .               *((V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))**2
     .                +(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))**2)
     .   *(fc50(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc81(z,yt,xt))
       aux1=aux1
     .  +VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .    *(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))/MCH2(I)
     .   *((V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .       *(fc50(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc81(z,yt,xt))
     . +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .   *(U(I,2,1)*U(J,2,1)+U(I,2,2)*U(J,2,2))*fc60(z,yt,xt))
!     .                          +asf(MUL)/4d0/Pi*fc91(z,yt,xt)))
     .  -VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .    *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))/MCH2(I)
     .   *((V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))
     .       *(fc50(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc81(z,yt,xt))
     . +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .   *(U(I,2,1)*U(J,2,2)-U(I,2,2)*U(J,2,1))*fc60(z,yt,xt))
!     .                          +asf(MUL)/4d0/Pi*fc91(z,yt,xt)))
       aux2=aux2
     .  +VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2/MCH2(I)
     .               *((V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))
     .                *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .                -(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .                *(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2)))
     .   *(fc50(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc81(z,yt,xt))
       CSCHAR=CSCHAR
     . +VVc*(V(J,1,1)*U(I,2,1)-V(J,1,2)*U(I,2,2))/MCH2(I)
     .   *mmu/vdq**2/(1d0+epst3)
     .      *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                   *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi
     .                   *(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *((V(J,1,1)*U(I,2,1)-V(J,1,2)*U(I,2,2))
     .    *(fc50(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc121(z,yt,xt))
     .  +dsqrt(MCH2(J)/MCH2(I))*(U(J,2,1)*V(I,1,1)-U(J,2,2)*V(I,1,2))
     .    *(fc60(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc131(z,yt,xt)))
     . +VVc*(V(J,1,2)*U(I,2,1)+V(J,1,1)*U(I,2,2))/MCH2(I)
     .   *mmu/vdq**2/(1d0+epst3)
     .      *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                   *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi
     .                   *(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *((V(J,1,2)*U(I,2,1)+V(J,1,1)*U(I,2,2))
     .    *(fc50(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc121(z,yt,xt))
     .  +dsqrt(MCH2(J)/MCH2(I))*(-U(J,2,2)*V(I,1,1)-U(J,2,1)*V(I,1,2))
     .    *(fc60(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc131(z,yt,xt)))
       CSCHARI=CSCHARI
     . -VVc*(V(J,1,2)*U(I,2,1)+V(J,1,1)*U(I,2,2))/MCH2(I)
     .   *mmu/vdq**2/(1d0+epst3)
     .      *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                   *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi
     .                   *(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *((V(J,1,1)*U(I,2,1)-V(J,1,2)*U(I,2,2))
     .    *(fc50(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc121(z,yt,xt))
     .  +dsqrt(MCH2(J)/MCH2(I))*(U(J,2,1)*V(I,1,1)-U(J,2,2)*V(I,1,2))
     .    *(fc60(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc131(z,yt,xt)))
     . +VVc*(V(J,1,1)*U(I,2,1)-V(J,1,2)*U(I,2,2))/MCH2(I)
     .   *mmu/vdq**2/(1d0+epst3)
     .      *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                   *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi
     .                   *(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *((V(J,1,2)*U(I,2,1)+V(J,1,1)*U(I,2,2))
     .    *(fc50(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc121(z,yt,xt))
     .  +dsqrt(MCH2(J)/MCH2(I))*(-U(J,2,2)*V(I,1,1)-U(J,2,1)*V(I,1,2))
     .    *(fc60(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc131(z,yt,xt)))
       CPCHAR=CPCHAR
     . -VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                    *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     .         /MCH2(I)*mmu/vdq**2/(1d0+epst3)
     .     *(V(J,1,1)*U(I,2,1)-V(J,1,2)*U(I,2,2))
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi*(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *((V(J,1,1)*U(I,2,1)-V(J,1,2)*U(I,2,2))
     .     *(fc50(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc121(z,yt,xt))
     .  -dsqrt(MCH2(J)/MCH2(I))*(U(J,2,1)*V(I,1,1)-U(J,2,2)*V(I,1,2))
     . *(fc60(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc131(z,yt,xt)))
     . -VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                    *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     .         /MCH2(I)*mmu/vdq**2/(1d0+epst3)
     .     *(V(J,1,2)*U(I,2,1)+V(J,1,1)*U(I,2,2))
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi*(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *((V(J,1,2)*U(I,2,1)+V(J,1,1)*U(I,2,2))
     .     *(fc50(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc121(z,yt,xt))
     .  -dsqrt(MCH2(J)/MCH2(I))*(-U(J,2,2)*V(I,1,1)-U(J,2,1)*V(I,1,2))
     . *(fc60(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc131(z,yt,xt)))
       CPCHARI=CPCHARI
     . +VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                    *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     .         /MCH2(I)*mmu/vdq**2/(1d0+epst3)
     .     *(V(J,1,2)*U(I,2,1)+V(J,1,1)*U(I,2,2))
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi*(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *((V(J,1,1)*U(I,2,1)-V(J,1,2)*U(I,2,2))
     .     *(fc50(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc121(z,yt,xt))
     .  -dsqrt(MCH2(J)/MCH2(I))*(U(J,2,1)*V(I,1,1)-U(J,2,2)*V(I,1,2))
     . *(fc60(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc131(z,yt,xt)))
     . -VVc*(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                    *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     .         /MCH2(I)*mmu/vdq**2/(1d0+epst3)
     .     *(V(J,1,1)*U(I,2,1)-V(J,1,2)*U(I,2,2))
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi*(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *((V(J,1,2)*U(I,2,1)+V(J,1,1)*U(I,2,2))
     .     *(fc50(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc121(z,yt,xt))
     .  -dsqrt(MCH2(J)/MCH2(I))*(-U(J,2,2)*V(I,1,1)-U(J,2,1)*V(I,1,2))
     . *(fc60(z,yt,xt)+asf(dsqrt(MSU2(1)))/4d0/Pi*fc131(z,yt,xt)))
      ENDDO
      ENDDO

      DO I=1,2
      DO J=1,2
      DO K=1,2
      DO L=1,2                                            ! Quartic squark coupling
      DO M=1,2
       z=MCH2(J)/MCH2(I)                       ! stops
       yt=MST2(K)/MCH2(I)
       zt=MST2(L)/MCH2(I)
       xt=MSNE2/MCH2(I)
       aux=aux-asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     .           *((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *((V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .     +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .       *(U(I,2,1)*U(J,2,1)+U(I,2,2)*U(J,2,2))*fc100(z,yt,zt,xt))
     .   -(CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *((V(I,1,1)*V(J,1,2)-V(I,1,2)*V(J,1,1))*fc90(z,yt,zt,xt)
     .     +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .       *(U(I,2,2)*U(J,2,1)-U(I,2,2)*U(J,2,1))*fc100(z,yt,zt,xt)))
     .        -asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     .           *((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *((V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .     +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .       *(U(I,2,1)*U(J,2,1)+U(I,2,2)*U(J,2,2))*fc100(z,yt,zt,xt))
     .   +(CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *((V(I,1,1)*V(J,1,2)-V(I,1,2)*V(J,1,1))*fc90(z,yt,zt,xt)
     .     +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .       *(U(I,2,2)*U(J,2,1)-U(I,2,2)*U(J,2,1))*fc100(z,yt,zt,xt)))
     . -asf(dsqrt(MST2(M)))/6d0/Pi/MW**2*MST2(M)/MCH2(I)
     .           *((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,1)+U(I,1,2)*U(J,1,2))
     .     -(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc50(z,yt,zt))
     .  +(CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,2)-U(I,1,2)*U(J,1,1))
     .     -(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))*fc50(z,yt,zt)))
     . -asf(dsqrt(MST2(M)))/6d0/Pi/MW**2*MST2(M)/MCH2(I)
     .           *((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,1)+U(I,1,2)*U(J,1,2))
     .     -(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc50(z,yt,zt))
     .  -(CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,2)-U(I,1,2)*U(J,1,1))
     .     -(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))*fc50(z,yt,zt)))
       auxe=auxe-asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     .           *((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .   -(CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *(V(I,1,1)*V(J,1,2)-V(I,1,2)*V(J,1,1))*fc90(z,yt,zt,xt))
     .        -asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     .           *((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .   +(CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *(V(I,1,1)*V(J,1,2)-V(I,1,2)*V(J,1,1))*fc90(z,yt,zt,xt))
     . -asf(dsqrt(MST2(M)))/6d0/Pi/MW**2*MST2(M)/MCH2(I)
     .           *((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,1)+U(I,1,2)*U(J,1,2))
     .     -(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc50(z,yt,zt))
     .  +(CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,2)-U(I,1,2)*U(J,1,1))
     .     -(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))*fc50(z,yt,zt)))
     . -asf(dsqrt(MST2(M)))/6d0/Pi/MW**2*MST2(M)/MCH2(I)
     .           *((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,1)+U(I,1,2)*U(J,1,2))
     .     -(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc50(z,yt,zt))
     .  -(CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,2)-U(I,1,2)*U(J,1,1))
     .     -(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))*fc50(z,yt,zt)))
       aux1=aux1+asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     .           *((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *((V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .     +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .       *(U(I,2,1)*U(J,2,1)+U(I,2,2)*U(J,2,2))*fc100(z,yt,zt,xt))
     .   -(CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *((V(I,1,1)*V(J,1,2)-V(I,1,2)*V(J,1,1))*fc90(z,yt,zt,xt)
     .     +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .       *(U(I,2,2)*U(J,2,1)-U(I,2,2)*U(J,2,1))*fc100(z,yt,zt,xt)))
     .        -asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     .           *((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *((V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .     +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .       *(U(I,2,1)*U(J,2,1)+U(I,2,2)*U(J,2,2))*fc100(z,yt,zt,xt))
     .   +(CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *((V(I,1,1)*V(J,1,2)-V(I,1,2)*V(J,1,1))*fc90(z,yt,zt,xt)
     .     +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .       *(U(I,2,2)*U(J,2,1)-U(I,2,2)*U(J,2,1))*fc100(z,yt,zt,xt)))
     . +asf(dsqrt(MST2(M)))/6d0/Pi/MW**2*MST2(M)/MCH2(I)
     .           *((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,1)+U(I,1,2)*U(J,1,2))
     .     -(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc50(z,yt,zt))
     .  +(CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,2)-U(I,1,2)*U(J,1,1))
     .     -(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))*fc50(z,yt,zt)))
     . -asf(dsqrt(MST2(M)))/6d0/Pi/MW**2*MST2(M)/MCH2(I)
     .           *((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,1)+U(I,1,2)*U(J,1,2))
     .     -(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc50(z,yt,zt))
     .  -(CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,2)-U(I,1,2)*U(J,1,1))
     .     -(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))*fc50(z,yt,zt)))
       aux2=aux2+asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     .           *((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .   -(CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *(V(I,1,1)*V(J,1,2)-V(I,1,2)*V(J,1,1))*fc90(z,yt,zt,xt))
     .        -asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     .           *((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .   +(CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *(V(I,1,1)*V(J,1,2)-V(I,1,2)*V(J,1,1))*fc90(z,yt,zt,xt))
     . +asf(dsqrt(MST2(M)))/6d0/Pi/MW**2*MST2(M)/MCH2(I)
     .           *((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,1)+U(I,1,2)*U(J,1,2))
     .     -(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc50(z,yt,zt))
     .  +(CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,2)-U(I,1,2)*U(J,1,1))
     .     -(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))*fc50(z,yt,zt)))
     . -asf(dsqrt(MST2(M)))/6d0/Pi/MW**2*MST2(M)/MCH2(I)
     .           *((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .             -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .             -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .            +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .             -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .            *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .             -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,2)*CCT(I,L,1)-CCT(J,K,1)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,1)+U(I,1,2)*U(J,1,2))
     .     -(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc50(z,yt,zt))
     .  -(CCT(J,K,1)*CCT(I,L,1)+CCT(J,K,2)*CCT(I,L,2))
     .  *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,zt)
     .     *(U(I,1,1)*U(J,1,2)-U(I,1,2)*U(J,1,1))
     .     -(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))*fc50(z,yt,zt)))
       CSCHAR=CSCHAR-asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     . *(1d0+asf(dsqrt(MST2(M)))/4d0/Pi*(1d0+2d0*dlog(MST2(M)/MGL**2)))
     .       *mmu/vdq**2/(1d0+epst3)
     .          *((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .            -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .            -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .           +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .            -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .            -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,1)*UT(L,1,1)*U(I,2,1)
     .            -CCT(J,K,1)*UT(L,1,2)*U(I,2,2)
     .            +CCT(J,K,2)*UT(L,1,2)*U(I,2,1)
     .            +CCT(J,K,2)*UT(L,1,1)*U(I,2,2))
     .   *((U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .       +(V(I,1,1)*U(J,2,1)-V(I,1,2)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt))
     .   -(CCT(J,K,2)*UT(L,1,1)*U(I,2,1)
     .            -CCT(J,K,2)*UT(L,1,2)*U(I,2,2)
     .            -CCT(J,K,1)*UT(L,1,2)*U(I,2,1)
     .            -CCT(J,K,1)*UT(L,1,1)*U(I,2,2))
     .   *((U(I,2,1)*V(J,1,2)+U(I,2,2)*V(J,1,1))*fc90(z,yt,zt,xt)
     .       +(V(I,1,2)*U(J,2,1)+V(I,1,1)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt)))
     .   -asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     . *(1d0+asf(dsqrt(MST2(M)))/4d0/Pi*(1d0+2d0*dlog(MST2(M)/MGL**2)))
     .       *mmu/vdq**2/(1d0+epst3)
     .          *((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .            -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .            -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .           -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .            -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .            -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,2)*UT(L,1,1)*U(I,2,1)
     .            -CCT(J,K,2)*UT(L,1,2)*U(I,2,2)
     .            -CCT(J,K,1)*UT(L,1,2)*U(I,2,1)
     .            -CCT(J,K,1)*UT(L,1,1)*U(I,2,2))
     .   *((U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .       +(V(I,1,1)*U(J,2,1)-V(I,1,2)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt))
     .   +(CCT(J,K,1)*UT(L,1,1)*U(I,2,1)
     .            -CCT(J,K,1)*UT(L,1,2)*U(I,2,2)
     .            +CCT(J,K,2)*UT(L,1,2)*U(I,2,1)
     .            +CCT(J,K,2)*UT(L,1,1)*U(I,2,2))
     .   *((U(I,2,1)*V(J,1,2)+U(I,2,2)*V(J,1,1))*fc90(z,yt,zt,xt)
     .       +(V(I,1,2)*U(J,2,1)+V(I,1,1)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt)))
       CSCHARI=CSCHARI+asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     . *(1d0+asf(dsqrt(MST2(M)))/4d0/Pi*(1d0+2d0*dlog(MST2(M)/MGL**2)))
     .       *mmu/vdq**2/(1d0+epst3)
     .          *((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .            -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .            -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .           -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .            -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .            -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,1)*UT(L,1,1)*U(I,2,1)
     .            -CCT(J,K,1)*UT(L,1,2)*U(I,2,2)
     .            +CCT(J,K,2)*UT(L,1,2)*U(I,2,1)
     .            +CCT(J,K,2)*UT(L,1,1)*U(I,2,2))
     .   *((U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .       +(V(I,1,1)*U(J,2,1)-V(I,1,2)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt))
     .   -(CCT(J,K,2)*UT(L,1,1)*U(I,2,1)
     .            -CCT(J,K,2)*UT(L,1,2)*U(I,2,2)
     .            -CCT(J,K,1)*UT(L,1,2)*U(I,2,1)
     .            -CCT(J,K,1)*UT(L,1,1)*U(I,2,2))
     .   *((U(I,2,1)*V(J,1,2)+U(I,2,2)*V(J,1,1))*fc90(z,yt,zt,xt)
     .       +(V(I,1,2)*U(J,2,1)+V(I,1,1)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt)))
     .   -asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     . *(1d0+asf(dsqrt(MST2(M)))/4d0/Pi*(1d0+2d0*dlog(MST2(M)/MGL**2)))
     .       *mmu/vdq**2/(1d0+epst3)
     .          *((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .            -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .            -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .           +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .            -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .            -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,2)*UT(L,1,1)*U(I,2,1)
     .            -CCT(J,K,2)*UT(L,1,2)*U(I,2,2)
     .            -CCT(J,K,1)*UT(L,1,2)*U(I,2,1)
     .            -CCT(J,K,1)*UT(L,1,1)*U(I,2,2))
     .   *((U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .       +(V(I,1,1)*U(J,2,1)-V(I,1,2)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt))
     .   +(CCT(J,K,1)*UT(L,1,1)*U(I,2,1)
     .            -CCT(J,K,1)*UT(L,1,2)*U(I,2,2)
     .            +CCT(J,K,2)*UT(L,1,2)*U(I,2,1)
     .            +CCT(J,K,2)*UT(L,1,1)*U(I,2,2))
     .   *((U(I,2,1)*V(J,1,2)+U(I,2,2)*V(J,1,1))*fc90(z,yt,zt,xt)
     .       +(V(I,1,2)*U(J,2,1)+V(I,1,1)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt)))
       CPCHAR=CPCHAR+asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     . *(1d0+asf(dsqrt(MST2(M)))/4d0/Pi*(1d0+2d0*dlog(MST2(M)/MGL**2)))
     .         *mmu/vdq**2/(1d0+epst3)
     .          *((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .            -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .            -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .           +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .            -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .            -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,1)*UT(L,1,1)*U(I,2,1)
     .         -CCT(J,K,1)*UT(L,1,2)*U(I,2,2)
     .         +CCT(J,K,2)*UT(L,1,2)*U(I,2,1)
     .         +CCT(J,K,2)*UT(L,1,1)*U(I,2,2))
     .    *((U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .       -(V(I,1,1)*U(J,2,1)-V(I,1,2)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt))
     .   -(CCT(J,K,2)*UT(L,1,1)*U(I,2,1)
     .         -CCT(J,K,2)*UT(L,1,2)*U(I,2,2)
     .         -CCT(J,K,1)*UT(L,1,2)*U(I,2,1)
     .         -CCT(J,K,1)*UT(L,1,1)*U(I,2,2))
     .    *((U(I,2,2)*V(J,1,1)+U(I,2,1)*V(J,1,2))*fc90(z,yt,zt,xt)
     .       -(V(I,1,2)*U(J,2,1)+V(I,1,1)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt)))
     . +asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     . *(1d0+asf(dsqrt(MST2(M)))/4d0/Pi*(1d0+2d0*dlog(MST2(M)/MGL**2)))
     .         *mmu/vdq**2/(1d0+epst3)
     .          *((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .            -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .            -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .           -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .            -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .            -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,2)*UT(L,1,1)*U(I,2,1)
     .         -CCT(J,K,2)*UT(L,1,2)*U(I,2,2)
     .         -CCT(J,K,1)*UT(L,1,2)*U(I,2,1)
     .         -CCT(J,K,1)*UT(L,1,1)*U(I,2,2))
     .    *((U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .       -(V(I,1,1)*U(J,2,1)-V(I,1,2)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt))
     .   +(CCT(J,K,1)*UT(L,1,1)*U(I,2,1)
     .         -CCT(J,K,1)*UT(L,1,2)*U(I,2,2)
     .         +CCT(J,K,2)*UT(L,1,2)*U(I,2,1)
     .         +CCT(J,K,2)*UT(L,1,1)*U(I,2,2))
     .    *((U(I,2,2)*V(J,1,1)+U(I,2,1)*V(J,1,2))*fc90(z,yt,zt,xt)
     .       -(V(I,1,2)*U(J,2,1)+V(I,1,1)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt)))
       CPCHARI=CPCHARI-asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     . *(1d0+asf(dsqrt(MST2(M)))/4d0/Pi*(1d0+2d0*dlog(MST2(M)/MGL**2)))
     .         *mmu/vdq**2/(1d0+epst3)
     .          *((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .            -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .            -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .           -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .            -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .            -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,1)*UT(L,1,1)*U(I,2,1)
     .         -CCT(J,K,1)*UT(L,1,2)*U(I,2,2)
     .         +CCT(J,K,2)*UT(L,1,2)*U(I,2,1)
     .         +CCT(J,K,2)*UT(L,1,1)*U(I,2,2))
     .    *((U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .       -(V(I,1,1)*U(J,2,1)-V(I,1,2)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt))
     .   -(CCT(J,K,2)*UT(L,1,1)*U(I,2,1)
     .         -CCT(J,K,2)*UT(L,1,2)*U(I,2,2)
     .         -CCT(J,K,1)*UT(L,1,2)*U(I,2,1)
     .         -CCT(J,K,1)*UT(L,1,1)*U(I,2,2))
     .    *((U(I,2,2)*V(J,1,1)+U(I,2,1)*V(J,1,2))*fc90(z,yt,zt,xt)
     .       -(V(I,1,2)*U(J,2,1)+V(I,1,1)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt)))
     . +asf(dsqrt(MST2(M)))/3d0/Pi*MST2(M)/MCH2(I)**2
     . *(1d0+asf(dsqrt(MST2(M)))/4d0/Pi*(1d0+2d0*dlog(MST2(M)/MGL**2)))
     .         *mmu/vdq**2/(1d0+epst3)
     .          *((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .            -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .            -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .           +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .            -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .           *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .            -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .  *((CCT(J,K,2)*UT(L,1,1)*U(I,2,1)
     .         -CCT(J,K,2)*UT(L,1,2)*U(I,2,2)
     .         -CCT(J,K,1)*UT(L,1,2)*U(I,2,1)
     .         -CCT(J,K,1)*UT(L,1,1)*U(I,2,2))
     .    *((U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))*fc90(z,yt,zt,xt)
     .       -(V(I,1,1)*U(J,2,1)-V(I,1,2)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt))
     .   +(CCT(J,K,1)*UT(L,1,1)*U(I,2,1)
     .         -CCT(J,K,1)*UT(L,1,2)*U(I,2,2)
     .         +CCT(J,K,2)*UT(L,1,2)*U(I,2,1)
     .         +CCT(J,K,2)*UT(L,1,1)*U(I,2,2))
     .    *((U(I,2,2)*V(J,1,1)+U(I,2,1)*V(J,1,2))*fc90(z,yt,zt,xt)
     .       -(V(I,1,2)*U(J,2,1)+V(I,1,1)*U(J,2,2))
     .        *dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,zt,xt)))
       z=MST2(L)/MCH2(I)
       yt=MST2(J)/MCH2(I)
       zt=MST2(K)/MCH2(I)
       aux=aux-asf(dsqrt(MST2(M)))/6d0/Pi*MST2(M)/MCH2(I)/MW**2
     . *(CCT(I,K,1)*CCT(I,L,1)+CCT(I,K,2)*CCT(I,L,2))*fc50(z,yt,zt)
     .      *((UT(K,1,1)*UT(J,1,1)+UT(K,1,2)*UT(J,1,2))
     .       *((UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .       +(UT(K,1,2)*UT(J,1,1)-UT(K,1,1)*UT(J,1,2))
     .       *(-(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .      +((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,1)*UT(L,1,1)+UT(J,1,2)*UT(L,1,2))
     .      +((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,2)*UT(L,1,1)-UT(J,1,1)*UT(L,1,2)))
     .     -asf(dsqrt(MST2(M)))/6d0/Pi*MST2(M)/MCH2(I)/MW**2
     . *(CCT(I,K,2)*CCT(I,L,1)-CCT(I,K,1)*CCT(I,L,2))*fc50(z,yt,zt)
     .      *((UT(K,1,2)*UT(J,1,1)-UT(K,1,1)*UT(J,1,2))
     .       *((UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .       -(UT(K,1,1)*UT(J,1,1)+UT(K,1,2)*UT(J,1,2))
     .       *(-(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .      +((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,1)*UT(L,1,1)+UT(J,1,2)*UT(L,1,2))
     .      -((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,2)*UT(L,1,1)-UT(J,1,1)*UT(L,1,2)))
       auxe=auxe-asf(dsqrt(MST2(M)))/6d0/Pi*MST2(M)/MCH2(I)/MW**2
     . *(CCT(I,K,1)*CCT(I,L,1)+CCT(I,K,2)*CCT(I,L,2))*fc50(z,yt,zt)
     .      *((UT(K,1,1)*UT(J,1,1)+UT(K,1,2)*UT(J,1,2))
     .       *((UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .       +(UT(K,1,2)*UT(J,1,1)-UT(K,1,1)*UT(J,1,2))
     .       *(-(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .      +((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,1)*UT(L,1,1)+UT(J,1,2)*UT(L,1,2))
     .      +((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,2)*UT(L,1,1)-UT(J,1,1)*UT(L,1,2)))
     .     -asf(dsqrt(MST2(M)))/6d0/Pi*MST2(M)/MCH2(I)/MW**2
     . *(CCT(I,K,2)*CCT(I,L,1)-CCT(I,K,1)*CCT(I,L,2))*fc50(z,yt,zt)
     .      *((UT(K,1,2)*UT(J,1,1)-UT(K,1,1)*UT(J,1,2))
     .       *((UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .       -(UT(K,1,1)*UT(J,1,1)+UT(K,1,2)*UT(J,1,2))
     .       *(-(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .      +((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,1)*UT(L,1,1)+UT(J,1,2)*UT(L,1,2))
     .      -((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,2)*UT(L,1,1)-UT(J,1,1)*UT(L,1,2)))
       aux1=aux1-asf(dsqrt(MST2(M)))/6d0/Pi*MST2(M)/MCH2(I)/MW**2
     . *(CCT(I,K,2)*CCT(I,L,1)-CCT(I,K,1)*CCT(I,L,2))*fc50(z,yt,zt)
     .      *((UT(K,1,1)*UT(J,1,1)+UT(K,1,2)*UT(J,1,2))
     .       *((UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .       +(UT(K,1,2)*UT(J,1,1)-UT(K,1,1)*UT(J,1,2))
     .       *(-(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .      +((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,1)*UT(L,1,1)+UT(J,1,2)*UT(L,1,2))
     .      +((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,2)*UT(L,1,1)-UT(J,1,1)*UT(L,1,2)))
     .     +asf(dsqrt(MST2(M)))/6d0/Pi*MST2(M)/MCH2(I)/MW**2
     . *(CCT(I,K,1)*CCT(I,L,1)+CCT(I,K,2)*CCT(I,L,2))*fc50(z,yt,zt)
     .      *((UT(K,1,2)*UT(J,1,1)-UT(K,1,1)*UT(J,1,2))
     .       *((UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .       -(UT(K,1,1)*UT(J,1,1)+UT(K,1,2)*UT(J,1,2))
     .       *(-(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .      +((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,1)*UT(L,1,1)+UT(J,1,2)*UT(L,1,2))
     .      -((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,2)*UT(L,1,1)-UT(J,1,1)*UT(L,1,2)))
       aux2=aux2-asf(dsqrt(MST2(M)))/6d0/Pi*MST2(M)/MCH2(I)/MW**2
     . *(CCT(I,K,2)*CCT(I,L,1)-CCT(I,K,1)*CCT(I,L,2))*fc50(z,yt,zt)
     .      *((UT(K,1,1)*UT(J,1,1)+UT(K,1,2)*UT(J,1,2))
     .       *((UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .       +(UT(K,1,2)*UT(J,1,1)-UT(K,1,1)*UT(J,1,2))
     .       *(-(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .      +((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,1)*UT(L,1,1)+UT(J,1,2)*UT(L,1,2))
     .      +((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,2)*UT(L,1,1)-UT(J,1,1)*UT(L,1,2)))
     .     +asf(dsqrt(MST2(M)))/6d0/Pi*MST2(M)/MCH2(I)/MW**2
     . *(CCT(I,K,1)*CCT(I,L,1)+CCT(I,K,2)*CCT(I,L,2))*fc50(z,yt,zt)
     .      *((UT(K,1,2)*UT(J,1,1)-UT(K,1,1)*UT(J,1,2))
     .       *((UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .       -(UT(K,1,1)*UT(J,1,1)+UT(K,1,2)*UT(J,1,2))
     .       *(-(UT(J,1,2)*UT(M,1,1)-UT(J,1,1)*UT(M,1,2)
     .         -UT(J,2,2)*UT(M,2,1)+UT(J,2,1)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,1)+UT(M,1,2)*UT(L,1,2)
     .         -UT(M,2,1)*UT(L,2,1)-UT(M,2,2)*UT(L,2,2))
     .        +(UT(J,1,1)*UT(M,1,1)+UT(J,1,2)*UT(M,1,2)
     .         -UT(J,2,1)*UT(M,2,1)-UT(J,2,2)*UT(M,2,2))
     .        *(UT(M,1,1)*UT(L,1,2)-UT(M,1,2)*UT(L,1,1)
     .         -UT(M,2,1)*UT(L,2,2)-UT(M,2,2)*UT(L,2,1)))
     .      +((UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       -(UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,1)*UT(L,1,1)+UT(J,1,2)*UT(L,1,2))
     .      -((UT(K,1,1)*UT(M,1,1)+UT(K,1,2)*UT(M,1,2)
     .        -UT(K,2,1)*UT(M,2,1)-UT(K,2,2)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,1)+UT(M,1,2)*UT(J,1,2)
     .        -UT(M,2,1)*UT(J,2,1)-UT(M,2,2)*UT(J,2,2))
     .       +(UT(K,1,2)*UT(M,1,1)-UT(K,1,1)*UT(M,1,2)
     .        -UT(K,2,2)*UT(M,2,1)+UT(K,2,1)*UT(M,2,2))
     .       *(UT(M,1,1)*UT(J,1,2)-UT(M,1,2)*UT(J,1,1)
     .        -UT(M,2,1)*UT(J,2,2)+UT(M,2,2)*UT(J,2,1)))
     .        *(UT(J,1,2)*UT(L,1,1)-UT(J,1,1)*UT(L,1,2)))
      ENDDO
      ENDDO
      ENDDO
       z=MCH2(J)/MCH2(I)                          ! scharms
       yt=MSU2(1)/MCH2(I)
       xt=MSNE2/MCH2(I)
       aux=aux-asf(dsqrt(MSU2(1)))/3d0/Pi*VVc*MSU2(1)/MCH2(I)**2
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .  *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .  *((V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc90(z,yt,yt,xt)
     .   +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .    *(U(I,2,1)*U(J,2,1)+U(I,2,2)*U(J,2,2))*fc100(z,yt,yt,xt))
     .    -asf(dsqrt(MSU2(1)))/3d0/Pi*VVc*MSU2(1)/MCH2(I)**2
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .  *(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))
     .  *((V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))*fc90(z,yt,yt,xt)
     .   +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .    *(U(I,2,1)*U(J,2,2)-U(I,2,2)*U(J,2,1))*fc100(z,yt,yt,xt))
     . -asf(dsqrt(MSU2(1)))/MW**2/6d0/Pi*VVc*MSU2(1)/MCH2(I)
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .  *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .   *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,yt)
     .     *(U(I,1,1)*U(J,1,1)+U(I,1,2)*U(J,1,2))
     .          -(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc50(z,yt,yt))
     . -asf(dsqrt(MSU2(1)))/MW**2/6d0/Pi*VVc*MSU2(1)/MCH2(I)
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .  *(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))
     .   *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,yt)
     .     *(U(I,1,1)*U(J,1,2)-U(I,1,2)*U(J,1,1))
     .        -(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))*fc50(z,yt,yt))
     . -asf(dsqrt(MSU2(1)))/6d0/Pi*VVc*MSU2(1)/MCH2(I)/MW**2
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .         *(V(I,1,1)**2+V(I,1,2)**2)*fc50(yt,yt,yt)
       auxe=auxe-asf(dsqrt(MSU2(1)))/3d0/Pi*VVc*MSU2(1)/MCH2(I)**2
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .  *(V(I,1,1)**2+V(I,1,2)**2)*(V(J,1,1)**2+V(J,1,2)**2)
     .      *fc90(z,yt,yt,xt)
     . -asf(dsqrt(MSU2(1)))/MW**2/6d0/Pi*VVc*MSU2(1)/MCH2(I)
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .             *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .  *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .    *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,yt)
     .      *(U(I,1,1)*U(J,1,1)+U(I,1,2)*U(J,1,2))
     .      -(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc50(z,yt,yt))
     . -asf(dsqrt(MSU2(1)))/MW**2/6d0/Pi*VVc*MSU2(1)/MCH2(I)
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .             *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .  *(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))
     .    *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,yt)
     .      *(U(I,1,1)*U(J,1,2)-U(I,1,2)*U(J,1,1))
     .      -(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))*fc50(z,yt,yt))
     . -asf(dsqrt(MSU2(1)))/6d0/Pi*VVc*MSU2(1)/MCH2(I)/MW**2
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .          *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .   *(V(I,1,1)**2+V(I,1,2)**2)*fc50(yt,yt,yt)
       aux1=aux1-asf(dsqrt(MSU2(1)))/3d0/Pi*VVc*MSU2(1)/MCH2(I)**2
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .  *(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))
     .  *((V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc90(z,yt,yt,xt)
     .   +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .    *(U(I,2,1)*U(J,2,1)+U(I,2,2)*U(J,2,2))*fc100(z,yt,yt,xt))
     .    +asf(dsqrt(MSU2(1)))/3d0/Pi*VVc*MSU2(1)/MCH2(I)**2
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .  *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .  *((V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))*fc90(z,yt,yt,xt)
     .   +2d0*(mmu/vdq)**2/g2q*dsqrt(MCH2(J)/MCH2(I))
     .    *(U(I,2,1)*U(J,2,2)-U(I,2,2)*U(J,2,1))*fc100(z,yt,yt,xt))
     . -asf(dsqrt(MSU2(1)))/MW**2/6d0/Pi*VVc*MSU2(1)/MCH2(I)
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .  *(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))
     .   *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,yt)
     .     *(U(I,1,1)*U(J,1,1)+U(I,1,2)*U(J,1,2))
     .          -(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc50(z,yt,yt))
     . +asf(dsqrt(MSU2(1)))/MW**2/6d0/Pi*VVc*MSU2(1)/MCH2(I)
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .               *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .  *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .   *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,yt)
     .     *(U(I,1,1)*U(J,1,2)-U(I,1,2)*U(J,1,1))
     .        -(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))*fc50(z,yt,yt))
       aux2=aux2
     . -asf(dsqrt(MSU2(1)))/MW**2/6d0/Pi*VVc*MSU2(1)/MCH2(I)
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .             *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .  *(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))
     .    *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,yt)
     .      *(U(I,1,1)*U(J,1,1)+U(I,1,2)*U(J,1,2))
     .      -(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))*fc50(z,yt,yt))
     . +asf(dsqrt(MSU2(1)))/MW**2/6d0/Pi*VVc*MSU2(1)/MCH2(I)
     .  *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .             *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))**2
     .  *(V(I,1,1)*V(J,1,1)+V(I,1,2)*V(J,1,2))
     .    *(2d0*dsqrt(MCH2(J)/MCH2(I))*fc60(z,yt,yt)
     .      *(U(I,1,1)*U(J,1,2)-U(I,1,2)*U(J,1,1))
     .      -(V(I,1,2)*V(J,1,1)-V(I,1,1)*V(J,1,2))*fc50(z,yt,yt))
       CSCHAR=CSCHAR-asf(dsqrt(MSU2(1)))/3d0/Pi*VVc/(1d0+epst3)
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi*(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                     *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *mmu/vdq**2*MSU2(1)/MCH2(I)**2
     .  *(U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))
     .  *((U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))*fc90(z,yt,yt,xt)
     .       +dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,yt,xt)
     .   *(V(I,1,1)*U(J,2,1)-V(I,1,2)*U(J,2,2)))
     .     -asf(dsqrt(MSU2(1)))/3d0/Pi*VVc/(1d0+epst3)
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi*(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                     *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *mmu/vdq**2*MSU2(1)/MCH2(I)**2
     .  *(U(I,2,2)*V(J,1,1)+U(I,2,1)*V(J,1,2))
     .  *((U(I,2,2)*V(J,1,1)+U(I,2,1)*V(J,1,2))*fc90(z,yt,yt,xt)
     .       +dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,yt,xt)
     .   *(V(I,1,1)*U(J,2,2)+V(I,1,2)*U(J,2,1)))
       CSCHARI=CSCHARI+asf(dsqrt(MSU2(1)))/3d0/Pi*VVc/(1d0+epst3)
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi*(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                     *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *mmu/vdq**2*MSU2(1)/MCH2(I)**2
     .  *(U(I,2,2)*V(J,1,1)+U(I,2,1)*V(J,1,2))
     .  *((U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))*fc90(z,yt,yt,xt)
     .       +dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,yt,xt)
     .   *(V(I,1,1)*U(J,2,1)-V(I,1,2)*U(J,2,2)))
     .     -asf(dsqrt(MSU2(1)))/3d0/Pi*VVc/(1d0+epst3)
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi*(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                     *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *mmu/vdq**2*MSU2(1)/MCH2(I)**2
     .  *(U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))
     .  *((U(I,2,2)*V(J,1,1)+U(I,2,1)*V(J,1,2))*fc90(z,yt,yt,xt)
     .       +dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,yt,xt)
     .   *(V(I,1,1)*U(J,2,2)+V(I,1,2)*U(J,2,1)))
       CPCHAR=CPCHAR+asf(dsqrt(MSU2(1)))/3d0/Pi*VVc/(1d0+epst3)
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi*(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *MSU2(1)/MCH2(I)**2*mmu/vdq**2
     .   *(U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))
     .   *((U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))*fc90(z,yt,yt,xt)
     .       -dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,yt,xt)
     .    *(V(I,1,1)*U(J,2,1)-V(I,1,2)*U(J,2,2)))
     .   +asf(dsqrt(MSU2(1)))/3d0/Pi*VVc/(1d0+epst3)
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi*(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *MSU2(1)/MCH2(I)**2*mmu/vdq**2
     .   *(U(I,2,2)*V(J,1,1)+U(I,2,1)*V(J,1,2))
     .   *((U(I,2,2)*V(J,1,1)+U(I,2,1)*V(J,1,2))*fc90(z,yt,yt,xt)
     .       -dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,yt,xt)
     .    *(V(I,1,1)*U(J,2,2)+V(I,1,2)*U(J,2,1)))
       CPCHARI=CPCHARI-asf(dsqrt(MSU2(1)))/3d0/Pi*VVc/(1d0+epst3)
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi*(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *MSU2(1)/MCH2(I)**2*mmu/vdq**2
     .   *(U(I,2,2)*V(J,1,1)+U(I,2,1)*V(J,1,2))
     .   *((U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))*fc90(z,yt,yt,xt)
     .       -dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,yt,xt)
     .    *(V(I,1,1)*U(J,2,1)-V(I,1,2)*U(J,2,2)))
     .   +asf(dsqrt(MSU2(1)))/3d0/Pi*VVc/(1d0+epst3)
     . *(1d0+asf(dsqrt(MSU2(1)))/4d0/Pi*(1d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *(1d0-asf(dsqrt(MSU2(1)))/4d0/Pi
     .                *(7d0/2d0+2d0*dlog(MSU2(1)/MGL**2)))
     . *MSU2(1)/MCH2(I)**2*mmu/vdq**2
     .   *(U(I,2,1)*V(J,1,1)-U(I,2,2)*V(J,1,2))
     .   *((U(I,2,2)*V(J,1,1)+U(I,2,1)*V(J,1,2))*fc90(z,yt,yt,xt)
     .       -dsqrt(MCH2(J)/MCH2(I))*fc100(z,yt,yt,xt)
     .    *(V(I,1,1)*U(J,2,2)+V(I,1,2)*U(J,2,1)))
      ENDDO
      ENDDO

      C10eCHAR=CACHAR+auxe*MW**2/4d0
      CACHAR=CACHAR+aux*MW**2/4d0
      CSCHAR=CSCHAR*MW**2/2d0/g2q
      CPCHAR=CPCHAR*MW**2/2d0/g2q
      C10eCHARI=CACHARI+aux2*MW**2/4d0
      CACHARI=CACHARI+aux1*MW**2/4d0
      CSCHARI=CSCHARI*MW**2/2d0/g2q
      CPCHARI=CPCHARI*MW**2/2d0/g2q

*         - Higgs penguin from integrated SUSY loops at large tanB [1]
      aux=0d0
      auxe=0d0
      do i=1,5
      aux=aux-((XH(i,1)-XH(i,2)*tanb)*sigRLbs
     .                        -XH(i,4)*(cosb+sinb*tanb)*sigRLbsI)
     . *XH(i,2)*sgn(MH0(i)-M_Bs**2)/(cosb
     . *dsqrt((MH0(i)-M_Bs**2)**2+MH0(i)*WIDTH(i)**2))
      auxe=auxe-((XH(i,1)-XH(i,2)*tanb)*sigRLbsI
     .                         +XH(i,4)*(cosb+sinb*tanb)*sigRLbs)
     . *XH(i,2)*sgn(MH0(i)-M_Bs**2)/(cosb
     . *dsqrt((MH0(i)-M_Bs**2)**2+MH0(i)*WIDTH(i)**2))
      enddo
      CSHP=CSHP+gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*aux
      CSHPI=CSHPI+gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*auxe

      aux=0d0
      auxe=0d0
      do i=1,5
      aux=aux-(-XH(i,4)*(cosb+sinb*tanb)*sigRLbs
     .                -(XH(i,1)-XH(i,2)*tanb)*sigRLbsI)
     .  *XH(i,4)*tanb*sgn(MH0(i)-M_Bs**2)/
     . dsqrt((MH0(i)-M_Bs**2)**2+MH0(i)*WIDTH(i)**2)
      auxe=auxe-((XH(i,1)-XH(i,2)*tanb)*sigRLbs
     .                -XH(i,4)*(cosb+sinb*tanb)*sigRLbsI)
     .  *XH(i,4)*tanb*sgn(MH0(i)-M_Bs**2)/
     . dsqrt((MH0(i)-M_Bs**2)**2+MH0(i)*WIDTH(i)**2)
      enddo
      CPHP=CPHP-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*aux
      CPHPI=CPHPI-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*auxe

*         - Summary
      I=1                              ! 0: SM; 1: NMSSM
      CA=CASM
      CS=0d0
      CP=0d0
      CAI=0d0
      CSI=0d0
      CPI=0d0

      IF(I.eq.1)then
       CA=CA+CAH+CACHAR
       CS=CSH+CSCHAR+CSHP
       CP=CPH+CPCHAR+CPHP
       CAI=CACHARI
       CSI=CSCHARI+CSHPI
       CPI=CPCHARI+CPHPI
      ENDIF

*  Error estimate: 10% of BSM contributions (NLO), 2.2% for SM [18] (non-param+mt)
      DCA=0.022d0*dabs(CASM)
      DCS=0d0
      DCP=0d0
      DCAI=0d0
      DCSI=0d0
      DCPI=0d0

      IF(I.eq.1)then
       DCA=DCA+0.1d0*(dabs(CAH)+dabs(CACHAR))
       DCS=0.1d0*(dabs(CSH)+dabs(CSCHAR)+dabs(CSHP))
       DCP=0.1d0*(dabs(CPH)+dabs(CPCHAR)+dabs(CPHP))
       DCAI=0.1d0*(dabs(CACHARI))
       DCSI=0.1d0*(dabs(CSCHARI)+dabs(CSHPI))
       DCPI=0.1d0*(dabs(CPCHARI)+dabs(CPHPI))
      ENDIF

*  Error estimate: 10% of BSM contributions (NLO), 2.2% for SM [18] (non-param+mt)
      DCA=0.022d0*dabs(CASM)
      DCS=0d0
      DCP=0d0

      IF(I.eq.1)then
       DCA=DCA+0.1d0*(dabs(CAH)+dabs(CACHAR))
       DCS=0.1d0*(dabs(CSH)+dabs(CSCHAR)+dabs(CSHP))
       DCP=0.1d0*(dabs(CPH)+dabs(CPCHAR)+dabs(CPHP))
      ENDIF

*	 2) Branching ratio BR[Bs --> mu+mu-]
      aux=GF**4*MW**4/32d0/Pi**5*M_Bs**5*tau_Bs*fBs**2    ! prefactor
     .    *VtsVtb2*dsqrt(1d0-4d0*mmu**2/M_Bs**2)

      BRBMUMU=aux*
     .      ((1d0-4d0*mmu**2/M_Bs**2)/(1d0+MS/Mb)**2*(CS**2+CSI**2)
     .       +(CP/(1d0+ms/mb)+2d0*mmu/M_Bs**2*CA)**2
     .       +(CPI/(1d0+ms/mb)+2d0*mmu/M_Bs**2*CAI)**2)

*  Tagged BR, cf. 1204.1737
      ADG=((CP/(1d0+ms/mb)+2d0*mmu/M_Bs**2*CA)**2
     .     -(CPI/(1d0+ms/mb)+2d0*mmu/M_Bs**2*CAI)**2
     .     -(1d0-4d0*mmu**2/M_Bs**2)/(1d0+MS/Mb)**2*CS**2
     .     +(1d0-4d0*mmu**2/M_Bs**2)/(1d0+MS/Mb)**2*CSI**2)
     .    /((1d0-4d0*mmu**2/M_Bs**2)/(1d0+MS/Mb)**2*(CS**2+CSI**2)
     .       +(CP/(1d0+ms/mb)+2d0*mmu/M_Bs**2*CA)**2
     .       +(CPI/(1d0+ms/mb)+2d0*mmu/M_Bs**2*CAI)**2)
      BRBMUMU=BRBMUMU*(1d0+0.088d0*ADG)/(1d0-0.088d0**2)    ! ys=0.088+/-0.014

*  Error bars:
      BRBMUMUMAX=Max(
     . ((CP+DCP)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CA+DCA))**2,
     . ((CP+DCP)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CA-DCA))**2,
     . ((CP-DCP)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CA+DCA))**2,
     . ((CP-DCP)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CA-DCA))**2)
     . +Max(
     . ((CPI+DCPI)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CAI+DCAI))**2,
     . ((CPI+DCPI)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CAI-DCAI))**2,
     . ((CPI-DCPI)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CAI+DCAI))**2,
     . ((CPI-DCPI)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CAI-DCAI))**2)
      BRBMUMUMAX=aux*(1d0+dsqrt(dtau_Bs**2        ! Total Bs width
     .               +(dVtsVtb2/VtsVtb2)**2       ! CKM
     .               +(2d0*dfBs)**2))             ! hadronic parameter
     . *((1d0-4d0*mmu**2/M_Bs**2)/(1d0+MS/Mb)**2*(dabs(CS)+DCS)**2
     .                       +BRBMUMUMAX)         ! higher orders
      BRBMUMUMAX=BRBMUMUMAX*(1d0+0.116d0*ADG+0.028d0* ! ys
     . dabs((ADG*(1d0-0.088d0**2)-2d0*0.088d0)/(1d0-0.088d0**2)))
     . /(1d0-0.088d0**2)

      IF(dabs(CP/(1d0+ms/mb)+2d0*mmu/M_Bs**2*CA).gt.
     .   (DCP/(1d0+ms/mb)+2d0*mmu/M_Bs**2*DCA))THEN
      BRBMUMUMIN=min(
     . ((CP+DCP)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CA+DCA))**2,
     . ((CP+DCP)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CA-DCA))**2,
     . ((CP-DCP)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CA+DCA))**2,
     . ((CP-DCP)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CA-DCA))**2)
       ELSE
       BRBMUMUMIN=-(dabs(CP/(1d0+ms/mb)+2d0*mmu/M_Bs**2*CA)
     .  -(DCP/(1d0+ms/mb)+2d0*mmu/M_Bs**2*DCA))**2
       ENDIF
      IF(dabs(CPI/(1d0+ms/mb)+2d0*mmu/M_Bs**2*CAI).gt.
     .   (DCPI/(1d0+ms/mb)+2d0*mmu/M_Bs**2*DCAI))THEN
      BRBMUMUMIN=BRBMUMUMIN+min(
     . ((CPI+DCPI)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CAI+DCAI))**2,
     . ((CPI+DCPI)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CAI-DCAI))**2,
     . ((CPI-DCPI)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CAI+DCAI))**2,
     . ((CPI-DCPI)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CAI-DCAI))**2)
      ELSE
      BRBMUMUMIN=BRBMUMUMIN-
     . (dabs(CPI/(1d0+ms/mb)+2d0*mmu/M_Bs**2*CAI)
     .  -(DCPI/(1d0+ms/mb)+2d0*mmu/M_Bs**2*DCAI))**2
      ENDIF
      IF(dabs(CS).gt.DCS)THEN
      BRBMUMUMIN=BRBMUMUMIN
     . +(1d0-4d0*mmu**2/M_Bs**2)/(1d0+MS/Mb)**2*(dabs(CS)-DCS)**2
      ELSE
      BRBMUMUMIN=BRBMUMUMIN-(1d0-4d0*mmu**2/M_Bs**2)
     .  /(1d0+MS/Mb)**2*(dabs(CS)-DCS)**2
      ENDIF
      IF(dabs(CSI).gt.DCSI)THEN
      BRBMUMUMIN=BRBMUMUMIN
     . +(1d0-4d0*mmu**2/M_Bs**2)/(1d0+MS/Mb)**2*(dabs(CSI)-DCSI)**2
      ELSE
      BRBMUMUMIN=Max(0d0,BRBMUMUMIN-(1d0-4d0*mmu**2/M_Bs**2)
     .  /(1d0+MS/Mb)**2*(dabs(CSI)-DCSI)**2)
      ENDIF
      BRBMUMUMIN=aux*(1d0-dsqrt(dtau_Bs**2    ! Total Bs width
     .               +(dVtsVtb2/VtsVtb2)**2       ! CKM
     .               +(2d0*dfBs)**2))             ! hadronic parameter
     .               *BRBMUMUMIN                  ! higher orders
      BRBMUMUMIN=Max(0d0,BRBMUMUMIN*(1d0+0.116d0*ADG-0.028d0* ! ys
     . dabs((ADG*(1d0-0.088d0**2)-2d0*0.088d0)/(1d0-0.088d0**2)))
     . /(1d0-0.088d0**2))
!      print*,'BRBmumu',BRBMUMUMIN,BRBMUMU,BRBMUMUMAX

*  Comparison with experimental data (source [1',7'])
      prob(35)=0d0

      IF(BRBMUMUmin.GE.BRBMUMUexpMax)
     .     PROB(35)=BRBMUMUmin/BRBMUMUexpMax-1d0
      IF(BRBMUMUmax.LE.BRBMUMUexpMin)
     .     PROB(35)=BRBMUMUmax/BRBMUMUexpMin-1d0

!      csqb=csqb+4.d0*(BRBMUMU-(BRBMUMUexpmin+BRBMUMUexpMax)/2d0)**2
!     c /((BRBMUMUMax-BRBMUMUmin)**2+(BRBMUMUexpMax-BRBMUMUexpmin)**2)

*	 3) Branching ratio BR[Bd --> mu+mu-]

      CSHP=-mmu*tanb**2/4d0/MW**2*(fh30(xt,z)            ! EW Higgs Penguin
     .      +asc0/4d0/Pi*(fh31(xt,z)
     .                       +8d0*xt*fh30p(xt,z)*dlog((sc0/MT0)**2)))
      CPHP=-CSHP

      aux=0d0
      do i=1,5
      aux=aux+XH(i,4)**2*MH0(i)
     .          *sgn(MH0(i)-M_Bd**2)/
     . dsqrt((MH0(i)-M_Bd**2)**2+MH0(i)*WIDTH(i)**2)
      enddo
      CPHP=CPHP*aux

      aux=0d0
      do i=1,5
      aux=aux+XH(i,2)**2*MH0(i)
     .         *sgn(MH0(i)-M_Bd**2)/
     . dsqrt((MH0(i)-M_Bd**2)**2+MH0(i)*WIDTH(i)**2)
      enddo
      CSHP=CSHP*aux

*         - Higgs penguin from integrated SUSY loops at large tanB [1]
      aux=0d0
      auxe=0d0
      do i=1,5
      aux=aux-((XH(i,1)-XH(i,2)*tanb)*sigRLbd
     .                        -XH(i,4)*(cosb+sinb*tanb)*sigRLbdI)
     . *XH(i,2)*sgn(MH0(i)-M_Bd**2)/(cosb
     . *dsqrt((MH0(i)-M_Bd**2)**2+MH0(i)*WIDTH(i)**2))
      auxe=auxe-((XH(i,1)-XH(i,2)*tanb)*sigRLbdI
     .                         +XH(i,4)*(cosb+sinb*tanb)*sigRLbd)
     . *XH(i,2)*sgn(MH0(i)-M_Bd**2)/(cosb
     . *dsqrt((MH0(i)-M_Bd**2)**2+MH0(i)*WIDTH(i)**2))
      enddo
      CSHP=CSHP+gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*aux
      CSHPI=CSHPI+gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*auxe

      aux=0d0
      auxe=0d0
      do i=1,5
      aux=aux-(-XH(i,4)*(cosb+sinb*tanb)*sigRLbd
     .                -(XH(i,1)-XH(i,2)*tanb)*sigRLbdI)
     .  *XH(i,4)*tanb*sgn(MH0(i)-M_Bd**2)/
     . dsqrt((MH0(i)-M_Bd**2)**2+MH0(i)*WIDTH(i)**2)
      auxe=auxe-((XH(i,1)-XH(i,2)*tanb)*sigRLbd
     .                -XH(i,4)*(cosb+sinb*tanb)*sigRLbdI)
     .  *XH(i,4)*tanb*sgn(MH0(i)-M_Bd**2)/
     . dsqrt((MH0(i)-M_Bd**2)**2+MH0(i)*WIDTH(i)**2)
      enddo
      CPHP=CPHP-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*aux
      CPHPI=CPHPI-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*auxe

*         - Summary
      I=1                              ! 0: SM; 1: NMSSM
      CA=CASM
      CS=0d0
      CP=0d0
      CAI=0d0
      CSI=0d0
      CPI=0d0

      IF(I.eq.1)then
       CA=CA+CAH+CACHAR
       CS=CSH+CSCHAR+CSHP
       CP=CPH+CPCHAR+CPHP
       CAI=CACHARI
       CSI=CSCHARI+CSHPI
       CPI=CPCHARI+CPHPI
      ENDIF

*  Error estimate: 10% of BSM contributions (NLO), 2.2% for SM [18] (non-param+mt)
      DCA=0.022d0*dabs(CASM)
      DCS=0d0
      DCP=0d0
      DCAI=0d0
      DCSI=0d0
      DCPI=0d0

      IF(I.eq.1)then
       DCA=DCA+0.1d0*(dabs(CAH)+dabs(CACHAR))
       DCS=0.1d0*(dabs(CSH)+dabs(CSCHAR)+dabs(CSHP))
       DCP=0.1d0*(dabs(CPH)+dabs(CPCHAR)+dabs(CPHP))
       DCAI=0.1d0*(dabs(CACHARI))
       DCSI=0.1d0*(dabs(CSCHARI)+dabs(CSHPI))
       DCPI=0.1d0*(dabs(CPCHARI)+dabs(CPHPI))
      ENDIF

      aux=GF**4*MW**4/32d0/Pi**5*M_Bd**5*tau_Bd*fBd**2    ! prefactor
     .    *VtdVtb2*dsqrt(1d0-4d0*mmu**2/M_Bd**2)

      BRBdMUMU=aux*
     .      ((1d0-4d0*mmu**2/M_Bd**2)/(1d0+Md/Mb)**2*(CS**2+CSI**2)
     .       +(CP/(1d0+md/mb)+2d0*mmu/M_Bd**2*CA)**2
     .       +(CPI/(1d0+md/mb)+2d0*mmu/M_Bd**2*CAI)**2)

*  Error bars:
      BRBdMUMUMAX=Max(
     . ((CP+DCP)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CA+DCA))**2,
     . ((CP+DCP)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CA-DCA))**2,
     . ((CP-DCP)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CA+DCA))**2,
     . ((CP-DCP)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CA-DCA))**2)
     . +Max(
     . ((CPI+DCPI)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CAI+DCAI))**2,
     . ((CPI+DCPI)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CAI-DCAI))**2,
     . ((CPI-DCPI)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CAI+DCAI))**2,
     . ((CPI-DCPI)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CAI-DCAI))**2)
      BRBdMUMUMAX=aux*(1d0+dsqrt(dtau_Bd**2       ! Total Bd width
     .               +(dVtdVtb2/VtdVtb2)**2       ! CKM
     .               +(2d0*dfBd)**2))             ! hadronic parameter
     . *((1d0-4d0*mmu**2/M_Bd**2)/(1d0+Md/Mb)**2*(dabs(CS)+DCS)**2
     .                       +BRBdMUMUMAX)        ! higher orders

      IF(dabs(CP/(1d0+md/mb)+2d0*mmu/M_Bd**2*CA).gt.
     .   (DCP/(1d0+md/mb)+2d0*mmu/M_Bd**2*DCA))THEN
      BRBdMUMUMIN=min(
     . ((CP+DCP)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CA+DCA))**2,
     . ((CP+DCP)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CA-DCA))**2,
     . ((CP-DCP)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CA+DCA))**2,
     . ((CP-DCP)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CA-DCA))**2)
       ELSE
       BRBdMUMUMIN=-(dabs(CP/(1d0+md/mb)+2d0*mmu/M_Bd**2*CA)
     .  -(DCP/(1d0+md/mb)+2d0*mmu/M_Bd**2*DCA))**2
       ENDIF
      IF(dabs(CPI/(1d0+md/mb)+2d0*mmu/M_Bd**2*CAI).gt.
     .   (DCPI/(1d0+md/mb)+2d0*mmu/M_Bd**2*DCAI))THEN
      BRBdMUMUMIN=BRBdMUMUMIN+min(
     . ((CPI+DCPI)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CAI+DCAI))**2,
     . ((CPI+DCPI)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CAI-DCAI))**2,
     . ((CPI-DCPI)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CAI+DCAI))**2,
     . ((CPI-DCPI)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CAI-DCAI))**2)
      ELSE
      BRBdMUMUMIN=BRBdMUMUMIN-
     . (dabs(CPI/(1d0+md/mb)+2d0*mmu/M_Bd**2*CAI)
     .  -(DCPI/(1d0+md/mb)+2d0*mmu/M_Bd**2*DCAI))**2
      ENDIF
      IF(dabs(CS).gt.DCS)THEN
      BRBdMUMUMIN=BRBdMUMUMIN
     . +(1d0-4d0*mmu**2/M_Bd**2)/(1d0+md/Mb)**2*(dabs(CS)-DCS)**2
      ELSE
      BRBdMUMUMIN=BRBdMUMUMIN-(1d0-4d0*mmu**2/M_Bd**2)
     .  /(1d0+md/Mb)**2*(dabs(CS)-DCS)**2
      ENDIF
      IF(dabs(CSI).gt.DCSI)THEN
      BRBdMUMUMIN=BRBdMUMUMIN
     . +(1d0-4d0*mmu**2/M_Bd**2)/(1d0+md/Mb)**2*(dabs(CSI)-DCSI)**2
      ELSE
      BRBdMUMUMIN=Max(0d0,BRBdMUMUMIN-(1d0-4d0*mmu**2/M_Bd**2)
     .  /(1d0+md/Mb)**2*(dabs(CSI)-DCSI)**2)
      ENDIF
      BRBdMUMUMIN=aux*(1d0-dsqrt(dtau_Bd**2       ! Total Bd width
     .               +(dVtdVtb2/VtdVtb2)**2       ! CKM
     .               +(2d0*dfBd)**2))             ! hadronic parameter
     .               *BRBdMUMUMIN                 ! higher orders
!      print*,'BRBdmumu',BRBdMUMUMIN,BRBdMUMU,BRBdMUMUMAX

*  Comparison with experimental data (source [1',7'])
      PROB(56)=0d0
      IF(BRBdMUMUmin.GE.BRBdMUMUexpMax)
     .     PROB(56)=BRBdMUMUmin/BRBdMUMUexpMax-1d0
      IF(BRBdMUMUmax.LE.BRBdMUMUexpMin)
     .     PROB(56)=BRBdMUMUmax/BRBdMUMUexpMin-1d0

!      csqb=csqb+4.d0*(BRBdMUMU-(BRBdMUMUexpmin+BRBdMUMUexpMax)/2d0)**2
!     c /((BRBdMUMUMax-BRBdMUMUmin)**2+(BRBdMUMUexpMax-BRBdMUMUexpmin)**2)


      return
      end
