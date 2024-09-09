      SUBROUTINE BOTTOMONIUM(PROB)

***********************************************************************
*   Subroutine to check bottomonium physics driven by a light CP-odd Higgs:
*      PROB(38) =/= 0  excluded by Upsilon(1S) -> H/A gamma
*      PROB(39) =/= 0  excluded by eta_b(1S) mass measurement
***********************************************************************

      IMPLICIT NONE


      DOUBLE PRECISION PROB(*),PI,SQR2
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION ALEM0,MY,BRYMUMU,C,ZZ,AP,GYEE
      DOUBLE PRECISION ALPHAS,ALSMY,DELTA,UPSILONTAU,UPSILONMU
      DOUBLE PRECISION M0,RETA,GAM2,XX,YY,F,YMAX,FMAX,D,MEMAX,UU,VV
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5),MBQM,MA,MH,RMAX,CB(5),CL(5)
      DOUBLE PRECISION DDCOS,FLYDIST,UPSILONHAD,UPSILONINV

      COMMON/ALEM0/ALEM0
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/CB/CB,CL

      PI=4d0*DATAN(1d0)
      SQR2=DSQRT(2d0)

* Light H Physics

      MH=SMASS(1)

* Test on Upsilon(1S) -> H gamma (from CLEO + BABAR)

      MY=9.46d0 ! Upsilon(1S) mass
      MBQM=4.9d0 ! b quark mass in quark models
      ALSMY=ALPHAS(MY,2) ! alpha_s at MY, 2 loop calculation
      BRYMUMU=2.48d-2 ! BR(Upsilon(1S) -> mu mu)

      IF(MH.LT.MY)THEN

c     Check whether Higgs long-lived

      FLYDIST=(1d0+(MH/MY)**2)/(1d0-(MH/MY)**2)*1.97327d-16/WIDTH(1)

      IF(FLYDIST.lt.1d0)then              ! Higgs short-lived

      ZZ=1d0-MH**2/MY**2 ! energy fraction of the photon
      AP=6d0*ZZ**.2d0 ! Nason function for QCD corrections
      C=1d0+4d0*ALSMY/(3d0*PI)*(4d0-AP) ! QCD corrections
      DELTA=1.2d0**2/MBQM**2 ! function for rel. corrections
      C=C* ! relativistic corrections (for MH<~8.8 GeV)
     .  (MY**2-MH**2)**2/(4d0*MBQM**2-MH**2)**2*(1d0-
     .  DELTA/3d0*(36d0*MBQM**2+MH**2)/(4d0*MBQM**2-MH**2))

      C=MAX(C,1d-6)

c     Limits on leptonic decays (CLEO+BABAR)
      RMAX=0d0
      RMAX=SQR2*PI*ALEM0*UPSILONTAU(MH)/(GF*MBQM**2*ZZ*C*BRYMUMU)
      IF(RMAX.NE.0d0)PROB(38)=DDIM(CB(1)**2*BRLL(1)/RMAX,1d0)

      RMAX=0d0
      RMAX=SQR2*PI*ALEM0*UPSILONMU(MH)/(GF*MBQM**2*ZZ*C*BRYMUMU)
      IF(RMAX.NE.0d0)PROB(38)=PROB(38)+DDIM(CB(1)**2*BRMM(1)/RMAX,1d0)

c     Limits on hadronic decays (BABAR)
      RMAX=SQR2*PI*ALEM0*UPSILONHAD(MH)/(GF*MBQM**2*ZZ*C*BRYMUMU)
      IF(RMAX.NE.0d0)PROB(38)=PROB(38)+DDIM(CB(1)**2*BRJJ(1)/RMAX,1d0)

c     Limits on invisible decays (CLEO+BABAR)
      RMAX=SQR2*PI*ALEM0*UPSILONINV(MH)/(GF*MBQM**2*ZZ*C*BRYMUMU)
      IF(RMAX.NE.0d0)PROB(38)=PROB(38)+DDIM(CB(1)**2*BRNEU(1,1,1)
     .                                                      /RMAX,1d0)

      ELSE                               ! Higgs long-lived

c     Limits on invisible decays (CLEO+BABAR)
      RMAX=SQR2*PI*ALEM0*UPSILONINV(MH)/(GF*MBQM**2*ZZ*C*BRYMUMU)
      IF(RMAX.NE.0d0)PROB(38)=PROB(38)+DDIM(CB(1)**2/RMAX,1d0)

      ENDIF

      ENDIF

* Light A Physics

      MA=PMASS(1)

* Test on Upsilon(1S) -> A gamma (from CLEO + BABAR)

      IF(MA.LT.MY)THEN

      FLYDIST=(1d0+(MA/MY)**2)/(1d0-(MA/MY)**2)*1.97327d-16/WIDTH(4)

      IF(FLYDIST.lt.1d0)then              ! Higgs short-lived

      ZZ=1d0-MA**2/MY**2 ! energy fraction of the photon
      AP=6d0*ZZ**.2d0 ! Nason function for QCD corrections
      C=1d0+4d0*ALSMY/(3d0*PI)*(4d0-AP) ! QCD corrections
      DELTA=1.2d0**2/MBQM**2 ! function for rel. corrections
      C=C* ! relativistic corrections (for MA<~8.8 GeV)
     .  (MY**2-MA**2)**2/(4d0*MBQM**2-MA**2)**2*(1d0-
     .  DELTA/3d0*(36d0*MBQM**2+MA**2)/(4d0*MBQM**2-MA**2))

      C=MAX(C,1d-6)

c     Limits on leptonic decays (CLEO+BABAR)
      RMAX=0d0
      RMAX=SQR2*PI*ALEM0*UPSILONTAU(MA)/(GF*MBQM**2*ZZ*C*BRYMUMU)
      IF(RMAX.NE.0d0)PROB(38)=PROB(38)+DDIM(CB(4)**2*BRLL(4)/RMAX,1d0)

      RMAX=0d0
      RMAX=SQR2*PI*ALEM0*UPSILONMU(MA)/(GF*MBQM**2*ZZ*C*BRYMUMU)
      IF(RMAX.NE.0d0)PROB(38)=PROB(38)+DDIM(CB(4)**2*BRMM(4)/RMAX,1d0)

c     Limits on hadronic decays (BABAR)
      RMAX=SQR2*PI*ALEM0*UPSILONHAD(MA)/(GF*MBQM**2*ZZ*C*BRYMUMU)
      IF(RMAX.NE.0d0)PROB(38)=PROB(38)+DDIM(CB(4)**2*BRJJ(4)/RMAX,1d0)

c     Limits on invisible decays (CLEO+BABAR)
      RMAX=SQR2*PI*ALEM0*UPSILONINV(MA)/(GF*MBQM**2*ZZ*C*BRYMUMU)
      IF(RMAX.NE.0d0)PROB(38)=PROB(38)+DDIM(CB(4)**2*BRNEU(4,1,1)
     .                                                      /RMAX,1d0)

      ELSE                               ! Higgs long-lived

c     Limits on invisible decays (CLEO+BABAR)
      RMAX=SQR2*PI*ALEM0*UPSILONINV(MA)/(GF*MBQM**2*ZZ*C*BRYMUMU)
      IF(RMAX.NE.0d0)PROB(38)=PROB(38)+DDIM(CB(4)**2/RMAX,1d0)

      ENDIF

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
       YY=UU*DDCOS(DACOS(VV)/3d0 + 4d0*PI/3d0)
       F=XX*YY*(1d0+GAM2/(XX+YY)**2)
       FMAX=MAX(FMAX,F)
      ENDIF

      RMAX=8d0*PI*MZ**2/(3d0*(g1+g2)/2d0*RETA*M0**3)*FMAX
      PROB(39)=DDIM(CB(4)**2/RMAX,1d0)

      END

*********************************************************************

      double precision function UPSILONHAD(x)
      
      implicit none
      double precision x
      UPSILONHAD=1d0
c Limit from Y(1S) -> gamma (A1->had) (1307.5306)
      IF(x.gt.0.5d0.and.x.lt.1.46d0)
     . UPSILONHAD=10d0*dexp(dlog(2d-7)+(x-0.5d0)/0.96d0*dlog(112d0))
      IF(x.gt.1.46d0.and.x.lt.6d0)
     . UPSILONHAD=10d0*dexp(dlog(224d-7)+(x-1.46d0)/4.54d0*dlog(30d0))
      IF(x.gt.6d0.and.x.lt.9d0)UPSILONHAD=
     . 10d0*dexp(dlog(6.72d-3)+(x-6d0)/3d0*dlog(0.3d0/6.72d-3))
c Limit from Y(2,3S) -> gamma (A1->had) (1108.3549)
c  converting to Y(1S) by factor BR(Y(1S)->mu+mu-)/BR(Y(2,3S)->mu+mu-)
      IF(x.gt.0.27d0.and.x.lt.1.2d0)
     . UPSILONHAD=min(UPSILONHAD,2.48d0/2.18d0*
     .        10d0*dexp(dlog(8d-8)+(x-0.27d0)/0.73d0*dlog(4d2)))
      IF(x.gt.1d0.and.x.lt.1.2d0)
     . UPSILONHAD=min(UPSILONHAD,2.48d0/2.18d0*
     .        10d0*dexp(dlog(2.4d-6)+(x-1d0)/0.2d0*dlog(2d-1)))
      IF(x.gt.1.2d0.and.x.lt.4.4d0)
     . UPSILONHAD=min(UPSILONHAD,2.48d0/2.18d0*
     .        10d0*dexp(dlog(4.8d-7)+(x-1.2d0)/3.2d0*dlog(1d1)))
      IF(x.gt.4.4d0.and.x.lt.7d0)
     . UPSILONHAD=min(UPSILONHAD,2.48d0/2.18d0*
     .        10d0*dexp(dlog(4.8d-6)+(x-4.4d0)/2.6d0*dlog(2.5d0)))
      return
      end

*********************************************************************

      double precision function UPSILONINV(x)
      
      implicit none
      double precision x,EG,MY
      MY=9.46d0 ! Upsilon(1S) mass
      UPSILONINV=1d0
c Limit from Y(1S) -> gamma A1inv (CLEO, Phys.Rev.D 51 (1995) 2053)
      EG=(MY**2-x**2)/(2d0*MY)
      IF(EG.gt.1d0.and.EG.lt.1.5d0)UPSILONINV=
     .      10d0*dexp(dlog(8d-5)+(EG-1d0)/0.5d0*dlog(5d0/8d0))
      IF(EG.gt.1.5d0.and.EG.lt.2d0)UPSILONINV=
     .      10d0*dexp(dlog(5d-5)+(EG-1.5d0)/0.5d0*dlog(1d-1))
      IF(EG.gt.2d0.and.EG.lt.4.7d0)UPSILONINV=
     .      10d0*dexp(dlog(5d-6)+(EG-2d0)/2.7d0*dlog(1.5d0/5d0))
c Limit from Y(1S) -> gamma A1inv (BABAR, 1007.4646)
      IF(x.lt.7d0)UPSILONINV=min(UPSILONINV,4d-6)
      IF(x.gt.7d0.and.x.lt.9d0)UPSILONINV=min(UPSILONINV,
     .        10d0*dexp(dlog(4d-7)+(x-7d0)/2d0*dlog(1d1)))
      return
      end

*********************************************************************
*********************************************************************

      SUBROUTINE BSG(PAR,PROB)

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

      INTEGER I,J,K,L,M,P
      DOUBLE PRECISION PAR(*),PROB(*)
      DOUBLE PRECISION Pi,aux,auxe,asf,H2,BB0,BB1,BB0L,sgn,sp2
      DOUBLE PRECISION TANB,COSB,SINB,au,ad,SST,SSB,AAT,AAB,AAS
      DOUBLE PRECISION ST(2),RST(2,2),SB(2),RSB(2,2),MTL,MTR,MBL
      DOUBLE PRECISION SL(2),RSL(2,2),SSL
      DOUBLE PRECISION XT,YT,ZT,mu,Mch2,gg1,gg2,scR,scR2,runmass
      DOUBLE PRECISION asmt,asmh,asmsusy,asc0,asmc,asmb,runmb
      DOUBLE PRECISION scb,sc0,MT0,MTH,mbh,MQU,MD,MS0,MC0,MB0,mmu,Vbdg
      DOUBLE PRECISION VVc,VVu,Vtb2,Vcs2,Vud2,VtdVtb2,VtsVtb2,Vub2,Vbsg
      DOUBLE PRECISION VVtu,VVtuim,VVtud,VVtudim,dVub,dVtdVtb2,dVtsVtb2
      DOUBLE PRECISION BBs,dBBs,BBd,dBBd,fB,dfB,fBs,dfBs,fBd,dfBd,tauB0
      DOUBLE PRECISION M_Bs,tau_Bs,dtau_Bs,M_Bd,tau_Bd,dtau_Bd,MBu,tauB
      DOUBLE PRECISION DMdexpmin,DMdexpMax,DMsexpmin,DMsexpMax,
     .       BRBMUMUexpMax,BRBMUMUexpMin,BRBTAUNUexpmin,
     .       BRBTAUNUexpMax,BRSGexpmin,BRSGexpMax
      DOUBLE PRECISION BRBdMUMUexpMax,BRBdMUMUexpMin
      DOUBLE PRECISION sigRLbs,sigRLbd,sigLRbs,sigLRbd,ADG
      DOUBLE PRECISION MCHH(2),CCT(2,2),CCD(2),etaH,D0B,D2B,etaS
      DOUBLE PRECISION CVLLSM,CVLLHIG,C1SLLHIG,C2SLLHIG,C1LRHIG,C2LRHIG
      DOUBLE PRECISION CVLLCHAR,C1SLLCHAR,C2SLLCHAR,C1LRCHAR,C2LRCHAR
      DOUBLE PRECISION C2LRDPH,C1SLLDPH,CVLL,C1SLL,C2SLL,C1LR,C2LR
      DOUBLE PRECISION PRLH_tb(2),PRLH_ts(2),PLRH_tb(2),PLRH_ts(2)
      DOUBLE PRECISION BVLL,B1SLL,B2SLL,B1LR,B2LR,PVLL,P1SLL,P2SLL
      DOUBLE PRECISION P1LR,P2LR,rh,S0
      DOUBLE PRECISION CcR,CcL,RD_taulSM,RD_taulexpmin,RD_taulexpmax
      DOUBLE PRECISION RDs_taulSM,RDs_taulexpmin,RDs_taulexpmax
      DOUBLE PRECISION C70SM,C80SM,C11SM,C41SM,C71SM,C81SM
      DOUBLE PRECISION C41HIG,C70HIG,C80HIG,dC7HIG,dC8HIG
      DOUBLE PRECISION C7HIG,C8HIG,C71HIG,C81HIG
      DOUBLE PRECISION C7CHARS,C8CHARS,C7CHAR,C8CHAR,C41CHAR,C71CHAR
      DOUBLE PRECISION C81CHAR,C70S0,C80S0,eta0
      DOUBLE PRECISION C70BSM,C80BSM,C71BSM,C81BSM,C41BSM
      DOUBLE PRECISION DC7BSM,DC8BSM,DC7BSM_b,DC8BSM_b
      DOUBLE PRECISION C70,C80,C11,C41,C71,C81,eta
      DOUBLE PRECISION FF1,FF2,FF3,FG1,FG2,FG3,ffh,fgh,esm,eh
      DOUBLE PRECISION echi,esfe,H17,H18,H27,H28,Q11,Q21,Q31,Q41
      DOUBLE PRECISION gg7,gg8,Delt7,Delt8,vdq
      DOUBLE PRECISION aa(8),bb(4),hh(8),h8(4)
      DOUBLE PRECISION m0011(8),m0021(8),m0031(8),m0041(8),m0051(8)
      DOUBLE PRECISION m0061(8),m0071(8),m0081(8),m0034(8),m0044(8)
      DOUBLE PRECISION m0054(8),m0064(8),m0074(8),m0084(8),m1012(8)
      DOUBLE PRECISION m1022(8),m1032(8),m1042(8),m1052(8),m1062(8)
      DOUBLE PRECISION m1072(8),m1082(8),m1112(8),m1122(8),m1132(8)
      DOUBLE PRECISION m1142(8),m1152(8),m1162(8),m1172(8),m1182(8)
      DOUBLE PRECISION C10b,C20b,C30b,C40b,C50b,C60b,C70b,C80b
      DOUBLE PRECISION C11b,C21b,C31b,C41b,C51b,C61b,C71b,C81b
      DOUBLE PRECISION C7EMb,MBp,MB_kin,MC_scb,z,EPSEW,af,afim,bf,bfim
      DOUBLE PRECISION ff11,ff12,ff17,ff18,ff22,ff27,ff28,ff47,ff48
      DOUBLE PRECISION ff77,ff78,ff88,ff17i,ff27i,ff18i,ff28i
      DOUBLE PRECISION K17,K17I,K27,K27I,K37,K37I,K47,K47I,K57,K57I
      DOUBLE PRECISION K67,K67I,K77,K78,K78I,BSGPERT,KC7BSM,KC8BSM
      DOUBLE PRECISION delt,delt2,lndelt,lndeltp,NP17,NP78,NP88,NP77
      DOUBLE PRECISION lambd1,lambd2,rho1,rho2,MCNP,HQET,CCSL,BRSL
      DOUBLE PRECISION dBPERT,BRDGexpmin,BRDGexpMax
      DOUBLE PRECISION CASM,CAH,CSH,CPH,CACHAR,CSCHAR,CPCHAR
      DOUBLE PRECISION CSHP,CPHP,CA,CS,CP,DCA,DCS,DCP
      DOUBLE PRECISION fh20,fh21,fh20p,fh70,fh71,fh70p,fh30,fh31,fh30p
      DOUBLE PRECISION fc41,fc51,fc50,fc81,fc60,fc121,fc131,fc90,fc100
      DOUBLE PRECISION fc40,fc31,fc30
      DOUBLE PRECISION C700SM,C800SM,C90SM,C100SM,C9H,C10H,C9CHAR
      DOUBLE PRECISION C10CHAR,C9eCHAR,C10eCHAR,C900,C1000,C9e00,C10e00
      DOUBLE PRECISION fhD0,fhD1,DeD1,hc30,hc31,qc51,YSM,WSM
      DOUBLE PRECISION C700,C800,C710,C810,RR7,RR8,RR9,RR10
      DOUBLE PRECISION BRBSmm,BRBSmmmin,BRBSmmmax,ALEMMB
      DOUBLE PRECISION BRBSee,BRBSeemin,BRBSeemax,Prefac
      DOUBLE PRECISION mb1S,CQ12,CQ22,IntpropS1,IntpropS2,IntpropS3
      DOUBLE PRECISION IntpropP1,IntpropP2,IntpropP3,IntpropSh1
      DOUBLE PRECISION IntpropSh2,IntpropSh3,IntpropPh1,IntpropPh2
      DOUBLE PRECISION IntpropPh3
      DOUBLE PRECISION BRBShmm,BRBShmmmin,BRBShmmmax
      DOUBLE PRECISION BRBShee,BRBSheemin,BRBSheemax
      DOUBLE PRECISION BRBSllminexp,BRBSllMaxexp,BRBShllminexp,
     . BRBShllMaxexp
      DOUBLE PRECISION BRBptoKpll,BRBptoKpllmin,BRBptoKpllMax
      DOUBLE PRECISION FL_BKll,FL_BKllmin,FL_BKllMax,S4_BKll,S4_BKllmin
      DOUBLE PRECISION S4_BKllMax,S5_BKll,S5_BKllmin,S5_BKllMax
      DOUBLE PRECISION BRBptoKpllminexp,BRBptoKpllMaxexp
      DOUBLE PRECISION CLSM,fh10,fh11,fh61,CLH,CLeCHAR,CLtauCHAR,RRL
      DOUBLE PRECISION BRBpKpnunuexpMax,BRBKsnunuexpMax,BRBXsnunuexpMax
      DOUBLE PRECISION T1,T2,T3,LO4B,VVtc,NLO4B

      DOUBLE PRECISION ALEM0
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION MS,MC,MBNP,MB,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2)
      DOUBLE PRECISION PCOMP(2,2),CMASS
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      DOUBLE PRECISION QSTSB
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
      DOUBLE PRECISION HTQ,HBQ,MTOPQ,MBOTQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION eps0,epst0,epst1,epst2,epst3,epsts,epstb
      DOUBLE PRECISION epsY32,epsY31,epsY23,epsY13,epscs,epscb

      COMMON/ALEM0/ALEM0
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/SMSPEC/MS,MC,MBNP,MB,MT,MTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,N
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/STSBSCALE/QSTSB
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
      COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/EPSCOUP/eps0,epst0,epst1,epst2,epst3,epsts,epstb,
     .               epsY32,epsY31,epsY23,epsY13,epscs,epscb

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
       asmh=asf(CMASS)

*   M_top/bot at the Charged Higgs Mass Scale:
      mth=mt*(asmh/asmt)**(4d0/7d0)/(1d0+4d0/(3d0*pi)*asmt)
      mbh=runmb(CMASS)

*   Susy scale (squark masses):
       asmsusy=asf(dsqrt(QSTSB))

*   m_top(m_top)(MSbar)
      MT0=MT/(1d0+4d0/(3d0*pi)*asmt+11d0/pi**2*asmt**2) !165d0

*	 2) EW gauge couplings at the SUSY scale:
      gg1=dsqrt(g1q)
      gg2=dsqrt(g2q)

*	 3) SUSY parameters:
*   Trig. Functions of Beta
      TANB=PAR(3)
      sinb=tanb/dsqrt(1d0+tanb**2)
      cosb=sinb/tanb
      au=1d0/tanb
      ad=-tanb

*   Sfermions
      AAT=PAR(12)
      AAB=PAR(13)    
      AAS=PAR(13)
  
*   Stop Masses and Mixing Angles
      ST(1)=MST1
      ST(2)=MST2
      SST=DSQRT(1d0-CST**2)
      RST(1,1)=CST
      RST(2,2)=CST
      RST(1,2)=SST
      RST(2,1)=-SST

*   Sbottom Masses and Mixing Angles
      SB(1)=MSB1
      SB(2)=MSB2
      SSB=DSQRT(1d0-CSB**2)
      RSB(1,1)=CSB
      RSB(2,2)=CSB
      RSB(1,2)=SSB
      RSB(2,1)=-SSB

*   Stau Masses and Mixing Angles
      SL(1)=MSL1
      SL(2)=MSL2
      SSL=DSQRT(1d0-CSL**2)
      RSL(1,1)=CSL
      RSL(2,2)=CSL
      RSL(1,2)=SSL
      RSL(2,1)=-SSL

*   Charginos
      mu=PAR(4)
      Mch2=PAR(21)

*	 4) Arguments of loop-functions

!      do i=1,2
!       xsqc(i)=(MUL/MCH(i))**2         !(m_Q/m_ch(i))^2
!       do k=1,2
!        xstc(k,i)=(ST(k)/MCH(i))**2       !(m_St(k)/m_ch(i))^2
!       enddo
!      enddo

!      xQg=(MUL/MGL)**2              !(m_Q/m_gl)^2
!      xBg1=(SB(1)/MGL)**2             !(m_Sb(1)/m_gl)^2
!      xBg2=(SB(2)/MGL)**2             !(m_Sb(2)/m_gl)^2
!      xTg1=(ST(1)/MGL)**2             !(m_St(1)/m_gl)^2
!      xTg2=(ST(2)/MGL)**2             !(m_St(2)/m_gl)^2

!      do i=1,5
!       do k=1,2
!        xTneu(k,i)=(ST(k)/MNEU(i))**2       !(m_St(k)/m_neu(i))^2
!        xBneu(k,i)=(SB(k)/MNEU(i))**2       !(m_Sb(k)/m_neu(i))^2
!       enddo
!      enddo
      tauB=1.638d-12/(6.58211915d-25)            ! [1']
      tauB0=1.519d-12/(6.58211915d-25)            ! [1']
      mBu=5.27925d0                              ! [6']

      M_Bs=5.36677d0                       ! [6']
      M_Bd=5.27958d0                       ! [6']

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

*   Muon mass
      mmu=0.10566d0  !GeV

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
      epst2=asf(scR)/(3d0*pi)*(BB1(0d0,MGL**2,MDL**2,scR2)
     .                        +BB1(0d0,MGL**2,MDR**2,scR2)
     . -2d0*MGL*(AAS-mu*tanbq)*BB0L(MGL**2,MDL**2,MDR**2,scR2))

      epst3=0d0
      do k=1,2
      epst3=epst3+asf(scR)/(3d0*pi)*(BB1(0d0,MGL**2,SB(k)**2,scR2)
     . -2d0*MGL/MB0*RSB(k,1)*RSB(k,2)*BB0(0d0,MGL**2,SB(k)**2,scR2))
      enddo

*	 * Neutralino / squark contribution
      aux=0d0
      do i=1,5
       aux=aux+2d0/H2Q*MNEU(i)                      
     .   *((gg1/3d0*N(i,1)-gg2*N(i,2))/dsqrt(2d0)*N(i,4)
     .                     *BB0(0d0,MNEU(i)**2,MDL**2,scR2)
     .     +DSQRT(2d0)/3d0*gg1*N(i,1)*N(i,4)
     .                     *BB0(0d0,MNEU(i)**2,MDR**2,scR2))
     .  +1d0/2d0*(gg1/3d0*N(i,1)-gg2*N(i,2))**2
     .                     *BB1(0d0,MNEU(i)**2,MDL**2,scR2)
     .  +2d0/9d0*gg1**2*N(i,1)**2
     .                     *BB1(0d0,MNEU(i)**2,MDR**2,scR2)
      enddo
      epst2=epst2+aux/(32d0*pi**2)

      aux=0d0
      do i=1,5
      do k=1,2
      aux=aux+2d0*MNEU(i)/MB0                    
     . *(RSB(k,1)*(gg1/3d0*N(i,1)-gg2*N(i,2))/dsqrt(2d0)
     .                                     +HBQ*RSB(k,2)*N(i,4))
     . *(dsqrt(2d0)*gg1/3d0*RSB(k,2)*N(i,1)+HBQ*RSB(k,1)*N(i,4))
     .                    *BB0(0d0,MNEU(i)**2,SB(k)**2,scR2)
     . +((dsqrt(2d0)*gg1/3d0*RSB(k,2)*N(i,1)+HBQ*RSB(k,1)*N(i,4))**2
     .     +(RSB(k,1)*(gg1/3d0*N(i,1)-gg2*N(i,2))/dsqrt(2d0)
     .                                     +HBQ*RSB(k,2)*N(i,4))**2)
     .                    *BB1(0d0,MNEU(i)**2,SB(k)**2,scR2)
      enddo
      enddo
      epst3=epst3+1d0/(32d0*pi**2)*aux

*	 * Chargino / squark contribution
      aux=0d0
      do i=1,2
      aux=aux-2d0*MCH(i)*gg2/H2Q
     .     *U(i,2)*V(i,1)*BB0(0d0,MCH(i)**2,MUL**2,scR2)
!     .  +(MC/H1Q)**2*V(i,2)**2*BB1(0d0,MCH(i)**2,MUR**2,scR2)
     .  +(gg2**2*V(i,1)**2) !+(MS0/H2Q)**2*U(i,2)**2)
     .                        *BB1(0d0,MCH(i)**2,MUL**2,scR2)
      enddo
      epst1=epst2+aux*Vud2/(32d0*pi**2)
      epst2=epst2+aux*Vcs2/(32d0*pi**2)

      aux=0d0
      do i=1,2
      do k=1,2
      aux=aux-2d0*MCH(i)/H2Q
     .  *RST(k,1)*U(i,2)*(gg2*RST(k,1)*V(i,1)-HTQ*RST(k,2)*V(i,2))
     .                            *BB0(0d0,MCH(i)**2,ST(k)**2,scR2)
     .  +(HBQ**2*U(i,2)**2*RST(k,1)**2
     .         +(gg2*RST(k,1)*V(i,1)-HTQ*RST(k,2)*V(i,2))**2)
     .       *BB1(0d0,MCH(i)**2,ST(k)**2,scR2)
      enddo
      enddo
      epst3=epst3+Vtb2/(32d0*pi**2)*aux

*        2) Corrections to the neutral Higgs couplings
*     * epsilon_Y^(JI) as in eq. (5.1) in [1] (up to a factor yt^2
*   and a factor tanb), with (3.7) for lambda_0^(JI) and (3.53)
*   for the CKM matrix elements in terms of V^eff
      epsY32=0d0
      do i=1,2
      do k=1,2
      epsY32=epsY32
     .  +(gg2*RST(k,1)*V(i,1)-HTQ*RST(k,2)*V(i,2))**2
     .                     *BB1(0d0,MCH(i)**2,ST(k)**2,scR2)
     .  -2d0*MCH(i)/H2Q*U(i,2)*RST(k,1)
     .      *(gg2*RST(k,1)*V(i,1)-HTQ*RST(k,2)*V(i,2))
     .                     *BB0(0d0,MCH(i)**2,ST(k)**2,scR2)
!     .  +(MS0/H2Q)**2*U(i,2)**2*RST(k,1)**2
!     .                     *BB1(0d0,MCH(i)**2,ST(k)**2,scR2)
      enddo
      epsY32=epsY32
     . +VVc*((gg2*V(i,1))**2*BB1(0d0,MCH(i)**2,MUL**2,scR2)
     .    -2d0*gg2*MCH(i)/H2Q*U(i,2)*V(i,1)
     .                      *BB0(0d0,MCH(i)**2,MUL**2,scR2))
!     .    +(MC/H1Q)**2*V(i,2)**2
!     .                      *BB1(0d0,MCH(i)**2,MUR**2,scR2)
!     .    +(MS0/H2Q)**2*U(i,2)**2
!     .                      *BB1(0d0,MCH(i)**2,MUL**2,scR2)
      enddo
      epsY32=epsY32/(32d0*pi**2)
      
      epsY31=0d0
      do i=1,2
      do k=1,2
      epsY31=epsY31
     .   +(gg2*RST(k,1)*V(i,1)-HTQ*RST(k,2)*V(i,2))**2
     .                *BB1(0d0,MCH(i)**2,ST(k)**2,scR2)
     .   -2d0*MCH(i)/H2Q*U(i,2)*RST(k,1)
     .       *(gg2*RST(k,1)*V(i,1)-HTQ*RST(k,2)*V(i,2))
     .                *BB0(0d0,MCH(i)**2,ST(k)**2,scR2)
!     .      +(MD/H2Q)**2*U(i,2)**2*RST(k,1)**2
!     .                *BB1(0d0,MCH(i)**2,ST(k)**2,scR2)
      enddo
      epsY31=epsY31
     . +VVu*((gg2*V(i,1))**2*BB1(0d0,MCH(i)**2,MUL**2,scR2)
     .    -2d0*gg2*MCH(i)/H2Q*U(i,2)*V(i,1)
     .                      *BB0(0d0,MCH(i)**2,MUL**2,scR2))
!     .    +(MQU/H1Q)**2*V(i,2)**2
!     .           *BB1(0d0,MCH(i)**2,MUR**2,scR2)
!     .    +(MD/H2Q)**2*U(i,2)**2
!     .          *BB1(0d0,MCH(i)**2,MUL**2,scR2))
      enddo
      epsY31=epsY31/(32d0*pi**2)

      epsY13=0d0
      do i=1,2
      do k=1,2
      epsY13=epsY13-2d0/H2Q*U(i,2)*MCH(i)*RST(k,1)
     .       *(gg2*RST(k,1)*V(i,1)-HTQ*RST(k,2)*V(i,2))
     .                       *BB0(0d0,MCH(i)**2,ST(k)**2,scR2)
     .  +(gg2*RST(k,1)*V(i,1)-HTQ*RST(k,2)*V(i,2))**2
     .                       *BB1(0d0,MCH(i)**2,ST(k)**2,scR2)
     .  +HBQ/H2Q*MB0*U(i,2)**2*RST(k,1)**2
     .                       *BB1(0d0,MCH(i)**2,ST(k)**2,scR2)
      enddo
      epsY13=epsY13
     . +VVu*(gg2**2*V(i,1)**2*BB1(0d0,MCH(i)**2,MUL**2,scR2)
!     .    +(MQU/H1Q)**2*V(i,2)**2*BB1(0d0,MCH(i)**2,MUR**2,scR2)
     .    -2d0*gg2/H2Q*U(i,2)*V(i,1)*MCH(i)
     .                       *BB0(0d0,MCH(i)**2,MUL**2,scR2)
     .    +HBQ/H2Q*MB0*U(i,2)**2
     .                       *BB1(0d0,MCH(i)**2,MUL**2,scR2))
      enddo
      epsY13=epsY13/(32d0*pi**2)

      epsY23=0d0
      do i=1,2
      do k=1,2
      epsY23=epsY23-2d0/H2Q*U(i,2)*MCH(i)*RST(k,1)
     .       *(gg2*RST(k,1)*V(i,1)-HTQ*RST(k,2)*V(i,2))
     .                        *BB0(0d0,MCH(i)**2,ST(k)**2,scR2)
     .  +(gg2*RST(k,1)*V(i,1)-HTQ*RST(k,2)*V(i,2))**2
     .                        *BB1(0d0,MCH(i)**2,ST(k)**2,scR2)
     .  +HBQ/H2Q*MB0*U(i,2)**2*RST(k,1)**2
     .                        *BB1(0d0,MCH(i)**2,ST(k)**2,scR2)
      enddo
      epsY23=epsY23
     . +VVc*(gg2**2*V(i,1)**2*BB1(0d0,MCH(i)**2,MUL**2,scR2)
!     .   +(MC0/H1Q)**2*V(i,2)**2*BB1(0d0,MCH(i)**2,MUR**2,scR2)
     .   -2d0*gg2/H2Q*U(i,2)*V(i,1)*MCH(i)
     .                        *BB0(0d0,MCH(i)**2,MUL**2,scR2)
     .   +HBQ/H2Q*MB0*U(i,2)**2
     .                        *BB1(0d0,MCH(i)**2,MUL**2,scR2))
      enddo
      epsY23=epsY23/(32d0*pi**2)

*   * Correction to the CKM matrix: epst0 (*tanb), Eq.(3.53) in [1]
      epst0=epst3-Vtb2*epsY31
      eps0=epst3-Vtb2*epsY32

*   * Couplings [X^s_RL]^JI as in eqs. (3.55) and (3.56), but
*       WITHOUT the S dependent mixing angles
      sigRLbs=MB0*epsY32/(H1Q*(1d0+eps0)*(1d0+epst3))
      
      sigRLbd=MB0*epsY31/(H1Q*(1d0+epst0)*(1d0+epst3))

      sigLRbs=MS0*epsY23/(H1Q*(1d0+eps0)*(1d0+epst3))
     .  *(1d0+epst3+(epst2-epst3)*epsY32/epsY23)/(1d0+epst2)

      sigLRbd=MD*epsY13/(H1Q*(1d0+epst0)*(1d0+epst3))
     .  *(1d0+epst3+(epst1-epst3)*epsY31/epsY13)/(1d0+epst1)

*       3) Corrections to the charged Higgs couplings
*          Epsilon_t as in [1] section 5.3
*    * eps_t(s) - gluino/squark contribution
       epsts=0d0
       do i=1,2
       epsts=epsts-asmsusy*2d0/(3d0*pi)*mu/MGL*RST(i,2)**2
     .                      *H2(ST(i)**2/MGL**2,MDL**2/MGL**2)
       enddo

*               - neutralino/chargino/squark contribution
       MTL=dsqrt(ST(1)**2*RST(1,1)**2+ST(2)**2*RST(2,1)**2)
       MTR=dsqrt(ST(1)**2*RST(1,2)**2+ST(2)**2*RST(2,2)**2)
       aux=0d0
       do i=1,2
       do j=1,5
       aux=aux+MNEU(j)/MCH(i)
     . *(U(i,2)*(gg1*N(j,1)+gg2*N(j,2))/dsqrt(2d0)-gg2*U(i,1)*N(j,4))
     . *(gg2*V(i,1)*N(j,3)*H2(MNEU(j)**2/MCH(i)**2,MTL**2/MCH(i)**2)
     .  +2d0*dsqrt(2d0)/3d0*gg1*V(i,1)*N(j,1)
     .                    *H2(MNEU(j)**2/MCH(i)**2,MTR**2/MCH(i)**2)
     .  -(gg1*N(j,1)/3d0-gg2*N(j,2))*V(i,2)/dsqrt(2d0)
     .                    *H2(MNEU(j)**2/MCH(i)**2,MDL**2/MCH(i)**2))
       enddo
       enddo
       epsts=epsts+aux/16d0/Pi**2

*    * eps_t(b) - gluino/squark contribution
       epstb=0d0
       do i=1,2
       do k=1,2
       epstb=epstb-asmsusy*2d0/(3d0*pi)*mu/MGL*(RST(i,2)*RSB(k,1))**2
     .                      *H2(ST(i)**2/MGL**2,SB(k)**2/MGL**2)
       enddo
       enddo

*               - neutralino/squark contribution
       aux=0d0
       do i=1,5
       do j=1,2
       do k=1,2
       aux=aux+N(i,4)*N(i,3)/MNEU(i)*(RST(j,1)*RSB(k,1))**2
     .                 *H2(ST(j)**2/MNEU(i)**2,SB(k)**2/MNEU(i)**2)
       enddo
       enddo
       enddo
        epstb=epstb-HBQ**2/(16d0*pi**2)*AAB*aux

*               - neutralino/chargino/squark contribution
       MBL=dsqrt(SB(1)**2*RST(1,1)**2+ST(2)**2*RST(2,1)**2)
       aux=0d0
       do i=1,2
       do j=1,5
       aux=aux+MNEU(j)/MCH(i)
     . *(U(i,2)*(gg1*N(j,1)+gg2*N(j,2))/dsqrt(2d0)-gg2*U(i,1)*N(j,4))
     . *(gg2*V(i,1)*N(j,3)*H2(MNEU(j)**2/MCH(i)**2,MTL**2/MCH(i)**2)
     .  +2d0*dsqrt(2d0)/3d0*gg1*V(i,2)*N(j,1)
     .                    *H2(MNEU(j)**2/MCH(i)**2,MTR**2/MCH(i)**2)
     .  -(gg1*N(j,1)/3d0-gg2*N(j,2))*V(i,2)/dsqrt(2d0)
     .                    *H2(MNEU(j)**2/MCH(i)**2,MBL**2/MCH(i)**2))
       enddo
       enddo
       epstb=epstb+aux/16d0/Pi**2

*    * eps_c(s) - gluino/squark contribution
       epscs=-asmsusy*2d0/(3d0*pi)*mu/MGL
     .                      *H2(MUL**2/MGL**2,MDL**2/MGL**2)

*               - neutralino/chargino/squark contribution
       aux=0d0
       do i=1,2
       do j=1,5
       aux=aux+MNEU(j)/MCH(i)
     . *(U(i,2)*(gg1*N(j,1)+gg2*N(j,2))/dsqrt(2d0)-gg2*U(i,1)*N(j,4))
     . *(gg2*V(i,1)*N(j,3)*H2(MNEU(j)**2/MCH(i)**2,MUL**2/MCH(i)**2)
     .  +2d0*dsqrt(2d0)/3d0*gg1*V(i,2)*N(j,1)
     .                    *H2(MNEU(j)**2/MCH(i)**2,MUR**2/MCH(i)**2)
     .  -(gg1*N(j,1)/3d0-gg2*N(j,2))*V(i,2)/dsqrt(2d0)
     .                    *H2(MNEU(j)**2/MCH(i)**2,MDL**2/MCH(i)**2))
       enddo
       enddo
       epscs=epscs+aux/16d0/Pi**2

*    * eps_c(b) - gluino/squark contribution
       epscb=0d0
       do k=1,2
       epscb=epscb-asmsusy*2d0/(3d0*pi)*mu/MGL*RSB(k,1)**2
     .                      *H2(MUL**2/MGL**2,SB(k)**2/MGL**2)
       enddo

*               - neutralino/squark contribution
       aux=0d0
       do i=1,5
       do k=1,2
       aux=aux+N(i,4)*N(i,3)/MNEU(i)*RSB(k,1)**2
     .                 *H2(MUL**2/MNEU(i)**2,SB(k)**2/MNEU(i)**2)
       enddo
       enddo
        epscb=epscb-HBQ**2/(16d0*pi**2)*AAB*aux

*               - neutralino/chargino/squark contribution
       aux=0d0
       do i=1,2
       do j=1,5
       aux=aux+MNEU(j)/MCH(i)
     . *(U(i,2)*(gg1*N(j,1)+gg2*N(j,2))/dsqrt(2d0)-gg2*U(i,1)*N(j,4))
     . *(gg2*V(i,1)*N(j,3)*H2(MNEU(j)**2/MCH(i)**2,MUL**2/MCH(i)**2)
     .  +2d0*dsqrt(2d0)/3d0*gg1*V(i,2)*N(j,1)
     .                    *H2(MNEU(j)**2/MCH(i)**2,MUR**2/MCH(i)**2)
     .  -(gg1*N(j,1)/3d0-gg2*N(j,2))*V(i,2)/dsqrt(2d0)
     .                    *H2(MNEU(j)**2/MCH(i)**2,MBL**2/MCH(i)**2))
       enddo
       enddo
       epscb=epscb+aux/16d0/Pi**2


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
*  Coefficients at CMASS
      MCHH(1)=MW
      MCHH(2)=CMASS

      aux=mth/dsqrt(h1q**2+h2q**2)
      PRLH_tb(1)=-aux
      PRLH_tb(2)=aux/tanb*(1d0-epstb*(1d0/tanb+tanb))
      PRLH_ts(1)=-aux
      PRLH_ts(2)=aux/tanb*(1d0-(1d0/tanb+tanb)
     .            *(epsts-epsY32*(epstb-epsts)/(1d0+eps0)))

      aux=mbh/dsqrt(h1q**2+h2q**2)
      PLRH_tb(1)=aux
      PLRH_tb(2)=aux*((1d0/tanb+tanb)/(1d0+epst3)-1d0/tanb)
      aux=runmass(0.095d0,CMASS)/dsqrt(h1q**2+h2q**2)
      PLRH_ts(1)=aux
      PLRH_ts(2)=aux*((1d0/tanb+tanb)/(1d0+epst2)
     .    *(1d0-epsY23/(1d0+epst2)/(1d0+eps0)
     .         -epsY32*(epst2-epst3)/(1d0+epst3)/(1d0+eps0))-1d0/tanb)

      CVLLHIG=-g2q/2d0*PRLH_tb(2)*PRLH_ts(2)
     .                          *mth**2*D0B(MW,MCHH(2),mth,mth)
     . +(PRLH_tb(2)*PRLH_ts(2))**2*D2B(MCHH(2),MCHH(2),mth,mth)/8d0
     . +PRLH_tb(2)*PRLH_ts(2)*PRLH_tb(1)*PRLH_ts(1)
     .                            *D2B(MCHH(1),MCHH(2),mth,mth)/4d0
      CVLLHIG=CVLLHIG/(GF*MW)**2

      C1SLLHIG=0d0
      DO I=1,2
      DO J=1,2
       C1SLLHIG=C1SLLHIG+PLRH_tb(J)*PLRH_tb(I)*PRLH_ts(I)*PRLH_ts(J)
     .                   *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)/2d0
      ENDDO
      ENDDO
      C1SLLHIG=C1SLLHIG/(GF*MW)**2

      aux=(1d0+epst3)/(1d0+epst2)
      C1LRHIG=0d0
      DO I=1,2
      DO J=1,2
       C1LRHIG=C1LRHIG+PRLH_tb(J)*PLRH_tb(I)*PRLH_ts(I)*PLRH_ts(J)
     . *(D2B(MCHH(I),MCHH(J),mth,mth)-aux*D2B(MCHH(I),MCHH(J),mth,0d0))
     .                                             /4d0
      ENDDO
      ENDDO
      C1LRHIG=C1LRHIG/(GF*MW)**2

      aux=(1d0+epst3)/(1d0+epst2)
      C2LRHIG=0d0
      DO I=1,2
       C2LRHIG=C2LRHIG-g2q*PLRH_tb(I)*PLRH_ts(I)/2d0
     .    *(D2B(MCHH(I),MW,mth,mth)
     .      -2d0*aux*D2B(MCHH(I),MW,mth,0d0)
     .      +aux**2*D2B(MCHH(I),MW,0d0,0d0))
      DO J=1,2
       C2LRHIG=C2LRHIG+PLRH_tb(J)*PRLH_tb(I)*PRLH_ts(I)*PLRH_ts(J)
     .                     *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)
      ENDDO
      ENDDO
      C2LRHIG=C2LRHIG/(GF*MW)**2

*  Running to sc0 [24] App.C
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

*          - Chargino / squark Boxes [1] App.(A.4.ii)
*  Coefficients at scR
      do j=1,2
       CCD(J)=HBQ/(1d0+epst3)*U(j,2)                      ! MB0/H2Q
      do k=1,2
       CCT(j,k)=gg2*V(j,1)*RST(k,1)-HTQ*V(j,2)*RST(k,2)   ! MT0/H1Q
      enddo
      enddo
      
      aux=0d0                   !3rd family contribution
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+CCT(j,k)*CCT(i,l)*CCT(i,k)*CCT(j,l)
     .                   *D2B(MCH(i),MCH(j),ST(k),ST(l))
      enddo                     !3rd family/1st-2nd interference
       aux=aux-2d0*g2q*V(j,1)*V(i,1)*CCT(i,k)*CCT(j,k)
     .                   *D2B(MCH(i),MCH(j),ST(k),MUL)
      enddo                     !1st-2nd family contribution
       aux=aux+(g2q*V(j,1)*V(i,1))**2*D2B(MCH(i),MCH(j),MUL,MUL)
      enddo
      enddo
      CVLLCHAR=aux/8d0/(GF*MW)**2

      aux=0d0                   !3rd family contribution
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+CCD(j)*CCD(i)*RST(k,1)*RST(l,1)*CCT(i,k)*CCT(j,l)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),ST(l))
      enddo                     !3rd family/1st-2nd interference
       aux=aux-2d0*dsqrt(g2q)*CCD(i)*V(j,1)*CCD(j)*RST(k,1)*CCT(i,k)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),MUL)
      enddo                     !1st-2nd family contribution
       aux=aux+g2q*V(j,1)*V(i,1)*CCD(j)*CCD(i)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),MUL,MUL)
      enddo
      enddo
      C1SLLCHAR=-aux/4d0/(GF*MW)**2
      C2SLLCHAR=aux/16d0/(GF*MW)**2

      aux=0d0                   !3rd family contribution
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+CCD(j)*CCD(i)*RST(l,1)**2*CCT(j,k)*CCT(i,k)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),ST(l))
      enddo                     !3rd family/1st-2nd interference
       aux=aux
     . -CCD(i)*CCD(j)*(CCT(i,k)*CCT(j,k)+g2q*V(i,1)*V(j,1)*RST(k,1)**2)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),MUL)
      enddo                     !1st-2nd family contribution
       aux=aux+g2q*V(j,1)*V(i,1)*CCD(j)*CCD(i)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),MUL,MUL)
      enddo
      enddo
      C1LRCHAR=-aux/2d0/(GF*MW)**2*MS0/MB0
      
      aux=0d0                   !3rd family contribution
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+CCD(j)**2*RST(k,1)*RST(l,1)*CCT(i,k)*CCT(i,l)
     .                   *D2B(MCH(i),MCH(j),ST(k),ST(l))
      enddo                     !3rd family/1st-2nd interference
       aux=aux-2d0*dsqrt(g2q)*V(i,1)*CCD(j)**2*RST(k,1)*CCT(i,k)
     .                   *D2B(MCH(i),MCH(j),ST(k),MUL)
      enddo                     !1st-2nd family contribution
       aux=aux+g2q*V(i,1)**2*CCD(j)**2*D2B(MCH(i),MCH(j),MUL,MUL)
      enddo
      enddo
      C2LRCHAR=-aux/2d0/(GF*MW)**2*MS0/MB0

*  Running to sc0 [24] App.C
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

*          - Double Penguin contributions [1] Eqs.(6.12)-(6.22)
      aux=0d0
      do i=1,3
      aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)**2
     . *sgn(SMASS(i)**2-M_Bs**2)/
     . dsqrt((SMASS(i)**2-M_Bs**2)**2+(SMASS(i)*WIDTH(i))**2)
      enddo
      do i=1,2
      aux=aux+(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)**2
     . *sgn(PMASS(i)**2-M_Bs**2)/
     . dsqrt((PMASS(i)**2-M_Bs**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo

      C2LRDPH=-(4d0*pi/(GF*MW))**2*sigRLbs*sigLRbs*aux

      aux=0d0
      do i=1,3
      aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)**2
     . *sgn(SMASS(i)**2-M_Bs**2)/
     . dsqrt((SMASS(i)**2-M_Bs**2)**2+(SMASS(i)*WIDTH(i))**2)
      enddo
      do i=1,2
      aux=aux-(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)**2
     . *sgn(PMASS(i)**2-M_Bs**2)/
     . dsqrt((PMASS(i)**2-M_Bs**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo

      C1SLLDPH=-(4d0*pi/(GF*MW))**2*sigRLbs**2*aux/2d0

*          - Summary
      CVLL=CVLLSM 

      I=1                              ! 0: SM; 1: NMSSM

      IF(I.eq.1)then
      CVLL=CVLL+CVLLHIG+CVLLCHAR

      C1SLL=C1SLLHIG+C1SLLCHAR+C1SLLDPH
      C2SLL=C2SLLHIG+C2SLLCHAR

      C1LR=C1LRHIG+C1LRCHAR
      C2LR=C2LRHIG+C2LRCHAR+C2LRDPH
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
      
      DMs=GF**2*MW**2/(24d0*pi**2)*M_Bs*fBs**2*VtsVtb2/(6.58211915d-13)
     .                  *dabs(aux)

*          - Error estimate
*      First: 2 sigma error bars from lattice Bag parameters: 
*      (2sigma, added quadratically)

      DMsMax=(dBBs*0.551d0*CVLL)**2

      auxe=-5d0/8d0*(1.0153d0*eta**(-0.6315d0)
     .                                   -0.0153d0*eta**(0.7184d0)
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(4.8177d0-5.2272d0*eta)
     .              +eta**(0.7184d0)*(0.3371d0+0.0724d0*eta)))*C1SLL
     . -5d0/8d0*(1.9325d0*(eta**(-0.6315d0)-eta**(0.7184d0))
     .      +asmb/4d0/pi*(eta**(-0.6315d0)*(9.1696d0-38.8778d0*eta)
     .            +eta**(0.7184d0)*(42.5021d0-12.7939d0*eta)))*C2SLL
     
      DMsMax=DMsMax+4d0*(auxe*0.041d0*1.42d0)**2

      auxe=-3d0/2d0*(0.0081d0*(eta**(0.7184d0)-eta**(-0.6315d0))
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(0.0531d0+0.0415d0*eta)
     .              -eta**(0.7184d0)*(0.0566d0+0.0380d0*eta)))*C1SLL
     . -3d0/2d0*(1.0153d0*eta**(0.7184d0)-0.0153d0*eta**(-0.6315d0)
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(0.1011d0+0.3083d0*eta)
     .             +eta**(0.7184d0)*(-7.1314d0+6.7219d0*eta)))*C2SLL

      DMsMax=DMsMax+4d0*(auxe*0.095d0*1.42d0)**2

      auxe=-(eta**(3d0/23d0)+asmb/4d0/Pi*(0.9250d0*eta**(-24d0/23d0)
     .          +eta**(3d0/23d0)*(-2.0994d0+1.1744d0*eta)))/2d0*C1LR
     . -(asmb/4d0/Pi*1.3875d0*(eta**(26d0/23d0)-eta**(-24d0/23d0)))
     .                                                     *C2LR/2d0

      IF(aux*auxe.gt.0d0)then
       DMsMax=DMsMax+4d0*(auxe*0.212d0*1.42d0)**2
       DMsmin=DMsMax+4d0*(auxe*0.067d0*1.42d0)**2
      ELSE
       DMsMax=DMsMax+4d0*(auxe*0.067d0*1.42d0)**2
       DMsmin=DMsMax+4d0*(auxe*0.212d0*1.42d0)**2
      ENDIF
      
      auxe=3d0/4d0*(2d0/3d0*(eta**(3d0/23d0)-eta**(-24d0/23d0))
     .      +asmb/4d0/Pi*(eta**(3d0/23d0)*(-11.7329d0+0.7829*eta)
     .           +eta**(-24d0/23d0)*(-5.3048d0+16.2548d0*eta)))*C1LR
     . +3d0/4d0*(eta**(-24d0/23d0)+asmb/4d0/Pi*(eta**(-24d0/23d0)
     .     *(7.9572d0-8.8822d0*eta)+0.9250d0*eta**(26d0/23d0)))*C2LR

      IF(aux*auxe.gt.0d0)then
       DMsMax=DMsMax+4d0*(auxe*0.054d0*1.42d0)**2
       DMsmin=DMsmin+4d0*(auxe*0.073d0*1.42d0)**2
      ELSE
       DMsMax=DMsMax+4d0*(auxe*0.073d0*1.42d0)**2
       DMsmin=DMsmin+4d0*(auxe*0.054d0*1.42d0)**2
      ENDIF

      DMsMax=(GF**2*MW**2/(24d0*pi**2)*M_Bs*fBs**2*VtsVtb2
     .            /6.58211915d-13)**2*DMsMax
      DMsmin=(GF**2*MW**2/(24d0*pi**2)*M_Bs*fBs**2*VtsVtb2
     .            /6.58211915d-13)**2*DMsmin

*      Second: Uncertainty from matching: 1% SM + 30% on each BSM contribution:
*      NB: SM error usually neglected in the literature
      aux=0.01d0*dabs(PVLL*CVLLSM)                               ! 1% SM

      IF(I.eq.1)then                                             ! 30% BSM
       aux=aux+0.3d0*(dabs(PVLL*CVLLHIG+P1SLL*C1SLLHIG+P2SLL*C2SLLHIG
     . +P1LR*C1LRHIG+P2LR*C2LRHIG)+dabs(PVLL*CVLLCHAR+P1SLL*C1SLLCHAR
     . +P2SLL*C2SLLCHAR+P1LR*C1LRCHAR+P2LR*C2LRCHAR)
     . +dabs(P1SLL*C1SLLDPH+P2LR*C2LRDPH))
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

!      csqb=csqb+4.d0*(DMs-(DMsexpmin+DMsexpMax)/2.d0)**2
!     c /((DMsMax-DMsmin)**2+(DMsexpMax-DMsexpmin)**2)

*	 3) Wilson coefficients for DMd at the matching scale
      
*          - Charged Higgs Boxes [1] App.(A.4.i)
*  Coefficients at CMASS
      MCHH(1)=MW
      MCHH(2)=CMASS

      aux=mth/dsqrt(h1q**2+h2q**2)
      PRLH_tb(1)=-aux
      PRLH_tb(2)=aux/tanb*(1d0-epstb*(1d0/tanb+tanb))
      PRLH_ts(1)=-aux
      PRLH_ts(2)=aux/tanb*(1d0-(1d0/tanb+tanb)
     .            *(epsts-epsY31*(epstb-epsts)/(1d0+epst0)))

      aux=mbh/dsqrt(h1q**2+h2q**2)
      PLRH_tb(1)=aux
      PLRH_tb(2)=aux*((1d0/tanb+tanb)/(1d0+epst3)-1d0/tanb)
      aux=runmass(0.005d0,CMASS)/dsqrt(h1q**2+h2q**2)
      PLRH_ts(1)=aux
      PLRH_ts(2)=aux*((1d0/tanb+tanb)/(1d0+epst1)
     .    *(1d0-epsY13/(1d0+epst1)/(1d0+epst0)
     .         -epsY31*(epst1-epst3)/(1d0+epst3)/(1d0+epst0))-1d0/tanb)

      CVLLHIG=-g2q/2d0*PRLH_tb(2)*PRLH_ts(2)
     .                          *mth**2*D0B(MW,MCHH(2),mth,mth)
     . +(PRLH_tb(2)*PRLH_ts(2))**2*D2B(MCHH(2),MCHH(2),mth,mth)/8d0
     . +PRLH_tb(2)*PRLH_ts(2)*PRLH_tb(1)*PRLH_ts(1)
     .                            *D2B(MCHH(1),MCHH(2),mth,mth)/4d0
      CVLLHIG=CVLLHIG/(GF*MW)**2

      C1SLLHIG=0d0
      DO I=1,2
      DO J=1,2
       C1SLLHIG=C1SLLHIG+PLRH_tb(J)*PLRH_tb(I)*PRLH_ts(I)*PRLH_ts(J)
     .                   *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)/2d0
      ENDDO
      ENDDO
      C1SLLHIG=C1SLLHIG/(GF*MW)**2

      C1LRHIG=0d0
      DO I=1,2
      DO J=1,2
       C1LRHIG=C1LRHIG+PRLH_tb(J)*PLRH_tb(I)*PRLH_ts(I)*PLRH_ts(J)
     .                          *D2B(MCHH(I),MCHH(J),mth,mth)/4d0
      ENDDO
      ENDDO
      C1LRHIG=C1LRHIG/(GF*MW)**2

      aux=(1d0+epst3)/(1d0+epst1)
      C2LRHIG=0d0
      DO I=1,2
       C2LRHIG=C2LRHIG-g2q*PLRH_tb(I)*PLRH_ts(I)/2d0
     .    *(D2B(MCHH(I),MW,mth,mth)
     .      -2d0*aux*D2B(MCHH(I),MW,mth,0d0)
     .      +aux**2*D2B(MCHH(I),MW,0d0,0d0))
      DO J=1,2
       C2LRHIG=C2LRHIG+PLRH_tb(J)*PRLH_tb(I)*PRLH_ts(I)*PLRH_ts(J)
     .                     *mth**2*D0B(MCHH(I),MCHH(J),mth,mth)
      ENDDO
      ENDDO
      C2LRHIG=C2LRHIG/(GF*MW)**2

*  Running to sc0 [24] App.C
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

*          - Chargino / squark Boxes [1] App.(A.4.ii)
*  Coefficients at scR
      do j=1,2
       CCD(J)=HBQ/(1d0+epst3)*U(j,2)
      do k=1,2
       CCT(j,k)=gg2*V(j,1)*RST(k,1)-HTQ*V(j,2)*RST(k,2)
      enddo
      enddo
      
      aux=0d0                   !3rd family contribution
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+CCT(j,k)*CCT(i,l)*CCT(i,k)*CCT(j,l)
     .                   *D2B(MCH(i),MCH(j),ST(k),ST(l))
      enddo                     !3rd family/1st-2nd interference
       aux=aux-2d0*g2q*V(j,1)*V(i,1)*CCT(i,k)*CCT(j,k)
     .                   *D2B(MCH(i),MCH(j),ST(k),MUL)
      enddo                     !1st-2nd family contribution
       aux=aux+(g2q*V(j,1)*V(i,1))**2*D2B(MCH(i),MCH(j),MUL,MUL)
      enddo
      enddo
      CVLLCHAR=aux/8d0/(GF*MW)**2

      aux=0d0                   !3rd family contribution
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+CCD(j)*CCD(i)*RST(k,1)*RST(l,1)*CCT(i,k)*CCT(j,l)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),ST(l))
      enddo                     !3rd family/1st-2nd interference
       aux=aux-2d0*dsqrt(g2q)*CCD(i)*V(j,1)*CCD(j)*RST(k,1)*CCT(i,k)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),MUL)
      enddo                     !1st-2nd family contribution
       aux=aux+g2q*V(j,1)*V(i,1)*CCD(j)*CCD(i)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),MUL,MUL)
      enddo
      enddo
      C1SLLCHAR=-aux/4d0/(GF*MW)**2
      C2SLLCHAR=aux/16d0/(GF*MW)**2

      aux=0d0                   !3rd family contribution
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+CCD(j)*CCD(i)*RST(l,1)**2*CCT(j,k)*CCT(i,k)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),ST(l))
      enddo                     !3rd family/1st-2nd interference
       aux=aux
     . -CCD(i)*CCD(j)*(CCT(i,k)*CCT(j,k)+g2q*V(i,1)*V(j,1)*RST(k,1)**2)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),ST(k),MUL)
      enddo                     !1st-2nd family contribution
       aux=aux+g2q*V(j,1)*V(i,1)*CCD(j)*CCD(i)
     .             *MCH(i)*MCH(j)*D0B(MCH(i),MCH(j),MUL,MUL)
      enddo
      enddo
      C1LRCHAR=-aux/2d0/(GF*MW)**2*MD/MB0
      
      aux=0d0                   !3rd family contribution
      do i=1,2
      do j=1,2
      do k=1,2
      do l=1,2
       aux=aux+CCD(j)**2*RST(k,1)*RST(l,1)*CCT(i,k)*CCT(i,l)
     .                   *D2B(MCH(i),MCH(j),ST(k),ST(l))
      enddo                     !3rd family/1st-2nd interference
       aux=aux-2d0*dsqrt(g2q)*V(i,1)*CCD(j)**2*RST(k,1)*CCT(i,k)
     .                   *D2B(MCH(i),MCH(j),ST(k),MUL)
      enddo                     !1st-2nd family contribution
       aux=aux+g2q*V(i,1)**2*CCD(j)**2*D2B(MCH(i),MCH(j),MUL,MUL)
      enddo
      enddo
      C2LRCHAR=-aux/2d0/(GF*MW)**2*MD/MB0

*  Running to sc0 [24] App.C
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

*          - Double Penguin contributions [1] Eqs.(6.12)-(6.22)
      aux=0d0
      do i=1,3
      aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)**2
     . *sgn(SMASS(i)**2-M_Bd**2)/
     . dsqrt((SMASS(i)**2-M_Bd**2)**2+(SMASS(i)*WIDTH(i))**2)
      enddo
      do i=1,2
      aux=aux+(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)**2
     . *sgn(PMASS(i)**2-M_Bd**2)/
     . dsqrt((PMASS(i)**2-M_Bd**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo

      C2LRDPH=-(4d0*pi/(GF*MW))**2*sigRLbd*sigLRbd*aux

      aux=0d0
      do i=1,3
      aux=aux+(SCOMP(i,1)-SCOMP(i,2)*tanb)**2
     . *sgn(SMASS(i)**2-M_Bd**2)/
     . dsqrt((SMASS(i)**2-M_Bd**2)**2+(SMASS(i)*WIDTH(i))**2)
      enddo
      do i=1,2
      aux=aux-(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)**2
     . *sgn(PMASS(i)**2-M_Bd**2)/
     . dsqrt((PMASS(i)**2-M_Bd**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo

      C1SLLDPH=-(4d0*pi/(GF*MW))**2*sigRLbd**2*aux/2d0

*          - Summary
      CVLL=CVLLSM 

      I=1                              ! 0: SM; 1: NMSSM

      IF(I.eq.1)then
      CVLL=CVLL+CVLLHIG+CVLLCHAR

      C1SLL=C1SLLHIG+C1SLLCHAR+C1SLLDPH
      C2SLL=C2SLLHIG+C2SLLCHAR

      C1LR=C1LRHIG+C1LRCHAR
      C2LR=C2LRHIG+C2LRCHAR+C2LRDPH
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
      
      DMd=GF**2*MW**2/(24d0*pi**2)*M_Bd*fBd**2*VtdVtb2/(6.58211915d-13)
     .                  *dabs(aux)

*          - Error estimate
*      First, error bars from uncertainties on lattice Bag parameters: 
*      (2sigma, added quadratically)

      DMdMax=(dBBd*0.551d0*CVLL)**2

      auxe=-5d0/8d0*(1.0153d0*eta**(-0.6315d0)
     .                                   -0.0153d0*eta**(0.7184d0)
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(4.8177d0-5.2272d0*eta)
     .              +eta**(0.7184d0)*(0.3371d0+0.0724d0*eta)))*C1SLL
     . -5d0/8d0*(1.9325d0*(eta**(-0.6315d0)-eta**(0.7184d0))
     .      +asmb/4d0/pi*(eta**(-0.6315d0)*(9.1696d0-38.8778d0*eta)
     .            +eta**(0.7184d0)*(42.5021d0-12.7939d0*eta)))*C2SLL
     
      DMdMax=DMdMax+4d0*(auxe*0.045d0*1.44d0)**2

      auxe=-3d0/2d0*(0.0081d0*(eta**(0.7184d0)-eta**(-0.6315d0))
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(0.0531d0+0.0415d0*eta)
     .              -eta**(0.7184d0)*(0.0566d0+0.0380d0*eta)))*C1SLL
     . -3d0/2d0*(1.0153d0*eta**(0.7184d0)-0.0153d0*eta**(-0.6315d0)
     .       +asmb/4d0/pi*(eta**(-0.6315d0)*(0.1011d0+0.3083d0*eta)
     .             +eta**(0.7184d0)*(-7.1314d0+6.7219d0*eta)))*C2SLL

      DMdMax=DMdMax+4d0*(auxe*0.103d0*1.44d0)**2

      auxe=-(eta**(3d0/23d0)+asmb/4d0/Pi*(0.9250d0*eta**(-24d0/23d0)
     .          +eta**(3d0/23d0)*(-2.0994d0+1.1744d0*eta)))/2d0*C1LR
     . -(asmb/4d0/Pi*1.3875d0*(eta**(26d0/23d0)-eta**(-24d0/23d0)))
     .                                                     *C2LR/2d0

      IF(aux*auxe.gt.0d0)then
       DMdMax=DMdMax+4d0*(auxe*0.204d0*1.44d0)**2
       DMdmin=DMdMax+4d0*(auxe*0.072d0*1.44d0)**2
      ELSE
       DMdMax=DMdMax+4d0*(auxe*0.072d0*1.44d0)**2
       DMdmin=DMdMax+4d0*(auxe*0.204d0*1.44d0)**2
      ENDIF
      
      auxe=3d0/4d0*(2d0/3d0*(eta**(3d0/23d0)-eta**(-24d0/23d0))
     .      +asmb/4d0/Pi*(eta**(3d0/23d0)*(-11.7329d0+0.7829*eta)
     .           +eta**(-24d0/23d0)*(-5.3048d0+16.2548d0*eta)))*C1LR
     . +3d0/4d0*(eta**(-24d0/23d0)+asmb/4d0/Pi*(eta**(-24d0/23d0)
     .     *(7.9572d0-8.8822d0*eta)+0.9250d0*eta**(26d0/23d0)))*C2LR

      IF(aux*auxe.gt.0d0)then
       DMdMax=DMdMax+4d0*(auxe*0.058d0*1.44d0)**2
       DMdmin=DMdmin+4d0*(auxe*0.076d0*1.44d0)**2
      ELSE
       DMdMax=DMdMax+4d0*(auxe*0.076d0*1.44d0)**2
       DMdmin=DMdmin+4d0*(auxe*0.058d0*1.44d0)**2
      ENDIF

      DMdMax=(GF**2*MW**2/(24d0*pi**2)*M_Bd*fBd**2*VtdVtb2
     .            /6.58211915d-13)**2*DMdMax
      DMdmin=(GF**2*MW**2/(24d0*pi**2)*M_Bd*fBd**2*VtdVtb2
     .            /6.58211915d-13)**2*DMdmin

*      Second, uncertainties from matching: 1% SM + 30% on each BSM contribution:
*      NB: SM error usually neglected in the literature
      aux=0.01d0*dabs(PVLL*CVLLSM)                                ! 1% SM

      IF(I.eq.1)then
       aux=aux+0.3d0*(dabs(PVLL*CVLLHIG+P1SLL*C1SLLHIG+P2SLL*C2SLLHIG
     . +P1LR*C1LRHIG+P2LR*C2LRHIG)+dabs(PVLL*CVLLCHAR+P1SLL*C1SLLCHAR
     . +P2SLL*C2SLLCHAR+P1LR*C1LRCHAR+P2LR*C2LRCHAR)
     . +dabs(P1SLL*C1SLLDPH+P2LR*C2LRDPH))
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
*             (0.5015ps-1 < DMd=0.5015ps-1 < 0.5095ps-1)
      prob(34)=0d0

      IF(DMdmin.GE.DMdexpMax)
     .     PROB(34)=DMdmin/DMdexpMax-1d0
      IF(DMdmax.LE.DMdexpmin)
     .     PROB(34)=DMdmax/DMdexpMin-1d0

!      csqb=csqb+4.d0*(DMd-(DMdexpmin+DMdexpMax)/2.d0)**2
!     c /((DMdMax-DMdmin)**2+(DMdexpMax-DMdexpmin)**2)
!       DMd=(P1SLL*C1SLLDPH+P2LR*C2LRDPH)/BBd

*	IV- Charged B decays: transitions mediated by a charged Higgs
*   1) BR[B+ -> tau+ nu_tau] following [20]

      rh=(1d0-(mBu/CMASS)**2*((1d0+tanb**2)/(1d0+epst0)-1d0))**2
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

!      csqb=csqb+4.d0*(BRBtaunu-(BRBtaunuexpmin+BRBtaunuexpMax)/2.d0)**2
!     c/((BRBtaunuMax-BRBtaunumin)**2+(BRBtaunuexpMax-BRBtaunuexpmin)**2)

*   2) RD=BR[B+ -> D tau+ nu_tau]/BR[B+ -> D l+ nu_l] following [23]
      RD_taulSM=0.297d0

      CcR=-runmb(mb)*mtau/CMASS**2*((1d0+tanb**2)/(1d0+eps0)-1d0)

      CcL=+runmass(1.25d0,mb)*mtau/CMASS**2*(1d0-(1d0/tanb+tanb)
     .                     *(epscb-epsY32/(1d0+eps0)*(epscs-epscb)))

      rh=1d0+1.5d0*(CcR+CcL)+(CcR+CcL)**2

      RD_taul=RD_taulSM*rh

*  Error estimate
      RD_taulMax=RD_taul+2d0*0.017d0                                ! 2sigma SM
     .  +0.3d0*(dabs(CcR)+dabs(CcL))*RD_taul*dabs(1.5d0+CcR+CcL)    !30% NP
      RD_taulmin=RD_taul-2d0*0.017d0                                ! 2sigma SM
     .  -0.3d0*(dabs(CcR)+dabs(CcL))*RD_taul*dabs(1.5d0+CcR+CcL)    !30% NP
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

      rh=1d0+0.12d0*(CcR-CcL)+0.05d0*(CcR-CcL)**2

      RDs_taul=RDs_taulSM*rh

*  Error estimate
      RDs_taulMax=RDs_taul+2d0*0.003d0                              ! 2sigma SM
     .  +0.3d0*(dabs(CcR)+dabs(CcL))*RD_taul                        !30% NP
     .               *dabs(0.12d0+0.05d0*(CcR-CcL))
      RDs_taulmin=RDs_taul-2d0*0.003d0                              ! 2sigma SM
     .  -0.3d0*(dabs(CcR)+dabs(CcL))*RD_taul*                       !30% NP
     .               *dabs(0.12d0+0.05d0*(CcR-CcL))
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
      yt=(MTH/CMASS)**2                                          !(m_t/m_H+)^2

      C70HIG=au**2/3d0*ff1(yt)-au*ad*ff2(yt)             ! LO
      C80HIG=au**2/3d0*fg1(yt)-au*ad*fg2(yt)

      C41HIG=au**2*eh(yt)                                ! NLO
      C71HIG=ffh(yt,tanb)-4d0/9d0*C41HIG
      C81HIG=fgh(yt,tanb)-1d0/6d0*C41HIG

      dC7HIG=-tanb*Vtb2*(epsts+(epsts-epstb)*epsY32/(1d0+eps0)  ! large tanB [1]
     .                                  +epst3/tanb/(1d0+epst3)
     . -Vtb2*epst3/(1d0+epst3)*(epsts+(epsts-epstb)*epsY32/(1d0+eps0))
     .                                                       )*ff2(yt)
      dC8HIG=-tanb*Vtb2*(epsts+(epsts-epstb)*epsY32/(1d0+eps0)  ! large tanB [1]
     .                                  +epst3/tanb/(1d0+epst3)
     . -Vtb2*epst3/(1d0+epst3)*(epsts+(epsts-epstb)*epsY32/(1d0+eps0))
     .                                                       )*fg2(yt)

       etaH=asmh/asc0                     ! Evolution from M_Higgs to sc0

       C7HIG=etaH**(16d0/21d0)*(C70HIG+dC7HIG)
     .       +8d0/3d0*(etaH**(2d0/3d0)-etaH**(16d0/21d0))
     .       *(C80HIG+dC8HIG)
       C8HIG=etaH**(2d0/3d0)*(C80HIG+dC8HIG)

*          - Chargino/Squark Contributions [4,5,6]
       do j=1,2
        CCD(J)=U(J,2)*MW/(dsqrt(2d0)*COSB*MCH(J))/(1d0+epst3)
       do k=1,2
        CCT(j,k)=V(j,1)*RST(k,1)
     . *(1d0-asf(ST(k))/(4d0*Pi)*(7d0/3d0+2d0*dlog((ST(K)/MGL)**2)))
     .          -V(j,2)*RST(k,2)*HTQ/dsqrt(g2)
     . *(1d0+asf(ST(k))/(4d0*Pi)*(1d0+2d0*dlog((ST(K)/MGL)**2)))
       enddo
       enddo

       C7CHARS=0d0
       C8CHARS=0d0
      
       do j=1,2                            ! LO Scharm/charginos
        C7CHARS=C7CHARS-VVc*(
     .          2d0/3d0*V(J,1)**2*(MW/MUL)**2*FF1((MUL/MCH(J))**2)
     .   *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))**2
     .         +CCD(J)*V(J,1)*FF3((MUL/MCH(J))**2)
     .   *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))
     .   *(1d0+asf(MUL)/(4d0*Pi)*(1d0+2d0*dlog((MUL/MGL)**2))))
        C8CHARS=C8CHARS-VVc*(
     .          2d0/3d0*V(J,1)**2*(MW/MUL)**2*FG1((MUL/MCH(J))**2)
     .   *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))**2
     .       +CCD(J)*V(J,1)*FG3((MUL/MCH(J))**2)
     .   *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))
     .   *(1d0+asf(MUL)/(4d0*Pi)*(1d0+2d0*dlog((MUL/MGL)**2))))
       enddo
       etaS=asf(MUL)/asc0                  ! Running between MUL and sc0
       C7CHAR=etaS**(16d0/21d0)*C7CHARS+
     .    8d0/3d0*(etaS**(14d0/21d0)-etaS**(16d0/21d0))*C8CHARS
       C8CHAR=etaS**(14d0/21d0)*C8CHARS

       do k=1,2                            ! LO Stop/charginos
       C7CHARS=0d0
       C8CHARS=0d0
       do j=1,2
         C7CHARS=C7CHARS-2d0/3d0*CCT(J,K)**2*(MW/ST(K))**2
     .                                    *FF1((ST(K)/MCH(J))**2)
     .  -CCD(J)*CCT(J,K)*RST(K,1)*FF3((ST(K)/MCH(J))**2)
     . *(1d0+asf(ST(k))/(4d0*Pi)*(1d0+2d0*dlog((ST(K)/MGL)**2)))
         C8CHARS=C8CHARS-2d0/3d0*CCT(J,K)**2*(MW/ST(K))**2
     .                                    *FG1((ST(K)/MCH(J))**2)
     .  -CCD(J)*CCT(J,K)*RST(K,1)*FG3((ST(K)/MCH(J))**2)
     . *(1d0+asf(ST(k))/(4d0*Pi)*(1d0+2d0*dlog((ST(K)/MGL)**2)))
       enddo
       etaS= asf(ST(k))/asc0                ! Running between ST(K) and sc0
       C7CHAR=C7CHAR+etaS**(16d0/21d0)*C7CHARS+
     .    8d0/3d0*(etaS**(14d0/21d0)-etaS**(16d0/21d0))*C8CHARS
       C8CHAR=C8CHAR+etaS**(14d0/21d0)*C8CHARS
       enddo

       C41CHAR=0d0
       C71CHAR=0d0
       C81CHAR=0d0

       do j=1,2                            ! NLO Scharm/chargino
       C41CHAR=C41CHAR+VVc*(
     .            V(J,1)**2*(MW/MUL)**2*echi((MUL/MCH(J))**2)
     .   *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))**2)
       C71CHAR=C71CHAR+VVc*(
     .             V(J,1)**2*(MW/MCH(J))**2*H17((MUL/MCH(J))**2)
     .   *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))**2
     .         -CCD(J)*V(J,1)*H27((MUL/MCH(J))**2)
     .   *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))
     .   *(1d0+asf(MUL)/(4d0*Pi)*(1d0+2d0*dlog((MUL/MGL)**2))))
       C81CHAR=C81CHAR+VVc*(
     .             V(J,1)**2*(MW/MCH(J))**2*H18((MUL/MCH(J))**2)
     .   *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))**2
     .         -CCD(J)*V(J,1)*H28((MUL/MCH(J))**2)
     .   *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))
     .   *(1d0+asf(MUL)/(4d0*Pi)*(1d0+2d0*dlog((MUL/MGL)**2))))
       enddo

       do j=1,2                            ! NLO Scharm quartic coupling
       C41CHAR=C41CHAR+VVc*(V(J,1)**2*(MW/MUL)**2
     .               *esfe((MUL/MCH(J))**2,(MUL/MCH(J))**2)
     .   *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))**2)
       C71CHAR=C71CHAR+VVc/2d0*(V(J,1)**2*(MW/MUL)**2
     .   *(Q11((MUL/MCH(J))**2,(MUL/MCH(J))**2)
     .           -2d0/3d0*Q21((MUL/MCH(J))**2,(MUL/MCH(J))**2))
     .   *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))**2
     .         -CCD(J)*V(J,1)*(MCH(J)/MUL)**2
     .   *(Q31((MUL/MCH(J))**2,(MUL/MCH(J))**2)
     .           -2d0/3d0*Q41((MUL/MCH(J))**2,(MUL/MCH(J))**2))
     .   *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))
     .   *(1d0+asf(MUL)/(4d0*Pi)*(1d0+2d0*dlog((MUL/MGL)**2))))
       C81CHAR=C81CHAR-VVc/2d0*(V(J,1)**2*(MW/MUL)**2
     .                        *Q21((MUL/MCH(J))**2,(MUL/MCH(J))**2)
     .   *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))**2
     .         -CCD(J)*V(J,1)*(MCH(J)/MUL)**2
     .                        *Q41((MUL/MCH(J))**2,(MUL/MCH(J))**2)
     .   *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))
     .   *(1d0+asf(MUL)/(4d0*Pi)*(1d0+2d0*dlog((MUL/MGL)**2))))
       enddo

       do j=1,2                            ! NLO Stop/chargino
       do k=1,2
       C41CHAR=C41CHAR
     .              +CCT(j,k)**2*(MW/ST(K))**2*echi((ST(K)/MCH(J))**2)
       C71CHAR=C71CHAR
     .              +CCT(j,k)**2*(MW/MCH(J))**2*H17((ST(K)/MCH(J))**2)
     .              -CCD(J)*RST(k,1)*CCT(J,K)*H27((ST(K)/MCH(J))**2)
     .   *(1d0+asf(ST(K))/(4d0*Pi)*(1d0+2d0*dlog((ST(K)/MGL)**2)))
       C81CHAR=C81CHAR
     .              +CCT(j,k)**2*(MW/MCH(J))**2*H18((ST(K)/MCH(J))**2)
     .              -CCD(J)*RST(k,1)*CCT(J,K)*H28((ST(K)/MCH(J))**2)
     .   *(1d0+asf(ST(K))/(4d0*Pi)*(1d0+2d0*dlog((ST(K)/MGL)**2)))
       enddo
       enddo

       do j=1,2                            ! NLO Stop quartic coupling
       do k=1,2
       do l=1,2
       do m=1,2
       C41CHAR=C41CHAR+(MW/ST(l))**2*CCT(k,J)*CCT(m,j)
     .                 *(RST(k,1)*RST(l,1)-RST(k,2)*RST(l,2))
     .                 *(RST(l,1)*RST(m,1)-RST(l,2)*RST(m,2))
     .              *esfe((ST(k)/MCH(J))**2,(ST(m)/MCH(J))**2)
       C71CHAR=C71CHAR+(MW/ST(l))**2/2d0
     .                 *(RST(k,1)*RST(l,1)-RST(k,2)*RST(l,2))
     .                 *(RST(l,1)*RST(m,1)-RST(l,2)*RST(m,2))
     . *(CCT(k,J)*CCT(m,j)*(Q11((ST(k)/MCH(J))**2,(ST(m)/MCH(J))**2)
     .               -2d0/3d0*Q21((ST(k)/MCH(J))**2,(ST(m)/MCH(J))**2))
     .  -CCT(k,j)*CCD(j)*RST(m,1)*(MCH(j)/ST(l))**2
     .   *(1d0+asf(ST(K))/(4d0*Pi)*(1d0+2d0*dlog((ST(K)/MGL)**2)))
     .   *(Q31((ST(k)/MCH(J))**2,(ST(m)/MCH(J))**2)
     .               -2d0/3d0*Q41((ST(k)/MCH(J))**2,(ST(m)/MCH(J))**2)))
       C81CHAR=C81CHAR-(MW/ST(l))**2/2d0
     .                 *(RST(k,1)*RST(l,1)-RST(k,2)*RST(l,2))
     .                 *(RST(l,1)*RST(m,1)-RST(l,2)*RST(m,2))
     . *(CCT(k,J)*CCT(m,j)*Q21((ST(k)/MCH(J))**2,(ST(m)/MCH(J))**2)
     .  -CCT(k,j)*CCD(j)*RST(m,1)*(MCH(j)/ST(l))**2
     .   *(1d0+asf(ST(K))/(4d0*Pi)*(1d0+2d0*dlog((ST(K)/MGL)**2)))
     .              *Q41((ST(k)/MCH(J))**2,(ST(m)/MCH(J))**2))
       enddo
       enddo
       enddo
       enddo

*          - Higgs Penguin Contribution [1]

      C70S0=0d0
      C80S0=0d0

      do i=1,3
       aux=(SCOMP(i,1)-SCOMP(i,2)*tanb)
     .   *(SCOMP(i,2)+SCOMP(i,1)*epst3/tanb)/SMASS(i)**2
       aux=1d0/18d0*MW**2/g2*MB0/(H2Q*(1d0+epst3))*sigRLbs*aux

       IF(SMASS(i).ge.mb)THEN                    ! Running to sc0
        eta0=asf(SMASS(i))/asc0
       ELSE
        eta0=asmb/asc0
       ENDIF

       IF(SMASS(i).gt.mt0)then
        C80S0=C80S0-3d0*eta0**(14d0/21d0)*aux
        C70S0=C70S0+eta0**(16d0/21d0)*aux
     . -3d0*8d0/3d0*(eta0**(14d0/21d0)-eta0**(16d0/21d0))*aux
       ELSE
        C80S0=C80S0-3d0*eta0**(14d0/23d0)*aux
        C70S0=C70S0+eta0**(16d0/23d0)*aux
     . -3d0*8d0/3d0*(eta0**(14d0/23d0)-eta0**(16d0/23d0))*aux
       ENDIF
      enddo

      do i=1,2
       aux=-(PCOMP(i,1)*cosb+PCOMP(i,1)*sinb*tanb)
     .       *(PCOMP(i,1)*sinb-PCOMP(i,1)*cosb*epst3/tanb)/PMASS(i)**2
       aux=1d0/18d0*MW**2/g2*MB0/(H2Q*(1d0+epst3))*sigRLbs*aux

       IF(PMASS(i).ge.mb)THEN                    ! Running to sc0
        eta0=asf(PMASS(i))/asc0
       ELSE
        eta0=asmb/asc0
       ENDIF

       IF(PMASS(i).gt.mt0)then
        C80S0=C80S0-3d0*eta0**(14d0/21d0)*aux
        C70S0=C70S0+eta0**(16d0/21d0)*aux
     . -3d0*8d0/3d0*(eta0**(14d0/21d0)-eta0**(16d0/21d0))*aux
       ELSE
        C80S0=C80S0-3d0*eta0**(14d0/23d0)*aux
        C70S0=C70S0+eta0**(16d0/23d0)*aux
     . -3d0*8d0/3d0*(eta0**(14d0/23d0)-eta0**(16d0/23d0))*aux
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
      ELSE                         ! collects all BSM contributions
       C70BSM=C7HIG+C7CHAR+C70S0
       C80BSM=C8HIG+C8CHAR+C80S0
       C41BSM=C41HIG+C41CHAR
       C71BSM=C71HIG+C71CHAR
       C81BSM=C81HIG+C81CHAR
      ENDIF
                 
*  DC7BSM and DC8BSM are (conservative) error estimates:
*  Dominant errors: HIG: Order (alphas**2) ~ 10%,
*  CHAR: Order(alphas**2) ~ 10%, PENGUIN: Order(alphas) ~ 30%
       DC7BSM=1d-1*DABS(C7HIG)+1d-1*DABS(C7CHAR)+3d-1*DABS(C70S0)
     .      +asc0/4d0/Pi*(3d-1*DABS(C71HIG)+3d-1*DABS(C71CHAR))
       DC8BSM=1d-1*DABS(C8HIG)+1d-1*DABS(C8CHAR)+3d-1*DABS(C80S0)
     .      +asc0/4d0/Pi*(3d-1*DABS(C81HIG)+3d-1*DABS(C81CHAR))

*  Finally, the total Wilson coefficients at sc0
       C70=C70SM+C70BSM
       C80=C80SM+C80BSM
       C11=C11SM
       C41=C41SM+C41BSM
       C71=C71SM+C71BSM
       C81=C81SM+C81BSM

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
*      Summation of the Electroweak Corrections, the SM Contribution
*      is set to 0.0071.
      EPSew=0.0071d0
     .     +ALEMMZ*(C7EMb/asmb-(eta**(aa(2))*C70BSM
     .   +8d0/3d0*(eta**(aa(1))-eta**(aa(2)))*C80BSM)*dlog(MZ/scb)/pi)

*   BSM error estimate from DC7BSM and DC8BSM:
      DC7BSM_b=eta**(aa(2))*DABS(DC7BSM)
     .   +8d0/3d0*DABS(eta**(aa(1))-eta**(aa(2)))*DABS(DC8BSM)
      DC8BSM_b=eta**(aa(1))*DABS(DC8BSM)

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
     .                 *(C10b+asmb/4d0/Pi*C11b)**2
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
     .      -2d0*VVtu*0.0136471d0
     .                 *(C10b+asmb/4d0/Pi*C11b)*(C80b+asmb/4d0/Pi*C81b)
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
     .      +2d0*0.0818824d0*VVtu
     .                 *(C20b+asmb/4d0/Pi*C21b)*(C80b+asmb/4d0/Pi*C81b)
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

      BSGPERT=C70b**2                              ! NNLO expansion [12]
     . +asmb/2d0/Pi*(C70b*C71b
     .              +ff11*C10b**2+ff12*C10b*C20b+ff22*C20b**2
     .              +K17*C10b*C70b+K27*C20b*C70b+K37*C30b*C70b
     .              +K47*C40b*C70b+K57*C50b*C70b+K67*C60b*C70b
     .              +K77*C70b**2+K78*C70b*C80b
     .              +ff18*C10b*C80b+ff28*C20b*C80b+ff48*C40b*C80b
     .              +ff88*C80b**2)
     . +(asmb/4d0/Pi)**2*(C71b**2
     .    +2d0*(2d0*ff11*C10b*C11b+ff12*(C10b*C21b+C11b*C20b)
     .         +2d0*ff22*C20b*C21b
     .         +K17*(C10b*C71b+C11b*C70b)
     .         +K27*(C20b*C71b+C21b*C70b)
     .         +K37*(C30b*C71b+C31b*C70b)
     .         +K47*(C40b*C71b+C41b*C70b)
     .         +K57*(C50b*C71b+C51b*C70b)
     .         +K67*(C60b*C71b+C61b*C70b)
     .         +2d0*K77*(C70b*C71b)
     .         +K78*(C70b*C81b+C71b*C80b)
     .         +ff18*(C10b*C81b+C11b*C80b)
     .         +ff28*(C20b*C81b+C21b*C80b)
     .         +ff48*(C40b*C81b+C41b*C80b)
     .         +2d0*ff88*(C80b*C81b)))
     . +2d0*EPSEW*C70b+0d-4/aux
     . +LO4B
     . +6.256900288d-3 !6.382d-3 6.7d-3                ! Mimicking the missing SM NNLO
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

      HQET=NP17*(C20b-1d0/6d0*C10b)*C70b+NP78*C70b*C80b
     .     +NP77*C70b**2+NP88*C80b**2

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
     .                                     *eta**(aa(2))      ! dC70b/dC70BSM
     . +(asmb/2d0/Pi*C70b+(asmb/4d0/Pi)**2*(2d0*C71b          ! d/dC71b
     . +K17*C10b+K27*C20b+K37*C30b+K47*C40b+K57*C50b+K67*C60b
     . +2d0*K77*C70b+K78*C80b))
     .            *37208d0/4761d0*eta**(aa(2))*(eta-1d0)      ! dC71b/dC70BSM
     . +2d0*C70b                                              ! d/dEPSew
     .         *((88d0/575d0*eta**(16d0/23d0)                 ! dEPSew/dC70BSM 
     .            -40d0/69d0*eta**(-7d0/23d0)
     .            +32d0/75d0*eta**(-9d0/23d0))/asmb
     .            -dlog(MZ/scb)/pi*eta**(aa(2)))*ALEMMZ
     . -1.266116316d-3                                        ! NNLO

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
     .          -eta**(aa(2))*(6698884d0/357075d0*eta-297664d0/14283d0))
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

*         - Error estimate:
*     First: Variations of BRSG from BSM uncertainties in DC7,8BSM:
      BRSGmax=BRSG+aux*(DABS(KC7BSM)*DC7BSM+DABS(KC8BSM)*DC8BSM)
      BRSGmin=BRSG-aux*(DABS(KC7BSM)*DC7BSM+DABS(KC8BSM)*DC8BSM)
*     Second: Add the SM + CKM + NP uncertainties +/-0.23d-4 [11] *2sigma
      BRSGmax=BRSGmax+2d0*0.23d-4
      BRSGmin=Max(BRSGmin-2d0*0.23d-4,0d0)
!      print*,'BRBsg',BRSGMIN,BRSG,BRSGMAX,aux*KC7BSM,aux*KC8BSM,aux
*     Comparison with experimental data:
      prob(32)=0d0

      IF(BRSGmin.GE.BRSGexpMax)
     .     PROB(32)=BRSGmin/BRSGexpMax-1d0
      IF(BRSGmax.LE.BRSGexpmin)
     .     PROB(32)=BRSGmax/BRSGexpmin-1d0

!      csqb=csqb+4.d0*(BRSG-(BRSGexpmin+BRSGexpMax)/2.d0)**2
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

      BSGPERT=C70b**2                              ! NNLO expansion [12]
     . +asmb/2d0/Pi*(C70b*C71b
     .              +ff11*C10b**2+ff12*C10b*C20b+ff22*C20b**2
     .              +K17*C10b*C70b+K27*C20b*C70b+K37*C30b*C70b
     .              +K47*C40b*C70b+K57*C50b*C70b+K67*C60b*C70b
     .              +K77*C70b**2+K78*C70b*C80b
     .              +ff18*C10b*C80b+ff28*C20b*C80b+ff48*C40b*C80b
     .              +ff88*C80b**2)
     . +(asmb/4d0/Pi)**2*(C71b**2
     .    +2d0*(2d0*ff11*C10b*C11b+ff12*(C10b*C21b+C11b*C20b)
     .         +2d0*ff22*C20b*C21b
     .         +K17*(C10b*C71b+C11b*C70b)
     .         +K27*(C20b*C71b+C21b*C70b)
     .         +K37*(C30b*C71b+C31b*C70b)
     .         +K47*(C40b*C71b+C41b*C70b)
     .         +K57*(C50b*C71b+C51b*C70b)
     .         +K67*(C60b*C71b+C61b*C70b)
     .         +2d0*K77*(C70b*C71b)
     .         +K78*(C70b*C81b+C71b*C80b)
     .         +ff18*(C10b*C81b+C11b*C80b)
     .         +ff28*(C20b*C81b+C21b*C80b)
     .         +ff48*(C40b*C81b+C41b*C80b)
     .         +2d0*ff88*(C80b*C81b)))
     . -LO4B+dBPERT
     . +2d0*EPSEW*C70b+0d-4/aux
     . +6.382d-3                              ! Mimicking the missing SM NNLO
     . -1.266116316d-3*C70BSM+1.419489589d-5*C80BSM         ! Mimicking the missing NNLO[SM] BSM
                                              !  coefficients (from [14])
      Vbdg=(0.2077d0)**2                               ! (V_ts.V_tb/V_cb)^2
      aux=6d0*ALEM0/(pi*CCSL)*Vbdg*BRSL                !prefactor
      BRDG=aux*(BSGPERT+HQET)

*         - New-Physics coefficients (linearized)
      KC7BSM=(2d0*C70b+asmb/2d0/Pi*(C71b+K17*C10b             ! d/dC70b
     . +K27*C20b+K37*C30b+K47*C40b+K57*C50b+K67*C60b+2d0*K77*C70b
     . +K78*C80b)+2d0*(asmb/4d0/Pi)**2*(K17*C11b+K27*C21b+K37*C31b
     . +K47*C41b+K57*C51b+K67*C61b+2d0*K77*C71b+K78*C81b)+2d0*EPSEW
     . +NP17*(C20b-1d0/6d0*C10b)+NP78*C80b+2d0*NP77*C70b)
     .                                     *eta**(aa(2))      ! dC70b/dC70BSM
     . +(asmb/2d0/Pi*C70b+(asmb/4d0/Pi)**2*(2d0*C71b          ! d/dC71b
     . +K17*C10b+K27*C20b+K37*C30b+K47*C40b+K57*C50b+K67*C60b
     . +2d0*K77*C70b+K78*C80b))
     .            *37208d0/4761d0*eta**(aa(2))*(eta-1d0)      ! dC71b/dC70BSM
     . +2d0*C70b                                              ! d/dEPSew
     .         *((88d0/575d0*eta**(16d0/23d0)                 ! dEPSew/dC70BSM 
     .            -40d0/69d0*eta**(-7d0/23d0)
     .            +32d0/75d0*eta**(-9d0/23d0))/asmb
     .            -dlog(MZ/scb)/pi*eta**(aa(2)))*ALEMMZ
     . -1.266116316d-3                                        ! NNLO

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
     .          -eta**(aa(2))*(6698884d0/357075d0*eta-297664d0/14283d0))
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

*         - Error estimate:
*     First: Variations of BRDG from BSM uncertainties in DC7,8BSM:
      BRDGmax=BRDG+aux*(DABS(KC7BSM)*DC7BSM+DABS(KC8BSM)*DC8BSM)
      BRDGmin=BRDG-aux*(DABS(KC7BSM)*DC7BSM+DABS(KC8BSM)*DC8BSM)
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
      yt=(MT0/CMASS)**2
      z=(CMASS/MW)**2

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
      do i=1,2
      aux=aux+PCOMP(i,1)**2*PMASS(i)**2
     .          *sgn(PMASS(i)**2-M_Bs**2)/
     . dsqrt((PMASS(i)**2-M_Bs**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo
      CPHP=CPHP*aux

      aux=0d0
      do i=1,3
      aux=aux+SCOMP(i,2)**2*SMASS(i)**2
     .         *sgn(SMASS(i)**2-M_Bs**2)/
     . dsqrt((SMASS(i)**2-M_Bs**2)**2+(SMASS(i)*WIDTH(i))**2)
      enddo
      CSHP=CSHP*aux

*         - Chargino / Squark [19]
      CACHAR=0d0

      DO I=1,2
      DO J=1,2
      DO K=1,2                                            ! Z-penguin
       z=(MCH(J)/MCH(I))**2
       yt=(ST(K)/MCH(I))**2
       CACHAR=CACHAR+CCT(J,K)*CCT(I,K)          ! 1 stop / 2 charginos
     .     *(2d0*MCH(J)/MCH(I)*U(J,1)*U(I,1)
     .              *(fc30(z,yt)+asf(ST(K))/4d0/Pi*fc31(z,yt))
     .      -V(J,1)*V(I,1)*(fc40(z,yt)+asf(ST(K))/4d0/Pi*fc41(z,yt)))

       z=(ST(J)/MCH(I))**2
       yt=(ST(K)/MCH(I))**2                     ! 2 stops / 1 chargino
       CACHAR=CACHAR+CCT(I,J)*CCT(I,K)*RST(K,1)*RST(J,1)
     .              *(fc40(yt,z)+asf(ST(K))/4d0/Pi*fc51(yt,z))
      ENDDO

       z=(MCH(J)/MCH(I))**2
       yt=(MUL/MCH(I))**2
       CACHAR=CACHAR                            ! 1 scharm / 2charginos
     .  +VVc*(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .    *V(I,1)*V(J,1)*(2d0*MCH(J)/MCH(I)*U(J,1)*U(I,1)
     .              *(fc30(z,yt)+asf(MUL)/4d0/Pi*fc31(z,yt))
     .      -V(J,1)*V(I,1)*(fc40(z,yt)+asf(MUL)/4d0/Pi*fc41(z,yt)))
      ENDDO

       z=(MUL/MCH(I))**2                       ! 2 scharms / 1 chargino
       CACHAR=CACHAR
     .  +VVc*(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .    *V(I,1)**2*(fc40(z,z)+asf(MUL)/4d0/Pi*fc51(z,z))
      ENDDO
      CACHAR=CACHAR/8d0

      vdq=H2Q
      aux=0d0
      auxe=0d0
      CSCHAR=0d0
      CPCHAR=0d0
      DO I=1,2
      DO J=1,2
      DO K=1,2                                            ! Box
       z=(MCH(J)/MCH(I))**2
       yt=(ST(K)/MCH(I))**2
       xt=(MNL/MCH(I))**2
       aux=aux+CCT(J,K)*CCT(I,K)/MCH(I)**2*(      ! 1 stop / 1 sneutrino / 2 charginos
     .          V(I,1)*V(J,1)*(fc50(z,yt,xt)
     .                          +asf(ST(K))/4d0/Pi*fc81(z,yt,xt))
     . +2d0*(mmu/vdq)**2/g2q*MCH(J)/MCH(I)*U(I,2)*U(J,2)*fc60(z,yt,xt))
!     .                          +asf(ST(K))/4d0/Pi*fc91(z,yt,xt))
       auxe=auxe+CCT(J,K)*CCT(I,K)/MCH(I)**2*(
     .          V(I,1)*V(J,1)*(fc50(z,yt,xt)
     .                          +asf(ST(K))/4d0/Pi*fc81(z,yt,xt)))
       CSCHAR=CSCHAR+CCT(J,K)*RST(K,1)*U(I,2)/MCH(I)**2
     . *(1d0+asf(ST(K))/4d0/Pi*(1d0+2d0*dlog((ST(K)/MGL)**2)))
     .     *mmu/vdq**2/(1d0+epst3)
     . *(V(J,1)*U(I,2)*(fc50(z,yt,xt)+asf(ST(K))/4d0/Pi*fc121(z,yt,xt))
     .  +MCH(J)/MCH(I)*U(J,2)*V(I,1)
     .                *(fc60(z,yt,xt)+asf(ST(K))/4d0/Pi*fc131(z,yt,xt)))
       CPCHAR=CPCHAR-CCT(J,K)*RST(K,1)*U(I,2)/MCH(I)**2
     . *(1d0+asf(ST(K))/4d0/Pi*(1d0+2d0*dlog((ST(K)/MGL)**2)))
     .     *mmu/vdq**2/(1d0+epst3)
     . *(V(J,1)*U(I,2)*(fc50(z,yt,xt)+asf(ST(K))/4d0/Pi*fc121(z,yt,xt))
     .  -MCH(J)/MCH(I)*U(J,2)*V(I,1)
     .                *(fc60(z,yt,xt)+asf(ST(K))/4d0/Pi*fc131(z,yt,xt)))
      ENDDO
       z=(MCH(J)/MCH(I))**2
       yt=(MUL/MCH(I))**2
       xt=(MNL/MCH(I))**2
       aux=aux                                    ! 1 scharm / 1 sneutrino / 2 charginos
     .  +VVc*(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .    *V(I,1)*V(J,1)/MCH(I)**2*(V(I,1)*V(J,1)*(fc50(z,yt,xt)
     .                          +asf(MUL)/4d0/Pi*fc81(z,yt,xt))
     . +2d0*(mmu/vdq)**2/g2q*MCH(J)/MCH(I)*U(I,2)*U(J,2)*fc60(z,yt,xt))
!     .                          +asf(MUL)/4d0/Pi*fc91(z,yt,xt)))
       auxe=auxe+VVc
     .    *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .               *(V(I,1)*V(J,1))**2/MCH(I)**2
     .             *(fc50(z,yt,xt)+asf(MUL)/4d0/Pi*fc81(z,yt,xt))
       CSCHAR=CSCHAR
     . +VVc*V(J,1)*U(I,2)/MCH(I)**2*mmu/vdq**2/(1d0+epst3)
     .      *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))
     . *(1d0+asf(MUL)/4d0/Pi*(1d0+2d0*dlog((MUL/MGL)**2)))
     . *(V(J,1)*U(I,2)*(fc50(z,yt,xt)+asf(MUL)/4d0/Pi*fc121(z,yt,xt))
     .  +MCH(J)/MCH(I)*U(J,2)*V(I,1)
     .                *(fc60(z,yt,xt)+asf(MUL)/4d0/Pi*fc131(z,yt,xt)))
       CPCHAR=CPCHAR
     . -VVc*(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))
     .         *V(J,1)*U(I,2)/MCH(I)**2*mmu/vdq**2/(1d0+epst3)
     . *(1d0+asf(MUL)/4d0/Pi*(1d0+2d0*dlog((MUL/MGL)**2)))
     . *(V(J,1)*U(I,2)*(fc50(z,yt,xt)+asf(MUL)/4d0/Pi*fc121(z,yt,xt))
     .  -MCH(J)/MCH(I)*U(J,2)*V(I,1)
     .                *(fc60(z,yt,xt)+asf(MUL)/4d0/Pi*fc131(z,yt,xt)))
      ENDDO
      ENDDO

      DO I=1,2
      DO J=1,2
      DO K=1,2
      DO L=1,2                                            ! Quartic squark coupling
      DO M=1,2
       z=(MCH(J)/MCH(I))**2                       ! stops
       yt=(ST(K)/MCH(I))**2
       zt=(ST(L)/MCH(I))**2
       xt=(MNL/MCH(I))**2
       aux=aux-asf(ST(M))/3d0/Pi*CCT(J,K)*CCT(I,L)*ST(M)**2/MCH(I)**4
     .           *(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .           *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .  *(V(I,1)*V(J,1)*fc90(z,yt,zt,xt)+2d0*(mmu/vdq)**2/g2q
     .          *MCH(J)/MCH(I)*U(I,2)*U(J,2)*fc100(z,yt,zt,xt))
     . -asf(ST(M))/6d0/Pi/MW**2*CCT(J,K)*CCT(I,L)*ST(M)**2/MCH(I)**2
     .           *(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .           *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .  *(2d0*MCH(J)/MCH(I)*U(I,1)*U(J,1)*fc60(z,yt,zt)
     .          -V(I,1)*V(J,1)*fc50(z,yt,zt))
       auxe=auxe-asf(ST(M))/3d0/Pi*CCT(J,K)*CCT(I,L)*ST(M)**2/MCH(I)**4
     .           *(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .           *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .               *V(I,1)*V(J,1)*fc90(z,yt,zt,xt)
     . -asf(ST(M))/6d0/Pi/MW**2*CCT(J,K)*CCT(I,L)*ST(M)**2/MCH(I)**2
     .           *(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .           *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .  *(2d0*MCH(J)/MCH(I)*U(I,1)*U(J,1)*fc60(z,yt,zt)
     .          -V(I,1)*V(J,1)*fc50(z,yt,zt))
       CSCHAR=CSCHAR-asf(ST(M))/3d0/Pi*CCT(J,K)*RST(L,1)*U(I,2)
     . *(1d0+asf(ST(M))/4d0/Pi*(1d0+2d0*dlog((ST(M)/MGL)**2)))
     .       *mmu/vdq**2/(1d0+epst3)
     . *ST(M)**2/MCH(I)**4*(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .           *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .  *(U(I,2)*V(J,1)*fc90(z,yt,zt,xt)
     .       +MCH(J)/MCH(I)*V(I,1)*U(J,2)*fc100(z,yt,zt,xt))
       CPCHAR=CPCHAR+asf(ST(M))/3d0/Pi*CCT(J,K)*RST(L,1)*U(I,2)
     . *(1d0+asf(ST(M))/4d0/Pi*(1d0+2d0*dlog((ST(M)/MGL)**2)))
     .         *mmu/vdq**2/(1d0+epst3)
     . *ST(M)**2/MCH(I)**4*(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .           *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .  *(U(I,2)*V(J,1)*fc90(z,yt,zt,xt)
     .       -MCH(J)/MCH(I)*V(I,1)*U(J,2)*fc100(z,yt,zt,xt))
       z=(ST(L)/MCH(I))**2
       yt=(ST(J)/MCH(I))**2
       zt=(ST(K)/MCH(I))**2
       aux=aux-asf(ST(M))/6d0/Pi*CCT(I,K)*CCT(I,L)*ST(M)**2/MCH(I)**2
     . *fc50(z,yt,zt)/MW**2*(RST(K,1)*RST(J,1)
     .      *(RST(J,1)*RST(M,1)-RST(J,2)*RST(M,2))
     .      *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .      +(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .      *(RST(M,1)*RST(J,1)-RST(M,2)*RST(J,2))
     .                       *RST(J,1)*RST(L,1))
       auxe=auxe-asf(ST(M))/6d0/Pi*CCT(I,K)*CCT(I,L)*ST(M)**2
     . /MCH(I)**2*fc50(z,yt,zt)/MW**2*(RST(K,1)*RST(J,1)
     .      *(RST(J,1)*RST(M,1)-RST(J,2)*RST(M,2))
     .      *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .      +(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .      *(RST(M,1)*RST(J,1)-RST(M,2)*RST(J,2))
     .                       *RST(J,1)*RST(L,1))
      ENDDO
      ENDDO
      ENDDO
       z=(MCH(J)/MCH(I))**2                       ! scharms
       yt=(MUL/MCH(I))**2
       xt=(MNL/MCH(I))**2
       aux=aux-asf(MUL)/3d0/Pi*VVc*MUL**2/MCH(I)**4*V(I,1)*V(J,1)
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .  *(V(I,1)*V(J,1)*fc90(z,yt,yt,xt)+2d0*(mmu/vdq)**2/g2q
     .         *MCH(J)/MCH(I)*U(I,2)*U(J,2)*fc100(z,yt,yt,xt))
     . -asf(MUL)/MW**2/6d0/Pi*VVc*MUL**2/MCH(I)**2*V(I,1)*V(J,1)
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .  *(2d0*MCH(J)/MCH(I)*U(I,1)*U(J,1)*fc60(z,yt,yt)
     .          -V(I,1)*V(J,1)*fc50(z,yt,yt))
     . -asf(MUL)/6d0/Pi*VVc*MUL**2/MCH(I)**2/MW**2*V(I,1)**2
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .         *fc50(yt,yt,yt)
       auxe=auxe-asf(MUL)/3d0/Pi*VVc*MUL**2/MCH(I)**4*V(I,1)*V(J,1)
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .                    *V(I,1)*V(J,1)*fc90(z,yt,yt,xt)
     . -asf(MUL)/MW**2/6d0/Pi*VVc*MUL**2/MCH(I)**2*V(I,1)*V(J,1)
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .  *(2d0*MCH(J)/MCH(I)*U(I,1)*U(J,1)*fc60(z,yt,yt)
     .          -V(I,1)*V(J,1)*fc50(z,yt,yt))
     . -asf(MUL)/6d0/Pi*VVc*MUL**2/MCH(I)**2/MW**2*V(I,1)**2
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .         *fc50(yt,yt,yt)
       CSCHAR=CSCHAR-asf(MUL)/3d0/Pi*VVc*U(I,2)*V(J,1)/(1d0+epst3)
     . *(1d0+asf(MUL)/4d0/Pi*(1d0+2d0*dlog((MUL/MGL)**2)))
     . *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))
     . *mmu/vdq**2*MUL**2/MCH(I)**4*(U(I,2)*V(J,1)*fc90(z,yt,yt,xt)
     .       +MCH(J)/MCH(I)*V(I,1)*U(J,2)*fc100(z,yt,yt,xt))
       CPCHAR=CPCHAR+asf(MUL)/3d0/Pi*VVc*U(I,2)*V(J,1)/(1d0+epst3)
     . *(1d0+asf(MUL)/4d0/Pi*(1d0+2d0*dlog((MUL/MGL)**2)))
     . *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))
     . *MUL**2/MCH(I)**4*mmu/vdq**2*(U(I,2)*V(J,1)*fc90(z,yt,yt,xt)
     .       -MCH(J)/MCH(I)*V(I,1)*U(J,2)*fc100(z,yt,yt,xt))
      ENDDO
      ENDDO

      C10eCHAR=CACHAR+auxe*MW**2/4d0
      CACHAR=CACHAR+aux*MW**2/4d0
      CSCHAR=CSCHAR*MW**2/2d0/g2q
      CPCHAR=CPCHAR*MW**2/2d0/g2q

*         - Higgs penguin from integrated SUSY loops at large tanB [1]
      aux=0d0
      do i=1,3
      aux=aux-(SCOMP(i,1)-SCOMP(i,2)*tanb)*SCOMP(i,2)
     . *sgn(SMASS(i)**2-M_Bs**2)/(cosb
     . *dsqrt((SMASS(i)**2-M_Bs**2)**2+(SMASS(i)*WIDTH(i))**2))
      enddo
      CSHP=CSHP+gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs*aux

      aux=0d0
      do i=1,2
      aux=aux-(-PCOMP(i,1)*cosb-PCOMP(i,1)*sinb*tanb)
     .  *PCOMP(i,1)*tanb*sgn(PMASS(i)**2-M_Bs**2)/
     . dsqrt((PMASS(i)**2-M_Bs**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo
      CPHP=CPHP-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs*aux

*         - Summary
      I=1                              ! 0: SM; 1: NMSSM
      CA=CASM
      CS=0d0
      CP=0d0

      IF(I.eq.1)then
       CA=CA+CAH+CACHAR
       CS=CSH+CSCHAR+CSHP
       CP=CPH+CPCHAR+CPHP
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
     .      ((1d0-4d0*mmu**2/M_Bs**2)/(1d0+MS/Mb)**2*CS**2
     .       +(CP/(1d0+ms/mb)+2d0*mmu/M_Bs**2*CA)**2)

*  Tagged BR, cf. 1204.1737
      ADG=((CP/(1d0+ms/mb)+2d0*mmu/M_Bs**2*CA)**2
     .     -(1d0-4d0*mmu**2/M_Bs**2)/(1d0+MS/Mb)**2*CS**2)
     .    /((1d0-4d0*mmu**2/M_Bs**2)/(1d0+MS/Mb)**2*CS**2
     .       +(CP/(1d0+ms/mb)+2d0*mmu/M_Bs**2*CA)**2)
      BRBMUMU=BRBMUMU*(1d0+0.088d0*ADG)/(1d0-0.088d0**2)    ! ys=0.088+/-0.014

*  Error bars:
      BRBMUMUMAX=Max(
     . ((CP+DCP)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CA+DCA))**2,
     . ((CP+DCP)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CA-DCA))**2,
     . ((CP-DCP)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CA+DCA))**2,
     . ((CP-DCP)/(1d0+ms/mb)+2d0*mmu/M_Bs**2*(CA-DCA))**2)
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
      IF(dabs(CS).gt.DCS)THEN
      BRBMUMUMIN=BRBMUMUMIN
     . +(1d0-4d0*mmu**2/M_Bs**2)/(1d0+MS/Mb)**2*(dabs(CS)-DCS)**2
      ELSE
      BRBMUMUMIN=Max(0d0,BRBMUMUMIN-(1d0-4d0*mmu**2/M_Bs**2)
     .  /(1d0+MS/Mb)**2*(dabs(CS)-DCS)**2)
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

!      csqb=csqb+4.d0*(BRBMUMU-(BRBMUMUexpmin+BRBMUMUexpMax)/2.d0)**2
!     c /((BRBMUMUMax-BRBMUMUmin)**2+(BRBMUMUexpMax-BRBMUMUexpmin)**2)

*	 3) Branching ratio BR[Bd --> mu+mu-]

      CSHP=-mmu*tanb**2/4d0/MW**2*(fh30(xt,z)            !  EW Higgs Penguin
     .      +asc0/4d0/Pi*(fh31(xt,z)
     .                       +8d0*xt*fh30p(xt,z)*dlog((sc0/MT0)**2)))
      CPHP=-CSHP

      aux=0d0
      do i=1,2
      aux=aux+PCOMP(i,1)**2*PMASS(i)**2
     .          *sgn(PMASS(i)**2-M_Bd**2)/
     . dsqrt((PMASS(i)**2-M_Bd**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo
      CPHP=CPHP*aux

      aux=0d0
      do i=1,3
      aux=aux+SCOMP(i,2)**2*SMASS(i)**2
     .         *sgn(SMASS(i)**2-M_Bd**2)/
     . dsqrt((SMASS(i)**2-M_Bd**2)**2+(SMASS(i)*WIDTH(i))**2)
      enddo
      CSHP=CSHP*aux

*         - Higgs penguin from integrated SUSY loops at large tanB [1]
      aux=0d0
      do i=1,3
      aux=aux-(SCOMP(i,1)-SCOMP(i,2)*tanb)*SCOMP(i,2)
     . *sgn(SMASS(i)**2-M_Bd**2)/(cosb
     . *dsqrt((SMASS(i)**2-M_Bd**2)**2+(SMASS(i)*WIDTH(i))**2))
      enddo
      CSHP=CSHP+gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbd*aux

      aux=0d0
      do i=1,2
      aux=aux-(-PCOMP(i,1)*cosb-PCOMP(i,1)*sinb*tanb)
     .  *PCOMP(i,1)*tanb*sgn(PMASS(i)**2-M_Bd**2)/
     . dsqrt((PMASS(i)**2-M_Bd**2)**2+(PMASS(i)*WIDTH(3+i))**2)
      enddo
      CPHP=CPHP-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbd*aux

*         - Summary
      I=1                              ! 0: SM; 1: NMSSM
      CA=CASM
      CS=0d0
      CP=0d0

      IF(I.eq.1)then
       CA=CA+CAH+CACHAR
       CS=CSH+CSCHAR+CSHP
       CP=CPH+CPCHAR+CPHP
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

      aux=GF**4*MW**4/32d0/Pi**5*M_Bd**5*tau_Bd*fBd**2    ! prefactor
     .    *VtdVtb2*dsqrt(1d0-4d0*mmu**2/M_Bd**2)

      BRBdMUMU=aux*
     .      ((1d0-4d0*mmu**2/M_Bd**2)/(1d0+Md/Mb)**2*CS**2
     .       +(CP/(1d0+md/mb)+2d0*mmu/M_Bd**2*CA)**2)

*  Error bars:
      BRBdMUMUMAX=Max(
     . ((CP+DCP)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CA+DCA))**2,
     . ((CP+DCP)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CA-DCA))**2,
     . ((CP-DCP)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CA+DCA))**2,
     . ((CP-DCP)/(1d0+md/mb)+2d0*mmu/M_Bd**2*(CA-DCA))**2)
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
      IF(dabs(CS).gt.DCS)THEN
      BRBdMUMUMIN=BRBdMUMUMIN
     . +(1d0-4d0*mmu**2/M_Bd**2)/(1d0+Md/Mb)**2*(dabs(CS)-DCS)**2
      ELSE
      BRBdMUMUMIN=Max(0d0,BRBdMUMUMIN-(1d0-4d0*mmu**2/M_Bd**2)
     .  /(1d0+Md/Mb)**2*(dabs(CS)-DCS)**2)
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

!      csqb=csqb+4.d0*(BRBdMUMU-(BRBdMUMUexpmin+BRBdMUMUexpMax)/2.d0)**2
!     c /((BRBdMUMUMax-BRBdMUMUmin)**2+(BRBdMUMUexpMax-BRBdMUMUexpmin)**2)


*	VII- Delt B = 1 - BR[B --> Xs l+l-]

*	 1) Wilson coefficients

      sc0=120d0
      eta=asc0/asf(sc0)

*         - LO SM-contribution
      xt=(MT0*(asf(sc0)/asf(MT0))**(12d0/23d0)
     .   *(1d0+7462d0/1587d0*(asf(sc0)-asf(MT0))/4d0/Pi)/MW)**2
      C90SM=YSM(xt)/S2TW+WSM(xt)+4d0/9d0*(1d0-dlog((sc0/MW)**2/xt))

      C100SM=-YSM(xt)/S2TW

      C700SM=0d0
      do i=1,8
       C700SM=C700SM+hh(i)*eta**(aa(i))
      enddo
      C700SM=C700SM+eta**(aa(2))*C70SM
     .   +8d0/3d0*(eta**(aa(1))-eta**(aa(2)))*C80SM

      C800SM=0d0
      do i=1,4
       C800SM=C800SM+h8(i)*eta**(bb(i))
      enddo
      C800SM=C800SM+eta**(aa(1))*(C80SM+313063d0/363036d0)

*         - Charged-Higgs [19]
      yt=(MT0*(asf(sc0)/asf(MT0))**(12d0/23d0)
     .   *(1d0+7462d0/1587d0*(asf(sc0)-asf(MT0))/4d0/Pi)/CMASS)**2
      z=(CMASS/MW)**2

      C9H=(4d0-1d0/S2TW)*(MT0/MW/tanb)**2*(fh20(yt)       ! Z-penguin
     .      +asc0/4d0/Pi*(fh21(yt)+8d0*fh20p(yt)*dlog((sc0/MT0)**2)))
     . -1d0/tanb**2*(fhD0(yt)
     .            +asc0/4d0/Pi*(fhD1(yt)+DeD1(yt)*dlog((sc0/MT0)**2)))

      C10H=CAH/S2TW

*         - Chargino / Squark [19]
      C9CHAR=0d0

      DO I=1,2
      DO J=1,2
      DO K=1,2                                            ! Z-penguin
       z=(MCH(J)/MCH(I))**2
       yt=(ST(K)/MCH(I))**2
       C9CHAR=C9CHAR+CCT(J,K)*CCT(I,K)          ! 1 stop / 2 charginos
     .     *(2d0*MCH(J)/MCH(I)*U(J,1)*U(I,1)
     .              *(fc30(z,yt)+asf(ST(K))/4d0/Pi*fc31(z,yt))
     .      -V(J,1)*V(I,1)*(fc40(z,yt)+asf(ST(K))/4d0/Pi*fc41(z,yt)))

       z=(ST(J)/MCH(I))**2
       yt=(ST(K)/MCH(I))**2                     ! 2 stops / 1 chargino
       C9CHAR=C9CHAR+CCT(I,J)*CCT(I,K)*RST(K,1)*RST(J,1)
     .              *(fc40(yt,z)+asf(ST(K))/4d0/Pi*fc51(yt,z))
      ENDDO

       z=(MCH(J)/MCH(I))**2
       yt=(MUL/MCH(I))**2
       C9CHAR=C9CHAR                            ! 1 scharm / 2charginos
     .  +VVc*(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .    *V(I,1)*V(J,1)*(2d0*MCH(J)/MCH(I)*U(J,1)*U(I,1)
     .              *(fc30(z,yt)+asf(MUL)/4d0/Pi*fc31(z,yt))
     .      -V(J,1)*V(I,1)*(fc40(z,yt)+asf(MUL)/4d0/Pi*fc41(z,yt)))
      ENDDO

       z=(MUL/MCH(I))**2                       ! 2 scharms / 1 chargino
       C9CHAR=C9CHAR
     .  +VVc*(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .    *V(I,1)**2*(fc40(z,z)+asf(MUL)/4d0/Pi*fc51(z,z))
      ENDDO

      aux=0d0
      DO I=1,2
      DO J=1,2
      DO K=1,2
      DO L=1,2                                            ! Quartic squark coupling
      DO M=1,2
       z=(MCH(J)/MCH(I))**2                       ! stops
       yt=(ST(K)/MCH(I))**2
       zt=(ST(L)/MCH(I))**2
       xt=(MNL/MCH(I))**2
       C9CHAR=C9CHAR
     . -asf(ST(M))/3d0/Pi*CCT(J,K)*CCT(I,L)*ST(M)**2/MCH(I)**2
     .           *(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .           *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .  *(2d0*MCH(J)/MCH(I)*U(I,1)*U(J,1)*fc60(z,yt,zt)
     .          -V(I,1)*V(J,1)*fc50(z,yt,zt))
       z=(ST(L)/MCH(I))**2
       yt=(ST(J)/MCH(I))**2
       zt=(ST(K)/MCH(I))**2
       C9CHAR=C9CHAR-asf(ST(M))/3d0/Pi*CCT(I,K)*CCT(I,L)
     . *ST(M)**2/MCH(I)**2*fc50(z,yt,zt)*(RST(K,1)*RST(J,1)
     .      *(RST(J,1)*RST(M,1)-RST(J,2)*RST(M,2))
     .      *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .      +(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .      *(RST(M,1)*RST(J,1)-RST(M,2)*RST(J,2))
     .                       *RST(J,1)*RST(L,1))
      ENDDO
      ENDDO
      ENDDO
       z=(MCH(J)/MCH(I))**2                       ! scharms
       yt=(MUL/MCH(I))**2
       xt=(MNL/MCH(I))**2
       C9CHAR=C9CHAR-asf(MUL)/3d0/Pi*VVc*MUL**2/MCH(I)**2
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .  *(V(I,1)*V(J,1)*(2d0*MCH(J)/MCH(I)*U(I,1)*U(J,1)*fc60(z,yt,zt)
     .    -V(I,1)*V(J,1)*fc50(z,yt,zt))
     .    +2d0**V(I,1)**2*fc50(yt,yt,yt))
      ENDDO
      ENDDO
      C9CHAR=C9CHAR*(4d0-1d0/S2TW)/8d0

      auxe=0d0
      aux=0d0
      DO I=1,2
      DO J=1,2
      DO K=1,2                                            ! Box
       z=(MCH(J)/MCH(I))**2
       yt=(ST(K)/MCH(I))**2
       xt=(MNL/MCH(I))**2
       aux=aux+CCT(J,K)*CCT(I,K)/MCH(I)**2*(      ! 1 stop / 1 sneutrino / 2 charginos
     .          V(I,1)*V(J,1)*(fc50(z,yt,xt)
     .                          +asf(ST(K))/4d0/Pi*fc81(z,yt,xt))
     . -2d0*(mmu/vdq)**2/g2q*MCH(J)/MCH(I)*U(I,2)*U(J,2)*fc60(z,yt,xt))
!     .                          +asf(ST(K))/4d0/Pi*fc91(z,yt,xt)))
       auxe=auxe+CCT(J,K)*CCT(I,K)/MCH(I)**2*V(I,1)*V(J,1)
     .          *(fc50(z,yt,xt)+asf(ST(K))/4d0/Pi*fc81(z,yt,xt))
      ENDDO
       z=(MCH(J)/MCH(I))**2
       yt=(MUL/MCH(I))**2
       xt=(MNL/MCH(I))**2
       aux=aux                                    ! 1 scharm / 1 sneutrino / 2 charginos
     .  +VVc*(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .    *V(I,1)*V(J,1)/MCH(I)**2*(V(I,1)*V(J,1)*(fc50(z,yt,xt)
     .                          +asf(MUL)/4d0/Pi*fc81(z,yt,xt))
     . -2d0*(mmu/vdq)**2/g2q*MCH(J)/MCH(I)*U(I,2)*U(J,2)*fc60(z,yt,xt))
!     .                          +asf(MUL)/4d0/Pi*fc91(z,yt,xt)))
       auxe=auxe+VVc*(V(I,1)*V(J,1))**2/MCH(I)**2
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .    *(fc50(z,yt,xt)+asf(MUL)/4d0/Pi*fc81(z,yt,xt))
      ENDDO
      ENDDO

      DO I=1,2
      DO J=1,2
      DO K=1,2
      DO L=1,2                                            ! Quartic squark coupling
      DO M=1,2
       z=(MCH(J)/MCH(I))**2                       ! stops
       yt=(ST(K)/MCH(I))**2
       zt=(ST(L)/MCH(I))**2
       xt=(MNL/MCH(I))**2
       aux=aux-asf(ST(M))/3d0/Pi*CCT(J,K)*CCT(I,L)*ST(M)**2/MCH(I)**4
     .           *(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .           *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .  *(V(I,1)*V(J,1)*fc90(z,yt,zt,xt)-2d0*(mmu/vdq)**2/g2q
     .          *MCH(J)/MCH(I)*U(I,2)*U(J,2)*fc100(z,yt,zt,xt))
       auxe=auxe-asf(ST(M))/3d0/Pi*CCT(J,K)*CCT(I,L)*ST(M)**2/MCH(I)**4
     .           *(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .           *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .                      *V(I,1)*V(J,1)*fc90(z,yt,zt,xt)
      ENDDO
      ENDDO
      ENDDO
       z=(MCH(J)/MCH(I))**2                       ! scharms
       yt=(MUL/MCH(I))**2
       xt=(MNL/MCH(I))**2
       aux=aux-asf(MUL)/3d0/Pi*VVc*MUL**2/MCH(I)**4*V(I,1)*V(J,1)
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .  *(V(I,1)*V(J,1)*fc90(z,yt,yt,xt)-2d0*(mmu/vdq)**2/g2q
     .         *MCH(J)/MCH(I)*U(I,2)*U(J,2)*fc100(z,yt,yt,xt))
       auxe=auxe-asf(MUL)/3d0/Pi*VVc*MUL**2/MCH(I)**4*V(I,1)*V(J,1)
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .           *V(I,1)*V(J,1)*fc90(z,yt,yt,xt)
      ENDDO
      ENDDO
      C9eCHAR=C9CHAR+auxe*MW**2/4d0/S2TW
      C9CHAR=C9CHAR+aux*MW**2/4d0/S2TW

      aux=0d0
      DO I=1,2
      DO K=1,2                                            ! D-term
       z=(ST(K)/MCH(I))**2                         ! stops
       aux=aux+CCT(I,K)**2/MCH(I)**2*(hc30(z)+asf(ST(K))/4d0/Pi*hc31(z))
      ENDDO
       z=(MUL/MCH(I))**2                           ! scharms
       aux=aux+VVc/MCH(I)**2*(hc30(z)+asf(MUL)/4d0/Pi*hc31(z))
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
      ENDDO

      DO I=1,2
      DO K=1,2
      DO L=1,2                                            ! Quartic squark coupling
      DO M=1,2
       yt=(ST(K)/MCH(I))**2                       ! stops
       zt=(ST(L)/MCH(I))**2
       aux=aux+asf(ST(M))/4d0/Pi*CCT(I,K)*CCT(I,L)*ST(M)**2/MCH(I)**4
     .           *(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .           *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .                   *qc51(yt,zt)
      ENDDO
      ENDDO
      ENDDO
       yt=(MUL/MCH(I))**2                        ! scharms
       aux=aux+asf(MUL)/4d0/Pi*VVc*MUL**2/MCH(I)**4*V(I,1)**2
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .                    *qc51(yt,yt)
      ENDDO
      C9CHAR=C9CHAR-aux*MW**2
      C9eCHAR=C9eCHAR-aux*MW**2

      C10CHAR=CACHAR/S2TW
      C10eCHAR=C10eCHAR/S2TW

      C900=C9H+C9CHAR
      C1000=C10H+C10CHAR
      C9e00=C9H+C9eCHAR
      C10e00=C10H+C10eCHAR

*  - C_7,8 from B-->Xs gamma:      (neglecting the imaginary parts, as for C_9,10)
      C700=eta**(aa(2))*C70BSM
     .   +8d0/3d0*(eta**(aa(1))-eta**(aa(2)))*C80BSM

      C800=eta**(aa(1))*C80BSM

      C710=0d0
      C810=0d0
      DO I=1,8
       C710=C710+eta**(aa(i)+1d0)*m0074(i)*C41BSM
       C810=C810+eta**(aa(i)+1d0)*m0084(i)*C41BSM
      ENDDO
      C710=C710+eta**(aa(2)+1d0)*C71BSM
     .         +8d0/3d0*eta*(eta**(aa(1))-eta**(aa(2)))*C81BSM
     .         +37208d0/4761d0*eta**(aa(2))*(eta-1d0)*C70BSM
     .         +(eta**(aa(1))*(256868d0/14283d0*eta-7164416d0/357075d0)
     .          -eta**(aa(2))*(6698884d0/357075d0*eta-297664d0/14283d0))
     .                                               *C80BSM

      C810=C810+eta**(aa(1))*(eta*C81BSM+6.7441d0*(eta-1d0)*C80BSM)

*	 2) Integrated rate between 1 and 6 GeV^2 for mu+mu- [21]

*  - Scalar coefficients (Higgs penguins) [1,19]
      mb1S=4.68d0
      CQ12=(CSH+CSCHAR)**2*IntpropS1(mmu,mb1S)
      CQ22=(CPH+CPCHAR)**2*IntpropP1(mmu,mb1S)

      aux=0d0
      do i=1,3
      aux=aux-(SCOMP(i,1)-SCOMP(i,2)*tanb)*SCOMP(i,2)/cosb
     .              *IntpropS2(mmu,mb1S,SMASS(i))
      enddo
      aux=gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs*aux
      CQ12=CQ12+2d0*(CSH+CSCHAR)*aux

      aux=0d0
      do i=1,2
      aux=aux-(-PCOMP(i,1)*cosb-PCOMP(i,1)*sinb*tanb)
     .  *PCOMP(i,1)*tanb*IntpropP2(mmu,mb1S,PMASS(i))
      enddo
      CQ22=CQ22-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs*aux
     .                                      *2d0*(CPH+CPCHAR)

      CSHP=-mmu*tanb**2/4d0/MW**2*(fh30(xt,z)            !  EW Higgs Penguin
     .      +asc0/4d0/Pi*(fh31(xt,z)
     .                       +8d0*xt*fh30p(xt,z)*dlog((sc0/MT0)**2)))
      CPHP=-CSHP

      aux=0d0
      do i=1,3
      aux=aux+SCOMP(i,2)**2*SMASS(i)**2
     .         *IntpropS2(mmu,mb1S,SMASS(i))
      enddo
      CQ12=CQ12+2d0*(CSH+CSCHAR)*aux*CSHP

      aux=0d0
      do i=1,2
      aux=aux+PCOMP(i,1)**2*PMASS(i)**2
     .          *IntpropP2(mmu,mb1S,PMASS(i))
      enddo
      CQ22=CQ22+2d0*(CPH+CPCHAR)*aux*CPHP

      aux=0d0
      do i=1,3
      do j=1,3
      aux=aux+(-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs
     .         *(SCOMP(i,1)-SCOMP(i,2)*tanb)*SCOMP(i,2)/cosb
     .         +SCOMP(i,2)**2*SMASS(i)**2*CSHP)
     .  *(-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs
     .         *(SCOMP(j,1)-SCOMP(j,2)*tanb)*SCOMP(j,2)/cosb
     .         +SCOMP(j,2)**2*SMASS(j)**2*CSHP)
     .       *IntpropS3(mmu,mb1S,SMASS(i),SMASS(j),WIDTH(i)*WIDTH(j))
      enddo
      enddo
      CQ12=CQ12+aux

      aux=0d0
      do i=1,2
      do j=1,2
      aux=aux
     .      +(-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs*
     .        (-PCOMP(i,1)*cosb-PCOMP(i,1)*sinb*tanb)*PCOMP(i,1)*tanb
     .          +PCOMP(i,1)**2*PMASS(i)**2*CPHP)
     .      *(-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs*
     .        (-PCOMP(j,1)*cosb-PCOMP(j,1)*sinb*tanb)*PCOMP(j,1)*tanb
     .          +PCOMP(j,1)**2*PMASS(j)**2*CPHP)
     .  *IntpropP3(mmu,mb1S,PMASS(i),PMASS(j),WIDTH(3+i)*WIDTH(3+j))
      enddo
      enddo
      CQ22=CQ22+aux

      CQ12=CQ12*(mb1S/S2TW)**2
      CQ22=CQ22*(mb1S/S2TW)**2

*  - SM + 'vector' BSM : Eq. (B.35)
      RR7=(C700+asf(sc0)/4d0/Pi*C710)/C700SM+1d0
      RR8=(C800+asf(sc0)/4d0/Pi*C810)/C800SM+1d0
      RR9=C900/C90SM+1d0
      RR10=C1000/C100SM+1d0

      BRBSmm=-0.0604983d-7*RR7-0.0157368d-7*RR8
     . +2.86075d-7*RR9-0.55849d-7*RR10+0.0685402d-7*RR7*RR8
     . -0.853217d-7*RR7*RR9+0.0153118d-7*RR7*RR10-0.09808d-7*RR8*RR9
     . +0.00185203d-7*RR8*RR10-0.107301d-7*RR9*RR10+0.287891d-7*RR7**2
     . +0.003817d-7*RR8**2+1.48796d-7*RR9**2+10.8441d-7*RR10**2
     . +2.31673d-7

*  - Prefactor Eq. (4.6), normalized by its value in [21]
      scb=5d0
      asmb=asf(scb)
      eta=asf(sc0)/asmb
      ALEMMB=asmb*23d0/8d0/(460d0*Pi*(1d0-eta)+327d0*asmb*dlog(eta))
     . *(-69d0*Pi+dsqrt(3d0*Pi*(1587d0*Pi+7360d0*Pi*ALEMMZ/asmb
     .   *(1d0-eta)+5232d0*ALEMMZ*dlog(eta))))
      aux=BRSL*4d0/CCSL*Vbsg/(1d0-2.41307d0*asmb/Pi
     .  -(21.29553d0+4.625050d0*dlog((scb/mbp)**2))*(asmb/Pi)**2
     .  +12d0/23d0*ALEMMB/asmb*(1d0-1d0/eta)
     .  +(lambd1-9d0*lambd2)/2d0/mbp**2)
      Prefac=aux/(0.1051d0*4d0/0.574d0*0.9621d0
     .  /(1d0-2.41307d0*0.217d0/Pi
     .  -(21.29553d0+4.625050d0*dlog((scb/4.68d0)**2))*(0.217d0/Pi)**2
     .  +12d0/23d0*ALEMMB/0.217d0*(1d0-1d0/eta)
     .  +(-0.362d0-9d0*0.12d0)/2d0/4.68d0**2))
      BRBSmm=BRBSmm*Prefac

*  - 'Scalar' BSM [22]
      BRBSmm=BRBSmm+BRSL*4d0/CCSL*Vbsg*(ALEMMB/(4d0*pi))**2
     .               *3d0/2d0*(cQ12+cQ22)

*  - Uncertainty
      BRBSmmmax=BRBSmm+0.18d-6                               ! SM ; Eq. (5.14) [21]
     . +0.254081d-8*dabs(RR7-1d0)+0.0357906d-8*dabs(RR8-1d0) ! 'vector' BSM, 10%
     . +4.77807d-8*dabs(RR9-1d0)+21.0396d-8*dabs(RR10-1d0)   ! linearizing Eq.(B35)
     . +BRSL*4d0/CCSL*Vbsg*(ALEMMB/(4d0*pi))**2*0.3d0*       ! 'scalar' BSM 30%
     .               *3d0/2d0*(cQ12+cQ22)
      BRBSmmmin=BRBSmm-0.18d-6                               ! SM ; Eq. (5.14)
     . -0.254081d-8*dabs(RR7-1d0)-0.0357906d-8*dabs(RR8-1d0) ! 'vector' BSM, 10%
     . -4.77807d-8*dabs(RR9-1d0)-21.0396d-8*dabs(RR10-1d0)   ! linearizing Eq.(B35)
     . -BRSL*4d0/CCSL*Vbsg*(ALEMMB/(4d0*pi))**2*0.3d0*       ! 'scalar' BSM 30%
     .               *3d0/2d0*(cQ12+cQ22)

*	 3) Integrated rate between 1 and 6 GeV^2 for e+e- [21]

*  - SM + 'vector' BSM : Eq. (B.32)
      RR7=(C700+asf(sc0)/4d0/Pi*C710)/C700SM+1d0
      RR8=(C800+asf(sc0)/4d0/Pi*C810)/C800SM+1d0
      RR9=C9e00/C90SM+1d0
      RR10=C10e00/C100SM+1d0

      BRBSee=-0.114406d-7*RR7-0.0198921d-7*RR8
     . +3.00071d-7*RR9-0.55849d-7*RR10+0.067386d-7*RR7*RR8
     . -0.868076d-7*RR7*RR9+0.0153118d-7*RR7*RR10-0.0992099d-7*RR8*RR9
     . +0.00185203d-7*RR8*RR10-0.107301d-7*RR9*RR10+0.280302d-7*RR7**2
     . +0.00377311d-7*RR8**2+1.52738d-7*RR9**2+11.1181d-7*RR10**2
     . +2.44915d-7

      BRBSee=BRBSee*Prefac

*  - Uncertainty
      BRBSeemax=BRBSee+0.20d-6                               ! SM ; Eq. (5.13) [21]
     . +0.33918d-8*dabs(RR7-1d0)+0.0423178d-8*dabs(RR8-1d0)  ! 'vector' BSM, 10%
     . +4.98088d-8*dabs(RR9-1d0)+21.5876d-8*dabs(RR10-1d0)   ! linearizing Eq.(B32)
      BRBSeemin=BRBSee-0.20d-6                               ! SM ; Eq. (5.13)
     . -0.33918d-8*dabs(RR7-1d0)-0.0423178d-8*dabs(RR8-1d0)  ! 'vector' BSM, 10%
     . -4.98088d-8*dabs(RR9-1d0)-21.5876d-8*dabs(RR10-1d0)   ! linearizing Eq.(B32)

*	 4) Integrated rate beyond 14.4 GeV^2 for mu+mu- [21]

*  - Scalar coefficients (Higgs penguins) [1,19]
      CQ12=(CSH+CSCHAR)**2*IntpropSh1(mmu,mb1S)
      CQ22=(CPH+CPCHAR)**2*IntpropPh1(mmu,mb1S)

      aux=0d0
      do i=1,3
      aux=aux-(SCOMP(i,1)-SCOMP(i,2)*tanb)*SCOMP(i,2)/cosb
     .              *IntpropSh2(mmu,mb1S,SMASS(i))
      enddo
      aux=gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs*aux
      CQ12=CQ12+2d0*(CSH+CSCHAR)*aux

      aux=0d0
      do i=1,2
      aux=aux-(-PCOMP(i,1)*cosb-PCOMP(i,1)*sinb*tanb)
     .  *PCOMP(i,1)*tanb*IntpropPh2(mmu,mb1S,PMASS(i))
      enddo
      CQ22=CQ22-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs*aux
     .                                      *2d0*(CPH+CPCHAR)

      aux=0d0
      do i=1,3
      aux=aux+SCOMP(i,2)**2*SMASS(i)**2
     .         *IntpropSh2(mmu,mb1S,SMASS(i))
      enddo
      CQ12=CQ12+2d0*(CSH+CSCHAR)*aux*CSHP

      aux=0d0
      do i=1,2
      aux=aux+PCOMP(i,1)**2*PMASS(i)**2
     .          *IntpropPh2(mmu,mb1S,SMASS(i))
      enddo
      CQ22=CQ22+2d0*(CPH+CPCHAR)*aux*CPHP

      aux=0d0
      do i=1,3
      do j=1,3
      aux=aux+(-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs
     .         *(SCOMP(i,1)-SCOMP(i,2)*tanb)*SCOMP(i,2)/cosb
     .         +SCOMP(i,2)**2*SMASS(i)**2*CSHP)
     .  *(-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs
     .         *(SCOMP(j,1)-SCOMP(j,2)*tanb)*SCOMP(j,2)/cosb
     .         +SCOMP(j,2)**2*SMASS(j)**2*CSHP)
     .       *IntpropSh3(mmu,mb1S,SMASS(i),SMASS(j),WIDTH(i)*WIDTH(j))
      enddo
      enddo
      CQ12=CQ12+aux

      aux=0d0
      do i=1,2
      do j=1,2
      aux=aux
     .      +(-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs*
     .        (-PCOMP(i,1)*cosb-PCOMP(i,1)*sinb*tanb)*PCOMP(i,1)*tanb
     .          +PCOMP(i,1)**2*PMASS(i)**2*CPHP)
     .      *(-gg2/(2d0*MW)*(pi/(GF*MW))**2*mmu/MB0*sigRLbs*
     .        (-PCOMP(j,1)*cosb-PCOMP(j,1)*sinb*tanb)*PCOMP(j,1)*tanb
     .          +PCOMP(j,1)**2*PMASS(j)**2*CPHP)
     .  *IntpropPh3(mmu,mb1S,PMASS(i),PMASS(j),WIDTH(3+i)*WIDTH(3+j))
      enddo
      enddo
      CQ22=CQ22+aux

      CQ12=CQ12*(mb1S/S2TW)**2
      CQ22=CQ22*(mb1S/S2TW)**2

*  - SM + 'vector' BSM : Eq. (B.38)
      RR7=(C700+asf(sc0)/4d0/Pi*C710)/C700SM+1d0
      RR8=(C800+asf(sc0)/4d0/Pi*C810)/C800SM+1d0
      RR9=C900/C90SM+1d0
      RR10=C1000/C100SM+1d0

      BRBShmm=-0.0871863d-7*RR7-0.00943852d-7*RR8
     . +0.594393d-7*RR9-0.0806142d-7*RR10+0.000835527d-7*RR7*RR8
     . -0.0601984d-7*RR7*RR9+0.00111614d-7*RR7*RR10
     . -0.00716282d-7*RR8*RR9+0.000119004d-7*RR8*RR10
     . -0.0168936 d-7*RR9*RR10+0.00370104d-7*RR7**2
     . +0.0000421485d-7*RR8**2+0.234333d-7*RR9**2+1.66583d-7*RR10**2
     . +0.292268d-7

      BRBShmm=BRBShmm*Prefac

*  - 'Scalar' BSM [22]
      BRBShmm=BRBShmm+BRSL*4d0/CCSL*Vbsg*(ALEMMB/(4d0*pi))**2
     .               *3d0/2d0*(cQ12+cQ22)

*  - Uncertainty
      BRBShmmmax=BRBShmm+0.14d-6                             ! SM ; Eq. (5.15) [21]
     . +0.138031d-8*dabs(RR7-1d0)+0.0155625d-8*dabs(RR8-1d0) ! 'vector' BSM, 10%
     . +0.978804d-8*dabs(RR9-1d0)+2.53114d-8*dabs(RR10-1d0)  ! linearizing Eq.(B38)
     . +BRSL*4d0/CCSL*Vbsg*(ALEMMB/(4d0*pi))**2*0.3d0*       ! 'scalar' BSM 30%
     .               *3d0/2d0*(cQ12+cQ22)
      BRBShmmmin=BRBShmm-0.14d-6                             ! SM ; Eq. (5.15)
     . -0.138031d-8*dabs(RR7-1d0)-0.0155625d-8*dabs(RR8-1d0) ! 'vector' BSM, 10%
     . -0.978804d-8*dabs(RR9-1d0)-2.53114d-8*dabs(RR10-1d0)  ! linearizing Eq.(B38)
     . -BRSL*4d0/CCSL*Vbsg*(ALEMMB/(4d0*pi))**2*0.3d0*       ! 'scalar' BSM 30%
     .               *3d0/2d0*(cQ12+cQ22)

*	 5) Integrated rate beyond 14.4 GeV^2 for e+e- [21]

*  - SM + 'vector' BSM : Eq. (B.37)
      RR7=(C700+asf(sc0)/4d0/Pi*C710)/C700SM+1d0
      RR8=(C800+asf(sc0)/4d0/Pi*C810)/C800SM+1d0
      RR9=C9e00/C90SM+1d0
      RR10=C10e00/C100SM+1d0

      BRBShee=-0.0723471d-7*RR7-0.00827793d-7*RR8
     . +0.511715d-7*RR9-0.0806142d-7*RR10+0.000709678d-7*RR7*RR8
     . -0.0516424d-7*RR7*RR9+0.00111614d-7*RR7*RR10
     . -0.00651216d-7*RR8*RR9+0.000119004d-7*RR8*RR10
     . -0.0168936d-7*RR9*RR10+0.00287361d-7*RR7**2
     . +0.0000373632d-7*RR8**2+0.211548d-7*RR9**2+1.50748d-7*RR10**2
     . +0.200589d-7

      BRBShee=BRBShee*Prefac

*  - Uncertainty
      BRBSheemax=BRBShee+0.14d-6                              ! SM ; Eq. (5.15) [21]
     . +0.116416d-8*dabs(RR7-1d0)+0.0423178d-8*dabs(RR8-1d0)  ! 'vector' BSM, 10%
     . +0.859763d-8*dabs(RR9-1d0)+2.91869d-8*dabs(RR10-1d0)   ! linearizing Eq.(B37)
      BRBSheemin=BRBShee-0.14d-6                              ! SM ; Eq. (5.15)
     . -0.116416d-8*dabs(RR7-1d0)-0.0423178d-8*dabs(RR8-1d0)  ! 'vector' BSM, 10%
     . -0.859763d-8*dabs(RR9-1d0)-2.91869d-8*dabs(RR10-1d0)   ! linearizing Eq.(B37)

*	 6) Summary and comparison with experiment

      BRBSll=(BRBSmm+BRBSee)/2d0
      BRBSllmin=(BRBSmmmin+BRBSeemin)/2d0
      BRBSllMax=(BRBSmmMax+BRBSeeMax)/2d0

*  - Comparison with experimental source
      prob(40)=0d0
      IF(BRBSllmin.GE.BRBSllMaxexp)
     .     PROB(40)=BRBSllmin/BRBSllMaxexp-1d0
      IF(BRBSllMax.LE.BRBSllminexp)
     .     PROB(40)=BRBSllmax/BRBSllminexp-1d0

!      csqb=csqb+4.d0*(BRBSll-(BRBSllMaxexp+BRBSllminexp)/2.d0)**2
!     c /((BRBSllMax-BRBSllmin)**2+(BRBSllMaxexp-BRBSllminexp)**2)

      BRBShll=(BRBShmm+BRBShee)/2d0
      BRBShllmin=(BRBShmmmin+BRBSheemin)/2d0
      BRBShllMax=(BRBShmmMax+BRBSheeMax)/2d0

*  - Comparison with experimental source
      IF(BRBShllmin.GE.BRBShllMaxexp)
     .     PROB(40)=PROB(40)+BRBShllmin/BRBShllMaxexp-1d0
      IF(BRBShllMax.LE.BRBShllminexp)
     .     PROB(40)=PROB(40)+BRBShllmax/BRBShllminexp-1d0

!      csqb=csqb+4.d0*(BRBShll-(BRBShllMaxexp+BRBShllminexp)/2.d0)**2
!     c /((BRBShllMax-BRBShllmin)**2+(BRBShllMaxexp-BRBShllminexp)**2)

!      print*,PROB(40)
!      print*,BRBSllmin,BRBSll,BRBSllMax,BRBSllMinexp,BRBSllMaxexp
!      print*,BRBShllmin,BRBShll,BRBShllMax,BRBShllMinexp,BRBShllMaxexp

*	 7) The B -> K* l+l- system (following [27])

      BRBptoKpll=1.11d-7+0.22d-7*C700+0.27d-7*C9e00-0.27d-7*C10e00     ! Eq.(7)
      BRBptoKpllmin=BRBptoKpll-(1.11d-7*0.3d0
     .      +0.22d-8*dabs(C700)+0.27d-8*(dabs(C9e00)+dabs(C10e00)))
      BRBptoKpllMax=BRBptoKpll+(1.11d-7*0.3d0
     .      +0.22d-8*dabs(C700)+0.27d-8*(dabs(C9e00)+dabs(C10e00)))
     
*  - Comparison with experimental source
c      IF(BRBptoKpllmin.GE.BRBptoKpllMaxexp)
c     .     PROB(40)=PROB(40)+BRBptoKpllmin/BRBptoKpllMaxexp
c      IF(BRBptoKpllMax.LE.BRBptoKpllminexp)
c     .     PROB(40)=PROB(40)+1d0-BRBptoKpllmax/BRBptoKpllminexp
      
      FL_BKll=0.77d0+0.25d0*C700+0.05d0*C9e00                          ! Eq.(4)
      FL_BKllmin=FL_BKll
     .          -(0.08d0+0.25d-1*dabs(C700)+0.05d-1*dabs(C9e00))
      FL_BKllMax=FL_BKll
     .          +(0.08d0+0.25d-1*dabs(C700)+0.05d-1*dabs(C9e00))
      
      S4_BKll=0.29d0                                                   ! Eq.(5)
      S4_BKllmin=0.29d0-0.08d0
      S4_BKllMax=0.29d0+0.08d0
      
      S5_BKll=-0.14d0-0.59d0*C700-0.09d0*C9e00                         ! Eq.(6)
      S5_BKllmin=S5_BKll-(0.04d0
     .                   +0.59d0*dabs(C700)+0.09d0*dabs(C9e00)) 
      S5_BKllMax=S5_BKll+(0.04d0
     .                   +0.59d0*dabs(C700)+0.09d0*dabs(C9e00)) 

c      print*,BRBptoKpllmin,BRBptoKpll,BRBptoKpllMax


*	VIII- Delt B = 1 - BR[B --> Xs nu nu]

*	 1) Wilson coefficients [19]

      sc0=MT0

*         - SM-contribution
      xt=(MT0/MW)**2
      CLSM=(fh10(xt)+asf(sc0)/4d0/Pi*fh11(xt))/4d0         ! Z-penguin
     . -(fh20(xt)+asf(sc0)/4d0/Pi*fh61(xt))                ! Box
      CLSM=0.994d0*CLSM                                    ! mt(mt) correction

*         - Charged-Higgs
      yt=(MT0/CMASS)**2
      CLH=-(MT0/MW/tanb)**2*(fh20(yt)+asc0/4d0/Pi*fh21(yt))! Z-penguin
     .                                                /8d0

*         - Chargino / Squark [19]
      CLeCHAR=0d0

      DO I=1,2
      DO J=1,2
      DO K=1,2                                            ! Z-penguin
       z=(MCH(J)/MCH(I))**2
       yt=(ST(K)/MCH(I))**2
       CLeCHAR=CLeCHAR+CCT(J,K)*CCT(I,K)        ! 1 stop / 2 charginos
     .     *(2d0*MCH(J)/MCH(I)*U(J,1)*U(I,1)
     .              *(fc30(z,yt)+asf(ST(K))/4d0/Pi*fc31(z,yt))
     .      -V(J,1)*V(I,1)*(fc40(z,yt)+asf(ST(K))/4d0/Pi*fc41(z,yt)))

       z=(ST(J)/MCH(I))**2
       yt=(ST(K)/MCH(I))**2                     ! 2 stops / 1 chargino
       CLeCHAR=CLeCHAR+CCT(I,J)*CCT(I,K)*RST(K,1)*RST(J,1)
     .              *(fc40(yt,z)+asf(ST(K))/4d0/Pi*fc51(yt,z))
      ENDDO

       z=(MCH(J)/MCH(I))**2
       yt=(MUL/MCH(I))**2
       CLeCHAR=CLeCHAR                          ! 1 scharm / 2charginos
     .  +VVc*(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .    *V(I,1)*V(J,1)*(2d0*MCH(J)/MCH(I)*U(J,1)*U(I,1)
     .              *(fc30(z,yt)+asf(MUL)/4d0/Pi*fc31(z,yt))
     .      -V(J,1)*V(I,1)*(fc40(z,yt)+asf(MUL)/4d0/Pi*fc41(z,yt)))
      ENDDO

       z=(MUL/MCH(I))**2                       ! 2 scharms / 1 chargino
       CLeCHAR=CLeCHAR
     .  +VVc*(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .    *V(I,1)**2*(fc40(z,z)+asf(MUL)/4d0/Pi*fc51(z,z))
      ENDDO
      CLeCHAR=-CLeCHAR/8d0
      CLtauCHAR=CLeCHAR

      aux=0d0
      
      DO I=1,2                                            ! Box
       xt=(MLL/MCH(I))**2
      DO J=1,2
       z=(MCH(J)/MCH(I))**2
      DO K=1,2
       yt=(ST(K)/MCH(I))**2                     ! 2 stops / 2 charginos

       aux=aux+U(J,1)*U(I,1)*CCT(J,K)*CCT(I,K)/MCH(I)**2
     .   *(fc50(z,yt,xt)+asf(ST(K))/4d0/Pi*fc81(z,yt,xt))
       
      ENDDO
       yt=(MUL/MCH(I))**2                       ! 2 scharms / 2 charginos
       
       aux=aux+VVc*U(J,1)*U(I,1)*V(J,1)*V(I,1)/MCH(I)**2
     . *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))
     .   *(fc50(z,yt,xt)+asf(MUL)/4d0/Pi*fc81(z,yt,xt))
      
      ENDDO
      ENDDO

      CLeCHAR=CLeCHAR+MW**2/4d0*aux

      aux=0d0
      
      DO I=1,2                                            ! Box
      DO J=1,2
       z=(MCH(J)/MCH(I))**2
      DO L=1,2
       xt=(SL(L)/MCH(I))**2
      DO K=1,2
       yt=(ST(K)/MCH(I))**2                     ! 2 stops / 2 charginos

       aux=aux+(U(j,1)*RSL(L,1)-MTAU/gg2/H2Q*U(j,2)*RSL(L,2))
     .      *(U(i,1)*RSL(L,1)-MTAU/gg2/H2Q*U(i,2)*RSL(L,2))
     .   *CCT(J,K)*CCT(I,K)/MCH(I)**2
     .   *(fc50(z,yt,xt)+asf(ST(K))/4d0/Pi*fc81(z,yt,xt))
       
      ENDDO
       yt=(MUL/MCH(I))**2                       ! 2 scharms / 2 charginos
       
       aux=aux+VVc*V(J,1)*V(I,1)/MCH(I)**2
     .    *(U(j,1)*RSL(L,1)-MTAU/gg2/H2Q*U(j,2)*RSL(L,2))
     .    *(U(i,1)*RSL(L,1)-MTAU/gg2/H2Q*U(i,2)*RSL(L,2))
     . *(1d0-asf(MUL)/(4d0*Pi)*(7d0/3d0+2d0*dlog((MUL/MGL)**2)))
     .   *(fc50(z,yt,xt)+asf(MUL)/4d0/Pi*fc81(z,yt,xt))

      ENDDO
      ENDDO
      ENDDO

      CLtauCHAR=CLtauCHAR+MW**2/4d0*aux

      aux=0d0

      DO I=1,2
      DO J=1,2
      DO K=1,2                                            ! Quartic squark coupling
      DO L=1,2                                            ! Z-penguin
      DO M=1,2
       z=(MCH(J)/MCH(I))**2                       ! stops
       yt=(ST(K)/MCH(I))**2
       zt=(ST(L)/MCH(I))**2
       
       aux=aux+asf(ST(M))/6d0/Pi*CCT(J,K)*CCT(I,L)*ST(M)**2/MCH(I)**2
     .           *(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .           *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .  *(2d0*MCH(J)/MCH(I)*U(I,1)*U(J,1)*fc60(z,yt,zt)
     .    -V(I,1)*V(J,1)*fc50(z,yt,zt))
     
       z=(ST(L)/MCH(I))**2
       yt=(ST(J)/MCH(I))**2
       zt=(ST(K)/MCH(I))**2
       
       aux=aux+asf(ST(M))/6d0/Pi*CCT(I,K)*CCT(I,L)*ST(M)**2/MCH(I)**2
     . *fc50(z,yt,zt)*(RST(K,1)*RST(J,1)
     .      *(RST(J,1)*RST(M,1)-RST(J,2)*RST(M,2))
     .      *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .      +(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .      *(RST(M,1)*RST(J,1)-RST(M,2)*RST(J,2))
     .                       *RST(J,1)*RST(L,1))
      ENDDO
      ENDDO
      ENDDO
       z=(MCH(J)/MCH(I))**2                       ! scharms
       yt=(MUL/MCH(I))**2
       aux=aux
     . +asf(MUL)/6d0/Pi*VVc*MUL**2/MCH(I)**2*V(I,1)*V(J,1)
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .  *(2d0*MCH(J)/MCH(I)*U(I,1)*U(J,1)*fc60(z,yt,yt)
     .          -V(I,1)*V(J,1)*fc50(z,yt,yt))
     . +asf(MUL)/3d0/Pi*VVc*MUL**2/MCH(I)**2*V(I,1)**2
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .         *fc50(yt,yt,yt)
      ENDDO
      ENDDO

      CLeCHAR=CLeCHAR+aux/4d0
      CLtauCHAR=CLtauCHAR+aux/4d0
      
      aux=0d0

      DO I=1,2
      DO J=1,2
      DO K=1,2                                            ! Quartic squark coupling
      DO L=1,2                                            ! Box
      DO M=1,2
       z=(MCH(J)/MCH(I))**2                       ! stops
       yt=(ST(K)/MCH(I))**2
       zt=(ST(L)/MCH(I))**2
       xt=(MLL/MCH(I))**2
       
       aux=aux-asf(ST(M))/3d0/Pi*CCT(J,K)*CCT(I,L)*U(J,1)*U(I,1)
     .           *(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .           *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .  *fc90(z,yt,zt,xt)/ST(M)**2

      ENDDO
      ENDDO
      ENDDO
       z=(MCH(J)/MCH(I))**2                       ! scharms
       yt=(MUL/MCH(I))**2
       xt=(MLL/MCH(I))**2
       aux=aux
     . -asf(MUL)/3d0/Pi*VVc*V(I,1)*V(J,1)*U(J,1)*U(I,1)
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .  *fc90(z,yt,yt,xt)/MUL**2
      ENDDO
      ENDDO
      
      CLeCHAR=CLeCHAR+MW**2/4d0*aux
      
      aux=0d0

      DO I=1,2
      DO J=1,2
       z=(MCH(J)/MCH(I))**2
      DO P=1,2
       xt=(SL(P)/MCH(I))**2
      DO K=1,2                                            ! Quartic squark coupling
      DO L=1,2                                            ! Box
      DO M=1,2
       yt=(ST(K)/MCH(I))**2                       ! stops
       zt=(ST(L)/MCH(I))**2
       
       aux=aux-asf(ST(M))/3d0/Pi*CCT(J,K)*CCT(I,L)
     .    *(U(j,1)*RSL(P,1)-MTAU/gg2/H2Q*U(j,2)*RSL(P,2))
     .    *(U(i,1)*RSL(P,1)-MTAU/gg2/H2Q*U(i,2)*RSL(P,2))
     .           *(RST(K,1)*RST(M,1)-RST(K,2)*RST(M,2))
     .           *(RST(M,1)*RST(L,1)-RST(M,2)*RST(L,2))
     .  *fc90(z,yt,zt,xt)/ST(M)**2

      ENDDO
      ENDDO
      ENDDO
       z=(MCH(J)/MCH(I))**2                       ! scharms
       yt=(MUL/MCH(I))**2
       aux=aux
     . -asf(MUL)/3d0/Pi*VVc*V(I,1)*V(J,1)
     .    *(U(j,1)*RSL(P,1)-MTAU/gg2/H2Q*U(j,2)*RSL(P,2))
     .    *(U(i,1)*RSL(P,1)-MTAU/gg2/H2Q*U(i,2)*RSL(P,2))
     .  *(1d0-asf(MUL)/4d0/Pi*(7d0/2d0+2d0*dlog((MUL/MGL)**2)))**2
     .  *fc90(z,yt,yt,xt)/MUL**2
      ENDDO
      ENDDO
      ENDDO
      
      CLtauCHAR=CLtauCHAR+MW**2/4d0*aux
      
*   Summary
      RRL=2d0/3d0*(CLSM+CLH+CLeCHAR)**2
     .    +1d0/3d0*(CLSM+CLH+CLtauCHAR)**2
      RRL=RRL/CLSM**2

*	 2) The branching ratio BR[B -> Xs nu nu]

*  SM estimate [28]
      BRBXsnunu=2.9d-5
*  Incorporating BSM
      BRBXsnunu=BRBXsnunu*RRL
*  Error estimate
      BRBXsnunumin=BRBXsnunu-2d0*0.3d-5                ! SM / CKM / Form factor uncertainty
     .   -0.2d0*(2d0/3d0*(dabs(CLH)+dabs(CLeCHAR))     ! 10% New Physics
     .          +1d0/3d0*(dabs(CLH)+dabs(CLtauCHAR)))/dabs(CLSM)*2.9d-5
      BRBXsnunuMax=BRBXsnunu+2d0*0.3d-5                ! SM / CKM / Form factor uncertainty
     .   +0.2d0*(2d0/3d0*(dabs(CLH)+dabs(CLeCHAR))     ! 10% New Physics
     .          +1d0/3d0*(dabs(CLH)+dabs(CLtauCHAR)))/dabs(CLSM)*2.9d-5
!      print*,BRBXsnunumin,BRBXsnunu,BRBXsnunuMax
*  Comparison with experiment
      PROB(57)=0d0
      IF(BRBXsnunumin.GE.BRBXsnunuexpMax)
     .     PROB(57)=BRBXsnunumin/BRBXsnunuexpMax-1d0

*	 3) The branching ratio BR[B+ -> K+ nu nu]

*  SM estimate [28]
      BRBpKpnunu=3.98d-6
*  Incorporating BSM
      BRBpKpnunu=BRBpKpnunu*RRL
*  Error estimate
      BRBpKpnunumin=BRBpKpnunu-2d0*0.47d-6             ! SM / CKM / Form factor uncertainty
     .   -0.2d0*(2d0/3d0*(dabs(CLH)+dabs(CLeCHAR))     ! 10% New Physics
     .          +1d0/3d0*(dabs(CLH)+dabs(CLtauCHAR)))/dabs(CLSM)*3.98d-6
      BRBpKpnunuMax=BRBpKpnunu+2d0*0.47d-6             ! SM / CKM / Form factor uncertainty
     .   +0.2d0*(2d0/3d0*(dabs(CLH)+dabs(CLeCHAR))     ! 10% New Physics
     .          +1d0/3d0*(dabs(CLH)+dabs(CLtauCHAR)))/dabs(CLSM)*3.98d-6
!      print*,BRBpKpnunumin,BRBpKpnunu,BRBpKpnunuMax
*  Comparison with experiment
      IF(BRBpKpnunumin.GE.BRBpKpnunuexpMax)
     .     PROB(57)=PROB(57)+BRBpKpnunumin/BRBpKpnunuexpMax-1d0

*	 4) The branching ratio BR[B0 -> K0* nu nu]

*  SM estimate [28]
      BRBKsnunu=9.19d-6
*  Incorporating BSM
      BRBKsnunu=BRBKsnunu*RRL
*  Error estimate
      BRBKsnunumin=BRBKsnunu-2d0*0.99d-6               ! SM / CKM / Form factor uncertainty
     .   -0.2d0*(2d0/3d0*(dabs(CLH)+dabs(CLeCHAR))     ! 10% New Physics
     .          +1d0/3d0*(dabs(CLH)+dabs(CLtauCHAR)))/dabs(CLSM)*9.19d-6
      BRBKsnunuMax=BRBKsnunu+2d0*0.99d-6               ! SM / CKM / Form factor uncertainty
     .   +0.2d0*(2d0/3d0*(dabs(CLH)+dabs(CLeCHAR))     ! 10% New Physics
     .          +1d0/3d0*(dabs(CLH)+dabs(CLtauCHAR)))/dabs(CLSM)*9.19d-6
!      print*,BRBKsnunumin,BRBKsnunu,BRBKsnunuMax
*  Comparison with experiment
      IF(BRBKsnunumin.GE.BRBKsnunuexpMax)
     .     PROB(57)=PROB(57)+BRBKsnunumin/BRBKsnunuexpMax-1d0

      return
      END

*******************************************************************

      double precision function BB0L(x,y,z,t)

      implicit none
      double precision x,y,z,t,BB0
      if(dabs(y-z).lt.1d-5)then
      BB0L=1d0/2d0/x
      if(dabs(y-x).gt.1d-5)BB0L=(-x+y+x*dlog(x/y))/(x-y)**2
      else
      BB0L=(BB0(0d0,x,y,t)-BB0(0d0,x,z,t))/(y-z)
      endif
      BB0L=-BB0L
      return
      end

*********************************************************************

      double precision function echi(x)
      
      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1d-3)then
      d=1d0/(x-1d0)**3
      dg=dlog(x)/(x-1d0)**4
      echi=x*(11d0-7d0*x+2d0*x**2)*d/18d0-x*dg/3d0
      else
      echi=1d0/12d0
      endif
      return
      end

*********************************************************************

      double precision function esfe(x,y)
      
      implicit none
      double precision x,y
      if(dabs(x-1d0).gt.1d-3.and.dabs(y-1d0).gt.1d-3
     .                             .and.dabs(x-y).gt.1d-3)then
      esfe=4d0/9d0/(x-y)*(dlog(x)/(x-1d0)**4-dlog(y)/(y-1d0)**4)
     . +(4d0*x**2*y**2-14d0*x*y**2+22d0*y**2-14d0*x**2*y+52d0*x*y
     .  -62d0*y+22d0*x**2-62d0*x+52d0)/27d0/(x-1d0)**3/(y-1d0)**3
      else
      if(dabs(x-1d0).gt.1d-3)THEN
      esfe=(25d0-48d0*x+36d0*x**2-16d0*x**3+3d0*x**4+12d0*dlog(x))
     . /27d0/(x-1d0)**5
      endif
      if(dabs(y-1d0).gt.1d-3)THEN
      esfe=(25d0-48d0*y+36d0*y**2-16d0*y**3+3d0*y**4+12d0*dlog(y))
     . /27d0/(y-1d0)**5
      endif
      if(dabs(y-x).lt.1d-3.and.dabs(x-1d0).gt.1d-3)THEN
      esfe=4d0*(3d0+10d0*x-18d0*x**2+6d0*x**3-x**4+12d0*x*dlog(x))
     .                       /27d0/x/(1d0-x**5)**5
      endif
      if(dabs(x-1d0).lt.1d-3.and.dabs(y-1d0).lt.1d-3)esfe=4d0/45d0
      endif
      return
      end

*********************************************************************

      double precision function H17(x)
      
      implicit none
      double precision x,Sp2
      if(dabs(x-1d0).gt.1d-3)then
      H17=(24d0*x**3+52d0*x**2-32d0*x)*Sp2(1d0-1d0/x)/9d0/(1d0-x)**4
     . +(-189d0*x**3-783d0*x**2+425d0*x+43d0)*dlog(x)/81d0/(1d0-x)**5
     . +(-1030d0*x**3-1899d0*x**2+1332d0*x+85d0)/243d0/(1d0-x)**4
      else
      H17=-119d0/1620d0
      endif
      return
      end

*********************************************************************

      double precision function H27(x)
      
      implicit none
      double precision x,Sp2
      if(dabs(x-1d0).gt.1d-3)then
      H27=(112d0*x**2-48d0*x)*Sp2(1d0-1d0/x)/9d0/(1d0-x)**3
     . +(12d0*x**3-176d0*x**2+64d0*x+16d0)*dlog(x)/9d0/(1d0-x)**4
     . +(-170d0*x**2+66d0*x+20d0)/9d0/(1d0-x)**3
      else
      H27=-1d0/81d0
      endif
      return
      end

*********************************************************************

      double precision function H18(x)
      
      implicit none
      double precision x,Sp2
      if(dabs(x-1d0).gt.1d-3)then
      H18=(-9d0*x**3-46d0*x**2-49d0*x)*Sp2(1d0-1d0/x)/12d0/(1d0-x)**4
     . +(81d0*x**3+594d0*x**2+1270d0*x+71d0)*dlog(x)/108d0/(1d0-x)**5
     . +(923d0*x**3+3042d0*x**2+6921d0*x+1210d0)/648d0/(1d0-x)**4
      else
      H18=11d0/135d0
      endif
      return
      end

*********************************************************************

      double precision function H28(x)
      
      implicit none
      double precision x,Sp2
      if(dabs(x-1d0).gt.1d-3)then
      H28=(-16d0*x**2-12d0*x)*Sp2(1d0-1d0/x)/3d0/(1d0-x)**3
     . +(52d0*x**2+109d0*x+7d0)*dlog(x)/6d0/(1d0-x)**4
     . +(95d0*x**2+180d0*x+61d0)/12d0/(1d0-x)**3
      else
      H28=29d0/54d0
      endif
      return
      end

*********************************************************************

      double precision function Q11(x,y)
      
      implicit none
      double precision x,y
      if(dabs(x-1d0).gt.1d-1.and.dabs(y-1d0).gt.1d-1
     .                             .and.dabs(x-y).gt.1d-3)then
      Q11=4d0/3d0/(x-y)
     .    *(x**2*dlog(x)/(x-1d0)**4-y**2*dlog(y)/(y-1d0)**4)
     . +(4d0*x**2*y**2+10d0*x*y**2-2d0*y**2+10d0*x**2*y-44d0*x*y
     .  +10d0*y-2d0*x**2+10d0*x+4d0)/9d0/(x-1d0)**3/(y-1d0)**3
      else
      if(dabs(x-1d0).gt.1d-1)THEN
      Q11=(-1d0+8d0*x-8d0*x**3+x**4+12d0*x**2*dlog(x))
     . /9d0/(x-1d0)**5
      endif
      if(dabs(y-1d0).gt.1d-1)THEN
      Q11=(-1d0+8d0*y-8d0*y**3+y**4+12d0*y**2*dlog(y))
     . /9d0/(y-1d0)**5
      endif
      if(dabs(y-x).lt.1d-1.and.dabs(x-1d0).gt.1d-1)THEN
      Q11=4d0*(1d0+9d0*x-9d0*x**2+6d0*x*(1d0+x)*dlog(x))
     .                       /9d0/(1d0-x)**5
      endif
      if(dabs(x-1d0).lt.1d-1.and.dabs(y-1d0).lt.1d-1)Q11=2d0/45d0
      endif
      return
      end

*********************************************************************

      double precision function Q21(x,y)
      
      implicit none
      double precision x,y
      if(dabs(x-1d0).gt.1d-1.and.dabs(y-1d0).gt.1d-1
     .                             .and.dabs(x-y).gt.1d-1)then
      Q21=4d0/3d0/(x-y)*(x*dlog(x)/(x-1d0)**4-y*dlog(y)/(y-1d0)**4)
     . +(-2d0*x**2*y**2+10d0*x*y**2+4d0*y**2+10d0*x**2*y-20d0*x*y
     .  -14d0*y+4d0*x**2-14d0*x+22d0)/9d0/(x-1d0)**3/(y-1d0)**3
      else
      if(dabs(x-1d0).gt.1d-1)THEN
      Q21=(3d0+10d0*x-18d0*x**2+6d0*x**3-x**4+12d0*x*dlog(x))
     . /9d0/(x-1d0)**5
      endif
      if(dabs(y-1d0).gt.1d-1)THEN
      Q21=(3d0+10d0*y-18d0*y**2+6d0*y**3-y**4+12d0*y*dlog(y))
     . /9d0/(y-1d0)**5
      endif
      if(dabs(y-x).lt.1d-1.and.dabs(x-1d0).gt.1d-1)THEN
      Q21=2d0*(17d0-9d0*x-9d0*x**2+x**3+6d0*x*(1d0+3d0*x)*dlog(x))
     .                       /9d0/(1d0-x)**5
      endif
      if(dabs(x-1d0).lt.1d-1.and.dabs(y-1d0).lt.1d-1)Q21=-1d0/15d0
      endif
      return
      end

*********************************************************************

      double precision function Q31(x,y)
      
      implicit none
      double precision x,y
      if(dabs(x-1d0).gt.1d-3.and.dabs(y-1d0).gt.1d-3
     .                             .and.dabs(x-y).gt.1d-3)then
      Q31=8d0/3d0/(x-y)
     .    *(-x**2*dlog(x)/(x-1d0)**3+y**2*dlog(y)/(y-1d0)**3)
     . +(-12d0*x*y+4d0*y+4d0*x+4d0)/3d0/(x-1d0)**2/(y-1d0)**2
      else
      if(dabs(x-1d0).gt.1d-3)THEN
      Q31=4d0*(1d0-6d0*x+3d0*x**2+2d0*x**3-6d0*x**2*dlog(x))
     . /9d0/(x-1d0)**4
      endif
      if(dabs(y-1d0).gt.1d-3)THEN
      Q31=4d0*(1d0-6d0*y+3d0*y**2+2d0*y**3-6d0*y**2*dlog(y))
     . /9d0/(y-1d0)**4
      endif
      if(dabs(y-x).lt.1d-3.and.dabs(x-1d0).gt.1d-3)THEN
      Q31=4d0*(1d0+4d0*x-5d0*x**2+2d0*x*(2d0+x)*dlog(x))
     .                       /3d0/(1d0-x)**4
      endif
      if(dabs(x-1d0).lt.1d-3.and.dabs(y-1d0).lt.1d-3)Q31=2d0/9d0
      endif
      return
      end

*********************************************************************

      double precision function Q41(x,y)
      
      implicit none
      double precision x,y
      if(dabs(x-1d0).gt.1d-3.and.dabs(y-1d0).gt.1d-3
     .                             .and.dabs(x-y).gt.1d-3)then
      Q41=8d0/3d0/(x-y)*(-x*dlog(x)/(x-1d0)**3+y*dlog(y)/(y-1d0)**3)
     . +(-4d0*x*y-4d0*y-4d0*x+12d0)/3d0/(x-1d0)**2/(y-1d0)**2
      else
      if(dabs(x-1d0).gt.1d-3)THEN
      Q41=-4d0*(2d0+3d0*x-6d0*x**2+x**3+6d0*x*dlog(x))
     . /9d0/(x-1d0)**4
      endif
      if(dabs(y-1d0).gt.1d-3)THEN
      Q41=-4d0*(2d0+3d0*y-6d0*y**2+y**3+6d0*y*dlog(y))
     . /9d0/(y-1d0)**4
      endif
      if(dabs(y-x).lt.1d-3.and.dabs(x-1d0).gt.1d-3)THEN
      Q41=4d0*(5d0-4d0*x-x**2+2d0*x*(1d0+2d0*x)*dlog(x))
     .                       /3d0/(1d0-x)**4
      endif
      if(dabs(x-1d0).lt.1d-3.and.dabs(y-1d0).lt.1d-3)Q41=-2d0/9d0
      endif
      return
      end

*********************************************************************

      double precision function gg7(x)
      
      implicit none
      double precision x,d,dg,sp2
      if(dabs(x-1d0).gt.1d-3)then
      d=Sp2(1d0-1d0/x)/(x-1d0)**4
      dg=dlog(x)/(x-1d0)**5
      gg7=x*(-16d0*x**3-122d0*x**2+80d0*x-8d0)*d/9d0
     .  +(-102d0*x**5-588d0*x**4-2262d0*x**3+3244d0*x**2
     .    -1364d0*x+208d0)*dg/81d0
     .  +(6d0*x**2+46d0*x-28d0)*x**2*dlog(x)**2/3d0/(x-1d0)**5
     .  +(1646d0*x**4+12205d0*x**3-10740d0*x**2+2509d0*x-436d0)/486d0
     .      /(x-1d0)**4
      else
      gg7=-3451d0/9720d0
      endif
      return
      end

*********************************************************************

      double precision function gg8(x)
      
      implicit none
      double precision x,d,dg,sp2
      if(dabs(x-1d0).gt.1d-3)then
      d=Sp2(1d0-1d0/x)/(x-1d0)**4
      dg=dlog(x)/(x-1d0)**5
      gg8=x*(-4d0*x**3+40d0*x**2+41d0*x+1d0)*d/6d0
     .  +(-210d0*x**5+1086d0*x**4+4893d0*x**3+2857d0*x**2
     .    -1994d0*x+280d0)*dg/216d0
     .  +(-17d0*x-31d0)*x**2*dlog(x)**2/2d0/(x-1d0)**5
     .  +(737d0*x**4-14102d0*x**3-28209d0*x**2+610d0*x-508d0)/1296d0
     .      /(x-1d0)**4
      else
      gg8=-9821d0/12960d0
      endif
      return
      end

*********************************************************************

      double precision function Delt7(x)
      
      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1d-3)then
      d=1d0/(x-1d0)**4
      dg=dlog(x)/(x-1d0)**5
      Delt7=(82d0*x**4+383d0*x**3+1086d0*x**2-1111d0*x+208d0)*d/81d0
     .  +2d0*x**2*(-3d0*x**2-23d0*x+14d0)*dg/3d0
      else
      Delt7=1369d0/810d0
      endif
      return
      end

*********************************************************************

      double precision function Delt8(x)
      
      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1d-3)then
      d=1d0/(x-1d0)**4
      dg=dlog(x)/(x-1d0)**5
      Delt8=(77d0*x**4-398d0*x**3-1509d0*x**2-902d0*x+140d0)*d/108d0
     .  +x**2*(17d0*x+31d0)*dg/2d0
      else
      Delt8=869d0/1080d0
      endif
      return
      end

*********************************************************************

      double precision function fh10(x)
      
      implicit none
      double precision x
      if(dabs(x-1d0).gt.1d-3)then
      fh10=x/2d0*((x-6d0)/(x-1d0)+(3d0*x+2d0)*dlog(x)/(x-1d0)**2)
      else
      fh10=3d0/4d0
      endif
      return
      end

*********************************************************************

      double precision function fh11(x)
      
      implicit none
      double precision x,Sp2
      if(dabs(x-1d0).gt.1d-3)then
      fh11=4d0*x/3d0*((29d0+7d0*x+4d0*x**2)/(x-1d0)**2
     .               -(23d0+14d0*x+3d0*x**2)*dlog(x)/(x-1d0)**3
     .               -3d0*(4d0+x**2)*Sp2(1d0-1d0/x)/(x-1d0)**2)
      else
      fh11=35d0/9d0
      endif
      return
      end

*********************************************************************

      double precision function fh61(x)
      
      implicit none
      double precision x,Sp2
      if(dabs(x-1d0).gt.1d-3)then
      fh61=2d0*x/3d0*((29d0+3d0*x)/(x-1d0)**2
     .               -(25d0+7d0*x)*dlog(x)/(x-1d0)**3
     .               -12d0*Sp2(1d0-1d0/x)/(x-1d0)**2)
      else
      fh61=11d0/9d0
      endif
      return
      end

*********************************************************************

      double precision function fh20(x)
      
      implicit none
      double precision x
      if(dabs(x-1d0).gt.1d-3)then
      fh20=-x/(x-1d0)+x*dlog(x)/(x-1d0)**2
      else
      fh20=-1d0/2d0
      endif
      return
      end

*********************************************************************

      double precision function fh20p(x)
      
      implicit none
      double precision x
      if(dabs(x-1d0).gt.1d-3)then
      fh20p=-x*(x-3d0)/(x-1d0)**2-2d0*x*dlog(x)/(x-1d0)**3
      else
      fh20p=-2d0/3d0
      endif
      return
      end

*********************************************************************

      double precision function fh21(x)
      
      implicit none
      double precision x,Sp2
      if(dabs(x-1d0).gt.1d-3)then
      fh21=-32d0*x*(x-3d0)/3d0/(x-1d0)**2
     .     +8d0*x*(3d0*x-11d0)*dlog(x)/3d0/(x-1d0)**3
     .     +8d0*x*(x-2d0)*Sp2(1d0-1d0/x)/(x-1d0)**2
      else
      fh21=26d0/9d0
      endif
      return
      end

*********************************************************************

      double precision function fh70(x,z)
      
      implicit none
      double precision x,z,fh20
      if(dabs(x-z).gt.1d-3)then
      if(dabs(x-1d0).gt.1d-3)then
      if(dabs(z-1d0).gt.1d-3)then
       fh70=x*dlog(x)/(x-1d0)/(x-z)+x*dlog(z)/(z-1d0)/(z-x)
      else
       fh70=fh20(x)
      endif
      else
      if(dabs(z-1d0).gt.1d-3)then
       fh70=fh20(z)/z
      else
       fh70=-1d0/2d0
      endif
      endif
      else
      if(dabs(x-1d0).gt.1d-3)then
       fh70=1d0/(x-1d0)-x*dlog(x)/(x-1d0)**2
      else
       fh70=-1d0/2d0
      endif
      endif
      return
      end

*********************************************************************

      double precision function fh70p(x,z)
      
      implicit none
      double precision x,z
      if(dabs(x-z).gt.1d-3)then
      if(dabs(x-1d0).gt.1d-3)then
      if(dabs(z-1d0).gt.1d-3)then
       fh70p=(z-x**2)*dlog(x)/(x-1d0)**2/(x-z)**2
     .      +z*dlog(z)/(x-z)**2/(z-1d0)
     .      +1d0/(x-1d0)/(x-z)
      else
       fh70p=2d0/(x-1d0)**2-(x+1d0)*dlog(x)/(x-1d0)**3
      endif
      else
      if(dabs(z-1d0).gt.1d-3)then
       fh70p=-(1d0+z)/2d0/(z-1d0)**2+z*dlog(z)/(z-1)**3
      else
       fh70p=-1d0/6d0
      endif
      endif
      else
      if(dabs(x-1d0).gt.1d-3)then
       fh70p=-(1d0+x)/2d0/x/(x-1d0)**2+dlog(x)/(x-1d0)**3
      else
       fh70p=-1d0/6d0
      endif
      endif
      return
      end

*********************************************************************

      double precision function fh71(x,z)
      
      implicit none
      double precision x,z,Sp2,Pi
      Pi=4d0*datan(1d0)
      if(dabs(x-z).gt.1d-3)then
      if(dabs(x-1d0).gt.1d-3)then
      if(dabs(z-1d0).gt.1d-3)then
       fh71=32d0*x/3d0/(x-1d0)/(x-z)+4d0*Pi**2/3d0*x/z
     .      -8d0*x*dlog(x)*(x**2-7d0*z+3d0*x*(1d0+z))/3d0
     .                               /(x-z)**2/(x-1d0)**2
     .      -8d0*x*dlog(z)*(3d0*x-7d0*z)/3d0/(x-z)**2/(z-1d0)
     .      -8d0*x*Sp2(1d0-1d0/x)/(z-1d0)
     .      +8d0*x*Sp2(1d0-z/x)/z/(z-1d0)
      else
       fh71=4d0/3d0*x*(Pi**2+2d0*(11d0-3d0*x)/(x-1d0)**2
     .     -2d0*(10d0-5d0*x+3d0*x**2)*dlog(x)/(x-1d0)**3
     .     -6d0*Sp2(1d0-1d0/x))
      endif
      else
      if(dabs(z-1d0).gt.1d-3)then
       fh71=4d0/3d0*(Pi**2/z+2d0*(1d0-5d0*z)/(z-1d0)**2
     . +2d0*(7d0*z-3d0)*dlog(z)/(z-1d0)**3+6d0*Sp2(1d0-z)/z/(z-1d0))
      else
       fh71=4d0/9d0*(3d0*Pi**2-31d0)
      endif
      endif
      else
      if(dabs(x-1d0).gt.1d-3)then
       fh71=4d0/3d0*(Pi**2+2d0*(x-5d0)/(x-1d0)**2
     .        +2d0*x*(7d0-3d0*x)*dlog(x)/(x-1d0)**3
     .        -6d0*x*Sp2(1d0-1d0/x)/(x-1d0))
      else
       fh71=4d0/9d0*(3d0*Pi**2-31d0)
      endif
      endif
      return
      end

*********************************************************************

      double precision function fh30(x,z)  ! (x* def. in [19])
      
      implicit none
      double precision x,z
      if(dabs(x-z).gt.1d-3)then
      if(dabs(x-1d0).gt.1d-3)then
      if(dabs(z-1d0).gt.1d-3)then
       fh30=x*(x*dlog(x)/(x-1d0)/(x-z)+z*dlog(z)/(z-1d0)/(z-x))
      else
       fh30=x*(-1d0/(x-1d0)+x*dlog(x)/(x-1d0)**2)
      endif
      else
      if(dabs(z-1d0).gt.1d-3)then
       fh30=-1d0/(z-1d0)+z*dlog(z)/(z-1d0)**2
      else
       fh30=1d0/2d0
      endif
      endif
      else
      if(dabs(x-1d0).gt.1d-3)then
       fh30=x*(1d0/(x-1d0)-dlog(x)/(x-1d0)**2)
      else
       fh30=1d0/2d0
      endif
      endif
      return
      end

*********************************************************************

      double precision function fc30(x,z)  ! (x* def. in [19])
      
      implicit none
      double precision x,z
      if(dabs(x-z).gt.1d-3)then
      if(dabs(x-1d0).gt.1d-3)then
      if(dabs(z-1d0).gt.1d-3)then
       fc30=x*dlog(x)/(x-1d0)/(x-z)+z*dlog(z)/(z-1d0)/(z-x)
      else
       fc30=-1d0/(x-1d0)+x*dlog(x)/(x-1d0)**2
      endif
      else
      if(dabs(z-1d0).gt.1d-3)then
       fc30=-1d0/(z-1d0)+z*dlog(z)/(z-1d0)**2
      else
       fc30=1d0/2d0
      endif
      endif
      else
      if(dabs(x-1d0).gt.1d-3)then
       fc30=1d0/(x-1d0)-dlog(x)/(x-1d0)**2
      else
       fc30=1d0/2d0
      endif
      endif
      return
      end

*********************************************************************

      double precision function fh30p(x,z)
      
      implicit none
      double precision x,z
      if(dabs(x-z).gt.1d-3)then
      if(dabs(x-1d0).gt.1d-3)then
      if(dabs(z-1d0).gt.1d-3)then
       fh30p=-x*(x-2d0*z+x*z)*dlog(x)/(x-1d0)**2/(x-z)**2
     .   +x/(x-1d0)/(x-z)
     .   +z**2*dlog(z)/(x-z)**2/(z-1d0)
      else
       fh30p=(1d0+x)/(x-1d0)**2-2d0*x*dlog(x)/(x-1d0)**3
      endif
      else
      if(dabs(z-1d0).gt.1d-3)then
       fh30p=(1d0-3d0*z)/2d0/(z-1d0)**2+z**2*dlog(z)/(z-1d0)**3
      else
       fh30p=1d0/3d0
      endif
      endif
      else
      if(dabs(x-1d0).gt.1d-3)then
       fh30p=(x-3d0)/2d0/(x-1d0)**2+dlog(x)/(x-1d0)**3
      else
       fh30p=1d0/3d0
      endif
      endif
      return
      end

*********************************************************************

      double precision function fh31(x,z)                 ! f(1)14 [19]
      
      implicit none
      double precision x,z,Sp2
      if(dabs(x-z).gt.1d-3)then
      if(dabs(x-1d0).gt.1d-3)then
      if(dabs(z-1d0).gt.1d-3)then
       fh31=32d0*x**2/3d0/(x-1d0)/(x-z)
     .   -8d0*x**2*(7d0*x*(1d0+z)-11d0*z-3d0*x**2)*dlog(x)
     .                               /3d0/(x-1d0)**2/(x-z)**2
     .   -8d0*x*z*(3d0*x-7d0*z)*dlog(z)/3d0/(x-z)**2/(z-1d0)
     .   -8d0*x/(z-1d0)*Sp2(1d0-1d0/x)+8d0*x/(z-1d0)*Sp2(1d0-z/x)
      else
       fh31=8d0/3d0*x*((7d0+x)/(x-1d0)**2
     .                -(3d0+5d0*x)*dlog(x)/(x-1d0)**3)
      endif
      else
      if(dabs(z-1d0).gt.1d-3)then
       fh31=8d0/3d0*((-5d0+14d0*z-9d0*z**2)/(z-1d0)**3
     .      +z*(7d0*z-3d0)*dlog(z)/(z-1d0)**3+3d0*Sp2(1d0-z)/(z-1d0))
      else
       fh31=-4d0/9d0
      endif
      endif
      else
      if(dabs(x-1d0).gt.1d-3)then
       fh31=-8d0*x/3d0*((9d0-5d0*x)/(x-1d0)**2
     .    +(3d0*x-7d0)*dlog(x)/(x-1d0)**3+3d0*Sp2(1d0-1d0/x)/(x-1d0))
      else
       fh31=-4d0/9d0
      endif
      endif
      return
      end

*********************************************************************

      double precision function fc40(x,z)
      
      implicit none
      double precision x,z
      if(dabs(x-z).gt.1d-3)then
      if(dabs(x-1d0).gt.1d-3)then
      if(dabs(z-1d0).gt.1d-3)then
       fc40=x**2*dlog(x)/(x-1d0)/(x-z)+z**2*dlog(z)/(z-1d0)/(z-x)
      else
       fc40=-1d0/(x-1d0)+x**2*dlog(x)/(x-1d0)**2
      endif
      else
      if(dabs(z-1d0).gt.1d-3)then
       fc40=-1d0/(z-1d0)+z**2*dlog(z)/(z-1d0)**2
      else
       fc40=3d0/2d0
      endif
      endif
      else
      if(dabs(x-1d0).gt.1d-3)then
       fc40=x/(x-1d0)+x*(x-2d0)*dlog(x)/(x-1d0)**2
      else
       fc40=3d0/2d0
      endif
      endif
      return
      end

*********************************************************************

      double precision function fc41(x,z)
      
      implicit none
      double precision x,z,Sp2
      if(dabs(x-z).gt.1d-3)then
      if(dabs(x-1d0).gt.1d-3)then
      if(dabs(z-1d0).gt.1d-3)then
       fc41=-59d0*x/6d0/(x-z)-z*(59d0-3d0*z)/6d0/(z-1d0)/(x-z)
     .   +4d0*x*(7d0*x**2-3d0*x*z+3d0*z**2)*dlog(x)
     .                        /3d0/(x-1d0)/(x-z)**2+2d0*dlog(z)**2
     .   +4d0*z**2*dlog(z)*(x*(18d0-11d0*z)-z*(11d0-4d0*z))
     .                        /3d0/(x-z)**2/(z-1d0)**2
     .   +4d0*(1d0+z**2)/(z-1d0)/(x-1d0)*Sp2(1d0-1d0/z)
     .   +4d0*(x**2+z**2)/(x-1d0)/(x-z)*Sp2(1d0-x/z)
      else
       fc41=(-47d0+50d0*x-59d0*x**2)/6d0/(x-1d0)**2
     . +4d0*x/3d0*(3d0-3d0*x+7d0*x**2)*dlog(x)/(x-1d0)**3
     . +4d0*(1d0+x**2)*Sp2(1d0-x)/(x-1d0)**2
      endif
      else
      if(dabs(z-1d0).gt.1d-3)then
       fc41=(-3d0+94d0*z+21d0*z**2)/6d0/(z-1d0)**2+2d0*dlog(z)**2
     . +4d0/3d0*(-3d0+3d0*z-21d0*z**2+7d0*z**3)/(z-1d0)**3*dlog(z)
     . -4d0*(-1d0+2d0*z+z**2)/(z-1d0)**2*Sp2(1d0-1d0/z)
      else
       fc41=-13d0/18d0
      endif
      endif
      else
      if(dabs(x-1d0).gt.1d-3)then
       fc41=(-59d0+50d0*x-47d0*x**2)/6d0/(x-1d0)**2+2d0*dlog(x)**2
     . +4d0*x*(18d0-18d0*x+7d0*x**2)*dlog(x)/3d0/(x-1d0)**3
     . +4d0*(1d0+x**2)/(x-1d0)**2*Sp2(1d0-1d0/x)
      else
       fc41=-13d0/18d0
      endif
      endif
      return
      end

*********************************************************************

      double precision function fc31(x,z)
      
      implicit none
      double precision x,z,Sp2
      if(dabs(x-z).gt.1d-3)then
      if(dabs(x-1d0).gt.1d-3)then
      if(dabs(z-1d0).gt.1d-3)then
       fc31=-28d0*z/3d0/(x-z)/(z-1d0)
     .     +2d0*x*(11d0*x+3d0*z)*dlog(x)/3d0/(x-1d0)/(x-z)**2
     .     +2d0*z*(x*(25d0-11d0*z)-z*(11d0+3d0*z))*dlog(z)
     .                        /3d0/(x-z)**2/(z-1d0)**2
     .     +4d0*(1d0+z)/(z-1d0)/(x-1d0)*Sp2(1d0-1d0/z)
     .     +4d0*(x+z)/(x-1d0)/(x-z)*Sp2(1d0-x/z)
      else
       fc31=2d0/3d0/(x-1d0)**3*(2d0*(-4d0+x+3d0*x**2)
     .      +x*(3d0+11d0*x)*dlog(x)+6d0*(x**2-1d0)*Sp2(1d0-x))
      endif
      else
      if(dabs(z-1d0).gt.1d-3)then
       fc31=-2d0/3d0/(z-1d0)**3*(11d0+6d0*z-17d0*z**2
     .  +(6d0+25d0*z-3d0*z**2)*dlog(z)+12d0*(z-1d0)*z*Sp2(1d0-1d0/z))
      else
       fc31=1d0/9d0
      endif
      endif
      else
      if(dabs(x-1d0).gt.1d-3)then
       fc31=2d0/3d0/(x-1d0)**3*(6d0+2d0*x-8d0*x**2
     .   +(11d0+3d0*x)*dlog(x)+6d0*(x**2-1d0)*Sp2(1d0-1d0/x))
      else
       fc31=1d0/9d0
      endif
      endif
      return
      end

*********************************************************************

      double precision function fc5a(x,z)
      
      implicit none
      double precision x,z,Sp2
      if(dabs(x-z).gt.1d-3)then
      if(dabs(x-1d0).gt.1d-3)then
      if(dabs(z-1d0).gt.1d-3)then
       fc5a=4d0*x*(1d0+x*(12d0+z)-z-6d0*x**2)*dlog(x)
     .                                        /3d0/(x-1d0)**2/(x-z)
     . -2d0*(1d0+6d0*x**2*(z-1d0)-3d0*x**3*(z-1d0)+x*(3d0*z-4d0))
     .                   *dlog(x)**2/3d0/(x-1d0)**2/(x-z)/(z-1d0)
     . +4d0*z*(3d0*x**2*(z-1d0)+x*z*(3d0-2d0*z)+z**2*(z-2d0))
     .                 *Sp2(1d0-x/z)/3d0/(x-1d0)/(x-z)**2/(z-1d0)
     . +4d0*(1d0-3d0*x-x**2*(3d0-6d0*z)-x**3)*Sp2(1d0-1d0/x)
     .                              /3d0/(x-1d0)/(x-z)/(z-1d0)
      else
       fc5a=0d0
      endif
      else
      if(dabs(z-1d0).gt.1d-3)then
       fc5a=0d0
      else
       fc5a=0d0
      endif
      endif
      else
      if(dabs(x-1d0).gt.1d-3)then
       fc5a=0d0
      else
       fc5a=0d0
      endif
      endif
      return
      end

*********************************************************************

      double precision function fc51(x,z)
      
      implicit none
      double precision x,z,Sp2,fc5a
      if(dabs(x-z).gt.1d-3)then
      if(dabs(x-1d0).gt.1d-3)then
      if(dabs(z-1d0).gt.1d-3)then
       fc51=-(83d0+27d0*x*(z-1d0)-27d0*z)/6d0/(x-1d0)/(z-1d0)
     . -fc5a(x,z)-fc5a(z,x)+4d0*dlog(x)*(
     . (x-2d0)*x**2*dlog(x)/(x-1d0)**2/(x-z)+(x*z-x-z)/(x-1d0)/(z-1d0)
     . -(z-2d0)*z**2/(x-z)/(z-1d0)**2*dlog(z))
      else
       fc51=-1d0/6d0/(x-1d0)**3*(-23d0+121d0*x-117d0*x**2+19d0*x**3
     . +4d0*dlog(x)*(9d0-22d0*x+47d0*x**2-20d0*x**3)
     . +4d0*dlog(x)**2*(1d0-3d0*x+6d0*x**2-3d0*x**3)
     . +24d0*x**2*Sp2(1d0-x)+8d0*(1d0-3d0*x+5d0*x**3)*Sp2(1d0-1d0/x))
      endif
      else
      if(dabs(z-1d0).gt.1d-3)then
       fc51=-1d0/6d0/(z-1d0)**3*(-23d0+121d0*z-117d0*z**2+19d0*z**3
     . -8d0*z*(2d0-16d0*z+7d0*z**2)*dlog(z)
     . +4d0*(1d0-3d0*z-6d0*z**2+3d0*z**3)*dlog(z)**2
     . +24d0*z**2*Sp2(1d0-z)+8d0*(1d0-3d0*z+5d0*z**3)*Sp2(1d0-1d0/z))
      else
       fc51=-17d0/18d0
      endif
      endif
      else
      if(dabs(x-1d0).gt.1d-3)then
       fc51=1d0/6d0/(x-1d0)**3/x*(7d0*x*(13d0-7d0*x-9d0*x**2+3d0*x**3)
     . +8d0*x**2*(37d0-33d0*x+10d0*x**2)*dlog(x)
     . +12d0*(x-1d0)**3*x*dlog(x)**2
     . -24d0*x*(1d0+x-3d0*x**2+x**3)*Sp2(1d0-1d0/x))
      else
       fc51=-17d0/18d0
      endif
      endif
      return
      end

*********************************************************************

      double precision function fc50(x,y,z)
      
      implicit none
      double precision x,y,z
      if(dabs(x-y).gt.1d-3)then
       if(dabs(x-z).gt.1d-3)then
        if(dabs(y-z).gt.1d-3)then
         if(dabs(x-1d0).gt.1d-3)then
          if(dabs(y-1d0).gt.1d-3)then
           if(dabs(z-1d0).gt.1d-3)then
            fc50=x**2*dlog(x)/(x-1d0)/(x-y)/(x-z)
     .          +y**2*dlog(y)/(y-1d0)/(y-x)/(y-z)
     .          +z**2*dlog(z)/(z-1d0)/(z-y)/(z-x)
           else
            fc50=x**2*dlog(x)/(x-1d0)**2/(x-y)
     .          +y**2*dlog(y)/(y-1d0)**2/(y-x)
     .          +1d0/(x-1d0)/(y-1d0)
           endif
          else
           if(dabs(z-1d0).gt.1d-3)then
            fc50=x**2*dlog(x)/(x-1d0)**2/(x-z)
     .          +z**2*dlog(z)/(z-1d0)**2/(z-x)
     .          +1d0/(x-1d0)/(z-1d0)
           else
            fc50=(-1d0+4d0*x-3d0*x**2+2d0*x**2*dlog(x))
     .           /2d0/(x-1d0)**3
           endif
          endif
         else
          if(dabs(y-1d0).gt.1d-3)then
           if(dabs(z-1d0).gt.1d-3)then
            fc50=z**2*dlog(z)/(z-1d0)**2/(z-y)
     .          +y**2*dlog(y)/(y-1d0)**2/(y-z)
     .          +1d0/(z-1d0)/(y-1d0)
           else
            fc50=(-1d0+4d0*y-3d0*y**2+2d0*y**2*dlog(y))
     .           /2d0/(y-1d0)**3
           endif
          else
           if(dabs(z-1d0).gt.1d-3)then
            fc50=(-1d0+4d0*z-3d0*z**2+2d0*z**2*dlog(z))
     .           /2d0/(z-1d0)**3
           else
            fc50=1d0/3d0
           endif
          endif
         endif
        else
         if(dabs(x-y).gt.1d-3)then
          if(dabs(x-1d0).gt.1d-3)then
           if(dabs(y-1d0).gt.1d-3)then
            fc50=x**2*dlog(x)/(x-1d0)/(x-y)**2
     .         -y*(x*y-2d0*x+y)*dlog(y)/(x-y)**2/(y-1d0)**2
     .         -y/(x-y)/(y-1d0)
           else
            fc50=(-1d0+4d0*x-3d0*x**2+2d0*x**2*dlog(x))
     .          /2d0/(x-1d0)**3
           endif
          else
           if(dabs(y-1d0).gt.1d-3)then
            fc50=(1d0+y)/(y-1d0)**2-2d0*y*dlog(y)/(y-1d0)**3
           else
            fc50=1d0/3d0
           endif
          endif
         else
          if(dabs(x-1d0).gt.1d-3)then
           fc50=(3d0-4d0*x+x**2+2d0*dlog(x))/2d0/(x-1d0)**3
          else
           fc50=1d0/3d0
          endif
         endif
        endif
       else
        if(dabs(x-y).gt.1d-3)then
         if(dabs(x-1d0).gt.1d-3)then
          if(dabs(y-1d0).gt.1d-3)then
           fc50=-x*(x*y+x-2d0*y)*dlog(x)/(x-1d0)**2/(x-y)**2
     .        +y**2*dlog(y)/(x-y)**2/(y-1d0)
     .        +x/(x-y)/(x-1d0)
          else
           fc50=(1d0+x)/(x-1d0)**2-2d0*x*dlog(x)/(x-1d0)**3
          endif
         else
          if(dabs(y-1d0).gt.1d-3)then
           fc50=(-1d0+4d0*y-3d0*y**2+2d0*y**2*dlog(y))
     .        /2d0/(y-1d0)**3
          else
            fc50=1d0/3d0
          endif
         endif
        else
         if(dabs(x-1d0).gt.1d-3)then
          fc50=(3d0-4d0*x+x**2+2d0*dlog(x))/2d0/(x-1d0)**3
         else
          fc50=1d0/3d0
         endif
        endif
       endif
      else
       if(dabs(x-z).gt.1d-3)then
        if(dabs(x-1d0).gt.1d-3)then
         if(dabs(z-1d0).gt.1d-3)then
          fc50=-x*(x*z-2d0*z+x)*dlog(x)/(x-1d0)**2/(x-z)**2
     .        +x/(x-z)/(x-1d0)
     .        +z**2*dlog(z)/(x-z)**2/(z-1d0)
         else
          fc50=(1d0+x)/(x-1d0)**2-2d0*x*dlog(x)/(x-1d0)**3
         endif
        else
         if(dabs(z-1d0).gt.1d-3)then
          fc50=(-1d0+4d0*z-3d0*z**2+2d0*z**2*dlog(z))
     .           /2d0/(z-1d0)**3
         else
          fc50=1d0/3d0
         endif
        endif
       else
        if(dabs(x-1d0).gt.1d-3)then
         fc50=(3d0-4d0*x+x**2+2d0*dlog(x))/2d0/(x-1d0)**3
        else
         fc50=1d0/3d0
        endif
       endif
      endif
      return
      end

*********************************************************************

      double precision function fc60(x,y,z)
      
      implicit none
      double precision x,y,z
      if(dabs(x-y).gt.1d-3)then
       if(dabs(x-z).gt.1d-3)then
        if(dabs(y-z).gt.1d-3)then
         if(dabs(x-1d0).gt.1d-3)then
          if(dabs(y-1d0).gt.1d-3)then
           if(dabs(z-1d0).gt.1d-3)then
            fc60=x*dlog(x)/(x-1d0)/(x-y)/(x-z)
     .          +y*dlog(y)/(y-1d0)/(y-x)/(y-z)
     .          +z*dlog(z)/(z-1d0)/(z-y)/(z-x)
           else
            fc60=x*dlog(x)/(x-1d0)**2/(x-y)+y*dlog(y)/(y-x)/(y-1d0)**2
     .           +1d0/(y-1d0)/(x-1d0)
           endif
          else
           if(dabs(z-1d0).gt.1d-3)then
            fc60=x*dlog(x)/(x-1d0)**2/(x-z)+z*dlog(z)/(z-x)/(z-1d0)**2
     .           +1d0/(z-1d0)/(x-1d0)
           else
            fc60=(1d0-x**2+2d0*x*dlog(x))/2d0/(x-1d0)**3
           endif
          endif
         else
          if(dabs(y-1d0).gt.1d-3)then
           if(dabs(z-1d0).gt.1d-3)then
            fc60=y*dlog(y)/(y-1d0)**2/(y-z)+z*dlog(z)/(z-y)/(z-1d0)**2
     .           +1d0/(z-1d0)/(y-1d0)
           else
            fc60=(1d0-y**2+2d0*y*dlog(y))/2d0/(y-1d0)**3
           endif
          else
           if(dabs(z-1d0).gt.1d-3)then
            fc60=(1d0-z**2+2d0*z*dlog(z))/2d0/(z-1d0)**3
           else
            fc60=-1d0/6d0
           endif
          endif
         endif
        else
         if(dabs(x-y).gt.1d-3)then
          if(dabs(x-1d0).gt.1d-3)then
           if(dabs(y-1d0).gt.1d-3)then
            fc60=x*dlog(x)/(x-1d0)/(x-y)**2
     .         +(x-y**2)*dlog(y)/(y-1d0)**2/(x-y)**2
     .         -1d0/(x-y)/(y-1d0)
           else
            fc60=(1d0-x**2+2d0*x*dlog(x))/2d0/(x-1d0)**3
           endif
          else
           if(dabs(y-1d0).gt.1d-3)then
            fc60=2d0/(y-1d0)**2-(1d0+y)*dlog(y)/(y-1d0)**3
           else
            fc60=-1d0/6d0
           endif
          endif
         else
          if(dabs(x-1d0).gt.1d-3)then
           fc60=(1d0-x**2+2d0*x*dlog(x))/2d0/(x-1d0)**3/x
          else
           fc60=-1d0/6d0
          endif
         endif
        endif
       else
        if(dabs(x-y).gt.1d-3)then
         if(dabs(x-1d0).gt.1d-3)then
          if(dabs(y-1d0).gt.1d-3)then
           fc60=y*dlog(y)/(y-1d0)/(x-y)**2
     .     +(y-x**2)*dlog(x)/(x-1d0)**2/(x-y)**2
     .     +1d0/(x-1d0)/(x-y)
          else
           fc60=2d0/(x-1d0)**2-(1d0+x)*dlog(x)/(x-1d0)**3
          endif
         else
          if(dabs(y-1d0).gt.1d-3)then
           fc60=(1d0-y**2+2d0*y*dlog(y))/2d0/(y-1d0)**3
          else
           fc60=-1d0/6d0
          endif
         endif
        else
         if(dabs(x-1d0).gt.1d-3)then
          fc60=(1d0-x**2+2d0*x*dlog(x))/2d0/(x-1d0)**3/x
         else
          fc60=-1d0/6d0
         endif
        endif
       endif
      else
       if(dabs(x-z).gt.1d-3)then
        if(dabs(x-1d0).gt.1d-3)then
         if(dabs(z-1d0).gt.1d-3)then
          fc60=z*dlog(z)/(z-1d0)/(x-z)**2
     .       +(z-x**2)*dlog(x)/(x-1d0)**2/(x-z)**2
     .       +1d0/(x-1d0)/(x-z)
         else
          fc60=2d0/(x-1d0)**2-(1d0+x)*dlog(x)/(x-1d0)**3
         endif
        else
         if(dabs(z-1d0).gt.1d-3)then
          fc60=(1d0-z**2+2d0*z*dlog(z))/2d0/(z-1d0)**3
         else
          fc60=-1d0/6d0
         endif
        endif
       else
        if(dabs(x-1d0).gt.1d-3)then
         fc60=(1d0-x**2+2d0*x*dlog(x))/2d0/(x-1d0)**3/x
        else
         fc60=-1d0/6d0
        endif
       endif
      endif
      return
      end

*********************************************************************

      double precision function fc81(x,y,z)
      
      implicit none
      double precision x,y,z,Sp2
      if(dabs(x-y).gt.1d-3)then
       if(dabs(x-z).gt.1d-3)then
        if(dabs(y-z).gt.1d-3)then
         if(dabs(x-1d0).gt.1d-3)then
          if(dabs(y-1d0).gt.1d-3)then
           if(dabs(z-1d0).gt.1d-3)then
            fc81=-28d0*y**2/3d0/(x-y)/(y-1d0)/(y-z)
     .          +4d0*x*(7d0*x**2-3d0*x*y+3d0*y**2)*dlog(x)
     .           /3d0/(x-1d0)/(x-y)**2/(x-z)
     .          +4d0*z*(7d0*z**2-3d0*z*y+3d0*y**2)*dlog(z)
     .           /3d0/(z-1d0)/(z-y)**2/(z-x)
     .          -4d0*y**2*(x*(4d0*y**2+18d0*z-11d0*y*(1d0+z))
     .                    +y*(3d0*y**2-11d0*z+4d0*y*(1d0+z)))
     .             *dlog(y)/3d0/(x-y)**2/(y-1d0)**2/(y-z)**2
     .          -4d0*(1d0+y**2)/(x-1d0)/(y-1d0)/(z-1d0)*Sp2(1d0-1d0/y)
     .          +4d0*(x**2+y**2)*Sp2(1d0-x/y)/(x-1d0)/(x-y)/(x-z)
     .          +4d0*(z**2+y**2)*Sp2(1d0-z/y)/(z-1d0)/(z-y)/(z-x)
           else
            fc81=4d0/3d0*(-7d0*y**2/(x-y)/(y-1d0)**2
     .          +(-7d0+3d0*y-3d0*y**2)/(x-1d0)/(y-1d0)**2
     .       +x*(7d0*x**2-3d0*x*y+3d0*y**2)*dlog(x)/(x-1d0)**2/(x-y)**2
     .          -y**2*(2d0*x*(-9d0+2d0*y)+y*(11d0+3d0*y))*dlog(y)
     .                 /(x-y)**2/(y-1d0)**3
     .          +3d0*(x**2+y**2)*Sp2(1d0-x/y)/(x-1d0)**2/(x-y)
     .          +3d0*(-(x-1d0)*(1d0+y**2)*dlog(y)
     .        +(y*(-1d0-2d0*y+y**2)+x*(-1d0+2d0*y+y**2))*Sp2(1d0-1d0/y))
     .                 /(x-1d0)**2/(y-1d0)**2)
           endif
          else
           if(dabs(z-1d0).gt.1d-3)then
            fc81=2d0/3d0*(-11d0/(x-1d0)/(z-1d0)
     .          +14d0*(-z+2d0*x*z-x)/(x-1d0)**2/(z-1d0)**2
     .          +2d0*x*(3d0-3d0*x+7d0*x**2)*dlog(x)/(x-1d0)**3/(x-z)
     .          +2d0*z*(3d0-3d0*z+7d0*z**2)*dlog(z)/(z-1d0)**3/(z-x)
     .          +6d0*(1d0+x**2)*Sp2(1d0-x)/(x-1d0)**2/(x-z)
     .          +6d0*(1d0+z**2)*Sp2(1d0-z)/(z-1d0)**2/(z-x))
           else
            fc81=2d0/9d0/(x-1d0)**4*(32d0-63d0*x+72d0*x**2-41d0*x**3
     .                    +6d0*x*(3d0-3d0*x+7d0*x**2)*dlog(x)
     .                    +18d0*(-1d0+x-x**2+x**3)*Sp2(1d0-x))
           endif
          endif
         else
          if(dabs(y-1d0).gt.1d-3)then
           if(dabs(z-1d0).gt.1d-3)then
            fc81=4d0/3d0*(7d0*y**2/(y-1d0)**2/(y-z)
     .          +(-7d0+3d0*y-3d0*y**2)/(y-1d0)**2/(z-1d0)
     .          -y**2*(3d0*y**2-18d0*z+y*(11d0+4d0*z))*dlog(y)
     .             /(y-1d0)**3/(y-z)**2
     .          +z*(3d0*y**2-3d0*y*z+7d0*z**2)*dlog(z)
     .             /(y-z)**2/(z-1d0)**2
     .          +3d0*((1d0+y**2)*(1d0-z)*dlog(y)
     .          +(y**3+y**2*(z-2d0)-z+y*(2d0*z-1d0))*Sp2(1d0-1d0/y))
     .             /(y-1d0)**2/(z-1d0)**2
     .          +3d0*(y**2+z**2)*Sp2(1d0-z/y)/(z-1d0)**2/(z-y))
           else
            fc81=-2d0/3d0/(y-1d0)**4*(-4d0+33d0*y-12d0*y**2-17d0*y**3
     .                   -3d0*(1d0-5d0*y-11d0*y**2+y**3)*dlog(y)
     .                   +12d0*(y-1d0)*y**2*Sp2(1d0-1d0/y))
           endif
          else
           if(dabs(z-1d0).gt.1d-3)then
            fc81=2d0/9d0/(z-1d0)**4*(32d0-63d0*z+72d0*z**2-41d0*z**3
     .           +6d0*z*(3d0-3d0*z+7d0*z**2)*dlog(z)
     .           +18d0*(-1d0+z-z**2+z**3)*Sp2(1d0-z))
           else
            fc81=1d0/9d0
           endif
          endif
         endif
        else
         if(dabs(x-y).gt.1d-3)then
          if(dabs(x-1d0).gt.1d-3)then
           if(dabs(y-1d0).gt.1d-3)then
            fc81=4d0/3d0*(-6d0*y/(y-1d0)/(y-x)
     .         +x*(7d0*x**2-3d0*x*y+3d0*y**2)*dlog(x)/(x-1d0)/(x-y)**3
     .         +y*(x*(29d0-15d0*y)+y*(y-15d0))/2d0/(y-1d0)**2/(x-y)**2
     .         +(x*y*(-18d0+7d0*y-3d0*y**2)+y**2*(7d0-3d0*y+3d0*y**2)
     .            +x**2*(18d0-18d0*y+7d0*y**2))*dlog(y)*y
     .              /(y-1d0)**3/(y-x)**3
     .         +3d0*(x**2+y**2)*Sp2(1d0-x/y)/(x-1d0)/(x-y)**2
     .         -3d0*(1d0+y**2)*Sp2(1d0-1d0/y)/(x-1d0)/(y-1d0)**2)
           else
            fc81=2d0/9d0/(x-1d0)**4*(32d0-63d0*x+72d0*x**2-41d0*x**3
     .           +6d0*x*(3d0-3d0*x+7d0*x**2)*dlog(x)
     .           +18d0*(-1d0+x-x**2+x**3)*Sp2(1d0-x))
           endif
          else
           if(dabs(y-1d0).gt.1d-3)then
            fc81=2d0/3d0/(y-1d0)**4*(14d0-3d0*y+6d0*y**2-17d0*y**3
     .            +6d0*(1d0+5d0*y+y**2)*dlog(y)
     .            +12d0*y*(-1d0+y**2)*Sp2(1d0-1d0/y))
           else
            fc81=1d0/9d0
           endif
          endif
         else
          if(dabs(x-1d0).gt.1d-3)then
           fc81=-2d0/9d0/(x-1d0)**4*(41d0-72d0*x+63d0*x**2-32d0*x**3
     .         +6d0*(7d0-3d0*x+3d0*x**2)*dlog(x)
     .         +18d0*(-1d0+x-x**2+x**3)*Sp2(1d0-1d0/x))
          else
           fc81=1d0/9d0
          endif
         endif
        endif
       else
        if(dabs(x-y).gt.1d-3)then
         if(dabs(x-1d0).gt.1d-3)then
          if(dabs(y-1d0).gt.1d-3)then
           fc81=4d0/3d0*(7d0*y**2/(x-y)**2/(y-1d0)
     .               +(7d0*x**2-3d0*x*y+3d0*y**2)/(x-1d0)/(x-y)**2
     .               -(3d0*x**2*(-7d0+y)*y+3d0*x*y**2-3d0*y**3
     .                +x**3*(7d0+11d0*y))*dlog(x)/(x-1d0)**2/(x-y)**3
     .               +y**2*(y*(4d0+3d0*y)+x*(-18d0+11d0*y))*dlog(y)
     .                 /(x-y)**3/(y-1d0)**2
     .               -3d0/(x-1d0)**2/(x-y)**2
     .                *((x-1d0)*(x**2+y**2)*dlog(x/y)
     .                +(2d0*x*(y-1d0)*y+x**2*(1d0+y)-y**2*(1d0+y))
     .                  *Sp2(1d0-x/y))
     .                -3d0*(1d0+y**2)*Sp2(1d0-1d0/y)/(x-1d0)**2/(y-1d0))
          else
           fc81=-2d0/3d0/(x-1d0)**4*(17d0-6d0*x+3d0*x**2-14d0*x**3
     .      +6d0*x*(1d0+5d0*x+x**2)*dlog(x)+12d0*(x**2-1d0)*Sp2(1d0-x))
          endif
         else
          if(dabs(y-1d0).gt.1d-3)then
           fc81=-2d0/3d0/(y-1d0)**4*(-4d0+33d0*y-12d0*y**2-17d0*y**3
     .           -3d0*(1d0-5d0*y-11d0*y**2+y**3)*dlog(y)
     .           +12d0*(y-1d0)*y**2*Sp2(1d0-1d0/y))
          else
           fc81=1d0/9d0
          endif
         endif
        else
         if(dabs(x-1d0).gt.1d-3)then
          fc81=-2d0/9d0/(x-1d0)**4*(41d0-72d0*x+63d0*x**2-32d0*x**3
     .         +6d0*(7d0-3d0*x+3d0*x**2)*dlog(x)
     .         +18d0*(-1d0+x-x**2+x**3)*Sp2(1d0-1d0/x))
         else
          fc81=1d0/9d0
         endif
        endif
       endif
      else
       if(dabs(x-z).gt.1d-3)then
        if(dabs(x-1d0).gt.1d-3)then
         if(dabs(z-1d0).gt.1d-3)then
          fc81=2d0/3d0*(-12d0*x/(x-1d0)/(x-z)
     .           -14d0*x*(x-2d0*z+x*z)/(x-1d0)**2/(x-z)**2
     .           +6d0*x*dlog(x)/(x-1d0)/(x-z)
     .           +x/(x-1d0)/(x-z)+2d0*x/(x-1d0)**3/(x-z)**3
     .            *(15d0*z**2+3d0*x**3*(1d0+z)-12d0*x*z*(1d0+z)
     .              +x**2*(4d0-5d0*z+4d0*z**2))*dlog(x)
     .           +2d0*z*(3d0*x**2-3d0*x*z+7d0*z**2)*dlog(z)
     .              /(z-1d0)/(z-x)**3
     .           -6d0*(1d0+x**2)*Sp2(1d0-1d0/x)/(x-1d0)**2/(z-1d0)
     .           +6d0*(x**2+z**2)*Sp2(1d0-z/x)/(x-z)**2/(z-1d0))
         else
          fc81=2d0/3d0/(x-1d0)**4*(14d0-3d0*x+6d0*x**2-17d0*x**3
     .           +6d0*(1d0+5d0*x+x**2)*dlog(x)
     .           +12d0*x*(x**2-1d0)*Sp2(1d0-1d0/x))
         endif
        else
         if(dabs(z-1d0).gt.1d-3)then
          fc81=2d0/9d0/(z-1d0)**4*(32d0-63d0*z+72d0*z**2-41d0*z**3
     .            +6d0*z*(3d0-3d0*z+7d0*z**2)*dlog(z)
     .            +18d0*(-1d0+z-z**2+z**3)*Sp2(1d0-z))
         else
          fc81=1d0/9d0
         endif
        endif
       else
        if(dabs(x-1d0).gt.1d-3)then
         fc81=-2d0/9d0/(x-1d0)**4*(41d0-72d0*x+63d0*x**2-32d0*x**3
     .         +6d0*(7d0-3d0*x+3d0*x**2)*dlog(x)
     .         +18d0*(-1d0+x-x**2+x**3)*Sp2(1d0-1d0/x))
        else
         fc81=1d0/9d0
        endif
       endif
      endif
      return
      end

*********************************************************************

      double precision function fc121(x,y,z)
      
      implicit none
      double precision x,y,z
      if(dabs(x-y).gt.1d-3)then
       if(dabs(x-z).gt.1d-3)then
        if(dabs(y-z).gt.1d-3)then
         if(dabs(x-1d0).gt.1d-3)then
          if(dabs(y-1d0).gt.1d-3)then
           if(dabs(z-1d0).gt.1d-3)then
            fc121=-28d0*y**2/3d0/(x-y)/(y-1d0)/(y-z)
     .        +4d0*x**2*(6d0*x+y)*dlog(x)/3d0/(x-1d0)/(x-y)**2/(x-z)
     .        +4d0*z**2*(6d0*z+y)*dlog(z)/3d0/(z-1d0)/(z-y)**2/(z-x)
     .        -4d0*y**2*dlog(y)/3d0/(x-y)**2/(y-1d0)**2/(y-z)**2
     .         *(x*(6d0*y**2+20d0*z-13d0*y*(1d0+z))
     .          +y*(y**2-13d0*z+6d0*y*(1d0+z)))
           else
            fc121=4d0/3d0*(-7d0*y**2/(x-y)/(y-1d0)**2
     .          +(6d0+y)/(1d0-x)/(y-1d0)**2
     .          +x**2*(6d0*x+y)*dlog(x)/(x-1d0)**2/(x-y)**2
     .          -y**2*(y*(13d0+y)+x*(-20d0+6d0*y))*dlog(y)
     .                            /(x-y)**2/(y-1d0)**3)
           endif
          else
           if(dabs(z-1d0).gt.1d-3)then
            fc121=4d0/3d0*(5d0/2d0/(x-1d0)/(z-1d0)
     .           +7d0*(-z+x*(2d0*z-1d0))/(x-1d0)**2/(z-1d0)**2
     .           +x**2*(1d0+6d0*x)*dlog(x)/(x-1d0)**3/(x-z)
     .           +z**2*(1d0+6d0*z)*dlog(z)/(z-1d0)**3/(z-x))
           else
            fc121=2d0*((x-1d0)*(-11d0+37d0*x-68d0*x**2)
     .           +6d0*x**2*(1d0+6d0*x)*dlog(x))/9d0/(x-1d0)**4
           endif
          endif
         else
          if(dabs(y-1d0).gt.1d-3)then
           if(dabs(z-1d0).gt.1d-3)then
            fc121=4d0/3d0*(7d0*y**2/(y-1d0)**2/(y-z)
     .                  -(6d0+y)/(y-1d0)**2/(z-1d0)
     .                  -y**2*(y**2-20d0*z+y*(13d0+6d0*z))*dlog(y)
     .                        /(y-1d0)**3/(y-z)**2
     .                  +z**2*(y+6d0*z)*dlog(z)/(y-z)**2/(z-1d0)**2)
           else
            fc121=2d0/3d0/(y-1d0)**4*(6d0-37d0*y+14d0*y**2+17d0*y**3
     .                                -2d0*y**2*(20d0+y)*dlog(y))
           endif
          else
           if(dabs(z-1d0).gt.1d-3)then
            fc121=2d0/9d0/(z-1d0)**4*(11d0-48d0*z+105d0*z**2-68d0*z**3
     .     +6d0*z**2*(1d0+6d0*z)*dlog(z))
           else
            fc121=17d0/9d0
           endif
          endif
         endif
        else
         if(dabs(x-y).gt.1d-3)then
          if(dabs(x-1d0).gt.1d-3)then
           if(dabs(y-1d0).gt.1d-3)then
            fc121=2d0/3d0/(x-1d0)/(x-y)**3/(y-1d0)**3*(
     .              2d0*x**2*(y-1d0)**3*(6d0*x+y)*dlog(x)
     .               -(x-1d0)*y*((y-1d0)*(4d0*x*(13d0-6d0*y)*y
     .               +y**2*(-19d0+5d0*y)+x**2*(-33d0+19d0*y))
     .             +2d0*(y**2*(6d0+y)+x*y*(-18d0+3d0*y+y**2)
     .              +x**2*(19d0-18d0*y+6d0*y**2))*dlog(y)))
           else
            fc121=2d0*(11d0-48d0*x+105d0*x**2-68d0*x**3
     .     +6d0*x**2*(1d0+6d0*x)*dlog(x))/9d0/(x-1d0)**4
           endif
          else
           if(dabs(y-1d0).gt.1d-3)then
            fc121=2d0*(12d0+23d0*y-40d0*y**2+5d0*y**3
     .      +2d0*y*(19d0+2d0*y)*dlog(y))/3d0/(y-1d0)**4
           else
            fc121=17d0/9d0
           endif
          endif
         else
          if(dabs(x-1d0).gt.1d-3)then
           fc121=-2d0/9d0/(x-1d0)**4*(68d0-105d0*x+48d0*x**2-11d0*x**3
     .           +6d0*(6d0+x)*dlog(x))
          else
           fc121=17d0/9d0
          endif
         endif
        endif
       else
        if(dabs(x-y).gt.1d-3)then
         if(dabs(x-1d0).gt.1d-3)then
          if(dabs(y-1d0).gt.1d-3)then
           fc121=4d0/3d0/(x-y)**3*(
     .        7d0*(x-y)*y**2/(y-1d0)
     .       +x*((x-1d0)*(6d0*x**2-5d0*x*y-y**2)
     .       -(x*(-18d0+y)*y-2d0*y**2+x**2*(6d0+13d0*y))*dlog(x))
     .                   /(x-1d0)**2
     .       +y**2*(y*(6d0+y)+x*(-20d0+13d0*y))*dlog(y)/(y-1d0)**2)
          else
           fc121=2d0*(5d0-40d0*x+23d0*x**2+12d0*x**3
     .      -2d0*x*(2d0+19d0*x)*dlog(x))/3d0/(x-1d0)**4
          endif
         else
          if(dabs(y-1d0).gt.1d-3)then
           fc121=2d0*(6d0-37d0*y+14d0*y**2+17d0*y**3
     .     -2d0*y**2*(20d0+y)*dlog(y))/3d0/(y-1d0)**4
          else
           fc121=17d0/9d0
          endif
         endif
        else
         if(dabs(x-1d0).gt.1d-3)then
          fc121=-2d0/9d0/(x-1d0)**4*(68d0-105d0*x+48d0*x**2-11d0*x**3
     .           +6d0*(6d0+x)*dlog(x))
         else
          fc121=17d0/9d0
         endif
        endif
       endif
      else
       if(dabs(x-z).gt.1d-3)then
        if(dabs(x-1d0).gt.1d-3)then
         if(dabs(z-1d0).gt.1d-3)then
          fc121=2d0/3d0/(x-1d0)**3/(x-z)**3/(z-1d0)*(
     .     2d0*x*(z-1d0)*(19d0*z**2+x**3*(1d0+z)
     .     -18d0*x*z*(1d0+z)+3d0*x**2*(2d0+z+2d0*z**2))*dlog(x)
     .     +(x-1d0)*(x*(z-1d0)*(5d0*x**3-33d0*z**2+x*z*(52d0+19d0*z)
     .      -x**2*(19d0+24d0*z))-2d0*(x-1d0)**2*z**2*(x+6d0*z)*dlog(z)))
         else
          fc121=2d0*(12d0+23d0*x-40d0*x**2+5d0*x**3
     .        +2d0*x*(19d0+2d0*x)*dlog(x))/3d0/(x-1d0)**4
         endif
        else
         if(dabs(z-1d0).gt.1d-3)then
          fc121=2d0*(11d0-48d0*z+105d0*z**2-68d0*z**3
     .    +6d0*z**2*(1d0+6d0*z)*dlog(z))/9d0/(z-1d0)**4
         else
          fc121=17d0/9d0
         endif
        endif
       else
        if(dabs(x-1d0).gt.1d-3)then
         fc121=-2d0/9d0/(x-1d0)**4*(68d0-105d0*x+48d0*x**2-11d0*x**3
     .           +6d0*(6d0+x)*dlog(x))
        else
         fc121=17d0/9d0
        endif
       endif
      endif
      return
      end

*********************************************************************

      double precision function fc131(x,y,z)
      
      implicit none
      double precision x,y,z
      if(dabs(x-y).gt.1d-3)then
       if(dabs(x-z).gt.1d-3)then
        if(dabs(y-z).gt.1d-3)then
         if(dabs(x-1d0).gt.1d-3)then
          if(dabs(y-1d0).gt.1d-3)then
           if(dabs(z-1d0).gt.1d-3)then
            fc131=-28d0*y/3d0/(x-y)/(y-1d0)/(y-z)
     .        +4d0*x*(6d0*x+y)*dlog(x)/3d0/(x-1d0)/(x-y)**2/(x-z)
     .        +4d0*z*(6d0*z+y)*dlog(z)/3d0/(z-1d0)/(z-y)**2/(z-x)
     .        +4d0*y*dlog(y)/3d0/(x-y)**2/(y-1d0)**2/(y-z)**2
     .         *(x*(y**2-13d0*z+6d0*y*(1d0+z))
     .          +y*(y-8d0*y**2+6d0*z+y*z))
           else
            fc131=4d0/3d0*(-7d0*y/(x-y)/(y-1d0)**2
     .            +(6d0+y)/(1d0-x)/(y-1d0)**2
     .            +x*(6d0*x+y)*dlog(x)/(x-1d0)**2/(x-y)**2
     .            +y*(x*(13d0+y)-2d0*y*(3d0+4d0*y))
     .                             *dlog(y)/(x-y)**2/(y-1d0)**3)
           endif
          else
           if(dabs(z-1d0).gt.1d-3)then
            fc131=4d0/3d0*(5d0/2d0/(x-1d0)/(z-1d0)
     .           +7d0*(x*z-1d0)/(x-1d0)**2/(z-1d0)**2
     .           +x*(1d0+6d0*x)*dlog(x)/(x-1d0)**3/(x-z)
     .           +z*(1d0+6d0*z)*dlog(z)/(z-1d0)**3/(z-x))
           else
            fc131=2d0/9d0/(x-1d0)**4*(
     .         (x-1d0)*(4d0-35d0*x-11d0*x**2)+6d0*x*(1d0+6d0*x)*dlog(x))
           endif
          endif
         else
          if(dabs(y-1d0).gt.1d-3)then
           if(dabs(z-1d0).gt.1d-3)then
            fc131=4d0/3d0*(7d0*y/(y-1d0)**2/(y-z)
     .    -(6d0+y)/(y-1d0)**2/(z-1d0)
     .    -y*(8d0*y**2-y*(z-6d0)-13d0*z)*dlog(y)/(y-1d0)**3/(y-z)**2
     .    +z*(y+6d0*z)*dlog(z)/(y-z)**2/(z-1d0)**2)
           else
            fc131=2d0/3d0/(y-1d0)**4*(
     .  -6d0-29d0*y+34d0*y**2+y**3-2d0*y*(13d0+8d0*y)*dlog(y))
           endif
          else
           if(dabs(z-1d0).gt.1d-3)then
            fc131=2d0/9d0/(z-1d0)**4*(
     .  -4d0+39d0*z-24d0*z**2-11d0*z**3+6d0*z*(1d0+6d0*z)*dlog(z))
           else
            fc131=-5d0/9d0
           endif
          endif
         endif
        else
         if(dabs(x-y).gt.1d-3)then
          if(dabs(x-1d0).gt.1d-3)then
           if(dabs(y-1d0).gt.1d-3)then
            fc131=-2d0/3d0/(x-1d0)/(x-y)**3/(y-1d0)**3*(
     .    -2d0*x*(y-1d0)**3*(6d0*x+y)*dlog(x)
     .    +(x-1d0)*((y-1d0)*(4d0*x*y*(6d0+y)+x**2*(-19d0+5d0*y)
     .                       -y**2*(5d0+9d0*y))
     .             +2d0*(x**2*(6d0+y)+y**3*(6d0+y)
     .                     +x*y*(1d0-21d0*y+6d0*y**2))*dlog(y)))
           else
            fc131=-2d0*(4d0-39d0*x+24d0*x**2+11d0*x**3
     .                 -6d0*x*(1d0+6d0*x)*dlog(x))/9d0/(x-1d0)**4
           endif
          else
           if(dabs(y-1d0).gt.1d-3)then
            fc131=2d0*(31d0-20d0*y-22d0*y**2
     .      +2d0*(6d0+14d0*y+y**2)*dlog(y))/3d0/(y-1d0)**4
           else
            fc131=-5d0/9d0
           endif
          endif
         else
          if(dabs(x-1d0).gt.1d-3)then
           fc131=-2d0/9d0/x/(x-1d0)**4*(11d0+24d0*x-39d0*x**2+4d0*x**3
     .                                  +6d0*x*(6d0+x)*dlog(x))
          else
           fc131=-5d0/9d0
           endif
         endif
        endif
       else
        if(dabs(x-y).gt.1d-3)then
         if(dabs(x-1d0).gt.1d-3)then
          if(dabs(y-1d0).gt.1d-3)then
           fc131=4d0/3d0/(x-y)**3*(
     .        7d0*(x-y)*y/(y-1d0)
     .       +((x-1d0)*(6d0*x**2-5d0*x*y-y**2)
     .        +(-6d0*x**3+13d0*x*y-8d0*x**2*y+y**2)*dlog(x))
     .                   /(x-1d0)**2
     .       +y*(x*(6d0*y-13d0)+y*(8d0*y-1d0))*dlog(y)/(y-1d0)**2)
          else
           fc131=-2d0*(11d0+20d0*x-31d0*x**2
     .          +2d0*(1d0+14d0*x+6d0*x**2)*dlog(x))/3d0/(x-1d0)**4
          endif
         else
          if(dabs(y-1d0).gt.1d-3)then
           fc131=2d0*(-6d0-29d0*y+34d0*y**2+y**3
     .          -2d0*y*(13d0+8d0*y)*dlog(y))/3d0/(y-1d0)**4
          else
           fc131=-5d0/9d0
          endif
         endif
        else
         if(dabs(x-1d0).gt.1d-3)then
          fc131=-2d0/9d0/x/(x-1d0)**4*(11d0+24d0*x-39d0*x**2+4d0*x**3
     .                                  +6d0*x*(6d0+x)*dlog(x))
         else
          fc131=-5d0/9d0
         endif
        endif
       endif
      else
       if(dabs(x-z).gt.1d-3)then
        if(dabs(x-1d0).gt.1d-3)then
         if(dabs(z-1d0).gt.1d-3)then
          fc131=2d0/3d0/(x-1d0)**3/(x-z)**3/(z-1d0)*(
     .    2d0*(x**4*(z-1d0)-21d0*x**2*(z-1d0)*z+6d0*(z-1d0)*z**2
     .        +6d0*x**3*(z**2-1d0)+x*z*(z**2-1d0))*dlog(x)
     .     -(x-1d0)*((z-1d0)*(9d0*x**3+x**2*(5d0-4d0*z)+19d0*z**2
     .        -x*z*(24d0+5d0*z))+2d0*(x-1d0)**2*z*(x+6d0*z)*dlog(z)))
         else
          fc131=2d0*(31d0-20d0*x-11d0*x**2
     .         +2d0*(6d0+14d0*x+x**2)*dlog(x))/3d0/(x-1d0)**4
         endif
        else
         if(dabs(z-1d0).gt.1d-3)then
          fc131=2d0*(-4d0+39d0*z-24d0*z**2-11d0*z**3
     .         +6d0*z*(1d0+6d0*z)*dlog(z))/9d0/(z-1d0)**4
         else
          fc131=-5d0/9d0
         endif
        endif
       else
        if(dabs(x-1d0).gt.1d-3)then
         fc131=-2d0/9d0/x/(x-1d0)**4*(11d0+24d0*x-39d0*x**2+4d0*x**3
     .                                  +6d0*x*(6d0+x)*dlog(x))
        else
         fc131=-5d0/9d0
        endif
       endif
      endif
      return
      end

*********************************************************************

      double precision function fc90(x,y,z,t)
      
      implicit none
      integer D,D2
      double precision x,y,z,t,x1,x2,x3,x4,aux

      D=0
      if(dabs(x-1d0).gt.1d-3)then
       x1=x
       if(dabs(y-1d0).gt.1d-3)then
        x2=y
        if(dabs(z-1d0).gt.1d-3)then
         x3=z
         if(dabs(t-1d0).gt.1d-3)then
          x4=t
         else
          D=1
         endif
        else
         D=1
         if(dabs(t-1d0).gt.1d-3)then
          x3=t
         else
          D=2
         endif
        endif
       else
        D=1
        if(dabs(z-1d0).gt.1d-3)then
         x2=z
         if(dabs(t-1d0).gt.1d-3)then
          x3=t
         else
          D=2
         endif
        else
         D=2
         if(dabs(t-1d0).gt.1d-3)then
          x2=t
         else
          D=3
         endif
        endif
       endif
      else
       D=1
       if(dabs(y-1d0).gt.1d-3)then
        x1=y
        if(dabs(z-1d0).gt.1d-3)then
         x2=z
         if(dabs(t-1d0).gt.1d-3)then
          x3=t
         else
          D=2
         endif
        else
         D=2
         if(dabs(t-1d0).gt.1d-3)then
          x2=t
         else
          D=3
         endif
        endif
       else
        D=2
        if(dabs(z-1d0).gt.1d-3)then
         x1=z
         if(dabs(t-1d0).gt.1d-3)then
          x2=t
         else
          D=3
         endif
        else
         D=3
         if(dabs(t-1d0).gt.1d-3)then
          x1=t
         else
          D=4
         endif
        endif
       endif
      endif

      if(D.eq.4)fc90=-1d0/12d0

      if(D.eq.3)then
       fc90=(-1d0+6d0*x1-3d0*x1**2-2d0*x1**3+6d0*x1**2*dlog(x1))
     .         /6d0/(x1-1d0)**4
      endif

      if(D.eq.2)then
       if(dabs(x1-x2).gt.1d-3)then
        fc90=-(1d0+x1+x2-3d0*x1*x2)/2d0/(x1-1d0)**2/(x2-1d0)**2
     .       +x1**2*dlog(x1)/(x1-1d0)**3/(x1-x2)
     .       +x2**2*dlog(x2)/(x2-1d0)**3/(x2-x1)
       else
        fc90=-(1d0+4d0*x1-5d0*x1**2+2d0*x1*(2d0+x1)*dlog(x1))
     .          /2d0/(x1-1d0)**4
       endif
      endif

      D2=0
      if(D.eq.1)then
       if(dabs(x1-x2).gt.1d-3)then
        if(dabs(x1-x3).gt.1d-3)then
         if(dabs(x2-x3).gt.1d-3)then
          fc90=x1**2*dlog(x1)/(x1-1d0)**2/(x1-x2)/(x1-x3)
     . +x2**2*dlog(x2)/(x2-1d0)**2/(x2-x1)/(x2-x3)
     . +x3**2*dlog(x3)/(x3-x1)/(x3-x2)/(x3-1d0)**2
     . -1d0/(x1-1d0)/(x2-1d0)/(x3-1d0)
         else
          D2=1
         endif
        else
         D2=1
         aux=x1
         x1=x2
         x2=aux
        endif
       else
         D2=1
        if(dabs(x1-x3).gt.1d-3)then
         aux=x1
         x1=x3
         x2=aux
        else
         D2=2
        endif
       endif

       if(D2.eq.2)fc90=(5d0-4d0*x1-x1**2+(2d0+4d0*x1)*dlog(x1))
     .                    /2d0/(x1-1d0)**4
       if(D2.eq.1)fc90=x1**2*dlog(x1)/(x1-1d0)**2/(x1-x2)**2
     .                -2d0*x2**2*dlog(x2)/(x2-1d0)**3/(x2-x1)
     .    +x2*(x2-x1+(x2-2d0*x1)*dlog(x2))/(x1-x2)**2/(x2-1d0)**2
     .    -1d0/(x1-1d0)/(x2-1d0)**2
      endif

      D2=0
      if(D.eq.0)then
       if(dabs(x1-x2).gt.1d-3)then
        if(dabs(x1-x3).gt.1d-3)then
         if(dabs(x1-x4).gt.1d-3)then
          if(dabs(x2-x3).gt.1d-3)then
           if(dabs(x2-x4).gt.1d-3)then
            if(dabs(x3-x4).gt.1d-3)then
             fc90=x1**2*dlog(x1)/(x1-1d0)/(x1-x2)/(x1-x3)/(x1-x4)
     .         +x2**2*dlog(x2)/(x2-1d0)/(x2-x1)/(x2-x3)/(x2-x4)
     .         +x3**2*dlog(x3)/(x3-1d0)/(x3-x2)/(x3-x1)/(x3-x4)
     .         +x4**2*dlog(x4)/(x4-1d0)/(x4-x2)/(x4-x3)/(x4-x1)
            else
             D2=1
            endif
           else
            D2=1
            aux=x2
            x2=x3
            x3=aux
           endif
          else
           D2=1
           if(dabs(x2-x4).gt.1d-3)then
            aux=x2
            x2=x4
            x3=aux
           else
            D2=2
           endif
          endif
         else
          D2=1
          if(dabs(x2-x3).gt.1d-3)then
           aux=x1
           x1=x2
           x2=x3
           x3=aux
          else
           D2=3
          endif
         endif
        else
         D2=1
         if(dabs(x1-x4).gt.1d-3)then
          if(dabs(x2-x4).gt.1d-3)then
           aux=x1
           x1=x2
           x2=x4
           x3=aux
          else
           D2=3
          endif
         else
          D2=2
          aux=x1
          x1=x2
          x2=aux
         endif
        endif
       else
        D2=1
        if(dabs(x1-x3).gt.1d-3)then
         if(dabs(x1-x4).gt.1d-3)then
          if(dabs(x3-x4).gt.1d-3)then
           aux=x1
           x1=x3
           x2=x4
           x3=aux
          else
           D2=3
           x2=x3
          endif
         else
          D2=2
          aux=x1
          x1=x3
          x2=aux
         endif
        else
         if(dabs(x1-x4).gt.1d-3)then
          D2=2
          aux=x1
          x1=x4
          x2=aux
         else
          D2=4
         endif
        endif
       endif

      IF(D2.eq.4)fc90=-(2d0+3d0*x1-6d0*x1**2+x1**3+6d0*x1*dlog(x1))
     .                /6d0/(x1-1d0)**4/x1

      IF(D2.eq.3)fc90=-x1*(-(x1-1d0)*(x1-x2)+(x1**2-2d0*x2+x1*x2)
     .                  *dlog(x1))/(x1-1d0)**2/(x1-x2)**3
     .            +x2*((x1-x2)*(x2-1d0)+(x1*(x2-2d0)+x2**2)*dlog(x2))
     .                   /(x2-1d0)**2/(x1-x2)**3

      IF(D2.eq.2)fc90=(2d0*x1**2*(x2-1d0)**3*dlog(x1)
     .  -(x1-1d0)*(-(x2-1d0)*(-x1**2*(x2-3d0)-4d0*x1*x2+x2**2*(1d0+x2))
     .  +2d0*(x1**2+x1*(x2-3d0)*x2**2+x2**3)*dlog(x2)))
     .     /2d0/(x1-1d0)/(x1-x2)**3/(x2-1d0)**3

      IF(D2.eq.1)fc90=x1**2*dlog(x1)/(x1-1d0)/(x1-x2)/(x1-x3)**2
     .       +x2**2*dlog(x2)/(x2-1d0)/(x2-x1)/(x2-x3)**2
     .       +x3*((x1-x3)*(x2-x3)*(x3-1d0)+(x1*(x2*(x3-2d0)+x3) 
     .                                     +x3*(x2-x3**2))*dlog(x3))
     .        /(x1-x3)**2/(x2-x3)**2/(x3-1d0)**2

      endif

      return
      end

*********************************************************************

      double precision function fc100(x,y,z,t)
      
      implicit none
      integer D,D2
      double precision x,y,z,t,x1,x2,x3,x4,aux

      D=0
      if(dabs(x-1d0).gt.1d-3)then
       x1=x
       if(dabs(y-1d0).gt.1d-3)then
        x2=y
        if(dabs(z-1d0).gt.1d-3)then
         x3=z
         if(dabs(t-1d0).gt.1d-3)then
          x4=t
         else
          D=1
         endif
        else
         D=1
         if(dabs(t-1d0).gt.1d-3)then
          x3=t
         else
          D=2
         endif
        endif
       else
        D=1
        if(dabs(z-1d0).gt.1d-3)then
         x2=z
         if(dabs(t-1d0).gt.1d-3)then
          x3=t
         else
          D=2
         endif
        else
         D=2
         if(dabs(t-1d0).gt.1d-3)then
          x2=t
         else
          D=3
         endif
        endif
       endif
      else
       D=1
       if(dabs(y-1d0).gt.1d-3)then
        x1=y
        if(dabs(z-1d0).gt.1d-3)then
         x2=z
         if(dabs(t-1d0).gt.1d-3)then
          x3=t
         else
          D=2
         endif
        else
         D=2
         if(dabs(t-1d0).gt.1d-3)then
          x2=t
         else
          D=3
         endif
        endif
       else
        D=2
        if(dabs(z-1d0).gt.1d-3)then
         x1=z
         if(dabs(t-1d0).gt.1d-3)then
          x2=t
         else
          D=3
         endif
        else
         D=3
         if(dabs(t-1d0).gt.1d-3)then
          x1=t
         else
          D=4
         endif
        endif
       endif
      endif

      if(D.eq.4)fc100=1d0/12d0

      if(D.eq.3)then
       fc100=(2d0+3d0*x1-6d0*x1**2+x1**3+6d0*x1*dlog(x1))
     .       /6d0/(x1-1d0)**4
      endif

      if(D.eq.2)then
       if(dabs(x1-x2).gt.1d-3)then
        fc100=(-3d0+x1+x2+x1*x2)/2d0/(x1-1d0)**2/(x2-1d0)**2
     .       +x1*dlog(x1)/(x1-1d0)**3/(x1-x2)
     .       +x2*dlog(x2)/(x2-1d0)**3/(x2-x1)
       else
        fc100=(-5d0+4d0*x1+x1**2-2d0*(1d0+2d0*x1)*dlog(x1))
     .          /2d0/(x1-1d0)**4
       endif
      endif

      D2=0
      if(D.eq.1)then
       if(dabs(x1-x2).gt.1d-3)then
        if(dabs(x1-x3).gt.1d-3)then
         if(dabs(x2-x3).gt.1d-3)then
          fc100=x1*dlog(x1)/(x1-1d0)**2/(x1-x2)/(x1-x3)
     . +x2*dlog(x2)/(x2-1d0)**2/(x2-x1)/(x2-x3)
     . +x3*dlog(x3)/(x3-x1)/(x3-x2)/(x3-1d0)**2
     . -1d0/(x1-1d0)/(x2-1d0)/(x3-1d0)
         else
          D2=1
         endif
        else
         D2=1
         aux=x1
         x1=x2
         x2=aux
        endif
       else
         D2=1
        if(dabs(x1-x3).gt.1d-3)then
         aux=x1
         x1=x3
         x2=aux
        else
         D2=2
        endif
       endif

       if(D2.eq.2)fc100=(1d0+4d0*x1-5d0*x1**2+2d0*x1*(2d0+x1)*dlog(x1))
     .                 /2d0/(x1-1d0)**4/x1
       if(D2.eq.1)fc100=x1*dlog(x1)/(x1-1d0)**2/(x1-x2)**2
     .                -2d0*x2*dlog(x2)/(x2-1d0)**3/(x2-x1)
     .    +x2*(-x1+x2-x1*dlog(x2))/(x1-x2)**2/(x2-1d0)**2
     .    -1d0/(x1-1d0)/(x2-1d0)**2
      endif

      D2=0
      if(D.eq.0)then
       if(dabs(x1-x2).gt.1d-3)then
        if(dabs(x1-x3).gt.1d-3)then
         if(dabs(x1-x4).gt.1d-3)then
          if(dabs(x2-x3).gt.1d-3)then
           if(dabs(x2-x4).gt.1d-3)then
            if(dabs(x3-x4).gt.1d-3)then
             fc100=x1*dlog(x1)/(x1-1d0)/(x1-x2)/(x1-x3)/(x1-x4)
     .         +x2*dlog(x2)/(x2-1d0)/(x2-x1)/(x2-x3)/(x2-x4)
     .         +x3*dlog(x3)/(x3-1d0)/(x3-x2)/(x3-x1)/(x3-x4)
     .         +x4*dlog(x4)/(x4-1d0)/(x4-x2)/(x4-x3)/(x4-x1)
            else
             D2=1
            endif
           else
            D2=1
            aux=x2
            x2=x3
            x3=aux
           endif
          else
           D2=1
           if(dabs(x2-x4).gt.1d-3)then
            aux=x2
            x2=x4
            x3=aux
           else
            D2=2
           endif
          endif
         else
          D2=1
          if(dabs(x2-x3).gt.1d-3)then
           aux=x1
           x1=x2
           x2=x3
           x3=aux
          else
           D2=3
          endif
         endif
        else
         D2=1
         if(dabs(x1-x4).gt.1d-3)then
          if(dabs(x2-x4).gt.1d-3)then
           aux=x1
           x1=x2
           x2=x4
           x3=aux
          else
           D2=3
          endif
         else
          D2=2
          aux=x1
          x1=x2
          x2=aux
         endif
        endif
       else
        D2=1
        if(dabs(x1-x3).gt.1d-3)then
         if(dabs(x1-x4).gt.1d-3)then
          if(dabs(x3-x4).gt.1d-3)then
           aux=x1
           x1=x3
           x2=x4
           x3=aux
          else
           D2=3
           x2=x3
          endif
         else
          D2=2
          aux=x1
          x1=x3
          x2=aux
         endif
        else
         if(dabs(x1-x4).gt.1d-3)then
          D2=2
          aux=x1
          x1=x4
          x2=aux
         else
          D2=4
         endif
        endif
       endif

      IF(D2.eq.4)fc100=(1d0-6d0*x1+3d0*x1**2+2d0*x1**3
     .                -6d0*x1**2*dlog(x1))/6d0/(x1-1d0)**4/x1**2

      IF(D2.eq.3)fc100=((x1-1d0)*(x1-x2)+(x1-2d0*x1**2+x2)
     .                  *dlog(x1))/(x1-1d0)**2/(x1-x2)**3
     .            +((x1-x2)*(x2-1d0)-(x1+x2-2d0*x2**2)*dlog(x2))
     .                   /(x2-1d0)**2/(x1-x2)**3

      IF(D2.eq.2)fc100=(2d0*x1*x2*(x2-1d0)**3*dlog(x1)
     .    +(x1-1d0)*((x2-1d0)*(x1**2*(x2+1d0)-4d0*x1*x2**2+
     .    x2**2*(3d0*x2-1d0))-2d0*x2*(x1+x1**2-3d0*x1*x2+x2**2)
     .      *dlog(x2)))/2d0/(x1-1d0)/(x1-x2)**3/(x2-1d0)**3/x2

      IF(D2.eq.1)fc100=x1*dlog(x1)/(x1-1d0)/(x1-x2)/(x1-x3)**2
     .       +x2*dlog(x2)/(x2-1d0)/(x2-x1)/(x2-x3)**2
     .       +((x1-x3)*(x2-x3)*(x3-1d0)+(x1*(x3**2-x2) 
     .                         +x3**2*(1d0+x2-2d0*x3))*dlog(x3))
     .        /(x1-x3)**2/(x2-x3)**2/(x3-1d0)**2

      endif

      return
      end

*********************************************************************

      double precision function fhD0(x)
      
      implicit none
      double precision x
      if(dabs(x-1d0).gt.1d-3)then
      fhD0=x*(47d0*x**2-79d0*x+38d0)/108d0/(x-1d0)**3
     . -x*(3d0*x**3-6d0*x+4d0)*dlog(x)/18d0/(x-1d0)**4
      else
      fhD0=1d0/24d0
      endif
      return
      end

*********************************************************************

      double precision function fhD1(x)
      
      implicit none
      double precision x,Sp2
      if(dabs(x-1d0).gt.1d-3)then
      fhD1=x/81d0*((380d0*x**3-528d0*x**2+72d0*x+128d0)*Sp2(1d0-1d0/x)
     .                           /(x-1d0)**4
     . +(596d0*x**3-672d0*x**2+64d0*x+204d0)*dlog(x)/(x-1d0)**5
     . +(-6175d0*x**3+9138d0*x**2-3927d0*x-764d0)/9d0/(x-1d0)**4)
      else
      fhD1=-169d0/540d0
      endif
      return
      end

*********************************************************************

      double precision function DeD1(x)
      
      implicit none
      double precision x
      if(dabs(x-1d0).gt.1d-3)then
      DeD1=x/81d0*((-352d0*x**3-972d0*x**2+1944d0*x-1052d0)
     .                   /3d0/(x-1d0)**4
     . +(432d0*x**3-456d0*x**2+40d0*x+128d0)*dlog(x)/(x-1d0)**5)
      else
      DeD1=-22d0/135d0
      endif
      return
      end

*********************************************************************

      double precision function hc30(x)
      
      implicit none
      double precision x
      if(dabs(x-1d0).gt.1d-3)then
      hc30=(-6d0*x**3+9d0*x**2-2d0)/9d0/(x-1d0)**4*dlog(x)
     . +(52d0*x**2-101d0*x+43d0)/54d0/(x-1d0)**3
      else
      hc30=-7d0/36d0
      endif
      return
      end

*********************************************************************

      double precision function hc31(x)
      
      implicit none
      double precision x,Sp2
      if(dabs(x-1d0).gt.1d-3)then
      hc31=(32d0*x**3+120d0*x**2-384d0*x+128d0)*Sp2(1d0-1d0/x)
     .                           /81d0/(x-1d0)**4
     . +(-108d0*x**4+1058d0*x**3-898d0*x**2-1098d0*x+710d0)*dlog(x)
     .                           /81d0/(x-1d0)**5
     . +(-304d0*x**3-13686d0*x**2+29076d0*x-12062d0)/729d0/(x-1d0)**4
      else
      hc31=-127d0/810d0
      endif
      return
      end

*********************************************************************

      double precision function qc51(x,y)
      
      implicit none
      double precision x,y
      if(dabs(x-y).gt.1d-3)then
      if(dabs(x-1d0).gt.1d-3)then
      if(dabs(y-1d0).gt.1d-3)then
       qc51=4d0/27d0/(x-y)*((6d0*x**3-9d0*x**2+2d0)*dlog(x)/(x-1d0)**4
     .                  -(6d0*y**3-9d0*y**2+2d0)*dlog(y)/(y-1d0)**4)
     . +(104d0*x**2*y**2-202d0*x*y**2+86d0*y**2-202d0*x**2*y
     .   +380d0*x*y-154d0*y+86d0*x**2-154d0*x+56d0)
     .    /81d0/(x-1d0)**3/(y-1d0)**3
      else
       qc51=(65d0-204d0*x+180d0*x**2-20d0*x**3-21d0*x**4
     .  +12d0*(2d0-9d0*x**2+6d0*x**3)*dlog(x))/81d0/(x-1d0)**5
      endif
      else
      if(dabs(y-1d0).gt.1d-3)then
       qc51=(65d0-204d0*y+180d0*y**2-20d0*y**3-21d0*y**4
     .  +12d0*(2d0-9d0*y**2+6d0*y**3)*dlog(y))/81d0/(y-1d0)**5
      else
       qc51=-4d0/135d0
      endif
      endif
      else
      if(dabs(x-1d0).gt.1d-3)then
       qc51=-8d0/81d0/(x-1d0)**5*(3d0+4d0*x-45d0*x**2+60d0*x**3
     . -22d0*x**4+3d0*x*(4d0-9d0*x+3d0*x**3)*dlog(x))
      else
       qc51=-4d0/135d0
      endif
      endif
      return
      end

*********************************************************************

      double precision function YSM(x)
      
      implicit none
      double precision x
      if(dabs(x-1d0).gt.1d-3)then
      YSM=3d0*x**2/8d0/(x-1d0)**2*dlog(x)+x/8d0*(x-4d0)/(x-1d0)
      else
      YSM=5d0/16d0
      endif
      return
      end

*********************************************************************

      double precision function WSM(x)
      
      implicit none
      double precision x
      if(dabs(x-1d0).gt.1d-3)then
      WSM=(-32d0*x**4+38d0*x**3+15d0*x**2-18d0*x)*dlog(x)
     .                         /18d0/(x-1d0)**4
     .  +(-18d0*x**4+163d0*x**3-259d0*x**2+108d0*x)/36d0/(x-1d0)**3
      else
      WSM=-173d0/216d0
      endif
      return
      end

*********************************************************************
*           1/mb^4.Int[(q^2-4*mmu^2).(1-q^2/mb^2).Sqrt[1-4*mmu^2/q^2].d(q^2),q^2=1..6]
      double precision function IntpropS1(mmu,mb)  
      
      implicit none
      double precision mmu,mb
      IntpropS1=2d0*dsqrt(9d0-6d0*mmu**2)/mb**8
     .   *(mb**4*(3d0-5d0*mmu**2)-4d0*mb**2*(6d0-7d0*mmu**2+mmu**4)
     .     +3d0*(18d0-18d0*mmu**2+mmu**4+mmu**6))
     . -dsqrt(1d0-4d0*mmu**2)/(12d0*mb**8)
     .   *(3d0+mb**4*(6d0-60d0*mmu**2)-18d0*mmu**2+6d0*mmu**4
     .     +36d0*mmu**6-8d0*mb**2*(1d0-7d0*mmu**2+6d0*mmu**4))
     . +2d0*mmu**4/mb**4*(3d0-4d0*(mmu/mb)**2+3d0*(mmu/mb)**4)
     .   *dlog(2d0*(3d0-mmu**2+dsqrt(9d0-6d0*mmu**2))
     .             /(1d0-2d0*mmu**2+dsqrt(1d0-4d0*mmu**2)))
      return
      end

*********************************************************************
*           1/mb^4.Int[q^2.(1-q^2/mb^2).Sqrt[1-4*mmu^2/q^2].d(q^2),q^2=1..6]
      double precision function IntpropP1(mmu,mb)  
      
      implicit none
      double precision mmu,mb
      IntpropP1=-2d0*dsqrt(9d0-6d0*mmu**2)/mb**8
     .   *(-54d0+mb**4*(-3d0+mmu**2)+6d0*mmu**2
     .    +5d0*mmu**4+5d0*mmu**6-4d0*mb**2*(-6d0+mmu**2+mmu**4))
     . +dsqrt(1d0-4d0*mmu**2)/(12d0*mb**8)
     .   *(-3d0+2d0*mmu**2+10d0*mmu**4+60d0*mmu**6
     .    +6d0*mb**4*(-1d0+2d0*mmu**2)
     .    -8d0*mb**2*(-1d0+mmu**2+6d0*mmu**4))
     . -2d0*mmu**4/mb**4*(1d0-4d0*(mmu/mb)**2+5d0*(mmu/mb)**4)
     .   *dlog(2d0*(3d0-mmu**2+dsqrt(9d0-6d0*mmu**2))
     .             /(1d0-2d0*mmu**2+dsqrt(1d0-4d0*mmu**2)))
      return
      end

*********************************************************************
*           1/mb^4.Int[(q^2-4*mmu^2)/(q^2-mA^2).(1-q^2/mb^2).Sqrt[1-4*mmu^2/q^2].d(q^2),q^2=1..6]
      double precision function IntpropS2(mmu,mb,mA)  
      
      implicit none
      double precision mmu,mb,mA,IntpropS1,mAp,Intp
      IF(mA.gt.30d0)then
       IntpropS2=-IntpropS1(mmu,mb)/mA**2
      else
       mAp=mA
       IF(dabs(mA-1d0).lt.1d-3.and.mA.le.1d0)mAp=0.999d0
       IF(dabs(mA-1d0).lt.1d-3.and.mA.gt.1d0)mAp=1.001d0
       IF(dabs(mA-dsqrt(6d0)).lt.1d-3.and.mA.le.dsqrt(6d0))
     .  mAp=2.44849d0
       IF(dabs(mA-dsqrt(6d0)).lt.1d-3.and.mA.gt.dsqrt(6d0))
     .  mAp=2.45049d0
       Intp=3.5513958263617518d0-0.4200909856255438d0*mAp**2
     . +0.01042285250455391d0*mAp**4+mAp**2
     .  *(1d0-0.09131419387829645*mAp**2+0.002084570500910778*mAp**4)
     .                 *dlog(dabs((mAp**2-6d0)/(mAp**2-1d0)))
       IntpropS2=Intp*IntpropS1(mmu,mb)/IntpropS1(0d0,4.68d0)/mb**4
      endif
      return
      end

*********************************************************************
*           1/mb^4.Int[q^2/(q^2-mA^2).(1-q^2/mb^2).Sqrt[1-4*mmu^2/q^2].d(q^2),q^2=1..6]
      double precision function IntpropP2(mmu,mb,mA)  
      
      implicit none
      double precision mmu,mb,mA,IntpropP1,mAp,Intp
      IF(mA.gt.30d0)then
       IntpropP2=-IntpropP1(mmu,mb)/mA**2
      else
       mAp=mA
       IF(dabs(mA-1d0).lt.1d-3.and.mA.le.1d0)mAp=0.999d0
       IF(dabs(mA-1d0).lt.1d-3.and.mA.gt.1d0)mAp=1.001d0
       IF(dabs(mA-dsqrt(6d0)).lt.1d-3.and.mA.le.dsqrt(6d0))
     .  mAp=2.44849d0
       IF(dabs(mA-dsqrt(6d0)).lt.1d-3.and.mA.gt.dsqrt(6d0))
     .  mAp=2.45049d0
       Intp=3.5513958263617518d0-0.4200909856255438d0*mAp**2
     . +0.01042285250455391d0*mAp**4+mAp**2
     .  *(1d0-0.09131419387829645*mAp**2+0.002084570500910778*mAp**4)
     .                 *dlog(dabs((mAp**2-6d0)/(mAp**2-1d0)))
       IntpropP2=Intp*IntpropP1(mmu,mb)/IntpropP1(0d0,4.68d0)/mb**4
      endif
      return
      end

*********************************************************************
*           1/mb^4.Int[(q^2-4*mmu^2)/(q^2-mA^2)/(q^2-mB^2).(1-q^2/mb^2).Sqrt[1-4*mmu^2/q^2].d(q^2),q^2=1..6]
      double precision function IntpropS3(mmu,mb,mA1,mA2,G12)  
      
      implicit none
      double precision mmu,mb,mA1,mA2,G12,IntpropS1,m1,m2,Intp
      IF(min(mA1,mA2).gt.30d0)then
       IntpropS3=IntpropS1(mmu,mb)/mA1**2/mA2**2
      else
       m1=min(mA1,mA2)
       m2=Max(mA1,mA2)
       IF(dabs(m1-1d0).lt.1d-3.and.m1.le.1d0)m1=0.999d0
       IF(dabs(m1-1d0).lt.1d-3.and.m1.gt.1d0)m1=1.001d0
       IF(dabs(m1-dsqrt(6d0)).lt.1d-3.and.m1.le.dsqrt(6d0))
     .    m1=2.44849d0
       IF(dabs(m1-dsqrt(6d0)).lt.1d-3.and.m1.gt.dsqrt(6d0))
     .    m1=2.45049d0
       IF(dabs(m2-1d0).lt.1d-3.and.m2.le.1d0)m2=0.999d0
       IF(dabs(m2-1d0).lt.1d-3.and.m2.gt.1d0)m2=1.001d0
       IF(dabs(m2-dsqrt(6d0)).lt.1d-3.and.m2.le.dsqrt(6d0))
     .    m2=2.44849d0
       IF(dabs(m2-dsqrt(6d0)).lt.1d-3.and.m2.gt.dsqrt(6d0))
     .    m2=2.45049d0
       If(dabs(m1-m2).gt.1d-3)then
        Intp=-0.42009098478536167d0+0.010422852483708186d0
     .                                     *(m1**2+m2**2)
     . +(m1**2*(1d0-0.09131419387829645d0*m1**2
     .                         +0.002084570500910778d0*m1**4)
     .           *dlog(dabs((m1**2-6d0)/(m1**2-1d0)))
     .  -m2**2*(1d0-0.09131419387829645d0*m2**2
     .                         +0.002084570500910778d0*m2**4)
     .           *dlog(dabs((m2**2-6d0)/(m2**2-1d0))))
     .                                     /(m1**2-m2**2)
        IntpropS3=Intp*IntpropS1(mmu,mb)/IntpropS1(0d0,4.68d0)/mb**4
       else                  ! G12=GamA*GamB
        Intp=-0.42009098478536155d0+0.02084570496741637d0*m1**2
     .  +m1*(-1d0+0.09131419387829645d0*m1**2-0.002084570500910778d0
     .                                                        *m1**4
     .   +G12*(-0.09131419387829645d0+0.006253711502732334d0*m1**2))
     .      *(datan((1d0-m1**2)/m1/dsqrt(G12))
     .         -datan((6d0-m1**2)/m1/dsqrt(G12)))/dsqrt(G12)
     .  +(0.5d0-(0.09131419387829645d0+0.001042285250455389d0*G12)
     .          *m1**2+0.003126855751366167d0*m1**4)
     .   *dlog(dabs((36d0+(-12d0+G12)*m1**2+m1**4)
     .             /(1d0+(-2d0+G12)*m1**2+m1**4)))
        IntpropS3=Intp*IntpropS1(mmu,mb)/IntpropS1(0d0,4.68d0)/mb**4
       endif
      endif
      return
      end

*********************************************************************
*           1/mb^4.Int[q^2/(q^2-mA^2)/(q^2-mB^2).(1-q^2/mb^2).Sqrt[1-4*mmu^2/q^2].d(q^2),q^2=1..6]
      double precision function IntpropP3(mmu,mb,mA1,mA2,G12)  
      
      implicit none
      double precision mmu,mb,mA1,mA2,G12,IntpropP1,m1,m2,Intp
      IF(min(mA1,mA2).gt.30d0)then
       IntpropP3=IntpropP1(mmu,mb)/mA1**2/mA2**2
      else
       m1=min(mA1,mA2)
       m2=Max(mA1,mA2)
       IF(dabs(m1-1d0).lt.1d-3.and.m1.le.1d0)m1=0.999d0
       IF(dabs(m1-1d0).lt.1d-3.and.m1.gt.1d0)m1=1.001d0
       IF(dabs(m1-dsqrt(6d0)).lt.1d-3.and.m1.le.dsqrt(6d0))
     .    m1=2.44849d0
       IF(dabs(m1-dsqrt(6d0)).lt.1d-3.and.m1.gt.dsqrt(6d0))
     .    m1=2.45049d0
       IF(dabs(m2-1d0).lt.1d-3.and.m2.le.1d0)m2=0.999d0
       IF(dabs(m2-1d0).lt.1d-3.and.m2.gt.1d0)m2=1.001d0
       IF(dabs(m2-dsqrt(6d0)).lt.1d-3.and.m2.le.dsqrt(6d0))
     .    m2=2.44849d0
       IF(dabs(m2-dsqrt(6d0)).lt.1d-3.and.m2.gt.dsqrt(6d0))
     .    m2=2.45049d0
       If(dabs(m1-m2).gt.1d-3)then
        Intp=-0.42009098478536167d0+0.010422852483708186d0
     .                                     *(m1**2+m2**2)
     . +(m1**2*(1d0-0.09131419387829645d0*m1**2
     .                         +0.002084570500910778d0*m1**4)
     .           *dlog(dabs((m1**2-6d0)/(m1**2-1d0)))
     .  -m2**2*(1d0-0.09131419387829645d0*m2**2
     .                         +0.002084570500910778d0*m2**4)
     .           *dlog(dabs((m2**2-6d0)/(m2**2-1d0))))
     .                                     /(m1**2-m2**2)
        IntpropP3=Intp*IntpropP1(mmu,mb)/IntpropP1(0d0,4.68d0)/mb**4
       else                  ! G12=GamA*GamB
        Intp=-0.42009098478536155d0+0.02084570496741637d0*m1**2
     .  +m1*(-1d0+0.09131419387829645d0*m1**2-0.002084570500910778d0
     .                                                        *m1**4
     .   +G12*(-0.09131419387829645d0+0.006253711502732334d0*m1**2))
     .      *(datan((1d0-m1**2)/m1/dsqrt(G12))
     .         -datan((6d0-m1**2)/m1/dsqrt(G12)))/dsqrt(G12)
     .  +(0.5d0-(0.09131419387829645d0+0.001042285250455389d0*G12)
     .          *m1**2+0.003126855751366167d0*m1**4)
     .   *dlog(dabs((36d0+(-12d0+G12)*m1**2+m1**4)
     .             /(1d0+(-2d0+G12)*m1**2+m1**4)))
        IntpropP3=Intp*IntpropP1(mmu,mb)/IntpropP1(0d0,4.68d0)/mb**4
       endif
      endif
      return
      end

*********************************************************************
*           1/mb^4.Int[(q^2-4*mmu^2).(1-q^2/mb^2).Sqrt[1-4*mmu^2/q^2].d(q^2),q^2=1..6]
      double precision function IntpropSh1(mmu,mb)  
      
      implicit none
      double precision mmu,mb
      IntpropSh1=dsqrt(14.4d0*(14.4d0-4d0*mmu**2))/mb**8
     .   *(-746.496d0+311.04d0*mmu**2-7.2d0*mmu**4-3d0*mmu**6
     .          +mb**4*(-7.2d0+5d0*mmu**2)
     .             +mb**2*(138.24d0-67.2d0*mmu**2-4d0*mmu**4))
     . +dsqrt(mb**2*(mb**2-4d0*mmu**2))/(12d0*mb**2)
     .   *(1d0-22d0*(mmu/mb)**2-42d0*(mmu/mb)**4+36d0*(mmu/mb)**6)
     . -2d0*mmu**4/mb**4*(3d0-4d0*(mmu/mb)**2+3d0*(mmu/mb)**4)
     .   *dlog((14.4d0-2d0*mmu**2+dsqrt(14.4d0*(14.4d0-4d0*mmu**2)))
     .           /(mb**2-2d0*mmu**2+dsqrt(mb**2*(mb**2-4d0*mmu**2))))
      return
      end

*********************************************************************
*           1/mb^4.Int[q^2.(1-q^2/mb^2).Sqrt[1-4*mmu^2/q^2].d(q^2),q^2=14.4..mb^2]
      double precision function IntpropPh1(mmu,mb)  
      
      implicit none
      double precision mmu,mb
      IntpropPh1=dsqrt(14.4d0*(14.4d0-4d0*mmu**2))/mb**8
     .   *(-746.496d0+34.56d0*mmu**2+12d0*mmu**4+5d0*mmu**6
     . +mb**4*(-7.2d0+mmu**2)+mb**2*(138.24d0-9.6d0*mmu**2-4d0*mmu**4))
     . +dsqrt(mb**2*(mb**2-4d0*mmu**2))/(12d0*mb**2)
     .   *(1d0-6d0*(mmu/mb)**2+38d0*(mmu/mb)**4-60d0*(mmu/mb)**6)
     . +2d0*mmu**4/mb**4*(1d0-4d0*(mmu/mb)**2+5d0*(mmu/mb)**4)
     .   *dlog((14.4d0-2d0*mmu**2+dsqrt(14.4d0*(14.4d0-4d0*mmu**2)))
     .           /(mb**2-2d0*mmu**2+dsqrt(mb**2*(mb**2-4d0*mmu**2))))
      return
      end

*********************************************************************
*           1/mb^4.Int[(q^2-4*mmu^2)/(q^2-mA^2).(1-q^2/mb^2).Sqrt[1-4*mmu^2/q^2].d(q^2),q^2=14.4..mb^2]
      double precision function IntpropSh2(mmu,mb,mA)  
      
      implicit none
      double precision mmu,mb,mA,IntpropSh1,mAp,Intp
      IF(mA.gt.30d0)then
       IntpropSh2=-IntpropSh1(mmu,mb)/mA**2
      else
       mAp=mA
       IF(dabs(mA-dsqrt(14.4d0)).lt.1d-3.and.mA.le.dsqrt(14.4d0))
     .        mAp=3.79373d0
       IF(dabs(mA-dsqrt(14.4d0)).lt.1d-3.and.mA.gt.dsqrt(14.4d0))
     .        mAp=3.79573d0
       IF(dabs(mA-4.68d0).lt.1d-3.and.mA.le.4.68d0)mAp=4.679d0
       IF(dabs(mA-4.68d0).lt.1d-3.and.mA.gt.4.68d0)mAp=4.681d0
       Intp=0.29342423289097846d0-0.40120387688455267d0*mAp**2
     . +0.01563928169475445d0*mAp**4
     . +mA**2*(-1d0+0.09131419387829645d0*mAp**2
     .                     -0.002084570500910778d0*mAp**4)
     .  *dlog(dabs((14.4d0- mAp**2)/(21.902399992497596d0-mAp**2)))
       IntpropSh2=Intp*IntpropSh1(mmu,mb)/IntpropSh1(0d0,4.68d0)/mb**4
      endif
      return
      end

*********************************************************************
*           1/mb^4.Int[q^2/(q^2-mA^2).(1-q^2/mb^2).Sqrt[1-4*mmu^2/q^2].d(q^2),q^2=14.4..mb^2]
      double precision function IntpropPh2(mmu,mb,mA)  
      
      implicit none
      double precision mmu,mb,mA,IntpropPh1,mAp,Intp
      IF(mA.gt.30d0)then
       IntpropPh2=-IntpropPh1(mmu,mb)/mA**2
      else
       mAp=mA
       IF(dabs(mA-dsqrt(14.4d0)).lt.1d-3.and.mA.le.dsqrt(14.4d0))
     .        mAp=3.79373d0
       IF(dabs(mA-dsqrt(14.4d0)).lt.1d-3.and.mA.gt.dsqrt(14.4d0))
     .        mAp=3.79573d0
       IF(dabs(mA-4.68d0).lt.1d-3.and.mA.le.4.68d0)mAp=4.679d0
       IF(dabs(mA-4.68d0).lt.1d-3.and.mA.gt.4.68d0)mAp=4.681d0
       Intp=0.29342423289097846d0-0.40120387688455267d0*mAp**2
     . +0.01563928169475445d0*mAp**4
     . +mA**2*(-1d0+0.09131419387829645d0*mAp**2
     .                     -0.002084570500910778d0*mAp**4)
     .  *dlog(dabs((14.4d0- mAp**2)/(21.902399992497596d0-mAp**2)))
       IntpropPh2=Intp*IntpropPh1(mmu,mb)/IntpropPh1(0d0,4.68d0)/mb**4
      endif
      return
      end

*********************************************************************
*           1/mb^4.Int[(q^2-4*mmu^2)/(q^2-mA^2)/(q^2-mB^2).(1-q^2/mb^2).Sqrt[1-4*mmu^2/q^2].d(q^2),q^2=14.4..mb^2]
      double precision function IntpropSh3(mmu,mb,mA1,mA2,G12)  
      
      implicit none
      double precision mmu,mb,mA1,mA2,G12,IntpropSh1,m1,m2,Intp
      IF(min(mA1,mA2).gt.30d0)then
       IntpropSh3=IntpropSh1(mmu,mb)/mA1**2/mA2**2
      else
       m1=min(mA1,mA2)
       m2=Max(mA1,mA2)
       IF(dabs(m1-dsqrt(14.4d0)).lt.1d-3.and.m1.le.dsqrt(14.4d0))
     .        m1=3.79373d0
       IF(dabs(m1-dsqrt(14.4d0)).lt.1d-3.and.m1.gt.dsqrt(14.4d0))
     .        m1=3.79573d0
       IF(dabs(m1-4.68d0).lt.1d-3.and.m1.le.4.68d0)m1=4.679d0
       IF(dabs(m1-4.68d0).lt.1d-3.and.m1.gt.4.68d0)m1=4.681d0
       IF(dabs(m2-dsqrt(14.4d0)).lt.1d-3.and.m2.le.dsqrt(14.4d0))
     .        m2=3.79373d0
       IF(dabs(m2-dsqrt(14.4d0)).lt.1d-3.and.m2.gt.dsqrt(14.4d0))
     .        m2=3.79573d0
       IF(dabs(m2-4.68d0).lt.1d-3.and.m2.le.4.68d0)m2=4.679d0
       IF(dabs(m2-4.68d0).lt.1d-3.and.m2.gt.4.68d0)m2=4.681d0
       If(dabs(m1-m2).gt.1d-3)then
        Intp=-0.40120387768696064d0+0.015639281726033016d0*(m1**2+m2**2)
     .   +(m1**2*(1d0-0.09131419387829645d0*m1**2
     .                          +0.002084570500910778d0*m1**4)
     .    *dlog(dabs((4.68d0**2-m1**2)/(14.4d0-m1**2)))
     .    -m2**2*(1d0-0.09131419387829645d0*m2**2
     .                          +0.002084570500910778d0*m2**4)
     .    *dlog(dabs((4.68d0**2-m2**2)/(14.4d0-m2**2))))/(m1**2-m2**2)
        IntpropSh3=Intp*IntpropSh1(mmu,mb)/IntpropSh1(0d0,4.68d0)/mb**4
       else                  ! G12=GamA*GamB
        Intp=-0.4012038768845528d0+0.0312785633895089d0*m1**2
     .  +m1*(1d0-0.0913142d0*m1**2+0.00208457d0*m1**4
     .      +G12*(0.0913142d0-0.00625371d0*m1**2))
     .      *(datan((4.68d0**2-m1**2)/m1/dsqrt(G12))
     .         -datan((14.4d0-m1**2)/m1/dsqrt(G12)))/dsqrt(G12)
     .  +(0.5d0-(0.0913142d0+0.00104229d0*G12)*m1**2
     .          +0.00312686d0*m1**4)
     .   *dlog(dabs((4.68d0**4+(-2d0*4.68**2+G12)*m1**2+m1**4)
     .             /(14.4d0**2+(-2d0*14.4d0+G12)*m1**2+m1**4)))
        IntpropSh3=Intp*IntpropSh1(mmu,mb)/IntpropSh1(0d0,4.68d0)/mb**4
       endif
      endif
      return
      end

*********************************************************************
*           1/mb^4.Int[q^2/(q^2-mA^2)/(q^2-mB^2).(1-q^2/mb^2).Sqrt[1-4*mmu^2/q^2].d(q^2),q^2=14.4..mb^2]
      double precision function IntpropPh3(mmu,mb,mA1,mA2,G12)  
      
      implicit none
      double precision mmu,mb,mA1,mA2,G12,IntpropPh1,m1,m2,Intp
      IF(min(mA1,mA2).gt.30d0)then
       IntpropPh3=IntpropPh1(mmu,mb)/mA1**2/mA2**2
      else
       m1=min(mA1,mA2)
       m2=Max(mA1,mA2)
       IF(dabs(m1-dsqrt(14.4d0)).lt.1d-3.and.m1.le.dsqrt(14.4d0))
     .        m1=3.79373d0
       IF(dabs(m1-dsqrt(14.4d0)).lt.1d-3.and.m1.gt.dsqrt(14.4d0))
     .        m1=3.79573d0
       IF(dabs(m1-4.68d0).lt.1d-3.and.m1.le.4.68d0)m1=4.679d0
       IF(dabs(m1-4.68d0).lt.1d-3.and.m1.gt.4.68d0)m1=4.681d0
       IF(dabs(m2-dsqrt(14.4d0)).lt.1d-3.and.m2.le.dsqrt(14.4d0))
     .        m2=3.79373d0
       IF(dabs(m2-dsqrt(14.4d0)).lt.1d-3.and.m2.gt.dsqrt(14.4d0))
     .        m2=3.79573d0
       IF(dabs(m2-4.68d0).lt.1d-3.and.m2.le.4.68d0)m2=4.679d0
       IF(dabs(m2-4.68d0).lt.1d-3.and.m2.gt.4.68d0)m2=4.681d0
       If(dabs(m1-m2).gt.1d-3)then
        Intp=-0.40120387768696064d0+0.015639281726033016d0*(m1**2+m2**2)
     .   +(m1**2*(1d0-0.09131419387829645d0*m1**2
     .                          +0.002084570500910778d0*m1**4)
     .    *dlog(dabs((4.68d0**2-m1**2)/(14.4d0-m1**2)))
     .    -m2**2*(1d0-0.09131419387829645d0*m2**2
     .                          +0.002084570500910778d0*m2**4)
     .    *dlog(dabs((4.68d0**2-m2**2)/(14.4d0-m2**2))))/(m1**2-m2**2)
        IntpropPh3=Intp*IntpropPh1(mmu,mb)/IntpropPh1(0d0,4.68d0)/mb**4
       else                  ! G12=GamA*GamB
        Intp=-0.4012038768845528d0+0.0312785633895089d0*m1**2
     .  +m1*(1d0-0.0913142d0*m1**2+0.00208457d0*m1**4
     .      +G12*(0.0913142d0-0.00625371d0*m1**2))
     .      *(datan((4.68d0**2-m1**2)/m1/dsqrt(G12))
     .         -datan((14.4d0-m1**2)/m1/dsqrt(G12)))/dsqrt(G12)
     .  +(0.5d0-(0.0913142d0+0.00104229d0*G12)*m1**2
     .          +0.00312686d0*m1**4)
     .   *dlog(dabs((4.68d0**4+(-2d0*4.68**2+G12)*m1**2+m1**4)
     .             /(14.4d0**2+(-2d0*14.4d0+G12)*m1**2+m1**4)))
        IntpropPh3=Intp*IntpropPh1(mmu,mb)/IntpropPh1(0d0,4.68d0)/mb**4
       endif
      endif
      return
      end

*********************************************************************

      double precision function D0B(x,y,z,t)
      
      implicit none
      integer D
      double precision x,y,z,t,x1,x2,x3,x4,aux

      x1=x**2
      x2=y**2
      x3=z**2
      x4=t**2

      D=0
      if(dabs(x1-x2).gt.1d-3)then
       if(dabs(x1-x3).gt.1d-3)then
        if(dabs(x1-x4).gt.1d-3)then
         if(dabs(x2-x3).gt.1d-3)then
          if(dabs(x2-x4).gt.1d-3)then
           if(dabs(x3-x4).gt.1d-3)then
            D=0
           else
            D=1
           endif
          else
           if(dabs(x3-x4).gt.1d-3)then
            D=1
            aux=x2
            x2=x3
            x3=aux
           else
            D=2
           endif
          endif
         else
          if(dabs(x2-x4).gt.1d-3)then
           D=1
           aux=x3
           x3=x4
           x4=aux
          else
           D=2
          endif
         endif
        else
         if(dabs(x2-x3).gt.1d-3)then
          D=1
          aux=x1
          x1=x3
          x3=aux
         else
          D=4
         endif
        endif
       else
        if(dabs(x1-x4).gt.1d-3)then
         if(dabs(x2-x4).gt.1d-3)then
          D=1
          aux=x1
          x1=x4
          x4=aux
         else
          D=4
         endif
        else
         D=2
         aux=x1
         x1=x2
         x2=aux
        endif
       endif
      else
       if(dabs(x1-x3).gt.1d-3)then
        if(dabs(x1-x4).gt.1d-3)then
         if(dabs(x3-x4).gt.1d-3)then
          D=1
          aux=x1
          x1=x3
          x3=aux
          aux=x2
          x2=x4
          x4=aux
         else
          D=4
          aux=x2
          x2=x3
          x3=aux
         endif
        else
         D=2
         aux=x1
         x1=x3
         x3=x2
         x2=aux
        endif
       else
        if(dabs(x1-x4).gt.1d-3)then
         D=2
         aux=x1
         x1=x4
         x4=x2
         x2=aux
        else
         D=3
        endif
       endif
      endif

      if(D.eq.3)then
       if(dabs(x1).gt.1d-3)then
        D0B=-1d0/6d0/x1**2
       else
        D0B=0d0
       endif
      endif

      if(D.eq.2)then
       if(dabs(x1).gt.1d-3)then
        if(dabs(x2).gt.1d-3)then
         D0B=-((x1+x2)/(x1-x2)**2
     .         +2d0*x1*x2*dlog(x2/x1)/(x1-x2)**3)/2d0/x2
        else
         D0B=0d0
        endif
       else
        if(dabs(x2).gt.1d-3)then
         D0B=-1d0/2d0/x2
        else
         D0B=0d0
        endif
       endif
      endif

      if(D.eq.4)then
       if(dabs(x1).gt.1d-3)then
        if(dabs(x2).gt.1d-3)then
         D0B=2d0/(x1-x2)**2+(x1+x2)*dlog(x2/x1)/(x1-x2)**3
        else
         D0B=0d0
        endif
       else
        D0B=0d0
       endif
      endif

      if(D.eq.1)then
       if(dabs(x1).gt.1d-3)then
        if(dabs(x2).gt.1d-3)then
         if(dabs(x3).gt.1d-3)then
          D0B=x2*dlog(x2/x1)/(x2-x1)/(x2-x3)**2+1d0/(x1-x3)/(x2-x3)
     .       +(x1*x2-x3**2)*dlog(x3/x1)/(x1-x3)**2/(x2-x3)**2
         else
          D0B=0d0
         endif
        else
         if(dabs(x3).gt.1d-3)then
          D0B=1d0/x3/(x3-x1)+dlog(x1/x3)/(x1-x3)**2
         else
          D0B=0d0
         endif
        endif
       else
        if(dabs(x3).gt.1d-3)then
         if(dabs(x2).gt.1d-3)then
          D0B=1d0/x3/(x3-x2)+dlog(x2/x3)/(x2-x3)**2
         else
          D0B=0d0
         endif
        else
         D0B=0d0
        endif
       endif
      endif

      if(D.eq.0)then
       if(dabs(x1).gt.1d-3)then
        if(dabs(x2).gt.1d-3)then
         if(dabs(x3).gt.1d-3)then
          if(dabs(x4).gt.1d-3)then
           D0B=x2*dlog(x2/x1)/(x2-x1)/(x2-x3)/(x2-x4)
     .        +x3*dlog(x3/x1)/(x3-x2)/(x3-x1)/(x3-x4)
     .        +x4*dlog(x4/x1)/(x4-x2)/(x4-x3)/(x4-x1)
          else
           D0B=dlog(x2/x1)/(x2-x1)/(x2-x3)+dlog(x3/x1)/(x3-x2)/(x3-x1)
          endif
         else
          if(dabs(x4).gt.1d-3)then
           D0B=dlog(x4/x1)/(x4-x2)/(x4-x1)+dlog(x2/x1)/(x2-x1)/(x2-x4)
          else
           D0B=0d0
          endif
         endif
        else
         if(dabs(x3).gt.1d-3)then
          if(dabs(x4).gt.1d-3)then
           D0B=dlog(x3/x1)/(x3-x1)/(x3-x4)+dlog(x4/x1)/(x4-x3)/(x4-x1)
          else
           D0B=0d0
          endif
         else
          D0B=0d0
         endif
        endif
       else
        if(dabs(x2).gt.1d-3)then
         if(dabs(x3).gt.1d-3)then
          if(dabs(x4).gt.1d-3)then
           D0B=dlog(x2)/(x2-x3)/(x2-x4)+dlog(x3)/(x3-x2)/(x3-x4)
     .        +dlog(x4)/(x4-x2)/(x4-x3)
          else
           D0B=0d0
          endif
         else
          D0B=0d0
         endif
        else
         D0B=0d0
        endif
       endif
      endif

      return
      end

*********************************************************************

      double precision function D2B(x,y,z,t)
      
      implicit none
      integer D
      double precision x,y,z,t,x1,x2,x3,x4,aux

      x1=x**2
      x2=y**2
      x3=z**2
      x4=t**2

      D=0
      if(dabs(x1-x2).gt.1d-3)then
       if(dabs(x1-x3).gt.1d-3)then
        if(dabs(x1-x4).gt.1d-3)then
         if(dabs(x2-x3).gt.1d-3)then
          if(dabs(x2-x4).gt.1d-3)then
           if(dabs(x3-x4).gt.1d-3)then
            D=0
           else
            D=1
           endif
          else
           if(dabs(x3-x4).gt.1d-3)then
            D=1
            aux=x2
            x2=x3
            x3=aux
           else
            D=2
           endif
          endif
         else
          if(dabs(x2-x4).gt.1d-3)then
           D=1
           aux=x3
           x3=x4
           x4=aux
          else
           D=2
          endif
         endif
        else
         if(dabs(x2-x3).gt.1d-3)then
          D=1
          aux=x1
          x1=x3
          x3=aux
         else
          D=4
         endif
        endif
       else
        if(dabs(x1-x4).gt.1d-3)then
         if(dabs(x2-x4).gt.1d-3)then
          D=1
          aux=x1
          x1=x4
          x4=aux
         else
          D=4
         endif
        else
         D=2
         aux=x1
         x1=x2
         x2=aux
        endif
       endif
      else
       if(dabs(x1-x3).gt.1d-3)then
        if(dabs(x1-x4).gt.1d-3)then
         if(dabs(x3-x4).gt.1d-3)then
          D=1
          aux=x1
          x1=x3
          x3=aux
          aux=x2
          x2=x4
          x4=aux
         else
          D=4
          aux=x2
          x2=x3
          x3=aux
         endif
        else
         D=2
         aux=x1
         x1=x3
         x3=x2
         x2=aux
        endif
       else
        if(dabs(x1-x4).gt.1d-3)then
         D=2
         aux=x1
         x1=x4
         x4=x2
         x2=aux
        else
         D=3
        endif
       endif
      endif

      if(D.eq.3)then
       if(dabs(x1).gt.1d-3)then
        D2B=1d0/3d0/x1
       else
        D2B=0d0
       endif
      endif

      if(D.eq.2)then
       if(dabs(x1).gt.1d-3)then
        if(dabs(x2).gt.1d-3)then
         D2B=-((3*x1-x2)/(x1-x2)**2
     .         +2d0*x1**2*dlog(x2/x1)/(x1-x2)**3)/2d0
        else
         D2B=0d0
        endif
       else
        if(dabs(x2).gt.1d-3)then
         D2B=1d0/2d0/x2
        else
         D2B=0d0
        endif
       endif
      endif

      if(D.eq.4)then
       if(dabs(x1).gt.1d-3)then
        if(dabs(x2).gt.1d-3)then
         D2B=(x1+x2)/(x1-x2)**2+2d0*x1*x2*dlog(x2/x1)/(x1-x2)**3
        else
         D2B=1d0/x1
        endif
       else
        if(dabs(x2).gt.1d-3)then
         D2B=1d0/x2
        else
         D2B=0d0
        endif
       endif
      endif

      if(D.eq.1)then
       if(dabs(x1).gt.1d-3)then
        if(dabs(x2).gt.1d-3)then
         if(dabs(x3).gt.1d-3)then
          D2B=-x2**2*dlog(x2/x1)/(x1-x2)/(x2-x3)**2+x3/(x1-x3)/(x2-x3)
     .   -x3/(x1-x3)**2/(x2-x3)**2*(x2*x3+x1*x3-2d0*x1*x2)*dlog(x3/x1)
         else
          D2B=dlog(x1/x2)/(x1-x2)
         endif
        else
         if(dabs(x3).gt.1d-3)then
          D2B=-1d0/(x1-x3)-x1*dlog(x3/x1)/(x1-x3)**2
         else
          D2B=0d0
         endif
        endif
       else
        if(dabs(x3).gt.1d-3)then
         if(dabs(x2).gt.1d-3)then
          D2B=-1d0/(x2-x3)-x1*dlog(x3/x2)/(x2-x3)**2
         else
          D2B=1d0/x3
         endif
        else
         D2B=0d0
        endif
       endif
      endif

      if(D.eq.0)then
       if(dabs(x1).gt.1d-3)then
        if(dabs(x2).gt.1d-3)then
         if(dabs(x3).gt.1d-3)then
          if(dabs(x4).gt.1d-3)then
           D2B=x2**2*dlog(x2/x1)/(x2-x1)/(x2-x3)/(x2-x4)
     .        +x3**2*dlog(x3/x1)/(x3-x2)/(x3-x1)/(x3-x4)
     .        +x4**2*dlog(x4/x1)/(x4-x2)/(x4-x3)/(x4-x1)
          else
           D2B=x2*dlog(x2/x1)/(x2-x1)/(x2-x3)
     .        +x3*dlog(x3/x1)/(x3-x2)/(x3-x1)
          endif
         else
          if(dabs(x4).gt.1d-3)then
           D2B=x4*dlog(x4/x1)/(x4-x2)/(x4-x1)
     .        +x2*dlog(x2/x1)/(x2-x1)/(x2-x4)
          else
           D2B=dlog(x1/x2)/(x1-x2)
          endif
         endif
        else
         if(dabs(x3).gt.1d-3)then
          if(dabs(x4).gt.1d-3)then
           D2B=x3*dlog(x3/x1)/(x3-x1)/(x3-x4)
     .        +x4*dlog(x4/x1)/(x4-x3)/(x4-x1)
          else
           D2B=dlog(x1/x3)/(x1-x3)
          endif
         else
          if(dabs(x4).gt.1d-3)then
           D2B=dlog(x1/x4)/(x1-x4)
          else
           D2B=0d0
          endif
         endif
        endif
       else
        if(dabs(x2).gt.1d-3)then
         if(dabs(x3).gt.1d-3)then
          if(dabs(x4).gt.1d-3)then
           D2B=x2*dlog(x2)/(x2-x3)/(x2-x4)
     . +x3*dlog(x3)/(x3-x2)/(x3-x4)+x4*dlog(x4)/(x4-x2)/(x4-x3)
          else
           D2B=dlog(x2/x3)/(x2-x3)
          endif
         else
          if(dabs(x4).gt.1d-3)then
           D2B=dlog(x2/x4)/(x2-x4)
          else
           D2B=0d0
          endif
         endif
        else
         if(min(x3,x4).gt.1d-3)then
          D2B=dlog(x3/x4)/(x3-x4)
         else
          D2B=0d0
         endif
        endif
       endif
      endif

      return
      end

**********************************************************************

* MSbar 5/6-flavour evolution for the Strong Coupling Constant
      double precision function asf(x)
      
*       ALPHA_S (x) - NB it uses 5 flavors for x<mt and 6 if x>=mt
      implicit none
      double precision x,pi,asc,fn,b0,b1,vvv,b0t,b1t,ast
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      pi=4d0*datan(1d0)
      asc=ALSMZ
      
      fn=5d0
      b0=11d0-2d0*fn/3d0
      b1=102d0-38d0*fn/3d0
      vvv=1d0-b0*asc/(2d0*pi)*dlog(MZ/x)
      asf=asc/vvv*(1d0-b1/b0*asc/(4d0*pi*vvv)*dlog(vvv))

      if(x.gt.MT)then
       vvv=1d0-b0*asc/(2d0*pi)*dlog(MZ/MT)
       ast=asc/vvv*(1d0-b1/b0*asc/(4d0*pi*vvv)*dlog(vvv))
       b0t=b0-2d0/3d0
       b1t=b1-38d0/3d0
       vvv=1d0-b0t*ast/(2d0*pi)*dlog(MT/x)
       asf=ast/vvv*(1d0-b1t/b0t*ast/(4d0*pi*vvv)*dlog(vvv))
      endif
      return
      end
      
****************************************************************

      double precision function f(x)
      
      implicit none
      integer i
      double precision b(12),x,z,cCc,sum
      z=-dlog(1d0-x)
      b(1)=-.5d0
      b(2)=1d0/6d0
      b(3)=0d0
      b(4)=-1d0/30d0
      b(5)=0d0
      b(6)=1d0/42d0
      b(7)=0d0
      b(8)=-1d0/30d0
      b(9)=0d0
      b(10)=5d0/66d0
      b(11)=0d0
      b(12)=-691d0/2730d0
      cCc=z
      sum=z
      do i=1,12
      cCc=cCc*z/(i+1d0)
      sum=sum+b(i)*cCc
      enddo
      f=sum
      end
            
****************************************************************
* Dilogarithm
      double precision function sp2(x)
*       Li_2(x)
      implicit none
      double precision x,pi,f,bla
      external f
      pi=4d0*datan(1d0)
      bla=0d0
      if(x.ge.-1d0.and.x.le.0.5d0)then
       bla=f(x)
      else if(x.gt.0.5d0.and.x.lt.1d0)then
       bla=-f(1d0-x)+pi**2/6d0-dlog(x)*dlog(1d0-x)
      else if(x.lt.-1d0)then
       bla=-f(1d0/x)-pi**2/6d0-.5d0*(dlog(-x))**2
*      else if(x.gt.-500d0.and.x.lt.-10d0)then
*       bla=1.50616d0+0.153765d0*x-0.0000484249d0*x**2
*     .      -2.69934d-8*x**3
*     .      -1.97807d0*dlog(dabs(x))-0.0245271d0*x*dlog(x)
      else if(x.ge.1d0)then
      bla=0d0
*       write(6,*)'error in dilog',x
      endif
      sp2=bla
      return
      end

****************************************************************

      double precision function sgn(x)
      
      implicit none
      double precision x
      if(x.ge.0d0)then
      sgn=1d0
      else
      sgn=-1d0
      endif
      end
      
********************************************************************

      double precision function ff1(x)

      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1d-3)then
      d=1d0/(x-1d0)**3
      dg=dlog(x)/(x-1d0)**4
      ff1=x*(7d0-5d0*x-8d0*x**2)*d/24d0
     .      +x**2*(3d0*x-2d0)*dg/4d0
      else
      ff1=-5d0/48d0
      endif
      return
      end

*********************************************************************

      double precision function fg1(x)

      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1d-2)then
      d=1d0/(x-1d0)**3
      dg=dlog(x)/(x-1d0)**4
      fg1=x*(2d0+5d0*x-x**2)*d/8d0-3d0*x**2*dg/4d0
      else
      fg1=-1d0/16d0
      endif
      return
      end

*********************************************************************

      double precision function ff2(x)

      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1d-2)then
      d=1d0/(x-1d0)**2
      dg=dlog(x)/(x-1d0)**3
      ff2=x*(3d0-5d0*x)*d/12d0+x*(3d0*x-2d0)*dg/6d0
      else
      ff2=-7d0/36d0
      endif
      return
      end

*********************************************************************

      double precision function fg2(x)

      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1d-2)then
      d=1d0/(x-1d0)**2
      dg=dlog(x)/(x-1d0)**3
      fg2=x*(3d0-x)*d/4d0-x*dg/2d0
      else
      fg2=-1d0/6d0
      endif
      return
      end

*********************************************************************

      double precision function ff3(x)

      implicit none
      double precision x,ff1,ff2
      ff3=2d0/3d0*(1d0-1d0/x)*ff1(x)+ff2(x)+23d0/36d0
      return
      end
*********************************************************************

      double precision function fg3(x)

      implicit none
      double precision x,fg1,fg2
      fg3=2d0/3d0*(1d0-1d0/x)*fg1(x)+fg2(x)+1d0/3d0
      return
      end

*********************************************************************

      double precision function esm(x)
      
      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1d-3)then
      d=1d0/(x-1d0)**3
      dg=dlog(x)/(x-1d0)**4
      esm=x*(-18d0+11d0*x+x**2)*d/12d0
     .      +x**2*(15d0-16d0*x+4d0*x**2)*dg/6d0
      esm=esm-2d0/3d0*dlog(x)
      else
      esm=43d0/72d0
      endif
      return
      end

*********************************************************************

      double precision function eh(x)
      
      implicit none
      double precision x,d,dg
      if(dabs(x-1d0).gt.1d-3)then
      d=1d0/(x-1d0)**3
      dg=dlog(x)/(x-1d0)**4
      eh=x*(16d0-29d0*x+7d0*x**2)*d/36d0+x*(3d0*x-2d0)
     .       *dg/6d0
      else
      eh=1d0/8d0
      endif
      return
      end
*********************************************************************

      double precision function af(x)
      
      implicit none
      double precision x,x2,x3,lnx,pi
      pi=4d0*datan(1d0)
      x2=x**2
      x3=x**3
      lnx=dlog(x)
      af=16d0/9d0*((5d0/2d0-1d0/3d0*pi**2-3d0*1.2021d0
     .     +(5d0/2d0-3d0/4d0*pi**2)*lnx+lnx**2/4d0
     .       +lnx**3/12d0)*x
     .    +(7d0/4d0+2d0/3d0*pi**2-1d0/2d0*pi**2*lnx-lnx**2/4d0
     .     +lnx**3/12d0)*x2
     .    +(-7d0/6d0-pi**2/4d0+2d0*lnx-3d0/4d0*lnx**2)*x3)
      return
      end

*********************************************************************

      double precision function afim(x)
      
      implicit none
      double precision x,x2,x3,lnx,pi
      pi=4d0*datan(1d0)
      x2=x**2
      x3=x**3
      lnx=dlog(x)
      afim=16d0/9d0*pi*((2d0-pi**2/6d0+lnx/2d0+lnx**2/2d0)*x
     .  +(1d0/2d0-pi**2/6d0-lnx+lnx**2/2d0)*x2+x3)
      return
      end

*********************************************************************

      double precision function bf(x)
      
      implicit none
      double precision x,x2,x3,lnx,pi
      pi=4d0*datan(1d0)
      x2=x**2
      x3=x**3
      lnx=dlog(x)
      bf=-8d0/9d0*((-3d0+pi**2/6d0-lnx)*x
     .    +(1d0/2d0+pi**2-2d0*lnx-lnx**2/2d0)*x2
     .    +(-25d0/12d0-pi**2/9d0-19d0/18d0*lnx+2d0*lnx**2)*x3
     .    -2d0/3d0*pi**2*x**(3d0/2d0))
      return
      end

*********************************************************************

      double precision function bfim(x)
      
      implicit none
      double precision x,x2,x3,lnx,pi
      pi=4d0*datan(1d0)
      x2=x**2
      x3=x**3
      lnx=dlog(x)
      bfim=-8d0/9d0*pi*(-x+(1d0-2d0*lnx)*x2
     .       +(-10d0/9d0+4d0/3d0*lnx)*x3)
      return
      end
            
*********************************************************************

      double precision function ffh(x,tanb)
      
      implicit none
      double precision x,x2,x3,x4,x5,dd3,dd4,dd5,lnx,ln2x,spx
      double precision sp2,sum1,sum2,tanb
      x2=x**2
      x3=x**3
      x4=x**4
      x5=x**5
      dd3=(x-1d0)**3
      dd4=(x-1d0)**4
      dd5=(x-1d0)**5
      lnx=dlog(x)
      ln2x=lnx**2
      spx=sp2(1d0-1d0/x)
      sum1=4d0*(-3d0+7d0*x-2d0*x2)/(3d0*dd3)*spx+
     .       (8d0-14d0*x-3d0*x2)/(3d0*dd4)*ln2x+
     .       2d0*(-3d0-x+12d0*x2-2d0*x3)/(3d0*dd4)*lnx
     .       +(7d0-13d0*x+2d0*x2)/dd3
      sum2=x*(18d0-37d0*x+8d0*x2)/dd4*spx
     .      +x*(-14d0+23d0*x+3d0*x2)/dd5*ln2x
     .      +(-50d0+251d0*x-174d0*x2-192d0*x3+21d0*x4)
     .      /(9d0*dd5)*lnx+
     .      (797d0-5436d0*x+7569d0*x2-1202d0*x3)/(108d0*dd4)
      ffh=-4d0/3d0*x*sum1+2d0/9d0*x/(tanb**2)*sum2
      return
      end

***********************************************************************

      double precision function fgh(x,tanb)
      
      implicit none
      double precision x,x2,x3,x4,dd2,dd3,dd4,dd5,lnx,ln2x,sp2
      double precision spx,sum1,sum2,tanb
      x2=x**2
      x3=x**3
      x4=x**4      
      dd2=(x-1d0)**2
      dd3=(x-1d0)**3
      dd4=(x-1d0)**4
      dd5=(x-1d0)**5
      lnx=dlog(x)
      ln2x=lnx**2
      spx=sp2(1d0-1d0/x)
      sum1=(-36d0+25d0*x-17d0*x**2)/(2d0*dd3)*spx+(19d0+17d0*x)
     .     /(dd4)*ln2x+(-3d0-187d0*x+12d0*x2-14d0*x3)/(4d0*dd4)*lnx
     .      +3d0*(143d0-44d0*x+29d0*x2)/(8d0*dd3)
      sum2=x*(30d0-17d0*x+13d0*x2)/dd4*spx
     .      -x*(31d0+17d0*x)/dd5*ln2x
     .      +(-226d0+817d0*x+1353d0*x2+318d0*x3+42d0*x4)
     .       /(36d0*dd5)*lnx
     .      +(1130d0-18153d0*x+7650d0*x2-4451d0*x3)/(216d0*dd4)
      fgh=-x/3d0*sum1+x/(6d0*tanb**2)*sum2
      return
      end

**********************************************************************

* For effective Yukawa Couplings
      double precision function h2(x,y)
      
      implicit none
      double precision x,y
      if(dabs(x-y).gt.1d-2)then
        if(dabs(x-1d0).lt.1d-2)then
             h2=1d0/(-1d0+y)-y*dlog(y)/(-1d0+y)**2
         else
            if(dabs(y-1d0).lt.1d-2)then
             h2=1d0/(-1d0+x) - x*dlog(x)/(-1d0+x)**2
            else
             h2=x*dlog(x)/((1d0-x)*(x-y))+
     .      y*dlog(y)/((1d0-y)*(-x+y))
            endif
         endif
      else
         if(dabs(x-1d0).lt.1d-2)then
            h2=-1d0/2d0
         else
            h2=1d0/(1d0 - x) + dlog(x)/(-1d0 + x)**2
         endif
      endif
      return
      end

*******************************************************************

      double precision function BB0(x,y,z,t)

      implicit none
      double precision x,y,z,t,NMB0
      BB0=1d0-NMB0(x,y,z,t)
      return
      end

*******************************************************************

      double precision function BB1(x,y,z,t)

      implicit none
      double precision x,y,z,t,BB0
      if(dabs(y-z).lt.1d-5)then
      BB1=1d0/2d0*BB0(x,y,z,t)
      else
      BB1=1d0/2d0*BB0(x,y,z,t)-(y+z)/(4d0*(y-z))
     .      +y*z/(2d0*(y-z)**2)*dlog(y/z)
      endif
      return
      end

*******************************************************************

      double precision function S0(x)

      implicit none
      double precision x
      S0=x*(1d0/4d0+9d0/(4d0*(1d0-x))-3d0/(2d0*(1d0-x)**2)
     .       -3d0*x**2*dlog(x)/(2d0*(1d0-x)**3))
      return
      end
      
*******************************************************************

      double precision function RUNMass(y,x)

      implicit none
      double precision x,y
      DOUBLE PRECISION PI,nf,ge0,ge1,b0,b1,aux,aux1,asf
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      PI=4d0*DATAN(1d0)

*       * x: renormalization scale
*       * y: quark running mass at 2 GeV;
*      Running Masses at 2 GeV (PDG 2006):
*      - Mu= 1.5 to 3.0 MeV
*      - Md= 3. to 7. MeV
*      - Ms= 95 +/- 25 MeV
*      - Mc= 1.25 +/- 0.09 GeV
*      - Mb= 4.20 +/- 0.07 GeV

      aux=y

       nf=4d0
       ge0=8d0
       ge1=404d0/3d0-38d0/3d0*nf
       b0=11d0-2d0*nf/3d0
       b1=102d0-38d0/3d0*nf
       aux1=2d0
      
       if(x.gt.MBP)then
        aux=aux*(asf(MBP)/asf(MC))**(ge0/(2d0*b0))
     .      *(1d0+ge0/(8d0*pi*b0)*(ge1/ge0-b1/b0)
     .           *(asf(MBP)-asf(aux1)))

        nf=5d0
        ge0=8d0
        ge1=404d0/3d0-38d0/3d0*nf
        b0=11d0-2d0*nf/3d0
        b1=102d0-38d0/3d0*nf
        aux1=MBP
      
        if(x.gt.MT)then
         aux=aux*(asf(MT)/asf(MBP))**(ge0/(2d0*b0))
     .      *(1d0+ge0/(8d0*pi*b0)*(ge1/ge0-b1/b0)
     .           *(asf(MT)-asf(aux1)))

         nf=6d0
         ge0=8d0
         ge1=404d0/3d0-38d0/3d0*nf
         b0=11d0-2d0*nf/3d0
         b1=102d0-38d0/3d0*nf
         aux1=MT
        endif
      
       endif
      
      RUNMass=aux*(asf(x)/asf(aux1))**(ge0/(2d0*b0))
     .      *(1d0+ge0/(8d0*pi*b0)*(ge1/ge0-b1/b0)
     .           *(asf(x)-asf(aux1)))
      
      return
      end
