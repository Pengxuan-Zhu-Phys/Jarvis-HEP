      SUBROUTINE MAGNMU_CPV(PROB)

*      Computation of the Muon Anomalous Moment
*
*       - Literature:
*      [1]:  F. Jegerlehner, "Essentials of the Muon g-2.",
*     Acta Phys.Polon.B38:3021,2007, hep-ph/0703125
*      [2]:  J. Bijnens and J. Prades, "The hadronic light-by-light
*      contribution to the muon anomalous magnetic moment: Where
*      DO we stand?", Mod. Phys. Lett. A 22 (2007) 767,
*      arXiv:hep-ph/0702170
*      [3]:  A. Czarnecki, W. J. Marciano, A. Vainshtein," Refinements in
*     electroweak contributions to the muon anomalous magnetic moment."
*     Phys.Rev.D67:073006,2003, Erratum-ibid.D73:119901,2006, hep-ph/0212229
*      [4]:  S. P. Martin, J. D. Wells
*     "Muon anomalous magnetic dipole moment in supersymmetric theories."
*     Phys.Rev.D64:035003,2001, hep-ph/0103067
*      [5]:  J. P. Leveille, "The Second Order Weak Correction to (G-2) of the
*     Muon in Arbitrary Gauge Models."
*     Nucl.Phys.B137:63,1978
*      [6]:  J. F. Gunion, D. Hooper, B. McElrath
*     "Light neutralino dark matter in the NMSSM"
*     Phys.Rev.D73:015011,2006, hep-ph/0509024
*      [7]:  G. Degrassi, G.F. Giudice, "QED logarithms in the electroweak
*     corrections to the muon anomalous magnetic moment."
*     Phys.Rev.D58:053007,1998, hep-ph/9803384
*      [8]:  S. Heinemeyer, D. Stockinger, G. Weiglein,
*     "Electroweak and supersymmetric two-loop corrections to (g-2)(mu)."
*     Nucl.Phys.B699:103-123,2004, hep-ph/0405255
*      [9]:  K. Cheung, C.-H. Chou, O.C.W. Kong, "Muon anomalous magnetic
*     moment, two Higgs doublet model, and supersymmetry."
*     Phys.Rev.D64:111301,2001, hep-ph/0103183
*      [10]: A. Arhrib, S. Baek, "Two loop Barr-Zee type
*     contributions to (g-2)(muon) in the MSSM."
*     Phys.Rev.D65:075002,2002, hep-ph/0104225
*      [11]: D. Stockinger, "The Muon Magnetic Moment and Supersymmetry."
*     J.Phys.G34:R45-R92,2007, hep-ph/0609168

      IMPLICIT NONE

      INTEGER I,K
      DOUBLE PRECISION PROB(*)
      DOUBLE PRECISION aux,Ymu,Pi,QNP,mtq,mbq,csqmag
      DOUBLE PRECISION amuexp,amuqed,damuqed,amuhadlo,amuhadnlo
      DOUBLE PRECISION amulbl,amuweaklo,amufermnlo,amubosnlo,damusm
      DOUBLE PRECISION errexp,errqed,errhadlo,errhadnlo
      DOUBLE PRECISION errlbl,errfermnlo,errbosnlo,errsm
      DOUBLE PRECISION amu2Lbos,amu2LHf,amu2LCh,amu2LSF,amu2L,amuerr,
     .   amuChar,amuNeutr,amuHiggs,amu1L
      DOUBLE PRECISION delmagmu,damumin,damumax,amuthmax,amuthmin
      DOUBLE PRECISION fc1,fc2,fn1,fn2,fhs,fhp,fhc,runmb,asf,fS,fPS,fSF
      DOUBLE PRECISION lambdHT(5,2),lambdHB(5,2),lambdHL(5,2)

      DOUBLE PRECISION ALEM0
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION mt,mb,mtau,mmuon,mel,MS,MC,MBP,MPI,MSTRANGE
      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION MSU2(2),MSD2(2),MSE2(2),MSNE2,MSMU2(2),
     . UMU(2,2,2)
      DOUBLE PRECISION MCH2(2),U(2,2,2),V(2,2,2)
      DOUBLE PRECISION MNEU(5),N(5,5,2)
      DOUBLE PRECISION MHC2,XC(2,2),MH2(5),XH(5,5),MA2
      DOUBLE PRECISION MST2(2),UT(2,2,2),MSB2(2),UB(2,2,2),MSL2(2),
     . UTAU(2,2,2),MSNT2
      DOUBLE PRECISION GRHSTST(5,2,2),GRHSBSB(5,2,2),GRHSLSL(5,2,2),
     . GRHSUSU(5,2,2),GRHSDSD(5,2,2),GRHSESE(5,2,2),GRHSNSN(5)
      DOUBLE PRECISION GIHSTST(5,2,2),GIHSBSB(5,2,2),GIHSLSL(5,2,2)
      DOUBLE PRECISION GRHCSTSB(2,2),GRHCSNSL(2),GRHCSUSD(2,2),
     . GRHCSNSE(2),GIHCSTSB(2,2),GIHCSNSL(2)
      DOUBLE PRECISION COH0CH(5,2,2,2),COH0NEU(5,5,5,2),
     . COHPNEUCHM(2,5,2,2),COHMNEUCHP(2,5,2,2)

      COMMON/ALEM0/ALEM0
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/SMFERM/mt,mb,mtau,mmuon,mel,MS,MC,MBP,MPI,MSTRANGE
      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/SFERM1SPEC/MSU2,MSD2,MSE2,MSNE2,MSMU2,UMU
      COMMON/CHASPEC/MCH2,U,V
      COMMON/NEUSPEC/MNEU,N
      COMMON/HISPEC/MHC2,XC,MH2,XH,MA2
      COMMON/SFERM3SPEC/MST2,UT,MSB2,UB,MSL2,UTAU,MSNT2
      COMMON/HISFCOUP/GRHSTST,GRHSBSB,GRHSLSL,GRHSUSU,GRHSDSD,
     . GRHSESE,GRHSNSN,GIHSTST,GIHSBSB,GIHSLSL,GRHCSTSB,GRHCSNSL,
     . GRHCSUSD,GRHCSNSE,GIHCSTSB,GIHCSNSL
      COMMON/HINOCOUP/COH0CH,COH0NEU,COHPNEUCHM,COHMNEUCHP
      COMMON/MAGMU/delmagmu,damumin,damumax,amuthmax,amuthmin


      pi=4d0*datan(1d0)
      Ymu=MMUON/vdq

*   A: Evaluation of the Discrepancy between Experiment and SM Contributions
*             (except for 2-loop bosonic: (-2.2 +/- 0.2)d-10: see below)

*       -> World Average after Fermilab: 0.00116592061(41)
      amuexp=11659206.1d-10
      errexp=(4.1d-10)**2
      
*       -> QED Contributions (4 loops):
      amuqed=11658471.8113d-10
      errqed=(0.0162d-10)**2
      
*      -> Difference Exp.-QED:
      damuqed=amuexp-amuqed
      
*      -> LO Hadronic contr. from [1]:
      amuhadlo=692.1d-10
      errhadlo=(5.6d-10)**2
      
*      -> NLO Hadronic contr. from [1]:
      amuhadnlo=-10.03d-10
      errhadnlo=(.22d-10)**2
      
*      -> Light-by-Light contr. from [2]:
      amulbl=11.0d-10
      errlbl=(4.0d-10)**2
      
*      -> Weak SM contr., 1 Loop (error negl.)
      amuweaklo=19.482d-10
      
*      -> Weak SM fermionic 2 loop contr. from [3]
      amufermnlo=-1.512d-10
      errfermnlo=(.1d-10)**2
      
*      -> Weak SM bosonic 2 loop contr., e.g. from [8]
*       (Serves just to estimate pure SM contribution.
*      It will be deduced from the 2-loop bosonic contr. below
*      in order to avoid double counting.)
      amubosnlo=GF*MMUON**2*ALEMMZ/(72d0*dsqrt(2d0)*pi**3)
     .      *(107d0+23d0*(1d0-4d0*s2TW)**2)*dlog(MMUON/MW)
      errbosnlo=(0.2d-10)**2
      
*      -> Difference Exp.-SM:
      damusm=damuqed-amuhadlo-amuhadnlo-amulbl-amuweaklo-amufermnlo
     .    -amubosnlo

*      -> 1-sigma error of this difference:
      errsm=dsqrt(errexp+errqed+errhadlo+errhadnlo+errlbl+errfermnlo
     .     +errbosnlo)

*      -> 2-sigma lower bound to the required SUSY contr.:
      damumin=damusm-2d0*errsm

*      -> 2-sigma upper bound to the required SUSY contr.:
      damumax=damusm+2d0*errsm


*   B: 1-loop contributions
*       - 1-Loop Chargino/SNeutrino Contribution (cf. [4]):
      aux=0d0
      DO k=1,2
      aux=aux+MMUON/(12d0*MSNE2)
     .     *(g2q*(V(k,1,1)**2+V(k,1,2)**2)
     .      +Ymu**2*(U(k,2,1)**2+U(k,2,2)**2))*Fc1(MCH2(k)/MSNE2)
     .     +2d0*dsqrt(MCH2(k))/(3d0*MSNE2)*dsqrt(g2q)*Ymu
     .       *(-V(k,1,1)*U(k,2,1)+V(k,1,2)*U(k,2,2))
     .                                        *Fc2(MCH2(k)/MSNE2)
      ENDDO
      amuChar=MMUON/(16d0*pi**2)*aux

*       - 1-Loop Neutralino/SMuon Contribution (cf. [4]):
      aux=0d0
      DO i=1,5
       DO k=1,2
       aux=aux-MMUON/(12d0*MSMU2(k))*((1d0/dsqrt(2d0)
     .  *(dsqrt(g1q)*N(i,1,1)+dsqrt(g2q)*N(i,2,1))*UMU(k,1,1)
     . +(dsqrt(g1q)*N(i,1,2)+dsqrt(g2q)*N(i,2,2))*UMU(k,1,2)/dsqrt(2d0)
     .     -Ymu*(N(i,4,1)*UMU(k,2,1)+N(i,4,2)*UMU(k,2,2)))**2
     . +((dsqrt(g1q)*N(i,1,1)+dsqrt(g2q)*N(i,2,1))*UMU(k,1,2)/dsqrt(2d0)
     .  -(dsqrt(g1q)*N(i,1,2)+dsqrt(g2q)*N(i,2,2))*UMU(k,1,1)/dsqrt(2d0)
     .     -Ymu*(N(i,4,1)*UMU(k,2,2)-N(i,4,2)*UMU(k,2,1)))**2
     .    +(dsqrt(2d0*g1q)*(N(i,1,1)*UMU(k,2,1)-N(i,1,2)*UMU(k,2,2))
     .     +Ymu*(N(i,4,1)*UMU(k,1,1)-N(i,4,2)*UMU(k,1,2)))**2
     .    +(dsqrt(2d0*g1q)*(N(i,1,1)*UMU(k,2,2)+N(i,1,2)*UMU(k,2,1))
     .     +Ymu*(N(i,4,2)*UMU(k,1,1)+N(i,4,1)*UMU(k,1,2)))**2)
     .                                       *fn1(MNEU(i)/MSMU2(k))
     .    +dsqrt(MNEU(i))/(3d0*MSMU2(k))*((1d0/dsqrt(2d0)
     .  *(dsqrt(g1q)*N(i,1,1)+dsqrt(g2q)*N(i,2,1))*UMU(k,1,1)
     . +(dsqrt(g1q)*N(i,1,2)+dsqrt(g2q)*N(i,2,2))*UMU(k,1,2)/dsqrt(2d0)
     .     -Ymu*(N(i,4,1)*UMU(k,2,1)+N(i,4,2)*UMU(k,2,2)))
     .      *(dsqrt(2d0*g1q)*(N(i,1,1)*UMU(k,2,1)-N(i,1,2)*UMU(k,2,2))
     .     +Ymu*(N(i,4,1)*UMU(k,1,1)-N(i,4,2)*UMU(k,1,2)))
     . -((dsqrt(g1q)*N(i,1,1)+dsqrt(g2q)*N(i,2,1))*UMU(k,1,2)/dsqrt(2d0)
     .  -(dsqrt(g1q)*N(i,1,2)+dsqrt(g2q)*N(i,2,2))*UMU(k,1,1)/dsqrt(2d0)
     .     -Ymu*(N(i,4,1)*UMU(k,2,2)-N(i,4,2)*UMU(k,2,1)))
     .  *(dsqrt(2d0*g1q)*(N(i,1,1)*UMU(k,2,2)+N(i,1,2)*UMU(k,2,1))
     .     +Ymu*(N(i,4,2)*UMU(k,1,1)+N(i,4,1)*UMU(k,1,2)))
     .                                       )*fn2(MNEU(i)/MSMU2(k))
       ENDDO
      ENDDO
      amuNeutr=MMUON/(16d0*pi**2)*aux

*       - 1-Loop Higgs/Muon Contribution (cf. [5,6]):
      aux=0d0
      DO i=1,5
      aux=aux+XH(i,2)**2*fhs(dsqrt(MH2(i))/MMUON)
     .       -XH(i,4)**2*sinb**2*fhp(dsqrt(MH2(i))/MMUON)
      ENDDO
      aux=aux+sinb**2*fhc(dsqrt(MHC2)/MMUON)
      amuHiggs=(Ymu/(4d0*pi))**2*aux
      
      amu1L=amuChar+amuNeutr+amuHiggs


*   C: 2-loop contributions

*       - Large Logarithms: 1-loop-like 2-loop contributions (cf. [7]):
      QNP=dsqrt(Max(MCH2(2),MSMU2(2)))
      amu1L=amu1L*(1d0-4d0*ALEM0/pi*dlog(QNP/MMUON))

*       - 2-loop bosonic contributions (Leading Log) (cf. [8]):
      aux=0d0
      DO i=1,5
      aux=aux+XH(i,2)*(XH(i,2)*cosb-XH(i,1)*sinb)/MH2(i)
      ENDDO
      aux=(cosb**2-sinb**2)*MZ**2/cosb*aux
      aux=(98d0+9d0*aux+23d0*(1d0-4d0*s2TW)**2)/30d0
      amu2Lbos=5d0/(24d0*dsqrt(2d0)*pi**3)*GF*MMUON**2*ALEMMZ
     .   *(aux*dlog(MMUON**2/MW**2))
      
*       - 2-loop Higgs/Fermion contribution (cf. [9]):
      mtq=MT/(1d0+4d0/(3d0*pi)*asf(MT)+11d0/pi**2*asf(MT)**2)
      mbq=runmb(MB)
      aux=0d0
      DO i=1,5
      aux=aux+4d0/3d0*XH(i,1)*XH(i,2)/sinb*FS(MTQ**2/MH2(i))
     .       +1d0/3d0*XH(i,2)**2/cosb*FS(MBQ**2/MH2(i))
     .       +XH(i,2)**2/cosb*FS(mtau**2/MH2(i))
      aux=aux+XH(i,4)**2*(4d0/3d0*FPS(MTQ**2/MH2(i))
     .      +1d0/3d0*FPS(MBQ**2/MH2(i))*tanb**2
     .      +FPS(mtau**2/MH2(i))*tanb**2)*cosb
      ENDDO
      aux=aux/cosb
      amu2LHf=GF*MMUON**2*ALEM0/(4d0*dsqrt(2d0)*pi**3)*aux
      
*       - 2-loop Sfermion contribution (cf. [10]):
      aux=0d0
      DO i=1,5
      DO k=1,2
      lambdHT(i,k)=gRHSTST(i,k,k)*2d0*MW/dsqrt(g2)
      lambdHB(i,k)=gRHSBSB(i,k,k)*2d0*MW/dsqrt(g2)
      lambdHL(i,k)=gRHSLSL(i,k,k)*2d0*MW/dsqrt(g2)
      
      aux=aux+XH(i,2)*(4d0/3d0*lambdHT(i,k)
     .      *FSF(MST2(k)/MH2(i))/MST2(k)
     .     +1d0/3d0*lambdHB(i,k)
     .      *FSF(MSB2(k)/MH2(i))/MSB2(k)
     .     +lambdHL(i,k)
     .      *FSF(MSL2(k)/MH2(i))/MSL2(k))

      lambdHT(i,k)=gIHSTST(i,k,k)*2d0*MW/dsqrt(g2)
      lambdHB(i,k)=gIHSBSB(i,k,k)*2d0*MW/dsqrt(g2)
      lambdHL(i,k)=gIHSLSL(i,k,k)*2d0*MW/dsqrt(g2)
      
      aux=aux-XH(i,4)*sinb*(4d0/3d0*lambdHT(i,k)
     .      *FSF(MST2(k)/MH2(i))/MST2(k)
     .     +1d0/3d0*lambdHB(i,k)
     .      *FSF(MSB2(k)/MH2(i))/MSB2(k)
     .     +lambdHL(i,k)
     .      *FSF(MSL2(k)/MH2(i))/MSL2(k))
      ENDDO
      ENDDO
      amu2LSF=GF*MMUON**2*ALEMMZ/(4d0*dsqrt(2d0)*pi**3*cosb)*aux

*       - 2-loop Chargino contribution (cf. [11]):
      aux=0d0
      DO k=1,2
      DO i=1,5
      aux=aux+XH(i,2)/cosb*MW/dsqrt(MCH2(k))
     .  *COH0CH(i,k,k,1)*fS(MCH2(k)/MH2(i))
      aux=aux+XH(i,4)*tanb*MW/dsqrt(MCH2(k))
     .  *COH0CH(i,k,k,2)*fPS(MCH2(k)/MH2(i))
      ENDDO
      ENDDO
      amu2LCh=GF*MMUON**2*dsqrt(2d0)*ALEMMZ/(8d0*pi**3)*aux


*  D: Conclusion and errors  (cf. [11])
*      -> Pure SUSY contribution to the discrepancy (G-2)mu:
*       (Here the bosonic SM 2-loop contr. is deduced)
      amu2L=amu1L+amu2Lbos+amu2LHf+amu2LSF+amu2LCh-amubosnlo
*      -> Theoretical errors (added linearly):
      amuerr=2.8d-10+.02d0*dabs(amu1L)+.3d0*dabs(amu2L-amu1L)
*       -> Theoretical Central Value:
      delmagmu=amu2L
*       -> Theoretical Error Bounds:
      amuthmax=amu2L+amuerr
      amuthmin=amu2L-amuerr
      
*       - Comparison with the required SUSY contr.:
      IF(damumin-amuthmax.ge.0d0)THEN
      PROB(37)=1d0-amuthmax/damumin
      ENDIF
      IF(damuMax-amuthmin.le.0d0)THEN
      PROB(37)=1d0-amuthmin/damuMax
      ENDIF

      csqmag=(delmagmu-damusm)**2/(errsm**2+amuerr**2/4.d0)

      RETURN
      END
