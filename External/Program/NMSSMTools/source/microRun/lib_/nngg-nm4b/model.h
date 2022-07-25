*     LanHEP output produced at Mon Apr  6 00:14:31 2015
*     Model named 'NMSSM(../spect.slha)+hgg'

      double precision Sqrt2, pi, degree, hbar_c2,bogus
      parameter (Sqrt2=1.41421356237309504880168872421D0)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      parameter (degree = pi/180D0)
      parameter (hbar_c2 = 3.8937966D8)
      parameter (bogus = -1D123)
      double complex cI
      parameter (cI = (0D0, 1D0))

      double precision Divergence
      common /renorm/ Divergence


      double precision EE, MZ, MW, GG, alfSMZ, MbMb, McMc, Mtp
      double precision Ml, Mq, CW, SW, C2W, LamQCD, rd, mu, Lambda
      double precision Kappa, aLambda, aKappa, muP, At, Ab, Al
      double precision Ml2, Ml3, Mr2, Mr3, Mq2, Mq3, Mu2, Mu3
      double precision Md2, Md3, tB, Mh1, Mh2, Mh3, Mha, Mhb, MHc
      double precision MNE1, MNE2, MNE3, MNE4, MNE5, MC1, MC2
      double precision MSG, MSuL, MSuR, MSdL, MSdR, MScL, MScR
      double precision MSsL, MSsR, MSt1, MSt2, MSb1, MSb2, QSUSY
      double precision MSeL, MSeR, MSmL, MSmR, MSl1, MSl2, MSne
      double precision MSnm, MSnl, Zl11, Zl12, Zl21, Zl22, Zh11
      double precision Zh12, Zh13, Zh21, Zh22, Zh23, Zh31, Zh32
      double precision Zh33, Za11, Za12, Za13, Za21, Za22, Za23
      double precision Zu11, Zu12, Zu21, Zu22, Zv11, Zv12, Zv21
      double precision Zv22, Zt11, Zt12, Zt21, Zt22, Zb11, Zb12
      double precision Zb21, Zb22, la1, la2, la3, la4, la5, la6
      double precision la7, la1s, la2s, la3s, la4s, la5s, la6s
      double precision la7s, la8s, aa1, aa2, aa3, aa4, aa5, aa6
      double precision mB1, mB2, X, dMb, xif, xis, MM3, MSP, sb
      double precision cb, t2b, xvev, Pa12, Pa22, Pa11, Pa21, Td3
      double precision Q, Mt, Mb, Mc, xH2, yH2, dMd, Td2, Au, Ad
      double precision fiuu, fidd, ficc, fiss, Zuu11, Zuu12, Zuu21
      double precision Zuu22, Zdd11, Zdd12, Zdd21, Zdd22, Zcc11
      double precision Zcc12, Zcc21, Zcc22, Zss11, Zss12, Zss21
      double precision Zss22, Mele, Mmu, MbMM, MtMM

      double complex Zn11, Zn12, Zn13, Zn14, Zn15, Zn21, Zn22
      double complex Zn23, Zn24, Zn25, Zn31, Zn32, Zn33, Zn34
      double complex Zn35, Zn41, Zn42, Zn43, Zn44, Zn45, Zn51
      double complex Zn52, Zn53, Zn54, Zn55, NMM11, NMM12, NMM13
      double complex NMM14, NMM15, NMM22, NMM23, NMM24, NMM25
      double complex NMM33, NMM34, NMM35, NMM44, NMM45, NMM55
      double complex MG1I, MG2I

      double precision AAABR(541)
      double complex AAABC(344)

      double precision quuMass, qudMass, lpdMass, neuMass, chaMass
      double precision sluMass, sldMass, sleMass, squMass, sqvMass
      double precision sqdMass, sqeMass, hisMass, hiaMass, MTR001
      double precision hisW, hiaW
      double precision MTR002, MTR003, MTR004, MTR005, MTR006
      double precision MTR007, MTR008, MTR009, MTR010, MTR011
      double precision MTR012, MTR013, MTR014, MTR015, MTR016
      double precision MTR017, MTR018, MTR019, MTR020, MTR021
      double precision MTR022, MTR023, MTR024, MTR025, MTR026
      double precision MTR027, MTR028, MTR029, MTR030, MTR031
      double precision MTR032, MTR033, MTR034, MTR035, MTR036
      double precision MTR037, MTR038, MTR039, MTR040, MTR041
      double precision MTR042, MTR043, MTR044, MTR045, MTR046
      double precision MTR047, MTR048, MTR049, MTR050, MTR051
      double precision MTR052, MTR053, MTR054, MTR061, MTR062
      double precision MTR067, MTR068, MTR069, MTR070, MTR071
      double precision MTR072, MTR077, MTR078, MTR084, MTR085
      double precision MTR086, MTR087, MTR088, MTR089, MTR090
      double precision MTR091, MTR092, MTR093, MTR094, MTR095
      double precision MTR096, MTR097, MTR098, MTR103, MTR104
      double precision MTR105, MTR106, MTR107, MTR108, MTR109
      double precision MTR133, MTR134, MTR141, MTR142, MTR143
      double precision MTR144, MTR145, MTR146, MTR147, MTR148
      double precision MTR149, MTR150, MTR151, MTR152, MTR153
      double precision MTR154, MTR155, MTR156, MTR157, MTR158
      double precision MTR159, MTR160, MTR161, MTR162, MTR163
      double precision MTR164, MTR165, MTR166, MTR167, MTR168
      double precision MTR169, MTR170, MTR171, MTR172, MTR173
      double precision MTR174, MTR175, MTR176, MTR177, MTR178
      double precision MTR179, MTR180, MTR181, MTR182, MTR183
      double precision MTR184, MTR185, MTR186, MTR187, MTR188
      double precision MTR189, MTR190, MTR191, MTR192, MTR193
      double complex MTR055, MTR056, MTR057, MTR058, MTR059, MTR060
      double complex MTR063, MTR064, MTR065, MTR066, MTR073, MTR074
      double complex MTR075, MTR076, MTR079, MTR080, MTR081, MTR082
      double complex MTR083, MTR099, MTR100, MTR101, MTR102, MTR110
      double complex MTR111, MTR112, MTR113, MTR114, MTR115, MTR116
      double complex MTR117, MTR118, MTR119, MTR120, MTR121, MTR122
      double complex MTR123, MTR124, MTR125, MTR126, MTR127, MTR128
      double complex MTR129, MTR130, MTR131, MTR132, MTR135, MTR136
      double complex MTR137, MTR138, MTR139, MTR140

      dimension
     & quuMass(3), qudMass(3), lpdMass(3), neuMass(5), chaMass(2), 
     & sluMass(3), sldMass(3), sleMass(3), squMass(3), sqvMass(3), 
     & sqdMass(3), sqeMass(3), hisMass(3), hiaMass(2), MTR001(3), 
     & hisW(3), hiaW(2),
     & MTR002(3), MTR003(3), MTR004(3), MTR005(3), MTR006(3), 
     & MTR007(3), MTR008(3), MTR009(3), MTR010(2), MTR011(3), 
     & MTR012(3), MTR013(2,3), MTR014(3,3), MTR015(3,3), 
     & MTR016(3,3), MTR017(3,3), MTR018(3), MTR019(3), MTR020(3), 
     & MTR021(2,3), MTR022(3,3), MTR023(3,3), MTR024(3), 
     & MTR025(3), MTR026(3,3), MTR027(3), MTR028(2,3), MTR029(3,3), 
     & MTR030(3,3), MTR031(3,3), MTR032(3), MTR033(2,3), 
     & MTR034(2,2,3), MTR035(3,3,3), MTR036(2), MTR037(3), 
     & MTR038(3), MTR039(3), MTR040(3), MTR041(3), MTR042(3), 
     & MTR043(3), MTR044(3), MTR045(3), MTR046(3), MTR047(3), 
     & MTR048(3), MTR049(3), MTR050(3), MTR051(3), MTR052(3), 
     & MTR053(2,3), MTR054(3), MTR055(2,2), MTR056(2,2), 
     & MTR057(2,2,2), MTR058(2,2,2), MTR059(2,2,3), MTR060(2,2,3), 
     & MTR061(2,3), MTR062(2,3), MTR063(2,5), MTR064(2,5), 
     & MTR065(2,5), MTR066(2,5), MTR067(2,3), MTR068(2,3), 
     & MTR069(2,3), MTR070(2,3), MTR071(3,2), MTR072(3,3), 
     & MTR073(5,3), MTR074(5,3), MTR075(5,3), MTR076(5,3), 
     & MTR077(2,3), MTR078(2,3), MTR079(5,3), MTR080(5,3), 
     & MTR081(5,3), MTR082(5,3), MTR083(5,3), MTR084(3), 
     & MTR085(3,2), MTR086(3,3), MTR087(3), MTR088(3), MTR089(3), 
     & MTR090(3), MTR091(2,3), MTR092(2,3), MTR093(2,3), 
     & MTR094(2,3), MTR095(3), MTR096(3), MTR097(3), MTR098(3), 
     & MTR099(5,3), MTR100(5,3), MTR101(5,3), MTR102(5,3), 
     & MTR103(3), MTR104(3,2), MTR105(3,3), MTR106(3), MTR107(3), 
     & MTR108(3), MTR109(3), MTR110(5,3), MTR111(5,3), MTR112(5,3), 
     & MTR113(5,3), MTR114(5,3), MTR115(5,3), MTR116(5,3), 
     & MTR117(5,3), MTR118(5,3), MTR119(5,3), MTR120(5,3), 
     & MTR121(5,3), MTR122(5,3), MTR123(2,5), MTR124(2,5), 
     & MTR125(2,5), MTR126(2,5), MTR127(5,5), MTR128(5,5), 
     & MTR129(5,5,2), MTR130(5,5,2), MTR131(5,5,3), MTR132(5,5,3), 
     & MTR133(2,2), MTR134(2,2), MTR135(2,5), MTR136(2,5), 
     & MTR137(2,5), MTR138(2,5), MTR139(5,5), MTR140(5,5), 
     & MTR141(2), MTR142(2), MTR143(3), MTR144(3), MTR145(3), 
     & MTR146(3), MTR147(3), MTR148(3), MTR149(3), MTR150(3), 
     & MTR151(3), MTR152(3), MTR153(3), MTR154(3), MTR155(3), 
     & MTR156(3), MTR157(3), MTR158(3), MTR159(3), MTR160(3), 
     & MTR161(3), MTR162(3), MTR163(3), MTR164(3), MTR165(3), 
     & MTR166(3), MTR167(3), MTR168(3), MTR169(3), MTR170(3), 
     & MTR171(3), MTR172(3), MTR173(3), MTR174(3), MTR175(3), 
     & MTR176(3), MTR177(3), MTR178(3), MTR179(3), MTR180(3), 
     & MTR181(3), MTR182(3), MTR183(3), MTR184(3), MTR185(3), 
     & MTR186(3), MTR187(3), MTR188(3), MTR189(3), MTR190(2,2), 
     & MTR191(2,2), MTR192(3,3), MTR193(3,3)

      common /mdl_mtrces/
     & quuMass, qudMass, lpdMass, neuMass, chaMass, sluMass, 
     & sldMass, sleMass, squMass, sqvMass, sqdMass, sqeMass, 
     & hisMass, hiaMass, hisW, hiaW, MTR001, MTR002, MTR003, MTR004, 
     & MTR005, MTR006, MTR007, MTR008, MTR009, MTR010, MTR011, 
     & MTR012, MTR013, MTR014, MTR015, MTR016, MTR017, MTR018, 
     & MTR019, MTR020, MTR021, MTR022, MTR023, MTR024, MTR025, 
     & MTR026, MTR027, MTR028, MTR029, MTR030, MTR031, MTR032, 
     & MTR033, MTR034, MTR035, MTR036, MTR037, MTR038, MTR039, 
     & MTR040, MTR041, MTR042, MTR043, MTR044, MTR045, MTR046, 
     & MTR047, MTR048, MTR049, MTR050, MTR051, MTR052, MTR053, 
     & MTR054, MTR055, MTR056, MTR057, MTR058, MTR059, MTR060, 
     & MTR061, MTR062, MTR063, MTR064, MTR065, MTR066, MTR067, 
     & MTR068, MTR069, MTR070, MTR071, MTR072, MTR073, MTR074, 
     & MTR075, MTR076, MTR077, MTR078, MTR079, MTR080, MTR081, 
     & MTR082, MTR083, MTR084, MTR085, MTR086, MTR087, MTR088, 
     & MTR089, MTR090, MTR091, MTR092, MTR093, MTR094, MTR095, 
     & MTR096, MTR097, MTR098, MTR099, MTR100, MTR101, MTR102, 
     & MTR103, MTR104, MTR105, MTR106, MTR107, MTR108, MTR109, 
     & MTR110, MTR111, MTR112, MTR113, MTR114, MTR115, MTR116, 
     & MTR117, MTR118, MTR119, MTR120, MTR121, MTR122, MTR123, 
     & MTR124, MTR125, MTR126, MTR127, MTR128, MTR129, MTR130, 
     & MTR131, MTR132, MTR133, MTR134, MTR135, MTR136, MTR137, 
     & MTR138, MTR139, MTR140, MTR141, MTR142, MTR143, MTR144, 
     & MTR145, MTR146, MTR147, MTR148, MTR149, MTR150, MTR151, 
     & MTR152, MTR153, MTR154, MTR155, MTR156, MTR157, MTR158, 
     & MTR159, MTR160, MTR161, MTR162, MTR163, MTR164, MTR165, 
     & MTR166, MTR167, MTR168, MTR169, MTR170, MTR171, MTR172, 
     & MTR173, MTR174, MTR175, MTR176, MTR177, MTR178, MTR179, 
     & MTR180, MTR181, MTR182, MTR183, MTR184, MTR185, MTR186, 
     & MTR187, MTR188, MTR189, MTR190, MTR191, MTR192, MTR193

      common /mdl_para/
     &    EE, MZ, MW, GG, alfSMZ, MbMb, McMc, Mtp, Ml, Mq, CW,
     &    SW, C2W, LamQCD, rd, mu, Lambda, Kappa, aLambda, aKappa,
     &    muP, At, Ab, Al, Ml2, Ml3, Mr2, Mr3, Mq2, Mq3, Mu2, Mu3,
     &    Md2, Md3, tB, Mh1, Mh2, Mh3, Mha, Mhb, MHc, MNE1, MNE2,
     &    MNE3, MNE4, MNE5, MC1, MC2, MSG, MSuL, MSuR, MSdL, MSdR,
     &    MScL, MScR, MSsL, MSsR, MSt1, MSt2, MSb1, MSb2, QSUSY,
     &    MSeL, MSeR, MSmL, MSmR, MSl1, MSl2, MSne, MSnm, MSnl,
     &    Zl11, Zl12, Zl21, Zl22, Zh11, Zh12, Zh13, Zh21, Zh22,
     &    Zh23, Zh31, Zh32, Zh33, Za11, Za12, Za13, Za21, Za22,
     &    Za23, Zn11, Zn12, Zn13, Zn14, Zn15, Zn21, Zn22, Zn23,
     &    Zn24, Zn25, Zn31, Zn32, Zn33, Zn34, Zn35, Zn41, Zn42,
     &    Zn43, Zn44, Zn45, Zn51, Zn52, Zn53, Zn54, Zn55, Zu11,
     &    Zu12, Zu21, Zu22, Zv11, Zv12, Zv21, Zv22, Zt11, Zt12,
     &    Zt21, Zt22, Zb11, Zb12, Zb21, Zb22, la1, la2, la3, la4,
     &    la5, la6, la7, la1s, la2s, la3s, la4s, la5s, la6s, la7s,
     &    la8s, aa1, aa2, aa3, aa4, aa5, aa6, mB1, mB2, X, dMb, xif,
     &    xis, MM3, MSP, sb, cb, t2b, xvev, Pa12, Pa22, Pa11, Pa21,
     &    Td3, Q, Mt, Mb, Mc, xH2, yH2, dMd, Td2, Au, Ad, fiuu,
     &    fidd, ficc, fiss, Zuu11, Zuu12, Zuu21, Zuu22, Zdd11,
     &    Zdd12, Zdd21, Zdd22, Zcc11, Zcc12, Zcc21, Zcc22, Zss11,
     &    Zss12, Zss21, Zss22, Mele, Mmu, MbMM, MtMM, NMM11, NMM12,
     &    NMM13, NMM14, NMM15, NMM22, NMM23, NMM24, NMM25, NMM33,
     &    NMM34, NMM35, NMM44, NMM45, NMM55, MG1I, MG2I, AAABR, AAABC

