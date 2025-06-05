Install["bin/gm2calc.mx"];

GM2CalcSetFlags[
    loopOrder -> 2,
    tanBetaResummation -> True,
    forceOutput -> False,
    runningCouplings -> True];

GM2CalcSetSMParameters[
    alpha0 -> 0.00729735,
    alphaMZ -> 0.0077552,
    alphaS -> 0.1184,
    MhSM -> 125.09,
    MW -> 80.385,
    MZ -> 91.1876,
    MT -> 173.34,
    mcmc -> 1.28,
    mu2GeV -> 0.0022,
    mbmb -> 4.18,
    ms2GeV -> 0.096,
    md2GeV -> 0.0047,
    ML -> 1.777,
    MM -> 0.1056583715,
    ME -> 0.000510998928,
    Mv1 -> 0,
    Mv2 -> 0,
    Mv3 -> 0,
    CKM -> IdentityMatrix[3] ];

(* calculate amu using the gauge basis input parameters *)
Print[{amu, Damu} /. GM2CalcAmuTHDMGaugeBasis[
    yukawaType        -> 2,
    lambda            -> { 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1 },
    TB                -> 3,
    m122              -> 200^2,
    zetau             -> 0,
    zetad             -> 0,
    zetal             -> 0,
    Deltau            -> 0 IdentityMatrix[3],
    Deltad            -> 0 IdentityMatrix[3],
    Deltal            -> 0 IdentityMatrix[3],
    Piu               -> 0 IdentityMatrix[3],
    Pid               -> 0 IdentityMatrix[3],
    Pil               -> 0 IdentityMatrix[3]]
];

(* calculate amu using the mass basis input parameters *)
Print[{amu, Damu} /. GM2CalcAmuTHDMMassBasis[
    yukawaType        -> 2,
    Mhh               -> { 125, 400 },
    MAh               -> 420,
    MHp               -> 440,
    sinBetaMinusAlpha -> 0.999,
    lambda6           -> 0,
    lambda7           -> 0,
    TB                -> 3,
    m122              -> 200^2,
    zetau             -> 0,
    zetad             -> 0,
    zetal             -> 0,
    Deltau            -> 0 IdentityMatrix[3],
    Deltad            -> 0 IdentityMatrix[3],
    Deltal            -> 0 IdentityMatrix[3],
    Piu               -> 0 IdentityMatrix[3],
    Pid               -> 0 IdentityMatrix[3],
    Pil               -> 0 IdentityMatrix[3]]
];
