Install["bin/gm2calc.mx"];

GM2CalcSetFlags[
    loopOrder -> 2,
    tanBetaResummation -> True,
    forceOutput -> False];

GM2CalcSetSMParameters[
    alphaMZ -> 0.0077552,   (* 1L *)
    alpha0 -> 0.00729735,   (* 2L *)
    alphaS -> 0.1184,       (* 2L *)
    MW -> 80.385,           (* 1L *)
    MZ -> 91.1876,          (* 1L *)
    MT -> 173.34,           (* 2L *)
    mbmb -> 4.18,           (* 2L *)
    ML -> 1.777,            (* 2L *)
    MM -> 0.1056583715];    (* 1L *)

(* calculate amu using the GM2Calc parameters *)
Print[{amu, Damu} /. GM2CalcAmuGM2CalcScheme[
    MAh    -> 1500,    (* 2L *)
    TB     -> 10,      (* 1L *)
    Mu     -> 350,     (* 1L *)
    MassB  -> 150,     (* 1L *)
    MassWB -> 300,     (* 1L *)
    MassG  -> 1000,    (* 2L *)
    mq2    -> 500^2 IdentityMatrix[3], (* 2L *)
    ml2    -> 500^2 IdentityMatrix[3], (* 1L *)
    mu2    -> 500^2 IdentityMatrix[3], (* 2L *)
    md2    -> 500^2 IdentityMatrix[3], (* 2L *)
    me2    -> 500^2 IdentityMatrix[3], (* 2L *)
    Au     -> 0 IdentityMatrix[3],     (* 2L *)
    Ad     -> 0 IdentityMatrix[3],     (* 2L *)
    Ae     -> 0 IdentityMatrix[3],     (* 1L *)
    Q      -> 454.7]                   (* 2L *)
];
