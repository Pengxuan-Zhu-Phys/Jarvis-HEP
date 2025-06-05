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

(* calculate amu using the SLHA parameters *)
Print[{amu, Damu} /. GM2CalcAmuSLHAScheme[
    (* pole masses *)
    MSvmL  -> 5.18860573 10^2,                    (* 1L *)
    MSm    -> {5.05095249 10^2, 5.25187016 10^2}, (* 1L *)
    MChi   -> {2.01611468 10^2, 4.10040273 10^2, -5.16529941 10^2, 5.45628749 10^2},  (* 1L *)
    MCha   -> {4.0998989 10^2, 5.46057190 10^2},  (* 1L *)
    MAh    -> 1.5 10^3,                           (* 2L *)
    TB     -> 40,                                 (* 1L, DR-bar scheme *)
    Mu     -> 500,                                (* initial guess *)
    MassB  -> 200,                                (* initial guess *)
    MassWB -> 400,                                (* initial guess *)
    MassG  -> 2000,                               (* 2L *)
    mq2    -> 7000^2 IdentityMatrix[3],           (* 2L *)
    ml2    -> 500^2 IdentityMatrix[3],            (* 2L *)
    mu2    -> 7000^2 IdentityMatrix[3],           (* 2L *)
    md2    -> 7000^2 IdentityMatrix[3],           (* 2L *)
    me2    -> 500^2 IdentityMatrix[3],            (* 2L *)
    Au     -> 0 IdentityMatrix[3],                (* 2L *)
    Ad     -> 0 IdentityMatrix[3],                (* 2L *)
    Ae     -> 0 IdentityMatrix[3],                (* 1L, DR-bar scheme *)
    Q      -> 1000]                               (* 2L *)
];
