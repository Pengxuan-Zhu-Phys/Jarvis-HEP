:Evaluate: BeginPackage["GM2Calc`"]

:Evaluate: loopOrder::usage =
    "Loop order of the calculation. Possible values: 0, 1, 2";

:Evaluate: tanBetaResummation::usage =
    "Enable/disable tan(beta) resummation. Possible values: True, False";

:Evaluate: forceOutput::usage =
    "Enforce output, even in case an error has occurred.";

:Evaluate: runningCouplings::usage =
    "Use running couplings in the THDM.";

:Evaluate: alphaMZ::usage =
    "Electromagnetic coupling in the MS-bar scheme at the scale MZ."

:Evaluate: alpha0::usage =
    "Electromagnetic coupling in the Thomson limit."

:Evaluate: alphaS::usage =
    "Strong coupling."

:Evaluate: MW::usage =
    "W-boson pole mass."

:Evaluate: MZ::usage =
    "Z-boson pole mass."

:Evaluate: MT::usage =
    "Top quark pole mass."

:Evaluate: mcmc::usage =
    "Charm quark MS-bar mass mc at the scale mc."

:Evaluate: mu2GeV::usage =
    "Up quark MS-bar mass mu at the scale 2 GeV."

:Evaluate: mbmb::usage =
    "Bottom quark MS-bar mass mb at the scale mb."

:Evaluate: ms2GeV::usage =
    "Strange quark MS-bar mass ms at the scale 2 GeV."

:Evaluate: md2GeV::usage =
    "Down quark MS-bar mass ms at the scale 2 GeV."

:Evaluate: mbMZ::usage =
    "Bottom quark DR-bar mass at the scale MZ."

:Evaluate: ML::usage =
    "Tau lepton pole mass."

:Evaluate: MM::usage =
    "Muon pole mass."

:Evaluate: ME::usage =
    "Electron pole mass."

:Evaluate: Mv1::usage =
    "Lightest neutrino pole mass."

:Evaluate: Mv2::usage =
    "2nd lightest neutrino pole mass."

:Evaluate: Mv3::usage =
    "Heaviest neutrino pole mass."

:Evaluate: CKM::usage =
    "CKM matrix."

:Evaluate: amu::usage =
    "Calculated value of the anomalous magnetic moment of the muon, amu = (g-2)/2.";

:Evaluate: amu1L::usage =
    "Calculated value of the anomalous magnetic moment of the muon, amu = (g-2)/2 (1-loop contribution only).";

:Evaluate: amu2L::usage =
    "Calculated value of the anomalous magnetic moment of the muon, amu = (g-2)/2 (2-loop contribution only).";

:Evaluate: amu2LF::usage =
    "Calculated value of the anomalous magnetic moment of the muon, amu = (g-2)/2 (2-loop fermionic contribution only).";

:Evaluate: amu2LB::usage =
    "Calculated value of the anomalous magnetic moment of the muon, amu = (g-2)/2 (2-loop bosonic contribution only).";

:Evaluate: Damu::usage =
    "Uncertainty of the calculated value of the anomalous magnetic moment of the muon.";

:Evaluate: UM::usage =
    "Mixing matrix of the negatively charged charginos.";

:Evaluate: UP::usage =
    "Mixing matrix of the positively charged charginos.";

:Evaluate: ZN::usage =
    "Mixing matrix of the neutralinos.";

:Evaluate: USt::usage =
    "Mixing matrix of the stops.";

:Evaluate: USb::usage =
    "Mixing matrix of the sbottoms.";

:Evaluate: UStau::usage =
    "Mixing matrix of the staus.";

:Evaluate: USm::usage =
    "Mixing matrix of the smuons.";

:Evaluate: MSveL::usage =
    "Electron-sneutrino mass.";

:Evaluate: MSvmL::usage =
    "Muon-sneutrino mass.";

:Evaluate: MSvtL::usage =
    "Tau-sneutrino mass.";

:Evaluate: MSe::usage =
    "Selectron masses.";

:Evaluate: MSm::usage =
    "Smuon masses.";

:Evaluate: MStau::usage =
    "Stau masses.";

:Evaluate: MSu::usage =
    "Sup masses.";

:Evaluate: MSd::usage =
    "Sdown masses.";

:Evaluate: MSc::usage =
    "Scharm masses.";

:Evaluate: MSs::usage =
    "Sstrange masses.";

:Evaluate: MSt::usage =
    "Stop masses.";

:Evaluate: MSb::usage =
    "Sbottom masses.";

:Evaluate: MChi::usage =
    "Neutralino masses.";

:Evaluate: MCha::usage =
    "Chargino masses.";

:Evaluate: MhSM::usage =
    "Standard Model Higgs boson mass.";

:Evaluate: Mhh::usage =
    "CP-even Higgs boson masses.";

:Evaluate: MAh::usage =
    "CP-odd Higgs boson mass.";

:Evaluate: MHp::usage =
    "charged Higgs boson mass.";

:Evaluate: Mhh::usage =
    "CP-even Higgs bosons masses.";

:Evaluate: TB::usage =
    "tan(beta) = v2/v1 (= vu/vd in the MSSM).";

:Evaluate: Mu::usage =
    "Superpotential mu-parameter.";

:Evaluate: MassB::usage =
    "Soft-breaking bino mass parameter.";

:Evaluate: MassWB::usage =
    "Soft-breaking wino mass parameter.";

:Evaluate: MassG::usage =
    "Soft-breaking gluino mass parameter.";

:Evaluate: mq2::usage =
    "Soft-breaking squared left-handed squark mass parameters.";

:Evaluate: mu2::usage =
    "Soft-breaking squared right-handed up-type squark mass parameters.";

:Evaluate: md2::usage =
    "Soft-breaking squared right-handed down-type squark mass parameters.";

:Evaluate: ml2::usage =
    "Soft-breaking squared left-handed slepton mass parameters.";

:Evaluate: me2::usage =
    "Soft-breaking squared right-handed down-type slepton mass parameters.";

:Evaluate: Au::usage =
    "Soft-breaking trilinear couplings between Higgs bosons and up-type squarks.";

:Evaluate: Ad::usage =
    "Soft-breaking trilinear couplings between Higgs bosons and down-type squarks.";

:Evaluate: Ae::usage =
    "Soft-breaking trilinear couplings between Higgs bosons and down-type sleptons.";

:Evaluate: Q::usage =
    "Renormalization scale.";

:Evaluate: Yu::usage =
    "Yukawa couplings of up-type quarks.";

:Evaluate: Yd::usage =
    "Yukawa couplings of down-type quarks.";

:Evaluate: Ye::usage =
    "Yukawa couplings of down-type leptons.";

:Evaluate: yukawaType::usage =
    "Yukawa type of the Two-Higgs Doublet Model.";

:Evaluate: sinBetaMinusAlpha::usage =
    "Mixing angle sin(beta - alpha) of the Higgs sector in the Two-Higgs Doublet Model.";

:Evaluate: lambda::usage =
    "Lagrangian parameters { lambda_1, ... , lambda_7 } in the Two-Higgs Doublet Model.";

:Evaluate: lambda6::usage =
    "Lagrangian parameter lambda_6 in the Two-Higgs Doublet Model.";

:Evaluate: lambda7::usage =
    "Lagrangian parameter lambda_7 in the Two-Higgs Doublet Model.";

:Evaluate: m122::usage =
    "Lagrangian parameter m_{12}^2 in the Two-Higgs Doublet Model.";

:Evaluate: zetau::usage =
    "Alignment parameter zeta_u in the Two-Higgs Doublet Model.";

:Evaluate: zetad::usage =
    "Alignment parameter zeta_d in the Two-Higgs Doublet Model.";

:Evaluate: zetal::usage =
    "Alignment parameter zeta_l in the Two-Higgs Doublet Model.";

:Evaluate: Deltau::usage =
    "Yukawa coupling Delta_u in the general Two-Higgs Doublet Model.";

:Evaluate: Deltad::usage =
    "Yukawa coupling Delta_d in the general Two-Higgs Doublet Model.";

:Evaluate: Deltal::usage =
    "Yukawa coupling Delta_l in the general Two-Higgs Doublet Model.";

:Evaluate: Piu::usage =
    "Yukawa coupling Pi_u in the general Two-Higgs Doublet Model.";

:Evaluate: Pid::usage =
    "Yukawa coupling Pi_d in the general Two-Higgs Doublet Model.";

:Evaluate: Pil::usage =
    "Yukawa coupling Pi_l in the general Two-Higgs Doublet Model.";

:Evaluate: GM2CalcSetFlags::usage =
    "GM2CalcSetFlags sets the configuration flags for GM2Calc.  Available flags are: {loopOrder, tanBetaResummation, forceOutput, runningCouplings}.  Unset flags are set to their default values, see Options[GM2CalcSetFlags].  Use GM2CalcGetFlags[] to retrieve the flags currently set."

:Evaluate: GM2CalcGetFlags::usage =
    "GM2CalcGetFlags returns the current configuration flags for GM2Calc."

:Evaluate: GM2CalcSetSMParameters::usage =
    "GM2CalcSetSMParameters sets the Standard Model parameters input parameters. Unset parameters are set to their default values, see Options[GM2CalcSetSMParameters].  Use GM2CalcGetSMParameters[] to retrieve the current values of the Standard Model parameters."

:Evaluate: GM2CalcGetSMParameters::usage =
    "GM2CalcGetSMParameters returns the Standard Model parameters."

:Evaluate: GM2CalcAmuSLHAScheme::usage =
    "GM2CalcAmuSLHAScheme calculates amu and its uncertainty in the MSSM using the given SLHA parameters (SLHA interface).  Unset SLHA parameters are set to zero.  See Options[GM2CalcAmuSLHAScheme] for all SLHA parameters and their default values."

:Evaluate: GM2CalcAmuGM2CalcScheme::usage =
    "GM2CalcAmuGM2CalcScheme calculates amu and its uncertainty in the MSSM using the given parameters in the GM2Calc-specific renormalization scheme (GM2Calc interface).  Unset parameters are set to zero.  See Options[GM2CalcAmuGM2CalcScheme] for all input parameters in the GM2Calc scheme and their default values."

:Evaluate: GM2CalcAmuTHDMGaugeBasis::usage =
    "GM2CalcAmuTHDMGaugeBasis calculates amu and its uncertainty in the Two-Higgs Doublet Model using the given input parameters in the gauge bassis.  Unset input parameters are set to zero.  See Options[GM2CalcAmuTHDMGaugeBasis] for all parameters and their default values."

:Evaluate: GM2CalcAmuTHDMMassBasis::usage =
    "GM2CalcAmuTHDMMassBasis calculates amu and its uncertainty in the Two-Higgs Doublet Model using the given input parameters.  Unset input parameters are set to zero.  See Options[GM2CalcAmuTHDMMassBasis] for all parameters and their default values."

:Evaluate: GM2CalcAmuSLHAScheme::error = "`1`";
:Evaluate: GM2CalcAmuSLHAScheme::warning = "`1`";

:Evaluate: GM2CalcAmuGM2CalcScheme::error = "`1`";
:Evaluate: GM2CalcAmuGM2CalcScheme::warning = "`1`";

:Evaluate: GM2CalcAmuTHDMMassBasis::error = "`1`";
:Evaluate: GM2CalcAmuTHDMGaugeBasis::error = "`1`";

:Evaluate: Begin["`Private`"]

:Begin:
:Function: GM2CalcSetFlags
:Pattern: GM2CalcSetFlags[OptionsPattern[]]
:Arguments: {
   OptionValue[loopOrder],
   Boole[OptionValue[tanBetaResummation]],
   Boole[OptionValue[forceOutput]],
   Boole[OptionValue[runningCouplings]] }
:ArgumentTypes: {Integer, Integer, Integer, Integer}
:ReturnType: Integer
:End:

:Evaluate: Options[GM2CalcSetFlags] = {
   loopOrder -> 2,
   tanBetaResummation -> True,
   forceOutput -> False,
   runningCouplings -> False }

:Evaluate: GM2CalcSetFlags::wronglooporder = "Unsupported loop order: `1`";

:Begin:
:Function: GM2CalcGetFlags
:Pattern: GM2CalcGetFlags[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: GM2CalcSetSMParameters
:Pattern: GM2CalcSetSMParameters[OptionsPattern[]]
:Arguments: {
   N @ OptionValue[alpha0],
   N @ OptionValue[alphaMZ],
   N @ OptionValue[alphaS],
   N @ OptionValue[MhSM],
   N @ OptionValue[MW],
   N @ OptionValue[MZ],
   N @ OptionValue[MT],
   N @ OptionValue[mcmc],
   N @ OptionValue[mu2GeV],
   N @ OptionValue[mbmb],
   N @ OptionValue[ms2GeV],
   N @ OptionValue[md2GeV],
   N @ OptionValue[Mv1],
   N @ OptionValue[Mv2],
   N @ OptionValue[Mv3],
   N @ OptionValue[ML],
   N @ OptionValue[MM],
   N @ OptionValue[ME],
   Re @ N @ OptionValue[CKM][[1,1]],
   Re @ N @ OptionValue[CKM][[1,2]],
   Re @ N @ OptionValue[CKM][[1,3]],
   Re @ N @ OptionValue[CKM][[2,1]],
   Re @ N @ OptionValue[CKM][[2,2]],
   Re @ N @ OptionValue[CKM][[2,3]],
   Re @ N @ OptionValue[CKM][[3,1]],
   Re @ N @ OptionValue[CKM][[3,2]],
   Re @ N @ OptionValue[CKM][[3,3]],
   Im @ N @ OptionValue[CKM][[1,1]],
   Im @ N @ OptionValue[CKM][[1,2]],
   Im @ N @ OptionValue[CKM][[1,3]],
   Im @ N @ OptionValue[CKM][[2,1]],
   Im @ N @ OptionValue[CKM][[2,2]],
   Im @ N @ OptionValue[CKM][[2,3]],
   Im @ N @ OptionValue[CKM][[3,1]],
   Im @ N @ OptionValue[CKM][[3,2]],
   Im @ N @ OptionValue[CKM][[3,3]] }
:ArgumentTypes: {
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real }
:ReturnType: Integer
:End:

:Evaluate: Options[GM2CalcSetSMParameters] = {
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
    CKM -> {{1,0,0}, {0,1,0}, {0,0,1}} }

:Begin:
:Function: GM2CalcGetSMParameters
:Pattern: GM2CalcGetSMParameters[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: GM2CalcAmuSLHAScheme
:Pattern: GM2CalcAmuSLHAScheme[OptionsPattern[]]
:Arguments: {
    N @ OptionValue[MSvmL],
    N @ OptionValue[MSm][[1]],
    N @ OptionValue[MSm][[2]],
    N @ OptionValue[MChi][[1]],
    N @ OptionValue[MChi][[2]],
    N @ OptionValue[MChi][[3]],
    N @ OptionValue[MChi][[4]],
    N @ OptionValue[MCha][[1]],
    N @ OptionValue[MCha][[2]],
    N @ OptionValue[MAh],
    N @ OptionValue[TB],
    N @ OptionValue[Mu],
    N @ OptionValue[MassB],
    N @ OptionValue[MassWB],
    N @ OptionValue[MassG],
    N @ OptionValue[mq2][[1,1]],
    N @ OptionValue[mq2][[2,2]],
    N @ OptionValue[mq2][[3,3]],
    N @ OptionValue[ml2][[1,1]],
    N @ OptionValue[ml2][[2,2]],
    N @ OptionValue[ml2][[3,3]],
    N @ OptionValue[mu2][[1,1]],
    N @ OptionValue[mu2][[2,2]],
    N @ OptionValue[mu2][[3,3]],
    N @ OptionValue[md2][[1,1]],
    N @ OptionValue[md2][[2,2]],
    N @ OptionValue[md2][[3,3]],
    N @ OptionValue[me2][[1,1]],
    N @ OptionValue[me2][[2,2]],
    N @ OptionValue[me2][[3,3]],
    N @ OptionValue[Au][[3,3]],
    N @ OptionValue[Ad][[3,3]],
    N @ OptionValue[Ae][[2,2]],
    N @ OptionValue[Ae][[3,3]],
    N @ OptionValue[Q] }
:ArgumentTypes: {
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real }
:ReturnType: Manual
:End:

:Evaluate: Options[GM2CalcAmuSLHAScheme] = {
    MSvmL     -> 0,
    MSm       -> {0, 0},
    MChi      -> {0, 0, 0, 0},
    MCha      -> {0, 0},
    MAh       -> 0,
    TB        -> 0,
    Mu        -> 0,
    MassB     -> 0,
    MassWB    -> 0,
    MassG     -> 0,
    mq2       -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    ml2       -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    mu2       -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    md2       -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    me2       -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Au        -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Ad        -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Ae        -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Q         -> 0 }

:Begin:
:Function: GM2CalcAmuGM2CalcScheme
:Pattern: GM2CalcAmuGM2CalcScheme[OptionsPattern[]]
:Arguments: {
   N @ OptionValue[MAh],
   N @ OptionValue[TB],
   N @ OptionValue[Mu],
   N @ OptionValue[MassB],
   N @ OptionValue[MassWB],
   N @ OptionValue[MassG],
   N @ OptionValue[mq2][[1,1]],
   N @ OptionValue[mq2][[2,2]],
   N @ OptionValue[mq2][[3,3]],
   N @ OptionValue[ml2][[1,1]],
   N @ OptionValue[ml2][[2,2]],
   N @ OptionValue[ml2][[3,3]],
   N @ OptionValue[mu2][[1,1]],
   N @ OptionValue[mu2][[2,2]],
   N @ OptionValue[mu2][[3,3]],
   N @ OptionValue[md2][[1,1]],
   N @ OptionValue[md2][[2,2]],
   N @ OptionValue[md2][[3,3]],
   N @ OptionValue[me2][[1,1]],
   N @ OptionValue[me2][[2,2]],
   N @ OptionValue[me2][[3,3]],
   N @ OptionValue[Au][[3,3]],
   N @ OptionValue[Ad][[3,3]],
   N @ OptionValue[Ae][[2,2]],
   N @ OptionValue[Ae][[3,3]],
   N @ OptionValue[Q] }
:ArgumentTypes: {
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real }
:ReturnType: Manual
:End:

:Evaluate: Options[GM2CalcAmuGM2CalcScheme] = {
    MAh    -> 0,
    TB     -> 0,
    Mu     -> 0,
    MassB  -> 0,
    MassWB -> 0,
    MassG  -> 0,
    mq2    -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    ml2    -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    mu2    -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    md2    -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    me2    -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Au     -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Ad     -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Ae     -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Q      -> 0 }

:Begin:
:Function: GM2CalcAmuTHDMGaugeBasis
:Pattern: GM2CalcAmuTHDMGaugeBasis[OptionsPattern[]]
:Arguments: {
   N @ OptionValue[yukawaType],
   N @ OptionValue[lambda][[1]],
   N @ OptionValue[lambda][[2]],
   N @ OptionValue[lambda][[3]],
   N @ OptionValue[lambda][[4]],
   N @ OptionValue[lambda][[5]],
   N @ OptionValue[lambda][[6]],
   N @ OptionValue[lambda][[7]],
   N @ OptionValue[TB],
   N @ OptionValue[m122],
   N @ OptionValue[zetau],
   N @ OptionValue[zetad],
   N @ OptionValue[zetal],
   Re @ N @ OptionValue[Deltau][[1,1]],
   Re @ N @ OptionValue[Deltau][[1,2]],
   Re @ N @ OptionValue[Deltau][[1,3]],
   Re @ N @ OptionValue[Deltau][[2,1]],
   Re @ N @ OptionValue[Deltau][[2,2]],
   Re @ N @ OptionValue[Deltau][[2,3]],
   Re @ N @ OptionValue[Deltau][[3,1]],
   Re @ N @ OptionValue[Deltau][[3,2]],
   Re @ N @ OptionValue[Deltau][[3,3]],
   Re @ N @ OptionValue[Deltad][[1,1]],
   Re @ N @ OptionValue[Deltad][[1,2]],
   Re @ N @ OptionValue[Deltad][[1,3]],
   Re @ N @ OptionValue[Deltad][[2,1]],
   Re @ N @ OptionValue[Deltad][[2,2]],
   Re @ N @ OptionValue[Deltad][[2,3]],
   Re @ N @ OptionValue[Deltad][[3,1]],
   Re @ N @ OptionValue[Deltad][[3,2]],
   Re @ N @ OptionValue[Deltad][[3,3]],
   Re @ N @ OptionValue[Deltal][[1,1]],
   Re @ N @ OptionValue[Deltal][[1,2]],
   Re @ N @ OptionValue[Deltal][[1,3]],
   Re @ N @ OptionValue[Deltal][[2,1]],
   Re @ N @ OptionValue[Deltal][[2,2]],
   Re @ N @ OptionValue[Deltal][[2,3]],
   Re @ N @ OptionValue[Deltal][[3,1]],
   Re @ N @ OptionValue[Deltal][[3,2]],
   Re @ N @ OptionValue[Deltal][[3,3]],
   Re @ N @ OptionValue[Piu][[1,1]],
   Re @ N @ OptionValue[Piu][[1,2]],
   Re @ N @ OptionValue[Piu][[1,3]],
   Re @ N @ OptionValue[Piu][[2,1]],
   Re @ N @ OptionValue[Piu][[2,2]],
   Re @ N @ OptionValue[Piu][[2,3]],
   Re @ N @ OptionValue[Piu][[3,1]],
   Re @ N @ OptionValue[Piu][[3,2]],
   Re @ N @ OptionValue[Piu][[3,3]],
   Re @ N @ OptionValue[Pid][[1,1]],
   Re @ N @ OptionValue[Pid][[1,2]],
   Re @ N @ OptionValue[Pid][[1,3]],
   Re @ N @ OptionValue[Pid][[2,1]],
   Re @ N @ OptionValue[Pid][[2,2]],
   Re @ N @ OptionValue[Pid][[2,3]],
   Re @ N @ OptionValue[Pid][[3,1]],
   Re @ N @ OptionValue[Pid][[3,2]],
   Re @ N @ OptionValue[Pid][[3,3]],
   Re @ N @ OptionValue[Pil][[1,1]],
   Re @ N @ OptionValue[Pil][[1,2]],
   Re @ N @ OptionValue[Pil][[1,3]],
   Re @ N @ OptionValue[Pil][[2,1]],
   Re @ N @ OptionValue[Pil][[2,2]],
   Re @ N @ OptionValue[Pil][[2,3]],
   Re @ N @ OptionValue[Pil][[3,1]],
   Re @ N @ OptionValue[Pil][[3,2]],
   Re @ N @ OptionValue[Pil][[3,3]] }
:ArgumentTypes: {
   Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real }
:ReturnType: Manual
:End:

:Evaluate: Options[GM2CalcAmuTHDMGaugeBasis] = {
    yukawaType -> 0,
    lambda     -> { 0, 0, 0, 0, 0, 0, 0 },
    TB         -> 0,
    m122       -> 0,
    zetau      -> 0,
    zetad      -> 0,
    zetal      -> 0,
    Deltau     -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Deltad     -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Deltal     -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Piu        -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Pid        -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Pil        -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}} }

:Begin:
:Function: GM2CalcAmuTHDMMassBasis
:Pattern: GM2CalcAmuTHDMMassBasis[OptionsPattern[]]
:Arguments: {
   N @ OptionValue[yukawaType],
   N @ OptionValue[Mhh][[1]],
   N @ OptionValue[Mhh][[2]],
   N @ OptionValue[MAh],
   N @ OptionValue[MHp],
   N @ OptionValue[sinBetaMinusAlpha],
   N @ OptionValue[lambda6],
   N @ OptionValue[lambda7],
   N @ OptionValue[TB],
   N @ OptionValue[m122],
   N @ OptionValue[zetau],
   N @ OptionValue[zetad],
   N @ OptionValue[zetal],
   Re @ N @ OptionValue[Deltau][[1,1]],
   Re @ N @ OptionValue[Deltau][[1,2]],
   Re @ N @ OptionValue[Deltau][[1,3]],
   Re @ N @ OptionValue[Deltau][[2,1]],
   Re @ N @ OptionValue[Deltau][[2,2]],
   Re @ N @ OptionValue[Deltau][[2,3]],
   Re @ N @ OptionValue[Deltau][[3,1]],
   Re @ N @ OptionValue[Deltau][[3,2]],
   Re @ N @ OptionValue[Deltau][[3,3]],
   Re @ N @ OptionValue[Deltad][[1,1]],
   Re @ N @ OptionValue[Deltad][[1,2]],
   Re @ N @ OptionValue[Deltad][[1,3]],
   Re @ N @ OptionValue[Deltad][[2,1]],
   Re @ N @ OptionValue[Deltad][[2,2]],
   Re @ N @ OptionValue[Deltad][[2,3]],
   Re @ N @ OptionValue[Deltad][[3,1]],
   Re @ N @ OptionValue[Deltad][[3,2]],
   Re @ N @ OptionValue[Deltad][[3,3]],
   Re @ N @ OptionValue[Deltal][[1,1]],
   Re @ N @ OptionValue[Deltal][[1,2]],
   Re @ N @ OptionValue[Deltal][[1,3]],
   Re @ N @ OptionValue[Deltal][[2,1]],
   Re @ N @ OptionValue[Deltal][[2,2]],
   Re @ N @ OptionValue[Deltal][[2,3]],
   Re @ N @ OptionValue[Deltal][[3,1]],
   Re @ N @ OptionValue[Deltal][[3,2]],
   Re @ N @ OptionValue[Deltal][[3,3]],
   Re @ N @ OptionValue[Piu][[1,1]],
   Re @ N @ OptionValue[Piu][[1,2]],
   Re @ N @ OptionValue[Piu][[1,3]],
   Re @ N @ OptionValue[Piu][[2,1]],
   Re @ N @ OptionValue[Piu][[2,2]],
   Re @ N @ OptionValue[Piu][[2,3]],
   Re @ N @ OptionValue[Piu][[3,1]],
   Re @ N @ OptionValue[Piu][[3,2]],
   Re @ N @ OptionValue[Piu][[3,3]],
   Re @ N @ OptionValue[Pid][[1,1]],
   Re @ N @ OptionValue[Pid][[1,2]],
   Re @ N @ OptionValue[Pid][[1,3]],
   Re @ N @ OptionValue[Pid][[2,1]],
   Re @ N @ OptionValue[Pid][[2,2]],
   Re @ N @ OptionValue[Pid][[2,3]],
   Re @ N @ OptionValue[Pid][[3,1]],
   Re @ N @ OptionValue[Pid][[3,2]],
   Re @ N @ OptionValue[Pid][[3,3]],
   Re @ N @ OptionValue[Pil][[1,1]],
   Re @ N @ OptionValue[Pil][[1,2]],
   Re @ N @ OptionValue[Pil][[1,3]],
   Re @ N @ OptionValue[Pil][[2,1]],
   Re @ N @ OptionValue[Pil][[2,2]],
   Re @ N @ OptionValue[Pil][[2,3]],
   Re @ N @ OptionValue[Pil][[3,1]],
   Re @ N @ OptionValue[Pil][[3,2]],
   Re @ N @ OptionValue[Pil][[3,3]] }
:ArgumentTypes: {
   Integer, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
   Real, Real, Real, Real, Real, Real, Real }
:ReturnType: Manual
:End:

:Evaluate: Options[GM2CalcAmuTHDMMassBasis] = {
    yukawaType        -> 0,
    Mhh               -> { 0, 0 },
    MAh               -> 0,
    MHp               -> 0,
    sinBetaMinusAlpha -> 0,
    lambda6           -> 0,
    lambda7           -> 0,
    TB                -> 0,
    m122              -> 0,
    zetau             -> 0,
    zetad             -> 0,
    zetal             -> 0,
    Deltau            -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Deltad            -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Deltal            -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Piu               -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Pid               -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    Pil               -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}} }

:Evaluate: End[]

:Evaluate: EndPackage[]


#include "gm2calc/gm2_1loop.h"
#include "gm2calc/gm2_2loop.h"
#include "gm2calc/gm2_uncertainty.h"
#include "gm2calc/gm2_version.h"
#include "gm2calc/MSSMNoFV_onshell.h"
#include "gm2calc/SM.h"
#include "gm2calc/THDM.h"
#include "gm2_uncertainty_helpers.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mathlink.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NELEMS(x) (sizeof(x) / sizeof((x)[0]))

#define MLPutRule(link,s)                       \
 do {                                           \
    MLPutFunction(link, "Rule", 2);             \
    MLPutSymbol(link, s);                       \
 } while (0)

#define MLPutRuleToReal(link,v,s)               \
 do {                                           \
    MLPutRule(link, s);                         \
    MLPutReal(link, v);                         \
 } while (0)

#define MLPutRuleToInteger(link,v,s)            \
 do {                                           \
    MLPutRule(link, s);                         \
    MLPutInteger(link, v);                      \
 } while (0)

#define MLPutRuleToString(link,v,s)             \
 do {                                           \
    MLPutRule(link, s);                         \
    MLPutString(link, v);                       \
 } while (0)

#define MLPutRealMatrix(link,v,dim1,dim2)                  \
  do {                                                     \
    long dims[] = { dim1, dim2 };                          \
    MLPutDoubleArray(link, v, dims, NULL, NELEMS(dims));   \
  } while (0)

/* macros which access matrix/vector elements through GM2Calc interface
   functions */

#define MLPutRealVectorInterface(link,M,dim)            \
   do {                                                 \
      double M[dim];                                    \
      for (unsigned i = 0; i < dim; i++) {              \
         M[i] = gm2calc_mssmnofv_get_##M(model, i);     \
      }                                                 \
      MLPutRealList(link, M, dim);                      \
   } while (0)

#define MLPutRealMatrixInterface(link,M,dim1,dim2)           \
   do {                                                      \
      double M[dim1][dim2];                                  \
      for (unsigned i = 0; i < dim1; i++) {                  \
         for (unsigned k = 0; k < dim2; k++) {               \
            M[i][k] = gm2calc_mssmnofv_get_##M(model, i, k); \
         }                                                   \
      }                                                      \
      MLPutRealMatrix(link, (double*)M, dim1, dim2);         \
   } while (0)

#define MLPutComplexMatrixInterface(link,M,dim1,dim2)             \
   do {                                                           \
      MLPutFunction(link, "List", dim1);                          \
      for (unsigned i = 0; i < dim1; i++) {                       \
         MLPutFunction(link, "List", dim2);                       \
         for (unsigned k = 0; k < dim2; k++) {                    \
            double re = 0., im = 0.;                              \
            re = gm2calc_mssmnofv_get_##M(model, i, k, &im);      \
            MLPutComplex(link, re, im);                           \
         }                                                        \
      }                                                           \
   } while (0)

#define MLPutComplexMatrixTHDM(link,M,dim1,dim2)                  \
   do {                                                           \
      MLPutFunction(link, "List", dim1);                          \
      for (unsigned i = 0; i < dim1; i++) {                       \
         MLPutFunction(link, "List", dim2);                       \
         for (unsigned k = 0; k < dim2; k++) {                    \
            double re = M##_real[i][k];                           \
            double im = M##_imag[i][k];                           \
            MLPutComplex(link, re, im);                           \
         }                                                        \
      }                                                           \
   } while (0)

#define MLPutRuleToRealVectorInterface(link,M,name,dim) \
 do {                                                   \
    MLPutFunction(link, "Rule", 2);                     \
    MLPutSymbol(link, name);                            \
    MLPutRealVectorInterface(link,M,dim);               \
 } while (0)

#define MLPutRuleToRealMatrixInterface(link,v,name,dim1,dim2)   \
 do {                                                           \
    MLPutFunction(link, "Rule", 2);                             \
    MLPutSymbol(link, name);                                    \
    MLPutRealMatrixInterface(link,v,dim1,dim2);                 \
 } while (0)

#define MLPutRuleToComplexMatrixInterface(link,M,name,dim1,dim2)        \
 do {                                                                   \
    MLPutFunction(link, "Rule", 2);                                     \
    MLPutSymbol(link, name);                                            \
    MLPutComplexMatrixInterface(link,M,dim1,dim2);                      \
 } while (0)

#define MLPutRuleToComplexMatrixTHDM(link,M,name,dim1,dim2)             \
 do {                                                                   \
    MLPutFunction(link, "Rule", 2);                                     \
    MLPutSymbol(link, name);                                            \
    MLPutComplexMatrixTHDM(link,M,dim1,dim2);                           \
 } while (0)

/* global configuration flags */
struct Config_flags {
   int loopOrder;
   int tanBetaResummation;
   int forceOutput;
   int runningCouplings;
} config_flags = {
   .loopOrder = 2,
   .tanBetaResummation = 1,
   .forceOutput = 0,
   .runningCouplings = 1
};

/* Standard Model parameters */
gm2calc_SM sm;

/* SUSY parameters for SLHA interface */
struct SLHA_parameters {
   double MSvmL;
   double MSm[2];
   double MChi[4];
   double MCha[2];
   double MAh;
   double TB;
   double Mu;
   double MassB;
   double MassWB;
   double MassG;
   double mq2[3][3];
   double ml2[3][3];
   double mu2[3][3];
   double md2[3][3];
   double me2[3][3];
   double Au[3][3];
   double Ad[3][3];
   double Ae[3][3];
   double Q;
};

/* SUSY parameters for GM2Calc interface */
struct GM2Calc_parameters {
   double MAh;
   double TB;
   double Mu;
   double MassB;
   double MassWB;
   double MassG;
   double mq2[3][3];
   double ml2[3][3];
   double mu2[3][3];
   double md2[3][3];
   double me2[3][3];
   double Au[3][3];
   double Ad[3][3];
   double Ae[3][3];
   double Q;
};

/******************************************************************/

void initialize_slha_parameters(struct SLHA_parameters* pars)
{
   pars->MSvmL     = 0.;
   memset(pars->MSm, 0, sizeof(pars->MSm[0]) * 2);
   memset(pars->MChi, 0, sizeof(pars->MChi[0]) * 4);
   memset(pars->MCha, 0, sizeof(pars->MCha[0]) * 2);
   pars->MAh       = 0.;
   pars->TB        = 0.;
   pars->Mu        = 0.;
   pars->MassB     = 0.;
   pars->MassWB    = 0.;
   pars->MassG     = 0.;
   memset(pars->mq2, 0, sizeof(pars->mq2[0][0]) * 3 * 3);
   memset(pars->ml2, 0, sizeof(pars->ml2[0][0]) * 3 * 3);
   memset(pars->mu2, 0, sizeof(pars->mu2[0][0]) * 3 * 3);
   memset(pars->md2, 0, sizeof(pars->md2[0][0]) * 3 * 3);
   memset(pars->me2, 0, sizeof(pars->me2[0][0]) * 3 * 3);
   memset(pars->Au, 0, sizeof(pars->Au[0][0]) * 3 * 3);
   memset(pars->Ad, 0, sizeof(pars->Ad[0][0]) * 3 * 3);
   memset(pars->Ae, 0, sizeof(pars->Ae[0][0]) * 3 * 3);
   pars->Q         = 0.;
}

/******************************************************************/

void initialize_gm2calc_parameters(struct GM2Calc_parameters* pars)
{
   pars->MAh     = 0.;
   pars->TB      = 0.;
   pars->Mu      = 0.;
   pars->MassB   = 0.;
   pars->MassWB  = 0.;
   pars->MassG   = 0.;
   memset(pars->mq2, 0, sizeof(pars->mq2[0][0]) * 3 * 3);
   memset(pars->ml2, 0, sizeof(pars->ml2[0][0]) * 3 * 3);
   memset(pars->mu2, 0, sizeof(pars->mu2[0][0]) * 3 * 3);
   memset(pars->md2, 0, sizeof(pars->md2[0][0]) * 3 * 3);
   memset(pars->me2, 0, sizeof(pars->me2[0][0]) * 3 * 3);
   memset(pars->Au, 0, sizeof(pars->Au[0][0]) * 3 * 3);
   memset(pars->Ad, 0, sizeof(pars->Ad[0][0]) * 3 * 3);
   memset(pars->Ae, 0, sizeof(pars->Ae[0][0]) * 3 * 3);
   pars->Q       = 0.;
}

/******************************************************************/

double sqr(double x) { return x*x; }

/******************************************************************/

void MLPutComplex(MLINK link, double re, double im)
{
   if (im == 0.) {
      MLPutReal(link, re);
   } else {
      MLPutFunction(link, "Complex", 2);
      MLPutReal(link, re);
      MLPutReal(link, im);
   }
}

/******************************************************************/

void put_error_message(const char* function_name,
                       const char* message_tag,
                       const char* message_str)
{
   MLPutFunction(stdlink, "CompoundExpression", 2);
   MLPutFunction(stdlink, "Message", 2);
   MLPutFunction(stdlink, "MessageName", 2);
   MLPutSymbol(stdlink, function_name);
   MLPutString(stdlink, message_tag);
   MLPutString(stdlink, message_str);
}

/******************************************************************/

gm2calc_error setup_model_slha_scheme(MSSMNoFV_onshell* model,
                                      const struct SLHA_parameters* pars)
{
   /* fill SM parameters */
   gm2calc_mssmnofv_set_alpha_MZ(model, sm.alpha_em_mz);
   gm2calc_mssmnofv_set_alpha_thompson(model, sm.alpha_em_0);
   gm2calc_mssmnofv_set_g3(model, sqrt(4 * M_PI * sm.alpha_s_mz));
   gm2calc_mssmnofv_set_MT_pole(model, sm.mu[2]);
   gm2calc_mssmnofv_set_MB_running(model, sm.md[2]);
   gm2calc_mssmnofv_set_MM_pole(model, sm.ml[1]);
   gm2calc_mssmnofv_set_ML_pole(model, sm.ml[2]);
   gm2calc_mssmnofv_set_MW_pole(model, sm.mw);
   gm2calc_mssmnofv_set_MZ_pole(model, sm.mz);

   /* fill pole masses */
   gm2calc_mssmnofv_set_MSvmL_pole(model, pars->MSvmL);
   gm2calc_mssmnofv_set_MSm_pole(model, 0, pars->MSm[0]);
   gm2calc_mssmnofv_set_MSm_pole(model, 1, pars->MSm[1]);
   gm2calc_mssmnofv_set_MChi_pole(model, 0, pars->MChi[0]);
   gm2calc_mssmnofv_set_MChi_pole(model, 1, pars->MChi[1]);
   gm2calc_mssmnofv_set_MChi_pole(model, 2, pars->MChi[2]);
   gm2calc_mssmnofv_set_MChi_pole(model, 3, pars->MChi[3]);
   gm2calc_mssmnofv_set_MCha_pole(model, 0, pars->MCha[0]);
   gm2calc_mssmnofv_set_MCha_pole(model, 1, pars->MCha[1]);
   gm2calc_mssmnofv_set_MAh_pole(model, pars->MAh);

   /* fill DR-bar parameters */
   gm2calc_mssmnofv_set_TB(model, pars->TB);
   gm2calc_mssmnofv_set_Mu(model, pars->Mu);
   gm2calc_mssmnofv_set_MassB(model, pars->MassB);
   gm2calc_mssmnofv_set_MassWB(model, pars->MassWB);
   gm2calc_mssmnofv_set_MassG(model, pars->MassG);
   for (unsigned i = 0; i < 3; i++) {
      for (unsigned k = 0; k < 3; k++) {
         gm2calc_mssmnofv_set_mq2(model, i, k, pars->mq2[i][k]);
         gm2calc_mssmnofv_set_ml2(model, i, k, pars->ml2[i][k]);
         gm2calc_mssmnofv_set_mu2(model, i, k, pars->mu2[i][k]);
         gm2calc_mssmnofv_set_md2(model, i, k, pars->md2[i][k]);
         gm2calc_mssmnofv_set_me2(model, i, k, pars->me2[i][k]);
         gm2calc_mssmnofv_set_Au(model, i, k, pars->Au[i][k]);
         gm2calc_mssmnofv_set_Ad(model, i, k, pars->Ad[i][k]);
         gm2calc_mssmnofv_set_Ae(model, i, k, pars->Ae[i][k]);
      }
   }
   gm2calc_mssmnofv_set_scale(model, pars->Q);

   /* convert DR-bar parameters to on-shell */
   const gm2calc_error error = gm2calc_mssmnofv_convert_to_onshell(model);

   return error;
}

/******************************************************************/

gm2calc_error setup_model_gm2calc_scheme(MSSMNoFV_onshell* model,
                                         const struct GM2Calc_parameters* pars)
{
   /* fill SM parameters */
   gm2calc_mssmnofv_set_alpha_MZ(model, sm.alpha_em_mz);
   gm2calc_mssmnofv_set_alpha_thompson(model, sm.alpha_em_0);
   gm2calc_mssmnofv_set_g3(model, sqrt(4 * M_PI * sm.alpha_s_mz));
   gm2calc_mssmnofv_set_MT_pole(model, sm.mu[2]);
   gm2calc_mssmnofv_set_MB_running(model, sm.md[2]);
   gm2calc_mssmnofv_set_MM_pole(model, sm.ml[1]);
   gm2calc_mssmnofv_set_ML_pole(model, sm.ml[2]);
   gm2calc_mssmnofv_set_MW_pole(model, sm.mw);
   gm2calc_mssmnofv_set_MZ_pole(model, sm.mz);

   /* fill DR-bar parameters */
   gm2calc_mssmnofv_set_TB(model, pars->TB);

   /* fill on-shell parameters */
   gm2calc_mssmnofv_set_Mu(model, pars->Mu);
   gm2calc_mssmnofv_set_MassB(model, pars->MassB);
   gm2calc_mssmnofv_set_MassWB(model, pars->MassWB);
   gm2calc_mssmnofv_set_MassG(model, pars->MassG);
   gm2calc_mssmnofv_set_MAh_pole(model, pars->MAh);
   gm2calc_mssmnofv_set_scale(model, pars->Q);

   /* fill DR-bar parameters */
   for (unsigned i = 0; i < 3; i++) {
      for (unsigned k = 0; k < 3; k++) {
         gm2calc_mssmnofv_set_mq2(model, i, k, pars->mq2[i][k]);
         gm2calc_mssmnofv_set_ml2(model, i, k, pars->ml2[i][k]);
         gm2calc_mssmnofv_set_mu2(model, i, k, pars->mu2[i][k]);
         gm2calc_mssmnofv_set_md2(model, i, k, pars->md2[i][k]);
         gm2calc_mssmnofv_set_me2(model, i, k, pars->me2[i][k]);
         gm2calc_mssmnofv_set_Au(model, i, k, pars->Au[i][k]);
         gm2calc_mssmnofv_set_Ad(model, i, k, pars->Ad[i][k]);
         gm2calc_mssmnofv_set_Ae(model, i, k, pars->Ae[i][k]);
      }
   }

   /* convert DR-bar parameters to on-shell */
   const gm2calc_error error = gm2calc_mssmnofv_calculate_masses(model);

   return error;
}

/******************************************************************/

void calculate_amu_mssmnofv(MSSMNoFV_onshell* model, double* amu1L, double* amu2L)
{
   if (config_flags.loopOrder > 0) {
      if (config_flags.tanBetaResummation) {
         *amu1L = gm2calc_mssmnofv_calculate_amu_1loop(model);
      } else {
         *amu1L = gm2calc_mssmnofv_calculate_amu_1loop_non_tan_beta_resummed(model);
      }
   }
   if (config_flags.loopOrder > 1) {
      if (config_flags.tanBetaResummation) {
         *amu2L = gm2calc_mssmnofv_calculate_amu_2loop(model);
      } else {
         *amu2L = gm2calc_mssmnofv_calculate_amu_2loop_non_tan_beta_resummed(model);
      }
   }
}

/******************************************************************/

double calculate_uncertainty_mssmnofv(MSSMNoFV_onshell* model, double amu1L, double amu2L)
{
   double delta_amu = 0.;

   if (config_flags.loopOrder == 0) {
      delta_amu = gm2calc_mssmnofv_calculate_uncertainty_amu_0loop_amu1L(model, amu1L);
   } else if (config_flags.loopOrder == 1) {
      delta_amu = gm2calc_mssmnofv_calculate_uncertainty_amu_1loop_amu2L(model, amu2L);
   } else if (config_flags.loopOrder > 1) {
      delta_amu = gm2calc_mssmnofv_calculate_uncertainty_amu_2loop(model);
   }

   return delta_amu;
}

/******************************************************************/

void calculate_amu_thdm(gm2calc_THDM* model, double* amu1L, double* amu2LF, double* amu2LB)
{
   if (config_flags.loopOrder > 0) {
      *amu1L = gm2calc_thdm_calculate_amu_1loop(model);
   }
   if (config_flags.loopOrder > 1) {
      *amu2LF = gm2calc_thdm_calculate_amu_2loop_fermionic(model);
      *amu2LB = gm2calc_thdm_calculate_amu_2loop_bosonic(model);
   }
}

/******************************************************************/

double calculate_uncertainty_thdm(gm2calc_THDM* model, double amu1L, double amu2L)
{
   double delta_amu = 0.;

   if (config_flags.loopOrder == 0) {
      delta_amu = gm2calc_thdm_calculate_uncertainty_amu_0loop_amu1L_amu2L(model, amu1L, amu2L);
   } else if (config_flags.loopOrder == 1) {
      delta_amu = gm2calc_thdm_calculate_uncertainty_amu_1loop_amu1L_amu2L(model, amu1L, amu2L);
   } else if (config_flags.loopOrder > 1) {
      delta_amu = gm2calc_thdm_calculate_uncertainty_amu_2loop_amu1L_amu2L(model, amu1L, amu2L);
   }

   return delta_amu;
}

/******************************************************************/

static void print_package()
{
   static int do_print = 1;

   if (do_print) {
      printf("===========================\n");
      printf("GM2Calc " GM2CALC_VERSION "\n");
      printf("http://gm2calc.hepforge.org\n");
      printf("===========================\n");

      do_print = 0;
   }
}

/******************************************************************/

int GM2CalcSetFlags(int loopOrder_, int tanBetaResummation_, int forceOutput_,
                    int runningCouplings_)
{
   char loop_order_str[12];

   print_package();

   if (loopOrder_ < 0 || loopOrder_ > 2) {
      snprintf(loop_order_str, sizeof(loop_order_str), "%d", loopOrder_);
      put_error_message("GM2CalcSetFlags", "wronglooporder", loop_order_str);
   }

   config_flags.loopOrder = loopOrder_;
   config_flags.tanBetaResummation = tanBetaResummation_;
   config_flags.forceOutput = forceOutput_;
   config_flags.runningCouplings = runningCouplings_;

   return 0;
}

/******************************************************************/

void GM2CalcGetFlags(void)
{
   MLPutFunction(stdlink, "List", 4);
   MLPutRuleToInteger(stdlink, config_flags.loopOrder, "loopOrder");
   MLPutRuleToString(stdlink, config_flags.tanBetaResummation ? "True" : "False", "tanBetaResummation");
   MLPutRuleToString(stdlink, config_flags.forceOutput ? "True" : "False", "forceOutput");
   MLPutRuleToString(stdlink, config_flags.runningCouplings ? "True" : "False", "runningCouplings");
   MLEndPacket(stdlink);
}

/******************************************************************/

int GM2CalcSetSMParameters(
   double alpha0_,
   double alphaMZ_,
   double alphaS_,
   double MhSM_,
   double MW_,
   double MZ_,
   double MT_,
   double mcmc_,
   double mu2GeV_,
   double mbmb_,
   double ms2GeV_,
   double md2GeV_,
   double Mv1_,
   double Mv2_,
   double Mv3_,
   double ML_,
   double MM_,
   double ME_,
   double CKM_real_11_,
   double CKM_real_12_,
   double CKM_real_13_,
   double CKM_real_21_,
   double CKM_real_22_,
   double CKM_real_23_,
   double CKM_real_31_,
   double CKM_real_32_,
   double CKM_real_33_,
   double CKM_imag_11_,
   double CKM_imag_12_,
   double CKM_imag_13_,
   double CKM_imag_21_,
   double CKM_imag_22_,
   double CKM_imag_23_,
   double CKM_imag_31_,
   double CKM_imag_32_,
   double CKM_imag_33_)
{
   sm.alpha_em_0 = alpha0_;
   sm.alpha_em_mz = alphaMZ_;
   sm.alpha_s_mz = alphaS_;
   sm.mh = MhSM_;
   sm.mw = MW_;
   sm.mz = MZ_;
   sm.mu[2] = MT_;
   sm.mu[1] = mcmc_;
   sm.mu[0] = mu2GeV_;
   sm.md[2] = mbmb_;
   sm.md[1] = ms2GeV_;
   sm.md[0] = md2GeV_;
   sm.mv[2] = Mv3_;
   sm.mv[1] = Mv2_;
   sm.mv[0] = Mv1_;
   sm.ml[2] = ML_;
   sm.ml[1] = MM_;
   sm.ml[0] = ME_;
   sm.ckm_real[0][0] = CKM_real_11_;
   sm.ckm_real[0][1] = CKM_real_12_;
   sm.ckm_real[0][2] = CKM_real_13_;
   sm.ckm_real[1][0] = CKM_real_21_;
   sm.ckm_real[1][1] = CKM_real_22_;
   sm.ckm_real[1][2] = CKM_real_23_;
   sm.ckm_real[2][0] = CKM_real_31_;
   sm.ckm_real[2][1] = CKM_real_32_;
   sm.ckm_real[2][2] = CKM_real_33_;
   sm.ckm_imag[0][0] = CKM_imag_11_;
   sm.ckm_imag[0][1] = CKM_imag_12_;
   sm.ckm_imag[0][2] = CKM_imag_13_;
   sm.ckm_imag[1][0] = CKM_imag_21_;
   sm.ckm_imag[1][1] = CKM_imag_22_;
   sm.ckm_imag[1][2] = CKM_imag_23_;
   sm.ckm_imag[2][0] = CKM_imag_31_;
   sm.ckm_imag[2][1] = CKM_imag_32_;
   sm.ckm_imag[2][2] = CKM_imag_33_;

   return 0;
}

/******************************************************************/

void GM2CalcGetSMParameters(void)
{
   MLPutFunction(stdlink, "List", 19);

   MLPutRuleToReal(stdlink, sm.alpha_em_0, "alpha0");
   MLPutRuleToReal(stdlink, sm.alpha_em_mz, "alphaMZ");
   MLPutRuleToReal(stdlink, sm.alpha_s_mz, "alphaS");
   MLPutRuleToReal(stdlink, sm.mh, "MhSM");
   MLPutRuleToReal(stdlink, sm.mw, "MW");
   MLPutRuleToReal(stdlink, sm.mz, "MZ");
   MLPutRuleToReal(stdlink, sm.mu[2], "MT");
   MLPutRuleToReal(stdlink, sm.mu[1], "mcmc");
   MLPutRuleToReal(stdlink, sm.mu[0], "mu2GeV");
   MLPutRuleToReal(stdlink, sm.md[2], "mbmb");
   MLPutRuleToReal(stdlink, sm.md[1], "ms2GeV");
   MLPutRuleToReal(stdlink, sm.md[0], "md2GeV");
   MLPutRuleToReal(stdlink, sm.mv[2], "Mv1");
   MLPutRuleToReal(stdlink, sm.mv[1], "Mv2");
   MLPutRuleToReal(stdlink, sm.mv[0], "Mv3");
   MLPutRuleToReal(stdlink, sm.ml[2], "ML");
   MLPutRuleToReal(stdlink, sm.ml[1], "MM");
   MLPutRuleToReal(stdlink, sm.ml[0], "ME");
   MLPutRuleToComplexMatrixTHDM(stdlink, sm.ckm, "CKM", 3, 3);

   MLEndPacket(stdlink);
}

/******************************************************************/

void create_result_list(MSSMNoFV_onshell* model)
{
   double amu1L = 0, amu2L = 0;
   calculate_amu_mssmnofv(model, &amu1L, &amu2L);
   const double amu = amu1L + amu2L;
   const double Damu = calculate_uncertainty_mssmnofv(model, amu1L, amu2L);

   MLPutFunction(stdlink, "List", 47);
   /* amu [4] */
   MLPutRuleToReal(stdlink, amu, "amu");
   MLPutRuleToReal(stdlink, amu1L, "amu1L");
   MLPutRuleToReal(stdlink, amu2L, "amu2L");
   MLPutRuleToReal(stdlink, Damu, "Damu");
   /* couplings [3] */
   MLPutRuleToReal(stdlink, sqr(gm2calc_mssmnofv_get_EL(model))/(4.*M_PI), "alphaMZ");
   MLPutRuleToReal(stdlink, sqr(gm2calc_mssmnofv_get_EL0(model))/(4.*M_PI), "alpha0");
   MLPutRuleToReal(stdlink, sqr(gm2calc_mssmnofv_get_g3(model))/(4.*M_PI), "alphaS");
   /* on-shell masses and parameters [40] */
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MM(model), "MM");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MT(model), "MT");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MBMB(model), "mbmb");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MB(model), "mbMZ");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_ML(model), "ML");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MW(model), "MW");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MZ(model), "MZ");
   MLPutRuleToRealVectorInterface(stdlink, MSm, "MSm", 2);
   MLPutRuleToRealMatrixInterface(stdlink, USm, "USm", 2, 2);
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MSvmL(model), "MSvmL");
   MLPutRuleToRealVectorInterface(stdlink, MStau, "MStau", 2);
   MLPutRuleToRealMatrixInterface(stdlink, UStau, "UStau", 2, 2);
   MLPutRuleToRealVectorInterface(stdlink, MSt, "MSt", 2);
   MLPutRuleToRealMatrixInterface(stdlink, USt, "USt", 2, 2);
   MLPutRuleToRealVectorInterface(stdlink, MSb, "MSb", 2);
   MLPutRuleToRealMatrixInterface(stdlink, USb, "USb", 2, 2);
   MLPutRuleToRealVectorInterface(stdlink, MCha, "MCha", 2);
   MLPutRuleToComplexMatrixInterface(stdlink, UM, "UM", 2, 2);
   MLPutRuleToComplexMatrixInterface(stdlink, UP, "UP", 2, 2);
   MLPutRuleToRealVectorInterface(stdlink, MChi, "MChi", 4);
   MLPutRuleToComplexMatrixInterface(stdlink, ZN, "ZN", 4, 4);
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MAh(model), "MAh");
   MLPutRuleToRealVectorInterface(stdlink, Mhh, "Mhh", 2);
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_TB(model), "TB");
   MLPutRuleToRealMatrixInterface(stdlink, Yu, "Yu", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, Yd, "Yd", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, Ye, "Ye", 3, 3);
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_Mu(model), "Mu");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MassB(model), "MassB");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MassWB(model), "MassWB");
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_MassG(model), "MassG");
   MLPutRuleToRealMatrixInterface(stdlink, mq2, "mq2", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, ml2, "ml2", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, mu2, "mu2", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, md2, "md2", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, me2, "me2", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, Au, "Au", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, Ad, "Ad", 3, 3);
   MLPutRuleToRealMatrixInterface(stdlink, Ae, "Ae", 3, 3);
   MLPutRuleToReal(stdlink, gm2calc_mssmnofv_get_scale(model), "Q");
   MLEndPacket(stdlink);
}

/******************************************************************/

void create_error_output(void)
{
   MLPutFunction(stdlink, "List", 0);
   MLEndPacket(stdlink);
}

/******************************************************************/

void GM2CalcAmuSLHAScheme(
   double MSvmL_,
   double MSm_1_,
   double MSm_2_,
   double MChi_1_,
   double MChi_2_,
   double MChi_3_,
   double MChi_4_,
   double MCha_1_,
   double MCha_2_,
   double MAh_,
   double TB_,
   double Mu_,
   double MassB_,
   double MassWB_,
   double MassG_,
   double mq2_11_,
   double mq2_22_,
   double mq2_33_,
   double ml2_11_,
   double ml2_22_,
   double ml2_33_,
   double mu2_11_,
   double mu2_22_,
   double mu2_33_,
   double md2_11_,
   double md2_22_,
   double md2_33_,
   double me2_11_,
   double me2_22_,
   double me2_33_,
   double Au_33_,
   double Ad_33_,
   double Ae_22_,
   double Ae_33_,
   double Q_
)
{
   struct SLHA_parameters pars;
   initialize_slha_parameters(&pars);

   pars.MSvmL     = MSvmL_;
   pars.MSm[0]    = MSm_1_;
   pars.MSm[1]    = MSm_2_;
   pars.MChi[0]   = MChi_1_;
   pars.MChi[1]   = MChi_2_;
   pars.MChi[2]   = MChi_3_;
   pars.MChi[3]   = MChi_4_;
   pars.MCha[0]   = MCha_1_;
   pars.MCha[1]   = MCha_2_;
   pars.MAh       = MAh_;
   pars.TB        = TB_;
   pars.Mu        = Mu_;
   pars.MassB     = MassB_;
   pars.MassWB    = MassWB_;
   pars.MassG     = MassG_;
   pars.mq2[0][0] = mq2_11_;
   pars.mq2[1][1] = mq2_22_;
   pars.mq2[2][2] = mq2_33_;
   pars.ml2[0][0] = ml2_11_;
   pars.ml2[1][1] = ml2_22_;
   pars.ml2[2][2] = ml2_33_;
   pars.mu2[0][0] = mu2_11_;
   pars.mu2[1][1] = mu2_22_;
   pars.mu2[2][2] = mu2_33_;
   pars.md2[0][0] = md2_11_;
   pars.md2[1][1] = md2_22_;
   pars.md2[2][2] = md2_33_;
   pars.me2[0][0] = me2_11_;
   pars.me2[1][1] = me2_22_;
   pars.me2[2][2] = me2_33_;
   pars.Au[2][2]  = Au_33_;
   pars.Ad[2][2]  = Ad_33_;
   pars.Ae[1][1]  = Ae_22_;
   pars.Ae[2][2]  = Ae_33_;
   pars.Q         = Q_;

   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();
   const gm2calc_error error = setup_model_slha_scheme(model, &pars);

   if (gm2calc_mssmnofv_have_warning(model)) {
      char msg[400];
      gm2calc_mssmnofv_get_warnings(model, msg, sizeof(msg));
      put_error_message("GM2CalcAmuSLHAScheme", "warning", msg);
   }

   if (error != gm2calc_NoError) {
      put_error_message("GM2CalcAmuSLHAScheme", "error",
                        gm2calc_error_str(error));
   }

   if (gm2calc_mssmnofv_have_problem(model)) {
      char msg[400];
      gm2calc_mssmnofv_get_problems(model, msg, sizeof(msg));
      put_error_message("GM2CalcAmuSLHAScheme", "error", msg);
   }

   if (error == gm2calc_NoError || config_flags.forceOutput) {
      create_result_list(model);
   } else {
      create_error_output();
   }

   gm2calc_mssmnofv_free(model);
}

/******************************************************************/

void GM2CalcAmuGM2CalcScheme(
   double MAh_,
   double TB_,
   double Mu_,
   double MassB_,
   double MassWB_,
   double MassG_,
   double mq2_11_,
   double mq2_22_,
   double mq2_33_,
   double ml2_11_,
   double ml2_22_,
   double ml2_33_,
   double mu2_11_,
   double mu2_22_,
   double mu2_33_,
   double md2_11_,
   double md2_22_,
   double md2_33_,
   double me2_11_,
   double me2_22_,
   double me2_33_,
   double Au_33_,
   double Ad_33_,
   double Ae_22_,
   double Ae_33_,
   double Q_)
{
   struct GM2Calc_parameters pars;
   initialize_gm2calc_parameters(&pars);

   pars.MAh       = MAh_;
   pars.TB        = TB_;
   pars.Mu        = Mu_;
   pars.MassB     = MassB_;
   pars.MassWB    = MassWB_;
   pars.MassG     = MassG_;
   pars.mq2[0][0] = mq2_11_;
   pars.mq2[1][1] = mq2_22_;
   pars.mq2[2][2] = mq2_33_;
   pars.ml2[0][0] = ml2_11_;
   pars.ml2[1][1] = ml2_22_;
   pars.ml2[2][2] = ml2_33_;
   pars.mu2[0][0] = mu2_11_;
   pars.mu2[1][1] = mu2_22_;
   pars.mu2[2][2] = mu2_33_;
   pars.md2[0][0] = md2_11_;
   pars.md2[1][1] = md2_22_;
   pars.md2[2][2] = md2_33_;
   pars.me2[0][0] = me2_11_;
   pars.me2[1][1] = me2_22_;
   pars.me2[2][2] = me2_33_;
   pars.Au[2][2]  = Au_33_;
   pars.Ad[2][2]  = Ad_33_;
   pars.Ae[1][1]  = Ae_22_;
   pars.Ae[2][2]  = Ae_33_;
   pars.Q         = Q_;

   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();
   const gm2calc_error error = setup_model_gm2calc_scheme(model, &pars);

   if (gm2calc_mssmnofv_have_warning(model)) {
      char msg[400];
      gm2calc_mssmnofv_get_warnings(model, msg, sizeof(msg));
      put_error_message("GM2CalcAmuGM2CalcScheme", "warning", msg);
   }

   if (error != gm2calc_NoError) {
      put_error_message("GM2CalcAmuGM2CalcScheme", "error",
                        gm2calc_error_str(error));
   }

   if (gm2calc_mssmnofv_have_problem(model)) {
      char msg[400];
      gm2calc_mssmnofv_get_problems(model, msg, sizeof(msg));
      put_error_message("GM2CalcAmuGM2CalcScheme", "error", msg);
   }

   if (error == gm2calc_NoError || config_flags.forceOutput) {
      create_result_list(model);
   } else {
      create_error_output();
   }

   gm2calc_mssmnofv_free(model);
}

/******************************************************************/

void GM2CalcAmuTHDMGaugeBasis(
   int yukawa_type_,
   double lambda_1_,
   double lambda_2_,
   double lambda_3_,
   double lambda_4_,
   double lambda_5_,
   double lambda_6_,
   double lambda_7_,
   double TB_,
   double m122_,
   double zeta_u_,
   double zeta_d_,
   double zeta_l_,
   double Deltau_11_,
   double Deltau_12_,
   double Deltau_13_,
   double Deltau_21_,
   double Deltau_22_,
   double Deltau_23_,
   double Deltau_31_,
   double Deltau_32_,
   double Deltau_33_,
   double Deltad_11_,
   double Deltad_12_,
   double Deltad_13_,
   double Deltad_21_,
   double Deltad_22_,
   double Deltad_23_,
   double Deltad_31_,
   double Deltad_32_,
   double Deltad_33_,
   double Deltal_11_,
   double Deltal_12_,
   double Deltal_13_,
   double Deltal_21_,
   double Deltal_22_,
   double Deltal_23_,
   double Deltal_31_,
   double Deltal_32_,
   double Deltal_33_,
   double Piu_11_,
   double Piu_12_,
   double Piu_13_,
   double Piu_21_,
   double Piu_22_,
   double Piu_23_,
   double Piu_31_,
   double Piu_32_,
   double Piu_33_,
   double Pid_11_,
   double Pid_12_,
   double Pid_13_,
   double Pid_21_,
   double Pid_22_,
   double Pid_23_,
   double Pid_31_,
   double Pid_32_,
   double Pid_33_,
   double Pil_11_,
   double Pil_12_,
   double Pil_13_,
   double Pil_21_,
   double Pil_22_,
   double Pil_23_,
   double Pil_31_,
   double Pil_32_,
   double Pil_33_)
{
   if (yukawa_type_ < 1 || yukawa_type_ > 6) {
      put_error_message("GM2CalcAmuTHDMGaugeBasis", "error",
                        "yukawaType must be between 1 and 6.");
      create_error_output();
      return;
   }

   gm2calc_THDM_gauge_basis basis;
   basis.yukawa_type = int_to_c_yukawa_type(yukawa_type_);
   basis.lambda[0] = lambda_1_;
   basis.lambda[1] = lambda_2_;
   basis.lambda[2] = lambda_3_;
   basis.lambda[3] = lambda_4_;
   basis.lambda[4] = lambda_5_;
   basis.lambda[5] = lambda_6_;
   basis.lambda[6] = lambda_7_;
   basis.tan_beta = TB_;
   basis.m122 = m122_;
   basis.zeta_u = zeta_u_;
   basis.zeta_d = zeta_d_;
   basis.zeta_l = zeta_l_;
   basis.Delta_u[0][0] = Deltau_11_;
   basis.Delta_u[0][1] = Deltau_12_;
   basis.Delta_u[0][2] = Deltau_13_;
   basis.Delta_u[1][0] = Deltau_21_;
   basis.Delta_u[1][1] = Deltau_22_;
   basis.Delta_u[1][2] = Deltau_23_;
   basis.Delta_u[2][0] = Deltau_31_;
   basis.Delta_u[2][1] = Deltau_32_;
   basis.Delta_u[2][2] = Deltau_33_;
   basis.Delta_d[0][0] = Deltad_11_;
   basis.Delta_d[0][1] = Deltad_12_;
   basis.Delta_d[0][2] = Deltad_13_;
   basis.Delta_d[1][0] = Deltad_21_;
   basis.Delta_d[1][1] = Deltad_22_;
   basis.Delta_d[1][2] = Deltad_23_;
   basis.Delta_d[2][0] = Deltad_31_;
   basis.Delta_d[2][1] = Deltad_32_;
   basis.Delta_d[2][2] = Deltad_33_;
   basis.Delta_l[0][0] = Deltal_11_;
   basis.Delta_l[0][1] = Deltal_12_;
   basis.Delta_l[0][2] = Deltal_13_;
   basis.Delta_l[1][0] = Deltal_21_;
   basis.Delta_l[1][1] = Deltal_22_;
   basis.Delta_l[1][2] = Deltal_23_;
   basis.Delta_l[2][0] = Deltal_31_;
   basis.Delta_l[2][1] = Deltal_32_;
   basis.Delta_l[2][2] = Deltal_33_;
   basis.Pi_u[0][0] = Piu_11_;
   basis.Pi_u[0][1] = Piu_12_;
   basis.Pi_u[0][2] = Piu_13_;
   basis.Pi_u[1][0] = Piu_21_;
   basis.Pi_u[1][1] = Piu_22_;
   basis.Pi_u[1][2] = Piu_23_;
   basis.Pi_u[2][0] = Piu_31_;
   basis.Pi_u[2][1] = Piu_32_;
   basis.Pi_u[2][2] = Piu_33_;
   basis.Pi_d[0][0] = Pid_11_;
   basis.Pi_d[0][1] = Pid_12_;
   basis.Pi_d[0][2] = Pid_13_;
   basis.Pi_d[1][0] = Pid_21_;
   basis.Pi_d[1][1] = Pid_22_;
   basis.Pi_d[1][2] = Pid_23_;
   basis.Pi_d[2][0] = Pid_31_;
   basis.Pi_d[2][1] = Pid_32_;
   basis.Pi_d[2][2] = Pid_33_;
   basis.Pi_l[0][0] = Pil_11_;
   basis.Pi_l[0][1] = Pil_12_;
   basis.Pi_l[0][2] = Pil_13_;
   basis.Pi_l[1][0] = Pil_21_;
   basis.Pi_l[1][1] = Pil_22_;
   basis.Pi_l[1][2] = Pil_23_;
   basis.Pi_l[2][0] = Pil_31_;
   basis.Pi_l[2][1] = Pil_32_;
   basis.Pi_l[2][2] = Pil_33_;

   gm2calc_THDM_config config;
   config.force_output = config_flags.forceOutput;
   config.running_couplings = config_flags.runningCouplings;

   gm2calc_THDM* model = 0;
   gm2calc_error error = gm2calc_thdm_new_with_gauge_basis(&model, &basis, &sm, &config);

   if (error == gm2calc_NoError) {
      double amu1L = 0, amu2LF = 0, amu2LB = 0;
      calculate_amu_thdm(model, &amu1L, &amu2LF, &amu2LB);
      const double damu = calculate_uncertainty_thdm(model, amu1L, amu2LF + amu2LB);

      MLPutFunction(stdlink, "List", 5);
      MLPutRuleToReal(stdlink, amu1L + amu2LF + amu2LB, "amu");
      MLPutRuleToReal(stdlink, amu1L, "amu1L");
      MLPutRuleToReal(stdlink, amu2LF, "amu2LF");
      MLPutRuleToReal(stdlink, amu2LB, "amu2LB");
      MLPutRuleToReal(stdlink, damu, "Damu");
      MLEndPacket(stdlink);
   } else {
      put_error_message("GM2CalcAmuTHDMMassBasis", "error",
                        gm2calc_error_str(error));
      create_error_output();
   }

   gm2calc_thdm_free(model);
}

/******************************************************************/

void GM2CalcAmuTHDMMassBasis(
   int yukawa_type_,
   double Mhh_1_,
   double Mhh_2_,
   double MAh_,
   double MHp_,
   double sin_beta_minus_alpha_,
   double lambda_6_,
   double lambda_7_,
   double TB_,
   double m122_,
   double zeta_u_,
   double zeta_d_,
   double zeta_l_,
   double Deltau_11_,
   double Deltau_12_,
   double Deltau_13_,
   double Deltau_21_,
   double Deltau_22_,
   double Deltau_23_,
   double Deltau_31_,
   double Deltau_32_,
   double Deltau_33_,
   double Deltad_11_,
   double Deltad_12_,
   double Deltad_13_,
   double Deltad_21_,
   double Deltad_22_,
   double Deltad_23_,
   double Deltad_31_,
   double Deltad_32_,
   double Deltad_33_,
   double Deltal_11_,
   double Deltal_12_,
   double Deltal_13_,
   double Deltal_21_,
   double Deltal_22_,
   double Deltal_23_,
   double Deltal_31_,
   double Deltal_32_,
   double Deltal_33_,
   double Piu_11_,
   double Piu_12_,
   double Piu_13_,
   double Piu_21_,
   double Piu_22_,
   double Piu_23_,
   double Piu_31_,
   double Piu_32_,
   double Piu_33_,
   double Pid_11_,
   double Pid_12_,
   double Pid_13_,
   double Pid_21_,
   double Pid_22_,
   double Pid_23_,
   double Pid_31_,
   double Pid_32_,
   double Pid_33_,
   double Pil_11_,
   double Pil_12_,
   double Pil_13_,
   double Pil_21_,
   double Pil_22_,
   double Pil_23_,
   double Pil_31_,
   double Pil_32_,
   double Pil_33_)
{
   if (yukawa_type_ < 1 || yukawa_type_ > 6) {
      put_error_message("GM2CalcAmuTHDMMassBasis", "error",
                        "yukawaType must be between 1 and 6.");
      create_error_output();
      return;
   }

   gm2calc_THDM_mass_basis basis;
   basis.yukawa_type = int_to_c_yukawa_type(yukawa_type_);
   basis.mh = Mhh_1_;
   basis.mH = Mhh_2_;
   basis.mA = MAh_;
   basis.mHp = MHp_;
   basis.sin_beta_minus_alpha = sin_beta_minus_alpha_;
   basis.lambda_6 = lambda_6_;
   basis.lambda_7 = lambda_7_;
   basis.tan_beta = TB_;
   basis.m122 = m122_;
   basis.zeta_u = zeta_u_;
   basis.zeta_d = zeta_d_;
   basis.zeta_l = zeta_l_;
   basis.Delta_u[0][0] = Deltau_11_;
   basis.Delta_u[0][1] = Deltau_12_;
   basis.Delta_u[0][2] = Deltau_13_;
   basis.Delta_u[1][0] = Deltau_21_;
   basis.Delta_u[1][1] = Deltau_22_;
   basis.Delta_u[1][2] = Deltau_23_;
   basis.Delta_u[2][0] = Deltau_31_;
   basis.Delta_u[2][1] = Deltau_32_;
   basis.Delta_u[2][2] = Deltau_33_;
   basis.Delta_d[0][0] = Deltad_11_;
   basis.Delta_d[0][1] = Deltad_12_;
   basis.Delta_d[0][2] = Deltad_13_;
   basis.Delta_d[1][0] = Deltad_21_;
   basis.Delta_d[1][1] = Deltad_22_;
   basis.Delta_d[1][2] = Deltad_23_;
   basis.Delta_d[2][0] = Deltad_31_;
   basis.Delta_d[2][1] = Deltad_32_;
   basis.Delta_d[2][2] = Deltad_33_;
   basis.Delta_l[0][0] = Deltal_11_;
   basis.Delta_l[0][1] = Deltal_12_;
   basis.Delta_l[0][2] = Deltal_13_;
   basis.Delta_l[1][0] = Deltal_21_;
   basis.Delta_l[1][1] = Deltal_22_;
   basis.Delta_l[1][2] = Deltal_23_;
   basis.Delta_l[2][0] = Deltal_31_;
   basis.Delta_l[2][1] = Deltal_32_;
   basis.Delta_l[2][2] = Deltal_33_;
   basis.Pi_u[0][0] = Piu_11_;
   basis.Pi_u[0][1] = Piu_12_;
   basis.Pi_u[0][2] = Piu_13_;
   basis.Pi_u[1][0] = Piu_21_;
   basis.Pi_u[1][1] = Piu_22_;
   basis.Pi_u[1][2] = Piu_23_;
   basis.Pi_u[2][0] = Piu_31_;
   basis.Pi_u[2][1] = Piu_32_;
   basis.Pi_u[2][2] = Piu_33_;
   basis.Pi_d[0][0] = Pid_11_;
   basis.Pi_d[0][1] = Pid_12_;
   basis.Pi_d[0][2] = Pid_13_;
   basis.Pi_d[1][0] = Pid_21_;
   basis.Pi_d[1][1] = Pid_22_;
   basis.Pi_d[1][2] = Pid_23_;
   basis.Pi_d[2][0] = Pid_31_;
   basis.Pi_d[2][1] = Pid_32_;
   basis.Pi_d[2][2] = Pid_33_;
   basis.Pi_l[0][0] = Pil_11_;
   basis.Pi_l[0][1] = Pil_12_;
   basis.Pi_l[0][2] = Pil_13_;
   basis.Pi_l[1][0] = Pil_21_;
   basis.Pi_l[1][1] = Pil_22_;
   basis.Pi_l[1][2] = Pil_23_;
   basis.Pi_l[2][0] = Pil_31_;
   basis.Pi_l[2][1] = Pil_32_;
   basis.Pi_l[2][2] = Pil_33_;

   gm2calc_THDM_config config;
   config.force_output = config_flags.forceOutput;
   config.running_couplings = config_flags.runningCouplings;

   gm2calc_THDM* model = 0;
   gm2calc_error error = gm2calc_thdm_new_with_mass_basis(&model, &basis, &sm, &config);

   if (error == gm2calc_NoError) {
      double amu1L = 0, amu2LF = 0, amu2LB = 0;
      calculate_amu_thdm(model, &amu1L, &amu2LF, &amu2LB);
      const double damu = calculate_uncertainty_thdm(model, amu1L, amu2LF + amu2LB);

      MLPutFunction(stdlink, "List", 5);
      MLPutRuleToReal(stdlink, amu1L + amu2LF + amu2LB, "amu");
      MLPutRuleToReal(stdlink, amu1L, "amu1L");
      MLPutRuleToReal(stdlink, amu2LF, "amu2LF");
      MLPutRuleToReal(stdlink, amu2LB, "amu2LB");
      MLPutRuleToReal(stdlink, damu, "Damu");
      MLEndPacket(stdlink);
   } else {
      put_error_message("GM2CalcAmuTHDMMassBasis", "error",
                        gm2calc_error_str(error));
      create_error_output();
   }

   gm2calc_thdm_free(model);
}

/******************************************************************/

int main(int argc, char *argv[])
{
   gm2calc_sm_set_to_default(&sm);
   return MLMain(argc, argv);
}
