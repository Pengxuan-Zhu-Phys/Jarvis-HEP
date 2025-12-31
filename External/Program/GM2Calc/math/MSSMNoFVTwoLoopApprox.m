(***************************************************************************************
* Mathematica code for the one-loop and two-loop approximations for 
* the MSSM contributions to the muon anomalous magnetic moment 
*
*  22/10/2013
*  Helvecio Fargnoli, Christoph Gnendiger, Sebastian Passehr,
*  Dominik Stoeckinger, Hyejung Stoeckinger-Kim
*
* This code provides the following definitions:
*
* amu1Lapprox: the one-loop contributions to muon (g-2) 
*              in mass-insertion approximation, see Sec. 5 of Ref. [1]. 
* amu2LFSfapprox: the leading logarithmic approximation of the 
*              fermion/sfermion two-loop contributions, see Sec. 5 of Ref. [1]. 
*
* If you use amu2LFSfapprox, please cite the following references:
*
*   [1]  Helvecio Fargnoli, Christoph Gnendiger, Sebastian Passehr, Dominik St\"ockinger, Hyejung St\"ockinger-Kim
*        % ``Two-Loop Corrections to the Muon Magnetic Moment from Fermion/Sfermion
*        % Loops in the MSSM: Full Results,''
*        [arXiv:1311.1775].
*   [2]  H.~Fargnoli, C.~Gnendiger, S.~Passehr, D.~St\"ockinger and H.~St\"ockinger-Kim,
*        % ``Non-decoupling two-loop corrections to (g-2)_\mu\ from fermion/sfermion loops in the MSSM,''
*        Physics Letters B 726 (2013), pp. 717-724, 
*        [arXiv:1309.0980 [hep-ph]].
*
* Comments: [1] defines and introduces the formulas contained in the present code;
*           [2] first introduces the fermion/sfermion-loop contributions and
*               explains the leading logarithmic behaviour.
***************************************************************************************
*
* Remarks:
* -The 1st and 2nd generation Yukawa couplings are neglected. 
* -The non-logarithmic numerical constants approximate the magnitude of the 
*   additional non-logarithmic contributions.
* -The parameters are similar to the notation of the FeynArts MSSM model file, i.e.
*   MM, ML, MT, MB, MW, MZ, CW, SW: masses of muon, tau, top, bottom, W-boson, Z-boson, cos/sin of theta_W.
*   M1, M2, MUE, TB, CB, SB: masses of bino, wino, superpotential mu-parameter, tanbeta, cosbeta, sinbeta.
*   MSU[i], MSD[i], MSQ[i], MSE[i], MSL[i]: soft susy breaking mass parameters for the 
*                        squark and slepton doublets and singlets of generation i (i=1,2,3).
*   MSf[1,1,2]: muon-sneutrino mass.
*   AlfaMZ: running fine structure constant at scale MZ.
*   MUDIM2: square of the dimensional regularization/reduction scale.
*   
*
***************************************************************************************)

(*-------------------DEFINITIONS------------------------*) 
g1Coupling = (2*Sqrt[AlfaMZ]*Sqrt[Pi])/CW;
g2Coupling = (2*Sqrt[AlfaMZ]*Sqrt[Pi])/SW;
Fa[x_, y_] := -((G3[x] - G3[y])/(x - y)); 
Fb[x_, y_] := -((G4[x] - G4[y])/(x - y));
G3[x_] := (1/(2*(x - 1)^3))*((x - 1)*(x - 3) + 2*Log[x]); 
G4[x_] := (1/(2*(x - 1)^3))*((x - 1)*(x + 1) - 2*x*Log[x]);
LogScale = Min[Abs[MUE], Abs[M1], Abs[M2], Abs[MSL[2]], Abs[MSE[2]]];
(* ---------------------------------------------------- *)

(*---------------One- and two-loop results--------------*)
amu1Lapprox = (
	       (g1Coupling^2*MM^2*MUE*TB*Fb[MSL[2]^2/M1^2, MSE[2]^2/M1^2])/(8*M1^3*Pi^2) 
	       - (g1Coupling^2*M1*MM^2*MUE*TB*Fb[M1^2/MSE[2]^2, MUE^2/MSE[2]^2])/(8*Pi^2*MSE[2]^4) 
	       + (g2Coupling^2*M2*MM^2*MUE*TB*Fa[M2^2/MSf[1, 1, 2]^2,MUE^2/MSf[1, 1, 2]^2])/(8*Pi^2*MSf[1, 1, 2]^4) 
	       + (g1Coupling^2*M1*MM^2*MUE*TB*Fb[M1^2/MSL[2]^2, MUE^2/MSL[2]^2])/(16*Pi^2*MSL[2]^4) 
	       - (g2Coupling^2*M2*MM^2*MUE*TB*Fb[M2^2/MSL[2]^2, MUE^2/MSL[2]^2])/(16*Pi^2*MSL[2]^4)
	       );
(*------------------------------------------------------*)
amu2LFSfapprox = (
	       (
		((0.03 + Deltag1 + DeltaTanBeta)*g1Coupling^2*MM^2*MUE*TB*
		 Fb[MSL[2]^2/M1^2, MSE[2]^2/M1^2])/(8*M1^3*Pi^2)
		) - 
	       (
		((0.04 + Deltag1 + DeltaTanBeta + DeltaYukBinoHiggsino + DeltaYukHiggsino)*
		 g1Coupling^2*M1*MM^2*MUE*TB*Fb[M1^2/MSE[2]^2,MUE^2/MSE[2]^2])/(8*Pi^2*MSE[2]^4)
		) + 
	       (
		((0.015 + Deltag2 + DeltaTanBeta + DeltaYukHiggsino + DeltaYukWinoHiggsino)*
		 g2Coupling^2*M2*MM^2*MUE*TB*Fa[M2^2/MSf[1, 1, 2]^2,MUE^2/MSf[1, 1, 2]^2])/(8*Pi^2*MSf[1, 1, 2]^4)
		) + 
	       (
		((0.015 + Deltag1 + DeltaTanBeta + DeltaYukBinoHiggsino + DeltaYukHiggsino)*
		 g1Coupling^2*M1*MM^2*MUE*TB*Fb[M1^2/MSL[2]^2,MUE^2/MSL[2]^2])/(16*Pi^2*MSL[2]^4)
		) - 
	       (
		((0.015 + Deltag2 + DeltaTanBeta + DeltaYukHiggsino + DeltaYukWinoHiggsino)*
		 g2Coupling^2*M2*MM^2*MUE*TB*Fb[M2^2/MSL[2]^2,MUE^2/MSL[2]^2])/(16*Pi^2*MSL[2]^4))
	       );
(*------------------------------------------------------*)

(*--------------The shifts Delta------------------------*)
Deltag1 = (AlfaMZ*((Log[MSD[1]/LogScale]+Log[MSD[2]/LogScale])/3 + Log[MSD[3]/LogScale]/3 + 
		   Log[MSE[3]/LogScale] + Log[MSL[3]/LogScale]/2 + Log[MSQ[1]/LogScale]/6 + Log[MSQ[2]/LogScale]/6 +
		   Log[MSQ[3]/LogScale]/6 + 4*(Log[MSU[1]/LogScale]+Log[MSU[2]/LogScale])/3 + 
		   (4*Log[MSU[3]/LogScale])/3))/(3*CW^2*Pi);
Deltag2 = (AlfaMZ*(Log[MSL[3]/LogScale]/2 + 3*Log[MSQ[1]/LogScale]/2 + 3*Log[MSQ[2]/LogScale]/2 +
		   (3*Log[MSQ[3]/LogScale])/2))/(3*Pi*SW^2);
DeltaYukHiggsino = ((6*AlfaMZ*MB^2*Pi*Log[MSD[3]/LogScale])/
		    (CB^2*MW^2*SW^2) + (2*AlfaMZ*ML^2*Pi*Log[MSE[3]/LogScale])/
		    (CB^2*MW^2*SW^2) + (2*AlfaMZ*ML^2*Pi*Log[MSL[3]/LogScale])/
		    (CB^2*MW^2*SW^2) + 3*((2*AlfaMZ*MB^2*Pi)/(CB^2*MW^2*SW^2) + 
					  (2*AlfaMZ*MT^2*Pi)/(MW^2*SB^2*SW^2))*Log[MSQ[3]/LogScale] + 
		    (6*AlfaMZ*MT^2*Pi*Log[MSU[3]/LogScale])/(MW^2*SB^2*SW^2))/(32*Pi^2);

DeltaYukBinoHiggsino = ((4*AlfaMZ*MT^2*Pi*Log[MSQ[3]/LogScale])/
			(MW^2*SB^2*SW^2) - (16*AlfaMZ*MT^2*Pi*Log[MSU[3]/LogScale])/
			(MW^2*SB^2*SW^2))/(16*Pi^2);

DeltaYukWinoHiggsino = (
			(-3*AlfaMZ*MT^2*Log[MSQ[3]/LogScale])/
			(4*MW^2*Pi*SB^2*SW^2)
			);

DeltaTanBeta = -(((-((6 AlfaMZ MB^2 Pi)/(CB^2 MW^2 SW^2)) - (
    2 AlfaMZ ML^2 Pi)/(CB^2 MW^2 SW^2) + (6 AlfaMZ MT^2 Pi)/(
    MW^2 SB^2 SW^2)) Log[Sqrt[MUDIM2]/LogScale])/(16 Pi^2));
(*------------------------------------------------------*)


