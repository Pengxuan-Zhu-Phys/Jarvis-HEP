(*

  This file contains the fermionic non-Standard Model 2-loop
  contributions to the anomalous magnetic moment of the muon in the
  Two-Higgs Doublet Model from [arxiv:1607.06292].

 *)

(* Eq (60), arxiv:1607.06292 *)
FlHp[ms2_, mf2_] :=
    Module[{xl = mf2/ms2},
           (
               xl + xl (xl - 1) (PolyLog[2, 1 - 1/xl] - Pi^2/6)
               + (xl - 1/2) Log[xl]
           )
    ]

(* Eq (61), arxiv:1607.06292
 
  Note: There is a misprint in Eq (61), arxiv:1607.06292v2: There
  should be no Phi function in the 2nd line of (61).
 *)
FdHp[ms2_, md2_, mu2_, qd_, qu_] :=
    FdHp[mu2/ms2, md2/ms2, qd, qu]

FdHp[xu_, xd_, qd_, qu_] :=
    Module[{cbar = (xu - qu) xu - (xd + qd) xd,
            y = (xu - xd)^2 - 2 (xu + xd) + 1,
            c = (xu - xd)^2 - qu xu + qd xd,
            s = (qu + qd)/4},
           (
               - (xu - xd)
               + (cbar/y - c (xu - xd)/y) Phi[xd, xu, 1]
               + c (PolyLog[2, 1 - xd/xu] - Log[xu]*Log[xd/xu]/2)
               + (s + xd) Log[xd] + (s - xu) Log[xu]
           )
    ]

(* Eq (62), arxiv:1607.06292 *)
FuHp[ms2_, md2_, mu2_, qd_, qu_] :=
    FuHp[mu2/ms2, md2/ms2, qd, qu]

FuHp[xu_, xd_, qd_, qu_] :=
    Module[{y = (xu - xd)^2 - 2 (xu + xd) + 1},
           (
               FdHp[xu, xd, 2 + qd, 2 + qu]
               - 4/3 (xu - xd - 1)/y Phi[xd, xu, 1]
               - 1/3 (Log[xd]^2 - Log[xu]^2)
           )
    ]
