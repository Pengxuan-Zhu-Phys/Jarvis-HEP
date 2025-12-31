(* 1-loop contribution to (g-2)/2 in the THDM

  cosab = Cos[beta - alpha]

  vev ~ 245 GeV is the SM-like Higgs vacuum expectation value

  zetaL = \zeta_f = Y^f_2 / cos(\beta)  - \sqrt{2} tan(\beta) M_f / v,
  zetaL is related to Eq.(17) from arxiv:1607.06292 by
  y^h_f = M_f / v sin(\beta-\alpha) + \xi^f cos(\beta-\alpha) / \sqrt{2}
  y^H_f = M_f / v cos(\beta-\alpha) - \xi^f sin(\beta-\alpha) / \sqrt{2}
  y^A_{d,l} = I \xi^{d,l}_A / \sqrt{2}
  y^A_u = -I \xi^u_A / \sqrt{2}

  The Y^f are related to the Yukawa part of the Lagrangian:

     -L = Y^u_1 Q.H_1 u + Y^u_2 Q.H_2 u + Y^d_1 Q.H_1 d + Y^d_2 Q.H_2 d + Y^e_1 L.H_1 e + Y^e_2 L.H_d e

 *)

(* coupling of down-type leptons i and j to h *)
yh[i_, j_] := ml[[2]]/vev*Sqrt[1 - cosab^2]*KroneckerDelta[i, j] + cosab/Sqrt[2] xiL[[i, j]]


(* coupling of down-type leptons i and j to H *)
yH[i_, j_] := ml[[2]]/vev*cosab*KroneckerDelta[i, j] - Sqrt[1 - cosab^2]/Sqrt[2] xiL[[i, j]]


(* coupling of down-type leptons i and j to H^\pm *)
yHp[i_, j_] := xiL[[i, j]]


(* coupling of down-type leptons i and j to A *)
yA[i_, j_] := I/Sqrt[2] xiL[[i, j]]


(* h, H and A contribution *)
Aloop1S[li_, ml_, mphi2_, y_] :=
    F1C[ml[[li]]^2/mphi2]/12*(Abs[y[li, 2]]^2 + Abs[y[2, li]]^2) +
    2/3*ml[[li]]/ml[[2]]*F2C[ml[[li]]^2/mphi2]*Conjugate[y[li, 2]]*Conjugate[y[2, li]]


(* H^\pm contribution *)
Aloop1Hp[li_, mnu_, mphi2_, y_] :=
    Abs[y[li, 2]]^2 (F1N[mnu[[2]]^2/mphi2] + F1N[mnu[[li]]^2/mphi2])/24


Prefactor1L[mphi2_] := ml[[2]]^2/mphi2 / (4 Pi)^2


amu1L :=
    Sum[
        Prefactor1L[mh^2] Aloop1S[li, ml, mh^2, yh] +
        Prefactor1L[mH^2] Aloop1S[li, ml, mH^2, yH] +
        Prefactor1L[mA^2] Aloop1S[li, ml, mA^2, yA] +
        Prefactor1L[mHp^2] Aloop1Hp[li, mnu, mHp^2, yHp]
      , {li, 1, 3}
    ]
