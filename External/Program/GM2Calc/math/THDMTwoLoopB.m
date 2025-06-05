(*

  This file contains the bosonic non-Standard Model 2-loop contributions
  to the anomalous magnetic moment of the muon in the Two-Higgs Doublet
  Model from [arxiv:1607.06292].

  amu2LBNonYuk : non-Yukawa contributions from Eq. (71)

  amu2LBYuk : Yukawa contributions from Eq. (52)

  amu2LBEW : electro-weak contributions from Eq. (49)

 *)

expandAmu = {
    CW2 -> CW^2,
    SW2 -> 1 - CW2,
    xA -> MA0^2 / MZ^2,
    xHp -> MHp^2 / MZ^2,
    xhSM -> MH^2 / MZ^2,
    xH -> MHH^2 / MZ^2,
    f1 -> 7/2 - 25/(2*CW2) + 4*CW2 - 4*CW2^2,
    f2 -> 2*(17 - 24*CW2 + 56*CW2^2 - 128*CW2^3 + 64*CW2^4),
    f3 -> (25 - 32*CW2 + 4*CW2^2)/(CW2*SW2),
    f4 -> 13/2 - 15*CW2 + 10*CW2^2,
    f5 -> CW2*(5 - 16*CW2 + 8*CW2^2)/SW2,
    f6 -> (7 - 14*CW2 + 4*CW2^2)/(4*CW2*SW2),
    f7 -> 1 - 6*CW2 + 4*CW2^2,
    f8 -> (13 - 20*CW2 + 4*CW2^2)/(CW2*SW2),
    f9 -> 7 - 12*CW2 + 8*CW2^2
}


(* Eq.(72), arxiv:1607.06292 *)
T0[u_, w_] :=
    9/CW2^2*(u - w)*(CW2*(u - w)*(u + 2*w) - (u - w)^3 + CW2^2*w)/(CW2^2 + (u - w)^2 - 2*CW2*(u + w))*Phi[u, w, CW2]


(* Eq.(73), arxiv:1607.06292 *)
T1[u_, w_] :=
    9/CW2^2*(u - w)*(CW2*w - (u - w)^2)*PolyLog[2, 1 - u/w]


(* Eq.(74), arxiv:1607.06292 *)
T2[u_, w_, sgn_] :=
    Log[u]*(
      + (6*u^2 + CW2*(u - xHp) + 2*CW2^2*(u - xHp))/(2*(u - w))
      + f6*(u - xHp)^2*(3*CW2^2 + 3*CW2*(u - xHp) + (u - xHp)^2)/(CW2*(u - w))
      + sgn*f7*3*u^2*(u - xHp)/((xA - xH)*(u - w))
      - f8*3*u*(u - xHp)^2/(2*(u - w))
      - f9*3*u*(u - xHp)/(2*(u - w))
    )

T2p[args__] := T2[args, +1]

T2m[args__] := T2[args, -1]


(* Eq.(75), arxiv:1607.06292 *)
T4[u_, w_] :=
    (u - w)*Log[u]/4*f5*(xA*(3 + 2*xH) - xA^2 + 3*xH - xH^2 - 3)


(* Eq.(76), arxiv:1607.06292 *)
T5[u_, w_] :=
    Log[u]*(
      3/2*u + f6/CW2*((u - w)^3 + 3*CW2*(u - w)^2 + 3*CW2^2*(u - w))
      - 3/2*f8*u*(u - w) - CW2/2 - CW2^2)


(* Eq.(77), arxiv:1607.06292 *)
T6[u_, w_] :=
    9/2*((u - w)*(u^2 - 2*u*w + w*(w - CW2))/CW2^2*Log[u/w]*Log[w/CW2]
         + Log[CW2]/CW2*(2*u^2 + u*(CW2 - 4*w) - w*(CW2 - 2*w)))


(* Eq.(78), arxiv:1607.06292 *)
(* Note: overall factor -1/2 is missing in arxiv:1607.06292v2 *)
T7[u_, w_] :=
    Module[{ s1 = u + w - 1 + Sqrt[1 + (u - w)^2 - 2*(u + w)] },
           -1/2*f5*(2*(u + w) - (u - w)^2 - 1)*Log[s1/(2*Sqrt[u*w])]*(u + w - 1 - 4*u*w/s1)
    ]


(* Eq.(79), arxiv:1607.06292 *)
T8[u_, w_] :=
    Module[{ s2 = u + w - CW2 + Sqrt[(u + w - CW2)^2 - 4*u*w] },
           2*f6*(4*u*w - (u + w - CW2)^2)*Log[s2/(2*Sqrt[u*w])]*((u + w)/CW2 - 4*u*w/(CW2*s2) - 1)
    ]


(* Eq (71), arxiv:1607:06292 *)
amu2LBNonYuk = (AL/(24*Pi*CW2*(1 - CW2)) MM/MZ)^2 (
    + (xA - xH)/(xA - xHp)*T2p[xA, xH]
    + T2m[xH, xHp]
    + (xA - xH)/(xA - xHp)*T4[xA, xHp]
    + T4[xH, xA]
    + T5[xHp, xH]
    + T5[xHp, xA]
    + T2p[xHp, xH]
    + T2p[xHp, xA]
    + T6[xA, xHp]
    + T6[xH, xHp]
    + T7[xA, xH]
    + T7[xHp, xHp]*(1 - 2*CW2)^2
    + T8[xA, xHp]
    + T8[xH, xHp]
    - 16/3*CW2*(1 - CW2)*(1 + 8*CW2 - 8*CW2^2)
    + 8*CW2^2*(1 - CW2)^2/(5*xHp)
    + f2*xHp
    - f3*xHp^2
    + f1*(xA^2 + xH^2)
    + f3*xHp*(xA + xH)
    + f4*(xA + xH)
    - f5*xA*xH
    + T1[xA, xHp]
    + T1[xH, xHp]
    + T0[xA, xHp]
    + T0[xH, xHp]
) //. expandAmu


(* Eq.(99), arxiv:1607.06292 *)
b[u_, w_] := AL*Pi/(CW2*(-1 + CW2))*(u + 2*w)


(* Eq.(100), arxiv:1607.06292 *)
Fm0[u_, w_] :=
    1/(AL*Pi) * CW2*(-1 + CW2)/(u + 2*w) * YF1[u, w]


(* Eq.(101), arxiv:1607.06292 *)
Fmp[u_, w_] :=
    (-9*(-1 + CW2))/(AL*Pi) * (T9[u, w]/2 + T10[u, w])


(* Eq.(121), arxiv:1607.06292 *)
YFZ[u_] :=
    Module[{
        z1 = 3*(17 - 48*CW2 + 32*CW2^2), (* Eq.(122) *)
        z2 = 5 - 12*CW2 + 8*CW2^2,       (* Eq.(123) *)
        z3 = 3*(1 - 3*CW2 + 2*CW2^2)     (* Eq.(124) *)
        },
        (
            + z1*u*PolyLog[2, 1 - u]
            + z2/(2*u^2)*(6*(-4 + u)*u + Pi^2*(4 + 3*u) + 6*u*(4 + u)*Log[u]
                          - 6*(4 + 3*u)*PolyLog[2, 1 - u] + 6*u*(2 + u)*Phi[u, 1, 1])
            + z3*u*(6 + Pi^2*(-4 + u)*u + 3*Log[u]*(4 + (-4 + u)*u*Log[u])
                    + 12*(-4 + u)*u*PolyLog[2, 1 - u] + 6*(-2 + u)*Phi[u, 1, 1])
        )
    ]


(* Eq.(125), arxiv:1607.06292 *)
YFW[u_] := (
    - 57/2*CW2 - 4*CW2^3*Pi^2/u^2 + 3*CW2^2*(32 - 3*Pi^2)/(4*u)
    + 3*(16*CW2^3 + 9*CW2^2*u + 12*CW2*u^2 - 19*u^3)*PolyLog[2, 1 - u/CW2]/(2*u^2)
    + 3*CW2*(16*CW2 + 19*u)*(Log[CW2/u])/(2*u)
    + 3*(4*CW2^2 - 50*CW2*u + 19*u^2)*Phi[u, CW2, CW2]/(2*(4*CW2-u)*u)
)

(* Eq.(102), arxiv:1607.06292 *)
YF1[u_, w_] := (
    - 72*CW2*(-1 + CW2)*(u + 2*w)/u - 36*CW2*(-1 + CW2)*(u + 2*w)/u*Log[w]
    + 9*(-8*CW2^2 - 3*u + 2*CW2*(4 + u))*(u + 2*w)/(2*(u-1)*u)*Log[u]
    - 9*(3 - 10*CW2 + 8*CW2^2)*w*(u + 2*w)/((4*w-1)*(u-1))*Phi[w, w, 1]
    + 9*(8*CW2^2 + 3*u - 2*CW2*(4 + u))*w*(u + 2*w)/((4*w-u)*(u-1)*u^2)*Phi[u, w, w]
)


(* Eq.(105), arxiv:1607.06292 *)
YF2[u_] :=
    Module[{
        f0 = 3/4*CW2^2*(-640 + 576*CW2 + 7*Pi^2),       (* Eq.(106) *)
        f1 = 96*CW2^3*(11 - 53*CW2 + 36*CW2^2),         (* Eq.(107) *)
        f2 = -3/4*CW2*(-66*CW2 - 48*CW2^2 + 672*CW2^3), (* Eq.(108) *)
        f3 = -3/4*CW2*(109 - 430*CW2 + 120*CW2^2),      (* Eq.(109) *)
        f4 = 96*CW2^3*(-11 + 9*CW2),                    (* Eq.(110) *)
        f5 = 45/2*CW2^2 + 192*CW2^3,                    (* Eq.(111) *)
        f6 = 3/4*CW2*(157 + 90*CW2),                    (* Eq.(112) *)
        f7 = -3/4*(18 + 61*CW2),                        (* Eq.(113) *)
        f8 = -7 + 61*CW2 - 162*CW2^2 + 96*CW2^3,        (* Eq.(114) *)
        f9 = 1 - 5*CW2 + 10*CW2^2,                      (* Eq.(115) *)
        f10 = -1728*CW2^4*(-1 + CW2),                   (* Eq.(116) *)
        f11 = 3*CW2^3*(-899 + 768*CW2),                 (* Eq.(117) *)
        f12 = 387*CW2^2 - 363*CW2^3,                    (* Eq.(118) *)
        f13 = 9/2*CW2*(57 + 106*CW2),                   (* Eq.(119) *)
        f14 = -15/2*(7 + 45*CW2)                        (* Eq.(120) *)
        },
        (
            + YFW[u]
            + YFZ[u]
            + 8*CW2^3*Pi^2/u^2 + f0/u + 393/8*CW2
            + (f1/u + f2 + f3*u)*Log[CW2]/((4*CW2-1)*(4*CW2-u))
            + (f4/u + f5 + f6*u + f7*u^2)*Log[u]/((u-1)*(4*CW2-u))
            - 3/2*(32*CW2^3/u^2 + 21*CW2^2/u + 15*CW2 - 35*u)*PolyLog[2, 1 - u/CW2]
            + (f8 + f9*u)*9*CW2*(-3 + 4*CW2)/2*Phi[CW2, CW2, 1]/((4*CW2-1)^2*(u-1))
            + (f10/u^2 + f11/u + f12 + f13*u + f14*u^2 + 105/2*u^3)*Phi[u, CW2, CW2]/
              ((4*CW2-u)^2*(u-1))
        )
    ]


(* Eq.(126), arxiv:1607.06292 *)
YF3[u_, w_] :=
    Module[{
        (* Eq.(127) *)
        a1 = -9*CW2*u^3 + 9*CW2*u^2*(3*CW2+w) + 27*CW2^2*u*(w-CW2) + 9*(CW2^4 - 4*CW2^3*w + 3*CW2^2*w^2),
        (* Eq.(128) *)
        a2 = 9*CW2^2*w/2 - 9*u^2*(5*CW2 + w) + u*(36*CW2^2 + 153*CW2*w/4) + 9*u^3,
        (* Eq.(129) *)
        a3 = 9*CW2*u^2 - 9/2*CW2*u*(4*CW2 + w),
        (* Eq.(130) *)
        a4 = -9/2*u^2*w*(2*CW2^2 + 9*CW2*w + 2*w^2) + 9/8*u*w*(32*CW2^3 + 13*CW2^2*w + 35*CW2*w^2) + 9*u^3*w^2,
        (* Eq.(131) *)
        a5 = -9*u^3*(CW2 + w) - 9*u*(3*CW2^3 + 2*CW2*w^2) + 9*u^2*(3*CW2^2 + 4*CW2*w + w^2) + 9/2*CW2^2*(2*CW2^2 - 6*CW2*w + w^2),
        (* Eq.(132) *)
        a6 = -9*u^4*(9*CW2 + w) + u*(81*CW2^3*w - 225*CW2^4) + 9*CW2^4*(w - CW2) - 9/2*u^2*(3*CW2^3 + 37*CW2^2*w) + u^3*(198*CW2^2 + 72*CW2*w) + 9*u^5,
        (* Eq.(133) *)
        a7 = -9*CW2*u^4 + 18*CW2*u^3*(2*CW2 + w) + 36*u*(CW2^4 - 2*CW2^3*w) - 9*CW2*u^2*(6*CW2^2 - CW2*w + w^2) - 9*CW2*(CW2 - 3*w)*(CW2^3 - 2*CW2^2*w + CW2*w^2)
        },
        (
            + 9*u*(2*CW2 - u + w)/w
            + (a1*(Log[u] - Log[CW2]) + 9*CW2^2*(CW2^2 - 4*CW2*w + 3*w^2)*Log[CW2])*(Log[w] - Log[CW2])/(2*w^2*(CW2-w))
            + a2*Log[u]/(w*(4*CW2-u))
            + a3*Log[w]/(w*(CW2-w))
            + a4*Log[CW2]/(w^2*(4*CW2-u)*(CW2-w))
            + a5/(CW2*w^2)*PolyLog[2, 1 - u/CW2]
            + a6/(u*CW2*(4*CW2-u)^2*(CW2-w))*Phi[u, CW2, CW2]
            + a7/(w^2*(CW2-w)*(CW2^2-2*CW2*(u+w)+(u-w)^2))*Phi[u, w, CW2]
        )
    ]


(* Eq.(103), arxiv:1607.06292 *)
T9[u_, w_] := (
    - 2*(CW2^2*w + CW2*(u^2 + u*w - 2*w^2) - (u-w)^3)*Phi[u, w, CW2]/
      ((CW2 - w)*(CW2^2 - 2*CW2*(u+w) + (u-w)^2))
    + 2*CW2^2*(u^2 - 4*u*w + 2*w^2)*Phi[u, w, w]/(w^2*(w-CW2)*(u-4*w))
    - 2*(CW2*u*(u-2*w) + w*(u-w)^2)*PolyLog[2, 1 - u/w]/w^2
)


(* Eq.(104), arxiv:1607.06292 *)
T10[u_, w_] := (
    (u^2 - CW2*w - 2*u*w + w^2)/(2*(CW2-w))*Log[w/u]*Log[w/CW2]
    + CW2*(CW2 + 2*u - 2*w)/(2*(CW2-w))*Log[w/CW2]
    + CW2*u/w*Log[w/u]
    + CW2/w*(w-u)
)


(* Eq (52), arxiv:1607:06292 *)
amu2LBYuk =
    Module[{
        (* Eq.(91), arxiv:1607.06292 *)
        a000 = b[xhSM, xHp]*Fm0[xhSM, xHp],
        (* Eq.(92), arxiv:1607.06292 *)
        a0z0 = -b[xH, 0]*(Fm0[xH, xHp] + Fmp[xH, xHp]),
        (* Eq.(93), arxiv:1607.06292 *)
        a500 = Fm0[xhSM, xHp],
        (* Eq.(94), arxiv:1607.06292 *)
        a5z0 = -1/2 (Fm0[xH, xHp] + Fmp[xH, xHp]),
        (* Eq.(95), arxiv:1607.06292 *)
        a001 = (+ b[xH, 0]*Fm0[xH, xHp]
                - b[xhSM, 0]*Fm0[xhSM, xHp]),
        (* Eq.(96), arxiv:1607.06292 *)
        a0z1 = (
            - (
                + b[xH, xHp]*(
                    + Fm0[xH, xHp]
                    + Fmp[xH, xHp]
                  )
                - YF3[xH, xHp]
                (* SM contributions with opposite sign: *)
                - b[xhSM, xHp]*(
                    + Fm0[xhSM, xHp]
                    + Fmp[xhSM, xHp]
                  )
                + YF3[xhSM, xHp]
            ) + YF2[xH] ),
        (* Eq.(97), arxiv:1607.06292 *)
        a501 = (+ Fm0[xH, xHp]/2
                - Fm0[xhSM, xHp]/2),
        (* Eq.(98), arxiv:1607.06292 *)
        a5z1 = (- Fm0[xH, xHp]
                - Fmp[xH, xHp]
                + Fm0[xhSM, xHp]
                + Fmp[xhSM, xHp])
        },
        (AL/(24*Pi*CW2*(1 - CW2)) MM/MZ)^2 (
            + a000
            + a0z0*(TB - 1/TB)*ZetaL
            + a500*Lambda5
            + a5z0*(TB - 1/TB)*Lambda567*ZetaL
            + (
                + a001*(TB - 1/TB)
                + a0z1*ZetaL
                + a501*(TB - 1/TB)*Lambda567
                + a5z1*Lambda5*ZetaL
            )*aeps
        ) //. expandAmu
    ]


(* expansion rules of 2-loop EW bosonic contributions *)
expandEW = {
    pi2 -> Pi^2,
    CW2 -> CW^2,
    CW4 -> CW2^2,
    CW6 -> CW2^3,
    CW8 -> CW2^4,
    CW10 -> CW2^5,
    CW12 -> CW2^6,
    CW14 -> CW2^7,
    xh2 -> xh^2,
    xh3 -> xh^3,
    xh4 -> xh^4,
    xh5 -> xh^5,
    xh6 -> xh^6,
    l32 -> l3^2,
    lc2 -> lc^2,
    lh2 -> lh^2,
    xh -> (MH/MZ)^2,
    s0 -> Sqrt[1 - 4*CW2],
    s1 -> Sqrt[-4*xh + xh2],
    s2 -> Sqrt[-((4*CW2 - xh)*xh)],
    lh -> Log[xh],
    lc -> Log[CW2],
    l1 -> Log[(2 - s1 - xh)/2],
    l2 -> Log[(-s1 + xh)/2],
    l3 -> Log[(1 - s0)/2],
    l4 -> Log[(2 - (s2 + xh)/CW2)/2],
    l5 -> Log[(-s2 + xh)/(2 CW2)],
    li1 -> PolyLog[2, (xh - s1)/2],
    li2 -> PolyLog[2, (2 - xh - s1)/2],
    li3 -> PolyLog[2, 1 - xh/CW2],
    li4 -> PolyLog[2, 1 - xh],
    li5 -> PolyLog[2, (1 - s0)/2],
    li6 -> PolyLog[2, (2 - (xh + s2)/CW2)/2],
    li7 -> PolyLog[2, (xh - s2)/(2 CW2)],
    prefEW -> AL^2 MM^2 ZetaL aeps / (MZ^2 Pi^2 CW^4 (1 - CW^2)^2 4608 (-1 + 4*CW^2) (-1 + xh) (4*CW^2 - xh)^2)
}


(* Eq (49), arxiv:1607:06292 *)
amu2LBEW = prefEW (
-80*(-6*li4 + pi2) + (36864*CW14*(6*l32 - 3*lc2 - 12*li5 + pi2))/s0 -
 (6912*CW12*(78*l32 - 39*lc2 - 156*li5 + 13*pi2 + 32*s0 + 16*lc*s0))/s0 +
 32*CW2*(-120 + 120*lh - 66*li4 + 11*pi2 + 60*l1*l2*s1 - 60*li1*s1 -
   60*li2*s1 + 10*pi2*s1) -
 (64*CW10*(-6570*l32 + 3285*lc2 + 13140*li5 - 1095*pi2 - 3744*s0 -
    1668*lc*s0 + 336*lh*s0 - 48*li3*s0 - 576*li4*s0 + 104*pi2*s0 +
    192*l1*l2*s0*s1 - 192*li1*s0*s1 - 192*li2*s0*s1 + 32*pi2*s0*s1))/s0 -
 4*CW4*(-3744 + 4704*lh + 1632*li4 - 272*pi2 + 2592*l1*l2*s1 - 2592*li1*s1 -
   2592*li2*s1 + 432*pi2*s1 - 1386*l4*l5*s2 + 1386*li6*s2 + 1386*li7*s2 -
   231*pi2*s2) + (4*CW6*(3024*l32 - 1512*lc2 - 6048*li5 + 504*pi2 - 3348*s0 -
    2244*lc*s0 + 4356*lh*s0 - 96*li3*s0 + 11136*li4*s0 - 1816*pi2*s0 +
    2304*l1*l2*s0*s1 - 2304*li1*s0*s1 - 2304*li2*s0*s1 + 384*pi2*s0*s1 -
    6222*l4*l5*s0*s2 + 6222*li6*s0*s2 + 6222*li7*s0*s2 - 1037*pi2*s0*s2))/
  s0 + (16*CW8*(-7596*l32 + 3798*lc2 + 15192*li5 - 1266*pi2 + 852*s0 +
    1224*lc*s0 - 564*lh*s0 + 48*li3*s0 - 4416*li4*s0 + 704*pi2*s0 +
    576*l1*l2*s0*s1 - 576*li1*s0*s1 - 576*li2*s0*s1 + 96*pi2*s0*s1 +
    678*l4*l5*s0*s2 - 678*li6*s0*s2 - 678*li7*s0*s2 + 113*pi2*s0*s2))/s0 +
 (-1280*CW4*(-6*li4 + pi2) + 8192*CW6*(-6*li4 + pi2) -
   768*CW10*(-4*li3 + 64*li4 - 10*pi2 + 90*l4*l5*s2 - 90*li6*s2 - 90*li7*s2 +
     15*pi2*s2) + 256*CW8*(336*li4 + 54*(l4*l5 - li6 - li7)*s2 +
     pi2*(-56 + 9*s2)) + 1024*CW12*(-12*li3 + 54*(l4*l5 - li6 - li7)*s2 +
     pi2*(2 + 9*s2)))/xh^2 + (2048*CW12*(108 + 54*lc - 54*lh + 6*li3 - pi2) +
   640*CW2*(-6*li4 + pi2) - 64*CW4*(-120 + 120*lh - 354*li4 + 59*pi2 +
     60*l1*l2*s1 - 60*li1*s1 - 60*li2*s1 + 10*pi2*s1) -
   16*CW8*(-9024 - 1920*lc + 7296*lh - 48*li3 - 192*li4 + 40*pi2 +
     2688*l1*l2*s1 - 2688*li1*s1 - 2688*li2*s1 + 448*pi2*s1 - 6594*l4*l5*s2 +
     6594*li6*s2 + 6594*li7*s2 - 1099*pi2*s2) +
   4*CW6*(-12288 + 12288*lh - 7680*li4 + 1280*pi2 + 6144*l1*l2*s1 -
     6144*li1*s1 - 6144*li2*s1 + 1024*pi2*s1 - 5442*l4*l5*s2 + 5442*li6*s2 +
     5442*li7*s2 - 907*pi2*s2) + 1024*CW10*(-330 - 147*lc + 195*lh - 6*li3 +
     12*li4 - pi2 + 24*l1*l2*s1 - 24*li1*s1 - 24*li2*s1 + 4*pi2*s1 -
     72*l4*l5*s2 + 72*li6*s2 + 72*li7*s2 - 12*pi2*s2))/xh +
 ((-14592*CW12*(6*l32 - 3*lc2 - 12*li5 + pi2))/s0 -
   20*(-24 + 24*lh + 6*li4 - pi2 + 12*l1*l2*s1 - 12*li1*s1 - 12*li2*s1 +
     2*pi2*s1) - (64*CW10*(-3762*l32 + 1881*lc2 + 7524*li5 - 627*pi2 -
      1824*s0 - 684*lc*s0 - 384*lh*s0 - 768*li4*s0 + 768*l1*l2*s0*s1 -
      768*li1*s0*s1 - 768*li2*s0*s1 + 128*pi2*s0*s1))/s0 +
   (64*CW8*(-3114*l32 + 1557*lc2 + 6228*li5 - 519*pi2 - 2853*s0 - 768*lc*s0 -
      291*lh*s0 - 24*li3*s0 - 1632*li4*s0 + 58*pi2*s0 + 1440*l1*l2*s0*s1 -
      1440*li1*s0*s1 - 1440*li2*s0*s1 + 240*pi2*s0*s1))/s0 +
   2*CW2*(864 + 96*lh + 1824*li4 - 304*pi2 + 288*l1*l2*s1 - 288*li1*s1 -
     288*li2*s1 + 48*pi2*s1 + 270*l4*l5*s2 - 270*li6*s2 - 270*li7*s2 +
     45*pi2*s2) + (4*CW4*(-1512*l32 + 756*lc2 + 3024*li5 - 252*pi2 -
      5190*s0 - 345*lc*s0 + 3081*lh*s0 - 804*li3*s0 - 6576*li4*s0 +
      818*pi2*s0 + 2496*l1*l2*s0*s1 - 2496*li1*s0*s1 - 2496*li2*s0*s1 +
      416*pi2*s0*s1 - 198*l4*l5*s0*s2 + 198*li6*s0*s2 + 198*li7*s0*s2 -
      33*pi2*s0*s2))/s0 - (16*CW6*(-3690*l32 + 1845*lc2 + 7380*li5 -
      615*pi2 - 3939*s0 - 780*lc*s0 + 1158*lh*s0 - 828*li3*s0 - 4848*li4*s0 +
      348*pi2*s0 + 3360*l1*l2*s0*s1 - 3360*li1*s0*s1 - 3360*li2*s0*s1 +
      560*pi2*s0*s1 + 342*l4*l5*s0*s2 - 342*li6*s0*s2 - 342*li7*s0*s2 +
      57*pi2*s0*s2))/s0)*xh +
 ((384*CW10*(6*l32 - 3*lc2 - 12*li5 + pi2 - 48*s0 - 96*lh*s0 - 96*lh2*s0 -
      512*li4*s0 - 32*pi2*s0 + 144*l1*l2*s0*s1 - 144*li1*s0*s1 -
      144*li2*s0*s1 + 24*pi2*s0*s1))/s0 +
   (12*CW6*(1734*l32 - 867*lc2 - 3468*li5 + 289*pi2 + 1368*s0 - 232*lc*s0 -
      964*lh*s0 - 2688*lh2*s0 - 1072*li3*s0 - 10688*li4*s0 - 936*pi2*s0 +
      384*l1*l2*s0*s1 - 384*li1*s0*s1 - 384*li2*s0*s1 + 64*pi2*s0*s1))/s0 -
   (16*CW8*(1206*l32 - 603*lc2 - 2412*li5 + 201*pi2 - 960*s0 + 72*lc*s0 -
      3264*lh*s0 - 4032*lh2*s0 - 19968*li4*s0 - 1344*pi2*s0 +
      4512*l1*l2*s0*s1 - 4512*li1*s0*s1 - 4512*li2*s0*s1 + 752*pi2*s0*s1))/
    s0 + CW2*(3867 + 426*lc - 2106*lh + 1572*li3 + 5568*li4 - 384*pi2 +
     (756*l32)/s0 - (378*lc2)/s0 - (1512*li5)/s0 + (126*pi2)/s0 -
     4032*l1*l2*s1 + 4032*li1*s1 + 4032*li2*s1 - 672*pi2*s1 - 420*l4*l5*s2 +
     420*li6*s2 + 420*li7*s2 - 70*pi2*s2) +
   4*(-150 + 90*lh - 90*li4 + 15*pi2 + 30*l1*l2*s1 - 30*li1*s1 - 30*li2*s1 +
     5*pi2*s1 - 48*l4*l5*s2 + 48*li6*s2 + 48*li7*s2 - 8*pi2*s2) +
   (6*CW4*(-1122*l32 + 561*lc2 + 2244*li5 - 187*pi2 - 1774*s0 - 48*lc*s0 +
      478*lh*s0 + 768*lh2*s0 - 512*li3*s0 - 224*li4*s0 + 372*pi2*s0 +
      2784*l1*l2*s0*s1 - 2784*li1*s0*s1 - 2784*li2*s0*s1 + 464*pi2*s0*s1 +
      792*l4*l5*s0*s2 - 792*li6*s0*s2 - 792*li7*s0*s2 + 132*pi2*s0*s2))/s0)*
  xh2 + (-3072*CW10*(-15*lh2 - 60*li4 - 5*pi2 + 6*l1*l2*s1 - 6*li1*s1 -
     6*li2*s1 + pi2*s1) + (2*CW4*(342*l32 - 171*lc2 - 684*li5 + 57*pi2 +
      3366*s0 + 834*lc*s0 + 6444*lh*s0 + 5184*lh2*s0 + 3144*li3*s0 +
      29184*li4*s0 + 1728*pi2*s0 - 8256*l1*l2*s0*s1 + 8256*li1*s0*s1 +
      8256*li2*s0*s1 - 1376*pi2*s0*s1))/s0 +
   (48*CW8*(30*l32 - 15*lc2 - 60*li5 + 5*pi2 + 192*s0 + 384*lh*s0 -
      1296*lh2*s0 - 4672*li4*s0 - 432*pi2*s0 + 96*l1*l2*s0*s1 -
      96*li1*s0*s1 - 96*li2*s0*s1 + 16*pi2*s0*s1))/s0 +
   (4*CW6*(-450*l32 + 225*lc2 + 900*li5 - 75*pi2 - 3936*s0 - 180*lc*s0 -
      7680*lh*s0 + 2016*lh2*s0 - 1920*li4*s0 + 672*pi2*s0 +
      7296*l1*l2*s0*s1 - 7296*li1*s0*s1 - 7296*li2*s0*s1 + 1216*pi2*s0*s1))/
    s0 + 4*(-6 - 15*lh - 48*li3 - 102*li4 + 102*l1*l2*s1 - 102*li1*s1 -
     102*li2*s1 + 17*pi2*s1 + 48*l4*l5*s2 - 48*li6*s2 - 48*li7*s2 +
     8*pi2*s2) - (CW2*(108*l32 - 54*lc2 - 216*li5 + 18*pi2 + 747*s0 +
      426*lc*s0 + 1350*lh*s0 + 2304*lh2*s0 + 804*li3*s0 + 9696*li4*s0 +
      768*pi2*s0 - 672*l1*l2*s0*s1 + 672*li1*s0*s1 + 672*li2*s0*s1 -
      112*pi2*s0*s1 + 768*l4*l5*s0*s2 - 768*li6*s0*s2 - 768*li7*s0*s2 +
      128*pi2*s0*s2))/s0)*xh3 + (-3072*CW10*(3*lh2 + 12*li4 + pi2) +
   768*CW8*(-9*lh2 - 36*li4 - 3*pi2 + 12*l1*l2*s1 - 12*li1*s1 - 12*li2*s1 +
     2*pi2*s1) - 24*(-6 - 12*lh - 12*lh2 - 8*li3 - 65*li4 - 4*pi2 +
     18*l1*l2*s1 - 18*li1*s1 - 18*li2*s1 + 3*pi2*s1) +
   48*CW4*(42 + 84*lh - 312*lh2 - 1136*li4 - 104*pi2 + 42*l1*l2*s1 -
     42*li1*s1 - 42*li2*s1 + 7*pi2*s1) -
   192*CW6*(6 + 12*lh - 156*lh2 - 608*li4 - 52*pi2 + 66*l1*l2*s1 -
     66*li1*s1 - 66*li2*s1 + 11*pi2*s1) +
   24*CW2*(-42 - 84*lh + 36*lh2 - 32*li3 + 28*li4 + 12*pi2 + 78*l1*l2*s1 -
     78*li1*s1 - 78*li2*s1 + 13*pi2*s1))*xh4 +
 (1536*CW8*(3*lh2 + 12*li4 + pi2) + 24*(-15*lh2 - 60*li4 - 5*pi2 +
     6*l1*l2*s1 - 6*li1*s1 - 6*li2*s1 + pi2*s1) +
   336*CW4*(-3*lh2 - 12*li4 - pi2 + 6*l1*l2*s1 - 6*li1*s1 - 6*li2*s1 +
     pi2*s1) - 192*CW6*(27*lh2 + 108*li4 + 9*pi2 + 6*l1*l2*s1 - 6*li1*s1 -
     6*li2*s1 + pi2*s1) - 24*CW2*(-81*lh2 - 324*li4 - 27*pi2 + 42*l1*l2*s1 -
     42*li1*s1 - 42*li2*s1 + 7*pi2*s1))*xh5 +
 (24*(3*lh2 + 12*li4 + pi2) - 168*CW2*(3*lh2 + 12*li4 + pi2) +
   336*CW4*(3*lh2 + 12*li4 + pi2) - 192*CW6*(3*lh2 + 12*li4 + pi2))*xh6
) //. expandEW


Null
