(* Loop functions for (g-2)_mu *)

Li2[x_] := PolyLog[2, x]

(* arxiv:1003.5820, Eq.(13) *)
F1C[0] := 4

F1C[1] := 1

F1C[Infinity] := 0

F1C[x_] := (2/(1 - x)^4)*(2 + 3*x - 6*x^2 + x^3 + 6*x*Log[x])

(* arxiv:1003.5820, Eq.(14) *)
F2C[1] := 1

F2C[Infinity] := 0

F2C[x_] := (3/(2*(1 - x)^3))*(-3 + 4*x - x^2 - 2*Log[x])

(* arxiv:1003.5820, Eq.(37) *)
F3C[1] := 1

F3C[Infinity] := 0

F3C[x_] := (4/(141*(1 - x)^4))*(
    (1 - x)*(151*x^2 - 335*x + 592)
    + 6*(21*x^3 - 108*x^2 - 93*x + 50)*Log[x]
    - 54*x*(x^2 - 2*x - 2)*Log[x]^2
    - 108*x*(x^2 - 2*x + 12)*Li2[1 - x]
)

(* arxiv:1003.5820, Eq.(38) *)
F4C[1] := 1

F4C[Infinity] := 0

F4C[x_] := (-9/(122*(1 - x)^3))*(
    8*(x^2 - 3*x + 2)
    + (11*x^2 - 40*x + 5)*Log[x]
    - 2*(x^2 - 2*x - 2)*Log[x]^2
    - 4*(x^2 - 2*x + 9)*Li2[1 - x]
)

(* arxiv:1003.5820, Eq.(15) *)
F1N[0] := 2

F1N[1] := 1

F1N[Infinity] := 0

F1N[x_] := (2/(1 - x)^4)*(1 - 6*x + 3*x^2 + 2*x^3 - 6*x^2*Log[x])

(* arxiv:1003.5820, Eq.(16) *)
F2N[0] := 3

F2N[1] := 1

F2N[Infinity] := 0

F2N[x_] := (3/(1 - x)^3)*(1 - x^2 + 2*x*Log[x])

(* arxiv:1003.5820, Eq.(39) *)
F3N[0] := 8/105

F3N[1] := 1

F3N[Infinity] := 0

F3N[x_] := (4/(105*(1 - x)^4))*(
    (1 - x)*(-97*x^2 - 529*x + 2)
    + 6*x^2*(13*x + 81)*Log[x]
    + 108*x*(7*x + 4)*Li2[1 - x]
)

(* arxiv:1003.5820, Eq.(40) *)
F4N[0] := (-3*(-9 + Pi^2))/4

F4N[1] := 1

F4N[Infinity] := 0

F4N[x_] := (-(9/4/(1 - x)^3))*(
    (x + 3)*(x*Log[x] + x - 1)
    + (6*x + 2)*Li2[1 - x]
)

(* arxiv:1311.1775, Eq.(6.3a) *)
Fa[1, 1] := 1/4

Fa[x_, x_] := (2 + 3 x - 6 x^2 + x^3 + 6 x Log[x])/(2 (-1 + x)^4 x)

Fa[x_, y_] := -(G3[x] - G3[y])/(x - y)

(* arxiv:1311.1775, Eq.(6.3b) *)
Fb[1, 1] := 1/12

Fb[x_, x_] := (-5 + 4 x + x^2 - 2 Log[x] - 4 x Log[x])/(2 (-1 + x)^4)

Fb[x_, y_] := -(G4[x] - G4[y])/(x - y)

(* arxiv:1311.1775, Eq.(6.4a) *)
G3[1] := 1/3

G3[x_] := 1/(2 (x - 1)^3) ((x - 1)(x - 3) + 2 Log[x])

(* arxiv:1311.1775, Eq.(6.4b) *)
G4[1] := 1/6

G4[x_] := 1/(2 (x - 1)^3) ((x - 1)(x + 1) - 2 x Log[x])

(* I[a,b,c] with squared arguments *)
I2abc[a_, a_, a_] := 1/(2 a)

I2abc[a_, a_, c_] := (a - c - c*Log[a/c])/(a - c)^2

I2abc[a_, c_, a_] := I2abc[a, a, c]

I2abc[c_, a_, a_] := I2abc[a, a, c]

I2abc[a_, b_, c_] :=
    (a b Log[a/b] + b c Log[b/c] + c a Log[c/a])/((a - b) (b - c) (a - c))

(* I[a,b,c] with non-squared arguments *)
Iabc[a_, b_, c_] := I2abc[a^2, b^2, c^2]

(* arxiv:hep-ph/0609168, Eq.(70) *)
fPS[0] := 0

fPS[1/4] := 2 Log[2]

fPS[z_] := Re[Module[{y = Sqrt[1 - 4z]}, 2z/y (PolyLog[2, 1 - (1-y)/(2z)] - PolyLog[2, 1 - (1+y)/(2z)])]]

(* arxiv:hep-ph/0609168, Eq.(71) *)
fS[0] := 0

fS[z_] := (2z - 1) fPS[z] - 2z(2 + Log[z])

(* arxiv:hep-ph/0609168, Eq.(72) *)
fsferm[0] := 0

fsferm[z_] := z/2 (2 + Log[z] - fPS[z])

(* arxiv:1607.06292, Eq.(60), with extra global prefactor factor z *)
fCl[0] := 0

fCl[z_] := z (z + z (z - 1) (PolyLog[2, 1 - 1/z] - Pi^2/6) + (z - 1/2) Log[z]);

(* loop function for fermionic 2-loop Barr-Zee diagram with Z boson and pseudoscalar mediator *)
FPZ[1/4, 1/4] := (-1 - 2*Log[2])/3

FPZ[x_, x_] := -2 x (fPS[x] + Log[x])/(-1 + 4*x)

FPZ[x_, y_] := (y fPS[x] - x fPS[y])/(x - y)

(* loop function for fermionic 2-loop Barr-Zee diagram with Z boson and scalar mediator *)
FSZ[1/4, 1/4] := (-1 + Log[16])/3

FSZ[x_, x_] := (2 x (1 - 4 x + 2 x fPS[x] + Log[x] - 2 x Log[x]))/(-1 + 4 x)

FSZ[x_, y_] := (y fS[x] - x fS[y])/(x - y)

(* loop function for leptonic 2-loop Barr-Zee diagram with W boson and scalar mediator *)
FCWl[x_, x_] := (-3*x + 12*x^2 + Pi^2*x^2 - 2*Pi^2*x^3 - 6*x^2*Log[1 - (-1 + x)/x] +
  6*x^2*Log[x] + 6*x^2*PolyLog[2, 1 - 1/x] - 6*x^3*PolyLog[2, 1 - 1/x] -
  12*x^2*PolyLog[2, (-1 + x)/x] + 18*x^3*PolyLog[2, (-1 + x)/x])/6

FCWl[x_, y_] := (y fCl[x] - x fCl[y])/(x - y)

(* arxiv:1502.04199, Eq.(25) *)
(* Module[{x}, w/2 Integrate[(2x(1-x)-1)/(w-x(1-x)) Log[w/(x(1-x))], {x,0,1}]] *)
F1[0] := 0

F1[1/4] := -1/2

F1[w_] := Re[-(w*(12*Sqrt[1 - 4*w] + Log[4096]*Log[-1 + Sqrt[1 - 4*w]] -
    w*Log[16777216]*Log[-1 + Sqrt[1 - 4*w]] + 6*Log[(1 - Sqrt[1 - 4*w])^(-1)]*
     Log[-1 + Sqrt[1 - 4*w]] - 12*w*Log[(1 - Sqrt[1 - 4*w])^(-1)]*
     Log[-1 + Sqrt[1 - 4*w]] - Log[4096]*Log[1 + Sqrt[1 - 4*w]] +
    w*Log[16777216]*Log[1 + Sqrt[1 - 4*w]] - 6*Log[(1 - Sqrt[1 - 4*w])^(-1)]*
     Log[1 + Sqrt[1 - 4*w]] + 12*w*Log[(1 - Sqrt[1 - 4*w])^(-1)]*
     Log[1 + Sqrt[1 - 4*w]] - 6*Log[-1 + Sqrt[1 - 4*w]]*
     Log[1 + Sqrt[1 - 4*w]] + 12*w*Log[-1 + Sqrt[1 - 4*w]]*
     Log[1 + Sqrt[1 - 4*w]] + 6*Log[1 + Sqrt[1 - 4*w]]^2 -
    12*w*Log[1 + Sqrt[1 - 4*w]]^2 + 6*Sqrt[1 - 4*w]*Log[w] +
    6*Log[-1 + Sqrt[1 - 4*w]]*Log[w] - 12*w*Log[-1 + Sqrt[1 - 4*w]]*Log[w] -
    6*Log[1 + Sqrt[1 - 4*w]]*Log[w] + 12*w*Log[1 + Sqrt[1 - 4*w]]*Log[w] +
    (6 - 12*w)*PolyLog[2, (-1 + Sqrt[1 - 4*w])/(1 + Sqrt[1 - 4*w])] +
    6*(-1 + 2*w)*PolyLog[2, (1 + Sqrt[1 - 4*w])/(-1 + Sqrt[1 - 4*w])]))/
    (6*Sqrt[1 - 4*w])]

(* arxiv:1502.04199, Eq.(26) *)
(* Module[{x}, w/2 Integrate[1/(w-x(1-x)) Log[w/(x(1-x))], {x,0,1}]] *)
F1t[0] := 0

F1t[1/4] := Log[2]

F1t[w_] := Re[-((w*((Log[-1 + Sqrt[1 - 4*w]] - Log[1 + Sqrt[1 - 4*w]])*
     (Log[1 + Sqrt[1 - 4*w]] - Log[(1 + Sqrt[1 - 4*w])/w] - Log[w]) +
    PolyLog[2, -(1 + Sqrt[1 - 4*w] - 2*w)/(2*w)] -
    PolyLog[2, (-1 + Sqrt[1 - 4*w] + 2*w)/(2*w)]))/Sqrt[1 - 4*w])]

(* arxiv:1502.04199, Eq.(27) *)
(* Module[{x}, 1/2 Integrate[x(x-1)/(w-x(1-x)) Log[w/(x(1-x))], {x,0,1}]] *)
F2[1/4] := 1 - Log[4]

F2[w_] := Re[(2*Sqrt[1 - 4*w] - 2*w*Log[4]*Log[-1 + Sqrt[1 - 4*w]] -
  2*w*Log[(1 - Sqrt[1 - 4*w])^(-1)]*Log[-1 + Sqrt[1 - 4*w]] +
  w*Log[16]*Log[1 + Sqrt[1 - 4*w]] + 2*w*Log[(1 - Sqrt[1 - 4*w])^(-1)]*
   Log[1 + Sqrt[1 - 4*w]] + 2*w*Log[-1 + Sqrt[1 - 4*w]]*
   Log[1 + Sqrt[1 - 4*w]] - 2*w*Log[1 + Sqrt[1 - 4*w]]^2 +
  Sqrt[1 - 4*w]*Log[w] - 2*w*Log[-1 + Sqrt[1 - 4*w]]*Log[w] +
  2*w*Log[1 + Sqrt[1 - 4*w]]*Log[w] -
  2*w*PolyLog[2, (-1 + Sqrt[1 - 4*w])/(1 + Sqrt[1 - 4*w])] +
  2*w*PolyLog[2, (1 + Sqrt[1 - 4*w])/(-1 + Sqrt[1 - 4*w])])/(2*Sqrt[1 - 4*w])]

(* arxiv:1502.04199, Eq.(28) *)
(* Module[{x}, 1/2 Integrate[(x w(3x(4x-1)+10)-x(1-x))/(w-x(1-x)) Log[w/(x(1-x))], {x,0,1}]] *)
(* F3[z_] := (17/2 - 15 z) fPS[z]/2 + (1/2 + 75/10 z) (2 + Log[z]) *)

F3[1/4] := 19/4

F3[w_] := Re[(4*Sqrt[-1 + 4*w] + 60*w*Sqrt[-1 + 4*w] + 2*Sqrt[-1 + 4*w]*Log[w] +
  30*w*Sqrt[-1 + 4*w]*Log[w] + (38*I)*w*Log[2]*Log[-1 - I*Sqrt[-1 + 4*w]] -
  (24*I)*w^2*Log[2]*Log[-1 - I*Sqrt[-1 + 4*w]] -
  38*w*Sqrt[-1 + 4*w]*Log[2]*Log[-1 - I*Sqrt[-1 + 4*w]] +
  24*w^2*Sqrt[-1 + 4*w]*Log[2]*Log[-1 - I*Sqrt[-1 + 4*w]] +
  (19*I)*w*Log[w]*Log[-1 - I*Sqrt[-1 + 4*w]] -
  (12*I)*w^2*Log[w]*Log[-1 - I*Sqrt[-1 + 4*w]] -
  19*w*Sqrt[-1 + 4*w]*Log[w]*Log[-1 - I*Sqrt[-1 + 4*w]] +
  12*w^2*Sqrt[-1 + 4*w]*Log[w]*Log[-1 - I*Sqrt[-1 + 4*w]] -
  (19*I)*w*Log[2]*Log[1 - I*Sqrt[-1 + 4*w]] +
  (12*I)*w^2*Log[2]*Log[1 - I*Sqrt[-1 + 4*w]] +
  19*w*Sqrt[-1 + 4*w]*Log[2]*Log[1 - I*Sqrt[-1 + 4*w]] -
  12*w^2*Sqrt[-1 + 4*w]*Log[2]*Log[1 - I*Sqrt[-1 + 4*w]] -
  (19*I)*w*Log[w]*Log[1 - I*Sqrt[-1 + 4*w]] +
  (12*I)*w^2*Log[w]*Log[1 - I*Sqrt[-1 + 4*w]] +
  19*w*Sqrt[-1 + 4*w]*Log[w]*Log[1 - I*Sqrt[-1 + 4*w]] -
  12*w^2*Sqrt[-1 + 4*w]*Log[w]*Log[1 - I*Sqrt[-1 + 4*w]] -
  (15*I)*w*Log[2]*Log[-1 + I*Sqrt[-1 + 4*w]] +
  (48*I)*w^2*Log[2]*Log[-1 + I*Sqrt[-1 + 4*w]] -
  19*w*Sqrt[-1 + 4*w]*Log[2]*Log[-1 + I*Sqrt[-1 + 4*w]] +
  12*w^2*Sqrt[-1 + 4*w]*Log[2]*Log[-1 + I*Sqrt[-1 + 4*w]] -
  (15*I)*w*Log[w]*Log[-1 + I*Sqrt[-1 + 4*w]] +
  (48*I)*w^2*Log[w]*Log[-1 + I*Sqrt[-1 + 4*w]] -
  19*w*Sqrt[-1 + 4*w]*Log[w]*Log[-1 + I*Sqrt[-1 + 4*w]] +
  12*w^2*Sqrt[-1 + 4*w]*Log[w]*Log[-1 + I*Sqrt[-1 + 4*w]] +
  (13*I)*w*Log[2]*Log[1 + I*Sqrt[-1 + 4*w]] -
  (66*I)*w^2*Log[2]*Log[1 + I*Sqrt[-1 + 4*w]] +
  19*w*Sqrt[-1 + 4*w]*Log[2]*Log[1 + I*Sqrt[-1 + 4*w]] -
  12*w^2*Sqrt[-1 + 4*w]*Log[2]*Log[1 + I*Sqrt[-1 + 4*w]] +
  (15*I)*w*Log[w]*Log[1 + I*Sqrt[-1 + 4*w]] -
  (48*I)*w^2*Log[w]*Log[1 + I*Sqrt[-1 + 4*w]] +
  19*w*Sqrt[-1 + 4*w]*Log[w]*Log[1 + I*Sqrt[-1 + 4*w]] -
  12*w^2*Sqrt[-1 + 4*w]*Log[w]*Log[1 + I*Sqrt[-1 + 4*w]] -
  (19*I)*w*Log[-1 - I*Sqrt[-1 + 4*w]]*Log[1 + I*Sqrt[-1 + 4*w]] +
  (12*I)*w^2*Log[-1 - I*Sqrt[-1 + 4*w]]*Log[1 + I*Sqrt[-1 + 4*w]] +
  19*w*Sqrt[-1 + 4*w]*Log[-1 - I*Sqrt[-1 + 4*w]]*Log[1 + I*Sqrt[-1 + 4*w]] -
  12*w^2*Sqrt[-1 + 4*w]*Log[-1 - I*Sqrt[-1 + 4*w]]*
   Log[1 + I*Sqrt[-1 + 4*w]] + (19*I)*w*Log[1 - I*Sqrt[-1 + 4*w]]*
   Log[1 + I*Sqrt[-1 + 4*w]] - (12*I)*w^2*Log[1 - I*Sqrt[-1 + 4*w]]*
   Log[1 + I*Sqrt[-1 + 4*w]] - 19*w*Sqrt[-1 + 4*w]*Log[1 - I*Sqrt[-1 + 4*w]]*
   Log[1 + I*Sqrt[-1 + 4*w]] + 12*w^2*Sqrt[-1 + 4*w]*
   Log[1 - I*Sqrt[-1 + 4*w]]*Log[1 + I*Sqrt[-1 + 4*w]] +
  (15*I)*w*Log[-1 + I*Sqrt[-1 + 4*w]]*Log[1 + I*Sqrt[-1 + 4*w]] -
  (48*I)*w^2*Log[-1 + I*Sqrt[-1 + 4*w]]*Log[1 + I*Sqrt[-1 + 4*w]] +
  19*w*Sqrt[-1 + 4*w]*Log[-1 + I*Sqrt[-1 + 4*w]]*Log[1 + I*Sqrt[-1 + 4*w]] -
  12*w^2*Sqrt[-1 + 4*w]*Log[-1 + I*Sqrt[-1 + 4*w]]*
   Log[1 + I*Sqrt[-1 + 4*w]] - (15*I)*w*Log[1 + I*Sqrt[-1 + 4*w]]^2 +
  (48*I)*w^2*Log[1 + I*Sqrt[-1 + 4*w]]^2 - 19*w*Sqrt[-1 + 4*w]*
   Log[1 + I*Sqrt[-1 + 4*w]]^2 + 12*w^2*Sqrt[-1 + 4*w]*
   Log[1 + I*Sqrt[-1 + 4*w]]^2 + (19*I)*w*Log[-1 - I*Sqrt[-1 + 4*w]]*
   Log[I/(I + Sqrt[-1 + 4*w])] - (12*I)*w^2*Log[-1 - I*Sqrt[-1 + 4*w]]*
   Log[I/(I + Sqrt[-1 + 4*w])] - 19*w*Sqrt[-1 + 4*w]*
   Log[-1 - I*Sqrt[-1 + 4*w]]*Log[I/(I + Sqrt[-1 + 4*w])] +
  12*w^2*Sqrt[-1 + 4*w]*Log[-1 - I*Sqrt[-1 + 4*w]]*
   Log[I/(I + Sqrt[-1 + 4*w])] - (2*I)*w*Log[1 + I*Sqrt[-1 + 4*w]]*
   Log[I/(I + Sqrt[-1 + 4*w])] - (18*I)*w^2*Log[1 + I*Sqrt[-1 + 4*w]]*
   Log[I/(I + Sqrt[-1 + 4*w])] - (19*I)*w*Log[1 - I*Sqrt[-1 + 4*w]]*
   Log[(2*I)/(I + Sqrt[-1 + 4*w])] + (12*I)*w^2*Log[1 - I*Sqrt[-1 + 4*w]]*
   Log[(2*I)/(I + Sqrt[-1 + 4*w])] + 19*w*Sqrt[-1 + 4*w]*
   Log[1 - I*Sqrt[-1 + 4*w]]*Log[(2*I)/(I + Sqrt[-1 + 4*w])] -
  12*w^2*Sqrt[-1 + 4*w]*Log[1 - I*Sqrt[-1 + 4*w]]*
   Log[(2*I)/(I + Sqrt[-1 + 4*w])] - (15*I)*w*Log[-1 + I*Sqrt[-1 + 4*w]]*
   Log[(2*I)/(I + Sqrt[-1 + 4*w])] + (48*I)*w^2*Log[-1 + I*Sqrt[-1 + 4*w]]*
   Log[(2*I)/(I + Sqrt[-1 + 4*w])] - 19*w*Sqrt[-1 + 4*w]*
   Log[-1 + I*Sqrt[-1 + 4*w]]*Log[(2*I)/(I + Sqrt[-1 + 4*w])] +
  12*w^2*Sqrt[-1 + 4*w]*Log[-1 + I*Sqrt[-1 + 4*w]]*
   Log[(2*I)/(I + Sqrt[-1 + 4*w])] + (17*I)*w*Log[1 + I*Sqrt[-1 + 4*w]]*
   Log[(2*I)/(I + Sqrt[-1 + 4*w])] - (30*I)*w^2*Log[1 + I*Sqrt[-1 + 4*w]]*
   Log[(2*I)/(I + Sqrt[-1 + 4*w])] + 19*w*Sqrt[-1 + 4*w]*
   Log[1 + I*Sqrt[-1 + 4*w]]*Log[(2*I)/(I + Sqrt[-1 + 4*w])] -
  12*w^2*Sqrt[-1 + 4*w]*Log[1 + I*Sqrt[-1 + 4*w]]*
   Log[(2*I)/(I + Sqrt[-1 + 4*w])] - (2*I)*w*(-17 + 30*w)*
   PolyLog[2, (-I + Sqrt[-1 + 4*w])/(I + Sqrt[-1 + 4*w])] +
  (2*I)*w*(-17 + 30*w)*PolyLog[2, (I + Sqrt[-1 + 4*w])/
  (-I + Sqrt[-1 + 4*w])])/(4*Sqrt[-1 + 4*w])]

fCSd[0, xd_, qu_, qd_] :=
    -1/12*(xd*(-12*xd + 2*Pi^2*qd*xd + 2*Pi^2*xd^2 - 3*(qd + qu + 4*xd)*Log[xd] + 6*xd*(qd + xd)*Log[xd]^2 - 3*qd*Log[xu] - 3*qu*Log[xu] + 12*xd*(qd + xd)*PolyLog[2, 1 - xd]))

fCSd[xu_, 0, qu_, qd_] := 0

fCSd[1/4, 1/4, qu_, qd_] := -1/2*((qd + qu)*Log[2])

fCSd[xu_, xu_, qu_, qd_] :=
    ((qd + qu)*xu*(Log[xu] + (Sqrt[1 - 4*xu]*xu*(Pi^2 + 6*Log[(1 - Sqrt[1 - 4*xu])/2]^2 - 3*Log[xu]^2 - 12*PolyLog[2, (1 - Sqrt[1 - 4*xu])/2]))/(-3 + 12*xu)))/2

(* calculate Phi[xd, xu, 1]/y with y = (xu - xd)^2 - 2*(xu + xd) + 1, properly handle the case y = 0 *)
PhiOverY[xu_, xd_] :=
    Which[PossibleZeroQ[xu - (1 - 2 Sqrt[xd] + xd)],
          -Log[-1 + Sqrt[xd]]/Sqrt[xd] + Log[xd]/(2*(-1 + Sqrt[xd])),
          PossibleZeroQ[xu - (1 + 2 Sqrt[xd] + xd)],
          Log[1 + Sqrt[xd]]/Sqrt[xd] - Log[xd]/(2*(1 + Sqrt[xd])),
          True,
          Phi[xd, xu, 1]/((xu - xd)^2 - 2*(xu + xd) + 1)
    ]

fCSd[xu_, xd_, qu_, qd_] :=
    Module[{s, c, cbar, lxu, lxd},
           s = 1/4*(qu + qd);
           c = (xu - xd)^2 - qu*xu + qd*xd;
           cbar = (xu - qu)*xu - (xd + qd)*xd;
           lxu = Log[xu];
           lxd = Log[xd];
           xd*(-(xu - xd) + (cbar - c*(xu - xd))*PhiOverY[xu, xd] + c*(Li2[1 - xd/xu] - 1/2*lxu*(lxd - lxu)) + (s + xd)*lxd + (s - xu)*lxu)
    ]

fCSu[xu_, xd_, qu_, qd_] :=
    Module[{lxu = Log[xu], lxd = Log[xd]},
           xu*(fCSd[xu, xd, qu + 2, qd + 2]/xd - 4/3*(xu - xd - 1)*PhiOverY[xu, xd] - 1/3*(lxd + lxu)*(lxd - lxu))
    ]

FCWd[xu_, xd_, yu_, yd_, qu_, qd_] :=
    (yd fCSd[xu, xd, qu, qd] - xd fCSd[yu, yd, qu, qd])/(xd - yd)

FCWu[xu_, xd_, yu_, yd_, qu_, qd_] :=
    (yu fCSu[xu, xd, qu, qd] - xu fCSu[yu, yd, qu, qd])/(xu - yu)

(* arxiv:1502.04199, Eq.(29) *)
G[wa_, wb_, x_] := Log[(wa*x+wb*(1-x))/(x(1-x))]/(x(1-x)-wa*x-wb*(1-x))

(* Integral of x^n G[wa,wb,x] over {x,0,1} *)
Gn[wa_, wb_, n_] := NIntegrate[x^n G[wa, wb, x], {x,0,1}]

(* Phi(x,y,z) from Davydychev and Tausk, Nucl. Phys. B397 (1993) 23 *)
(* Note: The parameters x,y,z are interpreted as squared masses *)
(* Formulation as in arxiv:1607.06292, Eq.(68) *)
Phi[x_, y_, z_] := LambdaK[x,y,z]/2 (
    + 2 Log[alphaPlus[x,y,z]] Log[alphaMinus[x,y,z]]
    - Log[x/z] Log[y/z]
    - 2 PolyLog[2, alphaPlus[x,y,z]]
    - 2 PolyLog[2, alphaMinus[x,y,z]]
    + Pi^2/3
)

(* arxiv:1607.06292, Eq.(69) *)
LambdaK[x_, y_, z_] := Sqrt[x^2 + y^2 + z^2 - 2 x y - 2 y z - 2 z x]

(* arxiv:1607.06292, Eq.(70) *)
alphaPlus[x_, y_, z_] := (z + x - y - LambdaK[x,y,z]) / (2 z)

(* arxiv:1607.06292, Eq.(70) *)
alphaMinus[x_, y_, z_] := (z - x + y - LambdaK[x,y,z]) / (2 z)
