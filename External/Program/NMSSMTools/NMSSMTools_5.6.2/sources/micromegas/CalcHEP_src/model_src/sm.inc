% Standard Model - unitary and t'Hooft-Feynman gauges.


external_func(alphaQCD,1).
external_func(McEff,1).
external_func(MbEff,1).
external_func(MtEff,1).
external_func(MqEff,2). 
external_func(initQCD5,4).


option ReduceGamma5=0.
let g5=gamma5.
use sm_tex.

parameter  EE  = 0.31333 : 'Electromagnetic coupling constant (<->1/128)',
	   GG  = 1.117   : 'Strong coupling constant (Z point)  (PDG-94)',
	   SW  = 0.4740  : 'sin of the Weinberg angle 0.474 - "on-shell",481 - "MS-bar" )'.


parameter  CW  = sqrt(1-SW**2) : 'cos of the Weinberg angle'.

% Masses 

parameter Q=100: 'Scale of effective running masses',
            MW=80.385: 'W boson mass',
            GF=EE**2/(2*SW*MW)**2/Sqrt2 :' experimental value 1.166E-5 [1/GeV^2]',  
            MZ=MW/CW:  'Z boson mass',
            Mtp  =172.5:'Top quark pole mass',
            McMc =1.23:     'Mc(Mc)  MS-BAR',
            MbMb =4.25:     'Mb(Mb)  MS-BAR',
            alphaSMZ=0.1184:'Srtong alpha(MZ)', 
            LamQCD=initQCD5(alphaSMZ,McMc,MbMb,Mtp),
            Mb=MbEff(Q),
            Mc=McEff(Q),
            Ms=MqEff(0.096,Q): 's-quark effective mass  via 2MeV running one' .  


do_if gauge_fixing==Feynman.

vector  
	A/A: (photon, gauge),
	Z/Z:('Z boson', mass MZ, width wZ = auto, gauge),
	G/G: (gluon, color c8, gauge),
        'W+'/'W-': ('W boson', mass MW, width wW = auto, gauge).

do_else.

vector  
	A/A: (photon, gauge),
	Z/Z:('Z boson', mass MZ, width wZ = auto),
	G/G: (gluon, color c8, gauge),
        'W+'/'W-': ('W boson', mass MW, width wW = auto).

end_if.

spinor 		ne/Ne:(neutrino,left),        e:(electron),
		nm/Nm:('mu-neutrino',left),   m:(muon),
		nl/Nl:('tau-neutrino',left),  l:('tau-lepton', mass Ml  = 1.777).

spinor		u:('u-quark',color c3),
		d:('d-quark',color c3),
		c:('c-quark',color c3, mass Mc),
		s:('s-quark',color c3, mass Ms),
		t:('t-quark',color c3, mass Mtp, width wt = auto),
		b:('b-quark',color c3, mass Mb).


do_if CKMdim==1.
 let  Vub=0, Vcb=0, Vtd=0, Vts=0, Vtb=1, Vud=1, Vus=0, Vcs=1, Vcd=0.

do_else.

parameter  s12 = 0.221   : 'Parameter of C-K-M matrix (PDG-94)',
           c12  = sqrt(1-s12**2).

do_if CKMdim==2.

  parameter  Vud = c12  ,
             Vus = s12  ,
             Vcs = Vud  ,
             Vcd = -Vus .
  let        Vub = 0, Vcb = 0, Vtd = 0, Vts = 0, Vtb = 1.
  OrthMatrix({{Vud,Vus}, {Vcd,Vcs}}).

do_else.
   
  parameter  s23 = 0.040   : 'Parameter of C-K-M matrix (PDG-94)',
             s13 = 0.0035  : 'Parameter of C-K-M matrix (PDG-94)',
             c23  = sqrt(1-s23**2),
             c13  = sqrt(1-s13**2).
		  
parameter  Vud = c12*c13                ,
           Vus = s12*c13                ,
           Vub = s13     	        ,
           Vcd = (-s12*c23-c12*s23*s13) ,
           Vcs = (c12*c23-s12*s23*s13)  ,
           Vcb = s23*c13                ,
           Vtd = (s12*s23-c12*c23*s13) 	,
           Vts = (-c12*s23-s12*c23*s13)	,
           Vtb = c23*c13  	        .

OrthMatrix( { {Vud,Vus,Vub}, {Vcd,Vcs,Vcb}, {Vtd,Vts,Vtb}} ).

end_if.
end_if.

scalar h/h:(Higgs, mass Mh = 125, width wh =auto).
									
let l1={ne,e}, L1={Ne,E}.
let l2={nm,m}, L2={Nm,M}.
let l3={nl,l}, L3={Nl,L}.

let q1={u,d}, Q1={U,D}, q1a={u,Vud*d+Vus*s+Vub*b}, Q1a={U,Vud*D+Vus*S+Vub*B}.
let q2={c,s}, Q2={C,S}, q2a={c,Vcd*d+Vcs*s+Vcb*b}, Q2a={C,Vcd*D+Vcs*S+Vcb*B}. 
let q3={t,b}, Q3={T,B}, q3a={t,Vtd*d+Vts*s+Vtb*b}, Q3a={T,Vtd*D+Vts*S+Vtb*B}.


%================== U(1) gauge field 
let B1= -SW*Z+CW*A.

%================== SU(2) gauge field  and Faddeev-Popov ghosts 
let  W = { ('W+'+'W-')/Sqrt2,     i*('W+'-'W-')/Sqrt2,     CW*Z+SW*A        }.
do_if gauge_fixing==Feynman.
let  gh= { ('W+.c'+'W-.c')/Sqrt2, i*('W+.c'-'W-.c')/Sqrt2, CW*'Z.c'+SW*'A.c'}.
let  Gh= { ('W+.C'+'W-.C')/Sqrt2, i*('W+.C'-'W-.C')/Sqrt2, CW*'Z.C'+SW*'A.C'}. 
end_if.

%================== SU(2) and U(1) gauge couplings and VEV
let g2=EE/SW, g1=EE/CW.
let vevh=2*MW/EE*SW.      


%================== Self-interaction of gauge bosons

lterm -F**2/4  where  F=deriv^mu*B1^nu-deriv^nu*B1^mu.
lterm -F**2/4  where  F=deriv^mu*G^nu^a-deriv^nu*G^mu^a+i*GG*f_SU3^a^b^c*G^mu^b*G^nu^c.
lterm -F**2/4  where  F=deriv^mu*W^nu^a-deriv^nu*W^mu^a+g2*eps^a^b^c*W^mu^b*W^nu^c.


%================== left fermion interaction with gauge fields

lterm  	i*anti(psi)*gamma*(1-g5)/2*(deriv-i*g2*tau*W/2-i*Y*g1*B1)*psi
		where 
			psi=l1,  Y=-1/2;
			psi=l2,  Y=-1/2;
			psi=l3,  Y=-1/2;
			psi=q1a, Y= 1/6;
			psi=q2a, Y= 1/6;
			psi=q3a, Y= 1/6.

%================== right fermion interaction with gauge fields

lterm  	i*anti(psi)*gamma*(1+g5)/2*(deriv - i*Y*g1*B1)*psi
		where 
			psi=e,Y= -1;
			psi=m,Y= -1;
			psi=l,Y= -1;
			psi=u, Y=  2/3;
			psi=c, Y=  2/3;
			psi=t, Y=  2/3;
			psi=d, Y= -1/3;
			psi=s, Y= -1/3;
			psi=b, Y= -1/3.

%================== quark-gluon interaction

lterm  GG*anti(psi)*lambda*gamma*G*psi where
	psi=q1; psi=q2; psi=q3.


%=================  Higgs doublet
do_if gauge_fixing==Feynman.
    let shd = {  i*'W+.f',  (vev(vevh)+h-i*'Z.f')/Sqrt2 }.   % scalar Higgs doublet
do_else.
    let shd = { 0,          (vev(vevh)+h)/Sqrt2         }. 
end_if.

let  shD = anti(shd).                                        % scalar Higgs anti-doublet. 

%===  Generalized dirivative of Higgs doublet 

let Dshd^mu^a = (deriv^mu-i*g1/2*B1^mu)*shd^a -i*g2/2*tau^a^b^c*W^mu^c*shd^b.
let DshD^mu^a = (deriv^mu+i*g1/2*B1^mu)*shD^a +i*g2/2*tau^b^a^c*W^mu^c*shD^b.
    
%===  Lagrangian 
lterm DshD*Dshd - Mh**2/vevh**2/2*(shd*shD-vevh**2/2)**2.


%================ Yukawa terms  
%=== for bottom quarks 
lterm  -M*Sqrt2/vevh*anti(pl)*(1+g5)/2*pr*shd +  AddHermConj 
    where
	M=Vud*0,  pl=q1a, pr=d;          % 0 stands for Md 
	M=Vus*Ms, pl=q1a, pr=s;
	M=Vub*Mb, pl=q1a, pr=b;
	M=Vcd*0,  pl=q2a, pr=d;
	M=Vcs*Ms, pl=q2a, pr=s;
	M=Vcb*Mb, pl=q2a, pr=b;
	M=Vtd*0,  pl=q3a, pr=d;
	M=Vts*Ms, pl=q3a, pr=s;
	M=Vtb*Mb, pl=q3a, pr=b.

%===  for up quarks
   let eps2 = {{0,-1},{1,0}}.  % eps2^{12} =1 LanHEP reads matrices by columns 
   lterm   -M*Sqrt2/vevh*anti(pl)*(1+g5)/2*eps2*pr*shD +  AddHermConj
   where
%	M=0 ,  pl=q1a, pr=u;
	M=Mc,  pl=q2a, pr=c;
	M=Mtp, pl=q3a, pr=t.

%=== for charged leptons 
lterm  -M*Sqrt2/vevh*anti(pl)*(1+g5)/2*pr*shd +  AddHermConj 
    where
%	M=Mm,   pl=l2,  pr=m;
	M=Ml,  pl=l3,  pr=l.
	
%====================== gauge fixing terms and corresponding ghost terms

%=== for SU(3)
lterm  -1/2*(deriv*G)**2.
lterm -i*GG*f_SU3*ccghost(G)*G^mu*deriv^mu*ghost(G).


%=== for SU(2)*U(1) 

do_if gauge_fixing==Feynman.

lterm  -1/2*(deriv*A)**2  -1/2*(2*(deriv*'W+'+vevh*g2/2*'W+.f')*(deriv*'W-'+vevh*g2/2*'W-.f') +(deriv*Z+ vevh*g2/2/CW*'Z.f')**2).
                
lterm -g2*eps*deriv*Gh*gh*W.

lterm vevh*g2**2/4*(-h*('W-.C'*'W+.c' + 'W+.C'*'W-.c' + 'Z.C'*'Z.c'/CW**2)  +i*'Z.f'*('W-.C'*'W+.c' - 'W+.C'*'W-.c') ).

lterm -i*g2**2*vevh/4/CW*(
	'W+.f'*('W-.C'*'Z.c'*(1-2*SW**2)+'W-.c'*'Z.C' + 2*CW*SW*'W-.C'*'A.c')
      -	'W-.f'*('W+.C'*'Z.c'*(1-2*SW**2)+'W+.c'*'Z.C' + 2*CW*SW*'W+.C'*'A.c')  
                    ).
end_if.

do_if hgg==On.

external_func(alphaQCD,1).
external_func(lGGhiggs,2).
external_func(lAAhiggs,2).

parameter LAAh=-cabs(lAAhiggs(Mh,str(h))),
          LGGh=-cabs(lGGhiggs(Mh,str(h))),
          aQCDh =alphaQCD(Mh)/acos(-1),
          RQCDh=sqrt(1+149/12*aQCDh+68.6482*aQCDh**2-212.447*aQCDh**3).

lterm  LAAh*(shd*shD-vevh**2/2)/vevh*F**2  where 
             F=deriv^mu*A^nu-deriv^nu*A^mu.

lterm  LGGh*RQCDh*(shd*shD-vevh**2/2) /vevh*F**2 where
             F=deriv^mu*G^nu^a-deriv^nu*G^mu^a.

do_if h4G==On.
parameter  Maux=1.
vector  G1t/G1T: (G1T, mass Maux, color c8, aux ('!*')).
vector  G2t/G2T: (G2T, mass Maux, color c8, aux ('!*')).
scalar  x1/X1:   (x1,  mass Maux, aux ("!*")).  


lterm 1/vevh*LGGh*RQCDh*(shd*shD-vevh**2/2)*(F12*F11+F11*F12)
        where
        F11=deriv^mu*G^nu^a-deriv^nu*G^mu^a, F12=i*GG*f_SU3^a^b^c*G^mu^b*G^nu^c.

lterm   1/vevh*LGGh*RQCDh*(shd*shD-vevh**2/2)*(GG/2*f_SU3^a^b^c*G^n^a*G^m^b*'G1t.t'^m^n^c).
lterm                                           GG/2*f_SU3^a^b^c*G^n^a*G^m^b*'G1T.t'^m^n^c. 
lterm   1/vevh*LGGh*RQCDh*(shd*shD-vevh**2/2-vevh*h)*x1*Maux.
lterm   (GG/2*f_SU3^a^b^c*G^n^a*G^m^b*'G2t.t'^m^n^c)*(X1*Maux).
lterm    GG/2*f_SU3^a^b^c*G^n^a*G^m^b*'G2T.t'^m^n^c.


end_if. 

end_if.

SetAngle(1-SW**2=CW**2).
SetEM(A,EE).

