      SUBROUTINE NS_CHARGINO

************************************************************************
*
*     This subroutine computes the chargino decays
*
*     It is a generalisation of the corresponding routine from
*     SDECAY: A Fortran code for the decays of the supersymmetric
*             particles in the MSSM
*     by M. Muhlleitner (Karlsruhe, Inst. Technol.),
*	 A. Djouadi (Orsay, LPT & CERN, Theory Division),
*	 Y. Mambrini (Orsay, LPT),
*     Comput.Phys.Commun.168:46-70 (2005), hep-ph/0311167.
*     SDECAY should be cited whenever NMSDECAY is used.
*
************************************************************************

      IMPLICIT NONE

      INTEGER I,J,jsign,k,GRFLAG

      DOUBLE PRECISION amuv,lamv,amuvdiv
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION ql(5,2),qr(5,2),ol(5,2),or(5,2)
      DOUBLE PRECISION opl(2,2),opr(2,2),onl(5,5),onr(5,5)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION alup(2,2),aldo(2,2)
      DOUBLE PRECISION alstor(2,2),akstor(2,2),alsbot(2,2),aksbot(2,2)
      DOUBLE PRECISION hchichi(3,5,5),achichi(2,5,5),hchachaR(3,2,2),
     .         hchachaL(3,2,2),achachaR(2,2,2),achachaL(2,2,2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION asupr,asupl,asdownr,asdownl,aser,asel,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION qscal,multilim,flagmulti,flagqcd,flagloop
      DOUBLE PRECISION charsel(2),charser(2),charsmu1(2),charsmu2(2),
     .         charstau1(2),charstau2(2),charsne1(2),charsne2(2),
     .         charsnm1(2),charsnm2(2),charsnt1(2),charsnt2(2),
     .         charsupl(2),charsupr(2),charsdownl(2),charsdownr(2),
     .         charst1(2),charst2(2),charsb1(2),charsb2(2),
     .         qcdcharst1(2),qcdcharst2(2),qcdcharsb1(2),qcdcharsb2(2),
     .         qcdcharsupl(2),qcdcharsupr(2),qcdcharsdownl(2),
     .         qcdcharsdownr(2),charwneut(2,5),charhcneut(2,5),
     .         char2zchic1,char2Hchic1(3),char2Achic1(2)
      DOUBLE PRECISION xchar1el,xchar1mu
      DOUBLE PRECISION xchitau(2,5),xchiel(2,5),xchiup(2,5),xchimu(2,5),
     .         xchich(2,5),xchitop(2,5),xgluiupdb(2),xgluichsb(2),
     .         xgluitopbb(2),xchar1tau,xchar1nue,xchar1numu,
     .         xchar1nutau,xchar1up,xchar1dow,xchar1ch,xchar1str,
     .         xchar1top,xchar1bot
      DOUBLE PRECISION chartot2(2),chartot(2),chartotmulti(2)
      DOUBLE PRECISION chartot2lo(2),chartot2nlo(2)
      DOUBLE PRECISION brcharsel(2),brcharser(2),brcharsmu1(2),
     .         brcharsmu2(2),brcharstau1(2),brcharstau2(2),
     .         brcharsne1(2),brcharsne2(2),brcharsnm1(2),brcharsnm2(2),
     .         brcharsnt1(2),brcharsnt2(2),brcharsupl(2),brcharsupr(2),
     .         brcharsdownl(2),brcharsdownr(2),brcharst1(2),
     .         brcharst2(2),brcharsb1(2),brcharsb2(2),brcharwneut(2,5),
     .         brcharhcneut(2,5),brcharzchic,brcharHchic(3),
     .         brcharAchic(2),brntaunut(2,5),brnelnue(2,5),
     .         brnmunumu(2,5),brnupdb(2,5),brnchsb(2,5),brntopbb(2,5),
     .         brglupdb(2),brglchsb(2),brgltopbb(2),brchee,brchmumu,
     .         brchtautau,brchnene,brchnmunmu,brchntauntau,brchupup,
     .         brchdodo,brchchch,brchstst,brchtoptop,brchbotbot
      DOUBLE PRECISION NS_lamb
      DOUBLE PRECISION NS_realicorr,NS_gltchar,NS_grtchar
      DOUBLE PRECISION NS_glbchar,NS_grbchar
      DOUBLE PRECISION NS_gamtop1,NS_gamtop2,NS_gamglui1,NS_gamglui2,
     .         NS_gamglui3,NS_gam11,NS_gam12,NS_gamvirt,NS_gamreal,
     .         NS_gamcfdec
      DOUBLE PRECISION NS_ftotqcd
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION charWgra(2),charHCgra(2)
      DOUBLE PRECISION brcharWgra(2),brcharHCgra(2),brneutGAMgra(5),
     .         brneutZgra(5),brneutHgra(5,3),brneutAgra(5,2),
     .         brgluiGLUgra,brselEgra,brserEgra,brsmu1MUgra,
     .         brsmu2MUgra,brstau1TAUgra,brstau2TAUgra,brsneNEgra,
     .         brsnmNMgra,brsntNTgra,brsulUgra,brsurUgra,brsdlDgra,
     .         brsdrDgra,brst1Tgra,brst2Tgra,brsb1Bgra,brsb2Bgra
      DOUBLE PRECISION KNG(5),KNZ(5),KNH(5,3),KNA(5,2),KCW(2),KCH(2)
      DOUBLE PRECISION M32,CGR,MPL

      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_refscale/qscal
      COMMON/CHARGINO_2GAMMA/charsel,charser,charsmu1,charsmu2,
     .         charstau1,charstau2,charsne1,charsne2,charsnm1,
     .         charsnm2,charsnt1,charsnt2,charsupl,charsupr,
     .         charsdownl,charsdownr,charst1,charst2,charsb1,
     .         charsb2,qcdcharst1,qcdcharst2,qcdcharsb1,qcdcharsb2,
     .         qcdcharsupl,qcdcharsupr,qcdcharsdownl,qcdcharsdownr,
     .         charwneut,charhcneut,char2zchic1,char2Hchic1,char2Achic1
      COMMON/CHARGINO_3GAMMA/xchar1el,xchar1mu,xchitau,xchiel,xchiup,
     .         xchimu,xchich,xchitop,xgluiupdb,xgluichsb,
     .         xgluitopbb,xchar1tau,xchar1nue,xchar1numu,
     .         xchar1nutau,xchar1up,xchar1dow,xchar1ch,xchar1str,
     .         xchar1top,xchar1bot
      COMMON/CHARGINO_WIDTH/chartot2,chartot,chartotmulti
      COMMON/CHARGINO_BR_2BD/brcharsel,brcharser,brcharsmu1,
     .         brcharsmu2,brcharstau1,brcharstau2,
     .         brcharsne1,brcharsne2,brcharsnm1,brcharsnm2,
     .         brcharsnt1,brcharsnt2,brcharsupl,brcharsupr,
     .         brcharsdownl,brcharsdownr,brcharst1,
     .         brcharst2,brcharsb1,brcharsb2,brcharwneut,
     .         brcharhcneut,brcharzchic,brcharHchic,
     .         brcharAchic
      COMMON/CHARGINO_BR_3BD/brntaunut,brnelnue,brnmunumu,
     .         brnupdb,brnchsb,brntopbb,
     .         brglupdb,brglchsb,brgltopbb,
     .         brchee,brchmumu,brchtautau,brchnene,
     .         brchnmunmu,brchntauntau,brchupup,brchdodo,
     .         brchchch,brchstst,brchtoptop,brchbotbot
      COMMON/GRAVITINO/brcharWgra,brcharHCgra,brneutGAMgra,
     .         brneutZgra,brneutHgra,brneutAgra,
     .         brgluiGLUgra,brselEgra,brserEgra,brsmu1MUgra,
     .         brsmu2MUgra,brstau1TAUgra,brstau2TAUgra,brsneNEgra,
     .         brsnmNMgra,brsntNTgra,brsulUgra,brsurUgra,brsdlDgra,
     .         brsdrDgra,brst1Tgra,brst2Tgra,brsb1Bgra,brsb2Bgra
      COMMON/GRAVCOUP/KNG,KNZ,KNH,KNA,KCW,KCH
      COMMON/M32/M32,CGR,MPL,GRFLAG
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asupr,asupl,asdownr,asdownl,aser,asel,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/NS_coup4/opl,opr,onl,onr
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup7/alup,aldo
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_coupNMSSM/hchichi,achichi,hchachaR,hchachaL,achachaR,
     .         achachaL
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor
      COMMON/NS_multilim/multilim
      COMMON/NS_FLAGS/flagmulti,flagqcd,flagloop

      EXTERNAL NS_lamb
      EXTERNAL NS_realicorr,NS_gltchar,NS_grtchar
      EXTERNAL NS_glbchar,NS_grbchar
      EXTERNAL NS_gamtop1,NS_gamtop2,NS_gamglui1,NS_gamglui2,
     .         NS_gamglui3,NS_gam11,NS_gam12,NS_gamvirt,NS_gamreal,
     .         NS_gamcfdec
      EXTERNAL NS_ftotqcd

*  Initialization

      do i=1,2
         chartot(i) = 0d0
         chartot2(i) = 0d0
         chartot2lo(i) = 0d0
         chartot2nlo(i) = 0d0
         chartotmulti(i) = 0d0
         charsel(i) = 0d0
         charser(i) = 0d0
         charsmu1(i) = 0d0
         charsmu2(i) = 0d0
         charstau1(i) = 0d0
         charstau2(i) = 0d0
         charsne1(i) = 0d0
         charsne2(i) = 0d0
         charsnm1(i) = 0d0
         charsnm2(i) = 0d0
         charsnt1(i) = 0d0
         charsnt2(i) = 0d0
         charsupl(i) = 0d0
         charsupr(i) = 0d0
         charsdownl(i) = 0d0
         charsdownr(i) = 0d0
         charst1(i) = 0d0
         charst2(i) = 0d0
         charsb1(i) = 0d0
         charsb2(i) = 0d0
         charWgra(i) = 0d0
         charHCgra(i) = 0d0
         qcdcharsupl(i)   = 0d0
         qcdcharsupr(i)   = 0d0
         qcdcharsdownl(i) = 0d0
         qcdcharsdownr(i) = 0d0
         qcdcharst1(i)    = 0d0
         qcdcharst2(i)    = 0d0
         qcdcharsb1(i)    = 0d0
         qcdcharsb2(i)    = 0d0
         xgluiupdb(i)     = 0d0
         xgluichsb(i)     = 0d0
         xgluitopbb(i)    = 0d0
         do j=1,5
            charwneut(i,j) = 0d0
            charhcneut(i,j) = 0d0
            xchitau(i,j) = 0d0
            xchiel(i,j) = 0d0
            xchiup(i,j) = 0d0
            xchimu(i,j) = 0d0
            xchich(i,j) = 0d0
            xchitop(i,j) = 0d0
         enddo
      enddo
      char2zchic1 = 0d0
      do i=1,3
         char2Hchic1(i) = 0d0
      enddo
      do i=1,2
         char2Achic1(i) = 0d0
      enddo
      xchar1tau     = 0d0
      xchar1nue     = 0d0
      xchar1numu    = 0d0
      xchar1nutau   = 0d0
      xchar1up      = 0d0
      xchar1dow     = 0d0
      xchar1ch      = 0d0
      xchar1str     = 0d0
      xchar1top     = 0d0
      xchar1bot     = 0d0
**********************************************************************
c             BEGIN 2-BODY DECAYS
***********************************************************************
c chi+_1/chi+_2 --> selectron_L+ + neutrino_e

      do i=1,2
         if(asel.le.amchar(i)) then
            charsel(i)=g2s/2d0*(ale(1,i)**2)*
     .        (amchar(i)**2-asel**2)*NS_lamb(0d0,asel/amchar(i))
     .        /(16*pi*amchar(i))
         else
            charsel(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> selectron_R+ + neutrino_e

      do i=1,2
         if(aser.le.amchar(i)) then
            charser(i)=g2s/2d0*(ale(2,i)**2)*
     .        (amchar(i)**2-aser**2)*NS_lamb(0d0,aser/amchar(i))
     .        /(16*pi*amchar(i))
         else
            charser(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> smu_1+ + neutrino_mu

      do i=1,2
         if (asmu1.le.amchar(i)) then
            charsmu1(i)=g2s/2d0*(almu(1,i)**2)*
     .        (amchar(i)**2-asmu1**2)*NS_lamb(0d0,asmu1/amchar(i))
     .        /(16*pi*amchar(i))
         else
            charsmu1(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> smu_2+ + neutrino_mu

      do i=1,2
         if (asmu2.le.amchar(i)) then
            charsmu2(i)=g2s/2d0*(almu(2,i)**2)*
     .        (amchar(i)**2-asmu2**2)*NS_lamb(0d0,asmu2/amchar(i))
     .        /(16*pi*amchar(i))
         else
            charsmu2(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> stau_1+ + neutrino_tau

      do i=1,2
         if (astau1.le.amchar(i)) then
            charstau1(i)=g2s/2d0*(altau(1,i)**2)*
     .        (amchar(i)**2-astau1**2)*NS_lamb(0d0,astau1/amchar(i))
     .        /(16*pi*amchar(i))
         else
            charstau1(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> stau_2+ + neutrino_tau

      do i=1,2
         if (astau2.le.amchar(i)) then
            charstau2(i)=g2s/2d0*((altau(2,i)**2)*
     .        (amchar(i)**2-astau2**2))*NS_lamb(0d0,astau2/amchar(i))
     .        /(16*pi*amchar(i))
         else
            charstau2(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> sneutrino_e1 + e+

      do i=1,2
         if(asne1.le.amchar(i)) then
            charsne1(i)=g2s/2d0*alsne(1,i)**2*
     .        (amchar(i)**2-asne1**2)*NS_lamb(0d0,asne1/amchar(i))
     .        /(16*pi*amchar(i))
         else
            charsne1(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> sneutrino_e2 + e+

      do i=1,2
         if(asne2.le.amchar(i)) then
            charsne2(i)=g2s/2d0*alsne(2,i)**2*
     .        (amchar(i)**2-asne2**2)*NS_lamb(0d0,asne2/amchar(i))
     .        /(16*pi*amchar(i))
         else
            charsne2(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> sneutrino_mu1 + mu+
      do i=1,2
         if ((asnmu1+mmuon).le.amchar(i)) then
            charsnm1(i)=g2s/2d0*(
     .          (alsnm(1,i)**2+blsnm(1,i)**2)*
     .          (amchar(i)**2-asnmu1**2+mmuon**2)
     .          +4d0*alsnm(1,i)*blsnm(1,i)*xmchar(i)*mmuon)*
     .          NS_lamb(mmuon/amchar(i),asnmu1/amchar(i))
     .          /(16*pi*amchar(i))
         else
            charsnm1(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> sneutrino_mu2 + mu+

      do i=1,2
         if ((asnmu2+mmuon).le.amchar(i)) then
            charsnm2(i)=g2s/2d0*(
     .          (alsnm(2,i)**2+blsnm(2,i)**2)*
     .          (amchar(i)**2-asnmu2**2+mmuon**2)
     .          +4d0*alsnm(2,i)*blsnm(2,i)*xmchar(i)*mmuon)*
     .          NS_lamb(mmuon/amchar(i),asnmu2/amchar(i))
     .          /(16*pi*amchar(i))
         else
            charsnm2(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> sneutrino_tau1 + tau+
      do i=1,2
         if ((asntau1+mtau).le.amchar(i)) then
            charsnt1(i)=g2s/2d0*(
     .          (alsnt(1,i)**2+blsnt(1,i)**2)*
     .          (amchar(i)**2-asntau1**2+mtau**2)
     .          +4d0*alsnt(1,i)*blsnt(1,i)*xmchar(i)*mtau)*
     .          NS_lamb(mtau/amchar(i),asntau1/amchar(i))
     .          /(16*pi*amchar(i))
         else
            charsnt1(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> sneutrino_tau2 + tau+

      do i=1,2
         if ((asntau2+mtau).le.amchar(i)) then
            charsnt2(i)=g2s/2d0*(
     .          (alsnt(2,i)**2+blsnt(2,i)**2)*
     .          (amchar(i)**2-asntau2**2+mtau**2)
     .          +4d0*alsnt(2,i)*blsnt(2,i)*xmchar(i)*mtau)*
     .          NS_lamb(mtau/amchar(i),asntau2/amchar(i))
     .          /(16*pi*amchar(i))
         else
            charsnt2(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> stop1 + bbar

      do i=1,2
         if ((ast1+mbp).le.amchar(i)) then
                      charst1(i)=3d0*g2s/2d0*(
     .           (alstor(1,i)**2+akstor(1,i)**2)*
     .           (amchar(i)**2-ast1**2+mbp**2)
     .           +4d0*alstor(1,i)*akstor(1,i)*xmchar(i)*mbp)*
     .          NS_lamb(mbp/amchar(i),ast1/amchar(i))/(16*pi*amchar(i))
         else
            charst1(i)=0d0
         endif
      enddo
C -- QCD corrections --
      if(flagqcd.eq.1d0) then
      do i=1,2
         if ((ast1+mbp).le.amchar(i)) then
            jsign = 0
            amuvdiv = amuv/qscal
            qcdcharst1(i)= 0d0
            qcdcharst1(i)= -g2s/24d0/pi**2/amchar(i)*g3s/(4d0*pi)*
     .           3d0/2d0*
     .           ((akstor(1,i)*NS_gltchar(1,i,amuv,amuvdiv,lamv)
     .           +alstor(1,i)*NS_grtchar(1,i,amuv,amuvdiv,lamv))*
     .           (-ast1**2+mbp**2+amchar(i)**2)
     .           +2d0*(akstor(1,i)*NS_grtchar(1,i,amuv,amuvdiv,lamv)
     .           +alstor(1,i)*NS_gltchar(1,i,amuv,amuvdiv,lamv))*
     .           mbp*xmchar(i))*NS_lamb(mbp/amchar(i),ast1/amchar(i))
     .           +g2s/(6d0*pi**2*amchar(i))*g3s/(4d0*pi)*
     .           3d0/2d0*(-1d0)*
     .           NS_realicorr(mbp,amchar(i),ast1,lamv,2,jsign,1,i,1)
         else
            qcdcharst1(i)= 0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> stop2 + bbar

      do i=1,2
         if ((ast2+mbp).le.amchar(i)) then
            charst2(i)=3d0*g2s/2d0*(
     .           (alstor(2,i)**2+akstor(2,i)**2)*
     .           (amchar(i)**2-ast2**2+mbp**2)
     .           +4d0*alstor(2,i)*akstor(2,i)*xmchar(i)*mbp)*
     .          NS_lamb(mbp/amchar(i),ast2/amchar(i))/(16*pi*amchar(i))
         else
            charst2(i)=0d0
         endif
      enddo

c -- QCD corrections --
      if(flagqcd.eq.1d0) then
      do i=1,2
         if ((ast2+mbp).le.amchar(i)) then
            jsign = 0
            amuvdiv = amuv/qscal
            qcdcharst2(i)=0d0
            qcdcharst2(i)= -g2s/24d0/pi**2/amchar(i)*g3s/(4d0*pi)*
     .           3d0/2d0*
     .           ((akstor(2,i)*NS_gltchar(2,i,amuv,amuvdiv,lamv)
     .            +alstor(2,i)*NS_grtchar(2,i,amuv,amuvdiv,lamv))*
     .           (-ast2**2+mbp**2+amchar(i)**2)
     .           +2d0*(akstor(2,i)*NS_grtchar(2,i,amuv,amuvdiv,lamv)
     .                 +alstor(2,i)*NS_gltchar(2,i,amuv,amuvdiv,lamv))*
     .           mbp*xmchar(i))*NS_lamb(mbp/amchar(i),ast2/amchar(i))
     .           +g2s/(6d0*pi**2*amchar(i))*g3s/(4d0*pi)*
     .           3d0/2d0*(-1d0)*
     .           NS_realicorr(mbp,amchar(i),ast2,lamv,2,jsign,2,i,1)
         else
            qcdcharst2(i)= 0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> sbottom1* + top

      do i=1,2
         if ((asb1+mt).le.amchar(i)) then
            charsb1(i)=3d0*g2s/2d0*(
     .           (alsbot(1,i)**2+aksbot(1,i)**2)*
     .           (amchar(i)**2-asb1**2+mt**2)
     .           +4d0*alsbot(1,i)*aksbot(1,i)*xmchar(i)*mt)*
     .           NS_lamb(mt/amchar(i),asb1/amchar(i))/(16*pi*amchar(i))
         else
            charsb1(i)=0d0
         endif
      enddo

c -- QCD corrections --
      if(flagqcd.eq.1d0) then
      do i=1,2
         if ((asb1+mt).le.amchar(i)) then
            jsign = 0
            amuvdiv = amuv/qscal
            qcdcharsb1(i)= 0d0
            qcdcharsb1(i)= -g2s/24d0/pi**2/amchar(i)*g3s/(4d0*pi)*
     .           3d0/2d0*
     .           ((aksbot(1,i)*NS_glbchar(1,i,amuv,amuvdiv,lamv)
     .            +alsbot(1,i)*NS_grbchar(1,i,amuv,amuvdiv,lamv))*
     .           (-asb1**2+mt**2+amchar(i)**2)
     .           +2d0*(aksbot(1,i)*NS_grbchar(1,i,amuv,amuvdiv,lamv)
     .                 +alsbot(1,i)*NS_glbchar(1,i,amuv,amuvdiv,lamv))*
     .           mt*xmchar(i))*NS_lamb(mt/amchar(i),asb1/amchar(i))
     .           +g2s/(6d0*pi**2*amchar(i))*g3s/(4d0*pi)*
     .           3d0/2d0*(-1d0)*
     .           NS_realicorr(mt,amchar(i),asb1,lamv,2,jsign,1,i,2)
         else
            qcdcharsb1(i)= 0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> sbottom2* + top

      do i=1,2
        if ((asb2+mt).le.amchar(i)) then
            charsb2(i)=3d0*g2s/2d0*(
     .           (alsbot(2,i)**2+aksbot(2,i)**2)*
     .           (amchar(i)**2-asb2**2+mt**2)
     .           +4d0*alsbot(2,i)*aksbot(2,i)*xmchar(i)*mt)*
     .           NS_lamb(mt/amchar(i),asb2/amchar(i))/(16*pi*amchar(i))
         else
            charsb2(i)=0d0
         endif
      enddo

c -- QCD corrections --
      if(flagqcd.eq.1d0) then
      do i=1,2
         if ((asb2+mt).le.amchar(i)) then
            jsign = 0
            amuvdiv = amuv/qscal
            qcdcharsb2(i)=0d0
            qcdcharsb2(i)= -g2s/24d0/pi**2/amchar(i)*g3s/(4d0*pi)*
     .           3d0/2d0*
     .           ((aksbot(2,i)*NS_glbchar(2,i,amuv,amuvdiv,lamv)
     .            +alsbot(2,i)*NS_grbchar(2,i,amuv,amuvdiv,lamv))*
     .           (-asb2**2+mt**2+amchar(i)**2)
     .           +2d0*(aksbot(2,i)*NS_grbchar(2,i,amuv,amuvdiv,lamv)
     .                 +alsbot(2,i)*NS_glbchar(2,i,amuv,amuvdiv,lamv))*
     .           mt*xmchar(i))*NS_lamb(mt/amchar(i),asb2/amchar(i))
     .           +g2s/(6d0*pi**2*amchar(i))*g3s/(4d0*pi)*
     .           3d0/2d0*(-1d0)*
     .           NS_realicorr(mt,amchar(i),asb2,lamv,2,jsign,2,i,2)
         else
            qcdcharsb2(i)= 0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> supl + downbar

      do i=1,2
         if (asupl.le.amchar(i)) then
            charsupl(i)=3d0*g2s/2d0*alup(1,i)**2*
     .           (amchar(i)**2-asupl**2)*
     .           NS_lamb(0d0,asupl/amchar(i))/(16*pi*amchar(i))
         else
            charsupl(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,2
         if(amchar(i).ge.asupl) then
            qcdcharsupl(i)=0d0
            qcdcharsupl(i)=4d0/3d0*g3s/(4d0*pi)/pi*charsupl(i)*
     .       NS_ftotqcd(asupl**2/amchar(i)**2,mgluino**2/amchar(i)**2)
         else
            qcdcharsupl(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> supr + downbar

      do i=1,2
         if (asupr.le.amchar(i)) then
            charsupr(i)=3d0*g2s/2d0*alup(2,i)**2*
     .           (amchar(i)**2-asupr**2)*
     .           NS_lamb(0d0,asupr/amchar(i))/(16*pi*amchar(i))
         else
            charsupr(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,2
         if(amchar(i).ge.asupr) then
            qcdcharsupr(i)=0d0
            qcdcharsupr(i)=4d0/3d0*g3s/(4d0*pi)/pi*charsupr(i)*
     .       NS_ftotqcd(asupr**2/amchar(i)**2,mgluino**2/amchar(i)**2)
         else
            qcdcharsupr(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> sdownl* + up

      do i=1,2
         if(asdownl.le.amchar(i)) then
            charsdownl(i)=3d0*g2s/2d0*aldo(1,i)**2*
     .           (amchar(i)**2-asdownl**2)*
     .           NS_lamb(0d0,asdownl/amchar(i))/(16*pi*amchar(i))
         else
            charsdownl(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,2
         if(amchar(i).ge.asdownl) then
            qcdcharsdownl(i)=0d0
            qcdcharsdownl(i)=4d0/3d0*g3s/(4d0*pi)/pi*charsdownl(i)*
     .      NS_ftotqcd(asdownl**2/amchar(i)**2,mgluino**2/amchar(i)**2)
         else
            qcdcharsdownl(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> sdownr* + up

      do i=1,2
         if(asdownr.le.amchar(i)) then
            charsdownr(i)=3d0*g2s/2d0*aldo(2,i)**2*
     .           (amchar(i)**2-asdownr**2)*
     .           NS_lamb(0d0,asdownr/amchar(i))/(16*pi*amchar(i))
         else
            charsdownr(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,2
         if(amchar(i).ge.asdownr) then
            qcdcharsdownr(i)=0d0
            qcdcharsdownr(i)=4d0/3d0*g3s/(4d0*pi)/pi*charsdownr(i)*
     .      NS_ftotqcd(asdownr**2/amchar(i)**2,mgluino**2/amchar(i)**2)
         else
            qcdcharsdownr(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c chi+_2 --> chi+_1 + Z

      if ((mz+amchar(1)).le.amchar(2)) then
         char2zchic1=g2s/2d0/cw**2*(-12d0*xmchar(1)*xmchar(2)*
     .        opl(1,2)*opr(1,2)+(opl(1,2)**2+opr(1,2)**2)*
     .        ((amchar(1)**2+amchar(2)**2-mz**2)+
     .        (amchar(2)**2+mz**2-amchar(1)**2)*
     .        (amchar(2)**2-amchar(1)**2-mz**2)/mz**2))
     .        *NS_lamb(amchar(1)/amchar(2),mz/amchar(2))
     .        /(16d0*pi*amchar(2))
      else
         char2zchic1=0d0
      endif
c -------------------------------------------------------------------- c
c chi+_2 --> chi+_1 + H(K) char2Hchic1(3) [NMSSM!!!]
      do K=1,3
         if ((SMASS(K)+amchar(1)).le.amchar(2)) then
         char2Hchic1(K)=g2s/2d0*(
     .        (hchachaL(K,1,2)**2+hchachaR(K,1,2)**2)*
     .        (amchar(1)**2+amchar(2)**2-SMASS(K)**2)+
     .        4d0*hchachaL(K,1,2)*hchachaR(K,1,2)*
     .        xmchar(1)*xmchar(2))*
     .        NS_lamb(amchar(1)/amchar(2),SMASS(K)/amchar(2))
     .        /(16d0*pi*amchar(2))
      else
         char2Hchic1(K)=0d0
      endif
      enddo
c -------------------------------------------------------------------- c
c chi+_2 --> chi+_1 + A(K)     char2Achic1(2) [NMSSM!!!]

      do K=1,2
      if ((PMASS(K)+amchar(1)).le.amchar(2)) then
         char2Achic1(K)=g2s/2d0*(
     .        (achachaL(K,1,2)**2+achachaR(K,1,2)**2)*
     .        (amchar(1)**2+amchar(2)**2-PMASS(K)**2)+
     .        4d0*achachaL(K,1,2)*achachaR(K,1,2)*
     .        xmchar(1)*xmchar(2))*
     .        NS_lamb(amchar(1)/amchar(2),PMASS(K)/amchar(2))
     .        /(16d0*pi*amchar(2))
      else
         char2Achic1(K)=0d0
      endif
      enddo
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + H+ [NMSSM!!!]

      do i=1,2
         do j=1,5
            if ((cmass+amneut(j)).le.amchar(i)) then
               charhcneut(i,j)=g2s/2d0*((ql(j,i)**2+qr(j,i)**2)*
     .              (amneut(j)**2+amchar(i)**2-cmass**2)+4d0*ql(j,i)*
     .              qr(j,i)*xmneut(j)*xmchar(i))*
     .              NS_lamb(amneut(j)/amchar(i),cmass/amchar(i))
     .              /(16d0*pi*amchar(i))
            else
               charhcneut(i,j)=0d0
            endif
         enddo
      enddo

c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + W+ [NMSSM!!!]

      do i=1,2
         do j=1,5
            if ((mw+amneut(j)).le.amchar(i)) then
               charwneut(i,j)=g2s/2d0*(
     .              -12d0*xmneut(j)*xmchar(i)*ol(j,i)*
     .              or(j,i)+(ol(j,i)**2+or(j,i)**2)*((amneut(j)**2+
     .              amchar(i)**2-mw**2)+
     .              (amchar(i)**2+mw**2-amneut(j)**2)*
     .              (amchar(i)**2-amneut(j)**2-mw**2)/mw**2))*
     .              NS_lamb(amneut(j)/amchar(i),mw/amchar(i))
     .              /(16d0*pi*amchar(i))
            else
               charwneut(i,j)=0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> W+ + gravitino

      if(GRFLAG.EQ.1)then
        do i=1,2
          if ((mw+M32).le.amchar(i)) then
            charWgra(i) = KCW(i)*amchar(i)**5/(48d0*PI*MPL**2*M32**2)
     .                    *(1d0-mw**2/amchar(i)**2)**4
          else
            charWgra(i) = 0d0
          endif
        enddo
      endif
c -------------------------------------------------------------------- c
c chi+_1/chi+_2 --> HC + gravitino

      if(GRFLAG.EQ.1)then
        do i=1,2
          if ((CMASS+M32).le.amchar(i)) then
            charHCgra(i) = KCH(i)*amchar(i)**5/(48d0*PI*MPL**2*M32**2)
     .                     *(1d0-CMASS**2/amchar(i)**2)**4
          else
            charHCgra(i) = 0d0
          endif
        enddo
      endif
c -------------------------------------------------------------------- c
c  ------------------ TWO BODY DECAY WIDTHS -------------------------- c
c -------------------------------------------------------------------- c

      do i=1,2
         chartot2lo(i) = charsel(i)+charser(i)+charsmu1(i)+charsmu2(i)
     .     +charstau1(i)+charstau2(i)+charsne1(i)+charsne2(i)
     .     +charsnm1(i)+charsnm2(i)+charsnt1(i)+charsnt2(i)
     .     +2d0*(charsupl(i)+charsupr(i)+2d0*charsdownl(i)
     .     +charsdownr(i))+charst1(i)+charst2(i)+charsb1(i)
     .     +charsb2(i)+charwneut(i,1)+charwneut(i,2)+charwneut(i,3)
     .     +charwneut(i,4)+charwneut(i,5)+charhcneut(i,1)
     .     +charhcneut(i,2)+charhcneut(i,3)+charhcneut(i,4)
     .     +charhcneut(i,5)+charWgra(i)+charHCgra(i)
      enddo
      chartot2lo(2) = chartot2lo(2)+char2zchic1+
     .     char2Hchic1(1)+char2Hchic1(2)+char2Hchic1(3)+
     .     char2Achic1(1)+char2Achic1(2)
c
c     ERROR RELATED TO FLAGQCD IS CORRECTED
         do i=1,2
            chartot2nlo(i) = chartot2lo(i)+qcdcharst1(i)+qcdcharst2(i)+
     .           qcdcharsb1(i)+qcdcharsb2(i)+2d0*(qcdcharsupl(i)+
     .           qcdcharsupr(i)+qcdcharsdownl(i)+qcdcharsdownr(i))
         enddo
         do i=1,2
            chartot2(i) = chartot2nlo(i)
         enddo

c -------------------------------------------------------------------- c
C ------------------------- THREE BODY DECAYS ------------------------ c
c -------------------------------------------------------------------- c
c
      if(flagmulti.eq.1d0) then

      CALL NS_xintegchipm

      do i=1,2
            chartotmulti(i)=0d0
            do j=1,5
               chartotmulti(i) = chartotmulti(i)+xchitau(i,j)+
     .              xchiel(i,j)+xchiup(i,j)+xchimu(i,j)+xchich(i,j)+
     .              xchitop(i,j)

            enddo
            chartotmulti(i) = chartotmulti(i)+xgluiupdb(i)+
     .           xgluichsb(i)+xgluitopbb(i)

            if(i.eq.2) then
               chartotmulti(i) = chartotmulti(i)+xchar1tau+xchar1nue
     .              +xchar1numu+xchar1nutau+xchar1up+xchar1dow
     .              +xchar1ch+xchar1str+xchar1top+xchar1bot
     .              +xchar1el+xchar1mu

            endif
      enddo
      endif
c ---------- Consider 3-body decays only if > multilim ---------------c
      do i = 1,2
         if (chartotmulti(i).lt.multilim*chartot2(i))then
            chartotmulti(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c ------------------------ Total widths -------------------------------c
c -------------------------------------------------------------------- c

      do i=1,2
              chartot(i)=chartot2(i)+chartotmulti(i)
      enddo
c -------------------------------------------------------------------- c
c -------------------- Chargino branching ratios ------------------- c
c -------------------------------------------------------------------- c
c                       -- 2-body decays --
      if(flagqcd.eq.1d0) then
         do i=1,2
            charsupl(i)   = charsupl(i)+qcdcharsupl(i)
            charsupr(i)   = charsupr(i)+qcdcharsupr(i)
            charsdownl(i) = charsdownl(i)+qcdcharsdownl(i)
            charsdownr(i) = charsdownr(i)+qcdcharsdownr(i)
            charst1(i)    = charst1(i)+qcdcharst1(i)
            charst2(i)    = charst2(i)+qcdcharst2(i)
            charsb1(i)    = charsb1(i)+qcdcharsb1(i)
            charsb2(i)    = charsb2(i)+qcdcharsb2(i)
         enddo
      endif

      do i=1,2
       if(chartot(i).ne.0d0)then

         brcharsel(i) = charsel(i)/chartot(i)
         brcharser(i) = charser(i)/chartot(i)
         brcharsmu1(i) = charsmu1(i)/chartot(i)
         brcharsmu2(i) = charsmu2(i)/chartot(i)
         brcharstau1(i) = charstau1(i)/chartot(i)
         brcharstau2(i) = charstau2(i)/chartot(i)
         brcharsne1(i) = charsne1(i)/chartot(i)
         brcharsne2(i) = charsne2(i)/chartot(i)
         brcharsnm1(i) = charsnm1(i)/chartot(i)
         brcharsnm2(i) = charsnm2(i)/chartot(i)
         brcharsnt1(i) = charsnt1(i)/chartot(i)
         brcharsnt2(i) = charsnt2(i)/chartot(i)
         brcharsupl(i) = charsupl(i)/chartot(i)
         brcharsupr(i) = charsupr(i)/chartot(i)
         brcharsdownl(i) = charsdownl(i)/chartot(i)
         brcharsdownr(i) = charsdownr(i)/chartot(i)
         brcharst1(i) = charst1(i)/chartot(i)
         brcharst2(i) = charst2(i)/chartot(i)
         brcharsb1(i) = charsb1(i)/chartot(i)
         brcharsb2(i) = charsb2(i)/chartot(i)
         brcharWgra(i) = charWgra(i)/chartot(i)
         brcharHCgra(i) = charHCgra(i)/chartot(i)
         do j=1,5
            brcharwneut(i,j) = charwneut(i,j)/chartot(i)
            brcharhcneut(i,j) = charhcneut(i,j)/chartot(i)
         enddo

       else

         brcharsel(i) = 0d0
         brcharser(i) = 0d0
         brcharsmu1(i) = 0d0
         brcharsmu2(i) = 0d0
         brcharstau1(i) = 0d0
         brcharstau2(i) = 0d0
         brcharsne1(i) = 0d0
         brcharsne2(i) = 0d0
         brcharsnm1(i) = 0d0
         brcharsnm2(i) = 0d0
         brcharsnt1(i) = 0d0
         brcharsnt2(i) = 0d0
         brcharsupl(i) = 0d0
         brcharsupr(i) = 0d0
         brcharsdownl(i) = 0d0
         brcharsdownr(i) = 0d0
         brcharst1(i) = 0d0
         brcharst2(i) = 0d0
         brcharsb1(i) = 0d0
         brcharsb2(i) = 0d0
         brcharWgra(i) = 0d0
         brcharHCgra(i) = 0d0
         do j=1,5
            brcharwneut(i,j) = 0d0
            brcharhcneut(i,j) = 0d0
         enddo

       endif
      enddo

      if(chartot(2).ne.0d0)then

       brcharzchic  = char2zchic1/chartot(2)
       do K=1,3
          brcharHchic(K) = char2Hchic1(K)/chartot(2)
       enddo
       do K=1,2
          brcharAchic(K) = char2Achic1(K)/chartot(2)
       enddo

      else

       brcharzchic  = 0d0
       do K=1,3
          brcharHchic(K) = 0d0
       enddo
       do K=1,2
          brcharAchic(K) = 0d0
       enddo

      endif
c
c                     -- 3-body decays --
c
      do i=1,2
       if(chartotmulti(i).ne.0d0)then

         do j=1,5
            brntaunut(i,j) = xchitau(i,j)/chartot(i)
            brnelnue(i,j)  = xchiel(i,j)/chartot(i)
            brnmunumu(i,j) = xchimu(i,j)/chartot(i)
            brnupdb(i,j)   = xchiup(i,j)/chartot(i)
            brnchsb(i,j)   = xchich(i,j)/chartot(i)
            brntopbb(i,j)  = xchitop(i,j)/chartot(i)
         enddo
         brglupdb(i)  = xgluiupdb(i)/chartot(i)
         brglchsb(i)  = xgluichsb(i)/chartot(i)
         brgltopbb(i) = xgluitopbb(i)/chartot(i)

       else

         do j=1,5
            brntaunut(i,j) = 0d0
            brnelnue(i,j)  = 0d0
            brnmunumu(i,j) = 0d0
            brnupdb(i,j)   = 0d0
            brnchsb(i,j)   = 0d0
            brntopbb(i,j)  = 0d0
         enddo
         brglupdb(i)  = 0d0
         brglchsb(i)  = 0d0
         brgltopbb(i) = 0d0

       endif
      enddo

      if(chartotmulti(2).ne.0d0)then

       brchee	   = xchar1el/chartot(2)
       brchmumu     = xchar1mu/chartot(2)
       brchtautau   = xchar1tau/chartot(2)
       brchnene     = xchar1nue/chartot(2)
       brchnmunmu   = xchar1numu/chartot(2)
       brchntauntau = xchar1nutau/chartot(2)
       brchupup     = xchar1up/chartot(2)
       brchdodo     = xchar1dow/chartot(2)
       brchchch     = xchar1ch/chartot(2)
       brchstst     = xchar1str/chartot(2)
       brchtoptop   = xchar1top/chartot(2)
       brchbotbot   = xchar1bot/chartot(2)

      else

       brchee	   = 0d0
       brchmumu     = 0d0
       brchtautau   = 0d0
       brchnene     = 0d0
       brchnmunmu   = 0d0
       brchntauntau = 0d0
       brchupup     = 0d0
       brchdodo     = 0d0
       brchchch     = 0d0
       brchstst     = 0d0
       brchtoptop   = 0d0
       brchbotbot   = 0d0

      endif

      end

c ==================================================================== c
c                          chargino 3-body decays                      c
c ==================================================================== c

      SUBROUTINE NS_xintegchipm
*
      IMPLICIT NONE
      INTEGER nx1t,ny1t,ni,nj,i
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION NS_ay,NS_by,NS_ax,NS_bx
      DOUBLE PRECISION NS_chipmtau,NS_chipmmu,NS_chipmel,NS_chipmup,
     .       NS_chipmtop,NS_charel,NS_charmu,NS_chartau,NS_charnue,
     .       NS_charnum,NS_charnutau,NS_charup,NS_chardow,NS_charbot,
     .       NS_chartop,NS_gluiupdb,NS_gluitopbb
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION chartot2(2),chartot(2),chartotmulti(2)
      DOUBLE PRECISION xchar1el,xchar1mu
      DOUBLE PRECISION xchitau(2,5),xchiel(2,5),xchiup(2,5),xchimu(2,5),
     .          xchich(2,5),xchitop(2,5),xgluiupdb(2),xgluichsb(2),
     .          xgluitopbb(2),xchar1tau,xchar1nue,xchar1numu,
     .          xchar1nutau,xchar1up,xchar1dow,xchar1ch,xchar1str,
     .          xchar1top,xchar1bot
      DOUBLE PRECISION xmu1,xmu2,xmu3,sum1
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_nx1/nx1t,ny1t
      COMMON/NS_indices/ni,nj
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/CHARGINO_3GAMMA/xchar1el,xchar1mu,xchitau,xchiel,xchiup,
     .          xchimu,xchich,xchitop,xgluiupdb,xgluichsb,
     .          xgluitopbb,xchar1tau,xchar1nue,xchar1numu,
     .          xchar1nutau,xchar1up,xchar1dow,xchar1ch,xchar1str,
     .          xchar1top,xchar1bot
      COMMON/CHARGINO_WIDTH/chartot2,chartot,chartotmulti
*
      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_chipmtau,NS_chipmmu,NS_chipmel,NS_chipmup,
     .       NS_chipmtop,NS_charel,NS_charmu,NS_chartau,NS_charnue,
     .       NS_charnum,NS_charnutau,NS_charup,NS_chardow,NS_charbot,
     .       NS_chartop,NS_gluiupdb,NS_gluitopbb

c -------------------------------------------------------------------- c
c ------------------------- neutralino e+ nu_e ----------------------- c
c -------------------------------------------------------------------- c
       do ni = 1,2
         do nj = 1,5
            xmu1=0d0
            xmu2=0d0
            xmu3=amneut(nj)**2/amchar(ni)**2
            if(amchar(ni).gt.amneut(nj)) then
               CALL NS_integ2(NS_chipmel,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xchiel(ni,nj) = sum1/64d0/(2d0*pi)**3*amchar(ni)
            else
               xchiel(ni,nj) = 0d0
            endif
         enddo
       enddo
c -------------------------------------------------------------------- c
c ----------------------- neutralino mu+ nu_mu ----------------------- c
c -------------------------------------------------------------------- c
       do ni = 1,2
         do nj = 1,5
            xmu1=0d0
            xmu2=0d0
            xmu3=amneut(nj)**2/amchar(ni)**2
            if(amchar(ni).gt.amneut(nj)) then
               CALL NS_integ2(NS_chipmmu,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xchimu(ni,nj) = sum1/64d0/(2d0*pi)**3*amchar(ni)
            else
               xchimu(ni,nj) = 0d0
            endif
         enddo
       enddo
c -------------------------------------------------------------------- c
c ------------------------ neutralino tau nutau ---------------------- c
c -------------------------------------------------------------------- c
      do ni = 1,2
         do nj = 1,5
            xmu1=0d0
            xmu2=mtau**2/amchar(ni)**2
            xmu3=amneut(nj)**2/amchar(ni)**2
            if(amchar(ni).gt.(mtau+amneut(nj))) then
               CALL NS_integ2(NS_chipmtau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xchitau(ni,nj) = sum1/64d0/(2d0*pi)**3*amchar(ni)
            else
               xchitau(ni,nj) = 0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c ----------------------- neutralino up downbar ---------------------- c
c -------------------------------------------------------------------- c
      do ni = 1,2
         do nj = 1,5
            xmu1=0d0
            xmu2=0d0
            xmu3=amneut(nj)**2/amchar(ni)**2
            if(amchar(ni).gt.amneut(nj)) then
               CALL NS_integ2(NS_chipmup,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xchiup(ni,nj) = sum1/64d0/(2d0*pi)**3*amchar(ni)*3d0
            else
               xchiup(ni,nj) = 0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c -------------------- neutralino charm strangebar ------------------- c
c -------------------------------------------------------------------- c
      do i=1,5
         xchich(1,i) = xchiup(1,i)
         xchich(2,i) = xchiup(2,i)
      enddo
c -------------------------------------------------------------------- c
c --------------------- neutralino top bottombar --------------------- c
c -------------------------------------------------------------------- c
      do ni = 1,2
         do nj = 1,5
            xmu1=mt**2/amchar(ni)**2
            xmu2=mbp**2/amchar(ni)**2
            xmu3=amneut(nj)**2/amchar(ni)**2
            if(amchar(ni).gt.(amneut(nj)+mt+mbp)) then
               CALL NS_integ2(NS_chipmtop,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xchitop(ni,nj) = sum1/64d0/(2d0*pi)**3*amchar(ni)*3d0
            else
               xchitop(ni,nj) = 0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c --------------------------- chargino e+ e- ------------------------- c
c -------------------------------------------------------------------- c
      xmu1=0d0
      xmu2=0d0
      xmu3=amchar(1)**2/amchar(2)**2
      if(amchar(2).gt.amchar(1)) then
         CALL NS_integ2(NS_charel,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .               xmu3,nx1t,ny1t,sum1)
         xchar1el = sum1/64d0/(2d0*pi)**3*amchar(2)
      else
         xchar1el = 0d0
      endif
c -------------------------------------------------------------------- c
c ------------------------- chargino mu+ mu- ------------------------- c
c -------------------------------------------------------------------- c
      xmu1=0d0
      xmu2=0d0
      xmu3=amchar(1)**2/amchar(2)**2
      if(amchar(2).gt.amchar(1)) then
         CALL NS_integ2(NS_charmu,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .               xmu3,nx1t,ny1t,sum1)
         xchar1mu = sum1/64d0/(2d0*pi)**3*amchar(2)
      else
         xchar1mu = 0d0
      endif
c -------------------------------------------------------------------- c
c ------------------------- chargino tau+ tau- ----------------------- c
c -------------------------------------------------------------------- c
      xmu1=mtau**2/amchar(2)**2
      xmu2=mtau**2/amchar(2)**2
      xmu3=amchar(1)**2/amchar(2)**2
      if(amchar(2).gt.(amchar(1)+2d0*mtau)) then
         CALL NS_integ2(NS_chartau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .               xmu3,nx1t,ny1t,sum1)
         xchar1tau = sum1/64d0/(2d0*pi)**3*amchar(2)
      else
         xchar1tau = 0d0
      endif
c -------------------------------------------------------------------- c
c ---------------------- chargino nubar_e nu_e ----------------------- c
c -------------------------------------------------------------------- c
      xmu1=0d0
      xmu2=0d0
      xmu3=amchar(1)**2/amchar(2)**2
      if(amchar(2).gt.amchar(1)) then
         CALL NS_integ2(NS_charnue,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .               xmu3,nx1t,ny1t,sum1)
         xchar1nue = sum1/64d0/(2d0*pi)**3*amchar(2)
      else
         xchar1nue = 0d0
      endif
c -------------------------------------------------------------------- c
c --------------------- chargino nubar_mu nu_mu ---------------------- c
c -------------------------------------------------------------------- c
      xmu1=0d0
      xmu2=0d0
      xmu3=amchar(1)**2/amchar(2)**2
      if(amchar(2).gt.amchar(1)) then
         CALL NS_integ2(NS_charnum,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .               xmu3,nx1t,ny1t,sum1)
         xchar1numu = sum1/64d0/(2d0*pi)**3*amchar(2)
      else
         xchar1numu = 0d0
      endif
c -------------------------------------------------------------------- c
c -------------------- chargino nubar_tau nu_tau --------------------- c
c -------------------------------------------------------------------- c
      xmu1=0d0
      xmu2=0d0
      xmu3=amchar(1)**2/amchar(2)**2
      if(amchar(2).gt.amchar(1)) then
         CALL NS_integ2(NS_charnutau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .               xmu2,xmu3,nx1t,ny1t,sum1)
         xchar1nutau = sum1/64d0/(2d0*pi)**3*amchar(2)
      else
         xchar1nutau = 0d0
      endif
c -------------------------------------------------------------------- c
c ------------------------ chargino up upbar ------------------------- c
c -------------------------------------------------------------------- c
      xmu1=0d0
      xmu2=0d0
      xmu3=amchar(1)**2/amchar(2)**2
      if(amchar(2).gt.amchar(1)) then
         CALL NS_integ2(NS_charup,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .               xmu3,nx1t,ny1t,sum1)
         xchar1up = 3d0*sum1/64d0/(2d0*pi)**3*amchar(2)
      else
         xchar1up = 0d0
      endif
c -------------------------------------------------------------------- c
c -------------------- chargino charm charmbar ----------------------- c
c -------------------------------------------------------------------- c
      xchar1ch = xchar1up
c -------------------------------------------------------------------- c
c ---------------------- chargino down downbar ----------------------- c
c -------------------------------------------------------------------- c
      xmu1=0d0
      xmu2=0d0
      xmu3=amchar(1)**2/amchar(2)**2
      if(amchar(2).gt.amchar(1)) then
         CALL NS_integ2(NS_chardow,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .               xmu3,nx1t,ny1t,sum1)
         xchar1dow = 3d0*sum1/64d0/(2d0*pi)**3*amchar(2)
      else
         xchar1dow = 0d0
      endif
c -------------------------------------------------------------------- c
c ------------------- chargino strange strangebar -------------------- c
c -------------------------------------------------------------------- c
      xchar1str = xchar1dow
c -------------------------------------------------------------------- c
c ----------------------- chargino top topbar ------------------------ c
c -------------------------------------------------------------------- c
      xmu1=mt**2/amchar(2)**2
      xmu2=mt**2/amchar(2)**2
      xmu3=amchar(1)**2/amchar(2)**2
      if(amchar(2).gt.(amchar(1)+2d0*mt)) then
         CALL NS_integ2(NS_chartop,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .               xmu3,nx1t,ny1t,sum1)
         xchar1top = 3d0*sum1/64d0/(2d0*pi)**3*amchar(2)
      else
         xchar1top = 0d0
      endif
c -------------------------------------------------------------------- c
c -------------------- chargino bottom bottombar --------------------- c
c -------------------------------------------------------------------- c
      xmu1=mbp**2/amchar(2)**2
      xmu2=mbp**2/amchar(2)**2
      xmu3=amchar(1)**2/amchar(2)**2
      if(amchar(2).gt.(amchar(1)+2d0*mbp)) then
         CALL NS_integ2(NS_charbot,NS_ax,NS_bx,NS_ay,NS_by,xmu1,xmu2,
     .               xmu3,nx1t,ny1t,sum1)

         xchar1bot = 3d0*sum1/64d0/(2d0*pi)**3*amchar(2)
      else
         xchar1bot = 0d0
      endif
c -------------------------------------------------------------------- c
c ------------------------- gluino up downbar ------------------------ c
c -------------------------------------------------------------------- c
      do ni=1,2
         xmu1=0d0
         xmu2=0d0
         xmu3=mgl**2/amchar(ni)**2

         if(amchar(ni).gt.dabs(mgl)) then
            CALL NS_integ2(NS_gluiupdb,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                  xmu2,xmu3,nx1t,ny1t,sum1)
            xgluiupdb(ni) = sum1/64d0/pi**3*amchar(ni)
         else
            xgluiupdb(ni) = 0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c ---------------------- gluino charm strangebar --------------------- c
c -------------------------------------------------------------------- c
      xgluichsb(1) = xgluiupdb(1)
      xgluichsb(2) = xgluiupdb(2)
c -------------------------------------------------------------------- c
c ----------------------- gluino top bottombar ----------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2
         xmu1=mt**2/amchar(ni)**2
         xmu2=mbp**2/amchar(ni)**2
         xmu3=mgl**2/amchar(ni)**2
         if(amchar(ni).gt.(dabs(mgl)+mt+mbp)) then
            CALL NS_integ2(NS_gluitopbb,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                  xmu2,xmu3,nx1t,ny1t,sum1)
            xgluitopbb(ni) = sum1/64d0/pi**3*amchar(ni)
         else
            xgluitopbb(ni) = 0d0
         endif
      enddo
      end
c
c ==================================================================== c
c =========================== neutralino e+ nue ====================== c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_chipmel(x1,x2)
c
      IMPLICIT NONE
      INTEGER ni,nj,i,k
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsel(2),dselb(2)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION ql(5,2),qr(5,2),ol(5,2),or(5,2),ble(2,2)
      DOUBLE PRECISION awff,vwff
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION xmuneut1,x3,X1,X2,y1,y2,y3,xmusel(2),
     .xmusnel(2),chipmsel,xmuw,dw,chipmw,chipmwsel,chipmwsnel
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
c
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/NS_coup18/awff,vwff
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
c
      xmuneut1 = amneut(nj)**2/amchar(ni)**2
c
      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c
      do i=1,2
         ble(1,i) = 0d0
         ble(2,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c                           sfermion exchange
c -------------------------------------------------------------------- c
      xmusel(1)  = ase1**2/amchar(ni)**2
      xmusel(2)  = ase2**2/amchar(ni)**2
      xmusnel(1) = asne1**2/amchar(ni)**2
      xmusnel(2) = asne2**2/amchar(ni)**2
c
      dsel(1)  = 1d0-x1-xmusel(1)
      dsel(2)  = 1d0-x1-xmusel(2)
      dselb(1) = 1d0-x2-xmusnel(1)
      dselb(2) = 1d0-x2-xmusnel(2)
c
      chipmsel=0d0
c
      if(amneut(nj).le.amchar(ni)) then
         do i=1,2
            do k=1,2
               if (xmusel(k).gt.1d0.and.xmusel(i).gt.1d0)then
               chipmsel=chipmsel+g2s**2/dsel(k)/dsel(i)*x1*y1*
     .              (ale(k,ni)*ale(i,ni)+ble(k,ni)*ble(i,ni))*
     .              (ae(k,nj)*ae(i,nj)+be(k,nj)*be(i,nj))
               endif
               if (xmusnel(k).gt.1d0.and.xmusnel(i).gt.1d0)then
               chipmsel=chipmsel
     .              +g2s**2/dselb(k)/dselb(i)*x2*y2*
     .              (alsne(k,ni)*alsne(i,ni)+blsne(k,ni)*blsne(i,ni))*
     .              (ane(k,nj)*ane(i,nj)+bne(k,nj)*bne(i,nj))
               endif
               if (xmusnel(k).gt.1d0.and.xmusel(i).gt.1d0)then
               chipmsel=chipmsel
     .              +g2s**2/dsel(i)/dselb(k)*(
     .              (ane(k,nj)*ale(i,ni)*blsne(k,ni)*be(i,nj)+
     .              alsne(k,ni)*ae(i,nj)*bne(k,nj)*ble(i,ni))*
     .              (-x1*y1-x2*y2+x3*y3)
     .              +2d0*(alsne(k,ni)*ane(k,nj)*ale(i,ni)*ae(i,nj)
     .              +blsne(k,ni)*bne(k,nj)*ble(i,ni)*be(i,nj))*
     .              xmneut(nj)/xmchar(ni)*y3)
               endif
c
            enddo
         enddo
      else
         chipmsel=0d0
      endif
c -------------------------------------------------------------------- c
c                               W+ exchange
c -------------------------------------------------------------------- c
      xmuw  = mw**2/amchar(ni)**2
      dw    = y3-xmuw
      chipmw=0d0
      if(amneut(nj).le.amchar(ni)) then
         if ((mw+amneut(nj)).gt.amchar(ni))then
         chipmw=g2s**2*16d0/dw**2*vwff**2*
     .        (ol(nj,ni)**2*x2*y2+or(nj,ni)**2*x1*y1
     .        -2d0*xmneut(nj)/xmchar(ni)*ol(nj,ni)*or(nj,ni)*y3)
         endif
      else
         chipmw=0d0
      endif
c -------------------------------------------------------------------- c
c                      W+ charged sfermion interference
c -------------------------------------------------------------------- c
      chipmwsel=0d0
c
      if(amneut(nj).le.amchar(ni)) then
         do i=1,2
         if ((mw+amneut(nj)).gt.amchar(ni).and.xmusel(i).gt.1d0)then
            chipmwsel=chipmwsel-g2s**2*4d0/dsel(i)/dw*(
     .           2d0*ale(i,ni)*ae(i,nj)*or(nj,ni)*vwff*x1*y1
     .           -2d0*xmneut(nj)/xmchar(ni)*ale(i,ni)*ae(i,nj)*
     .           ol(nj,ni)*vwff*y3)
         endif
         enddo
      else
         chipmwsel=0d0
      endif
c
c -------------------------------------------------------------------- c
c                         W+ sneutrino interference
c -------------------------------------------------------------------- c
      chipmwsnel=0d0
c
        if(amneut(nj).le.amchar(ni)) then
           do i=1,2
        if ((mw+amneut(nj)).gt.amchar(ni).and.xmusnel(i).gt.1d0)then
              chipmwsnel=chipmwsnel+g2s**2*4d0/dselb(i)/dw*(
     .             2d0*alsne(i,ni)*ane(i,nj)*ol(nj,ni)*vwff*x2*y2
     .             -2d0*xmneut(nj)/xmchar(ni)*alsne(i,ni)*ane(i,nj)*
     .             or(nj,ni)*vwff*y3)
        endif
           enddo
        else
           chipmwsnel=0d0
        endif
c -------------------------------------------------------------------- c
      NS_chipmel = chipmsel+chipmw+chipmwsel+chipmwsnel
c
      end
c
c ==================================================================== c
c =========================== neutralino mu+ numu ==================== c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_chipmmu(x1,x2)
c
      IMPLICIT NONE
      INTEGER ni,nj,i,k
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsel(2),dselb(2)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION ql(5,2),qr(5,2),ol(5,2),or(5,2),blmu(2,2)
      DOUBLE PRECISION awff,vwff
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION xmuneut1,x3,X1,X2,y1,y2,y3,xmusel(2),
     .xmusnel(2),chipmsel,xmuw,dw,chipmw,chipmwsel,chipmwsnel
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
c
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/NS_coup18/awff,vwff
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
c
      xmuneut1 = amneut(nj)**2/amchar(ni)**2
c
      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c
      do i=1,2
         blmu(1,i) = 0d0
         blmu(2,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c                           sfermion exchange
c -------------------------------------------------------------------- c
      xmusel(1)  = asmu1**2/amchar(ni)**2
      xmusel(2)  = asmu2**2/amchar(ni)**2
      xmusnel(1) = asnmu1**2/amchar(ni)**2
      xmusnel(2) = asnmu2**2/amchar(ni)**2
c
      dsel(1)  = 1d0-x1-xmusel(1)
      dsel(2)  = 1d0-x1-xmusel(2)
      dselb(1) = 1d0-x2-xmusnel(1)
      dselb(2) = 1d0-x2-xmusnel(2)
c
      chipmsel=0d0
c
      if(amneut(nj).le.amchar(ni)) then
         do i=1,2
            do k=1,2
               if (xmusel(k).gt.1d0.and.xmusel(i).gt.1d0)then
               chipmsel=chipmsel+g2s**2/dsel(k)/dsel(i)*x1*y1*
     .              (almu(k,ni)*almu(i,ni)+blmu(k,ni)*blmu(i,ni))*
     .              (amu(k,nj)*amu(i,nj)+bmu(k,nj)*bmu(i,nj))
               endif
               if (xmusnel(k).gt.1d0.and.xmusnel(i).gt.1d0)then
               chipmsel=chipmsel
     .              +g2s**2/dselb(k)/dselb(i)*x2*y2*
     .              (alsnm(k,ni)*alsnm(i,ni)+blsnm(k,ni)*blsnm(i,ni))*
     .              (anmu(k,nj)*anmu(i,nj)+bnmu(k,nj)*bnmu(i,nj))
               endif
               if (xmusnel(k).gt.1d0.and.xmusel(i).gt.1d0)then
               chipmsel=chipmsel
     .              +g2s**2/dsel(i)/dselb(k)*(
     .              (anmu(k,nj)*almu(i,ni)*blsnm(k,ni)*bmu(i,nj)+
     .              alsnm(k,ni)*amu(i,nj)*bnmu(k,nj)*blmu(i,ni))*
     .              (-x1*y1-x2*y2+x3*y3)
     .              +2d0*(alsnm(k,ni)*anmu(k,nj)*almu(i,ni)*amu(i,nj)
     .              +blsnm(k,ni)*bnmu(k,nj)*blmu(i,ni)*bmu(i,nj))*
     .              xmneut(nj)/xmchar(ni)*y3)
               endif
c
            enddo
         enddo
      else
         chipmsel=0d0
      endif
c -------------------------------------------------------------------- c
c                               W+ exchange
c -------------------------------------------------------------------- c
      xmuw  = mw**2/amchar(ni)**2
      dw    = y3-xmuw
      chipmw=0d0
      if(amneut(nj).le.amchar(ni)) then
         if ((mw+amneut(nj)).gt.amchar(ni))then
         chipmw=g2s**2*16d0/dw**2*vwff**2*
     .        (ol(nj,ni)**2*x2*y2+or(nj,ni)**2*x1*y1
     .        -2d0*xmneut(nj)/xmchar(ni)*ol(nj,ni)*or(nj,ni)*y3)
         endif
      else
         chipmw=0d0
      endif
c -------------------------------------------------------------------- c
c                      W+ charged sfermion interference
c -------------------------------------------------------------------- c
      chipmwsel=0d0
c
      if(amneut(nj).le.amchar(ni)) then
         do i=1,2
         if ((mw+amneut(nj)).gt.amchar(ni).and.xmusel(i).gt.1d0)then
            chipmwsel=chipmwsel-g2s**2*4d0/dsel(i)/dw*(
     .           2d0*almu(i,ni)*amu(i,nj)*or(nj,ni)*vwff*x1*y1
     .           -2d0*xmneut(nj)/xmchar(ni)*almu(i,ni)*amu(i,nj)*
     .           ol(nj,ni)*vwff*y3)
         endif
         enddo
      else
         chipmwsel=0d0
      endif
c
c -------------------------------------------------------------------- c
c                         W+ sneutrino interference
c -------------------------------------------------------------------- c
      chipmwsnel=0d0
c
        if(amneut(nj).le.amchar(ni)) then
           do i=1,2
        if ((mw+amneut(nj)).gt.amchar(ni).and.xmusnel(i).gt.1d0)then
              chipmwsnel=chipmwsnel+g2s**2*4d0/dselb(i)/dw*(
     .             2d0*alsnm(i,ni)*anmu(i,nj)*ol(nj,ni)*vwff*x2*y2
     .             -2d0*xmneut(nj)/xmchar(ni)*alsnm(i,ni)*anmu(i,nj)*
     .             or(nj,ni)*vwff*y3)
        endif
           enddo
        else
           chipmwsnel=0d0
        endif
c -------------------------------------------------------------------- c
      NS_chipmmu = chipmsel+chipmw+chipmwsel+chipmwsnel
c
      end
c
c ==================================================================== c
c ======================= neutralino tau+ nutau ====================== c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_chipmtau(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,i,K
      DOUBLE PRECISION dsto(2),dstob(2),bltau(2,2)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION ql(5,2),qr(5,2),ol(5,2),or(5,2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION achtop,vchtop,achtau,vchtau,awff,vwff
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION amne,xmustau(2),xmusn(2),xmutau,xmun
     .,xmuneut1,x3,x1,x2,y3,uh,th,chipmstau,xmuw,dw,chipmw,rh,sh,rk,
     .xmuch,dh,chipmh,chipmwstau,chipmhstau,chipmwh
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/NS_coup15/achtop,vchtop,achtau,vchtau
      COMMON/NS_coup18/awff,vwff
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
c
c --- neutrino mass ---
      amne = 0d0
c
      xmustau(1) = astau1**2/amchar(ni)**2
      xmustau(2) = astau2**2/amchar(ni)**2
      xmusn(1)   = asntau1**2/amchar(ni)**2
      xmusn(2)   = asntau2**2/amchar(ni)**2
      xmutau   = mtau**2/amchar(ni)**2
      xmun     = amne**2/amchar(ni)**2
      xmuneut1 = amneut(nj)**2/amchar(ni)**2
c
      x3 = 2d0-x1-x2
      y3 = 1d0+xmuneut1-x3
c
      do i=1,2
         bltau(1,i)   = 0d0
         bltau(2,i)   = 0d0
      enddo
c
      uh = 1d0-x1+xmun
      th = 1d0-x2+xmutau
c -------------------------------------------------------------------- c
c                         stau/sneutrino_tau exchange
c -------------------------------------------------------------------- c
      dsto(1)  = 1-x1-xmustau(1)+xmun
      dsto(2)  = 1-x1-xmustau(2)+xmun
      dstob(1) = 1-x2-xmusn(1)+xmutau
      dstob(2) = 1-x2-xmusn(2)+xmutau
c
      chipmstau=0d0
c
      if((amneut(nj)+mtau).le.amchar(ni)) then
         do i=1,2
            do k=1,2
       if (xmustau(k).gt.1d0.and.xmustau(i).gt.1d0) then
               chipmstau=chipmstau
     .          +g2s**2/dsto(k)/dsto(i)*(
     .           (altau(i,ni)*bltau(k,ni)+bltau(i,ni)*altau(k,ni))*
     .           (atau(i,nj)*btau(k,nj)+btau(i,nj)*atau(k,nj))*
     .           xmneut(nj)/xmchar(ni)*dsqrt(xmun*xmutau)*(-4d0)+
     .           (altau(i,ni)*altau(k,ni)+bltau(i,ni)*bltau(k,ni))*
     .           (atau(i,nj)*btau(k,nj)+btau(i,nj)*atau(k,nj))*
     .           xmneut(nj)/xmchar(ni)*dsqrt(xmutau)*2d0*
     .           (uh-xmun-1d0)+
     .           (altau(i,ni)*bltau(k,ni)+bltau(i,ni)*altau(k,ni))*
     .           (atau(i,nj)*atau(k,nj)+btau(i,nj)*btau(k,nj))*
     .           dsqrt(xmun)*2d0*(uh-xmutau-xmuneut1)+
     .           (altau(i,ni)*altau(k,ni)+bltau(i,ni)*bltau(k,ni))*
     .           (atau(i,nj)*atau(k,nj)+btau(i,nj)*btau(k,nj))*
     .           (-uh**2+uh*(1d0+xmuneut1+xmun+xmutau)-
     .           (xmuneut1+xmutau)*(1d0+xmun)))
       endif
      if (xmusn(k).gt.1d0.and.xmusn(i).gt.1d0) then
               chipmstau=chipmstau
     .           +g2s**2/dstob(k)/dstob(i)*(
     .           (alsnt(i,ni)*blsnt(k,ni)+blsnt(i,ni)*alsnt(k,ni))*
     .           (antau(i,nj)*bntau(k,nj)+bntau(i,nj)*antau(k,nj))*
     .           xmneut(nj)/xmchar(ni)*dsqrt(xmutau*xmun)*(-4d0)+
     .           (alsnt(i,ni)*alsnt(k,ni)+blsnt(i,ni)*blsnt(k,ni))*
     .           (antau(i,nj)*bntau(k,nj)+bntau(i,nj)*antau(k,nj))*
     .           xmneut(nj)/xmchar(ni)*dsqrt(xmun)*2d0*
     .           (th-xmutau-1d0)+
     .           (alsnt(i,ni)*blsnt(k,ni)+blsnt(i,ni)*alsnt(k,ni))*
     .           (antau(i,nj)*antau(k,nj)+bntau(i,nj)*bntau(k,nj))*
     .           dsqrt(xmutau)*2d0*(th-xmun-xmuneut1)+
     .           (alsnt(i,ni)*alsnt(k,ni)+blsnt(i,ni)*blsnt(k,ni))*
     .           (antau(i,nj)*antau(k,nj)+bntau(i,nj)*bntau(k,nj))*
     .           (-th**2+th*(1d0+xmuneut1+xmun+xmutau)-(xmuneut1+xmun)*
     .           (1d0+xmutau)))
       endif
       if (xmustau(k).gt.1d0.and.xmusn(i).gt.1d0)then
          chipmstau=chipmstau
     .           -2d0*g2s**2/dsto(k)/dstob(i)*(
     .           (blsnt(i,ni)*bltau(k,ni)*antau(i,nj)*atau(k,nj)
     .           +alsnt(i,ni)*altau(k,ni)*bntau(i,nj)*btau(k,nj))*
     .           xmneut(nj)/xmchar(ni)*dsqrt(xmun*xmutau)*(-2d0)+
     .           (alsnt(i,ni)*bltau(k,ni)*antau(i,nj)*atau(k,nj)
     .           +blsnt(i,ni)*altau(k,ni)*bntau(i,nj)*btau(k,nj))*
     .           xmneut(nj)/xmchar(ni)*dsqrt(xmun)*
     .           (th-xmutau-1d0)+
     .           (blsnt(i,ni)*altau(k,ni)*antau(i,nj)*atau(k,nj)
     .           +alsnt(i,ni)*bltau(k,ni)*bntau(i,nj)*btau(k,nj))*
     .           xmneut(nj)/xmchar(ni)*dsqrt(xmutau)*
     .           (uh-xmun-1d0)+
     .           (alsnt(i,ni)*altau(k,ni)*antau(i,nj)*atau(k,nj)
     .           +blsnt(i,ni)*bltau(k,ni)*bntau(i,nj)*btau(k,nj))*
     .           xmneut(nj)/xmchar(ni)*(uh+th-xmuneut1-1d0)+
     .           (altau(k,ni)*blsnt(i,ni)*atau(k,nj)*bntau(i,nj)
     .           +bltau(k,ni)*alsnt(i,ni)*btau(k,nj)*antau(i,nj))*
     .           dsqrt(xmun*xmutau)*(uh+th-xmun-xmutau)+
     .           (altau(k,ni)*alsnt(i,ni)*atau(k,nj)*bntau(i,nj)
     .           +bltau(k,ni)*blsnt(i,ni)*btau(k,nj)*antau(i,nj))*
     .           dsqrt(xmun)*(uh-xmutau-xmuneut1)+
     .           (bltau(k,ni)*blsnt(i,ni)*atau(k,nj)*bntau(i,nj)
     .           +altau(k,ni)*alsnt(i,ni)*btau(k,nj)*antau(i,nj))*
     .           dsqrt(xmutau)*(th-xmun-xmuneut1)+
     .           (bltau(k,ni)*alsnt(i,ni)*atau(k,nj)*bntau(i,nj)
     .           +altau(k,ni)*blsnt(i,ni)*btau(k,nj)*antau(i,nj))*
     .           (uh*th-xmun*xmutau-xmuneut1))
         endif
            enddo
         enddo
      else
         chipmstau=0d0
      endif
c
c -------------------------------------------------------------------- c
c                              W exchange
c -------------------------------------------------------------------- c
      xmuw = mw**2/amchar(ni)**2
      dw   = y3-xmuw
c
      chipmw = 0d0
c
      rh = xmuneut1+xmun+xmutau-th-uh+1d0
      sh = (xmuneut1-th-uh+1d0)*(xmun+xmutau)+4d0*xmun*xmutau
      rk = xmuneut1*(xmun+xmutau-th-uh+4d0)+xmun+xmutau-uh-th
c
      if((amneut(nj)+mtau).le.amchar(ni)) then
         if ((mw+amneut(nj)).gt.amchar(ni))then
         chipmw=chipmw+g2s**2/dw**2*(
     .    ol(nj,ni)*or(nj,ni)*2d0*vwff**2*
     .    xmneut(nj)/xmchar(ni)*(8d0/xmuw**2*rh*sh-16d0/xmuw*sh
     .    -16d0*(xmuneut1-uh-th+1d0))
     .    +(ol(nj,ni)**2+or(nj,ni)**2)*2d0*vwff**2*
     .    (-2d0/xmuw**2*rk*sh+8d0/xmuw*(xmuneut1*(2d0*xmun*xmutau+
     .    2d0*(xmun+xmutau)-xmun*th-xmutau*uh)+2d0*xmun*xmutau
     .    -xmun*uh-xmutau*th)+4d0*(xmuneut1*(uh+th-xmun-xmutau-2d0)+
     .    (xmun+xmutau)*(uh+th-1d0)-2d0*xmun*xmutau+th*(-th+1d0)+
     .    uh*(-uh+1d0)))+
     .    (ol(nj,ni)**2-or(nj,ni)**2)*vwff**2*8d0*(
     .    xmuneut1*(xmun-xmutau+th-uh)+(xmun+xmutau)*(th-uh)-xmun+
     .    xmutau+th*(-th+1d0)+uh*(uh-1d0))
     .    )
         endif
      else
         chipmw=0d0
      endif
c -------------------------------------------------------------------- c
c                              H+ exchange
c -------------------------------------------------------------------- c
      xmuch = cmass**2/amchar(ni)**2
      dh    = y3-xmuch
      chipmh=0d0
      if ((cmass+amneut(nj)).gt.amchar(ni)) then
      if((amneut(nj)+mtau).le.amchar(ni)) then
         chipmh=g2s**2/dh**2*(
     .    ql(nj,ni)*qr(nj,ni)*(
     .    (vchtau**2-achtau**2)*
     .    xmneut(nj)/xmchar(ni)*dsqrt(xmun*xmutau)*(-16d0)+
     .    (vchtau**2+achtau**2)*
     .    xmneut(nj)/xmchar(ni)*8d0*(1d0+xmuneut1-th-uh) )
     .    +(ql(nj,ni)**2+qr(nj,ni)**2)*(
     .    (vchtau**2-achtau**2)*
     .    dsqrt(xmun*xmutau)*4d0*(xmun+xmutau-th-uh)+
     .    (vchtau**2+achtau**2)*
     .    2d0*(xmuneut1*(uh+th-xmun-xmutau)+(xmun+xmutau)*(uh+th-1d0)
     .    -(th+uh)**2+uh+th)) )
      else
         chipmh=0d0
      endif
      endif
c -------------------------------------------------------------------- c
c                    W+ stau/sneutrino_tau interference
c -------------------------------------------------------------------- c
      chipmwstau=0d0
c
      if((amneut(nj)+mtau).le.amchar(ni)) then
         do i=1,2
         if ((mw+amneut(nj)).gt.amchar(ni).and.xmustau(i).gt.1.0d0)then
            chipmwstau=chipmwstau
     .       +g2s**2/dsto(i)/dw*2d0*vwff*(
     .      bltau(i,ni)*btau(i,nj)*or(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*dsqrt(xmun*xmutau)*(
     .      1/xmuw*(-4d0)*(1d0+xmuneut1+xmun+xmutau-uh-th)+16d0) +
     .      altau(i,ni)*atau(i,nj)*ol(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*(2d0/xmuw*((xmuneut1+1d0-uh-th)*
     .      (xmun+xmutau)+4d0*xmun*xmutau)+4d0*(1d0+xmuneut1-
     .      uh-th)) +
     .      bltau(i,ni)*btau(i,nj)*ol(nj,ni)*dsqrt(xmun*xmutau)*
     .      (2d0/xmuw*(xmuneut1*(xmun+xmutau-th-uh+4d0)+xmun+xmutau
     .      -th-uh)+4d0*(xmun+xmutau-th-uh)) +
     .      altau(i,ni)*atau(i,nj)*or(nj,ni)*
     .      (2d0/xmuw*(xmuneut1*(-2d0*xmun*xmutau+xmun*th
     .      -2d0*xmun+xmutau*uh-2d0*xmutau)+xmun*(-2d0*xmutau+uh)
     .      +xmutau*th)+4d0*(xmuneut1*(xmun-uh+1d0)+xmun*(xmutau-uh)
     .      +xmutau*(1d0-uh)+uh**2-uh)) +
     .      atau(i,nj)*bltau(i,ni)*ol(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*dsqrt(xmun)*(2d0/xmuw*
     .      (xmuneut1*(1d0+xmun-uh)+1d0+xmun*(xmutau-uh)+xmutau*
     .      (xmutau-th-2d0*uh+3d0)+th*(uh-1d0)+uh*(uh-2d0))+
     .      4d0*(1d0+xmutau-th)) +
     .      btau(i,nj)*altau(i,ni)*or(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*dsqrt(xmutau)*
     .      (2d0/xmuw*(xmuneut1*(-1d0-xmun+uh)-1d0+xmutau*(-1d0+uh)
     .      +uh*(2d0-th-uh)+th+xmun*(-xmun-xmutau+th+2d0*uh-2d0))+
     .      8d0*(1d0+xmun-uh)) +
     .      bltau(i,ni)*atau(i,nj)*or(nj,ni)*
     .      dsqrt(xmun)*((-2d0)/xmuw*(xmuneut1-uh+xmutau)*
     .      (xmuneut1+xmun+xmutau-th-uh+1d0)+8d0*(xmuneut1+xmutau
     .      -uh)) +
     .      altau(i,ni)*btau(i,nj)*ol(nj,ni)*
     .      dsqrt(xmutau)*(2d0/xmuw*(xmuneut1*(xmuneut1+
     .      3d0*xmun-th-2d0*uh+1d0)+xmun*(xmun+xmutau-th-2d0*uh)+
     .      xmutau*(1d0-uh)+uh*th+uh**2-uh)+4d0*(xmuneut1+xmun-th)) )
            endif
           if ((mw+amneut(nj)).gt.amchar(ni).and.xmusn(i).gt.1.0d0)then
            chipmwstau=chipmwstau
     .      -g2s**2/dstob(i)/dw*2d0*vwff*(
     .      blsnt(i,ni)*bntau(i,nj)*ol(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*dsqrt(xmun*xmutau)*(
     .      1/xmuw*(-4d0)*(1d0+xmuneut1+xmun+xmutau-uh-th)+16d0) +
     .      alsnt(i,ni)*antau(i,nj)*or(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*(2d0/xmuw*((xmun+xmutau)*(xmuneut1
     .      +1d0-th-uh)+4d0*xmun*xmutau)+4d0*(1d0+xmuneut1-uh
     .      -th)) +
     .      blsnt(i,ni)*bntau(i,nj)*or(nj,ni)*dsqrt(xmun*xmutau)*
     .      (2d0/xmuw*(xmuneut1*(xmun+xmutau-th-uh+4d0)+xmun+xmutau
     .      -th-uh)+4d0*(xmun+xmutau-th-uh)) +
     .      alsnt(i,ni)*antau(i,nj)*ol(nj,ni)*
     .      (2d0/xmuw*(xmuneut1*(-2d0*xmun*xmutau+xmun*th-2d0*xmun
     .      +xmutau*uh-2d0*xmutau)+xmun*(-2d0*xmutau+uh)+xmutau*th)
     .      +4d0*(xmuneut1*(xmutau-th+1d0)+xmun*(xmutau-th+1d0)+
     .      th*(-xmutau+th-1d0))) +
     .      antau(i,nj)*blsnt(i,ni)*or(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*dsqrt(xmutau)*
     .      (2d0/xmuw*(xmuneut1*(1d0+xmutau-th)+1d0+xmun*(xmun+
     .      xmutau-2d0*th-uh+3d0)+th*(th-xmutau+uh-2d0)-uh)
     .      +4d0*(1d0+xmun-uh)) +
     .      alsnt(i,ni)*bntau(i,nj)*ol(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*dsqrt(xmun)*
     .      (2d0/xmuw*(xmuneut1*(-1d0-xmutau+th)-1d0+xmun*(-xmutau+
     .      th-1d0)+xmutau*(2d0*th+uh-2d0-xmutau)-th*(th+uh)+2d0*
     .      th+uh)+8d0*(1d0+xmutau-th)) +
     .      blsnt(i,ni)*antau(i,nj)*ol(nj,ni)*dsqrt(xmutau)*(
     .      2d0/xmuw*(xmuneut1*(-xmuneut1-2d0*xmun-xmutau+2d0*th+uh
     .      -1d0)+xmun*(-xmun-xmutau+2d0*th+uh-1d0)+th*(xmutau-th-
     .      uh+1))+8d0*(xmuneut1+xmun-th)) +
     .      alsnt(i,ni)*bntau(i,nj)*or(nj,ni)*dsqrt(xmun)*(
     .      2d0/xmuw*(xmuneut1*(xmuneut1+3d0*xmutau-uh-2d0*th+1d0)
     .      +xmun*(xmutau-th+1d0)+xmutau*(xmutau-2d0*th-uh)+uh*th+
     .      th**2-th)+4d0*(xmuneut1+xmutau-uh)) )
            endif
         enddo
      else
         chipmwstau=0d0
      endif
c -------------------------------------------------------------------- c
c                     H+ stau/sneutrino_tau interference
c -------------------------------------------------------------------- c
      chipmhstau = 0d0
c
      if((amneut(nj)+mtau).le.amchar(ni)) then
         do i=1,2
      if ((cmass+amneut(nj)).gt.amchar(ni).and.xmustau(i).gt.1.0d0)then
            chipmhstau=chipmhstau
     .       -g2s**2/dh/dsto(i)*2d0*vchtau*(
     .       btau(i,nj)*bltau(i,ni)*qr(nj,ni)*(-2d0)*
     .       dsqrt(xmun)*(xmuneut1+xmutau-uh) +
     .       atau(i,nj)*altau(i,ni)*ql(nj,ni)*2d0*
     .       dsqrt(xmutau)*(xmuneut1+xmun-th) +
     .       btau(i,nj)*bltau(i,ni)*ql(nj,ni)*2d0*
     .       xmneut(nj)/xmchar(ni)*dsqrt(xmun)*(1d0+xmutau-th) +
     .       atau(i,nj)*altau(i,ni)*qr(nj,ni)*(-2d0)*
     .       xmneut(nj)/xmchar(ni)*dsqrt(xmutau)*(1d0+xmun-uh) +
     .       btau(i,nj)*altau(i,ni)*ql(nj,ni)*2d0*
     .       xmneut(nj)/xmchar(ni)*(1d0+xmuneut1-uh-th) +
     .       atau(i,nj)*bltau(i,ni)*qr(nj,ni)*(-4d0)*
     .       xmneut(nj)/xmchar(ni)*dsqrt(xmun*xmutau) +
     .       atau(i,nj)*bltau(i,ni)*ql(nj,ni)*2d0*
     .       dsqrt(xmun*xmutau)*(xmun+xmutau-th-uh) +
     .       altau(i,ni)*btau(i,nj)*qr(nj,ni)*2d0*
     .       (-uh**2-uh*th+uh*(1d0+xmun+xmutau+xmuneut1)-xmutau-
     .       xmuneut1*xmun) )
            endif
      if ((cmass+amneut(nj)).gt.amchar(ni).and.xmusn(i).gt.1.0d0)then
               chipmhstau=chipmhstau
     .       +2d0*g2s**2/dh/dstob(i)*2d0*vchtau*(
     .       bntau(i,nj)*blsnt(i,ni)*ql(nj,ni)*
     .       dsqrt(xmun)*(uh-xmuneut1-xmutau) +
     .       antau(i,nj)*alsnt(i,ni)*qr(nj,ni)*
     .       dsqrt(xmutau)*(-th+xmun+xmuneut1) +
     .       bntau(i,nj)*blsnt(i,ni)*qr(nj,ni)*
     .       dsqrt(xmun)*xmneut(nj)/xmchar(ni)*(-th+1d0+xmutau) +
     .       antau(i,nj)*alsnt(i,ni)*ql(nj,ni)*
     .       dsqrt(xmutau)*xmneut(nj)/xmchar(ni)*(uh-1d0-xmun) +
     .       bntau(i,nj)*alsnt(i,ni)*qr(nj,ni)*
     .       2d0*xmneut(nj)/xmchar(ni)*dsqrt(xmun*xmutau) +
     .       antau(i,nj)*blsnt(i,ni)*qr(nj,ni)*
     .       (uh*th+th**2-th*(1d0+xmun+xmutau+xmuneut1)+xmun+
     .       xmutau*xmuneut1) +
     .       antau(i,nj)*blsnt(i,ni)*ql(nj,ni)*
     .       xmneut(nj)/xmchar(ni)*(uh+th-xmuneut1-1d0) +
     .       bntau(i,nj)*alsnt(i,ni)*ql(nj,ni)*
     .       dsqrt(xmun*xmutau)*(uh+th-xmun-xmutau) )
            endif
         enddo
      else
         chipmhstau=0d0
      endif
c
c -------------------------------------------------------------------- c
c 	               interference W+ H-
c -------------------------------------------------------------------- c
      chipmwh=0d0
c
      if((amneut(nj)+mtau).le.amchar(ni)) then
      if ((mw+amneut(nj)).gt.amchar(ni).and.
     .(cmass+amneut(nj)).gt.amchar(ni))then
         chipmwh=chipmwh-4d0*g2s**2/dh/dw*vwff*vchtau*(
     .    (ol(nj,ni)*ql(nj,ni)+or(nj,ni)*qr(nj,ni))*
     .    xmneut(nj)/xmchar(ni)*dsqrt(xmutau)*(
     .    2d0/xmuw*(-2d0-xmuneut1*(xmun+xmutau-th-uh+2d0)-xmun*
     .    (2d0*xmun+xmutau-3d0*th-3d0*uh+5d0)-xmutau*(xmun-th-uh
     .    +1d0)-(th+uh)**2+3d0*(th+uh))+4d0*(-uh+xmun+1d0))
     .    +(ol(nj,ni)*qr(nj,ni)+or(nj,ni)*ql(nj,ni))*dsqrt(xmutau)*(
     .    2d0/xmuw*(2d0*xmuneut1**2-xmuneut1*(-5d0*xmun-xmutau+
     .    3d0*(th+uh)-2d0)-xmun*(-2d0*xmun-xmutau+3d0*(th+uh)-1d0)-
     .    xmutau*(-xmun+uh+th-1d0)+(th+uh)**2-(th+uh))
     .    -4d0*(xmuneut1+xmun-th) ) )
         endif
      else
         chipmwh=0d0
      endif
c
c -------------------------------------------------------------------- c
c
      NS_chipmtau = chipmstau+chipmw+chipmh+chipmwstau+chipmhstau+
     .              chipmwh
c
      end
c
c ==================================================================== c
c ========================= neutralino up downbar ==================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_chipmup(x1,x2)
c
      IMPLICIT NONE
      INTEGER ni,nj,i,K
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsdown(2),dsup(2),bldo(2,2),blup(2,2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION aup(2,5),bup(2,5),ado(2,5),bdo(2,5)
      DOUBLE PRECISION alup(2,2),aldo(2,2)
      DOUBLE PRECISION ql(5,2),qr(5,2),ol(5,2),or(5,2)
      DOUBLE PRECISION awff,vwff
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION xmuneut1,x3,x1,x2,Y1,Y2,Y3,xmusd(2),
     .chipmsdown,xmusup(2),chipmsup,chipmsdownsup,xmuw,dw
     .,chipmw,chipmwsdown,chipmwsup
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
c
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup7/alup,aldo
      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/NS_coup10/aup,bup,ado,bdo
      COMMON/NS_coup18/awff,vwff
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
c
      xmuneut1 = amneut(nj)**2/amchar(ni)**2
c
      x3=2d0-x1-x2
      y1=1d0-xmuneut1-x1
      y2=1d0-xmuneut1-x2
      y3=1d0+xmuneut1-x3
c
      do i=1,2
         bldo(1,i) = 0d0
         bldo(2,i) = 0d0
         blup(1,i) = 0d0
         blup(2,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c                            sdown exchange
c -------------------------------------------------------------------- c
      xmusd(1) = asdown1**2/amchar(ni)**2
      xmusd(2) = asdown2**2/amchar(ni)**2
c
      dsdown(1)=1d0-x1-xmusd(1)
      dsdown(2)=1d0-x1-xmusd(2)
c
      chipmsdown=0d0
c
      if(amneut(nj).le.amchar(ni)) then
         do i=1,2
            do k=1,2
               if (xmusd(k).gt.1.0d0.and.xmusd(i).gt.1.0d0)then
               chipmsdown=chipmsdown+g2s**2/dsdown(k)/dsdown(i)*x1*y1*
     .              (aldo(i,ni)*aldo(k,ni)+bldo(i,ni)*bldo(k,ni))
     .              *(ado(i,nj)*ado(k,nj)+bdo(i,nj)*bdo(k,nj))
               endif
            enddo
         enddo
      else
         chipmsdown=0d0
      endif
c -------------------------------------------------------------------- c
c                             sup exchange
c -------------------------------------------------------------------- c
      xmusup(1) = asup1**2/amchar(ni)**2
      xmusup(2) = asup2**2/amchar(ni)**2
c
      dsup(1)=1d0-x2-xmusup(1)
      dsup(2)=1d0-x2-xmusup(2)
c
      chipmsup=0d0
c
      if(amneut(nj).le.amchar(ni)) then
         do i=1,2
            do k=1,2
               if (xmusup(k).gt.1.0d0.and.xmusup(i).gt.1.0d0)then
               chipmsup=chipmsup+g2s**2/dsup(i)/dsup(k)*x2*y2*
     .              (alup(i,ni)*alup(k,ni)+blup(i,ni)*blup(k,ni))
     .              *(aup(i,nj)*aup(k,nj)+bup(i,nj)*bup(k,nj))
               endif
            enddo
         enddo
      else
         chipmsup=0d0
      endif
c -------------------------------------------------------------------- c
c                        sup sdown interference
c -------------------------------------------------------------------- c
      chipmsdownsup=0d0
c
      if(amneut(nj).le.amchar(ni)) then
         do i=1,2
            do k=1,2
               if (xmusup(k).gt.1.0d0.and.xmusd(i).gt.1.0d0)then
               chipmsdownsup=chipmsdownsup+
     .              g2s**2/dsup(k)/dsdown(i)*(
     .              (aup(k,nj)*aldo(i,ni)*blup(k,ni)*bdo(i,nj)+
     .              alup(k,ni)*ado(i,nj)*bup(k,nj)*bldo(i,ni))*
     .              (-x1*y1-x2*y2+x3*y3)
     .              +2d0*(alup(k,ni)*aup(k,nj)*aldo(i,ni)*ado(i,nj)
     .              +blup(k,ni)*bup(k,nj)*bldo(i,ni)*bdo(i,nj))
     .              *xmneut(nj)/xmchar(ni)*y3)
               endif
            enddo
         enddo
      else
         chipmsdownsup=0d0
      endif
c
c -------------------------------------------------------------------- c
c                             W+ exchange
c -------------------------------------------------------------------- c
      xmuw = mw**2/amchar(ni)**2
      dw   = y3-xmuw
c
      chipmw=0d0
c
      if(amneut(nj).le.amchar(ni)) then
         if ((mw+amneut(nj)).gt.amchar(ni))then
         chipmw=g2s**2*16d0/dw**2*vwff**2*
     .        (or(nj,ni)**2*x1*y1+ol(nj,ni)**2*x2*y2
     .        -2d0*ol(nj,ni)*or(nj,ni)*xmneut(nj)/xmchar(ni)*y3)
         endif
      else
         chipmw=0d0
      endif
c -------------------------------------------------------------------- c
c                         W+ sdown interference
c -------------------------------------------------------------------- c
      chipmwsdown=0d0
c
      if(amneut(nj).le.amchar(ni)) then
         do i=1,2
       if ((mw+amneut(nj)).gt.amchar(ni).and.xmusd(i).gt.1.0d0)then
            chipmwsdown=chipmwsdown-g2s**2*4d0/dsdown(i)/dw*
     .           (2d0*vwff*or(nj,ni)*aldo(i,ni)*ado(i,nj)*x1*y1
     .           -2d0*vwff*ol(nj,ni)*aldo(i,ni)*ado(i,nj)*
     .           xmneut(nj)/xmchar(ni)*y3)
       endif
         enddo
      else
         chipmwsdown=0d0
      endif
c
c -------------------------------------------------------------------- c
c                          W+ sup interference
c -------------------------------------------------------------------- c
      chipmwsup=0d0
c
      if(amneut(nj).le.amchar(ni)) then
         do i=1,2
       if ((mw+amneut(nj)).gt.amchar(ni).and.xmusup(i).gt.1.0d0)then
            chipmwsup=chipmwsup+g2s**2*4d0/dsup(i)/dw*(
     .           (2d0*vwff*ol(nj,ni)*alup(i,ni)*aup(i,nj)*x2*y2
     .           -2d0*vwff*or(nj,ni)*alup(i,ni)*aup(i,nj)*
     .           xmneut(nj)/xmchar(ni)*y3))
       endif
         enddo
      else
         chipmwsup=0d0
      endif
c
c -------------------------------------------------------------------- c
      NS_chipmup = chipmsdown+chipmsup+chipmsdownsup+chipmw+
     .             chipmwsdown+chipmwsup
c
      end
c
c ==================================================================== c
c ======================= neutralino top bottombar =================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_chipmtop(x1,x2)
c
      IMPLICIT NONE
      INTEGER ni,nj,i,k
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION ql(5,2),qr(5,2),ol(5,2),or(5,2)
      DOUBLE PRECISION dsb(2),dst(2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2),alstor(2,2),akstor(2,2)
      DOUBLE PRECISION abot(2,5),bbot(2,5),atopr(2,5),btopr(2,5)
      DOUBLE PRECISION achtop,vchtop,achtau,vchtau
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION VCH,ach,xmuneut1,xmusbot(2),xmustop(2),
     .xmut,xmub,x3,x1,x2,y3,uh,th,chipmsfer,
     .xmuw,dw,chipmw,rh,sh,rk,xmuch,dh,chipmh,chipmwsbot,chipmhsbot,
     .chipmwh,gmst(2),gmsb(2)
      DOUBLE PRECISION awff,vwff
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
c
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/NS_coup15/achtop,vchtop,achtau,vchtau
      COMMON/NS_coup18/awff,vwff
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor
      COMMON/NS_neutsbotbot/abot,bbot
      COMMON/NS_neutstoptop/atopr,btopr
c
      vch = vchtop
      ach = achtop
      gmst(1) = ast1
      gmst(2) = ast2
      gmsb(1) = asb1
      gmsb(2) = asb2
      xmuneut1 = amneut(nj)**2/amchar(ni)**2
      xmusbot(1) = asb1**2/amchar(ni)**2
      xmusbot(2) = asb2**2/amchar(ni)**2
      xmustop(1) = ast1**2/amchar(ni)**2
      xmustop(2) = ast2**2/amchar(ni)**2
      xmut     = mt**2/amchar(ni)**2
      xmub     = mbp**2/amchar(ni)**2
c
      x3 = 2d0-x1-x2
      y3 = 1d0+xmuneut1-x3
c
      uh = 1d0-x1+xmut
      th = 1d0-x2+xmub
c -------------------------------------------------------------------- c
c                      sbottom and stop exchange
c -------------------------------------------------------------------- c
      dsb(1) = 1d0-x1-xmusbot(1)+xmut
      dsb(2) = 1d0-x1-xmusbot(2)+xmut
      dst(1) = 1d0-x2-xmustop(1)+xmub
      dst(2) = 1d0-x2-xmustop(2)+xmub
c
      chipmsfer=0d0
c
      if ((amneut(nj)+mbp+mt).le.amchar(ni)) then
         do i=1,2
            do k=1,2
       if ((gmsb(k)+mt).gt.amchar(ni).and.
     .(gmsb(i)+mt).gt.amchar(ni))then
               chipmsfer=chipmsfer
     .           +g2s**2/dsb(k)/dsb(i)*(
     .           (alsbot(i,ni)*aksbot(k,ni)+aksbot(i,ni)*alsbot(k,ni))*
     .           (abot(i,nj)*bbot(k,nj)+bbot(i,nj)*abot(k,nj))*
     .           xmneut(nj)/xmchar(ni)*dsqrt(xmut*xmub)*(-4d0)+
     .           (alsbot(i,ni)*alsbot(k,ni)+aksbot(i,ni)*aksbot(k,ni))*
     .           (abot(i,nj)*bbot(k,nj)+bbot(i,nj)*abot(k,nj))*
     .           xmneut(nj)/xmchar(ni)*dsqrt(xmub)*2d0*
     .           (uh-xmut-1d0)+
     .           (alsbot(i,ni)*aksbot(k,ni)+aksbot(i,ni)*alsbot(k,ni))*
     .           (abot(i,nj)*abot(k,nj)+bbot(i,nj)*bbot(k,nj))*
     .           dsqrt(xmut)*2d0*(uh-xmub-xmuneut1)+
     .           (alsbot(i,ni)*alsbot(k,ni)+aksbot(i,ni)*aksbot(k,ni))*
     .           (abot(i,nj)*abot(k,nj)+bbot(i,nj)*bbot(k,nj))*
     .           (-uh**2+uh*(1d0+xmuneut1+xmut+xmub)-
     .           (xmuneut1+xmub)*(1d0+xmut)))
       endif
       if (xmustop(k).gt.1d0.and.xmustop(i).gt.1d0)then
               chipmsfer=chipmsfer
     .           +g2s**2/dst(k)/dst(i)*(
     .           (alstor(i,ni)*akstor(k,ni)+akstor(i,ni)*alstor(k,ni))*
     .           (atopr(i,nj)*btopr(k,nj)+btopr(i,nj)*atopr(k,nj))*
     .           xmneut(nj)/xmchar(ni)*dsqrt(xmub*xmut)*(-4d0)+
     .           (alstor(i,ni)*alstor(k,ni)+akstor(i,ni)*akstor(k,ni))*
     .           (atopr(i,nj)*btopr(k,nj)+btopr(i,nj)*atopr(k,nj))*
     .           xmneut(nj)/xmchar(ni)*dsqrt(xmut)*2d0*
     .           (th-xmub-1d0)+
     .           (alstor(i,ni)*akstor(k,ni)+akstor(i,ni)*alstor(k,ni))*
     .           (atopr(i,nj)*atopr(k,nj)+btopr(i,nj)*btopr(k,nj))*
     .           dsqrt(xmub)*2d0*(th-xmut-xmuneut1)+
     .           (alstor(i,ni)*alstor(k,ni)+akstor(i,ni)*akstor(k,ni))*
     .           (atopr(i,nj)*atopr(k,nj)+btopr(i,nj)*btopr(k,nj))*
     .           (-th**2+th*(1d0+xmuneut1+xmut+xmub)-(xmuneut1+xmut)*
     .           (1d0+xmub)))
        endif
        if ((gmsb(k)+mt).gt.amchar(ni).and.xmustop(i).gt.1d0)then
                chipmsfer=chipmsfer
     .           -2d0*g2s**2/dsb(k)/dst(i)*(
     .           (akstor(i,ni)*aksbot(k,ni)*atopr(i,nj)*abot(k,nj)
     .           +alstor(i,ni)*alsbot(k,ni)*btopr(i,nj)*bbot(k,nj))*
     .           xmneut(nj)/xmchar(ni)*dsqrt(xmut*xmub)*(-2d0)+
     .           (alstor(i,ni)*aksbot(k,ni)*atopr(i,nj)*abot(k,nj)
     .           +akstor(i,ni)*alsbot(k,ni)*btopr(i,nj)*bbot(k,nj))*
     .           xmneut(nj)/xmchar(ni)*dsqrt(xmut)*
     .           (th-xmub-1d0)+
     .           (akstor(i,ni)*alsbot(k,ni)*atopr(i,nj)*abot(k,nj)
     .           +alstor(i,ni)*aksbot(k,ni)*btopr(i,nj)*bbot(k,nj))*
     .           xmneut(nj)/xmchar(ni)*dsqrt(xmub)*
     .           (uh-xmut-1d0)+
     .           (alstor(i,ni)*alsbot(k,ni)*atopr(i,nj)*abot(k,nj)
     .           +akstor(i,ni)*aksbot(k,ni)*btopr(i,nj)*bbot(k,nj))*
     .           xmneut(nj)/xmchar(ni)*(uh+th-xmuneut1-1d0)+
     .           (alsbot(k,ni)*akstor(i,ni)*abot(k,nj)*btopr(i,nj)
     .           +aksbot(k,ni)*alstor(i,ni)*bbot(k,nj)*atopr(i,nj))*
     .           dsqrt(xmut*xmub)*(uh+th-xmut-xmub)+
     .           (alsbot(k,ni)*alstor(i,ni)*abot(k,nj)*btopr(i,nj)
     .           +aksbot(k,ni)*akstor(i,ni)*bbot(k,nj)*atopr(i,nj))*
     .           dsqrt(xmut)*(uh-xmub-xmuneut1)+
     .           (aksbot(k,ni)*akstor(i,ni)*abot(k,nj)*btopr(i,nj)
     .           +alsbot(k,ni)*alstor(i,ni)*bbot(k,nj)*atopr(i,nj))*
     .           dsqrt(xmub)*(th-xmut-xmuneut1)+
     .           (aksbot(k,ni)*alstor(i,ni)*abot(k,nj)*btopr(i,nj)
     .           +alsbot(k,ni)*akstor(i,ni)*bbot(k,nj)*atopr(i,nj))*
     .           (uh*th-xmut*xmub-xmuneut1))
       endif
            enddo
         enddo
      else
         chipmsfer=0d0
      endif
c
c -------------------------------------------------------------------- c
c                            W exchange
c -------------------------------------------------------------------- c
      xmuw = mw**2/amchar(ni)**2
      dw   = y3-xmuw
c
      chipmw = 0d0
c
      rh = xmuneut1+xmut+xmub-th-uh+1d0
      sh = (xmuneut1-th-uh+1d0)*(xmut+xmub)+4d0*xmut*xmub
      rk = xmuneut1*(xmut+xmub-th-uh+4d0)+xmut+xmub-uh-th
c
      if ((amneut(nj)+mt+mbp).le.amchar(ni)) then
         if ((mw+amneut(nj)).gt.amchar(ni))then
         chipmw=chipmw+g2s**2/dw**2*(
     .    ol(nj,ni)*or(nj,ni)*2d0*vwff**2*
     .    xmneut(nj)/xmchar(ni)*(8d0/xmuw**2*rh*sh-16d0/xmuw*sh
     .    -16d0*(xmuneut1-uh-th+1d0))
     .    +(ol(nj,ni)**2+or(nj,ni)**2)*2d0*vwff**2*
     .    (-2d0/xmuw**2*rk*sh+8d0/xmuw*(xmuneut1*(2d0*xmut*xmub+
     .    2d0*(xmut+xmub)-xmut*th-xmub*uh)+2d0*xmut*xmub
     .    -xmut*uh-xmub*th)+4d0*(xmuneut1*(uh+th-xmut-xmub-2d0)+
     .    (xmut+xmub)*(uh+th-1d0)-2d0*xmut*xmub+th*(-th+1d0)+
     .    uh*(-uh+1d0)))+
     .    (ol(nj,ni)**2-or(nj,ni)**2)*vwff**2*8d0*(
     .    xmuneut1*(xmut-xmub+th-uh)+(xmut+xmub)*(th-uh)-xmut+
     .    xmub+th*(-th+1d0)+uh*(uh-1d0))
     .    )
         endif
      else
         chipmw=0d0
      endif
c -------------------------------------------------------------------- c
c                            H+ exchange
c -------------------------------------------------------------------- c
      xmuch = cmass**2/amchar(ni)**2
      dh    = y3-xmuch
      chipmh=0d0
c
      if((amneut(nj)+mt+mbp).le.amchar(ni)) then
         if ((cmass+amneut(nj)).gt.amchar(ni))then
         chipmh=g2s**2/dh**2*(
     .    ql(nj,ni)*qr(nj,ni)*(
     .    (vch**2-ach**2)*
     .    xmneut(nj)/xmchar(ni)*dsqrt(xmut*xmub)*(-16d0)+
     .    (vch**2+ach**2)*
     .    xmneut(nj)/xmchar(ni)*8d0*(1d0+xmuneut1-th-uh) )
     .    +(ql(nj,ni)**2+qr(nj,ni)**2)*(
     .    (vch**2-ach**2)*
     .    dsqrt(xmut*xmub)*4d0*(xmut+xmub-th-uh)+
     .    (vch**2+ach**2)*
     .    2d0*(xmuneut1*(uh+th-xmut-xmub)+(xmut+xmub)*(uh+th-1d0)
     .    -(th+uh)**2+uh+th)) )
         endif
      else
         chipmh=0d0
      endif
c -------------------------------------------------------------------- c
c                     W+ sbottom/stop interference
c -------------------------------------------------------------------- c
      chipmwsbot = 0d0
c
      if((amneut(nj)+mt+mbp).le.amchar(ni)) then
         do i=1,2
        if ((gmsb(i)+mt).gt.amchar(ni).and.
     .(mw+amneut(nj)).gt.amchar(ni))then
            chipmwsbot=chipmwsbot
     .       +g2s**2/dsb(i)/dw*2d0*vwff*(
     .      aksbot(i,ni)*bbot(i,nj)*or(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*dsqrt(xmut*xmub)*(
     .      1/xmuw*(-4d0)*(1d0+xmuneut1+xmut+xmub-uh-th)+16d0) +
     .      alsbot(i,ni)*abot(i,nj)*ol(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*(2d0/xmuw*((xmuneut1+1d0-uh-th)*
     .      (xmut+xmub)+4d0*xmut*xmub)+4d0*(1d0+xmuneut1-
     .      uh-th)) +
     .      aksbot(i,ni)*bbot(i,nj)*ol(nj,ni)*dsqrt(xmut*xmub)*
     .      (2d0/xmuw*(xmuneut1*(xmut+xmub-th-uh+4d0)+xmut+xmub
     .      -th-uh)+4d0*(xmut+xmub-th-uh)) +
     .      alsbot(i,ni)*abot(i,nj)*or(nj,ni)*
     .      (2d0/xmuw*(xmuneut1*(-2d0*xmut*xmub+xmut*th
     .      -2d0*xmut+xmub*uh-2d0*xmub)+xmut*(-2d0*xmub+uh)
     .      +xmub*th)+4d0*(xmuneut1*(xmut-uh+1d0)+xmut*(xmub-uh)
     .      +xmub*(1d0-uh)+uh**2-uh)) +
     .      abot(i,nj)*aksbot(i,ni)*ol(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*dsqrt(xmut)*(2d0/xmuw*
     .      (xmuneut1*(1d0+xmut-uh)+1d0+xmut*(xmub-uh)+xmub*
     .      (xmub-th-2d0*uh+3d0)+th*(uh-1d0)+uh*(uh-2d0))+
     .      4d0*(1d0+xmub-th)) +
     .      bbot(i,nj)*alsbot(i,ni)*or(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*dsqrt(xmub)*
     .      (2d0/xmuw*(xmuneut1*(-1d0-xmut+uh)-1d0+xmub*(-1d0+uh)
     .      +uh*(2d0-th-uh)+th+xmut*(-xmut-xmub+th+2d0*uh-2d0))+
     .      8d0*(1d0+xmut-uh)) +
     .      aksbot(i,ni)*abot(i,nj)*or(nj,ni)*
     .      dsqrt(xmut)*((-2d0)/xmuw*(xmuneut1-uh+xmub)*
     .      (xmuneut1+xmut+xmub-th-uh+1d0)+8d0*(xmuneut1+xmub
     .      -uh)) +
     .      alsbot(i,ni)*bbot(i,nj)*ol(nj,ni)*
     .      dsqrt(xmub)*(2d0/xmuw*(xmuneut1*(xmuneut1+
     .      3d0*xmut-th-2d0*uh+1d0)+xmut*(xmut+xmub-th-2d0*uh)+
     .      xmub*(1d0-uh)+uh*th+uh**2-uh)+4d0*(xmuneut1+xmut-th)) )
            endif
            if (xmustop(i).gt.1d0.and.
     .(mw+amneut(nj)).gt.amchar(ni))then
            chipmwsbot=chipmwsbot
     .      -g2s**2/dst(i)/dw*2d0*vwff*(
     .      akstor(i,ni)*btopr(i,nj)*ol(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*dsqrt(xmut*xmub)*(
     .      1/xmuw*(-4d0)*(1d0+xmuneut1+xmut+xmub-uh-th)+16d0) +
     .      alstor(i,ni)*atopr(i,nj)*or(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*(2d0/xmuw*((xmut+xmub)*(xmuneut1
     .      +1d0-th-uh)+4d0*xmut*xmub)+4d0*(1d0+xmuneut1-uh
     .      -th)) +
     .      akstor(i,ni)*btopr(i,nj)*or(nj,ni)*dsqrt(xmut*xmub)*
     .      (2d0/xmuw*(xmuneut1*(xmut+xmub-th-uh+4d0)+xmut+xmub
     .      -th-uh)+4d0*(xmut+xmub-th-uh)) +
     .      alstor(i,ni)*atopr(i,nj)*ol(nj,ni)*
     .      (2d0/xmuw*(xmuneut1*(-2d0*xmut*xmub+xmut*th-2d0*xmut
     .      +xmub*uh-2d0*xmub)+xmut*(-2d0*xmub+uh)+xmub*th)
     .      +4d0*(xmuneut1*(xmub-th+1d0)+xmut*(xmub-th+1d0)+
     .      th*(-xmub+th-1d0))) +
     .      atopr(i,nj)*akstor(i,ni)*or(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*dsqrt(xmub)*
     .      (2d0/xmuw*(xmuneut1*(1d0+xmub-th)+1d0+xmut*(xmut+
     .      xmub-2d0*th-uh+3d0)+th*(th-xmub+uh-2d0)-uh)
     .      +4d0*(1d0+xmut-uh)) +
     .      alstor(i,ni)*btopr(i,nj)*ol(nj,ni)*
     .      xmneut(nj)/xmchar(ni)*dsqrt(xmut)*
     .      (2d0/xmuw*(xmuneut1*(-1d0-xmub+th)-1d0+xmut*(-xmub+
     .      th-1d0)+xmub*(2d0*th+uh-2d0-xmub)-th*(th+uh)+2d0*
     .      th+uh)+8d0*(1d0+xmub-th)) +
     .      akstor(i,ni)*atopr(i,nj)*ol(nj,ni)*dsqrt(xmub)*(
     .      2d0/xmuw*(xmuneut1*(-xmuneut1-2d0*xmut-xmub+2d0*th+uh
     .      -1d0)+xmut*(-xmut-xmub+2d0*th+uh-1d0)+th*(xmub-th-
     .      uh+1))+8d0*(xmuneut1+xmut-th)) +
     .      alstor(i,ni)*btopr(i,nj)*or(nj,ni)*dsqrt(xmut)*(
     .      2d0/xmuw*(xmuneut1*(xmuneut1+3d0*xmub-uh-2d0*th+1d0)
     .      +xmut*(xmub-th+1d0)+xmub*(xmub-2d0*th-uh)+uh*th+
     .      th**2-th)+4d0*(xmuneut1+xmub-uh)) )
            endif
         enddo
      else
         chipmwsbot=0d0
      endif
c
c -------------------------------------------------------------------- c
c                        H+ sbottom/stop interference
c -------------------------------------------------------------------- c
      chipmhsbot = 0d0
c
      if((amneut(nj)+mt+mbp).le.amchar(ni)) then
         do i=1,2
       if ((gmsb(i)+mt).gt.amchar(ni).and.
     .(cmass+amneut(nj)).gt.amchar(ni))then
            chipmhsbot=chipmhsbot
     .       -g2s**2/dh/dsb(i)*(
     .       (abot(i,nj)*alsbot(i,ni)*ql(nj,ni)*(vch+ach)
     .       +bbot(i,nj)*aksbot(i,ni)*qr(nj,ni)*(vch-ach))*(-2d0)*
     .       dsqrt(xmut)*(xmuneut1+xmub-uh) +
     .       (abot(i,nj)*alsbot(i,ni)*ql(nj,ni)*(vch-ach)
     .       +bbot(i,nj)*aksbot(i,ni)*qr(nj,ni)*(vch+ach))*2d0*
     .       dsqrt(xmub)*(xmuneut1+xmut-th) +
     .       (abot(i,nj)*alsbot(i,ni)*qr(nj,ni)*(vch+ach)
     .       +bbot(i,nj)*aksbot(i,ni)*ql(nj,ni)*(vch-ach))*2d0*
     .       xmneut(nj)/xmchar(ni)*dsqrt(xmut)*(1d0+xmub-th) +
     .       (abot(i,nj)*alsbot(i,ni)*qr(nj,ni)*(vch-ach)
     .       +bbot(i,nj)*aksbot(i,ni)*ql(nj,ni)*(vch+ach))*(-2d0)*
     .       xmneut(nj)/xmchar(ni)*dsqrt(xmub)*(1d0+xmut-uh) +
     .       (bbot(i,nj)*alsbot(i,ni)*ql(nj,ni)*(vch-ach)
     .       +abot(i,nj)*aksbot(i,ni)*qr(nj,ni)*(vch+ach))*2d0*
     .       xmneut(nj)/xmchar(ni)*(1d0+xmuneut1-uh-th) +
     .       (abot(i,nj)*aksbot(i,ni)*qr(nj,ni)*(vch-ach)
     .       +bbot(i,nj)*alsbot(i,ni)*ql(nj,ni)*(vch+ach))*(-4d0)*
     .       xmneut(nj)/xmchar(ni)*dsqrt(xmut*xmub) +
     .       (abot(i,nj)*aksbot(i,ni)*ql(nj,ni)*(vch-ach)
     .       +bbot(i,nj)*alsbot(i,ni)*qr(nj,ni)*(vch+ach))*2d0*
     .       dsqrt(xmut*xmub)*(xmut+xmub-th-uh) +
     .       (alsbot(i,ni)*bbot(i,nj)*qr(nj,ni)*(vch-ach)
     .       +aksbot(i,ni)*abot(i,nj)*ql(nj,ni)*(vch+ach))*2d0*
     .       (-uh**2-uh*th+uh*(1d0+xmut+xmub+xmuneut1)-xmub-
     .       xmuneut1*xmut) )
            endif
             if (xmustop(i).gt.1d0.and.
     .(cmass+amneut(nj)).gt.amchar(ni))then
            chipmhsbot=chipmhsbot
     .       +2d0*g2s**2/dh/dst(i)*(
     .       (btopr(i,nj)*akstor(i,ni)*ql(nj,ni)*(vch-ach)
     .       +atopr(i,nj)*alstor(i,ni)*qr(nj,ni)*(vch+ach))*
     .       dsqrt(xmut)*(uh-xmuneut1-xmub) +
     .       (atopr(i,nj)*alstor(i,ni)*qr(nj,ni)*(vch-ach)
     .       +btopr(i,nj)*akstor(i,ni)*ql(nj,ni)*(vch+ach))*
     .       dsqrt(xmub)*(-th+xmut+xmuneut1) +
     .       (atopr(i,nj)*alstor(i,ni)*ql(nj,ni)*(vch+ach)
     .       +btopr(i,nj)*akstor(i,ni)*qr(nj,ni)*(vch-ach))*
     .       dsqrt(xmut)*xmneut(nj)/xmchar(ni)*(-th+1d0+xmub) +
     .       (atopr(i,nj)*alstor(i,ni)*ql(nj,ni)*(vch-ach)
     .       +btopr(i,nj)*akstor(i,ni)*qr(nj,ni)*(vch+ach))*
     .       dsqrt(xmub)*xmneut(nj)/xmchar(ni)*(uh-1d0-xmut) +
     .       (btopr(i,nj)*alstor(i,ni)*qr(nj,ni)*(vch-ach)
     .       +atopr(i,nj)*akstor(i,ni)*ql(nj,ni)*(vch+ach))*
     .       2d0*xmneut(nj)/xmchar(ni)*dsqrt(xmut*xmub) +
     .       (atopr(i,nj)*akstor(i,ni)*qr(nj,ni)*(vch-ach)
     .       +btopr(i,nj)*alstor(i,ni)*ql(nj,ni)*(vch+ach))*
     .       (uh*th+th**2-th*(1d0+xmut+xmub+xmuneut1)+xmut+
     .       xmub*xmuneut1) +
     .       (atopr(i,nj)*akstor(i,ni)*ql(nj,ni)*(vch-ach)
     .       +btopr(i,nj)*alstor(i,ni)*qr(nj,ni)*(vch+ach))*
     .       xmneut(nj)/xmchar(ni)*(uh+th-xmuneut1-1d0) +
     .       (btopr(i,nj)*alstor(i,ni)*ql(nj,ni)*(vch-ach)
     .       +atopr(i,nj)*akstor(i,ni)*qr(nj,ni)*(vch+ach))*
     .       dsqrt(xmut*xmub)*(uh+th-xmut-xmub) )
         endif
         enddo
      else
         chipmhsbot=0d0
      endif
c
c -------------------------------------------------------------------- c
c 	               interference W+ H-
c -------------------------------------------------------------------- c
      chipmwh=0d0
c
      if((amneut(nj)+mt+mbp).le.amchar(ni)) then
       if ((cmass+amneut(nj)).gt.amchar(ni).and.
     .(mw+amneut(nj)).gt.amchar(ni))then
         chipmwh=chipmwh-2d0*g2s**2/dh/dw*vwff*(
     .    (ol(nj,ni)*ql(nj,ni)+or(nj,ni)*qr(nj,ni))*(
     .    (vch-ach)*xmneut(nj)/xmchar(ni)*dsqrt(xmub)*(
     .    2d0/xmuw*(-2d0-xmuneut1*(xmut+xmub-th-uh+2d0)-xmut*
     .    (2d0*xmut+xmub-3d0*th-3d0*uh+5d0)-xmub*(xmut-th-uh
     .    +1d0)-(th+uh)**2+3d0*(th+uh))+4d0*(-uh+xmut+1d0)) +
     .    (vch+ach)*xmneut(nj)/xmchar(ni)*dsqrt(xmut)*(
     .    2d0/xmuw*(2d0+xmuneut1*(xmut+xmub-th-uh+2d0)+xmut*
     .    (-th-uh+1d0+xmub)+xmub*(xmut+2d0*xmub-3d0*(th+uh)+5d0)
     .    +(th+uh)**2-3d0*(th+uh))+4d0*(th-xmub-1d0) ) )
     .    +(ol(nj,ni)*qr(nj,ni)+or(nj,ni)*ql(nj,ni))*(
     .    (vch-ach)*dsqrt(xmub)*(
     .    2d0/xmuw*(2d0*xmuneut1**2-xmuneut1*(-5d0*xmut-xmub+
     .    3d0*(th+uh)-2d0)-xmut*(-2d0*xmut-xmub+3d0*(th+uh)-1d0)-
     .    xmub*(-xmut+uh+th-1d0)+(th+uh)**2-(th+uh))
     .    -4d0*(xmuneut1+xmut-th) )+
     .    (vch+ach)*dsqrt(xmut)*(
     .    2d0/xmuw*(-2d0*xmuneut1**2+xmuneut1*(-xmut-5d0*xmub+3d0*
     .    (th+uh)-2d0)+xmut*(-xmub+th+uh-1d0)+xmub*(-xmut-2d0*xmub+
     .    3d0*(th+uh)-1d0)-(th+uh)**2+(th+uh)) +4d0*(xmuneut1+xmub
     .    -uh) ) ) )
       endif
      else
         chipmwh=0d0
      endif
c
      NS_chipmtop = chipmsfer+chipmw+chipmh+chipmwsbot+chipmhsbot+
     .              chipmwh
c
      end
c
c ==================================================================== c
c =========================== chargino1 e+ e- ======================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_charel(x1,x2)
c
      IMPLICIT NONE
      INTEGER ni,nj,i,k
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsel(2)
      DOUBLE PRECISION opl(2,2),opr(2,2),onl(5,5),onr(5,5)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                   azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION xmuneut1,x1,x2,x3,y1,y2,y3,xmusnel(2),
     .charsnel,xmuz,dz,charz,charzsnel
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
c
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup4/opl,opr,onl,onr
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_pi/PI,SQR2
c
      xmuneut1 = amchar(1)**2/amchar(2)**2
c
      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                          sneutrino exchange
c -------------------------------------------------------------------- c
      xmusnel(1) = asne1**2/amchar(2)**2
      xmusnel(2) = asne2**2/amchar(2)**2
c
      dsel(1)  = 1d0-x2-xmusnel(1)
      dsel(2)  = 1d0-x2-xmusnel(2)
c
      charsnel=0d0
c
      if(amchar(1).le.amchar(2)) then
         do i=1,2
            do k=1,2
               if (xmusnel(i).gt.1d0.and.xmusnel(k).gt.1d0)then
               charsnel=charsnel+g2s**2/dsel(k)/dsel(i)*x2*y2*
     .         (alsne(k,2)*alsne(i,2)+blsne(k,2)*blsne(i,2))*
     .         (alsne(k,1)*alsne(i,1)+blsne(k,1)*blsne(i,1))
               endif
            enddo
         enddo
      else
         charsnel=0d0
      endif
c
c -------------------------------------------------------------------- c
c                             Z exchange
c -------------------------------------------------------------------- c
c Note: azztautau,vzztautau = azzelel,vzzelel
      xmuz  = mz**2/amchar(2)**2
      dz    = y3-xmuz
c
         charz=0d0
c
      if(amchar(1).le.amchar(2)) then
        if((mz+amchar(1)).gt.amchar(2))then
         charz=g2s**2*4d0/dz**2/cw**2*
     .        (((azztautau+vzztautau)**2*opl(1,2)**2
     .        + (azztautau-vzztautau)**2*opr(1,2)**2)*x2*y2
     .        +((azztautau+vzztautau)**2*opr(1,2)**2
     .        + (azztautau-vzztautau)**2*opl(1,2)**2)*x1*y1
     .        -4d0*xmchar(1)/xmchar(2)*opl(1,2)*opr(1,2)
     .        *(azztautau**2+vzztautau**2)*y3 )
        endif
      else
         charz=0d0
      endif
c -------------------------------------------------------------------- c
c                      Z-sneutrino interference
c -------------------------------------------------------------------- c
      charzsnel=0d0
c
      if(amchar(1).le.amchar(2)) then
         do i=1,2
            if (xmusnel(i).gt.1d0
     ..and.(mz+amchar(1)).gt.amchar(2))then
            charzsnel=charzsnel+g2s**2*4d0/dsel(i)/dz/cw
     .           *((alsne(i,2)*alsne(i,1)*opl(1,2)*
     .           (azztautau+vzztautau)+
     .           blsne(i,2)*blsne(i,1)*opr(1,2)*
     .           (-azztautau+vzztautau))*x2*y2
     .           -(alsne(i,1)*alsne(i,2)*opr(1,2)*
     .           (azztautau+vzztautau)
     .           +blsne(i,1)*blsne(i,2)*opl(1,2)*
     .           (-azztautau+vzztautau)
     .           )*xmchar(1)/xmchar(2)*y3)
            endif
         enddo
      else
         charzsnel=0d0
      endif
c
c -------------------------------------------------------------------- c
      NS_charel = charsnel+charz+charzsnel

      end
c
c ==================================================================== c
c =========================== chargino1 mu+ mu- ====================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_charmu(x1,x2)
c
      IMPLICIT NONE
      INTEGER ni,nj,i,k
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsel(2)
      DOUBLE PRECISION opl(2,2),opr(2,2),onl(5,5),onr(5,5)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                   azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION xmuneut1,x1,x2,x3,y1,y2,y3,xmusnel(2),
     .charsnel,xmuz,dz,charz,charzsnel
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
c
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup4/opl,opr,onl,onr
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_pi/PI,SQR2
c
      xmuneut1 = amchar(1)**2/amchar(2)**2
c
      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                          sneutrino exchange
c -------------------------------------------------------------------- c
      xmusnel(1) = asnmu1**2/amchar(2)**2
      xmusnel(2) = asnmu2**2/amchar(2)**2
c
      dsel(1)  = 1d0-x2-xmusnel(1)
      dsel(2)  = 1d0-x2-xmusnel(2)
c
      charsnel=0d0
c
      if(amchar(1).le.amchar(2)) then
         do i=1,2
            do k=1,2
               if (xmusnel(i).gt.1d0.and.xmusnel(k).gt.1d0)then
               charsnel=charsnel+g2s**2/dsel(k)/dsel(i)*x2*y2*
     .         (alsnm(k,2)*alsnm(i,2)+blsnm(k,2)*blsnm(i,2))*
     .         (alsnm(k,1)*alsnm(i,1)+blsnm(k,1)*blsnm(i,1))
               endif
            enddo
         enddo
      else
         charsnel=0d0
      endif
c
c -------------------------------------------------------------------- c
c                             Z exchange
c -------------------------------------------------------------------- c
c Note: azztautau,vzztautau = azzelel,vzzelel
      xmuz  = mz**2/amchar(2)**2
      dz    = y3-xmuz
c
         charz=0d0
c
      if(amchar(1).le.amchar(2)) then
        if((mz+amchar(1)).gt.amchar(2))then
         charz=g2s**2*4d0/dz**2/cw**2*
     .        (((azztautau+vzztautau)**2*opl(1,2)**2
     .        + (azztautau-vzztautau)**2*opr(1,2)**2)*x2*y2
     .        +((azztautau+vzztautau)**2*opr(1,2)**2
     .        + (azztautau-vzztautau)**2*opl(1,2)**2)*x1*y1
     .        -4d0*xmchar(1)/xmchar(2)*opl(1,2)*opr(1,2)
     .        *(azztautau**2+vzztautau**2)*y3 )
        endif
      else
         charz=0d0
      endif
c -------------------------------------------------------------------- c
c                      Z-sneutrino_el interference
c -------------------------------------------------------------------- c
      charzsnel=0d0
c
      if(amchar(1).le.amchar(2)) then
         do i=1,2
            if (xmusnel(i).gt.1d0.and.(mz+amchar(1)).gt.amchar(2))then
            charzsnel=charzsnel+g2s**2*4d0/dsel(i)/dz/cw
     .           *((alsnm(i,2)*alsnm(i,1)*opl(1,2)*
     .           (azztautau+vzztautau)+
     .           blsnm(i,2)*blsnm(i,1)*opr(1,2)*
     .           (-azztautau+vzztautau))*x2*y2
     .           -(alsnm(i,1)*alsnm(i,2)*opr(1,2)*
     .           (azztautau+vzztautau)
     .           +blsnm(i,1)*blsnm(i,2)*opl(1,2)*
     .           (-azztautau+vzztautau)
     .           )*xmchar(1)/xmchar(2)*y3)
            endif
         enddo
      else
         charzsnel=0d0
      endif
c
c -------------------------------------------------------------------- c
      NS_charmu = charsnel+charz+charzsnel

      end
c ==================================================================== c
c ========================== chargino1 tau+ tau- ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_chartau(x1,x2)
c
      IMPLICIT NONE
      INTEGER i,j,k
      DOUBLE PRECISION dsto(2),sgn(2)
      DOUBLE PRECISION opl(2,2),opr(2,2),onl(5,5),onr(5,5)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION P(2,3)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION hchichi(3,5,5),achichi(2,5,5)
      DOUBLE PRECISION hchachaR(3,2,2),hchachaL(3,2,2),achachaR(2,2,2),
     .     achachaL(2,2,2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION vzz,azz,xmuneut1,xmusn(2),xmutau,
     .x1,x2,x3,y3,uh,th,charsntau,xmuz,dz,charztau,rh,sh,rk,
     .dhl(3),dhh(3),dha(2),dhna(2),charhtau(3,3),
     .charzsntau,hv1,hv2,hv3,hv4,hv5,hv6,hv7,hv8,charhsntau(3),
     .charhasntau(2),charzh(3),charza(2),charhatau(2,2)
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
c
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_coup4/opl,opr,onl,onr
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_CPodd_MIX/P
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_coupNMSSM/hchichi,achichi,hchachaR,hchachaL,achachaR,
     .     achachaL
c
      do i=1,2
         sgn(i) = 1d0
         if(xmchar(i).ge.0d0) then
            sgn(i) = 1d0
         elseif(xmchar(i).lt.0d0) then
            sgn(i) = -1d0
         endif
      enddo
      vzz = vzztautau
      azz = azztautau
c
      xmuneut1 = amchar(1)**2/amchar(2)**2
      xmusn(1)   = asntau1**2/amchar(2)**2
      xmusn(2)   = asntau2**2/amchar(2)**2
      xmutau   = mtau**2/amchar(2)**2
c
      x3 = 2d0-x1-x2
      y3 = 1d0+xmuneut1-x3
c
      uh = 1d0-x1+xmutau
      th = 1d0-x2+xmutau
c -------------------------------------------------------------------- c
c                        sneutrino_tau exchange
c -------------------------------------------------------------------- c
      dsto(1) = 1d0-x2-xmusn(1)+xmutau
      dsto(2) = 1d0-x2-xmusn(2)+xmutau
c
      charsntau=0d0
c
      if((amchar(1)+2d0*mtau).le.amchar(2)) then
         do i=1,2
            do k=1,2
      if (xmusn(k).gt.1d0.and.xmusn(i).gt.1d0)then
               charsntau=charsntau
     .           +g2s**2/dsto(k)/dsto(i)*(
     .           (alsnt(i,2)*blsnt(k,2)+blsnt(i,2)*alsnt(k,2))*
     .           (alsnt(i,1)*blsnt(k,1)+blsnt(i,1)*alsnt(k,1))*
     .           xmchar(1)/xmchar(2)*xmutau*(-4d0)+
     .           (alsnt(i,2)*alsnt(k,2)+blsnt(i,2)*blsnt(k,2))*
     .           (alsnt(i,1)*blsnt(k,1)+blsnt(i,1)*alsnt(k,1))*
     .           xmchar(1)/xmchar(2)*dsqrt(xmutau)*2d0*
     .           (th-xmutau-1d0)+
     .           (alsnt(i,2)*blsnt(k,2)+blsnt(i,2)*alsnt(k,2))*
     .           (alsnt(i,1)*alsnt(k,1)+blsnt(i,1)*blsnt(k,1))*
     .           dsqrt(xmutau)*2d0*(th-xmutau-xmuneut1)+
     .           (alsnt(i,2)*alsnt(k,2)+blsnt(i,2)*blsnt(k,2))*
     .           (alsnt(i,1)*alsnt(k,1)+blsnt(i,1)*blsnt(k,1))*
     .           (-th**2+th*(1d0+xmuneut1+2d0*xmutau)-
     .           (xmuneut1+xmutau)*(1d0+xmutau)))
       endif
            enddo
         enddo
      else
         charsntau=0d0
      endif
c
c -------------------------------------------------------------------- c
c                            Z exchange
c -------------------------------------------------------------------- c
      xmuz = mz**2/amchar(2)**2
      dz   = y3-xmuz
c
      charztau = 0d0
c
      rh = xmuneut1+2d0*xmutau-th-uh+1d0
      sh = (xmuneut1-th-uh+1d0)*2d0*xmutau+4d0*xmutau**2
      rk = xmuneut1*(2d0*xmutau-th-uh+4d0)+2d0*xmutau-uh-th
c
      if((amchar(1)+2d0*mtau).le.amchar(2)) then
         if ((mz+amchar(1)).gt.amchar(2))then
         charztau=g2s**2/dz**2/cw**2*(
     .    opl(1,2)*opr(1,2)*(vzz**2-azz**2)*
     .    xmchar(1)/xmchar(2)*xmutau*(-16d0/xmuz**2*rh**2+
     .    32d0/xmuz*rh-64d0)+
     .    opl(1,2)*opr(1,2)*(vzz**2+azz**2)*
     .    xmchar(1)/xmchar(2)*(8d0/xmuz**2*rh*sh-16d0/xmuz*sh
     .    -16d0*(xmuneut1-uh-th+1d0))+
     .    (opl(1,2)**2+opr(1,2)**2)*(vzz**2-azz**2)*
     .    xmutau*(4d0/xmuz**2*rh*rk-8d0/xmuz*rk+8d0*(uh+th-
     .    2d0*xmutau))
     .    +(opl(1,2)**2+opr(1,2)**2)*(vzz**2+azz**2)*
     .    (-2d0/xmuz**2*rk*sh+8d0/xmuz*(xmuneut1*(2d0*xmutau**2+
     .    4d0*xmutau-xmutau*(th+uh))+2d0*xmutau**2-xmutau*(uh+th))+
     .    4d0*(xmuneut1*(uh+th-2d0*xmutau-2d0)+2d0*xmutau*(uh+th-
     .    1d0)-2d0*xmutau**2+th*(-th+1d0)+uh*(-uh+1d0)))+
     .    (opl(1,2)**2-opr(1,2)**2)*vzz*azz*8d0*(
     .    xmuneut1*(th-uh)+2d0*xmutau*(th-uh)+th*(-th+1d0)+uh*(uh-
     .    1d0)))
         endif
      else
         charztau=0d0
      endif
c --------------------------------------------------------------------c
c       NMSSM CP EVEN Higgs exchange + Interference H(i)-H(j)
c --------------------------------------------------------------------c
      do i=1,3
         do j=1,3
c
      dhl(i)   = y3-smass(i)**2/amchar(2)**2
      dhh(j)   = y3-smass(j)**2/amchar(2)**2
c
         charhtau(i,j) =0d0
      if ((smass(i)+amchar(1)).gt.amchar(2).and.
     .(smass(j)+amchar(1)).gt.amchar(2))then
      if((amchar(1)+2d0*mtau).le.amchar(2)) then
             charhtau(i,j) = g2s**2/dhh(j)/dhl(i)
     .    *scaltau*S(i,2)*scaltau*S(j,2)*(
     .    (hchachaL(i,1,2)*hchachaR(j,1,2)+
     .    hchachaR(i,1,2)*hchachaL(j,1,2))*(
     .    xmchar(1)/xmchar(2)*xmutau*(-4d0)+
     .    xmchar(1)/xmchar(2)*2d0*(1d0+xmuneut1-uh-th))+
     .    (hchachaL(i,1,2)*hchachaL(j,1,2)+
     .    hchachaR(i,1,2)*hchachaR(j,1,2))*(
     .    xmutau*2d0*(2d0*xmutau-th-uh)+
     .    (xmuneut1*(uh+th-2d0*xmutau)+2d0*xmutau*(uh+th-1d0)
     .    +th+uh-(th+uh)**2)))
       else
         charhtau(i,j) =0d0
       endif
       endif
c
       enddo
      enddo
c -------------------------------------------------------------------- c
c     NMSSM CP ODD Higgs exchange+Interference A(i)-A(j)
c -------------------------------------------------------------------- c
      do i=1,2
         do j=1,2
      dha(i)   = y3-pmass(i)**2/amchar(2)**2
      dhna(j)  = y3-pmass(j)**2/amchar(2)**2
      charhatau(i,j)=0d0
      if ((pmass(i)+amchar(1)).gt.amchar(2).and.
     .(pmass(j)+amchar(1)).gt.amchar(2))then
      if((amchar(1)+2d0*mtau).le.amchar(2)) then

         charhatau(i,j)=g2s**2/dha(i)/dhna(j)*(scaltau*P(i,2))*
     .    (scaltau*P(j,2))*(
     .    achachaL(i,1,2)*achachaR(j,1,2)
     .    *xmchar(1)/xmchar(2)*xmutau*8d0+
     .    achachaL(i,1,2)*achachaR(j,1,2)*xmchar(1)/xmchar(2)*4d0*
     .    (1d0+xmuneut1-th-uh)+
     .    (achachaL(i,1,2)*achachaL(j,1,2)
     .    +achachaR(i,1,2)*achachaR(j,1,2))*xmutau*2d0*
     .    (-2d0*xmutau+th+uh)+
     .    (achachaL(i,1,2)*achachaL(j,1,2)+achachaR(i,1,2)*
     .    achachaR(j,1,2))*(xmuneut1*
     .    (uh+th-2d0*xmutau)+2d0*xmutau*(uh+th-1d0)
     .     -(th+uh)**2+uh+th))
      else
         charhatau(i,j)=0d0
      endif
      endif
c
        enddo
      enddo
c -------------------------------------------------------------------- c
c                    Z sneutrino_tau interference
c -------------------------------------------------------------------- c
      charzsntau=0d0
c
      if((amchar(1)+2d0*mtau).le.amchar(2)) then
         do i=1,2
            hv1 = alsnt(i,2)*alsnt(i,1)*opr(2,1)*(vzz-azz)
     .           +blsnt(i,2)*blsnt(i,1)*opl(2,1)*(vzz+azz)
            hv2 = alsnt(i,2)*alsnt(i,1)*opr(2,1)*(vzz+azz)
     .           +blsnt(i,2)*blsnt(i,1)*opl(2,1)*(vzz-azz)
            hv3 = alsnt(i,2)*alsnt(i,1)*opl(2,1)*(vzz-azz)
     .           +blsnt(i,2)*blsnt(i,1)*opr(2,1)*(vzz+azz)
            hv4 = alsnt(i,2)*alsnt(i,1)*opl(2,1)*(vzz+azz)
     .           +blsnt(i,2)*blsnt(i,1)*opr(2,1)*(vzz-azz)
            hv5 = alsnt(i,1)*blsnt(i,2)*opr(2,1)*(vzz+azz)
     .           +blsnt(i,1)*alsnt(i,2)*opl(2,1)*(vzz-azz)
            hv6 = alsnt(i,1)*blsnt(i,2)*opr(2,1)*(vzz-azz)
     .           +blsnt(i,1)*alsnt(i,2)*opl(2,1)*(vzz+azz)
            hv7 = alsnt(i,1)*blsnt(i,2)*opl(2,1)*(vzz+azz)
     .           +blsnt(i,1)*alsnt(i,2)*opr(2,1)*(vzz-azz)
            hv8 = alsnt(i,1)*blsnt(i,2)*opl(2,1)*(vzz-azz)
     .           +blsnt(i,1)*alsnt(i,2)*opr(2,1)*(vzz+azz)
            if (xmusn(i).gt.1d0.and.(mz+amchar(1)).gt.amchar(2))then
            charzsntau=charzsntau
     .      -g2s**2/dsto(i)/dz/cw*(
     .      hv1*xmchar(1)/xmchar(2)*xmutau*(
     .      1/xmuz*(-4d0)*(1d0+xmuneut1+2d0*xmutau-uh-th)+16d0) +
     .      hv2*xmchar(1)/xmchar(2)*(2d0/xmuz*(2d0*xmutau*(xmuneut1+
     .      1d0-th-uh)+4d0*xmutau**2)+4d0*(1d0+xmuneut1-uh-th)) +
     .      hv3*xmutau*(2d0/xmuz*(xmuneut1*(2d0*xmutau-th-uh+4d0)+
     .      2d0*xmutau-th-uh)+4d0*(2d0*xmutau-th-uh)) +
     .      hv4*(2d0/xmuz*(xmuneut1*(-2d0*xmutau**2+xmutau*th
     .      -2d0*xmutau+xmutau*uh-2d0*xmutau)+xmutau*(-2d0*xmutau
     .      +uh)+xmutau*th)+4d0*(xmuneut1*(xmutau-th+1d0)+xmutau*
     .      (xmutau-th)+xmutau*(1d0-th)+th**2-th)) +
     .      hv5*xmchar(1)/xmchar(2)*dsqrt(xmutau)*
     .      (2d0/xmuz*(xmuneut1*(1d0+xmutau-th)+1d0+xmutau*
     .      (2d0*xmutau-2d0*th-uh+3d0)+th*(th-xmutau+uh-2d0)-uh)+
     .      4d0*(1d0+xmutau-uh)) +
     .      hv6*xmchar(1)/xmchar(2)*dsqrt(xmutau)*
     .      (2d0/xmuz*(xmuneut1*(-1d0-xmutau+th)-1d0+xmutau*
     .      (-xmutau+th-1d0)+xmutau*(2d0*th+uh-2d0-xmutau)
     .      -th*(th+uh)+2d0*th+uh)+8d0*(1d0+xmutau-th)) +
     .      hv7*dsqrt(xmutau)*(2d0/xmuz*(xmuneut1*(-xmuneut1
     .      -3d0*xmutau+2d0*th+uh-1d0)+xmutau*(-2d0*xmutau+2d0*th+
     .      uh-1d0)+th*(xmutau-th-uh+1))+8d0*(xmuneut1+xmutau-th))+
     .      hv8*dsqrt(xmutau)*(2d0/xmuz*(xmuneut1*(xmuneut1+
     .      3d0*xmutau-uh-2d0*th+1d0)+xmutau*(xmutau-th+1d0)+
     .      xmutau*(xmutau-2d0*th-uh)+uh*th+th**2-th)+4d0*(xmuneut1
     .      +xmutau-uh)) )
            endif
         enddo
      else
         charzsntau=0d0
      endif
c
c -------------------------------------------------------------------- c
c                     Higgs-sneutrino_tau interference
c -------------------------------------------------------------------- c
      do j=1,3
         dhl(j)   = y3-smass(j)**2/amchar(2)**2
c
         charhsntau(j)=0d0
c
       if ((amchar(1)+2d0*mtau).le.amchar(2)) then
         do i=1,2
       if ((smass(j)+amchar(1)).gt.amchar(2).and.xmusn(i).gt.1d0)then
            charhsntau(j)=charhsntau(j)
     .    +2d0*g2s**2/dhl(j)/dsto(i)*(scaltau*S(j,2)/dsqrt(2d0))*(
     .       (alsnt(i,1)*alsnt(i,2)*hchachaR(j,1,2)+
     .        blsnt(i,2)*blsnt(i,1)*hchachaL(j,1,2))*(
     .       dsqrt(xmutau)*(uh-xmuneut1-xmutau) +
     .       dsqrt(xmutau)*(-th+xmutau+xmuneut1) ) +
     .       (alsnt(i,1)*alsnt(i,2)*hchachaL(j,1,2)+
     .        blsnt(i,2)*blsnt(i,1)*hchachaR(j,1,2))*(
     .       dsqrt(xmutau)*xmchar(1)/xmchar(2)*(-th+1d0+xmutau) +
     .       dsqrt(xmutau)*xmchar(1)/xmchar(2)*(uh-1d0-xmutau) ) +
     .       (alsnt(i,1)*blsnt(i,2)*hchachaL(j,1,2)
     .       +alsnt(i,2)*blsnt(i,1)*hchachaR(j,1,2))*(
     .       2d0*xmutau*xmchar(1)/xmchar(2) +
     .       xmchar(1)/xmchar(2)*(uh+th-xmuneut1-1d0) ) +
     .       (alsnt(i,1)*blsnt(i,2)*hchachaR(j,1,2)
     .       +alsnt(i,2)*blsnt(i,1)*hchachaL(j,1,2))*(
     .       (uh*th+th**2-th*(1d0+2d0*xmutau+xmuneut1)+xmutau+
     .       xmutau*xmuneut1) +
     .       xmutau*(uh+th-2d0*xmutau)) )
       endif
        enddo
      else
         charhsntau(j)=0d0
      endif
      enddo
c -------------------------------------------------------------------- c
c                     A-sneutrino_tau interference
c -------------------------------------------------------------------- c
      do j=1,2
c
      dha(j) = y3-pmass(j)**2/amchar(2)**2
      charhasntau(j)=0d0
c
      if ((amchar(1)+2d0*mtau).le.amchar(2)) then
         do i=1,2
      if ((pmass(j)+amchar(1)).gt.amchar(2).and.xmusn(i).gt.1d0)then
            charhasntau(j)=charhasntau(j)
     .       +2d0*g2s**2/dha(j)/dsto(i)
     .           *(scaltau*(-P(j,2))/dsqrt(2d0))*(
     .       (alsnt(i,1)*alsnt(i,2)*achachaR(j,1,2)-
     .        blsnt(i,2)*blsnt(i,1)*achachaL(j,1,2))*(
     .       dsqrt(xmutau)*sgn(2)*(uh-xmuneut1-xmutau) +
     .       dsqrt(xmutau)*sgn(2)*(-th+xmutau+xmuneut1)*(-1d0) ) +
     .       (alsnt(i,1)*alsnt(i,2)*achachaL(j,1,2)-
     .        blsnt(i,2)*blsnt(i,1)*achachaR(j,1,2))*(
     .       dsqrt(xmutau)*sgn(2)*xmchar(1)/xmchar(2)*
     .       (-th+1d0+xmutau) +
     .       dsqrt(xmutau)*sgn(2)*xmchar(1)/xmchar(2)*
     .       (uh-1d0-xmutau)*(-1d0)
     .       ) +
     .       (alsnt(i,1)*blsnt(i,2)*achachaL(j,1,2)
     .       -alsnt(i,2)*blsnt(i,1)*achachaR(j,1,2))*(
     .       2d0*xmutau*xmchar(1)/xmchar(2) +
     .       xmchar(1)/xmchar(2)*(uh+th-xmuneut1-1d0)*(-1d0) ) +
     .       (alsnt(i,1)*blsnt(i,2)*achachaR(j,1,2)
     .       -alsnt(i,2)*blsnt(i,1)*achachaL(j,1,2))*(
     .       (uh*th+th**2-th*(1d0+2d0*xmutau+xmuneut1)+xmutau+
     .       xmutau*xmuneut1)*(-1d0) +
     .       xmutau*(uh+th-2d0*xmutau)) )
       endif
         enddo
      else
         charhasntau(j)=0d0
      endif
      enddo
c -------------------------------------------------------------------- c
c 	                interference Z and H/h/NH/A/NA
c -------------------------------------------------------------------- c
      do i=1,3
c
      dhl(i) = y3-smass(i)**2/amchar(2)**2
      charzh(i)=0d0
      if ((smass(i)+amchar(1)).gt.amchar(2).and.
     .(mz+amchar(1)).gt.amchar(2))then
      if ((amchar(1)+2d0*mtau).le.amchar(2)) then
         charzh(i)=
     .  -2d0*g2s**2/cw/dhl(i)/dz*vzz*(scaltau*S(i,2))/dsqrt(2d0)*(
     .    xmchar(1)/xmchar(2)*dsqrt(xmutau)*
     .    (hchachaL(i,1,2)*opl(1,2)+hchachaR(i,1,2)*opr(1,2))*(
     .    4d0*(th-xmutau-1d0) - 4d0*(uh-xmutau-1d0) )
     .    +dsqrt(xmutau)*(hchachaL(i,1,2)*opr(1,2)
     .      +hchachaR(i,1,2)*opl(1,2))*(
     .    4d0*(xmuneut1+xmutau-uh)-4d0*(xmuneut1+xmutau-th)) )
      else
         charzh(i)=0d0
      endif
      endif
c
      enddo
      do i=1,2
c
       dha(i)   = y3-pmass(i)**2/amchar(2)**2
       charza(i)=0d0
      if ((pmass(i)+amchar(1)).gt.amchar(2).and.
     .(mz+amchar(1)).gt.amchar(2))then
      if ((amchar(1)+2d0*mtau).le.amchar(2)) then
         charza(i)=
     .    -2d0*g2s**2/cw/dha(i)/dz*azz*(-scaltau*P(i,2))/dsqrt(2d0)*(
     .    xmchar(1)/xmchar(2)*dsqrt(xmutau)*
     .    (achachaL(i,1,2)*opl(1,2)+achachaR(i,1,2)*opr(1,2))*(
     .    4d0/xmuz*(2d0+xmuneut1*(2d0*xmutau-th-uh+2d0)+xmutau*
     .    (-uh-th+1d0+xmutau)+xmutau*(3d0*xmutau-3d0*(th+uh)+5d0)
     .    +(th+uh)**2-3d0*(th+uh))+
     .    4d0*(th-xmutau-1d0) +4d0*(uh-xmutau-1d0) )
     .    +dsqrt(xmutau)
     .     *(achachaL(i,1,2)*opr(1,2)+achachaR(i,1,2)*opl(1,2))*(
     .    4d0/xmuz*(-2d0*xmuneut1**2+xmuneut1*(-6d0*xmutau+
     .    3d0*(th+uh)-2d0)+xmutau*(-xmutau+th+uh-1d0)+xmutau*
     .    (-3d0*xmutau+3d0*(uh+th)-1d0)-(th+uh)**2+(th+uh))+
     .    4d0*(xmuneut1+xmutau-uh)+4d0*(xmuneut1+xmutau-th)) )
      else
       charza(i)=0d0
      endif
      endif
c
      enddo
      NS_chartau = charsntau+charztau+charzsntau
     .             +charhasntau(1)+charhasntau(2)
     .             +charza(1)+charza(2)
        do i=1,3
          NS_chartau = NS_chartau+charhsntau(i)
          NS_chartau = NS_chartau+charzh(i)
        do j=1,3
          NS_chartau = NS_chartau+charhtau(i,j)
        enddo
        enddo
        do i=1,2
           do j=1,2
        NS_chartau = NS_chartau+charhatau(i,j)
        enddo
        enddo
        end
c
c ==================================================================== c
c ======================= chargino1 nu_e nubar_e ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_charnue(x1,x2)
c
      IMPLICIT NONE
      INTEGER I,k
      DOUBLE PRECISION dsel(2),ble(2,2)
      DOUBLE PRECISION opl(2,2),opr(2,2),onl(5,5),onr(5,5)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION xmuneut1,X1,X2,X3,charsel,charz,charzsel,y1,
     .y2,y3,xmusel(2),xmuz,dz
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
c
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_coup4/opl,opr,onl,onr
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
c
      xmuneut1 = amchar(1)**2/amchar(2)**2
c
      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c
      do i=1,2
         ble(1,i) = 0d0
         ble(2,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c                          selectron exchange
c -------------------------------------------------------------------- c
      xmusel(1) = ase1**2/amchar(2)**2
      xmusel(2) = ase2**2/amchar(2)**2
c
      dsel(1)  = 1d0-x1-xmusel(1)
      dsel(2)  = 1d0-x1-xmusel(2)
c
      charsel=0d0
c
      if (amchar(1).le.amchar(2)) then
         do i=1,2
            do k=1,2
      if (xmusel(k).gt.1.0d0.and.xmusel(i).gt.1.0d0)then
               charsel=charsel+g2s**2/dsel(k)/dsel(i)*x1*y1*
     .         (ale(k,2)*ale(i,2)+ble(k,2)*ble(i,2))*
     .         (ale(k,1)*ale(i,1)+ble(k,1)*ble(i,1))
      endif
            enddo
         enddo
      else
         charsel=0d0
      endif
c -------------------------------------------------------------------- c
c                             Z exchange
c -------------------------------------------------------------------- c
      xmuz  = mz**2/amchar(2)**2
      dz    = y3-xmuz
c
      charz=0d0
c
      if (amchar(1).le.amchar(2)) then
         if ((mz+amchar(1)).gt.amchar(2))then
         charz=charz+g2s**2*4d0/dz**2/cw**2*
     .        (((azzneutneut+vzzneutneut)**2*opl(1,2)**2
     .        + (azzneutneut-vzzneutneut)**2*opr(1,2)**2)*x2*y2
     .        +((azzneutneut+vzzneutneut)**2*opr(1,2)**2
     .        + (azzneutneut-vzzneutneut)**2*opl(1,2)**2)*x1*y1
     .        -4d0*xmchar(1)/xmchar(2)*opl(1,2)*opr(1,2)
     .        *(azzneutneut**2+vzzneutneut**2)*y3 )
         endif
      else
         charz=0d0
      endif
c -------------------------------------------------------------------- c
c                      Z-selectron interference
c -------------------------------------------------------------------- c
      charzsel=0d0
c
      if(amchar(1).le.amchar(2)) then
         do i=1,2
            if (xmusel(i).gt.1.0d0.and.(mz+amchar(1)).gt.amchar(2))then
            charzsel=charzsel-g2s**2*4d0/dsel(i)/dz/cw
     .           *((ale(i,2)*ale(i,1)*opr(1,2)*
     .           (azzneutneut+vzzneutneut)+
     .           ble(i,2)*ble(i,1)*opl(1,2)*
     .           (-azzneutneut+vzzneutneut))*x1*y1
     .           -(ale(i,1)*ale(i,2)*opl(1,2)*
     .           (azzneutneut+vzzneutneut)
     .           +ble(i,1)*ble(i,2)*opr(1,2)*
     .           (-azzneutneut+vzzneutneut)
     .           )*xmchar(1)/xmchar(2)*y3)
            endif
         enddo
      else
         charzsel=0d0
      endif
      NS_charnue = charsel+charz+charzsel
      end
c
c ==================================================================== c
c ======================= chargino1 nu_m nubar_m ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_charnum(x1,x2)
c
      IMPLICIT NONE
      INTEGER I,k
      DOUBLE PRECISION dsel(2),blmu(2,2)
      DOUBLE PRECISION opl(2,2),opr(2,2),onl(5,5),onr(5,5)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION xmuneut1,X1,X2,X3,charsel,charz,charzsel,y1,
     .y2,y3,xmusel(2),xmuz,dz
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
c
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_coup4/opl,opr,onl,onr
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
c
      xmuneut1 = amchar(1)**2/amchar(2)**2
c
      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c
      do i=1,2
         blmu(1,i) = 0d0
         blmu(2,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c                          selectron exchange
c -------------------------------------------------------------------- c
      xmusel(1) = asmu1**2/amchar(2)**2
      xmusel(2) = asmu2**2/amchar(2)**2
c
      dsel(1)  = 1d0-x1-xmusel(1)
      dsel(2)  = 1d0-x1-xmusel(2)
c
      charsel=0d0
c
      if (amchar(1).le.amchar(2)) then
         do i=1,2
            do k=1,2
      if (xmusel(k).gt.1.0d0.and.xmusel(i).gt.1.0d0)then
               charsel=charsel+g2s**2/dsel(k)/dsel(i)*x1*y1*
     .         (almu(k,2)*almu(i,2)+blmu(k,2)*blmu(i,2))*
     .         (almu(k,1)*almu(i,1)+blmu(k,1)*blmu(i,1))
      endif
            enddo
         enddo
      else
         charsel=0d0
      endif
c -------------------------------------------------------------------- c
c                             Z exchange
c -------------------------------------------------------------------- c
      xmuz  = mz**2/amchar(2)**2
      dz    = y3-xmuz
c
      charz=0d0
c
      if (amchar(1).le.amchar(2)) then
         if ((mz+amchar(1)).gt.amchar(2))then
         charz=charz+g2s**2*4d0/dz**2/cw**2*
     .        (((azzneutneut+vzzneutneut)**2*opl(1,2)**2
     .        + (azzneutneut-vzzneutneut)**2*opr(1,2)**2)*x2*y2
     .        +((azzneutneut+vzzneutneut)**2*opr(1,2)**2
     .        + (azzneutneut-vzzneutneut)**2*opl(1,2)**2)*x1*y1
     .        -4d0*xmchar(1)/xmchar(2)*opl(1,2)*opr(1,2)
     .        *(azzneutneut**2+vzzneutneut**2)*y3 )
         endif
      else
         charz=0d0
      endif
c -------------------------------------------------------------------- c
c                      Z-selectron interference
c -------------------------------------------------------------------- c
      charzsel=0d0
c
      if(amchar(1).le.amchar(2)) then
         do i=1,2
            if (xmusel(i).gt.1.0d0.and.(mz+amchar(1)).gt.amchar(2))then
            charzsel=charzsel-g2s**2*4d0/dsel(i)/dz/cw
     .           *((almu(i,2)*almu(i,1)*opr(1,2)*
     .           (azzneutneut+vzzneutneut)+
     .           blmu(i,2)*blmu(i,1)*opl(1,2)*
     .           (-azzneutneut+vzzneutneut))*x1*y1
     .           -(almu(i,1)*almu(i,2)*opl(1,2)*
     .           (azzneutneut+vzzneutneut)
     .           +blmu(i,1)*blmu(i,2)*opr(1,2)*
     .           (-azzneutneut+vzzneutneut)
     .           )*xmchar(1)/xmchar(2)*y3)
            endif
         enddo
      else
         charzsel=0d0
      endif
      NS_charnum = charsel+charz+charzsel
      end
c
c ==================================================================== c
c ===================== chargino1 nu_tau nubar_tau =================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_charnutau(x1,x2)
c
      IMPLICIT NONE
      INTEGER i,k
      DOUBLE PRECISION blto(2,2),dsl(2)
      DOUBLE PRECISION opl(2,2),opr(2,2),onl(5,5),onr(5,5)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
c
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION xmuneut1,x1,x2,x3,y1,y2,y3,xmusl(2),
     .charstau,xmuz,dz,charz,charzstau
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
c
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_coup4/opl,opr,onl,onr
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
c
      xmuneut1 = amchar(1)**2/amchar(2)**2
c
      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c
      do i=1,2
         blto(1,i) = 0d0
         blto(2,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c                           stau exchange
c -------------------------------------------------------------------- c
      xmusl(1) = astau1**2/amchar(2)**2
      xmusl(2) = astau2**2/amchar(2)**2
c
      dsl(1)  = 1d0-x1-xmusl(1)
      dsl(2)  = 1d0-x1-xmusl(2)
c
      charstau=0d0
c
      if (amchar(1).le.amchar(2)) then
         do i=1,2
            do k=1,2
      if (xmusl(k).gt.1.0d0.and.xmusl(i).gt.1.0d0)then
               charstau=charstau+g2s**2/dsl(k)/dsl(i)*x1*y1*
     .         (altau(k,2)*altau(i,2)+blto(k,2)*blto(i,2))*
     .         (altau(k,1)*altau(i,1)+blto(k,1)*blto(i,1))
      endif
            enddo
         enddo
      else
         charstau=0d0
      endif
c
c -------------------------------------------------------------------- c
c                             Z exchange
c -------------------------------------------------------------------- c
      xmuz  = mz**2/amchar(2)**2
      dz    = y3-xmuz
c
      charz=0d0
c
      if (amchar(1).le.amchar(2)) then
         if ((mz+amchar(1)).gt.amchar(2))then
         charz=charz+g2s**2*4d0/dz**2/cw**2*
     .        (((azzneutneut+vzzneutneut)**2*opl(1,2)**2
     .        + (azzneutneut-vzzneutneut)**2*opr(1,2)**2)*x2*y2
     .        +((azzneutneut+vzzneutneut)**2*opr(1,2)**2
     .        + (azzneutneut-vzzneutneut)**2*opl(1,2)**2)*x1*y1
     .        -4d0*xmchar(1)/xmchar(2)*opl(1,2)*opr(1,2)
     .        *(azzneutneut**2+vzzneutneut**2)*y3 )
         endif
      else
         charz=0d0
      endif
c -------------------------------------------------------------------- c
c                        Z-stau interference
c -------------------------------------------------------------------- c
      charzstau=0d0
c
      if(amchar(1).le.amchar(2)) then
         do i=1,2
            if (xmusl(i).gt.1.0d0.and.(mz+amchar(1)).gt.amchar(2))then
            charzstau=charzstau-g2s**2*4d0/dsl(i)/dz/cw
     .           *((altau(i,2)*altau(i,1)*opr(1,2)*
     .           (azzneutneut+vzzneutneut)+
     .           blto(i,2)*blto(i,1)*opl(1,2)*
     .           (-azzneutneut+vzzneutneut))*x1*y1
     .           -(altau(i,1)*altau(i,2)*opl(1,2)*
     .           (azzneutneut+vzzneutneut)
     .           +blto(i,1)*blto(i,2)*opr(1,2)*
     .           (-azzneutneut+vzzneutneut)
     .           )*xmchar(1)/xmchar(2)*y3)
            endif
         enddo
      else
         charzstau=0d0
      endif
c
      NS_charnutau = charstau+charz+charzstau
c
      end
c ==================================================================== c
c ======================== chargino1 up upbar ======================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_charup(x1,x2)
c
      IMPLICIT NONE
      INTEGER i,k
      DOUBLE PRECISION dsd(2),bldo(2,2)
      DOUBLE PRECISION opl(2,2),opr(2,2),onl(5,5),onr(5,5)
      DOUBLE PRECISION alup(2,2),aldo(2,2)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION xmuneut1,X1,x2,x3,y1,y2,y3,xmusd(2)
     .,charsdow,xmuz,dz,charz,charzsdow
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
c
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_coup4/opl,opr,onl,onr
      COMMON/NS_coup7/alup,aldo
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
c
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
c
      xmuneut1 = amchar(1)**2/amchar(2)**2
c
      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c
      do i=1,2
         bldo(1,i) = 0d0
         bldo(2,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c                            sdown exchange
c -------------------------------------------------------------------- c
      xmusd(1) = asdown1**2/amchar(2)**2
      xmusd(2) = asdown2**2/amchar(2)**2
c
      dsd(1)  = 1d0-x1-xmusd(1)
      dsd(2)  = 1d0-x1-xmusd(2)
c
      charsdow=0d0
c
      if (amchar(1).le.amchar(2)) then
         do i=1,2
            do k=1,2
               if (xmusd(k).gt.1d0.and.xmusd(i).gt.1d0)then
               charsdow=charsdow+g2s**2/dsd(k)/dsd(i)*x1*y1*
     .         (aldo(k,2)*aldo(i,2)+bldo(k,2)*bldo(i,2))*
     .         (aldo(k,1)*aldo(i,1)+bldo(k,1)*bldo(i,1))
               endif
            enddo
         enddo
      else
         charsdow=0d0
      endif
c
c -------------------------------------------------------------------- c
c                             Z exchange
c -------------------------------------------------------------------- c
      xmuz = mz**2/amchar(2)**2
      dz   = y3-xmuz
c
      charz=0d0
c
      if (amchar(1).le.amchar(2)) then
         if ((mz+amchar(1)).gt.amchar(2))then
         charz=charz+g2s**2*4d0/dz**2/cw**2*
     .        (((azztoptop+vzztoptop)**2*opl(1,2)**2
     .        + (azztoptop-vzztoptop)**2*opr(1,2)**2)*x2*y2
     .        +((azztoptop+vzztoptop)**2*opr(1,2)**2
     .        + (azztoptop-vzztoptop)**2*opl(1,2)**2)*x1*y1
     .        -4d0*dsqrt(xmuneut1)*opl(1,2)*opr(1,2)
     .        *(azztoptop**2+vzztoptop**2)*y3 )
         endif
      else
         charz=0d0
      endif
c -------------------------------------------------------------------- c
c                        Z-sdown interference
c -------------------------------------------------------------------- c
      charzsdow=0d0
c
      if(amchar(1).le.amchar(2)) then
         do i=1,2
             if (xmusd(i).gt.1d0.and.(mz+amchar(1)).gt.amchar(2))then
            charzsdow=charzsdow-g2s**2*4d0/dsd(i)/dz/cw
     .           *((aldo(i,2)*aldo(i,1)*
     .           opr(1,2)*(azztoptop+vzztoptop)+
     .           bldo(i,2)*bldo(i,1)*opl(1,2)*
     .           (-azztoptop+vzztoptop))*x1*y1
     .           -(aldo(i,1)*aldo(i,2)*opl(1,2)*(azztoptop+vzztoptop)
     .           +bldo(i,1)*bldo(i,2)*opr(1,2)*(-azztoptop+vzztoptop)
     .           )*dsqrt(xmuneut1)*y3)
            endif
         enddo
      else
         charzsdow=0d0
      endif
c
      NS_charup = charsdow+charz+charzsdow
c
      end
c ==================================================================== c
c ====================== chargino1 down downbar ====================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_chardow(x1,x2)
c
      IMPLICIT NONE
      INTEGER i,k
      DOUBLE PRECISION dsup(2),blup(2,2)
      DOUBLE PRECISION opl(2,2),opr(2,2),onl(5,5),onr(5,5)
      DOUBLE PRECISION alup(2,2),aldo(2,2)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION xmuneut1,x1,x2,x3,y1,y2,y3,xmusup(2),
     .charsup,xmuz,dz,charz,charzsup
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
c
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_coup4/opl,opr,onl,onr
      COMMON/NS_coup7/alup,aldo
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
c
      xmuneut1 = amchar(1)**2/amchar(2)**2
c
      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c
      do i=1,2
         blup(1,i) = 0d0
         blup(2,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c                            sup exchange
c -------------------------------------------------------------------- c
      xmusup(1) = asup1**2/amchar(2)**2
      xmusup(2) = asup2**2/amchar(2)**2
c
      dsup(1)  = 1-x2-xmusup(1)
      dsup(2)  = 1-x2-xmusup(2)
c
      charsup=0d0
c
      if (amchar(1).le.amchar(2)) then
         do i=1,2
            do k=1,2
      if (xmusup(k).gt.1d0.and.xmusup(i).gt.1d0)then
               charsup=charsup+g2s**2/dsup(k)/dsup(i)*x2*y2*
     .         (alup(k,2)*alup(i,2)+blup(k,2)*blup(i,2))*
     .         (alup(k,1)*alup(i,1)+blup(k,1)*blup(i,1))
      endif
            enddo
         enddo
      else
         charsup=0d0
      endif
c
c -------------------------------------------------------------------- c
c                             Z exchange
c -------------------------------------------------------------------- c
      xmuz = mz**2/amchar(2)**2
      dz   = y3-xmuz
c
      charz=0d0
c
      if (amchar(1).le.amchar(2)) then
         if((mz+amchar(1)).gt.amchar(2))then
         charz=charz+g2s**2*4d0/dz**2/cw**2*
     .        (((azzbotbot+vzzbotbot)**2*opl(1,2)**2
     .        + (azzbotbot-vzzbotbot)**2*opr(1,2)**2)*x2*y2
     .        +((azzbotbot+vzzbotbot)**2*opr(1,2)**2
     .        + (azzbotbot-vzzbotbot)**2*opl(1,2)**2)*x1*y1
     .        -4d0*dsqrt(xmuneut1)*opl(1,2)*opr(1,2)
     .        *(azzbotbot**2+vzzbotbot**2)*y3 )
         endif
      else
         charz=0d0
      endif
c -------------------------------------------------------------------- c
c                         Z-sup interference
c -------------------------------------------------------------------- c
      charzsup=0d0
c
      if(amchar(1).le.amchar(2)) then
         do i=1,2
      if (xmusup(i).gt.1d0.and.(mz+amchar(1)).gt.amchar(2))then
            charzsup=charzsup+g2s**2*4d0/dsup(i)/dz/cw
     .           *((alup(i,2)*alup(i,1)*opl(1,2)*(azzbotbot+vzzbotbot)+
     .           blup(i,2)*blup(i,1)*opr(1,2)*
     .           (-azzbotbot+vzzbotbot))*x2*y2
     .           -(alup(i,1)*alup(i,2)*opr(1,2)*(azzbotbot+vzzbotbot)
     .           +blup(i,1)*blup(i,2)*opl(1,2)*(-azzbotbot+vzzbotbot)
     .           )*dsqrt(xmuneut1)*y3)
            endif
         enddo
      else
         charzsup=0d0
      endif
c
      NS_chardow = charsup+charz+charzsup
c
      end
c ==================================================================== c
c ==================== chargino1 bottom bottombar ==================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_charbot(x1,x2)
c
      IMPLICIT NONE
      INTEGER ni,nj,i,J,k
      DOUBLE PRECISION dstop(2),sgn(2)
      DOUBLE PRECISION opl(2,2),opr(2,2),onl(5,5),onr(5,5)
      DOUBLE PRECISION alstor(2,2),akstor(2,2)
      DOUBLE PRECISION alstor1(2,2),akstor1(2,2),
     .alstor2(2,2),akstor2(2,2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION Hbbr(3),Abbr(2)
      DOUBLE PRECISION hchichi(3,5,5),achichi(2,5,5),
     .hchachaR(3,2,2),hchachaL(3,2,2),achachaR(2,2,2),
     .     achachaL(2,2,2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION xmust(2),xmub,UH,TH,X1,X2,X3,
     .vzz,azz,charstop,charzbot,y3,xmuneut1,xmuz,dz,
     .rh,sh,rk,charzstop,hv1,hv2,hv3,hv4,hv5,hv6,hv7,hv8,
     .charhstop(3),charhastop(2),charzh(3),charza(2),charhabot(2,2)
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION charhbot(3,3),dhl(3),dhh(3),dha(2),dhna(2)
c
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup4/opl,opr,onl,onr
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coupNMSSM/hchichi,achichi,hchachaR,hchachaL,achachaR,
     .     achachaL
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_charstopbot/alstor,akstor
      COMMON/NS_phibotbot/Hbbr,Abbr
c
      do i=1,2
         do j=1,2
            alstor1(i,j)=alstor(i,j)
            akstor1(i,j)=akstor(i,j)
            alstor2(i,j)=alstor(i,j)
            akstor2(i,j)=akstor(i,j)
         enddo
      enddo
c
      do i=1,2
         sgn(i) = 1d0
         if(xmchar(i).ge.0d0) then
            sgn(i) = 1d0
         elseif(xmchar(i).lt.0d0) then
            sgn(i) = -1d0
         endif
      enddo
c
      xmuneut1 = amchar(1)**2/amchar(2)**2
      xmust(1)   = ast1**2/amchar(2)**2
      xmust(2)   = ast2**2/amchar(2)**2
      xmub     = mbp**2/amchar(2)**2
c
      x3 = 2d0-x1-x2
      y3 = 1d0+xmuneut1-x3
c
      uh = 1d0-x1+xmub
      th = 1d0-x2+xmub
c
      vzz = vzzbotbot
      azz = azzbotbot
c -------------------------------------------------------------------- c
c                            stop exchange
c -------------------------------------------------------------------- c
      dstop(1)  = 1d0-x2-xmust(1)+xmub
      dstop(2)  = 1d0-x2-xmust(2)+xmub
c
      charstop=0d0
c
      if((amchar(1)+2d0*mbp).le.amchar(2)) then
         do i=1,2
            do k=1,2
               if (xmust(k).gt.1d0.and.xmust(i).gt.1d0)then
               charstop=charstop
     .           +g2s**2/dstop(k)/dstop(i)*(
     .           (alstor1(i,2)*akstor1(k,2)+akstor1(i,2)*alstor1(k,2))*
     .           (alstor2(i,1)*akstor2(k,1)+akstor2(i,1)*alstor2(k,1))*
     .           xmchar(1)/xmchar(2)*xmub*(-4d0)+
     .           (alstor1(i,2)*alstor1(k,2)+akstor1(i,2)*akstor1(k,2))*
     .           (alstor2(i,1)*akstor2(k,1)+akstor2(i,1)*alstor2(k,1))*
     .           xmchar(1)/xmchar(2)*dsqrt(xmub)*2d0*
     .           (th-xmub-1d0)+
     .           (alstor1(i,2)*akstor1(k,2)+akstor1(i,2)*alstor1(k,2))*
     .           (alstor2(i,1)*alstor2(k,1)+akstor2(i,1)*akstor2(k,1))*
     .           dsqrt(xmub)*2d0*(th-xmub-xmuneut1)+
     .           (alstor1(i,2)*alstor1(k,2)+akstor1(i,2)*akstor1(k,2))*
     .           (alstor2(i,1)*alstor2(k,1)+akstor2(i,1)*akstor2(k,1))*
     .           (-th**2+th*(1d0+xmuneut1+2d0*xmub)-
     .           (xmuneut1+xmub)*(1d0+xmub)))
               endif
            enddo
         enddo
      else
         charstop=0d0
      endif

c -------------------------------------------------------------------- c
c                            Z exchange
c -------------------------------------------------------------------- c
      xmuz = mz**2/amchar(2)**2
      dz   = y3-xmuz

      charzbot = 0d0

      rh = xmuneut1+2d0*xmub-th-uh+1d0
      sh = (xmuneut1-th-uh+1d0)*2d0*xmub+4d0*xmub**2
      rk = xmuneut1*(2d0*xmub-th-uh+4d0)+2d0*xmub-uh-th

      if((amchar(1)+2d0*mbp).le.amchar(2)) then
         if((mz+amchar(1)).gt.amchar(2))then
         charzbot=charzbot+g2s**2/dz**2/cw**2*(
     .    opl(1,2)*opr(1,2)*(vzz**2-azz**2)*
     .    xmchar(1)/xmchar(2)*xmub*(-16d0/xmuz**2*rh**2+
     .    32d0/xmuz*rh-64d0)+
     .    opl(1,2)*opr(1,2)*(vzz**2+azz**2)*
     .    xmchar(1)/xmchar(2)*(8d0/xmuz**2*rh*sh-16d0/xmuz*sh
     .    -16d0*(xmuneut1-uh-th+1d0))+
     .    (opl(1,2)**2+opr(1,2)**2)*(vzz**2-azz**2)*
     .    xmub*(4d0/xmuz**2*rh*rk-8d0/xmuz*rk+8d0*(uh+th-
     .    2d0*xmub))
     .    +(opl(1,2)**2+opr(1,2)**2)*(vzz**2+azz**2)*
     .    (-2d0/xmuz**2*rk*sh+8d0/xmuz*(xmuneut1*(2d0*xmub**2+
     .    4d0*xmub-xmub*(th+uh))+2d0*xmub**2-xmub*(uh+th))+
     .    4d0*(xmuneut1*(uh+th-2d0*xmub-2d0)+2d0*xmub*(uh+th-
     .    1d0)-2d0*xmub**2+th*(-th+1d0)+uh*(-uh+1d0)))+
     .    (opl(1,2)**2-opr(1,2)**2)*vzz*azz*8d0*(
     .    xmuneut1*(th-uh)+2d0*xmub*(th-uh)+th*(-th+1d0)+uh*(uh-
     .    1d0)))
         endif
      else
         charzbot=0d0
      endif
c -------------------------------------------------------------------- c
c         NMSSM CP ODD Higgs exchange+ Interference A(i)-A(j)
c -------------------------------------------------------------------- c
      do i=1,2
         do j=1,2

       dha(i)  = y3-pmass(i)**2/amchar(2)**2
       dhna(j)  = y3-pmass(j)**2/amchar(2)**2
       charhabot(i,j)=0d0
       if (pmass(i)+amchar(1).gt.amchar(2).and.
     .(pmass(j)+amchar(1)).gt.amchar(2))then
      if((amchar(1)+2d0*mbp).le.amchar(2)) then

         charhabot(i,j)=g2s**2/dha(i)/dhna(j)*Abbr(i)*Abbr(j)*(
     .    achachaL(i,1,2)*achachaR(j,1,2)
     .    *xmchar(1)/xmchar(2)*xmub*8d0+
     .    achachaL(i,1,2)*achachaR(j,1,2)*xmchar(1)/xmchar(2)*4d0*
     .    (1d0+xmuneut1-th-uh)+
     .    (achachaL(i,1,2)*achachaL(j,1,2)+achachaR(i,1,2)*
     .    achachaR(j,1,2))*xmub*2d0*
     .     (-2d0*xmub+th+uh)+
     .    (achachaL(i,1,2)*achachaL(j,1,2)+achachaR(i,1,2)*
     .     achachaR(j,1,2))
     .     *(xmuneut1*(uh+th-2d0*xmub)
     .    +2d0*xmub*(uh+th-1d0)-(th+uh)**2+uh+th))
      else
         charhabot(i,j)=0d0
      endif
      endif

      enddo
      enddo
c -------------------------------------------------------------------- c
c            NMSSM CP EVEN Higgs exchange + Interference H(i)-H(j)
c -------------------------------------------------------------------- c
        do i=1,3
        do j=1,3

      dhl(i)   = y3-smass(i)**2/amchar(2)**2
      dhh(j)   = y3-smass(j)**2/amchar(2)**2

         charhbot(i,j) =0d0
         if (smass(i)+amchar(1).gt.amchar(2)
     ..and.(smass(j)+amchar(1)).gt.amchar(2))then
      if((amchar(1)+2d0*mbp).le.amchar(2)) then

        charhbot(i,j) = g2s**2/dhh(j)/dhl(i)*hbbr(i)*hbbr(j)*(
     .   (hchachaL(i,1,2)*hchachaR(j,1,2)+
     .    hchachaR(i,1,2)*hchachaL(j,1,2))*(
     .   xmchar(1)/xmchar(2)*xmub*(-4d0) +
     .   xmchar(1)/xmchar(2)*2d0*(1d0+xmuneut1-uh-th)) +
     .   (hchachaL(i,1,2)*hchachaL(j,1,2)+
     .    hchachaR(i,1,2)*hchachaR(j,1,2))*(
     .   xmub*2d0*(2d0*xmub-th-uh)+
     .   (xmuneut1*(uh+th-2d0*xmub)+2d0*xmub*(uh+th-1d0)
     .   +th+uh-(th+uh)**2)))
      else
         charhbot(i,j) =0d0
      endif
      endif
        enddo
        enddo
c -------------------------------------------------------------------- c
c                       Z stop interference
c -------------------------------------------------------------------- c
      charzstop=0d0

      if((amchar(1)+2d0*mbp).le.amchar(2)) then
         do i=1,2
            hv1 = alstor1(i,2)*alstor2(i,1)*opr(2,1)*(vzz-azz)
     .           +akstor1(i,2)*akstor2(i,1)*opl(2,1)*(vzz+azz)
            hv2 = alstor1(i,2)*alstor2(i,1)*opr(2,1)*(vzz+azz)
     .           +akstor1(i,2)*akstor2(i,1)*opl(2,1)*(vzz-azz)
            hv3 = alstor1(i,2)*alstor2(i,1)*opl(2,1)*(vzz-azz)
     .           +akstor1(i,2)*akstor2(i,1)*opr(2,1)*(vzz+azz)
            hv4 = alstor1(i,2)*alstor2(i,1)*opl(2,1)*(vzz+azz)
     .           +akstor1(i,2)*akstor2(i,1)*opr(2,1)*(vzz-azz)

            hv5 = alstor2(i,1)*akstor1(i,2)*opr(2,1)*(vzz+azz)
     .           +akstor2(i,1)*alstor1(i,2)*opl(2,1)*(vzz-azz)
            hv6 = alstor2(i,1)*akstor1(i,2)*opr(2,1)*(vzz-azz)
     .           +akstor2(i,1)*alstor1(i,2)*opl(2,1)*(vzz+azz)
            hv7 = alstor2(i,1)*akstor1(i,2)*opl(2,1)*(vzz+azz)
     .           +akstor2(i,1)*alstor1(i,2)*opr(2,1)*(vzz-azz)
            hv8 = alstor2(i,1)*akstor1(i,2)*opl(2,1)*(vzz-azz)
     .           +akstor2(i,1)*alstor1(i,2)*opr(2,1)*(vzz+azz)
            if (xmust(i).gt.1d0.and.(mz+amchar(1)).gt.amchar(2))then
            charzstop=charzstop
     .      -g2s**2/dstop(i)/dz/cw*(
     .      hv1*xmchar(1)/xmchar(2)*xmub*(
     .      1/xmuz*(-4d0)*(1d0+xmuneut1+2d0*xmub-uh-th)+16d0) +
     .      hv2*xmchar(1)/xmchar(2)*(2d0/xmuz*(2d0*xmub*(xmuneut1+
     .      1d0-th-uh)+4d0*xmub**2)+4d0*(1d0+xmuneut1-uh-th)) +
     .      hv3*xmub*(2d0/xmuz*(xmuneut1*(2d0*xmub-th-uh+4d0)+
     .      2d0*xmub-th-uh)+4d0*(2d0*xmub-th-uh)) +
     .      hv4*(2d0/xmuz*(xmuneut1*(-2d0*xmub**2+xmub*th
     .      -2d0*xmub+xmub*uh-2d0*xmub)+xmub*(-2d0*xmub
     .      +uh)+xmub*th)+4d0*(xmuneut1*(xmub-th+1d0)+xmub*
     .      (xmub-th)+xmub*(1d0-th)+th**2-th)) +
     .      hv5*xmchar(1)/xmchar(2)*dsqrt(xmub)*
     .      (2d0/xmuz*(xmuneut1*(1d0+xmub-th)+1d0+xmub*
     .      (2d0*xmub-2d0*th-uh+3d0)+th*(th-xmub+uh-2d0)-uh)+
     .      4d0*(1d0+xmub-uh)) +
     .      hv6*xmchar(1)/xmchar(2)*dsqrt(xmub)*
     .      (2d0/xmuz*(xmuneut1*(-1d0-xmub+th)-1d0+xmub*
     .      (-xmub+th-1d0)+xmub*(2d0*th+uh-2d0-xmub)
     .      -th*(th+uh)+2d0*th+uh)+8d0*(1d0+xmub-th)) +
     .      hv7*dsqrt(xmub)*(2d0/xmuz*(xmuneut1*(-xmuneut1
     .      -3d0*xmub+2d0*th+uh-1d0)+xmub*(-2d0*xmub+2d0*th+
     .      uh-1d0)+th*(xmub-th-uh+1d0))+8d0*(xmuneut1+xmub-th)) +
     .      hv8*dsqrt(xmub)*(2d0/xmuz*(xmuneut1*(xmuneut1+
     .      3d0*xmub-uh-2d0*th+1d0)+xmub*(xmub-th+1d0)+
     .      xmub*(xmub-2d0*th-uh)+uh*th+th**2-th)+4d0*(xmuneut1
     .      +xmub-uh)) )
            endif
         enddo
      else
         charzstop=0d0
      endif
c -------------------------------------------------------------------- c
c                         Higgs-stop interference
c -------------------------------------------------------------------- c
      do j=1,3
      dhl(j)   = y3-smass(j)**2/amchar(2)**2
      charhstop(j)=0d0

      if ((amchar(1)+2d0*mbp).le.amchar(2)) then
         do i=1,2
      if ((smass(j)+amchar(1)).gt.amchar(2).and.xmust(i).gt.1d0)then
            charhstop(j)=charhstop(j)
     .       +2d0*g2s**2/dhl(j)/dstop(i)*(Hbbr(j)/dsqrt(2d0))*(
     .       (alstor2(i,1)*alstor1(i,2)*hchachaR(j,1,2)+
     .        akstor1(i,2)*akstor2(i,1)*hchachaL(j,1,2))*(
     .       dsqrt(xmub)*(uh-xmuneut1-xmub) +
     .       dsqrt(xmub)*(-th+xmub+xmuneut1) ) +
     .       (alstor2(i,1)*alstor1(i,2)*hchachaL(j,1,2)+
     .        akstor1(i,2)*akstor2(i,1)*hchachaR(j,1,2))*(
     .       dsqrt(xmub)*xmchar(1)/xmchar(2)*(-th+1d0+xmub) +
     .       dsqrt(xmub)*xmchar(1)/xmchar(2)*(uh-1d0-xmub) ) +
     .       (alstor2(i,1)*akstor1(i,2)*hchachaL(j,1,2)
     .       +alstor1(i,2)*akstor2(i,1)*hchachaR(j,1,2))*(
     .       2d0*xmub*xmchar(1)/xmchar(2) +
     .       xmchar(1)/xmchar(2)*(uh+th-xmuneut1-1d0) ) +
     .       (alstor2(i,1)*akstor1(i,2)*hchachaR(j,1,2)
     .       +alstor1(i,2)*akstor2(i,1)*hchachaL(j,1,2))*(
     .       (uh*th+th**2-th*(1d0+2d0*xmub+xmuneut1)+xmub+
     .       xmub*xmuneut1) +
     .       xmub*(uh+th-2d0*xmub)) )
            endif
         enddo
      else
         charhstop(j)=0d0
      endif

      enddo
c -------------------------------------------------------------------- c
c                       A-stop interference
c -------------------------------------------------------------------- c
      do j=1,2
      dha(j) = y3-pmass(j)**2/amchar(2)**2

      charhastop(j)=0d0

      if((amchar(1)+2d0*mbp).le.amchar(2)) then

         do i=1,2
      if (pmass(j)+amchar(1).gt.amchar(2).and.xmust(i).gt.1d0)then
            charhastop(j)=charhastop(j)
     .       +2d0*g2s**2/dha(j)/dstop(i)*(Abbr(j)/dsqrt(2d0))*(
     .       (alstor2(i,1)*alstor1(i,2)*achachaR(j,1,2)-
     .        akstor1(i,2)*akstor2(i,1)*achachaL(j,1,2))*(
     .       dsqrt(xmub)*(uh-xmuneut1-xmub) +
     .       dsqrt(xmub)*(-th+xmub+xmuneut1)*(-1d0) ) +
     .       (alstor2(i,1)*alstor1(i,2)*achachaL(j,1,2)-
     .        akstor1(i,2)*akstor2(i,1)*achachaR(j,1,2))*(
     .       dsqrt(xmub)*xmchar(1)/xmchar(2)*(-th+1d0+xmub) +
     .       dsqrt(xmub)*xmchar(1)/xmchar(2)*(uh-1d0-xmub)*(-1d0)
     .       ) +
     .       (alstor2(i,1)*akstor1(i,2)*achachaL(j,1,2)
     .       -alstor1(i,2)*akstor2(i,1)*achachaR(j,1,2))*(
     .       2d0*xmub*xmchar(1)/xmchar(2) +
     .       xmchar(1)/xmchar(2)*(uh+th-xmuneut1-1d0)*(-1d0) ) +
     .       (alstor2(i,1)*akstor1(i,2)*achachaR(j,1,2)
     .       -alstor1(i,2)*akstor2(i,1)*achachaL(j,1,2))*(
     .       (uh*th+th**2-th*(1d0+2d0*xmub+xmuneut1)+xmub+
     .       xmub*xmuneut1)*(-1d0) +
     .       xmub*(uh+th-2d0*xmub)) )
            endif
         enddo
      else
         charhastop(j)=0d0
      endif

      enddo
c -------------------------------------------------------------------- c
c 	                interference Z and H/h/NH,A,NA
c -------------------------------------------------------------------- c
      do i=1,3
      dhl(i)   = y3-smass(i)**2/amchar(2)**2

      charzh(i)=0d0

      if((amchar(1)+2d0*mbp).le.amchar(2)) then
         if ((smass(i)+amchar(1)).gt.amchar(2).and.
     .(mz+amchar(1)).gt.amchar(2))then
         charzh(i)=
     .    -2d0*g2s**2/cw/dhl(i)/dz*vzz*Hbbr(i)/dsqrt(2d0)*(
     .    xmchar(1)/xmchar(2)*dsqrt(xmub)*
     .    (hchachaL(i,1,2)*opl(1,2)+hchachaR(i,1,2)*opr(1,2))*(
     .    4d0*(th-xmub-1d0) - 4d0*(uh-xmub-1d0) )
     .    +dsqrt(xmub)*(hchachaL(i,1,2)*opr(1,2)
     .     +hchachaR(i,1,2)*opl(1,2))*(
     .    4d0*(xmuneut1+xmub-uh)-4d0*(xmuneut1+xmub-th)) )
         endif
      else
         charzh(i)=0d0
      endif

      enddo
*
      do i=1,2
      dha(i)   = y3-pmass(i)**2/amchar(2)**2

      charza(i)=0d0

      if((amchar(1)+2d0*mbp).le.amchar(2)) then
         if ((pmass(i)+amchar(1)).gt.amchar(2).and.
     .(mz+amchar(1)).gt.amchar(2))then
         charza(i)=
     .    -2d0*g2s**2/cw/dha(i)/dz*azz*Abbr(i)/dsqrt(2d0)*(
     .    xmchar(1)/xmchar(2)*dsqrt(xmub)*
     .    (achachaL(i,1,2)*opl(1,2)+achachaR(i,1,2)*opr(1,2))*(
     .    4d0/xmuz*(2d0+xmuneut1*(2d0*xmub-th-uh+2d0)+xmub*
     .    (-uh-th+1d0+xmub)+xmub*(3d0*xmub-3d0*(th+uh)+5d0)
     .    +(th+uh)**2-3d0*(th+uh))+
     .    4d0*(th-xmub-1d0) +4d0*(uh-xmub-1d0))
     .    +dsqrt(xmub)*(achachaL(i,1,2)*opr(1,2)
     .     +achachaR(i,1,2)*opl(1,2))*(
     .    4d0/xmuz*(-2d0*xmuneut1**2+xmuneut1*(-6d0*xmub+
     .    3d0*(th+uh)-2d0)+xmub*(-xmub+th+uh-1d0)+xmub*
     .    (-3d0*xmub+3d0*(uh+th)-1d0)-(th+uh)**2+(th+uh))+
     .    4d0*(xmuneut1+xmub-uh)+4d0*(xmuneut1+xmub-th)) )
         endif
      else
         charza(i)=0d0
      endif

      enddo

      NS_charbot = charstop+charzbot+charzstop
        do i=1,2
          NS_charbot=NS_charbot+charhastop(i)+charza(i)
        do j=1,2
          NS_charbot=NS_charbot+charhabot(i,j)
        enddo
        enddo

        do i=1,3
          NS_charbot = NS_charbot+charhstop(i)
          NS_charbot = NS_charbot+charzh(i)
        do j=1,3
          NS_charbot = NS_charbot+charhbot(i,j)
        enddo
        enddo

      end
c ==================================================================== c
c ======================= chargino1 top topbar ======================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_chartop(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,i,j,k
      DOUBLE PRECISION dsbot(2)
      DOUBLE PRECISION opl(2,2),opr(2,2),onl(5,5),onr(5,5)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2),alsbot1(2,2),
     . aksbot1(2,2),alsbot2(2,2),aksbot2(2,2)
      DOUBLE PRECISION hchichi(3,5,5),achichi(2,5,5)
      DOUBLE PRECISION hchachaR(3,2,2),hchachaL(3,2,2),achachaR(2,2,2),
     .     achachaL(2,2,2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION Httr(3),Attr(2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION
     .UH,TH,X1,X2,X3,xmuneut1,xmusb(2),xmut,y3,vzz,azz,charsbot
     .,xmuz,dz,charztop,rh,sh,rk,dhl(3),dhh(3),
     .charhtop(3,3),dha(2),dhna(2),charzsbot,
     .hv1,hv2,hv3,hv4,hv5,hv6,hv7,hv8,charhsbot(3),
     .charhasbot(2),charzh(3),charza(2),charhatop(2,2),gmsb(2)
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup4/opl,opr,onl,onr
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coupNMSSM/hchichi,achichi,hchachaR,hchachaL,achachaR,
     .     achachaL
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_phitoptop/Httr,Attr
*
      do i=1,2
         do j=1,2
            alsbot1(i,j)=alsbot(i,j)
            aksbot1(i,j)=aksbot(i,j)
            alsbot2(i,j)=alsbot(i,j)
            aksbot2(i,j)=aksbot(i,j)
         enddo
      enddo
      gmsb(1) = asb1
      gmsb(2) = asb2
      xmuneut1 = amchar(1)**2/amchar(2)**2
      xmusb(1)   = asb1**2/amchar(2)**2
      xmusb(2)   = asb2**2/amchar(2)**2
      xmut     = mt**2/amchar(2)**2

      x3 = 2d0-x1-x2
      y3 = 1d0+xmuneut1-x3

      uh = 1d0-x1+xmut
      th = 1d0-x2+xmut

      vzz = vzztoptop
      azz = azztoptop
c -------------------------------------------------------------------- c
c                          sbottom exchange
c -------------------------------------------------------------------- c
      dsbot(1) = 1-x1-xmusb(1)+xmut
      dsbot(2) = 1-x1-xmusb(2)+xmut

      charsbot=0d0

      if((amchar(1)+2d0*mt).le.amchar(2)) then
         do i=1,2
            do k=1,2
      if((gmsb(k)+mt).gt.amchar(2).and.(gmsb(i)+mt).gt.amchar(2))then
               charsbot=charsbot
     .          +g2s**2/dsbot(k)/dsbot(i)*(
     .           (alsbot1(i,2)*aksbot1(k,2)+aksbot1(i,2)*alsbot1(k,2))*
     .           (alsbot2(i,1)*aksbot2(k,1)+aksbot2(i,1)*alsbot2(k,1))*
     .           xmchar(1)/xmchar(2)*xmut*(-4d0)+
     .           (alsbot1(i,2)*alsbot1(k,2)+aksbot1(i,2)*aksbot1(k,2))*
     .           (alsbot2(i,1)*aksbot2(k,1)+aksbot2(i,1)*alsbot2(k,1))*
     .           xmchar(1)/xmchar(2)*dsqrt(xmut)*2d0*
     .           (uh-xmut-1d0)+
     .           (alsbot1(i,2)*aksbot1(k,2)+aksbot1(i,2)*alsbot1(k,2))*
     .           (alsbot2(i,1)*alsbot2(k,1)+aksbot2(i,1)*aksbot2(k,1))*
     .           dsqrt(xmut)*2d0*(uh-xmut-xmuneut1)+
     .           (alsbot1(i,2)*alsbot1(k,2)+aksbot1(i,2)*aksbot1(k,2))*
     .           (alsbot2(i,1)*alsbot2(k,1)+aksbot2(i,1)*aksbot2(k,1))*
     .           (-uh**2+uh*(1d0+xmuneut1+2d0*xmut)-(xmuneut1+xmut)*
     .           (1d0+xmut)))
      endif
            enddo
         enddo
      else
         charsbot=0d0
      endif

c -------------------------------------------------------------------- c
c                            Z exchange
c -------------------------------------------------------------------- c
      xmuz = mz**2/amchar(2)**2
      dz   = y3-xmuz

      charztop=0d0

      rh = xmuneut1+2d0*xmut-th-uh+1d0
      sh = (xmuneut1-th-uh+1d0)*2d0*xmut+4d0*xmut**2
      rk = xmuneut1*(2d0*xmut-th-uh+4d0)+2d0*xmut-uh-th

      if((amchar(1)+2d0*mt).le.amchar(2)) then
         if((mz+amchar(1)).gt.amchar(2))then
         charztop=charztop+g2s**2/dz**2/cw**2*(
     .    opl(1,2)*opr(1,2)*(vzz**2-azz**2)*
     .    xmchar(1)/xmchar(2)*xmut*(-16d0/xmuz**2*rh**2+
     .    32d0/xmuz*rh-64d0)+
     .    opl(1,2)*opr(1,2)*(vzz**2+azz**2)*
     .    xmchar(1)/xmchar(2)*(8d0/xmuz**2*rh*sh-16d0/xmuz*sh
     .    -16d0*(xmuneut1-uh-th+1d0))+
     .    (opl(1,2)**2+opr(1,2)**2)*(vzz**2-azz**2)*
     .    xmut*(4d0/xmuz**2*rh*rk-8d0/xmuz*rk+8d0*(uh+th-2d0*xmut))
     .    +(opl(1,2)**2+opr(1,2)**2)*(vzz**2+azz**2)*
     .    (-2d0/xmuz**2*rk*sh+8d0/xmuz*(xmuneut1*(2d0*xmut**2+
     .    4d0*xmut-xmut*(th+uh))+2d0*xmut**2-xmut*(uh+th))+4d0*(
     .    xmuneut1*(uh+th-2d0*xmut-2d0)+2d0*xmut*(uh+th-1d0)
     .    -2d0*xmut**2+th*(-th+1d0)+uh*(-uh+1d0)))+
     .    (opl(1,2)**2-opr(1,2)**2)*vzz*azz*8d0*(
     .    xmuneut1*(th-uh)+2d0*xmut*(th-uh)+th*(-th+1d0)+uh*(uh-1d0))
     .    )
         endif
      else
         charztop=0d0
      endif
c -------------------------------------------------------------------- c
c    NMSSM CP ODD Higgs exchange+ Interference A(i)-A(j)
c -------------------------------------------------------------------- c
      do i=1,2
         do j=1,2
      dha(i)   = y3-pmass(i)**2/amchar(2)**2
      dhna(j)  = y3-pmass(j)**2/amchar(2)**2
      charhatop(i,j)=0d0
      if ((pmass(i)+amchar(1)).gt.amchar(2).and.
     .(pmass(j)+amchar(1)).gt.amchar(2))then
      if((amchar(1)+2d0*mt).le.amchar(2)) then

         charhatop(i,j)=g2s**2/dha(i)/dhna(j)*Attr(i)*Attr(j)/2d0*(
     .    achachaL(i,1,2)*achachaR(j,1,2)*(
     .    xmchar(1)/xmchar(2)*xmut*16d0+
     .    xmchar(1)/xmchar(2)*8d0*(1d0+xmuneut1-th-uh) )+
     .    (achachaL(i,1,2)*achachaL(j,1,2)+
     .        achachaR(i,1,2)*achachaR(j,1,2))*(
     .    xmut*(-4d0)*(2d0*xmut-th-uh)+
     .    2d0*(xmuneut1*(uh+th-2d0*xmut)+2d0*xmut*(uh+th-1d0)
     .    -(th+uh)**2+uh+th)) )
      else
         charhatop(i,j)=0d0
      endif
      endif

         enddo
      enddo
c -------------------------------------------------------------------- c
c         NMSSM CP EVEN Higgs exchange + Interference H(i)-H(j)
c -------------------------------------------------------------------- c
        do i=1,3
        do j=1,3

      dhl(i)   = y3-smass(i)**2/amchar(2)**2
      dhh(j)   = y3-smass(j)**2/amchar(2)**2

         charhtop(i,j) =0d0
      if ((smass(i)+amchar(1)).gt.amchar(2).and.
     .(smass(j)+amchar(1)).gt.amchar(2))then
      if((amchar(1)+2d0*mt).le.amchar(2)) then

        charhtop(i,j) = g2s**2/dhh(j)/dhl(i)*Httr(i)*Httr(j)*(
     .   (hchachaL(i,1,2)*hchachaR(j,1,2)+
     .    hchachaR(i,1,2)*hchachaL(j,1,2))*(
     .   xmchar(1)/xmchar(2)*xmut*(-4d0)+
     .   xmchar(1)/xmchar(2)*2d0*(1d0+xmuneut1-uh-th) )+
     .   (hchachaL(i,1,2)*hchachaL(j,1,2)+
     .    hchachaR(i,1,2)*hchachaR(j,1,2))*(
     .   xmut*2d0*(2d0*xmut-th-uh)+
     .   (xmuneut1*(uh+th-2d0*xmut)+2d0*xmut*(uh+th-1d0)+th+uh
     .   -(th+uh)**2)) )
      else
         charhtop(i,j)=0d0
      endif
      endif

        enddo
        enddo
c -------------------------------------------------------------------- c
c                      Z-sbottom interference
c -------------------------------------------------------------------- c
      charzsbot=0d0

      if((amchar(1)+2d0*mt).le.amchar(2)) then
         do i=1,2
            hv1 = alsbot1(i,2)*alsbot2(i,1)*opl(2,1)*(vzz-azz)
     .           +aksbot1(i,2)*aksbot2(i,1)*opr(2,1)*(vzz+azz)
            hv2 = alsbot1(i,2)*alsbot2(i,1)*opl(2,1)*(vzz+azz)
     .           +aksbot1(i,2)*aksbot2(i,1)*opr(2,1)*(vzz-azz)
            hv3 = alsbot1(i,2)*alsbot2(i,1)*opr(2,1)*(vzz-azz)
     .           +aksbot1(i,2)*aksbot2(i,1)*opl(2,1)*(vzz+azz)
            hv4 = alsbot1(i,2)*alsbot2(i,1)*opr(2,1)*(vzz+azz)
     .           +aksbot1(i,2)*aksbot2(i,1)*opl(2,1)*(vzz-azz)

            hv5 = alsbot2(i,1)*aksbot1(i,2)*opl(2,1)*(vzz+azz)
     .           +aksbot2(i,1)*alsbot1(i,2)*opr(2,1)*(vzz-azz)
            hv6 = alsbot2(i,1)*aksbot1(i,2)*opl(2,1)*(vzz-azz)
     .           +aksbot2(i,1)*alsbot1(i,2)*opr(2,1)*(vzz+azz)
            hv7 = alsbot2(i,1)*aksbot1(i,2)*opr(2,1)*(vzz+azz)
     .           +aksbot2(i,1)*alsbot1(i,2)*opl(2,1)*(vzz-azz)
            hv8 = alsbot2(i,1)*aksbot1(i,2)*opr(2,1)*(vzz-azz)
     .           +aksbot2(i,1)*alsbot1(i,2)*opl(2,1)*(vzz+azz)
            if((gmsb(i)+mt).gt.amchar(2).and.(mz+amchar(1)).gt.amchar(2)
     .)then
            charzsbot=charzsbot
     .      +g2s**2/dsbot(i)/dz/cw*(
     .      hv1*xmchar(1)/xmchar(2)*xmut*(
     .      1/xmuz*(-4d0)*(1d0+xmuneut1+2d0*xmut-uh-th)+16d0) +
     .      hv2*xmchar(1)/xmchar(2)*(2d0/xmuz*((xmuneut1+1d0-uh-th)*
     .      2d0*xmut+4d0*xmut**2)+4d0*(1d0+xmuneut1-uh-th)) +
     .      hv3*xmut*(2d0/xmuz*(xmuneut1*(2d0*xmut-th-uh+4d0)+2d0*
     .      xmut-th-uh)+4d0*(2d0*xmut-th-uh)) +
     .      hv4*(2d0/xmuz*(xmuneut1*(-2d0*xmut**2+xmut*th-2d0*xmut+
     .      xmut*uh-2d0*xmut)+xmut*(-2d0*xmut+uh)+xmut*th)+4d0*(
     .      xmuneut1*(xmut-uh+1d0)+xmut*(xmut-uh)+xmut*(1d0-uh)+uh**2
     .      -uh)) +
     .      hv5*xmchar(1)/xmchar(2)*dsqrt(xmut)*(2d0/xmuz*
     .      (xmuneut1*(1d0+xmut-uh)+1d0+xmut*(xmut-uh)+xmut*
     .      (xmut-th-2d0*uh+3d0)+th*(uh-1d0)+uh*(uh-2d0))+
     .      4d0*(1d0+xmut-th)) +
     .      hv6*xmchar(1)/xmchar(2)*dsqrt(xmut)*
     .      (2d0/xmuz*(xmuneut1*(-1d0-xmut+uh)-1d0+xmut*(-1d0+uh)
     .      +uh*(2d0-th-uh)+th+xmut*(-2d0*xmut+th+2d0*uh-2d0))+
     .      8d0*(1d0+xmut-uh)) +
     .      hv7*dsqrt(xmut)*(
     .      (-2d0)/xmuz*(xmuneut1-uh+xmut)*(xmuneut1+2d0*xmut-th-uh
     .      +1d0)+8d0*(xmuneut1+xmut-uh)) +
     .      hv8*dsqrt(xmut)*(
     .      2d0/xmuz*(xmuneut1*(xmuneut1+3d0*xmut-th-2d0*uh+1d0)+
     .      xmut*(2d0*xmut-th-2d0*uh)+xmut*(1d0-uh)+uh*th+uh**2-uh)+
     .      4d0*(xmuneut1+xmut-th)))
            endif
         enddo
      else
         charzsbot=0d0
      endif

c -------------------------------------------------------------------- c
c                      Higgs-sbottom interference
c -------------------------------------------------------------------- c

        do j=1,3

         dhl(j) = y3-smass(j)**2/amchar(2)**2

      charhsbot(j)=0d0

      if((amchar(1)+2d0*mt).le.amchar(2)) then
         do i=1,2
      if ((smass(j)+amchar(1)).gt.amchar(2).
     .and.(gmsb(i)+mt).gt.amchar(2))then
            charhsbot(j)=charhsbot(j)
     .       +2d0*g2s**2/dhl(j)/dsbot(i)*(Httr(j)/dsqrt(2d0))*(
     .       (alsbot1(i,2)*alsbot2(i,1)*hchachaL(j,1,2)+
     .        aksbot1(i,2)*aksbot2(i,1)*hchachaR(j,1,2))*(
     .       dsqrt(xmut)*(xmuneut1+xmut-uh) +
     .       dsqrt(xmut)*(-xmuneut1-xmut+th)) +
     .       (alsbot1(i,2)*alsbot2(i,1)*hchachaR(j,1,2)+
     .        aksbot1(i,2)*aksbot2(i,1)*hchachaL(j,1,2))*(
     .       xmchar(1)/xmchar(2)*dsqrt(xmut)*(-1d0-xmut+th) +
     .       xmchar(1)/xmchar(2)*dsqrt(xmut)*(1d0+xmut-uh) )
     .       +(alsbot2(i,1)*aksbot1(i,2)*hchachaR(j,1,2)+
     .         alsbot1(i,2)*aksbot2(i,1)*hchachaL(j,1,2))*(
     .       xmchar(1)/xmchar(2)*(-1d0-xmuneut1+uh+th) +
     .       2d0*xmchar(1)/xmchar(2)*xmut ) +
     .       (alsbot2(i,1)*aksbot1(i,2)*hchachaL(j,1,2)+
     .        alsbot1(i,2)*aksbot2(i,1)*hchachaR(j,1,2))*(
     .       xmut*(-2d0*xmut+th+uh) +
     .       (uh**2+uh*th-uh*(1d0+2d0*xmut+xmuneut1)+xmut+xmuneut1*
     .       xmut)) )
            endif
         enddo
      else
         charhsbot(j)=0d0
      endif

      enddo
c -------------------------------------------------------------------- c
c                       A-sbottom interference
c -------------------------------------------------------------------- c
      do j=1,2

      dha(j)   = y3-pmass(j)**2/amchar(2)**2

      charhasbot(j)=0d0

      if((amchar(1)+2d0*mt).le.amchar(2)) then

         do i=1,2
            if ((pmass(j)+amchar(1)).gt.amchar(2).and.
     .(gmsb(i)+mt).gt.amchar(2))then
            charhasbot(j)=charhasbot(j)
     .       +2d0*g2s**2/dha(j)/dsbot(i)*(Attr(j)/dsqrt(2d0))*(
     .       (alsbot1(i,2)*alsbot2(i,1)*achachaL(j,1,2)-
     .        aksbot1(i,2)*aksbot2(i,1)*achachaR(j,1,2))*(
     .       dsqrt(xmut)*(xmuneut1+xmut-uh) -
     .       dsqrt(xmut)*(-xmuneut1-xmut+th)) +
     .       (alsbot1(i,2)*alsbot2(i,1)*achachaR(j,1,2)-
     .        aksbot1(i,2)*aksbot2(i,1)*achachaL(j,1,2))*(
     .       xmchar(1)/xmchar(2)*dsqrt(xmut)*(-1d0-xmut+th) -
     .       xmchar(1)/xmchar(2)*dsqrt(xmut)*(1d0+xmut-uh) )
     .       +(alsbot2(i,1)*aksbot1(i,2)*achachaR(j,1,2)-
     .         alsbot1(i,2)*aksbot2(i,1)*achachaL(j,1,2))*(
     .       xmchar(1)/xmchar(2)*(-1d0-xmuneut1+uh+th) -
     .       2d0*xmchar(1)/xmchar(2)*xmut ) +
     .       (-alsbot2(i,1)*aksbot1(i,2)*achachaL(j,1,2)+
     .        alsbot1(i,2)*aksbot2(i,1)*achachaR(j,1,2))*(
     .       xmut*(-2d0*xmut+th+uh) -
     .       (uh**2+uh*th-uh*(1d0+2d0*xmut+xmuneut1)+xmut+xmuneut1*
     .       xmut)) )
            endif
         enddo
      else
         charhasbot(j)=0d0
      endif

      enddo
c -------------------------------------------------------------------- c
c 	                interference Z and H/h/nh,A,NA
c -------------------------------------------------------------------- c
      do i=1,3
         dhl(i)   = y3-smass(i)**2/amchar(2)**2

      charzh(i)=0d0
      if ((smass(i)+amchar(1)).gt.amchar(2).and.
     .(mz+amchar(1)).gt.amchar(2))then
      if((amchar(1)+2d0*mt).le.amchar(2)) then
         charzh(i)=
     .    -2d0*g2s**2/cw/dhl(i)/dz*vzz*(Httr(i)/dsqrt(2d0))*(
     .    xmchar(1)/xmchar(2)*dsqrt(xmut)*
     .    (hchachaL(i,1,2)*opl(1,2)+hchachaR(i,1,2)*opr(1,2))*(
     .    4d0*(th-xmut-1d0)-4d0*(uh-xmut-1d0))
     .    +dsqrt(xmut)*(hchachaL(i,1,2)*opr(1,2)
     .      +hchachaR(i,1,2)*opl(1,2))*(
     .    4d0*(xmuneut1+xmut-uh)-4d0*(xmuneut1+xmut-th)) )
      else
         charzh(i)=0d0
      endif
      endif
      enddo

      do i=1,2
         dha(i)   = y3-pmass(i)**2/amchar(2)**2

      charza(i)=0d0
      if ((pmass(i)+amchar(2)).gt.amchar(2).and.
     .(mz+amchar(1)).gt.amchar(2))then
      if((amchar(1)+2d0*mt).le.amchar(2)) then
         charza(i)=
     .    -2d0*g2s**2/cw/dha(i)/dz*azz*(Attr(i)/dsqrt(2d0))*(
     .    xmchar(1)/xmchar(2)*dsqrt(xmut)*
     .    (achachaL(i,1,2)*opl(1,2)+achachaR(i,1,2)*opr(1,2))*(
     .    4d0/xmuz*(2d0+xmuneut1*(2d0*xmut-th-uh+2d0)+xmut*
     .    (-uh-th+1d0+xmut)+xmut*(3d0*xmut-3d0*(th+uh)+5d0)
     .    +(th+uh)**2-3d0*(th+uh))+
     .    4d0*(th-xmut-1d0)+4d0*(uh-xmut-1d0))
     .    +dsqrt(xmut)*(achachaL(i,1,2)*opr(1,2)
     .      +achachaR(i,1,2)*opl(1,2))*(
     .    4d0/xmuz*(-2d0*xmuneut1**2+xmuneut1*(-6d0*xmut+
     .    3d0*(th+uh)-2d0)+xmut*(-xmut+th+uh-1d0)+xmut*
     .    (-3d0*xmut+3d0*(uh+th)-1d0)-(th+uh)**2+(th+uh))+
     .    4d0*(xmuneut1+xmut-uh)+4d0*(xmuneut1+xmut-th)) )
      else
         charza(i)=0d0
      endif
      endif
      enddo

      NS_chartop = charsbot+charztop+charzsbot

      do i=1,2
        NS_chartop=NS_chartop+charhasbot(i)+charza(i)
          do j=1,2
          NS_chartop=NS_chartop+charhatop(i,j)
        enddo
        enddo

        do i=1,3
          NS_chartop = NS_chartop+charhsbot(i)
          NS_chartop = NS_chartop+charzh(i)
        do j=1,3
          NS_chartop = NS_chartop+charhtop(i,j)
        enddo
        enddo

      end
c ==================================================================== c
c ========================= gluino up downbar ======================== c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_gluiupdb(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,i,k
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION gul(2),gur(2),gdl(2),gdr(2),gtl(2),gtr(2),
     . gbl(2),gbr(2)
      DOUBLE PRECISION dsdow(2),dsup(2)
      DOUBLE PRECISION alup(2,2),aldo(2,2)
      DOUBLE PRECISION bldo(2,2),blup(2,2)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION gs21,gs22,xmuneut1,x1,x2,x3,y1,y2,y3
     .,xmusd(2),xmusu(2)
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_coup7/alup,aldo
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2

c --- the running strong coupling ---

      gs21 = g3s
      gs22 = gs21

      xmuneut1 = mgl**2/amchar(ni)**2

      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3

      do i=1,2
         bldo(1,i) = 0d0
         bldo(2,i) = 0d0
         blup(1,i) = 0d0
         blup(2,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c                           sdown/sup exchange
c -------------------------------------------------------------------- c
      xmusd(1) = asdown1/amchar(ni)**2
      xmusd(2) = asdown2/amchar(ni)**2
      xmusu(1) = asup1/amchar(ni)**2
      xmusu(2) = asup2/amchar(ni)**2

      dsdow(1) = 1-x1-xmusd(1)
      dsdow(2) = 1-x1-xmusd(2)
      dsup(1)  = 1-x2-xmusu(1)
      dsup(2)  = 1-x2-xmusu(2)

      NS_gluiupdb=0d0

      if (dabs(mgl).le.amchar(ni)) then
         do i=1,2
            do k=1,2
               if (xmusd(k).gt.1d0.and.xmusd(i).gt.1d0)then
               NS_gluiupdb=NS_gluiupdb
     .            +g2s*gs21/dsdow(i)/dsdow(k)*x1*y1
     .            *(aldo(i,ni)*aldo(k,ni)+bldo(i,ni)*bldo(k,ni))
     .            *(gdr(i)*gdr(k)+gdl(i)*gdl(k))
               endif
               if (xmusu(k).gt.1d0.and.xmusu(i).gt.1d0)then
                  NS_gluiupdb=NS_gluiupdb
     .            +g2s*gs22/dsup(i)/dsup(k)*x2*y2
     .            *(alup(i,ni)*alup(k,ni)+blup(i,ni)*blup(k,ni))
     .            *(gur(i)*gur(k)+gul(i)*gul(k))
               endif
               if (xmusu(k).gt.1d0.and.xmusd(i).gt.1d0)then
                  NS_gluiupdb=NS_gluiupdb
     .            +g2s*dsqrt(gs21)*dsqrt(gs22)/dsdow(i)/dsup(k)*(
     .            (gul(k)*gdr(i)*aldo(i,ni)*blup(k,ni)+
     .             gdl(i)*gur(k)*alup(k,ni)*bldo(i,ni))*
     .            (-x1*y1-x2*y2+x3*y3) +
     .            2d0*(gul(k)*gdl(i)*alup(k,ni)*aldo(i,ni)+
     .                  gur(k)*gdr(i)*blup(k,ni)*bldo(i,ni))*
     .            mgl/xmchar(ni)*y3 )
                endif
            enddo
         enddo
      else
         NS_gluiupdb=0d0
      endif

      end
c ==================================================================== c
c ======================== gluino top bottombar ====================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_gluitopbb(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,i,k
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION gul(2),gur(2),gdl(2),gdr(2),gtl(2),gtr(2),
     . gbl(2),gbr(2)
      DOUBLE PRECISION dsbot(2),dstop(2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2),alstor(2,2),akstor(2,2)
      DOUBLE PRECISION  MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION alup(2,2),aldo(2,2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION gs21,gs22,xmuneut1,xmut,xmub,xmusb(2),
     .xmust(2),X1,X2,uh,th,db1,db2,db3,db4,ab1,ab2,ab3,ab4,
     .dt1,dt2,dt3,dt4,at1,at2,at3,at4,gmsb(2)
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_coup7/alup,aldo
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor

c --- the running strong coupling ---
      gmsb(1) = asb1
      gmsb(2) = asb2
      gs21 = g3s
      gs22 = gs21
      xmuneut1 = mgl**2/amchar(ni)**2
      xmut     = mt**2/amchar(ni)**2
      xmub     = mbp**2/amchar(ni)**2
      xmusb(1)   = asb1**2/amchar(ni)**2
      xmusb(2)   = asb2**2/amchar(ni)**2
      xmust(1)   = ast1**2/amchar(ni)**2
      xmust(2)   = ast2**2/amchar(ni)**2
c -------------------------------------------------------------------- c
c                         sbottom/stop exchange
c -------------------------------------------------------------------- c
      dsbot(1)  = 1d0-x1-xmusb(1)+xmut
      dsbot(2)  = 1d0-x1-xmusb(2)+xmut
      dstop(1)  = 1d0-x2-xmust(1)+xmub
      dstop(2)  = 1d0-x2-xmust(2)+xmub

      NS_gluitopbb=0d0

      uh = 1d0-x1+xmut
      th = 1d0-x2+xmub

      if ((dabs(mgl)+mt+mbp).le.amchar(ni)) then
         do i=1,2
            do k=1,2
              db1 = alsbot(k,ni)*aksbot(i,ni)+alsbot(i,ni)*aksbot(k,ni)
              db2 = alsbot(k,ni)*alsbot(i,ni)+aksbot(k,ni)*aksbot(i,ni)
              db3 = gbl(k)*gbr(i)+gbl(i)*gbr(k)
              db4 = gbl(i)*gbl(k)+gbr(i)*gbr(k)

               ab1 = db1*db3
               ab2 = db2*db3
               ab3 = db1*db4
               ab4 = db2*db4

              dt1 = alstor(k,ni)*akstor(i,ni)+alstor(i,ni)*akstor(k,ni)
              dt2 = alstor(k,ni)*alstor(i,ni)+akstor(k,ni)*akstor(i,ni)
              dt3 = gtl(k)*gtr(i)+gtl(i)*gtr(k)
              dt4 = gtl(i)*gtl(k)+gtr(i)*gtr(k)

               at1 = dt1*dt3
               at2 = dt2*dt3
               at3 = dt1*dt4
               at4 = dt2*dt4
        if ((gmsb(k)+mt).gt.amchar(ni).and.
     .(gmsb(i)+mt).gt.amchar(ni))then
               NS_gluitopbb=NS_gluitopbb+g2s*gs21/dsbot(i)/dsbot(k)*
     .          (ab1*(-4d0)*mgl/xmchar(ni)*dsqrt(xmub*xmut)+
     .           ab2*2d0*mgl/xmchar(ni)*dsqrt(xmub)*(-xmut-1d0
     .           +uh) +
     .           ab3*2d0*dsqrt(xmut)*(-xmub-xmuneut1+uh)+
     .           ab4*(-uh**2+uh*(1d0+xmuneut1+xmub+xmut)-(xmuneut1
     .           +xmub)*(1d0+xmut)) )
         endif
         if ((xmust(i)).gt.1d0.and.(xmust(k)).gt.1d0)then
            NS_gluitopbb=NS_gluitopbb
     .           +g2s*gs22/dstop(i)/dstop(k)*
     .          (at1*(-4d0)*mgl/xmchar(ni)*dsqrt(xmub*xmut)+
     .           at2*2d0*mgl/xmchar(ni)*dsqrt(xmut)*(-xmub-1d0
     .           +th) +
     .           at3*2d0*dsqrt(xmub)*(-xmut-xmuneut1+th)+
     .           at4*(-th**2+th*(1d0+xmuneut1+xmub+xmut)-(xmuneut1
     .           +xmut)*(1d0+xmub)) )
         endif
         if ((gmsb(i)+mt).gt.amchar(ni).and.(xmust(k)).gt.1d0)then
            NS_gluitopbb=NS_gluitopbb
     .           +(-2d0)*g2s*dsqrt(gs21)*dsqrt(gs22)
     .           /dstop(k)/dsbot(i)*(
     .           (alstor(k,ni)*aksbot(i,ni)*gtl(k)*gbr(i)+
     .            alsbot(i,ni)*akstor(k,ni)*gbl(i)*gtr(k))*
     .            dsqrt(xmub*xmut)*(uh+th-xmub-xmut)+
     .           (alstor(k,ni)*alsbot(i,ni)*gtl(k)*gbr(i)+
     .            aksbot(i,ni)*akstor(k,ni)*gbl(i)*gtr(k))*dsqrt(xmub)*
     .           (th-xmuneut1-xmut) +
     .           (akstor(k,ni)*alsbot(i,ni)*gtl(k)*gbl(i)+
     .            aksbot(i,ni)*alstor(k,ni)*gbr(i)*gtr(k))*
     .           mgl/xmchar(ni)*dsqrt(xmub)*(uh-xmut-1d0) +
     .           (akstor(k,ni)*aksbot(i,ni)*gtl(k)*gbl(i)+
     .            alsbot(i,ni)*alstor(k,ni)*gbr(i)*gtr(k))*(-2d0)*
     .           mgl/xmchar(ni)*dsqrt(xmub*xmut) +
     .           (akstor(k,ni)*aksbot(i,ni)*gtl(k)*gbr(i)+
     .            alsbot(i,ni)*alstor(k,ni)*gbl(i)*gtr(k))*dsqrt(xmut)*
     .           (uh-xmuneut1-xmub) +
     .           (akstor(k,ni)*alsbot(i,ni)*gtl(k)*gbr(i)+
     .            aksbot(i,ni)*alstor(k,ni)*gbl(i)*gtr(k))*
     .           (uh*th-xmuneut1-xmut*xmub) +
     .           (alstor(k,ni)*alsbot(i,ni)*gtl(k)*gbl(i)+
     .            aksbot(i,ni)*akstor(k,ni)*gbr(i)*gtr(k))*
     .           mgl/xmchar(ni)*(uh+th-xmuneut1-1d0) +
     .           (alstor(k,ni)*aksbot(i,ni)*gtl(k)*gbl(i)+
     .            alsbot(i,ni)*akstor(k,ni)*gbr(i)*gtr(k))*
     .           dsqrt(xmut*xmuneut1)*(th-xmub-1d0) )
            endif
            enddo
         enddo
      else
         NS_gluitopbb=0d0
      endif


      end
