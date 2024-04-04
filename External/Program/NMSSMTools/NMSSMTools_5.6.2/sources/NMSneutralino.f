      SUBROUTINE NS_NEUTRALINO

************************************************************************
*
*     This subroutine computes the neutralino decays
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

      INTEGER I,J,ni,nj,jsign,k,GRFLAG

      DOUBLE PRECISION amuv,lamv,amuvdiv
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION ql(5,2),qr(5,2),or(5,2),ol(5,2)
      DOUBLE PRECISION opl(2,2),opr(2,2),onl(5,5),onr(5,5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION aup(2,5),bup(2,5),ado(2,5),bdo(2,5)
      DOUBLE PRECISION hchichi(3,5,5),achichi(2,5,5)
      DOUBLE PRECISION hchachaR(3,2,2),hchachaL(3,2,2),achachaR(2,2,2),
     .         achachaL(2,2,2)
      DOUBLE PRECISION abot(2,5),bbot(2,5),atopr(2,5),btopr(2,5)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION Q2
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION neutwchar(5,2),neutzneut(5,5),neuthcchar(5,2),
     .         neutHneut(5,5,3),neutAneut(5,5,2), neutsnel(5),
     .         neutsell(5),neutselr(5),neutsnmu(5),neutsmu1(5),
     .         neutsmu2(5),neutstau1(5),neutstau2(5),neutsn1(5),
     .         neutsn2(5),neutst1(5),neutst2(5),neutsupl(5),
     .         neutsupr(5),neutsdownl(5),neutsdownr(5),neutsb1(5),
     .         neutsb2(5)
      DOUBLE PRECISION qcdneutst1(5),qcdneutst2(5),qcdneutsb1(5),
     .          qcdneutsb2(5)
      DOUBLE PRECISION qcdneutsupl(5),qcdneutsupr(5),qcdneutsdownl(5),
     .         qcdneutsdownr(5)
      DOUBLE PRECISION xneutel(5,5),xneutmu(5,5),xneuttau(5,5),
     .         xneutnue(5,5),xneutnumu(5,5),xneutnutau(5,5),
     .         xneutup(5,5),xneutdow(5,5),xneutst(5,5),xneutch(5,5),
     .         xneutbot(5,5),xneuttop(5,5),xgluinoup(5),
     .         xgluinodo(5),xgluinoch(5),xgluinost(5),xgluinobot(5),
     .         xgluinotop(5),xchelne(5,2),xchmunmu(5,2),
     .         xchtauntau(5,2),xchubdow(5,2),xchcbs(5,2),xchtbb(5,2)
      DOUBLE PRECISION neuttot2(5),neuttot(5),neuttotmulti(5),
     .         neuttotrad(5)
      DOUBLE PRECISION neuttot2lo(5),neuttot2nlo(5)
      DOUBLE PRECISION brneutst1(5),brneutst2(5),brneutsb1(5),
     .         brneutsb2(5),
     .         brneutsupl(5),brneutsupr(5),brneutsdownl(5),
     .         brneutsdownr(5),brneutsnel(5),brneutsn1(5),
     .         brneutsn2(5),brneutsell(5),brneutselr(5),
     .         brneutsnmu(5),brneutsmu1(5),brneutsmu2(5),
     .         brneutstau1(5),brneutstau2(5),brneutwchar(5,2),
     .         brneuthcchar(5,2),brneutzneut(5,5),
     .         brneutHneut(5,5,3),brneutAneut(5,5,2)
      DOUBLE PRECISION brneutup(5,5),brneutdow(5,5),brneutch(5,5),
     .         brneutst(5,5),brneutbot(5,5),brneuttop(5,5),
     .         brneutel(5,5),brneutmu(5,5),brneuttau(5,5),
     .         brneutnue(5,5),brneutnumu(5,5),brneutnutau(5,5),
     .         brchubd(5,2),brchcbs(5,2),brchtbb(5,2),brchelne(5,2),
     .         brchmunmu(5,2),brchtauntau(5,2),brglup(5),brgldo(5),
     .         brglch(5),brglst(5),brgltop(5),brglbot(5),
     .         brnraddec(5,5)
      DOUBLE PRECISION nraddec(5,5)
      DOUBLE PRECISION NS_lamb,NS_realicorr
      DOUBLE PRECISION NS_glbneut,NS_grbneut,NS_gltneut,NS_grtneut
      DOUBLE PRECISION NS_ftotqcd,ninjphoton
      DOUBLE PRECISION multilim,flagmulti,flagqcd,flagloop
      DOUBLE PRECISION neutGAMgra(5),neutZgra(5),neutHgra(5,3),
     .         neutAgra(5,2)
      DOUBLE PRECISION brcharWgra(2),brcharHCgra(2),brneutGAMgra(5),
     .         brneutZgra(5),brneutHgra(5,3),brneutAgra(5,2),
     .         brgluiGLUgra,brselEgra,brserEgra,brsmu1MUgra,
     .         brsmu2MUgra,brstau1TAUgra,brstau2TAUgra,brsneNEgra,
     .         brsnmNMgra,brsntNTgra,brsulUgra,brsurUgra,brsdlDgra,
     .         brsdrDgra,brst1Tgra,brst2Tgra,brsb1Bgra,brsb2Bgra
      DOUBLE PRECISION KNG(5),KNZ(5),KNH(5,3),KNA(5,2),KCW(2),KCH(2)
      DOUBLE PRECISION M32,CGR,MPL

      COMMON/RENSCALE/Q2
      COMMON/NEUTRALINO_2GAMMA/neutwchar,neutzneut,neuthcchar,
     .         neutHneut,neutAneut,neutsnel,neutsell,neutselr,
     .         neutsnmu,neutsmu1,neutsmu2,neutstau1,
     .         neutstau2,neutsn1,neutsn2,neutst1,
     .         neutst2,neutsupl,neutsupr,neutsdownl,
     .         neutsdownr,neutsb1,neutsb2
      COMMON/NEUTRALINO_RAD/nraddec
      COMMON/NEUTRALINO_3GAMMA/xneutel,xneutmu,xneuttau,
     .         xneutnue,xneutnumu,xneutnutau,
     .         xneutup,xneutdow,xneutst,xneutch,
     .         xneutbot,xneuttop,xgluinoup,
     .         xgluinodo,xgluinoch,xgluinost,xgluinobot,
     .         xgluinotop,xchelne,xchmunmu,
     .         xchtauntau,xchubdow,xchcbs,xchtbb
      COMMON/NEUTRALINO_WIDTH/neuttot2,neuttot,neuttotmulti,
     .         neuttotrad
      COMMON/NEUTRALINO_BR_2BD/brneutst1,brneutst2,brneutsb1,brneutsb2,
     .         brneutsupl,brneutsupr,brneutsdownl,brneutsdownr,
     .         brneutsnel,brneutsn1,brneutsn2,brneutsell,brneutselr,
     .         brneutsnmu,brneutsmu1,brneutsmu2,
     .         brneutstau1,brneutstau2,brneutwchar,brneuthcchar,
     .         brneutzneut,brneutHneut,brneutAneut,brnraddec
      COMMON/NEUTRALINO_BR_3BD/brneutup,brneutdow,brneutch,brneutst,
     .         brneutbot,brneuttop,brneutel,brneutmu,brneuttau,
     .         brneutnue,brneutnumu,brneutnutau,brchubd,brchcbs,
     .         brchtbb,brchelne,brchmunmu,brchtauntau,brglup,brgldo,
     .         brglch,brglst,brgltop,brglbot
      COMMON/GRAVITINO/brcharWgra,brcharHCgra,brneutGAMgra,
     .         brneutZgra,brneutHgra,brneutAgra,
     .         brgluiGLUgra,brselEgra,brserEgra,brsmu1MUgra,
     .         brsmu2MUgra,brstau1TAUgra,brstau2TAUgra,brsneNEgra,
     .         brsnmNMgra,brsntNTgra,brsulUgra,brsurUgra,brsdlDgra,
     .         brsdrDgra,brst1Tgra,brst2Tgra,brsb1Bgra,brsb2Bgra
      COMMON/GRAVCOUP/KNG,KNZ,KNH,KNA,KCW,KCH
      COMMON/M32/M32,CGR,MPL,GRFLAG
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/NS_coup4/opl,opr,onl,onr
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup10/aup,bup,ado,bdo
      COMMON/NS_coupNMSSM/hchichi,achichi,hchachaR,hchachaL,achachaR,
     .         achachaL
      COMMON/NS_loopdecij/ni,nj
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_neutsbotbot/abot,bbot
      COMMON/NS_neutstoptop/atopr,btopr
      COMMON/NS_multilim/multilim
      COMMON/NS_FLAGS/flagmulti,flagqcd,flagloop

      EXTERNAL NS_lamb
      EXTERNAL NS_realicorr,NS_grtneut,NS_gltneut,NS_glbneut,NS_grbneut
      EXTERNAL NS_ftotqcd

*  Initialization

      do I=1,5
         neuttot(i)     = 0d0
         neuttot2(i)    = 0d0
         neuttot2lo(i)  = 0d0
         neuttot2nlo(i) = 0d0
         neuttotmulti(i)= 0d0
         neuttotrad(i)  = 0d0

         neutst1(i)      = 0d0
         neutst2(i)      = 0d0
         neutsb1(i)      = 0d0
         neutsb2(i)      = 0d0
         neutsupl(i)     = 0d0
         neutsupr(i)     = 0d0
         neutsdownl(i)   = 0d0
         neutsdownr(i)   = 0d0
         neutsnel(i)     = 0d0
         neutsnmu(i)     = 0d0
         neutsn1(i)      = 0d0
         neutsn2(i)      = 0d0
         neutsell(i)     = 0d0
         neutselr(i)     = 0d0
         neutsmu1(i)    = 0d0
         neutsmu2(i)    = 0d0
         neutstau1(i)    = 0d0
         neutstau2(i)    = 0d0
         neuthcchar(i,1) = 0d0
         neuthcchar(i,2) = 0d0
         neutwchar(i,1)  = 0d0
         neutwchar(i,2)  = 0d0
         neutGAMgra(i)   = 0d0
         neutZgra(i)     = 0d0
         do j=1,3
           neutHgra(i,j) = 0d0 
         enddo
         do j=1,2
           neutAgra(i,j) = 0d0 
         enddo

         do j=1,5
            nraddec(i,j) = 0d0

            neutzneut(i,j)  = 0d0
            do k=1,3
               neutHneut(i,j,k) = 0d0
            enddo
            do k=1,2
               neutAneut(i,j,k) = 0d0
            enddo
            xneutel(i,j)    = 0d0
            xneutmu(i,j)    = 0d0
            xneuttau(i,j)   = 0d0
            xneutnue(i,j)   = 0d0
            xneutnumu(i,j)  = 0d0
            xneutnutau(i,j) = 0d0
            xneutup(i,j)    = 0d0
            xneutdow(i,j)   = 0d0
            xneutst(i,j)    = 0d0
            xneutch(i,j)    = 0d0
            xneutbot(i,j)   = 0d0
            xneuttop(i,j)   = 0d0
         enddo

         do j = 1,2
            xchelne(i,j)    = 0d0
            xchmunmu(i,j)   = 0d0
            xchtauntau(i,j) = 0d0
            xchubdow(i,j)   = 0d0
            xchcbs(i,j)     = 0d0
            xchtbb(i,j)     = 0d0
         enddo

         qcdneutst1(i)    = 0d0
         qcdneutst2(i)    = 0d0
         qcdneutsb1(i)    = 0d0
         qcdneutsb2(i)    = 0d0
         qcdneutsupl(i)   = 0d0
         qcdneutsupr(i)   = 0d0
         qcdneutsdownl(i) = 0d0
         qcdneutsdownr(i) = 0d0
         xgluinoup(i)     = 0d0
         xgluinodo(i)     = 0d0
         xgluinoch(i)     = 0d0
         xgluinost(i)     = 0d0
         xgluinotop(i)    = 0d0
         xgluinobot(i)    = 0d0
      enddo

      ninjphoton       = 0d0
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> chi+_1/chi+_2 + W
      do i=1,5
         do j=1,2
            if(amneut(i).ge.(mw+amchar(j))) then
               neutwchar(i,j) = g2s/2d0*(
     .              -12d0*xmneut(i)*xmchar(j)*ol(i,j)*or(i,j)
     .              +(ol(i,j)**2+or(i,j)**2)*
     .              ((amchar(j)**2+amneut(i)**2-mw**2)
     .              +(amneut(i)**2+mw**2-amchar(j)**2)
     .              *(amneut(i)**2-amchar(j)**2-mw**2)/mw**2))
     .              *NS_lamb(amchar(j)/amneut(i),mw/amneut(i))
     .              /(16*pi*amneut(i))
            else
               neutwchar(i,j)=0d0
            endif
         enddo
      enddo
C -------------------------------------------------------------------- c
c  chi0_2/chi0_3/chi0_4/chi0_5 --> chi0_1/chi0_2/chi0_3/chi0_4 + Z

      do i=1,5
         do j=1,5
            if(j.ge.i) then
               neutzneut(i,j) = 0d0
            else
               if(amneut(i).ge.(amneut(j)+mz)) then
                  neutzneut(i,j)=g2s/2d0*(
     .                 -12*xmneut(i)*xmneut(j)*onl(i,j)*onr(i,j)
     .                 +(onl(i,j)**2+onr(i,j)**2)*
     .                 ((amneut(i)**2+amneut(j)**2-mz**2)
     .                 +(amneut(i)**2+mz**2-amneut(j)**2)
     .                 *(amneut(i)**2-amneut(j)**2-mz**2)/mz**2))
     .                 *NS_lamb(amneut(j)/amneut(i),mz/amneut(i))
     .                 /(16*pi*amneut(i))
               else
                  neutzneut(i,j)=0d0
               endif
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> chi+_1/chi+_2 + H-
      do i=1,5
         do j=1,2
            if(amneut(i).ge.(amchar(j)+cmass)) then
               neuthcchar(i,j) = g2s/2d0*
     .              ((ql(i,j)**2+qr(i,j)**2)*
     .              (amchar(j)**2+amneut(i)**2-cmass**2)
     .              +4d0*ql(i,j)*qr(i,j)*xmchar(j)*xmneut(i))
     .              *NS_lamb(amchar(j)/amneut(i),cmass/amneut(i))
     .              /(16*pi*amneut(i))
            else
               neuthcchar(i,j) =0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c  chi0_2/chi0_3/chi0_4/chi0_5 --> chi0_1/chi0_2/chi0_3/chi0_4 + H(k)
c                                         neutHneut(i,j,3) [NMSSM!!!]
      do i=1,5
         do j=1,5
            do k=1,3
               if(j.ge.i) then
                  neutHneut(i,j,k) = 0d0
               else
                  if(amneut(i).ge.(amneut(j)+SMASS(k))) then
                     neutHneut(i,j,k)=g2s/2d0*
     .                    4d0*(2d0*hchichi(k,i,j)**2*
     .                    (amneut(i)**2+amneut(j)**2-SMASS(k)**2)
     .                    +4d0*hchichi(k,i,j)**2*xmneut(i)*xmneut(j))
     .               *NS_lamb(amneut(j)/amneut(i),SMASS(k)/amneut(i))
     .                    /(16*pi*amneut(i))
                  else
                     neutHneut(i,j,k)=0d0
                  endif
               endif
            enddo
         enddo
      enddo
c -------------------------------------------------------------------- c
c  chi0_2/chi0_3/chi0_4/chi0_5 --> chi0_1/chi0_2/chi0_3/chi0_4 + A(k)
c                                        neutAneut(i,j,2) [NMSSM!!!]
      do i=1,5
         do j=1,5
            do k=1,2
               if(j.ge.i) then
                  neutAneut(i,j,k) = 0d0
               else
                  if(amneut(i).ge.(amneut(j)+PMASS(k))) then
                     neutAneut(i,j,k)=g2s/2d0*
     .                    4d0*(2d0*achichi(k,i,j)**2*
     .                    (amneut(i)**2+amneut(j)**2-PMASS(k)**2)
     .                    -4d0*achichi(k,i,j)**2*xmneut(i)*xmneut(j))
     .                    *NS_lamb(amneut(j)/amneut(i),
     .                    PMASS(k)/amneut(i))/(16*pi*amneut(i))
                  else
                     neutAneut(i,j,k)=0d0
                  endif
               endif
            enddo
         enddo
      enddo
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> sneutrino_el neutrino_el
      do i=1,5
         if(amneut(i).ge.asne1) then
            neutsnel(i)=g2s/2d0*ane(1,i)**2*
     .           (amneut(i)**2-asne1**2)*NS_lamb(0d0,asne1/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutsnel(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> selectron_l electron
      do i=1,5
         if(amneut(i).ge.ase1) then
            neutsell(i)=g2s/2d0*(ae(1,i)**2+be(1,i)**2)*
     .           (amneut(i)**2-ase1**2)*NS_lamb(0d0,ase1/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutsell(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> selectron_r electron
      do i=1,5
         if(amneut(i).ge.ase2) then
            neutselr(i)=g2s/2d0*(ae(2,i)**2+be(2,i)**2)*
     .           (amneut(i)**2-ase2**2)*NS_lamb(0d0,ase2/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutselr(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> sneutrino_mu neutrino_mu
      do i=1,5
         if(amneut(i).ge.asnmu1) then
            neutsnmu(i)=g2s/2d0*anmu(1,i)**2*
     .           (amneut(i)**2-asnmu1**2)*NS_lamb(0d0,asnmu1/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutsnmu(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> smuon_1 muon
      do i=1,5
         if(amneut(i).ge.asmu1) then
            neutsmu1(i)=g2s/2d0*(amu(1,i)**2+bmu(1,i)**2)*
     .           (amneut(i)**2-asmu1**2)*NS_lamb(0d0,asmu1/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutsmu1(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> smuon_2 muon
      do i=1,5
         if(amneut(i).ge.asmu2) then
            neutsmu2(i)=g2s/2d0*(amu(2,i)**2+bmu(2,i)**2)*
     .           (amneut(i)**2-asmu2**2)*NS_lamb(0d0,asmu2/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutsmu2(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> stau1 tau

      do i=1,5
         if(amneut(i).ge.(astau1+amtau)) then
            neutstau1(i)=g2s/2d0*((atau(1,i)**2+btau(1,i)**2)*
     .           (amneut(i)**2-astau1**2+amtau**2)
     .           +4d0*amtau*xmneut(i)*atau(1,i)*btau(1,i)
     .           )*NS_lamb(amtau/amneut(i),astau1/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutstau1(i)=0d0
         endif
      enddo
C -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> stau2 tau
      do i=1,5
         if(amneut(i).ge.(astau2+amtau)) then
            neutstau2(i)=g2s/2d0*((atau(2,i)**2+btau(2,i)**2)*
     .           (amneut(i)**2-astau2**2+amtau**2)
     .           +4d0*amtau*xmneut(i)*atau(2,i)*btau(2,i)
     .           )*NS_lamb(amtau/amneut(i),astau2/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutstau2(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> sneutrino_tau1 neutrino_tau
      do i=1,5
         if(amneut(i).ge.asntau1) then
            neutsn1(i)=g2s/2d0*(antau(1,i)**2+bntau(1,i)**2)*
     .           (amneut(i)**2-asntau1**2)*
     .           NS_lamb(0d0,asntau1/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutsn1(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> sneutrino_tau2 neutrino_tau
      do i=1,5
         if(amneut(i).ge.asntau2) then
            neutsn2(i)=g2s/2d0*(antau(2,i)**2+bntau(2,i)**2)*
     .           (amneut(i)**2-asntau2**2)*
     .           NS_lamb(0d0,asntau2/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutsn2(i)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> supl + up
      do i=1,5
         if(amneut(i).ge.asup1) then
            neutsupl(i)=3d0*g2s/2d0*(aup(1,i)**2+bup(1,i)**2)
     .           *(amneut(i)**2-asup1**2)
     .           *NS_lamb(0d0,asup1/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutsupl(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,5
         if(amneut(i).ge.asup1) then
            qcdneutsupl(i)=4d0/3d0*g3s/(4d0*pi)/pi*
     .        neutsupl(i)*
     .        NS_ftotqcd(asup1**2/amneut(i)**2,mgluino**2/amneut(i)**2)
         else
            qcdneutsupl(i)=0d0
         endif
      enddo
      endif
C -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> supr + up

      do i=1,5
         if(amneut(i).ge.asup2) then
            neutsupr(i)=3d0*g2s/2d0*(aup(2,i)**2+bup(2,i)**2)
     .           *(amneut(i)**2-asup2**2)
     .           *NS_lamb(0d0,asup2/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutsupr(i)=0d0
         endif
      enddo

c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,5
         if(amneut(i).ge.asup2) then
            qcdneutsupr(i)=4d0/3d0*g3s/(4d0*pi)/pi*
     .        neutsupr(i)*
     .        NS_ftotqcd(asup2**2/amneut(i)**2,mgluino**2/amneut(i)**2)
         else
            qcdneutsupr(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> sdownl + down
      do i=1,5
         if(amneut(i).ge.asdown1) then
            neutsdownl(i)=3d0*g2s/2d0*(ado(1,i)**2+bdo(1,i)**2)
     .           *(amneut(i)**2-asdown1**2)
     .           *NS_lamb(0d0,asdown1/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutsdownl(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,5
         if(amneut(i).ge.asdown1) then
            qcdneutsdownl(i)=4d0/3d0*g3s/(4d0*pi)/pi*
     .      neutsdownl(i)*
     .      NS_ftotqcd(asdown1**2/amneut(i)**2,mgluino**2/amneut(i)**2)
         else
            qcdneutsdownl(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> sdownr + down
      do i=1,5
         if(amneut(i).ge.asdown2) then
            neutsdownr(i)=3d0*g2s/2d0*(ado(2,i)**2+bdo(2,i)**2)
     .           *(amneut(i)**2-asdown2**2)
     .           *NS_lamb(0d0,asdown2/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutsdownr(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,5
         if(amneut(i).ge.asdown2) then
            qcdneutsdownr(i)=4d0/3d0*g3s/(4d0*pi)/pi*
     .      neutsdownr(i)*
     .      NS_ftotqcd(asdown2**2/amneut(i)**2,mgluino**2/amneut(i)**2)
         else
            qcdneutsdownr(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> stop1 + top
      do i=1,5
         if(amneut(i).ge.(ast1+amt)) then
            neutst1(i)=3d0*g2s/2d0*((atopr(1,i)**2+btopr(1,i)**2)
     .           *(amneut(i)**2-ast1**2+amt**2)
     .           +4d0*amt*xmneut(i)*atopr(1,i)*btopr(1,i)
     .           )*NS_lamb(amt/amneut(i),ast1/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutst1(i)=0d0
         endif
      enddo

c -- QCD corrections --
      if(flagqcd.eq.1d0) then
      amuvdiv = amuv/dsqrt(q2)

      do i=1,5
         if(amneut(i).ge.(ast1+amt)) then
            if(xmneut(i).le.0d0) then
               jsign = 1
            else
               jsign = 0
            endif

            qcdneutst1(i)= -g2s/24d0/pi**2/amneut(i)*g3s/(4d0*pi)*
     .           3d0/2d0*
     .           ((btopr(1,i)*NS_gltneut(1,i,amuv,amuvdiv,lamv)+
     .             atopr(1,i)*NS_grtneut(1,i,amuv,amuvdiv,lamv))*
     .            (-ast1**2+amt**2+amneut(i)**2)
     .           +2d0*(btopr(1,i)*NS_grtneut(1,i,amuv,amuvdiv,lamv)
     .                 +atopr(1,i)*NS_gltneut(1,i,amuv,amuvdiv,lamv))*
     .           amt*xmneut(i))*NS_lamb(amt/amneut(i),ast1/amneut(i))
     .           +g2s/(6d0*pi**2*amneut(i))*g3s/(4d0*pi)*
     .           3d0/2d0*(-1d0)*
     .           NS_realicorr(amt,amneut(i),ast1,lamv,1,jsign,1,i,1)
         else
            qcdneutst1(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> stop2 + top
      do i=1,5
         if(amneut(i).ge.(ast2+amt)) then
            neutst2(i)=3d0*g2s/2d0*((atopr(2,i)**2+btopr(2,i)**2)*
     .           (amneut(i)**2-ast2**2+amt**2)
     .           +4d0*amt*xmneut(i)*atopr(2,i)*btopr(2,i)
     .           )*NS_lamb(amt/amneut(i),ast2/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutst2(i)=0d0
         endif
      enddo

c -- QCD corrections --
      if(flagqcd.eq.1d0) then
      do i=1,5
         if(amneut(i).ge.(ast2+amt)) then
            if(xmneut(i).le.0d0) then
               jsign = 1
            else
               jsign = 0
            endif

            qcdneutst2(i)= -g2s/24d0/pi**2/amneut(i)*g3s/(4d0*pi)*
     .           3d0/2d0*
     .           ((btopr(2,i)*NS_gltneut(2,i,amuv,amuvdiv,lamv)+
     .             atopr(2,i)*NS_grtneut(2,i,amuv,amuvdiv,lamv))*
     .            (-ast2**2+amt**2+amneut(i)**2)
     .           +2d0*(btopr(2,i)*NS_grtneut(2,i,amuv,amuvdiv,lamv)
     .                 +atopr(2,i)*NS_gltneut(2,i,amuv,amuvdiv,lamv))*
     .           amt*xmneut(i))*NS_lamb(amt/amneut(i),ast2/amneut(i))
     .           +g2s/(6d0*pi**2*amneut(i))*g3s/(4d0*pi)*
     .           3d0/2d0*(-1d0)*
     .           NS_realicorr(amt,amneut(i),ast2,lamv,1,jsign,2,i,1)
         else
            qcdneutst2(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> sbottom1 + bottom
      do i=1,5
         if(amneut(i).ge.(asb1+amb)) then
            neutsb1(i)=3d0*g2s/2d0*((abot(1,i)**2+bbot(1,i)**2)*
     .           (amneut(i)**2-asb1**2+amb**2)
     .           +4d0*amb*xmneut(i)*abot(1,i)*bbot(1,i)
     .           )*NS_lamb(amb/amneut(i),asb1/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutsb1(i)=0d0
         endif
      enddo

c -- QCD corrections --
      if(flagqcd.eq.1d0) then
      do i=1,5
         if(amneut(i).ge.(asb1+amb)) then
            if(xmneut(i).le.0d0) then
               jsign = 1
            else
               jsign = 0
            endif

            qcdneutsb1(i)= -g2s/24d0/pi**2/amneut(i)*g3s/(4d0*pi)*
     .           3d0/2d0*
     .           ((bbot(1,i)*NS_glbneut(1,i,amuv,amuvdiv,lamv)+
     .             abot(1,i)*NS_grbneut(1,i,amuv,amuvdiv,lamv))*
     .            (-asb1**2+amb**2+amneut(i)**2)
     .           +2d0*(bbot(1,i)*NS_grbneut(1,i,amuv,amuvdiv,lamv)
     .                 +abot(1,i)*NS_glbneut(1,i,amuv,amuvdiv,lamv))*
     .           amb*xmneut(i))*NS_lamb(amb/amneut(i),asb1/amneut(i))
     .           +g2s/(6d0*pi**2*amneut(i))*g3s/(4d0*pi)*
     .           3d0/2d0*(-1d0)*
     .           NS_realicorr(amb,amneut(i),asb1,lamv,1,jsign,1,i,2)
         else
            qcdneutsb1(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> sbottom2 + bottom
      do i=1,5
         if(amneut(i).ge.(asb2+amb)) then
            neutsb2(i)=3d0*g2s/2d0*((abot(2,i)**2+bbot(2,i)**2)*
     .           (amneut(i)**2-asb2**2+amb**2)
     .           +4d0*amb*xmneut(i)*abot(2,i)*bbot(2,i)
     .           )*NS_lamb(amb/amneut(i),asb2/amneut(i))
     .           /(16*pi*amneut(i))
         else
            neutsb2(i)=0d0
         endif
      enddo
c -- QCD corrections --
      if(flagqcd.eq.1d0) then
      do i=1,5
         if(amneut(i).ge.(asb2+amb)) then
            if(xmneut(i).le.0d0) then
               jsign = 1
            else
               jsign = 0
            endif

            qcdneutsb2(i)= -g2s/24d0/pi**2/amneut(i)*g3s/(4d0*pi)*
     .           3d0/2d0*
     .           ((bbot(2,i)*NS_glbneut(2,i,amuv,amuvdiv,lamv)+
     .             abot(2,i)*NS_grbneut(2,i,amuv,amuvdiv,lamv))*
     .            (-asb2**2+amb**2+amneut(i)**2)
     .           +2d0*(bbot(2,i)*NS_grbneut(2,i,amuv,amuvdiv,lamv)
     .                 +abot(2,i)*NS_glbneut(2,i,amuv,amuvdiv,lamv))*
     .           amb*xmneut(i))*NS_lamb(amb/amneut(i),asb2/amneut(i))
     .           +g2s/(6d0*pi**2*amneut(i))*g3s/(4d0*pi)*
     .           3d0/2d0*(-1d0)*
     .           NS_realicorr(amb,amneut(i),asb2,lamv,1,jsign,2,i,2)
         else
            qcdneutsb2(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> gamma + gravitino

      if(GRFLAG.EQ.1)then
        do i=1,5
          if (M32.le.amneut(i)) then
            neutGAMgra(i) = KNG(i)*amneut(i)**5/(48d0*PI*MPL**2*M32**2)
          else
            neutGAMgra(i) = 0d0
          endif
        enddo
      ENDif
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> Z + gravitino

      if(GRFLAG.EQ.1)then
        do i=1,5
          if ((MZ+M32).le.amneut(i)) then
            neutZgra(i) = KNZ(i)*amneut(i)**5/(48d0*PI*MPL**2*M32**2)
     .                    *(1d0-mz**2/amneut(i)**2)**4
          else
            neutZgra(i) = 0d0
          endif
        enddo
      ENDif
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> H(j) + gravitino

      if(GRFLAG.EQ.1)then
        do i=1,5
          do j=1,3
            if ((SMASS(j)+M32).le.amneut(i)) then
              neutHgra(i,j) = KNH(i,j)*amneut(i)**5
     .                        /(48d0*PI*MPL**2*M32**2)
     .                        *(1d0-SMASS(j)**2/amneut(i)**2)**4
            else
              neutHgra(i,j) = 0d0
            endif
          enddo
        enddo
      ENDif
c -------------------------------------------------------------------- c
c  chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 --> A(j) + gravitino

      if(GRFLAG.EQ.1)then
        do i=1,5
          do j=1,2
            if ((PMASS(j)+M32).le.amneut(i)) then
              neutAgra(i,j) = KNA(i,j)*amneut(i)**5
     .                        /(48d0*PI*MPL**2*M32**2)
     .                        *(1d0-PMASS(j)**2/amneut(i)**2)**4
            else
              neutAgra(i,j) = 0d0
            endif
          enddo
        enddo
      ENDif
c -------------------------------------------------------------------- c
c                  SUM 2 BODY DECAYS
c -------------------------------------------------------------------- c
      do i=1,5
         neuttot2lo(i) = 2d0*(neutst1(i)+neutst2(i)+neutsb1(i)+
     .           neutsb2(i)+2d0*neutsupl(i)+2d0*neutsupr(i)+
     .           2d0*neutsdownl(i)+2d0*neutsdownr(i)+
     .           neutsnel(i)+neutsn1(i)+neutsn2(i)+neutsell(i)+
     .           neutselr(i)+neutstau1(i)+neutstau2(i)+
     .           neutsnmu(i)+neutsmu1(i)+neutsmu2(i)+
     .           neuthcchar(i,1)+neuthcchar(i,2)+neutwchar(i,1)+
     .           neutwchar(i,2))+neutGAMgra(i)+neutZgra(i)
         do j=1,3
            neuttot2lo(i) = neuttot2lo(i)+neutHgra(i,j)
         enddo
         do j=1,2
            neuttot2lo(i) = neuttot2lo(i)+neutAgra(i,j)
         enddo
         do j=1,5
            neuttot2lo(i) = neuttot2lo(i)+neutzneut(i,j)
            do k=1,3
               neuttot2lo(i) = neuttot2lo(i)+neutHneut(i,j,k)
            enddo
            do k=1,2
               neuttot2lo(i) = neuttot2lo(i)+neutAneut(i,j,k)
            enddo
         enddo
      enddo

         do i=1,5
            neuttot2nlo(i) = 2d0*(qcdneutst1(i)+qcdneutst2(i)+
     .           qcdneutsb1(i)+qcdneutsb2(i)+2d0*qcdneutsupl(i)+2d0*
     .           qcdneutsupr(i)+2d0*qcdneutsdownl(i)+2d0*
     .           qcdneutsdownr(i)) + neuttot2lo(i)
         enddo

         do i=1,5
            neuttot2(i) = neuttot2nlo(i)
         enddo
c--------------------------------------------------------------------- c
c ---------------- 3-body decays and 3-body total widths -------------
c--------------------------------------------------------------------- c
      if(flagmulti.eq.1d0) then

      CALL NS_xintegneut

      do i=1,5
            neuttotmulti(i) = 0d0
            do j=1,5
               neuttotmulti(i) = neuttotmulti(i)+xneutel(i,j)+
     .              xneutmu(i,j)+xneuttau(i,j)+xneutnue(i,j)+
     .              xneutnumu(i,j)+xneutnutau(i,j)+xneutup(i,j)+
     .              xneutdow(i,j)+xneutst(i,j)+xneutch(i,j)+
     .              xneutbot(i,j)+xneuttop(i,j)

            enddo
            do j=1,2
               neuttotmulti(i) = neuttotmulti(i)+2d0*(xchelne(i,j)+
     .              xchmunmu(i,j)+xchtauntau(i,j)+
     .              xchubdow(i,j)+xchcbs(i,j)+xchtbb(i,j))

            enddo
            neuttotmulti(i) = neuttotmulti(i)+xgluinoup(i)+
     .           xgluinodo(i)+xgluinoch(i)+xgluinost(i)+
     .           xgluinotop(i)+xgluinobot(i)

      enddo
      endif
c ---- Add 3-body decays if BR>multilim -------------------------------c
      do i = 1,5
         if (neuttotmulti(i).lt.multilim*neuttot2(i))then
            neuttotmulti(i)=0d0
         endif
      enddo
c--------------------------------------------------------------------- c
c -- loop decays neutralino_i -> neutralino_j + photon
c--------------------------------------------------------------------- c
      if(flagloop.eq.1d0) then
      do ni = 1,5
            do nj = 1,5
               if(nj.ge.ni) then
                  nraddec(ni,nj)=0d0
               elseif(nj.lt.ni) then
                  CALL NS_neutraddecay(ninjphoton)
                  nraddec(ni,nj)=ninjphoton
c ---- Add loop decays only if BR>multilim --------------- ------------c
                  if(nraddec(ni,nj).lt.multilim*neuttot2(ni)) then
                    nraddec(ni,nj)=0d0
                  endif
               endif
            enddo
      enddo

      do i=1,5
            neuttotrad(i) = 0d0
            do j=1,5
               neuttotrad(i) = neuttotrad(i)+nraddec(i,j)
            enddo
      enddo
      endif
c--------------------------------------------------------------------- c
c ---------------------------- total widths -------------------------- c
c--------------------------------------------------------------------- c
      do i=1,5
       neuttot(i) = neuttot2(i)+neuttotmulti(i)+neuttotrad(i)

      enddo
c--------------------------------------------------------------------- c
c ------ neutralino branching ratios for 2-body decays --------------- c
c--------------------------------------------------------------------- c
      do i=1,5
            neutst1(i)    = neutst1(i)+qcdneutst1(i)
            neutst2(i)    = neutst2(i)+qcdneutst2(i)
            neutsb1(i)    = neutsb1(i)+qcdneutsb1(i)
            neutsb2(i)    = neutsb2(i)+qcdneutsb2(i)
            neutsupl(i)   = neutsupl(i)+qcdneutsupl(i)
            neutsupl(i)   = neutsupl(i)+qcdneutsupl(i)
            neutsupr(i)   = neutsupr(i)+qcdneutsupr(i)
            neutsupr(i)   = neutsupr(i)+qcdneutsupr(i)
            neutsdownl(i) = neutsdownl(i)+qcdneutsdownl(i)
            neutsdownl(i) = neutsdownl(i)+qcdneutsdownl(i)
            neutsdownr(i) = neutsdownr(i)+qcdneutsdownr(i)
            neutsdownr(i) = neutsdownr(i)+qcdneutsdownr(i)
      enddo

      do i=1,5
       if(neuttot(i).ne.0d0)then

         brneutst1(i)    = neutst1(i)/neuttot(i)
         brneutst2(i)    = neutst2(i)/neuttot(i)
         brneutsb1(i)    = neutsb1(i)/neuttot(i)
         brneutsb2(i)    = neutsb2(i)/neuttot(i)
         brneutsupl(i)   = neutsupl(i)/neuttot(i)
         brneutsupr(i)   = neutsupr(i)/neuttot(i)
         brneutsdownl(i) = neutsdownl(i)/neuttot(i)
         brneutsdownr(i) = neutsdownr(i)/neuttot(i)
         brneutsnel(i)   = neutsnel(i)/neuttot(i)
         brneutsnmu(i)   = neutsnmu(i)/neuttot(i)
         brneutsn1(i)    = neutsn1(i)/neuttot(i)
         brneutsn2(i)    = neutsn2(i)/neuttot(i)
         brneutsell(i)   = neutsell(i)/neuttot(i)
         brneutselr(i)   = neutselr(i)/neuttot(i)
         brneutsmu1(i)   = neutsmu1(i)/neuttot(i)
         brneutsmu2(i)   = neutsmu2(i)/neuttot(i)
         brneutstau1(i)  = neutstau1(i)/neuttot(i)
         brneutstau2(i)  = neutstau2(i)/neuttot(i)
         brneutGAMgra(i) = neutGAMgra(i)/neuttot(i)
         brneutZgra(i)   = neutZgra(i)/neuttot(i)
         do j=1,3
            brneutHgra(i,j)  = neutHgra(i,j)/neuttot(i)
         enddo
         do j=1,2
            brneutAgra(i,j)  = neutAgra(i,j)/neuttot(i)
            brneutwchar(i,j)  = neutwchar(i,j)/neuttot(i)
            brneuthcchar(i,j) = neuthcchar(i,j)/neuttot(i)
         enddo
         do j=1,5
            brneutzneut(i,j)  = neutzneut(i,j)/neuttot(i)
            do k=1,3
               brneutHneut(i,j,k) = neutHneut(i,j,k)/neuttot(i)
            enddo
            do k=1,2
               brneutAneut(i,j,k) = neutAneut(i,j,k)/neuttot(i)
            enddo
         enddo

       else

         brneutst1(i)    = 0d0
         brneutst2(i)    = 0d0
         brneutsb1(i)    = 0d0
         brneutsb2(i)    = 0d0
         brneutsupl(i)   = 0d0
         brneutsupr(i)   = 0d0
         brneutsdownl(i) = 0d0
         brneutsdownr(i) = 0d0
         brneutsnel(i)   = 0d0
         brneutsnmu(i)   = 0d0
         brneutsn1(i)    = 0d0
         brneutsn2(i)    = 0d0
         brneutsell(i)   = 0d0
         brneutselr(i)   = 0d0
         brneutsmu1(i)   = 0d0
         brneutsmu2(i)   = 0d0
         brneutstau1(i)  = 0d0
         brneutstau2(i)  = 0d0
         brneutGAMgra(i) = 0d0
         brneutZgra(i)   = 0d0
         do j=1,3
            brneutHgra(i,j) = 0d0
         enddo
         do j=1,2
            brneutAgra(i,j)   = 0d0
            brneutwchar(i,j)  = 0d0
            brneuthcchar(i,j) = 0d0
         enddo
         do j=1,5
            brneutzneut(i,j) = 0d0
            do k=1,3
               brneutHneut(i,j,k) = 0d0
            enddo
            do k=1,2
               brneutAneut(i,j,k) = 0d0
            enddo
         enddo

       endif
      enddo

c--------------------------------------------------------------------- c
c ------ neutralino branching ratios for 3-body and loop decays ------ c
c--------------------------------------------------------------------- c
      do i=1,5

       if(neuttotmulti(i).ne.0d0)then

         do j=1,5
            brneutup(i,j)    = xneutup(i,j)/neuttot(i)
            brneutdow(i,j)   = xneutdow(i,j)/neuttot(i)
            brneutch(i,j)    = xneutch(i,j)/neuttot(i)
            brneutst(i,j)    = xneutst(i,j)/neuttot(i)
            brneutbot(i,j)   = xneutbot(i,j)/neuttot(i)
            brneuttop(i,j)   = xneuttop(i,j)/neuttot(i)
            brneutel(i,j)    = xneutel(i,j)/neuttot(i)
            brneutmu(i,j)    = xneutmu(i,j)/neuttot(i)
            brneuttau(i,j)   = xneuttau(i,j)/neuttot(i)
            brneutnue(i,j)   = xneutnue(i,j)/neuttot(i)
            brneutnumu(i,j)  = xneutnumu(i,j)/neuttot(i)
            brneutnutau(i,j) = xneutnutau(i,j)/neuttot(i)
         enddo
         do j=1,2
            brchubd(i,j)     = xchubdow(i,j)/neuttot(i)
            brchcbs(i,j)     = xchcbs(i,j)/neuttot(i)
            brchtbb(i,j)     = xchtbb(i,j)/neuttot(i)
            brchelne(i,j)    = xchelne(i,j)/neuttot(i)
            brchmunmu(i,j)   = xchmunmu(i,j)/neuttot(i)
            brchtauntau(i,j) = xchtauntau(i,j)/neuttot(i)
         enddo
         brglup(i)  = xgluinoup(i)/neuttot(i)
         brgldo(i)  = xgluinodo(i)/neuttot(i)
         brglch(i)  = xgluinoch(i)/neuttot(i)
         brglst(i)  = xgluinost(i)/neuttot(i)
         brgltop(i) = xgluinotop(i)/neuttot(i)
         brglbot(i) = xgluinobot(i)/neuttot(i)

       else

         do j=1,5
            brneutup(i,j)    = 0d0
            brneutdow(i,j)   = 0d0
            brneutch(i,j)    = 0d0
            brneutst(i,j)    = 0d0
            brneutbot(i,j)   = 0d0
            brneuttop(i,j)   = 0d0
            brneutel(i,j)    = 0d0
            brneutmu(i,j)    = 0d0
            brneuttau(i,j)   = 0d0
            brneutnue(i,j)   = 0d0
            brneutnumu(i,j)  = 0d0
            brneutnutau(i,j) = 0d0
         enddo
         do j=1,2
            brchubd(i,j)     = 0d0
            brchcbs(i,j)     = 0d0
            brchtbb(i,j)     = 0d0
            brchelne(i,j)    = 0d0
            brchmunmu(i,j)   = 0d0
            brchtauntau(i,j) = 0d0
         enddo
         brglup(i)  = 0d0
         brgldo(i)  = 0d0
         brglch(i)  = 0d0
         brglst(i)  = 0d0
         brgltop(i) = 0d0
         brglbot(i) = 0d0

       endif

       if(neuttot(i).ne.0d0)then

         do j=1,5
            brnraddec(i,j)   = nraddec(i,j)/neuttot(i)
         enddo

       else

         do j=1,5
            brnraddec(i,j) = 0d0
         enddo

       endif
 
      enddo

      end
c ==================================================================== c
c          Radiative decays neutralino_i -> neutralino_j photon        c
c ==================================================================== c
      SUBROUTINE NS_neutraddecay(ninjphoton)
*
      IMPLICIT NONE
      INTEGER Ni,nj
      DOUBLE PRECISION ninjphoton
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION eps(5),gefgh(2)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION NS_gcoupefgh,NS_gcoupabcd,NS_gcoupabcd0
      DOUBLE PRECISION gabcd,gabcd0,gijgamma
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_loopdecij/ni,nj
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      EXTERNAL NS_gcoupefgh,NS_gcoupabcd,NS_gcoupabcd0
*
      if(xmneut(nj).ge.0d0) then
         eps(nj) = 1d0
      else
         eps(nj) = -1d0
      endif

       if(amneut(ni).gt.amneut(nj)) then
            gefgh(1) = NS_gcoupefgh(xmneut(ni),xmneut(nj),xmchar(1),ni,
     .                              nj,1)
            gefgh(2) = NS_gcoupefgh(xmneut(ni),xmneut(nj),xmchar(2),ni,
     .                              nj,2)
         gabcd = NS_gcoupabcd(xmneut(ni),xmneut(nj),ni,nj)
         gabcd0 = NS_gcoupabcd0(xmneut(ni),xmneut(nj),ni,nj)
         gijgamma = g2s*dsqrt(g2s)*sw*amneut(ni)*eps(nj)/8d0/pi**2*
     .        (gefgh(1)+gefgh(2)+gabcd+gabcd0)
        else
            gijgamma = 0d0
        endif

      ninjphoton = gijgamma**2*(amneut(ni)**2-amneut(nj)**2)**3/
     .             8d0/pi/amneut(ni)**5

       end
c -------------------------------------------------------------------- c
      DOUBLE PRECISION FUNCTION NS_gcoupefgh(xmni,xmnj,xmchar,i,j,k)
*
      IMPLICIT NONE
      INTEGER i,j,k
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION ql(5,2),qr(5,2),or(5,2),ol(5,2)
      DOUBLE PRECISION XMNI,xmnj,amni,amnj,epsi,epsj,fl,fr,gl,gr,
     .amchar,xmchar
      COMPLEX*16 NS_iint,NS_i2int,NS_jint,NS_kint
*
      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      EXTERNAL NS_iint,NS_i2int,NS_jint,NS_kint

      if(xmni.le.0d0) then
         epsi = -1d0
      else
         epsi = 1d0
      endif
      if(xmnj.le.0d0) then
         epsj = -1d0
      else
         epsj = 1d0
      endif

      amni = dabs(xmni)
      amnj = dabs(xmnj)
      amchar = dabs(xmchar)
      fl = ol(i,k)
      fr = or(i,k)
      gl = ol(j,k)
      gr = or(j,k)

      NS_gcoupefgh = (gl*fl-gr*fr)*(epsi*amni*
     .     dreal((NS_i2int(amni,amnj,amchar,mw)
     .  -NS_jint(amni,amnj,amchar,mw)-NS_kint(amni,amnj,amchar,mw)))
     .     + epsj*amnj*dreal((NS_jint(amni,amnj,amchar,mw)
     .                       -NS_kint(amni,amnj,amchar,mw))) )
     .     + 2d0*amchar*(gl*fr-gr*fl)*
     .       dreal(NS_jint(amni,amnj,amchar,mw))

      return

      end

c -------------------------------------------------------------------- c
      DOUBLE PRECISION FUNCTION NS_gcoupabcd(xmni,xmnj,i,j)
*
      IMPLICIT NONE
      INTEGER I,J,K
      DOUBLE PRECISION abot(2,5),bbot(2,5)
      DOUBLE PRECISION atopr(2,5),btopr(2,5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION epsi,epsj,xmni,xmnj,amni,amnj,beta,sbeta,cbeta
     .,zip,zim,zjp,zjm,ef,cf,mfer,mbos,gl,gr,fl,fr,glfrgrfl,glflgrfr
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMPLEX*16 NS_iint,NS_i2int,NS_jint,NS_kint
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_neutsbotbot/abot,bbot
      COMMON/NS_neutstoptop/atopr,btopr

      EXTERNAL NS_iint,NS_i2int,NS_jint,NS_kint

      if(xmni.le.0d0) then
         epsi = -1d0
      else
         epsi = 1d0
      endif

      if(xmnj.le.0d0) then
         epsj = -1d0
      else
         epsj = 1d0
      endif

      amni = dabs(xmni)
      amnj = dabs(xmnj)

      beta  = datan(tanbq)
      sbeta = DDSIN(beta)
      cbeta = DDCOS(beta)

      zip = zz(i,2)+zz(i,1)*sw/cw
      zim = zz(i,2)-zz(i,1)*sw/cw
      zjp = zz(j,2)+zz(j,1)*sw/cw
      zjm = zz(j,2)-zz(j,1)*sw/cw
      NS_gcoupabcd  = 0d0

      do k=1,6,1
         if(k.eq.1) then
            ef   = 2d0/3d0
            cf   = 3d0
            mfer = amt
            mbos = ast1
            gl   = -dsqrt(2d0)*atopr(1,j)
            gr   = -dsqrt(2d0)*btopr(1,j)
            fl   = -dsqrt(2d0)*btopr(1,i)
            fr   = -dsqrt(2d0)*atopr(1,i)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.2) then
            ef   = 2d0/3d0
            cf   = 3d0
            mfer = amt
            mbos = ast2
            gl   = -dsqrt(2d0)*atopr(2,j)
            gr   = -dsqrt(2d0)*btopr(2,j)
            fl   = -dsqrt(2d0)*btopr(2,i)
            fr   = -dsqrt(2d0)*atopr(2,i)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.3) then
            ef   = -1d0/3d0
            cf   = 3d0
            mfer = amb
            mbos = asb1
            gl   = -dsqrt(2d0)*abot(1,j)
            gr   = -dsqrt(2d0)*bbot(1,j)
            fl   = -dsqrt(2d0)*bbot(1,i)
            fr   = -dsqrt(2d0)*abot(1,i)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.4) then
            ef   = -1d0/3d0
            cf   = 3d0
            mfer = amb
            mbos = asb2
            gl   = -dsqrt(2d0)*abot(2,j)
            gr   = -dsqrt(2d0)*bbot(2,j)
            fl   = -dsqrt(2d0)*bbot(2,i)
            fr   = -dsqrt(2d0)*abot(2,i)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.5) then
            ef   = -1d0
            cf   = 1d0
            mfer = amtau
            mbos = astau1
            gl   = -dsqrt(2d0)*atau(1,j)
            gr   = -dsqrt(2d0)*btau(1,j)
            fl   = -dsqrt(2d0)*btau(1,i)
            fr   = -dsqrt(2d0)*atau(1,i)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.6) then
            ef   = -1d0
            cf   = 1d0
            mfer = amtau
            mbos = astau2
            gl   = -dsqrt(2d0)*atau(2,j)
            gr   = -dsqrt(2d0)*btau(2,j)
            fl   = -dsqrt(2d0)*btau(2,i)
            fr   = -dsqrt(2d0)*atau(2,i)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         endif

         NS_gcoupabcd = NS_gcoupabcd -1d0/4d0*ef*cf*(
     .        glfrgrfl*(epsi*amni*dreal((NS_i2int(amni,amnj,mfer,mbos)-
     .                             NS_kint(amni,amnj,mfer,mbos)))
     .                 -epsj*amnj*dreal(NS_kint(amni,amnj,mfer,mbos)) )
     .        +mfer*glflgrfr*dreal(NS_iint(amni,amnj,mfer,mbos)) )
      enddo

      do k=1,2
         ef   = 1d0
         cf   = 1d0
         mfer = amchar(k)
         mbos = cmass

         glfrgrfl = 2d0*cbeta**2*(zz(j,4)*vv(k,1)+1d0/dsqrt(2d0)*
     .       zjp*vv(k,2))*(zz(i,4)*vv(k,1)+1d0/dsqrt(2d0)*zip*vv(k,2))
     .        -2d0*sbeta**2*(zz(j,3)*uu(k,1)-1d0/dsqrt(2d0)*zjp*
     .        uu(k,2))*(zz(i,3)*uu(k,1)-1d0/dsqrt(2d0)*zip*uu(k,2))
         glflgrfr = DDSIN(2d0*beta)*(uu(k,1)*vv(k,1)*(zz(j,4)*zz(i,3)-
     .        zz(j,3)*zz(i,4))-1d0/dsqrt(2d0)*zip*(zz(j,3)*uu(k,1)*
     .        vv(k,2)+zz(j,4)*uu(k,2)*vv(k,1))+1d0/dsqrt(2d0)*zjp*(
     .        zz(i,3)*uu(k,1)*vv(k,2)+zz(i,4)*uu(k,2)*vv(k,1)))

         NS_gcoupabcd = NS_gcoupabcd -1d0/4d0*ef*cf*(
     .        glfrgrfl*(epsi*amni*dreal((NS_i2int(amni,amnj,mfer,mbos)-
     .                             NS_kint(amni,amnj,mfer,mbos)))
     .                 -epsj*amnj*dreal(NS_kint(amni,amnj,mfer,mbos)) )
     .        +mfer*glflgrfr*dreal(NS_iint(amni,amnj,mfer,mbos)) )
      enddo

      do k=1,2
         ef   = 1d0
         cf   = 1d0
         mfer = amchar(k)
         mbos = mw

         glfrgrfl = 2d0*sbeta**2*(zz(j,4)*vv(k,1)+1d0/dsqrt(2d0)*
     .       zjp*vv(k,2))*(zz(i,4)*vv(k,1)+1d0/dsqrt(2d0)*zip*vv(k,2))
     .        -2d0*cbeta**2*(zz(j,3)*uu(k,1)-1d0/dsqrt(2d0)*zjp*
     .        uu(k,2))*(zz(i,3)*uu(k,1)-1d0/dsqrt(2d0)*zip*uu(k,2))
         glflgrfr = -DDSIN(2d0*beta)*(uu(k,1)*vv(k,1)*(zz(j,4)*zz(i,3)-
     .        zz(j,3)*zz(i,4))-1d0/dsqrt(2d0)*zip*(zz(j,3)*uu(k,1)*
     .        vv(k,2)+zz(j,4)*uu(k,2)*vv(k,1))+1d0/dsqrt(2d0)*zjp*(
     .        zz(i,3)*uu(k,1)*vv(k,2)+zz(i,4)*uu(k,2)*vv(k,1)))

         NS_gcoupabcd = NS_gcoupabcd -1d0/4d0*ef*cf*(
     .        glfrgrfl*(epsi*amni*dreal((NS_i2int(amni,amnj,mfer,mbos)-
     .                                   NS_kint(amni,amnj,mfer,mbos)))
     .                 -epsj*amnj*dreal(NS_kint(amni,amnj,mfer,mbos)) )
     .        +mfer*glflgrfr*dreal(NS_iint(amni,amnj,mfer,mbos)) )
      enddo

      return

      end
c -------------------------------------------------------------------- c
      DOUBLE PRECISION FUNCTION NS_gcoupabcd0(xmni,xmnj,i,j)
*
      IMPLICIT NONE
      INTEGER I,J,K
      DOUBLE PRECISION aup(2,5),bup(2,5),ado(2,5),bdo(2,5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      DOUBLE PRECISION epsi,epsj,xmni,xmnj,amni,amnj,
     .ef,cf,mfer,mbos,gl,gr,fl,fr,glfrgrfl,glflgrfr,
     .mu,md,me
      COMPLEX*16 NS_i2int0,NS_kint0,NS_iint
*
      COMMON/NS_pi/PI,SQR2
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup10/aup,bup,ado,bdo

      EXTERNAL NS_i2int0,NS_kint0,NS_iint

      if(xmni.le.0d0) then
         epsi = -1d0
      else
         epsi = 1d0
      endif

      if(xmnj.le.0d0) then
         epsj = -1d0
      else
         epsj = 1d0
      endif

      amni = dabs(xmni)
      amnj = dabs(xmnj)

      mu = 1d-2
      md = 1d-2
      me = 1d-2

      NS_gcoupabcd0  = 0d0

      do k=1,6,1
         if(k.eq.1) then
            ef   = 2d0/3d0
            cf   = 3d0
            mfer = mu
            mbos = asup1
            gl   = -dsqrt(2d0)*aup(1,j)
            gr   = -dsqrt(2d0)*bup(1,j)
            fl   = -dsqrt(2d0)*bup(1,i)
            fr   = -dsqrt(2d0)*aup(1,i)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.2) then
            ef   = 2d0/3d0
            cf   = 3d0
            mfer = mu
            mbos = asup2
            gl   = -dsqrt(2d0)*aup(2,j)
            gr   = -dsqrt(2d0)*bup(2,j)
            fl   = -dsqrt(2d0)*bup(2,i)
            fr   = -dsqrt(2d0)*aup(2,i)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.3) then
            ef   = -1d0/3d0
            cf   = 3d0
            mfer = md
            mbos = asdown1
            gl   = -dsqrt(2d0)*ado(1,j)
            gr   = -dsqrt(2d0)*bdo(1,j)
            fl   = -dsqrt(2d0)*bdo(1,i)
            fr   = -dsqrt(2d0)*ado(1,i)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.4) then
            ef   = -1d0/3d0
            cf   = 3d0
            mfer = md
            mbos = asdown2
            gl   = -dsqrt(2d0)*ado(2,j)
            gr   = -dsqrt(2d0)*bdo(2,j)
            fl   = -dsqrt(2d0)*bdo(2,i)
            fr   = -dsqrt(2d0)*ado(2,i)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.5) then
            ef   = -1d0
            cf   = 1d0
            mfer = me
            mbos = asmu1
            gl   = -dsqrt(2d0)*ae(1,j)
            gr   = -dsqrt(2d0)*be(1,j)
            fl   = -dsqrt(2d0)*be(1,i)
            fr   = -dsqrt(2d0)*ae(1,i)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         elseif(k.eq.6) then
            ef   = -1d0
            cf   = 1d0
            mfer = me
            mbos = asmu2
            gl   = -dsqrt(2d0)*ae(2,j)
            gr   = -dsqrt(2d0)*be(2,j)
            fl   = -dsqrt(2d0)*be(2,i)
            fr   = -dsqrt(2d0)*ae(2,i)
            glfrgrfl = gl*fr-gr*fl
            glflgrfr = gl*fl-gr*fr
         endif

         NS_gcoupabcd0 = NS_gcoupabcd0 -1d0/4d0*ef*cf*(
     .        glfrgrfl*(epsi*amni*dreal((NS_i2int0(amni,amnj,mbos)-
     .                             NS_kint0(amni,amnj,mbos)))
     .                 -epsj*amnj*dreal(NS_kint0(amni,amnj,mbos)) ) )
      enddo

      NS_gcoupabcd0 = 2d0*NS_gcoupabcd0

      return

      end
c ==================================================================== c
c                   neutralino_1/2/3/4/5 3-body decays
c ==================================================================== c
      SUBROUTINE NS_xintegneut
*
      IMPLICIT NONE
      INTEGER ni,nj,I
      INTEGER nx1t,ny1t
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION xneutel(5,5),xneutmu(5,5),xneuttau(5,5),
     .          xneutnue(5,5),xneutnumu(5,5),xneutnutau(5,5),
     .          xneutup(5,5),xneutdow(5,5),xneutst(5,5),xneutch(5,5),
     .          xneutbot(5,5),xneuttop(5,5),xgluinoup(5),
     .          xgluinodo(5),xgluinoch(5),xgluinost(5),xgluinobot(5),
     .          xgluinotop(5),xchelne(5,2),xchmunmu(5,2),
     .          xchtauntau(5,2),xchubdow(5,2),xchcbs(5,2),xchtbb(5,2)

      DOUBLE PRECISION NS_ay,NS_by,NS_ax,NS_bx
      DOUBLE PRECISION NS_neutel,NS_neuttau,NS_neutnue,NS_neutnutau
      DOUBLE PRECISION NS_neutup,NS_neutdow,NS_neutbot,NS_neuttop
      DOUBLE PRECISION NS_chelne,NS_chtauntau,NS_chubd,NS_chtbb
      DOUBLE PRECISION NS_xgluibot,NS_xgluido,NS_xgluiup,NS_xgluitop
      DOUBLE PRECISION xmu1,xmu2,xmu3,sum1
*
      COMMON/NEUTRALINO_3GAMMA/xneutel,xneutmu,xneuttau,
     .          xneutnue,xneutnumu,xneutnutau,
     .          xneutup,xneutdow,xneutst,xneutch,
     .          xneutbot,xneuttop,xgluinoup,
     .          xgluinodo,xgluinoch,xgluinost,xgluinobot,
     .          xgluinotop,xchelne,xchmunmu,
     .          xchtauntau,xchubdow,xchcbs,xchtbb

      COMMON/NS_indices/ni,nj
      COMMON/NS_nx1/nx1t,ny1t
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
*
      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_neutel,NS_neuttau,NS_neutnue,NS_neutnutau
      EXTERNAL NS_neutup,NS_neutdow,NS_neutbot,NS_neuttop
      EXTERNAL NS_chelne,NS_chtauntau,NS_chubd,NS_chtbb
      EXTERNAL NS_xgluibot,NS_xgluido,NS_xgluiup,NS_xgluitop
      EXTERNAL NS_chmunmu,NS_neutmu,NS_neutnumu
c -------------------------------------------------------------------- c
c --------------------------- gluino up upbar ------------------------ c
c -------------------------------------------------------------------- c
      do ni=1,5
         xmu1=0d0
         xmu2=0d0
         xmu3=mgluino**2/amneut(ni)**2

         if(amneut(ni).gt.dabs(mgluino)) then
            CALL NS_integ2(NS_xgluiup,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                  xmu2,xmu3,nx1t,ny1t,sum1)
            xgluinoup(ni) = sum1/64d0/pi**3*amneut(ni)
         else
            xgluinoup(ni) = 0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c ------------------------- gluino down downbar ---------------------- c
c -------------------------------------------------------------------- c
      do ni = 1,5
         xmu1=0d0
         xmu2=0d0
         xmu3=mgluino**2/amneut(ni)**2

         if(amneut(ni).gt.dabs(mgluino)) then
            CALL NS_integ2(NS_xgluido,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                  xmu2,xmu3,nx1t,ny1t,sum1)
            xgluinodo(ni) = sum1/64d0/pi**3*amneut(ni)
         else
            xgluinodo(ni) = 0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c ------------------------ gluino charm charmbar --------------------- c
c -------------------------------------------------------------------- c
      do i=1,5
         xgluinoch(i) = xgluinoup(i)
      enddo
c -------------------------------------------------------------------- c
c ---------------------- gluino strange strangebar ------------------- c
c -------------------------------------------------------------------- c
      do i=1,5
         xgluinost(i) = xgluinodo(i)
      enddo
c -------------------------------------------------------------------- c
c ------------------------- gluino top topbar ------------------------ c
c -------------------------------------------------------------------- c
      do ni=1,5
         xmu1=amt**2/amneut(ni)**2
         xmu2=xmu1
         xmu3=mgluino**2/amneut(ni)**2

         if(amneut(ni).gt.(2d0*amt+dabs(mgluino))) then
            CALL NS_integ2(NS_xgluitop,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                  xmu2,xmu3,nx1t,ny1t,sum1)
            xgluinotop(ni) = sum1/64d0/pi**3*amneut(ni)
         else
            xgluinotop(ni) = 0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c --------------------------- gluino b bbar -------------------------- c
c -------------------------------------------------------------------- c
      do ni=1,5
         xmu1=amb**2/amneut(ni)**2
         xmu2=xmu1
         xmu3=mgluino**2/amneut(ni)**2

         if(amneut(ni).gt.(2d0*amb+dabs(mgluino))) then
            CALL NS_integ2(NS_xgluibot,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                  xmu2,xmu3,nx1t,ny1t,sum1)
            xgluinobot(ni) = sum1/64d0/pi**3*amneut(ni)
         else
            xgluinobot(ni) = 0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c --------------------------- neutralino e+ e- ----------------------- c
c -------------------------------------------------------------------- c
      do ni = 2,5,1
         do nj = 1,5
            xmu1=0d0
            xmu2=0d0
            xmu3=amneut(nj)**2/amneut(ni)**2

            if(nj.ge.ni) then
               xneutel(ni,nj) = 0d0
            else
               if(amneut(ni).gt.amneut(nj)) then
                  CALL NS_integ2(NS_neutel,NS_ax,NS_bx,NS_ay,NS_by,
     .                        xmu1,xmu2,xmu3,nx1t,ny1t,sum1)
                  xneutel(ni,nj) = sum1/64d0/(2*pi)**3*amneut(ni)
               else
                  xneutel(ni,nj) = 0d0
               endif
            endif
         enddo
      enddo

      do i=1,5
         xneutel(1,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c ------------------------ neutralino mu+ mu- ------------------------ c
c -------------------------------------------------------------------- c
      do ni = 2,5,1
         do nj = 1,5
            xmu1=0d0
            xmu2=0d0
            xmu3=amneut(nj)**2/amneut(ni)**2

            if(nj.ge.ni) then
               xneutmu(ni,nj) = 0d0
            else
               if(amneut(ni).gt.amneut(nj)) then
                  CALL NS_integ2(NS_neutmu,NS_ax,NS_bx,NS_ay,NS_by,
     .                        xmu1,xmu2,xmu3,nx1t,ny1t,sum1)
                  xneutmu(ni,nj) = sum1/64d0/(2*pi)**3*amneut(ni)
               else
                  xneutmu(ni,nj) = 0d0
               endif
            endif
         enddo
      enddo

      do i=1,5
         xneutmu(1,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c ------------------------ neutralino tau+ tau- ---------------------- c
c -------------------------------------------------------------------- c
      do ni = 2,5,1
         do nj = 1,5
            xmu1=amtau**2/amneut(ni)**2
            xmu2=xmu1
            xmu3=amneut(nj)**2/amneut(ni)**2

            if(nj.ge.ni) then
               xneuttau(ni,nj) = 0d0
            else
               if(amneut(ni).gt.(2d0*amtau+amneut(nj))) then
                  CALL NS_integ2(NS_neuttau,NS_ax,NS_bx,NS_ay,NS_by,
     .                        xmu1,xmu2,xmu3,nx1t,ny1t,sum1)
                  xneuttau(ni,nj) = sum1/64d0/(2*pi)**3*amneut(ni)
               else
                  xneuttau(ni,nj) = 0d0
               endif
            endif
         enddo
      enddo

      do i=1,5
         xneuttau(1,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c ----------------------- neutralino nu_e nu_ebar -------------------- c
c -------------------------------------------------------------------- c
      do ni = 2,5,1
         do nj = 1,5
            xmu1=0d0
            xmu2=0d0
            xmu3=amneut(nj)**2/amneut(ni)**2

            if(nj.ge.ni) then
               xneutnue(ni,nj) = 0d0
            else
               if(amneut(ni).gt.amneut(nj)) then
                  CALL NS_integ2(NS_neutnue,NS_ax,NS_bx,NS_ay,NS_by,
     .                        xmu1,xmu2,xmu3,nx1t,ny1t,sum1)
                  xneutnue(ni,nj) = sum1/64d0/(2*pi)**3*amneut(ni)
               else
                  xneutnue(ni,nj) = 0d0
               endif
            endif
         enddo
      enddo

      do i=1,5
         xneutnue(1,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c --------------------- neutralino nu_mu nu_mubar -------------------- c
c -------------------------------------------------------------------- c
      do ni = 2,5,1
         do nj = 1,5
            xmu1=0d0
            xmu2=0d0
            xmu3=amneut(nj)**2/amneut(ni)**2

            if(nj.ge.ni) then
               xneutnumu(ni,nj) = 0d0
            else
               if(amneut(ni).gt.amneut(nj)) then
                  CALL NS_integ2(NS_neutnumu,NS_ax,NS_bx,NS_ay,NS_by,
     .                        xmu1,xmu2,xmu3,nx1t,ny1t,sum1)
                  xneutnumu(ni,nj) = sum1/64d0/(2*pi)**3*amneut(ni)
               else
                  xneutnumu(ni,nj) = 0d0
               endif
            endif
         enddo
      enddo

      do i=1,5
         xneutnumu(1,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c -------------------- neutralino nu_tau nu_taubar ------------------- c
c -------------------------------------------------------------------- c
      do ni = 2,5,1
         do nj = 1,5
            xmu1=0d0
            xmu2=0d0
            xmu3=amneut(nj)**2/amneut(ni)**2

            if(nj.ge.ni) then
               xneutnutau(ni,nj) = 0d0
            else
               if(amneut(ni).gt.amneut(nj)) then
                  CALL NS_integ2(NS_neutnutau,NS_ax,NS_bx,NS_ay,NS_by,
     .                        xmu1,xmu2,xmu3,nx1t,ny1t,sum1)
                  xneutnutau(ni,nj) = sum1/64d0/(2*pi)**3*amneut(ni)
               else
                  xneutnutau(ni,nj) = 0d0
               endif
            endif
         enddo
      enddo

      do i=1,5
         xneutnutau(1,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c ----------------------- neutralino up upbar ------------------------ c
c -------------------------------------------------------------------- c
      do ni = 2,5,1
         do nj = 1,5
            xmu1=0d0
            xmu2=0d0
            xmu3=amneut(nj)**2/amneut(ni)**2

            if(nj.ge.ni) then
               xneutup(ni,nj) = 0d0
            else
               if(amneut(ni).gt.amneut(nj)) then
                  CALL NS_integ2(NS_neutup,NS_ax,NS_bx,NS_ay,NS_by,
     .                        xmu1,xmu2,xmu3,nx1t,ny1t,sum1)
                  xneutup(ni,nj) = sum1/64d0/(2*pi)**3*amneut(ni)
                  xneutup(ni,nj) = 3d0*xneutup(ni,nj)
               else
                  xneutup(ni,nj) = 0d0
               endif
            endif
         enddo
      enddo

      do i=1,5
         xneutup(1,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c -------------------- neutralino down downbar ----------------------- c
c -------------------------------------------------------------------- c
      do ni = 2,5,1
         do nj = 1,5
            xmu1=0d0
            xmu2=0d0
            xmu3=amneut(nj)**2/amneut(ni)**2

            if(nj.ge.ni) then
               xneutdow(ni,nj) = 0d0
            else
               if(amneut(ni).gt.amneut(nj)) then
                  CALL NS_integ2(NS_neutdow,NS_ax,NS_bx,NS_ay,NS_by,
     .                        xmu1,xmu2,xmu3,nx1t,ny1t,sum1)
                  xneutdow(ni,nj) = sum1/64d0/(2*pi)**3*amneut(ni)
                  xneutdow(ni,nj) = 3d0*xneutdow(ni,nj)
               else
                  xneutdow(ni,nj) = 0d0
               endif
            endif
         enddo
      enddo

      do i=1,5
         xneutdow(1,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c --------------------- neutralino charm charmbar -------------------- c
c -------------------------------------------------------------------- c
      do ni = 1,5
         xneutch(1,ni) = xneutup(1,ni)
         xneutch(2,ni) = xneutup(2,ni)
         xneutch(3,ni) = xneutup(3,ni)
         xneutch(4,ni) = xneutup(4,ni)
         xneutch(5,ni) = xneutup(5,ni)
      enddo
c -------------------------------------------------------------------- c
c ----------------- neutralino strange strangebar -------------------- c
c -------------------------------------------------------------------- c
      do ni = 1,5
         xneutst(1,ni) = xneutdow(1,ni)
         xneutst(2,ni) = xneutdow(2,ni)
         xneutst(3,ni) = xneutdow(3,ni)
         xneutst(4,ni) = xneutdow(4,ni)
         xneutst(5,ni) = xneutdow(5,ni)
      enddo
c -------------------------------------------------------------------- c
c ----------------------- neutralino top topbar ---------------------- c
c -------------------------------------------------------------------- c
      do ni = 2,5,1
         do nj = 1,5
            xmu1=amt**2/amneut(ni)**2
            xmu2=amt**2/amneut(ni)**2
            xmu3=amneut(nj)**2/amneut(ni)**2

            if(nj.ge.ni) then
               xneuttop(ni,nj) = 0d0
            else
               if(amneut(ni).gt.(2d0*amt+amneut(nj))) then
                  CALL NS_integ2(NS_neuttop,NS_ax,NS_bx,NS_ay,NS_by,
     .                        xmu1,xmu2,xmu3,nx1t,ny1t,sum1)
                  xneuttop(ni,nj) = sum1/64d0/(2*pi)**3*amneut(ni)
                  xneuttop(ni,nj) = 3d0*xneuttop(ni,nj)
               else
                  xneuttop(ni,nj) = 0d0
               endif
            endif
         enddo
      enddo

      do i=1,5
         xneuttop(1,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c -------------------------- neutralino b bbar ----------------------- c
c -------------------------------------------------------------------- c
      do ni = 2,5,1
         do nj = 1,5
            xmu1=amb**2/amneut(ni)**2
            xmu2=xmu1
            xmu3=amneut(nj)**2/amneut(ni)**2

            if(nj.ge.ni) then
               xneutbot(ni,nj) = 0d0
            else
               if(amneut(ni).gt.(2d0*amb+amneut(nj))) then
                  CALL NS_integ2(NS_neutbot,NS_ax,NS_bx,NS_ay,NS_by,
     .                        xmu1,xmu2,xmu3,nx1t,ny1t,sum1)
                  xneutbot(ni,nj) = sum1/64d0/(2*pi)**3*amneut(ni)
                  xneutbot(ni,nj) = 3d0*xneutbot(ni,nj)
               else
                  xneutbot(ni,nj) = 0d0
               endif
            endif
         enddo
      enddo

      do i=1,5
         xneutbot(1,i) = 0d0
      enddo
c -------------------------------------------------------------------- c
c ------------------------ chargino- e+ nuebar ----------------------- c
c -------------------------------------------------------------------- c
      do ni = 1,5
         do nj = 1,2
            xmu1=0d0
            xmu2=0d0
            xmu3=amchar(nj)**2/amneut(ni)**2

            if(amneut(ni).gt.amchar(nj)) then
               CALL NS_integ2(NS_chelne,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xchelne(ni,nj) = sum1/64d0/(2*pi)**3*amneut(ni)
            else
               xchelne(ni,nj) = 0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c ---------------------- chargino- mu+ numubar ----------------------- c
c -------------------------------------------------------------------- c
      do ni = 1,5
         do nj = 1,2
            xmu1=0d0
            xmu2=0d0
            xmu3=amchar(nj)**2/amneut(ni)**2

            if(amneut(ni).gt.amchar(nj)) then
               CALL NS_integ2(NS_chmunmu,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xchmunmu(ni,nj) = sum1/64d0/(2*pi)**3*amneut(ni)
            else
               xchmunmu(ni,nj) = 0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c --------------------- chargino- tau+ nutaubar ---------------------- c
c -------------------------------------------------------------------- c
      do ni = 1,5
         do nj = 1,2
            xmu1=0d0
            xmu2=amtau**2/amneut(ni)**2
            xmu3=amchar(nj)**2/amneut(ni)**2

            if(amneut(ni).gt.(amchar(nj)+amtau)) then
               CALL NS_integ2(NS_chtauntau,NS_ax,NS_bx,NS_ay,NS_by,
     .                     xmu1,xmu2,xmu3,nx1t,ny1t,sum1)
               xchtauntau(ni,nj) = sum1/64d0/(2*pi)**3*amneut(ni)
            else
               xchtauntau(ni,nj) = 0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c ------------------------ chargino- up downbar ---------------------- c
c -------------------------------------------------------------------- c
      do ni = 1,5
         do nj = 1,2
            xmu1=0d0
            xmu2=0d0
            xmu3=amchar(nj)**2/amneut(ni)**2

            if(amneut(ni).gt.amchar(nj)) then
               CALL NS_integ2(NS_chubd,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xchubdow(ni,nj) = 3d0*sum1/64d0/(2*pi)**3*amneut(ni)
            else
               xchubdow(ni,nj) = 0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c --------------------- chargino- charm strangebar ------------------- c
c -------------------------------------------------------------------- c
      do ni=1,5
         do nj=1,2
            xchcbs(ni,nj) = xchubdow(ni,nj)
         enddo
      enddo
c -------------------------------------------------------------------- c
c ---------------------- chargino- top bottombar --------------------- c
c -------------------------------------------------------------------- c
      do ni = 1,5
         do nj = 1,2
            xmu1=amt**2/amneut(ni)**2
            xmu2=amb**2/amneut(ni)**2
            xmu3=amchar(nj)**2/amneut(ni)**2

            if(amneut(ni).gt.(amchar(nj)+amt+amb)) then
               CALL NS_integ2(NS_chtbb,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xchtbb(ni,nj) = 3d0*sum1/64d0/(2*pi)**3*amneut(ni)
            else
               xchtbb(ni,nj) = 0d0
            endif
         enddo
      enddo
      end
c ==================================================================== c
c ============================ gluino up upbar ======================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_xgluiup(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,i,k
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION
     .gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),gdr(2),gdl(2)
      DOUBLE PRECISION dsup(2),dsupb(2)
      DOUBLE PRECISION aup(2,5),bup(2,5),ado(2,5),bdo(2,5)

      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION xmuneut1,x1,x2,x3,y1,y2,y3,xmusup(2)
*
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_coup10/aup,bup,ado,bdo
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_indices/ni,nj
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar

      xmuneut1 = mgluino**2/amneut(ni)**2

      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                               sup exchange
c -------------------------------------------------------------------- c
      xmusup(1) = asup1**2/amneut(ni)**2
      xmusup(2) = asup2**2/amneut(ni)**2

      dsup(1)  = 1d0-x1-xmusup(1)
      dsup(2)  = 1d0-x1-xmusup(2)
      dsupb(1) = 1d0-x2-xmusup(1)
      dsupb(2) = 1d0-x2-xmusup(2)

      NS_xgluiup=0d0

      if (dabs(mgluino).le.amneut(ni)) then
         do i=1,2
            do k=1,2
          if (xmusup(k).gt.1.0d0.and.xmusup(i).gt.1.0d0)then
               NS_xgluiup=NS_xgluiup
     .              +g2s*g3s/dsup(i)/dsup(k)*x1*y1
     .              *(aup(i,ni)*aup(k,ni)+bup(i,ni)*bup(k,ni))
     .              *(gur(i)*gur(k)+gul(i)*gul(k))
     .              +g2s*g3s/dsupb(k)/dsupb(i)*x2*y2
     .              *(aup(i,ni)*aup(k,ni)+bup(i,ni)*bup(k,ni))*
     .               (gur(i)*gur(k)+gul(i)*gul(k))
     .              +g2s*dsqrt(g3s)*dsqrt(g3s)/dsupb(i)/dsup(k)
     .              *( ( aup(i,ni)*bup(k,ni)*gur(i)*gul(k)
     .                 + gur(k)*gul(i)*aup(k,ni)*bup(i,ni))
     .                 *(-x1*y1-x2*y2+x3*y3)
     .               +2d0*mgluino/xmneut(ni)*y3*
     .               ( aup(i,ni)*aup(k,ni)*gul(k)*gul(i)
     .                +gur(i)*gur(k)*bup(k,ni)*bup(i,ni)))
          endif
            enddo
         enddo
      else
         NS_xgluiup=0d0
      endif

      end
c ==================================================================== c
c =========================== gluino down downbar ==================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_xgluido(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,i,k
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION
     .gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),gdr(2),gdl(2)
      DOUBLE PRECISION dsdo(2),dsdob(2)
      DOUBLE PRECISION aup(2,5),bup(2,5),ado(2,5),bdo(2,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION xmuneut1,x1,x2,x3,y1,y2,y3,xmusdo(2)
*
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/NS_coup10/aup,bup,ado,bdo
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS

      xmuneut1 = mgluino**2/amneut(ni)**2

      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                             sdown exchange
c -------------------------------------------------------------------- c
      xmusdo(1) = asdown1**2/amneut(ni)**2
      xmusdo(2) = asdown2**2/amneut(ni)**2

      dsdo(1)  = 1d0-x1-xmusdo(1)
      dsdo(2)  = 1d0-x1-xmusdo(2)
      dsdob(1) = 1d0-x2-xmusdo(1)
      dsdob(2) = 1d0-x2-xmusdo(2)

      NS_xgluido=0d0

      if (dabs(mgluino).le.amneut(ni)) then
         do i=1,2
            do k=1,2
         if (xmusdo(k).gt.1.0d0.and.xmusdo(i).gt.1.0d0)then
               NS_xgluido=NS_xgluido
     .              +g2s*g3s/dsdo(i)/dsdo(k)*x1*y1
     .              *(ado(i,ni)*ado(k,ni)+bdo(i,ni)*bdo(k,ni))
     .              *(gdr(i)*gdr(k)+gdl(i)*gdl(k))
     .              +g2s*dsdob(k)/dsdob(i)*x2*y2
     .              *(ado(i,ni)*ado(k,ni)+bdo(i,ni)*bdo(k,ni))*
     .               (gdr(i)*gdr(k)+gdl(i)*gdl(k))
     .              +g2s*dsqrt(g3s)*dsqrt(g3s)/dsdob(i)/dsdo(k)
     .              *( ( ado(i,ni)*bdo(k,ni)*gdr(i)*gdl(k)
     .                 + gdr(k)*gdl(i)*ado(k,ni)*bdo(i,ni))
     .                 *(-x1*y1-x2*y2+x3*y3)
     .               +2d0*mgluino/xmneut(ni)*y3*
     .               ( ado(i,ni)*ado(k,ni)*gdl(k)*gdl(i)
     .                +gdr(i)*gdr(k)*bdo(k,ni)*bdo(i,ni)))
          endif
            enddo
         enddo
      else
         NS_xgluido=0d0
      endif

      end
c ==================================================================== c
c ========================== gluino top topbar ======================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_xgluitop(x1,x2)
      IMPLICIT NONE
*
      INTEGER I,K,NI,NJ
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION
     .gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),gdr(2),gdl(2)
      DOUBLE PRECISION dst(2),dstb(2),sgn(5)
      DOUBLE PRECISION atopr(2,5),btopr(2,5)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      DOUBLE PRECISION xmuneut1,xmust(2),xmut,uh,th,x1,x2,gmst(2)
*
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_neutstoptop/atopr,btopr
      do i=1,5
         sgn(i) = 1d0
         if(xmneut(i).gt.0d0) then
            sgn(i) = 1d0
         elseif(xmneut(i).lt.0d0) then
            sgn(i) = -1d0
         endif
      enddo

      xmuneut1 = mgluino**2/amneut(ni)**2

c -------------------------------------------------------------------- c
c                              stop exchange
c -------------------------------------------------------------------- c
      gmst(1) = ast1
      gmst(2) = ast2
      xmust(1)   = ast1**2/amneut(ni)**2
      xmust(2)   = ast2**2/amneut(ni)**2
      xmut     = amt**2/amneut(ni)**2

      dst(1)  = 1d0-x1-xmust(1)+xmut
      dst(2)  = 1d0-x1-xmust(2)+xmut
      dstb(1) = 1d0-x2-xmust(1)+xmut
      dstb(2) = 1d0-x2-xmust(2)+xmut

      NS_xgluitop = 0d0

      uh = 1d0-x1+xmut
      th = 1d0-x2+xmut

      if ((dabs(mgluino)+2d0*amt).le.amneut(ni)) then
         do i=1,2
            do k=1,2
      if ((gmst(k)+amt).gt.amneut(ni).and.(gmst(i)+amt).gt.amneut(ni))
     .then
               NS_xgluitop=NS_xgluitop
     .              +g2s*g3s/dst(k)/dst(i)*(
     .               (atopr(i,ni)*btopr(k,ni)+btopr(i,ni)*atopr(k,ni))*
     .               (gtr(i)*gtr(k)+gtl(i)*gtl(k))*
     .               2d0*dsqrt(xmut)*sgn(ni)*(uh-xmut-xmuneut1)+
     .               (atopr(i,ni)*atopr(k,ni)+btopr(i,ni)*btopr(k,ni))*
     .               (gtr(i)*gtr(k)+gtl(i)*gtl(k))*
     .               (xmut*(-xmut-xmuneut1+2d0*uh-1d0)+uh*(-uh+
     .                xmuneut1+1d0)-xmuneut1)+
     .               (atopr(k,ni)*btopr(i,ni)+btopr(k,ni)*atopr(i,ni))*
     .               (gtr(k)*gtl(i)+gtr(i)*gtl(k))*mgluino/xmneut(ni)*
     .               xmut*(-4d0)+
     .               (atopr(k,ni)*atopr(i,ni)+btopr(k,ni)*btopr(i,ni))*
     .               (gtr(k)*gtl(i)+gtr(i)*gtl(k))*mgluino/xmneut(ni)*
     .               dsqrt(xmut)*sgn(ni)*2d0*(uh-xmut-1d0))
     .              +g2s*g3s/dstb(k)/dstb(i)*(
     .               (atopr(i,ni)*btopr(k,ni)+btopr(i,ni)*atopr(k,ni))*
     .               (gtr(i)*gtl(k)+gtr(k)*gtl(i))*mgluino/xmneut(ni)*
     .               xmut*(-4d0)+
     .               (atopr(i,ni)*atopr(k,ni)+btopr(i,ni)*btopr(k,ni))*
     .               (gtr(i)*gtl(k)+gtr(k)*gtl(i))*mgluino/xmneut(ni)*
     .               dsqrt(xmut)*sgn(ni)*2d0*(th-xmut-1d0)+
     .               (atopr(i,ni)*btopr(k,ni)+atopr(k,ni)*btopr(i,ni))*
     .               (gtr(i)*gtr(k)+gtl(k)*gtl(i))*
     .               dsqrt(xmut)*sgn(ni)*2d0*(th-xmuneut1-xmut)+
     .               (atopr(i,ni)*atopr(k,ni)+btopr(k,ni)*btopr(i,ni))*
     .               (gtr(i)*gtr(k)+gtl(k)*gtl(i))*(xmut*(2d0*th-1d0
     .               -xmut-xmuneut1)+th*(-th+xmuneut1+1d0)-xmuneut1))
     .              -2d0*g2s*dsqrt(g3s)*dsqrt(g3s)/dst(k)/dstb(i)*(
     .               (gtr(k)*gtr(i)*atopr(k,ni)*atopr(i,ni)
     .               +gtl(k)*gtl(i)*btopr(k,ni)*btopr(i,ni))*
     .               mgluino/xmneut(ni)*xmut*(-2d0)+
     .               (gtr(k)*gtr(i)*atopr(k,ni)*btopr(i,ni)
     .               +gtl(k)*gtl(i)*btopr(k,ni)*atopr(i,ni))*
     .               mgluino/xmneut(ni)*dsqrt(xmut)*sgn(ni)*
     .               (th-xmut-1d0)+
     .               (gtr(k)*gtr(i)*atopr(i,ni)*btopr(k,ni)+
     .                gtl(k)*gtl(i)*btopr(i,ni)*atopr(k,ni))*
     .               mgluino/xmneut(ni)*dsqrt(xmut)*sgn(ni)*
     .               (uh-xmut-1d0)+
     .               (gtr(k)*gtr(i)*btopr(i,ni)*btopr(k,ni)+
     .                gtl(k)*gtl(i)*atopr(i,ni)*atopr(k,ni))*
     .               mgluino/xmneut(ni)*(uh+th-xmuneut1-1d0)+
     .               (gtr(k)*gtl(i)*atopr(i,ni)*btopr(k,ni)+
     .                gtl(k)*gtr(i)*atopr(k,ni)*btopr(i,ni))*
     .               xmut*(uh+th-2d0*xmut)+
     .               (gtr(k)*gtl(i)*btopr(i,ni)*btopr(k,ni)+
     .               gtl(k)*gtr(i)*atopr(k,ni)*atopr(i,ni))*
     .               dsqrt(xmut)*sgn(ni)*(uh-xmuneut1-xmut)+
     .               (gtr(k)*gtl(i)*atopr(i,ni)*atopr(k,ni)+
     .               gtl(k)*gtr(i)*btopr(k,ni)*btopr(i,ni))*
     .               dsqrt(xmut)*sgn(ni)*(th-xmut-xmuneut1)+
     .               (gtr(k)*gtl(i)*atopr(k,ni)*btopr(i,ni)+
     .                gtl(k)*gtr(i)*atopr(i,ni)*btopr(k,ni))*
     .               (uh*th-xmut**2-xmuneut1))
               endif
            enddo
         enddo
      else
         NS_xgluitop=0d0
      endif

      end
c ==================================================================== c
c =========================== gluino b bbar ========================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_xgluibot(x1,x2)
*
      IMPLICIT NONE
      INTEGER I,K,NI,NJ
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION
     .gtr(2),gtl(2),gbr(2),gbl(2),gur(2),gul(2),gdr(2),gdl(2)
      DOUBLE PRECISION dsbo(2),dsbob(2),sgn(5)
      DOUBLE PRECISION abot(2,5),bbot(2,5)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION xmuneut1,xmusb(2),xmub,uh,th,x1,x2
*
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_neutsbotbot/abot,bbot
*
      do i=1,5
         sgn(i) = 1d0
         if(xmneut(i).gt.0d0) then
            sgn(i) = 1d0
         elseif(xmneut(i).lt.0d0) then
            sgn(i) = -1d0
         endif
      enddo

      xmuneut1 = mgluino**2/amneut(ni)**2
c -------------------------------------------------------------------- c
c 			    sbottom exchange
c -------------------------------------------------------------------- c
      xmusb(1)   = asb1**2/amneut(ni)**2
      xmusb(2)   = asb2**2/amneut(ni)**2
      xmub     = amb**2/amneut(ni)**2

      dsbo(1)  = 1d0-x1-xmusb(1)+xmub
      dsbo(2)  = 1d0-x1-xmusb(2)+xmub
      dsbob(1) = 1d0-x2-xmusb(1)+xmub
      dsbob(2) = 1d0-x2-xmusb(2)+xmub

      NS_xgluibot = 0d0

      uh = 1d0-x1+xmub
      th = 1d0-x2+xmub

      if ((dabs(mgluino)+2d0*amb).le.amneut(ni)) then
         do i=1,2
            do k=1,2
      if (xmusb(k).gt.1d0.and.xmusb(i).gt.1d0)then
               NS_xgluibot=NS_xgluibot
     .              +g2s*g3s/dsbo(k)/dsbo(i)*(
     .               (abot(i,ni)*bbot(k,ni)+bbot(i,ni)*abot(k,ni))*
     .               (gbr(i)*gbr(k)+gbl(i)*gbl(k))*
     .               2d0*dsqrt(xmub)*sgn(ni)*(uh-xmub-xmuneut1)+
     .               (abot(i,ni)*abot(k,ni)+bbot(i,ni)*bbot(k,ni))*
     .               (gbr(i)*gbr(k)+gbl(i)*gbl(k))*
     .               (xmub*(-xmub-xmuneut1+2d0*uh-1d0)+uh*(-uh+
     .                xmuneut1+1d0)-xmuneut1)+
     .               (abot(k,ni)*bbot(i,ni)+bbot(k,ni)*abot(i,ni))*
     .               (gbr(k)*gbl(i)+gbr(i)*gbl(k))*mgluino/xmneut(ni)*
     .               xmub*(-4d0)+
     .               (abot(k,ni)*abot(i,ni)+bbot(k,ni)*bbot(i,ni))*
     .               (gbr(k)*gbl(i)+gbr(i)*gbl(k))*mgluino/xmneut(ni)*
     .               dsqrt(xmub)*sgn(ni)*2d0*(uh-xmub-1d0))
     .              +g2s*g3s/dsbob(k)/dsbob(i)*(
     .               (abot(i,ni)*bbot(k,ni)+bbot(i,ni)*abot(k,ni))*
     .               (gbr(i)*gbl(k)+gbr(k)*gbl(i))*mgluino/xmneut(ni)*
     .               xmub*(-4d0)+
     .               (abot(i,ni)*abot(k,ni)+bbot(i,ni)*bbot(k,ni))*
     .               (gbr(i)*gbl(k)+gbr(k)*gbl(i))*mgluino/xmneut(ni)*
     .               dsqrt(xmub)*sgn(ni)*2d0*(th-xmub-1d0)+
     .               (abot(i,ni)*bbot(k,ni)+abot(k,ni)*bbot(i,ni))*
     .               (gbr(i)*gbr(k)+gbl(k)*gbl(i))*
     .               dsqrt(xmub)*sgn(ni)*2d0*(th-xmuneut1-xmub)+
     .               (abot(i,ni)*abot(k,ni)+bbot(k,ni)*bbot(i,ni))*
     .               (gbr(i)*gbr(k)+gbl(k)*gbl(i))*(xmub*(2d0*th-1d0
     .               -xmub-xmuneut1)+th*(-th+xmuneut1+1d0)-xmuneut1))
     .              -2d0*g2s*dsqrt(g3s)*dsqrt(g3s)/dsbo(k)/dsbob(i)*(
     .               (gbr(k)*gbr(i)*abot(k,ni)*abot(i,ni)
     .               +gbl(k)*gbl(i)*bbot(k,ni)*bbot(i,ni))*
     .               mgluino/xmneut(ni)*xmub*(-2d0)+
     .               (gbr(k)*gbr(i)*abot(k,ni)*bbot(i,ni)
     .               +gbl(k)*gbl(i)*bbot(k,ni)*abot(i,ni))*
     .               mgluino/xmneut(ni)*dsqrt(xmub)*
     .               sgn(ni)*(th-xmub-1d0)+
     .               (gbr(k)*gbr(i)*abot(i,ni)*bbot(k,ni)+
     .                gbl(k)*gbl(i)*bbot(i,ni)*abot(k,ni))*
     .               mgluino/xmneut(ni)*dsqrt(xmub)*sgn(ni)*
     .               (uh-xmub-1d0)+
     .               (gbr(k)*gbr(i)*bbot(i,ni)*bbot(k,ni)+
     .                gbl(k)*gbl(i)*abot(i,ni)*abot(k,ni))*
     .               mgluino/xmneut(ni)*(uh+th-xmuneut1-1d0)+
     .               (gbr(k)*gbl(i)*abot(i,ni)*bbot(k,ni)+
     .                gbl(k)*gbr(i)*abot(k,ni)*bbot(i,ni))*
     .               xmub*(uh+th-2d0*xmub)+
     .               (gbr(k)*gbl(i)*bbot(i,ni)*bbot(k,ni)+
     .               gbl(k)*gbr(i)*abot(k,ni)*abot(i,ni))*
     .               dsqrt(xmub)*sgn(ni)*(uh-xmuneut1-xmub)+
     .               (gbr(k)*gbl(i)*abot(i,ni)*abot(k,ni)+
     .               gbl(k)*gbr(i)*bbot(k,ni)*bbot(i,ni))*
     .               dsqrt(xmub)*sgn(ni)*(th-xmub-xmuneut1)+
     .               (gbr(k)*gbl(i)*abot(k,ni)*bbot(i,ni)+
     .                gbl(k)*gbr(i)*abot(i,ni)*bbot(k,ni))*
     .               (uh*th-xmub**2-xmuneut1))
      endif
            enddo
         enddo
      else
         NS_xgluibot=0d0
      endif


      end
c ==================================================================== c
c =========================  neutralino up upbar ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_neutup(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsup(2),dsupb(2)
      DOUBLE PRECISION opl(2,2),opr(2,2),oppl(5,5),oppr(5,5)
      DOUBLE PRECISION aup(2,5),bup(2,5),ado(2,5),bdo(2,5)
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION x1,x2,x3,y1,y2,y3,xmusup(2),xneutsup,
     .xmuz,dz,xneutzup,xneutzsup,xmuneut1
*
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup4/opl,opr,oppl,oppr
      COMMON/NS_coup10/aup,bup,ado,bdo
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS

      xmuneut1 = amneut(nj)**2/amneut(ni)**2

      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                              sup exchange
c -------------------------------------------------------------------- c
      xmusup(1) = asup1**2/amneut(ni)**2
      xmusup(2) = asup2**2/amneut(ni)**2

      dsup(1)  = 1d0-x1-xmusup(1)
      dsup(2)  = 1d0-x1-xmusup(2)
      dsupb(1) = 1d0-x2-xmusup(1)
      dsupb(2) = 1d0-x2-xmusup(2)

      xneutsup=0d0

      if (amneut(nj).le.amneut(ni)) then
         do i=1,2
            do k=1,2
       if (xmusup(k).gt.1d0.and.xmusup(i).gt.1d0)then
               xneutsup=xneutsup
     .              +g2s**2/dsup(i)/dsup(k)*x1*y1
     .              *(aup(i,ni)*aup(k,ni)+bup(i,ni)*bup(k,ni))
     .              *(aup(i,nj)*aup(k,nj)+bup(i,nj)*bup(k,nj))
     .              +g2s**2/dsupb(k)/dsupb(i)*x2*y2
     .              *(aup(i,ni)*aup(k,ni)+bup(i,ni)*bup(k,ni))
     .              *(aup(i,nj)*aup(k,nj)+bup(i,nj)*bup(k,nj))
     .              +g2s**2/dsupb(k)/dsup(i)
     .              *( ( aup(i,ni)*bup(k,ni)*aup(k,nj)*bup(i,nj)
     .                 + aup(i,nj)*bup(k,nj)*aup(k,ni)*bup(i,ni))
     .                 *(-x1*y1-x2*y2+x3*y3)
     .               +2d0*xmneut(nj)/xmneut(ni)*y3*
     .               ( aup(i,ni)*aup(k,ni)*aup(k,nj)*aup(i,nj)
     .                +bup(i,nj)*bup(k,nj)*bup(k,ni)*bup(i,ni)))
       endif
            enddo
         enddo
      else
         xneutsup=0d0
      endif

c -------------------------------------------------------------------- c
c                                Z exchange
c -------------------------------------------------------------------- c
      xmuz = mz**2/amneut(ni)**2
      dz   = y3-xmuz

      xneutzup=0d0

      if (amneut(nj).le.amneut(ni)) then
         if ((mz+amneut(nj)).gt.amneut(ni))then
         xneutzup=g2s**2*4d0/dz**2*
     .        (((azztoptop+vzztoptop)**2*oppl(ni,nj)**2
     .        +(azztoptop-vzztoptop)**2*oppr(ni,nj)**2)*x2*y2
     .        +((azztoptop+vzztoptop)**2*oppr(ni,nj)**2
     .        +(azztoptop-vzztoptop)**2*oppl(ni,nj)**2)*x1*y1
     .        -4d0*xmneut(nj)/xmneut(ni)*oppl(ni,nj)*oppr(ni,nj)
     .        *(azztoptop**2+vzztoptop**2)*y3 )
         endif
      else
         xneutzup=0d0
      endif
c -------------------------------------------------------------------- c
c                            Z-sup interference
c -------------------------------------------------------------------- c
      xneutzsup=0d0

      if (amneut(nj).le.amneut(ni)) then
         do i=1,2
      if ((mz+amneut(nj)).gt.amneut(ni).and.xmusup(i).gt.1d0)then
            xneutzsup=xneutzsup-g2s**2*4d0/dsup(i)/dz
     .    *((aup(i,ni)*aup(i,nj)*oppr(ni,nj)*(azztoptop+vzztoptop)+
     .     bup(i,ni)*bup(i,nj)*oppl(ni,nj)*(-azztoptop+vzztoptop))*x1*y1
     .    -(aup(i,nj)*aup(i,ni)*oppl(ni,nj)*(azztoptop+vzztoptop)
     .    +bup(i,nj)*bup(i,ni)*oppr(ni,nj)*(-azztoptop+vzztoptop)
     .     )*xmneut(nj)/xmneut(ni)*y3)
     .    +g2s**2*4d0/dsupb(i)/dz
     .    *((aup(i,ni)*aup(i,nj)*oppl(ni,nj)*(azztoptop+vzztoptop)+
     .     bup(i,ni)*bup(i,nj)*oppr(ni,nj)*(-azztoptop+vzztoptop))*x2*y2
     .    -(aup(i,nj)*aup(i,ni)*oppr(ni,nj)*(azztoptop+vzztoptop)
     .    +bup(i,nj)*bup(i,ni)*oppl(ni,nj)*(-azztoptop+vzztoptop)
     .     )*xmneut(nj)/xmneut(ni)*y3)
      endif
         enddo
      else
         xneutzsup=0d0
      endif

c -------------------------------------------------------------------- c

      NS_neutup = xneutsup+xneutzup+xneutzsup

      end
c ==================================================================== c
c ======================  neutralino down downbar ==================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_neutdow(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsdo(2),dsdob(2)
      DOUBLE PRECISION opl(2,2),opr(2,2),oppl(5,5),oppr(5,5)
      DOUBLE PRECISION aup(2,5),bup(2,5),ado(2,5),bdo(2,5)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION x1,x2,x3,y1,y2,y3,xmusdo(2),xneutsdow,
     .xmuz,dz,xneutzdow,xneutzsdow,xmuneut1
*
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup4/opl,opr,oppl,oppr
      COMMON/NS_coup10/aup,bup,ado,bdo
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS

      xmuneut1 = amneut(nj)**2/amneut(ni)**2

      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                            sdown exchange
c -------------------------------------------------------------------- c
      xmusdo(1) = asdown1**2/amneut(ni)**2
      xmusdo(2) = asdown2**2/amneut(ni)**2

      dsdo(1)  = 1d0-x1-xmusdo(1)
      dsdo(2)  = 1d0-x1-xmusdo(2)
      dsdob(1) = 1d0-x2-xmusdo(1)
      dsdob(2) = 1d0-x2-xmusdo(2)

      xneutsdow=0d0

      if (amneut(nj).le.amneut(ni)) then
         do i=1,2
            do k=1,2
      if (xmusdo(k).gt.1d0.and.xmusdo(i).gt.1d0)then
               xneutsdow=xneutsdow
     .              +g2s**2/dsdo(i)/dsdo(k)*x1*y1
     .              *(ado(i,ni)*ado(k,ni)+bdo(i,ni)*bdo(k,ni))
     .              *(ado(i,nj)*ado(k,nj)+bdo(i,nj)*bdo(k,nj))
     .              +g2s**2/dsdob(k)/dsdob(i)*x2*y2
     .              *(ado(i,ni)*ado(k,ni)+bdo(i,ni)*bdo(k,ni))*
     .               (ado(i,nj)*ado(k,nj)+bdo(i,nj)*bdo(k,nj))
     .              +g2s**2/dsdob(k)/dsdo(i)
     .              *( ( ado(i,ni)*bdo(k,ni)*ado(k,nj)*bdo(i,nj)
     .                 + ado(i,nj)*bdo(k,nj)*ado(k,ni)*bdo(i,ni))
     .                *(-x1*y1-x2*y2+x3*y3)
     .               + 2d0*xmneut(nj)/xmneut(ni)*y3*
     .               ( ado(i,ni)*ado(k,ni)*ado(k,nj)*ado(i,nj)
     .               + bdo(i,nj)*bdo(k,nj)*bdo(k,ni)*bdo(i,ni)))
       endif
            enddo
         enddo
      else
         xneutsdow=0d0
      endif

c -------------------------------------------------------------------- c
c                               Z exchange
c -------------------------------------------------------------------- c
      xmuz = mz**2/amneut(ni)**2
      dz   = y3-xmuz
      xneutzdow=0d0
      if (amneut(nj).le.amneut(ni)) then
      if((mz+amneut(nj)).gt.amneut(ni))then
         xneutzdow=g2s**2*4d0/dz**2*
     .        (((azzbotbot+vzzbotbot)**2*oppl(ni,nj)**2
     .        +(azzbotbot-vzzbotbot)**2*oppr(ni,nj)**2)*x2*y2
     .        +((azzbotbot+vzzbotbot)**2*oppr(ni,nj)**2
     .        +(azzbotbot-vzzbotbot)**2*oppl(ni,nj)**2)*x1*y1
     .        -4d0*xmneut(nj)/xmneut(ni)*oppl(ni,nj)*oppr(ni,nj)
     .        *(azzbotbot**2+vzzbotbot**2)*y3 )
      endif
      else
         xneutzdow=0d0
      endif
c -------------------------------------------------------------------- c
c                           Z-sdown interference
c -------------------------------------------------------------------- c
      xneutzsdow=0d0

      if (amneut(nj).le.amneut(ni)) then
         do i=1,2
      if ((mz+amneut(nj)).gt.amneut(ni).and.xmusdo(i).gt.1d0)then
            xneutzsdow=xneutzsdow-g2s**2*4d0/dsdo(i)/dz
     .           *((ado(i,ni)*ado(i,nj)*oppr(ni,nj)*
     .           (azzbotbot+vzzbotbot)+
     .           bdo(i,ni)*bdo(i,nj)*oppl(ni,nj)*
     .           (-azzbotbot+vzzbotbot))*x1*y1
     .           -(ado(i,nj)*ado(i,ni)*oppl(ni,nj)*
     .            (azzbotbot+vzzbotbot)
     .           +bdo(i,nj)*bdo(i,ni)*oppr(ni,nj)*
     .            (-azzbotbot+vzzbotbot)
     .           )*xmneut(nj)/xmneut(ni)*y3)
     .           +g2s**2*4d0/dsdob(i)/dz
     .           *((ado(i,ni)*ado(i,nj)*oppl(ni,nj)*
     .           (azzbotbot+vzzbotbot)+
     .           bdo(i,ni)*bdo(i,nj)*oppr(ni,nj)*
     .           (-azzbotbot+vzzbotbot))*x2*y2
     .           -(ado(i,nj)*ado(i,ni)*oppr(ni,nj)*
     .            (azzbotbot+vzzbotbot)
     .           +bdo(i,nj)*bdo(i,ni)*oppl(ni,nj)*
     .            (-azzbotbot+vzzbotbot)
     .           )*xmneut(nj)/xmneut(ni)*y3)
      endif
         enddo
      else
         xneutzsdow=0d0
      endif


      NS_neutdow = xneutsdow+xneutzdow+xneutzsdow

      end
c ==================================================================== c
c =======================  neutralino top topbar ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_neuttop(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i,J
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION atopr(2,5),btopr(2,5)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION opl(2,2),opr(2,2),oppl(5,5),oppr(5,5)
      DOUBLE PRECISION sgn(5)
      DOUBLE PRECISION dst(2),dstb(2)
      DOUBLE PRECISION hchichi(3,5,5),achichi(2,5,5)
      DOUBLE PRECISION hchachaR(3,2,2),hchachaL(3,2,2),achachaR(2,2,2),
     .     achachaL(2,2,2)
      DOUBLE PRECISION Httr(3),Attr(2)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION X1,x3,y3,x2,xmuneut1,xmust(2),xmut,uh,th,
     .vzz,azz,xneutstop,xmuz,dz,xneutztop,rh,sh,rk,xneutzstop
      DOUBLE PRECISION xneuthl(3,3),dhl(3),dhh(3),xneuta(2,2),
     .da(2),dna(2),xneuthlstop(3),xneutastop(2),xneutza(2),gmst(2)
*
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup4/opl,opr,oppl,oppr
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_coupNMSSM/hchichi,achichi,hchachaR,hchachaL,achachaR,
     .     achachaL
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/NS_neutstoptop/atopr,btopr
      COMMON/NS_phitoptop/Httr,Attr
*
      do i=1,5
         sgn(i) = 1d0
         if(xmneut(i).ge.0d0) then
            sgn(i) = 1d0
         elseif(xmneut(i).lt.0d0) then
            sgn(i) = -1d0
         endif
      enddo
      gmst(1)=ast1
      gmst(2)=ast2
      xmuneut1 = amneut(nj)**2/amneut(ni)**2
      xmust(1)   = ast1**2/amneut(ni)**2
      xmust(2)   = ast2**2/amneut(ni)**2
      xmut     = amt**2/amneut(ni)**2

      x3 = 2d0-x1-x2
      y3 = 1d0+xmuneut1-x3

      uh = 1d0-x1+xmut
      th = 1d0-x2+xmut

      vzz = vzztoptop
      azz = azztoptop
c -------------------------------------------------------------------- c
c                            stop exchange
c -------------------------------------------------------------------- c
      dst(1)  = 1d0-x1-xmust(1)+xmut
      dst(2)  = 1d0-x1-xmust(2)+xmut
      dstb(1) = 1d0-x2-xmust(1)+xmut
      dstb(2) = 1d0-x2-xmust(2)+xmut

      xneutstop=0d0

      if ((amneut(nj)+2d0*amt).le.amneut(ni)) then
         do i=1,2
            do k=1,2
      if ((gmst(k)+amt).gt.amneut(ni).and.(gmst(i)+amt).gt.amneut(ni)
     .)then
               xneutstop=xneutstop
     .          +g2s**2/dst(k)/dst(i)*(
     .           (atopr(i,ni)*btopr(k,ni)+btopr(i,ni)*atopr(k,ni))*
     .           (atopr(i,nj)*btopr(k,nj)+btopr(i,nj)*atopr(k,nj))*
     .           xmneut(nj)/xmneut(ni)*xmut*(-4d0)+
     .           (atopr(i,ni)*atopr(k,ni)+btopr(i,ni)*btopr(k,ni))*
     .           (atopr(i,nj)*btopr(k,nj)+btopr(i,nj)*atopr(k,nj))*
     .           xmneut(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*2d0*
     .           (uh-xmut-1d0)+
     .           (atopr(i,ni)*btopr(k,ni)+btopr(i,ni)*atopr(k,ni))*
     .           (atopr(i,nj)*atopr(k,nj)+btopr(i,nj)*btopr(k,nj))*
     .           dsqrt(xmut)*sgn(ni)*2d0*(uh-xmut-xmuneut1)+
     .           (atopr(i,ni)*atopr(k,ni)+btopr(i,ni)*btopr(k,ni))*
     .           (atopr(i,nj)*atopr(k,nj)+btopr(i,nj)*btopr(k,nj))*
     .           (-uh**2+uh*(1d0+xmuneut1+2d0*xmut)-(xmuneut1+xmut)*
     .           (1d0+xmut)))
     .           +g2s**2/dstb(k)/dstb(i)*(
     .           (atopr(i,ni)*btopr(k,ni)+btopr(i,ni)*atopr(k,ni))*
     .           (atopr(i,nj)*btopr(k,nj)+btopr(i,nj)*atopr(k,nj))*
     .           xmneut(nj)/xmneut(ni)*xmut*(-4d0)+
     .           (atopr(i,ni)*atopr(k,ni)+btopr(i,ni)*btopr(k,ni))*
     .           (atopr(i,nj)*btopr(k,nj)+btopr(i,nj)*atopr(k,nj))*
     .           xmneut(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*2d0*
     .           (th-xmut-1d0)+
     .           (atopr(i,ni)*btopr(k,ni)+btopr(i,ni)*atopr(k,ni))*
     .           (atopr(i,nj)*atopr(k,nj)+btopr(i,nj)*btopr(k,nj))*
     .           dsqrt(xmut)*sgn(ni)*2d0*(th-xmut-xmuneut1)+
     .           (atopr(i,ni)*atopr(k,ni)+btopr(i,ni)*btopr(k,ni))*
     .           (atopr(i,nj)*atopr(k,nj)+btopr(i,nj)*btopr(k,nj))*
     .           (-th**2+th*(1d0+xmuneut1+2d0*xmut)-(xmuneut1+xmut)*
     .           (1d0+xmut)))
     .           -2d0*g2s**2/dst(k)/dstb(i)*(
     .           (btopr(i,ni)*btopr(k,ni)*atopr(i,nj)*atopr(k,nj)
     .           +atopr(i,ni)*atopr(k,ni)*btopr(i,nj)*btopr(k,nj))*
     .           xmneut(nj)/xmneut(ni)*xmut*(-2d0)+
     .           (atopr(i,ni)*btopr(k,ni)*atopr(i,nj)*atopr(k,nj)
     .           +atopr(k,ni)*btopr(i,ni)*btopr(i,nj)*btopr(k,nj))*
     .           xmneut(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*
     .           (th-xmut-1d0)+
     .           (atopr(k,ni)*btopr(i,ni)*atopr(i,nj)*atopr(k,nj)
     .           +atopr(i,ni)*btopr(k,ni)*btopr(i,nj)*btopr(k,nj))*
     .           xmneut(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*
     .           (uh-xmut-1d0)+
     .           (atopr(k,ni)*atopr(i,ni)*atopr(i,nj)*atopr(k,nj)
     .           +btopr(i,ni)*btopr(k,ni)*btopr(i,nj)*btopr(k,nj))*
     .           xmneut(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0)+
     .           (atopr(k,ni)*btopr(i,ni)*atopr(k,nj)*btopr(i,nj)
     .           +btopr(k,ni)*atopr(i,ni)*atopr(i,nj)*btopr(k,nj))*
     .           xmut*(uh+th-2d0*xmut)+
     .           (atopr(k,ni)*atopr(i,ni)*atopr(k,nj)*btopr(i,nj)
     .           +btopr(k,ni)*btopr(i,ni)*atopr(i,nj)*btopr(k,nj))*
     .           dsqrt(xmut)*sgn(ni)*(uh-xmut-xmuneut1)+
     .           (btopr(k,ni)*btopr(i,ni)*atopr(k,nj)*btopr(i,nj)
     .           +atopr(k,ni)*atopr(i,ni)*atopr(i,nj)*btopr(k,nj))*
     .           dsqrt(xmut)*sgn(ni)*(th-xmut-xmuneut1)+
     .           (btopr(k,ni)*atopr(i,ni)*atopr(k,nj)*btopr(i,nj)
     .           +atopr(k,ni)*btopr(i,ni)*atopr(i,nj)*btopr(k,nj))*
     .           (uh*th-xmut**2-xmuneut1))
       endif
            enddo
         enddo
      else
         xneutstop=0d0
      endif

c -------------------------------------------------------------------- c
c 	                     Z exchange
c -------------------------------------------------------------------- c
      xmuz = mz**2/amneut(ni)**2
      dz   = y3-xmuz

      xneutztop=0d0

      rh = xmuneut1+2d0*xmut-th-uh+1d0
      sh = (xmuneut1-th-uh+1d0)*2d0*xmut+4d0*xmut**2
      rk = xmuneut1*(2d0*xmut-th-uh+4d0)+2d0*xmut-uh-th

      if ((amneut(nj)+2d0*amt).le.amneut(ni)) then
         if ((amneut(nj)+mz).gt.amneut(ni))then
         xneutztop=xneutztop+g2s**2/dz**2*(
     .    oppl(ni,nj)*oppr(ni,nj)*(vzz**2-azz**2)*
     .    xmneut(nj)/xmneut(ni)*xmut*(-16d0/xmuz**2*rh**2+
     .    32d0/xmuz*rh-64d0)+
     .    oppl(ni,nj)*oppr(ni,nj)*(vzz**2+azz**2)*
     .    xmneut(nj)/xmneut(ni)*(8d0/xmuz**2*rh*sh-16d0/xmuz*sh
     .    -16d0*(xmuneut1-uh-th+1d0))+
     .    (oppl(ni,nj)**2+oppr(ni,nj)**2)*(vzz**2-azz**2)*
     .    xmut*(4d0/xmuz**2*rh*rk-8d0/xmuz*rk+8d0*(uh+th-2d0*xmut))
     .    +(oppl(ni,nj)**2+oppr(ni,nj)**2)*(vzz**2+azz**2)*
     .    (-2d0/xmuz**2*rk*sh+8d0/xmuz*(xmuneut1*(2d0*xmut**2+
     .    4d0*xmut-xmut*(th+uh))+2d0*xmut**2-xmut*(uh+th))+4d0*(
     .    xmuneut1*(uh+th-2d0*xmut-2d0)+2d0*xmut*(uh+th-1d0)
     .    -2d0*xmut**2+th*(-th+1d0)+uh*(-uh+1d0)))+
     .    (oppl(ni,nj)**2-oppr(ni,nj)**2)*vzz*azz*8d0*(
     .    xmuneut1*(th-uh)+2d0*xmut*(th-uh)+th*(-th+1d0)+uh*(uh-1d0))
     .    )
        endif
      else
         xneutztop=0d0
      endif
c -------------------------------------------------------------------- c
c     NMSSM CP EVEN Higgs exchange + Interference H(i)-H(j)
c -------------------------------------------------------------------- c
       do i=1,3
       do j=1,3

        dhl(i)   = y3-SMASS(i)**2/amneut(ni)**2
        dhh(j)   = y3-SMASS(j)**2/amneut(ni)**2

      xneuthl(i,j)=0d0

      if ((amneut(nj)+2d0*amt).le.amneut(ni)) then
         if ((SMASS(i)+amneut(nj)).gt.amneut(ni).and.
     .(SMASS(j)+amneut(nj)).gt.amneut(ni))then
         xneuthl(i,j)=g2s**2/dhl(i)/dhh(j)*Httr(i)*Httr(j)*
     .    hchichi(i,ni,nj)*hchichi(j,ni,nj)*(
     .    xmneut(nj)/xmneut(ni)*xmut*(-32d0)+
     .    xmneut(nj)/xmneut(ni)*16d0*(1d0+xmuneut1-th-uh)+
     .    xmut*16d0*(2d0*xmut-th-uh)+
     .    8d0*(xmuneut1*(uh+th-2d0*xmut)+2d0*xmut*(uh+th-1d0)
     .    -(th+uh)**2+uh+th))
      endif
      else
         xneuthl(i,j)=0d0
      endif


      enddo
      enddo
c -------------------------------------------------------------------- c
c            NMSSM CP ODD Higgs exchange + Interference A(i)-A(j)
c -------------------------------------------------------------------- c
      do i=1,2
      do j=1,2

      da(i)    = y3-PMASS(i)**2/amneut(ni)**2
      dna(j)   = y3-PMASS(j)**2/amneut(ni)**2

      xneuta(i,j)=0d0

      if ((amneut(nj)+2d0*amt).le.amneut(ni)) then
      if ((PMASS(i)+amneut(nj)).gt.amneut(ni).and.
     .(PMASS(j)+amneut(nj)).gt.amneut(ni))then
         xneuta(i,j)=g2s**2/da(i)/dna(j)*Attr(i)*Attr(j)*
     .    achichi(i,ni,nj)*achichi(j,ni,nj)*(
     .    xmneut(nj)/xmneut(ni)*xmut*(-32d0)+
     .    xmneut(nj)/xmneut(ni)*(-16d0)*(1d0+xmuneut1-th-uh)+
     .    xmut*(-16d0)*(2d0*xmut-th-uh)+
     .    8d0*(xmuneut1*(uh+th-2d0*xmut)+2d0*xmut*(uh+th-1d0)
     .    -(th+uh)**2+uh+th))
      endif
      else
         xneuta(i,j)=0d0
      endif

      enddo
      enddo

c -------------------------------------------------------------------- c
c    	                 interference Z-stop
c -------------------------------------------------------------------- c
      xneutzstop=0d0

      if ((amneut(nj)+2d0*amt).le.amneut(ni)) then
         do i=1,2
       if ((gmst(i)+amt).gt.amneut(ni).and.(mz+amneut(nj)).gt.
     .amneut(ni))then
            xneutzstop=xneutzstop
     .      +g2s**2/dst(i)/dz*oppl(ni,nj)*(
     .      ((atopr(i,ni)*atopr(i,nj)-btopr(i,ni)*btopr(i,nj))*vzz-
     .       (atopr(i,ni)*atopr(i,nj)+btopr(i,ni)*btopr(i,nj))*azz)*
     .      xmneut(nj)/xmneut(ni)*xmut*(
     .      1/xmuz*(-4d0)*(1d0+xmuneut1+2d0*xmut-uh-th)+16d0) +
     .      ((atopr(i,ni)*atopr(i,nj)-btopr(i,ni)*btopr(i,nj))*vzz+
     .       (atopr(i,ni)*atopr(i,nj)+btopr(i,ni)*btopr(i,nj))*azz)*
     .      xmneut(nj)/xmneut(ni)*(2d0/xmuz*((xmuneut1+1d0-uh-th)*
     .      2d0*xmut+4d0*xmut**2)+4d0*(1d0+xmuneut1-uh-th)) +
     .      ((btopr(i,ni)*btopr(i,nj)-atopr(i,ni)*atopr(i,nj))*vzz+
     .       (btopr(i,ni)*btopr(i,nj)+atopr(i,ni)*atopr(i,nj))*azz)*
     .      xmut*(2d0/xmuz*(xmuneut1*(2d0*xmut-th-uh+4d0)+2d0*xmut
     .      -th-uh)+4d0*(2d0*xmut-th-uh)) +
     .      ((btopr(i,ni)*btopr(i,nj)-atopr(i,ni)*atopr(i,nj))*vzz-
     .       (btopr(i,ni)*btopr(i,nj)+atopr(i,ni)*atopr(i,nj))*azz)*
     .      (2d0/xmuz*(xmuneut1*(-2d0*xmut**2+xmut*th-2d0*xmut+xmut*
     .       uh-2d0*xmut)+xmut*(-2d0*xmut+uh)+xmut*th)+4d0*(
     .      xmuneut1*(xmut-uh+1d0)+xmut*(xmut-uh)+xmut*(1d0-uh)+uh**2
     .      -uh)) +
     .      (atopr(i,nj)*btopr(i,ni)*(vzz+azz)-
     .       atopr(i,ni)*btopr(i,nj)*(vzz-azz))*
     .      xmneut(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*(2d0/xmuz*
     .      (xmuneut1*(1d0+xmut-uh)+1d0+xmut*(xmut-uh)+xmut*
     .      (xmut-th-2d0*uh+3d0)+th*(uh-1d0)+uh*(uh-2d0))+
     .      4d0*(1d0+xmut-th)) +
     .      (atopr(i,nj)*btopr(i,ni)*(vzz-azz)-
     .       atopr(i,ni)*btopr(i,nj)*(vzz+azz))*
     .      xmneut(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*
     .      (2d0/xmuz*(xmuneut1*(-1d0-xmut+uh)-1d0+xmut*(-1d0+uh)
     .      +uh*(2d0-th-uh)+th+xmut*(-2d0*xmut+th+2d0*uh-2d0))+
     .      8d0*(1d0+xmut-uh)) +
     .      (atopr(i,ni)*btopr(i,nj)*(vzz-azz)-
     .       atopr(i,nj)*btopr(i,ni)*(vzz+azz))*dsqrt(xmut)*sgn(ni)*(
     .      (-2d0)/xmuz*(xmuneut1-uh+xmut)*(xmuneut1+2d0*xmut-th-uh
     .      +1d0)+8d0*(xmuneut1+xmut-uh)) +
     .      (atopr(i,ni)*btopr(i,nj)*(vzz+azz)-
     .       atopr(i,nj)*btopr(i,ni)*(vzz-azz))*dsqrt(xmut)*sgn(ni)*(
     .      2d0/xmuz*(xmuneut1*(xmuneut1+3d0*xmut-th-2d0*uh+1d0)+
     .      xmut*(2d0*xmut-th-2d0*uh)+xmut*(1d0-uh)+uh*th+uh**2-uh)+
     .      4d0*(xmuneut1+xmut-th)) )
     .      -g2s**2/dstb(i)/dz*oppl(ni,nj)*(
     .      ((-atopr(i,ni)*atopr(i,nj)+btopr(i,ni)*btopr(i,nj))*vzz+
     .       (atopr(i,ni)*atopr(i,nj)+btopr(i,ni)*btopr(i,nj))*azz)*
     .      xmneut(nj)/xmneut(ni)*xmut*(
     .      1/xmuz*(-4d0)*(1d0+xmuneut1+2d0*xmut-uh-th)+16d0) +
     .      ((-atopr(i,ni)*atopr(i,nj)+btopr(i,ni)*btopr(i,nj))*vzz-
     .       (atopr(i,ni)*atopr(i,nj)+btopr(i,ni)*btopr(i,nj))*azz)*
     .      xmneut(nj)/xmneut(ni)*(2d0/xmuz*(2d0*xmut*(xmuneut1+1d0
     .      -th-uh)+4d0*xmut**2)+4d0*(1d0+xmuneut1-uh-th)) +
     .      ((atopr(i,ni)*atopr(i,nj)-btopr(i,ni)*btopr(i,nj))*vzz-
     .       (atopr(i,ni)*atopr(i,nj)+btopr(i,ni)*btopr(i,nj))*azz)*
     .      xmut*(2d0/xmuz*(xmuneut1*(2d0*xmut-th-uh+4d0)+2d0*xmut
     .      -th-uh)+4d0*(2d0*xmut-th-uh)) +
     .      ((atopr(i,ni)*atopr(i,nj)-btopr(i,ni)*btopr(i,nj))*vzz+
     .       (atopr(i,ni)*atopr(i,nj)+btopr(i,ni)*btopr(i,nj))*azz)*
     .      (2d0/xmuz*(xmuneut1*(-2d0*xmut**2+xmut*th-2d0*xmut+xmut*
     .       uh-2d0*xmut)+xmut*(-2d0*xmut+uh)+xmut*th)+4d0*(
     .      xmuneut1*(xmut-th+1d0)+xmut*(xmut-th)+xmut*(1d0-th)+th**2
     .      -th)) +
     .      (atopr(i,nj)*btopr(i,ni)*(-vzz-azz)+
     .       atopr(i,ni)*btopr(i,nj)*(vzz-azz))*
     .      xmneut(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*
     .      (2d0/xmuz*(xmuneut1*(1d0+xmut-th)+1d0+xmut*(2d0*xmut
     .      -2d0*th-uh+3d0)+th*(th-xmut+uh-2d0)-uh)+4d0*(1d0+xmut
     .      -uh)) +
     .      (atopr(i,nj)*btopr(i,ni)*(-vzz+azz)+
     .       atopr(i,ni)*btopr(i,nj)*(vzz+azz))*
     .      xmneut(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*
     .      (2d0/xmuz*(xmuneut1*(-1d0-xmut+th)-1d0+xmut*(-xmut+th
     .      -1d0)+xmut*(2d0*th+uh-2d0-xmut)-th*(th+uh)+2d0*th+uh)+
     .      8d0*(1d0+xmut-th)) +
     .      (atopr(i,ni)*btopr(i,nj)*(-vzz+azz)+
     .       atopr(i,nj)*btopr(i,ni)*(vzz+azz))*dsqrt(xmut)*sgn(ni)*(
     .      2d0/xmuz*(xmuneut1*(-xmuneut1-3d0*xmut+2d0*th+uh-1d0)
     .      +xmut*(-2d0*xmut+2d0*th+uh-1d0)+th*(xmut-th-uh+1d0))
     .      +8d0*(xmuneut1+xmut-th)) +
     .      (atopr(i,ni)*btopr(i,nj)*(-vzz-azz)+
     .       atopr(i,nj)*btopr(i,ni)*(vzz-azz))*dsqrt(xmut)*sgn(ni)*(
     .      2d0/xmuz*(xmuneut1*(xmuneut1+3d0*xmut-uh-2d0*th+1d0)+
     .      xmut*(xmut-th+1d0)+xmut*(xmut-2d0*th-uh)+uh*th+th**2-th)+
     .      4d0*(xmuneut1+xmut-uh)) )
            endif
         enddo
      else
         xneutzstop=0d0
      endif

c -------------------------------------------------------------------- c
c                        interference H(j)-stop
c -------------------------------------------------------------------- c
      do j=1,3

        dhl(j)   = y3-SMASS(j)**2/amneut(ni)**2
        xneuthlstop(j)=0d0

      if ((amneut(nj)+2d0*amt).le.amneut(ni)) then
         do i=1,2
            if ((SMASS(j)+amneut(nj)).gt.amneut(ni).and.
     .(gmst(i)+amt).gt.amneut(ni))then
            xneuthlstop(j)=xneuthlstop(j)
     .       +2d0*g2s**2/dhl(j)/dst(i)*(Httr(j)/dsqrt(2d0))*
     .       (-2d0)*hchichi(j,ni,nj)*(
     .       (atopr(i,nj)*atopr(i,ni)+btopr(i,ni)*btopr(i,nj))*(
     .       dsqrt(xmut)*sgn(ni)*(xmuneut1+xmut-uh) +
     .       dsqrt(xmut)*sgn(ni)*(-xmuneut1-xmut+th) +
     .       xmneut(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*(-1d0-xmut+th) +
     .       xmneut(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*(1d0+xmut-uh) )
     .       +(atopr(i,nj)*btopr(i,ni)+atopr(i,ni)*btopr(i,nj))*(
     .       xmneut(nj)/xmneut(ni)*(-1d0-xmuneut1+uh+th) +
     .       2d0*xmneut(nj)/xmneut(ni)*xmut +
     .       xmut*(-2d0*xmut+th+uh) +
     .       (uh**2+uh*th-uh*(1d0+2d0*xmut+xmuneut1)+xmut+xmuneut1*
     .       xmut)) )
     .       +2d0*g2s**2/dhl(j)/dstb(i)*(Httr(j)/dsqrt(2d0))*
     .       (-2d0)*hchichi(j,ni,nj)*(
     .       (atopr(i,nj)*atopr(i,ni)+btopr(i,ni)*btopr(i,nj))*(
     .       dsqrt(xmut)*sgn(ni)*(uh-xmuneut1-xmut) +
     .       dsqrt(xmut)*sgn(ni)*(-th+xmut+xmuneut1) +
     .       dsqrt(xmut)*sgn(ni)*xmneut(nj)/xmneut(ni)*(-th+1d0+xmut) +
     .       dsqrt(xmut)*sgn(ni)*xmneut(nj)/xmneut(ni)*(uh-1d0-xmut)) +
     .       (atopr(i,nj)*btopr(i,ni)+atopr(i,ni)*btopr(i,nj))*(
     .       2d0*xmut*xmneut(nj)/xmneut(ni) +
     .       (uh*th+th**2-th*(1d0+2d0*xmut+xmuneut1)+xmut+
     .       xmut*xmuneut1) +
     .       xmneut(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0) +
     .       xmut*(uh+th-2d0*xmut)) )
       endif
         enddo
      else
         xneuthlstop(j)=0d0
      endif

      enddo
c -------------------------------------------------------------------- c
c 	                 interference A(j)-stop
c -------------------------------------------------------------------- c
      do j=1,2
      da(j)    = y3-PMASS(j)**2/amneut(ni)**2
      xneutastop(j)=0d0

      if ((amneut(nj)+2d0*amt).le.amneut(ni)) then
         do i=1,2
      if ((PMASS(j)+amneut(nj)).gt.amneut(ni).and.
     . (gmst(i)+amt).gt.amneut(ni))then
            xneutastop(j)=xneutastop(j)
     .       +2d0*g2s**2/da(j)/dst(i)*(Attr(j)/dsqrt(2d0))*
     .       2d0*achichi(j,ni,nj)*(
     .       (atopr(i,nj)*atopr(i,ni)+btopr(i,ni)*btopr(i,nj))*(
     .       dsqrt(xmut)*sgn(ni)*(xmuneut1+xmut-uh) +
     .       dsqrt(xmut)*sgn(ni)*(-xmuneut1-xmut+th)*(-1d0) +
     .       xmneut(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*(-1d0-xmut+th)*
     .       (-1d0) +
     .       xmneut(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*(1d0+xmut-uh) )
     .       +(atopr(i,nj)*btopr(i,ni)+atopr(i,ni)*btopr(i,nj))*(
     .       xmneut(nj)/xmneut(ni)*(-1d0-xmuneut1+uh+th)*(-1d0) +
     .       2d0*xmneut(nj)/xmneut(ni)*xmut +
     .       xmut*(-2d0*xmut+th+uh)*(-1d0) +
     .       (uh**2+uh*th-uh*(1d0+2d0*xmut+xmuneut1)+xmut+xmuneut1*
     .       xmut)) )
     .       +2d0*g2s**2/da(j)/dstb(i)*(Attr(j)/dsqrt(2d0))*
     .       2d0*achichi(j,ni,nj)*(
     .       (atopr(i,nj)*atopr(i,ni)+btopr(i,ni)*btopr(i,nj))*(
     .       dsqrt(xmut)*sgn(ni)*(uh-xmuneut1-xmut)*(-1d0) +
     .       dsqrt(xmut)*sgn(ni)*(-th+xmut+xmuneut1) +
     .       dsqrt(xmut)*sgn(ni)*xmneut(nj)/xmneut(ni)*(-th+1d0+xmut) +
     .       dsqrt(xmut)*sgn(ni)*xmneut(nj)/xmneut(ni)*(uh-1d0-xmut)*
     .       (-1d0)) +
     .       (atopr(i,nj)*btopr(i,ni)+atopr(i,ni)*btopr(i,nj))*(
     .       2d0*xmut*xmneut(nj)/xmneut(ni) +
     .       (uh*th+th**2-th*(1d0+2d0*xmut+xmuneut1)+xmut+
     .       xmut*xmuneut1) +
     .       xmneut(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0)*(-1d0) +
     .       xmut*(uh+th-2d0*xmut)*(-1d0)) )
       endif
         enddo
      else
         xneutastop(j)=0d0
      endif

      enddo
c -------------------------------------------------------------------- c
c 	                interference Z and A(i)
C ---------------- INTERFERENCE TERMS Z-h/H ARE MISSING
c -------------------------------------------------------------------- c
      do i=1,2
      da(i)    = y3-PMASS(i)**2/amneut(ni)**2
      xneutza(i)=0d0

      if ((amneut(nj)+amt+amt).le.amneut(ni)) then
      if ((PMASS(i)+amneut(nj)).gt.amneut(ni).and.
     .(amneut(nj)+mz).gt.amneut(ni))then
         xneutza(i)=xneutza(i)-4d0*g2s**2/da(i)/dz*azz*
     .    Attr(i)/dsqrt(2d0)*2d0*achichi(i,ni,nj)*oppl(ni,nj)*(
     .    xmneut(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*(
     .    4d0/xmuz*(2d0+xmuneut1*(2d0*xmut-th-uh+2d0)+xmut*(-uh-th
     .    +1d0+xmut)+xmut*(3d0*xmut-3d0*(th+uh)+5d0)+(th+uh)**2
     .    -3d0*(th+uh))+
     .    4d0*(th-xmut-1d0)+
     .    4d0*(uh-xmut-1d0))
     .    -dsqrt(xmut)*sgn(ni)*(
     .    4d0/xmuz*(-2d0*xmuneut1**2+xmuneut1*(-6d0*xmut+
     .    3d0*(th+uh)-2d0)+xmut*(-xmut+th+uh-1d0)+xmut*(-3d0*xmut
     .    +3d0*(uh+th)-1d0)-(th+uh)**2+(th+uh))+
     .    4d0*(xmuneut1+xmut-uh)+
     .    4d0*(xmuneut1+xmut-th)) )
      endif
      else
         xneutza(i)=0d0
      endif

      enddo

      NS_neuttop=xneutztop+xneutstop+xneutzstop

      do i=1,3
         NS_neuttop=NS_neuttop+xneuthlstop(i)
         do j=1,3
            NS_neuttop=NS_neuttop+xneuthl(i,j)
         enddo
      enddo

      do i=1,2
         NS_neuttop=NS_neuttop+xneutastop(i)+ xneutza(i)
            do j=1,2
              NS_neuttop=NS_neuttop+xneuta(i,j)
            enddo
      enddo
      end
c ==================================================================== c
c =========================  neutralino b bbar ======================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_neutbot(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i,J
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION abot(2,5),bbot(2,5)
      DOUBLE PRECISION opl(2,2),opr(2,2),oppl(5,5),oppr(5,5)
      DOUBLE PRECISION sgn(5)
      DOUBLE PRECISION dsbo(2),dsbob(2)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION hchichi(3,5,5),achichi(2,5,5)
      DOUBLE PRECISION hchachaR(3,2,2),hchachaL(3,2,2),achachaR(2,2,2),
     .     achachaL(2,2,2)
      DOUBLE PRECISION Hbbr(3),Abbr(2)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION X1,x3,y3,x2,xmuneut1,xmusb(2),xmub,
     .uh,th,vzz,azz,xneutsbot,
     .xmuz,dz,xneutzbot,rh,sh,rk,xneutzsbot
      DOUBLE PRECISION xneuthl(3,3),dhl(3),dhh(3),xneuta(2,2),
     .da(2),dna(2),xneuthlsbot(3),xneutza(2),xneutasbot(2),gmsb(2)
*
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup4/opl,opr,oppl,oppr
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_coupNMSSM/hchichi,achichi,hchachaR,hchachaL,achachaR,
     .     achachaL
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/NS_neutsbotbot/abot,bbot
      COMMON/NS_phibotbot/Hbbr,Abbr
*
      do i=1,5
         sgn(i) = 1d0
         if(xmneut(i).ge.0d0) then
            sgn(i) = 1d0
         elseif(xmneut(i).lt.0d0) then
            sgn(i) = -1d0
         endif
      enddo

      xmuneut1 = amneut(nj)**2/amneut(ni)**2
      xmub     = amb**2/amneut(ni)**2

      x3 = 2d0-x1-x2
      y3 = 1d0+xmuneut1-x3

      uh = 1d0-x1+xmub
      th = 1d0-x2+xmub

      vzz = vzzbotbot
      azz = azzbotbot
c -------------------------------------------------------------------- c
c 	                    sbottom exchange
c -------------------------------------------------------------------- c
      gmsb(1)=asb1
      gmsb(2)=asb2
      xmusb(1) = asb1**2/amneut(ni)**2
      xmusb(2) = asb2**2/amneut(ni)**2

      dsbo(1)  = 1d0-x1-xmusb(1)+xmub
      dsbo(2)  = 1d0-x1-xmusb(2)+xmub
      dsbob(1) = 1d0-x2-xmusb(1)+xmub
      dsbob(2) = 1d0-x2-xmusb(2)+xmub

      xneutsbot=0d0

      if((amneut(nj)+2d0*amb).le.amneut(ni)) then
         do i=1,2
            do k=1,2
      if (xmusb(k).gt.1d0.and.xmusb(i).gt.1d0)then
               xneutsbot=xneutsbot
     .          +g2s**2/dsbo(k)/dsbo(i)*(
     .           (abot(i,ni)*bbot(k,ni)+bbot(i,ni)*abot(k,ni))*
     .           (abot(i,nj)*bbot(k,nj)+bbot(i,nj)*abot(k,nj))*
     .           xmneut(nj)/xmneut(ni)*xmub*(-4d0)+
     .           (abot(i,ni)*abot(k,ni)+bbot(i,ni)*bbot(k,ni))*
     .           (abot(i,nj)*bbot(k,nj)+bbot(i,nj)*abot(k,nj))*
     .           xmneut(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*2d0*
     .           (uh-xmub-1d0)+
     .           (abot(i,ni)*bbot(k,ni)+bbot(i,ni)*abot(k,ni))*
     .           (abot(i,nj)*abot(k,nj)+bbot(i,nj)*bbot(k,nj))*
     .           dsqrt(xmub)*sgn(ni)*2d0*(uh-xmub-xmuneut1)+
     .           (abot(i,ni)*abot(k,ni)+bbot(i,ni)*bbot(k,ni))*
     .           (abot(i,nj)*abot(k,nj)+bbot(i,nj)*bbot(k,nj))*
     .           (-uh**2+uh*(1d0+xmuneut1+2d0*xmub)-(xmuneut1+xmub)*
     .           (1d0+xmub)))
     .           +g2s**2/dsbob(k)/dsbob(i)*(
     .           (abot(i,ni)*bbot(k,ni)+bbot(i,ni)*abot(k,ni))*
     .           (abot(i,nj)*bbot(k,nj)+bbot(i,nj)*abot(k,nj))*
     .           xmneut(nj)/xmneut(ni)*xmub*(-4d0)+
     .           (abot(i,ni)*abot(k,ni)+bbot(i,ni)*bbot(k,ni))*
     .           (abot(i,nj)*bbot(k,nj)+bbot(i,nj)*abot(k,nj))*
     .           xmneut(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*2d0*
     .           (th-xmub-1d0)+
     .           (abot(i,ni)*bbot(k,ni)+bbot(i,ni)*abot(k,ni))*
     .           (abot(i,nj)*abot(k,nj)+bbot(i,nj)*bbot(k,nj))*
     .           dsqrt(xmub)*sgn(ni)*2d0*(th-xmub-xmuneut1)+
     .           (abot(i,ni)*abot(k,ni)+bbot(i,ni)*bbot(k,ni))*
     .           (abot(i,nj)*abot(k,nj)+bbot(i,nj)*bbot(k,nj))*
     .           (-th**2+th*(1d0+xmuneut1+2d0*xmub)-(xmuneut1+xmub)*
     .           (1d0+xmub)))
     .           -2d0*g2s**2/dsbo(k)/dsbob(i)*(
     .           (bbot(i,ni)*bbot(k,ni)*abot(i,nj)*abot(k,nj)
     .           +abot(i,ni)*abot(k,ni)*bbot(i,nj)*bbot(k,nj))*
     .           xmneut(nj)/xmneut(ni)*xmub*(-2d0)+
     .           (abot(i,ni)*bbot(k,ni)*abot(i,nj)*abot(k,nj)
     .           +abot(k,ni)*bbot(i,ni)*bbot(i,nj)*bbot(k,nj))*
     .           xmneut(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*
     .           (th-xmub-1d0)+
     .           (abot(k,ni)*bbot(i,ni)*abot(i,nj)*abot(k,nj)
     .           +abot(i,ni)*bbot(k,ni)*bbot(i,nj)*bbot(k,nj))*
     .           xmneut(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*
     .           (uh-xmub-1d0)+
     .           (abot(k,ni)*abot(i,ni)*abot(i,nj)*abot(k,nj)
     .           +bbot(i,ni)*bbot(k,ni)*bbot(i,nj)*bbot(k,nj))*
     .           xmneut(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0)+
     .           (abot(k,ni)*bbot(i,ni)*abot(k,nj)*bbot(i,nj)
     .           +bbot(k,ni)*abot(i,ni)*abot(i,nj)*bbot(k,nj))*
     .           xmub*(uh+th-2d0*xmub)+
     .           (abot(k,ni)*abot(i,ni)*abot(k,nj)*bbot(i,nj)
     .           +bbot(k,ni)*bbot(i,ni)*abot(i,nj)*bbot(k,nj))*
     .           dsqrt(xmub)*sgn(ni)*(uh-xmub-xmuneut1)+
     .           (bbot(k,ni)*bbot(i,ni)*abot(k,nj)*bbot(i,nj)
     .           +abot(k,ni)*abot(i,ni)*abot(i,nj)*bbot(k,nj))*
     .           dsqrt(xmub)*sgn(ni)*(th-xmub-xmuneut1)+
     .           (bbot(k,ni)*abot(i,ni)*abot(k,nj)*bbot(i,nj)
     .           +abot(k,ni)*bbot(i,ni)*abot(i,nj)*bbot(k,nj))*
     .           (uh*th-xmub**2-xmuneut1))
      endif
            enddo
         enddo
      else
         xneutsbot=0d0
      endif

c -------------------------------------------------------------------- c
c 			        Z exchange
c -------------------------------------------------------------------- c
      xmuz = mz**2/amneut(ni)**2
      dz   = y3-xmuz

      xneutzbot=0d0

      rh = xmuneut1+2d0*xmub-th-uh+1d0
      sh = (xmuneut1-th-uh+1d0)*2d0*xmub+4d0*xmub**2
      rk = xmuneut1*(2d0*xmub-th-uh+4d0)+2d0*xmub-uh-th

      if((amneut(nj)+2d0*amb).le.amneut(ni)) then
         if ((mz+amneut(nj)).gt.amneut(ni))then
         xneutzbot=xneutzbot+g2s**2/dz**2*(
     .    oppl(ni,nj)*oppr(ni,nj)*(vzz**2-azz**2)*
     .    xmneut(nj)/xmneut(ni)*xmub*(-16d0/xmuz**2*rh**2+
     .    32d0/xmuz*rh-64d0)+
     .    oppl(ni,nj)*oppr(ni,nj)*(vzz**2+azz**2)*
     .    xmneut(nj)/xmneut(ni)*(8d0/xmuz**2*rh*sh-16d0/xmuz*sh
     .    -16d0*(xmuneut1-uh-th+1d0))+
     .    (oppl(ni,nj)**2+oppr(ni,nj)**2)*(vzz**2-azz**2)*
     .    xmub*(4d0/xmuz**2*rh*rk-8d0/xmuz*rk+8d0*(uh+th-2d0*xmub))
     .    +(oppl(ni,nj)**2+oppr(ni,nj)**2)*(vzz**2+azz**2)*
     .    (-2d0/xmuz**2*rk*sh+8d0/xmuz*(xmuneut1*(2d0*xmub**2+
     .    4d0*xmub-xmub*(th+uh))+2d0*xmub**2-xmub*(uh+th))+4d0*(
     .    xmuneut1*(uh+th-2d0*xmub-2d0)+2d0*xmub*(uh+th-1d0)
     .    -2d0*xmub**2+th*(-th+1d0)+uh*(-uh+1d0)))+
     .    (oppl(ni,nj)**2-oppr(ni,nj)**2)*vzz*azz*8d0*(
     .    xmuneut1*(th-uh)+2d0*xmub*(th-uh)+th*(-th+1d0)+uh*(uh-1d0))
     .    )
         endif
      else
         xneutzbot=0d0
      endif
c -------------------------------------------------------------------- c
c      NMSSM CP EVEN Higgs exchange + Interference H(i)-H(j)
c -------------------------------------------------------------------- c
        do i=1,3
        do j=1,3

      dhl(i)   = y3-smass(i)**2/amneut(ni)**2
      dhh(j)   = y3-smass(j)**2/amneut(ni)**2

      xneuthl(i,j)=0d0

      if((amneut(nj)+2d0*amb).le.amneut(ni)) then
         if((SMASS(i)+amneut(nj)).gt.amneut(ni).and.
     .(SMASS(j)+amneut(nj)).gt.amneut(ni))then
         xneuthl(i,j)=g2s**2/dhl(i)/dhh(j)*Hbbr(i)*Hbbr(j)
     .    *hchichi(i,ni,nj)*hchichi(j,ni,nj)*(
     .    xmneut(nj)/xmneut(ni)*xmub*(-32d0)+
     .    xmneut(nj)/xmneut(ni)*16d0*(1d0+xmuneut1-th-uh)+
     .    xmub*16d0*(2d0*xmub-th-uh)+
     .    8d0*(xmuneut1*(uh+th-2d0*xmub)+2d0*xmub*(uh+th-1d0)
     .    -(th+uh)**2+uh+th))
      endif
      else
         xneuthl(i,j)=0d0
      endif

      enddo
      enddo
c -------------------------------------------------------------------- c
c       NMSSM CP ODD Higgs exchange + Interference A(i)-A(j)
c -------------------------------------------------------------------- c
       do i=1,2
       do j=1,2

      da(i)    = y3-pmass(i)**2/amneut(ni)**2
      dna(j)   = y3-pmass(j)**2/amneut(ni)**2

      xneuta(i,j)=0d0

      if((amneut(nj)+2d0*amb).le.amneut(ni)) then
       if ((PMASS(i)+amneut(nj)).gt.amneut(ni).and.
     .(PMASS(j)+amneut(nj)).gt.amneut(ni))then
         xneuta(i,j)=g2s**2/da(i)/dna(j)*Abbr(i)*Abbr(j)
     .    *achichi(i,ni,nj)*achichi(j,ni,nj)*(
     .    xmneut(nj)/xmneut(ni)*xmub*(-32d0)+
     .    xmneut(nj)/xmneut(ni)*(-16d0)*(1d0+xmuneut1-th-uh)+
     .    xmub*(-16d0)*(2d0*xmub-th-uh)+
     .    8d0*(xmuneut1*(uh+th-2d0*xmub)+2d0*xmub*(uh+th-1d0)
     .    -(th+uh)**2+uh+th))
      endif
      else
         xneuta(i,j)=0d0
      endif

      enddo
      enddo

c -------------------------------------------------------------------- c
c 	                 interference Z-sbottom
c -------------------------------------------------------------------- c
      xneutzsbot=0d0

      if((amneut(nj)+2d0*amb).le.amneut(ni)) then
         do i=1,2
      if (xmusb(i).gt.1d0.and.(mz+amneut(nj)).gt.amneut(ni))then
            xneutzsbot=xneutzsbot
     .      +g2s**2/dsbo(i)/dz*oppl(ni,nj)*(
     .      ((abot(i,ni)*abot(i,nj)-bbot(i,ni)*bbot(i,nj))*vzz-
     .       (abot(i,ni)*abot(i,nj)+bbot(i,ni)*bbot(i,nj))*azz)*
     .      xmneut(nj)/xmneut(ni)*xmub*(
     .      1/xmuz*(-4d0)*(1d0+xmuneut1+2d0*xmub-uh-th)+16d0) +
     .      ((abot(i,ni)*abot(i,nj)-bbot(i,ni)*bbot(i,nj))*vzz+
     .       (abot(i,ni)*abot(i,nj)+bbot(i,ni)*bbot(i,nj))*azz)*
     .      xmneut(nj)/xmneut(ni)*(2d0/xmuz*((xmuneut1+1d0-uh-th)*
     .      2d0*xmub+4d0*xmub**2)+4d0*(1d0+xmuneut1-uh-th)) +
     .      ((bbot(i,ni)*bbot(i,nj)-abot(i,ni)*abot(i,nj))*vzz+
     .       (bbot(i,ni)*bbot(i,nj)+abot(i,ni)*abot(i,nj))*azz)*
     .      xmub*(2d0/xmuz*(xmuneut1*(2d0*xmub-th-uh+4d0)+2d0*xmub
     .      -th-uh)+4d0*(2d0*xmub-th-uh)) +
     .      ((bbot(i,ni)*bbot(i,nj)-abot(i,ni)*abot(i,nj))*vzz-
     .       (bbot(i,ni)*bbot(i,nj)+abot(i,ni)*abot(i,nj))*azz)*
     .      (2d0/xmuz*(xmuneut1*(-2d0*xmub**2+xmub*th-2d0*xmub+xmub*
     .       uh-2d0*xmub)+xmub*(-2d0*xmub+uh)+xmub*th)+4d0*(
     .      xmuneut1*(xmub-uh+1d0)+xmub*(xmub-uh)+xmub*(1d0-uh)+uh**2
     .      -uh)) +
     .      (abot(i,nj)*bbot(i,ni)*(vzz+azz)-
     .       abot(i,ni)*bbot(i,nj)*(vzz-azz))*
     .      xmneut(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*
     .      (2d0/xmuz*(xmuneut1*(1d0+xmub-uh)+1d0+xmub*(xmub-uh)+
     .      xmub*(xmub-th-2d0*uh+3d0)+th*(uh-1d0)+uh*(uh-2d0))+
     .      4d0*(1d0+xmub-th)) +
     .      (abot(i,nj)*bbot(i,ni)*(vzz-azz)-
     .       abot(i,ni)*bbot(i,nj)*(vzz+azz))*
     .      xmneut(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*(2d0/xmuz*
     .      (xmuneut1*(-1d0-xmub+uh)-1d0+xmub*(-1d0+uh)+uh*(2d0-th
     .      -uh)+th+xmub*(-2d0*xmub+th+2d0*uh-2d0))+8d0*(1d0+xmub-
     .      uh)) +
     .      (abot(i,ni)*bbot(i,nj)*(vzz-azz)-
     .       abot(i,nj)*bbot(i,ni)*(vzz+azz))*dsqrt(xmub)*sgn(ni)*(
     .      (-2d0)/xmuz*(xmuneut1-uh+xmub)*(xmuneut1+2d0*xmub-th-uh
     .      +1d0)+8d0*(xmuneut1+xmub-uh)) +
     .      (abot(i,ni)*bbot(i,nj)*(vzz+azz)-
     .       abot(i,nj)*bbot(i,ni)*(vzz-azz))*dsqrt(xmub)*sgn(ni)*(
     .      2d0/xmuz*(xmuneut1*(xmuneut1+3d0*xmub-th-2d0*uh+1d0)+
     .      xmub*(2d0*xmub-th-2d0*uh)+xmub*(1d0-uh)+uh*th+uh**2-uh)+
     .      4d0*(xmuneut1+xmub-th)) )
     .      -g2s**2/dsbob(i)/dz*oppl(ni,nj)*(
     .      ((-abot(i,ni)*abot(i,nj)+bbot(i,ni)*bbot(i,nj))*vzz+
     .       (abot(i,ni)*abot(i,nj)+bbot(i,ni)*bbot(i,nj))*azz)*
     .      xmneut(nj)/xmneut(ni)*xmub*(
     .      1/xmuz*(-4d0)*(1d0+xmuneut1+2d0*xmub-uh-th)+16d0) +
     .      ((-abot(i,ni)*abot(i,nj)+bbot(i,ni)*bbot(i,nj))*vzz-
     .       (abot(i,ni)*abot(i,nj)+bbot(i,ni)*bbot(i,nj))*azz)*
     .      xmneut(nj)/xmneut(ni)*(2d0/xmuz*(2d0*xmub*(xmuneut1+1d0
     .      -th-uh)+4d0*xmub**2)+4d0*(1d0+xmuneut1-uh-th)) +
     .      ((abot(i,ni)*abot(i,nj)-bbot(i,ni)*bbot(i,nj))*vzz-
     .       (abot(i,ni)*abot(i,nj)+bbot(i,ni)*bbot(i,nj))*azz)*
     .      xmub*(2d0/xmuz*(xmuneut1*(2d0*xmub-th-uh+4d0)+2d0*xmub
     .      -th-uh)+4d0*(2d0*xmub-th-uh)) +
     .      ((abot(i,ni)*abot(i,nj)-bbot(i,ni)*bbot(i,nj))*vzz+
     .       (abot(i,ni)*abot(i,nj)+bbot(i,ni)*bbot(i,nj))*azz)*
     .      (2d0/xmuz*(xmuneut1*(-2d0*xmub**2+xmub*th-2d0*xmub+xmub*
     .       uh-2d0*xmub)+xmub*(-2d0*xmub+uh)+xmub*th)+4d0*(
     .      xmuneut1*(xmub-th+1d0)+xmub*(xmub-th)+xmub*(1d0-th)+th**2
     .      -th)) +
     .      (abot(i,nj)*bbot(i,ni)*(-vzz-azz)+
     .       abot(i,ni)*bbot(i,nj)*(vzz-azz))*
     .      xmneut(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*(2d0/xmuz*
     .      (xmuneut1*(1d0+xmub-th)+1d0+xmub*(2d0*xmub-2d0*th-uh+
     .      3d0)+th*(th-xmub+uh-2d0)-uh)+4d0*(1d0+xmub-uh)) +
     .      (abot(i,nj)*bbot(i,ni)*(-vzz+azz)+
     .       abot(i,ni)*bbot(i,nj)*(vzz+azz))*
     .      xmneut(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*(2d0/xmuz*
     .      (xmuneut1*(-1d0-xmub+th)-1d0+xmub*(-xmub+th-1d0)+xmub*
     .      (2d0*th+uh-2d0-xmub)-th*(th+uh)+2d0*th+uh)+8d0*(1d0
     .      +xmub-th)) +
     .      (abot(i,ni)*bbot(i,nj)*(-vzz+azz)+
     .       abot(i,nj)*bbot(i,ni)*(vzz+azz))*dsqrt(xmub)*sgn(ni)*(
     .      2d0/xmuz*(xmuneut1*(-xmuneut1-3d0*xmub+2d0*th+uh-1d0)
     .      +xmub*(-2d0*xmub+2d0*th+uh-1d0)+th*(xmub-th-uh+1))
     .      +8d0*(xmuneut1+xmub-th)) +
     .      (abot(i,ni)*bbot(i,nj)*(-vzz-azz)+
     .       abot(i,nj)*bbot(i,ni)*(vzz-azz))*dsqrt(xmub)*sgn(ni)*(
     .      2d0/xmuz*(xmuneut1*(xmuneut1+3d0*xmub-uh-2d0*th+1d0)+
     .      xmub*(xmub-th+1d0)+xmub*(xmub-2d0*th-uh)+uh*th+th**2-th)+
     .      4d0*(xmuneut1+xmub-uh)) )
      endif
         enddo
      else
         xneutzsbot=0d0
      endif

c -------------------------------------------------------------------- c
c 			interference H(j)-sbottom
c -------------------------------------------------------------------- c
      do j=1,3
         dhl(j)   = y3-SMASS(j)**2/amneut(ni)**2
      xneuthlsbot(j)=0d0


      if((amneut(nj)+2d0*amb).le.amneut(ni)) then
         do i=1,2
      if ((smass(j)+amneut(nj)).gt.amneut(ni).and.xmusb(i).gt.1d0)then
            xneuthlsbot(j)=xneuthlsbot(j)
     .       +2d0*g2s**2/dhl(j)/dsbo(i)*(Hbbr(j)/dsqrt(2d0))*
     .       (-2d0)*hchichi(j,ni,nj)*(
     .       (abot(i,nj)*abot(i,ni)+bbot(i,ni)*bbot(i,nj))*(
     .       dsqrt(xmub)*sgn(ni)*(xmuneut1+xmub-uh) +
     .       dsqrt(xmub)*sgn(ni)*(-xmuneut1-xmub+th) +
     .       xmneut(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*(-1d0-xmub+th)+
     .       xmneut(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*(1d0+xmub-uh) )
     .       +(abot(i,nj)*bbot(i,ni)+abot(i,ni)*bbot(i,nj))*(
     .       xmneut(nj)/xmneut(ni)*(-1d0-xmuneut1+uh+th) +
     .       2d0*xmneut(nj)/xmneut(ni)*xmub +
     .       xmub*(-2d0*xmub+th+uh) +
     .       (uh**2+uh*th-uh*(1d0+2d0*xmub+xmuneut1)+xmub+xmuneut1*
     .       xmub)) )
     .       +2d0*g2s**2/dhl(j)/dsbob(i)*(Hbbr(j)/dsqrt(2d0))*
     .       (-2d0)*hchichi(j,ni,nj)*(
     .       (abot(i,nj)*abot(i,ni)+bbot(i,ni)*bbot(i,nj))*(
     .       dsqrt(xmub)*sgn(ni)*(uh-xmuneut1-xmub) +
     .       dsqrt(xmub)*sgn(ni)*(-th+xmub+xmuneut1) +
     .       dsqrt(xmub)*sgn(ni)*xmneut(nj)/xmneut(ni)*(-th+1d0+xmub) +
     .       dsqrt(xmub)*sgn(ni)*xmneut(nj)/xmneut(ni)*(uh-1d0-xmub)) +
     .       (abot(i,nj)*bbot(i,ni)+abot(i,ni)*bbot(i,nj))*(
     .       2d0*xmub*xmneut(nj)/xmneut(ni) +
     .       (uh*th+th**2-th*(1d0+2d0*xmub+xmuneut1)+xmub+
     .       xmub*xmuneut1) +
     .       xmneut(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0) +
     .       xmub*(uh+th-2d0*xmub)) )
      endif
         enddo
      else
         xneuthlsbot(j)=0d0
      endif


      enddo
c -------------------------------------------------------------------- c
c                        interference A(j)-sbottom
c -------------------------------------------------------------------- c
      do j=1,2
      da(j)    = y3-pmass(j)**2/amneut(ni)**2
      xneutasbot(j)=0d0


      if((amneut(nj)+2d0*amb).le.amneut(ni)) then
         do i=1,2
       if ((pmass(j)+amneut(nj)).gt.amneut(ni).and.
     .xmusb(i).gt.1d0)then
            xneutasbot(j)=xneutasbot(j)
     .       +2d0*g2s**2/da(j)/dsbo(i)*(Abbr(j)/dsqrt(2d0))*
     .       2d0*achichi(j,ni,nj)*(
     .       (abot(i,nj)*abot(i,ni)+bbot(i,ni)*bbot(i,nj))*(
     .       dsqrt(xmub)*sgn(ni)*(xmuneut1+xmub-uh) +
     .       dsqrt(xmub)*sgn(ni)*(-xmuneut1-xmub+th)*(-1d0) +
     .       xmneut(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*(-1d0-xmub+th)*
     .       (-1d0) +
     .       xmneut(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*(1d0+xmub-uh) )
     .       +(abot(i,nj)*bbot(i,ni)+abot(i,ni)*bbot(i,nj))*(
     .       xmneut(nj)/xmneut(ni)*(-1d0-xmuneut1+uh+th)*(-1d0) +
     .       2d0*xmneut(nj)/xmneut(ni)*xmub +
     .       xmub*(-2d0*xmub+th+uh)*(-1d0) +
     .       (uh**2+uh*th-uh*(1d0+2d0*xmub+xmuneut1)+xmub+xmuneut1*
     .       xmub)) )
     .       +2d0*g2s**2/da(j)/dsbob(i)*(Abbr(j)/dsqrt(2d0))*
     .       2d0*achichi(j,ni,nj)*(
     .       (abot(i,nj)*abot(i,ni)+bbot(i,ni)*bbot(i,nj))*(
     .       dsqrt(xmub)*sgn(ni)*(uh-xmuneut1-xmub)*(-1d0) +
     .       dsqrt(xmub)*sgn(ni)*(-th+xmub+xmuneut1) +
     .       dsqrt(xmub)*sgn(ni)*xmneut(nj)/xmneut(ni)*(-th+1d0+xmub) +
     .       dsqrt(xmub)*sgn(ni)*xmneut(nj)/xmneut(ni)*(uh-1d0-xmub)*
     .       (-1d0)) +
     .       (abot(i,nj)*bbot(i,ni)+abot(i,ni)*bbot(i,nj))*(
     .       2d0*xmub*xmneut(nj)/xmneut(ni) +
     .       (uh*th+th**2-th*(1d0+2d0*xmub+xmuneut1)+xmub+
     .       xmub*xmuneut1) +
     .       xmneut(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0)*(-1d0) +
     .       xmub*(uh+th-2d0*xmub)*(-1d0)) )
      endif
         enddo
      else
         xneutasbot(j)=0d0

      endif
      enddo

c -------------------------------------------------------------------- c
c 	                interference Z and A(i)
C -----------------INTERFERENCE TERMS Z-h/H ARE MISSING
c -------------------------------------------------------------------- c
      do i=1,2
      da(i)= y3-pmass(i)**2/amneut(ni)**2

      xneutza(i)=0d0

      if((amneut(nj)+2d0*amb).le.amneut(ni)) then
      if ((pmass(i)+amneut(nj)).gt.amneut(ni).and.(mz+amneut(nj))
     ..gt.amneut(ni))then
         xneutza(i)=-4d0*g2s**2/da(i)/dz*azz*
     .    Abbr(i)/dsqrt(2d0)*2d0*achichi(i,ni,nj)*oppl(ni,nj)*(
     .    xmneut(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*(
     .    4d0/xmuz*(2d0+xmuneut1*(2d0*xmub-th-uh+2d0)+xmub*(-uh-th
     .    +1d0+xmub)+xmub*(3d0*xmub-3d0*(th+uh)+5d0)+(th+uh)**2
     .    -3d0*(th+uh))+
     .    4d0*(th-xmub-1d0)+
     .    4d0*(uh-xmub-1d0))
     .    -dsqrt(xmub)*sgn(ni)*(
     .    4d0/xmuz*(-2d0*xmuneut1**2+xmuneut1*(-6d0*xmub+
     .    3d0*(th+uh)-2d0)+xmub*(-xmub+th+uh-1d0)+xmub*(-3d0*xmub
     .    +3d0*(uh+th)-1d0)-(th+uh)**2+(th+uh))+
     .    4d0*(xmuneut1+xmub-uh)+
     .    4d0*(xmuneut1+xmub-th)) )
      endif
      else
         xneutza(i)=0d0
      endif

      enddo

      NS_neutbot=xneutzbot+xneutsbot+xneutzsbot

      do i=1,3
         NS_neutbot= NS_neutbot+xneuthlsbot(i)
         do j=1,3
            NS_neutbot=NS_neutbot+xneuthl(i,j)
         enddo
      enddo

      do i=1,2
         NS_neutbot=NS_neutbot+xneutasbot(i)+xneutza(i)
         do j=1,2
           NS_neutbot=NS_neutbot+xneuta(i,j)
         enddo
      enddo
      end
c ==================================================================== c
c =========================  neutralino e+ e- ======================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_neutel(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION opl(2,2),opr(2,2),oppl(5,5),oppr(5,5)
      DOUBLE PRECISION dsel(2),dselb(2)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      DOUBLE PRECISION X1,x3,y3,y1,y2,x2,xmuneut1,xmusel(2),
     .xneutsel,xmuz,dz,xneutzel,xneutzsel
*
      COMMON/NS_indices/ni,nj
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup4/opl,opr,oppl,oppr
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS

      xmuneut1 = amneut(nj)**2/amneut(ni)**2

      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                          selectron exchange
c -------------------------------------------------------------------- c
      xmusel(1) = ase1**2/amneut(ni)**2
      xmusel(2) = ase2**2/amneut(ni)**2

      dsel(1)  = 1-x1-xmusel(1)
      dsel(2)  = 1-x1-xmusel(2)
      dselb(1) = 1-x2-xmusel(1)
      dselb(2) = 1-x2-xmusel(2)

      xneutsel=0d0

      if (amneut(nj).le.amneut(ni)) then
         do i=1,2
            do k=1,2
      if (xmusel(i).gt.1d0.and.xmusel(k).gt.1d0)then
               xneutsel=xneutsel
     .              +g2s**2/dsel(k)/dsel(i)*x1*y1*
     .               (ae(i,ni)*ae(k,ni)+be(i,ni)*be(k,ni))*
     .               (ae(i,nj)*ae(k,nj)+be(i,nj)*be(k,nj))
     .              +g2s**2/dselb(k)/dselb(i)*x2*y2*
     .               (ae(i,ni)*ae(k,ni)+be(i,ni)*be(k,ni))*
     .               (ae(i,nj)*ae(k,nj)+be(i,nj)*be(k,nj))
     .              +g2s**2/dselb(k)/dsel(i)*
     .               ( ( ae(i,ni)*be(k,ni)*ae(k,nj)*be(i,nj)
     .                  +ae(i,nj)*be(k,nj)*ae(k,ni)*be(i,ni))
     .                 *(-x1*y1-x2*y2+x3*y3)
     .                +2d0*xmneut(nj)/xmneut(ni)*y3*
     .                 ( ae(i,ni)*ae(k,ni)*ae(k,nj)*ae(i,nj)
     .                  +be(i,nj)*be(k,nj)*be(k,ni)*be(i,ni)))
      endif
            enddo
         enddo
      else
         xneutsel=0d0
      endif

c -------------------------------------------------------------------- c
c                              Z exchange
c -------------------------------------------------------------------- c
      xmuz = mz**2/amneut(ni)**2
      dz   = y3-xmuz

      xneutzel=0d0

      if (amneut(nj).le.amneut(ni)) then
      if ((mz+amneut(nj)).gt.amneut(ni))then
         xneutzel=g2s**2*4d0/dz**2*
     .        (((azztautau+vzztautau)**2*oppl(ni,nj)**2
     .        + (azztautau-vzztautau)**2*oppr(ni,nj)**2)*x2*y2
     .        +((azztautau+vzztautau)**2*oppr(ni,nj)**2
     .        + (azztautau-vzztautau)**2*oppl(ni,nj)**2)*x1*y1
     .        -4d0*xmneut(nj)/xmneut(ni)*oppl(ni,nj)*oppr(ni,nj)
     .        *(azztautau**2+vzztautau**2)*y3 )
      endif
      else
         xneutzel=0d0
      endif
c -------------------------------------------------------------------- c
c                        Z-selectron interference
c -------------------------------------------------------------------- c
      xneutzsel=0d0

      if (amneut(nj).le.amneut(ni)) then
         do i=1,2
      if (xmusel(i).gt.1d0.and.(mz+amneut(nj)).gt.amneut(ni))then
           xneutzsel=xneutzsel-g2s**2*4d0/dsel(i)/dz
     .      *((ae(i,ni)*ae(i,nj)*oppr(ni,nj)*(azztautau+vzztautau)+
     .         be(i,ni)*be(i,nj)*oppl(ni,nj)*
     .        (-azztautau+vzztautau))*x1*y1
     .       -(ae(i,nj)*ae(i,ni)*oppl(ni,nj)*(azztautau+vzztautau)
     .        +be(i,nj)*be(i,ni)*oppr(ni,nj)*(-azztautau+vzztautau)
     .          )*xmneut(nj)/xmneut(ni)*y3)
     .        +g2s**2*4d0/dselb(i)/dz
     .      *((ae(i,ni)*ae(i,nj)*oppl(ni,nj)*(azztautau+vzztautau)+
     .         be(i,ni)*be(i,nj)*oppr(ni,nj)*
     .         (-azztautau+vzztautau))*x2*y2
     .        -(ae(i,nj)*ae(i,ni)*oppr(ni,nj)*(azztautau+vzztautau)
     .         +be(i,nj)*be(i,ni)*oppl(ni,nj)*
     .          (-azztautau+vzztautau)
     .          )*xmneut(nj)/xmneut(ni)*y3)
      endif
         enddo
      else
         xneutzsel=0d0
      endif


      NS_neutel = xneutsel+xneutzel+xneutzsel

      end
c ==================================================================== c
c =========================  neutralino mu+ mu- ====================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_neutmu(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION opl(2,2),opr(2,2),oppl(5,5),oppr(5,5)
      DOUBLE PRECISION dsel(2),dselb(2)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      DOUBLE PRECISION X1,x3,y3,y1,y2,x2,xmuneut1,xmusel(2),
     .xneutsel,xmuz,dz,xneutzel,xneutzsel
*
      COMMON/NS_indices/ni,nj
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup4/opl,opr,oppl,oppr
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS

      xmuneut1 = amneut(nj)**2/amneut(ni)**2

      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                          smon exchange
c -------------------------------------------------------------------- c
      xmusel(1) = asmu1**2/amneut(ni)**2
      xmusel(2) = asmu2**2/amneut(ni)**2

      dsel(1)  = 1-x1-xmusel(1)
      dsel(2)  = 1-x1-xmusel(2)
      dselb(1) = 1-x2-xmusel(1)
      dselb(2) = 1-x2-xmusel(2)

      xneutsel=0d0

      if (amneut(nj).le.amneut(ni)) then
         do i=1,2
            do k=1,2
      if (xmusel(i).gt.1d0.and.xmusel(k).gt.1d0)then
               xneutsel=xneutsel
     .              +g2s**2/dsel(k)/dsel(i)*x1*y1*
     .               (amu(i,ni)*amu(k,ni)+bmu(i,ni)*bmu(k,ni))*
     .               (amu(i,nj)*amu(k,nj)+bmu(i,nj)*bmu(k,nj))
     .              +g2s**2/dselb(k)/dselb(i)*x2*y2*
     .               (amu(i,ni)*amu(k,ni)+bmu(i,ni)*bmu(k,ni))*
     .               (amu(i,nj)*amu(k,nj)+bmu(i,nj)*bmu(k,nj))
     .              +g2s**2/dselb(k)/dsel(i)*
     .               ( ( amu(i,ni)*bmu(k,ni)*amu(k,nj)*bmu(i,nj)
     .                  +amu(i,nj)*bmu(k,nj)*amu(k,ni)*bmu(i,ni))
     .                 *(-x1*y1-x2*y2+x3*y3)
     .                +2d0*xmneut(nj)/xmneut(ni)*y3*
     .                 ( amu(i,ni)*amu(k,ni)*amu(k,nj)*amu(i,nj)
     .                  +bmu(i,nj)*bmu(k,nj)*bmu(k,ni)*bmu(i,ni)))
      endif
            enddo
         enddo
      else
         xneutsel=0d0
      endif

c -------------------------------------------------------------------- c
c                              Z exchange
c -------------------------------------------------------------------- c
      xmuz = mz**2/amneut(ni)**2
      dz   = y3-xmuz

      xneutzel=0d0

      if (amneut(nj).le.amneut(ni)) then
      if ((mz+amneut(nj)).gt.amneut(ni))then
         xneutzel=g2s**2*4d0/dz**2*
     .        (((azztautau+vzztautau)**2*oppl(ni,nj)**2
     .        + (azztautau-vzztautau)**2*oppr(ni,nj)**2)*x2*y2
     .        +((azztautau+vzztautau)**2*oppr(ni,nj)**2
     .        + (azztautau-vzztautau)**2*oppl(ni,nj)**2)*x1*y1
     .        -4d0*xmneut(nj)/xmneut(ni)*oppl(ni,nj)*oppr(ni,nj)
     .        *(azztautau**2+vzztautau**2)*y3 )
      endif
      else
         xneutzel=0d0
      endif
c -------------------------------------------------------------------- c
c                        Z-smuon interference
c -------------------------------------------------------------------- c
      xneutzsel=0d0

      if (amneut(nj).le.amneut(ni)) then
         do i=1,2
      if (xmusel(i).gt.1d0.and.(mz+amneut(nj)).gt.amneut(ni))then
           xneutzsel=xneutzsel-g2s**2*4d0/dsel(i)/dz
     .      *((amu(i,ni)*amu(i,nj)*oppr(ni,nj)*(azztautau+vzztautau)+
     .         bmu(i,ni)*bmu(i,nj)*oppl(ni,nj)*
     .        (-azztautau+vzztautau))*x1*y1
     .       -(amu(i,nj)*amu(i,ni)*oppl(ni,nj)*(azztautau+vzztautau)
     .        +bmu(i,nj)*bmu(i,ni)*oppr(ni,nj)*(-azztautau+vzztautau)
     .          )*xmneut(nj)/xmneut(ni)*y3)
     .        +g2s**2*4d0/dselb(i)/dz
     .      *((amu(i,ni)*amu(i,nj)*oppl(ni,nj)*(azztautau+vzztautau)+
     .         bmu(i,ni)*bmu(i,nj)*oppr(ni,nj)*
     .         (-azztautau+vzztautau))*x2*y2
     .        -(amu(i,nj)*amu(i,ni)*oppr(ni,nj)*(azztautau+vzztautau)
     .         +bmu(i,nj)*bmu(i,ni)*oppl(ni,nj)*
     .          (-azztautau+vzztautau)
     .          )*xmneut(nj)/xmneut(ni)*y3)
      endif
         enddo
      else
         xneutzsel=0d0
      endif


      NS_neutmu = xneutsel+xneutzel+xneutzsel

      end
c ==================================================================== c
c =======================  neutralino tau+ tau- ====================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_neuttau(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i,j
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION opl(2,2),opr(2,2),oppl(5,5),oppr(5,5)
      DOUBLE PRECISION dsto(2),dstob(2)
      DOUBLE PRECISION sgn(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION hchichi(3,5,5),achichi(2,5,5),
     . hchachaR(3,2,2),hchachaL(3,2,2),achachaR(2,2,2),
     .     achachaL(2,2,2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),
     .  CMASS,P(2,3)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION vzz,azz,xmuneut1,xmustau(2),
     .xmutau,x1,x2,x3,y3,uh,th,xneutstau,xmuz,dz,xneutztau,rh,sh,rk,
     .xneutzstau
      DOUBLE PRECISION xneuthl(3,3),dhl(3),dhh(3),da(2),dna(2),
     .xneuta(2,2),xneuthlstau(3),xneutastau(2),xneutza(2)
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
*
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup4/opl,opr,oppl,oppr
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coupNMSSM/hchichi,achichi,hchachaR,hchachaL,achachaR,
     .     achachaL
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_CPodd_MIX/P

      do i=1,5
         sgn(i) = 1d0
         if(xmneut(i).ge.0d0) then
            sgn(i) = 1d0
         elseif(xmneut(i).lt.0d0) then
            sgn(i) = -1d0
         endif
      enddo

      vzz = vzztautau
      azz = azztautau

      xmuneut1 = amneut(nj)**2/amneut(ni)**2
      xmustau(1) = astau1**2/amneut(ni)**2
      xmustau(2) = astau2**2/amneut(ni)**2
      xmutau   = amtau**2/amneut(ni)**2

      x3 = 2d0-x1-x2
      y3 = 1d0+xmuneut1-x3

      Uh = 1d0-x1+xmutau
      th = 1d0-x2+xmutau
c -------------------------------------------------------------------- c
c                              stau exchange
c -------------------------------------------------------------------- c
      dsto(1)  = 1d0-x1-xmustau(1)+xmutau
      dsto(2)  = 1d0-x1-xmustau(2)+xmutau
      dstob(1) = 1d0-x2-xmustau(1)+xmutau
      dstob(2) = 1d0-x2-xmustau(2)+xmutau

      xneutstau = 0d0

      if ((amneut(nj)+2d0*amtau).le.amneut(ni)) then
         do i=1,2
            do k=1,2
      if (xmustau(i).gt.1d0.and.xmustau(k).gt.1d0)then
               xneutstau=xneutstau
     .          +g2s**2/dsto(k)/dsto(i)*(
     .           (atau(i,ni)*btau(k,ni)+btau(i,ni)*atau(k,ni))*
     .           (atau(i,nj)*btau(k,nj)+btau(i,nj)*atau(k,nj))*
     .           xmneut(nj)/xmneut(ni)*xmutau*(-4d0)+
     .           (atau(i,ni)*atau(k,ni)+btau(i,ni)*btau(k,ni))*
     .           (atau(i,nj)*btau(k,nj)+btau(i,nj)*atau(k,nj))*
     .           xmneut(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*2d0*
     .           (uh-xmutau-1d0)+
     .           (atau(i,ni)*btau(k,ni)+btau(i,ni)*atau(k,ni))*
     .           (atau(i,nj)*atau(k,nj)+btau(i,nj)*btau(k,nj))*
     .           dsqrt(xmutau)*sgn(ni)*2d0*(uh-xmutau-xmuneut1)+
     .           (atau(i,ni)*atau(k,ni)+btau(i,ni)*btau(k,ni))*
     .           (atau(i,nj)*atau(k,nj)+btau(i,nj)*btau(k,nj))*
     .           (-uh**2+uh*(1d0+xmuneut1+2d0*xmutau)-
     .           (xmuneut1+xmutau)*(1d0+xmutau)))
     .           +g2s**2/dstob(k)/dstob(i)*(
     .           (atau(i,ni)*btau(k,ni)+btau(i,ni)*atau(k,ni))*
     .           (atau(i,nj)*btau(k,nj)+btau(i,nj)*atau(k,nj))*
     .           xmneut(nj)/xmneut(ni)*xmutau*(-4d0)+
     .           (atau(i,ni)*atau(k,ni)+btau(i,ni)*btau(k,ni))*
     .           (atau(i,nj)*btau(k,nj)+btau(i,nj)*atau(k,nj))*
     .           xmneut(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*2d0*
     .           (th-xmutau-1d0)+
     .           (atau(i,ni)*btau(k,ni)+btau(i,ni)*atau(k,ni))*
     .           (atau(i,nj)*atau(k,nj)+btau(i,nj)*btau(k,nj))*
     .           dsqrt(xmutau)*sgn(ni)*2d0*(th-xmutau-xmuneut1)+
     .           (atau(i,ni)*atau(k,ni)+btau(i,ni)*btau(k,ni))*
     .           (atau(i,nj)*atau(k,nj)+btau(i,nj)*btau(k,nj))*
     .           (-th**2+th*(1d0+xmuneut1+2d0*xmutau)-
     .           (xmuneut1+xmutau)*(1d0+xmutau)))
     .           -2d0*g2s**2/dsto(k)/dstob(i)*(
     .           (btau(i,ni)*btau(k,ni)*atau(i,nj)*atau(k,nj)
     .           +atau(i,ni)*atau(k,ni)*btau(i,nj)*btau(k,nj))*
     .           xmneut(nj)/xmneut(ni)*xmutau*(-2d0)+
     .           (atau(i,ni)*btau(k,ni)*atau(i,nj)*atau(k,nj)
     .           +atau(k,ni)*btau(i,ni)*btau(i,nj)*btau(k,nj))*
     .           xmneut(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*
     .           (th-xmutau-1d0)+
     .           (atau(k,ni)*btau(i,ni)*atau(i,nj)*atau(k,nj)
     .           +atau(i,ni)*btau(k,ni)*btau(i,nj)*btau(k,nj))*
     .           xmneut(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*
     .           (uh-xmutau-1d0)+
     .           (atau(k,ni)*atau(i,ni)*atau(i,nj)*atau(k,nj)
     .           +btau(i,ni)*btau(k,ni)*btau(i,nj)*btau(k,nj))*
     .           xmneut(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0)+
     .           (atau(k,ni)*btau(i,ni)*atau(k,nj)*btau(i,nj)
     .           +btau(k,ni)*atau(i,ni)*atau(i,nj)*btau(k,nj))*
     .           xmutau*(uh+th-2d0*xmutau)+
     .           (atau(k,ni)*atau(i,ni)*atau(k,nj)*btau(i,nj)
     .           +btau(k,ni)*btau(i,ni)*atau(i,nj)*btau(k,nj))*
     .           dsqrt(xmutau)*sgn(ni)*(uh-xmutau-xmuneut1)+
     .           (btau(k,ni)*btau(i,ni)*atau(k,nj)*btau(i,nj)
     .           +atau(k,ni)*atau(i,ni)*atau(i,nj)*btau(k,nj))*
     .           dsqrt(xmutau)*sgn(ni)*(th-xmutau-xmuneut1)+
     .           (btau(k,ni)*atau(i,ni)*atau(k,nj)*btau(i,nj)
     .           +atau(k,ni)*btau(i,ni)*atau(i,nj)*btau(k,nj))*
     .           (uh*th-xmutau**2-xmuneut1))
      endif
            enddo
         enddo
      else
         xneutstau=0d0
      endif

c -------------------------------------------------------------------- c
c 	                         Z exchange
c -------------------------------------------------------------------- c
      xmuz = mz**2/amneut(ni)**2
      dz   = y3-xmuz

      xneutztau=0d0

      rh = xmuneut1+2d0*xmutau-th-uh+1d0
      sh = (xmuneut1-th-uh+1d0)*2d0*xmutau+4d0*xmutau**2
      rk = xmuneut1*(2d0*xmutau-th-uh+4d0)+2d0*xmutau-uh-th

      if ((amneut(nj)+2d0*amtau).le.amneut(ni)) then
      if ((mz+amneut(nj)).gt.amneut(ni))then
         xneutztau=xneutztau+g2s**2/dz**2*(
     .    oppl(ni,nj)*oppr(ni,nj)*(vzz**2-azz**2)*
     .    xmneut(nj)/xmneut(ni)*xmutau*(-16d0/xmuz**2*rh**2+
     .    32d0/xmuz*rh-64d0)+
     .    oppl(ni,nj)*oppr(ni,nj)*(vzz**2+azz**2)*
     .    xmneut(nj)/xmneut(ni)*(8d0/xmuz**2*rh*sh-16d0/xmuz*sh
     .    -16d0*(xmuneut1-uh-th+1d0))+
     .    (oppl(ni,nj)**2+oppr(ni,nj)**2)*(vzz**2-azz**2)*
     .    xmutau*(4d0/xmuz**2*rh*rk-8d0/xmuz*rk+8d0*(uh+th-
     .    2d0*xmutau))
     .    +(oppl(ni,nj)**2+oppr(ni,nj)**2)*(vzz**2+azz**2)*
     .    (-2d0/xmuz**2*rk*sh+8d0/xmuz*(xmuneut1*(2d0*xmutau**2+
     .    4d0*xmutau-xmutau*(th+uh))+2d0*xmutau**2-xmutau*(uh+th))+
     .    4d0*(xmuneut1*(uh+th-2d0*xmutau-2d0)+2d0*xmutau*(uh+th-
     .    1d0)-2d0*xmutau**2+th*(-th+1d0)+uh*(-uh+1d0)))+
     .    (oppl(ni,nj)**2-oppr(ni,nj)**2)*vzz*azz*8d0*(
     .    xmuneut1*(th-uh)+2d0*xmutau*(th-uh)+th*(-th+1d0)+uh*(uh-
     .    1d0)))
      endif
      else
         xneutztau=0d0
      endif
c -------------------------------------------------------------------- c
c     NMSSM CP EVEN Higgs exchange + Interference H(i)-H(j)
c -------------------------------------------------------------------- c
        do i=1,3
        do j=1,3
      dhl(i)   = y3-smass(i)**2/amneut(ni)**2
      dhh(j)   = y3-smass(j)**2/amneut(ni)**2

      xneuthl(i,j)=0d0

      if((amneut(nj)+2d0*amtau).le.amneut(ni)) then
       if ((smass(i)+amneut(nj)).gt.amneut(ni).and.
     .(smass(j)+amneut(nj)).gt.amneut(ni))then
         xneuthl(i,j)=g2s**2/dhl(i)/dhh(j)*(scaltau*S(i,2))*
     .    (scaltau*S(j,2))*hchichi(i,ni,nj)*hchichi(j,ni,nj)*(
     .    xmneut(nj)/xmneut(ni)*xmutau*(-32d0)+
     .    xmneut(nj)/xmneut(ni)*16d0*(1d0+xmuneut1-th-uh)+
     .    xmutau*16d0*(2d0*xmutau-th-uh)+
     .    8d0*(xmuneut1*(uh+th-2d0*xmutau)+2d0*xmutau*(uh+th-1d0)
     .    -(th+uh)**2+uh+th))
       endif
      else
         xneuthl(i,j)=0d0
      endif

      enddo
      enddo
c -------------------------------------------------------------------- c
c     NMSSM CP ODD Higgs exchange + Interference A(i)-A(j)
c -------------------------------------------------------------------- c
         do i=1,2
         do j=1,2

      da(i)   = y3-pmass(i)**2/amneut(ni)**2
      dna(j)  = y3-pmass(j)**2/amneut(ni)**2

      xneuta(i,j)=0d0

      if((amneut(nj)+2d0*amtau).le.amneut(ni)) then
       if ((pmass(i)+amneut(nj)).gt.amneut(ni)
     ..and.(pmass(j)+amneut(nj)).gt.amneut(ni))then
      xneuta(i,j)=g2s**2/da(i)/dna(j)*(scaltau*P(i,2))*
     .    (scaltau*P(j,2))*achichi(i,ni,nj)*achichi(j,ni,nj)*(
     .    xmneut(nj)/xmneut(ni)*xmutau*(-32d0)+
     .    xmneut(nj)/xmneut(ni)*(-16d0)*(1d0+xmuneut1-th-uh)+
     .    xmutau*(-16d0)*(2d0*xmutau-th-uh)+
     .    8d0*(xmuneut1*(uh+th-2d0*xmutau)+2d0*xmutau*(uh+th-1d0)
     .    -(th+uh)**2+uh+th))
      endif
      else
         xneuta(i,j)=0d0
      endif

      enddo
      enddo

c -------------------------------------------------------------------- c
c                           Z-stau interference
c -------------------------------------------------------------------- c
      xneutzstau=0d0

      if((amneut(nj)+2d0*amtau).le.amneut(ni)) then
         do i=1,2
            if (xmustau(i).gt.1d0.and.(mz+amneut(nj)).gt.amneut(ni))
     .then
            xneutzstau=xneutzstau
     .      +g2s**2/dsto(i)/dz*oppl(ni,nj)*(
     .      ((atau(i,ni)*atau(i,nj)-btau(i,ni)*btau(i,nj))*vzz-
     .       (atau(i,ni)*atau(i,nj)+btau(i,ni)*btau(i,nj))*azz)*
     .      xmneut(nj)/xmneut(ni)*xmutau*(
     .      1/xmuz*(-4d0)*(1d0+xmuneut1+2d0*xmutau-uh-th)+16d0) +
     .      ((atau(i,ni)*atau(i,nj)-btau(i,ni)*btau(i,nj))*vzz+
     .       (atau(i,ni)*atau(i,nj)+btau(i,ni)*btau(i,nj))*azz)*
     .      xmneut(nj)/xmneut(ni)*(2d0/xmuz*((xmuneut1+1d0-uh-th)*
     .      2d0*xmutau+4d0*xmutau**2)+4d0*(1d0+xmuneut1-uh-th)) +
     .      ((btau(i,ni)*btau(i,nj)-atau(i,ni)*atau(i,nj))*vzz+
     .       (btau(i,ni)*btau(i,nj)+atau(i,ni)*atau(i,nj))*azz)*
     .      xmutau*(2d0/xmuz*(xmuneut1*(2d0*xmutau-th-uh+4d0)+
     .      2d0*xmutau-th-uh)+4d0*(2d0*xmutau-th-uh)) +
     .      ((btau(i,ni)*btau(i,nj)-atau(i,ni)*atau(i,nj))*vzz-
     .       (btau(i,ni)*btau(i,nj)+atau(i,ni)*atau(i,nj))*azz)*
     .      (2d0/xmuz*(xmuneut1*(-2d0*xmutau**2+xmutau*th-
     .      2d0*xmutau+xmutau*uh-2d0*xmutau)+xmutau*(-2d0*xmutau+uh)
     .      +xmutau*th)+4d0*(xmuneut1*(xmutau-uh+1d0)+xmutau*
     .      (xmutau-uh)+xmutau*(1d0-uh)+uh**2-uh)) +
     .      (atau(i,nj)*btau(i,ni)*(vzz+azz)-
     .       atau(i,ni)*btau(i,nj)*(vzz-azz))*
     .      xmneut(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*(2d0/xmuz*
     .      (xmuneut1*(1d0+xmutau-uh)+1d0+xmutau*(xmutau-uh)+xmutau*
     .      (xmutau-th-2d0*uh+3d0)+th*(uh-1d0)+uh*(uh-2d0))+
     .      4d0*(1d0+xmutau-th)) +
     .      (atau(i,nj)*btau(i,ni)*(vzz-azz)-
     .       atau(i,ni)*btau(i,nj)*(vzz+azz))*
     .      xmneut(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*
     .      (2d0/xmuz*(xmuneut1*(-1d0-xmutau+uh)-1d0+xmutau*(-1d0
     .      +uh)+uh*(2d0-th-uh)+th+xmutau*(-2d0*xmutau+th+2d0*uh
     .      -2d0))+8d0*(1d0+xmutau-uh)) +
     .      (atau(i,ni)*btau(i,nj)*(vzz-azz)-
     .       atau(i,nj)*btau(i,ni)*(vzz+azz))*dsqrt(xmutau)*sgn(ni)*(
     .      (-2d0)/xmuz*(xmuneut1-uh+xmutau)*(xmuneut1+2d0*xmutau-th
     .      -uh+1d0)+8d0*(xmuneut1+xmutau-uh)) +
     .      (atau(i,ni)*btau(i,nj)*(vzz+azz)-
     .       atau(i,nj)*btau(i,ni)*(vzz-azz))*dsqrt(xmutau)*sgn(ni)*(
     .      2d0/xmuz*(xmuneut1*(xmuneut1+3d0*xmutau-th-2d0*uh+1d0)+
     .      xmutau*(2d0*xmutau-th-2d0*uh)+xmutau*(1d0-uh)+uh*th+
     .      uh**2-uh)+4d0*(xmuneut1+xmutau-th)) )
     .      -g2s**2/dstob(i)/dz*oppl(ni,nj)*(
     .      ((-atau(i,ni)*atau(i,nj)+btau(i,ni)*btau(i,nj))*vzz+
     .       (atau(i,ni)*atau(i,nj)+btau(i,ni)*btau(i,nj))*azz)*
     .      xmneut(nj)/xmneut(ni)*xmutau*(
     .      1/xmuz*(-4d0)*(1d0+xmuneut1+2d0*xmutau-uh-th)+16d0) +
     .      ((-atau(i,ni)*atau(i,nj)+btau(i,ni)*btau(i,nj))*vzz-
     .       (atau(i,ni)*atau(i,nj)+btau(i,ni)*btau(i,nj))*azz)*
     .      xmneut(nj)/xmneut(ni)*(2d0/xmuz*(2d0*xmutau*(xmuneut1+1d0
     .      -th-uh)+4d0*xmutau**2)+4d0*(1d0+xmuneut1-uh-th)) +
     .      ((atau(i,ni)*atau(i,nj)-btau(i,ni)*btau(i,nj))*vzz-
     .       (atau(i,ni)*atau(i,nj)+btau(i,ni)*btau(i,nj))*azz)*
     .      xmutau*(2d0/xmuz*(xmuneut1*(2d0*xmutau-th-uh+4d0)+
     .      2d0*xmutau-th-uh)+4d0*(2d0*xmutau-th-uh)) +
     .      ((atau(i,ni)*atau(i,nj)-btau(i,ni)*btau(i,nj))*vzz+
     .       (atau(i,ni)*atau(i,nj)+btau(i,ni)*btau(i,nj))*azz)*
     .      (2d0/xmuz*(xmuneut1*(-2d0*xmutau**2+xmutau*th
     .      -2d0*xmutau+xmutau*uh-2d0*xmutau)+xmutau*(-2d0*xmutau
     .      +uh)+xmutau*th)+4d0*(xmuneut1*(xmutau-th+1d0)+xmutau*
     .      (xmutau-th)+xmutau*(1d0-th)+th**2-th)) +
     .      (atau(i,nj)*btau(i,ni)*(-vzz-azz)+
     .       atau(i,ni)*btau(i,nj)*(vzz-azz))*
     .      xmneut(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*
     .      (2d0/xmuz*(xmuneut1*(1d0+xmutau-th)+1d0+xmutau*
     .      (2d0*xmutau-2d0*th-uh+3d0)+th*(th-xmutau+uh-2d0)-uh)+
     .      4d0*(1d0+xmutau-uh)) +
     .      (atau(i,nj)*btau(i,ni)*(-vzz+azz)+
     .       atau(i,ni)*btau(i,nj)*(vzz+azz))*
     .      xmneut(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*
     .      (2d0/xmuz*(xmuneut1*(-1d0-xmutau+th)-1d0+xmutau*
     .      (-xmutau+th-1d0)+xmutau*(2d0*th+uh-2d0-xmutau)
     .      -th*(th+uh)+2d0*th+uh)+
     .      8d0*(1d0+xmutau-th)) +
     .      (atau(i,ni)*btau(i,nj)*(-vzz+azz)+
     .       atau(i,nj)*btau(i,ni)*(vzz+azz))*dsqrt(xmutau)*sgn(ni)*(
     .      2d0/xmuz*(xmuneut1*(-xmuneut1-3d0*xmutau+2d0*th+uh-1d0)
     .      +xmutau*(-2d0*xmutau+2d0*th+uh-1d0)+th*(xmutau-th-uh+1))
     .      +8d0*(xmuneut1+xmutau-th)) +
     .      (atau(i,ni)*btau(i,nj)*(-vzz-azz)+
     .       atau(i,nj)*btau(i,ni)*(vzz-azz))*dsqrt(xmutau)*sgn(ni)*(
     .      2d0/xmuz*(xmuneut1*(xmuneut1+3d0*xmutau-uh-2d0*th+1d0)+
     .      xmutau*(xmutau-th+1d0)+xmutau*(xmutau-2d0*th-uh)+
     .      uh*th+th**2-th)+4d0*(xmuneut1+xmutau-uh)) )
      endif
         enddo
      else
         xneutzstau=0d0
      endif

c -------------------------------------------------------------------- c
c                           H(j)-stau interference
c -------------------------------------------------------------------- c
      do j=1,3
      dhl(j)   = y3-smass(j)**2/amneut(ni)**2
      xneuthlstau(j)=0d0

      if ((amneut(nj)+2d0*amtau).le.amneut(ni)) then
         do i=1,2
       if ((smass(j)+amneut(nj)).gt.amneut(ni)
     ..and.xmustau(i).gt.1d0)then
            xneuthlstau(j)=xneuthlstau(j)
     .      +2d0*g2s**2/dhl(j)/dsto(i)*(scaltau*(S(j,2))/dsqrt(2d0))*
     .       (-2d0)*hchichi(j,ni,nj)*(
     .       (atau(i,nj)*atau(i,ni)+btau(i,ni)*btau(i,nj))*(
     .       dsqrt(xmutau)*sgn(ni)*(xmuneut1+xmutau-uh) +
     .       dsqrt(xmutau)*sgn(ni)*(-xmuneut1-xmutau+th) +
     .       xmneut(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*
     .       (-1d0-xmutau+th) +
     .       xmneut(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*
     .       (1d0+xmutau-uh) )
     .       +(atau(i,nj)*btau(i,ni)+atau(i,ni)*btau(i,nj))*(
     .       xmneut(nj)/xmneut(ni)*(-1d0-xmuneut1+uh+th) +
     .       2d0*xmneut(nj)/xmneut(ni)*xmutau +
     .       xmutau*(-2d0*xmutau+th+uh) +
     .       (uh**2+uh*th-uh*(1d0+2d0*xmutau+xmuneut1)+xmutau+
     .       xmuneut1*xmutau)) )
     .     +2d0*g2s**2/dhl(j)/dstob(i)*(scaltau*(S(j,2))/dsqrt(2d0))*
     .       (-2d0)*hchichi(j,ni,nj)*(
     .       (atau(i,nj)*atau(i,ni)+btau(i,ni)*btau(i,nj))*(
     .       dsqrt(xmutau)*sgn(ni)*(uh-xmuneut1-xmutau) +
     .       dsqrt(xmutau)*sgn(ni)*(-th+xmutau+xmuneut1) +
     .       dsqrt(xmutau)*sgn(ni)*xmneut(nj)/xmneut(ni)*(-th+1d0+
     .       xmutau) +
     .       dsqrt(xmutau)*sgn(ni)*xmneut(nj)/xmneut(ni)*(uh-1d0-
     .       xmutau)) +
     .       (atau(i,nj)*btau(i,ni)+atau(i,ni)*btau(i,nj))*(
     .       2d0*xmutau*xmneut(nj)/xmneut(ni) +
     .       (uh*th+th**2-th*(1d0+2d0*xmutau+xmuneut1)+xmutau+
     .       xmutau*xmuneut1) +
     .       xmneut(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0) +
     .       xmutau*(uh+th-2d0*xmutau)) )
            endif
         enddo
      else
         xneuthlstau(j)=0d0
      endif

      enddo
c -------------------------------------------------------------------- c
c 			   A(j)-stau interference
c -------------------------------------------------------------------- c
      do j=1,2
      da(j)   = y3-pmass(j)**2/amneut(ni)**2
      xneutastau(j)=0d0

      if ((amneut(nj)+2d0*amtau).le.amneut(ni)) then
         do i=1,2
      if ((pmass(j)+amneut(nj)).gt.amneut(ni).and.
     .xmustau(i).gt.1d0)then
            xneutastau(j)=xneutastau(j)
     .      +2d0*g2s**2/da(j)/dsto(i)*((-scaltau*p(j,2))/dsqrt(2d0))*
     .       2d0*achichi(j,ni,nj)*(
     .       (atau(i,nj)*atau(i,ni)+btau(i,ni)*btau(i,nj))*(
     .       dsqrt(xmutau)*sgn(ni)*(xmuneut1+xmutau-uh) +
     .       dsqrt(xmutau)*sgn(ni)*(-xmuneut1-xmutau+th)*(-1d0) +
     .       xmneut(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*(-1d0-xmutau
     .       +th)*(-1d0) +
     .       xmneut(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*(1d0+xmutau
     .       -uh) )
     .       +(atau(i,nj)*btau(i,ni)+atau(i,ni)*btau(i,nj))*(
     .       xmneut(nj)/xmneut(ni)*(-1d0-xmuneut1+uh+th)*(-1d0) +
     .       2d0*xmneut(nj)/xmneut(ni)*xmutau +
     .       xmutau*(-2d0*xmutau+th+uh)*(-1d0) +
     .       (uh**2+uh*th-uh*(1d0+2d0*xmutau+xmuneut1)+xmutau+
     .       xmuneut1*xmutau)) )
     .     +2d0*g2s**2/da(j)/dstob(i)*((-scaltau*p(j,2))/dsqrt(2d0))*
     .       2d0*achichi(j,ni,nj)*(
     .       (atau(i,nj)*atau(i,ni)+btau(i,ni)*btau(i,nj))*(
     .       dsqrt(xmutau)*sgn(ni)*(uh-xmuneut1-xmutau)*(-1d0) +
     .       dsqrt(xmutau)*sgn(ni)*(-th+xmutau+xmuneut1) +
     .       dsqrt(xmutau)*sgn(ni)*xmneut(nj)/xmneut(ni)*(-th+1d0+
     .       xmutau) +
     .       dsqrt(xmutau)*sgn(ni)*xmneut(nj)/xmneut(ni)*(uh-1d0-
     .       xmutau)*(-1d0)) +
     .       (atau(i,nj)*btau(i,ni)+atau(i,ni)*btau(i,nj))*(
     .       2d0*xmutau*xmneut(nj)/xmneut(ni) +
     .       (uh*th+th**2-th*(1d0+2d0*xmutau+xmuneut1)+xmutau+
     .       xmutau*xmuneut1) +
     .       xmneut(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0)*(-1d0) +
     .       xmutau*(uh+th-2d0*xmutau)*(-1d0)) )
      endif
         enddo
      else
         xneutastau(j)=0d0
      endif

      enddo

c -------------------------------------------------------------------- c
c 	                interference Z and A(i)
C  -------------------Z/h/H INTERFERENCE IS MISSING
c -------------------------------------------------------------------- c
        do i=1,2
      da(i)   = y3-pmass(i)**2/amneut(ni)**2

      xneutza(i)=0d0

      if ((amneut(nj)+2d0*amtau).le.amneut(ni)) then
         if ((pmass(i)+amneut(nj)).gt.amneut(ni).and.
     .(mz+amneut(nj)).gt.amneut(ni))then
         xneutza(i)=xneutza(i)-4d0*g2s**2/da(i)/dz*azz*
     .    (-scaltau*P(i,2))/dsqrt(2d0)*2d0*achichi(i,ni,nj)
     .    *oppl(ni,nj)*(
     .    xmneut(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*(
     .    4d0/xmuz*(2d0+xmuneut1*(2d0*xmutau-th-uh+2d0)+xmutau*
     .    (-uh-th+1d0+xmutau)+xmutau*(3d0*xmutau-3d0*(th+uh)+5d0)
     .    +(th+uh)**2-3d0*(th+uh))+
     .    4d0*(th-xmutau-1d0)+4d0*(uh-xmutau-1d0))
     .    -dsqrt(xmutau)*sgn(ni)*(
     .    4d0/xmuz*(-2d0*xmuneut1**2+xmuneut1*(-6d0*xmutau+
     .    3d0*(th+uh)-2d0)+xmutau*(-xmutau+th+uh-1d0)+xmutau*
     .    (-3d0*xmutau+3d0*(uh+th)-1d0)-(th+uh)**2+(th+uh))+
     .    4d0*(xmuneut1+xmutau-uh)+4d0*(xmuneut1+xmutau-th)) )
      endif
      else
         xneutza(i)=0d0
      endif

      enddo

      NS_neuttau = xneutstau+xneutztau+xneutzstau

      do i=1,3
         NS_neuttau=NS_neuttau+xneuthlstau(i)
         do j=1,3
            NS_neuttau=NS_neuttau+xneuthl(i,j)
         enddo
      enddo

      do i=1,2
         NS_neuttau=NS_neuttau+xneutastau(i)+xneutza(i)
         do J=1,2
            NS_neuttau=NS_neuttau+xneuta(i,j)
         enddo
      enddo

      end
c ==================================================================== c
c =====================  neutralino nu_e nu_ebar ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_neutnue(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION opl(2,2),opr(2,2),oppl(5,5),oppr(5,5)
      DOUBLE PRECISION dsnl(2),dsnlb(2)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION xmuneut1,x1,x2,x3,y3,y2,y1,xmusnl(2),
     .xneutsnl,dz,xneutznl,xneutzsnl,xmuz
*
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup4/opl,opr,oppl,oppr
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
*
      xmuneut1 = amneut(nj)**2/amneut(ni)**2
      xmuz=mz**2/amneut(ni)**2
      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                           snu_e exchange
c -------------------------------------------------------------------- c
      xmusnl(1) = asne1**2/amneut(ni)**2
      xmusnl(2) = asne2**2/amneut(ni)**2

      dsnl(1)  = 1d0-x1-xmusnl(1)
      dsnl(2)  = 1d0-x1-xmusnl(2)
      dsnlb(1) = 1d0-x2-xmusnl(1)
      dsnlb(2) = 1d0-x2-xmusnl(2)

      xneutsnl=0d0

      if (amneut(nj).le.amneut(ni)) then
         do i=1,2
            do k=1,2
       if (xmusnl(k).gt.1d0.and.xmusnl(i).gt.1d0)then
               xneutsnl=xneutsnl
     .              +g2s**2/dsnl(k)/dsnl(i)*x1*y1*
     .               (ane(i,ni)*ane(k,ni)+bne(i,ni)*bne(k,ni))*
     .               (ane(i,nj)*ane(k,nj)+bne(i,nj)*bne(k,nj))
     .              +g2s**2/dsnlb(k)/dsnlb(i)*x2*y2*
     .               (ane(i,ni)*ane(k,ni)+bne(i,ni)*bne(k,ni))*
     .               (ane(i,nj)*ane(k,nj)+bne(i,nj)*bne(k,nj))
     .              +g2s**2/dsnlb(k)/dsnl(i)*
     .               ( ( ane(i,ni)*bne(k,ni)*ane(k,nj)*bne(i,nj)
     .                  +ane(i,nj)*bne(k,nj)*ane(k,ni)*bne(i,ni))
     .                 *(-x1*y1-x2*y2+x3*y3)
     .                +2d0*xmneut(nj)/xmneut(ni)*y3*
     .                ( ane(i,ni)*ane(k,ni)*ane(k,nj)*ane(i,nj)
     .                 +bne(i,nj)*bne(k,nj)*bne(k,ni)*bne(i,ni)))
       endif
            enddo
         enddo
      else
         xneutsnl=0d0
      endif

c -------------------------------------------------------------------- c
c                             Z exchange
c -------------------------------------------------------------------- c
      dz   = y3-mz**2/amneut(ni)**2

      xneutznl=0d0

      if (amneut(nj).le.amneut(ni)) then
      if ((mz+amneut(nj)).gt.amneut(ni))then
         xneutznl=g2s**2*4d0/dz**2*
     .        (((azzneutneut+vzzneutneut)**2*oppl(ni,nj)**2
     .        + (azzneutneut-vzzneutneut)**2*oppr(ni,nj)**2)*x2*y2
     .        +((azzneutneut+vzzneutneut)**2*oppr(ni,nj)**2
     .        + (azzneutneut-vzzneutneut)**2*oppl(ni,nj)**2)*x1*y1
     .        -4d0*xmneut(nj)/xmneut(ni)*oppl(ni,nj)*oppr(ni,nj)
     .        *(azzneutneut**2+vzzneutneut**2)*y3 )
      endif
      else
         xneutznl=0d0
      endif
c -------------------------------------------------------------------- c
c                          Z-snue interference
c -------------------------------------------------------------------- c
      xneutzsnl=0d0

      if (amneut(nj).le.amneut(ni)) then
         do i=1,2
      if (xmusnl(i).gt.1d0.and.(mz+amneut(nj)).gt.amneut(ni))then
            xneutzsnl=xneutzsnl-g2s**2*4d0/dsnl(i)/dz*
     .           ((ane(i,ni)*ane(i,nj)*oppr(ni,nj)*
     .           (azzneutneut+vzzneutneut)+
     .           bne(i,ni)*bne(i,nj)*oppl(ni,nj)*
     .           (-azzneutneut+vzzneutneut))*x1*y1
     .           -(ane(i,nj)*ane(i,ni)*oppl(ni,nj)*
     .           (azzneutneut+vzzneutneut)
     .           +bne(i,nj)*bne(i,ni)*oppr(ni,nj)*
     .           (-azzneutneut+vzzneutneut))*xmneut(nj)/xmneut(ni)*y3)
     .           +g2s**2*4d0/dsnlb(i)/dz
     .           *((ane(i,ni)*ane(i,nj)*oppl(ni,nj)*
     .           (azzneutneut+vzzneutneut)+bne(i,ni)*bne(i,nj)*
     .           oppr(ni,nj)*(-azzneutneut+vzzneutneut))*x2*y2
     .           -(ane(i,nj)*ane(i,ni)*oppr(ni,nj)*
     .           (azzneutneut+vzzneutneut)
     .           +bne(i,nj)*bne(i,ni)*oppl(ni,nj)*
     .           (-azzneutneut+vzzneutneut))*xmneut(nj)/xmneut(ni)*y3)
      endif
         enddo
      else
         xneutzsnl=0d0
      endif

      NS_neutnue = xneutsnl+xneutznl+xneutzsnl

      end

c ==================================================================== c
c =====================  neutralino nu_mu nu_mubar =================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_neutnumu(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION opl(2,2),opr(2,2),oppl(5,5),oppr(5,5)
      DOUBLE PRECISION dsnl(2),dsnlb(2)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION xmuneut1,x1,x2,x3,y3,y2,y1,xmusnl(2),
     .xneutsnl,dz,xneutznl,xneutzsnl,xmuz
*
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup4/opl,opr,oppl,oppr
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
*
      xmuneut1 = amneut(nj)**2/amneut(ni)**2
      xmuz=mz**2/amneut(ni)**2
      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                           snu_mu exchange
c -------------------------------------------------------------------- c
      xmusnl(1) = asnmu1**2/amneut(ni)**2
      xmusnl(2) = asnmu2**2/amneut(ni)**2

      dsnl(1)  = 1d0-x1-xmusnl(1)
      dsnl(2)  = 1d0-x1-xmusnl(2)
      dsnlb(1) = 1d0-x2-xmusnl(1)
      dsnlb(2) = 1d0-x2-xmusnl(2)

      xneutsnl=0d0

      if (amneut(nj).le.amneut(ni)) then
         do i=1,2
            do k=1,2
       if (xmusnl(k).gt.1d0.and.xmusnl(i).gt.1d0)then
               xneutsnl=xneutsnl
     .              +g2s**2/dsnl(k)/dsnl(i)*x1*y1*
     .               (anmu(i,ni)*anmu(k,ni)+bnmu(i,ni)*bnmu(k,ni))*
     .               (anmu(i,nj)*anmu(k,nj)+bnmu(i,nj)*bnmu(k,nj))
     .              +g2s**2/dsnlb(k)/dsnlb(i)*x2*y2*
     .               (anmu(i,ni)*anmu(k,ni)+bnmu(i,ni)*bnmu(k,ni))*
     .               (anmu(i,nj)*anmu(k,nj)+bnmu(i,nj)*bnmu(k,nj))
     .              +g2s**2/dsnlb(k)/dsnl(i)*
     .               ( ( anmu(i,ni)*bnmu(k,ni)*anmu(k,nj)*bnmu(i,nj)
     .                  +anmu(i,nj)*bnmu(k,nj)*anmu(k,ni)*bnmu(i,ni))
     .                 *(-x1*y1-x2*y2+x3*y3)
     .                +2d0*xmneut(nj)/xmneut(ni)*y3*
     .                ( anmu(i,ni)*anmu(k,ni)*anmu(k,nj)*anmu(i,nj)
     .                 +bnmu(i,nj)*bnmu(k,nj)*bnmu(k,ni)*bnmu(i,ni)))
       endif
            enddo
         enddo
      else
         xneutsnl=0d0
      endif

c -------------------------------------------------------------------- c
c                             Z exchange
c -------------------------------------------------------------------- c
      dz   = y3-mz**2/amneut(ni)**2

      xneutznl=0d0

      if (amneut(nj).le.amneut(ni)) then
      if ((mz+amneut(nj)).gt.amneut(ni))then
         xneutznl=g2s**2*4d0/dz**2*
     .        (((azzneutneut+vzzneutneut)**2*oppl(ni,nj)**2
     .        + (azzneutneut-vzzneutneut)**2*oppr(ni,nj)**2)*x2*y2
     .        +((azzneutneut+vzzneutneut)**2*oppr(ni,nj)**2
     .        + (azzneutneut-vzzneutneut)**2*oppl(ni,nj)**2)*x1*y1
     .        -4d0*xmneut(nj)/xmneut(ni)*oppl(ni,nj)*oppr(ni,nj)
     .        *(azzneutneut**2+vzzneutneut**2)*y3 )
      endif
      else
         xneutznl=0d0
      endif
c -------------------------------------------------------------------- c
c                          Z-snumu interference
c -------------------------------------------------------------------- c
      xneutzsnl=0d0

      if (amneut(nj).le.amneut(ni)) then
         do i=1,2
      if (xmusnl(i).gt.1d0.and.(mz+amneut(nj)).gt.amneut(ni))then
            xneutzsnl=xneutzsnl-g2s**2*4d0/dsnl(i)/dz*
     .           ((anmu(i,ni)*anmu(i,nj)*oppr(ni,nj)*
     .           (azzneutneut+vzzneutneut)+
     .           bnmu(i,ni)*bnmu(i,nj)*oppl(ni,nj)*
     .           (-azzneutneut+vzzneutneut))*x1*y1
     .           -(anmu(i,nj)*anmu(i,ni)*oppl(ni,nj)*
     .           (azzneutneut+vzzneutneut)
     .           +bnmu(i,nj)*bnmu(i,ni)*oppr(ni,nj)*
     .           (-azzneutneut+vzzneutneut))*xmneut(nj)/xmneut(ni)*y3)
     .           +g2s**2*4d0/dsnlb(i)/dz
     .           *((anmu(i,ni)*anmu(i,nj)*oppl(ni,nj)*
     .           (azzneutneut+vzzneutneut)+bnmu(i,ni)*bnmu(i,nj)*
     .           oppr(ni,nj)*(-azzneutneut+vzzneutneut))*x2*y2
     .           -(anmu(i,nj)*anmu(i,ni)*oppr(ni,nj)*
     .           (azzneutneut+vzzneutneut)
     .           +bnmu(i,nj)*bnmu(i,ni)*oppl(ni,nj)*
     .           (-azzneutneut+vzzneutneut))*xmneut(nj)/xmneut(ni)*y3)
      endif
         enddo
      else
         xneutzsnl=0d0
      endif

      NS_neutnumu = xneutsnl+xneutznl+xneutzsnl

      end

c ==================================================================== c
c ===================  neutralino nu_tau nu_taubar =================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_neutnutau(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION opl(2,2),opr(2,2),oppl(5,5),oppr(5,5)
      DOUBLE PRECISION dsnt(2),dsntb(2)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION xmuneut1,x1,x2,x3,y3,y2,y1,xmusnt(2),
     .xneutsnt,dz,xneutznt,xneutzsnt
*
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup4/opl,opr,oppl,oppr
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2

      xmuneut1 = amneut(nj)**2/amneut(ni)**2

      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                           snu_tau exchange
c -------------------------------------------------------------------- c
      xmusnt(1) = asntau1**2/amneut(ni)**2
      xmusnt(2) = asntau2**2/amneut(ni)**2

      dsnt(1)  = 1d0-x1-xmusnt(1)
      dsnt(2)  = 1d0-x1-xmusnt(2)
      dsntb(1) = 1d0-x2-xmusnt(1)
      dsntb(2) = 1d0-x2-xmusnt(2)

      xneutsnt=0d0

      if (amneut(nj).le.amneut(ni)) then

         do i=1,2
            do k=1,2
      if (xmusnt(i).gt.1d0.and.xmusnt(k).gt.1d0)then
               xneutsnt=xneutsnt
     .           +g2s**2/dsnt(k)/dsnt(i)*x1*y1*
     .            (antau(i,ni)*antau(k,ni)+bntau(i,ni)*bntau(k,ni))*
     .            (antau(i,nj)*antau(k,nj)+bntau(i,nj)*bntau(k,nj))
     .           +g2s**2/dsntb(k)/dsntb(i)*x2*y2*
     .            (antau(i,ni)*antau(k,ni)+bntau(i,ni)*bntau(k,ni))*
     .            (antau(i,nj)*antau(k,nj)+bntau(i,nj)*bntau(k,nj))
     .           +g2s**2/dsntb(k)/dsnt(i)*
     .            ( ( antau(i,ni)*bntau(k,ni)*antau(k,nj)*bntau(i,nj)
     .              + antau(i,nj)*bntau(k,nj)*antau(k,ni)*bntau(i,ni))
     .              *(-x1*y1-x2*y2+x3*y3)
     .            +2d0*xmneut(nj)/xmneut(ni)*y3*
     .            ( antau(i,ni)*antau(k,ni)*antau(k,nj)*antau(i,nj)
     .             +bntau(i,nj)*bntau(k,nj)*bntau(k,ni)*bntau(i,ni)))
      endif
            enddo
         enddo
      else
         xneutsnt=0d0
      endif

c -------------------------------------------------------------------- c
c                             Z exchange
c -------------------------------------------------------------------- c

      dz   = y3-mz**2/amneut(ni)**2

      xneutznt=0d0

      if (amneut(nj).le.amneut(ni)) then
      if ((mz+amneut(nj)).gt.amneut(ni))then
         xneutznt=g2s**2*4d0/dz**2*
     .        (((azzneutneut+vzzneutneut)**2*oppl(ni,nj)**2
     .        + (azzneutneut-vzzneutneut)**2*oppr(ni,nj)**2)*x2*y2
     .        +((azzneutneut+vzzneutneut)**2*oppr(ni,nj)**2
     .        + (azzneutneut-vzzneutneut)**2*oppl(ni,nj)**2)*x1*y1
     .        -4d0*xmneut(nj)/xmneut(ni)*oppl(ni,nj)*oppr(ni,nj)
     .        *(azzneutneut**2+vzzneutneut**2)*y3 )
      endif
      else
         xneutznt=0d0
      endif
c -------------------------------------------------------------------- c
c                      Z-snu tau interference
c -------------------------------------------------------------------- c
      xneutzsnt=0d0

      if (amneut(nj).le.amneut(ni)) then
         do i=1,2
      if ((mz+amneut(nj)).gt.amneut(ni)
     ..and.xmusnt(i).gt.1d0)then
            xneutzsnt=xneutzsnt-g2s**2*4d0/dsnt(i)/dz
     .           *((antau(i,ni)*antau(i,nj)*oppr(ni,nj)
     .           *(azzneutneut+vzzneutneut)+
     .           bntau(i,ni)*bntau(i,nj)*oppl(ni,nj)*
     .           (-azzneutneut+vzzneutneut))*x1*y1
     .           -(antau(i,nj)*antau(i,ni)*oppl(ni,nj)*
     .           (azzneutneut+vzzneutneut)
     .           +bntau(i,nj)*bntau(i,ni)*oppr(ni,nj)*
     .           (-azzneutneut+vzzneutneut)
     .           )*xmneut(nj)/xmneut(ni)*y3)
     .           +g2s**2*4d0/dsntb(i)/dz
     .           *((antau(i,ni)*antau(i,nj)*oppl(ni,nj)*
     .           (azzneutneut+vzzneutneut)+
     .           bntau(i,ni)*bntau(i,nj)*oppr(ni,nj)*
     .           (-azzneutneut+vzzneutneut))*x2*y2
     .           -(antau(i,nj)*antau(i,ni)*oppr(ni,nj)*
     .           (azzneutneut+vzzneutneut)
     .           +bntau(i,nj)*bntau(i,ni)*oppl(ni,nj)*
     .           (-azzneutneut+vzzneutneut)
     .           )*xmneut(nj)/xmneut(ni)*y3)
      endif
         enddo
      else
         xneutzsnt=0d0
      endif

      NS_neutnutau = xneutsnt+xneutznt+xneutzsnt

      end
c ==================================================================== c
c =========================  chargino- e+ nu_e ======================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_chelne(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsel(2),dsne(2)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION ql(5,2),qr(5,2),or(5,2),ol(5,2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION awff,vwff
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION xmuneut1,x1,x2,x3,y3,y2,y1,xmusel(2),
     .xmusn(2),xneutsf,dw,xneutwsf,xneutwel,xmuw
*
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup18/awff,vwff
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS

      xmuneut1 = amchar(nj)**2/amneut(ni)**2

      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                        sfermion exchange
c -------------------------------------------------------------------- c
      xmusel(1) = ase1**2/amneut(ni)**2
      xmusel(2) = ase2**2/amneut(ni)**2
      xmusn(1)  = asne1**2/amneut(ni)**2
      xmusn(2)  = asne2**2/amneut(ni)**2

      dsel(1) = 1d0-x2-xmusel(1)
      dsel(2) = 1d0-x2-xmusel(2)
      dsne(1) = 1d0-x1-xmusn(1)
      dsne(2) = 1d0-x1-xmusn(2)

      xneutsf=0d0

      if(amchar(nj).le.amneut(ni)) then
         do i=1,2
            do k=1,2
       if (xmusel(i).gt.1d0.and.xmusel(k).gt.1d0)then
               xneutsf=xneutsf
     .              +g2s**2/dsel(k)/dsel(i)*x2*y2*
     .              (ae(i,ni)*ae(k,ni)+be(i,ni)*be(k,ni))*
     .              (ale(i,nj)*ale(k,nj))
       endif
       if (xmusn(i).gt.1d0.and.xmusn(k).gt.1d0)then
               xneutsf=xneutsf
     .              +g2s**2/dsne(k)/dsne(i)*x1*y1*
     .              (ane(i,ni)*ane(k,ni)+bne(i,ni)*bne(k,ni))*
     .              (alsne(i,nj)*alsne(k,nj)+blsne(i,nj)*blsne(k,nj))
       endif
       if (xmusn(i).gt.1d0.and.xmusel(k).gt.1d0)then
               xneutsf=xneutsf
     .              +g2s**2/dsne(i)/dsel(k)*
     .              (( ane(i,ni)*be(k,ni)*ale(k,nj)*blsne(i,nj))
     .              *(-x1*y1-x2*y2+x3*y3)
     .              +2d0*xmchar(nj)/xmneut(ni)*y3*
     .              ( ae(k,ni)*ane(i,ni)*alsne(i,nj)*ale(k,nj)))
       endif
            enddo
         enddo
      else
         xneutsf=0d0
      endif

c -------------------------------------------------------------------- c
c                             W+ exchange
c -------------------------------------------------------------------- c
      xmuw=mw**2/amneut(ni)**2
      dw   = y3-mw**2/amneut(ni)**2
      xneutwel=0d0
      if(amchar(nj).le.amneut(ni)) then
         if ((mw+amchar(nj)).gt.amneut(ni))then
         xneutwel=g2s**2*4d0/dw**2*
     .        (4d0*vwff**2*or(ni,nj)**2*x2*y2
     .        +4d0*vwff**2*ol(ni,nj)**2*x1*y1
     .        -8d0*xmchar(nj)/xmneut(ni)*ol(ni,nj)*or(ni,nj)*
     .        vwff**2*y3 )
         endif
      else
         xneutwel=0d0
      endif

c -------------------------------------------------------------------- c
c                        W+ - sfermion interference
c -------------------------------------------------------------------- c
      xneutwsf=0d0

      if(amchar(nj).le.amneut(ni)) then
         do i=1,2
      if ((mw+amchar(nj)).gt.amneut(ni).and.
     .xmusel(i).gt.1d0)then
            xneutwsf=xneutwsf-g2s**2*4d0/dsel(i)/dw
     .      *(ae(i,ni)*ale(i,nj)*or(ni,nj)*2d0*vwff*x2*y2
     .       -ale(i,nj)*ae(i,ni)*ol(ni,nj)*2d0*vwff*
     .        xmchar(nj)/xmneut(ni)*y3)
      endif
      if ((mw+amchar(nj)).gt.amneut(ni).and.
     .xmusn(i).gt.1d0)then
           xneutwsf=xneutwsf
     .       +g2s**2*4d0/dsne(i)/dw
     .      *(ane(i,ni)*alsne(i,nj)*ol(ni,nj)*2d0*vwff*x1*y1
     .       -alsne(i,nj)*ane(i,ni)*or(ni,nj)*2d0*vwff*
     .        xmchar(nj)/xmneut(ni)*y3)
      endif
         enddo
      else
         xneutwsf=0d0
      endif

      NS_chelne = xneutsf+xneutwel+xneutwsf


      end
c ==================================================================== c
c =========================  chargino- mu+ nu_mu ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_chmunmu(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsel(2),dsne(2)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION ql(5,2),qr(5,2),or(5,2),ol(5,2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION awff,vwff
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION xmuneut1,x1,x2,x3,y3,y2,y1,xmusel(2),
     .xmusn(2),xneutsf,dw,xneutwsf,xneutwel,xmuw
*
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup18/awff,vwff
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS

      xmuneut1 = amchar(nj)**2/amneut(ni)**2

      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                        sfermion exchange
c -------------------------------------------------------------------- c
      xmusel(1) = asmu1**2/amneut(ni)**2
      xmusel(2) = asmu2**2/amneut(ni)**2
      xmusn(1)  = asnmu1**2/amneut(ni)**2
      xmusn(2)  = asnmu2**2/amneut(ni)**2

      dsel(1) = 1d0-x2-xmusel(1)
      dsel(2) = 1d0-x2-xmusel(2)
      dsne(1) = 1d0-x1-xmusn(1)
      dsne(2) = 1d0-x1-xmusn(2)

      xneutsf=0d0

      if(amchar(nj).le.amneut(ni)) then
         do i=1,2
            do k=1,2
       if (xmusel(i).gt.1d0.and.xmusel(k).gt.1d0)then
               xneutsf=xneutsf
     .              +g2s**2/dsel(k)/dsel(i)*x2*y2*
     .              (amu(i,ni)*amu(k,ni)+bmu(i,ni)*bmu(k,ni))*
     .              (almu(i,nj)*almu(k,nj))
       endif
       if (xmusn(i).gt.1d0.and.xmusn(k).gt.1d0)then
               xneutsf=xneutsf
     .              +g2s**2/dsne(k)/dsne(i)*x1*y1*
     .              (anmu(i,ni)*anmu(k,ni)+bnmu(i,ni)*bnmu(k,ni))*
     .              (alsnm(i,nj)*alsnm(k,nj)+blsnm(i,nj)*blsnm(k,nj))
       endif
       if (xmusn(i).gt.1d0.and.xmusel(k).gt.1d0)then
               xneutsf=xneutsf
     .              +g2s**2/dsne(i)/dsel(k)*
     .              (( anmu(i,ni)*bmu(k,ni)*almu(k,nj)*blsnm(i,nj))
     .              *(-x1*y1-x2*y2+x3*y3)
     .              +2d0*xmchar(nj)/xmneut(ni)*y3*
     .              ( amu(k,ni)*anmu(i,ni)*alsnm(i,nj)*almu(k,nj)))
       endif
            enddo
         enddo
      else
         xneutsf=0d0
      endif

c -------------------------------------------------------------------- c
c                             W+ exchange
c -------------------------------------------------------------------- c
      xmuw=mw**2/amneut(ni)**2
      dw   = y3-mw**2/amneut(ni)**2
      xneutwel=0d0
      if(amchar(nj).le.amneut(ni)) then
         if ((mw+amchar(nj)).gt.amneut(ni))then
         xneutwel=g2s**2*4d0/dw**2*
     .        (4d0*vwff**2*or(ni,nj)**2*x2*y2
     .        +4d0*vwff**2*ol(ni,nj)**2*x1*y1
     .        -8d0*xmchar(nj)/xmneut(ni)*ol(ni,nj)*or(ni,nj)*
     .        vwff**2*y3 )
         endif
      else
         xneutwel=0d0
      endif

c -------------------------------------------------------------------- c
c                        W+ - sfermion interference
c -------------------------------------------------------------------- c
      xneutwsf=0d0

      if(amchar(nj).le.amneut(ni)) then
         do i=1,2
      if ((mw+amchar(nj)).gt.amneut(ni).and.
     .xmusel(i).gt.1d0)then
            xneutwsf=xneutwsf-g2s**2*4d0/dsel(i)/dw
     .      *(amu(i,ni)*almu(i,nj)*or(ni,nj)*2d0*vwff*x2*y2
     .       -almu(i,nj)*amu(i,ni)*ol(ni,nj)*2d0*vwff*
     .        xmchar(nj)/xmneut(ni)*y3)
      endif
      if ((mw+amchar(nj)).gt.amneut(ni).and.
     .xmusn(i).gt.1d0)then
           xneutwsf=xneutwsf
     .       +g2s**2*4d0/dsne(i)/dw
     .      *(anmu(i,ni)*alsnm(i,nj)*ol(ni,nj)*2d0*vwff*x1*y1
     .       -alsnm(i,nj)*anmu(i,ni)*or(ni,nj)*2d0*vwff*
     .        xmchar(nj)/xmneut(ni)*y3)
      endif
         enddo
      else
         xneutwsf=0d0
      endif

      NS_chmunmu = xneutsf+xneutwel+xneutwsf


      end
c ==================================================================== c
c =======================  chargino- tau+ nu_tau ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_chtauntau(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dstau(2),dsntau(2),sgn(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION ql(5,2),qr(5,2),or(5,2),ol(5,2)
      DOUBLE PRECISION hchichi(3,5,5),achichi(2,5,5)
      DOUBLE PRECISION hchachaR(3,2,2),hchachaL(3,2,2),achachaR(2,2,2),
     .     achachaL(2,2,2)
      DOUBLE PRECISION achtop,vchtop,achtau,vchtau
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION awff,vwff
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION xmuneut1,xmutau,xmustau(2),
     .xmusnt(2),x1,x2,x3,y3,uh,th,xneutsf,xmuw,dw,xneutwtau,
     .dh,xneuthtau,xneutwsf,xneuthsf,xneuthw
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
*
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup15/achtop,vchtop,achtau,vchtau
      COMMON/NS_coup18/awff,vwff
      COMMON/NS_coupNMSSM/hchichi,achichi,hchachaR,hchachaL,achachaR,
     .     achachaL
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_indices/ni,nj
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW

      do i=1,5
         sgn(i) = 1d0
         if(xmneut(i).gt.0d0) then
            sgn(i) = 1d0
         elseif(xmneut(i).lt.0d0) then
            sgn(i) = -1d0
         endif
      enddo

      xmuneut1 = amchar(nj)**2/amneut(ni)**2
      xmutau   = amtau**2/amneut(ni)**2
      xmustau(1) = astau1**2/amneut(ni)**2
      xmustau(2) = astau2**2/amneut(ni)**2
      xmusnt(1)  = asntau1**2/amneut(ni)**2
      xmusnt(2)  = asntau2**2/amneut(ni)**2

      x3 = 2d0-x1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                        sfermion exchange
c -------------------------------------------------------------------- c
      dsntau(1) = 1d0-x1-xmusnt(1)
      dsntau(2) = 1d0-x1-xmusnt(2)
      dstau(1)  = 1d0-x2-xmustau(1)+xmutau
      dstau(2)  = 1d0-x2-xmustau(2)+xmutau

      uh = 1d0-x1
      th = 1d0-x2+xmutau

      xneutsf=0d0

      if((amchar(nj)+amtau).le.amneut(ni)) then
         do i=1,2
            do k=1,2
       if (xmustau(i).gt.1d0.and.xmustau(k).gt.1d0 )then
               xneutsf=xneutsf
     .          +g2s**2/dstau(k)/dstau(i)*(
     .          (atau(i,ni)*atau(k,ni)+btau(i,ni)*btau(k,ni))*
     .          (altau(i,nj)*altau(k,nj))*
     .          (xmuneut1*(th-xmutau-1d0)+xmutau*(
     .          +th)+th*(-th+1d0))+
     .          (atau(i,ni)*btau(k,ni)+btau(i,ni)*atau(k,ni))*
     .          (altau(i,nj)*altau(k,nj))*
     .          dsqrt(xmutau)*sgn(ni)*2d0*(th-xmuneut1) )
       endif
       if (xmusnt(i).gt.1d0.and.xmusnt(k).gt.1d0 )then
               xneutsf=xneutsf
     .          +g2s**2/dsntau(k)/dsntau(i)*(
     .          (antau(i,ni)*antau(k,ni)+bntau(i,ni)*bntau(k,ni))*
     .          (alsnt(i,nj)*alsnt(k,nj)+blsnt(i,nj)*blsnt(k,nj))*
     .          (xmuneut1*(+uh-1d0)+xmutau*
     .          (uh-1d0)+uh*(-uh+1d0))+
     .          (antau(i,ni)*antau(k,ni)+bntau(i,ni)*bntau(k,ni))*
     .          (alsnt(i,nj)*blsnt(k,nj)+blsnt(i,nj)*alsnt(k,nj))*
     .          xmchar(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*2d0*
     .          (uh-1d0))
       endif
        if (xmusnt(k).gt.1d0.and.xmustau(i).gt.1d0 )then
               xneutsf=xneutsf
     .          -2d0*g2s**2/dsntau(k)/dstau(i)*(
     .          (alsnt(k,nj)*altau(i,nj)*atau(i,ni)*antau(k,ni))*
     .          xmchar(nj)/xmneut(ni)*(uh+th-1d0-xmuneut1)+
     .          (blsnt(k,nj)*altau(i,nj)*btau(i,ni)*antau(k,ni))*
     .          (uh*th-xmuneut1)+
     .          (alsnt(k,nj)*altau(i,nj)*btau(i,ni)*antau(k,ni))*
     .          xmchar(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*
     .          (uh-1d0)+
     .          (blsnt(k,nj)*altau(i,nj)*atau(i,ni)*antau(k,ni))*
     .          dsqrt(xmutau)*sgn(ni)*(-xmuneut1+th))
      endif
             enddo
         enddo
      else
         xneutsf=0d0
      endif

c -------------------------------------------------------------------- c
c                             W+ exchange
c -------------------------------------------------------------------- c

      xmuw = mw**2/amneut(ni)**2
      dw   = y3-xmuw
      xneutwtau=0d0
      if((amchar(nj)+amtau).le.amneut(ni)) then
       if ((mw+amchar(nj)).gt.amneut(ni))then
         xneutwtau=g2s**2*16d0/dw**2*vwff**2*(
     .        (ol(ni,nj)**2+or(ni,nj)**2)*(-xmuneut1)
     .        +ol(ni,nj)**2*(xmuneut1*uh+uh*xmutau
     .         -xmutau-uh**2+uh)
     .        +or(ni,nj)**2*(xmuneut1*(th-xmutau)+th*xmutau
     .         -th**2+th)
     .        +2d0*xmchar(nj)/xmneut(ni)*ol(ni,nj)*or(ni,nj)*
     .         (uh+th-1d0-xmuneut1)
     .        +1d0/xmuw**2*xmchar(nj)/xmneut(ni)*ol(ni,nj)*
     .         or(ni,nj)*
     .         (xmuneut1**2*xmutau+xmuneut1*xmutau**2
     .          +(th+uh-1d0)*(
     .          -2d0*xmuneut1*xmutau-xmutau**2)
     .          +xmutau*(2d0*uh*th-2d0*th-2d0*uh
     .          +1d0+uh**2+th**2))
     .        +1d0/4d0/xmuw**2*(ol(ni,nj)**2+or(ni,nj)**2)*
     .         (xmuneut1**2*xmutau*(th+uh-xmutau-4d0)
     .         +xmuneut1*(xmutau*(6d0*(uh+th)-th**2-uh**2-
     .          4d0-2d0*uh*th-2d0*xmutau)
     .          +xmutau**2*(uh+th))
     .         +xmutau**2*(uh+th-1d0)
     .         +xmutau*(-th**2-uh**2+uh+th-2d0*uh*th))
     .        +2d0/xmuw*xmchar(nj)/xmneut(ni)*ol(ni,nj)*or(ni,nj)*
     .         (-xmuneut1*xmutau+xmutau*(uh+th-1d0))
     .        +1d0/xmuw*(ol(ni,nj)**2+or(ni,nj)**2)*(
     .         xmuneut1*(2d0*xmutau
     .         -xmutau*uh)-xmutau*th) )
       endif
      else
         xneutwtau=0d0
      endif
c -------------------------------------------------------------------- c
c                            H+ exchange
c -------------------------------------------------------------------- c
      dh    = y3-cmass**2/amneut(ni)**2
      xneuthtau=0d0
      if((amchar(nj)+amtau).le.amneut(ni)) then
         if ((cmass+amchar(nj)).gt.amneut(ni))then
         xneuthtau=2d0*g2s**2/dh**2*(
     .     (ql(ni,nj)**2+qr(ni,nj)**2)*(vchtau**2+achtau**2)*
     .     ((1d0+xmuneut1+xmutau)*(uh+th)-(uh+th)**2-(1d0
     .     +xmuneut1)*xmutau)+
     .     4d0*ql(ni,nj)*qr(ni,nj)*xmchar(nj)/xmneut(ni)*
     .     (vchtau**2+achtau**2)*(1d0+xmuneut1-uh-th) )
          endif
      else
         xneuthtau=0d0

      endif
c -------------------------------------------------------------------- c
c                        W+ - sfermion interference
c -------------------------------------------------------------------- c
      xneutwsf=0d0

      if((amchar(nj)+amtau).le.amneut(ni)) then
         do i=1,2
      if ((mw+amchar(nj)).gt.amneut(ni).and.
     .xmusnt(i).gt.1d0)then
            xneutwsf=xneutwsf+g2s**2*4d0/dsntau(i)/dw*vwff*
     .      (2d0*antau(i,ni)*alsnt(i,nj)*ol(ni,nj)*
     .       (xmuneut1*(+uh-1d0)+uh*xmutau
     .        -xmutau+uh*(1d0-uh))+
     .       2d0*antau(i,ni)*alsnt(i,nj)*or(ni,nj)*
     .       xmchar(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0)+
     .       4d0*antau(i,ni)*blsnt(i,nj)*ol(ni,nj)*
     .       dsqrt(xmutau)*sgn(ni)*xmchar(nj)/xmneut(ni)*(uh-1d0)+
     .       2d0*antau(i,ni)*blsnt(i,nj)*or(ni,nj)*dsqrt(xmutau)*
     .       sgn(ni)*(th-xmuneut1)+
     .       1d0/xmuw*antau(i,ni)*alsnt(i,nj)*ol(ni,nj)*
     .       (xmutau*(xmuneut1*(2d0-uh)-th))+
     .       1d0/xmuw*antau(i,ni)*alsnt(i,nj)*or(ni,nj)*
     .       xmchar(nj)/xmneut(ni)*(xmuneut1*(-xmutau)
     .       +xmutau*(th+uh-1d0))+
     .       1d0/xmuw*antau(i,ni)*blsnt(i,nj)*ol(ni,nj)*
     .       dsqrt(xmutau)*sgn(ni)*xmchar(nj)/xmneut(ni)*(xmuneut1*
     .       (-uh+1d0)+xmutau*(-
     .       uh+1d0)+uh*th-th+uh**2-2d0*uh+1d0)+
     .       1d0/xmuw*antau(i,ni)*blsnt(i,nj)*or(ni,nj)*
     .       dsqrt(xmutau)*sgn(ni)*(-xmuneut1**2+xmuneut1*(+
     .       th+2d0*uh-1d0)+xmutau*(+uh
     .       -1d0)+uh*(-th-uh+1d0)))
            endif
            if ((mw+amchar(nj)).gt.amneut(ni).and.
     .xmustau(i).gt.1d0)then
            xneutwsf=xneutwsf
     .       -g2s**2*4d0/dstau(i)/dw*vwff*(
     .       2d0*atau(i,ni)*altau(i,nj)*or(ni,nj)*(xmuneut1*(-xmutau
     .       +th-1d0)+th*(xmutau)+th*(1d0-th))
     .       +2d0*atau(i,ni)*altau(i,nj)*ol(ni,nj)*
     .       xmchar(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0)+
     .       2d0*btau(i,ni)*altau(i,nj)*ol(ni,nj)*xmchar(nj)/xmneut(ni)*
     .       dsqrt(xmutau)*sgn(ni)*(uh-1d0)+
     .       4d0*btau(i,ni)*altau(i,nj)*or(ni,nj)*dsqrt(xmutau)*sgn(ni)*
     .       (th-xmuneut1)+
     .       1d0/xmuw*atau(i,ni)*altau(i,nj)*or(ni,nj)*(xmutau*
     .       (xmuneut1*(-uh+2d0)-th))+
     .       1d0/xmuw*atau(i,ni)*altau(i,nj)*ol(ni,nj)*
     .       xmchar(nj)/xmneut(ni)*(xmuneut1*(-xmutau)
     .       +xmutau*(uh+th-1d0))+
     .       1d0/xmuw*btau(i,ni)*altau(i,nj)*ol(ni,nj)*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*(xmuneut1*
     .       (-xmutau+th-1d0)+xmutau*(th)
     .       +uh*(-th+1d0)+th*(2d0-th)-1d0)+
     .       1d0/xmuw*btau(i,ni)*altau(i,nj)*or(ni,nj)*
     .       dsqrt(xmutau)*sgn(ni)*(
     .       xmuneut1**2+xmuneut1*(xmutau-2d0*th-uh+1d0)
     .       +xmutau*(-th)+th*(th+uh-1d0)))
      endif
         enddo
      else
         xneutwsf=0d0
      endif

c -------------------------------------------------------------------- c
c                        H+ sfermion interference
c -------------------------------------------------------------------- c
      xneuthsf = 0d0

      if((amchar(nj)+amtau).le.amneut(ni)) then
         do i=1,2
       if ((cmass+amchar(nj)).gt.amneut(ni).and.
     .xmusnt(i).gt.1d0)then
            xneuthsf=xneuthsf+2d0*g2s**2/dh/dsntau(i)*(
     .       (alsnt(i,nj)*bntau(i,ni)*qr(ni,nj)*(vchtau+achtau)
     .       +antau(i,ni)*blsnt(i,nj)*ql(ni,nj)*(vchtau-achtau))*
     .       xmchar(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0)+
     .       (alsnt(i,nj)*bntau(i,ni)*ql(ni,nj)*(vchtau+achtau)
     .       +antau(i,ni)*blsnt(i,nj)*qr(ni,nj)*(vchtau-achtau))*
     .       (xmuneut1*(-uh)+uh*(-xmutau)+xmutau+uh*
     .        (th+uh-1d0))+
     .       (alsnt(i,nj)*antau(i,ni)*qr(ni,nj)*(vchtau-achtau)
     .       +bntau(i,ni)*blsnt(i,nj)*ql(ni,nj)*(vchtau+achtau))*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*(-uh+1d0)+
     .       (alsnt(i,nj)*antau(i,ni)*ql(ni,nj)*(vchtau-achtau)
     .       +bntau(i,ni)*blsnt(i,nj)*qr(ni,nj)*(vchtau+achtau))*
     .       dsqrt(xmutau)*sgn(ni)*(-xmuneut1+th))
       endif
       if ((cmass+amchar(nj)).gt.amneut(ni).and.
     .xmustau(i).gt.1d0)then
          xneuthsf=xneuthsf
     .       +2d0*g2s**2/dh/dstau(i)*(
     .       (altau(i,nj)*btau(i,ni)*ql(ni,nj)*(vchtau-achtau))*
     .       xmchar(nj)/xmneut(ni)*(uh+th-1d0-xmuneut1)+
     .       (altau(i,nj)*btau(i,ni)*qr(ni,nj)*(vchtau-achtau))*
     .       (xmuneut1*(xmutau-th)+th*(-xmutau+th+uh-1d0))+
     .       (altau(i,nj)*atau(i,ni)*ql(ni,nj)*(vchtau-achtau))*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*(uh-1d0)+
     .       (altau(i,nj)*atau(i,ni)*qr(ni,nj)*(vchtau-achtau))*
     .       dsqrt(xmutau)*sgn(ni)*(xmuneut1-th))
          endif
         enddo
      else
         xneuthsf=0d0
      endif

c -------------------------------------------------------------------- c
c                           H+ W- interference
c -------------------------------------------------------------------- c
      xneuthw=0d0
      if((amchar(nj)+amtau).le.amneut(ni)) then
         if ((cmass+amchar(nj)).gt.amneut(ni).and.
     .(mw+amchar(nj)).gt.amneut(ni))then
            xneuthw=-g2s**2/dw/dh*vwff*(
     .      +(ql(ni,nj)*or(ni,nj)+qr(ni,nj)*ol(ni,nj))*(vchtau-achtau)*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*8d0*
     .       (-1d0+uh)
     .      +(ql(ni,nj)*ol(ni,nj)+qr(ni,nj)*or(ni,nj))*(vchtau-achtau)*
     .       dsqrt(xmutau)*sgn(ni)*8d0*(xmuneut1-th)+
     .      1d0/xmuw*(
     .      (ql(ni,nj)*or(ni,nj)+qr(ni,nj)*ol(ni,nj))*(vchtau-achtau)*
     .      xmchar(nj)/xmneut(ni)*dsqrt(xmutau)*sgn(ni)*4d0*(
     .      xmuneut1*(-uh-th+xmutau+2d0)+xmutau*(-uh-th+1d0)
     .      +uh*(uh-3d0
     .      +2d0*th)+th*(th-3d0)+2d0)+
     .      (ql(ni,nj)*ol(ni,nj)+qr(ni,nj)*or(ni,nj))*(vchtau-achtau)*
     .      dsqrt(xmutau)*sgn(ni)*4d0*(xmuneut1*(-2d0*xmuneut1
     .      -xmutau+3d0*(uh+th)-2d0)+xmutau*(th+uh-1d0
     .      )+uh*(-uh+1d0
     .      -2d0*th)+th*(1d0-th)) ) )
      endif
      else
         xneuthw=0d0
      endif


      NS_chtauntau = xneutsf+xneutwtau+xneuthtau+xneutwsf+xneuthsf+
     .               xneuthw

      end
c ==================================================================== c
c =======================  chargino- up downbar ====================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_chubd(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i
      DOUBLE PRECISION aup(2,5),bup(2,5),ado(2,5),bdo(2,5)
      DOUBLE PRECISION alup(2,2),aldo(2,2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION dsd(2),dsu(2)
      DOUBLE PRECISION ql(5,2),qr(5,2),or(5,2),ol(5,2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION awff,vwff
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION xmuneut1,xmusd(2),xmusu(2),
     .xneutsf,xmuw,dw,xneutw,xneutwsf,x1,x2,x3,y1,y2,y3
*
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_indices/ni,nj
      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/NS_coup7/alup,aldo
      COMMON/NS_coup10/aup,bup,ado,bdo
      COMMON/NS_coup18/awff,vwff
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar

      xmuneut1 = amchar(nj)**2/amneut(ni)**2

      x3 = 2d0-x1-x2
      y1 = 1d0-xmuneut1-x1
      y2 = 1d0-xmuneut1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                          sfermion exchange
c -------------------------------------------------------------------- c
      xmusd(1) = asdown1**2/amneut(ni)**2
      xmusd(2) = asdown2**2/amneut(ni)**2
      xmusu(1) = asup1**2/amneut(ni)**2
      xmusu(2) = asup2**2/amneut(ni)**2

      dsd(1) = 1d0-x2-xmusd(1)
      dsd(2) = 1d0-x2-xmusd(2)
      dsu(1) = 1d0-x1-xmusu(1)
      dsu(2) = 1d0-x1-xmusu(2)

      xneutsf=0d0

      if(amchar(nj).le.amneut(ni)) then
         do i=1,2
            do k=1,2
               if (xmusd(i).gt.1d0.and.xmusd(k).gt.1d0)then
               xneutsf=xneutsf
     .              +g2s**2/dsd(k)/dsd(i)*x2*y2*
     .              (ado(i,ni)*ado(k,ni)+bdo(i,ni)*bdo(k,ni))*
     .              (aldo(i,nj)*aldo(k,nj))
               endif
               if (xmusu(i).gt.1d0.and.xmusu(k).gt.1d0)then
               xneutsf=xneutsf
     .              +g2s**2/dsu(k)/dsu(i)*x1*y1*
     .              (aup(i,ni)*aup(k,ni)+bup(i,ni)*bup(k,ni))*
     .              (alup(i,nj)*alup(k,nj))
               endif
               if (xmusu(i).gt.1d0.and.xmusd(k).gt.1d0)then
               xneutsf=xneutsf
     .              +g2s**2/dsu(i)/dsd(k)*
     .              (2d0*xmchar(nj)/xmneut(ni)*y3*
     .              ( ado(k,ni)*aup(i,ni)*alup(i,nj)*aldo(k,nj)))
               endif
            enddo
         enddo
      else
         xneutsf=0d0
      endif

c -------------------------------------------------------------------- c
c                             W+ exchange
c -------------------------------------------------------------------- c
      xmuw = mw**2/amneut(ni)**2
      dw   = y3-xmuw
      xneutw=0d0
      if(amchar(nj).le.amneut(ni)) then
         if ((mw+amchar(nj)).gt.amneut(ni))then
         xneutw=g2s**2*4d0/dw**2*
     .        (4d0*vwff**2*or(ni,nj)**2*x2*y2
     .        +4d0*vwff**2*ol(ni,nj)**2*x1*y1
     .        -8d0*xmchar(nj)/xmneut(ni)*
     .        ol(ni,nj)*or(ni,nj)*vwff**2*y3 )
         endif
      else
         xneutw=0d0
      endif
c -------------------------------------------------------------------- c
c                        W+ - sfermion interference
c -------------------------------------------------------------------- c
      xneutwsf=0d0

      if(amchar(nj).le.amneut(ni)) then
         do i=1,2
            if ((mw+amchar(nj)).gt.amneut(ni).and.xmusd(i).gt.1d0)then
            xneutwsf=xneutwsf-g2s**2*4d0/dsd(i)/dw
     .           *(ado(i,ni)*aldo(i,nj)*or(ni,nj)*2d0*vwff*x2*y2
     .            -aldo(i,nj)*ado(i,ni)*ol(ni,nj)*2d0*vwff*
     .            xmchar(nj)/xmneut(ni)*y3)
            endif
            if ((mw+amchar(nj)).gt.amneut(ni).and.xmusu(i).gt.1d0)then
               xneutwsf=xneutwsf
     .           +g2s**2*4d0/dsu(i)/dw
     .           *(aup(i,ni)*alup(i,nj)*ol(ni,nj)*2d0*vwff*x1*y1
     .            -alup(i,nj)*aup(i,ni)*or(ni,nj)*2d0*vwff*
     .            xmchar(nj)/xmneut(ni)*y3)
            endif
         enddo
      else
         xneutwsf=0d0
      endif

      NS_chubd = xneutsf+xneutw+xneutwsf

      end
c ==================================================================== c
c =====================  chargino- top bottombar ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_chtbb(x1,x2)
*
      IMPLICIT NONE
      INTEGER k,ni,nj,i
      DOUBLE PRECISION atopr(2,5),btopr(2,5)
      DOUBLE PRECISION abot(2,5),bbot(2,5)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2),alstor(2,2),akstor(2,2)
      DOUBLE PRECISION ql(5,2),qr(5,2),or(5,2),ol(5,2)
      DOUBLE PRECISION dsbot(2),dstop(2),sgn(5)
      DOUBLE PRECISION vchtopr,achtopr,achtau,vchtau
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION awff,vwff
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION xmuneut1,xmut,xmub,uh,th,x1,x2,x3,y3,
     .xmusbot(2),xmustop(2),xneutsf,xneuthsf,xmuw,dw,
     .xneutw,dh,xneuth,xneutwsf,xneuthw,gmsb(2),gmst(2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
*
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/NS_coup15/achtopr,vchtopr,achtau,vchtau
      COMMON/NS_coup18/awff,vwff
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_indices/ni,nj
      COMMON/NS_neutstoptop/atopr,btopr
      COMMON/NS_neutsbotbot/abot,bbot
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor

      do i=1,5
         sgn(i) = 1d0
         if(xmneut(i).gt.0d0) then
            sgn(i) = 1d0
         elseif(xmneut(i).lt.0d0) then
            sgn(i) = -1d0
         endif
      enddo

      xmuneut1 = amchar(nj)**2/amneut(ni)**2
      xmut     = amt**2/amneut(ni)**2
      xmub     = amb**2/amneut(ni)**2

      uh = 1d0-x1+xmut
      th = 1d0-x2+xmub

      x3 = 2d0-x1-x2
      y3 = 1d0+xmuneut1-x3
c -------------------------------------------------------------------- c
c                          sfermion exchange
c -------------------------------------------------------------------- c
      gmsb(1)=asb1
      gmsb(2)=asb2
      gmst(1)=ast1
      gmst(2)=ast2
      xmusbot(1) = asb1**2/amneut(ni)**2
      xmusbot(2) = asb2**2/amneut(ni)**2
      xmustop(1) = ast1**2/amneut(ni)**2
      xmustop(2) = ast2**2/amneut(ni)**2

      dsbot(1) = 1-x2-xmusbot(1)+xmub
      dsbot(2) = 1-x2-xmusbot(2)+xmub
      dstop(1) = 1-x1-xmustop(1)+xmut
      dstop(2) = 1-x1-xmustop(2)+xmut

      xneutsf=0d0

      if((amchar(nj)+amt+amb).le.amneut(ni)) then
         do i=1,2
            do k=1,2
       if ((gmsb(i)+amt).gt.amneut(ni).and.(
     .gmsb(k)+amt).gt.amneut(ni))then
               xneutsf=xneutsf
     .          +g2s**2/dsbot(k)/dsbot(i)*(
     .          (abot(i,ni)*abot(k,ni)+bbot(i,ni)*bbot(k,ni))*
     .          (alsbot(i,nj)*alsbot(k,nj)+aksbot(i,nj)*aksbot(k,nj))*
     .          (xmuneut1*(th-xmub-1d0)+xmut*(th-1d0)+xmub*(-xmut
     .          +th)+th*(-th+1d0))+
     .          (abot(i,ni)*abot(k,ni)+bbot(i,ni)*bbot(k,ni))*
     .          (alsbot(i,nj)*aksbot(k,nj)+aksbot(i,nj)*alsbot(k,nj))*
     .          xmchar(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*2d0*
     .          (th-xmub-1d0)
     .          +(abot(i,ni)*bbot(k,ni)+bbot(i,ni)*abot(k,ni))*
     .          (alsbot(i,nj)*alsbot(k,nj)+aksbot(i,nj)*aksbot(k,nj))*
     .          dsqrt(xmub)*sgn(ni)*2d0*(th-xmut-xmuneut1)+
     .          (abot(i,ni)*bbot(k,ni)+bbot(i,ni)*abot(k,ni))*
     .          (alsbot(i,nj)*aksbot(k,nj)+aksbot(i,nj)*alsbot(k,nj))*
     .          xmchar(nj)/xmneut(ni)*dsqrt(xmut*xmub)*(-4d0) )
       endif
       if (xmustop(i).gt.1d0.and.xmustop(k).gt.1d0)then
               xneutsf=xneutsf
     .          +g2s**2/dstop(k)/dstop(i)*(
     .          (atopr(i,ni)*atopr(k,ni)+btopr(i,ni)*btopr(k,ni))*
     .          (alstor(i,nj)*alstor(k,nj)+akstor(i,nj)*akstor(k,nj))*
     .          (xmuneut1*(-xmut+uh-1d0)+xmut*(uh-xmub)+xmub*
     .          (uh-1d0)+uh*(-uh+1d0))+
     .          (atopr(i,ni)*atopr(k,ni)+btopr(i,ni)*btopr(k,ni))*
     .          (alstor(i,nj)*akstor(k,nj)+akstor(i,nj)*alstor(k,nj))*
     .          xmchar(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*2d0*
     .          (uh-xmut-1d0)
     .          +(atopr(i,ni)*btopr(k,ni)+btopr(i,ni)*atopr(k,ni))*
     .          (alstor(i,nj)*alstor(k,nj)+akstor(i,nj)*akstor(k,nj))*
     .          dsqrt(xmut)*sgn(ni)*2d0*(uh-xmub-xmuneut1)+
     .          (atopr(i,ni)*btopr(k,ni)+btopr(i,ni)*atopr(k,ni))*
     .          (alstor(i,nj)*akstor(k,nj)+akstor(i,nj)*alstor(k,nj))*
     .          dsqrt(xmut*xmub)*xmchar(nj)/xmneut(ni)*(-4d0))
       endif
       if ((gmsb(i)+amt).gt.amneut(ni).and.xmustop(k).gt.1d0)then
               xneutsf=xneutsf
     .          -2d0*g2s**2/dstop(k)/dsbot(i)*(
     .          (alstor(k,nj)*alsbot(i,nj)*abot(i,ni)*atopr(k,ni)+
     .           akstor(k,nj)*aksbot(i,nj)*bbot(i,ni)*btopr(k,ni))*
     .          xmchar(nj)/xmneut(ni)*(uh+th-1d0-xmuneut1)+
     .          (alstor(k,nj)*aksbot(i,nj)*abot(i,ni)*btopr(k,ni)+
     .           akstor(k,nj)*alsbot(i,nj)*bbot(i,ni)*atopr(k,ni))*
     .          (uh*th-xmut*xmub-xmuneut1)+
     .          (alstor(k,nj)*alsbot(i,nj)*abot(i,ni)*btopr(k,ni)+
     .           akstor(k,nj)*aksbot(i,nj)*bbot(i,ni)*atopr(k,ni))*
     .          xmchar(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*
     .          (th-xmub-1d0)+
     .          (alstor(k,nj)*alsbot(i,nj)*bbot(i,ni)*atopr(k,ni)+
     .           akstor(k,nj)*aksbot(i,nj)*abot(i,ni)*btopr(k,ni))*
     .          xmchar(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*
     .          (uh-xmut-1d0)+
     .          (alstor(k,nj)*aksbot(i,nj)*abot(i,ni)*atopr(k,ni)+
     .           akstor(k,nj)*alsbot(i,nj)*bbot(i,ni)*btopr(k,ni))*
     .          dsqrt(xmut)*sgn(ni)*(-xmuneut1-xmub+uh)+
     .          (alstor(k,nj)*aksbot(i,nj)*bbot(i,ni)*btopr(k,ni)+
     .           akstor(k,nj)*alsbot(i,nj)*abot(i,ni)*atopr(k,ni))*
     .          dsqrt(xmub)*sgn(ni)*(-xmuneut1-xmut-th)+
     .          (alstor(k,nj)*alsbot(i,nj)*bbot(i,ni)*btopr(k,ni)+
     .           akstor(k,nj)*aksbot(i,nj)*abot(i,ni)*atopr(k,ni))*
     .          xmchar(nj)/xmneut(ni)*dsqrt(xmut*xmub)*(-2d0)+
     .          (alstor(k,nj)*aksbot(i,nj)*bbot(i,ni)*atopr(k,ni)+
     .           akstor(k,nj)*alsbot(i,nj)*abot(i,ni)*btopr(k,ni))*
     .          dsqrt(xmut*xmub)*(uh+th-xmut-xmub) )
       endif
             enddo
         enddo
      else
         xneutsf=0d0
      endif

c -------------------------------------------------------------------- c
c                               W+ exchange
c -------------------------------------------------------------------- c
      xmuw = mw**2/amneut(ni)**2
      dw   = y3-xmuw
      xneutw=0d0
      if((amchar(nj)+amt+amb).le.amneut(ni)) then
         if ((mw+amchar(nj)).gt.amneut(ni))then
         xneutw=g2s**2*16d0/dw**2*vwff**2*(
     .        (ol(ni,nj)**2+or(ni,nj)**2)*(-xmuneut1-xmub*xmut)
     .        +ol(ni,nj)**2*(xmuneut1*(uh-xmut)+uh*(xmub+xmut)
     .         -xmub-uh**2+uh)
     .        +or(ni,nj)**2*(xmuneut1*(th-xmub)+th*(xmub+xmut)
     .         -xmut-th**2+th)
     .        +2d0*xmchar(nj)/xmneut(ni)*ol(ni,nj)*or(ni,nj)*
     .         (uh+th-1d0-xmuneut1)
     .        +1d0/xmuw**2*xmchar(nj)/xmneut(ni)*ol(ni,nj)*
     .         or(ni,nj)*
     .         (xmuneut1**2*(xmub+xmut)+xmuneut1*(xmub**2+xmut**2
     .          +6d0*xmub*xmut)+(th+uh-1d0)*(-2d0*xmuneut1*xmut
     .          -2d0*xmuneut1*xmub-xmub**2-xmut**2-6d0*xmut*
     .          xmub)+(xmut+xmub)*(2d0*uh*th-2d0*th-2d0*uh
     .          +1d0+uh**2+th**2+4d0*xmub*xmut))
     .        +1d0/4d0/xmuw**2*(ol(ni,nj)**2+or(ni,nj)**2)*
     .         (xmuneut1**2*(xmut+xmub)*(th+uh-xmut-xmub-4d0)
     .         +xmuneut1*((xmub+xmut)*(6d0*(uh+th)-th**2-uh**2-
     .          4d0-2d0*uh*th-2d0*xmub-2d0*xmut)+6d0*xmub*xmut*
     .          (uh+th-2d0)-4d0*xmub*xmut*(xmub+xmut+1d0)
     .          +(xmub**2+xmut**2)*(uh+th))
     .         +xmut**2*(-4d0*xmub+uh+th-1d0)
     .         +xmub**2*(-4d0*xmut+uh+th-1d0)
     .         +xmub*xmut*(6d0*th+6d0*uh-2d0)
     .         +(xmut+xmub)*(-th**2-uh**2+uh+th-2d0*uh*th))
     .        +2d0/xmuw*xmchar(nj)/xmneut(ni)*ol(ni,nj)*or(ni,nj)*
     .         (-xmuneut1*(xmub+xmut)-4d0*xmub*xmut+
     .          (xmub+xmut)*(uh+th-1d0))
     .        +1d0/xmuw*(ol(ni,nj)**2+or(ni,nj)**2)*(
     .         xmuneut1*(2d0*xmub*xmut+2d0*xmub+2d0*xmut
     .         -xmut*th-xmub*uh)-xmub*th-xmut*uh+2d0*xmub*xmut) )
         endif
      else
         xneutw=0d0
      endif
c -------------------------------------------------------------------- c
c                              H+ exchange
c -------------------------------------------------------------------- c
      dh    = y3-cmass**2/amneut(ni)**2
      xneuth=0d0
      if((amchar(nj)+amt+amb).le.amneut(ni)) then
         if ((cmass+amchar(nj)).gt.amneut(ni))then
         xneuth=2d0*g2s**2/dh**2*(
     .     (ql(ni,nj)**2+qr(ni,nj)**2)*(vchtopr**2+achtopr**2)*
     .     ((1d0+xmuneut1+xmut+xmub)*(uh+th)-(uh+th)**2-(1d0
     .     +xmuneut1)*(xmut+xmub))+
     .     (ql(ni,nj)**2+qr(ni,nj)**2)*(vchtopr**2-achtopr**2)*(-2d0)*
     .     dsqrt(xmut*xmub)*(uh+th-xmut-xmub)+
     .     4d0*ql(ni,nj)*qr(ni,nj)*xmchar(nj)/xmneut(ni)*
     .     (vchtopr**2+achtopr**2)*(1d0+xmuneut1-uh-th)+
     .     4d0*ql(ni,nj)*qr(ni,nj)*xmchar(nj)/xmneut(ni)*
     .     (vchtopr**2-achtopr**2)*(-2d0)*dsqrt(xmut*xmub) )
         endif
      else
         xneuth=0d0
      endif

c -------------------------------------------------------------------- c
c                         W+ - sfermion interference
c -------------------------------------------------------------------- c
      xneutwsf=0d0

      if((amchar(nj)+amt+amb).le.amneut(ni)) then
         do i=1,2
            if ((mw+amchar(nj)).gt.amneut(ni).and.xmustop(i).gt.1d0)
     .then
            xneutwsf=xneutwsf+g2s**2*4d0/dstop(i)/dw*vwff*
     .      (2d0*atopr(i,ni)*alstor(i,nj)*ol(ni,nj)*
     .       (xmuneut1*(-xmut+uh-1d0)+uh*(xmut+xmub)-xmut*xmub
     .        -xmub+uh*(1d0-uh))+
     .       2d0*atopr(i,ni)*alstor(i,nj)*or(ni,nj)*
     .       xmchar(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0)+
     .       2d0*btopr(i,ni)*alstor(i,nj)*or(ni,nj)*
     .       dsqrt(xmut)*sgn(ni)*xmchar(nj)/xmneut(ni)*(th-xmub-1d0)+
     .       4d0*atopr(i,ni)*akstor(i,nj)*ol(ni,nj)*
     .       dsqrt(xmub)*sgn(ni)*xmchar(nj)/xmneut(ni)*(uh-xmut-1d0)+
     .       4d0*btopr(i,ni)*alstor(i,nj)*ol(ni,nj)*
     .       dsqrt(xmut)*sgn(ni)*(uh-xmuneut1-xmub)+
     .       2d0*atopr(i,ni)*akstor(i,nj)*or(ni,nj)*
     .       dsqrt(xmub)*sgn(ni)*(th-xmuneut1-xmut)+
     .       2d0*btopr(i,ni)*akstor(i,nj)*or(ni,nj)*dsqrt(xmut*xmub)*
     .       (uh+th-xmut-xmub)+
     .       8d0*btopr(i,ni)*akstor(i,nj)*ol(ni,nj)*(-1d0)*
     .       dsqrt(xmut*xmub)*xmchar(nj)/xmneut(ni)+
     .       1d0/xmuw*atopr(i,ni)*alstor(i,nj)*ol(ni,nj)*
     .       (xmut*(xmuneut1*(xmub-th+2d0)+xmub-uh)+
     .        xmub*(xmuneut1*(xmut-uh+2d0)+xmut-th))+
     .       1d0/xmuw*atopr(i,ni)*alstor(i,nj)*or(ni,nj)*
     .       xmchar(nj)/xmneut(ni)*(xmuneut1*(-xmut-xmub)+xmut*(th+uh
     .       -1d0-2d0*xmub)+xmub*(th+uh-1d0-2d0*xmut))+
     .       1d0/xmuw*btopr(i,ni)*alstor(i,nj)*or(ni,nj)*
     .       dsqrt(xmut)*sgn(ni)*
     .       xmchar(nj)/xmneut(ni)*(xmuneut1*(-xmut+uh-1d0)+xmut*
     .       (-xmub+uh)+xmub*(-xmub+th+2d0*uh-3d0)+th*(1d0-uh)
     .       +uh*(2d0-uh)-1d0)+
     .       1d0/xmuw*atopr(i,ni)*akstor(i,nj)*ol(ni,nj)*
     .       dsqrt(xmub)*sgn(ni)*
     .       xmchar(nj)/xmneut(ni)*(xmuneut1*(xmut-uh+1d0)+xmut*(
     .       xmut-th-2d0*uh+2d0)+xmub*(xmut-uh+1d0)+uh*th-th+uh**2
     .       -2d0*uh+1d0)+
     .       1d0/xmuw*btopr(i,ni)*alstor(i,nj)*ol(ni,nj)*
     .       dsqrt(xmut)*sgn(ni)*(xmuneut1**2+xmuneut1*(xmut+2d0*xmub
     .       -th-2d0*uh+1d0)+xmut*(xmub-uh)+xmub*(xmub-th-2d0*uh+
     .       1d0)+uh*(th+uh-1d0))+
     .       1d0/xmuw*atopr(i,ni)*akstor(i,nj)*or(ni,nj)*
     .       dsqrt(xmub)*sgn(ni)*(-xmuneut1**2+xmuneut1*(-3d0*xmut+th+
     .       2d0*uh-1d0)+xmut*(-xmut+th+2d0*uh)+xmub*(-xmut+uh-1d0)
     .       +uh*(-th-uh+1d0))+
     .       2d0/xmuw*btopr(i,ni)*akstor(i,nj)*ol(ni,nj)*
     .       dsqrt(xmut*xmub)*xmchar(nj)/xmneut(ni)*(xmuneut1+xmut+
     .       xmub-uh-th+1d0)+
     .       1d0/xmuw*btopr(i,ni)*akstor(i,nj)*or(ni,nj)*
     .       dsqrt(xmut*xmub)*(xmuneut1*(-xmut-xmub+th+uh-4d0)
     .       -xmut-xmub+th+uh))
            endif
        if ((mw+amchar(nj)).gt.amneut(ni).and.
     .(gmsb(i)+amt).gt.amneut(ni))then
            xneutwsf=xneutwsf
     .       -g2s**2*4d0/dsbot(i)/dw*vwff*(
     .       2d0*abot(i,ni)*alsbot(i,nj)*or(ni,nj)*(xmuneut1*(-xmub
     .       +th-1d0)+th*(xmut+xmub)-xmut*xmub-xmut+th*(1d0-th))
     .       +2d0*abot(i,ni)*alsbot(i,nj)*ol(ni,nj)*
     .       xmchar(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0)+
     .       4d0*abot(i,ni)*aksbot(i,nj)*or(ni,nj)*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*(th-xmub-1d0)+
     .       2d0*bbot(i,ni)*alsbot(i,nj)*ol(ni,nj)*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*(uh-xmut-1d0)+
     .       2d0*abot(i,ni)*aksbot(i,nj)*ol(ni,nj)*dsqrt(xmut)*sgn(ni)*
     .       (uh-xmub-xmuneut1)+
     .       4d0*bbot(i,ni)*alsbot(i,nj)*or(ni,nj)*dsqrt(xmub)*sgn(ni)*
     .       (th-xmut-xmuneut1)+
     .       2d0*bbot(i,ni)*aksbot(i,nj)*ol(ni,nj)*dsqrt(xmub*xmut)*
     .       (uh+th-xmut-xmub)+
     .       8d0*bbot(i,ni)*aksbot(i,nj)*or(ni,nj)*dsqrt(xmub*xmut)*
     .       xmchar(nj)/xmneut(ni)*(-1d0)+
     .       1d0/xmuw*abot(i,ni)*alsbot(i,nj)*or(ni,nj)*(xmub*
     .       (xmuneut1*(xmut-uh+2d0)+xmut-th)+xmut*(xmuneut1*
     .       (xmub-th+2d0)+xmub-uh))+
     .       1d0/xmuw*abot(i,ni)*alsbot(i,nj)*ol(ni,nj)*
     .       xmchar(nj)/xmneut(ni)*(xmuneut1*(-xmut-xmub)+xmut*(uh
     .       +th-1d0-2d0*xmub)+xmub*(uh+th-1d0-2d0*xmut))+
     .       1d0/xmuw*abot(i,ni)*aksbot(i,nj)*or(ni,nj)*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*(xmuneut1*
     .       (xmub-th+1d0)+xmub*(xmub-2d0*th-uh+2d0)+xmut*(xmub-th
     .       +1d0)+uh*th-uh+th**2-2d0*th+1d0)+
     .       1d0/xmuw*bbot(i,ni)*alsbot(i,nj)*ol(ni,nj)*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*(xmuneut1*(-xmub+
     .       th-1d0)+xmub*(-xmut+th)+xmut*(-xmut+2d0*th+uh-3d0)+
     .       uh*(-th+1d0)+th*(2d0-th)-1d0)+
     .       1d0/xmuw*abot(i,ni)*aksbot(i,nj)*ol(ni,nj)*
     .       dsqrt(xmut)*sgn(ni)*(-xmuneut1**2+xmuneut1*(-3d0*xmub+
     .       2d0*th+uh-1d0)+xmub*(-xmub+2d0*th+uh)+xmut*(-xmub+th-
     .       1d0)+th*(-th-uh+1d0))+
     .       1d0/xmuw*bbot(i,ni)*alsbot(i,nj)*or(ni,nj)*
     .       dsqrt(xmub)*sgn(ni)*(xmuneut1**2+xmuneut1*(2d0*xmut+xmub
     .       -2d0*th-uh+1d0)+xmub*(xmut-th)+xmut*(xmut-2d0*th-uh+
     .       1d0)+th*(th+uh-1d0))+
     .       1d0/xmuw*bbot(i,ni)*aksbot(i,nj)*ol(ni,nj)*
     .       dsqrt(xmut*xmub)*(xmuneut1*(-xmut-xmub+th+uh-4d0)+uh
     .       +th-xmut-xmub)+
     .       2d0/xmuw*bbot(i,ni)*aksbot(i,nj)*or(ni,nj)*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmut*xmub)*(xmuneut1+xmut
     .       +xmub-th-uh+1d0))
            endif
         enddo
      else
         xneutwsf=0d0
      endif

c -------------------------------------------------------------------- c
c                          H+ sfermion interference
c -------------------------------------------------------------------- c
      xneuthsf = 0d0


      if((amchar(nj)+amt+amb).le.amneut(ni)) then
         do i=1,2
      if((cmass+amchar(nj)).gt.amneut(ni).and.xmustop(i).gt.1d0)then
            xneuthsf=xneuthsf+2d0*g2s**2/dh/dstop(i)*(
     .       (alstor(i,nj)*btopr(i,ni)*qr(ni,nj)*(vchtopr+achtopr)
     .       +atopr(i,ni)*akstor(i,nj)*ql(ni,nj)*(vchtopr-achtopr))*
     .       xmchar(nj)/xmneut(ni)*(uh+th-xmuneut1-1d0)+
     .       (alstor(i,nj)*btopr(i,ni)*ql(ni,nj)*(vchtopr+achtopr)
     .       +atopr(i,ni)*akstor(i,nj)*qr(ni,nj)*(vchtopr-achtopr))*
     .       (xmuneut1*(xmut-uh)+uh*(-xmut-xmub)+xmub+uh*
     .        (th+uh-1d0))+
     .       (alstor(i,nj)*atopr(i,ni)*qr(ni,nj)*(vchtopr+achtopr)
     .       +btopr(i,ni)*akstor(i,nj)*ql(ni,nj)*(vchtopr-achtopr))*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*(-xmub+th-1d0)+
     .       (alstor(i,nj)*atopr(i,ni)*qr(ni,nj)*(vchtopr-achtopr)
     .       +btopr(i,ni)*akstor(i,nj)*ql(ni,nj)*(vchtopr+achtopr))*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*(xmut-uh+1d0)+
     .       (alstor(i,nj)*atopr(i,ni)*ql(ni,nj)*(vchtopr+achtopr)
     .       +btopr(i,ni)*akstor(i,nj)*qr(ni,nj)*(vchtopr-achtopr))*
     .       dsqrt(xmut)*sgn(ni)*(xmuneut1+xmub-uh)+
     .       (alstor(i,nj)*atopr(i,ni)*ql(ni,nj)*(vchtopr-achtopr)
     .       +btopr(i,ni)*akstor(i,nj)*qr(ni,nj)*(vchtopr+achtopr))*
     .       dsqrt(xmub)*sgn(ni)*(-xmuneut1-xmut+th)+
     .       (alstor(i,nj)*btopr(i,ni)*qr(ni,nj)*(vchtopr-achtopr)
     .       +atopr(i,ni)*akstor(i,nj)*ql(ni,nj)*(vchtopr+achtopr))*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmut*xmub)*2d0+
     .       (alstor(i,nj)*btopr(i,ni)*ql(ni,nj)*(vchtopr-achtopr)
     .       +atopr(i,ni)*akstor(i,nj)*qr(ni,nj)*(vchtopr+achtopr))*
     .       dsqrt(xmut*xmub)*(-xmut-xmub+uh+th) )
            endif
       if((cmass+amchar(nj)).gt.amneut(ni).and.
     .(gmsb(i)+amt).gt.amneut(ni))then
            xneuthsf=xneuthsf
     .       +2d0*g2s**2/dh/dsbot(i)*(
     .       (alsbot(i,nj)*bbot(i,ni)*ql(ni,nj)*(vchtopr-achtopr)
     .       +aksbot(i,nj)*abot(i,ni)*qr(ni,nj)*(vchtopr+achtopr))*
     .       xmchar(nj)/xmneut(ni)*(uh+th-1d0-xmuneut1)+
     .       (alsbot(i,nj)*bbot(i,ni)*qr(ni,nj)*(vchtopr-achtopr)
     .       +aksbot(i,nj)*abot(i,ni)*ql(ni,nj)*(vchtopr+achtopr))*
     .       (xmuneut1*(xmub-th)+xmut*(1d0-th)+th*
     .       (-xmub+th+uh-1d0))+
     .       (alsbot(i,nj)*abot(i,ni)*ql(ni,nj)*(vchtopr+achtopr)
     .       +aksbot(i,nj)*bbot(i,ni)*qr(ni,nj)*(vchtopr-achtopr))*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*(xmub-th+1d0)+
     .       (alsbot(i,nj)*abot(i,ni)*ql(ni,nj)*(vchtopr-achtopr)
     .       +aksbot(i,nj)*bbot(i,ni)*qr(ni,nj)*(vchtopr+achtopr))*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*(uh-xmut-1d0)+
     .       (alsbot(i,nj)*abot(i,ni)*qr(ni,nj)*(vchtopr+achtopr)
     .       +aksbot(i,nj)*bbot(i,ni)*ql(ni,nj)*(vchtopr-achtopr))*
     .       dsqrt(xmut)*sgn(ni)*(uh-xmuneut1-xmub)+
     .       (alsbot(i,nj)*abot(i,ni)*qr(ni,nj)*(vchtopr-achtopr)
     .       +aksbot(i,nj)*bbot(i,ni)*ql(ni,nj)*(vchtopr+achtopr))*
     .       dsqrt(xmub)*sgn(ni)*(xmut+xmuneut1-th)+
     .       (alsbot(i,nj)*bbot(i,ni)*ql(ni,nj)*(vchtopr+achtopr)
     .       +aksbot(i,nj)*abot(i,ni)*qr(ni,nj)*(vchtopr-achtopr))*
     .       xmchar(nj)/xmneut(ni)*dsqrt(xmut*xmub)*2d0+
     .       (alsbot(i,nj)*bbot(i,ni)*qr(ni,nj)*(vchtopr+achtopr)
     .       +aksbot(i,nj)*abot(i,ni)*ql(ni,nj)*(vchtopr-achtopr))*
     .       dsqrt(xmut*xmub)*(uh+th-xmut-xmub) )
      endif
         enddo
      else
         xneuthsf=0d0
      endif

c -------------------------------------------------------------------- c
c                           H+ W- interference
c -------------------------------------------------------------------- c
      xneuthw=0d0
      if((amchar(nj)+amt+amb).le.amneut(ni)) then
      if ((cmass+amchar(nj)).gt.amneut(ni).and.
     .(mw+amchar(nj)).gt.amneut(ni))then
            xneuthw=-g2s**2/dw/dh*vwff*(
     .     (ql(ni,nj)*or(ni,nj)+qr(ni,nj)*ol(ni,nj))*(vchtopr+achtopr)*
     .     xmchar(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*8d0*(1d0-th+xmub)
     .     +(ql(ni,nj)*or(ni,nj)+qr(ni,nj)*ol(ni,nj))*(vchtopr-achtopr)*
     .     xmchar(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*8d0*(-1d0+uh
     .     -xmut)
     .     +(ql(ni,nj)*ol(ni,nj)+qr(ni,nj)*or(ni,nj))*(vchtopr+achtopr)*
     .     dsqrt(xmut)*sgn(ni)*8d0*(uh-xmuneut1-xmub)
     .     +(ql(ni,nj)*ol(ni,nj)+qr(ni,nj)*or(ni,nj))*(vchtopr-achtopr)*
     .     dsqrt(xmub)*sgn(ni)*8d0*(xmut+xmuneut1-th)+
     .     1d0/xmuw*(
     .     (ql(ni,nj)*or(ni,nj)+qr(ni,nj)*ol(ni,nj))*(vchtopr+achtopr)*
     .     xmchar(nj)/xmneut(ni)*dsqrt(xmut)*sgn(ni)*4d0*(
     .     xmuneut1*(uh+th-xmut-xmub-2d0)+xmut*(uh+th-2d0*xmub
     .     -1d0)+xmub*(-2d0*xmub+3d0*(th+uh)-5d0)+th*(3d0
     .     -2d0*uh-th)+uh*(3d0-uh)-2d0)+
     .     (ql(ni,nj)*or(ni,nj)+qr(ni,nj)*ol(ni,nj))*(vchtopr-achtopr)*
     .     xmchar(nj)/xmneut(ni)*dsqrt(xmub)*sgn(ni)*4d0*(
     .     xmuneut1*(-uh-th+xmut+xmub+2d0)+xmub*(-uh-th+1d0+
     .     2d0*xmut)+xmut*(2d0*xmut-3d0*(th+uh)+5d0)+uh*(uh-3d0
     .     +2d0*th)+th*(th-3d0)+2d0)+
     .     (ql(ni,nj)*ol(ni,nj)+qr(ni,nj)*or(ni,nj))*(vchtopr+achtopr)*
     .     dsqrt(xmut)*sgn(ni)*4d0*(xmuneut1*(2d0*xmuneut1+xmut+
     .     5d0*xmub-3d0*(th+uh)+2d0)+xmut*(2d0*xmub-th-uh+1d0)+
     .     xmub*(2d0*xmub-3d0*(uh+th)+1d0)+th*(th+2d0*uh-1d0)+uh*(
     .     uh-1d0))+
     .     (ql(ni,nj)*ol(ni,nj)+qr(ni,nj)*or(ni,nj))*(vchtopr-achtopr)*
     .     dsqrt(xmub)*sgn(ni)*4d0*(xmuneut1*(-2d0*xmuneut1-5d0*xmut
     .     -xmub+3d0*(uh+th)-2d0)+xmub*(th+uh-1d0-2d0*xmut)+
     .     xmut*(-2d0*xmut+3d0*(th+uh)-1d0)+uh*(-uh+1d0-2d0*th)+
     .     th*(1d0-th))) )
      endif
      else
         xneuthw=0d0
      endif

      NS_chtbb = xneutsf+xneutw+xneuth+xneutwsf+xneuthsf+xneuthw

      end
c-----------------------------------------------------------------------
