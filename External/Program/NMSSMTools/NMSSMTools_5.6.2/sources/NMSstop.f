      SUBROUTINE NS_STOP

************************************************************************
*
*     This subroutine computes the top squark decays
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

      INTEGER I,J,jsign,nj,k,GRFLAG

      DOUBLE PRECISION amuv,lamv,amuvdiv
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION alp, nf
      DOUBLE PRECISION amsq,scalmur
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .         cmm,smm,cum,sum,cdm,sdm,cem,sem,cnm,snm
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION rmtc,rmbc,rmtauc
      DOUBLE PRECISION amurefer
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION atopr(2,5),btopr(2,5)
      DOUBLE PRECISION alstor(2,2),akstor(2,2)
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2),gzmumu(2,2)
      DOUBLE PRECISION gwtb(2,2),gwntau(2,2),gwnmu(2,2)
      DOUBLE PRECISION gctbr(2,2)
      DOUBLE PRECISION Hstaustaur(3,2,2),Astaustaur(2,2,2)
      DOUBLE PRECISION Hsbotsbotr(3,2,2),Asbotsbotr(2,2,2)
      DOUBLE PRECISION Hstopstopr(3,2,2),Astopstopr(3,2,2)
      DOUBLE PRECISION gmsb(2)
      DOUBLE PRECISION delta11c,delta12c,delta13c,delta14c,delta15c,
     .         delta1H(3),delta2H(3),delta3H(3),delta4H(3),delta5H(3),
     .         delta1A(2),delta2A(2),delta3A(2),delta4A(2),delta5A(2),
     .         delta21c,delta22c,delta23c,delta24c,delta25c,
     .         adel1,adel2,adel3,adel4,adel5,
     .         bdel1,bdel2,bdel3,bdel4,bdel5
      DOUBLE PRECISION del1,del2,del3,del4,del5
      DOUBLE PRECISION st1neutt(5),st2neutt(5),st1charb(2),st2charb(2),
     .         st1hcsb(2),st2hcsb(2),st1wsb(2),st2wsb(2),
     .         qcdst1neut(5),qcdst2neut(5),qcdst1charb(2),
     .         qcdst2charb(2),qcdst1hcsb(2),qcdst2hcsb(2),
     .         qcdst1wsb(2),qcdst2wsb(2),
     .         st1glui,st2glui,qcdst1glui,qcdst2glui,
     .         st2H(3),qcdst2H(3),st2A(2),qcdst2A(2),
     .         st2ztop,qcdst2ztop,st1Tgra,st2Tgra
      DOUBLE PRECISION stoptot(2),stoptot2(2),stoptotmulti(2),
     .         stoptotrad(2),stoptot2lo(2),stoptot2nlo(2)
      DOUBLE PRECISION gamma,gammaup,gammagluino,gamma2
      DOUBLE PRECISION sigmato,sigmasn,sigmaw,sigmah,ainter,
     .         sigmae,sigmaq
      DOUBLE PRECISION xintegstopw(2,5),xintegstoph(2,5),
     .         xintegststau(2,2),xintegstsntau(2,2),xintegstsmu(2,2),
     .         xintegstsnmu(2),xintegstsel(2,2),xintegstsnel(2),
     .         xintegstbsbst(2,2),xintegstbbsbt(2,2),
     .         xintegsttausbnu(2,2),xintegstelsbnu(2,2),
     .         xintegstupsbdow(2,2),xintegst2st1tt,xintegst2st1startt,
     .         xintegst2st1bb,xintegst2st1uu,xintegst2st1dd,
     .         xintegst2st1ee,xintegst2st1nunu,xintegst2st1tautau
      DOUBLE PRECISION brst1neutt(5),brst2neutt(5),brst1charb(2),
     .         brst2charb(2),brst1hcsb(2),brst2hcsb(2),brst1wsb(2),
     .         brst2wsb(2),brst1glui,brst2glui,brst2H(3),brst2A(2),
     .         brst2ztop,brgamma,brgammaup,brgammagluino
      DOUBLE PRECISION brstopw(2,5),brstoph(2,5),brststau(2,2),
     .         brstsntau(2,2),brstsmu(2,2),brstsnmu(2),brstsel(2,2),
     .         brstsnel(2),brstbsbst(2,2),brstbbsbt(2,2),
     .         brsttausbnu(2,2),brstelsbnu(2,2),brstupsbdow(2,2),
     .         brst2st1tt,brst2st1startt,brst2st1bb,brst2st1uu,
     .         brst2st1dd,brst2st1ee,brst2st1nunu,brst2st1tautau
      DOUBLE PRECISION NS_lamb,resum
      DOUBLE PRECISION
     .         NS_gltneut,NS_grtneut,NS_gltchar,NS_grtchar,NS_corrreali
      DOUBLE PRECISION NS_gamtop1,NS_gamtop2,NS_gamglui1,NS_gamglui2,
     .         NS_gamglui3,NS_gam11,NS_gam12,NS_gamvirt,NS_gamreal,
     .         NS_gamcfdec
      DOUBLE PRECISION NS_gvirtgl,NS_gvirtmix,NS_stopsbot1719,
     .         NS_dcounterhc,NS_realcorr
      DOUBLE PRECISION NS_topneut1719,NS_dcounterneut,NS_gvirtmixdiv
      DOUBLE PRECISION NS_gluonvertex,NS_gluinoWvertex,NS_gluinoZvertex
      DOUBLE PRECISION NS_wavefuncvertex,NS_quarkmixW,NS_quarkmixZ
      DOUBLE PRECISION NS_realgluonem
      DOUBLE PRECISION NS_hisakakob1,NS_hisakakob2
      DOUBLE PRECISION flagmulti,flagqcd,flagloop,multilim
      DOUBLE PRECISION brcharWgra(2),brcharHCgra(2),brneutGAMgra(5),
     .         brneutZgra(5),brneutHgra(5,3),brneutAgra(5,2),
     .         brgluiGLUgra,brselEgra,brserEgra,brsmu1MUgra,
     .         brsmu2MUgra,brstau1TAUgra,brstau2TAUgra,brsneNEgra,
     .         brsnmNMgra,brsntNTgra,brsulUgra,brsurUgra,brsdlDgra,
     .         brsdrDgra,brst1Tgra,brst2Tgra,brsb1Bgra,brsb2Bgra
      DOUBLE PRECISION KNG(5),KNZ(5),KNH(5,3),KNA(5,2),KCW(2),KCH(2)
      DOUBLE PRECISION M32,CGR,MPL
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/STOP_WIDTH/stoptot,stoptot2,stoptotmulti,stoptotrad
      COMMON/STOP_BR_2BD/brst1neutt,brst2neutt,brst1charb,
     .         brst2charb,brst1hcsb,brst2hcsb,brst1wsb,
     .         brst2wsb,brst1glui,brst2glui,brst2H,brst2A,
     .         brst2ztop,brgamma,brgammaup,brgammagluino
      COMMON/STOP_BR_3BD/brstopw,brstoph,brststau,
     .         brstsntau,brstsmu,brstsnmu,brstsel,
     .         brstsnel,brstbsbst,brstbbsbt,
     .         brsttausbnu,brstelsbnu,brstupsbdow,
     .         brst2st1tt,brst2st1startt,brst2st1bb,brst2st1uu,
     .         brst2st1dd,brst2st1ee,brst2st1nunu,brst2st1tautau
      COMMON/GRAVITINO/brcharWgra,brcharHCgra,brneutGAMgra,
     .         brneutZgra,brneutHgra,brneutAgra,
     .         brgluiGLUgra,brselEgra,brserEgra,brsmu1MUgra,
     .         brsmu2MUgra,brstau1TAUgra,brstau2TAUgra,brsneNEgra,
     .         brsnmNMgra,brsntNTgra,brsulUgra,brsurUgra,brsdlDgra,
     .         brsdrDgra,brst1Tgra,brst2Tgra,brsb1Bgra,brsb2Bgra
      COMMON/GRAVCOUP/KNG,KNZ,KNH,KNA,KCW,KCH
      COMMON/M32/M32,CGR,MPL,GRFLAG
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_sfmixang/thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .         cmm,smm,cum,sum,cdm,sdm,cem,sem,cnm,snm
      COMMON/NS_runmcalc/rmtc,rmbc,rmtauc
      COMMON/NS_refscale/amurefer
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_coup19/gztt,gzbb,gztautau,gzmumu
      COMMON/NS_coup20/gwtb,gwntau,gwnmu
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_HIGGSSTST/Hstopstopr,Astopstopr
      COMMON/NS_HIGGSBTBT/Hsbotsbotr,Asbotsbotr
      COMMON/NS_HIGGSTAUTAU/Hstaustaur,Astaustaur
      COMMON/NS_FLAGS/flagmulti,flagqcd,flagloop
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_charstopbot/alstor,akstor
      COMMON/NS_neutstoptop/atopr,btopr
      COMMON/NS_hcsbotstop/gctbr
      COMMON/NS_multilim/multilim

      EXTERNAL NS_lamb
      EXTERNAL NS_gltneut,NS_grtneut,NS_corrreali,NS_gltchar,NS_grtchar
      EXTERNAL NS_gamtop1,NS_gamtop2,NS_gamglui1,NS_gamglui2,
     .         NS_gamglui3,NS_gam11,NS_gam12,NS_gamvirt,NS_gamreal,
     .         NS_gamcfdec
      EXTERNAL NS_gvirtgl,NS_gvirtmix,NS_stopsbot1719,NS_dcounterhc,
     .         NS_realcorr
      EXTERNAL NS_topneut1719,NS_dcounterneut,NS_gvirtmixdiv
      EXTERNAL NS_gluonvertex,NS_gluinoWvertex,NS_gluinoZvertex
      EXTERNAL NS_wavefuncvertex,NS_quarkmixW,NS_quarkmixZ
      EXTERNAL NS_realgluonem
      EXTERNAL NS_hisakakob1,NS_hisakakob2
      EXTERNAL resum

*  Initialization

      do i=1,2
         stoptot(i)    = 0d0
         stoptot2(i)   = 0d0
         stoptot2lo(i) = 0d0
         stoptot2nlo(i)= 0d0
         stoptotmulti(i) = 0d0
         stoptotrad(i) = 0d0

         st1charb(i) = 0d0
         st1hcsb(i)  = 0d0
         st1wsb(i)   = 0d0
         st2charb(i) = 0d0
         st2hcsb(i)  = 0d0
         st2wsb(i)   = 0d0
         qcdst1charb(i) = 0d0
         qcdst1hcsb(i)  = 0d0
         qcdst1wsb(i)   = 0d0
         qcdst2charb(i) = 0d0
         qcdst2hcsb(i)  = 0d0
         qcdst2wsb(i)   = 0d0
         xintegstsnel(i) = 0d0
         xintegstsnmu(i) = 0d0

         do j=1,5
            xintegstopw(i,j) = 0d0
            xintegstoph(i,j) = 0d0
         enddo

         do j=1,2
            xintegststau(i,j)    = 0d0
            xintegstsntau(i,j)   = 0d0
            xintegstsmu(i,j)     = 0d0
            xintegstsel(i,j)     = 0d0
            xintegstbsbst(i,j)   = 0d0
            xintegstbbsbt(i,j)   = 0d0
            xintegsttausbnu(i,j) = 0d0
            xintegstelsbnu(i,j)  = 0d0
            xintegstupsbdow(i,j) = 0d0
         enddo
      enddo
*
      do i=1,5
         st1neutt(i) = 0d0
         st2neutt(i) = 0d0
         qcdst1neut(i) = 0d0
         qcdst2neut(i) = 0d0
      enddo
      do i=1,3
         st2H(i)=0d0
         qcdst2H(i)=0d0
      enddo
      do i=1,2
         st2A(i)=0d0
         qcdst2A(i)=0d0
      enddo
      st1glui = 0d0
      st2glui = 0d0
      st2ztop = 0d0
      qcdst1glui = 0d0
      qcdst2glui = 0d0
      qcdst2ztop = 0d0
      xintegst2st1tt     = 0d0
      xintegst2st1bb     = 0d0
      xintegst2st1uu     = 0d0
      xintegst2st1dd     = 0d0
      xintegst2st1ee     = 0d0
      xintegst2st1nunu   = 0d0
      xintegst2st1tautau = 0d0
      xintegst2st1startt = 0d0
      gamma       = 0d0
      gammaup     = 0d0
      gammagluino = 0d0
      gamma2=0d0
*
      sigmato= 0d0
      sigmasn= 0d0
      sigmaw= 0d0
      sigmah= 0d0
      ainter= 0d0
      sigmae= 0d0
      sigmaq= 0d0
*
      do i=1,5
         brst1neutt(i)=0d0
         brst2neutt(i)=0d0
      enddo
      brst1glui=0d0
      brst2glui=0d0
      brst2ztop=0d0
      do i=1,2
         brst1charb(i)=0d0
         brst2charb(i)=0d0
         brst1hcsb(i)=0d0
         brst2hcsb(i)=0d0
         brst1wsb(i)=0d0
         brst2wsb(i)=0d0
         do j=1,5
            brstopw(i,j)=0d0
            brstoph(i,j)=0d0
         enddo
         do j=1,2
            brststau(i,j)=0d0
            brstsntau(i,j)=0d0
            brstsmu(i,j)=0d0
            brstsel(i,j)=0d0
            brstbsbst(i,j)=0d0
            brstbbsbt(i,j)=0d0
            brsttausbnu(i,j)=0d0
            brstelsbnu(i,j)=0d0
            brstupsbdow(i,j)=0d0
         enddo
         brstsnel(i)=0d0
         brstsnmu(i)=0d0
      enddo
      do i=1,3
         brst2H(i)=0d0
      enddo
      do i=1,2
         brst2A(i)=0d0
      enddo
      brst2st1tt=0d0
      brst2st1startt=0d0
      brst2st1bb=0d0
      brst2st1uu=0d0
      brst2st1dd=0d0
      brst2st1ee=0d0
      brst2st1nunu=0d0
      brst2st1tautau=0d0
      brgamma=0d0
      brgammaup=0d0
      brgammagluino=0d0
*
      do i=1,2
         gmsb(i)=0d0
      enddo
      amsq=0d0
      do i=1,3
         delta1H(i)=0d0
         delta2H(i)=0d0
         delta3H(i)=0d0
         delta4H(i)=0d0
         delta5H(i)=0d0
      enddo
      do i=1,2
         delta1A(i)=0d0
         delta2A(i)=0d0
         delta3A(i)=0d0
         delta4A(i)=0d0
         delta5A(i)=0d0
      enddo
      delta11c=0d0
      delta12c=0d0
      delta13c=0d0
      delta14c=0d0
      delta15c=0d0
      delta21c=0d0
      delta22c=0d0
      delta23c=0d0
      delta24c=0d0
      delta25c=0d0
      adel1=0d0
      adel2=0d0
      adel3=0d0
      adel4=0d0
      adel5=0d0
      bdel1=0d0
      bdel2=0d0
      bdel3=0d0
      bdel4=0d0
      bdel5=0d0

      st1Tgra=0d0
      st2Tgra=0d0
c ------------------------------------------------------------------ c
      gmsb(1) = asb1
      gmsb(2) = asb2
c -------------------------------------------------------------------- c
c For QCD corrections: the fixed scale is amurefer = Q, where
c the couplings are defined:
      amuvdiv = amuv/amurefer
c -------------------------------------------------------------------- c
c  stop1 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + top
      do i=1,5
         if(ast1.gt.(amneut(i)+amt)) then
            st1neutt(i)=g2s*((atopr(1,i)**2+btopr(1,i)**2)*(ast1**2-
     .           amt**2-amneut(i)**2)-4d0*atopr(1,i)*btopr(1,i)*
     .           amt*xmneut(i))*NS_lamb(amt/ast1,amneut(i)/ast1)
     .           /(16*pi*ast1)

         else
            st1neutt(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do j=1,5
         if(ast1.gt.(amneut(j)+amt)) then
            if(xmneut(j).le.0d0) then
               jsign = 1
            else
               jsign = 0
            endif

            qcdst1neut(j) = -g2s/24d0/pi**2/ast1*gs2/(4d0*pi)*
     .           ((btopr(1,j)*NS_gltneut(1,j,amuv,amuvdiv,lamv)+
     .             atopr(1,j)*NS_grtneut(1,j,amuv,amuvdiv,lamv))*
     .            (ast1**2-amt**2-amneut(j)**2)
     .           -2d0*(btopr(1,j)*NS_grtneut(1,j,amuv,amuvdiv,lamv)
     .                 +atopr(1,j)*NS_gltneut(1,j,amuv,amuvdiv,lamv))*
     .           amt*xmneut(j))*NS_lamb(amt/ast1,amneut(j)/ast1)
     .           +g2s/(6d0*pi**2*ast1)*gs2/(4d0*pi)*
     .           NS_corrreali(amt,amneut(j),ast1,lamv,1,jsign,1,j,1)

         else
            qcdst1neut(j) = 0d0
         endif
      enddo
      endif

c  ---stop2 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + top
      do i=1,5
         if(ast2.gt.(amneut(i)+amt)) then
            st2neutt(i)=g2s*((atopr(2,i)**2+btopr(2,i)**2)*(ast2**2-
     .           amt**2-amneut(i)**2)-4d0*atopr(2,i)*btopr(2,i)*
     .           amt*xmneut(i))*NS_lamb(amt/ast2,amneut(i)/ast2)
     .           /(16*pi*ast2)
         else
            st2neutt(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do j=1,5
         if(ast2.gt.(amneut(j)+amt)) then
            if(xmneut(j).le.0d0) then
               jsign = 1
            else
               jsign = 0
            endif

            qcdst2neut(j) = -g2s/24d0/pi**2/ast2*gs2/(4d0*pi)*
     .           ((btopr(2,j)*NS_gltneut(2,j,amuv,amuvdiv,lamv)
     .            +atopr(2,j)*NS_grtneut(2,j,amuv,amuvdiv,lamv))*
     .           (ast2**2-amt**2-amneut(j)**2)
     .           -2d0*(btopr(2,j)*NS_grtneut(2,j,amuv,amuvdiv,lamv)
     .                 +atopr(2,j)*NS_gltneut(2,j,amuv,amuvdiv,lamv))*
     .           amt*xmneut(j))*NS_lamb(amt/ast2,amneut(j)/ast2)
     .           +g2s/(6d0*pi**2*ast2)*gs2/(4d0*pi)*
     .           NS_corrreali(amt,amneut(j),ast2,lamv,1,jsign,2,j,1)
         else
            qcdst2neut(j) = 0d0
         endif
      enddo
      endif

c  stop1 --> chi+_1/chi+_2 + bottom
       do i=1,2
         if(ast1.gt.(amchar(i)+amb)) then
            st1charb(i)=g2s*((alstor(1,i)**2+akstor(1,i)**2)*
     .           (ast1**2-amb**2-amchar(i)**2)
     .           -4d0*alstor(1,i)*akstor(1,i)*
     .           amb*xmchar(i))*NS_lamb(amb/ast1,amchar(i)/ast1)
     .           /(16*pi*ast1)
         else
            st1charb(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do j=1,2
         if(ast1.gt.(amchar(j)+amb)) then
            jsign = 0

            qcdst1charb(j) = -g2s/24d0/pi**2/ast1*gs2/(4d0*pi)*
     .           ((akstor(1,j)*NS_gltchar(1,j,amuv,amuvdiv,lamv)
     .            +alstor(1,j)*NS_grtchar(1,j,amuv,amuvdiv,lamv))*
     .           (ast1**2-amb**2-amchar(j)**2)
     .           -2d0*(akstor(1,j)*NS_grtchar(1,j,amuv,amuvdiv,lamv)
     .                 +alstor(1,j)*NS_gltchar(1,j,amuv,amuvdiv,lamv))*
     .           amb*xmchar(j))*NS_lamb(amb/ast1,amchar(j)/ast1)
     .           +g2s/(6d0*pi**2*ast1)*gs2/(4d0*pi)*
     .           NS_corrreali(amb,amchar(j),ast1,lamv,2,jsign,1,j,1)

         else
            qcdst1charb(j) = 0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  stop2 --> chi+_1/chi+_2 + bottom
      do i=1,2
         if(ast2.gt.(amchar(i)+amb)) then
            st2charb(i)=g2s*((alstor(2,i)**2+akstor(2,i)**2)*
     .           (ast2**2-amb**2-amchar(i)**2)
     .           -4d0*alstor(2,i)*akstor(2,i)*
     .           amb*xmchar(i))*NS_lamb(amb/ast2,amchar(i)/ast2)
     .           /(16*pi*ast2)
         else
            st2charb(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do j=1,2
         if(ast2.gt.(amchar(j)+amb)) then
            jsign = 0

            qcdst2charb(j) = -g2s/24d0/pi**2/ast2*gs2/(4d0*pi)*
     .           ((akstor(2,j)*NS_gltchar(2,j,amuv,amuvdiv,lamv)
     .            +alstor(2,j)*NS_grtchar(2,j,amuv,amuvdiv,lamv))*
     .           (ast2**2-amb**2-amchar(j)**2)
     .           -2d0*(akstor(2,j)* NS_grtchar(2,j,amuv,amuvdiv,lamv)
     .                 +alstor(2,j)*NS_gltchar(2,j,amuv,amuvdiv,lamv))*
     .           amb*xmchar(j))*NS_lamb(amb/ast2,amchar(j)/ast2)
     .           +g2s/(6d0*pi**2*ast2)*gs2/(4d0*pi)*
     .      NS_corrreali(amb,amchar(j),ast2,lamv,2,jsign,2,j,1)
         else
            qcdst2charb(j) = 0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  stop1 --> gluino + top
      if(ast1.gt.(mgluino+amt)) then
         st1glui = 8d0*gs2*(ast1**2-amt**2-mgluino**2+4d0*amt*
     .        mgluino*DDSIN(thet)*DDCOS(thet))*
     .        NS_lamb(amt/ast1,mgluino/ast1)/(16d0*pi*ast1)/3d0
      else
         st1glui = 0d0
      endif
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      if(ast1.gt.(mgluino+amt)) then
         amsq    = 2d0*(asup1+asup2+asdown1+asdown2)/8d0
         scalmur = amurefer
         alp     = gs2/(4d0*pi)
         nf      = 6d0
         qcdst1glui = 8d0*pi*alp/3d0/amt**2*st1glui*
     .        NS_gamtop1(ast1,ast2,amt,mgluino,thet,1,amuv,lamv) +
     .        16d0*pi*alp**2*NS_lamb(amt/ast1,mgluino/ast1)/
     .        (9d0*amt**2*ast1)*
     .        NS_gamtop2(ast1,ast2,amt,mgluino,thet,1,amuv) +
     .        4d0*pi*alp/mgluino**2*st1glui*(nf-2d0)*
     .        NS_gamglui1(ast1,ast2,amsq,amt,mgluino,amuv) +
     .        2d0*pi*alp/mgluino**2*st1glui*
     .        NS_gamglui2(ast1,ast2,amt,thet,asb1,asb2,amb,theb,
     .                    mgluino,1,amuv) +
     .        4d0*pi*alp*3d0/mgluino**2*st1glui*
     .        NS_gamglui3(mgluino,amuv,lamv) +
     .        8d0*4d0/3d0*pi*alp*st1glui*
     .        NS_gam11(ast1,ast2,amt,mgluino,thet,1,amuv,lamv) +
     .        8d0*16d0/9d0*pi*alp**2/ast1*
     .        NS_lamb(amt/ast1,mgluino/ast1)*
     .        NS_gam12(ast1,ast2,amt,mgluino,thet,1,amuv,lamv,scalmur) +
     .        alp**2*NS_lamb(amt/ast1,mgluino/ast1)/ast1*
     .        NS_gamvirt(ast1,ast2,amt,mgluino,thet,1,amuv,lamv) +
     .        alp**2*NS_gamreal(ast1,amt,mgluino,thet,1,lamv) +
     .        alp/(4d0*pi)*st1glui*
     .        NS_gamcfdec(ast1,ast2,amt,asb1,asb2,amb,mgluino,amsq,amuv,
     .        scalmur)
      else
         qcdst1glui = 0d0
      endif
      endif
c -------------------------------------------------------------------- c
c  stop2 --> gluino + top

      if(ast2.gt.(mgluino+amt)) then
         st2glui = 8d0*gs2*((ast2**2-amt**2-mgluino**2)-4d0*amt*
     .        mgluino*DDSIN(thet)*DDCOS(thet))*
     .        NS_lamb(amt/ast2,mgluino/ast2)/(16d0*pi*ast2)/3d0
      else
         st2glui = 0d0
      endif
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      if(ast2.gt.(mgluino+amt)) then
         amsq    = 2d0*(asup1+asup2+asdown1+asdown2)/8d0
         scalmur = amurefer
         alp     = gs2/(4d0*pi)
         nf      = 6d0
         qcdst2glui = 8d0*pi*alp/3d0/amt**2*st2glui*
     .        NS_gamtop1(ast2,ast1,amt,mgluino,thet,2,amuv,lamv) +
     .        16d0*pi*alp**2*NS_lamb(amt/ast2,mgluino/ast2)/
     .        (9d0*amt**2*ast2)*
     .        NS_gamtop2(ast2,ast1,amt,mgluino,thet,2,amuv) +
     .        4d0*pi*alp/mgluino**2*st2glui*(nf-2d0)*
     .        NS_gamglui1(ast2,ast1,amsq,amt,mgluino,amuv) +
     .        2d0*pi*alp/mgluino**2*st2glui*
     .        NS_gamglui2(ast2,ast1,amt,thet,asb2,asb1,amb,theb,
     .                    mgluino,2,amuv) +
     .        4d0*pi*alp*3d0/mgluino**2*st2glui*
     .        NS_gamglui3(mgluino,amuv,lamv) +
     .        8d0*4d0/3d0*pi*alp*st2glui*
     .        NS_gam11(ast2,ast1,amt,mgluino,thet,2,amuv,lamv) +
     .        8d0*16d0/9d0*pi*alp**2/ast2*
     .        NS_lamb(amt/ast2,mgluino/ast2)*
     .        NS_gam12(ast2,ast1,amt,mgluino,thet,2,amuv,lamv,scalmur) +
     .        alp**2*NS_lamb(amt/ast2,mgluino/ast2)/ast2*
     .        NS_gamvirt(ast2,ast1,amt,mgluino,thet,2,amuv,lamv) +
     .        alp**2*NS_gamreal(ast2,amt,mgluino,thet,2,lamv) +
     .        alp/(4d0*pi)*st2glui*
     .        NS_gamcfdec(ast2,ast1,amt,asb2,asb1,amb,mgluino,amsq,amuv,
     .        scalmur)
      else
         qcdst2glui = 0d0
      endif
      endif
c -------------------------------------------------------------------- c
c  stop1 --> H+ + sbottom1/2
      do i=1,2
         if(ast1.gt.(gmsb(i)+cmass)) then
            st1hcsb(i)=g2s*amw**2*gctbr(1,i)**2*
     .           NS_lamb(gmsb(i)/ast1,cmass/ast1)/(16d0*pi*ast1)
         else
            st1hcsb(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do nj=1,2
      alp = gs2/(4d0*pi)
      if(ast1.gt.(gmsb(nj)+cmass)) then
         delta11c = -dsqrt(2d0)*amw**2*gctbr(1,nj)*
     .        NS_gvirtgl(ast1,cmass,gmsb(nj),lamv,amuv)
         delta12c =  -dsqrt(2d0)*amw**2*gctbr(1,3-nj)*
     .        1d0/(gmsb(3-nj)**2-gmsb(nj)**2)*
     .        NS_gvirtmix(asb1,asb2,gmsb(nj),mgluino,rmbc,theb,amuv)
     .        -dsqrt(2d0)*amw**2*gctbr(2,nj)*
     .        1d0/(ast2**2-ast1**2)*
     .        NS_gvirtmix(ast1,ast2,ast1,mgluino,rmtc,thet,amuv)
         delta13c = NS_stopsbot1719(amuv,1,nj)
         delta14c = NS_dcounterhc(ast1,rmtc,thet,1,gmsb(nj),rmbc,
     .                         theb,nj,mgluino,amuv,amuvdiv,lamv,1,nj)
         delta15c = NS_realcorr(cmass,ast1,gmsb(nj),lamv,6,0,1,nj,ast1)
         qcdst1hcsb(nj) = -g2s*amw**2/(24d0*dsqrt(2d0)*pi*amw**2*
     .        ast1)*alp/pi*NS_lamb(gmsb(nj)/ast1,cmass/ast1)*
     .        gctbr(1,nj)*(delta11c+delta12c+delta13c+delta14c+delta15c)
      else
         qcdst1hcsb(nj) = 0d0
      endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  stop2 --> H(K) + stop1
      do K=1,3
         if(ast2.gt.(ast1+SMASS(K))) then
            st2H(K)=g2s*amz**4/amw**2*Hstopstopr(K,2,1)**2*
     .           NS_lamb(ast1/ast2,SMASS(K)/ast2)/(16d0*pi*ast2)
         else
            st2H(K)=0d0
         endif
      enddo

c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do K=1,3
         if(ast2.gt.(ast1+SMASS(K))) then
            alp = gs2/(4d0*pi)
            delta1H(K) = -dsqrt(2d0)*amz**2*Hstopstopr(K,2,1)*
     .           NS_gvirtgl(ast2,MAX(1d0,SMASS(K)),ast1,lamv,amuv)
            delta2H(K) =  -dsqrt(2d0)*amz**2*Hstopstopr(K,2,2)*
     .        1d0/(ast2**2-ast1**2)*
     .        NS_gvirtmix(ast1,ast2,ast1,mgluino,rmtc,thet,amuv)
     .        -dsqrt(2d0)*amz**2*Hstopstopr(K,1,1)*
     .        1d0/(ast1**2-ast2**2)*
     .        NS_gvirtmix(ast1,ast2,ast2,mgluino,rmtc,thet,amuv)
             delta3H(K) = NS_topneut1719(K,amuv)
             delta4H(K) = NS_dcounterneut(ast1,ast2,rmtc,thet,mgluino,
     .           amuv,amuvdiv,lamv,1,K)
             delta5H(K) = NS_realcorr(MAX(1d0,SMASS(K)),ast2,ast1,
     .           lamv,K,1,2,1,ast2)

            qcdst2H(K) = -g2s*amz**2/(24d0*dsqrt(2d0)*pi*amw**2*ast2)*
     .           alp/pi*NS_lamb(ast1/ast2,SMASS(K)/ast2)*
     .           Hstopstopr(K,2,1)*
     .           (delta1H(K)+delta2H(K)+delta3H(K)+delta4H(K)+
     .           delta5H(K))
         else
            qcdst2H(K) = 0d0
         endif
      enddo
      endif

c -------------------------------------------------------------------- c
c  stop2 --> A(K) + stop1

      do K=1,2
         if(ast2.gt.(ast1+PMASS(K))) then
            st2A(K)=g2s*amz**4/amw**2*Astopstopr(K,2,1)**2*
     .           NS_lamb(ast1/ast2,PMASS(K)/ast2)/(16d0*pi*ast2)
         else
            st2A(K)=0d0
         endif
      enddo

c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do K=1,2
         if(ast2.gt.(ast1+PMASS(K))) then
            alp = gs2/(4d0*pi)
            delta1A(K) = dsqrt(2d0)*amz**2*Astopstopr(K,2,1)*
     .           NS_gvirtgl(ast2,MAX(1d0,PMASS(K)),ast1,lamv,amuv)
            delta2A(K) = 0d0
            delta3A(K) = NS_topneut1719(3+k,amuv)
            delta4A(K) = NS_dcounterneut(ast1,ast2,rmtc,thet,mgluino,
     .           amuv,amuvdiv,lamv,1,K+3)
            delta5A(K) = NS_realcorr(MAX(1d0,PMASS(K)),ast2,
     .           ast1,lamv,K+3,1,2,1,ast2)
            qcdst2A(K) = g2s*amz**2/(24d0*dsqrt(2d0)*pi*amw**2*ast2)*
     .           alp/pi*NS_lamb(ast1/ast2,PMASS(K)/ast2)*
     .           Astopstopr(K,2,1)*
     .           (delta1A(K)+delta2A(K)+delta3A(K)+delta4A(K)
     .           +delta5A(K))
         else
            qcdst2A(K) = 0d0
         endif
      enddo
      endif

c -------------------------------------------------------------------- c
c  stop2 --> H+ + sbottom1/2
      do i=1,2
         if(ast2.gt.(gmsb(i)+CMASS)) then
            st2hcsb(i)=g2s*amw**2*gctbr(2,i)**2*
     .           NS_lamb(gmsb(i)/ast2,CMASS/ast2)/(16d0*pi*ast2)
         else
            st2hcsb(i)=0d0
         endif
      enddo

c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do nj=1,2
      if(ast2.gt.(gmsb(nj)+cmass)) then
         alp = gs2/(4d0*pi)
         delta21c = -dsqrt(2d0)*amw**2*gctbr(2,nj)*
     .        NS_gvirtgl(ast2,cmass,gmsb(nj),lamv,amuv)
         delta22c =  -dsqrt(2d0)*amw**2*gctbr(2,3-nj)*
     .        1d0/(gmsb(3-nj)**2-gmsb(nj)**2)*
     .        NS_gvirtmix(asb1,asb2,gmsb(nj),mgluino,rmbc,theb,amuv)
     .        -dsqrt(2d0)*amw**2*gctbr(1,nj)*
     .        1d0/(ast1**2-ast2**2)*
     .        NS_gvirtmix(ast1,ast2,ast2,mgluino,rmtc,thet,amuv)
         delta23c = NS_stopsbot1719(amuv,2,nj)
         delta24c = NS_dcounterhc(ast2,rmtc,thet,2,gmsb(nj),rmbc,
     .                         theb,nj,mgluino,amuv,amuvdiv,lamv,2,nj)
         delta25c = NS_realcorr(cmass,ast2,gmsb(nj),lamv,6,0,2,nj,ast2)

         qcdst2hcsb(nj) = -g2s*amw**2/(24d0*dsqrt(2d0)*pi*amw**2*
     .        ast2)*alp/pi*NS_lamb(gmsb(nj)/ast2,cmass/ast2)*
     .        gctbr(2,nj)*(delta21c+delta22c+delta23c+delta24c+delta25c)

      else
         qcdst2hcsb(nj) = 0d0
      endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  stop2 --> Z + stop1
      if(ast2.gt.(ast1+mz)) then
         st2ztop=g2s/64d0/pi/cw**2/mz**2*ast2**3*gztt(2,1)**2*
     .           NS_lamb(ast1/ast2,mz/ast2)**3
      else
         st2ztop=0d0
      endif

c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      if(ast2.gt.(ast1+mz)) then
         alp = gs2/(4d0*pi)
         del1 = -alp/3d0/pi*gztt(2,1)/2d0/cw*
     .        NS_gluonvertex(ast2,ast1,mz,lamv,amuv)
         del2 = -alp/3d0/pi/cw*NS_gluinoZvertex(ast2,ast1,mz,lamv,
     .        amuv,mgluino,amt,1d0/2d0,2d0/3d0,sw,thet)
         del3 = alp/pi*NS_wavefuncvertex(ast2,ast1,amt,amt,
     .     thet,thet,1d0,1d0,2,1,mgluino,lamv,amuv)
         del4 = alp/pi*NS_quarkmixZ(ast2,thet,1d0/2d0,2d0/3d0,ast1,
     .     ast2,amt,mgluino,amuv)
         del5 = NS_realgluonem(ast2,ast1,mz,lamv)
         qcdst2ztop =  g2s/16d0/pi/mz**2*ast2**3*(gztt(2,1)/2d0/cw)*
     .        NS_lamb(ast1/ast2,mz/ast2)**3*(2d0*del1+2d0*del2
     .        +2d0*del3+2d0*del4) +
     .        g2s/3d0/pi**2/ast2*alp*(gztt(2,1)/(2d0*cw))**2*del5
      else
         qcdst2ztop = 0d0
      endif
      endif
c -------------------------------------------------------------------- c
c  stop1 --> W+ + sbottom1/2
      do i=1,2
         if(ast1.gt.(gmsb(i)+mw)) then
            st1wsb(i)=g2s/32d0/pi/mw**2*ast1**3*gwtb(1,i)**2*
     .                NS_lamb(gmsb(i)/ast1,mw/ast1)**3
         else
            st1wsb(i)=0d0
         endif
      enddo
c -- QCD corrections --
      if(flagqcd.eq.1d0) then
      do i=1,2
         if(ast1.gt.(gmsb(i)+mw)) then
            alp = gs2/(4d0*pi)
            adel1 = -alp/3d0/pi*gwtb(1,i)/dsqrt(2d0)*
     .           NS_gluonvertex(ast1,gmsb(i),mw,lamv,amuv)
            adel2 = -dsqrt(2d0)/3d0*alp/pi*
     .           NS_gluinoWvertex(ast1,gmsb(i),mw,lamv,amuv,mgluino,
     .           amt,amb,thet,theb,1,i)
            adel3 = alp/pi*NS_wavefuncvertex(ast1,gmsb(i),amt,amb,
     .           thet,theb,2d0,1d0,1,i,mgluino,lamv,amuv)
            adel4 = alp/pi*NS_quarkmixW(ast1,thet,1d0/2d0,2d0/3d0,
     .             ast1,ast2,amt,gmsb(i),theb,-1d0/2d0,-1d0/3d0,
     .             asb1,asb2,amb,1,i,mgluino,amuv)
            adel5 = NS_realgluonem(ast1,gmsb(i),mw,lamv)
            qcdst1wsb(i) = g2s/16d0/pi/mw**2*ast1**3*
     .           (gwtb(1,i)/dsqrt(2d0))*
     .           NS_lamb(gmsb(i)/ast1,mw/ast1)**3*(2d0*adel1
     .           +2d0*adel2+2d0*adel3+2d0*adel4) +
     .         g2s/3d0/pi**2/ast1*alp*(gwtb(1,i)/dsqrt(2d0))**2*adel5
         else
            qcdst1wsb(i) = 0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  stop2 --> W+ + sbottom1/2
      do i=1,2
         if(ast2.gt.(gmsb(i)+mw)) then
            st2wsb(i)=g2s/32d0/pi/mw**2*ast2**3*gwtb(2,i)**2*
     .                NS_lamb(gmsb(i)/ast2,mw/ast2)**3
         else
            st2wsb(i)=0d0
         endif
      enddo
c -- QCD corrections --
      if(flagqcd.eq.1d0) then
      do i=1,2
         if(ast2.gt.(gmsb(i)+mw)) then
            alp = gs2/(4d0*pi)
            bdel1 = -alp/3d0/pi*gwtb(2,i)/dsqrt(2d0)*
     .           NS_gluonvertex(ast2,gmsb(i),mw,lamv,amuv)
            bdel2 = -dsqrt(2d0)/3d0*alp/pi*
     .           NS_gluinoWvertex(ast2,gmsb(i),mw,lamv,amuv,mgluino,
     .           amt,amb,thet,theb,2,i)
            bdel3 = alp/pi*NS_wavefuncvertex(ast2,gmsb(i),amt,amb,
     .           thet,theb,2d0,1d0,2,i,mgluino,lamv,amuv)
            bdel4 = alp/pi*NS_quarkmixW(ast2,thet,1d0/2d0,2d0/3d0,
     .             ast1,ast2,amt,gmsb(i),theb,-1d0/2d0,-1d0/3d0,
     .             asb1,asb2,amb,2,i,mgluino,amuv)
            bdel5 = NS_realgluonem(ast2,gmsb(i),mw,lamv)
            qcdst2wsb(i) = g2s/16d0/pi/mw**2*ast2**3*
     .           (gwtb(2,i)/dsqrt(2d0))*
     .           NS_lamb(gmsb(i)/ast2,mw/ast2)**3*(2d0*bdel1
     .           +2d0*bdel2+2d0*bdel3+2d0*bdel4) +
     .          g2s/3d0/pi**2/ast2*alp*(gwtb(2,i)/dsqrt(2d0))**2*bdel5
         else
            qcdst2wsb(i) = 0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  stop1 --> top + gravitino

      IF(GRFLAG.EQ.1)then
        if ((amt+M32).le.ast1) then
          st1Tgra = ast1**5/(48d0*PI*MPL**2*M32**2)
     .              *(1d0-amt**2/ast1**2)**4
        else
          st1Tgra = 0d0
        endif
      ENDIF
c -------------------------------------------------------------------- c
c  stop2 --> top + gravitino

      IF(GRFLAG.EQ.1)then
        if ((amt+M32).le.ast2) then
          st2Tgra = ast2**5/(48d0*PI*MPL**2*M32**2)
     .              *(1d0-amt**2/ast2**2)**4
        else
          st2Tgra = 0d0
        endif
      ENDIF
c ------------------------------------------------------------------- c
c ---------------- 2-body decays and 2-body total widths ------------ c
c ------------------------------------------------------------------- c

      stoptot2lo(1)=st1neutt(1)+st1neutt(2)+st1neutt(3)+st1neutt(4)+
     .            st1neutt(5)+st1charb(1)+st1charb(2)+st1glui+
     .            st1hcsb(1)+st1hcsb(2)+st1wsb(1)+st1wsb(2)+st1Tgra

      stoptot2lo(2)=st2neutt(1)+st2neutt(2)+st2neutt(3)+st2neutt(4)+
     .            st2neutt(5)+st2charb(1)+st2charb(2)+st2glui+
     .            st2H(1)+st2H(2)+st2H(3)+st2A(2)+st2A(1)+st2hcsb(1)+
     .            st2hcsb(2)+st2wsb(1)+st2wsb(2)+st2ztop+st2Tgra

c UE: Resum if qcdcorr < -tree:
      do i=1,5
!      If(qcdst1neut(i).lt.-st1neutt(i))
!     .write(0,23)"Warning: large negative rad. corrs. to st1->t+chi0_",i
        qcdst1neut(i)=resum(st1neutt(i),qcdst1neut(i))
!      If(qcdst2neut(i).lt.-st2neutt(i))
!     .write(0,23)"Warning: large negative rad. corrs. to st2->t+chi0_",i
        qcdst2neut(i)=resum(st2neutt(i),qcdst2neut(i))
      enddo
      do i=1,2
!      If(qcdst1charb(i).lt.-st1charb(i))
!     .write(0,23)"Warning: large negative rad. corrs. to st1->b+char_",i
        qcdst1charb(i)=resum(st1charb(i),qcdst1charb(i))
!      If(qcdst2charb(i).lt.-st2charb(i))
!     .write(0,23)"Warning: large negative rad. corrs. to st2->b+char_",i
        qcdst2charb(i)=resum(st2charb(i),qcdst2charb(i))
!      If(qcdst1hcsb(i).lt.-st1hcsb(i))
!     .write(0,23)"Warning: large negative rad. corrs. to st1->hc+sb_",i
        qcdst1hcsb(i)=resum(st1hcsb(i),qcdst1hcsb(i))
!      If(qcdst2hcsb(i).lt.-st2hcsb(i))
!     .write(0,23)"Warning: large negative rad. corrs. to st2->hc+sb_",i
        qcdst2hcsb(i)=resum(st2hcsb(i),qcdst2hcsb(i))
!      If(qcdst1wsb(i).lt.-st1wsb(i))
!     .write(0,23)"Warning: large negative rad. corrs. to st1->W+sb_",i
        qcdst1wsb(i)=resum(st1wsb(i),qcdst1wsb(i))
!      If(qcdst2wsb(i).lt.-st2wsb(i))
!     .write(0,23)"Warning: large negative rad. corrs. to st2->W+sb_",i
        qcdst2wsb(i)=resum(st2wsb(i),qcdst2wsb(i))
!      If(qcdst2A(i).lt.-st2A(i))
!     .write(0,23)"Warning: large negative rad. corrs. to st2->st1+A_",i
        qcdst2A(i)=resum(st2A(i),qcdst2A(i))
      enddo
!      If(qcdst1glui.lt.-st1glui)
!     .write(0,*)"Warning: large negative rad. corrs. to st1->glui+t"
      qcdst1glui=resum(st1glui,qcdst1glui)
!      If(qcdst2glui.lt.-st2glui)
!     .write(0,*)"Warning: large negative rad. corrs. to st2->glui+t"
      qcdst2glui=resum(st2glui,qcdst2glui)
      do i=1,3
!      If(qcdst2H(i).lt.-st2H(i))
!     .write(0,23)"Warning: large negative rad. corrs. to st2->st1+H_",i
        qcdst2H(i)=resum(st2H(i),qcdst2H(i))
      enddo
!      If(qcdst2ztop.lt.-st2ztop)
!     .write(0,*)"Warning: large negative rad. corrs. to st2->st1+Z"
      qcdst2ztop=resum(st2ztop,qcdst2ztop)
c UE: End resummation
23     format(A,I1)
c
      stoptot2nlo(1) = stoptot2lo(1) + qcdst1neut(1)+qcdst1neut(2)+
     .            qcdst1neut(3)+qcdst1neut(4)+qcdst1neut(5)+
     .            qcdst1charb(1)+
     .            qcdst1charb(2)+qcdst1glui+qcdst1hcsb(1)+qcdst1hcsb(2)+
     .            qcdst1wsb(1)+qcdst1wsb(2)
      stoptot2nlo(2) = stoptot2lo(2) + qcdst2neut(1)+qcdst2neut(2)+
     .            qcdst2neut(3)+qcdst2neut(4)+qcdst2neut(5)+
     .            qcdst2charb(1)+
     .            qcdst2charb(2)+qcdst2glui+qcdst2hcsb(1)+qcdst2hcsb(2)+
     .            qcdst2H(1)+qcdst2H(2)+qcdst2H(3)+
     .            qcdst2A(2)+qcdst2A(1)+
     .            qcdst2ztop+qcdst2wsb(1)+
     .            qcdst2wsb(2)

      if(flagqcd.eq.0d0) then
         do i=1,2
            stoptot2(i) = stoptot2lo(i)
         enddo
      elseif(flagqcd.eq.1d0) then
         do i=1,2
            stoptot2(i) = stoptot2nlo(i)
         enddo
      endif
c -------------------------------------------------------------------- c
c ---------------- 3-body decays and 3-body total widths ------------- c
c -------------------------------------------------------------------- c
      if(flagmulti.eq.1d0) then

         CALL NS_xintegstop(xintegstopw,xintegstoph,xintegststau,
     .     xintegstsntau,xintegstsmu,xintegstsnmu,xintegstsel,
     .     xintegstsnel,xintegstbsbst,xintegstbbsbt,xintegsttausbnu,
     .     xintegstelsbnu,xintegstupsbdow,xintegst2st1tt,xintegst2st1bb,
     .     xintegst2st1uu,xintegst2st1dd,xintegst2st1ee,
     .     xintegst2st1nunu,xintegst2st1tautau,xintegst2st1startt)

      do i=1,2
            stoptotmulti(i) = 0d0
            do j=1,5
               stoptotmulti(i) = stoptotmulti(i)+xintegstopw(i,j)+
     .              xintegstoph(i,j)
            enddo
            do j=1,2
               stoptotmulti(i) = stoptotmulti(i)+xintegststau(i,j)+
     .              xintegstsntau(i,j)+xintegstsel(i,j)+
     .              xintegstsmu(i,j)+xintegstbsbst(i,j)+
     .              xintegstbbsbt(i,j)+xintegsttausbnu(i,j)+
     .              2d0*xintegstelsbnu(i,j)+2d0*xintegstupsbdow(i,j)
            enddo
            stoptotmulti(i) = stoptotmulti(i)+xintegstsnel(i)+
     .              xintegstsnmu(i)
            if(i.eq.2) then
               stoptotmulti(i) = stoptotmulti(i)+
     .              xintegst2st1tt+xintegst2st1bb+
     .              2d0*xintegst2st1uu+2d0*xintegst2st1dd+
     .              2d0*xintegst2st1ee+3d0*xintegst2st1nunu+
     .              xintegst2st1tautau+xintegst2st1startt
            endif

      enddo
      endif
c Check if numerically relevant:
      do i=1,2
      if (stoptotmulti(i).lt.multilim*stoptot2(i))then
         stoptotmulti(i)=0d0
      endif
      enddo

c -------------------------------------------------------------------- c
c --------------- loop decays and the total widths ------------------- c
c -------------------------------------------------------------------- c
c Only if no treelevel decays:
      if(stoptot2(1).eq.0d0.and.flagloop.eq.1d0) then
         CALL NS_hikasakob1(gamma,gammaup,gammagluino)
         stoptotrad(1) = gamma+gammaup+gammagluino
      endif

c      if(stoptot2(2).eq.0d0.and.flagloop.eq.1d0) then
c         CALL NS_hikasakob2(gamma2)
c      endif

      do i =1,2
         stoptot(i)=stoptot2(i)+stoptotmulti(i)+stoptotrad(i)
      enddo

c ---------------------- stop branching ratios ------------------- c
      if(flagqcd.eq.1d0) then
         do i=1,5
            st1neutt(i) = st1neutt(i)+qcdst1neut(i)
            st2neutt(i) = st2neutt(i)+qcdst2neut(i)
         enddo
         do i=1,2
            st1charb(i) = st1charb(i) + qcdst1charb(i)
            st2charb(i) = st2charb(i) + qcdst2charb(i)
            st1hcsb(i)  = st1hcsb(i) + qcdst1hcsb(i)
            st2hcsb(i)  = st2hcsb(i) + qcdst2hcsb(i)
            st1wsb(i)   = st1wsb(i) + qcdst1wsb(i)
            st2wsb(i)   = st2wsb(i) + qcdst2wsb(i)
         enddo
         st1glui = st1glui + qcdst1glui
         st2glui = st2glui + qcdst2glui
         do K=1,3
            st2H(K)=st2H(K)+qcdst2H(K)
         enddo
         do K=1,2
            st2A(K)=st2A(K)+qcdst2A(K)
         enddo
        st2ztop = st2ztop + qcdst2ztop
      endif

      if(stoptot(1).ne.0d0)then

       do i=1,5
         brst1neutt(i) = st1neutt(i)/stoptot(1)
       enddo
       do i=1,2
         brst1charb(i) = st1charb(i)/stoptot(1)
         brst1hcsb(i)  = st1hcsb(i)/stoptot(1)
         brst1wsb(i)   = st1wsb(i)/stoptot(1)
       enddo
       brst1glui = st1glui/stoptot(1)
       brst1Tgra = st1Tgra/stoptot(1)

      else

       do i=1,5
         brst1neutt(i) = 0d0
       enddo
       do i=1,2
         brst1charb(i) = 0d0
         brst1hcsb(i)  = 0d0
         brst1wsb(i)   = 0d0
       enddo
       brst1glui = 0d0
       brst1Tgra = 0d0

      endif

      if(stoptot(2).ne.0d0)then

       do i=1,5
         brst2neutt(i) = st2neutt(i)/stoptot(2)
       enddo
       do i=1,2
         brst2charb(i) = st2charb(i)/stoptot(2)
         brst2hcsb(i)  = st2hcsb(i)/stoptot(2)
         brst2wsb(i)   = st2wsb(i)/stoptot(2)
       enddo
       brst2glui = st2glui/stoptot(2)
       do K=1,3
         brst2H(K) = st2H(K)/stoptot(2)
       enddo
       do K=1,2
         brst2A(K) = st2A(K)/stoptot(2)
       enddo
       brst2ztop = st2ztop/stoptot(2)
       brst2Tgra = st2Tgra/stoptot(2)

      else

       do i=1,5
         brst2neutt(i) = 0d0
       enddo
       do i=1,2
         brst2charb(i) = 0d0
         brst2hcsb(i)  = 0d0
         brst2wsb(i)   = 0d0
       enddo
       brst2glui = 0d0
       do K=1,3
         brst2H(K) = 0d0
       enddo
       do K=1,2
         brst2A(K) = 0d0
       enddo
       brst2ztop = 0d0
       brst2Tgra = 0d0

      endif

c ------------------------ 3-body and loop BRs ---------------------- c
      do i=1,2
         if (stoptotmulti(i).ne.0d0)then

            do j=1,5
               brstopw(i,j) = xintegstopw(i,j)/stoptot(i)
               brstoph(i,j) = xintegstoph(i,j)/stoptot(i)
            enddo
            do j=1,2
               brststau(i,j)    = xintegststau(i,j)/stoptot(i)
               brstsntau(i,j)   = xintegstsntau(i,j)/stoptot(i)
               brstsmu(i,j)     = xintegstsmu(i,j)/stoptot(i)
               brstsel(i,j)     = xintegstsel(i,j)/stoptot(i)
               brstbsbst(i,j)   = xintegstbsbst(i,j)/stoptot(i)
               brstbbsbt(i,j)   = xintegstbbsbt(i,j)/stoptot(i)
               brsttausbnu(i,j) = xintegsttausbnu(i,j)/stoptot(i)
               brstelsbnu(i,j)  = xintegstelsbnu(i,j)/stoptot(i)
               brstupsbdow(i,j) = xintegstupsbdow(i,j)/stoptot(i)
            enddo
            brstsnel(i)    = xintegstsnel(i)/stoptot(i)
            brstsnmu(i)    = xintegstsnmu(i)/stoptot(i)

	 else

            do j=1,5
               brstopw(i,j) = 0d0
               brstoph(i,j) = 0d0
            enddo
            do j=1,2
               brststau(i,j)    = 0d0
               brstsntau(i,j)   = 0d0
               brstsmu(i,j)     = 0d0
               brstsel(i,j)     = 0d0
               brstbsbst(i,j)   = 0d0
               brstbbsbt(i,j)   = 0d0
               brsttausbnu(i,j) = 0d0
               brstelsbnu(i,j)  = 0d0
               brstupsbdow(i,j) = 0d0
            enddo
            brstsnel(i)    = 0d0
            brstsnmu(i)    = 0d0

         endif
      enddo

      if (stoptotmulti(2).ne.0d0)then

         brst2st1tt     = xintegst2st1tt/stoptot(2)
         brst2st1startt = xintegst2st1startt/stoptot(2)
         brst2st1bb     = xintegst2st1bb/stoptot(2)
         brst2st1uu     = xintegst2st1uu/stoptot(2)
         brst2st1dd     = xintegst2st1dd/stoptot(2)
         brst2st1ee     = xintegst2st1ee/stoptot(2)
         brst2st1nunu   = xintegst2st1nunu/stoptot(2)
         brst2st1tautau = xintegst2st1tautau/stoptot(2)

      else

         brst2st1tt     = 0d0
         brst2st1startt = 0d0
         brst2st1bb     = 0d0
         brst2st1uu     = 0d0
         brst2st1dd     = 0d0
         brst2st1ee     = 0d0
         brst2st1nunu   = 0d0
         brst2st1tautau = 0d0

      endif

      if(stoptot(1).ne.0d0 .and. stoptot2(1).eq.0d0
     .  .and. flagloop.eq.1d0)then

         brgamma       = gamma/stoptot(1)
         brgammaup     = gammaup/stoptot(1)
         brgammagluino = gammagluino/stoptot(1)

      else

         brgamma       = 0d0
         brgammaup     = 0d0
         brgammagluino = 0d0

      endif

      END

c ------------- ADDITIONAL SUBROUTINES ------------------------------- c
c               contain radiative and 3 body decays                    c
c -------------------------------------------------------------------- c
c ==================================================================== c
c  The decay stop1 -> charm photino
c  K.Hikasa and M.Kobayashi, Phys.Rev.D36 (1987) 724
c ==================================================================== c
      SUBROUTINE NS_hikasakob1(gamma,gammaup,gammagluino)
*
      IMPLICIT NONE
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION u(2,2),v(2,2),z(5,5),zp(5,5)
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION sdmuq,sdmhd2
      DOUBLE PRECISION sdatop,sdabot,sdatau,sdmu
      DOUBLE PRECISION sdmsq,sdmbr,sdmdr
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION rmtc,rmbc,rmtauc
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cmm,smm,cum,sum,cdm,sdm,cem,sem,cnm,snm
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION beta,amh12,amsc1,xktb,xkcb,amx,dl,dr,eps,
     .     f11c,xkub,gamma,gammaup,gammagluino
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_mixmat/u,v,z,zp
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_sfmixang/thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cmm,smm,cum,sum,cdm,sdm,cem,sem,cnm,snm
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_hikasakob/sdmuq,sdmhd2
      COMMON/NS_trilin_mu/sdatop,sdabot,sdatau,sdmu
      COMMON/NS_hikasakob02/sdmsq,sdmbr,sdmdr
      COMMON/NS_runmcalc/rmtc,rmbc,rmtauc
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_tanb/tanbeta_Z

      beta=datan(tanbeta_Z)
      amh12=-sdmhd2
      amsc1=dsqrt(sdmuq**2+
     .      (0.5d0-2d0/3d0*sw**2)*amz**2*DDCOS(2d0*beta))

      xktb = 0.99915d0
      xkcb = 0.041d0
      xkub = 0.0036d0
      amx  = 2d16

      dl=-g2s/16d0/pi**2*log((amx/mw)**2)*xktb*xkcb*scalb**2*
     .    (sdmuq**2+sdmbr**2+sdabot**2+amh12)

      dr=
     .g2s/16d0/pi**2*log((amx/mw)**2)*xktb*xkcb*scalb**2*amt*sdabot

      eps=-(dl*ct+dr*st)/(amsc1**2-ast1**2)

      f11c=-(2d0/3d0)*dsqrt(2d0)*sw*zp(1,1)-dsqrt(2d0)*(0.5d0-
     .     2d0/3d0*sw**2)*(zp(1,2)/cw)
c -- the decay stop1 -> neutralino1 charm --
      if(ast1.gt.amneut(1)) then
         gamma=g2s/16d0/pi*eps**2*f11c**2*ast1*
     .        (1d0-amneut(1)**2/ast1**2)**2
      else
         gamma=0d0
      endif
c -- the decay stop1 -> neutralino1 up --
      if(ast1.gt.amneut(1)) then
         gammaup=gamma*(xkub/xkcb)**2
      else
         gammaup=0d0
      endif
c -- the decay stop1 -> gluino charm
      if(ast1.gt.mgluino) then
         gammagluino = 2d0/3d0*gs2/4d0/pi*eps**2*ast1*
     .        (1d0-mgluino**2/ast1**2)**2
      else
         gammagluino=0d0
      endif
      end

c ==================================================================== c
c  The decay stop2 -> charm photino
c  adapted from K.Hikasa and M.Kobayashi, Phys.Rev.D36 (1987) 724
c ==================================================================== c
      SUBROUTINE NS_hikasakob2(gamma2)
*
      IMPLICIT NONE
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION u(2,2),v(2,2),z(5,5),zp(5,5)
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION sdmuq,sdmhd2
      DOUBLE PRECISION sdatop,sdabot,sdatau,sdmu
      DOUBLE PRECISION sdmsq,sdmbr,sdmdr
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cmm,smm,cum,sum,cdm,sdm,cem,sem,cnm,snm
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION beta,amh12,amsc1,xktb,xkcb,amx,dl,dr,eps,
     .     f11c,gamma2
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/NS_pi/PI,SQR2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_mixmat/u,v,z,zp
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_sfmixang/thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cmm,smm,cum,sum,cdm,sdm,cem,sem,cnm,snm
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_hikasakob/sdmuq,sdmhd2
      COMMON/NS_trilin_mu/sdatop,sdabot,sdatau,sdmu
      COMMON/NS_hikasakob02/sdmsq,sdmbr,sdmdr
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_tanb/tanbeta_Z

      beta=datan(tanbeta_Z)
      amh12=sdmhd2
      amsc1=dsqrt(sdmuq**2+(0.5d0-2d0/3d0*sw**2)*amz**2*
     .     DDCOS(2d0*beta))

      xktb = 0.99915d0
      xkcb = 0.041d0
      amx  = 2d16

      dl=-g2s/16d0/pi**2*log((amx/mw)**2)*xktb*xkcb*scalb**2*
     .   (sdmuq**2+sdmbr**2+sdabot**2+amh12)

      dr=
     .g2s/16d0/pi**2*log((amx/mw)**2)*xktb*xkcb*scalb**2*amt*sdabot

      eps=-(-dl*st+dr*ct)/(amsc1**2-ast2**2)

      f11c=-(2d0/3d0)*dsqrt(2d0)*sw*zp(1,1)-dsqrt(2d0)*(0.5d0-
     .     2d0/3d0*sw**2)*(zp(1,2)/cw)
c -- the decay stop2 -> neutralino1 charm --
      if(ast2.gt.amneut(1)) then
         gamma2=g2s/16d0/pi*eps**2*f11c**2*ast2*(1-amneut(1)**2
     .        /ast2**2)**2
      else
         gamma2=0d0
      endif

      end
c ==================================================================== c
c                           stop 3-body decays                         c
c ==================================================================== c
      SUBROUTINE NS_xintegstop(xintegstopw,xintegstoph,xintegststau,
     .     xintegstsntau,xintegstsmu,xintegstsnmu,xintegstsel,
     .     xintegstsnel,xintegstbsbst,xintegstbbsbt,xintegsttausbnu,
     .     xintegstelsbnu,xintegstupsbdow,xintegst2st1tt,xintegst2st1bb,
     .     xintegst2st1uu,xintegst2st1dd,xintegst2st1ee,
     .     xintegst2st1nunu,xintegst2st1tautau,xintegst2st1startt)
*
      IMPLICIT NONE
      INTEGER nx1t,ny1t,ni,nj
      DOUBLE PRECISION xintegstopw(2,5),xintegstoph(2,5),
     .     xintegststau(2,2),xintegstsntau(2,2),xintegstsmu(2,2),
     .     xintegstsnmu(2),xintegstsel(2,2),xintegstsnel(2),
     .     xintegstbsbst(2,2),xintegstbbsbt(2,2),
     .     xintegsttausbnu(2,2),xintegstelsbnu(2,2),
     .     xintegstupsbdow(2,2)
      DOUBLE PRECISION xintegst2st1tt,xintegst2st1bb,
     .     xintegst2st1uu,xintegst2st1dd,xintegst2st1ee,
     .     xintegst2st1nunu,xintegst2st1tautau,xintegst2st1startt
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION gmsb(2),gmst(2),gmstau(2),gmsmu(2),gmsel(2),
     .          gmsnt(2),gmsnm(2),gmsne(2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),amch
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION NS_ay,NS_by,NS_ax,NS_bx
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION XMU1,XMU2,XMU3,sum1
*
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,amch
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_nx1/nx1t,ny1t
      COMMON/NS_indices/ni,nj
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_stbchiw,NS_stbchihc,NS_stbnustau,NS_stbsnutau,
     .         NS_stbnusmu,NS_stbsnumu,NS_stbnusel,NS_stbsnuel,
     .         NS_stbsbstart,NS_stbbsbt,NS_sttausbnu,NS_stelsbnu,
     .         NS_st2st1tt,NS_st2st1bb,NS_st2st1uu,NS_st2st1dd,
     .         NS_st2st1ee,NS_st2st1nunu,NS_st2st1tautau,NS_st2st1startt

      gmst(1) = ast1
      gmst(2) = ast2
      gmsb(1) = asb1
      gmsb(2) = asb2

      gmstau(1) = astau1
      gmstau(2) = astau2
      gmsmu(1)  = asmu1
      gmsmu(2)  = asmu2
      gmsel(1)  = ase1
      gmsel(2)  = ase2
      gmsnt(1)  = asntau1
      gmsnt(2)  = asntau2
      gmsnm(1)  = asnmu1
      gmsnm(2)  = asnmu2
      gmsne(1)  = asne1
      gmsne(2)  = asne2
c -------------------------------------------------------------------- c
c -------------------------- neutralino_j b W+ ----------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2
         do nj=1,5
            xmu1=amb**2/gmst(ni)**2
            xmu2=amneut(nj)**2/gmst(ni)**2
            xmu3=mw**2/gmst(ni)**2

            if(gmst(ni).gt.(amneut(nj)+mw+amb)) then
               CALL NS_integ2(NS_stbchiw,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xintegstopw(ni,nj)=g2s**2/32d0/(2d0*pi)**3*gmst(ni)*sum1
                               else
               xintegstopw(ni,nj)=0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c -------------------------- neutralino_j b H+ ----------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2
         do nj=1,5
            xmu1=amb**2/gmst(ni)**2
            xmu2=amneut(nj)**2/gmst(ni)**2
            xmu3=amch**2/gmst(ni)**2

            if(gmst(ni).gt.(amneut(nj)+amch+amb)) then
               CALL NS_integ2(NS_stbchihc,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xintegstoph(ni,nj)=g2s**2/32d0/(2d0*pi)**3*gmst(ni)*sum1
            else
               xintegstoph(ni,nj)=0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c ------------------------- b stau neutrino_tau ---------------------- c
c -------------------------------------------------------------------- c

      do ni=1,2
         do nj=1,2
            xmu1=amb**2/gmst(ni)**2
            xmu2=0d0
            xmu3=gmstau(nj)**2/gmst(ni)**2

            if(gmst(ni).gt.(amb+gmstau(nj))) then
               CALL NS_integ2(NS_stbnustau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xintegststau(ni,nj)=g2s**2/32d0/(2d0*pi)**3*gmst(ni)*sum1
            else
               xintegststau(ni,nj)=0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c ------------------------- b sneutrino_tau tau ---------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2
         do nj=1,2
            xmu1=amb**2/gmst(ni)**2
            xmu2=amtau**2/gmst(ni)**2
            xmu3=gmsnt(nj)**2/gmst(ni)**2

            if(gmst(ni).gt.(gmsnt(nj)+amb+amtau)) then
               CALL NS_integ2(NS_stbsnutau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xintegstsntau(ni,nj)=g2s**2/32d0/(2d0*pi)**3*gmst(ni)*
     .                              sum1
            else
               xintegstsntau(ni,nj)=0d0
            endif
         enddo
      enddo

c -------------------------------------------------------------------- c
c ----------------------- b smuon neutrino_mu ------------------------ c
c -------------------------------------------------------------------- c
      do ni=1,2
         do nj=1,2
            xmu1=amb**2/gmst(ni)**2
            xmu2=0d0
            xmu3=gmsmu(nj)**2/gmst(ni)**2

            if(gmst(ni).gt.(gmsmu(nj)+amb)) then
               CALL NS_integ2(NS_stbnusmu,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xintegstsmu(ni,nj)=g2s**2/32d0/(2d0*pi)**3*gmst(ni)*sum1
            else
               xintegstsmu(ni,nj)=0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c ----------------------- b sneutrino_mu muon ------------------------ c
c -------------------------------------------------------------------- c
      do ni=1,2
         xmu1=amb**2/gmst(ni)**2
         xmu2=0d0
         xmu3=gmsnm(1)**2/gmst(ni)**2

         if(gmst(ni).gt.(gmsnm(1)+amb)) then
            CALL NS_integ2(NS_stbsnumu,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                  xmu2,xmu3,nx1t,ny1t,sum1)
            xintegstsnmu(ni)=g2s**2/32d0/(2d0*pi)**3*gmst(ni)*sum1
         else
            xintegstsnmu(ni)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c ----------------------- b selectron neutrino_e --------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2
         do nj=1,2
            xmu1=amb**2/gmst(ni)**2
            xmu2=0d0
            xmu3=gmsel(nj)**2/gmst(ni)**2

            if(gmst(ni).gt.(gmsel(nj)+amb)) then
               CALL NS_integ2(NS_stbnusel,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xintegstsel(ni,nj)=g2s**2/32d0/(2d0*pi)**3*gmst(ni)*sum1
            else
               xintegstsel(ni,nj)=0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c ----------------------- b sneutrino_e electron --------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2
         xmu1=amb**2/gmst(ni)**2
         xmu2=0d0
         xmu3=gmsne(1)**2/gmst(ni)**2

         if(gmst(ni).gt.(gmsne(1)+amb)) then
            CALL NS_integ2(NS_stbsnuel,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                  xmu2,xmu3,nx1t,ny1t,sum1)
            xintegstsnel(ni)=g2s**2/32d0/(2d0*pi)**3*gmst(ni)*sum1
         else
            xintegstsnel(ni)=0d0
         endif
      enddo
c -------------------------------------------------------------------- c
c ------------------------- sbottom_1/2* b top ----------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2
         do nj=1,2
            xmu1=amb**2/gmst(ni)**2
            xmu2=amt**2/gmst(ni)**2
            xmu3=gmsb(nj)**2/gmst(ni)**2

            if(gmst(ni).gt.(amt+amb+gmsb(nj))) then
               CALL NS_integ2(NS_stbsbstart,NS_ax,NS_bx,NS_ay,NS_by,
     .                     xmu1,xmu2,xmu3,nx1t,ny1t,sum1)
              xintegstbsbst(ni,nj)=1d0/32d0/(2d0*pi)**3*gmst(ni)*sum1

            else
               xintegstbsbst(ni,nj)=0d0
            endif

         enddo
      enddo
c -------------------------------------------------------------------- c
c ------------------------ sbottom_1/2 bbar top ---------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2
         do nj=1,2
            xmu1=amb**2/gmst(ni)**2
            xmu2=amt**2/gmst(ni)**2
            xmu3=gmsb(nj)**2/gmst(ni)**2

            if(gmst(ni).gt.(amt+amb+gmsb(nj))) then
               CALL NS_integ2(NS_stbbsbt,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)
               xintegstbbsbt(ni,nj)=1d0/32d0/(2d0*pi)**3*gmst(ni)*sum1
            else
               xintegstbbsbt(ni,nj)=0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c ---------------------- sbottom_1/2 tau+ nu_tau --------------------- c
c -------------------------------------------------------------------- c
       do ni=1,2
         do nj=1,2
            xmu1=0d0

            xmu2=amtau**2/gmst(ni)**2
            xmu3=gmsb(nj)**2/gmst(ni)**2

            if(gmst(ni).gt.(amtau+gmsb(nj))) then
               CALL NS_integ2(NS_sttausbnu,NS_ax,NS_bx,NS_ay,NS_by,
     .                     xmu1,xmu2,xmu3,nx1t,ny1t,sum1)
              xintegsttausbnu(ni,nj)=g2s**2/32d0/(2d0*pi)**3*gmst(ni)*
     .                               sum1
            else
               xintegsttausbnu(ni,nj)=0d0
            endif

         enddo
      enddo
c -------------------------------------------------------------------- c
c --------------------- sbottom_1/2 electron+ nu_e ------------------- c
c -------------------------------------------------------------------- c
      do ni=1,2
         do nj=1,2
            xmu1=0d0
            xmu2=0d0
            xmu3=gmsb(nj)**2/gmst(ni)**2

            if(gmst(ni).gt.gmsb(nj)) then
               CALL NS_integ2(NS_stelsbnu,NS_ax,NS_bx,NS_ay,NS_by,
     .                     xmu1,xmu2,xmu3,nx1t,ny1t,sum1)
               xintegstelsbnu(ni,nj)=g2s**2/32d0/(2d0*pi)**3*gmst(ni)*
     .                               sum1
            else
               xintegstelsbnu(ni,nj)=0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c ----------------------- sbottom_1/2 up downbar --------------------- c
c -------------------------------------------------------------------- c
       do ni=1,2
         do nj=1,2
            if(gmst(ni).gt.gmsb(nj)) then
               xintegstupsbdow(ni,nj)=3d0*xintegstelsbnu(ni,nj)
            else
               xintegstupsbdow(ni,nj)=0d0
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c ---------------------------- stop1* t t ---------------------------- c
c -------------------------------------------------------------------- c
      xmu1=amt**2/gmst(2)**2
      xmu2=amt**2/gmst(2)**2
      xmu3=gmst(1)**2/gmst(2)**2

      if(gmst(2).gt.(gmst(1)+2d0*amt)) then
         CALL NS_integ2(NS_st2st1startt,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .               xmu2,xmu3,nx1t,ny1t,sum1)
         xintegst2st1startt=1d0/32d0/(2d0*pi)**3*gmst(2)*sum1
      else
         xintegst2st1startt=0d0
      endif

c-------------------------------------------------------------------- c
c ---------------------------- stop1 t tbar -------------------------- c
c -------------------------------------------------------------------- c
      xmu1=amt**2/gmst(2)**2
      xmu2=amt**2/gmst(2)**2
      xmu3=gmst(1)**2/gmst(2)**2

      if(gmst(2).gt.(gmst(1)+2d0*amt)) then
         CALL NS_integ2(NS_st2st1tt,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .               xmu2,xmu3,nx1t,ny1t,sum1)
         xintegst2st1tt=1d0/32d0/(2d0*pi)**3*gmst(2)*sum1
      else
         xintegst2st1tt=0d0
      endif
c -------------------------------------------------------------------- c
c ---------------------------- stop1 b bbar -------------------------- c
c -------------------------------------------------------------------- c
      xmu1=amb**2/gmst(2)**2
      xmu2=amb**2/gmst(2)**2
      xmu3=gmst(1)**2/gmst(2)**2

      if((gmst(1)+2d0*amb).lt.gmst(2)) then
         CALL NS_integ2(NS_st2st1bb,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .               xmu2,xmu3,nx1t,ny1t,sum1)
         xintegst2st1bb=g2s**2/32d0/(2d0*pi)**3*gmst(2)*sum1*3d0
      else
         xintegst2st1bb=0d0
      endif

c -------------------------------------------------------------------- c
c --------------------------- stop1 up upbar ------------------------- c
c -------------------------------------------------------------------- c
      xmu1=0d0
      xmu2=0d0
      xmu3=gmst(1)**2/gmst(2)**2

      if(gmst(2).gt.gmst(1)) then
         CALL NS_integ2(NS_st2st1uu,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .               xmu2,xmu3,nx1t,ny1t,sum1)
         xintegst2st1uu=g2s**2/32d0/(2d0*pi)**3*gmst(2)*sum1*3d0
      else
         xintegst2st1uu=0d0
      endif
c -------------------------------------------------------------------- c
c ------------------------- stop1 down downbar ----------------------- c
c -------------------------------------------------------------------- c
      xmu1=0d0
      xmu2=0d0
      xmu3=gmst(1)**2/gmst(2)**2

      if(gmst(2).gt.gmst(1)) then
         CALL NS_integ2(NS_st2st1dd,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .               xmu2,xmu3,nx1t,ny1t,sum1)
         xintegst2st1dd=g2s**2/32d0/(2d0*pi)**3*gmst(2)*sum1*3d0
      else
         xintegst2st1dd=0d0
      endif
c -------------------------------------------------------------------- c
c ----------------------------- stop1 e+ e- -------------------------- c
c -------------------------------------------------------------------- c
      xmu1=0d0
      xmu2=0d0
      xmu3=gmst(1)**2/gmst(2)**2

      if(gmst(2).gt.gmst(1)) then
         CALL NS_integ2(NS_st2st1ee,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .               xmu2,xmu3,nx1t,ny1t,sum1)
         xintegst2st1ee=g2s**2/32d0/(2d0*pi)**3*gmst(2)*sum1
      else
         xintegst2st1ee=0d0
      endif
c -------------------------------------------------------------------- c
c ---------------------------- stop1 nu nubar ------------------------ c
c -------------------------------------------------------------------- c
      xmu1=0d0
      xmu2=0d0
      xmu3=gmst(1)**2/gmst(2)**2

      if(gmst(2).gt.gmst(1)) then
         CALL NS_integ2(NS_st2st1nunu,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .               xmu2,xmu3,nx1t,ny1t,sum1)
         xintegst2st1nunu=g2s**2/32d0/(2d0*pi)**3*gmst(2)*sum1
      else
         xintegst2st1nunu=0d0
      endif
c -------------------------------------------------------------------- c
c --------------------------- stop1 tau+ tau- ------------------------ c
c -------------------------------------------------------------------- c
      xmu1=amtau**2/gmst(2)**2
      xmu2=amtau**2/gmst(2)**2
      xmu3=gmst(1)**2/gmst(2)**2

      if(gmst(2).gt.(gmst(1)+2d0*amtau)) then
         CALL NS_integ2(NS_st2st1tautau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .               xmu2,xmu3,nx1t,ny1t,sum1)
         xintegst2st1tautau=g2s**2/32d0/(2d0*pi)**3*gmst(2)*sum1
      else
         xintegst2st1tautau=0d0
      endif
      end
c ==================================================================== c
c ========================== neutralino_j b W+ ======================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_stbchiw(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,k,l,i
      DOUBLE PRECISION gmst(2),xmuchar(2),dsb(2),dchi(2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION gwtb(2,2),gwntau(2,2),gwnmu(2,2)
      DOUBLE PRECISION ql(5,2),qr(5,2),ol(5,2),or(5,2)
      DOUBLE PRECISION alstor(2,2),akstor(2,2)
      DOUBLE PRECISION abot(2,5),bbot(2,5)
      DOUBLE PRECISION atopr(2,5),btopr(2,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
c ------------- additional (local) variables -------------------------- c
      DOUBLE PRECISION vw,xmuw,xmut,xmub,xmuneut,xmusb(2),x3,
     .x1,x2,y1,y2,y3,dt,stbchiwbb,stbchiwtt,stbchiwchichi,stbchiwchib,
     .stbchiwbt,stbchiwchit,gmsb(2)
C ---------------------------------------------------------------------
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_indices/ni,nj
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/NS_coup20/gwtb,gwntau,gwnmu
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
*
      COMMON/NS_neutsbotbot/abot,bbot
      COMMON/NS_neutstoptop/atopr,btopr
      COMMON/NS_charstopbot/alstor,akstor
*
      gmst(1)=ast1
      gmst(2)=ast2

      vw=1d0/dsqrt(2d0)
      gmsb(1)=asb1
      gmsb(2)=asb2
      xmuw       = mw**2/gmst(ni)**2
      xmut       = amt**2/gmst(ni)**2
      xmub       = amb**2/gmst(ni)**2
      xmuneut    = amneut(nj)**2/gmst(ni)**2
      xmusb(1)     = asb1**2/gmst(ni)**2
      xmusb(2)     = asb2**2/gmst(ni)**2
      xmuchar(1) = amchar(1)**2/gmst(ni)**2
      xmuchar(2) = amchar(2)**2/gmst(ni)**2

      x3=2d0-x1-x2
      y1=(1d0+xmuw-xmuneut-xmub-x3)/2d0
      y2=(1d0-xmuw+xmuneut-xmub-x2)/2d0
      y3=(1d0-xmuw-xmuneut+xmub-x1)/2d0

      dt      = 1d0-x2+xmuneut-xmut
      dsb(1)  = 1d0-x3+xmuw-xmusb(1)
      dsb(2)  = 1d0-x3+xmuw-xmusb(2)
      dchi(1) = 1d0-x1+xmub-xmuchar(1)
      dchi(2) = 1d0-x1+xmub-xmuchar(2)
c -------------------------------------------------------------------- c
c                           sbottom exchange
c -------------------------------------------------------------------- c
      stbchiwbb=0d0

      do k=1,2
         do l=1,2
            if ((gmsb(k)+mw).gt.gmst(ni).and.(gmsb(l)+mw).gt.gmst(ni))
     .then
            stbchiwbb=stbchiwbb+gwtb(ni,k)*gwtb(ni,l)/dsb(k)/dsb(l)*
     .           1d0/2d0*8d0*(
     .           (abot(k,nj)*abot(l,nj)+bbot(k,nj)*bbot(l,nj))*
     .           (-y1*(2d0*y1+xmuneut+xmub)+y1/xmuw*(y2+y3)**2) +
     .           dsqrt(xmub)*xmneut(nj)/gmst(ni)*
     .           (abot(k,nj)*bbot(l,nj)+abot(l,nj)*bbot(k,nj))*
     .           (2d0*y1+xmub+xmuneut-1d0/xmuw*(y2+y3)**2) )
       endif
         enddo
      enddo

c -------------------------------------------------------------------- c
c                             top exchange
c -------------------------------------------------------------------- c
      stbchiwtt=0d0
      if ((amt+amneut(nj)).gt.gmst(ni))then
      stbchiwtt=vw**2/dt**2*(
     .     dsqrt(xmut)*xmneut(nj)/gmst(ni)*atopr(ni,nj)*btopr(ni,nj)*
     .     (-4d0)*(xmub+3d0*y2+2d0*y2**2/xmuw) +
     .     2d0*btopr(ni,nj)**2*xmut*(y1+2d0*y2*y3/xmuw) +
     .     2d0*atopr(ni,nj)**2*(y1*(xmub-xmuw+4d0*y2)+2d0*y3*xmub+
     .     4d0*y2*y3+1d0/xmuw*(4d0*y2**2*y1-2d0*y2*y3*xmub)) )
      endif
c -------------------------------------------------------------------- c
c                           chargino exchange
c -------------------------------------------------------------------- c
      stbchiwchichi=0d0

      do k=1,2
         do l=1,2
      if ((amchar(k)+mb).gt.gmst(ni).and.(amchar(l)+mb).gt.gmst(ni))then
            stbchiwchichi=stbchiwchichi+2d0/dchi(k)/dchi(l)*(
     .           (alstor(ni,k)*alstor(ni,l)*ol(nj,k)*ol(nj,l)+
     .            akstor(ni,k)*akstor(ni,l)*or(nj,k)*or(nj,l))*
     .           (4d0*y3*(y1+y2+y1*y3/xmuw)+y1*(xmuneut-xmuw)
     .            +2d0*xmuneut*y2*(1-y3/xmuw)) +
     .            xmchar(k)*xmchar(l)/gmst(ni)**2*
     .           (alstor(ni,k)*alstor(ni,l)*or(nj,k)*or(nj,l)+
     .            akstor(ni,k)*akstor(ni,l)*ol(nj,k)*ol(nj,l))*
     .           (y1+2d0/xmuw*y2*y3) +
     .           (-3d0)*(y1+y2)*xmneut(nj)/gmst(ni)*(
     .           (alstor(ni,k)*alstor(ni,l)*or(nj,k)*ol(nj,l)+
     .            akstor(ni,k)*akstor(ni,l)*ol(nj,k)*or(nj,l))*
     .            xmchar(k)/gmst(ni) +
     .           (alstor(ni,k)*alstor(ni,l)*or(nj,l)*ol(nj,k)+
     .            akstor(ni,k)*akstor(ni,l)*ol(nj,l)*or(nj,k))*
     .             xmchar(l)/gmst(ni))+
     .           (-1d0)*(xmuneut+2d0*y3**2/xmuw+3d0*y3)*
     .           dsqrt(xmub)*(
     .           (alstor(ni,k)*akstor(ni,l)*or(nj,k)*or(nj,l)+
     .            akstor(ni,k)*alstor(ni,l)*ol(nj,k)*ol(nj,l))*
     .            xmchar(k)/gmst(ni)+
     .           (alstor(ni,k)*akstor(ni,l)*ol(nj,l)*ol(nj,k)+
     .            akstor(ni,k)*alstor(ni,l)*or(nj,l)*or(nj,k))*
     .           xmchar(l)/gmst(ni)) +
     .           3d0*dsqrt(xmub)*xmneut(nj)/gmst(ni)*
     .           xmchar(k)/gmst(ni)*xmchar(l)/gmst(ni)*
     .           (alstor(ni,k)*akstor(ni,l)*ol(nj,l)*or(nj,k)+
     .            akstor(ni,k)*alstor(ni,l)*or(nj,l)*ol(nj,k)) +
     .           3d0*dsqrt(xmub)*xmneut(nj)/gmst(ni)*
     .           (alstor(ni,k)*akstor(ni,l)*ol(nj,k)*or(nj,l)+
     .            akstor(ni,k)*alstor(ni,l)*or(nj,k)*ol(nj,l))*
     .           (xmuw+xmuneut+2d0*y3) )
         endif
         enddo
      enddo


c -------------------------------------------------------------------- c
c                    chargino sbottom interference
c -------------------------------------------------------------------- c
      stbchiwchib=0d0

      do k=1,2
         do i=1,2
      if ((amchar(k)+mb).gt.gmst(ni).and.(gmsb(i)+mw).gt.gmst(ni))then
            stbchiwchib=stbchiwchib+4d0*gwtb(ni,i)/dsqrt(2d0)/
     .           dsb(i)/dchi(k)*(
     .           xmchar(k)/gmst(ni)*xmneut(nj)/gmst(ni)*
     .           (alstor(ni,k)*abot(i,nj)*or(nj,k)+
     .            akstor(ni,k)*bbot(i,nj)*ol(nj,k))*
     .           (y1-y2/xmuw*(y2+y3)+xmub) +
     .           (alstor(ni,k)*abot(i,nj)*ol(nj,k)+
     .            akstor(ni,k)*bbot(i,nj)*or(nj,k))*
     .           ((y2+y3)*(xmuneut*y2-2d0*y1*y3)/xmuw+y1*(2d0*y1+y2
     .            -y3+xmuneut)+xmuneut*y2-xmub*(xmuneut+y3)) +
     .           dsqrt(xmub)*( xmneut(nj)/gmst(ni)*
     .           (abot(i,nj)*akstor(ni,k)*or(nj,k)
     .           +bbot(i,nj)*alstor(ni,k)*ol(nj,k)) +
     .           (abot(i,nj)*akstor(ni,k)*ol(nj,k)
     .           +bbot(i,nj)*alstor(ni,k)*or(nj,k))*xmchar(k)/gmst(ni))*
     .           (1d0/xmuw*y3*(y2+y3)-xmuneut-y1) )
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c                       top sbottom interference
c -------------------------------------------------------------------- c
      stbchiwbt=0d0

      do i=1,2
        if((gmsb(i)+mw).gt.gmst(ni).and.
     .(amt+amneut(nj)).gt.gmst(ni))then
         stbchiwbt=stbchiwbt+4d0*vw*gwtb(ni,i)/dsqrt(2d0)/dt/dsb(i)*(
     .        dsqrt(xmut)*xmneut(nj)/gmst(ni)*btopr(ni,nj)*abot(i,nj)*
     .        (y1+xmub-y2/xmuw*(y2+y3)) +
     .        atopr(ni,nj)*abot(i,nj)*(y1*y2*(1d0+2d0*(y2+y3)/xmuw)+
     .        xmuneut*y2-y1*y3-2d0*y1**2-y1*xmub+xmub*(xmuneut-y3)+
     .        1d0/xmuw*(-xmub*y2*y3-xmub*y3**2)) +
     .        dsqrt(xmub)*dsqrt(xmut)*btopr(ni,nj)*bbot(i,nj)*(-y1
     .        -xmuneut+1d0/xmuw*y3*(y2+y3)) +
     .        dsqrt(xmub)*xmneut(nj)/gmst(ni)*atopr(ni,nj)*bbot(i,nj)*
     .        (y1+xmub-1d0/xmuw*y2*(y2+y3)) )
         endif
      enddo

c -------------------------------------------------------------------- c
c                      chargino top interference
c -------------------------------------------------------------------- c
      stbchiwchit=0d0

      do i=1,2
      if
     .((amchar(i)+mb).gt.gmst(ni).and.(amt+amneut(nj)).gt.gmst(ni))then
         stbchiwchit=stbchiwchit+vw/dt/dchi(i)*(
     .        xmchar(i)/gmst(ni)*xmneut(nj)/gmst(ni)*
     .        atopr(ni,nj)*alstor(ni,i)*or(nj,i)*
     .        (-6d0*y2-4d0*y2**2/xmuw-2d0*xmub) +
     .        atopr(ni,nj)*alstor(ni,i)*ol(nj,i)*2d0*(
     .        y1*(2d0*y3+2d0*y2+4d0*y1-xmuw)+y2*(4d0*y3+xmuneut)
     .        -2d0*y2*(2d0*y1*y3-xmuneut*y2)/xmuw+
     .        2d0*xmub/xmuw*y3**2-xmub*xmuneut+xmub*y3) +
     .        dsqrt(xmut)*xmneut(nj)/gmst(ni)*
     .        btopr(ni,nj)*alstor(ni,i)*ol(nj,i)*
     .        (-6d0)*(y1+y2) +
     .        dsqrt(xmut)*xmchar(i)/gmst(ni)*alstor(ni,i)*btopr(ni,nj)
     .        *or(nj,i)
     .        *(2d0*y1+4d0/xmuw*y2*y3) +
     .        dsqrt(xmub*xmut)*btopr(ni,nj)*akstor(ni,i)*or(nj,i)*
     .        (-6d0*y3-4d0/xmuw*y3**2-2d0*xmuneut) +
     .        6d0*dsqrt(xmut*xmub)*xmchar(i)/gmst(ni)
     .        *xmneut(nj)/gmst(ni)*
     .        btopr(ni,nj)*akstor(ni,i)*ol(nj,i) +
     .        dsqrt(xmub)*xmchar(i)/gmst(ni)*atopr(ni,nj)*akstor(ni,i)
     .        *ol(nj,i)
     .        *(-6d0)*(y1+y3) +
     .        dsqrt(xmub)*xmneut(nj)/gmst(ni)*
     .        atopr(ni,nj)*akstor(ni,i)*or(nj,i)*( 6d0*xmuw+6d0*y3+
     .        6d0*y2+2d0*y1+4d0/xmuw*y2*y3) )
         endif
      enddo

      NS_stbchiw=stbchiwbb+stbchiwtt+stbchiwchichi+2d0*stbchiwchib+
     .           2d0*stbchiwbt+2d0*stbchiwchit

      end

c ==================================================================== c
c ========================== neutralino_j b H+ ======================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_stbchihc(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,k,l
      DOUBLE PRECISION db(2),dchi(2),gmst(2)
      DOUBLE PRECISION gctbr(2,2)
      DOUBLE PRECISION ql(5,2),qr(5,2),ol(5,2),or(5,2)
      DOUBLE PRECISION alstor(2,2),akstor(2,2)
      DOUBLE PRECISION abot(2,5),bbot(2,5)
      DOUBLE PRECISION atopr(2,5),btopr(2,5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),amch
      DOUBLE PRECISION chtbrunr,chtbrunl
c ------------- additional (local) variables -------------------------- c
      DOUBLE PRECISION xmuch,xmut,xmub,xmuw,xmuneut,xmusb(2),
     .xmuchar(2),x3,x1,x2,y1,y2,dt,stbchihcbb,stbchihctt,stbchihchichi,
     .stbchihcchib,stbchihcbt,stbchihcchit,gmsb(2)
C ---------------------------------------------------------------------
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_indices/ni,nj
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,amch
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_coup3/ql,qr,ol,or
      COMMON/NS_hcsbotstop/gctbr
      COMMON/NS_neutsbotbot/abot,bbot
      COMMON/NS_neutstoptop/atopr,btopr
      COMMON/NS_higgschudb/chtbrunr,chtbrunl
      COMMON/NS_charstopbot/alstor,akstor

      gmst(1)=ast1
      gmst(2)=ast2
      gmsb(1)=asb1
      gmsb(2)=asb2
      xmuch      = amch**2/gmst(ni)**2
      xmut       = amt**2/gmst(ni)**2
      xmub       = amb**2/gmst(ni)**2
      xmuw       = amw**2/gmst(ni)**2
      xmuneut    = amneut(nj)**2/gmst(ni)**2
      xmusb(1)     = asb1**2/gmst(ni)**2
      xmusb(2)     = asb2**2/gmst(ni)**2
      xmuchar(1) = amchar(1)**2/gmst(ni)**2
      xmuchar(2) = amchar(2)**2/gmst(ni)**2

      x3=2d0-x1-x2
      y1=(1d0+xmuch-xmuneut-xmub-x3)/2d0
      y2=(1d0-xmuch+xmuneut-xmub-x2)/2d0

      db(1)   = 1d0-x3+xmuch-xmusb(1)
      db(2)   = 1d0-x3+xmuch-xmusb(2)
      dt      = 1d0-x2+xmuneut-xmut
      dchi(1) = 1d0-x1+xmub-xmuchar(1)
      dchi(2) = 1d0-x1+xmub-xmuchar(2)
c -------------------------------------------------------------------- c
c                           sbottom exchange
c -------------------------------------------------------------------- c
      stbchihcbb=0d0

      do k=1,2
         do l=1,2
      if((gmsb(k)+amch).gt.gmst(ni).and.(gmsb(l)+amch).gt.gmst(ni))then
           stbchihcbb= stbchihcbb+gctbr(ni,k)*gctbr(ni,l)*xmuw/db(k)
     .           /db(l)*(
     .           (abot(k,nj)*abot(l,nj)+bbot(k,nj)*bbot(l,nj))*
     .           2d0*y1 +
     .           (abot(k,nj)*bbot(l,nj)+bbot(k,nj)*abot(l,nj))*
     .           dsqrt(xmub)*xmneut(nj)/gmst(ni)*(-2d0) )
       endif
        enddo
      enddo

c -------------------------------------------------------------------- c
c                              top exchange
c -------------------------------------------------------------------- c
      stbchihctt=0d0
      if ((amt+amneut(nj)).gt.gmst(ni))then
      stbchihctt=1d0/dt**2*(
     .    dsqrt(xmut)*xmneut(nj)/gmst(ni)*
     .    (chtbrunr**2+chtbrunl**2)*atopr(ni,nj)*btopr(ni,nj)*
     .    2d0*(2d0*y1-x1) +
     .    (chtbrunr**2*atopr(ni,nj)**2+chtbrunl**2*btopr(ni,nj)**2)*
     .    (-2d0*(xmuneut*x1-xmuneut*y1+y1)+x1*x2) +
     .    (chtbrunr**2*btopr(ni,nj)**2+chtbrunl**2*atopr(ni,nj)**2)*
     .    xmut*2d0*y1 +
     .    chtbrunr*chtbrunl*(atopr(ni,nj)**2+btopr(ni,nj)**2)*
     .    dsqrt(xmub*xmut)*2d0*(x2-2d0*xmuneut) +
     .    chtbrunr*chtbrunl*atopr(ni,nj)*btopr(ni,nj)*
     .    dsqrt(xmub)*xmneut(nj)/gmst(ni)*4d0*(x2-xmuneut-xmut-1d0) )
      endif
c -------------------------------------------------------------------- c
c                           chargino exchange
c -------------------------------------------------------------------- c
      stbchihchichi=0d0

               do k=1,2
         do l=1,2
      if ((amchar(k)+mb).gt.gmst(ni).and.(amchar(l)+mb).gt.gmst(ni))
     .then
            stbchihchichi=stbchihchichi+1d0/dchi(k)/dchi(l)*(
     .          (alstor(ni,k)*alstor(ni,l)*ql(nj,k)*ql(nj,l)+
     .           akstor(ni,k)*akstor(ni,l)*qr(nj,k)*qr(nj,l))*
     .          (x1*x2+2d0*y1*(xmub-1d0)-2d0*x2*xmub) +
     .          xmchar(k)/gmst(ni)*xmneut(nj)/gmst(ni)*
     .          (alstor(ni,k)*alstor(ni,l)*ql(nj,l)*qr(nj,k)+
     .           akstor(ni,k)*akstor(ni,l)*ql(nj,k)*qr(nj,l))*
     .          2d0*(y1+y2) +
     .          xmchar(l)/gmst(ni)*xmneut(nj)/gmst(ni)*
     .          (alstor(ni,k)*alstor(ni,l)*ql(nj,k)*qr(nj,l)+
     .           akstor(ni,k)*akstor(ni,l)*ql(nj,l)*qr(nj,k))*
     .          2d0*(y1+y2) +
     .          xmchar(k)*xmchar(l)/gmst(ni)**2*
     .          (alstor(ni,k)*alstor(ni,l)*qr(nj,k)*qr(nj,l)+
     .           akstor(ni,k)*akstor(ni,l)*ql(nj,k)*ql(nj,l))*2d0*y1 +
     .          dsqrt(xmub)*(xmchar(k)/gmst(ni)*
     .          (alstor(ni,l)*akstor(ni,k)*ql(nj,k)*ql(nj,l)+
     .           akstor(ni,l)*alstor(ni,k)*qr(nj,k)*qr(nj,l)) +
     .          (alstor(ni,l)*akstor(ni,k)*qr(nj,k)*qr(nj,l)+
     .           akstor(ni,l)*alstor(ni,k)*ql(nj,k)*ql(nj,l))*
     .          xmchar(l)/gmst(ni) )*(2d0*y1-x2) +
     .          dsqrt(xmub)*xmchar(k)*xmchar(l)/gmst(ni)**2*
     .          xmneut(nj)/gmst(ni)*
     .          (alstor(ni,l)*akstor(ni,k)*ql(nj,k)*qr(nj,l)+
     .           akstor(ni,l)*alstor(ni,k)*qr(nj,k)*ql(nj,l))*(-2d0) +
     .          dsqrt(xmub)*xmneut(nj)/gmst(ni)*
     .          (alstor(ni,l)*akstor(ni,k)*qr(nj,k)*ql(nj,l)+
     .           akstor(ni,l)*alstor(ni,k)*ql(nj,k)*qr(nj,l))*
     .          2d0*(-1d0-xmub+x1) )
         endif
         enddo
      enddo


c -------------------------------------------------------------------- c
c                    chargino sbottom interference
c -------------------------------------------------------------------- c
      stbchihcchib=0d0

      do k=1,2
         do l=1,2

      if ((amchar(l)+mb).gt.gmst(ni).and.(gmsb(k)+amch).gt.gmst(ni))then
            stbchihcchib=stbchihcchib+
     .          gctbr(ni,k)*dsqrt(xmuw)/db(k)/dchi(l)*(
     .          (alstor(ni,l)*abot(k,nj)*ql(nj,l)+
     .           akstor(ni,l)*bbot(k,nj)*qr(nj,l))*
     .          xmneut(nj)/gmst(ni)*(x1-2d0*xmub) +
     .          (akstor(ni,l)*bbot(k,nj)*ql(nj,l)+
     .           alstor(ni,l)*abot(k,nj)*qr(nj,l))*xmchar(l)/gmst(ni)*
     .          2d0*y1 +
     .          (alstor(ni,l)*bbot(k,nj)*ql(nj,l)+
     .           akstor(ni,l)*abot(k,nj)*qr(nj,l))*dsqrt(xmub)*
     .          (2d0*y1-x2) +
     .          (akstor(ni,l)*abot(k,nj)*ql(nj,l)+
     .           alstor(ni,l)*bbot(k,nj)*qr(nj,l))*dsqrt(xmub)*
     .          xmchar(l)/gmst(ni)*xmneut(nj)/gmst(ni)*(-2d0) )
            endif
         enddo
      enddo

c -------------------------------------------------------------------- c
c                       top sbottom interference
c -------------------------------------------------------------------- c
      stbchihcbt=0d0

      do k=1,2
         if((amt+amneut(nj)).gt.gmst(ni).and.(gmsb(k)+amch).gt.gmst(ni))
     .then
         stbchihcbt=stbchihcbt-gctbr(ni,k)*dsqrt(xmuw)/db(k)/dt*(
     .        dsqrt(xmut)*
     .        (atopr(ni,nj)*abot(k,nj)*chtbrunl+
     .         btopr(ni,nj)*bbot(k,nj)*chtbrunr)*2d0*y1 +
     .        (btopr(ni,nj)*abot(k,nj)*chtbrunl+
     .         atopr(ni,nj)*bbot(k,nj)*chtbrunr)*xmneut(nj)/gmst(ni)*
     .        (2d0*y1-x1) +
     .        (btopr(ni,nj)*bbot(k,nj)*chtbrunl+
     .         atopr(ni,nj)*abot(k,nj)*chtbrunr)*dsqrt(xmub)*
     .        (x2-2d0*xmuneut) +
     .        (btopr(ni,nj)*abot(k,nj)*chtbrunr+
     .         atopr(ni,nj)*bbot(k,nj)*chtbrunl)*dsqrt(xmub*xmut)*
     .        xmneut(nj)/gmst(ni)*(-2d0) )
         endif
      enddo

c -------------------------------------------------------------------- c
c                      chargino top interference
c -------------------------------------------------------------------- c
      stbchihcchit=0d0

      do k=1,2
      if
     .((amt+amneut(nj)).gt.gmst(ni).and.(amchar(k)+mb).gt.gmst(ni))then
         stbchihcchit=stbchihcchit-1d0/dchi(k)/dt*(
     .        dsqrt(xmut)*xmneut(nj)/gmst(ni)*
     .        (alstor(ni,k)*atopr(ni,nj)*ql(nj,k)*chtbrunl+
     .         akstor(ni,k)*btopr(ni,nj)*qr(nj,k)*chtbrunr)*
     .        (x1-2d0*xmub) +
     .        (alstor(ni,k)*btopr(ni,nj)*ql(nj,k)*chtbrunl+
     .         akstor(ni,k)*atopr(ni,nj)*qr(nj,k)*chtbrunr)*
     .        (xmuneut*(x1-2d0*xmub)-x2*(x1-xmub)+2d0*y1) +
     .        dsqrt(xmut)*xmchar(k)/gmst(ni)*
     .        (akstor(ni,k)*btopr(ni,nj)*ql(nj,k)*chtbrunr+
     .         alstor(ni,k)*atopr(ni,nj)*qr(nj,k)*chtbrunl)*2d0*y1
     .        +xmchar(k)/gmst(ni)*xmneut(nj)/gmst(ni)*
     .        (akstor(ni,k)*atopr(ni,nj)*ql(nj,k)*chtbrunr+
     .         alstor(ni,k)*btopr(ni,nj)*qr(nj,k)*chtbrunl)*
     .        (2d0*y1-x1) +
     .        (akstor(ni,k)*btopr(ni,nj)*ql(nj,k)*chtbrunl+
     .         alstor(ni,k)*atopr(ni,nj)*qr(nj,k)*chtbrunr)*
     .        dsqrt(xmub)*xmchar(k)/gmst(ni)*(x2-2d0*xmuneut) +
     .        (alstor(ni,k)*btopr(ni,nj)*ql(nj,k)*chtbrunr+
     .         akstor(ni,k)*atopr(ni,nj)*qr(nj,k)*chtbrunl)*
     .        dsqrt(xmub*xmut)*(2d0*y1-x2) +
     .        (akstor(ni,k)*atopr(ni,nj)*ql(nj,k)*chtbrunl+
     .         alstor(ni,k)*btopr(ni,nj)*qr(nj,k)*chtbrunr)*
     .        dsqrt(xmut*xmub)*xmchar(k)/gmst(ni)*
     .        xmneut(nj)/gmst(ni)*(-2d0) +
     .        (alstor(ni,k)*atopr(ni,nj)*ql(nj,k)*chtbrunr+
     .         akstor(ni,k)*btopr(ni,nj)*qr(nj,k)*chtbrunl)*
     .        dsqrt(xmub)*xmneut(nj)/gmst(ni)*
     .        (xmuch-xmub-xmuneut+1d0) )
         endif
      enddo

      NS_stbchihc=stbchihcbb+stbchihctt+stbchihchichi+2d0*stbchihcchib
     .            +2d0*stbchihcbt+2d0*stbchihcchit
      end
C ---------------------------------------------------------------------C
c ==================================================================== c
c ========================= b stau neutrino_tau ====================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_stbnustau(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,i,k,l
      DOUBLE PRECISION dchi(2),gmst(2),xmuchar(2),xmustau(2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION alstor(2,2),akstor(2,2),bltau(2,2)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION xmub,x1,x2
*
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_indices/ni,nj
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_charstopbot/alstor,akstor

      gmst(1)=ast1
      gmst(2)=ast2

      xmuchar(1) = amchar(1)**2/gmst(ni)**2
      xmuchar(2) = amchar(2)**2/gmst(ni)**2
      xmustau(1) = astau1**2/gmst(ni)**2
      xmustau(2) = astau2**2/gmst(ni)**2
      xmub       = amb**2/gmst(ni)**2

      dchi(1)=1d0-x1-xmuchar(1)+xmub
      dchi(2)=1d0-x1-xmuchar(2)+xmub

      do i=1,2
         bltau(1,i)=0d0
         bltau(2,i)=0d0
      enddo

      NS_stbnustau=0d0

      do k=1,2
         do l=1,2
       if ((amchar(k)+mb).gt.gmst(ni).and.(amchar(l)+mb).gt.gmst(ni))
     .then
           NS_stbnustau=NS_stbnustau+1d0/dchi(k)/dchi(l)*(
     .        (alstor(ni,k)*alstor(ni,l)*altau(nj,k)*altau(nj,l)
     .        +akstor(ni,k)*akstor(ni,l)*bltau(nj,k)*bltau(nj,l))*
     .        xmchar(k)*xmchar(l)/gmst(ni)**2*(x1+x2-1d0+xmustau(nj)
     .         -xmub)
     .       +(alstor(ni,k)*alstor(ni,l)*bltau(nj,k)*bltau(nj,l)
     .        +akstor(ni,k)*akstor(ni,l)*altau(nj,k)*altau(nj,l))*
     .        ((1d0-x1)*(1d0-x2)-xmustau(nj)+xmub*(xmustau(nj)+x1-x2
     .         -xmub))
     .       +(alstor(ni,k)*akstor(ni,l)*bltau(nj,k)*bltau(nj,l)
     .        +akstor(ni,k)*alstor(ni,l)*altau(nj,k)*altau(nj,l))*
     .        dsqrt(xmub)*xmchar(l)/gmst(ni)*(-1d0-xmub+xmustau(nj)+x1)
     .       +(alstor(ni,k)*akstor(ni,l)*altau(nj,k)*altau(nj,l)
     .        +akstor(ni,k)*alstor(ni,l)*bltau(nj,k)*bltau(nj,l))*
     .        dsqrt(xmub)*xmchar(k)/gmst(ni)
     .       *(-1d0-xmub+xmustau(nj)+x1))
            endif
         enddo
      enddo

      end
c ==================================================================== c
c ========================= b sneutrino_tau tau ====================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_stbsnutau(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,k,l
*
      DOUBLE PRECISION dchi(2),gmst(2),xmusn(2),xmuchar(2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION alstor(2,2),akstor(2,2)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION xmub,xmutau,x1,x2
*
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_indices/ni,nj
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_charstopbot/alstor,akstor
*
      gmst(1)=ast1
      gmst(2)=ast2

      xmuchar(1)= amchar(1)**2/gmst(ni)**2
      xmuchar(2)= amchar(2)**2/gmst(ni)**2
      xmusn(1)  = asntau1**2/gmst(ni)**2
      xmusn(2)  = asntau2**2/gmst(ni)**2
      xmub      = amb**2/gmst(ni)**2
      xmutau    = amtau**2/gmst(ni)**2

      dchi(1)=1d0-x1-xmuchar(1)+xmub
      dchi(2)=1d0-x1-xmuchar(2)+xmub

      NS_stbsnutau=0d0

      do k=1,2
         do l=1,2
       if ((amchar(k)+mb).gt.gmst(ni).and.(amchar(l)+mb).gt.gmst(ni))
     .then
            NS_stbsnutau=NS_stbsnutau+1d0/dchi(k)/dchi(l)*(
     .        (alstor(ni,k)*alstor(ni,l)*blsnt(nj,k)*blsnt(nj,l)+
     .         akstor(ni,k)*akstor(ni,l)*alsnt(nj,k)*alsnt(nj,l))*
     .        xmchar(k)*xmchar(l)/gmst(ni)**2*(x1+x2-1d0+xmusn(nj)
     .        -xmub-xmutau)
     .       +(alstor(ni,k)*alstor(ni,l)*alsnt(nj,k)*alsnt(nj,l)+
     .         akstor(ni,k)*akstor(ni,l)*blsnt(nj,k)*blsnt(nj,l))*
     .        ((1d0-x1)*(1d0-x2)-xmusn(nj)+xmub*(xmusn(nj)+x1-x2
     .        -xmub-xmutau)+xmutau)
     .       +(alstor(ni,k)*alstor(ni,l)*alsnt(nj,k)*blsnt(nj,l)+
     .         akstor(ni,k)*akstor(ni,l)*blsnt(nj,k)*alsnt(nj,l))*
     .        dsqrt(xmutau)*xmchar(l)/gmst(ni)*(-2d0*xmub+x1)
     .       +(alstor(ni,k)*akstor(ni,l)*alsnt(nj,k)*alsnt(nj,l)+
     .         akstor(ni,k)*alstor(ni,l)*blsnt(nj,k)*blsnt(nj,l))*
     .        dsqrt(xmub)*xmchar(l)/gmst(ni)*
     .        (-1d0-xmub-xmutau+xmusn(nj)+x1)
     .       +(alstor(ni,k)*alstor(ni,l)*blsnt(nj,k)*alsnt(nj,l)+
     .         akstor(ni,k)*akstor(ni,l)*alsnt(nj,k)*blsnt(nj,l))*
     .        dsqrt(xmutau)*xmchar(k)/gmst(ni)*(-2d0*xmub+x1)
     .       +(alstor(ni,k)*akstor(ni,l)*blsnt(nj,k)*blsnt(nj,l)+
     .         akstor(ni,k)*alstor(ni,l)*alsnt(nj,k)*alsnt(nj,l))*
     .        dsqrt(xmub)*xmchar(k)/gmst(ni)*
     .        (-1d0-xmub-xmutau+xmusn(nj)+x1)
     .       +(alstor(ni,k)*akstor(ni,l)*alsnt(nj,k)*blsnt(nj,l)+
     .         akstor(ni,k)*alstor(ni,l)*blsnt(nj,k)*alsnt(nj,l))*
     .        dsqrt(xmub*xmutau)*(-2d0-2d0*xmub+2d0*x1)
     .       +(alstor(ni,k)*akstor(ni,l)*blsnt(nj,k)*alsnt(nj,l)+
     .         akstor(ni,k)*alstor(ni,l)*alsnt(nj,k)*blsnt(nj,l))*
     .       dsqrt(xmub*xmutau)*xmchar(k)*xmchar(l)/gmst(ni)**2*(-2d0))
            endif
         enddo
      enddo

      end
c ==================================================================== c
c ======================== b smuon neutrino_mu ======================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_stbnusmu(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,i,k,l
      DOUBLE PRECISION dchi(2),gmst(2),xmusel(2),xmuchar(2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION alstor(2,2),akstor(2,2),blmu(2,2)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION xmub,x1,x2
*
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_indices/ni,nj
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_charstopbot/alstor,akstor

      gmst(1)=ast1
      gmst(2)=ast2

      xmuchar(1) = amchar(1)**2/gmst(ni)**2
      xmuchar(2) = amchar(2)**2/gmst(ni)**2
      xmusel(1)  = asmu1**2/gmst(ni)**2
      xmusel(2)  = asmu2**2/gmst(ni)**2
      xmub       = amb**2/gmst(ni)**2

      dchi(1)=1d0-x1-xmuchar(1)+xmub
      dchi(2)=1d0-x1-xmuchar(2)+xmub

      do i=1,2
         blmu(1,i)=0d0
         blmu(2,i)=0d0
      enddo

      NS_stbnusmu=0d0

      do k=1,2
         do l=1,2
       if ((amchar(k)+mb).gt.gmst(ni).and.(amchar(l)+mb).gt.gmst(ni))
     .then
            NS_stbnusmu=NS_stbnusmu+1d0/dchi(k)/dchi(l)*(
     .        (alstor(ni,k)*alstor(ni,l)*almu(nj,k)*almu(nj,l)
     .        +akstor(ni,k)*akstor(ni,l)*blmu(nj,k)*blmu(nj,l))*
     .         xmchar(k)*xmchar(l)/gmst(ni)**2*(x1+x2-1d0+xmusel(nj)
     .         -xmub)
     .       +(alstor(ni,k)*alstor(ni,l)*blmu(nj,k)*blmu(nj,l)
     .        +akstor(ni,k)*akstor(ni,l)*almu(nj,k)*almu(nj,l))*
     .        ((1d0-x1)*(1d0-x2)-xmusel(nj)+xmub*(xmusel(nj)+x1-x2
     .         -xmub))
     .       +(alstor(ni,k)*akstor(ni,l)*blmu(nj,k)*blmu(nj,l)
     .        +akstor(ni,k)*alstor(ni,l)*almu(nj,k)*almu(nj,l))*
     .        dsqrt(xmub)*xmchar(l)/gmst(ni)*(-1d0-xmub+xmusel(nj)+x1)
     .       +(alstor(ni,k)*akstor(ni,l)*almu(nj,k)*almu(nj,l)
     .        +akstor(ni,k)*alstor(ni,l)*blmu(nj,k)*blmu(nj,l))*
     .        dsqrt(xmub)*xmchar(k)/gmst(ni)
     .       *(-1d0-xmub+xmusel(nj)+x1) )
            endif
         enddo
      enddo

      end
c ==================================================================== c
c ======================= b sneutrino_mu muon ======================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_stbsnumu(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,L,K
      DOUBLE PRECISION dchi(2),gmst(2),xmusn(2),xmuchar(2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION alstor(2,2),akstor(2,2)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION xmub,x1,x2
*
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_indices/ni,nj
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_charstopbot/alstor,akstor

      gmst(1)=ast1
      gmst(2)=ast2

      xmuchar(1)= amchar(1)**2/gmst(ni)**2
      xmuchar(2)= amchar(2)**2/gmst(ni)**2
      xmusn(1)  = asnmu1**2/gmst(ni)**2
      xmusn(2)  = asnmu2**2/gmst(ni)**2
      xmub      = amb**2/gmst(ni)**2

      dchi(1)=1d0-x1-xmuchar(1)+xmub
      dchi(2)=1d0-x1-xmuchar(2)+xmub

      NS_stbsnumu=0d0

      do k=1,2
         do l=1,2
       if ((amchar(k)+mb).gt.gmst(ni).and.(amchar(l)+mb).gt.gmst(ni))
     .then
            NS_stbsnumu=NS_stbsnumu+1d0/dchi(k)/dchi(l)*(
     .        (alstor(ni,k)*alstor(ni,l)*blsnm(1,k)*blsnm(1,l)+
     .         akstor(ni,k)*akstor(ni,l)*alsnm(1,k)*alsnm(1,l))*
     .        xmchar(k)*xmchar(l)/gmst(ni)**2*(x1+x2-1d0+xmusn(1)
     .        -xmub)
     .       +(alstor(ni,k)*alstor(ni,l)*alsnm(1,k)*alsnm(1,l)+
     .         akstor(ni,k)*akstor(ni,l)*blsnm(1,k)*blsnm(1,l))*
     .        ((1d0-x1)*(1d0-x2)-xmusn(1)+xmub*(xmusn(1)+x1-x2
     .        -xmub))
     .       +(alstor(ni,k)*akstor(ni,l)*alsnm(1,k)*alsnm(1,l)+
     .         akstor(ni,k)*alstor(ni,l)*blsnm(1,k)*blsnm(1,l))*
     .        dsqrt(xmub)*xmchar(l)/gmst(ni)*(-1d0-xmub+xmusn(1)+x1)
     .       +(alstor(ni,k)*akstor(ni,l)*blsnm(1,k)*blsnm(1,l)+
     .         akstor(ni,k)*alstor(ni,l)*alsnm(1,k)*alsnm(1,l))*
     .        dsqrt(xmub)*xmchar(k)/gmst(ni)*(-1d0-xmub+xmusn(1)+x1) )
            endif
         enddo
      enddo

      end
c ==================================================================== c
c ======================== b selectron neutrino_e ==================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_stbnusel(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,i,k,l
      DOUBLE PRECISION dchi(2),gmst(2),xmusel(2),xmuchar(2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION alstor(2,2),akstor(2,2),ble(2,2)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION xmub,x1,x2
*
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_indices/ni,nj
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_charstopbot/alstor,akstor

      gmst(1)=ast1
      gmst(2)=ast2

      xmuchar(1) = amchar(1)**2/gmst(ni)**2
      xmuchar(2) = amchar(2)**2/gmst(ni)**2
      xmusel(1)  = ase1**2/gmst(ni)**2
      xmusel(2)  = ase2**2/gmst(ni)**2
      xmub       = amb**2/gmst(ni)**2

      dchi(1)=1d0-x1-xmuchar(1)+xmub
      dchi(2)=1d0-x1-xmuchar(2)+xmub

      do i=1,2
         ble(1,i)=0d0
         ble(2,i)=0d0
      enddo

      NS_stbnusel=0d0

      do k=1,2
         do l=1,2
       if ((amchar(k)+mb).gt.gmst(ni).and.(amchar(l)+mb).gt.gmst(ni))
     .then
            NS_stbnusel=NS_stbnusel+1d0/dchi(k)/dchi(l)*(
     .        (alstor(ni,k)*alstor(ni,l)*ale(nj,k)*ale(nj,l)
     .        +akstor(ni,k)*akstor(ni,l)*ble(nj,k)*ble(nj,l))*
     .         xmchar(k)*xmchar(l)/gmst(ni)**2*(x1+x2-1d0+xmusel(nj)
     .         -xmub)
     .       +(alstor(ni,k)*alstor(ni,l)*ble(nj,k)*ble(nj,l)
     .        +akstor(ni,k)*akstor(ni,l)*ale(nj,k)*ale(nj,l))*
     .        ((1d0-x1)*(1d0-x2)-xmusel(nj)+xmub*(xmusel(nj)+x1-x2
     .         -xmub))
     .       +(alstor(ni,k)*akstor(ni,l)*ble(nj,k)*ble(nj,l)
     .        +akstor(ni,k)*alstor(ni,l)*ale(nj,k)*ale(nj,l))*
     .        dsqrt(xmub)*xmchar(l)/gmst(ni)*(-1d0-xmub+xmusel(nj)+x1)
     .       +(alstor(ni,k)*akstor(ni,l)*ale(nj,k)*ale(nj,l)
     .        +akstor(ni,k)*alstor(ni,l)*ble(nj,k)*ble(nj,l))*
     .        dsqrt(xmub)*xmchar(k)/gmst(ni)
     .       *(-1d0-xmub+xmusel(nj)+x1) )
            endif
         enddo
      enddo

      end
c ==================================================================== c
c ======================= b sneutrino_e electron ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_stbsnuel(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,L,K
      DOUBLE PRECISION dchi(2),gmst(2),gmsn(2),xmusn(2),xmuchar(2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION alstor(2,2),akstor(2,2)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION xmub,x1,x2
*
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_indices/ni,nj
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_charstopbot/alstor,akstor

      gmst(1)=ast1
      gmst(2)=ast2

      xmuchar(1)= amchar(1)**2/gmst(ni)**2
      xmuchar(2)= amchar(2)**2/gmst(ni)**2
      xmusn(1)  = asne1**2/gmst(ni)**2
      xmusn(2)  = asne2**2/gmst(ni)**2
      xmub      = amb**2/gmst(ni)**2

      dchi(1)=1d0-x1-xmuchar(1)+xmub
      dchi(2)=1d0-x1-xmuchar(2)+xmub

      NS_stbsnuel=0d0

      do k=1,2
         do l=1,2
       if ((amchar(k)+mb).gt.gmst(ni).and.(amchar(l)+mb).gt.gmst(ni))
     .then
            NS_stbsnuel=NS_stbsnuel+1d0/dchi(k)/dchi(l)*(
     .        (alstor(ni,k)*alstor(ni,l)*blsne(1,k)*blsne(1,l)+
     .         akstor(ni,k)*akstor(ni,l)*alsne(1,k)*alsne(1,l))*
     .        xmchar(k)*xmchar(l)/gmst(ni)**2*(x1+x2-1d0+xmusn(1)
     .        -xmub)
     .       +(alstor(ni,k)*alstor(ni,l)*alsne(1,k)*alsne(1,l)+
     .         akstor(ni,k)*akstor(ni,l)*blsne(1,k)*blsne(1,l))*
     .        ((1d0-x1)*(1d0-x2)-xmusn(1)+xmub*(xmusn(1)+x1-x2
     .        -xmub))
     .       +(alstor(ni,k)*akstor(ni,l)*alsne(1,k)*alsne(1,l)+
     .         akstor(ni,k)*alstor(ni,l)*blsne(1,k)*blsne(1,l))*
     .        dsqrt(xmub)*xmchar(l)/gmst(ni)*(-1d0-xmub+xmusn(1)+x1)
     .       +(alstor(ni,k)*akstor(ni,l)*blsne(1,k)*blsne(1,l)+
     .         akstor(ni,k)*alstor(ni,l)*alsne(1,k)*alsne(1,l))*
     .        dsqrt(xmub)*xmchar(k)/gmst(ni)*(-1d0-xmub+xmusn(1)+x1) )
            endif
         enddo
      enddo

      end
c ==================================================================== c
c ========================== b sbottom_1/2* top ====================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_stbsbstart(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj,l,k
      DOUBLE PRECISION neutneut
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION dchi(2),gmst(2),xmuchar(2),gmsb(2),
     .xmusb(2),dneut(5),xmuneut(5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2),alstor(2,2),akstor(2,2)
      DOUBLE PRECISION abot(2,5),bbot(2,5)
      DOUBLE PRECISION atopr(2,5),btopr(2,5)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION
     .gul(2),gur(2),gdl(2),gdr(2),gtl(2),gtr(2),gbl(2),gbr(2)
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION xmut,XMUB,xmugl,x1,x2,dgl,charchar,gluiglui
     .,charneut
*
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_indices/ni,nj
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor
      COMMON/NS_neutsbotbot/abot,bbot
      COMMON/NS_neutstoptop/atopr,btopr
*
      gmst(1)=ast1
      gmst(2)=ast2
      gmsb(1)=asb1
      gmsb(2)=asb2
      do k=1,5
         xmuneut(k) = amneut(k)**2/gmst(ni)**2
      enddo
      xmuchar(1) = amchar(1)**2/gmst(ni)**2
      xmuchar(2) = amchar(2)**2/gmst(ni)**2
      xmusb(1)   = gmsb(1)**2/gmst(ni)**2
      xmusb(2)   = gmsb(2)**2/gmst(ni)**2
      xmut       = amt**2/gmst(ni)**2
      xmub       = amb**2/gmst(ni)**2
      xmugl      = mgluino**2/gmst(ni)**2

      dchi(1)=1d0-x1-xmuchar(1)+xmub
      dchi(2)=1d0-x1-xmuchar(2)+xmub

      do k=1,5
         dneut(k)=1d0-x2-xmuneut(k)+xmut
      enddo

      dgl = 1d0-x2-xmugl+xmut
c -------------------------------------------------------------------- c
c                           chargino exchange
c -------------------------------------------------------------------- c
      charchar=0d0

      do k=1,2
         do l=1,2

       if ((amchar(k)+mb).gt.gmst(ni).and.(amchar(l)+mb).gt.gmst(ni))
     .then
            charchar=charchar+3d0*g2s**2/dchi(k)/dchi(l)*(
     .           xmchar(k)*xmchar(l)/gmst(ni)**2*
     .           (-1d0+x1+x2-xmub-xmut+xmusb(nj))*
     .           (alstor(ni,k)*alstor(ni,l)*alsbot(nj,k)*alsbot(nj,l)+
     .            akstor(ni,k)*akstor(ni,l)*aksbot(nj,k)*aksbot(nj,l))
     .           +(xmub*(-x2+xmusb(nj)-xmub+x1-xmut)-x2+x1*x2-x1+1d0
     .             -xmusb(nj)+xmut)*
     .           (alstor(ni,k)*alstor(ni,l)*aksbot(nj,k)*aksbot(nj,l)+
     .            akstor(ni,k)*akstor(ni,l)*alsbot(nj,k)*alsbot(nj,l))
     .           +dsqrt(xmut)*xmchar(k)/gmst(ni)*(x1-2d0*xmub)*
     .           (alstor(ni,k)*alstor(ni,l)*alsbot(nj,k)*aksbot(nj,l)+
     .            akstor(ni,k)*akstor(ni,l)*alsbot(nj,l)*aksbot(nj,k))
     .           +dsqrt(xmut)*xmchar(l)/gmst(ni)*(x1-2d0*xmub)*
     .           (alstor(ni,k)*alstor(ni,l)*alsbot(nj,l)*aksbot(nj,k)+
     .            akstor(ni,k)*akstor(ni,l)*alsbot(nj,k)*aksbot(nj,l))
     .           +dsqrt(xmub*xmut)*(-2d0*xmub+2d0*x1-2d0)*
     .           (alstor(ni,k)*akstor(ni,l)*alsbot(nj,l)*aksbot(nj,k)+
     .            alstor(ni,l)*akstor(ni,k)*alsbot(nj,k)*aksbot(nj,l))
     .           +dsqrt(xmub)*xmchar(l)/gmst(ni)*
     .           (x1-xmub+xmusb(nj)-xmut-1d0)*
     .           (alstor(ni,k)*akstor(ni,l)*aksbot(nj,k)*aksbot(nj,l)+
     .            alstor(ni,l)*akstor(ni,k)*alsbot(nj,l)*alsbot(nj,k))
     .           +dsqrt(xmub)*xmchar(k)/gmst(ni)*
     .           (x1-xmub+xmusb(nj)-xmut-1d0)*
     .           (alstor(ni,k)*akstor(ni,l)*alsbot(nj,k)*alsbot(nj,l)+
     .            alstor(ni,l)*akstor(ni,k)*aksbot(nj,k)*aksbot(nj,l))
     .           +dsqrt(xmub*xmut)*xmchar(l)/gmst(ni)*
     .           xmchar(k)/gmst(ni)*(-2d0)*
     .           (alstor(ni,k)*akstor(ni,l)*alsbot(nj,k)*aksbot(nj,l)+
     .            alstor(ni,l)*akstor(ni,k)*alsbot(nj,l)*aksbot(nj,k)) )
            endif
         enddo
      enddo

c -------------------------------------------------------------------- c
c                          neutralino exchange
c -------------------------------------------------------------------- c
      neutneut=0d0

      do k=1,5
         do l=1,5
            if((amt+amneut(l)).gt.gmst(ni).and.
     .(amt+amneut(k)).gt.gmst(ni))then
            neutneut=neutneut+3d0*g2s**2/dneut(k)/dneut(l)*(
     .           xmneut(k)*xmneut(l)/gmst(ni)**2*
     .           (-1d0+x1+x2-xmub-xmut+xmusb(nj))*
     .           (atopr(ni,k)*atopr(ni,l)*abot(nj,k)*abot(nj,l)+
     .            btopr(ni,k)*btopr(ni,l)*bbot(nj,k)*bbot(nj,l))
     .           +(xmut*(x2+xmusb(nj)-xmub-x1-xmut)-x2+x1*x2-x1+1d0
     .             -xmusb(nj)+xmub)*
     .           (atopr(ni,k)*atopr(ni,l)*bbot(nj,k)*bbot(nj,l)+
     .            btopr(ni,k)*btopr(ni,l)*abot(nj,k)*abot(nj,l))
     .           +xmneut(k)/gmst(ni)*dsqrt(xmub)*(x2-2d0*xmut)*
     .           (atopr(ni,k)*atopr(ni,l)*abot(nj,k)*bbot(nj,l)+
     .            btopr(ni,k)*btopr(ni,l)*abot(nj,l)*bbot(nj,k))
     .           +xmneut(l)/gmst(ni)*dsqrt(xmub)*(x2-2d0*xmut)*
     .           (atopr(ni,k)*atopr(ni,l)*abot(nj,l)*bbot(nj,k)+
     .            btopr(ni,k)*btopr(ni,l)*abot(nj,k)*bbot(nj,l))
     .           +dsqrt(xmub*xmut)*(-2d0*xmut+2d0*x2-2d0)*
     .           (atopr(ni,k)*btopr(ni,l)*abot(nj,l)*bbot(nj,k)+
     .            atopr(ni,l)*btopr(ni,k)*abot(nj,k)*bbot(nj,l))
     .           +dsqrt(xmut)*xmneut(l)/gmst(ni)*
     .           (x2-xmub+xmusb(nj)-xmut-1d0)*
     .           (atopr(ni,k)*btopr(ni,l)*bbot(nj,k)*bbot(nj,l)+
     .            atopr(ni,l)*btopr(ni,k)*abot(nj,l)*abot(nj,k))
     .           +dsqrt(xmut)*xmneut(k)/gmst(ni)*
     .           (x2-xmub+xmusb(nj)-xmut-1d0)*
     .           (atopr(ni,k)*btopr(ni,l)*abot(nj,k)*abot(nj,l)+
     .            atopr(ni,l)*btopr(ni,k)*bbot(nj,k)*bbot(nj,l))
     .           +dsqrt(xmub*xmut)*xmneut(k)*xmneut(l)/gmst(ni)**2*
     .           (-2d0)*
     .           (atopr(ni,k)*btopr(ni,l)*abot(nj,k)*bbot(nj,l)+
     .            atopr(ni,l)*btopr(ni,k)*abot(nj,l)*bbot(nj,k)) )
            endif
         enddo
      enddo

c -------------------------------------------------------------------- c
c                            gluino exchange
c -------------------------------------------------------------------- c
      gluiglui=0d0
      if ((amt+mgluino).gt.gmst(ni))then
      gluiglui= 2d0/3d0*gs2**2*4d0/dgl**2*(
     .          xmugl*dsqrt(xmub*xmut)*(-4d0)*
     .          gtr(ni)*gtl(ni)*gbr(nj)*gbl(nj)
     .         +dsqrt(xmugl*xmub)*2d0*(x2-2d0*xmut)*
     .          gbr(nj)*gbl(nj)*(gtr(ni)**2+gtl(ni)**2)
     .         +dsqrt(xmub*xmut)*4d0*(x2-xmut-1d0)*
     .          gbr(nj)*gbl(nj)*gtr(ni)*gtl(ni)
     .         +dsqrt(xmut*xmugl)*(-2d0)*(1d0-x2+xmub+xmut-xmusb(nj))*
     .          gtr(ni)*gtl(ni)*(gbr(nj)**2+gbl(nj)**2)
     .         +xmugl*(x1+x2-xmut-xmub-1d0+xmusb(nj))*
     .          (gbr(nj)**2*gtr(ni)**2+gbl(nj)**2*gtl(ni)**2)
     .         +(x2*xmut-xmut*xmub+x1*x2-xmut**2-x1*xmut+xmut*xmusb(nj)
     .           -x1-xmusb(nj)+1d0+xmub-x2)*
     .          (gbr(nj)**2*gtl(ni)**2+gbl(nj)**2*gtr(ni)**2) )
      endif
c -------------------------------------------------------------------- c
c                    chargino neutralino interference
c -------------------------------------------------------------------- c
      charneut=0d0

      do k=1,2
         do l=1,5
      if((amchar(k)+mb).gt.gmst(ni).and.(amt+amneut(l)).gt.gmst(ni))then
            charneut=charneut+2d0*3d0*g2s**2/dchi(k)/dneut(l)*(
     .           xmchar(k)*xmneut(l)/gmst(ni)**2*
     .           (-1d0+x1+x2-xmub-xmut+xmusb(nj))*
     .           (atopr(ni,l)*abot(nj,l)*alstor(ni,k)*alsbot(nj,k)+
     .            btopr(ni,l)*bbot(nj,l)*akstor(ni,k)*aksbot(nj,k))
     .           +(x2-x1*x2+x1-1d0-xmut-xmub-2d0*xmub*xmut+x1*xmut
     .            +x2*xmub+xmusb(nj))*
     .           (atopr(ni,l)*bbot(nj,l)*alsbot(nj,k)*akstor(ni,k)+
     .            btopr(ni,l)*abot(nj,l)*alstor(ni,k)*aksbot(nj,k))
     .           +xmchar(k)/gmst(ni)*dsqrt(xmub)*(x2-2d0*xmut)*
     .           (atopr(ni,l)*bbot(nj,l)*alstor(ni,k)*alsbot(nj,k)+
     .            btopr(ni,l)*abot(nj,l)*akstor(ni,k)*aksbot(nj,k))
     .           +xmneut(l)/gmst(ni)*dsqrt(xmut)*(x1-2d0*xmub)*
     .           (atopr(ni,l)*abot(nj,l)*alstor(ni,k)*aksbot(nj,k)+
     .            btopr(ni,l)*bbot(nj,l)*akstor(ni,k)*alsbot(nj,k))
     .           +dsqrt(xmub*xmut)*(1d0-xmub-xmut+xmusb(nj))*
     .           (atopr(ni,l)*bbot(nj,l)*alstor(ni,k)*aksbot(nj,k)+
     .            btopr(ni,l)*abot(nj,l)*akstor(ni,k)*alsbot(nj,k))
     .           +dsqrt(xmut)*xmchar(k)/gmst(ni)*
     .           (x2-xmub+xmusb(nj)-xmut-1d0)*
     .           (atopr(ni,l)*bbot(nj,l)*akstor(ni,k)*aksbot(nj,k)+
     .            btopr(ni,l)*abot(nj,l)*alstor(ni,k)*alsbot(nj,k))
     .           +dsqrt(xmub)*xmneut(l)/gmst(ni)*
     .           (x1-xmub+xmusb(nj)-xmut-1d0)*
     .           (atopr(ni,l)*abot(nj,l)*akstor(ni,k)*alsbot(nj,k)+
     .            btopr(ni,l)*bbot(nj,l)*alstor(ni,k)*aksbot(nj,k))
     .           +dsqrt(xmub*xmut)*xmchar(k)*xmneut(l)/gmst(ni)**2*
     .           (-2d0)*
     .           (atopr(ni,l)*abot(nj,l)*akstor(ni,k)*aksbot(nj,k)+
     .            btopr(ni,l)*bbot(nj,l)*alstor(ni,k)*alsbot(nj,k)) )
            endif
         enddo
      enddo

      NS_stbsbstart=charchar+neutneut+gluiglui+charneut
      end
c ==================================================================== c
c ======================== sbottom_1/2 bbar top ====================== c
c ==================================================================== c
       DOUBLE PRECISION FUNCTION NS_stbbsbt(x1,x2)
*
       IMPLICIT NONE
       INTEGER ni,nj,l,k,i
       DOUBLE PRECISION gmst(2),gmsb(2),xmusb(2),dw(2),dch(2),gctbr(2,2)
     .,dneut(5),xmuneut(5)
       DOUBLE PRECISION
     .gul(2),gur(2),gdl(2),gdr(2),gtl(2),gtr(2),gbl(2),gbr(2)
       DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
       DOUBLE PRECISION abot(2,5),bbot(2,5)
       DOUBLE PRECISION atopr(2,5),btopr(2,5)
       DOUBLE PRECISION gwtb(2,2),gwntau(2,2),gwnmu(2,2)
       DOUBLE PRECISION scalb,scalt,scaltau,gs2
       DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
       DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

       DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),amch
       DOUBLE PRECISION AMZ,AMW
       DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
       DOUBLE PRECISION chtbrunr,chtbrunl,xmut,xmub,xmuw,
     .xmuch,xmugl,x2,x1,x3,dgl,stsbotww,stsbothh,stsbotgl,stsbotneut,
     .stsbothw,stsbotwneut,stsbothcneut
*
      COMMON/NS_indices/ni,nj
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,amch
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_coup20/gwtb,gwntau,gwnmu
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
*
      COMMON/NS_hcsbotstop/gctbr
      COMMON/NS_higgschudb/chtbrunr,chtbrunl
      COMMON/NS_neutstoptop/atopr,btopr
      COMMON/NS_neutsbotbot/abot,bbot
*
      gmst(1)=ast1
      gmst(2)=ast2
      gmsb(1)=asb1
      gmsb(2)=asb2
      xmusb(1)   = gmsb(1)**2/gmst(ni)**2
      xmusb(2)   = gmsb(2)**2/gmst(ni)**2
      xmut       = amt**2/gmst(ni)**2
      xmub       = amb**2/gmst(ni)**2
      xmuw       = mw**2/gmst(ni)**2
      xmuch      = amch**2/gmst(ni)**2
      xmugl      = mgluino**2/gmst(ni)**2
      do i=1,5
         xmuneut(i) = amneut(i)**2/gmst(ni)**2
      enddo

      x3 = 2d0-x1-x2

      dw(1)   = 1d0-x3+xmusb(1)-xmuw
      dw(2)   = 1d0-x3+xmusb(2)-xmuw
      dch(1)  = 1d0-x3+xmusb(1)-xmuch
      dch(2)  = 1d0-x3+xmusb(2)-xmuch
      dgl     = 1d0-x2+xmut-xmugl
      do i=1,5
         dneut(i) = 1d0-x2+xmut-xmuneut(i)
      enddo
c -------------------------------------------------------------------- c
c                              W exchange
c -------------------------------------------------------------------- c
      stsbotww=0d0
      if (mw+gmsb(nj).gt.gmst(ni))then
      stsbotww=3d0*g2s**2/4d0/dw(nj)**2*gwtb(ni,nj)**2*(
     .     4d0*(1d0+x1*x2-x1-x2-xmusb(nj))
     .   + xmut*(3d0-3d0*x1+x2+xmusb(nj)-xmut+2d0*xmub)
     .   + xmub*(3d0-3d0*x2+x1+xmusb(nj)-xmub)
     .   + xmut/xmuw**2*(1d0-xmusb(nj))**2*
     .     (-1d0+x1+x2+xmusb(nj)-xmut+2d0*xmub)
     .   + xmub/xmuw**2*(1d0-xmusb(nj))**2*
     .     (-1d0+x1+x2+xmusb(nj)-xmub)
     .   + 2d0*xmut/xmuw*(1d0-xmusb(nj))*
     .     (-x1+x2-xmut-1d0+xmusb(nj)+2d0*xmub)
     .   + 2d0*xmub/xmuw*(1d0-xmusb(nj))*
     .     (-x2+x1-xmub-1d0+xmusb(nj)) )
      endif
c -------------------------------------------------------------------- c
c                             H+ exchange
c -------------------------------------------------------------------- c
      stsbothh=0d0
      if (amch+gmsb(nj).gt.gmst(ni))then
      stsbothh=3d0*g2s**2*gctbr(ni,nj)**2*amw**2/gmst(ni)**2
     .     /dch(nj)**2*(
     .     (chtbrunr**2+chtbrunl**2)*(x1+x2+xmusb(nj)-1d0-xmut-xmub)
     .     -4d0*chtbrunr*chtbrunl*dsqrt(xmub*xmut) )
      endif
c -------------------------------------------------------------------- c
c                           gluino exchange
c -------------------------------------------------------------------- c
      stsbotgl=0d0
      if ((amt+mgluino).gt.gmst(ni))then
      stsbotgl=gs2*gs2*2d0/3d0/dgl**2*4d0*(
     .     gtr(ni)*gtl(ni)*gbr(nj)*gbl(nj)*dsqrt(xmub*xmut)*(
     .     -4d0*xmugl-4d0*(1d0+xmut-x2) ) +
     .     gtr(ni)*gtl(ni)*(gbr(nj)**2+gbl(nj)**2)*dsqrt(xmut*xmugl)*
     .     2d0*(-xmut-xmub+xmusb(nj)+x2-1d0) +
     .     gbr(nj)*gbl(nj)*(gtr(ni)**2+gtl(ni)**2)*dsqrt(xmub*xmugl)*
     .     2d0*(x2-2d0*xmut) +
     .     (gtr(ni)**2*gbl(nj)**2+gtl(ni)**2*gbr(nj)**2)*xmugl*
     .     (-1d0-xmut-xmub+xmusb(nj)+x1+x2) +
     .     (gtr(ni)**2*gbr(nj)**2+gtl(ni)**2*gbl(nj)**2)*
     .     (1d0+x1*x2-x1-x2-xmusb(nj)+xmub+xmut*(-xmub+xmusb(nj)+x2
     .      -x1-xmut)) )
      endif
c -------------------------------------------------------------------- c
c                        neutralino exchange
c -------------------------------------------------------------------- c
      stsbotneut = 0d0

      do k=1,5
         do l=1,5
      if((amt+amneut(k)).gt.gmst(ni).and.
     .(amt+amneut(l)).gt.gmst(ni))then
            stsbotneut=stsbotneut+3d0*g2s**2/dneut(k)/dneut(l)*(
     .          (abot(nj,k)*abot(nj,l)*atopr(ni,k)*atopr(ni,l)+
     .           bbot(nj,k)*bbot(nj,l)*btopr(ni,k)*btopr(ni,l))*
     .          ((1d0-x1)*(1d0-x2)-xmusb(nj)+xmut*(-x1+x2+xmusb(nj)
     .           -xmut-xmub)+xmub)
     .         +(abot(nj,k)*abot(nj,l)*btopr(ni,k)*btopr(ni,l)+
     .           bbot(nj,k)*bbot(nj,l)*atopr(ni,k)*atopr(ni,l))*
     .           xmneut(k)*xmneut(l)/gmst(ni)**2*(x1+x2-1d0+xmusb(nj)
     .           -xmut-xmub)
     .         +(abot(nj,k)*bbot(nj,l)*atopr(ni,k)*btopr(ni,l)+
     .           bbot(nj,k)*abot(nj,l)*btopr(ni,k)*atopr(ni,l))*
     .           2d0*dsqrt(xmub*xmut)*(-1d0+x2-xmut)
     .         +dsqrt(xmub)*(x2-2d0*xmut)*( xmneut(k)/gmst(ni)*
     .          (abot(nj,k)*bbot(nj,l)*btopr(ni,k)*btopr(ni,l)+
     .           bbot(nj,k)*abot(nj,l)*atopr(ni,k)*atopr(ni,l))
     .         +xmneut(l)/gmst(ni)*
     .          (abot(nj,k)*bbot(nj,l)*atopr(ni,k)*atopr(ni,l)+
     .           bbot(nj,k)*abot(nj,l)*btopr(ni,k)*btopr(ni,l)) )
     .         +dsqrt(xmut)*(-1d0+x2+xmusb(nj)-xmub-xmut)*
     .         (xmneut(k)/gmst(ni)*
     .          (abot(nj,k)*abot(nj,l)*btopr(ni,k)*atopr(ni,l)+
     .           bbot(nj,k)*bbot(nj,l)*atopr(ni,k)*btopr(ni,l))
     .         +xmneut(l)/gmst(ni)*
     .          (abot(nj,k)*abot(nj,l)*atopr(ni,k)*btopr(ni,l)+
     .           bbot(nj,k)*bbot(nj,l)*btopr(ni,k)*atopr(ni,l)) )
     .         +(abot(nj,l)*bbot(nj,k)*atopr(ni,k)*btopr(ni,l)+
     .           bbot(nj,l)*abot(nj,k)*btopr(ni,k)*atopr(ni,l))*
     .           (-2d0)*dsqrt(xmub*xmut)*
     .             xmneut(k)*xmneut(l)/gmst(ni)**2 )
            endif
         enddo
      enddo

c -------------------------------------------------------------------- c
c                           H+ W interference
c -------------------------------------------------------------------- c
      stsbothw=0d0
      if ((amch+gmsb(nj)).gt.gmst(ni).and.(mw+gmsb(nj)).gt.gmst(ni))
     .then
      stsbothw=-3d0*g2s**2*gctbr(ni,nj)*gwtb(ni,nj)*amw/gmst(ni)
     .     /dw(nj)/dch(nj)*
     .     (dsqrt(xmut)*chtbrunl*(xmusb(nj)+xmub-xmut+x2-x1-1d0)
     .     +dsqrt(xmub)*chtbrunr*(-xmusb(nj)+xmub-xmut+x2-x1+1d0)
     .     +dsqrt(xmut)/xmuw*chtbrunl*(xmusb(nj)-1d0)*(xmut-xmub-x1-x2
     .      +1d0-xmusb(nj))
     .     +dsqrt(xmub)/xmuw*chtbrunr*(xmusb(nj)-1d0)*(xmut-xmub+x1+x2
     .      -1d0+xmusb(nj)) )
      endif
c -------------------------------------------------------------------- c
c                       neutralino W interference
c -------------------------------------------------------------------- c
      stsbotwneut = 0d0

      do l=1,5
          if((mw+gmsb(nj)).gt.gmst(ni).and.(amt+amneut(l)).gt.gmst(ni))
     .then
         stsbotwneut=stsbotwneut
     .     +3d0*g2s**2/dw(nj)/dneut(l)*gwtb(ni,nj)*(
     .      bbot(nj,l)*btopr(ni,l)*dsqrt(xmut*xmub)*(1d0/xmuw*
     .      (xmusb(nj)*(xmusb(nj)+xmut-xmub-2d0)-xmut+xmub+1d0) +
     .      2d0*x2-xmusb(nj)-xmut+xmub-3d0) +
     .      abot(nj,l)*atopr(ni,l)*(1d0/xmuw*(1d0-xmusb(nj))*(xmut*
     .      (xmusb(nj)+x2-xmut+xmub-1d0)-x2*xmub) + xmut*(1d0+x2+
     .      xmusb(nj)+xmub-2d0*x1-xmut)-2d0*xmusb(nj)+2d0*x1*x2+2d0
     .      -2d0*x1-2d0*x2+xmub*(2d0-x2)) +
     .      bbot(nj,l)*atopr(ni,l)*dsqrt(xmub)*xmneut(l)/gmst(ni)*
     .      (1d0/xmuw*
     .      (xmusb(nj)-1d0)*(xmusb(nj)+x1+x2-1d0-xmub+xmut)+x2-x1
     .      -xmusb(nj)-xmut+xmub+1d0) +
     .      abot(nj,l)*btopr(ni,l)*dsqrt(xmut)*xmneut(l)/gmst(ni)*
     .      (1d0/xmuw*
     .      (1d0-xmusb(nj))*(xmusb(nj)+x1+x2-1d0+xmub-xmut)+x2-x1
     .      +xmusb(nj)-xmut+xmub-1d0) )
         endif
      enddo
c -------------------------------------------------------------------- c
c                    neutralino Higgs interference
c -------------------------------------------------------------------- c
      stsbothcneut = 0d0

      do l=1,5
      if ((amt+amneut(l)).gt.gmst(ni).and.(gmsb(nj)+amch).gt.gmst(ni))
     .then
         stsbothcneut=stsbothcneut
     .     +3d0*g2s**2*2d0/dch(nj)/dneut(l)*
     .      (-gctbr(ni,nj))*amw/gmst(ni)*(
     .      (chtbrunr*btopr(ni,l)*abot(nj,l)+chtbrunl*atopr(ni,l)*
     .       bbot(nj,l))*dsqrt(xmut*xmub)*xmneut(l)/gmst(ni)*(-2d0) +
     .      (chtbrunr*btopr(ni,l)*bbot(nj,l)+chtbrunl*atopr(ni,l)*
     .       abot(nj,l))*dsqrt(xmut)*(-1d0-xmut-xmub+xmusb(nj)+x2) +
     .      (chtbrunr*atopr(ni,l)*bbot(nj,l)+chtbrunl*btopr(ni,l)*
     .       abot(nj,l))*xmneut(l)/gmst(ni)*(-1d0-xmut-xmub+xmusb(nj)
     .       +x1+x2) +
     .      (chtbrunr*atopr(ni,l)*abot(nj,l)+chtbrunl*btopr(ni,l)*
     .       bbot(nj,l))*dsqrt(xmub)*(-2d0*xmut+x2) )
         endif
      enddo
      NS_stbbsbt = stsbotww+stsbothh+stsbothw+stsbotgl+stsbotneut+
     .             stsbotwneut+stsbothcneut
      end
c ==================================================================== c
c ====================== sbottom_1/2 tau+ nu_tau ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_sttausbnu(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj
      DOUBLE PRECISION gmst(2),gmsb(2),xmusb(2),dw(2),dch(2),gctbr(2,2)
      DOUBLE PRECISION gwtb(2,2),gwntau(2,2),gwnmu(2,2)
      DOUBLE PRECISION
     .gul(2),gur(2),gdl(2),gdr(2),gtl(2),gtr(2),gbl(2),gbr(2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),amch
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION chctaunur,chctaunul
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION xmutau,xmuw,xmuch,x1,x2,x3,ststauww,ststauhh
     .,ststauhw
*
      COMMON/NS_indices/ni,nj
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,amch
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_coup14/chctaunur,chctaunul
      COMMON/NS_coup20/gwtb,gwntau,gwnmu
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_hcsbotstop/gctbr
*
      gmst(1)=ast1
      gmst(2)=ast2
      gmsb(1)=asb1
      gmsb(2)=asb2
      xmusb(1) = gmsb(1)**2/gmst(ni)**2
      xmusb(2) = gmsb(2)**2/gmst(ni)**2
      xmutau   = amtau**2/gmst(ni)**2
      xmuw     = mw**2/gmst(ni)**2
      xmuch    = amch**2/gmst(ni)**2

      x3 = 2d0-x1-x2

      dw(1)  = 1d0-x3+xmusb(1)-xmuw
      dw(2)  = 1d0-x3+xmusb(2)-xmuw
      dch(1) = 1d0-x3+xmusb(1)-xmuch
      dch(2) = 1d0-x3+xmusb(2)-xmuch
c -------------------------------------------------------------------- c
c                             W exchange
c -------------------------------------------------------------------- c
      ststauww=0d0
      if((gmsb(nj)+mw).gt.gmst(ni))then
      ststauww=1d0/4d0/dw(nj)**2*gwtb(ni,nj)**2*(
     .     4d0*(1d0+x1*x2-x1-x2-xmusb(nj)) +
     .     xmutau*(3d0-3d0*x1+x2+xmusb(nj)-xmutau)
     .     +xmutau/xmuw**2*(1d0-xmusb(nj))**2*(-1d0+xmusb(nj)+x1+x2
     .      -xmutau)
     .     +2d0*xmutau/xmuw*(1d0-xmusb(nj))*(-x1+x2-xmutau-1d0+
     .      xmusb(nj)) )
      endif
c -------------------------------------------------------------------- c
c                             H+ exchange
c -------------------------------------------------------------------- c
      ststauhh=0d0
      if((gmsb(nj)+amch).gt.gmst(ni))then
      ststauhh=gctbr(ni,nj)**2*amw**2/gmst(ni)**2/dch(nj)**2*
     .     chctaunur**2*(x1+x2+xmusb(nj)-1d0-xmutau)
      endif
c -------------------------------------------------------------------- c
c                          H+ W interference
c -------------------------------------------------------------------- c
      ststauhw=0d0
      if((gmsb(nj)+amch).gt.gmst(ni).and.(gmsb(nj)+mw).gt.gmst(ni))then
      ststauhw=-gctbr(ni,nj)*gwtb(ni,nj)*amw/gmst(ni)/dw(nj)/dch(nj)*
     .      dsqrt(xmutau)*chctaunur*(
     .      1d0+x1-x2+xmutau-xmusb(nj)+
     .      1d0/xmuw*(xmusb(nj)-1d0)*(-xmutau+x1+x2-1d0+xmusb(nj)))
      endif
      NS_sttausbnu = ststauww+ststauhh+ststauhw

      end
c ==================================================================== c
c ======================== sbottom_1/2 e+ nu_e ======================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_stelsbnu(x1,x2)
*
      IMPLICIT NONE
      INTEGER ni,nj
      DOUBLE PRECISION gmst(2),gmsb(2),xmusb(2),dw(2)
      DOUBLE PRECISION gwtb(2,2),gwntau(2,2),gwnmu(2,2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION XMUW,X1,X2,X3,stselww
*
      COMMON/NS_indices/ni,nj
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_coup20/gwtb,gwntau,gwnmu
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
*
      gmst(1)=ast1
      gmst(2)=ast2
      gmsb(1)=asb1
      gmsb(2)=asb2

      xmusb(1) = gmsb(1)**2/gmst(ni)**2
      xmusb(2) = gmsb(2)**2/gmst(ni)**2
      xmuw     = mw**2/gmst(ni)**2

      x3 = 2d0-x1-x2

      dw(1) = 1d0-x3+xmusb(1)-xmuw
      dw(2) = 1d0-x3+xmusb(2)-xmuw
c -------------------------------------------------------------------- c
c                             W exchange
c -------------------------------------------------------------------- c
      stselww=0d0
      if ((gmsb(nj)+mw).gt.gmst(ni))then
      stselww=1d0/4d0/dw(nj)**2*gwtb(ni,nj)**2*
     .        4d0*(1d0+x1*x2-x1-x2-xmusb(nj))
      endif
c -------------------------------------------------------------------- c
      NS_stelsbnu = stselww

      end
c ==================================================================== c
c =========================== stop1* top top ========================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_st2st1startt(x1,x2)
*
      IMPLICIT NONE
      INTEGER l,k,i
      DOUBLE PRECISION atopr(2,5),btopr(2,5)
      DOUBLE PRECISION dneut(5),xmuneut(5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION
     .gul(2),gur(2),gdl(2),gdr(2),gtl(2),gtr(2),gbl(2),gbr(2)
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION xmut,xmugl,xmust1,x1,x2,dgl,st2st1neut
     .,st2st1gg
*
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_neutstoptop/atopr,btopr

      xmut   = amt**2/ast2**2
      do i=1,5
         xmuneut(i) = amneut(i)**2/ast2**2
      enddo
      xmugl  = mgluino**2/ast2**2
      xmust1 = ast1**2/ast2**2

      do i=1,5
         dneut(i) = 1d0-x1+xmut-xmuneut(i)
      enddo
      dgl = 1d0-x1+xmut-xmugl
c -------------------------------------------------------------------- c
c                          neutralino exchange
c -------------------------------------------------------------------- c
      st2st1neut=0d0

      do k=1,5
         do l=1,5
       if((amt+amneut(k)).gt.ast2.and.(amt+amneut(l)).gt.ast2)
     .then
            st2st1neut=st2st1neut+1d0/dneut(k)/dneut(l)*(
     .          (btopr(1,k)*btopr(1,l)*atopr(2,k)*atopr(2,l)+
     .           atopr(1,k)*atopr(1,l)*btopr(2,k)*btopr(2,l))*
     .          ((1d0-x1)*(1d0-x2)-xmust1+xmut*(x1-x2+xmust1
     .           -2d0*xmut)+xmut)
     .         +(btopr(1,k)*btopr(1,l)*btopr(2,k)*btopr(2,l)+
     .           atopr(1,k)*atopr(1,l)*atopr(2,k)*atopr(2,l))*
     .           xmneut(k)*xmneut(l)/ast2**2*(x1+x2-1d0+xmust1
     .           -2d0*xmut)
     .         +(btopr(1,k)*atopr(1,l)*atopr(2,k)*btopr(2,l)+
     .           atopr(1,k)*btopr(1,l)*btopr(2,k)*atopr(2,l))*
     .           2d0*xmut*(-1d0+x1-xmut)
     .         +dsqrt(xmut)*(x1-2d0*xmut)*(xmneut(k)/ast2*
     .          (btopr(1,k)*atopr(1,l)*btopr(2,k)*btopr(2,l)+
     .           atopr(1,k)*btopr(1,l)*atopr(2,k)*atopr(2,l))
     .         +xmneut(l)/ast2*
     .          (btopr(1,k)*atopr(1,l)*atopr(2,k)*atopr(2,l)+
     .           atopr(1,k)*btopr(1,l)*btopr(2,k)*btopr(2,l)) )
     .         +dsqrt(xmut)*(-1d0+x1+xmust1-2d0*xmut)*
     .         (xmneut(k)/ast2*
     .          (btopr(1,k)*btopr(1,l)*btopr(2,k)*atopr(2,l)+
     .           atopr(1,k)*atopr(1,l)*atopr(2,k)*btopr(2,l))
     .         +xmneut(l)/ast2*
     .          (btopr(1,k)*btopr(1,l)*atopr(2,k)*btopr(2,l)+
     .           atopr(1,k)*atopr(1,l)*btopr(2,k)*atopr(2,l)) )
     .         +(btopr(1,l)*atopr(1,k)*atopr(2,k)*btopr(2,l)+
     .           atopr(1,l)*btopr(1,k)*btopr(2,k)*atopr(2,l))*
     .           (-2d0)*xmut*xmneut(k)*xmneut(l)/ast2**2 )
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c                             gluino exchange
c -------------------------------------------------------------------- c
      st2st1gg=0d0
      if((amt+mgluino).gt.ast2)then
      st2st1gg=1d0/dgl**2*4d0*( -4d0*dsqrt(xmut*xmugl)*xmut*
     .     (gtr(2)*gtr(1)+gtl(1)*gtl(2))*(gtl(1)*gtr(2)+gtr(1)*gtl(2))
     .    +2d0*dsqrt(xmut*xmugl)*(
     .     (gtl(1)*gtr(1)*(gtr(2)**2+gtl(2)**2)+gtr(2)*gtl(2)*
     .     (gtl(1)**2+gtr(1)**2))*x1 +
     .     gtr(2)*gtl(2)*(gtl(1)**2+gtr(1)**2)*(xmust1-1d0))
     .    +(-2d0)*xmut*xmugl*(gtl(1)*gtl(2)+gtr(2)*gtr(1))**2
     .    +xmut*((gtl(1)**2*gtr(2)**2+gtr(1)**2*gtl(2)**2)*(1d0+x1
     .     -x2+xmust1)+4d0*gtl(1)*gtr(2)*gtr(1)*gtl(2)*(x1-1d0))
     .    +xmut**2*(-2d0)*(gtl(1)*gtr(2)+gtr(1)*gtl(2))**2
     .    +xmugl*(gtl(1)**2*gtl(2)**2+gtr(2)**2*gtr(1)**2)*(x1+x2+
     .     xmust1-1d0)
     .    +(gtl(1)**2*gtr(2)**2+gtr(1)**2*gtl(2)**2)*(x1*x2-x1-x2
     .     -xmust1+1d0) )
      endif
c -------------------------------------------------------------------- c
      NS_st2st1startt=3d0*g2s**2*st2st1neut+
     .                2d0/3d0*gs2**2*st2st1gg


      end
c ==================================================================== c
c ========================== stop1 top topbar ======================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_st2st1tt(x1,x2)
*
      IMPLICIT NONE
      INTEGER I,J,K,L
      DOUBLE PRECISION atopr(2,5),btopr(2,5)
      DOUBLE PRECISION Hstopstopr(3,2,2),Astopstopr(3,2,2)
      DOUBLE PRECISION dneut(5),xmuneut(5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2),gzmumu(2,2)
      DOUBLE PRECISION
     .gul(2),gur(2),gdl(2),gdr(2),gtl(2),gtr(2),gbl(2),gbr(2)
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION Httr(3),Attr(2)
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION P(2,3)
      DOUBLE PRECISION xmuh(3),xmua(2),dh(3),da(2)
      DOUBLE PRECISION aztt11,aztt12,aztt21,aztt22,xmut,xmuz,
     .xmugl,xmust1,x3,x1,x2,dgl,dz,vzz,azz,st2st1zz,st2st1hk,
     .st2st1aa,st2st1higgs,st2st1neut,st2st1gg,st2st1neutz,
     .st2st1hneut,st2st1hz
*
      COMMON/NS_HIGGSSTST/Hstopstopr,Astopstopr
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau,gzmumu
      COMMON/NS_coup21/gtr,gtl,gbr,gbl,gur,gul,gdr,gdl
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_CPodd_MIX/P
      COMMON/NS_MZMWscaleQ/AMZ,AMW
*
      COMMON/NS_neutstoptop/atopr,btopr
      COMMON/NS_phitoptop/Httr,Attr

      aztt11=gztt(1,1)
      aztt12=gztt(1,2)
      aztt21=gztt(2,1)
      aztt22=gztt(2,2)

c --- several definitions ---
      xmut   = amt**2/ast2**2
      do i=1,5
         xmuneut(i) = amneut(i)**2/ast2**2
      enddo
      xmugl  = mgluino**2/ast2**2
      xmuz   = mz**2/ast2**2
      xmust1 = ast1**2/ast2**2
      do i =1,3
      xmuh(i) = SMASS(i)**2/ast2**2
      enddo
      do i =1,2
      xmua(i) = PMASS(i)**2/ast2**2
      enddo

      x3  = 2d0-x1-x2

      do i=1,5
         dneut(i) = 1d0-x1+xmut-xmuneut(i)
      enddo
      dgl = 1d0-x1+xmut-xmugl
      dz  = 1d0-x3+xmust1-xmuz

      do i =1,3
      dh(i) = 1d0-x3+xmust1-xmuh(i)
      enddo
      do i =1,2
      da(i) = 1d0-x3+xmust1-xmua(i)
      enddo
      vzz = vzztoptop
      azz = azztoptop
c -------------------------------------------------------------------- c
c                              Z exchange
c -------------------------------------------------------------------- c
      st2st1zz=0d0
      if ((ast1+mz).gt.ast2)then
      st2st1zz = 1d0/4d0/cw**2/dz**2*aztt12**2*(
     .      8d0*xmut*azz**2*(-x1-x2+xmust1+3d0)
     .     -8d0*(azz**2+vzz**2)*(x1+x2-1d0+xmust1-x1*x2)
     .     -16d0*azz**2*(1d0-xmust1)**2*xmut/xmuz
     .     +8d0*azz**2*xmut/xmuz**2*(1d0-xmust1)**2*
     .     (x1+x2+xmust1-1d0))
      endif
c -------------------------------------------------------------------- c
c                              Higgs exchange
c -------------------------------------------------------------------- c
      st2st1hk=0.0d0
      st2st1aa=0.0d0

             do I = 1,3
                do J = 1,3
      if ((smass(i)+ast1).gt.ast2.and.(smass(j)+ast1).gt.ast2)then
           st2st1hk = st2st1hk+(2d0/ast2**2/dh(I)/dh(J)*
     .(Httr(I)/dsqrt(2d0))*(Httr(J)/dsqrt(2d0))*
     .    amz**4/amw**2*Hstopstopr(I,1,2)*Hstopstopr(J,1,2))*
     .    (-1d0+x1+x2+xmust1-4d0*xmut)
      endif
                enddo
             enddo


             do I = 1,2
                do J = 1,2
       if ((pmass(i)+ast1).gt.ast2.and.(pmass(j)+ast1).gt.ast2)then
            st2st1aa = st2st1aa+2d0/ast2**2/da(I)/da(J)*
     .(Attr(I)/dsqrt(2d0))*(Attr(J)/dsqrt(2d0))*
     .     amz**4/amw**2*Astopstopr(I,1,2)*Astopstopr(J,1,2)
     .*(-1d0+x1+x2+xmust1)
                   endif
                enddo
             enddo
      st2st1higgs = st2st1hk+st2st1aa
c -------------------------------------------------------------------- c
c                          neutralino exchange
c -------------------------------------------------------------------- c
      st2st1neut=0d0

      do k=1,5
         do l=1,5
            if ((amt+amneut(k)).gt.ast2
     ..and.(amt+amneut(l)).gt.ast2)then
            st2st1neut=st2st1neut+1d0/dneut(k)/dneut(l)*(
     .          (atopr(1,k)*atopr(1,l)*atopr(2,k)*atopr(2,l)+
     .           btopr(1,k)*btopr(1,l)*btopr(2,k)*btopr(2,l))*
     .          ((1d0-x1)*(1d0-x2)-xmust1+xmut*(x1-x2+xmust1
     .           -2d0*xmut)+xmut)
     .         +(atopr(1,k)*atopr(1,l)*btopr(2,k)*btopr(2,l)+
     .           btopr(1,k)*btopr(1,l)*atopr(2,k)*atopr(2,l))*
     .           xmneut(k)*xmneut(l)/ast2**2*(x1+x2-1d0+xmust1
     .           -2d0*xmut)
     .         +(atopr(1,k)*btopr(1,l)*atopr(2,k)*btopr(2,l)+
     .           btopr(1,k)*atopr(1,l)*btopr(2,k)*atopr(2,l))*
     .           2d0*xmut*(-1d0+x1-xmut)
     .         +dsqrt(xmut)*(x1-2d0*xmut)*(xmneut(k)/ast2*
     .          (atopr(1,k)*btopr(1,l)*btopr(2,k)*btopr(2,l)+
     .           btopr(1,k)*atopr(1,l)*atopr(2,k)*atopr(2,l))
     .         +xmneut(l)/ast2*
     .          (atopr(1,k)*btopr(1,l)*atopr(2,k)*atopr(2,l)+
     .           btopr(1,k)*atopr(1,l)*btopr(2,k)*btopr(2,l)) )
     .         +dsqrt(xmut)*(-1d0+x1+xmust1-2d0*xmut)*
     .         (xmneut(k)/ast2*
     .          (atopr(1,k)*atopr(1,l)*btopr(2,k)*atopr(2,l)+
     .           btopr(1,k)*btopr(1,l)*atopr(2,k)*btopr(2,l))
     .         +xmneut(l)/ast2*
     .          (atopr(1,k)*atopr(1,l)*atopr(2,k)*btopr(2,l)+
     .           btopr(1,k)*btopr(1,l)*btopr(2,k)*atopr(2,l)) )
     .         +(atopr(1,l)*btopr(1,k)*atopr(2,k)*btopr(2,l)+
     .           btopr(1,l)*atopr(1,k)*btopr(2,k)*atopr(2,l))*
     .           (-2d0)*xmut*xmneut(k)*xmneut(l)/ast2**2 )
            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c                             gluino exchange
c -------------------------------------------------------------------- c
      st2st1gg=0d0
      if ((amt+mgluino).gt.ast2)then
      st2st1gg=1d0/dgl**2*4d0*( -4d0*dsqrt(xmut*xmugl)*xmut*
     .     (gtr(2)*gtl(1)+gtr(1)*gtl(2))*(gtr(1)*gtr(2)+gtl(1)*gtl(2))
     .    +2d0*dsqrt(xmut*xmugl)*(
     .     (gtr(1)*gtl(1)*(gtr(2)**2+gtl(2)**2)+gtr(2)*gtl(2)*
     .     (gtr(1)**2+gtl(1)**2))*x1 +
     .     gtr(2)*gtl(2)*(gtr(1)**2+gtl(1)**2)*(xmust1-1d0))
     .    +(-2d0)*xmut*xmugl*(gtr(1)*gtl(2)+gtr(2)*gtl(1))**2
     .    +xmut*((gtr(1)**2*gtr(2)**2+gtl(1)**2*gtl(2)**2)*(1d0+x1
     .     -x2+xmust1)+4d0*gtr(1)*gtr(2)*gtl(1)*gtl(2)*(x1-1d0))
     .    +xmut**2*(-2d0)*(gtr(1)*gtr(2)+gtl(1)*gtl(2))**2
     .    +xmugl*(gtr(1)**2*gtl(2)**2+gtr(2)**2*gtl(1)**2)*(x1+x2+
     .     xmust1-1d0)
     .    +(gtr(1)**2*gtr(2)**2+gtl(1)**2*gtl(2)**2)*(x1*x2-x1-x2
     .     -xmust1+1d0) )
      endif
c -------------------------------------------------------------------- c
c                         neutralino Z interference
c -------------------------------------------------------------------- c
      st2st1neutz=0d0

      do l=1,5
        if ((amt+amneut(l)).gt.ast2.and.(mz+ast1).gt.ast2)then
         st2st1neutz=st2st1neutz+aztt12/cw/dz/dneut(l)*(
     .      xmut*(atopr(1,l)*atopr(2,l)*(vzz-azz)+
     .            btopr(1,l)*btopr(2,l)*(vzz+azz))*(
     .      1d0/xmuz*(xmust1*(xmust1-2d0)+1d0)+
     .      2d0*x1-xmust1-3d0)
     .     +(atopr(1,l)*atopr(2,l)*(vzz+azz)+
     .       btopr(1,l)*btopr(2,l)*(vzz-azz))*(
     .      1d0/xmuz*(1d0-xmust1)*(xmut*(xmust1+x1-1d0)-x1*xmut)
     .      +xmut*(1d0+x1+xmust1-2d0*x2)-2d0*xmust1+2d0*x1*x2+2d0
     .      -2d0*(x1+x2)+xmut*(-x1+2d0) )
     .     +dsqrt(xmut)*xmneut(l)/ast2*
     .      (atopr(1,l)*btopr(2,l)*(vzz-azz)+
     .       btopr(1,l)*atopr(2,l)*(vzz+azz))*(
     .      1d0/xmuz*(xmust1-1d0)*(xmust1+x1+x2-1d0)+x1-x2-xmust1
     .      +1d0 )
     .     +dsqrt(xmut)*xmneut(l)/ast2*
     .      (atopr(1,l)*btopr(2,l)*(vzz+azz)+
     .       btopr(1,l)*atopr(2,l)*(vzz-azz))*(
     .      1d0/xmuz*(1d0-xmust1)*(xmust1+x1+x2-1d0)+x1-x2+xmust1
     .      -1d0 ) )
         endif
      enddo
c -------------------------------------------------------------------- c
c                       neutralino Higgs interference
c -------------------------------------------------------------------- c
      st2st1hneut=0d0

       do l=1,5
           do I = 1,3
              if ((amt+amneut(l)).gt.ast2.and.
     .(smass(i)+ast1).gt.ast2)then
         st2st1hneut=st2st1hneut-2d0*(Httr(I)/dsqrt(2d0)/dneut(l)
     .        /dh(I)/ast2*amz**2/amw*Hstopstopr(I,1,2))
     .        *((atopr(1,l)*btopr(2,l)+atopr(2,l)*btopr(1,l))*
     .        xmneut(l)/ast2*(x1+x2+xmust1-1d0-4d0*xmut)
     .        +(atopr(1,l)*atopr(2,l)+btopr(2,l)*btopr(1,l))*
     .        dsqrt(xmut)*(2d0*x1-1d0+xmust1-4d0*xmut) )
              endif
          enddo
          enddo
          do l=1,5
           do I = 1,2
        if ((amt+amneut(l)).gt.ast2.and.
     .(pmass(i)+ast1).gt.ast2)then
              st2st1hneut=st2st1hneut+
     .        2d0*Attr(I)/dsqrt(2d0)/dneut(l)/da(I)/ast2*amz**2/amw*
     .        (-Astopstopr(I,1,2))*(
     .        (atopr(1,l)*atopr(2,l)-btopr(2,l)*btopr(1,l))*
     .        dsqrt(xmut)*(1d0-xmust1) +
     .        (atopr(1,l)*btopr(2,l)-atopr(2,l)*btopr(1,l))*
     .        xmneut(l)/ast2*(1d0-x1-x2-xmust1) )
              endif
              enddo
          enddo
c -------------------------------------------------------------------- c
c                          Higgs Z interference
c -------------------------------------------------------------------- c
          st2st1hz=0.0d0
          do I =1,3
        if((smass(i)+ast1).gt.ast2.and.(mz+ast1).gt.ast2)then
      st2st1hz=st2st1hz-2d0/2d0/cw*
     .   (aztt12*Httr(I)/dsqrt(2d0)*amz**2/amw*Hstopstopr(I,1,2)
     .   /ast2/dz/dh(I))*2d0*dsqrt(xmut)*vzz*2d0*(x1-x2)
             endif
          enddo
          do I =1,2
        if((pmass(i)+ast1).gt.ast2.and.(mz+ast1).gt.ast2)then
      st2st1hz=st2st1hz+2d0/2d0/cw*
     .    aztt12*Attr(I)/dsqrt(2d0)*amz**2/amw*
     .(-Astopstopr(I,1,2))/ast2/dz/da(I)*
     .   (2d0*dsqrt(xmut)*azz*(2d0/xmuz*(1d0+(xmust1-1d0)*(x1
     .   +x2)+xmust1**2-2d0*xmust1)+2d0-2d0*xmust1) )
            endif
         enddo
c -------------------------------------------------------------------- c
      NS_st2st1tt=3d0*g2s**2*(st2st1zz+st2st1higgs+st2st1neut+
     .            st2st1neutz+st2st1hneut+st2st1hz) +
     .            gs2**2*2d0/3d0*st2st1gg

      end
c ==================================================================== c
c ======================= stop1 bottom bottombar ===================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_st2st1bb(x1,x2)
*
      IMPLICIT NONE
      INTEGER i,j,k,l
      DOUBLE PRECISION dchi(2),xmuchar(2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION alstor(2,2),akstor(2,2),alstor1(2,2),
     .        blstor1(2,2),alstor2(2,2),blstor2(2,2)
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2),gzmumu(2,2)
      DOUBLE PRECISION Hstopstopr(3,2,2),Astopstopr(3,2,2)
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION Hbbr(3),Abbr(2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION P(2,3)
      DOUBLE PRECISION xmuh(3),xmua(2),dh(3),da(2)
      DOUBLE PRECISION aztt11,aztt12,aztt21,aztt22,xmub,xmuz,
     .xmust1,x3,x1,x2,dz,vzz,azz,st2st1zz,st2st1hk,
     .st2st1aa,st2st1higgs,st2st1chi,st2st1chiz,
     .st2st1hz,st2st1hchi
*
      COMMON/NS_HIGGSSTST/Hstopstopr,Astopstopr
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau,gzmumu
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_CPodd_MIX/P
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_charstopbot/alstor,akstor
      COMMON/NS_phibotbot/Hbbr,Abbr

      aztt11=gztt(1,1)
      aztt12=gztt(1,2)
      aztt21=gztt(2,1)
      aztt22=gztt(2,2)
      do i=1,2
         do j=1,2
            alstor1(i,j) = alstor(i,j)
            blstor1(i,j) = akstor(i,j)
            alstor2(i,j) = alstor(i,j)
            blstor2(i,j) = akstor(i,j)
         enddo
      enddo

c --- several definitions ---
      xmub       = amb**2/ast2**2
      xmuchar(1) = amchar(1)**2/ast2**2
      xmuchar(2) = amchar(2)**2/ast2**2
      xmuz       = mz**2/ast2**2
      xmust1     = ast1**2/ast2**2

      do i =1,3
      xmuh(i) = SMASS(i)**2/ast2**2
      enddo
      do i =1,2
      xmua(i) = PMASS(i)**2/ast2**2
      enddo

      x3  = 2d0-x1-x2

      dchi(1) = 1d0-x1+xmub-xmuchar(1)
      dchi(2) = 1d0-x1+xmub-xmuchar(2)
      dz      = 1d0-x3+xmust1-xmuz

      do i =1,3
      dh(i) = 1d0-x3+xmust1-xmuh(i)
      enddo
      do i =1,2
      da(i) = 1d0-x3+xmust1-xmua(i)
      enddo

      vzz = vzzbotbot
      azz = azzbotbot
c -------------------------------------------------------------------- c
c                               Z exchange
c -------------------------------------------------------------------- c
      st2st1zz=0d0
       if ((ast1+mz).gt.ast2)then
      st2st1zz = 1d0/4d0/cw**2/dz**2*aztt12**2*(
     .      8d0*xmub*azz**2*(-x1-x2+xmust1+3d0)
     .     -8d0*(azz**2+vzz**2)*(x1+x2-1d0+xmust1-x1*x2)
     .     -16d0*azz**2*(1d0-xmust1)**2*xmub/xmuz
     .     +8d0*azz**2*xmub/xmuz**2*(1d0-xmust1)**2*
     .     (x1+x2+xmust1-1d0))
      endif
c -------------------------------------------------------------------- c
c                              Higgs exchange
c -------------------------------------------------------------------- c
      st2st1hk=0.0d0
      st2st1aa=0.0d0
                do I = 1,3
                   do J = 1,3
      if ((smass(i)+ast1).gt.ast2.and.(smass(j)+ast1).gt.ast2)then
           st2st1hk = st2st1hk+(2d0/ast2**2/dh(I)/dh(J)*
     .(Hbbr(I)/dsqrt(2d0))*(Hbbr(J)/dsqrt(2d0))*
     .    amz**4/amw**2*Hstopstopr(I,1,2)*Hstopstopr(J,1,2))*
     .    (-1d0+x1+x2+xmust1-4d0*xmub)
           endif
                   enddo
                enddo

             do I = 1,2
                do J = 1,2
       if ((pmass(i)+ast1).gt.ast2.and.(pmass(j)+ast1).gt.ast2)then
            st2st1aa = st2st1aa+2d0/ast2**2/da(I)/da(J)*
     .(Abbr(I)/dsqrt(2d0))*(Abbr(J)/dsqrt(2d0))*
     .     amz**4/amw**2*Astopstopr(I,1,2)*Astopstopr(J,1,2)
     .*(-1d0+x1+x2+xmust1)
           endif
                enddo
             enddo
              st2st1higgs = st2st1hk+st2st1aa
c -------------------------------------------------------------------- c
c                          chargino exchange
c -------------------------------------------------------------------- c
      st2st1chi=0d0

      do k=1,2
         do l=1,2
      If((amchar(k)+mb).gt.ast2.and.(amchar(l)+mb).gt.ast2)
     .then
            st2st1chi=st2st1chi+1d0/dchi(k)/dchi(l)*(
     .          (alstor2(1,k)*alstor2(1,l)*alstor1(2,k)*alstor1(2,l)+
     .           blstor2(1,k)*blstor2(1,l)*blstor1(2,k)*blstor1(2,l))*
     .          ((1d0-x1)*(1d0-x2)-xmust1+xmub*(x1-x2+xmust1
     .           -2d0*xmub)+xmub)
     .         +(alstor2(1,k)*alstor2(1,l)*blstor1(2,k)*blstor1(2,l)+
     .           blstor2(1,k)*blstor2(1,l)*alstor1(2,k)*alstor1(2,l))*
     .           xmchar(k)*xmchar(l)/ast2**2*(x1+x2-1d0+xmust1
     .           -2d0*xmub)
     .         +(alstor2(1,k)*blstor2(1,l)*alstor1(2,k)*blstor1(2,l)+
     .           blstor2(1,k)*alstor2(1,l)*blstor1(2,k)*alstor1(2,l))*
     .           2d0*xmub*(-1d0+x1-xmub)
     .         +dsqrt(xmub)*(x1-2d0*xmub)*(xmchar(k)/ast2*
     .          (alstor2(1,k)*blstor2(1,l)*blstor1(2,k)*blstor1(2,l)+
     .           blstor2(1,k)*alstor2(1,l)*alstor1(2,k)*alstor1(2,l))
     .         +xmchar(l)/ast2*
     .          (alstor2(1,k)*blstor2(1,l)*alstor1(2,k)*alstor1(2,l)+
     .           blstor2(1,k)*alstor2(1,l)*blstor1(2,k)*blstor1(2,l)) )
     .         +dsqrt(xmub)*(-1d0+x1+xmust1-2d0*xmub)*
     .         (xmchar(k)/ast2*
     .          (alstor2(1,k)*alstor2(1,l)*blstor1(2,k)*alstor1(2,l)+
     .           blstor2(1,k)*blstor2(1,l)*alstor1(2,k)*blstor1(2,l))
     .         +xmchar(l)/ast2*
     .          (alstor2(1,k)*alstor2(1,l)*alstor1(2,k)*blstor1(2,l)+
     .           blstor2(1,k)*blstor2(1,l)*blstor1(2,k)*alstor1(2,l)) )
     .         +(alstor2(1,l)*blstor2(1,k)*alstor1(2,k)*blstor1(2,l)+
     .           blstor2(1,l)*alstor2(1,k)*blstor1(2,k)*alstor1(2,l))*
     .           (-2d0)*xmub*xmchar(k)/ast2*xmchar(l)/ast2 )

            endif
         enddo
      enddo
c -------------------------------------------------------------------- c
c                         chargino Z interference
c -------------------------------------------------------------------- c
      st2st1chiz=0d0

      do l=1,2

        if((amchar(l)+mb).gt.ast2.and.(ast1+mz).gt.ast2)then
         st2st1chiz=st2st1chiz+aztt12/cw/dz/dchi(l)*(
     .      xmub*(alstor2(1,l)*alstor1(2,l)*(vzz-azz)+
     .            blstor2(1,l)*blstor1(2,l)*(vzz+azz))*(
     .      1d0/xmuz*(xmust1*(xmust1-2d0)+1d0)+
     .      2d0*x1-xmust1-3d0)
     .     +(alstor2(1,l)*alstor1(2,l)*(vzz+azz)+
     .       blstor2(1,l)*blstor1(2,l)*(vzz-azz))*(
     .      1d0/xmuz*(1d0-xmust1)*(xmub*(xmust1+x1-1d0)-x1*xmub)
     .      +xmub*(1d0+x1+xmust1-2d0*x2)-2d0*xmust1+2d0*x1*x2+2d0
     .      -2d0*(x1+x2)+xmub*(-x1+2d0) )
     .     +dsqrt(xmub)*xmchar(l)/ast2*
     .      (alstor2(1,l)*blstor1(2,l)*(vzz-azz)+
     .       blstor2(1,l)*alstor1(2,l)*(vzz+azz))*(
     .      1d0/xmuz*(xmust1-1d0)*(xmust1+x1+x2-1d0)+x1-x2-xmust1
     .      +1d0 )
     .     +dsqrt(xmub)*xmchar(l)/ast2*
     .      (alstor2(1,l)*blstor1(2,l)*(vzz+azz)+
     .       blstor2(1,l)*alstor1(2,l)*(vzz-azz))*(
     .      1d0/xmuz*(1d0-xmust1)*(xmust1+x1+x2-1d0)+x1-x2+xmust1
     .      -1d0 ) )
         endif
      enddo
c -------------------------------------------------------------------- c
c                       chargino Higgs interference
c -------------------------------------------------------------------- c
      st2st1hchi=0d0

       do l=1,2
          do I = 1,3
             if ((amchar(l)+mb).gt.ast2.and.
     .(smass(i)+ast1).gt.ast2)then
         st2st1hchi=st2st1hchi-2d0*(Hbbr(I)/dsqrt(2d0)/dchi(l)
     .        /dh(I)/ast2*amz**2/amw*Hstopstopr(I,1,2))*(
     .        (alstor2(1,l)*blstor1(2,l)+alstor1(2,l)*blstor2(1,l))*
     .        xmchar(l)/ast2*(x1+x2+xmust1-1d0-4d0*xmub)
     .        +(alstor2(1,l)*alstor1(2,l)+blstor1(2,l)*blstor2(1,l))*
     .        dsqrt(xmub)*(2d0*x1-1d0+xmust1-4d0*xmub) )
         endif
         enddo
         do J = 1,2
             if ((amchar(l)+mb).gt.ast2.and.
     .(pmass(j)+ast1).gt.ast2)then
               st2st1hchi=st2st1hchi
     .        +2d0*Abbr(J)/dsqrt(2d0)/dchi(l)/da(J)/ast2*amz**2/amw*
     .        (-Astopstopr(J,1,2))*(
     .        (alstor2(1,l)*alstor1(2,l)-blstor1(2,l)*blstor2(1,l))*
     .        dsqrt(xmub)*(1d0-xmust1) +
     .        (alstor2(1,l)*blstor1(2,l)-alstor1(2,l)*blstor2(1,l))*
     .        xmchar(l)/ast2*(1d0-x1-x2-xmust1) )
               endif
         enddo
       enddo

c -------------------------------------------------------------------- c
c                         Z Higgs interference
c -------------------------------------------------------------------- c
       st2st1hz=0.0d0

       do I =1,3
         if ((ast1+mz).gt.ast2.and.
     .(smass(i)+ast1).gt.ast2)then
      st2st1hz=st2st1hz-2d0/2d0/cw*
     .   (aztt12*Hbbr(I)/dsqrt(2d0)*amz**2/amw*Hstopstopr(I,1,2)
     .   /ast2/dz/dh(I))*2d0*dsqrt(xmub)*vzz*2d0*(x1-x2)
          endif
       enddo
       do I =1,2
          if ((ast1+mz).gt.ast2.and.
     .(pmass(i)+ast1).gt.ast2)then
          st2st1hz=st2st1hz
     .   +2d0/2d0/cw*
     .    aztt12*Abbr(I)/dsqrt(2d0)*amz**2/amw*
     .   (-Astopstopr(I,1,2))/ast2/dz/da(I)*
     .   (2d0*dsqrt(xmub)*azz*(2d0/xmuz*(1d0+(xmust1-1d0)*(x1
     .   +x2)+xmust1**2-2d0*xmust1)+2d0-2d0*xmust1) )
          endif
          enddo
c -------------------------------------------------------------------- c

      NS_st2st1bb=st2st1zz+st2st1higgs+st2st1chi+st2st1chiz+st2st1hchi+
     .            st2st1hz

      end
c ==================================================================== c
c =========================== stop1 up upbar ========================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_st2st1uu(x1,x2)
*
      IMPLICIT NONE
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2),gzmumu(2,2)
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION aztt11,aztt12,aztt21,aztt22,xmuz,xmust1,x3,x1,x2
     .,dz,azz,vzz,st2st1zz
*
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau,gzmumu
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      aztt11=gztt(1,1)
      aztt12=gztt(1,2)
      aztt21=gztt(2,1)
      aztt22=gztt(2,2)

      xmuz   = mz**2/ast2**2
      xmust1 = ast1**2/ast2**2

      x3  = 2d0-x1-x2
      dz  = 1d0-x3+xmust1-xmuz

      vzz = vzztoptop
      azz = azztoptop
c -------------------------------------------------------------------- c
c                              Z exchange
c -------------------------------------------------------------------- c
      st2st1zz=0d0
      if ((ast1+mz).gt.ast2)then
      st2st1zz = 1d0/4d0/cw**2/dz**2*aztt12**2*
     .           2d0*(vzz**2+azz**2)*4d0*(1d0+x1*x2-x1-x2-xmust1)
      endif
c ----------------------------------------------------------------- c
      NS_st2st1uu=st2st1zz
      end
c ==================================================================== c
c ========================= stop1 down downbar ======================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_st2st1dd(x1,x2)
*
      IMPLICIT NONE
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2),gzmumu(2,2)
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION aztt11,aztt12,aztt21,aztt22,xmuz,xmust1,x3,x1,x2
     .,dz,azz,vzz,st2st1zz
*
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau,gzmumu
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      aztt11=gztt(1,1)
      aztt12=gztt(1,2)
      aztt21=gztt(2,1)
      aztt22=gztt(2,2)

      xmuz   = mz**2/ast2**2
      xmust1 = ast1**2/ast2**2

      x3  = 2d0-x1-x2
      dz  = 1d0-x3+xmust1-xmuz

      vzz = vzzbotbot
      azz = azzbotbot
c -------------------------------------------------------------------- c
c                              Z exchange
c -------------------------------------------------------------------- c
      st2st1zz =0d0
      if ((ast1+mz).gt.ast2)then
      st2st1zz = 1d0/4d0/cw**2/dz**2*aztt12**2*
     .     2d0*(vzz**2+azz**2)*4d0*(1d0+x1*x2-x1-x2-xmust1)
      endif

c -------------------------------------------------------------------- c
      NS_st2st1dd=st2st1zz

      end
c ==================================================================== c
c ============================= stop1 e+ e- ========================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_st2st1ee(x1,x2)
*
      IMPLICIT NONE
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2),gzmumu(2,2)
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION aztt11,aztt12,aztt21,aztt22,xmuz,xmust1,x3,x1,x2
     .,dz,azz,vzz,st2st1zz
*
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau,gzmumu
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

*
      aztt11=gztt(1,1)
      aztt12=gztt(1,2)
      aztt21=gztt(2,1)
      aztt22=gztt(2,2)

      xmuz   = mz**2/ast2**2
      xmust1 = ast1**2/ast2**2

      x3  = 2d0-x1-x2
      dz  = 1d0-x3+xmust1-xmuz

      vzz = vzztautau
      azz = azztautau
c -------------------------------------------------------------------- c
c                              Z exchange
c -------------------------------------------------------------------- c
      st2st1zz=0d0
      if ((ast1+mz).gt.ast2)then
      st2st1zz = 1d0/4d0/cw**2/dz**2*aztt12**2*
     .     2d0*(vzz**2+azz**2)*4d0*(1d0+x1*x2-x1-x2-xmust1)
      endif
c -------------------------------------------------------------------- c
      NS_st2st1ee=st2st1zz
      end
c ==================================================================== c
c =========================== stop1 nu nubar ========================= c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_st2st1nunu(x1,x2)
*
      IMPLICIT NONE
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2),gzmumu(2,2)
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION aztt11,aztt12,aztt21,aztt22,xmuz,xmust1,x3,x1,x2
     .,dz,azz,vzz,st2st1zz

*
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau,gzmumu
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      aztt11=gztt(1,1)
      aztt12=gztt(1,2)
      aztt21=gztt(2,1)
      aztt22=gztt(2,2)

      xmuz   = mz**2/ast2**2
      xmust1 = ast1**2/ast2**2

      x3  = 2d0-x1-x2
      dz  = 1d0-x3+xmust1-xmuz

      vzz = vzzneutneut
      azz = azzneutneut
c -------------------------------------------------------------------- c
c                              Z exchange
c -------------------------------------------------------------------- c
      st2st1zz =0d0
      if ((ast1+mz).gt.ast2)then
      st2st1zz = 1d0/4d0/cw**2/dz**2*aztt12**2*
     .     2d0*(vzz**2+azz**2)*4d0*(1d0+x1*x2-x1-x2-xmust1)
      endif

c -------------------------------------------------------------------- c

      NS_st2st1nunu=st2st1zz
      end
c ==================================================================== c
c =========================== stop1 tau+ tau- ======================== c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_st2st1tautau(x1,x2)
*
      IMPLICIT NONE
      INTEGER i,j
*
      DOUBLE PRECISION Hstopstopr(3,2,2),Astopstopr(3,2,2)
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2),gzmumu(2,2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION tanbeta_Z!,sw2calc,swcalc
      DOUBLE PRECISION P(2,3)
      DOUBLE PRECISION xmuh(3),xmua(2),dh(3),da(2)
      DOUBLE PRECISION aztt11,aztt12,aztt21,aztt22,xmuz,xmutau,x3,x1,x2,
     .dz,azz,vzz,st2st1zz,bet,b,xmust1,st2st1hk,st2st1aa,st2st1higgs,
     .st2st1hz
      DOUBLE PRECISION AMZ,AMW
*
      COMMON/NS_HIGGSSTST/Hstopstopr,Astopstopr
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup17/azztoptop,vzztoptop,azztautau,vzztautau,
     .                 azzneutneut,vzzneutneut,azzbotbot,vzzbotbot
      COMMON/NS_coup19/gztt,gzbb,gztautau,gzmumu
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_CPodd_MIX/P
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_tanb/tanbeta_Z
      COMMON/NS_MZMWscaleQ/AMZ,AMW

      aztt11=gztt(1,1)
      aztt12=gztt(1,2)
      aztt21=gztt(2,1)
      aztt22=gztt(2,2)

      bet=datan(tanbeta_Z)
      b  =bet

      xmuz   = mz**2/ast2**2
      xmutau = amtau**2/ast2**2
      xmust1 = ast1**2/ast2**2
      do i =1,3
      xmuh(i) = SMASS(i)**2/ast2**2
      enddo
      do i =1,2
      xmua(i) = PMASS(i)**2/ast2**2
      enddo

      x3  = 2d0-x1-x2
      dz  = 1d0-x3+xmust1-xmuz

      do i =1,3
      dh(i) = 1d0-x3+xmust1-xmuh(i)
      enddo
      do i =1,2
      da(i) = 1d0-x3+xmust1-xmua(i)
      enddo

      vzz = vzztautau
      azz = azztautau
c -------------------------------------------------------------------- c
c                              Z exchange
c -------------------------------------------------------------------- c
      st2st1zz =0d0
      if ((ast1+mz).gt.ast2)then
      st2st1zz = 1d0/4d0/cw**2/dz**2*aztt12**2*(
     .   1d0/xmuz**2*(4d0*(vzz**2-azz**2)*xmutau*(1d0-(x1+x2)*(1d0
     .   -xmust1)**2+xmust1*(-xmust1**2+3d0*xmust1-3d0))+
     .   2d0*(vzz**2+azz**2)*(1d0-xmust1)**2*(2d0*xmutau*(xmust1
     .   +x1+x2-1d0)) )
     .   +1d0/xmuz*(8d0*(vzz**2-azz**2)*xmutau*(1d0-xmust1)**2+
     .   4d0*(vzz**2+azz**2)*(1d0-xmust1)*2d0*xmutau*(xmust1-1d0))
     .   +4d0*(vzz**2-azz**2)*xmutau*(x1+x2-xmust1-3d0)
     .   +2d0*(vzz**2+azz**2)*(4d0*(1d0+x1*x2-x1-x2-xmust1)+
     .   xmutau*(-2d0*x1-2d0*x2+2d0*xmust1+6d0)) )
      endif
c -------------------------------------------------------------------- c
c                             Higgs exchange
c -------------------------------------------------------------------- c
       st2st1hk=0.0d0
       st2st1aa=0.0d0
       do I =1,3
         do J =1,3

      if ((smass(i)+ast1).gt.ast2.and.(smass(j)+ast1).gt.ast2)then
      st2st1hk = st2st1hk+(2d0/ast2**2/dh(I)/dh(J)*
     .(-scaltau/dsqrt(2d0)*(-S(i,2)))*(-scaltau/dsqrt(2d0)*(-S(J,2)))*
     .    amz**4/amw**2*Hstopstopr(I,1,2)*Hstopstopr(J,1,2))*
     .    (-1d0+x1+x2+xmust1-4d0*xmutau)
      endif
          enddo
       enddo
       do I =1,2
         do J =1,2

      if ((pmass(i)+ast1).gt.ast2.and.(pmass(j)+ast1).gt.ast2)then
      st2st1aa = st2st1aa+
     .     2d0/ast2**2/da(I)/da(J)
     .*(scaltau/dsqrt(2d0)*P(I,2))*(scaltau/dsqrt(2d0)*P(J,2))*
     .     amz**4/amw**2*Astopstopr(I,1,2)*Astopstopr(J,1,2)
     . *(-1d0+x1+x2+xmust1)
            endif
         enddo
       enddo
      st2st1higgs = st2st1hk+st2st1aa
c -------------------------------------------------------------------- c
c                           Higgs Z interference
c -------------------------------------------------------------------- c
      st2st1hz=0.0d0

      do I = 1,3

      if ((smass(i)+ast1).gt.ast2.and.((ast1+mz).gt.ast2))then
      st2st1hz= st2st1hz-1d0/cw*
     .   (aztt12*(-scaltau/dsqrt(2d0)*(-S(i,2)))
     .*amz**2/amw*Hstopstopr(I,1,2)
     .   /ast2/dz/dh(I))
     . *2d0*dsqrt(xmutau)*vzz*2d0*(x1-x2)
           endif
      enddo
      do J = 1,2

      if ((pmass(j)+ast1).gt.ast2.and.((ast1+mz).gt.ast2))then
      st2st1hz= st2st1hz
     .   + 1d0/cw*
     .    aztt12*(-scaltau/dsqrt(2d0)*(-P(J,2)))*amz**2/amw*
     .   (-Astopstopr(J,1,2))/ast2/dz/da(j)*
     .   (2d0*dsqrt(xmutau)*azz*(2d0/xmuz*(1d0+(xmust1-1d0)*(x1
     .   +x2)+xmust1**2-2d0*xmust1)+2d0-2d0*xmust1) )
         endif
      enddo
c -------------------------------------------------------------------- c
      NS_st2st1tautau=st2st1zz+st2st1higgs+st2st1hz
      end
c =======================================================================c
