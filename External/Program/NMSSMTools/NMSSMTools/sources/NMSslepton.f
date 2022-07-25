      SUBROUTINE NS_SLEPTON

************************************************************************
*
*     This subroutine computes the slepton decays
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

      INTEGER I,J,GRFLAG

      DOUBLE PRECISION NS_lamb,PI,SQR2,multilim
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION sw,cw,tw,flagmulti,flagqcd,flagloop
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .     blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION gztt(2,2),gzbb(2,2),gztautau(2,2),gzmumu(2,2)
      DOUBLE PRECISION gwtb(2,2),gwntau(2,2),gwnmu(2,2)
      DOUBLE PRECISION gcsntaur(2,2)
      DOUBLE PRECISION Hstaustaur(3,2,2),Astaustaur(2,2,2)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION asne2,asnmu2,asntau2
      DOUBLE PRECISION gmsntau(2),gmstau(2)
      DOUBLE PRECISION sellneute(5),selrneute(5),sellcharnue(2),
     .         selrcharnue(2),snellneut(5),snellchar(2),
     .         smu1neutmu(5),smu2neutmu(5),smu1charnumu(2),
     .         smu2charnumu(2),snmu1neut(5),snmu1char(2),
     .         stau1neut(5),stau2neut(5),stau1char(2),
     .         stau2char(2),stau1hcsn(2),stau2hcsn(2),
     .         stau1wsn(2),stau2wsn(2),stau2ztau,stau2H(3),
     .         stau2A(2),sntauneut(5),sntauchar(2),
     .         sntau1hcstau(2),sntau1wstau(2)
      DOUBLE PRECISION selltot,selltot2,selltot3,selrtot,selrtot2,
     .         selrtot3,smu1tot,smu1tot2,smu1tot3,smu2tot,smu2tot2,
     .         smu2tot3,stau1tot2,stau2tot,stau2tot2,stau2tot3,
     .         snelltot,snelltot2,snelltot3,snmu1tot,snmu1tot2,
     .         snmu1tot3,sntautot,sntautot2,sntautot3
      DOUBLE PRECISION brsellneute(5),brselrneute(5),brsellcharnue(2),
     .         brselrcharnue(2),brsnellneut(5),brsnellchar(5),
     .         brsmu1neutmu(5),brsmu2neutmu(5),brsmu1charnumu(2),
     .         brsmu2charnumu(2),brsnmu1neut(5),brsnmu1char(5),
     .         brstau1neut(5),brstau2neut(5),brstau1char(2),
     .         brstau2char(2),brstau1hcsn(2),brstau2hcsn(2),
     .         brstau1wsn(2),brstau2wsn(2),brstau2ztau,brstau2H(3),
     .         brstau2A(2),brsntauneut(5),brsntauchar(2),
     .         brsntau1hcstau(2),brsntau1wstau(2)
      DOUBLE PRECISION brsellstau1star,brsellstau1,
     .         brsellstau1nutau,brselrstau1star,brselrstau1,
     .         brselrstau1nutau,brsnestau1star,brsnestau1,
     .         brsnestau1nutau,brsmu1stau1star,brsmu1stau1,
     .         brsmu1stau1nutau,brsmu2stau1star,brsmu2stau1,
     .         brsmu2stau1nutau,brsnmustau1star,brsnmustau1,
     .         brsnmustau1nutau,brstau2stau1star,brstau2stau1,
     .         brstau2stau1nn,brsntaustau1star,brsntaustau1,
     .         brsntaustau1nutau
      DOUBLE PRECISION selEgra,serEgra,smu1MUgra,smu2MUgra,
     .         stau1TAUgra,stau2TAUgra,sneNEgra,snmNMgra,sntNTgra
      DOUBLE PRECISION xintegsellstau1star,xintegsellstau1,
     .         xintegsellstau1nutau
      DOUBLE PRECISION xintegselrstau1star,xintegselrstau1,
     .         xintegselrstau1nutau
       DOUBLE PRECISION xintegsnestau1star,xintegsnestau1,
     .         xintegsnestau1nutau
      DOUBLE PRECISION xintegsmu1stau1star,xintegsmu1stau1,
     .         xintegsmu1stau1nutau
      DOUBLE PRECISION xintegsmu2stau1star,xintegsmu2stau1,
     .         xintegsmu2stau1nutau
       DOUBLE PRECISION xintegsnmustau1star,xintegsnmustau1,
     .         xintegsnmustau1nutau
      DOUBLE PRECISION xintegstau2stau1star,xintegstau2stau1,
     .         xintegstau2stau1nn
      DOUBLE PRECISION xintegsntaustau1star,xintegsntaustau1,
     .         xintegsntaustau1nutau
      DOUBLE PRECISION brcharWgra(2),brcharHCgra(2),brneutGAMgra(5),
     .         brneutZgra(5),brneutHgra(5,3),brneutAgra(5,2),
     .         brgluiGLUgra,brselEgra,brserEgra,brsmu1MUgra,
     .         brsmu2MUgra,brstau1TAUgra,brstau2TAUgra,brsneNEgra,
     .         brsnmNMgra,brsntNTgra,brsulUgra,brsurUgra,brsdlDgra,
     .         brsdrDgra,brst1Tgra,brst2Tgra,brsb1Bgra,brsb2Bgra
      DOUBLE PRECISION KNG(5),KNZ(5),KNH(5,3),KNA(5,2),KCW(2),KCH(2)
      DOUBLE PRECISION M32,CGR,MPL

      COMMON/SLEPTON_2GAMMA/sellneute,selrneute,sellcharnue,selrcharnue,
     .         snellneut,snellchar,smu1neutmu,smu2neutmu,smu1charnumu,
     .         smu2charnumu,snmu1neut,snmu1char,stau1neut,stau2neut,
     .         stau1char,stau2char,stau1hcsn,stau2hcsn,stau1wsn,
     .         stau2wsn,stau2ztau,stau2H,stau2A,sntauneut,sntauchar,
     .         sntau1hcstau,sntau1wstau
      COMMON/SLEPTON_WIDTH/selltot,selltot2,selltot3,selrtot,
     .         selrtot2,selrtot3,smu1tot,smu1tot2,smu1tot3,smu2tot,
     .         smu2tot2,smu2tot3,stau1tot2,stau2tot,stau2tot2,
     .         stau2tot3,snelltot,snelltot2,snelltot3,snmu1tot,
     .         snmu1tot2,snmu1tot3,sntautot2,sntautot3,sntautot
      COMMON/SLEPTON_BR_2BD/brsellneute,brselrneute,brsellcharnue,
     .         brselrcharnue,brsnellneut,brsnellchar,brsmu1neutmu,
     .         brsmu2neutmu,brsmu1charnumu,brsmu2charnumu,brsnmu1neut,
     .         brsnmu1char,brstau1neut,brstau2neut,brstau1char,
     .         brstau2char,brstau1hcsn,brstau2hcsn,brstau1wsn,
     .         brstau2wsn,brstau2ztau,brstau2H,brstau2A,brsntauneut,
     .         brsntauchar,brsntau1hcstau,brsntau1wstau
      COMMON/SLEPTON_BR_3BD/brsellstau1star,brsellstau1,
     .         brsellstau1nutau,brselrstau1star,brselrstau1,
     .         brselrstau1nutau,brsnestau1star,brsnestau1,
     .         brsnestau1nutau,brsmu1stau1star,brsmu1stau1,
     .         brsmu1stau1nutau,brsmu2stau1star,brsmu2stau1,
     .         brsmu2stau1nutau,brsnmustau1star,brsnmustau1,
     .         brsnmustau1nutau,brstau2stau1star,brstau2stau1,
     .         brstau2stau1nn,brsntaustau1star,brsntaustau1,
     .         brsntaustau1nutau
      COMMON/GRAVITINO/brcharWgra,brcharHCgra,brneutGAMgra,
     .         brneutZgra,brneutHgra,brneutAgra,
     .         brgluiGLUgra,brselEgra,brserEgra,brsmu1MUgra,
     .         brsmu2MUgra,brstau1TAUgra,brstau2TAUgra,brsneNEgra,
     .         brsnmNMgra,brsntNTgra,brsulUgra,brsurUgra,brsdlDgra,
     .         brsdrDgra,brst1Tgra,brst2Tgra,brsb1Bgra,brsb2Bgra
      COMMON/GRAVCOUP/KNG,KNZ,KNH,KNA,KCW,KCH
      COMMON/M32/M32,CGR,MPL,GRFLAG
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_rhsneutr/asne2,asnmu2,asntau2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
      COMMON/NS_coup19/gztt,gzbb,gztautau,gzmumu
      COMMON/NS_coup20/gwtb,gwntau,gwnmu
      COMMON/NS_HIGGSTAUTAU/Hstaustaur,Astaustaur
      COMMON/NS_hcstausntau/gcsntaur
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_multilim/multilim
      COMMON/NS_FLAGS/flagmulti,flagqcd,flagloop

      EXTERNAL NS_lamb

*  Initialization

      selltot  = 0d0
      selltot2 = 0d0
      selltot3 = 0d0
      selrtot = 0d0
      selrtot2 = 0d0
      selrtot3 = 0d0
      do i=1,2
         sellcharnue(i) = 0d0
         selrcharnue(i) = 0d0
      enddo
      do j=1,5
         sellneute(j) = 0d0
         selrneute(j) = 0d0
      enddo

      snelltot =0d0
      snelltot2 =0d0
      snelltot3 =0d0
      do i=1,2
         snellchar(i) = 0d0
      enddo
      do j=1,5
         snellneut(j) = 0d0
      enddo

      smu1tot  = 0d0
      smu1tot2 = 0d0
      smu1tot3 = 0d0
      smu2tot = 0d0
      smu2tot2 = 0d0
      smu2tot3 = 0d0
      do i=1,2
         smu1charnumu(i) = 0d0
         smu2charnumu(i) = 0d0
      enddo
      do j=1,5
        smu1neutmu(j) = 0d0
        smu2neutmu(j) = 0d0
      enddo

      snmu1tot =0d0
      snmu1tot2 =0d0
      snmu1tot3 =0d0
      do i=1,2
         snmu1char(i) = 0d0
      enddo
      do j=1,5
         snmu1neut(j) = 0d0
      enddo

      stau1tot2 = 0d0
      stau2tot = 0d0
      stau2tot2 = 0d0
      stau2tot3 = 0d0
      do j=1,5
         stau1neut(j) = 0d0
         stau2neut(j) = 0d0
      enddo
      do i=1,2
         stau1char(i) = 0d0
         stau1hcsn(i) = 0d0
         stau1wsn(i)  = 0d0
         stau2char(i) = 0d0
         stau2hcsn(i) = 0d0
         stau2wsn(i)  = 0d0
      enddo
      do i=1,3
      stau2H(i) = 0d0
      enddo
      do i=1,2
      stau2A(i) = 0d0
      enddo
      stau2ztau = 0d0

      sntautot=0d0
      sntautot2=0d0
      sntautot3=0d0
      do i=1,2
         sntauchar(i)    = 0d0
         sntau1hcstau(i) = 0d0
         sntau1wstau(i)  = 0d0
      enddo
      do j=1,5
         sntauneut(j) = 0d0
      enddo

      selEgra  = 0d0
      serEgra  = 0d0
      smu1MUgra = 0d0
      smu2MUgra = 0d0
      stau1TAUgra = 0d0
      stau2TAUgra = 0d0
      sneNEgra  = 0d0
      snmNMgra  = 0d0
      sntNTgra  = 0d0

c ==================================================================== c
c                        selectron 2-body decays                       c
c ==================================================================== c

c  selectron_L --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + e-

      do i=1,5
         if(ase1.gt.amneut(i)) then
            sellneute(i)=g2s*ae(1,i)**2*
     .           (ase1**2-amneut(i)**2)*
     .           NS_lamb(0d0,amneut(i)/ase1)/(16*pi*ase1)
         else
            sellneute(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  selectron_L --> chi-_1/chi-_2 + neutrino_e

      do i=1,2
         if(ase1.gt.amchar(i)) then
            sellcharnue(i)=g2s*ale(1,i)**2*(ase1**2-amchar(i)**2)*
     .                     NS_lamb(0d0,amchar(i)/ase1)/(16*pi*ase1)
         else
            sellcharnue(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  selectron_R --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + e-

      do i=1,5
         if(ase2.gt.amneut(i)) then
            selrneute(i)=g2s*be(2,i)**2*
     .           (ase2**2-amneut(i)**2)*
     .           NS_lamb(0d0,amneut(i)/ase2)/(16*pi*ase2)
         else
            selrneute(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  selectron_R --> chi-_1/chi-_2 + neutrino_e

      do i=1,2
         if(ase2.gt.amchar(i)) then
            selrcharnue(i)=g2s*ale(2,i)**2*(ase2**2-amchar(i)**2)*
     .                     NS_lamb(0d0,amchar(i)/ase2)/(16*pi*ase2)
         else
            selrcharnue(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  selectron_L --> e- + gravitino

      if(GRFLAG.EQ.1)then
        if (M32.le.ase1) then
          selEgra = ase1**5/(48d0*PI*MPL**2*M32**2)
        else
          selEgra = 0d0
        endif
      endif
c -------------------------------------------------------------------- c
c  selectron_R --> e- + gravitino

      if(GRFLAG.EQ.1)then
        if (M32.le.ase2) then
          serEgra = ase2**5/(48d0*PI*MPL**2*M32**2)
        else
          serEgra = 0d0
        endif
      endif
c -------------------------------------------------------------------- c
c                  SUM 2 BODY DECAYS
c -------------------------------------------------------------------- c

      do i=1,5
        selltot2=selltot2+sellneute(i)
        selrtot2=selrtot2+selrneute(i)
      enddo

      do i=1,2
        selltot2=selltot2+sellcharnue(i)
        selrtot2=selrtot2+selrcharnue(i)
      enddo

      selltot2=selltot2+selEgra
      selrtot2=selrtot2+serEgra

c--------------------------------------------------------------------- c
c ----- selectron_L 3-body decays and 3-body total widths ------------ c
c--------------------------------------------------------------------- c
c Only: selectron_L -> stau_1* + e + tau      (xintegsellstau1star)
c       selectron_L -> stau_1 + e + tau       (xintegsellstau1)
c       selectron_L -> stau_1 + nu_e + nu_tau (xintegsellstau1nutau)

      if(flagmulti.eq.1d0) then
        CALL NS_xintegsell(xintegsellstau1star,xintegsellstau1,
     .    xintegsellstau1nutau)
        selltot3=xintegsellstau1star+xintegsellstau1
     .    +xintegsellstau1nutau
        if(selltot3.lt.multilim*selltot2)then
           selltot3=0d0
        endif
      endif
      selltot=selltot2+selltot3

c--------------------------------------------------------------------- c
c ----- selectron_R 3-body decays and 3-body total widths ------------ c
c--------------------------------------------------------------------- c
c Only: selectron_R -> stau_1* + e + tau      (xintegselrstaustarltau)
c       selectron_R -> stau_1 + e + tau       (xintegselrstaultau)
c       selectron_R -> stau_1 + nu_e + nu_tau (xintegselrstau1nutau)

      if(flagmulti.eq.1d0) then
        CALL NS_xintegselr(xintegselrstau1star,xintegselrstau1,
     .    xintegselrstau1nutau)
        selrtot3=xintegselrstau1star+xintegselrstau1
     .    +xintegselrstau1nutau
        if(selrtot3.lt.multilim*selrtot2)then
           selrtot3=0d0
        endif
      endif
      selrtot=selrtot2+selrtot3

c------------------------------------------------ c
c ----- selectron_L branching ratios ------------ c
c------------------------------------------------ c

      if(selltot.ne.0d0)then

       do i=1,5
         brsellneute(i) = sellneute(i)/selltot
       enddo
       do i=1,2
         brsellcharnue(i) = sellcharnue(i)/selltot
       enddo
       brselEgra = selEgra/selltot

      else

       do i=1,5
         brsellneute(i) = 0d0
       enddo
       do i=1,2
         brsellcharnue(i) = 0d0
       enddo
       brselEgra = 0d0

      endif

      if (selltot3.ne.0d0)then

         brsellstau1star  = xintegsellstau1star/selltot
         brsellstau1      = xintegsellstau1/selltot
         brsellstau1nutau = xintegsellstau1nutau/selltot

      else

         brsellstau1star  = 0d0
         brsellstau1      = 0d0
         brsellstau1nutau = 0d0

      endif

c------------------------------------------------ c
c ----- selectron_R branching ratios ------------ c
c------------------------------------------------ c

      if(selrtot.ne.0d0)then

       do i=1,5
         brselrneute(i) = selrneute(i)/selrtot
       enddo
       do i=1,2
         brselrcharnue(i) = selrcharnue(i)/selrtot
       enddo
       brserEgra = serEgra/selrtot

      else

       do i=1,5
         brselrneute(i) = 0d0
       enddo
       do i=1,2
         brselrcharnue(i) = 0d0
       enddo
       brserEgra = 0d0

      endif

      if(selrtot3.ne.0d0) then

         brselrstau1star  = xintegselrstau1star/selrtot
         brselrstau1      = xintegselrstau1/selrtot
         brselrstau1nutau = xintegselrstau1nutau/selrtot

      else

         brselrstau1star  = 0d0
         brselrstau1      = 0d0
         brselrstau1nutau = 0d0

      endif

c ==================================================================== c
c                       sneutrino_el 2-body decays                     c
c ==================================================================== c
c  sneutrino_el --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + neutrino_e

      do i=1,5
         if(asne1.gt.amneut(i)) then
            snellneut(i)=g2s*(ane(1,i)**2+bne(2,i)**2)*
     .        (asne1**2-amneut(i)**2)*NS_lamb(0d0,amneut(i)/asne1)
     .        /(16*pi*asne1)
         else
            snellneut(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  sneutrino_el --> chi+_1/chi+_2 + e-

      do i=1,2
         if(asne1.gt.amchar(i)) then
           snellchar(i)=g2s*alsne(1,i)**2*(asne1**2-amchar(i)**2)*
     .           NS_lamb(0d0,amchar(i)/asne1)/(16*pi*asne1)
         else
            snellchar(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  sneutrino_el --> neutrino_e + gravitino

      if(GRFLAG.EQ.1)then
        if (M32.le.asne1) then
          sneNEgra = asne1**5/(48d0*PI*MPL**2*M32**2)
        else
          sneNEgra = 0d0
        endif
      endif

c -------------------------------------------------------------------- c
c                  SUM 2 BODY DECAYS
c -------------------------------------------------------------------- c

      do i=1,5
         snelltot2=snelltot2+snellneut(i)
      enddo
      do i=1,2
         snelltot2=snelltot2+snellchar(i)
      enddo
      snelltot2=snelltot2+sneNEgra

c--------------------------------------------------------------------- c
c -------- sneutrino_el 3-body decays and 3-body total widths -------- c
c--------------------------------------------------------------------- c
c Only: sne -> stau1* + nu_e + tau* (xintegsnestau1star)
c       sne -> stau1 + nu_e + tau   (xintegsnestau1)
c       sne -> stau1 + e + nu_tau   (xintegsnestau1nutau)
      if(flagmulti.eq.1d0) then
        CALL NS_xintegsne(xintegsnestau1star,xintegsnestau1,
     .    xintegsnestau1nutau)
        snelltot3=xintegsnestau1star+xintegsnestau1+
     .    xintegsnestau1nutau
        if (snelltot3.lt.multilim*snelltot2) then
           snelltot3=0d0
        endif
      endif
      snelltot=snelltot2+snelltot3

c------------------------------------------------ c
c ----- sneutrino_el branching ratios ----------- c
c------------------------------------------------ c

      if(snelltot.ne.0d0)then

       do i=1,5
         brsnellneut(i) = snellneut(i)/snelltot
       enddo
       do i=1,2
         brsnellchar(i) = snellchar(i)/snelltot
       enddo
       brsneNEgra = sneNEgra/snelltot

      else

       do i=1,5
         brsnellneut(i) = 0d0
       enddo
       do i=1,2
         brsnellchar(i) = 0d0
       enddo
       brsneNEgra = 0d0

      endif

      if (snelltot3.ne.0d0) then

         brsnestau1star  = xintegsnestau1star/snelltot
         brsnestau1      = xintegsnestau1/snelltot
         brsnestau1nutau = xintegsnestau1nutau/snelltot

      else

        brsnestau1star  = 0d0
        brsnestau1      = 0d0
        brsnestau1nutau = 0d0

      endif

c ==================================================================== c
c                        smuon 2-body decays                           c
c ==================================================================== c

c  smuon_1 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + mu-

      do i=1,5
         if(asmu1.gt.amneut(i)) then
            smu1neutmu(i)=g2s*amu(1,i)**2*
     .           (asmu1**2-amneut(i)**2)*
     .           NS_lamb(0d0,amneut(i)/asmu1)/(16*pi*asmu1)
         else
            smu1neutmu(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  smuon_1 --> chi-_1/chi-_2 + neutrino_mu

      do i=1,2
         if(asmu1.gt.amchar(i)) then
            smu1charnumu(i)=g2s*almu(1,i)**2*(asmu1**2-amchar(i)**2)*
     .                     NS_lamb(0d0,amchar(i)/asmu1)/(16*pi*asmu1)
         else
            smu1charnumu(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  smuon_2 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + mu-

      do i=1,5
         if(asmu2.gt.amneut(i)) then
            smu2neutmu(i)=g2s*bmu(2,i)**2*
     .           (asmu2**2-amneut(i)**2)*
     .           NS_lamb(0d0,amneut(i)/asmu2)/(16*pi*asmu2)
         else
            smu2neutmu(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  smuon_2 --> chi-_1/chi-_2 + neutrino_mu

      do i=1,2
         if(asmu2.gt.amchar(i)) then
            smu2charnumu(i)=g2s*almu(2,i)**2*(asmu2**2-amchar(i)**2)*
     .                     NS_lamb(0d0,amchar(i)/asmu2)/(16*pi*asmu2)
         else
            smu2charnumu(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  smuon_1 --> mu- + gravitino

      if(GRFLAG.EQ.1)then
        if (M32.le.asmu1) then
          smu1MUgra = asmu1**5/(48d0*PI*MPL**2*M32**2)
        else
          smu1MUgra = 0d0
        endif
      endif
c -------------------------------------------------------------------- c
c  smuon_2 --> mu- + gravitino

      if(GRFLAG.EQ.1)then
        if (M32.le.asmu2) then
          smu2MUgra = asmu2**5/(48d0*PI*MPL**2*M32**2)
        else
          smu2MUgra = 0d0
        endif
      endif
c -------------------------------------------------------------------- c
c                  SUM 2 BODY DECAYS
c -------------------------------------------------------------------- c

      do i=1,5
        smu1tot2=smu1tot2+smu1neutmu(i)
        smu2tot2=smu2tot2+smu2neutmu(i)
      enddo

      do i=1,2
        smu1tot2=smu1tot2+smu1charnumu(i)
        smu2tot2=smu2tot2+smu2charnumu(i)
      enddo

      smu1tot2=smu1tot2+smu1MUgra
      smu2tot2=smu2tot2+smu2MUgra

c--------------------------------------------------------------------- c
c ----- smuon_1 3-body decays and 3-body total widths ---------------- c
c--------------------------------------------------------------------- c
c Only: smuon_1 -> stau_1* + mu + tau      (xintegsmu1stau1star)
c       smuon_1 -> stau_1 + mu + tau       (xintegsmu1stau1)
c       smuon_1 -> stau_1 + nu_mu + nu_tau (xintegsmu1stau1nutau)

      if(flagmulti.eq.1d0) then
        CALL NS_xintegsmu1(xintegsmu1stau1star,xintegsmu1stau1,
     .    xintegsmu1stau1nutau)
        smu1tot3=xintegsmu1stau1star+xintegsmu1stau1
     .    +xintegsmu1stau1nutau
        if(smu1tot3.lt.multilim*smu1tot2)then
           smu1tot3=0d0
        endif
      endif
      smu1tot=smu1tot2+smu1tot3

c--------------------------------------------------------------------- c
c ----- smuon_2 3-body decays and 3-body total widths ---------------- c
c--------------------------------------------------------------------- c
c Only: smuon_2 -> stau_1* + mu + tau      (xintegsmu2staustarltau)
c       smuon_2 -> stau_1 + mu + tau       (xintegsmu2staultau)
c       smuon_2 -> stau_1 + nu_mu + nu_tau (xintegsmu2stau1nutau)

      if(flagmulti.eq.1d0) then
        CALL NS_xintegsmu2(xintegsmu2stau1star,xintegsmu2stau1,
     .    xintegsmu2stau1nutau)
        smu2tot3=xintegsmu2stau1star+xintegsmu2stau1
     .    +xintegsmu2stau1nutau
        if(smu2tot3.lt.multilim*smu2tot2)then
           smu2tot3=0d0
        endif
      endif
      smu2tot=smu2tot2+smu2tot3

c------------------------------------------------ c
c ----- smuon_1 branching ratios ---------------- c
c------------------------------------------------ c

      if(smu1tot.ne.0d0)then

       do i=1,5
         brsmu1neutmu(i) = smu1neutmu(i)/smu1tot
       enddo
       do i=1,2
         brsmu1charnumu(i) = smu1charnumu(i)/smu1tot
       enddo
       brsmu1MUgra = smu1MUgra/smu1tot
 
      else

       do i=1,5
         brsmu1neutmu(i)= 0d0
       enddo
       do i=1,2
         brsmu1charnumu(i) = 0d0
       enddo
       brsmu1MUgra = 0d0

      endif

      if (smu1tot3.ne.0d0)then

         brsmu1stau1star  = xintegsmu1stau1star/smu1tot
         brsmu1stau1      = xintegsmu1stau1/smu1tot
         brsmu1stau1nutau = xintegsmu1stau1nutau/smu1tot
 
      else

         brsmu1stau1star  = 0d0
         brsmu1stau1      = 0d0
         brsmu1stau1nutau = 0d0

      endif

c------------------------------------------------ c
c ----- smuon_2 branching ratios ---------------- c
c------------------------------------------------ c

      if(smu2tot.ne.0d0)then

       do i=1,5
         brsmu2neutmu(i) = smu2neutmu(i)/smu2tot
       enddo
       do i=1,2
         brsmu2charnumu(i) = smu2charnumu(i)/smu2tot
       enddo
       brsmu2MUgra = smu2MUgra/smu2tot
 
      else

       do i=1,5
         brsmu2neutmu(i) = 0d0
       enddo
       do i=1,2
         brsmu2charnumu(i) = 0d0
       enddo
       brsmu2MUgra = 0d0
 
      endif

      if(smu2tot3.ne.0d0)then

         brsmu2stau1star  = xintegsmu2stau1star/smu2tot
         brsmu2stau1      = xintegsmu2stau1/smu2tot
         brsmu2stau1nutau = xintegsmu2stau1nutau/smu2tot
 
      else

         brsmu2stau1star  = 0d0
         brsmu2stau1      = 0d0
         brsmu2stau1nutau = 0d0

      endif

c ==================================================================== c
c                       sneutrino_mu1 2-body decays                    c
c ==================================================================== c
c  sneutrino_mu1 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + neutrino_mu

      do i=1,5
         if(asnmu1.gt.amneut(i)) then
            snmu1neut(i)=g2s*(anmu(1,i)**2+bnmu(2,i)**2)*
     .        (asnmu1**2-amneut(i)**2)*NS_lamb(0d0,amneut(i)/asnmu1)
     .        /(16*pi*asnmu1)
         else
            snmu1neut(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  sneutrino_mu1 --> chi+_1/chi+_2 + mu-

      do i=1,2
         if(asnmu1.gt.amchar(i)) then
           snmu1char(i)=g2s*alsnm(1,i)**2*(asnmu1**2-amchar(i)**2)*
     .           NS_lamb(0d0,amchar(i)/asnmu1)/(16*pi*asnmu1)
         else
            snmu1char(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  sneutrino_mu1 --> neutrino_mu + gravitino

      if(GRFLAG.EQ.1)then
        if (M32.le.asnmu1) then
          snmNMgra = asnmu1**5/(48d0*PI*MPL**2*M32**2)
        else
          snmNMgra = 0d0
        endif
      endif

c -------------------------------------------------------------------- c
c                  SUM 2 BODY DECAYS
c -------------------------------------------------------------------- c

      do i=1,5
         snmu1tot2=snmu1tot2+snmu1neut(i)
      enddo
      do i=1,2
         snmu1tot2=snmu1tot2+snmu1char(i)
      enddo
      snmu1tot2=snmu1tot2+snmNMgra

c--------------------------------------------------------------------- c
c -------- sneutrino_mu1 3-body decays and 3-body total widths ------- c
c--------------------------------------------------------------------- c
c Only: snmu -> stau1* + nu_mu + tau  (xintegsnmustau1star)
c       snmu -> stau1 + nu_mu + tau   (xintegsnmustau1)
c       snmu -> stau1 + mu + nu_tau   (xintegsnmustau1nutau)
      if(flagmulti.eq.1d0) then
        CALL NS_xintegsnmu(xintegsnmustau1star,xintegsnmustau1,
     .    xintegsnmustau1nutau)
        snmu1tot3=xintegsnmustau1star+xintegsnmustau1+
     .    xintegsnmustau1nutau
        if (snmu1tot3.lt.multilim*snmu1tot2) then
           snmu1tot3=0d0
        endif
      endif
      snmu1tot=snmu1tot2+snmu1tot3

c------------------------------------------------ c
c ----- sneutrino_mu1 branching ratios ---------- c
c------------------------------------------------ c

      if(snmu1tot.ne.0d0)then

       do i=1,5
         brsnmu1neut(i) = snmu1neut(i)/snmu1tot
       enddo
       do i=1,2
         brsnmu1char(i) = snmu1char(i)/snmu1tot
       enddo
       brsnmNMgra = sneNEgra/snmu1tot

      else

       do i=1,5
         brsnmu1neut(i) = 0d0
       enddo
       do i=1,2
         brsnmu1char(i) = 0d0
       enddo
       brsnmNMgra = 0d0

      endif

      if (snmu1tot3.ne.0d0) then

         brsnmustau1star  = xintegsnmustau1star/snmu1tot
         brsnmustau1      = xintegsnmustau1/snmu1tot
         brsnmustau1nutau = xintegsnmustau1nutau/snmu1tot

      else

          brsnmustau1star  = 0d0
          brsnmustau1      = 0d0
          brsnmustau1nutau = 0d0

      endif

c ==================================================================== c
c                       stau1/2 2-body decays                          c
c ==================================================================== c

      gmsntau(1) = asntau1
      gmsntau(2) = asntau2

c -------------------------------------------------------------------- c
c  stau1 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + tau

      do i=1,5
         if(astau1.gt.(amneut(i)+mtau)) then
            stau1neut(i)=g2s*((atau(1,i)**2+btau(1,i)**2)*
     .           (astau1**2-amneut(i)**2-mtau**2)
     .           -4*atau(1,i)*btau(1,i)*mtau*xmneut(i)
     .           )*NS_lamb(mtau/astau1,amneut(i)/astau1)
     .        /(16d0*pi*astau1)
         else
            stau1neut(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  stau2 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + tau

      do i=1,5
         if(astau2.gt.(amneut(i)+mtau)) then
            stau2neut(i)=g2s*((atau(2,i)**2+btau(2,i)**2)*
     .           (astau2**2-amneut(i)**2-mtau**2)
     .           -4*atau(2,i)*btau(2,i)*mtau*xmneut(i)
     .           )*NS_lamb(mtau/astau2,amneut(i)/astau2)
     .           /(16d0*pi*astau2)
         else
            stau2neut(i)=0d0
         endif
      enddo

c ----------------------------------------------------------------- c
c  stau1 --> chi+_1/chi+_2 + nu_tau
      do i=1,2
         if(astau1.gt.amchar(i)) then
             stau1char(i)=g2s*altau(1,i)**2*(astau1**2-amchar(i)**2)*
     .        NS_lamb(0d0,amchar(i)/astau1)/(16*pi*astau1)
         else
            stau1char(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  stau2 --> chi-_1/chi-_2 + nu_tau
      do i=1,2
         if(astau2.gt.amchar(i)) then
            stau2char(i)=g2s*altau(2,i)**2*(astau2**2-amchar(i)**2)*
     .        NS_lamb(0d0,amchar(i)/astau2)/(16*pi*astau2)
         else
            stau2char(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  stau1 --> H- + sneutrino_tau1/2
      do i=1,2
         if(astau1.gt.(gmsntau(i)+cmass)) then
            stau1hcsn(i)=g2s*mw**2*gcsntaur(i,1)**2*
     .           NS_lamb(gmsntau(i)/astau1,cmass/astau1)
     .           /(16d0*pi*astau1)
         else
            stau1hcsn(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  stau2 --> H- + sneutrino_tau1/2
      do i=1,2
         if(astau2.gt.(gmsntau(i)+cmass)) then
            stau2hcsn(i)=g2s*mw**2*gcsntaur(i,2)**2*
     .           NS_lamb(gmsntau(i)/astau2,cmass/astau2)
     .           /(16d0*pi*astau2)
         else
            stau2hcsn(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  stau2 --> h(i) + stau1 [stau2H(3)]
      do I=1,3
      if(astau2.gt.(astau1+SMASS(I))) then
         stau2H(i)=g2s*mz**4/mw**2*Hstaustaur(I,2,1)**2*
     .         NS_lamb(astau1/astau2,SMASS(I)/astau2)/(16d0*pi*astau2)
      else
         stau2H(I)=0d0
      endif
      enddo

c -------------------------------------------------------------------- c
c  stau2 --> A(i) + stau1 [stau2A(2)]
      do I=1,2
      if(astau2.gt.(astau1+PMASS(I))) then
            stau2A(I)=g2s*mz**4/mw**2*Astaustaur(I,2,1)**2*
     .         NS_lamb(astau1/astau2,PMASS(I)/astau2)/(16d0*pi*astau2)
         else
           stau2A(I)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  stau2 --> Z + stau1

      if(astau2.gt.(astau1+mz)) then
        stau2ztau=g2s/64d0/pi/cw**2/mz**2*astau2**3*gztautau(2,1)**2*
     .           NS_lamb(astau1/astau2,mz/astau2)**3
      else
         stau2ztau=0d0
      endif

c -------------------------------------------------------------------- c
c  stau1 --> W- + sntau1/2
      do i=1,2
         if(astau1.gt.(gmsntau(i)+mw)) then
            stau1wsn(i)=g2s/32d0/pi/mw**2*astau1**3*gwntau(i,1)**2*
     .                NS_lamb(gmsntau(i)/astau1,mw/astau1)**3
         else
            stau1wsn(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  stau2 --> W- + sntau1/2
      do i=1,2
         if(astau2.gt.(gmsntau(i)+mw)) then
            stau2wsn(i)=g2s/32d0/pi/mw**2*astau2**3*gwntau(i,2)**2*
     .                NS_lamb(gmsntau(i)/astau2,mw/astau2)**3
         else
            stau2wsn(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  stau1 --> tau + gravitino

      if(GRFLAG.EQ.1)then
        if ((mtau+M32).le.astau1) then
          stau1TAUgra = astau1**5/(48d0*PI*MPL**2*M32**2)
     .                *(1d0-mtau**2/astau1**2)**4
        else
          stau1TAUgra = 0d0
        endif
      endif

c -------------------------------------------------------------------- c
c  stau2 --> tau + gravitino

      if(GRFLAG.EQ.1)then
        if ((mtau+M32).le.astau2) then
          stau2TAUgra = astau2**5/(48d0*PI*MPL**2*M32**2)
     .                *(1d0-mtau**2/astau2**2)**4
        else
          stau2TAUgra = 0d0
        endif
      endif

c -------------------------------------------------------------------- c
c                  SUM 2 BODY DECAYS
c -------------------------------------------------------------------- c

      stau1tot2=stau1neut(1)+stau1neut(2)+stau1neut(3)+stau1neut(4)+
     .          stau1neut(5)+
     .          stau1char(1)+stau1char(2)+stau1hcsn(1)+stau1hcsn(2)+
     .          stau1wsn(1)+stau1wsn(2)+stau1TAUgra

      stau2tot2=stau2neut(1)+stau2neut(2)+stau2neut(3)+stau2neut(4)+
     .          stau2neut(5)+
     .          stau2char(1)+stau2char(2)+stau2hcsn(1)+stau2hcsn(2)+
     .          stau2wsn(1)+stau2wsn(2)+stau2H(1)+stau2H(2)+stau2H(3)+
     .          stau2A(1)+stau2A(2)+stau2ztau+stau2TAUgra

c--------------------------------------------------------------------- c
c ----- stau_2 3-body decays and 3-body total widths ----------------- c
c--------------------------------------------------------------------- c
c Only: stau_2 -> stau_1* + tau + tau*     (xintegstau2stau1star)
c       stau_2 -> stau_1 + tau + tau       (xintegstau2stau1)
c       stau_2 -> stau_1 + nu_tau + nu_tau (xintegstau2stau1nn)
      if(flagmulti.eq.1d0) then
        CALL NS_xintegstau2(xintegstau2stau1star,xintegstau2stau1,
     .    xintegstau2stau1nn)
        stau2tot3=xintegstau2stau1star+xintegstau2stau1
     .    +xintegstau2stau1nn
        if (stau2tot3.lt.multilim*stau2tot2)then
           stau2tot3=0d0
        endif
      endif
      stau2tot=stau2tot2+stau2tot3

c-------------------------------------------- c
c ----- stau1/2 branching ratios ------------ c
c-------------------------------------------- c

      if(stau1tot2.ne.0d0)then

       do i=1,5
         brstau1neut(i)=stau1neut(i)/stau1tot2
       enddo
       do i=1,2
         brstau1char(i)=stau1char(i)/stau1tot2
         brstau1hcsn(i)=stau1hcsn(i)/stau1tot2
         brstau1wsn(i) =stau1wsn(i)/stau1tot2
       enddo
       brstau1TAUgra=stau1TAUgra/stau1tot2

      else

       do i=1,5
         brstau1neut(i)=0d0
       enddo
       do i=1,2
         brstau1char(i)=0d0
         brstau1hcsn(i)=0d0
         brstau1wsn(i) =0d0
       enddo
       brstau1TAUgra=0d0

      endif

      if(stau2tot.ne.0d0)then

       do i=1,5
         brstau2neut(i)=stau2neut(i)/stau2tot
       enddo
       do i=1,2
         brstau2char(i)=stau2char(i)/stau2tot
         brstau2hcsn(i)=stau2hcsn(i)/stau2tot
         brstau2wsn(i) =stau2wsn(i)/stau2tot
       enddo
       do I=1,3
         brstau2H(I)=stau2H(I)/stau2tot
       enddo
       do I=1,2
         brstau2A(I)=stau2A(I)/stau2tot
       enddo
       brstau2ztau=stau2ztau/stau2tot
       brstau2TAUgra=stau2TAUgra/stau2tot

      else

       do i=1,5
         brstau2neut(i)=0d0
       enddo
       do i=1,2
         brstau2char(i)=0d0
         brstau2hcsn(i)=0d0
         brstau2wsn(i) =0d0
       enddo
       do I=1,3
         brstau2H(I) =0d0
       enddo
       do I=1,2
         brstau2A(I) =0d0
       enddo
       brstau2ztau=0d0
       brstau2TAUgra=0d0

      endif

      if (stau2tot3.ne.0d0)then

         brstau2stau1star=xintegstau2stau1star/stau2tot
         brstau2stau1    =xintegstau2stau1/stau2tot
         brstau2stau1nn  =xintegstau2stau1nn/stau2tot

      else

         brstau2stau1star=0d0
         brstau2stau1    =0d0
         brstau2stau1nn  =0d0

      endif

c ==================================================================== c
c                       sneutrino_tau 2-body decays                    c
c ==================================================================== c

      gmstau(1) = astau1
      gmstau(2) = astau2
c -------------------------------------------------------------------- c
c  sneutrino_tau1 --> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + neutrino_tau
      do i=1,5
         if(asntau1.gt.amneut(i)) then
            sntauneut(i)=g2s*(antau(1,i)**2+bntau(1,i)**2)*
     .        (asntau1**2-amneut(i)**2)*
     .           NS_lamb(0d0,amneut(i)/asntau1)
     .        /(16*pi*asntau1)
         else
            sntauneut(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  sneutrino_tau1 --> chi+_1/chi+_2 + tau-
      do i=1,2
         if(asntau1.gt.(amchar(i)+mtau)) then
            sntauchar(i)=g2s*((alsnt(1,i)**2+blsnt(1,i)**2)*
     .           (asntau1**2-amchar(i)**2-mtau**2)
     .           -4d0*alsnt(1,i)*blsnt(1,i)*xmchar(i)*mtau)*
     .           NS_lamb(mtau/asntau1,amchar(i)/asntau1)
     .           /(16*pi*asntau1)
         else
            sntauchar(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  sneutrino_tau1 --> H- + stau1/2
c      CALL NS_hcstausntau(gcsntaur)
      do i=1,2
         if(asntau1.gt.(gmstau(i)+cmass)) then
           sntau1hcstau(i)=g2s*mw**2*gcsntaur(1,i)**2*
     .      NS_lamb(gmstau(i)/asntau1,cmass/asntau1)/(16d0*pi*asntau1)
         else
            sntau1hcstau(i)=0d0
         endif
      end do

c -------------------------------------------------------------------- c
c  sneutrino_tau1 --> W- + stau1/2
      do i=1,2
         if(asntau1.gt.(gmstau(i)+mw)) then
            sntau1wstau(i)=g2s/32d0/pi/mw**2*asntau1**3*
     .                gwntau(1,i)**2*
     .                NS_lamb(gmstau(i)/asntau1,mw/asntau1)**3
         else
            sntau1wstau(i)=0d0
         endif
      enddo

c -------------------------------------------------------------------- c
c  sneutrino_tau1 --> neutrino_tau + gravitino

      if(GRFLAG.EQ.1)then
        if (M32.le.asntau1) then
          sntNTgra = astau2**5/(48d0*PI*MPL**2*M32**2)
        else
          sntNTgra = 0d0
        endif
      endif

c -------------------------------------------------------------------- c
c                  SUM 2 BODY DECAYS
c -------------------------------------------------------------------- c

      sntautot2=sntauneut(1)+sntauneut(2)+sntauneut(3)+
     .          sntauneut(4)+sntauneut(5)+sntauchar(1)+sntauchar(2)+
     .          sntau1hcstau(1)+sntau1hcstau(2)+
     .          sntau1wstau(1)+sntau1wstau(2)+sntNTgra

c--------------------------------------------------------------------- c
c ----- sneutrino_tau 3-body decays and 3-body total widths ---------- c
c--------------------------------------------------------------------- c
c Only: sntau -> stau1* + nu_tau + tau  (xintegsntaustau1star)
c       sntau -> stau1 + nu_tau + tau   (xintegsntaustau1)
c       sntau -> stau1 + tau + nu_tau   (xintegsntaustau1nutau)
      if(flagmulti.eq.1d0) then
        CALL NS_xintegsntau(xintegsntaustau1star,xintegsntaustau1,
     .    xintegsntaustau1nutau)
        sntautot3=xintegsntaustau1star+xintegsntaustau1+
     .    xintegsntaustau1nutau
        if (sntautot3.lt.multilim*sntautot2)then
           sntautot3=0d0
        endif
      endif
      sntautot=sntautot2+sntautot3

c---------------------------------------------------- c
c ----- sneutrino_tau branching ratios -------------- c
c---------------------------------------------------- c

      if(sntautot.ne.0d0)then

      do i=1,5
         brsntauneut(i) = sntauneut(i)/sntautot
      enddo
      do i=1,2
         brsntauchar(i)    = sntauchar(i)/sntautot
         brsntau1wstau(i)  = sntau1wstau(i)/sntautot
         brsntau1hcstau(i) = sntau1hcstau(i)/sntautot
      enddo
      brsntNTgra = sntNTgra/sntautot

      else

      do i=1,5
         brsntauneut(i) = 0d0
      enddo
      do i=1,2
         brsntauchar(i)    = 0d0
         brsntau1wstau(i)  = 0d0
         brsntau1hcstau(i) = 0d0
      enddo
      brsntNTgra = 0d0

      endif

      if (sntautot3.ne.0d0)then

         brsntaustau1star  = xintegsntaustau1star/sntautot
         brsntaustau1	   = xintegsntaustau1/sntautot
         brsntaustau1nutau = xintegsntaustau1nutau/sntautot

      else

         brsntaustau1star  = 0d0
         brsntaustau1	   = 0d0
         brsntaustau1nutau = 0d0

      endif

      end

c ==================================================================== c
c                   selectron_L 3-body decays                          c
c ==================================================================== c

      SUBROUTINE NS_xintegsell(xintegsellstau1star,xintegsellstau1,
     .  xintegsellstau1nutau)
*
      INTEGER nx1t,ny1t
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION xmu1,xmu2,xmu3,sum,PI,SQR2
      DOUBLE PRECISION xintegsellstau1star,xintegsellstau1,
     .  xintegsellstau1nutau
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_nx1/nx1t,ny1t
*
      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_sellstau1star,NS_sellstau1,NS_sellstau1nutau

c -------------------------------------------------------------------- c
c ---------------- selectron_L -> stau_1* ---------------------------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=(mtau/ase1)**2
      xmu3=(astau1/ase1)**2

      if(ase1.gt.(astau1+mtau)) then
         CALL NS_integ2(NS_sellstau1star,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsellstau1star=g2s**2/32d0/(2d0*pi)**3*ase1*sum
      else
         xintegsellstau1star=0d0
      endif

c -------------------------------------------------------------------- c
c ---------------- selectron_L -> stau_1 ----------------------------- c
c -------------------------------------------------------------------- c

      if(ase1.gt.(astau1+mtau)) then
         CALL NS_integ2(NS_sellstau1,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsellstau1=g2s**2/32d0/(2d0*pi)**3*ase1*sum
      else
         xintegsellstau1=0d0
      endif

c -------------------------------------------------------------------- c
c --------------- selectron_L -> stau_1, chargino exchange ----------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=0d0
      xmu3=(astau1/ase1)**2

      if(ase1.gt.astau1) then
         CALL NS_integ2(NS_sellstau1nutau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsellstau1nutau=g2s**2/32d0/(2d0*pi)**3*ase1*sum
      else
         xintegsellstau1nutau=0d0
      endif

      end

c ==================================================================== c
c     selectron_L 3-body decay into stau_1*, neutralino exchange       c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_sellstau1star(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
*
      xtau=(mtau/ase1)**2
      do i=1,5
        xneut(i)=(MNEU(i)/ase1)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau1=(astau1/ase1)**2

      NS_sellstau1star=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.ase1).and.
     .       (dabs(mneu(l)).gt.ase1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_sellstau1star=NS_sellstau1star+1d0/dneut(k)/dneut(l)*
     .        ae(1,k)*ae(1,l)*(btau(1,k)*btau(1,l)*
     .        ((1d0-x1)*(1d0-x2)-xstau1+xtau)
     .        +(atau(1,k)*atau(1,l))*
     .        MNEU(k)*MNEU(l)/ase1**2*(x1+x2-1d0+xstau1-xtau)
     .        +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .        +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c     selectron_L 3-body decay into stau_1, neutralino exchange        c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_sellstau1(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau

      xtau=(mtau/ase1)**2
      do i=1,5
        xneut(i)=(MNEU(i)/ase1)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau1=(astau1/ase1)**2

      NS_sellstau1=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.ase1).and.
     .       (dabs(mneu(l)).gt.ase1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_sellstau1=NS_sellstau1+1d0/(dneut(k)*dneut(l))*
     .        ae(1,k)*ae(1,l)*(atau(1,k)*atau(1,l)*
     .        ((1d0-x1)*(1d0-x2)-xstau1+xtau)
     .        +(btau(1,k)*btau(1,l))*
     .        MNEU(k)*MNEU(l)/ase1**2*(x1+x2-1d0+xstau1-xtau)
     .        +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .        +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c     selectron_L 3-body decay into stau_1, chargino exchange          c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_sellstau1nutau(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xstau1,dchar(2),xchar(2)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .   blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt

      do i=1,2
        xchar(i)=(MCH(i)/ase1)**2
        dchar(i)=1d0-x1-xchar(i)
      enddo

      xstau1=(astau1/ase1)**2

      NS_sellstau1nutau=0d0

      do k=1,2
        do l=1,2
          if((dabs(mch(k)).gt.ase1).and.(dabs(mch(l)).gt.ase1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_sellstau1nutau=NS_sellstau1nutau+1d0/(dchar(k)*dchar(l))*
     .        ale(1,k)*ale(1,l)*altau(1,k)*altau(1,l)*
     .        ((1d0-x1)*(1d0-x2)-xstau1)

           endif
        enddo
      enddo

      end

c ==================================================================== c
c                        selectron_R 3-body decays                     c
c ==================================================================== c

      SUBROUTINE NS_xintegselr(xintegselrstau1star,xintegselrstau1,
     .  xintegselrstau1nutau)
*
      INTEGER nx1t,ny1t
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION xmu1,xmu2,xmu3,sum,PI,SQR2
      DOUBLE PRECISION xintegselrstau1star,xintegselrstau1,
     .  xintegselrstau1nutau
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_nx1/nx1t,ny1t
*
      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_selrstau1star,NS_selrstau1,NS_selrstau1nutau

c -------------------------------------------------------------------- c
c ---------------- selectron_R -> stau_1* ---------------------------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=(mtau/ase2)**2
      xmu3=(astau1/ase2)**2

      if(ase2.gt.(astau1+mtau)) then
         CALL NS_integ2(NS_selrstau1star,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegselrstau1star=g2s**2/32d0/(2d0*pi)**3*ase2*sum
      else
         xintegselrstau1star=0d0
      endif

c -------------------------------------------------------------------- c
c ---------------- selectron_R -> stau_1 ----------------------------- c
c -------------------------------------------------------------------- c

      if(ase2.gt.(astau1+mtau)) then
         CALL NS_integ2(NS_selrstau1,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegselrstau1=g2s**2/32d0/(2d0*pi)**3*ase2*sum
      else
         xintegselrstau1=0d0
      endif

c -------------------------------------------------------------------- c
c --------------- selectron_R -> stau_1, chargino exchange ----------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=0d0
      xmu3=(astau1/ase2)**2

      if(ase2.gt.astau1) then
         CALL NS_integ2(NS_selrstau1nutau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegselrstau1nutau=g2s**2/32d0/(2d0*pi)**3*ase2*sum
      else
         xintegselrstau1nutau=0d0
      endif

      end

c ==================================================================== c
c     selectron_R 3-body decay into stau_1*, neutralino exchange       c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_selrstau1star(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
*
      xtau=(mtau/ase2)**2
      do i=1,5
        xneut(i)=(MNEU(i)/ase2)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau=(astau1/ase2)**2

      NS_selrstau1star=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.ase2).and.
     .       (dabs(mneu(l)).gt.ase2)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_selrstau1star=NS_selrstau1star+1d0/(dneut(k)*dneut(l))*
     .        be(2,k)*be(2,l)*(btau(1,k)*btau(1,l)*
     .        dabs(mneu(k)*mneu(l))/ase2**2*(x1+x2-1d0+xstau-xtau)
     .        +atau(1,k)*atau(1,l)*((1d0-x1)*(1d0-x2)-xstau+xtau)
     .        +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .        +btau(1,k)*atau(1,l)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c     selectron_R 3-body decay into stau_1, neutralino exchange        c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_selrstau1(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau

      xtau=(mtau/ase2)**2
      do i=1,5
        xneut(i)=(MNEU(i)/ase2)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau=(astau1/ase2)**2

      NS_selrstau1=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.ase2).and.
     .       (dabs(mneu(l)).gt.ase2)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_selrstau1=NS_selrstau1+1d0/(dneut(k)*dneut(l))*
     .        be(2,k)*be(2,l)*(atau(1,k)*atau(1,l)*
     .        dabs(mneu(k)*mneu(l))/ase2**2*(x1+x2-1d0+xstau-xtau)
     .        +btau(1,k)*btau(1,l)*((1d0-x1)*(1d0-x2)-xstau+xtau)
     .        +btau(1,k)*atau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .        +atau(1,k)*btau(1,l)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c     selectron_R 3-body decay into stau_1, chargino exchange          c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_selrstau1nutau(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xstau1,dchar(2),xchar(2)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .   blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt

      do i=1,2
        xchar(i)=(MCH(i)/ase2)**2
        dchar(i)=1d0-x1-xchar(i)
      enddo

      xstau1=(astau1/ase2)**2

      NS_selrstau1nutau=0d0

      do k=1,2
        do l=1,2
          if((dabs(mch(k)).gt.ase2).and.(dabs(mch(l)).gt.ase2)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_selrstau1nutau=NS_selrstau1nutau+1d0/(dchar(k)*dchar(l))*
     .        ale(2,k)*ale(2,l)*altau(1,k)*altau(1,l)*
     .        ((1d0-x1)*(1d0-x2)-xstau1)

           endif
        enddo
      enddo

      end

c ==================================================================== c
c                          sne 3-body decays                           c
c ==================================================================== c

      SUBROUTINE NS_xintegsne(xintegsnestau1star,xintegsnestau1,
     .  xintegsnestau1nutau)
*
      INTEGER nx1t,ny1t
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION xmu1,xmu2,xmu3,sum,PI,SQR2
      DOUBLE PRECISION xintegsnestau1star,xintegsnestau1,
     .  xintegsnestau1nutau
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_nx1/nx1t,ny1t
*
      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_snestau1star,NS_snestau1,NS_snestau1char

c -------------------------------------------------------------------- c
c --------------------- sne -> stau_1* ------------------------------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=mtau**2/asne1**2
      xmu3=astau1**2/asne1**2

      if(asne1.gt.(astau1+mtau)) then
         CALL NS_integ2(NS_snestau1star,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsnestau1star=g2s**2/32d0/(2d0*pi)**3*asne1*sum
      else
         xintegsnestau1star=0d0
      endif

c -------------------------------------------------------------------- c
c --------------------- sne -> stau_1 -------------------------------- c
c -------------------------------------------------------------------- c

      if(asne1.gt.(astau1+mtau)) then
         CALL NS_integ2(NS_snestau1,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsnestau1=g2s**2/32d0/(2d0*pi)**3*asne1*sum
      else
         xintegsnestau1=0d0
      endif

c -------------------------------------------------------------------- c
c --------------------- sne -> stau_1, chargino exchange ------------- c
c -------------------------------------------------------------------- c

      xmu1=mtau**2/asne1**2
      xmu2=0d0
      xmu3=astau1**2/asne1**2

      if(asne1.gt.astau1) then
         CALL NS_integ2(NS_snestau1char,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsnestau1nutau=g2s**2/32d0/(2d0*pi)**3*asne1*sum
      else
         xintegsnestau1nutau=0d0
      endif

       end

c ==================================================================== c
c         snu_e 3-body decay into stau_1*, neutralino exchange         c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_snestau1star(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
*
      xtau=(mtau/asne1)**2
      do i=1,5
        xneut(i)=(MNEU(i)/asne1)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau1=(astau1/asne1)**2

      NS_snestau1star=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.asne1).and.
     .       (dabs(mneu(l)).gt.asne1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_snestau1star=NS_snestau1star+1d0/dneut(k)/dneut(l)*
     .        ane(1,k)*ane(1,l)*(btau(1,k)*btau(1,l)
     .        *((1d0-x1)*(1d0-x2)-xstau1+xtau)
     .        +(atau(1,k)*atau(1,l))*
     .        (MNEU(k)*MNEU(l))/asne1**2*(x1+x2-1d0+xstau1-xtau)
     .        +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .        +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c        snu_e 3-body decay into stau_1, neutralino exchange           c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_snestau1(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau

      xtau=(mtau/asne1)**2
      do i=1,5
        xneut(i)=(MNEU(i)/asne1)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau1=(astau1/asne1)**2

      NS_snestau1=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.asne1).and.
     .       (dabs(mneu(l)).gt.asne1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_snestau1=NS_snestau1+1d0/(dneut(k)*dneut(l))*
     .        ane(1,k)*ane(1,l)*(atau(1,k)*atau(1,l)
     .        *((1d0-x1)*(1d0-x2)-xstau1+xtau)
     .        +btau(1,k)*btau(1,l)*
     .        (MNEU(k)*MNEU(l))/asne1**2*(x1+x2-1d0+xstau1-xtau)
     .        +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .        +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c        snu_e 3-body decay into stau_1, chargino exchange             c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_snestau1char(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xstau1,dchar(2),xchar(2)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .   blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt

      do i=1,2
        xchar(i)=(MCH(i)/asne1)**2
        dchar(i)=1d0-x1-xchar(i)
      enddo
      xstau1=(astau1/asne1)**2

      NS_snestau1char=0d0

      do k=1,2
        do l=1,2
          if((dabs(mch(k)).gt.asne1).and.(dabs(mch(l)).gt.asne1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_snestau1char=NS_snestau1char+1d0/(dchar(k)*dchar(l))*
     .       (altau(1,k)*altau(1,l))
     .       *alsne(1,k)*alsne(1,l)*
     .       ((1d0-x1)*(1d0-x2)-xstau1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c                   smuon_1 3-body decays                              c
c ==================================================================== c

      SUBROUTINE NS_xintegsmu1(xintegsmu1stau1star,xintegsmu1stau1,
     .  xintegsmu1stau1nutau)
*
      INTEGER nx1t,ny1t
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION xmu1,xmu2,xmu3,sum,PI,SQR2
      DOUBLE PRECISION xintegsmu1stau1star,xintegsmu1stau1,
     .  xintegsmu1stau1nutau
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_nx1/nx1t,ny1t
*
      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_smu1stau1star,NS_smu1stau1,NS_smu1stau1nutau

c -------------------------------------------------------------------- c
c ---------------- smuon_1 -> stau_1* -------------------------------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=(mtau/asmu1)**2
      xmu3=(astau1/asmu1)**2

      if(asmu1.gt.(astau1+mtau)) then
         CALL NS_integ2(NS_smu1stau1star,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsmu1stau1star=g2s**2/32d0/(2d0*pi)**3*asmu1*sum
      else
         xintegsmu1stau1star=0d0
      endif

c -------------------------------------------------------------------- c
c ---------------- smuon_1 -> stau_1 --------------------------------- c
c -------------------------------------------------------------------- c

      if(asmu1.gt.(astau1+mtau)) then
         CALL NS_integ2(NS_smu1stau1,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsmu1stau1=g2s**2/32d0/(2d0*pi)**3*asmu1*sum
      else
         xintegsmu1stau1=0d0
      endif

c -------------------------------------------------------------------- c
c --------------- smuon_1 -> stau_1, chargino exchange --------------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=0d0
      xmu3=(astau1/asmu1)**2

      if(asmu1.gt.astau1) then
         CALL NS_integ2(NS_smu1stau1nutau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsmu1stau1nutau=g2s**2/32d0/(2d0*pi)**3*asmu1*sum
      else
         xintegsmu1stau1nutau=0d0
      endif

      end

c ==================================================================== c
c     smuon_1 3-body decay into stau_1*, neutralino exchange           c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_smu1stau1star(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
*
      xtau=(mtau/asmu1)**2
      do i=1,5
        xneut(i)=(MNEU(i)/asmu1)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau1=(astau1/asmu1)**2

      NS_smu1stau1star=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.asmu1).and.
     .       (dabs(mneu(l)).gt.asmu1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_smu1stau1star=NS_smu1stau1star+1d0/dneut(k)/dneut(l)*
     .        amu(1,k)*amu(1,l)*(btau(1,k)*btau(1,l)*
     .        ((1d0-x1)*(1d0-x2)-xstau1+xtau)
     .        +(atau(1,k)*atau(1,l))*
     .        MNEU(k)*MNEU(l)/asmu1**2*(x1+x2-1d0+xstau1-xtau)
     .        +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .        +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c     smuon_1 3-body decay into stau_1, neutralino exchange            c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_smu1stau1(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau

      xtau=(mtau/asmu1)**2
      do i=1,5
        xneut(i)=(MNEU(i)/asmu1)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau1=(astau1/asmu1)**2

      NS_smu1stau1=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.asmu1).and.
     .       (dabs(mneu(l)).gt.asmu1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_smu1stau1=NS_smu1stau1+1d0/(dneut(k)*dneut(l))*
     .        amu(1,k)*amu(1,l)*(atau(1,k)*atau(1,l)*
     .        ((1d0-x1)*(1d0-x2)-xstau1+xtau)
     .        +(btau(1,k)*btau(1,l))*
     .        MNEU(k)*MNEU(l)/asmu1**2*(x1+x2-1d0+xstau1-xtau)
     .        +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .        +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c     smuon_1 3-body decay into stau_1, chargino exchange              c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_smu1stau1nutau(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xstau1,dchar(2),xchar(2)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .   blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt

      do i=1,2
        xchar(i)=(MCH(i)/asmu1)**2
        dchar(i)=1d0-x1-xchar(i)
      enddo

      xstau1=(astau1/asmu1)**2

      NS_smu1stau1nutau=0d0

      do k=1,2
        do l=1,2
          if((dabs(mch(k)).gt.asmu1).and.(dabs(mch(l)).gt.asmu1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_smu1stau1nutau=NS_smu1stau1nutau+1d0/(dchar(k)*dchar(l))*
     .        almu(1,k)*almu(1,l)*altau(1,k)*altau(1,l)*
     .        ((1d0-x1)*(1d0-x2)-xstau1)

           endif
        enddo
      enddo

      end

c ==================================================================== c
c                        smuon_2 3-body decays                         c
c ==================================================================== c

      SUBROUTINE NS_xintegsmu2(xintegsmu2stau1star,xintegsmu2stau1,
     .  xintegsmu2stau1nutau)
*
      INTEGER nx1t,ny1t
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION xmu1,xmu2,xmu3,sum,PI,SQR2
      DOUBLE PRECISION xintegsmu2stau1star,xintegsmu2stau1,
     .  xintegsmu2stau1nutau
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_nx1/nx1t,ny1t
*
      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_smu2stau1star,NS_smu2stau1,NS_smu2stau1nutau

c -------------------------------------------------------------------- c
c ---------------- smuon_2 -> stau_1* -------------------------------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=(mtau/asmu2)**2
      xmu3=(astau1/asmu2)**2

      if(asmu2.gt.(astau1+mtau)) then
         CALL NS_integ2(NS_smu2stau1star,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsmu2stau1star=g2s**2/32d0/(2d0*pi)**3*asmu2*sum
      else
         xintegsmu2stau1star=0d0
      endif

c -------------------------------------------------------------------- c
c ---------------- smuon_2 -> stau_1 --------------------------------- c
c -------------------------------------------------------------------- c

      if(asmu2.gt.(astau1+mtau)) then
         CALL NS_integ2(NS_smu2stau1,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsmu2stau1=g2s**2/32d0/(2d0*pi)**3*asmu2*sum
      else
         xintegsmu2stau1=0d0
      endif

c -------------------------------------------------------------------- c
c --------------- smuon_2 -> stau_1, chargino exchange --------------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=0d0
      xmu3=astau1**2/asmu2**2

      if(asmu2.gt.astau1) then
         CALL NS_integ2(NS_smu2stau1nutau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsmu2stau1nutau=g2s**2/32d0/(2d0*pi)**3*asmu2*sum
      else
         xintegsmu2stau1nutau=0d0
      endif

      end

c ==================================================================== c
c     smuon_2 3-body decay into stau_1*, neutralino exchange           c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_smu2stau1star(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
*
      xtau=(mtau/asmu2)**2
      do i=1,5
        xneut(i)=(MNEU(i)/asmu2)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau=(astau1/asmu2)**2

      NS_smu2stau1star=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.asmu2).and.
     .       (dabs(mneu(l)).gt.asmu2)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_smu2stau1star=NS_smu2stau1star+1d0/(dneut(k)*dneut(l))*
     .        bmu(2,k)*bmu(2,l)*(btau(1,k)*btau(1,l)*
     .        dabs(mneu(k)*mneu(l))/asmu2**2*(x1+x2-1d0+xstau-xtau)
     .        +atau(1,k)*atau(1,l)*((1d0-x1)*(1d0-x2)-xstau+xtau)
     .        +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .        +btau(1,k)*atau(1,l)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c     smuon_2 3-body decay into stau_1, neutralino exchange            c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_smu2stau1(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau

      xtau=(mtau/asmu2)**2
      do i=1,5
        xneut(i)=(MNEU(i)/asmu2)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau=(astau1/asmu2)**2

      NS_smu2stau1=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.asmu2).and.
     .       (dabs(mneu(l)).gt.asmu2)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_smu2stau1=NS_smu2stau1+1d0/(dneut(k)*dneut(l))*
     .        bmu(2,k)*bmu(2,l)*(atau(1,k)*atau(1,l)*
     .        dabs(mneu(k)*mneu(l))/asmu2**2*(x1+x2-1d0+xstau-xtau)
     .        +btau(1,k)*btau(1,l)*((1d0-x1)*(1d0-x2)-xstau+xtau)
     .        +btau(1,k)*atau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .        +atau(1,k)*btau(1,l)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c     smuon_2 3-body decay into stau_1, chargino exchange              c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_smu2stau1nutau(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xstau1,dchar(2),xchar(2)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .   blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt

      do i=1,2
        xchar(i)=(MCH(i)/asmu2)**2
        dchar(i)=1d0-x1-xchar(i)
      enddo

      xstau1=(astau1/asmu2)**2

      NS_smu2stau1nutau=0d0

      do k=1,2
        do l=1,2
          if((dabs(mch(k)).gt.asmu2).and.(dabs(mch(l)).gt.asmu2)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_smu2stau1nutau=NS_smu2stau1nutau+1d0/(dchar(k)*dchar(l))*
     .        almu(2,k)*almu(2,l)*altau(1,k)*altau(1,l)*
     .        ((1d0-x1)*(1d0-x2)-xstau1)

           endif
        enddo
      enddo

      end

c ==================================================================== c
c                          snmu 3-body decays                          c
c ==================================================================== c

      SUBROUTINE NS_xintegsnmu(xintegsnmustau1star,xintegsnmustau1,
     .  xintegsnmustau1nutau)
*
      INTEGER nx1t,ny1t
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION xmu1,xmu2,xmu3,sum,PI,SQR2
      DOUBLE PRECISION xintegsnmustau1star,xintegsnmustau1,
     .  xintegsnmustau1nutau
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_nx1/nx1t,ny1t
*
      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_snmustau1star,NS_snmustau1,NS_snmustau1char

c -------------------------------------------------------------------- c
c --------------------- snmu -> stau_1* ------------------------------ c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=mtau**2/asnmu1**2
      xmu3=astau1**2/asnmu1**2

      if(asnmu1.gt.(astau1+mtau)) then
         CALL NS_integ2(NS_snmustau1star,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsnmustau1star=g2s**2/32d0/(2d0*pi)**3*asnmu1*sum
      else
         xintegsnmustau1star=0d0
      endif

c -------------------------------------------------------------------- c
c --------------------- snmu -> stau_1 ------------------------------- c
c -------------------------------------------------------------------- c

      if(asnmu1.gt.(astau1+mtau)) then
         CALL NS_integ2(NS_snmustau1,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsnmustau1=g2s**2/32d0/(2d0*pi)**3*asnmu1*sum
      else
         xintegsnmustau1=0d0
      endif

c -------------------------------------------------------------------- c
c --------------------- snmu -> stau_1, chargino exchange ------------ c
c -------------------------------------------------------------------- c

      xmu1=mtau**2/asnmu1**2
      xmu2=0d0
      xmu3=astau1**2/asnmu1**2

      if(asnmu1.gt.astau1) then
         CALL NS_integ2(NS_snmustau1char,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsnmustau1nutau=g2s**2/32d0/(2d0*pi)**3*asnmu1*sum
      else
         xintegsnmustau1nutau=0d0
      endif

       end

c ==================================================================== c
c         snu_mu 3-body decay into stau_1*, neutralino exchange        c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_snmustau1star(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
*
      xtau=(mtau/asnmu1)**2
      do i=1,5
        xneut(i)=(MNEU(i)/asnmu1)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau1=(astau1/asnmu1)**2

      NS_snmustau1star=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.asnmu1).and.
     .       (dabs(mneu(l)).gt.asnmu1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_snmustau1star=NS_snmustau1star+1d0/dneut(k)/dneut(l)*
     .        anmu(1,k)*anmu(1,l)*(btau(1,k)*btau(1,l)
     .        *((1d0-x1)*(1d0-x2)-xstau1+xtau)
     .        +(atau(1,k)*atau(1,l))*
     .        (MNEU(k)*MNEU(l))/asnmu1**2*(x1+x2-1d0+xstau1-xtau)
     .        +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .        +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c        snu_mu 3-body decay into stau_1, neutralino exchange          c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_snmustau1(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau

      xtau=(mtau/asnmu1)**2
      do i=1,5
        xneut(i)=(MNEU(i)/asnmu1)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau1=(astau1/asnmu1)**2

      NS_snmustau1=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.asnmu1).and.
     .       (dabs(mneu(l)).gt.asnmu1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_snmustau1=NS_snmustau1+1d0/(dneut(k)*dneut(l))*
     .        anmu(1,k)*anmu(1,l)*(atau(1,k)*atau(1,l)
     .        *((1d0-x1)*(1d0-x2)-xstau1+xtau)
     .        +btau(1,k)*btau(1,l)*
     .        (MNEU(k)*MNEU(l))/asnmu1**2*(x1+x2-1d0+xstau1-xtau)
     .        +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .        +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c        snu_mu 3-body decay into stau_1, chargino exchange            c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_snmustau1char(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xstau1,dchar(2),xchar(2)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .   blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt

      do i=1,2
        xchar(i)=(MCH(i)/asnmu1)**2
        dchar(i)=1d0-x1-xchar(i)
      enddo
      xstau1=(astau1/asnmu1)**2

      NS_snmustau1char=0d0

      do k=1,2
        do l=1,2
          if((dabs(mch(k)).gt.asnmu1).and.(dabs(mch(l)).gt.asnmu1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_snmustau1char=NS_snmustau1char+1d0/(dchar(k)*dchar(l))*
     .       (altau(1,k)*altau(1,l))
     .       *alsnm(1,k)*alsnm(1,l)*
     .       ((1d0-x1)*(1d0-x2)-xstau1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c                        stau_2 3-body decays                          c
c ==================================================================== c

      SUBROUTINE NS_xintegstau2(xintegstau2stau1star,xintegstau2stau1,
     .  xintegstau2stau1nn)
*
      INTEGER nx1t,ny1t
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION xmu1,xmu2,xmu3,sum,PI,SQR2
      DOUBLE PRECISION xintegstau2stau1star,xintegstau2stau1,
     .  xintegstau2stau1nn
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_nx1/nx1t,ny1t
*
      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_stau2stau1star,NS_stau2stau1,NS_stau2stau1nn

c -------------------------------------------------------------------- c
c --------------------- stau_2 -> stau_1* ---------------------------- c
c -------------------------------------------------------------------- c
      xmu1=mtau**2/astau2**2
      xmu2=mtau**2/astau2**2
      xmu3=astau1**2/astau2**2

      if(astau2.gt.(astau1+2d0*mtau)) then
         CALL NS_integ2(NS_stau2stau1star,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegstau2stau1star=g2s**2/32d0/(2d0*pi)**3*astau2*sum
      else
         xintegstau2stau1star=0d0
      endif

c -------------------------------------------------------------------- c
c --------------------- stau_2 -> stau_1 ----------------------------- c
c -------------------------------------------------------------------- c

      if(astau2.gt.(astau1+2d0*mtau)) then
         CALL NS_integ2(NS_stau2stau1,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegstau2stau1=g2s**2/32d0/(2d0*pi)**3*astau2*sum
      else
         xintegstau2stau1=0d0
      endif

c -------------------------------------------------------------------- c
c --------------------- stau_2 -> stau_1, chargino exchange ---------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=0d0
      xmu3=astau1**2/astau2**2

      if(astau2.gt.astau1) then
         CALL NS_integ2(NS_stau2stau1nn,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegstau2stau1nn=g2s**2/32d0/(2d0*pi)**3*astau2*sum
      else
         xintegstau2stau1nn=0d0
      endif

      end

c ==================================================================== c
c        stau_2 3-body decay into stau_1*, neutralino exchange         c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_stau2stau1star(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau
*
      xtau=(mtau/astau2)**2
      do i=1,5
        xneut(i)=(MNEU(i)/astau2)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau1=(astau1/astau2)**2

      NS_stau2stau1star=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.astau2-mtau).and.
     .       (dabs(mneu(l)).gt.astau2-mtau)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_stau2stau1star=NS_stau2stau1star+1d0/dneut(k)/dneut(l)*(
     .        (btau(1,k)*btau(1,l)*atau(2,k)*atau(2,l)+
     .        atau(1,k)*atau(1,l)*btau(2,k)*btau(2,l))*
     .        ((1d0-x1)*(1d0-x2)-xstau1+xtau*(x1-x2+xstau1
     .        -2d0*xtau)+xtau)
     .        +(btau(1,k)*btau(1,l)*btau(2,k)*btau(2,l)+
     .        atau(1,k)*atau(1,l)*atau(2,k)*atau(2,l))*
     .        MNEU(k)*MNEU(l)/astau2**2*(x1+x2-1d0+xstau1
     .        -2d0*xtau)
     .        +(btau(1,k)*atau(1,l)*atau(2,k)*btau(2,l)+
     .        atau(1,k)*btau(1,l)*btau(2,k)*atau(2,l))*
     .        2d0*xtau*(-1d0+x1-xtau)
     .        +dsqrt(xtau)*(x1-2d0*xtau)*(MNEU(k)/astau2*
     .        (btau(1,k)*atau(1,l)*btau(2,k)*btau(2,l)+
     .        atau(1,k)*btau(1,l)*atau(2,k)*atau(2,l))
     .        +MNEU(l)/astau2*
     .        (btau(1,k)*atau(1,l)*atau(2,k)*atau(2,l)+
     .        atau(1,k)*btau(1,l)*btau(2,k)*btau(2,l)) )
     .        +dsqrt(xtau)*(-1d0+x1+xstau1-2d0*xtau)*
     .        (MNEU(k)/astau2*
     .        (btau(1,k)*btau(1,l)*btau(2,k)*atau(2,l)+
     .         atau(1,k)*atau(1,l)*atau(2,k)*btau(2,l))
     .         +MNEU(l)/astau2*
     .         (btau(1,k)*btau(1,l)*atau(2,k)*btau(2,l)+
     .         atau(1,k)*atau(1,l)*btau(2,k)*atau(2,l)) )
     .         +(btau(1,l)*atau(1,k)*atau(2,k)*btau(2,l)+
     .         atau(1,l)*btau(1,k)*btau(2,k)*atau(2,l))*
     .         (-2d0)*xtau*MNEU(k)*MNEU(l)/astau2**2 )
          endif
        enddo
      enddo

      end

c ==================================================================== c
c        stau_2 3-body decay into stau_1, neutralino exchange          c
c ==================================================================== c
      DOUBLE PRECISION FUNCTION NS_stau2stau1(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau

      xtau=(mtau/astau2)**2
      do i=1,5
        xneut(i)=(MNEU(i)/astau2)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau1=(astau1/astau2)**2

      NS_stau2stau1=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.astau2-mtau).and.
     .       (dabs(mneu(l)).gt.astau2-mtau)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_stau2stau1=NS_stau2stau1+1d0/(dneut(k)*dneut(l))*(
     .        (atau(1,k)*atau(1,l)*atau(2,k)*atau(2,l)+
     .        btau(1,k)*btau(1,l)*btau(2,k)*btau(2,l))*
     .        ((1d0-x1)*(1d0-x2)-xstau1+xtau*(x1-x2+xstau1
     .        -2d0*xtau)+xtau)
     .        +(atau(1,k)*atau(1,l)*btau(2,k)*btau(2,l)+
     .        btau(1,k)*btau(1,l)*atau(2,k)*atau(2,l))*
     .        MNEU(k)*MNEU(l)/astau2**2*(x1+x2-1d0+xstau1
     .        -2d0*xtau)
     .        +(atau(1,k)*btau(1,l)*atau(2,k)*btau(2,l)+
     .        btau(1,k)*atau(1,l)*btau(2,k)*atau(2,l))*
     .        2d0*xtau*(-1d0+x1-xtau)
     .        +dsqrt(xtau)*(x1-2d0*xtau)*(MNEU(k)/astau2*
     .        (atau(1,k)*btau(1,l)*btau(2,k)*btau(2,l)+
     .        btau(1,k)*atau(1,l)*atau(2,k)*atau(2,l))
     .        +MNEU(l)/astau2*
     .        (atau(1,k)*btau(1,l)*atau(2,k)*atau(2,l)+
     .        btau(1,k)*atau(1,l)*btau(2,k)*btau(2,l)) )
     .        +dsqrt(xtau)*(-1d0+x1+xstau1-2d0*xtau)*
     .        (MNEU(k)/astau2*
     .        (atau(1,k)*atau(1,l)*btau(2,k)*atau(2,l)+
     .        btau(1,k)*btau(1,l)*atau(2,k)*btau(2,l))
     .        +MNEU(l)/astau2*
     .        (atau(1,k)*atau(1,l)*atau(2,k)*btau(2,l)+
     .        btau(1,k)*btau(1,l)*btau(2,k)*atau(2,l)) )
     .        +(atau(1,l)*btau(1,k)*atau(2,k)*btau(2,l)+
     .        btau(1,l)*atau(1,k)*btau(2,k)*atau(2,l))*
     .        (-2d0)*xtau*MNEU(k)*MNEU(l)/astau2**2 )
          endif
        enddo
      enddo

      end
c ==================================================================== c
c        stau_2 3-body decay into stau_1, chargino exchange            c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_stau2stau1nn(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xstau1,dchar(2),xchar(2)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .   blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt

      do i=1,2
        xchar(i)=(MCH(i)/astau2)**2
        dchar(i)=1d0-x1-xchar(i)
      enddo

      xstau1=(astau1/astau2)**2

      NS_stau2stau1nn=0d0

      do k=1,2
        do l=1,2
          if((dabs(mch(k)).gt.astau2).and.(dabs(mch(l)).gt.astau2)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_stau2stau1nn=NS_stau2stau1nn+1d0/(dchar(k)*dchar(l))*
     .        altau(2,k)*altau(2,l)*altau(1,k)*altau(1,l)*
     .        ((1d0-x1)*(1d0-x2)-xstau1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c                          sntau 3-body decays                         c
c ==================================================================== c

      SUBROUTINE NS_xintegsntau(xintegsntaustau1star,xintegsntaustau1,
     .  xintegsntaustau1nutau)
*
      INTEGER nx1t,ny1t
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION xmu1,xmu2,xmu3,sum,PI,SQR2
      DOUBLE PRECISION xintegsntaustau1star,xintegsntaustau1,
     .  xintegsntaustau1nutau
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_pi/PI,SQR2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_nx1/nx1t,ny1t

      EXTERNAL NS_ay,NS_by,NS_ax,NS_bx
      EXTERNAL NS_sntaustau1star,NS_sntaustau1,NS_sntaustau1nutau

c -------------------------------------------------------------------- c
c --------------------- sntau -> stau_1* ----------------------------- c
c -------------------------------------------------------------------- c

      xmu1=0d0
      xmu2=mtau**2/asntau1**2
      xmu3=astau1**2/asntau1**2

      if(asntau1.gt.(astau1+mtau)) then
        CALL NS_integ2(NS_sntaustau1star,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsntaustau1star=g2s**2/32d0/(2d0*pi)**3*asntau1*sum
      else
         xintegsntaustau1star=0d0
      endif

c -------------------------------------------------------------------- c
c --------------------- sntau -> stau_1 ------------------------------ c
c -------------------------------------------------------------------- c

      if(asntau1.gt.(astau1+mtau)) then
         CALL NS_integ2(NS_sntaustau1,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsntaustau1=g2s**2/32d0/(2d0*pi)**3*asntau1*sum
      else
         xintegsntaustau1=0d0
      endif

c -------------------------------------------------------------------- c
c --------------------- sntau -> stau_1, chargino exchange ----------- c
c -------------------------------------------------------------------- c

      xmu1=mtau**2/asntau1**2
      xmu2=0d0
      xmu3=astau1**2/asntau1**2

      if(asntau1.gt.mtau+astau1) then
         CALL NS_integ2(NS_sntaustau1nutau,NS_ax,NS_bx,NS_ay,NS_by,xmu1,
     .        xmu2,xmu3,nx1t,ny1t,sum)
         xintegsntaustau1nutau=g2s**2/32d0/(2d0*pi)**3*asntau1*sum
      else
         xintegsntaustau1nutau=0d0
      endif
c -------------------------------------------------------------------- c

      end

c ==================================================================== c
c        sntau 3-body decay into stau_1*, neutralino exchange         c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_sntaustau1star(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau

      xtau=(mtau/asntau1)**2
      do i=1,5
        xneut(i)=(MNEU(i)/asntau1)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau1=(astau1/asntau1)**2

      NS_sntaustau1star=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.asntau1).and.
     .       (dabs(mneu(l)).gt.asntau1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_sntaustau1star=NS_sntaustau1star+1d0/dneut(k)/dneut(l)*
     .        antau(1,k)*antau(1,l)*(btau(1,k)*btau(1,l)
     .        *((1d0-x1)*(1d0-x2)-xstau1+xtau)
     .        +(atau(1,k)*atau(1,l))*
     .        (MNEU(k)*MNEU(l))/asntau1**2*(x1+x2-1d0+xstau1-xtau)
     .        +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .        +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c        sntau 3-body decay into stau_1, neutralino exchange          c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_sntaustau1(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau1,dneut(5),xneut(5)
      DOUBLE PRECISION ae(2,5),be(2,5),amu(2,5),bmu(2,5),atau(2,5),
     .     btau(2,5),ane(2,5),bne(2,5),anmu(2,5),bnmu(2,5),
     .     antau(2,5),bntau(2,5)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup8/ae,be,amu,bmu,atau,btau,ane,bne,anmu,bnmu,antau,
     .bntau

      xtau=(mtau/asntau1)**2
      do i=1,5
        xneut(i)=(MNEU(i)/asntau1)**2
        dneut(i)=1d0-x1+xtau-xneut(i)
      enddo
      xstau1=(astau1/asntau1)**2

      NS_sntaustau1=0d0

      do k=1,5
        do l=1,5
          if((dabs(mneu(k)).gt.asntau1).and.
     .       (dabs(mneu(l)).gt.asntau1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_sntaustau1=NS_sntaustau1+1d0/(dneut(k)*dneut(l))*
     .        antau(1,k)*antau(1,l)*(atau(1,k)*atau(1,l)
     .        *((1d0-x1)*(1d0-x2)-xstau1+xtau)
     .        +(btau(1,k)*btau(1,l))*
     .        (MNEU(k)*MNEU(l))/asntau1**2*(x1+x2-1d0+xstau1-xtau)
     .        +atau(1,k)*btau(1,l)*dsqrt(xneut(l)*xtau)*x1
     .        +atau(1,l)*btau(1,k)*dsqrt(xneut(k)*xtau)*x1)
          endif
        enddo
      enddo

      end

c ==================================================================== c
c        sntau 3-body decay into stau_1, chargino exchange             c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_sntaustau1nutau(x1,x2)
*
      INTEGER i,k,l
      DOUBLE PRECISION x1,x2,xtau,xstau1,dchar(2),xchar(2)
      DOUBLE PRECISION ale(2,2),almu(2,2),altau(2,2),alsne(2,2),
     .   blsne(2,2),alsnm(2,2),blsnm(2,2),alsnt(2,2),blsnt(2,2)
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
*
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_coup5/ale,almu,altau,alsne,blsne,alsnm,blsnm,alsnt,blsnt

      xtau=(mtau/asntau1)**2
      do i=1,2
        xchar(i)=(MCH(i)/asntau1)**2
        dchar(i)=1d0-x1+xtau-xchar(i)
      enddo
      xstau1=(astau1/asntau1)**2

      NS_sntaustau1nutau=0d0

      do k=1,2
        do l=1,2
          if((dabs(mch(k)).gt.asntau1).and.
     .       (dabs(mch(l)).gt.asntau1)) then
c Otherwise: 2-body decays are possible
c -> singularity from on-shell intermediate state
            NS_sntaustau1nutau=NS_sntaustau1nutau+1d0/(dchar(k)*
     .        dchar(l))*alsnt(1,k)*alsnt(1,l)*altau(1,k)*altau(1,l)*
     .        ((1d0-x1)*(1d0-x2)-xstau1+xtau*(xstau1+x1-x2-xtau))
          endif
        enddo
      enddo

      end

