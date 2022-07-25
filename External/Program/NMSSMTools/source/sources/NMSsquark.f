      SUBROUTINE NS_SQUARKS

************************************************************************
*
*     This subroutine computes the u, d, c & s squark decays
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

      INTEGER I,GRFLAG

      DOUBLE PRECISION NS_lamb
      DOUBLE PRECISION amuv,lamv
      DOUBLE PRECISION PI,SQR2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION flagmulti,flagqcd,flagloop
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      DOUBLE PRECISION SCALb,SCALt,scaltau,gs2
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION alup(2,2),aldo(2,2)
      DOUBLE PRECISION aup(2,5),bup(2,5),ado(2,5),bdo(2,5)
      DOUBLE PRECISION thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .         cmm,smm,cum,sum,cdm,sdm,cem,sem,cnm,snm
      DOUBLE PRECISION scala,alp,ca,cf,amsq,rval
      DOUBLE PRECISION amurefer
      DOUBLE PRECISION suplneutup(5),suprneutup(5),suplchardow(2),
     .         suprchardow(2),qcdsuplneutup(5),qcdsuprneutup(5),
     .         qcdsuplchardow(2),qcdsuprchardow(2)
      DOUBLE PRECISION sdowlneutdow(5),sdowlcharup(2),sdowrneutdow(5),
     .         sdowrcharup(2),qcdsdowlneutdow(5),qcdsdowlcharup(2),
     .         qcdsdowrneutdow(5),qcdsdowrcharup(2)
      DOUBLE PRECISION sulUgra,surUgra,sdlDgra,sdrDgra
      DOUBLE PRECISION supltot2,suprtot2,sdowltot2,sdowrtot2
      DOUBLE PRECISION brsuplnup(5),brsuplcdow(2),brsuplglui,
     .         brsuprnup(5),brsuprcdow(2),brsuprglui,
     .         brsdowlndow(5),brsdowlchup(2),brsdowlglui,
     .         brsdowrndow(5),brsdowrchup(2),brsdowrglui
      DOUBLE PRECISION supltot2lo,suprtot2lo,supltot2nlo,suprtot2nlo,
     .         sdowltot2lo,sdowrtot2lo,sdowltot2nlo,sdowrtot2nlo
      DOUBLE PRECISION suplglui,qcdsuplglui,qcdsuprglui,
     .         sdowlglui,sdowrglui,qcdsdowlglui,qcdsdowrglui,suprglui
      DOUBLE PRECISION NS_ftotqcd,NS_gama,NS_gamfcap,NS_gamf,
     .         NS_gamglui2,NS_gamrendec,resum
      DOUBLE PRECISION brcharWgra(2),brcharHCgra(2),brneutGAMgra(5),
     .         brneutZgra(5),brneutHgra(5,3),brneutAgra(5,2),
     .         brgluiGLUgra,brselEgra,brserEgra,brsmu1MUgra,
     .         brsmu2MUgra,brstau1TAUgra,brstau2TAUgra,brsneNEgra,
     .         brsnmNMgra,brsntNTgra,brsulUgra,brsurUgra,brsdlDgra,
     .         brsdrDgra,brst1Tgra,brst2Tgra,brsb1Bgra,brsb2Bgra
      DOUBLE PRECISION KNG(5),KNZ(5),KNH(5,3),KNA(5,2),KCW(2),KCH(2)
      DOUBLE PRECISION M32,CGR,MPL

      COMMON/SQUARK_2GAMMA/suplneutup,suprneutup,suplchardow,
     .         suprchardow,qcdsuplneutup,qcdsuprneutup,
     .         qcdsuplchardow,qcdsuprchardow,
     .         sdowlneutdow,sdowlcharup,sdowrneutdow,
     .         sdowrcharup,qcdsdowlneutdow,qcdsdowlcharup,
     .         qcdsdowrneutdow,qcdsdowrcharup
      COMMON/SQUARK_WIDTH/supltot2,suprtot2,sdowltot2,sdowrtot2
      COMMON/SQUARK_BR_2BD/brsuplnup,brsuplcdow,brsuplglui,
     .         brsuprnup,brsuprcdow,brsuprglui,
     .         brsdowlndow,brsdowlchup,brsdowlglui,
     .         brsdowrndow,brsdowrchup,brsdowrglui
      COMMON/GRAVITINO/brcharWgra,brcharHCgra,brneutGAMgra,
     .         brneutZgra,brneutHgra,brneutAgra,
     .         brgluiGLUgra,brselEgra,brserEgra,brsmu1MUgra,
     .         brsmu2MUgra,brstau1TAUgra,brstau2TAUgra,brsneNEgra,
     .         brsnmNMgra,brsntNTgra,brsulUgra,brsurUgra,brsdlDgra,
     .         brsdrDgra,brst1Tgra,brst2Tgra,brsb1Bgra,brsb2Bgra
      COMMON/GRAVCOUP/KNG,KNZ,KNH,KNA,KCW,KCH
      COMMON/M32/M32,CGR,MPL,GRFLAG
      COMMON/NS_pi/PI,SQR2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/NS_FLAGS/flagmulti,flagqcd,flagloop
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu
      COMMON/NS_sfmixang/thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .         cmm,smm,cum,sum,cdm,sdm,cem,sem,cnm,snm
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_coup7/alup,aldo
      COMMON/NS_coup10/aup,bup,ado,bdo
      COMMON/NS_refscale/amurefer
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_qcdscales/amuv,lamv

      EXTERNAL NS_lamb
      EXTERNAL NS_ftotqcd
      EXTERNAL NS_gama,NS_gamfcap,NS_gamf,NS_gamglui2,NS_gamrendec
      EXTERNAL resum

*  Initialization

      supltot2 = 0d0
      suprtot2 = 0d0
      supltot2lo = 0d0
      suprtot2lo = 0d0
      supltot2nlo = 0d0
      suprtot2nlo = 0d0
      do i=1,5
         suplneutup(i) = 0d0
         suprneutup(i) = 0d0
         qcdsuplneutup(i) = 0d0
         qcdsuprneutup(i) = 0d0
      enddo
      do i=1,2
         suplchardow(i) = 0d0
         suprchardow(i) = 0d0
         qcdsuplchardow(i) = 0d0
         qcdsuprchardow(i) = 0d0
      enddo
      suplglui = 0d0
      suprglui = 0d0
      qcdsuplglui = 0d0
      qcdsuprglui = 0d0

      sdowltot2 = 0d0
      sdowrtot2 = 0d0
      sdowltot2lo = 0d0
      sdowrtot2lo = 0d0
      sdowltot2nlo = 0d0
      sdowrtot2nlo = 0d0
      do i=1,5
         sdowlneutdow(i) = 0d0
         sdowrneutdow(i) = 0d0
         qcdsdowlneutdow(i) = 0d0
         qcdsdowrneutdow(i) = 0d0
      enddo
      do i=1,2
         sdowlcharup(i) = 0d0
         sdowrcharup(i) = 0d0
         qcdsdowlcharup(i) = 0d0
         qcdsdowrcharup(i) = 0d0
      enddo
      sdowlglui = 0d0
      sdowrglui = 0d0
      qcdsdowlglui = 0d0
      qcdsdowrglui = 0d0

      sulUgra = 0d0
      surUgra = 0d0
      sdlDgra = 0d0
      sdrDgra = 0d0

c -------------------------------------------------------------------- c
c  supl -> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + up
      do i=1,5
         if(amneut(i).le.asup1) then
            suplneutup(i)=g2s*(aup(1,i)**2+bup(1,i)**2)*(asup1**2
     .           -amneut(i)**2)*NS_lamb(0d0,amneut(i)/asup1)
     .           /(16*pi*asup1)
         else
            suplneutup(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,5
         if(amneut(i).le.asup1) then
            qcdsuplneutup(i)=4d0/3d0*gs2/(4d0*pi)/pi*
     .           suplneutup(i)*
     .           NS_ftotqcd(amneut(i)**2/asup1**2,mgluino**2/asup1**2)
         else
            qcdsuplneutup(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  supl -> chi+_1/chi+_2 + down
      do i=1,2
         if (amchar(i).le.asup1) then
            suplchardow(i)=g2s*alup(1,i)**2*(asup1**2-amchar(i)**2)*
     .           NS_lamb(0d0,amchar(i)/asup1)/(16*pi*asup1)
         else
            suplchardow(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,2
         if(amchar(i).le.asup1) then
            qcdsuplchardow(i)=4d0/3d0*gs2/(4d0*pi)/pi*
     .           suplchardow(i)*
     .           NS_ftotqcd(amchar(i)**2/asup1**2,mgluino**2/asup1**2)
         else
            qcdsuplchardow(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  supr -> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + up
      do i=1,5
         if(amneut(i).le.asup2) then
            suprneutup(i)=g2s*(aup(2,i)**2+bup(2,i)**2)*(asup2**2
     .           -amneut(i)**2)*NS_lamb(0d0,amneut(i)/asup2)
     .           /(16*pi*asup2)
         else
            suprneutup(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,5
         if(amneut(i).le.asup2) then
            qcdsuprneutup(i)=4d0/3d0*gs2/(4d0*pi)/pi*
     .           suprneutup(i)*
     .           NS_ftotqcd(amneut(i)**2/asup2**2,mgluino**2/asup2**2)
         else
            qcdsuprneutup(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  supr -> chi+_1/chi+_2 + down
      do i=1,2
         if (amchar(i).le.asup2) then
            suprchardow(i)=g2s*alup(2,i)**2*(asup2**2-amchar(i)**2)*
     .           NS_lamb(0d0,amchar(i)/asup2)/(16*pi*asup2)
         else
            suprchardow(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,2
         if(amchar(i).le.asup2) then
            qcdsuprchardow(i)=4d0/3d0*gs2/(4d0*pi)/pi*
     .           suprchardow(i)*
     .           NS_ftotqcd(amchar(i)**2/asup2**2,mgluino**2/asup2**2)
         else
            qcdsuprchardow(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  supl --> gluino + up
      if(asup1.gt.mgluino) then
         suplglui = 8d0*gs2*(asup1**2-mgluino**2)*
     .              NS_lamb(0d0,mgluino/asup1)/(16d0*pi*asup1)/3d0
      else
         suplglui = 0d0
      endif
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      if(asup1.gt.mgluino) then
         scala = amurefer
         alp   = gs2/(4d0*pi)
         ca    = 3d0
         cf    = 4d0/3d0
         amsq  = 2d0*(asup1+asup2+asdown1+asdown2)/8d0
         rval  = mgluino**2/amsq**2

         qcdsuplglui = suplglui*alp/pi*( ca*NS_gama(rval) +
     .        cf*NS_gamfcap(rval) + 4d0*NS_gamf(rval) +
     .        2d0*pi**2/mgluino**2*
     .        NS_gamglui2(ast1,ast2,amt,thet,asb1,asb2,amb,theb,mgluino,
     .        1,amsq) +
     .        NS_gamrendec(amsq,ast1,ast2,amt,asb1,asb2,mgluino,scala) )
      else
         qcdsuplglui = 0d0
      endif
      endif
c -------------------------------------------------------------------- c
c  supr --> gluino + up

      if(asup2.gt.mgluino) then
         suprglui = 8d0*gs2*(asup2**2-mgluino**2)*
     .              NS_lamb(0d0,mgluino/asup2)/(16d0*pi*asup2)/3d0
      else
         suprglui = 0d0
      endif

c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      if(asup2.gt.mgluino) then
         scala = amurefer
         alp   = gs2/(4d0*pi)
         ca    = 3d0
         cf    = 4d0/3d0
         amsq  = 2d0*(asup1+asup2+asdown1+asdown2)/8d0
         rval  = mgluino**2/amsq**2

         qcdsuprglui = suprglui*alp/pi*( ca*NS_gama(rval) +
     .        cf*NS_gamfcap(rval) + 4d0*NS_gamf(rval) +
     .        2d0*pi**2/mgluino**2*
     .        NS_gamglui2(ast2,ast1,amt,thet,asb2,asb1,amb,theb,mgluino,
     .        2,amsq) +
     .        NS_gamrendec(amsq,ast1,ast2,amt,asb1,asb2,mgluino,scala) )
      else
         qcdsuprglui = 0d0
      endif
      endif
c -------------------------------------------------------------------- c
c  supL --> up + gravitino

      IF(GRFLAG.EQ.1)THEN
        if (M32.le.asup1) then
          sulUgra = asup1**5/(48d0*PI*MPL**2*M32**2)
        else
          sulUgra = 0d0
        endif
      ENDIF
c -------------------------------------------------------------------- c
c  supR --> up + gravitino

      IF(GRFLAG.EQ.1)THEN
        if (M32.le.asup2) then
          surUgra = asup2**5/(48d0*PI*MPL**2*M32**2)
        else
          surUgra = 0d0
        endif
      ENDIF
c -------------------------------------------------------------------- c
c                  SUM 2 BODY DECAYS
c -------------------------------------------------------------------- c
       supltot2lo = suplneutup(1)+suplneutup(2)+suplneutup(3)
     .           + suplneutup(4) + suplneutup(5)
     .           + suplchardow(1)+suplchardow(2)+suplglui+sulUgra

      suprtot2lo = suprneutup(1)+suprneutup(2)+suprneutup(3)
     .           + suprneutup(4)+suprneutup(5)
     .           + suprchardow(1)+suprchardow(2) +suprglui+surUgra

      supltot2nlo = supltot2lo + qcdsuplneutup(1)+qcdsuplneutup(2)
     .            + qcdsuplneutup(3)+qcdsuplneutup(4)+qcdsuplneutup(5)
     .            +qcdsuplchardow(1)+ qcdsuplchardow(2)+qcdsuplglui

      suprtot2nlo = suprtot2lo + qcdsuprneutup(1)+qcdsuprneutup(2)
     .            + qcdsuprneutup(3)+qcdsuprneutup(4)+qcdsuprneutup(5)
     .            + qcdsuprchardow(1)+ qcdsuprchardow(2)+qcdsuprglui

      if(flagqcd.eq.0d0) then
         supltot2 = supltot2lo
         suprtot2 = suprtot2lo
      elseif(flagqcd.eq.1d0) then
         supltot2 = supltot2nlo
         suprtot2 = suprtot2nlo
      endif

      if(flagqcd.eq.1d0) then
         do i=1,5
            suplneutup(i) = suplneutup(i)+qcdsuplneutup(i)
            suprneutup(i) = suprneutup(i)+qcdsuprneutup(i)
         enddo
         do i=1,2
            suplchardow(i) = suplchardow(i)+qcdsuplchardow(i)
            suprchardow(i) = suprchardow(i)+qcdsuprchardow(i)
         enddo
         suplglui = suplglui+qcdsuplglui
         suprglui = suprglui+qcdsuprglui
      endif

c---------------------------------------------------- c
c ----- sup_L/R branching ratios -------------------- c
c---------------------------------------------------- c

      if(supltot2.ne.0d0)then

       do i=1,5
         brsuplnup(i) = suplneutup(i)/supltot2
       enddo
       do i=1,2
         brsuplcdow(i) = suplchardow(i)/supltot2
       enddo
       brsuplglui = suplglui/supltot2
       brsulUgra = sulUgra/supltot2

      else

       do i=1,5
         brsuplnup(i) = 0d0
       enddo
       do i=1,2
         brsuplcdow(i) = 0d0
       enddo
       brsuplglui = 0d0
       brsulUgra = 0d0

      endif

      if(suprtot2.ne.0d0)then

       do i=1,5
         brsuprnup(i) = suprneutup(i)/suprtot2
       enddo
       do i=1,2
         brsuprcdow(i) = suprchardow(i)/suprtot2
       enddo
       brsuprglui = suprglui/suprtot2
       brsurUgra = surUgra/suprtot2

      else

       do i=1,5
         brsuprnup(i) = 0d0
       enddo
       do i=1,2
         brsuprcdow(i) = 0d0
       enddo
       brsuprglui = 0d0
       brsurUgra = 0d0

      endif

c -------------------------------------------------------------------- c
c  sdownl -> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + down

      do i=1,5
         if(amneut(i).le.asdown1) then
            sdowlneutdow(i)=g2s*((ado(1,i)**2+bdo(1,i)**2)*
     .           (asdown1**2-amneut(i)**2)
     .           )*NS_lamb(0d0,amneut(i)/asdown1)
     .           /(16*pi*asdown1)
         else
            sdowlneutdow(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,5
         if(amneut(i).le.asdown1) then
            qcdsdowlneutdow(i)=4d0/3d0*gs2/(4d0*pi)/pi*
     .         sdowlneutdow(i)*
     .         NS_ftotqcd(amneut(i)**2/asdown1**2,mgluino**2/asdown1**2)
         else
            qcdsdowlneutdow(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  sdownl -> chi-_1/chi-_2 up

      do i=1,2
         if(amchar(i).le.asdown1) then
            sdowlcharup(i)=g2s*aldo(1,i)**2*
     .           (asdown1**2-amchar(i)**2)*
     .           NS_lamb(0d0,amchar(i)/asdown1)/(16*pi*asdown1)
         else
            sdowlcharup(i)=0d0
         endif
      enddo
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,2
         if(amchar(i).le.asdown1) then
            qcdsdowlcharup(i)=4d0/3d0*gs2/(4d0*pi)/pi*
     .         sdowlcharup(i)*
     .         NS_ftotqcd(amchar(i)**2/asdown1**2,mgluino**2/asdown1**2)
         else
            qcdsdowlcharup(i)=0d0
         endif
      end do
      endif
c -------------------------------------------------------------------- c
c  sdownr -> chi0_1/chi0_2/chi0_3/chi0_4/chi0_5 + down

      do i=1,5
         if(amneut(i).le.asdown2) then
            sdowrneutdow(i)=g2s*((ado(2,i)**2+bdo(2,i)**2)*
     .           (asdown2**2-amneut(i)**2)
     .           )*NS_lamb(0d0,amneut(i)/asdown2)
     .           /(16*pi*asdown2)
         else
            sdowrneutdow(i)=0d0
         endif
      enddo

c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,5
         if(amneut(i).le.asdown2) then
            qcdsdowrneutdow(i)=4d0/3d0*gs2/(4d0*pi)/pi*
     .         sdowrneutdow(i)*
     .         NS_ftotqcd(amneut(i)**2/asdown2**2,mgluino**2/asdown2**2)
         else
            qcdsdowrneutdow(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  sdownr -> chi-_1/chi-_2 up

      do i=1,2
         if(amchar(i).le.asdown2) then
            sdowrcharup(i)=g2s*aldo(2,i)**2*
     .           (asdown2**2-amchar(i)**2)*
     .           NS_lamb(0d0,amchar(i)/asdown2)/(16*pi*asdown2)
         else
            sdowrcharup(i)=0d0
         endif
      enddo

c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      do i=1,2
         if(amchar(i).le.asdown2) then
            qcdsdowrcharup(i)=4d0/3d0*gs2/(4d0*pi)/pi*
     .         sdowrcharup(i)*
     .         NS_ftotqcd(amchar(i)**2/asdown2**2,mgluino**2/asdown2**2)
         else
            qcdsdowrcharup(i)=0d0
         endif
      enddo
      endif
c -------------------------------------------------------------------- c
c  sdownl --> gluino + down

      if(asdown1.gt.mgluino) then
         sdowlglui = 8d0*gs2*(asdown1**2-mgluino**2)*
     .            NS_lamb(0d0,mgluino/asdown1)/(16d0*pi*asdown1)/3d0
      else
         sdowlglui = 0d0
      endif
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      if(asdown1.gt.mgluino) then
         scala = amurefer
         alp   = gs2/(4d0*pi)
         ca    = 3d0
         cf    = 4d0/3d0
         amsq  = 2d0*(asup1+asup2+asdown1+asdown2)/8d0
         rval  = mgluino**2/amsq**2

         qcdsdowlglui = sdowlglui*alp/pi*( ca*NS_gama(rval) +
     .        cf*NS_gamfcap(rval) + 4d0*NS_gamf(rval) +
     .        2d0*pi**2/mgluino**2*
     .        NS_gamglui2(ast1,ast2,amt,thet,asb1,asb2,amb,theb,mgluino,
     .        1,amsq) +
     .        NS_gamrendec(amsq,ast1,ast2,amt,asb1,asb2,mgluino,scala) )
      else
         qcdsdowlglui = 0d0
      endif
      endif
c -------------------------------------------------------------------- c
c  sdownr --> gluino + down

      if(asdown2.gt.mgluino) then
         sdowrglui = 8d0*gs2*(asdown2**2-mgluino**2)*
     .             NS_lamb(0d0,mgluino/asdown2)/(16d0*pi*asdown2)/3d0
      else
         sdowrglui = 0d0
      endif
c --- QCD corrections ---
      if(flagqcd.eq.1d0) then
      if(asdown2.gt.mgluino) then
         scala = amurefer
         alp   = gs2/(4d0*pi)
         ca    = 3d0
         cf    = 4d0/3d0
         amsq  = 2d0*(asup1+asup2+asdown1+asdown2)/8d0
         rval  = mgluino**2/amsq**2

         qcdsdowrglui = sdowrglui*alp/pi*( ca*NS_gama(rval) +
     .        cf*NS_gamfcap(rval) + 4d0*NS_gamf(rval) +
     .        2d0*pi**2/mgluino**2*
     .        NS_gamglui2(ast2,ast1,amt,thet,asb2,asb1,amb,theb,mgluino,
     .        2,amsq) +
     .        NS_gamrendec(amsq,ast1,ast2,amt,asb1,asb2,mgluino,scala) )
      else
         qcdsdowrglui = 0d0
      endif
      endif
c -------------------------------------------------------------------- c
c  sdownL --> down + gravitino

      IF(GRFLAG.EQ.1)THEN
        if (M32.le.asdown1) then
          sdlDgra = asdown1**5/(48d0*PI*MPL**2*M32**2)
        else
          sdlDgra = 0d0
        endif
      ENDIF
c -------------------------------------------------------------------- c
c  sdownR --> down + gravitino

      IF(GRFLAG.EQ.1)THEN
        if (M32.le.asdown2) then
          sdrDgra = asdown2**5/(48d0*PI*MPL**2*M32**2)
        else
          sdrDgra = 0d0
        endif
      ENDIF
c -------------------------------------------------------------------- c
c                  SUM 2 BODY DECAYS
c -------------------------------------------------------------------- c
      sdowltot2lo = sdowlneutdow(1)+sdowlneutdow(2)+sdowlneutdow(3)+
     .              sdowlneutdow(4)+sdowlneutdow(5)+
     .              sdowlcharup(1)+sdowlcharup(2)+sdowlglui+sdlDgra

      sdowrtot2lo = sdowrneutdow(1)+sdowrneutdow(2)+sdowrneutdow(3)+
     .              sdowrneutdow(4)+sdowrneutdow(5)+
     .              sdowrcharup(1)+sdowrcharup(2)+sdowrglui+sdrDgra

      sdowltot2nlo = sdowltot2lo + qcdsdowlneutdow(1)+
     .               qcdsdowlneutdow(2)+qcdsdowlneutdow(3)+
     .               qcdsdowlneutdow(4)+qcdsdowlneutdow(5)+
     .               qcdsdowlcharup(1)+qcdsdowlcharup(2)+qcdsdowlglui

      sdowrtot2nlo = sdowrtot2lo + qcdsdowrneutdow(1)+
     .               qcdsdowrneutdow(2)+qcdsdowrneutdow(3)+
     .               qcdsdowrneutdow(4)+qcdsdowrneutdow(5)+
     .               qcdsdowrcharup(1)+qcdsdowrcharup(2)+qcdsdowrglui

      if(flagqcd.eq.0d0) then
         sdowltot2 = sdowltot2lo
         sdowrtot2 = sdowrtot2lo
      elseif(flagqcd.eq.1d0) then
         sdowltot2 = sdowltot2nlo
         sdowrtot2 = sdowrtot2nlo
      endif

      if(flagqcd.eq.1d0) then
c UE: Resum if qcdcorr < -tree:
         do i=1,5
!      If(qcdsdowlneutdow(i).lt.-sdowlneutdow(i))
!     .write(0,23)"Warning: large negative rad. corrs. to sdl->d+chi0_",i
       qcdsdowlneutdow(i)=resum(sdowlneutdow(i),qcdsdowlneutdow(i))
            sdowlneutdow(i) = sdowlneutdow(i)+qcdsdowlneutdow(i)
!      If(qcdsdowrneutdow(i).lt.-sdowrneutdow(i))
!     .write(0,23)"Warning: large negative rad. corrs. to sdr->d+chi0_",i
       qcdsdowrneutdow(i)=resum(sdowrneutdow(i),qcdsdowrneutdow(i))
            sdowrneutdow(i) = sdowrneutdow(i)+qcdsdowrneutdow(i)
         enddo
         do i=1,2
!      If(qcdsdowlcharup(i).lt.-sdowlcharup(i))
!     .write(0,23)"Warning: large negative rad. corrs. to sdl->u+char_",i
       qcdsdowlcharup(i)=resum(sdowlcharup(i),qcdsdowlcharup(i))
            sdowlcharup(i) = sdowlcharup(i)+qcdsdowlcharup(i)
!      If(qcdsdowrcharup(i).lt.-sdowrcharup(i))
!     .write(0,23)"Warning: large negative rad. corrs. to sdr->u+char_",i
       qcdsdowrcharup(i)=resum(sdowrcharup(i),qcdsdowrcharup(i))
            sdowrcharup(i) = sdowrcharup(i)+qcdsdowrcharup(i)
         enddo
!      If(qcdsdowlglui.lt.-sdowlglui)
!     .write(0,23)"Warning: large negative rad. corrs. to sdl->d+gluino"
       qcdsdowlglui=resum(sdowlglui,qcdsdowlglui)
         sdowlglui = sdowlglui+qcdsdowlglui
!      If(qcdsdowrglui.lt.-sdowrglui)
!     .write(0,23)"Warning: large negative rad. corrs. to sdr->d+gluino"
       qcdsdowrglui=resum(sdowrglui,qcdsdowrglui)
         sdowrglui = sdowrglui+qcdsdowrglui
      endif
23     format(A,I1)

c---------------------------------------------------- c
c ----- sdown_L/R branching ratios ------------------ c
c---------------------------------------------------- c

      if(sdowltot2.ne.0d0)then

       do i=1,5
         brsdowlndow(i) = sdowlneutdow(i)/sdowltot2
       enddo
       do i=1,2
         brsdowlchup(i) = sdowlcharup(i)/sdowltot2
       enddo
       brsdowlglui = sdowlglui/sdowltot2
       brsdlDgra = sdlDgra/sdowltot2

      else

       do i=1,5
         brsdowlndow(i) = 0d0
       enddo
       do i=1,2
         brsdowlchup(i) = 0d0
       enddo
       brsdowlglui = 0d0
       brsdlDgra = 0d0

      endif

      if(sdowrtot2.ne.0d0)then

       do i=1,5
         brsdowrndow(i) = sdowrneutdow(i)/sdowrtot2
       enddo
       do i=1,2
         brsdowrchup(i) = sdowrcharup(i)/sdowrtot2
       enddo
       brsdowrglui = sdowrglui/sdowrtot2
       brsdrDgra = sdrDgra/sdowrtot2

      else

       do i=1,5
         brsdowrndow(i) = 0d0
       enddo
       do i=1,2
         brsdowrchup(i) = 0d0
       enddo
       brsdowrglui = 0d0
       brsdrDgra = 0d0

      endif

      END

