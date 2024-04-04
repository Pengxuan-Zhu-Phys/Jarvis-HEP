      double precision function fhf2f(itb, imu, iml, imr, imbr,
     $        iat, iab, img2, img3, ima,
     $   imt, imb, itp)
     
      implicit real*8(a-z)

      double precision itb,imu,iml,imr,imbr,
     $                    iat,iab,img2,img3,ima,imt,imb
      integer itp


      double precision MSt1, MSt2, Mgl, MT, MB, MW, MZ, MA
     $               , stt, ctt, stb, ctb  
     $               , MSb1, MSb2, Mue, PI, sw2, sw, cw
     $               , cf, el, gs, a, as, gf
     $               , tb, b, c2b, sb, cb, pref, eps, eins
     $               , msusytl, msusytr, msusybl, msusybr, mlrt, mlrb
     $               , x2, delmst, msusytaul, msusytaur
      complex*16 cspen, i, res, res1, res2, res3, res4, res5, res6
      integer r, s, t, dr1l
      double precision MSmuLtot, MSmuRtot, MSmuneut

      common/masses/MSt1, MSt2, MSb1, MSb2, Mgl, Mue, delmst
      common/input/msusytl, msusytr, msusybl, msusybr, mlrt, mlrb,
     $             msusytaul, msusytaur
      common/prec/tb, b, c2b, sb, cb, MZ, MW, MA, sw2, sw, cw, MT, MB, 
     $             gf, as, el, a, gs, stb, cf, stt, eps, i, eins, pi
      common /Sbottomshift/ dr1l
      common /SmuonSector/ MSmuLtot, MSmuRtot, MSmuneut

      double precision xmh12, xmh22, xma, xsa, xca
      common/xhiggs/ xmh12, xmh22, xma, xsa, xca
c -------------------------------------------------------------------

      real*8 umix(1:2,1:2),vmix(1:2,1:2),nmix(1:4,1:4)
      real*8 mcha(1:2),mne(1:4)
      real*8 mma,ttb,mmt,mmm,mmue,mmsusy,au,ad,mh1,mh2,mh12,mh22,l1,l2
      real*8 mtlr,mblr
      double precision mlhlle1, mlhlle2, 
     $       mlhllediag1, mlhllediag2, mhhllediag1, mhhllediag2
      double precision allle1, allle2
      integer ii,j,k,ic,pri,naeh,selec,selec2,selec3,selec4,
     $        selec5,selec6, sps, lhbms
      integer lleselec1,msbarselec,mtmsbarselec, msbarselecsave, mpgut
      double precision delrholimit, delrholimitnew, delrhores,
     $                 delrho1loop, delrho2loopgluon, delrho2loopgluino,
     $                 prefdelrhogluino,
     $                 delrho2loopmt2yuk, delrho2loopmt2yuksm
      double precision msusytrnew
      double precision yepsilon,ymuee,ypi,ymz,ymw,mgf,yas,yalphasmz,
     $                 ycf,yeps,yeins,ymbb
      complex*16 yi
      double precision mst1msbar, mst2msbar, sttmsbar, qt1, qt2, qts,
     $                 mst1os, mst2os, sttos,
     $                 msb1msbar, msb2msbar, stbmsbar, qb1, qb2, qbs,
     $                 mbs1os, msb2os, stbos
      double precision xmsusytl, xmsusytr, xmsusybl, xmsusybr, 
     $                 xmtlr, xmblr
      double precision lhseol, hhseol, xhseol, p1seol, p2seol, p1p2seol,
     $                 mlheff1, mhheff1, aleff1, alefff1,
     $                 mlheff2, mhheff2, aleff2, alefff2
      double precision mhptree, mhp1l
      double precision msusytlmsbar, msusytrmsbar, msusybrmsbar, 
     $                 atmsbar, abmsbar, xtmsbar, xbmsbar, qmsbar

c      external fapr,gundu,gundutl,ftest

      common /smpara1/ ymw,ymz,mmt
      common /smpara2/ yepsilon,ymuee,ypi,ygf,ycf,yeps,yi,yeins,
     $                 yas,yalphasmz,ymbb
      common /susypara/ ttb,mma,mmm,mmue,au,ad
      common /param/ ssw2,ssw,ccw2,ccw,ppi,elec2,elec,mmz,mmz2,mmw2,mmw,
     &               beta,alpha
      common /singl/ epsilon,muee,lambda
      common /susyset/ mu,mm,mp
      common /mass/ mel,mmu,mta,mup,mdn,mch,mst,mbb,mbb2,mtt,mtt2,
     &              melsl,mmusl,mtasl,mupsl,mvesl,mvmsl,mvtsl,
     &              mdnsl,mstsl,mchsl,mtsl,mbsl,mhh,mlh,maa,mhp,
     &              melsr,mmusr,mtasr,mupsr,mvesr,mvmsr,mvtsr,
     &              mdnsr,mstsr,mchsr,mtsr,mbsr, mcha,mne
      common /mixing/ umix,vmix,nmix
      common /fangle/ ang1,ang2,ang3,ang4,ang5,ang6,ang7,ang8,ang9,
     &                ang10,ang11,ang12
      common /abreak/mssupq,mssdnq,mssdnl
      common /break/ mq2,mu2,mb2,md2,mf2,mfd2
      common /chargedhiggs/ mhptree, mhp1l
      common / err/ ic
      common /renpara/xo,zo,mgll
      common /print/pri,naeh,selec,selec2,selec4,selec5,selec6
      common /lle/ mtlr, mblr, mlhlle1, mlhlle2,
     $             mlhllediag1, mlhllediag2, mhhllediag1, mhhllediag2
      common /lle2/ allle1, allle2
      common /mhapp/ mlhapp, mhhapp, mh2app, mh1app
      common /zeromom/ mlheff1, mhheff1, alefff1, 
     $                 mlheff2, mhheff2, alefff2
      common /lleselec/ lleselec1
      common /msbarselec/ msbarselec, mtmsbarselec
      common /gutrel/ mpgut
      integer msbariter
      common /msbarselec2/ msbariter
      double precision mtpole, mbpole
      common /polemasses/ mtpole, mbpole
      integer alphatsq
      common /alphat2/alphatsq

      double precision mudim
      common /msbar/ mudim


      integer delmbresum
      double precision dmb
      double precision msb1dmb, msb2dmb, stbdmb, tsbdmb
      common /deltambresum/dmb, msb1dmb, msb2dmb, stbdmb, tsbdmb, 
     $                     delmbresum
      integer error


      double precision jtb,jmu,jml,jmr,jat,jma,jmt,jmbr,
     $            jab,jmg3,jmg2,jmb,jmhh,jmhl
      common/dsav/ jtb,jmu,jml,jmr,jat,jma,jmt,jmbr,
     $  jab,jmg3,jmg2,jmb,jmhh,jmhl
      save
      
      if(itb.ne.jtb .or. imu.ne.jmu .or. iml.ne.jml
     $ .or. imr.ne.jmr .or. iat.ne.jat .or. ima.ne.jma
     $ .or. imt.ne.jmt .or. imb.ne.jmb .or. imbr.ne.jmbr
     $ .or. jab.ne.iab .or. img2.ne.jmg2 .or. img3.ne.jmg3) then
            
      ttb=itb
      mmue=imu
      msusytl=iml
      msusytr=imr
      msusybr=imbr
      mtlr=iat
      mblr=iab
      mbpole=imb
      mmt=imt
      mgl=img3
      mma=ima
      mmm=img2
      
      jtb=itb
      jmu=imu
      jml=iml
      jmr=imr
      jmbr=imbr
      jat=iat
      jab=iab
      jma=ima
      jmg2=img2
      jmg3=img3
      jmt=imt
      jmb=imb
      
      call FeynHiggsFast()
      endif
            
      if(itp.eq.1) fhf2f=jmhl
      if(itp.eq.2) fhf2f=jmhh
      if(itp.eq.3) fhf2f=sin(alefff2)
      if(itp.eq.4) fhf2f=cos(alefff2)
      if(itp.eq.5) fhf2f=real(mhp1l)
      return
      end

      subroutine FeynHiggsFast()

c --------------------------------------------------------------
c
c     FeynHiggsFast
c     =============
c      
c       Calculation of the masses of the neutral CP-even
c       Higgs bosons in the MSSM
c       
c       Authors: Sven Heinemeyer (one-, two-loop part, new renormalization)
c                Andreas Dabelstein (one-loop part)
c                Markus Frank (new renormalization)
c       
c       Based on hep-ph/9803277, hep-ph/9807423, hep-ph/9812472,
c                hep-ph/9903404, hep-ph/9910283
c       by S. Heinemeyer, W. Hollik, G. Weiglein
c       and on hep-ph/0001002
c       by M. Carena, H. Haber, S. Heinemeyer, W. Hollik,
c          C. Wagner and G. Weiglein
c       new non-log O(alpha_t^2) corrections taken from hep-ph/0112177
c       by A. Brignole, G. Degrassi, P. Slavich and F. Zwirner
c
c       new renormalization implemented based on hep-ph/0202166
c       by M. Frank, S. Heinemeyer, W. Hollik and G. Weiglein
c      
c       In case of problems or questions,
c       contact Sven Heinemeyer
c       email: Sven.Heinemeyer@physik.uni-muenchen.de
c       
c       FeynHiggs homepage:
c       http://www.feynhiggs.de
c
c --------------------------------------------------------------


      implicit real*8(a-z)
c -------------------------------------------------------------------
c varcom.h
c
      double precision MSt1, MSt2, Mgl, MT, MB, MW, MZ, MA
     $               , stt, ctt, stb, ctb  
     $               , MSb1, MSb2, Mue, PI, sw2, sw, cw
     $               , cf, el, gs, a, as, gf
     $               , tb, b, c2b, sb, cb, pref, eps, eins
     $               , msusytl, msusytr, msusybl, msusybr, mlrt, mlrb
     $               , x2, delmst, msusytaul, msusytaur
      complex*16 cspen, i, res, res1, res2, res3, res4, res5, res6
      integer r, s, t, dr1l
      double precision MSmuLtot, MSmuRtot, MSmuneut

      common/masses/MSt1, MSt2, MSb1, MSb2, Mgl, Mue, delmst
      common/input/msusytl, msusytr, msusybl, msusybr, mlrt, mlrb,
     $             msusytaul, msusytaur
      common/prec/tb, b, c2b, sb, cb, MZ, MW, MA, sw2, sw, cw, MT, MB, 
     $             gf, as, el, a, gs, stb, cf, stt, eps, i, eins, pi
      common /Sbottomshift/ dr1l
      common /SmuonSector/ MSmuLtot, MSmuRtot, MSmuneut

      double precision xmh12, xmh22, xma, xsa, xca
      common/xhiggs/ xmh12, xmh22, xma, xsa, xca
c -------------------------------------------------------------------

      real*8 umix(1:2,1:2),vmix(1:2,1:2),nmix(1:4,1:4)
      real*8 mcha(1:2),mne(1:4)
      real*8 mma,ttb,mmt,mmm,mmue,mmsusy,au,ad,mh1,mh2,mh12,mh22,l1,l2
      real*8 mtlr,mblr
      double precision mlhlle1, mlhlle2, 
     $       mlhllediag1, mlhllediag2, mhhllediag1, mhhllediag2
      double precision allle1, allle2
      integer ii,j,k,ic,pri,naeh,selec,selec2,selec3,selec4,
     $        selec5,selec6, sps, lhbms
      integer lleselec1,msbarselec,mtmsbarselec, msbarselecsave, mpgut
      double precision delrholimit, delrholimitnew, delrhores,
     $                 delrho1loop, delrho2loopgluon, delrho2loopgluino,
     $                 prefdelrhogluino,
     $                 delrho2loopmt2yuk, delrho2loopmt2yuksm
      double precision msusytrnew
      double precision yepsilon,ymuee,ypi,ymz,ymw,mgf,yas,yalphasmz,
     $                 ycf,yeps,yeins,ymbb
      complex*16 yi
      double precision mst1msbar, mst2msbar, sttmsbar, qt1, qt2, qts,
     $                 mst1os, mst2os, sttos,
     $                 msb1msbar, msb2msbar, stbmsbar, qb1, qb2, qbs,
     $                 mbs1os, msb2os, stbos
      double precision xmsusytl, xmsusytr, xmsusybl, xmsusybr, 
     $                 xmtlr, xmblr
      double precision lhseol, hhseol, xhseol, p1seol, p2seol, p1p2seol,
     $                 mlheff1, mhheff1, aleff1, alefff1,
     $                 mlheff2, mhheff2, aleff2, alefff2
      double precision mhptree, mhp1l
      double precision msusytlmsbar, msusytrmsbar, msusybrmsbar, 
     $                 atmsbar, abmsbar, xtmsbar, xbmsbar, qmsbar

c      external fapr,gundu,gundutl,ftest

      common /smpara1/ ymw,ymz,mmt
      common /smpara2/ yepsilon,ymuee,ypi,ygf,ycf,yeps,yi,yeins,
     $                 yas,yalphasmz,ymbb
      common /susypara/ ttb,mma,mmm,mmue,au,ad
      common /param/ ssw2,ssw,ccw2,ccw,ppi,elec2,elec,mmz,mmz2,mmw2,mmw,
     &               beta,alpha
      common /singl/ epsilon,muee,lambda
      common /susyset/ mu,mm,mp
      common /mass/ mel,mmu,mta,mup,mdn,mch,mst,mbb,mbb2,mtt,mtt2,
     &              melsl,mmusl,mtasl,mupsl,mvesl,mvmsl,mvtsl,
     &              mdnsl,mstsl,mchsl,mtsl,mbsl,mhh,mlh,maa,mhp,
     &              melsr,mmusr,mtasr,mupsr,mvesr,mvmsr,mvtsr,
     &              mdnsr,mstsr,mchsr,mtsr,mbsr, mcha,mne
      common /mixing/ umix,vmix,nmix
      common /fangle/ ang1,ang2,ang3,ang4,ang5,ang6,ang7,ang8,ang9,
     &                ang10,ang11,ang12
      common /abreak/mssupq,mssdnq,mssdnl
      common /break/ mq2,mu2,mb2,md2,mf2,mfd2
      common /chargedhiggs/ mhptree, mhp1l
      common / err/ ic
      common /renpara/xo,zo,mgll
      common /print/pri,naeh,selec,selec2,selec4,selec5,selec6
      common /lle/ mtlr, mblr, mlhlle1, mlhlle2,
     $             mlhllediag1, mlhllediag2, mhhllediag1, mhhllediag2
      common /lle2/ allle1, allle2
      common /mhapp/ mlhapp, mhhapp, mh2app, mh1app
      common /zeromom/ mlheff1, mhheff1, alefff1, 
     $                 mlheff2, mhheff2, alefff2
      common /lleselec/ lleselec1
      common /msbarselec/ msbarselec, mtmsbarselec
      common /gutrel/ mpgut
      integer msbariter
      common /msbarselec2/ msbariter
      double precision mtpole, mbpole
      common /polemasses/ mtpole, mbpole
      integer alphatsq
      common /alphat2/alphatsq

      double precision mudim
      common /msbar/ mudim


      integer delmbresum
      double precision dmb
      double precision msb1dmb, msb2dmb, stbdmb, tsbdmb
      common /deltambresum/dmb, msb1dmb, msb2dmb, stbdmb, tsbdmb, 
     $                     delmbresum
      integer error
      double precision jtb,jmu,jml,jmr,jat,jma,jmt,jmbr,
     $            jab,jmg3,jmg2,jmb,jmhh,jmhl
      common/dsav/ jtb,jmu,jml,jmr,jat,jma,jmt,jmbr,
     $  jab,jmg3,jmg2,jmb,jmhh,jmhl

      save

c ----------------------------------------------------------------
c setting the SM parameters 
c ----------------------------------------------------------------
      epsilon = 1.01d0
      muee     = 1.d0      

      pi = 3.14159265897d0

      mz  = 91.187d0
      mw =  80.451d0
      gf = 1.16639d-5
      alphasmz = 0.118d0
      cf = 4d0/3d0
      eps = 1d-10
      i = (0d0, 1d0)
      eins = 1d0

      delrholimit = 2d-3


c passing over the SM parameters
      yepsilon = epsilon
      ymuee = muee
      ypi = pi
      ymz = mz
      ymw = mw
      ygf = gf
      yas = as
      yalphasmz = alphasmz
      ycf = cf
      yeps = eps
      yi = i
      yeins = eins
      ymbb = mbb

c -----------------------------------------------------------------
c
c  switches
c  ========
c
c  selec2/ii = 1: full 1-loop + 2-loop QCD
c              2: same as 1, but in addition with
c                           mt = 166.5 at 2-loop
c              3: same as 2, but in addition with
c                           Yukawa term added for light Higgs
c
c  selec4 = 1: mt = pole mass at 2-loop in the stop mass matrix
c           2: mt = running mass
c           this this an option only if selec2/ii > 1,
c           otherwise selec4 = 1 is automatically chosen.
c
c  selec = 1: top/stop sector only in 1-loop calculation
c          2: top/stop + bottom/sbottom sector only
c          3: full 1-loop calculation
c
c  selec5 = 1: full result for the two-loop SEs is used
c           2: LL-Expansion for the two-loop SEs is used
c
c  selec6 = 2: the short formula from the LLE is also used in order to
c              to derive the Higgs masses 
c              You get two values for the Higgs masses then: 
c              the 'complete' one and the 'short' one
c         = 1: the short formula is not used, you get only one result
c
c  lleselec1 = 1: the LLE contribution is evaluated with the short
c                 (and published) formula
c            = 2: the LLE contribution is evaluated with the long
c                 (and unpublished) formula -> LLExpansionP2.f
c
c  msbarselec = 1: input parameters are on-shell parameter
c             = 2: input parameters are MSbar parameters
c                  (they are then transformed into the on-shell scheme
c                   in order to use them for FeynHiggs)
c
c  alphatsq = 1: alpha_t^2 corrections in LL approximation (RGiEP)
c           = 2:           full calculation (EP 2-loop)
c
c  delmbresum = 1: Delta mb corrections are not included in mh calculation
c             = 2: Delta mb corrections are included   
c
c -----------------------------------------------------------------

c     single point calculation with online parameter input


      

 90   continue

      sps = 0
      lhbms = 0
      selec5 = 1
      selec6 = 1
      lleselec1 = 1
      mpgut = 1
      mtmsbarselec = 1
      msbarselec = 2
      msbarselecsave = msbarselec
      ii=4
      selec3=1
      selec4=1
      delmbresum=1
      print1=1
      print2=1
      print3=1
      delrholimit=0
      
      msbariter = 0
      if (selec3.ne.1) selec4 = 1
      if (delrholimit.eq.0d0) delrholimit = 2d-3

      if (ii.le.3) alphatsq = 1
      if (ii.eq.4) alphatsq = 2
      mtselec = 0


c --> Input of parameters

      if (selec3.eq.1) then
c --> Input: OS, Msusy, Xt, Xb

c      read(*,*) ttb, msusytl, msusytr, msusybr, mtlr, mblr, 
c     $          mmt, mbpole, mgl, mmue, mmm, mma, mudim, selec
      if (ttb.eq.0d0) goto 90
      mudim=1
      selec=3
      mmsusy = msusytl
      msusybl = msusytl
      if (msusytr.eq.0d0) msusytr = mmsusy
      if (msusybr.eq.0d0) msusybr = mmsusy
      if (mmt.eq.0d0) mmt = 175d0
      if (mmt.eq.1d0) mmt = 174.3d0
      if (mbpole.eq.0d0) mbpole = 4.5d0
      if (mbpole.eq.1d0) mbpole = 4.25d0
      mbb = 2.97d0
      mtpole = mmt
      if (mgl.eq.0d0) mgl = 500d0
      if (mgl.eq.1d0) mgl = mmsusy
      if (mmue.eq.1d0) mmue = 200d0
      if (mmue.eq.2d0) mmue = mmsusy
      if (mmm.eq.1d0) mmm = 400d0
      if (mmm.eq.2d0) mmm = mmsusy
      if (mudim.eq.0d0) mudim = mmt/2d0
      if (mudim.eq.1d0) mudim = mmt
      if (mudim.eq.2d0) mudim = mmt*2d0
      

      else
         write(*,*) 'internal error for selec3... stopping'
         stop
      endif




c --> Evaluating the OS parametes (if necessary...)

      if (selec3.eq.1) then
c --> Input: OS, Msusy, Xt, Xb

      au = mtlr + mmue/ttb
      if (mblr.eq.1d0) then
         ad = au
         mblr = ad - mmue * ttb
      endif
      ad = mblr + mmue * ttb
      mlrt = mtlr
      mlrb = mblr

      elseif (selec3.eq.2) then
c --> Input: OS, MSt1, MSt2, stt, MSb1, MSb2, stb

      call def4b(mst1,mst2,stt,mmt,2d0/3d0,
     $           ttb,msusytl,msusytr,mtlr,1)
      call def4b(msb1,msb2,stb,mbb,-1d0/3d0,
     $           ttb,msusybl,msusybr,mblr,2)
      if ((msusytl.eq.0d0).or.(msusybl.eq.0d0)) then
         write(*,*) 'inconsistency in the sfermion sector'
         write(*,*) '... skipping calculation'
         mh12 = 119.9999d0
         goto 100
      endif
      au = mtlr + mmue/ttb
      ad = mblr + mmue*ttb
      mlrt = mtlr
      mlrb = mblr
      write(*,*)
      write(*,*) 'Transition to unphysical parameters completed'

      elseif (selec3.eq.3) then
c --> Input: MSbar, MSt1, MSt2, stt, MSb1, MSb2, stb

      mst1msbar = mst1
      mst2msbar = mst2
      sttmsbar = stt
      mtos = mmt
      msb1msbar = msb1
      msb2msbar = msb2
      stbmsbar = stb
      mbos = mbb
c      write(*,*) 'vor MSbartoOnShellPPtopbot'
      call MSbartoOnShellPPtopbot(mst1msbar,mst2msbar,sttmsbar,
     $     qt1, qt2, qts,
     $     msb1msbar, msb2msbar, stbmsbar,
     $     qb1, qb2, qbs,
     $     mtos, mbos, mgl, yalphasmz, ymz,
     $     mst1os, mst2os, sttos,
     $     msb1os, msb2os, stbos)
      mst1 = mst1os
      mst2 = mst2os
      stt = sttos
      msb1 = msb1os
      msb2 = msb2os
      stb = stbos

      call def4b(mst1,mst2,stt,mmt,2d0/3d0,
     $           ttb,msusytl,msusytr,mtlr,1)
      call def4b(msb1,msb2,stb,mbb,-1d0/3d0,
     $           ttb,msusybl,msusybr,mblr,2)
      if ((msusytl.eq.0d0).or.(msusybl.eq.0d0)) then
         write(*,*) 'inconsistency in the sfermion sector'
         write(*,*) '... skipping calculation'
         goto 100
      endif

      au = mtlr + mmue/ttb
      ad = mblr + mmue*ttb
c      al = mllr + mmue*ttb
      mlrt = mtlr
      mlrb = mblr
c      write(*,*) 'At, Ab:', real(au), real(ad)
      write(*,*) 'Transition to unphysical OS parameters completed'


      elseif (selec3.eq.4) then
c --> Input: MSbar, MSt1, MSt2, stt, Msusy_bot_R, Xb

      mst1msbar = mst1
      mst2msbar = mst2
      sttmsbar = stt
      mtos = mmt
c      write(*,*) 'vor MSbartoOnShellPP'
      call MSbartoOnShellPP(mst1msbar,mst2msbar,sttmsbar,
     $     qt1, qt2, qts,
     $     mtos, mgl, alphasmz, ymz,
     $     mst1os, mst2os, sttos)
      mst1 = mst1os
      mst2 = mst2os
      stt = sttos

      call def5(ttb,mmt,ymw,ymz,mst1,mst2,stt, msusybr, mblr,
     $          xmsusytl,xmsusytr,xmsusybl,xmsusybr,xmtlr,xmblr)
c      write(*,*) 'Results of def5:'
c      write(*,*) real(xmsusytl), real(xmsusytr), real(xmsusybl),
c     $           real(xmsusybr), real(xmtlr), real(xmblr)
      msusytl = xmsusytl
      msusytr = xmsusytr
      msusybl = xmsusybl
      msusybr = xmsusybr
      if (msusybr.eq.1d0) msusybr = msusybl
      mtlr = xmtlr
      mblr = xmblr
      if (msusytl.eq.0d0) then 
         write(*,*) 'Stop masses and the mixing angle do not fit'
     $              // ' together'
         goto 100
      endif
      au = mtlr + mmue/ttb
      ad = mblr + mmue*ttb
      if (mblr.eq.1d0) ad = au
      mblr = ad - mmue*ttb
      mlrt = mtlr
      mlrb = mblr
      write(*,*) 'Transition to unphysical OS parameters completed'


      else
         write(*,*) 'internal error for selec3... stopping'
         stop
      endif
            



      ymbb = mbb

      if (print1.eq.0) then
      write(*,*)
      write(*,*) "Your OS parameters:"
      write(*,*) "MT, Msusy(top-left), Msusy(top-right), Xt, At"
      write(*,*) "MB, Msusy(bot-left), Msusy(bot-rigth), Xb, Ab"
      write(*,*) "tb, Mgl, Mu, M2, MA"
      write(*,*) real(mmt), real(msusytl), real(msusytr), real(mtlr),
     $                                                    real(au)
      write(*,*) real(mbb), real(msusybl), real(msusybr), real(mblr),
     $                                                    real(ad)
      write(*,*) real(ttb), real(mgl), real(mmue), real(mmm), real(mma)
      write(*,*)
      endif


      selec2 = ii
      dr1l = 0
c      write(*,*) "-------------------------------------------------"
c      write(*,*) "-------------------------------------------------"
c      write(*,*) "Performing the calculation..."
       write(*,*) 
     $ "FeynHiggsFast 1.2.2 (Aug 09 2002) (implemented Apr 29 2003)"  
      call feynhiggssub(mh1,mh2,mh12,mh22)
c      write(*,*) "-------------------------------------------------"
c      write(*,*) "-------------------------------------------------"
c      write(*,*) mst1, mst2
      if (   (real(mst1).le.0d0).or.(real(mst2).le.0d0).
     $    or.(real(msb1).le.0d0).or.(real(msb2).le.0d0)) then
         write(*,*) "negative entry in sfermion mass matrix"
         write(*,*) "MSt1, MSt2: ", real(mst1), real(mst2)
         write(*,*) "MSb1, MSb2: ", real(msb1), real(msb2)
         goto 100
      endif
      if (   (real(mst1).le.65d0).or.(real(mst2).le.65d0).
     $    or.(real(msb1).le.65d0).or.(real(msb2).le.65d0)) then
         write(*,*) "WARNING - WARNING - WARNING - WARNING"
         write(*,*) "experimental excluded sfermion masses:"
         write(*,*) "MSt1, MSt2: ", real(mst1), real(mst2)
         write(*,*) "MSb1, MSb2: ", real(msb1), real(msb2)
         write(*,*) "WARNING - WARNING - WARNING - WARNING"
      else
         if (print2.eq.0) then
         write(*,*) "The Sfermion masses:"
         write(*,*) "MSt1, MSt2: ", real(mst1), real(mst2)
         write(*,*) "MSb1, MSb2: ", real(msb1), real(msb2)
         ctt = dsqrt(1d0 - stt**2)
         ctb = dsqrt(1d0 - stb**2)
         write(*,*) "cos(mixing angles): ", real(ctt), real(ctb)
         write(*,*) "sin(mixing angles): ", real(stt), real(stb)
         endif
      endif
      if (print3.eq.0) then
         write(*,*) "Chargino Masses:", real(mcha(1)), real(mcha(2))
         write(*,*) "The mixing matrices:"
         write(*,*) "             U                             V"
         write(*,*) real(umix(1,1)), real(umix(1,2)),
     $              real(vmix(1,1)), real(vmix(1,2))
         write(*,*) real(umix(2,1)), real(umix(2,2)),
     $              real(vmix(2,1)), real(vmix(2,2))
         write(*,*) "Neutralino Masses:"
         write(*,*) real(mne(1)), real(mne(2)), 
     $              real(mne(3)), real(mne(4))
         write(*,*) "The mixing matrix:"
         write(*,*) real(nmix(1,1)), real(nmix(1,2)),
     $              real(nmix(1,3)), real(nmix(1,4))
         write(*,*) real(nmix(2,1)), real(nmix(2,2)),
     $              real(nmix(2,3)), real(nmix(2,4))
         write(*,*) real(nmix(3,1)), real(nmix(3,2)),
     $              real(nmix(3,3)), real(nmix(3,4))
         write(*,*) real(nmix(4,1)), real(nmix(4,2)),
     $              real(nmix(4,3)), real(nmix(4,4))
      endif

c      write(*,*) "----------------------------------------------------"
      write(*,*) "----------------------------------------------------"
      write(*,*) "The results:  light Higgs     heavy Higgs     alpha"
      write(*,*) "----------------------------------------------------"
      write(*,*) "mh-tree :   ", real(mlh), real(mhh), real(alpha)
      write(*,*) "----------------------------------------------------"
      write(*,*) "mh-1loop:"
      if ((mh1.ne.119.9999d0).and.(mh2.ne.119.9999d0)) then
         write(*,*) " --> BEST  :", real(mh1), real(mh2), real(alefff1)
      else 
         write(*,*) " --> BEST  : Higgs sector not ok at 1-loop"
      endif
      write(*,*) "Yuk-approx :", real(mlhapp), real(mhhapp)
      if (((dabs(mlhapp - mh1).ge.10d0).or.
     $     (dabs(mhhapp - mh2).ge.10d0)).and.
     $    (mh1.ne.119.9999d0).and.(mh2.ne.119.9999d0)) then
         write(*,*) "WARNING - WARNING - WARNING - WARNING"
         write(*,*) "possible numerical instability detected"
         write(*,*) "WARNING - WARNING - WARNING - WARNING"
      endif
      write(*,*) "----------------------------------------------------"
      write(*,*) "mh-2loop:"
      if ((mh12.ne.119.9999d0).and.(mh22.ne.119.9999d0)) then
         write(*,*) " --> BEST  :", real(mh12),real(mh22), real(alefff2)
      else 
         write(*,*) " --> BEST  : Higgs sector not ok at 2-loop"
      endif
      jmhl=mh12
      jmhh=mh22
      write(*,*) "Yuk-approx :", real(mh1app), real(mh2app)
      if (((dabs(mh1app - mh12).ge.10d0).or.
     $     (dabs(mh2app - mh22).ge.10d0)).and.
     $    (mh1.ne.119.9999d0).and.(mh2.ne.119.9999d0)) then
         write(*,*) "WARNING - WARNING - WARNING - WARNING"
         write(*,*) "possible numerical instability detected"
         write(*,*) "WARNING - WARNING - WARNING - WARNING"
      endif
      write(*,*) "----------------------------------------------------"
c      write(*,*) "----------------------------------------------------"
      write(*,*) "charged Higgs, tree:", real(mhptree)
      write(*,*) "             1-loop:", real(mhp1l)


 99   continue

      write(*,*) "----------------------------------------------------"
c      write(*,*) "----------------------------------------------------"
      call delrho(mst1, mst2, stt, msb1, msb2, stb, mgl, mmt, mbb, 
     $            gf, as, cf, gs, el, mz, mw,
     $            mlheff2, mhheff2, mma, alefff2, datan(ttb), 
     $            delrho1loop, delrho2loopgluon, delrho2loopgluino,
     $            delrho2loopmt2yuk, delrho2loopmt2yuksm)
      delrhores = delrho1loop + delrho2loopgluon + delrho2loopgluino +
     $            (delrho2loopmt2yuk - delrho2loopmt2yuksm)
      if (dabs(delrhores).gt.delrholimit) then
         write(*,*) 'WARNING: Delta rho > experimental limit'
      endif
c      write(*,*) 'Delta rho 1-loop         : ', delrho1loop
c      write(*,*) 'Delta rho 2-loop (gluon) : ', delrho2loopgluon
c      write(*,*) 'Delta rho 2-loop (gluino): ', delrho2loopgluino
c      write(*,*) 'Delta rho 2-loop (MT4eff): ', 
c     $               (delrho2loopmt2yuk - delrho2loopmt2yuksm)
      write(*,*) 'Delta rho total          : ', delrhores
      write(*,*) "----------------------------------------------------"


c      write(*,*)
c      write(*,*) "collected WARNINGS:"
      if (   (real(mst1).le.65d0).or.(real(mst2).le.65d0).
     $    or.(real(msb1).le.65d0).or.(real(msb2).le.65d0)) then
         write(*,*) "experimental excluded sfermion masses"
      endif
      if ((mh1.eq.119.9999d0).or.(mh2.eq.119.9999d0)) then
         write(*,*) " --> BEST  : Higgs sector not ok at 1-loop"
      endif
c      if (mlheff1.eq.0d0) then
c         write(*,*) "zero mom.  : Higgs sector not ok at 1-loop"
c      endif
      if (((dabs(mlhapp - mh1).ge.10d0).or.
     $     (dabs(mhhapp - mh2).ge.10d0)).and.
     $    (mh1.ne.119.9999d0).and.(mh2.ne.119.9999d0)) then
         write(*,*) "possible numerical instability detected"
      endif
      if ((mh12.eq.119.9999d0).or.(mh22.eq.119.9999d0)) then
         write(*,*) " --> BEST  : Higgs sector not ok at 2-loop"
      endif
c      if (mlheff2.eq.0d0) then
c         write(*,*) "zero mom.  : Higgs sector not ok at 2-loop"
c      endif
      if (((dabs(mh1app - mh12).ge.10d0).or.
     $     (dabs(mh2app - mh22).ge.10d0)).and.
     $    (mh1.ne.119.9999d0).and.(mh2.ne.119.9999d0)) then
         write(*,*) "possible numerical instability detected"
      endif
      if (dabs(delrhores).gt.delrholimit) then
         write(*,*) 'Delta rho > experimental limit'
      endif
c      write(*,*)
c      write(*,*) "----------------------------------------------------"
c      write(*,*) "----------------------------------------------------"
      error = 0
c      write(*,*) "General WARNINGS:"
c      write(*,*)
      if (dmb.gt.0.8d0) then
         write(*,*) "Delta mb > 0.8"
         error = 1
      endif
      if (dmb.lt.-0.8d0) then
         write(*,*) "Delta mb < -0.8"
         error = 1
      endif
      if (dabs(mtlr)/dsqrt(msusytl*msusytr).ge.2.3d0) then
         write(*,*) "|Xt|/Msusy > 2.3"
         error = 1
      endif
      if (error.eq.1) then 
         write(*,*) "possible implication: numerical instability"
      endif
      if (sps.eq.11) then
         write(*,*) "WARNING: SPS1b implemented only as approximation"
      endif
c      write(*,*) "----------------------------------------------------"
c      write(*,*) "----------------------------------------------------"
      write(*,*)



 100  continue
      
      end
