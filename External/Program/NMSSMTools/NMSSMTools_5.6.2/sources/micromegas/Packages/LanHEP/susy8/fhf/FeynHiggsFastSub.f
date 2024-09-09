
      subroutine feynhiggssub(mh1,mh2,mh12,mh22)

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
      complex*16 P1se, P1se1, P2se, P2se1, P1P2se, P1P2se1
      complex*16 p1setl, p2setl, p1p2setl
      complex*16 LLExpansionP2MTrun, LLExpansionP2MTrun1 
      real*8 umix(1:2,1:2),vmix(1:2,1:2),nmix(1:4,1:4)
      real*8 mcha(1:2),mne(1:4)
      real*8 mma,ttb,mmt,mmm,mmue,mmsusy,au,ad,mh1,mh2,mh12,mh22
      real*8 l1, l2, l3, l4
c      double precision mlhlle1, mlhlle2, 
c     $       mlhllediag1, mlhllediag2, mhhllediag1, mhhllediag2
      integer ii,j,k,ic,pri,naeh,selec,selec2,selec4,selec5,selec6
      integer lleselec1,mpgut
      double precision m1zl, m1ol, m1tl, m2zl, m2ol, m2tl, 
     $                 mlhapp, mhhapp, mh2app, mh1app, 
     $                 azl, aol, atl,
     $                 mp12, mp22, mp1p22, mp12ol, mp22ol, mp1p22ol,
     $                 mp12tl, mp22tl, mp1p22tl,
     $                 mdiag1, mdiag2, mixang,
     $                 renlh, renhh, renxh, renlhsub, renhhsub, renxhsub
c      double precision allle1, allle2
      double precision alphasmz, alphas, mtrun, xttilde, vvv, ttt, 
     $                 delmlhsq, ms1, ms2, fac, alem
      double precision yepsilon,ymuee,ypi,ymz,ymw,mgf,yas,yalphasmz,
     $                 ycf,yeps,yeins,ymbb
      double precision lhseol, hhseol, xhseol, p1seol, p2seol, p1p2seol,
     $                 mlheff1sq, mhheff1sq, mlheff2sq, mhheff2sq,
     $                 mlheff1, mhheff1, aleff1, alefff1,
     $                 mlheff2, mhheff2, aleff2, alefff2
      complex*16 yi
      double precision mhptree, mhp1l

      common /smpara1/ ymw,ymz,mmt
      common /smpara2/ yepsilon,ymuee,ypi,ygf,ycf,yeps,yi,yeins,
     $                 yas,yalphasmz,ymbb
      common /susypara/ ttb,mma,mmm,mmue,au,ad
      common /gutrel/ mpgut
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
      common / err/ ic
      common /renpara/xo,zo,mgll
      common /print/pri,naeh,selec,selec2,selec4,selec5,selec6
      common /selftl/p1setl,p2setl,p1p2setl
c      common /lle/ mtlr, mblr, mlhlle1, mlhlle2,
c     $             mlhllediag1, mlhllediag2, mhhllediag1, mhhllediag2
c      common /lle2/ allle1, allle2
      common /lleselec/ lleselec1
      common /mhapp/ mlhapp, mhhapp, mh2app, mh1app
      common /zeromom/ mlheff1, mhheff1, alefff1, 
     $                 mlheff2, mhheff2, alefff2
      common /chargedhiggs/ mhptree, mhp1l
      double precision mtpole, mbpole
      common /polemasses/ mtpole, mbpole
      integer alphatsq
      common /alphat2/alphatsq
      double precision delp1, delp2, delp1p2


      integer delmbresum
      double precision dmb
      double precision msb1dmb, msb2dmb, stbdmb, tsbdmb
      common /deltambresum/dmb, msb1dmb, msb2dmb, stbdmb, tsbdmb, 
     $                     delmbresum

c passing over the SM parameters
      epsilon = yepsilon
      muee = ymuee
      pi = ypi
      mz = ymz
      mw = ymw
      gf = ygf
      alphasmz = yalphasmz
      cf = ycf
      eps = yeps
      i = yi
      eins = yeins
      mbb = ymbb

      ppi = pi

      mmz  = mz
      mz2 = mmz**2
      mmz2 = mmz**2
      mmw =  mw
      mw2 = mmw**2
      mmw2 = mmw**2
      cw2 = (mmw/mmz)**2
      ccw2 = (mmw/mmz)**2
      cw = dsqrt(ccw2)
      ccw = dsqrt(ccw2)
      sw2 = 1.d0 - ccw2
      ssw2 = 1.d0 - ccw2
      sw = dsqrt(ssw2)
      ssw = dsqrt(ssw2)
      el = dsqrt(8d0*gf*MMW**2*(1-MMW**2/MMZ**2)/dsqrt(2d0))
      elec = el
      elec2 = el**2
      alem = el**2/(4d0 * pi) 
      gmue = gf

      mtt=mmt
      mtt2 = mtt**2
      mt=mmt
      mb  = mbb
      mbb2 = mbb**2

      alphas = alphasmz/(1d0 + (11d0 - 10d0/3d0)/(4d0*ppi) * alphasmz
     $                       * dlog(mmt**2/mmz**2))
      as = alphas
      gs = dsqrt(4*ppi*as)
      call mtrunning(mmt, mtrun, alphas)

      maa=mma
      ma = mma

      mlheff1 = 0d0
      mlheff2 = 0d0


c softbreaking parameter

c-sh  Uebergabe
      mq = msusytl
      mq2 = msusytl**2
      mu2 = msusytr**2
      mb2 = msusybl**2  ! should be equal to mq2
c      mb2 = msusytl**2  ! should be equal to mq2
      md2 = msusybr**2
      mf2  = msusytl**2
      mfd2 = mf2
c      write(*,*) 'Msusy_top/bot_L:', real(msusytl), real(msusybl)

c --> passing over the parameters to Dabelstein's program
      mm=mmm
c --> mue convention:
c     Here the mue convention is fixed. Internally FeynHiggs works with
c     the MtLR = At + mue CTb. Due to the '-' sign below, externally
c     i.e. in the front-end one can work with the MtLR = At - mue CTb
c     convention.
      mu = -mmue
      mue = -mmue
      if (mpgut.eq.1) then
         mp = 5.d0/3.d0 * mm * (ssw/ccw)**2
c         write(*,*) 'M1:', real(mp)
      elseif (mpgut.ne.2) then
         write(*,*) 'M1 is not defined !!!'
         stop
      endif
      mgll = mgl

      mel = 0.51d-3
      mmu = 0.1057d0
      mta = 1.777d0
      mup = 0.0415d0
      mdn = 0.04151d0
      mch = 1.50d0
      mst = 0.150d0

c-sh  Uebergabe
      beta = datan(ttb)
      c2b = dcos(2d0*datan(ttb))
      cb = dcos(datan(ttb))
      sb = dsin(datan(ttb))

      mlh2 = 0.5d0*(maa**2+mmz**2 - dsqrt((maa**2+mmz**2)**2 - 4.d0*
     &              mmz**2*maa**2*dcos(2.d0*beta)**2))
      mlh = dsqrt(mlh2)
      Mhh2 = 0.5d0*(maa**2+mmz**2 + dsqrt((maa**2+mmz**2)**2 -
     &              4.d0*mmz**2*maa**2*dcos(2.d0*beta)**2))
      mhh = dsqrt(mhh2)
      mhp = dsqrt(maa**2 + mmw**2)
      mhptree = mhp

csh   new formula: eq (7.5) from S.H.'s PhD thesis
      alpha = datan((-(mma**2 + mmz**2) * sb * cb)/
     $              (mmz**2 * cb**2 + mma**2 * sb**2 - mlh**2))

c-sh  Uebergabe
      mssupq=au
c-sh  Uebergabe
      mssdnq=ad
      mssdnl =   mssdnq


      call mix
      call def2
      if (delmbresum.eq.2) then
         dmb = deltambnoew(alphas, mst1, mst2, msb1, msb2,
     $                         au, ad, mgl, mmt, mbb, mmue, ttb)
         write(*,*) 'evaluating Delta mb for Higgs masses:', real(dmb)
         call def4(msusybl, msusybr, ad - mmue * ttb, mbb/(1d0+dmb), 
     $             ttb, msb1dmb, msb2dmb, stbdmb, 2)
         tsbdmb = dasin(stbdmb)
         write(*,*) 'MSb1, MSb2, stb:'
         write(*,*) real(msb1dmb), real(msb2dmb), real(stbdmb)
         if ((msb1dmb.le.0d0).or.(msb2dmb.le.0d0)) then
            msb1dmb = msb1
            msb2dmb = msb2
            stbdmb = stb
            tsbdmb = dasin(stbdmb)
            write(*,*) 'WARNING: unphysical sbottom masses'
            write(*,*) ' --> using sbottom masses without Delta mb'//
     $                 'corrections'
            write(*,*) ' ==> results inconsistent!'
         endif
      endif

      if ((mst1.le.0d0).or.(mst2.le.0d0).or.
     $    (msb2.le.0d0).or.(msb2.le.0d0)) then
         mh1 = 119.9999d0
         mh2 = 119.9999d0
         mh12 = 119.9999d0
         mh22 = 119.9999d0
         goto 999
      endif


c --> find approximation for light Higgs mass
      p1se = p1se1()
      p2se = p2se1()
      p1p2se = p1p2se1()
c      write(*,*) 'self-energies at one loop:'
c      write(*,*) real(p1se), real(p2se), real(p1p2se)
      mlhapp = (.5d0 * (ma**2 + mz**2 - p2se - p1se) -
     $     (.25d0 * ((ma**2 + mz**2)**2 + (p2se - p1se)**2) -
     $      ma**2 * mz**2 * (dcos(2d0*beta))**2 +
     $      .5d0 * (p1se - p2se) * dcos(2d0*beta) * (ma**2 - mz**2) +
     $     p1p2se * dsin(2d0*beta) * (ma**2 + mz**2) +
     $     p1p2se**2 )**.5d0 )**.5d0
      mhhapp = (.5d0 * (ma**2 + mz**2 - p2se - p1se) +
     $     (.25d0 * ((ma**2 + mz**2)**2 + (p2se - p1se)**2) -
     $      ma**2 * mz**2 * (dcos(2d0*beta))**2 +
     $      .5d0 * (p1se - p2se) * dcos(2d0*beta) * (ma**2 - mz**2) +
     $     p1p2se * dsin(2d0*beta) * (ma**2 + mz**2) +
     $     p1p2se**2 )**.5d0 )**.5d0
c      write(*,*) "light higgs approximation: ", real(mlhapp)



c      write(*,*) 'vor mt-Aenderung'
      mtold = mt
c --> new top mass for two-loop contribution
      if (selec2.ge.2) then
c      write(*,*) 'mt = ', real(mt), real(mtt), real(mmt)
      write(*,*) "using running mt for two-loop contribution:", mtrun
      mt = mtrun
      if (selec4.eq.2) then
         call def2
         write(*,*) '... also for mt in Stop mass matrix'
      endif
      endif
      if ((mst1.le.0d0).or.(mst2.le.0d0).or.
     $    (msb2.le.0d0).or.(msb2.le.0d0)) then
c         write(*,*) 'Problem in Sfermion sector:'
c         write(*,*) 'MSt1, MSt2:', real(mst1), real(mst2)
c         write(*,*) 'MSb1, MSb2:', real(msb1), real(msb2)
         mh1 = 119.9999d0
         mh2 = 119.9999d0
         mh12 = 119.9999d0
         mh22 = 119.9999d0
         goto 999
      endif
c --> end of new parameter definition

c      write(*,*) 'vor 2-loop Berechnung'
c     calculation of the 2-loop contribution
      if (selec5.eq.1) then
         p1setl = 0d0
         p1p2setl = 0d0
      if (msusytl.eq.msusytr) then
         ms2 = dsqrt(msusytl**2 + mtrun**2) 
      else
         ms2 = msfkt(msusytl, msusytr, mtrun)
      endif                                             
      LLExpansionP2MTrun = LLExpansionP2MTrun1(mtrun,ms2)      
           p2setl = cf*el**2*as/
     $            (3d0 * 2d0**8 * sw**2 * pi**3 * sb**2)*
     $            mtrun**2/mw**2 * dreal(llexpansionp2mtrun)

c$$$         write(*,*) 'self-energies at two loop:'
c$$$         write(*,*) p1setl, p2setl, p1p2setl
c$$$         write(*,*) real(p1setl), real(p2setl), real(p1p2setl)
c$$$         write(*,*) real(mt), real(mst1), real(mst2), real(stt)
c$$$         write(*,*) real(cf), real(el), real(gs), real(sw), real(pi),
c$$$     $              real(sb), real(mw), real(mgl), real(eps)
c$$$         write(*,*) sb**2, (real(p2setl) * sb**2 *
c$$$     $   (1d0 + (4d0*MZ**2*(1d0 - 2d0*sb**2)*(1d0 - sb**2))/MA**2))
      else
         p1setl = 0d0
         p1p2setl = 0d0
         p2setl = 0d0
      endif


c----------------------------------------------------------------

      delmlhsq = 0d0
      if (selec2.ge.3) then
c     including the leading Yukawa term for the light Higgs mass
c     this term is taken from Carena, Espinoza, Quiros, Wagner
c     Nucl. Phys. B461 (1996) 407

      if (mh12.ne.119.9999d0) then

      write(*,*) 'including two-loop Yukawa term'

      if (alphatsq.eq.1) then
         write(*,*) "... using LL approximation (RGiEP)"
      if (selec4.eq.2) then
         write(*,*) "... also with running mt in stop mass matrix"
      endif

      xttilde = (((MSt2**2 - MSt1**2)/(4d0*mtrun**2) * 
     $            (2d0 * stt * dsqrt(1d0 - stt**2))**2 )**2 *
     $           (2d0 - (MSt2**2 + MSt1**2)/(MSt2**2 - MSt1**2) *
     $                  dlog(MSt2**2/MSt1**2)) +
     $           (MSt2**2 - MSt1**2)/(2d0 * mtrun**2) *
     $            (2d0 * stt * dsqrt(1d0 - stt**2))**2 *
     $                  dlog(MSt2**2/MSt1**2) )
      ttt = .5d0 * dlog(MSt2**2 * MSt1**2/mtrun**4)
      vvv = 174.1d0
      delmlhsq = ( 3d0/(4d0 * ppi**2) * mtrun**4/vvv**2 *
     $            ( 1d0/(16d0 * ppi**2) * 3d0/2d0 * mtrun**2/vvv**2 *
     $             ( xttilde * ttt + ttt**2)))

      p2setl = p2setl -1d0/sb**2 * delmlhsq

c --> new routine from hep-ph/0112177
      else
         write(*,*) "... using full calculation (EP 2-loop)"
         call BDSZHiggs(mtrun**2, mma**2, msusybl**2,
     $                  mst1**2, mst2**2, stt, dsqrt(1d0-stt**2),
     $                  mmt**2, mue, ttb, 246.218d0**2, 1, 
     $                  delp1, delp2, delp1p2)
         p1setl = p1setl - delp1
         p2setl = p2setl - delp2
         p1p2setl = p1p2setl - delp1p2
         delmlhsq = sb**2 * delp2
      endif
c --> end of new routine from hep-ph/0112177

      endif
      endif

c----------------------------------------------------------------
      

c --> two-loop approximations for the Higgs masses
c      write(*,*) 'pure top:', real(p1se), real(p2se), real(p1p2se)
      mh1app = (.5d0 * (ma**2 + mz**2 
     $     - (p2se + p2setl) - (p1se + p1setl)) -
     $     (.25d0 * ((ma**2 + mz**2)**2 
     $     + ((p2se + p2setl) - (p1se + p1setl))**2) -
     $      ma**2 * mz**2 * (dcos(2d0*beta))**2 +
     $      .5d0 * ((p1se + p1setl) 
     $     - (p2se + p2setl)) * dcos(2d0*beta) * (ma**2 - mz**2) +
     $     (p1p2se + p1p2setl) * dsin(2d0*beta) * (ma**2 + mz**2) +
     $     (p1p2se + p1p2setl)**2 )**.5d0 )**.5d0
      mh2app = (.5d0 * (ma**2 + mz**2 
     $     - (p2se + p2setl) - (p1se + p1setl)) +
     $     (.25d0 * ((ma**2 + mz**2)**2 
     $     + ((p2se + p2setl) - (p1se + p1setl))**2) -
     $      ma**2 * mz**2 * (dcos(2d0*beta))**2 +
     $      .5d0 * ((p1se + p1setl) 
     $     - (p2se + p2setl)) * dcos(2d0*beta) * (ma**2 - mz**2) +
     $     (p1p2se + p1p2setl) * dsin(2d0*beta) * (ma**2 + mz**2) +
     $     (p1p2se + p1p2setl)**2 )**.5d0 )**.5d0

  




c --> setting back the top mass to its original value for next calculation
      mt = mtold
c      write(*,*) 'mt = ', real(mt), real(mtt), real(mmt)
      if (selec4.eq.2) then
c --> setting back the sfermion masses
         call def2
      endif
c --> end of setting back the parameters

      

c --> charged Higgs at one loop
      mhp1l = dsqrt(ma**2 + mw**2 + 3d0/(4d0 * pi) * alem/(sw2 * mw2) *
     $            (  2d0 * mt**2 * mb**2/(sb**2 * cb**2) 
     $             - mw2 * (mt**2/sb**2 + mb**2/cb**2) 
     $            + 2d0/3d0 * mw2**2) * dlog(msusytl/mt)
     $            + mw2/(6d0 * pi) * 15d0 * alem/cw2 * dlog(msusytl/mw))


c --> neutral Higgs masses at two loop
      call mlhren(0.01d0, lhseol)
      lhseol = lhseol + mlh**2
      call mhhren(0.01d0, hhseol)
      hhseol = hhseol + mhh**2
      call mxhren(0.01d0, xhseol)
      xhseol = xhseol 
      p1seol =   dcos(alpha)**2 * hhseol 
     $         + dsin(alpha)**2 * lhseol
     $         - 2d0 * dsin(alpha) * dcos(alpha) * xhseol
      p2seol =   dsin(alpha)**2 * hhseol
     $         + dcos(alpha)**2 * lhseol
     $         + 2d0 * dsin(alpha) * dcos(alpha) * xhseol
      p1p2seol =   dsin(alpha) * dcos(alpha) * (hhseol - lhseol)
     $           + (dcos(alpha)**2 - dsin(alpha)**2) * xhseol

      mlheff1sq = (.5d0 * (ma**2 + mz**2 
     $                   - p2seol - p1seol) -
     $     (.25d0 * ((ma**2 + mz**2)**2 
     $               + (p2seol - p1seol)**2) -
     $      ma**2 * mz**2 * (dcos(2d0*beta))**2 +
     $      .5d0 * (p1seol - p2seol) * dcos(2d0*beta) 
     $                           * (ma**2 - mz**2) +
     $     p1p2seol * dsin(2d0*beta) * (ma**2 + mz**2) +
     $     p1p2seol**2 )**.5d0 )
      if (mlheff1sq.gt.0d0) then
         mlheff1 = (.5d0 * (ma**2 + mz**2 
     $                   - p2seol - p1seol) -
     $     (.25d0 * ((ma**2 + mz**2)**2 
     $               + (p2seol - p1seol)**2) -
     $      ma**2 * mz**2 * (dcos(2d0*beta))**2 +
     $      .5d0 * (p1seol - p2seol) * dcos(2d0*beta) 
     $                           * (ma**2 - mz**2) +
     $     p1p2seol * dsin(2d0*beta) * (ma**2 + mz**2) +
     $     p1p2seol**2 )**.5d0 )**.5d0
      else
         write(*,*) 'WARNING: error at effective light Higgs' //
     $              ' mass at 1 loop'
         mlheff1 = 0d0
      endif
      mhheff1sq = (.5d0 * (ma**2 + mz**2 
     $                   - p2seol - p1seol) +
     $     (.25d0 * ((ma**2 + mz**2)**2 
     $               + (p2seol - p1seol)**2) -
     $      ma**2 * mz**2 * (dcos(2d0*beta))**2 +
     $      .5d0 * (p1seol - p2seol) * dcos(2d0*beta) 
     $                           * (ma**2 - mz**2) +
     $     p1p2seol * dsin(2d0*beta) * (ma**2 + mz**2) +
     $     p1p2seol**2 )**.5d0 )
      if (mhheff1sq.gt.0d0) then
         mhheff1 = (.5d0 * (ma**2 + mz**2 
     $                   - p2seol - p1seol) +
     $     (.25d0 * ((ma**2 + mz**2)**2 
     $               + (p2seol - p1seol)**2) -
     $      ma**2 * mz**2 * (dcos(2d0*beta))**2 +
     $      .5d0 * (p1seol - p2seol) * dcos(2d0*beta) 
     $                           * (ma**2 - mz**2) +
     $     p1p2seol * dsin(2d0*beta) * (ma**2 + mz**2) +
     $     p1p2seol**2 )**.5d0 )**.5d0
      else
         write(*,*) 'WARNING: error at effective heavy Higgs' //
     $              ' mass at 1 loop'
         mhheff1 = 0d0
      endif
c      write(*,*) 'after mhheff1'
      if ((mlheff1.ne.0d0).and.(mhheff1.ne.0d0)) then
         aleff1 = datan((-(ma**2 + mz**2) * sb * cb - p1p2seol)/
     $            (mz**2 * cb**2 + ma**2 * sb**2 - p1seol - mlheff1**2))
      else
         aleff1 = 0d0
      endif
      alefff1 = aleff1


      mlheff2 = (.5d0 * (ma**2 + mz**2 
     $                   - (p2seol+p2setl) - (p1seol+p1setl)) -
     $     (.25d0 * ((ma**2 + mz**2)**2 
     $               + ((p2seol+p2setl) - (p1seol+p1setl))**2) -
     $      ma**2 * mz**2 * (dcos(2d0*beta))**2 +
     $      .5d0 * ((p1seol+p1setl) - (p2seol+p2setl)) * dcos(2d0*beta) 
     $                           * (ma**2 - mz**2) +
     $     (p1p2seol+p1p2setl) * dsin(2d0*beta) * (ma**2 + mz**2) +
     $     (p1p2seol+p1p2setl)**2 )**.5d0 )**.5d0
      mhheff2 = (.5d0 * (ma**2 + mz**2 
     $                   - (p2seol+p2setl) - (p1seol+p1setl)) +
     $     (.25d0 * ((ma**2 + mz**2)**2 
     $               + ((p2seol+p2setl) - (p1seol+p1setl))**2) -
     $      ma**2 * mz**2 * (dcos(2d0*beta))**2 +
     $      .5d0 * ((p1seol+p1setl) - (p2seol+p2setl)) * dcos(2d0*beta) 
     $                           * (ma**2 - mz**2) +
     $     (p1p2seol+p1p2setl) * dsin(2d0*beta) * (ma**2 + mz**2) +
     $     (p1p2seol+p1p2setl)**2 )**.5d0 )**.5d0
      aleff2 = datan((-(ma**2 + mz**2) * sb * cb 
     $                 - (p1p2seol+dreal(p1p2setl)))/
     $               (mz**2 * cb**2 + ma**2 * sb**2 
     $                 - (p1seol+dreal(p1setl)) - mlheff2**2))
      alefff2 = aleff2

      write(*,*) 'effective result:', real(ma), real(mz), real(beta)
      write(*,*) real(mlheff1), real(mhheff1), real(aleff1)
      write(*,*) real(mlheff2), real(mhheff2), real(aleff2)
      mh1 = mlheff1
      mh2 = mhheff1
      mh12 = mlheff2
      mh22 = mhheff2

c ----------------------------------------------------------------



 999  continue





      end

c ===========================================================

      subroutine mlhren (x, mout)
c
      implicit double precision (a-z)
c
      real*8 umix(1:2,1:2),vmix(1:2,1:2),nmix(1:4,1:4)
      real*8 mcha(1:2),mne(1:4)
      integer pri,naeh,selec,selec2,selec4,selec5,selec6
c
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
      common /renpara/xo,zo,mgll
      common /print/pri,naeh,selec,selec2,selec4,selec5,selec6
      common /sigmaasave/siasave, dsiasave
c
      q2 = x**2
c      write(*,*)'q2=',q2
c
c renorm.light scalar higgs mass
c
      call sigmaa  (maa**2 ,sab,sas,saf,sac,sat)
      if (selec.eq.4) then
         selfa = (sab + sas + saf      )
      elseif (selec.eq.3) then
         selfa = (sab + sas + saf + sac)
      elseif (selec.eq.2) then
         selfa = (      sas + saf      )
      elseif (selec.eq.1) then
         selfa = (      sas + saf      )
      else
         write(*,*) "Error in mlhren: selec out or range"
      endif
c      siasave = selfa
c MF 200801: selfa is now really the A-Selfenergy and not
c Dabelstein's combination of sigmaa and dsigmaa
c
c
      call tadlh (0.d0, sltb,slts,sltf,sltc,sltt)
c      write(*,*)'tadlh:',sltb,slts,sltf,sltc,sltt
      if (selec.eq.4) then
         tadlt = sltb + slts + sltf
      elseif (selec.eq.3) then
         tadlt = sltb + slts + sltf + sltc
      elseif (selec.eq.2) then
         tadlt =        slts + sltf         
      elseif (selec.eq.1) then
         tadlt =        slts + sltf         
c         write(*,*) 'tadlt:', slts, sltf
      else
         write(*,*) "Error in mlhren: selec out or range"
      endif
c     tadlt = slts + sltf
c     tadlt = sltt
c
      call tadhh (0.d0, shtb,shts,shtf,shtc,shtt)
c      write(*,*)'tadhh:',shtb,shts,shtf,shtc,shtt
      if (selec.eq.4) then
         tadht = shtb + shts + shtf
      elseif (selec.eq.3) then
         tadht = shtb + shts + shtf + shtc
      elseif (selec.eq.2) then
         tadht =        shts + shtf         
      elseif (selec.eq.1) then
         tadht =        shts + shtf         
      else
         write(*,*) "Error in mlhren: selec out or range"
      endif
c     tadht = shts + shtf
c     tadht = shtt
c
      call sigmaz (mmz**2 ,szb,szs,szf,szc,szt)
c      write(*,*)'sigmaz:',szb,szs,szf,szc,szt
      if (selec.eq.4) then
         selfz = szb + szs + szf
      elseif (selec.eq.3) then
         selfz = szb + szs + szf + szc
      elseif (selec.eq.2) then
         selfz =       szs + szf
      elseif (selec.eq.1) then
         selfz =       szs + szf
      else
         write(*,*) "Error in mlhren: selec out or range"
      endif


      dtb = (dzh0() - dzhh())/(2.0D0*dcos(2.0D0*alpha))
c MF 200801: tan(beta) Counterterm dtb contains the missing mudim-
c terms necessary for MSbar.

      deltamh = dcos(beta-alpha)**2 * selfa +
     &          elec/(2.d0*ssw*mmw) * dsin(beta-alpha)**2*dcos(beta-
     &          alpha) * tadht -
     &          elec/(2.d0*ssw*mmw) * dsin(beta-alpha) * (1.d0 +
     &          dcos(beta-alpha)**2) * tadlt +
     &          dsin(alpha+beta)**2 * selfz -
     &          dtb * dsin(2*beta) *
     &          (maa**2 * dsin(beta-alpha) * dcos(beta-alpha) -
     &           mmz**2 * dsin(alpha+beta) * dcos(alpha+beta))
c MF 200801: deltamh is now written in a general form including dtb
c      write(*,*) "nach deltamh",deltamh
c

      call sigmalh (q2, slhb,slhs,slhf,slhc,slht)
c      write(*,*)'sigmalh:',slhb,slhs,slhf,slhc,slht
      if (selec.eq.4) then
         selfh =  slhb + slhs + slhf
      elseif (selec.eq.3) then
         selfh =  slhb + slhs + slhf + slhc
      elseif (selec.eq.2) then
         selfh =         slhs + slhf 
      elseif (selec.eq.1) then
         selfh =         slhs + slhf 
      else
         write(*,*) "Error in mlhren: selec out or range"
      endif

      xxlh = selfh + dzh0()*(q2 - mlh**2) - deltamh
c MF 210801: Counterterm dZHH contains the missing mudim-terms
c necessary for MSbar.

c
      mout = xxlh + q2 - mlh**2
c
      
      
      return
      end
c
c ============================================================
c
      subroutine mhhren (x, mout)
c
      implicit double precision (a-z)
c
      real*8 umix(1:2,1:2),vmix(1:2,1:2),nmix(1:4,1:4)
      real*8 mcha(1:2),mne(1:4)
      integer pri,naeh,selec,selec2,selec4,selec5,selec6
c
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
      common /renpara/xo,zo,mgll
      common /print/pri,naeh,selec,selec2,selec4,selec5,selec6
c
      q2 = x**2
c
c renorm. heavy scalar higgs mass
c
      call sigmaa  (maa**2 ,sab,sas,saf,sac,sat)
      if (selec.eq.4) then
         selfa = (sab + sas + saf      )
      elseif (selec.eq.3) then
         selfa = (sab + sas + saf + sac)
      elseif (selec.eq.2) then
         selfa = (      sas + saf      )
      elseif (selec.eq.1) then
         selfa = (      sas + saf      )
      else
         write(*,*) "Error in mhhren: selec out or range"
      endif
c MF 200801: selfa is now really the A-Selfenergy and not
c Dabelstein's combination of sigmaa and dsigmaa
c
      call tadlh (0.d0, sltb,slts,sltf,sltc,sltt)
      if (selec.eq.4) then
         tadlt = sltb + slts + sltf
      elseif (selec.eq.3) then
         tadlt = sltb + slts + sltf + sltc
      elseif (selec.eq.2) then
         tadlt =        slts + sltf         
      elseif (selec.eq.1) then
         tadlt =        slts + sltf         
      else
         write(*,*) "Error in mhhren: selec out or range"
      endif

      call tadhh (0.d0, shtb,shts,shtf,shtc,shtt)
      if (selec.eq.4) then
         tadht = shtb + shts + shtf
      elseif (selec.eq.3) then
         tadht = shtb + shts + shtf + shtc
      elseif (selec.eq.2) then
         tadht =        shts + shtf         
      elseif (selec.eq.1) then
         tadht =        shts + shtf         
      else
         write(*,*) "Error in mhhren: selec out or range"
      endif

      call sigmaz (mmz**2 ,szb,szs,szf,szc,szt)
      if (selec.eq.4) then
         selfz = szb + szs + szf
      elseif (selec.eq.3) then
         selfz = szb + szs + szf + szc
      elseif (selec.eq.2) then
         selfz =       szs + szf
      elseif (selec.eq.1) then
         selfz =       szs + szf
      else
         write(*,*) "Error in mhhren: selec out or range"
      endif


      dtb = (dzh0() - dzhh())/(2.0D0*dcos(2.0D0*alpha))
c MF 200801: tan(beta) Counterterm dtb contains the missing mudim-
c terms necessary for MSbar.


      deltamh = dsin(beta-alpha)**2 * selfa -
     &          elec/(2.d0*ssw*mmw) * dcos(beta-alpha) * (1.d0 +
     &          dsin(beta-alpha)**2 ) * tadht +
     &          elec/(2.d0*ssw*mmw) * dcos(beta-alpha)**2 *
     &          dsin(beta-alpha) * tadlt +
     &          dcos(beta+alpha)**2 * selfz +
     &          dtb * dsin(2*beta) *
     &          (maa**2 * dsin(beta-alpha) * dcos(beta-alpha) -
     &           mmz**2 * dsin(alpha+beta) * dcos(alpha+beta))
c MF 200801: deltamh is now written in a general form including dtb
c
      call sigmahh (q2, slhb,slhs,slhf,slhc,slht)
      if (selec.eq.4) then
         selfH =  slhb + slhs + slhf
      elseif (selec.eq.3) then
         selfH =  slhb + slhs + slhf + slhc
      elseif (selec.eq.2) then
         selfH =         slhs + slhf 
      elseif (selec.eq.1) then
         selfH =         slhs + slhf 
c         write(*,*) 'subroutine mhhren:', real(slhs), real(slhf) 
      else
         write(*,*) "Error in mhhren: selec out or range"
      endif

      xxlh = selfH + dzhh()*(q2 - mhh**2) - deltamh
c MF 210801: Counterterm dZHH contains the missing mudim-terms
c necessary for MSbar.

c
      mout = xxlh + q2 - mhh**2
c      write(*,*) 'mhhren:', real(dsqrt(q2)), real(xxlh), real(mhh**2)

      return
      end
c
c ====================================================================
c
      subroutine mxhren (x, mout)
c
      implicit double precision (a-z)
c
      real*8 umix(1:2,1:2),vmix(1:2,1:2),nmix(1:4,1:4)
      real*8 mcha(1:2),mne(1:4)
      integer pri,naeh,selec,selec2,selec4,selec5,selec6
c
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
      common /renpara/xo,zo,mgll
      common /print/pri,naeh,selec,selec2,selec4,selec5,selec6

      q2 = x**2
c
c renorm. scalar higgs mixing
c
      call sigmaa  (maa**2 ,sab,sas,saf,sac,sat)
      if (selec.eq.4) then
         selfa = (sab + sas + saf      )
      elseif (selec.eq.3) then
         selfa = (sab + sas + saf + sac)
      elseif (selec.eq.2) then
         selfa = (      sas + saf      )
      elseif (selec.eq.1) then
         selfa = (      sas + saf      )
      else
         write(*,*) "Error in mxhren: selec out or range"
      endif
c MF 200801: selfa is now really the A0-Selfenergy and not
c Dabelstein's combination of sigmaa and dsigmaa
c
      call tadlh (0.d0, sltb,slts,sltf,sltc,sltt)
      if (selec.eq.4) then
         tadlt = sltb + slts + sltf
      elseif (selec.eq.3) then
         tadlt = sltb + slts + sltf + sltc
      elseif (selec.eq.2) then
         tadlt =        slts + sltf         
      elseif (selec.eq.1) then
         tadlt =        slts + sltf         
      else
         write(*,*) "Error in mlhren: selec out or range"
      endif

      call tadhh (0.d0, shtb,shts,shtf,shtc,shtt)
      if (selec.eq.4) then
         tadht = shtb + shts + shtf
      elseif (selec.eq.3) then
         tadht = shtb + shts + shtf + shtc
      elseif (selec.eq.2) then
         tadht =        shts + shtf         
      elseif (selec.eq.1) then
         tadht =        shts + shtf         
      else
         write(*,*) "Error in mlhren: selec out or range"
      endif

      call sigmaz (mmz**2 ,szb,szs,szf,szc,szt)
      if (selec.eq.4) then
         selfz = szb + szs + szf
      elseif (selec.eq.3) then
         selfz = szb + szs + szf + szc
      elseif (selec.eq.2) then
         selfz =       szs + szf
      elseif (selec.eq.1) then
         selfz =       szs + szf
      else
         write(*,*) "Error in mlhren: selec out or range"
      endif

      dtb = (dzh0() - dzhh())/(2.0D0*dcos(2.0D0*alpha))
c MF 200801: tan(beta) Counterterm dtb contains the missing mudim-
c terms necessary for MSbar.


      deltamh = dsin(alpha-beta)*dcos(alpha-beta) * selfa -
     &          elec/(2.d0*ssw*mmw) * dsin(beta-alpha)**3 * tadht -
     &          elec/(2.d0*ssw*mmw) * dcos(beta-alpha)**3 * tadlt -
     &          dsin(beta+alpha)*dcos(beta+alpha) * selfz -
     &          dtb * dsin(beta) * dcos(beta) *
     &          (maa**2 * (dcos(beta-alpha)**2-dsin(beta-alpha)**2) +
     &           mmz**2 * (dcos(alpha+beta)**2-dsin(alpha+beta)**2))
c MF 200801: deltamh is now written in a general form including dtb

      call sigmaxh (q2, slhb,slhs,slhf,slhc,slht)
      if (selec.eq.4) then
         selfx =  slhb + slhs + slhf
      elseif (selec.eq.3) then
         selfx =  slhb + slhs + slhf + slhc
      elseif (selec.eq.2) then
         selfx =         slhs + slhf 
      elseif (selec.eq.1) then
         selfx =         slhs + slhf 
      else
         write(*,*) "Error in mlhren: selec out or range"
      endif

      xxlh = selfx + dsin(2*alpha)*dtb*(q2 - (mhh**2 + mlh**2)/2.0D0) -
     &   deltamh
c MF 210801: the field renormalization (the term prop. to q2) contains
c the missing mudim-terms necessary for MSbar. Note: dtb is used for
c convenience only, it has simply the same analytic form as dZHh.

c
      mout = xxlh
c
      return
      end

c----------------------------------------------------------------
      double precision function msfkt(msusytl, msusytr, mtrun)
 
      double precision msusytl, msusytr, mtrun
 
c$$$       msfkt = dsqrt(dsqrt(msusytl**2 * msusytr**2
c$$$     $                    + mt**2 * (msusytl**2 + msusytr**2) + mt**4
c$$$     $                    + mt**2 * mtlr**2) - mt**2)
      msfkt = dsqrt( dsqrt(  msusytl**2 * msusytr**2
     $         + mtrun**2 * (msusytl**2 + msusytr**2) + mtrun**4))
 
      end
 
 
c----------------------------------------------------------------
c
c --> MSbartoOnShellPPtopbot.f
c
c----------------------------------------------------------------
      
      subroutine MSbartoOnShellPPtopbot(mst1msbar,mst2msbar,sttmsbar,
     $     qt1, qt2, qts,
     $     msb1msbar, msb2msbar, stbmsbar,
     $     qb1, qb2, qbs,
     $     mtos, mbos, xmgl, asmz, mz,
     $     mst1os, mst2os, sttos,
     $     msb1os, msb2os, stbos)

      complex*16 ShiftMSbarOSMSt1sq, ShiftMSbarOSMSt1sq1
      complex*16 ShiftMSbarOSMSt2sq, ShiftMSbarOSMSt2sq1
      complex*16 ShiftMSbarOSstt, ShiftMSbarOSstt1
      complex*16 ShiftMSbarOSMSb1sq, ShiftMSbarOSMSb1sq1
      complex*16 ShiftMSbarOSMSb2sq, ShiftMSbarOSMSb2sq1
      complex*16 ShiftMSbarOSstb, ShiftMSbarOSstb1

c -------------------------------------------------------------------
c varcom2.h
c
      double precision mst1, mst2, stt, mgl, mt, cf, gs, eps, pi, mue
      complex*16 i

      common /msbartoos/ mst1, mst2, stt, mgl, mt, cf, gs, eps, pi,
     $	                 mue, i
c -------------------------------------------------------------------
c -------------------------------------------------------------------
c varcom3.h
c
      double precision msb1, msb2, stb, mb

      common /msbartoos3/ msb1, msb2, stb, mb
c -------------------------------------------------------------------

      double precision mst1msbar, mst2msbar, sttmsbar, qt1, qt2, qts,
     $                 mst1os,    mst2os,    sttos
      double precision msb1msbar, msb2msbar, stbmsbar, qb1, qb2, qbs,
     $                 msb1os,    msb2os,    stbos
      double precision mst1org, mst2org, sttorg,
     $                 msb1org, msb2org, stborg
      double precision mtos, mbos, asmz, as, mz, xmgl
      double precision shiftt1, shiftt2, shifttx
      double precision shiftb1, shiftb2, shiftbx
      integer itert, iterb
      integer msbariter
      common /msbarselec2/ msbariter
 
c$$$      write(*,*) 'variables in MSbartoOnShellPP:'
c$$$      write(*,*) real(mst1msbar), real(mst2msbar), real(sttmsbar),
c$$$     $           real(qt1), real(qt2), real(qts), real(mtos), 
c$$$     $           real(msb1msbar), real(msb2msbar), real(stbmsbar),
c$$$     $           real(qb1), real(qb2), real(qbs), real(mbos), 
c$$$     $           real(xmgl), real(asmz), real(mz)
      eps = 1d-6
      cf = 4d0/3d0
      pi = 3.14159265358979d0
      i = (0d0, 1d0)

      mst1 = mst1msbar
      mst2 = mst2msbar
      stt = sttmsbar
      mt = mtos
      msb1 = msb1msbar
      msb2 = msb2msbar
      stb = stbmsbar
      mb = mbos
      call alphasmt(asmz,mt,mz,as)
      gs = dsqrt(4d0 * pi * as)
      mgl = xmgl
      itert = 0
      iterb = 0

      mst1org = mst1
      mst2org = mst2
      sttorg = stt
      msb1org = msb1
      msb2org = msb2
      stborg = stb

      write(*,*) 'MSbar to OS: Stop sector evaluation'
 80   continue

      mue = qt1
      ShiftMSbarOSMSt1sq = ShiftMSbarOSMSt1sq1()
      shiftt1 = dreal(ShiftMSbarOSMSt1sq)
c      write(*,*) 'shiftt1:', real(shiftt1)
      mue = qt2
      ShiftMSbarOSMSt2sq = ShiftMSbarOSMSt2sq1()
      shiftt2 = dreal(ShiftMSbarOSMSt2sq)
c      write(*,*) 'shiftt2:', real(shiftt2)
      mue = qts
      ShiftMSbarOSstt = ShiftMSbarOSstt1()
      shifttx = dreal(ShiftMSbarOSstt)
c      write(*,*) 'shifttx:', real(shifttx)

      if ((mst1org**2-shiftt1).ge.0d0) then
         mst1os = dsqrt(mst1org**2 - shiftt1)
      else
         mst1os = 0d0
         itert = 51
         write(*,*) 'no MSt1 OS'
      endif
      if ((mst2org**2-shiftt2).ge.0d0) then
         mst2os = dsqrt(mst2org**2 - shiftt2)
      else
         mst2os = 0d0
         itert = 51
         write(*,*) 'no MSt2 OS'
      endif
      sttos = sttorg - shifttx

      if (dabs(sttos).gt.1d0) then
      write(*,*) "MSbar to OnShell: |stt| > 1"
         if ((sttos-1d0).le.1d-2) then
            sttos = 1d0-1d-5
         elseif (dabs((-sttos-1d0)).le.1d-2) then
            sttos = -1d0+1d-5
         else
            sttos = 2d0
         endif
      endif

      if ((itert.le.50).and.(msbariter.eq.1)) then
c$$$         if ((dabs((mst1os - mst1)/mst1os).ge.1d-2).or.
c$$$     $       (dabs((mst2os - mst2)/mst2os).ge.1d-2).or.
c$$$     $       (dabs((sttos - stt)/sttos).ge.1d-2)) then
         if ((dabs((mst1os - mst1)).ge.1d0).or.
     $       (dabs((mst2os - mst2)).ge.1d0).or.
     $       (dabs((sttos - stt)).ge.1d-3)) then
            itert = itert + 1
            mst1 = mst1os
            mst2 = mst2os
            stt = sttos
            write(*,*) 
            write(*,*) "MSbar to OnShell: one more iteration:", itert
            write(*,*)
            write(*,*) "MSt1, MSt2, stt"
            write(*,*) real(mst1msbar), real(mst2msbar), real(sttmsbar)
            write(*,*) real(mst1os), real(mst2os), real(sttos)
            write(*,*) real(shiftt1), real(shiftt2), real(shifttx)
            goto 80
         endif
         elseif ((mst1os.ne.0d0).and.(mst2os.ne.0d0).and.
     $           (msbariter.eq.1)) then
         write(*,*) 'WARNING: More than 50 iterations in OS Stop calc.'
      endif


      write(*,*) 'MSbar to OS: Sbottom sector evaluation'
 85   continue

      mue = qb1
      ShiftMSbarOSMSb1sq = ShiftMSbarOSMSb1sq1()
      shiftb1 = dreal(ShiftMSbarOSMSb1sq)
      mue = qb2
      ShiftMSbarOSMSb2sq = ShiftMSbarOSMSb2sq1()
      shiftb2 = dreal(ShiftMSbarOSMSb2sq)
      mue = qbs
      ShiftMSbarOSstb = ShiftMSbarOSstb1()
      shiftbx = dreal(ShiftMSbarOSstb)

      msb1os = dsqrt(msb1org**2 - shiftb1)
      msb2os = dsqrt(msb2org**2 - shiftb2)
      stbos = stborg - shiftbx
c$$$      if (dabs(stbos).gt.1d0) then
c$$$         if ((stbos-1d0).le.1d-5) then
c$$$            stbos = 1d0-1d-5
c$$$         elseif (dabs((-stbos-1d0)).le.1d-5) then
c$$$            stbos = -1d0+1d-5
c$$$         else
c$$$            stbos = 2d0
c$$$         endif
c$$$      endif

      if ((iterb.le.50).and.(msbariter.eq.1)) then
         if ((dabs((msb1os - msb1)/msb1os).ge.1d-2).or.
     $       (dabs((msb2os - msb2)/msb2os).ge.1d-2).or.
     $       (dabs((stbos - stb)/stbos).ge.1d-2)) then
            iterb = iterb + 1
            msb1 = msb1os
            msb2 = msb2os
            stb = stbos
            write(*,*) 
            write(*,*) "MSbar to OnShell: one more iteration:", iterb
            write(*,*)
            write(*,*) "MSb1, MSb2, stb"
            write(*,*) real(msb1msbar), real(msb2msbar), real(stbmsbar)
            write(*,*) real(msb1os), real(msb2os), real(stbos)
            write(*,*) real(shiftb1), real(shiftb2), real(shiftbx)
            goto 85
         endif
         elseif ((msb1os.ne.0d0).and.(msb2os.ne.0d0).and.
     $           (msbariter.eq.1)) then
         write(*,*) 'WARNING: More than 20 iterations in OS Sbot calc.'
      endif



      write(*,*) 
      write(*,*) "MSbar to OnShell"
      write(*,*) "================"
      write(*,*)
      write(*,*) "MSt1, MSt2, stt"
      write(*,*) real(mst1msbar), real(mst2msbar), real(sttmsbar)
      write(*,*) real(mst1os), real(mst2os), real(sttos)
      write(*,*) real(shiftt1), real(shiftt2), real(shifttx)
      write(*,*) 
      write(*,*) "MSb1, MSb2, stb"
      write(*,*) real(msb1msbar), real(msb2msbar), real(stbmsbar)
      write(*,*) real(msb1os), real(msb2os), real(stbos)
      write(*,*) real(shiftb1), real(shiftb2), real(shiftbx)
      write(*,*) 

 100  continue

      end

c----------------------------------------------------------------
c
c --> MSbartoOnShellPP.f
c
c----------------------------------------------------------------
      
      subroutine MSbartoOnShellPP(mst1msbar,mst2msbar,sttmsbar,
     $     q1, q2, qs,
     $     mtos, xmgl, asmz, mz,
     $     mst1os, mst2os, sttos)

      complex*16 ShiftMSbarOSMSt1sq, ShiftMSbarOSMSt1sq1
      complex*16 ShiftMSbarOSMSt2sq, ShiftMSbarOSMSt2sq1
      complex*16 ShiftMSbarOSstt, ShiftMSbarOSstt1

c -------------------------------------------------------------------
c varcom2.h
c
      double precision mst1, mst2, stt, mgl, mt, cf, gs, eps, pi, mue
      complex*16 i

      common /msbartoos/ mst1, mst2, stt, mgl, mt, cf, gs, eps, pi,
     $	                 mue, i
c -------------------------------------------------------------------

      double precision mst1msbar, mst2msbar, sttmsbar, q1, q2, qs,
     $                 mst1os,    mst2os,    sttos
      double precision mtos, asmz, as, mz, xmgl
      double precision shift1, shift2, shiftx
 
c      write(*,*) 'variables in MSbartoOnShellPP:'
c      write(*,*) real(mst1msbar), real(mst2msbar), real(sttmsbar),
c     $           real(q1), real(q2), real(qs), real(mtos), real(xmgl),
c     $           real(asmz), real(mz)
      eps = 1d-6
      cf = 4d0/3d0
      pi = 3.14159265358979d0
      i = (0d0, 1d0)

      mst1 = mst1msbar
      mst2 = mst2msbar
      stt = sttmsbar
      mt = mtos
      call alphasmt(asmz,mt,mz,as)
      gs = dsqrt(4d0 * pi * as)
      mgl = xmgl

      mue = q1
      ShiftMSbarOSMSt1sq = ShiftMSbarOSMSt1sq1()
      shift1 = dreal(ShiftMSbarOSMSt1sq)
      mue = q2
      ShiftMSbarOSMSt2sq = ShiftMSbarOSMSt2sq1()
      shift2 = dreal(ShiftMSbarOSMSt2sq)
      mue = qs
      ShiftMSbarOSstt = ShiftMSbarOSstt1()
      shiftx = dreal(ShiftMSbarOSstt)

      mst1os = dsqrt(mst1**2 - shift1)
      mst2os = dsqrt(mst2**2 - shift2)
      sttos = stt - shiftx



      write(*,*) 
      write(*,*) "MSbar to OnShell"
      write(*,*) "================"
      write(*,*)
      write(*,*) "MSt1, MSt2, stt"
      write(*,*) real(mst1msbar), real(mst2msbar), real(sttmsbar)
      write(*,*) real(mst1os), real(mst2os), real(sttos)
      write(*,*) real(shift1), real(shift2), real(shiftx)
      write(*,*) 

 100  continue

      end



