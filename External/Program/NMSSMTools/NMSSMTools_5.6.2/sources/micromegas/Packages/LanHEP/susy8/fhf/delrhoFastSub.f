      subroutine delrho(mst1, mst2, stt, msb1, msb2, stb, mgl, 
     $           mt, mb, xgf, as, cf, gs, el, mz, mw, 
     $           mh12eff, mh22eff, ma, aleff, beta,
     $           res1loop, res2loopgluon, res2loopgluino,
     $           res2loopmt2yuk, res2loopmt2yuksm)

      double precision mst1, mst2, stt, ctt, mt, 
     $                 msb1, msb2, stb, ctb, mb, mgl, 
     $                 gf, as, pi, res, xgf, cf, gs, el, mw, mz, pref
      double precision oneloop, twoloop, xyf0, xyf1,
     $                 res1loop, res2loopgluon, res2loopgluino,
     $                 res2loopmt2yuk, res2loopmt2yuksm
      complex*16 delrhogluino

      double precision xmh12, xmh22, xma, xsa, xca,
     $                 mh12eff, mh22eff, ma, aleff, beta
      common/xhiggs/ xmh12, xmh22, xma, xsa, xca

      integer delmbresum
      double precision dmb, mbbdmb
      double precision msb1dmb, msb2dmb, stbdmb, tsbdmb
      double precision msb1save, msb2save, stbsave
      common /deltambresum/dmb, msb1dmb, msb2dmb, stbdmb, tsbdmb, 
     $                     delmbresum
      mbbdmb = mbb/(1d0 + dmb)

      if (delmbresum.eq.2) then
         msb1save = msb1
         msb2save = msb2
         stbsave = stb
         msb1 = msb1dmb
         msb2 = msb2dmb
         stb = stbdmb
      endif

c$$$      write(*,*) 'paramters in delrho:'
c$$$      write(*,*) 'mst1, mst2, stt, msb1, msb2, stb, xgf, as, mt, mb'
c$$$      write(*,*) 'mgl, cf, gs, el, mz, mw, mh, mH, MA, al, be'
c$$$      write(*,*) real(mst1), real(mst2), real(stt), real(msb1), 
c$$$     $           real(msb2), real(stb), real(xgf), real(as),
c$$$     $           real(mt), real(mb), real(mgl), real(cf), real(gs),
c$$$     $           real(el), real(mz), real(mw), real(mh12eff),
c$$$     $           real(mh22eff), real(ma), real(aleff), real(beta)

      pi = 3.14159265358979d0
      gf = xgf
      ctt = dsqrt(1d0 - stt**2)
      ctb = dsqrt(1d0 - stb**2)
      pref = (3d0 * cf * gs**2 * el**2 * mz**2)/
     $       (2**6 * mw**2 * (mw**2 - mz**2) * pi**4)
      prefmt2yuk = (3d0 * el**4 * mt**4 * mz**4)/
     $             (2**12 * mw**4 * (mw**2 - mz**2)**2 * pi**4)
      xmh12 = mh12eff
      xmh22 = mh22eff
      xsa = dsin(aleff)
      xca = dcos(aleff)
      xma = ma
 
      oneloop = 3d0 * gf/(8d0 * dsqrt(2d0) * pi**2) *
     $          ( - stt**2 * ctt**2 * xyf0(mst1**2, mst2**2)
     $            - stb**2 * ctb**2 * xyf0(msb1**2, msb2**2)
     $            + ctt**2 * ctb**2 * xyf0(mst1**2, msb1**2)
     $            + ctt**2 * stb**2 * xyf0(mst1**2, msb2**2)
     $            + stt**2 * ctb**2 * xyf0(mst2**2, msb1**2)
     $            + stt**2 * stb**2 * xyf0(mst2**2, msb2**2) )

      twoloop = gf * as/(4d0 * dsqrt(2d0) * pi**3) *
     $          ( - stt**2 * ctt**2 * xyf1(mst1**2, mst2**2)
     $            - stb**2 * ctb**2 * xyf1(msb1**2, msb2**2)
     $            + ctt**2 * ctb**2 * xyf1(mst1**2, msb1**2)
     $            + ctt**2 * stb**2 * xyf1(mst1**2, msb2**2)
     $            + stt**2 * ctb**2 * xyf1(mst2**2, msb1**2)
     $            + stt**2 * stb**2 * xyf1(mst2**2, msb2**2) )

      delrhogluino = 0d0

c      write(*,*) 'DelrhoSub:'
      res1loop = oneloop
      res2loopgluon = twoloop
      res2loopgluino = dreal(delrhoGluino) * pref


      res2loopmt2yuk = 0d0
      res2loopmt2yuksm = 0d0
c      write(*,*) 'MT2Yuk:', res2loopmt2yuk, res2loopmt2yuksm, prefmt2yuk




      if (delmbresum.eq.2) then
         msb1 = msb1save
         msb2 = msb2save
         stb = stbsave
         ctb = dsqrt(1d0 - stb**2)
      endif

      end


c-------------------------------------------------------------------

      double precision function xyf0(x,y)

      double precision x, y

      if (x.ne.y) then
c         write(*,*) 'xyf0:', x, y
      xyf0 = x + y - (2d0 * x * y)/(x - y) * dlog(x/y)
      else
      xyf0 = 0d0
      endif

      end


c-------------------------------------------------------------------

      double precision function xyf1(x2,y2)

      double precision x2, y2
      complex*16 ff, x, y, cspen
      
      x = (dsqrt(x2) - (0d0,1d0) * 10d-10)**2
      y = (dsqrt(y2) - (0d0,1d0) * 10d-10)**2

      if (x.ne.y) then
c         write(*,*) 'xyf1:', x, y
      ff = x + y 
     $   - (2d0 * x * y)/(x - y) * cdlog(x/y) * (2d0 + x/y * cdlog(x/y))
     $   + (x + y)*x**2/(x - y)**2 * (cdlog(x/y))**2 
     $   - 2d0 * (x - y) * cspen(1d0 - x/y)

      xyf1 = dreal(ff)
      else
      xyf1 = 0d0
      endif

      end


c-------------------------------------------------------------------

      double precision function switchoff(m)

      double precision m,so

      so = (m - 750d0)/125d0
      if (so.le.0d0) so = 0d0
      
      switchoff = so**2

      end

c-------------------------------------------------------------------




      
