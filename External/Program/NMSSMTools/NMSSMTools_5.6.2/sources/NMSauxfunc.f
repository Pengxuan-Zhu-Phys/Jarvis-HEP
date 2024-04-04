************************************************************************
*
*     These subroutines are generalisations of corresponding routines
*     from SDECAY: A Fortran code for the decays of the supersymmetric
*                  particles in the MSSM
*     by M. Muhlleitner (Karlsruhe, Inst. Technol.),
*	 A. Djouadi (Orsay, LPT & CERN, Theory Division),
*	 Y. Mambrini (Orsay, LPT),
*     Comput.Phys.Commun.168:46-70 (2005), hep-ph/0311167.
*     SDECAY should be cited whenever NMSDECAY is used.
*
************************************************************************

c ------------------------------------------------------------------c
c ------------------ FUNCTION resum ------------------------------- c
c ------------------------------------------------------------------c
c --- x: tree level, y: rad. corr., returns y or y_modif. --------- c

      DOUBLE PRECISION FUNCTION resum(x,y)
      IMPLICIT NONE
      DOUBLE PRECISION x,y

      if(x.ne.0d0) then
        if(y.le.-x) then
          resum=x*y/(x-y)
        else
          resum=y
        endif
      else
        resum=0d0
      endif
      return
      end

c -------------------------------------------------------------------- c
c ------------------ the FUNCTION lambda ----------------------------- c
c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_lamb(x,y)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)

      NS_lamb=dsqrt((1d0-x**2-y**2)**2-4d0*x**2*y**2)

      return
      end

c -------------------------------------------------------------------- c
c ---------------------- the integration limits ---------------------- c
c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_ax(xmu1)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      NS_ax=2d0*dsqrt(xmu1)
      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_BX(xmu1,xmu2,xmu3)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      NS_BX=1d0+xmu1-(dsqrt(xmu2)+dsqrt(xmu3))**2
      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_ay(xmu1,xmu2,xmu3,x1)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)

      a = 1d0-x1+xmu1
      b = (x1-2d0)*(x1-1d0-xmu2+xmu3-xmu1)

      delta = (4d0*xmu1-x1**2)*
     .        (4d0*xmu2*xmu3-(x1-1d0+xmu2+xmu3-xmu1)**2)

      if (delta.lt.0d0) then
         NS_ay=0d0
      else
         r1=(b-dsqrt(delta))/(2d0*a)
         r2=(b+dsqrt(delta))/(2d0*a)
         NS_ay=r1
      endif

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_BY(xmu1,xmu2,xmu3,x1)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)

      a = 1d0-x1+xmu1
      b = (x1-2d0)*(x1-1d0-xmu2+xmu3-xmu1)

      delta = (4d0*xmu1-x1**2)*
     .        (4d0*xmu2*xmu3-(x1-1d0+xmu2+xmu3-xmu1)**2)

      if (delta.lt.0d0) then
         NS_BY=0d0
      else
         r1=(b-dsqrt(delta))/(2d0*a)
         r2=(b+dsqrt(delta))/(2d0*a)
         NS_BY=r2
      endif

      end

c -------------------------------------------------------------------- c
c ----------------------- the integration routine -------------------- c
c -------------------------------------------------------------------- c

      SUBROUTINE NS_INTEG2(F,NS_AX,NS_BX,NS_AY,NS_BY,xmu1,xmu2,xmu3,
     .                     NX,NY,SUM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DOUBLE PRECISION NS_AY,NS_BY,NS_AX,NS_BX
      DIMENSION RX(1000),WX(1000),RY(1000),WY(1000)
      EXTERNAL F,NS_AY,NS_BY,NS_AX,NS_BX

      AXX=NS_AX(xmu1)
      BXX=NS_BX(xmu1,xmu2,xmu3)

      CALL NS_GAUS(NX,AXX,BXX,RX,WX)

      SX=0d0
      do  1 K=1,NX
         X=RX(K)
         AYX=NS_AY(xmu1,xmu2,xmu3,X)
         BYX=NS_BY(xmu1,xmu2,xmu3,X)
         CALL NS_GAUS(NY,AYX,BYX,RY,WY)
         SY=0d0
         do 2 J=1,NY
            SY=SY+WY(J)*F(X,RY(J))
 2       CONTINUE
         SX=SX+WX(K)*SY
 1    CONTINUE

      SUM=SX
      end

c -------------------------------------------------------------------- c

      SUBROUTINE NS_GAUS(N,A,B,X,W)
C     GAUSS-LEGENDRE INTEGRATION FROM A TO B (WEIGHT = 1.)
C     CALCULATES GAUSSIAN POINTS X AND WEIGHTS W
C                      FOR N=4,6,8,12,16,24,32 ;
C     if N IS DIFFERENT FROM THESE VALUES THE PROGRAM DECOMPOSES
C     THE INTERVAL [A,B] AND N INTO CORRESPONDING PEACES
C     N MUST BE EVEN AND >= 4,if IT IS NOT,IT IS CHANGED !
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION XG(16,7),DG(16,7),YG(16),EG(16),X(1000),W(1000),NI(7),
     .          NG(7)
      DATA NBEG,NA,NI,NG /9*0,4,6,8,12,16,24,32/
      DATA XG/ .43056815579702629d0 , .16999052179242813d0 , 14*0d0,
     X         .46623475710157601d0 , .33060469323313226d0 ,
     Y         .11930959304159845d0 ,                        13*0d0,
     Z         .48014492824876812d0 , .39833323870681337d0 ,
     1         .26276620495816449d0 , .9171732124782490 d-1, 12*0d0,
     2         .49078031712335963d0 , .45205862818523743d0 ,
     3         .38495133709715234d0 , .29365897714330872d0 ,
     4         .18391574949909010d0 , .62616704255734458d-1, 10*0d0,
     5         .49470046749582497d0 , .47228751153661629d0 ,
     6         .43281560119391587d0 , .37770220417750152d0 ,
     7         .30893812220132187d0 , .22900838882861369d0 ,
     8         .14080177538962946d0 , .47506254918818720d-1,  8*0d0,
     9         .49759360999851068d0 , .48736427798565475d0 ,
     A         .46913727600136638d0 , .44320776350220052d0 ,
     B         .41000099298695146d0 , .37006209578927718d0 ,
     C         .32404682596848778d0 , .27271073569441977d0 ,
     D         .21689675381302257d0 , .15752133984808169d0 ,
     E         .95559433736808150d-1, .32028446431302813d-1, 20*0d0/
      DATA YG/ .49863193092474078d0 , .49280575577263417d0 ,
     G         .48238112779375322d0 , .46745303796886984d0 ,
     H         .44816057788302606d0 , .42468380686628499d0 ,
     I         .39724189798397120d0 , .36609105937014484d0 ,
     J         .33152213346510760d0 , .29385787862038116d0 ,
     K         .25344995446611470d0 , .21067563806531767d0 ,
     L         .16593430114106382d0 , .11964368112606854d0 ,
     M         .72235980791398250d-1, .24153832843869158d-1/
      DATA DG/ .17392742256872693d0 , .32607257743127307d0 , 14*0d0,
     X         .85662246189585173d-1, .18038078652406930d0 ,
     Y         .23395696728634552d0 ,                        13*0d0,
     Z         .50614268145188130d-1, .11119051722668724d0 ,
     1         .15685332293894364d0 , .18134189168918099d0 , 12*0d0,
     2         .23587668193255914d-1, .53469662997659215d-1,
     3         .8003916427167311 d-1, .10158371336153296d0 ,
     4         .11674626826917740d0 , .12457352290670139d0 , 10*0d0,
     5         .13576229705877047d-1, .31126761969323946d-1,
     6         .47579255841246392d-1, .62314485627766936d-1,
     7         .74797994408288370d-1, .84578259697501270d-1,
     8         .91301707522461790d-1, .94725305227534250d-1,  8*0d0,
     9         .61706148999935998d-2, .14265694314466832d-1,
     A         .22138719408709903d-1, .29649292457718390d-1,
     B         .36673240705540153d-1, .43095080765976638d-1,
     C         .48809326052056944d-1, .53722135057982817d-1,
     D         .57752834026862801d-1, .60835236463901696d-1,
     E         .62918728173414148d-1, .63969097673376078d-1, 20*0d0/
      DATA EG/ .35093050047350483d-2, .8137197365452835 d-2,
     G         .12696032654631030d-1, .17136931456510717d-1,
     H         .21417949011113340d-1, .25499029631188088d-1,
     I         .29342046739267774d-1, .32911111388180923d-1,
     J         .36172897054424253d-1, .39096947893535153d-1,
     K         .41655962113473378d-1, .43826046502201906d-1,
     L         .45586939347881942d-1, .46922199540402283d-1,
     M         .47819360039637430d-1, .48270044257363900d-1/
C
      if(NBEG.EQ.0) then
      NBEG=1
      do 10 I=1,16
      XG(I,7)=YG(I)
   10 DG(I,7)=EG(I)
      end if
C
C     N MUST BE EVEN AND >= 4
C
      NN=(N/2)*2
      if(NN.LT.4) NN=4
      if(NN.NE.N) then
      WRITE (*,*) N,' GAUSS-POINTS WERE NOT POSSIBLE'
      N=NN
      WRITE (*,*) ' INSTEAD NOW ',N,' POINTS ARE USED'
      end if
      if(NA.NE.N) then
      NA=N
      NR=NA
      do 11 L=7,1,-1
      NI(L)=NR/NG(L)
      NR=NR-NG(L)*NI(L)
      if(NR.EQ.2) then
      NI(L)=NI(L)-1
      NR=NR+NG(L)
      end if
 11   CONTINUE
      if(NR.NE.0) WRITE (*,*) 'ERROR IN GAUSS: NR=',NR
      end if
C
      DELP=(B-A)/N
      if(DELP.EQ.0d0) GO TO 15
      A1=A
      I0=0
      do 14 L1=7,1,-1
      NIN=NI(L1)
      if(NIN.EQ.0) GO TO 14
      DEL=DELP*NG(L1)
      M   = NG(L1)/2
      do 13 K=1,NIN
      ABM =A1+DEL*0.5d0
      do 12 J=1,M
      I   = M+J
      L   = M+1-J
      J1=J+I0
      I1=I+I0
      X(J1)=ABM-DEL*XG(J,L1)
      X(I1)=ABM+DEL*XG(L,L1)
      W(J1)=    DEL*DG(J,L1)
   12 W(I1)=    DEL*DG(L,L1)
      I0=I0+NG(L1)
   13 A1=A1+DEL
   14 CONTINUE
C
C     TEST
C
      if(I0.NE.N) WRITE (*,*) 'ERROR IN GAUSS :',I0,'.NE.',N,
     +                        ' A1=',A1,' B=',B
      RETURN
   15 do 16 L=1,N
      X(L)=A
   16 W(L)=0d0
      RETURN
      end

c -------------------------------------------------------------------- c
c ---------------- Help FUNCTIONs for radiative decays --------------- c
c -------------------------------------------------------------------- c

c -------- Integrals needed in the radiative gluino decays:  --------- c
c -------- gluino -> neutralino_i gluon.                     --------- c
c -------- Taken from: Haber/Wyler Nucl.Phys.B323 (1989) 267 --------- c

      DOUBLE COMPLEX FUNCTION NS_iint(mj,mi,m,mc)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)

      DOUBLE COMPLEX NS_ccspen,lami,lamj,ctmp1,ctmp2,ctmp3,ctmp4,tmpa,
     .               tmpb,m2s,mc2s

      EXTERNAL NS_ccspen

      eps = 1d-10

      m2s  = dcmplx(m**2,-m**2*eps)
      mc2s = dcmplx(mc**2,-mc**2*eps)

      tmpa = dcmplx((m2s+mc2s-mi**2)**2-4d0*m2s*mc2s)
      tmpb = dcmplx((m2s+mc2s-mj**2)**2-4d0*m2s*mc2s)

      lami = cdsqrt(tmpa)
      lamj = cdsqrt(tmpb)

      ctmp1 = NS_ccspen( (mj**2+m2s-mc2s+lamj)/(2d0*m2s) )
      ctmp2 = NS_ccspen( (mj**2+m2s-mc2s-lamj)/(2d0*m2s) )
      ctmp3 = NS_ccspen( (mi**2+m2s-mc2s+lami)/(2d0*m2s) )
      ctmp4 = NS_ccspen( (mi**2+m2s-mc2s-lami)/(2d0*m2s) )

      if(mc.gt.1d4.or.m.gt.1d4) then
         NS_iint = 1d0/(mc2s-m2s)*(1d0-
     .                        mc2s/(mc2s-m2s)*cdlog(mc2s/m2s))
      else
         NS_iint = -1d0/(mj**2-mi**2)*(ctmp1+ctmp2-ctmp3-ctmp4)
      endif

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE COMPLEX FUNCTION NS_jint(mj,mi,m,mc)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)

      DOUBLE COMPLEX lami,lamj,tmpa,tmpb,ctmp1,ctmp2,ctmp3,ctmp4,m2s,
     .               mc2s,NS_iint

      EXTERNAL NS_iint

      eps = 1d-10

      m2s  = dcmplx(m**2,-m**2*eps)
      mc2s = dcmplx(mc**2,-mc**2*eps)

      tmpa = dcmplx((m2s+mc2s-mi**2)**2-4d0*m2s*mc2s)
      tmpb = dcmplx((m2s+mc2s-mj**2)**2-4d0*m2s*mc2s)

      lami = cdsqrt(tmpa)
      lamj = cdsqrt(tmpb)

      ctmp1 = dcmplx(m2s+mc2s-mj**2+lamj)
      ctmp2 = dcmplx(m2s+mc2s-mi**2+lami)

      ctmp3 = (cdlog(ctmp1/dcmplx(2d0*m*mc)))**2
      ctmp4 = (cdlog(ctmp2/dcmplx(2d0*m*mc)))**2

      if(mc.gt.1d4.or.m.gt.1d4) then
         NS_jint = 1d0/(mc2s-m2s)*cdlog(m2s/mc2s)
     .        - NS_iint(mj,mi,m,mc)
      else
         NS_jint = 1d0/(mj**2-mi**2)*( ctmp3 - ctmp4 )
     .        - NS_iint(mj,mi,m,mc)
      endif

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE COMPLEX FUNCTION NS_i2int(mj,mi,m,mc)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)

      DOUBLE COMPLEX lami,lamj,tmpa,tmpb,m2s,mc2s,ctmp3,ctmp4,ctmp5,
     .               ctmp6

      eps = 1d-10

      m2s  = dcmplx(m**2,-m**2*eps)
      mc2s = dcmplx(mc**2,-mc**2*eps)

      tmpa = dcmplx((m2s+mc2s-mi**2)**2-4d0*m2s*mc2s)
      tmpb = dcmplx((m2s+mc2s-mj**2)**2-4d0*m2s*mc2s)

      lami = cdsqrt(tmpa)
      lamj = cdsqrt(tmpb)

      ctmp3 = m2s+mc2s-dcmplx(mj**2)-lamj
      ctmp4 = m2s+mc2s-dcmplx(mj**2)+lamj
      ctmp5 = m2s+mc2s-dcmplx(mi**2)-lami
      ctmp6 = m2s+mc2s-dcmplx(mi**2)+lami

      if(mc.gt.1d4.or.m.gt.1d4) then
         NS_i2int = -1d0/(mc2s-m2s)**2*((m2s+mc2s)/2d0-
     .        m2s*mc2s/(mc2s-m2s)*cdlog(mc2s/m2s))
      else
         NS_i2int = (mc2s-m2s)/(2d0*mi**2*mj**2)*cdlog(m2s/mc2s) +
     .        1d0/(mj**2-mi**2)*(
     .        lamj/(2d0*mj**2)*cdlog(ctmp3/ctmp4) -
     .        lami/(2d0*mi**2)*cdlog(ctmp5/ctmp6) )
      endif

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE COMPLEX FUNCTION NS_kint(mj,mi,m,mc)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION m,mc,mj,mi
      DOUBLE COMPLEX NS_iint,NS_jint,NS_i2int,m2s,mc2s

      EXTERNAL NS_iint,NS_jint,NS_i2int

      eps = 1d-10

      m2s  = dcmplx(m**2,-m**2*eps)
      mc2s = dcmplx(mc**2,-mc**2*eps)

      if(mc.gt.1d4.or.m.gt.1d4) then
         NS_kint = 1d0/2d0*NS_i2int(mj,mi,m,mc)
      else
         NS_kint = -1d0/(mj**2-mi**2)*(1d0+m2s*NS_iint(mj,mi,m,mc)+
     .        mc2s*NS_jint(mj,mi,m,mc)-mj**2*NS_i2int(mj,mi,m,mc))
      endif

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE COMPLEX FUNCTION NS_iinthelp(mj,mi,m,mc)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE COMPLEX NS_ccspen,lami,lamj,ctmp1,ctmp2,ctmp3,ctmp4,tmpa,
     .               tmpb,m2s,mc2s

      EXTERNAL NS_ccspen

      eps = 1d-10

      m2s  = dcmplx(m**2,-m**2*eps)
      mc2s = dcmplx(mc**2,-mc**2*eps)

      tmpa = dcmplx((m2s+mc2s-mi**2)**2-4d0*m2s*mc2s)
      tmpb = dcmplx((m2s+mc2s-mj**2)**2-4d0*m2s*mc2s)

      lami = cdsqrt(tmpa)
      lamj = cdsqrt(tmpb)

      ctmp1 = NS_ccspen( (mj**2+m2s-mc2s+lamj)/(2d0*m2s) )
      ctmp2 = NS_ccspen( (mj**2+m2s-mc2s-lamj)/(2d0*m2s) )
      ctmp3 = NS_ccspen( (mi**2+m2s-mc2s+lami)/(2d0*m2s) )
      ctmp4 = NS_ccspen( (mi**2+m2s-mc2s-lami)/(2d0*m2s) )

      NS_iinthelp = -1d0/(mj**2-mi**2)*(ctmp1+ctmp2-ctmp3-ctmp4)

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE COMPLEX FUNCTION NS_jint0(mj,mi,mc)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE COMPLEX NS_iinthelp

      EXTERNAL NS_iinthelp

      NS_jint0 = NS_iinthelp(mj,mi,mc,0d0)

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE COMPLEX FUNCTION NS_i2int0(mj,mi,mc)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE COMPLEX mc2s

      eps = 1d-10

      mc2s = dcmplx(mc**2,-mc**2*eps)

      if(mc.gt.1d4) then
         NS_i2int0 = -1d0/2d0/mc2s
      else
         NS_i2int0 = -mc2s/mi**2/mj**2*cdlog(mc2s) + 1d0/(mj**2-mi**2)*
     .        ((mc2s-mi**2)/mi**2*cdlog(mc2s-mi**2) -
     .        (mc2s-mj**2)/mj**2*cdlog(mc2s-mj**2) )
      endif

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE COMPLEX FUNCTION NS_kint0(mj,mi,mc)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE COMPLEX NS_jint0,NS_i2int0,mc2s

      EXTERNAL NS_jint0,NS_i2int0

      eps = 1d-10

      mc2s = dcmplx(mc**2,-mc**2*eps)

      if(mc.gt.1d4) then
         NS_kint0 = 1d0/2d0*NS_i2int0(mj,mi,mc)
      else
         NS_kint0 = -1d0/(mj**2-mi**2)*(1d0+mc2s*NS_jint0(mj,mi,mc)
     .        -mj**2*NS_i2int0(mj,mi,mc))
      endif

      return

      end

c -------------------------------------------------------------------- c
c -------------- Help FUNCTIONs for the QCD corrections -------------- c
c -------------------------------------------------------------------- c

c -------------------------------------------------------------------- c
c ------------------ A.Bartl et al., hep-ph/9710286 ------------------ c
c -------------------------------------------------------------------- c
c -- The vertex corrections -- c

      DOUBLE PRECISION FUNCTION NS_gluonvertex(ami,amj,amv,lamv,amuv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_C1_lam,SD_C2_lam,SD_B02
      DOUBLE COMPLEX SD_C0_lam

      c1Den = - dreal(SD_C0_lam(ami,amj,amv,lamv)) -
     .    NS_C1_lam(ami,amj,amv,ami,lamv,amj,amuv,lamv)
      c2Den = SD_C2_lam(ami,amj,amv,ami,lamv,amj,amuv,lamv)

      NS_gluonvertex = SD_B02(ami**2,lamv,ami,amuv**2) +
     .     SD_B02(amj**2,lamv,amj,amuv**2) - 2d0*
     .     (ami**2+amj**2-amv**2)*
     .     (dreal(SD_C0_lam(ami,amj,amv,lamv)) + C1Den + C2Den)

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gluinoZvertex(ami,amj,amv,lamv,amuv,
     .     mgluino,amq,iq3L,eq,sw,thetasq)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION iq3L,NS_C1,SD_C2,SD_B02
      DOUBLE PRECISION DDCOS,DDSIN
      DOUBLE COMPLEX SD_C03

      c1Den = - dreal(SD_C03(ami**2,amv**2,amj**2,mgluino,amq,amq)) -
     .    NS_C1(ami,amv,amj,mgluino,amq,amq,amuv)
      c2Den = SD_C2(ami,amv,amj,mgluino,amq,amq,amuv)

      NS_gluinoZvertex = iq3L*(2d0*mgluino**2*
     .     dreal(SD_C03(ami**2,amv**2,amj**2,mgluino,amq,amq))
     .     + ami**2*C1Den + amj**2*C2Den
     .     +(mgluino**2-amq**2)*(C1Den+C2Den) +
     .     SD_B02(amv**2,amq,amq,amuv**2))*DDSIN(2d0*thetasq) +
     .     2d0*mgluino*amq*(iq3L-2d0*eq*sw**2)*
     .     ( dreal(SD_C03(ami**2,amv**2,amj**2,mgluino,amq,amq))
     .     + C1Den + C2Den )*DDCOS(2d0*thetasq)

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gluinoWvertex(ami,amj,amv,lamv,amuv,
     .     mgluino,amq,amqp,thsq,thsqp,i,j)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)
      DOUBLE PRECISION NS_C1,SD_C2,SD_B02
      DOUBLE PRECISION DDCOS,DDSIN
      DOUBLE COMPLEX SD_C03
      DIMENSION r(2,2),rp(2,2)

      r(1,1) = DDCOS(thsq)
      r(1,2) = DDSIN(thsq)
      r(2,1) = -DDSIN(thsq)
      r(2,2) = DDCOS(thsq)

      rp(1,1) = DDCOS(thsqp)
      rp(1,2) = DDSIN(thsqp)
      rp(2,1) = -DDSIN(thsqp)
      rp(2,2) = DDCOS(thsqp)

      c1Den = - dreal(SD_C03(ami**2,amv**2,amj**2,mgluino,amq,amqp)) -
     .    NS_C1(ami,amv,amj,mgluino,amq,amqp,amuv)
      c2Den = SD_C2(ami,amv,amj,mgluino,amq,amqp,amuv)

      NS_gluinoWvertex = mgluino*(
     .     dreal(SD_C03(ami**2,amv**2,amj**2,mgluino,amq,amqp)) +
     .     C1Den + C2Den )*
     .     ( amq*r(i,2)*rp(j,1) + amqp*r(i,1)*rp(j,2) )
     .     -(ami**2*C1Den + amj**2*C2Den + mgluino**2*(2d0*
     .       dreal(SD_C03(ami**2,amv**2,amj**2,mgluino,amq,amqp))
     .      + C1Den + C2Den )
     .      + SD_B02(amv**2,amq,amqp,amuv**2) )*r(i,1)*rp(j,1)
     .     -amq*amqp*( C1Den + C2Den )*r(i,2)*rp(j,2)

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_wavefuncvertex(amsqi,amsqpj,amq,amqp,
     .     thetasq,thetasqp,vecindex,quarkindex,ii,jj,mgluino,lamv,amuv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)
      DOUBLE PRECISION NS_deltaZnngluon,NS_deltaZnngluino,
     .         NS_deltaZnnpgluino,NS_deltaZnnpsquark

      DIMENSION gvqqp(2,2)
      DIMENSION gztt(2,2),gzbb(2,2),gztautau(2,2),gzmumu(2,2)
      DIMENSION gwtb(2,2),gwntau(2,2),gwnmu(2,2)
      DIMENSION gmst(2),gmsb(2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1

      COMMON/NS_coup19/gztt,gzbb,gztautau,gzmumu
      COMMON/NS_coup20/gwtb,gwntau,gwnmu
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1


      EXTERNAL NS_deltaZnngluon,NS_deltaZnngluino,NS_deltaZnnpgluino,
     .         NS_deltaZnnpsquark

      gmst(1) = ast1
      gmst(2) = ast2
      gmsb(1) = asb1
      gmsb(2) = asb2

      if(ii.eq.1) then
         ik = 2
      elseif(ii.eq.2) then
         ik = 1
      endif

      if(jj.eq.1) then
         il = 2
      elseif(jj.eq.2) then
         il = 1
      endif

      if(quarkindex.eq.1d0) then
         if(vecindex.eq.1d0) then
            do i=1,2
               do j=1,2
                  gvqqp(i,j) = 1d0/2d0/cw*gztt(i,j)
                  amsqk  = gmst(ik)
                  amsqpl = gmst(il)
                  amsq1  = gmst(1)
                  amsq2  = gmst(2)
                  amsqp1 = gmst(1)
                  amsqp2 = gmst(2)
               enddo
            enddo
         elseif(vecindex.eq.2d0) then
            do i=1,2
               do j=1,2
                  gvqqp(i,j) = 1d0/dsqrt(2d0)*gwtb(i,j)
                  amsqk  = gmst(ik)
                  amsqpl = gmsb(il)
                  amsq1  = gmst(1)
                  amsq2  = gmst(2)
                  amsqp1 = gmsb(1)
                  amsqp2 = gmsb(2)
               enddo
            enddo
         endif
      elseif(quarkindex.eq.2d0) then
         if(vecindex.eq.1d0) then
            do i=1,2
               do j=1,2
                  gvqqp(i,j) = 1d0/2d0/cw*gzbb(i,j)
                  amsqk  = gmsb(ik)
                  amsqpl = gmsb(il)
                  amsq1  = gmsb(1)
                  amsq2  = gmsb(2)
                  amsqp1 = gmsb(1)
                  amsqp2 = gmsb(2)
               enddo
            enddo
         elseif(vecindex.eq.2d0) then
            do i=1,2
               do j=1,2
                  gvqqp(i,j) = 1d0/dsqrt(2d0)*gwtb(j,i)
                  amsqk  = gmsb(ik)
                  amsqpl = gmst(il)
                  amsq1  = gmsb(1)
                  amsq2  = gmsb(2)
                  amsqp1 = gmst(1)
                  amsqp2 = gmst(2)
               enddo
            enddo
         endif
      endif

      NS_wavefuncvertex =1d0/2d0*(NS_deltaZnngluon(amsqi,lamv,amuv)+
     .                NS_deltaZnngluon(amsqpj,lamv,amuv))*gvqqp(ii,jj)
     .  +  1d0/2d0*(NS_deltaZnngluino(amsqi,mgluino,amq,amuv,dble(ii),
     .                thetasq) +
     .                NS_deltaZnngluino(amsqpj,mgluino,amqp,amuv,
     .                dble(jj),thetasqp) )*gvqqp(ii,jj)
     .    +NS_deltaZnnpgluino(amsqi,amsqk,mgluino,amq,amuv,thetasq)*
     .     gvqqp(ik,jj)
     .    +NS_deltaZnnpgluino(amsqpj,amsqpl,mgluino,amqp,amuv,thetasqp)*
     .     gvqqp(ii,il)
     .    +NS_deltaZnnpsquark(amsq1,amsq2,amsqi,amsqk,amuv,thetasq)*
     .     gvqqp(ik,jj)
     .    +NS_deltaZnnpsquark(amsqp1,amsqp2,amsqpj,amsqpl,amuv,thetasqp)
     .     *gvqqp(ii,il)

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_quarkmixZ(amsq,thetasq,iq3L,eq,amsq1,
     .     amsq2,amq,mgluino,amuv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION iq3L,NS_A01,SD_B02
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/NS_weinberg/sw,cw,tw

      deltathetasqsq = 1d0/6d0*DDSIN(4d0*thetasq)/(amsq1**2-amsq2**2)
     .     *( NS_A01(amsq2**2,amuv**2) - NS_A01(amsq1**2,amuv**2) )

      v11 = 4d0*(iq3L*DDCOS(thetasq)**2-eq*sw**2)
      v22 = 4d0*(iq3L*DDSIN(thetasq)**2-eq*sw**2)

      deltathetasqgl = 1d0/3d0*mgluino*amq/iq3L/(amsq1**2-amsq2**2)*
     .     ( SD_B02(amsq2**2,mgluino,amq,amuv**2)*v11 -
     .       SD_B02(amsq1**2,mgluino,amq,amuv**2)*v22 )

      NS_quarkmixZ = -1d0/cw*iq3L*DDCOS(2d0*thetasq)*(
     .     deltathetasqsq + deltathetasqgl )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_quarkmixW(amsq,thetasq,iq3L,eq,amsq1,
     .     amsq2,amq,amsqp,thetasqp,iqp3L,eqp,amsqp1,amsqp2,amqp,ii,jj,
     .     mgluino,amuv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)
      DOUBLE PRECISION iq3L,iqp3L,NS_A01,SD_B02
      DOUBLE PRECISION DDCOS,DDSIN
      DIMENSION cijwsq(2,2),cijwsqp(2,2)
      COMMON/NS_weinberg/sw,cw,tw


      cijwsq(1,1) = -DDSIN(thetasq)*DDCOS(thetasqp)
      cijwsq(1,2) =  DDSIN(thetasq)*DDSIN(thetasqp)
      cijwsq(2,1) = -DDCOS(thetasq)*DDCOS(thetasqp)
      cijwsq(2,2) =  DDCOS(thetasq)*DDSIN(thetasqp)

      cijwsqp(1,1) = -DDSIN(thetasqp)*DDCOS(thetasq)
      cijwsqp(1,2) = -DDCOS(thetasq)*DDCOS(thetasqp)
      cijwsqp(2,1) =  DDSIN(thetasq)*DDSIN(thetasqp)
      cijwsqp(2,2) =  DDCOS(thetasqp)*DDSIN(thetasq)

c --------------------------------------

      deltathetasqsq = 1d0/6d0*DDSIN(4d0*thetasq)/(amsq1**2-amsq2**2)
     .     *( NS_A01(amsq2**2,amuv**2) - NS_A01(amsq1**2,amuv**2) )

      deltathetasqpsq = 1d0/6d0*DDSIN(4d0*thetasqp)
     .     /(amsqp1**2-amsqp2**2)*
     .     ( NS_A01(amsqp2**2,amuv**2) - NS_A01(amsqp1**2,amuv**2) )

c --------------------------------------

      v11sq = 4d0*(iq3L*DDCOS(thetasq)**2-eq*sw**2)
      v22sq = 4d0*(iq3L*DDSIN(thetasq)**2-eq*sw**2)

      deltathetasqgl = 1d0/3d0*mgluino*amq/iq3L/(amsq1**2-amsq2**2)*
     .     ( SD_B02(amsq2**2,mgluino,amq,amuv**2)*v11sq -
     .       SD_B02(amsq1**2,mgluino,amq,amuv**2)*v22sq )

c --------------------------------------

      v11sqp = 4d0*(iqp3L*DDCOS(thetasqp)**2-eqp*sw**2)
      v22sqp = 4d0*(iqp3L*DDSIN(thetasqp)**2-eqp*sw**2)

      deltathetasqpgl = 1d0/3d0*mgluino*amqp/iqp3L
     .     /(amsqp1**2-amsqp2**2)*
     .     ( SD_B02(amsqp2**2,mgluino,amqp,amuv**2)*v11sqp -
     .       SD_B02(amsqp1**2,mgluino,amqp,amuv**2)*v22sqp )

c --------------------------------------

      NS_quarkmixW = 1d0/dsqrt(2d0)*(
     .     cijwsq(ii,jj)*(deltathetasqsq+deltathetasqgl) +
     .     cijwsqp(ii,jj)*(deltathetasqpsq+deltathetasqpgl) )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_realgluonem(amsqi,amsqpj,amv,lamv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)

      DOUBLE COMPLEX NS_ccspen,NS_kappa

      EXTERNAL NS_kappa,NS_ccspen

      kap = dreal(NS_kappa(amsqi**2,amsqpj**2,amv**2,0d0))

      b0 = (amsqi**2-amsqpj**2-amv**2+kap)/2d0/amsqpj/amv
      b1 = (amsqi**2-amsqpj**2+amv**2-kap)/2d0/amsqi/amv
      b2 = (amsqi**2+amsqpj**2-amv**2-kap)/2d0/amsqi/amsqpj

      lb0 = dreal(cdlog(dcmplx(b0)))
      lb1 = dreal(cdlog(dcmplx(b1)))
      lb2 = dreal(cdlog(dcmplx(b2)))
      lb12 = dreal(cdlog(dcmplx(b1/b2)))
      lb02 = dreal(cdlog(dcmplx(b0/b2)))

      hint = 1d0/4d0/amsqi**2*(kap/2d0*(amsqi**2+amsqpj**2+amv**2)+
     .     2d0*amsqi**2*amsqpj**2*lb2 +
     .     2d0*amsqi**2*amv**2*lb1 +
     .     2d0*amsqpj**2*amv**2*lb0 )

      hint0 = 1d0/4d0/amsqi**2*(-2d0*amsqpj**2*lb2-2d0*amv**2*lb1-
     .     kap)

      hint1 = 1d0/4d0/amsqi**2*(-2d0*amsqi**2*lb2-2d0*amv**2*lb0-
     .     kap)

      hint00 = 1d0/4d0/amsqi**4*(
     .     kap*dlog(kap**2/(lamv*amsqi*amsqpj*amv)) - kap -
     .     (amsqpj**2-amv**2)*lb12 - amsqi**2*lb0 )

      hint11 = 1d0/4d0/amsqi**2/amsqpj**2*(
     .     kap*dlog(kap**2/(lamv*amsqi*amsqpj*amv)) - kap -
     .     (amsqi**2-amv**2)*lb02 - amsqpj**2*lb1 )

      hint01 = dreal(1d0/4d0/amsqi**2*(
     .     -2d0*dlog((lamv*amsqi*amsqpj*amv)/kap**2)*lb2 +
     .     2d0*lb2**2 - lb0**2 - lb1**2 +
     .     2d0*NS_ccspen(dcmplx(1d0-b2**2)) -
     .     NS_ccspen(dcmplx(1-b0**2))
     .     - NS_ccspen(dcmplx(1-b1**2)) ) )

      NS_realgluonem = 2d0*hint - kap**2/amv**2*( hint0+hint1+
     .     amsqi**2*hint00+amsqpj**2*hint11+
     .     (amsqi**2+amsqpj**2-amv**2)*hint01 )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_deltaZnngluon(amsq,lamv,amuv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION SD_B02,SD_BP02

      NS_deltaZnngluon = 2d0/3d0*(
     .     SD_B02(amsq**2,0d0,amsq,amuv**2) +
     .     2d0*amsq**2*SD_BP02(amsq**2,lamv,amsq,amuv**2) )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_deltaZnngluino(amsq,mgluino,amq,
     .     amuv,n,thetasq)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION SD_B02,SD_BP02
      DOUBLE PRECISION DDCOS,DDSIN

      NS_deltaZnngluino = -2d0/3d0*(
     .     SD_B02(amsq**2,mgluino,amq,amuv**2) +
     .     (amsq**2-amq**2-mgluino**2)*
     .     SD_BP02(amsq**2,mgluino,amq,amuv**2) -
     .     2d0*amq*mgluino*(-1d0)**n*DDSIN(2d0*thetasq)*
     .     SD_BP02(amsq**2,mgluino,amq,amuv**2) )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_deltaZnnpgluino(amsq,amsqp,mgluino,
     .     amq,amuv,thetasq)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION SD_B02
      DOUBLE PRECISION DDCOS,DDSIN

      NS_deltaZnnpgluino = - 1d0/(amsq**2-amsqp**2)*4d0/3d0*
     .     mgluino*amq*DDCOS(2d0*thetasq)*
     .     SD_B02(amsq**2,mgluino,amq,amuv**2)

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_deltaZnnpsquark(amsq1,amsq2,amsq,
     .     amsqp,amuv,thetasq)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_A01
      DOUBLE PRECISION DDCOS,DDSIN

      NS_deltaZnnpsquark = - 1d0/(amsq**2-amsqp**2)*1d0/6d0*
     .     DDSIN(4d0*thetasq)*( NS_A01(amsq2**2,amuv**2) -
     .     NS_A01(amsq1**2,amuv**2) )

      return

      end

c -------------------------------------------------------------------- c
c ----- A.Arhrib, A.Djouadi, W.Hollik, C.Juenger, hep-ph/9702426 ----- c
c -------------------------------------------------------------------- c
c -- Virtual corrections for the decays squark_i -> Higgs squark_j' -- c

      DOUBLE PRECISION FUNCTION NS_gvirtgl(ami,amhi,amj,lamv,amuv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION SD_B02
      DOUBLE COMPLEX SD_C0_lam

      amisq = ami**2
      amjsq = amj**2
      amhsq = amhi**2
      amusq = amuv**2

      NS_gvirtgl = SD_B02(amisq,lamv,ami,amusq) +
     .             SD_B02(amjsq,lamv,amj,amusq) -
     .             SD_B02(amhsq,ami,amj,amusq)  +
     .             2d0*(amisq+amjsq-amhsq)*
     .             dreal(SD_C0_lam(amj,ami,amhi,lamv))

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gvirtmix(am1,am2,amij,amgl,amq,theq,
     .                                   amuv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_A01,SD_B02
      DOUBLE PRECISION DDCOS,DDSIN

      c2q = DDCOS(2d0*theq)
      s2q = DDSIN(2d0*theq)

      am1sq  = am1**2
      am2sq  = am2**2
      amijsq = amij**2
      amusq  = amuv**2

      NS_gvirtmix = c2q*s2q*(NS_A01(am2sq,amusq)-
     .           NS_A01(am1sq,amusq))+
     .           4d0*c2q*amq*amgl*SD_B02(amijsq,amgl,amq,amusq)

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gvirtmixdiv(am1,am2,amij,amgl,amq,
     .                                         theq,amuv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION SD_B02_DIV
      DOUBLE PRECISION DDCOS,DDSIN

      c2q = DDCOS(2d0*theq)
      s2q = DDSIN(2d0*theq)

      am1sq  = am1**2
      am2sq  = am2**2
      amijsq = amij**2
      amusq  = amuv**2

      NS_gvirtmixdiv = c2q*s2q*(am2sq*dlog(amusq)-am1sq*dlog(amusq) )
     .     +4d0*c2q*amq*amgl*SD_B02_DIV(amijsq,amgl,amq,amusq)

      return

      end

c -------------------------------------------------------------------- c
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
c -------------------------------------------------------------------- c
C REF:::PRD57,5860
c -------------------------------------------------------------------- c
C**********************************************************************C
C**********************************************************************C
      DOUBLE PRECISION FUNCTION NS_topneut1719(nh,amuv)
*
      IMPLICIT NONE
*
      INTEGER nh
      DOUBLE COMPLEX SD_C03
*
      DOUBLE PRECISION m(2,2,2)
      DOUBLE PRECISION amusq,amuv,aml,amh,amnh,ama,amna,
     .coupphi,squarktopneut,amq,amar,amch,
     .ast1sq,ast2sq,amlsq,amhsq,amnhsq,amasq,amnasq,v1,v2,a1,a2,gluinoex
      DOUBLE PRECISION sdthet,sdtheb,sdthel,sdthem,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION Hstopstopr(3,2,2),Astopstopr(3,2,2)
      DOUBLE PRECISION Httr(3),Attr(2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION SD_B02
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION runmt,runmb,rmtauc
*
      COMMON/NS_sfmixang/sdthet,sdtheb,sdthel,sdthem,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_HIGGSSTST/Hstopstopr,Astopstopr
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1

      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_phitoptop/Httr,Attr
*
      CALL NS_smatrix(m)
*
      aml = MAX(1d0,SMASS(1))
      amh = SMASS(2)
      amnh = SMASS(3)
      ama = MAX(1d0,PMASS(1))
      amar = ama
      amna = PMASS(2)
      amch = CMASS
      amusq  = amuv**2
      amq    = runmt
      ast1sq = ast1**2
      ast2sq = ast2**2
      amlsq  = aml**2
      amhsq  = amh**2
      amnhsq = amnh**2
      amasq  = ama**2
      amnasq = amna**2
*
      v1 = 1d0/2d0*(ct-st)
      v2 = -1d0/2d0*(ct+st)
      a1 = -v2
      a2 = v1
C------------------EQ.17------------------------
      if(nh.eq.1) then
         coupphi = m(2,1,1)*Hstopstopr(1,1,1)*m(1,1,1)*
     .        SD_B02(aml**2,ast1,ast1,amusq) +
     .        m(2,1,1)*Hstopstopr(1,1,2)*m(2,1,1)*
     .        SD_B02(aml**2,ast1,ast2,amusq) +
     .        m(2,2,1)*Hstopstopr(1,1,2)*m(1,1,1)*
     .        SD_B02(aml**2,ast2,ast1,amusq) +
     .        m(2,2,1)*Hstopstopr(1,2,2)*m(2,1,1)*
     .        SD_B02(aml**2,ast2,ast2,amusq)
      elseif(nh.eq.2) then
         coupphi = m(2,1,1)*Hstopstopr(2,1,1)*m(1,1,1)*
     .        SD_B02(amh**2,ast1,ast1,amusq) +
     .        m(2,1,1)*Hstopstopr(2,1,2)*m(2,1,1)*
     .        SD_B02(amh**2,ast1,ast2,amusq) +
     .        m(2,2,1)*Hstopstopr(2,1,2)*m(1,1,1)*
     .        SD_B02(amh**2,ast2,ast1,amusq) +
     .        m(2,2,1)*Hstopstopr(2,2,2)*m(2,1,1)*
     .        SD_B02(amh**2,ast2,ast2,amusq)
       elseif(nh.eq.3) then
         coupphi = m(2,1,1)*Hstopstopr(3,1,1)*m(1,1,1)*
     .        SD_B02(amnh**2,ast1,ast1,amusq) +
     .        m(2,1,1)*Hstopstopr(3,1,2)*m(2,1,1)*
     .        SD_B02(amnh**2,ast1,ast2,amusq) +
     .        m(2,2,1)*Hstopstopr(3,1,2)*m(1,1,1)*
     .        SD_B02(amnh**2,ast2,ast1,amusq) +
     .        m(2,2,1)*Hstopstopr(3,2,2)*m(2,1,1)*
     .        SD_B02(amnh**2,ast2,ast2,amusq)
       elseif(nh.eq.4) then
         coupphi = m(2,1,1)*Astopstopr(1,1,1)*m(1,1,1)*
     .        SD_B02(ama**2,ast1,ast2,amusq) -
     .        m(2,1,1)*Astopstopr(1,1,2)*m(2,1,1)*
     .        SD_B02(ama**2,ast1,ast2,amusq)-
     .        m(2,2,1)*Astopstopr(1,1,2)*m(1,1,1)*
     .        SD_B02(ama**2,ast2,ast1,amusq)-
     .        m(2,2,1)*Astopstopr(1,2,2)*m(2,1,1)*
     .        SD_B02(ama**2,ast2,ast2,amusq)
      elseif(nh.eq.5) then
         coupphi = m(2,1,1)*Astopstopr(2,1,1)*m(1,1,1)*
     .        SD_B02(amna**2,ast1,ast2,amusq) -
     .        m(2,1,1)*Astopstopr(2,1,2)*m(2,1,1)*
     .        SD_B02(amna**2,ast1,ast2,amusq)-
     .        m(2,2,1)*Astopstopr(2,1,2)*m(1,1,1)*
     .        SD_B02(amna**2,ast2,ast1,amusq)-
     .        m(2,2,1)*Astopstopr(2,2,2)*m(2,1,1)*
     .        SD_B02(amna**2,ast2,ast2,amusq)
      endif

      coupphi = dsqrt(2d0)*AMZ**2*coupphi

c the result Eq.(17) in the paper

            squarktopneut = coupphi
c --------------------------------------------
c the result Eq.(19) in the paper

      if(nh.eq.1) then
         gluinoex = 4d0*AMW*Httr(1)*(
     .        (amq*(v2*v1+a2*a1)+mgluino*(v2*v1-a2*a1))*
     .        (SD_B02(ast2**2,mgluino,amq,amusq)+
     .         SD_B02(ast1**2,mgluino,amq,amusq)) +
     .        2d0*amq*(v2*v1+a2*a1)*
     .        SD_B02(aml**2,amq,amq,amusq) +
     .        (-mgluino*(v2*v1-a2*a1)*(aml**2-4d0*amq**2)
     .         -(v2*v1+a2*a1)*(ast1**2*amq+ast2**2*amq-(mgluino**2+
     .         amq**2)*2d0*amq) )*
     .        dreal(SD_C03(ast2sq,amlsq,ast1sq,mgluino,amq,amq)) )
      elseif(nh.eq.2) then
         gluinoex = 4d0*AMW*Httr(2)*(
     .        (amq*(v2*v1+a2*a1)+mgluino*(v2*v1-a2*a1))*
     .        (SD_B02(ast2**2,mgluino,amq,amusq)+
     .         SD_B02(ast1**2,mgluino,amq,amusq)) +
     .        2d0*amq*(v2*v1+a2*a1)*
     .        SD_B02(amh**2,amq,amq,amusq) +
     .        (-mgluino*(v2*v1-a2*a1)*(amh**2-4d0*amq**2)
     .         -(v2*v1+a2*a1)*(ast1**2*amq+ast2**2*amq-(mgluino**2+
     .         amq**2)*2d0*amq) )*
     .        dreal(SD_C03(ast2sq,amhsq,ast1sq,mgluino,amq,amq)) )
          elseif(nh.eq.3) then
         gluinoex = 4d0*AMW*Httr(3)*(
     .        (amq*(v2*v1+a2*a1)+mgluino*(v2*v1-a2*a1))*
     .        (SD_B02(ast2**2,mgluino,amq,amusq)+
     .         SD_B02(ast1**2,mgluino,amq,amusq)) +
     .        2d0*amq*(v2*v1+a2*a1)*
     .        SD_B02(amnh**2,amq,amq,amusq) +
     .        (-mgluino*(v2*v1-a2*a1)*(amnh**2-4d0*amq**2)
     .         -(v2*v1+a2*a1)*(ast1**2*amq+ast2**2*amq-(mgluino**2+
     .         amq**2)*2d0*amq) )*
     .        dreal(SD_C03(ast2sq,amnhsq,ast1sq,mgluino,amq,amq)) )
         elseif(nh.eq.4) then
         gluinoex = -4d0*AMW*Attr(1)*(
     .        (-amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(ast2**2,mgluino,amq,amusq)+
     .        (amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(ast1**2,mgluino,amq,amusq) +
     .        (mgluino*(v2*a1-a2*v1)*ama**2
     .         +(v2*a1+a2*v1)*(ast1**2*amq-ast2**2*amq) )*
     .        dreal(SD_C03(ast2sq,amasq,ast1sq,mgluino,amq,amq)) )
      elseif(nh.eq.5) then
         gluinoex = -4d0*AMW*Attr(2)*(
     .        (-amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(ast2**2,mgluino,amq,amusq)+
     .        (amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(ast1**2,mgluino,amq,amusq) +
     .        (mgluino*(v2*a1-a2*v1)*amna**2
     .         +(v2*a1+a2*v1)*(ast1**2*amq-ast2**2*amq) )*
     .        dreal(SD_C03(ast2sq,amnasq,ast1sq,mgluino,amq,amq)) )
      endif

      NS_topneut1719 = squarktopneut + gluinoex

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_botneut1719(nh,amuv)
      IMPLICIT NONE
*
      INTEGER nh
      DOUBLE COMPLEX SD_C03
*
      DOUBLE PRECISION m(2,2,2)
      DOUBLE PRECISION amusq,amuv,aml,amh,amnh,ama,amna,
     .coupphi,squarkbotneut,amq,amar,amch,
     .asb1sq,asb2sq,amlsq,amhsq,amnhsq,amasq,amnasq,v1,v2,a1,a2,gluinoex
      DOUBLE PRECISION sdthet,sdtheb,sdthel,sdthem,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION Hsbotsbotr(3,2,2),Asbotsbotr(2,2,2)
      DOUBLE PRECISION Hbbr(3),Abbr(2)
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION SD_B02
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
***
      COMMON/NS_sfmixang/sdthet,sdtheb,sdthel,sdthem,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_HIGGSBTBT/Hsbotsbotr,Asbotsbotr
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_phibotbot/Hbbr,Abbr
*
      aml = MAX(1d0,SMASS(1))
      amh = SMASS(2)
      amnh = SMASS(3)
      ama = MAX(1d0,PMASS(1))
      amar = ama
      amna = PMASS(2)
      amch = CMASS
      CALL NS_smatrix(m)
      amusq  = amuv**2
c -----------------------------

      if(nh.eq.1) then
         coupphi = m(2,1,2)*Hsbotsbotr(1,1,1)*m(1,1,2)*
     .             SD_B02(aml**2,asb1,asb1,amusq) +
     .             m(2,1,2)*Hsbotsbotr(1,1,2)*m(2,1,2)*
     .             SD_B02(aml**2,asb1,asb2,amusq) +
     .             m(2,2,2)*Hsbotsbotr(1,1,2)*m(1,1,2)*
     .             SD_B02(aml**2,asb2,asb1,amusq) +
     .             m(2,2,2)*Hsbotsbotr(1,2,2)*m(2,1,2)*
     .             SD_B02(aml**2,asb2,asb2,amusq)
      elseif(nh.eq.2) then
         coupphi = m(2,1,2)*Hsbotsbotr(2,1,1)*m(1,1,2)*
     .             SD_B02(amh**2,asb1,asb1,amusq) +
     .             m(2,1,2)*Hsbotsbotr(2,1,2)*m(2,1,2)*
     .             SD_B02(amh**2,asb1,asb2,amusq) +
     .             m(2,2,2)*Hsbotsbotr(2,1,2)*m(1,1,2)*
     .             SD_B02(amh**2,asb2,asb1,amusq) +
     .             m(2,2,2)*Hsbotsbotr(2,2,2)*m(2,1,2)*
     .             SD_B02(amh**2,asb2,asb2,amusq)
      elseif(nh.eq.3) then
         coupphi = m(2,1,2)*Hsbotsbotr(3,1,1)*m(1,1,2)*
     .             SD_B02(amnh**2,asb1,asb1,amusq) +
     .             m(2,1,2)*Hsbotsbotr(3,1,2)*m(2,1,2)*
     .             SD_B02(amnh**2,asb1,asb2,amusq) +
     .             m(2,2,2)*Hsbotsbotr(3,1,2)*m(1,1,2)*
     .             SD_B02(amnh**2,asb2,asb1,amusq) +
     .             m(2,2,2)*Hsbotsbotr(3,2,2)*m(2,1,2)*
     .             SD_B02(amnh**2,asb2,asb2,amusq)
      elseif(nh.eq.4) then
         coupphi = m(2,1,2)*Asbotsbotr(1,1,1)*m(1,1,2)*
     .             SD_B02(ama**2,asb1,asb2,amusq) -
     .              m(2,1,2)*Asbotsbotr(1,1,2)*m(2,1,2)*
     .             SD_B02(ama**2,asb1,asb2,amusq)-
     .             m(2,2,2)*Asbotsbotr(1,1,2)*m(1,1,2)*
     .             SD_B02(ama**2,asb2,asb1,amusq)-
     .             m(2,2,2)*Asbotsbotr(1,2,2)*m(2,1,2)*
     .             SD_B02(ama**2,asb2,asb1,amusq)
      elseif(nh.eq.5) then
         coupphi = m(2,1,2)*Asbotsbotr(2,1,1)*m(1,1,2)*
     .             SD_B02(amna**2,asb1,asb2,amusq) -
     .              m(2,1,2)*Asbotsbotr(2,1,2)*m(2,1,2)*
     .             SD_B02(amna**2,asb1,asb2,amusq)-
     .             m(2,2,2)*Asbotsbotr(2,1,2)*m(1,1,2)*
     .             SD_B02(amna**2,asb2,asb1,amusq)-
     .             m(2,2,2)*Asbotsbotr(2,2,2)*m(2,1,2)*
     .             SD_B02(amna**2,asb2,asb1,amusq)

      endif

      coupphi = dsqrt(2d0)*AMZ**2*coupphi

c the result Eq.(17) in the paper

      squarkbotneut = coupphi
c --------------------------------------------
      amq    = runmb
      asb1sq = asb1**2
      asb2sq = asb2**2
      amlsq  = aml**2
      amhsq  = amh**2
      amnhsq  = amnh**2
      amasq  = ama**2
      amnasq  = amna**2
**
      v1 = 1d0/2d0*(cb-sb)
      v2 = -1d0/2d0*(cb+sb)
      a1 = -v2
      a2 = v1

c the result Eq.(19) in the paper

      if(nh.eq.1) then
         gluinoex = 4d0*AMW*Hbbr(1)*(
     .        (amq*(v2*v1+a2*a1)+mgluino*(v2*v1-a2*a1))*
     .        (SD_B02(asb2**2,mgluino,amq,amusq)+
     .         SD_B02(asb1**2,mgluino,amq,amusq)) +
     .        2d0*amq*(v2*v1+a2*a1)*
     .        SD_B02(aml**2,amq,amq,amusq) +
     .        (-mgluino*(v2*v1-a2*a1)*(aml**2-4d0*amq**2)
     .         -(v2*v1+a2*a1)*(asb1**2*amq+asb2**2*amq-(mgluino**2+
     .         amq**2)*2d0*amq) )*
     .        dreal(SD_C03(asb2sq,amlsq,asb1sq,mgluino,amq,amq)) )
      elseif(nh.eq.2) then
         gluinoex = 4d0*AMW*Hbbr(2)*(
     .        (amq*(v2*v1+a2*a1)+mgluino*(v2*v1-a2*a1))*
     .        (SD_B02(asb2**2,mgluino,amq,amusq)+
     .         SD_B02(asb1**2,mgluino,amq,amusq)) +
     .        2d0*amq*(v2*v1+a2*a1)*
     .        SD_B02(amh**2,amq,amq,amusq) +
     .        (-mgluino*(v2*v1-a2*a1)*(amh**2-4d0*amq**2)
     .         -(v2*v1+a2*a1)*(asb1**2*amq+asb2**2*amq-(mgluino**2+
     .         amq**2)*2d0*amq) )*
     .        dreal(SD_C03(asb2sq,amhsq,asb1sq,mgluino,amq,amq)) )
      elseif(nh.eq.3) then
         gluinoex = 4d0*AMW*Hbbr(3)*(
     .        (amq*(v2*v1+a2*a1)+mgluino*(v2*v1-a2*a1))*
     .        (SD_B02(asb2**2,mgluino,amq,amusq)+
     .         SD_B02(asb1**2,mgluino,amq,amusq)) +
     .        2d0*amq*(v2*v1+a2*a1)*
     .        SD_B02(amnh**2,amq,amq,amusq) +
     .        (-mgluino*(v2*v1-a2*a1)*(amnh**2-4d0*amq**2)
     .         -(v2*v1+a2*a1)*(asb1**2*amq+asb2**2*amq-(mgluino**2+
     .         amq**2)*2d0*amq) )*
     .        dreal(SD_C03(asb2sq,amnhsq,asb1sq,mgluino,amq,amq)) )
               elseif(nh.eq.4) then
         gluinoex = -4d0*AMW*Abbr(1)*(
     .        (-amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(asb2**2,mgluino,amq,amusq)+
     .        (amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(asb1**2,mgluino,amq,amusq) +
     .        (mgluino*(v2*a1-a2*v1)*ama**2
     .         +(v2*a1+a2*v1)*(asb1**2*amq-asb2**2*amq) )*
     .        dreal(SD_C03(asb2sq,amasq,asb1sq,mgluino,amq,amq)) )
      elseif(nh.eq.5) then
         gluinoex = -4d0*AMW*Abbr(2)*(
     .        (-amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(asb2**2,mgluino,amq,amusq)+
     .        (amq*(v2*a1+a2*v1)-mgluino*(v2*a1-a2*v1))*
     .         SD_B02(asb1**2,mgluino,amq,amusq) +
     .        (mgluino*(v2*a1-a2*v1)*amna**2
     .         +(v2*a1+a2*v1)*(asb1**2*amq-asb2**2*amq) )*
     .        dreal(SD_C03(asb2sq,amnasq,asb1sq,mgluino,amq,amq)) )
      endif

      NS_botneut1719 = squarkbotneut + gluinoex

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_stopsbot1719(amuv,ni,nj)
*
      IMPLICIT NONE
      DOUBLE COMPLEX SD_C03
*
      DOUBLE PRECISION m(2,2,2),gctbr(2,2),gmst(2),vq(2),aq(2),vqp(2),
     .          aqp(2),gmsb(2)
      INTEGER ni,nj
      DOUBLE PRECISION amusq,amuv,
     .coupphi,squarkstopsbot,gluinoex,amq,amqp,amchsq,vs,as
      DOUBLE PRECISION sdthet,sdtheb,sdthel,sdthem,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION chtbrunr,chtbrunl
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION ama,aml,amh,amch,amar,amna,amnh
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION SD_B02
*
      COMMON/NS_sfmixang/sdthet,sdtheb,sdthel,sdthem,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      COMMON/NS_hcsbotstop/gctbr
      COMMON/NS_higgschudb/chtbrunr,chtbrunl
**
      aml = MAX(1d0,SMASS(1))
      amh = SMASS(2)
      amnh = SMASS(3)
      ama = MAX(1d0,PMASS(1))
      amar = ama
      amna = PMASS(2)
      amch = CMASS
*
      gmst(1) = ast1
      gmst(2) = ast2
      gmsb(1) = asb1
      gmsb(2) = asb2
*
      CALL NS_smatrix(m)
      amusq = amuv**2
c  --------------

      coupphi = m(ni,1,1)*gctbr(1,1)*m(1,nj,2)*
     .     SD_B02(amch**2,ast1,asb1,amusq) +
     .     m(ni,1,1)*gctbr(1,2)*m(2,nj,2)*
     .     SD_B02(amch**2,ast1,asb2,amusq) +
     .     m(ni,2,1)*gctbr(2,1)*m(1,nj,2)*
     .     SD_B02(amch**2,ast2,asb1,amusq) +
     .     m(ni,2,1)*gctbr(2,2)*m(2,nj,2)*
     .     SD_B02(amch**2,ast2,asb2,amusq)

      coupphi = dsqrt(2d0)*AMW**2*coupphi

c the result Eq.(17) in the paper

      squarkstopsbot = coupphi

c --------------------------------------------

      amq  =  runmt
      amqp =  runmb

      amchsq = amch**2

      vq(1) = 1d0/2d0*(ct-st)
      vq(2) = -1d0/2d0*(ct+st)
      aq(1) = -vq(2)
      aq(2) = vq(1)

      vqp(1) = 1d0/2d0*(cb-sb)
      vqp(2) = -1d0/2d0*(cb+sb)
      aqp(1) = -vqp(2)
      aqp(2) = vqp(1)

      vs = 2d0*dsqrt(2d0)*AMW*(chtbrunr+chtbrunl)
      as = 2d0*dsqrt(2d0)*AMW*(chtbrunr-chtbrunl)

c the result Eq.(19) in the paper

      gluinoex = (vs*(amq*(vq(ni)*vqp(nj)+aq(ni)*aqp(nj))+
     .     mgluino*(vq(ni)*vqp(nj)-aq(ni)*aqp(nj)))
     .     -as*(amq*(vq(ni)*aqp(nj)+aq(ni)*vqp(nj))+
     .     mgluino*(vq(ni)*aqp(nj)-aq(ni)*vqp(nj))))*
     .     SD_B02(gmst(ni)**2,mgluino,amq,amusq) +
     .     (vs*(amqp*(vq(ni)*vqp(nj)+aq(ni)*aqp(nj))+
     .     mgluino*(vq(ni)*vqp(nj)-aq(ni)*aqp(nj))) +
     .     as*(amqp*(vq(ni)*aqp(nj)+aq(ni)*vqp(nj)) -
     .     mgluino*(vq(ni)*aqp(nj)-aq(ni)*vqp(nj))))*
     .     SD_B02(gmsb(nj)**2,mgluino,amqp,amusq) +
     .     (vs*(amq*(vq(ni)*vqp(nj)+aq(ni)*aqp(nj))+
     .     amqp*(vq(ni)*vqp(nj)+aq(ni)*aqp(nj))) -
     .     as*(amq*(vq(ni)*aqp(nj)+aq(ni)*vqp(nj))-
     .     amqp*(vq(ni)*aqp(nj)+aq(ni)*vqp(nj))))*
     .     SD_B02(amchsq,amq,amqp,amusq) +
     .     (as*mgluino*(vq(ni)*aqp(nj)-aq(ni)*vqp(nj))*
     .      (amch**2-(amq-amqp)**2) -
     .      vs*mgluino*(vq(ni)*vqp(nj)-aq(ni)*aqp(nj))*
     .      (amch**2-(amq+amqp)**2) +
     .      as*(vq(ni)*aqp(nj)+aq(ni)*vqp(nj))*(gmsb(nj)**2*amq
     .      -gmst(ni)**2*amqp-(mgluino**2-amq*amqp)*(amq-amqp)) -
     .      vs*(vq(ni)*vqp(nj)+aq(ni)*aqp(nj))*(gmsb(nj)**2*amq
     .      +gmst(ni)**2*amqp-(mgluino**2+amq*amqp)*(amq+amqp)) )*
     . dreal(SD_C03(gmst(ni)**2,amchsq,gmsb(nj)**2,mgluino,amq,amqp))


      NS_stopsbot1719 = squarkstopsbot + gluinoex
      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_sbotstop1719(amuv,ni,nj)
*
      IMPLICIT NONE
      INTEGER ni,nj
      DOUBLE COMPLEX SD_C03
*
      DOUBLE PRECISION m(2,2,2),gctbr(2,2),gmsb(2),vq(2),aq(2),vqp(2),
     .          aqp(2),gmst(2)
      DOUBLE PRECISION amusq,coupphi,amuv,
     .gluinoex,amq,amqp,amchsq,vs,as
      DOUBLE PRECISION squarksbotstop
      DOUBLE PRECISION ama,aml,amh,amch,amar,amna,amnh
      DOUBLE PRECISION sdthet,sdtheb,sdthel,sdthem,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION chtbrunr,chtbrunl
      DOUBLE PRECISION SD_B02
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
*
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      COMMON/NS_sfmixang/sdthet,sdtheb,sdthel,sdthem,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_hcsbotstop/gctbr
      COMMON/NS_higgschudb/chtbrunr,chtbrunl
**
      aml = MAX(1d0,SMASS(1))
      amh = SMASS(2)
      amnh = SMASS(3)
      ama = MAX(1d0,PMASS(1))
      amar = ama
      amna = PMASS(2)
      amch = CMASS
*
      gmsb(1) = asb1
      gmsb(2) = asb2
**
      gmst(1) = ast1
      gmst(2) = ast2
*
      CALL NS_smatrix(m)
      amusq = amuv**2
**
      coupphi = m(nj,1,1)*gctbr(1,1)*m(1,ni,2)*
     .     SD_B02(amch**2,ast1,asb1,amusq) +
     .     m(nj,1,1)*gctbr(1,2)*m(2,ni,2)*
     .     SD_B02(amch**2,ast1,asb2,amusq) +
     .     m(nj,2,1)*gctbr(2,1)*m(1,ni,2)*
     .     SD_B02(amch**2,ast2,asb1,amusq) +
     .     m(nj,2,1)*gctbr(2,2)*m(2,ni,2)*
     .     SD_B02(amch**2,ast2,asb2,amusq)

      coupphi = dsqrt(2d0)*AMW**2*coupphi

c the result Eq.(17) in the paper

      squarksbotstop = coupphi

c --------------------------------------------

      amq  = runmt
      amqp = runmb

      amchsq = amch**2

      vqp(1) = 1d0/2d0*(cb-sb)
      vqp(2) = -1d0/2d0*(cb+sb)
      aqp(1) = -vqp(2)
      aqp(2) = vqp(1)

      vq(1) = 1d0/2d0*(ct-st)
      vq(2) = -1d0/2d0*(ct+st)
      aq(1) = -vq(2)
      aq(2) = vq(1)

      vs = 2d0*dsqrt(2d0)*AMW*(chtbrunr+chtbrunl)
      as = 2d0*dsqrt(2d0)*AMW*(chtbrunr-chtbrunl)

c the result Eq.(19) in the paper

      gluinoex = (vs*(amq*(vq(nj)*vqp(ni)+aq(nj)*aqp(ni))+
     .     mgluino*(vq(nj)*vqp(ni)-aq(nj)*aqp(ni)))
     .     -as*(amq*(vq(nj)*aqp(ni)+aq(nj)*vqp(ni))+
     .     mgluino*(vq(nj)*aqp(ni)-aq(nj)*vqp(ni))))*
     .     SD_B02(gmst(nj)**2,mgluino,amq,amusq) +
     .     (vs*(amqp*(vq(nj)*vqp(ni)+aq(nj)*aqp(ni))+
     .     mgluino*(vq(nj)*vqp(ni)-aq(nj)*aqp(ni))) +
     .     as*(amqp*(vq(nj)*aqp(ni)+aq(nj)*vqp(ni)) -
     .     mgluino*(vq(nj)*aqp(ni)-aq(nj)*vqp(ni))))*
     .     SD_B02(gmsb(ni)**2,mgluino,amqp,amusq) +
     .     (vs*(amq*(vq(nj)*vqp(ni)+aq(nj)*aqp(ni))+
     .     amqp*(vq(nj)*vqp(ni)+aq(nj)*aqp(ni))) -
     .     as*(amq*(vq(nj)*aqp(ni)+aq(nj)*vqp(ni))-
     .     amqp*(vq(nj)*aqp(ni)+aq(nj)*vqp(ni))))*
     .     SD_B02(amchsq,amq,amqp,amusq) +
     .     (as*mgluino*(vq(nj)*aqp(ni)-aq(nj)*vqp(ni))*
     .      (amch**2-(amq-amqp)**2) -
     .      vs*mgluino*(vq(nj)*vqp(ni)-aq(nj)*aqp(ni))*
     .      (amch**2-(amq+amqp)**2) +
     .      as*(vq(nj)*aqp(ni)+aq(nj)*vqp(ni))*(gmsb(ni)**2*amq
     .      -gmst(nj)**2*amqp-(mgluino**2-amq*amqp)*(amq-amqp)) -
     .      vs*(vq(nj)*vqp(ni)+aq(nj)*aqp(ni))*(gmsb(ni)**2*amq
     .      +gmst(nj)**2*amqp-(mgluino**2+amq*amqp)*(amq+amqp)) )*
     . dreal(SD_C03(gmst(nj)**2,amchsq,gmsb(ni)**2,mgluino,amq,amqp))

      NS_sbotstop1719 = squarksbotstop + gluinoex

      return

      end

c -------------------------------------------------------------------- c

      subroutine NS_smatrix(m)

      DOUBLE PRECISION m(2,2,2)
      DOUBLE PRECISION sdthet,sdtheb,sdthel,sdthem,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn

      COMMON/NS_sfmixang/sdthet,sdtheb,sdthel,sdthem,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn

      m(1,1,1) = ct**2-st**2
      m(1,2,1) = -2d0*st*ct
      m(2,1,1) = m(1,2,1)
      m(2,2,1) = -m(1,1,1)

      m(1,1,2) = cb**2-sb**2
      m(1,2,2) = -2d0*sb*cb
      m(2,1,2) = m(1,2,2)
      m(2,2,2) = -m(1,1,2)

      end
c -------------------------------------------------------------------- c

c -------------------------------------------------------------------- c
c -- the counterterm for squark2 -> h/H/A squark1 --

      DOUBLE PRECISION FUNCTION NS_dcounterneut(amsq1,amsq2,amq,theq,
     . mgluino,amuv,amuvdiv,lamv,nq,nh)
*
      IMPLICIT NONE
      INTEGER nq,nh
      DOUBLE PRECISION lamv,amq,theq,amuv,amuvdiv,amqq,amsq1,amsq2
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION dtlttr(3,2,2),dtattr(2,2,2),
     .datlttr(3,2,2),datattr(2,2,2),
     .dthlttr(3,2,2),dthattr(2,2,2)
      DOUBLE PRECISION dblbbr(3,2,2),dbabbr(2,2,2),
     .dablbbr(3,2,2),dababbr(2,2,2),
     .dthlbbr(3,2,2),dthabbr(2,2,2)
      DOUBLE PRECISION gqqr(2,2),dmqqr(2,2),daqr(2,2),dthr(2,2)
      DOUBLE PRECISION mgluino
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION NS_deltaz,NS_deltamqdiv,
     .NS_deltaAq,NS_deltathdiv
      DOUBLE PRECISION Hstopstopr(3,2,2),Astopstopr(3,2,2)
      DOUBLE PRECISION Hsbotsbotr(3,2,2),Asbotsbotr(2,2,2)
      DOUBLE PRECISION runmt,runmb,rmtauc
*
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_HIGGSSTST/Hstopstopr,Astopstopr
      COMMON/NS_HIGGSBTBT/Hsbotsbotr,Asbotsbotr
      COMMON/NS_higgsst1st2deriv/
     .dtlttr,dtattr,datlttr,datattr,dthlttr,dthattr
      COMMON/NS_higgssb1sb2deriv/dblbbr,dbabbr,
     .dablbbr,dababbr,dthlbbr,dthabbr

c --- the running mass ---

      if(nq.eq.1) then
         amqq = runmt
      elseif(nq.eq.2) then
         amqq = runmb
      endif
C
c --- some definitions ---
C ----- FOR OTHER HIGGS--CORRECTIONS ARE ADDEd--------------

      if(nq.eq.1.and.nh.eq.1) then
         gqqr(2,1)  = Hstopstopr(1,1,2)
         dmqqr(2,1) = dtlttr(1,2,1)
         daqr(2,1)  = datlttr(1,2,1)
         dthr(2,1)  = dthlttr(1,2,1)
      elseif(nq.eq.1.and.nh.eq.2) then
         gqqr(2,1)  = Hstopstopr(2,1,2)
         dmqqr(2,1) = dtlttr(2,2,1)
         daqr(2,1)  = datlttr(2,2,1)
         dthr(2,1)  = dthlttr(2,2,1)
C --NMSSM
      elseif(nq.eq.1.and.nh.eq.3) then
         gqqr(2,1)  = Hstopstopr(3,1,2)
         dmqqr(2,1) = dtlttr(3,2,1)
         daqr(2,1)  = datlttr(3,2,1)
         dthr(2,1)  = dthlttr(3,2,1)
      elseif(nq.eq.1.and.nh.eq.4) then
         gqqr(2,1)  = -Astopstopr(1,1,2)
         dmqqr(2,1) = -dtattr(1,2,1)
         daqr(2,1)  = -datattr(1,2,1)
         dthr(2,1)  = -dthattr(1,2,1)
      elseif(nq.eq.1.and.nh.eq.5) then
         gqqr(2,1)  = -Astopstopr(2,1,2)
         dmqqr(2,1) = -dtattr(2,2,1)
         daqr(2,1)  = -datattr(2,2,1)
         dthr(2,1)  = -dthattr(2,2,1)
C --NMSSM
C! NEGATIVE SIGN IN FRON OF dblbbr(1,2,1) AND dablbbr(1,2,1) HAVE BEEN ADDED
C! FOR CONSISTENCY WITH SDECAY FORMULA. SIMILARLY THERE SHOULD NOT BE ANY
C! NEGATIVE SIGN BEFORE dblbbr(2,2,1) AND dablbbr(2,2,1)

      elseif(nq.eq.2.and.nh.eq.1) then
         gqqr(2,1)  = Hsbotsbotr(1,1,2)
         dmqqr(2,1) = -dblbbr(1,2,1)
         daqr(2,1)  = -dablbbr(1,2,1)
         dthr(2,1)  = dthlbbr(1,2,1)
      elseif(nq.eq.2.and.nh.eq.2) then
         gqqr(2,1)  = Hsbotsbotr(2,1,2)
         dmqqr(2,1) = dblbbr(2,2,1)
         daqr(2,1)  = dablbbr(2,2,1)
         dthr(2,1)  = dthlbbr(2,2,1)
C --NMSSM
      elseif(nq.eq.2.and.nh.eq.3) then
         gqqr(2,1)  = Hsbotsbotr(3,1,2)
         dmqqr(2,1) = dblbbr(3,2,1)
         daqr(2,1)  = dablbbr(3,2,1)
         dthr(2,1)  = dthlbbr(3,2,1)
      elseif(nq.eq.2.and.nh.eq.4) then
         gqqr(2,1)  = -Asbotsbotr(1,1,2)
         dmqqr(2,1) = -dbabbr(1,2,1)
         daqr(2,1)  = -dababbr(1,2,1)
         dthr(2,1)  = -dthabbr(1,2,1)
      elseif(nq.eq.2.and.nh.eq.5) then
         gqqr(2,1)  = -Asbotsbotr(2,1,2)
         dmqqr(2,1) = -dbabbr(2,2,1)
         daqr(2,1)  = -dababbr(2,2,1)
         dthr(2,1)  = -dthabbr(2,2,1)
C --NMSSM
      endif

      NS_dcounterneut = gqqr(2,1)/2d0*(
     .     NS_deltaz(amsq2,mgluino,amq,theq,amuv,lamv,2)
     .     +NS_deltaz(amsq1,mgluino,amq,theq,amuv,lamv,1) )
     .     +dmqqr(2,1)*amqq*
     .     NS_deltamqdiv(amsq1,amsq2,mgluino,amq,theq,amuvdiv,lamv) +
     .     daqr(2,1)*
     .     NS_deltaAq(amsq2,amsq1,amsq2,mgluino,amq,theq,amuvdiv,lamv,
     .                nq) +
     .   dthr(2,1)*NS_deltathdiv(amsq1,amsq2,mgluino,amq,theq,amuvdiv)

      NS_dcounterneut = -dsqrt(2d0)*amz**2*NS_dcounterneut

      return

      end
c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_deltaz(amsq,mgluino,amq,theq,amuv,
     .                                 lamv,ni)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      integer ni
      DOUBLE PRECISION SD_B02,SD_BP02
      DOUBLE PRECISION DDCOS,DDSIN

      amusq = amuv**2

      NS_deltaz = 2d0*(
     .     (amq**2+mgluino**2-amsq**2)*
     .     SD_BP02(amsq**2,amq,mgluino,amusq) -
     .     SD_B02(amsq**2,amq,mgluino,amusq)
     .     + SD_B02(amsq**2,lamv,amsq,amusq) +
     .     2d0*(-1d0)**ni*DDSIN(2d0*theq)*mgluino*amq*
     .     SD_BP02(amsq**2,amq,mgluino,amusq)
     .     +2d0*amsq**2*
     .     SD_BP02(amsq**2,lamv,amsq,amusq) )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_deltaAq(amdec,amsq1,amsq2,mgluino,
     .                                  amq,theq,amuv,lamv,nq)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      integer nq
      DOUBLE PRECISION NS_deltamqdiv,NS_deltathdiv,NS_deltamsqdiv
      DOUBLE PRECISION DDCOS,DDSIN

c --- the running masses ---

      amqq = amq

c --------------------------

      NS_deltaAq = (amsq1**2-amsq2**2)/(2d0*amqq)*(
     .     2d0*DDCOS(2d0*theq)*
     .     NS_deltathdiv(amsq1,amsq2,mgluino,amq,theq,amuv)
     .     -DDSIN(2d0*theq)*
     .     NS_deltamqdiv(amsq1,amsq2,mgluino,amq,theq,amuv,lamv) )
     .     +DDSIN(2d0*theq)*
     .     NS_deltamsqdiv(amsq1,mgluino,amq,theq,amuv,lamv,1,nq)/amqq
     .     -DDSIN(2d0*theq)*
     .     NS_deltamsqdiv(amsq2,mgluino,amq,theq,amuv,lamv,2,nq)/amqq

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_deltathdiv(amsq1,amsq2,mgluino,amq,
     .                                        theq,amuv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION SD_B02_DIV
      DOUBLE PRECISION DDCOS,DDSIN

      amusq = amuv**2

      NS_deltathdiv = 1d0/(amsq1**2-amsq2**2)*(4d0*mgluino*amq*
     .     DDCOS(2d0*theq)/2d0*(
     .     SD_B02_DIV(amsq1**2,amq,mgluino,amusq)+
     .     SD_B02_DIV(amsq2**2,amq,mgluino,amusq) ) +
     .     DDCOS(2d0*theq)*DDSIN(2d0*theq)*(amsq2**2*dlog(amusq)
     .     -amsq1**2*log(amusq) ) )

      return

      end


c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_deltamqdiv(amsq1,amsq2,mgluino,amq,
     .                                        theq,amuv,lamv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_B1_DIV,SD_B02_DIV
      DOUBLE PRECISION DDCOS,DDSIN

      amusq = amuv**2

      NS_deltamqdiv = -(2d0*NS_B1_DIV(amq**2,amq,lamv,amusq) +
     .     4d0*SD_B02_DIV(amq**2,amq,lamv,amusq) +
     .     NS_B1_DIV(amq**2,mgluino,amsq1,amusq) +
     .     NS_B1_DIV(amq**2,mgluino,amsq2,amusq) ) +
     .     DDSIN(2d0*theq)*mgluino/amq*(
     .     SD_B02_DIV(amq**2,mgluino,amsq1,amusq) -
     .     SD_B02_DIV(amq**2,mgluino,amsq2,amusq) )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_deltamsqdiv(amsq,mgluino,amq,theq,
     .                                         amuv,lamv,i,nq)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)
      INTEGER NQ
      DIMENSION gmsq(2)
      DOUBLE PRECISION NS_A01_DIV,SD_B02_DIV
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1

      amusq = amuv**2
      ip = 3 - i

      if(nq.eq.1) then
         gmsq(1) = ast1
         gmsq(2) = ast2
      elseif(nq.eq.2) then
         gmsq(1) = asb1
         gmsq(2) = asb2
      endif

      NS_deltamsqdiv = -(amq**2+mgluino**2-amsq**2)*
     .     SD_B02_DIV(amsq**2,amq,mgluino,amusq) - 2d0*amsq**2*
     .     SD_B02_DIV(amsq**2,lamv,amsq,amusq) -
     .     NS_A01_DIV(mgluino**2,amusq) - NS_A01_DIV(amq**2,amusq) +
     .     1d0/2d0*((1d0+DDCOS(2d0*theq)**2)*
     .     NS_A01_DIV(amsq**2,amusq)
     .     + DDSIN(2d0*theq)**2*NS_A01_DIV(gmsq(ip)**2,amusq) ) -
     .     2d0*(-1d0)**i*DDSIN(2d0*theq)*mgluino*amq*
     .     SD_B02_DIV(amsq**2,amq,mgluino,amusq)

      return

      end

c -------------------------------------------------------------------- c
c ----------------------- The real corrections ----------------------- c

      DOUBLE PRECISION FUNCTION NS_realcorr(mphi,msq,msqp,lamv,nh,nq,
     .                                 ni,nj,scala)

      IMPLICIT NONE
      integer nh,nq,ni,nj
      DOUBLE PRECISION gctbr(2,2),lamv
      DOUBLE PRECISION msq,mphi,msqp,
     .GPHI,KAP,B0,B1,B2,lb0,lb1,lb2,scala
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION Hstopstopr(3,2,2),Astopstopr(3,2,2)
      DOUBLE PRECISION Hsbotsbotr(3,2,2),Asbotsbotr(2,2,2)
      DOUBLE COMPLEX NS_ccspen,NS_kappa
*
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_HIGGSSTST/Hstopstopr,Astopstopr
      COMMON/NS_HIGGSBTBT/Hsbotsbotr,Asbotsbotr
      COMMON/NS_hcsbotstop/gctbr

c --- some definitions ---

      if(nh.eq.6.and.nq.eq.0) then
         gphi = -dsqrt(2d0)*amw**2*gctbr(ni,nj)
      endif

      if(nh.eq.1.and.nq.eq.1) then
         gphi = -dsqrt(2d0)*amz**2*Hstopstopr(1,1,2)
      elseif(nh.eq.1.and.nq.eq.2) then
         gphi = -dsqrt(2d0)*amz**2*Hsbotsbotr(1,1,2)
      elseif(nh.eq.2.and.nq.eq.1) then
         gphi = -dsqrt(2d0)*amz**2*Hstopstopr(2,1,2)
      elseif(nh.eq.2.and.nq.eq.2) then
         gphi = -dsqrt(2d0)*amz**2*Hsbotsbotr(2,1,2)
c --NMSSM
      elseif(nh.eq.4.and.nq.eq.1) then
         gphi = dsqrt(2d0)*amz**2*Astopstopr(1,1,2)
      elseif(nh.eq.4.and.nq.eq.2) then
         gphi = dsqrt(2d0)*amz**2*Asbotsbotr(1,1,2)
      elseif(nh.eq.3.and.nq.eq.1) then
         gphi = -dsqrt(2d0)*amz**2*Hstopstopr(3,1,2)
      elseif(nh.eq.3.and.nq.eq.2) then
         gphi = -dsqrt(2d0)*amz**2*Hsbotsbotr(3,1,2)
      elseif(nh.eq.5.and.nq.eq.1) then
         gphi = dsqrt(2d0)*amz**2*Astopstopr(2,1,2)
      elseif(nh.eq.5.and.nq.eq.2) then
         gphi = dsqrt(2d0)*amz**2*Asbotsbotr(2,1,2)
C --NMSSM
      endif

      kap = dreal(NS_kappa(msq**2,mphi**2,msqp**2,0d0))

      b0 = (msq**2-mphi**2-msqp**2+kap)/(2d0*mphi*msqp)

      b1 = (msq**2-mphi**2+msqp**2-kap)/(2d0*msq*msqp)

      b2 = (mphi**2+msq**2-msqp**2-kap)/(2d0*mphi*msq)

c -----------------------

      lb0 = dreal(cdlog(dcmplx(b0)))
      lb1 = dreal(cdlog(dcmplx(b1)))
      lb2 = dreal(cdlog(dcmplx(b2)))

      NS_realcorr = 2d0*gphi/kap*(
     .     (mphi**2-msq**2-msqp**2)*(
     .     -2d0*dlog((lamv*mphi*msq*msqp)/kap**2)*lb1
     .     +2d0*lb1**2-lb0**2-lb2**2
     .     +2d0*dreal(NS_ccspen(dcmplx(1d0-b1**2)))
     .     -dreal(NS_ccspen(dcmplx(1d0-b0**2)))
     .     -dreal(NS_ccspen(dcmplx(1d0-b2**2))) )
     .     +2d0*kap*dlog((lamv*mphi*msq*msqp)/kap**2)
     .     +4d0*kap+(2d0*mphi**2+msq**2+msqp**2)*lb1
     .     +(mphi**2+2d0*msqp**2)*lb2+(mphi**2+2d0*msq**2)*lb0 )

      return

      end
c -------------------------------------------------------------------- c
c -------------------------------------------------------------------- c
c ---------- A.Djouadi, W.Hollik, C.Juenger, hep-ph/9609419 ---------- c
c -------------------------------------------------------------------- c

c --- QCD corrections to the light squark decays --- c

      DOUBLE PRECISION FUNCTION NS_ftotqcd(kap,gam)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)

      DOUBLE COMPLEX NS_ccspen,funci

      pi = 4d0*datan(1d0)

      if(kap*gam.lt.1d0) then
         funci = NS_ccspen(dcmplx((gam-1d0)/(gam*kap-1d0))) -
     .        NS_ccspen(dcmplx(kap*(gam-1d0)/(gam*kap-1d0))) -
     .        NS_ccspen(dcmplx((gam+kap-2d0)/(gam*kap-1d0))) +
     .        NS_ccspen(dcmplx(kap*(gam+kap-2d0)/(gam*kap-1d0)))
      elseif(kap*gam.ge.1d0) then
         funci = -NS_ccspen(dcmplx((gam*kap-1d0)/(gam-1d0))) +
     .        NS_ccspen(dcmplx((gam*kap-1d0)/(gam+kap-2d0))) +
     .        NS_ccspen(dcmplx((gam*kap-1d0)/(kap*(gam-1d0)))) -
     .        NS_ccspen(dcmplx((gam*kap-1d0)/(kap*(gam+kap-2d0)))) -
     .        dlog(kap)*cdlog(dcmplx((gam+kap-2d0)/(gam-1d0)))
      endif

      NS_ftotqcd = -1d0/8d0*( (4d0*gam**2-27d0*gam+25d0)/
     .(gam-1d0)
     .     + (3d0*kap-5d0)/(kap-1d0) ) - pi**2/3d0
     .     - 2d0*dreal(NS_ccspen(dcmplx(kap)))
     .     - 1d0/2d0*(gam**2-1d0)*
     .     dreal(cdlog(dcmplx((gam-1d0)/gam)))
     .     + (3d0*gam**2-4d0*gam+2d0)/
     .     (4d0*(1d0-gam)**2)*dlog(gam) -3d0/2d0*dlog(1d0-kap) +
     .     3d0/4d0*(kap**2-4d0*kap)/(kap-1d0)**2*dlog(kap)
     .     -dlog(kap)*dlog(1d0-kap) + dsqrt(kap*gam)*
     .     (1d0/kap*dlog(1d0-kap)+1d0/(1d0-kap)*(gam*dlog(gam)
     .     -(gam-1d0)*dreal(cdlog(dcmplx(gam-1d0))) )
     .     + (kap+gam-2d0)/(1d0-kap)**2*dreal(funci) )

      return

      end

c -------------------------------------------------------------------- c
C SEE REF PRD55,6975 QCD CORRECTIONS TO SCALAR QURAK DECAYS BY A. DJOUADI
C ET.AL
c --- Heavy squark decays                       --- c
c --- Virtual corrections for the decays        --- c
c --- squark_i -> chargino_j/neutralino_j quark --- c

      DOUBLE PRECISION FUNCTION NS_gltneut(ni,nj,amusc,amuscdiv,lamsc)
*
      IMPLICIT NONE
      INTEGER ni,nj,k,idec,I,J
      DOUBLE PRECISION gmst(2),atopr(2,5),btopr(2,5),vt(2),
     . at(2),del(2,2)
      DOUBLE PRECISION fnt1(2,5),fnt2(2,5),fct1(2,2),fct2(2,2),
     .     fnt1ik(2,5,2),fnt2ik(2,5,2),fnt3ik(2,5,2),fnt4ik(2,5,2),
     .     fnt5ik(2,5,2),fnt6ik(2,5,2),fnt7ik(2,5,2),
     .     fct1ik(2,2,2),fct2ik(2,2,2),fct3ik(2,2,2),fct4ik(2,2,2),
     .     fct5ik(2,2,2),fct6ik(2,2,2),fct7ik(2,2,2)
      DOUBLE PRECISION amuv,lamv,Z(5,5)
      DOUBLE PRECISION thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION amusc,lamsc,amuvdiv,amuscdiv,BET,amqt
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION SD_B02,NS_delztr,NS_delzst,NS_delthdiv,
     ,NS_A01,NS_delmtdiv
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_sfmixang/thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1

      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_neutstoptop/atopr,btopr

        do i =1,5
           do j=1,5
              Z(i,j) = ZZ(i,j)
           enddo
        enddo
      bet = datan(TANBETA_Z)
      gmst(1) = ast1
      gmst(2) = ast2

      idec = ni

      amuv    = amusc
      lamv    = lamsc
      amuvdiv = amuscdiv

      vt(1) = 1d0/2d0*(DDCOS(thet)-DDSIN(thet))
      vt(2) = -1d0/2d0*(DDCOS(thet)+DDSIN(thet))
      at(1) = -vt(2)
      at(2) = vt(1)

      del(1,1) = 1d0
      del(1,2) = 0d0
      del(2,1) = 0d0
      del(2,2) = 1d0

      CALL NS_ftFUNCTIONs(fnt1,fnt2,fnt1ik,fnt2ik,fnt3ik,fnt4ik,fnt5ik,
     .     fnt6ik,fnt7ik,fct1,fct2,fct1ik,fct2ik,fct3ik,fct4ik,fct5ik,
     .     fct6ik,fct7ik)

c --- the running couplings ---

      amqt = amt

      NS_gltneut = 0d0

      do k=1,2
         NS_gltneut = NS_gltneut - 2d0*(
     .        atopr(k,nj)*(
     .        (vt(k)*vt(ni)+at(k)*at(ni))*fnt4ik(ni,nj,k) -
     .        (at(k)*vt(ni)+vt(k)*at(ni))*fnt5ik(ni,nj,k) +
     .        (vt(k)*vt(ni)-at(k)*at(ni))*fnt6ik(ni,nj,k) -
     .        (at(k)*vt(ni)-vt(k)*at(ni))*fnt7ik(ni,nj,k) ) +
     .        btopr(k,nj)*(
     .        (vt(k)*vt(ni)+at(k)*at(ni))*fnt1ik(ni,nj,k) -
     .        (at(k)*vt(ni)+vt(k)*at(ni))*fnt1ik(ni,nj,k) +
     .        (vt(k)*vt(ni)-at(k)*at(ni))*fnt2ik(ni,nj,k) -
     .        (at(k)*vt(ni)-vt(k)*at(ni))*fnt3ik(ni,nj,k)) )
      enddo

      NS_gltneut = NS_gltneut + btopr(ni,nj)*fnt1(ni,nj) +
     .     atopr(ni,nj)*fnt2(ni,nj)

      NS_gltneut = NS_gltneut + (-1d0)**ni*(del(1,ni)*btopr(2,nj)+
     .     del(2,ni)*btopr(1,nj))/(ast1**2-ast2**2)*(
     .     4d0*amqt*mgluino*DDCOS(2d0*thet)*
     .     SD_B02(gmst(ni)**2,amqt,mgluino,amuv**2) +
     .     DDCOS(2d0*thet)*DDSIN(2d0*thet)*
     .     (NS_A01(ast2**2,amuv**2)-NS_A01(ast1**2,amuv**2)))

      if(ni.eq.1) then
         NS_gltneut = NS_gltneut + 1d0/2d0*btopr(1,nj)*(
     .        NS_delztr(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast1,thet,amuv,lamv,1) ) -
     .        1d0/(dsqrt(2d0)*amw*DDSIN(bet))*z(nj,4)*DDCOS(thet)*
     .       NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv)*
     .        runmt +
     .        runmt/(dsqrt(2d0)*amw*DDSIN(bet))*z(nj,4)*DDSIN(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) -
     .        (-dsqrt(2d0))*sw*(2d0/3d0*zp(nj,1)-2d0/3d0*sw/cw*
     .        zp(nj,2))*DDCOS(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv)

      elseif(ni.eq.2) then
         NS_gltneut = NS_gltneut + 1d0/2d0*btopr(2,nj)*(
     .        NS_delztr(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast2,thet,amuv,lamv,2) ) -
     .        1d0/(dsqrt(2d0)*amw*DDSIN(bet))*z(nj,4)*(-DDSIN(thet))*
     .        NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv)*
     .        runmt +
     .        runmt/(dsqrt(2d0)*amw*DDSIN(bet))*z(nj,4)*DDCOS(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) -
     .        (-dsqrt(2d0))*sw*(2d0/3d0*zp(nj,1)-2d0/3d0*sw/cw*
     .        zp(nj,2))*(-DDSIN(thet))*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv)

      endif
      NS_gltneut = (-1d0)*NS_gltneut


      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_glbneut(ni,nj,amusc,amuscdiv,lamsc)

      IMPLICIT NONE
      integer ni,nj,k,idec,I,J

      DOUBLE PRECISION gmsb(2),abot(2,5),bbot(2,5),vb(2),ab(2),del(2,2)
      DOUBLE PRECISION fnb1(2,5),fnb2(2,5),fcb1(2,2),fcb2(2,2),
     .     fnb1ik(2,5,2),fnb2ik(2,5,2),fnb3ik(2,5,2),fnb4ik(2,5,2),
     .     fnb5ik(2,5,2),fnb6ik(2,5,2),fnb7ik(2,5,2),
     .     fcb1ik(2,2,2),fcb2ik(2,2,2),fcb3ik(2,2,2),fcb4ik(2,2,2),
     .     fcb5ik(2,2,2),fcb6ik(2,2,2),fcb7ik(2,2,2)
      DOUBLE PRECISION amuv,lamv,z(5,5)
      DOUBLE PRECISION thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION amusc,lamsc,amuvdiv,amuscdiv,BET,amqb
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION SD_B02,NS_delztr,NS_delzst,NS_delthdiv,
     .NS_A01,NS_delmtdiv
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/NS_sfmixang/thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_neutsbotbot/abot,bbot

      bet=datan(TANBETA_z)
       do i =1,5
           do j=1,5
             Z(i,j) = ZZ(i,j)
           enddo
      enddo

C     -------------------------------------
      gmsb(1) = asb1
      gmsb(2) = asb2

      idec = ni

      amuv    = amusc
      lamv    = lamsc
      amuvdiv = amuscdiv

      vb(1) = 1d0/2d0*(DDCOS(theb)-DDSIN(theb))
      vb(2) = -1d0/2d0*(DDCOS(theb)+DDSIN(theb))
      ab(1) = -vb(2)
      ab(2) = vb(1)

      del(1,1) = 1d0
      del(1,2) = 0d0
      del(2,1) = 0d0
      del(2,2) = 1d0

      CALL NS_fbFUNCTIONs(fnb1,fnb2,fnb1ik,fnb2ik,fnb3ik,fnb4ik,fnb5ik,
     .     fnb6ik,fnb7ik,fcb1,fcb2,fcb1ik,fcb2ik,fcb3ik,fcb4ik,fcb5ik,
     .     fcb6ik,fcb7ik)

c --- the running couplings ---

      amqb = amb

c -------------

      NS_glbneut = 0d0

      do k=1,2
         NS_glbneut = NS_glbneut - 2d0*(
     .        abot(k,nj)*(
     .        (vb(k)*vb(ni)+ab(k)*ab(ni))*fnb4ik(ni,nj,k) -
     .        (ab(k)*vb(ni)+vb(k)*ab(ni))*fnb5ik(ni,nj,k) +
     .        (vb(k)*vb(ni)-ab(k)*ab(ni))*fnb6ik(ni,nj,k) -
     .        (ab(k)*vb(ni)-vb(k)*ab(ni))*fnb7ik(ni,nj,k) ) +
     .        bbot(k,nj)*(
     .        (vb(k)*vb(ni)+ab(k)*ab(ni))*fnb1ik(ni,nj,k) -
     .        (ab(k)*vb(ni)+vb(k)*ab(ni))*fnb1ik(ni,nj,k) +
     .        (vb(k)*vb(ni)-ab(k)*ab(ni))*fnb2ik(ni,nj,k) -
     .        (ab(k)*vb(ni)-vb(k)*ab(ni))*fnb3ik(ni,nj,k)) )
      enddo

      NS_glbneut = NS_glbneut + bbot(ni,nj)*fnb1(ni,nj) +
     .     abot(ni,nj)*fnb2(ni,nj)

      NS_glbneut = NS_glbneut + (-1d0)**ni*(del(1,ni)*bbot(2,nj)+
     .     del(2,ni)*bbot(1,nj))/(asb1**2-asb2**2)*(
     .     4d0*amqb*mgluino*DDCOS(2d0*theb)*
     .     SD_B02(gmsb(ni)**2,amqb,mgluino,amuv**2) +
     .     DDCOS(2d0*theb)*DDSIN(2d0*theb)*
     .     (NS_A01(asb2**2,amuv**2)-NS_A01(asb1**2,amuv**2)))

      if(ni.eq.1) then
         NS_glbneut = NS_glbneut + 1d0/2d0*bbot(1,nj)*(
     .        NS_delztr(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb1,theb,amuv,lamv,1) ) -
     .        1d0/(dsqrt(2d0)*amw*DDCOS(bet))*z(nj,3)*DDCOS(theb)*
     .        NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv)*
     .        runmb +
     .        runmb/(dsqrt(2d0)*amw*DDCOS(bet))*z(nj,3)*DDSIN(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) -
     .        (-dsqrt(2d0))*sw*(-1d0/3d0*zp(nj,1)+1d0/3d0*sw/cw*
     .        zp(nj,2))*DDCOS(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv)
      elseif(ni.eq.2) then
         NS_glbneut = NS_glbneut + 1d0/2d0*bbot(2,nj)*(
     .        NS_delztr(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb2,theb,amuv,lamv,2) ) -
     .        1d0/(dsqrt(2d0)*amw*DDCOS(bet))*z(nj,3)*(-DDSIN(theb))*
     .        NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv)*
     .        runmb +
     .        runmb/(dsqrt(2d0)*amw*DDCOS(bet))*z(nj,3)*DDCOS(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) -
     .        (-dsqrt(2d0))*sw*(-1d0/3d0*zp(nj,1)+1d0/3d0*sw/cw*
     .        zp(nj,2))*(-DDSIN(theb))*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv)
      endif

      NS_glbneut = (-1d0)*NS_glbneut

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_grtneut(ni,nj,amusc,amuscdiv,lamsc)
*
      IMPLICIT NONE
      integer ni,nj,k,idec,I,J
      DOUBLE PRECISION gmst(2),atopr(2,5),btopr(2,5),vt(2),
     .at(2),del(2,2)
      DOUBLE PRECISION fnt1(2,5),fnt2(2,5),fct1(2,2),fct2(2,2),
     .     fnt1ik(2,5,2),fnt2ik(2,5,2),fnt3ik(2,5,2),fnt4ik(2,5,2),
     .     fnt5ik(2,5,2),fnt6ik(2,5,2),fnt7ik(2,5,2),
     .     fct1ik(2,2,2),fct2ik(2,2,2),fct3ik(2,2,2),fct4ik(2,2,2),
     .     fct5ik(2,2,2),fct6ik(2,2,2),fct7ik(2,2,2)

      DOUBLE PRECISION amuv,lamv,z(5,5)
      DOUBLE PRECISION thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION amusc,lamsc,amuvdiv,amuscdiv,BET,amqt
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION SD_B02,NS_delzst,NS_delthdiv,
     ,NS_A01,NS_delmtdiv,NS_delztl
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/NS_sfmixang/thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_neutstoptop/atopr,btopr

      bet=datan(TANBETA_Z)
      do i =1,5
           do j=1,5
              Z(i,j) = ZZ(i,j)
           enddo
      enddo
      gmst(1) = ast1
      gmst(2) = ast2

      idec = ni

      amuv    = amusc
      lamv    = lamsc
      amuvdiv = amuscdiv

      vt(1) = 1d0/2d0*(DDCOS(thet)-DDSIN(thet))
      vt(2) = -1d0/2d0*(DDCOS(thet)+DDSIN(thet))
      at(1) = -vt(2)
      at(2) = vt(1)

      del(1,1) = 1d0
      del(1,2) = 0d0
      del(2,1) = 0d0
      del(2,2) = 1d0

      CALL NS_ftFUNCTIONs(fnt1,fnt2,fnt1ik,fnt2ik,fnt3ik,fnt4ik,fnt5ik,
     .     fnt6ik,fnt7ik,fct1,fct2,fct1ik,fct2ik,fct3ik,fct4ik,fct5ik,
     .     fct6ik,fct7ik)

      amqt = amt

      NS_grtneut = 0d0

      do k=1,2
         NS_grtneut = NS_grtneut - 2d0*(
     .        btopr(k,nj)*(
     .        (vt(k)*vt(ni)+at(k)*at(ni))*fnt4ik(ni,nj,k) +
     .        (at(k)*vt(ni)+vt(k)*at(ni))*fnt5ik(ni,nj,k) +
     .        (vt(k)*vt(ni)-at(k)*at(ni))*fnt6ik(ni,nj,k) +
     .        (at(k)*vt(ni)-vt(k)*at(ni))*fnt7ik(ni,nj,k) ) +
     .        atopr(k,nj)*(
     .        (vt(k)*vt(ni)+at(k)*at(ni))*fnt1ik(ni,nj,k) +
     .        (at(k)*vt(ni)+vt(k)*at(ni))*fnt1ik(ni,nj,k) +
     .        (vt(k)*vt(ni)-at(k)*at(ni))*fnt2ik(ni,nj,k) +
     .        (at(k)*vt(ni)-vt(k)*at(ni))*fnt3ik(ni,nj,k)) )
      enddo

      NS_grtneut = NS_grtneut + atopr(ni,nj)*fnt1(ni,nj) +
     .         btopr(ni,nj)*fnt2(ni,nj)

      NS_grtneut = NS_grtneut + (-1d0)**ni*(del(1,ni)*atopr(2,nj)+
     .     del(2,ni)*atopr(1,nj))/(ast1**2-ast2**2)*(
     .     4d0*amqt*mgluino*DDCOS(2d0*thet)*
     .     SD_B02(gmst(ni)**2,amqt,mgluino,amuv**2) +
     .     DDCOS(2d0*thet)*DDSIN(2d0*thet)*
     .     (NS_A01(ast2**2,amuv**2)-NS_A01(ast1**2,amuv**2)))

        if(ni.eq.1) then

         NS_grtneut = NS_grtneut + 1d0/2d0*atopr(1,nj)*(
     .        NS_delztl(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast1,thet,amuv,lamv,1) ) -
     .        1d0/(dsqrt(2d0)*amw*DDSIN(bet))*z(nj,4)*DDSIN(thet)*
     .        NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv)*
     .        runmt -
     .        runmt/(dsqrt(2d0)*amw*DDSIN(bet))*z(nj,4)*DDCOS(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) +
     .        dsqrt(2d0)*sw*(2d0/3d0*zp(nj,1)+(1d0/2d0-
     .         2d0/3d0*sw**2)*1d0/sw/cw*zp(nj,2))*DDSIN(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv)

      elseif(ni.eq.2) then
         NS_grtneut = NS_grtneut + 1d0/2d0*atopr(2,nj)*(
     .        NS_delztl(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast2,thet,amuv,lamv,2) ) -
     .        1d0/(dsqrt(2d0)*amw*DDSIN(bet))*z(nj,4)*DDCOS(thet)*
     .        NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv)*
     .        runmt -
     .        runmt/(dsqrt(2d0)*amw*DDSIN(bet))*z(nj,4)*(-DDSIN(thet))*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) +
     .        dsqrt(2d0)*sw*(2d0/3d0*zp(nj,1)+(1d0/2d0-
     .         2d0/3d0*sw**2)*1d0/sw/cw*zp(nj,2))*DDCOS(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv)

      endif
      NS_grtneut = (-1d0)*NS_grtneut

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_grbneut(ni,nj,amusc,amuscdiv,lamsc)
*
      IMPLICIT NONE
      INTEGER ni,nj,k,idec,I,J
      DOUBLE PRECISION gmsb(2),abot(2,5),bbot(2,5),vb(2),ab(2),del(2,2)
      DOUBLE PRECISION fnb1(2,5),fnb2(2,5),fcb1(2,2),fcb2(2,2),
     .     fnb1ik(2,5,2),fnb2ik(2,5,2),fnb3ik(2,5,2),fnb4ik(2,5,2),
     .     fnb5ik(2,5,2),fnb6ik(2,5,2),fnb7ik(2,5,2),
     .     fcb1ik(2,2,2),fcb2ik(2,2,2),fcb3ik(2,2,2),fcb4ik(2,2,2),
     .     fcb5ik(2,2,2),fcb6ik(2,2,2),fcb7ik(2,2,2)
      DOUBLE PRECISION amuv,lamv,z(5,5)
      DOUBLE PRECISION thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION amusc,lamsc,amuvdiv,amuscdiv,BET,amqb
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION SD_B02,NS_delzst,NS_delthdiv,
     ,NS_A01,NS_delmtdiv,NS_delztl
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/NS_sfmixang/thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_neutsbotbot/abot,bbot

      bet=datan(TANBETA_Z)
      do i =1,5
           do j=1,5
              Z(i,j) = ZZ(i,j)
           enddo
      enddo

      gmsb(1) = asb1
      gmsb(2) = asb2

      idec = ni

      amuv = amusc
      lamv = lamsc
      amuvdiv = amuscdiv

      vb(1) = 1d0/2d0*(DDCOS(theb)-DDSIN(theb))
      vb(2) = -1d0/2d0*(DDCOS(theb)+DDSIN(theb))
      ab(1) = -vb(2)
      ab(2) = vb(1)

      del(1,1) = 1d0
      del(1,2) = 0d0
      del(2,1) = 0d0
      del(2,2) = 1d0

      CALL NS_fbFUNCTIONs(fnb1,fnb2,fnb1ik,fnb2ik,fnb3ik,fnb4ik,fnb5ik,
     .     fnb6ik,fnb7ik,fcb1,fcb2,fcb1ik,fcb2ik,fcb3ik,fcb4ik,fcb5ik,
     .     fcb6ik,fcb7ik)

      amqb = amb
      NS_grbneut= 0d0

      do k=1,2
        NS_grbneut=NS_grbneut- 2d0*(
     .        bbot(k,nj)*(
     .        (vb(k)*vb(ni)+ab(k)*ab(ni))*fnb4ik(ni,nj,k) +
     .        (ab(k)*vb(ni)+vb(k)*ab(ni))*fnb5ik(ni,nj,k) +
     .        (vb(k)*vb(ni)-ab(k)*ab(ni))*fnb6ik(ni,nj,k) +
     .        (ab(k)*vb(ni)-vb(k)*ab(ni))*fnb7ik(ni,nj,k) ) +
     .        abot(k,nj)*(
     .        (vb(k)*vb(ni)+ab(k)*ab(ni))*fnb1ik(ni,nj,k) +
     .        (ab(k)*vb(ni)+vb(k)*ab(ni))*fnb1ik(ni,nj,k) +
     .        (vb(k)*vb(ni)-ab(k)*ab(ni))*fnb2ik(ni,nj,k) +
     .        (ab(k)*vb(ni)-vb(k)*ab(ni))*fnb3ik(ni,nj,k)) )
      enddo

      NS_grbneut=NS_grbneut+ abot(ni,nj)*fnb1(ni,nj) +
     .         bbot(ni,nj)*fnb2(ni,nj)

      NS_grbneut=NS_grbneut+ (-1d0)**ni*(del(1,ni)*abot(2,nj)+
     .     del(2,ni)*abot(1,nj))/(asb1**2-asb2**2)*(
     .     4d0*amqb*mgluino*DDCOS(2d0*theb)*
     .     SD_B02(gmsb(ni)**2,amqb,mgluino,amuv**2) +
     .     DDCOS(2d0*theb)*DDSIN(2d0*theb)*
     .     (NS_A01(asb2**2,amuv**2)-NS_A01(asb1**2,amuv**2)))

      if(ni.eq.1) then
        NS_grbneut=NS_grbneut+ 1d0/2d0*abot(1,nj)*(
     .        NS_delztl(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb1,theb,amuv,lamv,1) ) -
     .        1d0/(dsqrt(2d0)*amw*DDCOS(bet))*z(nj,3)*DDSIN(theb)*
     .        NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv)*
     .        runmb -
     .        runmb/(dsqrt(2d0)*amw*DDCOS(bet))*z(nj,3)*DDCOS(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) +
     .        dsqrt(2d0)*sw*(-1d0/3d0*zp(nj,1)+(-1d0/2d0+
     .         1d0/3d0*sw**2)*1d0/sw/cw*zp(nj,2))*DDSIN(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv)
      elseif(ni.eq.2) then
        NS_grbneut=NS_grbneut+ 1d0/2d0*abot(2,nj)*(
     .        NS_delztl(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb2,theb,amuv,lamv,2) ) -
     .        1d0/(dsqrt(2d0)*amw*DDCOS(bet))*z(nj,3)*DDCOS(theb)*
     .        NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv)*
     .        runmb -
     .        runmb/(dsqrt(2d0)*amw*DDCOS(bet))*z(nj,3)*(-DDSIN(theb))*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) +
     .        dsqrt(2d0)*sw*(-1d0/3d0*zp(nj,1)+(-1d0/2d0+
     .         1d0/3d0*sw**2)*1d0/sw/cw*zp(nj,2))*DDCOS(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv)
      endif

      NS_grbneut= (-1d0)*NS_grbneut

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gltchar(ni,nj,amusc,amuscdiv,lamsc)
*
      IMPLICIT NONE
      integer ni,nj,k,idec,I,J
*
      DOUBLE PRECISION gmst(2),vt(2),at(2),vb(2),ab(2),del(2,2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2),alstor(2,2),akstor(2,2)
      DOUBLE PRECISION fnt1(2,5),fnt2(2,5),fct1(2,2),fct2(2,2),
     .     fnt1ik(2,5,2),fnt2ik(2,5,2),fnt3ik(2,5,2),fnt4ik(2,5,2),
     .     fnt5ik(2,5,2),fnt6ik(2,5,2),fnt7ik(2,5,2),
     .     fct1ik(2,2,2),fct2ik(2,2,2),fct3ik(2,2,2),fct4ik(2,2,2),
     .     fct5ik(2,2,2),fct6ik(2,2,2),fct7ik(2,2,2)
      DOUBLE PRECISION amuv,lamv,z(5,5)
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION amusc,lamsc,amuvdiv,amuscdiv,BET,amqt,amqb
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION SD_B02,NS_delztr,NS_delzst,NS_delthdiv,
     ,NS_A01,NS_delmtdiv
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/NS_sfmixang/thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1

      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor

       bet=datan(TANBETA_Z)
       do i =1,5
           do j=1,5
              Z(i,j) = ZZ(i,j)
           enddo
       enddo
*
      amuv    = amusc
      lamv    = lamsc
      amuvdiv = amuscdiv

      idec = ni

      gmst(1) = ast1
      gmst(2) = ast2

      vt(1) = 1d0/2d0*(DDCOS(thet)-DDSIN(thet))
      vt(2) = -1d0/2d0*(DDCOS(thet)+DDSIN(thet))
      at(1) = -vt(2)
      at(2) = vt(1)

      vb(1) = 1d0/2d0*(DDCOS(theb)-DDSIN(theb))
      vb(2) = -1d0/2d0*(DDCOS(theb)+DDSIN(theb))
      ab(1) = -vb(2)
      ab(2) = vb(1)

      del(1,1) = 1d0
      del(1,2) = 0d0
      del(2,1) = 0d0
      del(2,2) = 1d0

      CALL NS_ftFUNCTIONs(fnt1,fnt2,fnt1ik,fnt2ik,fnt3ik,fnt4ik,fnt5ik,
     .     fnt6ik,fnt7ik,fct1,fct2,fct1ik,fct2ik,fct3ik,fct4ik,fct5ik,
     .     fct6ik,fct7ik)

      amqt = amt
      amqb = amb

      NS_gltchar = 0d0

      do k=1,2
         NS_gltchar = NS_gltchar -2d0*(
     .        alsbot(k,nj)*(
     .        (vb(k)*vt(ni)+ab(k)*at(ni))*fct4ik(ni,nj,k) -
     .        (ab(k)*vt(ni)+vb(k)*at(ni))*fct5ik(ni,nj,k) +
     .        (vb(k)*vt(ni)-ab(k)*at(ni))*fct6ik(ni,nj,k) -
     .        (ab(k)*vt(ni)-vb(k)*at(ni))*fct7ik(ni,nj,k) ) +
     .        aksbot(k,nj)*(
     .        (vb(k)*vt(ni)+ab(k)*at(ni))*fct1ik(ni,nj,k) -
     .        (ab(k)*vt(ni)+vb(k)*at(ni))*fct1ik(ni,nj,k) +
     .        (vb(k)*vt(ni)-ab(k)*at(ni))*fct2ik(ni,nj,k) -
     .        (ab(k)*vt(ni)-vb(k)*at(ni))*fct3ik(ni,nj,k)) )
      enddo


     .
      NS_gltchar = NS_gltchar +
     .         akstor(ni,nj)*fct1(ni,nj) + alstor(ni,nj)*fct2(ni,nj)

      NS_gltchar = NS_gltchar + (-1d0)**ni*(del(1,ni)*akstor(2,nj)+
     .     del(2,ni)*akstor(1,nj))/(ast1**2-ast2**2)*(
     .     4d0*amqt*mgluino*DDCOS(2d0*thet)*
     .     SD_B02(gmst(ni)**2,amqt,mgluino,amuv**2) +
     .     DDCOS(2d0*thet)*DDSIN(2d0*thet)*
     .     (NS_A01(ast2**2,amuv**2)-NS_A01(ast1**2,amuv**2)))

         if(ni.eq.1) then
         NS_gltchar = NS_gltchar + 1d0/2d0*akstor(1,nj)*(
     .        NS_delztr(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast1,thet,amuv,lamv,1)  +
     .      2d0*NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv))
     .        - runmb*uu(nj,2)/dsqrt(2d0)/amw/DDCOS(bet)*DDSIN(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv)

      elseif(ni.eq.2) then
         NS_gltchar = NS_gltchar + 1d0/2d0*akstor(2,nj)*(
     .        NS_delztr(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast2,thet,amuv,lamv,2)  +
     .      2d0*NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv))
     .        - runmb*uu(nj,2)/dsqrt(2d0)/amw/DDCOS(bet)*DDCOS(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv)
      endif

      NS_gltchar = (-1d0)*NS_gltchar

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_glbchar(ni,nj,amusc,amuscdiv,lamsc)
*
      IMPLICIT NONE
      integer ni,nj,k,idec,I,J
      DOUBLE PRECISION gmsb(2),vt(2),at(2),vb(2),ab(2),del(2,2)
      DOUBLE PRECISION alstor(2,2),akstor(2,2),aksbot(2,2),alsbot(2,2)
      DOUBLE PRECISION fnb1(2,5),fnb2(2,5),fcb1(2,2),fcb2(2,2),
     .     fnb1ik(2,5,2),fnb2ik(2,5,2),fnb3ik(2,5,2),fnb4ik(2,5,2),
     .     fnb5ik(2,5,2),fnb6ik(2,5,2),fnb7ik(2,5,2),
     .     fcb1ik(2,2,2),fcb2ik(2,2,2),fcb3ik(2,2,2),fcb4ik(2,2,2),
     .     fcb5ik(2,2,2),fcb6ik(2,2,2),fcb7ik(2,2,2)
      DOUBLE PRECISION amuv,lamv,z(5,5)
**
      DOUBLE PRECISION thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION amusc,lamsc,amuvdiv,amuscdiv,BET,amqt,amqb
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION SD_B02,NS_delztr,NS_delzst,NS_delthdiv,
     ,NS_A01,NS_delmtdiv
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/NS_sfmixang/thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor

      bet=datan(TANBETA_Z)
c
      do i =1,5
           do j=1,5
              Z(i,j) = ZZ(i,j)
           enddo
      enddo

      amuv = amusc
      lamv = lamsc
      amuvdiv = amuscdiv

      idec = ni

      gmsb(1) = asb1
      gmsb(2) = asb2

      vt(1) = 1d0/2d0*(DDCOS(thet)-DDSIN(thet))
      vt(2) = -1d0/2d0*(DDCOS(thet)+DDSIN(thet))
      at(1) = -vt(2)
      at(2) = vt(1)

      vb(1) = 1d0/2d0*(DDCOS(theb)-DDSIN(theb))
      vb(2) = -1d0/2d0*(DDCOS(theb)+DDSIN(theb))
      ab(1) = -vb(2)
      ab(2) = vb(1)

      del(1,1) = 1d0
      del(1,2) = 0d0
      del(2,1) = 0d0
      del(2,2) = 1d0

      CALL NS_fbFUNCTIONs(fnb1,fnb2,fnb1ik,fnb2ik,fnb3ik,fnb4ik,fnb5ik,
     .     fnb6ik,fnb7ik,fcb1,fcb2,fcb1ik,fcb2ik,fcb3ik,fcb4ik,fcb5ik,
     .     fcb6ik,fcb7ik)

      amqt = amt
      amqb = amb

      NS_glbchar = 0d0

      do k=1,2
         NS_glbchar = NS_glbchar - 2d0*(
     .        alstor(k,nj)*(
     .        (vt(k)*vb(ni)+at(k)*ab(ni))*fcb4ik(ni,nj,k) -
     .        (at(k)*vb(ni)+vt(k)*ab(ni))*fcb5ik(ni,nj,k) +
     .        (vt(k)*vb(ni)-at(k)*ab(ni))*fcb6ik(ni,nj,k) -
     .        (at(k)*vb(ni)-vt(k)*ab(ni))*fcb7ik(ni,nj,k) ) +
     .        akstor(k,nj)*(
     .        (vt(k)*vb(ni)+at(k)*ab(ni))*fcb1ik(ni,nj,k) -
     .        (at(k)*vb(ni)+vt(k)*ab(ni))*fcb1ik(ni,nj,k) +
     .        (vt(k)*vb(ni)-at(k)*ab(ni))*fcb2ik(ni,nj,k) -
     .        (at(k)*vb(ni)-vt(k)*ab(ni))*fcb3ik(ni,nj,k)) )
      enddo

      NS_glbchar = NS_glbchar +
     .         aksbot(ni,nj)*fcb1(ni,nj) + alsbot(ni,nj)*fcb2(ni,nj)

      NS_glbchar = NS_glbchar + (-1d0)**ni*(del(1,ni)*aksbot(2,nj)+
     .     del(2,ni)*aksbot(1,nj))/(asb1**2-asb2**2)*(
     .     4d0*amqb*mgluino*DDCOS(2d0*theb)*
     .     SD_B02(gmsb(ni)**2,amqb,mgluino,amuv**2) +
     .     DDCOS(2d0*theb)*DDSIN(2d0*theb)*
     .     (NS_A01(asb2**2,amuv**2)-NS_A01(asb1**2,amuv**2)))

      if(ni.eq.1) then
         NS_glbchar = NS_glbchar + 1d0/2d0*aksbot(1,nj)*(
     .        NS_delztr(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb1,theb,amuv,lamv,1)  +
     .      2d0*NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv))
     .        - runmt*vv(nj,2)/dsqrt(2d0)/amw/DDSIN(bet)*DDSIN(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv)
      elseif(ni.eq.2) then
         NS_glbchar = NS_glbchar + 1d0/2d0*aksbot(2,nj)*(
     .        NS_delztr(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb2,theb,amuv,lamv,2)  +
     .      2d0*NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv))
     .        - runmt*vv(nj,2)/dsqrt(2d0)/amw/DDSIN(bet)*DDCOS(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv)
      endif

      NS_glbchar = (-1d0)*NS_glbchar

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_grtchar(ni,nj,amusc,amuscdiv,lamsc)
*
      IMPLICIT NONE
      integer ni,nj,k,idec,I,J
*
      DOUBLE PRECISION gmst(2),vt(2),at(2),vb(2),ab(2),del(2,2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2),alstor(2,2),akstor(2,2)
      DOUBLE PRECISION fnt1(2,5),fnt2(2,5),fct1(2,2),fct2(2,2),
     .     fnt1ik(2,5,2),fnt2ik(2,5,2),fnt3ik(2,5,2),fnt4ik(2,5,2),
     .     fnt5ik(2,5,2),fnt6ik(2,5,2),fnt7ik(2,5,2),
     .     fct1ik(2,2,2),fct2ik(2,2,2),fct3ik(2,2,2),fct4ik(2,2,2),
     .     fct5ik(2,2,2),fct6ik(2,2,2),fct7ik(2,2,2)
      DOUBLE PRECISION amuv,lamv,z(5,5)
      DOUBLE PRECISION thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION amusc,lamsc,amuvdiv,amuscdiv,BET,amqt,amqb
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION SD_B02,NS_delzst,NS_delthdiv,
     ,NS_A01,NS_delmtdiv,NS_delztl
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/NS_sfmixang/thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor

      bet=datan(TANBETA_Z)
        do i =1,5
           do j=1,5
              Z(i,j) = ZZ(i,j)
           enddo
        enddo
*
      amuv = amusc
      lamv = lamsc
      amuvdiv = amuscdiv

      idec = ni

      gmst(1) = ast1
      gmst(2) = ast2

      vt(1) = 1d0/2d0*(DDCOS(thet)-DDSIN(thet))
      vt(2) = -1d0/2d0*(DDCOS(thet)+DDSIN(thet))
      at(1) = -vt(2)
      at(2) = vt(1)

      vb(1) = 1d0/2d0*(DDCOS(theb)-DDSIN(theb))
      vb(2) = -1d0/2d0*(DDCOS(theb)+DDSIN(theb))
      ab(1) = -vb(2)
      ab(2) = vb(1)

      del(1,1) = 1d0
      del(1,2) = 0d0
      del(2,1) = 0d0
      del(2,2) = 1d0

      CALL NS_ftFUNCTIONs(fnt1,fnt2,fnt1ik,fnt2ik,fnt3ik,fnt4ik,fnt5ik,
     .     fnt6ik,fnt7ik,fct1,fct2,fct1ik,fct2ik,fct3ik,fct4ik,fct5ik,
     .     fct6ik,fct7ik)
      amqt = amt
      amqb = amb

      NS_grtchar = 0d0

      do k=1,2
         NS_grtchar = NS_grtchar -2d0*(
     .        aksbot(k,nj)*(
     .        (vb(k)*vt(ni)+ab(k)*at(ni))*fct4ik(ni,nj,k) +
     .        (ab(k)*vt(ni)+vb(k)*at(ni))*fct5ik(ni,nj,k) +
     .        (vb(k)*vt(ni)-ab(k)*at(ni))*fct6ik(ni,nj,k) +
     .        (ab(k)*vt(ni)-vb(k)*at(ni))*fct7ik(ni,nj,k) ) +
     .        alsbot(k,nj)*(
     .        (vb(k)*vt(ni)+ab(k)*at(ni))*fct1ik(ni,nj,k) +
     .        (ab(k)*vt(ni)+vb(k)*at(ni))*fct1ik(ni,nj,k) +
     .        (vb(k)*vt(ni)-ab(k)*at(ni))*fct2ik(ni,nj,k) +
     .        (ab(k)*vt(ni)-vb(k)*at(ni))*fct3ik(ni,nj,k)) )
      enddo

      NS_grtchar = NS_grtchar +
     .         alstor(ni,nj)*fct1(ni,nj) + akstor(ni,nj)*fct2(ni,nj)

      NS_grtchar = NS_grtchar + (-1d0)**ni*(del(1,ni)*alstor(2,nj)+
     .     del(2,ni)*alstor(1,nj))/(ast1**2-ast2**2)*(
     .     4d0*amqt*mgluino*DDCOS(2d0*thet)*
     .     SD_B02(gmst(ni)**2,amqt,mgluino,amuv**2) +
     .     DDCOS(2d0*thet)*DDSIN(2d0*thet)*
     .     (NS_A01(ast2**2,amuv**2)-NS_A01(ast1**2,amuv**2)))

      if(ni.eq.1) then
         NS_grtchar = NS_grtchar + 1d0/2d0*alstor(1,nj)*(
     .        NS_delztl(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast1,thet,amuv,lamv,1) ) +
     .        1d0/(dsqrt(2d0)*amw*DDSIN(bet))*vv(nj,2)*DDSIN(thet)*
     .        NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv)*
     .        runmt +
     .        vv(nj,1)*DDSIN(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) +
     .        runmt*vv(nj,2)/dsqrt(2d0)/amw/DDSIN(bet)*DDCOS(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv)
      elseif(ni.eq.2) then
         NS_grtchar = NS_grtchar + 1d0/2d0*alstor(2,nj)*(
     .        NS_delztl(amqb,mgluino,asb1,asb2,theb,amuv,lamv) +
     .        NS_delzst(amqt,mgluino,ast2,thet,amuv,lamv,2) ) +
     .        1d0/(dsqrt(2d0)*amw*DDSIN(bet))*vv(nj,2)*DDCOS(thet)*
     .        NS_delmtdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv,lamv)*
     .        runmt +
     .        vv(nj,1)*DDCOS(thet)*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv) +
     .        runmt*vv(nj,2)/dsqrt(2d0)/amw/DDSIN(bet)*(-DDSIN(thet))*
     .        NS_delthdiv(amqt,mgluino,ast1,ast2,thet,amuvdiv)
      endif

      NS_grtchar = (-1d0)*NS_grtchar
      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_grbchar(ni,nj,amusc,amuscdiv,lamsc)
*
      IMPLICIT NONE
      integer ni,nj,k,idec,I,J
*
      DOUBLE PRECISION gmsb(2),vt(2),at(2),vb(2),ab(2),del(2,2)
      DOUBLE PRECISION alsbot(2,2),aksbot(2,2),alstor(2,2),akstor(2,2)
      DOUBLE PRECISION fnb1(2,5),fnb2(2,5),fcb1(2,2),fcb2(2,2),
     .     fnb1ik(2,5,2),fnb2ik(2,5,2),fnb3ik(2,5,2),fnb4ik(2,5,2),
     .     fnb5ik(2,5,2),fnb6ik(2,5,2),fnb7ik(2,5,2),
     .     fcb1ik(2,2,2),fcb2ik(2,2,2),fcb3ik(2,2,2),fcb4ik(2,2,2),
     .     fcb5ik(2,2,2),fcb6ik(2,2,2),fcb7ik(2,2,2)
      DOUBLE PRECISION amuv,lamv,z(5,5)
      DOUBLE PRECISION uu(2,2),vv(2,2),zz(5,5),zp(5,5)
      DOUBLE PRECISION tanbeta_Z
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION lamsc,amuvdiv,amusc,amuscdiv,BET,amqt,amqb
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION SD_B02,NS_delzst,NS_delthdiv,
     ,NS_A01,NS_delmtdiv,NS_delztl
      DOUBLE PRECISION thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/NS_sfmixang/thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_mixmat/uu,vv,zz,zp
      COMMON/NS_tanb/tanbeta_Z
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor

      bet=datan(TANBETA_Z)
        do i =1,5
           do j=1,5
              Z(i,j) = ZZ(i,j)
           enddo
        enddo

      amuv    = amusc
      lamv    = lamsc
      amuvdiv = amuscdiv

      idec = ni

      gmsb(1) = asb1
      gmsb(2) = asb2

      vt(1) = 1d0/2d0*(DDCOS(thet)-DDSIN(thet))
      vt(2) = -1d0/2d0*(DDCOS(thet)+DDSIN(thet))
      at(1) = -vt(2)
      at(2) = vt(1)

      vb(1) = 1d0/2d0*(DDCOS(theb)-DDSIN(theb))
      vb(2) = -1d0/2d0*(DDCOS(theb)+DDSIN(theb))
      ab(1) = -vb(2)
      ab(2) = vb(1)

      del(1,1) = 1d0
      del(1,2) = 0d0
      del(2,1) = 0d0
      del(2,2) = 1d0

      CALL NS_fbFUNCTIONs(fnb1,fnb2,fnb1ik,fnb2ik,fnb3ik,fnb4ik,fnb5ik,
     .     fnb6ik,fnb7ik,fcb1,fcb2,fcb1ik,fcb2ik,fcb3ik,fcb4ik,fcb5ik,
     .     fcb6ik,fcb7ik)

      amqt = amt
      amqb = amb

      NS_grbchar = 0d0

      do k=1,2
         NS_grbchar = NS_grbchar - 2d0*(
     .        akstor(k,nj)*(
     .        (vt(k)*vb(ni)+at(k)*ab(ni))*fcb4ik(ni,nj,k) +
     .        (at(k)*vb(ni)+vt(k)*ab(ni))*fcb5ik(ni,nj,k) +
     .        (vt(k)*vb(ni)-at(k)*ab(ni))*fcb6ik(ni,nj,k) +
     .        (at(k)*vb(ni)-vt(k)*ab(ni))*fcb7ik(ni,nj,k) ) +
     .        alstor(k,nj)*(
     .        (vt(k)*vb(ni)+at(k)*ab(ni))*fcb1ik(ni,nj,k) +
     .        (at(k)*vb(ni)+vt(k)*ab(ni))*fcb1ik(ni,nj,k) +
     .        (vt(k)*vb(ni)-at(k)*ab(ni))*fcb2ik(ni,nj,k) +
     .        (at(k)*vb(ni)-vt(k)*ab(ni))*fcb3ik(ni,nj,k)) )
      enddo

      NS_grbchar = NS_grbchar +
     .         alsbot(ni,nj)*fcb1(ni,nj) + aksbot(ni,nj)*fcb2(ni,nj)

      NS_grbchar = NS_grbchar + (-1d0)**ni*(del(1,ni)*alsbot(2,nj)+
     .     del(2,ni)*alsbot(1,nj))/(asb1**2-asb2**2)*(
     .     4d0*amqb*mgluino*DDCOS(2d0*theb)*
     .     SD_B02(gmsb(ni)**2,amqb,mgluino,amuv**2) +
     .     DDCOS(2d0*theb)*DDSIN(2d0*theb)*
     .     (NS_A01(asb2**2,amuv**2)-NS_A01(asb1**2,amuv**2)))

      if(ni.eq.1) then
         NS_grbchar = NS_grbchar + 1d0/2d0*alsbot(1,nj)*(
     .        NS_delztl(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb1,theb,amuv,lamv,1) ) +
     .        1d0/(dsqrt(2d0)*amw*DDCOS(bet))*uu(nj,2)*DDSIN(theb)*
     .        NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv)*
     .        runmb +
     .        uu(nj,1)*DDSIN(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) +
     .        runmb*uu(nj,2)/dsqrt(2d0)/amw/DDCOS(bet)*DDCOS(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv)
      elseif(ni.eq.2) then
         NS_grbchar = NS_grbchar + 1d0/2d0*alsbot(2,nj)*(
     .        NS_delztl(amqt,mgluino,ast1,ast2,thet,amuv,lamv) +
     .        NS_delzst(amqb,mgluino,asb2,theb,amuv,lamv,2) ) +
     .        1d0/(dsqrt(2d0)*amw*DDCOS(bet))*uu(nj,2)*DDCOS(theb)*
     .        NS_delmtdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv,lamv)*
     .        runmb +
     .        uu(nj,1)*DDCOS(theb)*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv) +
     .        runmb*uu(nj,2)/dsqrt(2d0)/amw/DDCOS(bet)*(-DDSIN(theb))*
     .        NS_delthdiv(amqb,mgluino,asb1,asb2,theb,amuvdiv)
      endif

      NS_grbchar = (-1d0)*NS_grbchar

      return

      end

c -------------------------------------------------------------------- c

      subroutine NS_ftFUNCTIONs(fnt1,fnt2,fnt1ik,fnt2ik,fnt3ik,fnt4ik,
     .     fnt5ik,fnt6ik,fnt7ik,fct1,fct2,fct1ik,fct2ik,fct3ik,fct4ik,
     .     fct5ik,fct6ik,fct7ik)
*
      IMPLICIT NONE
      integer k,idec,I,J
      DOUBLE COMPLEX SD_C03,SD_C0_lam
*
      DOUBLE PRECISION fnt1(2,5),fnt2(2,5),fct1(2,2),fct2(2,2),
     .     fnt1ik(2,5,2),fnt2ik(2,5,2),fnt3ik(2,5,2),fnt4ik(2,5,2),
     .     fnt5ik(2,5,2),fnt6ik(2,5,2),fnt7ik(2,5,2),
     .     fct1ik(2,2,2),fct2ik(2,2,2),fct3ik(2,2,2),fct4ik(2,2,2),
     .     fct5ik(2,2,2),fct6ik(2,2,2),fct7ik(2,2,2)
      DOUBLE PRECISION gmst(2),gmsb(2)
      DOUBLE PRECISION amqt,amqb,amuv,lamv
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION SD_B02,NS_C1_lam,SD_C2_lam,NS_C1,SD_C2
*
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_decindex/idec
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
*
      gmst(1) = ast1
      gmst(2) = ast2
      gmsb(1) = asb1
      gmsb(2) = asb2

      amqt = amt
      amqb = amb

      do i=1,2
         do j=1,5
            fnt1(i,j) = SD_B02(gmst(i)**2,lamv,gmst(i),amuv**2) +
     .        2d0*amqt**2*dreal(SD_C0_lam(amqt,gmst(i),amneut(j),lamv))
     .         -2d0*gmst(i)**2*(
     .         NS_C1_lam(amqt,gmst(i),amneut(j),amqt,lamv,gmst(i),amuv,
     .                lamv)-
     .         SD_C2_lam(amqt,gmst(i),amneut(j),amqt,lamv,gmst(i),amuv,
     .                lamv))
     .         +2d0*amneut(j)**2*
     .         NS_C1_lam(amqt,gmst(i),amneut(j),amqt,lamv,gmst(i),amuv,
     .                lamv)

            fnt2(i,j) = -2d0*amqt*xmneut(j)*(
     .         dreal(SD_C0_lam(amqt,gmst(i),amneut(j),lamv)) +
     .         NS_C1_lam(amqt,gmst(i),amneut(j),amqt,lamv,gmst(i),amuv,
     .                lamv) )

            do k=1,2
               fnt1ik(i,j,k) = mgluino*xmneut(j)*(
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .                             gmst(k),mgluino,amqt)) +
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) )
               fnt2ik(i,j,k) = xmneut(j)*(amqt*(
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) +
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) + amqt*
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) )
               fnt3ik(i,j,k) = xmneut(j)*(-amqt*(
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) +
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) + amqt*
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) )
               fnt4ik(i,j,k) = mgluino*(amqt*
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) + amqt*(
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) )
               fnt5ik(i,j,k) = mgluino*(amqt*
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) - amqt*(
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) )
               fnt6ik(i,j,k) = gmst(k)**2*
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) + amqt*amqt*(
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) +
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) + amqt**2*(
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) + amneut(j)**2*
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) +
     .              SD_B02(gmst(i)**2,mgluino,amqt,amuv**2)
               fnt7ik(i,j,k) = gmst(k)**2*
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) - amqt*amqt*(
     .              dreal(SD_C03(amqt**2,gmst(i)**2,amneut(j)**2,
     .              gmst(k),mgluino,amqt)) +
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) + amqt**2*(
     .              NS_C1(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) ) + amneut(j)**2*
     .              SD_C2(amqt,gmst(i),amneut(j),gmst(k),mgluino,
     .                      amqt,amuv) +
     .              SD_B02(gmst(i)**2,mgluino,amqt,amuv**2)
            enddo
         enddo
      enddo

      do i=1,2
         do j=1,2
            fct1(i,j) = SD_B02(gmst(i)**2,lamv,gmst(i),amuv**2) +
     .       2d0*amqb**2*dreal(SD_C0_lam(amqb,gmst(i),amchar(j),lamv))
     .         -2d0*gmst(i)**2*(
     .         NS_C1_lam(amqb,gmst(i),amchar(j),amqb,lamv,gmst(i),amuv,
     .                lamv) -
     .         SD_C2_lam(amqb,gmst(i),amchar(j),amqb,lamv,gmst(i),amuv,
     .                lamv))
     .         +2d0*amchar(j)**2*
     .         NS_C1_lam(amqb,gmst(i),amchar(j),amqb,lamv,gmst(i),amuv,
     .                lamv)

            fct2(i,j) = -2d0*amqb*xmchar(j)*(
     .         dreal(SD_C0_lam(amqb,gmst(i),amchar(j),lamv)) +
     .         NS_C1_lam(amqb,gmst(i),amchar(j),amqb,lamv,gmst(i),amuv,
     .                lamv) )

            do k=1,2
               fct1ik(i,j,k) = mgluino*xmchar(j)*(
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) +
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) )
               fct2ik(i,j,k) = xmchar(j)*(amqb*(
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) +
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) + amqt*
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) )
               fct3ik(i,j,k) = xmchar(j)*(-amqb*(
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) +
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) + amqt*
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) )
               fct4ik(i,j,k) = mgluino*(amqt*
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) + amqb*(
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) )
               fct5ik(i,j,k) = mgluino*(amqt*
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) - amqb*(
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) )
               fct6ik(i,j,k) = gmsb(k)**2*
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) + amqt*amqb*(
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) +
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) + amqb**2*(
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) + amchar(j)**2*
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) +
     .              SD_B02(gmst(i)**2,mgluino,amqt,amuv**2)
               fct7ik(i,j,k) = gmsb(k)**2*
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) - amqt*amqb*(
     .              dreal(SD_C03(amqb**2,gmst(i)**2,amchar(j)**2,
     .              gmsb(k),mgluino,amqt)) +
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) + amqb**2*(
     .              NS_C1(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) -
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) ) + amchar(j)**2*
     .              SD_C2(amqb,gmst(i),amchar(j),gmsb(k),mgluino,
     .                      amqt,amuv) +
     .              SD_B02(gmst(i)**2,mgluino,amqt,amuv**2)
            enddo
         enddo
      enddo

      end

C------------------------------------------------------------------------C

      subroutine NS_fbFUNCTIONs(fnb1,fnb2,fnb1ik,fnb2ik,fnb3ik,fnb4ik,
     .     fnb5ik,fnb6ik,fnb7ik,fcb1,fcb2,fcb1ik,fcb2ik,fcb3ik,fcb4ik,
     .     fcb5ik,fcb6ik,fcb7ik)
*
      IMPLICIT NONE
      integer k,I,J,idec
      DOUBLE COMPLEX SD_C03,SD_C0_lam
*
      DOUBLE PRECISION fnb1(2,5),fnb2(2,5),fcb1(2,2),fcb2(2,2),
     .     fnb1ik(2,5,2),fnb2ik(2,5,2),fnb3ik(2,5,2),fnb4ik(2,5,2),
     .     fnb5ik(2,5,2),fnb6ik(2,5,2),fnb7ik(2,5,2),
     .     fcb1ik(2,2,2),fcb2ik(2,2,2),fcb3ik(2,2,2),fcb4ik(2,2,2),
     .     fcb5ik(2,2,2),fcb6ik(2,2,2),fcb7ik(2,2,2)
      DOUBLE PRECISION gmst(2),gmsb(2)
      DOUBLE PRECISION amuv,lamv
      DOUBLE PRECISION amqt,amqb
      DOUBLE PRECISION SD_B02,NS_C1_lam,SD_C2_lam,NS_C1,SD_C2
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION runmt,runmb,rmtauc
*
      COMMON/NS_decindex/idec
      COMMON/NS_qcdscales/amuv,lamv
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SMSPEC/MS,MC,MB,AMB,AMT,AMTAU,MMUON,MZ,MW
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
*
      gmst(1) = ast1
      gmst(2) = ast2
      gmsb(1) = asb1
      gmsb(2) = asb2

      amqt = amt
      amqb = amb

      do i=1,2
         do j=1,5
            fnb1(i,j) = SD_B02(gmsb(i)**2,lamv,gmsb(i),amuv**2) +
     .       2d0*amqb**2*dreal(SD_C0_lam(amqb,gmsb(i),amneut(j),lamv))
     .         -2d0*gmsb(i)**2*(
     .         NS_C1_lam(amqb,gmsb(i),amneut(j),amqb,lamv,gmsb(i),amuv,
     .                lamv)-
     .         SD_C2_lam(amqb,gmsb(i),amneut(j),amqb,lamv,gmsb(i),amuv,
     .                lamv))
     .         +2d0*amneut(j)**2*
     .         NS_C1_lam(amqb,gmsb(i),amneut(j),amqb,lamv,gmsb(i),amuv,
     .                lamv)

            fnb2(i,j) = -2d0*amqb*xmneut(j)*(
     .         dreal(SD_C0_lam(amqb,gmsb(i),amneut(j),lamv)) +
     .         NS_C1_lam(amqb,gmsb(i),amneut(j),amqb,lamv,gmsb(i),amuv,
     .                lamv) )

            do k=1,2
               fnb1ik(i,j,k) = mgluino*xmneut(j)*(
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .                             gmsb(k),mgluino,amqb)) +
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) )
               fnb2ik(i,j,k) = xmneut(j)*(amqb*(
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) +
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) + amqb*
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) )
               fnb3ik(i,j,k) = xmneut(j)*(-amqb*(
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) +
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) + amqb*
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) )
               fnb4ik(i,j,k) = mgluino*(amqb*
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) + amqb*(
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) )
               fnb5ik(i,j,k) = mgluino*(amqb*
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) - amqb*(
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) )
               fnb6ik(i,j,k) = gmsb(k)**2*
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) + amqb*amqb*(
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) +
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) + amqb**2*(
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) + amneut(j)**2*
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) +
     .              SD_B02(gmsb(i)**2,mgluino,amqb,amuv**2)
               fnb7ik(i,j,k) = gmsb(k)**2*
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) - amqb*amqb*(
     .              dreal(SD_C03(amqb**2,gmsb(i)**2,amneut(j)**2,
     .              gmsb(k),mgluino,amqb)) +
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) + amqb**2*(
     .              NS_C1(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) ) + amneut(j)**2*
     .              SD_C2(amqb,gmsb(i),amneut(j),gmsb(k),mgluino,
     .                      amqb,amuv) +
     .              SD_B02(gmsb(i)**2,mgluino,amqb,amuv**2)
            enddo
         enddo
      enddo

      do i=1,2
         do j=1,2
            fcb1(i,j) = SD_B02(gmsb(i)**2,lamv,gmsb(i),amuv**2) +
     .       2d0*amqt**2*dreal(SD_C0_lam(amqt,gmsb(i),amchar(j),lamv))
     .         -2d0*gmsb(i)**2*(
     .         NS_C1_lam(amqt,gmsb(i),amchar(j),amqt,lamv,gmsb(i),amuv,
     .                lamv) -
     .         SD_C2_lam(amqt,gmsb(i),amchar(j),amqt,lamv,gmsb(i),amuv,
     .                lamv))
     .         +2d0*amchar(j)**2*
     .         NS_C1_lam(amqt,gmsb(i),amchar(j),amqt,lamv,gmsb(i),amuv,
     .                lamv)

            fcb2(i,j) = -2d0*amqt*xmchar(j)*(
     .         dreal(SD_C0_lam(amqt,gmsb(i),amchar(j),lamv)) +
     .         NS_C1_lam(amqt,gmsb(i),amchar(j),amqt,lamv,gmsb(i),amuv,
     .                lamv) )

            do k=1,2
               fcb1ik(i,j,k) = mgluino*xmchar(j)*(
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) +
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) )
               fcb2ik(i,j,k) = xmchar(j)*(amqt*(
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) +
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) + amqb*
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) )
               fcb3ik(i,j,k) = xmchar(j)*(-amqt*(
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) +
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) + amqb*
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) )
               fcb4ik(i,j,k) = mgluino*(amqb*
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) + amqt*(
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) )
               fcb5ik(i,j,k) = mgluino*(amqb*
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) - amqt*(
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) )
               fcb6ik(i,j,k) = gmst(k)**2*
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) + amqb*amqt*(
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) +
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) + amqt**2*(
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) + amchar(j)**2*
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) +
     .              SD_B02(gmsb(i)**2,mgluino,amqb,amuv**2)
               fcb7ik(i,j,k) = gmst(k)**2*
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) - amqb*amqt*(
     .              dreal(SD_C03(amqt**2,gmsb(i)**2,amchar(j)**2,
     .              gmst(k),mgluino,amqb)) +
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) + amqt**2*(
     .              NS_C1(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) -
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) ) + amchar(j)**2*
     .              SD_C2(amqt,gmsb(i),amchar(j),gmst(k),mgluino,
     .                      amqb,amuv) +
     .              SD_B02(gmsb(i)**2,mgluino,amqb,amuv**2)
            enddo
         enddo
      enddo

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_delmtdiv(amq,mgluino,ast1,ast2,thet,
     .                                   amuv,lamv)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_sigmardiv,NS_sigmaldiv,NS_sigmasdiv

      NS_delmtdiv = 1d0/2d0*(
     .     NS_sigmardiv(amq,mgluino,ast1,ast2,thet,amuv,lamv) +
     .     NS_sigmaldiv(amq,mgluino,ast1,ast2,thet,amuv,lamv) ) +
     .     NS_sigmasdiv(amq,mgluino,ast1,ast2,thet,amuv,lamv)

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_delztr(amq,mgluino,ast1,ast2,thet,
     .                                 amuv,lamv)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_sigmar,NS_sigmalp,NS_sigmarp,NS_sigmasp

      NS_delztr = -NS_sigmar(amq,mgluino,ast1,ast2,thet,amuv,lamv) -
     .     amq**2*(NS_sigmalp(amq,mgluino,ast1,ast2,thet,amuv,lamv) +
     .     NS_sigmarp(amq,mgluino,ast1,ast2,thet,amuv,lamv) +
     .     2d0*NS_sigmasp(amq,mgluino,ast1,ast2,thet,amuv,lamv) )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_delztl(amq,mgluino,ast1,ast2,thet,
     .                                 amuv,lamv)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_sigmal,NS_sigmalp,NS_sigmarp,NS_sigmasp

      NS_delztl = -NS_sigmal(amq,mgluino,ast1,ast2,thet,amuv,lamv) -
     .     amq**2*(NS_sigmalp(amq,mgluino,ast1,ast2,thet,amuv,lamv) +
     .     NS_sigmarp(amq,mgluino,ast1,ast2,thet,amuv,lamv) +
     .     2d0*NS_sigmasp(amq,mgluino,ast1,ast2,thet,amuv,lamv) )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_delzst(amqt,mgluino,amsq,thet,
     .                                 amuv,lamv,ni)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_B1,SD_BP1,SD_BP02,SD_B02
      DOUBLE PRECISION DDCOS,DDSIN
      integer ni

      NS_delzst = 2d0*(-2d0*NS_B1(amsq**2,amsq,lamv,amuv**2)
     .     -2d0*amsq**2*SD_BP1(amsq**2,amsq,lamv,amuv**2) +
     .     (amqt**2+mgluino**2-amsq**2)*
     .     SD_BP02(amsq**2,amqt,mgluino,amuv**2) -
     .     SD_B02(amsq**2,amqt,mgluino,amuv**2) +
     .     (-1d0)**ni*2d0*DDSIN(2d0*thet)*amqt*mgluino*
     .     SD_BP02(amsq**2,amqt,mgluino,amuv**2) )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_sigmar(amq,mgluino,ast1,ast2,thet,
     .                                 amuv,lamv)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_B1
      DOUBLE PRECISION DDCOS,DDSIN

      NS_sigmar = -(2d0*NS_B1(amq**2,amq,lamv,amuv**2)+
     .     (1d0-DDCOS(2d0*thet))*NS_B1(amq**2,mgluino,ast1,amuv**2)
     .    +(1d0+DDCOS(2d0*thet))*NS_B1(amq**2,mgluino,ast2,amuv**2)
     .     )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_sigmardiv(amq,mgluino,ast1,ast2,thet,
     .                                    amuv,lamv)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_B1_DIV
      DOUBLE PRECISION DDCOS,DDSIN

      NS_sigmardiv = -(2d0*NS_B1_DIV(amq**2,amq,lamv,amuv**2)+
     .     (1d0-DDCOS(2d0*thet))*NS_B1_DIV(amq**2,mgluino,ast1,amuv**2)
     .    +(1d0+DDCOS(2d0*thet))*NS_B1_DIV(amq**2,mgluino,ast2,amuv**2)
     .     )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_sigmarp(amq,mgluino,ast1,ast2,thet,
     .                                 amuv,lamv)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION SD_BP1
      DOUBLE PRECISION DDCOS,DDSIN

      NS_sigmarp = -(2d0*SD_BP1(amq**2,amq,lamv,amuv**2)+
     .     (1d0-DDCOS(2d0*thet))*SD_BP1(amq**2,mgluino,ast1,amuv**2)
     .    +(1d0+DDCOS(2d0*thet))*SD_BP1(amq**2,mgluino,ast2,amuv**2)
     .     )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_sigmal(amq,mgluino,ast1,ast2,thet,
     .                                 amuv,lamv)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_B1
      DOUBLE PRECISION DDCOS,DDSIN

      NS_sigmal = -(2d0*NS_B1(amq**2,amq,lamv,amuv**2)+
     .     (1d0+DDCOS(2d0*thet))*NS_B1(amq**2,mgluino,ast1,amuv**2)
     .    +(1d0-DDCOS(2d0*thet))*NS_B1(amq**2,mgluino,ast2,amuv**2)
     .     )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_sigmaldiv(amq,mgluino,ast1,ast2,thet,
     .                                    amuv,lamv)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_B1_DIV
      DOUBLE PRECISION DDCOS,DDSIN

      NS_sigmaldiv = -(2d0*NS_B1_DIV(amq**2,amq,lamv,amuv**2)+
     .     (1d0+DDCOS(2d0*thet))*NS_B1_DIV(amq**2,mgluino,ast1,amuv**2)
     .    +(1d0-DDCOS(2d0*thet))*NS_B1_DIV(amq**2,mgluino,ast2,amuv**2)
     .     )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_sigmalp(amq,mgluino,ast1,ast2,thet,
     .                                 amuv,lamv)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION SD_BP1
      DOUBLE PRECISION DDCOS,DDSIN

      NS_sigmalp = -(2d0*SD_BP1(amq**2,amq,lamv,amuv**2)+
     .     (1d0+DDCOS(2d0*thet))*SD_BP1(amq**2,mgluino,ast1,amuv**2)
     .    +(1d0-DDCOS(2d0*thet))*SD_BP1(amq**2,mgluino,ast2,amuv**2)
     .     )

      return

      end


c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_sigmasdiv(amq,mgluino,ast1,ast2,thet,
     .                                    amuv,lamv)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION SD_B02_DIV
      DOUBLE PRECISION DDCOS,DDSIN

      NS_sigmasdiv = -(4d0*SD_B02_DIV(amq**2,amq,lamv,amuv**2)+
     .     mgluino/amq*DDSIN(2d0*thet)*(
     .     SD_B02_DIV(amq**2,mgluino,ast1,amuv**2)
     .    -SD_B02_DIV(amq**2,mgluino,ast2,amuv**2)
     .     ) )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_sigmasp(amq,mgluino,ast1,ast2,thet,
     .                                 amuv,lamv)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION SD_BP02
      DOUBLE PRECISION DDCOS,DDSIN

      NS_sigmasp = -(4d0*SD_BP02(amq**2,amq,lamv,amuv**2)+
     .     mgluino/amq*DDSIN(2d0*thet)*(
     .     SD_BP02(amq**2,mgluino,ast1,amuv**2)
     .    -SD_BP02(amq**2,mgluino,ast2,amuv**2)
     .     ) )

      return

      end


c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_delthdiv(amqt,mgluino,ast1,ast2,
     .                                   thet,amuv)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION SD_B02_DIV
      DOUBLE PRECISION DDCOS,DDSIN

      NS_delthdiv = 1d0/(ast1**2-ast2**2)*(4d0*amqt*mgluino*
     .     DDCOS(2d0*thet)*SD_B02_DIV(ast2**2,amqt,mgluino,amuv**2) +
     .     DDCOS(2d0*thet)*DDSIN(2d0*thet)*
     .     (ast2**2-ast1**2)*dlog(amuv**2) )

      return

      end

c -------------------------------------------------------------------- c
c ----------------------- The real corrections ----------------------- c

      DOUBLE PRECISION FUNCTION NS_corrreali(amq,mcharneut,amsti,lamv,
     .     icharneut,isign,ni,nj,idec)
*
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)
      INTEGER ni,nj
      DIMENSION atopr(2,5),btopr(2,5),alstor(2,2),akstor(2,2),
     .     abot(2,5),bbot(2,5),alsbot(2,2),aksbot(2,2)
      DOUBLE COMPLEX NS_ccspen,NS_kappa
*
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_neutstoptop/atopr,btopr
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor
      COMMON/NS_neutsbotbot/abot,bbot
      EXTERNAL NS_ccspen,NS_kappa

      kap = dreal(NS_kappa(amsti**2,amq**2,mcharneut**2,0d0))

      b0 = (amsti**2-amq**2-mcharneut**2+kap)/(2d0*amq*mcharneut)
      b1 = (amsti**2-amq**2+mcharneut**2-kap)/(2d0*amsti*mcharneut)
      b2 = (amsti**2+amq**2-mcharneut**2-kap)/(2d0*amsti*amq)

      lb0 = dreal(cdlog(dcmplx(b0)))
      lb1 = dreal(cdlog(dcmplx(b1)))
      lb2 = dreal(cdlog(dcmplx(b2)))
      lb12 = dreal(cdlog(dcmplx(b1/b2)))
      lb02 = dreal(cdlog(dcmplx(b0/b2)))

      hi00 = 1d0/4d0/amsti**4*(kap*
     .     dlog(kap**2/(lamv*amsti*amq*mcharneut)) - kap -
     .     (amq**2-mcharneut**2)*lb12 - amsti**2*lb0 )
      hi11 = 1d0/4d0/(amq**2*amsti**2)*(kap*
     .     dlog(kap**2/(lamv*amsti*amq*mcharneut)) - kap -
     .     (amsti**2-mcharneut**2)*lb02 - amq**2*lb1)
      hi01 = dreal(1d0/(4d0*amsti**2)*(-2d0*
     .     dlog((lamv*amsti*amq*mcharneut)/kap**2)*lb2 +
     .     2d0*lb2**2 - lb0**2 - lb1**2 +
     .     2d0*NS_ccspen(dcmplx(1d0-b2**2)) -
     .     NS_ccspen(dcmplx(1-b0**2))
     .     - NS_ccspen(dcmplx(1-b1**2)) ) )
      hi = 1d0/(4d0*amsti**2)*(kap/2d0*
     .     (amsti**2+amq**2+mcharneut**2) + 2d0*amsti**2*amq**2*
     .     lb2 + 2d0*amsti**2*mcharneut**2*lb1 +
     .     2d0*amq**2*mcharneut**2*lb0)
      hi0 = 1d0/(4d0*amsti**2)*(-2d0*amq**2*lb2
     .     -2d0*mcharneut**2*lb1-kap)
      hi1 = 1d0/(4d0*amsti**2)*(-2d0*amsti**2*lb2
     .     -2d0*mcharneut**2*lb0-kap)
      hi10 = 1d0/(4d0*amsti**2)*(amsti**4*lb2-
     .     mcharneut**2*(2d0*amq**2-2d0*amsti**2+mcharneut**2)*
     .     lb0 - kap/4d0*(amq**2-3d0*amsti**2+
     .     5d0*mcharneut**2) )

      if(icharneut.eq.1.and.idec.eq.1) then
         cli = - btopr(ni,nj)
         cri = - atopr(ni,nj)
      elseif(icharneut.eq.2.and.idec.eq.1) then
         cli = - akstor(ni,nj)
         cri = - alstor(ni,nj)
      elseif(icharneut.eq.1.and.idec.eq.2) then
         cli = - bbot(ni,nj)
         cri = - abot(ni,nj)
      elseif(icharneut.eq.2.and.idec.eq.2) then
         cli = - aksbot(ni,nj)
         cri = - alsbot(ni,nj)
      endif

      if(isign.eq.1) then
         epschi = -1d0
      elseif(isign.eq.0) then
         epschi = 1d0
      endif

      NS_corrreali = 8d0*cli*cri*amq*mcharneut*epschi*((amsti**2+amq**2
     .     -mcharneut**2)*hi01+amsti**2*hi00+amq**2*hi11+hi0+hi1)
     .     +(cli**2+cri**2)*(2d0*(amq**2+mcharneut**2-amsti**2)*
     .     (amsti**2*hi00+amq**2*hi11+hi0+hi1) + 2d0*(amq**4-
     .     (mcharneut**2-amsti**2)**2)*hi01-hi-hi10)

      return

      end

c -------------------------------------------------------------------- c
c ----------------------- The real corrections ----------------------- c
c ----------- for the processes gaugino -> squark + quark ------------ c

      DOUBLE PRECISION FUNCTION NS_realicorr(amq,mcharneut,amsti,lamv,
     .     icharneut,isign,ni,nj,idec)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)
*
      INTEGER ni,nj
      DOUBLE PRECISION atopr(2,5),btopr(2,5),alstor(2,2),akstor(2,2),
     .     abot(2,5),bbot(2,5),alsbot(2,2),aksbot(2,2)

      DOUBLE COMPLEX NS_ccspen,NS_kappa
*
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/NS_runmcalc/rmtc,rmbc,rmtauc
      COMMON/NS_neutstoptop/atopr,btopr
      COMMON/NS_charsbottop/alsbot,aksbot
      COMMON/NS_charstopbot/alstor,akstor
      COMMON/NS_neutsbotbot/abot,bbot
*
      EXTERNAL NS_kappa

      kap = dreal(NS_kappa(amsti**2,amq**2,mcharneut**2,0d0))

      b0 = (mcharneut**2-amq**2-amsti**2+kap)/(2d0*amq*amsti)
      b1 = (mcharneut**2-amq**2+amsti**2-kap)/(2d0*amsti*mcharneut)
      b2 = (mcharneut**2+amq**2-amsti**2-kap)/(2d0*mcharneut*amq)

      lb0 = dreal(cdlog(dcmplx(b0)))
      lb1 = dreal(cdlog(dcmplx(b1)))
      lb2 = dreal(cdlog(dcmplx(b2)))
      lb12 = dreal(cdlog(dcmplx(b1/b2)))
      lb01 = dreal(cdlog(dcmplx(b0/b1)))
      lb02 = dreal(cdlog(dcmplx(b0/b2)))

      hi11 = 1d0/4d0/(amq**2*mcharneut**2)*(kap*
     .     dlog(kap**2/(lamv*amsti*amq*mcharneut)) - kap -
     .     (mcharneut**2-amsti**2)*lb02 - amq**2*lb1)
      hi22 = 1d0/4d0/(mcharneut**2*amsti**2)*(kap*
     .     dlog(kap**2/(lamv*amsti*amq*mcharneut)) - kap -
     .     (mcharneut**2-amq**2)*lb01 - amsti**2*lb2)
      hi21 = dreal(1d0/(4d0*mcharneut**2)*(-2d0*
     .     dlog((lamv*amsti*amq*mcharneut)/kap**2)*lb0 +
     .     2d0*lb0**2 - lb2**2 - lb1**2 +
     .     2d0*NS_ccspen(dcmplx(1d0-b0**2)) -
     .     NS_ccspen(dcmplx(1-b2**2))
     .     - NS_ccspen(dcmplx(1-b1**2)) ) )
      hi = 1d0/(4d0*mcharneut**2)*(kap/2d0*
     .     (mcharneut**2+amq**2+amsti**2) + 2d0*mcharneut**2*amq**2*
     .     lb2 + 2d0*amsti**2*mcharneut**2*lb1 +
     .     2d0*amq**2*amsti**2*lb0)
      hi2 = 1d0/(4d0*mcharneut**2)*(-2d0*mcharneut**2*lb1
     .     -2d0*amq**2*lb0-kap)
      hi1 = 1d0/(4d0*mcharneut**2)*(-2d0*mcharneut**2*lb2
     .     -2d0*amsti**2*lb0-kap)
      hi12 = 1d0/(4d0*mcharneut**2)*(amsti**4*lb0-
     .     mcharneut**2*(2d0*amq**2-2d0*amsti**2+mcharneut**2)*
     .     lb2 - kap/4d0*(amq**2-3d0*amsti**2+
     .     5d0*mcharneut**2) )

      if(icharneut.eq.1.and.idec.eq.1) then
         cli = - btopr(ni,nj)
         cri = - atopr(ni,nj)
      elseif(icharneut.eq.2.and.idec.eq.1) then
         cli = - akstor(ni,nj)
         cri = - alstor(ni,nj)
      elseif(icharneut.eq.1.and.idec.eq.2) then
         cli = - bbot(ni,nj)
         cri = - abot(ni,nj)
      elseif(icharneut.eq.2.and.idec.eq.2) then
         cli = - aksbot(ni,nj)
         cri = - alsbot(ni,nj)
      endif

      if(isign.eq.1) then
         epschi = -1d0
      elseif(isign.eq.0) then
         epschi = 1d0
      endif

      NS_realicorr = 8d0*cli*cri*amq*mcharneut*epschi*((amsti**2+amq**2
     .     -mcharneut**2)*hi21+amsti**2*hi22+amq**2*hi11+hi2+hi1)
     .     +(cli**2+cri**2)*(2d0*(amq**2+mcharneut**2-amsti**2)*
     .     (amsti**2*hi22+amq**2*hi11+hi2+hi1) + 2d0*(amq**4-
     .     (mcharneut**2-amsti**2)**2)*hi21-hi-hi12)

      return

      end

c -------------------------------------------------------------------- c
c ---------- Beenakker, Hoepker and Zerwas, hep-ph/9602378 ----------- c
c -------------------------------------------------------------------- c

c - QCD corrections to the decays light squark -> light quark + gluino c
      DOUBLE PRECISION FUNCTION NS_gama(r)
*
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE COMPLEX NS_ccspen
*
      pi = 4d0*datan(1d0)

      NS_gama = 3d0/(r-1d0)*dreal(NS_ccspen(dcmplx(1d0-r))) -
     .     r/(r-1d0)*dreal(NS_ccspen(dcmplx(-r)))
     .     + (5d0*r-6d0)/(12d0*(r-1d0))*pi**2 +
     .     59d0/24d0 + r/(4d0*(r-1d0)) +
     .     ( (3d0+r)/(2d0*(r-1d0))*dlog(r) - 2d0)*
     .     dlog(dabs(1d0-r)) +
     .     ( (r*(5d0*r-6d0))/(4d0*(r-1d0)**2) -
     .     r/(r-1d0)*dlog(1d0+r) )*dlog(r)

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gamfcap(r)
*
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE COMPLEX NS_ccspen
*
      pi = 4d0*datan(1d0)

      NS_gamfcap = -2d0/(r-1d0)*dreal(NS_ccspen(dcmplx(1d0-r))) +
     .     (2d0*r)/(r-1d0)*dreal(NS_ccspen(dcmplx(-r))) +
     .     (4d0-3d0*r)/(6d0*(r-1d0))*pi**2 + 5d0/2d0 - r/2d0 +
     .     ( r- r**2/2d0 - (r+1d0)/(r-1d0)*dlog(r) )*
     .     dlog(dabs(1d0-r)) +
     .     ( 2d0*r/(r-1d0)*dlog(1d0+r)-r+r**2/2d0)*dlog(r)

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gamf(r)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)

      NS_gamf = -3d0/(4d0*r) +
     .     ((r-1d0)*(r+3d0))/(4d0*r**2)*dlog(dabs(1d0-r))

      return

      end

c -------------------------------------------------------------------- c
c maggie changed with respect to the paper 26/3/03

      DOUBLE PRECISION FUNCTION NS_gamrendec(amsq,amst1,amst2,amt,amsb1,
     .     amsb2,amgl,scala)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)

      mur = scala

      NS_gamrendec = 3d0/4d0*dlog(mur**2/amsq**2) - 1d0/4d0

      NS_gamrendec =  NS_gamrendec + 4d0/12d0*dlog(mur**2/amsq**2)
     .     + 1d0/24d0*dlog(mur**2/amsb1**2) + 1d0/24d0*
     .     dlog(mur**2/amsb2**2)
     .     + 1d0/24d0*dlog(mur**2/amst1**2) + 1d0/24d0*
     .     dlog(mur**2/amst2**2) + 1d0/6d0*dlog(mur**2/amt**2)
     .     + 1d0/2d0*dlog(mur**2/amgl**2)
      return

      end

c end maggie changed

c -------------------------------------------------------------------- c
c ------ Beenakker, Hoepker, Plehn and Zerwas, hep-ph/9610313 -------- c
c -------------------------------------------------------------------- c

c -------- QCD corrections to the decay stop1/2 -> top gluino -------- c
c -------- and sbottom1/2 -> bottom gluino, gluino -> stop1/2 top ---- c
c -------- gluino -> sbottom1/2 bottom ------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gamtop1(amst1,amst2,amt,amgl,thet,
     .                                     ival,amuv,lamv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)
      DOUBLE PRECISION NS_A01,SD_B02,SD_BP02
      DOUBLE PRECISION isign
      DOUBLE PRECISION DDCOS,DDSIN

      pi = 4d0*datan(1d0)

      if(ival.eq.1) then
         isign = 1d0
      elseif(ival.eq.2) then
         isign = -1d0
      endif

      sig2t = isign*amt*amgl*DDSIN(2d0*thet)

      NS_gamtop1 = 2d0*NS_A01(amt**2,amuv**2) -2d0*amt**2 +
     .   2d0*NS_A01(amgl**2,amuv**2) - NS_A01(amst1**2,amuv**2)
     .   - NS_A01(amst2**2,amuv**2) +
     .   (amt**2+amst1**2-amgl**2)*SD_B02(amt**2,amgl,amst1,amuv**2)
     .  +(amt**2+amst2**2-amgl**2)*SD_B02(amt**2,amgl,amst2,amuv**2)
     .   -4d0*amt**2*sig2t*(SD_BP02(amt**2,amgl,amst1,amuv**2)
     .   - SD_BP02(amt**2,amgl,amst2,amuv**2) ) + 2d0*amt**2*(
     .   (amgl**2+amt**2-amst1**2)*SD_BP02(amt**2,amgl,amst1,amuv**2)
     .  +(amgl**2+amt**2-amst2**2)*SD_BP02(amt**2,amgl,amst2,amuv**2)
     .  -4d0*amt**2*SD_BP02(amt**2,lamv,amt,amuv**2) )

      NS_gamtop1 = -1d0/16d0/pi**2*NS_gamtop1

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gamtop2(amst1,amst2,amt,amgl,thet,
     .     ival,amuv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)
      DOUBLE PRECISION NS_A01,SD_B02
      DOUBLE PRECISION isign
      DOUBLE PRECISION DDCOS,DDSIN

      pi = 4d0*datan(1d0)

      if(ival.eq.1) then
         isign = 1d0
      elseif(ival.eq.2) then
         isign = -1d0
      endif

      NS_gamtop2 = NS_A01(amst2**2,amuv**2) -
     .  NS_A01(amst1**2,amuv**2)
     .  -(amgl**2+amt**2-amst1**2)*SD_B02(amt**2,amgl,amst1,amuv**2)
     .  +(amgl**2+amt**2-amst2**2)*SD_B02(amt**2,amgl,amst2,amuv**2)

      NS_gamtop2 = (isign*DDCOS(2d0*thet))**2*(amgl**2+amt**2-amst1**2)*
     .     NS_gamtop2

      NS_gamtop2 = -1d0/16d0/pi**2*NS_gamtop2

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gamglui1(amst1,amst2,amsq,amt,amgl,
     .     amuv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_A01,SD_BP02

      pi = 4d0*datan(1d0)

      NS_gamglui1 = -NS_A01(amsq**2,amuv**2)+(amsq**2+amgl**2)*
     .     SD_B02(amgl**2,amsq,0d0,amuv**2) + 2d0*amgl**2*
     .     (amgl**2-amsq**2)*SD_BP02(amgl**2,amsq,0d0,amuv**2)

      NS_gamglui1 = -1d0/16d0/pi**2*NS_gamglui1

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gamglui2(amst1,amst2,amt,thet,amsb1,
     .     amsb2,amb,theb,amgl,ival,amuv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)
      DOUBLE PRECISION isign
      DOUBLE PRECISION NS_A01,SD_B02,SD_BP02
      DOUBLE PRECISION DDCOS,DDSIN

      pi = 4d0*datan(1d0)

      if(ival.eq.1) then
         isign = 1d0
      elseif(ival.eq.2) then
         isign = -1d0
      endif

      sig2t = isign*amt*amgl*DDSIN(2d0*thet)
      sig2b = isign*amb*amgl*DDSIN(2d0*theb)

      NS_gamglui2 = 2d0*NS_A01(amt**2,amuv**2) -
     .     NS_A01(amst1**2,amuv**2) - NS_A01(amst2**2,amuv**2)
     .     + (amst1**2+amgl**2-amt**2)*
     .     SD_B02(amgl**2,amst1,amt,amuv**2) +
     .     (amst2**2+amgl**2-amt**2)*
     .     SD_B02(amgl**2,amst2,amt,amuv**2) -
     .     4d0*amgl**2*sig2t*( SD_BP02(amgl**2,amst1,amt,amuv**2)
     .     - SD_BP02(amgl**2,amst2,amt,amuv**2) )
     .     + 2d0*amgl**2*( (amgl**2+amt**2-amst1**2)*
     .     SD_BP02(amgl**2,amst1,amt,amuv**2) +
     .     (amgl**2+amt**2-amst2**2)*
     .     SD_BP02(amgl**2,amst2,amt,amuv**2) )

c maggie changed with respect to the paper 26/3/03
      NS_gamglui2 = NS_gamglui2 +
     .     2d0*NS_A01(amb**2,amuv**2) -
     .     NS_A01(amsb1**2,amuv**2) - NS_A01(amsb2**2,amuv**2)
     .     + (amsb1**2+amgl**2-amb**2)*
     .     SD_B02(amgl**2,amsb1,amb,amuv**2) +
     .     (amsb2**2+amgl**2-amb**2)*
     .     SD_B02(amgl**2,amsb2,amb,amuv**2) -
     .     4d0*amgl**2*sig2b*( SD_BP02(amgl**2,amsb1,amb,amuv**2)
     .     - SD_BP02(amgl**2,amsb2,amb,amuv**2) )
     .     + 2d0*amgl**2*( (amgl**2+amb**2-amsb1**2)*
     .     SD_BP02(amgl**2,amsb1,amb,amuv**2) +
     .     (amgl**2+amb**2-amsb2**2)*
     .     SD_BP02(amgl**2,amsb2,amb,amuv**2) )
c end maggie changed

      NS_gamglui2 = -1d0/16d0/pi**2*NS_gamglui2

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gamglui3(amgl,amuv,lamv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_A01,SD_BP02

      pi = 4d0*datan(1d0)

      NS_gamglui3 =  NS_A01(amgl**2,amuv**2) - amgl**2 -
     .     4d0*amgl**4*SD_BP02(amgl**2,lamv,amgl,amuv**2)

      NS_gamglui3 = -1d0/16d0/pi**2*NS_gamglui3

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gam11(amst1,amst2,amt,amgl,thet,ival,
     .     amuv,lamv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)
      DOUBLE PRECISION isign
      DOUBLE PRECISION SD_B02,SD_BP02
      DOUBLE PRECISION DDCOS,DDSIN

      pi = 4d0*datan(1d0)

      if(ival.eq.1) then
         isign = 1d0
      elseif(ival.eq.2) then
         isign = -1d0
      endif

      sig2t = isign*amt*amgl*DDSIN(2d0*thet)

      NS_gam11 = SD_B02(amst1**2,amgl,amt,amuv**2) -
     .     SD_B02(amst1**2,lamv,amst1,amuv**2) +
     .     2d0*sig2t*SD_BP02(amst1**2,amgl,amt,amuv**2) -
     .     (amgl**2+amt**2-amst1**2)*
     .     SD_BP02(amst1**2,amgl,amt,amuv**2) - 2d0*amst1**2*
     .     SD_BP02(amst1**2,lamv,amst1,amuv**2)

      NS_gam11 = -1d0/16d0/pi**2*NS_gam11

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gam12(amst1,amst2,amt,amgl,thet,ival,
     .     amuv,lamv,scala)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)
      DOUBLE PRECISION isign
      DOUBLE PRECISION NS_A01,SD_B02,NS_A01_DIV,SD_B02_DIV
      DOUBLE PRECISION DDCOS,DDSIN

      pi = 4d0*datan(1d0)

      if(ival.eq.1) then
         isign = 1d0
      elseif(ival.eq.2) then
         isign = -1d0
      endif

      sig2t = isign*amt*amgl*DDSIN(2d0*thet)

      NS_gam12 = NS_A01(amst2**2,amuv**2) - NS_A01(amst1**2,amuv**2)
     .     + 4d0*amt**2*amgl**2/sig2t*
     .     SD_B02(amst1**2,amgl,amt,amuv**2)

      NS_gam12 = NS_gam12 -
     .     ( NS_A01_DIV(amst2**2,(amuv/scala)**2)
     .     - NS_A01_DIV(amst1**2,(amuv/scala)**2)
     .     + 4d0*amt**2*amgl**2/sig2t*
     .       SD_B02_DIV(amst1**2,amgl,amt,(amuv/scala)**2) )

      NS_gam12 = 1d0/(amst1**2-amst2**2)*sig2t*DDCOS(2d0*thet)**2*
     .     NS_gam12

      NS_gam12 = -1d0/16d0/pi**2*NS_gam12

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gamvirt(amst1,amst2,amtop,amgl,thet,
     .     ival,amuv,lamv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)
      DOUBLE PRECISION isign
      DOUBLE PRECISION DDCOS,DDSIN

      DIMENSION fffunc(3),fafunc(3)

      COMMON/NS_qcdscales/amuvv,lamvv
      COMMON/NS_relmasses/mst1,mst2,mgl,mtop

      amuvv = amuv
      lamvv = lamv
      mst1  = amst1
      mst2  = amst2
      mgl   = amgl
      mtop  = amtop

      pi = 4d0*datan(1d0)

      if(ival.eq.1) then
         isign = 1d0
      elseif(ival.eq.2) then
         isign = -1d0
      endif

      sig2t = isign*amtop*amgl*DDSIN(2d0*thet)

      CALL NS_fifaFUNCTIONs(fffunc,fafunc,fbfunc)

      NS_gamvirt = 64d0/9d0*pi*
     .     ( fffunc(1) + sig2t*fffunc(2) + sig2t**2*fffunc(3) ) +
     .     8d0*pi*
     .     (fafunc(1) + sig2t*fafunc(2) + sig2t**2*fafunc(3)) +
     .     16d0/3d0*pi*
     .     (-(amgl**2+mtop**2-amst1**2)+2d0*sig2t)*fbfunc

      return

      end

c -------------------------------------------------------------------- c
c --- this FUNCTION is for the stop1/2, sbottom1/2 decays ---

      DOUBLE PRECISION FUNCTION NS_gamreal(amst1,amt,amgl,thet,ival,
     .                                     lamv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)
      DOUBLE PRECISION ist1gl,ist1t,iglgl,ist1st1,itt,igl,ist1,it,
     .     itgl,iglst1
      DOUBLE PRECISION isign
      DOUBLE PRECISION DDCOS,DDSIN
      DOUBLE COMPLEX NS_kappa,NS_ccspen

      pi = 4d0*datan(1d0)

      if(ival.eq.1) then
         isign = 1d0
      elseif(ival.eq.2) then
         isign = -1d0
      endif

      sig2t = isign*amt*amgl*DDSIN(2d0*thet)

      m0 = amst1
      m1 = amt
      m2 = amgl

      kap = dreal(NS_kappa(m0**2,m1**2,m2**2,0d0))

      b0 = (m0**2-m1**2-m2**2+kap)/(2d0*m1*m2)
      b1 = (m0**2-m1**2+m2**2-kap)/(2d0*m0*m2)
      b2 = (m0**2+m1**2-m2**2-kap)/(2d0*m0*m1)

      ist1gl = dreal(
     .     1d0/(4d0*m0**2)*(-2d0*dlog((lamv*m0*m1*m2)/kap**2)*
     .     dlog(b1) + 2d0*(dlog(b1))**2 - (dlog(b0))**2 -
     .     (dlog(b2))**2 + 2d0*NS_ccspen(dcmplx(1d0-b1**2))
     .     - NS_ccspen(dcmplx(1d0-b0**2)) -
     .     NS_ccspen(dcmplx(1d0-b2**2)) ) )

      ist1t = dreal(
     .     1d0/(4d0*m0**2)*(-2d0*dlog((lamv*m0*m1*m2)/kap**2)*
     .     dlog(b2) + 2d0*(dlog(b2))**2 - (dlog(b0))**2 -
     .     (dlog(b1))**2 + 2d0*NS_ccspen(dcmplx(1d0-b2**2))
     .     - NS_ccspen(dcmplx(1d0-b0**2)) -
     .     NS_ccspen(dcmplx(1d0-b1**2)) ) )

      iglgl = 1d0/(4d0*m2**2*m0**2)*(kap*dlog(kap**2/(lamv*m0*m1*m2))
     .     -kap-(m0**2-m1**2)*dlog(b0/b1)-m2**2*dlog(b2) )

      ist1st1 = 1d0/(4d0*m0**4)*(kap*dlog(kap**2/(lamv*m0*m1*m2))
     .     -kap-(m1**2-m2**2)*dlog(b1/b2)-m0**2*dlog(b0))

      itt = 1d0/(4d0*m1**2*m0**2)*(kap*dlog(kap**2/(lamv*m0*m1*m2))
     .     -kap-(m0**2-m2**2)*dlog(b0/b2)-m1**2*dlog(b1) )

      igl = 1d0/(4d0*m0**2)*(-2d0*m0**2*dlog(b1)-2d0*m1**2*dlog(b0)
     .     -kap)

      ist1 = 1d0/(4d0*m0**2)*(-2d0*m1**2*dlog(b2)-
     .     2d0*m2**2*dlog(b1)-kap)

      it = 1d0/(4d0*m0**2)*(-2d0*m0**2*dlog(b2)-2d0*m2**2*dlog(b0)
     .     -kap)

      itgl = -1d0/(4d0*m0**2)*(-m2**4*dlog(b0)+m0**2*(2d0*m1**2
     .     -2d0*m2**2)*dlog(b2) + m0**4*dlog(b2) +
     .     kap/4d0*(m1**2+5d0*m0**2-3d0*m2**2) )

      iglst1 = 1d0/(4d0*m0**2)*(m0**4*dlog(b1)-m1**2*(2d0*m2**2
     .     -2d0*m0**2+m1**2)*dlog(b0) -kap/4d0*(m2**2-3d0*m0**2
     .     +5d0*m1**2) )

      NS_gamreal = 8d0/pi/amst1*(-(amgl**2+amt**2-amst1**2)+
     .     2d0*sig2t)*(-(amst1**2-amt**2)*ist1gl+amt**2*ist1t
     .     -amgl**2*iglgl-igl) +
     .     32d0/9d0/pi/amst1*(-(amgl**2+amt**2-amst1**2)+
     .     2d0*sig2t)*(-amst1**2*ist1st1-amt**2*itt-
     .     (amt**2+amst1**2-amgl**2)*ist1t-ist1-it) +
     .     4d0/3d0/pi/amst1*(4d0/3d0*itgl-3d0*iglst1)

      return

      end

c -------------------------------------------------------------------- c
c --- this FUNCTION is for the gluino decays ---

      DOUBLE PRECISION FUNCTION NS_gamrealgl(amst1,amt,amgl,thet,ival,
     .                                       lamv)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)
      DOUBLE PRECISION ist1gl,ist1t,iglgl,ist1st1,itt,igl,ist1,it,
     .     itgl,iglst1
      DOUBLE PRECISION isign
      DOUBLE PRECISION DDCOS,DDSIN
      DOUBLE COMPLEX NS_kappa,NS_ccspen

      pi = 4d0*datan(1d0)

      if(ival.eq.1) then
         isign = 1d0
      elseif(ival.eq.2) then
         isign = -1d0
      endif

      sig2t = isign*amt*amgl*DDSIN(2d0*thet)

      m0 = amgl
      m1 = amt
      m2 = amst1

      kap = dreal(NS_kappa(m0**2,m1**2,m2**2,0d0))

      b0 = (m0**2-m1**2-m2**2+kap)/(2d0*m1*m2)
      b1 = (m0**2-m1**2+m2**2-kap)/(2d0*m0*m2)
      b2 = (m0**2+m1**2-m2**2-kap)/(2d0*m0*m1)

      ist1gl = dreal(
     .     1d0/(4d0*m0**2)*(-2d0*dlog((lamv*m0*m1*m2)/kap**2)*
     .     dlog(b1) + 2d0*(dlog(b1))**2 - (dlog(b0))**2 -
     .     (dlog(b2))**2 + 2d0*NS_ccspen(dcmplx(1d0-b1**2))
     .     - NS_ccspen(dcmplx(1d0-b0**2)) -
     .     NS_ccspen(dcmplx(1d0-b2**2)) ) )

      ist1t = dreal(
     .     1d0/(4d0*m0**2)*(-2d0*dlog((lamv*m0*m1*m2)/kap**2)*
     .     dlog(b0) + 2d0*(dlog(b0))**2 - (dlog(b1))**2 -
     .     (dlog(b2))**2 + 2d0*NS_ccspen(dcmplx(1d0-b0**2))
     .     - NS_ccspen(dcmplx(1d0-b1**2)) -
     .     NS_ccspen(dcmplx(1d0-b2**2)) ) )

      iglgl = 1d0/(4d0*m0**4)*(kap*dlog(kap**2/(lamv*m0*m1*m2))
     .     -kap-(m1**2-m2**2)*dlog(b1/b2)-m0**2*dlog(b0) )

      ist1st1 = 1d0/(4d0*m2**2*m0**2)*(
     .     kap*dlog(kap**2/(lamv*m0*m1*m2))
     .     -kap-(m0**2-m1**2)*dlog(b0/b1)-m2**2*dlog(b2))

      itt = 1d0/(4d0*m1**2*m0**2)*(kap*dlog(kap**2/(lamv*m0*m1*m2))
     .     -kap-(m0**2-m2**2)*dlog(b0/b2)-m1**2*dlog(b1) )

      igl = 1d0/(4d0*m0**2)*(-2d0*m1**2*dlog(b2)-2d0*m2**2*dlog(b1)
     .     -kap)

      ist1 = 1d0/(4d0*m0**2)*(-2d0*m0**2*dlog(b1)
     .      -2d0*m1**2*dlog(b0)-kap)

      it = 1d0/(4d0*m0**2)*(-2d0*m0**2*dlog(b2)-2d0*m2**2*dlog(b0)
     .     -kap)

      itgl = -1d0/(4d0*m0**2)*(-m0**4*dlog(b2)+m2**2*(2d0*m1**2
     .     -2d0*m0**2)*dlog(b0) + m2**4*dlog(b0) +
     .     kap/4d0*(m1**2+5d0*m2**2-3d0*m0**2) )

      iglst1 = 1d0/(4d0*m0**2)*(m2**4*dlog(b1)-m1**2*(2d0*m0**2
     .     -2d0*m2**2+m1**2)*dlog(b2) -kap/4d0*(m0**2-3d0*m2**2
     .     +5d0*m1**2) )

      NS_gamrealgl = 8d0/pi/amst1*(-(amgl**2+amt**2-amst1**2)+
     .     2d0*sig2t)*(-(amst1**2-amt**2)*ist1gl+amt**2*ist1t
     .     -amgl**2*iglgl-igl) +
     .     32d0/9d0/pi/amst1*(-(amgl**2+amt**2-amst1**2)+
     .     2d0*sig2t)*(-amst1**2*ist1st1-amt**2*itt-
     .     (amt**2+amst1**2-amgl**2)*ist1t-ist1-it) +
     .     4d0/3d0/pi/amst1*(4d0/3d0*itgl-3d0*iglst1)

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_gamcfdec(amst1,amst2,amt,amsb1,amsb2,
     .     amb,amgl,amsq,amuv,scala)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)

      mur = scala

      NS_gamcfdec = -(-dlog(mur**2/amuv**2) )*3d0 - 1d0

      NS_gamcfdec = NS_gamcfdec + 8d0/3d0

c maggie changed with respect to the paper 26/3/03
      NS_gamcfdec = NS_gamcfdec + 4d0*( 4d0/12d0*dlog(mur**2/amsq**2)
     .     + 1d0/24d0*dlog(mur**2/amsb1**2) + 1d0/24d0*
     .     dlog(mur**2/amsb2**2)
     .     + 1d0/24d0*dlog(mur**2/amst1**2) + 1d0/24d0*
     .     dlog(mur**2/amst2**2) + 1d0/6d0*dlog(mur**2/amt**2)
     .     + 1d0/2d0*dlog(mur**2/amgl**2) )
c maggie changed with respect to the paper 26/3/03

      return

      end

c -------------------------------------------------------------------- c

      subroutine NS_fifaFUNCTIONs(fffunc,fafunc,fbfunc)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      IMPLICIT INTEGER (I-J)

      DOUBLE COMPLEX SD_C03,SD_C0_lam

      DIMENSION fffunc(3),fafunc(3)
      DOUBLE PRECISION SD_B02

      COMMON/NS_qcdscales/amuv,lamv
      COMMON/NS_relmasses/amst1,amst2,amgl,amt

      pi = 4d0*datan(1d0)

      fffunc(1) = 2d0*(amt**2+amgl**2)*
     .   SD_B02(amst1**2,amgl,amt,amuv**2) +
     .   (amst1**2+amt**2+amgl**2)*SD_B02(amst1**2,lamv,amst1,amuv**2)
     .   + 2d0*(amgl**2-amst1**2)*SD_B02(amt**2,lamv,amt,amuv**2)
     .   - 2d0*amt**2*SD_B02(amt**2,amgl,amst2,amuv**2)
     .   - 4d0*amgl**2*SD_B02(amgl**2,amt,amst1,amuv**2)
     .   + 4d0*amgl**2*(amst1**2-amgl**2)*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst1) )
     .   + 2d0*amt**2*(amst1**2+amst2**2-2d0*amt**2)*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst2) )

      fffunc(2) = -2d0*SD_B02(amst1**2,lamv,amst1,amuv**2)
     .   -2d0*SD_B02(amt**2,lamv,amt,amuv**2)
     .   -4d0*SD_B02(amst1**2,amgl,amt,amuv**2)
     .   +2d0*SD_B02(amt**2,amgl,amst1,amuv**2)
     .   +4d0*SD_B02(amgl**2,amt,amst1,amuv**2)
     .   +4d0*(amgl**2+amt**2-amst1**2)*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst1) )
     .   +2d0*(amst1**2-amst2**2)*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst2) )

      fffunc(3) = 1d0/(amgl**2*amt**2)*( 2d0*amt**2*(
     .   SD_B02(amt**2,amgl,amst2,amuv**2) -
     .   SD_B02(amt**2,amgl,amst1,amuv**2) ) +
     .   (amgl**2+amt**2-amst1**2)*( (amgl**2+amt**2-amst1**2) - 4d0*
     .   amt**2 )*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst1) )
     .   -((amgl**2+amt**2-amst1**2)*(amgl**2+amt**2-amst2**2) -
     .   2d0*amt**2*(amgl**2+amt**2-amst1**2) - 2d0*amt**2*
     .   (amgl**2+amt**2-amst2**2) )*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst2) ) )

      fafunc(1) = -2d0*(amgl**2+amt**2-amst1**2) +
     .   4d0*(amt**2-amst1**2)*SD_B02(amgl**2,lamv,amgl,amuv**2)
     .   +2d0*amt**2*( SD_B02(amt**2,amgl,amst2,amuv**2) -
     .   SD_B02(amt**2,amgl,amst1,amuv**2) ) + 4d0*amgl**2*
     .   SD_B02(amgl**2,amt,amst1,amuv**2)
     .   + 4d0*amgl**2*(amgl**2-amst1**2)*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst1) )
     .   -2d0*amt**2*(amst1**2+amst2**2-2d0*amt**2)*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst2) )

      fafunc(2) = 4d0 - 4d0*SD_B02(amgl**2,lamv,amgl,amuv**2)
     .   -4d0*SD_B02(amgl**2,amt,amst1,amuv**2)
     .   -4d0*(amgl**2+amt**2-amst1**2)*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst1) )
     .   -2d0*(amst1**2-amst2**2)*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst2) )

      fafunc(3) = 1d0/(amgl**2*amt**2)*(
     .   2d0*amt**2*( SD_B02(amt**2,amgl,amst1,amuv**2)
     .   - SD_B02(amt**2,amgl,amst2,amuv**2) )
     .   + (amgl**2+amt**2-amst1**2)*(4d0*amt**2-(amgl**2+amt**2-
     .   amst1**2))*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst1) )
     .   + ( (amgl**2+amt**2-amst1**2)*(amgl**2+amt**2-amst2**2)
     .   -2d0*amt**2*(amgl**2+amt**2-amst1**2) - 2d0*amt**2*
     .   (amgl**2+amt**2-amst2**2) )*
     .   dreal( SD_C03(amst1**2,amt**2,amgl**2,amt,amgl,amst2) ))

      fbfunc = 3d0*( (amt**2+amst1**2-amgl**2)*
     .   dreal( SD_C0_lam(amst1,amt,amgl,lamv) )
     .   - (amgl**2+amt**2-amst1**2)*
     .   dreal( SD_C0_lam(amt,amgl,amst1,lamv) )
     .   - (amst1**2+amgl**2-amt**2)*
     .   dreal( SD_C0_lam(amgl,amst1,amt,lamv) ) )
     .   - 8d0/3d0*(amt**2+amst1**2-amgl**2)*
     .   dreal( SD_C0_lam(amst1,amt,amgl,lamv) )

      do i=1,3,1
         fffunc(i) = -1d0/16d0/pi**2*fffunc(i)
         fafunc(i) = -1d0/16d0/pi**2*fafunc(i)
      enddo
         fbfunc    = -1d0/16d0/pi**2*fbfunc

      end

c -------------------------------------------------------------------- c
c ---------- The A FUNCTION for the higher order corrections --------- c

      DOUBLE PRECISION FUNCTION NS_A01(s,mu2)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)

      NS_A01 = s*(-dlog(s/mu2)+1d0)

      return

      end

c -------------------------------------------------------------------- c
c ---------------------- The divergent piece of A01 ------------------ c

      DOUBLE PRECISION FUNCTION NS_A01_DIV(s,mu2)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)

      NS_A01_DIV = s*dlog(mu2)

      return

      end

c -------------------------------------------------------------------- c
c -------- The FUNCTION B1 for the higher order corrections ---------- c

      DOUBLE PRECISION FUNCTION NS_B1(s,m1,m2,mu2)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_A01,SD_B02

      NS_B1 = 1d0/2d0/s*( NS_A01(m1**2,mu2)-NS_A01(m2**2,mu2)
     .     +(m2**2-m1**2-s)*SD_B02(s,m1,m2,mu2) )

      return

      end

c -------------------------------------------------------------------- c
c ----------------------- The divergent piece of B1 ------------------ c

      DOUBLE PRECISION FUNCTION NS_B1_DIV(s,m1,m2,mu2)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION SD_B02_DIV

      NS_B1_DIV = 1d0/2d0/s*( m1**2*log(mu2)-m2**2*log(mu2)
     .     +(m2**2-m1**2-s)*SD_B02_DIV(s,m1,m2,mu2) )

      return

      end

c -------------------------------------------------------------------- c
c ------------------- The derivative of B1: dB1/ ds ------------------ c

      DOUBLE PRECISION FUNCTION SD_BP1(s,m1,m2,mu2)

      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION NS_A01,SD_B02,SD_BP02

      SD_BP1 = -1d0/2d0/s**2*( NS_A01(m1**2,mu2)-
     .     NS_A01(m2**2,mu2)+(m2**2-m1**2)*SD_B02(s,m1,m2,mu2) )
     .     +1d0/2d0/s*(m2**2-m1**2)*SD_BP02(s,m1,m2,mu2)
     .     -1d0/2d0*SD_BP02(s,m1,m2,mu2)

      return

      end

c -------------------------------------------------------------------- c
c ------- The C FUNCTION for a small mass lambda: -------------------- c
c ------- C0(msqp**2,msq**2,mphi**2,msqp,lamv,msq) =     ------------- c
c ------- C0(msq**2,mphi**2,msqp**2,lamv,msq,msqp) =     ------------- c
c - 4.2.2003 M. Muehlleitner ----------------------------------------- c
c -------------------------------------------------------------------- c

      FUNCTION SD_C0_lam(msqp,msq,mphi,lamv)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION m1,m2,lamv,msq,msqp,mphi

      DOUBLE COMPLEX SD_C0_lam,NS_ccspen,dlxs,dlfc,ieps,ys,xs

      ieps = dcmplx(0d0,1d-17)
      pi   = 4d0*datan(1d0)

      m1 = msqp
      m2 = msq
      s  = mphi**2

      ys = cdsqrt(1d0 - (4d0*m1*m2)/(s-(m1-m2)**2+ieps))

      xs = (ys-1d0)/(ys+1d0)

      dlxs = cdlog(xs)
      dlfc = cdlog(1d0-xs**2)

      SD_C0_lam = xs/(m1*m2*(1d0-xs**2))*(
     .     dlxs*(-dlog(lamv**2/(m1*m2))-1d0/2d0*dlxs+
     .     2d0*dlfc ) +
     .     NS_ccspen(1d0-xs*m1/m2) +
     .     NS_ccspen(1d0-xs*m2/m1) + NS_ccspen(xs**2) +
     .     1d0/2d0*(dlog(m1/m2))**2 - pi**2/6d0 )

      return

      end

c -------------------------------------------------------------------- c
c --      The B and C FUNCTIONs for the higher order corrections ----- c
c --      Spence FUNCTION.                                       ----- c
c -- taken from hdecay.f Version 3.0,                            ----- c
c -- authors: A.Djouadi, J.Kalinowski and M.Spira                ----- c
c -------------------------------------------------------------------- c

      FUNCTION NS_kappa(a,b,c,d)

      DOUBLE PRECISION a,b,c,d
      DOUBLE COMPLEX NS_kappa,ieps

      ieps = dcmplx(0d0,1d-17)

      if(A.eq.B) then
         NS_KAPPA = cdsqrt((C*(C-4d0*A))*(1+IEPS*D))
      elseif(B.eq.C) then
         NS_KAPPA = cdsqrt((A*(A-4d0*B))*(1+IEPS*D))
      elseif(A.eq.C) then
         NS_KAPPA = cdsqrt((B*(B-4d0*A))*(1+IEPS*D))
      else
         NS_KAPPA = CDSQRT((A**2+B**2+C**2-2*(A*B+A*C+B*C))
     .        * (1+IEPS*D))
      endif

      end

************************************************************************
      FUNCTION SD_C03(P1,P2,P3,M1,M2,M3)
************************************************************************
*  SCALAR 3-POINT FUNCTION                                             *
*  P1,P2,P3 = SQUARED EXTERNAL MOMENTA  			       *
*----------------------------------------------------------------------*
*  5.12.96  M. SPIRA    					       *
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DOUBLE PRECISION M1,M2,M3
      DOUBLE PRECISION R(0:2)
      DOUBLE COMPLEX SD_C03,NS_CCSPEN,NS_ETA,IEPS,IM
      DOUBLE COMPLEX ALP(0:2),X(0:2,2),Y0(0:2),Y(0:2,2)
      DOUBLE COMPLEX CDUM,CX,CY
C     DOUBLE PRECISION NS_KAPPA
      DOUBLE COMPLEX NS_KAPPA
c maggie changed 17/2/03
      DOUBLE COMPLEX ALPHA
      EPS = 1d-8*(P1+P2+P3)
      IM = DCMPLX(0d0,1d0)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IEPS = DCMPLX(0d0,1d-17)
c     IEPS = DCMPLX(0d0,1d-20)
      PI = 4*DATAN(1d0)
      XX = 0d0
      if(P1.NE.0d0.OR.XX.NE.0d0)then
       Q10 = P1
      else
       Q10 = EPS
      endif
      if(P3.NE.0d0.OR.XX.NE.0d0)then
       Q20 = P3
      else
       Q20 = EPS
      endif
      if(P2.NE.0d0.OR.XX.NE.0d0)then
       Q21 = P2
      else
       Q21 = EPS
      endif
      R(0) = P2
      R(1) = P3
      R(2) = P1
      SM0 = M1**2
      SM1 = M2**2
      SM2 = M3**2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ALPHA  = NS_KAPPA(Q10,Q21,Q20,1d0)
      ALP(0) = NS_KAPPA(Q21,SM1,SM2,SIGN(1d0,Q21))
      ALP(1) = NS_KAPPA(Q20,SM2,SM0,SIGN(1d0,Q20))
      ALP(2) = NS_KAPPA(Q10,SM0,SM1,SIGN(1d0,Q10))
      X(0,1) = (Q21 - SM1 + SM2 + ALP(0))/2/Q21
      X(0,2) = (Q21 - SM1 + SM2 - ALP(0))/2/Q21
      X(1,1) = (Q20 - SM2 + SM0 + ALP(1))/2/Q20
      X(1,2) = (Q20 - SM2 + SM0 - ALP(1))/2/Q20
      X(2,1) = (Q10 - SM0 + SM1 + ALP(2))/2/Q10
      X(2,2) = (Q10 - SM0 + SM1 - ALP(2))/2/Q10
      Y0(0) = (Q21*(Q21-Q20-Q10+2*SM0-SM1-SM2) - (Q20-Q10)*(SM1-SM2)
     .      + ALPHA*(Q21-SM1+SM2))/2/ALPHA/Q21
      Y0(1) = (Q20*(Q20-Q10-Q21+2*SM1-SM2-SM0) - (Q10-Q21)*(SM2-SM0)
     .      + ALPHA*(Q20-SM2+SM0))/2/ALPHA/Q20
      Y0(2) = (Q10*(Q10-Q21-Q20+2*SM2-SM0-SM1) - (Q21-Q20)*(SM0-SM1)
     .      + ALPHA*(Q10-SM0+SM1))/2/ALPHA/Q10
      Y(0,1) = Y0(0) - X(0,1)
      Y(0,2) = Y0(0) - X(0,2)
      Y(1,1) = Y0(1) - X(1,1)
      Y(1,2) = Y0(1) - X(1,2)
      Y(2,1) = Y0(2) - X(2,1)
      Y(2,2) = Y0(2) - X(2,2)
      cDUM=0d0
      do I=0,2
       do J=1,2
        cDUM = CDUM+NS_CCSPEN((Y0(I)-1d0)/Y(I,J))
     .         -NS_CCSPEN(Y0(I)/Y(I,J))
        cX = NS_ETA(1d0-X(I,J),1d0/Y(I,J))
        if(CX.NE.DCMPLX(0d0,0d0))then
         cDUM = CDUM + CX*CDLOG((Y0(I)-1)/Y(I,J))
        endif
        cY = NS_ETA(-X(I,J),1d0/Y(I,J))
        if(CY.NE.DCMPLX(0d0,0d0))then
         cDUM = CDUM - CY*CDLOG(Y0(I)/Y(I,J))
        endif
       enddo
       cX = NS_ETA(-X(I,1),-X(I,2))
       if(CX.NE.DCMPLX(0d0,0d0))then
        cDUM = CDUM - CX*CDLOG((1d0-Y0(I))/(-Y0(I)))
       endif
       cY = NS_ETA(Y(I,1),Y(I,2))
       if(CY.NE.DCMPLX(0d0,0d0))then
        cDUM = CDUM + CY*CDLOG((1d0-Y0(I))/(-Y0(I)))
       endif
       A = -R(I)
       B = -DIMAG(Y(I,1)*Y(I,2))
       if(A.GT.0d0.AND.B.GT.0d0) then
        cDUM = CDUM + 2d0*PI*IM*CDLOG((1d0-Y0(I))/(-Y0(I)))
       endif
      enddo
      SD_C03 = CDUM/ALPHA
      RETURN
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      c
C        SUBROUTINE CALCULATING THE FINITE REAL PART OF THE            c
C          GENERAL MASSIVE TWO POINT FUNCTION                          c
C                                                                      c
C           SD_B02(P.P,M1,M2,MU**2)                                    c
C           SD_BP02(P.P,M1,M2,MU**2)                                   c
C                                                                      c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION SD_B02(s,m1,m2,mu2)

      IMPLICIT NONE

      DOUBLE PRECISION s,m1,m2,mu2,m12,m22
      DOUBLE COMPLEX zkappa,x1,x2

      m12 = m1**2
      m22 = m2**2

      if(s.eq.m22) then
         zkappa = cdsqrt(dcmplx(m12*(m12-4d0*s)))
      elseif(s.eq.m12) then
         zkappa = cdsqrt(dcmplx(m22*(m22-4d0*s)))
      elseif(m12.eq.m22) then
         zkappa = cdsqrt(dcmplx(s*(s-4d0*m12)))
      else
         zkappa=cdsqrt(dcmplx(s**2+m12**2+m22**2
     .        -2d0*(s*m12+s*m22+m12*m22)))
      endif

      if (s.eq.0d0) then
         if (m12.eq.m22) then
            SD_B02=-dlog(m12/mu2)
         else
            SD_B02=1d0 - m12/(m12-m22)*dlog(m12/mu2)
     .                 + m22/(m12-m22)*dlog(m22/mu2)
         endif
      else
         if ((m12.eq.0d0).and.(m22.eq.0d0)) then
            SD_B02=2d0 - dlog(s/mu2)
         elseif ((m12.eq.s).and.(m22.eq.0d0)) then
            SD_B02=2d0 - dlog(m12/mu2)
         elseif ((m22.eq.s).and.(m12.eq.0d0)) then
            SD_B02=2d0 - dlog(m22/mu2)
         elseif (m12.eq.0d0) then
            SD_B02=2d0 - (s-m22)/s*dlog( dabs(m22-s)/m22 )
     .                 - dlog(m22/mu2)
         elseif (m22.eq.0d0) then
            SD_B02=2d0 - (s-m12)/s*dlog( dabs(m12-s)/m12 )
     .                 - dlog(m12/mu2)
         else
            x1=dcmplx( (s-m22+m12+zkappa)/(2d0*s) )
            x2=dcmplx( (s-m22+m12-zkappa)/(2d0*s) )
            SD_B02=dreal( 2d0+ dlog(mu2/m22)
     .                       + x1*cdlog(1d0-1d0/x1)
     .                       + x2*cdlog(1d0-1d0/x2))
         endif
      endif

      return
      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION SD_BP02(s,m1,m2,mu2)

      IMPLICIT NONE

      DOUBLE PRECISION s,m1,m2,mu2,m12,m22
      DOUBLE COMPLEX zkappa,x1,x2

      m12 = m1**2
      m22 = m2**2

      if(s.eq.m22) then
         zkappa = cdsqrt(dcmplx(m12*(m12-4d0*s)))
      elseif(s.eq.m12) then
         zkappa = cdsqrt(dcmplx(m22*(m22-4d0*s)))
      elseif(m12.eq.m22) then
         zkappa = cdsqrt(dcmplx(s*(s-4d0*m12)))
      else
         zkappa=cdsqrt(dcmplx(s**2+m12**2+m22**2
     .        -2d0*(s*m12+s*m22+m12*m22)))
      endif

      if (s.eq.0d0) then
         if (m12.eq.m22) then
            SD_BP02=1d0/(6d0*m12)
         else
            SD_BP02=( (m12+m22)/2d0
     .        - m12*m22/(m12-m22)*dlog(m12/m22) )/(m12-m22)**2
         endif
      elseif ((s.eq.m12).and.(m22.eq.0d0)) then
         SD_BP02=( -1d0 + dlog(m12/mu2)/2d0 )/m12
      elseif ((s.eq.m22).and.(m12.eq.0d0)) then
         SD_BP02=( -1d0 + dlog(m22/mu2)/2d0 )/m22
      elseif (m22.eq.0d0) then
         if(m12.ge.s) then
            SD_BP02=( -1d0 - m12/s*dlog((m12-s)/m12) )/s
         elseif(m12.lt.s) then
            SD_BP02=( -1d0 - m12/s*dlog((-m12+s)/m12) )/s
         endif
      else
         x1=dcmplx( (s-m22+m12+zkappa)/(2d0*s) )
         x2=dcmplx( (s-m22+m12-zkappa)/(2d0*s) )
         SD_BP02=dreal( -1d0 + ( x1*(1d0-x1)*cdlog(1d0-1d0/x1)
     .                     - x2*(1d0-x2)*cdlog(1d0-1d0/x2) )
     .                                                  /(x1-x2) )/s
      endif

      return
      end

************************************************************************
        FUNCTION NS_ETA(C1,C2)
************************************************************************
*       COMPLEX ETA-FUNKTION                                           *
*----------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       *
************************************************************************
        IMPLICIT LOGICAL(A-Z)
        DOUBLE COMPLEX NS_ETA,C1,C2
        DOUBLE PRECISION PI,IM1,IM2,IM12

        PI = 4d0*DATAN(1d0)
        IM1 = DIMAG(C1)
        IM2 = DIMAG(C2)
        IM12 = DIMAG(C1*C2)

        if(IM1.LT.0d0.AND.IM2.LT.0d0.AND.IM12.GT.0d0) then
            NS_ETA = DCMPLX(0d0,2d0*PI)
        elseif (IM1.GT.0d0.AND.IM2.GT.0d0.AND.IM12.LT.0d0) then
            NS_ETA = DCMPLX(0d0,-2d0*PI)
        else
            NS_ETA = DCMPLX(0d0)
        endif
        end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        FUNCTION NS_CCSPEN(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK                      c
C----------------------------------------------------------------------C
C       20.07.83    LAST CHANGED 10.05.89        ANSGAR DENNER         c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        IMPLICIT INTEGER (I-N)
        DOUBLE COMPLEX NS_CCSPEN,W,SUM,Z,U
        DOUBLE PRECISION RZ,AZ,A1
        DOUBLE PRECISION B(9)/
     1   0.1666666666666666666666666667d0,
     2  -0.0333333333333333333333333333d0,
     3   0.0238095238095238095238095238d0,
     4  -0.0333333333333333333333333333d0,
     5   0.0757575757575757575757575758d0,
     6  -0.2531135531135531135531135531d0,
     7   1.1666666666666666666666666667d0,
     8  -7.09215686274509804d0         ,
     9  54.97117794486215539d0         /
C     BEACHTE:                 B(N)=B2N
C     B(1)=1./6.
C     B(2)=-1./30.
C     B(3)=1./42.
C     B(4)=-1./30.
C     B(5)=5./66.
C     B(6)=-691./2730.
C     B(7)=7./6.
C     B(8)=-3617./510.
C     B(9)=43867./798.
C     B(10)=-174611./330.
C     B(11)=854513./138.
C     PI=3.1415926535897932384
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...
C
      Z =Z*DCMPLX(1d0)
      RZ=DREAL(Z)
      AZ=CDABS(Z)
      A1=CDABS(1d0-Z)
C     if((SNGL(RZ) .EQ. 0d0) .AND. (SNGL(DIMAG(Z)) .EQ. 0d0)) then
C ---> CHANGED  10.5.89
      if(AZ .LT. 1d-20) then
        NS_CCSPEN=-CDLOG(1d0-Z)
        RETURN
      end if
      if((SNGL(RZ) .EQ. 1d0) .AND. (SNGL(DIMAG(Z)) .EQ. 0d0)) then
        NS_CCSPEN=1.64493406684822643d0
        RETURN
      end if
      if(RZ.GT.5d-1) GOTO 20
      if(AZ.GT.1d0) GOTO 10
      W=-CDLOG(1d0-Z)
      SUM=W-0.25d0*W*W
      U=W
      if(CDABS(U).LT.1d-10) GOTO 2
      do 1 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      if(CDABS(U*B(K)/SUM).LT.1d-20) GOTO 2
      SUM=SUM+U*B(K)
 1    CONTINUE
 2    NS_CCSPEN=SUM
      RETURN
10    W=-CDLOG(1d0-1d0/Z)
      SUM=W-0.25d0*W*W
      U=W
      if(CDABS(U).LT.1d-10) GOTO 12

      do 11 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      if(CDABS(B(K)*U/SUM).LT.1d-20) GOTO 12
      SUM=SUM+U*B(K)
11    CONTINUE
12    NS_CCSPEN=-SUM-1.64493406684822643d0-.5d0*CDLOG(-Z)**2
      RETURN
20    if(A1.GT.1d0) GOTO 30
      W=-CDLOG(Z)
      SUM=W-0.25d0*W*W
      U=W
      if(CDABS(U).LT.1d-10) GOTO 22
      do 21 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      if(CDABS(U*B(K)/SUM).LT.1d-20) GOTO 22
      SUM=SUM+U*B(K)
21    CONTINUE
22    NS_CCSPEN=-SUM+1.64493406684822643d0-CDLOG(Z)*CDLOG(1d0-Z)
      RETURN
30    W=CDLOG(1d0-1d0/Z)
      SUM=W-0.25d0*W*W
      U=W
      if(CDABS(U).LT.1d-10) GOTO 32
      do 31 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      if(CDABS(U*B(K)/SUM).LT.1d-20) GOTO 32
      SUM=SUM+U*B(K)
31    CONTINUE
32    NS_CCSPEN=SUM+3.28986813369645287d0
     *               +.5d0*CDLOG(Z-1d0)**2-CDLOG(Z)*CDLOG(1d0-Z)
      end

c ==================================================================== c
c                     The C11 and C12 FUNCTIONs                        c
c                    c_mu = p1_mu*C11 + p2_mu*C12                      c
c  11.2.03 M.Muehlleitner                                              c
c  p1,p2,p3 squared EXTERNAL momenta                                   c
c ==================================================================== c

      DOUBLE PRECISION FUNCTION NS_C1(p1,p2,p3,m1,m2,m3,amuv)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION m1,m2,m3
      DOUBLE COMPLEX SD_C03

      f1 = m2**2-m1**2-p1**2
      f2 = m3**2-m2**2-p3**2+p1**2

      den = p1**2*p2**2- (p3**2-p1**2-p2**2)**2/4d0

      r1 = 1d0/2d0*(
     .     f1*dreal(SD_C03(p1**2,p2**2,p3**2,m1,m2,m3)) +
     .     SD_B02(p3**2,m1,m3,amuv**2) -
     .     SD_B02(p2**2,m2,m3,amuv**2) )

      r2 = 1d0/2d0*(
     .     f2*dreal(SD_C03(p1**2,p2**2,p3**2,m1,m2,m3)) +
     .     SD_B02(p1**2,m1,m2,amuv**2) -
     .     SD_B02(p3**2,m1,m3,amuv**2) )

      NS_C1 = 1d0/den*( p2**2*r1 - (p3**2-p1**2-p2**2)/2d0*r2 )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION SD_C2(p1,p2,p3,m1,m2,m3,amuv)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION SD_B02
      DOUBLE PRECISION m1,m2,m3
      DOUBLE COMPLEX SD_C03

      f1 = m2**2-m1**2-p1**2
      f2 = m3**2-m2**2-p3**2+p1**2

      den = p1**2*p2**2 - (p3**2-p1**2-p2**2)**2/4d0

      r1 = 1d0/2d0*(
     .     f1*dreal(SD_C03(p1**2,p2**2,p3**2,m1,m2,m3)) +
     .     SD_B02(p3**2,m1,m3,amuv**2) -
     .     SD_B02(p2**2,m2,m3,amuv**2) )

      r2 = 1d0/2d0*(
     .     f2*dreal(SD_C03(p1**2,p2**2,p3**2,m1,m2,m3)) +
     .     SD_B02(p1**2,m1,m2,amuv**2) -
     .     SD_B02(p3**2,m1,m3,amuv**2) )

      SD_C2 = 1d0/den*( -(p3**2-p1**2-p2**2)/2d0*r1 + p1**2*r2 )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION NS_C1_lam(p1,p2,p3,m1,m2,m3,amuv,lamv)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION SD_B02
      DOUBLE PRECISION m1,m2,m3,lamv
      DOUBLE COMPLEX SD_C0_lam

      f1 = m2**2-m1**2-p1**2
      f2 = m3**2-m2**2-p3**2+p1**2

      den = p1**2*p2**2- (p3**2-p1**2-p2**2)**2/4d0

      r1 = 1d0/2d0*( f1*dreal(SD_C0_lam(p1,p2,p3,lamv)) +
     .     SD_B02(p3**2,m1,m3,amuv**2) -
     .     SD_B02(p2**2,m2,m3,amuv**2) )

      r2 = 1d0/2d0*( f2*dreal(SD_C0_lam(p1,p2,p3,lamv)) +
     .     SD_B02(p1**2,m1,m2,amuv**2) -
     .     SD_B02(p3**2,m1,m3,amuv**2) )

      NS_C1_lam = 1d0/den*( p2**2*r1 - (p3**2-p1**2-p2**2)/2d0*r2 )

      return

      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION SD_C2_lam(p1,p2,p3,m1,m2,m3,amuv,lamv)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION SD_B02
      DOUBLE PRECISION m1,m2,m3,lamv
      DOUBLE COMPLEX SD_C0_lam

      f1 = m2**2-m1**2-p1**2
      f2 = m3**2-m2**2-p3**2+p1**2

      den = p1**2*p2**2 - (p3**2-p1**2-p2**2)**2/4d0

      r1 = 1d0/2d0*( f1*dreal(SD_C0_lam(p1,p2,p3,lamv)) +
     .     SD_B02(p3**2,m1,m3,amuv**2) -
     .     SD_B02(p2**2,m2,m3,amuv**2) )

      r2 = 1d0/2d0*( f2*dreal(SD_C0_lam(p1,p2,p3,lamv)) +
     .     SD_B02(p1**2,m1,m2,amuv**2) -
     .     SD_B02(p3**2,m1,m3,amuv**2) )

      SD_C2_lam = 1d0/den*( -(p3**2-p1**2-p2**2)/2d0*r1 + p1**2*r2 )

      return

      end

c -------------------------------------------------------------------- c
c            The divergent pieces of the B FUNCTIONs                   c
c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION SD_B02_DIV(s,m1,m2,mu2)

      IMPLICIT NONE

      DOUBLE PRECISION s,m1,m2,mu2,m12,m22

      m12 = m1**2
      m22 = m2**2

      if (s.eq.0d0) then
         if (m12.eq.m22) then
            SD_B02_DIV=dlog(mu2)
         else
            SD_B02_DIV= + m12/(m12-m22)*dlog(mu2)
     .               - m22/(m12-m22)*dlog(mu2)
         endif
      else
         if ((m12.eq.0d0).and.(m22.eq.0d0)) then
            SD_B02_DIV= dlog(mu2)
         elseif ((m12.eq.s).and.(m22.eq.0d0)) then
            SD_B02_DIV= dlog(mu2)
         elseif ((m22.eq.s).and.(m12.eq.0d0)) then
            SD_B02_DIV= dlog(mu2)
         elseif (m12.eq.0d0) then
            SD_B02_DIV= dlog(mu2)
         elseif (m22.eq.0d0) then
            SD_B02_DIV= dlog(mu2)
         else
            SD_B02_DIV= dlog(mu2)
         endif
      endif

      return
      end

c -------------------------------------------------------------------- c

      DOUBLE PRECISION FUNCTION SD_BP02_DIV(s,m1,m2,mu2)

      IMPLICIT NONE

      DOUBLE PRECISION s,m1,m2,mu2,m12,m22

      m12 = m1**2
      m22 = m2**2

      if (s.eq.0d0) then
         if (m12.eq.m22) then
            SD_BP02_DIV=0d0
         else
            SD_BP02_DIV=0d0
         endif
      elseif ((s.eq.m12).and.(m22.eq.0d0)) then
         SD_BP02_DIV=( - dlog(mu2)/2d0 )/m12
      elseif ((s.eq.m22).and.(m12.eq.0d0)) then
         SD_BP02_DIV=( - dlog(mu2)/2d0 )/m22
      else
         SD_BP02_DIV= 0d0
      endif

      return
      end

c -------------------- Derivatives of couplings ---------------------- c
c -------------------------------------------------------------------- c
c                   H+ - stop1/2 - sbottom1/2 couplings                c
c -------------------------------------------------------------------- c
      subroutine NS_hcsbotstopderiv(gcdmtr,gcdmbr,gcdabr,gcdatr,
     .                      gcdthtr,gcdthbr)
      IMPLICIT NONE
*
      DOUBLE PRECISION chctbdt(2,2),chctbdb(2,2),chctbab(2,2),
     .     chctbat(2,2),
     .     chctbtt(2,2),gcdthtr(2,2),gcdthbr(2,2),chctbbb(2,2),
     .     gcdmtr(2,2),gcdmbr(2,2),gcdabr(2,2),gcdatr(2,2)
      DOUBLE PRECISION thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .cmm,smm,cum,sum,cdm,sdm,cem,sem,cnm,snm
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION au,ad,al,mu
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION b,tgbet
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION tmpt,s11t,s12t,s21t,s22t,tmpb,s11b,s12b,
     . s21b,s22b,s11ab,s12ab,s21ab,s22ab,s11at,s12at,s21at,s22at,
     .ctt,stt,cbb,sbb,s11,s12,s21,s22
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/NS_sfmixang/thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .cmm,smm,cum,sum,cdm,sdm,cem,sem,cnm,snm
      COMMON/NS_trilin_mu/au,ad,al,mu
      COMMON/NS_weinberg/sw,cw,tw
      COMMON/NS_tanb/tanbeta
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
**
      b     = datan(tanbeta)
      tgbet = tanbeta

c ---- derivative d/dmtop ----

      tmpt = 1d0/dsqrt(2d0)/amw/DDSIN(b)

      s11t = 1d0/dsqrt(2d0)*amw*(2d0*tmpt*scalt)*DDSIN(2d0*b)
      s12t = 0d0
      s21t = tmpt*DDSIN(b)*(au/tgbet+mu)
      s22t = dsqrt(2d0)*amw*scalb*tmpt

      chctbdt(1,1)=(-ct*cb*s11t-st*sb*s22t-ct*sb*s12t-st*cb*s21t)
      chctbdt(1,2)=(ct*sb*s11t-ct*cb*s12t+sb*st*s21t-st*cb*s22t)
      chctbdt(2,1)=(st*cb*s11t+st*sb*s12t-ct*cb*s21t-ct*sb*s22t)
      chctbdt(2,2)=(-st*sb*s11t+st*cb*s12t+ct*sb*s21t-ct*cb*s22t)

      gcdmtr(1,1)=chctbdt(1,1)/amw
      gcdmtr(1,2)=chctbdt(1,2)/amw
      gcdmtr(2,1)=chctbdt(2,1)/amw
      gcdmtr(2,2)=chctbdt(2,2)/amw

c ---- derivative d/dmbottom ----

      tmpb = 1d0/dsqrt(2d0)/amw/DDCOS(b)

      s11b = 1d0/dsqrt(2d0)*amw*(2d0*tmpb*scalb)*DDSIN(2d0*b)
      s12b = tmpb*DDCOS(b)*(ad*tgbet+mu)
      s21b = 0d0
      s22b = dsqrt(2d0)*amw*tmpb*scalt

      chctbdb(1,1)=(-ct*cb*s11b-st*sb*s22b-ct*sb*s12b-st*cb*s21b)
      chctbdb(1,2)=(ct*sb*s11b-ct*cb*s12b+sb*st*s21b-st*cb*s22b)
      chctbdb(2,1)=(st*cb*s11b+st*sb*s12b-ct*cb*s21b-ct*sb*s22b)
      chctbdb(2,2)=(-st*sb*s11b+st*cb*s12b+ct*sb*s21b-ct*cb*s22b)

      gcdmbr(1,1)=chctbdb(1,1)/amw
      gcdmbr(1,2)=chctbdb(1,2)/amw
      gcdmbr(2,1)=chctbdb(2,1)/amw
      gcdmbr(2,2)=chctbdb(2,2)/amw

c ---- derivative d/dAb ----

      s11ab = 0d0
      s12ab = scalb*DDCOS(b)*tgbet
      s21ab = 0d0
      s22ab = 0d0

      chctbab(1,1)=(-ct*cb*s11ab-st*sb*s22ab-ct*sb*s12ab-st*cb*s21ab)
      chctbab(1,2)=(ct*sb*s11ab-ct*cb*s12ab+sb*st*s21ab-st*cb*s22ab)
      chctbab(2,1)=(st*cb*s11ab+st*sb*s12ab-ct*cb*s21ab-ct*sb*s22ab)
      chctbab(2,2)=(-st*sb*s11ab+st*cb*s12ab+ct*sb*s21ab-ct*cb*s22ab)

      gcdabr(1,1)=chctbab(1,1)/amw
      gcdabr(1,2)=chctbab(1,2)/amw
      gcdabr(2,1)=chctbab(2,1)/amw
      gcdabr(2,2)=chctbab(2,2)/amw

c ---- derivative d/dAt ----

      s11at = 0d0
      s12at = 0d0
      s21at = scalt*DDSIN(b)*1d0/tgbet
      s22at = 0d0

      chctbat(1,1)=(-ct*cb*s11at-st*sb*s22at-ct*sb*s12at-st*cb*s21at)
      chctbat(1,2)=(ct*sb*s11at-ct*cb*s12at+sb*st*s21at-st*cb*s22at)
      chctbat(2,1)=(st*cb*s11at+st*sb*s12at-ct*cb*s21at-ct*sb*s22at)
      chctbat(2,2)=(-st*sb*s11at+st*cb*s12at+ct*sb*s21at-ct*cb*s22at)

      gcdatr(1,1)=chctbat(1,1)/amw
      gcdatr(1,2)=chctbat(1,2)/amw
      gcdatr(2,1)=chctbat(2,1)/amw
      gcdatr(2,2)=chctbat(2,2)/amw

c ---- derivative d/dtheta_t ----

      ctt = -DDSIN(thet)
      stt = DDCOS(thet)

      s11 = 1d0/dsqrt(2d0)*amw*(scalb**2+scalt**2-1d0)*DDSIN(2d0*b)
      s12 = scalb*DDCOS(b)*(ad*tgbet+mu)
      s21 = scalt*DDSIN(b)*(au/tgbet+mu)
      s22 = dsqrt(2d0)*amw*scalb*scalt

      chctbtt(1,1)=(-ctt*cb*s11-stt*sb*s22-ctt*sb*s12-stt*cb*s21)
      chctbtt(1,2)=(ctt*sb*s11-ctt*cb*s12+sb*stt*s21-stt*cb*s22)
      chctbtt(2,1)=(stt*cb*s11+stt*sb*s12-ctt*cb*s21-ctt*sb*s22)
      chctbtt(2,2)=(-stt*sb*s11+stt*cb*s12+ctt*sb*s21-ctt*cb*s22)

      gcdthtr(1,1)=chctbtt(1,1)/amw
      gcdthtr(1,2)=chctbtt(1,2)/amw
      gcdthtr(2,1)=chctbtt(2,1)/amw
      gcdthtr(2,2)=chctbtt(2,2)/amw

c ---- derivative d/dtheta_b ----

      cbb = -DDSIN(theb)
      sbb = DDCOS(theb)

      s11 = 1d0/dsqrt(2d0)*amw*(scalb**2+scalt**2-1d0)*DDSIN(2d0*b)
      s12 = scalb*DDCOS(b)*(ad*tgbet+mu)
      s21 = scalt*DDSIN(b)*(au/tgbet+mu)
      s22 = dsqrt(2d0)*amw*scalb*scalt

      chctbbb(1,1)=(-ct*cbb*s11-st*sbb*s22-ct*sbb*s12-st*cbb*s21)
      chctbbb(1,2)=(ct*sbb*s11-ct*cbb*s12+sbb*st*s21-st*cbb*s22)
      chctbbb(2,1)=(st*cbb*s11+st*sbb*s12-ct*cbb*s21-ct*sbb*s22)
      chctbbb(2,2)=(-st*sbb*s11+st*cbb*s12+ct*sbb*s21-ct*cbb*s22)

      gcdthbr(1,1)=chctbbb(1,1)/amw
      gcdthbr(1,2)=chctbbb(1,2)/amw
      gcdthbr(2,1)=chctbbb(2,1)/amw
      gcdthbr(2,2)=chctbbb(2,2)/amw

      end

c -------------------------------------------------------------------- c
c --------------------------- The counterterms ----------------------- c

      DOUBLE PRECISION FUNCTION NS_dcounterhc(amsq,amq,theq,ni,amsqp,
     .                   amqp,theqp,nj,mgluino,amuv,amuvdiv,lamv,ic,jc)

      IMPLICIT NONE
      INTEGER ic,jc,ni,nj
      DOUBLE PRECISION lamv
      DOUBLE PRECISION gctbr(2,2),gcdthtr(2,2),gcdthbr(2,2),gcdmtr(2,2),
     .          gcdmbr(2,2),gcdabr(2,2),gcdatr(2,2)
      DOUBLE PRECISION AMZ,AMW
      DOUBLE PRECISION mgluino
      DOUBLE PRECISION asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      DOUBLE PRECISION thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      DOUBLE PRECISION scalb,scalt,scaltau,gs2
      DOUBLE PRECISION NS_deltaz,NS_deltamqdiv,NS_deltaAq
     .,NS_deltathdiv
      DOUBLE PRECISION runmt,runmb,rmtauc
      DOUBLE PRECISION amtt,ambb,amsq,amq,theq,amuv,amsqp,amqp,
     .theqp,amuvdiv
*
      COMMON/NS_sfmixang/thet,theb,thel,them,ct,st,cb,sb,cl,sl,
     .     cm,sm,cu,su,cd,sd,ce,se,cn,sn
      COMMON/NS_scala/scalb,scalt,scaltau,gs2
      COMMON/SFSPEC/asup2,asup1,asdown2,asdown1,ase2,ase1,asne1,
     .     ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .     CST,CSB,CSL,asmu1,asmu2,asnmu1
      COMMON/NS_MZMWscaleQ/AMZ,AMW
      COMMON/NS_runmcalc/runmt,runmb,rmtauc
      COMMON/NS_hcsbotstop/gctbr
*
      CALL NS_hcsbotstopderiv(gcdmtr,gcdmbr,gcdabr,gcdatr,gcdthtr,
     .                        gcdthbr)
c --- the running mass ---
      amtt = runmt
      ambb = runmb
*
      NS_dcounterhc = 1d0/2d0*gctbr(ic,jc)*(
     .     NS_deltaz(amsq,mgluino,amq,theq,amuv,lamv,ni) +
     .     NS_deltaz(amsqp,mgluino,amqp,theqp,amuv,lamv,nj) ) +
     .     gcdmtr(ic,jc)*amtt*
     .     NS_deltamqdiv(ast1,ast2,mgluino,amtt,thet,amuvdiv,lamv) +
     .     gcdmbr(ic,jc)*ambb*
     .     NS_deltamqdiv(asb1,asb2,mgluino,ambb,theb,amuvdiv,lamv) +
     .     gcdatr(ic,jc)*
     .     NS_deltaAq(amsq,ast1,ast2,mgluino,amtt,thet,amuvdiv,lamv,1) +
     .     gcdabr(ic,jc)*
     .     NS_deltaAq(amsq,asb1,asb2,mgluino,ambb,theb,amuvdiv,lamv,2) +
     .     gcdthtr(ic,jc)*NS_deltathdiv(ast1,ast2,mgluino,amtt,thet,
     .                               amuvdiv) +
     .     gcdthbr(ic,jc)*NS_deltathdiv(asb1,asb2,mgluino,ambb,theb,
     .                               amuvdiv)

      NS_dcounterhc = -dsqrt(2d0)*amw**2*NS_dcounterhc

      return

      end
