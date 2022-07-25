* This file contains
*
* the main routine FULLHIG for the full 1- and 2-Loop corrections
*   from Degrassi/Slavich,
*   ``On the radiative corrections to the neutral Higgs boson masses in
*    the NMSSM,'' Nucl. Phys.  B {\bf 825} (2010) 119
*   [arXiv:0907.4682 [hep-ph]].
*
* SUBROUTINE mysort(msd,ZS)
* SUBROUTINE getdVB(g,gp,mzpole,mwpole,Q,dVB)
* SUBROUTINE effpot(lp,mt,mg,T1,T2,st,ct,q2,tanb,vv,l,xx,as,
*     .     DMS,DMP)
* SUBROUTINE makefuncs(mt,mg,T1,T2,s2t,c2t,q,tanb,At,mu,
*     .     F1t,F2t,F3t,Ft,FA)
* SUBROUTINE makederiv(mt,mg,T1,T2,s2t,c2t,q,
*     .     DT1,DT1T1,DT1t,DT1c2t,DT1T2,Dtt,Dc2t,Dc2tc2t,Dtc2t,Dcptmptt)
* SUBROUTINE diagonalize(mq,mw,mz,msql,msqr,Aq,mu,tb,iq,
*     .     msq2,sth,cth)
* SUBROUTINE squarks(mtsm,mbsm,mz,mw,mqs,mts,mbs,At,Ab,mu,tb,mg,as,
*     .     A0,vv,q2,mtmssm,mbmssm,mstop2,msbot2,cst,sst,csb,ssb,
*     .     errsqua,mtpole,asma,ast,v2)
* SUBROUTINE getdmt(mt,mstop2,sst,cst,mb,msbot2,ssb,csb,
*     .      vv,mg,A0,mu,tb,as,q2,dmt)
* SUBROUTINE getdmb(mt,mstop2,sst,cst,mb,msbot2,ssb,csb,
*     .      vv,mg,A0,mu,tb,as,q2,dd,eb)
* SUBROUTINE swap12(M)
* SUBROUTINE jacobi(a,np,d,v)
* SUBROUTINE HEigensystem(n, A,ldA, d, U,ldU, sort)
* SUBROUTINE gettadS(g,gp,ht,hb,htau,v1,v2,Q,tadS)
* SUBROUTINE getPiSS(g,gp,ht,hb,htau,v1,v2,p,Q,piSS)
* SUBROUTINE getPiPP(g,gp,ht,hb,htau,v1,v2,p,Q,piPP)
* SUBROUTINE getPiZZ(g,gp,ht,hb,htau,v1,v2,p,Q,piZZ)
* SUBROUTINE getPiWW(g,gp,ht,hb,htau,v1,v2,p,Q,piWW)
* SUBROUTINE treemasses(g,gp,ll,kk,ht,hb,htau,v1,v2,xx,M1,M2,
*     .     Ak,Al,At,Ab,Atau,mQ3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,mer,
*     .     Q,errmass)
*     defining the
*   COMMON/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
*     .     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue
* SUBROUTINE tree_charginos(g,ll,v1,v2,M2,xx,xmc,u,v)
* SUBROUTINE tree_neutralinos(g,gp,ll,kk,v1,v2,xx,M1,M2,xmn,Z)
* SUBROUTINE tree_higgses(g,gp,ll,kk,v1,v2,xx,Ak,Al,
*     .     mss,maa,mhc,RS,RP,RC,errhiggs)
* SUBROUTINE tree_sfermions(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
*     .     mQ3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,mer,mstop,msbot,
*     .     mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue,errsfer)
* SUBROUTINE diagsfe(n,hf,g,gp,ll,v1,v2,xx,Af,mL,mR,mass,R,error)
* SUBROUTINE scalarcouplings(g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
*     .    Al,Ak,At,Ab,Atau)
* SUBROUTINE coupl_s_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
*     .     Rt,Rb,Rtau,lsstt,lssbb,lsstata,lssntnt,lssuu,lssdd,
*     .     lssee,lssnn,lstt,lsbb,lstata,lsntnt,lsuu,lsdd,lsee,lsnn)
* SUBROUTINE coupl_s_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
*     .     lsshh,lssaa,lshh,lsaa,lsscc,lscc)
* SUBROUTINE coupl_s_ino(g,gp,ll,kk,NN,UU,VV,lsnene,lschch)
* SUBROUTINE coupl_p_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
*     .     Rt,Rb,Rtau,lpptt,lppbb,lpptata,lppntnt,lppuu,lppdd,
*     .     lppee,lppnn,lptt,lpbb,lptata,lpntnt,lpuu,lpdd,lpee,lpnn)
* SUBROUTINE coupl_p_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
*     .     lpphh,lppaa,lpah,lppcc,lpcc)
* SUBROUTINE coupl_p_ino(g,gp,ll,kk,NN,UU,VV,lpnene,lpchch)
* SUBROUTINE coupl_Z_ino(g,gp,NN,UU,VV,lznene,azchch,bzchch)
* SUBROUTINE coupl_W_ino(g,NN,UU,VV,awnech,bwnech)
* SUBROUTINE twlpyuk(mt,mb,A0,T1,T2,B1,B2,st,ct,sb,cb,q,l,xx,tanb,
*     .     vv,DMS,DMP)
* SUBROUTINE makefuncstb(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,
*     .     q,mu,vv,tanb,F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,
*     .     Ft,Fb,FA)
* SUBROUTINE makederivtb(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,
*     .     q,mu,vv,tanb)

      SUBROUTINE FULLHIG(mom,l,k,mu,tanb,mq3,mtr,mbr,
     .     mQ,mur,mdr,mL3,mtaur,mL,mer,At,Ab,Atau,Al,Ak,M1,M2,mg,Q,
     .     mss,maa,OS,OP,mhc,err)

      IMPLICIT NONE

      LOGICAL mom,errsqua
      INTEGER i,j,kk
      INTEGER err,errmass,IL
      DOUBLE PRECISION l,k,mu,tanb,mq3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,
     .     mer,At,Ab,Atau,Al,Ak,M1,M2,mg,q,mss(3),maa(3),OS(3,3),
     .     OP(3,3),mhc
      DOUBLE PRECISION g,gp,gb2,vv,cb,sb,v1,v2,slasf,asq,mtsm,mbsm,
     .     runt,runb,pi,mtmssm,mbmssm,mstop2(2),msbot2(2),cst,sst,
     .     csb,ssb,As,MStr(3,3),MPtr(3,3),DMS_2lt(3,3),
     .     DMP_2lt(3,3),DMS_2lb(3,3),DMP_2lb(3,3),MS(3,3),MP(3,3),
     .     ms2(3),ma2(3),xx,tadS(3),piSS(3,3),piPP(3,3),
     .     gold,ms2cor(3),ma2cor(3),
     .     MSfull(3,3),MPfull(3,3),pa,ps,vev(3),mz,mw,
     .     piZZ_MZ,piWW_MW,piWW_0,dVB,sq2,T1,T2,B1,B2,
     .     A0,DMS_2ly(3,3),DMP_2ly(3,3)
      DOUBLE PRECISION OS0(3,3),OP0(3,3),asma,ast
      DOUBLE PRECISION mZpole,mWpole,mtpole,mbmb,mtau,asmz,GF
* For charged Higgs:
      DOUBLE PRECISION mhcsqtree, mhcsqcorr,phc,piHpHm

      COMMON/SMINPUTS/mZpole,mWpole,mtpole,mbmb,mtau,asmz,GF

*      write(*,*) 'inputs_1',mom,l,k,mu,tanb,mq3,mtr,mbr,
*     .     mQ,mur,mdr,mL3,mtaur,mL,mer,At,Ab,Atau,Al,Ak,M1,M2,mg,Q
*      write(*,*) 'inputs_2',mZpole,mWpole,mtpole,mbmb,mtau,asmz,GF

      err=0

*     some preliminary quantities

      pi = 4d0*atan(1d0)
      sq2 = sqrt(2d0)
      cb = 1d0/sqrt(1d0+tanb**2)
      sb = tanb*cb

      asq = slasf(q,mtpole,asmz,mzpole)

*     determine the running couplings and quark/squark masses

      vv  = 1d0/sqrt(2d0*sq2*GF)    ! note: v ~ 174
      g = sqrt(2d0*mwpole**2/vv**2)
      gp = sqrt(2d0*(mzpole**2-mwpole**2)/vv**2)
      mz = mzpole
      mw = mwpole

      v1 = vv*cb
      v2 = vv*sb
      xx = mu/l

      A0 = l*xx*(Al+k*xx)/sb/cb      ! the would-be mA^2 in the MSSM limit

      mtsm = runt(q,mtpole,asmz,mzpole,vv,mbmb) ! SM, DRbar running masses at q
      mbsm = runb(q,mbmb,mtpole,asmz,mzpole)

* mod. by UE:
      IF(A0.le.1d4) A0=1d4
* asma=Alpha_s at the scale A0:
      asma = slasf(dsqrt(A0),mtpole,asmz,mzpole)
* ast=Alpha_s at the scale mtsm:
      ast = slasf(mtsm,mtpole,asmz,mzpole)
* (Used in SUBROUTINE squarks)
* End mod. by UE      

      do i=1,5

      CALL squarks(mtsm,mbsm,mz,mw,mq3,mtr,mbr,At,Ab,l*xx,tanb,mg,
     .      asq,A0,vv,q**2,mtmssm,mbmssm,mstop2,msbot2,cst,sst,
     .      csb,ssb,errsqua,mtpole,asma,ast,v2)

       IF(errsqua) THEN
         err=3
         RETURN
       ENDIF

      CALL treemasses(g,gp,l,k,mtmssm/v2,mbmssm/v1,mtau/v1,v1,v2,xx,
     .      M1,M2,Ak,Al,At,Ab,Atau,mQ3,mtr,mbr,mQ,mur,mdr,mL3,
     .      mtaur,mL,mer,Q,errmass)

       IF(errmass.ne.0) THEN
         err=errmass
         RETURN
       ENDIF
 
      CALL getPiZZ(g,gp,mtmssm/v2,mbmssm/v1,mtau/v1,
     .      v1,v2,mzpole,Q,piZZ_MZ)

      CALL getPiWW(g,gp,mtmssm/v2,mbmssm/v1,mtau/v1,
     .      v1,v2,mwpole,Q,piWW_MW)

      CALL getPiWW(g,gp,mtmssm/v2,mbmssm/v1,mtau/v1,
     .      v1,v2,0d0,Q,piWW_0)

      CALL getdVB(g,gp,mzpole,mwpole,Q,dVB)

       vv = 1d0/sqrt(2d0*sq2*GF)/sqrt(1d0-piWW_0/mwpole**2-dVB)
       v1 = vv*cb
       v2 = vv*sb

*UE,9.9.2014:
      IF(piWW_MW/mwpole**2.GE.0.0d0) THEN
       g = sqrt(2d0*mwpole**2/vv**2*(1d0+piWW_MW/mwpole**2))
      ELSE
       g = sqrt(2d0*mwpole**2/vv**2/(1d0-piWW_MW/mwpole**2))
      ENDIF
      IF(piZZ_MZ/mzpole**2.GE.0.0d0) THEN
       gp = sqrt(2d0*mzpole**2/vv**2*(1d0+piZZ_MZ/mzpole**2)-g**2)
      ELSE
       gp = sqrt(2d0*mzpole**2/vv**2/(1d0-piZZ_MZ/mzpole**2)-g**2)
      ENDIF

       mz = sqrt(g**2+gp**2)/sq2*vv
       mw = g/sq2*vv

      enddo

*     tree-level mass matrices

      As = Al+k*xx            ! shortcuts
      gb2 = (g**2+gp**2)/2d0

      MStr(1,1) = gb2*v1**2+l*xx*v2/v1*As
      MStr(1,2) = (2d0*l**2-gb2)*v1*v2-l*xx*As
      MStr(1,3) = 2d0*l**2*v1*xx-l*v2*(As+k*xx)
      MStr(2,2) = gb2*v2**2+l*xx*v1/v2*As
      MStr(2,3) = 2d0*l**2*v2*xx-l*v1*(As+k*xx)
      MStr(3,3) = l*Al*v1*v2/xx+k*xx*(Ak+4d0*k*xx)
      MStr(2,1) = MStr(1,2)
      MStr(3,1) = MStr(1,3)
      MStr(3,2) = MStr(2,3)

      gold = mz**2        ! gauge-fixing mass

      MPtr(1,1) = l*xx*v2/v1*As+cb**2*gold
      MPtr(1,2) = l*xx*As-cb*sb*gold
      MPtr(1,3) = l*v2*(As-3d0*k*xx)
      MPtr(2,2) = l*xx*v1/v2*As+sb**2*gold
      MPtr(2,3) = l*v1*(As-3d0*k*xx)
      MPtr(3,3) = 4d0*l*k*v1*v2+l*Al*v1*v2/xx-3d0*k*Ak*xx
      MPtr(2,1) = MPtr(1,2)
      MPtr(3,1) = MPtr(1,3)
      MPtr(3,2) = MPtr(2,3)

*     compute the top corrections at zero momentum

      T1 = mstop2(1)
      T2 = mstop2(2)


      CALL effpot(2,mtmssm,mg,T1,T2,sst,cst,q**2,
     .     tanb,vv,l,xx,asq,DMS_2lt,DMP_2lt)

*     compute the bottom corrections at zero momentum

      B1 = msbot2(1)
      B2 = msbot2(2)


      CALL effpot(2,mbmssm,mg,B1,B2,ssb,csb,q**2,
     .     1d0/tanb,vv,l,xx,asq,DMS_2lb,DMP_2lb)


      CALL swap12(DMS_2lb)

      CALL swap12(DMP_2lb)

*     compute the two-loop top-bot Yukawa corrections

 
      CALL twlpyuk(mtmssm,mbmssm,A0,T1,T2,B1,B2,sst,cst,ssb,csb,q**2,
     .     l,xx,tanb,vv,DMS_2ly,DMP_2ly)

*     compute the full one-loop at zero momentum


      CALL treemasses(g,gp,l,k,mtmssm/v2,mbmssm/v1,mtau/v1,v1,
     .     v2,xx,M1,M2,Ak,Al,At,Ab,Atau,mQ3,mtr,mbr,mQ,mur,mdr,mL3,
     .     mtaur,mL,mer,Q,errmass)

       IF(errmass.ne.0) THEN
         err=errmass
         RETURN
       ENDIF


      CALL scalarcouplings(g,gp,l,k,mtmssm/v2,mbmssm/v1,
     .     mtau/v1,v1,v2,xx,Al,Ak,At,Ab,Atau)


      CALL gettadS(g,gp,mtmssm/v2,mbmssm/v1,
     .     mtau/v1,v1,v2,Q,tadS)


      CALL getPiSS(g,gp,mtmssm/v2,mbmssm/v1,
     .     mtau/v1,v1,v2,0d0,Q,piSS)


      CALL getPiPP(g,gp,mtmssm/v2,mbmssm/v1,
     .     mtau/v1,v1,v2,0d0,Q,piPP)

*     put all together (except the momentum-dependent corrections)

      vev(1) = v1
      vev(2) = v2
      vev(3) = xx

      do i = 1,3
       do j = 1,3

          MS(i,j) = MStr(i,j) ! start with tree level
          MP(i,j) = MPtr(i,j)

          ! add the full two loop

             MS(i,j) = MS(i,j)+DMS_2lt(i,j)
             MP(i,j) = MP(i,j)+DMP_2lt(i,j)
             MS(i,j) = MS(i,j)+DMS_2lb(i,j)
             MP(i,j) = MP(i,j)+DMP_2lb(i,j)
             MS(i,j) = MS(i,j)+DMS_2ly(i,j)
             MP(i,j) = MP(i,j)+DMP_2ly(i,j)

          ! add zero-mom 1-loop
          ! (If mom=true: Mfull matrices are recalculated)

            MSfull(i,j) = MS(i,j)-piSS(i,j)
            MPfull(i,j) = MP(i,j)-piPP(i,j)

             if(i.eq.j) then

              MSfull(i,j) = MSfull(i,j)+tadS(i)/sqrt(2d0)/vev(i)
              MPfull(i,j) = MPfull(i,j)+tadS(i)/sqrt(2d0)/vev(i)

             endif

       enddo
      enddo

*     simple diagonalization
*    (check of pos. masses^2, and for the mixing matrices)

      CALL jacobi(MSfull,3,ms2,OS)

      CALL jacobi(MPfull,3,ma2,OP)

      CALL mysort(ms2,OS)

      CALL mysort(ma2,OP)

      IF(ms2(1).le.0d0) THEN
       mss(1)=ms2(1)
       err=1
       RETURN
      ENDIF

      IF(ma2(1).le.0d0) THEN
       maa(1)=ma2(1)
       err=2
       RETURN
      ENDIF

      if(mom) then    ! iterative procedure for the EXTERNAL mom

       do kk = 1,3          ! one round for each eigenstate

        ps = sqrt(ms2(kk))

        if(abs(ma2(kk)).le.1d-6) then
         pa = 0d0      ! otherwise the PV FUNCTIONs freak out
        else
         pa = sqrt(ma2(kk))
        endif
      
        IL=0
 500    continue
        IL=IL+1

        CALL getPiSS(g,gp,mtmssm/v2,mbmssm/v1,
     .         mtau/v1,v1,v2,ps,Q,piSS)
    
        CALL getPiPP(g,gp,mtmssm/v2,mbmssm/v1,
     .         mtau/v1,v1,v2,pa,Q,piPP)

        do i=1,3
         do j = 1,3

          MSfull(i,j) = MS(i,j)-piSS(i,j)
          MPfull(i,j) = MP(i,j)-piPP(i,j)

          if(i.eq.j) then

           MSfull(i,j) = MSfull(i,j)+tadS(i)/sqrt(2d0)/vev(i)
           MPfull(i,j) = MPfull(i,j)+tadS(i)/sqrt(2d0)/vev(i)

          endif
         enddo
        enddo

        CALL jacobi(MSfull,3,ms2cor,OS0)
    
        CALL jacobi(MPfull,3,ma2cor,OP0)

        CALL mysort(ms2cor,OS0)
    
        CALL mysort(ma2cor,OP0)

        if(abs(ps**2-abs(ms2cor(kk)))/ps**2
     .         +abs(pa**2-abs(ma2cor(kk)))/pa**2.gt.1d-4) then

         if(IL.ge.10) then
          err=4
          return
         endif

         ps = sqrt(abs(ms2cor(kk)))
         pa = sqrt(abs(ma2cor(kk)))
         goto 500

        else

         ms2(kk) = ms2cor(kk)
         ma2(kk) = ma2cor(kk)

        endif

       enddo
      endif

*     take the square roots and check that it's all right

      do i = 1,3
       if(ms2(i).ge.0d0) then
          mss(i) = sqrt(ms2(i))
       else
          mss(i) = -sqrt(abs(ms2(i)))
          err=1
          return
       endif
      enddo

      do i = 1,3
       if(ma2(i).ge.0d0) then
          maa(i) = sqrt(ma2(i))
       else
          maa(i) = -sqrt(abs(ma2(i)))
          err=2
          return
       endif
      enddo

* Charged Higgs mass (by courtesy of P. Slavich and K.H. Phan)

      mhcsqtree = (l*xx*As-l**2*v1*v2)/sb/cb + mw**2

* UE: Estimate 2-loop corrs. borrowed from the CP-odd Higgs:

      mhcsqtree = mhcsqtree +
     .     (DMP_2lt(1,2)+DMP_2lb(1,2)+DMP_2ly(1,2))/sb/cb

      mhcsqcorr = mhcsqtree

      phc = sqrt(dabs(mhcsqcorr))

      call getPiHpHm(g,gp,l,k,mtmssm/v2,mbmssm/v1,
     .        mtau/v1,v1,v2,xx,Al,At,Ab,Atau,phc,Q,piHpHm)

      mhcsqcorr = mhcsqtree - piHpHm 
     .        + sb**2*tadS(1)/v1/sqrt(2d0) + cb**2*tadS(2)/v2/sqrt(2d0)
     
      phc = sqrt(mhcsqcorr)
     
      if(mom) then

       IL=0
 600   continue
       IL=IL+1

       call getPiHpHm(g,gp,l,k,mtmssm/v2,mbmssm/v1,
     .           mtau/v1,v1,v2,xx,Al,At,Ab,Atau,phc,Q,piHpHm)

       mhcsqcorr = mhcsqtree - piHpHm
     .        + sb**2*tadS(1)/v1/sqrt(2d0) + cb**2*tadS(2)/v2/sqrt(2d0)

        if(abs(phc**2-abs(mhcsqcorr))/phc**2.gt.1d-4)then

         if(IL.ge.10) then
          err=4
          return
         endif

         phc = sqrt(abs(mhcsqcorr))
         goto 600

        endif
       endif

*     now take the square root (if allowed)

       if(mhcsqcorr.gt.0d0) then
          mhc = sqrt(mhcsqcorr)
       else
          mhc=-dsqrt(dabs(mhc))
          err=5
       endif

      end

*
***********************************************************************
*

      SUBROUTINE mysort(msd,ZS)

      IMPLICIT NONE

      DOUBLE PRECISION msd(3),ZS(3,3),mss(3),OS(3,3)
      INTEGER i,j,kk

*     order the disorder

      mss(1) = dmin1(msd(1),msd(2),msd(3))
      mss(3) = dmax1(msd(1),msd(2),msd(3))

      do i=1,3
       if(msd(i).gt.mss(1).and.msd(i).lt.mss(3)) then
          mss(2) = msd(i)
       endif
      enddo

      do i = 1,3
       do j = 1,3
          if(mss(i).eq.msd(j)) then
             do kk = 1,3
              OS(i,kk) = ZS(kk,j)
             enddo
          endif
       enddo
      enddo

      do i=1,3
       msd(i) = mss(i)
       do j=1,3
          ZS(i,j)=OS(i,j)
       enddo
      enddo

      end

*
***********************************************************************
*

      SUBROUTINE getdVB(g,gp,mzpole,mwpole,Q,dVB)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,mzpole,mwpole,Q,dVB
      DOUBLE PRECISION ch2,sh2,c2,s2,pi,rho

      pi = 4d0*atan(1d0)

      ch2 = g**2/(g**2+gp**2)
      sh2 = 1d0-ch2

      c2 = (mwpole/mzpole)**2
      s2 = 1d0-c2

      rho = c2/ch2

      dVB = 6d0+log(c2)/s2*(3.5d0-2.5d0*s2-sh2*(5d0-1.5d0*c2/ch2))
      dVB = dVB-4d0*log(mzpole**2/Q**2)
      dVB = g**2/16d0/pi**2*dVB
*     dVB = dVB*rho

      end

*
***********************************************************************
*

      SUBROUTINE effpot(lp,mt,mg,T1,T2,st,ct,q2,tanb,vv,l,xx,as,
     .     DMS,DMP)

      IMPLICIT NONE

      INTEGER i,j,lp
      DOUBLE PRECISION mt,mg,T1,T2,st,ct,q2,tanb,vv,l,xx,as,
     .     DMS(3,3),DMP(3,3)
      DOUBLE PRECISION c2t,s2t,At,mu,Xt,ht,sbe,pi,k
      DOUBLE PRECISION F1t,F2t,F3t,Ft,FA
      DOUBLE PRECISION DDSIN

      pi = 4d0*atan(1d0)

      if(lp.eq.1) then
       k = 3d0/(16d0*Pi**2)   ! one-loop factor
      elseif(lp.eq.2) then
       k = as/(16d0*Pi**3)    ! two-loop factor
      else
       k = 0d0
      endif

      s2t = 2d0*ct*st
      c2t = ct**2-st**2

      mu = l*xx
      Xt = (T1-T2)*s2t/2d0/mt
      At = Xt+mu/tanb

      sbe = DDSIN(datan(tanb))

      ht = mt/vv/sbe          ! v ~ 174

      if(lp.eq.1) then        !the usual one-loop FUNCTIONs

       Ft = T1*(log(T1/q2)-1d0)-T2*(log(T2/q2)-1d0)

       F1t = log(T1*T2/mt**4)

       F2t = log(T1/T2)

       F3t = 2d0-(T1+T2)/(T1-T2)*log(T1/T2)

       FA = Ft

      elseif(lp.eq.2) then

*UE: FA from makefuncs explodes for At=0
      if(dabs(At).le.1d-12) At=1d-12

      CALL makefuncs(mt,mg,T1,T2,s2t,c2t,q2,tanb,At,mu,
     .      F1t,F2t,F3t,Ft,FA)

      endif

*     now build up the results

      DMS(1,1) = .5d0*ht**2*mu**2*s2t**2*F3t
     .    +ht**2*tanb*mu*At/(T1-T2)*Ft

      DMS(1,2) =-ht**2*mu*mt*s2t*F2t-.5d0*ht**2*At*mu*s2t**2*F3t
     .    -ht**2*mu*At/(T1-T2)*Ft

      DMS(1,3) = .5d0*ht*l*mu*mt*s2t**2/tanb*F3t
     .    -ht*l*mt*(At-2d0*mu/tanb)/(T1-T2)*Ft

      DMS(2,2) = 2d0*ht**2*mt**2*F1t+2d0*ht**2*At*mt*s2t*F2t
     .    +.5d0*ht**2*At**2*s2t**2*F3t
     .    +ht**2/tanb*mu*At/(T1-T2)*Ft

      DMS(2,3) = -.5d0*ht*l*At*mt*s2t**2/tanb*F3t
     .    -ht*l*mt**2*s2t/tanb*F2t-ht*l*mt*At/(T1-T2)/tanb*Ft

      DMS(3,3) = .5d0*l**2*s2t**2*mt**2/tanb**2*F3t
     .     +l**2*mt**2/tanb*At/mu/(T1-T2)*Ft

      DMS(2,1) = DMS(1,2)
      DMS(3,1) = DMS(1,3)
      DMS(3,2) = DMS(2,3)

      DMP(1,1) = ht**2*mu*At/(T1-T2)*FA*tanb

      DMP(1,2) = ht**2*mu*At/(T1-T2)*FA

      DMP(1,3) = l*ht*mt*At/(T1-T2)*FA

      DMP(2,2) = ht**2*mu*At/(T1-T2)*FA/tanb

      DMP(2,3) = l*ht*mt*At/(T1-T2)*FA/tanb

      DMP(3,3) = l**2*mt**2*At/mu/(T1-T2)*FA/tanb

      DMP(2,1) = DMP(1,2)
      DMP(3,1) = DMP(1,3)
      DMP(3,2) = DMP(2,3)

      do i=1,3
       do j=1,3
          DMS(i,j) = k*DMS(i,j)
          DMP(i,j) = k*DMP(i,j)
       enddo
      enddo

      end

*
***********************************************************************
*

      SUBROUTINE makefuncs(mt,mg,T1,T2,s2t,c2t,q,tanb,At,mu,
     .     F1t,F2t,F3t,Ft,FA)

      IMPLICIT NONE

      DOUBLE PRECISION mt,mg,T1,T2,s2t,c2t,q,tanb,At,mu,
     .     F1t,F2t,F3t,Ft,FA

      DOUBLE PRECISION DT1,DT2,Dc2t,DT1T1,DT2T2,Dtt,Dc2tc2t,
     .     DT1t,DT2t,DT1T2,Dtc2t,DT1c2t,DT2c2t,Dcptmptt,
     .     Dtt_1,Dc2t_1,Dc2tc2t_1,Dtc2t_1,Dcptmptt_1,
     .     Dtt_2,Dc2t_2,Dc2tc2t_2,Dtc2t_2,Dcptmptt_2


      CALL makederiv(mt,mg,T1,T2,s2t,c2t,q,
     .     DT1,DT1T1,DT1t,DT1c2t,DT1T2,
     .     Dtt_1,Dc2t_1,Dc2tc2t_1,Dtc2t_1,Dcptmptt_1)


      CALL makederiv(mt,mg,T2,T1,-s2t,c2t,q,
     .     DT2,DT2T2,DT2t,DT2c2t,DT1T2,
     .     Dtt_2,Dc2t_2,Dc2tc2t_2,Dtc2t_2,Dcptmptt_2)

      Dtt = Dtt_1+Dtt_2
      Dc2t = Dc2t_1+Dc2t_2
      Dc2tc2t = Dc2tc2t_1+Dc2tc2t_2
      Dtc2t = Dtc2t_1+Dtc2t_2
      Dcptmptt = Dcptmptt_1+Dcptmptt_2

      F1t = Dtt+DT1T1+DT2T2+2d0*(DT1t+DT2t+DT1T2)

      F2t = DT1T1-DT2T2+DT1t-DT2t
     .     -4d0*c2t**2/(T1-T2)*(Dtc2t+DT1c2t+DT2c2t)

      F3t = DT1T1+DT2T2-2d0*DT1T2
     .    -2d0/(T1-T2)*(DT1-DT2)
     .    +16d0*c2t**2/(T1-T2)**2*(c2t**2*Dc2tc2t+2d0*Dc2t)
     .     -8d0*c2t**2/(T1-T2)*(DT1c2t-DT2c2t)

      Ft = DT1-DT2-4d0*c2t**2/(T1-T2)*Dc2t

      FA = Ft-2d0*mu/tanb/At/(T1-T2)/s2t**2*Dcptmptt

      end

*
***********************************************************************
*

      SUBROUTINE makederiv(mt,mg,T1,T2,s2t,c2t,q,
     .     DT1,DT1T1,DT1t,DT1c2t,DT1T2,Dtt,Dc2t,Dc2tc2t,Dtc2t,Dcptmptt)

      IMPLICIT NONE

      DOUBLE PRECISION mt,mg,T1,T2,s2t,c2t,q,
     .     DT1,DT1T1,DT1t,DT1c2t,DT1T2,Dtt,Dc2t,Dc2tc2t,Dtc2t,Dcptmptt
      DOUBLE PRECISION delt,phi,II,JJ

      DOUBLE PRECISION t,g,Logt,Logg,LogT1,LogT2,pphi,del,III

      t = mt**2
      g = mg**2

      Logt = Log(t/q)
      Logg = Log(g/q)
      LogT1 = Log(T1/q)
      LogT2 = Log(T2/q)
      pphi = phi(T1,g,t)
      del = delt(T1,g,t)
      III = II(q,T1,g,t)

      Dc2t = .5d0*JJ(q,T1,T1)-.5d0*JJ(q,T1,T2)+2d0*mg*mt/s2t*III

      Dc2tc2t = mg*mt/s2t**3*III

      Dcptmptt = -4d0*mg*mt*s2t*III

      DT1= -6d0*T1+2d0*mg*mt*s2t+4d0*t*(1d0-logt+logT1)
     .     +4d0*g*(1d0-logg+logT1)+((5d0-c2t**2)*T1-s2t**2*T2
     .     -4d0*mg*mt*s2t)*logT1+(-3d0+c2t**2)*T1*logT1**2
     .     +s2t**2*T2*logT1*logT2-(2d0*(g+t-T1)-2d0*mg*mt*s2t)
     .     *(logt*(logT1-logg)+logT1*logg)+(2d0/t*(del+2d0*g*t)
     .     -2d0*mg/mt*s2t*(g+t-T1))*pphi

      DT1T1= -(1d0+c2t**2)+4d0/T1*(g+t-mg*mt*s2t)-s2t**2*T2/T1
     .     *(1d0-logT2)+(3d0+c2t**2+8d0*g*t/del-4d0*mg*mt*s2t/del
     .     *(g+t-T1))*logT1-4d0*t/del/T1*(del-g*(g-t-T1)+mg*mt*s2t
     .     *(g-t+T1))*logt-4d0*g/del/T1*(del+t*(g-t+T1)-mg*mt*s2t
     .     *(g-t-T1))*logg+(-3d0+c2t**2)*logT1**2
     .     +2d0*(logt*(logT1-logg)+logT1*logg)-2d0/t/del*((g+t-T1)
     .     *(del-2d0*g*t)+4d0*mg**3*mt**3*s2t)*pphi

      DT1c2t= (T2*(1d0-logT2)-T1*(1d0-logT1))*LogT1
     .     -mg*mt/s2t*(1d0-2d0*logT1+logt*(logT1-logg)
     .     +logT1*logg-(g+t-T1)/t*pphi)

      DT1t= mg/mt*s2t+4d0*g/del*(T1-g-t+2d0*mg*mt*s2t)*logg
     .     +4d0/del*(2d0*g*t-mg*mt*s2t*(g+t-T1))*logt+2d0/del
     .     *(2d0*g*(g-t-T1)-mg/mt*s2t*(del-2d0*t*(t-g-T1)))*logT1
     .     +(-2d0+mg/mt*s2t)*(logt*(logT1-logg)+logT1*logg)
     .     +1d0/del/t*(mg/mt*s2t*(del*(T1-g-3d0*t)
     .     +2d0*t*((t-T1)**2-g**2))
     .     +2d0*(g-T1)**3+2d0*t*(del+(2d0*T1-t)*(g+T1)))*pphi

      DT1T2= s2t**2*logT1*logT2

      Dtt= -2d0-5d0/2d0*mg/mt**3*s2t*T1+6d0*logt**2
     .     +4d0*g/del*(g-t-T1+mt/mg*s2t*(t-g-T1))*logt
     .     -4d0*g/del*(g-t+T1+mg/mt*s2t*(t-g+T1))*logg
     .     +(8d0*g*T1/del+2d0*mg/mt**3*s2t*T1
     .      *(1d0-2d0*t/del*(t+g-T1)))*logT1
     .     -(2d0-mg*s2t*T1/2d0/mt**3)*logg*(logt-logT1)
     .     -(2d0+mg*s2t*T1/2d0/mt**3)*logt*logT1
     .     +mg*s2t/2d0/mt**3*(g+3d0*t)*logT1*(logt-logg)
     .     -2d0/del/t*(mg/mt**3*s2t*(del**2/4d0+t*(g-2d0*t+T1)*del
     .     +t**2*(g-t+T1)**2)-T1*(del+(g+t)*(2d0*t-T1))-(g-t)**3)*pphi

      Dtc2t = -mg/mt/s2t/2d0*(5d0*T1
     .     -4d0*T1*logT1+(g-3d0*t)*logT1*(logg-logt)
     .     +T1*(logt*(logT1-logg)+logg*logT1)+(del/t-2d0*(g-t+T1))*pphi)

      end

*
***********************************************************************
*

      DOUBLE PRECISION FUNCTION myA0(m,q)
      DOUBLE PRECISION m,q

      if(m.ne.0d0) then
       myA0 = m*(1d0-Log(m/q))
      else
       myA0 = 0d0
      endif

      end

*
***********************************************************************
*

      DOUBLE PRECISION FUNCTION myB0(p,m1,m2,q)

      IMPLICIT NONE
      DOUBLE PRECISION p, m1, m2
      DOUBLE PRECISION mudim2, divergence, lambda2, q
      DOUBLE PRECISION acc, eps, minacc
      DOUBLE COMPLEX x1, x2, y1, y2, r, be0
      DOUBLE COMPLEX Ieps, onePeps, oneMeps
      COMMON/cutoff/mudim2, divergence, lambda2
      PARAMETER (acc = 1d-12)
      PARAMETER (eps = 1d-20)
      PARAMETER (Ieps = (0d0,1d0)*eps)
      PARAMETER (onePeps = 1d0 + Ieps)
      PARAMETER (oneMeps = 1d0 - Ieps)

      DOUBLE COMPLEX fpv, xlogx
      EXTERNAL fpv, xlogx

      divergence = 0d0
      lambda2 = 0d0
      mudim2 = q
      minacc = acc*(m1 + m2)

* general case
      if(abs(p) .gt. minacc) then
  
      CALL roots(p, m1, m2, x1, x2, y1, y2, r)
        if(abs(y1) .gt. .5d0 .and. abs(y2) .gt. .5d0) then
          be0 = -log(m2/mudim2) -
     +      fpv(1, x1, y1) - fpv(1, x2, y2)
        else if(abs(x1) .lt. 10d0 .and. abs(x2) .lt. 10d0) then
          be0 = 2 - log(p*oneMeps/mudim2) +
     +      xlogx(-x1) + xlogx(-x2) - xlogx(y1) - xlogx(y2)
        else if(abs(x1) .gt. .5d0 .and. abs(x2) .gt. .5d0) then
          be0 = -log(m1/mudim2) -
     +      fpv(1, y1, x1) - fpv(1, y2, x2)
        else
          be0 = 1d100
        endif

* zero momentum
      else if(abs(m1 - m2) .gt. minacc) then
        x2 = oneMeps*m1/(m1 - m2)
        y2 = oneMeps*m2/(m2 - m1)
        if(abs(y2) .gt. .5d0) then
          be0 = -log(m2/mudim2) - fpv(1, x2, y2)
        else
          be0 = -log(m1/mudim2) - fpv(1, y2, x2)
        endif
      else
        be0 = -log(m2/mudim2)
      endif

      myB0 = dble(be0 + divergence)

      end

*
**********************************************************************
*

      DOUBLE PRECISION FUNCTION myB1(p,m1,m2,q)

      IMPLICIT NONE

      DOUBLE PRECISION p,m1,m2,q
      DOUBLE PRECISION myA0,myB0

      if(p.eq.0d0) then
       if(abs(m1-m2).le.1d-8) then
          myB1 = -Log(m1/q)/2d0
       else
          if(m1.eq.0d0) then
             myB1 = (1d0-2d0*Log(m2/q))/4d0
          elseif(m2.eq.0d0) then
             myB1 = (3d0-2d0*Log(m1/q))/4d0
          else
             myB1 = (1d0-Log(m2/q)+m1**2/(m1-m2)**2*Log(m2/m1)
     .            +(m1+m2)/(m1-m2)/2d0)/2d0
          endif
       endif
      else
       myB1 = (myA0(m2,q)-myA0(m1,q)+(p+m1-m2)*myB0(p,m1,m2,q))/2d0/p
      endif

      end

*
**********************************************************************
*

      DOUBLE PRECISION FUNCTION myF(q,m1,m2,mu2)

      IMPLICIT NONE
      DOUBLE PRECISION q,m1,m2,mu2,myA0,myB0

      myF = myA0(m1,mu2)-2d0*myA0(m2,mu2)
     .     -(2d0*q+2d0*m1-m2)*myB0(q,m1,m2,mu2)

      end

*
***********************************************************************
*

      DOUBLE PRECISION FUNCTION myG(q,m1,m2,mu2)

      IMPLICIT NONE
      DOUBLE PRECISION q,m1,m2,mu2,myA0,myB0

      if(q.eq.0d0.and.m1.eq.0d0.and.m2.eq.0d0) then
       myG = 0d0
      else
       myG = (q-m1-m2)*myB0(q,m1,m2,mu2)-myA0(m1,mu2)-myA0(m2,mu2)
      endif

      end

*
***********************************************************************
*

      DOUBLE PRECISION FUNCTION myB22(q,m1,m2,mu2)

      IMPLICIT NONE
      DOUBLE PRECISION q,m1,m2,mu2,myA0,myB0,myB1

      if(q.eq.0d0.and.m1.eq.0d0.and.m2.eq.0d0) then
       myB22 = 0d0
      else
       myB22 = ((myA0(m1,mu2)+myA0(m2,mu2))/2d0
     .      +(m1+m2-q/2d0)*myB0(q,m1,m2,mu2)
     .      +(m2-m1)*(myB1(q,m1,m2,mu2)-myB0(q,m1,m2,mu2)/2d0)
     .      +m1+m2-q/3d0)/6d0
      endif

      end

*
***********************************************************************
*

      DOUBLE PRECISION FUNCTION myB22T(q,m1,m2,mu2)

      IMPLICIT NONE
      DOUBLE PRECISION q,m1,m2,mu2,myA0,myB22

      myB22T = myB22(q,m1,m2,mu2)-myA0(m1,mu2)/4d0-myA0(m2,mu2)/4d0

      end

*
**********************************************************************
*

      DOUBLE PRECISION FUNCTION myH(q,m1,m2,mu2)

      IMPLICIT NONE
      DOUBLE PRECISION q,m1,m2,mu2,myG,myB22

      myH = 4d0*myB22(q,m1,m2,mu2)+myG(q,m1,m2,mu2)

      end

*
**********************************************************************
*

      DOUBLE PRECISION FUNCTION JJ(q,m1,m2)

      IMPLICIT NONE
      DOUBLE PRECISION q,m1,m2

      JJ = m1*m2*(Log(m1/q)-1d0)*(Log(m2/q)-1d0)

      end

*
**********************************************************************
*

      DOUBLE PRECISION FUNCTION II(q,m1,m2,m3)

      IMPLICIT NONE
      DOUBLE PRECISION q,m1,m2,m3,delt,phi

      II = (m1-m2-m3)/2d0*Log(m2/q)*Log(m3/q)
     .     +(m2-m1-m3)/2d0*Log(m1/q)*Log(m3/q)
     .     +(m3-m1-m2)/2d0*Log(m1/q)*Log(m2/q)
     .     +2d0*m1*log(m1/q)+2d0*m2*Log(m2/q)+2d0*m3*Log(m3/q)
     .     -2.5d0*(m1+m2+m3)-delt(m1,m2,m3)/2d0/m3*phi(m1,m2,m3)

      end

*
**********************************************************************
*

      FUNCTION phi(x,y,z)

*     from Davydychev and Tausk, Nucl. Phys. B397 (1993) 23

      IMPLICIT NONE
      DOUBLE PRECISION x,y,z,phi,pphi,myphi

      if(x.le.z.and.y.le.z) then
       pphi = myphi(x,y,z)
      elseif(z.le.x.and.y.le.x) then
       pphi = z/x*myphi(z,y,x)
      elseif(z.le.y.and.x.le.y) then
       pphi = z/y*myphi(z,x,y)
      endif

      phi = pphi

      end

      FUNCTION myphi(x,y,z)

      IMPLICIT NONE

      DOUBLE PRECISION x,y,z,myphi
      DOUBLE PRECISION u,v
      DOUBLE PRECISION Pi,SLLi2
      DOUBLE COMPLEX clam,cxp,cxm,SLCLI2,ccphi

      PARAMETER (pi = 3.1415926535897932384626433832795029d0)

*     auxiliary variables

      u = x/z
      v = y/z

      if(u.le.1d-8) then

       if(v.ne.1d0) then
          myphi = (log(u)*log(v)+2d0*SLLi2(1d0-v))/(1d0-v)
       else
          myphi = 2d0-log(u)
       endif

      elseif(v.le.1d-8) then

       if(u.ne.1d0) then
          myphi = (log(v)*log(u)+2d0*SLLi2(1d0-u))/(1d0-u)
       else
          myphi = 2d0-log(v)
       endif

      else

       if((1d0-u-v)**2.ge.4d0*u*v) then
          clam = DCMPLX(sqrt((1d0-u-v)**2-4d0*u*v),0d0)
       else
          clam = DCMPLX(0d0,sqrt(4d0*u*v-(1d0-u-v)**2))
       endif

       cxp = (1d0+(u-v)-clam)/2d0
       cxm = (1d0-(u-v)-clam)/2d0

*     phi FUNCTION from eq. (A4)

       ccphi = (2d0*log(cxp)*log(cxm)-log(u)*log(v)-
     .      2d0*(SLCLI2(cxp)+SLCLI2(cxm))+Pi**2/3d0)/clam
       myphi = DBLE(ccphi)

      endif

      end

*
***********************************************************************
*

      FUNCTION SLLi2(x)

      IMPLICIT NONE

      DOUBLE COMPLEX SLCLI2,z
      DOUBLE PRECISION x,SLLi2

      z = DCMPLX(x,0d0)
      SLLi2 = DBLE(SLCLI2(z))

      end

*
***********************************************************************
*

      DOUBLE COMPLEX FUNCTION SLCLI2(Z)

*     just CALL the Dilog routine

      DOUBLE COMPLEX Z,Dilog

      SLCLI2 = Dilog(Z)

      end

*
**********************************************************************
*
* Dilog.F
* complex dilogarithm
* this file is part of FeynHiggs
* last modified 20 Oct 05 th

      DOUBLE COMPLEX FUNCTION Dilog(z)
      IMPLICIT NONE
      DOUBLE COMPLEX z

      DOUBLE COMPLEX Dilogsum
      EXTERNAL Dilogsum

      DOUBLE PRECISION absz, abs1z
      DOUBLE COMPLEX t, mlogz

      DOUBLE PRECISION pi, zeta2
      PARAMETER (pi = 3.1415926535897932384626433832795029d0)
      PARAMETER (zeta2 = pi*pi/6d0)

      absz = abs(z)
      if( absz .lt. 1d-20 ) then
       Dilog = -log(1d0-z)
       return
      endif

      abs1z = abs(1d0-z)
      if( abs1z .lt. 1d-20 ) then
         Dilog = zeta2
         return
      endif

      if( DBLE(z) .gt. .5d0 ) then
         mlogz = -log(z)
         t = zeta2+mlogz*log(1d0-z)
         if( abs1z .gt. 1 ) then
            Dilog = Dilogsum(log(1d0-1d0/z))+zeta2 +
     .           .5d0*log(z-1d0)**2+t
         else
          Dilog = -Dilogsum(mlogz)+t
       endif
      else
       if( absz .gt. 1 ) then
          Dilog = -Dilogsum(-log(1d0-1d0/z))-zeta2-.5d0*log(-z)**2
       else
          Dilog = Dilogsum(-log(1d0-z))
       endif
      endif
      end

************************************************************************

      DOUBLE COMPLEX FUNCTION Dilogsum(w)
      IMPLICIT NONE
      DOUBLE COMPLEX w

      DOUBLE COMPLEX u, t
      INTEGER k

      DOUBLE PRECISION b2, b4, b6, b8, b10, b12, b14
      DOUBLE PRECISION b16, b18, b20, b22
      PARAMETER (b2 = 1d0/6d0)
      PARAMETER (b4 = -1d0/30d0)
      PARAMETER (b6 = 1d0/42d0)
      PARAMETER (b8 = -1d0/30d0)
      PARAMETER (b10 = 5d0/66d0)
      PARAMETER (b12 = -691d0/2730d0)
      PARAMETER (b14 = 7d0/6d0)
      PARAMETER (b16 = -3617d0/510d0)
      PARAMETER (b18 = 43867d0/798d0)
      PARAMETER (b20 = -174611d0/330d0)
      PARAMETER (b22 = 854513d0/138d0)

      DOUBLE PRECISION bernoulliB(11)
      data bernoulliB /b2, b4, b6, b8, b10, b12, b14,
     .     b16, b18, b20, b22/

      Dilogsum = w*(1d0-.25d0*w)
      if( abs(w) .lt. 1d-10 ) return

      u = w
      do k = 1, 11
       u = u*w**2/DBLE(2d0*k*(2d0*k+1d0))
       t = u*bernoulliB(k)
       Dilogsum = Dilogsum+t
       if( abs(t) .lt. 1d-16*abs(Dilogsum) ) return
      enddo

      end

      FUNCTION delt(x,y,z)
      DOUBLE PRECISION delt,x,y,z

      delt = x**2+y**2+z**2-2d0*(x*y+x*z+y*z)

      end

***********************************************************************
*     FUNCTION IN THE HALL-RATTAZZI-SARID TERM
***********************************************************************

      FUNCTION SLH2(x,y)

      IMPLICIT NONE

      DOUBLE PRECISION x,y,SLH2,eps

      eps = 1d-8

      if(abs(x-y).ge.eps) then
       if(abs(x-1d0).ge.eps.and.abs(y-1d0).ge.eps) then
          SLH2 = x*log(x)/(1d0-x)/(x-y)+y*log(y)/(1d0-y)/(y-x)
       elseif(abs(x-1d0).ge.eps.and.abs(y-1d0).lt.eps) then
          SLH2 = (-1d0+x-x*log(x))/(x-1d0)**2
       elseif(abs(x-1d0).lt.eps.and.abs(y-1d0).ge.eps) then
          SLH2 = (-1d0+y-y*log(y))/(y-1d0)**2
       else
          SLH2 = -.5d0
       endif
      else
       if(abs(x-1d0).ge.eps) then
          SLH2 = (1d0-x+log(x))/(x-1d0)**2
       else
          SLH2 = -.5d0
       endif
      endif

      end

*
***********************************************************************
*

      SUBROUTINE diagonalize(mq,mw,mz,msql,msqr,Aq,mu,tb,iq,
     .     msq2,sth,cth)

      IMPLICIT NONE

      INTEGER iq
      DOUBLE PRECISION mq,mw,mz,msql,msqr,Aq,mu,tb,msq2(2),sth,cth

      DOUBLE PRECISION mq2,mz2,mw2,c2b,xq,yq,zq,tth,dx,dy

      mq2 = mq**2
      mz2 = mz**2
      mw2 = mw**2
      c2b = (1d0-tb**2)/(1d0+tb**2)

      if(iq.eq.1) then
       zq = mq*(Aq-mu/tb)
       dx = 1d0/4d0*mz2*c2b
       dy = 1d0/12d0*(8d0*mw2-5d0*mz2)*c2b
      elseif(iq.eq.2) then
       zq = mq*(Aq-mu*tb)
       dx = -1d0/4d0*mz2*c2b
       dy = -1d0/12d0*(4d0*mw2-mz2)*c2b
      else
       write(*,*) 'ERROR: iq out of range'
      endif

      xq = mq2+1d0/2d0*(msql**2+msqr**2)+dx
      yq = 1d0/2d0*(msql**2-msqr**2)+dy

      msq2(1) = xq+sqrt(yq**2+zq**2)
      msq2(2) = xq-sqrt(yq**2+zq**2)

      if(zq.eq.0.and.yq.ge.0)then
       cth=1d0
       sth=0d0
      elseif(zq.eq.0.and.yq.lt.0)then
       cth=0d0
       sth=1d0
      else
       tth = 1d0/zq*(sqrt(yq**2+zq**2)-yq)
       cth = 1d0/sqrt(1d0+tth**2)
       sth = cth*tth
      endif

      end

*
***********************************************************************
*

      SUBROUTINE squarks(mtsm,mbsm,mz,mw,mqs,mts,mbs,At,Ab,mu,tb,mg,as,
     .     A0,vv,q2,mtmssm,mbmssm,mstop2,msbot2,cst,sst,csb,ssb,
     .     errsqua,mtpole,asma,ast,v2)

*     compute the quark running masses and the squark masses and mixing

      IMPLICIT NONE

      DOUBLE PRECISION mtsm,mbsm,mz,mw,mqs,mts,mbs,At,Ab,mu,tb,mg,as,
     .     A0,vv,q2,mtmssm,mbmssm,mstop2(2),msbot2(2),cst,sst,csb,ssb,
     .     mtpole,asma,ast,v2,htmt,pi,dla,sb2,fact,mtma,htma,dlqa

      LOGICAL errsqua

      DOUBLE PRECISION dmt,eb,dd,K,K0,eps

      INTEGER i,imax

      pi=4d0*atan(1d0)

      errsqua = .false.

*     first computation of the squarks using the SM running quark masses

      mtmssm = mtsm
      mbmssm = mbsm

      CALL diagonalize(mtmssm,mw,mz,mqs,mts,At,mu,tb,1,
     .     mstop2,sst,cst)

      CALL diagonalize(mbmssm,mw,mz,mqs,mbs,Ab,mu,tb,2,
     .     msbot2,ssb,csb)
      
      eps=1d-4
      imax=100
      i=0
      K=0d0

!      PRINT*,"          I","           K","                      DeltaK"
 1    I=I+1      
      K0=K

*     compute the threshold corrections to the quark masses
      CALL getdmt(mtmssm,mstop2,sst,cst,mbmssm,msbot2,ssb,csb,
     .      vv,mg,A0,mu,tb,as,q2,dmt)

      CALL getdmb(mtmssm,mstop2,sst,cst,mbmssm,msbot2,ssb,csb,
     .      vv,mg,A0,mu,tb,as,q2,dd,eb)

*     recompute the quark masses

       mtmssm = mtsm+dmt
      
* mod. by UE: resum large logs ~ht^2*ln(Q^2/mt^2) (not included in dmt):

       htmt=mtpole/v2/(1d0+4d0*ast/(3d0*pi)+10.9d0*(ast/pi)**2)
       dla=(asma/ast)**(1d0/7d0)
       sb2=tb**2/(1d0+tb**2)
       fact=1d0-9d0*sb2*htmt**2/(8d0*pi*ast)*(1d0-dla)

       mtma=mtmssm*fact**(-1d0/6d0)
       htma=htmt*dla**4*fact**(-0.5d0)
       dlqa=(as/asma)**(1d0/7d0)
       mtmssm=mtma*(1d0-9d0*htma**2/(8d0*pi*asma)*
     .    (1d0-dlqa))**(-1d0/6d0)

* end mod. by UE

       K = (1d0-dd)/(1d0+eb*tb) ! "resummation"
!      PRINT*,I,K,DABS(K0-K)/K
      if(DABS(K0-K)/K.LT.eps)i=imax

       mbmssm = K*mbsm

*     now recompute the squark masses and mixing

 
      CALL diagonalize(mtmssm,mw,mz,mqs,mts,At,mu,tb,1,
     .      mstop2,sst,cst)

 
      CALL diagonalize(mbmssm,mw,mz,mqs,mbs,Ab,mu,tb,2,
     .      msbot2,ssb,csb)

      if(.not.mstop2(2).ge.0d0) then
       errsqua = .true.
       return
      endif

      if(.not.msbot2(2).ge.0d0) then
       errsqua = .true.
       return
      endif

      if(I.LT.imax)goto 1

      if(DABS(K0-K)/K.GT.eps)errsqua = .true.

!      PRINT*,"err =",errsqua
!      PRINT*,""
 
      end

*
***********************************************************************
*

      SUBROUTINE getdmt(mt,mstop2,sst,cst,mb,msbot2,ssb,csb,
     .      vv,mg,A0,mu,tb,as,q2,dmt)

      IMPLICIT NONE

      DOUBLE PRECISION mt,mstop2(2),sst,cst,mb,msbot2(2),ssb,csb,
     .     vv,mg,A0,mu,tb,as,q2,dmt

      DOUBLE PRECISION pi,cbe,sbe,ht,hb,mt2,mb2,mu2,mg2,T1,T2,B1,B2,
     .     dmts,dmty,myB0,myB1,cbe2,sbe2,ht2,hb2

      pi=4d0*atan(1d0)

      cbe = 1d0/sqrt(1d0+tb**2)
      sbe = tb*cbe

      cbe2 = cbe**2
      sbe2 = sbe**2

      ht = mt/vv/sbe
      hb = mb/vv/cbe

      ht2 = ht**2
      hb2 = hb**2

      T1 = mstop2(1)
      T2 = mstop2(2)

      B1 = msbot2(1)
      B2 = msbot2(2)

      mt2 = mt**2
      mb2 = mb**2
      mu2 = mu**2
      mg2 = mg**2

      dmts = as/3d0/pi*mt*(
     .     myB1(mt2,mg2,T1,q2)+myB1(mt2,mg2,T2,q2)
     .     -2d0*sst*cst*mg/mt*(myB0(mt2,mg2,T1,q2)-myB0(mt2,mg2,T2,q2)))
      
      IF(A0.ge.0d0) THEN
       dmty = mt/32d0/pi**2*(
* mod. by UE:
* Scale q of myB1 replaced; large logs are resummed in SUBROUTINE
* squarks
     .     ht2*(
c Before 18/6/2014:
c     .       2d0*myB1(mt2,mt2,0d0,mt2)+myB1(mt2,mb2,0d0,mt2)
c     .       +2d0*cbe2*(myB1(mt2,mt2,A0,A0)-myB1(mt2,mt2,0d0,mt2))
c     .       +cbe2*(myB1(mt2,mb2,A0,A0)-myB1(mt2,mb2,0d0,mt2)
c After 18/6/2014 (the difference is added to runt()):
     .       +2d0*cbe2*myB1(mt2,mt2,A0,A0)
     .       +cbe2*myB1(mt2,mb2,A0,A0)
* end mod. by UE
     .       +myB1(mt2,mu2,T1,q2)+myB1(mt2,mu2,T2,q2))
     .     +hb2*(cbe2*myB1(mt2,mb2,0d0,q2)+sbe2*myB1(mt2,mb2,A0,q2))
     .     -2d0*ht*hb*sbe*cbe*mb/mt*
     .     (myB0(mt2,mb2,0d0,q2)-myB0(mt2,mb2,A0,q2))
     .     +(ht2*csb**2+hb2*ssb**2)*myB1(mt2,mu2,B1,q2)
     .     +(ht2*ssb**2+hb2*csb**2)*myB1(mt2,mu2,B2,q2)
     .     +2d0*ht*hb*ssb*csb*mu/mt*
     .     (myB0(mt2,mu2,B1,q2)-myB0(mt2,mu2,B2,q2)))
      ELSE
       dmty=0d0
      endif


      dmt = dmts+dmty

      end

*
***********************************************************************
*

      SUBROUTINE getdmb(mt,mstop2,sst,cst,mb,msbot2,ssb,csb,
     .      vv,mg,A0,mu,tb,as,q2,dd,eb)

      IMPLICIT NONE

      DOUBLE PRECISION mt,mstop2(2),sst,cst,mb,msbot2(2),ssb,csb,
     .     vv,mg,A0,mu,tb,as,q2,dd,eb

      DOUBLE PRECISION pi,cbe,sbe,ht,hb,mt2,mu2,mg2,T1,T2,B1,B2,
     .     myB0,myB1,ht2,hb2,At,Ab,dds,ebs,ddy,eby,
     .     lhiggs,hhiggs,higgsino

      pi=4d0*atan(1d0)

      cbe = 1d0/sqrt(1d0+tb**2)
      sbe = tb*cbe

      ht = mt/vv/sbe
      hb = mb/vv/cbe

      ht2 = ht**2
      hb2 = hb**2

      T1 = mstop2(1)
      T2 = mstop2(2)


      B1 = msbot2(1)
      B2 = msbot2(2)

      mt2 = mt**2
      mu2 = mu**2
      mg2 = mg**2

      At = sst*cst*(T1-T2)/mt+mu/tb
      Ab = ssb*csb*(B1-B2)/mb+mu*tb

*     the strong part

*UE.10.9.2014
      IF(DABS(B1-B2).GE.1d-9) THEN
        dds = as/3d0/pi*(
     .    2d0*Ab*mg/(B1-B2)*(myB0(0d0,mg2,B1,q2)-myB0(0d0,mg2,B2,q2))
     .    -myB1(0d0,mg2,B1,q2)-myB1(0d0,mg2,B2,q2))
        ebs = -as/3d0/pi*(
     .    2d0*mu*mg/(B1-B2)*(myB0(0d0,mg2,B1,q2)-myB0(0d0,mg2,B2,q2)))
      ELSE
        dds = as/3d0/pi*(
     .    -1d0*Ab*mg/B1)
        ebs = -as/3d0/pi*(
     .    -2d0*mu*mg/B1)
      ENDIF

*     the Yukawa part (in the large-tanB limit)

      lhiggs = -ht2/4d0*(5d0-6d0*Log(mt2/q2))

      IF(A0.gt.0d0) THEN
       hhiggs = 2d0*ht2*(mt2-A0+A0*Log(A0/q2)-mt2*Log(mt2/q2))/(mt2-A0)
     .    +hb2*mt2/2d0/(mt2-A0)**2*(mt2-A0+mt2*Log(A0/mt2))
     .    +3d0*hb**2/4d0*(1d0-2d0*Log(A0/q2))
      ELSE
       hhiggs=0d0
      ENDIF

*UE.10.9.2014
      IF(DABS(T1-T2).GE.1d-9) THEN
        higgsino = hb**2*(myB1(0d0,mu2,B1,q2)+myB1(0d0,mu2,B2,q2))
     .    +(hb2*cst**2+ht2*sst**2)*myB1(0d0,mu2,T1,q2)
     .    +(hb2*sst**2+ht2*cst**2)*myB1(0d0,mu2,T2,q2)
     .    -2d0*ht2*mu2/(T1-T2)*(myB0(0d0,mu2,T1,q2)-myB0(0d0,mu2,T2,q2))
        eby = ht**2/16d0/pi**2*At*mu/(T1-T2)*
     .    (T1/(T1-mu2)*Log(T1/mu2)-T2/(T2-mu2)*Log(T2/mu2))
      ELSE
        higgsino = hb**2*(myB1(0d0,mu2,B1,q2)+myB1(0d0,mu2,B2,q2))
     .    +(hb2*cst**2+ht2*sst**2)*myB1(0d0,mu2,T1,q2)
     .    +(hb2*sst**2+ht2*cst**2)*myB1(0d0,mu2,T2,q2)
     .    +2d0*ht2*mu2/T1
        eby = ht**2/16d0/pi**2*At*mu*(2D0/T1-1D0/(T1-MU2))
      ENDIF

      ddy = -(lhiggs+hhiggs+higgsino)/32d0/pi**2

*     all together

      dd = dds+ddy
      eb = ebs+eby

      end

*
***********************************************************************
*

      DOUBLE PRECISION FUNCTION runt(x,tc,asc,zm,vv,bmass)

*     compute the running (SM, MSbar) top mass

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION mt2,mb2,myB1

      pi=4d0*atan(1d0)

      if(x.le.tc)then
       fn=5d0
      else
       fn=6d0
      endif

      b0=11d0-2d0*fn/3d0
      b1=102d0-38d0*fn/3d0
      g0=8d0
      g1=404d0/3d0-40d0*fn/9d0

      asx=slasf(x,tc,asc,zm)
      ast=slasf(tc,tc,asc,zm)
      rrr=tc*(asx/ast)**(g0/(2d0*b0))

*     this is the relation between mpole/mtrun(mtrun)
*      pol1=1d0+4d0*ast/(3d0*pi)+8.243d0*(ast/pi)**2

*     this is the relation between mpole/mt(mpole)
      pol2= 1d0+4d0*ast/(3d0*pi) +10.9d0*(ast/pi)**2

c mod. by UE:
c The following contribution to dmt in getdmt() from the light Higgs/Goldstone
c is added as corrhig, using mt=tc/pol2 and ht2=mt2/vv**2/(1-cb**2):
c       dmty = mt/32d0/pi**2*(ht2*(1-cb**2)*(
c     .       2d0*myB1(mt2,mt2,0d0,mt2)+myB1(mt2,mb2,0d0,mt2)))
      mt2=(tc/pol2)**2
      mb2=bmass**2
      corrhig=1d0+mt2/32d0/pi**2/vv**2
     .      *(2d0*myB1(mt2,mt2,0d0,mt2)+myB1(mt2,mb2,0d0,mt2))

      corr=1d0+ast*g0/(4d0*pi*2d0*b0)*(-b1/b0+g1/g0)*(asx/ast-1d0)

      runt=rrr*corr/pol2*corrhig

      runt = runt*(1d0-asx/pi/3d0-(asx/pi)**2*55d0/144d0) ! shift to DRbar

      end

*
***********************************************************************
*

      DOUBLE PRECISION FUNCTION slasf(x,tc,asc,zm)

*     compute the running (SM, MSbar) alpha_s

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      pi = 4d0*atan(1d0)

      fn=5d0

      b0=11d0-2d0*fn/3d0
      b1=102d0-38d0*fn/3d0
      vvv=1d0-b0*asc/(2d0*pi)*log(zm/x)

      if(x.le.tc) then        ! five flavors

       slasf=asc/vvv*(1d0-b1/b0*asc/(4d0*pi*vvv)*log(vvv))

      else

       vvv=1d0-b0*asc/(2d0*pi)*log(zm/tc) ! first evolve up to q=mt

       ast=asc/vvv*(1d0-b1/b0*asc/(4d0*pi*vvv)*log(vvv))

       b0t=b0-2d0/3d0       ! six flavours
       b1t=b1-38d0/3d0
       vvv=1d0-b0t*ast/(2d0*pi)*log(tc/x) !     now evolve up to the scale >mt

       slasf=ast/vvv*(1d0-b1t/b0t*ast/(4d0*pi*vvv)*log(vvv))

      endif

      end

*
************************************************************************
*

      DOUBLE PRECISION FUNCTION runb(x,mbmb,mt,asmz,mz)

*     compute the running (SM, DRbar) bottom mass

      IMPLICIT NONE

      DOUBLE PRECISION x,mbmb,mt,asmz,mz
      DOUBLE PRECISION pi,b0,b1,g0,g1,asx,asb,slasf,ast,rrr,corr,mbmt
      INTEGER nf

      pi=4d0*atan(1d0)

      if(x.le.mt) then

       nf = 5             ! evolve from mb to x with nf=5
       b0 = 11d0-2d0*nf/3d0
       b1 = 102d0-38d0*nf/3d0

       g0 = 8d0
       g1 = 404d0/3d0-40d0*nf/9d0

       asx=slasf(x,mt,asmz,mz)
       asb=slasf(mbmb,mt,asmz,mz)

       rrr=mbmb*(asx/asb)**(g0/(2d0*b0))
       corr=1d0+asb*g0/(4d0*pi*2d0*b0)*(-b1/b0+g1/g0)*(asx/asb-1d0)

       runb=rrr*corr

      else

       nf = 5             ! first evolve from mb to mt with nf=5
       b0 = 11d0-2d0*nf/3d0
       b1 = 102d0-38d0*nf/3d0

       g0 = 8d0
       g1 = 404d0/3d0-40d0*nf/9d0

       ast=slasf(mt,mt,asmz,mz)
       asb=slasf(mbmb,mt,asmz,mz)

       rrr=mbmb*(ast/asb)**(g0/(2d0*b0))
       corr=1d0+asb*g0/(4d0*pi*2d0*b0)*(-b1/b0+g1/g0)*(ast/asb-1d0)

       mbmt=rrr*corr

       nf = 6             ! then evolve from mt to x with nf=6
       b0 = 11d0-2d0*nf/3d0
       b1 = 102d0-38d0*nf/3d0

       g0 = 8d0
       g1 = 404d0/3d0-40d0*nf/9d0

       asx=slasf(x,mt,asmz,mz)

       rrr=mbmt*(asx/ast)**(g0/(2d0*b0))
       corr=1d0+ast*g0/(4d0*pi*2d0*b0)*(-b1/b0+g1/g0)*(asx/ast-1d0)

       runb=rrr*corr

      endif

      runb = runb*(1d0-asx/pi/3d0-(asx/pi)**2*55d0/144d0) ! shift to DRbar

      end

*
***********************************************************************
*

      SUBROUTINE swap12(M)

      IMPLICIT NONE
      DOUBLE PRECISION M(3,3),temp

      temp = M(1,1)
      M(1,1) = M(2,2)
      M(2,2) = temp

      temp = M(1,3)
      M(1,3) = M(2,3)
      M(2,3) = temp

      temp = M(3,1)
      M(3,1) = M(3,2)
      M(3,2) = temp

      end

*
***********************************************************************
*

      SUBROUTINE jacobi(a,np,d,v)

*     just CALLs Tomas Hahn's diagonalization routine

      INTEGER np,i,j
      DOUBLE PRECISION a(np,np),d(np),v(np,np)
      DOUBLE COMPLEX M(np,np),U(np,np)

      do i=1,np              ! turn to complex
       do j=1,np
          M(i,j) = DCMPLX(a(i,j))
       enddo
      enddo


      CALL HEigensystem(np, M, np, d, U, np, 0)

      do i=1,np              ! back to real
       do j=1,np
          v(i,j) = DBLE(U(j,i)) ! the other jacobi had this convention
       enddo
      enddo

      end

* diagonalization of a Hermitian n-by-n matrix using the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 27 Sep 07 th
************************************************************************
** HEigensystem diagonalizes a Hermitian n-by-n matrix.
** Input: n, A = n-by-n matrix, Hermitian
** (only the upper triangle of A needs to be filled).
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^+ = U A U^-1 with U unitary.

      SUBROUTINE HEigensystem(n, A,ldA, d, U,ldU, sort)
      IMPLICIT NONE
      INTEGER n, ldA, ldU, sort
      DOUBLE COMPLEX A(ldA,*), U(ldU,*)
      DOUBLE PRECISION d(*)

      INTEGER p, q, j
      DOUBLE PRECISION red, off, thresh
      DOUBLE PRECISION delta, t, invc, s
      DOUBLE COMPLEX x, y, Apq
      DOUBLE PRECISION ev(2,16)

      INTEGER sweep
      COMMON /nsweeps/ sweep

      DOUBLE PRECISION sq
      DOUBLE COMPLEX c
      sq(c) = DBLE(c*DCONJG(c))

      if( n .gt. 16 ) then
        print *, "Dimension too large"
        d(1) = -999d0
        return
      endif

      do p = 1, n
        ev(1,p) = 0d0
        ev(2,p) = DBLE(A(p,p))
        d(p) = ev(2,p)
      enddo

      do p = 1, n
        do q = 1, n
          U(q,p) = 0
        enddo
        U(p,p) = 1
      enddo

      red = .04d0/n**4

      do sweep = 1, 50
        off = 0
        do q = 2, n
          do p = 1, q-1
            off = off+sq(A(p,q))
          enddo
        enddo
        if( off .lt. 2d0**(-103) ) goto 1

        thresh = 0
        if( sweep .lt. 4 ) thresh = off*red

        do q = 2, n
          do p = 1, q-1
            off = sq(A(p,q))
            if( sweep .gt. 4 .and. off .lt.
     .            2d0**(-103)*max(ev(2,p)**2, ev(2,q)**2) ) then
            A(p,q) = 0
            else
            if( off .gt. thresh ) then
              t = .5d0*(ev(2,p)-ev(2,q))
              t = 1d0/(t+sign(sqrt(t**2+off), t))

              delta = t*off
              ev(1,p) = ev(1,p)+delta
              ev(2,p) = d(p)+ev(1,p)
              ev(1,q) = ev(1,q)-delta
              ev(2,q) = d(q)+ev(1,q)

              invc = sqrt(delta*t+1d0)
              s = t/invc
              t = delta/(invc+1d0)

              Apq = A(p,q)

              do j = 1, p-1
                x = A(j,p)
                y = A(j,q)
                A(j,p) = x+s*(DCONJG(Apq)*y-t*x)
                A(j,q) = y-s*(Apq*x+t*y)
              enddo

              do j = p+1, q-1
                x = A(p,j)
                y = A(j,q)
                A(p,j) = x+s*(Apq*DCONJG(y)-t*x)
                A(j,q) = y-s*(Apq*DCONJG(x)+t*y)
              enddo

              do j = q+1, n
                x = A(p,j)
                y = A(q,j)
                A(p,j) = x+s*(Apq*y-t*x)
                A(q,j) = y-s*(DCONJG(Apq)*x+t*y)
              enddo

              A(p,q) = 0

              do j = 1, n
                x = U(p,j)
                y = U(q,j)
                U(p,j) = x+s*(Apq*y-t*x)
                U(q,j) = y-s*(DCONJG(Apq)*x+t*y)
              enddo
            endif
            endif
          enddo
        enddo

        do p = 1, n
          ev(1,p) = 0
          d(p) = ev(2,p)
        enddo
      enddo

      print *, "Bad convergence in HEigensystem"

1     if( sort .eq. 0 ) return

* sort the eigenvalues

      do p = 1, n-1
        j = p
        t = d(p)
        do q = p+1, n
          if( sort*(t-d(q)) .gt. 0 ) then
            j = q
            t = d(q)
          endif
        enddo

        if( j .ne. p ) then
          d(j) = d(p)
          d(p) = t
          do q = 1, n
            x = U(p,q)
            U(p,q) = U(j,q)
            U(j,q) = x
          enddo
        endif
      enddo
      end

*
***********************************************************************
*

      SUBROUTINE gettadS(g,gp,ht,hb,htau,v1,v2,Q,tadS)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ht,hb,htau,v1,v2,Q,tadS(3)

      DOUBLE PRECISION mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     .     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),
     .     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     .     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2),myA0,pi,
     .     mw2,mz2,mt,mb,mtau,gb2,sq2

      COMMON/coupl_s/lsstt,lssbb,lsstata,lssntnt,lssuu,lssdd,
     .     lssee,lssnn,lstt,lsbb,lstata,lsntnt,lsuu,lsdd,lsee,lsnn,
     .     lsshh,lssaa,lshh,lsaa,lsscc,lscc,lsnene,lschch

      DOUBLE PRECISION lsstt(3,3,2,2),lssbb(3,3,2,2),lsstata(3,3,2,2),
     .     lssntnt(3,3),lssuu(3,3,2,2),lssdd(3,3,2,2),lssee(3,3,2,2),
     .     lssnn(3,3),lstt(3,2,2),lsbb(3,2,2),lstata(3,2,2),lsntnt(3),
     .     lsuu(3,2,2),lsdd(3,2,2),lsee(3,2,2),lsnn(3),
     .     lsshh(3,3,3,3),lssaa(3,3,3,3),lshh(3,3,3),lsaa(3,3,3),
     .     lsscc(3,3,2,2),lscc(3,2,2),lschch(3,2,2),lsnene(3,5,5)

      INTEGER i,k

      COMMON/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     .     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue

      pi=4d0*atan(1d0)
      sq2=sqrt(2d0)

*     the gauge and fermion contributions

      gb2=(g**2+gp**2)/2d0

      mw2=g**2/2d0*(v1**2+v2**2)
      mz2=gb2*(v1**2+v2**2)

      mt=ht*v2
      mb=hb*v1
      mtau=htau*v1

      tadS(1)=0d0
     .     -6d0*sq2*hb*mb*myA0(mb**2,Q**2)
     .     -2d0*sq2*htau*mtau*myA0(mtau**2,Q**2)
     .     +3d0*v1/sq2*g**2*myA0(mw2,Q**2)+3d0*v1/sq2*gb2*myA0(mz2,Q**2)

      tadS(2)=0d0
     .     -6d0*sq2*ht*mt*myA0(mt**2,Q**2)
     .     +3d0*v2/sq2*g**2*myA0(mw2,Q**2)+3d0*v2/sq2*gb2*myA0(mz2,Q**2)

      tadS(3)=0d0

*     the sfermion contributions

*      CALL coupl_s_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
*     .     Rt,Rb,Rtau,lsstt,lssbb,lsstata,lssntnt,lssuu,lssdd,
*     .     lssee,lssnn,lstt,lsbb,lstata,lsntnt,lsuu,lsdd,lsee,lsnn)

      do i=1,3

       do k=1,2

          tadS(i)=tadS(i)
     .        +3d0*lstt(i,k,k)*myA0(mstop(k)**2,Q**2)
     .        +3d0*lsbb(i,k,k)*myA0(msbot(k)**2,Q**2)
     .        +lstata(i,k,k)*myA0(mstau(k)**2,Q**2)
     .        +6d0*lsuu(i,k,k)*myA0(msup(k)**2,Q**2)
     .        +6d0*lsdd(i,k,k)*myA0(msdown(k)**2,Q**2)
     .        +2d0*lsee(i,k,k)*myA0(msel(k)**2,Q**2)

       enddo

       tadS(i)=tadS(i)      ! add sneutrinos (no sum)
     .       +2d0*lsnn(i)*myA0(msnue**2,Q**2)
     .       +lsntnt(i)*myA0(msnutau**2,Q**2)

      enddo

*     the Higgs contributions

*      CALL coupl_s_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
*     .     lsshh,lssaa,lshh,lsaa,lsscc,lscc)

      do i=1,3

       do k=1,3           ! neutral Higgses

          tadS(i)=tadS(i)
     .        +lshh(i,k,k)*myA0(mhh(k)**2,Q**2)
     .        +lsaa(i,k,k)*myA0(maa(k)**2,Q**2)

       enddo

       do k=1,2           ! charged Higgses

          tadS(i)=tadS(i)
     .        +lscc(i,k,k)*myA0(mhc(k)**2,Q**2)

       enddo

      enddo

*     the chargino and neutralino contributions

*      CALL coupl_s_ino(g,gp,ll,kk,NN,UU,VV,lsnene,lschch)

      do i=1,3

       do k=1,5           ! neutralinos
          tadS(i)=tadS(i)
     .        -4d0*lsnene(i,k,k)*mne(k)*myA0(mne(k)**2,Q**2)
       enddo

       do k=1,2           ! charginos
          tadS(i)=tadS(i)
     .        -4d0*lschch(i,k,k)*mch(k)*myA0(mch(k)**2,Q**2)
       enddo

      enddo

      do i=1,3
       tadS(i)=tadS(i)/16d0/pi**2
      enddo

      end

*
***********************************************************************
*

      SUBROUTINE getPiSS(g,gp,ht,hb,htau,v1,v2,p,Q,piSS)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ht,hb,htau,v1,v2,p,Q,piSS(3,3)

      DOUBLE PRECISION mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     .     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),
     .     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     .     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2),
     .     myB0,myA0,myF,myG,pi,cb,sb,mw2,mz2,mt2,mb2,mtau2,gb2

      COMMON/coupl_s/lsstt,lssbb,lsstata,lssntnt,lssuu,lssdd,
     .     lssee,lssnn,lstt,lsbb,lstata,lsntnt,lsuu,lsdd,lsee,lsnn,
     .     lsshh,lssaa,lshh,lsaa,lsscc,lscc,lsnene,lschch

      DOUBLE PRECISION lsstt(3,3,2,2),lssbb(3,3,2,2),lsstata(3,3,2,2),
     .     lssntnt(3,3),lssuu(3,3,2,2),lssdd(3,3,2,2),lssee(3,3,2,2),
     .     lssnn(3,3),lstt(3,2,2),lsbb(3,2,2),lstata(3,2,2),lsntnt(3),
     .     lsuu(3,2,2),lsdd(3,2,2),lsee(3,2,2),lsnn(3),
     .     lsshh(3,3,3,3),lssaa(3,3,3,3),lshh(3,3,3),lsaa(3,3,3),
     .     lsscc(3,3,2,2),lscc(3,2,2),lschch(3,2,2),lsnene(3,5,5)

      INTEGER i,j,k,l

      COMMON/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     .     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue

      pi=4d0*atan(1d0)

      do i=1,3              ! initialize
       do j=1,3
          PiSS(i,j)=0d0
       enddo
      enddo

*     the gauge and fermion contributions

      cb=v1/sqrt(v1**2+v2**2)
      sb=v2/sqrt(v1**2+v2**2)
      gb2=(g**2+gp**2)/2d0

      mw2=g**2/2d0*(v1**2+v2**2)
      mz2=gb2*(v1**2+v2**2)

      mt2=ht**2*v2**2
      mb2=hb**2*v1**2
      mtau2=htau**2*v1**2

      PiSS(1,1) =
     .     3d0*hb**2*((p**2-4d0*mb2)*myB0(p**2,mb2,mb2,Q**2)
     .     -2d0*myA0(mb2,Q**2))
     .     +htau**2*((p**2-4d0*mtau2)*myB0(p**2,mtau2,mtau2,Q**2)
     .     -2d0*myA0(mtau2,Q**2))
     .     +7d0/2d0*cb**2*(g**2*mw2*myB0(p**2,mw2,mw2,Q**2)
     .     +gb2*mz2*myB0(p**2,mz2,mz2,Q**2))
     .     +2d0*g**2*myA0(mw2,Q**2)+2d0*gb2*myA0(mz2,Q**2)

      PiSS(1,2) =
     .     7d0/2d0*cb*sb*(g**2*mw2*myB0(p**2,mw2,mw2,Q**2)
     .     +gb2*mz2*myB0(p**2,mz2,mz2,Q**2))

      PiSS(2,1)=PiSS(1,2)

      PiSS(2,2)=3d0*ht**2*((p**2-4d0*mt2)*myB0(p**2,mt2,mt2,Q**2)
     .     -2d0*myA0(mt2,Q**2))
     .     +7d0/2d0*sb**2*(g**2*mw2*myB0(p**2,mw2,mw2,Q**2)
     .     +gb2*mz2*myB0(p**2,mz2,mz2,Q**2))
     .     +2d0*g**2*myA0(mw2,Q**2)+2d0*gb2*myA0(mz2,Q**2)

*     pseudoscalar-gauge contribution

      do i=1,2
       do j=1,2

          do k=1,3
             PiSS(i,j)=PiSS(i,j)+(-1d0)**(i+j)*
     .            gb2/2d0*RP(k,i)*RP(k,j)*myF(p**2,maa(k)**2,mz2,Q**2)
          enddo

          do k=1,2
             PiSS(i,j)=PiSS(i,j)+(-1d0)**(i+j)*
     .           g**2/2d0*RC(k,i)*RC(k,j)*myF(p**2,mhc(k)**2,mw2,Q**2)
          enddo

       enddo
      enddo

*     the sfermion contributions

*      CALL coupl_s_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
*     .     Rt,Rb,Rtau,lsstt,lssbb,lsstata,lssntnt,lssuu,lssdd,
*     .     lssee,lssnn,lstt,lsbb,lstata,lsntnt,lsuu,lsdd,lsee,lsnn)


      do i=1,3
       do j=1,3
          do k=1,2

             PiSS(i,j)=PiSS(i,j)
     .           +6d0*lsstt(i,j,k,k)*myA0(mstop(k)**2,Q**2)
     .           +6d0*lssbb(i,j,k,k)*myA0(msbot(k)**2,Q**2)
     .           +2d0*lsstata(i,j,k,k)*myA0(mstau(k)**2,Q**2)
     .           +12d0*lssuu(i,j,k,k)*myA0(msup(k)**2,Q**2)
     .           +12d0*lssdd(i,j,k,k)*myA0(msdown(k)**2,Q**2)
     .           +4d0*lssee(i,j,k,k)*myA0(msel(k)**2,Q**2)

             do l=1,2

              PiSS(i,j)=PiSS(i,j)
     .            +3d0*lstt(i,k,l)*lstt(j,l,k)
     .             *myB0(p**2,mstop(k)**2,mstop(l)**2,Q**2)
     .            +3d0*lsbb(i,k,l)*lsbb(j,l,k)
     .             *myB0(p**2,msbot(k)**2,msbot(l)**2,Q**2)
     .            +lstata(i,k,l)*lstata(j,l,k)
     .             *myB0(p**2,mstau(k)**2,mstau(l)**2,Q**2)
     .            +6d0*lsuu(i,k,l)*lsuu(j,l,k)
     .             *myB0(p**2,msup(k)**2,msup(l)**2,Q**2)
     .            +6d0*lsdd(i,k,l)*lsdd(j,l,k)
     .             *myB0(p**2,msdown(k)**2,msdown(l)**2,Q**2)
     .            +2d0*lsee(i,k,l)*lsee(j,l,k)
     .             *myB0(p**2,msel(k)**2,msel(l)**2,Q**2)
             enddo
          enddo

          PiSS(i,j)=PiSS(i,j) ! add sneutrinos (no sum)
     .        +4d0*lssnn(i,j)*myA0(msnue**2,Q**2)
     .        +2d0*lsnn(i)*lsnn(j)
     .         *myB0(p**2,msnue**2,msnue**2,Q**2)
     .        +2d0*lssntnt(i,j)*myA0(msnutau**2,Q**2)
     .        +lsntnt(i)*lsntnt(j)
     .         *myB0(p**2,msnutau**2,msnutau**2,Q**2)

       enddo
      enddo

*     the Higgs contributions

*      CALL coupl_s_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
*     .     lsshh,lssaa,lshh,lsaa,lsscc,lscc)

      do i=1,3
       do j=1,3

          do k=1,3        ! neutral Higgses

             PiSS(i,j)=PiSS(i,j)
     .           +2d0*lsshh(i,j,k,k)*myA0(mhh(k)**2,Q**2)
     .           +2d0*lssaa(i,j,k,k)*myA0(maa(k)**2,Q**2)

             do l=1,3

              PiSS(i,j)=PiSS(i,j)
     .            +2d0*lshh(i,k,l)*lshh(j,k,l)
     .             *myB0(p**2,mhh(k)**2,mhh(l)**2,Q**2)
     .            +2d0*lsaa(i,k,l)*lsaa(j,k,l)
     .             *myB0(p**2,maa(k)**2,maa(l)**2,Q**2)
             enddo
          enddo

          do k=1,2        ! charged Higgses

             PiSS(i,j)=PiSS(i,j)
     .           +2d0*lsscc(i,j,k,k)*myA0(mhc(k)**2,Q**2)

             do l=1,2

              PiSS(i,j)=PiSS(i,j)
     .            +lscc(i,k,l)*lscc(j,l,k)
     .             *myB0(p**2,mhc(k)**2,mhc(l)**2,Q**2)

             enddo

          enddo
       enddo
      enddo

*     the chargino and neutralino contributions

*     CALL coupl_s_ino(g,gp,ll,kk,NN,UU,VV,lsnene,lschch)

      do i=1,3
       do j=1,3

          do k=1,5        ! neutralinos
             do l=1,5

              PiSS(i,j)=PiSS(i,j)
     .            +4d0*lsnene(i,k,l)*lsnene(j,k,l)*(
     .             myG(p**2,mne(k)**2,mne(l)**2,Q**2)
     .            -2d0*mne(k)*mne(l)
     .             *myB0(p**2,mne(k)**2,mne(l)**2,Q**2))
             enddo
          enddo
          do k=1,2        ! charginos
             do l=1,2

              PiSS(i,j)=PiSS(i,j)
     .            +2d0*(lschch(i,k,l)*lschch(j,k,l)*
     .             myG(p**2,mch(k)**2,mch(l)**2,Q**2)
     .            -2d0*lschch(i,k,l)*lschch(j,l,k)*mch(k)*mch(l)
     .             *myB0(p**2,mch(k)**2,mch(l)**2,Q**2))
             enddo
          enddo

       enddo
      enddo

      do i=1,3
       do j=1,3
          PiSS(i,j)=PiSS(i,j)/16d0/pi**2
       enddo
      enddo

      end

*
***********************************************************************
*

      SUBROUTINE getPiPP(g,gp,ht,hb,htau,v1,v2,p,Q,piPP)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ht,hb,htau,v1,v2,p,Q,piPP(3,3)

      DOUBLE PRECISION mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     .     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),
     .     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     .     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2),
     .     myB0,myA0,myF,myG,pi,cb,sb,mw2,mz2,mt2,mb2,mtau2,gb2,ghost

      COMMON/coupl_p/lpptt,lppbb,lpptata,lppntnt,lppuu,lppdd,
     .     lppee,lppnn,lptt,lpbb,lptata,lpntnt,lpuu,lpdd,lpee,lpnn,
     .     lpphh,lppaa,lpah,lppcc,lpcc,lpnene,lpchch

      DOUBLE PRECISION lpptt(3,3,2,2),lppbb(3,3,2,2),lpptata(3,3,2,2),
     .     lppntnt(3,3),lppuu(3,3,2,2),lppdd(3,3,2,2),lppee(3,3,2,2),
     .     lppnn(3,3),lptt(3,2,2),lpbb(3,2,2),lptata(3,2,2),lpntnt(3),
     .     lpuu(3,2,2),lpdd(3,2,2),lpee(3,2,2),lpnn(3),
     .     lpphh(3,3,3,3),lppaa(3,3,3,3),lpah(3,3,3),
     .     lppcc(3,3,2,2),lpcc(3,2,2),lpchch(3,2,2),lpnene(3,5,5)

      INTEGER i,j,k,l

      COMMON/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     .     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue

      pi=4d0*atan(1d0)

      do i=1,3              ! initialize
       do j=1,3
          PiPP(i,j)=0d0
       enddo
      enddo

*     the gauge and fermion contributions

      cb=v1/sqrt(v1**2+v2**2)
      sb=v2/sqrt(v1**2+v2**2)
      gb2=(g**2+gp**2)/2d0

      mw2=g**2/2d0*(v1**2+v2**2)
      mz2=gb2*(v1**2+v2**2)

      mt2=ht**2*v2**2
      mb2=hb**2*v1**2
      mtau2=htau**2*v1**2

      ghost=g**2/2d0*mw2*myB0(p**2,mw2,mw2,Q**2)

      PiPP(1,1) =
     .     3d0*hb**2*(p**2*myB0(p**2,mb2,mb2,Q**2)-2d0*myA0(mb2,Q**2))
     .     +htau**2*(p**2*myB0(p**2,mtau2,mtau2,Q**2)
     .     -2d0*myA0(mtau2,Q**2))
     .     +2d0*g**2*myA0(mw2,Q**2)+2d0*gb2*myA0(mz2,Q**2)
     .     +cb**2*ghost

      PiPP(1,2) =
     .     -cb*sb*ghost

      PiPP(2,1)=PiPP(1,2)

      PiPP(2,2) =
     .     3d0*ht**2*(p**2*myB0(p**2,mt2,mt2,Q**2)-2d0*myA0(mt2,Q**2))
     .     +2d0*g**2*myA0(mw2,Q**2)+2d0*gb2*myA0(mz2,Q**2)
     .     +sb**2*ghost

*     scalar-Z contribution

      do i=1,2
       do j=1,2

          do k=1,3
             PiPP(i,j)=PiPP(i,j)+(-1d0)**(i+j)*
     .            gb2/2d0*RS(k,i)*RS(k,j)*myF(p**2,mhh(k)**2,mz2,Q**2)
          enddo

          do k=1,2
             PiPP(i,j)=PiPP(i,j) +
     .           g**2/2d0*RC(k,i)*RC(k,j)*myF(p**2,mhc(k)**2,mw2,Q**2)
          enddo

       enddo
      enddo

*     the sfermion contributions

*      CALL coupl_p_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
*     .     Rt,Rb,Rtau,lpptt,lppbb,lpptata,lppntnt,lppuu,lppdd,
*     .     lppee,lppnn,lptt,lpbb,lptata,lpntnt,lpuu,lpdd,lpee,lpnn)


      do i=1,3
       do j=1,3
          do k=1,2

             PiPP(i,j)=PiPP(i,j)
     .           +6d0*lpptt(i,j,k,k)*myA0(mstop(k)**2,Q**2)
     .           +6d0*lppbb(i,j,k,k)*myA0(msbot(k)**2,Q**2)
     .           +2d0*lpptata(i,j,k,k)*myA0(mstau(k)**2,Q**2)
     .           +12d0*lppuu(i,j,k,k)*myA0(msup(k)**2,Q**2)
     .           +12d0*lppdd(i,j,k,k)*myA0(msdown(k)**2,Q**2)
     .           +4d0*lppee(i,j,k,k)*myA0(msel(k)**2,Q**2)

             do l=1,2

              PiPP(i,j)=PiPP(i,j)
     .            -3d0*lptt(i,k,l)*lptt(j,l,k)
     .             *myB0(p**2,mstop(k)**2,mstop(l)**2,Q**2)
     .            -3d0*lpbb(i,k,l)*lpbb(j,l,k)
     .             *myB0(p**2,msbot(k)**2,msbot(l)**2,Q**2)
     .            -lptata(i,k,l)*lptata(j,l,k)
     .             *myB0(p**2,mstau(k)**2,mstau(l)**2,Q**2)
     .            -6d0*lpuu(i,k,l)*lpuu(j,l,k)
     .             *myB0(p**2,msup(k)**2,msup(l)**2,Q**2)
     .            -6d0*lpdd(i,k,l)*lpdd(j,l,k)
     .             *myB0(p**2,msdown(k)**2,msdown(l)**2,Q**2)
     .            -2d0*lpee(i,k,l)*lpee(j,l,k)
     .             *myB0(p**2,msel(k)**2,msel(l)**2,Q**2)
             enddo
          enddo

          PiPP(i,j)=PiPP(i,j) ! add sneutrinos (no sum)
     .        +4d0*lppnn(i,j)*myA0(msnue**2,Q**2)
     .        -2d0*lpnn(i)*lpnn(j)
     .         *myB0(p**2,msnue**2,msnue**2,Q**2)
     .        +2d0*lppntnt(i,j)*myA0(msnutau**2,Q**2)
     .        -lpntnt(i)*lpntnt(j)
     .         *myB0(p**2,msnutau**2,msnutau**2,Q**2)

       enddo
      enddo

*     the Higgs contributions

*      CALL coupl_p_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
*     .     lpphh,lppaa,lpah,lppcc,lpcc)

      do i=1,3
       do j=1,3

          do k=1,3       ! neutral Higgses

             PiPP(i,j)=PiPP(i,j)
     .           +2d0*lpphh(i,j,k,k)*myA0(mhh(k)**2,Q**2)
     .           +2d0*lppaa(i,j,k,k)*myA0(maa(k)**2,Q**2)

             do l=1,3

              PiPP(i,j)=PiPP(i,j)
     .            +lpah(i,k,l)*lpah(j,k,l)
     .             *myB0(p**2,maa(k)**2,mhh(l)**2,Q**2)

             enddo
          enddo

          do k=1,2        ! charged Higgses

             PiPP(i,j)=PiPP(i,j)
     .           +2d0*lppcc(i,j,k,k)*myA0(mhc(k)**2,Q**2)

             do l=1,2

              PiPP(i,j)=PiPP(i,j)
     .            -lpcc(i,k,l)*lpcc(j,l,k)
     .             *myB0(p**2,mhc(k)**2,mhc(l)**2,Q**2)

             enddo

          enddo
       enddo
      enddo

*     the chargino and neutralino contributions

*      CALL coupl_p_ino(g,gp,ll,kk,NN,UU,VV,lpnene,lpchch)

      do i=1,3
       do j=1,3

          do k=1,5        ! neutralinos
             do l=1,5

              PiPP(i,j)=PiPP(i,j)
     .            +4d0*lpnene(i,k,l)*lpnene(j,k,l)*(
     .             myG(p**2,mne(k)**2,mne(l)**2,Q**2)
     .             +2d0*mne(k)*mne(l)
     .             *myB0(p**2,mne(k)**2,mne(l)**2,Q**2))
             enddo
          enddo

          do k=1,2        ! charginos
             do l=1,2

              PiPP(i,j)=PiPP(i,j)
     .            +2d0*(lpchch(i,k,l)*lpchch(j,k,l)*
     .             myG(p**2,mch(k)**2,mch(l)**2,Q**2)
     .            +2d0*lpchch(i,k,l)*lpchch(j,l,k)*mch(k)*mch(l)
     .             *myB0(p**2,mch(k)**2,mch(l)**2,Q**2))
             enddo
          enddo

       enddo
      enddo

      do i=1,3
       do j=1,3
          PiPP(i,j)=PiPP(i,j)/16d0/pi**2
       enddo
      enddo

      end

*
***********************************************************************
*

      SUBROUTINE getPiZZ(g,gp,ht,hb,htau,v1,v2,p,Q,piZZ)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ht,hb,htau,v1,v2,p,Q,piZZ

      DOUBLE PRECISION mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     .     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),del(2,2),
     .     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     .     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2),
     .     myB0,myB22T,myH,pi,cb,sb,mw2,mz2,mt2,mb2,mtau2,gb2,
     .     sw2,cw2,guL,guR,gdL,gdR,geL,geR,gnu,lznene(5,5),
     .     azchch(2,2),bzchch(2,2)

      INTEGER i,j

      COMMON/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     .     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue

      pi=4d0*atan(1d0)

      PiZZ=0d0

      cb=v1/sqrt(v1**2+v2**2)
      sb=v2/sqrt(v1**2+v2**2)
      gb2=(g**2+gp**2)/2d0
      sw2=gp**2/2d0/gb2
      cw2=1d0-sw2

      mw2=g**2/2d0*(v1**2+v2**2)
      mz2=gb2*(v1**2+v2**2)

*     the higgs and gauge contributions

      do i=1,3            ! neutral scalars

       PiZZ=PiZZ          ! scalar-Z
     .      +2d0*gb2*mz2*(cb*RS(i,1)+sb*RS(i,2))**2
     .      *myB0(p**2,mhh(i)**2,mz2,Q**2)

       do j=1,3

          PiZZ=PiZZ
     .        -2d0*gb2*(RS(i,1)*RP(j,1)-RS(i,2)*RP(j,2))**2
     .         *myB22T(p**2,mhh(i)**2,maa(j)**2,Q**2)

       enddo
      enddo

      do i=1,2            ! charged scalars

       PiZZ=PiZZ
     .      -2d0*gb2*(cw2-sw2)**2
     .      *myB22T(p**2,mhc(i)**2,mhc(i)**2,Q**2)

      enddo

      PiZZ=PiZZ             ! pure gauge
     .     -4d0*gb2*cw2**2*(2d0*p**2+mw2-mz2*sw2**2/cw2)
     .     *myB0(p**2,mw2,mw2,Q**2)
     .     -16d0*gb2*cw2**2*myB22T(p**2,mw2,mw2,Q**2)

*     the fermion contributions

      guL= .5d0-2d0/3d0*sw2
      gdL=-.5d0+1d0/3d0*sw2
      gnu= .5d0
      geL=-.5d0+sw2
      guR= 2d0/3d0*sw2
      gdR=-1d0/3d0*sw2
      geR=-sw2

      mt2=ht**2*v2**2
      mb2=hb**2*v1**2
      mtau2=htau**2*v1**2

      PiZZ=PiZZ+2d0*gb2*(
     .     (6d0*(guL**2+guR**2)+6d0*(gdL**2+gdR**2)+2d0*(geL**2+geR**2)
     .     +3d0*gnu**2)*myH(p**2,0d0,0d0,Q**2)
     .     +3d0*((guL**2+guR**2)*myH(p**2,mt2,mt2,Q**2)
     .     -4d0*guL*guR*mt2*myB0(p**2,mt2,mt2,Q**2))
     .     +3d0*((gdL**2+gdR**2)*myH(p**2,mb2,mb2,Q**2)
     .     -4d0*gdL*gdR*mb2*myB0(p**2,mb2,mb2,Q**2))
     .     +((geL**2+geR**2)*myH(p**2,mtau2,mtau2,Q**2)
     .     -4d0*geL*geR*mtau2*myB0(p**2,mtau2,mtau2,Q**2)))

*     the sfermion contributions

      del(1,1)=1d0
      del(1,2)=0d0
      del(2,1)=0d0
      del(2,2)=1d0

      do i=1,2
       do j=1,2

          PiZZ=PiZZ-8d0*gb2*(
     .         6d0*(guL*del(i,1)*del(j,1)-guR*del(i,2)*del(j,2))**2
     .         *myB22T(p**2,msup(i)**2,msup(j)**2,Q**2)
     .         +6d0*(gdL*del(i,1)*del(j,1)-gdR*del(i,2)*del(j,2))**2
     .         *myB22T(p**2,msdown(i)**2,msdown(j)**2,Q**2)
     .         +2d0*(geL*del(i,1)*del(j,1)-geR*del(i,2)*del(j,2))**2
     .         *myB22T(p**2,msel(i)**2,msel(j)**2,Q**2)
     .         +3d0*(guL*Rt(i,1)*Rt(j,1)-guR*Rt(i,2)*Rt(j,2))**2
     .         *myB22T(p**2,mstop(i)**2,mstop(j)**2,Q**2)
     .         +3d0*(gdL*Rb(i,1)*Rb(j,1)-gdR*Rb(i,2)*Rb(j,2))**2
     .         *myB22T(p**2,msbot(i)**2,msbot(j)**2,Q**2)
     .         +(geL*Rtau(i,1)*Rtau(j,1)-geR*Rtau(i,2)*Rtau(j,2))**2
     .         *myB22T(p**2,mstau(i)**2,mstau(j)**2,Q**2))

       enddo
      enddo

      PiZZ=PiZZ-8d0*gb2*gnu**2*(
     .     2d0*myB22T(p**2,msnue**2,msnue**2,Q**2)
     .     +myB22T(p**2,msnutau**2,msnutau**2,Q**2))

*     the chargino and neutralino contributions


      CALL coupl_Z_ino(g,gp,NN,UU,VV,lznene,azchch,bzchch)

      do i=1,5            ! neutralinos
       do j=1,5

          PiZZ=PiZZ
     .         +4d0*lznene(i,j)**2*(myH(p**2,mne(i)**2,mne(j)**2,Q**2)
     .         -2d0*mne(i)*mne(j)*myB0(p**2,mne(i)**2,mne(j)**2,Q**2))

       enddo
      enddo

      do i=1,2            ! charginos
       do j=1,2

          PiZZ=PiZZ
     .         +(azchch(i,j)**2+bzchch(i,j)**2)
     .         *myH(p**2,mch(i)**2,mch(j)**2,Q**2)
     .         +4d0*azchch(i,j)*bzchch(i,j)*
     .         mch(i)*mch(j)*myB0(p**2,mch(i)**2,mch(j)**2,Q**2)

       enddo
      enddo

      PiZZ=PiZZ/16d0/pi**2

      end

*
***********************************************************************
*

      SUBROUTINE getPiWW(g,gp,ht,hb,htau,v1,v2,p,Q,piWW)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ht,hb,htau,v1,v2,p,Q,piWW

      DOUBLE PRECISION mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     .     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),del(2,2),
     .     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     .     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2),
     .     myB0,myB22T,myH,pi,cb,sb,mw2,mz2,mt2,mb2,mtau2,gb2,
     .     sw2,cw2,awnech(5,2),bwnech(5,2)

      INTEGER i,j

      COMMON/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     .     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue

      pi=4d0*atan(1d0)

      PiWW=0d0

      cb=v1/sqrt(v1**2+v2**2)
      sb=v2/sqrt(v1**2+v2**2)
      gb2=(g**2+gp**2)/2d0
      sw2=gp**2/2d0/gb2
      cw2=1d0-sw2

      mw2=g**2/2d0*(v1**2+v2**2)
      mz2=gb2*(v1**2+v2**2)

*     the higgs and gauge contributions

      do i=1,3

       PiWW=PiWW          ! scalar-W
     .       +g**2*mw2*(cb*RS(i,1)+sb*RS(i,2))**2
     .      *myB0(p**2,mhh(i)**2,mw2,Q**2)

       do j=1,2

          PiWW=PiWW       ! scalar-charged
     .        -g**2*(RS(i,1)*RC(j,1)-RS(i,2)*RC(j,2))**2
     .         *myB22T(p**2,mhh(i)**2,mhc(j)**2,Q**2)

          PiWW=PiWW       ! pseudo-charged
     .        -g**2*(RP(i,1)*RC(j,1)+RP(i,2)*RC(j,2))**2
     .         *myB22T(p**2,maa(i)**2,mhc(j)**2,Q**2)

       enddo
      enddo

      PiWW=PiWW             ! pure gauge
     .     -8d0*g**2*cw2*myB22T(p**2,mw2,mz2,Q**2)
     .     -g**2*((4d0*p**2+mw2+mz2)*cw2-mz2*sw2**2)
     .     *myB0(p**2,mz2,mw2,Q**2)
     .     -sw2*g**2*(8d0*myB22T(p**2,mw2,0d0,Q**2)
     .     +4d0*p**2*myB0(p**2,mw2,0d0,Q**2))

*     the fermion contributions

      mt2=ht**2*v2**2
      mb2=hb**2*v1**2
      mtau2=htau**2*v1**2

      PiWW=PiWW+g**2/2d0*(
     .     8d0*myH(p**2,0d0,0d0,Q**2)
     .     +3d0*myH(p**2,mt2,mb2,Q**2)
     .     +myH(p**2,mtau2,0d0,Q**2))

*     the sfermion contributions

      del(1,1)=1d0
      del(1,2)=0d0
      del(2,1)=0d0
      del(2,2)=1d0

      do i=1,2
       do j=1,2

          PiWW=PiWW-2d0*g**2*(
     .         6d0*(del(i,1)*del(j,1))**2
     .         *myB22T(p**2,msup(i)**2,msdown(j)**2,Q**2)
     .         +3d0*(Rt(i,1)*Rb(j,1))**2
     .         *myB22T(p**2,mstop(i)**2,msbot(j)**2,Q**2))

       enddo

       PiWW=PiWW-2d0*g**2*(
     .      2d0*del(i,1)**2*myB22T(p**2,msel(i)**2,msnue**2,Q**2)
     .      +Rtau(i,1)**2*myB22T(p**2,mstau(i)**2,msnutau**2,Q**2))

      enddo

*     the chargino/neutralino contribution


      CALL coupl_W_ino(g,NN,UU,VV,awnech,bwnech)

      do i=1,5
       do j=1,2

          PiWW=PiWW
     .         +(awnech(i,j)**2+bwnech(i,j)**2)
     .         *myH(p**2,mne(i)**2,mch(j)**2,Q**2)
     .         +4d0*awnech(i,j)*bwnech(i,j)*
     .         mne(i)*mch(j)*myB0(p**2,mne(i)**2,mch(j)**2,Q**2)

       enddo
      enddo

      PiWW=PiWW/16d0/pi**2

      end

*
***********************************************************************
*

      SUBROUTINE treemasses(g,gp,ll,kk,ht,hb,htau,v1,v2,xx,M1,M2,
     .     Ak,Al,At,Ab,Atau,mQ3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,mer,
     .     Q,errmass)

      DOUBLE PRECISION g,gp,ll,kk,ht,hb,htau,v1,v2,xx,M1,M2,
     .     Ak,Al,At,Ab,Atau,mQ3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,mer,Q

      INTEGER errmass,errhiggs

      DOUBLE PRECISION mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     .     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),
     .     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     .     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2)

      LOGICAL errsfer

      COMMON/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     .     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue

      COMMON/FOREFFPOT/mt,T1,T2,st,ct,Q2
      DOUBLE PRECISION mt,T1,T2,st,ct,Q2

      errmass=0

*     compute all the tree-level masses


      CALL tree_sfermions(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     .     mQ3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,mer,mstop,msbot,
     .     mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue,errsfer)

      IF(errsfer) THEN
       errmass=3
       RETURN
      ENDIF

      mt=ht*v2
      T1=mstop(1)**2        ! for the Higgs mass calculation
      T2=mstop(2)**2
      st=Rt(1,2)            ! always true???
      ct=Rt(1,1)
      Q2=Q**2

      errhiggs=0


      CALL tree_higgses(g,gp,ll,kk,v1,v2,xx,Ak,Al,mhh,maa,mhc,RS,RP,RC,
     .     errhiggs)

      IF(errhiggs.ne.0) THEN
       errmass=errhiggs
       RETURN
      ENDIF


      CALL tree_charginos(g,ll,v1,v2,M2,xx,mch,UU,VV)


      CALL tree_neutralinos(g,gp,ll,kk,v1,v2,xx,M1,M2,mne,NN)

*      write(*,*) 'scalar',mhh
*      write(*,*) RS(1,1),RS(1,2),RS(1,3)
*      write(*,*) RS(2,1),RS(2,2),RS(2,3)
*      write(*,*) RS(3,1),RS(3,2),RS(3,3)
*
*      write(*,*) 'pseudo',maa
*      write(*,*) RP(1,1),RP(1,2),RP(1,3)
*      write(*,*) RP(2,1),RP(2,2),RP(2,3)
*      write(*,*) RP(3,1),RP(3,2),RP(3,3)
*
*      write(*,*) 'charged',mhc
*
*      write(*,*) 'sup',msup,mstop
*      write(*,*) Rt(1,1),Rt(1,2)
*      write(*,*) Rt(2,1),Rt(2,2)
*
*      write(*,*) 'sdown',msdown,msbot
*      write(*,*) Rb(1,1),Rb(1,2)
*      write(*,*) Rb(2,1),Rb(2,2)
*
*      write(*,*) 'slep',msel,mstau
*      write(*,*) Rtau(1,1),Rtau(1,2)
*      write(*,*) Rtau(2,1),Rtau(2,2)
*
*      write(*,*) 'sneut',msnue,msnutau
*
*      write(*,*) 'charg',mch
*      write(*,*) UU(1,1),UU(1,2),VV(1,1),VV(1,2)
*      write(*,*) UU(2,1),UU(2,2),VV(2,1),VV(2,2)
*
*      write(*,*) 'neut',mne
*      write(*,*) NN(1,1),NN(1,2),NN(1,3),NN(1,4),NN(1,5)
*      write(*,*) NN(2,1),NN(2,2),NN(2,3),NN(2,4),NN(2,5)
*      write(*,*) NN(3,1),NN(3,2),NN(3,3),NN(3,4),NN(3,5)
*      write(*,*) NN(4,1),NN(4,2),NN(4,3),NN(4,4),NN(4,5)
*      write(*,*) NN(5,1),NN(5,2),NN(5,3),NN(5,4),NN(5,5)

      end

*
***********************************************************************
*

      SUBROUTINE tree_charginos(g,ll,v1,v2,M2,xx,xmc,u,v)

      IMPLICIT NONE

      DOUBLE PRECISION g,ll,v1,v2,M2,xx,xmc(2),u(2,2),v(2,2)
      DOUBLE PRECISION mc(2,2),xtx(2,2),xxt(2,2),ut(2,2),vt(2,2),
     .     bu(2),bv(2),umc(2),vmc(2),mtemp

      INTEGER i,j,k

      mc(1,1)=M2
      mc(1,2)=g*v2
      mc(2,1)=g*v1
      mc(2,2)=ll*xx

      do i=1,2
       do j=1,2
          xxt(i,j)=0d0
          xtx(i,j)=0d0
       enddo
      enddo

      do i=1,2
       do j=1,2
          do k=1,2
             xxt(i,j)=xxt(i,j)+mc(i,k)*mc(j,k)
             xtx(i,j)=xtx(i,j)+mc(k,i)*mc(k,j)
          enddo
       enddo
      enddo


      CALL jacobi(xxt,2,umc,ut)

      CALL jacobi(xtx,2,vmc,vt)

      if(abs(vmc(1)-umc(1)).gt.1d-6) then ! swap eigenstates
       do j=1,2
          bv(j)=vt(j,1)
          vt(j,1)=vt(j,2)
          vt(j,2)=bv(j)
       enddo
      endif

      do i=1,2
       do j=1,2
          u(i,j)=ut(j,i)
          v(i,j)=vt(j,i)
       enddo
      enddo

      xmc(1)=0d0
      xmc(2)=0d0

      do i=1,2
       do j=1,2
          xmc(1)=xmc(1)+u(1,i)*mc(i,j)*v(1,j)
          xmc(2)=xmc(2)+u(2,i)*mc(i,j)*v(2,j)
       enddo
      enddo

*     order the eigenstates

      if(abs(xmc(1)).gt.abs(xmc(2))) then
       mtemp=xmc(1)
       xmc(1)=xmc(2)
       xmc(2)=mtemp
       do j=1,2
          bu(j)=u(1,j)
          u(1,j)=u(2,j)
          u(2,j)=bu(j)
          bv(j)=v(1,j)
          v(1,j)=v(2,j)
          v(2,j)=bv(j)
       enddo
      endif

      end

*
***********************************************************************
*

      SUBROUTINE tree_neutralinos(g,gp,ll,kk,v1,v2,xx,M1,M2,xmn,Z)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ll,kk,v1,v2,xx,M1,M2,xmn(5),Z(5,5)

      DOUBLE PRECISION xn(5,5),mn(5),ymn(5),xx0,xx1,sq2,zx(5,5)

      INTEGER i,i1,j,idummy,iord(5),irem(5)


      sq2=sqrt(2d0)

      xn(1,1)=M1
      xn(1,2)=0d0
      xn(1,3)=-gp*v1/sq2
      xn(1,4)=gp*v2/sq2
      xn(1,5)=0d0
      xn(2,1)=0d0
      xn(2,2)=M2
      xn(2,3)=g*v1/sq2
      xn(2,4)=-g*v2/sq2
      xn(2,5)=0d0
      xn(3,1)=xn(1,3)
      xn(3,2)=xn(2,3)
      xn(3,3)=0d0
      xn(3,4)=-ll*xx
      xn(3,5)=-ll*v2
      xn(4,1)=xn(1,4)
      xn(4,2)=xn(2,4)
      xn(4,3)=xn(3,4)
      xn(4,4)=0d0
      xn(4,5)=-ll*v1
      xn(5,1)=xn(1,5)
      xn(5,2)=xn(2,5)
      xn(5,3)=xn(3,5)
      xn(5,4)=xn(4,5)
      xn(5,5)=2d0*kk*xx


      CALL jacobi(xn,5,ymn,zx)

*     ordering the disorder

      do i=1,5
       mn(i)=abs(ymn(i))
      enddo

      xx0=dmin1(mn(1),mn(2),mn(3),mn(4),mn(5))
      xx1=dmax1(mn(1),mn(2),mn(3),mn(4),mn(5))
      idummy=1
      do i=1,5
       if(mn(i).eq.xx0)then
          iord(1)=i
       elseif(mn(i).eq.xx1)then
          iord(5)=i
       else
          irem(idummy)=i
          idummy=idummy+1
       endif
      enddo

      xx0=dmin1(mn(irem(1)),mn(irem(2)),mn(irem(3)))
      xx1=dmax1(mn(irem(1)),mn(irem(2)),mn(irem(3)))

      do i=1,3
       if(mn(irem(i)).eq.xx0)then
          iord(2)=irem(i)
       elseif(mn(irem(i)).eq.xx1)then
          iord(4)=irem(i)
       else
          iord(3)=irem(i)
       endif
      enddo
*
      do j=1,5
       i=iord(j)
       xmn(j)=ymn(i)
       do i1=1,5
          z(j,i1)=zx(i1,i)    ! note that ZX ~ Z^T
       enddo
      enddo

      end

*
***********************************************************************
*

      SUBROUTINE tree_higgses(g,gp,ll,kk,v1,v2,xx,Ak,Al,
     .     mss,maa,mhc,RS,RP,RC,errhiggs)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ll,kk,v1,v2,xx,Ak,Al,
     .     mss(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2)

      INTEGER errhiggs

      DOUBLE PRECISION gb2,As,MS(3,3),MP(3,3),ZS(3,3),ZP(3,3),
     .     ms2(3),msd(3),ma2(3),mad(3)

      DOUBLE PRECISION cb,sb,gold

      INTEGER i,j,k

      COMMON/FOREFFPOT/mt,T1,T2,st,ct,Q2
      DOUBLE PRECISION mt,T1,T2,st,ct,Q2,DMS(3,3),DMP(3,3)

      gb2=(g**2+gp**2)/2d0
      As=Al+kk*xx

      MS(1,1)=gb2*v1**2+ll*xx*v2/v1*As
      MS(1,2)=(2d0*ll**2-gb2)*v1*v2-ll*xx*As
      MS(1,3)=2d0*ll**2*v1*xx-ll*v2*(As+kk*xx)
      MS(2,2)=gb2*v2**2+ll*xx*v1/v2*As
      MS(2,3)=2d0*ll**2*v2*xx-ll*v1*(As+kk*xx)
      MS(3,3)=ll*Al*v1*v2/xx+kk*xx*(Ak+4d0*kk*xx)
      MS(2,1)=MS(1,2)
      MS(3,1)=MS(1,3)
      MS(3,2)=MS(2,3)

      gold=gb2*(v1**2+v2**2)  ! gauge-fixing mass
      sb=v2/sqrt(v1**2+v2**2)
      cb=v1/sqrt(v1**2+v2**2)

      MP(1,1)=ll*xx*v2/v1*As+cb**2*gold
      MP(1,2)=ll*xx*As-cb*sb*gold
      MP(1,3)=ll*v2*(As-3d0*kk*xx)
      MP(2,2)=ll*xx*v1/v2*As+sb**2*gold
      MP(2,3)=ll*v1*(As-3d0*kk*xx)
      MP(3,3)=4d0*ll*kk*v1*v2+ll*Al*v1*v2/xx-3d0*kk*Ak*xx
      MP(2,1)=MP(1,2)
      MP(3,1)=MP(1,3)
      MP(3,2)=MP(2,3)


      CALL jacobi(MS,3,ms2,ZS)

      CALL jacobi(MP,3,ma2,ZP)

*     take the square roots and check that it's all right

      errhiggs=0

      do i=1,3
       if(ms2(i).ge.0d0) then
          msd(i)=sqrt(ms2(i))
       else
          errhiggs=1
       endif
      enddo

      do i=1,3
       if(ma2(i).ge.-1d-12) then ! allow for tiny nonzero Goldstone mass
          mad(i)=sqrt(abs(ma2(i)))
       else
          errhiggs=2
       endif
      enddo

*     if there was a tachionic mass, try again with corrected matrices

      if(errhiggs.ne.0) then

 
      CALL effpot(1,mt,0d0,T1,T2,st,ct,Q2,v2/v1,
     .      sqrt(v1**2+v2**2),ll,xx,0d0,DMS,DMP)

       do i=1,3
          do j=1,3
             MS(i,j)=MS(i,j)+DMS(i,j)
             MP(i,j)=MP(i,j)+DMP(i,j)
          enddo
       enddo

 
      CALL jacobi(MS,3,ms2,ZS)
 
      CALL jacobi(MP,3,ma2,ZP)

       errhiggs=0

       do i=1,3
          if(ms2(i).ge.0d0) then
             msd(i)=sqrt(ms2(i))
          else
* Mod. by UE:
             msd(i)=1d0
*             errhiggs=1
*             RETURN
* End mod. by UE
          endif
       enddo

       do i=1,3
          if(ma2(i).ge.-1d-12) then ! allow for tiny nonzero Goldstone mass
             mad(i)=sqrt(abs(ma2(i)))
          else
* Mod. by UE:
            mad(i)=1d0
*             errhiggs=2
*             RETURN
* End mod. by UE
          endif
       enddo

      endif

*     order the disorder

      mss(1)=dmin1(msd(1),msd(2),msd(3))
      mss(3)=dmax1(msd(1),msd(2),msd(3))

      do i=1,3
       if(msd(i).gt.mss(1).and.msd(i).lt.mss(3)) then
          mss(2)=msd(i)
       endif
      enddo

      do i=1,3
       do j=1,3
          if(mss(i).eq.msd(j)) then
             do k=1,3
              RS(i,k)=ZS(k,j)
             enddo
          endif
       enddo
      enddo

      maa(1)=dmin1(mad(1),mad(2),mad(3))
      maa(3)=dmax1(mad(1),mad(2),mad(3))

* mod. by UE, 11.9.2014:
      IF(maa(1).le.1d-3) maa(1)=1d-3

      do i=1,3
       if(mad(i).gt.maa(1).and.mad(i).lt.maa(3)) then
          maa(2)=mad(i)
       endif
      enddo

      do i=1,3
       do j=1,3
          if(maa(i).eq.mad(j)) then
             do k=1,3
              RP(i,k)=ZP(k,j)
             enddo
          endif
       enddo
      enddo

*.     add gauge-fixing mass for G0 (improve later)
*
*      maa(1)=sqrt(gb2*(v1**2+v2**2))

*     the charged Higgses

      mhc(1)=sqrt(g**2*(v1**2+v2**2)/2d0)

      IF(((ll*xx*As-ll**2*v1*v2)*(v1**2+v2**2)/v1/v2
     .     +g**2*(v1**2+v2**2)/2d0).ge.1d0) THEN
       mhc(2)=sqrt((ll*xx*As-ll**2*v1*v2)*(v1**2+v2**2)/v1/v2
     .        +g**2*(v1**2+v2**2)/2d0)
      ELSE
       mhc(2)=1d0
      ENDIF

      RC(1,1)=-v1/sqrt(v1**2+v2**2)
      RC(1,2)=v2/sqrt(v1**2+v2**2)
      RC(2,1)=RC(1,2)
      RC(2,2)=v1/sqrt(v1**2+v2**2)

      end

*
***********************************************************************
*

      SUBROUTINE tree_sfermions(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     .     mQ3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,mer,mstop,msbot,
     .     mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue,errsfer)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     .     mQ3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,mer,mstop(2),msbot(2),
     .     mstau(2),Rt(2,2),Rb(2,2),Rtau(2,2),msup(2),msdown(2),msel(2),
     .     msnutau,msnue

      LOGICAL errsfer,errt,errb,errta

      DOUBLE PRECISION msuL,msuR,msdL,msdR,mslL,mslR

*     first the ones that mix


      CALL diagsfe(1,ht,g,gp,ll,v1,v2,xx,At,mQ3,mtr,mstop,Rt,errt)

      CALL diagsfe(2,hb,g,gp,ll,v1,v2,xx,Ab,mQ3,mbr,msbot,Rb,errb)

      CALL diagsfe(3,htau,g,gp,ll,v1,v2,xx,Atau,mL3,mtaur,mstau,Rtau,
     .     errta)

      errsfer=errt.or.errb.or.errta

*     then the others

      msuL=sqrt(mQ**2+(g**2/2d0-gp**2/6d0)*(v1**2-v2**2)/2d0)
      msdL=sqrt(mQ**2+(-g**2/2d0-gp**2/6d0)*(v1**2-v2**2)/2d0)
      mslL=sqrt(mL**2+(-g**2/2d0+gp**2/2d0)*(v1**2-v2**2)/2d0)

      msuR=sqrt(mur**2+(2d0*gp**2/3d0)*(v1**2-v2**2)/2d0)
      msdR=sqrt(mdr**2+(-gp**2/3d0)*(v1**2-v2**2)/2d0)
      mslR=sqrt(mer**2+(-gp**2)*(v1**2-v2**2)/2d0)

      msnutau=sqrt(mL3**2+(g**2/2d0+gp**2/2d0)*(v1**2-v2**2)/2d0)
      msnue  =sqrt(mL**2+(g**2/2d0+gp**2/2d0)*(v1**2-v2**2)/2d0)

      msup(1)=msuL
      msup(2)=msuR
      msdown(1)=msdL
      msdown(2)=msdR
      msel(1)=mslL
      msel(2)=mslR

      end

*
***********************************************************************
*

      SUBROUTINE diagsfe(n,hf,g,gp,ll,v1,v2,xx,Af,mL,mR,mass,R,error)

      IMPLICIT NONE

      INTEGER n
      DOUBLE PRECISION hf,g,gp,ll,v1,v2,xx,Af,mL,mR,mass(2),R(2,2)
      LOGICAL error

      DOUBLE PRECISION dL,dR,X,mf,MS(2,2),RT(2,2),mass2(2)
      INTEGER i,j

      if(n.eq.1) then
       dL=( g**2/2d0-gp**2/6d0)*(v1**2-v2**2)/2d0
       dR=( 2d0*gp**2/3d0)*(v1**2-v2**2)/2d0
       mf=hf*v2
       X=Af-ll*xx*v1/v2
      elseif(n.eq.2) then
       dL=(-g**2/2d0-gp**2/6d0)*(v1**2-v2**2)/2d0
       dR=(  -gp**2/3d0)*(v1**2-v2**2)/2d0
       mf=hf*v1
       X=Af-ll*xx*v2/v1
      elseif(n.eq.3) then
       dL=(-g**2/2d0+gp**2/2d0)*(v1**2-v2**2)/2d0
       dR=(  -gp**2)*(v1**2-v2**2)/2d0
       mf=hf*v1
       X=Af-ll*xx*v2/v1
      endif

      MS(1,1)=mL**2+mf**2+dL
      MS(1,2)=mf*X
      MS(2,1)=MS(1,2)
      MS(2,2)=mR**2+mf**2+dR


      CALL jacobi(MS,2,mass2,Rt)

      error=.false.

      if(mass2(1).ge.0d0) then
       mass(1)=sqrt(mass2(1))
      else
       error=.true.
       mass(1)=-sqrt(abs(mass2(1)))
      endif

      if(mass2(2).ge.0d0) then
       mass(2)=sqrt(mass2(2))
      else
       error=.true.
       mass(2)=-sqrt(abs(mass2(2)))
      endif

      do i=1,2
       do j=1,2
          R(i,j)=Rt(j,i)
       enddo
      enddo

      end

*
***********************************************************************
*

      SUBROUTINE scalarcouplings(g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
     .    Al,Ak,At,Ab,Atau)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
     .    Al,Ak,At,Ab,Atau

      COMMON/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     .     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue

      DOUBLE PRECISION mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     .     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),
     .     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     .     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2)

      COMMON/coupl_s/lsstt,lssbb,lsstata,lssntnt,lssuu,lssdd,
     .     lssee,lssnn,lstt,lsbb,lstata,lsntnt,lsuu,lsdd,lsee,lsnn,
     .     lsshh,lssaa,lshh,lsaa,lsscc,lscc,lsnene,lschch

      DOUBLE PRECISION lsstt(3,3,2,2),lssbb(3,3,2,2),
     .     lsstata(3,3,2,2),lssntnt(3,3),lssuu(3,3,2,2),lssdd(3,3,2,2),
     .     lssee(3,3,2,2),lssnn(3,3),lstt(3,2,2),lsbb(3,2,2),
     .     lstata(3,2,2),lsntnt(3),lsuu(3,2,2),lsdd(3,2,2),lsee(3,2,2),
     .     lsnn(3),lsshh(3,3,3,3),lssaa(3,3,3,3),lshh(3,3,3),
     .     lsaa(3,3,3),lsscc(3,3,2,2),lscc(3,2,2),
     .     lsnene(3,5,5),lschch(3,2,2)

      COMMON/coupl_p/lpptt,lppbb,lpptata,lppntnt,lppuu,lppdd,
     .     lppee,lppnn,lptt,lpbb,lptata,lpntnt,lpuu,lpdd,lpee,lpnn,
     .     lpphh,lppaa,lpah,lppcc,lpcc,lpnene,lpchch

      DOUBLE PRECISION lpptt(3,3,2,2),lppbb(3,3,2,2),
     .     lpptata(3,3,2,2),lppntnt(3,3),lppuu(3,3,2,2),lppdd(3,3,2,2),
     .     lppee(3,3,2,2),lppnn(3,3),lptt(3,2,2),lpbb(3,2,2),
     .     lptata(3,2,2),lpntnt(3),lpuu(3,2,2),lpdd(3,2,2),lpee(3,2,2),
     .     lpnn(3),lpphh(3,3,3,3),lppaa(3,3,3,3),lpah(3,3,3),
     .     lppcc(3,3,2,2),lpcc(3,2,2),lpnene(3,5,5),lpchch(3,2,2)

*     fill the COMMON blocks

      CALL coupl_s_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     .     Rt,Rb,Rtau,lsstt,lssbb,lsstata,lssntnt,lssuu,lssdd,
     .     lssee,lssnn,lstt,lsbb,lstata,lsntnt,lsuu,lsdd,lsee,lsnn)


      CALL coupl_s_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
     .     lsshh,lssaa,lshh,lsaa,lsscc,lscc)


      CALL coupl_s_ino(g,gp,ll,kk,NN,UU,VV,lsnene,lschch)


      CALL coupl_p_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     .     Rt,Rb,Rtau,lpptt,lppbb,lpptata,lppntnt,lppuu,lppdd,
     .     lppee,lppnn,lptt,lpbb,lptata,lpntnt,lpuu,lpdd,lpee,lpnn)


      CALL coupl_p_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
     .     lpphh,lppaa,lpah,lppcc,lpcc)


      CALL coupl_p_ino(g,gp,ll,kk,NN,UU,VV,lpnene,lpchch)

      end

*
***********************************************************************
*

      SUBROUTINE coupl_s_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     .     Rt,Rb,Rtau,lsstt,lssbb,lsstata,lssntnt,lssuu,lssdd,
     .     lssee,lssnn,lstt,lsbb,lstata,lsntnt,lsuu,lsdd,lsee,lsnn)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     .     Rt(2,2),Rb(2,2),Rtau(2,2),lsstt(3,3,2,2),lssbb(3,3,2,2),
     .     lsstata(3,3,2,2),lssntnt(3,3),lssuu(3,3,2,2),lssdd(3,3,2,2),
     .     lssee(3,3,2,2),lssnn(3,3),lstt(3,2,2),lsbb(3,2,2),
     .     lstata(3,2,2),lsntnt(3),lsuu(3,2,2),lsdd(3,2,2),lsee(3,2,2),
     .     lsnn(3),lsstt_int(3,3,2,2),lssbb_int(3,3,2,2),
     .     lsstata_int(3,3,2,2),lstt_int(3,2,2),lsbb_int(3,2,2),
     .     lstata_int(3,2,2)

      DOUBLE PRECISION sq2,gb2,sw2,guL,guR,gdL,gdR,geL,geR,gnu

      INTEGER i,j,k,l,a,b

      sq2=sqrt(2d0)
      gb2=(g**2+gp**2)/2d0
      sw2=gp**2/2d0/gb2

      guL= .5d0-2d0/3d0*sw2
      gdL=-.5d0+1d0/3d0*sw2
      gnu= .5d0
      geL=-.5d0+sw2
      guR= 2d0/3d0*sw2
      gdR=-1d0/3d0*sw2
      geR=-sw2

*     FIRST TWO GENERATIONS

*     quartic, up squarks

      do i=1,3
       do j=1,3
          do k=1,2
             do l=1,2
              lssuu(i,j,k,l)=0d0
             enddo
          enddo
       enddo
      enddo

      lssuu(1,1,1,1)=gb2*guL/2d0
      lssuu(1,1,2,2)=gb2*guR/2d0
      lssuu(2,2,1,1)=-gb2*guL/2d0
      lssuu(2,2,2,2)=-gb2*guR/2d0

*     trilinear, up squarks

      do i=1,3
       do k=1,2
          do l=1,2
             lsuu(i,k,l)=0d0
          enddo
       enddo
      enddo

      lsuu(1,1,1)=sq2*gb2*guL*v1
      lsuu(1,2,2)=sq2*gb2*guR*v1
      lsuu(2,1,1)=-sq2*gb2*guL*v2
      lsuu(2,2,2)=-sq2*gb2*guR*v2

*     quartic, down squarks

      do i=1,3
       do j=1,3
          do k=1,2
             do l=1,2
              lssdd(i,j,k,l)=0d0
             enddo
          enddo
       enddo
      enddo

      lssdd(1,1,1,1)=gb2*gdL/2d0
      lssdd(1,1,2,2)=gb2*gdR/2d0
      lssdd(2,2,1,1)=-gb2*gdL/2d0
      lssdd(2,2,2,2)=-gb2*gdR/2d0

*     trilinear, down squarks

      do i=1,3
       do k=1,2
          do l=1,2
             lsdd(i,k,l)=0d0
          enddo
       enddo
      enddo

      lsdd(1,1,1)=sq2*gb2*gdL*v1
      lsdd(1,2,2)=sq2*gb2*gdR*v1
      lsdd(2,1,1)=-sq2*gb2*gdL*v2
      lsdd(2,2,2)=-sq2*gb2*gdR*v2

*     quartic, charged sleptons

      do i=1,3
       do j=1,3
          do k=1,2
             do l=1,2
              lssee(i,j,k,l)=0d0
             enddo
          enddo
       enddo
      enddo

      lssee(1,1,1,1)=gb2*geL/2d0
      lssee(1,1,2,2)=gb2*geR/2d0
      lssee(2,2,1,1)=-gb2*geL/2d0
      lssee(2,2,2,2)=-gb2*geR/2d0

*     trilinear, charged sleptons

      do i=1,3
       do k=1,2
          do l=1,2
             lsee(i,k,l)=0d0
          enddo
       enddo
      enddo

      lsee(1,1,1)=sq2*gb2*geL*v1
      lsee(1,2,2)=sq2*gb2*geR*v1
      lsee(2,1,1)=-sq2*gb2*geL*v2
      lsee(2,2,2)=-sq2*gb2*geR*v2

*     quartic, sneutrinos

      do i=1,3
       do j=1,3
          lssnn(i,j)=0d0
       enddo
      enddo

      lssnn(1,1)=gb2*gnu/2d0
      lssnn(2,2)=-gb2*gnu/2d0

*     trilinear, sneutrinos

      lsnn(1)=sq2*gb2*gnu*v1
      lsnn(2)=-sq2*gb2*gnu*v2
      lsnn(3)=0d0

*     THIRD GENERATION

*     quartic, top squarks

      do i=1,3
       do j=1,3
          do k=1,2
             do l=1,2
              lsstt_int(i,j,k,l)=lssuu(i,j,k,l)
             enddo
          enddo
       enddo
      enddo

      lsstt_int(2,2,1,1)=lsstt_int(2,2,1,1)+ht**2/2d0
      lsstt_int(2,2,2,2)=lsstt_int(2,2,2,2)+ht**2/2d0

      lsstt_int(1,3,1,2)=lsstt_int(1,3,1,2)-ht*ll/4d0
      lsstt_int(1,3,2,1)=lsstt_int(1,3,1,2)
      lsstt_int(3,1,1,2)=lsstt_int(1,3,1,2)
      lsstt_int(3,1,2,1)=lsstt_int(1,3,1,2)

*     trilinear, top squarks

      do i=1,3
       do k=1,2
          do l=1,2
             lstt_int(i,k,l)=lsuu(i,k,l)
          enddo
       enddo
      enddo

      lstt_int(2,1,1)=lstt_int(2,1,1)+sq2*ht**2*v2
      lstt_int(2,2,2)=lstt_int(2,2,2)+sq2*ht**2*v2

      lstt_int(1,1,2)=lstt_int(1,1,2)-ht*ll*xx/sq2
      lstt_int(2,1,2)=lstt_int(2,1,2)+ht*At/sq2
      lstt_int(3,1,2)=lstt_int(3,1,2)-ht*ll*v1/sq2

      lstt_int(1,2,1)=lstt_int(1,1,2)
      lstt_int(2,2,1)=lstt_int(2,1,2)
      lstt_int(3,2,1)=lstt_int(3,1,2)

*     quartic, bottom squarks

      do i=1,3
       do j=1,3
          do k=1,2
             do l=1,2
              lssbb_int(i,j,k,l)=lssdd(i,j,k,l)
             enddo
          enddo
       enddo
      enddo

      lssbb_int(1,1,1,1)=lssbb_int(1,1,1,1)+hb**2/2d0
      lssbb_int(1,1,2,2)=lssbb_int(1,1,2,2)+hb**2/2d0

      lssbb_int(2,3,1,2)=lssbb_int(2,3,1,2)-hb*ll/4d0
      lssbb_int(2,3,2,1)=lssbb_int(2,3,1,2)
      lssbb_int(3,2,1,2)=lssbb_int(2,3,1,2)
      lssbb_int(3,2,2,1)=lssbb_int(2,3,1,2)

*     trilinear, bottom squarks

      do i=1,3
       do k=1,2
          do l=1,2
             lsbb_int(i,k,l)=lsdd(i,k,l)
          enddo
       enddo
      enddo

      lsbb_int(1,1,1)=lsbb_int(1,1,1)+sq2*hb**2*v1
      lsbb_int(1,2,2)=lsbb_int(1,2,2)+sq2*hb**2*v1

      lsbb_int(1,1,2)=lsbb_int(1,1,2)+hb*Ab/sq2
      lsbb_int(2,1,2)=lsbb_int(2,1,2)-hb*ll*xx/sq2
      lsbb_int(3,1,2)=lsbb_int(3,1,2)-hb*ll*v2/sq2

      lsbb_int(1,2,1)=lsbb_int(1,1,2)
      lsbb_int(2,2,1)=lsbb_int(2,1,2)
      lsbb_int(3,2,1)=lsbb_int(3,1,2)

*     quartic, staus

      do i=1,3
       do j=1,3
          do k=1,2
             do l=1,2
               lsstata_int(i,j,k,l)=lssee(i,j,k,l)
             enddo
          enddo
       enddo
      enddo

      lsstata_int(1,1,1,1)=lsstata_int(1,1,1,1)+htau**2/2d0
      lsstata_int(1,1,2,2)=lsstata_int(1,1,2,2)+htau**2/2d0

      lsstata_int(2,3,1,2)=lsstata_int(2,3,1,2)-htau*ll/4d0
      lsstata_int(2,3,2,1)=lsstata_int(2,3,1,2)
      lsstata_int(3,2,1,2)=lsstata_int(2,3,1,2)
      lsstata_int(3,2,2,1)=lsstata_int(2,3,1,2)

*     trilinear, staus

      do i=1,3
       do k=1,2
          do l=1,2
             lstata_int(i,k,l)=lsee(i,k,l)
          enddo
       enddo
      enddo

      lstata_int(1,1,1)=lstata_int(1,1,1)+sq2*htau**2*v1
      lstata_int(1,2,2)=lstata_int(1,2,2)+sq2*htau**2*v1

      lstata_int(1,1,2)=lstata_int(1,1,2)+htau*Atau/sq2
      lstata_int(2,1,2)=lstata_int(2,1,2)-htau*ll*xx/sq2
      lstata_int(3,1,2)=lstata_int(3,1,2)-htau*ll*v2/sq2

      lstata_int(1,2,1)=lstata_int(1,1,2)
      lstata_int(2,2,1)=lstata_int(2,1,2)
      lstata_int(3,2,1)=lstata_int(3,1,2)

*     sneutrinos (same as first two generations)

      do i=1,3
       lsntnt(i)=lsnn(i)
       do j=1,3
          lssntnt(i,j)=lssnn(i,j)
       enddo
      enddo

*     now rotate the third-generation couplings

*     quartic
      do i=1,3
       do j=1,3
          do k=1,2
             do l=1,2
              lsstt(i,j,k,l)=0d0
              lssbb(i,j,k,l)=0d0
              lsstata(i,j,k,l)=0d0
              do a=1,2
                 do b=1,2
                  lsstt(i,j,k,l)=lsstt(i,j,k,l)
     .                   +Rt(k,a)*Rt(l,b)*lsstt_int(i,j,a,b)
                  lssbb(i,j,k,l)=lssbb(i,j,k,l)
     .                   +Rb(k,a)*Rb(l,b)*lssbb_int(i,j,a,b)
                  lsstata(i,j,k,l)=lsstata(i,j,k,l)
     .                   +Rtau(k,a)*Rtau(l,b)*lsstata_int(i,j,a,b)
                 enddo
              enddo
             enddo
          enddo
       enddo
      enddo

*     trilinear
      do i=1,3
       do k=1,2
          do l=1,2
             lstt(i,k,l)=0d0
             lsbb(i,k,l)=0d0
             lstata(i,k,l)=0d0
             do a=1,2
              do b=1,2
                 lstt(i,k,l)=lstt(i,k,l)
     .                +Rt(k,a)*Rt(l,b)*lstt_int(i,a,b)
                 lsbb(i,k,l)=lsbb(i,k,l)
     .                +Rb(k,a)*Rb(l,b)*lsbb_int(i,a,b)
                 lstata(i,k,l)=lstata(i,k,l)
     .                +Rtau(k,a)*Rtau(l,b)*lstata_int(i,a,b)
              enddo
             enddo
          enddo
       enddo
      enddo

      end

*
***********************************************************************
*

      SUBROUTINE coupl_s_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
     .     lsshh,lssaa,lshh,lsaa,lsscc,lscc)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ll,kk,v1,v2,xx,Al,Ak,RS(3,3),RP(3,3),
     .     lsshh(3,3,3,3),lssaa(3,3,3,3),lshh(3,3,3),lsaa(3,3,3),
     .     lsscc(3,3,2,2),lscc(3,2,2)

      DOUBLE PRECISION lssss(3,3,3,3),lsspp(3,3,3,3),
     .     lsss(3,3,3),lspp(3,3,3),gb2,sq2,c2b,s2b

      INTEGER i,j,k,l,a,b

      gb2=(g**2+gp**2)/2d0
      sq2=sqrt(2d0)

      do i=1,3              ! initialize
       do j=1,3
          do k=1,3
             lsss(i,j,k)=0d0
             lspp(i,j,k)=0d0
             lshh(i,j,k)=0d0
             lsaa(i,j,k)=0d0
             do l=1,3
              lssss(i,j,k,l)=0d0
              lsspp(i,j,k,l)=0d0
              lsshh(i,j,k,l)=0d0
              lssaa(i,j,k,l)=0d0
             enddo
          enddo
          do k=1,2
             do l=1,2
              lsscc(i,j,k,l)=0d0
             enddo
          enddo
       enddo
       do j=1,2
          do k=1,2
             lscc(i,j,k)=0d0
          enddo
       enddo
      enddo

*     quartic neutral couplings

      lssss(1,1,1,1)=gb2/16d0

      lssss(2,2,2,2)=gb2/16d0

      lssss(3,3,3,3)=kk**2/4d0

      lssss(1,1,2,2)=(2d0*ll**2-gb2)/48d0
      lssss(1,2,1,2)=lssss(1,1,2,2) ! symmetrize the couplings
      lssss(1,2,2,1)=lssss(1,1,2,2) ! (got a better proposal? ;-)
      lssss(2,1,1,2)=lssss(1,1,2,2)
      lssss(2,1,2,1)=lssss(1,1,2,2)
      lssss(2,2,1,1)=lssss(1,1,2,2)

      lssss(1,2,3,3)=-ll*kk/24d0
      lssss(2,1,3,3)=lssss(1,2,3,3)
      lssss(2,3,1,3)=lssss(1,2,3,3)
      lssss(1,3,2,3)=lssss(1,2,3,3)
      lssss(3,1,2,3)=lssss(1,2,3,3)
      lssss(3,2,1,3)=lssss(1,2,3,3)
      lssss(3,1,3,2)=lssss(1,2,3,3)
      lssss(3,2,3,1)=lssss(1,2,3,3)
      lssss(3,3,1,2)=lssss(1,2,3,3)
      lssss(3,3,2,1)=lssss(1,2,3,3)
      lssss(1,3,3,2)=lssss(1,2,3,3)
      lssss(2,3,3,1)=lssss(1,2,3,3)

      lssss(1,1,3,3)=ll**2/24d0
      lssss(1,3,1,3)=lssss(1,1,3,3)
      lssss(1,3,3,1)=lssss(1,1,3,3)
      lssss(3,1,1,3)=lssss(1,1,3,3)
      lssss(3,1,3,1)=lssss(1,1,3,3)
      lssss(3,3,1,1)=lssss(1,1,3,3)
      lssss(2,2,3,3)=lssss(1,1,3,3)
      lssss(2,3,2,3)=lssss(1,1,3,3)
      lssss(2,3,3,2)=lssss(1,1,3,3)
      lssss(3,2,2,3)=lssss(1,1,3,3)
      lssss(3,2,3,2)=lssss(1,1,3,3)
      lssss(3,3,2,2)=lssss(1,1,3,3)

      lsspp(1,1,1,1)=gb2/8d0
      lsspp(2,2,2,2)=lsspp(1,1,1,1)

      lsspp(1,1,2,2)=(2d0*ll**2-gb2)/8d0
      lsspp(2,2,1,1)=lsspp(1,1,2,2)

      lsspp(1,1,3,3)=ll**2/4d0
      lsspp(2,2,3,3)=lsspp(1,1,3,3)
      lsspp(3,3,1,1)=lsspp(1,1,3,3)
      lsspp(3,3,2,2)=lsspp(1,1,3,3)

      lsspp(1,2,3,3)=ll*kk/4d0
      lsspp(2,1,3,3)=lsspp(1,2,3,3)
      lsspp(3,3,1,2)=lsspp(1,2,3,3)
      lsspp(3,3,2,1)=lsspp(1,2,3,3)

      lsspp(1,3,2,3)=-ll*kk/4d0
      lsspp(3,1,2,3)=lsspp(1,3,2,3)
      lsspp(1,3,3,2)=lsspp(1,3,2,3)
      lsspp(3,1,3,2)=lsspp(1,3,2,3)
      lsspp(2,3,1,3)=lsspp(1,3,2,3)
      lsspp(3,2,1,3)=lsspp(1,3,2,3)
      lsspp(2,3,3,1)=lsspp(1,3,2,3)
      lsspp(3,2,3,1)=lsspp(1,3,2,3)

      lsspp(3,3,3,3)=kk**2/2d0

*     rotate the quartics

      do i=1,3
       do j=1,3
          do k=1,3
             do l=1,3
              do a=1,3
                 do b=1,3
                  lsshh(i,j,k,l)=lsshh(i,j,k,l)
     .                   +6d0*RS(k,a)*RS(l,b)*lssss(i,j,a,b)
                  lssaa(i,j,k,l)=lssaa(i,j,k,l)
     .                   +RP(k,a)*RP(l,b)*lsspp(i,j,a,b)
                 enddo
              enddo
             enddo
          enddo
       enddo
      enddo

*     trilinear neutral couplings

      lsss(1,1,1)=gb2*v1/2d0/sq2

      lsss(2,2,2)=gb2*v2/2d0/sq2

      lsss(1,2,2)=v1/6d0/sq2*(2d0*ll**2-gb2)
      lsss(2,1,2)=lsss(1,2,2)
      lsss(2,2,1)=lsss(1,2,2)

      lsss(1,1,2)=v2/6d0/sq2*(2d0*ll**2-gb2)
      lsss(1,2,1)=lsss(1,1,2)
      lsss(2,1,1)=lsss(1,1,2)

      lsss(3,1,1)=ll**2*xx/3d0/sq2
      lsss(1,1,3)=lsss(3,1,1)
      lsss(1,3,1)=lsss(3,1,1)
      lsss(3,2,2)=lsss(3,1,1)
      lsss(2,2,3)=lsss(3,1,1)
      lsss(2,3,2)=lsss(3,1,1)

      lsss(3,3,3)=kk*Ak/3d0/sq2+sq2*kk**2*xx

      lsss(1,3,3)=ll/3d0/sq2*(ll*v1-kk*v2)
      lsss(3,1,3)=lsss(1,3,3)
      lsss(3,3,1)=lsss(1,3,3)

      lsss(2,3,3)=ll/3d0/sq2*(ll*v2-kk*v1)
      lsss(3,2,3)=lsss(2,3,3)
      lsss(3,3,2)=lsss(2,3,3)

      lsss(1,2,3)=-ll*Al/6d0/sq2-ll*kk*xx/3d0/sq2
      lsss(1,3,2)=lsss(1,2,3)
      lsss(2,1,3)=lsss(1,2,3)
      lsss(2,3,1)=lsss(1,2,3)
      lsss(3,1,2)=lsss(1,2,3)
      lsss(3,2,1)=lsss(1,2,3)


      lspp(1,1,1)=gb2*v1/2d0/sq2

      lspp(2,2,2)=gb2*v2/2d0/sq2

      lspp(1,2,2)=v1/2d0/sq2*(2d0*ll**2-gb2)

      lspp(2,1,1)=v2/2d0/sq2*(2d0*ll**2-gb2)

      lspp(3,1,1)=ll**2*xx/sq2
      lspp(3,2,2)=lspp(3,1,1)

      lspp(3,3,3)=-kk*Ak/sq2+sq2*kk**2*xx

      lspp(1,3,3)=ll/sq2*(ll*v1+kk*v2)

      lspp(2,3,3)=ll/sq2*(ll*v2+kk*v1)

      lspp(3,1,3)=-ll*kk*v2/sq2
      lspp(3,3,1)=lspp(3,1,3)

      lspp(3,2,3)=-ll*kk*v1/sq2
      lspp(3,3,2)=lspp(3,2,3)

      lspp(1,2,3)=ll*Al/2d0/sq2-ll*kk*xx/sq2
      lspp(1,3,2)=lspp(1,2,3)
      lspp(2,1,3)=lspp(1,2,3)
      lspp(2,3,1)=lspp(1,2,3)

      lspp(3,1,2)=ll*Al/2d0/sq2+ll*kk*xx/sq2
      lspp(3,2,1)=lspp(3,1,2)

*     rotate the trilinears

      do i=1,3
       do k=1,3
          do l=1,3
             do a=1,3
              do b=1,3
                 lshh(i,k,l)=lshh(i,k,l)
     .                +3d0*RS(k,a)*RS(l,b)*lsss(i,a,b)
                 lsaa(i,k,l)=lsaa(i,k,l)
     .                +RP(k,a)*RP(l,b)*lspp(i,a,b)
              enddo
             enddo
          enddo
       enddo
      enddo

*     quartic charged couplings

      c2b=(v1**2-v2**2)/(v1**2+v2**2)
      s2b=2d0*v1*v2/(v1**2+v2**2)

      lsscc(1,1,1,1)=(g**2+gp**2*c2b)/8d0
      lsscc(2,2,2,2)=lsscc(1,1,1,1)

      lsscc(1,1,2,2)=(g**2-gp**2*c2b)/8d0
      lsscc(2,2,1,1)=lsscc(1,1,2,2)

      lsscc(1,2,1,1)=(2d0*ll**2-g**2)*s2b/8d0
      lsscc(2,1,1,1)=lsscc(1,2,1,1)

      lsscc(1,2,2,2)=-(2d0*ll**2-g**2)*s2b/8d0
      lsscc(2,1,2,2)=lsscc(1,2,2,2)

      lsscc(1,1,1,2)=-gp**2/8d0*s2b
      lsscc(1,1,2,1)=lsscc(1,1,1,2)

      lsscc(2,2,1,2)=gp**2/8d0*s2b
      lsscc(2,2,2,1)=lsscc(2,2,1,2)

      lsscc(1,2,1,2)=(2d0*ll**2-g**2)*c2b/8d0
      lsscc(1,2,2,1)=lsscc(1,2,1,2)
      lsscc(2,1,1,2)=lsscc(1,2,1,2)
      lsscc(2,1,2,1)=lsscc(1,2,1,2)

      lsscc(3,3,1,2)=-kk*ll/2d0*c2b
      lsscc(3,3,2,1)=lsscc(3,3,1,2)

      lsscc(3,3,1,1)=ll/2d0*(ll-kk*s2b)

      lsscc(3,3,2,2)=ll/2d0*(ll+kk*s2b)

*     trilinear charged couplings

      lscc(1,1,1)=(v1*(g**2+gp**2*c2b)+v2*(2d0*ll**2-g**2)*s2b)/2d0/sq2

      lscc(1,2,2)=(v1*(g**2-gp**2*c2b)-v2*(2d0*ll**2-g**2)*s2b)/2d0/sq2

      lscc(1,1,2)=(-v1*gp**2*s2b+v2*(2d0*ll**2-g**2)*c2b)/2d0/sq2
      lscc(1,2,1)=lscc(1,1,2)

      lscc(2,1,1)=(v2*(g**2-gp**2*c2b)+v1*(2d0*ll**2-g**2)*s2b)/2d0/sq2

      lscc(2,2,2)=(v2*(g**2+gp**2*c2b)-v1*(2d0*ll**2-g**2)*s2b)/2d0/sq2

      lscc(2,1,2)=(v2*gp**2*s2b+v1*(2d0*ll**2-g**2)*c2b)/2d0/sq2
      lscc(2,2,1)=lscc(2,1,2)

      lscc(3,1,1)=ll/sq2*(2d0*ll*xx-(Al+2d0*kk*xx)*s2b)

      lscc(3,2,2)=ll/sq2*(2d0*ll*xx+(Al+2d0*kk*xx)*s2b)

      lscc(3,1,2)=-ll/sq2*(Al+2d0*kk*xx)*c2b
      lscc(3,2,1)=lscc(3,1,2)

      end

*
***********************************************************************
*

      SUBROUTINE coupl_s_ino(g,gp,ll,kk,NN,UU,VV,lsnene,lschch)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ll,kk,NN(5,5),UU(2,2),VV(2,2),
     .     lsnene(3,5,5),lschch(3,2,2)

      DOUBLE PRECISION ls00(3,5,5),lspm(3,2,2),sq2

      INTEGER i,k,l,a,b

      sq2=sqrt(2d0)

*     neutralino couplings

      do i=1,3
       do k=1,5
          do l=1,5
             ls00(i,k,l)=0d0
             lsnene(i,k,l)=0d0
          enddo
       enddo
      enddo

      ls00(1,1,3) =-gp/4d0
      ls00(1,3,1)=ls00(1,1,3)

      ls00(2,1,4)=gp/4d0
      ls00(2,4,1)=ls00(2,1,4)

      ls00(1,2,3)=g/4d0
      ls00(1,3,2)=ls00(1,2,3)

      ls00(2,2,4)=-g/4d0
      ls00(2,4,2)=ls00(2,2,4)

      ls00(3,5,5)=kk/sq2

      ls00(1,4,5)=-ll/2d0/sq2
      ls00(1,5,4)=ls00(1,4,5)
      ls00(2,3,5)=ls00(1,4,5)
      ls00(2,5,3)=ls00(1,4,5)
      ls00(3,3,4)=ls00(1,4,5)
      ls00(3,4,3)=ls00(1,4,5)

      do i=1,3
       do k=1,5
          do l=1,5
             do a=1,5
              do b=1,5
                 lsnene(i,k,l)=lsnene(i,k,l)
     .                +NN(k,a)*NN(l,b)*ls00(i,a,b)
              enddo
             enddo
          enddo
       enddo
      enddo

*     chargino couplings

      do i=1,3
       do k=1,2
          do l=1,2
             lspm(i,k,l)=0d0
             lschch(i,k,l)=0d0
          enddo
       enddo
      enddo

      lspm(1,1,2)=g/sq2

      lspm(2,2,1)=g/sq2

      lspm(3,2,2)=ll/sq2

      do i=1,3
       do k=1,2
          do l=1,2
             do a=1,2
              do b=1,2
                 lschch(i,k,l)=lschch(i,k,l)
     .                +VV(k,a)*UU(l,b)*lspm(i,a,b)
              enddo
             enddo
          enddo
       enddo
      enddo

      end

*
***********************************************************************
*

      SUBROUTINE coupl_p_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     .     Rt,Rb,Rtau,lpptt,lppbb,lpptata,lppntnt,lppuu,lppdd,
     .     lppee,lppnn,lptt,lpbb,lptata,lpntnt,lpuu,lpdd,lpee,lpnn)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     .     Rt(2,2),Rb(2,2),Rtau(2,2),lpptt(3,3,2,2),lppbb(3,3,2,2),
     .     lpptata(3,3,2,2),lppntnt(3,3),lppuu(3,3,2,2),lppdd(3,3,2,2),
     .     lppee(3,3,2,2),lppnn(3,3),lptt(3,2,2),lpbb(3,2,2),
     .     lptata(3,2,2),lpntnt(3),lpuu(3,2,2),lpdd(3,2,2),lpee(3,2,2),
     .     lpnn(3),lpptt_int(3,3,2,2),lppbb_int(3,3,2,2),
     .     lpptata_int(3,3,2,2),lptt_int(3,2,2),lpbb_int(3,2,2),
     .     lptata_int(3,2,2)

      DOUBLE PRECISION sq2,gb2,sw2,guL,guR,gdL,gdR,geL,geR,gnu

      INTEGER i,j,k,l,a,b

      sq2=sqrt(2d0)
      gb2=(g**2+gp**2)/2d0
      sw2=gp**2/2d0/gb2

      guL= .5d0-2d0/3d0*sw2
      gdL=-.5d0+1d0/3d0*sw2
      gnu= .5d0
      geL=-.5d0+sw2
      guR= 2d0/3d0*sw2
      gdR=-1d0/3d0*sw2
      geR=-sw2

*     FIRST TWO GENERATIONS

*     quartic, up squarks

      do i=1,3
       do j=1,3
          do k=1,2
             do l=1,2
              lppuu(i,j,k,l)=0d0
             enddo
          enddo
       enddo
      enddo

      lppuu(1,1,1,1)=gb2*guL/2d0
      lppuu(1,1,2,2)=gb2*guR/2d0
      lppuu(2,2,1,1)=-gb2*guL/2d0
      lppuu(2,2,2,2)=-gb2*guR/2d0

*     trilinear, up squarks

      do i=1,3
       do k=1,2
          do l=1,2
             lpuu(i,k,l)=0d0
          enddo
       enddo
      enddo

*     quartic, down squarks

      do i=1,3
       do j=1,3
          do k=1,2
             do l=1,2
              lppdd(i,j,k,l)=0d0
             enddo
          enddo
       enddo
      enddo

      lppdd(1,1,1,1)=gb2*gdL/2d0
      lppdd(1,1,2,2)=gb2*gdR/2d0
      lppdd(2,2,1,1)=-gb2*gdL/2d0
      lppdd(2,2,2,2)=-gb2*gdR/2d0

*     trilinear, down squarks

      do i=1,3
       do k=1,2
          do l=1,2
             lpdd(i,k,l)=0d0
          enddo
       enddo
      enddo

*     quartic, charged sleptons

      do i=1,3
       do j=1,3
          do k=1,2
             do l=1,2
              lppee(i,j,k,l)=0d0
             enddo
          enddo
       enddo
      enddo

      lppee(1,1,1,1)=gb2*geL/2d0
      lppee(1,1,2,2)=gb2*geR/2d0
      lppee(2,2,1,1)=-gb2*geL/2d0
      lppee(2,2,2,2)=-gb2*geR/2d0

*     trilinear, charged sleptons

      do i=1,3
       do k=1,2
          do l=1,2
             lpee(i,k,l)=0d0
          enddo
       enddo
      enddo

*     quartic, sneutrinos

      do i=1,3
       do j=1,3
          lppnn(i,j)=0d0
       enddo
      enddo

      lppnn(1,1)=gb2*gnu/2d0
      lppnn(2,2)=-gb2*gnu/2d0

*     trilinear, sneutrinos

      lpnn(1)=0d0
      lpnn(2)=0d0
      lpnn(3)=0d0

*     THIRD GENERATION

*     quartic, top squarks

      do i=1,3
       do j=1,3
          do k=1,2
             do l=1,2
              lpptt_int(i,j,k,l)=lppuu(i,j,k,l)
             enddo
          enddo
       enddo
      enddo

      lpptt_int(2,2,1,1)=lpptt_int(2,2,1,1)+ht**2/2d0
      lpptt_int(2,2,2,2)=lpptt_int(2,2,2,2)+ht**2/2d0

      lpptt_int(1,3,1,2)=lpptt_int(1,3,1,2)+ht*ll/4d0
      lpptt_int(1,3,2,1)=lpptt_int(1,3,1,2)
      lpptt_int(3,1,1,2)=lpptt_int(1,3,1,2)
      lpptt_int(3,1,2,1)=lpptt_int(1,3,1,2)

*     trilinear, top squarks

      do i=1,3
       do k=1,2
          do l=1,2
             lptt_int(i,k,l)=0d0
          enddo
       enddo
      enddo

      lptt_int(1,1,2)=ht*ll*xx/sq2
      lptt_int(2,1,2)=ht*At/sq2
      lptt_int(3,1,2)=ht*ll*v1/sq2

      lptt_int(1,2,1)=-lptt_int(1,1,2)
      lptt_int(2,2,1)=-lptt_int(2,1,2)
      lptt_int(3,2,1)=-lptt_int(3,1,2)

*     quartic, bottom squarks

      do i=1,3
       do j=1,3
          do k=1,2
             do l=1,2
              lppbb_int(i,j,k,l)=lppdd(i,j,k,l)
             enddo
          enddo
       enddo
      enddo

      lppbb_int(1,1,1,1)=lppbb_int(1,1,1,1)+hb**2/2d0
      lppbb_int(1,1,2,2)=lppbb_int(1,1,2,2)+hb**2/2d0

      lppbb_int(2,3,1,2)=lppbb_int(2,3,1,2)+hb*ll/4d0
      lppbb_int(2,3,2,1)=lppbb_int(2,3,1,2)
      lppbb_int(3,2,1,2)=lppbb_int(2,3,1,2)
      lppbb_int(3,2,2,1)=lppbb_int(2,3,1,2)

*     trilinear, bottom squarks

      do i=1,3
       do k=1,2
          do l=1,2
             lpbb_int(i,k,l)=0d0
          enddo
       enddo
      enddo

      lpbb_int(1,1,2)=hb*Ab/sq2
      lpbb_int(2,1,2)=hb*ll*xx/sq2
      lpbb_int(3,1,2)=hb*ll*v2/sq2

      lpbb_int(1,2,1)=-lpbb_int(1,1,2)
      lpbb_int(2,2,1)=-lpbb_int(2,1,2)
      lpbb_int(3,2,1)=-lpbb_int(3,1,2)

*     quartic, staus

      do i=1,3
       do j=1,3
          do k=1,2
             do l=1,2
               lpptata_int(i,j,k,l)=lppee(i,j,k,l)
             enddo
          enddo
       enddo
      enddo

      lpptata_int(1,1,1,1)=lpptata_int(1,1,1,1)+htau**2/2d0
      lpptata_int(1,1,2,2)=lpptata_int(1,1,2,2)+htau**2/2d0

      lpptata_int(2,3,1,2)=lpptata_int(2,3,1,2)+htau*ll/4d0
      lpptata_int(2,3,2,1)=lpptata_int(2,3,1,2)
      lpptata_int(3,2,1,2)=lpptata_int(2,3,1,2)
      lpptata_int(3,2,2,1)=lpptata_int(2,3,1,2)

*     trilinear, staus

      do i=1,3
       do k=1,2
          do l=1,2
             lptata_int(i,k,l)=0d0
          enddo
       enddo
      enddo

      lptata_int(1,1,2)=htau*Atau/sq2
      lptata_int(2,1,2)=htau*ll*xx/sq2
      lptata_int(3,1,2)=htau*ll*v2/sq2

      lptata_int(1,2,1)=-lptata_int(1,1,2)
      lptata_int(2,2,1)=-lptata_int(2,1,2)
      lptata_int(3,2,1)=-lptata_int(3,1,2)

*     sneutrinos (same as first two generations)

      do i=1,3
       lpntnt(i)=lpnn(i)
       do j=1,3
          lppntnt(i,j)=lppnn(i,j)
       enddo
      enddo

*     now rotate the third-generation couplings

*     quartic
      do i=1,3
       do j=1,3
          do k=1,2
             do l=1,2
              lpptt(i,j,k,l)=0d0
              lppbb(i,j,k,l)=0d0
              lpptata(i,j,k,l)=0d0
              do a=1,2
                 do b=1,2
                  lpptt(i,j,k,l)=lpptt(i,j,k,l)
     .                   +Rt(k,a)*Rt(l,b)*lpptt_int(i,j,a,b)
                  lppbb(i,j,k,l)=lppbb(i,j,k,l)
     .                   +Rb(k,a)*Rb(l,b)*lppbb_int(i,j,a,b)
                  lpptata(i,j,k,l)=lpptata(i,j,k,l)
     .                   +Rtau(k,a)*Rtau(l,b)*lpptata_int(i,j,a,b)
                 enddo
              enddo
             enddo
          enddo
       enddo
      enddo

*     trilinear
      do i=1,3
       do k=1,2
          do l=1,2
             lptt(i,k,l)=0d0
             lpbb(i,k,l)=0d0
             lptata(i,k,l)=0d0
             do a=1,2
              do b=1,2
                 lptt(i,k,l)=lptt(i,k,l)
     .                +Rt(k,a)*Rt(l,b)*lptt_int(i,a,b)
                 lpbb(i,k,l)=lpbb(i,k,l)
     .                +Rb(k,a)*Rb(l,b)*lpbb_int(i,a,b)
                 lptata(i,k,l)=lptata(i,k,l)
     .                +Rtau(k,a)*Rtau(l,b)*lptata_int(i,a,b)
              enddo
             enddo
          enddo
       enddo
      enddo

      end

*
***********************************************************************
*

      SUBROUTINE coupl_p_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
     .     lpphh,lppaa,lpah,lppcc,lpcc)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ll,kk,v1,v2,xx,Al,Ak,RS(3,3),RP(3,3),
     .     lpphh(3,3,3,3),lppaa(3,3,3,3),lpah(3,3,3),
     .     lppcc(3,3,2,2),lpcc(3,2,2)

      DOUBLE PRECISION lpppp(3,3,3,3),lsspp(3,3,3,3),
     .     lspp(3,3,3),gb2,sq2,c2b,s2b

      INTEGER i,j,k,l,a,b

      gb2=(g**2+gp**2)/2d0
      sq2=sqrt(2d0)

      do i=1,3              ! initialize
       do j=1,3
          do k=1,3
             lspp(i,j,k)=0d0
             lpah(i,j,k)=0d0
             do l=1,3
              lpppp(i,j,k,l)=0d0
              lsspp(i,j,k,l)=0d0
              lpphh(i,j,k,l)=0d0
              lppaa(i,j,k,l)=0d0
             enddo
          enddo
          do k=1,2
             do l=1,2
              lppcc(i,j,k,l)=0d0
             enddo
          enddo
       enddo
       do j=1,2
          do k=1,2
             lpcc(i,j,k)=0d0
          enddo
       enddo
      enddo

*     quartic neutral couplings

      lpppp(1,1,1,1)=gb2/16d0

      lpppp(2,2,2,2)=gb2/16d0

      lpppp(3,3,3,3)=kk**2/4d0

      lpppp(1,1,2,2)=(2d0*ll**2-gb2)/48d0
      lpppp(1,2,1,2)=lpppp(1,1,2,2) ! symmetrize the couplings
      lpppp(1,2,2,1)=lpppp(1,1,2,2) ! (got a better proposal? ;-)
      lpppp(2,1,1,2)=lpppp(1,1,2,2)
      lpppp(2,1,2,1)=lpppp(1,1,2,2)
      lpppp(2,2,1,1)=lpppp(1,1,2,2)

      lpppp(1,2,3,3)=-ll*kk/24d0
      lpppp(2,1,3,3)=lpppp(1,2,3,3)
      lpppp(2,3,1,3)=lpppp(1,2,3,3)
      lpppp(1,3,2,3)=lpppp(1,2,3,3)
      lpppp(3,1,2,3)=lpppp(1,2,3,3)
      lpppp(3,2,1,3)=lpppp(1,2,3,3)
      lpppp(3,1,3,2)=lpppp(1,2,3,3)
      lpppp(3,2,3,1)=lpppp(1,2,3,3)
      lpppp(3,3,1,2)=lpppp(1,2,3,3)
      lpppp(3,3,2,1)=lpppp(1,2,3,3)
      lpppp(1,3,3,2)=lpppp(1,2,3,3)
      lpppp(2,3,3,1)=lpppp(1,2,3,3)

      lpppp(1,1,3,3)=ll**2/24d0
      lpppp(1,3,1,3)=lpppp(1,1,3,3)
      lpppp(1,3,3,1)=lpppp(1,1,3,3)
      lpppp(3,1,1,3)=lpppp(1,1,3,3)
      lpppp(3,1,3,1)=lpppp(1,1,3,3)
      lpppp(3,3,1,1)=lpppp(1,1,3,3)
      lpppp(2,2,3,3)=lpppp(1,1,3,3)
      lpppp(2,3,2,3)=lpppp(1,1,3,3)
      lpppp(2,3,3,2)=lpppp(1,1,3,3)
      lpppp(3,2,2,3)=lpppp(1,1,3,3)
      lpppp(3,2,3,2)=lpppp(1,1,3,3)
      lpppp(3,3,2,2)=lpppp(1,1,3,3)

      lsspp(1,1,1,1)=gb2/8d0
      lsspp(2,2,2,2)=lsspp(1,1,1,1)

      lsspp(1,1,2,2)=(2d0*ll**2-gb2)/8d0
      lsspp(2,2,1,1)=lsspp(1,1,2,2)

      lsspp(1,1,3,3)=ll**2/4d0
      lsspp(2,2,3,3)=lsspp(1,1,3,3)
      lsspp(3,3,1,1)=lsspp(1,1,3,3)
      lsspp(3,3,2,2)=lsspp(1,1,3,3)

      lsspp(1,2,3,3)=ll*kk/4d0
      lsspp(2,1,3,3)=lsspp(1,2,3,3)
      lsspp(3,3,1,2)=lsspp(1,2,3,3)
      lsspp(3,3,2,1)=lsspp(1,2,3,3)

      lsspp(1,3,2,3)=-ll*kk/4d0
      lsspp(3,1,2,3)=lsspp(1,3,2,3)
      lsspp(1,3,3,2)=lsspp(1,3,2,3)
      lsspp(3,1,3,2)=lsspp(1,3,2,3)
      lsspp(2,3,1,3)=lsspp(1,3,2,3)
      lsspp(3,2,1,3)=lsspp(1,3,2,3)
      lsspp(2,3,3,1)=lsspp(1,3,2,3)
      lsspp(3,2,3,1)=lsspp(1,3,2,3)

      lsspp(3,3,3,3)=kk**2/2d0

*     rotate the quartics

      do i=1,3
       do j=1,3
          do k=1,3
             do l=1,3
              do a=1,3
                 do b=1,3
                  lpphh(i,j,k,l)=lpphh(i,j,k,l)
     .                   +RS(k,a)*RS(l,b)*lsspp(a,b,i,j)
                  lppaa(i,j,k,l)=lppaa(i,j,k,l)
     .                   +6d0*RP(k,a)*RP(l,b)*lpppp(i,j,a,b)
                 enddo
              enddo
             enddo
          enddo
       enddo
      enddo

*     trilinear neutral couplings

      lspp(1,1,1)=gb2*v1/2d0/sq2

      lspp(2,2,2)=gb2*v2/2d0/sq2

      lspp(1,2,2)=v1/2d0/sq2*(2d0*ll**2-gb2)

      lspp(2,1,1)=v2/2d0/sq2*(2d0*ll**2-gb2)

      lspp(3,1,1)=ll**2*xx/sq2
      lspp(3,2,2)=lspp(3,1,1)

      lspp(3,3,3)=-kk*Ak/sq2+sq2*kk**2*xx

      lspp(1,3,3)=ll/sq2*(ll*v1+kk*v2)

      lspp(2,3,3)=ll/sq2*(ll*v2+kk*v1)

      lspp(3,1,3)=-ll*kk*v2/sq2
      lspp(3,3,1)=lspp(3,1,3)

      lspp(3,2,3)=-ll*kk*v1/sq2
      lspp(3,3,2)=lspp(3,2,3)

      lspp(1,2,3)=ll*Al/2d0/sq2-ll*kk*xx/sq2
      lspp(1,3,2)=lspp(1,2,3)
      lspp(2,1,3)=lspp(1,2,3)
      lspp(2,3,1)=lspp(1,2,3)

      lspp(3,1,2)=ll*Al/2d0/sq2+ll*kk*xx/sq2
      lspp(3,2,1)=lspp(3,1,2)

*     rotate the trilinears

      do i=1,3
       do k=1,3
          do l=1,3
             do a=1,3
              do b=1,3
                 lpah(i,k,l)=lpah(i,k,l)
     .                +2d0*RS(l,a)*RP(k,b)*lspp(a,b,i)
              enddo
             enddo
          enddo
       enddo
      enddo

*     quartic charged couplings

      c2b=(v1**2-v2**2)/(v1**2+v2**2)
      s2b=2d0*v1*v2/(v1**2+v2**2)

      lppcc(1,1,1,1)=(g**2+gp**2*c2b)/8d0
      lppcc(2,2,2,2)=lppcc(1,1,1,1)

      lppcc(1,1,2,2)=(g**2-gp**2*c2b)/8d0
      lppcc(2,2,1,1)=lppcc(1,1,2,2)

      lppcc(1,2,1,1)=-(2d0*ll**2-g**2)*s2b/8d0
      lppcc(2,1,1,1)=lppcc(1,2,1,1)

      lppcc(1,2,2,2)=(2d0*ll**2-g**2)*s2b/8d0
      lppcc(2,1,2,2)=lppcc(1,2,2,2)

      lppcc(1,1,1,2)=-gp**2/8d0*s2b
      lppcc(1,1,2,1)=lppcc(1,1,1,2)

      lppcc(2,2,1,2)=gp**2/8d0*s2b
      lppcc(2,2,2,1)=lppcc(2,2,1,2)

      lppcc(1,2,1,2)=-(2d0*ll**2-g**2)*s2b/8d0
      lppcc(1,2,2,1)=lppcc(1,2,1,2)
      lppcc(2,1,1,2)=lppcc(1,2,1,2)
      lppcc(2,1,2,1)=lppcc(1,2,1,2)

      lppcc(3,3,1,2)=kk*ll/2d0*c2b
      lppcc(3,3,2,1)=lppcc(3,3,1,2)

      lppcc(3,3,1,1)=ll/2d0*(ll+kk*s2b)

      lppcc(3,3,2,2)=ll/2d0*(ll-kk*s2b)

*     trilinear charged couplings

      lpcc(1,1,2)=v2*(2d0*ll**2-g**2)/2d0/sq2
      lpcc(1,2,1)=-lpcc(1,1,2)

      lpcc(2,1,2)=v1*(2d0*ll**2-g**2)/2d0/sq2
      lpcc(2,2,1)=-lpcc(2,1,2)

      lpcc(3,1,2)=ll/sq2*(Al-2d0*kk*xx)
      lpcc(3,2,1)=-lpcc(3,1,2)

      end

*
***********************************************************************
*

      SUBROUTINE coupl_p_ino(g,gp,ll,kk,NN,UU,VV,lpnene,lpchch)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ll,kk,NN(5,5),UU(2,2),VV(2,2),
     .     lpnene(3,5,5),lpchch(3,2,2)

      DOUBLE PRECISION lp00(3,5,5),lppm(3,2,2),sq2

      INTEGER i,k,l,a,b

      sq2=sqrt(2d0)

*     neutralino couplings

      do i=1,3
       do k=1,5
          do l=1,5
             lp00(i,k,l)=0d0
             lpnene(i,k,l)=0d0
          enddo
       enddo
      enddo

      lp00(1,1,3)=gp/4d0
      lp00(1,3,1)=lp00(1,1,3)

      lp00(2,1,4)=-gp/4d0
      lp00(2,4,1)=lp00(2,1,4)

      lp00(1,2,3)=-g/4d0
      lp00(1,3,2)=lp00(1,2,3)

      lp00(2,2,4)=g/4d0
      lp00(2,4,2)=lp00(2,2,4)

      lp00(3,5,5)=kk/sq2

      lp00(1,4,5)=-ll/2d0/sq2
      lp00(1,5,4)=lp00(1,4,5)
      lp00(2,3,5)=lp00(1,4,5)
      lp00(2,5,3)=lp00(1,4,5)
      lp00(3,3,4)=lp00(1,4,5)
      lp00(3,4,3)=lp00(1,4,5)

      do i=1,3
       do k=1,5
          do l=1,5
             do a=1,5
              do b=1,5
                 lpnene(i,k,l)=lpnene(i,k,l)
     .                +NN(k,a)*NN(l,b)*lp00(i,a,b)
              enddo
             enddo
          enddo
       enddo
      enddo

*     chargino couplings

      do i=1,3
       do k=1,2
          do l=1,2
             lppm(i,k,l)=0d0
             lpchch(i,k,l)=0d0
          enddo
       enddo
      enddo

      lppm(1,1,2)=-g/sq2

      lppm(2,2,1)=-g/sq2

      lppm(3,2,2)=ll/sq2

      do i=1,3
       do k=1,2
          do l=1,2
             do a=1,2
              do b=1,2
                 lpchch(i,k,l)=lpchch(i,k,l)
     .                +VV(k,a)*UU(l,b)*lppm(i,a,b)
              enddo
             enddo
          enddo
       enddo
      enddo

      end

*
***********************************************************************
*

      SUBROUTINE coupl_Z_ino(g,gp,NN,UU,VV,lznene,azchch,bzchch)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,NN(5,5),UU(2,2),VV(2,2),
     .     lznene(5,5),azchch(2,2),bzchch(2,2)

      DOUBLE PRECISION lz00(5,5),lzpm(2,2),sq2,gb2,c2w

      INTEGER i,j,k,l

      sq2=sqrt(2d0)
      gb2=(g**2+gp**2)/2d0
      c2w=(g**2-gp**2)/(g**2+gp**2)

*     neutralino couplings

      do i=1,5
       do j=1,5
          lz00(i,j)=0d0
          lznene(i,j)=0d0
       enddo
      enddo

      lz00(3,3)=sqrt(gb2)/2d0/sq2

      lz00(4,4)=-lz00(3,3)

      do i=1,5
       do j=1,5
          do k=1,5
             do l=1,5
              lznene(i,j)=lznene(i,j)
     .             +NN(i,k)*NN(j,l)*lz00(k,l)
             enddo
          enddo
       enddo
      enddo

*     chargino couplings

      do i=1,2
       do j=1,2
          lzpm(i,j)=0d0
          azchch(i,j)=0d0
          bzchch(i,j)=0d0
       enddo
      enddo

      lzpm(1,1)=g**2/sqrt(gb2)/sq2

      lzpm(2,2)=sqrt(gb2)/sq2*c2w

      do i=1,2
       do j=1,2
          do k=1,2
             do l=1,2
              azchch(i,j)=azchch(i,j)
     .             +VV(i,k)*VV(j,l)*lzpm(k,l)
              bzchch(i,j)=bzchch(i,j)
     .             +UU(i,k)*UU(j,l)*lzpm(k,l)
             enddo
          enddo
       enddo
      enddo

      end

*
***********************************************************************
*

      SUBROUTINE coupl_W_ino(g,NN,UU,VV,awnech,bwnech)

      IMPLICIT NONE

      DOUBLE PRECISION g,NN(5,5),UU(2,2),VV(2,2),
     .     awnech(5,2),bwnech(5,2)

      DOUBLE PRECISION aw0p(5,2),bw0p(5,2),sq2

      INTEGER i,j,k,l

      sq2=sqrt(2d0)

      do i=1,5
       do j=1,2
          aw0p(i,j)=0d0
          bw0p(i,j)=0d0
          awnech(i,j)=0d0
          bwnech(i,j)=0d0
       enddo
      enddo

      aw0p(2,1)=-g
      aw0p(4,2)=g/sq2

      bw0p(2,1)=-g
      bw0p(3,2)=-g/sq2

      do i=1,5
       do j=1,2
          do k=1,5
             do l=1,2
              awnech(i,j)=awnech(i,j)
     .             +NN(i,k)*VV(j,l)*aw0p(k,l)
              bwnech(i,j)=bwnech(i,j)
     .             +NN(i,k)*UU(j,l)*bw0p(k,l)
             enddo
          enddo
       enddo
      enddo

      end

*
***********************************************************************
*

      SUBROUTINE twlpyuk(mt,mb,A0,T1,T2,B1,B2,st,ct,sb,cb,q,l,xx,tanb,
     .     vv,DMS,DMP)

      IMPLICIT NONE

      INTEGER i,j
      DOUBLE PRECISION mt,mb,A0,T1,T2,B1,B2,st,ct,sb,cb,q,l,xx,tanb,vv,
     .     DMS(3,3),DMP(3,3)
      DOUBLE PRECISION c2t,s2t,c2b,s2b,At,Ab,Xt,Xb,t,b,cbe,sbe,ht,hb,
     .     pi,k,mu,DMA
      DOUBLE PRECISION F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,FA!1,FA2,FA3
      DOUBLE PRECISION DDCOS,DDSIN

      pi = 4d0*atan(1d0)

      t = mt**2
      b = mb**2
      mu = l*xx

      s2t = 2d0*ct*st
      s2b = 2d0*cb*sb
      c2t = ct**2-st**2
      c2b = cb**2-sb**2

      Xt = (T1-T2)*s2t/2d0/mt
      Xb = (B1-B2)*s2b/2d0/mb
      At = Xt+mu/tanb
      Ab = Xb+mu*tanb

      sbe = DDSIN(datan(tanb))
      cbe = DDCOS(datan(tanb))

      ht = mt/vv/sbe
      hb = mb/vv/cbe

      k = 3d0/(16d0*Pi**2)**2


      CALL makefuncstb(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb,
     .     F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,FA)

      DMS(1,1) = .5d0*ht**2*mu**2*s2t**2*F3t
     .    +2d0*hb**2*mb**2*F1b+2d0*hb**2*Ab*mb*s2b*F2b
     .    +.5d0*hb**2*Ab**2*s2b**2*F3b
     .    -2d0*hb*ht*mb*mu*s2t*F4t-ht*hb*mu*Ab*s2t*s2b*F5
     .    +ht**2*tanb*mu*At/(T1-T2)*Ft+hb**2*tanb*mu*Ab/(B1-B2)*Fb

      DMS(1,2) = -ht**2*mu*mt*s2t*F2t-.5d0*ht**2*At*mu*s2t**2*F3t
     .    -hb**2*mu*mb*s2b*F2b-.5d0*hb**2*Ab*mu*s2b**2*F3b
     .    +ht*hb*mb*At*s2t*F4t+hb*ht*mt*Ab*s2b*F4b
     .    +.5d0*ht*hb*s2t*s2b*(At*Ab+mu**2)*F5
     .    +2d0*ht*hb*mt*mb*F6
     .    -ht**2*mu*At/(T1-T2)*Ft-hb**2*mu*Ab/(B1-B2)*Fb

      DMS(2,2) = .5d0*hb**2*mu**2*s2b**2*F3b
     .    +2d0*ht**2*mt**2*F1t+2d0*ht**2*At*mt*s2t*F2t
     .    +.5d0*ht**2*At**2*s2t**2*F3t
     .    -2d0*ht*hb*mt*mu*s2b*F4b-hb*ht*mu*At*s2b*s2t*F5
     .    +ht**2/tanb*mu*At/(T1-T2)*Ft+hb**2/tanb*mu*Ab/(B1-B2)*Fb

      DMS(1,3) = 0d0
*     .     .5d0*ht*s2t**2*l*mu*mt/tanb*F3t
*     .    -ht*l*mt*(At-2d0*mu/tanb)/(T1-T2)*Ft
*     .    -.5d0*hb*s2b**2*l*Ab*mb*tanb*F3b
*     .    -hb*s2b*l*b*tanb*F2b-hb*l*mb*Ab*tanb/(B1-B2)*Fb
*     .    -l*ht*b*s2t*F4t-.5d0*ht*l*mb*s2t*s2b*(Ab-mu*tanb)*F5

      DMS(2,3) = 0d0
*     .    -.5d0*ht*s2t**2*l*At*mt/tanb*F3t-ht*s2t*l*t/tanb*F2t
*     .    -ht*l*mt*At/tanb/(T1-T2)*Ft
*     .    +.5d0*hb*s2b**2*l*mu*mb*tanb*F3b
*     .    -hb*l*mb*(Ab-2d0*mu*tanb)/(B1-B2)*Fb
*     .    -l*hb*t*s2b*F4b-.5d0*hb*l*mt*s2t*s2b*(At-mu/tanb)*F5

      DMS(3,3) = 0d0
*     .     .5d0*l**2*s2t**2*t/tanb**2*F3t
*     .    +l**2*t/tanb*At/mu/(T1-T2)*Ft
*     .    +.5d0*l**2*s2b**2*b*tanb**2*F3b
*     .    +l**2*b*tanb*Ab/mu/(B1-B2)*Fb
*     .    +l**2*mb*mt*s2t*s2b*F5

      DMS(2,1) = DMS(1,2)
      DMS(3,1) = DMS(1,3)
      DMS(3,2) = DMS(2,3)

      DMA = ht**2*mu*At/(T1-T2)*Ft+hb**2*mu*Ab/(B1-B2)*Fb -2d0*ht*hb*FA

      DMP(1,1) = DMA*tanb

      DMP(1,2) = DMA

      DMP(2,2) = DMA/tanb

      DMP(1,3) = 0d0
*     .     l*ht*mt*At/(T1-T2)*Ft
*     .    +l*hb*mb*Ab/(B1-B2)*Fb*tanb
*     .    -2d0*l*hb*mt*FA2

      DMP(2,3) = 0d0
*     .     l*ht*mt*At/(T1-T2)*Ft/tanb
*     .    +l*hb*mb*Ab/(B1-B2)*Fb
*     .    -2d0*l*ht*mb*FA2

      DMP(3,3) = 0d0
*     .     l**2*t*At/mu/(T1-T2)*Ft/tanb
*     .    +l**2*b*Ab/mu/(B1-B2)*Fb*tanb
*     .    -l**2*FA3

      DMP(2,1) = DMP(1,2)
      DMP(3,1) = DMP(1,3)
      DMP(3,2) = DMP(2,3)

      do i=1,3
       do j=1,3
        DMS(i,j) = k*DMS(i,j)
        DMP(i,j) = k*DMP(i,j)
       enddo
      enddo

      end

*
***********************************************************************
*

      SUBROUTINE makefuncstb(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,
     .     q,mu,vv,tanb,F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,
     .     Ft,Fb,FA)

      IMPLICIT NONE

      DOUBLE PRECISION t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb,
     .     F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,FA!1,FA2,FA3

      DOUBLE PRECISION D1t,DT1,DT2,Dc2t,DT1T1,DT2T2,Dtt,Dc2tc2t,DT1t,
     .     DT2t,DT1T2,Dtc2t,DT1c2t,DT2c2t,Dtb,DT1b,DT2b,DB1t,DB2t,DT1B1,
     .     DT2B1,DT1B2,DT2B2,Dbc2t,DB1c2t,DB2c2t,DT1c2b,DT2c2b,Dc2tc2b,
     .     Dcptpb,Dcpttptb,Dcpbptt,Dcptptb,Dcptmptt,Dcpbmptb,
     .     Dspbmptbspbptt,Dsptmpttsptptb,Dsptmpttspbmptb

      COMMON/listderivtb/D1t,DT1,DT2,Dc2t,DT1T1,DT2T2,
     .     Dtt,Dc2tc2t,DT1t,DT2t,DT1T2,
     .     Dtc2t,DT1c2t,DT2c2t,Dtb,DT1b,DT2b,DB1t,DB2t,DT1B1,DT2B1,
     .     DT1B2,DT2B2,Dbc2t,DB1c2t,DB2c2t,DT1c2b,DT2c2b,Dc2tc2b,
     .     Dcptpb,Dcpttptb,Dcpbptt,Dcptptb,Dcptmptt,Dcpbmptb,
     .     Dspbmptbspbptt,Dsptmpttsptptb,Dsptmpttspbmptb

      DOUBLE PRECISION D1b,DB1,DB2,Dc2b,DB1B1,DB2B2,Dbb,Dc2bc2b,DB1b,
     .     DB2b,DB1B2,Dbc2b,DB1c2b,DB2c2b,Dtc2b

      DOUBLE PRECISION Xt,Xb,At,Ab


      CALL makederivtb(b,t,A0,B1,B2,T1,T2,s2b,c2b,s2t,c2t,
     .     q,-mu,vv,1d0/tanb)
*     note: the potential was computed with the opposite convention for mu

      D1b = D1t
      DB1 = DT1
      DB2 = DT2
      Dc2b = Dc2t
      DB1B1 = DT1T1
      DB2B2 = DT2T2
      Dbb = Dtt
      Dc2bc2b = Dc2tc2t
      DB1b = DT1t
      DB2b = DT2t
      DB1B2 = DT1T2
      Dbc2b = Dtc2t
      DB1c2b = DT1c2t
      DB2c2b = DT2c2t
      Dtc2b = Dbc2t


      CALL makederivtb(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,-mu,vv,tanb)

      F1t = Dtt+DT1T1+DT2T2+2d0*(DT1t+DT2t+DT1T2)

      F2t = DT1T1-DT2T2+DT1t-DT2t
     .    -4d0*c2t**2/(T1-T2)*(Dtc2t+DT1c2t+DT2c2t)

      F3t = DT1T1+DT2T2-2d0*DT1T2
     .    -2d0/(T1-T2)*(DT1-DT2)
     .    +16d0*c2t**2/(T1-T2)**2*(c2t**2*Dc2tc2t+2d0*Dc2t)
     .    -8d0*c2t**2/(T1-T2)*(DT1c2t-DT2c2t)

      F4t = DT1b+DT1B1+DT1B2-DT2b-DT2B1-DT2B2
     .    -4d0*c2t**2/(T1-T2)*(DB1c2t+DB2c2t+Dbc2t)

      F1b = Dbb+DB1B1+DB2B2+2d0*(DB1b+DB2b+DB1B2)

      F2b = DB1B1-DB2B2+DB1b-DB2b
     .    -4d0*c2b**2/(B1-B2)*(Dbc2b+DB1c2b+DB2c2b)

      F3b = DB1B1+DB2B2-2d0*DB1B2
     .    -2d0/(B1-B2)*(DB1-DB2)
     .    +16d0*c2b**2/(B1-B2)**2*(c2b**2*Dc2bc2b+2d0*Dc2b)
     .    -8d0*c2b**2/(B1-B2)*(DB1c2b-DB2c2b)

      F4b = DB1t+DT1B1-DT1B2-DB2t+DT2B1-DT2B2
     .    -4d0*c2b**2/(B1-B2)*(DT1c2b+DT2c2b+Dtc2b)

      F5  = DT1B1-DT1B2-DT2B1+DT2B2
     .    +16d0*c2t**2*c2b**2/(T1-T2)/(B1-B2)*Dc2tc2b
     .    -4d0*c2t**2/(T1-T2)*(DB1c2t-DB2c2t)
     .    -4d0*c2b**2/(B1-B2)*(DT1c2b-DT2c2b)

      F6 = Dtb+DT1b+DT2b+DB1t+DB2t
     .    +DT1B1+DT1B2+DT2B1+DT2B2

      Ft = DT1-DT2-4d0*c2t**2/(T1-T2)*Dc2t

      Fb = DB1-DB2-4d0*c2b**2/(B1-B2)*Dc2b

      Xt = (T1-T2)*s2t/2d0/sqrt(t)
      Xb = (B1-B2)*s2b/2d0/sqrt(b)

      At = Xt+mu/tanb
      Ab = Xb+mu*tanb

      FA = Dcptpb/4d0/Sqrt(b)/Sqrt(t)
     .    +4d0*Sqrt(t)*Sqrt(b)*(At*Ab-mu**2)**2
     .     /s2t**2/s2b**2/(T1-T2)**2/(B1-B2)**2*Dcpttptb
     .    +Sqrt(t)/Sqrt(b)/s2t**2/(T1-T2)**2
     .     *(At**2*Dcpbptt+mu**2/tanb**2*Dcptmptt)
     .    +Sqrt(b)/Sqrt(t)/s2b**2/(B1-B2)**2
     .     *(Ab**2*Dcptptb+mu**2*tanb**2*Dcpbmptb)
     .    +2d0*mu/s2t/s2b/(T1-T2)/(B1-B2)
     .     *(At*tanb*Dspbmptbspbptt+Ab/tanb*Dsptmpttsptptb
     .    -mu*Dsptmpttspbmptb)

*      FA2 = 4d0*Sqrt(t)*Sqrt(b)*(At*Ab-mu**2)*(Ab/tanb+At*tanb-2d0*mu)
*     .     /s2t**2/s2b**2/(T1-T2)**2/(B1-B2)**2*Dcpttptb
*     .    +Sqrt(t)/Sqrt(b)/tanb/s2t**2/(T1-T2)**2
*     .     *(At*Dcpbptt+mu/tanb*Dcptmptt)
*     .    +Sqrt(b)/Sqrt(t)*tanb/s2b**2/(B1-B2)**2
*     .     *(Ab*Dcptptb+mu*tanb*Dcpbmptb)
*     .    +1d0/s2t/s2b/(T1-T2)/(B1-B2)
*     .     *((At*tanb+mu)*Dspbmptbspbptt+(Ab/tanb+mu)*Dsptmpttsptptb
*     .    -2d0*mu*Dsptmpttspbmptb)
*
*      FA3 = 8d0*t*b*(Ab/tanb+At*tanb-2d0*mu)**2
*     .     /s2t**2/s2b**2/(T1-T2)**2/(B1-B2)**2*Dcpttptb
*     .    +2d0*t/tanb**2/s2t**2/(T1-T2)**2*(Dcpbptt+Dcptmptt)
*     .    +2d0*b*tanb**2/s2b**2/(B1-B2)**2*(Dcptptb+Dcpbmptb)
*     .    +4d0*Sqrt(t)*Sqrt(b)/s2t/s2b/(T1-T2)/(B1-B2)
*     .     *(Dspbmptbspbptt+Dsptmpttsptptb-Dsptmpttspbmptb)

      end

*
***********************************************************************
*

      SUBROUTINE makederivtb(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,
     .     q,mu,vv,tanb)

      IMPLICIT DOUBLE PRECISION (t)
      IMPLICIT CHARACTER (a-s,u-z)

      DOUBLE PRECISION t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb
      DOUBLE PRECISION ht,hb,mt,mb,Xt,Xb,Yt,Yb,sbe,cbe,mu2,Nc
      DOUBLE PRECISION Delt,phi,SLLi2
      EXTERNAL delt,phi,SLLi2

      DOUBLE PRECISION D1t,DT1,DT2,Dc2t,DT1T1,DT2T2,Dtt,Dc2tc2t,DT1t,
     .     DT2t,DT1T2,Dtc2t,DT1c2t,DT2c2t,Dtb,DT1b,DT2b,DB1t,DB2t,DT1B1,
     .     DT2B1,DT1B2,DT2B2,Dbc2t,DB1c2t,DB2c2t,DT1c2b,DT2c2b,Dc2tc2b,
     .     Dcptpb,Dcpttptb,Dcpbptt,Dcptptb,Dcptmptt,Dcpbmptb,
     .     Dspbmptbspbptt,Dsptmpttsptptb,Dsptmpttspbmptb

      COMMON/listderivtb/D1t,DT1,DT2,Dc2t,DT1T1,DT2T2,
     .     Dtt,Dc2tc2t,DT1t,DT2t,DT1T2,
     .     Dtc2t,DT1c2t,DT2c2t,Dtb,DT1b,DT2b,DB1t,DB2t,DT1B1,DT2B1,
     .     DT1B2,DT2B2,Dbc2t,DB1c2t,DB2c2t,DT1c2b,DT2c2b,Dc2tc2b,
     .     Dcptpb,Dcpttptb,Dcpbptt,Dcptptb,Dcptmptt,Dcpbmptb,
     .     Dspbmptbspbptt,Dsptmpttsptptb,Dsptmpttspbmptb

      DOUBLE PRECISION  Logt,Logb,Logmu2,LogA0,LogT1,LogT2,LogB1,LogB2,
     .     phimu2tT1,phimu2T2T,phiB1tmu2,phiB2tmu2,phiT1bmu2,phiT2bmu2,
     .     phiA0T1T1,phiA0T2T2,phiA0T1T2,phiA0B1B1,phiA0B2B2,
     .     phiA0B1T1,phiA0B2T1,phiA0B1T2,phiA0B2T2,phiA0tt,phiA0bt,
     .     deltmu2tT1,deltmu2T2T,deltB1tmu2,deltB2tmu2,
     .     deltT1bmu2,delT2Tbmu2,deltA0tt,deltA0bt,
     .     deltA0T1T1,deltA0T2T2,deltA0T1T2,deltA0B1B1,deltA0B2B2,
     .     deltA0B1T1,deltA0B2T1,deltA0B1T2,deltA0B2T2,
     .     Li2T1T2,Li2T1B1,Li2T1B2,Li2T2B1,Li2T2B2,Li2bt

      DOUBLE PRECISION DDCOS,DDSIN

      Nc = 3d0

      mt = dsqrt(t)
      mb = dsqrt(b)

      sbe = DDSIN(datan(tanb))
      cbe = DDCOS(datan(tanb))

      ht = mt/vv/sbe
      hb = mb/vv/cbe

      Xt = (T1-T2)*s2t/2d0/mt
      Xb = (B1-B2)*s2b/2d0/mb
      Yt  = Xt-mu/sbe/cbe
      Yb  = Xb-mu/sbe/cbe

      mu2 = mu**2

      Logt = Log(t/q)
      Logb = Log(b/q)
      Logmu2 = Log(mu2/q)
      LogA0 = Log(A0/q)
      LogT1 = Log(T1/q)
      LogT2 = Log(T2/q)
      LogB1 = Log(B1/q)
      LogB2 = Log(B2/q)
      phimu2tT1 = phi(mu2,t,T1)
      phimu2T2T = phi(mu2,t,T2)
      phiB1tmu2 = phi(B1,t,mu2)
      phiB2tmu2 = phi(B2,t,mu2)
      phiT1bmu2 = phi(T1,b,mu2)
      phiT2bmu2 = phi(T2,b,mu2)
      phiA0T1T1 = phi(A0,T1,T1)
      phiA0T1T2 = phi(A0,T1,T2)
      phiA0T2T2 = phi(A0,T2,T2)
      phiA0B1B1 = phi(A0,B1,B1)
      phiA0B2B2 = phi(A0,B2,B2)
      phiA0B1T1 = phi(A0,B1,T1)
      phiA0B1T2 = phi(A0,B1,T2)
      phiA0B2T1 = phi(A0,B2,T1)
      phiA0B2T2 = phi(A0,B2,T2)
      phiA0tt = phi(A0,t,t)
      phiA0bt = phi(A0,b,t)
      deltmu2tT1 = delt(mu2,t,T1)
      deltmu2T2T = delt(mu2,t,T2)
      deltB1tmu2 = delt(B1,t,mu2)
      deltB2tmu2 = delt(B2,t,mu2)
      deltT1bmu2 = delt(T1,b,mu2)
      delT2Tbmu2 = delt(T2,b,mu2)
      deltA0T1T1 = delt(A0,T1,T1)
      deltA0T1T2 = delt(A0,T1,T2)
      deltA0T2T2 = delt(A0,T2,T2)
      deltA0B1B1 = delt(A0,B1,B1)
      deltA0B2B2 = delt(A0,B2,B2)
      deltA0B1T1 = delt(A0,B1,T1)
      deltA0B1T2 = delt(A0,B1,T2)
      deltA0B2T1 = delt(A0,B2,T1)
      deltA0B2T2 = delt(A0,B2,T2)
      deltA0tt = delt(A0,t,t)
      deltA0bt = delt(A0,b,t)
      Li2T1T2 = SLLi2(1d0-T1/T2)
      Li2T1B1 = SLLi2(1d0-T1/B1)
      Li2T2B1 = SLLi2(1d0-T2/B1)
      Li2T1B2 = SLLi2(1d0-T1/B2)
      Li2T2B2 = SLLi2(1d0-T2/B2)
      Li2bt = SLLi2(1d0-b/t)

      tmp1 = 2d0-LogA0-Logb+2d0*Logt
      tmp2 = 2d0-LogB1-Logmu2+2d0*Logt
      tmp3 = 2d0-LogB2-Logmu2+2d0*Logt
      tmp4 = 2d0-Logb-Logmu2+2d0*LogT1
      tmp5 = 2d0-Logb-Logmu2+2d0*LogT2
      tmp6 = 0.25d0*hb**2/c2t-0.25d0*ht**2/c2t
      tmp7 = -(0.25d0*hb**2/c2t)+0.25d0*ht**2/c2t
      tmp8 = c2t**2-0.5d0*((-1d0+Nc)*s2t**2)
      tmp9 = cbe**2*ht**2+hb**2*sbe**2
      tmp10 = cbe**2*hb**2+ht**2*sbe**2
      tmp11 = -(0.5d0*(cbe**2*ht**2*mb*mt*s2b*s2t))-
     .   0.5d0*(hb**2*mb*mt*s2b*s2t*sbe**2)
      tmp12 = 0.5d0*(cbe**2*ht**2*mb*mt*s2b*s2t)+
     .   0.5d0*(hb**2*mb*mt*s2b*s2t*sbe**2)
      tmp13 = -(0.5d0*(cbe**2*hb**2*mb*mt*s2b*s2t))-
     .   0.5d0*(ht**2*mb*mt*s2b*s2t*sbe**2)
      tmp14 = 0.5d0*(cbe**2*hb**2*mb*mt*s2b*s2t)+
     .   0.5d0*(ht**2*mb*mt*s2b*s2t*sbe**2)
      tmp15 = 2d0/t**2-(2d0*(-1d0+Logt))/t**2
      tmp16 = 2d0/t**2-(2d0*Logt)/t**2
      tmp17 = 2d0/T1**2-(2d0*LogT1)/T1**2
      tmp18 = 2d0/T2**2-(2d0*LogT2)/T2**2
      tmp19 = -(0.0625d0*
     .      (B1*(-1d0+LogB1)*(-1d0+LogT1)*T1)/(c2b*c2t))+
     .   0.0625d0*(B2*(-1d0+LogB2)*(-1d0+LogT1)*T1)/(c2b*c2t)+
     .   0.0625d0*(B1*(-1d0+LogB1)*(-1d0+LogT2)*T2)/(c2b*c2t)-
     .   0.0625d0*(B2*(-1d0+LogB2)*(-1d0+LogT2)*T2)/(c2b*c2t)
      tmp20 = 0.25d0*((1d0+c2b)*(1d0+c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0+c2b)*s2t*Xb)/mt
      tmp21 = 0.0625d0*b/(c2b*c2t)+0.125d0*(mb*mt)/(s2b*s2t)+
     .   0.0625d0*t/(c2b*c2t)-0.125d0*(mb*Xb)/(c2t*s2b)+
     .   0.125d0*(mt*Xb)/(c2b*s2t)-0.0625d0*Xb**2/(c2b*c2t)
      tmp22 = -(0.0625d0*b/(c2b*c2t))-
     .   0.125d0*(mb*mt)/(s2b*s2t)-0.0625d0*t/(c2b*c2t)+
     .   0.125d0*(mb*Xb)/(c2t*s2b)-0.125d0*(mt*Xb)/(c2b*s2t)+
     .   0.0625d0*Xb**2/(c2b*c2t)
      tmp23 = 0.25d0*Xb/c2b-0.25d0*Xt/c2b
      tmp24 = -(0.25d0*Xb/c2b)+0.25d0*Xt/c2b
      tmp25 = 0.125d0*Xb/c2t**3-0.125d0*Xt/c2t**3
      tmp26 = -(0.125d0*Xb/c2t**3)+0.125d0*Xt/c2t**3
      tmp27 = 0.25d0*Xb/c2t-0.25d0*Xt/c2t
      tmp28 = -(0.25d0*Xb/c2t)+0.25d0*Xt/c2t
      tmp29 = 0.25d0*((1d0+c2b)*(1d0+c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0+c2b)*s2t*Xt)/mt
      tmp30 = (-2d0*mt*Xt)/s2t-Xt**2
      tmp31 = (2d0*mt*Xt)/s2t-Xt**2
      tmp32 = 0.0625d0*b/(c2b*c2t)+0.125d0*(mb*mt)/(s2b*s2t)+
     .   0.0625d0*t/(c2b*c2t)+0.125d0*(mb*Xt)/(c2t*s2b)-
     .   0.125d0*(mt*Xt)/(c2b*s2t)-0.0625d0*Xt**2/(c2b*c2t)
      tmp33 = -(0.0625d0*b/(c2b*c2t))-
     .   0.125d0*(mb*mt)/(s2b*s2t)-0.0625d0*t/(c2b*c2t)-
     .   0.125d0*(mb*Xt)/(c2t*s2b)+0.125d0*(mt*Xt)/(c2b*s2t)+
     .   0.0625d0*Xt**2/(c2b*c2t)
      tmp34 = 4d0*t-4d0*mt*s2t*Xt+s2t**2*Xt**2
      tmp35 = 4d0*t+4d0*mt*s2t*Xt+s2t**2*Xt**2
      tmp36 = 0.25d0*((1d0+c2b)*(1d0+c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0+c2b)*s2t*Yb)/mt
      tmp37 = 0.0625d0*b/(c2b*c2t)+0.125d0*(mb*mt)/(s2b*s2t)+
     .   0.0625d0*t/(c2b*c2t)-0.125d0*(mb*Yb)/(c2t*s2b)+
     .   0.125d0*(mt*Yb)/(c2b*s2t)-0.0625d0*Yb**2/(c2b*c2t)
      tmp38 = -(0.0625d0*b/(c2b*c2t))-
     .   0.125d0*(mb*mt)/(s2b*s2t)-0.0625d0*t/(c2b*c2t)+
     .   0.125d0*(mb*Yb)/(c2t*s2b)-0.125d0*(mt*Yb)/(c2b*s2t)+
     .   0.0625d0*Yb**2/(c2b*c2t)
      tmp39 = 0.25d0*Yb/c2b-0.25d0*Yt/c2b
      tmp40 = -(0.25d0*Yb/c2b)+0.25d0*Yt/c2b
      tmp41 = 0.125d0*Yb/c2t**3-0.125d0*Yt/c2t**3
      tmp42 = -(0.125d0*Yb/c2t**3)+0.125d0*Yt/c2t**3
      tmp43 = 0.25d0*Yb/c2t-0.25d0*Yt/c2t
      tmp44 = -(0.25d0*Yb/c2t)+0.25d0*Yt/c2t
      tmp45 = 0.25d0*((1d0+c2b)*(1d0+c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0+c2b)*s2t*Yt)/mt
      tmp46 = (-2d0*mt*Yt)/s2t-Yt**2
      tmp47 = (2d0*mt*Yt)/s2t-Yt**2
      tmp48 = 0.0625d0*b/(c2b*c2t)+0.125d0*(mb*mt)/(s2b*s2t)+
     .   0.0625d0*t/(c2b*c2t)+0.125d0*(mb*Yt)/(c2t*s2b)-
     .   0.125d0*(mt*Yt)/(c2b*s2t)-0.0625d0*Yt**2/(c2b*c2t)
      tmp49 = -(0.0625d0*b/(c2b*c2t))-
     .   0.125d0*(mb*mt)/(s2b*s2t)-0.0625d0*t/(c2b*c2t)-
     .   0.125d0*(mb*Yt)/(c2t*s2b)+0.125d0*(mt*Yt)/(c2b*s2t)+
     .   0.0625d0*Yt**2/(c2b*c2t)
      tmp50 = 4d0*t-4d0*mt*s2t*Yt+s2t**2*Yt**2
      tmp51 = 4d0*t+4d0*mt*s2t*Yt+s2t**2*Yt**2
      tmp52 = 0.0625d0*(1d0-c2b)/c2t**3-0.0625d0*(1d0+c2b)/c2t**3
      tmp53 = -(0.0625d0*(1d0-c2b)/c2t**3)
     .      +0.0625d0*(1d0+c2b)/c2t**3
      tmp54 = 0.125d0*(1d0-c2b)/c2t-0.125d0*(1d0+c2b)/c2t
      tmp55 = -(0.125d0*(1d0-c2b)/c2t)+0.125d0*(1d0+c2b)/c2t
      tmp56 = 0.25d0*((1d0+c2b)*(1d0-c2t))
     .      +0.25d0*((1d0-c2b)*(1d0+c2t))
      tmp57 = 0.125d0*(1d0-c2t)/c2b-0.125d0*(1d0+c2t)/c2b
      tmp58 = -(0.125d0*(1d0-c2t)/c2b)+0.125d0*(1d0+c2t)/c2b
      tmp59 = 0.25d0*((1d0-c2b)*(1d0-c2t))
     .      +0.25d0*((1d0+c2b)*(1d0+c2t))
      tmp60 = 0.5d0*((1d0+c2b)*hb**2)+0.5d0*((1d0-c2b)*ht**2)
      tmp61 = 0.5d0*((1d0-c2b)*hb**2)+0.5d0*((1d0+c2b)*ht**2)
      tmp62 = 0.5d0*((1d0+c2t)*hb**2)+0.5d0*((1d0-c2t)*ht**2)
      tmp63 = 0.5d0*((1d0-c2t)*hb**2)+0.5d0*((1d0+c2t)*ht**2)
      tmp64 = (A0-b)/b+LogA0-Logb-t/b
      tmp65 = (-LogA0+Logt)*(A0-t)+(-LogA0+Logt)*t
      tmp66 = (A0-b)*(-LogA0+Logb)+(-LogA0-Logb+2d0*Logt)*t
      tmp67 = (-LogB1+Logmu2)*(B1-mu2)+
     .   (-LogB1-Logmu2+2d0*Logt)*t
      tmp68 = (-LogB2+Logmu2)*(B2-mu2)+
     .   (-LogB2-Logmu2+2d0*Logt)*t
      tmp69 = b*(-LogA0+2d0*Logb-Logt)+(LogA0-Logt)*(-A0+t)
      tmp70 = (-LogA0+Logt)*t+(LogA0-Logt)*(-A0+t)
      tmp71 = B1*(2d0*LogB1-Logmu2-Logt)+
     .   (Logmu2-Logt)*(-mu2+t)
      tmp72 = B2*(2d0*LogB2-Logmu2-Logt)+
     .   (Logmu2-Logt)*(-mu2+t)
      tmp73 = Logmu2-Logt-B1/t-(-mu2+t)/t
      tmp74 = Logmu2-Logt-B2/t-(-mu2+t)/t
      tmp75 = 2d0*B1*LogB1-2d0*B2*LogB2-
     .   0.5d0*(deltB1tmu2*phiB1tmu2)/mu2+
     .   0.5d0*(deltB2tmu2*phiB2tmu2)/mu2+
     .   0.5d0*(Logmu2*Logt*(B1-mu2-t))-
     .   0.5d0*(Logmu2*Logt*(B2-mu2-t))+
     .   0.5d0*(LogB1*Logt*(-B1+mu2-t))-
     .   0.5d0*(LogB2*Logt*(-B2+mu2-t))+
     .   0.5d0*(LogB1*Logmu2*(-B1-mu2+t))-
     .   0.5d0*(LogB2*Logmu2*(-B2-mu2+t))-2.5d0*(B1+mu2+t)+
     .   2.5d0*(B2+mu2+t)
      tmp76 = -5d0*b+4d0*b*Logb+Logt**2*(b-t)-5d0*t-
     .   2d0*Li2bt*(-b+t)+Logt*(-2d0*b*Logb+4d0*t)
      tmp77 = -6d0+4d0*Logt+(2d0-2d0*Logt)*Logt+(4d0*t-2d0*Logt*t)/t
      tmp78 = -Logb+Logmu2-(b-mu2)/b-T1/b
      tmp79 = 2d0*(-1d0+LogT1)*T1+2d0*(-1d0+LogT1)**2*T1
      tmp80 = (-LogA0+LogT1)*(A0-T1)+(-LogA0+LogT1)*T1
      tmp81 = (A0-B1)*(-LogA0+LogB1)+
     .   (-LogA0-LogB1+2d0*LogT1)*T1
      tmp82 = (A0-B2)*(-LogA0+LogB2)+
     .   (-LogA0-LogB2+2d0*LogT1)*T1
      tmp83 = (-Logb+Logmu2)*(b-mu2)+
     .   (-Logb-Logmu2+2d0*LogT1)*T1
      tmp84 = (Logmu2-Logt)*(-mu2+t)+
     .   (-Logmu2-Logt+2d0*LogT1)*T1
      tmp85 = B1*(-LogA0+2d0*LogB1-LogT1)+
     .   (LogA0-LogT1)*(-A0+T1)
      tmp86 = B2*(-LogA0+2d0*LogB2-LogT1)+
     .   (LogA0-LogT1)*(-A0+T1)
      tmp87 = (-LogA0+LogT1)*T1+(LogA0-LogT1)*(-A0+T1)
      tmp88 = 2d0*A0*LogA0+2d0*B1*LogB1+2d0*LogT1*T1+
     .   0.5d0*(LogB1*LogT1*(A0-B1-T1))+
     .   0.5d0*(LogA0*LogT1*(-A0+B1-T1))-
     .   0.5d0*(deltA0B1T1*phiA0B1T1)/T1+
     .   0.5d0*(LogA0*LogB1*(-A0-B1+T1))-2.5d0*(A0+B1+T1)
      tmp89 = 2d0*A0*LogA0+2d0*B2*LogB2+2d0*LogT1*T1+
     .   0.5d0*(LogB2*LogT1*(A0-B2-T1))+
     .   0.5d0*(LogA0*LogT1*(-A0+B2-T1))-
     .   0.5d0*(deltA0B2T1*phiA0B2T1)/T1+
     .   0.5d0*(LogA0*LogB2*(-A0-B2+T1))-2.5d0*(A0+B2+T1)
      tmp90 = b*(2d0*Logb-Logmu2-LogT1)+
     .   (Logmu2-LogT1)*(-mu2+T1)
      tmp91 = (-Logmu2+2d0*Logt-LogT1)*t+
     .   (Logmu2-LogT1)*(-mu2+T1)
      tmp92 = 2d0*b*Logb+2d0*Logmu2*mu2+2d0*LogT1*T1-
     .   0.5d0*(deltT1bmu2*phiT1bmu2)/mu2+
     .   0.5d0*(Logmu2*LogT1*(b-mu2-T1))+
     .   0.5d0*(Logb*LogT1*(-b+mu2-T1))+
     .   0.5d0*(Logb*Logmu2*(-b-mu2+T1))-2.5d0*(b+mu2+T1)
      tmp93 = 2d0*A0*LogA0-A0*LogA0*LogT1+4d0*LogT1*T1+
     .   0.5d0*(LogT1**2*(A0-2d0*T1))-
     .   0.5d0*(deltA0T1T1*phiA0T1T1)/T1-2.5d0*(A0+2d0*T1)
      tmp94 = -10d0*T1+4d0*LogT1*T1+LogT1*(4d0*T1-2d0*LogT1*T1)
      tmp95 = -6d0+4d0*LogT1+(2d0-2d0*LogT1)*LogT1+
     .   (4d0*T1-2d0*LogT1*T1)/T1
      tmp96 = (-LogA0+2d0*LogT1-LogT2)*T1+
     .   (-LogA0+LogT2)*(A0-T2)
      tmp97 = -Logb+Logmu2-(b-mu2)/b-T2/b
      tmp98 = 2d0*(-1d0+LogT2)*T2+2d0*(-1d0+LogT2)**2*T2
      tmp99 = (-LogA0+LogT2)*(A0-T2)+(-LogA0+LogT2)*T2
      tmp100 = (A0-B1)*(-LogA0+LogB1)+
     .   (-LogA0-LogB1+2d0*LogT2)*T2
      tmp101 = (A0-B2)*(-LogA0+LogB2)+
     .   (-LogA0-LogB2+2d0*LogT2)*T2
      tmp102 = (-Logb+Logmu2)*(b-mu2)+
     .   (-Logb-Logmu2+2d0*LogT2)*T2
      tmp103 = (Logmu2-Logt)*(-mu2+t)+
     .   (-Logmu2-Logt+2d0*LogT2)*T2
      tmp104 = (LogA0-LogT1)*(-A0+T1)+
     .   (-LogA0-LogT1+2d0*LogT2)*T2
      tmp105 = B1*(-LogA0+2d0*LogB1-LogT2)+
     .   (LogA0-LogT2)*(-A0+T2)
      tmp106 = B2*(-LogA0+2d0*LogB2-LogT2)+
     .   (LogA0-LogT2)*(-A0+T2)
      tmp107 = (-LogA0+LogT2)*T2+(LogA0-LogT2)*(-A0+T2)
      tmp108 = 2d0*A0*LogA0+2d0*B1*LogB1+2d0*LogT2*T2+
     .   0.5d0*(LogB1*LogT2*(A0-B1-T2))+
     .   0.5d0*(LogA0*LogT2*(-A0+B1-T2))-
     .   0.5d0*(deltA0B1T2*phiA0B1T2)/T2+
     .   0.5d0*(LogA0*LogB1*(-A0-B1+T2))-2.5d0*(A0+B1+T2)
      tmp109 = 2d0*A0*LogA0+2d0*B2*LogB2+2d0*LogT2*T2+
     .   0.5d0*(LogB2*LogT2*(A0-B2-T2))+
     .   0.5d0*(LogA0*LogT2*(-A0+B2-T2))-
     .   0.5d0*(deltA0B2T2*phiA0B2T2)/T2+
     .   0.5d0*(LogA0*LogB2*(-A0-B2+T2))-2.5d0*(A0+B2+T2)
      tmp110 = b*(2d0*Logb-Logmu2-LogT2)+
     .   (Logmu2-LogT2)*(-mu2+T2)
      tmp111 = (-Logmu2+2d0*Logt-LogT2)*t+
     .   (Logmu2-LogT2)*(-mu2+T2)
      tmp112 = 2d0*b*Logb+2d0*Logmu2*mu2+2d0*LogT2*T2-
     .   0.5d0*(delT2Tbmu2*phiT2bmu2)/mu2+
     .   0.5d0*(Logmu2*LogT2*(b-mu2-T2))+
     .   0.5d0*(Logb*LogT2*(-b+mu2-T2))+
     .   0.5d0*(Logb*Logmu2*(-b-mu2+T2))-2.5d0*(b+mu2+T2)
      tmp113 = 2d0*LogT1*T1-2d0*LogT2*T2-
     .   0.5d0*(deltT1bmu2*phiT1bmu2)/mu2+
     .   0.5d0*(delT2Tbmu2*phiT2bmu2)/mu2+
     .   0.5d0*(Logmu2*LogT1*(b-mu2-T1))+
     .   0.5d0*(Logb*LogT1*(-b+mu2-T1))+
     .   0.5d0*(Logb*Logmu2*(-b-mu2+T1))-2.5d0*(b+mu2+T1)-
     .   0.5d0*(Logmu2*LogT2*(b-mu2-T2))-
     .   0.5d0*(Logb*LogT2*(-b+mu2-T2))-
     .   0.5d0*(Logb*Logmu2*(-b-mu2+T2))+2.5d0*(b+mu2+T2)
      tmp114 = 2d0*A0*LogA0-A0*LogA0*LogT2+4d0*LogT2*T2+
     .   0.5d0*(LogT2**2*(A0-2d0*T2))-
     .   0.5d0*(deltA0T2T2*phiA0T2T2)/T2-2.5d0*(A0+2d0*T2)
      tmp115 = -10d0*T2+4d0*LogT2*T2+LogT2*(4d0*T2-2d0*LogT2*T2)
      tmp116 = -6d0+4d0*LogT2+(2d0-2d0*LogT2)*LogT2+
     .   (4d0*T2-2d0*LogT2*T2)/T2
      tmp117 = 0.25d0*((1d0-c2b)*(1d0+c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0-c2b)*s2t*Xb)/mt
      tmp118 = 0.25d0*((1d0-c2b)*(1d0-c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0-c2b)*s2t*Xb)/mt
      tmp119 = 0.25d0*((1d0+c2b)*(1d0-c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0+c2b)*s2t*Xb)/mt
      tmp120 = 0.25d0*(b*(1d0+c2b)*(1d0-c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0+c2t)*t)+
     .   0.5d0*((1d0-c2t)*mb*s2b*Xb)-0.5d0*((1d0-c2b)*mt*s2t*Xb)+
     .   0.25d0*((1d0-c2b)*(1d0-c2t)*Xb**2)
      tmp121 = 0.25d0*(b*(1d0-c2b)*(1d0-c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0+c2t)*t)-
     .   0.5d0*((1d0-c2t)*mb*s2b*Xb)-0.5d0*((1d0+c2b)*mt*s2t*Xb)+
     .   0.25d0*((1d0+c2b)*(1d0-c2t)*Xb**2)
      tmp122 = -(0.125d0*(b*(1d0+c2b))/c2t)+
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0-c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Xb)/c2t+0.25d0*((1d0-c2b)*mt*Xb)/s2t-
     .   0.125d0*((1d0-c2b)*Xb**2)/c2t
      tmp123 = 0.125d0*(b*(1d0+c2b))/c2t-
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0-c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Xb)/c2t-0.25d0*((1d0-c2b)*mt*Xb)/s2t+
     .   0.125d0*((1d0-c2b)*Xb**2)/c2t
      tmp124 = -(0.125d0*(b*(1d0-c2b))/c2t)-
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0+c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Xb)/c2t+0.25d0*((1d0+c2b)*mt*Xb)/s2t-
     .   0.125d0*((1d0+c2b)*Xb**2)/c2t
      tmp125 = 0.125d0*(b*(1d0-c2b))/c2t+
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0+c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Xb)/c2t-0.25d0*((1d0+c2b)*mt*Xb)/s2t+
     .   0.125d0*((1d0+c2b)*Xb**2)/c2t
      tmp126 = 0.25d0*(b*(1d0+c2b)*(1d0+c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0-c2t)*t)+
     .   0.5d0*((1d0+c2t)*mb*s2b*Xb)+0.5d0*((1d0-c2b)*mt*s2t*Xb)+
     .   0.25d0*((1d0-c2b)*(1d0+c2t)*Xb**2)
      tmp127 = 0.25d0*(b*(1d0-c2b)*(1d0+c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0-c2t)*t)-
     .   0.5d0*((1d0+c2t)*mb*s2b*Xb)+0.5d0*((1d0+c2b)*mt*s2t*Xb)+
     .   0.25d0*((1d0+c2b)*(1d0+c2t)*Xb**2)
      tmp128 = 0.5d0*((1d0+c2b)*Xb)+0.5d0*((1d0-c2b)*Xt)
      tmp129 = 0.5d0*((1d0-c2b)*Xb)+0.5d0*((1d0+c2b)*Xt)
      tmp130 = 0.5d0*((1d0+c2t)*Xb)+0.5d0*((1d0-c2t)*Xt)
      tmp131 = 0.5d0*((1d0-c2t)*Xb)+0.5d0*((1d0+c2t)*Xt)
      tmp132 = 0.25d0*((1d0-c2b)*(1d0-c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0-c2b)*s2t*Xt)/mt
      tmp133 = 0.25d0*((1d0-c2b)*(1d0+c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0-c2b)*s2t*Xt)/mt
      tmp134 = 0.25d0*((1d0+c2b)*(1d0-c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0+c2b)*s2t*Xt)/mt
      tmp135 = 0.25d0*(b*(1d0+c2b)*(1d0-c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0+c2t)*t)-
     .   0.5d0*((1d0-c2t)*mb*s2b*Xt)+0.5d0*((1d0-c2b)*mt*s2t*Xt)+
     .   0.25d0*((1d0-c2b)*(1d0-c2t)*Xt**2)
      tmp136 = 0.25d0*(b*(1d0-c2b)*(1d0-c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0+c2t)*t)+
     .   0.5d0*((1d0-c2t)*mb*s2b*Xt)+0.5d0*((1d0+c2b)*mt*s2t*Xt)+
     .   0.25d0*((1d0+c2b)*(1d0-c2t)*Xt**2)
      tmp137 = -(0.125d0*(b*(1d0+c2b))/c2t)+
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0-c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Xt)/c2t-0.25d0*((1d0-c2b)*mt*Xt)/s2t-
     .   0.125d0*((1d0-c2b)*Xt**2)/c2t
      tmp138 = 0.125d0*(b*(1d0+c2b))/c2t-
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0-c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Xt)/c2t+0.25d0*((1d0-c2b)*mt*Xt)/s2t+
     .   0.125d0*((1d0-c2b)*Xt**2)/c2t
      tmp139 = -(0.125d0*(b*(1d0-c2b))/c2t)-
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0+c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Xt)/c2t-0.25d0*((1d0+c2b)*mt*Xt)/s2t-
     .   0.125d0*((1d0+c2b)*Xt**2)/c2t
      tmp140 = 0.125d0*(b*(1d0-c2b))/c2t+
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0+c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Xt)/c2t+0.25d0*((1d0+c2b)*mt*Xt)/s2t+
     .   0.125d0*((1d0+c2b)*Xt**2)/c2t
      tmp141 = 0.25d0*(b*(1d0+c2b)*(1d0+c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0-c2t)*t)-
     .   0.5d0*((1d0+c2t)*mb*s2b*Xt)-0.5d0*((1d0-c2b)*mt*s2t*Xt)+
     .   0.25d0*((1d0-c2b)*(1d0+c2t)*Xt**2)
      tmp142 = 0.25d0*(b*(1d0-c2b)*(1d0+c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0-c2t)*t)+
     .   0.5d0*((1d0+c2t)*mb*s2b*Xt)-0.5d0*((1d0+c2b)*mt*s2t*Xt)+
     .   0.25d0*((1d0+c2b)*(1d0+c2t)*Xt**2)
      tmp143 = 0.25d0*((1d0-c2b)*(1d0+c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0-c2b)*s2t*Yb)/mt
      tmp144 = 0.25d0*((1d0-c2b)*(1d0-c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0-c2b)*s2t*Yb)/mt
      tmp145 = 0.25d0*((1d0+c2b)*(1d0-c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0+c2b)*s2t*Yb)/mt
      tmp146 = 0.25d0*(b*(1d0+c2b)*(1d0-c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0+c2t)*t)+
     .   0.5d0*((1d0-c2t)*mb*s2b*Yb)-0.5d0*((1d0-c2b)*mt*s2t*Yb)+
     .   0.25d0*((1d0-c2b)*(1d0-c2t)*Yb**2)
      tmp147 = 0.25d0*(b*(1d0-c2b)*(1d0-c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0+c2t)*t)-
     .   0.5d0*((1d0-c2t)*mb*s2b*Yb)-0.5d0*((1d0+c2b)*mt*s2t*Yb)+
     .   0.25d0*((1d0+c2b)*(1d0-c2t)*Yb**2)
      tmp148 = -(0.125d0*(b*(1d0+c2b))/c2t)+
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0-c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Yb)/c2t+0.25d0*((1d0-c2b)*mt*Yb)/s2t-
     .   0.125d0*((1d0-c2b)*Yb**2)/c2t
      tmp149 = 0.125d0*(b*(1d0+c2b))/c2t-
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0-c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Yb)/c2t-0.25d0*((1d0-c2b)*mt*Yb)/s2t+
     .   0.125d0*((1d0-c2b)*Yb**2)/c2t
      tmp150 = -(0.125d0*(b*(1d0-c2b))/c2t)-
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0+c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Yb)/c2t+0.25d0*((1d0+c2b)*mt*Yb)/s2t-
     .   0.125d0*((1d0+c2b)*Yb**2)/c2t
      tmp151 = 0.125d0*(b*(1d0-c2b))/c2t+
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0+c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Yb)/c2t-0.25d0*((1d0+c2b)*mt*Yb)/s2t+
     .   0.125d0*((1d0+c2b)*Yb**2)/c2t
      tmp152 = 0.25d0*(b*(1d0+c2b)*(1d0+c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0-c2t)*t)+
     .   0.5d0*((1d0+c2t)*mb*s2b*Yb)+0.5d0*((1d0-c2b)*mt*s2t*Yb)+
     .   0.25d0*((1d0-c2b)*(1d0+c2t)*Yb**2)
      tmp153 = 0.25d0*(b*(1d0-c2b)*(1d0+c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0-c2t)*t)-
     .   0.5d0*((1d0+c2t)*mb*s2b*Yb)+0.5d0*((1d0+c2b)*mt*s2t*Yb)+
     .   0.25d0*((1d0+c2b)*(1d0+c2t)*Yb**2)
      tmp154 = 0.5d0*((1d0+c2b)*Yb)+0.5d0*((1d0-c2b)*Yt)
      tmp155 = 0.5d0*((1d0-c2b)*Yb)+0.5d0*((1d0+c2b)*Yt)
      tmp156 = 0.5d0*((1d0+c2t)*Yb)+0.5d0*((1d0-c2t)*Yt)
      tmp157 = 0.5d0*((1d0-c2t)*Yb)+0.5d0*((1d0+c2t)*Yt)
      tmp158 = 0.25d0*((1d0-c2b)*(1d0-c2t))+
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0-c2b)*s2t*Yt)/mt
      tmp159 = 0.25d0*((1d0-c2b)*(1d0+c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt+0.25d0*((1d0-c2b)*s2t*Yt)/mt
      tmp160 = 0.25d0*((1d0+c2b)*(1d0-c2t))-
     .   0.25d0*(mb*s2b*s2t)/mt-0.25d0*((1d0+c2b)*s2t*Yt)/mt
      tmp161 = 0.25d0*(b*(1d0+c2b)*(1d0-c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0+c2t)*t)-
     .   0.5d0*((1d0-c2t)*mb*s2b*Yt)+0.5d0*((1d0-c2b)*mt*s2t*Yt)+
     .   0.25d0*((1d0-c2b)*(1d0-c2t)*Yt**2)
      tmp162 = 0.25d0*(b*(1d0-c2b)*(1d0-c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0+c2t)*t)+
     .   0.5d0*((1d0-c2t)*mb*s2b*Yt)+0.5d0*((1d0+c2b)*mt*s2t*Yt)+
     .   0.25d0*((1d0+c2b)*(1d0-c2t)*Yt**2)
      tmp163 = -(0.125d0*(b*(1d0+c2b))/c2t)+
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0-c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Yt)/c2t-0.25d0*((1d0-c2b)*mt*Yt)/s2t-
     .   0.125d0*((1d0-c2b)*Yt**2)/c2t
      tmp164 = 0.125d0*(b*(1d0+c2b))/c2t-
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0-c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Yt)/c2t+0.25d0*((1d0-c2b)*mt*Yt)/s2t+
     .   0.125d0*((1d0-c2b)*Yt**2)/c2t
      tmp165 = -(0.125d0*(b*(1d0-c2b))/c2t)-
     .   0.25d0*(mb*mt*s2b)/s2t+0.125d0*((1d0+c2b)*t)/c2t-
     .   0.25d0*(mb*s2b*Yt)/c2t-0.25d0*((1d0+c2b)*mt*Yt)/s2t-
     .   0.125d0*((1d0+c2b)*Yt**2)/c2t
      tmp166 = 0.125d0*(b*(1d0-c2b))/c2t+
     .   0.25d0*(mb*mt*s2b)/s2t-0.125d0*((1d0+c2b)*t)/c2t+
     .   0.25d0*(mb*s2b*Yt)/c2t+0.25d0*((1d0+c2b)*mt*Yt)/s2t+
     .   0.125d0*((1d0+c2b)*Yt**2)/c2t
      tmp167 = 0.25d0*(b*(1d0+c2b)*(1d0+c2t))+
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0-c2b)*(1d0-c2t)*t)-
     .   0.5d0*((1d0+c2t)*mb*s2b*Yt)-0.5d0*((1d0-c2b)*mt*s2t*Yt)+
     .   0.25d0*((1d0-c2b)*(1d0+c2t)*Yt**2)
      tmp168 = 0.25d0*(b*(1d0-c2b)*(1d0+c2t))-
     .   0.5d0*(mb*mt*s2b*s2t)+0.25d0*((1d0+c2b)*(1d0-c2t)*t)+
     .   0.5d0*((1d0+c2t)*mb*s2b*Yt)-0.5d0*((1d0+c2b)*mt*s2t*Yt)+
     .   0.25d0*((1d0+c2b)*(1d0+c2t)*Yt**2)
      tmp169 = -Li2T1B1-0.5d0*(-LogB1+LogT1)**2
      tmp170 = -Li2T1B2-0.5d0*(-LogB2+LogT1)**2
      tmp171 = -Li2T1T2-0.5d0*(LogT1-LogT2)**2
      tmp172 = -Li2T2B1-0.5d0*(-LogB1+LogT2)**2
      tmp173 = -Li2T2B2-0.5d0*(-LogB2+LogT2)**2
      tmp174 = -A0+(A0-t)**2/b-t
      tmp175 = -A0+(A0-t)**2/t-t
      tmp176 = -A0+(A0-T1)**2/B1-T1
      tmp177 = -A0+(A0-T1)**2/B2-T1
      tmp178 = -A0+(A0-T1)**2/T1-T1
      tmp179 = -A0-T1+(A0-T1)**2/T2
      tmp180 = -A0+(A0-T2)**2/B1-T2
      tmp181 = -A0+(A0-T2)**2/B2-T2
      tmp182 = -A0+(A0-T2)**2/T2-T2
      tmp183 = cbe**2*hb**2*tmp119+ht**2*sbe**2*tmp133-
     .   cbe*hb*ht*sbe*((mb*tmp56)/mt-0.5d0*(s2b*s2t)-
     .      0.5d0*(s2b*tmp130)/mt)
      tmp184 = cbe**2*hb**2*tmp118+ht**2*sbe**2*tmp29-
     .   cbe*hb*ht*sbe*((mb*tmp59)/mt+0.5d0*(s2b*s2t)+
     .      0.5d0*(s2b*tmp130)/mt)
      tmp185 = cbe**2*hb**2*tmp127+ht**2*sbe**2*tmp135-
     .   cbe*hb*ht*sbe*(mb*s2t*tmp128-mt*s2b*tmp130+
     .      2d0*mb*mt*tmp56-0.5d0*(s2b*s2t*(b+t))-
     .      0.5d0*(s2b*s2t*Xb*Xt))
      tmp186 = cbe**2*hb**2*tmp126+ht**2*sbe**2*tmp136-
     .   cbe*hb*ht*sbe*(mb*s2t*tmp129+mt*s2b*tmp130+
     .      2d0*mb*mt*tmp59+0.5d0*(s2b*s2t*(b+t))+
     .      0.5d0*(s2b*s2t*Xb*Xt))
      tmp187 = ht**2*sbe**2*tmp132+cbe**2*hb**2*tmp20-
     .   cbe*hb*ht*sbe*((mb*tmp59)/mt+0.5d0*(s2b*s2t)-
     .      0.5d0*(s2b*tmp131)/mt)
      tmp188 = cbe**2*hb**2*tmp117+ht**2*sbe**2*tmp134-
     .   cbe*hb*ht*sbe*((mb*tmp56)/mt-0.5d0*(s2b*s2t)+
     .      0.5d0*(s2b*tmp131)/mt)
      tmp189 = cbe**2*hb**2*tmp121+ht**2*sbe**2*tmp141-
     .   cbe*hb*ht*sbe*(-(mb*s2t*tmp128)-mt*s2b*tmp131+
     .      2d0*mb*mt*tmp59+0.5d0*(s2b*s2t*(b+t))+
     .      0.5d0*(s2b*s2t*Xb*Xt))
      tmp190 = cbe**2*hb**2*tmp120+ht**2*sbe**2*tmp142-
     .   cbe*hb*ht*sbe*(-(mb*s2t*tmp129)+mt*s2b*tmp131+
     .      2d0*mb*mt*tmp56-0.5d0*(s2b*s2t*(b+t))-
     .      0.5d0*(s2b*s2t*Xb*Xt))
      tmp191 = hb**2*sbe**2*tmp145+cbe**2*ht**2*tmp159+
     .   cbe*hb*ht*sbe*((mb*tmp56)/mt-0.5d0*(s2b*s2t)-
     .      0.5d0*(s2b*tmp156)/mt)
      tmp192 = hb**2*sbe**2*tmp144+cbe**2*ht**2*tmp45+
     .   cbe*hb*ht*sbe*((mb*tmp59)/mt+0.5d0*(s2b*s2t)+
     .      0.5d0*(s2b*tmp156)/mt)
      tmp193 = hb**2*sbe**2*tmp153+cbe**2*ht**2*tmp161+
     .   cbe*hb*ht*sbe*(mb*s2t*tmp154-mt*s2b*tmp156+
     .      2d0*mb*mt*tmp56-0.5d0*(s2b*s2t*(b+t))-
     .      0.5d0*(s2b*s2t*Yb*Yt))
      tmp194 = hb**2*sbe**2*tmp152+cbe**2*ht**2*tmp162+
     .   cbe*hb*ht*sbe*(mb*s2t*tmp155+mt*s2b*tmp156+
     .      2d0*mb*mt*tmp59+0.5d0*(s2b*s2t*(b+t))+
     .      0.5d0*(s2b*s2t*Yb*Yt))
      tmp195 = cbe**2*ht**2*tmp158+hb**2*sbe**2*tmp36+
     .   cbe*hb*ht*sbe*((mb*tmp59)/mt+0.5d0*(s2b*s2t)-
     .      0.5d0*(s2b*tmp157)/mt)
      tmp196 = hb**2*sbe**2*tmp143+cbe**2*ht**2*tmp160+
     .   cbe*hb*ht*sbe*((mb*tmp56)/mt-0.5d0*(s2b*s2t)+
     .      0.5d0*(s2b*tmp157)/mt)
      tmp197 = hb**2*sbe**2*tmp147+cbe**2*ht**2*tmp167+
     .   cbe*hb*ht*sbe*(-(mb*s2t*tmp154)-mt*s2b*tmp157+
     .      2d0*mb*mt*tmp59+0.5d0*(s2b*s2t*(b+t))+
     .      0.5d0*(s2b*s2t*Yb*Yt))
      tmp198 = hb**2*sbe**2*tmp146+cbe**2*ht**2*tmp168+
     .   cbe*hb*ht*sbe*(-(mb*s2t*tmp155)+mt*s2b*tmp157+
     .      2d0*mb*mt*tmp56-0.5d0*(s2b*s2t*(b+t))-
     .      0.5d0*(s2b*s2t*Yb*Yt))
      tmp199 = cbe**2*hb**2*tmp125+ht**2*sbe**2*tmp137-
     .   cbe*hb*ht*sbe*(-(mt*s2b*tmp27)+2d0*mb*mt*tmp54+
     .      0.25d0*(s2b*(b+t))/s2t-0.5d0*(mb*tmp128)/s2t+
     .      0.25d0*(s2b*Xb*Xt)/s2t)
      tmp200 = cbe**2*hb**2*tmp123+ht**2*sbe**2*tmp139-
     .   cbe*hb*ht*sbe*(mt*s2b*tmp27+2d0*mb*mt*tmp55-
     .      0.25d0*(s2b*(b+t))/s2t-0.5d0*(mb*tmp129)/s2t-
     .      0.25d0*(s2b*Xb*Xt)/s2t)
      tmp201 = cbe**2*hb**2*tmp124+ht**2*sbe**2*tmp138-
     .   cbe*hb*ht*sbe*(-(mt*s2b*tmp28)+2d0*mb*mt*tmp55-
     .      0.25d0*(s2b*(b+t))/s2t+0.5d0*(mb*tmp128)/s2t-
     .      0.25d0*(s2b*Xb*Xt)/s2t)
      tmp202 = cbe**2*hb**2*tmp122+ht**2*sbe**2*tmp140-
     .   cbe*hb*ht*sbe*(mt*s2b*tmp28+2d0*mb*mt*tmp54+
     .      0.25d0*(s2b*(b+t))/s2t+0.5d0*(mb*tmp129)/s2t+
     .      0.25d0*(s2b*Xb*Xt)/s2t)
      tmp203 = hb**2*sbe**2*tmp151+cbe**2*ht**2*tmp163+
     .   cbe*hb*ht*sbe*(-(mt*s2b*tmp43)+2d0*mb*mt*tmp54+
     .      0.25d0*(s2b*(b+t))/s2t-0.5d0*(mb*tmp154)/s2t+
     .      0.25d0*(s2b*Yb*Yt)/s2t)
      tmp204 = hb**2*sbe**2*tmp149+cbe**2*ht**2*tmp165+
     .   cbe*hb*ht*sbe*(mt*s2b*tmp43+2d0*mb*mt*tmp55-
     .      0.25d0*(s2b*(b+t))/s2t-0.5d0*(mb*tmp155)/s2t-
     .      0.25d0*(s2b*Yb*Yt)/s2t)
      tmp205 = hb**2*sbe**2*tmp150+cbe**2*ht**2*tmp164+
     .   cbe*hb*ht*sbe*(-(mt*s2b*tmp44)+2d0*mb*mt*tmp55-
     .      0.25d0*(s2b*(b+t))/s2t+0.5d0*(mb*tmp154)/s2t-
     .      0.25d0*(s2b*Yb*Yt)/s2t)
      tmp206 = hb**2*sbe**2*tmp148+cbe**2*ht**2*tmp166+
     .   cbe*hb*ht*sbe*(mt*s2b*tmp44+2d0*mb*mt*tmp54+
     .      0.25d0*(s2b*(b+t))/s2t+0.5d0*(mb*tmp155)/s2t+
     .      0.25d0*(s2b*Yb*Yt)/s2t)
      tmp207 = (b**2*(Logb-Logt))/((1d0-b/t)**2*t**4)+
     .   b/((1d0-b/t)*t**3)+(2d0*b*(Logb-Logt))/((1d0-b/t)*t**3)
      tmp208 = (-2d0*b*(Logb-Logt))/((1d0-b/t)*t**2)+
     .   (-2d0-2d0*Logb)/t+(2d0*Logt)/t-
     .   (2d0*(Logb-Logt))/((1d0-b/t)*t)+
     .   (2d0*b*(Logb-Logt)*(-b+t))/((1d0-b/t)**2*t**3)+
     .   (2d0*(-b+t))/((1d0-b/t)*t**2)+
     .   (2d0*(Logb-Logt)*(-b+t))/((1d0-b/t)*t**2)
      tmp209 = -1d0+2d0*Li2bt+4d0*Logb+(-2d0-2d0*Logb)*Logt+
     .   Logt**2-(2d0*(Logb-Logt)*(-b+t))/((1d0-b/t)*t)
      tmp210 = -5d0-2d0*Li2bt+4d0*Logt-Logt**2+
     .   (2d0*Logt*(b-t))/t+
     .   (2d0*b*(Logb-Logt)*(-b+t))/((1d0-b/t)*t**2)+
     .   (-2d0*b*Logb+4d0*t)/t
      tmp211 = -1d0+4d0*LogB1+(-2d0-2d0*LogB1)*LogT1+
     .   LogT1**2-(2d0*(LogB1-LogT1)*(-B1+T1))/
     .    ((1d0-B1/T1)*T1)+2d0*tmp169
      tmp212 = -1d0+4d0*LogB2+(-2d0-2d0*LogB2)*LogT1+
     .   LogT1**2-(2d0*(LogB2-LogT1)*(-B2+T1))/
     .    ((1d0-B2/T1)*T1)+2d0*tmp170
      tmp213 = -5d0+4d0*LogT1-LogT1**2+
     .   (2d0*LogT1*(B1-T1))/T1+
     .   (2d0*B1*(LogB1-LogT1)*(-B1+T1))/((1d0-B1/T1)*T1**2)+
     .   (-2d0*B1*LogB1+4d0*T1)/T1-2d0*tmp169
      tmp214 = -5d0+4d0*LogT1-LogT1**2+
     .   (2d0*LogT1*(B2-T1))/T1+
     .   (2d0*B2*(LogB2-LogT1)*(-B2+T1))/((1d0-B2/T1)*T1**2)+
     .   (-2d0*B2*LogB2+4d0*T1)/T1-2d0*tmp170
      tmp215 = -1d0+4d0*LogB1+(-2d0-2d0*LogB1)*LogT2+
     .   LogT2**2-(2d0*(LogB1-LogT2)*(-B1+T2))/
     .    ((1d0-B1/T2)*T2)+2d0*tmp172
      tmp216 = -1d0+4d0*LogB2+(-2d0-2d0*LogB2)*LogT2+
     .   LogT2**2-(2d0*(LogB2-LogT2)*(-B2+T2))/
     .    ((1d0-B2/T2)*T2)+2d0*tmp173
      tmp217 = -5d0+4d0*LogT2-LogT2**2+
     .   (2d0*LogT2*(B1-T2))/T2+
     .   (2d0*B1*(LogB1-LogT2)*(-B1+T2))/((1d0-B1/T2)*T2**2)+
     .   (-2d0*B1*LogB1+4d0*T2)/T2-2d0*tmp172
      tmp218 = -5d0+4d0*LogT2-LogT2**2+
     .   (2d0*LogT2*(B2-T2))/T2+
     .   (2d0*B2*(LogB2-LogT2)*(-B2+T2))/((1d0-B2/T2)*T2**2)+
     .   (-2d0*B2*LogB2+4d0*T2)/T2-2d0*tmp173
      tmp219 = -1d0+LogT1**2+LogT1*(-2d0-2d0*LogT2)+
     .   4d0*LogT2-(2d0*(-LogT1+LogT2)*(T1-T2))/
     .    (T1*(1d0-T2/T1))+2d0*tmp171
      tmp220 = -5d0+4d0*LogT1-LogT1**2+
     .   (2d0*LogT1*(-T1+T2))/T1+(4d0*T1-2d0*LogT2*T2)/T1+
     .   (2d0*(-LogT1+LogT2)*(T1-T2)*T2)/(T1**2*(1d0-T2/T1))-
     .   2d0*tmp171
      tmp221 = -0.5d0+2d0*LogT1-
     .   (phiA0T1T2*(-A0+T1-T2))/T2-0.5d0*(LogA0*LogT1)+
     .   0.5d0*(LogA0*LogT2)-0.5d0*(LogT1*LogT2)+
     .   0.5d0*(LogT2*(A0-T1-T2))/T1+
     .   0.5d0*(LogA0*(-A0-T1+T2))/T1-
     .   0.5d0*(deltA0T1T2*((phiA0T1T2*(A0-T1+T2))/deltA0T1T2+
     .       (T2*tmp96)/(deltA0T1T2*T1)))/T2
      tmp222 = -0.5d0+2d0*LogT2-
     .   (phiA0B1T2*(-A0-B1+T2))/T2+0.5d0*(LogA0*LogB1)-
     .   0.5d0*(LogA0*LogT2)-0.5d0*(LogB1*LogT2)+
     .   0.5d0*(LogB1*(A0-B1-T2))/T2+
     .   0.5d0*(LogA0*(-A0+B1-T2))/T2-
     .   0.5d0*(deltA0B1T2*((B1*phiA0B1T2*(A0+B1-T2))/
     .        (deltA0B1T2*T2)+(B1*tmp100)/(deltA0B1T2*T2)))/B1
      tmp223 = -0.5d0+2d0*LogT2-
     .   (phiA0B2T2*(-A0-B2+T2))/T2+0.5d0*(LogA0*LogB2)-
     .   0.5d0*(LogA0*LogT2)-0.5d0*(LogB2*LogT2)+
     .   0.5d0*(LogB2*(A0-B2-T2))/T2+
     .   0.5d0*(LogA0*(-A0+B2-T2))/T2-
     .   0.5d0*(deltA0B2T2*((B2*phiA0B2T2*(A0+B2-T2))/
     .        (deltA0B2T2*T2)+(B2*tmp101)/(deltA0B2T2*T2)))/B2
      tmp224 = -0.5d0+2d0*LogT2-
     .   (phiT2bmu2*(-b-mu2+T2))/mu2+0.5d0*(Logb*Logmu2)-
     .   0.5d0*(Logb*LogT2)-0.5d0*(Logmu2*LogT2)+
     .   0.5d0*(Logmu2*(b-mu2-T2))/T2+
     .   0.5d0*(Logb*(-b+mu2-T2))/T2-
     .   0.5d0*(delT2Tbmu2*((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .       (mu2*tmp102)/(delT2Tbmu2*T2)))/mu2
      tmp225 = 0.5d0-2d0*LogT2+
     .   (phiT2bmu2*(-b-mu2+T2))/mu2-0.5d0*(Logb*Logmu2)+
     .   0.5d0*(Logb*LogT2)+0.5d0*(Logmu2*LogT2)-
     .   0.5d0*(Logmu2*(b-mu2-T2))/T2-
     .   0.5d0*(Logb*(-b+mu2-T2))/T2+
     .   0.5d0*(delT2Tbmu2*((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .       (mu2*tmp102)/(delT2Tbmu2*T2)))/mu2
      tmp226 = -(phiT2bmu2/delT2Tbmu2)-
     .   (2d0*phiT2bmu2*(b+mu2-T2)*(-b-mu2+T2))/
     .    delT2Tbmu2**2-(mu2*tmp102)/(delT2Tbmu2*T2**2)-
     .   (2d0*mu2*(-b-mu2+T2)*tmp102)/(delT2Tbmu2**2*T2)+
     .   ((b+mu2-T2)*((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .      (mu2*tmp102)/(delT2Tbmu2*T2)))/delT2Tbmu2+
     .   (mu2*tmp5)/(delT2Tbmu2*T2)
      tmp227 = phiT2bmu2/delT2Tbmu2-
     .   (2d0*phiT2bmu2*(b-mu2-T2)*(b+mu2-T2))/
     .    delT2Tbmu2**2-(2d0*mu2*(b-mu2-T2)*tmp102)/
     .    (delT2Tbmu2**2*T2)+
     .   ((b+mu2-T2)*((phiT2bmu2*(-b+mu2+T2))/delT2Tbmu2+
     .      (mu2*tmp110)/(b*delT2Tbmu2)))/delT2Tbmu2+
     .   (mu2*tmp97)/(delT2Tbmu2*T2)
      tmp228 = -0.5d0+2d0*Logt-
     .   (phimu2T2T*(-mu2+t-T2))/T2-0.5d0*(Logmu2*Logt)+
     .   0.5d0*(Logmu2*LogT2)-0.5d0*(Logt*LogT2)+
     .   0.5d0*(LogT2*(mu2-t-T2))/t+
     .   0.5d0*(Logmu2*(-mu2-t+T2))/t-
     .   0.5d0*(deltmu2T2T*((mu2*phimu2T2T*(mu2-t+T2))/
     .        (deltmu2T2T*T2)+(mu2*tmp111)/(deltmu2T2T*t)))/mu2
      tmp229 = -0.5d0+2d0*Logt-(phiA0bt*(-A0-b+t))/t+
     .   0.5d0*(LogA0*Logb)-0.5d0*(LogA0*Logt)-0.5d0*(Logb*Logt)+
     .   0.5d0*(Logb*(A0-b-t))/t+0.5d0*(LogA0*(-A0+b-t))/t-
     .   0.5d0*(deltA0bt*((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .       (b*tmp66)/(deltA0bt*t)))/b
      tmp230 = -((b*phiA0bt)/(deltA0bt*t))-
     .   (2d0*b*phiA0bt*(A0+b-t)*(-A0-b+t))/
     .    (deltA0bt**2*t)+(b*tmp1)/(deltA0bt*t)-
     .   (b*tmp66)/(deltA0bt*t**2)-
     .   (2d0*b*(-A0-b+t)*tmp66)/(deltA0bt**2*t)+
     .   ((A0+b-t)*((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .      (b*tmp66)/(deltA0bt*t)))/deltA0bt
      tmp231 = -0.5d0+2d0*Logt-
     .   (phiB1tmu2*(-B1-mu2+t))/mu2+0.5d0*(LogB1*Logmu2)-
     .   0.5d0*(LogB1*Logt)-0.5d0*(Logmu2*Logt)+
     .   0.5d0*(Logmu2*(B1-mu2-t))/t+
     .   0.5d0*(LogB1*(-B1+mu2-t))/t-
     .   0.5d0*(deltB1tmu2*((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .       (mu2*tmp67)/(deltB1tmu2*t)))/mu2
      tmp232 = -(phiB1tmu2/deltB1tmu2)-
     .   (2d0*phiB1tmu2*(B1+mu2-t)*(-B1-mu2+t))/
     .    deltB1tmu2**2+(mu2*tmp2)/(deltB1tmu2*t)-
     .   (mu2*tmp67)/(deltB1tmu2*t**2)-
     .   (2d0*mu2*(-B1-mu2+t)*tmp67)/(deltB1tmu2**2*t)+
     .   ((B1+mu2-t)*((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .      (mu2*tmp67)/(deltB1tmu2*t)))/deltB1tmu2
      tmp233 = phiB1tmu2/deltB1tmu2-
     .   (2d0*phiB1tmu2*(-B1-mu2+t)*(-B1+mu2+t))/
     .    deltB1tmu2**2+((-B1+mu2+t)*
     .      ((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .      (mu2*tmp67)/(deltB1tmu2*t)))/deltB1tmu2-
     .   (2d0*mu2*(-B1-mu2+t)*tmp71)/(B1*deltB1tmu2**2)+
     .   (mu2*tmp73)/(B1*deltB1tmu2)
      tmp234 = -0.5d0+2d0*Logt-
     .   (phiB2tmu2*(-B2-mu2+t))/mu2+0.5d0*(LogB2*Logmu2)-
     .   0.5d0*(LogB2*Logt)-0.5d0*(Logmu2*Logt)+
     .   0.5d0*(Logmu2*(B2-mu2-t))/t+
     .   0.5d0*(LogB2*(-B2+mu2-t))/t-
     .   0.5d0*(deltB2tmu2*((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .       (mu2*tmp68)/(deltB2tmu2*t)))/mu2
      tmp235 = -((phiB1tmu2*(-B1-mu2+t))/mu2)+
     .   (phiB2tmu2*(-B2-mu2+t))/mu2+0.5d0*(LogB1*Logmu2)-
     .   0.5d0*(LogB2*Logmu2)-0.5d0*(LogB1*Logt)+
     .   0.5d0*(LogB2*Logt)+0.5d0*(Logmu2*(B1-mu2-t))/t-
     .   0.5d0*(Logmu2*(B2-mu2-t))/t+
     .   0.5d0*(LogB1*(-B1+mu2-t))/t-
     .   0.5d0*(LogB2*(-B2+mu2-t))/t-
     .   0.5d0*(deltB1tmu2*((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .       (mu2*tmp67)/(deltB1tmu2*t)))/mu2+
     .   0.5d0*(deltB2tmu2*((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .       (mu2*tmp68)/(deltB2tmu2*t)))/mu2
      tmp236 = -(phiB2tmu2/deltB2tmu2)-
     .   (2d0*phiB2tmu2*(B2+mu2-t)*(-B2-mu2+t))/
     .    deltB2tmu2**2+(mu2*tmp3)/(deltB2tmu2*t)-
     .   (mu2*tmp68)/(deltB2tmu2*t**2)-
     .   (2d0*mu2*(-B2-mu2+t)*tmp68)/(deltB2tmu2**2*t)+
     .   ((B2+mu2-t)*((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .      (mu2*tmp68)/(deltB2tmu2*t)))/deltB2tmu2
      tmp237 = phiB2tmu2/deltB2tmu2-
     .   (2d0*phiB2tmu2*(-B2-mu2+t)*(-B2+mu2+t))/
     .    deltB2tmu2**2+((-B2+mu2+t)*
     .      ((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .      (mu2*tmp68)/(deltB2tmu2*t)))/deltB2tmu2-
     .   (2d0*mu2*(-B2-mu2+t)*tmp72)/(B2*deltB2tmu2**2)+
     .   (mu2*tmp74)/(B2*deltB2tmu2)
      tmp238 = -0.5d0+2d0*LogT1-
     .   (phiA0B1T1*(-A0-B1+T1))/T1+0.5d0*(LogA0*LogB1)-
     .   0.5d0*(LogA0*LogT1)-0.5d0*(LogB1*LogT1)+
     .   0.5d0*(LogB1*(A0-B1-T1))/T1+
     .   0.5d0*(LogA0*(-A0+B1-T1))/T1-
     .   0.5d0*(deltA0B1T1*((B1*phiA0B1T1*(A0+B1-T1))/
     .        (deltA0B1T1*T1)+(B1*tmp81)/(deltA0B1T1*T1)))/B1
      tmp239 = -0.5d0+2d0*LogT1-
     .   (phiA0B2T1*(-A0-B2+T1))/T1+0.5d0*(LogA0*LogB2)-
     .   0.5d0*(LogA0*LogT1)-0.5d0*(LogB2*LogT1)+
     .   0.5d0*(LogB2*(A0-B2-T1))/T1+
     .   0.5d0*(LogA0*(-A0+B2-T1))/T1-
     .   0.5d0*(deltA0B2T1*((B2*phiA0B2T1*(A0+B2-T1))/
     .        (deltA0B2T1*T1)+(B2*tmp82)/(deltA0B2T1*T1)))/B2
      tmp240 = -0.5d0+2d0*LogT1-
     .   (phiT1bmu2*(-b-mu2+T1))/mu2+0.5d0*(Logb*Logmu2)-
     .   0.5d0*(Logb*LogT1)-0.5d0*(Logmu2*LogT1)+
     .   0.5d0*(Logmu2*(b-mu2-T1))/T1+
     .   0.5d0*(Logb*(-b+mu2-T1))/T1-
     .   0.5d0*(deltT1bmu2*((phiT1bmu2*(b+mu2-T1))/deltT1bmu2+
     .       (mu2*tmp83)/(deltT1bmu2*T1)))/mu2
      tmp241 = -(phiT1bmu2/deltT1bmu2)-
     .   (2d0*phiT1bmu2*(b+mu2-T1)*(-b-mu2+T1))/
     .    deltT1bmu2**2+(mu2*tmp4)/(deltT1bmu2*T1)-
     .   (mu2*tmp83)/(deltT1bmu2*T1**2)-
     .   (2d0*mu2*(-b-mu2+T1)*tmp83)/(deltT1bmu2**2*T1)+
     .   ((b+mu2-T1)*((phiT1bmu2*(b+mu2-T1))/deltT1bmu2+
     .      (mu2*tmp83)/(deltT1bmu2*T1)))/deltT1bmu2
      tmp242 = phiT1bmu2/deltT1bmu2-
     .   (2d0*phiT1bmu2*(b-mu2-T1)*(b+mu2-T1))/
     .    deltT1bmu2**2+(mu2*tmp78)/(deltT1bmu2*T1)-
     .   (2d0*mu2*(b-mu2-T1)*tmp83)/(deltT1bmu2**2*T1)+
     .   ((b+mu2-T1)*((phiT1bmu2*(-b+mu2+T1))/deltT1bmu2+
     .      (mu2*tmp90)/(b*deltT1bmu2)))/deltT1bmu2
      tmp243 = -0.5d0+2d0*Logt-
     .   (phimu2tT1*(-mu2+t-T1))/T1-0.5d0*(Logmu2*Logt)+
     .   0.5d0*(Logmu2*LogT1)-0.5d0*(Logt*LogT1)+
     .   0.5d0*(LogT1*(mu2-t-T1))/t+
     .   0.5d0*(Logmu2*(-mu2-t+T1))/t-
     .   0.5d0*(deltmu2tT1*((mu2*phimu2tT1*(mu2-t+T1))/
     .        (deltmu2tT1*T1)+(mu2*tmp91)/(deltmu2tT1*t)))/mu2
      tmp244 = 2d0/T1-Logb/T1-Logmu2/T1-
     .   0.5d0*(Logmu2*(b-mu2-T1))/T1**2-
     .   0.5d0*(Logb*(-b+mu2-T1))/T1**2-
     .   0.5d0*(2d0*phiT1bmu2+deltT1bmu2*tmp241+
     .       4d0*(-b-mu2+T1)*
     .      ((phiT1bmu2*(b+mu2-T1))/deltT1bmu2+
     .        (mu2*tmp83)/(deltT1bmu2*T1)))/mu2
      D1t = -(tmp109*tmp195)-tmp108*tmp196-
     .   2d0*hb*ht*mt*mu*s2b*tmp235-(hb*ht*mu*s2b*tmp75)/mt+
     .   tmp61*(-(B1*(-1d0+LogB1))+2d0*B1*LogB1-
     .      B1*(-1d0+LogB1)*(-1d0+Logt)+(-1d0+Logmu2)*mu2+
     .      2d0*Logmu2*mu2+(-1d0+Logmu2)*(-1d0+Logt)*mu2+
     .      2d0*Logt*t-(B1-mu2-t)*tmp231-
     .      0.5d0*(deltB1tmu2*phiB1tmu2)/mu2+
     .      0.5d0*(Logmu2*Logt*(B1-mu2-t))+
     .      0.5d0*(LogB1*Logt*(-B1+mu2-t))+
     .      0.5d0*(LogB1*Logmu2*(-B1-mu2+t))-2.5d0*(B1+mu2+t)
     .      )+tmp60*(-(B2*(-1d0+LogB2))+2d0*B2*LogB2-
     .      B2*(-1d0+LogB2)*(-1d0+Logt)+(-1d0+Logmu2)*mu2+
     .      2d0*Logmu2*mu2+(-1d0+Logmu2)*(-1d0+Logt)*mu2+
     .      2d0*Logt*t-(B2-mu2-t)*tmp234-
     .      0.5d0*(deltB2tmu2*phiB2tmu2)/mu2+
     .      0.5d0*(Logmu2*Logt*(B2-mu2-t))+
     .      0.5d0*(LogB2*Logt*(-B2+mu2-t))+
     .      0.5d0*(LogB2*Logmu2*(-B2-mu2+t))-2.5d0*(B2+mu2+t)
     .      )-0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-
     .      5d0*T1+LogT1*(-2d0*B2*LogB2+4d0*T1)-
     .      2d0*(-B2+T1)*tmp170)*tmp183)-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      tmp184)-0.5d0*((-5d0*B2+4d0*B2*LogB2+
     .      LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      tmp187)-0.5d0*((-5d0*B1+4d0*B1*LogB1+
     .      LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      tmp188)-4d0*cbe*hb*ht*mb*mt*sbe*
     .    (0.5d0-2d0*Logt+(phiA0bt*(-A0-b+t))/t-
     .      0.5d0*(LogA0*Logb)+0.5d0*(LogA0*Logt)+
     .      0.5d0*(Logb*Logt)-0.5d0*(Logb*(A0-b-t))/t-
     .      0.5d0*(LogA0*(-A0+b-t))/t+0.5d0*tmp210+
     .      0.5d0*(deltA0bt*((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .          (b*tmp66)/(deltA0bt*t)))/b)-
     .   (2d0*cbe*hb*ht*mb*sbe*
     .      (-2d0*A0*LogA0-2d0*b*Logb-2d0*Logt*t-
     .      0.5d0*(Logb*Logt*(A0-b-t))-
     .      0.5d0*(LogA0*Logt*(-A0+b-t))+
     .      0.5d0*(deltA0bt*phiA0bt)/t-
     .      0.5d0*(LogA0*Logb*(-A0-b+t))+2.5d0*(A0+b+t)+
     .      0.5d0*tmp76))/mt+
     .   tmp10*(b*(-1d0+Logb)+b*(-1d0+Logb)*(-1d0+Logt)+
     .      0.5d0*((b+t)*tmp210)+0.5d0*tmp76)
      D1t = D1t-tmp192*tmp88-tmp191*tmp89+
     .   tmp9*(-(A0*(-1d0+LogA0))+2d0*A0*LogA0+b*(-1d0+Logb)+
     .      2d0*b*Logb-A0*(-1d0+LogA0)*(-1d0+Logt)+
     .      b*(-1d0+Logb)*(-1d0+Logt)+2d0*Logt*t-
     .      (A0-b-t)*tmp229+0.5d0*(Logb*Logt*(A0-b-t))+
     .      0.5d0*(LogA0*Logt*(-A0+b-t))-
     .      0.5d0*(deltA0bt*phiA0bt)/t+
     .      0.5d0*(LogA0*Logb*(-A0-b+t))-2.5d0*(A0+b+t))+
     .   ht**2*(2d0*(-1d0+Logmu2)*mu2+4d0*Logmu2*mu2+
     .      2d0*(-1d0+Logmu2)*(-1d0+Logt)*mu2+4d0*Logt*t-
     .      (-1d0+LogT1)*T1-(-1d0+Logt)*(-1d0+LogT1)*T1+
     .      2d0*LogT1*T1-(-1d0+LogT2)*T2-
     .      (-1d0+Logt)*(-1d0+LogT2)*T2+2d0*LogT2*T2-
     .      (-mu2-t+T2)*tmp228-(-mu2-t+T1)*tmp243+
     .      sbe**2*(-10d0*t+2d0*(-1d0+Logt)*t+2d0*(-1d0+Logt)**2*t+
     .       4d0*Logt*t+Logt*(4d0*t-2d0*Logt*t)+t*tmp77)+
     .      0.5d0*(Logt*LogT1*(mu2-t-T1))+
     .      0.5d0*(Logmu2*LogT1*(-mu2+t-T1))-
     .      0.5d0*(deltmu2tT1*phimu2tT1)/T1+
     .      0.5d0*(Logmu2*Logt*(-mu2-t+T1))-
     .      2.5d0*(mu2+t+T1)+0.5d0*(Logt*LogT2*(mu2-t-T2))+
     .      0.5d0*(Logmu2*LogT2*(-mu2+t-T2))-
     .      0.5d0*(deltmu2T2T*phimu2T2T)/T2+
     .      0.5d0*(Logmu2*Logt*(-mu2-t+T2))-
     .      2.5d0*(mu2+t+T2)+
     .      cbe**2*(-2d0*A0*(-1d0+LogA0)-
     .       2d0*A0*(-1d0+LogA0)*(-1d0+Logt)+2d0*(-1d0+Logt)*t+
     .       2d0*(-1d0+Logt)**2*t+
     .       2d0*(2d0*A0*LogA0-A0*LogA0*Logt+4d0*Logt*t+
     .          0.5d0*(Logt**2*(A0-2d0*t))-
     .          0.5d0*(deltA0tt*phiA0tt)/t-2.5d0*(A0+2d0*t))-
     .       (A0-2d0*t)*(-1d0+4d0*Logt-Logt**2-(A0*LogA0)/t+
     .          (2d0*A0*phiA0tt)/t+(Logt*(A0-2d0*t))/t+
     .          0.5d0*(deltA0tt*phiA0tt)/t**2-
     .          0.5d0*(deltA0tt*
     .            ((A0*phiA0tt)/deltA0tt+
     .              (phiA0tt*tmp175)/deltA0tt+
     .              tmp65/deltA0tt+tmp70/deltA0tt))/t)))-
     .   0.25d0*(ht**2*(cbe**2*tmp114*(4d0-(2d0*s2t*Yt)/mt)+
     .      cbe**2*tmp93*(4d0+(2d0*s2t*Yt)/mt)+
     .      0.5d0*(sbe**2*tmp115*(4d0-(2d0*s2t*Xt)/mt))+
     .      0.5d0*(sbe**2*tmp94*(4d0+(2d0*s2t*Xt)/mt))))
      DT1 = -(tmp194*tmp238)-tmp193*tmp239-
     .   2d0*hb*ht*mb*mu*s2t*tmp240+
     .   hb**2*(0.25d0*(B1*(1d0-c2b)*(1d0+c2t)*(-1d0+LogB1))+
     .      0.25d0*(B2*(1d0+c2b)*(1d0+c2t)*(-1d0+LogB2))+
     .      sbe**2*(0.5d0*(A0*(1d0+c2t)*(-1d0+LogA0))+
     .       0.5d0*(A0*(1d0+c2t)*(-1d0+LogA0)*(-1d0+LogT1)))+
     .      0.25d0*(B1*(1d0-c2b)*(1d0+c2t)*(-1d0+LogB1)*
     .       (-1d0+LogT1))+
     .      0.25d0*(B2*(1d0+c2b)*(1d0+c2t)*(-1d0+LogB2)*(-1d0+LogT1))
     .      )+tmp62*(-(b*(-1d0+Logb))-2d0*b*Logb-
     .      b*(-1d0+Logb)*(-1d0+LogT1)-(-1d0+Logmu2)*mu2-
     .      2d0*Logmu2*mu2-(-1d0+Logmu2)*(-1d0+LogT1)*mu2-
     .      2d0*LogT1*T1-(-b-mu2+T1)*tmp240+
     .      0.5d0*(deltT1bmu2*phiT1bmu2)/mu2-
     .      0.5d0*(Logmu2*LogT1*(b-mu2-T1))-
     .      0.5d0*(Logb*LogT1*(-b+mu2-T1))-
     .      0.5d0*(Logb*Logmu2*(-b-mu2+T1))+2.5d0*(b+mu2+T1))
     .    -0.5d0*(tmp186*tmp213)-0.5d0*(tmp185*tmp214)+
     .   ht**2*((-1d0+LogT2)*T2*tmp8+
     .      (-1d0+LogT1)*(-1d0+LogT2)*T2*tmp8+
     .      cbe**2*(A0*(-1d0+LogA0)*(1d0+0.5d0*(1d0-c2t))+
     .       A0*(-1d0+LogA0)*(-1d0+LogT1)*(1d0+0.5d0*(1d0-c2t)))+
     .      0.25d0*(B1*(1d0+c2b)*(1d0-c2t)*(-1d0+LogB1))+
     .      0.25d0*(B2*(1d0-c2b)*(1d0-c2t)*(-1d0+LogB2))+
     .      0.25d0*(B1*(1d0+c2b)*(1d0-c2t)*(-1d0+LogB1)*
     .       (-1d0+LogT1))+
     .      0.25d0*(B2*(1d0-c2b)*(1d0-c2t)*(-1d0+LogB2)*
     .       (-1d0+LogT1))+0.25d0*((1d0+Nc)*s2t**2*tmp79))+
     .   ht**2*(-((-1d0+Logmu2)*mu2)-2d0*Logmu2*mu2-
     .      (-1d0+Logmu2)*(-1d0+LogT1)*mu2-(-1d0+Logt)*t-
     .      2d0*Logt*t-(-1d0+Logt)*(-1d0+LogT1)*t-2d0*LogT1*T1-
     .      0.5d0*(Logt*LogT1*(mu2-t-T1))-
     .      0.5d0*(Logmu2*LogT1*(-mu2+t-T1))+
     .      0.5d0*(deltmu2tT1*phimu2tT1)/T1-
     .      0.5d0*(Logmu2*Logt*(-mu2-t+T1))+
     .      2.5d0*(mu2+t+T1)-
     .      (-mu2-t+T1)*
     .       (-0.5d0+2d0*LogT1-(phimu2tT1*(-mu2-t+T1))/T1+
     .       0.5d0*(Logmu2*Logt)-0.5d0*(Logmu2*LogT1)-
     .       0.5d0*(Logt*LogT1)+0.5d0*(Logt*(mu2-t-T1))/T1+
     .       0.5d0*(Logmu2*(-mu2+t-T1))/T1-
     .       0.5d0*(deltmu2tT1*
     .           ((mu2*phimu2tT1*(mu2+t-T1))/
     .            (deltmu2tT1*T1)+(mu2*tmp84)/(deltmu2tT1*T1)
     .             ))/mu2))
      DT1 = DT1-0.25d0*
     .    (ht**2*((1d0+c2t**2)*sbe**2*tmp220*Xt**2+
     .      2d0*(1d0+c2t**2)*cbe**2*tmp221*Yt**2+
     .      cbe**2*tmp51*
     .       (-1d0+4d0*LogT1-LogT1**2-(A0*LogA0)/T1+
     .         (2d0*A0*phiA0T1T1)/T1+(LogT1*(A0-2d0*T1))/T1+
     .         0.5d0*(deltA0T1T1*phiA0T1T1)/T1**2-
     .         0.5d0*(deltA0T1T1*
     .             ((A0*phiA0T1T1)/deltA0T1T1+
     .             (phiA0T1T1*tmp178)/deltA0T1T1+
     .             tmp80/deltA0T1T1+tmp87/deltA0T1T1))/T1)+
     .      0.5d0*(sbe**2*tmp35*tmp95)))
      DT2 = -(tmp198*tmp222)-tmp197*tmp223-
     .   2d0*hb*ht*mb*mu*s2t*tmp225+
     .   hb**2*(0.25d0*(B1*(1d0-c2b)*(1d0-c2t)*(-1d0+LogB1))+
     .      0.25d0*(B2*(1d0+c2b)*(1d0-c2t)*(-1d0+LogB2))+
     .      sbe**2*(0.5d0*(A0*(1d0-c2t)*(-1d0+LogA0))+
     .       0.5d0*(A0*(1d0-c2t)*(-1d0+LogA0)*(-1d0+LogT2)))+
     .      0.25d0*(B1*(1d0-c2b)*(1d0-c2t)*(-1d0+LogB1)*
     .       (-1d0+LogT2))+
     .      0.25d0*(B2*(1d0+c2b)*(1d0-c2t)*(-1d0+LogB2)*(-1d0+LogT2))
     .      )+tmp63*(-(b*(-1d0+Logb))-2d0*b*Logb-
     .      b*(-1d0+Logb)*(-1d0+LogT2)-(-1d0+Logmu2)*mu2-
     .      2d0*Logmu2*mu2-(-1d0+Logmu2)*(-1d0+LogT2)*mu2-
     .      2d0*LogT2*T2-(-b-mu2+T2)*tmp224+
     .      0.5d0*(delT2Tbmu2*phiT2bmu2)/mu2-
     .      0.5d0*(Logmu2*LogT2*(b-mu2-T2))-
     .      0.5d0*(Logb*LogT2*(-b+mu2-T2))-
     .      0.5d0*(Logb*Logmu2*(-b-mu2+T2))+2.5d0*(b+mu2+T2))
     .    +ht**2*(-((-1d0+Logmu2)*mu2)-2d0*Logmu2*mu2-
     .      (-1d0+Logmu2)*(-1d0+LogT2)*mu2-(-1d0+Logt)*t-
     .      2d0*Logt*t-(-1d0+Logt)*(-1d0+LogT2)*t-2d0*LogT2*T2-
     .      0.5d0*(Logt*LogT2*(mu2-t-T2))-
     .      0.5d0*(Logmu2*LogT2*(-mu2+t-T2))+
     .      0.5d0*(deltmu2T2T*phimu2T2T)/T2-
     .      0.5d0*(Logmu2*Logt*(-mu2-t+T2))+
     .      2.5d0*(mu2+t+T2)-
     .      (-mu2-t+T2)*
     .       (-0.5d0+2d0*LogT2-(phimu2T2T*(-mu2-t+T2))/T2+
     .       0.5d0*(Logmu2*Logt)-0.5d0*(Logmu2*LogT2)-
     .       0.5d0*(Logt*LogT2)+0.5d0*(Logt*(mu2-t-T2))/T2+
     .       0.5d0*(Logmu2*(-mu2+t-T2))/T2-
     .       0.5d0*(deltmu2T2T*
     .           ((mu2*phimu2T2T*(mu2+t-T2))/
     .            (deltmu2T2T*T2)+
     .             (mu2*tmp103)/(deltmu2T2T*T2)))/mu2))-
     .   0.5d0*(tmp190*tmp217)-0.5d0*(tmp189*tmp218)+
     .   ht**2*((-1d0+LogT1)*T1*tmp8+
     .      (-1d0+LogT1)*(-1d0+LogT2)*T1*tmp8+
     .      cbe**2*(A0*(-1d0+LogA0)*(1d0+0.5d0*(1d0+c2t))+
     .       A0*(-1d0+LogA0)*(-1d0+LogT2)*(1d0+0.5d0*(1d0+c2t)))+
     .      0.25d0*(B1*(1d0+c2b)*(1d0+c2t)*(-1d0+LogB1))+
     .      0.25d0*(B2*(1d0-c2b)*(1d0+c2t)*(-1d0+LogB2))+
     .      0.25d0*(B1*(1d0+c2b)*(1d0+c2t)*(-1d0+LogB1)*
     .       (-1d0+LogT2))+
     .      0.25d0*(B2*(1d0-c2b)*(1d0+c2t)*(-1d0+LogB2)*
     .       (-1d0+LogT2))+0.25d0*((1d0+Nc)*s2t**2*tmp98))
      DT2 = DT2-0.25d0*
     .    (ht**2*((1d0+c2t**2)*sbe**2*tmp219*Xt**2+
     .      2d0*(1d0+c2t**2)*cbe**2*Yt**2*
     .       (-0.5d0+2d0*LogT2-(phiA0T1T2*(-A0-T1+T2))/T2+
     .         0.5d0*(LogA0*LogT1)-0.5d0*(LogA0*LogT2)-
     .         0.5d0*(LogT1*LogT2)+
     .         0.5d0*(deltA0T1T2*phiA0T1T2)/T2**2+
     .         0.5d0*(LogT1*(A0-T1-T2))/T2+
     .         0.5d0*(LogA0*(-A0+T1-T2))/T2-
     .         0.5d0*(deltA0T1T2*
     .             (tmp104/deltA0T1T2+
     .             (phiA0T1T2*tmp179)/deltA0T1T2))/T2)+
     .      0.5d0*(sbe**2*tmp116*tmp34)+
     .      cbe**2*tmp50*
     .       (-1d0+4d0*LogT2-LogT2**2-(A0*LogA0)/T2+
     .         (2d0*A0*phiA0T2T2)/T2+(LogT2*(A0-2d0*T2))/T2+
     .         0.5d0*(deltA0T2T2*phiA0T2T2)/T2**2-
     .         0.5d0*(deltA0T2T2*
     .             ((A0*phiA0T2T2)/deltA0T2T2+
     .             tmp107/deltA0T2T2+
     .             (phiA0T2T2*tmp182)/deltA0T2T2+
     .             tmp99/deltA0T2T2))/T2)))
      Dc2t = (hb*ht*mb*mu*tmp113)/s2t-tmp109*tmp205-
     .   tmp108*tmp206+(b*(-1d0+Logb)*(-1d0+Logmu2)*mu2-
     .      b*(-1d0+Logb)*(-1d0+LogT2)*T2-
     .      (-1d0+Logmu2)*(-1d0+LogT2)*mu2*T2-
     .      (-b-mu2+T2)*tmp112)*tmp7-tmp204*tmp88-
     .   tmp203*tmp89+tmp6*
     .    (b*(-1d0+Logb)*(-1d0+Logmu2)*mu2-
     .      b*(-1d0+Logb)*(-1d0+LogT1)*T1-
     .      (-1d0+Logmu2)*(-1d0+LogT1)*mu2*T1-
     .      (-b-mu2+T1)*tmp92)+
     .   hb**2*(0.125d0*(B1*(1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT1)*T1)/
     .      c2t+0.125d0*(B2*(1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT1)*
     .        T1)/c2t+sbe**2*
     .       (0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT1)*T1)/c2t-
     .       0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT2)*T2)/c2t)-
     .      0.125d0*(B1*(1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT2)*T2)/
     .      c2t-0.125d0*(B2*(1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT2)*
     .        T2)/c2t)+
     .   ht**2*((-1d0+LogT1)*(-1d0+LogT2)*T1*T2*
     .       (1d0+0.5d0*(-1d0+Nc))-
     .      0.125d0*(B1*(1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT1)*T1)/
     .      c2t-0.125d0*(B2*(1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT1)*
     .        T1)/c2t+cbe**2*
     .       (-(0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT1)*T1)/c2t)+
     .       0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT2)*T2)/c2t)+
     .      0.125d0*(B1*(1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT2)*T2)/
     .      c2t+0.125d0*(B2*(1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT2)*
     .        T2)/c2t-0.25d0*
     .       ((1d0+Nc)*((-1d0+LogT1)**2*T1**2+
     .         (-1d0+LogT2)**2*T2**2)))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      tmp199)-0.5d0*((-5d0*B1+4d0*B1*LogB1+
     .      LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      tmp200)-0.5d0*((-5d0*B2+4d0*B2*LogB2+
     .      LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      tmp201)-0.5d0*((-5d0*B1+4d0*B1*LogB1+
     .      LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      tmp202)
      Dc2t = Dc2t-0.25d0*
     .    (ht**2*(cbe**2*tmp114*tmp47+cbe**2*tmp46*tmp93+
     .      sbe**2*(-5d0*T1-5d0*T2+4d0*LogT2*T2+
     .         LogT1**2*(-T1+T2)+LogT1*(4d0*T1-2d0*LogT2*T2)-
     .         2d0*(T1-T2)*tmp171)*Xt**2+
     .      2d0*cbe**2*Yt**2*
     .       (2d0*A0*LogA0+2d0*LogT1*T1+2d0*LogT2*T2+
     .         0.5d0*(LogT1*LogT2*(A0-T1-T2))+
     .         0.5d0*(LogA0*LogT2*(-A0+T1-T2))-
     .         0.5d0*(deltA0T1T2*phiA0T1T2)/T2+
     .         0.5d0*(LogA0*LogT1*(-A0-T1+T2))-
     .         2.5d0*(A0+T1+T2))+0.5d0*(sbe**2*tmp115*tmp31)+
     .      0.5d0*(sbe**2*tmp30*tmp94)))
      DT1T1 = -2d0*hb*ht*mb*mu*s2t*tmp244+
     .   hb**2*(0.25d0*(B1*(1d0-c2b)*(1d0+c2t)*(-1d0+LogB1))/T1+
     .      0.25d0*(B2*(1d0+c2b)*(1d0+c2t)*(-1d0+LogB2))/T1+
     .      0.5d0*(A0*(1d0+c2t)*(-1d0+LogA0)*sbe**2)/T1)+
     .   ht**2*(((-1d0+LogT2)*T2*tmp8)/T1+
     .      (A0*cbe**2*(-1d0+LogA0)*(1d0+0.5d0*(1d0-c2t)))/T1+
     .      0.25d0*(B1*(1d0+c2b)*(1d0-c2t)*(-1d0+LogB1))/T1+
     .      0.25d0*(B2*(1d0-c2b)*(1d0-c2t)*(-1d0+LogB2))/T1+
     .      0.25d0*((1d0+Nc)*s2t**2*
     .       (8d0*(-1d0+LogT1)+2d0*(-1d0+LogT1)**2+
     .         (2d0/T1**2-(2d0*(-1d0+LogT1))/T1**2)*T1**2)))-
     .   0.5d0*(((4d0*B2*(LogB2-LogT1))/((1d0-B2/T1)*T1**2)+8d0/T1-
     .      (4d0*LogT1)/T1-
     .      2d0*((B2**2*(LogB2-LogT1))/((1d0-B2/T1)**2*T1**4)+
     .         B2/((1d0-B2/T1)*T1**3)+
     .         (2d0*B2*(LogB2-LogT1))/((1d0-B2/T1)*T1**3))*
     .       (-B2+T1)-(-2d0*B2*LogB2+4d0*T1)/T1**2+
     .      (B2-T1)*tmp17)*tmp185)-
     .   0.5d0*(((4d0*B1*(LogB1-LogT1))/((1d0-B1/T1)*T1**2)+8d0/T1-
     .      (4d0*LogT1)/T1-
     .      2d0*((B1**2*(LogB1-LogT1))/((1d0-B1/T1)**2*T1**4)+
     .         B1/((1d0-B1/T1)*T1**3)+
     .         (2d0*B1*(LogB1-LogT1))/((1d0-B1/T1)*T1**3))*
     .       (-B1+T1)-(-2d0*B1*LogB1+4d0*T1)/T1**2+
     .      (B1-T1)*tmp17)*tmp186)-
     .   tmp194*(2d0/T1-LogA0/T1-LogB1/T1-
     .      0.5d0*(LogB1*(A0-B1-T1))/T1**2-
     .      0.5d0*(LogA0*(-A0+B1-T1))/T1**2-
     .      0.5d0*((2d0*B1*phiA0B1T1)/T1+
     .        4d0*(-A0-B1+T1)*
     .         ((B1*phiA0B1T1*(A0+B1-T1))/(deltA0B1T1*T1)+
     .           (B1*tmp81)/(deltA0B1T1*T1))+
     .        deltA0B1T1*
     .         ((B1*(2d0-LogA0-LogB1+2d0*LogT1))/
     .            (deltA0B1T1*T1)-
     .           (B1*phiA0B1T1)/(deltA0B1T1*T1)-
     .           (2d0*B1*phiA0B1T1*(A0+B1-T1)*(-A0-B1+T1))/
     .            (deltA0B1T1**2*T1)-
     .           (B1*tmp81)/(deltA0B1T1*T1**2)-
     .           (2d0*B1*(-A0-B1+T1)*tmp81)/
     .            (deltA0B1T1**2*T1)+
     .           ((A0+B1-T1)*
     .            ((B1*phiA0B1T1*(A0+B1-T1))/
     .               (deltA0B1T1*T1)+
     .              (B1*tmp81)/(deltA0B1T1*T1)))/deltA0B1T1))/
     .      B1)
      DT1T1 = DT1T1-tmp193*
     .    (2d0/T1-LogA0/T1-LogB2/T1-
     .      0.5d0*(LogB2*(A0-B2-T1))/T1**2-
     .      0.5d0*(LogA0*(-A0+B2-T1))/T1**2-
     .      0.5d0*((2d0*B2*phiA0B2T1)/T1+
     .        4d0*(-A0-B2+T1)*
     .         ((B2*phiA0B2T1*(A0+B2-T1))/(deltA0B2T1*T1)+
     .           (B2*tmp82)/(deltA0B2T1*T1))+
     .        deltA0B2T1*
     .         ((B2*(2d0-LogA0-LogB2+2d0*LogT1))/
     .            (deltA0B2T1*T1)-
     .           (B2*phiA0B2T1)/(deltA0B2T1*T1)-
     .           (2d0*B2*phiA0B2T1*(A0+B2-T1)*(-A0-B2+T1))/
     .            (deltA0B2T1**2*T1)-
     .           (B2*tmp82)/(deltA0B2T1*T1**2)-
     .           (2d0*B2*(-A0-B2+T1)*tmp82)/
     .            (deltA0B2T1**2*T1)+
     .           ((A0+B2-T1)*
     .            ((B2*phiA0B2T1*(A0+B2-T1))/
     .               (deltA0B2T1*T1)+
     .              (B2*tmp82)/(deltA0B2T1*T1)))/deltA0B2T1))/
     .      B2)+tmp62*(-((b*(-1d0+Logb))/T1)-
     .      ((-1d0+Logmu2)*mu2)/T1-(-b-mu2+T1)*tmp244+
     .      2d0*(0.5d0-2d0*LogT1+(phiT1bmu2*(-b-mu2+T1))/mu2-
     .       0.5d0*(Logb*Logmu2)+0.5d0*(Logb*LogT1)+
     .       0.5d0*(Logmu2*LogT1)-
     .       0.5d0*(Logmu2*(b-mu2-T1))/T1-
     .       0.5d0*(Logb*(-b+mu2-T1))/T1+
     .       0.5d0*(deltT1bmu2*
     .           ((phiT1bmu2*(b+mu2-T1))/deltT1bmu2+
     .             (mu2*tmp83)/(deltT1bmu2*T1)))/mu2))
      DT1T1 = DT1T1+ht**2*
     .    (-(((-1d0+Logmu2)*mu2)/T1)-((-1d0+Logt)*t)/T1+
     .      2d0*(0.5d0-2d0*LogT1+(phimu2tT1*(-mu2-t+T1))/T1-
     .       0.5d0*(Logmu2*Logt)+0.5d0*(Logmu2*LogT1)+
     .       0.5d0*(Logt*LogT1)-0.5d0*(Logt*(mu2-t-T1))/T1-
     .       0.5d0*(Logmu2*(-mu2+t-T1))/T1+
     .       0.5d0*(deltmu2tT1*
     .           ((mu2*phimu2tT1*(mu2+t-T1))/
     .            (deltmu2tT1*T1)+(mu2*tmp84)/(deltmu2tT1*T1)
     .             ))/mu2)-
     .      (-mu2-t+T1)*
     .       (2d0/T1-Logmu2/T1-Logt/T1-
     .       0.5d0*(Logt*(mu2-t-T1))/T1**2-
     .       0.5d0*(Logmu2*(-mu2+t-T1))/T1**2-
     .       0.5d0*((2d0*mu2*phimu2tT1)/T1+
     .           4d0*(-mu2-t+T1)*
     .            ((mu2*phimu2tT1*(mu2+t-T1))/
     .             (deltmu2tT1*T1)+
     .            (mu2*tmp84)/(deltmu2tT1*T1))+
     .           deltmu2tT1*
     .            (((2d0-Logmu2-Logt+2d0*LogT1)*mu2)/
     .             (deltmu2tT1*T1)-
     .            (mu2*phimu2tT1)/(deltmu2tT1*T1)-
     .            (2d0*mu2*phimu2tT1*(mu2+t-T1)*
     .               (-mu2-t+T1))/(deltmu2tT1**2*T1)-
     .            (mu2*tmp84)/(deltmu2tT1*T1**2)-
     .            (2d0*mu2*(-mu2-t+T1)*tmp84)/
     .             (deltmu2tT1**2*T1)+
     .            ((mu2+t-T1)*
     .               ((mu2*phimu2tT1*(mu2+t-T1))/
     .                  (deltmu2tT1*T1)+
     .                 (mu2*tmp84)/(deltmu2tT1*T1)))/deltmu2tT1
     .            ))/mu2))
      DT1T1 = DT1T1-0.25d0*
     .    (ht**2*((1d0+c2t**2)*sbe**2*
     .       (8d0/T1-(4d0*LogT1)/T1-(4d0*T1-2d0*LogT2*T2)/T1**2+
     .         (4d0*(-LogT1+LogT2)*T2)/(T1**2*(1d0-T2/T1))-
     .         2d0*(T1-T2)*
     .          (((-LogT1+LogT2)*T2**2)/
     .             (T1**4*(1d0-T2/T1)**2)+
     .            T2/(T1**3*(1d0-T2/T1))+
     .            (2d0*(-LogT1+LogT2)*T2)/(T1**3*(1d0-T2/T1)))+
     .         (-T1+T2)*tmp17)*Xt**2+
     .      0.5d0*(sbe**2*(4d0/T1+(2d0*(2d0-2d0*LogT1))/T1-
     .           (2d0*LogT1)/T1-(4d0*T1-2d0*LogT1*T1)/T1**2)*tmp35)
     .       +cbe**2*tmp51*
     .       (-((deltA0T1T1*phiA0T1T1)/T1**3)+
     .         (A0*LogA0)/T1**2+4d0/T1-(4d0*LogT1)/T1+
     .         (-4d0*A0*phiA0T1T1+
     .            deltA0T1T1*
     .             ((A0*phiA0T1T1)/deltA0T1T1+
     .             (phiA0T1T1*tmp178)/deltA0T1T1+
     .             tmp80/deltA0T1T1+tmp87/deltA0T1T1))/T1**2+
     .           0.5d0*((A0-2d0*T1)*tmp17)-
     .         0.5d0*(-8d0*A0*
     .            ((A0*phiA0T1T1)/deltA0T1T1+
     .              (phiA0T1T1*tmp178)/deltA0T1T1+
     .              tmp80/deltA0T1T1+tmp87/deltA0T1T1)+
     .             deltA0T1T1*
     .            ((4d0*A0**2*phiA0T1T1)/deltA0T1T1**2+
     .              (phiA0T1T1*
     .                 (-1d0-(A0-T1)**2/T1**2-
     .                   (2d0*(A0-T1))/T1))/deltA0T1T1+
     .              (1d0+(A0-T1)/T1)/deltA0T1T1+
     .              (1d0-(-A0+T1)/T1)/deltA0T1T1+
     .              (4d0*A0*phiA0T1T1*tmp178)/deltA0T1T1**2+
     .              (4d0*A0*tmp80)/deltA0T1T1**2+
     .              (4d0*A0*tmp87)/deltA0T1T1**2+
     .              (A0*
     .                 ((A0*phiA0T1T1)/deltA0T1T1+
     .                   (phiA0T1T1*tmp178)/deltA0T1T1+
     .                   tmp80/deltA0T1T1+tmp87/deltA0T1T1))/
     .               deltA0T1T1+
     .              (tmp178*
     .                 ((A0*phiA0T1T1)/deltA0T1T1+
     .                   (phiA0T1T1*tmp178)/deltA0T1T1+
     .                   tmp80/deltA0T1T1+tmp87/deltA0T1T1))/
     .               deltA0T1T1))/T1)+
     .      2d0*(1d0+c2t**2)*cbe**2*Yt**2*
     .       (2d0/T1-LogA0/T1-LogT2/T1-
     .         0.5d0*(LogT2*(A0-T1-T2))/T1**2-
     .         0.5d0*(LogA0*(-A0-T1+T2))/T1**2-
     .         0.5d0*(2d0*phiA0T1T2+
     .             4d0*(-A0+T1-T2)*
     .            ((phiA0T1T2*(A0-T1+T2))/deltA0T1T2+
     .              (T2*tmp96)/(deltA0T1T2*T1))+
     .             deltA0T1T2*
     .            (-(phiA0T1T2/deltA0T1T2)+
     .              ((2d0-LogA0+2d0*LogT1-LogT2)*T2)/
     .               (deltA0T1T2*T1)-
     .              (2d0*phiA0T1T2*(-A0+T1-T2)*
     .                 (A0-T1+T2))/deltA0T1T2**2-
     .              (T2*tmp96)/(deltA0T1T2*T1**2)-
     .              (2d0*(-A0+T1-T2)*T2*tmp96)/
     .               (deltA0T1T2**2*T1)+
     .              ((A0-T1+T2)*
     .                 ((phiA0T1T2*(A0-T1+T2))/
     .                  deltA0T1T2+(T2*tmp96)/(deltA0T1T2*T1)
     .                   ))/deltA0T1T2))/T2)))
      DT2T2 = hb**2*(0.25d0*
     .       (B1*(1d0-c2b)*(1d0-c2t)*(-1d0+LogB1))/T2+
     .      0.25d0*(B2*(1d0+c2b)*(1d0-c2t)*(-1d0+LogB2))/T2+
     .      0.5d0*(A0*(1d0-c2t)*(-1d0+LogA0)*sbe**2)/T2)+
     .   ht**2*(-(((-1d0+Logmu2)*mu2)/T2)-((-1d0+Logt)*t)/T2+
     .      2d0*(0.5d0-2d0*LogT2+(phimu2T2T*(-mu2-t+T2))/T2-
     .       0.5d0*(Logmu2*Logt)+0.5d0*(Logmu2*LogT2)+
     .       0.5d0*(Logt*LogT2)-0.5d0*(Logt*(mu2-t-T2))/T2-
     .       0.5d0*(Logmu2*(-mu2+t-T2))/T2+
     .       0.5d0*(deltmu2T2T*
     .           ((mu2*phimu2T2T*(mu2+t-T2))/
     .            (deltmu2T2T*T2)+
     .             (mu2*tmp103)/(deltmu2T2T*T2)))/mu2)-
     .      (-mu2-t+T2)*
     .       (2d0/T2-Logmu2/T2-Logt/T2-
     .       0.5d0*(Logt*(mu2-t-T2))/T2**2-
     .       0.5d0*(Logmu2*(-mu2+t-T2))/T2**2-
     .       0.5d0*((2d0*mu2*phimu2T2T)/T2+
     .           4d0*(-mu2-t+T2)*
     .            ((mu2*phimu2T2T*(mu2+t-T2))/
     .             (deltmu2T2T*T2)+
     .            (mu2*tmp103)/(deltmu2T2T*T2))+
     .           deltmu2T2T*
     .            (((2d0-Logmu2-Logt+2d0*LogT2)*mu2)/
     .             (deltmu2T2T*T2)-
     .            (mu2*phimu2T2T)/(deltmu2T2T*T2)-
     .            (2d0*mu2*phimu2T2T*(mu2+t-T2)*
     .               (-mu2-t+T2))/(deltmu2T2T**2*T2)-
     .            (mu2*tmp103)/(deltmu2T2T*T2**2)-
     .            (2d0*mu2*(-mu2-t+T2)*tmp103)/
     .             (deltmu2T2T**2*T2)+
     .            ((mu2+t-T2)*
     .               ((mu2*phimu2T2T*(mu2+t-T2))/
     .                  (deltmu2T2T*T2)+
     .                 (mu2*tmp103)/(deltmu2T2T*T2)))/
     .             deltmu2T2T))/mu2))-
     .   0.5d0*(((4d0*B2*(LogB2-LogT2))/((1d0-B2/T2)*T2**2)+8d0/T2-
     .      (4d0*LogT2)/T2-
     .      2d0*((B2**2*(LogB2-LogT2))/((1d0-B2/T2)**2*T2**4)+
     .         B2/((1d0-B2/T2)*T2**3)+
     .         (2d0*B2*(LogB2-LogT2))/((1d0-B2/T2)*T2**3))*
     .       (-B2+T2)-(-2d0*B2*LogB2+4d0*T2)/T2**2+
     .      (B2-T2)*tmp18)*tmp189)-
     .   0.5d0*(((4d0*B1*(LogB1-LogT2))/((1d0-B1/T2)*T2**2)+8d0/T2-
     .      (4d0*LogT2)/T2-
     .      2d0*((B1**2*(LogB1-LogT2))/((1d0-B1/T2)**2*T2**4)+
     .         B1/((1d0-B1/T2)*T2**3)+
     .         (2d0*B1*(LogB1-LogT2))/((1d0-B1/T2)*T2**3))*
     .       (-B1+T2)-(-2d0*B1*LogB1+4d0*T2)/T2**2+
     .      (B1-T2)*tmp18)*tmp190)
      DT2T2 = DT2T2-tmp198*
     .    (2d0/T2-LogA0/T2-LogB1/T2-
     .      0.5d0*(LogB1*(A0-B1-T2))/T2**2-
     .      0.5d0*(LogA0*(-A0+B1-T2))/T2**2-
     .      0.5d0*((2d0*B1*phiA0B1T2)/T2+
     .        4d0*(-A0-B1+T2)*
     .         ((B1*phiA0B1T2*(A0+B1-T2))/(deltA0B1T2*T2)+
     .           (B1*tmp100)/(deltA0B1T2*T2))+
     .        deltA0B1T2*
     .         ((B1*(2d0-LogA0-LogB1+2d0*LogT2))/
     .            (deltA0B1T2*T2)-
     .           (B1*phiA0B1T2)/(deltA0B1T2*T2)-
     .           (2d0*B1*phiA0B1T2*(A0+B1-T2)*(-A0-B1+T2))/
     .            (deltA0B1T2**2*T2)-
     .           (B1*tmp100)/(deltA0B1T2*T2**2)-
     .           (2d0*B1*(-A0-B1+T2)*tmp100)/
     .            (deltA0B1T2**2*T2)+
     .           ((A0+B1-T2)*
     .            ((B1*phiA0B1T2*(A0+B1-T2))/
     .               (deltA0B1T2*T2)+
     .              (B1*tmp100)/(deltA0B1T2*T2)))/deltA0B1T2))/
     .      B1)-tmp197*
     .    (2d0/T2-LogA0/T2-LogB2/T2-
     .      0.5d0*(LogB2*(A0-B2-T2))/T2**2-
     .      0.5d0*(LogA0*(-A0+B2-T2))/T2**2-
     .      0.5d0*((2d0*B2*phiA0B2T2)/T2+
     .        4d0*(-A0-B2+T2)*
     .         ((B2*phiA0B2T2*(A0+B2-T2))/(deltA0B2T2*T2)+
     .           (B2*tmp101)/(deltA0B2T2*T2))+
     .        deltA0B2T2*
     .         ((B2*(2d0-LogA0-LogB2+2d0*LogT2))/
     .            (deltA0B2T2*T2)-
     .           (B2*phiA0B2T2)/(deltA0B2T2*T2)-
     .           (2d0*B2*phiA0B2T2*(A0+B2-T2)*(-A0-B2+T2))/
     .            (deltA0B2T2**2*T2)-
     .           (B2*tmp101)/(deltA0B2T2*T2**2)-
     .           (2d0*B2*(-A0-B2+T2)*tmp101)/
     .            (deltA0B2T2**2*T2)+
     .           ((A0+B2-T2)*
     .            ((B2*phiA0B2T2*(A0+B2-T2))/
     .               (deltA0B2T2*T2)+
     .              (B2*tmp101)/(deltA0B2T2*T2)))/deltA0B2T2))/
     .      B2)+tmp63*(-((b*(-1d0+Logb))/T2)-
     .      ((-1d0+Logmu2)*mu2)/T2+2d0*tmp225-
     .      (-b-mu2+T2)*
     .       (2d0/T2-Logb/T2-Logmu2/T2-
     .       0.5d0*(Logmu2*(b-mu2-T2))/T2**2-
     .       0.5d0*(Logb*(-b+mu2-T2))/T2**2-
     .       0.5d0*(2d0*phiT2bmu2+
     .           4d0*(-b-mu2+T2)*
     .            ((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .            (mu2*tmp102)/(delT2Tbmu2*T2))+
     .           delT2Tbmu2*tmp226)/mu2))-
     .   2d0*hb*ht*mb*mu*s2t*
     .    (-2d0/T2+Logb/T2+Logmu2/T2+
     .      0.5d0*(Logmu2*(b-mu2-T2))/T2**2+
     .      0.5d0*(Logb*(-b+mu2-T2))/T2**2+
     .      0.5d0*(2d0*phiT2bmu2+
     .        4d0*(-b-mu2+T2)*
     .         ((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .           (mu2*tmp102)/(delT2Tbmu2*T2))+
     .        delT2Tbmu2*tmp226)/mu2)
      DT2T2 = DT2T2+ht**2*
     .    (((-1d0+LogT1)*T1*tmp8)/T2+
     .      (A0*cbe**2*(-1d0+LogA0)*(1d0+0.5d0*(1d0+c2t)))/T2+
     .      0.25d0*(B1*(1d0+c2b)*(1d0+c2t)*(-1d0+LogB1))/T2+
     .      0.25d0*(B2*(1d0-c2b)*(1d0+c2t)*(-1d0+LogB2))/T2+
     .      0.25d0*((1d0+Nc)*s2t**2*
     .       (8d0*(-1d0+LogT2)+2d0*(-1d0+LogT2)**2+
     .         (2d0/T2**2-(2d0*(-1d0+LogT2))/T2**2)*T2**2)))
      DT2T2 = DT2T2-0.25d0*
     .    (ht**2*((1d0+c2t**2)*sbe**2*
     .       (4d0/T2-(2d0*LogT1)/T2+
     .         (4d0*(-LogT1+LogT2))/(T1*(1d0-T2/T1))-
     .         2d0*(T1-T2)*
     .          ((-LogT1+LogT2)/(T1**2*(1d0-T2/T1)**2)+
     .            1/(T1*T2*(1d0-T2/T1))))*Xt**2+
     .      2d0*(1d0+c2t**2)*cbe**2*Yt**2*
     .       (-((deltA0T1T2*phiA0T1T2)/T2**3)+2d0/T2-
     .         LogA0/T2-LogT1/T2+
     .         (2d0*phiA0T1T2*(-A0-T1+T2)+
     .            deltA0T1T2*
     .             (tmp104/deltA0T1T2+
     .             (phiA0T1T2*tmp179)/deltA0T1T2))/T2**2-
     .         0.5d0*(LogT1*(A0-T1-T2))/T2**2-
     .         0.5d0*(LogA0*(-A0+T1-T2))/T2**2-
     .         0.5d0*(2d0*phiA0T1T2+
     .             4d0*(-A0-T1+T2)*
     .            (tmp104/deltA0T1T2+
     .              (phiA0T1T2*tmp179)/deltA0T1T2)+
     .             deltA0T1T2*
     .            ((2d0-LogA0-LogT1+2d0*LogT2)/deltA0T1T2-
     .              (phiA0T1T2*(A0-T1)**2)/
     .               (deltA0T1T2*T2**2)-
     .              (2d0*(-A0-T1+T2)*tmp104)/deltA0T1T2**2-
     .              (2d0*phiA0T1T2*(-A0-T1+T2)*tmp179)/
     .               deltA0T1T2**2+
     .              (tmp179*
     .                 (tmp104/deltA0T1T2+
     .                   (phiA0T1T2*tmp179)/deltA0T1T2))/
     .               deltA0T1T2))/T2)+
     .      0.5d0*(sbe**2*(4d0/T2+(2d0*(2d0-2d0*LogT2))/T2-
     .           (2d0*LogT2)/T2-(4d0*T2-2d0*LogT2*T2)/T2**2)*tmp34)
     .       +cbe**2*tmp50*
     .       (-((deltA0T2T2*phiA0T2T2)/T2**3)+
     .         (A0*LogA0)/T2**2+4d0/T2-(4d0*LogT2)/T2+
     .         (-4d0*A0*phiA0T2T2+
     .            deltA0T2T2*
     .             ((A0*phiA0T2T2)/deltA0T2T2+
     .             tmp107/deltA0T2T2+
     .             (phiA0T2T2*tmp182)/deltA0T2T2+
     .             tmp99/deltA0T2T2))/T2**2+
     .         0.5d0*((A0-2d0*T2)*tmp18)-
     .         0.5d0*(-8d0*A0*
     .            ((A0*phiA0T2T2)/deltA0T2T2+
     .              tmp107/deltA0T2T2+
     .              (phiA0T2T2*tmp182)/deltA0T2T2+
     .              tmp99/deltA0T2T2)+
     .             deltA0T2T2*
     .            ((4d0*A0**2*phiA0T2T2)/deltA0T2T2**2+
     .              (phiA0T2T2*
     .                 (-1d0-(A0-T2)**2/T2**2-
     .                   (2d0*(A0-T2))/T2))/deltA0T2T2+
     .              (1d0+(A0-T2)/T2)/deltA0T2T2+
     .              (1d0-(-A0+T2)/T2)/deltA0T2T2+
     .              (4d0*A0*tmp107)/deltA0T2T2**2+
     .              (4d0*A0*phiA0T2T2*tmp182)/deltA0T2T2**2+
     .              (4d0*A0*tmp99)/deltA0T2T2**2+
     .              (A0*
     .                 ((A0*phiA0T2T2)/deltA0T2T2+
     .                   tmp107/deltA0T2T2+
     .                   (phiA0T2T2*tmp182)/deltA0T2T2+
     .                   tmp99/deltA0T2T2))/deltA0T2T2+
     .              (tmp182*
     .                 ((A0*phiA0T2T2)/deltA0T2T2+
     .                   tmp107/deltA0T2T2+
     .                   (phiA0T2T2*tmp182)/deltA0T2T2+
     .                   tmp99/deltA0T2T2))/deltA0T2T2))/T2)))
      tmp245 = (2d0*(-1d0+Logmu2)*mu2)/t-((-1d0+LogT1)*T1)/t-
     .   ((-1d0+LogT2)*T2)/t+2d0*tmp228+2d0*tmp243+
     .   sbe**2*(8d0*(-1d0+Logt)+2d0*(-1d0+Logt)**2+
     .      t*(4d0/t+(2d0*(2d0-2d0*Logt))/t-(2d0*Logt)/t-
     .       (4d0*t-2d0*Logt*t)/t**2)+t**2*tmp15+2d0*tmp77)-
     .   (-mu2-t+T2)*(2d0/t-Logmu2/t-LogT2/t-
     .      0.5d0*(LogT2*(mu2-t-T2))/t**2-
     .      0.5d0*(Logmu2*(-mu2-t+T2))/t**2-
     .      0.5d0*((2d0*mu2*phimu2T2T)/T2+
     .        4d0*(-mu2+t-T2)*
     .         ((mu2*phimu2T2T*(mu2-t+T2))/(deltmu2T2T*T2)+
     .           (mu2*tmp111)/(deltmu2T2T*t))+
     .        deltmu2T2T*
     .         (((2d0-Logmu2+2d0*Logt-LogT2)*mu2)/
     .            (deltmu2T2T*t)-
     .           (mu2*phimu2T2T)/(deltmu2T2T*T2)-
     .           (2d0*mu2*phimu2T2T*(-mu2+t-T2)*
     .            (mu2-t+T2))/(deltmu2T2T**2*T2)-
     .           (mu2*tmp111)/(deltmu2T2T*t**2)-
     .           (2d0*mu2*(-mu2+t-T2)*tmp111)/
     .            (deltmu2T2T**2*t)+
     .           ((mu2-t+T2)*
     .            ((mu2*phimu2T2T*(mu2-t+T2))/
     .               (deltmu2T2T*T2)+
     .              (mu2*tmp111)/(deltmu2T2T*t)))/deltmu2T2T))/
     .      mu2)+cbe**2*
     .    (8d0*(-1d0+Logt)+2d0*(-1d0+Logt)**2-
     .      (2d0*A0*(-1d0+LogA0))/t+t**2*tmp15+
     .      4d0*(-1d0+4d0*Logt-Logt**2-(A0*LogA0)/t+
     .       (2d0*A0*phiA0tt)/t+(Logt*(A0-2d0*t))/t+
     .       0.5d0*(deltA0tt*phiA0tt)/t**2-
     .       0.5d0*(deltA0tt*
     .           ((A0*phiA0tt)/deltA0tt+
     .             (phiA0tt*tmp175)/deltA0tt+tmp65/deltA0tt+
     .             tmp70/deltA0tt))/t)-
     .      (A0-2d0*t)*(-((deltA0tt*phiA0tt)/t**3)+
     .       (A0*LogA0)/t**2+4d0/t-(4d0*Logt)/t+
     .       (-4d0*A0*phiA0tt+
     .          deltA0tt*
     .           ((A0*phiA0tt)/deltA0tt+
     .             (phiA0tt*tmp175)/deltA0tt+tmp65/deltA0tt+
     .             tmp70/deltA0tt))/t**2+
     .       0.5d0*((A0-2d0*t)*tmp16)-
     .       0.5d0*(-8d0*A0*((A0*phiA0tt)/deltA0tt+
     .            (phiA0tt*tmp175)/deltA0tt+tmp65/deltA0tt+
     .            tmp70/deltA0tt)+
     .           deltA0tt*
     .            ((4d0*A0**2*phiA0tt)/deltA0tt**2+
     .            (phiA0tt*
     .               (-1d0-(A0-t)**2/t**2-(2d0*(A0-t))/t))/
     .             deltA0tt+(1d0+(A0-t)/t)/deltA0tt+
     .            (1d0-(-A0+t)/t)/deltA0tt+
     .            (4d0*A0*phiA0tt*tmp175)/deltA0tt**2+
     .            (4d0*A0*tmp65)/deltA0tt**2+
     .            (4d0*A0*tmp70)/deltA0tt**2+
     .            (A0*((A0*phiA0tt)/deltA0tt+
     .                 (phiA0tt*tmp175)/deltA0tt+
     .                 tmp65/deltA0tt+tmp70/deltA0tt))/
     .             deltA0tt+
     .            (tmp175*
     .               ((A0*phiA0tt)/deltA0tt+
     .                 (phiA0tt*tmp175)/deltA0tt+
     .                 tmp65/deltA0tt+tmp70/deltA0tt))/
     .             deltA0tt))/t))
      tmp245 = tmp245-
     .   (-mu2-t+T1)*(2d0/t-Logmu2/t-LogT1/t-
     .      0.5d0*(LogT1*(mu2-t-T1))/t**2-
     .      0.5d0*(Logmu2*(-mu2-t+T1))/t**2-
     .      0.5d0*((2d0*mu2*phimu2tT1)/T1+
     .        4d0*(-mu2+t-T1)*
     .         ((mu2*phimu2tT1*(mu2-t+T1))/(deltmu2tT1*T1)+
     .           (mu2*tmp91)/(deltmu2tT1*t))+
     .        deltmu2tT1*
     .         (((2d0-Logmu2+2d0*Logt-LogT1)*mu2)/
     .            (deltmu2tT1*t)-
     .           (mu2*phimu2tT1)/(deltmu2tT1*T1)-
     .           (2d0*mu2*phimu2tT1*(-mu2+t-T1)*
     .            (mu2-t+T1))/(deltmu2tT1**2*T1)-
     .           (mu2*tmp91)/(deltmu2tT1*t**2)-
     .           (2d0*mu2*(-mu2+t-T1)*tmp91)/
     .            (deltmu2tT1**2*t)+
     .           ((mu2-t+T1)*
     .            ((mu2*phimu2tT1*(mu2-t+T1))/
     .               (deltmu2tT1*T1)+
     .              (mu2*tmp91)/(deltmu2tT1*t)))/deltmu2tT1))/
     .      mu2)
      Dtt = ht**2*tmp245+
     .   tmp10*(-5d0-2d0*Li2bt+4d0*Logt-Logt**2+
     .      (b*(-1d0+Logb))/t+(2d0*Logt*(b-t))/t+
     .      (2d0*b*(Logb-Logt)*(-b+t))/((1d0-b/t)*t**2)+
     .      (-2d0*b*Logb+4d0*t)/t+
     .      0.5d0*((b+t)*((4d0*b*(Logb-Logt))/((1d0-b/t)*t**2)+
     .         8d0/t-(4d0*Logt)/t-(-2d0*b*Logb+4d0*t)/t**2+
     .         (b-t)*tmp16-2d0*(-b+t)*tmp207)))+
     .   tmp61*(-((B1*(-1d0+LogB1))/t)+((-1d0+Logmu2)*mu2)/t+
     .      2d0*tmp231-(B1-mu2-t)*
     .       (2d0/t-LogB1/t-Logmu2/t-
     .       0.5d0*(Logmu2*(B1-mu2-t))/t**2-
     .       0.5d0*(LogB1*(-B1+mu2-t))/t**2-
     .       0.5d0*(2d0*phiB1tmu2+deltB1tmu2*tmp232+
     .           4d0*(-B1-mu2+t)*
     .            ((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .            (mu2*tmp67)/(deltB1tmu2*t)))/mu2))+
     .   tmp60*(-((B2*(-1d0+LogB2))/t)+((-1d0+Logmu2)*mu2)/t+
     .      2d0*tmp234-(B2-mu2-t)*
     .       (2d0/t-LogB2/t-Logmu2/t-
     .       0.5d0*(Logmu2*(B2-mu2-t))/t**2-
     .       0.5d0*(LogB2*(-B2+mu2-t))/t**2-
     .       0.5d0*(2d0*phiB2tmu2+deltB2tmu2*tmp236+
     .           4d0*(-B2-mu2+t)*
     .            ((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .            (mu2*tmp68)/(deltB2tmu2*t)))/mu2))-
     .   2d0*hb*ht*mu*s2b*(tmp235/mt+
     .      mt*(-(LogB1/t)+LogB2/t-
     .       0.5d0*(Logmu2*(B1-mu2-t))/t**2+
     .       0.5d0*(Logmu2*(B2-mu2-t))/t**2-
     .       0.5d0*(LogB1*(-B1+mu2-t))/t**2+
     .       0.5d0*(LogB2*(-B2+mu2-t))/t**2-
     .       0.5d0*(2d0*phiB1tmu2+deltB1tmu2*tmp232+
     .           4d0*(-B1-mu2+t)*
     .            ((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .            (mu2*tmp67)/(deltB1tmu2*t)))/mu2+
     .       0.5d0*(2d0*phiB2tmu2+deltB2tmu2*tmp236+
     .           4d0*(-B2-mu2+t)*
     .            ((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .            (mu2*tmp68)/(deltB2tmu2*t)))/mu2)-
     .      0.25d0*tmp75/t**1.5d0)
      Dtt = Dtt+tmp9*
     .    (-((A0*(-1d0+LogA0))/t)+(b*(-1d0+Logb))/t+2d0*tmp229-
     .      (A0-b-t)*(2d0/t-LogA0/t-Logb/t-
     .       0.5d0*(Logb*(A0-b-t))/t**2-
     .       0.5d0*(LogA0*(-A0+b-t))/t**2-
     .       0.5d0*((2d0*b*phiA0bt)/t+deltA0bt*tmp230+
     .           4d0*(-A0-b+t)*
     .            ((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .            (b*tmp66)/(deltA0bt*t)))/b))-
     .   4d0*cbe*hb*ht*mb*sbe*
     .    ((0.5d0-2d0*Logt+(phiA0bt*(-A0-b+t))/t-
     .       0.5d0*(LogA0*Logb)+0.5d0*(LogA0*Logt)+
     .       0.5d0*(Logb*Logt)-0.5d0*(Logb*(A0-b-t))/t-
     .       0.5d0*(LogA0*(-A0+b-t))/t+0.5d0*tmp210+
     .       0.5d0*(deltA0bt*
     .           ((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .             (b*tmp66)/(deltA0bt*t)))/b)/mt+
     .      mt*(-2d0/t+LogA0/t+Logb/t+
     .       0.5d0*(Logb*(A0-b-t))/t**2+
     .       0.5d0*(LogA0*(-A0+b-t))/t**2+
     .       0.5d0*((4d0*b*(Logb-Logt))/((1d0-b/t)*t**2)+8d0/t-
     .          (4d0*Logt)/t-(-2d0*b*Logb+4d0*t)/t**2+
     .          (b-t)*tmp16-2d0*(-b+t)*tmp207)+
     .       0.5d0*((2d0*b*phiA0bt)/t+deltA0bt*tmp230+
     .           4d0*(-A0-b+t)*
     .            ((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .            (b*tmp66)/(deltA0bt*t)))/b)-
     .      0.25d0*(-2d0*A0*LogA0-2d0*b*Logb-2d0*Logt*t-
     .        0.5d0*(Logb*Logt*(A0-b-t))-
     .        0.5d0*(LogA0*Logt*(-A0+b-t))+
     .        0.5d0*(deltA0bt*phiA0bt)/t-
     .        0.5d0*(LogA0*Logb*(-A0-b+t))+
     .        2.5d0*(A0+b+t)+0.5d0*tmp76)/t**1.5d0)+
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (cbe*hb*ht*sbe*
     .       (0.25d0*(s2b*tmp130)/t**1.5d0
     .        -0.5d0*(mb*tmp56)/t**1.5d0)-
     .        cbe**2*hb**2*
     .       (0.125d0*(mb*s2b*s2t)/t**1.5d0-
     .         0.125d0*((1d0+c2b)*s2t*Xb)/t**1.5d0)-
     .      ht**2*sbe**2*
     .       (0.125d0*(mb*s2b*s2t)/t**1.5d0-
     .         0.125d0*((1d0-c2b)*s2t*Xt)/t**1.5d0)))
      Dtt = Dtt+tmp89*
     .    (-(cbe*hb*ht*sbe*
     .       (0.25d0*(s2b*tmp156)/t**1.5d0
     .        -0.5d0*(mb*tmp56)/t**1.5d0))-
     .      hb**2*sbe**2*
     .       (0.125d0*(mb*s2b*s2t)/t**1.5d0-
     .       0.125d0*((1d0+c2b)*s2t*Yb)/t**1.5d0)-
     .      cbe**2*ht**2*(0.125d0*(mb*s2b*s2t)/t**1.5d0-
     .       0.125d0*((1d0-c2b)*s2t*Yt)/t**1.5d0))+
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (cbe*hb*ht*sbe*
     .       (0.25d0*(s2b*tmp131)/t**1.5d0
     .        -0.5d0*(mb*tmp59)/t**1.5d0)-
     .        cbe**2*hb**2*
     .       (-(0.125d0*(mb*s2b*s2t)/t**1.5d0)+
     .         0.125d0*((1d0+c2b)*s2t*Xb)/t**1.5d0)-
     .      ht**2*sbe**2*
     .       (-(0.125d0*(mb*s2b*s2t)/t**1.5d0)+
     .         0.125d0*((1d0-c2b)*s2t*Xt)/t**1.5d0)))+
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (cbe*hb*ht*sbe*
     .       (-(0.25d0*(s2b*tmp130)/t**1.5d0)-
     .         0.5d0*(mb*tmp59)/t**1.5d0)-
     .      cbe**2*hb**2*
     .       (-(0.125d0*(mb*s2b*s2t)/t**1.5d0)-
     .         0.125d0*((1d0-c2b)*s2t*Xb)/t**1.5d0)-
     .      ht**2*sbe**2*
     .       (-(0.125d0*(mb*s2b*s2t)/t**1.5d0)-
     .         0.125d0*((1d0+c2b)*s2t*Xt)/t**1.5d0)))+
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (cbe*hb*ht*sbe*
     .       (-(0.25d0*(s2b*tmp131)/t**1.5d0)-
     .         0.5d0*(mb*tmp56)/t**1.5d0)-
     .      cbe**2*hb**2*
     .       (0.125d0*(mb*s2b*s2t)/t**1.5d0+
     .         0.125d0*((1d0-c2b)*s2t*Xb)/t**1.5d0)-
     .      ht**2*sbe**2*
     .       (0.125d0*(mb*s2b*s2t)/t**1.5d0+
     .         0.125d0*((1d0+c2b)*s2t*Xt)/t**1.5d0)))-
     .   0.25d0*(ht**2*((cbe**2*s2t*tmp114*Yt)/t**1.5d0-
     .      (cbe**2*s2t*tmp93*Yt)/t**1.5d0+
     .      0.5d0*(s2t*sbe**2*tmp115*Xt)/t**1.5d0-
     .      0.5d0*(s2t*sbe**2*tmp94*Xt)/t**1.5d0))
      Dtt = Dtt+tmp109*
     .    (-(cbe*hb*ht*sbe*
     .       (0.25d0*(s2b*tmp157)/t**1.5d0
     .        -0.5d0*(mb*tmp59)/t**1.5d0))-
     .      hb**2*sbe**2*
     .       (-(0.125d0*(mb*s2b*s2t)/t**1.5d0)+
     .       0.125d0*((1d0+c2b)*s2t*Yb)/t**1.5d0)-
     .      cbe**2*ht**2*(-(0.125d0*(mb*s2b*s2t)/t**1.5d0)+
     .       0.125d0*((1d0-c2b)*s2t*Yt)/t**1.5d0))+
     .   tmp88*(-(cbe*hb*ht*sbe*
     .       (-(0.25d0*(s2b*tmp156)/t**1.5d0)-
     .         0.5d0*(mb*tmp59)/t**1.5d0))-
     .      hb**2*sbe**2*(-(0.125d0*(mb*s2b*s2t)/t**1.5d0)-
     .       0.125d0*((1d0-c2b)*s2t*Yb)/t**1.5d0)-
     .      cbe**2*ht**2*(-(0.125d0*(mb*s2b*s2t)/t**1.5d0)-
     .       0.125d0*((1d0+c2b)*s2t*Yt)/t**1.5d0))+
     .   tmp108*(-(cbe*hb*ht*sbe*
     .       (-(0.25d0*(s2b*tmp157)/t**1.5d0)-
     .         0.5d0*(mb*tmp56)/t**1.5d0))-
     .      hb**2*sbe**2*(0.125d0*(mb*s2b*s2t)/t**1.5d0+
     .       0.125d0*((1d0-c2b)*s2t*Yb)/t**1.5d0)-
     .      cbe**2*ht**2*(0.125d0*(mb*s2b*s2t)/t**1.5d0+
     .       0.125d0*((1d0+c2b)*s2t*Yt)/t**1.5d0))
      Dc2tc2t = (b*(-1d0+Logb)*(-1d0+Logmu2)*mu2-
     .      b*(-1d0+Logb)*(-1d0+LogT2)*T2-
     .      (-1d0+Logmu2)*(-1d0+LogT2)*mu2*T2-
     .      (-b-mu2+T2)*tmp112)*
     .    (0.125d0*hb**2/c2t**3-0.125d0*ht**2/c2t**3)+
     .   (b*(-1d0+Logb)*(-1d0+Logmu2)*mu2-
     .      b*(-1d0+Logb)*(-1d0+LogT1)*T1-
     .      (-1d0+Logmu2)*(-1d0+LogT1)*mu2*T1-
     .      (-b-mu2+T1)*tmp92)*
     .    (-(0.125d0*hb**2/c2t**3)+0.125d0*ht**2/c2t**3)+
     .   ht**2*(0.0625d0*(B1*(1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT1)*T1)/
     .      c2t**3+0.0625d0*
     .       (B2*(1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT1)*T1)/c2t**3+
     .      cbe**2*(0.125d0*(A0*(-1d0+LogA0)*(-1d0+LogT1)*T1)/
     .         c2t**3-0.125d0*
     .        (A0*(-1d0+LogA0)*(-1d0+LogT2)*T2)/c2t**3)-
     .      0.0625d0*(B1*(1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT2)*T2)/
     .      c2t**3-0.0625d0*
     .       (B2*(1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT2)*T2)/c2t**3)+
     .   hb**2*(-(0.0625d0*(B1*(1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT1)*
     .          T1)/c2t**3)-
     .      0.0625d0*(B2*(1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT1)*T1)/
     .      c2t**3+sbe**2*
     .       (-(0.125d0*(A0*(-1d0+LogA0)*(-1d0+LogT1)*T1)/c2t**3)+
     .       0.125d0*(A0*(-1d0+LogA0)*(-1d0+LogT2)*T2)/c2t**3)+
     .      0.0625d0*(B1*(1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT2)*T2)/
     .      c2t**3+0.0625d0*
     .       (B2*(1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT2)*T2)/c2t**3)+
     .   0.5d0*(hb*ht*mb*mu*tmp113)/s2t**3+
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (-(cbe**2*hb**2*
     .         (0.0625d0*(b*(1d0-c2b))/c2t**3-
     .           0.125d0*(mb*mt*s2b)/s2t**3-
     .           0.0625d0*((1d0+c2b)*t)/c2t**3-
     .           0.125d0*(mb*s2b*Xb)/c2t**3+
     .           0.125d0*((1d0+c2b)*mt*Xb)/s2t**3+
     .           0.0625d0*((1d0+c2b)*Xb**2)/c2t**3))+
     .      cbe*hb*ht*sbe*
     .       (-(mt*s2b*tmp25)+2d0*mb*mt*tmp52-
     .         0.125d0*(s2b*(b+t))/s2t**3+
     .         0.25d0*(mb*tmp128)/s2t**3-0.125d0*(s2b*Xb*Xt)/s2t**3
     .         )-ht**2*sbe**2*
     .       (-(0.0625d0*(b*(1d0+c2b))/c2t**3)-
     .         0.125d0*(mb*mt*s2b)/s2t**3+
     .         0.0625d0*((1d0-c2b)*t)/c2t**3+
     .         0.125d0*(mb*s2b*Xt)/c2t**3+
     .         0.125d0*((1d0-c2b)*mt*Xt)/s2t**3-
     .         0.0625d0*((1d0-c2b)*Xt**2)/c2t**3)))
      Dc2tc2t = Dc2tc2t+
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (-(cbe**2*hb**2*
     .         (-(0.0625d0*(b*(1d0-c2b))/c2t**3)+
     .           0.125d0*(mb*mt*s2b)/s2t**3+
     .           0.0625d0*((1d0+c2b)*t)/c2t**3+
     .           0.125d0*(mb*s2b*Xb)/c2t**3-
     .           0.125d0*((1d0+c2b)*mt*Xb)/s2t**3-
     .           0.0625d0*((1d0+c2b)*Xb**2)/c2t**3))+
     .      cbe*hb*ht*sbe*
     .       (-(mt*s2b*tmp26)+2d0*mb*mt*tmp53+
     .         0.125d0*(s2b*(b+t))/s2t**3-
     .         0.25d0*(mb*tmp128)/s2t**3+0.125d0*(s2b*Xb*Xt)/s2t**3
     .         )-ht**2*sbe**2*
     .       (0.0625d0*(b*(1d0+c2b))/c2t**3+
     .         0.125d0*(mb*mt*s2b)/s2t**3-
     .         0.0625d0*((1d0-c2b)*t)/c2t**3-
     .         0.125d0*(mb*s2b*Xt)/c2t**3-
     .         0.125d0*((1d0-c2b)*mt*Xt)/s2t**3+
     .         0.0625d0*((1d0-c2b)*Xt**2)/c2t**3)))+
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (-(cbe**2*hb**2*
     .         (0.0625d0*(b*(1d0+c2b))/c2t**3+
     .           0.125d0*(mb*mt*s2b)/s2t**3-
     .           0.0625d0*((1d0-c2b)*t)/c2t**3+
     .           0.125d0*(mb*s2b*Xb)/c2t**3+
     .           0.125d0*((1d0-c2b)*mt*Xb)/s2t**3+
     .           0.0625d0*((1d0-c2b)*Xb**2)/c2t**3))+
     .      cbe*hb*ht*sbe*
     .       (mt*s2b*tmp25+2d0*mb*mt*tmp53+
     .         0.125d0*(s2b*(b+t))/s2t**3+
     .         0.25d0*(mb*tmp129)/s2t**3+0.125d0*(s2b*Xb*Xt)/s2t**3
     .         )-ht**2*sbe**2*
     .       (-(0.0625d0*(b*(1d0-c2b))/c2t**3)+
     .         0.125d0*(mb*mt*s2b)/s2t**3+
     .         0.0625d0*((1d0+c2b)*t)/c2t**3-
     .         0.125d0*(mb*s2b*Xt)/c2t**3+
     .         0.125d0*((1d0+c2b)*mt*Xt)/s2t**3-
     .         0.0625d0*((1d0+c2b)*Xt**2)/c2t**3)))
      Dc2tc2t = Dc2tc2t+
     .   tmp109*(-(hb**2*sbe**2*
     .       (0.0625d0*(b*(1d0-c2b))/c2t**3-
     .         0.125d0*(mb*mt*s2b)/s2t**3-
     .         0.0625d0*((1d0+c2b)*t)/c2t**3-
     .         0.125d0*(mb*s2b*Yb)/c2t**3+
     .         0.125d0*((1d0+c2b)*mt*Yb)/s2t**3+
     .         0.0625d0*((1d0+c2b)*Yb**2)/c2t**3))-
     .      cbe*hb*ht*sbe*(-(mt*s2b*tmp41)+2d0*mb*mt*tmp52-
     .       0.125d0*(s2b*(b+t))/s2t**3+
     .       0.25d0*(mb*tmp154)/s2t**3-0.125d0*(s2b*Yb*Yt)/s2t**3)-
     .      cbe**2*ht**2*
     .       (-(0.0625d0*(b*(1d0+c2b))/c2t**3)-
     .       0.125d0*(mb*mt*s2b)/s2t**3+
     .       0.0625d0*((1d0-c2b)*t)/c2t**3+
     .       0.125d0*(mb*s2b*Yt)/c2t**3+
     .       0.125d0*((1d0-c2b)*mt*Yt)/s2t**3-
     .       0.0625d0*((1d0-c2b)*Yt**2)/c2t**3))-
     .   0.25d0*(ht**2*((cbe**2*mt*tmp114*Yt)/s2t**3-
     .      (cbe**2*mt*tmp93*Yt)/s2t**3+
     .      0.5d0*(mt*sbe**2*tmp115*Xt)/s2t**3-
     .      0.5d0*(mt*sbe**2*tmp94*Xt)/s2t**3))+
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (-(cbe**2*hb**2*
     .         (-(0.0625d0*(b*(1d0+c2b))/c2t**3)-
     .           0.125d0*(mb*mt*s2b)/s2t**3+
     .           0.0625d0*((1d0-c2b)*t)/c2t**3-
     .           0.125d0*(mb*s2b*Xb)/c2t**3-
     .           0.125d0*((1d0-c2b)*mt*Xb)/s2t**3-
     .           0.0625d0*((1d0-c2b)*Xb**2)/c2t**3))+
     .      cbe*hb*ht*sbe*
     .       (mt*s2b*tmp26+2d0*mb*mt*tmp52-
     .         0.125d0*(s2b*(b+t))/s2t**3-
     .         0.25d0*(mb*tmp129)/s2t**3-0.125d0*(s2b*Xb*Xt)/s2t**3
     .         )-ht**2*sbe**2*
     .       (0.0625d0*(b*(1d0-c2b))/c2t**3-
     .         0.125d0*(mb*mt*s2b)/s2t**3-
     .         0.0625d0*((1d0+c2b)*t)/c2t**3+
     .         0.125d0*(mb*s2b*Xt)/c2t**3-
     .         0.125d0*((1d0+c2b)*mt*Xt)/s2t**3+
     .         0.0625d0*((1d0+c2b)*Xt**2)/c2t**3)))
      Dc2tc2t = Dc2tc2t+
     .   tmp89*(-(hb**2*sbe**2*
     .       (-(0.0625d0*(b*(1d0-c2b))/c2t**3)+
     .         0.125d0*(mb*mt*s2b)/s2t**3+
     .         0.0625d0*((1d0+c2b)*t)/c2t**3+
     .         0.125d0*(mb*s2b*Yb)/c2t**3-
     .         0.125d0*((1d0+c2b)*mt*Yb)/s2t**3-
     .         0.0625d0*((1d0+c2b)*Yb**2)/c2t**3))-
     .      cbe*hb*ht*sbe*(-(mt*s2b*tmp42)+2d0*mb*mt*tmp53+
     .       0.125d0*(s2b*(b+t))/s2t**3-
     .       0.25d0*(mb*tmp154)/s2t**3+0.125d0*(s2b*Yb*Yt)/s2t**3)-
     .      cbe**2*ht**2*
     .       (0.0625d0*(b*(1d0+c2b))/c2t**3+
     .       0.125d0*(mb*mt*s2b)/s2t**3-
     .       0.0625d0*((1d0-c2b)*t)/c2t**3-
     .       0.125d0*(mb*s2b*Yt)/c2t**3-
     .       0.125d0*((1d0-c2b)*mt*Yt)/s2t**3+
     .       0.0625d0*((1d0-c2b)*Yt**2)/c2t**3))+
     .   tmp108*(-(hb**2*sbe**2*
     .       (0.0625d0*(b*(1d0+c2b))/c2t**3+
     .         0.125d0*(mb*mt*s2b)/s2t**3-
     .         0.0625d0*((1d0-c2b)*t)/c2t**3+
     .         0.125d0*(mb*s2b*Yb)/c2t**3+
     .         0.125d0*((1d0-c2b)*mt*Yb)/s2t**3+
     .         0.0625d0*((1d0-c2b)*Yb**2)/c2t**3))-
     .      cbe*hb*ht*sbe*(mt*s2b*tmp41+2d0*mb*mt*tmp53+
     .       0.125d0*(s2b*(b+t))/s2t**3+
     .       0.25d0*(mb*tmp155)/s2t**3+0.125d0*(s2b*Yb*Yt)/s2t**3)-
     .      cbe**2*ht**2*
     .       (-(0.0625d0*(b*(1d0-c2b))/c2t**3)+
     .       0.125d0*(mb*mt*s2b)/s2t**3+
     .       0.0625d0*((1d0+c2b)*t)/c2t**3-
     .       0.125d0*(mb*s2b*Yt)/c2t**3+
     .       0.125d0*((1d0+c2b)*mt*Yt)/s2t**3-
     .       0.0625d0*((1d0+c2b)*Yt**2)/c2t**3))+
     .   tmp88*(-(hb**2*sbe**2*
     .       (-(0.0625d0*(b*(1d0+c2b))/c2t**3)-
     .         0.125d0*(mb*mt*s2b)/s2t**3+
     .         0.0625d0*((1d0-c2b)*t)/c2t**3-
     .         0.125d0*(mb*s2b*Yb)/c2t**3-
     .         0.125d0*((1d0-c2b)*mt*Yb)/s2t**3-
     .         0.0625d0*((1d0-c2b)*Yb**2)/c2t**3))-
     .      cbe*hb*ht*sbe*(mt*s2b*tmp42+2d0*mb*mt*tmp52-
     .       0.125d0*(s2b*(b+t))/s2t**3-
     .       0.25d0*(mb*tmp155)/s2t**3-0.125d0*(s2b*Yb*Yt)/s2t**3)-
     .      cbe**2*ht**2*
     .       (0.0625d0*(b*(1d0-c2b))/c2t**3-
     .       0.125d0*(mb*mt*s2b)/s2t**3-
     .       0.0625d0*((1d0+c2b)*t)/c2t**3+
     .       0.125d0*(mb*s2b*Yt)/c2t**3-
     .       0.125d0*((1d0+c2b)*mt*Yt)/s2t**3+
     .       0.0625d0*((1d0+c2b)*Yt**2)/c2t**3))
      DT1t = -(tmp192*tmp238)-tmp191*tmp239-
     .   0.5d0*(tmp184*tmp213)-0.5d0*(tmp183*tmp214)+
     .   ht**2*(1d0-3d0*Logt+Logmu2*Logt-
     .      (-1d0+Logt)*(-1d0+LogT1)+LogT1-Logmu2*LogT1+
     .      (phimu2tT1*(-mu2+t-T1))/T1-
     .      (phimu2tT1*(-mu2-t+T1))/T1-
     .      0.5d0*(LogT1*(mu2-t-T1))/t+
     .      0.5d0*(Logt*(mu2-t-T1))/T1+
     .      0.5d0*(Logmu2*(-mu2+t-T1))/T1-
     .      0.5d0*(Logmu2*(-mu2-t+T1))/t-
     .      0.5d0*(deltmu2tT1*
     .        ((mu2*phimu2tT1*(mu2+t-T1))/(deltmu2tT1*T1)+
     .          (mu2*tmp84)/(deltmu2tT1*T1)))/mu2+
     .      0.5d0*(deltmu2tT1*
     .        ((mu2*phimu2tT1*(mu2-t+T1))/(deltmu2tT1*T1)+
     .          (mu2*tmp91)/(deltmu2tT1*t)))/mu2-
     .      (-mu2-t+T1)*
     .       (phimu2tT1/T1-
     .       ((-mu2+t-T1)*
     .          ((mu2*phimu2tT1*(mu2+t-T1))/
     .             (deltmu2tT1*T1)+(mu2*tmp84)/(deltmu2tT1*T1))
     .          )/mu2-((-mu2-t+T1)*
     .          ((mu2*phimu2tT1*(mu2-t+T1))/
     .             (deltmu2tT1*T1)+(mu2*tmp91)/(deltmu2tT1*t)))
     .         /mu2+0.5d0*Logmu2/t-0.5d0*LogT1/t+
     .       0.5d0*Logmu2/T1-0.5d0*Logt/T1+
     .       0.5d0*(mu2-t-T1)/(t*T1)-
     .       0.5d0*(deltmu2tT1*
     .           ((mu2*phimu2tT1)/(deltmu2tT1*T1)-
     .             (2d0*mu2*phimu2tT1*(-mu2+t-T1)*
     .              (mu2+t-T1))/(deltmu2tT1**2*T1)+
     .             (mu2*(Logmu2-Logt-(-mu2+t)/t-T1/t))/
     .            (deltmu2tT1*T1)-
     .             (2d0*mu2*(-mu2+t-T1)*tmp84)/
     .            (deltmu2tT1**2*T1)+
     .             ((mu2+t-T1)*
     .              ((mu2*phimu2tT1*(mu2-t+T1))/
     .                 (deltmu2tT1*T1)+
     .                (mu2*tmp91)/(deltmu2tT1*t)))/deltmu2tT1))
     .          /mu2))-
     .   0.25d0*(ht**2*(cbe**2*(4d0+(2d0*s2t*Yt)/mt)*
     .       (-1d0+4d0*LogT1-LogT1**2-(A0*LogA0)/T1+
     .         (2d0*A0*phiA0T1T1)/T1+(LogT1*(A0-2d0*T1))/T1+
     .         0.5d0*(deltA0T1T1*phiA0T1T1)/T1**2-
     .         0.5d0*(deltA0T1T1*
     .             ((A0*phiA0T1T1)/deltA0T1T1+
     .             (phiA0T1T1*tmp178)/deltA0T1T1+
     .             tmp80/deltA0T1T1+tmp87/deltA0T1T1))/T1)+
     .      0.5d0*(sbe**2*tmp95*(4d0+(2d0*s2t*Xt)/mt))))
      DT2t = -(tmp196*tmp222)-tmp195*tmp223+
     .   ht**2*(1d0-3d0*Logt+Logmu2*Logt-
     .      (-1d0+Logt)*(-1d0+LogT2)+LogT2-Logmu2*LogT2+
     .      (phimu2T2T*(-mu2+t-T2))/T2-
     .      (phimu2T2T*(-mu2-t+T2))/T2-
     .      0.5d0*(LogT2*(mu2-t-T2))/t+
     .      0.5d0*(Logt*(mu2-t-T2))/T2+
     .      0.5d0*(Logmu2*(-mu2+t-T2))/T2-
     .      0.5d0*(Logmu2*(-mu2-t+T2))/t-
     .      0.5d0*(deltmu2T2T*
     .        ((mu2*phimu2T2T*(mu2+t-T2))/(deltmu2T2T*T2)+
     .          (mu2*tmp103)/(deltmu2T2T*T2)))/mu2+
     .      0.5d0*(deltmu2T2T*
     .        ((mu2*phimu2T2T*(mu2-t+T2))/(deltmu2T2T*T2)+
     .          (mu2*tmp111)/(deltmu2T2T*t)))/mu2-
     .      (-mu2-t+T2)*
     .       (phimu2T2T/T2-
     .       ((-mu2+t-T2)*
     .          ((mu2*phimu2T2T*(mu2+t-T2))/
     .             (deltmu2T2T*T2)+(mu2*tmp103)/(deltmu2T2T*T2)
     .            ))/mu2-
     .       ((-mu2-t+T2)*
     .          ((mu2*phimu2T2T*(mu2-t+T2))/
     .             (deltmu2T2T*T2)+(mu2*tmp111)/(deltmu2T2T*t))
     .          )/mu2+0.5d0*Logmu2/t-0.5d0*LogT2/t+
     .       0.5d0*Logmu2/T2-0.5d0*Logt/T2+
     .       0.5d0*(mu2-t-T2)/(t*T2)-
     .       0.5d0*(deltmu2T2T*
     .           ((mu2*phimu2T2T)/(deltmu2T2T*T2)-
     .             (2d0*mu2*phimu2T2T*(-mu2+t-T2)*
     .              (mu2+t-T2))/(deltmu2T2T**2*T2)+
     .             (mu2*(Logmu2-Logt-(-mu2+t)/t-T2/t))/
     .            (deltmu2T2T*T2)-
     .             (2d0*mu2*(-mu2+t-T2)*tmp103)/
     .            (deltmu2T2T**2*T2)+
     .             ((mu2+t-T2)*
     .              ((mu2*phimu2T2T*(mu2-t+T2))/
     .                 (deltmu2T2T*T2)+
     .                (mu2*tmp111)/(deltmu2T2T*t)))/deltmu2T2T)
     .           )/mu2))-0.5d0*(tmp188*tmp217)-
     .   0.5d0*(tmp187*tmp218)-
     .   0.25d0*(ht**2*(cbe**2*(4d0-(2d0*s2t*Yt)/mt)*
     .       (-1d0+4d0*LogT2-LogT2**2-(A0*LogA0)/T2+
     .         (2d0*A0*phiA0T2T2)/T2+(LogT2*(A0-2d0*T2))/T2+
     .         0.5d0*(deltA0T2T2*phiA0T2T2)/T2**2-
     .         0.5d0*(deltA0T2T2*
     .             ((A0*phiA0T2T2)/deltA0T2T2+
     .             tmp107/deltA0T2T2+
     .             (phiA0T2T2*tmp182)/deltA0T2T2+
     .             tmp99/deltA0T2T2))/T2)+
     .      0.5d0*(sbe**2*tmp116*(4d0-(2d0*s2t*Xt)/mt))))
      DT1T2 = ht**2*(c2t**2+(-1d0+LogT1)*tmp8+
     .      (-1d0+LogT2)*tmp8+(-1d0+LogT1)*(-1d0+LogT2)*tmp8-
     .      0.5d0*((-1d0+Nc)*s2t**2))-
     .   0.25d0*(ht**2*((1d0+c2t**2)*sbe**2*
     .       ((2d0*LogT1)/T1+(-2d0-2d0*LogT2)/T1+
     .         (2d0*(-LogT1+LogT2)*(T1-T2)*T2)/
     .          (T1**3*(1d0-T2/T1)**2)-
     .         (2d0*(-LogT1+LogT2))/(T1*(1d0-T2/T1))+
     .         (2d0*(T1-T2))/(T1**2*(1d0-T2/T1))+
     .         (2d0*(-LogT1+LogT2)*(T1-T2))/
     .          (T1**2*(1d0-T2/T1))-
     .         (2d0*(-LogT1+LogT2)*T2)/(T1**2*(1d0-T2/T1)))*Xt**2
     .       +2d0*(1d0+c2t**2)*cbe**2*Yt**2*
     .       ((phiA0T1T2*(-A0+T1-T2))/T2**2+phiA0T1T2/T2-
     .         ((-A0+T1-T2)*
     .            (tmp104/deltA0T1T2+
     .            (phiA0T1T2*tmp179)/deltA0T1T2))/T2-
     .         ((-A0-T1+T2)*
     .            ((phiA0T1T2*(A0-T1+T2))/deltA0T1T2+
     .            (T2*tmp96)/(deltA0T1T2*T1)))/T2+
     .         0.5d0*LogA0/T1-0.5d0*LogT2/T1+0.5d0*LogA0/T2-
     .         0.5d0*LogT1/T2+0.5d0*(A0-T1-T2)/(T1*T2)+
     .         0.5d0*(deltA0T1T2*
     .             ((phiA0T1T2*(A0-T1+T2))/deltA0T1T2+
     .             (T2*tmp96)/(deltA0T1T2*T1)))/T2**2-
     .         0.5d0*(deltA0T1T2*
     .             (phiA0T1T2/deltA0T1T2+
     .             ((LogA0-LogT2-T1/T2+(A0-T2)/T2)*T2)/
     .              (deltA0T1T2*T1)-
     .             (2d0*phiA0T1T2*(-A0-T1+T2)*
     .                (A0-T1+T2))/deltA0T1T2**2+
     .             ((A0-T1+T2)*
     .                (tmp104/deltA0T1T2+
     .                  (phiA0T1T2*tmp179)/deltA0T1T2))/
     .              deltA0T1T2+tmp96/(deltA0T1T2*T1)-
     .             (2d0*T2*(-A0-T1+T2)*tmp96)/
     .              (deltA0T1T2**2*T1)))/T2)))
      Dtc2t = -(0.5d0*((-5d0*B2+4d0*B2*LogB2+
     .        LogT1**2*(B2-T1)-5d0*T1+
     .        LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (-(cbe*hb*ht*sbe*
     .           ((mb*tmp54)/mt+0.25d0*s2b/s2t-
     .             0.5d0*(s2b*tmp27)/mt))+
     .        cbe**2*hb**2*
     .         (-(0.125d0*(1d0+c2b)/c2t)+
     .           0.125d0*(mb*s2b)/(mt*s2t)-
     .           0.125d0*((1d0+c2b)*Xb)/(mt*s2t))+
     .        ht**2*sbe**2*
     .         (0.125d0*(1d0-c2b)/c2t+0.125d0*(mb*s2b)/(mt*s2t)-
     .           0.125d0*((1d0-c2b)*Xt)/(mt*s2t)))))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (-(cbe*hb*ht*sbe*
     .         ((mb*tmp55)/mt-0.25d0*s2b/s2t-
     .           0.5d0*(s2b*tmp28)/mt))+
     .      cbe**2*hb**2*
     .       (0.125d0*(1d0+c2b)/c2t-0.125d0*(mb*s2b)/(mt*s2t)+
     .         0.125d0*((1d0+c2b)*Xb)/(mt*s2t))+
     .      ht**2*sbe**2*
     .       (-(0.125d0*(1d0-c2b)/c2t)-0.125d0*(mb*s2b)/(mt*s2t)+
     .         0.125d0*((1d0-c2b)*Xt)/(mt*s2t))))-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (-(cbe*hb*ht*sbe*
     .         ((mb*tmp55)/mt-0.25d0*s2b/s2t+
     .           0.5d0*(s2b*tmp27)/mt))+
     .      cbe**2*hb**2*
     .       (-(0.125d0*(1d0-c2b)/c2t)-0.125d0*(mb*s2b)/(mt*s2t)-
     .         0.125d0*((1d0-c2b)*Xb)/(mt*s2t))+
     .      ht**2*sbe**2*
     .       (0.125d0*(1d0+c2b)/c2t-0.125d0*(mb*s2b)/(mt*s2t)-
     .         0.125d0*((1d0+c2b)*Xt)/(mt*s2t))))-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (-(cbe*hb*ht*sbe*
     .         ((mb*tmp54)/mt+0.25d0*s2b/s2t+
     .           0.5d0*(s2b*tmp28)/mt))+
     .      cbe**2*hb**2*
     .       (0.125d0*(1d0-c2b)/c2t+0.125d0*(mb*s2b)/(mt*s2t)+
     .         0.125d0*((1d0-c2b)*Xb)/(mt*s2t))+
     .      ht**2*sbe**2*
     .       (-(0.125d0*(1d0+c2b)/c2t)+0.125d0*(mb*s2b)/(mt*s2t)+
     .         0.125d0*((1d0+c2b)*Xt)/(mt*s2t))))
      Dtc2t = Dtc2t-tmp89*
     .    (cbe*hb*ht*sbe*((mb*tmp54)/mt+0.25d0*s2b/s2t-
     .       0.5d0*(s2b*tmp43)/mt)+
     .      hb**2*sbe**2*(-(0.125d0*(1d0+c2b)/c2t)+
     .       0.125d0*(mb*s2b)/(mt*s2t)-
     .       0.125d0*((1d0+c2b)*Yb)/(mt*s2t))+
     .      cbe**2*ht**2*(0.125d0*(1d0-c2b)/c2t+
     .       0.125d0*(mb*s2b)/(mt*s2t)-
     .       0.125d0*((1d0-c2b)*Yt)/(mt*s2t)))-
     .   tmp109*(cbe*hb*ht*sbe*
     .       ((mb*tmp55)/mt-0.25d0*s2b/s2t-0.5d0*(s2b*tmp44)/mt)+
     .      hb**2*sbe**2*(0.125d0*(1d0+c2b)/c2t-
     .       0.125d0*(mb*s2b)/(mt*s2t)+
     .       0.125d0*((1d0+c2b)*Yb)/(mt*s2t))+
     .      cbe**2*ht**2*(-(0.125d0*(1d0-c2b)/c2t)-
     .       0.125d0*(mb*s2b)/(mt*s2t)+
     .       0.125d0*((1d0-c2b)*Yt)/(mt*s2t)))-
     .   tmp88*(cbe*hb*ht*sbe*
     .       ((mb*tmp55)/mt-0.25d0*s2b/s2t+0.5d0*(s2b*tmp43)/mt)+
     .      hb**2*sbe**2*(-(0.125d0*(1d0-c2b)/c2t)-
     .       0.125d0*(mb*s2b)/(mt*s2t)-
     .       0.125d0*((1d0-c2b)*Yb)/(mt*s2t))+
     .      cbe**2*ht**2*(0.125d0*(1d0+c2b)/c2t-
     .       0.125d0*(mb*s2b)/(mt*s2t)-
     .       0.125d0*((1d0+c2b)*Yt)/(mt*s2t)))-
     .   tmp108*(cbe*hb*ht*sbe*
     .       ((mb*tmp54)/mt+0.25d0*s2b/s2t+0.5d0*(s2b*tmp44)/mt)+
     .      hb**2*sbe**2*(0.125d0*(1d0-c2b)/c2t+
     .       0.125d0*(mb*s2b)/(mt*s2t)+
     .       0.125d0*((1d0-c2b)*Yb)/(mt*s2t))+
     .      cbe**2*ht**2*(-(0.125d0*(1d0+c2b)/c2t)+
     .       0.125d0*(mb*s2b)/(mt*s2t)+
     .       0.125d0*((1d0+c2b)*Yt)/(mt*s2t)))-
     .   0.25d0*(ht**2*((cbe**2*tmp114*Yt)/(mt*s2t)-
     .      (cbe**2*tmp93*Yt)/(mt*s2t)+
     .      0.5d0*(sbe**2*tmp115*Xt)/(mt*s2t)-
     .      0.5d0*(sbe**2*tmp94*Xt)/(mt*s2t)))
      DT1c2t = -(tmp204*tmp238)-tmp203*tmp239+
     .   (hb*ht*mb*mu*tmp240)/s2t+
     .   hb**2*(0.125d0*(B1*(1d0-c2b)*(-1d0+LogB1))/c2t+
     .      0.125d0*(B2*(1d0+c2b)*(-1d0+LogB2))/c2t+
     .      sbe**2*(0.25d0*(A0*(-1d0+LogA0))/c2t+
     .       0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT1))/c2t)+
     .      0.125d0*(B1*(1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT1))/c2t+
     .      0.125d0*(B2*(1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT1))/c2t)+
     .   tmp6*(-(b*(-1d0+Logb))-2d0*b*Logb-
     .      b*(-1d0+Logb)*(-1d0+LogT1)-(-1d0+Logmu2)*mu2-
     .      2d0*Logmu2*mu2-(-1d0+Logmu2)*(-1d0+LogT1)*mu2-
     .      2d0*LogT1*T1-(-b-mu2+T1)*tmp240+
     .      0.5d0*(deltT1bmu2*phiT1bmu2)/mu2-
     .      0.5d0*(Logmu2*LogT1*(b-mu2-T1))-
     .      0.5d0*(Logb*LogT1*(-b+mu2-T1))-
     .      0.5d0*(Logb*Logmu2*(-b-mu2+T1))+2.5d0*(b+mu2+T1))
     .    -0.5d0*(tmp200*tmp213)-0.5d0*(tmp199*tmp214)+
     .   ht**2*(-(0.125d0*(B1*(1d0+c2b)*(-1d0+LogB1))/c2t)-
     .      0.125d0*(B2*(1d0-c2b)*(-1d0+LogB2))/c2t+
     .      cbe**2*(-(0.25d0*(A0*(-1d0+LogA0))/c2t)-
     .       0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT1))/c2t)-
     .      0.125d0*(B1*(1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT1))/c2t-
     .      0.125d0*(B2*(1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT1))/c2t+
     .      (-1d0+LogT2)*T2*(1d0+0.5d0*(-1d0+Nc))+
     .      (-1d0+LogT1)*(-1d0+LogT2)*T2*(1d0+0.5d0*(-1d0+Nc))-
     .      0.25d0*((1d0+Nc)*tmp79))-
     .   0.25d0*(ht**2*(sbe**2*tmp220*Xt**2+
     .      2d0*cbe**2*tmp221*Yt**2+
     .      cbe**2*tmp46*
     .       (-1d0+4d0*LogT1-LogT1**2-(A0*LogA0)/T1+
     .         (2d0*A0*phiA0T1T1)/T1+(LogT1*(A0-2d0*T1))/T1+
     .         0.5d0*(deltA0T1T1*phiA0T1T1)/T1**2-
     .         0.5d0*(deltA0T1T1*
     .             ((A0*phiA0T1T1)/deltA0T1T1+
     .             (phiA0T1T1*tmp178)/deltA0T1T1+
     .             tmp80/deltA0T1T1+tmp87/deltA0T1T1))/T1)+
     .      0.5d0*(sbe**2*tmp30*tmp95)))
      DT2c2t = -(tmp206*tmp222)-tmp205*tmp223+
     .   (hb*ht*mb*mu*tmp225)/s2t+
     .   hb**2*(-(0.125d0*(B1*(1d0-c2b)*(-1d0+LogB1))/c2t)-
     .      0.125d0*(B2*(1d0+c2b)*(-1d0+LogB2))/c2t+
     .      sbe**2*(-(0.25d0*(A0*(-1d0+LogA0))/c2t)-
     .       0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT2))/c2t)-
     .      0.125d0*(B1*(1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT2))/c2t-
     .      0.125d0*(B2*(1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT2))/c2t)+
     .   tmp7*(-(b*(-1d0+Logb))-2d0*b*Logb-
     .      b*(-1d0+Logb)*(-1d0+LogT2)-(-1d0+Logmu2)*mu2-
     .      2d0*Logmu2*mu2-(-1d0+Logmu2)*(-1d0+LogT2)*mu2-
     .      2d0*LogT2*T2-(-b-mu2+T2)*tmp224+
     .      0.5d0*(delT2Tbmu2*phiT2bmu2)/mu2-
     .      0.5d0*(Logmu2*LogT2*(b-mu2-T2))-
     .      0.5d0*(Logb*LogT2*(-b+mu2-T2))-
     .      0.5d0*(Logb*Logmu2*(-b-mu2+T2))+2.5d0*(b+mu2+T2))
     .    -0.5d0*(tmp202*tmp217)-0.5d0*(tmp201*tmp218)+
     .   ht**2*(0.125d0*(B1*(1d0+c2b)*(-1d0+LogB1))/c2t+
     .      0.125d0*(B2*(1d0-c2b)*(-1d0+LogB2))/c2t+
     .      cbe**2*(0.25d0*(A0*(-1d0+LogA0))/c2t+
     .       0.25d0*(A0*(-1d0+LogA0)*(-1d0+LogT2))/c2t)+
     .      0.125d0*(B1*(1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT2))/c2t+
     .      0.125d0*(B2*(1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT2))/c2t+
     .      (-1d0+LogT1)*T1*(1d0+0.5d0*(-1d0+Nc))+
     .      (-1d0+LogT1)*(-1d0+LogT2)*T1*(1d0+0.5d0*(-1d0+Nc))-
     .      0.25d0*((1d0+Nc)*tmp98))-
     .   0.25d0*(ht**2*(sbe**2*tmp219*Xt**2+
     .      2d0*cbe**2*Yt**2*
     .       (-0.5d0+2d0*LogT2-(phiA0T1T2*(-A0-T1+T2))/T2+
     .         0.5d0*(LogA0*LogT1)-0.5d0*(LogA0*LogT2)-
     .         0.5d0*(LogT1*LogT2)+
     .         0.5d0*(deltA0T1T2*phiA0T1T2)/T2**2+
     .         0.5d0*(LogT1*(A0-T1-T2))/T2+
     .         0.5d0*(LogA0*(-A0+T1-T2))/T2-
     .         0.5d0*(deltA0T1T2*
     .             (tmp104/deltA0T1T2+
     .             (phiA0T1T2*tmp179)/deltA0T1T2))/T2)+
     .      0.5d0*(sbe**2*tmp116*tmp31)+
     .      cbe**2*tmp47*
     .       (-1d0+4d0*LogT2-LogT2**2-(A0*LogA0)/T2+
     .         (2d0*A0*phiA0T2T2)/T2+(LogT2*(A0-2d0*T2))/T2+
     .         0.5d0*(deltA0T2T2*phiA0T2T2)/T2**2-
     .         0.5d0*(deltA0T2T2*
     .             ((A0*phiA0T2T2)/deltA0T2T2+
     .             tmp107/deltA0T2T2+
     .             (phiA0T2T2*tmp182)/deltA0T2T2+
     .             tmp99/deltA0T2T2))/T2)))
      Dtb = tmp10*(-1d0+Logb+(-1d0+Logb)*(-1d0+Logt)+
     .      Logt+0.5d0*((b+t)*tmp208)
     .      +0.5d0*tmp209+0.5d0*tmp210)-
     .     tmp108*(-(0.125d0*(cbe**2*ht**2*s2b*s2t)/(mb*mt))-
     .      0.125d0*(hb**2*s2b*s2t*sbe**2)/(mb*mt)+
     .      0.5d0*(cbe*hb*ht*sbe*tmp56)/(mb*mt))-
     .   tmp109*(0.125d0*(cbe**2*ht**2*s2b*s2t)/(mb*mt)+
     .      0.125d0*(hb**2*s2b*s2t*sbe**2)/(mb*mt)+
     .      0.5d0*(cbe*hb*ht*sbe*tmp59)/(mb*mt))-
     .   (2d0*cbe*hb*ht*mt*sbe*
     .      (0.5d0-2d0*Logt+(phiA0bt*(-A0-b+t))/t-
     .      0.5d0*(LogA0*Logb)+0.5d0*(LogA0*Logt)+
     .      0.5d0*(Logb*Logt)-0.5d0*(Logb*(A0-b-t))/t-
     .      0.5d0*(LogA0*(-A0+b-t))/t+0.5d0*tmp210+
     .      0.5d0*(deltA0bt*
     .          ((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .            (b*tmp66)/(deltA0bt*t)))/b))/mb-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (-(0.125d0*(cbe**2*hb**2*s2b*s2t)/(mb*mt))-
     .      0.125d0*(ht**2*s2b*s2t*sbe**2)/(mb*mt)-
     .      0.5d0*(cbe*hb*ht*sbe*tmp56)/(mb*mt)))-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (-(0.125d0*(cbe**2*hb**2*s2b*s2t)/(mb*mt))-
     .      0.125d0*(ht**2*s2b*s2t*sbe**2)/(mb*mt)-
     .      0.5d0*(cbe*hb*ht*sbe*tmp56)/(mb*mt)))-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (0.125d0*(cbe**2*hb**2*s2b*s2t)/(mb*mt)+
     .      0.125d0*(ht**2*s2b*s2t*sbe**2)/(mb*mt)-
     .      0.5d0*(cbe*hb*ht*sbe*tmp59)/(mb*mt)))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (0.125d0*(cbe**2*hb**2*s2b*s2t)/(mb*mt)+
     .      0.125d0*(ht**2*s2b*s2t*sbe**2)/(mb*mt)-
     .      0.5d0*(cbe*hb*ht*sbe*tmp59)/(mb*mt)))
      Dtb = Dtb-tmp89*
     .    (-(0.125d0*(cbe**2*ht**2*s2b*s2t)/(mb*mt))-
     .      0.125d0*(hb**2*s2b*s2t*sbe**2)/(mb*mt)+
     .      0.5d0*(cbe*hb*ht*sbe*tmp56)/(mb*mt))-
     .   tmp88*(0.125d0*(cbe**2*ht**2*s2b*s2t)/(mb*mt)+
     .      0.125d0*(hb**2*s2b*s2t*sbe**2)/(mb*mt)+
     .      0.5d0*(cbe*hb*ht*sbe*tmp59)/(mb*mt))-
     .   (2d0*cbe*hb*ht*mb*sbe*
     .      (0.5d0-2d0*Logb+(phiA0bt*(-A0+b-t))/t+
     .      0.5d0*(LogA0*Logb)-0.5d0*(LogA0*Logt)+
     .      0.5d0*(Logb*Logt)-0.5d0*(Logt*(A0-b-t))/b-
     .      0.5d0*(deltA0bt*phiA0bt)/(b*t)-
     .      0.5d0*(LogA0*(-A0-b+t))/b+0.5d0*tmp209+
     .      0.5d0*(deltA0bt*
     .          ((b*phiA0bt*tmp174)/(deltA0bt*t)+
     .            tmp69/deltA0bt))/b))/mt-
     .   4d0*cbe*hb*ht*mb*mt*sbe*
     .    (-(phiA0bt/t)-(phiA0bt*(-A0-b+t))/(b*t)+
     .      ((-A0+b-t)*
     .       ((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .         (b*tmp66)/(deltA0bt*t)))/b+
     .      ((-A0-b+t)*
     .       ((b*phiA0bt*tmp174)/(deltA0bt*t)+tmp69/deltA0bt))/
     .       b-0.5d0*LogA0/b+0.5d0*Logt/b-0.5d0*LogA0/t+
     .      0.5d0*Logb/t-0.5d0*(A0-b-t)/(b*t)+0.5d0*tmp208-
     .      0.5d0*(deltA0bt*((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .          (b*tmp66)/(deltA0bt*t)))/b**2+
     .      0.5d0*(deltA0bt*((b*phiA0bt)/(deltA0bt*t)-
     .          (2d0*b*phiA0bt*(-A0+b-t)*(A0+b-t))/
     .           (deltA0bt**2*t)+(b*tmp64)/(deltA0bt*t)+
     .          tmp66/(deltA0bt*t)-
     .          (2d0*b*(-A0+b-t)*tmp66)/(deltA0bt**2*t)+
     .          ((A0+b-t)*
     .             ((b*phiA0bt*tmp174)/(deltA0bt*t)+
     .             tmp69/deltA0bt))/deltA0bt))/b)-
     .   (cbe*hb*ht*sbe*(-2d0*A0*LogA0-2d0*b*Logb-2d0*Logt*t-
     .      0.5d0*(Logb*Logt*(A0-b-t))-
     .      0.5d0*(LogA0*Logt*(-A0+b-t))+
     .      0.5d0*(deltA0bt*phiA0bt)/t-
     .      0.5d0*(LogA0*Logb*(-A0-b+t))+2.5d0*(A0+b+t)+
     .      0.5d0*tmp76))/(mb*mt)
      Dtb = Dtb+tmp9*
     .    (-2d0+3d0*Logb+(-1d0+Logb)*(-1d0+Logt)+3d0*Logt-
     .      Logb*Logt-(phiA0bt*(-A0+b-t))/t-
     .      (phiA0bt*(-A0-b+t))/t+
     .      0.5d0*(Logt*(A0-b-t))/b+
     .      0.5d0*(deltA0bt*phiA0bt)/(b*t)+
     .      0.5d0*(Logb*(A0-b-t))/t+
     .      0.5d0*(LogA0*(-A0+b-t))/t+
     .      0.5d0*(LogA0*(-A0-b+t))/b-
     .      0.5d0*(deltA0bt*((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .          (b*tmp66)/(deltA0bt*t)))/b-
     .      0.5d0*(deltA0bt*((b*phiA0bt*tmp174)/(deltA0bt*t)+
     .          tmp69/deltA0bt))/b-
     .      (A0-b-t)*(phiA0bt/t+
     .       (phiA0bt*(-A0-b+t))/(b*t)-
     .       ((-A0+b-t)*
     .          ((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .            (b*tmp66)/(deltA0bt*t)))/b-
     .       ((-A0-b+t)*
     .          ((b*phiA0bt*tmp174)/(deltA0bt*t)+
     .            tmp69/deltA0bt))/b+0.5d0*LogA0/b-
     .       0.5d0*Logt/b+0.5d0*LogA0/t-0.5d0*Logb/t+
     .       0.5d0*(A0-b-t)/(b*t)+
     .       0.5d0*(deltA0bt*
     .           ((b*phiA0bt*(A0+b-t))/(deltA0bt*t)+
     .             (b*tmp66)/(deltA0bt*t)))/b**2-
     .       0.5d0*(deltA0bt*
     .           ((b*phiA0bt)/(deltA0bt*t)-
     .             (2d0*b*phiA0bt*(-A0+b-t)*(A0+b-t))/
     .            (deltA0bt**2*t)+(b*tmp64)/(deltA0bt*t)+
     .             tmp66/(deltA0bt*t)-
     .             (2d0*b*(-A0+b-t)*tmp66)/(deltA0bt**2*t)+
     .             ((A0+b-t)*
     .              ((b*phiA0bt*tmp174)/(deltA0bt*t)+
     .                tmp69/deltA0bt))/deltA0bt))/b))
      DT1b = -((hb*ht*mu*s2t*tmp240)/mb)-
     .   2d0*hb*ht*mb*mu*s2t*
     .    (phiT1bmu2/mu2-
     .      ((b-mu2-T1)*
     .       ((phiT1bmu2*(b+mu2-T1))/deltT1bmu2+
     .         (mu2*tmp83)/(deltT1bmu2*T1)))/mu2-
     .      ((-b-mu2+T1)*
     .       ((phiT1bmu2*(-b+mu2+T1))/deltT1bmu2+
     .         (mu2*tmp90)/(b*deltT1bmu2)))/mu2+0.5d0*Logmu2/b-
     .      0.5d0*LogT1/b-0.5d0*Logb/T1+0.5d0*Logmu2/T1+
     .      0.5d0*(-b+mu2-T1)/(b*T1)-
     .      0.5d0*(deltT1bmu2*tmp242)/mu2)+
     .   tmp62*(1d0-3d0*Logb+Logb*Logmu2-
     .      (-1d0+Logb)*(-1d0+LogT1)+LogT1-Logmu2*LogT1+
     .      (phiT1bmu2*(b-mu2-T1))/mu2-
     .      (phiT1bmu2*(-b-mu2+T1))/mu2-
     .      0.5d0*(LogT1*(-b+mu2-T1))/b+
     .      0.5d0*(Logmu2*(b-mu2-T1))/T1+
     .      0.5d0*(Logb*(-b+mu2-T1))/T1-
     .      0.5d0*(Logmu2*(-b-mu2+T1))/b-
     .      (-b-mu2+T1)*
     .       (phiT1bmu2/mu2-
     .       ((b-mu2-T1)*
     .          ((phiT1bmu2*(b+mu2-T1))/deltT1bmu2+
     .            (mu2*tmp83)/(deltT1bmu2*T1)))/mu2-
     .       ((-b-mu2+T1)*
     .          ((phiT1bmu2*(-b+mu2+T1))/deltT1bmu2+
     .            (mu2*tmp90)/(b*deltT1bmu2)))/mu2+
     .       0.5d0*Logmu2/b-0.5d0*LogT1/b-0.5d0*Logb/T1+
     .       0.5d0*Logmu2/T1+0.5d0*(-b+mu2-T1)/(b*T1)-
     .       0.5d0*(deltT1bmu2*tmp242)/mu2)-
     .      0.5d0*(deltT1bmu2*
     .        ((phiT1bmu2*(b+mu2-T1))/deltT1bmu2+
     .          (mu2*tmp83)/(deltT1bmu2*T1)))/mu2+
     .      0.5d0*(deltT1bmu2*
     .        ((phiT1bmu2*(-b+mu2+T1))/deltT1bmu2+
     .          (mu2*tmp90)/(b*deltT1bmu2)))/mu2)-
     .   0.5d0*(tmp214*(-(cbe*hb*ht*sbe*
     .         ((mt*tmp56)/mb-0.5d0*(s2b*s2t)+
     .           0.5d0*(s2t*tmp128)/mb))+
     .      cbe**2*hb**2*
     .       (0.25d0*((1d0-c2b)*(1d0+c2t))-
     .         0.25d0*(mt*s2b*s2t)/mb-0.25d0*((1d0+c2t)*s2b*Xb)/mb)
     .       +ht**2*sbe**2*
     .       (0.25d0*((1d0+c2b)*(1d0-c2t))-
     .         0.25d0*(mt*s2b*s2t)/mb-0.25d0*((1d0-c2t)*s2b*Xt)/mb)
     .      ))
      DT1b = DT1b-tmp239*
     .    (cbe*hb*ht*sbe*((mt*tmp56)/mb-0.5d0*(s2b*s2t)+
     .       0.5d0*(s2t*tmp154)/mb)+
     .      hb**2*sbe**2*(0.25d0*((1d0-c2b)*(1d0+c2t))-
     .       0.25d0*(mt*s2b*s2t)/mb-0.25d0*((1d0+c2t)*s2b*Yb)/mb)+
     .      cbe**2*ht**2*
     .       (0.25d0*((1d0+c2b)*(1d0-c2t))-0.25d0*(mt*s2b*s2t)/mb-
     .       0.25d0*((1d0-c2t)*s2b*Yt)/mb))-
     .   tmp238*(cbe*hb*ht*sbe*
     .       ((mt*tmp59)/mb+0.5d0*(s2b*s2t)+0.5d0*(s2t*tmp155)/mb)+
     .      hb**2*sbe**2*
     .       (0.25d0*((1d0+c2b)*(1d0+c2t))+0.25d0*(mt*s2b*s2t)/mb+
     .       0.25d0*((1d0+c2t)*s2b*Yb)/mb)+
     .      cbe**2*ht**2*(0.25d0*((1d0-c2b)*(1d0-c2t))+
     .       0.25d0*(mt*s2b*s2t)/mb+0.25d0*((1d0-c2t)*s2b*Yt)/mb))-
     .     0.5d0*(tmp213*(-(cbe*hb*ht*sbe*
     .         ((mt*tmp59)/mb+0.5d0*(s2b*s2t)+
     .           0.5d0*(s2t*tmp129)/mb))+
     .      cbe**2*hb**2*
     .       (0.25d0*((1d0+c2b)*(1d0+c2t))+
     .         0.25d0*(mt*s2b*s2t)/mb+0.25d0*((1d0+c2t)*s2b*Xb)/mb)
     .       +ht**2*sbe**2*
     .       (0.25d0*((1d0-c2b)*(1d0-c2t))+
     .         0.25d0*(mt*s2b*s2t)/mb+0.25d0*((1d0-c2t)*s2b*Xt)/mb)
     .      ))
      DT2b = -((hb*ht*mu*s2t*tmp225)/mb)+
     .   tmp63*(1d0-3d0*Logb+Logb*Logmu2-
     .      (-1d0+Logb)*(-1d0+LogT2)+LogT2-Logmu2*LogT2+
     .      (phiT2bmu2*(b-mu2-T2))/mu2-
     .      (phiT2bmu2*(-b-mu2+T2))/mu2-
     .      0.5d0*(LogT2*(-b+mu2-T2))/b+
     .      0.5d0*(Logmu2*(b-mu2-T2))/T2+
     .      0.5d0*(Logb*(-b+mu2-T2))/T2-
     .      0.5d0*(Logmu2*(-b-mu2+T2))/b-
     .      0.5d0*(delT2Tbmu2*
     .        ((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .          (mu2*tmp102)/(delT2Tbmu2*T2)))/mu2+
     .      0.5d0*(delT2Tbmu2*
     .        ((phiT2bmu2*(-b+mu2+T2))/delT2Tbmu2+
     .          (mu2*tmp110)/(b*delT2Tbmu2)))/mu2-
     .      (-b-mu2+T2)*
     .       (phiT2bmu2/mu2-
     .       ((b-mu2-T2)*
     .          ((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .            (mu2*tmp102)/(delT2Tbmu2*T2)))/mu2-
     .       ((-b-mu2+T2)*
     .          ((phiT2bmu2*(-b+mu2+T2))/delT2Tbmu2+
     .            (mu2*tmp110)/(b*delT2Tbmu2)))/mu2+
     .       0.5d0*Logmu2/b-0.5d0*LogT2/b-0.5d0*Logb/T2+
     .       0.5d0*Logmu2/T2+0.5d0*(-b+mu2-T2)/(b*T2)-
     .       0.5d0*(delT2Tbmu2*tmp227)/mu2))-
     .   2d0*hb*ht*mb*mu*s2t*
     .    (-(phiT2bmu2/mu2)+
     .      ((b-mu2-T2)*
     .       ((phiT2bmu2*(b+mu2-T2))/delT2Tbmu2+
     .         (mu2*tmp102)/(delT2Tbmu2*T2)))/mu2+
     .      ((-b-mu2+T2)*
     .       ((phiT2bmu2*(-b+mu2+T2))/delT2Tbmu2+
     .         (mu2*tmp110)/(b*delT2Tbmu2)))/mu2-
     .      0.5d0*Logmu2/b+0.5d0*LogT2/b+0.5d0*Logb/T2-
     .      0.5d0*Logmu2/T2-0.5d0*(-b+mu2-T2)/(b*T2)+
     .      0.5d0*(delT2Tbmu2*tmp227)/mu2)-
     .   0.5d0*(tmp218*(-(cbe*hb*ht*sbe*
     .         ((mt*tmp59)/mb+0.5d0*(s2b*s2t)-
     .           0.5d0*(s2t*tmp128)/mb))+
     .      cbe**2*hb**2*
     .       (0.25d0*((1d0-c2b)*(1d0-c2t))+
     .         0.25d0*(mt*s2b*s2t)/mb-0.25d0*((1d0-c2t)*s2b*Xb)/mb)
     .       +ht**2*sbe**2*
     .       (0.25d0*((1d0+c2b)*(1d0+c2t))+
     .         0.25d0*(mt*s2b*s2t)/mb-0.25d0*((1d0+c2t)*s2b*Xt)/mb)
     .      ))
      DT2b = DT2b-tmp223*
     .    (cbe*hb*ht*sbe*((mt*tmp59)/mb+0.5d0*(s2b*s2t)-
     .       0.5d0*(s2t*tmp154)/mb)+
     .      hb**2*sbe**2*(0.25d0*((1d0-c2b)*(1d0-c2t))+
     .       0.25d0*(mt*s2b*s2t)/mb-0.25d0*((1d0-c2t)*s2b*Yb)/mb)+
     .      cbe**2*ht**2*
     .       (0.25d0*((1d0+c2b)*(1d0+c2t))+0.25d0*(mt*s2b*s2t)/mb-
     .       0.25d0*((1d0+c2t)*s2b*Yt)/mb))-
     .   tmp222*(cbe*hb*ht*sbe*
     .       ((mt*tmp56)/mb-0.5d0*(s2b*s2t)-0.5d0*(s2t*tmp155)/mb)+
     .      hb**2*sbe**2*
     .       (0.25d0*((1d0+c2b)*(1d0-c2t))-0.25d0*(mt*s2b*s2t)/mb+
     .       0.25d0*((1d0-c2t)*s2b*Yb)/mb)+
     .      cbe**2*ht**2*(0.25d0*((1d0-c2b)*(1d0+c2t))-
     .       0.25d0*(mt*s2b*s2t)/mb+0.25d0*((1d0+c2t)*s2b*Yt)/mb))-
     .     0.5d0*(tmp217*(-(cbe*hb*ht*sbe*
     .         ((mt*tmp56)/mb-0.5d0*(s2b*s2t)-
     .           0.5d0*(s2t*tmp129)/mb))+
     .      cbe**2*hb**2*
     .       (0.25d0*((1d0+c2b)*(1d0-c2t))-
     .         0.25d0*(mt*s2b*s2t)/mb+0.25d0*((1d0-c2t)*s2b*Xb)/mb)
     .       +ht**2*sbe**2*
     .       (0.25d0*((1d0-c2b)*(1d0+c2t))-
     .         0.25d0*(mt*s2b*s2t)/mb+0.25d0*((1d0+c2t)*s2b*Xt)/mb)
     .      ))
      DB1t = -(tmp196*(-0.5d0+2d0*LogB1-
     .      (phiA0B1T2*(-A0+B1-T2))/T2-0.5d0*(LogA0*LogB1)+
     .      0.5d0*(LogA0*LogT2)-0.5d0*(LogB1*LogT2)+
     .      0.5d0*(LogT2*(A0-B1-T2))/B1+
     .      0.5d0*(deltA0B1T2*phiA0B1T2)/(B1*T2)+
     .      0.5d0*(LogA0*(-A0-B1+T2))/B1-
     .      0.5d0*(deltA0B1T2*
     .          (tmp105/deltA0B1T2+
     .            (B1*phiA0B1T2*tmp180)/(deltA0B1T2*T2)))/B1))-
     .   0.5d0*(tmp184*tmp211)-0.5d0*(tmp188*tmp215)-
     .   2d0*hb*ht*mt*mu*s2b*
     .    (phiB1tmu2/mu2-
     .      ((B1-mu2-t)*
     .       ((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .         (mu2*tmp67)/(deltB1tmu2*t)))/mu2-
     .      ((-B1-mu2+t)*
     .       ((phiB1tmu2*(-B1+mu2+t))/deltB1tmu2+
     .         (mu2*tmp71)/(B1*deltB1tmu2)))/mu2+
     .      0.5d0*Logmu2/B1-0.5d0*Logt/B1-0.5d0*LogB1/t+
     .      0.5d0*Logmu2/t+0.5d0*(-B1+mu2-t)/(B1*t)-
     .      0.5d0*(deltB1tmu2*tmp233)/mu2)-
     .   (hb*ht*mu*s2b*(-0.5d0+2d0*LogB1-
     .      (phiB1tmu2*(B1-mu2-t))/mu2-
     .      0.5d0*(LogB1*Logmu2)-0.5d0*(LogB1*Logt)+
     .      0.5d0*(Logmu2*Logt)+0.5d0*(Logt*(-B1+mu2-t))/B1+
     .      0.5d0*(Logmu2*(-B1-mu2+t))/B1-
     .      0.5d0*(deltB1tmu2*
     .          ((phiB1tmu2*(-B1+mu2+t))/deltB1tmu2+
     .            (mu2*tmp71)/(B1*deltB1tmu2)))/mu2))/mt+
     .   tmp61*(1d0+LogB1-LogB1*Logmu2-
     .      (-1d0+LogB1)*(-1d0+Logt)-3d0*Logt+Logmu2*Logt-
     .      (phiB1tmu2*(B1-mu2-t))/mu2+
     .      (phiB1tmu2*(-B1-mu2+t))/mu2+
     .      0.5d0*(Logt*(-B1+mu2-t))/B1-
     .      0.5d0*(Logmu2*(B1-mu2-t))/t-
     .      0.5d0*(LogB1*(-B1+mu2-t))/t+
     .      0.5d0*(Logmu2*(-B1-mu2+t))/B1-
     .      (B1-mu2-t)*
     .       (phiB1tmu2/mu2-
     .       ((B1-mu2-t)*
     .          ((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .            (mu2*tmp67)/(deltB1tmu2*t)))/mu2-
     .       ((-B1-mu2+t)*
     .          ((phiB1tmu2*(-B1+mu2+t))/deltB1tmu2+
     .            (mu2*tmp71)/(B1*deltB1tmu2)))/mu2+
     .       0.5d0*Logmu2/B1-0.5d0*Logt/B1-0.5d0*LogB1/t+
     .       0.5d0*Logmu2/t+0.5d0*(-B1+mu2-t)/(B1*t)-
     .       0.5d0*(deltB1tmu2*tmp233)/mu2)+
     .      0.5d0*(deltB1tmu2*
     .        ((phiB1tmu2*(B1+mu2-t))/deltB1tmu2+
     .          (mu2*tmp67)/(deltB1tmu2*t)))/mu2-
     .      0.5d0*(deltB1tmu2*
     .        ((phiB1tmu2*(-B1+mu2+t))/deltB1tmu2+
     .          (mu2*tmp71)/(B1*deltB1tmu2)))/mu2)
      DB1t = DB1t-tmp192*
     .    (-0.5d0+2d0*LogB1-(phiA0B1T1*(-A0+B1-T1))/T1-
     .      0.5d0*(LogA0*LogB1)+0.5d0*(LogA0*LogT1)-
     .      0.5d0*(LogB1*LogT1)+0.5d0*(LogT1*(A0-B1-T1))/B1+
     .      0.5d0*(deltA0B1T1*phiA0B1T1)/(B1*T1)+
     .      0.5d0*(LogA0*(-A0-B1+T1))/B1-
     .      0.5d0*(deltA0B1T1*
     .        ((B1*phiA0B1T1*tmp176)/(deltA0B1T1*T1)+
     .          tmp85/deltA0B1T1))/B1)
      DB2t = -(tmp195*(-0.5d0+2d0*LogB2-
     .      (phiA0B2T2*(-A0+B2-T2))/T2-0.5d0*(LogA0*LogB2)+
     .      0.5d0*(LogA0*LogT2)-0.5d0*(LogB2*LogT2)+
     .      0.5d0*(LogT2*(A0-B2-T2))/B2+
     .      0.5d0*(deltA0B2T2*phiA0B2T2)/(B2*T2)+
     .      0.5d0*(LogA0*(-A0-B2+T2))/B2-
     .      0.5d0*(deltA0B2T2*
     .          (tmp106/deltA0B2T2+
     .            (B2*phiA0B2T2*tmp181)/(deltA0B2T2*T2)))/B2))-
     .   0.5d0*(tmp183*tmp212)-0.5d0*(tmp187*tmp216)-
     .   2d0*hb*ht*mt*mu*s2b*
     .    (-(phiB2tmu2/mu2)+
     .      ((B2-mu2-t)*
     .       ((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .         (mu2*tmp68)/(deltB2tmu2*t)))/mu2+
     .      ((-B2-mu2+t)*
     .       ((phiB2tmu2*(-B2+mu2+t))/deltB2tmu2+
     .         (mu2*tmp72)/(B2*deltB2tmu2)))/mu2-
     .      0.5d0*Logmu2/B2+0.5d0*Logt/B2+0.5d0*LogB2/t-
     .      0.5d0*Logmu2/t-0.5d0*(-B2+mu2-t)/(B2*t)+
     .      0.5d0*(deltB2tmu2*tmp237)/mu2)+
     .   tmp60*(1d0+LogB2-LogB2*Logmu2-
     .      (-1d0+LogB2)*(-1d0+Logt)-3d0*Logt+Logmu2*Logt-
     .      (phiB2tmu2*(B2-mu2-t))/mu2+
     .      (phiB2tmu2*(-B2-mu2+t))/mu2+
     .      0.5d0*(Logt*(-B2+mu2-t))/B2-
     .      0.5d0*(Logmu2*(B2-mu2-t))/t-
     .      0.5d0*(LogB2*(-B2+mu2-t))/t+
     .      0.5d0*(Logmu2*(-B2-mu2+t))/B2-
     .      (B2-mu2-t)*
     .       (phiB2tmu2/mu2-
     .       ((B2-mu2-t)*
     .          ((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .            (mu2*tmp68)/(deltB2tmu2*t)))/mu2-
     .       ((-B2-mu2+t)*
     .          ((phiB2tmu2*(-B2+mu2+t))/deltB2tmu2+
     .            (mu2*tmp72)/(B2*deltB2tmu2)))/mu2+
     .       0.5d0*Logmu2/B2-0.5d0*Logt/B2-0.5d0*LogB2/t+
     .       0.5d0*Logmu2/t+0.5d0*(-B2+mu2-t)/(B2*t)-
     .       0.5d0*(deltB2tmu2*tmp237)/mu2)+
     .      0.5d0*(deltB2tmu2*
     .        ((phiB2tmu2*(B2+mu2-t))/deltB2tmu2+
     .          (mu2*tmp68)/(deltB2tmu2*t)))/mu2-
     .      0.5d0*(deltB2tmu2*
     .        ((phiB2tmu2*(-B2+mu2+t))/deltB2tmu2+
     .          (mu2*tmp72)/(B2*deltB2tmu2)))/mu2)-
     .   (hb*ht*mu*s2b*(0.5d0-2d0*LogB2+
     .      (phiB2tmu2*(B2-mu2-t))/mu2+
     .      0.5d0*(LogB2*Logmu2)+0.5d0*(LogB2*Logt)-
     .      0.5d0*(Logmu2*Logt)-0.5d0*(Logt*(-B2+mu2-t))/B2-
     .      0.5d0*(Logmu2*(-B2-mu2+t))/B2+
     .      0.5d0*(deltB2tmu2*
     .          ((phiB2tmu2*(-B2+mu2+t))/deltB2tmu2+
     .            (mu2*tmp72)/(B2*deltB2tmu2)))/mu2))/mt
      DB2t = DB2t-tmp191*
     .    (-0.5d0+2d0*LogB2-(phiA0B2T1*(-A0+B2-T1))/T1-
     .      0.5d0*(LogA0*LogB2)+0.5d0*(LogA0*LogT1)-
     .      0.5d0*(LogB2*LogT1)+0.5d0*(LogT1*(A0-B2-T1))/B2+
     .      0.5d0*(deltA0B2T1*phiA0B2T1)/(B2*T1)+
     .      0.5d0*(LogA0*(-A0-B2+T1))/B2-
     .      0.5d0*(deltA0B2T1*
     .        ((B2*phiA0B2T1*tmp177)/(deltA0B2T1*T1)+
     .          tmp86/deltA0B2T1))/B2)
      DT1B1 = ht**2*(0.25d0*((1d0+c2b)*(1d0-c2t))+
     .      0.25d0*((1d0+c2b)*(1d0-c2t)*(-1d0+LogB1))+
     .      0.25d0*((1d0+c2b)*(1d0-c2t)*(-1d0+LogT1))+
     .      0.25d0*((1d0+c2b)*(1d0-c2t)*(-1d0+LogB1)*(-1d0+LogT1)))+
     .     hb**2*(0.25d0*((1d0-c2b)*(1d0+c2t))+
     .      0.25d0*((1d0-c2b)*(1d0+c2t)*(-1d0+LogB1))+
     .      0.25d0*((1d0-c2b)*(1d0+c2t)*(-1d0+LogT1))+
     .      0.25d0*((1d0-c2b)*(1d0+c2t)*(-1d0+LogB1)*(-1d0+LogT1)))-
     .     0.5d0*(((-2d0*B1*(LogB1-LogT1))/((1d0-B1/T1)*T1**2)+
     .      (-2d0-2d0*LogB1)/T1+(2d0*LogT1)/T1-
     .      (2d0*(LogB1-LogT1))/((1d0-B1/T1)*T1)+
     .      (2d0*B1*(LogB1-LogT1)*(-B1+T1))/
     .       ((1d0-B1/T1)**2*T1**3)+
     .      (2d0*(-B1+T1))/((1d0-B1/T1)*T1**2)+
     .      (2d0*(LogB1-LogT1)*(-B1+T1))/((1d0-B1/T1)*T1**2))*
     .      tmp186)-tmp194*
     .    (phiA0B1T1/T1+(phiA0B1T1*(-A0-B1+T1))/(B1*T1)-
     .      ((-A0+B1-T1)*
     .       ((B1*phiA0B1T1*(A0+B1-T1))/(deltA0B1T1*T1)+
     .         (B1*tmp81)/(deltA0B1T1*T1)))/B1-
     .      ((-A0-B1+T1)*
     .       ((B1*phiA0B1T1*tmp176)/(deltA0B1T1*T1)+
     .         tmp85/deltA0B1T1))/B1+0.5d0*LogA0/B1-
     .      0.5d0*LogT1/B1+0.5d0*LogA0/T1-0.5d0*LogB1/T1+
     .      0.5d0*(A0-B1-T1)/(B1*T1)+
     .      0.5d0*(deltA0B1T1*
     .        ((B1*phiA0B1T1*(A0+B1-T1))/(deltA0B1T1*T1)+
     .          (B1*tmp81)/(deltA0B1T1*T1)))/B1**2-
     .      0.5d0*(deltA0B1T1*
     .        ((B1*phiA0B1T1)/(deltA0B1T1*T1)-
     .          (2d0*B1*phiA0B1T1*(-A0+B1-T1)*(A0+B1-T1))/
     .           (deltA0B1T1**2*T1)+
     .          (B1*((A0-B1)/B1+LogA0-LogB1-T1/B1))/
     .           (deltA0B1T1*T1)+tmp81/(deltA0B1T1*T1)-
     .          (2d0*B1*(-A0+B1-T1)*tmp81)/
     .           (deltA0B1T1**2*T1)+
     .          ((A0+B1-T1)*
     .             ((B1*phiA0B1T1*tmp176)/(deltA0B1T1*T1)+
     .             tmp85/deltA0B1T1))/deltA0B1T1))/B1)
      DT2B1 = hb**2*(0.25d0*((1d0-c2b)*(1d0-c2t))+
     .      0.25d0*((1d0-c2b)*(1d0-c2t)*(-1d0+LogB1))+
     .      0.25d0*((1d0-c2b)*(1d0-c2t)*(-1d0+LogT2))+
     .      0.25d0*((1d0-c2b)*(1d0-c2t)*(-1d0+LogB1)*(-1d0+LogT2)))+
     .     ht**2*(0.25d0*((1d0+c2b)*(1d0+c2t))+
     .      0.25d0*((1d0+c2b)*(1d0+c2t)*(-1d0+LogB1))+
     .      0.25d0*((1d0+c2b)*(1d0+c2t)*(-1d0+LogT2))+
     .      0.25d0*((1d0+c2b)*(1d0+c2t)*(-1d0+LogB1)*(-1d0+LogT2)))-
     .     tmp198*(phiA0B1T2/T2+
     .      (phiA0B1T2*(-A0-B1+T2))/(B1*T2)-
     .      ((-A0+B1-T2)*
     .       ((B1*phiA0B1T2*(A0+B1-T2))/(deltA0B1T2*T2)+
     .         (B1*tmp100)/(deltA0B1T2*T2)))/B1-
     .      ((-A0-B1+T2)*
     .       (tmp105/deltA0B1T2+
     .         (B1*phiA0B1T2*tmp180)/(deltA0B1T2*T2)))/B1+
     .      0.5d0*LogA0/B1-0.5d0*LogT2/B1+0.5d0*LogA0/T2-
     .      0.5d0*LogB1/T2+0.5d0*(A0-B1-T2)/(B1*T2)+
     .      0.5d0*(deltA0B1T2*
     .        ((B1*phiA0B1T2*(A0+B1-T2))/(deltA0B1T2*T2)+
     .          (B1*tmp100)/(deltA0B1T2*T2)))/B1**2-
     .      0.5d0*(deltA0B1T2*
     .        ((B1*phiA0B1T2)/(deltA0B1T2*T2)-
     .          (2d0*B1*phiA0B1T2*(-A0+B1-T2)*(A0+B1-T2))/
     .           (deltA0B1T2**2*T2)+
     .          (B1*((A0-B1)/B1+LogA0-LogB1-T2/B1))/
     .           (deltA0B1T2*T2)+tmp100/(deltA0B1T2*T2)-
     .          (2d0*B1*(-A0+B1-T2)*tmp100)/
     .           (deltA0B1T2**2*T2)+
     .          ((A0+B1-T2)*
     .             (tmp105/deltA0B1T2+
     .             (B1*phiA0B1T2*tmp180)/(deltA0B1T2*T2)))/
     .           deltA0B1T2))/B1)-
     .   0.5d0*(((-2d0*B1*(LogB1-LogT2))/((1d0-B1/T2)*T2**2)+
     .      (-2d0-2d0*LogB1)/T2+(2d0*LogT2)/T2-
     .      (2d0*(LogB1-LogT2))/((1d0-B1/T2)*T2)+
     .      (2d0*B1*(LogB1-LogT2)*(-B1+T2))/
     .       ((1d0-B1/T2)**2*T2**3)+
     .      (2d0*(-B1+T2))/((1d0-B1/T2)*T2**2)+
     .      (2d0*(LogB1-LogT2)*(-B1+T2))/((1d0-B1/T2)*T2**2))*
     .      tmp190)
      DT1B2 = ht**2*(0.25d0*((1d0-c2b)*(1d0-c2t))+
     .      0.25d0*((1d0-c2b)*(1d0-c2t)*(-1d0+LogB2))+
     .      0.25d0*((1d0-c2b)*(1d0-c2t)*(-1d0+LogT1))+
     .      0.25d0*((1d0-c2b)*(1d0-c2t)*(-1d0+LogB2)*(-1d0+LogT1)))+
     .     hb**2*(0.25d0*((1d0+c2b)*(1d0+c2t))+
     .      0.25d0*((1d0+c2b)*(1d0+c2t)*(-1d0+LogB2))+
     .      0.25d0*((1d0+c2b)*(1d0+c2t)*(-1d0+LogT1))+
     .      0.25d0*((1d0+c2b)*(1d0+c2t)*(-1d0+LogB2)*(-1d0+LogT1)))-
     .     0.5d0*(((-2d0*B2*(LogB2-LogT1))/((1d0-B2/T1)*T1**2)+
     .      (-2d0-2d0*LogB2)/T1+(2d0*LogT1)/T1-
     .      (2d0*(LogB2-LogT1))/((1d0-B2/T1)*T1)+
     .      (2d0*B2*(LogB2-LogT1)*(-B2+T1))/
     .       ((1d0-B2/T1)**2*T1**3)+
     .      (2d0*(-B2+T1))/((1d0-B2/T1)*T1**2)+
     .      (2d0*(LogB2-LogT1)*(-B2+T1))/((1d0-B2/T1)*T1**2))*
     .      tmp185)-tmp193*
     .    (phiA0B2T1/T1+(phiA0B2T1*(-A0-B2+T1))/(B2*T1)-
     .      ((-A0+B2-T1)*
     .       ((B2*phiA0B2T1*(A0+B2-T1))/(deltA0B2T1*T1)+
     .         (B2*tmp82)/(deltA0B2T1*T1)))/B2-
     .      ((-A0-B2+T1)*
     .       ((B2*phiA0B2T1*tmp177)/(deltA0B2T1*T1)+
     .         tmp86/deltA0B2T1))/B2+0.5d0*LogA0/B2-
     .      0.5d0*LogT1/B2+0.5d0*LogA0/T1-0.5d0*LogB2/T1+
     .      0.5d0*(A0-B2-T1)/(B2*T1)+
     .      0.5d0*(deltA0B2T1*
     .        ((B2*phiA0B2T1*(A0+B2-T1))/(deltA0B2T1*T1)+
     .          (B2*tmp82)/(deltA0B2T1*T1)))/B2**2-
     .      0.5d0*(deltA0B2T1*
     .        ((B2*phiA0B2T1)/(deltA0B2T1*T1)-
     .          (2d0*B2*phiA0B2T1*(-A0+B2-T1)*(A0+B2-T1))/
     .           (deltA0B2T1**2*T1)+
     .          (B2*((A0-B2)/B2+LogA0-LogB2-T1/B2))/
     .           (deltA0B2T1*T1)+tmp82/(deltA0B2T1*T1)-
     .          (2d0*B2*(-A0+B2-T1)*tmp82)/
     .           (deltA0B2T1**2*T1)+
     .          ((A0+B2-T1)*
     .             ((B2*phiA0B2T1*tmp177)/(deltA0B2T1*T1)+
     .             tmp86/deltA0B2T1))/deltA0B2T1))/B2)
      DT2B2 = hb**2*(0.25d0*((1d0+c2b)*(1d0-c2t))+
     .      0.25d0*((1d0+c2b)*(1d0-c2t)*(-1d0+LogB2))+
     .      0.25d0*((1d0+c2b)*(1d0-c2t)*(-1d0+LogT2))+
     .      0.25d0*((1d0+c2b)*(1d0-c2t)*(-1d0+LogB2)*(-1d0+LogT2)))+
     .     ht**2*(0.25d0*((1d0-c2b)*(1d0+c2t))+
     .      0.25d0*((1d0-c2b)*(1d0+c2t)*(-1d0+LogB2))+
     .      0.25d0*((1d0-c2b)*(1d0+c2t)*(-1d0+LogT2))+
     .      0.25d0*((1d0-c2b)*(1d0+c2t)*(-1d0+LogB2)*(-1d0+LogT2)))-
     .     tmp197*(phiA0B2T2/T2+
     .      (phiA0B2T2*(-A0-B2+T2))/(B2*T2)-
     .      ((-A0+B2-T2)*
     .       ((B2*phiA0B2T2*(A0+B2-T2))/(deltA0B2T2*T2)+
     .         (B2*tmp101)/(deltA0B2T2*T2)))/B2-
     .      ((-A0-B2+T2)*
     .       (tmp106/deltA0B2T2+
     .         (B2*phiA0B2T2*tmp181)/(deltA0B2T2*T2)))/B2+
     .      0.5d0*LogA0/B2-0.5d0*LogT2/B2+0.5d0*LogA0/T2-
     .      0.5d0*LogB2/T2+0.5d0*(A0-B2-T2)/(B2*T2)+
     .      0.5d0*(deltA0B2T2*
     .        ((B2*phiA0B2T2*(A0+B2-T2))/(deltA0B2T2*T2)+
     .          (B2*tmp101)/(deltA0B2T2*T2)))/B2**2-
     .      0.5d0*(deltA0B2T2*
     .        ((B2*phiA0B2T2)/(deltA0B2T2*T2)-
     .          (2d0*B2*phiA0B2T2*(-A0+B2-T2)*(A0+B2-T2))/
     .           (deltA0B2T2**2*T2)+
     .          (B2*((A0-B2)/B2+LogA0-LogB2-T2/B2))/
     .           (deltA0B2T2*T2)+tmp101/(deltA0B2T2*T2)-
     .          (2d0*B2*(-A0+B2-T2)*tmp101)/
     .           (deltA0B2T2**2*T2)+
     .          ((A0+B2-T2)*
     .             (tmp106/deltA0B2T2+
     .             (B2*phiA0B2T2*tmp181)/(deltA0B2T2*T2)))/
     .           deltA0B2T2))/B2)-
     .   0.5d0*(((-2d0*B2*(LogB2-LogT2))/((1d0-B2/T2)*T2**2)+
     .      (-2d0-2d0*LogB2)/T2+(2d0*LogT2)/T2-
     .      (2d0*(LogB2-LogT2))/((1d0-B2/T2)*T2)+
     .      (2d0*B2*(LogB2-LogT2)*(-B2+T2))/
     .       ((1d0-B2/T2)**2*T2**3)+
     .      (2d0*(-B2+T2))/((1d0-B2/T2)*T2**2)+
     .      (2d0*(LogB2-LogT2)*(-B2+T2))/((1d0-B2/T2)*T2**2))*
     .      tmp189)
      Dbc2t = tmp7*(2d0*b*Logb+(-1d0+Logmu2)*mu2+
     .      (-1d0+Logb)*(-1d0+Logmu2)*mu2+2d0*Logmu2*mu2-
     .      (-1d0+LogT2)*T2-(-1d0+Logb)*(-1d0+LogT2)*T2+
     .      2d0*LogT2*T2-0.5d0*(delT2Tbmu2*phiT2bmu2)/mu2+
     .      0.5d0*(Logmu2*LogT2*(b-mu2-T2))+
     .      0.5d0*(Logb*LogT2*(-b+mu2-T2))+
     .      0.5d0*(Logb*Logmu2*(-b-mu2+T2))-
     .      2.5d0*(b+mu2+T2)-
     .      (-b-mu2+T2)*
     .       (-0.5d0+2d0*Logb-(phiT2bmu2*(b-mu2-T2))/mu2-
     .       0.5d0*(Logb*Logmu2)-0.5d0*(Logb*LogT2)+
     .       0.5d0*(Logmu2*LogT2)+
     .       0.5d0*(LogT2*(-b+mu2-T2))/b+
     .       0.5d0*(Logmu2*(-b-mu2+T2))/b-
     .       0.5d0*(delT2Tbmu2*
     .           ((phiT2bmu2*(-b+mu2+T2))/delT2Tbmu2+
     .             (mu2*tmp110)/(b*delT2Tbmu2)))/mu2))+
     .   0.5d0*(hb*ht*mu*tmp113)/(mb*s2t)+
     .   tmp6*(2d0*b*Logb+(-1d0+Logmu2)*mu2+
     .      (-1d0+Logb)*(-1d0+Logmu2)*mu2+2d0*Logmu2*mu2-
     .      (-1d0+LogT1)*T1-(-1d0+Logb)*(-1d0+LogT1)*T1+
     .      2d0*LogT1*T1-0.5d0*(deltT1bmu2*phiT1bmu2)/mu2+
     .      0.5d0*(Logmu2*LogT1*(b-mu2-T1))+
     .      0.5d0*(Logb*LogT1*(-b+mu2-T1))+
     .      0.5d0*(Logb*Logmu2*(-b-mu2+T1))-
     .      2.5d0*(b+mu2+T1)-
     .      (-b-mu2+T1)*
     .       (-0.5d0+2d0*Logb-(phiT1bmu2*(b-mu2-T1))/mu2-
     .       0.5d0*(Logb*Logmu2)-0.5d0*(Logb*LogT1)+
     .       0.5d0*(Logmu2*LogT1)+
     .       0.5d0*(LogT1*(-b+mu2-T1))/b+
     .       0.5d0*(Logmu2*(-b-mu2+T1))/b-
     .       0.5d0*(deltT1bmu2*
     .           ((phiT1bmu2*(-b+mu2+T1))/deltT1bmu2+
     .             (mu2*tmp90)/(b*deltT1bmu2)))/mu2))+
     .   (hb*ht*mb*mu*(-((phiT1bmu2*(b-mu2-T1))/mu2)+
     .      (phiT2bmu2*(b-mu2-T2))/mu2-0.5d0*(Logb*LogT1)+
     .      0.5d0*(Logmu2*LogT1)+0.5d0*(Logb*LogT2)-
     .      0.5d0*(Logmu2*LogT2)+0.5d0*(LogT1*(-b+mu2-T1))/b+
     .      0.5d0*(Logmu2*(-b-mu2+T1))/b-
     .      0.5d0*(LogT2*(-b+mu2-T2))/b-
     .      0.5d0*(Logmu2*(-b-mu2+T2))/b+
     .      0.5d0*(delT2Tbmu2*
     .          ((phiT2bmu2*(-b+mu2+T2))/delT2Tbmu2+
     .            (mu2*tmp110)/(b*delT2Tbmu2)))/mu2-
     .      0.5d0*(deltT1bmu2*
     .          ((phiT1bmu2*(-b+mu2+T1))/deltT1bmu2+
     .            (mu2*tmp90)/(b*deltT1bmu2)))/mu2))/s2t
      Dbc2t = Dbc2t-0.5d0*
     .    ((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (-(cbe*hb*ht*sbe*
     .         ((mt*tmp55)/mb-0.25d0*s2b/s2t-
     .           0.25d0*tmp129/(mb*s2t)))+
     .      cbe**2*hb**2*
     .       (0.125d0*(1d0+c2b)/c2t-0.125d0*(mt*s2b)/(mb*s2t)+
     .         0.125d0*(s2b*Xb)/(c2t*mb))+
     .      ht**2*sbe**2*
     .       (-(0.125d0*(1d0-c2b)/c2t)-0.125d0*(mt*s2b)/(mb*s2t)-
     .         0.125d0*(s2b*Xt)/(c2t*mb))))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (-(cbe*hb*ht*sbe*
     .         ((mt*tmp55)/mb-0.25d0*s2b/s2t+
     .           0.25d0*tmp128/(mb*s2t)))+
     .      cbe**2*hb**2*
     .       (-(0.125d0*(1d0-c2b)/c2t)-0.125d0*(mt*s2b)/(mb*s2t)+
     .         0.125d0*(s2b*Xb)/(c2t*mb))+
     .      ht**2*sbe**2*
     .       (0.125d0*(1d0+c2b)/c2t-0.125d0*(mt*s2b)/(mb*s2t)-
     .         0.125d0*(s2b*Xt)/(c2t*mb))))-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (-(cbe*hb*ht*sbe*
     .         ((mt*tmp54)/mb+0.25d0*s2b/s2t+
     .           0.25d0*tmp129/(mb*s2t)))+
     .      cbe**2*hb**2*
     .       (-(0.125d0*(1d0+c2b)/c2t)+0.125d0*(mt*s2b)/(mb*s2t)-
     .         0.125d0*(s2b*Xb)/(c2t*mb))+
     .      ht**2*sbe**2*
     .       (0.125d0*(1d0-c2b)/c2t+0.125d0*(mt*s2b)/(mb*s2t)+
     .         0.125d0*(s2b*Xt)/(c2t*mb))))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (-(cbe*hb*ht*sbe*
     .         ((mt*tmp54)/mb+0.25d0*s2b/s2t-
     .           0.25d0*tmp128/(mb*s2t)))+
     .      cbe**2*hb**2*
     .       (0.125d0*(1d0-c2b)/c2t+0.125d0*(mt*s2b)/(mb*s2t)-
     .         0.125d0*(s2b*Xb)/(c2t*mb))+
     .      ht**2*sbe**2*
     .       (-(0.125d0*(1d0+c2b)/c2t)+0.125d0*(mt*s2b)/(mb*s2t)+
     .         0.125d0*(s2b*Xt)/(c2t*mb))))
      Dbc2t = Dbc2t-tmp88*
     .    (cbe*hb*ht*sbe*((mt*tmp55)/mb-0.25d0*s2b/s2t-
     .       0.25d0*tmp155/(mb*s2t))+
     .      hb**2*sbe**2*(0.125d0*(1d0+c2b)/c2t-
     .       0.125d0*(mt*s2b)/(mb*s2t)+0.125d0*(s2b*Yb)/(c2t*mb))+
     .      cbe**2*ht**2*
     .       (-(0.125d0*(1d0-c2b)/c2t)-0.125d0*(mt*s2b)/(mb*s2t)-
     .       0.125d0*(s2b*Yt)/(c2t*mb)))-
     .   tmp109*(cbe*hb*ht*sbe*
     .       ((mt*tmp55)/mb-0.25d0*s2b/s2t+0.25d0*tmp154/(mb*s2t))+
     .      hb**2*sbe**2*
     .       (-(0.125d0*(1d0-c2b)/c2t)-0.125d0*(mt*s2b)/(mb*s2t)+
     .       0.125d0*(s2b*Yb)/(c2t*mb))+
     .      cbe**2*ht**2*(0.125d0*(1d0+c2b)/c2t-
     .       0.125d0*(mt*s2b)/(mb*s2t)-0.125d0*(s2b*Yt)/(c2t*mb)))-
     .     tmp108*(cbe*hb*ht*sbe*
     .       ((mt*tmp54)/mb+0.25d0*s2b/s2t+0.25d0*tmp155/(mb*s2t))+
     .      hb**2*sbe**2*
     .       (-(0.125d0*(1d0+c2b)/c2t)+0.125d0*(mt*s2b)/(mb*s2t)-
     .       0.125d0*(s2b*Yb)/(c2t*mb))+
     .      cbe**2*ht**2*(0.125d0*(1d0-c2b)/c2t+
     .       0.125d0*(mt*s2b)/(mb*s2t)+0.125d0*(s2b*Yt)/(c2t*mb)))-
     .     tmp89*(cbe*hb*ht*sbe*
     .       ((mt*tmp54)/mb+0.25d0*s2b/s2t-0.25d0*tmp154/(mb*s2t))+
     .      hb**2*sbe**2*
     .       (0.125d0*(1d0-c2b)/c2t+0.125d0*(mt*s2b)/(mb*s2t)-
     .       0.125d0*(s2b*Yb)/(c2t*mb))+
     .      cbe**2*ht**2*(-(0.125d0*(1d0+c2b)/c2t)+
     .       0.125d0*(mt*s2b)/(mb*s2t)+0.125d0*(s2b*Yt)/(c2t*mb)))
      DB1c2t = hb**2*(0.125d0*((1d0-c2b)*(-1d0+LogT1)*T1)/c2t+
     .      0.125d0*((1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT1)*T1)/c2t-
     .      0.125d0*((1d0-c2b)*(-1d0+LogT2)*T2)/c2t-
     .      0.125d0*((1d0-c2b)*(-1d0+LogB1)*(-1d0+LogT2)*T2)/c2t)+
     .   ht**2*(-(0.125d0*((1d0+c2b)*(-1d0+LogT1)*T1)/c2t)-
     .      0.125d0*((1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT1)*T1)/c2t+
     .      0.125d0*((1d0+c2b)*(-1d0+LogT2)*T2)/c2t+
     .      0.125d0*((1d0+c2b)*(-1d0+LogB1)*(-1d0+LogT2)*T2)/c2t)-
     .   tmp206*(-0.5d0+2d0*LogB1-(phiA0B1T2*(-A0+B1-T2))/T2-
     .      0.5d0*(LogA0*LogB1)+0.5d0*(LogA0*LogT2)-
     .      0.5d0*(LogB1*LogT2)+0.5d0*(LogT2*(A0-B1-T2))/B1+
     .      0.5d0*(deltA0B1T2*phiA0B1T2)/(B1*T2)+
     .      0.5d0*(LogA0*(-A0-B1+T2))/B1-
     .      0.5d0*(deltA0B1T2*
     .        (tmp105/deltA0B1T2+
     .          (B1*phiA0B1T2*tmp180)/(deltA0B1T2*T2)))/B1)-
     .   0.5d0*(tmp200*tmp211)-0.5d0*(tmp202*tmp215)-
     .   tmp204*(-0.5d0+2d0*LogB1-(phiA0B1T1*(-A0+B1-T1))/T1-
     .      0.5d0*(LogA0*LogB1)+0.5d0*(LogA0*LogT1)-
     .      0.5d0*(LogB1*LogT1)+0.5d0*(LogT1*(A0-B1-T1))/B1+
     .      0.5d0*(deltA0B1T1*phiA0B1T1)/(B1*T1)+
     .      0.5d0*(LogA0*(-A0-B1+T1))/B1-
     .      0.5d0*(deltA0B1T1*
     .        ((B1*phiA0B1T1*tmp176)/(deltA0B1T1*T1)+
     .          tmp85/deltA0B1T1))/B1)
      DB2c2t = ht**2*(-(0.125d0*
     .       ((1d0-c2b)*(-1d0+LogT1)*T1)/c2t)-
     .      0.125d0*((1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT1)*T1)/c2t+
     .      0.125d0*((1d0-c2b)*(-1d0+LogT2)*T2)/c2t+
     .      0.125d0*((1d0-c2b)*(-1d0+LogB2)*(-1d0+LogT2)*T2)/c2t)+
     .   hb**2*(0.125d0*((1d0+c2b)*(-1d0+LogT1)*T1)/c2t+
     .      0.125d0*((1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT1)*T1)/c2t-
     .      0.125d0*((1d0+c2b)*(-1d0+LogT2)*T2)/c2t-
     .      0.125d0*((1d0+c2b)*(-1d0+LogB2)*(-1d0+LogT2)*T2)/c2t)-
     .   tmp205*(-0.5d0+2d0*LogB2-(phiA0B2T2*(-A0+B2-T2))/T2-
     .      0.5d0*(LogA0*LogB2)+0.5d0*(LogA0*LogT2)-
     .      0.5d0*(LogB2*LogT2)+0.5d0*(LogT2*(A0-B2-T2))/B2+
     .      0.5d0*(deltA0B2T2*phiA0B2T2)/(B2*T2)+
     .      0.5d0*(LogA0*(-A0-B2+T2))/B2-
     .      0.5d0*(deltA0B2T2*
     .        (tmp106/deltA0B2T2+
     .          (B2*phiA0B2T2*tmp181)/(deltA0B2T2*T2)))/B2)-
     .   0.5d0*(tmp199*tmp212)-0.5d0*(tmp201*tmp216)-
     .   tmp203*(-0.5d0+2d0*LogB2-(phiA0B2T1*(-A0+B2-T1))/T1-
     .      0.5d0*(LogA0*LogB2)+0.5d0*(LogA0*LogT1)-
     .      0.5d0*(LogB2*LogT1)+0.5d0*(LogT1*(A0-B2-T1))/B2+
     .      0.5d0*(deltA0B2T1*phiA0B2T1)/(B2*T1)+
     .      0.5d0*(LogA0*(-A0-B2+T1))/B2-
     .      0.5d0*(deltA0B2T1*
     .        ((B2*phiA0B2T1*tmp177)/(deltA0B2T1*T1)+
     .          tmp86/deltA0B2T1))/B2)
      DT1c2b = ht**2*(0.125d0*(B1*(1d0-c2t)*(-1d0+LogB1))/c2b-
     .      0.125d0*(B2*(1d0-c2t)*(-1d0+LogB2))/c2b+
     .      0.125d0*(B1*(1d0-c2t)*(-1d0+LogB1)*(-1d0+LogT1))/c2b-
     .      0.125d0*(B2*(1d0-c2t)*(-1d0+LogB2)*(-1d0+LogT1))/c2b)+
     .   hb**2*(-(0.125d0*(B1*(1d0+c2t)*(-1d0+LogB1))/c2b)+
     .      0.125d0*(B2*(1d0+c2t)*(-1d0+LogB2))/c2b-
     .      0.125d0*(B1*(1d0+c2t)*(-1d0+LogB1)*(-1d0+LogT1))/c2b+
     .      0.125d0*(B2*(1d0+c2t)*(-1d0+LogB2)*(-1d0+LogT1))/c2b)-
     .   0.5d0*(tmp214*(cbe**2*hb**2*
     .       (-(0.125d0*(b*(1d0+c2t))/c2b)+
     .         0.25d0*(mb*mt*s2t)/s2b+0.125d0*((1d0-c2t)*t)/c2b+
     .         0.25d0*((1d0+c2t)*mb*Xb)/s2b+
     .         0.25d0*(mt*s2t*Xb)/c2b+0.125d0*((1d0+c2t)*Xb**2)/c2b
     .         )-cbe*hb*ht*sbe*
     .       (mb*s2t*tmp23+2d0*mb*mt*tmp57+
     .         0.25d0*(s2t*(b+t))/s2b+0.5d0*(mt*tmp130)/s2b+
     .         0.25d0*(s2t*Xb*Xt)/s2b)+
     .      ht**2*sbe**2*
     .       (0.125d0*(b*(1d0-c2t))/c2b+0.25d0*(mb*mt*s2t)/s2b-
     .         0.125d0*((1d0+c2t)*t)/c2b+
     .         0.25d0*((1d0-c2t)*mb*Xt)/s2b-
     .         0.25d0*(mt*s2t*Xt)/c2b-0.125d0*((1d0-c2t)*Xt**2)/c2b
     .         )))-0.5d0*
     .    (tmp213*(cbe**2*hb**2*
     .       (0.125d0*(b*(1d0+c2t))/c2b-0.25d0*(mb*mt*s2t)/s2b-
     .         0.125d0*((1d0-c2t)*t)/c2b-
     .         0.25d0*((1d0+c2t)*mb*Xb)/s2b-
     .         0.25d0*(mt*s2t*Xb)/c2b-0.125d0*((1d0+c2t)*Xb**2)/c2b
     .         )-cbe*hb*ht*sbe*
     .       (mb*s2t*tmp24+2d0*mb*mt*tmp58-
     .         0.25d0*(s2t*(b+t))/s2b-0.5d0*(mt*tmp130)/s2b-
     .         0.25d0*(s2t*Xb*Xt)/s2b)+
     .      ht**2*sbe**2*
     .       (-(0.125d0*(b*(1d0-c2t))/c2b)-
     .         0.25d0*(mb*mt*s2t)/s2b+0.125d0*((1d0+c2t)*t)/c2b-
     .         0.25d0*((1d0-c2t)*mb*Xt)/s2b+
     .         0.25d0*(mt*s2t*Xt)/c2b+0.125d0*((1d0-c2t)*Xt**2)/c2b
     .         )))
      DT1c2b = DT1c2b-
     .   tmp239*(hb**2*sbe**2*
     .       (-(0.125d0*(b*(1d0+c2t))/c2b)+0.25d0*(mb*mt*s2t)/s2b+
     .       0.125d0*((1d0-c2t)*t)/c2b+
     .       0.25d0*((1d0+c2t)*mb*Yb)/s2b+0.25d0*(mt*s2t*Yb)/c2b+
     .       0.125d0*((1d0+c2t)*Yb**2)/c2b)+
     .      cbe*hb*ht*sbe*(mb*s2t*tmp39+2d0*mb*mt*tmp57+
     .       0.25d0*(s2t*(b+t))/s2b+0.5d0*(mt*tmp156)/s2b+
     .       0.25d0*(s2t*Yb*Yt)/s2b)+
     .      cbe**2*ht**2*(0.125d0*(b*(1d0-c2t))/c2b+
     .       0.25d0*(mb*mt*s2t)/s2b-0.125d0*((1d0+c2t)*t)/c2b+
     .       0.25d0*((1d0-c2t)*mb*Yt)/s2b-0.25d0*(mt*s2t*Yt)/c2b-
     .       0.125d0*((1d0-c2t)*Yt**2)/c2b))-
     .   tmp238*(hb**2*sbe**2*
     .       (0.125d0*(b*(1d0+c2t))/c2b-0.25d0*(mb*mt*s2t)/s2b-
     .       0.125d0*((1d0-c2t)*t)/c2b-
     .       0.25d0*((1d0+c2t)*mb*Yb)/s2b-0.25d0*(mt*s2t*Yb)/c2b-
     .       0.125d0*((1d0+c2t)*Yb**2)/c2b)+
     .      cbe*hb*ht*sbe*(mb*s2t*tmp40+2d0*mb*mt*tmp58-
     .       0.25d0*(s2t*(b+t))/s2b-0.5d0*(mt*tmp156)/s2b-
     .       0.25d0*(s2t*Yb*Yt)/s2b)+
     .      cbe**2*ht**2*(-(0.125d0*(b*(1d0-c2t))/c2b)-
     .       0.25d0*(mb*mt*s2t)/s2b+0.125d0*((1d0+c2t)*t)/c2b-
     .       0.25d0*((1d0-c2t)*mb*Yt)/s2b+0.25d0*(mt*s2t*Yt)/c2b+
     .       0.125d0*((1d0-c2t)*Yt**2)/c2b))
      DT2c2b = hb**2*(-(0.125d0*
     .       (B1*(1d0-c2t)*(-1d0+LogB1))/c2b)+
     .      0.125d0*(B2*(1d0-c2t)*(-1d0+LogB2))/c2b-
     .      0.125d0*(B1*(1d0-c2t)*(-1d0+LogB1)*(-1d0+LogT2))/c2b+
     .      0.125d0*(B2*(1d0-c2t)*(-1d0+LogB2)*(-1d0+LogT2))/c2b)+
     .   ht**2*(0.125d0*(B1*(1d0+c2t)*(-1d0+LogB1))/c2b-
     .      0.125d0*(B2*(1d0+c2t)*(-1d0+LogB2))/c2b+
     .      0.125d0*(B1*(1d0+c2t)*(-1d0+LogB1)*(-1d0+LogT2))/c2b-
     .      0.125d0*(B2*(1d0+c2t)*(-1d0+LogB2)*(-1d0+LogT2))/c2b)-
     .   0.5d0*(tmp218*(cbe**2*hb**2*
     .       (-(0.125d0*(b*(1d0-c2t))/c2b)-
     .         0.25d0*(mb*mt*s2t)/s2b+0.125d0*((1d0+c2t)*t)/c2b+
     .         0.25d0*((1d0-c2t)*mb*Xb)/s2b-
     .         0.25d0*(mt*s2t*Xb)/c2b+0.125d0*((1d0-c2t)*Xb**2)/c2b
     .         )-cbe*hb*ht*sbe*
     .       (-(mb*s2t*tmp23)+2d0*mb*mt*tmp58-
     .         0.25d0*(s2t*(b+t))/s2b+0.5d0*(mt*tmp131)/s2b-
     .         0.25d0*(s2t*Xb*Xt)/s2b)+
     .      ht**2*sbe**2*
     .       (0.125d0*(b*(1d0+c2t))/c2b-0.25d0*(mb*mt*s2t)/s2b-
     .         0.125d0*((1d0-c2t)*t)/c2b+
     .         0.25d0*((1d0+c2t)*mb*Xt)/s2b+
     .         0.25d0*(mt*s2t*Xt)/c2b-0.125d0*((1d0+c2t)*Xt**2)/c2b
     .         )))-0.5d0*
     .    (tmp217*(cbe**2*hb**2*
     .       (0.125d0*(b*(1d0-c2t))/c2b+0.25d0*(mb*mt*s2t)/s2b-
     .         0.125d0*((1d0+c2t)*t)/c2b-
     .         0.25d0*((1d0-c2t)*mb*Xb)/s2b+
     .         0.25d0*(mt*s2t*Xb)/c2b-0.125d0*((1d0-c2t)*Xb**2)/c2b
     .         )-cbe*hb*ht*sbe*
     .       (-(mb*s2t*tmp24)+2d0*mb*mt*tmp57+
     .         0.25d0*(s2t*(b+t))/s2b-0.5d0*(mt*tmp131)/s2b+
     .         0.25d0*(s2t*Xb*Xt)/s2b)+
     .      ht**2*sbe**2*
     .       (-(0.125d0*(b*(1d0+c2t))/c2b)+
     .         0.25d0*(mb*mt*s2t)/s2b+0.125d0*((1d0-c2t)*t)/c2b-
     .         0.25d0*((1d0+c2t)*mb*Xt)/s2b-
     .         0.25d0*(mt*s2t*Xt)/c2b+0.125d0*((1d0+c2t)*Xt**2)/c2b
     .         )))
      DT2c2b = DT2c2b-
     .   tmp223*(hb**2*sbe**2*
     .       (-(0.125d0*(b*(1d0-c2t))/c2b)-0.25d0*(mb*mt*s2t)/s2b+
     .       0.125d0*((1d0+c2t)*t)/c2b+
     .       0.25d0*((1d0-c2t)*mb*Yb)/s2b-0.25d0*(mt*s2t*Yb)/c2b+
     .       0.125d0*((1d0-c2t)*Yb**2)/c2b)+
     .      cbe*hb*ht*sbe*(-(mb*s2t*tmp39)+2d0*mb*mt*tmp58-
     .       0.25d0*(s2t*(b+t))/s2b+0.5d0*(mt*tmp157)/s2b-
     .       0.25d0*(s2t*Yb*Yt)/s2b)+
     .      cbe**2*ht**2*(0.125d0*(b*(1d0+c2t))/c2b-
     .       0.25d0*(mb*mt*s2t)/s2b-0.125d0*((1d0-c2t)*t)/c2b+
     .       0.25d0*((1d0+c2t)*mb*Yt)/s2b+0.25d0*(mt*s2t*Yt)/c2b-
     .       0.125d0*((1d0+c2t)*Yt**2)/c2b))-
     .   tmp222*(hb**2*sbe**2*
     .       (0.125d0*(b*(1d0-c2t))/c2b+0.25d0*(mb*mt*s2t)/s2b-
     .       0.125d0*((1d0+c2t)*t)/c2b-
     .       0.25d0*((1d0-c2t)*mb*Yb)/s2b+0.25d0*(mt*s2t*Yb)/c2b-
     .       0.125d0*((1d0-c2t)*Yb**2)/c2b)+
     .      cbe*hb*ht*sbe*(-(mb*s2t*tmp40)+2d0*mb*mt*tmp57+
     .       0.25d0*(s2t*(b+t))/s2b-0.5d0*(mt*tmp157)/s2b+
     .       0.25d0*(s2t*Yb*Yt)/s2b)+
     .      cbe**2*ht**2*(-(0.125d0*(b*(1d0+c2t))/c2b)+
     .       0.25d0*(mb*mt*s2t)/s2b+0.125d0*((1d0-c2t)*t)/c2b-
     .       0.25d0*((1d0+c2t)*mb*Yt)/s2b-0.25d0*(mt*s2t*Yt)/c2b+
     .       0.125d0*((1d0+c2t)*Yt**2)/c2b))
      Dc2tc2b = hb**2*tmp19+ht**2*tmp19-
     .   tmp89*(hb**2*sbe**2*tmp38+cbe**2*ht**2*tmp49+
     .      cbe*hb*ht*sbe*(-(0.25d0*(mb*mt)/(c2b*c2t))-
     .       0.125d0*(b+t)/(s2b*s2t)-0.5d0*(mb*tmp39)/s2t+
     .       0.5d0*(mt*tmp43)/s2b-0.125d0*(Yb*Yt)/(s2b*s2t)))-
     .   tmp108*(hb**2*sbe**2*tmp38+cbe**2*ht**2*tmp49+
     .      cbe*hb*ht*sbe*(-(0.25d0*(mb*mt)/(c2b*c2t))-
     .       0.125d0*(b+t)/(s2b*s2t)+0.5d0*(mb*tmp40)/s2t-
     .       0.5d0*(mt*tmp44)/s2b-0.125d0*(Yb*Yt)/(s2b*s2t)))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (cbe**2*hb**2*tmp22+ht**2*sbe**2*tmp33-
     .      cbe*hb*ht*sbe*
     .       (-(0.25d0*(mb*mt)/(c2b*c2t))-
     .         0.125d0*(b+t)/(s2b*s2t)-0.5d0*(mb*tmp23)/s2t+
     .         0.5d0*(mt*tmp27)/s2b-0.125d0*(Xb*Xt)/(s2b*s2t))))-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (cbe**2*hb**2*tmp22+ht**2*sbe**2*tmp33-
     .      cbe*hb*ht*sbe*
     .       (-(0.25d0*(mb*mt)/(c2b*c2t))-
     .         0.125d0*(b+t)/(s2b*s2t)+0.5d0*(mb*tmp24)/s2t-
     .         0.5d0*(mt*tmp28)/s2b-0.125d0*(Xb*Xt)/(s2b*s2t))))-
     .   0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (cbe**2*hb**2*tmp21+ht**2*sbe**2*tmp32-
     .      cbe*hb*ht*sbe*
     .       (0.25d0*(mb*mt)/(c2b*c2t)+0.125d0*(b+t)/(s2b*s2t)-
     .         0.5d0*(mb*tmp24)/s2t-0.5d0*(mt*tmp27)/s2b+
     .         0.125d0*(Xb*Xt)/(s2b*s2t))))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (cbe**2*hb**2*tmp21+ht**2*sbe**2*tmp32-
     .      cbe*hb*ht*sbe*
     .       (0.25d0*(mb*mt)/(c2b*c2t)+0.125d0*(b+t)/(s2b*s2t)+
     .         0.5d0*(mb*tmp23)/s2t+0.5d0*(mt*tmp28)/s2b+
     .         0.125d0*(Xb*Xt)/(s2b*s2t))))
      Dc2tc2b = Dc2tc2b-
     .   tmp88*(hb**2*sbe**2*tmp37+cbe**2*ht**2*tmp48+
     .      cbe*hb*ht*sbe*(0.25d0*(mb*mt)/(c2b*c2t)+
     .       0.125d0*(b+t)/(s2b*s2t)-0.5d0*(mb*tmp40)/s2t-
     .       0.5d0*(mt*tmp43)/s2b+0.125d0*(Yb*Yt)/(s2b*s2t)))-
     .   tmp109*(hb**2*sbe**2*tmp37+cbe**2*ht**2*tmp48+
     .      cbe*hb*ht*sbe*(0.25d0*(mb*mt)/(c2b*c2t)+
     .       0.125d0*(b+t)/(s2b*s2t)+0.5d0*(mb*tmp39)/s2t+
     .       0.5d0*(mt*tmp44)/s2b+0.125d0*(Yb*Yt)/(s2b*s2t)))
      Dcptpb = -2d0*cbe*hb*ht*mb*mt*sbe*tmp108*tmp56+
     .   cbe*hb*ht*mb*mt*sbe*
     .    (-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*tmp56
     .    +cbe*hb*ht*mb*mt*sbe*
     .    (-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*tmp56
     .    -2d0*cbe*hb*ht*mb*mt*sbe*tmp109*tmp59+
     .   cbe*hb*ht*mb*mt*sbe*
     .    (-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*tmp59
     .    +cbe*hb*ht*mb*mt*sbe*
     .    (-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*tmp59
     .    -2d0*cbe*hb*ht*mb*mt*sbe*tmp59*tmp88-
     .   2d0*cbe*hb*ht*mb*mt*sbe*tmp56*tmp89-
     .   4d0*cbe*hb*ht*mb*mt*sbe*
     .    (-2d0*A0*LogA0-2d0*b*Logb-2d0*Logt*t-
     .      0.5d0*(Logb*Logt*(A0-b-t))-
     .      0.5d0*(LogA0*Logt*(-A0+b-t))+
     .      0.5d0*(deltA0bt*phiA0bt)/t-
     .      0.5d0*(LogA0*Logb*(-A0-b+t))+2.5d0*(A0+b+t)+
     .      0.5d0*tmp76)
      Dcpttptb = 0.25d0*(cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*Xb*
     .      Xt)-0.25d0*(cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*Xb*
     .      Xt)-0.25d0*(cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*Xb*
     .      Xt)+0.25d0*(cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*Xb*
     .      Xt)+0.5d0*(cbe*hb*ht*s2b*s2t*sbe*tmp108*Yb*Yt)-
     .   0.5d0*(cbe*hb*ht*s2b*s2t*sbe*tmp109*Yb*Yt)-
     .   0.5d0*(cbe*hb*ht*s2b*s2t*sbe*tmp88*Yb*Yt)+
     .   0.5d0*(cbe*hb*ht*s2b*s2t*sbe*tmp89*Yb*Yt)
      Dcpbptt = -2d0*hb*ht*mb*mu*s2t*tmp113-
     .   cbe*hb*ht*sbe*tmp89*(mb*s2t*tmp154-0.5d0*(b*s2b*s2t))-
     .   cbe*hb*ht*sbe*tmp108*
     .    (-(mb*s2t*tmp155)-0.5d0*(b*s2b*s2t))-
     .   cbe*hb*ht*sbe*tmp109*
     .    (-(mb*s2t*tmp154)+0.5d0*(b*s2b*s2t))-
     .   cbe*hb*ht*sbe*tmp88*(mb*s2t*tmp155+0.5d0*(b*s2b*s2t))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (mb*s2t*tmp128-0.5d0*(b*s2b*s2t)))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (-(mb*s2t*tmp129)-0.5d0*(b*s2b*s2t)))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (-(mb*s2t*tmp128)+0.5d0*(b*s2b*s2t)))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (mb*s2t*tmp129+0.5d0*(b*s2b*s2t)))
      Dcptptb = -2d0*hb*ht*mt*mu*s2b*tmp75-
     .   cbe*hb*ht*sbe*tmp89*
     .    (-(mt*s2b*tmp156)-0.5d0*(s2b*s2t*t))-
     .   cbe*hb*ht*sbe*tmp108*(mt*s2b*tmp157-0.5d0*(s2b*s2t*t))-
     .   cbe*hb*ht*sbe*tmp88*(mt*s2b*tmp156+0.5d0*(s2b*s2t*t))-
     .   cbe*hb*ht*sbe*tmp109*
     .    (-(mt*s2b*tmp157)+0.5d0*(s2b*s2t*t))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (-(mt*s2b*tmp130)-0.5d0*(s2b*s2t*t)))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (mt*s2b*tmp131-0.5d0*(s2b*s2t*t)))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169)*
     .      (mt*s2b*tmp130+0.5d0*(s2b*s2t*t)))+
     .   0.5d0*(cbe*hb*ht*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (-(mt*s2b*tmp131)+0.5d0*(s2b*s2t*t)))
      Dcptmptt = -(tmp109*
     .      (0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t)+
     .      hb**2*sbe**2*
     .       (0.5d0*(mb*mt*s2b*s2t)-0.5d0*((1d0+c2b)*mt*s2t*Yb))+
     .      cbe**2*ht**2*
     .       (0.5d0*(mb*mt*s2b*s2t)-0.5d0*((1d0-c2b)*mt*s2t*Yt))))-
     .     tmp89*(-(0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t))+
     .      hb**2*sbe**2*(-(0.5d0*(mb*mt*s2b*s2t))+
     .       0.5d0*((1d0+c2b)*mt*s2t*Yb))+
     .      cbe**2*ht**2*(-(0.5d0*(mb*mt*s2b*s2t))+
     .       0.5d0*((1d0-c2b)*mt*s2t*Yt)))-
     .   tmp108*(-(0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t))+
     .      hb**2*sbe**2*(-(0.5d0*(mb*mt*s2b*s2t))-
     .       0.5d0*((1d0-c2b)*mt*s2t*Yb))+
     .      cbe**2*ht**2*(-(0.5d0*(mb*mt*s2b*s2t))-
     .       0.5d0*((1d0+c2b)*mt*s2t*Yt)))-
     .   0.25d0*(ht**2*(-2d0*mt*s2t*sbe**2*tmp115*Xt+
     .      2d0*mt*s2t*sbe**2*tmp94*Xt-
     .      4d0*cbe**2*mt*s2t*tmp114*Yt+4d0*cbe**2*mt*s2t*tmp93*Yt)
     .      )-0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-
     .      5d0*T2+LogT2*(-2d0*B2*LogB2+4d0*T2)-
     .      2d0*(-B2+T2)*tmp173)*
     .      (-(0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t))+
     .      cbe**2*hb**2*
     .       (0.5d0*(mb*mt*s2b*s2t)-0.5d0*((1d0+c2b)*mt*s2t*Xb))+
     .      ht**2*sbe**2*
     .       (0.5d0*(mb*mt*s2b*s2t)-0.5d0*((1d0-c2b)*mt*s2t*Xt))))-
     .     0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t)+
     .      cbe**2*hb**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))+0.5d0*((1d0+c2b)*mt*s2t*Xb))
     .       +ht**2*sbe**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))+0.5d0*((1d0-c2b)*mt*s2t*Xt))
     .      ))-0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-
     .      5d0*T2+LogT2*(-2d0*B1*LogB1+4d0*T2)-
     .      2d0*(-B1+T2)*tmp172)*
     .      (0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t)+
     .      cbe**2*hb**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))-0.5d0*((1d0-c2b)*mt*s2t*Xb))
     .       +ht**2*sbe**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))-0.5d0*((1d0+c2b)*mt*s2t*Xt))
     .      ))-0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-
     .      5d0*T1+LogT1*(-2d0*B1*LogB1+4d0*T1)-
     .      2d0*(-B1+T1)*tmp169)*
     .      (-(0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t))+
     .      cbe**2*hb**2*
     .       (0.5d0*(mb*mt*s2b*s2t)+0.5d0*((1d0-c2b)*mt*s2t*Xb))+
     .      ht**2*sbe**2*
     .       (0.5d0*(mb*mt*s2b*s2t)+0.5d0*((1d0+c2b)*mt*s2t*Xt))))
      Dcptmptt = Dcptmptt-
     .   tmp88*(0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t)+
     .      hb**2*sbe**2*(0.5d0*(mb*mt*s2b*s2t)+
     .       0.5d0*((1d0-c2b)*mt*s2t*Yb))+
     .      cbe**2*ht**2*(0.5d0*(mb*mt*s2b*s2t)+
     .       0.5d0*((1d0+c2b)*mt*s2t*Yt)))
      Dcpbmptb = -(tmp89*
     .      (-(0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe))+
     .      hb**2*sbe**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))-0.5d0*((1d0+c2t)*mb*s2b*Yb))
     .       +cbe**2*ht**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))-0.5d0*((1d0-c2t)*mb*s2b*Yt))
     .      ))-0.25d0*(hb**2*
     .      (2d0*cbe**2*(-10d0*B1+4d0*B1*LogB1+
     .         LogB1*(4d0*B1-2d0*B1*LogB1))*mb*s2b*Xb-
     .      2d0*cbe**2*(-10d0*B2+4d0*B2*LogB2+
     .         LogB2*(4d0*B2-2d0*B2*LogB2))*mb*s2b*Xb+
     .      4d0*mb*s2b*sbe**2*Yb*
     .       (2d0*A0*LogA0+4d0*B1*LogB1-A0*LogA0*LogB1-
     .         2.5d0*(A0+2d0*B1)+0.5d0*((A0-2d0*B1)*LogB1**2)-
     .         0.5d0*(deltA0B1B1*phiA0B1B1)/B1)-
     .      4d0*mb*s2b*sbe**2*Yb*
     .       (2d0*A0*LogA0+4d0*B2*LogB2-A0*LogA0*LogB2-
     .         2.5d0*(A0+2d0*B2)+0.5d0*((A0-2d0*B2)*LogB2**2)-
     .         0.5d0*(deltA0B2B2*phiA0B2B2)/B2)))-
     .   0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170)*
     .      (0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe)+
     .      cbe**2*hb**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))-0.5d0*((1d0+c2t)*mb*s2b*Xb))
     .       +ht**2*sbe**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))-0.5d0*((1d0-c2t)*mb*s2b*Xt))
     .      ))-0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-
     .      5d0*T1+LogT1*(-2d0*B1*LogB1+4d0*T1)-
     .      2d0*(-B1+T1)*tmp169)*
     .      (-(0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe))+
     .      cbe**2*hb**2*
     .       (0.5d0*(mb*mt*s2b*s2t)+0.5d0*((1d0+c2t)*mb*s2b*Xb))+
     .      ht**2*sbe**2*
     .       (0.5d0*(mb*mt*s2b*s2t)+0.5d0*((1d0-c2t)*mb*s2b*Xt))))-
     .     0.5d0*((-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173)*
     .      (-(0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe))+
     .      cbe**2*hb**2*
     .       (0.5d0*(mb*mt*s2b*s2t)-0.5d0*((1d0-c2t)*mb*s2b*Xb))+
     .      ht**2*sbe**2*
     .       (0.5d0*(mb*mt*s2b*s2t)-0.5d0*((1d0+c2t)*mb*s2b*Xt))))-
     .     0.5d0*((-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172)*
     .      (0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe)+
     .      cbe**2*hb**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))+0.5d0*((1d0-c2t)*mb*s2b*Xb))
     .       +ht**2*sbe**2*
     .       (-(0.5d0*(mb*mt*s2b*s2t))+0.5d0*((1d0+c2t)*mb*s2b*Xt))
     .      ))
      Dcpbmptb = Dcpbmptb-
     .   tmp88*(0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe)+
     .      hb**2*sbe**2*(0.5d0*(mb*mt*s2b*s2t)+
     .       0.5d0*((1d0+c2t)*mb*s2b*Yb))+
     .      cbe**2*ht**2*(0.5d0*(mb*mt*s2b*s2t)+
     .       0.5d0*((1d0-c2t)*mb*s2b*Yt)))-
     .   tmp109*(0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe)+
     .      hb**2*sbe**2*(0.5d0*(mb*mt*s2b*s2t)-
     .       0.5d0*((1d0-c2t)*mb*s2b*Yb))+
     .      cbe**2*ht**2*(0.5d0*(mb*mt*s2b*s2t)-
     .       0.5d0*((1d0+c2t)*mb*s2b*Yt)))-
     .   tmp108*(-(0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe))+
     .      hb**2*sbe**2*(-(0.5d0*(mb*mt*s2b*s2t))+
     .       0.5d0*((1d0-c2t)*mb*s2b*Yb))+
     .      cbe**2*ht**2*(-(0.5d0*(mb*mt*s2b*s2t))+
     .       0.5d0*((1d0+c2t)*mb*s2b*Yt)))
      Dspbmptbspbptt =
     .  -(0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe*tmp108))+
     .   0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe*tmp109)-
     .   0.25d0*(b*cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169))+
     .   0.25d0*(b*cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170))+
     .   0.25d0*(b*cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172))-
     .   0.25d0*(b*cbe*hb*ht*s2b*s2t*sbe*
     .      (-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173))+
     .   0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe*tmp88)-
     .   0.5d0*(b*cbe*hb*ht*s2b*s2t*sbe*tmp89)
      Dsptmpttsptptb =
     .  -(0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t*tmp108))+
     .   0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t*tmp109)-
     .   0.25d0*(cbe*hb*ht*s2b*s2t*sbe*t*
     .      (-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169))+
     .   0.25d0*(cbe*hb*ht*s2b*s2t*sbe*t*
     .      (-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-5d0*T1+
     .      LogT1*(-2d0*B2*LogB2+4d0*T1)-2d0*(-B2+T1)*tmp170))+
     .   0.25d0*(cbe*hb*ht*s2b*s2t*sbe*t*
     .      (-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-5d0*T2+
     .      LogT2*(-2d0*B1*LogB1+4d0*T2)-2d0*(-B1+T2)*tmp172))-
     .   0.25d0*(cbe*hb*ht*s2b*s2t*sbe*t*
     .      (-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-5d0*T2+
     .      LogT2*(-2d0*B2*LogB2+4d0*T2)-2d0*(-B2+T2)*tmp173))+
     .   0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t*tmp88)-
     .   0.5d0*(cbe*hb*ht*s2b*s2t*sbe*t*tmp89)
      Dsptmpttspbmptb =
     .  -(tmp108*tmp11)-tmp109*tmp12-tmp12*tmp88-
     .   tmp11*tmp89-0.5d0*
     .    (tmp14*(-5d0*B1+4d0*B1*LogB1+LogT1**2*(B1-T1)-5d0*T1+
     .      LogT1*(-2d0*B1*LogB1+4d0*T1)-2d0*(-B1+T1)*tmp169))-
     .   0.5d0*(tmp13*(-5d0*B2+4d0*B2*LogB2+LogT1**2*(B2-T1)-
     .      5d0*T1+LogT1*(-2d0*B2*LogB2+4d0*T1)-
     .      2d0*(-B2+T1)*tmp170))-
     .   0.5d0*(tmp13*(-5d0*B1+4d0*B1*LogB1+LogT2**2*(B1-T2)-
     .      5d0*T2+LogT2*(-2d0*B1*LogB1+4d0*T2)-
     .      2d0*(-B1+T2)*tmp172))-
     .   0.5d0*(tmp14*(-5d0*B2+4d0*B2*LogB2+LogT2**2*(B2-T2)-
     .      5d0*T2+LogT2*(-2d0*B2*LogB2+4d0*T2)-
     .      2d0*(-B2+T2)*tmp173))

      end

*      
***********************************************************************
*

      SUBROUTINE getPiHpHm(g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
     $     Al,At,Ab,Atau,p,Q,piHpHm)
* by courtesy of P. Slavich and K.H. Phan

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
     $     Al,At,Ab,Atau,Q,p,piHpHm
      
      DOUBLE PRECISION mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     $     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),
     $     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     $     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2),
     $     myB0,myF,myA0,myG,pi,cb,sb,mw2,mz2,mt2,mb2,mtau2,
     $     gb2,sw2,cw2,aHpnech(5,2),bHpnech(5,2),
     $     lHpHmtt(2,2),lHpHmbb(2,2),lHpHmtata(2,2),lHpHmntnt,
     $     lHpHmuu(2,2),lHpHmdd(2,2),lHpHmee(2,2),lHpHmnn,
     $     lHptb(2,2),lHpntta(2),lHpud(2,2),lHpne(2),  
     $     lHpHmhh(3,3),lHpHmaa(3,3),lHphc(3,2),lHpac(3,2),lHpHmcc(2),
     $     higgs,weak,fermions,sfermions,sfermions3g,inos

      INTEGER i,j
      
      COMMON/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     $     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue
      
      pi = 4d0*atan(1d0)

      PiHpHm = 0d0

      cb = v1/sqrt(v1**2+v2**2)
      sb = v2/sqrt(v1**2+v2**2)
      gb2 = (g**2+gp**2)/2d0
      sw2 = gp**2/2d0/gb2
      cw2 = 1-sw2
      mw2 = g**2/2d0*(v1**2+v2**2)
      mz2 = gb2*(v1**2+v2**2)

*     the higgs+gauge and pure gauge contributions

      weak = 0d0

      do i = 1,3                

         weak = weak            ! scalar-W
     $        + g**2/4d0*(sb*RS(i,1)-cb*RS(i,2))**2
     $        *myF(p**2,mhh(i)**2,mw2,Q**2)
         
         weak = weak            ! pseudoscalar-W
     $        + g**2/4d0*(sb*RP(i,1)+cb*RP(i,2))**2
     $        *myF(p**2,maa(i)**2,mw2,Q**2)
         
      enddo

      weak = weak               ! pure gauge
     $     +g**2/4d0*(cw2-sw2)**2/cw2*myF(p**2,mhc(2)**2,mz2,Q**2)
     $     +sw2*g**2*myF(p**2,mhc(2)**2,0d0,Q**2)
     $     +2*g**2*myA0(mw2,Q**2)
     $     +g**2*(cw2-sw2)**2/cw2*myA0(mz2,Q**2)
      
*     the pure higgs contributions

      call coupl_Hp_hh(g,gp,ll,kk,v1,v2,xx,Al,RS,RP,
     $     lHpHmhh,lHpHmaa,lHphc,lHpac,lHpHmcc)

      higgs = 0d0

      do i = 1,3
         do j=1,2
            
            higgs = higgs       ! charged+scalar and charged+pseudo
     $           +lHphc(i,j)**2*myB0(p**2,mhh(i)**2,mhc(j)**2,Q**2)
     $           +lHpac(i,j)**2*myB0(p**2,maa(i)**2,mhc(j)**2,Q**2)
         enddo

         higgs = higgs          ! scalar and pseudoscalar bubbles
     $        +lHpHmhh(i,i)*myA0(mhh(i)**2,Q**2)
     $        +lHpHmaa(i,i)*myA0(maa(i)**2,Q**2)
      enddo

      higgs = higgs             ! charged bubbles (note the symmetry factor)
     $     +lHpHmcc(1)*myA0(mhc(1)**2,Q**2)
     $     +4*lHpHmcc(2)*myA0(mhc(2)**2,Q**2)
               
*     the fermion contributions

      mt2 = ht**2*v2**2
      mb2 = hb**2*v1**2
      mtau2 = htau**2*v1**2
      
      fermions = 
     $     +3*((ht**2*cb**2+hb**2*sb**2)*myG(p**2,mt2,mb2,Q**2)
     $     -2*ht*hb*Sqrt(mt2*mb2)*2*sb*cb*myB0(p**2,mt2,mb2,Q**2))
     $     +htau**2*sb**2*myG(p**2,mtau2,0d0,Q**2)

*     the sfermion contributions

      call coupl_Hp_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     Rt,Rb,Rtau,lHpHmtt,lHpHmbb,lHpHmtata,lHpHmntnt,lHpHmuu,
     $     lHpHmdd,lHpHmee,lHpHmnn,lHptb,lHpntta,lHpud,lHpne)

      sfermions = 0d0           ! first two generations
      
      do i=1,2
         do j = 1,2

            sfermions = sfermions ! up-down squarks
     $          +6*lHpud(i,j)**2*myB0(p**2,msup(i)**2,msdown(j)**2,Q**2)
         enddo
         
         sfermions = sfermions  ! up-down sleptons
     $        +2*lHpne(i)**2*myB0(p**2,msel(i)**2,msnue**2,Q**2)
         
         sfermions = sfermions  ! charged-sfermion bubble
     $        +6*lHpHmuu(i,i)*myA0(msup(i)**2,Q**2)
     $        +6*lHpHmdd(i,i)*myA0(msdown(i)**2,Q**2)
     $        +2*lHpHmee(i,i)*myA0(msel(i)**2,Q**2)

      enddo
      
      sfermions = sfermions     ! sneutrino bubble
     $     + 2*lHpHmnn*myA0(msnue**2,Q**2)
      
      sfermions3g = 0d0         ! third generation

      do i=1,2
         do j = 1,2

            sfermions3g = sfermions3g ! top-bottom squarks
     $          +3*lHptb(i,j)**2*myB0(p**2,mstop(i)**2,msbot(j)**2,Q**2)
         enddo
         
         sfermions3g = sfermions3g ! stau-sneutrino
     $        +lHpntta(i)**2*myB0(p**2,mstau(i)**2,msnutau**2,Q**2)
         
         sfermions3g = sfermions3g ! charged-sfermion bubble
     $        +3*lHpHmtt(i,i)*myA0(mstop(i)**2,Q**2)
     $        +3*lHpHmbb(i,i)*myA0(msbot(i)**2,Q**2)
     $        +lHpHmtata(i,i)*myA0(mstau(i)**2,Q**2)
         
      enddo

      sfermions3g = sfermions3g           ! sneutrino bubble       
     $     +lHpHmntnt*myA0(msnutau**2,Q**2)
      
*     the chargino/neutralino contribution

      inos = 0d0

      call coupl_Hp_ino(g,gp,ll,v1,v2,NN,UU,VV,aHpnech,bHpnech)

      do i = 1,5
         do j = 1,2
            
            inos = inos
     $           +(aHpnech(i,j)**2+bHpnech(i,j)**2)
     $           *myG(p**2,mne(i)**2,mch(j)**2,Q**2)
     $           -4*aHpnech(i,j)*bHpnech(i,j)*
     $           mne(i)*mch(j)*myB0(p**2,mne(i)**2,mch(j)**2,Q**2)
         enddo
      enddo

      PiHpHm = weak + higgs + fermions + sfermions + sfermions3g + inos

c      write(*,*) 'gauge =',weak
c      write(*,*) 'higgs =',higgs
c      write(*,*) 'fermions =',fermions
c      write(*,*) 'sfermions1+2g =',sfermions
c      write(*,*) 'sfermions3g =',sfermions3g
c      write(*,*) 'inos =',inos
c      write(*,*) 'total =',PiHpHm

      PiHpHm = PiHpHm/16d0/pi**2

      end
      
*      
***********************************************************************
*     

      SUBROUTINE coupl_Hp_ino(g,gp,ll,v1,v2,NN,UU,VV,aHpnech,bHpnech)
      
      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ll,v1,v2,NN(5,5),UU(2,2),VV(2,2),
     $     aHpnech(5,2),bHpnech(5,2)

      DOUBLE PRECISION cb,sb,sq2

      INTEGER i,j
      
      cb = v1/sqrt(v1**2+v2**2)
      sb = v2/v1*cb

      sq2 = sqrt(2d0)

      do i=1,5
         do j=1,2
            aHpnech(i,j) = 0d0
            bHpnech(i,j) = 0d0
         enddo
      enddo

      do i=1,5
         do j=1,2
            aHpnech(i,j) = aHpnech(i,j)
     $           +sb*(-g*NN(i,3)*UU(j,1)
     $           +(g*NN(i,2)+gp*NN(i,1))*UU(j,2)/sq2)
     $           -cb*ll*NN(i,5)*UU(j,2)
            
            bHpnech(i,j) = bHpnech(i,j)
     $           +cb*(-g*NN(i,4)*VV(j,1)
     $           -(g*NN(i,2)+gp*NN(i,1))*VV(j,2)/sq2)
     $           -sb*ll*NN(i,5)*VV(j,2)
         enddo
      enddo
         
      end

*      
***********************************************************************
*     

      SUBROUTINE coupl_Hp_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     Rt,Rb,Rtau,lHpHmtt,lHpHmbb,lHpHmtata,lHpHmntnt,lHpHmuu,
     $     lHpHmdd,lHpHmee,lHpHmnn,lHptb,lHpntta,lHpud,lHpne)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     Rt(2,2),Rb(2,2),Rtau(2,2),lHpHmtt(2,2),lHpHmbb(2,2),
     $     lHpHmtata(2,2),lHpHmntnt,lHpHmuu(2,2),lHpHmdd(2,2),
     $     lHpHmee(2,2),lHpHmnn,lHptb(2,2),lHpntta(2),lHpud(2,2),
     $     lHpne(2),lHpHmtt_int(2,2),lHpHmbb_int(2,2),
     $     lHpHmtata_int(2,2),lHptb_int(2,2),lHpntta_int(2)

      DOUBLE PRECISION YuL,YuR,YdL,YdR,YeL,YeR,Ynu,cb,sb,c2b

      INTEGER i,j,k,l

      cb = v1/sqrt(v1**2+v2**2)
      sb = v2/v1*cb
      c2b = cb**2-sb**2

      YuL = 1/3d0 
      YdL = 1/3d0
      Ynu = -1d0
      YeL = -1d0
      YuR = -4/3d0
      YdR = 2/3d0
      YeR = 2d0

*     FIRST TWO GENERATIONS

*     quartic

      do i=1,2
         do j=1,2
            lHpHmuu(i,j) = 0d0
            lHpHmdd(i,j) = 0d0
            lHpHmee(i,j) = 0d0
         enddo
      enddo

      lHpHmuu(1,1) = c2b/4d0*(g**2+YuL*gp**2)
      lHpHmuu(2,2) = c2b/4d0*YuR*gp**2
      
      lHpHmdd(1,1) = -c2b/4d0*(g**2-YdL*gp**2)
      lHpHmdd(2,2) = c2b/4d0*YdR*gp**2
      
      lHpHmee(1,1) = -c2b/4d0*(g**2-YeL*gp**2)
      lHpHmee(2,2) = c2b/4d0*YeR*gp**2

      lHpHmnn = c2b/4d0*(g**2+Ynu*gp**2)
      
*     trilinear

      do i=1,2
         do j=1,2
            lHpud(i,j) = 0d0
         enddo
            lHpne(i) = 0d0            
      enddo

      lHpud(1,1) = g**2/2d0*(v2*cb+v1*sb)
      
      lHpne(1) = g**2/2d0*(v2*cb+v1*sb)

*     THIRD GENERATION (add the Yukawa interactions)

*     quartic

      do i=1,2
         do j=1,2
            lHpHmtt_int(i,j) = lHpHmuu(i,j)
            lHpHmbb_int(i,j) = lHpHmdd(i,j)
            lHpHmtata_int(i,j) = lHpHmee(i,j)
         enddo
      enddo

      lHpHmtt_int(1,1) = lHpHmtt_int(1,1) + hb**2*sb**2
      lHpHmtt_int(2,2) = lHpHmtt_int(2,2) + ht**2*cb**2

      lHpHmbb_int(1,1) = lHpHmbb_int(1,1) + ht**2*cb**2
      lHpHmbb_int(2,2) = lHpHmbb_int(2,2) + hb**2*sb**2

      lHpHmtata_int(2,2) = lHpHmtata_int(2,2) + htau**2*sb**2

      lHpHmntnt = lHpHmnn + htau**2*sb**2

*     trilinear

      lHptb_int(1,1) = lHpud(1,1) - ht**2*v2*cb - hb**2*v1*sb
      lHptb_int(2,2) = lHpud(2,2) - ht*hb*(v1*cb + v2*sb)
      lHptb_int(1,2) = lHpud(1,2) - hb*(ll*xx*cb + Ab*sb)
      lHptb_int(2,1) = lHpud(2,1) - ht*(ll*xx*sb + At*cb)

      lHpntta_int(1) = lHpne(1) - htau**2*v1*sb
      lHpntta_int(2) = lHpne(2) - htau*(ll*xx*cb + Atau*sb)
      
*     now rotate the third-generation couplings

      do i = 1,2
         do j = 1,2
            lHpHmtt(i,j) = 0d0
            lHpHmbb(i,j) = 0d0
            lHpHmtata(i,j) = 0d0
            lHptb(i,j) = 0d0
            do k=1,2
               do l=1,2
                  lHpHmtt(i,j) = lHpHmtt(i,j)
     $                 +Rt(i,k)*Rt(j,l)*lHpHmtt_int(k,l)
                  lHpHmbb(i,j) = lHpHmbb(i,j)
     $                 +Rb(i,k)*Rb(j,l)*lHpHmbb_int(k,l)
                  lHpHmtata(i,j) = lHpHmtata(i,j)
     $                 +Rtau(i,k)*Rtau(j,l)*lHpHmtata_int(k,l)
                  lHptb(i,j) = lHptb(i,j)
     $                 +Rt(i,k)*Rb(j,l)*lHptb_int(k,l)
               enddo
            enddo
         enddo
      enddo
      
*     stau-sneutrino trilinear

      do i = 1,2
         lHpntta(i) = 0d0
         do j=1,2
            lHpntta(i) = lHpntta(i) + Rtau(i,j)*lHpntta_int(j)
         enddo
      enddo

      end

*      
***********************************************************************
*     

      SUBROUTINE coupl_Hp_hh(g,gp,ll,kk,v1,v2,xx,Al,RS,RP,
     $     lHpHmhh,lHpHmaa,lHphc,lHpac,lHpHmcc)

      IMPLICIT NONE

      DOUBLE PRECISION g,gp,ll,kk,v1,v2,xx,Al,RS(3,3),RP(3,3),
     $     lHpHmhh(3,3),lHpHmaa(3,3),lHphc(3,2),lHpac(3,2),lHpHmcc(2)

      DOUBLE PRECISION sq2,c2b,s2b,
     $     lHpHmss(3,3),lHpHmpp(3,3),lHpsc(3,2),lHppc(3,2)

      INTEGER i,j,k,l
      
      c2b = (v1**2-v2**2)/(v1**2+v2**2)
      s2b = 2*v1*v2/(v1**2+v2**2)
      sq2 = sqrt(2d0)

      do i=1,3                  ! initialize
         do j=1,3
            lHpHmhh(i,j) = 0d0
            lHpHmaa(i,j) = 0d0
            lHpHmss(i,j) = 0d0
            lHpHmpp(i,j) = 0d0
            if(j.lt.3) then
               lHphc(i,j) = 0d0
               lHpac(i,j) = 0d0
               lHpsc(i,j) = 0d0
               lHppc(i,j) = 0d0   
            endif
         enddo
      enddo

*     quartic charged-neutral couplings
 
      lHpHmss(1,1) = (g**2-gp**2*c2b)/8d0
      lHpHmss(1,2) = -(2*ll**2-g**2)*s2b/8d0
      lHpHmss(2,1) = lHpHmss(1,2)
      lHpHmss(2,2) = (g**2+gp**2*c2b)/8d0 
      lHpHmss(3,3) = ll*(ll+kk*s2b)/2d0

      lHpHmpp(1,1) = (g**2-gp**2*c2b)/8d0
      lHpHmpp(1,2) = (2*ll**2-g**2)*s2b/8d0
      lHpHmpp(2,1) = lHpHmpp(1,2)
      lHpHmpp(2,2) = (g**2+gp**2*c2b)/8d0 
      lHpHmpp(3,3) = ll*(ll-kk*s2b)/2d0
      
*     rotate the quartics

      do i=1,3              
         do j=1,3
            do k=1,3
               do l=1,3
                  lHpHmhh(i,j) = lHpHmhh(i,j)
     $                 +RS(i,k)*RS(j,l)*lHpHmss(k,l)
                  lHpHmaa(i,j) = lHpHmaa(i,j)
     $                 +RP(i,k)*RP(j,l)*lHpHmpp(k,l)
               enddo
            enddo
         enddo
      enddo
      
*     trilinear charged-neutral couplings

      lHpsc(1,1) = (-v1*gp**2*s2b+v2*(2*ll**2-g**2)*c2b)/2d0/sq2
      lHpsc(1,2) = (v1*(g**2-gp**2*c2b)-v2*(2*ll**2-g**2)*s2b)/2d0/sq2
      lHpsc(2,1) = (v2*gp**2*s2b+v1*(2*ll**2-g**2)*c2b)/2d0/sq2
      lHpsc(2,2) = (v2*(g**2+gp**2*c2b)-v1*(2*ll**2-g**2)*s2b)/2d0/sq2
      lHpsc(3,1) = -ll/sq2*(Al+2*kk*xx)*c2b
      lHpsc(3,2) = ll/sq2*(2*ll*xx+(Al+2*kk*xx)*s2b)

      lHppc(1,1) = v2*(2*ll**2-g**2)/2d0/sq2               
      lHppc(2,1) = v1*(2*ll**2-g**2)/2d0/sq2               
      lHppc(3,1) = ll/sq2*(Al-2*kk*xx)

*     rotate the trilinears

      do i=1,3              
         do j=1,2
            do k=1,3
               lHphc(i,j) = lHphc(i,j) + RS(i,k)*lHpsc(k,j)
               lHpac(i,j) = lHpac(i,j) + RP(i,k)*lHppc(k,j)
            enddo
         enddo
      enddo

*     quartic charged-charged couplings

c      lHpHmcc(1) = -(g**2+gp**2)*(c2b**2-s2b**2)/4d0
c      lHpHmcc(2) = (g**2+gp**2)*c2b**2/8d0   
      lHpHmcc(1) = -(g**2+gp**2)*(c2b**2-s2b**2)/4d0 + ll**2*c2b**2
      lHpHmcc(2) = (g**2+gp**2)*c2b**2/8d0 + ll**2/4d0*s2b**2  

      end

