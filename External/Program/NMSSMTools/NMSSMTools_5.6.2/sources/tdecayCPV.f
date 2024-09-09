        SUBROUTINE TDECAY_CPV()

c ==================================================================== c
c                           top 2-body decays                          c
c ==================================================================== c
c Incl. dominant rad. corrs.:
c top --> bottom + W+/-: 
c   1-loop QCD as in Li, Oakes, Yuan, PRD43 (1991) 3759 for mb -> 0
c   no Susy (<~ 4%)
c top --> bottom + H+/-:
c   1-loop QCD as in Czarnecki, Davidson, hep-ph/9301237, for mb -> 0
c      (BUT: with terms ~ mb*tanb/mt)
c   Susy: ~delta_mb/mb only, from Guasch et al., hep-ph/9507461
c      and hep-ph/0003109
c ==================================================================== c

      IMPLICIT NONE 

      INTEGER I,J

      DOUBLE PRECISION topneutrstop(5,2)
      DOUBLE PRECISION topbw,topbh
      DOUBLE PRECISION sw,cw,tw
      DOUBLE PRECISION gmst(2),lamb_funct,rmt,rmb
      DOUBLE PRECISION alsmt,runmb,sp,xmb,xmw,xmh
      DOUBLE PRECISION dtbwqcd,gp,gm,nb,nllr
      DOUBLE PRECISION PI,SQR2,MHC

      DOUBLE PRECISION MHC2,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      DOUBLE PRECISION DELMB,DELML,DEL1
      DOUBLE PRECISION MNEU(5),NEU(5,5,2)
      DOUBLE PRECISION MST2P(2),MSB2P(2),MSU2P(2),MSD2P(2)
      DOUBLE PRECISION CONESTtL(5,2,2),CONESTtR(5,2,2),
     . CONESBbL(5,2,2),CONESBbR(5,2,2),CONESLlL(5,2,2),CONESLlR(5,2,2),
     . CONESNnL(5,2),CONESUuL(5,2,2),CONESUuR(5,2,2),
     . CONESDdL(5,2,2),CONESDdR(5,2,2),CONESEeL(5,2,2),
     . CONESEeR(5,2,2)
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2)
      DOUBLE PRECISION toptot

      COMMON/HISPEC/MHC2,XC,MH0,XH,MA2
      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/SMFERM/mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/NEUSPEC/MNEU,NEU
      COMMON/SFERMPSPEC/MST2P,MSB2P,MSU2P,MSD2P
      COMMON/NEUSFfCOUP/CONESTtL,CONESTtR,CONESBbL,CONESBbR,CONESLlL,
     . CONESLlR,CONESNnL,CONESUuL,CONESUuR,CONESDdL,CONESDdR,CONESEeL,
     . CONESEeR
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/topwidth/toptot


      PI=4d0*DATAN(1d0)
      sqr2=dsqrt(2d0)
      MHC=dsqrt(MHC2)

* Parameters at M_top:

        ALSMT=ALSMZ/(1d0+23d0/(12d0*PI)*ALSMZ*DLOG(MT**2/MZ**2))
        rmt=MT/(1d0+4d0*ALSMT/(3d0*PI)+11d0*(ALSMT/PI)**2)
        rmb=RUNMB(MT)

        xmb=(rmb/rmt)**2
        xmw=(mw/rmt)**2
        xmh=(mhc/mt)**2

c -------------------------------------------------------------------- c
c top --> bottom + W+/-

      IF(mt.gt.(mb+MW)) THEN
c   QCD correction (mb -> 0):
        dtbwqcd=2d0*alsmt/(3d0*pi)*(
     .    2d0*xmw*(2d0*xmw-1d0)*(1d0+xmw)*dlog(xmw)/
     .           ((1d0-xmw)**2*(1d0+2d0*xmw))
     .    -(5d0+4d0*xmw)/(1d0+2d0*xmw)*dlog(1-xmw)
     .    +2d0*(sp(1d0-xmw)-sp(xmw))-pi**2
     .    +(5d0+9d0*xmw+6d0*xmw**2)/(2d0*(1d0-xmw)*(1d0+2d0*xmw)))
         topbw = GF/(8d0*pi*sqr2)*mt**3*lamb_funct(mb/mt,mw/mt)*
     .     ((1-xmb)**2+(1+xmb)*xmw-2d0*xmw**2)*(1d0+dtbwqcd)
      ELSE
         topbw = 0d0
      ENDIF

c -------------------------------------------------------------------- c
c top --> bottom + H+/-

      IF(mt.gt.(mb+MHC)) THEN
c  QCD corrections (mb -> 0):
        gp=(1-xmh)*(sp(1-xmh)+9d0/8d0-pi**2/3d0
     .    +dlog(xmh)*dlog(1-xmh)/2d0+xmh*dlog(xmh)/(xmh-1)/2d0
     .    +(1/xmh/2d0-5d0/4d0)*dlog(1-xmh)+3d0/8d0*dlog(xmb))
        gm=-3d0/4d0*(1-xmh)*dlog(xmb)
c For tree level + Susy corrections:
        nb=(1+xmb-xmh)*(1/tanb**2+xmb*tanb**2)+4d0*xmb
        nllr=(1d0+xmb-xmh)*xmb*tanb**2+2d0*xmb
c Susy corrections: DELMB
c
         topbh = GF/(8d0*pi*sqr2)*mt**3*lamb_funct(mb/mt,mhc/mt)*(nb
     .     -2d0*nllr*DELMB
     .     +4d0*alsmt/(3d0*pi)*(2d0*(1/tanb**2+xmb*tanb**2)*gp
     .     +(1/tanb**2-xmb*tanb**2)*gm))
      ELSE
         topbh = 0d0
      ENDIF

c -------------------------------------------------------------------- c
c top --> stop(j) + chi^0_i from Sdecay:

      DO i=1,5
         DO j=1,2
               topneutrstop(i,j) = 0d0
         ENDDO
      ENDDO

       IF(dsqrt(mneu(1))+dsqrt(MST2P(1)).lt.mt) THEN

         cw=DSQRT(G2/(G1+G2))
         sw=DSQRT(G1/(G1+G2))
         tw=sw/cw

        gmst(1) = dsqrt(MST2P(1))
        gmst(2) = dsqrt(MST2P(2))

      DO i=1,5
         DO j=1,2
            IF(mt.gt.(dsqrt(mneu(i))+gmst(j))) THEN
               topneutrstop(i,j) = 1d0/32d0/pi/mt*(
     .    (CONESTtL(I,J,1)*CONESTtR(I,J,1)
     .     -CONESTtL(I,J,2)*CONESTtR(I,J,2))*4d0*mt*dsqrt(mneu(i))
     .    +(CONESTtL(I,J,1)**2+CONESTtR(I,J,1)**2
     .     +CONESTtL(I,J,2)**2+CONESTtR(I,J,2)**2)*
     .              (mt**2-gmst(j)**2+mneu(i)) )*
     .              lamb_funct(dsqrt(mneu(i))/mt,gmst(j)/mt)
            ELSE
               topneutrstop(i,j) = 0d0
            ENDIF
         ENDDO
      ENDDO
      
      ENDIF
c -------------------------------------------------------------------- c

         toptot = topbw+topbh
         DO i=1,5
            DO j=1,2
               toptot = toptot + topneutrstop(i,j)
            ENDDO
         ENDDO
         
         brtopbw = topbw/toptot
         brtopbh = topbh/toptot

         DO i=1,5
            DO j=1,2
               brtopneutrstop(i,j) = topneutrstop(i,j)/toptot
            ENDDO
         ENDDO

      END
